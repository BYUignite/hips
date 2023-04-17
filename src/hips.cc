#include "hips.h"

#include <iostream>
#include <cmath>
#include <iomanip>
#include "randomGenerator.h"
#include <fstream>
#include <sstream>

using namespace std;

///////////////////////////////////////////////////////////////////////////////
/** Constructor
  */

hips::hips(int     nLevels_, 
           double  domainLength_, 
           double  tau0_, 
           double  C_param_, 
           int     forceTurb_,
           int     nVar_,
           vector<double> &ScHips_,
           shared_ptr<Cantera::Solution> cantSol) : 
    nLevels(nLevels_), domainLength(domainLength_), tau0(tau0_),
    C_param(C_param_), forceTurb(forceTurb_),       ScHips(ScHips_),   
    nVar(nVar_),       LrandSet(true),              cvodeD(cantSol) { 

    varData = vector<vector<double> * > (nVar);                  // or varData.resize(nVar)
   
    gas = cantSol->thermo(); 
    nsp = gas->nSpecies();                                           // read the size of species

    //-------------------------- Set number of parcels, level lengthscales, timescales, and rates, and i_plus && i_batchelor
     
    iEta = nLevels - 3;        // Kolmogorov level; if nLevels = 7, then 0, 1, 2, 3, (4), 5, 6; iEta=4 is the lowest swap level: swap grandchildren of iEta=4 at level 6.
           
    int maxSc = 1.0;

    for (int i=0; i<ScHips.size(); i++)
        maxSc = ScHips[i]>maxSc ? ScHips[i] : maxSc;
    
    if (maxSc > 1.0)
        nLevels += ceil(log(maxSc)/log(4));       // Changing number of levels!
    
    Nm1 = nLevels - 1;
    Nm2 = nLevels - 2;
    Nm3 = nLevels - 3;
    
    // -------------------------- 
    
    nparcels = static_cast<int>(pow(2, Nm1));
    parcelTimes.resize(nparcels,0);
    i_batchelor.resize(nVar,0);
    
    vector<double> levelLengths(nLevels);      // including all levels, but last 2 don't count:
    vector<double> levelTaus(nLevels);         // smallest scale is 2 levels up from bottom
    levelRates   = vector<double>(nLevels);

    for (int i=0; i<nLevels; i++) {
        levelLengths[i] = domainLength * pow(Afac,i);
        levelTaus[i] = tau0 * pow(levelLengths[i]/domainLength, 2.0/3.0) / C_param;
        levelRates[i] = 1.0/levelTaus[i] * pow(2.0,i);
    }
    
    LScHips = ScHips.size() > 0 ? true : false;
    if (LScHips) {                             // correct levels for high Sc (levels > Kolmogorov)
        for (int i=iEta+1; i<nLevels; i++) {
            levelTaus[i] = tau0 *
            pow(levelLengths[iEta]/domainLength, 2.0/3.0) /
            C_param;
            levelRates[i] = 1.0/levelTaus[i] * pow(2.0,i);
        }
    }

    //-------------------------------------------------

    eddyRate_total = 0.0;
    for (int i=0; i<=Nm3; i++)
        eddyRate_total += levelRates[i];
    
    eddyRate_inertial = 0.0;
    for (int i=0; i<=iEta; i++)
        eddyRate_inertial += levelRates[i];
    
    //-------------------
    
    i_plus.resize(nVar);
    for (int k=0; k<nVar; k++) {
        if (ScHips[k] < 1.0)
            i_batchelor[k] = iEta + 1.5*log(ScHips[k])/log(4);
        else if (ScHips[k] > 1.0)
            i_batchelor[k] = iEta +     log(ScHips[k])/log(4);
        else
            i_batchelor[k] = iEta;

        i_plus[k] = ceil(i_batchelor[k]);
    }
    
    //------------------- Set the parcel addresses (index array)
    
    pLoc.resize(nparcels);
    for (int i=0; i<nparcels; i++)
        pLoc[i] = i;
    
    //-------------------------------------

    rand = new randomGenerator(seed);
} 

///////////////////////////////////////////////////////////////////////////////
/** The HiPS solver
 * @param tRun \input simulation run time
 *
 * Sample tee (time of next eddy event)
 * Select and swap subtrees (at current time, not at tee, so that we know whos involved)
 * If the eddy event is at the parcel/micromixing level:
 *     React involved parcels from their current time to tee
 *     Mix the involved parcels (micromixing)
 */

void hips::calculateSolution(const double tRun) {
    //cout<<"cal test 1"<<endl;
    int    nEddies = 0;                // number of eddy events
    int    fileCounter = 0;            // number of data files written
    int    iLevel;                     // tree level of EE with top at iLevel=0
    int    iTree;                      // one of two subtrees involved in swap at iLevel
    double dtEE;                       // time increment to next eddy event 

    time = 0.0;                        // initialize simulation time

    sample_hips_eddy(dtEE, iLevel);    // get first EE at time 0+dtEE
    nEddies++;
   //cout<<"cal test 2"<<endl;

    while (time+dtEE<=tRun) {

        time += dtEE;
        selectAndSwapTwoSubtrees(iLevel, iTree);
        // cout<<"cal test 3"<<endl;
        advanceHips(iLevel, iTree);    // reaction and micromixing (if needed) to t=time
      // cout<<"cal test 4"<<endl;


        sample_hips_eddy(dtEE, iLevel);
       // cout<<"cal test 5"<<endl;

        nEddies++;

        if(nEddies %200 == 0) writeData(++fileCounter, time);
    }
    time = tRun;
    iLevel = 0; iTree  = 0;
    reactParcels_LevelTree(iLevel, iTree);      // react all parcels up to end time
}


///////////////////////////////////////////////////////////////////////////////
/** Samples stochastic eddy events on the hips tree: time and level.
*   @param dtEE   \output time increment to next eddy event (EE)
*   @param iLevel \output tree level of EE
*/

void hips::sample_hips_eddy(double &dtEE, int &iLevel) {

    static double c1 = 1.0 - pow(2.0, 5.0/3.0*(iEta+1));
    static double c2 = pow(2.0, Nm2) - pow(2.0, iEta+1);
    static double c3 = pow(2.0, iEta+1);

    //--------------- time to next eddy

    double r = rand->getRand();
    dtEE = -log(r)/eddyRate_total;

    //----------------- get eddy level

    r = rand->getRand();

    if ( r <= eddyRate_inertial/eddyRate_total) {     // inertial region
        r = rand->getRand();
        iLevel = ceil(3.0/5.0*log2(1.0-r*c1) - 1.0);
        if (iLevel < 0)    iLevel = 0;
        if (iLevel > iEta) iLevel = iEta;
    }

    else {                                            // "Batchelor" region
        r = rand->getRand();
        iLevel = ceil(log2(r*c2 + c3) - 1.0);
        if(iLevel < iEta+1) iLevel = iEta+1;
        if(iLevel > Nm3) iLevel = Nm3;
    }

    return;
}

///////////////////////////////////////////////////////////////////////////////
/** Function performs eddy events: parcel swaps.
    @param iLevel \input  level of the tree for the base of the swap.
    @param iTree  \output which subtree on the level is selected

    Randomly select a node on iLevel.
    Go down two levels and select nodes 0q and 1r, where q, r are randomly 0 or 1
    Find the starting index of the Q-tree and R-tree to swap and the number of parcels.
    Then swap the cells.

    For a 6 level tree: 0, 1, 2, 3, 4, 5:
    If iLevel = 1, then suppose i=1, 0q = 00 and 1r = 11:
    Then we are swaping 0100** with 0111** or (01|00|**) with (01|11|**)
       or i0qs with i1rs, where i = 01; 0q = 00; 1r = 11; and s = **

    We use bitwise shifts for easy powers of 2.
    The swap is done by adding or subtracting a value (shift),
        which should be equivalent to flipping the swapping the two 0q bits and 1r bits.
                                                                                                              Level
                                                                                                            ---------
                                                    *                                                           0
                                                 /     \
                                              /           \
                                           /                 \
                                        /                       \
                                     /                             \
                                  /                                   \
                               /                                         \
                            /                                               \
                           *                                                (*)  01|0000                        1
                          / \                                               / \
                        /     \                                           /     \
                      /         \                                       /         \
                    /             \                                   /             \
                  /                 \                               /                 \
                /                     \                           /                     \
               *                       *                         *                       *                      2
              / \                     / \                       / \                     / \
            /     \                 /     \                   /     \                 /     \
          /         \             /         \               /         \             /         \
         *           *           *           *            [*] 00|**    *           *          [*] 11|**         3
        / \         / \         / \         / \           / \         / \         / \         / \
       /   \       /   \       /   \       /   \         /   \       /   \       /   \       /   \
      *     *     *     *     *     *     *     *       *     *     *     *     *     *     *     *             4
     / \   / \   / \   / \   / \   / \   / \   / \     / \   / \   / \   / \   / \   / \   / \   / \
    00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15   16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31           5
                                                      ^^^^^^^^^^^                         ^^^^^^^^^^^

*/

void hips::selectAndSwapTwoSubtrees(const int iLevel, int &iTree) {

    iTree = rand->getRandInt((1 << iLevel)-1);
    int zero_q = rand->getRandInt(1);                                    // 0q where q is 0 or 1
    int one_r  = 2 + rand->getRandInt(1);                                // 1r where r is 0 or 1

    int Qstart = (zero_q << (Nm3-iLevel)) + (iTree << (Nm1-iLevel));     // starting index of Q parcels
    int Rstart = (one_r  << (Nm3-iLevel)) + (iTree << (Nm1-iLevel));     // starting index of R parcels
    int nPswap = 1 << (Nm3-iLevel);                                      // number of parcels that will be swapped

    int Qend = Qstart + nPswap;                                          // inclusive indices are Qstart to Qend-1
    int Rend = Rstart + nPswap;                                          // inclusive indices are Rstart to Rend-1
    vector<int> aa(pLoc.begin()+Qstart, pLoc.begin()+Qend);
    copy(pLoc.begin()+Rstart, pLoc.begin()+Rend, pLoc.begin()+Qstart);   // python: pLoc[Qstart:Qend]=pLoc[Rstart:Rend]
    copy(aa.begin(), aa.end(), pLoc.begin()+Rstart);                     // python: pLoc[Rstart:Rend]=aa
}


///////////////////////////////////////////////////////////////////////////////
/** React and mix parcels that are involved in a micromixing process.
  * This is determined by the level and the tree within that level.
  * Might not do anything if the eddy doesn't cause micromixing.
  * @param iLevel \input level that the eddy event occurred.
  * @param iTree  \input root note of the eddy event at iLevel.
  */

void hips::advanceHips(const int iLevel, const int iTree) {
     
    if (forceTurb==2 && iLevel==0)              // forcing for statistically stationary
        forceProfile();

    bool rxnDone = false;                       // react all variables once
    for (int k=0; k<nVar; k++) {                //   upon finding first variable needing micromixing
        if ( (iLevel >= i_plus[k]) || 
             (iLevel==i_plus[k]-1 && rand->getRand() <= i_plus[k]-i_batchelor[k]) ) {
                if(!rxnDone) {
                   reactParcels_LevelTree(iLevel, iTree);
                   rxnDone = true;
             }
            mixAcrossLevelTree(k, iLevel, iTree);
        }
    }

    }

///////////////////////////////////////////////////////////////////////////////
/** React parcels that are involved in a micromixing process.
  * Parcels might react for different amounts of time depending on when they last
  * reacted, which is stored in the parcelTimes array.
  * @param iLevel \input level that the eddy event occurred.
  * @param iTree  \input root note of the eddy event at iLevel.
  */

void hips::reactParcels_LevelTree(const int iLevel, const int iTree) {
    
    int nP  = 1 << (Nm1 - iLevel);
    int istart = iTree << nP;
    int iend = istart + nP;
    int ime;
    double dt;

    for (int i=istart; i<iend; i++) {
        ime = pLoc[i];

        setState(ime);
        dt = time-parcelTimes[ime];
        cvodeD.integrate(0.01);
        // cvodeD.integrate(dt);

        varData[0][0][ime] = gas->temperature();
        for (int k=0; k<nsp; k++)
            varData[k+1][0][ime] = gas->massFraction(k);

//          varData[0][0][ime] = gas->temperature(); 
//          varData[1][0][ime] = gas->massFraction(0);
//          varData[2][0][ime] = gas->massFraction(1);
//          varData[3][0][ime] = gas->massFraction(2);
//          varData[4][0][ime] = gas->massFraction(3);
//          varData[5][0][ime] = gas->massFraction(4);
//          varData[6][0][ime] = gas->massFraction(5);
//
        // todo: set the reactor initial condition with state in parcel ime
        // todo: integrate the reactor for Δt = time - parcelTimes[ime]
        // todo: set the parcel state based on the end time of the reactor
        // todo: set parcelTimes[ime] = time
    }

}

///////////////////////////////////////////////////////////////////////////////
/** Mix parcels uniformly (using average) at given level and tree
  * Mixing at levels above the lowest enables low Sc variables.
  * @param kVar   \input  variable index to mix (normally a transported var as determined by caller)
  * @param iLevel \input  grandchildren of this iLevel will be mixed
  * @param iTree  \input  at the given iLevel, mix only this subtree
  *
  * Example: For a 5 level tree, we have levels 0, 1, 2, 3, 4 (top to bottom).
  *   
  *          Let iLevel = 2, then nPmix = 2 and we are mixing pairs of parcels at the base of the tree.
  *          iTree can be 0, 1, 2, or 3.
  *          if iTree = 0 we will mix parcels (0,1) and (2,3) as the left subtree and right subtree:
  *                (istart=1, iend=2) and (istart=2, iend=4)
  *          if iTree = 1 we will mix parcels (4,5) and (6,7) with (istart=4,iend=6), (istart=6, iend=8)
  *          if iTree = 2 we will mix parcels (8,9) and (10,11) with (istart=8,iend=10), (istart=10, iend=12)
  *          if iTree = 3 we will mix parcels (12,13) and (14,15) with (istart=12,iend=14), (istart=14, iend=16)
  *         
  *          Let iLevel=1 then we will be mixing groups of four parcels.
  *          iTree can be 0 or 1
  *          if iTree = 0 we will mix parcels (0,1,2,3) and (4,5,6,7) with (istart=0,iend=4), (istart=4, iend=8)
  *          if iTree = 1 we will mix parcels (8,9,10,11) and (12,13,14,15) with (istart=8,iend=12), (istart=12, iend=16)
  *
  * recall: 3 << 4 means 3*2^4 (or 3 = 000011 and 3<<4 = 110000 = 48), that is, we shift the bits left 4 places.
  *          
  * NOTE: BE CAREFUL WITH MIXING SOME SCALARS, LIKE MASS FRACTIONS; CURRENT CODE ASSUMES ALL PARCELS HAVE SAME DENSITY (mixing Yi directly)
  */

void hips::mixAcrossLevelTree(int kVar, const int iLevel, const int iTree) {
    
    int istart;
    int iend;

    int nPmix = 1 << (nLevels - iLevel - 2);   // number of parcels mixed together

    int ime;

    //---------- Mix left branch of iTree

    istart = iTree << (Nm1-iLevel);  
    iend = istart + nPmix;

    double s = 0;                   // initialize sum to 0
    for (int i=istart; i<iend; i++) {
        ime = pLoc[i];
        s += varData[kVar][0][ime];
    }
    for (int i=istart; i<iend; i++) {
        ime = pLoc[i];
        varData[kVar][0][ime] = s / nPmix; 
    }

    //--------- Mix right branch of iTree

    istart = iend;
    iend = istart + nPmix;

    s = 0;                   // initialize sum to 0
    for (int i=istart; i<iend; i++) {
        ime = pLoc[i];
        s += varData[kVar][0][ime];
    }
    for (int i=istart; i<iend; i++) {
        ime = pLoc[i];
        varData[kVar][0][ime] = s / nPmix;  
    }   
}

///////////////////////////////////////////////////////////////////////////////
/** Force hips profile to get statistically stationary
  * Meant to demonstrate forcing for simple scalars.
  * Such as, for a mixture fraction variable that varies between 0 and 1.
  * Code here is setup to force the left half of parcels to average 0,
  * and the right half of parcels to average 1.
  */

void hips::forceProfile() {

    for (int k=0; k<varData.size(); k++) {
        double s;

        //---------- force the left half of parcels to average 0

        s = 0.0;
        for (int i=0; i<nparcels>>1; i++)
            s += varData[k][0][pLoc[i]];
        
        s /= (nparcels>>1);
        for (int i=0; i<nparcels>>1; i++)
            varData[k][0][pLoc[i]] += (-s - 0.0);
         
        //---------- force the right half of parcels to average 1

        s = 0.0;
        for (int i=nparcels>>1; i<nparcels; i++)
            s += varData[k][0][pLoc[i]];

        s /= (nparcels>>1);
        for (int i=nparcels>>1; i<nparcels; i++)
            varData[k][0][pLoc[i]] += (-s + 1.0);   
    }
}

///////////////////////////////////////////////////////////////////////////////
/** Write a data file
  * @param ifile \input index used in the file name (sequential incrementing: Data_00001.dat, Data_00002.dat, etc.)
  * @param outputTime \input time of data file, written in file.
  */

void hips::writeData(const int ifile, const double outputTime) {

    stringstream ss1;
    string       s1;
    ss1.clear(); ss1<<setfill('0') << setw(5) << ifile; ss1 >> s1;
    string fname = "../data/Data_" + s1 + ".dat";

    ofstream ofile(fname.c_str());

    ofile << "# time = " << outputTime;
    ofile <<"\n# Grid Poins = " << nparcels;

    for (int i=0; i<nVar; i++){
        ofile << setw(14) <<  "variable_"+std::to_string(i);

    }

    ofile << scientific;
    ofile << setprecision(10);

    for (int i=0; i<nparcels; i++) {
        ofile<<endl; 
 
        for (int k=0; k<nVar; k++)
            ofile<<setw(19) <<varData[k][0][pLoc[i]];
    }

    ofile.close();
}

///////////////////////////////////////////////////////////////////////////////
/** Set the Cantera gas state.
  * @param ipt \input parcel to set the state for.
  */

void hips::setState(const int &ipt){

    vector<double> yi(nsp);
    double T = varData[0][0][ipt]; 
    
    for (int k=0; k<nsp; k++)
        yi[k] = varData[k+1][0][ipt];

    gas->setState_TPY(T, OneAtm, &yi[0]);
    //cout<<"test in state is "<<endl;
    
}

///////////////////////////////////////////////////////////////////////////////
