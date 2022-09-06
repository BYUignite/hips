/**
 * @file micromixer.cc
 * @brief Source file for class \ref micromixer
 */

#include "micromixer.h"
#include "domain.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <numeric> //accumulate

#include "interp_linear.h"


///////////////////////////////////////////////////////////////////////////////
/** micromixer constructor function
 */

void micromixer::setNominalStepSize() {
    dtStepNominal = domn->pram->diffCFL * domn->solv->tMix;
}


micromixer::micromixer() {
    cvode = new cvodeDriver();
    nsteps = 0;
   // if(LcvodeSet) delete cvode;

}

///////////////////////////////////////////////////////////////////////////////
/** micromixer initialization function
 *
 * @param p_domn  \input set domain pointer with.
 */

void micromixer::init(domain *p_domn) {
    domn = p_domn;

    bool LincludeRhsMix = (domn->pram->Lsolver=="SEMI-IMPLICIT") ? true : false;
    cvode->init(domn, LincludeRhsMix);
}

///////////////////////////////////////////////////////////////////////////////

/** Mix parcels uniformly (using average) at given level
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

void micromixer::mixAcrossLevelTree(const int kVar, const int iLevel, const int iTree) {

    int istart;
    int iend;

    int nPmix = 1<<(domn->pram->nLevels - iLevel - 2);     // number of parcels mixed together
                                                           // also = # to skip in loop
                                  
    int ime;            // index of first cell (mapped to istart)
    double s;           // initialize sum to 0

    bool LdoChi = (domn->v[kVar]->var_name == "mixf" || 
                   domn->v[kVar]->var_name == "mixf_000") ? true : false;

    // TODO: generalize this to compute chi for any scalar, not just mixf, mixf_0

    //--------- Mix left branch of iTree

    //istart = iTree << (Nm1-iLevel);  // same thing

    istart = iTree << (domn->solv->Nm1-iLevel);  // same thing
    iend   = istart + nPmix;

    s = 0;            // initialize sum to 0
    for(int i=istart; i<iend; i++) {
        ime = domn->solv->pLoc[i];
        s += domn->v[kVar]->d[ime];
    }
    for(int i=istart; i<iend; i++) {
        ime = domn->solv->pLoc[i];

        if(LdoChi)
            domn->chi->setVar(domn->v[kVar]->d[ime], s/nPmix, ime);

        domn->v[kVar]->d[ime] = s / nPmix;
       

    }

    //--------- Mix right branch of iTree

    istart = iend;        // same thing
    iend   = istart + nPmix;

    s = 0;            // initialize sum to 0
    for(int i=istart; i<iend; i++) {
        ime = domn->solv->pLoc[i];
        s += domn->v[kVar]->d[ime];
    }
    for(int i=istart; i<iend; i++) {
        ime = domn->solv->pLoc[i];

        if(LdoChi)
            domn->chi->setVar(domn->v[kVar]->d[ime], s/nPmix, ime);

        domn->v[kVar]->d[ime] = s / nPmix;
    }
}

///////////////////////////////////////////////////////////////////////////////////////

/** Advance ODT solution: diffusion and reaction
  * @param tstart  \input start time for mixing
  * @param tend    \input end   time for mixing
  * @param iLevel  \input root node of the eddy (default value -1, see header file)
  *                       in case of iLevel=-1, no simple mixing is done, only advancement
  */

void micromixer::advanceHips(const double p_tstart, const double p_tend, const int iLevel) {

 tstart = p_tstart;
    tend   = p_tend;

    cout << endl << "here 1" << endl; //doldb 
    if(domn->pram->forceHips==2 && iLevel==0)   // forcing for statistically stationary
        forceProfile();
    cout << endl << "here 2" << endl; //doldb 


    if(domn->pram->LsimpleMix){
        for(int k=0; k<domn->v.size(); k++){
    cout << endl << "here 3 " << domn->v[k]->var_name << " " << domn->v.size() << " " << endl; //doldb 
            if(!domn->v[k]->L_transported)
                continue;
            if(iLevel == domn->v[k]->i_plus-1 &&                                       // mix scalar across level i_minus
               domn->rand->getRand() <= domn->v[k]->i_plus-domn->v[k]->i_batchelor)    // with probability i_plus - i_batchelor
                
                    mixAcrossLevelTree(k, iLevel,domn->solv-> iTree);  
            else if(iLevel >= domn->v[k]->i_plus ){                                    // scalars at or below i_plus are fully mixed
               
                mixAcrossLevelTree(k, iLevel,domn->solv-> iTree);                              //    the initial condition has uniform mixing below Batchelor
            }
        }
    }

    //if(domn->pram->LsimpleMix){{{
    //    for(int k=0; k<domn->v.size(); k++)
    //        if(domn->v.at(k)->L_transported){
    //            for(int i=0; i<domn->ngrd; i+=2){
    //                int ime = domn->solv->pLoc[i];
    //                int inb = domn->solv->pLoc[i+1];
    //                double val = (domn->v[k]->d[ime] + domn->v[k]->d[inb])/2.0;
    //                domn->v[k]->d[ime] = val;
    //                domn->v[k]->d[inb] = val;
    //            }
    //        }}}}

    setNominalStepSize();

    for(time=tstart; time<tend; time+=dt) {

        setStepSize();
        if(domn->pram->Lsolver=="EXPLICIT")
            advanceHipsSingleStep_Explicit();
        else if(domn->pram->Lsolver=="SEMI-IMPLICIT")
            advanceHipsSingleStep_SemiImplicit();
        else if(domn->pram->Lsolver=="STRANG")
            advanceHipsSingleStep_StrangSplit();

        domn->io->dumpDomainIfNeeded();

    }
   
}


///////////////////////////////////////////////////////////////////////////////

/** Force hips profile to get statistically stationary
  */

void micromixer::forceProfile() {

    for(int k=0; k<domn->v.size(); k++){
        if(!domn->v[k]->L_transported)
            continue;

        double s;
        //---------- force the left half of parcels to average 0

        s=0.0;
        for(int i=0; i<domn->ngrd>>1; i++)
            s += domn->v[k]->d[domn->solv->pLoc[i]];
        s /= (domn->ngrd>>1);
        for(int i=0; i<domn->ngrd>>1; i++)
            domn->v[k]->d[domn->solv->pLoc[i]] += (-s - 0.0);

        //---------- force the right half of parcels to average 1

        s=0.0;
        for(int i=domn->ngrd>>1; i<domn->ngrd; i++)
            s += domn->v[k]->d[domn->solv->pLoc[i]];
        s /= (domn->ngrd>>1);
        for(int i=domn->ngrd>>1; i<domn->ngrd; i++)
            domn->v[k]->d[domn->solv->pLoc[i]] += (-s + 1.0);
    }
}


////////////////////////////////////////////////////////////////////////////////
/** Set time step size.
 *  This is based on a diffusive (or other) timescale.
 *  This is a uniform step size for the given integration period.
 *  The actual step size will be rest based on a data dump or tend.
 */

void micromixer::setStepSize() {

    dt = dtStepNominal;

    if (time+dt > domn->io->dumpTimes.at(domn->io->iNextDumpTime)) {
        dt = domn->io->dumpTimes.at(domn->io->iNextDumpTime) - time;
        domn->io->LdoDump = true;
    }

    if (time + dt > tend) {
        dt = (tend - time)*(1.0+1.0E-8);
        domn->io->LdoDump = false;
    }
}

///////////////////////////////////////////////////////////////////////////////
/** Advance ODT solution: diffusion and reaction
 */

void micromixer::advanceHipsSingleStep_Explicit(){

 


    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported) {
            domn->v.at(k)->getRhsMix(gf, dxc);
            domn->v.at(k)->getRhsSrc();
        }

    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            for(int i=0; i < domn->ngrd; i++)
                domn->v.at(k)->d.at(i) = domn->v.at(k)->d.at(i) + dt*( domn->v.at(k)->rhsMix.at(i) + domn->v.at(k)->rhsSrc.at(i));

  


}

///////////////////////////////////////////////////////////////////////////////
/** Advance HiPS solution: diffusion and reaction; Some terms are implicit, others explicit.
 *  Nominally the mixing terms are explicit. Calling the cvode driver.
 *  First order.
 *  dphi/dt =  D(phi_0) + S(phi) : solving from t0 to t1.
 *  Here, D() is the diffusive term, and S() is the (stiff) source term.
 *  We solve the whole RHS implicitly, but the D(phi_0) is fixed at time 0.
 */

void micromixer::advanceHipsSingleStep_SemiImplicit() {

    if(domn->pram->Lsolver!="SEMI-IMPLICIT")
        return;

    //--------------- Set the explicit (mixing) terms

    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            domn->v.at(k)->getRhsMix(gf, dxc);

    //--------------- Perform the implicit integration on each cell

    for(int i=0; i<domn->ngrd; i++)
        cvode->integrateCell(i, dt);

    //---------------

  }

///////////////////////////////////////////////////////////////////////////////
/** Advance HiPS solution: diffusion and reaction; Some terms are implicit, others explicit.
 *  Nominally the mixing terms are explicit. Calling the cvode driver.
 *  First order.
 *  dphi/dt =  D(phi_0) + S(phi) : solving from t0 to t1.
 *  Here, D() is the diffusive term, and S() is the (stiff) source term.
 *  We solve the whole RHS implicitly, but the D(phi_0) is fixed at time 0.
 */

void micromixer::advanceHipsSingleStep_StrangSplit() {

    if(domn->pram->Lsolver!="STRANG")
        return;

    //--------------- First step: phi_1 = phi_0 + 0.5*dt*D(phi_0)

 
 
  

    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            domn->v.at(k)->getRhsMix(gf, dxc);

    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            for(int i=0; i < domn->ngrd; i++)
                domn->v.at(k)->d.at(i) = domn->v.at(k)->d.at(i) + 0.5*dt*domn->v.at(k)->rhsMix.at(i);



    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            domn->v.at(k)->getRhsSrc();

    for(int i=0; i<domn->ngrd; i++)
        cvode->integrateCell(i, dt);



  
    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            domn->v.at(k)->getRhsMix(gf, dxc);

    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            for(int i=0; i < domn->ngrd; i++)
                domn->v.at(k)->d.at(i) = domn->v.at(k)->d.at(i) + 0.5*dt*domn->v.at(k)->rhsMix.at(i);


    //-------------------------

}
///////////////////////////////////////////////////////////////////////////////////



