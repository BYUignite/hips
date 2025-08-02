
/// \file hips.cc
/// \brief Implementation of the HiPS (Hierarchical Parcel Swapping) model.
///
/// This file contains the core implementation of the HiPS model, used for 
/// simulating turbulent mixing. It supports optional functionality for chemical 
/// reactions, which can be enabled by defining the `REACTIONS_ENABLED` macro.
///
/// Dependencies:
/// - YAML-CPP: Used for reading and parsing configuration files.
/// - Batch reactor classes (`batchReactor_cvode.h`, `batchReactor_cantera.h`): 
///   Included only when `REACTIONS_ENABLED` is defined, enabling reaction modeling.
/// - Standard C++ libraries: Utilized for input/output, mathematical computations, 
///   and data handling.
///
/// \note Define `REACTIONS_ENABLED` at compile time to enable chemical reaction functionality.
////////////////////////////////////////////////////////////////////////////////////////////////

#include "hips.h"

#ifdef REACTIONS_ENABLED
#include "batchReactor_cvode.h"
#include "batchReactor_cantera.h"
#endif

#include "randomGenerator.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
hips::hips(int nLevels_, 
           double domainLength_, 
           double tau0_, 
           double C_param_, 
           int forceTurb_,
           int nVar_,
           vector<double> &ScHips_,
           bool performReaction_,
           shared_ptr<void> vcantSol,
           int seed,
           int realization_): 
    nLevels(nLevels_), 
    domainLength(domainLength_), 
    tau0(tau0_),
    C_param(C_param_), 
    forceTurb(forceTurb_),       
    ScHips(ScHips_),   
    nVar(nVar_),                       
    LrandSet(true),              
    rand(seed),
    performReaction(performReaction_),
    realization(realization_){

    #ifdef REACTIONS_ENABLED
    if(performReaction) {
        shared_ptr<Cantera::Solution> cantSol = static_pointer_cast<Cantera::Solution>(vcantSol);
        gas = cantSol->thermo(); 
        nsp = gas->nSpecies();

        bRxr = make_shared<batchReactor_cvode>(cantSol);                                // By default, use batchReactor_cvode

        // Uncomment the following line to switch to batchReactor_cantera
        // bRxr = make_unique<batchReactor_cantera>(cantSol);
    }
    #endif

    // Resize vectors to the number of variables
    varData.resize(nVar);
    varName.resize(nVar); 

    set_tree(nLevels, domainLength, tau0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

hips::hips(double C_param_, 
           int forceTurb_,
           int nVar_,
           vector<double> &ScHips_,
           bool performReaction_,
           shared_ptr<void> vcantSol,
           int seed,
           int realization_): 
    C_param(C_param_), 
    forceTurb(forceTurb_),       
    nVar(nVar_),                       
    ScHips(ScHips_),   
    LrandSet(true),              
    rand(seed),
    performReaction(performReaction_),
    realization(realization_){

    #ifdef REACTIONS_ENABLED
    // Initialize Cantera thermo phase and species count.
    if(performReaction) {
        shared_ptr<Cantera::Solution> cantSol = static_pointer_cast<Cantera::Solution>(vcantSol);
        gas = cantSol->thermo(); 
        nsp = gas->nSpecies();

        // Set up the default batch reactor (cvode).
       // bRxr = make_shared<batchReactor_cvode>(cantSol);

        // Uncomment the following line to switch to batchReactor_cantera.
         bRxr = make_shared<batchReactor_cantera>(cantSol);
    }
    #endif

    // Resize vectors to accommodate the number of variables.
    varData.resize(nVar);
    varName.resize(nVar);        
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void hips::set_tree(int nLevels_, double domainLength_, double tau0_){    
 
    nLevels= nLevels_; 
    domainLength = domainLength_; 
    tau0 = tau0_; 

    if (nLevels == -1)  
        nLevels = nL; 

    iEta = nLevels - 3;                                // Kolmogorov level; if nLevels = 7, then 0, 1, 2, 3, (4), 5, 6; iEta=4 is the lowest swap level: swap grandchildren of iEta=4 at level 6.
           
    int maxSc = 1.0;

    for (int i=0; i<ScHips.size(); i++)
        maxSc = ScHips[i]>maxSc ? ScHips[i] : maxSc;
    
    if (maxSc > 1.0)
        nLevels += ceil(log(maxSc)/log(4));            // Changing number of levels!
    
    Nm1 = nLevels - 1;
    Nm2 = nLevels - 2;
    Nm3 = nLevels - 3;
    
    // -------------------------- 
    
    nparcels = static_cast<int>(pow(2, Nm1));
    parcelTimes.resize(nparcels,0);
    i_batchelor.resize(nVar,0);
    
    vector<double> levelLengths(nLevels);              // Including all levels, but last 2 don't count:
    vector<double> levelTaus(nLevels);                 // Smallest scale is 2 levels up from bottom
    levelRates   = vector<double>(nLevels);

    for (int i=0; i<nLevels; i++) {
        levelLengths[i] = domainLength * pow(Afac,i);
       //levelLengths[i] = domainLength * pow(Anew,i);

        levelTaus[i] = tau0 * pow(levelLengths[i]/domainLength, 2.0/3.0) / C_param;
        levelRates[i] = 1.0/levelTaus[i] * pow(2.0,i);
    }

    LScHips = ScHips.size() > 0 ? true : false;
    if (LScHips) {                                     // Ccorrect levels for high Sc (levels > Kolmogorov)
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
            i_batchelor[k] = iEta + log(ScHips[k])/log(4);
        else
            i_batchelor[k] = iEta;

        i_plus[k] = ceil(i_batchelor[k]);
    }
    
    //------------------- Set the parcel addresses (index array)

    varRho.resize(nparcels);
    Temp.resize(nparcels);
 
    pLoc.resize(nparcels);
    for (int i=0; i<nparcels; i++)
        pLoc[i] = i;
} 

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void hips::set_tree(double Re_, double domainLength_, double tau0_, std::string approach_) {
    Re = Re_;
    domainLength = domainLength_;
    tau0 = tau0_;
    approach = approach_;

    double baseLevelEstimate = (3.0 / 4) * log(1 / Re) / log(Afac);                               // Calculate the base tree level estimate (non-integer)
    int baseLevel;

    if (approach == "rounding") {
        baseLevel = round(baseLevelEstimate);                                                     // Round the base level to the nearest integer
    } 
    else if (approach == "probability") {
        baseLevel = ceil(baseLevelEstimate);                                                      // Ceil the base level to the nearest integer
        int previousLevel = baseLevel - 1;
        Prob = baseLevelEstimate - previousLevel;                                                 // Calculate the probability
    } 
    else if (approach == "micromixing") {
        baseLevel = ceil(baseLevelEstimate);                                                      // Ceil the base level to the nearest integer
        lStar = std::pow(Re, -3.0 / 4);                                                           // Calculate lStar based on Re
    } 
    else if (approach == "dynamic_A") {
        baseLevel = round(baseLevelEstimate);                                                     // Round the base level to the nearest integer
        Anew = exp(-log(Re) / ((4.0 / 3.0) * baseLevel));                                         // Calculate the new value of parameter A
    } 
    else {
        throw std::invalid_argument("Invalid approach specified");                                // Handle invalid approach case if needed
    }

    nL = baseLevel + 3;                                                                           // Set the number of levels for the binary tree structure

    //----------------------------------------------------------------------
    nLevels = nL;
    iEta = nLevels - 3;  // Kolmogorov level 

    int maxSc = 1;
    for (const auto &sc : ScHips)
        maxSc = std::max(maxSc, static_cast<int>(sc));
    if (maxSc > 1.0)
        nLevels += ceil(log(maxSc) / log(4));                           

    Nm1 = nLevels - 1;
    Nm2 = nLevels - 2;
    Nm3 = nLevels - 3;

    nparcels = static_cast<int>(pow(2, Nm1));
    parcelTimes.resize(nparcels, 0);
    i_batchelor.resize(nVar, 0);

    std::vector<double> levelLengths(nLevels);                                // Including all levels, but last 2 don't count
    std::vector<double> levelTaus(nLevels);                                   // Smallest scale is 2 levels up from bottom
    levelRates.resize(nLevels);

    for (int i = 0; i < nLevels; ++i) {
        levelLengths[i] = domainLength * pow((approach == "dynamic_A" ? Anew : Afac), i);

        levelTaus[i] = tau0 * pow(levelLengths[i] / domainLength, 2.0 / 3.0) / C_param;
        levelRates[i] = 1.0 / levelTaus[i] * pow(2.0, i);
    }

    if (approach == "micromixing") {                                          // Adjust rates for micromixing model
        levelTaus[Nm3] = tau0 * pow(lStar / domainLength, 2.0 / 3.0) / C_param;
        levelRates[Nm3] = 1.0 / levelTaus[Nm3] * pow(2.0, Nm3);
    }

    if (approach == "probability") {                                          // Adjust final mixing rate based on probability
        levelRates[Nm3] = levelRates[nL - 3] * Prob;
    }

    LScHips = !ScHips.empty();                                               // Correct levels for high Sc (levels > Kolmogorov)
    if (LScHips) {
        for (int i = iEta + 1; i < nLevels; ++i) {
            levelTaus[i] = tau0 * pow(levelLengths[iEta] / domainLength, 2.0 / 3.0) / C_param;
            levelRates[i] = 1.0 / levelTaus[i] * pow(2.0, i);
        }
    }

    //-----------------------------------------------------

    eddyRate_total = 0.0;
    for (int i = 0; i <= Nm3; ++i)
        eddyRate_total += levelRates[i];

    eddyRate_inertial = 0.0;
    for (int i = 0; i <= iEta; ++i)
        eddyRate_inertial += levelRates[i];

    //-----------------------------------------------------

    i_plus.resize(nVar);
    for (int k = 0; k < nVar; ++k) {
        if (ScHips[k] < 1.0)
            i_batchelor[k] = iEta + 1.5 * log(ScHips[k]) / log(4);
        else if (ScHips[k] > 1.0)
            i_batchelor[k] = iEta + log(ScHips[k]) / log(4);
        else
            i_batchelor[k] = iEta;
        i_plus[k] = ceil(i_batchelor[k]);
    }

    //-----------------------------------------------------

    varRho.resize(nparcels);
    Temp.resize(nparcels);
    pLoc.resize(nparcels);
    for (int i = 0; i < nparcels; ++i)
        pLoc[i] = i;
}

////////////////////////////////////////////////////////////////////////////////////
/// \brief Assigns variables, their corresponding weights, and names to the parcels in the HiPS tree.
///
/// This function projects the provided variables onto the parcels within the HiPS tree structure, 
/// using their corresponding weights. It ensures that each parcel in the tree is associated with 
/// the correct variable and weight, along with a descriptive name for better identification.
///
/// \param v         Vector of variables to be assigned to the parcels in the HiPS tree.
/// \param w         Vector of weights corresponding to each variable or parcel.
/// \param varN      String representing the name of the variable being assigned.
///
/// \note The size of `v` and `w` must match to ensure a one-to-one correspondence between 
///       variables and their weights. The function does not perform size validation internally.
///
/// \warning Ensure that the weights in `w` are normalized or appropriately scaled, as they 
///          directly influence the projection and subsequent simulations.
////////////////////////////////////////////////////////////////////////////////////

void hips::set_varData(std::vector<double> &v, std::vector<double> &w, const std::string &varN) {  
    
    varData[currentIndex] = std::make_shared<std::vector<double>>(projection(v, w));
    varName[currentIndex] = varN;

    currentIndex++; 
}

//////////////////////////////////////////////////////////////////////////////////// 
/// \brief Assigns variables, weights, names, and densities to the parcels in the HiPS tree.
///
/// This overloaded function assigns the specified variables, along with their associated weights, 
/// names, and densities, to the parcels within the HiPS tree structure. The function incorporates 
/// particle density during the projection process, ensuring a more accurate representation of 
/// parcel properties in the simulation.
///
/// \param v         Vector of variables to be assigned to the HiPS tree.
/// \param w         Vector of weights corresponding to each variable or parcel.
/// \param varN      String representing the name of the variable being assigned.
/// \param rho       Vector of densities corresponding to each flow particle.
///
/// \note This function is specifically overloaded to account for particle density, which enhances 
///       the accuracy of the parcel projection. Ensure that the size of `v`, `w`, and `rho` are consistent.
///
/// \warning The values in `rho` should be physically meaningful and consistent with the simulation's 
///          requirements. Improper density values may lead to inaccuracies or instabilities in the simulation.
////////////////////////////////////////////////////////////////////////////////////

void hips::set_varData(std::vector<double> &v, std::vector<double> &w, const std::string &varN, const std::vector<double> &rho) {
    
    //old std::pair<std::vector<double>, std::vector<double>> results = projection(v, w, rho);
    //old std::vector<double> vh = results.first;
    //old std::vector<double> rho_h = results.second;
    //old
    //old varData[currentIndex] = std::make_shared<std::vector<double>>(projection(v, w));            
    //old varRho = std::vector<double>(rho_h);
    //old varName[currentIndex] = varN;
    //old
    //old currentIndex++; 

    std::pair<std::vector<double>, std::vector<double>> results = projection(v, w, rho);

    varData[currentIndex] = std::make_shared<std::vector<double>>(results.first);
    varRho = results.second;
    varName[currentIndex] = varN;
 
    currentIndex++; 
}

//////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Projects values from flow particles onto HiPS parcels assuming constant density.
///
/// This function maps the values of flow particles onto HiPS parcels under the assumption of constant density.
/// The projection is performed using the following equation:
/// \f[
/// \sum_{i=0}^{\text{Number of Flow Particles (FP)}} (\phi_{\text{FP}} \, \mathrm{d}x_{\text{FP}})_{i} = 
/// \sum_{j=0}^{\text{Number of HiPS Parcels (HP)}} (\phi_{\text{HP}} \, \mathrm{d}x_{\text{HP}})_{j}
/// \f]
/// This ensures conservation of properties such as mass or concentration during the projection.
///
/// \param vcfd          Vector of variables from flow particles to be mapped to HiPS parcels.
/// \param weight        Vector of weights, with one weight assigned to each flow particle.
/// \return              Vector of projected values for HiPS parcels.
///
/// \note The function assumes constant density throughout the domain. For cases with varying density, 
///       use an appropriate overloaded function or method.
/// 
/// \warning Ensure that the `vcfd` and `weight` vectors have matching sizes, as any discrepancy 
///          may result in undefined behavior or incorrect projections.
///////////////////////////////////////////////////////////////////////////////

std::vector<double> hips::projection(std::vector<double> &vcfd, std::vector<double> &weight) {
    
    xc = setGridCfd(weight);                               // Populate the physical domain for flow particles
    xh = setGridHips(nparcels);                            // Populate the physical domain for hips parcels

    int nc = xc.size() - 1;                 
    int nh = xh.size() - 1; 

    std::vector<double> vh(nh, 0.0);
    int jprev = 0;

    for(int i = 0; i < nh; i++) {
        for(int j = jprev + 1; j <= nc; ++j) {
            if(xc[j] <= xh[i + 1]) {
                double d1 = xc[j] - xc[j - 1];
                double d2 = xc[j] - xh[i];
                // calculation of shortest distance
                // handling if distance is zero but substarction gives comutational error
                double d  = std::min(d1, d2) < 1E-15 ? 0 : std::min(d1, d2);

                vh[i] += vcfd[j - 1] * d;
            } 
            else {
                double d1 = xh[i + 1] - xc[j - 1];
                double d2 = xh[i + 1] - xh[i];
                // calculation of shortest distance
                // handling if distance is zero but substarction gives comutational error
                double d  = std::min(d1, d2) < 1E-15 ? 0 : std::min(d1, d2);

                vh[i] += vcfd[j - 1] * d;
                jprev = j - 1;
                break;
            }
        }
        vh[i] /= (xh[i + 1] - xh[i]);   
    }
    return vh;
}

///////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Projects the values from flow particles onto HiPS parcels, accounting for particle density.
///
/// This function calculates the projection of values from flow particles onto HiPS parcels, 
/// incorporating the density of each particle into the computation. The resulting projection 
/// ensures conservation of both the property values and the density. The function returns a 
/// pair of vectors: one representing the projected values on the HiPS parcels and another for 
/// the densities.
///
/// \param vcfd         Vector of variables from flow particles to be projected onto HiPS parcels.
/// \param weight       Vector of weights corresponding to each flow particle.
/// \param density      Vector of densities associated with each flow particle.
///
/// \return A pair of vectors:
///         - First vector: Projected values on the HiPS parcels.
///         - Second vector: Densities on the HiPS parcels.
///
/// \note This is an overloaded version of the projection function, designed to include particle density 
///       in the computation. Ensure that the `vcfd`, `weight`, and `density` vectors are of equal size 
///       for consistency.
///
/// \warning Discrepancies in the size of the input vectors (`vcfd`, `weight`, and `density`) may lead 
///          to undefined behavior or inaccurate projections. Verify the inputs before calling this function.
///////////////////////////////////////////////////////////////////////////////////////////////////

std::pair<std::vector<double>, std::vector<double>> hips::projection(std::vector<double> &vcfd, 
                                                                     std::vector<double> &weight, 
                                                                     const std::vector<double> &density) {
    
    xc = setGridCfd(weight);
    xh = setGridHips(nparcels);

    // Initialize variables
    int nc = xc.size() - 1;
    int nh = xh.size() - 1;

    std::vector<double> vh(nh, 0.0);
    std::vector<double> rho_h(nh, 0.0);
    int jprev = 0;

    for (int i = 0; i < nh; i++) {

        for (int j = jprev + 1; j <= nc; ++j) {
            if (xc[j] <= xh[i + 1]) {
                double d = std::min(xc[j] - xc[j - 1], xc[j] - xh[i]);
                rho_h[i] += density[j - 1] * d;
                vh[i] += density[j - 1] * vcfd[j - 1] * d;
            } 
            else {
                double d = std::min(xh[i + 1] - xc[j - 1], xh[i + 1] - xh[i]);
                rho_h[i] += density[j - 1] * d;
                vh[i] += density[j - 1] * vcfd[j - 1] * d;
                jprev = j - 1;
                break;
            }
        }

        // ------------------------ Normalize results

        rho_h[i] /= (xh[i + 1] - xh[i]);
        vh[i] /= rho_h[i] * (xh[i + 1] - xh[i]);
    }
    return {vh, rho_h};
}

/////////////////////////////////////////////////////////////////////////////////
/// \brief Generates a physical domain for flow particles based on their weights.
///
/// This function creates a grid of positions for flow particles, where each particle occupies a portion 
/// of the domain proportional to its weight. The total length of the domain is assumed to be 1, and the 
/// sum of all portions equals 1. The resulting vector represents the positions of particles along the domain.
///
/// \param w         Vector of weights, where each weight determines the portion of the domain occupied 
///                  by a particle.
/// \return          A vector of grid positions for the flow particles.
///
/// \note The function assumes that the weights in `w` are normalized or properly scaled such that the 
///       total sum matches the domain length of 1. If the weights are not normalized, the resulting grid 
///       may not represent a valid physical domain.
///
/// \warning Ensure that the input weight vector `w` is non-empty and contains positive values. Zero or 
///          negative weights may lead to undefined behavior or invalid domain generation.
/////////////////////////////////////////////////////////////////////////////////

std::vector<double> hips::setGridCfd(std::vector<double> &w) {

    double sumw = 0.0;
    for(int i=0; i<w.size(); i++)
        sumw += w[i];
   
    std::vector<double> pos;                               // Initializing a vector to hold the grid positions
    double posL = 0.0;                                     // Initializing the starting position

    int i = 0;

    while (i <= w.size()) {                               // Generate the grid positions based on the weights
        pos.push_back(posL);                              // Add the current position to the grid
        posL += w[i]/sumw;                                // Move to the next position by adding the corresponding weight
        i++;                                              
    }
    return pos;                                           // Return the generated grid positions
}

///////////////////////////////////////////////////////////////////////////////
/// \brief Generates a physical domain for HiPS parcels.
///
/// This function creates a grid of positions for HiPS parcels, where each parcel occupies 
/// an equal portion of the physical domain. The total size of the domain corresponds to the 
/// size specified in the `setGridCfd()` function, ensuring consistency between the HiPS 
/// and flow particle domains.
///
/// \param N         The number of grid points for the HiPS parcels.
/// \return          A vector representing the grid positions for the HiPS parcels.
///
/// \note The function assumes that the physical domain is evenly divided among the parcels. 
///       Ensure that the number of grid points (`N`) is consistent with the physical domain size 
///       defined in the simulation setup.
///
/// \warning If `N` is less than or equal to zero, the function may produce an empty or invalid grid. 
///          Validate the input to avoid unexpected behavior.
///////////////////////////////////////////////////////////////////////////////

std::vector<double> hips::setGridHips(int N){

    std::vector<double> xh(N + 1);                               // Initialize a vector to hold the grid points
    double step = 1.0 / N;                                       // Calculate the step size

    for(int i = 0; i <= N; i++)                                 // Populate the grid with evenly spaced points
        xh[i] = i * step;
        
    return xh;                                                  // Return the generated grid
}

///////////////////////////////////////////////////////////////////////////////////
/// \brief Runs the HiPS simulation, advancing the solution using eddy events.
///
/// This function performs the core HiPS loop: sampling eddy events, performing 
/// subtree swaps, advancing parcels, and optionally triggering reactions.
/// It runs until the specified simulation time (`tRun`) is reached and writes 
/// data periodically based on either eddy count or elapsed simulation time.
///
/// ### Key operations:
/// - Samples the next eddy event time (`dtEE`)
/// - Selects and swaps subtrees at a given level
/// - Applies micromixing and reactions (if enabled)
/// - Writes output data either:
///     - Every `outputIntervalEddy` eddy events (if enabled), or
///     - Every `outputIntervalTime` seconds (if enabled)
/// - At the end of the simulation, calls `saveAllParameters()` to store 
///   input and configuration data in `../post/parameters.dat`.
///
/// \param tRun              Total simulation run time (in seconds).
/// \param shouldWriteData   Flag to enable/disable periodic data writing.
///
/// \note To control output frequency, use:
///       - `setOutputIntervalEddy(int interval)`
///       - `setOutputIntervalTime(double interval)`
///       - If neither is called, the default behavior is writing every 1000 eddy events.
///
/// \note Output files are saved using `writeData(realization, ...)` and include the realization index.
///
/// \note At the end of the run, `saveAllParameters()` is automatically called 
///       to document simulation settings.
///
/// \warning Long simulations may generate many output files. Adjust output intervals or 
///          disable writing (`shouldWriteData = false`) to manage storage needs.
///
/// \see hips::writeData(), hips::saveAllParameters()
///////////////////////////////////////////////////////////////////////////////////

void hips::calculateSolution(const double tRun, bool shouldWriteData) {
    
    unsigned long long nEddies = 0;                               // Number of eddy events
    int fileCounter = 0;                                          // Number of data files written
    int iLevel;                                                   // Tree level of EE with top at iLevel=0
    int iTree;                                                    // One of two subtrees involved in swap at iLevel                                           
    time = 0.0;                                                   // Initialize simulation time
    int lastEddyOutput = 0;                                       // Track last eddy-based output event

    // Apply default values if user hasn't set them
    if (!useEddyBasedWriting && !useTimeBasedWriting) {
        outputIntervalEddy = DEFAULT_EDDY_INTERVAL;
        useEddyBasedWriting = true;  // Default to eddy-based writing
    }

    sample_hips_eddy(dtEE, iLevel);                               // Get first EE at time 0+dtEE
    nEddies++;
    eddyCounter = 0;                                              // Reset eddy counter at start
    lastOutputTime = 0.0;                                         // Reset last output time

    while (time + dtEE <= tRun) {
        time += dtEE;
        selectAndSwapTwoSubtrees(iLevel, iTree);
        advanceHips(iLevel, iTree);                              // Reaction and micromixing (if needed) to t=time

        sample_hips_eddy(dtEE, iLevel);
        nEddies++;
        eddyCounter++;  

        //  Only check the selected mode (set in example code or by default)
        bool writeByEddy = (useEddyBasedWriting && eddyCounter >= lastEddyOutput + outputIntervalEddy);
        bool writeByTime = (useTimeBasedWriting && time - lastOutputTime >= outputIntervalTime);

        if (shouldWriteData) {
            if (writeByEddy) {  //  Only write if using eddy-based writing
                writeData(realization, ++fileCounter, time);
                lastEddyOutput = eddyCounter;  //  Update last output event
            }
            else if (writeByTime) {  // Only write if using time-based writing
                 writeData(realization, ++fileCounter, time);
                lastOutputTime = time;
            }
        }
    }

    // Ensure the final time step completes
    time = tRun;
    iLevel = 0; 
    iTree = 0;

    if (performReaction)
        reactParcels_LevelTree(iLevel, iTree);                   // React all parcels up to end time
    saveAllParameters();
}

///////////////////////////////////////////////////////////////////////////////
/// \brief Samples stochastic eddy events on the HiPS tree, determining the time increment and tree level.
///
/// This function performs stochastic sampling to determine when (\f$\Delta t_{EE}\f$) and at what level (\f$i_{Level}\f$)
/// in the HiPS tree the next eddy event will occur. The time to the next eddy event is sampled based on the total eddy rate.
/// The tree level of the event is chosen depending on whether it occurs in the inertial or Batchelor region of turbulence.
///
/// \param dtEE         Time increment to the next eddy event (\f$\Delta t_{EE}\f$), sampled stochastically.
/// \param iLevel       Tree level (\f$i_{Level}\f$) at which the eddy event occurs, determined probabilistically.
///
/// \note The function distinguishes between events in the inertial and Batchelor regions based on turbulence properties.
///       Ensure that the HiPS tree is correctly initialized before calling this function.
///
/// \warning The stochastic nature of this function requires a properly seeded random generator to ensure reproducibility 
///          in simulations where determinism is necessary.
////////////////////////////////////////////////////////////////////////////////

void hips::sample_hips_eddy(double &dtEE, int &iLevel) {

    static double c1 = 1.0 - pow(2.0, 5.0/3.0*(iEta+1));
    static double c2 = pow(2.0, Nm2) - pow(2.0, iEta+1);
    static double c3 = pow(2.0, iEta+1);

    //--------------- time to next eddy

    double r = rand.getRand();
    dtEE = -log(r)/eddyRate_total;

    //----------------- get eddy level

    r = rand.getRand();

    if ( r <= eddyRate_inertial/eddyRate_total) {     // Inertial region
        r = rand.getRand();
        iLevel = ceil(3.0/5.0*log2(1.0-r*c1) - 1.0);
        if (iLevel < 0)    iLevel = 0;
        if (iLevel > iEta) iLevel = iEta;
    }

    else {                                            // "Batchelor" region
        r = rand.getRand();
        iLevel = ceil(log2(r*c2 + c3) - 1.0);
        if (iLevel < iEta+1) iLevel = iEta+1;
        if (iLevel > Nm3) iLevel = Nm3;
    }
    return;
}

///////////////////////////////////////////////////////////////////////////////
/// \brief Performs eddy events by swapping parcels within the HiPS tree.
///
/// This function executes parcel swaps by randomly selecting nodes at the specified level of the HiPS tree.
/// It identifies the starting indices of subtrees to be swapped, calculates the number of parcels to swap, 
/// and performs the swap operation efficiently using bitwise operations. The process mimics the hierarchical 
/// structure of turbulent mixing.
///
/// \param iLevel         Input level of the tree where the base of the swap occurs.
/// \param iTree          Output parameter indicating which subtree at the given level is selected for swapping.
///
/// The process is as follows:
/// - Randomly select a node on `iLevel`.
/// - Traverse two levels down to identify subtrees `0q` and `1r`, where `q` and `r` are random binary values (0 or 1).
/// - Determine the starting indices of the `Q-tree` and `R-tree` for the swap and compute the number of parcels.
/// - Swap the corresponding parcels between the subtrees.
///
/// ### Example for a 6-level tree:
/// - Tree levels: 0, 1, 2, 3, 4, 5.
/// - If `iLevel = 1`:
///   - Suppose the selected node is `i = 01`.
///   - Subtrees for swapping are `0q = 00` and `1r = 11`.
///   - Swapping involves parcels `0100**` with `0111**`, or `(01|00|**)` with `(01|11|**)`.
///   - In binary terms, the swap is equivalent to exchanging `i0qs` with `i1rs`, where:
///     - `i = 01`
///     - `0q = 00`
///     - `1r = 11`
///     - `s = **` (remaining bits).
///
/// ### Implementation:
/// - Bitwise operations are used for efficient calculations of powers of 2.
/// - The swap operation is performed by flipping the bits for `0q` and `1r`, which effectively swaps the subtrees.
///
/// ### Visual Representation of Tree:
/// ```
/// Level 0          * (root)
///                /     \
/// Level 1      *       (*)
///              / \       / \
/// Level 2    *   *     *   *
/// Level 3   *   [*]    *   [*]
/// Level 4  * *   * *   * *   * *
/// Level 5 00 01  02 03  ...  30 31
/// ```
/// - Subtrees `0q` and `1r` correspond to specific branches of the tree.
/// - Swapping occurs within highlighted sections of Level 3, identified by bit manipulation.
///
/// \warning Ensure that the input `iLevel` is within the valid range of tree levels and that the tree is 
///          properly initialized before invoking this function.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void hips::selectAndSwapTwoSubtrees(const int iLevel, int &iTree) {

    iTree = rand.getRandInt((1 << iLevel)-1);
    int zero_q = rand.getRandInt(1);                                    // 0q where q is 0 or 1
    int one_r  = 2 + rand.getRandInt(1);                                // 1r where r is 0 or 1

    int Qstart = (zero_q << (Nm3-iLevel)) + (iTree << (Nm1-iLevel));     // starting index of Q parcels
    int Rstart = (one_r  << (Nm3-iLevel)) + (iTree << (Nm1-iLevel));     // starting index of R parcels
    int nPswap = 1 << (Nm3-iLevel);                                      // number of parcels that will be swapped

    int Qend = Qstart + nPswap;                                          // inclusive indices are Qstart to Qend-1
    int Rend = Rstart + nPswap;                                          // inclusive indices are Rstart to Rend-1
    vector<int> aa(pLoc.begin()+Qstart, pLoc.begin()+Qend);
    copy(pLoc.begin()+Rstart, pLoc.begin()+Rend, pLoc.begin()+Qstart);   // python: pLoc[Qstart:Qend]=pLoc[Rstart:Rend]
    copy(aa.begin(), aa.end(), pLoc.begin()+Rstart);                     // python: pLoc[Rstart:Rend]=aa
}

/////////////////////////////////////////////////////////////////////////////////
/// \brief Advances the HiPS model by simulating micromixing and reactions at a specific tree level.
///
/// This function models the interaction of parcels within the HiPS tree structure at a given level, 
/// simulating the effects of micromixing and reactions. The process involves determining whether 
/// conditions for micromixing are met and, if so, triggering reactions and mixing for parcels under 
/// the specified node. The operation depends on turbulence forcing and micromixing thresholds.
///
/// \param iLevel    The level of the HiPS tree at which the eddy event occurs.
/// \param iTree     The root node of the eddy event at the specified level.
///
/// ### Key Operations:
/// - Identifies parcels within the specified tree level.
/// - Checks whether conditions for micromixing are satisfied.
/// - Applies turbulence forcing if enabled.
/// - Triggers reactions for the first variable that meets micromixing conditions.
/// - Applies mixing to subsequent variables without additional reactions.
///
/// \note Reactions are only performed for the first variable that satisfies the micromixing conditions. 
///       For all subsequent variables, only mixing operations are performed.
///
/// \warning Ensure that the HiPS tree structure is properly initialized and that `iLevel` and `iTree` 
///          correspond to valid levels and nodes within the tree to prevent undefined behavior.
//////////////////////////////////////////////////////////////////////////////////

void hips::advanceHips(const int iLevel, const int iTree) {

    if (forceTurb == 2 && iLevel == 0) {
        forceProfile();                                                  // Forcing for statistically stationary
    }

    bool rxnDone = false;                                               // React all variables once
    for (int k = 0; k < nVar; k++) {                                    // Upon finding first variable needing micromixing
        // Combined condition check with approach condition
        if ((iLevel >= i_plus[k]) || 
            (iLevel == i_plus[k] - 1 && rand.getRand() <= i_plus[k] - i_batchelor[k])) {
                if (!rxnDone && performReaction) {
                    reactParcels_LevelTree(iLevel, iTree);
                    rxnDone = true;
                }
                mixAcrossLevelTree(k, iLevel, iTree);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
/// \brief Retrieves the index of a variable by its name in the `varName` list.
///
/// This function searches the `varName` list for the specified variable name and determines its index.
/// It uses `std::find` and `std::distance` to locate the variable efficiently. If the variable name 
/// is not found, the function throws a `std::runtime_error`. This functionality is crucial for 
/// referencing variables dynamically within the HiPS model.
///
/// \param varName The name of the variable to search for (e.g., "enthalpy").
/// \return int The index of the variable in the `varName` list.
/// \throws std::runtime_error If the variable name is not found in the `varName` list.
///
/// ### Usage Example:
/// ```cpp
/// // Assume varName list is populated: {"temperature", "enthalpy", "density"}
/// int index = hips.get_varIndex("enthalpy");  // Returns 1
/// ```
///
/// \note This function is case-sensitive and assumes that the `varName` list has been populated,
///       typically via the `set_varData` method. Ensure that the variable names match exactly,
///       including case, to avoid errors.
///
/// \warning If the `varName` list is empty or not properly populated, the function may throw
///          unexpected errors. Verify the list contents before invoking this method.
///
/// \see set_varData
///////////////////////////////////////////////////////////////////////////////

int hips::getVariableIndex(const std::string &varName) const {

    auto it = std::find(this->varName.begin(), this->varName.end(), varName);
    if (it == this->varName.end()) {
        throw std::runtime_error("Error: Variable name '" + varName + "' not found.");
    }
    return std::distance(this->varName.begin(), it);
}

///////////////////////////////////////////////////////////////////////////////
/// \brief Simulates chemical reactions for parcels affected by a micromixing event.
///
/// This function performs chemical reactions for parcels involved in a micromixing process at a specific 
/// level and tree node within the HiPS structure. The reaction times are determined based on the last 
/// reaction time stored in `parcelTimes`. It dynamically retrieves the indices of relevant variables 
/// such as `enthalpy` and species mass fractions for accurate state updates.
///
/// \param iLevel         The tree level where the eddy event occurred.
/// \param iTree          The root node of the eddy event at the specified level.
///
/// \details
/// - `getVariableIndex("enthalpy")` is used to retrieve the index for `enthalpy`.
/// - Indices for species are dynamically retrieved using `gas->speciesName(i)`.
/// - The function updates parcel states, including `enthalpy` and species mass fractions, based on the reactions.
///
/// \throws std::runtime_error If a variable name is not found in the `varName` list.
///
/// \note 
/// - The `varName` list must be populated using `set_varData` before invoking this function.
/// - Reaction functionality is only available if `REACTIONS_ENABLED` is defined during compilation.
/// - This function relies on the HiPS model's proper initialization and an accurate setup of `parcelTimes`.
///
/// \warning Ensure that the necessary reaction data and variable names are correctly configured.
///          Missing or incorrectly configured `varName` entries may lead to runtime errors.
/// \see set_varData, getVariableIndex
///////////////////////////////////////////////////////////////////////////////

void hips::reactParcels_LevelTree(const int iLevel, const int iTree) {

    #ifdef REACTIONS_ENABLED

    int enthalpyIdx = getVariableIndex("enthalpy");  // Dynamically find enthalpy index
    int nP = 1 << (Nm1 - iLevel);
    int istart = iTree * nP;
    int iend = istart + nP;
    int ime;
    double dt;
    double h;
    std::vector<double> y(nsp);

    for (int i = istart; i < iend; i++) {
        ime = pLoc[i];
        dt = time - parcelTimes[ime];

        // Access enthalpy using its index
        h = (*varData[enthalpyIdx])[ime];

          // Access species using their indices
        for (int k = 0; k < nsp; k++) {
            int speciesIdx = getVariableIndex(gas->speciesName(k));
            y[k] = (*varData[speciesIdx])[ime];
        }
        if (performReaction) {
            bRxr->react(h, y, dt);
            varRho[ime] = bRxr->getDensity();
            Temp[ime] = bRxr->temperature;
        }
        // Update enthalpy and species
        (*varData[enthalpyIdx])[ime] = h;
        for (int k = 0; k < nsp; k++) {
            int speciesIdx = getVariableIndex(gas->speciesName(k));
            (*varData[speciesIdx])[ime] = y[k];
        }

        parcelTimes[ime] = time;
    }
    #endif
}

////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Uniformly mixes parcels at a specified level and subtree within the HiPS model.
///
/// This function performs uniform mixing of parcels based on their average values at a specified 
/// level and subtree of the HiPS tree. It is particularly suited for low Schmidt number (\f$Sc\f$) 
/// variables and supports mixing at levels above the lowest in the tree. Parcels are mixed in pairs 
/// or groups depending on the level and subtree index.
///
/// \param kVar   Index of the variable to be mixed (typically a transported variable determined by the caller).
/// \param iLevel The tree level whose grandchildren will be mixed.
/// \param iTree  The subtree at the given level where mixing will occur.
///
/// ### Process Overview:
/// - At each level, the function mixes parcels in pairs or larger groups as determined by the tree structure.
/// - The level (\p iLevel) defines the size of the groups to mix:
///   - Higher levels mix larger groups (e.g., groups of 4 parcels at Level 1).
///   - Lower levels mix smaller groups (e.g., pairs of parcels at Level 2).
/// - Mixing is performed uniformly by averaging the variable values (\p kVar) of parcels.
///
/// ### Example for a 5-Level Tree (Levels: 0 to 4):
/// - \p iLevel = 2:
///   - If \p iTree = 0, parcels (0,1) and (2,3) will be mixed.
///   - If \p iTree = 1, parcels (4,5) and (6,7) will be mixed.
///   - If \p iTree = 2, parcels (8,9) and (10,11) will be mixed.
///   - If \p iTree = 3, parcels (12,13) and (14,15) will be mixed.
/// - \p iLevel = 1:
///   - If \p iTree = 0, parcels (0,1,2,3) and (4,5,6,7) will be mixed.
///   - If \p iTree = 1, parcels (8,9,10,11) and (12,13,14,15) will be mixed.
///
/// ### Notes:
/// - **Density Assumption**: The function assumes all parcels have the same density. Be cautious when mixing scalars like mass fractions, as the current implementation directly mixes \f$Y_i\f$.
/// - **Index Calculation**: The function uses bitwise left shift (\p <<) to calculate the starting and ending indices for parcel mixing.
///
/// \warning Ensure that \p iLevel and \p iTree correspond to valid levels and subtrees in the HiPS structure to prevent undefined behavior.
////////////////////////////////////////////////////////////////////////////////////////////

void hips::mixAcrossLevelTree(int kVar, const int iLevel, const int iTree) {
    
    int istart;
    int iend;

    int nPmix = 1 << (nLevels - iLevel - 2);                   // Number of parcels mixed together

    int ime;

    //---------- Mix left branch of iTree

    istart = iTree << (Nm1-iLevel);  
    iend = istart + nPmix;

    double s = 0;                                               // Initialize sum to 0
    for (int i=istart; i<iend; i++) {
        ime = pLoc[i];
        s += (*varData[kVar])[ime];
    }
    for (int i=istart; i<iend; i++) {
        ime = pLoc[i];
        (*varData[kVar])[ime] = s / nPmix; 
    }

    //--------- Mix right branch of iTree

    istart = iend;
    iend = istart + nPmix;

    s = 0;                   // initialize sum to 0
    for (int i=istart; i<iend; i++) {
        ime = pLoc[i];
        s += (*varData[kVar])[ime];
    }
    for (int i=istart; i<iend; i++) {
        ime = pLoc[i];
        (*varData[kVar])[ime] = s / nPmix; 
    }
}

///////////////////////////////////////////////////////////////////////////
/// \brief Adjusts the HiPS profile to enforce statistical stationarity.
///
/// This function modifies the parcel values in the HiPS profile to achieve statistical stationarity. 
/// Specifically, it adjusts the average value of parcels in the left half of the profile to 0 and the 
/// average value in the right half to 1. This is particularly useful for simple scalar variables, such 
/// as a mixture fraction ranging between 0 and 1.
///
/// ### Process Overview:
/// - For each variable in the HiPS profile:
///   1. Compute the average value of parcels in the left half.
///   2. Compute the average value of parcels in the right half.
///   3. Adjust the parcel values to enforce the desired averages (0 for the left half, 1 for the right half).
///
/// \note 
/// - This function directly modifies the parcel data in the HiPS profile.
/// - It assumes that variable values are bounded and appropriate for normalization (e.g., scalars between 0 and 1).
///
/// \warning 
/// - Ensure that the HiPS profile is properly initialized and contains valid data before calling this function.
/// - Unbounded or inappropriate variable values may lead to unexpected results.
///////////////////////////////////////////////////////////////////////////
 
void hips::forceProfile() {
    // Loop through each variable in the HiPS profile
    for (int k = 0; k < varData.size(); k++) {
        double s;                                                    // Temporary variable for summation

        //---------- Force the left half of parcels to average 0 ----------

        for (int i = 0; i < nparcels >> 1; i++)
            s += (*varData[k])[pLoc[i]];                             // Calculate the sum of values in the left half of parcels
        
        s /= (nparcels >> 1); // Calculate the average of values in the left half of parcels
        
        for (int i = 0; i < nparcels >> 1; i++)
            (*varData[k])[pLoc[i]] += (-s - 0.0);                    // Adjust values in the left half of parcels to achieve an average of 0

        //---------- Force the right half of parcels to average 1 ----------
        s = 0.0;

        for (int i = nparcels >> 1; i < nparcels; i++)
            s += (*varData[k])[pLoc[i]];                             // Calculate the sum of values in the right half of parcels
        
        s /= (nparcels >> 1);                                        // Calculate the average of values in the right half of parcels
        
        for (int i = nparcels >> 1; i < nparcels; i++)
            (*varData[k])[pLoc[i]] += (-s + 1.0);                    // Adjust values in the right half of parcels to achieve an average of 1
    }
}

///////////////////////////////////////////////////////////////////////////////////
/// \brief Writes simulation data to a file for a specific realization, time, and file index.
///
/// This function saves the current state of simulation data to a file, organized by 
/// realization and time step. The output is stored in a subdirectory corresponding to 
/// the realization index (e.g., `rlz_00000`, `rlz_00001`), and the file is named sequentially 
/// (e.g., `Data_00001.dat`, `Data_00002.dat`) to maintain an ordered record of outputs. 
/// The simulation time is also written into the file for reference.
///
/// ### Key Operations:
/// 1. Creates a subdirectory named `rlz_XXXXX` based on the realization index.
/// 2. Constructs the output filename using the sequential file index.
/// 3. Writes the simulation variables (e.g., temperature, species mass fractions) in 
///    scientific format with high precision.
///
/// \param real        Realization index used to name the subdirectory (`rlz_XXXXX`).
/// \param ifile       Sequential index for naming the output file within the realization.
/// \param outputTime  Simulation time associated with the data, included in the file header.
///
/// \note 
/// - The function ensures the output directory for the specified realization is created if it does not exist.
/// - It writes all data with high precision for accurate post-processing.
///
/// \warning 
/// - If the function fails to create the directory or open the output file, it may throw 
///   an error or silently fail. Ensure file system permissions and disk space are sufficient.
///////////////////////////////////////////////////////////////////////////////////

void hips::writeData(int real, const int ifile, const double outputTime) {

    stringstream ss1, ss2;
    string s1, s2;

    // Prepare directory path
    ss1 << "../data/rlz_" << setfill('0') << setw(5) << real;
    ss1 >> s1;

    // Create directories if they don't exist
    #ifdef _WIN32
        system(("mkdir " + s1).c_str());
    #else
        system(("mkdir -p " + s1).c_str());
    #endif

    ss2 << "Data_" << setfill('0') << setw(5) << ifile << ".dat";
    ss2 >> s2;

    string fname = s1 + "/" + s2;
    ofstream ofile(fname.c_str());
    cout << endl << "writing data for time " << outputTime << " to file: " << fname.c_str();

    // Check if file opened successfully
    if (!ofile) {
        cerr << "Error: Unable to open file " << fname << " for writing!" << endl;
        return;
    }

    // Write metadata (header information)
    ofile << "# time = " << outputTime << "\n";
    ofile << "# Grid Points = " << nparcels << "\n";
        
    // Write column names (include temperature if reactions are enabled)
    if(performReaction)
        ofile << setw(19) << "# Temp";  // Include temperature column if reactions are enabled

    for (const auto& varN : varName) {
        ofile << setw(19) << "# " << varN;
    }
    ofile << endl;  // End of the header line

    // Set scientific notation and precision
    ofile << scientific;
    ofile << setprecision(10);

    // Write data
    for (int i = 0; i < nparcels; i++) {
        // Write temperature first if reactions are enabled
        if(performReaction)
            ofile << setw(19) << Temp[pLoc[i]];

        // Write variables
        for (int k = 0; k < nVar; k++) {
            ofile << setw(19) << (*varData[k])[pLoc[i]];
        }
        ofile << endl;
    }

    ofile.close();
   // cout << "Data successfully written to: " << fname << endl;
}
   
///////////////////////////////////////////////////////////////////////////////
/// \brief Projects HiPS parcel values back onto the flow particles.
///
/// This function reverses the projection process, redistributing the values stored in the HiPS 
/// parcels back to the flow particles. It ensures conservation of properties such as mass or 
/// concentration by maintaining consistency between the HiPS parcels and flow particles.
///
/// ### Conservation Principle:
/// The projection follows the equation:
/// \f[
/// \sum_{j=0}^{\text{Number of HP}} (\phi_{\text{HP}} \, \mathrm{d}x_{\text{HP}})_{j} = 
/// \sum_{i=0}^{\text{Number of FP}} (\phi_{\text{FP}} \, \mathrm{d}x_{\text{FP}})_{i}
/// \f]
/// where:
/// - \f$\phi_{\text{HP}}\f$: Values in HiPS parcels.
/// - \f$\mathrm{d}x_{\text{HP}}\f$: Differential volume elements for HiPS parcels.
/// - \f$\phi_{\text{FP}}\f$: Values in flow particles.
/// - \f$\mathrm{d}x_{\text{FP}}\f$: Differential volume elements for flow particles.
///
/// \param vh           Vector of values from HiPS parcels to be projected back.
/// \return             Vector of values redistributed onto the flow particles.
///
/// \note 
/// - This function is the reverse of the projection function, ensuring consistency in 
///   value transfers between HiPS parcels and flow particles.
/// - The input vector \p vh should be consistent with the HiPS parcel structure.
///
/// \warning 
/// - Ensure that the HiPS parcels have been properly populated with values before invoking this function.
/// - Mismatches in data sizes between HiPS parcels and flow particles may lead to unexpected results.
//////////////////////////////////////////////////////////////////////////////

//  cfd
//  i=    0         1         2         3
//   |    *    |    *    |    *    |    *    |

//  hips
//   | * | * | * | * | * | * | * | * | * | * |
//  j= 0   1   2   3   4   5   6   7   8   9

std::vector<double> hips::projection_back(std::vector<double> &vh) {

    int nh = xh.size() - 1;
    int nc = xc.size() - 1;

    std::vector<double> vc(nc, 0.0);
    int jprev = 0;

    for (int i = 0; i < nc; ++i) {
        for (int j = jprev + 1; j <= nh; ++j) {
            if (xh[j] <= xc[i + 1]) {
                double d1 = xh[j] - xh[j - 1];
                double d2 = xh[j] - xc[i];
                // calculation of shortest distance
                // handling if distance is zero but substarction gives comutational error
                double d  = std::min(d1, d2) < 1E-15 ? 0 : std::min(d1, d2);

                vc[i] += vh[j - 1] * d;
            } else {
                double d1 = xc[i + 1] - xh[j - 1];
                double d2 = xc[i + 1] - xc[i];
                // calculation of shortest distance
                // handling if distance is zero but substarction gives comutational error
                double d  = std::min(d1, d2) < 1E-15 ? 0 : std::min(d1, d2);

                vc[i] += vh[j - 1] * d;

                jprev = j - 1;
                break;
            }
        }
        vc[i] /= (xc[i + 1] - xc[i]);
    }
    return vc;
}

/////////////////////////////////////////////////////////////////////////////////////////
/// \brief Projects HiPS parcel values and densities back onto the flow particles.
///
/// This function reverses the projection process, redistributing both the values and densities stored 
/// in the HiPS parcels back to the flow particles. It ensures conservation of both the property values 
/// and densities, maintaining consistency between the HiPS parcels and flow particles.
///
/// ### Conservation Principle:
/// The projection follows the equation:
/// \f[
/// \sum_{j=0}^{\text{Number of HP}} (\phi_{\text{HP}} \, \rho_{\text{HP}} \, \mathrm{d}x_{\text{HP}})_{j} = 
/// \sum_{i=0}^{\text{Number of FP}} (\phi_{\text{FP}} \, \rho_{\text{FP}} \, \mathrm{d}x_{\text{FP}})_{i}
/// \f]
/// where:
/// - \f$\phi_{\text{HP}}\f$: Values in HiPS parcels.
/// - \f$\rho_{\text{HP}}\f$: Densities in HiPS parcels.
/// - \f$\mathrm{d}x_{\text{HP}}\f$: Differential volume elements for HiPS parcels.
/// - \f$\phi_{\text{FP}}\f$: Values in flow particles.
/// - \f$\rho_{\text{FP}}\f$: Densities in flow particles.
/// - \f$\mathrm{d}x_{\text{FP}}\f$: Differential volume elements for flow particles.
///
/// \param vh           Vector of values from HiPS parcels to be projected back.
/// \param rho_h        Vector of density values from HiPS parcels.
/// \return             A pair of vectors:
///                     - The first vector contains the values projected back onto the flow particles.
///                     - The second vector contains the densities redistributed to the flow particles.
///
/// \note 
/// - This function is the reverse of the projection function that includes density.
/// - The input vectors \p vh and \p rho_h should be consistent with the HiPS parcel structure and sizes.
///
/// \warning 
/// - Ensure that the HiPS parcels are populated with valid values and densities before invoking this function.
/// - Mismatches in data sizes between HiPS parcels and flow particles may lead to inaccurate results.
/////////////////////////////////////////////////////////////////////////////////////////////////////////

std::pair<std::vector<double>, std::vector<double>> hips::projection_back_with_density(std::vector<double> &vh, 
                                                                                       std::vector<double> &rho_h) {
    int nh = xh.size() - 1;
    int nc = xc.size() - 1;

    std::vector<double> vc(nc, 0.0);
    std::vector<double> rho_c(nc, 0.0);
    int jprev = 0;

    for (int i = 0; i < nc; ++i) {

        for (int j = jprev + 1; j <= nh; ++j) {
            if (xh[j] <= xc[i + 1]) {
                double d = std::min(xh[j] - xh[j - 1], xh[j] - xc[i]);
                rho_c[i] += rho_h[j - 1] * d;
                vc[i] += vh[j - 1] * rho_h[j - 1] * d;
            } else {
                double d = std::min(xc[i + 1] - xh[j - 1], xc[i + 1] - xc[i]);
                rho_c[i] += rho_h[j - 1] * d;
                vc[i] += vh[j - 1] * rho_h[j - 1] * d;

                jprev = j - 1;
                break;
            }
        }

        // Normalize the results
        vc[i] /= rho_c[i];
    }
    return {vc, rho_c};
}

/////////////////////////////////////////////////////////////////////////////////////////
/// \brief Retrieves the final data from the simulation.
///
/// This function returns the final state of the simulation, packaged as a vector of vectors. 
/// It is particularly useful when integrating HiPS as a subgrid model in CFD simulations, 
/// enabling seamless transfer of data for further analysis or post-processing.
///
/// \return A vector of vectors containing the final results, where:
///         - Each inner vector represents a specific variable or property.
///         - The outer vector contains all variables across parcels.
///
/// \note 
/// - This function is specifically designed for use when HiPS is employed as a subgrid model in CFD simulations.
/// - The structure of the returned data ensures compatibility with CFD solvers that require parcel-level data.
///
/// \warning 
/// - Ensure the simulation has reached completion before calling this function to avoid incomplete or inconsistent data.
/// - The returned data structure should be interpreted according to the simulation setup and variable ordering.
//////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> hips::get_varData(){

    std::vector<std::vector<double>> varDataProjections;               // Vector to store modified data projections

    for (int i = 0; i < varData.size(); i++)                           // Loop through each element of varData and project the data back
        varDataProjections.push_back(projection_back((*varData[i])));  // Project the data back and store it in vh

    return varDataProjections;                                         // Return the vector containing modified data projections
}
////////////////////////////////////////////////////////////////////////////////////
/// \brief Retrieves final simulation data, including both values and densities.
///
/// This function processes the HiPS data and projects both the values and densities back 
/// onto the flow particles. It returns a vector of pairs, where each pair contains the 
/// values and corresponding densities for a specific variable. This function is designed 
/// to support HiPS as a subgrid model in CFD simulations, ensuring compatibility with 
/// solvers that require both values and density information.
///
/// \return A vector of pairs:
///         - Each pair consists of two vectors:
///             - The first vector contains the final results for the values.
///             - The second vector contains the corresponding density results.
///
/// \note 
/// - This function is tailored for integrating HiPS as a subgrid model in CFD simulations, 
///   providing both value and density data for accurate modeling.
/// - Ensure that all HiPS parcels are properly initialized and contain valid values and densities 
///   before invoking this function.
///////////////////////////////////////////////////////////////////////////////////

std::vector<std::pair<std::vector<double>, std::vector<double>>> hips::get_varData_with_density() {

    std::vector<std::pair<std::vector<double>, std::vector<double>>> varDataProjections;  // Vector to store data projections with densities

    for (int i = 0; i < varData.size(); i++) {
        // Extract the value and density data
        std::vector<double> vh = (*varData[i]);  // Assuming varData[i][0] holds the values
        std::vector<double> rho_h = varRho;      // Assuming varData[i][1] holds the densities

        // Project the data back with density and store the result
        varDataProjections.push_back(projection_back_with_density(vh, rho_h));
    }

    return varDataProjections;  // Return the vector containing modified data projections with densities
}

/////////////////////////////////////////////////////////////////////////////////////
/// \brief Sets the interval (in number of eddy events) for writing simulation data.
///
/// Calling this function enables eddy-based writing and disables time-based writing.
/// If the user does not call this function or `setOutputIntervalTime()`, the default 
/// interval is set to 1000 eddy events.
///
/// \param interval The number of eddy events between data writes (e.g., `1000` writes data every 1000 eddies).
///
/// \note Calling this function automatically disables time-based writing (`setOutputIntervalTime()`).
///       To revert to default settings, the user must explicitly set a new interval or avoid calling this function.
///////////////////////////////////////////////////////////////////////////////////////////////////

void hips::setOutputIntervalEddy(int interval) {

    outputIntervalEddy = interval;
    useEddyBasedWriting = true;  ///< Enables eddy-based writing
    useTimeBasedWriting = false; ///< Disables time-based writing
}

/////////////////////////////////////////////////////////////////////////////////
/// \brief Sets the interval (in simulation time) for writing simulation data.
///
/// Calling this function enables time-based writing and disables eddy-based writing.
/// If the user does not call this function or `setOutputIntervalEddy()`, the default 
/// behavior is eddy-based writing every 1000 eddies.
///
/// \param interval The time interval (in seconds) between data writes (e.g., `0.1` writes data every 0.1s).
///
/// \note Calling this function automatically disables eddy-based writing (`setOutputIntervalEddy()`).
///       To revert to default settings, the user must explicitly set a new interval or avoid calling this function.
///////////////////////////////////////////////////////////////////////////////////////

void hips::setOutputIntervalTime(double interval) {

    outputIntervalTime = interval;
    useTimeBasedWriting = true;  ///< Enables time-based writing
    useEddyBasedWriting = false; ///< Disables eddy-based writing
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Saves all user-defined to a file.
///
/// This function writes both input parameters (provided by the user) for post-processing. The file is stored in the `post/` directory.
///
/// The file includes:
/// - **User-defined input parameters**, such as grid levels, domain size, turbulence settings, and variable names.
///
/// \note The parameters are saved to `../post/parameters.dat`. Ensure that the `post/` directory exists, 
///       or the function may fail to write the file.
///
/// \warning If the file cannot be openedx, an error message is printed, and no data is saved.
//////////////////////////////////////////////////////////////////////////////////////////////////

void hips::saveAllParameters() {

    std::string filepath = "../post/parameters.dat";  ///< Output file path for simulation parameters
    std::ofstream file(filepath);

    if (!file) {
        std::cerr << "Error: Could not open " << filepath << " for writing!\n";
        return;
    }

    // Write user-defined input parameters
    file << "nLevels " << nLevels << "\n";            ///< Number of hierarchical levels in the HiPS model
    file << "domainLength " << domainLength << "\n";  ///< Length of the computational domain
    file << "tau0 " << tau0 << "\n";                  ///< Reference eddy turnover time
    file << "C_param " << C_param << "\n";            ///< Model constant controlling turbulence behavior
    file << "forceTurb " << forceTurb << "\n";        ///< Flag for forced turbulence (1 = enabled, 0 = disabled)
    file << "nVar " << nVar << "\n";                  ///< Number of variables tracked in the simulation
    file << "performReaction " << performReaction << "\n";  ///< Flag indicating whether chemical reactions are simulated
    file << "realization " << realization << "\n";    ///< Current simulation realization (for multiple runs)

    // Write variable names
    if (!varName.empty()) {
        file << "varName ";
        for (const std::string &name : varName) {
            file << name << " ";  ///< Separate variable names by spaces
        }
        file << "\n";
    } else {
        file << "varName (undefined)\n";
    }

    // Write i_batchelor vector if it exists and has the correct size
    if (!i_batchelor.empty() && i_batchelor.size() == nVar) {
        file << "i_batchelor ";
        for (const auto &val : i_batchelor) {
            file << val << " ";  ///< Print each element separated by space
        }
        file << "\n";
    } else {
        file << "i_batchelor (undefined or size mismatch)\n";
    }

    file.close();
    std::cout << endl << "All parameters saved in: " << filepath << std::endl;
}
/////////////////////////////////////////////////////////////////////////////
