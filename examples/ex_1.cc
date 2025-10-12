
/// \file ex_1.cc
/// \brief Example demonstrating scalar mixing dynamics using HiPS.
/// 
/// This example illustrates the basic behavior of scalar mixing dynamics. 
/// Fluid parcels are initially divided into two groups (0 and 1) and then 
/// mix over time to form a uniform scalar field. The example uses the HiPS 
/// model to demonstrate how Schmidt numbers affect mixing dynamics.
/// 
/// ## Compilation:
/// ```bash
/// cd build
/// cmake ..
/// make
/// sudo make install
/// ```
/// 
/// ## Execution:
/// ```bash
/// cd ../run
/// ./ex_1.x
/// ```
/// 
/// ## Output:
/// The simulation produces scalar concentration profiles stored in the output file.

///////////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include "HiPS.h"

using namespace std;

// Initializes mixing fractions for fluid parcels
// - numParcels: Number of parcels in the simulation
// - Returns a vector with initialized mixing fractions (0 or 1)

vector<double> initializeMixingFractions(int numParcels) {
    vector<double> mixingFractions(numParcels);
    for (int i = 0; i < numParcels; i++) {
        // Assign 0.0 to the first half and 1.0 to the second half
        mixingFractions[i] = (i < numParcels / 2) ? 0.0 : 1.0;
    }
    return mixingFractions;
}

///////////////////////////////////////////////////////////////////////////////////////
int main() {
    // HiPS tree and simulation parameters
    int nLevels = 9;                           // Number of hierarchical levels
    double domainLength = 1.0;                  // Simulation domain length
    double tau0 = 1.0;                         // Initial time scale for the largest eddies
    double C_param = 0.5;                      // Eddy turnover rate multiplier
    double tRun = 300.0;                       // Total simulation time
    int forceTurb = 2;                         // Forcing parameter to impose turbulent profile
    vector<double> ScHips = {0.0625, 1.0, 16.0};    // Schmidt numbers for low and high diffusivity
    int numVariables = 3;                      // Number of scalar fields

    // Set up HiPS tree and calculate the number of parcels
    HiPS hips(nLevels, domainLength, tau0, C_param, forceTurb, numVariables, ScHips, false);
    int numParcels = hips.get_nparcels()
;

    // Initialize mixing fractions for each scalar variable
    vector<vector<double>> mixingFractions(numVariables);
    vector<double> weights(numParcels, 1.0 / numParcels);  // Uniform weights

    for (int i = 0; i < numVariables; ++i) {
        // Initialize mixing fractions for variable i
        mixingFractions[i] = initializeMixingFractions(numParcels);
        hips.set_varData(mixingFractions[i], weights, "mixf_0" + to_string(i));
    }

    // Set output interval in terms of time
    hips.setOutputIntervalTime(60.0);  // Save results every 100 seconds

    // Write initial condition
    hips.writeData(1, 0, 0.0 );

    // Run the simulation and calculate mixing dynamics
    hips.calculateSolution(tRun, true);

    return 0;
}

