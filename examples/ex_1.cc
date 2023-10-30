
#include "hips.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include <vector>
#include <string>
#include <iostream>

using namespace std;

/**
 * @brief Initialize mixing fractions.
 *
 * This function initializes mixing fractions for a given number of parcels.
 *
 * @param numParcels The number of parcels.
 * @return A vector containing initialized mixing fractions.
 */

vector<double> initializeMixingFractions(int numParcels) {
    vector<double> mixingFractions;
    for (int i = 0; i < numParcels; i++) {
        // Initialize mixing fractions to 0.0 for the first half and 1.0 for the second half
        mixingFractions.push_back(i < numParcels / 2 ? 0.0 : 1.0);
    }
    return mixingFractions;
}

///////////////////////////////////////////////////////////////////////////////

int main() {

    //---------- HiPS tree and constructor parameters

    int numLevels = 12;
    double domainLength = 1.0;
    double tau0 = 1.0;
    double C_param = 0.5;
    double tRun = 700.0;
    int forceTurb = 2;
    vector<double> ScHips(1, 1);

    //---------- Gas solution setup

    auto cantSol = Cantera::newSolution("gri30.yaml");
    auto gas = cantSol->thermo();
    size_t numSpecies = gas->nSpecies();
    int numVariables = 1;  

    //---------- HiPS tree creation

    hips HiPS(numLevels, domainLength, tau0, C_param, forceTurb, numVariables, ScHips, 
              cantSol, false);
    int numParcels = HiPS.nparcels;

    //---------- Initialize mixing fractions

    vector<double> initialMixingFractions = initializeMixingFractions(numParcels);

    //---------- Set state vectors in each parcel with initial mixing fractions

    HiPS.set_varData(initialMixingFractions, 0);

    //---------- Advance HiPS for mixing 

    HiPS.calculateSolution(tRun);

    return 0;
}

