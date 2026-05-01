
/// \file ex_2.cc
/// \brief Example demonstrating premixed combustion using HiPS and Cantera.
/// 
/// This example simulates turbulent premixed flames using a stoichiometric 
/// ethylene/air mixture. It initializes a HiPS domain with six levels, where 
/// 25% of the parcels are pre-combusted, and the remaining 75% contain fresh 
/// reactants. The simulation uses Cantera to handle chemical kinetics.
/// 
/// ## Compilation:
/// ```bash
/// cd build
/// cmake .. -DREACTIONS_ENABLED=ON
/// make
/// sudo make install
/// ```
/// 
/// ## Execution:
/// ```bash
/// cd ../run
/// ./ex_2.x
/// ```
/// 
/// ## Output:
/// The simulation produces temperature and species concentration profiles, 
/// saved in the output file.
/// 
/// ## Dependencies:
/// - HiPS library (compiled with reaction support)
/// - Cantera library (for chemical kinetics)

///////////////////////////////////////////////////////////////////////////////////////
#include "hips.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <memory>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////
int main() {
    // HiPS simulation parameters for premixed combustion
    int nLevels = 8;             // Number of hierarchical levels
    double domainLength = 0.01;  // Domain length scale
    double tau0 = 0.00005;      // Mixing timescale (fast mixing)
    double C_param = 0.5;        // Eddy rate multiplier
    double tRun = 0.002;             // Total simulation runtime
    bool forceTurb = false;           // No forced turbulence
    vector<double> ScHips(54, 1); // Schmidt number (unity for all species)

    // Initialize Cantera for combustion chemistry
    auto cantSol = Cantera::newSolution("gri30.yaml");
    auto gas = cantSol->thermo();
    size_t nsp = gas->nSpecies();
    int nVar = nsp + 1; // Number of variables (species + enthalpy)
    double P = 101325.0;
    int seed = 11;

    // Create HiPS instance with reaction support
    hips HiPS(nLevels, domainLength, tau0, C_param, forceTurb, nVar, ScHips, true, cantSol, P, seed);

    int nparcels = HiPS.get_nparcels();

    // Define variables for species mass fractions, temperature, and enthalpy
    vector<vector<double>> ysp(nsp, vector<double>(nparcels, 0));
    vector<double> h(nparcels);
    vector<double> rho(nparcels);

    vector<double> y0(nsp), y1(nsp); // Initial species mass fractions
    double T0 = 300.0, T1 = 300.0;   // Initial temperature
    double h0, h1;                   // Enthalpy values for reactants and products

    vector<string> variableNames(nsp + 1); // Variable names: enthalpy and species

    // Fraction of parcels that are pre-combusted
    double fracBurn = 0.25;

    // Initialize fresh reactants
    gas->setState_TPX(T0, Cantera::OneAtm, "C2H4:1, O2:3, N2:11.25");
    h0 = gas->enthalpy_mass();
    gas->getMassFractions(y0.data());

    // Assign fresh reactant properties to parcels
    for (int i = 0; i < nsp; i++) {
        for (int j = 0; j <= (1 - fracBurn) * nparcels; j++) {
            h[j] = h0;
            ysp[i][j] = y0[i];
            rho[j] = gas->density();
        }
    }

    // Initialize pre-combusted parcels
    gas->setState_TPX(T1, Cantera::OneAtm, "C2H4:1, O2:3, N2:11.25");
    h1 = gas->enthalpy_mass();
    gas->setState_HP(h1, gas->pressure());
    gas->equilibrate("HP");
    gas->getMassFractions(y1.data());

    // Assign pre-combusted properties to parcels
    for (int i = 0; i < nsp; i++) {
        for (int j = ((1 - fracBurn) * nparcels + 1); j < nparcels; j++) {
            h[j] = h1;
            ysp[i][j] = y1[i];
            rho[j] = gas->density();
        }
    }

    // Set initial conditions in HiPS

    variableNames[0] = "enthalpy";
    for (int i = 0; i < ysp.size(); i++) {
        variableNames[i + 1] = gas->speciesName(i);
    }

    vector<double> weight(nparcels, 1.0 / nparcels); // Uniform weights

    // Assign variables to HiPS
    HiPS.set_varData(h, weight, variableNames[0], rho);
    for (int k = 0; k < ysp.size(); k++) 
        HiPS.set_varData(ysp[k], weight, variableNames[k + 1], rho);

    // Set output interval in terms of time
    HiPS.setOutputIntervalTime(tRun/10);  // Save results every 100 seconds

    // Write initial condition
    HiPS.writeData(1, 0, 0.0 );
  
    // Run the combustion simulation
    HiPS.calculateSolution(tRun, true);

    return 0;
}

