#pragma once

#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include <memory>
#include <string>
#include <vector>

////////////////////////////////////////////////////////////////////////////////

class cgas {


    //////////////////// DATA MEMBERS //////////////////////

    public: 

    std::shared_ptr<Cantera::ThermoPhase> gas;

    //////////////////// MEMBER FUNCTIONS /////////////////

    int nSpecies() { return gas->nSpecies(); }
    double enthalpy_mass() { return gas->enthalpy_mass(); }
    double density() { return gas->density(); }
    double pressure() { return gas->pressure(); }
    void setState_TPX(double T, double P, std::string X) { gas->setState_TPX(T, P, X); }
    void setState_HP(double H, double P) { gas->setState_HP(H, P); }
    std::vector<double> getMassFractions() {
        std::vector<double> y(gas->nSpecies());
        gas->getMassFractions(y.data());
        return y;
    }
    void equilibrate(std::string XY) { gas->equilibrate(XY); }
    std::string speciesName(int k) { return gas->speciesName(k); }
    std::vector<std::string> speciesNames() {
        std::vector<std::string> names;
        names.reserve(gas->nSpecies());
        for (size_t k = 0; k < gas->nSpecies(); k++) {
            names.push_back(gas->speciesName(k));
        }
        return names;
    }

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////
    
    cgas(std::string mechanismName) {
        auto sol = Cantera::newSolution(mechanismName);
        gas = sol->thermo();
    }

   ~cgas(){}
};
