#pragma once

#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "integrator_cvode.h"

class batchReactor_cvode {

////////////////////////////// DATA MEMBERS /////////////////////////////

private:

    std::shared_ptr<Cantera::ThermoPhase> gas;
    std::shared_ptr<Cantera::Kinetics>    kin;
    std::unique_ptr<integrator_cvode>     integrator;

    int                                   nvar;

    double                                h_fixed;    // adiabatic h during integrate
    double                                P_fixed;    // pressure during integrate

/////////////////////////////// MEMBER FUNCTIONS//////////////////////////////////////////

public:

    batchReactor_cvode(std::shared_ptr<Cantera::Solution> cantSol);

    void react(double &h, std::vector<double> &y, const double tRun);

    int rhsf(const double t, const double *vars, double *dvarsdt);

    virtual ~batchReactor_cvode() {}

};
