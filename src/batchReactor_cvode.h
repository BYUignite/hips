#pragma once

#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "integrator_cvode.h"

#include <memory>
#include <vector>

////////////////////////////////////////////////////////////////////////////////

class batchReactor_cvode {

////////////////////////////// DATA MEMBERS /////////////////////////////

private:

    std::shared_ptr<Cantera::ThermoPhase> gas;        ///< Cantera thermo object
    std::shared_ptr<Cantera::Kinetics>    kin;        ///< Cantera kinetics object
    std::unique_ptr<integrator_cvode>     integrator; ///< cvode integrator wrappter

    int                                   nvar;       ///< number of variables/equations solved

    double                                h_fixed;    ///< adiabatic h during integrate
    double                                P_fixed;    ///< pressure during integrate

/////////////////////////////// MEMBER FUNCTIONS//////////////////////////////////////////

public:

    batchReactor_cvode(std::shared_ptr<Cantera::Solution> cantSol);

    void react(double &h, std::vector<double> &y, const double tRun);

    int rhsf(const double t, const double *vars, double *dvarsdt);  // dydt = rhsf

    virtual ~batchReactor_cvode() {}

};
