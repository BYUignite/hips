#pragma once

#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "integrator_cvode.h"

#include <memory>
#include <vector>

////////////////////////////////////////////////////////////////////////////////

class batchReactor {

////////////////////////////// DATA MEMBERS /////////////////////////////

protected:

    std::shared_ptr<Cantera::ThermoPhase> gas;        ///< Cantera thermo object
    std::shared_ptr<Cantera::Kinetics>    kin;        ///< Cantera kinetics object

    int                                   nvar;       ///< number of variables/equations solved

    double                                h_fixed;    ///< adiabatic h during integrate
    double                                P_fixed;    ///< pressure during integrate

/////////////////////////////// MEMBER FUNCTIONS//////////////////////////////////////////

public:

    batchReactor() {};

    virtual void react(double &h, std::vector<double> &y, const double tRun) = 0;

    virtual ~batchReactor() {}
};
