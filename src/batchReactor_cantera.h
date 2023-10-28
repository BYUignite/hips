#pragma once

#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/numerics/Integrator.h"

#include <memory>
#include <vector>

////////////////////////////////////////////////////////////////////////////////

class batchReactor_cantera : public Cantera::FuncEval {

    ////////////////////// DATA MEMBERS /////////////////////

private:

    std::shared_ptr<Cantera::ThermoPhase> gas;        ///< Cantera thermo object
    std::shared_ptr<Cantera::Kinetics>    kin;        ///< Cantera kinetics object
    std::unique_ptr<Cantera::Integrator>  integrator; ///< Cantera cvode wrapper

    int                                   nvar;       ///< number of variables/equations solved

    double                                h_fixed;    ///< adiabatic h during integrate
    double                                P_fixed;    ///< pressure during integrate


    ////////////////////// MEMBER FUNCTIONS /////////////////

public: 

    batchReactor_cantera(std::shared_ptr<Cantera::Solution> cantSol);

    void react(double &h, std::vector<double> &y, const double tRun);

    void eval(double t, double *vars, double *dvarsdt, double *not_used); // rhsf: dydt = rhsf

    size_t neq() { return nvar; }          // called by Cantera

    void getState(double* y) {             // called by cantera to set y
        gas->getMassFractions(y);
    }

    virtual ~batchReactor_cantera() {}

};
