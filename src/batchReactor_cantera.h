#pragma once

#include "batchReactor.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/numerics/Integrator.h"

#include <memory>
#include <vector>

////////////////////////////////////////////////////////////////////////////////

class batchReactor_cantera : public batchReactor, public Cantera::FuncEval {

////////////////////// DATA MEMBERS /////////////////////

std::unique_ptr<Cantera::Integrator>  integrator; ///< Cantera cvode wrapper

////////////////////// MEMBER FUNCTIONS /////////////////

public: 

    batchReactor_cantera(std::shared_ptr<Cantera::Solution> cantSol);

    virtual void react(double &h, std::vector<double> &y, const double tRun);

    void eval(double t, double *vars, double *dvarsdt, double *not_used); // rhsf: dydt = rhsf

    size_t neq() { return nvar; }          // called by Cantera

    void getState(double* y) {             // called by cantera to set y
        gas->getMassFractions(y);
    }

    virtual ~batchReactor_cantera() {}
};
