#pragma once

#include "batchReactor.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "integrator_cvode.h"

#include <memory>
#include <vector>

////////////////////////////////////////////////////////////////////////////////

class batchReactor_cvode : public batchReactor {

////////////////////////////// DATA MEMBERS /////////////////////////////
    
    std::unique_ptr<integrator_cvode>     integrator; ///< cvode integrator wrappter

/////////////////////////////// MEMBER FUNCTIONS//////////////////////////////////////////

public:

    batchReactor_cvode(std::shared_ptr<Cantera::Solution> cantSol);

    virtual void react(double &h, std::vector<double> &y, const double tRun);

    int rhsf(const double t, const double *vars, double *dvarsdt);  // dydt = rhsf

    virtual ~batchReactor_cvode() {}
};
