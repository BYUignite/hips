#pragma once

#include "BatchReactor.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "integrator_cvode.h"

#include <memory>
#include <vector>

////////////////////////////////////////////////////////////////////////////////

class BatchReactor_cvode : public BatchReactor {

////////////////////////////// DATA MEMBERS /////////////////////////////
    
    std::shared_ptr<integrator_cvode>     integrator;     ///< cvode integrator wrapper
    
/////////////////////////////// MEMBER FUNCTIONS//////////////////////////////////////////

public:

    BatchReactor_cvode(std::shared_ptr<Cantera::Solution> cantSol);

    virtual void react(double &h, std::vector<double> &y, const double tRun);

    int rhsf(const double t, const double *vars, double *dvarsdt);  // dydt = rhsf
    
    virtual ~BatchReactor_cvode() {}
};
