#include "batchReactor_cvode.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////

int rhsf_cvode(realtype t, N_Vector varsCV, N_Vector dvarsdtCV, void *user_data);

///////////////////////////////////////////////////////////////////////////////
/** Constructor
  */

batchReactor_cvode::batchReactor_cvode(std::shared_ptr<Cantera::Solution> cantSol) {

    gas = cantSol->thermo(); 
    kin = cantSol->kinetics(); 

    nvar = gas->nSpecies();

    vector<double> atol(nvar, 1E-12);
    double         rtol = 1E-4;

    integrator = make_unique<integrator_cvode>(rhsf_cvode, this, nvar, rtol, atol);
}
////////////////////////////////////////////////////////////////////////////////

void batchReactor_cvode::react(double &h, vector<double> &y, const double tRun) {

    h_fixed = h;
    P_fixed = gas->pressure();

    integrator->integrate(y, tRun);
}

////////////////////////////////////////////////////////////////////////////////

int batchReactor_cvode::rhsf(const double t, const double *vars, double *dvarsdt) {

    gas->setMassFractions_NoNorm(vars);
    gas->setState_HP(h_fixed, P_fixed);

    double rho = gas->density();
    vector<double> rr(nvar);
    kin->getNetProductionRates(&rr[0]);
    for (size_t k=0; k < gas->nSpecies(); k++)
        dvarsdt[k] = rr[k] * gas->molecularWeight(k) / rho;

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// CVODE interface; CVODE calls this function, which then calls user_data's rhsf

int rhsf_cvode(realtype t, N_Vector varsCV, N_Vector dvarsdtCV, void *user_data) {

    batchReactor_cvode *bRxr = static_cast<batchReactor_cvode *>(user_data);

    double *vars    = N_VGetArrayPointer(varsCV);
    double *dvarsdt = N_VGetArrayPointer(dvarsdtCV);

    int rv = bRxr->rhsf(t, vars, dvarsdt);

    return rv;
}
