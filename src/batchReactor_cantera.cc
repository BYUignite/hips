#include "batchReactor_cantera.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
/** Constructor
  */

batchReactor_cantera::batchReactor_cantera(std::shared_ptr<Cantera::Solution> cantSol) {

    gas = cantSol->thermo(); 
    kin = cantSol->kinetics(); 

    nvar = gas->nSpecies();

    integrator = unique_ptr<Cantera::Integrator>(Cantera::newIntegrator("CVODE"));
    integrator->setTolerances(1E-4, 1E-10);       // rtol, atol
    integrator->setMaxSteps(5000);
    integrator->initialize(0.0, *this);

}
////////////////////////////////////////////////////////////////////////////////

void batchReactor_cantera::react(double &h, vector<double> &y, const double tRun) {

    gas->setMassFractions(&y[0]);
    gas->setState_HP(h, gas->pressure());

    h_fixed = h;
    P_fixed = gas->pressure();

    integrator->reinitialize(0.0, *this);

    integrator->integrate(tRun);

    double *solution = integrator->solution();
    for (int k=0; k<nvar; k++)
        y[k] = solution[k];
}

////////////////////////////////////////////////////////////////////////////////

void batchReactor_cantera::eval(double t, double *vars, double *dvarsdt, double *not_used) {

    gas->setMassFractions_NoNorm(vars);
    gas->setState_HP(h_fixed, P_fixed);

    vector<double> rr(gas->nSpecies());

    double rho = gas->density();
    kin->getNetProductionRates(&rr[0]);
    for (size_t k=0; k < gas->nSpecies(); k++)
        dvarsdt[k] = rr[k] * gas->molecularWeight(k) / rho;
}
