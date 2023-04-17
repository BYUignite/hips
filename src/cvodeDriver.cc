
#include <iostream>
#include "cvodeDriver.h"

#include <vector>

using namespace Cantera;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

cvodeDriver::cvodeDriver(shared_ptr<Solution> sol) {
    gas = sol->thermo();
    kin = sol->kinetics();
    //cout <<"test 1 in cvodeD  "<<endl;
    Neq = gas->nSpecies();
    rr.resize(gas->nSpecies());

    integrator = unique_ptr<Integrator>(newIntegrator("CVODE"));
    integrator->initialize(0.0, *this);
}

////////////////////////////////////////////////////////////////////////////////

void cvodeDriver::integrate(double dt) {

    h_fixed = gas->enthalpy_mass();
    //cout <<"temperature    "<<gas->temperature()<<endl;
    P_fixed = gas->pressure();
  //  cout <<"test 2 in cvodeD  "<<endl;
    integrator->reinitialize(0.0, *this);

    integrator->integrate(dt);

    double *solution = integrator->solution();
    //doldb gas->setMassFractions_NoNorm(solution); // shouldn't use nonorm after integration is done
    gas->setMassFractions(solution);
    gas->setState_HP(h_fixed, P_fixed);

//   cout<<"gas of o2 "<<gas->massFraction(0) <<endl;
//cout<<"gas of  "<<gas->massFraction(1) <<endl;
//
//   cout<<"gas of CH4 "<<gas->massFraction(2) <<endl;
//cout<<"gas of  "<<gas->massFraction(3) <<endl;
//
//cout<<"gas of  "<<gas->massFraction(4) <<endl;
//
//   cout<<"gas of N2 "<<gas->massFraction(5) <<endl;
}

////////////////////////////////////////////////////////////////////////////////

void cvodeDriver::eval(double t, double *y, double *dydt, double *not_used) {

    double *massFracs = &y[0];
    double *dYdt = &dydt[0];
   
    gas->setMassFractions_NoNorm(massFracs);
    
    gas->setState_HP(h_fixed, P_fixed);
  
  // cout <<"test 3 in cvodeD  "<<endl;

    double rho = gas->density();
    kin->getNetProductionRates(&rr[0]);
    for (size_t k = 0; k < gas->nSpecies(); k++)
        dYdt[k] = rr[k] * gas->molecularWeight(k) / rho;
}
