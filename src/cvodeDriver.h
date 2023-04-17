#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/base/Solution.h"
#include "cantera/numerics/Integrator.h"
#include <memory>
#include <vector>

using namespace Cantera;

////////////////////////////////////////////////////////////////////////////////

class cvodeDriver : public Cantera::FuncEval {

    ////////////////////// DATA MEMBERS /////////////////////

    private:

    std::shared_ptr<ThermoPhase> gas;        // Cantera thermo object
    std::shared_ptr<Kinetics>    kin;        // Cantera kinetics object
    int                          Neq;        // # equations being solved
    double                       h_fixed;    // adiabatic h during integrate
    double                       P_fixed;    // pressure during integrate
    std::vector<double>          rr;         // reaction rate (kmol/m3*s)

    std::unique_ptr<Integrator>  integrator; // Cantera cvode wrapper

    ////////////////////// MEMBER FUNCTIONS /////////////////

    public: 

    void integrate(double dt);            // integrate to dt; gas holds soln

    void eval(double  t,                  // rhsf: dy/dt = rhsf
              double *y, 
              double *ydot, 
              double *not_used);

    size_t neq() { return Neq; }          // called by Cantera

    void getState(double* y) {            // called by cantera to set y
        gas->getMassFractions(&y[0]);
    }

    ////////////////////// CONSTRUCTOR FUNCTION /////////////

    cvodeDriver(shared_ptr<Solution> sol);

};
