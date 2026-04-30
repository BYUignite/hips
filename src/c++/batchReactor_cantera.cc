#include "batchReactor_cantera.h"

///////////////////////////////////////////////////////////////////////////////
/// \brief Constructor for batchReactor_cantera.
/// 
/// Initializes a batch reactor object using Cantera solution for simulations.
/// 
/// \param cantSol A shared pointer to a Cantera solution object.
/// 
///////////////////////////////////////////////////////////////////////

batchReactor_cantera::batchReactor_cantera(std::shared_ptr<Cantera::Solution> cantSol) {
    gas = cantSol->thermo(); 
    kin = cantSol->kinetics(); 

    // Get the number of species in the system
    nvar = gas->nSpecies();

    // Initialize the integrator with CVODE
    integrator.reset(Cantera::newIntegrator("CVODE"));        // Use reset() instead of make_unique
    integrator->setTolerances(1E-4, 1E-10);                   // Set tolerances for CVODE
    integrator->setMaxSteps(5000);                            // Set maximum number of steps
    integrator->initialize(0.0, *this);                       // Initialize the integrator
}

///////////////////////////////////////////////////////////////////////////////
/// \brief Simulates a reaction in the batch reactor.
/// 
/// \param h Heat of the reactor.
/// \param y Vector of species mass fractions.
/// \param tRun Time for the simulation.
///
///////////////////////////////////////////////////////////////////////////////

void batchReactor_cantera::react(double &h, std::vector<double> &y, const double tRun) {
    // Set mass fractions and state
    gas->setMassFractions(&y[0]);
    gas->setState_HP(h, gas->pressure());

    // Store fixed enthalpy and pressure
    h_fixed = h;
    P_fixed = gas->pressure();

    // Reinitialize the integrator
    integrator->reinitialize(0.0, *this);

    // Integrate over the given time
    integrator->integrate(tRun);

    // Get the solution and update species mass fractions
    double *solution = integrator->solution();
    for (int k = 0; k < nvar; k++)
        y[k] = solution[k];
}

///////////////////////////////////////////////////////////////////////////////
/// \brief Evaluates the reaction rates.
/// 
/// \param t Time.
/// \param vars Variables.
/// \param dvarsdt Derivatives of variables with respect to time.
/// \param not_used Unused parameter.
/////////////////////////////////////////////////////////////////////////////////

void batchReactor_cantera::eval(double t, double *vars, double *dvarsdt, double *not_used) {
    // Set mass fractions and state
    gas->setMassFractions_NoNorm(vars);
    gas->setState_HP(h_fixed, P_fixed);

    // Calculate density and reaction rates
    double rho = gas->density();
    std::vector<double> rr(gas->nSpecies());
    kin->getNetProductionRates(&rr[0]);
    for (size_t k = 0; k < gas->nSpecies(); k++)
        dvarsdt[k] = rr[k] * gas->molecularWeight(k) / rho;
}
