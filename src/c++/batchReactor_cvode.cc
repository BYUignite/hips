#include "batchReactor_cvode.h"

///////////////////////////////////////////////////////////////////////////////
/// \brief Callback function for the right-hand side of the ODE system used by CVODE.
///
/// This function calculates the right-hand side of the ODE system.
/// 
/// \param t Current time.
/// \param varsCV N_Vector containing the variables (species concentrations).
/// \param dvarsdtCV N_Vector containing the derivatives of the variables with respect to time.
/// \param user_data Pointer to the user-defined data (batchReactor_cvode instance).
/// \return 0 on success.
////////////////////////////////////////////////////////////////////////////////

int rhsf_cvode(sunrealtype t, N_Vector varsCV, N_Vector dvarsdtCV, void *user_data);

///////////////////////////////////////////////////////////////////////////////
/// \brief Constructor for the batchReactor_cvode class.
/// 
/// Initializes a batch reactor object using CVODE for simulations.
/// 
/// \param cantSol A shared pointer to a Cantera solution object.
/// 
///////////////////////////////////////////////////////////////////////////////////

batchReactor_cvode::batchReactor_cvode(std::shared_ptr<Cantera::Solution> cantSol) {

    gas = cantSol->thermo(); 
    kin = cantSol->kinetics(); 

    nvar = gas->nSpecies();

    // Set absolute and relative tolerances
    std::vector<double> atol(nvar, 1E-10);
    double rtol = 1E-4;

    // Initialize the integrator with CVODE
    integrator = std::make_shared<integrator_cvode>(rhsf_cvode, this, nvar, rtol, atol);
}

///////////////////////////////////////////////////////////////////////////////
/// \brief Simulates a reaction in the batch reactor.
/// 
/// \param h Heat of the reactor.
/// \param y Vector of species mass fractions.
/// \param tRun Time for the simulation.
///
////////////////////////////////////////////////////////////////////////////////////////

void batchReactor_cvode::react(double &h, std::vector<double> &y, const double tRun) {

    // Store fixed enthalpy and pressure
    h_fixed = h;
    P_fixed = gas->pressure();

    // Integrate over the given time
    integrator->integrate(y, tRun);
}

///////////////////////////////////////////////////////////////////////////////
/// \brief Computes the right-hand side of the ODE system.
/// 
/// This function is used by CVODE to compute the derivatives of the variables with respect to time.
/// 
/// \param t Current time.
/// \param vars Pointer to an array containing the variables (species concentrations).
/// \param dvarsdt Pointer to an array where the computed derivatives will be stored.
/// \return 0 on success.
////
/////////////////////////////////////////////////////////////////////////////////////

int batchReactor_cvode::rhsf(const double t, const double *vars, double *dvarsdt) {

    // Set mass fractions and state
    gas->setMassFractions_NoNorm(vars);
    gas->setState_HP(h_fixed, P_fixed);

    // Calculate density and reaction rates
    rho = gas->density();
    temperature = gas->temperature();

    std::vector<double> rr(nvar);
    kin->getNetProductionRates(&rr[0]);
    for (size_t k = 0; k < gas->nSpecies(); k++)
        dvarsdt[k] = rr[k] * gas->molecularWeight(k) / rho;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// CVODE interface; CVODE calls this function, which then calls user_data's rhsf

/// \brief Callback function for the right-hand side of the ODE system used by CVODE.
/// 
/// This function is called by CVODE, which in turn calls the rhsf method of the batchReactor_cvode instance.
/// 
/// \param t Current time.
/// \param varsCV N_Vector containing the variables (species concentrations).
/// \param dvarsdtCV N_Vector containing the derivatives of the variables with respect to time.
/// \param user_data Pointer to the user-defined data (batchReactor_cvode instance).
/// \return 0 on success.
/////////////////////////////////////////////////////////////////////////////////

int rhsf_cvode(sunrealtype t, N_Vector varsCV, N_Vector dvarsdtCV, void *user_data) {
 
    batchReactor_cvode *bRxr = static_cast<batchReactor_cvode *>(user_data);       // Cast user_data pointer to batchReactor_cvode pointer

    double *vars    = N_VGetArrayPointer_Serial(varsCV);                                 // Get array pointers for variables and derivatives

    double *dvarsdt = N_VGetArrayPointer_Serial(dvarsdtCV);                              // Call rhsf method of batchReactor_cvode instance
    
    int rv = bRxr->rhsf(t, vars, dvarsdt);

    return rv;
}
