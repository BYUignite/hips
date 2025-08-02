#pragma once

#ifdef  REACTIONS_ENABLED
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "batchReactor.h"
#endif

#include <iostream>
#include <vector>
#include <string>
#include <vector>
#include <memory>
#include "randomGenerator.h"

class hips {

    ////////////////////////////// DATA MEMBERS /////////////////////////////

public:
  
    int nparcels;                                                  ///< number of parcels
    int realization;                                               ///< number of realizations
    std::vector<std::shared_ptr<std::vector<double>>> varData;     ///< vector of pointers to vector
    std::vector<double> varRho;                                    ///< density
    std::vector<int> pLoc;                                         ///< parcel index array for fast implementation of swaps
 
#ifdef REACTIONS_ENABLED
    std::shared_ptr<Cantera::ThermoPhase> gas;                     ///< Shared pointer to a Cantera thermochemistry object
    std::shared_ptr<batchReactor> bRxr;                            ///< Unique pointer to the integrator object
#endif

    double domainLength;                                           ///< length of domain (m)
    double tau0;                                                   ///< integral timescale
    double C_param;                                                ///< Eddy frequency parameter
    std::vector<double> Temp;                                      ///< Vector containg temperature in each parcel;
    
private:
    
    int currentIndex = 0;                                          ///< member variable to keep track of current index of variables
    int nLevels;                                                   ///< number of tree levels
    int nLevels_;                                                  ///< number of tree levels?
    int forceTurb;                                                 ///< forcing function for statistically stationary: -1 = none, 1 = source term, 2 = dir
    int nVar;                                                      ///< number of parcel variables (e.g., h, ysp)
    int nsp;                                                       ///< number of species
    int Nm1;                                                       ///< nLevels - 1 
    int Nm2;                                                       ///< nLevels - 2                     
    int Nm3;                                                       ///< nLevels - 3 
    int iEta;                                                      ///< Kolmogorov level (needed for variable Sc scalars)
    int nL;                                                        ///< adjusted number of levels based on the Reynolds number

    bool LScHips;                                                  ///< hips schmidt number
    bool performReaction;                                          ///< flag indicating whether chemical reactions are performed in the simulation 
    bool LrandSet;                                                 ///< flag indicating new randomGen  --> allow deletion
    
    double time;                                                   ///< current simulation time
    double eddyRate_total;                                         ///< total rate of all eddies 0 through nLevels-3
    double eddyRate_inertial;                                      ///< total rate of all eddies 0 through iEta (= eddyRate_total if Sc=1) 
    double Afac = 0.5;                                             ///< level lengthscale reduction factor (0.5)
    double Re;                                                     ///< Reynolds number
    double dtEE;                                                   ///< time increment to next eddy event 
    double Prob;                                                   ///< probability value for probability-based solution
    double lStar;                                                  ///< length of the level associated with the Reynolds number 
    double Anew;                                                   ///< adjusted level lengthscale reduction factor for dynamic adjustment of reduction factor

    randomGenerator rand;                                               
    
    std::vector<int> i_plus;                                       ///< ceil(i_batchelor)
    std::vector<double> ScHips;                                    ///< vector containing Schmidt numbers related to each variable
    std::vector<std::string> varName;                              ///< vector containing the names of parcel variables
    std::vector<double> parcelTimes;                               ///< current times corresponding to the parcel states
    std::vector<double> levelRates;                                ///< list of eddy event rates at each level
    std::vector<double> i_batchelor;                               ///< Batchelor level for variable Sc scalars; NOTE: double, as in, between levels
    std::vector<double> xc;                                        ///< vector containing physical domain of flow particles
    std::vector<double> xh;                                        ///< vector containing physical domain of HiPS parcels
    
    std::string  approach;
    int outputIntervalEddy = 10;                                   ///< Default: write data every 10 eddy events
    double outputIntervalTime = 0.1;                               ///< Default: write data every 0.1s
    int eddyCounter = 0;                                           ///< Counter for eddy events
    double lastOutputTime = 0.0;                                   ///< Last time data was written
    bool useEddyBasedWriting = false;                              ///< Tracks if eddy writing is set
    bool useTimeBasedWriting = false;                              ///< Tracks if time writing is set
    const int DEFAULT_EDDY_INTERVAL = 1000;                        ///< Default: Write every 1000 eddies
    const double DEFAULT_TIME_INTERVAL = 0.1;                      ///< Default: Write every 0.1s
  


    ////////////////////////////// MEMBER FUNCTIONS /////////////////////////////

public:

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Sets up the HiPS tree using explicitly specified tree parameters.
    ///
    /// This method builds the binary HiPS tree using a user-defined number of levels and physical parameters.
    /// It also automatically adjusts the number of levels to account for high Schmidt numbers, ensuring that
    /// micromixing is properly resolved.
    ///
    /// \param nLevels_         Base number of levels in the HiPS tree.
    /// \param domainLength_    Domain size for determining eddy length scales.
    /// \param tau0_            Time scale of the smallest eddy (Kolmogorov scale).
    /// \param ScHips_          Vector of Schmidt numbers (one per variable).
    ///
    /// \note This function is used by the full constructor and can also be called manually after 
    ///       using the dynamic constructor.
    ///
    /// \warning The number of levels may be increased automatically for large Schmidt numbers 
    ///          to ensure accurate scalar mixing across scales.

    void set_tree(int nLevels_, double domainLength_, double tau0_);  // setting the HiPS tree based on the number of levels


    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Dynamically builds the HiPS tree based on Reynolds number and a selected strategy.
    ///
    /// This method configures the HiPS tree using a continuous Reynolds number and one of several 
    /// initialization strategies. It is intended for simulations where local turbulence conditions 
    /// change over time or space, requiring the tree to be updated dynamically.
    ///
    /// \param Re_              Reynolds number used to determine base tree level.
    /// \param domainLength_    Domain length for spatial scaling.
    /// \param tau0_            Base time scale for the largest eddy.
    /// \param ScHips_          Vector of Schmidt numbers (one per variable).
    /// \param approach_        Strategy to convert continuous Re to tree level:
    ///                         - "rounding"     → Round to nearest discrete level
    ///                         - "probability"  → Use probabilistic interpolation between levels
    ///                         - "micromixing"  → Use fixed level with adjusted mixing rate
    ///                         - "dynamic_A"    → Adjust geometric scale factor A to fit Re
    ///
    /// \note This method is ideal for Lagrangian simulations using grid cells with different Re values.
    ///       It supports runtime reconfiguration of the tree without reinitializing the hips object.
    ///
    /// \warning Ensure consistent `approach_` handling across the simulation to avoid inconsistencies.
    ///
    /// \see hips::set_tree(int, ...) for direct-level setup.
  
    void set_tree(double Re_, double domainLength_, double tau0_, std::string approach_ = "rounding");
  
////////////////////////////////////////////////////////////////////////////////////////////////////////

    void set_varData(std::vector<double> &v, std::vector<double> &w, const std::string &varN);      // passing all variables to vector of pointer 
    void set_varData(std::vector<double> &v, std::vector<double> &w, const std::string &varN,
                     const std::vector<double> &rho);                                               // passing all variables to vector of pointer

    std::vector<std::vector<double>> get_varData();                                                 // Retrieves modified data from the HiPS library and stores it in the provided vector. 
    std::vector<std::pair<std::vector<double>, std::vector<double>>> get_varData_with_density();
    
    void setOutputIntervalTime(double interval);
    void setOutputIntervalEddy(int interval);

    void calculateSolution(const double tRun, bool shouldWriteData =false);                         // Running simulations 

    void writeData(int real, const int ifile, const double outputTime);                                       // Writing the results for a user-defined number of eddies in the data folder.

private:

    std::vector<double> projection(std::vector<double> &vcfd, std::vector<double> &weight);         // Perform vector projection of flow particles onto hips parcels operation without density 
    
    std::pair<std::vector<double>, std::vector<double>>  projection(std::vector<double> &vcfd, std::vector<double> &weight,                     
                                   const std::vector<double> &density);                             // Perform vector projection flow particles onto hips parcels operation with density 
    
    std::vector<double> setGridHips(int N);                                                         // Set Hips grid with a specified number of grid points equal to number of parcels   
    std::vector<double> setGridCfd(std::vector<double> &w);                                         // Set CFD grid using provided weight vector
    std::vector<double>  projection_back(std::vector<double> &vb);                       
    std::pair<std::vector<double>, std::vector<double>> projection_back_with_density(std::vector<double> &vh, 
                                                                                       std::vector<double> &rho_h);
    



    void sample_hips_eddy(double &dt, int &iLevel);                                                 // Sample hips eddy with specified time step and level                                                
    void selectAndSwapTwoSubtrees(const int iLevel, int &iTree);                                    // Select and swap two subtrees in the level tree
    void advanceHips(const int iLevel, const int iTree);                                            // Advancing simulations to do mixing and reaction
   
    int getVariableIndex(const std::string &varName) const;

    void reactParcels_LevelTree(const int iLevel, const int iTree);                                 // Reacting parcels involved in micro-mixing
    void mixAcrossLevelTree(int kVar, const int iMixLevel, const int iTree);                        // Mixing paecels involved in micr0-mixing.
   
    void forceProfile();
    
    void saveAllParameters();                                                                       // Function to save ALL parameters
    
    /////////////////////////  COSTRUCTORS  ///////////////////////// 


public:

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Constructor for initializing a HiPS object without building the tree immediately.
///
/// This constructor is designed for simulations where the HiPS tree must be configured dynamically, 
/// such as grid-based simulations with cell-specific turbulence properties. The tree can be created 
/// or updated later using one of the `set_tree` functions.
/// \param C_param_         Eddy coefficient controlling mixing rate.
/// \param forceTurb_       Flag to enforce turbulence activation.
/// \param nVar_            Number of transported variables.
/// \param performReaction_ Enables chemical reactions if set to true.
/// \param cantSol          Cantera solution object (required when REACTIONS_ENABLED).
/// \param seed             Random seed (negative for random initialization).
/// \param realization_     Realization index for ensemble or parallel runs.
///
/// \note This constructor does not call any `set_tree()` function. The user is responsible 
///       for building the tree explicitly by calling either version of `set_tree(...)`.
///
/// \note See the full constructor for chemical integrator behavior under `REACTIONS_ENABLED`.
///
/// \see hips::set_tree() for deferred tree construction.
////////////////////////////////////////////////////////////////////////////////////////////////////
   
    hips(double C_param_,
         int forceTurb_,
         int nVar_,
         std::vector<double> &ScHips_,
         bool performReaction,
         std::shared_ptr<void> vcantSol = nullptr,
         int seed = 10,
         int realization_ = 1);

//////////////////////////////////////////////////////////////
/// \brief Constructor for initializing the full HiPS tree at the time of object creation.
///
/// This constructor is intended for simulations where the domain structure and turbulence 
/// properties are fixed throughout the run (e.g., standalone mixing or reacting cases).
/// The HiPS tree is initialized immediately using the provided parameters, and the 
/// structure remains constant throughout the simulation.
///
/// \param nLevels_         Number of levels in the HiPS binary tree (can be adjusted for high Sc).
/// \param domainLength_    Domain length for defining spatial scales.
/// \param tau0_            Characteristic time scale for the smallest eddy.
/// \param C_param_         Eddy coefficient controlling mixing rate.
/// \param forceTurb_       Flag to enforce turbulence activation.
/// \param nVar_            Number of transported variables.
/// \param ScHips_          Vector of Schmidt numbers (one per variable).
/// \param performReaction_ Enables chemical reactions if set to true.
/// \param cantSol          Cantera solution object (required when REACTIONS_ENABLED).
/// \param seed             Random seed (negative for random initialization).
/// \param realization_     Realization index for ensemble or parallel runs.
///
/// \note This constructor calls `set_tree(nLevels, domainLength, tau0, ScHips)` internally 
///       to fully build the tree at initialization time. This setup is optimal for cases 
///       where tree reconfiguration is not needed during runtime.
///
/// \note If `REACTIONS_ENABLED` is defined, the default chemical integrator is `batchReactor_cvode`.
///       To use `batchReactor_cantera` instead, uncomment the corresponding line in the constructor code.
///
/// \see hips::set_tree() for the internal tree setup logic.
/////////////////////////////////////////////////////////////////////////////////////////////////////////


    hips(int nLevels_,
         double domainLength_,
         double tau0_,
         double C_param_,
         int forceTurb_,
         int nVar_,
         std::vector<double> &ScHips_,
         bool performReaction,
         std::shared_ptr<void> vcantSol = nullptr,
         int seed = 10, 
         int realization_ = 1);

    void resetForNewRealization() { // Reset the number of Index to use for new realization
        currentIndex = 0;
    }
 
};
