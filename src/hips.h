#ifndef  HIPS_H
#define  HIPS_H

#include "cantera/thermo.h"
#include "cantera/base/Solution.h"
#include "cantera/numerics/Integrator.h"
#include "cvodeDriver.h"

#include <iostream>
#include <vector>
#include <string>
#include <map>

//#include "state.h"
#include "randomGenerator.h"

/////////////////////////////////////////////////////////////////////////////


class hips {


////////////////////////////// DATA MEMBERS /////////////////////////////

    public:

        int                                   nparcels; ///< number of parcels
        std::vector<std::vector<double > * >  varData;  ///< vector of pointers to vector
        std::shared_ptr<Cantera::ThermoPhase> gas;

    private:

        int                 nLevels;
        double              domainLength;        ///< length of domain (m) 
        double              tau0;                ///< integral timescale
        double              C_param;             ///< Eddy frequency parameter
        int                 forceTurb;           ///< forcing function for statistically stationary: -1 = none, 1 = source term, 2 = dir
        int                 nVar;                ///< number of parcel variables (e.g., h, ysp)
        int                 nsp;                 ///< number of species

        double              time;          
        int                 Nm1;                 ///< nLevels - 1 
        int                 Nm2;                 ///< nLevels - 2                     
        int                 Nm3;                 ///< nLevels - 3 
        int                 iEta;                ///< Kolmogorov level (needed for variable Sc scalars)
        double              eddyRate_total;      ///< total rate of all eddies 0 through nLevels-3
        double              eddyRate_inertial;   ///< tota    public:
        bool                LScHips;             ///< hips schmidt number
        std::vector<double> i_batchelor;         ///< Batchelor level for variable Sc scalars; NOTE: double, as in, between levels
        std::vector<int>    i_plus;              ///< ceil(i_batchelor)
        std::vector<double> parcelTimes;         ///< current times corresponding to the parcel states
      
        std::vector<double> levelRates;          ///< list of eddy event rates at each level
        std::vector<int>    pLoc;                ///< parcel index array for fast implementation of swaps
        std::vector<double> ScHips; 
        cvodeDriver         cvodeD;     

        bool                LrandSet;            ///< flag indicating new randomGen  --> allow deletion
        randomGenerator     *rand;

        const  int          seed = 30;           ///<  random number generator seed (negative to randomize it)
        const double        Afac = 0.5;          ///< level lengthscale reduction factor (0.5) 
               
  /////////////////////////////// MEMBER FUNCTIONS//////////////////////////////////////////

    public:
        
        hips(int             nLevels_,
             double          domainLength_,
             double          tau0_,
             double          C_param_,
             int             forceTurb_,
             int             nVar_,
             std::vector<double> &ScHips_,
             std::shared_ptr<Cantera::Solution> cantSol);

        virtual ~hips() { if(LrandSet) delete rand; }

        void calculateSolution(const double tRun);                            // Solution for hips
        void set_varData(std::vector<double> &v, int i) { varData[i] = &v; }  // passing all variables to vector of pointer     

    private:

        void sample_hips_eddy(double &dt, int &iLevel);
        void selectAndSwapTwoSubtrees(const int iLevel, int &iTree);
        void advanceHips(const int iLevel, const int iTree);
        void reactParcels_LevelTree(const int iLevel, const int iTree);
        void mixAcrossLevelTree( int kVar, const int iMixLevel, const int iTree);                       // const int iTree=-1);
        void forceProfile();
        void writeData(const int ifile, const double outputTime);

        void setState(const int &ipt);     // setting gas state in each parcel
};
#endif


