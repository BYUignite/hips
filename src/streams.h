#pragma once

#include "cantera/base/Solution.h"
#include "cantera/thermo.h"

#include <vector>
#include <memory>

////////////////////////////////////////////////////////////////////////////////

/** Class implementing streams for use in mixing and or reaction problems.
 *  This is writting in terms of mixture fraction with streams defined in an
 *  input file.  The class can implement products of complete combustion, or
 *  equilibrium (through the Cantera IdealGasPhase object, if desired).
 *  This class holds a pointer to a Cantera IdealGasPhase object (defined up front
 *  in main) that computes thermodynamic, kinetic, and transport data.
 *
 *  @author David O. Lignell
 */

class streams {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        double              h0;              ///< stream mixf=0 enthalpy (J/kg)
        double              h1;              ///< stream mixf=1 enthalpy (J/kg)
        std::vector<double> y0;              ///< stream mixf=0 composition vector
        std::vector<double> y1;              ///< stream mixf=1 composition vector

        double P;

        std::shared_ptr<Cantera::ThermoPhase> gas;

        double mixfStoic;       ///< stoichiometric mixture fraction
        int    nspc;            ///< number of species in gas mechanism
        double beta0;           ///< mixf = (beta-beta0) / (beta1-beta0)
        double beta1;           ///< mixf = (beta-beta0) / (beta1-beta0)

        std::vector<double> gCHON;           ///< gammas, as in beta = sum_i (y_i*gamma_i)

    //////////////////// MEMBER FUNCTIONS /////////////////

        void getProdOfCompleteComb(const double mixf,
                                   std::vector<double> &ypcc,
                                   double &hpcc,
                                   double &Tpcc);

        void getEquilibrium_HP(const double mixf,
                                     std::vector<double> &yeq,
                                     double &heq,
                                     double &Teq);

        void getEquilibrium_TP(const double mixf,
                                     double Teq,
                                     std::vector<double> &yeq,
                                     double &heq);

        void getMixingState(const double mixf,
                            std::vector<double> &ymix,
                            double &hmix,
                            double &Tmix);

        double getMixtureFraction(const double *y,
                                  const bool doBeta01=false);

    private:

        void setStoicMixf();
        std::vector<double> setElementMassFracs(const double *y);
        std::vector<double> setElementMoleFracs(const double *y);
        std::vector<double> getElementMoles(const double *x,
                                       double &nOnotFromO2,
                                       double &nHnotFromH2O,
                                       double &nCnotFromCO2);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        streams() {}
        streams(std::shared_ptr<Cantera::Solution> csol,
                const double _P,
                const double _h0, 
                const double _h1,
                const std::vector<double> &_y0, 
                const std::vector<double> &_y1);

        ~streams(){}

};



