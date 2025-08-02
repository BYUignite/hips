#include <catch2/catch_test_macros.hpp>
#include "hips.h"

#include <memory>
#include <cmath>

using namespace std;

//////////////////////////////////////////////////////////////////////////

TEST_CASE( "Test HiPS library" ) {

    int            nLevels = 9;                  // Number of hierarchical levels
    double         domainLength = 1.0;           // Simulation domain length
    double         tau0 = 1.0;                   // Initial time scale for the largest eddies
    double         C_param = 0.5;                // Eddy turnover rate multiplier
    double         tRun = 300.0;                 // Total simulation time
    int            forceTurb = 2;                // Forcing parameter to impose turbulent profile

    //////////////////////////////////////////////////////////////////////////

    SECTION("Test variable setup: projection and back projection") {

        vector<double> ScHips = {1.0}; // Schmidt numbers for low and high diffusivity
        int            nvars = 1;      // Number of scalar fields

        hips H(nLevels, domainLength, tau0, C_param, forceTurb, nvars, ScHips, false);

        vector<double> var = {0.985, 0.784, 0.508, 0.404, 0.913, 0.126, 0.561};
        vector<double> w   = {0.422, 6.845, 0.851, 7.806, 3.512, 4.258, 2.071};   // test: not a power of 2, not uniform, don't sum to 1

        double sumw = 0.0;                  // normalize w
        for(int i=0; i<var.size(); i++)
            sumw += w[i];
        for(int i=0; i<var.size(); i++)
            w[i] /= sumw;
        double sum1 = 0.0;                  // compute raw input sum
        for(int i=0; i<var.size(); i++)
            sum1 += var[i]*w[i];

        H.set_varData(var, w, "test");
        double sum2 = 0.0;                  // compute HiPS projected sum
        for(int i=0; i<H.nparcels; i++)
            sum2 += (*H.varData[0])[H.pLoc[i]];
        sum2 /= H.nparcels;

        REQUIRE( abs((sum1 - sum2)/sum1) < 1E-14 );   // results equal to within roundoff error

        //------------

        auto var2 = H.get_varData();
        double sum3 = 0.0;                 // compute back_projected sum
        for(int i=0; i<var2[0].size(); i++)
            sum3 += var2[0][i]*w[i];

        REQUIRE( abs((sum3 - sum2)/sum3) < 1E-14 );   // results equal to within roundoff error
    }

    //////////////////////////////////////////////////////////////////////////

    SECTION("Test overall solve") {

        vector<double> ScHips = {0.0625, 1.0, 16.0}; // Schmidt numbers for low and high diffusivity
        int            nvars = 3;                    // Number of scalar fields

        hips H(nLevels, domainLength, tau0, C_param, forceTurb, nvars, ScHips, false);

        //--------- initialize vars

        vector<vector<double>> vars(nvars, vector<double>(H.nparcels));
        vector<double> weights(H.nparcels, 1.0/H.nparcels);  // Uniform weights

        for (int k=0; k<nvars; k++) {
            for (int i=0; i<H.nparcels; i++)
                vars[k][i] = (i<H.nparcels) ? 0 : 1;
            H.set_varData(vars[k], weights, "mixf_0" + to_string(k));
        }
        //--------- solve

        H.calculateSolution(tRun, false);

        //--------- test

        REQUIRE( H.nparcels == (1 << (nLevels+1)) );                     // based on ScHips set above
        REQUIRE( abs((*H.varData[0])[H.pLoc[580]] - 0.3863866554595461) < 1E-14 );   // for given setup
        REQUIRE( abs((*H.varData[1])[H.pLoc[580]] - 1.0837252030932285) < 1E-14 );
        REQUIRE( abs((*H.varData[2])[H.pLoc[580]] - 0.9954210282686337) < 1E-14 );
    }
}
