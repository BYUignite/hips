#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "hips.h"

#include <memory>
#include <cmath>

using namespace std;
using Catch::Matchers::WithinRel;

//////////////////////////////////////////////////////////////////////////

TEST_CASE( "Test HiPS library" ) {

    int            nLevels = 9;                  // Number of hierarchical levels
    double         domainLength = 1.0;           // Simulation domain length
    double         tau0 = 1.0;                   // Initial time scale for the largest eddies
    double         C_param = 0.5;                // Eddy turnover rate multiplier
    double         tRun = 300.0;                 // Total simulation time

    //////////////////////////////////////////////////////////////////////////

    SECTION("Test variable setup: projection and back projection") {

        vector<double> ScHips = {1.0};    // Schmidt numbers for low and high diffusivity
        int            nvars = 1;         // Number of scalar fields
        bool           forceTurb = true;  // Forcing parameter to impose turbulent profile (not applicable here)

        hips H(nLevels, domainLength, tau0, C_param, forceTurb, nvars, ScHips, false);

        vector<double> var = {0.985, 0.784, 0.508, 0.404, 0.913, 0.126, 0.561};
        vector<double> w   = {0.422, 6.845, 0.851, 7.806, 3.512, 4.258, 2.071};   // test: not a power of 2, not uniform, don't sum to 1

        //--------- compute raw input sum

        double sumw = 0.0;                                  // normalize w
        for(int i=0; i<var.size(); i++) sumw += w[i];
        for(int i=0; i<var.size(); i++) w[i] /= sumw;
        double sum1 = 0.0;
        for(int i=0; i<var.size(); i++) sum1 += var[i]*w[i];

        //--------- compute sum of HiPS projected variable

        H.set_varData(var, w, "test");
        auto pLoc    = H.get_pLoc();
        auto varData = H.get_HipsVarData_ptr();

        double sum2 = 0.0;
        for(int i=0; i<H.get_nparcels(); i++)
            sum2 += (*varData[0])[pLoc[i]];
        sum2 /= H.get_nparcels();

        //------------ compute back_projected sum

        auto var2 = H.get_varData();
        double sum3 = 0.0;
        for(int i=0; i<var2[0].size(); i++)
            sum3 += var2[0][i]*w[i];

        //------------ test sum1=sum2 and sum2=sum3

        REQUIRE_THAT( sum1, WithinRel(sum2, 1E-14) );
        REQUIRE_THAT( sum2, WithinRel(sum3, 1E-14) );

    }

    //////////////////////////////////////////////////////////////////////////

    SECTION("Test overall solve with forcing") {

        vector<double> ScHips = {0.0625, 1.0, 16.0}; // Schmidt numbers for low and high diffusivity
        int            nvars = 3;                    // Number of scalar fields
        bool           forceTurb = true;             // Forcing parameter to impose turbulent profile

        hips H(nLevels, domainLength, tau0, C_param, forceTurb, nvars, ScHips, false);

        //--------- initialize vars

        vector<vector<double>> vars(nvars, vector<double>(H.get_nparcels()));
        vector<double> weights(H.get_nparcels(), 1.0/H.get_nparcels());  // Uniform weights

        for (int k=0; k<nvars; k++) {
            for (int i=0; i<H.get_nparcels(); i++)
                vars[k][i] = (i<H.get_nparcels()/2) ? 0 : 1;
            H.set_varData(vars[k], weights, "mixf_0" + to_string(k));
        }
        
        //--------- inititial raw conservation

        vector<double> sum1(nvars, 0.0);
        for(int i=0; i<H.get_nparcels(); i++)
            for(int k=0; k<nvars; k++)
                sum1[k] += vars[k][i]*weights[i];

        //--------- solve

        H.calculateSolution(tRun, false);

        //--------- recover variables, final conservation

        auto var2 = H.get_varData();

        vector<double> sum2(var2.size(), 0.0);
        for(int i=0; i<var2[0].size(); i++)
            for(int k=0; k<var2.size(); k++)
                sum2[k] += var2[k][i]*weights[i];

        //--------- test

        REQUIRE( H.get_nparcels() == (1 << (nLevels+1)) );        // based on ScHips set above
        REQUIRE_THAT(sum1[0],      WithinRel(sum2[0], 1E-14));
        REQUIRE_THAT(sum1[1],      WithinRel(sum2[1], 1E-14));
        REQUIRE_THAT(sum1[2],      WithinRel(sum2[2], 1E-14));
        REQUIRE_THAT(var2[0][849], WithinRel(1.7783400336213617, 1E-14));
        REQUIRE_THAT(var2[1][849], WithinRel(1.6986828219520567, 1E-14)); 
        REQUIRE_THAT(var2[2][849], WithinRel(1.2389359078700295, 1E-14));

    }

    //////////////////////////////////////////////////////////////////////////

}
