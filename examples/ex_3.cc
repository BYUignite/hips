
/// \file ex_3.cc
/// \brief Example showing possible structure of HiPS as a subgrid model
/// 
/// ## Compilation:
/// ```bash
/// cd build
/// cmake ..
/// make
/// sudo make install
/// ```
/// 
/// ## Execution:
/// ```bash
/// cd ../run
/// ./ex_3.x
/// ```

///////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include <memory>
#include "hips.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////
int main() {

    double C_param = 0.5;             // Eddy turnover rate multiplier
    bool forceTurb = true;            // Forcing parameter to impose turbulent profile
    int nVar = 1;                     // Number of scalar fields
    vector<double> ScHips = {1.0};    // Schmidt numbers for low and high diffusivity
    bool performReaction = false;

    hips HiPS(C_param, forceTurb, nVar, ScHips, performReaction);

    //--------- initialize pseudo CFD variable

    int nCFD = 4;                          // number of CFD grid cells

    vector<vector<double>> var(nCFD);      // var[iCell][iParticle]
    var[0] = vector<double>(20, 0.0);      // 20 particles in cell 0
    var[1] = vector<double>(40, 0.0);      // 40 particles in cell 1
    var[2] = vector<double>(60, 0.0);      // 60 particles in cell 2
    var[3] = vector<double>(80, 0.0);      // 80 particles in cell 3

    vector<vector<double>> w(nCFD);
    w[0] = vector<double>(20, 1.0/20);     // uniform particle weights in cell 0
    w[1] = vector<double>(40, 1.0/40);     // uniform particle weights in cell 1
    w[2] = vector<double>(60, 1.0/60);     // uniform particle weights in cell 2
    w[3] = vector<double>(80, 1.0/80);     // uniform particle weights in cell 3

    for(int i=0; i<nCFD; i++)                    // initialize variable values
        for(int k=var[i].size()/2; k<var[i].size(); k++)
            var[i][k] = 1.0;

    /////////////// USER CFD steps...
    
    /////////////// do HiPS mixing in each CFD cell:

    for(int iCell=0; iCell<nCFD; iCell++) {

        cout << endl << "running HiPS in cell " << iCell << endl;

        //---------- reset tree

        double tau0 = 1.0;                         // time scale for the largest eddies
        int nLevels = 9;                           // Number of levels
        double domainLength = 1.0;                 // domain length scale

        HiPS.set_tree(nLevels, domainLength, tau0);

        //---------- set variables (CFD particles to HiPS parcels)

        HiPS.set_varData(var[iCell], w[iCell], "mixf");
        
        //---------- advance hips

        double tRun = 300.0;
        HiPS.calculateSolution(tRun);

        //---------- get variables (HiPS parcels to CFD particles)

        var[iCell] = HiPS.get_varData()[0];
    }

    /////////////// USER CFD steps ...

    return 0;
}
