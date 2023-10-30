#include "hips.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

///////////////////////////////////////////////////////////////////////////////

int main() {
  
    //---------- setting  HiPS tree and constructor

    int            nLevels      = 9;           // # HiPS tree levels
    double         domainLength = 1.0;         // HiPS tree lengthscale
    double         tau0         = 1.0;         // HiPS tree timescale
    double         C_param      = 0.5;         // HiPS eddy rate multiplier
    double         tRun         = 400.0;       // tree evolution time
    int            forceTurb    = 0;           // Force turbulent profile
    vector<double> ScHips(7,1);

    //---------- setting gas solution

    //auto   cantSol = Cantera::newSolution("2S_CH4_BFER.yaml", "2S_CH4_BFER", "None");
    auto   cantSol = Cantera::newSolution("gri30.yaml");
    auto   gas = cantSol->thermo();
    size_t nsp = gas->nSpecies();
    int    nVar = nsp+1;           // 1 + nsp: h, y 

    hips HiPS(nLevels,
              domainLength,
              tau0, 
              C_param,
              forceTurb,
              nVar,
              ScHips,
              cantSol, 
              nullptr,
              true);

    int nparcels = HiPS.nparcels;

    //------------- declaring state vectors

    vector<vector <double > > ysp(nsp, vector<double> (nparcels, 0));
    vector<double>            h(nparcels); 

    vector<double>            y0(nsp);
    double                    T0;
    double                    h0;

    vector<double>            y1(nsp);
    double                    T1;
    double                    h1;

    //-------------- set state vectors in each parcel

    T0 = 300.0;
    gas->setState_TPX(T0, Cantera::OneAtm, "CH4:1");
    gas->getMassFractions(&y0[0]);
    h0 = gas->enthalpy_mass();

    T1 = 300;
    gas->setState_TPX(T1, Cantera::OneAtm, "O2:2, N2:7.52");
    gas->getMassFractions(&y1[0]);
    h1 = gas->enthalpy_mass();

    double mixf = 16/(16+29*9.52);  // stoichiometric
    vector<double> ymix(nsp);
    for(int k=0; k<nsp; k++)
        ymix[k] = y0[k]*(1-mixf) + y1[k]*mixf;
    double hmix = h0*(1-mixf) + h1*mixf;
    vector<double> yeq(nsp);
    gas->setMassFractions(&ymix[0]);
    gas->setState_HP(hmix, gas->pressure());
    gas->equilibrate("HP");

    double fracBurn = 0.5;
    for(int i=0; i<fracBurn*nparcels; i++) {
        h[i] = hmix;
        for (int k=0; k<nsp; k++)
            ysp[k][i] = gas->massFraction(k);
    }
    for(int i=fracBurn*nparcels; i<nparcels; i++) {
        h[i] = hmix;
        for (int k=0; k<nsp; k++)
            ysp[k][i] = ymix[k];
    }

    //----------------- set array of hips variables from state variables

    HiPS.set_varData(h, 0);
    for (int k=0; k<ysp.size(); k++) 
        HiPS.set_varData(ysp[k],k+1);

    //----------------- advancing HiPS to do rxn and mixing

    HiPS.calculateSolution(tRun);

    return 0;

}  
