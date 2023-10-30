#include "hips.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

/**
 * Calculates vh using the provided project_cfd_to_hips function.
 * @param xc The vector of cfd particle positions.
 * @param vc The vector of cfd particle values.
 * @param xh The vector of hips particle positions.
 * @return The calculated vh vector.
 */
vector<double> projectCfdonHiPs(const vector<double>& xc, const vector<double>& vc, const vector<double>& xh) {
    int jprev = 0;
    int nh = xh.size() - 1;
    int nc = xc.size() - 1;
    vector<double> vh(nh, 0.0);

    for (int i = 0; i < nh; ++i) {
        for (int j = jprev + 1; j <= nc; ++j) {
            if (xc[j] <= xh[i + 1]) {
                double Δ1 = xc[j] - xc[j - 1];
                double Δ2 = xc[j] - xh[i];
                double Δ = std::min(Δ1, Δ2);
                vh[i] += vc[j - 1] * Δ;
            } else {
                double Δ1 = xh[i + 1] - xc[j - 1];
                double Δ2 = xh[i + 1] - xh[i];
                double Δ = std::min(Δ1, Δ2);
                vh[i] += vc[j - 1] * Δ;
                jprev = j - 1;
                break;
            }
        }
        vh[i] /= (xh[i + 1] - xh[i]);
    }

    return vh;
}

/**
 * Calculates vc using the provided project_hips_to_cfd function.
 * @param vh The vector of vh values.
 * @param xh The vector of hips particle positions.
 * @param xc The vector of cfd particle positions.
 * @return The calculated vc vector.
 */
vector<double> projectHipsToCfd(const vector<double>& vh, const vector<double>& xh, const vector<double>& xc) {
    int jprev = 0;
    int nc = xc.size() - 1;
    vector<double> vc(nc, 0.0);

    for (int i = 0; i < nc; ++i) {
        for (int j = jprev + 1; j < xh.size(); ++j) {
            double Δ1 = xh[j] - xh[j - 1];
            double Δ2 = xh[j] - xc[i];
            double Δ = (Δ1 < Δ2) ? Δ1 : Δ2;

            if (xh[j] <= xc[i + 1]) {
                vc[i] += vh[j - 1] * Δ;
            } else {
                double Δ1 = xc[i + 1] - xh[j - 1];
                double Δ2 = xc[i + 1] - xc[i];
                Δ = (Δ1 < Δ2) ? Δ1 : Δ2;
                vc[i] += vh[j - 1] * Δ;
                jprev = j - 1;
                break;
            }
        }
        vc[i] /= (xc[i + 1] - xc[i]);
    }

    return vc;
}

int main() {
    
    vector<double> xc = {0.0, 0.1,0.2, 0.4, 0.6, 1.0};
    vector<double> vc = {1, 0.5, 1.7, 2, 1.5 };

    vector<double> xh(17, 0.0);
    for(int i = 0; i <= 16; i++) {
        xh[i] = i / 16.0;
    }

    vector<double> vh = projectCfdonHiPs(xc, vc, xh);

    // HiPS tree and constructor parameters
    int numLevels = log2(xh.size()) + 1; // Calculate numLevels based on xh
    double domainLength = 1.0;
    double tau0 = 1.0;
    double C_param = 0.5;
    double tRun = 20.0;
    int forceTurb = 0;
    vector<double> ScHips(1, 1);

    // Gas solution setup
    auto cantSol = Cantera::newSolution("gri30.yaml");
    auto gas = cantSol->thermo();
    size_t numSpecies = gas->nSpecies();


    int numVariables = 1;
    vector<double> modifiedData;

    // HiPS tree creation
    hips HiPS(numLevels, domainLength, tau0, C_param, forceTurb, numVariables, ScHips, cantSol, nullptr, false);

    // Set state vectors in each parcel using vh
    HiPS.set_varData(vh, 0);

    // Advance HiPS for mixing
    HiPS.calculateSolution(tRun);

    // retrieving data
    HiPS.get_varData(modifiedData, 0);

    
    return 0;
}
