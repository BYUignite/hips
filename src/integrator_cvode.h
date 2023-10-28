#pragma once

#include <vector>

#include <cvode/cvode.h>               // prototypes for CVODE fcts., consts.
#include <nvector/nvector_serial.h>    // access to serial N_Vector
#include <sunmatrix/sunmatrix_dense.h> // access to dense SUNMatrix
#include <sunlinsol/sunlinsol_dense.h> // access to dense SUNLinearSolver

////////////////////////////////////////////////////////////////////////////////

class integrator_cvode {

public: 

//////////////////// DATA MEMBERS ////////////////////

    SUNContext      sun;           // sundials object
    void           *cmem;         // cvode object
    N_Vector        vars;             // vector of variables being solved
    N_Vector        atol;          // vector atol (for each variable)
    realtype        rtol;          // scalar rtol
    unsigned        nvar;           // number of equations being solved
    SUNMatrix       J;             // matrix for linear solver
    SUNLinearSolver LS;            // linear solver
    int             rv;            // return value: checking status of calls

//////////////////// MEMBER FUNCTIONS ////////////////////

integrator_cvode(
                 int (*Func)(realtype, N_Vector, N_Vector, void*),
                 void *                _user_data,
                 const int             _nvar, 
                 const double          _rtol,
                 const std::vector<double> &_atol) :
    nvar(_nvar),
    rtol(_rtol) {

    rv = SUNContext_Create(NULL, &sun);

    vars = N_VNew_Serial(nvar, sun);
    atol = N_VNew_Serial(nvar, sun);
    for(int k=0; k<nvar; ++k)
        NV_Ith_S(atol, k) = _atol[k];

    cmem = CVodeCreate(CV_BDF, sun);
    rv   = CVodeSetUserData(cmem, _user_data);
    rv   = CVodeInit(cmem, Func, 0.0, vars);
    rv   = CVodeSVtolerances(cmem, rtol, atol);
    rv   = CVodeSetMaxNumSteps(cmem, 5000);

    J    = SUNDenseMatrix(nvar, nvar, sun);   // linear solver matrix J
    LS   = SUNLinSol_Dense(vars, J, sun);          // set linear solver
    rv   = CVodeSetLinearSolver(cmem, LS, J);  // associate matrix J and solver LS
}

//--------------

int integrate(std::vector<double> &y, const realtype dt) {

    for(int k=0; k<nvar; ++k)
        NV_Ith_S(vars, k) = y[k];

    realtype t;

    rv = CVodeReInit(cmem, 0.0, vars);
    rv = CVode(cmem, dt, vars, &t, CV_NORMAL);

    for(int k=0; k<nvar; ++k)
        y[k] = NV_Ith_S(vars,k);

    return rv;
}

//--------------

~integrator_cvode() {

    N_VDestroy(vars);
    N_VDestroy(atol);
    CVodeFree(&cmem);
    SUNLinSolFree(LS);
    SUNMatDestroy(J);
    SUNContext_Free(&sun);
}
};
