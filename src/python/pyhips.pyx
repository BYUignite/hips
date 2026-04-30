##
# @file pyrad.pyx
# \brief cython file for Python interface to c++ code
#

from libcpp cimport bool
from libcpp.memory  cimport shared_ptr
from libcpp.string  cimport string
from libcpp.utility cimport pair
from libcpp.vector  cimport vector
from cython.operator cimport dereference as deref

import numpy as np

#================================================================================

cdef extern from "../c++/hips.h":
    cdef cppclass hips:

        hips(double, bool, int, vector[double]&,
             bool, shared_ptr[void], int, int) except +

        hips(int, double, double, double, bool, int,
             vector[double]&, bool, shared_ptr[void], int, int) except +

        void set_tree(int, double, double)
        void set_tree(double, double, double, string)

        void set_varData(vector[double]&, vector[double]&, const string&)
        void set_varData(vector[double]&, vector[double]&,
                         const string&, const vector[double]&)

        vector[vector[double]] get_varData()
        pair[vector[vector[double]], vector[double]] get_varData_with_density()

        void setOutputIntervalTime(double)
        void setOutputIntervalEddy(int)

        void calculateSolution(double, bool)
        void writeData(int, int, double)
        int get_nparcels() const

        const vector[int]& get_pLoc() const
        const vector[shared_ptr[vector[double]]]& get_HipsVarData_ptr() const

#================================================================================

cdef vector[double] _as_vector_double(object values):
    cdef vector[double] out
    cdef object arr
    cdef Py_ssize_t i, n

    arr = np.asarray(values, dtype=np.float64).ravel()
    n = arr.shape[0]
    out.reserve(n)
    for i in range(n):
        out.push_back(arr[i])
    return out

#------------------

cdef vector[vector[double]] _copy_vector_vector_double(vector[vector[double]] values):
    return values

#------------------

cdef class pyhips:

    cdef hips *hipsptr

    #------------------

    def __cinit__(self,
                  double C_param,
                  bint forceTurb,
                  int nVar,
                  ScHips,
                  bint performReaction,
                  vcantSol=None,
                  int seed=10,
                  int realization=1,
                  nLevels=None,
                  domainLength=None,
                  tau0=None):
        cdef vector[double] sc_hips = _as_vector_double(ScHips)

        self.hipsptr = NULL
        if vcantSol is not None:
            raise NotImplementedError("Python wrapper does not yet support vcantSol")

        if nLevels is None and domainLength is None and tau0 is None:
            self.hipsptr = new hips(C_param, forceTurb, nVar, sc_hips,
                                    performReaction, shared_ptr[void](),
                                    seed, realization)
        elif nLevels is not None and domainLength is not None and tau0 is not None:
            self.hipsptr = new hips(nLevels, domainLength, tau0, C_param,
                                    forceTurb, nVar, sc_hips,
                                    performReaction, shared_ptr[void](),
                                    seed, realization)
        else:
            raise ValueError("Provide either all of nLevels/domainLength/tau0 or none of them")

    #------------------

    def __dealloc__(self):
        if self.hipsptr != NULL:
            del self.hipsptr

    #------------------

    def set_tree(self, level_or_re, double domainLength, double tau0, ReApproach="rounding"):
        cdef string re_approach
        if isinstance(level_or_re, int):
            self.hipsptr.set_tree(level_or_re, domainLength, tau0)
        else:
            re_approach = ReApproach.encode("utf-8")
            self.hipsptr.set_tree(level_or_re, domainLength, tau0, re_approach)

    #------------------

    def set_varData(self, v, w, varN, rho=None):
        cdef vector[double] v_vec = _as_vector_double(v)
        cdef vector[double] w_vec = _as_vector_double(w)
        cdef vector[double] rho_vec
        cdef string var_name = varN.encode("utf-8")

        if rho is None:
            self.hipsptr.set_varData(v_vec, w_vec, var_name)
        else:
            rho_vec = _as_vector_double(rho)
            self.hipsptr.set_varData(v_vec, w_vec, var_name, rho_vec)

    #------------------

    def get_varData(self):
        cdef vector[vector[double]] data = self.hipsptr.get_varData()
        cdef Py_ssize_t i

        return [np.array(data[i]) for i in range(data.size())]

    #------------------

    def get_varData_with_density(self):
        cdef pair[vector[vector[double]], vector[double]] data
        cdef Py_ssize_t i

        data = self.hipsptr.get_varData_with_density()
        return [np.array(data.first[i]) for i in range(data.first.size())], np.array(data.second)

    #------------------

    def setOutputIntervalTime(self, double interval):
        self.hipsptr.setOutputIntervalTime(interval)

    #------------------

    def setOutputIntervalEddy(self, int interval):
        self.hipsptr.setOutputIntervalEddy(interval)

    #------------------

    def calculateSolution(self, double tRun, bint shouldWriteData=False):
        self.hipsptr.calculateSolution(tRun, shouldWriteData)

    #------------------

    def writeData(self, int real, int ifile, double outputTime):
        self.hipsptr.writeData(real, ifile, outputTime)

    #------------------

    def get_nparcels(self):
        return self.hipsptr.get_nparcels()

    #------------------

    def get_pLoc(self):
        cdef const vector[int] *p_loc = &self.hipsptr.get_pLoc()
        cdef Py_ssize_t i

        return np.array([deref(p_loc)[i] for i in range(deref(p_loc).size())], dtype=np.int32)

    #------------------

    def get_HipsVarData_ptr(self):
        cdef const vector[shared_ptr[vector[double]]] *var_data = &self.hipsptr.get_HipsVarData_ptr()
        cdef Py_ssize_t i

        return [np.array(deref(deref(var_data)[i])) for i in range(deref(var_data).size())]

#================================================================================
