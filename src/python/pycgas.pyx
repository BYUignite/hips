from libcpp.string cimport string
from libcpp.vector cimport vector

#================================================================================

cdef extern from "../c++/cgas.h":
    cdef cppclass cgas:

        cgas(string) except +

        int nSpecies()
        double enthalpy_mass()
        double density()
        double pressure()
        void setState_TPX(double, double, string)
        void setState_HP(double, double)
        vector[double] getMassFractions()
        void equilibrate(string)
        string speciesName(int)
        vector[string] speciesNames()

#================================================================================

cdef class pycgas:

    cdef cgas *cgasptr

    #------------------

    def __cinit__(self, mechanismName):
        cdef string mechanism_name = mechanismName.encode("utf-8")

        self.cgasptr = new cgas(mechanism_name)

    #------------------

    def __dealloc__(self):
        if self.cgasptr != NULL:
            del self.cgasptr

    #------------------

    def nSpecies(self):
        return self.cgasptr.nSpecies()

    def enthalpy_mass(self):
        return self.cgasptr.enthalpy_mass()

    def density(self):
        return self.cgasptr.density()

    def pressure(self):
        return self.cgasptr.pressure()

    def setState_TPX(self, T, P, X):
        cdef string composition = X.encode("utf-8")
        self.cgasptr.setState_TPX(T, P, composition)

    def setState_HP(self, H, P):
        self.cgasptr.setState_HP(H, P)

    def getMassFractions(self):
        cdef vector[double] y = self.cgasptr.getMassFractions()
        return y

    def equilibrate(self, XY):
        cdef string equil_type = XY.encode("utf-8")
        self.cgasptr.equilibrate(equil_type)

    def speciesName(self, k):
        cdef string name = self.cgasptr.speciesName(k)
        return name.decode("utf-8")

    def speciesNames(self):
        cdef vector[string] names = self.cgasptr.speciesNames()
        return [names[k].decode("utf-8") for k in range(names.size())]
