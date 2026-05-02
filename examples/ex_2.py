# execution: python ex_2.py
# make sure the pyhips library path is in your PYTHONPATH environment variable
# make sure HiPS was built with reaction support and cantera is installed

###############################################################################

import pyhips
import pycgas
import numpy as np

#--------------------------- set parameters

nLevels = 8
domainLength = 0.01
tau0 = 0.00005
C_param = 0.5
tRun = 0.002
forceTurb = False
doReaction = True
seed = 11
realization = 1
mechanismName = "gri30.yaml"
P = 101325.

gas = pycgas.pycgas(mechanismName)
nsp = gas.nSpecies()
nVar = nsp + 1
ScHips = np.ones(nVar)


#--------------------------- setup the tree

hips = pyhips.pyhips(C_param, forceTurb, nVar, ScHips, doReaction, \
                     mechanismName, P, seed, realization, nLevels, domainLength, tau0)

#--------------------------- set the data

nparcels = hips.get_nparcels()
fracBurn = 0.25

ysp = np.zeros((nsp, nparcels))
h = np.zeros(nparcels)
rho = np.zeros(nparcels)

T0 = 300.0
T1 = 300.0
composition = "C2H4:1, O2:3, N2:11.25"

#---------- initialize fresh reactants

gas.setState_TPX(T0, P, composition)

freshEnd = int((1.0 - fracBurn) * nparcels)
h[:freshEnd+1] = gas.enthalpy_mass()
rho[:freshEnd+1] = gas.density()
for i in range(nsp):
    ysp[i, :freshEnd+1] = gas.getMassFractions()[i]

#---------- initialize combusted parcels

gas.setState_TPX(T1, P, composition)
gas.equilibrate("HP")

h[freshEnd + 1:] = gas.enthalpy_mass()
rho[freshEnd + 1:] = gas.density()
for i in range(nsp):
    ysp[i, freshEnd+1:] = gas.getMassFractions()[i]

variableNames = ["enthalpy"] + gas.speciesNames()
w = np.ones(nparcels) / nparcels

##############################################################

hips.set_varData(h, w, variableNames[0], rho)
for k in range(nsp):
    hips.set_varData(ysp[k], w, variableNames[k+1], rho)

#--------------------------- initial condition

print("Saving files to ../data/")

hips.writeData(1, 0, 0.0)

#--------------------------- run solution

hips.setOutputIntervalTime(tRun/10)

doWrite = True
hips.calculateSolution(tRun, doWrite)

print(len(hips.get_HipsVarData_ptr()))
