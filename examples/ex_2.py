# execution: python ex_2.py
# make sure the pyhips library path is in your PYTHONPATH environment variable
# make sure HiPS was built with reaction support and cantera is installed

###############################################################################

#doldb import cantera as ct
import pyhips
import numpy as np

import pickle

with open('data.pkl', 'rb') as f:
    h = pickle.load(f)
    rho = pickle.load(f)
    ysp = pickle.load(f)
    variableNames = pickle.load(f)
    w = pickle.load(f)

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

#doldb gas = ct.Solution(mechanismName)
nsp = 53 #doldb gas.n_species
nVar = nsp + 1
ScHips = np.ones(nVar)

#--------------------------- setup the tree

hips = pyhips.pyhips(C_param, forceTurb, nVar, ScHips, doReaction, \
                     mechanismName, P, seed, realization, nLevels, domainLength, tau0)

#--------------------------- set the data

nparcels = 128 #hips.get_nparcels()
fracBurn = 0.25

#doldb ysp = np.zeros((nsp, nparcels))
#doldb h = np.zeros(nparcels)
#doldb rho = np.zeros(nparcels)

T0 = 300.0
T1 = 300.0
composition = "C2H4:1, O2:3, N2:11.25"

#doldb #---------- initialize fresh reactants
#doldb 
#doldb gas.TPX = T0, ct.one_atm, composition
#doldb 
#doldb freshEnd = int((1.0 - fracBurn) * nparcels)
#doldb h[:freshEnd+1] = gas.enthalpy_mass
#doldb rho[:freshEnd+1] =gas.density
#doldb for i in range(nsp):
#doldb     ysp[i, :freshEnd+1] = gas.Y[i]
#doldb 
#doldb #---------- initialize combusted parcels
#doldb 
#doldb gas.TPX = T1, ct.one_atm, composition
#doldb gas.equilibrate("HP")
#doldb 
#doldb h[freshEnd + 1:] = gas.enthalpy_mass
#doldb rho[freshEnd + 1:] = gas.density
#doldb for i in range(nsp):
#doldb     ysp[i, freshEnd+1:] = gas.Y[i]
#doldb 
#doldb variableNames = ["enthalpy"] + gas.species_names
#doldb w = np.ones(nparcels) / nparcels

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
