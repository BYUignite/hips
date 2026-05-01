# execution: python ex_2.py
# make sure the pyhips library path is in your PYTHONPATH environment variable
# make sure HiPS was built with reaction support and cantera is installed

###############################################################################

import cantera as ct
import pyhips
import numpy as np
import pickle

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

gas = ct.Solution(mechanismName)
nsp = gas.n_species
nVar = nsp + 1
ScHips = np.ones(nVar)

#--------------------------- setup the tree

#doldb hips = pyhips.pyhips(C_param, forceTurb, nVar, ScHips, doReaction, \
#doldb                     mechanismName, seed, realization, nLevels, domainLength, tau0)

#--------------------------- set the data

nparcels = 128 #hips.get_nparcels()
fracBurn = 0.25

ysp = np.zeros((nsp, nparcels))
h = np.zeros(nparcels)
rho = np.zeros(nparcels)

T0 = 300.0
T1 = 300.0
composition = "C2H4:1, O2:3, N2:11.25"

#---------- initialize fresh reactants

gas.TPX = T0, ct.one_atm, composition

freshEnd = int((1.0 - fracBurn) * nparcels)
h[:freshEnd+1] = gas.enthalpy_mass
rho[:freshEnd+1] =gas.density
for i in range(nsp):
    ysp[i, :freshEnd+1] = gas.Y[i]

#---------- initialize combusted parcels

gas.TPX = T1, ct.one_atm, composition
gas.equilibrate("HP")

h[freshEnd + 1:] = gas.enthalpy_mass
rho[freshEnd + 1:] = gas.density
for i in range(nsp):
    ysp[i, freshEnd+1:] = gas.Y[i]

variableNames = ["enthalpy"] + gas.species_names
w = np.ones(nparcels) / nparcels

##############################################################

with open('data.pkl', 'wb') as f:
    pickle.dump(h, f)
    pickle.dump(rho, f)
    pickle.dump(ysp, f)
    pickle.dump(variableNames, f)
    pickle.dump(w, f)

##############################################################

#doldb hips.set_varData(h, w, variableNames[0], rho)
#doldb for k in range(nsp):
#doldb     hips.set_varData(ysp[k], w, variableNames[k+1], rho)
#doldb 
#doldb #--------------------------- initial condition
#doldb 
#doldb print("Saving files to ../data/")
#doldb 
#doldb hips.writeData(1, 0, 0.0)
#doldb 
#doldb #--------------------------- run solution
#doldb 
#doldb hips.setOutputIntervalTime(tRun/10)
#doldb 
#doldb doWrite = True
#doldb hips.calculateSolution(tRun, doWrite)
#doldb 
#doldb print(len(hips.get_HipsVarData_ptr()))
