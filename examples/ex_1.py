import pyhips
import numpy as np

#--------------------------- set parameters

nLevels = 9
domainLength = 1.0
tau0 = 1.0
C_param = 0.5
tRun = 300.0
forceTurb = 2
ScHips = np.array([0.065, 1.0, 16.0])
numVariables = 3
doReaction = False
writeData = True
seed = 10
realization = 1
vcantSol = None

#--------------------------- setup the tree

# Option 1 of 2
hips = pyhips.pyhips(C_param, forceTurb, numVariables, ScHips, doReaction, vcantSol, seed, realization)
hips.set_tree(nLevels, domainLength, tau0)

# Option 2 of 2 (same results)
#hips = pyhips.pyhips(C_param, forceTurb, numVariables, ScHips, doReaction, \
#                     vcantSol, seed, realization, nLevels, domainLength, tau0)

#--------------------------- set the data

mixf = np.zeros(hips.get_nparcels()); 
mixf[hips.get_nparcels() // 2:] = 1.0
w = np.ones(hips.get_nparcels())/hips.get_nparcels()

hips.set_varData(mixf, w, "mixf_00")
hips.set_varData(mixf, w, "mixf_01")
hips.set_varData(mixf, w, "mixf_02")

#--------------------------- initial condition

print("Saving files to ../data/")

hips.writeData(1, 0, 0.0)

#--------------------------- run solution

hips.setOutputIntervalTime(60.0)

hips.calculateSolution(tRun, writeData)

#=========================== recover the data (projected)

v = hips.get_varData()      # list of numpy arrays: v[0] is a numpy array of the first variable, etc.

np.savetxt("../data/v0.dat", v[0])  # will save variable 0 to file v0.dat
np.savetxt("../data/v1.dat", v[1])  # will save variable 1 to file v1.dat
np.savetxt("../data/v2.dat", v[2])  # will save variable 2 to file v2.dat

#---------- data (projected), along with density

v, rho = hips.get_varData_with_density()   # v is as for hips.get_varData above; rho is density

#---------- raw hips data (not projected)

vraw = hips.get_HipsVarData_ptr()

#---------- pLoc

pLoc = hips.get_pLoc()                     # integer index array (numpy) of parcel locations
