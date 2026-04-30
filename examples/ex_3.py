# execution: python ex_3.py
# make sure the pyhips library path is in your PYTHONPATH environment variable

###############################################################################

import pyhips
import numpy as np

############################ set parameters

C_param = 0.5
forceTurb = True
nVar = 1
ScHips = np.array([1.0])
performReaction = False

############################ setup the tree

hips = pyhips.pyhips(C_param, forceTurb, nVar, ScHips, performReaction)  # more parameters have default values

########## initialize pseudo CFD variables

nCFD = 4

var = [ np.zeros(20),            # 20 particles in cell 0
        np.zeros(40),            # 40 particles in cell 0
        np.zeros(60),            # 60 particles in cell 0
        np.zeros(80) ]           # 80 particles in cell 0

for i in range(len(var)): 
    var[i][len(var[i])//2:] = 1.0

w = [ np.full(20, 1.0/20),       # uniform particle weights in cell 0
      np.full(40, 1.0/40),       # uniform particle weights in cell 1
      np.full(60, 1.0/60),       # uniform particle weights in cell 2
      np.full(80, 1.0/80) ]      # uniform particle weights in cell 3

########## USER CFD steps...

########## do HiPS mixing in each CFD cell:

for iCell in range(nCFD):

    print(f"running HiPS in cell {iCell}")

    #---------- reset tree

    tau0 = 1.0
    nLevels = 9
    domainLength = 1.0

    hips.set_tree(nLevels, domainLength, tau0)

    #---------- set variables (CFD particles to HiPS parcels)

    hips.set_varData(var[iCell], w[iCell], "mixf")

    #---------- advance hips

    tRun = 300.0
    hips.calculateSolution(tRun)

    #---------- get variables (HiPS parcels to CFD particles)

    var[iCell] = hips.get_varData()[0]

########## USER CFD steps...

