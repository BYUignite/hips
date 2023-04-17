from __future__ import division
import glob
import numpy as np
from data_tools import get_nRlz, get_dataHeaderVars, get_inputFileParameter
from data_tools import get_data_realization, extrap1d, get_domainBounds
from data_tools import compute_pdf
from scipy.interpolate import UnivariateSpline
#--------------------------------------------------------------------------------------------

def hips_pdf(DI, LSR=False) :

   
    #--------------------------------------------------------------------------------------------

    dataFiles = glob.glob(DI['cdir']+"data_py/data_*.npy")
    ntimes    = len(dataFiles)
    varNames  = get_dataHeaderVars(DI)
    times     = get_inputFileParameter(DI, ("dumpTimes",))    
    nvar      = len(varNames)
    nLevels   = get_inputFileParameter(DI, ("params", "nLevels"))    
    nx = int(2**(nLevels-1))
    nrlz      = get_nRlz(DI, nx)
   

    #--------------------------------------------------------------------------------------------
    chi = np.zeros(0)


    for itime in range(ntimes):

        for irlz in range(nrlz):

            data          = get_data_realization(DI, itime, irlz, nx)
            chiColumn     = data[:,1]
            chi           = np.append(chi, chiColumn)
    #---------------------------------------------------------------------------------------------

    log_chi     = np.log10([ c for c in chi if c > 0 ])
    log_chi_min = -12
    log_chi_max = np.max(log_chi) 
    Plogchi,chibins     = compute_pdf(log_chi, log_chi_min , log_chi_max, nbins =100)
    #--------------------------------------------------------------------------------------------

    head = "chibins        "

    for i in range(1, len(varNames), 2):

        head     += "P10(log)_"+ varNames[i] + "   "

    fname    = DI['pdir'] + "PDF_Chi.dat"

    var       = np.column_stack((chibins,Plogchi))

    np.savetxt(fname, var, header=head, fmt="%15.8e ")

    
    




  
