
#------------------------------------------------------------------------------------------------------------
# function [ values ] = cntc_getparameters(ire, icp)
#
# used for retrieving various parameters from a contact problem needed for cntc_getcpresults
#
#  values      - values of the parameters obtained from CONTACT
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 5: m=*, wtd    - default icp=-1

import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_getparameters(ire, icp):

    # default: W/R contact, all patches
    if (not isinstance(ire, int)):
        ire =  1
    if (not isinstance(icp, int)):
        icp = -1

    lenarr = 10;
    tmparr = np.zeros(lenarr, dtype=c_double)

    cntc_dll.cntc_getparameters(c_int(ire), c_int(icp), c_int(lenarr), 
                                                                  tmparr.ctypes.data_as(POINTER(c_double)))

    values = {}
    values['veloc']  = tmparr[1-1]
    values['chi']    = tmparr[2-1]
    values['dq']     = tmparr[3-1]
    values['spinxo'] = tmparr[4-1]
    values['spinyo'] = tmparr[5-1]
    values['tau_c0'] = tmparr[6-1]

    return values

# end function cntc_getparameters

#------------------------------------------------------------------------------------------------------------
