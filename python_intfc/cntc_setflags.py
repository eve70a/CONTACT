
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_setflags(ire, icp, params, values)
#
# used for configuring various flags for a contact problem
#
#  lenflg         - length of params/values arrays
#  params(lenflg) - codes of the parameters to be communicated to CONTACT
#  values(lenflg) - values of the parameters to be communicated to CONTACT
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 5: m=*, wtd    - default icp=-1

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, POINTER

def cntc_setflags(ire, icp, params, values):

    if (not ire):
        ire = 1
    if (not icp):
        icp = -1; # default: W/R contact, all patches

    # convert params and values to NumPy ndarrays, allowing for scalar or list of scalars
    if (not isinstance(params, np.ndarray)):
        params = np.atleast_1d( np.array( params, dtype=c_int ) )
    if (not isinstance(values, np.ndarray)):
        values = np.atleast_1d( np.array( values, dtype=c_int ) )

    lenflg = len(params)
    if (len(values)!=lenflg):
        sys.exit('ERROR in cntc_setflags: invalid params/values provided.');

    cntc_dll.cntc_setflags(c_int(ire), c_int(icp), c_int(lenflg), params.ctypes.data_as(POINTER(c_int)),
                                               values.ctypes.data_as(POINTER(c_int)))

# end function cntc_setflags

#------------------------------------------------------------------------------------------------------------

