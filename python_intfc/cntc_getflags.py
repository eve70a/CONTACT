
#------------------------------------------------------------------------------------------------------------
# function [ values ] = cntc_getflags(ire, icp, params)
#
# used for retrieving various configuring flags from a contact problem
#
#  lenflg         - length of params/values arrays
#  params(lenflg) - codes of the parameters to be obtained from CONTACT
#  values(lenflg) - values of the parameters obtained from CONTACT
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 7: m=*, wtd or cp - default icp=-1

import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_getflags(ire, icp, params):

    # default: W/R contact, overall data
    if (not isinstance(ire, int)):
        ire =  1
    if (not isinstance(icp, int)):
        icp = -1

    # convert params to NumPy ndarray, allowing for scalar or list of scalars
    if (not isinstance(params, np.ndarray)):
        params = np.atleast_1d( np.array( params, dtype=c_int ) )

    lenarr = len(params)
    values = np.zeros(lenarr, dtype=c_double)

    cntc_dll.cntc_getflags(c_int(ire), c_int(icp), c_int(lenarr), params.ctypes.data_as(POINTER(c_double)),
                                                                  values.ctypes.data_as(POINTER(c_double)))

    return values

# end function cntc_getflags

#------------------------------------------------------------------------------------------------------------

