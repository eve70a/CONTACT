
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_setmetadata(ire, icp, params, values)
#
# used for configuring various metadata for a contact problem
#
#  lenmta         - length of params/values arrays
#  params(lenmta) - codes of the metadata to be communicated to CONTACT
#  values(lenmta) - values of the metadata to be communicated to CONTACT
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 5: m=*, wtd    - default icp=-1

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_setmetadata(ire, icp, params, values):

    if (not ire):
        ire = 1
    if (not icp):
        icp = -1  # default: W/R contact, all patches

    # convert params and values to NumPy ndarrays
    if (not isinstance(params, np.ndarray)):
        if (isinstance(params, list)):
            params = np.array( params, dtype=c_int )
        else:
            sys.exit('ERROR in cntc_setmetadata: invalid params provided.')
    if (not isinstance(values, np.ndarray)):
        if (isinstance(values, list)):
            values = np.array( values, dtype=c_double )
        else:
            sys.exit('ERROR in cntc_setmetadata: invalid values provided.')

    lenmta = len(params)

    if (len(values)!=lenmta):
        sys.exit('ERROR in cntc_setmetadata: invalid params/values provided.')

    cntc_dll.cntc_setmetadata(c_int(ire), c_int(icp), c_int(lenmta), 
                             params.ctypes.data_as(POINTER(c_int)), values.ctypes.data_as(POINTER(c_double)))

# end function cntc_setmetadata

#------------------------------------------------------------------------------------------------------------

