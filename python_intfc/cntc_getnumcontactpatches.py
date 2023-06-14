
#------------------------------------------------------------------------------------------------------------
# function [ npatch ] = cntc_getnumcontactpatches(ire)
#
# return the number of contact patches used in a w/r contact problem
#
#  npatch        - number of separate contact patches
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 2: m=1, wtd    - no icp needed

from python_intfc          import cntc_dll
from ctypes                import c_int

def cntc_getnumcontactpatches(ire=1):

    npatch = c_int()

    cntc_dll.cntc_getnumcontactpatches(c_int(ire), npatch)

    return npatch.value

# end function cntc_getnumcontactpatches

#------------------------------------------------------------------------------------------------------------

