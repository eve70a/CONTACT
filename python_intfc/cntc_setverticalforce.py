
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_setverticalforce(ire, fz)
#
# set the total vertical force between the contacting bodies for a w/r contact problem (module 1)
# Note: this function sets control digit N = 1
#
#  fz             - total vertical force between the two bodies [force]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 2: m=1, wtd    - no icp needed

from python_intfc          import cntc_dll
from ctypes                import c_int, c_double

def cntc_setverticalforce(ire, fz):

    if (not isinstance(ire, int)):
        ire = 1

    cntc_dll.cntc_setverticalforce(c_int(ire), c_double(fz));

# end function cntc_setverticalforce

#------------------------------------------------------------------------------------------------------------

