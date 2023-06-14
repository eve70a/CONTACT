
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_finalize(ire)
#
# Finalize calculations and clean-up for a result element
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 0: m=*, glob   - no icp needed

from python_intfc          import cntc_dll
from ctypes                import c_int

def cntc_finalize(ire=1):

    cntc_dll.cntc_finalize(c_int(ire))

# end function cntc_finalize

#------------------------------------------------------------------------------------------------------------

