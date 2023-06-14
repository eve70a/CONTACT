
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_finalizelast
#
# Finalize, clean-up, close output files
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 0: m=*, glob   - no icp needed

from python_intfc          import cntc_dll
from ctypes                import c_int

def cntc_finalizelast():

    cntc_dll.cntc_finalizelast()

# end function cntc_finalizelast

#------------------------------------------------------------------------------------------------------------

