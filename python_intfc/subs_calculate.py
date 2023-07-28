
#------------------------------------------------------------------------------------------------------------
# function [ ierror ] = subs_calculate(ire, icp, idebug)
#
# perform subsurface stress calculation for a contact problem
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 7: m=*, wtd or cp - default icp=-1

from python_intfc          import cntc_dll
from ctypes                import c_int, c_double
from .cntc_getmagicnumbers import cntc_getmagicnumbers

def subs_calculate(ire, icp, idebug=1):

    if (not isinstance(ire, int)):
        ire =  1
    if (not isinstance(icp, int)):
        icp = -1

    ierr = c_int()

    cntc_dll.subs_calculate(c_int(ire), c_int(icp), ierr)

    if (ierr.value):
        CNTC = cntc_getmagicnumbers()
        if (idebug>=1 and ierr.value==CNTC['err_allow']):
            print('subs_calculate: no valid license found for CONTACT library (%d).' % ierr.value)
        elif (idebug>=1 and ierr.value<0):
            print('subs_calculate: an error occurred in the CONTACT library (%d).' % ierr.value)
    # endif (ierr)

    return ierr.value

# end function subs_calculate

#------------------------------------------------------------------------------------------------------------

