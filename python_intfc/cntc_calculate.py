
#------------------------------------------------------------------------------------------------------------
# function [ ierror ] = cntc_calculate(ire, icp, idebug)
#
# perform actual CONTACT calculation for a contact problem
# in:   integer    idebug       - show warnings (2), errors (1) or hide all messages (0)
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 5: m=*, wtd    - default icp=-1

from python_intfc          import cntc_dll
from ctypes                import c_int
from .cntc_getmagicnumbers import cntc_getmagicnumbers

def cntc_calculate(ire=1, icp=-1, idebug=1):

    ierr = c_int(-999)

    try:
        cntc_dll.cntc_calculate(c_int(ire), c_int(icp), ierr)
    except Exception as e:
        print("An exception occurred..."+ str(e))

    if (ierr.value):
        CNTC = cntc_getmagicnumbers()
        if (idebug>=1 and ierr.value==CNTC['err_allow']):
            print('cntc_calculate: no valid license found for CONTACT library (%d).' % ierr.value)
        elif (idebug>=1 and ierr.value==CNTC['err_profil']):
            print('cntc_calculate: an error is found in the rail or wheel profile specification (%d).'
                                                                                     % ierr.value)
        elif (idebug>=1 and ierr.value==CNTC['err_broydn']):
            print('cntc_calculate: no solution found in Broyden algorithm (%d).'     % ierr.value)
        elif (idebug>=1 and ierr.value<0):
            print('cntc_calculate: an error occurred in the CONTACT library (%d).'   % ierr.value)
        elif (idebug>=2 and ierr.value>0):
            print('cntc_calculate: potential contact may be too small: there are %d' % ierr.value, 
                                                                      'points adjacent to the boundaries.' )
    # endif (ierr)

    return ierr.value

# end function cntc_calculate

#------------------------------------------------------------------------------------------------------------

