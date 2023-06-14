#------------------------------------------------------------------------------------------------------------
# function [ ifcver, ierror ] = cntc_initialize(ire, imodul, outpath, idebug)
#
# upon first call: initialize the addon internal data and initialize output channels,
#                  print version information;
# for each ire:   initialize and return the addon version number.
#
#  in:  integer    ire          - result element ID
#       integer    imodul       - module number 1=w/r contact, 3=basic contact
#       character  outpath(*)   - [optional] full path of output directory
#       integer    idebug       - [optional] show (1) or hide (0) error messages
#  out: integer    ifcver       - version of the CONTACT add-on
#       integer    ierror       - error flag
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 0: m=*, glob   - no icp needed

import sys
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, byref, pointer
from .cntc_getmagicnumbers import cntc_getmagicnumbers

def cntc_initialize(ire, imodul, outpath=' ', idebug=1):

    if (not imodul):
        sys.exit('cntc_initialize: please select module 1 or 3')

    # initialize the internal data of the library, open its output streams

    outpath_bytes = outpath.encode( 'utf-8' )

    ifcver        = c_int()
    ierror        = c_int()
    len_outpath   = c_int( len(outpath_bytes) )
    c_imodul      = c_int(imodul)
    p_ire         = pointer(c_int(ire))

    cntc_dll.cntc_initialize( p_ire, c_imodul, ifcver, ierror, outpath_bytes, len_outpath)

    # print('cntc_initialize: ifcver = ',ifcver.value,', ierror = ',ierror.value)

    if (ierror.value):
        CNTC = cntc_getmagicnumbers()
        if (idebug>=1 and ierror.value==CNTC['err_allow']):
            print('cntc_initialize: no license found or license invalid, check output-file (%d).' % ierror.value);
        elif (idebug>=1 and ierror.value<0):
            print('cntc_initialize: an error occurred in the CONTACT library (%d).' % ierror.value);
    # endif (ierror)

    return [ifcver.value, ierror.value]

# end function cntc_initialize

#------------------------------------------------------------------------------------------------------------

