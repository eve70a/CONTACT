#------------------------------------------------------------------------------------------------------------
# function [ CNTC, ifcver, ierror ] = cntc_initlibrary(outpath, expnam, idebug);
#
# load the library into Python, initialize its internal data and output channels
#
#  in:  character  outpath(*)   - [optional] full path of output directory
#       character  expnam(*)    - [optional] experiment name
#       integer    idebug       - [optional] show (1) or hide (0) error messages
#  out: integer    CNTC         - struct with 'magic numbers' for configuring CONTACT
#       integer    ifcver       - version of the CONTACT add-on
#       integer    ierror       - error flag
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 0: m=*, glob   - no icp needed

from python_intfc          import cntc_dll
from ctypes                import c_int, c_double
from .cntc_getmagicnumbers import cntc_getmagicnumbers

def cntc_initlibrary(outpath=' ', expnam=' ', idebug=1):

    # initialize the internal data of the library, open its output streams

    outpath_bytes = outpath.encode( 'utf-8' )
    expnam_bytes  = expnam.encode( 'utf-8' )

    ifcver        = c_int( )
    ierror        = c_int( )
    ioutput       = c_int( 0 )
    len_outpath   = c_int( len(outpath_bytes) )
    len_expnam    = c_int( len(expnam_bytes) )

    cntc_dll.cntc_initializefirst( ifcver, ierror, ioutput, outpath_bytes, expnam_bytes, 
                                   len_outpath, len_expnam )

    # print('initializeFirst: ifcver = ',ifcver.value,', ierror = ',ierror.value)

    # return a dict with 'magic numbers' for setting flags later on

    CNTC = cntc_getmagicnumbers()

    if (ierror.value==CNTC['err_allow']):
        print('cntc_initialize: no license found or license invalid, check output-file (%d).' % ierror.value);
    elif ierror:
        print('cntc_initlibrary: An error occurred in the CONTACT library (%d)' % ierror.value)
    elif (idebug>=2):
        print('Initialization ok')

    return [CNTC, ifcver.value, ierror.value]

# end function cntc_initlibrary

#------------------------------------------------------------------------------------------------------------
