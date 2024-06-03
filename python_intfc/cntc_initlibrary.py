#------------------------------------------------------------------------------------------------------------
# function [ CNTC, ifcver, ierror ] = cntc_initlibrary(wrkdir, outdir, expnam, idebug);
#
# load the library into Python, initialize its internal data and output channels
#
#  in:  character  wrkdir(*)    - [optional] effective working folder
#       character  outdir(*)    - [optional] output folder
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

def cntc_initlibrary(wrkdir=' ', outdir=' ', expnam=' ', idebug=1):

    # initialize the internal data of the library, open its output streams

    wrkdir_bytes  = wrkdir.encode( 'utf-8' )
    outdir_bytes  = outdir.encode( 'utf-8' )
    expnam_bytes  = expnam.encode( 'utf-8' )

    ifcver        = c_int( )
    ierror        = c_int( )
    ioutput       = c_int( 0 )
    len_wrkdir    = c_int( len(wrkdir_bytes) )
    len_outdir    = c_int( len(outdir_bytes) )
    len_expnam    = c_int( len(expnam_bytes) )

    cntc_dll.cntc_initializefirst_new( ifcver, ierror, ioutput, wrkdir_bytes, outdir_bytes, expnam_bytes, 
                                                                    len_wrkdir, len_outdir, len_expnam )

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
