#------------------------------------------------------------------------------------------------------------
# function [ ierror ] = cntc_readinpfile(ire, inp_type, fname)
#
#  read settings from inp-file
#
#  in:  integer    ire          - result element ID
#       integer    inp_type     - type of inp-file: CNTC_inp_spck, ...
#       character  fname(*)     - filename
#  out: integer    ierror       - error flag
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 2: m=1, wtd    - no icp needed

import sys
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, byref, pointer
from .cntc_getmagicnumbers import cntc_getmagicnumbers

def cntc_readinpfile(ire, inp_type, fname):

    if (not isinstance(ire, int)):
        ire = 1
    if (not isinstance(inp_type, int)):
        sys.exit('cntc_readinpfile: please provide type of input-file')

    fname = fname.encode( 'utf-8' )

    ierror        = c_int()
    len_fname     = c_int( len(fname) )
    c_inptype     = c_int(inp_type)
    p_ire         = pointer(c_int(ire))

    cntc_dll.cntc_readinpfile( p_ire, c_inptype, fname, len_fname, ierror)

    return ierror.value

# end function cntc_readinpfile

#------------------------------------------------------------------------------------------------------------

