
#------------------------------------------------------------------------------------------------------------
# function [ values ] = cntc_getparameters(ire, icp, itask)
#
# used for retrieving various parameters from a contact problem
#
#  itask       - selected group of parameters (1: cntc_getcpresults, 2: material, 3: friction)
#  values      - values of the parameters obtained from CONTACT
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 6: m=*, cp     - require icp>0, default 1

import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_getparameters(ire, icp, itask):

    # default: W/R contact, all patches
    if (not isinstance(ire, int)):
        ire =  1
    if (not isinstance(icp, int)):
        icp =  1
    if (not isinstance(itask, int)):
        itask = 1

    lenarr = 30;
    tmparr = np.zeros(lenarr, dtype=c_double)

    cntc_dll.cntc_getparameters(c_int(ire), c_int(icp), c_int(itask), c_int(lenarr), 
                                                                  tmparr.ctypes.data_as(POINTER(c_double)))

    values = {}
    if (itask==1):
        values['veloc' ] = tmparr[ 1-1]
        values['chi'   ] = tmparr[ 2-1]
        values['dq'    ] = tmparr[ 3-1]
        values['spinxo'] = tmparr[ 4-1]
        values['spinyo'] = tmparr[ 5-1]
        values['tau_c0'] = tmparr[ 6-1]
    elif (itask==2):
        values['gg1'   ] = tmparr[ 1-1]
        values['gg2'   ] = tmparr[ 2-1]
        values['gg'    ] = tmparr[ 3-1]
        values['poiss1'] = tmparr[ 4-1]
        values['poiss2'] = tmparr[ 5-1]
        values['poiss' ] = tmparr[ 6-1]
        values['ak'    ] = tmparr[ 7-1]
        values['flx1'  ] = tmparr[ 8-1]
        values['flx2'  ] = tmparr[ 9-1]
        values['flx3'  ] = tmparr[10-1]
        values['k0_mf' ] = tmparr[11-1]
        values['alfamf'] = tmparr[12-1]
        values['betamf'] = tmparr[13-1]
        values['k_eff' ] = tmparr[14-1]
        values['gg3'   ] = tmparr[15-1]
        values['laythk'] = tmparr[16-1]
        values['tau_c0'] = tmparr[17-1]
        values['k_tau' ] = tmparr[18-1]
        values['cdampn'] = tmparr[19-1]
        values['cdampt'] = tmparr[20-1]
        values['dfnmax'] = tmparr[21-1]
        values['dftmax'] = tmparr[22-1]
    elif (itask==3):
        values['lmeth' ] = round(tmparr[ 1-1])
        values['nvf'   ] = round(tmparr[ 2-1])
        values['memdst'] = tmparr[ 3-1]
        values['mem_s0'] = tmparr[ 4-1]
        values['alpha' ] = [ tmparr[ 5+3*i-1] for i in range(values['nvf']) ]
        values['fstat' ] = [ tmparr[ 6+3*i-1] for i in range(values['nvf']) ]
        values['fkin'  ] = [ tmparr[ 7+3*i-1] for i in range(values['nvf']) ]

    return values

# end function cntc_getparameters

#------------------------------------------------------------------------------------------------------------
