
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_setmaterialparameters(ire, icp, m_digit, rparam)
#
# set the material parameters for a contact problem
#
# M-digit = 0: rparam = [ poiss1, poiss2, g1,  g2                                 ];
#           1: rparam = [ poiss1, poiss2, g1,  g2,  fg1,   fg2,    vt1,    vt2    ];
#           2: rparam = [ poiss1, poiss2, g1,  g2,  flx,   k0_mf,  alfamf, betamf ];
#           3: rparam = [ poiss1, poiss2, g1,  g2,  k0_mf, alfamf, betamf         ];
#           4: rparam = [ poiss1, poiss2, g1,  g2,  g3,    laythk, tau_c0, k_tau  ];
# values >= 10 are used to configure the M2-digit, M2 = m_digit - 10
# M2-digit = 2: rparam = [ cdampn, cdampt, dfnmax, dftmax ];
#
# dimensions: poiss1, poiss2  [-],        g1, g2    [force/area],   
#             fg1, fg2        [-],        vt1, vt2  [time],
#             flx         [volume/force], k0_mf, alfamf, betamf [-], 
#             g3, tau_c0  [force/area],   laythk    [length],        k_tau [force/volume]
#             cdampn, cdampt  [time],     dfnmax, dftmax   [force/time]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 5: m=*, wtd    - default icp=-1

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_setmaterialparameters(ire, icp, m_digit, rparam):
    # default: W/R contact, all patches

    if (not isinstance(ire, int)):
        ire =  1
    if (not isinstance(icp, int)):
        icp = -1
    if (not isinstance(m_digit, int)):
        m_digit = 0 # default: linearly elastic material
    # convert rparam to NumPy ndarray
    if (not isinstance(rparam, np.ndarray)):
        if (isinstance(rparam, list)):
            rparam = np.array( rparam, dtype=c_double )
        else:
            sys.exit('ERROR in cntc_setmaterialparameters: invalid rparam provided.')

    nparam = len(rparam)

    cntc_dll.cntc_setmaterialparameters(c_int(ire), c_int(icp), c_int(m_digit), c_int(nparam),
                                                                 rparam.ctypes.data_as(POINTER(c_double)))

# end function cntc_setmaterialparameters

#------------------------------------------------------------------------------------------------------------

