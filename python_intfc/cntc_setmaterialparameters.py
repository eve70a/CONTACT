
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_setmaterialparameters(ire, icp, m_digit, rparam)
#
# set the M-digit and material parameters for a contact problem
#
# values < 10 are used to configure the M-digit
#    0: purely elastic material,           params = [ nu1, nu2, g1, g2 ]
#    1: visco-elastic material,            params = [ nu1, nu2, g1, g2, fg1,   fg2,    vt1,    vt2    ]
#    2: modified Fastsim, 1 flexibility    params = [ nu1, nu2, g1, g2, flx,   k0_mf,  alfamf, betamf ]
#    3: modified Fastsim, 3 flexibilities  params = [ nu1, nu2, g1, g2, k0_mf, alfamf, betamf         ]
#    4: elastic + elasto-plastic 3rd body  params = [ nu1, nu2, g1, g2, g3,    laythk, tau_c0, k_tau  ]
#    5: modified FaStrip                   params = [ nu1, nu2, g1, g2, k0_mf, alfamf, betamf ]
# dimensions: nu1, nu2        [-],        g1, g2    [force/area],   
#             fg1, fg2        [-],        vt1, vt2  [time],
#             flx         [volume/force], k0_mf, alfamf, betamf [-], 
#             g3, tau_c0   [force/area],  laythk    [length],        k_tau [force/volume]
#
# values in [20,29] are used to configure the M2-digit, M2 = imeth - 20
#   20: no force-based proportional damping   params = [ ]
#   22: force-based proportional damping      params = [ cdampn, cdampt, dfnmax, dftmax ]
#             cdampn, cdampt  [time],     dfnmax, dftmax   [force/time]
#
# values in [30,39] are used to configure the M3-digit, M3 = imeth - 30
#   30: no added flexibility                   params = [ ]
#   32: add flexibility at small displacement  params = [ uthrs, lpow0, lpow1 ]
#             uthrs [length],             lpow0, lpow1 [-]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2026 by Vtech CMCC.
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

