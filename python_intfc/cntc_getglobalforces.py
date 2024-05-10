
#------------------------------------------------------------------------------------------------------------
# function [ values ] = cntc_getglobalforces(ire, icp)
#
# return the overall forces for a w/r contact problem (module 1 only)
#
#  rvalues(lenarr)   - contact forces [force] and moments [force.length]
#
#  The following values are returned, if permitted by the length of rvalues:
#   1 - FX_TR    - total force on the output body, component in track longitudinal x-direction
#   2 - FY_TR    - total force on the output body, component in track lateral y-direction
#   3 - FZ_TR    - total force on the output body, component in track vertical z-direction
#   4 - MX_RR_TR - total moment on output body about rail profile marker, component in track x-direction
#   5 - MY_RR_TR - total moment on output body about rail profile marker, component in track y-direction
#   6 - MZ_RR_TR - total moment on output body about rail profile marker, component in track z-direction
#
#   7 - FX_WS    - total force on the output body, component in wheelset longitudinal x-direction
#   8 - FY_WS    - total force on the output body, component in wheelset lateral y-direction
#   9 - FZ_WS    - total force on the output body, component in wheelset vertical z-direction
#  10 - MX_RW_WS - total moment on output body about wheel profile marker, component in wheelset x-direction
#  11 - MY_RW_WS - total moment on output body about wheel profile marker, component in wheelset y-direction
#  12 - MZ_RW_WS - total moment on output body about wheel profile marker, component in wheelset z-direction
#                  Note that the 'output body' is the rail when using CONTACTs unit convention
#
#  13 - FX_RR    - total force on the output body, component in rail profile x-direction
#  14 - FY_RR    - total force on the output body, component in rail profile y-direction
#  15 - FZ_RR    - total force on the output body, component in rail profile z-direction
#  16 - MX_RR_TR - total moment on output body about rail profile marker, component in rail x-direction
#  17 - MY_RR_TR - total moment on output body about rail profile marker, component in rail y-direction
#  18 - MZ_RR_TR - total moment on output body about rail profile marker, component in rail z-direction
#
#  19 - FX_RW    - total force on the output body, component in wheel profile x-direction
#  20 - FY_RW    - total force on the output body, component in wheel profile y-direction
#  21 - FZ_RW    - total force on the output body, component in wheel profile z-direction
#  22 - MX_RW_WS - total moment on output body about wheel profile marker, component in wheel x-direction
#  23 - MY_RW_WS - total moment on output body about wheel profile marker, component in wheel y-direction
#  24 - MZ_RW_WS - total moment on output body about wheel profile marker, component in wheel z-direction
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 4: m=1, wtd/cp - default icp=-1

import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_getglobalforces(ire=1, icp=-1):

    lenarr = 30
    values = np.zeros(lenarr, dtype=c_double)

    cntc_dll.cntc_getglobalforces(c_int(ire), c_int(icp), c_int(lenarr), 
                                                                 values.ctypes.data_as(POINTER(c_double)))

    return values

# end function cntc_getglobalforces

#------------------------------------------------------------------------------------------------------------

