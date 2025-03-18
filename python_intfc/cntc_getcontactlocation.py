
#------------------------------------------------------------------------------------------------------------
# function [ rvalues ] = cntc_getcontactlocation(ire, icp)
#
# return the contact reference location for a contact problem (module 1 only)
#
#  rvalues(lenarr)   - contact location parameters. 
#
#  The following values are returned, if permitted by the length of rvalues:
#   1 - XCP_TR    - x-position of the contact reference point in track coordinates
#   2 - YCP_TR    - y-position of the contact reference point in track coordinates
#   3 - ZCP_TR    - z-position of the contact reference point in track coordinates
#   4 - DELTCP_TR - contact reference angle: rotation about track x-axis from the track positive z-axis to
#                   the contact positive n-axis, with sign according the right-hand rule
#
#   5 - XCP_R     - x-position of the contact reference point in rail profile coordinates
#   6 - YCP_R     - y-position of the contact reference point in rail profile coordinates
#   7 - ZCP_R     - z-position of the contact reference point in rail profile coordinates
#   8 - SCP_R     - s-parameter of the contact reference point, measured along the rail profile
#   9 - DELTCP_R  - rotation about rail x-axis from rail positive z-axis to contact positive n-axis
#
#  10 - XCP_W     - x-position of the contact reference point in wheel profile coordinates
#  11 - YCP_W     - y-position of the contact reference point in wheel profile coordinates
#  12 - ZCP_W     - z-position of the contact reference point in wheel profile coordinates
#  13 - SCP_W     - s-parameter of the contact reference point, measured along the wheel profile
#  14 - DELTCP_W  - rotation about wheel x-axis from wheel positive z-axis to contact positive n-axis
#
#  15 - XPN_TR    - x-position of the pressure center of gravity in track coordinates
#  16 - YPN_TR    - y-position of the pressure center of gravity in track coordinates
#
#  18 - DY_DEFL   - lateral rail shift according to massless rail deflection
#  19 - DZ_DEFL   - vertical rail shift according to massless rail deflection
#
#  21 - XW_TR     - x-position of wheel profile marker in track coordinates
#  22 - YW_TR     - y-position of wheel profile marker in track coordinates
#  23 - ZR_TR     - z-position of wheel profile marker in track coordinates
#  24 - ROLLW_TR  - roll angle of wheel profile marker in track coordinates
#  25 - YAWW_TR   - yaw angle of wheel profile marker in track coordinates
#
#  26 - YR_TR     - y-position of rail profile marker in track coordinates
#  27 - ZR_TR     - z-position of rail profile marker in track coordinates
#  28 - ROLLR_TR  - roll angle of rail profile marker in track coordinates
#
#  30 - XCP_WS    - x-position of the contact reference point in wheelset coordinates
#  31 - YCP_WS    - y-position of the contact reference point in wheelset coordinates
#  32 - ZCP_WS    - z-position of the contact reference point in wheelset coordinates
#  33 - DELTCP_WS - rotation about wheelset x-axis from wheelset positive z-axis to contact positive n-axis
#
#  34 - XW_WS     - x-position of wheel profile marker in wheelset coordinates
#  35 - YW_WS     - y-position of wheel profile marker in wheelset coordinates
#  36 - ZR_WS     - z-position of wheel profile marker in wheelset coordinates
#  37 - ROLLW_WS  - roll angle of wheel profile marker in wheelset coordinates
#  38 - YAWW_WS   - yaw angle of wheel profile marker in wheelset coordinates
#
#  The "contact reference point" is the origin of the contact local coordinate system. It is determined by
#  a heuristic rule and is centered within the contact patch in a weighted sense.
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 3: m=1, cp     - require icp>0, default 1

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_getcontactlocation(ire=1, icp=1):

    if (icp<=0):
        sys.exit('ERROR in cntc_getcontactlocation: not available for icp=%d' % icp)

    lenarr = 32
    values = np.zeros(lenarr, dtype=c_double)

    cntc_dll.cntc_getcontactlocation(c_int(ire), c_int(icp), c_int(lenarr),
                                                                 values.ctypes.data_as(POINTER(c_double)))

    return values

# end function cntc_getcontactlocation

#------------------------------------------------------------------------------------------------------------

