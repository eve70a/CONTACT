
#------------------------------------------------------------------------------------------------------------
# function [ rvalues ] = cntc_getwheelsetposition(ire)
#
# return the wheelset position a w/r contact problem (module 1 only)
#
#  rvalues(lenarr)   - wheelset position parameters. 
#
#  The following values are returned, if permitted by the length of rvalues:
#   1 - S_WS      - s-position of the wheelset center of mass along the track center line
#   2 - Y_WS      - lateral y-position of the wheelset center of mass in track coordinates
#   3 - Z_WS      - vertical z-position of the wheelset center of mass in track coordinates
#   4 - ROLL_WS   - wheelset roll angle with respect to the track plane
#   5 - YAW_WS    - wheelset yaw angle with respect to the track center line x_tr
#   6 - PITCH_WS  - wheelset pitch angle, i.e. rotation about the wheelset axle
#
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 2: m=1, wtd    - no icp needed

import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_getwheelsetposition(ire=1):

    lenarr = 6
    values = np.zeros(lenarr, dtype=c_double)

    cntc_dll.cntc_getwheelsetposition(c_int(ire), c_int(lenarr), values.ctypes.data_as(POINTER(c_double)))

    return values

# end function cntc_getwheelsetposition

#------------------------------------------------------------------------------------------------------------

