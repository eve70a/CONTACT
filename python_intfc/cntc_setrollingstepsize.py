
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_setrollingstepsize(ire, icp, chi, dq)
#
# set the rolling direction and step size for a contact problem
#
# in w/r contact,    icp = -1
#     chi            - ignored
#     dqrel          - rolling step size relative to grid size dx [-]
# in generic contact, icp > 0,
#     chi            - rolling direction [angle]
#     dq             - rolling step size [length]
#
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 5: m=*, wtd    - default icp=-1

from python_intfc          import cntc_dll
from ctypes                import c_int, c_double

def cntc_setrollingstepsize(ire, icp, chi, dq):

    if (not isinstance(ire, int)):
        ire =  1
    if (not isinstance(icp, int)):
        icp = -1; # default: W/R contact, all patches
        chi =  0;

    cntc_dll.cntc_setrollingstepsize(c_int(ire), c_int(icp), c_double(chi), c_double(dq))

# end function cntc_setrollingstepsize

#------------------------------------------------------------------------------------------------------------

