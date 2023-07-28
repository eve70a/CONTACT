
#------------------------------------------------------------------------------------------------------------
# function [ ] = cntc_setmaterialproperties(ire, icp, g1, nu1, g2, nu2)
#
# set the material properties for a contact problem.    Obsolete, replaced by cntc_setmaterialparameters.
#
#  g1, g2         - modulus of rigidity for body 1, 2 [force/area]
#  nu1, nu2       - Poisson's ratio for body 1, 2 [-]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 5: m=*, wtd    - default icp=-1

from python_intfc          import cntc_dll
from ctypes                import c_int, c_double

def cntc_setmaterialproperties(ire, icp, g1, nu1, g2, nu2):
    # default: W/R contact, all patches

    if (not isinstance(ire, int)):
        ire =  1
    if (not isinstance(icp, int)):
        icp = -1
    cntc_dll.cntc_setmaterialproperties(c_int(ire), c_int(icp), c_double(g1), c_double(nu1), 
                                                                c_double(g2), c_double(nu2))

# end function cntc_setmaterialproperties

#------------------------------------------------------------------------------------------------------------

