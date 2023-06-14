
#------------------------------------------------------------------------------------------------------------
# function [ sens ] = cntc_getsensitivities(ire, icp, lenout, lenin)
#
# return the sensitivities of the total forces for a contact problem
#
#  lenout, lenin      - requested number of outputs (forces) and inputs (creepages or shifts)
#  sens(lenout,lenin) - matrix of sensitivities
#
#  the inputs are ordered   1: pen, 2: cksi, 3: ceta, 4: cphi
#  in rolling, T=2,3, the units are pen [length], cksi, ceta [-],      cphi [angle/length]
#  in shifts,  T=1,   the units are pen [length], cksi, ceta [length], cphi [angle]
#
#  the outputs are ordered  1: fn,  2: fx,   3: fy,   4: mz
#  the units are fn [force], fx, fy [-], mz [force.length]
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 6: m=*, cp     - require icp>0, default 1

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER

def cntc_getsensitivities(ire, icp, lenout, lenin):
    print('function called for first time:', __name__)

    if (not ire):
        ire = 1
    if (not icp):
        icp = 1
    if (icp<=0):
        sys.exit('ERROR in cntc_getsensitivities: not available for icp=%d' % icp)
    if (lenout*lenin<=0):
        sys.exit('ERROR in cntc_getsensitivities: invalid lenout=%d, lenin=%d' % (lenout, lenin))

    lenarr = lenout * lenin
    sens = np.zeros(lenarr, dtype=c_double)

    cntc_dll.cntc_getsensitivities(c_int(ire), c_int(icp), c_int(lenout), c_int(lenin), 
                                                                 sens.ctypes.data_as(POINTER(c_double)))

    # lenout runs fastest in the library -> rows in numpy -> transpose to (lenout,lenin)-array
    sens = np.transpose( np.reshape(sens, (lenin, lenout)) )

    return sens

# end function cntc_getsensitivities

#------------------------------------------------------------------------------------------------------------

