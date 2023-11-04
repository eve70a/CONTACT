import os
import sys
import importlib
from   platform  import system
from   math      import pi
import numpy         as np
from   ctypes    import c_int, c_double

#------------------------------------------------------------------------------------------------------------
# Preparation: Locate the CONTACT library, load into Python environment
#------------------------------------------------------------------------------------------------------------

# The CONTACT dll is loaded automatically upon importing the python_intfc package
# print('Cntc_dll=', cntc.cntc_dll)

# Check if the Python interface to CONTACT can be found on the Python search path;
# If not: check environment variable CONTACTDIR, or set a default value

if (not importlib.util.find_spec('python_intfc')):
    # Check if CONTACT installation folder is set, if so, take it from there
    if (os.environ['CONTACTDIR']):
        CONTACTDIR = os.environ['CONTACTDIR']
        print('Found environment variable CONTACTDIR as', CONTACTDIR)
    else:
        # Default/fall-back value for CONTACT installation folder
        if (system() == 'Linux'):
            CONTACTDIR = '/v3/CMCC/contact'
        else:
            CONTACTDIR = 'C:\\CMCC\\contact'
            # CONTACTDIR = 'C:\Program Files\Vtech CMCC\contact_v23.2'
        print('Setting CONTACTDIR as', CONTACTDIR)
    sys.path.append( CONTACTDIR )

if (not importlib.util.find_spec('python_intfc')):
    sys.exit('ERROR: cannot find the Python interface to CONTACT.')

import python_intfc as cntc

# $Revision: 2447 $, $Date: 2023-11-04 15:03:17 +0100 (Sat, 04 Nov 2023) $

#------------------------------------------------------------------------------------------------------------
# Part 1: Initialize the CONTACT library, register problem "iwhe = 1"
#------------------------------------------------------------------------------------------------------------

# outpath = 'c:\\temp'
outpath = ' '
expnam  = ' '
idebug  = 1
[CNTC, ifcver, ierror] = cntc.initlibrary(outpath, expnam, idebug)

if ierror:
   sys.exit('An error occurred, ierror = %d' % ierror)

# set global flags: level of debug output

params = np.array([CNTC['if_idebug'] ], dtype=c_int) 
values = np.array([      1           ], dtype=c_int) 

cntc.setglobalflags(params, values)

# initialize one result element using 'module 3' = basic contact

iwhe   = 1   # "wheel number"
icp    = 1   # "contact problem on wheel"
imodul = 3   # basic contact
[ifcver, ierror] = cntc.initialize(iwhe, imodul)

#------------------------------------------------------------------------------------------------------------
# Part 2: Configure the main flags & control digits of the contact problem
#------------------------------------------------------------------------------------------------------------
 
flags = np.zeros( 10, dtype=c_int )
values = np.zeros( 10, dtype=c_int )
                                              # CONTACT unit convention: [mm], [mm/s], [N], acting on body (1)
flags[0] = CNTC['if_units' ]; values[0] = CNTC['un_cntc']; 
flags[1] = CNTC['if_wrtinp']; values[1] = 0;  # no .inp-file needed
flags[2] = CNTC['ic_matfil']; values[2] = 0;  # no .mat-file needed
flags[3] = CNTC['ic_output']; values[3] = 3;  # O=1, min. output to .out-file
flags[4] = CNTC['ic_flow'  ]; values[4] = 4;  # W=1, a little progress output
flags[5] = CNTC['ic_norm'  ]; values[5] = 1;  # N=1, FN prescribed
flags[6] = CNTC['ic_tang'  ]; values[6] = 3;  # T=3, steady state rolling
flags[7] = CNTC['ic_pvtime']; values[7] = 2;  # P=2, no previous time

cntc.setflags(iwhe, icp, flags, values)

# material data

gg    = 82000   # [N/mm^2]
poiss = 0.28    # [-]

cntc.setmaterialproperties(iwhe, icp, gg, poiss, gg, poiss) 

# Hertzian discretization, automatic potential contact area, ellipticity a/b = 0.5

ipotcn = -2     # ellipticity
mx     = 48 
my     = 44 
a1     = 0.0008 # [1/mm]
aob    = 0.5    # [-]
scale  = 1.1 

params = [mx, my, a1, aob, scale]

cntc.sethertzcontact(iwhe, icp, ipotcn, params)

# rolling velocity

veloc  = 10000  # [mm/s]

cntc.setreferencevelocity(iwhe, icp, veloc)

# total normal force specified, see N-digit above

fn     = 106700 # [N]

cntc.setnormalforce(iwhe, icp, fn)

# configuration of iteration accuracy

g_digit =   0
maxgs   =  80
maxin   = 200
maxnr   =   5
maxout  =   1
eps     =  1e-6
iparam  = [maxgs, maxin, maxnr, maxout]
rparam  = [eps]

cntc.setsolverflags(iwhe, icp, g_digit, iparam, rparam)

#------------------------------------------------------------------------------------------------------------
# Part 3: perform loops for all the cases to be computed
#------------------------------------------------------------------------------------------------------------

# creepages: a list of values for creating a creep-force curve

cksi_m0  = [ 0.00001, 0.0004, 0.0008, 0.0012, 0.0016, 0.0020, 0.0024, 0.0028,
              0.0032, 0.0036, 0.0040, 0.0044, 0.0048, 0.0052, 0.0056, 0.0060,
              0.0080, 0.0100, 0.0200, 0.0300, 0.0400, 0.0500, 0.0750, 0.1000,
              0.1250, 0.1500, 0.1750, 0.2000, 0.2250, 0.2500 ]  # [-]
cksi_m3  = [ 0.00001, 0.0010, 0.0020, 0.0030, 0.0040, 0.0050, 0.0060, 0.0070,
              0.0080, 0.0100, 0.0120, 0.0140, 0.0160, 0.0180, 0.0200, 0.0220,
              0.0250, 0.0300, 0.0350, 0.0400, 0.0500, 0.0650, 0.0800, 0.1000,
              0.1250, 0.1500, 0.1750, 0.2000, 0.2250, 0.2500 ]  # [-]
cksi_m4  = [ 0.00001, 0.0010, 0.0020, 0.0030, 0.0040, 0.0050, 0.0060, 0.0070,
              0.0080, 0.0090, 0.0100, 0.0110, 0.0120, 0.0130, 0.0140, 0.0150,
              0.0160, 0.0180, 0.0200, 0.0300, 0.0400, 0.0500, 0.0750, 0.1000,
              0.1250, 0.1500, 0.1750, 0.2000, 0.2250, 0.2500 ]  # [-]
ceta   = 0  # [-]
cphi   = 0  # [rad/mm]

#  icurv = 1:   Original Kalker theory,     Coulomb friction
#          2:   Extended CONTACT,           exponential falling friction
#          3:   Original Fastsim algorithm, Coulomb friction
#          4:   Modified Fastsim algorithm, exponential falling friction

for icurv in range(1, 5):

    # friction data

    if (icurv==1 or icurv==3):
        imeth = 0       # L-digit 0: Coulomb friction
        fstat = 0.33    # [-]
        fkin  = fstat 
        cntc.setfrictionmethod(iwhe, icp, imeth, [fstat,fkin]);
    elif (icurv==2 or icurv==4):
        imeth = 4       # L-digit 4: exponential falling friction
        if (icurv==2):
            fkin   = 0.14   # [-]
            fexp1  = 0.19
            sabsh1 = 1250   # [mm/s]
        else:
            fkin   = 0.1368 # [-]
            fexp1  = 0.2232
            sabsh1 =  990.0 # [mm/s]
        memdst = 0.003     # [mm]
        mem_s0 = 1.000     # [mm/s]
        params = [fkin, fexp1, sabsh1, 0, 0, memdst, mem_s0]
        cntc.setfrictionmethod(iwhe, icp, imeth, params)

    # interfacial layer

    if (icurv == 1):
        imeth  = 0      # M-digit 0: linear elastic half-spaces, no interfacial layer
        nparam = 0 
        params = [] 
        cntc.setinterfaciallayer(iwhe, icp, imeth, params) 
    elif (icurv == 2):
        imeth  = 4      # M-digit 4: elasto-plastic interfacial layer
        gg3    = 8200   # [N/mm2]
        h3     = 1.25   # [mm]
        params = [gg3, h3, 0, 0] 
        nparam = len(params) 
        cntc.setinterfaciallayer(iwhe, icp, imeth, params) 
    else:
        imeth  = 3      # M-digit 3: modified Fastsim
        if (icurv == 3):
            k0_mf  = 1.0    # [-]
            alfamf = 1.0    # [-]
            betamf = 1.0    # [-]
        else:
            k0_mf  = 0.54   # [-]
            alfamf = 0.02   # [-]
            betamf = 0.62   # [-]
        params = [k0_mf, alfamf, betamf] 
        nparam = len(params) 
        cntc.setinterfaciallayer(iwhe, icp, imeth, params) 

    # Traction bound: for CONTACT or Fastsim

    if (icurv <= 2):
        flags = [ CNTC['ic_bound'] ]; values = [ 0 ]    # B=0, elastic half-space
    else:
        flags = [ CNTC['ic_bound'] ]; values = [ 3 ]    # B=3, parabolic traction bound
    cntc.setflags(iwhe, icp, flags, values)

    if (icurv==1 or icurv==3):
        cksi  = cksi_m0
    elif (icurv==2):
        cksi  = cksi_m4
    else:
        cksi  = cksi_m3

    # loop over all creepages cksi provided

    fx_list = np.ndarray( (len(cksi), 4), dtype=c_double )

    for iksi in range(0, len(cksi)):

        # set creepages according to next value from cksi

        cntc.setcreepages(iwhe, icp, cksi[iksi], ceta, cphi)

        # compute the contact problem

        ierror = cntc.calculate(iwhe, icp)
        if (ierror!=0): 
            break

        # get forces on upper body (1) (CONTACT unit convention)

        fn, fx, fy, mz = cntc.getcontactforces(iwhe, icp)

        fx_list[iksi][icurv-1] = fx

    # end for iksi
# end for icurv

# Cleanup

cntc.finalize(iwhe)
cntc.closelibrary

print('Done')

#------------------------------------------------------------------------------------------------------------
