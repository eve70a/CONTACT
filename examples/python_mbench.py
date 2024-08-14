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
            # CONTACTDIR = 'C:\Program Files\Vtech CMCC\contact_v24.1'
        print('Setting CONTACTDIR as', CONTACTDIR)
    sys.path.append( CONTACTDIR )

if (not importlib.util.find_spec('python_intfc')):
    sys.exit('ERROR: cannot find the Python interface to CONTACT.')

import python_intfc as cntc

# $Revision: 2593 $, $Date: 2024-08-14 15:18:54 +0200 (Wed, 14 Aug 2024) $

#------------------------------------------------------------------------------------------------------------
# Part 1: Initialize the CONTACT library, register two problems 
#         "iwhe = 1" for left wheel, "iwhe = 2" for right wheel
#------------------------------------------------------------------------------------------------------------

# outdir = 'c:\\temp'
wrkdir  = ' '
outdir  = ' '
expnam  = ' '
idebug  = 1
[CNTC, ifcver, ierror] = cntc.initlibrary(wrkdir, outdir, expnam, idebug)

if ierror:
   sys.exit('An error occurred, ierror = %d' % ierror)

# set global flags: level of debug output

params = np.array([CNTC['if_idebug'] ], dtype=c_int) 
values = np.array([      1           ], dtype=c_int) 
lenflg = len(params)

cntc.setglobalflags(params, values)

# initialize two result elements using 'module 1'

imodul = 1;  # w/r contact
for iwhe in [1, 2]:
   [ifcver, ierror] = cntc.initialize(iwhe, imodul)

#-----------------------------------------------------------------------  
# Part 2: Configure the main flags & control digits of the contact problem
#-----------------------------------------------------------------------  
 
for iwhe in [1, 2]:

   flags = np.zeros( 10, dtype=c_int )
   values = np.zeros( 10, dtype=c_int )
                                           # CONTACT unit convention: [mm], [mm/s], [N], acting on body (1)
   flags[0] = CNTC['if_units' ]; values[0] = CNTC['un_cntc']; 
   flags[1] = CNTC['ic_config']; values[1] = iwhe-1; # C1 digit (CONFIG): 0 for left side, 1 for right side
   flags[2] = CNTC['ic_tang'  ]; values[2] = 3;      # T=3, steady state rolling
   flags[3] = CNTC['ic_pvtime']; values[3] = 2;      # P=2, no previous time
   flags[4] = CNTC['ic_norm'  ]; values[4] = 1;      # N=1, FZ prescribed
   flags[5] = CNTC['ic_discns']; values[5] = 2;      # D=2, planar contact
   flags[6] = CNTC['if_wrtinp']; values[6] = 0;      # no .inp-file needed
   flags[7] = CNTC['ic_matfil']; values[7] = 0;      # no .mat-file needed
   flags[8] = CNTC['ic_output']; values[8] = 3;      # O=1, min. output to .out-file
   flags[9] = CNTC['ic_flow'  ]; values[9] = 2;      # W=2, little progress output

   cntc.setflags(iwhe, [], flags, values)

   # material data

   gg    = 82000;  # [N/mm^2]
   poiss = 0.28;   # [-]

   cntc.setmaterialproperties(iwhe, [], gg, poiss, gg, poiss);

   # friction data

   if (iwhe <= 2):     # left wheel

      imeth = 0;      # L-digit 0: Coulomb friction
      fstat = 0.30;   # [-]
      params = fstat * np.ones( 2, dtype=c_double )
      cntc.setfrictionmethod(iwhe, [], imeth, params)

   elif (1):          # right wheel, main example

      imeth = 10;     # VL-digit 10: Coulomb friction, variable across rail
      nvf = 2;
      params = np.array( [ nvf, -20*pi/180, 0.20, 0.20, -10*pi/180, 0.30, 0.30 ], dtype=c_double )
      cntc.setfrictionmethod(iwhe, [], imeth, params);

 # else               # right wheel, falling friction

 #    imeth  = 14;    # L-digit 4: exponential falling friction
 #    memdst = 0.003; # [mm]
 #    mem_s0 = 1.000; # [mm/s]
 #    #        alpha [rad], fkin, fexp1, sabsh1, -, -
 #    table  = [ -20*pi/180, 0.09, 0.11, 1250, 0, 0, memdst, mem_s0;
 #               -10*pi/180, 0.14, 0.16, 1250, 0, 0, memdst, mem_s0 ];
 #    nvf    = size(table,1);
 #    params = [ nvf, reshape(table', 1, 8*nvf) ];
 #    cntc_setfrictionmethod(iwhe, [], imeth, length(params), params);
 # end

   # total vertical force specified, see N-digit above

   fz = 10000;   # [N]

   cntc.setverticalforce(iwhe, fz)

   # grid discretization

   dx     = 0.2;   # [mm]
   ds     = 0.2;   # [mm]
   a_sep  = pi/2;  # [rad]
   d_sep  = 8.0;   # [mm]
   d_comb = 4.0;   # [mm]

   params = np.array( [dx, ds, a_sep, d_sep, d_comb], dtype=c_double )

   cntc.setpotcontact(iwhe, [], -1, params)

   dqrel = 1;
   cntc.setrollingstepsize(iwhe, [], [], dqrel);

   # track dimensions & deviations at current side (Z=3)

   ztrack = 3;
   params = np.array( [ 14, 0, 1435, 0, 0, 0, 0, 0, 0, 0, 0 ], dtype=c_double )
   cntc.settrackdimensions(iwhe, ztrack, params)

   # rail profile

   itype = -1;
   sclfac = 1; smooth = 0; rparam = [sclfac, smooth];

   if (iwhe==1): # using Simpack prr-file
      mirror_y = 0; iparam = [itype, 0, mirror_y];
      cntc.setprofileinputfname(iwhe, 'MBench_UIC60_v3.prr', iparam, rparam)
   else:         # using Miniprof ban-file
      mirror_y = 1; iparam = [itype, 0, mirror_y];
      cntc.setprofileinputfname(iwhe, 'MBench_UIC60_v3.ban', iparam, rparam)

   # wheelset dimensions

   ewheel = 3;
   params = np.array( [1360, -70, 460], dtype=c_double )
   cntc.setwheelsetdimensions(iwhe, ewheel, params);

   # wheel profile

   mirror_y = 0; iparam = [itype, 0, mirror_y];
   cntc.setprofileinputfname(iwhe, 'MBench_S1002_v3.prw', iparam, rparam)

   ewheel = 1
   prm = np.zeros(12)
   cntc.setwheelsetflexibility(iwhe, ewheel, prm)

   # positions for subsurface stress calculation

#  iblk = 1
#  isubs = 1
#  nz = 11
#  zl = 1e-6
#  dz = 0.5
#  cntc.subs_addblock(iwhe, [], iblk, isubs, [], [], [nz, zl, dz])

   iblk  = 1
   isubs = 9
   xp = [  0. ]
   yp = [  0. ]
   zp = [ 0, 0.1, 0.2 ]
   cntc.subs_addblock(iwhe, [], iblk, isubs, xp, yp, zp)

# end for iwhe

#-----------------------------------------------------------------------  
# Part 3: perform loops for all the cases to be computed
#-----------------------------------------------------------------------  

# set wheelset positions according to the example

y_ws    = np.arange( 0,   10.5,    0.5 ) # final value 10.5 not included
yaw_ws  = np.arange( 0, 0.0252, 0.0012 )
roll_ws = np.array( [  0.00000000, -0.00002304, -0.00005049, -0.00008103, -0.00011280,
                      -0.00014570, -0.00018030, -0.00021680, -0.00025570, -0.00029770,
                      -0.00035540, -0.00047770, -0.00062720, -0.00437600, -0.00639300,
                      -0.00764200, -0.00860600, -0.00940800, -0.01010113, -0.01071386,
                      -0.01126431 ] )
vpitch  = np.array( [ -4.34811810, -4.34741340, -4.34657520, -4.34624400, -4.34591970,
                      -4.34556270, -4.34515030, -4.34466880, -4.34409300, -4.34337150,
                      -4.33536370, -4.33188640, -4.32937180, -4.27488340, -4.26356290,
                      -4.25757470, -4.25348570, -4.25032450, -4.24775610, -4.24556650,
                      -4.24363750 ] )

# initialize result structures similar to parse_out1 results

results = []
for iwhe in [1, 2]:
   results.append( dict() )
   results[iwhe-1]['ws_pos']   = { 'x':[],  'y':[],  'z':[],  'roll':[],  'yaw':[],  'pitch':[],
                                   'vx':[], 'vy':[], 'vz':[], 'vroll':[], 'vyaw':[], 'vpitch':[] }
   results[iwhe-1]['tot_forc'] = { 'fx_tr':[], 'fy_tr':[], 'fz_tr':[], 'fx_ws':[], 'fy_ws':[], 'fz_ws':[] }
   results[iwhe-1]['npatch']   = list()
   results[iwhe-1]['cp_pos']   = { 'xtr':[], 'ytr':[], 'ztr':[], 'delttr':[], 'yr':[], 'zr':[],
                                   'xw':[], 'yw':[], 'zw':[], 'ncon':[], 'nadh':[], 'nslip':[] }
   results[iwhe-1]['cp_creep'] = { 'pen':[], 'cksi':[], 'ceta':[], 'cphi':[], 'veloc':[] }
   results[iwhe-1]['cp_force'] = { 'fn':[], 'fx':[], 'fs':[], 'mz':[], 'elen':[], 'fric':[],
                                   'sigvm':[], 'vm_x':[], 'vm_y':[], 'vm_z':[] }

# loop over the lateral displacements y_ws

ncase = len(y_ws)
for iwhe in [1, 2]:  # wheel number

   for jcase in range(0, ncase):
      icase = jcase + 1

      ws_pos = np.array( [    0, y_ws[jcase], 0, roll_ws[jcase], yaw_ws[jcase],   0      ], dtype=c_double )
      ws_vel = np.array( [ 2000,     0,       0,      0,            0,     vpitch[jcase] ], dtype=c_double )

      ewheel = 2
      cntc.setwheelsetposition(iwhe, ewheel, ws_pos)
      cntc.setwheelsetvelocity(iwhe, ewheel, ws_vel)

      # compute the contact problem

      print('Starting case %d for wheel %d...' % (icase, iwhe))

      ierror = cntc.calculate(iwhe);
      if (ierror!=0): 
         break

      # compute subsurface stresses

      ierror = cntc.subs_calculate(iwhe)

      # get total forces on upper body (1) (global coordinates)

      values = cntc.getglobalforces(iwhe)
      results[iwhe-1]['tot_forc']['fx_tr'].append( values[1-1] )
      results[iwhe-1]['tot_forc']['fy_tr'].append( values[2-1] )
      results[iwhe-1]['tot_forc']['fz_tr'].append( values[3-1] )
      results[iwhe-1]['tot_forc']['fx_ws'].append( values[7-1] )
      results[iwhe-1]['tot_forc']['fy_ws'].append( values[8-1] )
      results[iwhe-1]['tot_forc']['fz_ws'].append( values[9-1] )

      # get number of contact patches

      npatch = cntc.getnumcontactpatches(iwhe)
      results[iwhe-1]['npatch'].append( npatch )

      # get detailed results per contact patch

      # add list for [icase-1]
      results[iwhe-1]['cp_pos'  ]['xtr'   ].append( [] )
      results[iwhe-1]['cp_pos'  ]['ytr'   ].append( [] )
      results[iwhe-1]['cp_pos'  ]['ztr'   ].append( [] )
      results[iwhe-1]['cp_pos'  ]['delttr'].append( [] )
      results[iwhe-1]['cp_pos'  ]['yr'    ].append( [] )
      results[iwhe-1]['cp_pos'  ]['zr'    ].append( [] )
      results[iwhe-1]['cp_pos'  ]['xw'    ].append( [] )
      results[iwhe-1]['cp_pos'  ]['yw'    ].append( [] )
      results[iwhe-1]['cp_pos'  ]['zw'    ].append( [] )
      results[iwhe-1]['cp_creep']['veloc' ].append( [] )
      results[iwhe-1]['cp_creep']['cksi'  ].append( [] )
      results[iwhe-1]['cp_creep']['ceta'  ].append( [] )
      results[iwhe-1]['cp_creep']['cphi'  ].append( [] )
      results[iwhe-1]['cp_creep']['pen'   ].append( [] )
      results[iwhe-1]['cp_force']['fn'    ].append( [] )
      results[iwhe-1]['cp_force']['fx'    ].append( [] )
      results[iwhe-1]['cp_force']['fs'    ].append( [] )
      results[iwhe-1]['cp_force']['mz'    ].append( [] )
      results[iwhe-1]['cp_force']['sigvm' ].append( [] )
      results[iwhe-1]['cp_force']['vm_x'  ].append( [] )
      results[iwhe-1]['cp_force']['vm_y'  ].append( [] )
      results[iwhe-1]['cp_force']['vm_z'  ].append( [] )

      for icp in range(1, npatch+1):

         # get contact reference location

         values = cntc.getcontactlocation(iwhe, icp)

         # add data for [icp-1]
         results[iwhe-1]['cp_pos']['xtr'   ][icase-1].append( values[1-1] )
         results[iwhe-1]['cp_pos']['ytr'   ][icase-1].append( values[2-1] )
         results[iwhe-1]['cp_pos']['ztr'   ][icase-1].append( values[3-1] )
         results[iwhe-1]['cp_pos']['delttr'][icase-1].append( values[4-1] )
         results[iwhe-1]['cp_pos']['yr'    ][icase-1].append( values[6-1] )
         results[iwhe-1]['cp_pos']['zr'    ][icase-1].append( values[7-1] )
         results[iwhe-1]['cp_pos']['xw'    ][icase-1].append( values[10-1] )
         results[iwhe-1]['cp_pos']['yw'    ][icase-1].append( values[11-1] )
         results[iwhe-1]['cp_pos']['zw'    ][icase-1].append( values[12-1] )

         # get reference velocity

         veloc = cntc.getreferencevelocity(iwhe, icp)
         results[iwhe-1]['cp_creep']['veloc'][icase-1].append( veloc )

         # get penetration and creepages 

         pen = cntc.getpenetration(iwhe, icp)
         cksi, ceta, cphi = cntc.getcreepages(iwhe, icp)

         results[iwhe-1]['cp_creep']['pen' ][icase-1].append( pen )
         results[iwhe-1]['cp_creep']['cksi'][icase-1].append( cksi )
         results[iwhe-1]['cp_creep']['ceta'][icase-1].append( ceta )
         results[iwhe-1]['cp_creep']['cphi'][icase-1].append( cphi )

         # get forces and moment in local coordinates

         fn, tx, ty, mz = cntc.getcontactforces(iwhe, icp)

         results[iwhe-1]['cp_force']['fn'][icase-1].append( fn )
         results[iwhe-1]['cp_force']['fx'][icase-1].append( tx )
         results[iwhe-1]['cp_force']['fs'][icase-1].append( ty )
         results[iwhe-1]['cp_force']['mz'][icase-1].append( mz )

         # get maximum von mises stress

         nx, ny, nz = cntc.subs_getblocksize(iwhe)

         iblk = 1
         table = cntc.subs_getresults(iwhe, icp, iblk, [1,2,3,8,9])
    #    [vm_max, ii_max] = max(table(:,4));

    #    results[iwhe-1]['cp_force']['sigvm'](icase,icp) = vm_max;
    #    results[iwhe-1]['cp_force']['vm_x'](icase,icp)  = table(ii_max,1);
    #    results[iwhe-1]['cp_force']['vm_y'](icase,icp)  = table(ii_max,2);
    #    results[iwhe-1]['cp_force']['vm_z'](icase,icp)  = table(ii_max,3);

      # end for icp

      # get tractions at y_ws = 8mm (1st patch only)

      if (icase==2 and iwhe==1):
         mx, my = cntc.getnumelements(iwhe, 1)
         print('mx=',mx,', my=',my)
         pn, px, py = cntc.gettractions(iwhe, 1)
         print('Array pn: size=',pn.size,', shape=',pn.shape)

   # end for jcase
# end for iwhe

tcpu, twall = cntc.getcalculationtime(1)
print('ire=1: tcpu=',tcpu,', twall=',twall)

print('Done')

cntc.closelibrary()

#------------------------------------------------------------------------------------------------------------
