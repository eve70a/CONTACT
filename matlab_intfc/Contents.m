%
% Matlab interface for the CONTACT library (trunk).
%
% cntc_initlibrary              - load library into Matlab
% cntc_initialize               - initialize data-structures
% cntc_getmagicnumbers          - get 'magic numbers' for configuring CONTACT
% cntc_setglobalflags           - configuration of global parameters
%
% cntc_setflags                 - configuration of parameters per contact problem
% cntc_setmetadata              - configuration of various metadata in CONTACT
% cntc_setsolverflags           - configuration of solver parameters
%
% cntc_setprofileinputfname     - set a wheel or rail profile filename for a w/r problem
% cntc_setprofileinputvalues    - set a wheel or rail profile (table) for a w/r problem
% cntc_settrackdimensions       - set the track geometry for a w/r problem   New version
% *cntc_settrackdimensions_old  - set the track geometry for a w/r problem   Obsolete: use _new version
% cntc_settrackdimensions_new   - set the track geometry for a w/r problem
% cntc_setwheelsetdimensions    - set wheelset geometry specification
% cntc_setwheelsetposition      - set wheelset position specification
% cntc_setwheelsetvelocity      - set wheelset velocity specification
% cntc_setwheelsetflexibility   - set wheelset flexibility specification
% cntc_setverticalforce         - set total vertical force
%
% cntc_setmaterialparameters    - set material parameters
% cntc_settemperaturedata       - set heat related material parameters
% *cntc_setmaterialproperties   - set material parameters            Obsolete: use setmaterialparameters
% *cntc_setinterfaciallayer     - set third body layer parameters    Obsolete: use setmaterialparameters
% cntc_setfrictionmethod        - set friction parameters
% cntc_settimestep              - set time step size used in shift problems
% cntc_setreferencevelocity     - set rolling velocity
% cntc_setrollingstepsize       - set rolling direction and step size
%
% cntc_sethertzcontact          - set Hertzian problem specification
% cntc_setpotcontact            - set grid discretization, pot.contact area for non-Hertzian cases
%
% cntc_setpenetration           - set approach/penetration
% cntc_setnormalforce           - set total normal force
% cntc_setundeformeddistc       - set undeformed distance elementwise
%
% cntc_setcreepages             - set creepages
% cntc_setextrarigidslip        - set extra term of the rigid slip elementwise
% cntc_settangentialforces      - set total tangential forces
%
% cntc_calculate                - perform actual calculation
%
% cntc_getflags                 - retrieve configuration parameters for a contact problem
% cntc_getparameters            - retrieve configuration data for a contact problem
%
% cntc_getwheelsetposition      - get wheelset position parameters
% cntc_getwheelsetvelocity      - get wheelset velocity parameters
% *cntc_getprofilevalues        - retrieve wheel or rail profile after smoothing  Obsolete: use _new version
% cntc_getprofilevalues_new     - retrieve wheel or rail profile after smoothing
% cntc_getnumcontactpatches     - get the number of separate contact patches
% cntc_getcontactlocation       - get the location of one contact patch
% cntc_getglobalforces          - get total forces in track/wheelset coordinates
% cntc_getreferencevelocity     - get rolling velocity
% cntc_gethertzcontact          - get parameters from a Hertzian contact problem
% cntc_getpotcontact            - get parameters of the potential contact area
%
% cntc_getcpresults             - get results of contact patch in form used by plot3d
% cntc_getcontactforces         - get total forces per patch in local coordinates
% cntc_getpenetration           - get approach/penetration
% cntc_getcreepages             - get creepages
% cntc_getcontactpatchareas     - get size of contact area
% cntc_getmaximumpressure       - get maximum pressure
% cntc_getmaximumtraction       - get maximum shear stress
% cntc_getmaximumtemperature    - get maximum surface temperatures
% cntc_getsensitivities         - get sensitivities of forces w.r.t. creepages
%
% cntc_getelementdivision       - get elementwise adhesion/slip areas
% cntc_getfielddata             - get elementwise output values
% cntc_gettractions             - get elementwise tractions
% cntc_getdisplacements         - get elementwise displacement differences
% cntc_getmicroslip             - get elementwise micro-slip velocity
%
% cntc_getcalculationtime       - get calculation time used
% cntc_resetcalculationtime     - reset timers
%
% subs_addblock                 - define 3D grid for subsurface stress calculation
% subs_calculate                - perform subsurface stress calculation
% subs_getblocksize             - get number of points for one block
% subs_getresults               - get detailed outputs of subsurface stress calculation
%
% cntc_finalize                 - cleanup for one result element
% cntc_closelibrary             - finalize, cleanup and unload the library

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

