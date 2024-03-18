function [] = matlab_varprof(expnam)

% function [] = matlab_varprof(expnam)
%
% Test program for variable profiles: switches and crossings, wheel out-of-roundness
%    expnam = cross_brute, cross+wing, mbench_intrup, etc.

   EXP_CROSS_BRUTE         =  1; EXP_CROSS_LOCUS         =  2; EXP_CROSS_WING          =  3;
   EXP_CW_INTERRUPT        =  4; EXP_MBENCH_BRUTE        =  5; EXP_MBENCH_INTRUP       =  6;
   EXP_MBENCH_LOCUS        =  7; EXP_TWO_PATCHES         =  8; EXP_WING_BRUTE          =  9;
   EXP_WING_LOCUS          = 10; EXP_CHALMERS_FLAT_FZ125 = 11; EXP_ROUNDED_FLAT_D09    = 12;

   known_exp = strvcat('cross_brute',        'cross_locus',            'cross+wing',     ...
                       'cw_interrupt',       'mbench_brute',           'mbench_intrup',  ...
                       'mbench_locus',       'two_patches',            'wing_brute',     ...
                       'wing_locus',         'chalmers_flat_fz125',    'rounded_flat_d09');
   iexp = 0;
   for i = 1 : size(known_exp,1)
      if (strcmp(expnam, deblank(known_exp(i,:))))
         iexp = i;
      end
   end
   if (iexp==0)
      disp([' Unknown experiment ',expnam,', aborting.']);
      return
   end

   %---------------------------------------------------------------------------------------------------------
   % Initialization of the CONTACT library
   %---------------------------------------------------------------------------------------------------------

   [CNTC, ifcver, ierror] = cntc_initlibrary;

   idebug = 1;             % 1: just a bit of information from the library
   cntc_setglobalflags(CNTC.if_idebug, idebug);

   imodul =  1;
   iwhe   =  1;
   icp    = -1;
   [ifcver, ierror] = cntc_initialize(iwhe, imodul);

   % Configure result element, setting permanent data

   % set flags: configure control digits & output, defaults, may be changed per experiment

   clear flags values;
   flags( 1) = CNTC.if_units ; values( 1) = CNTC.un_cntc;
   flags( 2) = CNTC.ic_config; values( 2) = 1;    % C1: 0=left wheel, 1=right
   flags( 3) = CNTC.ic_tang  ; values( 3) = 3;    % T=3: steady state rolling
   flags( 4) = CNTC.ic_pvtime; values( 4) = 2;    % P=2: no previous time
   flags( 5) = CNTC.ic_discns; values( 5) = 2;    % D=2: planar contact
   flags( 6) = CNTC.if_wrtinp; values( 6) = 0;    %   0: no .inp-file needed
   flags( 7) = CNTC.ic_matfil; values( 7) = 0;    % A=0: no .mat-file needed
   flags( 8) = CNTC.ic_output; values( 8) = 3;    % O=1: min. output to .out-file
   flags( 9) = CNTC.ic_flow  ; values( 9) = 2;    % W=2: little progress output
   flags(10) = CNTC.ic_npomax; values(10) = 20000; % max #elements in pot.contact

   cntc_setflags(iwhe, [], flags, values);

   % G=0: set maxgs, maxin, maxnr, maxout, eps

   gdigit  = 0;
   values  = [ 999, 100, 30, 1 ];
   rvalues = [ 1d-5 ];
   cntc_setsolverflags(iwhe, [], gdigit, values, rvalues);

   % material data using M = 0 (fully elastic model)

   gg    = 82000;  % [N/mm^2]
   poiss = 0.28;   % [-]

   mdigit = 0;
   rparam = [ poiss, poiss, gg, gg ];
   cntc_setmaterialparameters(iwhe, [], mdigit, rparam);

   % friction data using L = 0: Coulomb friction

   ldigit = 0;      % L-digit 0: Coulomb friction
   fstat  = 0.30;   % [-]
   cntc_setfrictionmethod(iwhe, [], ldigit, [fstat,fstat]);

   % set rolling step size: ratio c = dq / dx

   dqrel = 1;
   cntc_setrollingstepsize(iwhe, [], [], dqrel);

   % set track dimensions & deviations

   ztrack = 3;
   params = [ -1, 750, 0, 0, 0, 0, 0, 0, 0, 0, 0 ];
   cntc_settrackdimensions_new(iwhe, ztrack, params);

   %---------------------------------------------------------------------------------------------------------
   % Further settings dependent on experiment name
   %---------------------------------------------------------------------------------------------------------

   if (iexp==EXP_CROSS_BRUTE | iexp==EXP_CROSS_LOCUS)

      clear flags values;
      if (iexp==EXP_CROSS_LOCUS)
         flags(1) = CNTC.ic_discns ; values(1) = 2;            % D=2: contact locus method (default)
      else
         flags(1) = CNTC.ic_discns ; values(1) = 5;            % D=5: brute force method
      end
      cntc_setflags(iwhe, [], flags, values);

      % set grid discretization

      ipotcn = -1;     % -1 = w/r contact with fixed grid sizes & combination of patches
      dx     = 0.2;    % [mm]
      ds     = 0.05;   % [mm]
      a_sep  = pi/2;   % [rad]
      d_sep  =  8.0;   % [mm]
      d_comb =  4.0;   % [mm]
      d_turn = 12.0;   % [mm]
      params = [dx, ds, a_sep, d_sep, d_comb, d_turn];

      cntc_setpotcontact(iwhe, [], ipotcn, params);

      % set variable rail profile for switch/crossing

      fname   = '../profiles/cross_nose.slcs';

      itype   = -1;     % using filename extension
      mirrory =  0;     % no mirroring
      mirrorz = -1;     % no mirroring
      sclfac  =  1;     % already in [mm], no scaling
      smooth  =  0;     % no smoothing
      iparam = [ itype, 0, mirrory, mirrorz ];
      rparam = [ sclfac, smooth ];

      cntc_setprofileinputfname(iwhe, fname, iparam, rparam);

      ncase = 1;
      s_ws(1:ncase)     =     0.0; 
      y_ws(1:ncase)     =     0.0; 
      z_ws(1:ncase)     =    -2.1275; 
      fz_ws(1:ncase)    =    -1.0; 
      pitch_ws(1:ncase) =     0.0; 
      yaw_ws(1:ncase)   =     0.0; 
      roll_ws(1:ncase)  =     0.0; 
      vpitch(1:ncase)   =    -4.34811810; 

   elseif (iexp==EXP_CROSS_WING)

      clear flags values;
      flags(1) = CNTC.ic_discns ; values(1) = 2;            % D=2: contact locus method (default)
      cntc_setflags(iwhe, [], flags, values);

      % set grid discretization

      ipotcn = -1;     % -1 = w/r contact with fixed grid sizes & combination of patches
      dx     =  0.2;   % [mm]
      ds     =  0.05;  % [mm]
      a_sep  = pi/2;   % [rad]
      d_sep  =  8.0;   % [mm]
      d_comb =  4.0;   % [mm]
      d_turn = 12.0;   % [mm]
      params = [dx, ds, a_sep, d_sep, d_comb, d_turn];

      cntc_setpotcontact(iwhe, [], ipotcn, params);

      % set variable rail profile for switch/crossing

      fname   = '../profiles/cross+wing.slcs';

      itype   = -1;     % using filename extension
      mirrory =  0;     % no mirroring
      mirrorz = -1;     % no mirroring
      sclfac  =  1;     % already in [mm], no scaling
      smooth  =  0;     % no smoothing
      iparam = [ itype, 0, mirrory, mirrorz ];
      rparam = [ sclfac, smooth ];

      cntc_setprofileinputfname(iwhe, fname, iparam, rparam);

      ncase = 1;
      s_ws(1:ncase)     =   100.0; 
      y_ws(1:ncase)     =     0.0; 
      z_ws(1:ncase)     =     2.5212; 
      fz_ws(1:ncase)    =    -1.0; 
      pitch_ws(1:ncase) =     0.0; 
      yaw_ws(1:ncase)   =     0.0; 
      roll_ws(1:ncase)  =     0.0; 
      vpitch(1:ncase)   =    -4.34811810; 

   elseif (iexp==EXP_CW_INTERRUPT)

      clear flags values;
      flags(1) = CNTC.ic_discns ; values(1) = 2;            % D=2: contact locus method (default)
      cntc_setflags(iwhe, [], flags, values);

      % set grid discretization

      ipotcn = -1;     % -1 = w/r contact with fixed grid sizes & combination of patches
      dx     =  0.2;   % [mm]
      ds     =  0.05;  % [mm]
      a_sep  = pi/2;   % [rad]
      d_sep  =  8.0;   % [mm]
      d_comb =  4.0;   % [mm]
      d_turn = 12.0;   % [mm]
      params = [dx, ds, a_sep, d_sep, d_comb, d_turn];

      cntc_setpotcontact(iwhe, [], ipotcn, params);

      % set variable rail profile for switch/crossing

      fname   = '../profiles/cross+wing_extd.slcs';

      itype   = -1;     % using filename extension
      mirrory =  0;     % no mirroring
      mirrorz = -1;     % no mirroring
      sclfac  =  1;     % already in [mm], no scaling
      smooth  =  0;     % no smoothing
      iparam = [ itype, 0, mirrory, mirrorz ];
      rparam = [ sclfac, smooth ];

      cntc_setprofileinputfname(iwhe, fname, iparam, rparam);

      ncase = 2;
      s_ws(1:ncase)     =       0.0; 
      y_ws(1:ncase)     = [    40.0,           0.0    ];
      z_ws(1:ncase)     = [    -1.6736,        2.5287 ];
      fz_ws(1:ncase)    = [    10.0e3,        -1.0    ];
      pitch_ws(1:ncase) =       0.0; 
      yaw_ws(1:ncase)   =       0.0; 
      roll_ws(1:ncase)  =       0.0; 
      vpitch(1:ncase)   =      -4.34811810; 

   elseif (iexp==EXP_MBENCH_BRUTE |  iexp==EXP_MBENCH_LOCUS)

      clear flags values;
      if (iexp==EXP_MBENCH_LOCUS)
         flags(1) = CNTC.ic_discns ; values(1) = 2;            % D=2: contact locus method (default)
      else
         flags(1) = CNTC.ic_discns ; values(1) = 5;            % D=5: brute force method
      end
      cntc_setflags(iwhe, [], flags, values);

      % set grid discretization

      ipotcn = -1;     % -1 = w/r contact with fixed grid sizes & combination of patches
      dx     =  0.2;   % [mm]
      ds     =  0.2;   % [mm]
      a_sep  = pi/2;   % [rad]
      d_sep  =  8.0;   % [mm]
      d_comb =  4.0;   % [mm]
      d_turn = 12.0;   % [mm]
      params = [dx, ds, a_sep, d_sep, d_comb, d_turn];

      cntc_setpotcontact(iwhe, [], ipotcn, params);

      % set variable rail profile for switch/crossing

      fname   = '../profiles/uk_crossing.slcs';

      itype   = -1;     % using filename extension
      mirrory =  0;     % no mirroring
      mirrorz = -1;     % no mirroring
      sclfac  =  1e3;   % data in [m], 1000x scaling
      smooth  =  0;     % no smoothing
      iparam = [ itype, 0, mirrory, mirrorz ];
      rparam = [ sclfac, smooth ];

      cntc_setprofileinputfname(iwhe, fname, iparam, rparam);

      ncase = 35;
      s_ws(1:ncase)     =     239.0; 
      y_ws(1:ncase)     = [     0.80 + 0.01 * [1:ncase]  ];
      z_ws(1:ncase)     =       0.0; 
      fz_ws(1:ncase)    =      10.0e3;
      pitch_ws(1:ncase) =       0.0; 
      yaw_ws(1:ncase)   =       0.0; 
      roll_ws(1:ncase)  =       0.0; 
      vpitch(1:ncase)   =      -4.34811810; 

   elseif (iexp==EXP_MBENCH_INTRUP)

      clear flags values;
      flags(1) = CNTC.ic_discns ; values(1) = 2;            % D=2: contact locus method (default)
      cntc_setflags(iwhe, [], flags, values);

      % set grid discretization

      ipotcn = -1;     % -1 = w/r contact with fixed grid sizes & combination of patches
      dx     =  0.2;   % [mm]
      ds     =  0.2;   % [mm]
      a_sep  = pi/2;   % [rad]
      d_sep  =  8.0;   % [mm]
      d_comb =  4.0;   % [mm]
      d_turn = 12.0;   % [mm]
      params = [dx, ds, a_sep, d_sep, d_comb, d_turn];

      cntc_setpotcontact(iwhe, [], ipotcn, params);

      % set variable rail profile for switch/crossing

      fname   = '../profiles/uk_interrupt_v2.slcs';

      itype   = -1;     % using filename extension
      mirrory =  0;     % no mirroring
      mirrorz = -1;     % no mirroring
      sclfac  =  1e3;   % data in [m], 1000x scaling
      smooth  =  0;     % no smoothing
      iparam = [ itype, 0, mirrory, mirrorz ];
      rparam = [ sclfac, smooth ];

      cntc_setprofileinputfname(iwhe, fname, iparam, rparam);

      ncase  = 43;
      s_ws(1:ncase)     = [  110.0,   111.0,   112.0,   113.0,   114.0,   115.0,   115.5,   116.0, ...
                             116.5,   117.0,   117.5,   118.0,   118.1,   118.2,   118.3,   118.4, ...
                             118.5,   118.6,   118.7,   118.8,   118.9,   119.0,   119.1,   119.2, ...
                             119.3,   119.4,   119.5,   119.6,   119.7,   119.8,   119.9,   120.0, ...
                             120.5,   121.0,   122.0,   123.0,   124.0,   125.0,   126.0,   127.0, ...
                             128.0,   129.0,   130.0 ];
      y_ws(1:ncase)     =     0.78; 
      z_ws(1:ncase)     =    -0.87; 
      fz_ws(1:ncase)    =    -1.0; 
      pitch_ws(1:ncase) =     0.0; 
      yaw_ws(1:ncase)   =     0.0; 
      roll_ws(1:ncase)  =     0.0; 
      vpitch(1:ncase)   =    -4.34811810; 

   elseif (iexp==EXP_TWO_PATCHES)

      clear flags values;
      flags(1) = CNTC.ic_discns ; values(1) = 5;         % D=5: brute force method
      cntc_setflags(iwhe, [], flags, values);

      % set grid discretization

      ipotcn = -1;     % -1 = w/r contact with fixed grid sizes & combination of patches
      dx     =  0.2;   % [mm]
      ds     =  0.2;   % [mm]
      a_sep  = pi/2;   % [rad]
      d_sep  =  8.0;   % [mm]
      d_comb =  4.0;   % [mm]
      d_turn = 12.0;   % [mm]
      params = [dx, ds, a_sep, d_sep, d_comb, d_turn];

      cntc_setpotcontact(iwhe, [], ipotcn, params);

      % set variable rail profile for switch/crossing

      fname   = '../profiles/uk_interrupt_v2.slcs';

      itype   = -1;     % using filename extension
      mirrory =  0;     % no mirroring
      mirrorz = -1;     % no mirroring
      sclfac  =  1e3;   % data in [m], 1000x scaling
      smooth  =  0;     % no smoothing
      iparam = [ itype, 0, mirrory, mirrorz ];
      rparam = [ sclfac, smooth ];

      cntc_setprofileinputfname(iwhe, fname, iparam, rparam);

      ncase = 1;
      s_ws(1:ncase)     =    160.0;  
      y_ws(1:ncase)     =      3.52; 
      z_ws(1:ncase)     =      0.0; 
      fz_ws(1:ncase)    =     10.0e3;
      pitch_ws(1:ncase) =      0.0; 
      yaw_ws(1:ncase)   =      0.0; 
      roll_ws(1:ncase)  =      0.0; 
      vpitch(1:ncase)   =     -4.34811810; 

   elseif (iexp==EXP_WING_BRUTE | iexp==EXP_WING_LOCUS)

      clear flags values;
      if (iexp==EXP_WING_LOCUS)
         flags(1) = CNTC.ic_discns ; values(1) = 2;            % D=2: contact locus method (default)
      else
         flags(1) = CNTC.ic_discns ; values(1) = 5;            % D=5: brute force method
      end
      cntc_setflags(iwhe, [], flags, values);

      % set grid discretization

      ipotcn = -1;     % -1 = w/r contact with fixed grid sizes & combination of patches
      dx     =  0.2;   % [mm]
      ds     =  0.05;  % [mm]
      a_sep  = pi/2;   % [rad]
      d_sep  =  8.0;   % [mm]
      d_comb =  4.0;   % [mm]
      d_turn = 12.0;   % [mm]
      params = [dx, ds, a_sep, d_sep, d_comb, d_turn];

      cntc_setpotcontact(iwhe, [], ipotcn, params);

      % set variable rail profile for switch/crossing

      fname   = '../profiles/wing_rail.slcs';

      itype   = -1;     % using filename extension
      mirrory =  0;     % no mirroring
      mirrorz = -1;     % no mirroring
      sclfac  =  1;     % already in [mm], no scaling
      smooth  =  0;     % no smoothing
      iparam = [ itype, 0, mirrory, mirrorz ];
      rparam = [ sclfac, smooth ];

      cntc_setprofileinputfname(iwhe, fname, iparam, rparam);

      ncase = 1
      s_ws(1:ncase)     =  100.0; 
      y_ws(1:ncase)     =    0.0; 
      z_ws(1:ncase)     =    0.6589; 
      fz_ws(1:ncase)    =   -1.0; 
      pitch_ws(1:ncase) =    0.0; 
      yaw_ws(1:ncase)   =    0.0; 
      roll_ws(1:ncase)  =    0.0; 
      vpitch(1:ncase)   =   -4.34811810; 

   elseif (iexp==EXP_CHALMERS_FLAT_FZ125)

      clear flags values;
      flags(1) = CNTC.ic_discns ; values(1) = 2;         % D=2: contact locus method (default)
      flags(2) = CNTC.ic_npomax ; values(2) = 20000;     % max #elements in pot.contact
      cntc_setflags(iwhe, [], flags, values);

      % set track dimensions & deviations

      ztrack = 3;
      params = [ 14,   0,   1435,   0.020,   0,   0,   0,   0,   0,   0,   0 ];
      cntc_settrackdimensions_new(iwhe, ztrack, params);

      % set grid discretization

      ipotcn = -1;     % -1 = w/r contact with fixed grid sizes & combination of patches
      dx     =  0.4;   % [mm]
      ds     =  0.4;   % [mm]
      a_sep  = pi/2;   % [rad]
      d_sep  =  8.0;   % [mm]
      d_comb =  4.0;   % [mm]
      d_turn = 12.0;   % [mm]
      params = [dx, ds, a_sep, d_sep, d_comb, d_turn];

      cntc_setpotcontact(iwhe, [], ipotcn, params);

      % set constant rail profile

      fname   = '../../examples/r300_wide.prr';

      itype   = -1;     % using filename extension
      mirrory =  0;     % no mirroring
      mirrorz = -1;     % no mirroring
      sclfac  =  1;     % data in [mm], no scaling
      smooth  =  0;     % no smoothing
      iparam = [ itype, 0, mirrory, mirrorz ];
      rparam = [ sclfac, smooth ];

      cntc_setprofileinputfname(iwhe, fname, iparam, rparam);

      ncase = 26;
      s_ws(1:ncase)     =       0.0; 
      y_ws(1:ncase)     =       0.0; 
      z_ws(1:ncase)     =       0.0; 
      fz_ws(1:ncase)    =     125.0e3;
      pitch_ws(1:ncase) = [   -24.0 - 1.0 * [1:ncase] ] * pi/180;
      yaw_ws(1:ncase)   =       0.0; 
      roll_ws(1:ncase)  =       0.0; 
      vpitch(1:ncase)   =      -4.08190679; 

   elseif (iexp==EXP_ROUNDED_FLAT_D09)

      clear flags values;
      flags(1) = CNTC.ic_discns ; values(1) = 2;         % D=2: contact locus method (default)
      flags(2) = CNTC.ic_npomax ; values(2) = 20000;     % max #elements in pot.contact
      cntc_setflags(iwhe, [], flags, values);

      % set track dimensions & deviations

      ztrack = 3;
      params = [ -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ];
      cntc_settrackdimensions_new(iwhe, ztrack, params);

      % set grid discretization

      ipotcn = -1;     % -1 = w/r contact with fixed grid sizes & combination of patches
      dx     =  0.2;   % [mm]
      ds     =  0.2;   % [mm]
      a_sep  = pi/2;   % [rad]
      d_sep  =  8.0;   % [mm]
      d_comb =  4.0;   % [mm]
      d_turn = 12.0;   % [mm]
      params = [dx, ds, a_sep, d_sep, d_comb, d_turn];

      cntc_setpotcontact(iwhe, [], ipotcn, params);

      % set constant rail profile

      fname   = '../profiles/circ_r300.prr';

      itype   = -1;     % using filename extension
      mirrory =  0;     % no mirroring
      mirrorz =  0;     % auto mirroring
      sclfac  =  1;     % data in [mm], no scaling
      smooth  =  0;     % no smoothing
      iparam = [ itype, 0, mirrory, mirrorz ];
      rparam = [ sclfac, smooth ];

      cntc_setprofileinputfname(iwhe, fname, iparam, rparam);

      ncase = 15;
      s_ws(1:ncase)     =       0.0; 
      y_ws(1:ncase)     =       0.0; 
      z_ws(1:ncase)     =       0.0; 
      fz_ws(1:ncase)    =      10.0e3;
      pitch_ws(1:ncase) = [    -8.0 + 1.0 * [1:ncase] ] * pi/180;
      yaw_ws(1:ncase)   =       0.0; 
      roll_ws(1:ncase)  =       0.0; 
      vpitch(1:ncase)   =      -4.44444444; 

   else

      disp(['ERROR: unknown experiment ',expnam,', aborting.'])
      return

   end

   % set wheelset dimensions

   if     (iexp==EXP_CHALMERS_FLAT_FZ125)
      rvalues(1:3) = [ 1360,   -70,   490 ];
   elseif (iexp==EXP_ROUNDED_FLAT_D09)
      rvalues(1:3) = [    0,     0,   450 ];
   else
      rvalues(1:3) = [ 1360,   -70,   460 ];
   end

   ewheel = 3;
   cntc_setwheelsetdimensions(iwhe, ewheel, rvalues);

   % set wheel profile

   if     (iexp==EXP_CHALMERS_FLAT_FZ125)
      fname   = '../../examples/S1002_flat.slcw';
      smooth    = 5;      % lambda smoothing
   elseif (iexp==EXP_ROUNDED_FLAT_D09)
      fname   = '../profiles/flat_d09.slcw';
      smooth    = 0;      % no smoothing
   else
      fname   = '../profiles/MBench_S1002_v3.prw';
      smooth    = 0;      % no smoothing
   end

   itype   = -1;     % using filename extension
   mirrory =  0;     % no y-mirroring
   mirrorz = -1;     % no z-mirroring
   sclfac  =  1;     % already in [mm], no scaling
   iparam = [ itype, 0, mirrory, mirrorz ];
   rparam = [ sclfac, smooth ];

   cntc_setprofileinputfname(iwhe, fname, iparam, rparam);

   %---------------------------------------------------------------------------------------------------------
   % Run cases, with outputs in out-file
   %---------------------------------------------------------------------------------------------------------

   for icase = 1 : ncase

      % adjust discretisation settings for cw_interrupt

      if (icase==2 & iexp==EXP_CW_INTERRUPT)
         clear flags values;
         flags(1) = CNTC.ic_discns ; values(1) = 5;          % D=5: brute force method
         cntc_setflags(iwhe, [], flags, values);
      end

      % set wheelset position & velocity

      rvalues = [ s_ws(icase), y_ws(icase), z_ws(icase), roll_ws(icase), yaw_ws(icase), pitch_ws(icase) ];
      ewheel  = 2;
      cntc_setwheelsetposition(iwhe, ewheel, rvalues);

      rvalues = [ 2000,   0,   0,   0,   0,   vpitch(icase) ];
      ewheel  = 2;
      cntc_setwheelsetvelocity(iwhe, ewheel, rvalues);

      if (fz_ws(icase)>0)
         cntc_setverticalforce(iwhe, fz_ws(icase));

         clear flags values;
         flags(1) = CNTC.ic_norm; values(1) = 1;            % N=1: vert.force prescribed
         cntc_setflags(iwhe, [], flags, values);
      else
         clear flags values;
         flags(1) = CNTC.ic_norm; values(1) = 0;            % N=0: approach prescribed
         cntc_setflags(iwhe, [], flags, values);
      end

      if (idebug>=1)
         disp(sprintf('  ...start processing case %d', icase));
      end

      % perform the actual calculation for this wheel

      ierr = cntc_calculate(iwhe, icp);

      % check/report error conditions

      if (ierr==CNTC.err_allow)
         disp('      ERROR: no valid license found')
         return;
      elseif (ierr==CNTC.err_profil)
         disp('      ERROR: the rail and/or wheel profile files could not be found or processed')
         return
      elseif (ierr<0)
         disp(sprintf('      ERROR: an error occured in the CONTACT library, ierr=%d.', ierr));
      end

      itask  = 1;  % 1=rail or wheel profile coordinates, 2=track coords
      itype  = 1;  % 0=rail, 1=wheel
      isampl = 0;  % 0=internal sampling, 1=constant ds_out
      prw = cntc_getprofilevalues(iwhe, itask, [itype,isampl], []);

      figure(2); hold on;
      plot(prw(:,1), prw(:,2), '-');

   end % icase

   % finalize, close CONTACT library

   cntc_finalize(iwhe)

end % function matlab_varprof

