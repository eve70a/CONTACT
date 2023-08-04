
%------------------------------------------------------------------------------------------------------------
% function [ CNTC ] = cntc_getmagicnumbers();
%
% get the struct with 'magic numbers' for configuring CONTACT
%
%  CNTC         - struct with 'magic numbers' for configuring CONTACT
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 0: m=*, glob   - no icp needed

function [ CNTC ] = cntc_getmagicnumbers();

   % return a struct with 'magic numbers' for setting flags later on

   CNTC.if_units  = 1933; % code for setting the units convention
                          % note that w/r profiles are given in [mm] in all cases
   CNTC.un_cntc   = 1934; % code for selecting CONTACT's unit convention,
                          %          length in [mm], veloc in [mm/s], force in [N], acting on body 1
   CNTC.un_spck   = 1935; % code for selecting SIMPACK's unit convention,
                          %          length in [m], veloc in [m/s], force in [N], acting on body 2
   CNTC.un_si     = 1936; % code for selecting SI unit convention,
                          %          length in [m], veloc in [m/s], force in [N], acting on body 1
   CNTC.un_imper  = 1937; % code for selecting imperial units (n.y.a.),
                          %          length in [in], veloc in [???], force in [lbf], on body 1

   CNTC.ic_config = 1967; % code for setting the control digit C1, CONFIG:
                          %   0: wheelset on track, left rail
                          %   1: wheelset on track, right rail
                          %   4: wheelset on rollers, left side
                          %   5: wheelset on rollers, right side
   CNTC.ic_pvtime = 1970; % code for setting the control digit P, PVTIME:
                          %   0: continuation, sequence of cases (n.y.a. in module 1)
                          %   2: initiation of contact, first case
   CNTC.ic_bound  = 1971; % code for setting the control digit B, BOUND:
                          %   0: compute traction bound from normal problem
                          %   2, 3: set parabolic/elliptical traction bound (Fastsim)
                          %   4: set elliptical traction bound according to the SDEC approach
                          %   5: use non-Hertzian traction bound according to the KPEC approach
                          %   6: use non-Hertzian traction bound according to modified ANALYN approach
   CNTC.ic_tang   = 1972; % code for setting the control digit T, TANG:
                          %   0: frictionless, no tangential problem
                          %   1: shift or rolling with (material-) fixed coordinates
                          %   2: transient rolling with (contact-) moving coordinates
                          %   3: steady rolling with direct method, moving coordinates
   CNTC.ic_norm   = 1973; % code for setting the control digits N1 and N3, NORM:
                          % module 1: N1:  0: vertical position Z_WS prescribed
                          %                1: vertical force FZ prescribed
                          % module 3: N3:  0: approach PEN prescribed
                          %                1: normal force FN prescribed
   CNTC.ic_force  = 1974; % code for setting the control digit F, FORCE:
                          %   0: creepages CKSI and CETA prescribed (default)
                          %   1: tangential force FX and creepage CETA prescribed
                          %   2: tangential forces FX and FY prescribed
   CNTC.ic_frclaw = 1976; % code for retrieving the control digit L, FRCLAW
   CNTC.ic_discns = 1977; % code for setting the control digit D, DISCNS (module 1):
                          %   2: planar contact, using creep calculation
                          %   3: planar contact, extended rigid slip calculation
                          %   4: conformal contact on curved surface
   CNTC.ic_inflcf = 1978; % code for setting the control digit C3, INFLCF:
                          %   2: analytical IF for halfspace, piecewise constant
                          %   3: analytical IF for halfspace, bilinear elements
                          %   4: IF for conformal contact, using angle correction
   CNTC.ic_mater  = 1979; % code for retrieving the control digit M, MATER
   CNTC.ic_xflow  = 1981; % code for setting extended debug print output (X, XFLOW)
                          % <=0: no additional debug output
                          %  >0: codeword PSFLRIN for profil, smooth, force, locate,
                          %      readln, inflcf, nmdbg
   CNTC.ic_heat   = 1982; % code for retrieving the control digit H, HEAT (cntc_getflags)
   CNTC.ic_iestim = 1983; % code for setting the control digit I, IESTIM:
                          %   0: start from zero initial estimate (default)
                          %   3: use element division and tractions of previous case (n.y.a. in module 1)
   CNTC.ic_output = 1984; % code for setting the O-control digit O, OUTPUT:
                          %   0-4: amount of output on the results of the problem
   CNTC.ic_flow   = 1985; % code for setting the W-control digit W, FLOW:
                          %   0-9: amount of output on the flow of the computations
   CNTC.ic_return = 1986; % code for setting the R-control digit R, RETURN:
                          %   0-1: perform actual calculation
                          %   2-3: perform checks but skip actual calculation
   CNTC.ic_matfil = 1987; % code for setting the A-control digit A, MATFIL:
                          %   0: the mat-file is not created
                          %   1: write detailed results for contact area
                          %   2: write detailed results for potential contact area
   CNTC.ic_sens   = 1988; % code for setting the control digit S2, SENS:
                          %   0: sensitivities are not computed/printed
                          %   2: compute sensitivities for normal problem
                          %   3: compute sensitivities for normal+tangential problems
   CNTC.ic_ifmeth = 1989; % code for setting IF_METH for the Blanco-approach (C3=4)
                          %   0: fast, constant curvature (default),
                          %   1: detailed, varying curvature
   CNTC.ic_ifvari = 1990; % code for setting VARIANT for the Blanco-approach (C3=4)
                          %   variants: 1--4, default 4.
   CNTC.ic_sbsout = 1991; % code for setting the control digit O_s, OUTPUT_SUBS:
                          %   extent of output on subsurf.stresses to out-file:
                          %   0: no results are printed
                          %   1: print maximum values of primary stress invariants
                          %   2: more maximum values: Tresca, principal stresses
                          %   3: not used
                          %   4: print detailed results: displacements, invariants
   CNTC.ic_sbsfil = 1992; % code for setting the control digit A_s, MATFIL_SUBS:
                          %   0: the subs-file is not created
                          %   1: displacements and stress-invariants are written
                          %   2: the full stress tensor is written as well
   CNTC.ic_npomax = 1993; % code for setting the max. #elements in the potential contact area

   CNTC.if_idebug = 2000; % code for level of print-output of the addon itself
                          %     0 = none, 1 = overview (default), 2 = full, 3+ = debug
   CNTC.if_licdbg = 2001; % code for tracing which license is used
                          %     0 = none, 1 = overview (default), 2 = full, 3+ = debug
   CNTC.if_wrtinp = 2002; % code for activating writing of a CONTACT input-file
                          %     0 = no, 1 = write .inp-file (default)
   CNTC.if_openmp = 2003; % >0: number of threads to use per contact patch in CONTACT
                          %     -1 = automatic (#cores), 0 = no multi-threading (default)
                          % note: per-patch multi-threading is disabled if multi-
                          %       threading is used to compute different contact
                          %       patches concurrently

   % codes for cntc_getFieldData:

   CNTC.fld_h      =   1; % for retrieving array h      [length]
   CNTC.fld_mu     =   2; % for retrieving array mu     [-]
   CNTC.fld_px     =   3; % for retrieving array px     [force/area]
   CNTC.fld_py     =   4; % for retrieving array py     [force/area]
   CNTC.fld_pn     =   5; % for retrieving array pn     [force/area]
   CNTC.fld_ux     =   7; % for retrieving array ux     [length]
   CNTC.fld_uy     =   8; % for retrieving array uy     [length]
   CNTC.fld_un     =   9; % for retrieving array un     [length]
   CNTC.fld_taucrt =  11; % for retrieving array taucrt [force/area]
   CNTC.fld_uplsx  =  12; % for retrieving array uplsx  [length]
   CNTC.fld_uplsy  =  13; % for retrieving array uplsy  [length]
   CNTC.fld_sx     =  15; % for retrieving array sx     [-]
   CNTC.fld_sy     =  16; % for retrieving array sy     [-]
   CNTC.fld_temp1  =  20; % for retrieving array temp1  [C]
   CNTC.fld_temp2  =  21; % for retrieving array temp2  [C]
   CNTC.fld_wx     =  22; % for retrieving array wx     [-]
   CNTC.fld_wy     =  23; % for retrieving array wy     [-]

   % in module 3, the calling environment may provide additional information to be stored in the mat-file:

   CNTC.mt_tim    = 2011; % code for setting the simulation time
   CNTC.mt_xr     = 2012; % code for setting xr, long. pos. of contact point on rail
   CNTC.mt_yr     = 2013; % code for setting yr, lat. pos. of contact point on rail
   CNTC.mt_xw     = 2014; % code for setting xw, long. pos. of contact point on wheel
   CNTC.mt_yw     = 2015; % code for setting yw, lat. pos. of contact point on wheel
   CNTC.mt_sw     = 2016; % code for setting sw, long. wheel position in track coordinates
   CNTC.mt_run    = 2021; % code for setting irun, run number
   CNTC.mt_axle   = 2022; % code for setting iax, axle number
   CNTC.mt_side   = 2023; % code for setting iside, side number

   CNTC.err_allow  =  -12; % no license found or invalid license
   CNTC.err_search =  -25; % failure in contact search
   CNTC.err_ftot   =  -26; % no solution in total force iteration
   CNTC.err_norm   =  -27; % no convergence in NORM
   CNTC.err_tang   =  -28; % no convergence in TANG
   CNTC.err_tol    =  -29; % iteration stopped with residual > tolerance
   CNTC.err_icp    =  -31; % invalid ire/icp combination
   CNTC.err_profil =  -32; % error reading w/r profile file(s)
   CNTC.err_frclaw =  -33; % invalid friction parameters
   CNTC.err_other  =  -99; % other error, unspecified
   CNTC.err_broydn =  CNTC.err_ftot; % obsolete, kept for backward compatibility

end % function cntc_getmagicnumbers

%------------------------------------------------------------------------------------------------------------

