function [npatch, ws_pos, tot_forc, cp_pos, cp_creep, cp_force] = parse_out1(fname, ire_out, idebug)
%
% function [ sol ] = parse_out1(fname, [ire_out], [idebug])
%    or
% function [ npatch, ws_pos, tot_forc, cp_pos, cp_creep, cp_force, cp_subs ] = ...
%                                                               parse_out1(fname, [ire_out], [idebug])
%
% Reads CONTACT output for w/r contact cases (module 1) for result element 'ire_out' from out-file fname
% and returns the overall values from it:
%   npatch   = number of contact patches per case
%   ws_pos   = (input) wheel-set positions, orientations and velocities
%              struct with x, y, z, roll, yaw, pitch, vx, ... vpitch
%   tot_forc = total forces in track coordinates & wheelset coordinates
%              struct with fx_tr, ... fz_ws, mx_tr, ... mz_ws
%   cp_pos   = contact reference positions
%              struct with xtr, ytr, ztr, delttr, yr, zr, xw, yw, zw, ncon, nadh, nslip, nplast, prev_icp,
%              using npatch rows
%   cp_creep = list of creepages   
%              struct with pen, cksi, ceta, cphi, veloc - npatch rows
%   cp_force = forces and moments in contact reference coordinates
%              struct with fn, fx, fs, mn, el.en, fricw, pmax, temp - npatch rows
%   cp_subs  = maximum subsurface stresses and locations
%              [sighyd, sigvm, sigtr, sigma1, sigma3, sigxx, sigyy, sigzz]
%

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

if (nargin<2 | isempty(ire_out))
   ire_out = -1; % result element number as used in CONTACT library for which output is wanted, -1 == all
end
if (nargin<3 | isempty(idebug))
   idebug=0; % display input when idebug>=2, sections found when idebug>=5
end

% open file

if (isempty(strfind(fname,'.out')) & isempty(strfind(fname,'.ref_out')))
   fname=[deblank(fname),'.out'];
end
if (~exist(fname,'file'))
   disp(sprintf('\nERROR: cannot find file %s\n',fname))
   npatch=[]; ws_pos=[]; tot_forc=[]; cp_pos=[]; cp_creep=[]; cp_force=[];
   return
end

% list the special lines that precede the output data

headers = strvcat('Case', ...
                  'WHEEL-RAIL CONTACT', ...
                  'WHEEL-ROLLER CONTACT', ...
                  'WHEEL-SET POSITION', ...
                  'FLEXIBLE WHEEL-SET DEVIATIONS', ...
                  'RAIL IRREGULARITY', ...
                  'MASSLESS RAIL DEFLECTION', ...
                  'TOTAL FORCES AND MOMENTS', ...
                  'AVERAGE CONTACT POSITION', ...
                  'FX(TR)      FY(TR)', ...
                  'ELAST.EN.   FRIC', ...
                  'DATA FOR CONTACT PATCH', ...
                  'CONTACT REFERENCE LOCATION', ...
                  'KINEMATIC CONSTANTS', ...
                  'TOTAL FORCES, TORSIONAL', ...
                  'CONTACT STATISTICS', ...
                  'ABSMAX SIGHYD');

subs_keys = strvcat('SIGHYD', 'SIGVM', 'SIGTR', 'SIGMA1', 'SIGMA3', 'SIGXX', 'SIGYY', 'SIGZZ');
nsubs_keys = size(subs_keys,1);

% read contents of file into memory

f = read_file(fname, idebug);
f = find_all_headers(f, headers, idebug);

% read first/second lines for checking that it is a .out-file.

iline = 0;
[s, iline] = read_line(f, iline, idebug);
if (~isempty(strfind(s, 'refresh the license')))
   [s, iline] = read_line(f, iline, idebug);
end
if (~isempty(strfind(s, 'w/r contact on result element')))
   % CONTACT library: case starts at line 1
   iline = 0;
else
   [s, iline] = read_line(f, iline, idebug);
   if (isempty(strfind(s, 'detailed investigation of 3D')))
      disp('ERROR: this file doesn''t look like a CONTACT .out-file.');
      return;
   end
end

% parse the contents of the file

icase  = 0;
max_patch = 10;
ire    = []; jcase  = []; tim    = []; npatch  = []; 
x_ws   = []; y_ws   = []; z_ws   = []; roll    = []; yaw    = []; pitch  = [];
vx_ws  = []; vy_ws  = []; vz_ws  = []; vroll   = []; vyaw   = []; vpitch = [];
dx_whl = []; dy_whl = []; dz_whl = []; drollw  = []; dyaww  = []; dpitchw = [];
vxwhl  = []; vywhl  = []; vzwhl  = []; vrollw  = []; vyaww  = []; vpitchw = [];
dyrail = []; dzrail = []; drrail = []; vyrail  = []; vzrail = []; vrrail = [];
kyrail = []; dydefl = []; fydefl = []; kzrail  = []; dzdefl = []; fzdefl = [];
fx_tr  = []; fy_tr  = []; fz_tr  = []; fx_ws   = []; fy_ws  = []; fz_ws  = [];
xcp_tr = []; ycp_tr = []; zcp_tr = []; delt_tr = []; ycp_r  = []; zcp_r  = [];
xcp_w  = []; ycp_w  = []; zcp_w  = [];
pen    = []; cksi   = []; ceta   = []; cphi    = []; veloc  = [];
fn_loc = []; fx_loc = []; fs_loc = []; mn_loc  = []; elen   = []; fric   = [];
temp1  = []; temp2  = []; pmax   = []; 
ncon   = []; nadh   = []; nslip  = []; nplast  = []; prevcp = [];
subs_max = [];

while (iline<f.nline)

   % skip lines until one of the headers is found

   [ s, iline, ihdr ] = next_header(f, iline, idebug);
   if (ihdr<0), break; end

   % Process the consecutive section

   if (~isempty(strfind(s, 'Case')))

      icase = icase + 1;
      if (idebug>=5)
         disp(sprintf('...Case %d starts at line %d',icase,iline));
      end

      % Case     1
      % Case  29246 for w/r contact on result element  2, t=      0.08706221

      tmp = sscanf(s,' Case %d for w/r contact on result element %d, t= %f');
      jcase(icase) = tmp(1);
      if (~isempty(strfind(s, 'result element')))
         ire(icase) = tmp(2);
      end
      if (~isempty(strfind(s, ', t=')))
         tim(icase) = tmp(3);
      end

      % initialize the patch number for module 3

      ipatch = 1;

      % initialize output-values

      npatch(icase)  = 0;
      x_ws(icase)    = NaN; y_ws(icase)    = NaN; z_ws(icase)    = NaN; 
      roll(icase)    = NaN; yaw(icase)     = NaN; pitch(icase)   = NaN;
      vx_ws(icase)   = NaN; vy_ws(icase)   = NaN; vz_ws(icase)   = NaN; 
      vroll(icase)   = NaN; vyaw(icase)    = NaN; vpitch(icase)  = NaN;
      dx_whl(icase)  = NaN; dy_whl(icase)  = NaN; dz_whl(icase)  = NaN; 
      drollw(icase)  = NaN; dyaww(icase)   = NaN; dpitchw(icase) = NaN;
      vxwhl(icase)   = NaN; vywhl(icase)   = NaN; vzwhl(icase)   = NaN;
      vrollw(icase)  = NaN; vyaww(icase)   = NaN; vpitchw(icase) = NaN;
      fx_tr(icase)   = NaN; fy_tr(icase)   = NaN; fz_tr(icase)   = NaN; 
      fx_ws(icase)   = NaN; fy_ws(icase)   = NaN; fz_ws(icase)   = NaN;
      mx_tr(icase)   = NaN; my_tr(icase)   = NaN; mz_tr(icase)   = NaN; 
      mx_ws(icase)   = NaN; my_ws(icase)   = NaN; mz_ws(icase)   = NaN;

      al = [1:max_patch];
      xcp_tr(al,icase) = NaN; ycp_tr(al,icase) = NaN; zcp_tr(al,icase) = NaN;
      delt_tr(al,icase)= NaN; ycp_r(al,icase)  = NaN; zcp_r(al,icase)  = NaN;
      xcp_w(al,icase)  = NaN; ycp_w(al,icase)  = NaN; zcp_w(al,icase)  = NaN;
      pen(al,icase)    = NaN; veloc(al,icase)  = NaN;
      cksi(al,icase)   = NaN; ceta(al,icase)   = NaN; cphi(al,icase)   = NaN;
      fn_loc(al,icase) = NaN; fx_loc(al,icase) = NaN; fs_loc(al,icase) = NaN; 
      mn_loc(al,icase) = NaN; elen(al,icase)   = NaN; fric(al,icase)   = NaN;
      temp1(al,icase)  = NaN; temp2(al,icase)  = NaN; pmax(al,icase)   = NaN;
      ncon(al,icase)   = NaN; nadh(al,icase)   = NaN; nslip(al,icase)  = NaN;
      nplast(al,icase) = NaN; prevcp(al,al,icase) = 0;

      if (~isempty(subs_max))
         subs_max(1:4, al, icase, 1:nsubs_keys) = NaN;
      end
   end

   if (~isempty(strfind(s, 'WHEEL-RAIL CONTACT')) | ~isempty(strfind(s, 'WHEEL-ROLLER CONTACT')))

      % WHEEL-RAIL CONTACT, LEFT WHEEL,  1 CONTACT PATCH

      if (idebug>=5), disp(sprintf('...Found number of patches at line %d',iline)); end
      npatch(icase) = sscanf(s, 'WHEEL-R%*s CONTACT, %*s WHEEL, %d');
   end

   if (~isempty(strfind(s, 'WHEEL-SET POSITION')))

      % WHEEL-SET POSITION AND VELOCITY
      %     X_WS        Y_WS        Z_WS        ROLL        YAW         PITCH
      %     0.000       0.000      0.2154       0.000       0.000       0.000    
      % 
      %     VX_WS       VY_WS       VZ_WS       VROLL       VYAW        VPITCH
      %     40.00       1.000E-02   0.000       0.000       0.000      -86.93    

      if (idebug>=5)
         disp(sprintf('...Found wheel-set position and velocity at line %d',iline));
      end
      [s, iline] = read_line(f, iline, idebug);
      [s, iline] = read_line(f, iline, idebug);
      tmp = sscanf(s, '%f %f %f %f %f %f');
      x_ws(icase)   = tmp(1); y_ws(icase)   = tmp(2); z_ws(icase)   = tmp(3);
      roll(icase)   = tmp(4); yaw(icase)    = tmp(5); pitch(icase)  = tmp(6);

      [s, iline] = read_line(f, iline, idebug);
      [s, iline] = read_line(f, iline, idebug);
      [s, iline] = read_line(f, iline, idebug);
      tmp = sscanf(s, '%f %f %f %f %f %f');
      vx_ws(icase)  = tmp(1); vy_ws(icase)  = tmp(2); vz_ws(icase)  = tmp(3);
      vroll(icase)  = tmp(4); vyaw(icase)   = tmp(5); vpitch(icase) = tmp(6);
   end

   if (~isempty(strfind(s, 'FLEXIBLE WHEEL-SET DEVIATIONS')))

      % FLEXIBLE WHEEL-SET DEVIATIONS
      %     DXWHL       DYWHL       DZWHL       DROLLW      DYAWW       DPITCHW
      %         0.000       0.000       0.000  -4.000E-02       0.000       0.000
      %     VXWHL       VYWHL       VZWHL       VROLLW      VYAWW       VPITCHW
      %         0.000       0.000       0.000       0.000       0.000       0.000

      if (idebug>=5)
         disp(sprintf('...Found flexible wheel-set deviations at line %d',iline));
      end
      [s, iline] = read_line(f, iline, idebug);
      [s, iline] = read_line(f, iline, idebug);
      tmp = sscanf(s, '%f %f %f %f %f %f');
      dx_whl(icase) = tmp(1); dy_whl(icase) = tmp(2); dz_whl(icase)  = tmp(3);
      drollw(icase) = tmp(4); dyaww(icase)  = tmp(5); dpitchw(icase) = tmp(6);

      [s, iline] = read_line(f, iline, idebug);
      [s, iline] = read_line(f, iline, idebug);
      tmp = sscanf(s, '%f %f %f %f %f %f');
      vxwhl(icase)  = tmp(1); vywhl(icase)  = tmp(2); vzwhl(icase)   = tmp(3);
      vrollw(icase) = tmp(4); vyaww(icase)  = tmp(5); vpitchw(icase) = tmp(6);
   end

   if (~isempty(strfind(s, 'RAIL IRREGULARITY')))

      % RAIL IRREGULARITY
      %     DYRAIL      DZRAIL      DROLLR      VYRAIL      VZRAIL      VROLLR
      %    -2.773E-02      0.2268       0.000      -213.1       0.000       0.000

      if (idebug>=5)
         disp(sprintf('...Found rail irregularities at line %d',iline));
      end
      [s, iline] = read_line(f, iline, idebug);
      [s, iline] = read_line(f, iline, idebug);
      tmp = sscanf(s, '%f %f %f %f %f %f');
      dyrail(icase) = tmp(1); dzrail(icase) = tmp(2); drrail(icase) = tmp(3);
      vyrail(icase) = tmp(4); vzrail(icase) = tmp(5); vrrail(icase) = tmp(6);
   end

   if (~isempty(strfind(s, 'MASSLESS RAIL DEFLECTION')))

      % MASSLESS RAIL DEFLECTION
      %     KY_RAIL     DY_DEFL     FY_RAIL     KZ_RAIL     DZ_DEFL     FZ_RAIL
      %      100.0000     -1.0000   1.162E+05      0.0000      0.0000   1.230E+05

      if (idebug>=5)
         disp(sprintf('...Found massless rail deflection at line %d',iline));
      end
      [s, iline] = read_line(f, iline, idebug);
      [s, iline] = read_line(f, iline, idebug);
      tmp = sscanf(s, '%f %f %f %f %f %f');
      kyrail(icase) = tmp(1); dydefl(icase) = tmp(2); fydefl(icase) = tmp(3);
      kzrail(icase) = tmp(4); dzdefl(icase) = tmp(5); fzdefl(icase) = tmp(6);
   end

   if (~isempty(strfind(s, 'TOTAL FORCES AND MOMENTS')) | ...
       ~isempty(strfind(s, 'FX(TR)      FY(TR)')) )

      % TOTAL FORCES AND MOMENTS ON RAIL | ROLLER
      %     FX(TR)      FY(TR)      FZ(TR)      FX(WS)      FY(WS)      FZ(WS)
      %       70.09      -40.92   1.000E+04       70.09      -40.92   1.000E+04
      %
      %     MX(TR)      MY(TR)      MZ(TR)      MX(WS)      MY(WS)      MZ(WS)
      %   5.170E+07   1.092E+06  -1.812E+07   5.544E+07   1.492E+06  -1.944E+07

      if (idebug>=5)
         disp(sprintf('...Found global forces and moments at line %d',iline));
      end
      % support for compact .out-files without 'TOTAL FORCES...'
      if (isempty(strfind(s, 'FX(TR)      FY(TR)')))
         [s, iline] = read_line(f, iline, idebug);
      end
      [s, iline] = read_line(f, iline, idebug);
      tmp = sscanf(s, '%f %f %f %f %f %f');
      fx_tr(icase) = tmp(1); fy_tr(icase) = tmp(2); fz_tr(icase) = tmp(3);
      fx_ws(icase) = tmp(4); fy_ws(icase) = tmp(5); fz_ws(icase) = tmp(6);

      % check for MX(TR) on following lines
      [s, iline] = read_line(f, iline, idebug);
      [s, iline] = read_line(f, iline, idebug);
      if (~isempty(strfind(s, 'MX(TR)      MY(TR)')))
         [s, iline] = read_line(f, iline, idebug);
         tmp = sscanf(s, '%f %f %f %f %f %f');
         mx_tr(icase) = tmp(1); my_tr(icase) = tmp(2); mz_tr(icase) = tmp(3);
         mx_ws(icase) = tmp(4); my_ws(icase) = tmp(5); mz_ws(icase) = tmp(6);
      else
         iline = iline - 2;
         [s, iline] = read_line(f, iline, idebug);
      end
   end

   if (~isempty(strfind(s, 'DATA FOR CONTACT PATCH')))

      % DATA FOR CONTACT PATCH  1 (NEW PATCH)
      % DATA FOR CONTACT PATCH  1 (PREV  1,  2)

      if (idebug>=5)
         disp(sprintf('...Found contact patch number at line %d',iline));
      end
      ipatch = sscanf(strtrim(s), '----- DATA FOR CONTACT PATCH %d');

      if (~isempty(strfind(s, 'NEW PATCH')))
         prevcp(ipatch,:,icase) = 0;
      else
         [tmp, nval] = sscanf(strtrim(s), '----- DATA FOR CONTACT PATCH %*d (PREV %d, %d, %d, %d)');
         prevcp(ipatch,1:nval,icase) = tmp;
      end

      % support for compact .out-files without 'WHEEL-RAIL CONTACT...'
      if (ipatch>npatch(icase))
         npatch(icase) = ipatch;
      end
   end

   if (~isempty(strfind(s, 'CONTACT REFERENCE LOCATION')))

      % CONTACT REFERENCE LOCATION
      %     XCP(TR)     YCP(TR)     ZCP(TR)    DELTCP(TR)   YCP(R)      ZCP(R)
      %      0.0000   -751.0869      0.1489      0.0312     -9.4445      0.1490
      %
      %     XCP(W)      YCP(W)      ZCP(W)
      %      0.0000      1.0869     -0.0491

      if (idebug>=5)
         disp(sprintf('...Found contact point location at line %d',iline));
      end
      [s, iline] = read_line(f, iline, idebug);
      [s, iline] = read_line(f, iline, idebug);
      tmp = sscanf(s, '%f %f %f %f %f %f');
      xcp_tr(ipatch,icase) = tmp(1); ycp_tr(ipatch,icase)  = tmp(2);
      zcp_tr(ipatch,icase) = tmp(3); delt_tr(ipatch,icase) = tmp(4);
      ycp_r(ipatch,icase)  = tmp(5); zcp_r(ipatch,icase)   = tmp(6);

      [s, iline] = read_line(f, iline, idebug);
      [s, iline] = read_line(f, iline, idebug);
      [s, iline] = read_line(f, iline, idebug);
      tmp = sscanf(s, '%f %f %f');
      xcp_w(ipatch,icase) = tmp(1); ycp_w(ipatch,icase) = tmp(2);
      zcp_w(ipatch,icase) = tmp(3);
   end

   if (~isempty(strfind(s, 'KINEMATIC CONSTANTS')))

      % KINEMATIC CONSTANTS
      %     CHI         DQ          VELOC       CKSI        CETA        CPHI
      %     0.000       1.000           1.000   7.282E-04   7.938E-04  -1.736E-04

      % get the creepages of a contact patch

      if (idebug>=5)
         disp(sprintf('...Found kinematic constants at line %d',iline));
      end
      [s1, iline] = read_line(f, iline, idebug);
      [s2, iline] = read_line(f, iline, idebug);
      tmp = sscanf(s2, '%f %f %f %f %f %f');
      veloc(ipatch,icase)  = tmp(3);
      cksi(ipatch,icase)   = tmp(4);
      ceta(ipatch,icase)   = tmp(5);
      cphi(ipatch,icase)   = tmp(6);
   end

   if (~isempty(strfind(s, 'TOTAL FORCES, TORSIONAL')) | ...
       ~isempty(strfind(s, 'ELAST.EN.   FRIC')) )

      % T=1, F=2:
      % TOTAL FORCES, TORSIONAL MOMENT, ELASTIC ENERGY, FRICTIONAL POWER
      %     FN          FX          FS          MN         ELAST.EN.  FRIC.POWER
      %     25.00      0.1885       0.000       0.000       4.736       0.000
      %     FN/G       SHIFT X      SHIFT Y    APPROACH       PMAX
      %     25.00      -2.708E-03  -1.106E-15  6.492E-03    12.34

      % get the creepages or forces

      if (idebug>=5)
         disp(sprintf('...Found contact patch forces at line %d',iline));
      end
      % support for compact .out-files without 'TOTAL FORCES...'
      if (isempty(strfind(s, 'ELAST.EN.   FRIC'))) 
         [s , iline] = read_line(f, iline, idebug);
      end
      [s1, iline] = read_line(f, iline, idebug);
      tmp = sscanf(s1, '%f %f %f %f %f %f');
      fn_loc(ipatch,icase)    = tmp(1);
      fx_loc(ipatch,icase)    = tmp(2); % storing absolute force
      fs_loc(ipatch,icase)    = tmp(3); % storing absolute force
      mn_loc(ipatch,icase)    = tmp(4);
      elen(ipatch,icase)      = tmp(5);
      fric(ipatch,icase)      = tmp(6);

      [s2, iline] = read_line(f, iline, idebug);
      [s3, iline] = read_line(f, iline, idebug);
      tmp = sscanf(s3, '%f %f %f %f %f %f');
      if (~isempty([strfind(s2, 'CREEP X'), strfind(s2, 'SHIFT X')]))
         cksi(ipatch,icase)   = tmp(2);
      % else
      %    fx_loc(ipatch,icase) = tmp(2); % storing relative force
      end
      if (~isempty([strfind(s2, 'CREEP Y'), strfind(s2, 'SHIFT Y')]))
         ceta(ipatch,icase)   = tmp(3);
      % else
      %    fs_loc(ipatch,icase) = tmp(3); % storing relative force
      end
      pen(ipatch,icase) = tmp(4);
      if (~isempty(strfind(s2, 'PMAX')))
         pmax(ipatch,icase)  = tmp(5);              % >= v22.1
      end
      if     (~isempty(strfind(s2, 'MAX(T1)')))     % <= v21.1
         temp1(ipatch,icase) = tmp(5);
         temp2(ipatch,icase) = tmp(6);
      elseif (~isempty(strfind(s2, 'MAX(T1,T2)')))  % >= v22.1
         temp1(ipatch,icase) = tmp(6);
      end
   end % TOTAL FORCES, TORSIONAL

   if (~isempty(strfind(s, 'CONTACT STATISTICS')))

      % CONTACT STATISTICS
      % N: NUMBER OF ELEMENTS IN REGION,  I: NUMBER OF ITERATIONS,
      % POT: POTENTIAL CONTACT AREA,  CON: CONTACT,  ADH: ADHESION
      %    NPOT   NCON   NADH  NSLIP  INORM  ITANG
      %   3731   1286    919    367      1      1
      % or
      %    NPOT   NCON   NADH  NSLIP  NPLAST INORM  ITANG
      %      30     26     14      8      4      1      1

      % get the number of elements in contact/adhesion/slip

      if (idebug>=5)
         disp(sprintf('...Found contact statistics at line %d',iline));
      end
      [s , iline] = read_line(f, iline, idebug);
      [s , iline] = read_line(f, iline, idebug);
      [s , iline] = read_line(f, iline, idebug);
      [s1, iline] = read_line(f, iline, idebug);
      tmp = sscanf(s1, '%f %f %f %f %f %f');
      ncon(ipatch,icase)    = tmp(2);
      nadh(ipatch,icase)    = tmp(3);
      nslip(ipatch,icase)   = tmp(4);
      if (~isempty(strfind(s, 'NPLAST')))
         nplast(ipatch,icase) = tmp(5);
      end
   end % STATISTICS

   if (~isempty(strfind(s, 'ABSMAX SIGHYD')))

      if (isempty(subs_max))
         subs_max(1:4, 1:max_patch, 1:icase, 1:nsubs_keys) = NaN;
      end

      % ABSMAX SIGHYD =    -0.853 AT (X,Y,Z) = (   0.000,   0.000,   0.000)
      % MAX     SIGVM =     0.665 AT (X,Y,Z) = (   0.000,   0.000,   0.400)
      % MAX     SIGTR =     0.665 AT (X,Y,Z) = (   0.000,   0.000,   0.400)
      % MAX    SIGMA1 =     0.008 AT (X,Y,Z) = (   0.000,   0.000,   1.250)
      % MIN    SIGMA3 =    -1.000 AT (X,Y,Z) = (   0.000,   0.000,   0.000)
      % ABSMAX  SIGXX =    -0.780 AT (X,Y,Z) = (   0.000,   0.000,   0.000)
      % ABSMAX  SIGYY =    -0.780 AT (X,Y,Z) = (   0.000,   0.000,   0.000)
      % ABSMAX  SIGZZ =    -1.000 AT (X,Y,Z) = (   0.000,   0.000,   0.000)
 
      % process all lines containing a maximum value

      while(iline<=f.nline & ~isempty(strfind(s, 'AT (X,Y,Z)')))

         % get the key number
         key = strtrim(s(8:14));
         ikey = 1;
         while(ikey<nsubs_keys & ~strcmp(key, deblank(subs_keys(ikey,:))))
             ikey = ikey + 1;
         end
         if (~strcmp(key, deblank(subs_keys(ikey,:))))
            disp(sprintf('Internal error: unknown key "%s".', key));
            break
         end

         % get the values
         tmp = sscanf(s(17:end),' %f AT (X,Y,Z) = ( %f, %f, %f)');

         % store in subs_max
         subs_max(1:4, ipatch, icase, ikey) = tmp;

         % read next line from file
         [ s, iline ] = read_line(f, iline, idebug);
      end

   end % ABSMAX, SIGHYD

end % while (iline<f.nline)

% select cases for result element ire_out, if provided

if (ire_out>0 & ~isempty(ire))
   ncase    = length(jcase);
   i_select = find(ire == ire_out);
   disp(sprintf('Total #cases = %d, selected %d cases for ire_out = %d', ncase, length(i_select), ire_out));

   npatch = npatch(i_select);
   ire    = ire   (i_select); jcase  = jcase (i_select); tim     = tim    (i_select);
   x_ws   = x_ws  (i_select); y_ws   = y_ws  (i_select); z_ws    = z_ws   (i_select); 
   roll   = roll  (i_select); yaw    = yaw   (i_select); pitch   = pitch  (i_select);
   vx_ws  = vx_ws (i_select); vy_ws  = vy_ws (i_select); vz_ws   = vz_ws  (i_select); 
   vroll  = vroll (i_select); vyaw   = vyaw  (i_select); vpitch  = vpitch (i_select);
   dx_whl = dx_whl(i_select); dy_whl = dy_whl(i_select); dz_whl  = dz_whl (i_select); 
   drollw = drollw(i_select); dyaww  = dyaww (i_select); dpitchw = dpitchw(i_select);
   vxwhl  = vxwhl (i_select); vywhl  = vywhl (i_select); vzwhl   = vzwhl  (i_select);
   vrollw = vrollw(i_select); vyaww  = vyaww (i_select); vpitchw = vpitchw(i_select);
   fx_tr  = fx_tr (i_select); fy_tr  = fy_tr (i_select); fz_tr   = fz_tr  (i_select); 
   fx_ws  = fx_ws (i_select); fy_ws  = fy_ws (i_select); fz_ws   = fz_ws  (i_select);
   mx_tr  = mx_tr (i_select); my_tr  = my_tr (i_select); mz_tr   = mz_tr  (i_select); 
   mx_ws  = mx_ws (i_select); my_ws  = my_ws (i_select); mz_ws   = mz_ws  (i_select);
 
   xcp_tr  = xcp_tr (:,i_select); ycp_tr  = ycp_tr(:,i_select); zcp_tr  = zcp_tr(:,i_select);
   delt_tr = delt_tr(:,i_select); ycp_r   = ycp_r (:,i_select); zcp_r   = zcp_r (:,i_select);
   xcp_w   = xcp_w  (:,i_select); ycp_w   = ycp_w (:,i_select); zcp_w   = zcp_w (:,i_select);
   pen     = pen    (:,i_select); veloc   = veloc (:,i_select);
   cksi    = cksi   (:,i_select); ceta    = ceta  (:,i_select); cphi    = cphi  (:,i_select);
   fn_loc  = fn_loc (:,i_select); fx_loc  = fx_loc(:,i_select); fs_loc  = fs_loc(:,i_select);
   mn_loc  = mn_loc (:,i_select); elen    = elen  (:,i_select); fric    = fric  (:,i_select);
   temp1   = temp1  (:,i_select); temp2   = temp2 (:,i_select); pmax    = pmax  (:,i_select);
   ncon    = ncon   (:,i_select); nadh    = nadh  (:,i_select); nslip   = nslip (:,i_select);
   nplast  = nplast (:,i_select); prevcp  = prevcp(:,:,i_select);
   if (~isempty(subs_max))
      subs_max = subs_max(:,:,i_select,:);
   end
end

% group arrays into structures

ws_pos   = struct('x',     x_ws,   'y',     y_ws,   'z',     z_ws,   ...
                  'roll',  roll,   'yaw',   yaw,    'pitch', pitch,  ...
                  'vx',    vx_ws,  'vy',    vy_ws,  'vz',    vz_ws,  ...
                  'vroll', vroll,  'vyaw',  vyaw,   'vpitch',vpitch);
if (any(~isnan(dx_whl)))
   ws_pos.dx_whl = dx_whl; ws_pos.dy_whl = dy_whl; ws_pos.dz_whl  = dz_whl; 
   ws_pos.drollw = drollw; ws_pos.dyaww  = dyaww;  ws_pos.dpitchw = dpitchw;
   ws_pos.vxwhl  = vxwhl;  ws_pos.vywhl  = vywhl;  ws_pos.vzwhl   = vzwhl; 
   ws_pos.vrollw = vrollw; ws_pos.vyaww  = vyaww;  ws_pos.vpitchw = vpitchw;
end
if (any(~isnan(dyrail)))
   ws_pos.dyrail = dyrail; ws_pos.dzrail = dzrail;  ws_pos.drrail  = drrail; 
   ws_pos.vyrail = vyrail; ws_pos.vzrail = vzrail;  ws_pos.vrrail  = vrrail; 
end
if (any(~isnan(kyrail)))
   ws_pos.kyrail = kyrail; ws_pos.dydefl = dydefl;  ws_pos.fydefl  = fydefl; 
   ws_pos.kzrail = kzrail; ws_pos.dzdefl = dzdefl;  ws_pos.fzdefl  = fzdefl; 
end
tot_forc = struct('fx_tr', fx_tr,  'fy_tr', fy_tr,  'fz_tr', fz_tr,  ...
                  'fx_ws', fx_ws,  'fy_ws', fy_ws,  'fz_ws', fz_ws,  ...
                  'mx_tr', mx_tr,  'my_tr', my_tr,  'mz_tr', mz_tr,  ...
                  'mx_ws', mx_ws,  'my_ws', my_ws,  'mz_ws', mz_ws);

max_patch = max(npatch);
al        = [1:max_patch];
cp_pos    = struct('xcp_tr', xcp_tr(al,:),  'ycp_tr',  ycp_tr(al,:),  ...
                   'zcp_tr', zcp_tr(al,:),  'delt_tr', delt_tr(al,:), ...
                   'ycp_r',  ycp_r(al,:),   'zcp_r',   zcp_r(al,:),   ...
                   'xcp_w',  xcp_w(al,:),   'ycp_w',   ycp_w(al,:),   ...
                   'zcp_w',  zcp_w(al,:),   'ncon',    ncon(al,:),    ...
                   'nadh',   nadh(al,:),    'nslip',   nslip(al,:),   ...
                   'nplast', nplast(al,:),  'prev_icp', prevcp(al,al,:));
cp_creep  = struct('pen',    pen(al,:),     'cksi',    cksi(al,:),    ...
                   'ceta',   ceta(al,:),    'cphi',    cphi(al,:),    ...
                   'veloc',  veloc(al,:));
cp_force  = struct('fn',     fn_loc(al,:),  'fx',      fx_loc(al,:),  ...
                   'fs',     fs_loc(al,:),  'mn',      mn_loc(al,:),  ...
                   'elen',   elen(al,:),    'fric',    fric(al,:));
if (~all(isnan(pmax(al,:))))
   cp_force.pmax  = pmax(al,:);
end
if (~all(isnan(temp1(al,:))) & all(isnan(temp2(al,:))))
   cp_force.temp  = temp1(al,:);  % >= v22.1
elseif (~all(isnan(temp1(al,:))))
   cp_force.temp1 = temp1(al,:);  % <= v21.1
   cp_force.temp2 = temp2(al,:);
end

if (isempty(subs_max))
   cp_subs   = [];
else
   cp_subs   = struct('sighyd', squeeze(subs_max(:,al,:,1)), ...
                      'sigvm',  squeeze(subs_max(:,al,:,2)), ...
                      'sigtr',  squeeze(subs_max(:,al,:,3)), ...
                      'sigma1', squeeze(subs_max(:,al,:,4)), ...
                      'sigma3', squeeze(subs_max(:,al,:,5)), ...
                      'sigxx',  squeeze(subs_max(:,al,:,6)), ...
                      'sigyy',  squeeze(subs_max(:,al,:,7)), ...
                      'sigzz',  squeeze(subs_max(:,al,:,8)));
end

% convert separate structs into one overall struct, if requested

if (nargout<=1)
   sol = struct( 'ws_pos',ws_pos, 'tot_forc',tot_forc, 'npatch',  npatch, ...
                 'cp_pos',cp_pos, 'cp_creep',cp_creep, 'cp_force',cp_force,...
                 'cp_subs',cp_subs);
   if (~isempty(ire)), sol.ire = ire; end
   if (~isempty(tim)), sol.tim = tim; end
   sol.icase = jcase;
   npatch = sol;
end

end % function parse_out1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read_file: Helper function for reading entire file into memory

function [ myfile ] = read_file(fname, idebug)

% function [ file_contents ] = read_line(fname, idebug)

myfile.fname = fname;
myfile.nline = 0;
myfile.contents = {};

f = fopen(fname);
while(~feof(f))
   myfile.nline = myfile.nline + 1; 
   myfile.contents{myfile.nline} = fgets(f);
end
fclose(f);

end % function read_file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find_all_headers: Helper function for locating headers in file

function [ f ] = find_all_headers(f, headers, idebug)

% function [ f ] = find_all_headers(f, headers, idebug)

f.headers = zeros(f.nline,1);

num_hdr = size(headers,1);
for ihdr = 1 : num_hdr
   ix = strfind(f.contents, deblank(headers(ihdr,:)));
   for iline = 1 : f.nline
      if (~isempty(ix{iline}))
         f.headers(iline) = ihdr;
      end
   end
end

end % function find_all_headers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read_line: Helper function for getting next line of input

function [ s, iline ] = read_line(f, iline, idebug)

% function [ s, iline ] = read_line(f, iline, idebug)

iline = iline + 1; 
s = f.contents{iline};

if (idebug>=2),
   % debug: print line without trailing \n (?)
   if (length(s)>=2 & double(s(end-1))==13)
      disp(sprintf('line %4d: %s', iline, s(1:end-2) ))
   elseif (length(s)>=1 & double(s(end))==10)
      disp(sprintf('line %4d: %s', iline, s(1:end-1) ))
   else
      disp(sprintf('line %4d: %s', iline, s))
   end
end

end % function read_line

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% next_header: Helper function for locating the next section

function [ s, iline, ihdr ] = next_header(f, iline, idebug)

% function [ s, iline, ihdr ] = next_header(f, iline, idebug)

iprev = iline;
inext = iline+1;
while(inext<f.nline & f.headers(inext)<=0)
   inext = inext + 1;
end

if (f.headers(inext)>0)
   if (idebug>=6)
      disp(sprintf('...Next_header: inext = %d, fast-forword to iline = %d',inext,iline+inext));
   end
   iline = inext;
   s     = strtrim(f.contents{iline});
   ihdr  = f.headers(iline);
else
   if (idebug>=6)
      disp(sprintf('...Next_header: no more headers'));
   end
   iline = f.nline;
   s     = [];
   ihdr  = -1;
end

if (idebug>=2),
   for jline = iprev+1 : iline
      sloc = f.contents{jline};
      % debug: print line without trailing \n
      if (length(sloc)>=2 & double(sloc(end-1))==13)
         disp(sprintf('line %4d: %s', jline, sloc(1:end-2) ))
      elseif (length(sloc)>=1 & double(sloc(end))==10)
         disp(sprintf('line %4d: %s', jline, sloc(1:end-1) ))
      else
         disp(sprintf('line %4d: %s', jline, sloc))
      end
   end
end

end % function next_header

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
