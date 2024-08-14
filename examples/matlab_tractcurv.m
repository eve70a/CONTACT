
%============================================================================================================
% Example program using CONTACT library for Matlab. See User guide, Sections 7.5 & 5.8.
% Calculation of creep--force curves.
%============================================================================================================
% Part 1: Initialize the CONTACT library, register problem "iwhe = 1"
%============================================================================================================

if (~exist('cntc_initlibrary.m','file'))
   % set location of CONTACT installation folder
   % contactdir = 'C:\Program Files\Vtech CMCC\contact_v24.1';
   contactdir = '..';
   addpath([contactdir, '\matlab_intfc']);
   addpath([contactdir, '\matlab']);
end

if (exist('cntc_initlibrary')~=2)
   disp('ERROR: cant find CONTACT library on the Matlab search path');
   return
end
[CNTC, ifcver, ierror] = cntc_initlibrary;

if (ierror<0)
   disp(sprintf('An error occurred, ierror = %d, check output-file.',ierror));
   return
end

idebug = 1;             % 1: just a bit of information from the library
cntc_setglobalflags(CNTC.if_idebug, idebug);

iwhe   = 1; % "wheel number"
icp    = 1; % "contact problem on wheel"
imodul = 3; % basic contact
[ifcver, ierror] = cntc_initialize(iwhe, imodul);

%============================================================================================================
% Part 2: Configure the main flags & control digits of the contact problem
%============================================================================================================
 
            % CONTACT unit convention: [mm], [mm/s], [N], acting on body (1)
flags = []; values = [];
flags(1) = CNTC.if_units ; values(1) = CNTC.un_cntc;
flags(2) = CNTC.if_wrtinp; values(2) = 0;    %   0, no .inp-file needed
flags(3) = CNTC.ic_matfil; values(3) = 0;    % A=0, no .mat-file needed
flags(4) = CNTC.ic_output; values(4) = 3;    % O=1, min. output to .out-file
flags(5) = CNTC.ic_flow  ; values(5) = 4;    % W=1, a little of progress output
flags(6) = CNTC.ic_norm  ; values(6) = 1;    % N=1, FN prescribed
flags(7) = CNTC.ic_tang  ; values(7) = 3;    % T=3, steady state rolling
flags(8) = CNTC.ic_pvtime; values(8) = 2;    % P=2, no previous time

cntc_setflags(iwhe, icp, flags, values);

% G=0: set maxgs, maxin, maxnr, maxout, eps

gdigit  = 0;
values  = [ 80, 200, 5, 1 ];
rvalues = [ 1d-6 ];
cntc_setsolverflags(iwhe, icp, gdigit, values, rvalues);

% Hertzian discretization, automatic potential contact area
% ellipticity a/b = 0.5

ipotcn = -2;   % ellipticity
mx  = 48;
my  = 44;
a1  = 0.0008;  % [1/mm]
aob = 0.5;     % [-]
scale = 1.1;

params = [mx, my, a1, aob, scale];
cntc_sethertzcontact(iwhe, icp, ipotcn, params);

% rolling velocity

veloc = 10000; % [mm/s]

cntc_setreferencevelocity(iwhe, icp, veloc);

% total normal force specified, see N-digit above

fn = 106700;   % [N]

cntc_setnormalforce(iwhe, icp, fn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 3: perform loops for all the cases to be computed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% creepages: a list of values for creating a creep-force curve

cksi_m0  = [ 0.00001, 0.0004, 0.0008, 0.0012, 0.0016, 0.0020, 0.0024, 0.0028, ...
              0.0032, 0.0036, 0.0040, 0.0044, 0.0048, 0.0052, 0.0056, 0.0060, ...
              0.0080, 0.0100, 0.0200, 0.0300, 0.0400, 0.0500, 0.0750, 0.1000, ...
              0.1250, 0.1500, 0.1750, 0.2000, 0.2250, 0.2500 ]; % [-]
cksi_m3  = [ 0.00001, 0.0010, 0.0020, 0.0030, 0.0040, 0.0050, 0.0060, 0.0070, ...
              0.0080, 0.0100, 0.0120, 0.0140, 0.0160, 0.0180, 0.0200, 0.0220, ...
              0.0250, 0.0300, 0.0350, 0.0400, 0.0500, 0.0650, 0.0800, 0.1000, ...
              0.1250, 0.1500, 0.1750, 0.2000, 0.2250, 0.2500 ]; % [-]
cksi_m4  = [ 0.00001, 0.0010, 0.0020, 0.0030, 0.0040, 0.0050, 0.0060, 0.0070, ...
              0.0080, 0.0090, 0.0100, 0.0110, 0.0120, 0.0130, 0.0140, 0.0150, ...
              0.0160, 0.0180, 0.0200, 0.0300, 0.0400, 0.0500, 0.0750, 0.1000, ...
              0.1250, 0.1500, 0.1750, 0.2000, 0.2250, 0.2500 ]; % [-]
ceta = 0; % [-]
cphi = 0; % [rad/mm]

%  icurv = 1:   Original Kalker theory,     Coulomb friction
%          2:   Extended CONTACT,           exponential falling friction
%          3:   Original Fastsim algorithm, Coulomb friction
%          4:   Modified Fastsim algorithm, exponential falling friction

for icurv = 1 : 4

   % friction data

   if (icurv == 1 | icurv == 3)
      cksi = cksi_m0;
      ldigit = 0;      % L-digit 0: Coulomb friction
      fstat = 0.33;   % [-]
      fkin  = fstat;
      cntc_setfrictionmethod(iwhe, icp, ldigit, [fstat,fkin]);
   elseif (icurv == 2 | icurv == 4)
      ldigit = 4;      % L-digit 4: exponential falling friction
      if (icurv == 2)
         cksi = cksi_m4;
         fkin   = 0.14;   % [-]
         fexp1  = 0.19;
         sabsh1 = 1250;   % [mm/s]
      else
         cksi = cksi_m3;
         fkin   = 0.1368; % [-]
         fexp1  = 0.2232;
         sabsh1 =  990.0; % [mm/s]
      end
      memdst = 0.003; % [mm]
      mem_s0 = 1.000; % [mm/s]
      params = [fkin, fexp1, sabsh1, 0, 0, memdst, mem_s0];
      cntc_setfrictionmethod(iwhe, icp, ldigit, params);
   end

   % material data, interfacial layer

   gg = 82000;    % [N/mm^2]
   nu = 0.28;     % [-]

   if (icurv == 1)

      mdigit = 0;     % M-digit 0: linear elastic half-spaces, no interfacial layer
      params = [ nu, nu, gg, gg ];

   elseif (icurv == 2)

      mdigit = 4;     % M-digit 4: elastic interfacial layer, no plasticity
      gg3    = 8200;  % [N/mm2]
      h3     = 1.25;  % [mm]
      tau_c0 = 0;     % [N/mm2]
      k_tau  = 0;     % [N/mm2/mm]
      params = [ nu, nu, gg, gg, gg3, h3, tau_c0, k_tau];

   else

      mdigit  = 3;     % M-digit 3: modified Fastsim, automatic flexbilities L1,L2,L3
      if (icurv == 3)
         k0_mf  = 1.0;   % [-]
         alfamf = 1.0;   % [-]
         betamf = 1.0;   % [-]
      else
         k0_mf  = 0.54;  % [-]
         alfamf = 0.02;  % [-]
         betamf = 0.62;  % [-]
      end
      params = [ nu, nu, gg, gg, k0_mf, alfamf, betamf];
   end

   nparam = length(params);
   cntc_setmaterialparameters(iwhe, icp, mdigit, params);

   % Traction bound: for CONTACT or Fastsim

   clear flags values;
   if (icurv <= 2)
      flags(1) = CNTC.ic_bound ; values(1) = 0;    % B=0, elastic half-space
   else
      flags(1) = CNTC.ic_bound ; values(1) = 3;    % B=3, parabolic traction bound
   end
   cntc_setflags(iwhe, icp, flags, values);

   for iksi = 1:length(cksi)

      % set creepages according to next value from cksi

      cntc_setcreepages(iwhe, icp, cksi(iksi), ceta, cphi);

      % compute the contact problem

      ierror = cntc_calculate(iwhe, icp);

      % get forces on upper body (1) (CONTACT unit convention)

      [fn, fx, fy, mz] = cntc_getcontactforces(iwhe, icp);

      fx_list(iksi,icurv) = fx;

   end
end

% Cleanup

cntc_finalize(iwhe);
cntc_closelibrary;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 4: plotting of the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

make_fig = 1;
if (make_fig)
   figure(1); clf;
   set(1, 'defaultaxesfontsize',14);
   set(1, 'defaulttextfontsize',14);
   set(1, 'defaultlinelinewidth',2);
   set(1, 'defaultaxeslinewidth',2);
   hold on;

   clear sol;
   if (exist('tractcurv.out'))
      sol = parse_out3('tractcurv.out');
   elseif (exist('tractcurv.ref_out'))
      sol = parse_out3('tractcurv.ref_out');
   end
   if (exist('sol'))
      c2 = reshape(sol.creep.cksi, 30, 4);
      f2 = reshape(sol.force.fx, 30, 4);
      fstat = [0.33, 0.33, 0.33, 0.36];
      l1 = plot(100*c2(:,1), -f2(:,1)*fstat(1), '--','color','b');
      set(gca,'colororderindex', 3);
      l4 = plot(100*c2(:,4), -f2(:,4)*fstat(4), '-');
      set(gca,'colororderindex', 2);
      l2 = plot(100*c2(:,2), -f2(:,2)*fstat(2), '-');

      l5 = plot(100*cksi_m0, -fx_list(:,1)/fn, '.', 'color', [0 .7 .7], 'markersize', 12);
      l8 = plot(100*cksi_m3, -fx_list(:,4)/fn, '.', 'color', [0 .7 .7], 'markersize', 12);
      l6 = plot(100*cksi_m4, -fx_list(:,2)/fn, '.', 'color', [0 .7 .7], 'markersize', 12);
      legend([l1,l2,l4,l6],'Original Kalker theory','Extended CONTACT model', ...
                    'modified Fastsim algorithm', 'CONTACT library', ...
                    'location','SouthEast', 'autoupdate','off');
   else
      l5=plot(100*cksi_m0, -fx_list(:,1)/fn, 'b--');
      set(gca,'colororderindex', 3);
      l8=plot(100*cksi_m3, -fx_list(:,4)/fn, '-');
      set(gca,'colororderindex', 2);
      l6=plot(100*cksi_m4, -fx_list(:,2)/fn, '-*');
      legend([l5,l6,l8],'Original Kalker theory','Extended CONTACT model', ...
                    'modified Fastsim algorithm', 'location','SouthEast', 'autoupdate','off');
   end

   axis([0 5 0 0.4501]); grid on;
   xlabel('Longitudinal creepage \xi [%]')
   ylabel('Traction coefficient F_x/F_n [-]');

   text(0.2, 0.42, 'reduced initial slope');
   plot([0.4, 0.3], [0.405, 0.25], 'k-');
   plot([0.154 0.57], [0.25 0.25], 'k-');
   text(4, 0.35, 'falling friction effect','horizontalalignment','center');
   plot([4 4], [0.33, 0.297], 'k-');
   text(4.8, 0.42, 'User guide, Section 5.9','horizontalalignment','right');
end

% $Revision: 2593 $, $Date: 2024-08-14 15:18:54 +0200 (Wed, 14 Aug 2024) $
