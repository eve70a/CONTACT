
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example program using CONTACT library for Matlab. See User guide, Sections 7.5 & 5.6.
% Calculation of subsurface stresses, cf. [Kalker1990], Figure 5.20.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1: Initialize the CONTACT library, register problem "iwhe = 1"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~exist('cntc_initlibrary.m','file'))
   % set location of CONTACT installation folder
   % contactdir = 'C:\Program Files\Vtech CMCC\contact_v23.1';
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2: Configure the main flags & control digits of the contact problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
            % CONTACT unit convention: [mm], [mm/s], [N], acting on body (1)
flags( 1) = CNTC.if_units ; values( 1) = CNTC.un_cntc;
flags( 2) = CNTC.if_wrtinp; values( 2) = 0;    % no .inp-file needed
flags( 3) = CNTC.ic_matfil; values( 3) = 0;    % no .mat-file needed
flags( 4) = CNTC.ic_output; values( 4) = 1;    % min. output to .out-file
flags( 5) = CNTC.ic_flow  ; values( 5) = 1;    % a little of progress output
flags( 6) = CNTC.ic_norm  ; values( 6) = 1;    % N=1, FN prescribed
flags( 7) = CNTC.ic_tang  ; values( 7) = 0;    % T=0, frictionless
flags( 8) = CNTC.ic_pvtime; values( 8) = 2;    % P=2, no previous time
flags( 9) = CNTC.ic_sbsout; values( 9) = 2;    % O_s=2, subs-results to out-file
flags(10) = CNTC.ic_sbsfil; values(10) = 1;    % A_s=1, write subs-file

cntc_setflags(iwhe, icp, flags, values);

% material data

gg = 1.00;     % [N/mm^2]
nu = 0.28;     % [-]

cntc_setmaterialproperties(iwhe, icp, gg, nu, gg, nu);

% Hertzian discretization, one element of 1 x 1 mm

ipotcn = -3;   % semi-axes aa, bb
mx  = 1;
my  = 1;
aa  = 0.5;     % [mm]
bb  = 0.5;     % [mm]
scale = 1.0;

params = [mx, my, aa, bb, scale];
cntc_sethertzcontact(iwhe, icp, ipotcn, params);

% rolling velocity

veloc = 10000; % [mm/s]

cntc_setreferencevelocity(iwhe, icp, veloc);

% total normal force specified, see N-digit above

fn = 1;   % [N]

cntc_setnormalforce(iwhe, icp, fn);

% coefficient of friction

imeth = 0;      % L-digit 0: Coulomb friction
fstat = 1.0;    % [-]
fkin  = fstat;
cntc_setfrictionmethod(iwhe, icp, imeth, [fstat,fkin]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set positions for evaluation of subsurface stresses

iblk  = 1;
isubs = 9;
xarr  = [ 0 ];
yarr  = [ 0 ];
zarr  = [ 0    0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9 ...
          1.0  1.25  1.667 2.5   5.0 ];

subs_addblock(iwhe, icp, iblk, isubs, xarr, yarr, zarr);

iblk  = 2;
xarr  = [ -5.0  -3.0  -2.0  -1.0  -0.8  -0.6  -0.4  -0.3  -0.2  -0.1   0.0 ...
           0.1   0.2   0.3   0.4   0.6   0.8   1.0   2.0   3.0   5.0 ];
yarr  = xarr;

subs_addblock(iwhe, icp, iblk, isubs, xarr, yarr, zarr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the contact problem for case 1: frictionless compression

ierror = cntc_calculate(iwhe, icp);

% compute subsurface stresses for case 1

ierror = subs_calculate(iwhe, icp);

% for block 1, get selected results of subsurface stress calculation

icol = [3, 7, 8]; % ['z', 'sighyd', 'sigvm']
table1a = subs_getresults(iwhe, icp, 1, icol);

% for block 2, get struct with all results of subsurface stress calculation

blk1b   = subs_getresults(iwhe, icp, 2, 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the contact problem for case 2: frictional rolling

clear flags values;
flags(1) = CNTC.ic_tang  ; values(1) = 3;    % T=3, steady state rolling
flags(2) = CNTC.ic_force ; values(2) = 1;    % F=1, FX prescribed
cntc_setflags(iwhe, icp, flags, values);

cntc_settangentialforces(iwhe, icp, 0.999, 0);

ierror = cntc_calculate(iwhe, icp);

% compute subsurface stresses for case 2

ierror = subs_calculate(iwhe, icp);

% for block 1, get selected results of subsurface stress calculation

icol = [3, 7, 8]; % ['z', 'sighyd', 'sigvm']
table2a = subs_getresults(iwhe, icp, 1, icol);

% for block 2, get struct with all results of subsurface stress calculation

blk2b   = subs_getresults(iwhe, icp, 2, 'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create pictures cf. example 'subsurf.inp', User guide Sec. 6.9

% plot hydrostatic & von Mises stress at (x,y)=(0,0)

sighyd1 = table1a(:,2);
sigvm1  = table1a(:,3);
sigvm2  = table2a(:,3);
difvm   = sqrt(sigvm2.^2 - sigvm1.^2);

figure(1); clf; hold on;
plot(table1a(:,1), -sighyd1, '-o');
plot(table1a(:,1),  sigvm1, '-*');
plot(table1a(:,1),  difvm, '--*');
axis([0 2.5 0 1.8]); grid; set(gca,'ytick',[0:0.2:1.8]);
set(gca,'dataaspectratio',[1.5 1 1]);
xlabel('z [mm]');
ylabel('stress [N/mm^2]');
legend('-\sigma_{hyd} for p_z=1, p_x=p_y=0', ...
        '\sigma_{vm} for p_z=1, p_x=p_y=0', ...
        '\sigma_{vm} for p_z=0, p_x=1, p_y=0');
text(0.95, 0.60, 'single element, dx=dy=1', ...
        'horizontalalignment','right', 'units','normalized');
text(0.95, 0.52, 'G=1, \nu=0.28', ...
        'horizontalalignment','right', 'units','normalized');

% plot displacements uz along Oxz-plane, y=0

figure(2); clf; hold on;
opt = plotstrs;
opt.field   = 'uz';
opt.yslc    = 0;
opt.typplot = 'contourf';
opt.cntrlvl = [0:0.02:0.08, 0.12:0.04:0.40];
plotstrs(blk2b, opt);
% grid on;
set(gca,'clim',[0 0.40]);
h = findobj(gcf,'Type','colorbar');
set(h,'ylim',[0 0.40], 'ytick',opt.cntrlvl);
title('');

% Cleanup

cntc_finalize(iwhe);
cntc_closelibrary;

% $Revision: 2372 $, $Date: 2023-06-15 16:58:45 +0200 (Thu, 15 Jun 2023) $
