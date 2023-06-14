
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1: Initialize the CONTACT library, register problem "iwhe = 1"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~exist('cntc_initlibrary.m','file'))
   % set location of CONTACT installation folder
   % contact_dir = 'c:\program files (x86)\KalkerSoftware\contact_v19.1';
   contact_dir = '..';
   addpath([contact_dir, '\matlab_intfc']);
   addpath([contact_dir, '\matlab']);
end

if (exist('cntc_initlibrary')~=2)
   disp('ERROR: cant find CONTACT library on the Matlab search path');
   return
end
[CNTC, ifcver, ierror] = cntc_initlibrary;

if (ierror<0)
   disp(sprintf('An error occurred, ierror = %d',ierror));
   disp(sprintf('(ierror=-12 means that no valid license is found.)'));
   return
end

idebug = 1;             % just a bit of information from the library
cntc_setglobalflags(1, CNTC.if_idebug, idebug);

iwhe   = 1; % "wheel number"
icp    = 1; % "contact problem on wheel"
imodul = 3; % basic contact
[ifcver, ierror] = cntc_initialize(iwhe, imodul);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2: Configure the main flags & control digits of the contact problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
            % CONTACT unit convention: [mm], [mm/s], [N], acting on body (1)
flags(1) = CNTC.if_units ; values(1) = CNTC.un_cntc;
flags(2) = CNTC.if_wrtinp; values(2) = 1;    % write .inp-file
flags(3) = CNTC.ic_matfil; values(3) = 2;    % write .mat-file
flags(4) = CNTC.ic_output; values(4) = 4;    % general output to .out-file
flags(5) = CNTC.ic_flow  ; values(5) = 4;    % lot of progress output
flags(6) = CNTC.ic_norm  ; values(6) = 1;    % N=1, FN prescribed
flags(7) = CNTC.ic_tang  ; values(7) = 3;    % T=3, steady state rolling
flags(8) = CNTC.ic_pvtime; values(8) = 2;    % P=2, no previous time

cntc_setflags(iwhe, icp, length(flags), flags, values);

% solver settings

cntc_setsolverflags(iwhe, icp, 0, 4, [999 100 30 1], 1, 1e-5);

% material data

gg = 82000;    % [N/mm^2]   (82000e6 N/m^2)
nu = 0.28;     % [-]

cntc_setmaterialproperties(iwhe, icp, gg, nu, gg, nu);

% temperature inputs

h_digit =   3;
bktemp =    0;     % [C]
heatcp =  450;     % [J/kg-C]
lambda =   50e-3;  % [W/mm-C]   (50 W/m-C)
dens   = 7.85e-6;  % [kg/mm^3]  (7.85e3 kg/m3)
params = [ bktemp, heatcp, lambda, dens, bktemp, heatcp, lambda, dens ];

cntc_settemperaturedata(iwhe, icp, h_digit, 8, params);

% non-Hertzian discretization, extended potential contact area

ipotcn = 1;    % lower-left + step sizes
mx  = 216;
my  =  45;
xl  =  -7.90;  % [mm]  (-0.0079 m)
yl  = -11.25;  % [mm]  (-0.01125 m)
dx  =   0.20;  % [mm]  ( 0.0002 m)
dy  =   0.50;  % [mm]  ( 0.0005 m)

params = [mx, my, xl, yl, dx, dy];
cntc_setpotcontact(iwhe, icp, ipotcn, 6, params);

% quadratic u.d.

b1 = 0.001015;  % [1/mm]  (1.015 /mm)
b3 = 0.0004215; % [1/mm]  (0.4215 /mm)
ibase = 1;
cntc_setundeformeddistc(iwhe, icp, ibase, 6, [b1,0,b3,0,0,0]);

% rolling velocity

veloc = 30000;  % [mm/s]  (30 m/s)

cntc_setreferencevelocity(iwhe, icp, veloc);

% rolling direction & stepsize

chi = pi;
dq  = dx;
cntc_setrollingstepsize(iwhe, icp, chi, dq);

% total normal force specified, see N-digit above

fn = 100000;   % [N]

cntc_setnormalforce(iwhe, icp, fn);

% creepages

cksi = -0.03333;
ceta =  0;
cphi =  0;     % [rad/mm]  (0 rad/m)

cntc_setcreepages(iwhe, icp, cksi, ceta, cphi);

% Calculate three cases: basic with temperature, friction reducing, friction increasing

tcentr = []; mu = []; px = []; srel = [];

for icase = 1 : 3

   % friction data

   if (icase==1)

      imeth = 0;      % L-digit 0: Coulomb friction
      fstat = 0.30;   % [-]
      fkin  = fstat;
      cntc_setfrictionmethod(iwhe, icp, imeth, 2, [fstat,fkin]);

   else

      imeth = 6;      % L-digit 6: Temperature-dependent friction
      fref  = 0.30;   % [-]
      tref  = 0;      % [C]
      if (icase==2)
         dfheat = -0.18; % [-]
      else
         dfheat =  0.18; % [-]
      end
      dtheat = 400;   % [C]
      memdst = 0;     % [mm]    (0 m)
      mem_s0 = 1;     % [mm/s]  (0.001 m/s)
      cntc_setfrictionmethod(iwhe, icp, imeth, 6, [fref, tref, dfheat, dtheat, memdst, mem_s0]);

   end

   % compute the contact problem

   ierror = cntc_calculate(iwhe, icp);

   % get maximum temperatures on both bodies

   [t1, t2] = cntc_getmaximumtemperature(iwhe, icp);
   disp(sprintf('Case %d: max. temperatures t1=%6.1f, t2=%6.1f', icase, t1, t2));

   % get grid for making pictures, select y=0

   [mx, my, xcl, ycl, dx, dy] = cntc_getpotcontact(iwhe, icp);
   s.x = xcl + [0:mx-1] * dx;
   s.y = ycl + [0:my-1] * dy;
   [~,iy0] = min(abs(s.y));

   % get fields for making pictures

   s.temp1 = cntc_getfielddata(iwhe, icp, CNTC.fld_temp1, mx, my);
   s.mu    = cntc_getfielddata(iwhe, icp, CNTC.fld_mu, mx, my);
   s.px    = cntc_getfielddata(iwhe, icp, CNTC.fld_px, mx, my);
   s.sx    = cntc_getfielddata(iwhe, icp, CNTC.fld_sx, mx, my);
   s.sy    = cntc_getfielddata(iwhe, icp, CNTC.fld_sy, mx, my);
   s.srel  = sqrt(s.sx.^2 + s.sy.^2);
   tcentr = [tcentr; s.temp1(iy0,:)];
   mu     = [mu;     s.mu(iy0,:)];
   px     = [px;     s.px(iy0,:)];
   srel   = [srel;   s.srel(iy0,:)];

   % figure(icase)
   % surf(s.x, s.y, s.temp1); view([0 90]); shading flat; colorbar;
end

% Cleanup

cntc_finalize(iwhe);
cntc_closelibrary;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 4: plotting of the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

make_fig = 1;
if (make_fig)

   a = 6;  % [mm]  (0.006 m)
   figure(6); clf;
   p = get(gcf,'position'); p(3:4) = [800 600]; set(gcf,'position',p);
   subplot(2,2,1); hold on;
   plot(s.x/a, px, '-');
   ax = axis; ax(1:2) = [-1.2 6]; axis(ax); grid on
   legend('Constant \mu', '\Delta\mu_{heat}=-0.18', '\Delta\mu_{heat}=+0.18');
   % legend('reference','decreasing \mu', 'increasing \mu');
   title('traction');
   ylabel('p_x [N/mm^2]');

   subplot(2,2,2);
   plot(s.x/a, srel, '-');
   axis([-1.2 6 0.025 0.040]); grid on
   title('relative slip');
   ylabel('s [-]');

   subplot(2,2,3); hold on;
   plot(s.x/a, tcentr, '-');
   plot([-1 -1], [175 225], 'k--');
   plot([ 1  1], [175 225], 'k--');
   axis([-1.2 6 0 250]); grid on
   title('wheel temperature');
   xlabel('x/a [-]'); ylabel('T_a [{}^\circ{}C]');

   subplot(2,2,4); hold on
   plot(s.x/a, mu, '-');
   plot([-1 -1], [0.375 0.425], 'k--');
   plot([ 1  1], [0.375 0.425], 'k--');
   ax = axis; ax(1:4) = [-1.2 6 0.2 0.42]; axis(ax); grid on
   title('friction coefficient');
   xlabel('x/a [-]'); ylabel('\mu [-]');

end

