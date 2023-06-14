
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Part 1: Initialize the CONTACT library, register problem "iwhe = 1"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Part 2: Configure the main flags & control digits of the contact problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 
            % SIMPACK unit convention: [m], [m/s], [N], acting on body (2)
flags(1) = CNTC.if_units ; values(1) = CNTC.un_spck;
flags(2) = CNTC.if_wrtinp; values(2) = 1;    % write .inp-file
flags(3) = CNTC.ic_matfil; values(3) = 0;    % write .mat-file
flags(4) = CNTC.ic_output; values(4) = 4;    % general output to .out-file
flags(5) = CNTC.ic_flow  ; values(5) = 4;    % lot of progress output
flags(6) = CNTC.ic_norm  ; values(6) = 1;    % N=1, FN prescribed
flags(7) = CNTC.ic_tang  ; values(7) = 1;    % T=3, steady state rolling
flags(8) = CNTC.ic_pvtime; values(8) = 2;    % P=2, no previous time

cntc_setflags(iwhe, icp, length(flags), flags, values);

% solver settings

cntc_setsolverflags(iwhe, icp, 0, 4, [999 100 30 1], 1, 1e-6);

% material data

gg = 82000e6;  % [N/m^2]
nu = 0.28;     % [-]

cntc_setmaterialproperties(iwhe, icp, gg, nu, gg, nu);

% Hertzian discretization, one element of 1 mm2

ipotcn = -3;   % semi-axes
mx  = 1;
my  = 1;
aa  = 0.0005;  % [m]
bb  = 0.0005;  % [m]
scale = 1.0;

params = [mx, my, aa, bb, scale];
cntc_sethertzcontact(iwhe, icp, ipotcn, 5, params);

% friction data

imeth = 0;     % L-digit 0: Coulomb friction
fstat = 0.70;  % [-]
fkin  = fstat;
cntc_setfrictionmethod(iwhe, icp, imeth, 2, [fstat,fkin]);

% total normal force specified, see N-digit above

fn  = 900;     % [N]
cntc_setnormalforce(iwhe, icp, fn);

% compute series of steps for four interfacial layers

nser  = 4;
ncase = [13, 17, 13, 13];
gg3   = 1e6  * [ 176.63, 106.25,  176.63, 20.238 ]; % [N/m2]
tauc0 = 1e6  * [ 560.00, 200.00,  400.00, 20.000 ]; % [N/m2]
ktau  = 1e9  * [   0.00, 435.00, -195.00, 20.400 ]; % [N/m3]
cksi  = 1e-3 * [ -0.035 -0.035 -0.030 -0.050*ones(1,10), ...
                 -0.040 -0.040 -0.040 -0.040 -0.036 0.000 0.036 0.000 -0.046 -0.044 -0.050*ones(1,7), ...
                 -0.025 -0.025 -0.050*ones(1,11), ...
                 -0.020 -0.040 -0.040 -0.050*ones(1,10) ]; % [m]
         
iofs  = 0;
shft  =  ones(20,nser) * NaN;
px    = zeros(20,nser);
fx    = zeros(20,nser);
ux    = zeros(20,nser);
uplx  = zeros(20,nser);

for iser = 1 : nser

   % interfacial layer

   imeth  = 4;         % M-digit 4: elasto-plastic interfacial layer
   h3     = 0.000020;  % [m]
   params = [gg3(iser), h3, tauc0(iser), ktau(iser)];
   nparam = length(params);
   cntc_setinterfaciallayer(iwhe, icp, imeth, nparam, params);

   % output data

   shft(1,iser) = 0;

   icase = 0;
   while(icase < ncase(iser))

      icase = icase + 1;

      % set shift for current series

      if (icase<=1)
         cntc_setflags(iwhe, icp, 1, CNTC.ic_pvtime, 2);  % P=2, no previous time
      else
         cntc_setflags(iwhe, icp, 1, CNTC.ic_pvtime, 0);  % P=0, full sequence
      end

      cntc_setcreepages(iwhe, icp, cksi(iofs+icase), 0, 0);

      % compute the contact problem

      ierror = cntc_calculate(iwhe, icp);

      % get outputs

      shft(1+icase,iser) = shft(icase,iser) - cksi(iofs+icase);
      
      [fni, fxi] = cntc_getcontactforces(iwhe, icp);
      fx(1+icase,iser)   = fxi;

      fld_px = cntc_getfielddata(iwhe, icp, CNTC.fld_px, mx, my);
      px(1+icase,iser)   = fld_px(1,1);
      fld_ux = cntc_getfielddata(iwhe, icp, CNTC.fld_ux, mx, my);
      ux(1+icase,iser)   = fld_ux(1,1);
      fld_uplx = cntc_getfielddata(iwhe, icp, CNTC.fld_uplsx, mx, my);
      uplx(1+icase,iser) = fld_uplx(1,1);
   end

   iofs = iofs + ncase(iser);
end


% Cleanup

cntc_finalize(iwhe);
cntc_closelibrary;

clear flags values gg nu ipotcn mx my aa bb scale params imeth fstat fkin fn cksi fni fxi
clear nparam imeth gg3 h3 ktau tauc0 fld_px fld_ux fld_uplx
clear CNTC ifcver idebug ierror imodul iwhe icp iser icase iofs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Part 4: plotting of the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

make_fig = 1;
if (make_fig)

   figure(1); clf; hold on;
   plot(1e6*shft, -px/1e6, '.-', 'markersize',12);
   plot(1e6*shft, -fx, ':');
   xlabel('applied displacement [\mu m]');
   ylabel('shear stress \tau [N/mm^2]');
   axis([0 600 0 620]); grid on
   legend('"magnetite"','"clay"','"sand"','"MoS2"', 'location', 'southeast');
   text(110, 255, '\rightarrow', 'rotation', 18);
   text(163, 155, '\rightarrow', 'rotation', -108);
   text(188, 140, '\rightarrow', 'rotation', 72);
   text(260, 315, '\rightarrow', 'rotation', 18);
   plot(70, 560, 'k.', 'markersize', 12);
   text(75, 520, '$(u_{c0}, \tau_{c0})$', 'interpreter', 'latex');
   plot([460 535 535], [368 368 398], 'k-', 'linewidth', 1);
   text(545, 360, '$k_u$', 'interpreter', 'latex');

   figure(2); clf; hold on;
   l = plot(shft, ux+uplx);
   set(l(4), 'linestyle',':');
   set(l(3), 'linestyle','--');
   xlabel('total shift [m]');
   ylabel('total displacement [m]');
   axis([0 6e-4 0 6e-4]); grid on;

   clear make_fig l;
end

