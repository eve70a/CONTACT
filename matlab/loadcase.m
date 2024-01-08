
function [ s1, s2, s3, s4, s5, s6, s7, s8, s9, s10 ] = loadcase( expnam, icase, ipatch )
%
% [ s1, s2, ... ] = loadcase( expnam, [icase], [ipatch] )
%
% Read CONTACT mat-file(s) for one case and return struct(s) with its(their) contents.
%
% expnam    = the experiment name used, i.e. the filename without .mat extension
% icase     = the case-number to be read
% ipatch    = optional contact patch number, <=0 means all available
%               0 = all available, returned as an array of structures [s1]
%              -1 = all available, returned using separate structures  s1, s2, ...
%             default: 0 when nargin<=1, -1 when nargin>1
% s1,...    = structure(s) with results of the computation

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

max_patch = 10;
s1 = []; s2 = []; s3 = []; s4 = []; s5 = [];
s6 = []; s7 = []; s8 = []; s9 = []; s10 = [];

if (nargin<2 | isempty(icase))
   icase = 1;
end
if (nargin<3 | isempty(ipatch))
   if (nargout<=1)
      ipatch = 0;
   else
      ipatch = -1;
   end
end

if (ipatch<0)
   use_array = 0;
else
   use_array = 1;
end

% check for the presence of output-files "experiment.<icase>[a-z].mat", ...

has_output = zeros(max_patch,1);
for jpatch = 1:max_patch
   if (icase<=9999)
      fname = [expnam, sprintf('.%04d%c.mat', icase, 96+jpatch)];
   else
      fname = [expnam, sprintf('.%06d%c.mat', icase, 96+jpatch)];
   end
   if (exist(fname,'file'))
      has_output(jpatch) = 1;
   end
end

% set default number of patches: all available

if (max(has_output)<=0)
   npatch = 1;
else
   npatch = max(find(has_output));
end

% set requested range of patches

if (ipatch<=0)
   ipatch = [1 : npatch];
end

% loop over the requested patches

for jpatch = ipatch

   % npatch may be reduced to value obtained from actual mat-files

   if (jpatch > npatch)
      continue
   end

   tmp = [];
   if (jpatch>1 & ~has_output(jpatch))

      % short-cut for [b-z]: no data found

      disp(sprintf('No output found for contact patch %d.',jpatch));
      continue

   elseif (has_output(jpatch))

      % load file "<icase>[a-z]" if present

      if (icase<=9999)
         fname = [expnam, sprintf('.%04d%c.mat', icase, 96+jpatch)];
      else
         fname = [expnam, sprintf('.%06d%c.mat', icase, 96+jpatch)];
      end
      tmp = load(fname, '-ascii');

   else

      % else, check for presence of file "experiment.0001.mat" or "experiment.000001.mat"

      if (icase<=9999)
         fname = [expnam, sprintf('.%04d.mat', icase)];
      else
         fname = [expnam, sprintf('.%06d.mat', icase)];
      end
      fname_ic = fname;

      if (exist(fname,'file'))

         tmp = load(fname, '-ascii');

      else

         % else, check for presence of file "experiment.mat"

         fname = [expnam,'.mat'];
         if (exist(fname,'file'))
            tmp = load(fname, '-ascii');
         else
            disp(['ERROR: cannot find file ',fname,' nor ',fname_ic]);
         end
      end
   end

   % create solution struct

   if (isempty(tmp))
      strc = [];
   else
      strc = fill_struct(tmp, fname);
   end

   % Each patch tells the actual npatch for the run in which it was created
   % warn & update npatch as needed

   if (~isempty(strc) & strc.meta.npatch>0 & strc.meta.npatch < npatch)
      disp(sprintf('Warning: found %d mat-files for run with %d contact patches. Left-over from earlier run?', ...
                npatch, strc.meta.npatch));
      npatch = min(npatch, strc.meta.npatch);
   end

   % put solution struct into output argument

   eval(sprintf('s%d = strc;', jpatch));

end

if (use_array)
   s1 = [s1, s2, s3, s4, s5, s6, s7, s8, s9, s10];
end

end % function loadcase

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [ sol ] = fill_struct(tmp, fname)

% check the file format

fmtmat = tmp(1,end);
if (fmtmat>2320)
   disp('ERROR: the mat-file is created with a newer version of CONTACT.');
   disp('       Please use the corresponding version of this Matlab script.');
   sol=[];
   return
end

numcol = size(tmp,2);
irow = 1;

sol = struct('config',[], 'h_digit',[], 'meta',struct(), 'mater',struct(), 'fric',struct(), ...
             'kincns',struct('t_digit',[]) );

% extract the metadata from the first three rows

if (fmtmat>=2320)
   sol.meta.tim      = tmp(irow,1);     % time: may be provided by calling program
   sol.meta.s_ws     = tmp(irow,2);     % s-position along track curve, wrt. inertial reference
   sol.meta.th_ws    = tmp(irow,3);
   sol.meta.y_r      = tmp(irow,4);     % position of rail marker wrt track reference
   sol.meta.z_r      = tmp(irow,5);
   sol.meta.roll_r   = tmp(irow,6);     % orientation of rail marker wrt track reference
   sol.meta.rnom_rol = tmp(irow,7);
   irow = irow + 1;

   sol.meta.x_w      = tmp(irow,1);     % position of wheel marker wrt track reference
   sol.meta.y_w      = tmp(irow,2);
   sol.meta.z_w      = tmp(irow,3);
   sol.meta.roll_w   = tmp(irow,4);     % orientation of wheel marker wrt track reference
   sol.meta.yaw_w    = tmp(irow,5);
   sol.meta.rnom_whl = tmp(irow,6);
   irow = irow + 1;

   sol.meta.xcp_r    = tmp(irow,1);     % position of contact ref.point wrt rail marker
   sol.meta.ycp_r    = tmp(irow,2);
   sol.meta.zcp_r    = tmp(irow,3);
   sol.meta.deltcp_r = tmp(irow,4);     % orientation of contact ref.point wrt rail marker
   sol.meta.xcp_w    = tmp(irow,5);     % position of contact ref.point wrt wheel marker
   sol.meta.ycp_w    = tmp(irow,6);
   sol.meta.zcp_w    = tmp(irow,7);
   sol.meta.npatch   = tmp(irow,8);
   sol.meta.ipatch   = tmp(irow,9);
   sol.meta.spinxo   = tmp(irow,10);
   sol.meta.spinyo   = tmp(irow,11);
   irow = irow + 1;
elseif (fmtmat>=1801)
   sol.meta.tim      = tmp(irow,1);     % time: may be provided by calling program
   sol.meta.s_ws     = tmp(irow,2);     % s-position along track curve, wrt. inertial reference
   sol.meta.x_w      = tmp(irow,3);     % position of wheel marker wrt track reference
   sol.meta.y_w      = tmp(irow,4);
   sol.meta.z_w      = tmp(irow,5);
   sol.meta.roll_w   = tmp(irow,6);     % orientation of wheel marker wrt track reference
   sol.meta.yaw_w    = tmp(irow,7);
   sol.meta.pitch_w  = 0;
   sol.meta.y_r      = tmp(irow,8);     % position of rail marker wrt track reference
   sol.meta.z_r      = tmp(irow,9);
   sol.meta.roll_r   = tmp(irow,10);    % orientation of rail marker wrt track reference
   irow = irow + 1;

   sol.meta.xcp_r    = tmp(irow,1);     % position of contact ref.point wrt rail marker
   sol.meta.ycp_r    = tmp(irow,2);
   sol.meta.zcp_r    = tmp(irow,3);
   sol.meta.deltcp_r = tmp(irow,4);     % orientation of contact ref.point wrt rail marker
   sol.meta.xcp_w    = tmp(irow,5);     % position of contact ref.point wrt wheel marker
   sol.meta.ycp_w    = tmp(irow,6);
   sol.meta.zcp_w    = tmp(irow,7);
   sol.meta.npatch   = tmp(irow,8);
   sol.meta.ipatch   = tmp(irow,9);
   sol.meta.rnom_whl = tmp(irow,10);
   sol.meta.rnom_rol = tmp(irow,11);
   irow = irow + 1;
else
   sol.meta.tim      = 0;
   sol.meta.s_ws     = 0;
   sol.meta.x_w      = 0;
   sol.meta.y_w      = 0;
   sol.meta.z_w      = 0;
   sol.meta.roll_w   = 0;
   sol.meta.yaw_w    = 0;
   sol.meta.pitch_w  = 0;
   sol.meta.y_r      = 0;
   sol.meta.z_r      = 0;
   sol.meta.roll_r   = 0;
   sol.meta.xcp_r    = 0;
   sol.meta.ycp_r    = 0;
   sol.meta.zcp_r    = 0;
   sol.meta.deltcp_r = 0;
   sol.meta.xcp_w    = 0;
   sol.meta.ycp_w    = 0;
   sol.meta.zcp_w    = 0;
   sol.meta.npatch   = 0;
   sol.meta.ipatch   = 0;
   sol.meta.rnom_whl = 460;
   sol.meta.rnom_rol = 800;
   if (fmtmat>=1401)
      sol.meta.tim      = tmp(irow,1);
      sol.meta.xcp_r    = tmp(irow,2);
      sol.meta.ycp_r    = tmp(irow,3);
      irow = irow + 1;
   end
end

% extract the grid discretisation parameters from the next row

sol.mx  = tmp(irow, 1); mx = sol.mx;
sol.my  = tmp(irow, 2); my = sol.my;
sol.xl  = tmp(irow, 3);
sol.yl  = tmp(irow, 4);
sol.dx  = tmp(irow, 5);
sol.dy  = tmp(irow, 6);
sol.kincns.chi = tmp(irow, 7);
sol.kincns.dq  = tmp(irow, 8);
if (fmtmat>=1901)
   sol.config  = tmp(irow, 9);
   sol.d_digit = tmp(irow, 10);
elseif (sol.meta.yrw<0)
   sol.config = 0; sol.d_digit = 0;
else
   sol.config = 1; sol.d_digit = 0;
end
if (fmtmat>=2310)
   sol.meta.ynom_whl = tmp(irow,11);
else
   sol.meta.ynom_whl = 0;
end

% x_offset, y_offset: shift, real world coordinates to be used for the
%              origin of the contact coordinate system

sol.x_offset = [];
sol.y_offset = [];

% construct the grid coordinates

sol.x = sol.xl + ([1:mx]-0.5)*sol.dx;
sol.y = sol.yl + ([1:my]-0.5)*sol.dy;

% extract the material constants from the next row

if (fmtmat<1220)
   disp('WARNING: the mat-file is created with an old version of CONTACT.');
   disp('         The material constants are not filled in.');
   sol.mater = struct('m_digit',0,'gg',[1 1],'poiss',[1,1]);
   sol.kincns.t_digit = 3;
else
   irow = irow + 1;
   sol.kincns.t_digit = tmp(irow,1);
   sol.mater.m_digit  = tmp(irow,2);
   sol.mater.gg       = tmp(irow,3:4);
   sol.mater.poiss    = tmp(irow,5:6); % note: poiss(2) n.a. when m_digit==2 or 3
   if (sol.mater.m_digit==1)
      sol.mater.fg      = tmp(irow,7:8);
      sol.mater.tc      = tmp(irow,9:10);
   elseif (sol.mater.m_digit==2 | sol.mater.m_digit==3)
      sol.mater.poiss(2)= sol.mater.poiss(1); % note: poiss(2) n.a. when m_digit==2 or 3
      sol.mater.flx     = tmp(irow,6:8);
      sol.mater.k0_mf   = tmp(irow,9);
      sol.mater.alfamf  = tmp(irow,10);
      sol.mater.betamf  = tmp(irow,11);
   elseif (sol.mater.m_digit==4)
      sol.mater.gg3     = tmp(irow,7);
      sol.mater.laythk  = tmp(irow,8);
      sol.mater.tau_c0  = tmp(irow,9);
      sol.mater.k_tau   = tmp(irow,10);
   end

   g1 = sol.mater.gg(1); nu1 = sol.mater.poiss(1);
   g2 = sol.mater.gg(2); nu2 = sol.mater.poiss(2);
   sol.mater.ga = 2.0 / (1.0/g1 + 1.0/g2);
   sol.mater.si = 0.5 * sol.mater.ga * (nu1/g1 + nu2/g2);
   sol.mater.ak = 0.25 * sol.mater.ga * ((1-2*nu1)/g1 - (1-2*nu2)/g2);

end
if (sol.mater.m_digit~=4)
   has_plast = 0;
else
   has_plast = (sol.mater.tau_c0>1e-10 & sol.mater.tau_c0<1e10);
end

% extract the friction law parameters from the next row

irow = irow + 1;
frclaw = tmp(irow, 1); sol.fric.frclaw = frclaw;
if (frclaw<0 | frclaw>6)
   error(sprintf('The value for the L-digit (%d) in %s is incorrect. Maybe the file is generated with an older CONTACT version?',frclaw,fname));
   return
end
if (fmtmat<1220)
   idum=0;
   sol.h_digit = 0;
else
   idum=1; % introduced dummy in 2nd column in version 12.2
   sol.h_digit = tmp(irow,2);
end
has_heat   = (sol.h_digit>=1);
sol.kincns.veloc  = tmp(irow, 2+idum);
if (frclaw == 0)
   fstat  = tmp(irow, 3+idum); sol.fric.fstat  = fstat;
   fkin   = tmp(irow, 4+idum); sol.fric.fkin   = fkin;
elseif (frclaw == 2)
   fkin   = tmp(irow, 3+idum); sol.fric.fkin   = fkin;
   flin1  = tmp(irow, 4+idum); sol.fric.flin1  = flin1;
   sabsh1 = tmp(irow, 5+idum); sol.fric.sabsh1 = sabsh1;
   flin2  = tmp(irow, 6+idum); sol.fric.flin1  = flin1;
   sabsh2 = tmp(irow, 7+idum); sol.fric.sabsh2 = sabsh2;
   fstat  = fkin + flin1 + flin2; sol.fric.fstat = fstat;
elseif (frclaw == 3)
   fkin   = tmp(irow, 3+idum); sol.fric.fkin   = fkin;
   frat1  = tmp(irow, 4+idum); sol.fric.frat1  = frat1;
   sabsh1 = tmp(irow, 5+idum); sol.fric.sabsh1 = sabsh1;
   frat2  = tmp(irow, 6+idum); sol.fric.frat1  = frat1;
   sabsh2 = tmp(irow, 7+idum); sol.fric.sabsh2 = sabsh2;
   fstat  = fkin + frat1 + frat2; sol.fric.fstat = fstat;
elseif (frclaw == 4)
   fkin   = tmp(irow, 3+idum); sol.fric.fkin   = fkin;
   fexp1  = tmp(irow, 4+idum); sol.fric.fexp1  = fexp1;
   sabsh1 = tmp(irow, 5+idum); sol.fric.sabsh1 = sabsh1;
   fexp2  = tmp(irow, 6+idum); sol.fric.fexp2  = fexp2;
   sabsh2 = tmp(irow, 7+idum); sol.fric.sabsh2 = sabsh2;
   fstat  = fkin + fexp1 + fexp2; sol.fric.fstat = fstat;
elseif (frclaw == 6)
   fref   = tmp(irow, 3+idum); sol.fric.fref   = fref;
   tref   = tmp(irow, 4+idum); sol.fric.tref   = tref;
   dfheat = tmp(irow, 5+idum); sol.fric.dfheat = dfheat;
   dtheat = tmp(irow, 6+idum); sol.fric.dtheat = dtheat;
   fstat  = fref;              sol.fric.fstat  = fstat;
   fkin   = fref+dfheat;       sol.fric.fkin   = fkin;
end
if (frclaw >=2 & frclaw <=4)
   sol.fric.memdst = tmp(irow, 8+idum);
   sol.fric.mem_s0 = tmp(irow, 9+idum);
elseif (frclaw==6)
   sol.fric.memdst = tmp(irow, 7+idum);
   sol.fric.mem_s0 = tmp(irow, 8+idum);
end

% handle case of empty contact area

exter_el = 0; 
adhes_el = 1; 
slip_el  = 2;
plast_el = 3;

irow = irow + 1;
iend = size(tmp,1);
ncon = iend - irow + 1;

if (ncon<=0)
   disp('ERROR: mat-file contains no elements in the contact area');
   sol.eldiv  =  ones(my,mx)*exter_el;
   sol.h      = zeros(my,mx);
   sol.mu     = zeros(my,mx);
   sol.pn     = zeros(my,mx);
   sol.px     = zeros(my,mx);
   sol.py     = zeros(my,mx);
   sol.un     = zeros(my,mx);
   sol.ux     = zeros(my,mx);
   sol.uy     = zeros(my,mx);
   sol.srel   = zeros(my,mx);
   sol.shft   = zeros(my,mx);
   sol.trcbnd = zeros(my,mx);
   return
end

% get the numbers of the elements that are inside the contact area

ir = tmp(irow:end,1);

% extract the data for all active elements from the table

sol.eldiv =  ones(mx*my,1)*exter_el;
sol.h     =  ones(mx*my,1)*2*max(tmp(irow:iend,3));
sol.mu    = zeros(mx*my,1);
sol.pn    = zeros(mx*my,1);
sol.px    = zeros(mx*my,1);
sol.py    = zeros(mx*my,1);
sol.un    = zeros(mx*my,1);
sol.ux    = zeros(mx*my,1);
sol.uy    = zeros(mx*my,1);
sol.srel  = zeros(mx*my,1);
if (has_plast)
   sol.taucrt = zeros(mx*my,1);
   sol.uplsx  = zeros(mx*my,1);
   sol.uplsy  = zeros(mx*my,1);
end
if (has_heat)
   sol.temp1 = zeros(mx*my,1);
   sol.temp2 = zeros(mx*my,1);
end

if (fmtmat<1220)
   idum=0;
else
   idum=1; % introduced mu in 4th column in version 12.2
end

sol.eldiv(ir) = tmp(irow:iend, 2);
sol.h(ir)     = tmp(irow:iend, 3);
if (fmtmat>=1220)
   sol.mu(ir)    = tmp(irow:iend, 4);
end
sol.pn(ir)    = tmp(irow:iend,3+idum+1);
sol.px(ir)    = tmp(irow:iend,3+idum+2);
sol.py(ir)    = tmp(irow:iend,3+idum+3);
sol.un(ir)    = tmp(irow:iend,3+idum+4);
sol.ux(ir)    = tmp(irow:iend,3+idum+5);
sol.uy(ir)    = tmp(irow:iend,3+idum+6);
sol.srel(ir)  = tmp(irow:iend,3+idum+7);
iofs = 3+idum+7;
if (has_plast)
   sol.taucrt(ir) = tmp(irow:iend,iofs+1);
   sol.uplsx(ir)  = tmp(irow:iend,iofs+2);
   sol.uplsy(ir)  = tmp(irow:iend,iofs+3);
   iofs = iofs + 3;
end
if (has_heat)
   sol.temp1(ir) = tmp(irow:iend,iofs+1);
   sol.temp2(ir) = tmp(irow:iend,iofs+2);
   iofs = iofs + 2;
end

sol.eldiv = reshape(sol.eldiv,mx,my)';
sol.h     = reshape(sol.h    ,mx,my)';
sol.mu    = reshape(sol.mu   ,mx,my)';
sol.pn    = reshape(sol.pn   ,mx,my)';
sol.px    = reshape(sol.px   ,mx,my)';
sol.py    = reshape(sol.py   ,mx,my)';
sol.un    = reshape(sol.un   ,mx,my)';
sol.ux    = reshape(sol.ux   ,mx,my)';
sol.uy    = reshape(sol.uy   ,mx,my)';
sol.srel  = reshape(sol.srel ,mx,my)';
if (has_plast)
   sol.taucrt = reshape(sol.taucrt,mx,my)';
   sol.uplsx  = reshape(sol.uplsx ,mx,my)';
   sol.uplsy  = reshape(sol.uplsy ,mx,my)';
end
if (has_heat)
   sol.temp1  = reshape(sol.temp1 ,mx,my)';
   sol.temp2  = reshape(sol.temp2 ,mx,my)';
end

% in rolling, compute the shift from the slip velocity
% in shifts, fill the shft-array by copying (dq==1)

if (sol.kincns.t_digit>=2)
   sol.shft = sol.srel * sol.kincns.dq;
else
   sol.shft = sol.srel;
end

% compute the traction bound

if (~has_plast)
   sol.trcbnd =     sol.mu .* sol.pn;
else
   sol.trcbnd = min(sol.mu .* sol.pn, sol.taucrt);
end

end % function fill_struct

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

