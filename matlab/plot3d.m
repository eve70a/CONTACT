
function [ opts3 ] = plot3d(sol, opt3, prr, prw, subs)

%
% function [ opts3 ] = plot3d([sol], [opt3], [prr], [prw], [subs])
%
% Make a plot of different quantities in the contact area.
%
%  sol    == structure with CONTACT results as returned by loadcase.
%            may be empty; may be array of sol structures (#patches>1).
%
%  opt3, opts3 == structure with options for plot3d.
%           A valid structure opts3 is returned when plot3d is called with
%           no input arguments.
%
%  prr, prw == rail and wheel profiles (structs) or profile filenames. 
%           prr may be a variable profile in slcs-format.
%
%  subs   == structure with CONTACT results as returned by loadstrs.
%            may be empty; may be array of subs structures (#patches>1, 1 block/patch).
%
%  opt3.field: 
%   'eldiv', 'eldiv_spy', 'eldiv_contour': different ways of presenting the element division;
%   'h':        the undeformed distance of the two bodies;
%   'mu', 'taucrt':   the actual local coefficient of friction and the critical shear stress
%               for the tangential plasticity model;
%   'pn', 'px', 'py': the normal and tangential tractions acting on body (1),
%               i.e. the rail in module 1;
%   'ptabs', 'ptarg': the magnitude and direction of the tangential tractions;
%   'ptvec':    a vector-plot of the tangential tractions;
%   'ptabs+vec': show magnitude and direction of tangential tractions in one plot;
%   'un', 'ux', 'uy': the displacement differences in normal and tangential directions;
%   'uplsx', 'uplsy': the components of the accumulated plastic deformation (if M=4);
%   'utabs+vec', 'upls+vec': elastic and plastic displacements in tangential directions,
%               magnitude and orientation;
%   'sx', 'sy': the components of the local shift (T=1) or relative micro-slip velocity (T=2,3);
%   'shft', 'sabs', 'srel': the magnitude of the shift or absolute or relative slip velocity;
%   'shft+vec', 'sabs+vec', 'srel+vec': show in one plot the magnitude and direction of the
%               local shift resp. the absolute or relative slip velocity;
%   'fricdens': the frictional power density;
%   'temp1', 'temp2': surface temperatures of bodies 1 and 2 (if H>=1).
%
%  Other options:
%    rw_surfc: 'none' for flat view, 'prr' for rail coordinates, 'prw' for wheel surface view,
%               or 'both' for track coordinates;
%    exterval: the value to be plotted at points of the exterior area;
%    typplot:  type of plot: 'surf' (default), 'contourf', 'rw_rear' or 'rw_side';
%    view:     the view direction, e.g. [30 30], or 'rail' for [90 -90], 'rdyn' for [180 90];
%    xrange:   the range of the x-axis to be displayed in 3D plots
%    xysteps:  one or two step sizes for sampling surfaces along the x- and y-axes
%    zrange:   the range of the z-axis to be displayed;
%    ixrange:  the selection of the elements [1,mx] in x-direction to be displayed in the plot;
%    iyrange:  the selection of the elements [1,my] in y-direction to be displayed in the plot;
%    numvecx:  the maximum number of vectors to be displayed in x-direction;
%    numvecy:  same as numvecx for y-direction;
%    xstretch: stretching of x-axis esp. for vectors on elongated contact patches (b>>a)
%    vecscale: manual scaling factor ([mm] per [N/mm2]) for vectors on vector-plots.
%               used also for quantities in 'rw_rear' and 'rw_side' views
%    veccolor: color specification for vectors on vector-plots. Default: 'b' for ptvec, 'k' for ptabs+vec;
%               used also for quantitites in 'rw_rear' and 'rw_side' views, default matlab_color(6)
%    vecwidth: line-width for vectors on vector-plots;
%    addeldiv: show contours of element divisions on top of results (1) or not (0)
%              (only in 2D plots, e.g. view [0 90])
%    eldivcol: color specification for element division, e.g. ['b';'m';'g']
%              or [0 0 0.56; 0.5 0 0; 0.8 0.25 0.42]
%    eldivwid: line-width for contours for element division;
%    colormap: this changes the colormap for 3D plots;
%    addplot:  clear (0) or do not clear (1) the figure before plotting.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% determine version information on the Matlab graphics system

old_graphics = verLessThan('matlab', '8.4.0'); % 8.4 == R2014b
new_graphics = ~old_graphics;

% construct a local struct "myopt" in which default values are filled in for
% all available options.

myopt = struct( ...
   'field',    'default', ...
   'rw_surfc', 'none', ...
   'exterval',  NaN, ...
   'typplot',  'surf', ...
   'view',     'default', ...
   'xrange',    [], ...
   'xysteps',   [], ...
   'zrange',    [], ...
   'ixrange',  'auto', ...
   'iyrange',  'auto', ...
   'numvecx',   15, ...
   'numvecy',   15, ...
   'xstretch',  1, ...
   'vecscale',  [], ...
   'veccolor',  [], ...
   'vecwidth',  [], ...
   'addeldiv',  [], ...
   'eldivcol',  [0 0 0.56; 0.5 0 0; 0.8 0.25 0.42], ...
   'eldivwid',  [], ...
   'colormap', 'parula', ...
   'addplot',    0  ... % clear plot (0) or add to exist.plot (1, experimental)
);
if (old_graphics)
   myopt.colormap = 'jet';
end

% If the user has not supplied any arguments,
%    return default options as first and only output argument

if (nargin<1)
   opts3 = myopt;
   return
end

% If the user has not supplied an opt3-struct, use the default

if (nargin<2 | isempty(opt3))
   opt3 = myopt;
end

if (nargin<3)
   prr = [];
end
if (nargin<4)
   prw = [];
end
if (nargin<5)
   subs = [];
end

% Check whether user-supplied opt3-struct contains unknown options

useropts = fieldnames(opt3);
ierror = 0;
for i = 1:length(useropts)
   if (~isfield(myopt, useropts{i}))
      ierror = ierror + 1;
      disp(sprintf('Unknown option "%s" will be ignored',useropts{i}));
      if (ierror==1)
         disp(sprintf('You may use "opt3=rmfield(opt3,''%s'');" to remove the option',useropts{i}));
      end
   end
end

% Overwrite all values in "myopt" with user-supplied values

myopts = fieldnames(myopt);
for i = 1:length(myopts)
   if (isfield(opt3, myopts{i}))
      myopt = setfield(myopt, myopts{i}, mygetfield(opt3,myopts{i}));
   end
end

% Set default field

if (strcmp(myopt.field, 'default'))
   if (isempty(sol) | sol(1).kincns.t_digit<=0)
      myopt.field = 'pn';
   else
      myopt.field = 'ptabs+vec';
   end
end

% Check whether requested field is valid

ok_fields=strvcat('eldiv', 'eldiv_spy', 'eldiv_contour', 'h', 'mu', ...
                  'pn', 'px', 'py', 'ptabs', 'ptarg', 'ptvec', 'ptabs+vec', ...
                  'un', 'ux', 'uy', 'utabs+vec', ...
                  'taucrt', 'uplsx', 'uplsy', 'upls+vec', ...
                  'sx', 'sy', 'shft', 'sabs', 'srel', ...
                  'shft_vec', 'sabs_vec', 'srel_vec', ...
                  'shft+vec', 'sabs+vec', 'srel+vec', 'fricdens', ...
                  'temp1', 'temp2');
found=0;
for istr=1:size(ok_fields,1)
   if (strcmp(myopt.field, deblank(ok_fields(istr,:))))
      found=1;
   end
end
user_added_fld = 0;
if (~found & isfield(sol, myopt.field))
   % disp(sprintf('Plotting user-added field "%s"', myopt.field));
   user_added_fld = 1;
   found = 1;
end
if (~found)
   disp(sprintf('Unknown field requested="%s"; Available options:',myopt.field));
   disp(ok_fields)
   return;
end

% Return if no solution is given (handy for w/r contact: call plot3d for fixed #patches)

if (isempty(sol))
   % disp('empty sol, returning')
   return;
end

% prepare the rail profile, if provided

if (ischar(prr))
   is_wheel = 0; mirror_y = 0; % note: Miniprof profiles usually need to be mirrored
   prr = read_profile(prr, is_wheel, mirror_y);
end
if (~isempty(prr) & ~isfield(prr,'nslc'))
   if (isempty(prr.ProfileY) | isempty(prr.ProfileZ))
      disp(sprintf('ERROR: rail profile does not contain any data.'));
      prr = [];
   elseif (prr.ProfileY(end) - prr.ProfileY(1) < 0)
      disp(sprintf('ERROR: rail profile Y-coordinates must be ascending.'));
      return
   end
end

% prepare the wheel profile, if provided

if (ischar(prw))
   is_wheel = 1; mirror_y = 0;
   prw = read_profile(prw, is_wheel, mirror_y);
end
if (~isempty(prw))
   if (isempty(prw.ProfileY) | isempty(prw.ProfileZ))
      disp(sprintf('ERROR: wheel profile does not contain any data.'));
      prw = [];
   end
   if (~isempty(prw) & prw.ProfileY(end) - prw.ProfileY(1) > 0)
      disp(sprintf('ERROR: wheel profile Y-coordinates must be descending.'));
      return
   end
end

% Check whether requested rail or wheel surface is available

if (~any(strcmp(myopt.rw_surfc,{'none','prr','prw','both'})))
   disp(sprintf('Unknown option rw_surfc="%s" will be ignored',myopt.rw_surfc));
   disp(sprintf('Available values are "none" (default), "prr", "prw", and "both".'));
   myopt.rw_surfc = 'none';
end

if (isempty(prw) & any(strcmp(myopt.rw_surfc,{'prw','both'})))
   disp(sprintf('No wheel profile provided; ignoring rw_surfc="%s". ', myopt.rw_surfc));
   if (strcmp(myopt.rw_surfc, 'both'))
      myopt.rw_surfc = 'prr';
   else
      myopt.rw_surfc = 'none';
   end
end
if (isempty(prr) & any(strcmp(myopt.rw_surfc,{'prr','both'})))
   disp(sprintf('No rail profile provided; ignoring rw_surfc="%s".', myopt.rw_surfc));
   if (strcmp(myopt.rw_surfc, 'both') & ~isempty(prw))
      myopt.rw_surfc = 'prw';
   else
      myopt.rw_surfc = 'none';
   end
end

if (~strcmp(myopt.rw_surfc,'none') & any(strcmp(myopt.field, {'eldiv_spy'})))
   disp(sprintf('Option "%s" for rw_surfc is not available for field "%s"', myopt.rw_surfc, myopt.field));
   myopt.rw_surfc = 'none';
end

% Check whether requested typplot is valid, set default when needed

if (~any(strcmp(myopt.typplot,{'surf','contourf','rw_rear','rw_side'})))
%    typplot:  type of plot: 'surf' (default), 'contourf', 'rw_rear' or 'rw_side';
   disp(sprintf('Unknown value "%s" for option typplot will be ignored.',myopt.typplot));
   disp(sprintf('Available values are "surf" (default), "contourf", "rw_rear", and "rw_side".'));
   myopt.typplot = 'surf';
end
if (any(strcmp(myopt.typplot,{'rw_rear','rw_side'})) & ~any(strcmp(myopt.rw_surfc,{'prr','prw','both'})))
   disp(sprintf('Option rw_surfc="%s" therefore option typplot="%s" will be ignored', myopt.rw_surfc, ...
                myopt.typplot));
   myopt.typplot = 'surf';
end

% Check whether requested view is valid, set default when needed

if (isempty(myopt.view)), myopt.view = 'default'; end
if (ischar(myopt.view))
   if (strcmp(myopt.view,'rail') & strcmp(myopt.typplot,'surf') & ~strcmp(myopt.rw_surfc,'none'))
      % 3d surface plot should not use 2d view
      myopt.view = 'default';
   end
   if (strcmp(myopt.view,'rail'))
      % wheel-rail view: 2D, rolling direction == plot y-direction
      myopt.view=[90 -90];
   elseif (strcmp(myopt.view,'rdyn'))
      % RecurDyn view: 2D, x- and y-axes reversed
      myopt.view=[180 90];
   elseif (strcmp(myopt.view,'default') & ~strcmp(myopt.rw_surfc,'none'))
      % default: 3D wheel-rail profile view
      if (is_left_side(sol))
         myopt.view=[115 30];
      else
         myopt.view=[65 30];
      end
   else
      % default: 2D or 3D view depending on field to be plotted
      if (~strcmp(myopt.view,'default'))
         disp(sprintf('Unknown option view="%s" will be ignored',myopt.view));
         opt.view = 'default';
      end
      if (strcmp(myopt.field,'h')        || ...
          strcmp(myopt.field,'fricdens') );
         myopt.view=[-37.5 30];
      else
         myopt.view=[0 90]; % traditional view, rolling to the right
         % myopt.view=[90 -90]; % rail view, rolling up in the picture
      end
   end 
end

if (isempty(myopt.colormap))
   myopt.colormap = 'default'; % matlab built-in
end

% fill in default for addeldiv (eldiv_contour)

if (isempty(myopt.addeldiv))
   if (strcmp(myopt.field,'eldiv') | strcmp(myopt.field,'eldiv_spy'))
      myopt.addeldiv = 0;
   else
      myopt.addeldiv = 1;
   end
end
if (ischar(myopt.addeldiv))
   if (strcmp(myopt.addeldiv,'n') | strcmp(myopt.addeldiv,'no'))
      myopt.addeldiv = 0;
   elseif (strcmp(myopt.addeldiv,'y') | strcmp(myopt.addeldiv,'yes'))
      myopt.addeldiv = 1;
   else
      disp(sprintf('Unknown value "%s" for addeldiv will be ignored', ...
                                                          myopt.addeldiv));
      myopt.addeldiv = 0;
   end
end
if (any(strcmp(myopt.typplot,{'rw_rear','rw_side'})))
   myopt.addeldiv = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot arrays of sol structures using recursion (esp. w/r contacts with multiple patches)

if (isstruct(sol) & length(sol)>1)
   if (strcmp(myopt.rw_surfc,'none') | (isempty(prr) & isempty(prw)))
      disp(sprintf('ERROR: plotting %d patches at once requires using rail or wheel profiles', length(sol)));
      return
   end

   for isol = 1 : length(sol)
      pmax(isol) = max(max(sol(isol).pn));
   end
   myopt.vecscale = get_vecscale( sol, myopt, pmax ); % fix vecscale -- same for all patches

   [~,jsol] = sort(pmax, 'descend');    % process highest->lowest -- axes based on first one
   for i = 1 : length(sol)
      plot3d( sol(jsol(i)), myopt, prr, prw );
      myopt.addplot = 1;
   end
   return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add profiles to the sol struct

if (any(strcmp(myopt.rw_surfc, {'prr','both'})))
   if (isfield(prr,'nslc')) % variable profile
      sol.slcs = prr;
      sol.prr = get_profile_slice(sol.slcs, sol.meta.s_ws+sol.meta.xcp_r);
   else
      sol.prr = prr;
   end
end
if (any(strcmp(myopt.rw_surfc, {'prw','both'})))
   sol.prw = prw;
end
clear prr prw;

% Add subsurface data to the sol struct

if (~isempty(subs))
   sol.subs = subs;
end

% Add defaults for missing components of sol - for convenience of C.library

if (~isfield(sol, 'x_offset')), sol.x_offset = []; end
if (~isfield(sol, 'y_offset')), sol.y_offset = []; end
if (~isfield(sol, 'kincns.t_digit')), sol.kincns.t_digit = 3; end
if (~isfield(sol, 'x'))
   sol.x = sol.xl + ([1:sol.mx]-0.5)*sol.dx;
end
if (~isfield(sol, 'y'))
   sol.y = sol.yl + ([1:sol.my]-0.5)*sol.dy;
end

% Check sizes of main solution arrays - for convenience of C. library

if (any(size(sol.pn)~=[sol.my sol.mx]))
   disp(sprintf(['ERROR: incorrect size for array pn (%d,%d), ', ...
                        'expecting (%d,%d)'], size(sol.pn), sol.my, sol.mx));
   return
end
if (any(size(sol.eldiv)~=[sol.my sol.mx]))
   disp(sprintf(['ERROR: incorrect size for array eldiv (%d,%d), ', ...
                        'expecting (%d,%d)'], size(sol.eldiv), sol.my, sol.mx));
   return
end

% Prefix for figure title, dependent on type of problem

pfxtit = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start actual processing to produce the requested plot
%

if (myopt.addplot<=0)
   clf;  % restart figure
elseif (myopt.addplot==2)
   cla;  % restart subplot?
end

% start by plotting the rail and/or wheel surfaces

if (myopt.addplot<=0 & ~strcmp(myopt.rw_surfc,'none'))
   show_profiles(sol, myopt);
elseif (isfield(sol,'slcs') & any(strcmp(myopt.rw_surfc,{'prr','both'})) & ...
                              any(strcmp(myopt.typplot,{'rw_side','rw_rear'})) )
   show_profiles(sol, myopt);
end
hold on;

% set element numbers ix, iy to be plotted: interval of [1:mx(my)]

[ix_plot, iy_plot] = plot_ranges(sol, myopt);

% make the plot, depending on the field that is requested

if (strcmp(myopt.field,'h'))

   show_scalar_field(sol, 'h', myopt);
   title ('Undeformed distance H');

end
if (strcmp(myopt.field,'mu'))

   show_scalar_field(sol, 'mu', myopt);
   title ('Actual local friction coefficient \mu(x)');

end
if (strcmp(myopt.field,'pn'))

   show_scalar_field(sol, 'pn', myopt);
   title ('Normal pressure P_n');

end
if (strcmp(myopt.field,'px'))

   show_scalar_field(sol, derived_data(sol, myopt, 'px'), myopt);
   title ([pfxtit, 'Tangential traction P_x']);

end
if (strcmp(myopt.field,'py'))

   show_scalar_field(sol, derived_data(sol, myopt, 'py'), myopt);
   title ([pfxtit, 'Tangential traction P_y']);

end
if (strcmp(myopt.field,'ptabs') | strcmp(myopt.field,'ptabs+vec'))

   pt = derived_data(sol, myopt, 'pt');
   show_scalar_field(sol, pt, myopt);
   title ([pfxtit, 'Magnitude of tangential tractions |P_t|']);

end
if (strcmp(myopt.field,'ptvec') | strcmp(myopt.field,'ptabs+vec'))

   [px, py] = derived_data(sol, myopt, 'px', 'py');
   show_vec_field(sol, px, py, myopt);
   title ([pfxtit, 'Tangential tractions p_t']);

end
if (strcmp(myopt.field,'ptarg'))

   [eldiv, ptarg] = derived_data(sol, myopt, 'eldiv', 'ptarg');
   % ptarg: range -180:180 deg
   % shift lower bound (-180) to "arg_low", to avoid wrap-around in the plot
   msk=(eldiv>0);
   cnt(1) = nnz(msk.*(ptarg> 135 | ptarg<-135));
   cnt(2) = nnz(msk.*(ptarg>-135 & ptarg<-45));
   cnt(3) = nnz(msk.*(ptarg> -45 & ptarg< 45));
   cnt(4) = nnz(msk.*(ptarg>  45 & ptarg<135));
   [mn,ix]=min(cnt); arg_low=-270+90*ix;
   ix=find(ptarg<arg_low); ptarg(ix)=ptarg(ix)+360;
   ix=find(ptarg>arg_low+360); ptarg(ix)=ptarg(ix)-360;

   show_scalar_field(sol, ptarg, myopt);
   if (isempty(myopt.zrange))
      ix = find(eldiv>=1);
      rng=[ min(min(ptarg(ix))), max(max(ptarg(ix))) ];
      rng=rng + 0.05*[-1 1]*(rng(2)-rng(1));
      set(gca,'clim', rng);
   end
   title ([pfxtit, 'Direction of tangential tractions arg(P_t)']);

end
if (strcmp(myopt.field,'un'))

   show_scalar_field(sol, 'un', myopt);
   title ('Normal displacement difference U_n');

end
if (strcmp(myopt.field,'ux'))

   show_scalar_field(sol, 'ux', myopt);
   title ([pfxtit, 'Tangential displacement difference U_x']);

end
if (strcmp(myopt.field,'uy'))

   show_scalar_field(sol, 'uy', myopt);
   title ([pfxtit, 'Tangential displacement difference U_y']);

end
if (strcmp(myopt.field,'utabs+vec'))

   ut = sqrt(sol.ux.^2 + sol.uy.^2);
   show_scalar_field(sol, ut, myopt);
   show_vec_field(sol, sol.ux, sol.uy, myopt);
   title ('Tangential displacements u_t');

end
if (strcmp(myopt.field,'uplsx'))

   show_scalar_field(sol, 'uplsx', myopt);
   title ('Plastic deformation u_{pl,x}');

end
if (strcmp(myopt.field,'uplsy'))

   show_scalar_field(sol, 'uplsy', myopt);
   title ('Plastic deformation u_{pl,y}');

end
if (strcmp(myopt.field,'upls+vec'))

   upls = sqrt(sol.uplsx.^2 + sol.uplsy.^2);
   show_scalar_field(sol, upls, myopt);
   show_vec_field(sol, sol.uplsx, sol.uplsy, myopt);
   title ('Plastic deformation u_{pl}');

end
if (strcmp(myopt.field,'taucrt'))

   show_scalar_field(sol, 'taucrt', myopt);
   title ('Tangential yield stress \tau_{crt}');

end

% compute sx, sy when needed

if (strcmp(myopt.field,'sx')       | strcmp(myopt.field,'sy') | ...
    strcmp(myopt.field,'shft_vec') | strcmp(myopt.field,'shft+vec') | ...
    strcmp(myopt.field,'sabs_vec') | strcmp(myopt.field,'sabs+vec') | ...
    strcmp(myopt.field,'srel_vec') | strcmp(myopt.field,'srel+vec'))

   if (isfield(sol,'sx') & isfield(sol,'sy'))
      sx = sol.sx;
      sy = sol.sy;
   else
      pxdir = sol.px ./ max(1d-9,sqrt(sol.px.^2 + sol.py.^2));
      pydir = sol.py ./ max(1d-9,sqrt(sol.px.^2 + sol.py.^2));
      if (sol.kincns.t_digit>=2)
         sx = -sol.srel .* pxdir;
         sy = -sol.srel .* pydir;
      else
         sx = -sol.shft .* pxdir;
         sy = -sol.shft .* pydir;
         if (strfind(myopt.field,'sabs'))
            sx = sx * sol.kincns.veloc;
            sy = sy * sol.kincns.veloc;
         end
      end
   end
end

if (strcmp(myopt.field,'sx'))

   show_scalar_field(sol, sx, myopt);
   if (sol.kincns.t_digit>=2)
      title ([pfxtit, 'Relative slip component s_x']);
   else
      title ([pfxtit, 'Shift component S_x']);
   end
end
if (strcmp(myopt.field,'sy'))

   show_scalar_field(sol, sy, myopt);
   if (sol.kincns.t_digit>=2)
      title ([pfxtit, 'Relative slip component s_y']);
   else
      title ([pfxtit, 'Shift component S_y']);
   end

end
if (strcmp(myopt.field,'shft') | strcmp(myopt.field,'shft+vec'))

   show_scalar_field(sol, 'shft', myopt);
   title ([pfxtit, 'Shift distance']);

end
if (strcmp(myopt.field,'sabs') | strcmp(myopt.field,'sabs+vec'))

   sabs = derived_data(sol, myopt, 'sabs');
   show_scalar_field(sol, sabs, myopt);
   title ([pfxtit, 'Absolute slip velocity']);

end
if (strcmp(myopt.field,'srel') | strcmp(myopt.field,'srel+vec'))

   show_scalar_field(sol, 'srel', myopt);
   title ([pfxtit, 'Relative slip velocity']);

end
if (strcmp(myopt.field,'shft_vec') | strcmp(myopt.field,'shft+vec') | ...
    strcmp(myopt.field,'sabs_vec') | strcmp(myopt.field,'sabs+vec') | ...
    strcmp(myopt.field,'srel_vec') | strcmp(myopt.field,'srel+vec'))

   show_vec_field(sol, sx, sy, myopt)

   if (strcmp(myopt.field,'shft_vec') | strcmp(myopt.field,'shft+vec'))
      title ([pfxtit, 'Tangential shift S_\tau']);
   elseif (strcmp(myopt.field,'sabs_vec') | strcmp(myopt.field,'sabs+vec'))
      title ([pfxtit, 'Absolute slip velocity s_a']);
   elseif (strcmp(myopt.field,'srel_vec') | strcmp(myopt.field,'srel+vec'))
      title ([pfxtit, 'Relative slip velocity s_r']);
   end

end
if (strcmp(myopt.field,'fricdens'))

   % unit [W / mm2] = ([mm/s]/1000) * [-]  * [ N/mm2],   [W] = [J/s] = [N.m/s]
   % note: srel is dimensionless, using V.dt=dq=1 in shifts (T=1)
   fricdens = sol.kincns.veloc/1000 * sol.srel .* sqrt(sol.px.^2 + sol.py.^2);
   show_scalar_field(sol, fricdens, myopt);
   title ([pfxtit, 'Frictional power density']);

end
if (strcmp(myopt.field,'temp1'))

   show_scalar_field(sol, 'temp1', myopt);
   title ('Surface temperature T^{(1)}');

end
if (strcmp(myopt.field,'temp2'))

   show_scalar_field(sol, 'temp2', myopt);
   title ('Surface temperature T^{(2)}');

end
if (user_added_fld & isfield(sol, myopt.field))

   show_scalar_field(sol, myopt.field, myopt);
   title (sprintf('User-added field "%s"', myopt.field));

end
if (strcmp(myopt.field,'eldiv'))

   myopt.zrange = [-0.1 2.1];
   eldiv = derived_data(sol, myopt, 'eldiv');
   ix = find(eldiv==3); eldiv(ix) = 1.5; % set Plast(3) between Adhes(1) and Slip(2)
   show_scalar_field(sol, eldiv, myopt);
   h = get(gcf,'children'); ix = strmatch('Colorbar', get(h,'tag'));
   if (length(ix)==1)
      jx = find(sol.eldiv==3);
      if (~isempty(jx))
         set(h(ix),'ytick',[0,1,1.5,2],'yticklabel',['exter';'stick';'plast';'slip ']);
      else
         set(h(ix),'ytick',[0,1,2],'yticklabel',['exter';'stick';'slip ']);
      end
      if (isnan(myopt.exterval))
         set(h(ix),'ylim',[0.9 2.1]);
      end
   end
   title ([pfxtit, 'The element divisions']);

end
if (strcmp(myopt.field,'eldiv_spy'))

   eldiv = derived_data(sol, myopt, 'eldiv');
   mask = (eldiv==1); spy(mask, 'b.', 6); % ones in the adhesion area
   mask = (eldiv==2); spy(mask, 'r.', 6); % ones in the slip area
   mask = (eldiv==3); spy(mask, 'g.', 6); % ones in the platicity area
   title ([pfxtit, 'The element divisions (red=slip, blue=adhesion)']);
   axis xy
   ix = get(gca,'xtick')'; iy = get(gca,'ytick')';
   xcentr = sol.xl+(ix+0.5)*sol.dx;
   ycentr = sol.yl+(iy+0.5)*sol.dy;
   if (~isempty(sol.x_offset)), xcentr = xcentr + sol.x_offset; end
   if (~isempty(sol.y_offset)), ycentr = ycentr + sol.y_offset; end
   set(gca,'xticklabel', num2str(xcentr, '%5.3f'));
   set(gca,'yticklabel', num2str(ycentr, '%5.3f'));
   set_xyzlabels(sol, myopt);
   if (~is_2d_view(gca, myopt.view))
      v=axis; axis([v -1e-6 1e-6]); axis equal; 
      set(gca,'ztick',[]);
   end
   view(myopt.view);
end
if ( strcmp(myopt.field,'eldiv_contour') | myopt.addeldiv )
   % eldiv_contour can be used in 2D view & on rail/wheel surface
   if (is_2d_view(gca,myopt.view) | ~strcmp(myopt.rw_surfc,'none') )
      show_eldiv_contour(sol, myopt);
   end
end
if (any(strcmp(myopt.rw_surfc,{'prr','both'})) & ~any(strcmp(myopt.typplot,{'rw_rear','rw_side'})))
   h = findobj(gcf, 'Tag','prr');
   clim = get(gca,'clim');
   set(h,'cdata', clim(1)*ones(size(get(h,'xdata'))));
end
if (any(strcmp(myopt.rw_surfc,{'prw','both'})) & ~any(strcmp(myopt.typplot,{'rw_rear','rw_side'})))
   h = findobj(gcf, 'Tag','prw');
   clim = get(gca,'clim');
   set(h,'cdata', clim(2)*ones(size(get(h,'xdata'))));
end

end % function plot3d

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ ] = show_profiles(sol, opt)
% 
% function [ ] = show_profiles(sol, opt)
%
% create 3d surf plot of rail and/or wheel profiles
%   sol     - structure with CONTACT results as returned by loadcase
%   opt     - structure with options as defined by plot3d.
%

   % set the x range for plotting

   if (isempty(opt.xrange))
      if (any(strcmp(opt.rw_surfc,{'prr','both'})))

         xmax = max( abs(sol.meta.xcp_r) + 0.40 * sol.mx * sol.dx, 1.20 * sol.mx * sol.dx);
         if (isfield(sol,'slcs')), xmax = 3 * xmax; end

      else

         xmax = max( abs(sol.meta.xcp_w) + 0.40 * sol.mx * sol.dx, 1.20 * sol.mx * sol.dx);
      end

      opt.xrange = [-xmax, xmax ];
   end

   % set target step-sizes for plotting

   if (isempty(opt.xysteps))
      xlen = opt.xrange(2) - opt.xrange(1);
      ylen = 0;
      if (any(strcmp(opt.rw_surfc,{'prr','both'})))
         ylen = max(ylen, max(sol.prr.ProfileY) - min(sol.prr.ProfileY));
      end
      if (any(strcmp(opt.rw_surfc,{'prw','both'})))
         ylen = max(ylen, max(sol.prw.ProfileY) - min(sol.prw.ProfileY));
      end
      opt.xysteps = max(xlen/4, ylen/15);
   end
   if (length(opt.xysteps)==1)
      opt.xysteps = [1 1] * opt.xysteps;
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % if requested, show the rail profile

   want_rail = 1;

   if (any(strcmp(opt.rw_surfc,{'prr','both'})) & strcmp(opt.typplot,'rw_rear'))

      % 2d rear view of y-z-plane

      % form profile curve at x_cp

      xval = sol.meta.xcp_r;
      [ xcurv, ycurv, zcurv ] = make_3d_surface( sol, opt, want_rail, [xval, xval], [], 1, ...
                                                                        [], [], []);

      plot(ycurv, zcurv, 'color',matlab_color(1));
      set(gca,'ydir','reverse');
      axis equal;
      hold on

   elseif (any(strcmp(opt.rw_surfc,{'prr','both'})) & strcmp(opt.typplot,'rw_side'))

      % 2d side view of x-z-plane

      % s-value with y(ip) closest to y_cp

      [~,ip] = min( abs(sol.prr.ProfileY - sol.meta.ycp_r) );
      sval = sol.prr.ProfileS(ip);

      % form profile curve at given x, s

      totlen = opt.xrange(2) - opt.xrange(1);
      [ xcurv, ycurv, zcurv ] = make_3d_surface( sol, opt, want_rail, opt.xrange, [], opt.xysteps(1)/10, ...
                                                                        [sval, sval], [], 1);

      plot(xcurv, zcurv, 'color',matlab_color(1));
      set(gca,'ydir','reverse');
      axis equal;
      hold on

   elseif (any(strcmp(opt.rw_surfc,{'prr','both'})))

      % 3d view of rail surface

      % form profile surface at given x, all s - fine sampling

      % disp('...full surface')
      [ xsurf, ysurf, zsurf ] = make_3d_surface( sol, opt, want_rail, opt.xrange, [], opt.xysteps(1)/10, ...
                                                                           [], [], []);

      % plot 90% transparent surface (alpha=0.1); color-value?

      csurf = -1 * ones(size(xsurf));
      l=surf(xsurf, ysurf, zsurf, csurf, 'EdgeAlpha',0.1, 'FaceAlpha',0.1);
      set(l,'tag','prr');

      % set view appropriate for rails with z positive downwards

      set(gca,'xdir','reverse', 'zdir','reverse');
      axis equal;
      hold on

      % determine transverse curves along the surface at fixed steps in x

      % disp('...transverse curves')
      [ xcurv, ycurv, zcurv ] = make_3d_surface( sol, opt, want_rail, opt.xrange, [], opt.xysteps(1), ...
                                                                           [], [], []);

      % plot transverse curves along the surface

      l = plot3(xcurv', ycurv', zcurv', 'color',[.5 .5 .5], 'linewidth',0.5);

      [~,imin] = min(abs(xcurv(:,1)));
      set(l(imin), 'color', [.3 .3 .3], 'linewidth',1);

      % determine longitudinal curves along the surface at fixed steps in s

      % disp('...longitudinal curves')
      [ xcurv, ycurv, zcurv ] = make_3d_surface( sol, opt, want_rail, opt.xrange, [], opt.xysteps(1)/10, ...
                                                                     [], [], opt.xysteps(2));

      % plot curves along surface in longitudinal direction

      clear l ip;
      l = plot3(xcurv, ycurv, zcurv, 'color',[.5 .5 .5], 'linewidth',0.5);
      imid = ceil(length(l)/2);
      set(l(imid), 'color', [.3 .3 .3], 'linewidth',1);

   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % if requested, show the wheel profile

   want_rail = 0;

   if (any(strcmp(opt.rw_surfc,{'prw','both'})) & strcmp(opt.typplot,'rw_rear'))

      % 2d rear view of y-z-plane

      if (sol.meta.npatch>1)
         xval = 0;      % use principal profile; xcp_w can be much different between patches
      else
         xval = sol.meta.xcp_w;
      end
      [ xcurv, ycurv, zcurv ] = make_3d_surface( sol, opt, want_rail, [xval, xval], [], 1, ...
                                                               [], [], []);

      plot(ycurv, zcurv, 'color',matlab_color(3));
      set(gca,'ydir','reverse');
      axis equal;
      hold on

   elseif (any(strcmp(opt.rw_surfc,{'prw','both'})) & strcmp(opt.typplot,'rw_side'))

      % 2d side view of x-z-plane

      % s-value with y(ip) closest to y_cp

      [~,ip] = min( abs(sol.prw.ProfileY - sol.meta.ycp_w) );
      sval = sol.prw.ProfileS(ip);

      % form wheel profile curve at given x, s

      [ xcurv, ycurv, zcurv ] = make_3d_surface( sol, opt, want_rail, opt.xrange, [], opt.xysteps(1)/10, ...
                                                                        [sval, sval], [], 1);

      plot(xcurv, zcurv, 'color',matlab_color(3));
      set(gca,'ydir','reverse');
      axis equal;
      hold on

   elseif (any(strcmp(opt.rw_surfc,{'prw','both'})))

      % 3d view of wheel surface

      % x-range a little larger than the potential contact area

      xmax = max( abs(sol.meta.xcp_w) + 0.40 * sol.mx * sol.dx, ...
                                            1.20 * sol.mx * sol.dx);
      if (isfield(sol,'slcs')), xmax = 3 * xmax; end

      % form profile surface at given x, all s - fine sampling

      [ xsurf, ysurf, zsurf ] = make_3d_surface( sol, opt, want_rail, opt.xrange, [], opt.xysteps(1)/10, ...
                                                                                [], [], []);

      % plot 90% transparent surface (alpha=0.1); color-value?

      csurf = -1 * ones(size(xsurf));
      l=surf(xsurf, ysurf, zsurf, csurf, 'EdgeAlpha',0.1, 'FaceAlpha',0.1);
      set(l,'tag','prw');

      % set view appropriate for rails with z positive downwards

      set(gca,'xdir','reverse', 'zdir','reverse');
      axis equal;
      hold on

      % determine transverse curves along the surface at fixed steps in x

      [ xcurv, ycurv, zcurv ] = make_3d_surface( sol, opt, want_rail, opt.xrange, [], opt.xysteps(1), ...
                                                                                [], [], []);

      % plot transverse curves along the surface

      l = plot3(xcurv', ycurv', zcurv', 'color',[.5 .5 .5], 'linewidth',0.5);

      [~,imin] = min(abs(xcurv(:,1)));
      set(l(imin), 'color', [.3 .3 .3], 'linewidth',1);

      % determine curves along the surface at fixed steps in s

      [ xcurv, ycurv, zcurv ] = make_3d_surface( sol, opt, want_rail, opt.xrange, [], opt.xysteps(1)/10, ...
                                                                          [], [], opt.xysteps(2));

      % plot curves along surface in longitudinal direction

      clear l ip;
      l = plot3(xcurv, ycurv, zcurv, 'color',[.5 .5 .5], 'linewidth',0.5);
      imid = ceil(length(l)/2);
      set(l(imid), 'color', [.3 .3 .3], 'linewidth',1);

   end

end % function show_profiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ xsurf, ysurf, zsurf ] = make_3d_surface( sol, opt, want_rail, ...
                                                              xrange, dx_true, dx_appx, ...
                                                              srange, ds_true, ds_appx, idebug)

% create 2D arrays (x,y,z) for a prismatic rail (prr), variable rail (slcs), or roller (nom_radius)
%  - the profile will be plotted for xrange = [xmin, xmax] in track coordinates
%  - fixed steps of size dx_true may be requested, or n-time fitting steps of approximately dx_appx
%  - the whole profile will be used if srange = [smin, smax] is left empty
%  - the profile will be sampled at sj-values in the profile close to uniform sj = [smin: ds: smax]
%  - all profile points will be used if both ds_true and ds_appx are left empty
%  - for variable profiles, x-positions and arc-lengths s are replaced by the (u,v)-parametrization

   use_intern   = 1;
   if (nargin< 5), dx_true = []; end
   if (nargin< 6), dx_appx = []; end
   if (nargin< 7), srange  = []; end
   if (nargin< 8), ds_true = []; end
   if (nargin< 9), ds_appx = []; end
   if (nargin<10 | isempty(idebug)), idebug  =  0; end

   if (want_rail)
      prf = sol.prr;
      nom_radius = sol.meta.rnom_rol;
   else
      prf = sol.prw;
      nom_radius = sol.meta.rnom_whl;
   end
   has_slcs   = (want_rail & isfield(sol,'slcs'));

   % determine longitudinal positions xi (track coords) for evaluation of surface

   if (~isempty(dx_true))

      % if 'final' dx_true is prescribed,
      %    set xi = { i * dx_true } for appropriate i

      i0 = ceil( xrange(1) / dx_true );
      i1 = floor( xrange(2) / dx_true );
      xi = [i0 : i1] * dx_true;

   elseif (~isempty(dx_appx))

      % if 'target' dx_appx is given, 
      %    divide [xmin,xmax] in odd #intervals of size dx_true close to dx_appx

      nintv   = floor( (xrange(2) - xrange(1)) / dx_appx );
      if (mod(nintv,2)==1)  % prefer odd #lines for symmetric [-xmax, xmax]
         nintv = nintv + 1;
      end
      dx = (xrange(2) - xrange(1)) / max(1,nintv);
      xi = xrange(1) + [0:nintv] * dx;

   else

      disp('Error: either dx_true or dx_appx must be given');
      return

   end

   % determine lateral positions sj for evaluation of surface

   if (has_slcs)
      profile_s = sol.slcs.vj;
   else
      profile_s = prf.ProfileS;
   end

   if (isempty(srange))
      srange = [ profile_s(1), profile_s(end) ];
   end

   if (~isempty(ds_true))

      % if 'final' ds_true is prescribed,
      %    set si = { j * ds_true } for appropriate j

      j0 = ceil( (srange(1)-profile_s(1)) / ds_true );
      j1 = floor( (srange(2)-profile_s(1)) / ds_true );
      sj = [j0 : j1] * ds_true;

   elseif (~isempty(ds_appx))

      % if 'target' ds_appx is given, 
      %    divide [smin,smax] in odd #intervals of size ds_true close to ds_appx

      if (has_slcs)
         len = 1.3 * max(max(sol.slcs.ysurf)) - min(min(sol.slcs.ysurf));
         nintv = floor ( len / ds_appx );
      else
         nintv = floor( (srange(2) - srange(1)) / ds_appx );
      end
      if (mod(nintv,2)==1)
         nintv = nintv + 1;  % prefer odd #lines for symmetric [-smax, smax]
      end
      ds = (srange(2) - srange(1)) / max(1,nintv);
      sj = srange(1) + [0:nintv] * ds;

   else

      % if neither is given, use all profile points in range srange

      j  = find( profile_s>=srange(1) & profile_s<=srange(2) );
      sj = profile_s(j);

   end

   % determine indices js closest to each sj

   js = zeros(size(sj));
   for j = 1 : length(sj)
      [~,js(j)] = min(abs(profile_s - sj(j)));
   end

   % form profile surface at given xi and js

   nx    = length(xi);
   ns    = length(js);

   if (want_rail & ~has_slcs & abs(nom_radius)<1e-3)

      % prismatic rail

      xsurf = xi' * ones(1,ns);
      ysurf = ones(nx,1) * prf.ProfileY(js)';
      zsurf = ones(nx,1) * prf.ProfileZ(js)';
      
   elseif (want_rail & myisfield(sol, 'slcs.spl2d') & use_intern & exist('eval_2dspline'))

      % variable rail profile, using internal 2D spline algorithm

      % form profile surface at given ui (==x_fc) and vj

      sf_min = min(sol.slcs.spl2d.ui);
      sf_max = max(sol.slcs.spl2d.ui);
      ui     = max(sf_min, min(sf_max, sol.meta.s_ws + xi));
      vj     = sj;

      [ ~, xsurf, ysurf, zsurf ] = eval_2dspline( sol.slcs.spl2d, ui, vj );

   elseif (want_rail & has_slcs)

      xsurf = xi' * ones(1,ns);
      ysurf = []; zsurf = [];
      for ix = 1 : nx
         x_cur = sol.meta.s_ws + xi(ix);
         tmp_prr = get_profile_slice(sol.slcs, x_cur);
         ysurf = [ysurf; tmp_prr.ProfileY(js)'];
         zsurf = [zsurf; tmp_prr.ProfileZ(js)'];
      end

   elseif (want_rail)

      % roller surface: circle x^2 + (rnom - z)^2 = (rnom - z(0,y))^2

      xsurf = xi' * ones(1,ns);
      ysurf = ones(nx,1) * prf.ProfileY(js)';
      r_y   = nom_radius - prf.ProfileZ(js)';
      zsurf = nom_radius - sqrt( max(0, (ones(nx,1)*r_y).^2 - xsurf.^2) );

   else

      % wheel surface: circle x^2 + (rnom + z)^2 = (rnom + z(0,y))^2

      xsurf = xi' * ones(1,ns);
      ysurf = ones(nx,1) * prf.ProfileY(js)';
      r_y   =  nom_radius + prf.ProfileZ(js)';
      zsurf = -nom_radius + sqrt( max(0, (ones(nx,1)*r_y).^2 - xsurf.^2) );

   end

   % profiles are given for right-side w/r combination
   % mirror y for left-side w/r combination

   if (is_left_side(sol))
      ysurf = -ysurf;
   end

   % transform to track coordinates if needed

   if (want_rail & strcmp(opt.rw_surfc,'both'))

      [xsurf, ysurf, zsurf] = rail_to_track_coords(sol, xsurf, ysurf, zsurf);

   elseif (strcmp(opt.rw_surfc,'both'))

      [xsurf, ysurf, zsurf] = wheel_to_track_coords(sol, xsurf, ysurf, zsurf);

   end

end % function make_3d_surface

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ ] = show_scalar_field(sol, field, opt)
% 
% function [ ] = show_scalar_field(sol, field, opt)
%
% generic routine for plotting a grid with cell-centered values
%   sol     - structure with CONTACT results as returned by loadcase
%   field   - the actual array to be plotted, or a field name within sol
%   opt     - structure with options as defined by plot3d.
%

   % delegate work to separate functions depending on type of plot

   if (strcmp(opt.typplot, 'rw_rear'))
      show_rw_rear_view(sol, field, opt);
      if (isfield(sol, 'subs'))
         show_subs_field(sol, opt);
      end
   elseif (strcmp(opt.typplot, 'rw_side'))
      show_rw_side_view(sol, field, opt);
   else
      show_3d_field(sol, field, opt);
   end

end % function show_scalar_field

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ ] = show_subs_field(sol, opt)
% 
% function [ ] = show_subs_field(sol, field, opt)
%
% add plot of subsurface stresses to rw_rear plot

   if (0==1)
      % use slice closest to xc=0
      [~,ix] = min(abs(sol.subs.x));
      fld = squeeze( sol.subs.sigvm(ix,:,:) );
   else
      fld = squeeze( max(abs(sol.subs.sigvm), [], 1) );
   end

   % plot the potential contact area

   [YC, ZC] = meshgrid(sol.subs.y, sol.subs.z);
   if (strcmp(opt.rw_surfc,'prr'))
      [ ~, ypln, zpln ] = cntc_to_rail_coords(sol, [], YC, ZC);
   elseif (strcmp(opt.rw_surfc,'prw'))
      [ ~, ypln, zpln ] = cntc_to_wheel_coords(sol, [], YC, ZC);
   elseif (strcmp(opt.rw_surfc,'both'))
      [ ~, ypln, zpln ] = cntc_to_track_coords(sol, [], YC, ZC);
   end

   [C, h] = contourf( ypln, zpln, fld' );
   set(h, 'linewidth', 0.01);

end % function show_subs_field

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ ] = show_3d_field(sol, field, opt)

   % determine version information on the Matlab graphics system

   old_graphics = verLessThan('matlab', '8.4.0'); % 8.4 == R2014b
   new_graphics = ~old_graphics;

   % get field to be plotted

   if (ischar(field))
      eval(['tmp = sol.',field,';']);     
   else
      tmp = field;
   end

   % set element numbers ix, iy to be plotted: interval of [1:mx(my)]

   [ix_plot, iy_plot] = plot_ranges(sol, opt);

   % select the range to be plotted

   tmp    = tmp(iy_plot,ix_plot);
   eldiv  = sol.eldiv(iy_plot,ix_plot);
   xcornr = sol.xl + [ix_plot(1)-1, ix_plot]*sol.dx;
   ycornr = sol.yl + [iy_plot(1)-1, iy_plot]*sol.dy;
   if (~isempty(sol.x_offset)), xcornr = xcornr + sol.x_offset; end
   if (~isempty(sol.y_offset)), ycornr = ycornr + sol.y_offset; end
   xcornr = opt.xstretch * xcornr;

   % fill exterior area with default values

   if (~isempty(opt.exterval))
      ii = find(eldiv==0);
      tmp(ii) = opt.exterval;
   end

   if (isempty(opt.exterval) | ~isnan(opt.exterval))

      % expand one row/column for shading flat for cell-centered values

      tmp = [tmp, tmp(:,end)];
      tmp = [tmp; tmp(end,:)];

   elseif (exist('repelem','builtin'))

      % with flat shading, Matlab doesn't show faces where one of the corners
      % has value NaN --> duplicate values, shift-duplicate corner coordinates

      tmp = repelem(tmp,2,2);
      xcornr = [xcornr(1), repelem(xcornr(2:end-1),2), xcornr(end)];
      ycornr = [ycornr(1), repelem(ycornr(2:end-1),2), ycornr(end)];

   else

      % fall-back for Matlab versions < R2015a

      % expand one row/column for shading flat for cell-centered values

      tmp = [tmp, tmp(:,end)];
      tmp = [tmp; tmp(end,:)];

      %  perform dilatation:
      %  at exterior elements that lie to the right or above of interior
      %  elements, the tmp-value must be non-NaN

      sz=size(tmp);
      % rows 1:end-1 ==> rows 2:end
      [i,j]=find(eldiv(1:end-1,:)~=0 & eldiv(2:end,:)==0);
      ii0=sub2ind(sz,i,j); ii1=sub2ind(sz,i+1,j  ); tmp(ii1)=tmp(ii0);

      % row end ==> row end+1
      [i,j]=find(eldiv(end,:)~=0); i = i + size(eldiv,1) - 1;
      ii0=sub2ind(sz,i,j); ii1=sub2ind(sz,i+1,j  ); tmp(ii1)=tmp(ii0);
                           ii1=sub2ind(sz,i+1,j+1); tmp(ii1)=tmp(ii0);
      % columns 1:end-1 ==> columns 2:end
      [i,j]=find(eldiv(:,1:end-1)~=0 & eldiv(:,2:end)==0);
      ii0=sub2ind(sz,i,j); ii1=sub2ind(sz,i  ,j+1); tmp(ii1)=tmp(ii0);

      % column end ==> column end+1
      [i,j]=find(eldiv(:,end)~=0); j = j + size(eldiv,2) - 1;
      ii0=sub2ind(sz,i,j); ii1=sub2ind(sz,i+1,j  ); tmp(ii1)=tmp(ii0);
                           ii1=sub2ind(sz,i+1,j+1); tmp(ii1)=tmp(ii0);

      % diagonally iend,jend -> iend+1,jend+1
      [i,j]=find(eldiv(1:end-1,1:end-1)~=0 & eldiv(2:end,2:end)==0);
      ii0=sub2ind(sz,i,j); ii1=sub2ind(sz,i+1,j+1); tmp(ii1)=tmp(ii0);
   end

   % mark lower-left corner

   if (~strcmp(field,'ptarg'))
      minval   = min(min(tmp));
      tmp(1,1) = -0.2*max(max(abs(tmp)));
   end

   % expand x,y-coordinates to 2d arrays

   nx = length(xcornr);
   ny = length(ycornr);
   xcornr = ones(ny,1) * xcornr;
   ycornr = ycornr' * ones(1,nx);

   % set z-coordinates for plotting

   if (strcmp(opt.rw_surfc,'none'))
      zcornr = tmp;
   elseif (strcmp(opt.rw_surfc,'prr'))
      [ xcornr, ycornr, zcornr ] = cntc_to_rail_coords(sol, xcornr, ycornr);
      ii = find(isnan(tmp)); zcornr(ii) = NaN;
   elseif (strcmp(opt.rw_surfc,'prw'))
      [ xcornr, ycornr, zcornr ] = cntc_to_wheel_coords(sol, xcornr, ycornr);
      ii = find(isnan(tmp)); zcornr(ii) = NaN;
   elseif (strcmp(opt.rw_surfc,'both'))
      [ xcornr, ycornr, zcornr ] = cntc_to_track_coords(sol, xcornr, ycornr);
      ii = find(isnan(tmp)); zcornr(ii) = NaN;
   end

   % make plot, adjust colormap, z-range

   if (strcmp(opt.colormap,'none') | strcmp(opt.colormap,'black'))
      mesh( xcornr, ycornr, zcornr, tmp );
      view(opt.view);
      colormap([0 0 0]);
   else
      if (strcmp(opt.typplot,'contourf'))
         contourf( xcornr, ycornr, tmp, 'ShowText','on' );
      else
         % surf( xcornr, ycornr, tmp );
         surf( xcornr, ycornr, zcornr, tmp );
      end
      view(opt.view);
      if (strcmp(opt.rw_surfc,'none'))
         shading faceted;
      else
         shading flat;
      end
      colormap(opt.colormap);
      if (~isempty(opt.zrange))
         % use prescribed clim
         set(gca,'clim',opt.zrange);
      else
         % avoid spoiling of clim by special value at point (1,1):
         cl=get(gca,'clim');
         if (minval<cl(2))      % ensure that clim is increasing
            cl(1)=max(minval,cl(1));
            set(gca,'clim',cl);
         end
      end
   end

   % adjust axis

   nw_ax = [ min(min(xcornr)), max(max(xcornr)), min(min(ycornr)), max(max(ycornr)) ];
   if (opt.addplot==1)  % allow growing of plot region only
      pv_ax = axis;
      nw_ax = [ min(pv_ax(1),nw_ax(1)), max(pv_ax(2),nw_ax(2)) ...
                min(pv_ax(3),nw_ax(3)), max(pv_ax(4),nw_ax(4)) ];
   end

   if (strcmp(opt.rw_surfc,'none'))
      axis( nw_ax );
   end

   % optionally override automatic ticks, set ticks at cell-centers

   if (0==1)
      ixstep = max(1, round( (ix_plot(end)-ix_plot(1)+1) / 5 ));
      iystep = max(1, round( (iy_plot(end)-iy_plot(1)+1) / 5 ));
      x_tick = sol.xl + (0.5+[ix_plot(1):ixstep:ix_plot(end)])*sol.dx;
      y_tick = sol.yl + (0.5+[iy_plot(1):iystep:iy_plot(end)])*sol.dy;
      set(gca,'xtick', x_tick * opt.xstretch, ...
                               'xticklabel', num2str(x_tick', '%5.3f'));
      set(gca,'ytick', y_tick, 'yticklabel', num2str(y_tick', '%5.3f'));
   elseif (opt.xstretch~=1)
%     set tick-labels to unstretched values
      xt = get(gca,'xtick');
      xv = xt / opt.xstretch;
%     #digits before decimal point, negative: #zeros after decimal point
      lg = ceil(log10(max(abs(xv)))); 
      if (lg>=2)
         fmt = '%2.0f';
      elseif (lg==1 | lg==0)
         fmt = '%3.1f';
      elseif (lg==-1)
         fmt = '%4.2f';
      else
         fmt = '%5.3f';
      end
      set(gca, 'xtick', xt', 'xticklabel', num2str(xv', fmt));
   end

   % set grid, labels

   grid on
   set_xyzlabels(sol, opt);

   % draw colorbar, set label, if not done so before
   % note: no colorbar on mesh plot 

   if (opt.addplot<=0 & ...
       ~(strcmp(opt.colormap,'none') | strcmp(opt.colormap,'black')))

      h = colorbar;
      if (strcmp(opt.field,'pn'))
         ylabel(h, 'Pressure [N/mm^2]');
      elseif (strcmp(opt.field,'ptabs') | ...
          strcmp(opt.field,'ptabs+vec') | ...
          strcmp(opt.field,'px') | ...
          strcmp(opt.field,'py'))
         ylabel(h, 'Traction [N/mm^2]');
      elseif (strcmp(opt.field,'un') | strcmp(opt.field,'ux') | ...
              strcmp(opt.field,'uy') | strcmp(opt.field,'h'))
         ylabel(h, 'Displacement [mm]');
      elseif (strcmp(opt.field,'ptarg'))
         ylabel(h, 'Direction [deg]');
      elseif (strcmp(opt.field,'shft') | (sol.kincns.t_digit<=1 & ...
              (strcmp(opt.field,'sx')  | strcmp(opt.field,'sy'))))
         ylabel(h, 'Shift distance [mm]');
      elseif (strcmp(opt.field,'uplsx') | strcmp(opt.field,'uplsy') | strcmp(opt.field,'upls+vec'))
         ylabel(h, 'Plastic deformation [mm]');
      elseif (strcmp(opt.field,'taucrt'))
         ylabel(h, 'Tangential yield stress [N/mm^2]');
      elseif (strcmp(opt.field,'srel') | ...
              (sol.kincns.t_digit>=2 & (strcmp(opt.field,'sx')  | strcmp(opt.field,'sy'))) )
         ylabel(h, 'Relative slip velocity [-]');
      elseif (strcmp(opt.field,'sabs'))
         ylabel(h, 'Absolute slip velocity [mm/s]');
      elseif (strcmp(opt.field,'fricdens'))
         ylabel(h, 'Frictional power density [W/mm^2]');
      elseif (strcmp(opt.field,'temp1') | strcmp(opt.field,'temp2'))
         ylabel(h, 'Temperature (increase) [{}^\circ{}C]');
      % provisional support for user-added fields:
      elseif (strcmp(opt.field,'wx') | strcmp(opt.field,'wy'))
         ylabel(h, 'Rigid slip velocity [-]');
      elseif (strcmp(opt.field,'kyield'))
         ylabel(h, 'Yield strength K [N/mm^2]');
      elseif (strcmp(opt.field,'rcf'))
         ylabel(h, 'RCF index [-]');
      end
   end

end % function show_3d_field

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ ] = show_vec_field(sol, vx, vy, opt)
% 
% function [ ] = show_vec_field(sol, vx, vy, opt)
%
% generic routine for a quiver plot for cell-centered values
%   sol     - structure with CONTACT results as returned by loadcase
%   vx, vy  - the vector data to be plotted
%   opt     - structure with options as defined by plot3d.
%

   if (any(strcmp(opt.typplot,{'rw_rear','rw_side'})))
      return
   end

   % determine whether there are just vectors in the plot or magnitudes as well

   has_mgn = ~isempty(strfind(opt.field,'+vec'));

   % add to existing plot if magnitudes are already there

   if (has_mgn)
      hold on;
   end

   % set the vector color and line width, using defaults depending on
   % field to be plotted

   if (has_mgn)
      col='k';
      wid=2;
   else
      col='b';
      wid=1;
   end
   if (~isempty(opt.veccolor))
      col=opt.veccolor;
   end
   if (~isempty(opt.vecwidth))
      wid=opt.vecwidth;
   end

   % draw vectors in maximally about numvecx x numvecy points

   [ix_plot, iy_plot] = plot_ranges(sol, opt);
   ixstep = max(1, round( (ix_plot(end)-ix_plot(1)+1) / opt.numvecx ));
   iystep = max(1, round( (iy_plot(end)-iy_plot(1)+1) / opt.numvecy ));
   ix_vec = [ix_plot(1) : ixstep : ix_plot(end)];
   iy_vec = [iy_plot(1) : iystep : iy_plot(end)];
   x  = ones(size(iy_vec'))*sol.x(ix_vec);
   y  = sol.y(iy_vec)'*ones(size(ix_vec));
   if (~isempty(sol.x_offset)), x = x + sol.x_offset; end
   if (~isempty(sol.y_offset)), y = y + sol.y_offset; end
   x = x * opt.xstretch;

   % select components for actual vectors

   vx = vx(iy_vec,ix_vec);
   vy = vy(iy_vec,ix_vec);

   % eliminate vectors for points outside the contact area

   ix = find(sol.eldiv(iy_vec,ix_vec)==0); x(ix)=NaN;

   % set z-coordinates for plotting

   if (strcmp(opt.rw_surfc,'none'))
      zval = 1e6 * sign(opt.view(2)); z = zval * ones(size(x));
      vz = zeros(size(vx));
   elseif (strcmp(opt.rw_surfc,'prr'))
      [ x, y, z, deltar ] = cntc_to_rail_coords(sol, x, y);
      vz = sin(deltar) * vy;
      vy = cos(deltar) * vy;
   elseif (strcmp(opt.rw_surfc,'prw'))
      [ x, y, z, deltaw ] = cntc_to_wheel_coords(sol, x, y);
      vz = sin(deltaw) * vy;
      vy = cos(deltaw) * vy;
   elseif (strcmp(opt.rw_surfc,'both'))
      [ x, y, z, deltar ] = cntc_to_track_coords(sol, x, y);
      vz = sin(deltar) * vy;
      vy = cos(deltar) * vy;
   end

   % draw the actual vectors using quiver

   if (isempty(opt.vecscale) | opt.vecscale<=0)
      h = quiver3(x, y, z, vx, vy, vz, col);
   else
      % manual scaling of vectors, needed for addplot option:
      scl = opt.vecscale;
      h = quiver3(x, y, z, scl*vx, scl*vy, scl*vz, 0, col);
   end
   set(h,'linewidth',wid);

   % set huge z-data (+ or -) for the vectors when there are magnitudes as well

   % if (has_mgn)
   %    set(h,'zdata',sign(opt.view(2))*1e5*ones(size(get(h,'vdata'))));
   %    set(h,'wdata',    zeros(size(get(h,'vdata'))));
   % end

   % change axis range in 3D plots

   if (~is_2d_view(gca, opt.view) & strcmp(opt.rw_surfc,'none')) 
      v=axis; axis([v -1e-6 1e-6]); axis equal; 
      set(gca,'ztick',[]);
   end
   view(opt.view);

   % set labels

   set_xyzlabels(sol, opt);

end % function show_vec_field

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ ] = show_eldiv_contour(sol, opt)
% 
% function [ ] = show_eldiv_contour(sol, opt)
%
% routine for plotting outlines of contact and adhesion areas 
%   sol     - structure with CONTACT results as returned by loadcase
%   opt     - structure with options as defined by plot3d.
%

   if (sol.mx<=1 | sol.my<=1)
      if (0==1)
         disp('WARNING: option "eldiv_contour" is not available for 2D grids');
      end
      return
   end

   % set the line width, using default depending on field to be plotted
   if (~isempty(opt.eldivwid))
      wid = opt.eldivwid;
   elseif (strcmp(opt.field,'eldiv_contour'))
      wid = 1;
   else
      wid = 2;
   end

   [ix_plot, iy_plot] = plot_ranges(sol, opt);
   xcentr = sol.x(ix_plot); ycentr = sol.y(iy_plot);
   if (~isempty(sol.x_offset)), xcentr = xcentr + sol.x_offset; end
   if (~isempty(sol.y_offset)), ycentr = ycentr + sol.y_offset; end
   xcentr = xcentr * opt.xstretch;

   eldiv = derived_data(sol, opt, 'eldiv');
   eldiv = mod(eldiv,10); % diffcase: difference is stored in tens-digit
   eldiv = eldiv(iy_plot,ix_plot);

   % expand selection 'eldiv' and x/ycentr when there are elements in
   % contact in first/last rows/columns
   if (any(eldiv(1,:)))
      eldiv = [zeros(1,size(eldiv,2)) ; eldiv];
      ycentr = [ ycentr(1)-sol.dy, ycentr ];
   end
   if (any(eldiv(end,:)))
      eldiv = [eldiv ; zeros(1,size(eldiv,2))];
      ycentr = [ ycentr, ycentr(end)+sol.dy ];
   end
   if (any(eldiv(:,1)))
      eldiv = [zeros(size(eldiv,1),1) , eldiv];
      xcentr = [ xcentr(1)-sol.dx, xcentr ];
   end
   if (any(eldiv(:,end)))
      eldiv = [eldiv , zeros(size(eldiv,1),1)];
      xcentr = [ xcentr, xcentr(end)+sol.dx ];
   end

   % plot double line around adhesion area
   l1 = []; l2 = []; l3 = [];
   mask = 1.0*(eldiv==1);
   col  = opt.eldivcol(1,:);
   zval = 1e5;
   if (nnz(mask)>0 & ~strcmp(opt.field,'pn'))
      ContourLvl = [0.45, 0.55]';
      c = contourc(xcentr, ycentr, mask, ContourLvl);
      nx = size(c,2); ix = 1;
      while(ix < nx);
         npnt = c(2,ix); c(:,ix) = NaN; ix = ix + 1 + npnt;
      end

      % set xyz-coordinates for plotting

      if (strcmp(opt.rw_surfc,'none'))
         xcrv = c(1,:); ycrv = c(2,:); zcrv = zval*[1;-1]*ones(1,nx);
      elseif (strcmp(opt.rw_surfc,'prr'))
         [ xcrv, ycrv, zcrv ] = cntc_to_rail_coords(sol, c(1,:), c(2,:));
      elseif (strcmp(opt.rw_surfc,'prw'))
         [ xcrv, ycrv, zcrv ] = cntc_to_wheel_coords(sol, c(1,:), c(2,:));
      elseif (strcmp(opt.rw_surfc,'both'))
         [ xcrv, ycrv, zcrv ] = cntc_to_track_coords(sol, c(1,:), c(2,:));
      end

      l1 = plot3(xcrv, ycrv, zcrv, 'color',col, 'linewidth',wid);
   end

   % plot double line around plasticity area
   mask = 1.0*(eldiv==3);
   col  = opt.eldivcol(3,:);
   zval = 1e5;
   if (nnz(mask)>0 & ~strcmp(opt.field,'pn'))
      ContourLvl = [0.45, 0.55]';
      c = contourc(xcentr, ycentr, mask, ContourLvl);
      nx = size(c,2); ix = 1;
      while(ix < nx);
         npnt = c(2,ix); c(:,ix) = NaN; ix = ix + 1 + npnt;
      end

      % set xyz-coordinates for plotting

      if (strcmp(opt.rw_surfc,'none'))
         xcrv = c(1,:); ycrv = c(2,:); zcrv = zval*[1;-1]*ones(1,nx);
      elseif (strcmp(opt.rw_surfc,'prr'))
         [ xcrv, ycrv, zcrv ] = cntc_to_rail_coords(sol, c(1,:), c(2,:));
      elseif (strcmp(opt.rw_surfc,'prw'))
         [ xcrv, ycrv, zcrv ] = cntc_to_wheel_coords(sol, c(1,:), c(2,:));
      elseif (strcmp(opt.rw_surfc,'both'))
         [ xcrv, ycrv, zcrv ] = cntc_to_track_coords(sol, c(1,:), c(2,:));
      end

      l3 = plot3(xcrv, ycrv, zcrv, 'color',col, 'linewidth',wid);
   end

   % plot single line around contact area
   mask = 1.0*(eldiv>=1);
   col = opt.eldivcol(2,:);
   if (nnz(mask)>0)
      ContourLvl = [0.5, 9.99]';
      c = contourc(xcentr, ycentr, mask, ContourLvl);
      nx = size(c,2); ix = 1;
      while(ix < nx);
         npnt = c(2,ix); c(:,ix) = NaN; ix = ix + 1 + npnt;
      end

      % set xyz-coordinates for plotting

      if (strcmp(opt.rw_surfc,'none'))
         xcrv = c(1,:); ycrv = c(2,:); zcrv = zval*[1;-1]*ones(1,nx);
      elseif (strcmp(opt.rw_surfc,'prr'))
         [ xcrv, ycrv, zcrv ] = cntc_to_rail_coords(sol, c(1,:), c(2,:));
      elseif (strcmp(opt.rw_surfc,'prw'))
         [ xcrv, ycrv, zcrv ] = cntc_to_wheel_coords(sol, c(1,:), c(2,:));
      elseif (strcmp(opt.rw_surfc,'both'))
         [ xcrv, ycrv, zcrv ] = cntc_to_track_coords(sol, c(1,:), c(2,:));
      end

      l2 = plot3(xcrv, ycrv, zcrv, 'color',col, 'linewidth',wid);
   end

   % adapt axis in order to fit the plotted range of the contact area
   if (strcmp(opt.rw_surfc,'none'))
      xcornr = sol.xl + [ ix_plot(1)-1, ix_plot ]*sol.dx;
      ycornr = sol.yl + [ iy_plot(1)-1, iy_plot ]*sol.dy;
      if (~isempty(sol.x_offset)), xcornr = xcornr + sol.x_offset; end
      if (~isempty(sol.y_offset)), ycornr = ycornr + sol.y_offset; end
      nw_ax = [ min(xcornr), max(xcornr), min(ycornr), max(ycornr) ];
      if (opt.addplot==1)  % allow growing of plot region only
         pv_ax = axis;
         nw_ax = [ min(pv_ax(1),nw_ax(1)), max(pv_ax(2),nw_ax(2)) ...
                   min(pv_ax(3),nw_ax(3)), max(pv_ax(4),nw_ax(4)) ];
      end
      axis( nw_ax );

      % in 3D plots, concentrate on plane Oxy using tiny z-range
      %    except when the eldiv is added to another 3D plot
      if (~is_2d_view(gca, opt.view) & strcmp(opt.field,'eldiv_contour'))
         v=axis; axis([v -1e-6 1e-6]); axis equal; 
         set(gca,'ztick',[]);
      end
   end

   % set view and labels
   view(opt.view);
   set_xyzlabels(sol, opt);

end % function show_eldiv_contour

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function show_rw_rear_view(sol, field, opt)

% show_rw_rear_view: plot w/r profiles, contact plane & contact results

   % plot the potential contact area

   yc = sol.y; zc = zeros(1,sol.my);
   if (strcmp(opt.rw_surfc,'prr'))
      [ ~, ypln, zpln ] = cntc_to_rail_coords(sol, [], yc, zc);
   elseif (strcmp(opt.rw_surfc,'prw'))
      [ ~, ypln, zpln ] = cntc_to_wheel_coords(sol, [], yc, zc);
   elseif (strcmp(opt.rw_surfc,'both'))
      [ ~, ypln, zpln ] = cntc_to_track_coords(sol, [], yc, zc);
   end

   plot(ypln, zpln, '.-', 'markersize',12, 'color',matlab_color(5));

   % get field to be plotted, replace values in exterior elements
   if (ischar(field))
      eval(['val = sol.',field,';']);     
   else
      val = field;
   end
   if (~isempty(opt.exterval))
      eldiv = derived_data(sol, opt, 'eldiv');
      adjac = eldiv;
      adjac(1:end-1,:) = max(adjac(1:end-1,:), eldiv(2:end  ,:));
      adjac(2:end  ,:) = max(adjac(2:end  ,:), eldiv(1:end-1,:));
      adjac(:,1:end-1) = max(adjac(:,1:end-1), eldiv(:,2:end  ));
      adjac(:,2:end  ) = max(adjac(:,2:end  ), eldiv(:,1:end-1));
      ii = find(adjac==0);
      val(ii) = opt.exterval;
   end

   % convert field to 2d view, aggregating values over columns
   % TODO: provide max, absmax, maxabs, sum
   val   = max(val, [], 2)';
   scale = get_vecscale( sol, opt, val );

   yc = sol.y; zc = -scale*val;
   if (strcmp(opt.rw_surfc,'prr'))
      [ ~, yval, zval ] = cntc_to_rail_coords(sol, [], yc, zc);
   elseif (strcmp(opt.rw_surfc,'prw'))
      [ ~, yval, zval ] = cntc_to_wheel_coords(sol, [], yc, zc);
   elseif (strcmp(opt.rw_surfc,'both'))
      [ ~, yval, zval ] = cntc_to_track_coords(sol, [], yc, zc);
   end

   col = matlab_color(6);
   if (~isempty(opt.veccolor))
      col = opt.veccolor;
   end
   plot( yval, zval, 'linewidth',1, 'color',col );
   plot( [ypln; yval], [zpln; zval], 'linewidth',1, 'color',col );

   % plot contact reference point, origin of contact local coordinate system

   axlen = 0.1 * (sol.y(end) - sol.y(1));
   yc = [0,0,1]*axlen; zc = [1,0,0]*axlen;
   if (strcmp(opt.rw_surfc,'prr'))
      [ ~, yval, zval, delta ] = cntc_to_rail_coords(sol, [], yc, zc);
   elseif (strcmp(opt.rw_surfc,'prw'))
      [ ~, yval, zval, delta ] = cntc_to_wheel_coords(sol, [], yc, zc);
   elseif (strcmp(opt.rw_surfc,'both'))
      [ ~, yval, zval, delta ] = cntc_to_track_coords(sol, [], yc, zc);
   end

   if (exist('plot_cm','file'))
      l  = plot_cm([yval(2),zval(2)], 0.3, axlen, delta*180/pi, 'k');
      set(l, 'Tag','CM');
   else
      l1 = plot(yval, zval, 'k');
      l2 = plot(yval(2), zval(2), 'k.', 'markersize',18);
      set([l1,l2], 'Tag','CM');
   end

   set_yzlabels(sol, opt);
   grid on;

end % function show_rw_rear_view

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function show_rw_side_view(sol, field, opt)

% show_rw_side_view: plot w/r profiles, contact plane & contact results

   % plot the potential contact area

   xc = sol.x; zc = zeros(1,sol.mx);
   if (strcmp(opt.rw_surfc,'prr'))
      [ xpln, ~, zpln ] = cntc_to_rail_coords(sol, xc, [], zc);
   elseif (strcmp(opt.rw_surfc,'prw'))
      [ xpln, ~, zpln ] = cntc_to_wheel_coords(sol, xc, [], zc);
   elseif (strcmp(opt.rw_surfc,'both'))
      [ xpln, ~, zpln ] = cntc_to_track_coords(sol, xc, [], zc);
   end

   plot(xpln, zpln, '.-', 'markersize',12, 'color',matlab_color(5));

   % get field to be plotted, replace values in exterior elements
   if (ischar(field))
      eval(['val = sol.',field,';']);     
   else
      val = field;
   end
   if (~isempty(opt.exterval))
      eldiv = derived_data(sol, opt, 'eldiv');
      adjac = eldiv;
      adjac(1:end-1,:) = max(adjac(1:end-1,:), eldiv(2:end  ,:));
      adjac(2:end  ,:) = max(adjac(2:end  ,:), eldiv(1:end-1,:));
      adjac(:,1:end-1) = max(adjac(:,1:end-1), eldiv(:,2:end  ));
      adjac(:,2:end  ) = max(adjac(:,2:end  ), eldiv(:,1:end-1));
      ii = find(adjac==0);
      val(ii) = opt.exterval;
   end

   % convert field to 2d view, aggregating values over rows
   % TODO: provide max, absmax, maxabs, sum

   val   = max(val, [], 1);
   scale = get_vecscale( sol, opt, val );

   xc = sol.x; zc = -scale*val;
   if (strcmp(opt.rw_surfc,'prr'))
      [ xval, ~, zval ] = cntc_to_rail_coords(sol, xc, [], zc);
   elseif (strcmp(opt.rw_surfc,'prw'))
      [ xval, ~, zval ] = cntc_to_wheel_coords(sol, xc, [], zc);
   elseif (strcmp(opt.rw_surfc,'both'))
      [ xval, ~, zval ] = cntc_to_track_coords(sol, xc, [], zc);
   end

   col = matlab_color(6);
   if (~isempty(opt.veccolor))
      col = opt.veccolor;
   end

   plot( xval, zval, 'linewidth',1, 'color',col );
   plot( [xpln; xval], [zpln; zval], 'linewidth',1, 'color',col );

   % plot contact reference point, origin of contact local coordinate system

   v = axis;
   axlen = 0.03 * (v(2) - v(1));
   xc = [0,0,1]*axlen; zc = [1,0,0]*axlen;
   if (strcmp(opt.rw_surfc,'prr'))
      [ xval, ~, zval ] = cntc_to_rail_coords(sol, xc, [], zc);
   elseif (strcmp(opt.rw_surfc,'prw'))
      [ xval, ~, zval ] = cntc_to_wheel_coords(sol, xc, [], zc);
   elseif (strcmp(opt.rw_surfc,'both'))
      [ xval, ~, zval ] = cntc_to_track_coords(sol, xc, [], zc);
   end

   if (exist('plot_cm','file'))
      roty = 0; % hack for 'grade' normal [nx,0,nz] instead of [0,0,1]
      if (isfield(sol.meta, 'roty')), roty = sol.meta.roty*180/pi; end
      l = plot_cm([xval(2),zval(2)], 0.15*axlen, axlen, -roty, 'k');
      set(l, 'Tag','CM');
   else
      l1 = plot(xval, zval, 'k');
      l2 = plot(xval(2), zval(2), 'k.', 'markersize',18);
      set([l1,l2], 'Tag','CM');
   end

   set_xyzlabels(sol, opt);
   grid on;

end % function show_rw_side_view

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ is2d ] = is_2d_view(ax, view);

% return 1 if the current plot uses 2D view and 0 in case of a 3D view
% 2D view: myopt.view = [ *, +/-90 ];  axis returns vector of 4 elements.

   is2d = (abs(90 - abs(view(2))) <= 1e-3);

end % function is_2d_view

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ ] = set_yzlabels(sol, opt)

   if (strcmp(opt.rw_surfc,'prr'))
      xlabel ('y_r [mm]');
      ylabel ('z_r [mm]');
   elseif (strcmp(opt.rw_surfc,'prw'))
      xlabel ('y_w [mm]');
      ylabel ('z_w [mm]');
   elseif (strcmp(opt.rw_surfc,'both'))
      xlabel ('y_{tr} [mm]');
      ylabel ('z_{tr} [mm]');
   end

end % function set_yzlabels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ ] = set_xyzlabels(sol, opt)

   if (strcmp(opt.typplot,'rw_rear'))

      if (strcmp(opt.rw_surfc,'prr'))
         xlabel ('y_r [mm]');
         ylabel ('z_r [mm]');
      elseif (strcmp(opt.rw_surfc,'prw'))
         xlabel ('y_w [mm]');
         ylabel ('z_w [mm]');
      elseif (strcmp(opt.rw_surfc,'both'))
         xlabel ('y_{tr} [mm]');
         ylabel ('z_{tr} [mm]');
      end

   elseif (strcmp(opt.typplot,'rw_side'))

      if (strcmp(opt.rw_surfc,'prr'))
         xlabel ('x_r [mm]');
         ylabel ('z_r [mm]');
      elseif (strcmp(opt.rw_surfc,'prw'))
         xlabel ('x_w [mm]');
         ylabel ('z_w [mm]');
      elseif (strcmp(opt.rw_surfc,'both'))
         xlabel ('x_{tr} [mm]');
         ylabel ('z_{tr} [mm]');
      end

   else

      if (strcmp(opt.rw_surfc,'prr'))
         xlabel ('x_r [mm]');
         ylabel ('y_r [mm]');
      elseif (strcmp(opt.rw_surfc,'prw'))
         xlabel ('x_w [mm]');
         ylabel ('y_w [mm]');
      elseif (strcmp(opt.rw_surfc,'both'))
         xlabel ('x_{tr} [mm]');
         ylabel ('y_{tr} [mm]');
      else
         xlabel ('x_c [mm]');
         if (sol.d_digit==4)
            ylabel ('s_c [mm]');
         else
            ylabel ('y_c [mm]');
         end
%     else
%        xlabel ('X-coordinate [mm]');
%        ylabel ('Y-coordinate [mm]');
      end
   end

end % function set_xyzlabels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ix_plot, iy_plot] = plot_ranges(sol, opt);

% set element numbers ix, iy to be plotted: interval of [1:mx(my)]

if (isempty(opt.ixrange) | ...
    (ischar(opt.ixrange) & strcmp(opt.ixrange,'all')) )
   ix_plot = [1, sol.mx];
elseif (ischar(opt.ixrange) & strcmp(opt.ixrange,'tight'))
   num_inter = sum(sol.eldiv,1);
   ix_plot = find(num_inter);
   ix_plot = [ix_plot(1)-1 : ix_plot(end)+1 ];
elseif (ischar(opt.ixrange) & strcmp(opt.ixrange,'auto'))
   num_inter = sum(sol.eldiv,1);
   ix_plot = find(num_inter);
   if (isempty(ix_plot))
      ix_plot = [1, sol.mx];
   else
      ix_plot = [ix_plot(1)-3 : ix_plot(end)+3 ];
   end
else
   ix_plot = opt.ixrange;
end
if (isempty(opt.iyrange) | ...
    (ischar(opt.iyrange) & strcmp(opt.iyrange,'all')) )
   iy_plot = [1, sol.my];
elseif (ischar(opt.iyrange) & strcmp(opt.iyrange,'tight'))
   num_inter = sum(sol.eldiv,2);
   iy_plot = find(num_inter);
   iy_plot = [iy_plot(1)-1 : iy_plot(end)+1 ];
elseif (ischar(opt.iyrange) & strcmp(opt.iyrange,'auto'))
   num_inter = sum(sol.eldiv,2);
   iy_plot = find(num_inter);
   if (isempty(iy_plot))
      iy_plot = [1, sol.my];
   else
      iy_plot = [iy_plot(1)-3 : iy_plot(end)+3 ];
   end
else
   iy_plot = opt.iyrange;
end
ix_plot = max(1,min(ix_plot)) : min(sol.mx,max(ix_plot));
iy_plot = max(1,min(iy_plot)) : min(sol.my,max(iy_plot));

end % function plot_ranges

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ f1, f2, f3, f4, f5, f6 ] = ...
                         derived_data(sol, opt, n1, n2, n3, n4, n5, n6)

px = sol.px; py = sol.py; eldiv = sol.eldiv;

nfield = nargin - 2;
for ifield = 1 : nfield
   % get name of requested derived field, by copying n1, n2, .. or n6
   eval(sprintf('nam = n%d;',ifield));

   % compute the derived field
   if (strcmp(nam,'eldiv'))
      field = eldiv;
   elseif (strcmp(nam,'px'))
      field = px;
   elseif (strcmp(nam,'py'))
      field = py;
   elseif (strcmp(nam,'pt'))
      field = sqrt(px.^2 + py.^2);
   elseif (strcmp(nam,'ptarg'))
      field = 180/pi*atan2(py, px);   % range [-180,180]
   elseif (strcmp(nam,'pxdir'))
      field = sol.px ./ max(1d-9,sqrt(sol.px.^2 + sol.py.^2));
   elseif (strcmp(nam,'sabs'))
      field = sol.kincns.veloc * sol.srel;
   else
      disp(sprintf('Error: unknown derived field %s.',nam));
   end

   % store derived field, by copying to f1, f2, .. or f6
   eval(sprintf('f%d = field;',ifield));
end

end % function derived_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ scale ] = get_vecscale( sol, opt, val )

   if (~isempty(opt.vecscale))
      scale = opt.vecscale;
   else
      scale = 10 / max(abs(val));
   end

end % function get_vecscale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ is_left ] = is_left_side(sol)

if (length(sol)>1)
   is_left = (sol(1).config==0 | sol(1).config==4);
else
   is_left = (sol.config==0 | sol.config==4);
end

end % function is_left_side

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ prr ] = get_profile_slice( slcs, s_i, make_plot )

if (nargin<3)
   make_plot = 0;
end
i0 = find(slcs.s <  s_i, 1, 'last');
i1 = find(slcs.s >= s_i, 1, 'first');
% disp([s_i, i0, i1])
if (i0 >= slcs.nslc)
   yi = slcs.ysurf(end,:); zi = slcs.zsurf(end,:);
   % disp(sprintf('s_i = %5.1f > s(end) = %5.1f, using slice %d', s_i, slcs.s(end), i1));
elseif (i1 <= 1)
   yi = slcs.ysurf(1,:); zi = slcs.zsurf(1,:);
   % disp(sprintf('s_i = %5.1f < s(1) = %5.1f, using slice %d', s_i, slcs.s(1), i0));
else
   fac0 = (slcs.s(i1) - s_i) / (slcs.s(i1) - slcs.s(i0));
   fac1 = (s_i - slcs.s(i0)) / (slcs.s(i1) - slcs.s(i0));
   if (fac0>0.99 & nnz(slcs.mask_j(i0,:))>nnz(slcs.mask_j(i1,:)))
      % disp(sprintf('Using longer slice i0=%d', i0));
      yi = slcs.ysurf(i0,:); zi = slcs.zsurf(i0,:);
   elseif (fac1>0.99 & nnz(slcs.mask_j(i1,:))>nnz(slcs.mask_j(i0,:)))
      % disp(sprintf('Using longer slice i1=%d', i1));
      yi = slcs.ysurf(i1,:); zi = slcs.zsurf(i1,:);
   else
      % disp(sprintf('Using %5.3f * slice %d + %5.3f * slice %d',fac0, i0, fac1, i1));
      yi = fac0 * slcs.ysurf(i0,:) + fac1 * slcs.ysurf(i1,:);
      zi = fac0 * slcs.zsurf(i0,:) + fac1 * slcs.zsurf(i1,:);
   end

   if (make_plot)
      tmpfig = gcf;
      figure(make_plot); clf; hold on;
      plot([slcs.ysurf([i0:i1],:);yi]', [slcs.zsurf([i0:i1],:);zi]');
      set(gca,'ydir','reverse'); grid on; axis equal;
      legend(sprintf('slice %d, s=%6.2f',i0,slcs.s(i0)), sprintf('slice %d, s=%6.2f',i1,slcs.s(i1)), ...
                sprintf('interpolated, s=%6.2f',s_i), 'location','southeast')
      figure(tmpfig);
   end
end

% set arc-length s, accounting for NaN's in 'missing parts' 
% s  = slcs.vj; ix = find(isnan(yi)); s(ix) = NaN;

prr = struct('ProfileY',yi', 'ProfileZ',zi', 'ProfileS',slcs.vj);

end % function get_profile_slice

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ xtr, ytr, ztr, delta_tr ] = cntc_to_track_coords(sol, xc, yc, zc);

% transform contact [xc,yc/s]-coordinates to track [xtr,ytr,ztr] coordinates

   if (isempty(xc))
      xc = zeros(size(yc));
   end
   if (isempty(yc))
      yc = zeros(size(zc));
   end
   if (nargin<4 | isempty(zc))
      zc = zeros(size(yc));
   end

   % reshape inputs to row vectors

   m = size(xc,1); n = size(xc,2);
   xc = reshape(xc, 1, m*n);
   yc = reshape(yc, 1, m*n);
   zc = reshape(zc, 1, m*n);
   coords = [xc; yc; zc];

   % form transformation matrices and vectors

   Rref = rotx(sol.meta.deltcp_r*180/pi);
   if (isfield(sol.meta, 'roty')) % hack for 'grade' normal [nx,0,nz] instead of [0,0,1]
      Rref = Rref * roty(sol.meta.roty*180/pi);
   end
   oref = [sol.meta.xcp_r; sol.meta.ycp_r; sol.meta.zcp_r];

   R_r  = rotx(sol.meta.roll_r*180/pi);
   o_r  = [0; sol.meta.y_r; sol.meta.z_r];

   if (sol.d_digit~=4)

      % change coordinates from contact to rail to track coordinates

      coords = o_r*ones(1,m*n) + R_r * oref*ones(1,m*n) + R_r * Rref * coords;

   else

      % mirror ProfileY for left-side w/r pairs

      prf_s = sol.prr.ProfileS; prf_y = sol.prr.ProfileY; prf_z = sol.prr.ProfileZ;
      if (is_left_side(sol))
         prf_s = -flipud(prf_s); prf_y = -flipud(prf_y); prf_z =  flipud(prf_z);
      end

      % find oref in the rail profile

      dst = sqrt( (prf_y-sol.meta.ycp_r).^2 + (prf_z-sol.meta.zcp_r).^2 );
      [~,ix] = min(dst);
      rg = [max(1,ix-5) : min(length(prf_y),ix+5)];

      sref = interp1(prf_y(rg), prf_s(rg), sol.meta.ycp_r);
      if (isnan(sref))
         disp(sprintf('ERROR: prr seems different from prr used in simulation (needs mirroring?).'))
         sref = prf_s(ix);
      end

      % find grid points in the rail profile

      sc = yc;
      sr = sref + sc;
      yr_at_s = interp1(prf_s, prf_y, sr);
      zr_at_s = interp1(prf_s, prf_z, sr);

      % determine n-vectors on the rail profile

      ds = 0.1;
      dy = (interp1(prf_s, prf_y, sr+ds) - yr_at_s) / ds;
      dz = (interp1(prf_s, prf_z, sr+ds) - zr_at_s) / ds;
      nvec = [zeros(1,m*n); -dz; dy] ./ ([1;1;1] * sqrt(dy.^2+dz.^2));

      % hack: slight raise to put the data above the surfaces
      % zc = zc - 0.3;

      coords = [xc+oref(1); yr_at_s; zr_at_s] + nvec .* ([1;1;1] * zc);

      % change coordinates from rail to track coordinates

      coords = o_r*ones(1,m*n) + R_r * coords;

   end

   % extract components & reshape to the original array sizes

   xtr = reshape(coords(1,:), m, n);
   ytr = reshape(coords(2,:), m, n);
   ztr = reshape(coords(3,:), m, n);
        % TODO: correction for \Delta\phi_rr
        % delta_tr = sol.meta.deltcp_r + sol.meta.roll_r;
   delta_tr = sol.meta.deltcp_r + sol.meta.roll_r;

end % function cntc_to_track_coords

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ xr, yr, zr, deltar ] = cntc_to_rail_coords(sol, xc, yc, zc);

% transform contact [xc,yc/s]-coordinates to rail [xr,yr,zr] coordinates

   if (nargin<3), yc = []; end
   if (nargin<4), zc = []; end
   sz = max( [ size(xc) ; size(yc) ; size(zc) ] );
   if (isempty(xc)), xc = zeros(sz); end
   if (isempty(yc)), yc = zeros(sz); end
   if (isempty(zc)), zc = zeros(sz); end

   % reshape inputs to row vectors

   m = size(xc,1); n = size(xc,2);
   xc = reshape(xc, 1, m*n);
   yc = reshape(yc, 1, m*n);
   zc = reshape(zc, 1, m*n);
   coords = [xc; yc; zc];

   % form transformation matrices and vectors

   Rref = rotx(sol.meta.deltcp_r*180/pi);
   if (isfield(sol.meta, 'roty')) % hack for 'grade' normal [nx,0,nz] instead of [0,0,1]
      Rref = Rref * roty(sol.meta.roty*180/pi);
   end

   % for variable rails, x_r-coordinates are defined by the slices-file rather than aligned with track x_tr

   has_slcs = isfield(sol,'slcs');
   if (has_slcs)
      oref = [sol.meta.s_ws+sol.meta.xcp_r; sol.meta.ycp_r; sol.meta.zcp_r];
   else
      oref = [              sol.meta.xcp_r; sol.meta.ycp_r; sol.meta.zcp_r];
   end

   if (sol.d_digit~=4)

      % change coordinates from contact to rail coordinates

      coords = oref*ones(1,m*n) + Rref * coords;

   else

      if (has_slcs), disp('Rail coords not supported for conformal on variable profile'); end

      % mirror ProfileY for left-side w/r pairs

      prf_s = sol.prr.ProfileS; prf_y = sol.prr.ProfileY; prf_z = sol.prr.ProfileZ;
      if (is_left_side(sol))
         prf_s = -flipud(prf_s); prf_y = -flipud(prf_y); prf_z =  flipud(prf_z);
      end

      % find oref in the rail profile

      dst = sqrt( (prf_y-sol.meta.ycp_r).^2 + (prf_z-sol.meta.zcp_r).^2 );
      [~,ix] = min(dst);
      rg = [max(1,ix-5), min(length(prf_y),ix+5)];

      sref = interp1(prf_y(rg), prf_s(rg), sol.meta.ycp_r);

      % find grid points in the rail profile

      sc = yc;
      sr = sref + sc;
      yr_at_s = interp1(prf_s, prf_y, sr);
      zr_at_s = interp1(prf_s, prf_z, sr);

      % determine n-vectors on the rail profile

      ds = 0.1;
      dy = (interp1(prf_s, prf_y, sr+ds) - yr_at_s) / ds;
      dz = (interp1(prf_s, prf_z, sr+ds) - zr_at_s) / ds;
      nvec = [zeros(1,m*n); -dz; dy] ./ ([1;1;1] * sqrt(dy.^2+dz.^2));

      coords = [xc; yr_at_s; zr_at_s] + nvec .* ([1;1;1] * zc);

   end

   % extract components & reshape to the original array sizes

   xr = reshape(coords(1,:), m, n);
   yr = reshape(coords(2,:), m, n);
   zr = reshape(coords(3,:), m, n);
   deltar = sol.meta.deltcp_r;

end % function cntc_to_rail_coords

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ xw, yw, zw, deltaw ] = cntc_to_wheel_coords(sol, xc, yc, zc);

% transform contact [xc,yc/s]-coordinates to wheel [xw,yw,zw] coordinates

   if (nargin<3), yc = []; end
   if (nargin<4), zc = []; end
   sz = max( [ size(xc) ; size(yc) ; size(zc) ] );
   if (isempty(xc)), xc = zeros(sz); end
   if (isempty(yc)), yc = zeros(sz); end
   if (isempty(zc)), zc = zeros(sz); end

   % reshape inputs to row vectors

   m = size(xc,1); n = size(xc,2);
   xc = reshape(xc, 1, m*n);
   yc = reshape(yc, 1, m*n);
   zc = reshape(zc, 1, m*n);
   coords = [xc; yc; zc];

   % form transformation matrices and vectors

   Rref = rotx(sol.meta.deltcp_r*180/pi);
   oref = [sol.meta.xcp_r; sol.meta.ycp_r; sol.meta.zcp_r];

   R_r  = rotx(sol.meta.roll_r*180/pi);
   o_r  = [0; sol.meta.y_r; sol.meta.z_r];

   R_w  = rotx(sol.meta.roll_w*180/pi) * rotz(sol.meta.yaw_w*180/pi);
   o_w  = [sol.meta.x_w; sol.meta.y_w; sol.meta.z_w];

   if (sol.d_digit~=4)

      % change coordinates from contact to rail to track to wheel coordinates

      xtr = o_r*ones(1,m*n) + R_r * oref*ones(1,m*n) + R_r * Rref * coords;
      coords = R_w' * ( xtr - o_w*ones(1,m*n) );

   else

      % mirror ProfileY for left-side w/r pairs

      prf_s = sol.prw.ProfileS; prf_y = sol.prw.ProfileY; prf_z = sol.prw.ProfileZ;
      if (is_left_side(sol))
         prf_s = -flipud(prf_s); prf_y = -flipud(prf_y); prf_z =  flipud(prf_z);
      end

      % change oref from rail to track to wheel coordinates

      wref = R_w' * ( o_r + R_r * oref - o_w );

      % find wref.y in the wheel profile

      dst = sqrt( (prf_y-wref(2)).^2 + (prf_z-wref(3)).^2 );
      [~,ix] = min(dst);
      rg = [max(1,ix-5), min(length(prf_y),ix+5)];

      sref = interp1(prf_y(rg), prf_s(rg), wref(2));

      % find grid points in the wheel profile

      sc = yc;
      sw = sref - sc;    % sw increasing to the left
      yw_at_s = interp1(prf_s, prf_y, sw);
      zw_at_s = interp1(prf_s, prf_z, sw);

      % determine n-vectors on the wheel profile

      ds = 0.1;
      dy = (interp1(prf_s, prf_y, sw-ds) - yw_at_s) / ds;
      dz = (interp1(prf_s, prf_z, sw-ds) - zw_at_s) / ds;
      nvec = [zeros(1,m*n); -dz; dy] ./ ([1;1;1] * sqrt(dy.^2+dz.^2));

      coords = [xc; yw_at_s; zw_at_s] + nvec .* ([1;1;1] * zc);
   end

   % extract coordinates, add wheel curvature

   fac_r = 1 / (2*sol.meta.rnom_whl);
   xw = coords(1,:);
   yw = coords(2,:);
   zw = coords(3,:) - fac_r * xw.^2;

   % reshape to the original array sizes

   xw = reshape(xw, m, n);
   yw = reshape(yw, m, n);
   zw = reshape(zw, m, n);

   deltaw = sol.meta.deltcp_r;

end % function cntc_to_wheel_coords

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ xtr, ytr, ztr ] = rail_to_track_coords(sol, xr, yr, zr);

% transform rail [xr,yr,zr]-coordinates to track [xtr,ytr,ztr] coordinates

   % reshape inputs to row vectors

   m = size(xr,1); n = size(xr,2);
   xr = reshape(xr, 1, m*n);
   yr = reshape(yr, 1, m*n);
   zr = reshape(zr, 1, m*n);
   coords = [xr; yr; zr];

   % form transformation matrices and vectors

   R_r  = rotx(sol.meta.roll_r*180/pi);

   % for variable rails, x_r-coordinates are defined by the slices-file rather than aligned with track x_tr

   has_slcs = isfield(sol,'slcs');
   if (has_slcs)
      o_r  = [-sol.meta.s_ws; sol.meta.y_r; sol.meta.z_r];
   else
      o_r  = [          0   ; sol.meta.y_r; sol.meta.z_r];
   end

   % change coordinates from rail to track coordinates

   coords = o_r * ones(1,m*n) + R_r * coords;

   % extract & reshape to the original array sizes

   xtr = reshape(coords(1,:), m, n);
   ytr = reshape(coords(2,:), m, n);
   ztr = reshape(coords(3,:), m, n);

end % function rail_to_track_coords

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ xtr, ytr, ztr ] = wheel_to_track_coords(sol, xw, yw, zw);

% transform wheel [xw,yw,zw]-coordinates to track [xtr,ytr,ztr] coordinates

   % reshape inputs to row vectors

   m = size(xw,1); n = size(xw,2);
   xw = reshape(xw, 1, m*n);
   yw = reshape(yw, 1, m*n);
   zw = reshape(zw, 1, m*n);
   coords = [xw; yw; zw];

   % form transformation matrices and vectors from rw/ws to tr coordinates

   R_w_tr = rotx(sol.meta.roll_w*180/pi) * rotz(sol.meta.yaw_w*180/pi);
   o_w_tr = [sol.meta.x_w; sol.meta.y_w; sol.meta.z_w];

   % change coordinates from wheel to track coordinates

   coords = o_w_tr * ones(1,m*n) + R_w_tr * coords;

   % extract & reshape to the original array sizes

   xtr = reshape(coords(1,:), m, n);
   ytr = reshape(coords(2,:), m, n);
   ztr = reshape(coords(3,:), m, n);

end % function wheel_to_track_coords

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ xr, yr, zr ] = wheel_to_rail_coords(sol, xw, yw, zw);

% transform wheel [xw,yw,zw]-coordinates to rail [xr,yr,zr] coordinates

   % reshape inputs to row vectors

   m = size(xw,1); n = size(xw,2);
   xw = reshape(xw, 1, m*n);
   yw = reshape(yw, 1, m*n);
   zw = reshape(zw, 1, m*n);
   coords = [xw; yw; zw];

   % form transformation matrices and vectors

   R_r  = rotx(sol.meta.roll_r*180/pi);
   o_r  = [0; sol.meta.y_r; sol.meta.z_r];

   R_w  = rotx(sol.meta.roll_w*180/pi) * rotz(sol.meta.yaw_w*180/pi);
   o_w  = [sol.meta.x_w; sol.meta.y_w; sol.meta.z_w];

   % change coordinates from wheel to track to rail coordinates

   coords = R_r' * (o_w-o_r)*ones(1,m*n) + R_r' * R_w * coords;

   % extract & reshape to the original array sizes

   xr = reshape(coords(1,:), m, n);
   yr = reshape(coords(2,:), m, n);
   zr = reshape(coords(3,:), m, n);

end % function wheel_to_rail_coords

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ rot ] = rotx(roll_deg)

  sn = sin(roll_deg*pi/180); 
  cs = cos(roll_deg*pi/180);
  rot = [ 1,  0,   0;
          0, cs, -sn;
          0, sn,  cs];

end % function rotx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ rot ] = roty(pitch_deg)

  sn = sin(pitch_deg*pi/180); 
  cs = cos(pitch_deg*pi/180);
  rot = [ cs,  0,  sn;
           0,  1,   0;
         -sn,  0,  cs];

end % function roty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ rot ] = rotz(yaw_deg)

  sn = sin(yaw_deg*pi/180); 
  cs = cos(yaw_deg*pi/180);
  rot = [ cs, -sn, 0;
          sn,  cs, 0;
           0,   0, 1];

end % function rotz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ rgb ] = matlab_color( num )

   matlab_colors = [
            0    0.4470    0.7410
       0.8500    0.3250    0.0980
       0.9290    0.6940    0.1250
       0.4940    0.1840    0.5560
       0.4660    0.6740    0.1880
       0.3010    0.7450    0.9330
       0.6350    0.0780    0.1840
   ];
   num_colors = size(matlab_colors,1);

   % reshape to column vector, clip, get rows

   num = reshape(num, prod(size(num)), 1);
   num = max(1, min(num_colors, num));
   rgb = matlab_colors(num, :);

end % function matlab_color

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

