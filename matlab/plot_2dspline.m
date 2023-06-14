
function [ opt ] = plot_2dspline( sol, slcs, opt )

% function [ opt ] = plot_2dspline( [sol], slcs, opt )
%
% create 3d surface plot of variable rail profile
%   sol     - optional results from contact calculation
%   slcs    - structure with variable profile as returned by read_slices
%   opt     - plot configuration; called without options, a struct is returned with default settings
%

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   if (nargin<1 | isempty(sol))
      sol = struct();
   end
   if (nargin<2 | isempty(slcs))
      slcs = struct();
   end
   if (nargin<3 | isempty(opt))
      opt = struct();
   end

   % fill default options

   if (~isfield(opt, 'xrange'))
      opt.xrange   = [];
   end
   if (~isfield(opt, 'urange'))
      opt.urange   = [];
   end
   if (~isfield(opt, 'xysteps'))
      opt.xysteps  = [];
   end
   if (~isfield(opt, 'show_slc'))   % slice numbers that are highlighted
      opt.show_slc = [];
   end
   if (~isfield(opt, 'slc_color'))  % color-spec for slices that are highlighted
      opt.slc_color = 'r';
   end
   if (~isempty(opt.slc_color) & ~any(size(opt.slc_color,2)==[1 3]))
      disp('ERROR: slc_color should be column vector with 1 or 3 entries per row');
      opt.slc_color = 'r';
   end
   if (~isfield(opt, 'show_feat'))  % feature numbers that are highlighted
      opt.show_feat = [];
   end
   if (~isfield(opt, 'feat_color')) % color-spec for features that are highlighted
      opt.feat_color = 'r';
   end
   if (~isempty(opt.feat_color) & ~any(size(opt.feat_color,2)==[1 3]))
      disp('ERROR: feat_color should be column vector with 1 or 3 entries per row');
      opt.feat_color = 'r';
   end
   if (~isfield(opt, 'coordsys'))
      opt.coordsys = 'fc';
   end
   if (~isfield(opt, 'view'))
      opt.view     = [110 20];
   end
   if (~isfield(opt, 'zoom'))
      opt.zoom     = [];
   end
   if (~isfield(opt, 'addplot'))
      opt.addplot  = 0;
   end

   if (isempty(fields(slcs)))   % called with no arguments: return default options
      return;
   end
   if (~isfield(slcs, 's'))
      disp('Error: 2nd argument does not look like a slcs-structure, no s-positions')
      return;
   end
   if (~isfield(slcs, 'spl2d'))
      disp('Error: slcs does not have 2d spline data');
      return;
   end

   if (~opt.addplot)
      clf;
   end

   % set target step-sizes for plotting

   if (isempty(opt.xrange))
      opt.xrange = [ min(slcs.s), max(slcs.s) ];
   end
   if (isempty(opt.urange))
      opt.urange = [ min(slcs.uj), max(slcs.uj) ];
   end
   if (isempty(opt.xysteps))
      xlen = opt.xrange(2) - opt.xrange(1);
      ylen = max(max(slcs.ysurf)) - min(min(slcs.ysurf));
      opt.xysteps = max(xlen/4, ylen/15);
   end
   if (isempty(opt.zoom))
      xlen = opt.xrange(2) - opt.xrange(1);
      ylen = max(max(slcs.ysurf)) - min(min(slcs.ysurf));
      opt.zoom = [ max(1, 0.5*xlen/ylen) 1 1 ]; % data aspect ratio, e.g. [1000 1 1]
   end
   if (length(opt.xysteps)==1)
      opt.xysteps = [1 1] * opt.xysteps;
   end

   % form profile surface at given x, all u - fine sampling

   [ xsurf, ysurf, zsurf ] = make_3d_surface( sol, slcs, opt, opt.xrange, [], opt.xysteps(1)/10, ...
                                                              opt.urange, [], []);

   % plot 90% transparent surface (alpha=0.1); color-value?

   csurf = -1 * ones(size(xsurf));
   l=surf(xsurf, ysurf, zsurf, csurf, 'EdgeAlpha',0.1, 'FaceAlpha',0.1);
   set(l,'tag','prr');

   % set view appropriate for rails with z positive downwards

   set(gca,'xdir','reverse', 'zdir','reverse');
   axis equal;
   view(opt.view);
   hold on

   if (strcmp(opt.coordsys,'fc'))
      xlabel('x_{fc} [mm]');
      ylabel('y_{r} [mm]');
      zlabel('z_{r} [mm]');
   end
   set(gca,'dataaspectratio',opt.zoom);

   % determine transverse curves along the surface at fixed steps in x

   [ xcurv, ycurv, zcurv ] = make_3d_surface( sol, slcs, opt, opt.xrange, [], opt.xysteps(1), ...
                                                                        [], [], []);

   % plot transverse curves along the surface

   l = plot3(xcurv', ycurv', zcurv', 'color',[.5 .5 .5], 'linewidth',0.5);

   [~,imin] = min(abs(xcurv(:,1)));
   set(l(imin), 'color', [.3 .3 .3], 'linewidth',1);

   % determine longitudinal curves along the surface at fixed steps in s

   [ xcurv, ycurv, zcurv ] = make_3d_surface( sol, slcs, opt, opt.xrange, [], opt.xysteps(1)/10, ...
                                                                  [], [], opt.xysteps(2));

   % plot curves along surface in longitudinal direction

   clear l ip;
   l = plot3(xcurv, ycurv, zcurv, 'color',[.5 .5 .5], 'linewidth',0.5);
   imid = ceil(length(l)/2);
   set(l(imid), 'color', [.3 .3 .3], 'linewidth',1);

   % plot transverse curves at selected slices

   if (~isempty(opt.show_slc))
      islc = unique( sort( max(1, min(slcs.nslc, opt.show_slc)) ));
      nslc = length(islc);
      ncol = size(opt.slc_color,1);
      ju   = [ max(1, round(opt.urange(1))) : min(slcs.npnt,round(opt.urange(2))) ];
      nu   = length(ju);
      for i = 1 : nslc
         jslc = islc(i);
         icol = mod((i-1), ncol) + 1;
         col  = opt.slc_color(icol,:);
         if (length(col)==1 & isnumeric(col)), col = matlab_color(col); end
         plot3( slcs.s(jslc)*ones(1,nu), slcs.ysurf(jslc,ju), slcs.zsurf(jslc,ju), 'color',col);
      end
   end

   % plot longitudinal curves at selected features

   if (~isempty(opt.show_feat))
      ifeat = unique( sort( max(1, min(slcs.nfeat, opt.show_feat)) ));
      nfeat = length(ifeat);
      ncol  = size(opt.feat_color,1);
      for i = 1 : nfeat
         ip   = slcs.iseg_p( ifeat(i) );
         [ xcurv, ycurv, zcurv ] = make_3d_surface( sol, slcs, opt, opt.xrange, [], opt.xysteps(1)/10, ...
                                                                  slcs.uj(ip)*[1 1], [], opt.xysteps(2));
         icol = mod((i-1), ncol) + 1;
         col  = opt.feat_color(icol,:);
         if (length(col)==1 & isnumeric(col)), col = matlab_color(col); end
         % plot3( slcs.s, slcs.ysurf(:,ip), slcs.zsurf(:,ip), 'color',col );
         plot3( xcurv, ycurv, zcurv, 'color',col);
      end
   end

end % function show_profiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ xsurf, ysurf, zsurf ] = make_3d_surface( sol, slcs, opt, ...
                                                              xrange, dx_true, dx_appx, ...
                                                              urange, du_true, du_appx)

% create 2D arrays (x,y,z) for a variable rail profile (slcs)
%  - the profile will be plotted for xrange = [xmin, xmax]
%  - fixed steps of size dx_true may be requested, or n-time fitting steps of approximately dx_appx
%  - the whole profile will be used if urange = [umin, umax] is left empty
%  - the profile will be sampled at uj-values in the profile close to uniform uj = [umin: du: umax]
%  - all profile points will be used if both du_true and du_appx are left empty

   if (nargin< 5), dx_true = []; end
   if (nargin< 6), dx_appx = []; end
   if (nargin< 7), urange  = []; end
   if (nargin< 8), du_true = []; end
   if (nargin< 9), du_appx = []; end

   % determine longitudinal positions xi for evaluation of surface

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

   % determine lateral positions uj for evaluation of surface

   if (isempty(urange))
      urange = [ slcs.uj(1), slcs.uj(end) ];
   end

   if (~isempty(du_true))

      % if 'final' du_true is prescribed,
      %    set uj = { j * du_true } for appropriate j

      j0 = ceil(  (urange(1)-slcs.uj(1)) / du_true );
      j1 = floor( (urange(2)-slcs.uj(1)) / du_true );
      uj = [j0 : j1] * du_true;

   elseif (~isempty(du_appx))

      % if 'target' du_appx is given, 
      %    divide [umin,umax] in odd #intervals of size du_true close to du_appx

      nintv = floor( (urange(2) - urange(1)) / du_appx );
      du = (urange(2) - urange(1)) / max(1,nintv);
      uj = urange(1) + [0:nintv] * du;

   else

      % if neither is given, use all profile points

      j  = find( slcs.uj>=urange(1) & slcs.uj<=urange(2) );
      uj = slcs.uj(j);

   end

   % form profile surface at given xi and uj

   [ ~, ysurf, zsurf ] = eval_2dspline( slcs.spl2d, xi, uj );

   nu    = length(uj);
   xsurf = xi' * ones(1,nu);

   % transform to track coordinates if needed

   if (strcmp(opt.coordsys, 'tr'))

      [xsurf, ysurf, zsurf] = rail_to_track_coords(sol, xsurf, ysurf, zsurf);

   end

end % function make_3d_surface

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

   Rrr  = rotx(sol.meta.rollrr*180/pi);
   orr  = [0; sol.meta.yrr; sol.meta.zrr];

   if (sol.d_digit~=4)

      % change coordinates from contact to rail to track coordinates

      coords = orr*ones(1,m*n) + Rrr * oref*ones(1,m*n) + Rrr * Rref * coords;

   else

      % find oref in the rail profile

      dst = sqrt( (sol.prr.ProfileY-sol.meta.ycp_r).^2 + ...
                  (sol.prr.ProfileZ-sol.meta.zcp_r).^2 );
      [~,ix] = min(dst);
      rg = [max(1,ix-5) : min(length(sol.prr.ProfileY),ix+5)];

      sref = interp1(sol.prr.ProfileY(rg), sol.prr.ProfileS(rg), sol.meta.ycp_r);
      if (isnan(sref))
         disp(sprintf('ERROR: prr seems different from prr used in simulation (needs mirroring?).'))
         sref = sol.prr.ProfileS(ix);
      end

      % find grid points in the rail profile

      sc = yc;
      sr = sref + sc;
      yr_at_s = interp1(sol.prr.ProfileS, sol.prr.ProfileY, sr);
      zr_at_s = interp1(sol.prr.ProfileS, sol.prr.ProfileZ, sr);

      % determine n-vectors on the rail profile

      ds = 0.1;
      dy = (interp1(sol.prr.ProfileS, sol.prr.ProfileY, sr+ds) - yr_at_s) / ds;
      dz = (interp1(sol.prr.ProfileS, sol.prr.ProfileZ, sr+ds) - zr_at_s) / ds;
      nvec = [zeros(1,m*n); -dz; dy] ./ ([1;1;1] * sqrt(dy.^2+dz.^2));

      % hack: slight raise to put the data above the surfaces
      % zc = zc - 0.3;

      coords = [xc+oref(1); yr_at_s; zr_at_s] + nvec .* ([1;1;1] * zc);

      % change coordinates from rail to track coordinates

      coords = orr*ones(1,m*n) + Rrr * coords;

   end

   % extract components & reshape to the original array sizes

   xtr = reshape(coords(1,:), m, n);
   ytr = reshape(coords(2,:), m, n);
   ztr = reshape(coords(3,:), m, n);
        % TODO: correction for \Delta\phi_rr
        % delta_tr = sol.meta.deltcp_r + sol.meta.rollrr;
   delta_tr = sol.meta.deltcp_r + sol.meta.rollrr;

   % mirror y-coordinates for left rail/roller

   if (is_left_side(sol))
      ytr = -ytr;
      delta_tr = -delta_tr;
   end
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
   oref = [sol.meta.xcp_r; sol.meta.ycp_r; sol.meta.zcp_r];

   if (sol.d_digit~=4)

      % change coordinates from contact to rail coordinates

      coords = oref*ones(1,m*n) + Rref * coords;

   else

      % find oref in the rail profile

      dst = sqrt( (sol.prr.ProfileY-sol.meta.ycp_r).^2 + ...
                  (sol.prr.ProfileZ-sol.meta.zcp_r).^2 );
      [~,ix] = min(dst);
      rg = [max(1,ix-5), min(length(sol.prr.ProfileY),ix+5)];

      sref = interp1(sol.prr.ProfileY(rg), sol.prr.ProfileS(rg), sol.meta.ycp_r);

      % find grid points in the rail profile

      sc = yc;
      sr = sref + sc;
      yr_at_s = interp1(sol.prr.ProfileS, sol.prr.ProfileY, sr);
      zr_at_s = interp1(sol.prr.ProfileS, sol.prr.ProfileZ, sr);

      % determine n-vectors on the rail profile

      ds = 0.1;
      dy = (interp1(sol.prr.ProfileS, sol.prr.ProfileY, sr+ds) - yr_at_s) / ds;
      dz = (interp1(sol.prr.ProfileS, sol.prr.ProfileZ, sr+ds) - zr_at_s) / ds;
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

function [ xtr, ytr, ztr ] = rail_to_track_coords(sol, xr, yr, zr);

% transform rail [xr,yr,zr]-coordinates to track [xtr,ytr,ztr] coordinates

   % reshape inputs to row vectors

   m = size(xr,1); n = size(xr,2);
   xr = reshape(xr, 1, m*n);
   yr = reshape(yr, 1, m*n);
   zr = reshape(zr, 1, m*n);
   coords = [xr; yr; zr];

   % form transformation matrices and vectors

   Rrr  = rotx(sol.meta.rollrr*180/pi);
   orr  = [0; sol.meta.yrr; sol.meta.zrr];

   % change coordinates from rail to track coordinates

   coords = orr * ones(1,m*n) + Rrr * coords;

   % extract & reshape to the original array sizes

   xtr = reshape(coords(1,:), m, n);
   ytr = reshape(coords(2,:), m, n);
   ztr = reshape(coords(3,:), m, n);

   % mirror y-coordinates for left rail/roller

   if (is_left_side(sol))
      ytr = -ytr;
   end
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

   % form transformation matrices and vectors

   Rrw  = rotx(sol.meta.rollrw*180/pi) * rotz(sol.meta.yawrw*180/pi);
   orw  = [sol.meta.xrw; sol.meta.yrw; sol.meta.zrw];

   % change coordinates from wheel to track coordinates

   coords = orw * ones(1,m*n) + Rrw * coords;

   % extract & reshape to the original array sizes

   xtr = reshape(coords(1,:), m, n);
   ytr = reshape(coords(2,:), m, n);
   ztr = reshape(coords(3,:), m, n);

   % mirror y-coordinates for left rail/roller

   if (is_left_side(sol))
      ytr = -ytr;
   end
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

   Rrr  = rotx(sol.meta.rollrr*180/pi);
   orr  = [0; sol.meta.yrr; sol.meta.zrr];

   Rrw  = rotx(sol.meta.rollrw*180/pi) * rotz(sol.meta.yawrw*180/pi);
   orw  = [sol.meta.xrw; sol.meta.yrw; sol.meta.zrw];

   % change coordinates from wheel to track to rail coordinates

   coords = Rrr' * (orw-orr)*ones(1,m*n) + Rrr' * Rrw * coords;

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

