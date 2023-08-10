
function [ myopt ] = plot_2dspline( sol, slcs, opt )

% function [ opt ] = plot_2dspline( [sol], slcs, opt )
%
% create 3d surface plot of variable rail profile
%   sol     - optional results from CONTACT calculation, used for transforming to track coordinates
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

   % construct a local struct "myopt" in which default values are filled in for all available options.

   myopt = struct( ...
      'urange',     [],  ...
      'vrange',     [],  ...
      'xysteps',    [],  ...
      'show_slc',   [],  ...   % slice numbers that are highlighted
      'slc_color',  'r', ...   % color-spec for slices that are highlighted
      'show_feat',  [],  ...   % feature numbers that are highlighted
      'feat_color', 'r', ...   % color-spec for features that are highlighted
      'coordsys',   'fc', ...
      'view',       [110 20], ...
      'zoom',       [], ...
      'addplot',    0   ...
   );

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

   % Check whether user-supplied opt3-struct contains unknown options

   useropts = fieldnames(opt);
   ierror = 0;
   for i = 1:length(useropts)
      if (~isfield(myopt, useropts{i}))
         ierror = ierror + 1;
         disp(sprintf('Unknown option "%s" will be ignored',useropts{i}));
         if (ierror==1)
            disp(sprintf('You may use "opt=rmfield(opt,''%s'');" to remove the option',useropts{i}));
         end
      end
   end

   % Overwrite all values in "myopt" with user-supplied values

   myopts = fieldnames(myopt);
   for i = 1:length(myopts)
      if (isfield(opt, myopts{i}))
         myopt = setfield(myopt, myopts{i}, mygetfield(opt,myopts{i}));
      end
   end

   if (~any(size(myopt.slc_color,2)==[1 3]))
      disp('ERROR: slc_color should be column vector with 1 or 3 entries per row');
      myopt.slc_color = 'r';
   end
   if (~any(size(myopt.feat_color,2)==[1 3]))
      disp('ERROR: feat_color should be column vector with 1 or 3 entries per row');
      myopt.feat_color = 'r';
   end

   % start plotting

   if (~myopt.addplot)
      clf;
   end

   % set target step-sizes for plotting

   if (isempty(myopt.urange))
      myopt.urange = [ min(slcs.spl2d.ui), max(slcs.spl2d.ui) ];
   end
   if (isempty(myopt.vrange))
      myopt.vrange = [ min(slcs.spl2d.vj), max(slcs.spl2d.vj) ];
   end
   if (isempty(myopt.xysteps))
      xlen = max(max(slcs.xsurf)) - min(min(slcs.xsurf));
      ylen = max(max(slcs.ysurf)) - min(min(slcs.ysurf));
      myopt.xysteps = max(xlen/4, ylen/15);
   end
   if (isempty(myopt.zoom))
      xlen = max(max(slcs.xsurf)) - min(min(slcs.xsurf));
      ylen = max(max(slcs.ysurf)) - min(min(slcs.ysurf));
      myopt.zoom = [ max(1, 0.5*xlen/ylen) 1 1 ]; % data aspect ratio, e.g. [1000 1 1]
   end
   if (length(myopt.xysteps)==1)
      myopt.xysteps = [1 1] * myopt.xysteps;
   end

   % form profile surface at given x, all u - fine sampling

   [ xsurf, ysurf, zsurf ] = make_3d_surface( sol, slcs, myopt, myopt.urange, [], myopt.xysteps(1)/10, ...
                                                              myopt.vrange, [], []);

   % plot 90% transparent surface (alpha=0.1); color-value?

   csurf = -1 * ones(size(xsurf));
   l=surf(xsurf, ysurf, zsurf, csurf, 'EdgeAlpha',0.1, 'FaceAlpha',0.1);
   set(l,'tag','prr');

   % set view appropriate for rails with z positive downwards

   set(gca,'xdir','reverse', 'zdir','reverse');
   axis equal;
   view(myopt.view);
   hold on

   if (strcmp(myopt.coordsys,'fc'))
      xlabel('x_{fc} [mm]');
      ylabel('y_{r} [mm]');
      zlabel('z_{r} [mm]');
   end
   set(gca,'dataaspectratio',myopt.zoom);

   % determine transverse curves along the surface at fixed steps in x

   [ xcurv, ycurv, zcurv ] = make_3d_surface( sol, slcs, myopt, myopt.urange, [], myopt.xysteps(1), ...
                                                                        [], [], []);

   % plot transverse curves along the surface

   l = plot3(xcurv', ycurv', zcurv', 'color',[.5 .5 .5], 'linewidth',0.5);

   [~,imin] = min(abs(xcurv(:,1)));
   set(l(imin), 'color', [.3 .3 .3], 'linewidth',1);

   % determine longitudinal curves along the surface at fixed steps in s

   [ xcurv, ycurv, zcurv ] = make_3d_surface( sol, slcs, myopt, ...
                                     myopt.urange, [], myopt.xysteps(1)/10, [], [], myopt.xysteps(2));

   % plot curves along surface in longitudinal direction

   clear l ip;
   l = plot3(xcurv, ycurv, zcurv, 'color',[.5 .5 .5], 'linewidth',0.5);
   imid = ceil(length(l)/2);
   set(l(imid), 'color', [.3 .3 .3], 'linewidth',1);

   % plot transverse curves at selected slices

   if (~isempty(myopt.show_slc))
      islc = unique( sort( max(1, min(slcs.nslc, myopt.show_slc)) ));
      nslc = length(islc);
      ncol = size(myopt.slc_color,1);
      jv   = [ max(1, round(myopt.vrange(1))) : min(slcs.npnt,round(myopt.vrange(2))) ];
      nv   = length(jv);
      for i = 1 : nslc
         jslc = islc(i);
         if (slcs.s(jslc)>=myopt.urange(1) & slcs.s(jslc)<=myopt.urange(end))
            icol = mod((i-1), ncol) + 1;
            col  = myopt.slc_color(icol,:);
            if (length(col)==1 & isnumeric(col)), col = matlab_color(col); end
            plot3( slcs.xsurf(jslc,jv), slcs.ysurf(jslc,jv), slcs.zsurf(jslc,jv), 'color',col);
         end
      end
   end

   % plot longitudinal curves at selected features

   if (~isempty(myopt.show_feat))
      ifeat = unique( sort( max(1, min(slcs.nfeat, myopt.show_feat)) ));
      nfeat = length(ifeat);
      ncol  = size(myopt.feat_color,1);
      for j = 1 : nfeat
         jp   = slcs.iseg_p( ifeat(j) );
         [ xcurv, ycurv, zcurv ] = make_3d_surface( sol, slcs, myopt, ...
                       myopt.urange, [], myopt.xysteps(1)/10, slcs.vj(jp)*[1 1], [], myopt.xysteps(2));
         jcol = mod((j-1), ncol) + 1;
         col  = myopt.feat_color(jcol,:);
         if (length(col)==1 & isnumeric(col)), col = matlab_color(col); end
         % plot3( slcs.xsurf(:,jp), slcs.ysurf(:,jp), slcs.zsurf(:,jp), 'color',col );
         plot3( xcurv, ycurv, zcurv, 'color',col);
      end
   end

end % function plot_2dspline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ xsurf, ysurf, zsurf ] = make_3d_surface( sol, slcs, opt, ...
                                                              urange, du_true, du_appx, ...
                                                              vrange, dv_true, dv_appx, idebug)

% create 2D arrays (x,y,z) for a variable rail profile (slcs)
%  - the profile will be plotted for urange = [umin, umax]
%  - fixed steps of size du_true may be requested, or n-time fitting steps of approximately du_appx
%  - the whole profile will be used if vrange = [vmin, vmax] is left empty
%  - the profile will be sampled at vj-values in the profile close to uniform vj = [vmin: dv: vmax]
%  - all profile points will be used if both dv_true and dv_appx are left empty

   if (nargin< 5), du_true = []; end
   if (nargin< 6), du_appx = []; end
   if (nargin< 7), vrange  = []; end
   if (nargin< 8), dv_true = []; end
   if (nargin< 9), dv_appx = []; end
   if (nargin<10), idebug  = 0;  end

   % determine longitudinal positions ui for evaluation of surface

   if (~isempty(du_true))

      % if 'final' du_true is prescribed,
      %    set ui = { i * du_true } for appropriate i

      i0 = ceil( urange(1) / du_true );
      i1 = floor( urange(2) / du_true );
      ui = [i0 : i1] * du_true;
      if (idebug>=1), disp(sprintf('Using du_true=%4.2f, i0=%d, i1=%d',du_true,i0,i1)); end

   elseif (~isempty(du_appx))

      % if 'target' du_appx is given, 
      %    diuide [xmin,xmax] in odd #intervals of size du_true close to du_appx

      nintv   = floor( (urange(2) - urange(1)) / du_appx );
      if (mod(nintv,2)==1)  % prefer odd #lines for symmetric [-xmax, xmax]
         nintv = nintv + 1;
      end
      du = (urange(2) - urange(1)) / max(1,nintv);
      ui = urange(1) + [0:nintv] * du;
      if (idebug>=1), disp(sprintf('Using du_appx=%4.2f, nintv=%d',du_appx,nintv)); end

   else

      disp('Error: either du_true or du_appx must be given');
      return

   end

   % determine lateral positions vj for evaluation of surface

   if (isempty(vrange))
      vrange = [ slcs.vj(1), slcs.vj(end) ];
   end

   if (~isempty(dv_true))

      % if 'final' dv_true is prescribed,
      %    set vj = { j * dv_true } for appropriate j

      j0 = ceil(  (vrange(1)-slcs.vj(1)) / dv_true );
      j1 = floor( (vrange(2)-slcs.vj(1)) / dv_true );
      vj = [j0 : j1] * dv_true;
      if (idebug>=1), disp(sprintf('Using dv_true=%4.2f, j0=%d, j1=%d',dv_true,j0,j1)); end

   elseif (~isempty(dv_appx))

      % if 'target' dv_appx is given, 
      %    diuide [vmin,vmax] in odd #intervals of size dv_true close to dv_appx

      nintv = floor( (vrange(2) - vrange(1)) / dv_appx );
      dv = (vrange(2) - vrange(1)) / max(1,nintv);
      vj = vrange(1) + [0:nintv] * dv;
      if (idebug>=1), disp(sprintf('Using dv_appx=%4.2f, nintv=%d',dv_appx,nintv)); end

   else

      % if neither is given, use all profile points in range vrange

      j  = find( slcs.vj>=vrange(1) & slcs.vj<=vrange(2) );
      vj = slcs.vj(j);

   end

   % form profile surface at given ui and vj

   [ ~, xsurf, ysurf, zsurf ] = eval_2dspline( slcs.spl2d, ui, vj );

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

