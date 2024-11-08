
function [ myopt ] = plot_2dspline( slcs, opt )

% function [ opt ] = plot_2dspline( slcs, opt )
%
% create 3d surface plot of variable rail profile
%   slcs    - structure with variable profile as returned by read_slices
%   opt     - plot configuration; called without options, a struct is returned with default settings
%
% plot options:
%   typplot:     type of plot: 'surf' (default), 'refl', 'topol'
%   field:       compute/show surface indicator: 'none', 'dydu', 'dzdu', 'dydv', 'dzdv'
%   urange:      range of surface parameter u for longitudinal direction
%   vrange:      range of surface parameter v for lateral direction
%   xysteps:     (approx) distance between drawn lines on surface
%   show_slc:    list of slice numbers that are highlighted
%   slc_color:   color-spec for slices that are highlighted
%   show_feat:   list of feature numbers that are highlighted
%   feat_color:  color-spec for features that are highlighted
%   view:        matlab view direction, e.g. [50 25] (azimuth, elevation), can be 'rail' or 'default'
%   refl_avec:   light source direction for reflection plot, e.g. [0,1,0] or [100, 30] (az, el)
%   refl_dth:    interval spacing d theta for reflection plot, e.g. 30 deg (default)
%   zoom:        data aspect ratio, e.g. [1000 1 1] (rail) or [0.1 1 1] (wheel, theta)
%   addplot:     clear (0) or do not clear (1) the figure before plotting.

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   if (nargin<1 | isempty(slcs))
      slcs = struct();
   end
   if (nargin<2 | isempty(opt))
      opt = struct();
   end

   % construct a local struct "myopt" in which default values are filled in for all available options.

   myopt = struct( ...
      'typplot',    'surf',  ...
      'field',      [],  ...
      'urange',     [],  ...
      'vrange',     [],  ...
      'xysteps',    [],  ...
      'show_slc',   [],  ...
      'slc_color',  'r', ...
      'show_feat',  [],  ...
      'feat_color', 'r', ...
      'view',       'default', ...
      'refl_avec',  [0,1,0], ...
      'refl_dth',   30, ...
      'zoom',       [], ...
      'addplot',    0   ...
   );

   if (isempty(fields(slcs)))   % called with no arguments: return default options
      return;
   end
   if (~isfield(slcs, 'u'))
      disp('Error: 2nd argument does not look like a slcs-structure, no u-positions')
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

   if (strcmp(myopt.view,'rail'))
      % wheel-rail view: 2D, rolling direction == plot y-direction
      myopt.view=[90 -90];
   elseif (strcmp(myopt.view,'default'))
      if (strcmp(myopt.typplot,'topol'))
         myopt.view=[90 -90];
      else
         myopt.view=[110 20];
      end
   end 

   % start plotting

   if (~myopt.addplot)
      clf;
   end

   if (max(myopt.urange)<slcs.spl2d.tui(4) | min(myopt.urange)>slcs.spl2d.tui(end-3))
      disp(sprintf('opt.urange=[%3.1f,%3.1f] outside spline range u=[%3.1f,%3.1f]', ...
           myopt.urange, slcs.spl2d.tui([4,end-3])));
      return;
   end
   if (max(myopt.vrange)<slcs.spl2d.tvj(4) | min(myopt.vrange)>slcs.spl2d.tvj(end-3))
      disp(sprintf('opt.vrange=[%3.1f,%3.1f] outside spline range v=[%3.1f,%3.1f]', ...
           myopt.vrange, slcs.spl2d.tvj([4,end-3])));
      return;
   end

   % set target step-sizes for plotting

   if (isempty(myopt.urange))
      xlen = max(max(slcs.xsurf)) - min(min(slcs.xsurf));
   else
      xlen = max(opt.urange) - min(opt.urange);
   end
   if (isempty(myopt.urange))
      myopt.urange = [ min(slcs.spl2d.ui), max(slcs.spl2d.ui) ];
   end
   if (isempty(myopt.vrange))
      myopt.vrange = [ min(slcs.spl2d.vj), max(slcs.spl2d.vj) ];
   end
   if (isempty(myopt.xysteps))
      ylen = max(slcs.spl2d.vj) - min(slcs.spl2d.vj);
      myopt.xysteps = [xlen/4, ylen/10]; % using 4/10 intervals == 5/11 lines
   end
   if (isempty(myopt.zoom))
      if (strcmp(myopt.typplot,'topol'))
         ylen = max(max(slcs.vj)) - min(min(slcs.vj));
      else
         ylen = max(max(slcs.ysurf)) - min(min(slcs.ysurf));
      end
      myopt.zoom = [ 0.5*xlen/ylen 1 1 ];       % wheel: data aspect ratio, e.g. [0.05 1 1]
      if (~slcs.is_wheel)
         myopt.zoom(1) = max(1, myopt.zoom(1)); % rail: data aspect ratio could be [1000 1 1]
      end
   end
   if (length(myopt.xysteps)==1)
      myopt.xysteps = [1 1] * myopt.xysteps;
   end

   % form profile surface at given x, all v - fine sampling

   [ xsurf, ysurf, zsurf, ui_samp, vj_samp ] = ...
              make_3d_surface( slcs, myopt, myopt.urange, [], myopt.xysteps(1)/20, myopt.vrange, [], []);

   % compute surface indicator at sampling positions [ui_samp] x [vj_samp]

   csurf = -1 * ones(size(xsurf));
   if (any(strcmp(myopt.field,{'dist','disty','distz'})))
      if (slcs.u_intpol==1 & isfield(slcs,'appx'))
                % recompute [xyz]surf for dist in topology plot
         [ ~, xintp, yintp, zintp ] = eval_2dspline( slcs.spl2d, ui_samp, vj_samp );
         [ ~, xappx, yappx, zappx ] = eval_2dspline( slcs.appx, ui_samp, vj_samp );
      elseif (slcs.u_intpol==2 & isfield(slcs,'intp'))
         [ ~, xappx, yappx, zappx ] = eval_2dspline( slcs.spl2d, ui_samp, vj_samp );
         [ ~, xintp, yintp, zintp ] = eval_2dspline( slcs.intp, ui_samp, vj_samp );
      end
      dist  = sqrt( (xintp-xappx).^2 + (yintp-yappx).^2 + (zintp-zappx).^2 );
      disty = yintp - yappx;
      distz = zintp - zappx;
      if     (strcmp(myopt.field,'dist'))
         csurf = dist;
      elseif (strcmp(myopt.field,'disty'));
         csurf = disty;
      elseif (strcmp(myopt.field,'distz'));
         csurf = distz;
      end
   elseif (any(strcmp(myopt.field,{'dydu','dzdu','dydv','dzdv','dzdx'})))
      [dxdu, dydu, dzdu] = eval_2dspline_deriv(slcs.spl2d, ui_samp, vj_samp, 1);
      [dxdv, dydv, dzdv] = eval_2dspline_deriv(slcs.spl2d, ui_samp, vj_samp, 2);
      if     (strcmp(myopt.field,'dydu'))
         csurf = dydu;
      elseif (strcmp(myopt.field,'dzdu'))
         csurf = dzdu;
      elseif (strcmp(myopt.field,'dydv'))
         csurf = dydv;
      elseif (strcmp(myopt.field,'dzdv'))
         csurf = dzdv;
      elseif (strcmp(myopt.field,'dzdx'))
         % fails on cross-sections with vertical slopes dy/dv=0
         csurf = dzdu - (dydu ./ dydv) .* dzdv;
      end
   end

   % plot 90% transparent surface (alpha=0.1); color-value according to field computed above

   if (any(strcmp(myopt.typplot,{'surf','topol'})))
      if (isempty(myopt.field) | strcmp(myopt.field,'none'))
         l=surf(xsurf, ysurf, zsurf, csurf, 'EdgeAlpha',0.1, 'FaceAlpha',0.1);
      else
         l=surf(xsurf, ysurf, zsurf, csurf, 'EdgeAlpha',0.1, 'FaceAlpha',1);
         h=colorbar; ylabel(h, myopt.field);
      end
   end

   % compute and plot reflection lines

   if (strcmp(myopt.typplot, 'refl'))
      % using 4 positions per knot span
      ix = 3 + find( slcs.spl2d.tui(4:end-3)>= myopt.urange(1) & slcs.spl2d.tui(4:end-3)<=myopt.urange(end) );
      ui_samp = refine_intv( [myopt.urange(1) ; slcs.spl2d.tui(ix) ; myopt.urange(end) ], 12);
      if (isempty(ix))
         disp(sprintf('ERROR: spline surface has no data for u=[%3.1f,%3.1f], available u=[%3.1f,%3.1f]',...
                myopt.urange, slcs.spl2d.tui(4), slcs.spl2d.tui(end-3)));
         return
      end
      ix = 3 + find( slcs.spl2d.tvj(4:end-3)>= myopt.vrange(1) & slcs.spl2d.tvj(4:end-3)<=myopt.vrange(end) );
      vj_samp = refine_intv( slcs.spl2d.tvj(ix), 12 );
      if (isempty(ix))
         disp(sprintf('ERROR: spline surface has no data for v=[%3.1f,%3.1f], available v=[%3.1f,%3.1f]',...
                myopt.vrange, slcs.spl2d.tvj(4), slcs.spl2d.tvj(end-3)));
         return
      end

      [xsurf, ysurf, zsurf, th_rfl] = make_reflec( slcs.spl2d, myopt.refl_avec, myopt.view, ui_samp, vj_samp);
      th_col = floor(th_rfl/myopt.refl_dth) * myopt.refl_dth;

      % tmp   = parula(12);
      % m_bin = repmat(tmp([3,12],:), n_bin/2, 1);
      n_bin = round(360 / myopt.refl_dth);
      tmp   = [0 0 0; 1 1 .6];
      m_bin = repmat(tmp, n_bin/2, 1);

      l=surf(xsurf, ysurf, zsurf, th_col, 'EdgeColor','none', 'FaceAlpha',0.5);
      colormap(m_bin);
   end

   set(l,'Tag','slcs');

   % set view appropriate for rails with z positive downwards

   if (~strcmp(myopt.typplot,'topol'))
      set(gca,'xdir','reverse', 'zdir','reverse');
   end
   axis equal;
   view(myopt.view);
   hold on

   if (strcmp(myopt.typplot,'topol'))
      xlabel('u [-]');
      ylabel('v [-]');
   elseif (slcs.is_wheel)
      xlabel('\theta_{w} [rad]');
      ylabel('y_{w} [mm]');
      zlabel('dr_{w} [mm]');
   else
      xlabel('u_{fc} [mm]');
      ylabel('y_{r} [mm]');
      zlabel('z_{r} [mm]');
   end
   set(gca,'dataaspectratio',myopt.zoom);

   % determine transverse curves along the surface at selected u (x) positions

   [ xcurv, ycurv, zcurv ] = make_3d_surface( slcs, myopt, myopt.urange, [], myopt.xysteps(1), ...
                                                                        myopt.vrange, [], []);
   % plot transverse curves along the surface

   l = plot3(xcurv', ycurv', zcurv', 'color',[.5 .5 .5], 'linewidth',0.5, 'Tag','transverse');

   % highlight curve closest to x==0
   % [~,imin] = min(abs(xcurv(:,1)));
   % set(l(imin), 'color', [.3 .3 .3], 'linewidth',1);

   % determine longitudinal curves along the surface at selected v (y) positions

   [ xcurv, ycurv, zcurv ] = make_3d_surface( slcs, myopt, ...
                                 myopt.urange, [], myopt.xysteps(1)/20, myopt.vrange, [], myopt.xysteps(2));

   % plot curves along surface in longitudinal direction

   clear l ip;
   l = plot3(xcurv, ycurv, zcurv, 'color',[.5 .5 .5], 'linewidth',0.5, 'Tag','longitudinal');

   % highlight curve closest to mid of vrange used
   % imid = ceil(length(l)/2);
   % set(l(imid), 'color', [.3 .3 .3], 'linewidth',1);

   % plot transverse curves at selected slices

   if (~isempty(myopt.show_slc))
      islc = unique( sort( max(1, min(slcs.nslc, myopt.show_slc)) ));
      nslc = length(islc);
      ncol = size(myopt.slc_color,1);
      jv   = [ max(1, round(myopt.vrange(1))) : min(slcs.npnt,round(myopt.vrange(2))) ];
      nv   = length(jv);
      for i = 1 : nslc
         jslc = islc(i);
         if (slcs.u(jslc)>=myopt.urange(1) & slcs.u(jslc)<=myopt.urange(end))
            icol = mod((i-1), ncol) + 1;
            col  = myopt.slc_color(icol,:);
            if (length(col)==1 & isnumeric(col)), col = matlab_color(col); end
            if (1==1)
               [ xcurv, ycurv, zcurv ] = make_3d_surface( slcs, myopt, slcs.u(jslc)*[1 1], [], 1, ...
                                                                myopt.vrange, [], myopt.xysteps(2)/20);
            else
               xcurv = slcs.xsurf(jslc,jv); ycurv = slcs.ysurf(jslc,jv); zcurv = slcs.zsurf(jslc,jv);
            end
            plot3( xcurv, ycurv, zcurv, 'color',col, 'linewidth',1, 'Tag','slice');
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
         [ xcurv, ycurv, zcurv ] = make_3d_surface( slcs, myopt, ...
                       myopt.urange, [], myopt.xysteps(1)/20, slcs.vj(jp)*[1 1], [], myopt.xysteps(2));
         jcol = mod((j-1), ncol) + 1;
         col  = myopt.feat_color(jcol,:);
         if (length(col)==1 & isnumeric(col)), col = matlab_color(col); end
         % plot3( slcs.xsurf(:,jp), slcs.ysurf(:,jp), slcs.zsurf(:,jp), 'color',col );
         plot3( xcurv, ycurv, zcurv, 'color',col, 'linewidth',1, 'Tag','feature');
      end
   end

end % function plot_2dspline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ xsurf, ysurf, zsurf, ui, vj ] = make_3d_surface( slcs, opt, ...
                                                              urange, du_true, du_appx, ...
                                                              vrange, dv_true, dv_appx, idebug)

% create 2D arrays (x,y,z) for a variable rail profile (slcs)
%  - the profile will be plotted for urange = [umin, umax]
%  - fixed steps of size du_true may be requested, or n-time fitting steps of approximately du_appx
%  - the whole profile will be used if vrange = [vmin, vmax] is left empty
%  - the profile will be sampled at vj-values in the profile close to uniform vj = [vmin: dv: vmax]
%  - all profile points will be used if both dv_true and dv_appx are left empty

   if (nargin<4), du_true = []; end
   if (nargin<5), du_appx = []; end
   if (nargin<6), vrange  = []; end
   if (nargin<7), dv_true = []; end
   if (nargin<8), dv_appx = []; end
   if (nargin<9), idebug  = 0;  end

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
      %    divide [xmin,xmax] in odd #intervals of size du_true close to du_appx

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
      %    divide [vmin,vmax] in odd #intervals of size dv_true close to dv_appx

      nintv = floor( (vrange(2) - vrange(1)) / dv_appx );
      dv = (vrange(2) - vrange(1)) / max(1,nintv);
      vj = vrange(1) + [0:nintv] * dv;
      if (idebug>=1), disp(sprintf('Using dv_appx=%4.2f, nintv=%d',dv_appx,nintv)); end

   else

      % if neither is given, use all profile points in range vrange

      j  = find( slcs.vj>=vrange(1) & slcs.vj<=vrange(2) );
      vj = slcs.vj(j);

   end

   ui = max(ui, slcs.spl2d.tui(4)); ui = min(ui, slcs.spl2d.tui(end-3)); % eliminate round-off error
   vj = max(vj, slcs.spl2d.tvj(4)); vj = min(vj, slcs.spl2d.tvj(end-3));
   if (size(ui,2)>1), ui = ui'; end     % ui: column vector nu x 1
   if (size(vj,1)>1), vj = vj'; end     % vj: row vector 1 x nv

   % form profile surface at given ui and vj

   [ ~, xsurf, ysurf, zsurf ] = eval_2dspline( slcs.spl2d, ui, vj );

   % for a topology plot, output (usurf, vsurf, 1) instead of [xyz]surf

   if (strcmp(opt.typplot,'topol'))
      xsurf = ui * ones(1,length(vj));
      ysurf = ones(length(ui),1) * vj;
      ix = find(~isnan(zsurf));
      zsurf(ix) = 1;
   end

end % function make_3d_surface

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ v_rfn ] = refine_intv( v_in, k_rfn )

% refine intervals v_in(1) : v_in(2), v_in(2) : v_in(3), ... into k_rfn>=2 pieces each

   n_in  = length( v_in );
   n_rfn = k_rfn * (n_in-1) + 1;
   v_rfn = zeros(n_rfn,1);
   f_rfn = [0 : 1/k_rfn : 1];

   for i_in = 1 : n_in-1
      i_out = k_rfn * (i_in-1) + [1 : k_rfn+1];
      dv_in = v_in(i_in+1) - v_in(i_in);
      v_rfn(i_out) = v_in(i_in) + f_rfn * dv_in;
   end
end % function refine_intv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
