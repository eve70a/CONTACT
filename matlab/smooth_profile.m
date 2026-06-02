
function [ p_out, spl ] = smooth_profile( p_in, lfilt_arg, lambda_arg, ds_sampl_arg, ikinks_arg, ...
                                          iaccel_arg, use_wgt_arg, use_bspline_arg, ds_bspl_arg, idebug_arg );

% function [ opt ]        = smooth_profile;
% function [ p_out, spl ] = smooth_profile( p_in, opt );
% function [ p_out, spl ] = smooth_profile( p_in, [l_filt], [lambda], [ds_sampl], [ikinks], [iaccel], 
%                                                           [use_wgt], [use_bspline], [ds_bspl], [idebug] );
%
% compute smoothed profile using parametric smoothing spline approximation with smoothness
% parameters l_filt [mm] or lambda.
%
% p_in may be a struct ('ProfileY', 'ProfileZ') or array [y_prf, z_prf]
%
% re-sample at step-size ds_sampl, if given, or use s-positions of input profile
%
% enable detection of kinks (-1, default), disable ([]), or set manually [kink1, kink2,...]
% enable detection of accelerations (-1, default), disable ([]) or set manually [accel1, accel2, ...]
% use uniform weighting (use_wgt<=0), by segment length (1, default), increase wgt for curvature (2)
% using PP-spline (use_bspline=0) or B-spline (use_bspline=1, default) with target step ds_bspl.
%
% default method:   opt.ikinks=-1, opt.use_wgt=1, opt.use_bspline=1, opt.ds_bspl=[]
% default settings: opt.l_filt=[], opt.ds_sampl=[]

% Copyright 2008-2026 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% construct a local struct "myopt" in which default values are filled in for all available options.

myopt = struct( ...
   'l_filt',       [], ...
   'lambda',       [], ...
   'ds_sampl',     [], ...
   'ikinks',       -1, ...
   'iaccel',       [], ...
   'use_wgt',       1, ...
   'use_bspline',   1, ...
   'ds_bspl',      [], ...
   'idebug',        0, ...
   'kink_high',      pi/6, ...  % kink detection, high/low dalpha, window +/-wid
   'kink_low',       pi/6 / 5, ...
   'kink_wid',        2.0, ...
   'corrc_cos_thrs',  0.5, ...  % correction in of summed chord length at cos(dalph)<thrs
   'corrc_fds_min',   0.1, ...
   'corrc_kmax',       10  ...
);

% If the user has not supplied any arguments, return default options in first output argument

if (nargin<=0)
   p_out = myopt;
   spl   = [];
   return
end

% If the user supplied an opt-struct, rename to opt
% else, start from the default, overwrite with user-options

if (nargin<=2 & isstruct(lfilt_arg))

   opt = lfilt_arg;

else

   opt = myopt;

   if (nargin>=2 & ~isempty(lfilt_arg))
      opt.l_filt   = lfilt_arg;
   end
   if (nargin>=3 & ~isempty(lambda_arg))
      opt.lambda   = lambda_arg;
   end
   if (nargin>=4 & ~isempty(ds_sampl_arg))
      opt.ds_sampl = ds_sampl_arg;
   end
   if (nargin>=5)
      opt.ikinks   = ikinks_arg;
   end
   if (nargin>=6)
      opt.iaccel   = iaccel_arg;
   end
   if (nargin>=7 & ~isempty(use_wgt_arg))
      opt.use_wgt  = use_wgt_arg;
   end
   if (nargin>=8 & ~isempty(use_bspline_arg))
      opt.use_bspline = use_bspline_arg;
   end
   if (nargin>=9)
      opt.ds_bspl = ds_bspl_arg;
   end
   if (nargin>=10 & ~isempty(idebug_arg))
      opt.idebug  = idebug_arg;
   end

end

is_wheel = 0;
if (isstruct(p_in))
   is_wheel = p_in.is_wheel;
end

if (opt.use_bspline)
   opt.use_deriv = 3;
else
   opt.use_deriv = 2;
end
opt.use_repl = 1;

% compute lambda or l_filt from each other

if     ( isempty(opt.l_filt) &  isempty(opt.lambda))
   disp('Either l_filt or lambda must be given.');
   return;
elseif (~isempty(opt.l_filt) & ~isempty(opt.lambda))
   disp('l_filt and lambda may not both be given.');
   return;
elseif (~isempty(opt.l_filt))
   if (opt.use_deriv==3)
      opt.lambda = (opt.l_filt / (2*pi))^6;
   else
      opt.lambda = (opt.l_filt / (2*pi))^4;
   end
else
   if (opt.use_deriv==3)
      opt.l_filt = 2*pi * opt.lambda^(1/6);
   else
      opt.l_filt = 2*pi * opt.lambda^(1/4);
   end
end
if (isempty(opt.ds_bspl))
   ds_bspl = max(0.001, opt.l_filt/4);
end

% obtain y, z and s (optional) from input profile

if (isstruct(p_in))
   y_in = p_in.ProfileY;
   z_in = p_in.ProfileZ;
   if (isfield(p_in, 'ProfileS'))
      s_in = p_in.ProfileS;
   else
      s_in = [];
   end
   if (size(y_in,2)>1), y_in = y_in'; end       % make column vectors
   if (size(z_in,2)>1), z_in = z_in'; end
   if (size(s_in,2)>1), s_in = s_in'; end
else
   if (size(p_in,2)>size(p_in,1))
      p_in = p_in';
   end
   y_in = p_in(:,1);
   z_in = p_in(:,2);
   s_in = [];
end

% compute arc-lengths s if not provided in input profile

if (isempty(s_in))
   s_in = make_arclength(y_in, z_in, opt);
end

if (~isstruct(p_in))
   p_in = struct('ProfileY',y_in, 'ProfileZ',z_in, 'ProfileS',s_in, 'is_wheel',is_wheel);
end

% compute weights on data
%  0: uniform
%  1: using w_i = ds_i
%  2: w_i = ds_i * f, with f = f(rho)

if (opt.use_wgt<=0)
   wgt = ones(size(y_in));
   if (~opt.use_bspline)
      wgt([1,end  ]) = wgt([1,end  ]) * 20;     % PP-spline: increase weight at end points
   end
else
   ds_in = diff(s_in);
   ds_in = [ds_in(1); (ds_in(1:end-1)+ds_in(2:end))/2; ds_in(end)];
   wgt   = ds_in;
   if (~opt.use_bspline)
      wgt([1,end  ]) = wgt([1,end  ]) * 20;     % PP-spline: increase weight at end points
   end

   if (opt.use_wgt>=2)
      % estimate curvatures needed for adaptive weighting, using moving window [-sw,sw]

      imeth = 0; s_windw = 4; k_hlf = 4;
      curv = estimate_curv(p_in, imeth, [s_windw, k_hlf]);

      ix = find(abs(curv)>0.1);
      fac = 50; wgt(ix) = fac * wgt(ix);
      ix = find(abs(curv)>0.07 & abs(curv)<=0.1);
      fac = 15; wgt(ix) = fac * wgt(ix);
      ix = find(abs(curv)>0.035 & abs(curv)<=0.07);
      fac = 8; wgt(ix) = fac * wgt(ix);
      ix = find(abs(curv)>0.020 & abs(curv)<=0.035);
      fac = 8; wgt(ix) = fac * wgt(ix);
   end
end

% ikinks=-1: autodetect 'true kinks' between successive segments

if (length(opt.ikinks)~=1 | opt.ikinks>0)
   ikinks_lst = opt.ikinks;
else
   if (~exist('detect_kinks'))
      ikinks_lst = [];
   else
      scale_z    = 1.0;
      make_plot  = 0;
      ikinks_lst = detect_kinks( y_in, z_in, opt.kink_high, opt.kink_low, opt.kink_wid, scale_z, ...
                                                                                is_wheel, make_plot);
   end

   % PP-spline: increase weight at kink points
   if (~isempty(ikinks_lst) & ~opt.use_bspline)
      disp('Increasing weight * 100 at break points')
      wgt(ikinks_lst) = 100 * wgt(ikinks_lst);
   end
end

% create spline representation

spl = make_spline( s_in, y_in, z_in, opt.lambda, wgt, ikinks_lst, opt.iaccel, opt.use_bspline, ...
                                                 opt.ds_bspl, opt.use_deriv, opt.use_repl, opt.idebug );

% evaluate spline at output positions

if (~isempty(opt.ds_sampl))
   is_min = ceil(min(spl.s)/opt.ds_sampl);
   is_max = floor(max(spl.s)/opt.ds_sampl);
   s_out = [is_min : is_max]' * opt.ds_sampl;
   % disp(sprintf('sampling at is=[%d:%d], s=[%5.2f,%5.2f]', is_min, is_max, s_out(1), s_out(end)));
elseif ( opt.use_bspline )
   s_out = s_in;
else
   s_out = spl.s;
end

[ ~, y_out, z_out ] = eval_spline( spl, s_out );

% create profile structure

if (isstruct(p_in))
   p_out = p_in;
else
   p_out = struct();
end
p_out.ProfileS = s_out;
p_out.ProfileY = y_out;
p_out.ProfileZ = z_out;
p_out.spl      = spl;

% estimate surface inclination
% rail:   in [-pi, pi], wheel:  in [0, 2*pi]

alph = atan2( diff(z_out), diff(y_out) );
if (is_wheel)
   ix = find(alph<0); alph(ix) = alph(ix) + 2*pi;
end
alph = [alph(1); 0.5*(alph(1:end-1)+alph(2:end)); alph(end)];
p_out.ProfileAngle = alph;

% reset curvature

p_out.ProfileCurvature = zeros(size(p_out.ProfileS));

end % function smooth_profile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
