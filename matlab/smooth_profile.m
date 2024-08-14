
function [ p_out, spl ] = smooth_profile( p_in, l_filt, lambda, ds_sampl, ikinks, iaccel, ...
                                                                use_wgt, use_bspline, ds_bspl, idebug );

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

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

if (nargin<2)
   l_filt = [];
end
if (nargin<3)
   lambda = [];
end
if (nargin<4)
   ds_sampl = [];
end
if (nargin<5)
   ikinks = -1;
end
if (nargin<6)
   iaccel = [];
end
if (nargin<7 | isempty(use_wgt))
   use_wgt = 1;
end
if (nargin<8 | isempty(use_bspline))
   use_bspline = 1;
end
if (nargin<9 | isempty(ds_bspl))
   ds_bspl = 2;
end
if (nargin<10 | isempty(idebug))
   idebug = 0;
end

is_wheel = 0;
if (isstruct(p_in))
   is_wheel = p_in.is_wheel;
end

if (use_bspline)
   use_deriv = 3;
else
   use_deriv = 2;
end
use_repl = 1;

% compute lambda when l_filt is given

if     ( isempty(l_filt) &  isempty(lambda))
   disp('Either L_filt or lambda must be given.');
   return;
elseif (~isempty(l_filt) & ~isempty(lambda))
   disp('L_filt and lambda may not both be given.');
   return;
elseif (~isempty(l_filt))
   if (use_deriv==3)
      lambda = l_filt^6 / (64*pi^6);
   else
      lambda = l_filt^4 / (16*pi^4);
   end
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
   if (size(y_in,2)>1), y_in = y_in'; end
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
   s_in = make_arclength(y_in, z_in);
end

if (~isstruct(p_in))
   p_in = struct('ProfileY',y_in, 'ProfileZ',z_in, 'ProfileS',s_in, 'is_wheel',is_wheel);
end

% estimate curvatures needed for adaptive weighting, using moving window [-sw,sw]

imeth = 0; s_windw = 4; k_hlf = 4;
curv = estimate_curv(p_in, imeth, [s_windw, k_hlf]);

% compute weights on data
%  0: uniform
%  1: using w_i = ds_i
%  2: w_i = ds_i * f, with f = f(rho)

if (use_wgt<=0)
   wgt = ones(size(y_in));
   if (~use_bspline)
      wgt([1,end  ]) = wgt([1,end  ]) * 20;     % PP-spline: increase weight at end points
   end
else
   ds_in = diff(s_in);
   ds_in = [ds_in(1); (ds_in(1:end-1)+ds_in(2:end))/2; ds_in(end)];
   wgt   = ds_in;
   if (~use_bspline)
      wgt([1,end  ]) = wgt([1,end  ]) * 20;     % PP-spline: increase weight at end points
   end

   if (use_wgt>=2)
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

% ikinks=-1: autodetect kinks between successive segments

if (length(ikinks)==1 & ikinks<0)
   dalph_thrs_high = pi/6; 
   dalph_thrs_low  = dalph_thrs_high / 5;
   dst_max         = 0.0;
   scale_z         = 1.0;
   is_wheel        = [];
   make_plot       = 0;
   ikinks = detect_kinks( s_in, y_in, z_in, dalph_thrs_high, dalph_thrs_low, dst_max, scale_z, ...
                                                                                is_wheel, make_plot);

   % PP-spline: increase weight at kink points
   if (~isempty(ikinks) & ~use_bspline)
      disp('Increasing weight * 100 at break points')
      wgt(ikinks) = 100 * wgt(ikinks);
   end
end

% create spline representation

spl = make_spline( s_in, y_in, z_in, lambda, wgt, ikinks, iaccel, use_bspline, ds_bspl, ...
                                                                  use_deriv, use_repl, idebug );

% evaluate spline at output positions

if (~isempty(ds_sampl))
   is_min = ceil(min(spl.s)/ds_sampl);
   is_max = floor(max(spl.s)/ds_sampl);
   s_out = [is_min : is_max]' * ds_sampl;
   % disp(sprintf('sampling at is=[%d:%d], s=[%5.2f,%5.2f]', is_min, is_max, s_out(1), s_out(end)));
elseif ( use_bspline )
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
