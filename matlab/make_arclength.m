
function [ s ] = make_arclength( y, z, opt )

% function [ s ] = make_arclength( y, z, opt )
%
% compute summed chord length parameterization { s_i } of curve { (y_i, z_i) }:
% cumulative sum of distances ds between successive points
% with correction for fluctuations wrt locally average heading
%
% opt = struct('cos_thrs',0.5, 'frac_ds_min',0.1, 'kmax',10);

% Copyright 2008-2026 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

if (nargin<3 | ~isstruct(opt))
   opt = struct;
end
if (~isfield(opt, 'cos_thrs'))
   opt.cos_thrs = 0.5;
end
if (~isfield(opt, 'frac_ds_min'))
   opt.frac_ds_min = 0.1;
end
if (~isfield(opt, 'kmax'))
   opt.kmax        =  10;
end

% convert to column vectors

if (size(y,2)>size(y,1)), y = y'; end
if (size(z,2)>size(z,1)), z = z'; end

% compute summed chord length without correction

dy = diff(y);
dz = diff(z);
ds = sqrt(dy.^2 + dz.^2);
s  = [0; cumsum(ds) ];

% check for larger angles dalpha between successive segments

cos_da = (dy(1:end-1) .* dy(2:end) + dz(1:end-1) .* dz(2:end)) ./ (ds(1:end-1) .* ds(2:end));
ix     = find(cos_da<opt.cos_thrs);

if (0==1)
   figure(15); clf; hold on;
   plot(cos_da);
   plot([0, length(dy)], opt.cos_thrs*[1 1], '--', 'color',matlab_color(3));
   grid on;
   xlabel('point $i$', 'interpreter','latex');
   ylabel('$\cos\delta\alpha_i$', 'interpreter','latex');
   title('Cosine of angles $\delta\alpha$ between segments', 'interpreter','latex');
end

% for wiggly signal, correct lengths using the average heading

if (~isempty(ix))
   % disp(sprintf('%d segments with change dalpha > %3.1f deg', length(ix), acos(opt.cos_thrs)*180/pi ));

   % determine average heading (ty_avg,tz_avg) per segment

   [ ty_avg, tz_avg ] = avg_heading( y, z, opt.kmax );

   % compute angle davg between segment and average heading

   cos_tavg = (dy .* ty_avg + dz .* tz_avg) ./ ds;

   % scale segment lengths ds to the length in the average heading direction

   ds = ds .* max(opt.frac_ds_min, cos_tavg);
   s  = [0; cumsum(ds) ];
end

% shift s-values to get s = 0 at y = 0
%  - interpolate [y, s] to y=0
%  - shift s := s - s0

is_monotonic = (all(diff(y)>0) | all(diff(y)<0));
ix_zeros     = find(abs(y)<1e-20);
if (isempty(ix_zeros))
   ix_zeros  = find(y(1:end-1).*y(2:end)<0);
end

if (is_monotonic & min(y)<0 & max(y)>0)
   sref = interp1(y, s, 0);
elseif (length(ix_zeros)==1)
   if (ix_zeros<length(y))
      ix   = ix_zeros + [ 0,1];
   else
      ix   = ix_zeros + [-1,0];
   end
   sref = interp1(y(ix), s(ix), 0);
else
   sref = 0;
end
s  = s - sref;

end % function make_arclength

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ ty_avg, tz_avg, y_avg, z_avg ] = avg_heading( y, z, kmax )

% function [ ty_avg, tz_avg, y_avg, z_avg ] = avg_heading( y, z, kmax )
% 
% for each segment i, determine average heading (ty_avg, tz_avg)
%
% using orthogonal projection to fit the best line through [ y(i+k), z(i+k) ], k = [-kmax:kmax+1]
% with center (y_avg,z_avg) and tangent (ty_avg,tz_avg)

n      = length(y)-1;
ty_avg = zeros(n,1);
tz_avg = zeros(n,1);

if (n<2*kmax+3)
   disp(sprintf('avg_heading: cannot use kmax=%d for n=%d segments', kmax, n));
   return;
end

% determine initial estimate: heading for window:  x(i+k+1) - x(i-k)

y_sta = y(max(  1,[1:n]-kmax));
y_end = y(min(n+1,[1:n]+1+kmax));
z_sta = z(max(  1,[1:n]-kmax));
z_end = z(min(n+1,[1:n]+1+kmax));
dy    = y_end - y_sta;
dz    = z_end - z_sta;
ds    = sqrt(dy.^2 + dz.^2);
ty_avg = dy ./ ds;
tz_avg = dz ./ ds;
y_avg  = y;
z_avg  = z;

% improve estimate using orthogonal projection

use_fast_meth = 1;

if (~use_fast_meth)
   % slow implementation showing full algorithm
 
   y_avg  = zeros(n,1);
   z_avg  = zeros(n,1);
   Syy  = zeros(n,1);
   Syz  = zeros(n,1);
   Szz  = zeros(n,1);

   for i = 1 : n
      j        = [ max(1,i-kmax) : min(n+1,i+kmax+1) ];
      y_avg(i) = mean( y(j) );
      z_avg(i) = mean( z(j) );

      dy       = y(j) - y_avg(i);
      dz       = z(j) - z_avg(i);
      Syy(i)   = mean( dy.^2 );  % == mean(y(j).^2) - y_avg(i)^2;
      Syz(i)   = mean( dy.*dz ); % == mean(y(j)*z(j)) - y_avg(i)*z_avg(i);
      Szz(i)   = mean( dz.^2 );  % == mean(z(j).^2) - z_avg(i)^2;
   end

else
   % fast algorithm using vector operations

   y_avg = conv( y(1:n), ones(2*kmax+2,1)/(2*kmax+2), 'same');
   z_avg = conv( z(1:n), ones(2*kmax+2,1)/(2*kmax+2), 'same');
   Syy   = conv( y(1:n).^2, ones(2*kmax+2,1)/(2*kmax+2), 'same') - y_avg.^2;
   Szz   = conv( z(1:n).^2, ones(2*kmax+2,1)/(2*kmax+2), 'same') - z_avg.^2;
   Syz   = conv( y(1:n).*z(1:n), ones(2*kmax+2,1)/(2*kmax+2), 'same') - y_avg.*z_avg;
   for i = 1 : kmax,
      y_avg(i) = mean(y(1:i+kmax+1));
      z_avg(i) = mean(z(1:i+kmax+1));
      Syy(i)   = mean(y(1:i+kmax+1).^2) - y_avg(i)^2;
      Szz(i)   = mean(z(1:i+kmax+1).^2) - z_avg(i)^2;
      Syz(i)   = mean(y(1:i+kmax+1).*z(1:i+kmax+1)) - y_avg(i)*z_avg(i);
   end
   for i = n-kmax : n
      y_avg(i) = mean(y(i-kmax:n+1));
      z_avg(i) = mean(z(i-kmax:n+1));
      Syy(i)   = mean(y(i-kmax:n+1).^2) - y_avg(i)^2;
      Szz(i)   = mean(z(i-kmax:n+1).^2) - z_avg(i)^2;
      Syz(i)   = mean(y(i-kmax:n+1).*z(i-kmax:n+1)) - y_avg(i)*z_avg(i);
   end
end

% average heading alph_avg:

alph_avg = 0.5 * atan2( 2*Syz, Syy-Szz );
ty_orth  = cos(alph_avg);
tz_orth  = sin(alph_avg);

% flip sign where t_orth goes against first estimate

ix = find( ty_orth.*ty_avg + tz_orth.*tz_avg < 0 );
ty_orth(ix) = -ty_orth(ix);
tz_orth(ix) = -tz_orth(ix);

% overwrite initial estimate

ty_avg = ty_orth;
tz_avg = tz_orth;

end % function avg_heading
