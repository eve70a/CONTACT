
function [ dst_max, l, n ] = plot_update(prf1, prf2, fac, plot_n, use_angl2, norm_dist, show_angl, ...
                                                        add_plot, force_spline)

% function [ dst_max, l, n ] = plot_update(prf1, prf2, fac, plot_n, use_angl2, norm_dist, show_angl,
%                                                       add_plot, force_spline)
%
% show difference of two profiles at magnification 'fac'
%
% fac       == magnification factor, negative/default: automatic setting
% plot_n    == connect profile and update with normal lines at points [1 : plot_n : end]
% use_angl2 == use surface angle of second profile
% norm_dist == plot tangential updates (<0), full updates (0) or normal updates (>0)
% show_angl == add markers at specified surface inclinations [ show_angl ] (deg)
% add_plot  == clear current figure (0) or add plot to figure (1)
% force_spline == assume different s parameterization even if #points is equal

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

if (nargin<3 | isempty(fac))
   fac       = -1;
end
if (nargin<4 | isempty(plot_n))
   plot_n    = 5;
end
if (nargin<5 | isempty(use_angl2))
   use_angl2 = 0;
end
if (nargin<6 | isempty(norm_dist))
   norm_dist = 1;
end
if (nargin<7);
   show_angl = [];
end
if (nargin<8 | isempty(add_plot))
   add_plot  = 0;
end
if (nargin<9 | isempty(force_spline))
   force_spline = 0;
end

% get (y,z)-data for the two profiles

y1_in = prf1.ProfileY;
z1_in = prf1.ProfileZ;
n1    = length(y1_in);
if (min(size(y1_in))>1)
   disp('ERROR in plot_update: y1_in should be an n x 1 vector');
   return;
end
if (size(y1_in,1)<size(y1_in,2))
   y1_in = y1_in'; z1_in = z1_in';  % make column vectors
end

y2    = prf2.ProfileY;
z2    = prf2.ProfileZ;
n2    = length(y2);
if (min(size(y2))>1)
   disp('ERROR in plot_update: y2 should be an n x 1 vector');
   return;
end
if (size(y2,1)<size(y2,2))
   y2 = y2'; z2 = z2';              % make column vectors
end

% different #points: resample using spline interpolation

if (n1~=n2 | force_spline)
   % form spline for profile 1
   lambda = 0; use_bspline = 0; wgt = []; ikinks = []; iaccel = [];
   spl  = make_spline([], y1_in, z1_in, lambda, wgt, ikinks, iaccel, use_bspline);

   % determine coordinates s1 in spline 1 for points (y2,z2)
   % figure(5); clf; hold on;
   % plot(y1_in, z1_in, '-o');
   % plot(y2(774:776), z2(774:776), '*');
   % set(gca,'ydir','reverse');
   % grid on;
   % axis equal;
   % axis([-0.5 0.5 1.4 2.2]);

   idebug = -2; ifig = -5;
   si   = map_pts_to_spline(spl, y2, z2, idebug, ifig);

   % determine (y1,z1) at s1-positions
   [~, y1_int, z1_int] = eval_spline(spl, si);
else
   y1_int = y1_in;
   z1_int = z1_in;
end

% tangent and inward normal on profile surface at positions used for prf2

if (use_angl2)
   dy = diff(y2);
   dz = diff(z2);
else
   dy = diff(y1_int);
   dz = diff(z1_int);
end
ds = sqrt(dy.^2 + dz.^2);

ty = dy ./ ds;
tz = dz ./ ds;
ny = -tz;
nz =  ty;

% move normals from segments to end-points

ty = [ty(1); (ty(1:end-1)+ty(2:end))/2; ty(end)];
tz = [tz(1); (tz(1:end-1)+tz(2:end))/2; tz(end)];
ny = [ny(1); (ny(1:end-1)+ny(2:end))/2; ny(end)];
nz = [nz(1); (nz(1:end-1)+nz(2:end))/2; nz(end)];
alph = atan2(tz, ty);

% dst_n, dst_t: signed distance in normal/tangential direction between input profiles

dst_y = (y2 - y1_int);
dst_z = (z2 - z1_int);
dst_n = ny .* dst_y + nz .* dst_z;
dst_t = ty .* dst_y + tz .* dst_z;

if (norm_dist>0)
   [dst_max, ixm] = max(abs(dst_n));
elseif (norm_dist==0)
   [dst_max, ixm] = max(dst_y.^2 + dst_z.^2);
else
   [dst_max, ixm] = max(abs(dst_t));
end

% fac<0: automatic scaling factor, fac * dst_max ~=~ 4% to 8% * tot_width

if (fac<0)
   fac_list  = [ -1e-6 0.00005 0.0001 0.00025 0.0005 0.001 0.0025 0.005 0.01 0.025 0.05 0.1 0.25 0.5 1.0;
                  1000   20000  10000    4000   2000  1000    400   200  100    40   20  10    4   2   1];
   tot_width = max(y1_in) - min(y1_in);
   ifac      = find(fac_list(1,:)<25*dst_max/tot_width, 1, 'last');
   fac       = fac_list(2,ifac);
end

% plot original profiles 1 + 2

if (~add_plot)
   clf; hold on;
end
set(gca,'ydir','reverse');
axis equal;
grid on;

l(1) = plot(y1_in, z1_in);
l(2) = plot(y2, z2, '--');

% plot magnified updates 2-1 and maximum update

if (norm_dist>0)        % updates in normal direction
   l(3) = plot(y1_int+fac*dst_n.*ny, z1_int+fac*dst_n.*nz);
   plot(y1_int(ixm)+fac*dst_n(ixm)*ny(ixm), z1_int(ixm)+fac*dst_n(ixm)*nz(ixm), '*');
elseif (norm_dist==0)   % full updates normal + tangential
   l(3) = plot(y1_int+fac*dst_y, z1_int+fac*dst_z);
   plot(y1_int(ixm)+fac*dst_y(ixm), z1_int(ixm)+fac*dst_z(ixm), '*');
else                    % updates in tangential direction, normal to curve
   l(3) = plot(y1_int+fac*dst_t.*ny, z1_int+fac*dst_t.*nz);
   plot(y1_int(ixm)+fac*dst_t(ixm)*ny(ixm), z1_int(ixm)+fac*dst_t(ixm)*nz(ixm), '*');
end

% connect 1 with magnified updates at each profile point i

n = []; j = 0;
if (plot_n)
   if (norm_dist>0)
      for i = [1 : plot_n : length(dst_n)]
         if (fac*abs(dst_n(i))>0.01)
            j = j + 1;
            n(j) = plot(y1_int(i)+[0 fac*dst_n(i)*ny(i)], z1_int(i)+[0 fac*dst_n(i)*nz(i)], ...
                           'color',matlab_color(6), 'linewidth',1);
         end
      end
   elseif (norm_dist==0)
      dst = sqrt(dst_y.^2+dst_z.^2);
      for i = [1 : plot_n : length(dst)]
         if (fac*dst(i)>0.01)
            j = j + 1;
            n(j) = plot(y1_int(i)+[0 fac*dst_y(i)], z1_int(i)+[0 fac*dst_z(i)], ...
                           'color',matlab_color(4), 'linewidth',1);
         end
      end
   else
      for i = [1 : plot_n : length(dst_t)]
         if (fac*abs(dst_t(i))>0.01)
            j = j + 1;
            n(j) = plot(y1_int(i)+[0 fac*dst_t(i)*ny(i)], z1_int(i)+[0 fac*dst_t(i)*nz(i)], ...
                           'color',matlab_color(5), 'linewidth',1);
         end
      end
   end
end

% pull smoothed profile to foreground

if (1==1)
   c = get(gca,'children');
   c = [c(end-1); c(1:end-2); c(end)];
   set(gca,'children',c);
end

% plot angle indicators

if (~isempty(show_angl))
   if (size(show_angl,1)>size(show_angl,2)), show_angl = show_angl'; end
   alph_deg = alph * 180/pi;
   for angl = show_angl
      % disp(sprintf('Show angle marker at %4d deg', angl));
      % take profile points at start of segments with zero crossing
      ix = find( (alph_deg(1:end-1)-angl).*(alph_deg(2:end)-angl) <= 0 );
      if (~isempty(ix))
         % disp(sprintf('  found at ix=%4d %4d %4d %4d',ix));
         for j = 1 : length(ix)
            yi = y1_int(ix(j));
            zi = z1_int(ix(j));
            nvec = [-sin(alph(ix(j))), cos(alph(ix(j)))]; % inward normal
            plot(yi+5*nvec(1)*[-1,1], zi+5*nvec(2)*[-1,1], '-','color',matlab_color(5),'linewidth',1);
            text(yi-5*nvec(1), zi-5*nvec(2), sprintf('$%d^\\circ$',angl), 'interpreter','latex', ...
                                'fontsize',10, 'horizontalalignment','center', 'verticalalignment','bottom');
         end
      end
   end
end

% print biggest updates on flange

if (0==1)
   ix = find( (y1_int>5 & y1_int<40) | (y1_int>125 & y1_int<138));
   for j = ix'
      if (abs(dst_n(j))==max(abs(dst_n([j-15:j+15]))))
         dst = fac * dst_n(j) + 1 * sign(dst_n(j));
         text(y1_int(j)+dst*ny(j), z1_int(j)+dst*nz(j), ...
                sprintf('%5.3f',dst_n(j)), 'horizontalalignment','center', 'verticalalignment','middle');
      end
   end
end

% show magnification factor

if (1==1) % northeast / northwest
   text(0.95, 0.93, sprintf('$%d\\times$',round(fac)), 'units','normalized', ...
                                        'horizontalalignment','right', 'interpreter','latex');
   text(0.05, 0.93, sprintf('max update %5.3f mm',dst_max), 'units','normalized', ...
                                        'horizontalalignment','left', 'interpreter','latex');
else % southeast / southwest
   text(0.95, 0.07, sprintf('$%d\\times$',round(fac)), 'units','normalized', ...
                                        'horizontalalignment','right', 'interpreter','latex');
   text(0.05, 0.07, sprintf('max update %5.3f mm',dst_max), 'units','normalized', ...
                                        'horizontalalignment','left', 'interpreter','latex');
end

