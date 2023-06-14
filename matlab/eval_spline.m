
function [ s_out, y_out, z_out ] = eval_spline( spl, s_in, y_in, idebug )

% function [ s_out, y_out, z_out ] = eval_spline( spl, s_in, y_in, idebug )
%
% evaluate parametric spline {s, a0y-a3y, a0z-a3z}
%  - if s_in is given, determine y(s_in), z(s_in)
%  - if y_in is given, determine corresponding s(y) and z(s(y))

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

if (nargin<2)
   s_in = [];
end
if (nargin<3)
   y_in = [];
end
if (nargin<4 | isempty(idebug))
   idebug = 0;
end
if ( (isempty(s_in) & isempty(y_in)) | (~isempty(s_in) & ~isempty(y_in)) )
   disp('ERROR: either s_in or y_in must be given');
   return;
end

if (~isempty(s_in))
   s_out = s_in;
else
   s_out = spline_get_s_at_y( spl, y_in, idebug );
end

if (~isempty(y_in))
   y_out = y_in;
else
   y_out = eval_1d_spline(spl.s, spl.ay0, spl.ay1, spl.ay2, spl.ay3, s_out, idebug);
end

if (nargout>=3)
   z_out = eval_1d_spline(spl.s, spl.az0, spl.az1, spl.az2, spl.az3, s_out, idebug);
end

end % function eval_spline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ s_out ] = spline_get_s_at_y(spl, yev, idebug)

if (nargin<3 | isempty(idebug))
   idebug = 0;
end

if (size(yev,2) > size(yev,1)), yev = yev'; end;

% make input y-values ascending
if (yev(end) > yev(1))
   is_ascending = 1;
else
   yev = flipud(yev);
   is_ascending = 0;
end
if (any(sort(yev)~=yev))
   disp('ERROR: input y_ev should be sorted (ascending or descending)');
   return
end

s  = spl.s;
a0 = spl.ay0;
a1 = spl.ay1;
a2 = spl.ay2;
a3 = spl.ay3;
np = length(a0);

tiny_a2 = 1e-10;
tiny_a3 = 1e-10;
s_out = zeros(size(yev));

for iout = 1 : length(yev)

   % find segment iseg [ a0(i), a0(i+1) ] containing yev

   if (0==1 | ~isfield(spl, 'top_ybrk'))
      iseg = locate_segment( np, a0, yev(iout) );
   else
      npart = length(spl.top_ybrk) - 1;
      itop  = locate_segment( npart+1, spl.top_ybrk, yev(iout) );
      if (itop<=0)
         iseg = 0;
         %disp(sprintf('yev(%2d)=%8.3f lies in top part %d, y=   (-inf,%7.3f]', iout, yev(iout), itop, ...
         %       spl.top_ybrk(itop+1)));
      elseif (itop>=npart+1)
         iseg = np+1;
         %disp(sprintf('yev(%2d)=%8.3f lies in top part %d, y=[%7.3f, +inf)', iout, yev(iout), itop, ...
         %       spl.top_ybrk(end) ));
      else
         %disp(sprintf('yev(%2d)=%8.3f lies in top part %d, y=[%7.3f,%7.3f]', iout, yev(iout), itop, ...
         %       spl.top_ybrk(itop), spl.top_ybrk(itop+1)));
         isec = spl.top_sec(itop);
         ip0 = spl.isec_uni_y(isec);
         ip1 = spl.isec_uni_y(isec+1);
         % disp(sprintf('yev(%2d)=%7.3f lies in section %d, y=[%7.3f,%7.3f]', iout, yev(iout), isec, ...
         %           a0(ip0), a0(ip1)));
         iseg   = ip0-1 + locate_segment( ip1-ip0+1, a0(ip0:ip1), yev(iout) );
      end
   end

   if (idebug>=1)
      if (iseg<=0)
         disp(sprintf('Searching y(%2d)=%8.3f: iseg=%4d, sseg= ( +/- inf,%8.3f], yseg=( +/- inf,%8.3f]',...
                      iout, yev(iout), iseg, s(1), a0(1)));
      elseif (iseg>=np)
         disp(sprintf('Searching y(%2d)=%8.3f: iseg=%4d, sseg= [%8.3f, +/- inf), yseg=[%8.3f, +/- inf)',...
                      iout, yev(iout), iseg, s(np), a0(np)));
      else
         disp(sprintf('Searching y(%2d)=%8.3f: iseg=%4d, sseg= [%8.3f,%8.3f], yseg=[%8.3f,%8.3f]',...
                   iout, yev(iout), iseg, s(iseg), s(iseg+1), a0(iseg), a0(iseg+1)));
      end
   end

   if (isempty(iseg) | iseg<=0)   % before start of spline range: yev < a0(1)

      s_out(iout) = s(1) - 1e-9;  % set s slightly before s(1)

   elseif (iseg>=np)           % after end of spline range: yev > a0(end)

      s_out(iout) = s(end) + 1e-9; % set s slightly after s(end)

   else

      % Solve cubic equation for segment:

      ppcoef = [a0(i), a1(iseg), a2(iseg), a3(iseg)];
      sseg   = s(iseg+[0,1]);
      [sout1, found] = solve_cubic_eq(ppcoef, sseg, yev(iout), idebug, iout, 1, iseg);

      if (found & idebug>=2)
         disp(sprintf('iout = %d: found sl =%7.2f', iout, sout1));
      elseif (~found & idebug>=3)
         disp(sprintf('iout = %d: no solution',iout));
      end

      s_out(iout) = sout1;

   end % yev in interior
end % iout

if (~is_ascending)
   s_out = flipud(s_out);
end

end % function spline_get_s_at_y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ iseg ] = locate_segment( np, vp, vi )

   ascending = (vp(end) >= vp(1));

   if (ascending)
      if (vi < vp(1))
         iseg = 0;
      elseif (vi > vp(end))
         iseg = np + 1;
      else
         % find segment iseg: vp(iseg) <= vi <= vp(iseg+1)
         iseg = find( vp<=vi, 1, 'last');
      end
   else
      if (vi > vp(1))
         iseg = 0;
      elseif (vi < vp(end))
         iseg = np + 1;
      else
         % find segment iseg: vp(iseg) >= vi >= vp(iseg+1)
         iseg = find( vp>=vi, 1, 'last');
      end
   end

end % function locate_segment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ y ] = eval_1d_spline(xsp, a0, a1, a2, a3, xev, idebug)

% evaluate plain (non-parametric) spline y=y(x)

iev_debug = 1;

y = zeros(size(xev));

% copy NaN inputs to output

ix = find(isnan(xev));
y(ix) = NaN;

% assume that xev is sorted

n = length(xsp) - 1;

% all x < x(1): constant extrapolation

i = find(xev<xsp(1));
y(i) = a0(1);

if (idebug>=5 & any(i==iev_debug))
   disp(sprintf('xev(%d) lies before start xsp(1)', iev_debug));
end

% all x(1) <= x <= x(n+1) : cubic polynomials

for isp = 1 : n
   i = find(xev>=xsp(isp) & xev<xsp(isp+1));
   xloc = xev(i) - xsp(isp);
   y(i) = a3(isp) * xloc.^3 + a2(isp) * xloc.^2 + a1(isp) * xloc + a0(isp);

   if (idebug>=5 & any(i==iev_debug))
      j = find(i==iev_debug);
      disp(sprintf('xev(%d) lies in segment isp %d, xloc=%8.3f', iev_debug, isp, xloc(j)));
   end
end

% all x > x(n+1): constant extrapolation

i = find(xev>=xsp(n+1));
xloc = xsp(n+1) - xsp(n);
y(i) = a3(n) * xloc.^3 + a2(n) * xloc.^2 + a1(n) * xloc + a0(n);

end % function eval_1d_spline
