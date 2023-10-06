
function [ dy_out, ddy_out, dz_out, ddz_out, kappa, rcurv ] = ...
                                eval_spline_deriv( spl, s_out, idebug )

% function [ dy_out, ddy_out, dz_out, ddz_out, kappa, rcurv ] = ...
%                               eval_spline_deriv( spl, s_out, idebug )
%
% evaluate derivatives for parametric spline {s, a0y-a3y, a0z-a3z}
%  - determine dy/ds, d2y/ds2, dz/ds and d2z/ds2

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

if (nargin<2 | isempty(s_out))
   s_out = spl.s;
end
if (nargin<3 | isempty(idebug))
   idebug = 0;
end

[dy_out, ddy_out] = eval_1d_spline_deriv(spl.s, spl.ay0, spl.ay1, ...
                                                spl.ay2, spl.ay3, s_out);

[dz_out, ddz_out] = eval_1d_spline_deriv(spl.s, spl.az0, spl.az1, ...
                                                spl.az2, spl.az3, s_out);

kappa = (dy_out.*ddz_out - dz_out.*ddy_out) ./ (dy_out.^2 + dz_out.^2).^(3/2);

k_min = 1e-6;
rcurv = zeros(size(kappa));
ix = find(abs(kappa)>=k_min); rcurv(ix) = 1 ./ kappa(ix);
ix = find(abs(kappa)< k_min); rcurv(ix) = sign(kappa(ix)) / k_min;

end % function eval_spline_deriv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ dy, ddy ] = eval_1d_spline_deriv(xsp, a0, a1, a2, a3, xev)

% evaluate plain (non-parametric) spline y=y(x)

dy  = zeros(size(xev));
ddy = zeros(size(xev));

% assume that xev is sorted

n = length(xsp) - 1;

% all x < x(1): constant extrapolation

i = find(xev<xsp(1));
dy(i)  = a1(1);
ddy(i) = a2(1);

% all x(1) <= x <= x(n+1) : cubic polynomials

for isp = 1 : n
   i = find(xev>=xsp(isp) & xev<xsp(isp+1));
   xloc = xev(i) - xsp(isp);
   dy(i)  = 3 * a3(isp) * xloc.^2 + 2 * a2(isp) * xloc    + a1(isp);
   ddy(i) = 6 * a3(isp) * xloc    + 2 * a2(isp);
end

% all x > x(n+1): constant extrapolation

i = find(xev>=xsp(n+1));
xloc = xsp(n+1) - xsp(n);
dy(i)  = 3 * a3(n) * xloc.^2 + 2 * a2(n) * xloc    + a1(n);
ddy(i) = 6 * a3(n) * xloc    + 2 * a2(n);

end % function eval_1d_spline_deriv

