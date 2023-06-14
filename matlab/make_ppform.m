
function [ spl ] = make_ppform( spl, D1, D2, D3 )

% function [ spl ] = make_ppform( spl, [D1], [D2], [D3] )
%
% determine PP-form: value of spline & derivatives at start of each segment

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   % compute D1, D2, D3 if not provided

   if (nargin<4 | isempty(D1) | isempty(D2) | isempty(D3))

      % average step sizes over 1/2/3 adjacent intervals

      nknot  = length(spl.tj);
      dtj_k1 =  spl.tj(2:end) - spl.tj(1:nknot-1);
      dtj_k2 = (spl.tj(3:end) - spl.tj(1:nknot-2)) / 2;
      dtj_k3 = (spl.tj(4:end) - spl.tj(1:nknot-3)) / 3;

      % inverse average step sizes, zero at zero step size

      dtj_inv1 = zeros(size(dtj_k1)); ix  = find(dtj_k1>0); dtj_inv1(ix) = 1./ dtj_k1(ix);
      dtj_inv2 = zeros(size(dtj_k2)); ix  = find(dtj_k2>0); dtj_inv2(ix) = 1./ dtj_k2(ix);
      dtj_inv3 = zeros(size(dtj_k3)); ix  = find(dtj_k3>0); dtj_inv3(ix) = 1./ dtj_k3(ix);

      % backward difference matrices using average step lengths 1/2/3

      D1   = spdiags(dtj_inv1, 0, nknot-1, nknot-2) - spdiags(dtj_inv1(2:end), -1, nknot-1, nknot-2);
      D2   = spdiags(dtj_inv2, 0, nknot-2, nknot-3) - spdiags(dtj_inv2(2:end), -1, nknot-2, nknot-3);
      D3   = spdiags(dtj_inv3, 0, nknot-3, nknot-4) - spdiags(dtj_inv3(2:end), -1, nknot-3, nknot-4);
   end

   % determine value of each basis-function at the start of each segment ksi_i

   s = max(spl.s, spl.tj(4));      % avoid warning for extended knots
   s = min(    s, spl.tj(end-3));  %  --> constant extrapolation outside basic interval
   [B1, B2, B3, B4] = eval_bspline_basisfnc( spl.tj, s );

   % determine PP-form: value of spline & derivatives at start of each segment

   spl.ay0 = B4 * spl.cy;
   spl.ay1 = B3 * D3 * spl.cy;
   spl.ay2 = B2 * D2 * D3 * spl.cy / 2;
   spl.ay3 = B1 * D1 * D2 * D3 * spl.cy / 6;
   spl.az0 = B4 * spl.cz;
   spl.az1 = B3 * D3 * spl.cz;
   spl.az2 = B2 * D2 * D3 * spl.cz / 2;
   spl.az3 = B1 * D1 * D2 * D3 * spl.cz / 6;

end % function make_ppform

