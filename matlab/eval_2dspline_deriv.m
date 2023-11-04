
function [ dx_out, dy_out, dz_out ] = eval_2dspline_deriv( spl2d, u_in, v_in, idir, idebug )

% function [ dx_out, dy_out, dz_out ] = eval_2dspline_deriv( spl2d, u_in, v_in, idir, idebug )
%
% evaluate 1st derivatives of parametric 2d spline {tui, tvj, cij_x, cij_y, cij_z} 
%    for parameter idir = 1 (u) or 2 (v)
%
% the inputs (u,v) are considered a list when of the same size, else a tensor grid

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   if (nargin<4 | isempty(idir))
      idir = 1;
   end
   if (nargin<5 | isempty(idebug))
      idebug = 0;
   end
   if (min(size(u_in))>1 | min(size(v_in))>1)
      disp('ERROR: u_in, v_in must be 1D arrays');
      return;
   end

   has_xij = isfield(spl2d, 'cij_x');

   if (all(size(u_in)==size(v_in)))
      use_list = 1; typ_eval = 'list';
   else
      use_list = 0; typ_eval = 'tens.grid';
   end
   if (size(u_in,2)>size(u_in,1))
      u_in = u_in';     % make column vector
   end
   if (size(v_in,1)>size(v_in,2))
      v_in = v_in';     % make row vector
   end
   k = 4;

   nknotu = length(spl2d.tui);
   nknotv = length(spl2d.tvj);
   nsplu  = nknotu - k;
   nsplv  = nknotv - k;

   if (idebug>=1)
      disp(sprintf('2d spline evaluation (u,v) --> (dx,dy,dz), u: %d x %d, v: %d x %d (%s), idir=%d', ...
                                                                    size(u_in), size(v_in), typ_eval, idir));
   end

   if (idir==1)

      % determine derivatives dx/du, dy/du, dz/du

      %  1. form noutv u-splines at output-locations v_in

      [~, ~, ~, B4] = eval_bspline_basisfnc( spl2d.tvj, v_in, k, idebug );

      % ci_y: [ nsplu, noutv ], B4: [ noutv, nsplv ], cij_y: [ nsplu, nsplv ], v_in: [ noutv ]
      if (has_xij)
         ci_x = (B4 * spl2d.cij_x')';
      end
      ci_y = (B4 * spl2d.cij_y')';
      ci_z = (B4 * spl2d.cij_z')';
      mask = (B4 * spl2d.mask_j')';

      %  2.a collocation matrix for 1st derivative of u-splines at u_in

      [~, ~, B3, B4] = eval_bspline_basisfnc( spl2d.tui, u_in, k, idebug );

      %  2.b backward difference matrix using average step lengths

      dtj_k3 = (spl2d.tui(4:end) - spl2d.tui(1:nknotu-3)) / 3;
      dtj_inv3 = zeros(size(dtj_k3)); ix  = find(dtj_k3>0); dtj_inv3(ix) = 1./ dtj_k3(ix);
      D3     = spdiags(dtj_inv3, 0, nknotu-3, nknotu-4) - spdiags(dtj_inv3(2:end), -1, nknotu-3, nknotu-4);

      %  2. evaluate 1st derivative of u-splines at u_in

      % y_out: [ noutu, noutv ], B3: [ noutu, nsplu+1 ], D3: [nsplu+1, nsplu], 
      %                                                  ci_y: [ nsplu, noutv ], u_in: [ noutu ]

      dy_out = B3 * D3 * ci_y;
      if (has_xij)
         dx_out = B3 * D3 * ci_x;
      else
         dx_out = ones(size(dy_out));    % half-parametric spline: x==u
      end
      dz_out = B3 * D3 * ci_z;
      mask   = B4 * mask;

   else

      % determine derivatives dx/dv, dy/dv, dz/dv

      %  1. form noutu v-splines at output-locations u_in

      [~, ~, ~, B4] = eval_bspline_basisfnc( spl2d.tui, u_in, k, idebug );

      % cj_y: [ noutu, nsplv ], B4: [ noutu, nsplu ], cij_y: [ nsplu, nsplv ], u_in: [ noutu ]

      if (has_xij)
         cj_x = B4 * spl2d.cij_x;
      end
      cj_y = B4 * spl2d.cij_y;
      cj_z = B4 * spl2d.cij_z;
      mask = B4 * spl2d.mask_j;

      %  2.a collocation matrix for 1st derivative of v-splines at v_in

      [~, ~, B3, B4] = eval_bspline_basisfnc( spl2d.tvj, v_in, k, idebug );

      %  2.b backward difference matrix using average step lengths

      dtj_k3 = (spl2d.tvj(4:end) - spl2d.tvj(1:nknotv-3))' / 3;
      dtj_inv3 = zeros(size(dtj_k3)); ix  = find(dtj_k3>0); dtj_inv3(ix) = 1./ dtj_k3(ix);
      D3     = spdiags(dtj_inv3, 0, nknotv-3, nknotv-4) - spdiags(dtj_inv3(2:end), -1, nknotv-3, nknotv-4);

      %  2. evaluate 1st derivative of v-splines at v_in

      % y_out: [ noutu, noutv ], B3: [ noutv, nsplv+1 ], D3: [nsplv+1, nsplv], 
      %                                                  cj_y: [ noutu, nsplv ], v_in: [ noutv ]

      dy_out = (B3 * D3 * cj_y')';
      if (has_xij)
         dx_out = (B3 * D3 * cj_x')';
      else
         dx_out = zeros(size(dy_out));    % half-parametric spline: x==u
      end
      dz_out = (B3 * D3 * cj_z')';
      mask   = (B4 * mask')';

   end

   % remove points that don't have full support

   tiny  = 1e-10;
   ix    = find(mask<1-tiny);
   dx_out(ix) = NaN;
   dy_out(ix) = NaN;
   dz_out(ix) = NaN;

   % output a list when both inputs are of same size

   if (use_list)
      dx_out = diag(dx_out);
      dy_out = diag(dy_out);
      dz_out = diag(dz_out);
      if (any(size(dy_out)~=size(u_in)))
         dx_out = dx_out';
         dy_out = dy_out';
         dz_out = dz_out';
      end
   end

end % function eval_2dspline_deriv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
