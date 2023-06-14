
function [ u_out, y_out, z_out ] = eval_2dspline( spl2d, x_in, u_in, y_in, idebug )

% function [ u_out, y_out, z_out ] = eval_2dspline( spl2d, x_in, u_in, y_in, idebug )
%
% evaluate parametric 2d spline {txi, tuj, cij_y, cij_z} 
%  - if u_in is given, determine [yz]_out(x_in, u_in) for given (x_in, u_in)
%  - if y_in is given, determine u_out and corresponding z_out(x_in, u_out)
%
% the inputs (x,u) or (x,y) are considered a list when of the same size, else a tensor grid

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   if (nargin<3)
      u_in = [];
   end
   if (nargin<4)
      y_in = [];
   end
   if (nargin<5 | isempty(idebug))
      idebug = 0;
   end
   if ( (isempty(u_in) & isempty(y_in)) | (~isempty(u_in) & ~isempty(y_in)) )
      disp('ERROR: either u_in or y_in must be given');
      return;
   end

   if (~isempty(u_in))

      % forward calculation u_in --> y_out, z_out

      % disp('u_in specified, forward (x,u) --> (y,z)')
      u_out = u_in;
      [ y_out, z_out ] = eval_2dspline_forward(spl2d, x_in, u_in, idebug);

   else

      % inverse calculation y_in --> u_out

      % disp('y_in specified, inverse (x,y) --> (u)')
      u_out = eval_2dspline_inverse( spl2d, x_in, y_in, idebug );

      % forward calculation u_out --> z_out

      % disp('y_in specified, forward (x,u) --> (z)')
      y_out = y_in;
      [ ~, z_out ] = eval_2dspline_forward(spl2d, x_in, u_out, idebug);

   end

end % function eval_2dspline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y_out, z_out] = eval_2dspline_forward(spl2d, x_in, u_in, idebug)

% function [y_out, z_out] = eval_2dspline_forward(spl2d, x_in, u_in, idebug)
%
% if x_in and u_in have the same size, evaluate 2d spline at pairs (x_in(i), u_in(i))
%                                else, evaluate 2d spline at tensor grid x_in x u_in.

   if (all(size(x_in)==size(u_in)))
      use_list = 1; typ_eval = 'list';
   else
      use_list = 0; typ_eval = 'tens.grid';
   end
   k = 4;

   if (idebug>=1)
      disp(sprintf('Forward spline evaluation (x,u) --> (y,z), x: %d x %d, u: %d x %d (%s)', ...
                                                                        size(x_in), size(u_in), typ_eval));
   end

   % determine collocation matrix for u-direction: evaluate each B-spline at each output location
   %  - check basic interval [tj(4),tj(nknot-k+1)]

   n0 = nnz(u_in < spl2d.tuj(4));
   n1 = nnz(u_in > spl2d.tuj(end-3));
   if (n0+n1>0)
      disp(sprintf('Warning: there are %d u-positions before and %d after spline u-range [%3.1f,%3.1f]', ...
                n0, n1, spl2d.tuj(4), spl2d.tuj(end-3)));
   end

   numnan = nnz(isnan(u_in));
   if (idebug>=3 | (idebug>=1 & numnan>0))
      disp(sprintf('array u_in has %d NaN-values', nnz(isnan(u_in))));
   end

   [~, ~, ~, Bmat] = eval_bspline_basisfnc( spl2d.tuj, u_in, k );

   % matrix-multply to get 1d spline coefficients for x-direction

   % ci_y: [ nsplx, nout ], Bmat: [ nout, nsplu ], cij_y: [ nsplx, nsplu ], u_in: [ nout ]
   ci_y = (Bmat * spl2d.cij_y')';
   ci_z = (Bmat * spl2d.cij_z')';
   mask = (Bmat * spl2d.mask_j')';

   % determine collocation matrix for x-direction: evaluate each B-spline at each output location

   n0 = nnz(x_in < spl2d.txi(4));
   n1 = nnz(x_in > spl2d.txi(end-3));
   if (n0+n1>0)
      disp(sprintf('Warning: there are %d x-positions before and %d after spline x-range [%3.1f,%3.1f]', ...
                n0, n1, spl2d.txi(4), spl2d.txi(end-3)));
   end

   [~, ~, ~, Bmat] = eval_bspline_basisfnc( spl2d.txi, x_in, k );

   % matrix-multply to get output values

   % y_out: [ noutx, noutu ], Bmat: [ noutx, nsplx ], ci_y: [ nsplx, noutu ], x_in: [ noutx ]
   y_out = Bmat * ci_y;
   z_out = Bmat * ci_z;
   mask  = Bmat * mask;

   % remove points that don't have full support

   tiny  = 1e-10;
   ix    = find(mask<1-tiny);
   y_out(ix) = NaN;
   z_out(ix) = NaN;

   % output a list when both inputs are of same size

   if (all(size(x_in)==size(u_in)))
      y_out = diag(y_out);
      z_out = diag(z_out);
      if (any(size(y_out)~=size(x_in)))
         y_out = y_out';
         z_out = z_out';
      end
   end

end % function eval_2dspline_forward

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ u_out ] = eval_2dspline_inverse(spl2d, x_in, y_in, idebug)

% function [ u_out ] = eval_2dspline_inverse(spl2d, x_in, y_in, idebug)
%
% if x_in and y_in have the same size, invert 2d spline at pairs (x_in(i), y_in(i))
%                                else, invert 2d spline at tensor grid x_in x y_in.

   tiny   = 1e-10 * max(abs(spl2d.xi));
   k      = 4;
   nknotx = length(spl2d.txi);
   nknotu = length(spl2d.tuj);
   nsplx  = nknotx - k;
   nsplu  = nknotu - k;

   if (all(size(x_in)==size(y_in)))
      use_list = 1; typ_eval = 'list';
      nx    = length(x_in);
      ny    = 1;
      u_out = ones(size(x_in));
   else
      use_list = 0; typ_eval = 'tens.grid';
      nx    = length(x_in);
      ny    = length(y_in);
      u_out = ones(nx, ny);
   end

   if (idebug>=1)
      disp(sprintf('Inverse spline evaluation (x,y) --> (u),   x: %d x %d, y: %d x %d (%s)', ...
                                                                    size(x_in), size(y_in), typ_eval));
   end

   % check basic interval in x-direction, knots [t(4), t(nknot-k+1)]

   if (any(x_in<spl2d.txi(k)) | any(x_in>spl2d.txi(nknotx-k+1)))
      disp(sprintf('ERROR: all input x_in must lie in range of spline x=[%3.1e,%3.1e]', ...
                spl2d.txi(k), spl2d.txi(nknotx-k+1)));
      return;
   end

   % determine collocation matrix for computing coefficients cj_y

   [~, ~, ~, Bmatx] = eval_bspline_basisfnc( spl2d.txi, x_in );

   % loop over output positions x_in

   for iout = 1 : nx

      % determine spline coefficients cj_y at position x_in(iout)

      % cj_y: [ nsplu, 1 ], Bmat: [ nx, nsplx ], cij_y: [ nsplx, nsplu ], x_in: [ nx ]

      cj_y = (Bmatx(iout,:) * spl2d.cij_y)';

      % n+k knots gives n basisfunctions and n spline coefficients; 
      % basisfunction j becomes nonzero at knot j until knot j+k.
      % values of y in segment j = [t_j,t_j+1) are defined by spline coefficients j-k+1:j
      % n+k knots gives n+k-1 segments. The first k-1 and last k-1 segments are outside the basic interval
      % this gives n-k+1 segments in the basic interval with numbers k to n

      % determine interval [ylow, yhig] for segments j, j = k to n

      tmp            = NaN * ones(4, nknotu);
      tmp(:,k:nsplu) = [ cj_y(1:nsplu-3)';  cj_y(2:nsplu-2)';  cj_y(3:nsplu-1)';  cj_y(4:nsplu)' ];

      ylow = min(tmp); yhig = max(tmp);

      % no inner loop when using list-input, size(x_in) == size(y_in)

      if (use_list)
         j0 = iout; j1 = iout;
      else
         j0 = 1; j1 = ny;
      end

      for jout = j0 : j1

         % determine segments jseg that may contain y(jout) at xout

         y    = y_in(jout);
         jseg = find( ylow(1:nsplu)<=y & y<=yhig(1:nsplu) );
         nseg = length(jseg);

         if (nseg<=0)
            % disp(sprintf('ERROR: no segments for (x,y) = (%6.1f,%6.1f)', x_in(iout), y));
            found = 0;
         else
            if (idebug>=2)
               disp(sprintf(['(i,j)=(%2d,%2d), (x,y)=(%6.1f,%6.1f): found %d possible segments jseg = ', ...
                   '%d, %d, %d, %d, %d, %d'], iout, jout, x_in(iout), y_in(jout), nseg, jseg(1:min(6,nseg))));
               if (idebug>=-4)

                  [~, ~, ~, Bmatu] = eval_bspline_basisfnc( spl2d.tuj, spl2d.uj );
                  % ytmp: [ nu, 1 ], Bmatu: [ nu, nu ], cj_y: [ nu, 1 ]
                  ytmp = Bmatu * cj_y;
                  tuj  = spl2d.tuj;
                  tgrev(1:nsplu) = (tuj(2:nsplu+1) + tuj(3:nsplu+2) + tuj(4:nsplu+3)) / 3;
                  tm   = (tgrev(1:nsplu-1) + tgrev(2:end)) / 2;

                  figure(7); clf; hold on;
                  tmpt  = reshape( [1;1]*tuj(1:nsplu+1), 1, 2*nsplu+2 );
                  tmpyl = reshape( [1;1]*ylow(1:nsplu), 1, 2*nsplu );
                  tmpyh = reshape( [1;1]*yhig(1:nsplu), 1, 2*nsplu );
                  plot(tmpt(2:end-1), tmpyl);
                  plot(tmpt(2:end-1), tmpyh);
                  plot(spl2d.uj([1,end]), y*[1 1], '--');
                  plot(spl2d.uj, ytmp);
                  plot([1;1]*(tuj(jseg)+tuj(jseg+1))/2, [ylow(jseg);yhig(jseg)], 'color',matlab_color(3))
                  plot(spl2d.tuj(nsplu+1), cj_y(end), '*');
                  for jj = jseg
                      text( (tuj(jj)+tuj(jj+1))/2, ylow(jj), sprintf('j=%d',jj), ...
                        'horizontalalignment','center', 'verticalalignment','top');
                  end
                  grid on
                  xlabel('u_j'); ylabel('y_r');
                  legend('y_{low}(u;x)', 'y_{hig}(u;x)', 'target y_{out}', 'actual y(u), x=x_{out}', ...
                                                                                'location','northwest')
               end
            end

            % try potential segments one-by-one until a solution is found

            jj = 1;
            found = 0;
            while(jj<=nseg & ~found)

               % spline segment: interval [u_j, u_{j+1}]

               useg = spl2d.tuj( jseg(jj)+[0,1] );

               % PP-spline coefficients for segment:

               coef  = get_ppcoef_1seg_sparse(spl2d, cj_y, jseg(jj));
               if (0==1)
                  coef0 = get_ppcoef_1seg_full(spl2d, cj_y, jseg(jj));
                  coef1 = get_ppcoef_1seg_sparse(spl2d, cj_y, jseg(jj));
                  if (max(abs(coef0-coef1)>1e-10))
                     disp('coef:')
                     disp([coef0; coef1])
                  end
               end

               % Solve cubic equation for segment:

               [u, found] = solve_cubic_eq(coef, useg, y, idebug, iout, jout, jseg(jj));
               if (found & idebug>=2)
                  disp(sprintf('jj = %d: found u =%7.2f', jj, u));
               elseif (~found & idebug>=3)
                  disp(sprintf('jj = %d: no solution',jj));
               end

               if (~found), jj = jj + 1; end
            end % jj

         end % nseg>0

         % store solution

         if (~found), 
            if (idebug>=2)
               disp(sprintf('No solution for (i,j)=(%d,%d), setting NaN', iout, jout));
            end
            u = NaN;
         end
         if (use_list)
            u_out(iout) = u;
         else
            u_out(iout,jout) = u;
         end

      end % for jout
   end % for iout

end % function eval_2dspline_inverse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ coef ] = get_ppcoef_1seg_full(spl2d, cj_y, jseg)

% function [ coef ] = get_ppcoef_1seg_full(spl2d, cj_y, jseg)
%
% determine PP-spline coefficients for segment jseg

   % average step sizes over 1/2/3 adjacent intervals

   nknot  = length(spl2d.tuj);
   dtj_k1 = (spl2d.tuj(2:end) - spl2d.tuj(1:nknot-1))';        
   dtj_k2 = (spl2d.tuj(3:end) - spl2d.tuj(1:nknot-2))' / 2;
   dtj_k3 = (spl2d.tuj(4:end) - spl2d.tuj(1:nknot-3))' / 3;

   % inverse average step sizes, zero at zero step size

   dtj_inv1 = zeros(size(dtj_k1)); ix  = find(dtj_k1>0); dtj_inv1(ix) = 1./ dtj_k1(ix);
   dtj_inv2 = zeros(size(dtj_k2)); ix  = find(dtj_k2>0); dtj_inv2(ix) = 1./ dtj_k2(ix);
   dtj_inv3 = zeros(size(dtj_k3)); ix  = find(dtj_k3>0); dtj_inv3(ix) = 1./ dtj_k3(ix);

   % backward difference matrices using average step lengths 1/2/3

   D1   = spdiags(dtj_inv1, 0, nknot-1, nknot-2) - spdiags(dtj_inv1(2:end), -1, nknot-1, nknot-2);
   D2   = spdiags(dtj_inv2, 0, nknot-2, nknot-3) - spdiags(dtj_inv2(2:end), -1, nknot-2, nknot-3);
   D3   = spdiags(dtj_inv3, 0, nknot-3, nknot-4) - spdiags(dtj_inv3(2:end), -1, nknot-3, nknot-4);

   % evaluate function and derivative values at position uj(jseg)

   [B1, B2, B3, B4] = eval_bspline_basisfnc( spl2d.tuj, spl2d.tuj(jseg)+1e-9 );

   ay0 = B4 * cj_y;
   ay1 = B3 * D3 * cj_y;
   ay2 = B2 * D2 * D3 * cj_y / 2;
   ay3 = B1 * D1 * D2 * D3 * cj_y / 6;
   coef = [ay0, ay1, ay2, ay3];

   if (0==1)
      tmp = [B4,0,0,0;  B3,0,0;  B2,0;  B1];
      disp(full(tmp(:,jseg:jseg+3)))

      figure(5); clf;
      spy(tmp);
      c = get(gca,'children');
      set(c,'marker','*');
      axis([jseg-1 jseg+4 0 5]);
      axis normal
      hold on
      plot((jseg+3)*[1 1], [0 5], '--');
      grid on;
      title(sprintf('non-zeros for j = %d',jseg))
   end

end % function get_ppcoef_1seg_full

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ coef ] = get_ppcoef_1seg_sparse(spl2d, cj_y, jseg)

% function [ coef ] = get_ppcoef_1seg_sparse(spl2d, cj_y, jseg)
%
% determine PP-spline coefficients for segment jseg

   % average step sizes over 1/2/3 adjacent intervals

   nknotu = length(spl2d.tuj);
   tuj    = spl2d.tuj(jseg-3:jseg+4);            % dtj_k3 uses (jseg+1)+3 - (jseg+1)
   dtj_k1 = (tuj(2:end) - tuj(1:end-1))';        % used for jseg  :jseg+1, last one * 0.0
   dtj_k2 = (tuj(3:end) - tuj(1:end-2))' / 2;    % used for jseg-1:jseg+1, last one * 0.0
   dtj_k3 = (tuj(4:end) - tuj(1:end-3))' / 3;    % used for jseg-2:jseg+1, last one * 0.0

   % inverse average step sizes, zero at zero step size

   dtj_inv1 = zeros(9,1); ix  = find(dtj_k1>0); dtj_inv1(ix) = 1./ dtj_k1(ix);
   dtj_inv2 = zeros(8,1); ix  = find(dtj_k2>0); dtj_inv2(ix) = 1./ dtj_k2(ix);
   dtj_inv3 = zeros(7,1); ix  = find(dtj_k3>0); dtj_inv3(ix) = 1./ dtj_k3(ix);
   % disp(dtj_inv3')

   % evaluate function and derivative values at position tuj(jseg)

   B1   = zeros(1,5);
   B2   = zeros(1,5);
   B3   = zeros(1,5);
   B4   = zeros(1,5);

   u = tuj(4)+1e-9;
   for jj = 4 : 4       % B(j;1) is nonzero for jseg only
      B1(1,jj)    = (u >= tuj(jj)  & u < tuj(jj+1));
   end
   for jj = 3 : 4       % B(j;2) uses B(j+1;1), nonzero for jseg-1 and jseg,
                        % using dt1(jseg:jseg), using 0 * dt1(jseg-1) and 0 * dt1(jseg+1)
      B2(1,jj) = ( (u-tuj(jj  )) .* B1(:,jj  ) * dtj_inv1(jj  )  +  ...
                   (tuj(jj+2)-u) .* B1(:,jj+1) * dtj_inv1(jj+1) ) / 1;
   end
   for jj = 2 : 4       % B(j;3) uses B(j+1;2), nonzero for jseg-2 to jseg, using dt2(jseg-1:jseg)
      B3(1,jj) = ( (u-tuj(jj  )) .* B2(:,jj  ) * dtj_inv2(jj  )  +  ...
                   (tuj(jj+3)-u) .* B2(:,jj+1) * dtj_inv2(jj+1) ) / 2;
   end
   for jj = 1 : 4       % B(j;4) uses B(j+1;3), nonzero for jseg-3 to jseg, using dt3(jseg-2:jseg)
      B4(1,jj) = ( (u-tuj(jj  )) .* B3(:,jj  ) * dtj_inv3(jj  )  +  ...
                   (tuj(jj+4)-u) .* B3(:,jj+1) * dtj_inv3(jj+1) ) / 3;
   end
   % disp(B1); disp(B2); disp(B3); disp(B4)

   d1_cy = zeros(4,1);
   d2_cy = zeros(4,1);
   d3_cy = zeros(4,1);
   d1_cy(2:4) = diff( cj_y(jseg-3:jseg)) .* dtj_inv3(2:4);  % using dt3(jseg-2:jseg)
   d2_cy(3:4) = diff( d1_cy(    2:   4)) .* dtj_inv2(3:4);  % using dt2(jseg-1:jseg)
   d3_cy(4:4) = diff( d2_cy(    3:   4)) .* dtj_inv1(4:4);  % using dt1(jseg  :jseg)

   ay0 = B4(1,1:4) *  cj_y(jseg-3:jseg);
   ay1 = B3(1,2:4) * d1_cy(     2:   4);
   ay2 = B2(1,3:4) * d2_cy(     3:   4) / 2;
   ay3 = B1(1,4:4) * d3_cy(     4:   4) / 6;
   coef = [ay0, ay1, ay2, ay3];

end % function get_ppcoef_1seg_sparse

