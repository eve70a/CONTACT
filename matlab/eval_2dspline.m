
function [ v_out, x_out, y_out, z_out ] = eval_2dspline( spl2d, u_in, v_in, y_in, idebug )

% function [ v_out, x_out, y_out, z_out ] = eval_2dspline( spl2d, u_in, v_in, y_in, idebug )
%
% evaluate parametric 2d spline {tui, tvj, cij_x, cij_y, cij_z} 
%  - if v_in is given, determine [xyz]_out(u_in, v_in) for given (u_in, v_in)
%  - if y_in is given, determine v_out and corresponding x_out,z_out(u_in, v_out)
%
% the inputs (u,v) or (u,y) are considered a list when of the same size, else a tensor grid

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   if (isstruct(spl2d) & isfield(spl2d, 'slc_file'))
      slcs = spl2d;
      spl2d = slcs.spl2d;
   end
   if (nargin<3)
      v_in = [];
   end
   if (nargin<4)
      y_in = [];
   end
   if (nargin<5 | isempty(idebug))
      idebug = 0;
   end
   if ( (isempty(v_in) & isempty(y_in)) | (~isempty(v_in) & ~isempty(y_in)) )
      disp('ERROR: either v_in or y_in must be given');
      return;
   end
   if (min(size(u_in))>1 | min(size(v_in))>1 | min(size(y_in))>1)
      disp('ERROR: u_in, v_in, y_in must be 1D arrays');
      return;
   end
   if (size(u_in,2)>size(u_in,1))
      u_in = u_in';     % make column vector
   end
   if (size(v_in,1)>size(v_in,2))
      v_in = v_in';     % make row vector
   end

   if (spl2d.use_cylindr)
      u_in = apply_wrap_around( spl2d, u_in );
   end

   if (~isempty(v_in))

      % forward calculation u_in, v_in --> x_out, y_out, z_out

      % disp('u_in specified, forward (u,v) --> (x,y,z)')
      v_out = v_in;
      [ x_out, y_out, z_out ] = eval_2dspline_forward(spl2d, u_in, v_in, idebug);

   else

      % inverse calculation u_in, y_in --> v_out

      % disp('y_in specified, inverse (u,y) --> (v)')
      v_out = eval_2dspline_inverse( spl2d, u_in, y_in, idebug );

      % for tensor grid nu x nv, reshape 2d v_out to 1d list
      nu = length(u_in); ny = length(y_in);
      if (nu~=ny)
         u_in  = repmat(u_in, ny, 1);
         v_out = reshape(v_out, nu*ny, 1);
      end

      % forward calculation v_out --> x_out, z_out

      % disp('y_in specified, forward (u,v) --> (x,z)')
      y_out = y_in;
      [ x_out, ~, z_out ] = eval_2dspline_forward(spl2d, u_in, v_out, idebug);
      if (nu~=ny)
         x_out = reshape(x_out, nu, ny);
         y_out = repmat(y_in, nu, 1);
         z_out = reshape(z_out, nu, ny);
      end

   end

end % function eval_2dspline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_out, y_out, z_out] = eval_2dspline_forward(spl2d, u_in, v_in, idebug)

% function [x_out, y_out, z_out] = eval_2dspline_forward(spl2d, u_in, v_in, idebug)
%
% if u_in and v_in have the same size, evaluate 2d spline at pairs (u_in(i), v_in(i))
%                                else, evaluate 2d spline at tensor grid u_in x v_in.

   has_xij = isfield(spl2d, 'cij_x');

   if (all(size(u_in)==size(v_in)))
      use_list = 1; typ_eval = 'list';
   else
      use_list = 0; typ_eval = 'tens.grid';
   end
   k = 4;

   if (idebug>=1)
      disp(sprintf('Forward spline evaluation (u,v) --> (x,y,z), u: %d x %d, v: %d x %d (%s)', ...
                                                                        size(u_in), size(v_in), typ_eval));
   end

   % determine collocation matrix for v-direction: evaluate each B-spline at each output location
   %  - check basic interval [tj(4),tj(nknot-k+1)]

   n0 = nnz(v_in < spl2d.tvj(4));
   n1 = nnz(v_in > spl2d.tvj(end-3));
   if (n0+n1>0)
      disp(sprintf('Warning: there are %d v-positions before and %d after spline v-range [%3.1f,%3.1f]', ...
                n0, n1, spl2d.tvj(4), spl2d.tvj(end-3)));
   end

   numnan = nnz(isnan(v_in));
   if (idebug>=3 | (idebug>=1 & numnan>0))
      disp(sprintf('array v_in has %d NaN-values', nnz(isnan(v_in))));
   end

   [~, ~, ~, Bmat] = eval_bspline_basisfnc( spl2d.tvj, v_in, k );
   % nknotu = length(spl2d.tui); nsplu  = nknotu - k;
   % nknotv = length(spl2d.tvj); nsplv  = nknotv - k;

   % matrix-multply to get 1d spline coefficients for u-direction

   % ci_y: [ nsplu, nout ], Bmat: [ nout, nsplv ], cij_y: [ nsplu, nsplv ], v_in: [ nout ]
   if (has_xij)
      ci_x = (Bmat * spl2d.cij_x')';
   end
   ci_y = (Bmat * spl2d.cij_y')';
   ci_z = (Bmat * spl2d.cij_z')';
   mask = (Bmat * spl2d.mask_j')';

   % determine collocation matrix for u-direction: evaluate each B-spline at each output location

   tiny  = 1e-10;
   n0 = nnz(u_in < spl2d.tui(4)-tiny);
   n1 = nnz(u_in > spl2d.tui(end-3)+tiny);
   if (n0+n1>0)
      disp(sprintf('Warning: there are %d u-positions before and %d after spline u-range [%3.1f,%3.1f]', ...
                n0, n1, spl2d.tui(4), spl2d.tui(end-3)));
   end

   [~, ~, ~, Bmat] = eval_bspline_basisfnc( spl2d.tui, u_in, k );

   % matrix-multply to get output values

   % y_out: [ noutu, noutv ], Bmat: [ noutu, nsplu ], ci_y: [ nsplu, noutv ], u_in: [ noutu ]
   if (has_xij)
      x_out = Bmat * ci_x;
   else
      x_out = u_in * ones(1,length(v_in));       % half-parametric spline: x==u
   end
   y_out = Bmat * ci_y;
   z_out = Bmat * ci_z;
   mask  = Bmat * mask;

   % remove points that don't have full support

   ix    = find(mask<1-tiny);
   x_out(ix) = NaN;
   y_out(ix) = NaN;
   z_out(ix) = NaN;

   % output a list when both inputs are of same size

   if (use_list)
      x_out = diag(x_out);
      y_out = diag(y_out);
      z_out = diag(z_out);
      if (any(size(y_out)~=size(u_in)))
         x_out = x_out';
         y_out = y_out';
         z_out = z_out';
      end
   end

end % function eval_2dspline_forward

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ v_out ] = eval_2dspline_inverse(spl2d, u_in, y_in, idebug)

% function [ v_out ] = eval_2dspline_inverse(spl2d, u_in, y_in, idebug)
%
% if u_in and y_in have the same size, invert 2d spline at pairs (u_in(i), y_in(i))
%                                else, invert 2d spline at tensor grid u_in x y_in.

   tiny   = 1e-10 * max(abs(spl2d.ui));
   k      = 4;
   nknotu = length(spl2d.tui);
   nknotv = length(spl2d.tvj);
   nsplu  = nknotu - k;
   nsplv  = nknotv - k;

   if (all(size(u_in)==size(y_in)))
      use_list = 1; typ_eval = 'list';
      nu    = length(u_in);
      ny    = 1;
      v_out = ones(size(u_in));
   else
      use_list = 0; typ_eval = 'tens.grid';
      nu    = length(u_in);
      ny    = length(y_in);
      v_out = ones(nu, ny);
   end

   if (idebug>=1)
      disp(sprintf('Inverse spline evaluation (u,y) --> (v),   u: %d x %d, y: %d x %d (%s)', ...
                                                                    size(u_in), size(y_in), typ_eval));
   end

   % check basic interval in u-direction, knots [t(4), t(nknot-k+1)]

   if (any(u_in<spl2d.tui(k)) | any(u_in>spl2d.tui(nknotu-k+1)))
      disp(sprintf('ERROR: all input u_in must lie in range of spline u=[%3.1e,%3.1e]', ...
                spl2d.tui(k), spl2d.tui(nknotu-k+1)));
      return;
   end

   % determine collocation matrix for computing coefficients cj_y

   [~, ~, ~, Bmatx] = eval_bspline_basisfnc( spl2d.tui, u_in );

   % loop over output positions u_in

   for iout = 1 : nu

      % determine spline coefficients cj_y at position u_in(iout)

      % cj_y: [ nsplv, 1 ], Bmat: [ nu, nsplu ], cij_y: [ nsplu, nsplv ], u_in: [ nu ]

      cj_y = (Bmatx(iout,:) * spl2d.cij_y)';

      % n+k knots gives n basisfunctions and n spline coefficients; 
      % basisfunction j becomes nonzero at knot j until knot j+k.
      % values of y in segment j = [t_j,t_j+1) are defined by spline coefficients j-k+1:j
      % n+k knots gives n+k-1 segments. The first k-1 and last k-1 segments are outside the basic interval
      % this gives n-k+1 segments in the basic interval with numbers k to n

      % determine interval [ylow, yhig] for segments j, j = k to n

      tmp            = NaN * ones(4, nknotv);
      tmp(:,k:nsplv) = [ cj_y(1:nsplv-3)';  cj_y(2:nsplv-2)';  cj_y(3:nsplv-1)';  cj_y(4:nsplv)' ];

      ylow = min(tmp); yhig = max(tmp);

      % no inner loop when using list-input, size(u_in) == size(y_in)

      if (use_list)
         j0 = iout; j1 = iout;
      else
         j0 = 1; j1 = ny;
      end

      for jout = j0 : j1

         % determine segments jseg that may contain y(jout) at xout

         y    = y_in(jout);
         jseg = find( ylow(1:nsplv)<=y & y<=yhig(1:nsplv) );
         nseg = length(jseg);

         if (nseg<=0)
            % disp(sprintf('ERROR: no segments for (u,y) = (%6.1f,%6.1f)', u_in(iout), y));
            found = 0;
         else
            if (idebug>=2)
               disp(sprintf(['(i,j)=(%2d,%2d), (u,y)=(%6.1f,%6.1f): found %d possible segments jseg = ', ...
                   '%d, %d, %d, %d, %d, %d'], iout, jout, u_in(iout), y_in(jout), nseg, jseg(1:min(6,nseg))));
               if (idebug>=-4)

                  [~, ~, ~, Bmatv] = eval_bspline_basisfnc( spl2d.tvj, spl2d.vj );
                  % ytmp: [ nv, 1 ], Bmatv: [ nv, nv ], cj_y: [ nv, 1 ]
                  ytmp = Bmatv * cj_y;
                  tvj  = spl2d.tvj;
                  tgrev(1:nsplv) = (tvj(2:nsplv+1) + tvj(3:nsplv+2) + tvj(4:nsplv+3)) / 3;
                  tm   = (tgrev(1:nsplv-1) + tgrev(2:end)) / 2;

                  figure(7); clf; hold on;
                  tmpt  = reshape( [1;1]*tvj(1:nsplv+1), 1, 2*nsplv+2 );
                  tmpyl = reshape( [1;1]*ylow(1:nsplv), 1, 2*nsplv );
                  tmpyh = reshape( [1;1]*yhig(1:nsplv), 1, 2*nsplv );
                  plot(tmpt(2:end-1), tmpyl);
                  plot(tmpt(2:end-1), tmpyh);
                  plot(spl2d.vj([1,end]), y*[1 1], '--');
                  plot(spl2d.vj, ytmp);
                  plot([1;1]*(tvj(jseg)+tvj(jseg+1))/2, [ylow(jseg);yhig(jseg)], 'color',matlab_color(3))
                  plot(spl2d.tvj(nsplv+1), cj_y(end), '*');
                  for jj = jseg
                      text( (tvj(jj)+tvj(jj+1))/2, ylow(jj), sprintf('j=%d',jj), ...
                        'horizontalalignment','center', 'verticalalignment','top');
                  end
                  grid on
                  xlabel('v_j'); ylabel('y_r');
                  legend('y_{low}(v;u)', 'y_{hig}(v;u)', 'target y_{out}', 'actual y(v), u=u_{out}', ...
                                                                                'location','northwest')
               end
            end

            % try potential segments one-by-one until a solution is found

            jj = 1;
            found = 0;
            while(jj<=nseg & ~found)

               % spline segment: interval [v_j, v_{j+1}]

               vseg = spl2d.tvj( jseg(jj)+[0,1] );

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

               [v, found] = solve_cubic_eq(coef, vseg, y, idebug, iout, jout, jseg(jj));
               if (found & idebug>=2)
                  disp(sprintf('jj = %d: found v =%7.2f', jj, v));
               elseif (~found & idebug>=3)
                  disp(sprintf('jj = %d: no solution',jj));
               end

               found = 1;
               if (~found), jj = jj + 1; end
            end % jj

         end % nseg>0

         % store solution

         if (~found), 
            if (idebug>=2)
               disp(sprintf('No solution for (i,j)=(%d,%d), setting NaN', iout, jout));
            end
            v = NaN;
         end
         if (use_list)
            v_out(iout) = v;
         else
            v_out(iout,jout) = v;
         end

      end % for jout
   end % for iout

end % function eval_2dspline_inverse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ th_ev ] = apply_wrap_around( spl2d, th_ev )

% - apply wrap-around th_ev --> [-pi,pi)
% - restrict th_ev \in spline range [th0,th1] to obtain constant extrapolation

   th_ev = mod(th_ev+pi, 2*pi) - pi;

   th0   = spl2d.tui(4);
   th1   = spl2d.tui(end-3);
   th_ev = max(th0, min(th1, th_ev));

end % function apply_wrap_around

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ coef ] = get_ppcoef_1seg_full(spl2d, cj_y, jseg)

% function [ coef ] = get_ppcoef_1seg_full(spl2d, cj_y, jseg)
%
% determine PP-spline coefficients for segment jseg

   % average step sizes over 1/2/3 adjacent intervals

   nknot  = length(spl2d.tvj);
   dtj_k1 = (spl2d.tvj(2:end) - spl2d.tvj(1:nknot-1))';        
   dtj_k2 = (spl2d.tvj(3:end) - spl2d.tvj(1:nknot-2))' / 2;
   dtj_k3 = (spl2d.tvj(4:end) - spl2d.tvj(1:nknot-3))' / 3;

   % inverse average step sizes, zero at zero step size

   dtj_inv1 = zeros(size(dtj_k1)); ix  = find(dtj_k1>0); dtj_inv1(ix) = 1./ dtj_k1(ix);
   dtj_inv2 = zeros(size(dtj_k2)); ix  = find(dtj_k2>0); dtj_inv2(ix) = 1./ dtj_k2(ix);
   dtj_inv3 = zeros(size(dtj_k3)); ix  = find(dtj_k3>0); dtj_inv3(ix) = 1./ dtj_k3(ix);

   % backward difference matrices using average step lengths 1/2/3

   D1   = spdiags(dtj_inv1, 0, nknot-1, nknot-2) - spdiags(dtj_inv1(2:end), -1, nknot-1, nknot-2);
   D2   = spdiags(dtj_inv2, 0, nknot-2, nknot-3) - spdiags(dtj_inv2(2:end), -1, nknot-2, nknot-3);
   D3   = spdiags(dtj_inv3, 0, nknot-3, nknot-4) - spdiags(dtj_inv3(2:end), -1, nknot-3, nknot-4);

   % evaluate function and derivative values at position vj(jseg)

   [B1, B2, B3, B4] = eval_bspline_basisfnc( spl2d.tvj, spl2d.tvj(jseg)+1e-9 );

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

   nknotv = length(spl2d.tvj);
   tvj    = spl2d.tvj(jseg-3:jseg+4);            % dtj_k3 uses (jseg+1)+3 - (jseg+1)
   dtj_k1 = (tvj(2:end) - tvj(1:end-1))';        % used for jseg  :jseg+1, last one * 0.0
   dtj_k2 = (tvj(3:end) - tvj(1:end-2))' / 2;    % used for jseg-1:jseg+1, last one * 0.0
   dtj_k3 = (tvj(4:end) - tvj(1:end-3))' / 3;    % used for jseg-2:jseg+1, last one * 0.0

   % inverse average step sizes, zero at zero step size

   dtj_inv1 = zeros(9,1); ix  = find(dtj_k1>0); dtj_inv1(ix) = 1./ dtj_k1(ix);
   dtj_inv2 = zeros(8,1); ix  = find(dtj_k2>0); dtj_inv2(ix) = 1./ dtj_k2(ix);
   dtj_inv3 = zeros(7,1); ix  = find(dtj_k3>0); dtj_inv3(ix) = 1./ dtj_k3(ix);
   % disp(dtj_inv3')

   % evaluate function and derivative values at position tvj(jseg)

   B1   = zeros(1,5);
   B2   = zeros(1,5);
   B3   = zeros(1,5);
   B4   = zeros(1,5);

   u = tvj(4)+1e-9;
   for jj = 4 : 4       % B(j;1) is nonzero for jseg only
      B1(1,jj)    = (u >= tvj(jj)  & u < tvj(jj+1));
   end
   for jj = 3 : 4       % B(j;2) uses B(j+1;1), nonzero for jseg-1 and jseg,
                        % using dt1(jseg:jseg), using 0 * dt1(jseg-1) and 0 * dt1(jseg+1)
      B2(1,jj) = ( (u-tvj(jj  )) .* B1(:,jj  ) * dtj_inv1(jj  )  +  ...
                   (tvj(jj+2)-u) .* B1(:,jj+1) * dtj_inv1(jj+1) ) / 1;
   end
   for jj = 2 : 4       % B(j;3) uses B(j+1;2), nonzero for jseg-2 to jseg, using dt2(jseg-1:jseg)
      B3(1,jj) = ( (u-tvj(jj  )) .* B2(:,jj  ) * dtj_inv2(jj  )  +  ...
                   (tvj(jj+3)-u) .* B2(:,jj+1) * dtj_inv2(jj+1) ) / 2;
   end
   for jj = 1 : 4       % B(j;4) uses B(j+1;3), nonzero for jseg-3 to jseg, using dt3(jseg-2:jseg)
      B4(1,jj) = ( (u-tvj(jj  )) .* B3(:,jj  ) * dtj_inv3(jj  )  +  ...
                   (tvj(jj+4)-u) .* B3(:,jj+1) * dtj_inv3(jj+1) ) / 3;
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

