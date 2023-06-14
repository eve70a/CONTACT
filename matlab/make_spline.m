
function [ spl ] = make_spline( s, y, z, lambda, wgt, ikinks, iaccel, use_bspline, ds_bspl, ...
                                                                      use_deriv, use_repl, idebug )

% function [ spl ] = make_spline( [s], y, [z], lambda, [wgt], [ikinks], [iaccel], [use_bspline], ...
%                                                           [ds_bspl], [use_deriv], [use_repl], [idebug] )
%
% compute parametric smoothing spline {s, a0y-a3y, a0z-a3z} (pp-form) for data {y(s),z(s)},
%   s         = arc-length, will be computed if not provided
%   lambda    = weight for smoothness, default 0 (no smoothing)
%   wgt       = weight for data-points, default all ones.
%   ikinks    = discontinuous 1st derivatives: sharp corners/kinks at list of indices [ i1, i2, ... ]
%   iaccel    = discontinuous 2nd derivatives: radius jumps, acceleration points [i1, i2, ...]
%   ds_bspl   = target step size for B-spline method
%   use_deriv = 2 or 3 for penalty on 2nd or 3rd derivative
%   use_repl  = 1=replicate knots at s(1) and s(end) or 0=extend knots, replicating intervals ds(1), ds(end)

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

if (nargin>=3 & isempty(y) & ~isempty(z))
   y      = ones(size(z));
end
if (nargin<3 | isempty(z))
   z      = ones(size(y));
end
if (nargin<4 | isempty(lambda))
   lambda = 0;  % default: no smoothing
end
if (nargin<5 | isempty(wgt))
   wgt = ones(size(y));
end
if (nargin<6)
   ikinks = []; % default: no kinks/sharp corners
end
if (nargin<7)
   iaccel = []; % default: no acceleration points
end
if (nargin<8 | isempty(use_bspline))
   use_bspline = 1;  % default: new B-spline method
end
if (nargin<9 | isempty(ds_bspl))
   ds_bspl = 2;
end
if (nargin<10 | isempty(use_deriv))
   use_deriv = 3;
   if (~use_bspline), use_deriv = 2; end
end
if (nargin<11 | isempty(use_repl))
   use_repl = 1;
end
if (nargin<12 | isempty(idebug))
   idebug = 0;
end

if (~use_bspline & ~isempty(iaccel))
   disp('Accelerations are ignored in the PP-form method.');
end
if (~use_bspline & use_deriv~=2)
   disp('The PP-spline uses use_deriv=2.');
   use_deriv = 2;
end

npnt = length(y);
if (size(s,1)==1), s = s'; end % make column vector
if (size(y,1)==1), y = y'; end
if (size(z,1)==1), z = z'; end
if (size(wgt,1)==1), wgt = wgt'; end
if (size(ikinks,1) > size(ikinks,2)) % make row vector
   ikinks = ikinks';
end
if (size(iaccel,1) > size(iaccel,2)) % make row vector
   iaccel = iaccel';
end

ix = find(iaccel<=0 | iaccel>npnt);
if (~isempty(ix))
   disp(sprintf('ERROR: %d accelerations out of range [1,%d]', length(ix), npnt));
   disp(iaccel')
   return;
end

if (isempty(s))
   s = make_arclength( y, z );
end
if (any(size(s)~=size(y)) | any(size(z)~=size(y)) | any(size(wgt)~=size(y)))
   disp('Arrays s, y, z, wgt have incompatible sizes');
   disp( [size(s); size(y); size(z); size(wgt)] );
   return;
end

if (use_bspline)
   spl = make_bspline( s, y, z, lambda, wgt, ikinks, iaccel, ds_bspl, use_deriv, use_repl, idebug );
else
   spl = make_pp_spline( s, y, z, lambda, wgt, ikinks, idebug );
end

end % function make_spline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ spl ] = make_bspline( s_in, y_in, z_in, lambda, wgt, ikinks, iaccel, ds_bspl, ...
                                                                use_deriv, use_repl, idebug )

% function [ spl ] = make_bspline( s_in, y_in, z_in, lambda, wgt, ikinks, iaccel, ds_bspl, ...
%                                                               use_deriv, use_repl, [idebug] )
%
% compute parametric cubic smoothing spline (y(s),z(s)) using least squares B-spline method with
% penalty on 2nd or 3rd derivative

if (nargin<11 | isempty(idebug))
   idebug = 1;
end
k  = 4;  % spline order 4, cubic splines

% check size of ds_bspl: too small value spoils the condition of matrix M

ds_min = 0.1 * lambda^(1/(2*use_deriv));
if (ds_bspl < ds_min)
   if (idebug>=1)
      disp(sprintf('Input ds_bspl = %6.1e too small, using ds_min = %6.3f', ds_bspl, ds_min));
   end
   ds_bspl = ds_min;
end

% determine appropriate knot-vector tj
%  - simple algorithm: just apply ds_bspl as much as possible, needs some smoothing
%  - advanced algorithm: no more than 1 knot per meas.segment, no smaller than ds_bspl, fewer at kinks, accel and boundaries

use_simple  = 1;      % default
lmb_min     = 1e-6;
if (lambda<lmb_min), use_simple = 0; end

[tj, keep_reduc] = make_knot_vector( s_in, ikinks, iaccel, ds_bspl, use_simple, use_repl, idebug );
nknot = length(tj);

% determine collocation matrix: evaluate each B-spline at each measurement location

nmeas = length(s_in);
[~, ~, ~, Bmat] = eval_bspline_basisfnc( tj, s_in );
if (idebug>=5 & nmeas<20), disp(full(Bmat)); end

% weight matrix W for input data

Wmat = spdiags(wgt, 0, nmeas, nmeas);
% disp(Wmat)

% average step sizes over 1/2/3 adjacent intervals

dtj_k1 =  tj(2:end) - tj(1:nknot-1);
dtj_k2 = (tj(3:end) - tj(1:nknot-2)) / 2;
dtj_k3 = (tj(4:end) - tj(1:nknot-3)) / 3;

% inverse average step sizes, zero at zero step size

dtj_inv1 = zeros(size(dtj_k1)); ix  = find(dtj_k1>0); dtj_inv1(ix) = 1./ dtj_k1(ix);
dtj_inv2 = zeros(size(dtj_k2)); ix  = find(dtj_k2>0); dtj_inv2(ix) = 1./ dtj_k2(ix);
dtj_inv3 = zeros(size(dtj_k3)); ix  = find(dtj_k3>0); dtj_inv3(ix) = 1./ dtj_k3(ix);

% backward difference matrices using average step lengths 1/2/3

D1   = spdiags(dtj_inv1, 0, nknot-1, nknot-2) - spdiags(dtj_inv1(2:end), -1, nknot-1, nknot-2);
D2   = spdiags(dtj_inv2, 0, nknot-2, nknot-3) - spdiags(dtj_inv2(2:end), -1, nknot-2, nknot-3);
D3   = spdiags(dtj_inv3, 0, nknot-3, nknot-4) - spdiags(dtj_inv3(2:end), -1, nknot-3, nknot-4);

% penalty matrix D^T * C * D for penalty on 2nd or 3rd derivative

if (use_deriv==2)
   Dmat = D2(3:nknot-k,2:nknot-k) * D3(2:nknot-k,1:nknot-k);
   Cmat =   2/3 * spdiags( dtj_k2(1:nknot-k),  0, nknot-k, nknot-k ) + ...
          + 1/6 * spdiags( dtj_k1(2:nknot-k), -1, nknot-k, nknot-k ) + ...
          + 1/6 * spdiags( dtj_k1(1:nknot-k),  1, nknot-k, nknot-k );
   Cmat = Cmat(3:nknot-k,3:nknot-k);
else
   Dmat = D1(4:nknot-k,3:nknot-k) * D2(3:nknot-k,2:nknot-k) * D3(2:nknot-k,1:nknot-k);
   Cmat = spdiags( dtj_k1, 0, nknot-1, nknot-1 );
   Cmat = Cmat(4:nknot-k,4:nknot-k);
end

% Pmat = Dmat' * Cmat * Dmat;
% dt = dtj_k1(1);
% disp(full(dt^5*Pmat(1:5,1:6)))
% disp(flipud(fliplr( full(dt^5*Pmat(end-4:end,end-5:end))) ))

if (idebug>=3)
   for i = 1:min(10,nknot-1)
      tmp = full( Dmat(i,max(1,i-3):min(i,nknot-4) ));
      disp(sprintf('i=%3d: d(i,:) = %10.3e %10.3e %10.3e %10.3e', i, tmp))
   end
end

% combine basisfunctions at segments that use a reduced order

if (any(~keep_reduc))
   if (idebug>=0)
      disp(sprintf('Total %d knots, %d knots removed, %d after reduction', nknot, nnz(~keep_reduc), ...
                        nnz(keep_reduc)))
   end
 
   Bmat = reduce_matrix(Bmat, tj, keep_reduc, idebug);
   Dmat = reduce_matrix(Dmat, tj, keep_reduc);
end

% determine spline coefficients, solving the least squares optimization equation

if (idebug>=2)
   if (0==1)
      disp('Bt*W*B:');
      tmp = Bmat'*Wmat*Bmat;
   else
      disp('B:');
      tmp = Bmat;
   end
   if (size(tmp,1)>20)
      disp(full(tmp(1:10, 1:10)))
      disp(full(tmp(end-9:end, end-9:end)))
   else
      disp(full(tmp))
   end
   % figure(6); spy(tmp);

   if (lambda>1e-12)
      for i = 1: min(10,nknot-1)
         tmp = full( Dmat(i,max(1,i-5):min(i-2,nknot-4) ));
         disp(sprintf('i=%3d: dr(i,:) = %10.3e %10.3e %10.3e %10.3e', i, tmp));
      end

      tmp = Dmat'*Cmat*Dmat;
      if (size(tmp,1)>20)
         disp('Dt*C*D (top-left):');
         disp(sprintf('%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n', ...
                                                                        full(tmp(1:10, 1:10))))
         disp('Dt*C*D (bottom-right):');
         disp(sprintf('%10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n', ...
                                                                        full(tmp(end-9:end, end-9:end))))
      else
         disp('Dt*C*D:');
         disp(full(tmp))
      end

      % disp(log10(diag(tmp)))

      disp('[Bt*W*B + lmb * Dt*C*D, rhs]:')
      tmp = [Bmat'*Wmat*Bmat + lambda * Dmat'*Cmat*Dmat, Bmat'*Wmat*[y_in,z_in]];
      if (size(tmp,1)>10)
         disp(full(tmp(1:10, [1:8,end-1:end])))
      else
         disp(full(tmp));
      end
   end

end

Mat = (Bmat'*Wmat*Bmat + lambda * Dmat'*Cmat*Dmat);
cnd = condest(Mat);
if (cnd>1e10)
   disp(sprintf('make_bspline: matrix condition %6.2e, spline may be unstable', cnd));
end

coef_y = Mat \ (Bmat' * Wmat * y_in);
coef_z = Mat \ (Bmat' * Wmat * z_in);

if (idebug>=3)
   disp('B-spline coefficients:')
   if (length(coef_y)>20)
      disp([ [1:10,length(coef_y)+[-9:0]]', coef_y([1:10,end-9:end]), coef_z([1:10,end-9:end])])
   else
      disp([ [1:length(coef_y)]', coef_y, coef_z])
   end
end

coef_y = expand_vector( coef_y, tj, keep_reduc );
coef_z = expand_vector( coef_z, tj, keep_reduc );

if (idebug>=3 & ~all(keep_reduc))
   disp('Expanded coefficients:')
   disp([ [1:10,length(coef_y)+[-9:0]]', coef_y([1:10,end-9:end]), coef_z([1:10,end-9:end])])
end

ksi = unique(tj);
spl = struct('s',ksi, 'tj',tj, 'cy',coef_y, 'cz',coef_z);
spl = make_ppform( spl, D1, D2, D3 );

end % function make_bspline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ is_ok ] = check_kink_accel( si, ikinks, iaccel )

% function [ is_ok ] = check_kink_accel( si, ikinks, iaccel )
%
% check requirements of B-spline method on kinks and accelerations

is_ok  = 0;

nmeas  = length(si);
nkink  = length(ikinks);
naccel = length(iaccel);
ikinks = sort(ikinks);
iaccel = sort(iaccel);
if (size(ikinks,1)>size(ikinks,2)), ikinks = ikinks'; end
if (size(iaccel,1)>size(iaccel,2)), iaccel = iaccel'; end

% perform checks on ikinks and iaccel:
%  - all kinks and accelerations must be different

if (length(unique([ikinks iaccel])) < length([ikinks iaccel]))
   disp('ERROR: check_kink_accel: ikinks and iaccel must all be different.')
   disp(ikinks);
   disp(iaccel);
   is_ok = -1;
   return;
end

%  - all kinks must lie in interior of measurement range

if (any(ikinks<=1 | ikinks>=nmeas))
   disp(sprintf('ERROR: check_kink_accel: ikinks must be >= 2 and <= nmeas = %d', nmeas))
   disp(ikinks);
   is_ok = -1;
   return;
end

%  - all accelerations must lie in interior of measurement range, away from boundaries

if (any(iaccel<=3 | iaccel>=nmeas-2))
   disp(sprintf('ERROR: check_kink_accel: iaccel must be >= 4 and <= nmeas-3 = %d', nmeas-3))
   disp(iaccel);
   is_ok = -1;
   return;
end

%  - all kinks must have 0 or >=2 points between them, accelerations and kinks need >= 2 points between them

iall  = [  1, nmeas,  ikinks,          iaccel           ];
itype = [  0,    0 ,  1*ones(1,nkink), 2*ones(1,naccel) ];
[iall, iperm] = sort(iall);
itype = itype(:,iperm);
idist = diff(iall);

if (any(idist<=2 & (itype(1:end-1)==2 | itype(2:end)==2) ))
   disp('ERROR: check_kink_accel: accelerations must have at least 2 points between them and');
   disp('                         kinks, boundaries and other accelerations');
   disp([iall; itype])
   is_ok = -1;
   return
end
if (any(idist==2 & itype(1:end-1)<=1 & itype(2:end)==1))
   disp(sprintf('ERROR: check_kink_accel: kinks must be adjacent or have >=2 points between them'));
   disp(ikinks);
   is_ok = -1;
   return;
end

is_ok = 1;

end % function check_kink_accel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ tj, keep_reduc ] = make_knot_vector( si, ikinks, iaccel, ds_bspl, use_simple, ...
                                                                                use_repl, idebug )

% function [ tj, keep_reduc ] = make_knot_vector( s_meas, ikinks, iaccel, ds_bspl, use_simple, ...
%                                                                               use_repl, [idebug] )
%
% determine knot-vector with repeated knots at kinks and accelerations
% reduce knots in intervals between double kinks or kinks adjacent to boundaries
%  - ds_bspl       = target step size
%  - use_simple    = true:  apply ds_bspl as much as possible
%                    false: at most 1 knot per data point
%  - use_repl      = knots at boundaries: replicated [1,1,1,1]*tj(1) or extended tj(1)+[-3,-2,-1,0]*dt1
%  - keep_reduc    = flag indicating that tj is/isnt kept in reduced knots t^*

if (nargin<6 | isempty(idebug))
   idebug = 1;
end
if (idebug>=3)
   make_plot = 6;
else
   make_plot = 0;
end

nmeas  = length(si);
ikinks = sort(ikinks);
iaccel = sort(iaccel);
if (size(ikinks,1)>size(ikinks,2)), ikinks = ikinks'; end
if (size(iaccel,1)>size(iaccel,2)), iaccel = iaccel'; end

is_ok = check_kink_accel( si, ikinks, iaccel );
if (is_ok<=0)
   disp('ERROR with kinks and/or accelerations')
   return;
end

if (use_simple)
   [ tj ] = make_knot_vector_simple(si, ikinks, iaccel, ds_bspl, use_repl, idebug);
else
   [ tj ] = make_knot_vector_advanced(si, ikinks, iaccel, ds_bspl, use_repl, idebug);
end
nknot = length(tj);

% filter knots for segments with reduced order

iall   = [1, ikinks, nmeas];
iclose = find(diff(iall)==1);
isegm  = iall(iclose);

keep_reduc = ones(size(tj));
for i = 1 : length(isegm)
   sr   = si(isegm(i));
   jsta = find(tj>sr, 1, 'first') - 4;
   keep_reduc(jsta+[1,2]) = 0;
end

% plot knot multiplicities & reduced knots

if (make_plot)
   ksi  = unique(tj);
   nksi = length(ksi);
   mult_j = zeros(nksi,1);
   mult_r = zeros(nksi,1);
   for iksi = 1 : nksi
      mult_j(iksi) = nnz(tj==ksi(iksi));
      mult_r(iksi) = nnz(tj==ksi(iksi) & keep_reduc);
   end

   figure(make_plot); clf; hold on;
   plot(si, zeros(size(si)), '-o');

   for i = 1 : length(si)
      text(si(i), -0.03, num2str(i), 'horizontalalignment','center', 'verticalalignment','top')
   end

   for m = 1 : 4
      ik = find(mult_j>=m);
      plot(ksi(ik), 0.1*m*ones(size(ik)), '*', 'color',matlab_color(2));
      ik = find(mult_r>=m);
      plot(ksi(ik), 0.1*m*ones(size(ik)), 'o', 'color',matlab_color(4));
   end
   v = axis; v(3:4) = [-0.2 1]; axis(v);
   grid on;
   legend('measurement s_i', 'knots t_j', 'reduced knots t^*_j', 'location','northwest');
   xlabel('s / t [mm]');
end

if (any(isnan(tj)))
   disp(sprintf('knot vector has %d NaN-values at positions %d, %d, %d, %d, %d, %d, %d, %d', ...
        nnz(isnan(tj)), find(isnan(tj))));
end
if (idebug>=3)
   disp(' tj,  keep_reduc:')
   disp([ tj , keep_reduc ]')
%  disp(tj(end-5:end)')
end

end % function make_knot_vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ tj ] = make_knot_vector_simple( si, ikinks, iaccel, ds_bspl, use_repl, idebug )

% function [ tj ] = make_knot_vector_simple( s_meas, ikinks, iaccel, ds_bspl, use_repl, [idebug] )
%
% simple version that just applies ds_bspl as much as possible -- needs lambda>0 to avoid singular systems

nmeas  = length(si);
nkink  = length(ikinks);
naccel = length(iaccel);

% define sections between kinks and acceleration points

iranges = sort([ 1 , ikinks, iaccel, nmeas ]);
nsec    = length(iranges) - 1;

tj   = [];
tiny = 1e-12;

% loop over sections

for isec = 1 : nsec

   % take next section

   sec     = [iranges(isec), iranges(isec+1)];

   % determine actual step dt \approx ds_bspl that fits integer multiple times

   len_sec = si(sec(2)) - si(sec(1));
   if (diff(sec)<=1)
      nintv   = 1;      % double kink or kink next to boundary
   else
      nintv   = max(1, round(len_sec / ds_bspl));
   end
   dt      = len_sec / nintv;

   % add knots for this section including start/end-points

   tj   = [tj; si(sec(1)) + [0:nintv]' * dt ];

   % accelerations will be repeated once, end of isec & start of isec+1
   % kinks will be repeated once more

   if (any(ikinks==sec(2)))
      tj = [ tj ; tj(end) ];
      if (idebug>=1)
         disp(sprintf('kink i = %3d: repeat tj = %7.3f', sec(2), tj(end)));
      end
   end
end

if (use_repl)
   % repeat first and last knots three times for cubic spline
   tj   = [ [1;1;1]*tj(1); tj; [1;1;1]*tj(end) ];
else
   % extend by repeating first and last steps three times for cubic spline
   dt0  = tj(2) - tj(1);
   dt1  = tj(end) - tj(end-1);
   tj   = [ tj(1)+[-3:0]'*dt0; tj(2:end-1); tj(end)+[0:3]'*dt1 ];
end

end % make_knot_vector_simple

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ tj ] = make_knot_vector_advanced( si, ikinks, iaccel, ds_bspl, use_repl, idebug )

% function [ tj ] = make_knot_vector_advanced( s_meas, ikinks, iaccel, ds_bspl, use_repl, [idebug] )
%
% advanced knot-vector with steps >= ds_bspl, using at most one knot per data-interval, 
%          with measurements skipped for not-a-knot conditions at boundaries, kinks and accelerations

% define sections between kinks and acceleration points

nmeas   = length(si);
iranges = sort([ 1 , ikinks, iaccel, nmeas ]);
nsec    = length(iranges) - 1;

tj   = [];
tiny = 1e-12;

% loop over sections

for isec = 1 : nsec

   % take next section

   sec     = [iranges(isec), iranges(isec+1)];

   % determine whether 2nd measurement and before-last measurements should be skipped

   start_kink = (any(ikinks==sec(1)));
   end_kink   = (any(ikinks==sec(2)));
   end_accel  = (any(iaccel==sec(2)));

   skip_2nd   = (isec==1 | start_kink);
   skip_2last = (isec==nsec | end_kink | end_accel);

   % determine step dt \approx ds_bspl that fits integer multiple times
   %    this will be the default for steps in this section

   len_sec = si(sec(2)) - si(sec(1));
   nintv   = round(len_sec / ds_bspl);
   dt      = len_sec / nintv;

   if (idebug>=1)
      disp(sprintf('section %d: ip=[%3d,%3d], s=[%8.3f,%8.3f], dt=%7.4f, nintv=%4d', ...
                   isec, sec, si(sec(1)), si(sec(2)), dt, nintv));
      if (isec==1),    disp('           first section - skip 2nd measurement'); end
      if (start_kink), disp('           kink at start - skip 2nd measurement'); end
      if (isec==nsec), disp('           last section  - skip before-last measurement'); end
      if (end_kink),   disp('           kink at end   - skip before-last measurement'); end
      if (end_accel),  disp('           accel at end  - skip before-last measurement'); end
   end

   % add knots for this section

   tj   = [tj; si(sec(1))];

   iprf = sec(1);
   while(iprf < sec(2))
      % disp(sprintf('iprf=%2d, tj(end)=%6.3f, si+dt=%6.3f', iprf, tj(end), si(iprf)+dt));

      % take step dt

      tj_new = tj(end) + dt;

      % enlarge if needed to go to next meas.segment

      tj_new = max(tj_new, si(iprf+1));

      % enlarge if needed to go to third meas.segment

      if (skip_2nd & iprf==sec(1) & sec(2)-sec(1)>=2)
         tj_new = max(tj_new, si(iprf+2));
      end

      % if we are inside the very last two meas.segments, jump to end of section

      if (skip_2last & iprf>=sec(2)-2)
         tj_new = si(sec(2));
      end

      % if step takes us inside the last meas.segment, jump to end of section

      if (skip_2last & tj_new >= si(sec(2)-1))
         tj_new = si(sec(2));
      end

      % if step takes us inside before-last meas.segment, go to beginning of this segment

      if (skip_2last & sec(2)-sec(1)>=2 & tj_new >= si(sec(2)-2) & tj_new <= si(sec(2)-1))
         tj_new = si(sec(2)-2);
      end

      % never go beyond end of section

      tj_new = min(tj_new, si(sec(2)));

      % add to list of knots

      tj = [tj; tj_new];

      % increment iprf to the last measurement that has been covered

      iprf = find(si<=tj(end)+tiny, 1, 'last');
   end

   % accelerations will be repeated once, end of isec / start of isec+1
   % kinks will be repeated once more

   if (end_kink)
      tj = [ tj ; tj(end) ];
      if (idebug>=1)
         disp(sprintf('kink i = %3d: repeat tj = %7.3f', sec(2), tj(end)));
      end
   end
end

if (use_repl)
   % repeat first and last knots three times for cubic spline
   tj   = [ [1;1;1]*tj(1); tj; [1;1;1]*tj(end) ];
else
   % extend by repeating first and last steps three times for cubic spline
   dt0  = tj(2) - tj(1);
   dt1  = tj(end) - tj(end-1);
   tj   = [ tj(1)+[-3:0]'*dt0; tj(2:end-1); tj(end)+[0:3]'*dt1 ];
end

end % make_knot_vector_advanced

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Rmat ] = reduce_matrix(Mat, tj, keep_reduc, idebug)

% function [ Rmat ] = reduce_matrix(Mat, tj, keep_reduc)
%
% Combine columns for segments using a reduced order, e.g. linear segments between double kinks
%   tj         == knot vector (1:nknot)
%   keep_reduc == flags for knots to be kept in reduced knot vector tr = t^* (1:nreduc)

if (nargin<4 | isempty(idebug))
   idebug = 0;
end

nknot  = length(tj);

jsta_reduc = find(keep_reduc(1:end-1) & ~keep_reduc(2:end));
if (idebug>=2 & ~isempty(jsta_reduc))
   disp('Using reduced order at knots jsta =');
   disp(jsta_reduc');
end

% select indices j_orig of knots that are kept in reduced knot vector

j_orig = find(keep_reduc);
tr     = tj(j_orig); 

% new number for knots that are kept

j_reduc = cumsum(keep_reduc) .* keep_reduc;
nreduc  = length(tr);

if (size(Mat,2)~=nknot-4)
   disp('Incorrect size Mat')
   disp([size(Mat,2), nknot-4])
   return
end

% combine columns of Mat(.,nknot-4) into Rmat(.,nreduc-4)

%  - copy columns that are kept as a whole

Rmat  = Mat(:,j_orig(1:nreduc-4));

%  - add columns that are deleted

for j0 = jsta_reduc
   Rmat(:,j_reduc(j0)  ) = Rmat(:,j_reduc(j0)  ) + 2/3 * Mat(:,j0+1) + 1/3 * Mat(:,j0+2);
   Rmat(:,j_reduc(j0+3)) = Rmat(:,j_reduc(j0+3)) + 1/3 * Mat(:,j0+1) + 2/3 * Mat(:,j0+2);
end


end % function reduce_matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ a ] = expand_vector( ra, tj, keep_reduc, idebug )

% function [ a ] = expand_vector( ra, tj, keep_reduc )
%    ra = coefficients for reduced basisfunctions a^*
%    tj = full knot vector t_j (1:nknot)
%    keep_reduc = flags for knots kept in reduced knot vector tr = t^* (1:nreduc)

if (nargin<4 | isempty(idebug))
   idebug = 0;
end

nknot   = length(tj);
nreduc  = nnz(keep_reduc);

% locate knots jsta_reduc at start of segment with reduced order

jsta_reduc = find(keep_reduc(1:end-1) & ~keep_reduc(2:end));
if (idebug>=2 & ~isempty(jsta_reduc))
   disp('Using reduced order at knots jsta =');
   disp(jsta_reduc');
end

% select indices j_orig of knots that are kept in reduced knot vector

j_orig = find(keep_reduc);

% new number for knots that are kept

j_reduc = cumsum(keep_reduc) .* keep_reduc;

% expand coefficients ra(1:nreduc-4) to coefficients a(1:nknot-4)

a = zeros(nknot-4,1);

for j = 1 : nknot-4
   if     (~keep_reduc(j) & any(jsta_reduc==j-1))

      a(j) = 2/3 * ra(j_reduc(j-1)) + 1/3 * ra(j_reduc(j+2));

   elseif (~keep_reduc(j) & any(jsta_reduc==j-2))

      a(j) = 1/3 * ra(j_reduc(j-2)) + 2/3 * ra(j_reduc(j+1));

   elseif (keep_reduc(j))

      a(j) = ra(j_reduc(j));

   end
end

end % function expand_vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ spl ] = make_pp_spline( s, y, z, lambda, wgt, ikinks, idebug )

% function [ spl ] = make_pp_spline( s, y, z, lambda, wgt, ikinks, idebug )
%
% compute parametric cubic smoothing spline (y(s),z(s)) using PP-form method

if (nargin<7 | isempty(idebug))
   idebug = 0;
end
use_knot = 1; % using free ends (0) or not-a-knot boundaries (1).

% arguments have been checked: s, y, z, wgt are all [npnt,1]

% combine measurements that lie too close together

ds_min = 1e-3 * max(diff(s));
[s, y, z, wgt] = combine_close_meas(s, y, z, wgt, ds_min);
% disp(sprintf('Ratio hmax / hmin = %3.1e', max(diff(s)) / min(diff(s)) ));

npnt = length(s);

% process separate sections between break-points

nsec    = length(ikinks) + 1;
iranges = [ 1 , ikinks , npnt ];
ay0 = []; ay1 = []; ay2 = []; ay3 = [];
az0 = []; az1 = []; az2 = []; az3 = [];

for isec = 1 : nsec
   rg = [iranges(isec) : iranges(isec+1)];
   if (idebug>=1)
      disp(sprintf('PPspline: make section %2d, rg=[%4d:%4d], s0=%10.6f', isec, rg(1), rg(end), s(rg(1)) ));
   end

   [a0, a1, a2, a3] = make_pp_spline_section( s(rg), y(rg), lambda, wgt(rg), use_knot, idebug );
   if (isec<nsec)       % a0,a2 have npnt values, a1, a3 have nseg values
      ay0 = [ay0; a0(1:end-1)]; ay1 = [ay1; a1];
      ay2 = [ay2; a2(1:end-1)]; ay3 = [ay3; a3];
   else
      ay0 = [ay0; a0]; ay1 = [ay1; a1];
      ay2 = [ay2; a2]; ay3 = [ay3; a3];
   end

   [a0, a1, a2, a3] = make_pp_spline_section( s(rg), z(rg), lambda, wgt(rg), use_knot, idebug );
   if (isec<nsec)
      az0 = [az0; a0(1:end-1)]; az1 = [az1; a1];
      az2 = [az2; a2(1:end-1)]; az3 = [az3; a3];
   else
      az0 = [az0; a0]; az1 = [az1; a1];
      az2 = [az2; a2]; az3 = [az3; a3];
   end
end

% add a1, a3 at end-point to make all arrays equal size
% v (s)  = a0 + a1     * s +   a2     * s^2 +   a3     * s^3
% v'(s1) =      a1(s0)     + 2*a2(s0) * ds  + 3*a3(s0) * ds^2

npnt = length(s);
ds   = s(npnt) - s(npnt-1);
ay1(npnt) = ay1(npnt-1) + 2 * ds * ay2(npnt-1) + 3 * ds^2 * ay3(npnt-1);
az1(npnt) = az1(npnt-1) + 2 * ds * az2(npnt-1) + 3 * ds^2 * az3(npnt-1);
ay3(npnt) = ay3(npnt-1);
az3(npnt) = az3(npnt-1);

% create structure with final result

spl = struct('y',y, 'z',z, 's',s, 'lambda',lambda, 'wgt',wgt, 'ikinks', ikinks, 'iaccel', [], ...
             'ay0',ay0, 'ay1',ay1, 'ay2',ay2, 'ay3',ay3, ...
             'az0',az0, 'az1',az1, 'az2',az2, 'az3',az3);

end % function make_pp_spline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [a0, a1, a2, a3] = make_pp_spline_section( s, y, lambda, wgt, use_knot, idebug )

% compute smoothing spline {a0--a3} for data {s,y} with
% mis-fit weight lambda and data weights wgt

if (nargin<6 | isempty(idebug))
   idebug = 0;
end

% data   has        npnt   points    numbered ipnt = 1..npnt
% spline has nseg = npnt-1 segments  numbered iseg = 1..nseg

npnt = length(s);
nseg = npnt - 1;
% disp(sprintf('%d points, %d intervals, %d equations',npnt,nseg,npnt));

if (size(s,2)>size(s,1)), s = s'; end
if (size(y,2)>size(y,1)), y = y'; end
a3 = zeros(nseg,1);
a2 = zeros(npnt,1);
a1 = zeros(nseg,1);
a0 = zeros(npnt,1);

% compute inverse weight matrix Sigma
% Note: should this be squared?

Sigma = spdiags(1./wgt, 0, npnt, npnt);

% set relative importance of smoothness of the result

kappa = 2/3 * lambda;

if (idebug>=5)
   disp(sprintf('kappa = %6.3f, Sigma =', kappa));
   disp(full(Sigma));
end

% compute distances between points

h   = diff(s);

if (idebug>=5)
   disp('h =');
   disp(h');
end

% build R- and Q-matrices: npnt equations, including boundary conditions
%  R * a2 = Q^T * a0

R   = sparse(npnt,npnt);
Q   = sparse(npnt,npnt);

if (~use_knot | nseg<=1)  % free boundary at start
   R(1, 1) = 2*h(1);
   Q(1, 1) = 0;
else            % not-a-knot boundary
   R(1, 1) = -h(2);
   R(1, 2) =  h(1) + h(2);
   R(1, 3) = -h(1);
   Q(1, 1) =  0;
   Q(2, 1) =  0;
   Q(3, 1) =  0;
end

% equations for interior points

for i = 2 : npnt-1
   R(i,i-1) =    h(i-1);
   R(i,i)   = 2*(h(i-1) + h(i));
   R(i,i+1) =             h(i);

   Q(i-1,i) =             3 / h(i-1);
   Q(i  ,i) = -3 / h(i) - 3 / h(i-1);
   Q(i+1,i) =  3 / h(i);
end

if (~use_knot | nseg<=1)  % free boundary at end
   R(npnt,npnt) = 1;
   Q(npnt,npnt) = 0;
else
   R(npnt,npnt-2) = -h(nseg);
   R(npnt,npnt-1) =  h(nseg-1) + h(nseg);
   R(npnt,npnt  ) = -h(nseg-1);
   Q(npnt-2,npnt) = 0;
   Q(npnt-1,npnt) = 0;
   Q(npnt  ,npnt) = 0;
end

% solve system for a2 and a0

A   = kappa * Q' * Sigma * Q + R;
rhs = Q' * y;

if (idebug>=5)
   disp('R =');
   disp(full(R));
   disp('Qt =');
   disp(full(Q'));
   disp('kappa * Q'' * Sigma * Q =');
   disp(full(kappa*Q'*Sigma*Q));
   disp('A =');
   disp(full(A));
   disp('rhs =');
   disp(rhs');
end

% scale free boundaries at start & end
% note: Q ~ 1/h, wgt ~ 1/h, therefore A ~ 1/h^3

mx_D         = max(abs(diag(A)));
A(1,1)       = mx_D * A(1,1);
A(npnt,npnt) = mx_D * A(npnt,npnt);

c   = condest(A);       % condition c grows with n^4 * (hmax/hmin)^2 * f(lambda)
if (idebug>=2)
   disp(sprintf('nmeas = %d, lambda = %4.0f: cond(A) = %3.1e', npnt, lambda, c));
end
if (c>1e10)
   disp(sprintf('Warning: cond(A) = %3.1e',c));
end
if (isnan(c))
   disp('h(1:6):')
   disp(h(1:6))
   disp('R(1;5,1:6):')
   disp(R(1:5,1:6))
   disp('Q(1;5,1:6):')
   disp(Q(1:5,1:6))
end

if (1==1)

   % direct implementation:

   % imid = round(npnt/2);
   % tmp = Q' * Q;
   % disp( R(imid,imid-3:imid+3) )

   a2   = A \ rhs;

else

   if (~use_knot)
      % eliminate boundary values b_1 = b_np = 0 from equations 2 and npnt-1 for LDLt decomposition
      A(2,1) = 0;
      A(npnt-1,end) = 0;

      % use L-D-Lt factorization

      [L,D,P] = ldl(A);    % Pt * A * P = L * D * Lt
      if (any(P~=eye(npnt)))
         disp('LDLt: Non-trivial permutation');
      end

      tmp1 = L \ rhs;
      tmp2 = D \ tmp1;
      a2   = L' \ tmp2;

   else

      [L,U,P] = lu(A);    % P * A = L * U
      if (any(P~=eye(npnt)))
         disp('LU: Non-trivial permutation');
      end

        % A * a2 = rhs -->  P * A * a2 = P * rhs --> L * U * a2 = P * rhs
        % U * a2 = tmp -->  L * tmp = P * rhs
      tmp = L \ (P * rhs);
      a2  = U \ tmp;

   end

end

a0 = y - kappa * Sigma * Q * a2;

% spline parameters

a3(1:nseg) = diff(a2(1:npnt)) ./ (3*h(1:nseg));
a1(1:nseg) = diff(a0(1:npnt)) ./ h(1:nseg) ...
                - a3(1:nseg).*h(1:nseg).^2 - a2(1:nseg).*h(1:nseg);

if (idebug>=5)
   disp('a0 =');
   disp(a0');
   disp('a1 =');
   disp(a1');
   disp('a2 =');
   disp(a2');
   disp('a3 =');
   disp(a3');
end

if (0==1 & is_z)
   ix = 465;
   grad_L1 = -inv(Sigma) * (y - a0);
   mat = Q;
   % mat(ix, ix+[-5:3])
   grad_L2 =  kappa * Q * inv(R) * Q' * a0;
   % [y(ix), a0(ix), grad_L1(ix), -grad_L2(ix)]

   figure(5); clf; hold on;
   % spy(abs(mat)>1e5);
   l = plot(grad_L1+grad_L2, '-*'); % set(l([2,4]),'linestyle','--');
end

end % function make_pp_spline_section

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ s, y, z, wgt ] = combine_close_meas( s, y, z, wgt, ds_min )

% function [ s, y, z, wgt ] = combine_close_meas( s, y, z, lambda, wgt, ikinks, idebug )
%
% combine measurements that lie too close together

iter = 0;
comb = ones(size(s));
ds   = diff(s);
ix   = find(ds < ds_min);

while(~isempty(ix))
   iter = iter + 1;
   disp(sprintf('It %d: merging %d measurements with neighbouring ones', iter, length(ix)));

   % replace s, y, z by weighted averages
   s(ix)   = (comb(ix) .* s(ix) + comb(ix+1) .* s(ix+1)) ./ (comb(ix) + comb(ix+1));
   y(ix)   = (comb(ix) .* y(ix) + comb(ix+1) .* y(ix+1)) ./ (comb(ix) + comb(ix+1));
   z(ix)   = (comb(ix) .* z(ix) + comb(ix+1) .* z(ix+1)) ./ (comb(ix) + comb(ix+1));
   s(ix+1) = []; y(ix+1) = []; z(ix+1) = [];

   % replace wgt by sum of weights
   wgt(ix) = (wgt(ix) + wgt(ix+1));
   wgt(ix+1) = [];

   % compute ds anew, iterate
   ds = diff(s);
   ix = find(ds < ds_min);
end % while

end % function combine_close_meas

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

