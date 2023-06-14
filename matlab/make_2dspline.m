
function [ spl2d ] = make_2dspline( xi, uj, yij, zij, mask_j, use_approx, idebug, show_fig, fig_ofs )

% function [ spl2d ] = make_2dspline( xi, uj, yij, zij, [mask_j], [use_approx], [idebug], 
%                                                                               [show_fig], [fig_ofs] )
%
% compute parametric 2D spline {txi, tuj, cij_y, cij_z} (b-form) for data {yij,zij},
%   xi      = data positions in longitudinal (x-) direction
%   uj      = (scaled) data positions in lateral (u-) direction
%   mask_j  = mask for inactive points (0) esp. for 'sections' with shorter slices
%   use_approx = select approximating spline (1, default) or interpolating spline (0) for x-direction

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   k      = 4;
   nmeasx = length(xi);
   nmeasu = length(uj);

   if (nargin<5 | isempty(mask_j))
      mask_j = ones(nmeasx,nmeasu);
   end
   if (nargin<6 | isempty(use_approx))
      use_approx = 1;
   end
   if (nargin<7 | isempty(idebug))
      idebug = 0;
   end
   if (nargin<8 | isempty(show_fig))
      show_fig = [];
   end
   if (nargin<9 | isempty(fig_ofs))
      fig_ofs = 0;
   end

   if (nmeasx<k | nmeasu<k)
      disp(sprintf('Error: spline needs >= k=%d data positions in x (nmeasx=%d) and u (nmeasu=%d)', ...
                        k, k, nmeasx, nmeasu));
      return;
   end

   if (idebug>=2)
      disp(' ');
      if (use_approx)
         disp(sprintf('make_2dspline: nmeasx = %d, nmeasu = %d, approximating spline', nmeasx, nmeasu));
      else
         disp(sprintf('make_2dspline: nmeasx = %d, nmeasu = %d, interpolating spline', nmeasx, nmeasu));
      end
   end

   if (size(xi,1)==1), xi = xi'; end         % make column vector xi
   if (size(uj,2)==1), uj = uj'; end         % make row vector uj
   if (size(yij,1)~=nmeasx), yij = yij'; end % make yij (nmeasx, nmeasu)
   if (size(zij,1)~=nmeasx), zij = zij'; end % make zij (nmeasx, nmeasu)

   if (any(size(yij)~=[nmeasx nmeasu]) | any(size(zij)~=[nmeasx nmeasu]))
      disp('Arrays xi, uj, yij, zij have incompatible sizes');
      disp( [size(xi); size(uj); size(yij); size(zij)] );
      return;
   end

   % determine master knot-vector tuj_full with knots at all measurement points + extension at start/end

   all_meas = 1; use_repl = 0;
   tuj_full = make_knot_vector_atmeas( uj, 'u', all_meas, use_repl, idebug );

   % determine master knot-vector txi_full with knots at all measurement points + extension at start/end

   all_meas = 1; use_repl = 0;
   txi_full = make_knot_vector_atmeas( xi, 'x', all_meas, use_repl, idebug );

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Phase 1: build B-splines per slice in lateral direction u -> (y,z)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % consider 'sections' of slices ix0:ix1 with same mask per slice

   ci_y = zeros(nmeasx, nmeasu+2);
   ci_z = zeros(nmeasx, nmeasu+2);

   ix1 = 0;
   while(ix1 < nmeasx)

      % select slices ix0:ix1 with the same mask per slice

      ix0 = ix1 + 1;
      ix1 = ix0;
      while (ix1<nmeasx & all(mask_j(ix1+1,:)==mask_j(ix0,:)))
         ix1 = ix1 + 1;
      end

      % process slices ix0:ix1 in one go with nx right hand sides

      j0  = find(mask_j(ix0,:), 1, 'first');
      j1  = find(mask_j(ix0,:), 1, 'last');
      nj  = j1 - j0 + 1;

      if (idebug>=1)
         disp(sprintf('compute  ci_[syz] for ix=[%3d,%3d], active ju=[%3d,%3d]', ix0,ix1, j0,j1));
      end

      % select local knot-vector tuj from tuj_full

      jof = [ 1+[-3:0], 3:nj-2, nj+[0:3] ];
      jk  = j0 + jof  + k-2;
      jk1 = j0 + 2    + k-2;
      jk2 = j0 + nj-1 + k-2;
      tuj = tuj_full(jk);

      % disp('knots tuj:')
      % disp(tuj)

      % determine collocation matrix for u-direction: evaluate each B-spline at each measurement location

      [~, ~, ~, Bmat] = eval_bspline_basisfnc( tuj, uj(j0:j1) );

      % determine 1d spline coefficients ci_y, ci_z, solving B [ci_y, ci_z] = [yij, zij]

      cnd = condest(Bmat);
      if (cnd>1e10)
         disp(sprintf('make_2dspline: Bmat(u) has size %3d x %3d, condition %6.2e', size(Bmat), cnd));
      end

      % ci_y: [ nmeasx, nmeasu ], Bmat: [ nmeasu, nmeasu ], yij: [ nmeasx, nmeasu ]

      % disp('inp_c_y:'); disp(yij(ix0:ix1,j0:j1))

      c_y = Bmat \ yij(ix0:ix1,j0:j1)';
      c_z = Bmat \ zij(ix0:ix1,j0:j1)';

      % disp('out_c_y:'); disp(c_y)

      % insert knots at not-a-knot positions

      [ t1, c_y ] = insert_knot( tuj, tuj_full(jk1), c_y, k, idebug );
      [ t1, c_z ] = insert_knot( tuj, tuj_full(jk1), c_z, k, idebug );

      [ t2, c_y ] = insert_knot( t1, tuj_full(jk2), c_y, k, idebug );
      [ t2, c_z ] = insert_knot( t1, tuj_full(jk2), c_z, k, idebug );

      % store results

      ci_y(ix0:ix1, j0:j1+2) = c_y';
      ci_z(ix0:ix1, j0:j1+2) = c_z';

   end

   if (any(show_fig==4))
      figure(4+fig_ofs); clf;
      surf(ci_z');
      view([90 -90]);
      shading flat;
      colorbar;
      xlabel('x');
      ylabel('u');
   end

   if (any(show_fig==5))
      figure(5+fig_ofs); clf;
      plot(uj, ci_z(5:6,2:end-1), '-o');
      xlabel('u');
      grid on;
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Phase 2: build B-splines for longitudinal direction x -> (cy,cz)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % consider 'regions' of points j0:j1 with same mask per point

   cij_y    = zeros(nmeasx+2, nmeasu+2);
   cij_z    = zeros(nmeasx+2, nmeasu+2);
   spl_mask = zeros(nmeasx+2, nmeasu+2);

   j1 = 0;
   while(j1 < nmeasu)

      % select 'interpolation paths' j0:j1 with the same mask per point

      j0 = j1 + 1;
      j1 = j0;
      while (j1<nmeasu & all(mask_j(:,j1+1)==mask_j(:,j0)))
         j1 = j1 + 1;
      end

      % determine corresponding columns jt0:jt1 in arrays ci_y, ci_z

      if (j0<=1 | nnz(mask_j(:,j0))>nnz(mask_j(:,j0-1)))
         jt0 = j0;       % add guard band
      else
         jt0 = j0 + 2;   % move 1st column to neighbour
      end
      if (j1>=nmeasu | nnz(mask_j(:,j1))>nnz(mask_j(:,j1+1)))
         jt1 = j1 + 2;   % add guard band
      else
         jt1 = j1;       % move last column to neighbour
      end

      % determine active ix0:ix1 and corresponding rows it0:it1 in cij_y, cij_z

      ix0  = find(mask_j(:,j0), 1, 'first');
      ix1  = find(mask_j(:,j0), 1, 'last');
      nx   = ix1 - ix0 + 1;
      it0  = ix0;
      it1  = ix1 + 2;

      if (idebug>=1)
         disp(sprintf('compute cij_[syz] for jt=[%3d,%3d], active ix=[%3d,%3d]', jt0,jt1, ix0,ix1));
      end

      if (use_approx)

         % approximating spline -- interpret (cy, cz) as control points for 2D spline

         % select interior 'control points' from ci_y, ci_z

         njt = jt1 - jt0 + 1;
         c_y = zeros(nx+2,njt);
         c_z = zeros(nx+2,njt);

         c_y(2:nx+1,:) = ci_y(ix0:ix1,jt0:jt1);
         c_z(2:nx+1,:) = ci_z(ix0:ix1,jt0:jt1);

         % add fantom control points at rows 1 and nx+2

         dt_ext = txi_full(3+ix0)   - txi_full(3+ix0-1);
         dt_int = txi_full(3+ix0+1) - txi_full(3+ix0);
         c_y(1,:)  = c_y(2,:) - dt_ext/dt_int * (c_y(3,:) - c_y(2,:));
         c_z(1,:)  = c_z(2,:) - dt_ext/dt_int * (c_z(3,:) - c_z(2,:));
         if (idebug>=1)
            disp(sprintf('start: t_ext = %5.2f, t_bnd = %5.2f, t_int = %5.2f', txi_full(3+ix0+[-1:1])));
            disp(sprintf('       c_ext = %5.2f *c_bnd + %5.2f *c_int', [1,0]+[1,-1]*dt_ext/dt_int))
         end

         dt_int = txi_full(3+ix1)   - txi_full(3+ix1-1);
         dt_ext = txi_full(3+ix1+1) - txi_full(3+ix1);
         c_y(nx+2,:)  = c_y(nx+1,:) + dt_ext/dt_int * (c_y(nx+1,:) - c_y(nx,:));
         c_z(nx+2,:)  = c_z(nx+1,:) + dt_ext/dt_int * (c_z(nx+1,:) - c_z(nx,:));
         if (idebug>=1)
            disp(sprintf('end:   t_int = %5.2f, t_bnd = %5.2f, t_ext = %5.2f', txi_full(3+ix1+[-1:1])));
            disp(sprintf('       c_ext = %5.2f *c_bnd + %5.2f *c_int', [1,0]+[1,-1]*dt_ext/dt_int))
         end

      else

         % interpolating spline -- take (cy, cz) as data points for 2D spline

         % select local knot-vector txi from txi_full
    
         iof = [ 1+[-3:0], 3:nx-2, nx+[0:3] ];
         ik  = ix0 + iof  + k-2;
         ik1 = ix0 + 2    + k-2;
         ik2 = ix0 + nx-1 + k-2;
         txi = txi_full(ik);
         % disp('knots txi_full:')
         % disp(txi_full')
         % disp('knots txi:')
         % disp(txi')
    
         % determine collocation matrix for x-direction: evaluate each B-spline at each measurement location
    
         [~, ~, ~, Bmat] = eval_bspline_basisfnc( txi, xi(ix0:ix1) );
    
         % determine 2d spline coefficients cij_y, cij_z, solving B [cij_y, cij_z] = [ci_y, ci_z]
    
         cnd = condest(Bmat);
         if (cnd>1e10)
            disp(sprintf('make_2dspline: Bmat(x) has size %3d x %3d, condition %6.2e', size(Bmat), cnd));
         end
    
         % cij_y: [ nmeasx, nmeasu+2 ], Bmat: [ nmeasx, nmeasx ], ci_y: [ nmeasx, nmeasu+2 ]
    
         c_y = Bmat \ ci_y(ix0:ix1,jt0:jt1);
         c_z = Bmat \ ci_z(ix0:ix1,jt0:jt1);
    
         % insert knots at not-a-knot positions
    
         [ t1, c_y ] = insert_knot( txi, txi_full(ik1), c_y, k, idebug );
         [ t1, c_z ] = insert_knot( txi, txi_full(ik1), c_z, k, idebug );
    
         [ t2, c_y ] = insert_knot( t1, txi_full(ik2), c_y, k, idebug );
         [ t2, c_z ] = insert_knot( t1, txi_full(ik2), c_z, k, idebug );

      end % use_approx

      % store results

      cij_y(ix0:ix1+2, jt0:jt1) = c_y;
      cij_z(ix0:ix1+2, jt0:jt1) = c_z;
      spl_mask(ix0:ix1+2, jt0:jt1) = 1;
   end

   spl2d = struct('xi',xi, 'uj',uj, 'yij',yij, 'zij',zij, 'mask_j',spl_mask, ...
                           'txi',txi_full, 'tuj',tuj_full, 'cij_y',cij_y, 'cij_z',cij_z);

   if (idebug>=1)

      % check correctness of interpolating spline at meas. positions

      [ ~, y_out, z_out ] = eval_2dspline( spl2d, xi, uj, [], idebug );

      disp(sprintf('Max diff(spline y - yij) = %3.1e, max(spline z - zij) = %3.1e', ...
                                        max(max(abs(y_out-yij))), max(max(abs(z_out-zij))) ))
   end

end % function make_2dspline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ tj ] = make_knot_vector_atmeas( xi, namcoor, all_meas, use_repl, idebug )

% function [ tj ] = make_knot_vector_atmeas( xi, namcoor, all_meas, use_repl, idebug )
% all_meas = 1 == knots at all measurements, 0 == skip 2nd and n-1th for not-a-knot b.c.
% use_repl = 1 == repeat 1st/last points, 0 == repeat 1st/last intervals

   nmeas  = length(xi);
   nspl   = nmeas;      % number of B-splines equal to number of measured values
   if (all_meas)
      nknot = nspl + 6; % cubic spline, order k=4: add 3 knots at both ends
   else
      nknot = nspl + 4; % not-a-knot b.c.: remove 2nd and n-1th meas. points
   end

   if (size(xi,1)>1)
      tj     = zeros(nknot,1);
   else
      tj     = zeros(1,nknot);
   end

   dt0    = xi(2) - xi(1);       % compute 1st and last intervals
   dt1    = xi(end) - xi(end-1);

   if (all_meas & use_repl)

      tj(1:4)           = xi(1);                % replicate x(1) and x(end) 4 times for order k=4, 
      tj(5:nknot-4)     = xi(2:nmeas-1);        % use all measurement points
      tj(nknot-3:nknot) = xi(nmeas);

   elseif (all_meas)

      tj(1:4)           = xi(1)    +[-3:0]*dt0; % replicate 1st and last intervals 3 times for order k=4, 
      tj(5:nknot-4)     = xi(2:nmeas-1);        % use all measurement points
      tj(nknot-3:nknot) = xi(nmeas)+[ 0:3]*dt1;

   elseif (use_repl)

      tj(1:4)           = xi(1);                % replicate x(1) and x(end) 4 times for order k=4, 
      tj(5:nknot-4)     = xi(3:nmeas-2);        % skip x(2) and x(end-1) for not-a-knot boundaries,
      tj(nknot-3:nknot) = xi(nmeas);            %    retain internal x(i) for interpolating spline

   else

      tj(1:4)           = xi(1)    +[-3:0]*dt0; % replicate first & last intervals 3 times for order k=4
      tj(5:nknot-4)     = xi(3:nmeas-2);        % skip x(2) and x(end-1) for not-a-knot boundaries
      tj(nknot-3:nknot) = xi(nmeas)+[ 0:3]*dt1;

   end

   if (idebug>=3)
      disp(sprintf(' knot vector t%c has %d knots, nmeas= %d:', namcoor, nknot, nmeas));
      disp(tj')
   %  disp(tj(end-5:end)')
   end

end % function make_knot_vector_atmeas

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ tjx, a_out ] = insert_knot( tj, tnew, a_in, k, idebug )

% function [ tjx, a_out ] = insert_knot( tj, tnew, a_in, k, [idebug] )
%
% insert one new knot tnew into knot vector tj and compute B-spline coefficients a_out 

   if (nargin<5 | isempty(idebug))
      idebug = 0;
   end

   if (size(tj,1)>1), tj = tj'; end % using row vector tj

   nknot = length(tj);
   nspl  = size(a_in, 1);
   ncol  = size(a_in, 2);
   if (nknot ~= nspl+k)
      disp(sprintf('ERROR: size of tj (%d) does not match with a_in (%d, %d)', nknot, nspl, k));
      return;
   end

   % determine position jnew for tnew, form extended knot vector

   jnew = find( tj > tnew, 1, 'first');

   tjx  = [tj(1:jnew-1) , tnew, tj(jnew:end) ];

   a_out            = zeros(nspl+1,ncol);
   for j = 1 : jnew - k
      a_out(j,:) = a_in(j,:);
   end
   for j = jnew-k+1 : jnew-1
      fj         = (tjx(jnew) - tj(j)) / (tj(j+3) - tj(j));
      a_out(j,:) = (1 - fj) * a_in(j-1,:) + fj * a_in(j,:);
   end
   for j = jnew : nspl+1
      a_out(j,:) = a_in(j-1,:);
   end

end % function insert_knot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

