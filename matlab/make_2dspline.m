
function [ spl2d ] = make_2dspline( ui, vj, xij, yij, zij, mask_j, use_approx, use_insert, ...
                                                             use_cylindr, idebug, show_fig, fig_ofs )

% function [ spl2d ] = make_2dspline( ui, vj, xij, yij, zij, [mask_j], [use_approx], [use_insert],
%                                                            [use_cylindr], [idebug], [show_fig], [fig_ofs] )
% can be used with ui==xi, xij==[]
%
% compute parametric 2D spline {tui, tvj, cij_x, cij_y, cij_z} (B-form) for data {xij,yij,zij},
%   ui      = parameter (data) positions in longitudinal (u-) direction
%   vj      = parameter (scaled data) positions in lateral (v-) direction
%   mask_j  = mask for inactive points (0) esp. for 'sections' with shorter slices
%   use_approx  = select approximating spline (1, default) or interpolating spline (0) for x-direction
%   use_insert  = insert knots at not-a-knot positions (1, default) or not (0)
%   use_cylindr = define spline on u \in [-pi,pi), using wrap-around in evaluation

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   k      = 4;
   nmeasu = length(ui);
   nmeasv = length(vj);

   if (nargin<6 | isempty(mask_j))
      mask_j = ones(nmeasu,nmeasv);
   end
   if (nargin<7 | isempty(use_approx))
      use_approx = 1;
   end
   if (nargin<8 | isempty(use_insert))
      use_insert = 1;   % option for pure tensor splines: do not(0)/do(1) fill in not-a-knot positions
   end
   if (nargin<9 | isempty(use_cylindr))
      use_cylindr = 0;
   end
   if (nargin<10 | isempty(idebug))
      idebug = 0;
   end
   if (nargin<11 | isempty(show_fig))
      show_fig = [];
   end
   if (nargin<12 | isempty(fig_ofs))
      fig_ofs = 0;
   end

   if (~all(mask_j))
      use_insert = 1;   % filling-in not-a-knots is required when a non-trivial mask is given
   end

   if (nmeasu<k | nmeasv<k)
      disp(sprintf('Error: spline needs >= k=%d data positions in u (nmeasu=%d) and v (nmeasv=%d)', ...
                        k, nmeasu, nmeasv));
      return;
   end

   if (idebug>=2)
      disp(' ');
      if (use_approx)
         disp(sprintf('make_2dspline: nmeasu = %d, nmeasv = %d, approximating spline', nmeasu, nmeasv));
      else
         disp(sprintf('make_2dspline: nmeasu = %d, nmeasv = %d, interpolating spline', nmeasu, nmeasv));
      end
   end

   has_xij = (~isempty(xij));
   has_yij = (~isempty(yij));
   has_zij = (~isempty(zij));

   if (size(ui,1)==1), ui = ui'; end         % make column vector ui
   if (size(vj,2)==1), vj = vj'; end         % make row vector vj
   if (has_xij & size(xij,1)~=nmeasu), xij = xij'; end % make xij (nmeasu, nmeasv)
   if (has_yij & size(yij,1)~=nmeasu), yij = yij'; end % make yij (nmeasu, nmeasv)
   if (has_zij & size(zij,1)~=nmeasu), zij = zij'; end % make zij (nmeasu, nmeasv)

   if ((has_xij & any(size(xij)~=[nmeasu nmeasv])) | ...
       (has_yij & any(size(yij)~=[nmeasu nmeasv])) | ...
       (has_zij & any(size(zij)~=[nmeasu nmeasv])))

      disp('Arrays ui, vj, xij, yij, zij have incompatible sizes');
      disp( [size(ui); size(vj); size(xij); size(yij); size(zij)] );
      return;
   end

   % developed originally with use_insert=1, n_ins=2; skip knot-insertion when use_insert=0.

   if (use_insert)
      n_ins = 2;
   else
      n_ins = 0;
   end
   n_hlf = n_ins / 2;

   % determine master knot-vector tui_full with knots at all measurement points + extension at start/end

   use_insert = 1;
   use_repl = 0;
   tui_full = make_knot_vector_atmeas( ui, 'u', use_insert, use_repl, idebug );

   % determine master knot-vector tvj_full with knots at all measurement points + extension at start/end

   use_repl = 0;
   tvj_full = make_knot_vector_atmeas( vj, 'v', use_insert, use_repl, idebug );

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Phase 1: build B-splines per slice in lateral direction v -> (y,z)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % consider contiguous 'sections' of slices i0:i1 with same mask per slice

   ci_x = zeros(nmeasu, nmeasv+n_ins);
   ci_y = zeros(nmeasu, nmeasv+n_ins);
   ci_z = zeros(nmeasu, nmeasv+n_ins);

   i1 = 0;
   while(i1 < nmeasu)

      % select contiguous slices i0:i1 with the same mask per slice

      i0 = i1 + 1;
      i1 = i0;
      while (i1<nmeasu & all(mask_j(i1+1,:)==mask_j(i0,:)))
         i1 = i1 + 1;
      end
      ni  = i1 - i0 + 1;

      % process slices i0:i1 in one go with ni right hand sides

      j0  = find(mask_j(i0,:), 1, 'first');
      j1  = find(mask_j(i0,:), 1, 'last');
      nj  = j1 - j0 + 1;

      if (idebug>=1)
         disp(sprintf('compute  ci_[xyz] for ix=[%3d,%3d], active jv=[%3d,%3d]', i0,i1, j0,j1));
      end

      % select local knot-vector tvj from tvj_full, skipping not-a-knot positions

      jof = [ 1+[-3:0], 3:nj-2, nj+[0:3] ];  % jof = numbering within selection
      jk  = j0 + jof  + k-2;                 % jk  = converted to overall numbering
      jk1 = j0 + 2    + k-2;
      jk2 = j0 + nj-1 + k-2;
      tvj = tvj_full(jk);

      % disp('knots tvj:')
      % disp(tvj)

      % determine collocation matrix for v-direction: evaluate each B-spline at each measurement location

      [~, ~, ~, Bmat] = eval_bspline_basisfnc( tvj, vj(j0:j1) );

      % determine 1d spline coefficients ci_x, ci_y, ci_z, solving B [ci_x, ci_y, ci_z] = [xij, yij, zij]

      cnd = condest(Bmat);
      if (cnd>1e10)
         disp(sprintf('make_2dspline: Bmat(v) has size %3d x %3d, condition %6.2e', size(Bmat), cnd));
      end

      % ci_y: [ nmeasu, nmeasv ], Bmat: [ nmeasv, nmeasv ], yij: [ nmeasu, nmeasv ]

      % disp('inp_c_y:'); disp(yij(i0:i1,j0:j1))

      if (has_xij), c_x = Bmat \ xij(i0:i1,j0:j1)'; end
      if (has_yij), c_y = Bmat \ yij(i0:i1,j0:j1)'; end
      if (has_zij), c_z = Bmat \ zij(i0:i1,j0:j1)'; end

      % disp('out_c_y:'); disp(c_y)

      if (use_insert)

         % insert knots at not-a-knot positions

         if (has_xij)
            [ t1, c_x ] = insert_knot( tvj, tvj_full(jk1), c_x, k, idebug );
            [ t2, c_x ] = insert_knot( t1, tvj_full(jk2), c_x, k, idebug );
         end
         if (has_yij)
            [ t1, c_y ] = insert_knot( tvj, tvj_full(jk1), c_y, k, idebug );
            [ t2, c_y ] = insert_knot( t1, tvj_full(jk2), c_y, k, idebug );
         end
         if (has_zij)
            [ t1, c_z ] = insert_knot( tvj, tvj_full(jk1), c_z, k, idebug );
            [ t2, c_z ] = insert_knot( t1, tvj_full(jk2), c_z, k, idebug );
         end

      end

      % store results

      if (has_xij), ci_x(i0:i1, j0:j1+n_ins) = c_x'; end
      if (has_yij), ci_y(i0:i1, j0:j1+n_ins) = c_y'; end
      if (has_zij), ci_z(i0:i1, j0:j1+n_ins) = c_z'; end

   end

   if (any(show_fig==4) & has_zij)
      figure(4+fig_ofs); clf;
      surf(ci_z');
      view([90 -90]);
      shading flat;
      colorbar;
      xlabel('u');
      ylabel('v');
   end

   if (any(show_fig==5) & has_zij)
      figure(5+fig_ofs); clf;
      plot(vj, ci_z(5:6,2:end-1), '-o');
      xlabel('v');
      grid on;
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Phase 2: build B-splines for longitudinal direction x -> (cy,cz)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % consider 'regions' of points j0:j1 with same mask per point

   cij_x    = zeros(nmeasu+n_ins, nmeasv+n_ins);
   cij_y    = zeros(nmeasu+n_ins, nmeasv+n_ins);
   cij_z    = zeros(nmeasu+n_ins, nmeasv+n_ins);
   spl_mask = zeros(nmeasu+n_ins, nmeasv+n_ins);

   j1 = 0;
   while(j1 < nmeasv)

      % select 'interpolation paths' j0:j1 with the same mask per point

      j0 = j1 + 1;
      j1 = j0;
      while (j1<nmeasv & all(mask_j(:,j1+1)==mask_j(:,j0)))
         j1 = j1 + 1;
      end

      % determine corresponding columns jt0:jt1 in arrays ci_x, ci_y, ci_z

      if (~use_insert)
         jt0 = j0;
         jt1 = j1;
      else
         jt0 = j0 + n_hlf;
         jt1 = j1 + n_hlf;

         if (j0<=1 | nnz(mask_j(:,j0))>nnz(mask_j(:,j0-1)))
            jt0 = jt0 - 1;     % add guard band
         else
            jt0 = jt0 + 1;     % move 1st column to neighbour
         end
         if (j1>=nmeasv | nnz(mask_j(:,j1))>nnz(mask_j(:,j1+1)))
            jt1 = jt1 + 1;     % add guard band
         else
            jt1 = jt1 - 1;     % move last column to neighbour
         end
      end
      njt = jt1 - jt0 + 1;

      % loop over contiguous regions i0:i1

      i_prev = 0;
      while(i_prev<nmeasu)

          % determine contiguous active i0:i1 and corresponding rows it0:it1 in cij_x, cij_y, cij_z

          i0   = i_prev + find(mask_j(i_prev+1:end,j0), 1, 'first');
          i1   = i0-1 + find(mask_j(i0+1:end,j0)==0, 1, 'first');
          if (isempty(i1)), i1 = nmeasu; end
          ni   = i1 - i0 + 1;
          it0  = i0;
          it1  = i1 + 2;

         if (idebug>=1)
            disp(sprintf('compute cij_[xyz] for jt=[%3d,%3d], active ix=[%3d,%3d]', jt0,jt1, i0,i1));
         end
    
         if (use_approx)
    
            % approximating spline -- interpret (cx, cy, cz) as control points for 2D spline
    
            c_x = zeros(ni+n_ins,njt);
            c_y = zeros(ni+n_ins,njt);
            c_z = zeros(ni+n_ins,njt);
    
            if (~use_insert)
    
               if (has_xij), c_x(1:ni,:) = ci_x(i0:i1,jt0:jt1); end
               if (has_yij), c_y(1:ni,:) = ci_y(i0:i1,jt0:jt1); end
               if (has_zij), c_z(1:ni,:) = ci_z(i0:i1,jt0:jt1); end
    
            else
    
               % select interior 'control points' from ci_x, ci_y, ci_z
    
               if (has_xij), c_x(2:ni+1,:) = ci_x(i0:i1,jt0:jt1); end
               if (has_yij), c_y(2:ni+1,:) = ci_y(i0:i1,jt0:jt1); end
               if (has_zij), c_z(2:ni+1,:) = ci_z(i0:i1,jt0:jt1); end
    
               % add fantom control points at rows 1 and ni+2
    
               dt_ext = tui_full(3+i0)   - tui_full(3+i0-1);
               dt_int = tui_full(3+i0+1) - tui_full(3+i0);
               c_x(1,:)  = c_x(2,:) - dt_ext/dt_int * (c_x(3,:) - c_x(2,:));
               c_y(1,:)  = c_y(2,:) - dt_ext/dt_int * (c_y(3,:) - c_y(2,:));
               c_z(1,:)  = c_z(2,:) - dt_ext/dt_int * (c_z(3,:) - c_z(2,:));
               if (idebug>=1)
                  disp(sprintf('start: t_ext = %5.2f, t_bnd = %5.2f, t_int = %5.2f', tui_full(3+i0+[-1:1])));
                  disp(sprintf('       c_ext = %5.2f *c_bnd + %5.2f *c_int', [1,0]+[1,-1]*dt_ext/dt_int))
               end
    
               dt_int = tui_full(3+i1)   - tui_full(3+i1-1);
               dt_ext = tui_full(3+i1+1) - tui_full(3+i1);
               c_x(ni+2,:)  = c_x(ni+1,:) + dt_ext/dt_int * (c_x(ni+1,:) - c_x(ni,:));
               c_y(ni+2,:)  = c_y(ni+1,:) + dt_ext/dt_int * (c_y(ni+1,:) - c_y(ni,:));
               c_z(ni+2,:)  = c_z(ni+1,:) + dt_ext/dt_int * (c_z(ni+1,:) - c_z(ni,:));
               if (idebug>=1)
                  disp(sprintf('end:   t_int = %5.2f, t_bnd = %5.2f, t_ext = %5.2f', tui_full(3+i1+[-1:1])));
                  disp(sprintf('       c_ext = %5.2f *c_bnd + %5.2f *c_int', [1,0]+[1,-1]*dt_ext/dt_int))
               end
            end
    
         else
    
            % interpolating spline -- take (cx, cy, cz) as data points for 2D spline
    
            % select local knot-vector tui from tui_full
       
            iof = [ 1+[-3:0], 3:ni-2, ni+[0:3] ]; % iof = numbering within selection
            ik  = i0 + iof  + k-2;               % ik  = converted to overall numbering
            ik1 = i0 + 2    + k-2;
            ik2 = i0 + ni-1 + k-2;
            tui = tui_full(ik);
            % disp('knots tui_full:')
            % disp(tui_full')
            % disp('knots tui:')
            % disp(tui')
       
            % determine collocation matrix for u-direction: evaluate each B-spline at each measurement location
       
            [~, ~, ~, Bmat] = eval_bspline_basisfnc( tui, ui(i0:i1) );
       
            % determine 2d spline coefficients cij_[xyz], solving B [cij_x, cij_y, cij_z] = [ci_x, ci_y, ci_z]
       
            cnd = condest(Bmat);
            if (cnd>1e10)
               disp(sprintf('make_2dspline: Bmat(u) has size %3d x %3d, condition %6.2e', size(Bmat), cnd));
            end
       
            % cij_y: [ nmeasu, nmeasv+2 ], Bmat: [ nmeasu, nmeasu ], ci_y: [ nmeasu, nmeasv+2 ]
       
            if (has_xij), c_x = Bmat \ ci_x(i0:i1,jt0:jt1); end
            if (has_yij), c_y = Bmat \ ci_y(i0:i1,jt0:jt1); end
            if (has_zij), c_z = Bmat \ ci_z(i0:i1,jt0:jt1); end
       
            if (use_insert)
    
               % insert knots at not-a-knot positions
    
               if (has_xij)
                  [ t1, c_x ] = insert_knot( tui, tui_full(ik1), c_x, k, idebug );
                  [ t2, c_x ] = insert_knot( t1, tui_full(ik2), c_x, k, idebug );
               end
               if (has_yij)
                  [ t1, c_y ] = insert_knot( tui, tui_full(ik1), c_y, k, idebug );
                  [ t2, c_y ] = insert_knot( t1, tui_full(ik2), c_y, k, idebug );
               end
               if (has_zij)
                  [ t1, c_z ] = insert_knot( tui, tui_full(ik1), c_z, k, idebug );
                  [ t2, c_z ] = insert_knot( t1, tui_full(ik2), c_z, k, idebug );
               end
            end
    
         end % use_approx
    
         % store results
    
         if (has_xij), cij_x(i0:i1+n_ins, jt0:jt1) = c_x; end
         if (has_yij), cij_y(i0:i1+n_ins, jt0:jt1) = c_y; end
         if (has_zij), cij_z(i0:i1+n_ins, jt0:jt1) = c_z; end
         spl_mask(i0:i1+n_ins, jt0:jt1) = 1;

         % prepare for next section [i0:i1]

         i_prev = i1;
         if (isempty(find(mask_j(i1+1:end,j0)))), i_prev = nmeasu; end

      end % while (contiguous i0:i1)
   end

   spl2d = struct('ui',ui, 'vj',vj, 'use_cylindr',use_cylindr, 'xij',xij, 'yij',yij, 'zij',zij, ...
                           'mask_j',spl_mask, 'tui',tui_full, 'tvj',tvj_full);

   if (~use_insert)
      spl2d.tui = spl2d.tui([1:k, k+2:end-k-1, end-k+1:end]);
      spl2d.tvj = spl2d.tvj([1:k, k+2:end-k-1, end-k+1:end]);
   end

   spl2d.tui_grev = conv(spl2d.tui(2:end-1), [1 1 1]/3, 'valid');
   spl2d.tvj_grev = conv(spl2d.tvj(2:end-1), [1 1 1]/3, 'valid');

   if (has_xij), spl2d.cij_x = cij_x; end
   if (has_yij), spl2d.cij_y = cij_y; end
   if (has_zij), spl2d.cij_z = cij_z; end

   % temporarily, copy intermediate coefficients to spl2d

   if (has_xij), spl2d.ci_x = ci_x; end
   if (has_yij), spl2d.ci_y = ci_y; end
   if (has_zij), spl2d.ci_z = ci_z; end

   if (idebug>=1)

      % check correctness of interpolating spline at meas. positions

      [ ~, x_out, y_out, z_out ] = eval_2dspline( spl2d, ui, vj, [], idebug );

      disp(sprintf('Max diff(spline y - yij) = %3.1e, max(spline z - zij) = %3.1e', ...
                                        max(max(abs(y_out-yij))), max(max(abs(z_out-zij))) ))
   end

end % function make_2dspline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ tj ] = make_knot_vector_atmeas( si, namcoor, use_insert, use_repl, idebug )

% function [ tj ] = make_knot_vector_atmeas( si, namcoor, use_insert, use_repl, idebug )
% use_insert = 1 == knots at all measurements, 0 == skip 2nd and n-1th for not-a-knot b.c.
% use_repl   = 1 == repeat 1st/last points, 0 == repeat 1st/last intervals

   nmeas  = length(si);
   nspl   = nmeas;      % number of B-splines equal to number of measured values
   if (use_insert)
      nknot = nspl + 6; % cubic spline, order k=4: add 3 knots at both ends
   else
      nknot = nspl + 4; % not-a-knot b.c.: remove 2nd and n-1th meas. points
   end

   if (size(si,1)>1)
      tj     = zeros(nknot,1);
   else
      tj     = zeros(1,nknot);
   end

   dt0    = si(2) - si(1);       % compute 1st and last intervals
   dt1    = si(end) - si(end-1);

   if (use_insert & use_repl)

      tj(1:4)           = si(1);                % replicate s(1) and s(end) 4 times for order k=4, 
      tj(5:nknot-4)     = si(2:nmeas-1);        % use all measurement points
      tj(nknot-3:nknot) = si(nmeas);

   elseif (use_insert)

      tj(1:4)           = si(1)    +[-3:0]*dt0; % replicate 1st and last intervals 3 times for order k=4, 
      tj(5:nknot-4)     = si(2:nmeas-1);        % use all measurement points
      tj(nknot-3:nknot) = si(nmeas)+[ 0:3]*dt1;

   elseif (use_repl)

      tj(1:4)           = si(1);                % replicate s(1) and s(end) 4 times for order k=4, 
      tj(5:nknot-4)     = si(3:nmeas-2);        % skip s(2) and s(end-1) for not-a-knot boundaries,
      tj(nknot-3:nknot) = si(nmeas);            %    retain internal s(i) for interpolating spline

   else

      tj(1:4)           = si(1)    +[-3:0]*dt0; % replicate first & last intervals 3 times for order k=4
      tj(5:nknot-4)     = si(3:nmeas-2);        % skip s(2) and s(end-1) for not-a-knot boundaries
      tj(nknot-3:nknot) = si(nmeas)+[ 0:3]*dt1;

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

