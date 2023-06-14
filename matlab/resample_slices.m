
function [ slcs, ierror ] = resample_slices( slcs, ds_max2d, idebug, show_fig, fig_ofs )

% function [ slcs, ierror ] = resample_slices( slcs, ds_max2d, [idebug], [show_fig], [fig_ofs] )
%
% resample the slices on the basis of the feature information, 
% using the same number of points per part in each slice.
%
% in:  slcs.prr    - profiles per slice, with different #points(islc) and arc-length L(islc)
%      slcs.s_feat - break points per slice [ s0, s1, .., sn ] corresponding to 'geometric features'
%      ds_max2d    - maximum step size ds after resampling
%
% s_feat(islc,:) = [ -1, -1, 5.0, 10.3, 20.0, 999.0, -1, -1 ]
%         -->  this has 8 breaks, 7 parts
%         -->  for this slice 'islc', breaks 1, 2, 7, 8 are non-existent
%         -->  break 3 == 5.0 means that 5 mm will be trimmed off from the start of the profile
%         -->  break 6 == 999.0 will be changed to s(end), the end-point of the profile
 
% method:
%  - for this slice, len_part(islc,:) = [ NaN, NaN, 5.3, 9.7, s(end)-20.0, NaN, NaN ]
%  - the maximum max(len_part) is computed over all slices,
%                         say max_len = [ 5.2, 5.6, 5.8, 9.7, 15.0, 10.0, 9.0 ]
%  - a target step size ds_max2d is given, say ds_max2d = 1.0
%  - the number of points per part is set accordingly, e.g. n_p = [ 6, 6, 6, 10, 15, 10, 9 ]
%  - each slice is resampled using its own step sizes ds, e.g. ds_p = [ NaN, NaN, 5.3/(6-1), ... ]

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

   if (nargin<3 | isempty(idebug))
      idebug   = 0;
   end
   if (nargin<4)
      show_fig = [];
   end
   if (nargin<5 | isempty(fig_ofs))
      fig_ofs = 0;
   end
   is_debug = 0; % slice for which information is printed
   if (idebug>=3 & ~isempty(show_fig))
      is_debug = 8;
   end
   ierror = 0;

   s_feat = slcs.s_feat;
   nfeat  = size(s_feat, 2);
   nslc   = size(s_feat, 1);

   % check size of s_feat array

   if (nslc~=slcs.nslc | nfeat<=1)
      disp(sprintf('Incorrect break points s_feat, needs at least 2 s-values for each of %d slices', ...
                                                                                slcs.nslc));
      disp(size(s_feat));
      return;
   end

   % check and adapt the feature information

   for is = 1 : nslc

      % - mark non-used breaks (features) for this slice by NaN-value

      ib  = find(s_feat(is,:) < 0);
      s_feat(is,ib) = NaN;

      % - determine active 'breaks' [ib0:ib1] for each slice

      ib0 = find(s_feat(is,:) >= 0, 1, 'first');
      ib1 = find(s_feat(is,:) >= 0, 1, 'last');
      slcs.slc_ib(is,1:2) = [ib0, ib1];

      % - shift s_feat from [0, L] as used in the parts-file to [s_0, s_end] as used in the profile

      s_feat(is,ib0:ib1) = s_feat(is,ib0:ib1) + slcs.prr(is).ProfileS(1);

      % - clip last value after s_end

      s_feat(is,ib1) = min(s_feat(is,ib1), slcs.prr(is).ProfileS(end));

      % - check that s_feat are in strictly increasing order and in range of profile

      if (any(diff(s_feat(is,ib0:ib1))<=0))
         disp(sprintf('Incorrect break points s_feat for slice %d, must be strictly increasing.',is))
         disp([ib0, ib1])
         disp(s_feat(is,:))
         m = min(diff(s_feat(is,ib0:ib1)));
         disp(sprintf('Minimum s_p - s_{p-1} = %3.1e', min(m)));
         return;
      end

      % - determine the length of each part in this slice

      len_part(is,ib0:ib1-1) = diff(s_feat(is,ib0:ib1));

   end % for is: adapt parts-info

   % determine the maximum length per part over all of the slices

   L_part = max( len_part );

   % check scaling

   tot_len = sum(L_part);
   if (tot_len < 5*ds_max2d)
      disp(sprintf('Warning: step ds_max2d = %5.2f may be too large for total length = %5.2f', ...
                ds_max2d, tot_len));
   end

   % determine number of points and step size du for each of the parts

   nseg_p = ceil( L_part / ds_max2d );

   % determine the first point number i_p for each of the parts after resampling

   i_p    = cumsum([1, nseg_p]);

   % set total number of points after resampling

   n_pnt  = i_p(end);

   % set the sampling positions vj

   vj  = [1 : n_pnt];

   if (idebug>=2 | n_pnt<=4)
      for ipart = 1 : nfeat-1
         disp(sprintf('part %d: %3d segments, points [%3d,%3d]', ipart, nseg_p(ipart), i_p(ipart+[0,1])));
      end
   end

   % create output arrays for resampled surface

   slcs.npnt   = n_pnt;
   slcs.iseg_p = i_p;
   slcs.vj     = vj;
   slcs.mask_j =      zeros(slcs.nslc, n_pnt);
   slcs.xsurf  = NaN * ones(slcs.nslc, n_pnt);
   slcs.ysurf  = NaN * ones(slcs.nslc, n_pnt);
   slcs.zsurf  = NaN * ones(slcs.nslc, n_pnt);

   % create mask array for use in spline computation

   for is = 1 : slcs.nslc
      ib0 = slcs.slc_ib(is,1);  % active features [ib0:ib1]
      ib1 = slcs.slc_ib(is,2);
      i0 = slcs.iseg_p( ib0 );  % active points [i0:i1]
      i1 = slcs.iseg_p( ib1 );
      slcs.mask_j(is, i0:i1) = 1;
   end

   % loop over slices, resample each slice

   for is = 1 : slcs.nslc

      slc = slcs.prr(is);
      if (idebug>=3 & is==is_debug)
         disp(sprintf(['Slice %d: s=[%5.1f,%5.1f], s_p={%5.1f,',...
                       '%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f}'], ...
                       is, slcs.prr(is).ProfileS([1,end]), s_feat(is,:)));
      end

      % - determine breaks [ib0:ib1] available for this slice

      ib0 = find(~isnan(s_feat(is,:)), 1, 'first');
      ib1 = find(~isnan(s_feat(is,:)), 1, 'last');

      % - set the positions s_j corresponding to u_j

      sj     = NaN * ones(n_pnt,1);
      for ipart = ib0 : ib1-1

         i0  = i_p(ipart);
         i1  = i_p(ipart+1);
         v_len = vj(i1) - vj(i0);
         s_ofs = s_feat(is,ipart);
         s_len = s_feat(is,ipart+1) - s_feat(is,ipart);
         sj(i0:i1) = s_ofs + s_len * (vj(i0:i1) - vj(i0)) / v_len;

         if (idebug>=3 & is==is_debug)
            disp(sprintf(' part %d: s=[%6.1f,%6.1f], sj={%6.1f,%6.1f ..%6.1f}', ...
                        ipart, s_feat(is,ipart+[0:1]), sj([i0,i0+1,i1]) ));
         end
      end

      % resample at selected sj-positions

      i0  = i_p(ib0);
      i1  = i_p(ib1);

      if (exist('make_spline'))

         % make spline for the slice using Vtech internal function, evaluate at positions sj

         lambda = 0; wgt = []; ikinks = []; iaccel = []; use_bspline = 0; ds_bspline = 1;
         spl = make_spline( slc.ProfileS, slc.ProfileY, slc.ProfileZ, lambda, wgt, ikinks, iaccel, ...
                                                                           use_bspline, ds_bspline);
         [~, yj, zj] = eval_spline(spl, sj(i0:i1));
         slcs.xsurf(is,i0:i1) = slcs.s(is);
         slcs.ysurf(is,i0:i1) = yj;
         slcs.zsurf(is,i0:i1) = zj;

      else
         % linear interpolation. Note: s_i(end) can be larger than ProfileS(end) by round-off error

         slcs.xsurf(is,:) = slcs.s(is);
         slcs.ysurf(is,:) = interp1( slc.ProfileS, slc.ProfileY, sj, 'linear', 'extrap' )';
         slcs.zsurf(is,:) = interp1( slc.ProfileS, slc.ProfileZ, sj, 'linear', 'extrap' )';
      end

      nnan = nnz( isnan(slcs.ysurf(is,i0:i1))+isnan(slcs.zsurf(is,i0:i1)) );
      if (nnan>0)
         disp(sprintf('resample_slices: slice %3d has %d nan-values',is, nnan));
         ierror = 1;
      end

      if (idebug>=3 & is==is_debug)
         if (any(show_fig==2))
            figure(2+fig_ofs); clf; hold on;
            plot(slcs.ysurf(is,:));
            plot(i_p, slcs.ysurf(is,i_p), 'o');
            grid on;
            xlabel('point j'); ylabel('s [mm]');
            title(sprintf('Slice %d: resampled y-coordinate',is));
         end

         if (any(show_fig==3))
            figure(3+fig_ofs); clf; hold on;
            plot(slcs.ysurf(is,:), slcs.zsurf(is,:));
            plot(slcs.ysurf(is,i_p), slcs.zsurf(is,i_p), 'o');
            set(gca,'ydir','reverse');
            grid on;
            xlabel('point j'); ylabel('s [mm]');
            title(sprintf('Slice %d: resampled profile',is));
         end
      end
   end

   if (any(show_fig==4))
      figure(4+fig_ofs); clf;
      surf(slcs.xsurf', slcs.ysurf', slcs.zsurf');
      set(gca, 'ydir','reverse', 'zdir','reverse');
      % v = axis; v(1:2) = slcs.s([5,end-4]); axis(v);
      v = get(gca,'dataaspectratio'); v(3) = v(2); set(gca,'dataaspectratio',v);
      shading flat; colorbar;
      xlabel('s_{fc} [mm]');
      ylabel('y_{r} [mm]');
      zlabel('z_{r} [mm]');
   end

end % function resample_slices
