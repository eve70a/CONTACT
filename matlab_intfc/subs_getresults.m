
%------------------------------------------------------------------------------------------------------------
% function [ table ] = subs_getresults(ire, icp, iblk, icol)
%   or
% function [ blk ] = subs_getresults(ire, icp, iblk, 'all')
%
% return selected columns of the results of the subsurface stress calculation for one block for a
% contact problem.
%
%  iblk              - number of elements in potential contact area
%  icol(ncol)        - requested columns of result data
%                       1-- 3: x,y,z        - positions of points for subsurface stress calculation
%                       4-- 6: ux,uy,uz     - elastic displacements
%                       7-- 9: sighyd,vm,tr - hydrostatic, von Mises and Tresca stresses
%                      10--12: sigma1,2,3   - principal stresses
%                      13--15: sigxx,yx,zx  - components of the full stress tensor
%                      16--18: sigxy,yy,zy
%                      19--21: sigxz,yz,zz
%  table(npnt,ncol)  - result data, with npnt = nx*ny*nz entries filled per column
%  blk               - structure with subsurface results, as used by loadstrs / plotstrs
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ table ] = subs_getresults(ire, icp, iblk, icol)
   global libname;

   totcol = 21;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in subs_getresults: not available for icp=%d',icp));
      return
   end

   if (nargin<3 | isempty(iblk))
      iblk = 1;
   end
   [ nx, ny, nz ] = subs_getblocksize(ire, icp, iblk);
   if (any([nx, ny, nz]<=0))
      disp(sprintf('ERROR in subs_getresults: no data for iblk = %d', iblk));
      disp([nx, ny, nz])
      table = [];
      return
   end

   make_struct = 0;
   if (nargin>=4 & ischar(icol))
      make_struct = 1;
      icol = [1:totcol];
   end
   if (nargin<4 | isempty(icol) | any(icol<=0) | any(icol>totcol))
      disp(sprintf('ERROR in subs_getresults: columns must be 1 <= icol <= %d',totcol));
      disp(icol)
      return
   end

   npnt = nx * ny * nz;
   ncol = length(icol);

   p_table = libpointer('doublePtr',zeros(npnt,ncol));

   calllib(libname,'subs_getresults', ire, icp, iblk, npnt, ncol, icol, p_table);

   table = p_table.value;
   % table = reshape(p_table.value, npnt, ncol);

   if (make_struct)
      blk = struct('nx',nx, 'ny',ny, 'nz',nz);

      % sort data such that x runs fastest, then y, then z
      [~, iperm] = sort(table(:,1)); table = table(iperm,:);
      [~, iperm] = sort(table(:,2)); table = table(iperm,:);
      [~, iperm] = sort(table(:,3)); table = table(iperm,:);

      % the data-lines are given with iz running fastest, ix slowest.
      blk.npoints = nx * ny * nz;
      blk.x       = reshape(table(:, 1), nx, ny, nz); blk.x = squeeze(blk.x(:,1,1));
      blk.y       = reshape(table(:, 2), nx, ny, nz); blk.y = squeeze(blk.y(1,:,1)); blk.y = blk.y';
      blk.z       = reshape(table(:, 3), nx, ny, nz); blk.z = squeeze(blk.z(1,1,:));
      blk.ux      = reshape(table(:, 4), nx, ny, nz);
      blk.uy      = reshape(table(:, 5), nx, ny, nz);
      blk.uz      = reshape(table(:, 6), nx, ny, nz);
      blk.sighyd  = reshape(table(:, 7), nx, ny, nz);
      blk.sigvm   = reshape(table(:, 8), nx, ny, nz);
      blk.sigtr   = reshape(table(:, 9), nx, ny, nz);
      blk.sigma1  = reshape(table(:,10), nx, ny, nz);
      blk.sigma2  = reshape(table(:,11), nx, ny, nz);
      blk.sigma3  = reshape(table(:,12), nx, ny, nz);
      blk.sigxx   = reshape(table(:,13), nx, ny, nz);
      blk.sigxy   = reshape(table(:,16), nx, ny, nz);
      blk.sigyy   = reshape(table(:,17), nx, ny, nz);
      blk.sigxz   = reshape(table(:,19), nx, ny, nz);
      blk.sigyz   = reshape(table(:,20), nx, ny, nz);
      blk.sigzz   = reshape(table(:,21), nx, ny, nz);
      table = blk;

   end % make_struct

end % subs_getresults

%------------------------------------------------------------------------------------------------------------

