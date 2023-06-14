
%------------------------------------------------------------------------------------------------------------
% function [ nx, ny, nz ] = subs_getblocksize(ire, icp, iblk)
%
% get the number of points in a block used for subsurface stress calculation
%
%  nx, ny, nz    - number of points used in x-, y- and z-directions
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ nx, ny, nz ] = subs_getblocksize(ire, icp, iblk)
   global libname;

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

   p_nx = libpointer('int32Ptr',-1);
   p_ny = libpointer('int32Ptr',-1);
   p_nz = libpointer('int32Ptr',-1);

   calllib(libname,'subs_getblocksize', ire, icp, iblk, p_nx, p_ny, p_nz);

   nx = double(p_nx.value);
   ny = double(p_ny.value);
   nz = double(p_nz.value);

end % subs_getblocksize

%------------------------------------------------------------------------------------------------------------
