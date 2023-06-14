
%------------------------------------------------------------------------------------------------------------
% function [ ] = subs_addblock(ire, icp, iblk, isubs, xparam, yparam, zparam)
%
% set the parameters for a block of points for the subsurface stress calculation for a contact problem
%  iblk           - block number; all blocks with this number or higher will be discarded
%  isubs          - type of block specification
%  x/y/zparam     - parameters describing x/y/z-coordinates of the block
%
%  isubs
%    1: xparam = [                 ]  yparam = [                 ]  zparam = [ NZ, ZL, DZ   ]
%    2: xparam = [ IXL, IXINC, IXH ]  yparam = [ IYL, IYINC, IYH ]  zparam = [ NZ, ZL, DZ   ]
%    3: xparam = [ IX(i), i=1:nx   ]  yparam = [ IY(j), j=1:ny   ]  zparam = [ NZ, ZL, DZ   ]
%    5: xparam = [                 ]  yparam = [                 ]  zparam = [ Z(k), k=1:nz ]
%    6: xparam = [ IXL, IXINC, IXH ]  yparam = [ IYL, IYINC, IYH ]  zparam = [ Z(k), k=1:nz ]
%    7: xparam = [ IX(i), i=1:nx   ]  yparam = [ IY(j), j=1:ny   ]  zparam = [ Z(k), k=1:nz ]
%    9: xparam = [ X(i),  i=1:nx   ]  yparam = [ Y(j),  j=1:ny   ]  zparam = [ Z(k), k=1:nz ]
%       ix, iy [-],    x, y, z, zl, dz [length]
%
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 7: m=*, wtd or cp - default icp=-1

function [ ] = subs_addblock(ire, icp, iblk, isubs, xparam, yparam, zparam)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = -1; % default: fall-back for all contact patches
   end
   if (nargin<3 | isempty(iblk))
      iblk = 1;
   end
   if (nargin<7 | isempty(isubs) | isempty(zparam) )
      disp('ERROR in subs_addblock: isubs, zparam are mandatory.');
      return
   end
   if (any(isubs==[1,5]))
      xparam = [ -999 ]; yparam = [ -999 ];
   elseif (isempty(xparam) | isempty(yparam))
      disp(sprintf('ERROR in subs_addblock: xparam, yparam are mandatory for isubs=%d.',isubs));
      return
   end

   npx = length(xparam);
   npy = length(yparam);
   npz = length(zparam);

   calllib(libname,'subs_addblock', ire, icp, iblk, isubs, npx, npy, npz, xparam, yparam, zparam);

end % subs_addblock

%------------------------------------------------------------------------------------------------------------

