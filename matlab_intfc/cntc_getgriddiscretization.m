
%------------------------------------------------------------------------------------------------------------
% function [ dx, dy ] = cntc_getgriddiscretization(ire, icp)
%
% get the grid discretization step sizes dx,dy for a contact problem
%
%  dx, dy        - grid discretization step sizes [length]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ dx, dy ] = cntc_getgriddiscretization(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_getgriddiscretization: not available for icp=%d',icp));
      return
   end

   p_dx = libpointer('doublePtr',-1);
   p_dy = libpointer('doublePtr',-1);

   calllib(libname,'cntc_getgriddiscretization', ire, icp, p_dx, p_dy);

   dx = p_dx.value;
   dy = p_dy.value;

end % cntc_getgriddiscretization

%------------------------------------------------------------------------------------------------------------

