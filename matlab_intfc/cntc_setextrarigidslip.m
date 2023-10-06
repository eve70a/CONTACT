
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setextrarigidslip(ire, icp, wx, wy)
%
% set the extra term of the rigid slip for all elements in the potential contact area for a contact problem
%
%  mx, my               - number of elements in potential contact area
%  wx(my,mx), wy(my,mx) - in rolling, T=2,3: extra relative rigid slip  [-]
%                         in shifts,  T=1:   extra rigid shift distance [length]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 1: m=3, cp     - require icp>0, default 1

function [ ] = cntc_setextrarigidslip(ire, icp, wx, wy)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_setextrarigidslip: not available for icp=%d',icp));
      return
   end
   if (nargin<4 | isempty(wx) | isempty(wy))
      disp('ERROR in cntc_setextrarigidslip: wx, wy are mandatory.');
      return
   end

   [ mx, my ] = cntc_getnumelements(ire, icp);

   if (any(size(wx)~=[my mx]) | any(size(wy)~=[my mx]))
      disp('ERROR in cntc_setextrarigidslip: arrays wx, wy must have size (my,mx).');
      return
   end

   wx = reshape(wx', lenarr, 1);
   wy = reshape(wy', lenarr, 1);

   calllib(libname,'cntc_setextrarigidslip', ire, icp, lenarr, wx, wy);

end % cntc_setextrarigidslip

%------------------------------------------------------------------------------------------------------------

