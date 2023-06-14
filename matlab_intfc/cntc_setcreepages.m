
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setcreepages(ire, icp, vx, vy, phi)
%
% set the kinematic constants (creepages) for a contact problem
% note: vx is ignored when the F-digit is 1 or 2, vy is ignored when F=1.
%
%  vx, vy, phi    - in rolling, T=2,3: long/lat/spin creepages [-, -, angle/length]
%                   in shifts,  T=1:   long/lat/spin shift     [length, length, angle]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 1: m=3, cp     - require icp>0, default 1

function [ ] = cntc_setcreepages(ire, icp, vx, vy, phi)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_setcreepages: not available for icp=%d',icp));
      return
   end
   if (nargin<5 | isempty(vx) | isempty(vy) | isempty(phi))
      disp('ERROR in cntc_setcreepages: creepages vx,vy,phi are mandatory.');
      return
   end

   calllib(libname,'cntc_setcreepages', ire, icp, vx, vy, phi);

end % cntc_setcreepages

%------------------------------------------------------------------------------------------------------------

