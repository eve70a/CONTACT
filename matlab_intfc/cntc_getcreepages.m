
%------------------------------------------------------------------------------------------------------------
% function [ vx, vy, phi ] = cntc_getcreepages(ire, icp)
%
% get the kinematic constants (creepages) for a contact problem
%
%  vx, vy, phi    - in rolling, T=2,3: long/lat/spin creepages [-, -, angle/length]
%                   in shifts,  T=1:   long/lat/spin shift [length, length, angle]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ vx, vy, phi ] = cntc_getcreepages(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_getcreepages: not available for icp=%d',icp));
      return
   end

   p_vx  = libpointer('doublePtr',-1);
   p_vy  = libpointer('doublePtr',-1);
   p_phi = libpointer('doublePtr',-1);

   calllib(libname,'cntc_getcreepages', ire, icp, p_vx, p_vy, p_phi);

   vx  = p_vx.value;
   vy  = p_vy.value;
   phi = p_phi.value;

end % cntc_getcreepages

%------------------------------------------------------------------------------------------------------------

