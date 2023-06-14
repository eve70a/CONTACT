
%------------------------------------------------------------------------------------------------------------
% function [ rvalues ] = cntc_getcontactlocation(ire, icp)
%
% return the contact reference location for a contact problem (module 1 only)
%
%  rvalues(lenarr)   - contact location parameters. 
%
%  The following values are returned, if permitted by the length of rvalues:
%   1 - XCP_TR    - x-position of the contact reference point in track coordinates
%   2 - YCP_TR    - y-position of the contact reference point in track coordinates
%   3 - ZCP_TR    - z-position of the contact reference point in track coordinates
%   4 - DELTCP_TR - contact reference angle: rotation about track x-axis from the track positive z-axis to the
%                   contact positive n-axis, with sign according the right-hand rule
%
%   5 - XCP_R     - x-position of the contact reference point in rail profile coordinates
%   6 - YCP_R     - y-position of the contact reference point in rail profile coordinates
%   7 - ZCP_R     - z-position of the contact reference point in rail profile coordinates
%   8 - SCP_R     - s-parameter of the contact reference point, measured along the rail profile
%   9 - DELTCP_R  - rotation about rail x-axis from rail positive z-axis to contact positive n-axis
%
%  10 - XCP_W     - x-position of the contact reference point in wheel profile coordinates
%  11 - YCP_W     - y-position of the contact reference point in wheel profile coordinates
%  12 - ZCP_W     - z-position of the contact reference point in wheel profile coordinates
%  13 - SCP_W     - s-parameter of the contact reference point, measured along the wheel profile
%  14 - DELTCP_W  - rotation about wheel x-axis from wheel positive z-axis to contact positive n-axis
%
%  15 - XPN_TR    - x-position of the pressure center of gravity in track coordinates
%  16 - YPN_TR    - y-position of the pressure center of gravity in track coordinates
%
%  21 - XW_TR     - x-position of wheel profile marker in track coordinates
%  22 - YW_TR     - y-position of wheel profile marker in track coordinates
%  23 - ZR_TR     - z-position of wheel profile marker in track coordinates
%  24 - ROLLW_TR  - roll angle of wheel profile marker in track coordinates
%  25 - YAWW_TR   - yaw angle of wheel profile marker in track coordinates
%
%  26 - YR_TR     - y-position of rail profile marker in track coordinates
%  27 - ZR_TR     - z-position of rail profile marker in track coordinates
%  28 - ROLLR_TR  - roll angle of rail profile marker in track coordinates
%
%  The "contact reference point" is the origin of the contact local coordinate system. It is determined by
%  a heuristic rule and is centered within the contact patch in a weighted sense.
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 3: m=1, cp     - require icp>0, default 1

function [ rvalues ] = cntc_getcontactlocation(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_getcontactlocation: not available for icp=%d',icp));
      return
   end

   lenarr = 28;
   p_values = libpointer('doublePtr',zeros(lenarr,1));

   calllib(libname,'cntc_getcontactlocation', ire, icp, lenarr, p_values);

   rvalues = p_values.value;

end % cntc_getcontactlocation

%------------------------------------------------------------------------------------------------------------

