
%------------------------------------------------------------------------------------------------------------
% function [ values ] = cntc_getglobalforces(ire, icp)
%
% return the overall forces for a w/r contact problem (module 1 only)
%
% rvalues(lenarr)   - contact forces [force] and moments [force.length]
%
% The following values are returned, if permitted by the length of rvalues:
%   1 - FX_TR    - total force on the output body, component in track longitudinal x-direction
%   2 - FY_TR    - total force on the output body, component in track lateral y-direction
%   3 - FZ_TR    - total force on the output body, component in track vertical z-direction
%   4 - MX_RR_TR - total moment on output body about rail profile marker, component in track x-direction
%   5 - MY_RR_TR - total moment on output body about rail profile marker, component in track y-direction
%   6 - MZ_RR_TR - total moment on output body about rail profile marker, component in track z-direction
%
%   7 - FX_WS    - total force on the output body, component in wheelset longitudinal x-direction
%   8 - FY_WS    - total force on the output body, component in wheelset lateral y-direction
%   9 - FZ_WS    - total force on the output body, component in wheelset vertical z-direction
%  10 - MX_RW_WS - total moment on output body about wheel profile marker, component in wheelset x-direction
%  11 - MY_RW_WS - total moment on output body about wheel profile marker, component in wheelset y-direction
%  12 - MZ_RW_WS - total moment on output body about wheel profile marker, component in wheelset z-direction
%                  Note that the 'output body' is the rail when using CONTACTs unit convention
%
%  13 - FX_RR    - total force on the output body, component in rail profile x-direction
%  14 - FY_RR    - total force on the output body, component in rail profile y-direction
%  15 - FZ_RR    - total force on the output body, component in rail profile z-direction
%  16 - MX_RR_RR - total moment on output body about rail profile marker, component in rail x-direction
%  17 - MY_RR_RR - total moment on output body about rail profile marker, component in rail y-direction
%  18 - MZ_RR_RR - total moment on output body about rail profile marker, component in rail z-direction
%
%  19 - FX_RW    - total force on the output body, component in wheel profile x-direction
%  20 - FY_RW    - total force on the output body, component in wheel profile y-direction
%  21 - FZ_RW    - total force on the output body, component in wheel profile z-direction
%  22 - MX_RW_RW - total moment on output body about wheel profile marker, component in wheel x-direction
%  23 - MY_RW_RW - total moment on output body about wheel profile marker, component in wheel y-direction
%  24 - MZ_RW_RW - total moment on output body about wheel profile marker, component in wheel z-direction
%
%  25 - X_AVG    - average contact position (minimum total moment), track longitudinal x-direction
%  26 - Y_AVG    - average contact position (minimum total moment), track lateral y-direction
%  27 - Z_AVG    - average contact position (minimum total moment), track vertical z-direction
%  28 - MX_AV_TR - total moment on output body about average contact position, comp. in track x-direction
%  29 - MY_AV_TR - total moment on output body about average contact position, comp. in track y-direction
%  30 - MZ_AV_TR - total moment on output body about average contact position, comp. in track z-direction
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 4: m=1, wtd/cp - default icp=-1

function [ values ] = cntc_getglobalforces(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = -1;   % default: sum over all contact patches
   end

   lenarr = 30;
   p_values = libpointer('doublePtr',zeros(lenarr,1));

   calllib(libname,'cntc_getglobalforces', ire, icp, lenarr, p_values);

   values = p_values.value;

end % cntc_getglobalforces

%------------------------------------------------------------------------------------------------------------

