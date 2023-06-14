
%------------------------------------------------------------------------------------------------------------
% function [ rvalues ] = cntc_getwheelsetvelocity(ire)
%
% return the wheelset velocity a w/r contact problem (module 1 only)
%
%  rvalues(lenarr)   - wheelset velocity parameters. 
%
%  The following values are returned, if permitted by the length of rvalues:
%   1 - VS_WS      - s-velocity of the wheelset center of mass along the track center line
%   2 - VY_WS      - velocity in lateral y-direction of the wheelset center of mass (track coordinates)
%   3 - VZ_WS      - velocity in vertical z-direction of the wheelset center of mass (track coordinates)
%   4 - VROLL_WS   - wheelset roll velocity, angular velocity with respect to the track plane
%   5 - VYAW_WS    - wheelset yaw velocity, angular velocity with respect to the track center line x_tr
%   6 - VPITCH_WS  - wheelset pitch velocity, i.e. angular velocity about the wheelset axle
%
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 2: m=1, wtd    - no icp needed

function [ rvalues ] = cntc_getwheelsetvelocity(ire)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end

   lenarr = 6;
   p_values = libpointer('doublePtr',zeros(lenarr,1));

   calllib(libname,'cntc_getwheelsetvelocity', ire, lenarr, p_values);

   rvalues = p_values.value;

end % cntc_getwheelsetvelocity

%------------------------------------------------------------------------------------------------------------

