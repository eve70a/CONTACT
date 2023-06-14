
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setwheelsetposition(ire, ewheel, params) 
%
% set the wheelset position state data for a wheel-rail contact problem
%  ewheel         - type of position specification (E-digit)
%  nparam         - number of parameters provided
%  params(nparam) - depending on method that is used
%
%  E=0  : keep wheelset position from previous specification, ignore params provided
%  E=1-8: new wheelset position   params = [s, y, z, roll, yaw, pitch]
%
%  dimensions:   s_ws, y_ws, z_ws [length],       roll, yaw, pitch [angle]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 2: m=1, wtd    - no icp needed

function [ ] = cntc_setwheelsetposition(ire, ewheel, params)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(ewheel))
      ewheel = 1;
   end
   if (nargin<3 | isempty(params))
      disp('ERROR in cntc_setwheelsetposition: invalid params provided.');
      return
   end

   calllib(libname,'cntc_setwheelsetposition', ire, ewheel, length(params), params);

end % cntc_setwheelsetposition

%------------------------------------------------------------------------------------------------------------

