
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setwheelsetvelocity(ire, ewheel, params) 
%
% continue the wheelset description for a wheel-rail contact problem, set wheelset velocity data
%  ewheel         - type of velocity specification (E-digit)
%  nparam         - number of parameters provided
%  params(nparam) - depending on method that is used
%
%  E=0-1: no new wheelset velocity    params = [ ]
%  E=2-5:    new wheelset velocity    params = [ vx, vy, vz, vroll, vyaw, vpitch ]
%            for roller-rigs vx is replaced by rpitch (C1=4,5)
%            position increments v * dt are used in transient shifts (T=1)
%
%  dimensions:     vx, vy, vz    : [veloc]   vroll, vyaw, vpitch, rpitch  : [ang.veloc]
%                  shft_sws--zws : [length]  shft_rol, shft_yaw, shft_pit : [angle]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 2: m=1, wtd    - no icp needed

function [ ] = cntc_setwheelsetvelocity(ire, ewheel, params)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(ewheel))
      ewheel = 2;
   end
   if (nargin<3 | isempty(params))
      disp('ERROR in cntc_setwheelsetvelocity: invalid params provided.');
      return
   end

   calllib(libname,'cntc_setwheelsetvelocity', ire, ewheel, length(params), params);

end % cntc_setwheelsetvelocity

%------------------------------------------------------------------------------------------------------------

