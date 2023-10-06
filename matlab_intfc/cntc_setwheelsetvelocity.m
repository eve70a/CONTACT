
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setwheelsetvelocity(ire, ewheel, params) 
%
% continue the wheelset description for a wheel-rail contact problem
%  ewheel         - type of velocity specification (E-digit)
%  nparam         - number of parameters provided
%  params(nparam) - depending on method that is used
%
%  E=0-1: keep velocity settings from previous specification, ignore params provided
%  E=2-8: new wheelset velocity   params = [vs, vy, vz, vroll, vyaw, vpitch]
%
%  dimensions:  vs_ws, vy_ws, vz_ws [veloc],     vroll, vyaw, vpitch [ang.veloc]
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

