
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setwheelsetdimensions(ire, ewheel, params) 
%
% set the wheelset description for a wheel-rail contact problem
%  ewheel         - control digit EWHEEL
%  nparam         - number of parameters provided
%  params(nparam) - depending on method that is used
%
%  E=1-2, 5-6: keep old wheelset dimensions, ignore params provided
%  E=3-4, 7-8: new wheelset geometry   params = [fbdist, fbpos, nomrad]
%
%  dimensions:  fbdist, fbpos, nomrad [length]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 2: m=1, wtd    - no icp needed

function [ ] = cntc_setwheelsetdimensions(ire, ewheel, params)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<3 | isempty(ewheel) | isempty(params))
      disp('ERROR in cntc_setwheelsetdimensions: ewheel and params are mandatory.');
      return
   end

   calllib(libname,'cntc_setwheelsetdimensions', ire, ewheel, length(params), params);

end % cntc_setwheelsetdimensions

%------------------------------------------------------------------------------------------------------------

