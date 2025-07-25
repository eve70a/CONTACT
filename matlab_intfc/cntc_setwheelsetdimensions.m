
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setwheelsetdimensions(ire, ewheel, params) 
%
% set the wheelset description for a wheel-rail contact problem
%  ewheel         - control digit EWHEEL
%  nparam         - number of parameters provided
%  params(nparam) - depending on method that is used
%
%  E=0-2, 4: no new geometry         params = []
%  E=3, 5:   new wheelset geometry   params = [fbdist, fbpos, nomrad, {ytape}]
%
%  dimensions:  fbdist, fbpos, nomrad, ytape [length]
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
   if (nargin<2 | isempty(ewheel))
      disp('ERROR in cntc_setwheelsetdimensions: ewheel is mandatory.');
      return
   end
   if (nargin<3 | isempty(params))
      params = [];
   end

   calllib(libname,'cntc_setwheelsetdimensions', ire, ewheel, length(params), params);

end % cntc_setwheelsetdimensions

%------------------------------------------------------------------------------------------------------------

