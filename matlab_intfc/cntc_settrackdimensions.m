
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_settrackdimensions(ire, ztrack, params) 
%
% set the track or roller-rig description for a wheel-rail contact problem
%
%   !!!  This is has become an alias for cntc_settrackdimensions_new.
%        This used to link to cntc_settrackdimensions_old in the previous version. !!!
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 2: m=1, wtd    - no icp needed

function [ ] = cntc_settrackdimensions(ire, ztrack, params)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(ztrack))
      disp('ERROR in cntc_settrackdimensions: ztrack is mandatory.');
      return
   end
   if (nargin<3 | isempty(params))
      disp('ERROR in cntc_settrackdimensions: invalid params provided.');
      return
   end

   calllib(libname,'cntc_settrackdimensions', ire, ztrack, length(params), params);

end % cntc_settrackdimensions

%------------------------------------------------------------------------------------------------------------

