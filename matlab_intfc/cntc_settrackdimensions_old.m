
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_settrackdimensions_old(ire, ztrack, params) 
%
% set the track or roller-rig description for a wheel-rail contact problem
%  ztrack    - control digit ZTRACK
%  params    - depending on method that is used
%
%    1: new design track dimensions   params = [gaugwd, gaught, cant, nomrad ]
%    2: new track deviations          params = [dyrail, dzrail, drollr, vyrail, vzrail, vrollr]
%    3: new dimensions & track deviations for current side of the track
%           params = [gaugwd, gaught, cant, nomrad, dyrail, dzrail, drollr, vyrail, vzrail, vrollr]
%
% dimensions: gaugwd, gaught, nomrad, dyrail, dzrail [length],  cant, drollr [angle],
%                                     vyrail, vzrail [veloc],         vrollr [ang.veloc]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 2: m=1, wtd    - no icp needed

function [ ] = cntc_settrackdimensions_old(ire, ztrack, params)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(ztrack))
      disp('ERROR in cntc_settrackdimensions_old: ztrack is mandatory.');
      return
   end
   if (nargin<3 | isempty(params))
      disp('ERROR in cntc_settrackdimensions_old: invalid params provided.');
      return
   end

   calllib(libname,'cntc_settrackdimensions_old', ire, ztrack, length(params), params);

end % cntc_settrackdimensions

%------------------------------------------------------------------------------------------------------------

