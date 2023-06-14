
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_settrackdimensions_new(ire, ztrack, params) 
%
% set the track or roller-rig description for a wheel-rail contact problem
%  ztrack    - control digit ZTRACK
%  params    - depending on method that is used
%
%    1: new design track dimensions   params = [gaught, gaugsq, gaugwd, cant, nomrad],   if gaught >  0,
%                                         or   [gaught, raily0, railz0, cant, nomrad],   if gaught <= 0.
%    2: new track deviations          params = [dyrail, dzrail, drollr, vyrail, vzrail, vrollr]
%    3: new dimensions & track deviations for current side of the track
%                                     params(1:5) cf. Z=1 followed by params(6:11) cf. Z=2.
%
% dimensions: gaught, gaugwd, raily0, railz0, nomrad, dyrail, dzrail [length],    cant, drollr [angle]
%                                                     vyrail, vzrail [veloc],           vrollr [ang.veloc]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 2: m=1, wtd    - no icp needed

function [ ] = cntc_settrackdimensions_new(ire, ztrack, params)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(ztrack))
      disp('ERROR in cntc_settrackdimensions_new: ztrack is mandatory.');
      return
   end
   if (nargin<3 | isempty(params))
      disp('ERROR in cntc_settrackdimensions_new: invalid params provided.');
      return
   end

   calllib(libname,'cntc_settrackdimensions_new', ire, ztrack, length(params), params);

end % cntc_settrackdimensions

%------------------------------------------------------------------------------------------------------------

