
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_settrackdimensions(ire, ztrack, params) 
%
% set the track or roller-rig description for a wheel-rail contact problem
%  ztrack    - control digit ZTRACK
%  params    - depending on method that is used
%
%    0: maintain track dimensions     params = [ ]
%    1: new design track dimensions   params = [gaught, gaugsq, gaugwd, cant, nomrad],   if gaught >  0,
%                                         or   [gaught, raily0, railz0, cant, nomrad],   if gaught <= 0.
%    2: new track deviations          params = [dyrail, dzrail, drollr, vyrail, vzrail, vrollr]
%    3: new dimensions & track deviations for current side of the track
%                                     params = params(1:5) cf. Z=1 followed by params(6:11) cf. Z=2;
%                                              additionally, [kyrail, fyrail, kzrail, fzrail] when F=3.
%
% dimensions: gaught, gaugwd, raily0, railz0, nomrad, dyrail, dzrail [length],    cant, drollr [angle]
%                                                     vyrail, vzrail [veloc],           vrollr [ang.veloc]
%                                                kyrail, kzrail [force/length], fyrail, fzrail [force]
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

