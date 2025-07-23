
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_settrackdimensions(ire, ztrack, params) 
%
% set the track or roller-rig description for a wheel-rail contact problem
%  ztrack - control digit ZTRACK      params - depending on method that is used
%    0: maintain track dimensions     params = [ ]
%    1: new design track dimensions   params = [gaught,  dummy, gaugwd, cant, nomrad, curv],  if gaught >  0,
%                                         or   [gaught, raily0, railz0, cant, nomrad, curv],  if gaught <= 0.
%    2: new track deviations          params = [dyrail, dzrail, drollr, vyrail, vzrail, vrollr]
%    3: new dimensions & track deviations for current side of the track
%                                     params = params(1:6) cf. Z=1 followed by params(7:12) cf. Z=2;
%  ztrack >= 30 is used to configure the massless rail model, F=3
%   32/33: "F=3, Z=2/3":              params = [kyrail, fyrail, kzrail, fzrail, {dystep0} ]
%
% dimensions: gaught, gaugwd, raily0, railz0, nomrad  [length],  cant  [angle],  curv [1/length]
%             dyrail, dzrail  [length],  drollr  [angle],  vyrail, vzrail  [veloc], vrollr  [ang.veloc]
%             kyrail, kzrail  [force/length],  fyrail, fzrail  [force],  dystep0  [length]
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
   if (nargin<3)
      params = [];
   end

   calllib(libname,'cntc_settrackdimensions', ire, ztrack, length(params), params);

end % cntc_settrackdimensions

%------------------------------------------------------------------------------------------------------------

