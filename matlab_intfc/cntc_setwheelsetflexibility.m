
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setwheelsetflexibility(ire, ewheel, params) 
%
% set the description of wheelset flexibilities for a wheel-rail contact problem
%    ewheel         - type of wheelset flexibilities specification (E-digit)
%    nparam         - number of parameters provided
%    params(nparam) - depending on method that is used
%
%    E=0  : keep wheelset flexibility from previous specification, ignore params provided
%    E=1-4: no wheelset flexibility               params = []
%    E=5,7: new wheelset flexibility parameters   params = [dxwhl, dywhl, dzwhl, drollw, dyaww, dpitchw,
%                                                           vxwhl, vywhl, vzwhl, vrollw, vyaww, vpitchw],
%    E=6,8: same as 5,7, with separate values for both sides of the wheelset
%
%    dimensions:   dxwhl, dywhl, dzwhl [length],  drollw, dyaww, dpitchw [angle]
%                  vxwhl, vywhl, vzwhl [veloc],   vrollw, vyaww, vpitchw [ang.veloc]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 2: m=1, wtd    - no icp needed

function [ ] = cntc_setwheelsetflexibility(ire, ewheel, params)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(ewheel))
      ewheel = 0;
   end
   if (nargin<3 | isempty(params))
      disp('ERROR in cntc_setwheelsetflexibility: invalid params provided.');
      return
   end

   calllib(libname,'cntc_setwheelsetflexibility', ire, ewheel, length(params), params);

end % cntc_setwheelsetflexibility

%------------------------------------------------------------------------------------------------------------

