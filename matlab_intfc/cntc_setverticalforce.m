
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setverticalforce(ire, fz)
%
% set the total vertical force between the contacting bodies for a w/r contact problem (module 1)
% Note: this function sets control digit N = 1
%
%  fz             - total vertical force between the two bodies [force]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 2: m=1, wtd    - no icp needed

function [ ] = cntc_setverticalforce(ire, fz)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(fz))
      disp('ERROR in cntc_setverticalforce: fz is mandatory.');
      return
   end

   calllib(libname,'cntc_setverticalforce', ire, fz);

end % cntc_setverticalforce

%------------------------------------------------------------------------------------------------------------

