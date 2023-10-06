
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setpenetration(ire, icp, pen)
%
% set the approach (penetration) of the bodies as a whole for a contact problem
% Note: this function sets control digit N = 0
%
%  pen            - penetration/approach of the two bodies [length]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 1: m=3, cp     - require icp>0, default 1

function [ ] = cntc_setpenetration(ire, icp, pen)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_setpenetration: not available for icp=%d',icp));
      return
   end
   if (nargin<3 | isempty(pen))
      disp('ERROR in cntc_setpenetration: pen is mandatory.');
      return
   end

   calllib(libname,'cntc_setpenetration', ire, icp, pen);

end % cntc_setpenetration

%------------------------------------------------------------------------------------------------------------

