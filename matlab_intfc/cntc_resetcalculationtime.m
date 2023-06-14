
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_resetcalculationtime(ire, icp)
%
% reset the accumulated cpu-time and wall-clock-time used for a contact problem
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ ] = cntc_resetcalculationtime(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_resetcalculationtime: not available for icp=%d',icp));
      return
   end

   calllib(libname,'cntc_resetcalculationtime', ire, icp);

end % cntc_resetcalculationtime

%------------------------------------------------------------------------------------------------------------

