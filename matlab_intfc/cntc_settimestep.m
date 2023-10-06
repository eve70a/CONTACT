
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_settimestep(ire, icp, dt)
%
% set the time step size dt for a contact problem, in particular for T = 0 or 1 (shifts)
%
%  dt          - time step size [time]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 1: m=3, cp     - require icp>0, default 1

function [ ] = cntc_settimestep(ire, icp, dt)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_settimestep: not available for icp=%d',icp));
      return
   end
   if (nargin<3 | isempty(dt))
      disp('ERROR in cntc_settimestep: dt is mandatory.');
      return
   end

   calllib(libname,'cntc_settimestep', ire, icp, dt);

end % cntc_settimestep

%------------------------------------------------------------------------------------------------------------

