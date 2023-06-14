
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setreferencevelocity(ire, icp, veloc)
%
% set the rolling velocity for a contact problem
%
%  veloc          - absolute rolling velocity [veloc]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 1: m=3, cp     - require icp>0, default 1

function [ ] = cntc_setreferencevelocity(ire, icp, veloc)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_setreferencevelocity: not available for icp=%d',icp));
      return
   end
   if (nargin<3 | isempty(veloc))
      disp('ERROR in cntc_setreferencevelocity: veloc is mandatory.');
      return
   end

   calllib(libname,'cntc_setreferencevelocity', ire, icp, veloc);

end % cntc_setreferencevelocity

%------------------------------------------------------------------------------------------------------------

