
%------------------------------------------------------------------------------------------------------------
% function [ veloc ] = cntc_getreferencevelocity(ire, icp)
%
% get the rolling velocity for a contact problem
%
%  veloc          - absolute rolling velocity [veloc]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ veloc ] = cntc_getreferencevelocity(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_getreferencevelocity: not available for icp=%d',icp));
      return
   end

   p_veloc = libpointer('doublePtr',-1);

   calllib(libname,'cntc_getreferencevelocity', ire, icp, p_veloc);

   veloc = p_veloc.value;

end % cntc_getreferencevelocity

%------------------------------------------------------------------------------------------------------------

