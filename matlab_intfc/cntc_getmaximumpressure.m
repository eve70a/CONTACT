
%------------------------------------------------------------------------------------------------------------
% function [ pmax ] = cntc_getmaximumpressure(ire, icp)
%
% return the maximum normal pressure in a contact problem
%
%  pnmax         - maximum pressure in contact patch [force/area]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ pmax ] = cntc_getmaximumpressure(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_getmaximumpressure: not available for icp=%d',icp));
      return
   end

   p_pmax = libpointer('doublePtr',-1);

   calllib(libname,'cntc_getmaximumpressure', ire, icp, p_pmax);

   pmax = p_pmax.value;

end % cntc_getmaximumpressure

%------------------------------------------------------------------------------------------------------------

