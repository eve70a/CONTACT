
%------------------------------------------------------------------------------------------------------------
% function [ ptmax ] = cntc_getmaximumtraction(ire, icp)
%
% return the maximum tangential traction in a contact problem
%
%  ptmax         - maximum traction |pt| in contact patch [force/area]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ ptmax ] = cntc_getmaximumtraction(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_getmaximumtraction: not available for icp=%d',icp));
      return
   end

   p_ptmax = libpointer('doublePtr',-1);

   calllib(libname,'cntc_getmaximumtraction', ire, icp, p_ptmax);

   ptmax = p_ptmax.value;

end % cntc_getmaximumtraction

%------------------------------------------------------------------------------------------------------------

