
%------------------------------------------------------------------------------------------------------------
% function [ pen ] = cntc_getpenetration(ire, icp)
%
% return the penetration (approach) for a contact problem
%
%  pen           - penetration [length]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ pen ] = cntc_getpenetration(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_getpenetration: not available for icp=%d',icp));
      return
   end

   p_pen = libpointer('doublePtr',-1);

   calllib(libname,'cntc_getpenetration', ire, icp, p_pen);

   pen = p_pen.value;

end % cntc_getpenetration

%------------------------------------------------------------------------------------------------------------

