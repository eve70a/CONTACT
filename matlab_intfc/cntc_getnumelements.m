
%------------------------------------------------------------------------------------------------------------
% function [ mx, my ] = cntc_getnumelements(ire, icp)
%
% return the number of elements in the potential contact area used for a contact problem,
%            length of tractions arrays
%
%  mx, my        - number of discretization elements in long/lat dirs
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ mx, my ] = cntc_getnumelements(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_getnumelements: not available for icp=%d',icp));
      return
   end

   p_mx = libpointer('int32Ptr',-1);
   p_my = libpointer('int32Ptr',-1);

   calllib(libname,'cntc_getnumelements', ire, icp, p_mx, p_my);

   mx = double(p_mx.value);
   my = double(p_my.value);

end % cntc_getnumelements

%------------------------------------------------------------------------------------------------------------

