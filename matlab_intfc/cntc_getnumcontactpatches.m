
%------------------------------------------------------------------------------------------------------------
% function [ npatch ] = cntc_getnumcontactpatches(ire)
%
% return the number of contact patches used in a w/r contact problem
%
%  npatch        - number of separate contact patches
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 2: m=1, wtd    - no icp needed

function [ npatch ] = cntc_getnumcontactpatches(ire)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end

   p_npatch = libpointer('int32Ptr',-1);

   calllib(libname,'cntc_getnumcontactpatches', ire, p_npatch);

   npatch = double(p_npatch.value);

end % cntc_getnumcontactpatches

%------------------------------------------------------------------------------------------------------------

