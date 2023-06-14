
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_finalize(ire)
%
% Finalize calculations and clean-up for a result element
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 0: m=*, glob   - no icp needed

function [ ] = cntc_finalize(ire)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end

   calllib(libname,'cntc_finalize', ire);

end % cntc_finalize

%------------------------------------------------------------------------------------------------------------

