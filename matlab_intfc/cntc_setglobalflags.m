
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setglobalflags(params, values)
%
% used for configuring flags that are the same for all contact problems
%
%  params - codes of the parameters to be communicated to CONTACT
%  values - values of the parameters to be communicated to CONTACT
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 0: m=*, glob   - no icp needed

function [ ] = cntc_setglobalflags(params, values)
   global libname;

   if (nargin<2 | isempty(params) | isempty(values))
      disp('ERROR in cntc_setglobalflags: params and values are mandatory.');
      return
   end
   if (length(params) ~= length(values))
      disp('ERROR in cntc_setglobalflags: params and values must provide same number of items.');
      return
   end

   calllib(libname,'cntc_setglobalflags', length(params), params, values);

end % cntc_setglobalflags

%------------------------------------------------------------------------------------------------------------

