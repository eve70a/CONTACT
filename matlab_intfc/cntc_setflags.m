
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setflags(ire, icp, params, values)
%
% used for configuring various flags for a contact problem
%
%  lenflg         - length of params/values arrays
%  params(lenflg) - codes of the parameters to be communicated to CONTACT
%  values(lenflg) - values of the parameters to be communicated to CONTACT
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 5: m=*, wtd    - default icp=-1

function [ ] = cntc_setflags(ire, icp, params, values)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = -1; % default: W/R contact, all patches
   end
   if (nargin<4 | isempty(params) | isempty(values))
      disp('ERROR in cntc_setflags: params and values are mandatory.');
      return
   end
   if (length(params)~=length(values))
      disp('ERROR in cntc_setflags: inconsistent length of params and values.');
      return
   end

   calllib(libname,'cntc_setflags', ire, icp, length(params), params, values);

end % cntc_setflags

%------------------------------------------------------------------------------------------------------------

