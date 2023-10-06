
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setmetadata(ire, icp, params, values)
%
% used for configuring various metadata for a contact problem
%
%  params   - codes of the metadata to be communicated to CONTACT
%  values   - values of the metadata to be communicated to CONTACT
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 1: m=3, cp     - require icp>0, default 1

function [ ] = cntc_setmetadata(ire, icp, params, values)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_setmetadata: not available for icp=%d',icp));
      return
   end
   if (nargin<4 | isempty(params) | isempty(values))
      disp('ERROR in cntc_setmetadata: params and values are mandatory.');
      return
   end
   if (length(params) ~= length(values))
      disp('ERROR in cntc_setmetadata: params and values must provide same number of items.');
      return
   end

   calllib(libname,'cntc_setmetadata', ire, icp, length(params), params, values);

end % cntc_setmetadata

%------------------------------------------------------------------------------------------------------------

