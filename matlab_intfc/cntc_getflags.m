
%------------------------------------------------------------------------------------------------------------
% function [ values ] = cntc_getflags(ire, icp, params)
%
% used for retrieving various configuring flags from a contact problem
%
%  lenflg         - length of params/values arrays
%  params(lenflg) - codes of the parameters to be obtained from CONTACT
%  values(lenflg) - values of the parameters obtained from CONTACT
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 5: m=*, wtd    - default icp=-1

function [ values ] = cntc_getflags(ire, icp, params)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = -1; % default: W/R contact, all patches
   end
   if (nargin<3 | isempty(params))
      disp('ERROR in cntc_getflags: params is mandatory.');
      return
   end

   lenarr = length(params);
   p_values = libpointer('int32Ptr',zeros(lenarr,1));

   calllib(libname,'cntc_getflags', ire, icp, lenarr, params, p_values);

   % convert to double for consistency with loadcase.m
   % note that for int32-arrays, NaN == 0
   values = double(p_values.value);

end % cntc_getflags

%------------------------------------------------------------------------------------------------------------

