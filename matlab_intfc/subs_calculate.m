
%------------------------------------------------------------------------------------------------------------
% function [ ierror ] = subs_calculate(ire, icp, idebug)
%
% perform subsurface stress calculation for a contact problem
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 7: m=*, wtd or cp - default icp=-1

function [ ierror ] = subs_calculate(ire, icp, idebug)

   global libname;
   CNTC_err_allow  = -12;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = -1; % default: W/R contact, compute for all patches
   end
   if (nargin<3 | isempty(idebug))
      idebug = 1;
   end

   p_ierr = libpointer('int32Ptr',-1);

   calllib(libname,'subs_calculate', ire, icp, p_ierr);

   ierror = double(p_ierr.value);
   if (idebug>=1 & ierror==CNTC_err_allow)
      disp(sprintf('subs_calculate: no valid license found for CONTACT library (%d).',ierror));
   elseif (idebug>=1 & ierror<0)
      disp(sprintf('subs_calculate: an error occurred in the CONTACT library (%d).',ierror));
   end

end % subs_calculate

%------------------------------------------------------------------------------------------------------------

