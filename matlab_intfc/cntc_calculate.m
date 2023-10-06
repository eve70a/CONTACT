
%------------------------------------------------------------------------------------------------------------
% function [ ierror ] = cntc_calculate(ire, icp, idebug)
%
% perform actual CONTACT calculation for a contact problem
% in:   integer    idebug       - show warnings (2), errors (1) or hide all messages (0)
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 5: m=*, wtd    - default icp=-1

function [ ierror ] = cntc_calculate(ire, icp, idebug)

   global libname;
   CNTC_err_allow  = -12;
   CNTC_err_profil = -32;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = -1; % default: W/R contact, all patches
   end
   if (nargin<3 | isempty(idebug))
      idebug = 1;
   end

   p_ierr = libpointer('int32Ptr',-1);

   calllib(libname,'cntc_calculate', ire, icp, p_ierr);

   ierror = double(p_ierr.value);
   if (idebug>=1 & ierror==CNTC_err_allow)
      disp(sprintf('cntc_calculate: no valid license found for CONTACT library (%d).',ierror));
   elseif (idebug>=1 & ierror==CNTC_err_profil)
      disp(sprintf('cntc_calculate: an error is found in the rail or wheel profile specification (%d).',ierror));
   elseif (idebug>=1 & ierror<0)
      disp(sprintf('cntc_calculate: an error occurred in the CONTACT library (%d).',ierror));
   elseif (idebug>=2 & ierror>0)
      disp(sprintf('cntc_calculate: potential contact may be too small: there are %d points adjacent to the boundaries.',ierror));
   end

end % cntc_calculate

%------------------------------------------------------------------------------------------------------------

