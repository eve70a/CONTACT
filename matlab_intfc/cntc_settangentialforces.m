
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_settangentialforces(ire, icp, fx, fy)
%
% set the total tangential forces a contact problem
% note: fx is ignored when the F-digit is 0, fy is ignored when F=0 or 1.
%
%  fx, fy       - total tangential forces relative to fstat*fn [-]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 1: m=3, cp     - require icp>0, default 1

function [ ] = cntc_settangentialforces(ire, icp, fx, fy)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (nargin<4 | isempty(fx) | isempty(fy))
      disp('ERROR in cntc_settangentialforces: total forces fx, fy are mandatory.');
      return
   end

   calllib(libname,'cntc_settangentialforces', ire, icp, fx, fy);

end % cntc_settangentialforces

%------------------------------------------------------------------------------------------------------------

