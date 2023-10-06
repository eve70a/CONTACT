
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setrollingstepsize(ire, icp, chi, dq)
%
% set the rolling direction and step size for a contact problem
%
% in w/r contact,    icp = -1
%     chi            - ignored
%     dqrel          - rolling step size relative to grid size dx [-]
% in generic contact, icp > 0,
%     chi            - rolling direction [angle]
%     dq             - rolling step size [length]
%
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 5: m=*, wtd    - default icp=-1

function [ ] = cntc_setrollingstepsize(ire, icp, chi, dq)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = -1;
      chi =  0;
   end
   if (nargin<4 | isempty(chi) | isempty(dq))
      disp('ERROR in cntc_setrollingstepsize: chi and dq are mandatory.');
      return
   end

   calllib(libname,'cntc_setrollingstepsize', ire, icp, chi, dq);

end % cntc_setrollingstepsize

%------------------------------------------------------------------------------------------------------------

