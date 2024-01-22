
%------------------------------------------------------------------------------------------------------------
% function [ rvalues ] = cntc_gethertzcontact(ire, icp)
%
% get the parameters from a Hertzian contact problem
%
% The following values are returned in rvalues:
%   1 - A1       - curvature in rolling direction [1/length]
%   2 - B1       - curvature in lateral direction [1/length]
%   3 - AA       - semi-axis in rolling direction [length]
%   4 - BB       - semi-axis in lateral direction [length]
%   5 - RHO      - effective radius of curvature, 2 / (A1 + B1) [length]
%   6 - CP       - effective semi-axis, sqrt(AA * BB) [length]
%   7 - SCALE    - potential contact scale factor [-]
%   8 - BNEG     - semi-axis of negative half-ellipse in lateral direction [length]
%   9 - BPOS     - semi-axis of positive half-ellipse in lateral direction [length]
%  10 - AOB      - ellipticity AA/BB [-]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 1: m=3, cp     - require icp>0, default 1

function [ rvalues ] = cntc_gethertzcontact(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_gethertzcontact: not available for icp=%d',icp));
      return
   end

   lenarr = 10;
   p_rvalues = libpointer('doublePtr',zeros(lenarr,1));

   calllib(libname,'cntc_gethertzcontact', ire, icp, lenarr, p_rvalues);

   rvalues = p_rvalues.value;

end % cntc_gethertzcontact

%------------------------------------------------------------------------------------------------------------

