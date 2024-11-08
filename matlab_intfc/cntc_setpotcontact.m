
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setpotcontact(ire, icp, ipotcn, params)
%
% set the parameters of the potential contact area for a contact problem
%  ipotcn   - type of specification for the potential contact area
%  params   - depending on method that is used
%
%  for w/r contact, icp = -1:
%    0: w/r contact, fixed grid sizes,  params = [ dx, ds, n.a. ]
%   -1: w/r contact, fixed grid sizes,  params = [ dx, ds, a_sep, d_sep, d_comb, [d_turn, g_miss] ]
%
%  for generic contact, icp > 0:
%    1: lower-left + grid sizes,        params = [ mx, my, xl , yl , dx , dy  ]
%    2: lower-left + upper right,       params = [ mx, my, xl , yl , xh , yh  ]
%    3: 1st center + grid sizes,        params = [ mx, my, xc1, yc1, dx , dy  ]
%    4: 1st center + last center,       params = [ mx, my, xc1, yc1, xcm, ycm ]
%
%  dimensions: mx, my [-],   a_sep [angle],
%              dx, ds, dy, d_sep, d_comb, d_turn, g_miss, xl, yl, xh, yh, xc1, yc1, xcm, ycm [length]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 5: m=*, wtd    - default icp=-1

function [ ] = cntc_setpotcontact(ire, icp, ipotcn, params)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = -1;
   end
   if (nargin<4 | isempty(ipotcn) | isempty(params))
      disp('ERROR in cntc_setpotcontact: ipotcn, params are mandatory.');
      return
   end

   calllib(libname,'cntc_setpotcontact', ire, icp, ipotcn, length(params), params);

end % cntc_setpotcontact

%------------------------------------------------------------------------------------------------------------

