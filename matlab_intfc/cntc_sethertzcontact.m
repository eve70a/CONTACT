
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_sethertzcontact(ire, icp, ipotcn, params)
%
% set the parameters for a Hertzian contact problem
%  ipotcn    - type of specification of the Hertzian geometry
%  params    - depending on method that is used
%
%   -6: SDEC approach, union of two half ellipses                  params = [ mx, my, aa , bneg, bpos, scale ]
%   -5: Hertzian rectangular contact, half-sizes prescribed,       params = [ mx, my, aa , bb , scale ]
%   -4: Hertzian rectangular contact, curv+half width prescribed,  params = [ mx, my, a1 , bb , scale ]
%   -3: Hertzian elliptical contact, semi-axes prescribed,         params = [ mx, my, aa , bb , scale ]
%   -2: Hertzian elliptical contact, ellipticity prescribed,       params = [ mx, my, a1 , aob, scale ]
%   -1: Hertzian elliptical contact, curvatures prescribed,        params = [ mx, my, a1 , b1 , scale ]
%
% dimensions:  mx, my, aob, scale [-],    a1, b1: [1/length],    aa, bneg, bpos, bb: [length]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 1: m=3, cp     - require icp>0, default 1

function [ ] = cntc_sethertzcontact(ire, icp, ipotcn, params)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_sethertzcontact: not available for icp=%d',icp));
      return
   end
   if (nargin<4 | isempty(ipotcn) | isempty(params))
      disp('ERROR in cntc_sethertzcontact: ipotcn, params are mandatory.');
      return
   end

   calllib(libname,'cntc_sethertzcontact', ire, icp, ipotcn, length(params), params);

end % cntc_sethertzcontact

%------------------------------------------------------------------------------------------------------------

