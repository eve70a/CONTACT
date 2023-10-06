
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setsolverflags(ire, icp, gdigit, iparam, rparam)
%
% set parameters for the iterative solution algorithms
%  gdigit            - G-digit, which solvers to use
%  iparams, rparam   - depending on method that is used
%
%  0: use default solvers      - iparam = [maxgs, maxin, maxnr, maxout],         rparam = [eps]
%  1: keep parameters          - iparam = [ ],                                   rparam = [ ]
%  2: always use ConvexGS      - iparam = [maxgs, maxin, maxnr, maxout, inislp], rparam = [eps, omegah, omegas, omgslp]
%  3: use SteadyGS if possible - iparam = [maxgs, maxin, maxnr, maxout, inislp], rparam = [eps, omegah, omegas, omgslp]
%  4: use default solvers      - iparam = [maxgs, maxin, maxnr, maxout, inislp], rparam = [eps, omgslp]
%  5: GDsteady if possible     - iparam = [maxgs, maxin, maxnr, maxout]          rparam = [eps, pow_s, omg_s]
%  6: flags sensitivities      - iparam = [mxsens],                              rparam = [epsens]
%
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 5: m=*, wtd    - default icp=-1

function [ ] = cntc_setsolverflags(ire, icp, gdigit, iparam, rparam)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = -1; % default: W/R contact, all patches
   end
   if (nargin<3 | isempty(gdigit))
      disp('ERROR in cntc_setsolverflags: gdigit is mandatory.');
      return
   end
   if (gdigit==1)
      nints = 0; iparam = 0; nreals = 0; rparam = 0;
   end
   if (nargin<5 | isempty(iparam) | isempty(rparam))
      disp('ERROR in cntc_setsolverflags: iparam, rparam are mandatory.');
      return
   end
   if (gdigit~=1)
      nints = length(iparam); nreals = length(rparam);
   end

   calllib(libname,'cntc_setsolverflags', ire, icp, gdigit, nints, iparam, nreals, rparam);

end % cntc_setsolverflags

%------------------------------------------------------------------------------------------------------------

