
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setmaterialproperties(ire, icp, g1, nu1, g2, nu2)
%
% set the material properties for a contact problem.    Obsolete, replaced by cntc_setmaterialparameters.
%
%  g1, g2         - modulus of rigidity for body 1, 2 [force/area]
%  nu1, nu2       - Poisson's ratio for body 1, 2 [-]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 5: m=*, wtd    - default icp=-1

function [ ] = cntc_setmaterialproperties(ire, icp, g1, nu1, g2, nu2)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = -1; % default: W/R contact, all patches
   end
   if (nargin<6 | isempty(g1) | isempty(nu1) | isempty(g2) | isempty(nu2))
      disp('ERROR in cntc_setmaterialproperties: g1, nu1, g2, nu2 are mandatory.');
      return
   end

   calllib(libname,'cntc_setmaterialproperties', ire, icp, g1, nu1, g2, nu2);

end % cntc_setmaterialproperties

%------------------------------------------------------------------------------------------------------------

