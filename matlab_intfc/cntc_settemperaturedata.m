
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_settemperaturedata(ire, icp, imeth, params)
%
% set parameters for the temperature calculation for a contact problem
%  imeth    - type of temperature model used (H-digit)
%  params   - depending on method that is used
%    0: no temperature calculation,     params = []
%    1: keep old parameters,            params = []
%    3: calculate temperature based on new parameters and steady rolling,
%       params = [bktemp1, heatcp1, lambda1, dens1, bktemp2, heatcp2, lambda2, dens2]
%
% dimensions:  bktemp: [C],  heatcp: [J/kg-C],  lambda: [W/length-C],  dens: [kg/length^3]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 5: m=*, wtd    - default icp=-1

function [ ] = cntc_settemperaturedata(ire, icp, imeth, params)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = -1; % default: W/R contact, all patches
   end
   if (nargin<3 | isempty(imeth))
      disp('ERROR in cntc_settemperaturedata: imeth is mandatory.');
      return
   end
   if (imeth<=1)
      params = [ 0 ];
   end
   if (nargin<4 | isempty(params))
      disp('ERROR in cntc_settemperaturedata: invalid params provided.');
      return
   end

   calllib(libname,'cntc_settemperaturedata', ire, icp, imeth, length(params), params);

end % cntc_settemperaturedata

%------------------------------------------------------------------------------------------------------------

