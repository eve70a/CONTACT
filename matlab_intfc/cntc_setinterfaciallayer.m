
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setinterfaciallayer(ire, icp, imeth, params)
%
% set parameters for the interfacial layer for a contact problem.   
%                                                       *Obsolete, replaced by cntc_setmaterialparameters.
%  imeth     - type of interfacial layer used
%  params    - depending on method that is used
%
%    0: clean interface, no layer     params = [ ]
%  2,3: Modified FASTSIM algorithm    params = [ k0_mf, alfamf, betamf ]
%    4: elasto-plastic layer          params = [ G3, laythk, tau_c0, k_tau ]
%    
%  dimensions: k0_mf, alfamf, betamf: [-],
%              G3, tau_c0: [force/area],  laythk: [length],  k_tau: [force/area/length]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 5: m=*, wtd    - default icp=-1

function [ ] = cntc_setinterfaciallayer(ire, icp, imeth, params)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = -1; % default: W/R contact, all patches
   end
   if (nargin<3 | isempty(imeth))
      disp('ERROR in cntc_setinterfaciallayer: imeth is mandatory.');
      return
   end
   if (imeth==0)
      params = [ 0 ];
   end
   if (nargin<4 | isempty(params))
      disp('ERROR in cntc_setinterfaciallayer: invalid params provided.');
      return
   end

   calllib(libname,'cntc_setinterfaciallayer', ire, icp, imeth, length(params), params);

end % cntc_setinterfaciallayer

%------------------------------------------------------------------------------------------------------------

