
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setfrictionmethod(ire, icp, imeth, params)
%
% set parameters for the friction law for a contact problem
%  imeth          - type of friction law used: 10 * V-digit + 1 * L-digit
%  params         - depending on method that is used
%
%  L = 0: Coulomb friction,             lparam = [fstat, fkin]
%      1: keep previous parameters,     lparam = [ ]
%      2: linear falling friction,      lparam = [fkin, flin1, sabsh1, flin2, sabsh2, memdst, mem_s0]
%      3: rational falling friction,    lparam = [fkin, frat1, sabsh1, frat2, sabsh2, memdst, mem_s0]
%      4: exponential falling friction, lparam = [fkin, fexp1, sabsh1, fexp2, sabsh2, memdst, mem_s0]
%      5: exponential falling friction, lparam = [fstat, polach_a, polach_b]
%      6: temperature dep. friction,    lparam = [fref, tref, dfheat, dtheat, memdst, mem_s0]
%
%  When V = 0, params = lparam(1:n), with n=2 for L=0, n=7 for L=2--4, n=3 for L=5, n=6 for L=6
%  When V = 1, params = [ nvf, ...
%                         alphvf(1),   lparam( 1 ,1:n) ],  ...
%                                   ...
%                         alphvf(nvf), lparam(nvf,1:n) ]
%  When V = 2, same as for V = 1, using svf instead of alphvf
%    
%  dimensions:  alphvf: [angle],  fstat, fkin, flin1,2, frat1,2, fexp1,2, polach_a, fref, dfheat: [-]
%               sabsh1,2, mem_s0: [veloc],  svf, memdst: [length],  polach_b [1/veloc],  tref, dtheat [C]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 5: m=*, wtd    - default icp=-1

function [ ] = cntc_setfrictionmethod(ire, icp, imeth, params)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = -1; % default: W/R contact, all patches
   end
   if (nargin<3 | isempty(imeth))
      disp('ERROR in cntc_setfrictionmethod: imeth is mandatory.');
      return
   end
   if (nargin<4)
      params = [];
   end

   calllib(libname,'cntc_setfrictionmethod', ire, icp, imeth, length(params), params);

end % cntc_setfrictionmethod

%------------------------------------------------------------------------------------------------------------

