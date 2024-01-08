
%------------------------------------------------------------------------------------------------------------
% function [ values ] = cntc_getparameters(ire, icp)
%
% used for retrieving various parameters from a contact problem needed for cntc_getcpresults
%
%  lenarr         - length of values array
%  values(lenarr) - values of the parameters obtained from CONTACT
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 5: m=*, wtd    - default icp=-1

function [ values ] = cntc_getparameters(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = -1; % default: W/R contact, all patches
   end

   lenarr = 10;
   p_values = libpointer('doublePtr',zeros(lenarr,1));

   calllib(libname,'cntc_getparameters', ire, icp, lenarr, p_values);

   tmp    = p_values.value;
   values = struct();
   values.veloc  = tmp(1);
   values.chi    = tmp(2);
   values.dq     = tmp(3);
   values.spinxo = tmp(4);
   values.spinyo = tmp(5);
   values.tau_c0 = tmp(6);

end % cntc_getparameters

%------------------------------------------------------------------------------------------------------------

