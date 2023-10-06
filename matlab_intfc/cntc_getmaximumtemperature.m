
%------------------------------------------------------------------------------------------------------------
% function [ t1max, t2max ] = cntc_getmaximumtemperature(ire, icp)
%
% return the maximum contact temperature in a contact problem
%
%  t1max, t2max  - maximum surface temperatures in bodies 1 and 2 in contact patch [C]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ t1max, t2max ] = cntc_getmaximumtemperature(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_getmaximumtemperature: not available for icp=%d',icp));
      return
   end

   p_t1max = libpointer('doublePtr',-1);
   p_t2max = libpointer('doublePtr',-1);

   calllib(libname,'cntc_getmaximumtemperature', ire, icp, p_t1max, p_t2max);

   t1max = p_t1max.value;
   t2max = p_t2max.value;

end % cntc_getmaximumtemperature

%------------------------------------------------------------------------------------------------------------

