
%------------------------------------------------------------------------------------------------------------
% function [ un, ux, uy ] = cntc_getdisplacements(ire, icp)
%
% return the displ.differences for all elements in the potential contact area for a contact problem
%  mx, my                            - number of elements in potential contact area
%  un(my,mx), ux(my,mx), uy(my,mx)   - displacement difference [length]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ un, ux, uy ] = cntc_getdisplacements(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_getdisplacements: not available for icp=%d',icp));
      return
   end

   [ mx, my ] = cntc_getnumelements(ire, icp);
   lenarr = mx * my;
   p_un = libpointer('doublePtr',zeros(lenarr,1));
   p_ux = libpointer('doublePtr',zeros(lenarr,1));
   p_uy = libpointer('doublePtr',zeros(lenarr,1));

   calllib(libname,'cntc_getdisplacements', ire, icp, lenarr, p_un, p_ux, p_uy);

   un = reshape(p_un.value, mx, my)';
   ux = reshape(p_ux.value, mx, my)';
   uy = reshape(p_uy.value, mx, my)';

end % cntc_getdisplacements

%------------------------------------------------------------------------------------------------------------

