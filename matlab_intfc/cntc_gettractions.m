
%------------------------------------------------------------------------------------------------------------
% function [ pn, px, py ] = cntc_gettractions(ire, icp)
%
% return the tractions for all elements in the potential contact area for a contact problem
%            note the order of the arguments, with pn (z-direction) occurring before px,py.
%
%  mx, my                            - number of elements in potential contact area
%  pn(my,mx), px(my,mx), py(my,mx)   - surface tractions for all elements of contact area, [force/area]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ pn, px, py ] = cntc_gettractions(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_gettractions: not available for icp=%d',icp));
      return
   end

   [ mx, my ] = cntc_getnumelements(ire, icp);

   if (mx*my <= 0)
      pn = [ 0 ];
      px = [ 0 ];
      py = [ 0 ];
   else
      lenarr = mx * my;
      p_pn = libpointer('doublePtr',zeros(lenarr,1));
      p_px = libpointer('doublePtr',zeros(lenarr,1));
      p_py = libpointer('doublePtr',zeros(lenarr,1));

      calllib(libname,'cntc_gettractions', ire, icp, lenarr, p_pn, p_px, p_py);

      pn = reshape(p_pn.value, mx, my)';
      px = reshape(p_px.value, mx, my)';
      py = reshape(p_py.value, mx, my)';
   end

end % cntc_gettractions

%------------------------------------------------------------------------------------------------------------

