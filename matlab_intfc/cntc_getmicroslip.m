
%------------------------------------------------------------------------------------------------------------
% function [ sx, sy ] = cntc_getmicroslip(ire, icp)
%
% return the relative micro-slip velocity for all elements in the potential contact area for
%            a contact problem
%
%  mx, my               - number of elements in potential contact area
%  sx(my,mx), sy(my,mx) - in rolling, T=2,3: relative micro-slip velocity [-]
%                         in shifts,  T=1:   shift distance [length]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ sx, sy ] = cntc_getmicroslip(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_getmicroslip: not available for icp=%d',icp));
      return
   end

   [ mx, my ] = cntc_getnumelements(ire, icp);

   if (mx*my <= 0)
      sx = [ 0 ];
      sy = [ 0 ];
   else
      lenarr = mx * my;
      p_sx = libpointer('doublePtr',zeros(lenarr,1));
      p_sy = libpointer('doublePtr',zeros(lenarr,1));

      calllib(libname,'cntc_getmicroslip', ire, icp, lenarr, p_sx, p_sy);

      sx = reshape(p_sx.value, mx, my)';
      sy = reshape(p_sy.value, mx, my)';
   end

end % cntc_getmicroslip

%------------------------------------------------------------------------------------------------------------

