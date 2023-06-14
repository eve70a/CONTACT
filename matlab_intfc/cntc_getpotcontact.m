
%------------------------------------------------------------------------------------------------------------
% function [ mx, my, xc1, yc1, dx, dy ] = cntc_getpotcontact(ire, icp)
%
% get the parameters of the potential contact area for a contact problem
%    3: first center + grid sizes,        params = [ mx, my, xc1, yc1, dx, dy ]
%
%  mx, my        - number of elements in x- and y-directions [-]
%  xc1, yc1      - position of first element center [length]
%  dx, dy        - grid discretization step sizes [length]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ mx, my, xc1, yc1, dx, dy ] = cntc_getpotcontact(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_getpotcontact: not available for icp=%d',icp));
      return
   end

   lenarr = 6;
   p_values = libpointer('doublePtr',zeros(lenarr,1));

   calllib(libname,'cntc_getpotcontact', ire, icp, lenarr, p_values);

   values = p_values.value;

   mx  = round(values(1));
   my  = round(values(2));
   xc1 = values(3);
   yc1 = values(4);
   dx  = values(5);
   dy  = values(6);

end % cntc_getpotcontact

%------------------------------------------------------------------------------------------------------------
