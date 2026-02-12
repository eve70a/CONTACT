
%------------------------------------------------------------------------------------------------------------
% function [ fld ] = cntc_getfielddata(ire, icp, ifld)
%
% return the values of field 'ifld' for all elements in the potential contact area for a contact problem
%
%  mx, my        - number of elements in potential contact area
%  fld(my,mx)    - output value of field 'ifld' for all elements of contact area
%
%                  ifld  =  CNTC_fld_h      =    1 ! for retrieving array h         [length]
%                           CNTC_fld_mu     =    2 ! for retrieving array mu        [-]
%                           CNTC_fld_px     =    3 ! for retrieving array px        [force/area]
%                           CNTC_fld_py     =    4 ! for retrieving array py        [force/area]
%                           CNTC_fld_pn     =    5 ! for retrieving array pn        [force/area]
%                           CNTC_fld_ux     =    7 ! for retrieving array ux        [length]
%                           CNTC_fld_uy     =    8 ! for retrieving array uy        [length]
%                           CNTC_fld_un     =    9 ! for retrieving array un        [length]
%                           CNTC_fld_taucrt =   11 ! for retrieving array taucrt    [force/area]
%                           CNTC_fld_uplsx  =   12 ! for retrieving array uplsx     [length]
%                           CNTC_fld_uplsy  =   13 ! for retrieving array uplsy     [length]
%                           CNTC_fld_sx     =   15 ! for retrieving array sx        [-]
%                           CNTC_fld_sy     =   16 ! for retrieving array sy        [-]
%                           CNTC_fld_temp1  =   20 ! for retrieving array temp1     [C]
%                           CNTC_fld_temp2  =   21 ! for retrieving array temp2     [C]
%                           CNTC_fld_wx     =   22 ! for retrieving array wx        [-]
%                           CNTC_fld_wy     =   23 ! for retrieving array wy        [-]
%
%  fld(my)       - output value of field 'ifld' for all rows of contact area
%                           CNTC_fld_crefy  =   32 ! for retrieving array cref_y    [length]
%                           CNTC_fld_crefz  =   33 ! for retrieving array cref_z    [length]
%                           CNTC_fld_cangl  =   34 ! for retrieving array cref_incl [angle]

%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2026 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ fld ] = cntc_getfielddata(ire, icp, ifld)
   global libname;

   CNTC = cntc_getmagicnumbers();

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
   if (nargin<3 | isempty(ifld))
      ifld = CNTC.fld_pn;
   end

   [ mx, my ] = cntc_getnumelements(ire, icp);

   if (mx*my <= 0)
      fld = [0];
   else
      if (ifld>=CNTC.fld_crefy)
         lenarr =      my;      % array (my) per row of c.grid
      else
         lenarr = mx * my;      % array (my,mx) for entire c.grid
      end
      p_val = libpointer('doublePtr',zeros(lenarr,1));

      calllib(libname,'cntc_getfielddata', ire, icp, ifld, lenarr, p_val);

      if (lenarr==my)
         fld = p_val.value;
      else
         fld = reshape(p_val.value, mx, my)';
      end
   end

end % cntc_getfielddata

%------------------------------------------------------------------------------------------------------------

