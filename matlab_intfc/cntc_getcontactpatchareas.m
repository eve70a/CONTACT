
%------------------------------------------------------------------------------------------------------------
% function [ carea, harea, sarea, parea ] = cntc_getcontactpatchareas(ire, icp)
%
% return the area of contact for a contact problem
%  carea   - area of contact patch [area]
%  harea   - area of adhesion area [area]
%  sarea   - area of slip area [area]
%  parea   - area of plasticity area [area]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ carea, harea, sarea, parea ] = cntc_getcontactpatchareas(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_getcontactpatchareas: not available for icp=%d',icp));
      return
   end

   p_carea  = libpointer('doublePtr',-1);
   p_harea  = libpointer('doublePtr',-1);
   p_sarea  = libpointer('doublePtr',-1);

   calllib(libname,'cntc_getcontactpatchareas', ire, icp, p_carea, p_harea, p_sarea);

   carea = p_carea.value;
   harea = p_harea.value;
   sarea = p_sarea.value;
   parea = carea - harea - sarea;

end % cntc_getcontactpatchareas

%------------------------------------------------------------------------------------------------------------

