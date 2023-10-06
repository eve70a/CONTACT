
%------------------------------------------------------------------------------------------------------------
% function [ fn, tx, ty, mz ] = cntc_getcontactforces(ire, icp)
%
% return the total forces and torsional moment for a contact problem in contact local coordinates
%
%  fn            - total normal force [force]
%  tx, ty        - total tangential forces [force]
%  mz            - total torsional moment [force.length]
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ fn, tx, ty, mz ] = cntc_getcontactforces(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_getcontactforces: not available for icp=%d',icp));
      return
   end

   p_fn = libpointer('doublePtr',-1);
   p_tx = libpointer('doublePtr',-1);
   p_ty = libpointer('doublePtr',-1);
   p_mz = libpointer('doublePtr',-1);

   calllib(libname,'cntc_getcontactforces', ire, icp, p_fn, p_tx, p_ty, p_mz);

   fn = p_fn.value;
   tx = p_tx.value;
   ty = p_ty.value;
   mz = p_mz.value;

end % cntc_getcontactforces

%------------------------------------------------------------------------------------------------------------

