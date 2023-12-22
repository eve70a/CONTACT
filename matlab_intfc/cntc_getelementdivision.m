
%------------------------------------------------------------------------------------------------------------
% function [ eldiv ] = cntc_getelementdivision(ire, icp)
%
% return flags for all elements in the potential contact area for a contact problem
% indicating whether the element is in Exerior (0), Adhesion (1), Slip (2) or Plasticity (3).
%
%  mx, my       - number of elements in potential contact area
%  eldiv(my,mx) - element division of contact area, 0=E, 1=H, 2=S, 3=P.
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ eldiv ] = cntc_getelementdivision(ire, icp)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_getelementdivision: not available for icp=%d',icp));
      return
   end

   [ mx, my ] = cntc_getnumelements(ire, icp);

   if (mx*my <= 0)
      eldiv = [ 0 ];
   else
      lenarr = mx * my;
      p_eldiv = libpointer('int32Ptr',zeros(lenarr,1));

      calllib(libname,'cntc_getelementdivision', ire, icp, lenarr, p_eldiv);

      % convert to double for consistency with loadcase.m
      % note that for int32-arrays, NaN == 0
      eldiv = double(reshape(p_eldiv.value, mx, my)');
   end

end % cntc_getelementdivision

%------------------------------------------------------------------------------------------------------------

