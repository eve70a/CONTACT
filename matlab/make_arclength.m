
function [ s ] = make_arclength( y, z )

% function [ s ] = make_arclength( y, z )
%
% compute arc-length parameterization { s_i } of curve { (y_i, z_i) }:
% cumulative sum of distances ds between successive points

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

dy = diff(y);
dz = diff(z);
ds = sqrt(dy.^2 + dz.^2);
s  = [0; cumsum(ds) ];

% interpolate [y, s] to y=0, shift s := s - s0

is_monotonic = (all(diff(y)>0) | all(diff(y)<0));
ix_zeros     = find(abs(y)<1e-20);
if (isempty(ix_zeros))
   ix_zeros  = find(y(1:end-1).*y(2:end)<0);
end

if (is_monotonic & min(y)<0 & max(y)>0)
   sref = interp1(y, s, 0);
elseif (length(ix_zeros)==1)
   if (ix_zeros<length(y))
      ix   = ix_zeros + [ 0,1];
   else
      ix   = ix_zeros + [-1,0];
   end
   sref = interp1(y(ix), s(ix), 0);
else
   sref = 0;
end
s  = s - sref;

end % function make_arclength

