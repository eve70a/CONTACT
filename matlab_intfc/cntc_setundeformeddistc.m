
%------------------------------------------------------------------------------------------------------------
% function [ ] = cntc_setundeformeddistc(ire, icp, ibase, prmudf)
%
% set the undeformed distance function through a formula or by element-wise specification
%
%  ibase          - type of undeformed distance specification
%  nparam         - number of parameters provided
%  prmudf(nparam) - parameters of undef.dist, depending on method that is used
%
%    1: quadratic function            6    params = [b1, b2, b3, b4, b5, b6]
%    2: circular-x, piecewise-lin-y   5+nn params = [nn, xm, rm, y1, dy1], [b(k), k=1..nn]
%    3: quadratic plus two sines      8    params = [b1, b2, b3, b4, b5, b6, b7, b8]
%    9: elementwise specification     npot params = [h(i), i=1..npot] - undeformed distance per elem. [length]
%                                                   note: positive values == separation between profiles
%
% when ibase=9, prmudf may be of size (my,mx) as well, with nparam=mx*my.
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 1: m=3, cp     - require icp>0, default 1

function [ ] = cntc_setundeformeddistc(ire, icp, ibase, prmudf)
   global libname;

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_setundeformeddistc: not available for icp=%d',icp));
      return
   end
   if (nargin<4 | isempty(ibase) | isempty(prmudf))
      disp('ERROR in cntc_setundeformeddistc: ibase, prmudf are mandatory.');
      return
   end

   % check size of prmudf-array

   nparam = prod(size(prmudf));
   if (ibase==9)
      [ mx, my ] = cntc_getnumelements(ire, icp);
      if (nparam~=mx*my)
         disp(sprintf('ERROR in cntc_setundeformeddistc: expecting nparam = mx*my = %d',mx*my));
         return
      end
      if (min(size(prmudf))>1 & any(size(prmudf)~=[my mx]))
         disp('ERROR in cntc_setundeformeddistc:');
         disp(sprintf('      expecting array prmudf with size (%d,%d), got size (%d,%d)',my,mx,size(prmudf)))
         return
      end
   end

   % if prmudf is 2-dimensional, reshape to 1-d array

   if (ibase==9 & min(size(prmudf))>1)
      prmudf = reshape(prmudf', nparam, 1);
   end

   calllib(libname,'cntc_setundeformeddistc', ire, icp, ibase, nparam, prmudf);

end % cntc_setundeformeddistc

%------------------------------------------------------------------------------------------------------------

