
function [ h ] = plot_arrow( pos0, vec, col, scale, width, pln )
%
% function [ h ] = plot_arrow( pos0, vec, col, scale, width, [pln] )
%
% pos0  = position of tail of arrow (may be 3-vector for plot3)
% vec   = vector from tail to head of arrow (may be a 3-vector for plot3)
% col   = color-spec for arrow; default='b';
% scale = scale factor for arrow head; default=1: 25% of arrow size
% width = width factor for arrow head; default=1: about 30 deg.
% pln   = binormal direction used in 3d plot to define plane for arrow-head
%

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

if (nargin<3 | isempty(col))
   col = 'b';
end
if (nargin<4 | isempty(scale))
   scale = 1;
end
if (nargin<5 | isempty(width))
   width = 1;
end
if (nargin<6 | isempty(pln))
   pln = [ 0, 0, 1 ];
end

% if col is just a single value, it is interpreted as a color index

if (isnumeric(col) & length(col)==1)
   col = matlab_color(col);
end

% convert arguments into column vectors

if (size(pos0,2)>1), pos0=pos0'; end
if (size(vec ,2)>1), vec =vec'; end
if (size(pln ,2)>1), pln =pln'; end
pln = pln / norm(pln);

% use plot3 instead of plot when pos is a 3-vector

is_3d = (length(pos0)==3);

% adapt binormal pln if in direction of vec

if (is_3d)
   vp = vec'*pln / norm(vec);
   if (abs(vp)>0.99)
      if (abs(pln'*[0;0;1])>0.99)
         pln = [0;1;0];
      else
         pln = [0;0;1];
      end
   end
end

% set line to be plotted

lin  = [pos0, pos0+vec];

% determine arrow head

tng  = vec;
if (is_3d)
   nrm  = cross(tng, pln);
else
   asp  = get(gca,'dataaspectratio');
   nrm  = [-tng(2)*asp(1); tng(1)*asp(2)] / sqrt(asp(1)*asp(2));
end

head = [ pos0+vec-0.25*scale*tng+0.125*width*scale*nrm,  ...
         pos0+vec, ...
         pos0+vec-0.25*scale*tng-0.125*width*scale*nrm];

if (~is_3d)
   h(1) = plot(lin(1,:), lin(2,:), 'color',col);
   hold on;
   h(2) = plot(head(1,:), head(2,:), 'color',col);
else
   nl = size(lin,2); nh = size(head,2);
   h(1) = plot3(lin(1,:), lin(2,:), lin(3,:), 'color',col);
   hold on;
   h(2) = plot3(head(1,:), head(2,:), head(3,:), 'color',col);
end

%stel norm = [ 1, 0 ]
%     aspc = [ 2, 1 ]
% ==> transpose = [ 0, 1 ] = half zo lang.
