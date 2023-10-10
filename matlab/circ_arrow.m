
function [ h ] = circ_arrow( pos_c, r, th0, th1, col, scale, width, rot )
%
% function [ h ] = circ_arrow( pos_c, r, th0, th1, col, scale, width, [rot] )
%
% pos_c = position of arrow center of rotation, 2-vector or 3-vector
% r     = circle radius;  negative r: reversed rotation direction.
%         vectorial: [rx,ry], separate radius for x and y-direction
% th0   = angle of tail of arrow [deg]
% th1   = angle of head of arrow
% col   = color-spec for arrow; default='b';
% scale = scale factor for arrow head; default=1: 25% of arrow size
% width = width factor for arrow head; default=1: about 30 deg.
% rot   = for 3d plots: rotation matrix from [x,y,0] to desired orientation
%         [0,0,1; 1,0,0; 0,1,0] for arc in yz-plane
%

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

if (nargin<5 | isempty(col))
   col = 1;
end
if (nargin<6 | isempty(scale))
   scale = 1;
end
if (nargin<7 | isempty(width))
   width = 1;
end
if (nargin<8 | isempty(rot))
   rot   = eye(3);
end

% if r is a scalar, copy to make a 2-vector [rx,ry]
if (max(size(r))==1)
   r = [r,r];
end
rx=r(1); ry=r(2);

% if col is just a single value, it is interpreted as a color index
if (isnumeric(col) & length(col)==1)
   col = matlab_color(col);
end

if (size(pos_c,2)>1), pos_c=pos_c'; end % make column vectors

th0 = th0 * pi/180;             % convert to radians
th1 = th1 * pi/180;             % convert to radians

% use plot3 instead of plot when pos is a 3-vector

is_3d = (length(pos_c)==3);
if (is_3d)
   pos_c = rot' * pos_c;
else
   pos_c = [pos_c; 0];
end

hold on;

% compute circular arc with center pos_c, radius [rx,ry,0] and angles theta (th)
th   = th0 + [0:30]/30 * (th1-th0);
lin  = pos_c*ones(size(th)) + [rx*cos(th); ry*sin(th); 0*th];

pos1 = lin(:,end);
% plot(pos_c(1), pos_c(2), 'b*');

% compute nice theta-value for getting the arrow head direction
th_h = th1 - min(1,0.25*scale) * (th1-th0);

% compute point on the curve for the back of the arrow head
pos_h = pos_c + [rx*cos(th_h); ry*sin(th_h); 0*th_h];

% compute tangent vector at th_h == orientation of head;
% make size of arrowhead proportional to length of curve (th1-th0)
tng  = [-ry*sin(th_h); rx*cos(th_h); 0] * (th1-th0);

% compute the normal direction
nrm  = [-tng(2); tng(1); 0];

% reverse direction when needed
if (th0>th1), nrm=-nrm; end

% compute arrow head points, on the inside closer to the curve
head = [ pos_h+0.10*width*scale*nrm,  ...
         pos1, ...
         pos_h-0.13*width*scale*nrm];

if (is_3d)
   lin  = rot * lin;
   head = rot * head;
   h(1) = plot3(lin(1,:), lin(2,:), lin(3,:), 'color',col);
   h(2) = plot3(head(1,:), head(2,:), head(3,:), 'color',col);
else
   h(1) = plot(lin(1,:), lin(2,:), 'color',col);
   h(2) = plot(head(1,:), head(2,:), 'color',col);
end
