
function [ h ] = plot_cm(O, r, veclen, theta, col, no_vec, no_head)

% function [ h ] = plot_cm(O, r, veclen, theta, col, no_vec, no_head)
%
% plot a CM marker
%  O       - position of CM wrt Matlab axes
%  r       - radius of symbol
%  veclen  - length of arrows - starting at r
%  theta   - rotation wrt coordinate system [deg]
%  col     - color
%  no_vec  - disable vectors
%  no_head - disable vector heads

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

if (nargin<1 | isempty(O))
   O = [0; 0];
end
if (size(O,2)>size(O,1))
   O = O';
end
if (nargin<2 | isempty(r))
   r = 1;
end
if (nargin<3 | isempty(veclen))
   veclen = 2*r;
end
if (nargin<4 | isempty(theta))
   theta = 0;
end
if (nargin<5 | isempty(col))
   col = 'b';
end
if (nargin<6 | isempty(no_vec))
   no_vec = 0;
end
if (nargin<7 | isempty(no_head))
   no_head = 0;
end

if (isnumeric(col) & length(col)==1)
   col = matlab_color(col);
end
if (no_head)
   scl_head = 0.001;
else
   scl_head = 1;
end

h = zeros(7,1);

hold on;

% put circle at origin O

th = [0:5:360]*pi/180;
sn = sin(th); cs=cos(th);
h(1) = plot(O(1)+r*cs, O(2)+r*sn,'-','color',col);

% put filled quarter circles in 2nd & 4th quadrants

th = (theta + [90:5:180])*pi/180;
sn = sin(th); cs=cos(th);
h(2) = fill(O(1)+[0,r*cs,0], O(2)+[0,r*sn,0], col, 'edgecolor','none');

th = (theta + [270:5:360])*pi/180;
sn = sin(th); cs=cos(th);
h(3) = fill(O(1)+[0,r*cs,0], O(2)+[0,r*sn,0], col, 'edgecolor','none');


if (~no_vec)
   % plot x-axis

   n = [ cos(theta*pi/180); sin(theta*pi/180)];
   O1 = O + r*n;
   h(4:5) = plot_arrow(O1, veclen*n, col, scl_head);

   % plot y-axis

   t = [ -n(2) ; n(1) ];
   O2 = O + r*t;
   h(6:7) = plot_arrow(O2, veclen*t, col, scl_head);
end

