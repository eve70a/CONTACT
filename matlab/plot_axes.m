
function [ h ] = plot_axes(O, scl, theta, zpos, axnams, col, mrksiz, ...
                           siz_o, fontsiz, txtcol, th_txt, fac_txt, ...
                           no_orie, no_head, left_hnd)

% function [ h ] = plot_axes(O, scl, theta, zpos, axnams, col, mrksiz, ...
%                            siz_o, fontsiz, txtcol, th_txt, fac_txt, ...
%                            no_orie, no_head, left_hnd)
%
% plot coordinate axes Oxyz when looked upon from either +inf or -inf
% at one coordinate axis.
%
% O        = position of origin on the current axes
% scl      = length of vectors on the current axes
% theta    = angle from Matlab x-axis to first coordinate axis [deg]
% zpos     = height of viewpoint on 3rd axis
% axnams   = labels to be put on axes x1,x2,x3
%                axnams = 'xyz', zpos <0: plot Oxy viewed from z = -inf
%                axnams = 'xyz', zpos>=0: plot Oxy viewed from z =  inf
%                axnams = 'yzx', zpos <0: plot Oyz viewed from x = -inf
%                axnams = 'zxy', zpos <0: plot Oxz viewed from y = -inf, theta=-90: z pointing upwards
%                axnams = 'zyx', zpos>=0: plot Oxz viewed from y =  inf, theta=-90: z pointing downwards
% col      = color-spec
% mrksiz   = marker size (dot in circle for origin, default 8)
% siz_o    = relative size of circle for origin (default 0.12)
% fontsiz  = font size
% txtcol   = color-spec for text labels
% th_txt   = rotation angle [deg] for text labels (default theta)
% fac_txt  = scale factor(s) for label positioning, can be 3x1 or 3x2
% no_orie  = option to suppress sense of rotation
% no_head  = option to suppress arrow heads
% left_hnd = option to use left-handed sense of rotation
%
% output h(1) = marker; h(2:4) = lines center; h(5:8) = arrows; h(9:12) = sense of rotation; h(13:15) = labels

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

if (nargin<1 | isempty(O))
   O = [0; 0];
end
if (size(O,2)>size(O,1))
   O = O';
end
if (nargin<2 | isempty(scl))
   scl = 1;
end
if (nargin<3 | isempty(theta))
   theta = 0;
end
if (nargin<4 | isempty(zpos))
   zpos = 1;
end
if (nargin<5 | isempty(axnams))
   axnams = 'xyz';
end
if (size(axnams,1)==1)
   axnams = axnams';
end
if (nargin<6 | isempty(col))
   col = 'b';
end
if (nargin<7 | isempty(mrksiz))
   mrksiz = 8;
end
if (nargin<8 | isempty(siz_o))
   siz_o   = 0.12;
end
if (nargin<9 | isempty(fontsiz))
   fontsiz = 16;
end
if (nargin<10 | isempty(txtcol))
   txtcol = 'k';
end
if (nargin<11 | isempty(th_txt))
   th_txt = theta;
end
if (nargin<12 | isempty(fac_txt))
   fac_txt = 1;
end
if (nargin<13 | isempty(no_orie))
   no_orie = 0;
end
if (nargin<14 | isempty(no_head))
   no_head = 0;
end
if (nargin<15 | isempty(left_hnd))
   left_hnd = 0;
end

% associate ivec with matlab x, jvec with matlab y, kvec out of the screen

% use left-handed convention if Matlab ydir is reversed

gca_ydir = get(gca,'ydir');
if (strcmp(gca_ydir,'reverse'))
   fac_y = -1;
else
   fac_y =  1;
end

% reverse sense of rotation when left_hnd is used

if (left_hnd)
   fac_y = -fac_y;
end

% reverse angles if z-position < 0

if (zpos >= 0)
   fac_ang =  1;
else
   fac_ang = -1;
end
theta  = fac_y * fac_ang * theta;
th_txt =         fac_ang * th_txt;

if (isnumeric(col) & length(col)==1)
   col = matlab_color(col);
end
if (isnumeric(txtcol) & length(txtcol)==1)
   txtcol = matlab_color(txtcol);
end
if (no_head)
   scl_head = 0.001;
else
   scl_head = 1;
end

% expand fac_txt to 3x2 matrix
if (size(fac_txt,2)==3), fac_txt = fac_txt'; end
if (size(fac_txt,1)==1), fac_txt = ones(3,1)*fac_txt; end
if (size(fac_txt,2)==1), fac_txt = fac_txt*ones(1,2); end

% text labels: pos(1) * ivec + pos(2) * jvec
pos_txt = scl * fac_txt .* [1.00, 0.25;   0.25, 1.00;   -0.10, -0.20];

% compute unit vectors ivec, jvec in Matlab x,y-coordinates

ivec = [ cos(theta*pi/180); sin(theta*pi/180) ];
jvec = fac_y * fac_ang * [ -ivec(2); ivec(1) ];

% initialize output-array: handles

h = zeros(14,1);

hold on;

th=[0:10:360]*pi/180;
sn=sin(th); cs=cos(th);

% put circle with cross or dot at origin O

r = siz_o*scl;
h(4) = plot(O(1)+r*cs, O(2)+r*sn,'-','color',col);

if (zpos>=0)
   % looking down x3-axis: dot
   h(1) = plot(O(1)+0*scl, O(2)+0*scl, '.','color',col, 'markersize',mrksiz);
else
   % looking up x3-axis: cross
   ij = (ivec+jvec) / sqrt(2);
   h(2) = plot(O(1)+r*[-1,1]*ij(1), O(2)+r*[-1,1]*ij(2), 'color',col);
   ij = (ivec-jvec) / sqrt(2);
   h(3) = plot(O(1)+r*[-1,1]*ij(1), O(2)+r*[-1,1]*ij(2), 'color',col);
end

intrp = 'tex';
if (strfind(reshape(axnams,1,prod(size(axnams))), '$'))
   intrp = 'latex';
end

% plot first axis ivec

O1 = O + r*ivec;
h(5:6) = plot_arrow(O1, (scl-r)*ivec, col, scl_head);

l1 = O + pos_txt(1,1)*ivec + pos_txt(1,2)*jvec;
h(13) = text(l1(1), l1(2), deblank(axnams(1,:)), 'rotation', th_txt, ...
                   'horizontalalignment','center', 'interpreter', intrp, ...
                   'fontsize',fontsiz, 'color', txtcol);

% plot second axis jvec

O2 = O + r*jvec;
h(7:8) = plot_arrow(O2, (scl-r)*jvec, col, scl_head);

l2 = O + pos_txt(2,1)*ivec + pos_txt(2,2)*jvec;
h(14) = text(l2(1), l2(2), deblank(axnams(2,:)), 'rotation', th_txt, ...
                   'horizontalalignment','center', 'interpreter', intrp, ...
                   'fontsize',fontsiz, 'color', txtcol);

% plot label for third axis kvec

l3 = O + pos_txt(3,1)*ivec + pos_txt(3,2)*jvec;
h(15) = text(l3(1), l3(2), deblank(axnams(3,:)), 'rotation',th_txt, ...
                   'horizontalalignment','center', 'interpreter', intrp, ...
                   'fontsize',fontsiz, 'color', txtcol);

% plot optional sense of positive rotation

if (no_orie<=0)
   Orot = O + 0.3*scl*ivec - 0.3*scl*jvec;
   f = 0.05*scl*[-1,1];
   h( 9) = plot(Orot(1)+ivec(1)*f, Orot(2)+ivec(2)*f, 'color',col);
   h(10) = plot(Orot(1)+jvec(1)*f, Orot(2)+jvec(2)*f, 'color',col);
   if (fac_ang*fac_y>0)
      h(11:12) = circ_arrow(Orot, 0.15*scl, theta-120, theta+150, col, ...
                                  0.7*scl_head, 0.8);
   else
      h(11:12) = circ_arrow(Orot, 0.15*scl, theta+120, theta-150, col, ...
                                  0.7*scl_head, 0.8);
   end
end

