
function [ h ] = plot_axes3(O, scl, rot, axnams, col, mrksiz, ...
                           fontsiz, txtcol, th_txt, fac_txt, ...
                           no_head, left_hnd)

% function [ h ] = plot_axes3(O, scl, rot, axnams, col, mrksiz, ...
%                            fontsiz, txtcol, th_txt, fac_txt, ...
%                            no_head, left_hnd)
%
% 3d plot of coordinate axes Oxyz rotated by [roll,yaw,pitch].
%
% O        = position of origin on the current axes
% scl      = length of vectors on the current axes
% rot      = Euler angles [r,y,p] in [deg] or 3x3 rotation matrix
% axnams   = labels to be put on axes x1,x2,x3
% col      = color-spec
% mrksiz   = marker size (dot in circle for origin, default 8)
% fontsiz  = font size
% txtcol   = color-spec for text labels
% th_txt   = rotation angle [deg] for text labels
% fac_txt  = scale factor(s) for label positioning, can be 3x1 or 3x3
% no_head  = option to suppress arrow heads
% left_hnd = option to use left-handed sense of rotation

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
if (nargin<3 | isempty(rot))
   rot = eye(3);
end
if (nargin<4 | isempty(axnams))
   axnams = 'xyz';
end
if (size(axnams,1)==1)
   axnams = axnams';
end
if (nargin<5 | isempty(col))
   col = 'b';
end
if (nargin<6 | isempty(mrksiz))
   mrksiz = 8;
end
if (nargin<7 | isempty(fontsiz))
   fontsiz = 16;
end
if (nargin<8 | isempty(txtcol))
   txtcol = 'k';
end
if (nargin<9 | isempty(th_txt))
   th_txt = 0;
end
if (nargin<10 | isempty(fac_txt))
   fac_txt = 1;
end
if (nargin<11 | isempty(no_head))
   no_head = 0;
end
if (nargin<12 | isempty(left_hnd))
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

% convert color-spec to [r g b]

if (isnumeric(col) & length(col)==1)
   col = matlab_color(col);
end
if (isnumeric(txtcol) & length(txtcol)==1)
   txtcol = matlab_color(txtcol);
end

% set size of arrow head

if (no_head)
   scl_head = 0.001;
else
   scl_head = 1;
end

% expand fac_txt to 3x2 matrix

if (size(fac_txt,2)==3), fac_txt = fac_txt'; end
if (size(fac_txt,1)==1), fac_txt = ones(3,1)*fac_txt; end
if (size(fac_txt,2)==1), fac_txt = fac_txt*ones(1,3); end

% text labels: pos(1) * ivec + pos(2) * jvec

pos_txt = scl * fac_txt .* [1.00, 0.25, 0;   0.25, 1.00, 0;   0.20, 0.20, 1.00];

% compute unit vectors ivec, jvec, kvec in Matlab x,y,z-coordinates

if (any(size(rot)~=[3 3]))
   roll  = rot(1);
   yaw   = rot(2);
   pitch = rot(3);
   rot   = rotx(roll) * rotz(yaw) * roty(pitch);
end
ivec = rot(:,1);
jvec = rot(:,2);
kvec = rot(:,3);

% initialize output-array: handles [dot, line,line,txt, line,line,txt, line,line,txt]

h = zeros(10,1);

% put dot at origin O

h(1) = plot3(O(1), O(2), O(3), '.','color',col, 'markersize',mrksiz);
hold on;

intrp = 'tex';
if (strfind(reshape(axnams,1,prod(size(axnams))), '$'))
   intrp = 'latex';
end

% plot first axis ivec

h(2:3) = plot_arrow(O, scl*ivec, col, scl_head, [], jvec);

l1 = O + pos_txt(1,1)*ivec + pos_txt(1,2)*jvec + pos_txt(1,3)*kvec;
h(4) = text(l1(1), l1(2), l1(3), deblank(axnams(1,:)), 'rotation', th_txt, ...
                   'horizontalalignment','center', 'interpreter', intrp, ...
                   'fontsize',fontsiz, 'color', txtcol);

% plot second axis jvec

h(5:6) = plot_arrow(O, scl*jvec, col, scl_head, [], ivec);

l2 = O + pos_txt(2,1)*ivec + pos_txt(2,2)*jvec + pos_txt(2,3)*kvec;
h(7) = text(l2(1), l2(2), l2(3), deblank(axnams(2,:)), 'rotation', th_txt, ...
                   'horizontalalignment','center', 'interpreter', intrp, ...
                   'fontsize',fontsiz, 'color', txtcol);

% plot third axis kvec

h(8:9) = plot_arrow(O, scl*kvec, col, scl_head, [], jvec);

l3 = O + pos_txt(3,1)*ivec + pos_txt(3,2)*jvec + pos_txt(3,3)*kvec;
h(10) = text(l3(1), l3(2), l3(3), deblank(axnams(3,:)), 'rotation',th_txt, ...
                   'horizontalalignment','center', 'interpreter', intrp, ...
                   'fontsize',fontsiz, 'color', txtcol);

end % function plot_axes3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ rot ] = rotx(roll_deg)

  sn = sin(roll_deg*pi/180); 
  cs = cos(roll_deg*pi/180);
  rot = [ 1,  0,   0;
          0, cs, -sn;
          0, sn,  cs];

end % function rotx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ rot ] = roty(pitch_deg)

  sn = sin(pitch_deg*pi/180); 
  cs = cos(pitch_deg*pi/180);
  rot = [ cs,  0,  sn;
           0,  1,   0;
         -sn,  0,  cs];

end % function roty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ rot ] = rotz(yaw_deg)

  sn = sin(yaw_deg*pi/180); 
  cs = cos(yaw_deg*pi/180);
  rot = [ cs, -sn, 0;
          sn,  cs, 0;
           0,   0, 1];

end % function rotz

