
print_fig = 1;
iver = 2;

% wheel position

rnom  = 460;
x_ws  =   0;
z_ws  = -rnom;

if (iver==1)
   th_ws = 0;
else
   th_ws = 50;
end

% wheel surface

thw = [0 : 360];
xw = rnom * sin(thw*pi/180);
zw = rnom * cos(thw*pi/180);

% wheel surface rotated/translated to track coordinates

Rot = [ cos(th_ws*pi/180), sin(th_ws*pi/180);
       -sin(th_ws*pi/180), cos(th_ws*pi/180)];
tmp = Rot * [xw; zw];
xw  = x_ws + tmp(1,:);
zw  = z_ws + tmp(2,:);

% start plotting

figure(1); clf; hold on;
axis equal;
axis([-550 550 -930 160]);
grid on;
set(gca,'ydir','reverse');
xlabel('x_{tr} [mm]');
ylabel('z_{tr} [mm]');

% rail plane

if (iver>=2)
   plot(1.1*rnom*[-1 -0.05 NaN 0.05 1], [0 0 0 0 0], 'color',matlab_color(1));

   % text(0, -100, 'O_{track}', 'horizontalalignment','right');
   % text(len, -20, 'y_{track}', 'verticalalignment', 'top');
   % text(-70, len-20, 'z_{track}', 'verticalalignment','top');

   plot([x_ws,xw(361-th_ws)], [z_ws+110,zw(361-th_ws)],':', 'color',matlab_color(1));
end

% wheel body coordinate system

if (iver<=2)
   l = plot_axes([x_ws,z_ws], 150, th_ws-90, 1, ['$z_{wcyl}$';'$x_{wcyl}$';'$O_{wcyl}$'], 'k', ...
                                [], [], [], [], th_ws, [1.2,1.2,4.5;2.0,1.2,2])
   % set(l([11:14]), 'linewidth', 1.2);
   delete(l([1,2,4,5,7,8,11:14]));
   if (iver==2), delete(l([6,9])); end
end

l = plot_cm([x_ws,z_ws], 18, 132, -th_ws, 'k');
set(l([1,4:end]), 'linewidth', 1.2);

% wheel outline

set(gca,'colororderindex',5);
plot(xw, zw)

% wheel point at theta=0

set(gca,'colororderindex',5);
plot(xw(1), zw(1), '.', 'markersize', 15)

set(gca,'colororderindex',5);
plot([x_ws,xw(1)], [z_ws,zw(1)],':');

text(xw(1)+1, zw(1)+5, '$\theta_{wc}=0$', 'verticalalignment','top', 'interpreter','latex');

p = [x_ws; z_ws] + 0.75*rnom * [sin(th_ws*pi/180); cos(th_ws*pi/180)];
text(p(1)+15, p(2), '$r_{wc}$', 'interpreter','latex', 'verticalalignment','bottom');

% plot track coordinate axes

if (iver>=2)
   len = 150;
   l = plot_axes([0,0], 150, -90, 1, ['$z_{tr}$';'$x_{tr}$';'        '], 'k', [], [], [], [], 0, 1.2);
   set(l([1,2,4,5,7,8,11:14]),'linewidth',1);
   % set(l, 'linewidth',1);
end

% - plot angle from x_m-axis (angle th_m) to -vertical (angle -90)

if (iver>=2)
   col = matlab_color(3);
   l = circ_arrow( [x_ws,z_ws], 0.5*rnom, 90, 90-th_ws, col, 0.7, 0.7 );
   set(l, 'linewidth', 1.5);

   p = [x_ws; z_ws] + 0.52*rnom * [sin(th_ws/2*pi/180); cos(th_ws/2*pi/180)];
   text(p(1), p(2), '$\theta_{ws}$', 'interpreter','latex','verticalalignment','top');
end

if (print_fig)
   axis off;
   print('-djpeg95', sprintf('wheel_theta%d.jpg',iver));
end


