
print_fig = 1;

addpath('../../matlab_intern');
set(0,'defaultaxeslinewidth',2);
set(0,'defaultlinelinewidth',2);
set(0,'defaulttextfontsize',15);
set(0,'defaultaxesfontsize',15);

cant_rai = 15 * pi/180;

RotCant = [cos(cant_rai), -sin(cant_rai); 
           sin(cant_rai),  cos(cant_rai)];

% rail profile in rail coordinates
% first point == center, top of rail

[yprf,zprf] = rail_profile();

rail_height = max(zprf) - min(zprf);

% rotate rail profile with -cant about rail origin

tmp  = RotCant' * [yprf; zprf];
yact = tmp(1,:);
zact = tmp(2,:);

% shift up to just touch the track plane

zact = zact - min(zact);
z_rr = zact(1);  % center of uncanted rail

% interpolate and insert new point at gauge measuring height

gauge_height = 20;
ixg  = find(zact>gauge_height, 1, 'first');
fac  = (gauge_height-zact(ixg-1)) / (zact(ixg)-zact(ixg-1));
yint = yact(ixg-1) + fac * (yact(ixg)-yact(ixg-1));
zint = zact(ixg-1) + fac * (zact(ixg)-zact(ixg-1));
yact = [yact(1:ixg-1), yint, yact(ixg:end)];
zact = [zact(1:ixg-1), zint, zact(ixg:end)];

% shift right to put gauge measuring point at half gauge width

gauge_width = 1435;
y_rr        = gauge_width/2 - yact(ixg);
yact = yact + y_rr;

figure(1); clf; hold on;
p=get(gcf,'position'); p(3:4)=[710 420]; set(gcf,'position',p);

% plot track plane

yvec = [1; 0];
zvec = [0; 1];
plane0 = - 900*yvec;
plane1 =   900*yvec;
plot([0, plane1(1)], [0, plane1(2)])                               % right half

col = [1,1,1] - 0.5 * ([1,1,1] - matlab_color(1));
plot([plane0(1), 0], [plane0(2), 0], 'color',col, 'linewidth',1.2) % left half

set(gca,'ydir','reverse')
axis equal
axis([-950 950 -400 600]);
grid on

% plot coordinate axes

len = 180;
l = plot_axes([0,0], len, 0, -1, '   ','k',[],[],[],[],[],0);
set(l(1:14), 'linewidth',1);

text(len*yvec(1),    len*yvec(2)-0, 'y_{track}', 'verticalalignment', 'top');
text(len*zvec(1)-70, len*zvec(2)-20, 'z_{track}', 'verticalalignment','top');

% plot right rail

tmp = [yact; zact];

set(gca,'colororderindex',2);
plot(tmp(1,:), tmp(2,:))

% plot left rail

tmp = [-yact; zact];

col = [1,1,1] - 0.5 * ([1,1,1] - matlab_color(2));
plot(tmp(1,:), tmp(2,:), 'color',col, 'linewidth',1.2);

% plot cant angles

tmp = RotCant' * [ 0, 0; -50, rail_height ];
tmp(:,3) = tmp(:,2) + [ 0; -50-rail_height];
plot( y_rr+tmp(1,:),  z_rr+tmp(2,:), 'k--', 'linewidth',1);
col = 0.5 * [1 1 1];
plot(-y_rr-tmp(1,:),  z_rr+tmp(2,:), '--', 'color',col, 'linewidth',1);

circ_arrow([ y_rr+tmp(1,2);z_rr+tmp(2,2)], rail_height+50, ...
                        -90, -90-cant_rai*180/pi, 5, 2, 1.3);
col = [1,1,1] - 0.5 * ([1,1,1] - matlab_color(5));
l = circ_arrow([-y_rr-tmp(1,2);z_rr+tmp(2,2)], rail_height+50, ...
                        -90, -90+cant_rai*180/pi, col, 2, 1.3);
set(l, 'linewidth',1.2);
text(y_rr, z_rr+rail_height+80, 'rail cant>0', 'horizontalalignment','center');

% plot right gauge measuring point / gauge stop

tmp = [ gauge_width/2+[0 0 -65 -65];
        gauge_height*[0 1 1 0] ];
plot(tmp(1,:), tmp(2,:), 'color',0.6*[1 1 1]);

% plot left gauge measuring point / gauge stop

tmp = [ -gauge_width/2-[0 0 -65 -65];
         gauge_height*[0 1 1 0] ];
plot(tmp(1,:), tmp(2,:), 'color',0.8*[1 1 1], 'linewidth',1.2);

% plot gauge width G/2 <--> G/2

gw = gauge_width;
plot( gw/2*[1 1], [-80 -140], 'k-', 'linewidth', 1);
plot(    0*[1 1], [-80 -140], 'k-', 'linewidth', 1);
text( gw/4, -110, '$G/2$', 'horizontalalignment','center', 'interpreter','latex');
l = plot_arrow([0.30*gw, -110], [ 0.2*gw, 0], 'k', 0.5); set(l,'linewidth',1);
plot([0 0.2*gw], -110*[1 1], 'k', 'linewidth',1);

col = 0.5 * [ 1 1 1 ];
plot(-gw/2*[1 1], [-80 -140], '-', 'color',col, 'linewidth', 1);
text(-0.235*gw, -110, '$-G/2$', 'color',col, 'horizontalalignment','center', 'interpreter','latex');
l = plot_arrow([-0.30*gw, -110], [-0.2*gw, 0], col, 0.5); set(l,'linewidth',0.8);
plot([0 -0.17*gw], -110*[1 1], 'color',col, 'linewidth',0.8);

if (print_fig)
   axis off
   set(gcf,'paperpositionmode','auto');
   print -djpeg95 track_coords.jpg
end
