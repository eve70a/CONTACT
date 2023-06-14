
print_fig = 0;

addpath('../../matlab_intern');
set(0,'defaultaxeslinewidth',2);
set(0,'defaultlinelinewidth',2);
set(0,'defaulttextfontsize',15);
set(0,'defaultaxesfontsize',15);

% track geometry parameters:

cant_rai = atan(1/20); % 2.86 deg

gauge_height = 14;
gauge_width  = 1435;

RotCant = [cos(cant_rai), -sin(cant_rai); 
           sin(cant_rai),  cos(cant_rai)];

% get rail profile in rail coordinates, the first point is the rail origin:

[yprf,zprf] = rail_profile();

% rotate rail profile with -cant about rail origin

tmp  = RotCant' * [yprf; zprf];
yact = tmp(1,:); zact = tmp(2,:);

% shift up to just touch the track plane

zact = zact - min(zact);
z_rr = zact(1);

% interpolate and insert new point at gauge measuring height

ixg  = find(zact>gauge_height, 1, 'first');
fac  = (gauge_height-zact(ixg-1)) / (zact(ixg)-zact(ixg-1));
yint = yact(ixg-1) + fac * (yact(ixg)-yact(ixg-1));
zint = zact(ixg-1) + fac * (zact(ixg)-zact(ixg-1));
yact = [yact(1:ixg-1), yint, yact(ixg:end)];
zact = [zact(1:ixg-1), zint, zact(ixg:end)];

% shift right to put gauge measuring point at half gauge width

y_rr = gauge_width/2 - yact(ixg);
yact = yact + y_rr;

disp(sprintf('Rail origin at (%7.3f,%6.3f)', y_rr, z_rr));

figure(1); clf; hold on;
set(gca,'ydir','reverse')
axis equal
axis([ 650 850 -30 170]);
grid on
ylabel('z_{track}');
% xlabel('y_{track}');
text(860, 165, 'y_{track}');

% plot track plane

yvec   = [1; 0];
zvec   = [0; 1];
plane0 = -100*yvec;
plane1 =  900*yvec;
plot([plane0(1),plane1(1)], [plane0(2),plane1(2)])

% plot track coordinate axes

if (0==1)
   len = 50;
   l = plot_axes([0,0], len, 0, -1, '   ','k',[],[],[],[],[],0);
   text(0, -20, 'O_{track}', 'horizontalalignment','right');
   text(len*yvec(1),    len*yvec(2)- 0, 'y_{track}', 'verticalalignment','top');
   text(len*zvec(1)-40, len*zvec(2)-20, 'z_{track}', 'verticalalignment','top');
end

% plot right rail

set(gca,'colororderindex',2);
plot(yact, zact);
plot(yact(1), zact(1), 'k.', 'markersize',10);

% plot gauge measuring point

plot(yact(ixg)+[-45 0 0], gauge_height*[1 1 -0.2], 'color',0.6*[1 1 1]);
text(yact(ixg)-15, -3, '0.5 g_{width}', 'horizontalalignment','center', ...
                                        'verticalalignment','bottom');
text(yact(ixg)-25, gauge_height+2, 'g_{height}', ...
                'horizontalalignment','right');

% plot rail coordinate axes

len = 50;
plot_cm([y_rr z_rr], len*0.09, len*0.91, -cant_rai*180/pi, 'k');

yrail = RotCant' * [1; 0];
zrail = RotCant' * [0; 1];

text(y_rr+5, -3, 'O_{rail}', 'horizontalalignment','center', ...
                                        'verticalalignment','bottom');
text(y_rr+len*yrail(1)-2, len*yrail(2)+ 0, 'y_{rail}', ...
                                        'verticalalignment','top');
text(y_rr+len*zrail(1)+5, len*zrail(2)-32, 'z_{rail}', ...
                                        'verticalalignment','top');

if (print_fig)
   set(gcf,'paperpositionmode','auto');
   print -djpeg95 rail_coords.jpg
end

