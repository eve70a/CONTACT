
print_fig = 0;

addpath('../../matlab_intern');
set(0,'defaultaxeslinewidth',2);
set(0,'defaultlinelinewidth',2);
set(0,'defaulttextfontsize',15);
set(0,'defaultaxesfontsize',15);

% track origin in Isis-coordinates

ytrk = 0; ztrk = 0;
elev_trk = 10 * pi/180;
cant_rai =  0 * pi/180;

RotElev = [cos(elev_trk), -sin(elev_trk); 
           sin(elev_trk),  cos(elev_trk)];
RotCant = [cos(cant_rai), -sin(cant_rai); 
           sin(cant_rai),  cos(cant_rai)];

% rail profile in rail coordinates

[yprf,zprf] = rail_profile();

% total height; width of rail head

rail_height = max(zprf) - min(zprf);
ix = find(zprf<rail_height/2);
rail_width  = max(yprf(ix)) - min(yprf(ix));

% rotate rail profile with -cant about rail origin

tmp  = RotCant' * [yprf; zprf];
yact = tmp(1,:); zact = tmp(2,:);

% shift up to just touch the track plane

zact       = zact - min(zact);
ztrk_orail = zact(1);

% interpolate and insert new point at gauge measuring height

gauge_height = 14;
ixg  = find(zact>gauge_height, 1, 'first');
fac  = (gauge_height-zact(ixg-1)) / (zact(ixg)-zact(ixg-1));
yint = yact(ixg-1) + fac * (yact(ixg)-yact(ixg-1));
zint = zact(ixg-1) + fac * (zact(ixg)-zact(ixg-1));
yact = [yact(1:ixg-1), yint, yact(ixg:end)];
zact = [zact(1:ixg-1), zint, zact(ixg:end)];

% shift right to put gauge measuring point at half gauge width

gauge_width = 1435;
ytrk_orail  = gauge_width/2 - yact(ixg);
yact = yact + ytrk_orail;

figure(1); clf; hold on;
set(gca,'ydir','reverse')
axis equal
axis([-950 950 -400 600]);
grid on

% plot Isys coordinate axes

if (0==1)
   len = 200;
   ang =   0;
   plot_axes([-800, 300], len, ang, -1, '   ', 'k');
   text(-750, 370, 'O_{isys}')
end

% plot track plane

yvec = RotElev * [1; 0];
zvec = RotElev * [0; 1];
plane0 = [ytrk;ztrk]- 900*yvec;
plane1 = [ytrk;ztrk]+ 900*yvec;
set(gca,'colororderindex',1);
plot([plane0(1),plane1(1)], [plane0(2),plane1(2)])

text(700, -100, 'track plane', 'horizontalalignment','right')
plot([360, 300], [-60, 50], 'k-','linewidth',1);

% plot right rail

tmp = RotElev * [yact; zact];

set(gca,'colororderindex',2);
plot(ytrk+tmp(1,:), ztrk+tmp(2,:))
% plot(ytrk+tmp(1,1), ztrk+tmp(2,1), 'k.', 'markersize',10);

dash = RotElev * RotCant' * [ [0 0] ; [-20 rail_height+50] ];
plot(ytrk+tmp(1,1)+dash(1,:), ztrk+tmp(2,1)+dash(2,:), 'k--', 'linewidth',1);

% plot left rail

tmp = RotElev * [-yact; zact];

set(gca,'colororderindex',2);
plot(ytrk+tmp(1,:), ztrk+tmp(2,:))
% plot(ytrk+tmp(1,1), ztrk+tmp(2,1), 'k.', 'markersize',10);

dash = RotElev * RotCant * [ [0 0] ; [-20 rail_height+50] ];
plot(ytrk+tmp(1,1)+dash(1,:), ztrk+tmp(2,1)+dash(2,:), 'k--', 'linewidth',1);

% plot gauge measuring point on right rail

if (0==1)
   tmp = RotElev * [ gauge_width/2+[0 0 -65 -65];
                     gauge_height*[0 1 1 0] ];
   plot(ytrk+tmp(1,:), ztrk+tmp(2,:), 'color',0.6*[1 1 1]);

% plot gauge measuring point on left rail

   tmp = RotElev * [ -gauge_width/2-[0 0 -65 -65];
                      gauge_height*[0 1 1 0] ];
   plot(ytrk+tmp(1,:), ztrk+tmp(2,:), 'color',0.6*[1 1 1]);
end

% plot gauge width

y2 = gauge_width/2;
y1 = y2 + rail_width;

tmp = RotElev * [ [  y2  y2  y2  y1  y1  y1  y1 -y1 -y1 -y1 -y1 -y2 -y2 -y2 ];
                  [ -30 -70 -50 -50 -30 -70 -50 -50 -70 -30 -50 -50 -70 -30 ] ];
plot(ytrk+tmp(1,:), ztrk+tmp(2,:), 'k-', 'linewidth',1);
ang = 180/pi * -elev_trk;
text(ytrk+mean(tmp(1,3:4))+13, ztrk+mean(tmp(2,3:4))-60, 'w', ...
        'rotation', ang, 'horizontalalignment','center');
text(ytrk+mean(tmp(1,7:8))+13, ztrk+mean(tmp(2,7:8))-60, 'G', ...
        'rotation', ang, 'horizontalalignment','center');

% plot rail center distance

tmp = RotElev * [ (gauge_width+rail_width)/2 * [1 -1] ;
                  (rail_height+30) * [ 1 1 ] ];
plot(ytrk+tmp(1,:), ztrk+tmp(2,:), 'k', 'linewidth', 1);
text(ytrk+mean(tmp(1,:))+13, ztrk+mean(tmp(2,:))-60, 'd', ...
        'rotation', ang, 'horizontalalignment','center');

plot(ytrk+tmp(1,1)+[0 -600], ztrk+tmp(2,1)*[1 1], 'k', 'linewidth',1);
circ_arrow(tmp(:,1), 500, -180, -180+elev_trk*180/pi, 1);
text(170, 260, '\phi', 'horizontalalignment','right')

text(ytrk+tmp(1,1)-10, ztrk+tmp(2,1)+50, 'horizontal plane', ...
                'horizontalalignment','right');

if (print_fig)
   axis off
   set(gcf,'paperpositionmode','auto');
   print -djpeg75 isys_coords.jpg
end
