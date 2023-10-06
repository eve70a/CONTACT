
print_fig = 0;

addpath('../../matlab_intern');
set(0,'defaultaxeslinewidth',2);
set(0,'defaultlinelinewidth',2);
set(0,'defaulttextfontsize',15);
set(0,'defaultaxesfontsize',15);

r_roller =  800;
r_axle   =   70;
rnom_ws  =  460;
ax_len   = 1900;
cm_len   =  150;

gauge_height = 14;
gauge_width  = 1435;

% roller profile in rail coordinates, first point is center of top of rail

[yprf, zprf] = roller_profile(r_roller, r_axle);

% track origin in Isis-coordinates

yvec    = [1; 0];
zvec    = [0; 1];

% shift profile down to just touch the axle

zprf = zprf - max(zprf) + r_roller - r_axle;

% interpolate and insert new point at gauge measuring height

ixg  = find(zprf>gauge_height, 1, 'first');
fac  = (gauge_height-zprf(ixg-1)) / (zprf(ixg)-zprf(ixg-1));
yint = yprf(ixg-1) + fac * (yprf(ixg)-yprf(ixg-1));
zint = zprf(ixg-1) + fac * (zprf(ixg)-zprf(ixg-1));
yprf = [yprf(1:ixg-1), yint, yprf(ixg:end)];
zprf = [zprf(1:ixg-1), zint, zprf(ixg:end)];

% shift right to put gauge measuring point at half gauge width

y_rr = gauge_width/2 - yprf(ixg);
yprf = yprf + y_rr;


figure(1); clf; hold on;
set(gca,'ydir','reverse')
axis equal
axis([-1050 1050 -400 1000]);
grid on

% plot track plane

plot([-ax_len/2-50 ax_len/2+50], [0 0]);

% plot roller axle

set(gca,'colororderindex',2);
plot(ax_len/2*[-1 1 1 -1 -1], r_roller+r_axle*[1 1 -1 -1 1]);

% plot dashed line for roller axle

plot( (ax_len/2+50)*[-1 1], r_roller*[1 1], 'k--');

% plot right rail

set(gca,'colororderindex',2);
plot(yprf(:), zprf(:))

l=plot_cm([yprf(1), zprf(1)], cm_len*0.16, cm_len*0.84, 0, 'k');
set(l,'linewidth',1);
text(yprf(1)+30, zprf(1)-110, 'O_{rrail}', 'horizontalalignment','center');

% plot right rail again, flipped w.r.t. axle

zflip = 2*r_roller - zprf;

set(gca,'colororderindex',2);
plot(yprf(:), zflip(:));

% plot gauge measuring stop

ygaug = gauge_width/2+[0 0 -65 -65];
zgaug = gauge_height *[0 1   1   0];
plot(ygaug(:), zgaug(:), 'color',0.6*[1 1 1]);

% plot left rail : shifted by dr to separate from the track plane

hgt = max(zprf);
dr  = 30;
zlft = -dr + (hgt+dr) * zprf/hgt;

set(gca,'colororderindex',2);
plot(-yprf(:), zlft(:))

l=plot_cm([-yprf(1), zprf(1)], cm_len*0.16, cm_len*0.84, 0,'k');
set(l,'linewidth',1);
text(-yprf(1), zlft(1)-110, 'O_{lrail}', 'horizontalalignment','center');

% plot left rail again, flipped w.r.t. axle

set(gca,'colororderindex',2);
plot(-yprf(:), zflip(:));

% plot gauge measuring stop

ygaug = [ygaug, ygaug(1)];
zgaug = [zgaug, zgaug(1)];
plot(-ygaug(:), -dr+zgaug(:), 'color',0.6*[1 1 1]);

% plot mist over left side

l = fill([-1000 0 0 -1000 -1000], [1100 1100 -400 -400 1100], [1 1 1]);
set(l, 'edgecolor','none', 'facealpha', 0.7);

% plot coordinate axes

len = 250;
l = plot_axes([0,0], len, 0, -1, '   ','k',[],[],[],[],[],0);
set(l, 'linewidth',1);

text(0, -110, 'O_{track}', 'horizontalalignment','right');
text(len*yvec(1), len*yvec(2)-20, 'y_{track}', ...
                                        'verticalalignment', 'top');
text(len*zvec(1)-70, len*zvec(2)-20, 'z_{track}', ...
                                        'verticalalignment','top');

% plot axle center of mass

l=plot_cm([0 r_roller], cm_len*0.16, cm_len*0.84, 0, 'k');
set(l,'linewidth',1);
text(-40, r_roller(1)-80, 'O_{rol}', 'horizontalalignment','right');

% plot nominal radius

plot_arrow( [-400, r_roller/2], [0,  r_roller/2], 5, 0.6);
plot_arrow( [-400, r_roller/2], [0, -r_roller/2], 5, 0.6);
text(-370, 450, 'r_{nom,r}')

if (print_fig)
   axis off
   set(gcf,'paperpositionmode','auto');
   print -djpeg95 rollerrig_coords.jpg
end
