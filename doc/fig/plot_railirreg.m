
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

% rail irregularity:

dphi_rai =  5 * pi/180;
dy_rai   = 15;
dz_rai   =  5;

RotCant = [cos(cant_rai), -sin(cant_rai); 
           sin(cant_rai),  cos(cant_rai)];
RotPhi  = [cos(dphi_rai), -sin(dphi_rai); 
           sin(dphi_rai),  cos(dphi_rai)];

% get rail profile in rail coordinates, the first point is the rail origin:

[yprf,zprf] = rail_profile();

% rotate rail profile with -cant about rail origin

tmp = RotCant' * [yprf; zprf];
ydsgn = tmp(1,:); zdsgn = tmp(2,:);

% shift up to just touch the track plane

zdsgn = zdsgn - min(zdsgn);
z_rr_dsgn = zdsgn(1);

% interpolate and insert new point at gauge measuring height

ixg  = find(zdsgn>gauge_height, 1, 'first');
fac  = (gauge_height-zdsgn(ixg-1)) / (zdsgn(ixg)-zdsgn(ixg-1));
yint = ydsgn(ixg-1) + fac * (ydsgn(ixg)-ydsgn(ixg-1));
zint = zdsgn(ixg-1) + fac * (zdsgn(ixg)-zdsgn(ixg-1));
ydsgn = [ydsgn(1:ixg-1), yint, ydsgn(ixg:end)];
zdsgn = [zdsgn(1:ixg-1), zint, zdsgn(ixg:end)];

% shift right to put gauge measuring point at half gauge width

y_rr_dsgn  = gauge_width/2 - ydsgn(ixg);
ydsgn = ydsgn + y_rr_dsgn;

disp(sprintf('Rail origin at (%7.3f,%6.3f)', y_rr_dsgn, z_rr_dsgn));

figure(2); clf; hold on;

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

% plot right rail at design position, with marker at rail origin

set(gca,'colororderindex',2);
plot(ydsgn, zdsgn, '--', 'linewidth', 1);

% plot rail coordinate axes for design position

len = 50;
h = plot_cm([y_rr_dsgn z_rr_dsgn], len*0.09, len*0.91, -cant_rai*180/pi, ...
        .5*[1 1 1]);
set(h, 'linewidth',1, 'linestyle','-');


% plot gauge measuring bar

plot(ydsgn(ixg)+[-45 0 0], gauge_height*[1 1 -0.2], 'color',0.6*[1 1 1]);
%text(ydsgn(ixg)-15, -3, '0.5 g_{width}', 'horizontalalignment','center', ...
%                                        'verticalalignment','bottom');
%text(ydsgn(ixg)-25, gauge_height+2, 'g_{height}', ...
%                'horizontalalignment','right');

% rotate & shift rail to actual position

tmp = RotPhi * [ydsgn-y_rr_dsgn; zdsgn-z_rr_dsgn];
yact = tmp(1,:) + y_rr_dsgn + dy_rai;
zact = tmp(2,:) + z_rr_dsgn + dz_rai;

y_rr_act = yact(1);
z_rr_act = zact(1);

% plot right rail at actual position

set(gca,'colororderindex',2);
plot(yact, zact, '-', 'linewidth', 2);

% plot rail coordinate axes for actual position

h = plot_cm([y_rr_act z_rr_act], len*0.09, len*0.91, ...
                (dphi_rai-cant_rai)*180/pi, 'k');

% yrail = RotPhi * RotCant' * [1; 0];
% zrail = RotPhi * RotCant' * [0; 1];
% text(y_rr_act+7, z_rr_act-6, 'O_{rrail}', ...
%                 'horizontalalignment','center', 'verticalalignment','bottom');
% text(y_rr_act+len*yrail(1)-2, z_rr_act+len*yrail(2)+0, ...
%                         'y_{rrail}', 'verticalalignment', 'top');
% text(y_rr_act+len*zrail(1)+5, z_rr_act+len*zrail(2)-32, ...
%                         'z_{rrail}', 'verticalalignment','top');

% plot displacements Delta y, Delta z

plot(y_rr_dsgn*[1 1], [-12 -6], 'k', 'linewidth', 1);
h1 = plot_arrow([y_rr_dsgn, -9], [dy_rai,0], 'k');
set(h1, 'linewidth',1);
text( y_rr_dsgn+dy_rai/2, -12, '\Delta{}y_{r}', ...
        'horizontalalignment','center', 'verticalalignment','bottom');

plot([830 836], [0 0], 'k', 'linewidth', 1);
h1 = plot_arrow([833,0], [0,dz_rai], 'k', 2, 2);
set(h1, 'linewidth', 1);
text( 836, 2, '\Delta{}z_{r}', 'verticalalignment','top');

if (print_fig)
   set(gcf,'paperpositionmode','auto');
   print -djpeg95 rail_irreg.jpg
end

