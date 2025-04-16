
print_fig = 0;

% wheelset geometry parameters

nom_radius  =  460;
tape_circle =   70;
flange_back = 1360;
flback_ypos =  -70;

% wheelset position w.r.t. track system:

y_cm = 0;
z_cm = -nom_radius;
roll = 0;
roll_deg = roll * 180/pi;

% the wheel profile is given here with origin at the tape circle line

prw = read_profile('../../examples/MBench_S1002_v3.prw');

% position of tape circle line in wheelset-coordinates

ywset_tapcrc = flange_back/2 + tape_circle;
zwset_tapcrc = nom_radius;

RotRoll = [cos(roll), -sin(roll); 
           sin(roll),  cos(roll)];

yvec = RotRoll * [1;0];
zvec = RotRoll * [0;1];

% position where we want to place the wheel profile marker

yrw_ws = flange_back/2 - flback_ypos;
zrw_ws = nom_radius;

% make closed wheel contour using nominal radius

yprf = [ prw.ProfileY               ;  flipud(prw.ProfileY)              ];
zprf = [ prw.ProfileZ + zwset_tapcrc; -flipud(prw.ProfileZ)-zwset_tapcrc ];
yprf = [ yprf ; yprf(1) ]';
zprf = [ zprf ; zprf(1) ]';

% shift flange back to desired position

yprf = yprf + ywset_tapcrc;

figure(2); clf; hold on;
set(gca,'ydir','reverse');
axis equal;
axis([-200 1000 -1000 200]);
grid on;

% plot wheel profile

set(gca,'colororderindex',5);
tmp = RotRoll * [yprf; zprf];
plot(y_cm+tmp(1,:), z_cm+tmp(2,:));

% plot axle

set(gca,'colororderindex',5);
r_axle = 70;
y_axle = [ 0 flange_back/2 NaN max(yprf)+[0 70] ];
y_axle = [ y_axle, fliplr(y_axle) ];
z_axle = r_axle * [ ones(1,5), -ones(1,5) ];
tmp = RotRoll * [y_axle; z_axle];
plot(y_cm+tmp(1,:), z_cm+tmp(2,:));

y_axle = [ -180 0 NaN 0 -180];
z_axle = r_axle * [ 1 1 NaN -1 -1 ];
tmp = RotRoll * [y_axle; z_axle];
col = [1,1,1] - 0.5 * ([1,1,1] - matlab_color(5));
plot(y_cm+tmp(1,:), z_cm+tmp(2,:), 'color',col, 'linewidth',1.2);

% plot wheel-set coordinate axes

len = 180;
l=plot_cm([y_cm z_cm], len*0.15, len*0.85, roll_deg, 'k');
set(l,'linewidth',1);

text(y_cm+50, z_cm-100, 'O_{wset}', 'horizontalalignment','center', 'rotation', -roll_deg);
text(y_cm+len+15, z_cm-45, 'y_{wset}', 'verticalalignment','top', 'rotation', -roll_deg);
text(y_cm-70, z_cm+len-40, 'z_{wset}', 'verticalalignment','top', 'rotation', -roll_deg);

% plot flange-back distance <--> d/2

tmp = RotRoll * [   flange_back/2 * [  0  0  0  .32 NaN .68  1   1  1 ] ;
                      -nom_radius * [ .55 .45 .5 .5 NaN .5 .5 .45 .55 ] ];
plot(y_cm+tmp(1,:), z_cm+tmp(2,:), 'k', 'linewidth',1);
ix=[4,6];
text(y_cm+mean(tmp(1,ix)), z_cm+mean(tmp(2,ix)), 'd_{flng}/2', 'rotation', -roll_deg, ...
                                                                'horizontalalignment','center');

% plot tape circle line

set(gca,'colororderindex',2);
tmp = RotRoll * [ywset_tapcrc*[1 1]; zwset_tapcrc*[1 -1]];
plot( y_cm+tmp(1,:), z_cm+tmp(2,:) );

% plot marker for flange back on axle

tmp = RotRoll * [ flange_back/2 + [ 0 0 0 ]; [-30 0 30] ];
plot(y_cm+tmp(1,:), z_cm+tmp(2,:), 'k-');

text(y_cm+tmp(1,2)-5, z_cm+tmp(2,2)-0.85*nom_radius, 'flange-', 'rotation', -roll_deg, ...
                                                                'horizontalalignment','right');
text(y_cm+tmp(1,2)-5, z_cm+tmp(2,2)-0.70*nom_radius, 'back', 'rotation', -roll_deg, ...
                                                                'horizontalalignment','right');

% plot cross-hairs for tape circle on axle

tmp = RotRoll * [ ywset_tapcrc + [   0 30 NaN   0  0 ] ;
                                 [   0  0 NaN -30 30 ] ];
plot(y_cm+tmp(1,:), z_cm+tmp(2,:), 'k-');

text(y_cm+tmp(1,1)+20, z_cm+tmp(2,1)-0.84*nom_radius, 'y profile', 'rotation', -roll_deg);
text(y_cm+tmp(1,1)+20, z_cm+tmp(2,1)-0.69*nom_radius, ' reference', 'rotation', -roll_deg);
text(y_cm+tmp(1,1)+10, z_cm+tmp(2,1)+0.6*nom_radius, 'r_{nom,w}+dr', 'rotation', -roll_deg);

% plot wheel origin

len = 150;
tmp = RotRoll * [yrw_ws*[1 1]; zrw_ws*[1 -1]];
l = plot_cm([y_cm+tmp(1,1) z_cm+tmp(2,1)], len*0.15, len*0.85, roll_deg, 'k');
set(l, 'linewidth',1);

text(y_cm+tmp(1,1)+20, z_cm+tmp(2,1)+60, 'O_{wheel}', 'rotation', -roll_deg);

% plot indicators for flback_ypos

tmp = RotRoll * [  yrw_ws + flback_ypos * [ 1 1 NaN 1 1 ] ;
                   zrw_ws +          [ -20 20 NaN 50 90 ] ];
plot(y_cm+tmp(1,:), z_cm+tmp(2,:), 'k', 'linewidth',1);
text(y_cm+tmp(1,end)+50, z_cm+tmp(2,end)+0, 'y_{fbpos}', 'rotation', -roll_deg, ...
                                                                'horizontalalignment','right');

if (print_fig)
   axis off
   set(gcf,'paperpositionmode','auto');
   print -djpeg95 wheelset_def.jpg
end

