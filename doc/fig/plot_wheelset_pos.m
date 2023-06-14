
print_fig = 0;

if (~exist('read_simpack'))
   addpath('../../matlab');
end
if (~exist('plot_arrow'))
   addpath('../../matlab_intern');
end

% wheelset geometry parameters

nom_radius  =  460;
tape_circle =   70;
flange_back = 1360;
flback_ypos =  -70;

% the wheel profile is given with origin at the tape circle line

prw = read_simpack('../../examples/MBench_S1002_v3.prw');

figure(2); clf; hold on;
set(gca,'ydir','reverse');
axis equal;
axis([-200 1000 -1000 200]);
grid on;

% plot track coordinate axes

len = 200;
l = plot_axes([0,0], len, 0, -1, '   ','k',[],[],[],[],[],0);
set(l, 'linewidth',1);

text(0, -100, 'O_{track}', 'horizontalalignment','right');
text(len, -20, 'y_{track}', 'verticalalignment', 'top');
text(-70, len-20, 'z_{track}', 'verticalalignment','top');

% plot wheelset at initial (reference) and actual positions

for ipos = 1 : 2

   % wheelset position w.r.t. track system:

   if (ipos==1)
      y_cm = 0;
      z_cm = -nom_radius;
      roll = 0;
   else
      y_cm = 50;
      z_cm = -nom_radius - 10;
      roll = 2.0 * pi/180;
   end
   roll_deg = roll * 180/pi;

   % position of tape circle line in wheelset-coordinates

   ywset_tapcrc = flange_back/2 + tape_circle;
   zwset_tapcrc = nom_radius;

   RotRoll = [cos(roll), -sin(roll); 
              sin(roll),  cos(roll)];

   yvec = RotRoll * [1;0];
   zvec = RotRoll * [0;1];

   % position where we want to pretend that the origin was located

   yrw_ws = flange_back/2 - flback_ypos;
   zrw_ws = nom_radius;

   % make closed wheel contour using nominal radius

   yprf = [ prw.ProfileY               ;  flipud(prw.ProfileY)              ];
   zprf = [ prw.ProfileZ + zwset_tapcrc; -flipud(prw.ProfileZ)-zwset_tapcrc ];
   yprf = [ yprf ; yprf(1) ]';
   zprf = [ zprf ; zprf(1) ]';

   % shift flange back to desired position

   yprf = yprf + ywset_tapcrc;

   % plot wheel profile

   set(gca,'colororderindex',5);
   tmp = RotRoll * [yprf; zprf];
   l = plot(y_cm+tmp(1,:), z_cm+tmp(2,:));
   if (ipos==1), set(l,'linewidth',1); end

   % plot axle

   set(gca,'colororderindex',5);
   r_axle = 70;
   y_axle = [ 0 flange_back/2 NaN max(yprf)+[0 70] ];
   y_axle = [ y_axle, fliplr(y_axle) ];
   z_axle = r_axle * [ ones(1,5), -ones(1,5) ];
   tmp = RotRoll * [y_axle; z_axle];
   l = plot(y_cm+tmp(1,:), z_cm+tmp(2,:));
   if (ipos==1), set(l,'linewidth',1); end

   y_axle = [ -180 0 NaN 0 -180];
   z_axle = r_axle * [ 1 1 NaN -1 -1 ];
   tmp = RotRoll * [y_axle; z_axle];
   col = [1,1,1] - 0.5 * ([1,1,1] - matlab_color(5));
   l = plot(y_cm+tmp(1,:), z_cm+tmp(2,:), 'color',col, 'linewidth',1.2);
   if (ipos==1), set(l,'linewidth',0.5); end

   % plot wheel-set coordinate axes

   col = [.6; 0] * [1 1 1];

   len = 180;
   l=plot_cm([y_cm z_cm], len*0.15, len*0.85, roll_deg, col(ipos,:));
   set(l,'linewidth',1);

   if (ipos==2)
      text(y_cm+50, z_cm-100, 'O_{wset}', 'horizontalalignment','center', ...
                                                   'rotation', -roll_deg);
      text(y_cm+len+15, z_cm-45, 'y_{wset}', 'verticalalignment','top', ...
                                                   'rotation', -roll_deg);
      text(y_cm-70, z_cm+len-40, 'z_{wset}', 'verticalalignment','top', ...
                                                   'rotation', -roll_deg);
   end

   % plot marker for flange back on axle

   tmp = RotRoll * [ flange_back/2 + [ 0 0 0 ]; [-30 0 30] ];
   plot(y_cm+tmp(1,:), z_cm+tmp(2,:), '-', 'color', col(ipos,:));

   % plot cross-hairs for tape circle on axle

   tmp = RotRoll * [ ywset_tapcrc + [   0 30 NaN   0  0 ] ;
                                    [   0  0 NaN -30 30 ] ];
   plot(y_cm+tmp(1,:), z_cm+tmp(2,:), '-', 'color', col(ipos,:));

   % plot wheel origin

   len = 150;
   tmp = RotRoll * [yrw_ws*[1 1]; zrw_ws*[1 -1]];
   l = plot_cm([y_cm+tmp(1,1) z_cm+tmp(2,1)], len*0.15, len*0.85, ...
                                                roll_deg, col(ipos,:));
   set(l, 'linewidth',1);

   if (ipos==2)
      text(y_cm+tmp(1,1)+20, z_cm+tmp(2,1)+60, 'O_{wheel}', ...
                                                'rotation', -roll_deg);
   end

end % ipos

if (print_fig)
   axis off
   set(gcf,'paperpositionmode','auto');
   print -djpeg95 wheelset_pos.jpg
end

