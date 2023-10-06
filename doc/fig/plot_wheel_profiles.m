
show_all = 1; % 9 profiles in one plot
show_one = 0; % nice picture of std profile
print_fig = 1;

v = ver; v = sscanf(v.Version, '%f');

if (v<=8.65)
   % Vortech82:
   fntsiz = 10;
   set(0,'defaulttextfontsize',fntsiz);
   set(0,'defaultaxesfontsize',fntsiz);
   set(0,'defaultlinelinewidth',2);
   set(0,'defaultaxeslinewidth',1.5);
else
   fntsiz = 15;
end

if (~exist('plot3d') | ~exist('plot_arrow'))
   addpath('../../../contact/matlab');
   addpath('../../../contact/matlab_intern');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read/create a variety of rail profiles

prw1 = read_profile('../../../../../contact/examples/MBench_S1002_v3.prw');

% create prw2 for left wheel

prw2 = prw1;
prw2.ProfileY = -prw2.ProfileY;

% change prw1 to remove cant & centered in y

% roll = 0.0248;
% rot  = [ cos(roll), -sin(roll); sin(roll), cos(roll) ];
% tmp  = (rot * [prw1.ProfileY, prw1.ProfileZ]')';
% prw1.ProfileY = tmp(:,1) - mean(tmp(:,1));
% prw1.ProfileZ = (tmp(:,2) + flipud(tmp(:,2)))/ 2;

prw3 = prw1;
prw3.ProfileY = prw3.ProfileY - min(prw3.ProfileY);
prw3.ProfileZ = prw3.ProfileZ + 5;

prw4 = read_miniprof('../../../../old-2017-2020/m18002-support-oldknow/2-miniprof-format/20050510-0891.whl');

prw5 = read_profile('../../../../../Projecten/Langeveld/m19001-evenaar-sommen/2-profielen/wheel1_p1.whl');

prw6 = read_profile('../../6-circ-profile/profiles/circ_r50.prw');
prw7 = read_profile('../../6-circ-profile/profiles/flat_ang200.prw');

prw8 = read_simpack('../../../../../Projecten/SouthObserv/1-log-profile/wheel_log_v1.prw');
% scale to 120 mm width, flip z positive downwards
prw8.ProfileY = (60/87) * prw8.ProfileY;
prw8.ProfileZ = -prw8.ProfileZ;

% create double-flanged crane wheel with three rounded corners

rcirc =  5;
th = [90 : -3 : 6] * pi/180;
ycirc = rcirc * cos(th);
zcirc = rcirc * sin(th);

hlf_wid = 140/2;
fl_wid  =  25;
fl_hgt  =  30;

% p3: field side, top of flange
p3_z = fl_hgt - rcirc; p3_y = hlf_wid - rcirc;           
% p2: rail center, top of flange
p2_z = fl_hgt - rcirc; p2_y = hlf_wid - fl_wid + rcirc;  
% p1: inside corner of flange
p1_z =          rcirc;
dz = (p2_z+zcirc(end)) - (p1_z-zcirc(end));
dy = dz * tan(6*pi/180);
p1_y = p2_y - dy - 2*ycirc(end);

prw9.ProfileY = [ 0, p1_y+ycirc, p2_y-fliplr(ycirc), p3_y+ycirc, hlf_wid ];
prw9.ProfileZ = [ 0, p1_z-zcirc, p2_z+fliplr(zcirc), p3_z+zcirc,  -10  ];
prw9.ProfileY = [ -fliplr(prw9.ProfileY), prw9.ProfileY(2:end) ];
prw9.ProfileZ = [  fliplr(prw9.ProfileZ), prw9.ProfileZ(2:end) ];

is_wheel = 1;
prw10 = read_profile('../../../../../contact/testbank/profiles/fr_ol.csv', is_wheel);
prw10.ProfileY = prw10.ProfileY - mean(prw10.ProfileY);
prw10.ProfileY = [prw10.ProfileY(12); prw10.ProfileY(12:end); prw10.ProfileY(end)];
prw10.ProfileZ = [prw10.ProfileZ(12)-6.5; prw10.ProfileZ(12:end); prw10.ProfileZ(end)-5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot a variety of rail profiles in one picture

if (show_all)
   figure(1); clf;
   % p = get(gcf,'position'); p(3) = 3.2*p(4); set(gcf,'position',p);
   nrow = 3; ncol = 4;
   p = get(gcf,'position'); p(3:4) = [1000 480]; set(gcf,'position',p);

   subplot(nrow,ncol,1);
   plot( prw1.ProfileY, prw1.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([-80 70 -20 50]); grid on
   plot_cm([0;0], 4, 10, 0, 'k');
   text(0.92, 0.18, 'a)', 'units','normalized', 'horizontalalignment','right')

   subplot(nrow,ncol,2);
   plot(prw2.ProfileY, prw2.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([-70 80 -20 50]); grid on
   plot_cm([0;0], 4, 10, 0, 'k');
   text(0.08, 0.18, 'b)', 'units','normalized')

   subplot(nrow,ncol,3);
   plot(prw3.ProfileY, prw3.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([-10 140 -20 50]); grid on
   plot_cm([0;0], 4, 10, 0, 'k');
   text(0.92, 0.18, 'c)', 'units','normalized', 'horizontalalignment','right')

   subplot(nrow,ncol,4);
   plot(prw4.ProfileY, prw4.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([0 150 -50 20]); grid on
   text(0.92, 0.18, 'd)', 'units','normalized', 'horizontalalignment','right')

   subplot(nrow,ncol,5);
   plot(prw6.ProfileY, prw6.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([-75 75 -55 15]); grid on
   text(0.08, 0.18, 'e)', 'units','normalized')

   subplot(nrow,ncol,6);
   plot(prw7.ProfileY, prw7.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([-75 75 -30 40]); grid on
   text(0.08, 0.18, 'f)', 'units','normalized')

   subplot(nrow,ncol,7);
   plot(prw8.ProfileY, prw8.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([-75 75 -35 35]); grid on
   text(0.08, 0.18, 'g)', 'units','normalized')

   subplot(nrow,ncol,8);
   plot(prw10.ProfileY, prw10.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([-15 15 -12 2]); grid on
   text(0.94, 0.18, 'h)', 'units','normalized', 'horizontalalignment','right')

   subplot(nrow,ncol,9);
   plot(prw9.ProfileY, prw9.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([-75 75 -20 50]); grid on
   text(0.25, 0.18, 'i)', 'units','normalized');

   if (print_fig)
      set(gcf,'paperpositionmode','auto');
      print -djpeg95 wheel_profiles.jpg
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot rail profile in same way as 'nice pictures' of 7-examples

if (show_one)
   figure(2); clf; hold on;
   plot( prw9.ProfileY, prw9.ProfileZ)
   plot([p1_y, p2_y, p3_y], [p1_z, p2_z, p3_z], '*');
   % plot_cm([0;0], 2.2, 8, 0, 'k');
   % text(3, 6, 'm_{rr(rr)}')
   % plot([p1_y+ycirc(end), p2_y-ycirc(end)], [p1_z-zcirc(end) p2_z+zcirc(end)], '*')
   % th = 84 * pi/180;
   % v = [cos(th); sin(th)];
   % plot(p1_y+ycirc(end)+[0 30]*v(1), p1_z-zcirc(end)+[0 30]*v(2))
   set(gca,'ydir','reverse');
   axis equal
   axis([ -70 70 -10 35]); grid on
   xlabel('y_{rr} [mm]');
   ylabel('z_{rr} [mm]');

   if (print_fig)
      set(gcf,'paperpositionmode','auto');
      print -djpeg95 rr_profile.jpg
   end
end

