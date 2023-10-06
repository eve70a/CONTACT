
show_all = 1; % 11 profiles in one plot
show_one = 0; % nice picture of std profile
print_fig = 0;

v = ver; v = sscanf(v.Version, '%f');

if (v<=8.65)
   % Vortech82:
   set(0,'defaulttextfontsize',10);
   set(0,'defaultaxesfontsize',10);
   set(0,'defaultlinelinewidth',2);
   set(0,'defaultaxeslinewidth',1.5);
   fntsiz = 10;
else
   fntsiz = 15;
end

if (~exist('plot3d') | ~exist('plot_arrow'))
   addpath('../../matlab');
   addpath('../../matlab_intern');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read/create a variety of rail profiles

prr1 = read_profile('../../examples/MBench_UIC60_v3.prr');

% create prr2 with enlarged cant angle

prr2 = prr1;
roll = -10*pi/180;
rot  = [ cos(roll), -sin(roll); sin(roll), cos(roll) ];
tmp  = (rot * [prr1.ProfileY, prr1.ProfileZ]')';
prr2.ProfileY = tmp(:,1);
prr2.ProfileZ = tmp(:,2);

% change prr1 to remove cant & centered in y

roll = 0.0248;
rot  = [ cos(roll), -sin(roll); sin(roll), cos(roll) ];
tmp  = (rot * [prr1.ProfileY, prr1.ProfileZ]')';
prr1.ProfileY = tmp(:,1) - mean(tmp(:,1));
prr1.ProfileZ = (tmp(:,2) + flipud(tmp(:,2)))/ 2;

prr3 = prr1;
prr3.ProfileY = prr3.ProfileY - min(prr3.ProfileY);
prr3.ProfileZ = prr3.ProfileZ + 5;

idebug = 0;
prr4 = read_miniprof('../../../Projecten/Sentient/2019-wear-calculation/1-get-inputs/avg_low2.ban', idebug);
is_wheel = 1; mirror_y = 0; mirror_z = 1;
prr4 = modify_profile(prr4, is_wheel, mirror_y, mirror_z, [], [], [], [], idebug);

prr5 = read_miniprof('../../../meldingen/old-2017-2020/m18002-support-oldknow/2-miniprof-format/20050901-0211.ban');
prr5.ProfileY = [ -27.6; prr5.ProfileY;  29.4 ];
prr5.ProfileZ = [ -36.3; prr5.ProfileZ; -39.4 ];

prr6 = read_miniprof('../../../meldingen/old-2013-2016/m13014-pract-melis/1-from-bjorn/Section_35_l.rail', idebug);
is_wheel = 0; mirror_y = 1; mirror_z = 0;
prr6 = modify_profile(prr6, is_wheel, mirror_y, mirror_z, [], [], [], [], idebug);

prr7 = read_profile('../../../meldingen/old-2013-2016/m16006-wr-module/6-circ-profile/profiles/circ_r50.prr');

prr8 = read_profile('../../../meldingen/old-2013-2016/m16006-wr-module/6-circ-profile/profiles/flat_ang000.prr');

prr9 = read_simpack('../../../meldingen/old-2013-2016/m16006-wr-module/5-deep-groove/profiles/groove_r1500.prr');
prr9.ProfileY = [ -45; -45; prr9.ProfileY;  45; 45 ];
prr9.ProfileZ = [  20; -15; prr9.ProfileZ; -15; 20 ];

mirror_y = 1; mirror_z = 1;
prr10 = read_profile('../../../Projecten/Sentient/2021-wear-calculation/2-guard-rails/rail/NY2005-LOW_100-ARA-B_guardrail.ban', [], mirror_y, mirror_z);

mirror_y = 1; mirror_z = 1;
prr11 = read_profile('../../../meldingen/old-2013-2016/m16006-wr-module/report/fig/Profile_53R1_full.prr', [], mirror_y, mirror_z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot a variety of rail profiles in one picture

if (show_all)
   figure(1); clf;
   nrow = 3; ncol = 4;
   p = get(gcf,'position'); p(3:4) = [1000 480]; set(gcf,'position',p);

   subplot(nrow,ncol,1);
   plot( prr1.ProfileY, prr1.ProfileZ)
   plot_cm([0;0], 4, 10, 0, 'k');
   set(gca,'ydir','reverse');
   axis equal
   axis([-55 55 -15 55]); grid on
   text(0.09, 0.14, 'a)', 'units','normalized')

   subplot(nrow,ncol,2);
   plot(prr2.ProfileY, prr2.ProfileZ)
   plot_cm([0;0], 4, 10, 0, 'k');
   set(gca,'ydir','reverse');
   axis equal
   axis([-55 55 -15 55]); grid on
   text(0.05, 0.14, 'b)', 'units','normalized')

   subplot(nrow,ncol,3);
   plot(prr3.ProfileY, prr3.ProfileZ)
   plot_cm([0;0], 4, 10, 0, 'k');
   set(gca,'ydir','reverse');
   axis equal
   axis([-15 95 -15 55]); grid on
   text(0.05, 0.14, 'c)', 'units','normalized')

   subplot(nrow,ncol,4);
   plot(prr4.ProfileY, prr4.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([-65 45 -15 55]); grid on
   text(0.09, 0.14, 'd)', 'units','normalized')
   text(0.94, 0.14, 'track center \rightarrow', 'horizontalalignment','right', ...
                                                                'units','normalized', 'fontsize',fntsiz);
   set(gca,'xtick',[-60:30:30]);

   subplot(nrow,ncol,5);
   plot(prr5.ProfileY, prr5.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([-55 55 -55 15]); grid on
   text(0.09, 0.14, 'e)', 'units','normalized')

   subplot(nrow,ncol,6);
   plot(prr6.ProfileY, prr6.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([-45 65 -15 55]); grid on
   text(0.05, 0.14, 'f)', 'units','normalized')
   text(0.94, 0.14, 'switch/stock rail', 'horizontalalignment','right', ...
                                                                'units','normalized', 'fontsize',fntsiz);
   set(gca,'xtick',[-60:30:60]);

   subplot(nrow,ncol,7);
   plot(prr10.ProfileY, prr10.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([-140 80 -30 110]); grid on
   set(gca,'xtick',[-120:60:60], 'ytick',[0:40:80]);
   text(0.05, 0.14, 'g)', 'units','normalized');
   text(0.94, 0.14, 'guard/stock rail', 'horizontalalignment','right', ...
                                                                'units','normalized', 'fontsize',fntsiz);

   subplot(nrow,ncol,8);
   plot(prr11.ProfileY, prr11.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([-112 42 -20 78]); grid on
   set(gca,'xtick',[-90:30:60], 'ytick',[0:30:60]);
   text(0.05, 0.14, 'h)', 'units','normalized');
   text(0.96, 0.14, 'tram rail', 'horizontalalignment','right', 'units','normalized', 'fontsize',fntsiz);

   subplot(nrow,ncol,9);
   plot(prr7.ProfileY, prr7.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([-55 55 -15 55]); grid on
   text(0.12, 0.14, 'i)', 'units','normalized')
   text(0.91, 0.14, 'circular arc', 'horizontalalignment','right', 'units','normalized', 'fontsize',fntsiz);

   subplot(nrow,ncol,10);
   plot(prr8.ProfileY, prr8.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([-55 55 -15 55]); grid on
   text(0.09, 0.14, 'j)', 'units','normalized')
   text(0.94, 0.14, 'flat', 'horizontalalignment','right', 'units','normalized', 'fontsize',fntsiz);

   subplot(nrow,ncol,11);
   plot(prr9.ProfileY, prr9.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([-55 55 -35 35]); grid on
   text(0.12, 0.14, 'k)', 'units','normalized')
   text(0.94, 0.14, 'deep groove', 'horizontalalignment','right', 'units','normalized', 'fontsize',fntsiz);

   if (print_fig)
      set(gcf,'paperpositionmode','auto');
      print -djpeg95 rail_profiles.jpg
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot rail profile in same way as 'nice pictures' of 7-examples

if (show_one)
   figure(2); clf; hold on;
   plot( prr1.ProfileY, prr1.ProfileZ)
   plot_cm([0;0], 2.2, 8, 0, 'k');
   text(3, 6, 'm_{rr(rr)}')
   set(gca,'ydir','reverse');
   axis equal
   axis([-60 60 -20 50]); grid on
   xlabel('y_{rr} [mm]');
   ylabel('z_{rr} [mm]');

   if (print_fig)
      set(gcf,'paperpositionmode','auto');
      print -djpeg95 rr_profile.jpg
   end
end

