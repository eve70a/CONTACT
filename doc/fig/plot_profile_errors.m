
show_all = 1; % all profiles in one plot
show_one = 0; % single picture of one profile
print_fig = 1;

v = ver('matlab'); v = sscanf(v.Version, '%f');

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

sentient19 = '../../../Projecten/Sentient/2019-wear-calculation/';
sentient20 = '../../../Projecten/Sentient/2020-wear-calculation/';
sentient21 = '../../../Projecten/Sentient/2021-wear-calculation/';
sentient22 = '../../../Projecten/Sentient/2022-wear-calculation/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read/create a variety of rail profiles showing challenging features

prf1 = read_profile([sentient22,'2-prog2-error/case1_SB/rail/current_left_rail.ban'], [], 1);
prf2 = read_profile([sentient19,'1-get-inputs/avg_HR.BAN']);
prf3 = read_profile([sentient22,'1-damage-nrs\1_CN_SHARP\wheel\CGEX_1419_L3.trc.trc.whl']);

prf4 = read_profile('../../testbank/profiles/Car7422_0001r.whl');
prf5 = read_profile('../../testbank/profiles/rail_left_iter288.ban', [], 1);
prf6 = read_profile('../../testbank/profiles/groove_r1025.prr');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot a variety of rail profiles in one picture

if (show_all)
   figure(1); clf;
   nrow = 2; ncol = 3;
   p = get(gcf,'position'); p(3:4) = [1000 480]; set(gcf,'position',p);

   % subplot axes: aspect ratio 5 : 4

   subplot(nrow,ncol,1);
   plot( prf1.ProfileY, prf1.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([-49.2 -46.7 31 33]); grid on
   text(0.09, 0.14, 'a)', 'units','normalized')

   subplot(nrow,ncol,2);
   plot( prf2.ProfileY, prf2.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([29.237 29.287 2.528 2.568]); grid on
   text(0.09, 0.14, 'b)', 'units','normalized')

   subplot(nrow,ncol,3);
   plot( prf3.ProfileY, prf3.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([0 40  2 34]); grid on
   text(0.09, 0.14, 'c)', 'units','normalized')

   subplot(nrow,ncol,4);
   plot( prf4.ProfileY, prf4.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([16 19.5 26.8 29.6]); grid on
   text(0.09, 0.14, 'd)', 'units','normalized')

   subplot(nrow,ncol,5);
   plot( prf5.ProfileY, prf5.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([-32 -29 14.5 16.9]); grid on
   text(0.09, 0.14, 'e)', 'units','normalized')

   subplot(nrow,ncol,6);
   plot( prf6.ProfileY, prf6.ProfileZ)
   set(gca,'ydir','reverse');
   axis equal
   axis([-18 -3 -15 -3]); grid on
   text(0.09, 0.14, 'f)', 'units','normalized')

   if (print_fig)
      set(gcf,'paperpositionmode','auto');
      print -djpeg95 profile_errors.jpg
   end
end

if (show_one>0)
   figure(2); clf; hold on;
   eval(sprintf('prf = prf%d;', show_one));

   plot( prf.ProfileY, prf.ProfileZ)
   % plot_cm([0;0], 4, 10, 0, 'k');
   set(gca,'ydir','reverse');
   axis equal
   axis([-55 55 -15 55]); grid on
   text(0.09, 0.14, 'a)', 'units','normalized')
end
