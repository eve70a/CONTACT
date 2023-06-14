
if (~exist('loadcase'))
   addpath('../matlab');
end
if (~exist('plot_arrow'))
   addpath('../matlab_intern');
end

expnam = strvcat('m3_tempdep', 'conv_temp', 'm1_tempdep', 'euro_tractcurv');
expnam = 'm3_tempdep';

pause_after_plot = 1 * (size(expnam,1)>1);
print_fig = 1;

%  determine version information on the Matlab graphics system

old_graphics = verLessThan('matlab', '8.4.0');
new_graphics = ~old_graphics;
set(0,'defaultfigurepaperpositionmode','auto');

% Ertz testcase met temperature-dependent friction

if (~isempty(strmatch('m3_tempdep',expnam)))
   ncase = 3;
   tcentr = []; mu = []; px = []; srel = [];
   for icase = 1 : ncase
      s = loadcase('m3_tempdep', icase);
      [~,iy0] = min(abs(s.y));
      tcentr = [tcentr; s.temp1(iy0,:)];
      mu     = [mu;     s.mu(iy0,:)];
      px     = [px;     s.px(iy0,:)];
      srel   = [srel;   s.srel(iy0,:)];
   end
   a = 6;

   figure(6); clf;
   p = get(gcf,'position'); p(3:4) = [800 600]; set(gcf,'position',p);
   subplot(2,2,1); hold on;
   plot(s.x/a, px, '-');
   ax = axis; ax(1:2) = [-1.2 6]; axis(ax); grid on
   legend('Constant \mu', '\Delta\mu_{heat}=-0.18', '\Delta\mu_{heat}=+0.18');
   % legend('reference','decreasing \mu', 'increasing \mu');
   title('traction');
   ylabel('p_x [N/mm^2]');

   subplot(2,2,2);
   plot(s.x/a, srel, '-');
   axis([-1.2 6 0.025 0.040]); grid on
   title('relative slip');
   ylabel('s [-]');

   subplot(2,2,3); hold on;
   plot(s.x/a, tcentr, '-');
   plot([-1 -1], [175 225], 'k--');
   plot([ 1  1], [175 225], 'k--');
   axis([-1.2 6 0 250]); grid on
   title('wheel temperature');
   xlabel('x/a [-]'); ylabel('T_a [{}^\circ{}C]');

   subplot(2,2,4); hold on
   plot(s.x/a, mu, '-');
   plot([-1 -1], [0.375 0.425], 'k--');
   plot([ 1  1], [0.375 0.425], 'k--');
   ax = axis; ax(1:4) = [-1.2 6 0.2 0.42]; axis(ax); grid on
   title('friction coefficient');
   xlabel('x/a [-]'); ylabel('\mu [-]');

   if (print_fig)
      print -djpeg95 figs/m3_tempdep.jpg
   end

   if (pause_after_plot), pause; end
   clear a ax icase iy0 l mu ncase px s srel tcentr
end

% Ertz testcase, kleine Delta T zorgt voor convergentie-problemen

if (~isempty(strmatch('conv_temp',expnam)))
   s1 = loadcase('conv_temp', 1);
   s2 = loadcase('conv_temp', 2);
   s3 = loadcase('conv_temp', 3);
   s4 = loadcase('conv_temp', 4);

   opt = plot3d;
   opt.field = 'mu';
   figure(1); clf;
   plot3d(s2, opt);

   figure(2); clf;
   plot3d(s4, opt);

   [~,iy0] = min(abs(s1.y));
   figure(3); clf; hold on;
   plot(s1.x, [s1.temp1(iy0,:); s2.temp1(iy0,:); ...
                                s3.temp1(iy0,:); s4.temp1(iy0,:)]);
   axis([-6 12 0 160]); grid on;
   xlabel('x_c [mm]');
   ylabel('Temperature');
   legend('\omega_{slp} = 1.0', '\omega_{slp} = 0.5', '\omega_{slp} = 0.2', ...
                'memdst=0.1\mu{}m', 'location','southeast');

   if (pause_after_plot), pause; end
   clear iy0 opt s1 s2 s3 s4;
end

% Roller on flat plane o.b.v. module 1, w/r profiles

if (~isempty(strmatch('m1_tempdep',expnam)))
   s1 = loadcase('m1_tempdep', 1);
   s2 = loadcase('m1_tempdep', 2);
   s3 = loadcase('m1_tempdep', 3);

   [~,iy0] = min(abs(s1.y));
   figure(1); clf; hold on;
   l = plot(s1.x, [s1.temp1(iy0,:); s2.temp1(iy0,:); s3.temp1(iy0,:)]);
   set(l(3),'linestyle','--');
   axis([-3 3 0 50]); grid on;
   xlabel('x_c [mm]');
   title('Surface temperature');

   figure(2); clf; hold on;
   l = plot(s1.x, [s1.mu(iy0,:); s2.mu(iy0,:); s3.mu(iy0,:)]);
   set(l(3),'linestyle','--');
   axis([-3 3 0.275 0.305]); grid on;
   xlabel('x_c [mm]');
   title('Coefficient of friction');

   if (print_fig), 
      figure(1); print -djpeg95 figs/m1_temp.jpg
      figure(2); print -djpeg95 figs/m1_frict.jpg
   end

   if (pause_after_plot), pause; end
   clear iy0 l s1 s2 s3;
end

% Traction curve cf. Eurosprinter testcase

if (~isempty(strmatch('euro_tractcurv',expnam)))
   fstat = 0.33;

   [creep, force]= parse_out3('euro_tractcurv.out');
   cksi =          reshape(creep(:,2), 30, 4);
   fx   = -fstat * reshape(force(:,2), 30, 4);
   t1   =          reshape(force(:,7), 30, 4);
   t2   =          reshape(force(:,8), 30, 4);

   figure(1); clf;
   l = plot(100*cksi, fx);
   set(l(1), 'color','b', 'linestyle','--');
   axis([0 5 0 0.4]); grid on;
   xlabel('Longitudinal creepage \xi [%]');
   ylabel('Traction coefficient F_x/F_n [-]');
   legend('Original Kalker theory', ...
           'Temp.dep., \Delta\mu= 0.19, \Delta{}T=1200', ...
           'Temp.dep., \Delta\mu=-0.19, \Delta{}T=1200', ...
           'Temp.dep., \Delta\mu=-0.19, \Delta{}T= 400', 'location','southeast')

   figure(2); clf; hold on;
   l = plot(100*cksi, t1);
   set(l(1), 'color','b', 'linestyle','--');
   axis([0 25 0 800]);
   grid on;
   xlabel('Longitudinal creepage \xi [%]');
   ylabel('Max. temperature [^\circ{}C]');
   legend('Constant friction, \mu=0.33', ...
           'Temp.dep., \Delta\mu= 0.19, \Delta{}T=1200', ...
           'Temp.dep., \Delta\mu=-0.19, \Delta{}T=1200', ...
           'Temp.dep., \Delta\mu=-0.19, \Delta{}T= 400', 'location','northwest')

   if (print_fig)
      figure(1); print -djpeg95 figs/euro_tract_temp.jpg
      figure(2); print -djpeg95 figs/euro_max_temp.jpg
   end

   if (pause_after_plot), pause; end
   clear cksi creep force fstat fx opt;
end

clear old_graphics new_graphics pause_after_plot print_fig expnam;

