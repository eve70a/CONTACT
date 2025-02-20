
if (~exist('loadcase'))
   addpath('../matlab');
end
if (~exist('plot_arrow'))
   addpath('../matlab_intern');
end

expnam = strvcat('carter_plast', 'catt2d_plast', 'catt3d_plast', 'mbench_a22_left', ...
                 'plastic_3bl', 'tractcurve_plast' );
% expnam = 'mbench_a22_left';

pause_after_plot = 1 * (size(expnam,1)>1);
print_figures = 0;

%  determine version information on the Matlab graphics system

old_graphics = verLessThan('matlab', '8.4.0');
new_graphics = ~old_graphics;
set(0,'defaultfigurepaperpositionmode','auto');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Carter_plast: memo 17.050, v2.01, Figs 10, 11: krel=0.980, 0.6, 0, -0.3 ?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if (~isempty(strmatch('carter_plast',expnam)))
   for icase = [ 2, 7, 9, 12]
      s = loadcase('carter_plast', icase);
      if (~isfield(s,'taucrt'))
         s.mater.tau_c0 = 60;
         s.taucrt = 60*ones(size(s.x));
         s.uplsx  = zeros(size(s.x));
      end
      ylim = repmat([ 6e-6, 1.4e-4, 6e-4, 7e-4], 1, 3);
      krel = repmat([ 0.98, 0.60, 0, -0.333], 1, 3);

      figure(icase); clf;
      subplot(2,2,[1 2]); hold on;
      plot(s.x, s.px, '-o');
      plot(s.x, s.mu.*s.pn, '--');
      plot(s.x, s.mater.tau_c0*ones(1,s.mx), '--');
      set(gca,'colororderindex',3);
      plot(s.x, s.taucrt, '.-', 'markersize',15);
      axis([-1.1 1.1 0 60]); grid on;
      c = get(gca,'children'); legend(c([4,1]), 'p_x', '\tau_c');
      % legend('p_x', '\mu p_n', '\tau_{c0}', '\tau_c');
      set(gca,'colororderindex',1);
      ii = find(s.eldiv==1);
      plot(s.x(ii), zeros(size(ii)), '.', 'markersize',15);
      ii = find(s.eldiv==2);
      plot(s.x(ii), zeros(size(ii)), '.', 'markersize',15);
      ii = find(s.eldiv==3);
      plot(s.x(ii), zeros(size(ii)), '.', 'markersize',15);
      title(sprintf('k_{rel} = %6.3f, k_{tau} = %3.1e', ...
                                            krel(icase), s.mater.k_tau))

      subplot(2,2,3); hold on;
      set(gca,'colororderindex',2);
      plot(s.x, s.uplsx, '-*');
      axis([-1.1 1.1 0 ylim(icase)]); grid on;
      legend('u_{pl,x}', 'location','northeast');

      subplot(2,2,4); hold on;
      set(gca,'colororderindex',2);
      plot(s.x, s.srel, '-o');
      axis([-1.1 1.1 0 14e-4]); grid on;
      legend('s', 'location','northeast');

      if (print_figures), 
         fname = sprintf('figs/carter_krel%03d_dx%02d.jpg', 1000*krel(icase), ...
                                100*s.dx);
         print('-djpeg95', fname);
      end
   end

   if (pause_after_plot), pause; end
   clear c fname icase ii krel s ylim;
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2D Cattaneo to Carter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if (~isempty(strmatch('catt2d_plast',expnam)))
   for icase = [ 1, 2, 5, 10, 25, 50]
      s = loadcase('catt2d_plast', icase);
      if (~isfield(s,'taucrt'))
         s.mater.tau_c0 = 60;
         s.taucrt = 60*ones(size(s.x));
         s.uplsx  = zeros(size(s.x));
      end

      figure(icase); clf;
      subplot(2,2,[1 2]); hold on;
      plot(s.x, s.px, '-o');
      plot(s.x, s.mu.*s.pn, '--');
      plot(s.x, s.mater.tau_c0*ones(1,s.mx), '--');
      set(gca,'colororderindex',3);
      plot(s.x, s.taucrt, '.-', 'markersize',15);
      axis([-1.1 1.4 0 60]); grid on;
      c = get(gca,'children'); legend(c([4,1]), 'p_x', '\tau_c');
      % legend('p_x', '\mu p_n', '\tau_{c0}', '\tau_c');
      set(gca,'colororderindex',1);
      ii = find(s.eldiv==1);
      plot(s.x(ii), zeros(size(ii)), '.', 'markersize',15);
      ii = find(s.eldiv==2);
      plot(s.x(ii), zeros(size(ii)), '.', 'markersize',15);
      ii = find(s.eldiv==3);
      plot(s.x(ii), zeros(size(ii)), '.', 'markersize',15);
      title(sprintf('t = %3d', icase))

      subplot(2,2,3); hold on;
      set(gca,'colororderindex',2);
      plot(s.x, s.uplsx, '-*');
      axis([-1.1 1.1 0 1.4e-4]); grid on;
      legend('u_{pl,x}', 'location','northeast');

      subplot(2,2,4); hold on;
      set(gca,'colororderindex',2);
      plot(s.x, s.srel, '-o');
      axis([-1.1 1.1 0 14e-4]); grid on;
      legend('s', 'location','northeast');

      if (print_figures), 
         fname = sprintf('figs/catt2d_t%03d.jpg', icase);
         print('-djpeg95', fname);
      end
   end

   if (pause_after_plot), pause; end
   clear c fname icase ii krel s;
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 3D Cattaneo to Carter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if (~isempty(strmatch('catt3d_plast',expnam)))
   s2  = loadcase('catt3d_plast', 2);
   s9  = loadcase('catt3d_plast', 10);
   s58 = loadcase('catt3d_plast', 59);
   opt = plot3d;
   opt.addplot = 1;

   [~,iy0] = min(abs(s2.y));

   figure(101); clf;    % fig number > max(case) used above
   p = get(gcf,'position'); p(2)=min(400,p(2)); p(3:4)=[800 600]; set(gcf,'position',p);

   subplot(2,3,4);
   plot(s2.x, -s2.px(iy0,:), '-', s2.x, s2.mu(iy0,:).*s2.pn(iy0,:), '--', ...
        s2.x, s2.taucrt(iy0,:), ':');
   axis([-4 4 0 0.0075]); grid on;
   ylabel('p_x [N/mm^2]');
   subplot(2,3,5);
   plot(s9.x, -s9.px(iy0,:), '-', s9.x, s9.mu(iy0,:).*s9.pn(iy0,:), '--', ...
        s9.x, s9.taucrt(iy0,:), ':');
   axis([-4 4 0 0.0075]); grid on;
   xlabel('x_c [mm]');
   subplot(2,3,6);
   plot(s58.x, -s58.px(iy0,:), '-', s58.x, s58.mu(iy0,:).*s58.pn(iy0,:), '--', ...
        s58.x, s58.taucrt(iy0,:), ':');
   axis([-4 4 0 0.0075]); grid on;

   subplot(2,3,1);
   plot3d(s2, opt);
   axis([-4 4 -4 4]);
   xlabel('');
   title('t = 0');
   subplot(2,3,2);
   plot3d(s9, opt);
   axis([-4 4 -4 4]);
   xlabel(''); ylabel('');
   title('t = 7');
   subplot(2,3,3);
   plot3d(s58, opt);
   axis([-4 4 -4 4]);
   xlabel(''); ylabel('');
   title('t = 56');

   colorbar;
   h = get(gcf,'children');
   p = [0.89 0.15 0.02 0.70;
        0.66 0.58 0.21 0.34;
        0.39 0.58 0.21 0.34;
        0.12 0.58 0.21 0.34;
        0.66 0.12 0.21 0.34;
        0.39 0.12 0.21 0.34;
        0.12 0.12 0.21 0.34];
   for i=1:7
      set(h(i), 'position', p(i,:));
   end

   if (print_figures), 
      print -djpeg95 figs/catt3d_tract.jpg
   end

   opt = plot3d;
   opt.field = 'eldiv';
   figure(2); clf;
   plot3d(s58, opt);

   if (print_figures), 
      print -djpeg95 figs/catt3d_eldiv.jpg
   end

   opt.field = 'taucrt';
   opt.eldivcol(1,:) = [ 0.44 0.75 0.50 ];
   figure(3); clf;
   plot3d(s58, opt);

   if (print_figures), 
      print -djpeg95 figs/catt3d_taucrt.jpg
   end

   opt.field = 'upls+vec';
   opt.eldivcol(1,:) = [ 0.44 0.75 0.50 ];
   figure(4); clf;
   plot3d(s58, opt);

   if (print_figures), 
      print -djpeg95 figs/catt3d_taucrt.jpg
   end

   if (pause_after_plot), pause; end
   clear h i iy0 opt p s2 s58 s9;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Module 1: Manchester benchmark
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if (~isempty(strmatch('mbench_a22_left',expnam)))
   s4 = loadcase('mbench_a22_left');
   opt = plot3d;
   opt.view = 'rail';
   opt.numvecx = 5;
   opt.vecwidth = 2;

   figure(102); clf;     % fig number > max(case) used above
   p = get(gcf,'position'); p(3)=2*p(4); set(gcf,'position',p);
   plot3d(s4, opt);
   v=get(gca,'dataaspectratio'); v(2)=v(1); set(gca,'dataaspectratio',v);

   if (print_figures), 
      print -djpeg95 figs/mbench_ptabs.jpg
   end

   [~,iy2] = min(abs(s4.y-2));
   [~,iy3] = min(abs(s4.y+2));
   figure(2); clf; hold on;
   plot(s4.x, -s4.px(iy2,:));
   plot(s4.x, s4.fric.fstat*s4.pn(iy2,:), '--');
   plot(s4.x, s4.taucrt(iy2,:), ':');
   plot(s4.x, -s4.px(iy3,:));
   plot(s4.x, s4.fric.fstat*s4.pn(iy3,:), '--');
   axis([-3 4 0 100]); grid on;
   xlabel('x_c [mm]');
   ylabel('traction [N/mm^2]');
   legend('-p_x,   y= 2', ' \mu{}p_n, y= 2', ' \tau_c,    y= 2', '-p_x,   y=-2', ' \mu{}p_n, y=-2');

   if (print_figures), 
      print -djpeg95 figs/mbench_slices.jpg
   end

   figure(103); clf;     % fig number > max(case) used above
   opt.field = 'eldiv';
   p = get(gcf,'position'); p(3)=2*p(4); set(gcf,'position',p);
   plot3d(s4, opt);


   if (pause_after_plot), pause; end
   clear iy2 iy3 s4 opt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Tractcurve_plast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if (~isempty(strmatch('tractcurve_plast',expnam)))
   [c, f] = parse_out3('tractcurve_plast', 0);
   cksi = reshape(c(:,2), 30, 8);
   fx = reshape(f(:,2), 30, 8);
   fstat = 0.33;
   tau_c0 = 150;

   figure(1); clf;
   plot(100*cksi(:,[1,3:8]), -fstat*fx(:,[1,3:8]), '.-', 'markersize',12);
   xlabel('Longitudinal creepage \xi [%]');
   ylabel('Traction coefficient F_x/F_n [-]');
   axis([0 8 0 0.35]); grid on
   legend('Elastic', 'k_{rel}= 0.1', 'k_{rel}= 0.03', 'k_{rel}= 0.01', ...
          'k_{rel}= 0.003', 'k_{rel}= 0', 'k_{rel}=-0.001', 'location','southeast');
   text(0.4, 0.025, sprintf('h^3 = 0, \\mu = %4.2f, \\tau_{c0} = %3.0f', fstat, tau_c0));

   if (print_figures), 
      print -djpeg95 figs/tractcurv_krel.jpg
   end

   if (pause_after_plot), pause; end
   clear c cksi f fstat fx
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clear old_graphics new_graphics pause_after_plot print_figures expnam;

