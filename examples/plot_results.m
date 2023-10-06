
if (~exist('loadcase'))
   addpath('../matlab');
end
if (~exist('plot_arrow'))
   addpath('../matlab_intern');
end

expnam = strvcat('cattaneo', 'carter2d', 'bentall', 'visc_cylindr', 'catt_to_cart', ...
                 'subsurf', 'spence35', 'mbench', 'kpec', 'tractcurv', 'conformal', ...
                 'ertz_temperature', 'plastic_3bl', 'fastsim', 'veldep_fric');
% expnam = strvcat('kpec');
% expnam = strvcat('plastic_3bl');
expnam = strvcat('mbench');

pause_after_plot = 1 * (size(expnam,1)>1);
print_figures = 0;

% determine version information on the Matlab graphics system

old_graphics = verLessThan('matlab', '8.4.0');
new_graphics = ~old_graphics;

if (1==1);
   set(0,'defaultlinelinewidth',2);
   set(0,'defaultaxeslinewidth',2);
end
if (~verLessThan('matlab', '8.7.0'))
   set(0,'defaultlegendautoupdate', 'off');
end
set(0,'defaulttextfontsize',13);
set(0,'defaultaxesfontsize',13);
set(0,'defaultfigurepaperpositionmode','auto');

% 5.1 : The Cattaneo shift problem

if (~isempty(strmatch('cattaneo',expnam)))
   figure(1); clf; hold on;
   s = loadcase('cattaneo',2); % the tangential problem concerns the
   opt2 = plot2d;              % second case in the actual input-file
   opt2.yslc = 0.0;
   opt2.facpt = -1.0;
   plot2d(s, opt2);
   if (print_figures), print -djpeg75 cattaneo_px_y0.jpg; end

   if (pause_after_plot), pause; end
   clear s opt2;
end
   
% 5.2 : The 2D Carter/Fromm problem

if (~isempty(strmatch('carter2d',expnam)))
   figure(1); clf; hold on;
   s = loadcase('carter2d');
   opt2 = plot2d;
   opt2.yslc = 0.0;
   opt2.facpt = 1.0;
   plot2d(s, opt2);
   if (print_figures), print -djpeg75 carter2d_px.jpg; end

   if (pause_after_plot), pause; end
   clear s opt2;
end

% 5.3 : The 2D Bentall-Johnson testcase

if (~isempty(strmatch('bentall',expnam)))
   sol = parse_out3('bentall');
   s5 = loadcase('bentall',5);
   s6 = loadcase('bentall',6);

   for icase = 5:6
      ifig = icase - 4;
      figure(ifig); clf; hold on
      tmp = loadcase('bentall', icase);
      sx = sign(-tmp.px) .* tmp.srel;
      scl = 10000;
      l1 = plot(tmp.x, -tmp.px, 'b-');
      l2 = plot(tmp.x, scl*sx, 'r--');
      l3 = plot(tmp.x, tmp.fric.fstat*tmp.pn, 'b--');
      l4 = plot([-1 1], scl*sol.creep.cksi(icase)*[1 1], 'm-.');
      set([l1,l2,l3,l4],'linewidth',2);
      grid on
      if (icase==5)
         axis([-0.11 0.11 -6 6])
      else
         axis([-0.11 0.11 -4 8])
      end
      xlabel('x-coordinate [mm]')
      ylabel('Traction q(x) [N/mm^2]');
      if (icase==5)
         lg_pos = 'NorthEast';
      else
         lg_pos = 'SouthWest';
      end
      ix = find(tmp.eldiv(2:end)~=tmp.eldiv(1:end-1)) +1;
      for i = ix
         plot((tmp.x(i)-0.5*tmp.dx)*[1 1], [0, -mean(tmp.px(i-1:i))], 'k--');
      end
      legend([l1,l2,l4], 'Traction q(x)', 'Rel.slip s_x \cdot 10^{4}', ...
                                'Creepage \cdot 10^{4}', 'location', lg_pos);

      if (print_figures)
         if (icase==5)
            print -djpeg75 bentall_px_fxpos.jpg
         elseif (icase==6)
            print -djpeg75 bentall_px_fxneg.jpg
         end
      end
   end
   if (pause_after_plot), pause; end
   clear sol i icase ifig ix l1 l2 l3 l4 lg_pos s5 s6 scl sx tmp;
end

% 5.4 : The viscoelastic two-cylinder Carter/Fromm problem 

if (~isempty(strmatch('visc_cylindr',expnam)))
   figure(11);
   set(11,'position',[220 590 730 500]);
   set(11,'defaulttextfontsize',11);
   set(11,'defaultaxesfontsize',11);
   clf; hold on;

%  material parameters used:
   G571    = 571.5;
   G1143   = 1143;
   fstat   =  0.3;
%  reference values for scaling, from elastic case with G=571.5:
   ah_ref  = 6.34;
   pn_ref  = 60.3;
   px_ref  = pn_ref * fstat;
   
   sol_vt0   =loadcase('visc_cylindr',1); sol_vt0.x = sol_vt0.x / ah_ref;
   sol_vt9   =loadcase('visc_cylindr',2); sol_vt9.x = sol_vt9.x / ah_ref;
   sol_vt45  =loadcase('visc_cylindr',3); sol_vt45.x = sol_vt45.x / ah_ref;
   sol_G571  =loadcase('visc_cylindr',5); sol_G571.x = sol_G571.x / ah_ref;
   sol_G1143 =loadcase('visc_cylindr',6); sol_G1143.x = sol_G1143.x / ah_ref;

   ix = [1:3:sol_G571.mx];
   set(gca,'colororderindex',1);
   plot(sol_vt0.x,       -sol_vt0.px      /px_ref);
   plot(sol_vt9.x,       -sol_vt9.px      /px_ref);
   plot(sol_vt45.x,      -sol_vt45.px     /px_ref);
   plot(sol_G1143.x,     -sol_G1143.px    /px_ref);
   plot(sol_G571.x(ix),  -sol_G571.px(ix) /px_ref,'.','markersize',14);
   plot(sol_G1143.x(ix), -sol_G1143.px(ix)/px_ref,'v','markersize',6,'linewidth',1.5);
   set(gca,'colororderindex',1);
   plot(sol_vt0.x,       sol_vt0.pn  /pn_ref,'--');
   plot(sol_vt9.x,       sol_vt9.pn  /pn_ref,'-.');
   plot(sol_vt45.x,      sol_vt45.pn /pn_ref,'-.');
   plot(sol_G1143.x,     sol_G1143.pn/pn_ref,'--');
   c = get(gca,'children');
   axis equal
   axis([-2.0 1.2 -0.1 1.7]);
   set(gca,'xtick',[-1 -0.5 0 0.5 1])
   set(gca,'ytick',[ 0  0.5 1 1.5 2])
   grid on
   xlabel('Scaled coordinate x / a_h [-]')
   ylabel('Scaled p / p_0 and -q / \mu p_0 [-]')
   legend('viscoelastic, \tau_C=0 s', ...
          'viscoelastic, \tau_C=0.009 s', ...
          'viscoelastic, \tau_C=0.045 s', ...
          'viscoelastic, \tau_C=2 s', ...
          'elastic, G=571.5 MPa', ...
          'elastic, G=1143 MPa','Location','NorthWest')
   title('Tractions for varying relaxation distance')
   if (print_figures), print -djpeg75 visc_cylindr.jpg; end

   if (pause_after_plot), pause; end
   clear sol_vt0 sol_vt9 sol_vt45 sol_G571 sol_G1143 
   clear G571 G1143 ah_ref pn_ref px_ref fstat ix c
end

% 5.5 : Instationary problems - from Cattaneo to Carter

if (~isempty(strmatch('catt_to_cart',expnam)))
   figure(1); clf; hold on;
   r3=loadcase('catt_to_cart',13);  % "r3" = rolling, 3 units == case 13
   s3=loadcase('catt_to_cart',70);  % "s3" = shift, 3 units == case 70
   opt=plot2d; opt.yslc=0; opt.facpt=-1;
   plot2d(r3,opt);
   if (new_graphics), set(gca,'ColorOrderIndex',4); end
   l=plot(s3.x(1:end-12), -s3.px(17,13:end), '-*');  % row 17 == centerline y=0.
   if (old_graphics), set(l,'color','r'); end
   if (print_figures), print -djpeg75 catt_distc_12dx.jpg; end

   figure(2); clf; hold on;
   opt=plot3d; opt.field='eldiv';
   if (old_graphics), opt.colormap='cool'; end
   % if (new_graphics), opt.exterval=NaN; end
   plot3d(r3,opt);
   title('');
   if (print_figures), print -djpeg75 catt_eldiv_3units.jpg; end

   if (pause_after_plot), pause; end
   clear r3 s3 opt l;
end

% 5.6 : The calculation of subsurface stresses

if (~isempty(strmatch('subsurf',expnam)))
   figure(9); 
   set(9,'defaulttextfontsize',10);
   set(9,'defaultaxesfontsize',10);
   clf; hold on;
   [s1a, s1b] = loadstrs('subsurf',1);
   [s2a, s2b] = loadstrs('subsurf',2);
   dif = diffstrs(s2a, s1a);
   plot(s1a.z,  -squeeze(s1a.sighyd(1,1,:)), '-o');
   plot(s1a.z,   squeeze(s1a.sigvm(1,1,:)), '-*')
   l=plot(dif.z, squeeze(dif.sigvm(1,1,:)), '--*');
   if (old_graphics), set(l,'color','r'); end
   axis([0 2.5 0 1.8]); grid; set(gca,'ytick',[0:0.2:1.8]);
   legend('-\sigma_{hyd} for p_z=1, p_x=p_y=0', ...
          '\sigma_{vm} for p_z=1, p_x=p_y=0', ...
          '\sigma_{vm} for p_z=0, p_x=1, p_y=0');
   text(0.95, 0.60, 'single element, dx=dy=1', ...
                'horizontalalignment','right', 'units','normalized');
   text(0.95, 0.52, 'G=1, \nu=0.28', ...
                'horizontalalignment','right', 'units','normalized');
   xlabel('z [mm]')
   ylabel('stress [N/mm^2]');
   % set(gca,'position',[0.2 0.1 0.5 0.8])
   if (old_graphics), set(gcf,'paperposition',[0.25 2.5 6 8]); end
   if (new_graphics), set(gca,'dataaspectratio',[1.5 1 1]); end
   if (print_figures), print -djpeg75 subsurf_zaxis.jpg; end

   figure(2); clf; hold on;
   opt = plotstrs; opt.yslc = 0; opt.field = 'uz';
   opt.typplot = 'contourf';
   opt.cntrlvl = [0:0.02:0.08, 0.12:0.04:0.40];
   plotstrs(s2b, opt);
   % grid on;
   set(gca,'clim',[0 0.40]);
   h = findobj(gcf,'Type','colorbar');
   set(h,'ylim',[0 0.40], 'ytick',opt.cntrlvl);
   title('');
   if (print_figures), print -djpeg75 subsurf_uz.jpg; end

   if (pause_after_plot), pause; end
   clear c dif h ix j l mx s s1a s1b s2a s2b opt;
end

% 5.6* : Dissimilar materials: Spence compression

if (~isempty(strmatch('spence35',expnam)))
   for i=1:35, 
      eval(sprintf('s%d=loadcase(''spence35'',i);',i)); 
      eval(sprintf('str%d=loadstrs(''spence35'',i);',i));
   end

   % retrieve displacements ux along centerline y=0, both surface and subsurfc.
   ix=23;
   usrf=[];
   usub=[];
   for i=1:35, 
      eval(sprintf('usrf=[usrf, s%d.ux(:,ix)];',i)); 
      % 1st index of ux: depends on z-coordinates used in subsurf.calc.
      eval(sprintf('usub=[usub, squeeze(str%d.ux(1,:,:))];',i)); 
   end

   % plot element-division, direction of tangential tractions
   figure(1); clf;
   opt3=plot3d; opt3.field='ptvec'; opt3.numvecx=35; opt3.numvecy=35;
   plot3d(s35,opt3);
   axis([-3.7 3.7 -3.7 3.7])
   % title('Spence compression, 35 x 35 elements, final stage');
   title('');
   set(gca,'box','on')
   if (print_figures), print -djpeg75 spence_ptvec_case35.jpg; end

   % plot tangential tractions along radial line, traction bound, Hertzian

   fn=0.4705; a=3.5; b=3.5;
   pmax=1.5/pi*fn/(a*b);
   xhz=[0:0.05:a]; phz=pmax*sqrt(1-(xhz/a).^2);
   figure(2); clf; hold on;
   plot(xhz, s35.fric.fstat*phz, 'r--');
   plot(s35.x, s35.fric.fstat*s35.pn(ix,:), 'b--');
   plot(s35.x, s35.px(ix,:), 'b-o');
   grid on
   xlabel('X-axis [mm]');
   ylabel('Surface traction P_x / G [-]');
   legend('Hertzian traction bound','Traction bound f_{stat}*p_n','Tangential traction p_x')
   % title('Tangential tractions for the Spence problem');
   title('');
   axis([0 3.7 -0.001 0.008]);
   l=text(3.6,0.0049,'Stage 35, F_n=0.4705*G','horizontalalignment','right');
   fac = max(max(s35.pn))/pmax;
   l=text(3.6,0.0042,sprintf('Max p_n=%4.2f*p_{hz}',fac),'horizontalalignment','right');
   set(gca,'box','on')
   if (print_figures), print -djpeg75 spence_px.jpg; end

   if (pause_after_plot), pause; end
   for i=1:35, eval(sprintf('clear s%d str%d',i,i)); end
   clear i ix usrf usub opt3 fn a b pmax xhz phz l fac;
end

% 5.7 : Manchester contact benchmark

if (~isempty(strmatch('mbench',expnam)))
   % "case 11": 5 mm, 12.0 mrad
   l11          = loadcase('mbench_a22_left',11);
   [r11a, r11b] = loadcase('mbench_a22_right',11);
   % "case 13": 6 mm, 14.4 mrad
   l13 = loadcase('mbench_a22_left',13);
   r13 = loadcase('mbench_a22_right',13);
   prr = read_profile('Mbench_UIC60_v3.prr');
   prw = read_profile('Mbench_S1002_v3.prw');
   opt = plot3d;
   opt.rw_surfc = 'prr';
   opt.vecwidth = 0.5;
   opt.eldivwid = 1;

   figure(1); clf;
   plot3d(l11, opt, prr);
   set(gca,'cameraviewangle',8)
   set(gca,'cameratarget',[10 -6 25])
   title('')
   h=xlabel('x_{lft}');
   set(h,'position', [2, 37, 42])
   h=ylabel('y_{lft}');
   set(h,'position',[-12 -8 43])
   text(0, 0.95, 'Case A-2.2, 5 mm, 12.0 mrad', 'units', 'normalized')
   text(0, 0.87, 'Left wheel', 'units', 'normalized');

   if (print_figures), print -djpeg75 mbench_a22_5mm_lft_prr.jpg; end

   figure(2); clf;
   if (~isempty(r11b))
      opta = opt; opta.addplot = 1;
      plot3d(r11b, opt, prr);
      plot3d(r11a, opta, prr);
   else
      plot3d(r11a, opt, prr);
   end
   set(gca,'cameraviewangle',8)
   set(gca,'cameratarget',[10 -6 25])
   h=xlabel('x_{rgt}');
   set(h,'position', [2, -49 42])
   h=ylabel('y_{rgt}');
   set(h,'position',[-12 -4 43])
   title('')
   text(0, 0.95, 'Case A-2.2, 5 mm, 12.0 mrad', 'units', 'normalized')
   text(0, 0.87, 'Right wheel', 'units', 'normalized');

   if (print_figures), print -djpeg75 mbench_a22_5mm_rgt_prr.jpg; end

   figure(3); clf;
   plot3d(l13, opt, prr);
   set(gca,'cameraviewangle',8)
   set(gca,'cameratarget',[10 -6 25])
   view([60 30])
   title('')
   h=xlabel('x_{lft}');
   set(h,'position', [5, -33, 42])
   h=ylabel('y_{lft}');
   set(h,'position',[-12 12 43])
   text(0, 0.95, 'Case A-2.2, 6 mm, 14.4 mrad', 'units', 'normalized')
   text(0, 0.87, 'Left wheel', 'units', 'normalized');

   if (print_figures), print -djpeg75 mbench_a22_6mm_lft_prr.jpg; end

   figure(4); clf;
   plot3d(r13, opt, prr);
   set(gca,'cameraviewangle',8)
   set(gca,'cameratarget',[10 -6 25])
   h=xlabel('x_{rgt}');
   set(h,'position', [2, -49 42])
   h=ylabel('y_{rgt}');
   set(h,'position',[-12 -4 43])
   title('')
   text(0, 0.95, 'Case A-2.2, 6 mm, 14.4 mrad', 'units', 'normalized')
   text(0, 0.87, 'Right wheel', 'units', 'normalized');

   if (print_figures), print -djpeg75 mbench_a22_6mm_rgt_prr.jpg; end

   opt = plot3d; opt.field='pn'; opt.rw_surfc='both'; opt.typplot='rw_rear';

   figure(5); clf;
   plot3d(l11, opt, prr, prw);
   axis([-800 -720 -25 10]);
   text(0.05, 0.92, 'Case A-2.2, 5 mm, 12.0 mrad', 'units', 'normalized')
   text(0.05, 0.78, 'Left wheel', 'units', 'normalized');

   figure(6); clf;
   plot3d(r11a, opt, prr, prw);
   axis([ 700  780 -15 20]);
   text(0.05, 0.92, 'Case A-2.2, 5 mm, 12.0 mrad', 'units', 'normalized')
   text(0.05, 0.78, 'Right wheel', 'units', 'normalized');

   if (pause_after_plot), pause; end
   clear h l11 l13 r11a r11b r13 opt opta prr prw
end

% 5.7 : Manchester contact benchmark, KPEC + ANALYN

if (~isempty(strmatch('kpec',expnam)))

   prr  = read_profile('Mbench_UIC60_v3.prr');
   prw  = read_profile('Mbench_S1002_v3.prw');
   sol  = parse_out1('mbench_a22_kpec.out');

   for icase = 1 : 3
      kpec = loadcase('mbench_a22_kpec', 0+icase);
      anly = loadcase('mbench_a22_kpec', 3+icase);
      cntc = loadcase('mbench_a22_kpec', 6+icase);

      % plot contact patches for CONTACT, KPEC and ANALYN

      bbox = cntc_boundbox(cntc);
      bbox = cntc_boundbox(kpec, bbox);
      bbox = cntc_boundbox(anly, bbox);
      bbox = bbox + 2*[ -cntc.dx  cntc.dx  -cntc.dy  cntc.dy ];
      len   = bbox(2) - bbox(1);
      wid   = bbox(4) - bbox(3);
      if (icase==1)
         nx = 1; ny = 3; xc_lb = 2.3; yc_lb = -6;
      elseif (icase==2)
         nx = 1; ny = 3; xc_lb = 2.9; yc_lb =  5;
      elseif (icase==3)
         nx = 2; ny = 2; xc_lb = 3.0; yc_lb =  5;
      else
         xc_lb = []; yc_lb = [];
         if (len>0.5*wid)
            nx = 3; ny = 1;
         else
            nx = 1; ny = 3;
         end
      end

      opt = plot3d; 
      opt.field = 'eldiv';
      opt.view  = 'rail';
      opt.addplot = 1;

      figure(5+icase); clf; 
      subplot(ny,nx,1);
      plot3d(cntc, opt);
      axis equal;
      axis(bbox);
      title(''); ylabel(''); text(xc_lb, yc_lb, 'CONTACT', 'horizontalalignment','center');
      text(-xc_lb, yc_lb, sprintf('y_{ws}=%4.1f', sol.ws_pos.y(icase)), 'horizontalalignment','center');

      subplot(ny,nx,nx*ny-1);
      plot3d(kpec, opt);
      axis equal;
      axis(bbox);
      title(''); ylabel(''); text(xc_lb, yc_lb, 'KPEC', 'horizontalalignment','center');

      subplot(ny,nx,nx*ny);
      plot3d(anly, opt);
      axis equal;
      axis(bbox);
      title(''); ylabel(''); text(xc_lb, yc_lb, 'ANALYN', 'horizontalalignment','center');

      if (print_figures)
         print('-djpeg95', sprintf('mbench_kpec_el%d.jpg', icase));
      end
   end % icase

   if (pause_after_plot), pause; end
   clear sol cntc kpec anly prr prw opt icase bbox len wid nx ny  xc_lb yc_lb
end

% 5.8 : Calculation of traction curves, Eurosprinter test-case

if (~isempty(strmatch('tractcurv',expnam)))
   figure(8); clf;
   set(8,'defaulttextfontsize',10);
   set(8,'defaultaxesfontsize', 8);
   set(8,'defaultaxeslinewidth',1);
   hold on;

   fstat = [0.33, 0.33, 0.33, 0.36];   % [-]
   sol  = parse_out3('tractcurv.out');
   cksi = 100 * reshape(sol.creep.cksi, 30, 4);
   fx   =       reshape(sol.force.fx, 30, 4); fx = -fx * diag(fstat);

   dir = 'meldingen/afgesloten/m11014-slope-tractcurv/4-euro-final/';
   file = 'engel-data.mat';
   if (exist(['../../' dir file]))
      file =['../../' dir file];
   elseif (exist(['../../../' dir file]))
      file =['../../../' dir file];
   end
   if (exist(file))
      load(file);
      % l0 = plot(xmeas, ymeas, '.', 'color', [0 .75 .75], 'markersize', 12);
      l0 = plot(xmeas,ymeas,'*', 'color', [0 .7 .7], 'linewidth',1 );
   else
      l0 = [];
   end

   l1=plot(cksi(:,1), fx(:,1), '--', 'color', 'b');
   % l3=plot(cksi(:,3), fx(:,3), ':', 'color', matlab_color(4));
   l4=plot(cksi(:,4), fx(:,4), '-', 'color', matlab_color(3));
   l2=plot(cksi(:,2), fx(:,2), '-','color', matlab_color(2));
   axis([0 5 0 0.4501]); grid on;
   xlabel('Longitudinal creepage \xi [%]')
   ylabel('Traction coefficient F_x/F_n [-]');

   if (exist('l0') & ~isempty(l0))
      h=legend([l0,l1,l2,l4],'Measurement data', 'Original Kalker theory', ...
                 'Extended CONTACT model', 'Modified Fastsim algorithm', ...
                 'location','SouthEast');
   else
      h=legend([l1,l2,l4],'Original Kalker theory','Extended CONTACT model', ...
                 'Modified Fastsim algorithm', 'location','SouthEast');
   end
   p = get(h, 'position'); p(2) = p(2)-0.01; set(h,'position',p);

   l6  = text(4.8, 0.135, '36 km/h = 10 m/s','horizontalalignment','right');
   l7  = text(0.3, 0.135, 'Data from Engel, Beck & Alders (1998):');
   l8  = text(0.3, 0.110, '"Verschleissreduzierende Rad-schlupfregelung mit');
   l9  = text(0.3, 0.085, 'hoher Kraftschlussausnutzung",');
   l10 = text(0.3, 0.060, 'Elektrische Bahnen 96, pp.201-209');
   l11 = text(0.3, 0.035, 'Siemens Eurosprinter 127 001 locomotive');
   set([l6,l7,l8,l9,l10,l11], 'fontsize', 8);

   % print "dashed" lines for original and reduced slope:
   %  - estimate slope df / dxi at Fx/Fn = 0.1
   %  - set y-values for the dashes, compute x-values
   %  - plot dashes, plot arrow, plot text

   [~,ix1] = min(abs(fx(:,1)-0.1)); dfdxi1  = fx(ix1,1) / cksi(ix1,1);
   [~,ix2] = min(abs(fx(:,2)-0.1)); dfdxi2  = fx(ix2,2) / cksi(ix2,2);
   flin   = [0, 0.07, NaN, 0.12, 0.18, NaN, 0.23, 0.29, NaN, 0.34, 0.40];
   xilin1 = flin / dfdxi1; xilin2 = flin / dfdxi2;

   plot(xilin1, flin, 'k-');
   plot(xilin2, flin, 'k-');
   plot_arrow([0.23,0.37],[0.35,0],'k',1,0.5);
   text(0.2,0.42,'reduced initial slope');

   % print "lines" for original and falling friction:
   %  - locate points at suitable creepage - plot straight line segments

   xi_a = 3.0; xi_b = 4.0; df = 0.004;
   [~,ix1a] = min(abs(cksi(:,1)-xi_a)); [~,ix1b] = min(abs(cksi(:,1)-xi_b));
   [~,ix2a] = min(abs(cksi(:,2)-xi_a)); [~,ix2b] = min(abs(cksi(:,2)-xi_b));
   plot(cksi([ix1a,ix1b]), fx([ix1a,ix1b],1), 'k');
   plot(cksi([ix2a,ix2b]), fx([ix2a,ix2b],2), 'k');
   fmean = mean( [ fx([ix1a,ix1b],1) , fx([ix2a,ix2b],2) ] );
   plot_arrow([mean([xi_a,xi_b]), fmean(1)-df], ...
                   [0,fmean(2)-fmean(1)+2*df],'k', 1.5, 2);
   text(mean([xi_a,xi_b]), 0.35, 'falling friction effect', ...
                                             'horizontalalignment','center');
   if (print_figures), print -djpeg95 eurosprinter6.jpg; end
   if (pause_after_plot), pause; end

   c = get(gca,'children'); delete(c(1:10))
   axis([0 25 0 0.35]);
   set(l6 , 'position', [23.5, 0.115]);
   set(l7 , 'position', [1.5, 0.115]);
   set(l8 , 'position', [1.5, 0.085]);
   set(l9 , 'position', [1.5, 0.065]);
   set(l10, 'position', [1.5, 0.035]);
   set(l11, 'position', [1.5, 0.015]);

   if (print_figures), print -djpeg95 eurosprinter5.jpg; end
   if (pause_after_plot), pause; end
   clear c h l l0 l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 l11 dir file p
   clear fstat sol cksi fx xmeas ymeas flin fmean df
   clear ix1 ix1a ix1b ix2 ix2a ix2b xi_a xi_b xilin1 xilin2 dfdxi1 dfdxi2 
end

% 5.9 Conformal contact situation

if (~isempty(strmatch('conformal',expnam)))
   opt=plot3d;
   opt.exterval=NaN;
   opt.view=[90 -90];
   opt.field='ptabs+vec';
   opt.zrange=[0 320];
   opt.numvecx=10;
   opt.numvecy=10;
   opt.addplot=1;

   t1 = loadcase('conformal',1);
   t3 = loadcase('conformal',3);
   t4 = loadcase('conformal',4);
   t5 = loadcase('conformal',5);
   dif3 = diffcase(t5,t3);
   dif4 = diffcase(t5,t4);
   a = 8.6; a = 0.5*ceil(a/0.5+1);
   b = 6.6; b = 0.5*ceil(b/0.5+1);
   disp(sprintf('Max pressures: Hertzian %5.1f, full conformal %5.1f',...
                max(max(t1.pn)), max(max(t5.pn))));

   figure(1); clf;
   subplot(1,2,1);
   plot3d(t1,opt); axis equal; axis([-a a -b b]);
   title(''); ylabel(''); colormap jet;
   xlabel('Rolling X-coordinate [mm]');
   text(a-0.5, -b+0.5, 'Hertzian', 'fontsize',15);

   subplot(1,2,2);
   plot3d(t5,opt); axis equal; axis([-a a -b b]);
   xlabel(''); title(''); colormap jet;
   l=ylabel('Lateral S-coordinate [mm]');
   p=get(l,'position'); p(2)=-6; set(l,'position',p);
   text(a-0.5, -b+0.5, 'Conformal', 'fontsize',15);

   c=colorbar;
   h=get(gcf,'children');
   set(h(3),'position',[0.09 0.13 0.36 0.80]);
   set(h(2),'position',[0.45 0.13 0.36 0.80]);
   set(h(1),'position',[0.82 0.24 0.03 0.58]);
   set(h(2),'xticklabel',[]);
   ylabel(c, 'Traction [N/mm^2]');

   if (print_figures), print -djpeg75 conform_ptabs1.jpg; end

   opt.zrange=[0 130];
   opt.vecscale=0.005;
   opt.veccolor='m';
   figure(2); clf;
   subplot(1,2,1);
   plot3d(dif3,opt); axis equal; axis([-a a -b b]);
   title(''); ylabel('')
   xlabel('Rolling X-coordinate [mm]');
   text(a-0.5, -b+0.5, 'Planar creepage', 'fontsize',15);

   subplot(1,2,2);
   plot3d(dif4,opt); axis equal; axis([-a a -b b]);
   xlabel(''); title('');
   l=ylabel('Lateral S-coordinate [mm]');
   p=get(l,'position'); p(2)=-6; set(l,'position',p);
   text(a-0.5, -b+0.5, 'Halfspace IF', 'fontsize',15);

   c=colorbar;
   h=get(gcf,'children');
   set(h(3),'position',[0.09 0.13 0.36 0.80]);
   set(h(2),'position',[0.45 0.13 0.36 0.80]);
   set(h(1),'position',[0.82 0.24 0.03 0.58]);
   set(h(2),'xticklabel',[]);
   ylabel(c, 'Traction [N/mm^2]');

   if (print_figures), print -djpeg75 conform_ptabs2.jpg; end
   if (pause_after_plot), pause; end
   clear b t1 t3 t4 t5 dif3 dif4 opt c h l p;
end

% 5.10: Temperature calculation

if (~isempty(strmatch('ertz_temperature',expnam)))
   for i=1:4, 
      eval(sprintf('s%d=loadcase(''ertz_temperature'',i);',i)); 
   end

   % plot surface temperature due to frictional heating, bktemp=0

   opt = plot3d;
   opt.field    = 'temp1';
   opt.typplot  = 'contourf';
   opt.exterval = [];
   opt.ixrange  = 'all';
   opt.view     = [0 90];
   m = colormap('thermal'); m = m(49:end,:);
   opt.eldivcol(2,:) = m(160,:);
   a = 5.88; 

   figure(1); clf; hold on
   plot3d(s1, opt);
   shading flat;
   title('case 1: sliding at 1 m/s');
   colormap(m);
   set(gca,'clim',[0 200]);

   if (print_figures), print -djpeg95 temperature_fld.jpg; end

   % plot temperatures Tr,Tw along centerline y=0 for conduction

   iy=23;
   figure(2); clf; 
   subplot(2,1,1); hold on;
   l = plot(s1.x/a, [s1.temp2(iy,:); s1.temp1(iy,:)]);
   set(l(2), 'linestyle','--');
   axis([-1.5 6 0 175]); grid on;
   legend('T^{(2)} (wheel)', 'T^{(1)} (rail)');
   title('case 1: sliding at 1 m/s');

   subplot(2,1,2); hold on;
   plot(s2.x/a, [s2.temp2(iy,:); s2.temp1(iy,:)]);
   axis([-1.5 6 0 155]); grid on;
   xlabel('Scaled coordinate x/a [-]');
   h=ylabel('Temperature (increase) [{}^\circ{}C]');
   p=get(h,'position'); p(2)=230; set(h,'position',p);
   legend('T^{(2)} (wheel)', 'T^{(1)} (rail)', 'location','east');
   title('case 2: different bulk temperatures');

   if (print_figures), print -djpeg95 temperature_tser.jpg; end
   if (pause_after_plot), pause; end
   clear a h i iy l m opt p s1 s2 s3 s4
end

% 5.11: Third body layer cf. Hou et al

if (~isempty(strmatch('plastic_3bl',expnam)))
   sol   = parse_out3('plastic_3bl');
   cksi  = -sol.creep.cksi;
   fstat =  0.70;
   fx    =  fstat * sol.force.fn .* sol.force.fx;

   nser = 4;
   ncase = [13, 17, 13, 13];
   shft = zeros(20,nser);
   px   = zeros(20,nser);
   iofs = 1;
   for iser = 1 : nser
      irg = [1:ncase(iser)];
      shft(1+irg,iser) = cumsum(cksi(iofs+irg));
      shft(1+irg(end)+1:20,iser) = NaN;
      px(1+irg,iser) = fx(iofs+irg);
      iofs = iofs + ncase(iser);
   end

   figure(1); clf; hold on;
   plot(1000*shft, px, '.-', 'markersize',12);
   xlabel('applied displacement [\mu m]');
   ylabel('shear stress \tau [N/mm^2]');
   axis([0 600 0 620]); grid on
   legend('"magnetite"','"clay"','"sand"','"MoS2"', 'location', 'southeast');
   text(110, 255, '\rightarrow', 'rotation', 18);
   text(163, 155, '\rightarrow', 'rotation', -108);
   text(188, 140, '\rightarrow', 'rotation', 72);
   text(260, 315, '\rightarrow', 'rotation', 18);
   plot(70, 560, 'k.', 'markersize', 12);
   text(75, 520, '$(u_{c0}, \tau_{c0})$', 'interpreter', 'latex');
   plot([460 535 535], [368 368 398], 'k-', 'linewidth', 1);
   text(545, 360, '$k_u$', 'interpreter', 'latex');

   if (print_figures), print -djpeg95 hou1997_example.jpg; end
   if (pause_after_plot), pause; end
   clear sol cksi fstat fx iofs irg iser ncase nser px shft
end

% 5.12: The use of the FASTSIM algorithm

if (~isempty(strmatch('fastsim',expnam)))
   s1 = loadcase('fastsim'); 
   s4 = loadcase('fastsim',4); 
   s5 = loadcase('fastsim',5); 
   s6 = loadcase('fastsim',6); 
   opt = plot3d;
   if (1==1)
      opt.field='ptabs+vec'; opt.zrange=[0 400];
      opt.numvecx = 22; opt.numvecy = 12;
      f = 2e5;
      nam_plot = 'eldiv';
   else
      opt.field='sabs+vec'; opt.zrange=[0 0.4]; opt.view=[0 90];
      opt.addeldiv=1; opt.eldivcol=[0.75 0 0.75; 0.5 0 0];
      opt.veccolor='r';
      f = 0.4;
      nam_plot = 'slip';
   end
   opt.exterval=NaN;

   figure(1); clf; hold on;
   plot3d(s1, opt); 
   shading flat
   axis equal
   axis([-8.8 8.8 -4.4 4.4]);
   l1=plot3([6.0, 7.2],[-1.6 -2.4],[f f],'m'); t1=text(6.8,-2.8,'Adhesion');
   l2=plot3([4.4, 5.2],[-3.2 -3.6],[f f],'m'); t2=text(5.4,-3.8,'Slip');
   set([l1,l2],'linewidth',2); set([t1,t2],'fontsize',14);
   title('');
   if (print_figures), print('-djpeg75',['fastsim_',nam_plot,'_cntc.jpg']); end

   figure(2); clf; hold on;
   plot3d(s4, opt); 
   shading flat
   axis equal
   axis([-8.8 8.8 -4.4 4.4]);
   l1=plot3([-7.2,-5.6],[-2.0 -0.8],[f f],'m'); t1=text(-8.4,-2.4,'Adhesion');
   l2=plot3([-6.0,-4.8],[-3.2 -2.4],[f f],'m'); t2=text(-6.8,-3.6,'Slip');
   l3=plot3([ 6.0, 7.2],[-1.6 -2.4],[f f],'m'); t3=text( 6.8,-2.8 ,'Adhesion');
   set([l1,l2,l3],'linewidth',2); set([t1,t2,t3],'fontsize',14);
   title('');
   if (print_figures), print('-djpeg75',['fastsim_',nam_plot,'_ellip.jpg']); end

   figure(3); clf; hold on;
   plot3d(s5, opt); 
   shading flat
   axis equal
   axis([-8.8 8.8 -4.4 4.4]);
   l1=plot3([-6.0,-4.8],[-3.2 -2.4],[f f],'m'); t1=text(-6.8,-3.6,'Slip');
   l2=plot3([ 6.0, 7.2],[-1.6 -2.4],[f f],'m'); t2=text( 6.8,-2.8 ,'Adhesion');
   set([l1,l2],'linewidth',2); set([t1,t2],'fontsize',14);
   title('');
   if (print_figures), print('-djpeg75',['fastsim_',nam_plot,'_parab.jpg']); end

   figure(4); clf; hold on;
   plot3d(s6, opt); 
   shading flat
   axis equal
   axis([-8.8 8.8 -4.4 4.4]);
   l1=plot3([ 6.0, 4.8],[-3.2 -2.4],[f f],'m'); t1=text( 6.0,-3.6,'Slip');
   l2=plot3([-6.0,-7.2],[-1.6 -3.0],[f f],'m'); t2=text(-8.6,-3.6 ,'Adhesion');
   set([l1,l2],'linewidth',2); set([t1,t2],'fontsize',14);
   title('');
   if (print_figures), print('-djpeg75',['fastsim_',nam_plot,'_revers.jpg']); end

   if (pause_after_plot), pause; end
   clear nam_plot s1 s4 s5 s6 opt c l1 l2 l3 t1 t2 t3 f;
end

% 5.12* : Velocity-dependent friction

if (~isempty(strmatch('veldep_fric',expnam)))
   for t=[1,2,4,52], 
      eval(sprintf('s%d=loadcase(''veldep_fric'',t);',t));
   end
   g = 82000;

   figure(1); clf; hold on
   % get point numbers for x=-4, x=4, select points for markers
   ix=find(abs(4-abs(s1.x))<0.001); 
   dp=(ix(2)-ix(1))/16; rg=round([ix(1):dp:ix(2)]);
   % plot traction bound for s1
   plot(s1.x, s1.trcbnd/g, '--')
   % plot tractions p, slip s
   if (new_graphics), set(gca,'ColorOrderIndex',3); end
   plot(s1.x, [s1.px; s2.px; s4.px]/g, '-');
   if (new_graphics), set(gca,'ColorOrderIndex',3); end
   plot(s1.x, [s1.srel; s2.srel; s4.srel]*s1.kincns.veloc/1e6, '-.');
   if (new_graphics), set(gca,'ColorOrderIndex',3); end
   plot(s1.x(rg), s1.px(rg)/g, 'o');
   plot(s1.x(rg), s2.px(rg)/g, '*');
   plot(s1.x(rg), s4.px(rg)/g, 'v');
   if (new_graphics), set(gca,'ColorOrderIndex',3); end
   plot(s1.x(rg), s1.srel(rg)*s1.kincns.veloc/1e6, 'o');
   plot(s1.x(rg), s2.srel(rg)*s1.kincns.veloc/1e6, '*');
   plot(s1.x(rg), s4.srel(rg)*s1.kincns.veloc/1e6, 'v');
   if (old_graphics)
      l = flipud(get(gca,'children'));
      set(l([3,6,9,12]),'color',[0 .5 0]);
      set(l([4,7,10,13]),'color','r');
   end
   % set labels, legend
   axis([-4.5 4.5 0 0.0020]); grid on
   set(gca,'ytick',0:0.0005:0.0020)
   xlabel('{\it x}-coordinate [mm]');
   ylabel('Tractions {\it p/G} [-], Slip |{\it s}| [km/s]')
   title(['Grid \delta{}x=',num2str(s1.dx),' mm, V=45 m/s']);
   [legh,objh,outh,outm] = legend('{\it g}, case 1', '{\it p}, case 1','{\it p}, case 2','{\it p}, case 3','{\it s}, case 1','{\it s}, case 2','{\it s}, case 3');
   set(objh([11,17]),'marker','o');
   set(objh([13,19]),'marker','*');
   set(objh([15,21]),'marker','v');

   if (print_figures)
      l=num2str(s1.fric.frclaw); grd=num2str(s1.dx); grd=[grd(1),grd(3:end)];
      eval(['print -djpeg75 veldep_l',l,'_dx',grd,'_mult.jpg']);
   end

   figure(2); clf; hold on
   % plot traction bound for s52
   plot(s52.x, s52.trcbnd/g, '--')
   % plot tractions p, slip s
   plot(s52.x, s52.px/g, '-');
   l=plot(s52.x, s52.srel*s52.kincns.veloc/1e6, '-.');
   if (old_graphics), set(l,'color','r'); end
   % set labels, legend
   axis([-4.5 4.5 0 0.0020]); grid on
   set(gca,'ytick',0:0.0005:0.0020)
   xlabel('x-coordinate [mm]');
   ylabel('Tractions p/G [-], Slip |s_a| [km/s]')
   text( 1.4, 0.00135, '\mu p_n(x)');
   text(-2.9, 0.0009, 'p_x(x)');
   text(-3.6, 0.00015, 's_a(x)');
   title(['Grid \delta{}x=',num2str(s1.dx),' mm, V=45 m/s']);

   if (print_figures)
      l=num2str(s1.fric.frclaw); grd=num2str(s1.dx); grd=[grd(1),grd(3:end)];
      t=num2str(52-3);
      eval(['print -djpeg75 veldep_l',l,'_dx',grd,'_t',t,'.jpg']);
   end

   if (pause_after_plot), pause; end
   clear t s1 s2 s4 s52 ix dp legh objh outh outm l g t grd rg;
end

clear old_graphics new_graphics pause_after_plot print_figures expnam;

