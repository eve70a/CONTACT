
show_fig = [ 4 5 ];
show_case = [ 10, 12, 14]; % -23, -25, -27 deg

sol = parse_out1('chalmers_flat_fz125.out');
s = loadcase('chalmers_flat_fz125', show_case(1));
sref = loadcase('chalmers_flat_fz125', 1);
sref.x_offset = sref.meta.xcp_r;
sref.y_offset = sref.meta.ycp_r;

mirror_y = 0; mirror_z = 0; scale_yz = 1;
prr = read_profile('../../examples/r300_wide.prr', [], mirror_y, mirror_z, scale_yz);
mirror_z = -1;
if (~exist('slcw'))
   slcw = read_profile('../../examples/S1002_flat.slcw', [], mirror_y, mirror_z, scale_yz);
end
clear mirror_y mirror_z scale_yz;

if (any(show_fig==1))
   opt = plot_2dspline;
   opt.urange    = [0 60]*pi/180;
   opt.typplot   = 'refl';
   opt.refl_avec = [0 1 0];
   opt.view      = [60 30];

   figure(1); clf;
   plot_2dspline(slcw, opt);
   hold on;
   th_w = -sol.ws_pos.pitch(show_case) + sol.cp_pos.xcp_w(show_case)/s.meta.rnom_whl;
   plot3(th_w, sol.cp_pos.ycp_w(show_case), sol.cp_pos.zcp_w(show_case), '.', 'markersize',12);
   if (slcw.u_method==1)
      title('Interpolating spline');
   else
      title('Approximating spline');
   end
   axis([0.2 0.7 -5 25 -5 5]);
end

if (any(show_fig==2))
   figure(2); clf; hold on;
   plot(sol.ws_pos.pitch*180/pi, sol.tot_forc.fz_tr/1e3, '-*');
   plot(sol.ws_pos.pitch(show_case)*180/pi, sol.tot_forc.fz_tr(show_case)/1e3, 'o');
   grid on;
   xlabel('\theta_{ws} [deg]');
   ylabel('F_z [kN]');
   title('Total force');
end

if (any(show_fig==3))
   figure(3); clf; hold on;
   plot(sol.ws_pos.pitch*180/pi, sol.cp_force.pmax, '-*');
   plot(sol.ws_pos.pitch(show_case)*180/pi, sol.cp_force.pmax(show_case), 'o');
   grid on;
   v = axis; v(3)=0; axis(v);
   xlabel('\theta_{ws} [deg]');
   ylabel('p_{max} [N/mm^2]');
   title('Max. pressure');
end

if (any(show_fig==4))
   opt = plot3d;
   opt.field = 'pn';
   opt.view  = 'rail';
   opt.addeldiv = 0;
   opt.addplot = 1;
   opt.zrange  = [0 2500];
 
   optel = opt;
   optel.field = 'eldiv_contour';
   optel.addplot = 1;

   iver = 2;
   if (iver==1)
      ncol = 6;
      nrow = 4;
      show_cases = [2 : 25];
   else
      ncol = 4;
      nrow = 2;
      show_cases = [4 : 2 : 18];
   end

   figure(4); clf;
   p = get(gcf,'position'); p(3) = 2.1*p(4); set(gcf,'position',p,'paperpositionmode','auto');
   for isub = 1 : length(show_cases)
      icase = show_cases(isub);
      s = loadcase('chalmers_flat_fz125', icase);
      s.x_offset = s.meta.xcp_r;
      s.y_offset = s.meta.ycp_r;

      irow = floor( (isub-1) / ncol ) + 1;
      icol = mod( (isub-1) , ncol ) + 1;
      subplot(nrow, ncol, isub);

      plot3d(s, opt);

      plot3d(sref, optel);
      c = get(gca,'children');
      delete(c(3:4));

      axis([-20 20 -15 20]);
      shading flat;
      % title(sprintf('$\\theta=%3d^\\circ$', round(sol.ws_pos.pitch(icase)*180/pi)), 'interpreter','latex');
      title('');
      if (icol==1), xlabel('x_r [mm]'); else, xlabel(''); end
      if (irow==nrow), ylabel('y_r [mm]'); else, ylabel(''); end
      text(-15, -13, sprintf('$\\theta_{ws}=%3d^\\circ$', round(sol.ws_pos.pitch(icase)*180/pi)), ...
                                                                                'interpreter','latex');
   end
end

if (any(show_fig==5))
   opt = plot3d;
   opt.field = 'pn';
   opt.view  = 'rail';
   opt.zrange  = [0 800];
   opt.addeldiv = 0;

   icase = 21;
   s = loadcase('chalmers_flat_fz125', icase);
   s.x_offset = s.meta.xcp_w;
   s.y_offset = s.meta.ycp_w;

   figure(5); clf;
   plot3d(s, opt);
   plot(s.meta.xcp_w, s.meta.ycp_w, '.', 'color',matlab_color(2), 'markersize',12);

   opt.field = 'eldiv_contour';
   opt.addplot = 1;
   plot3d(sref, opt);
   c = get(gca,'children');
   delete(c(3:4));

   % axis([-5 6 4 16]);
   title('');
   text(0.08, 0.92, sprintf('$\\theta=%3d^\\circ, F_z=%4.1f$ kN', round(sol.ws_pos.pitch(icase)*180/pi), ...
                sol.tot_forc.fz_tr(icase)/1e3), 'units','normalized', 'interpreter','latex');
   xlabel('x_w [mm]');
   ylabel('y_w [mm]');
end

if (any(show_fig==6))
   icase = 11;
   s   = loadcase('chalmers_flat_fz125', icase);

   opt = plot3d;
   opt.field    = 'pn';
   opt.typplot  = 'rw_rear';
   opt.rw_surfc = 'both';

   figure(6); clf;
   plot3d(s, opt, prr, slcw);
end

if (any(show_fig==7))
   ix = find(slcw.xsurf(:,1)*490>250, 1, 'first');
   iy = find(slcw.ysurf(ix,:)>=10, 1, 'last');

   du   = 0.05;
   u_ev = [0:du:40] * pi/180;
   y_ev = [8:1:12];
   [~,x_ev,~,z_ev] = eval_2dspline(slcw, u_ev, [], y_ev);
   ix = 1+round([23,25,27]/du);

   iver = 3;
   if (iver==1)
      fac = 180/pi;
   else
      fac = 490;
   end

   figure(7); clf; hold on;
   plot(x_ev*fac, z_ev);
   plot(x_ev(ix,3)*fac, z_ev(ix,3),'.', 'markersize',12, 'color',matlab_color(2));
   set(gca,'ydir','reverse');
   grid on;
   ylabel('z_{w} [mm]');
   legend({'$y_w=8$ mm', '$y_w=9$ mm', '$y_w=10$ mm', '$y_w=11$ mm', '$y_w=12$ mm'}, ...
                                'interpreter','latex', 'location','northwest');
   if (iver==1)
      axis([15 35 -0.35 -0.2]);
      xlabel('\theta_{wc} [deg]');
   elseif (iver==2)
      axis([130 300 -0.35 -0.2]);
      xlabel('r\theta_{wc} [mm]');
   else
      plot(slcw.xsurf(:,iy)*fac, slcw.zsurf(:,iy), ':.', 'markersize',9, 'linewidth',1, 'displayname','input data');
      axis([220 260 -0.32 -0.2]);
      xlabel('r\theta_{wc} [mm]');
      l = get(gca,'children'); delete(l([4,6]));
   end
   if (slcw.u_method==1)
      title('Interpolating spline');
   else
      title('Approximating spline');
   end
end
