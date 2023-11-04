
show_fig = [ 3 : 6   ];
show_case = [ 2, 5, 8]; % -23, -25, -27 deg

r_nom   = 450;

sol = parse_out1('rounded_flat_d09.out');
s = loadcase('rounded_flat_d09', show_case(1));

mirror_y = 0; mirror_z = -1; scale_yz = 1;
prr = read_profile('../profiles/circ_r300.prr', [], mirror_y, mirror_z, scale_yz);
if (~exist('slcw'))
   slcw = read_profile('../profiles/flat_d09.slcw', [], mirror_y, mirror_z, scale_yz);
end
clear mirror_y mirror_z scale_yz;

if (any(show_fig==1))
   opt = plot_2dspline;
   % opt.urange   = [-0 60]*pi/180;
   opt.typplot   = 'refl';
   opt.refl_avec = [105 30];
   opt.refl_dth  = 5;
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
   % axis([0.2 0.7 -5 25 -5 5]);
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
   opt.zrange  = [0 700];

   ncol = 4;
   nrow = 2;
   iver = 3;

   figure(4); clf;
   p = get(gcf,'position'); p(3) = 2.5*p(4); set(gcf,'position',p);
   set(gcf, 'paperpositionmode','auto');
   for ipic = 1 : 8
      if (iver==1)
         icase = ipic;
      elseif (iver==2)
         icase = 9 - ipic;
      elseif (iver==3)
         icase = 16 - ipic;
      end
      s = loadcase('rounded_flat_d09', icase);
      s.x_offset = s.meta.xcp_w;
      s.y_offset = s.meta.ycp_w;
      irow = floor( (ipic-1) / ncol ) + 1;
      icol = mod( (ipic-1) , ncol ) + 1;
      isub = icol + (irow-1)*ncol;
      subplot(nrow, ncol, isub);
      plot3d(s, opt);
      axis([-15 15 -24 24]);
      shading flat;
      title(sprintf('$\\theta_{ws}=%3d^\\circ$', round(sol.ws_pos.pitch(icase)*180/pi)), ...
                                                                        'interpreter','latex');
      xlabel('x_w [mm]');
      ylabel('y_w [mm]');
   end
end

if (any(show_fig==5))
   opt = plot3d;
   opt.field = 'pn';
   opt.view  = 'rail';
   opt.zrange  = [0 800];
   opt.addeldiv = 0;

   icase =  8;
   s = loadcase('rounded_flat_d09', icase);
   s.x_offset = s.meta.xcp_w;
   s.y_offset = s.meta.ycp_w;

   figure(5); clf;
   plot3d(s, opt);
   shading flat;
   plot(s.meta.xcp_w, s.meta.ycp_w, '.', 'color',matlab_color(2), 'markersize',12);
   axis([-4 4 -25 25]);
   title('');
   text(0.08, 0.92, sprintf('$\\theta_{ws}=%3d^\\circ, F_z=%4.1f$ kN', round(sol.ws_pos.pitch(icase)*180/pi), ...
                sol.tot_forc.fz_tr(icase)/1e3), 'units','normalized', 'interpreter','latex');
   xlabel('x_w [mm]');
   ylabel('y_w [mm]');
end

if (any(show_fig==6))
   s   = loadcase('rounded_flat_d09', 8);

   opt = plot3d;
   opt.field    = 'pn';
   opt.typplot  = 'rw_rear';
   opt.rw_surfc = 'both';
   opt.vecscale = 20 / 1000;

   figure(6); clf;
   plot3d(s, opt, prr, slcw);
   axis([-50 50 -28 20]);
   text(-45, -25, sprintf('$\\theta_{ws}=%3d^\\circ, F_z=%4.1f$ kN', round(sol.ws_pos.pitch(icase)*180/pi), ...
                sol.tot_forc.fz_tr(icase)/1e3), 'interpreter','latex', 'verticalalignment','top');
   text(-45, 5, '$r_{crown}=300$ mm', 'interpreter','latex', 'verticalalignment','top');
end

if (any(show_fig==7))
   ix = find(slcw.xsurf(:,1)>=0, 1, 'first');
   iy = find(slcw.ysurf(ix,:)>=0, 1, 'last');

   du   = 0.001;
   u_ev = [-0.15:du:0.15];
   y_ev = [-15:5:5];
   [~,x_ev,~,z_ev] = eval_2dspline(slcw, u_ev, [], y_ev);
   ix = 1+round( ([-7:0]*pi/180+0.15) / du);

   iver = 3;
   if (iver==1)
      fac = 1;
   else
      fac = r_nom;
   end

   figure(7); clf; hold on;
   plot(x_ev*fac, z_ev);
   plot(x_ev(ix,4)*fac, z_ev(ix,4),'.', 'markersize',12, 'color',matlab_color(2));
   set(gca,'ydir','reverse');
   grid on;
   ylabel('z_{w} [mm]');
   legend({'$y_w=-15$ mm', '$y_w=-10$ mm', '$y_w=-5$ mm', '$y_w=0$ mm', '$y_w=5$ mm'}, ...
                                'interpreter','latex', 'location','northwest');
   if (iver==1)
      axis([-0.15 0.15 -1.0  0.1]);
      xlabel('\theta_{wc} [rad]');
   elseif (iver==2)
      axis([-60 60 -1.0  0.1]);
      xlabel('r\theta_{wc} [mm]');
   else
      plot(slcw.xsurf(:,iy)*fac, slcw.zsurf(:,iy), ':.', 'markersize',9, 'linewidth',1, 'displayname','input data');
      axis([-60 -35 -0.10 0.02]);
      xlabel('r\theta_{wc} [mm]');
      l = get(gca,'children'); delete(l([2,3]));
   end
   if (slcw.u_method==1)
      title('Interpolating spline');
   else
      title('Approximating spline');
   end
end
