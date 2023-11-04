
print_fig = 1;
show_figs = [  6:13 ];
% set(0,'defaultlinelinewidth',2);
% set(0,'defaultaxeslinewidth',2);
% set(0,'defaulttextfontsize',15);
% set(0,'defaultaxesfontsize',15);

icase = 1;
s = loadcase('test_pictures', icase);

prr = read_profile('../../examples/MBench_UIC60_v3.prr');
prw = read_profile('../../examples/MBench_S1002_v3.prw');
nprr = length(prr.ProfileY);
nprw = length(prw.ProfileY);

% lower the wheel to get clear picture of interpen.region
 
s.meta.z_w  = s.meta.z_w + 1;

% set rail & wheel markers in track-coordinates

o_r   = [ s.meta.y_r ; s.meta.z_r ];
o_w   = [ s.meta.y_w ; s.meta.z_w ];

% compute rotation matrices

cs = cos(s.meta.roll_r); sn = sin(s.meta.roll_r); 
rot_r = [cs, -sn; sn, cs];             % between rr and tr

cs = cos(s.meta.roll_w); sn = sin(s.meta.roll_w);
rot_w = [cs, -sn; sn, cs];             % between ws/rw and tr

% rotate rail profile to track coordinates

prr_tr = o_r * ones(1,nprr) + ...
         rot_r * [ prr.ProfileY' ; prr.ProfileZ' ];

% rotate wheel profile to track coordinates

prw_tr = o_w * ones(1,nprw) + ...
         rot_w * [ prw.ProfileY' ; prw.ProfileZ' ];

% compute undeformed distance in interpen.region

y_ud  = 712:0.1:738;
zr_ud = interp1( prr_tr(1,:), prr_tr(2,:), y_ud);
zw_ud = interp1( prw_tr(1,:), prw_tr(2,:), y_ud);
z_ud  = min(0, zr_ud - zw_ud);
ix_ud = 1:15:length(y_ud); % y_ud(ix_ud) = 712:1.5:738;

% compute alpha at rail surface

tmp   = atan2( zr_ud(3:end)-zr_ud(1:end-2), y_ud(3:end)-y_ud(1:end-2) );
ar_ud = [tmp(1), tmp, tmp(end)];

% compute contact reference point

y_cref = sum(z_ud.*y_ud) / sum(z_ud) - 0.3;     % manual offset for presentation purpose
z_cref = interp1( prr_tr(1,:), prr_tr(2,:), y_cref);
o_cptr = [ y_cref ; z_cref ];

% compute contact reference angle

dy        = 0.1;
dz        = interp1( prr_tr(1,:), prr_tr(2,:), y_cref+dy) - z_cref;
deltcp_tr = atan2(dz, dy);

cs = cos(deltcp_tr); sn = sin(deltcp_tr); 
rot_cptr = [ cs, -sn ; sn, cs ];        % between cp and tr

% rotate rail profile to contact coordinates

prr_cp = rot_cptr' * (prr_tr - o_cptr*ones(1,nprr));
o_rcp = rot_cptr' * (o_r - o_cptr);
deltcp_r = deltcp_tr - s.meta.roll_r;

% rotate wheel profile to contact coordinates

prw_cp = rot_cptr' * (prw_tr - o_cptr*ones(1,nprw));
o_wcp = rot_cptr' * (o_w - o_cptr);
deltcp_w = deltcp_tr - s.meta.roll_w;

% compute interpenetration depth & contact point Q on wheel surface

pen = min(z_ud) * cos(deltcp_tr);
q_cptr = o_cptr + pen * [ sin(deltcp_tr) ; -cos(deltcp_tr) ];

% compute extent of interpenetration region

ix = find(zr_ud-zw_ud < 0, 1, 'first') - 1;
ysta = y_ud(ix); %  - 0.75;
zsta = interp1( prr_tr(1,:), prr_tr(2,:), ysta);

ix = find(zr_ud-zw_ud < 0, 1, 'last') + 1;
yend = y_ud(ix); %  + 0.75;
zend = interp1( prr_tr(1,:), prr_tr(2,:), yend);

% project [ysta;zsta] and [yend,zend] onto contact reference plane
% sp == planar (s,n)-coordinates

sp_sta_cp = rot_cptr' * ([ysta; zsta] - o_cptr); sp_sta_cp(2) = 0;
sp_end_cp = rot_cptr' * ([yend; zend] - o_cptr); sp_end_cp(2) = 0;
sp_sta_tr = o_cptr + rot_cptr * sp_sta_cp;
sp_end_tr = o_cptr + rot_cptr * sp_end_cp;

% determine sp-grid == planar (s,n)-coordinates

dsp = 1.0;
sp_cgrid = [ round(sp_sta_cp/dsp)-2 : round(sp_end_cp/dsp)+2 ] * dsp;
tmp = o_cptr + rot_cptr * [sp_cgrid; zeros(size(sp_cgrid))];
y_cgrid = tmp(1,:);
z_cgrid = tmp(2,:);

% interpolate rail and wheel to contact grid positions

nr_cgrid = interp1( prr_cp(1,:), prr_cp(2,:), sp_cgrid);
nw_cgrid = interp1( prw_cp(1,:), prw_cp(2,:), sp_cgrid);
nw_cgrid = max(nw_cgrid, nr_cgrid-0.15);

% compute profile sr- & sw-coordinates at contact reference

sr_cp = interp1( prr_tr(1,:), prr.ProfileS, o_cptr(1) );
sw_cp = interp1( prw_tr(1,:), prw.ProfileS, q_cptr(1) );

% determine sc-grid == curved (s,n)-coordinates

sc_sta = interp1( prr_tr(1,:), prr.ProfileS, ysta ) - sr_cp;
sc_end = interp1( prr_tr(1,:), prr.ProfileS, yend ) - sr_cp;

dsc = 0.5;
sc_grid = [ round(sc_sta/dsc)-10 : round(sc_end/dsc)+13 ] * dsc;
sr_grid = sr_cp + sc_grid;
sw_grid = sw_cp - sc_grid;

yr_sc = interp1( prr.ProfileS, prr_tr(1,:), sr_grid);
zr_sc = interp1( prr.ProfileS, prr_tr(2,:), sr_grid);
yw_sc = interp1( prw.ProfileS, prw_tr(1,:), sw_grid);
zw_sc = interp1( prw.ProfileS, prw_tr(2,:), sw_grid);
yc_sc = (yr_sc + yw_sc) / 2;
zc_sc = (zr_sc + zw_sc) / 2;

% determine origin for curved contact surface

[~,ix] = min(abs(sc_grid));
o_crvtr = [yc_sc(ix); zc_sc(ix)];

% determine normal vectors at sc_sta & sc_end

ix = find(sc_grid<sc_sta, 1, 'last');
sc_sta_tr = [yc_sc(ix) ; zc_sc(ix)];
nc_sta = [ diff(zc_sc(ix:ix+1)) ; -diff(yc_sc(ix:ix+1)) ];
nc_sta = nc_sta / norm(nc_sta);

ix = find(sc_grid<sc_end, 1, 'last');
sc_end_tr = [yc_sc(ix) ; zc_sc(ix)];
nc_end = [ diff(zc_sc(ix:ix+1)) ; -diff(yc_sc(ix:ix+1)) ];
nc_end = nc_end / norm(nc_end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figures 1--14: plot profiles, markers, vertical overlap, interpen. area in track coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ifig = 1 : 14
   if (any(show_figs==ifig))
      figure(ifig); clf; hold on;
      axis equal; grid on;
      set(gca,'ydir','reverse', 'ticklabelinterpreter','latex');
      xlabel('$y_{tr}$ [mm]', 'interpreter','latex');
      ylabel('$z_{tr}$ [mm]', 'interpreter','latex');
      if (ifig<=1)
         axis([690 810 -20 50]);
      else
         axis([710 740 -2.5 15]);
      end

      % plot rail profile + rail profile marker

      plot(prr_tr(1,:), prr_tr(2,:), '-', 'markersize', 12);

      if (ifig<=1)
         plot_cm(o_r, 2.2, 8, s.meta.roll_r*180/pi, 'k');
         text(o_r(1)+1, o_r(2)+6, '$\mathbf{m}_{r(tr)}$', 'interpreter','latex', ...
                                                        'rotation', -s.meta.roll_r*180/pi);
      end

      % plot wheel profile

      set(gca,'colororderindex',3);
      plot(prw_tr(1,:), prw_tr(2,:), '-', 'markersize', 12);

      if (ifig<=1)
         plot_cm(o_w, 2.2, 8, s.meta.roll_w*180/pi, 'k');
         text(o_w(1)+0, o_w(2)-2, '$\mathbf{m}_{w(tr)}$', 'interpreter','latex', ...
                                'rotation', -s.meta.roll_w*180/pi, 'verticalalignment','bottom');
      end

      % plot search grid

      if (ifig>=2 & ifig<=11)
         y_i = y_ud(ix_ud);
         z_i = zeros(size(y_i));
         plot(y_i, z_i, '+', 'color',matlab_color(2), 'linewidth',1);
         plot(y_i, z_i, 'k-', 'linewidth',1);
      end

      if (ifig>=2 & ifig<=3)
         text(y_i(1), -0.4, 'gap mesh $(x_{ij}, y_{ij})$', 'interpreter','latex', ...
                                                                        'verticalalignment','bottom');
      end

      % plot vertical search lines

      if (ifig>=3 & ifig<=5)
         plot([1;1]*y_i, [z_i; z_i+13], ':', 'color',matlab_color(2), 'linewidth',1);
      end
      if (ifig>=3 & ifig<=4)
         text(715.5, 3.5, '`rays''', 'interpreter','latex', 'verticalalignment','bottom');
      end

      % plot intersections at rail surface

      if (ifig>=4 & ifig<=5)
         for ix = ix_ud
            plot(y_ud(ix)+0.45*[-1 1], zr_ud(ix)*[1 1], '-', 'color',matlab_color(2));
         end
         if (ifig>=4 & ifig<=5)
            ix = ix_ud(6);
            text(y_ud(ix)-0.5, zr_ud(ix), 'rail $z^r_{ij}$', 'interpreter','latex', ...
                                'horizontalalignment','right', 'verticalalignment','bottom');
         end
      end

      % plot intersections at wheel surface

      if (ifig>=5 & ifig<=5)
         for ix = ix_ud
            plot(y_ud(ix)+0.45*[-1 1], zw_ud(ix)*[1 1], '-', 'color',matlab_color(4));
         end
         if (ifig>=5 & ifig<=5)
            ix = ix_ud(4);
            text(y_ud(ix)-0.5, zw_ud(ix), 'wheel $z^w_{ij}$', 'interpreter','latex', ...
                                'horizontalalignment','right', 'verticalalignment','bottom');
         end
      end

      % plot vertical undeformed distance

      if (ifig>=6 & ifig<=11)
         zw_tmp = max(zw_ud, zr_ud-0.2);
         plot([1;1]*y_ud(ix_ud), [zr_ud(ix_ud); zw_tmp(ix_ud)], 'color', matlab_color(2));
      end
      if (ifig>=6 & ifig<=7)
         ix = ix_ud(7);
         zm = (zr_ud(ix) + zw_tmp(ix))/2;
         plot(y_ud(ix)+[0,2.2], zm+[0,1], 'color', matlab_color(2), 'linewidth',1);
         text(y_ud(ix)+2.4, zm+1, 'gap $g_{ij}=z_{ij}^r-z_{ij}^w$', 'interpreter','latex');
      end

      % plot surface inclination alpha_r

      if (ifig>=6 & ifig<=9)
         for ix = ix_ud(5:13)
            dy = 0.45*cos(ar_ud(ix)); dz = 0.45*sin(ar_ud(ix));
            plot(y_ud(ix)+dy*[-1 1], zr_ud(ix)-0.1+dz*[-1 1], '-', 'color',matlab_color(2));
         end
      end
      if (ifig==6)
         ix = ix_ud(6);
         plot(y_ud(ix)+[-1.4 0], zr_ud(ix)*[1 1], 'k-', 'linewidth',1);
         text(y_ud(ix)-1.5, zr_ud(ix), '$\alpha^r_{ij}$', 'interpreter','latex', ...
                             'horizontalalignment','right', 'verticalalignment','bottom');
      end

      % plot bounds of interpenetration area on horizontal plane

      if (ifig>=7 & ifig<=10)
         plot([1;1]*ysta, [zsta;0], 'k--', 'linewidth',1);
         text(ysta, -1.2, '$y_{sta(tr)}$', 'interpreter','latex', 'horizontalalignment','center');

         plot([1;1]*yend, [zend;0], 'k--', 'linewidth',1);
         text(yend, -1.2, '$y_{end(tr)}$', 'interpreter','latex', 'horizontalalignment','center');
      end

      % plot average contact position

      if (ifig>=8 & ifig<=9)
         plot(o_cptr(1)*[1 1], [0, o_cptr(2)], 'k-', 'linewidth',1);
         text(o_cptr(1), -1.2, 'avg. $y_{ref}$', 'interpreter','latex', 'horizontalalignment','center');
      end
      if (ifig==8)
         plot(o_cptr(1)*[1 1], o_cptr(2)+[0 2], 'k-', 'linewidth',1);
      end
      if (ifig==9)
         plot(o_cptr(1)+[-2.5 0], o_cptr(2)*[1 1], 'k-', 'linewidth',1);
         text(o_cptr(1)-2.7, o_cptr(2), '$z_{ref}$', 'interpreter','latex', 'horizontalalignment','right', ...
                                                                'verticalalignment','bottom');
      end

      % plot reference angle

      if (ifig==9)
         y_ang = o_cptr(1) + 0 + [4, 0, 4]*cos(deltcp_tr);
         z_ang = o_cptr(2) + 4 + [0, 0, 4]*sin(deltcp_tr);
         plot(y_ang, z_ang, 'k-', 'linewidth',1);
         text(y_ang(3)+0.2, z_ang(2), 'avg. $\alpha_{ref}$', 'interpreter','latex', ...
                                        'horizontalalignment','left', 'verticalalignment','bottom');
      end

      % plot grid on tangent plane 

      if (ifig>=13)
         nvec = [-sin(deltcp_tr); cos(deltcp_tr)]; nvec = nvec / norm(nvec);
         dy = 0.3*nvec(1)*[-1 1];
         dz = 0.3*nvec(2)*[-1 1];
         for ix = 1 : length(y_cgrid)
            plot(y_cgrid(ix)+dy, z_cgrid(ix)+dz, 'k-', 'linewidth',1);
         end
      end
      if (ifig==13)
         ix = round(length(y_cgrid)/2);
         text(y_cgrid(ix)-0.5*nvec(1), z_cgrid(ix)-0.5*nvec(2), 'contact grid $(x_{ij},s_{ij})$', ...
                'interpreter','latex', 'rotation',-deltcp_tr*180/pi, ...
                'horizontalalignment','center', 'verticalalignment','bottom');
      end

      % plot tangent plane

      if (ifig>=10)
         y_pln = o_cptr(1) + 10*[-1, 1]*cos(deltcp_tr);
         z_pln = o_cptr(2) + 10*[-1, 1]*sin(deltcp_tr);
         plot(y_pln, z_pln, '-', 'color',matlab_color(5));
      end

      if (ifig>=10 & ifig<=12)
         nvec = [-sin(deltcp_tr); cos(deltcp_tr)]; nvec = nvec / norm(nvec);
         ix = round(length(y_cgrid)/2);
         text(y_cgrid(ix)-0.5*nvec(1), z_cgrid(ix)-0.5*nvec(2), 'tangent plane', ...
                'interpreter','latex', 'rotation',-deltcp_tr*180/pi, ...
                'horizontalalignment','center', 'verticalalignment','bottom');
      end

      % plot reference marker

      if (ifig>=20 & ifig<=12)
         plot_cm(o_cptr, 0.6, 2.0, deltcp_tr*180/pi, 'k');
         text(o_cptr(1)-1, o_cptr(2)-1.5, '$\mathbf{m}_{cp(tr)}$', 'interpreter','latex', ...
                        'rotation',-deltcp_tr*180/pi, 'horizontalalignment','center');
      end
      if (ifig>=10)
         plot(o_cptr(1), o_cptr(2), 'k.', 'markersize',15);
      end

      % plot bounds of interpenetration area on tangent plane

      if (ifig>=10 & ifig<=13)
         plot([ysta,sp_sta_tr(1)], [zsta,sp_sta_tr(2)], 'k--', 'linewidth',1);
         plot([yend,sp_end_tr(1)], [zend,sp_end_tr(2)], 'k--', 'linewidth',1);
      end
      if (ifig==12)
         text(sp_sta_tr(1)-1, sp_sta_tr(2)-1, '$s_{sta}$', 'interpreter','latex', ...
              'rotation',-deltcp_tr*180/pi, 'horizontalalignment','center');
         text(sp_end_tr(1)-1, sp_end_tr(2)-1, '$s_{end}$', 'interpreter','latex', ...
              'rotation',-deltcp_tr*180/pi, 'horizontalalignment','center');
      end

      % plot search lines in normal direction

      if (ifig==14)
         nvec = [-sin(deltcp_tr); cos(deltcp_tr)]; nvec = nvec / norm(nvec);
         dy = nvec(1)*[ 0.3 8];
         dz = nvec(2)*[ 0.3 8];
         for ix = 1 : length(y_cgrid)
            plot(y_cgrid(ix)+dy, z_cgrid(ix)+dz, 'r:', 'linewidth',1);
         end
      end

      if (print_fig)
         set(gcf,'paperpositionmode','auto');
         print('-djpeg95', sprintf('wr_vert_ud%d.jpg',ifig));
      end
   end % show_fig
end % for ifig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 21:23 conformal: plot curved reference plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ifig = 21 : 23
   if (any(show_figs==ifig))
      iver = ifig - 20;    % 1: with cm, 2: with [sn]-axes, 3: dotted cref

      figure(ifig); clf; hold on;
      axis equal; grid on;
      set(gca,'ydir','reverse', 'ticklabelinterpreter','latex');
      xlabel('$y_{tr}$ [mm]', 'interpreter','latex');
      ylabel('$z_{tr}$ [mm]', 'interpreter','latex');
      axis([710 740 -2.5 15]);

      % plot rail profile

      plot(prr_tr(1,:), prr_tr(2,:), '-');
      % plot(yr_sc, zr_sc, '.', 'markersize', 12);

      % plot wheel profile

      set(gca,'colororderindex',3);
      plot(prw_tr(1,:), prw_tr(2,:), '-');
      % plot(yw_sc, zw_sc, '.', 'markersize', 12);

      % plot curved contact reference plane

      set(gca,'colororderindex',5);
      if (iver<=2)
         plot(yc_sc, zc_sc, '-');
      else
         plot(yc_sc, zc_sc, ':');
      end
   
      if (iver==1)
         plot_cm(o_crvtr, 0.6, 2.0, deltcp_tr*180/pi, 'k');
         text(o_crvtr(1)-1, o_crvtr(2)-2, '$\mathbf{m}_{cp(tr)}$', 'interpreter','latex', ...
              'rotation',-deltcp_tr*180/pi, 'horizontalalignment','center');
      end

      % plot normal and tangential vectors at different places

      if (iver>=2)
         [~,iy1] = min(abs(yc_sc-721));
         [~,iy2] = min(abs(yc_sc-732));
         for iy = [iy1, iy2]
            dy = yc_sc(iy+1) - yc_sc(iy);
            dz = zc_sc(iy+1) - zc_sc(iy);
            tvec = [ dy; dz]; tvec = tvec / norm(tvec);
            nvec = [-dz; dy]; nvec = nvec / norm(nvec);
            delt = atan2(dz, dy) * 180/pi;
            l = plot_axes( [yc_sc(iy), zc_sc(iy)], 2, delt, -1, ['$s$';'$n$';'   '], 'k', 0, 0, ...
                                                                                [],[],[], 1.2, 1);
            set(l(5:8), 'linewidth',1);
         end
      end

      % plot start & end of interpenetration region

      plot( sc_sta_tr(1)+nc_sta(1)*[2 0], sc_sta_tr(2)+nc_sta(2)*[2 0], 'k--', 'linewidth',1);
      plot( sc_end_tr(1)+nc_end(1)*[2 0], sc_end_tr(2)+nc_end(2)*[2 0], 'k--', 'linewidth',1);

      delt_sta = atan2( nc_sta(1), -nc_sta(2) );
      delt_end = atan2( nc_end(1), -nc_end(2) );
      text(sc_sta_tr(1)+2*nc_sta(1)-1, sc_sta_tr(2)+2*nc_sta(2)-1, ...
           '$\tilde{s}_{sta}$', 'fontsize',17, 'interpreter','latex', ...
           'rotation', -delt_sta*180/pi, 'horizontalalignment','center');
      text(sc_end_tr(1)+2*nc_end(1)-1, sc_end_tr(2)+2*nc_end(2)-1, ...
           '$\tilde{s}_{end}$', 'fontsize',17, 'interpreter','latex', ...
           'rotation', -delt_end*180/pi, 'horizontalalignment','center');

      if (print_fig)
         set(gcf,'paperpositionmode','auto');
         print('-djpeg95', sprintf('wr_curved_ref%d.jpg',iver));
      end
   end % show_fig
end % ifig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 31-33: plot profiles in terms of contact coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ifig = 31 : 33
   if (any(show_figs==ifig))
      iver = ifig - 30;

      figure(ifig); clf; hold on;
      axis equal; grid on
      xlabel('$s_{cp}$ [mm]', 'interpreter','latex');
      ylabel('$n_{cp}$ [mm]', 'interpreter','latex');
      set(gca,'ydir','reverse', 'ticklabelinterpreter','latex');

      if (iver==1)
         axis([-40 80 -20 50]);
      else
         axis([-10 10 -3.667 8]);
      end

      % plot rail profile

      plot(prr_cp(1,:), prr_cp(2,:));
      plot_cm(o_rcp, 2.2, 8, -deltcp_r*180/pi, 'k');
      text(o_rcp(1)-1, o_rcp(2)+6, '$\mathbf{m}_{r(cp)}$', 'interpreter','latex', ...
                   'rotation', (deltcp_tr-s.meta.roll_r) *180/pi)

      % plot wheel profile

      plot(prw_cp(1,:), prw_cp(2,:), 'color',matlab_color(3));
      plot_cm(o_wcp, 2.2, 8, -deltcp_w*180/pi, 'k');
      text(o_wcp(1)+2, o_wcp(2)-7, '$\mathbf{m}_{w(cp)}$', 'interpreter','latex', ...
                   'rotation', (deltcp_tr-s.meta.roll_w) *180/pi)

      % plot contact plane

      y_pln = 0 + 25*cos(0)*[-1, 2];
      z_pln = 0 + 25*sin(0)*[-1, 2];
      plot(y_pln, z_pln, '-', 'color',matlab_color(5));

      % plot bounds of interpenetration area on tangent plane

      if (iver<=-10)
         plot(sp_sta_cp(1)*[1 1], [1,sp_sta_cp(2)], 'k--');
         plot(sp_end_cp(1)*[1 1], [1,sp_end_cp(2)], 'k--');
      end

      % plot vertical search lines

      if (iver==3)
         plot([1;1]*sp_grid, [-0.3;5]*ones(size(sp_grid)), ':', 'color',matlab_color(2), 'linewidth',1);
      end

      % plot contact grid

      if (iver>=2)
         plot(sp_grid, zeros(size(sp_grid)), 'k.', 'markersize', 9);
      end
      if (iver==2)
         text(-9, -0.3, 'contact grid', 'interpreter','latex', 'verticalalignment','bottom');
      end

      % plot step size ds

      if (iver==2)
         sp_lbl = mean(sp_grid(15:16));
         plot(sp_grid(15)*[1 1], [0,-0.8], 'k-', 'linewidth',1);
         plot(sp_grid(16)*[1 1], [0,-0.8], 'k-', 'linewidth',1);
         text(sp_lbl, -1, '$\delta{}s$', 'interpreter','latex', 'horizontalalignment','center', ...
                                   'verticalalignment','bottom')
      end

      % plot vertical undeformed distance

      if (iver>=3)
         plot([1;1]*sp_cgrid, [nr_cgrid;nw_cgrid], 'color',matlab_color(2));
      end
      if (iver>=3)
         ix = 7;
         n_i = [nr_cgrid(ix), nw_cgrid(ix)];
         plot(sp_grid(ix)+[0,2], mean(n_i)+[0,1.5], '-', 'color',matlab_color(2), 'linewidth',1);
         text(sp_grid(ix)+2.2, mean(n_i)+1.5, '$h_{ij} = n^r_{ij} - n^w_{ij}$', 'interpreter','latex');
      end

      % plot contact reference point

      if (iver==1)
         plot_cm([0;0], 2.2, 8, 0, 'k');
         text(0+1, 0-3, '$\mathbf{m}_{cp(cp)}$', 'interpreter','latex', 'horizontalalignment','center', ...
                              'verticalalignment','bottom');
      end
      if (iver==2)
         plot_cm([0;0], 0.4, 1.5, 0, 'k');
         text(0, 0-1, '$\mathbf{m}_{cp(cp)}$', 'interpreter','latex', 'horizontalalignment','center', ...
                        'verticalalignment','bottom');
      end
      if (iver>=3)
         plot(0, 0, 'k.', 'markersize',18);
      end

      if (print_fig)
         set(gcf,'paperpositionmode','auto');
         print('-djpeg95', sprintf('wr_norm_ud%d.jpg', iver));
      end
   end % show_fig
end % ifig

