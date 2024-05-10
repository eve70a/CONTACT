
icase      = 20;
show_wheel = 0;
show_patch = 1;
show_fig   = [ 1 2 3 4 5 6 ];

% expnam = 'cross_brute'; slc_file = '../profiles/cross_nose.slcs'; scale_yz = 1;
% expnam = 'cross_locus'; slc_file = '../profiles/cross_nose.slcs'; scale_yz = 1;
% expnam = 'wing_brute'; slc_file = '../profiles/wing_rail.slcs'; scale_yz = 1;
% expnam = 'cross+wing'; slc_file = '../profiles/cross+wing.slcs'; scale_yz = 1;
% expnam = 'cw_interrupt'; slc_file = '../profiles/cross+wing_extd.slcs'; scale_yz = 1;
% expnam = 'mbench_brute'; slc_file = '../profiles/uk_crossing.slcs'; scale_yz = 1000;
% expnam = 'mbench_locus'; slc_file = '../profiles/uk_crossing.slcs'; scale_yz = 1000;
  expnam = 'mbench_intrup'; slc_file = '../profiles/uk_interrupt_v2.slcs'; scale_yz = 1000;
% expnam = 'two_patches'; slc_file = '../profiles/uk_interrupt_v2.slcs'; scale_yz = 1000;

if (~exist('slcs') | ~strcmp(slcs.slc_file, slc_file))
   mirror_y = -1; mirror_z = -1; idebug   =  1;
   slcs  = read_profile(slc_file, [], mirror_y, mirror_z, scale_yz, idebug);
   clear mirror_y mirror_z idebug;
end

prw = read_profile('\cmcc\contact-cmcc\examples\MBench_S1002_v3.prw');

sol   = parse_out1([expnam,'.out']);
ncase = length(sol.npatch);
icase = min(icase, ncase);
s     = loadcase(expnam, icase);
if (strcmp(expnam,'cross_brute') | strcmp(expnam,'cross_locus'))
   s.meta.roty = 0.1;
end
show_patch = min(length(s), show_patch);

if (abs(s(show_patch).meta.roll_r)>1e-10)
   disp('Ignoring rail roll...')
end

% 1: plot rear view

if (any(show_fig==1))
   opt = plot3d;
   opt.typplot  = 'rw_rear';
   opt.rw_surfc = 'prr';
   opt.field    = 'pn';
   opt.xysteps  = 20;
   opt.vecscale = 20 / 2000;
   if (show_wheel), opt.rw_surfc = 'both'; end

   figure(1); clf;
   plot3d(s(show_patch), opt, slcs, prw);
end

% 2: plot side view

if (any(show_fig==2))
   opt = plot3d;
   opt.typplot  = 'rw_side';
   opt.rw_surfc = 'prr';
   opt.field    = 'pn';
   opt.xrange   = [-150, 150];
   opt.xysteps  = 20;
   opt.vecscale = 50 / 5000;
   if (show_wheel), opt.rw_surfc = 'both'; end

   figure(2);
   plot3d(s(show_patch), opt, slcs, prw);
end

% 3: plot 3d view

if (any(show_fig==3))
   opt = plot3d;
   opt.typplot  = 'surf';
   opt.rw_surfc = 'prr';
   opt.field    = 'pn';
   opt.xrange   = [-100, 100];
   opt.xysteps  = 20;

   if (1==1 & strcmp(expnam,'wing_brute'))
      opt.xrange   = [-60, 60];
      opt.xysteps  = [30, 8];
      opt.view     = [85, 30];
   end
   if (show_wheel), opt.rw_surfc = 'both'; end

   figure(3);
   plot3d(s, opt, slcs, prw);
end

% 4: plot pressures

if (any(show_fig==4))
   opt = plot3d;
   opt.field = 'pn';
   opt.view  = 'rail';

   % transform (xc,yc) coordinates to (xtr,ytr)
   s(show_patch).x_offset = s(show_patch).meta.xcp_r;
   s(show_patch).y_offset = s(show_patch).meta.ycp_r + s(show_patch).meta.y_r;

   figure(4);
   plot3d(s(show_patch), opt);
   set(gca,'dataaspectratio',[2 1 1])

   xlabel('x_{tr} [mm]'); ylabel('y_{tr} [mm]');
   % axis([-8 3 778 780.3]);
end

% 5: plot load transfer

if (any(show_fig==5) & any(strcmp(expnam, {'mbench_brute','mbench_locus'})))
   figure(5); clf; hold on;
   ix0 = find(sol.npatch>1, 1, 'first');
   ix1 = find(sol.npatch>1, 1, 'last');
   plot(sol.ws_pos.y(1:ix1), sol.cp_force.fn(1,1:ix1)/1e3, '-o');
   plot(sol.ws_pos.y(ix0:end), [sol.cp_force.fn(2,ix0:ix1), sol.cp_force.fn(1,ix1+1:end)]/1e3, '-*');
   plot(sol.ws_pos.y, sol.npatch);
   axis([0.80 1.15 0 11]);
   grid on;
   xlabel('y_{ws} [mm]');
   ylabel('F_n [kN]');
   legend('F_n, wing rail', 'F_n, crossing nose', '#patches', 'location','north');
end

% 6: plot load across interruption

if (any(show_fig==6) & strcmp(expnam, 'mbench_intrup'))
   figure(6); clf; hold on;
   plot(sol.ws_pos.x, sol.tot_forc.fz_tr/1e3, '-o');
   for iu = 1 : slcs.nslc
      if (slcs.u(iu)>110 & slcs.u(iu)<130)
         ic = find(sol.ws_pos.x>slcs.u(iu), 1, 'first');
         plot(slcs.u(iu)*[1 1], sol.tot_forc.fz_tr(ic)/1e3+0.05*[-1 1], '--', 'color',matlab_color(2))
         if (slcs.slc_ib(iu,1)<max(slcs.slc_ib([iu-1,iu+1],1)) | ...
             slcs.slc_ib(iu,2)>min(slcs.slc_ib([iu-1,iu+1],2)) )
            text(slcs.u(iu), sol.tot_forc.fz_tr(ic)/1e3+0.06, 'interruption', ...
                                'horizontalalignment','center', 'verticalalignment','bottom');
         end
      end
   end
   axis([110 130 12.4 13.0]);
   grid on;
   xlabel('s_{ws} [mm]');
   ylabel('F_n [kN]');
   legend('total force F_z', 'profile slices', 'location','northwest');
end

clear expnam slc_file show_fig show_wheel icase ncase show_patch opt iu ix0 ix1
