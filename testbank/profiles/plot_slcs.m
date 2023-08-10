
show_fig = [1 2 3];

if (any(show_fig==1))

   is_wheel = 0; mirror_y = 0; mirror_z = -1; scale_yz = 1;
   % slcs = read_profile('cross_nose.slcs', is_wheel, mirror_y, mirror_z, scale_yz);
   % slcs = read_profile('cross+wing.slcs', is_wheel, mirror_y, mirror_z, scale_yz);
   % slcs = read_profile('cross+wing_extd.slcs', is_wheel, mirror_y, mirror_z, scale_yz);
   % slcs = read_profile('wing_rail.slcs', is_wheel, mirror_y, mirror_z, scale_yz);

   % 1: settings for cross_nose, cross+wing, cross+wing_extd, wing_rail

   opt = plot_2dspline;
   opt.urange     = [ ];
   opt.zoom       = [1 1 1 ];
   opt.show_slc   = [ 1 : 99];
   opt.show_feat  = [1:10];
   opt.slc_color  = [1:6]';
   opt.feat_color = [2:7]';

   figure(1); clf; hold on;
   plot_2dspline([], slcs, opt);
end

if (any(show_fig==2))

   is_wheel = 0; mirror_y = 0; mirror_z = -1; scale_yz = 1; scale_yz = 1000;
   slcs = read_profile('uk_crossing.slcs', is_wheel, mirror_y, mirror_z, scale_yz);

   % 2: settings for uk_crossing

   opt = plot_2dspline;
   opt.urange     = [-2000 4000 ];
   opt.zoom       = [10 1 1 ];
   opt.show_slc   = [ 1 : 99];
   opt.show_feat  = [1:10];
   opt.slc_color  = [1:6]';
   opt.feat_color = [2:7]';

   figure(2); clf; hold on;
   plot_2dspline([], slcs, opt);

end

if (any(show_fig==3))

   is_wheel = 0; mirror_y = 0; mirror_z = -1; scale_yz = 1; scale_yz = 1000;
   slcs = read_profile('uk_interrupt_v1.slcs', is_wheel, mirror_y, mirror_z, scale_yz);

   % 3: settings for uk_interrupt_v1

   opt = plot_2dspline;
   opt.urange     = [ ];
   opt.zoom       = [ 1 1 1 ];
   opt.show_slc   = [ 1:99];
   opt.show_feat  = [ 1:99 ];
   opt.slc_color  = [1:6]';
   opt.feat_color = [2:7]';
   opt.view       = [60 20];

   figure(3); clf; hold on;
   plot_2dspline([], slcs, opt);
   title(slcs.slc_file, 'interpreter','none');

end
