
show_fig = [ 2 ];

% 1: UIC60 profile, canonical form

if (any(show_fig==1))

   is_wheel = 0; mirror_y = 0; mirror_z = -1; scale_yz = 1;
   prr = read_profile('MBench_UIC60_v3.prr', is_wheel, mirror_y, mirror_z, scale_yz);

   yofs = 20;
   zofs =  5;
   phi  = -1/20;

   rot = [cos(phi), -sin(phi); sin(phi), cos(phi)];
   tmp = rot * [ prr.ProfileY' ; prr.ProfileZ' ];
   y_r = yofs + tmp(1,:);
   z_r = zofs + tmp(2,:);

   figure(1); clf; hold on;
   plot(y_r, z_r);
   set(gca,'ydir','reverse');
   grid on
   axis equal
   plot_cm([0,0], 2, 10, 0, 'k');
   axis([-45 60 -15 50]);
   xlabel('y_r [mm]');
   ylabel('z_r [mm]');
   text(0.05, 0.9, '\leftarrow track center', 'units','normalized', 'horizontalalignment','left');
   text(0.95, 0.9, 'field side \rightarrow', 'units','normalized', 'horizontalalignment','right');
end

% 2: wheel in 'wrong' conventions on input
% 3: wheel with input conventions corrected

if (any(show_fig==2))

   is_wheel = 1; mirror_y = -1; mirror_z = -1; scale_yz = 1;
   prw = read_profile('Car7216_0001r.whl', is_wheel, mirror_y, mirror_z, scale_yz);

   % create profile with 'wrong' sense

   inch = 25.4;
   prw.ProfileY = -prw.ProfileY / inch;
   prw.ProfileZ = -prw.ProfileZ / inch;

   figure(2); clf; hold on;
   plot(prw.ProfileY, prw.ProfileZ);
   grid on
   axis equal
   axis([-5.5 0.5 -1.5 1.5]);
   xlabel('x [inch]');
   ylabel('y [inch]');

   % create profile with corrected conventions

   is_wheel = 1; mirror_y = -1; mirror_z = -1; scale_yz = 1;
   prw = read_profile('Car7216_0001r.whl', is_wheel, mirror_y, mirror_z, scale_yz);

   figure(3); clf; hold on;
   plot(prw.ProfileY, prw.ProfileZ);
   grid on
   axis equal
   set(gca,'ydir','reverse');
   plot_cm([0,0], 3, 14, 0, 'k');
   axis([-10 140 -35 35]);
   xlabel('y_w [mm]');
   ylabel('z_w [mm]');
   text(0.05, 0.9, '\leftarrow track center', 'units','normalized', 'horizontalalignment','left');
   text(0.95, 0.9, 'field side \rightarrow', 'units','normalized', 'horizontalalignment','right');
end

if (any(show_fig==5))

   % 5: Crossing profile, canonical form

   is_wheel = 0; mirror_y = 0; mirror_z = -1; scale_yz = 1000;
   prr = read_profile('uk_crossing/Crossing_48.txt', is_wheel, mirror_y, mirror_z, scale_yz);

   figure(5); clf; hold on;
   plot(prr.ProfileY, prr.ProfileZ);
   set(gca,'ydir','reverse');
   grid on
   axis equal
   plot_cm([0,0], 3, 12, 0, 'k');
   axis([-50 100 -20 60]);
   xlabel('y_r [mm]');
   ylabel('z_r [mm]');
   text(0.05, 0.9, '\leftarrow track center', 'units','normalized', 'horizontalalignment','left');
   text(0.95, 0.9, 'field side \rightarrow', 'units','normalized', 'horizontalalignment','right');
end
