
show_fig = [1 2 3 4];

s1 = loadcase('catt3d_plast',1);
s2 = loadcase('catt3d_plast',2);
s3 = loadcase('catt3d_plast',59);
dif = diffcase(s1, s2);
iy0 = 17;

if (any(show_fig==1))
   figure(1); clf; hold on;
   plot(s1.x, s1.px(iy0,:))
   plot(s2.x, s2.px(iy0,:))
   plot(s3.x, s3.px(iy0,:),'--')
   grid on;
   set(gca,'colororderindex',1);
   plot(s1.x, 0.001*s1.eldiv(iy0,:), '-*');
   plot(s2.x, 0.001*s2.eldiv(iy0,:), '--o');
   plot(s3.x, 0.001*s3.eldiv(iy0,:), '-.x');
   title('Px');
end

if (any(show_fig==2))
   figure(2); clf; hold on;
   plot(s1.x, s1.uplsx(iy0,:))
   plot(s2.x, s2.uplsx(iy0,:))
   plot(s3.x, s3.uplsx(iy0,:),'--')
   grid on;
   legend('SteadyGS','ConvexGS','Transient');
   title('Uplsx');
end

if (any(show_fig==3))
   figure(3); clf; hold on;
   plot(s1.x, s1.taucrt(iy0,:))
   plot(s2.x, s2.taucrt(iy0,:))
   plot(s3.x, s3.taucrt(iy0,:))
   grid on;
   title('Taucrt');
end

if (any(show_fig==4))
   figure(4);
   opt = plot3d;
   opt.field = 'taucrt';
   plot3d(dif, opt);
end
