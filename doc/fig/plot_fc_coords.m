
print_fig = 0;
use_3d = 0;

% centerline: straight [c0 - c1], circular [c1 - c3]
% left rail:  straight [l0 - l1], circular [l1 - l3]
% right rail: straight [r0 - r1], circular [r1 - r3]

wid  = 4;
th0  = 10;
th3  = 90;
c0   = [ 5; -35 ];
l0   = c0 - wid * [sin(th0*pi/180); -cos(th0*pi/180)];
r0   = c0 + wid * [sin(th0*pi/180); -cos(th0*pi/180)];

len1 = 30;
th1  = th0;
c1   = c0 + len1 * [cos(th1*pi/180);  sin(th1*pi/180) ];
l1   = c1 - wid  * [sin(th1*pi/180); -cos(th1*pi/180)];
r1   = c1 + wid  * [sin(th1*pi/180); -cos(th1*pi/180)];

r2   = 40;
th2  = [ th0:1:th3 ];
nth  = length(th2);
cm   = c1 - r2 * [ sin(th2(1)*pi/180); -cos(th2(1)*pi/180) ];
c2   = cm*ones(1,nth) + r2 * [ sin(th2*pi/180); -cos(th2*pi/180) ];
l2   = c2 - wid * [sin(th2*pi/180); -cos(th2*pi/180)];
r2   = c2 + wid * [sin(th2*pi/180); -cos(th2*pi/180)];

c3   = c2(:,end);

% introduce superelevation by roll angle phi

phi  =  10 * pi/180;
sphi =  sin(phi);
zl   = -wid * sphi;
zr   =  wid * sphi;

figure(1); clf; hold on;
axis equal;
if (use_3d)
   axis([-10 80 -50 20 -15 15]);
   view([-20 30]);
else
   axis([-15 80 -50 15]);
end

% plot track curve and rail curves

n = length(c2);
plot3([c0(1),c1(1),c2(1,:)], [c0(2),c1(2),c2(2,:)], zeros(1,n+2), '-', 'color',matlab_color(2));
hold on;
plot3([l0(1),l1(1),l2(1,:)], [l0(2),l1(2),l2(2,:)], [0,zl*ones(1,n+1)], 'k-', 'linewidth',1);
plot3([r0(1),r1(1),r2(1,:)], [r0(2),r1(2),r2(2,:)], [0,zr*ones(1,n+1)], 'k-', 'linewidth',1);

grid on;
xlabel('x_{isys} [m]');
ylabel('y_{isys} [m]');

% plot Isys axis system

if (use_3d)
   axnams = ['$x_{isys}$'; '$y_{isys}$'; '$z_{isys}$']; col = 1; th_txt = 0;
   rot    = [180,0,0];
   factxt = [1.4 1.4 1.1];
   l = plot_axes3([0,0,0], 10, rot, axnams, col, [],[],[], th_txt, factxt, 0 );
   set(l([4,7,10]), 'interpreter','latex');
else
   axnams = ['$x_{isys}$'; '$y_{isys}$'; '$z_{isys}$']; col = 1; th_txt = 0;
   l = plot_axes([0,0], 10, 0, -1, axnams, col, ...
                   [],[],[],[], th_txt, [1.1,1.5; 2,1; 3,2.0], 1 );
   set(l([6,9,10]), 'interpreter','latex');
end

% plot Fc axis system at s_fc=0

if (use_3d)
   axnams = ['$x_{fc}$'; '$y_{fc}$'; '$z_{fc}$']; col = 2; th_txt = 0;
   factxt = [1.3 1.3 1.2];
   rot = [180, -th0, 0];
   l = plot_axes3([c0;0], 10, rot, axnams, col, [],[],[], th_txt, factxt, 0 );
   set(l([4,7,10]), 'interpreter','latex');
   set(l(1),'color','k');
else
   axnams = ['$x_{fc}$'; '$y_{fc}$'; '$z_{fc}$']; col = 2; th_txt = 0;
   l = plot_axes(c0, 10, -th0, -1, axnams, col, ...
                   [],[],[],[], th_txt, [1.1,0.9; 1.5,1; 3,1.5], 1 );
   set(l([6,9,10]), 'interpreter','latex');
end

% plot Gamma(s)

cg = 0.7 * c1 + 0.3 * c0;
lg = cg - 10 * [sin(th1*pi/180); -cos(th1*pi/180)];
plot([cg(1),lg(1)], [cg(2),lg(2)], 'k-', 'linewidth',1);
text(lg(1), lg(2)+1, '$\boldmath{\Gamma}(s)$', 'interpreter','latex', ...
                'verticalalignment','bottom', 'horizontalalignment','center');

% plot arrow for track radius r_curv

ixa = find(th2==85);
l = plot_arrow(cm, c2(:,ixa)-cm, 'k', 0.25);
set(l, 'linewidth',1);
plot(cm(1), cm(2), 'k.', 'markersize',12);
cmm = (cm + c2(:,ixa))/2;
text(cmm(1), cmm(2)+1, '$r_{curv}$', 'interpreter','latex', ...
                'horizontalalignment','center', 'verticalalignment','bottom');

if (0==1)

   % plot two slices

   ix1  = find(th2==30);
   ix2  = find(th2==40);
   slc1 = c2(:,ix1)*[1 1 1 1 1] + ...
                   wid * [sin(th2(ix1)*pi/180); -cos(th2(ix1)*pi/180)] * [-2:2];
   slc2 = c2(:,ix2)*[1 1 1 1 1] + ...
                   wid * [sin(th2(ix2)*pi/180); -cos(th2(ix2)*pi/180)] * [-2:2];

   plot(slc1(1,:), slc1(2,:), '--', 'color',matlab_color(4));
   plot(slc2(1,:), slc2(2,:), '--', 'color',matlab_color(4));
   % plot(slc1(1,2:end-1), slc1(2,2:end-1), '.', 'markersize',9);
   % plot(slc2(1,2:end-1), slc2(2,2:end-1), '.', 'markersize',9);

   % plot arrows between two slices

   s1 = (slc1(:,1) + slc1(:,2)) / 2;
   s2 = (slc2(:,1) + slc2(:,2)) / 2;
   s3 = (slc1(:,4) + slc1(:,5)) / 2;
   s4 = (slc2(:,4) + slc2(:,5)) / 2;
   sl = (s1 + s2) / 2;
   sr = (s3 + s4) / 2;
   l1 = plot_arrow(sl, s1-sl, 'k', 1.8);
   l2 = plot_arrow(sl, s2-sl, 'k', 1.8);
   l3 = plot_arrow(sr, s3-sr, 'k', 1.4);
   l4 = plot_arrow(sr, s4-sr, 'k', 1.4);
   set([l1;l2;l3;l4], 'linewidth',1);
else

   % plot grid on fc-coordinate system

   for k = 1 : 7
      ix  = find(th2==12+3*k);
      slc_x(k,:) = c2(1,ix) + wid * sin(th2(ix)*pi/180) * [-2:2];
      slc_y(k,:) = c2(2,ix) - wid * cos(th2(ix)*pi/180) * [-2:2];
      slc_z(k,:) =          + wid * sphi * [-2:2];
   end
   plot3(slc_x, slc_y, slc_z, '-', 'linewidth',1, 'color',matlab_color(2));
   plot3(slc_x', slc_y', slc_z', '-', 'linewidth',1, 'color',matlab_color(2));
end

% plot cartesian tr grid as used in CONTACT

ixtr = find(th2==60);
slctr = [c2(:,ixtr);0] * [1 1 1 1 1] + ...
                wid * [sin(th2(ixtr)*pi/180); -cos(th2(ixtr)*pi/180); sphi] * [-2:2];

tx = cos(th2(ixtr)*pi/180);
ty = sin(th2(ixtr)*pi/180);
for k = -3 : 3
   plot3(slctr(1,:)+2*k*tx, slctr(2,:)+2*k*ty, slctr(3,:), '-', ...
                'linewidth',1, 'color',matlab_color(5));
end
for i = 1 : 5
   plot3(slctr(1,i)+2*[-3,3]*tx, slctr(2,i)+2*[-3,3]*ty, slctr(3,i)+[0 0], '-', ...
                'linewidth',1, 'color',matlab_color(5));
end

% plot reference marker for tr axis system at th_fc = th(ixtr)

if (use_3d)
   axnams = ['$x_{tr}$'; '$y_{tr}$'; '        ']; col = 5; th_txt = 0;
   rot = [180-phi*180/pi, -th2(ixtr), 0];
   factxt = [1.6 1.4 1.2; 0.5 1.4 1.2; 1 1.4 1.2];
   l = plot_axes3([c2(:,ixtr);0], 10, rot, axnams, col, [],[],[], th_txt, factxt, 0 );
   set(l(1),'color','k');
   set(l([4,7,10]), 'interpreter','latex');
else
   axnams = ['$x_{tr}$'; '$y_{tr}$'; '        ']; col = 5; th_txt = 0;
   l = plot_axes(c2(:,ixtr), 10, -th2(ixtr), -1, axnams, col, ...
                   [],[],[],[], th_txt, [1.2,0.0; 1.2,1.2; 3,2.2], 1 );
   set(l([6,9,10]), 'interpreter','latex');
end

if (print_fig)
   axis off
   if (use_3d)
      print -djpeg95 fc_coord_sys3.jpg
   else
      print -djpeg95 fc_coord_sys.jpg
   end
end
