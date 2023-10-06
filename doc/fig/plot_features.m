
print_fig = 1;

slcs = read_profile('../../testbank/profiles/uk_interrupt_v1.slcs', [], -1, -1, 1000);

figure(2); clf; hold on

for is = 6 : 12
   y = slcs.ysurf(is,:);
   z = slcs.zsurf(is,:);
   x = slcs.s(is) * ones(size(y));
   if (is<=7)
      plot3(x, y, z, 'color',matlab_color(2), 'linewidth',1);
   else
      plot3(x, y, z, 'color',matlab_color(6), 'linewidth',1);
   end
end

islc = [6:12];
for ift = 1 : slcs.nfeat
   ip = slcs.iseg_p(ift);
   plot3( slcs.s(islc), slcs.ysurf(islc,ip), slcs.zsurf(islc,ip), '.:', 'color',matlab_color(4), ...
                                                               'linewidth',1, 'markersize',9 );
end
clear x y z is islc ift ip;

plot_arrow([92,30,0], [15,0,0], 4);
text(103, 34, 0, '$s_{fc}$ (track curve)', 'interpreter','latex', 'fontsize',14);

rot = [0,0,1; 1,0,0; 0,1,0];
circ_arrow([80,41,24], 15, -180, -100, 2, [], [], rot);
text(80, 31, 18, '$s$ (chord length)', 'interpreter','latex', 'fontsize',14);

text(80, 16, 30, '$0$', 'interpreter','latex', 'fontsize',12);
text(80, 22-2, -1, sprintf('$%5.2f$',slcs.s_feat(6,end-1)), 'interpreter','latex', 'fontsize',12, ...
                                        'horizontalalignment','left', 'verticalalignment','top' );
text(80, 100, 6, sprintf('`$%4.0f$''',slcs.s_feat(6,end)), 'interpreter','latex', 'fontsize',12, ...
                                        'horizontalalignment','right' );

grid on;
xlabel('s_{fc} [mm]');
ylabel('y_r [mm]');
zlabel('z_r [mm]');
set(gca,'ydir','reverse', 'zdir','reverse');
set(gca,'dataaspectratio',[0.5 1 1]);
axis([80 140 -45 100 -10 55]);
view([-110 15]);

if (print_fig)
   print -djpeg95 uk_interrupt_v1.jpg
end
