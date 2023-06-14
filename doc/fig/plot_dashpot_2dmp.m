
print_fig=0;

set(0,'defaultlinelinewidth',2);
set(0,'defaultaxeslinewidth',2);
set(0,'defaulttextfontsize',15);
set(0,'defaultaxesfontsize',15);
set(0,'defaultfigurepaperpositionmode','auto');

figure(1); clf; hold on
axis equal

R=1;
th=0:360; sn=sin(th*pi/180); cs=cos(th*pi/180);

% in vierkant [0,1] x [0,1]:
x_spring=[0  0   -1  1   -1  1   -1  1   0   0 ]/2;
y_spring=[0 0.4 0.5 0.7 0.9 1.1 1.3 1.5 1.6 2.0]/2;

% in rechthoek [-0.25,0.25] x [0,1]:
x_damp_lrg  =[-1 -1 1 1 NaN -0.6 0.6 NaN 0 0 NaN  0  0]/4;
y_damp_lrg  =[ 2  1 1 2 NaN  1.5 1.5 NaN 1 0 NaN 1.5 3]/3;

% in rechthoek [-0.25,0.25] x [0,1]:
x_damp_sml  =[-1 -1 1 1 NaN 0 0 NaN  0  0]/4;
y_damp_sml  =[ 2  1 1 2 NaN 1 0 NaN 1.3 3]/3;

% in vierkant [0,1] x [0,1]:
x_par_spring_damp = [0.35+x_spring*0.4, NaN, -0.35+x_damp_lrg*0.8 NaN, ...
                    -0.35 0.35 NaN, 0 0    NaN, -0.35 0.35 NaN, 0    0];
x_par_spring_damp = (x_par_spring_damp+0.55) / 1.1;
y_par_spring_damp = [0.15+0.7*y_spring, NaN,  0.15+0.7*y_damp_lrg NaN,  ...
                     0.15 0.15 NaN, 0 0.15 NaN,  0.85 0.85 NaN, 0.85 1];

% horizontaal, in rechthoek [0, 1.6] x [-0.25,0.25]
y_damp_lrg(end-3)=0.1;
y_damp_lrg(end)=0.8;
x_ser_spring_damp = [y_damp_lrg, NaN, y_spring+0.7] - 0.1;
y_ser_spring_damp = [x_damp_lrg, NaN, x_spring/2];
y_damp_lrg(end-3)=0;
y_damp_lrg(end)=1;

% center positions of wheels 1, 2:
x1 = 0; y1 = R+0.4;
x2 = 3; y2 = y1;

% first wheel

plot(x1, y1, '.', 'markersize', 30)
plot(x1+R*cs, y1+R*sn, 'b');

% second wheel

plot(x2, y2, '.', 'markersize', 30);
plot(x2+R*cs, y2+R*sn, 'b');

% rail

plot([x1-1.5*R x2+1.5*R],  0.00*R*[1 1], 'r');
plot([x1-1.5*R x2+1.5*R], -0.30*R*[1 1], 'r');

text(x1, -0.5, 'normal', 'horizontalalignment','center');
text(x2, -0.5, 'tangential', 'horizontalalignment','center');

axis([x1-1.5*R x2+1.5*R -0.7 y1+1.2*R]);

if (0==1)
   % first version: parallel spring and damper

   % vertical: spring
   sclx = 0.4;
   scly = 0.6;
   plot(x1+sclx*x_spring, y1-R+scly*(y_spring-1), 'k');

   % horizontal: spring + damper
   sclx = 0.45;
   scly = 0.6;
   yc = y2 - R + sclx*(-0.65);
   plot(x2+scly*(y_par_spring_damp-0.5), y2-R+sclx*(x_par_spring_damp-1.15),'k');
   plot((x2-0.5*scly)*[1 1], [ 0  yc ], 'k');
   plot((x2+0.5*scly)*[1 1], [ yc 2.1*yc], 'k');

elseif (1==1)
   % second version: serial spring and damper

   % vertical: spring and two dampers
   ofsx = 0;
   sclx = 0.4;
   scly = 0.39;
   l1=plot_arrow([x1-0.4,y1-R-0.35], [0.8, 0.3], [0 .75 .75], 0.5);
   plot(x1+ofsx+sclx*x_spring, y1-R+scly*(y_spring-1), 'k');

   % unilateral contact allowing separation
   plot(x1+ofsx+sclx*[-0.8 0.8], y1-R-scly+[0.01 0.01], 'k');

   % dampers
   ofsx = 0.27;
   sclx = 0.15;
   scly = 0.41;
   l4 = plot(x1+ofsx+sclx*x_damp_sml, y1-R+0.03+scly*(y_damp_sml-1), 'k');

   ofsx = -0.27;
   l5 = plot(x1+ofsx+sclx*x_damp_sml, y1-R+0.03+scly*(y_damp_sml-1), 'k');

   % horizontal: spring + damper
   % in rechthoek [0,1.7] x [-0.25,0.25]
   sclx = 0.6;
   scly = 0.45;
   yc = (y2 - R)/2;
   l2=plot_arrow([x2-0.4,y2-R-0.35], [0.8, 0.3], [0 .75 .75], 0.5);
   plot(x2+sclx*(x_ser_spring_damp-0.8), yc+scly*y_ser_spring_damp,'k');
   plot((x2-0.8*sclx)*[1 1], [ 0  yc ], 'k');
   plot((x2+0.8*sclx)*[1 1], [ yc 2.5*yc], 'k');

   % friction: down-arrow:
   l3=fill(x2+sclx*-0.8+[-0.07 0.07 0], yc+[-0.1, -0.1, -0.2],'k');

   if (print_fig)
      axis off
      axis tight
      print -depsc cntc_fric_spring_dashpot2_var.eps
      print -djpeg95 cntc_fric_spring_dashpot2_var.jpg

      % set([l1;l2],'visible','off');
      % print -depsc cntc_fric_spring_dashpot2.eps
      % print -djpeg95 cntc_fric_spring_dashpot2.jpg

      % set([l3],'visible','off');
      % print -depsc cntc_spring_dashpot2.eps
      % print -djpeg95 cntc_spring_dashpot2.jpg

      % set([l1;l2],'visible','on');
      % print -depsc cntc_spring_dashpot2_var.eps
      % print -djpeg95 cntc_spring_dashpot2_var.jpg
   end

else
   % third version: just a damper

   % vertical: spring
   sclx = 0.4;
   scly = 0.4;
   l1=plot_arrow([x1-0.4,y1-R-0.35], [0.8, 0.3], [0 .75 .75], 0.5);
   plot(x1+sclx*x_spring, y1-R+scly*(y_spring-1), 'k');

   % horizontal: damper
   % in rechthoek [0,1.7] x [-0.25,0.25]
   sclx = 0.6;
   scly = 0.45;
   yc = (y2 - R)/2;
   l2=plot_arrow([x2-0.4,y2-R-0.35], [0.8, 0.3], [0 .75 .75], 0.5);
   plot(x2+sclx*(y_damp_lrg-0.5), yc+scly*x_damp_lrg,'k');
   plot((x2-0.5*sclx)*[1 1], [ 0  yc ], 'k');
   plot((x2+0.5*sclx)*[1 1], [ yc 2.3*yc], 'k');

   % friction: down-arrow:
   l3=fill(x2+sclx*-0.5+[-0.07 0.07 0], yc+[-0.1, -0.1, -0.2],'k');

   if (print_fig)
      axis off
      axis tight
      print -depsc cntc_fric_dashpot_var.eps
      print -djpeg95 cntc_fric_dashpot_var.jpg

      set([l1;l2],'visible','off');
      print -depsc cntc_fric_dashpot.eps
      print -djpeg95 cntc_fric_dashpot.jpg

      set([l3],'visible','off');
      print -depsc cntc_dashpot.eps
      print -djpeg95 cntc_dashpot.jpg

      set([l1;l2],'visible','on');
      print -depsc cntc_dashpot_var.eps
      print -djpeg95 cntc_dashpot_var.jpg
   end
end

