% produces figures:
% iversion = 1: just a friction element
% iversion = 2: friction element + damper
% iversion = 3: friction element + damper + spring
%
%     print -djpeg95 cntc_fric_spring_dashpot_var.jpg
%     print -djpeg95 cntc_fric_spring_dashpot.jpg
%     print -djpeg95 cntc_spring_dashpot.jpg
%     print -djpeg95 cntc_spring_dashpot_var.jpg
%

iversion = 2;
print_fig = 0;

figure(1); clf; hold on
axis equal

R=1;
th=0:360; sn=sin(th*pi/180); cs=cos(th*pi/180);

% hor.spring in vierkant [0,1] x [0,1]:
x_spring=[0  0   -1  1   -1  1   -1  1   0   0 ]/2;
y_spring=[0 0.4 0.5 0.7 0.9 1.1 1.3 1.5 1.6 2.0]/2;

% vert.damper in rechthoek [-0.25,0.25] x [0,1]:
x_damp_lrg  =[-1 -1 1 1 NaN -0.6 0.6 NaN 0 0 NaN  0  0]/4;
y_damp_lrg  =[ 2  1 1 2 NaN  1.5 1.5 NaN 1 0 NaN 1.5 3]/3;

% vert.damper in rechthoek [-0.25,0.25] x [0,1]:
x_damp_sml  =[-1 -1 1 1 NaN 0 0 NaN  0  0]/4;
y_damp_sml  =[ 2  1 1 2 NaN 1 0 NaN 1.3 3]/3;

% modified vert.damper in rechthoek [-0.25,0.25] x [0,1]:
x_damp_star  =[-1 -1 1 1 NaN 0 0 NaN  0  0 NaN -0.6 0.6 NaN -0.6 0.6 ]/4;
y_damp_star  =[ 2  1 1 2 NaN 1 0 NaN 1.2 3 NaN  1.35 1.8 NaN  1.8 1.35 ]/3;

% vert.damper + spring in vierkant [0,1] x [0,1]:
x_par_spring_damp = [0.35+x_spring*0.4, NaN, -0.35+x_damp_lrg*0.8 NaN, ...
                    -0.35 0.35 NaN, 0 0    NaN, -0.35 0.35 NaN, 0    0];
x_par_spring_damp = (x_par_spring_damp+0.55) / 1.1;
y_par_spring_damp = [0.15+0.7*y_spring, NaN,  0.15+0.7*y_damp_lrg NaN,  ...
                     0.15 0.15 NaN, 0 0.15 NaN,  0.85 0.85 NaN, 0.85 1];

% horizontaal, in rechthoek [0, 1.6] x [-0.25,0.25]
% nieuw: damp_star
y_damp_star(7)=0.1;
y_damp_star(10)=0.8;
% figure(2); clf; plot(x_damp_star, y_damp_star); figure(1);
x_ser_spring_damp = [y_damp_star, NaN, y_spring+0.7] - 0.1;
y_ser_spring_damp = [x_damp_star, NaN, x_spring/2];
y_damp_star(7)=0;
y_damp_star(10)=1;


% center positions of wheels 1, 2:
x1 = 0; y1 = R+0.4;
x2 = 3; y2 = y1;

% first wheel

plot(x1, y1, 'b.', 'markersize', 30)
plot(x1+R*cs, y1+R*sn, 'b');

% second wheel

plot(x2, y2, 'b.', 'markersize', 30);
plot(x2+R*cs, y2+R*sn, 'b');

% rail

plot([x1-1.5*R x2+1.5*R],  0.00*R*[1 1], 'r');
plot([x1-1.5*R x2+1.5*R], -0.30*R*[1 1], 'r');

text(x1, -0.5, 'normal', 'horizontalalignment','center');
text(x2, -0.5, 'tangential', 'horizontalalignment','center');

axis([x1-1.5*R x2+1.5*R -0.7 y1+1.2*R]);

if     (iversion==1)
   % first variant:  friction element

   if (1==1)
      % vertical: spring + separation element
      sclx = 0.4;
      scly = 0.39;
      l1=plot_arrow([x1-0.4,y1-R-0.35], [0.8, 0.3], [0 .75 .75], 0.5);
      plot(x1+sclx*x_spring, y1-R+scly*(y_spring-1), 'k');

      % normal: unilateral contact allowing separation
      ofsx = -0.1;
      plot(x1+ofsx+sclx*[-0.4 1.0], y1-R-scly+[0.01 0.01], 'k');

   else
      % vertical: spring and damper
      ofsx = -0.1;
      sclx = 0.4;
      scly = 0.39;
      l1=plot_arrow([x1-0.4,y1-R-0.35], [0.8, 0.3], [0 .75 .75], 0.5);
      plot(x1+ofsx+sclx*x_spring, y1-R+scly*(y_spring-1), 'k');

      % normal: unilateral contact allowing separation
      plot(x1+ofsx+sclx*[-0.4 1.0], y1-R-scly+[0.01 0.01], 'k');

      % normal: vertical damper
      ofsx = 0.2;
      sclx = 0.25;
      scly = 0.40;
      l4 = plot(x1+ofsx+sclx*x_damp_lrg, y1-R+0.01+scly*(y_damp_lrg-1), 'k');
   end

   % horizontal: damper
   % in rechthoek [0,1.7] x [-0.25,0.25]
   sclx = 0.5;
   scly = 0.45;
   yc = (y2 - R)/2;
   l2=plot_arrow([x2-0.4,y2-R-0.35], [0.8, 0.3], [0 .75 .75], 0.5);
   plot(x2+sclx*[-0.5 -0.5 0.5 0.5], [ 0  yc yc 2.3*yc], 'k');

   % friction: down-arrow:
   l3=fill(x2+sclx*-0.5+[-0.07 0.07 0], yc+[-0.1, -0.1, -0.2],'k');

   if (print_fig)
      set(gcf, 'paperpositionmode','auto');

      axis off
      axis tight
      set([l1;l2],'visible','off');
      print -djpeg95 cntc_fricelem.jpg

      set([l1;l2],'visible','on');
      print -djpeg95 cntc_fricelem_var.jpg
   end

elseif (iversion==2)
   % second variant: friction element + damper

   if (0==1)
      % vertical: just a spring
      sclx = 0.4;
      scly = 0.4;
      l1=plot_arrow([x1-0.4,y1-R-0.35], [0.8, 0.3], [0 .75 .75], 0.5);
      plot(x1+sclx*x_spring, y1-R+scly*(y_spring-1), 'k');
   else
      % vertical: spring and damper
      ofsx = -0.1;
      sclx = 0.4;
      scly = 0.39;
      l1=plot_arrow([x1-0.4,y1-R-0.35], [0.8, 0.3], [0 .75 .75], 0.5);
      plot(x1+ofsx+sclx*x_spring, y1-R+scly*(y_spring-1), 'k');

      % normal: unilateral contact allowing separation
      plot(x1+ofsx+sclx*[-0.4 1.0], y1-R-scly+[0.01 0.01], 'k');

      % normal: vertical damper
      ofsx = 0.2;
      sclx = 0.25;
      scly = 0.40;
      l4 = plot(x1+ofsx+sclx*x_damp_lrg, y1-R+0.01+scly*(y_damp_lrg-1), 'k');
   end

   % horizontal: damper
   % in rechthoek [0,1.7] x [-0.25,0.25]
   sclx = 0.6;
   scly = 0.45;
   yc = (y2 - R)/2;
   l2=plot_arrow([x2-0.4,y2-R-0.35], [0.8, 0.3], [0 .75 .75], 0.5);
   plot(x2+sclx*(y_damp_star-0.5), yc+scly*x_damp_star,'k');
   plot((x2-0.5*sclx)*[1 1], [ 0  yc ], 'k');
   plot((x2+0.5*sclx)*[1 1], [ yc 2.3*yc], 'k');

   % friction: down-arrow:
   l3=fill(x2+sclx*-0.5+[-0.07 0.07 0], yc+[-0.1, -0.1, -0.2],'k');

   if (print_fig)
      set(gcf, 'paperpositionmode','auto');

      axis off
      axis tight
      print -djpeg95 cntc_fric_dashpot_var.jpg

      set([l1;l2],'visible','off');
      print -djpeg95 cntc_fric_dashpot.jpg

      set([l3],'visible','off');
      print -djpeg95 cntc_dashpot.jpg

      set([l1;l2],'visible','on');
      print -djpeg95 cntc_dashpot_var.jpg
   end
elseif (iversion==3)
   % third variant: friction element + damper + spring

   % vertical: spring and damper
   ofsx = -0.1;
   sclx = 0.4;
   scly = 0.39;
   l1=plot_arrow([x1-0.4,y1-R-0.35], [0.8, 0.3], [0 .75 .75], 0.5);
   plot(x1+ofsx+sclx*x_spring, y1-R+scly*(y_spring-1), 'k');

   % normal: unilateral contact allowing separation
   plot(x1+ofsx+sclx*[-0.4 1.0], y1-R-scly+[0.01 0.01], 'k');

   % normal: vertical damper
   ofsx = 0.2;
   sclx = 0.25;
   scly = 0.40;
   l4 = plot(x1+ofsx+sclx*x_damp_lrg, y1-R+0.01+scly*(y_damp_lrg-1), 'k');

   % tangential: spring + damper
   % in rechthoek [0,1.7] x [-0.25,0.25]
   sclx = 0.6;
   scly = 0.45;
   yc = (y2 - R)/2;
   l2=plot_arrow([x2-0.4,y2-R-0.35], [0.8, 0.3], [0 .75 .75], 0.5);
   plot(x2+sclx*(x_ser_spring_damp-0.8), yc+scly*y_ser_spring_damp,'k');
   plot((x2-0.8*sclx)*[1 1], [ 0  yc ], 'k');
   plot((x2+0.8*sclx)*[1 1], [ yc 2.5*yc], 'k');

   % tangential: friction: down-arrow:
   l3=fill(x2+sclx*-0.8+[-0.07 0.07 0], yc+[-0.1, -0.1, -0.2],'k');

   if (print_fig)
      set(gcf, 'paperpositionmode','auto');

      axis off
      axis tight
      print -djpeg95 cntc_fric_spring_dashpot_var.jpg

      set([l1;l2],'visible','off');
      print -djpeg95 cntc_fric_spring_dashpot.jpg

      set([l3],'visible','off');
      print -djpeg95 cntc_spring_dashpot.jpg

      set([l1;l2],'visible','on');
      print -djpeg95 cntc_spring_dashpot_var.jpg
   end

end

