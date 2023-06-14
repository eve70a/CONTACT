
print_fig=0;

set(0,'defaultlinelinewidth',2);
set(0,'defaultaxeslinewidth',2);
set(0,'defaulttextfontsize',15);
set(0,'defaultaxesfontsize',15);
set(0,'defaultfigurepaperpositionmode','auto');

if (~exist('plot_arrow'))
   addpath('../../../../../contact/matlab');
end

figure(1); clf; hold on
axis equal

R=1;
th=0:360; sn=sin(th*pi/180); cs=cos(th*pi/180);

% in vierkant [-0.5,0.5] x [0,1]:
x_spring=[0  0   -1  1   -1  1   -1  1   0   0 ]/2;
y_spring=[0 0.4 0.5 0.7 0.9 1.1 1.3 1.5 1.6 2.0]/2;

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

if (1==1)
   % vertical: spring and damper

   % arrow for variable spring - behind of spring:
   l1=plot_arrow([x1-0.4,y1-R-0.35], [0.8, 0.3], [0 .75 .75], 0.5);

   % spring:
   ofsx = -0.1;
   sclx = 0.4;
   scly = 0.39;
   plot(x1+ofsx+sclx*x_spring, y1-R+scly*(y_spring-1), 'k');

   % unilateral contact allowing separation
   plot(x1+ofsx+sclx*[-0.4 0.4], y1-R-scly+[0.01 0.01], 'k');

   % damper
   ofsx = 0.2;
   sclx = 0.25;
   scly = 0.40;
   l4 = plot(x1+ofsx+sclx*x_damp, y1-R+0.01+scly*(y_damp-1), 'k');

   % horizontal: spring + friction element

   % arrow for variable spring - behind of spring:
   l2=plot_arrow([x2-0.3,y2-R-0.35], [0.65, 0.28], [0 .75 .75], 0.6);

   % in rechthoek [0,0.7]+[0.7,1.7] x [-0.25,0.25]
   sclx = 0.8;
   scly = 0.2;
   yc = (y2 - R)/2;
   plot(x2+sclx*(y_spring-0.5), yc+scly*x_spring,'k');
   % friction: down-arrow:
   plot(x2+sclx*[-0.5 -0.5], yc+[0, -0.1],'k');
   fill(x2+sclx*-0.5+[-0.07 0.07 0], yc+[-0.1, -0.1, -0.2],'k');
   plot((x2+sclx*0.5)*[1 1], [ yc 2.5*yc], 'k');

   if (print_fig)
      axis off
      axis tight
      print -depsc cntc_fric_spring_var.eps
      print -djpeg95 cntc_fric_spring_var.jpg

      delete([l1;l2]);
      print -depsc cntc_fric_spring.eps
      print -djpeg95 cntc_fric_spring.jpg
   end

end

