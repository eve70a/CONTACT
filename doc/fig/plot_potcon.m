
print_fig = 1;

set(0, 'defaultlinelinewidth', 1);
set(0, 'defaulttextfontsize', 9);

figure(1); clf; hold on;

% size of potential contact area:

xl=-10  ; dx=2.0; xh= 6; mx=round((xh-xl)/dx);
yl= -4.8; dy=1.6; yh= 4.8; my=round((yh-yl)/dy);
xi=xl-dx/2+[1:mx]*dx;
yj=yl-dy/2+[1:my]*dy;

% plot the elements:

for ix=0:mx
   plot(xl+ix*dx*[1 1], [yl,yh], 'b');
end
for iy=0:my
   plot([xl,xh], yl+iy*dy*[1 1], 'b');
end

% plot the element centers:

for iy=1:my
   for ix=1:mx,
      plot(xi(ix),yj(iy),'b.','markersize',12);
   end
end

% numbering of the elements

l=text(xi(1),yj(1)-dy,'ix=1','horizontalalignment','center','verticalalignment','top');
plot(xi(1)*[1 1], yj(1)+dy*[-0.5,-0.9],'k');
l=text(xi(mx),yj(1)-dy,'ix=mx','horizontalalignment','center','verticalalignment','top');
plot(xi(mx)*[1 1], yj(1)+dy*[-0.5,-0.9],'k');
l=text(xi(floor(mx/2))-0.5*dx, yj(1)-dy, 'x \rightarrow forward rolling', 'verticalalignment','top');

l=text(xi(1)-dx,yj(1),'iy=1','horizontalalignment','right','verticalalignment','middle');
plot(xi(1)+dx*[-0.5,-0.9], yj(1)*[1,1],'k');
l=text(xi(1)-dx,yj(my),'iy=my','horizontalalignment','right','verticalalignment','middle');
plot(xi(1)+dx*[-0.5,-0.9], yj(my)*[1,1],'k');
l=text(xi(1)-dx,yj(round(my/2)),'y \rightarrow','rotation',90,'horizontalalignment','right','verticalalignment','middle');

% element sizes dx, dy

plot(xi(7)+dx*[-.5 .5],yj(my)+dy*[1 1],'k');
plot(xi(7)-dx/2+[0.1 0 0.1], yj(my)+dy+[-0.1 0 0.1],'k');
plot(xi(7)+dx/2-[0.1 0 0.1], yj(my)+dy+[-0.1 0 0.1],'k');
text(xi(7), yj(my)+1.2*dy, '\delta{}x', 'horizontalalignment','center', 'verticalalignment','bottom')

plot(xi(mx)+dx*[1 1],yj(5)+dy*[-.5 .5],'k');
plot(xi(mx)+dx+[-0.1 0 0.1],yj(5)-dy/2+[0.1 0 0.1],'k');
plot(xi(mx)+dx+[-0.1 0 0.1],yj(5)+dy/2-[0.1 0 0.1],'k');
text(xi(mx)+1.2*dx, yj(5), '\delta{}y', 'horizontalalignment','left', 'verticalalignment','middle')

% element numbering

l=text(xi(1)+0.1*dx, yj(1)+0.05*dy, '1','verticalalignment','bottom');
l=text(xi(2)+0.1*dx, yj(1)+0.05*dy, '2','verticalalignment','bottom');
l=text(xi(1),        yj(2)+0.05*dy, 'mx+1','verticalalignment','bottom', 'horizontalalignment','center');
l=text(xi(mx),       yj(1)+0.05*dy, 'I=mx','verticalalignment','bottom', 'horizontalalignment','center');
l=text(xi(mx),       yj(my)+0.05*dy, 'I=npot','verticalalignment','bottom', 'horizontalalignment','center');
l=text(xi(mx),       yj(my-2)+0.05*dy, 'iy\cdot{}mx','verticalalignment','bottom', 'horizontalalignment','center');

% notations (xl,yl), (xc1,yc1)

th = 0:5:360;
r=0.12;
plot(xl+r*cos(th), yl+r*sin(th), 'r-');
plot(xl-[0 0.4*dx], yl-[0 0.4*dy], 'k-');
text(xl-dx/2, yl-dy/2, '(xl,yl)','horizontalalignment','right','verticalalignment','top')

plot(xl+dx/2+r*cos(th), yl+dy/2+r*sin(th), 'r-');
plot(xl+dx*[0.5 1.2], yl+dy*[0.5 -0.4], 'k-');
text(xl+dx, yl-dy/2, '(xc1,yc1)','horizontalalignment','left','verticalalignment','top')

plot(xh+r*cos(th), yh+r*sin(th), 'r-');
plot(xh+[0 dx/2], yh+[0 0.3*dy], 'k-');
text(xh+.6*dx, yh+dy/2, '(xh,yh)')

plot(xh-dx/2+r*cos(th), yh-dy/2+r*sin(th), 'r-');
plot(xh+dx*[-.5 .5], yh+dy*[-.5 -.5], 'k-');
text(xh+0.6*dx, yh-.5*dy, '(xcm,ycm)')

% plot a hypothetical actual contact area:

aa = 3.3*dx; bb = 2.1*dy;
xm = (xl + xh) / 2;
ym = (yl + yh) / 2;
th = [0:360] * pi/180;
plot(xm+aa*cos(th), ym+bb*sin(th), 'k-');

plot([xm+aa*cos(th(355)), 7], [ym+bb*sin(th(355)), -1], 'k-');
l=text(7.2, -1.1, 'actual contact area','horizontalalignment','left','verticalalignment','middle');
plot([6, 7], [-2.1 -2.4], 'b-');
l=text(7.2, -2.6, 'potential contact area','horizontalalignment','left','verticalalignment','middle');

% set axis, print figure

axis equal
axis([-11 11 -6 6]);
axis off

if (print_fig)
   set(gcf,'paperpositionmode','auto');
   print -djpeg95 potcon_numbering.jpg
end
