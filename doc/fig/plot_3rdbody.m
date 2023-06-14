
print_fig = 0;

set(0, 'defaultlinelinewidth', 1);
set(0, 'defaulttextfontsize', 12);

col1 = matlab_color(5); % rigid slip
col2 = matlab_color(2); % stress
col3 = matlab_color(1); % elastic displacement
col4 = matlab_color(3); % plastic displacement

% configure deformations:

uel_1 = 1.6;
upl_1 = 0.8; zpl_1 = 0.5;
utot_1 = upl_1 + uel_1;
uel_2 = 1.2;
uel_3 = 1.6;
upl_3 = 2.0;
utot_3 = upl_3 + uel_3;
utot = utot_1 + uel_2 + utot_3;
px   = 3.5;

% compute deformation profile

zdef = [0:0.1:3]; z0 = 9;
udef = (zdef-z0).^2 / z0^2;

upl_z = max(0, zpl_1-zdef) / zpl_1;

% define shape of upper/lower body

mid  = 10;
xsrf = 0:0.1:3; ysrf = (xsrf/3).^2; 
xsrf = [-fliplr(xsrf)-mid/2, xsrf+mid/2]; ysrf = [ fliplr(ysrf), ysrf];


figure(1); clf; hold on;
axis equal
axis([-13 12.5 -6 6]);
grid on;

% third body layer is horizontally centered
% plot contours of the layer

plot([-6 -5], [ 1  1], 'k--');
plot([-5  5], [ 1  1], 'k-');
plot([ 5  6], [ 1  1], 'k--');
plot([-6 -5], [-1 -1], 'k--');
plot([-5  5], [-1 -1], 'k-');
plot([ 5  6], [-1 -1], 'k--');
text(6, 0, 'interface layer (3), $G^{(3)}$', 'interpreter','latex');
% text(5, 0,'interface layer (3), $G^{(3)}\ll G^{(1)}$', 'interpreter','latex');
h= fill([-6 6 6 -6 -6], [1 1 -1 -1 1], matlab_color(1));
set(h,'facealpha',0.1, 'edgecolor','none');

% plot surface of upper & lower bodies

x1 = -utot_3/2 - utot_1; y1 = 3;
plot(x1+xsrf, y1+ysrf, 'k-');
text(6, 4.5, 'body (1), $G^{(1)}$', 'interpreter','latex');
h= fill(x1+xsrf, y1+ysrf, matlab_color(1));
set(h,'facealpha',0.1, 'edgecolor','none');

x2 =  utot_3/2 + uel_2; y2 = -3;
plot(x2+xsrf, y2-ysrf, 'k-');
text(6, -5.0, 'body (2), $G^{(2)}$', 'interpreter','latex');
% text(5, -5.0, 'body (2), $G^{(2)}>G^{(1)}$', 'interpreter','latex');
h= fill(x2+xsrf, y2-ysrf, matlab_color(1));
set(h,'facealpha',0.1, 'edgecolor','none');

% plot reference for elastic deformation in upper & lower bodies

plot( x1*[1 1], [ 5.5,  2.3], ':', 'color',col3);
plot( x2*[1 1], [-5.5, -2.3], ':', 'color',col3);

% plot deformation on upper & lower body

plot(x1+uel_1*udef+upl_1*upl_z, y1+zdef, '-', 'color',col3);
plot(x1+upl_1*upl_z, y1+zdef, '-', 'color',col4);
plot(x2-uel_2*udef, y2-zdef, '-', 'color',col3);

fac = (zpl_1-0.2)/zpl_1;
plot_arrow([x1, y1+0.2], [upl_1*fac, 0], col4, 2);
plot_arrow([x1+upl_1*fac, y1+0.2], [uel_1, 0], col3);
plot_arrow([x2, y2-0.2], [ -uel_2, 0], col3);
text(x1-0.2, y1+0.4, '$u_{pl}^{(1)}$', 'interpreter','latex', ...
        'horizontalalignment','right', 'verticalalignment','bottom');
text(x1+uel_1+0.5, y1+0.4, '$u_{el}^{(1)}$', 'interpreter','latex', ...
        'horizontalalignment','left', 'verticalalignment','bottom');
text(x2-uel_2-0.2, y2-0.4, '$u_{el}^{(2)}$', 'interpreter','latex', ...
        'horizontalalignment','right', 'verticalalignment','top');

% plot overall rigid shift w

plot( x2*[1 1 NaN 1 1], [ 0.5 1.5 NaN 4.5, 5.5], ':', 'color',col3);
plot_arrow([x2, 5], [ -utot, 0], col1, 0.2);
text(0, 5, '$w$', 'interpreter','latex', 'verticalalignment','bottom');

% plot reference for elastic deformation in 3rd body

x = upl_3/2 * [  1    1 NaN -1 -1   ];
y =           [ -1.7  0 NaN  0  1.7 ];
plot(x, y, ':', 'color',col3);

% plot elastic deformation in 3rd body

x = [ utot_3, utot_3, NaN, -utot_3, -utot_3 ] / 2;
y = [   -3,     -1,   NaN,    1,       3   ];
plot(x, y, '--', 'color',col3);
x = [ utot_3, upl_3, NaN, -upl_3, -utot_3 ] / 2;
y = [   -1,      0,   NaN,    0,       1  ];
plot(x, y, '-', 'color',col3);

fac = 0.8 / 1;
plot_arrow([-upl_3/2,  1-0.2], [-fac*uel_3/2, 0], col3, 1.5);
plot_arrow([ upl_3/2, -1+0.2], [ fac*uel_3/2, 0], col3, 1.5);

text(-utot_3/2,  1-0.2, '$u_{el}^{(3,top)}$', 'interpreter','latex', ...
        'horizontalalignment','right', 'verticalalignment','top');
text( utot_3/2, -1+0.2, '$u_{el}^{(3,bot)}$', 'interpreter','latex', ...
        'horizontalalignment','left', 'verticalalignment','bottom');

% plot plastic deformation in 3rd body

x = upl_3/2 * [ -1 -1 ];
y = 0.25 * [1 -1 ];
plot(x, y, 'color',col4);
plot_arrow([-upl_3/2,0], [upl_3,0], col4, 0.5);

text(-0.4, -0.1, '$u_{pl}^{(3)}$', 'interpreter','latex', ...
        'horizontalalignment','left', 'verticalalignment','middle');

% plot stresses on 1st, 3rd-top, 3rd-bot, 2nd bodies

plot_arrow([-utot_3/2,  3-0.2], [ px, 0], col2, 0.4);
plot_arrow([-utot_3/2,  1+0.2], [-px, 0], col2, 0.4);
plot_arrow([ utot_3/2, -1-0.2], [ px, 0], col2, 0.4);
plot_arrow([ utot_3/2, -3+0.2], [-px, 0], col2, 0.4);

text(x1+utot_1+px/2, y1-0.4, '$p^{(1)}$', 'interpreter','latex', ...
        'horizontalalignment','left', 'verticalalignment','top');
text(x1+utot_1-px/3, 1+0.3, '$p^{(3,top)}=-p^{(1)}$', 'interpreter','latex', ...
        'horizontalalignment','right', 'verticalalignment','bottom');
text(x2-uel_2+px/2, -1-0.4, '$p^{(3,bot)}=-p^{(3,top)}$', 'interpreter','latex', ...
        'horizontalalignment','left', 'verticalalignment','top');
text(x2-uel_2-px/3, y2+0.3, '$p^{(2)}=-p^{(3,bot)}$', 'interpreter','latex', ...
        'horizontalalignment','right', 'verticalalignment','bottom');

% plot legend

plot_arrow([-13, -2], [2, 0], col1, 0.7);
text(-10.3, -2, 'rigid displ.', 'interpreter','latex');
plot_arrow([-13, -3], [2, 0], col2, 0.7);
text(-10.3, -3, 'stress', 'interpreter','latex');
plot_arrow([-13, -4], [2, 0], col3, 0.7);
text(-10.3, -4, 'elastic displ.', 'interpreter','latex');
plot_arrow([-13, -5], [2, 0], col4, 0.7);
text(-10.3, -5, 'plastic displ.', 'interpreter','latex');

if (print_fig)
   axis off
   print -djpeg95 tractions_def2.jpg
end
