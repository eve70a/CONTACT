
% plastic_one_shift: shift along x-direction, then along path, then y-direction

print_fig = 0;
show_fig = [ 1:3 ];

expnam = 'plastic_one_shift';
ncase = 27;

col = matlab_color(1:5);

% initialize at zero
tauc  = zeros(1,ncase);
pn    = zeros(1,ncase);
px    = zeros(1,ncase);
py    = zeros(1,ncase);
uelx  = zeros(1,ncase);
uely  = zeros(1,ncase);
duplx = zeros(1,ncase);
duply = zeros(1,ncase);

for t = 1:ncase
   s = loadcase(expnam,t);
   time = t+1;
   tauc(time)  = s.taucrt;
   pn(time)    = s.pn;
   px(time)    = s.px;
   py(time)    = s.py;
   uelx(time)  = s.ux;
   uely(time)  = s.uy;
   uplx(time)  = s.uplsx;
   uply(time)  = s.uplsy;
   shft_abs(time) = s.shft;
end

duplx  = [0, diff(uplx)];
duply  = [0, diff(uply)];
utotx = uelx + uplx;
utoty = uely + uply;

pt_abs   = sqrt(px.^2    + py.^2);
uel_abs  = sqrt(uelx.^2  + uely.^2);
upl_abs  = sqrt(uplx.^2  + uply.^2);
utot_abs = sqrt(utotx.^2 + utoty.^2);
ustar    = cumsum( sqrt( duplx.^2 + duply.^2 ) );

pn(1)   = pn(2);
tauc(1) = tauc(2);

T = 0:ncase;

if (any(show_fig==1))
   figure(1); clf; hold on;
   plot(uplx+uelx, uply+uely, '-*', 'color',col(2,:));
   plot(uplx, uply, '-*', 'color',col(3,:));
   axis equal
   grid on;
   xlabel('u_x');
   ylabel('u_y');
   title('displacements path');
   legend('u_{tot}','u_{pl}','location','northwest');

   if (print_fig)
      print -djpeg95 'figs/one_el_path.png'
   end
end

if (any(show_fig==2))
   figure(2); clf; hold on;
   plot(T, px, '.-', 'color',col(2,:) );
   plot(T, py, '.-', 'color',col(4,:) );
   plot(T, pt_abs, '--*', 'color',col(1,:) );
   plot(T, tauc, '--', 'color',col(3,:) );
   axis([0 30 0 5.5]);
   grid on;
   legend('p_x','p_y','|p|','\tau_c', 'location','southeast');
   xlabel('time t');
   ylabel('traction [N/mm^2]');
   title('tractions');

   if (print_fig)
      print -djpeg95 'figs/one_el_tract.png'
   end
end

if (any(show_fig==3))
   figure(3); clf; hold on;
   plot(T, utot_abs+shft_abs, '-o', 'color',col(1,:));
   plot(T, uel_abs,           '-*', 'color',col(2,:));
   plot(T, ustar,             '-*', 'color',col(3,:));
   plot(T, shft_abs   ,       '-*', 'color',col(4,:));
   legend('total |u_{tot}|', ...
          'elastic |u_{el}|', ...
          'plastic |u_{pl}|', ...
          'slip |s|', ...
          'Location','NorthWest');
   axis([0 30 0 9]); grid on;
   xlabel('time t');
   ylabel('displacement');
   title('cumulative displacement');

   if (print_fig)
      print -djpeg95 'figs/one_el_displc.png'
   end
end

