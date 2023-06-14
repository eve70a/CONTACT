
print_fig = 1;
show_A = 1;
show_B = 0;
show_u_el = 0;

% define products to be plotted

n_curv = 4;
utot_c0 = zeros(n_curv,1); tau_c0 = zeros(n_curv,1);
utot_c1 = zeros(n_curv,1); tau_c1 = zeros(n_curv,1);
utot_c0(1) =  70; tau_c0(1) = 560; utot_c1(1) = 600; tau_c1(1) = 560; % magnetite
utot_c0(2) =  40; tau_c0(2) = 200; utot_c1(2) = 540; tau_c1(2) = 400; % clay
utot_c0(3) =  50; tau_c0(3) = 400; utot_c1(3) = 550; tau_c1(3) = 300; % sand
utot_c0(4) =  20; tau_c0(4) =  20; utot_c1(4) = 520; tau_c1(4) =  30; % MoS2
nams = strvcat('magnetite','clay','sand','MoS2');
% nams = strvcat('G=8, \tau_{c0}=560', 'G=5, \tau_{c0}=200', 'G=8, \tau_{c0}=400', 'G=1, \tau_{c0}=20');

G_el   = tau_c0 ./ utot_c0;
k_u    = (tau_c1 - tau_c0) ./ (utot_c1 - utot_c0);
k_tau  = G_el .* k_u ./ (G_el - k_u);

for ic = 1:n_curv
   disp(sprintf('%9s: utot_c0=%3d, tau_c0=%3d, G_el= %4.2f, k_u=%6.3f, k_tau=%8.4f', nams(ic,:), utot_c0(ic), tau_c0(ic), G_el(ic), k_u(ic), k_tau(ic)));
end

% define range of displacements for curves

u_tot  = [0:600]; 

% define unit-vectors for 1d to 2d conversion

n_utot  = length(u_tot);
e1 = ones(n_curv,1);
e2 = ones(1,n_utot);

% compute displacement - stress curves, tau = tau(u_tot)

tau    = (G_el*u_tot) .* (e1*u_tot <= utot_c0*e2) + ...
         (tau_c0*e2 + (k_u*e2) .* (e1*u_tot-utot_c0*e2)) .* (e1*u_tot > utot_c0*e2);
% tau    = [0:tau_c1];
% u_tot  = (tau/G).*(tau<=tau_c0) + (utot_c0+(tau-tau_c0)/k).*(tau>tau_c0);

% plot curves

figure(1); clf; hold on;
plot(u_tot, tau);
grid on;
axis([0 600 0 620]);
xlabel('total displacement u_{tot} [\mu m]');
ylabel('shear stress \tau [MPa]');
% legend(nams, 'location','south');

text(300, 515, 'k_u= 0', 'horizontalalignment','center');
text(580, 430, sprintf('k_u=%4.1f', k_u(2)), 'horizontalalignment','right');
text(580, 250, sprintf('k_u=%5.1f', k_u(3)), 'horizontalalignment','right');
text(580,  60, sprintf('k_u=%5.2f', k_u(4)), 'horizontalalignment','right');

% show decomposition of u_tot into upl + uel

% point A: clay, utot = 200

if (show_A)
   ic    = 2;    % clay
   utot_A = 200;
   tau_A = tau_c0(ic) + k_u(ic) * (utot_A - utot_c0(ic));
   uel_A = tau_A / G_el(ic);
   upl_A = utot_A - uel_A;

   tmp_y = [ tau_c0(ic) 330 ];
   tmp_x = tmp_y / G_el(ic);
   set(gca,'colororderindex',2);
   plot(tmp_x, tmp_y, '--');

   if (show_u_el)
      x = [0 0 0 uel_A*[1 1 1 1] utot_A*[1 1 1]];
      y = tau_A + 15*[-1 1 0 [0 -1 1 0] [0 -1 1]];
      text(50, 310, 'u_{el}', 'horizontalalignment','right');
   else
      x = [uel_A*[1 1 1 1] utot_A*[1 1 1]];
      y = tau_A + 15*[[0 -1 1 0] [0 -1 1]];
   end
   plot(x, y, 'k', 'linewidth',1);
   plot([ 130 120], [295 270], 'k', 'linewidth',1);
   text(150, 310, 'u_{pl}', 'horizontalalignment','center');
end

% point B: sand, utot = 300

if (show_B)
   ic    = 3;    % sand
   utot_B = 300;
   tau_B = tau_c0(ic) + k_u(ic) * (utot_B - utot_c0(ic));
   uel_B = tau_B / G_el(ic);
   upl_B = utot_B - uel_B;

   if (show_u_el)
      x = [0 0 0 uel_B*[1 1 1 1] utot_B*[1 1 1]];
      y = tau_B + 15*[[-1 1 0] [0 -1 1 0] [0 -1 1]];
   else
      x = [uel_B*[1 1 1 1] utot_B*[1 1 1]];
      y = tau_B + 15*[[0 -1 1 0] [0 -1 1]];
   end
   plot(x, y, 'k', 'linewidth',1);

   % text(50, 310, 'u_{el}', 'horizontalalignment','right');
   plot([ 120 130 ], [345 325 ], 'k', 'linewidth',1);
end

% add line for unloading at point A

tmp_y = tau_A + [ 0 -160 ];
tmp_x = utot_A + (tmp_y - tau_A) / G_el(2);
set(gca,'colororderindex',2);
plot(tmp_x, tmp_y, '--');

text(utot_A+15, tau_A, '\leftarrow unloading', 'rotation', 75, 'horizontalalignment','right');
text(240, 255, 'further loading \rightarrow', 'rotation', 17);

if (print_fig)
   print -djpeg95 hou1997_schm_utot.jpg
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define range for plot of tau_c

u_pl = [ 0 : 600 ];
tau_c = tau_c0*e2 + k_tau * u_pl;

figure(2); clf; hold on;
plot(u_pl, tau_c);
axis([0 600 0 620]);
grid on
xlabel('plastic displacement u^*_{pl} [\mu m]');
ylabel('yield limit \tau_c [MPa]');
% legend(nams, 'location','south');

text(300, 515, 'k_\tau= 0', 'horizontalalignment','center');
text(580, 470, sprintf('k_\\tau=%6.3f', k_tau(2)), 'horizontalalignment','right');
text(580, 250, sprintf('k_\\tau=%7.3f', k_tau(3)), 'horizontalalignment','right');
text(580,  60, sprintf('k_\\tau=%7.4f', k_tau(4)), 'horizontalalignment','right');

% plot position A

set(gca,'colororderindex',2);
plot([3 3], [tau_c0(2) 330], '--');

if (show_A)
   x = [0 0 0 upl_A*[1 1 1]];
   y = tau_A + 15*[-1 1 0 0 -1 1];
   plot(x, y, 'k', 'linewidth',1);
end

plot([0 60], tau_c0(2)*[1 1], 'k--', 'linewidth', 1);
text(65, tau_c0(2), '\tau_{c0}', 'verticalalignment','top');

% plot u_pl for position B

if (show_B)
   x = [0 0 0 upl_B*[1 1 1]];
   y = tau_B + 15*[-1 1 0 0 -1 1];
   plot(x, y, 'k', 'linewidth',1);
end

% add line for unloading at point A

tmp_y = tau_A + [ 0 -160 ];
tmp_x = upl_A + [ 0  0 ];
set(gca,'colororderindex',2);
plot(tmp_x, tmp_y, '--');

text(upl_A+15, tau_A, '\leftarrow unloading', 'rotation', 90, 'horizontalalignment','right');
text(188, 255, 'further loading \rightarrow', 'rotation', 18);

if (print_fig)
   print -djpeg95 hou1997_schm_tauc.jpg
end

