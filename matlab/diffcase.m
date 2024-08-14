
function [ dif ] = diffcase(sol1, sol2)

%
% function [ dif ] = diffcase(sol1, sol2)
%
% Compute the difference of the results for two calculations sol1 and sol2.
%
% dif         - output struct "sol1 - sol2"
% sol1, sol2  - structs with surface tractions as defined by loadcase.
%

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

dif = [];

% compare arrays of sol structures using recursion (esp. w/r contacts with multiple patches)

if (length(sol1)>1 | length(sol2)>1)
   if (length(sol1) ~= length(sol2))
      disp(sprintf('ERROR: sol1 has %d contact patches whereas sol2 has %d patches', ...
                        length(sol1), length(sol2)));
      return;
   end

   dif = diffcase(sol1(1), sol2(1));
   for isol = 2 : length(sol1)
      dif(isol) = diffcase(sol1(isol), sol2(isol));
   end
   return
end


% Check equality of the structures w.r.t. grids used

if (sol1.mx ~= sol2.mx | sol1.my ~= sol2.my)
   disp(sprintf('The two structs concern different grids: %dx%d, %dx%d',...
        sol1.mx, sol1.my, sol2.mx, sol2.my));
   return
end

grd1 = [sol1.xl, sol1.yl, sol1.dx, sol1.dy];
grd2 = [sol2.xl, sol2.yl, sol2.dx, sol2.dy];
tol  = 0.01 * max(sol1.dx,sol1.dy);
if (any( abs(grd1-grd2) > tol ))
   disp(sprintf(['The two structs concern different grids:\n',...
             '     grid 1: (xl,yl)=(%6.3f,%6.3f), (dx,dy)=(%6.3f,%6.3f)\n', ...
             '     grid 2: (xl,yl)=(%6.3f,%6.3f), (dx,dy)=(%6.3f,%6.3f)\n'], ...
             grd1, grd2));
   return
end

% Copy the administration from sol1

dif.d_digit = sol1.d_digit;
dif.mx      = sol1.mx;
dif.my      = sol1.my;
dif.xl      = sol1.xl;
dif.yl      = sol1.yl;
dif.dx      = sol1.dx;
dif.dy      = sol1.dy;
dif.x       = sol1.x;
dif.y       = sol1.y;

dif.kincns.t_digit = sol1.kincns.t_digit;
dif.kincns.chi     = sol1.kincns.chi;
dif.kincns.dq      = sol1.kincns.dq;
dif.kincns.veloc   = sol1.kincns.veloc;

dif.config  = sol1.config;
dif.h_digit = sol1.h_digit;
dif.mater   = sol1.mater;
dif.fric    = sol1.fric;

dif.x_offset  = sol1.x_offset;
dif.y_offset  = sol1.y_offset;

% copy element division from sol1 and encode the difference

dif.eldiv = sol1.eldiv + 10*(sol1.eldiv-sol2.eldiv);

% determine differences of solution variables

dif.h      = sol1.h      - sol2.h;
sol.mu     = sol1.mu     - sol2.mu;
dif.pn     = sol1.pn     - sol2.pn;
dif.px     = sol1.px     - sol2.px;
dif.py     = sol1.py     - sol2.py;
dif.un     = sol1.un     - sol2.un;
dif.ux     = sol1.ux     - sol2.ux;
dif.uy     = sol1.uy     - sol2.uy;
dif.srel   = sol1.srel   - sol2.srel;
dif.shft   = sol1.shft   - sol2.shft;
dif.trcbnd = sol1.trcbnd - sol2.trcbnd;
if (isfield(sol1,'taucrt') & isfield(sol2,'taucrt'))
   dif.taucrt = sol1.taucrt - sol2.taucrt;
   dif.uplsx  = sol1.uplsx  - sol2.uplsx;
   dif.uplsy  = sol1.uplsy  - sol2.uplsy;
end
if (isfield(sol1,'temp1') & isfield(sol2,'temp1'))
   dif.temp1 = sol1.temp1 - sol2.temp1;
   dif.temp2 = sol1.temp2 - sol2.temp2;
end

% determine the difference in magnitude pt

dif.pt     = sqrt(sol1.px.^2 + sol1.py.^2) - sqrt(sol2.px.^2 + sol2.py.^2);
