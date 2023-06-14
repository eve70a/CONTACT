
function [ u_out, found ] = solve_cubic_eq(ppcoef, useg, yev, idebug, iout, jout, jseg)

% function [ u_out, found ] = solve_cubic_eq(ppcoef, useg, yev, idebug, iout, jout, jseg)
%
% solve 3rd degree equation f(u) = y, returning one root in interval [useg(1),useg(2)]
% ppcoef = [a0, a1, a2, a3],  d = a0-yev, c = a1, b = a2, a = a3

   if (nargin<4 | isempty(idebug))
      idebug = 0;
   end
   if (nargin<5 | isempty(iout))
      iout   = 0;
   end
   if (nargin<6 | isempty(jout))
      jout   = 0;
   end
   if (nargin<7 | isempty(jseg))
      jseg   = 0;
   end
   if (jseg==-335), idebug = 5; end

   if (1==1)
      [ u_out1, found ] = solve_cubic_cardano(ppcoef, useg, yev, idebug, iout, jout, jseg);
      [ u_out2, found ] = solve_cubic_newton(ppcoef, useg, yev, idebug, iout, jout, jseg);
      if (abs(u_out1-u_out2)>1e-6)
         disp(sprintf('  jseg = %d: different solutions u_newt = %8.4f, u_card = %8.4f, diff = %3.1e', ...
                                                                 jseg, u_out2, u_out1, abs(u_out2-u_out1)));
      end
      u_out = u_out1;
   else
      [ u_out, found ] = solve_cubic_cardano(ppcoef, useg, yev, idebug, iout, jout, jseg);
      % [ u_out, found ] = solve_cubic_newton(ppcoef, useg, yev, idebug, iout, jout, jseg);
   end

end % function solve_cubic_eq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ u_out, found ] = solve_cubic_newton(ppcoef, useg, yev, idebug, iout, jout, jseg)

% function [ u_out, found ] = solve_cubic_newton(ppcoef, useg, yev, idebug, iout, jout, jseg)
%
% solve 3rd degree equation f(u) = y, returning one root in interval [useg(1),useg(2)]

   % solve f(u) = a3 ul^3 + a2 ul^2 + a1 ul + a0 - yev = 0, with ul the local coordinate in the segment
   % using Newton-Raphson, cannot jump over extremal values, needs appropriate u^0

   a0 = ppcoef(1);
   a1 = ppcoef(2);
   a2 = ppcoef(3);
   a3 = ppcoef(4);

   if (idebug>=3)
      ul = [-100:0.01:20];
      u4 = useg(2) - useg(1);
      f_u = a3 * ul.^3 +   a2 * ul.^2 + a1 * ul + a0 - yev;

      figure(12); clf; hold on;
      plot(ul, f_u);
      v = axis;
      plot([0 0 NaN u4 u4], [v(3:4) NaN v(3:4)], '--');
      plot(v(1:2), 0*[1 1], '--');
      grid on;
      xlabel(sprintf('u_l in seg %d', jseg))
      ylabel('y_l = f(u_l)')
      title(sprintf('Newton for (i,j)=(%d,%d)',iout,jout));
      l = legend(sprintf('cubic, jseg=%d',jseg));
      % set(l, 'autoupdate','off');
   end

   tiny_a2 = 1e-10;
   tiny_a3 = 1e-10;

   % check for local extremal values: 
   % zeros of quadratic f'(ul) = 3a ul^2 + 2b ul + c = 0

   discr = (2*a2)^2 - 4 * (3*a3) * a1;

   if (idebug>=3)
      disp(sprintf('  jseg = %d: a0 = %6.3e, a1 = %6.3e, a2 = %6.3e, a3 = %6.3e, D = %6.3e', ...
                                                                        jseg, a0, a1, a2, a3, discr));
   end

   u0 = 0;
   u4 = useg(2) - useg(1);

   if (abs(a2)<tiny_a2 & abs(a3)<tiny_a3)

      % segment function is linear, within roundoff precision; use initial estimate ul = 0

      if (idebug>=5), disp('  case 1'); end
      ul    = u0;
      f_u   = a0 - yev;
      df_du = a1;

   elseif (discr<=0)

      % no local extremal values: use initial estimate ul = 0

      if (idebug>=5), disp('  case 2'); end
      ul    = u0;
      f_u   = a0 - yev;
      df_du = a1;

   elseif (abs(a3)<tiny_a3)

      % segment function is quadratic, within roundoff precision

      if (idebug>=5), disp('  case 3'); end
      ul    = u0;
      f_u   = a0 - yev;
      df_du = a1;

      % check extremal value, should have opposite sign to a2

      u1    = -a1 / (2 * (2*a2));  % extremum
      f1    = a2 * u1^2 + a1 * u1 + a0 - yev;
      if (a2*f1 > 0)
         % disp('case 3: no solution')
         ul = NaN;
      end

   else

      u1 = ( -(2*a2) - sqrt(discr) ) / (2 * (3*a3));
      u2 = ( -(2*a2)               ) / (2 * (3*a3));  % inflection point
      u3 = ( -(2*a2) + sqrt(discr) ) / (2 * (3*a3));
      % swap u1,u3 when u1>u3
      if (u1>u3), tmp=u1; u1=u3; u3=tmp; end

      f0 = a0 - yev;
      f1 = a3 * u1^3 + a2 * u1^2 + a1 * u1 + a0 - yev;
      f2 = a3 * u2^3 + a2 * u2^2 + a1 * u2 + a0 - yev;
      f3 = a3 * u3^3 + a2 * u3^2 + a1 * u3 + a0 - yev;
      f4 = a3 * u4^3 + a2 * u4^2 + a1 * u4 + a0 - yev;
   
      if (idebug>=3)
         v = axis; dy = (v(4)-v(3)) / 5;
         plot(u0*[1 1], f0+dy*[-1,1],'b--'); text(u0, f0-1.1*dy, 's_0');
         plot(u1*[1 1], f1+dy*[-1,1],'r--'); text(u1, f1-1.1*dy, 's_1');
         plot(u2      , f2          ,'o'  ); text(u2, f2-1.1*dy, 's_2');
         plot(u3*[1 1], f3+dy*[-1,1],'r--'); text(u3, f3-1.1*dy, 's_3');
         plot(u4*[1 1], f4+dy*[-1,1],'b--'); text(u4, f4-1.1*dy, 's_4');
         % disp([a2, discr, f0, f1, f3, f4])
      end

      % there can be one zero before u1, one between u1 and u3, and one after u3
      % select initial estimate: start-point or end-point of segment jseg

      i1_in_u04 = ((u4-u1)*(u1-u0) > 0);
      i2_in_u04 = ((u4-u3)*(u3-u0) > 0);

      if (f1*f3>0 & a3*f1>0)
         % one zero, to the left of u1
         if (u1<=u0) % no solution within [u0,u4]
            if (idebug>=5), disp('  case 4.1.a'); end
            ul = u1+(u1-u0)+(u1-u2); f_u = a3 * ul.^3 +   a2 * ul.^2 + a1 * ul + a0 - yev;
         else
            if (idebug>=5), disp('  case 4.1.b'); end
            ul = u0; f_u = f0;
         end
      elseif (f1*f3>0)
         % one zero, to the right of u3
         if (u3>=u4) % no solution within [u0,u4]
            if (idebug>=5), disp('  case 4.2.a'); end
            ul = u3+(u3-u4)+(u3-u2); f_u = a3 * ul.^3 +   a2 * ul.^2 + a1 * ul + a0 - yev;
         else
            if (idebug>=5), disp('  case 4.2.b'); end
            ul = u4; f_u = f4;
         end
      elseif (i1_in_u04 & i2_in_u04)
         % three zeros, middle one between u0 and u4
         if (idebug>=5), disp('  case 4.3'); end
         ul = (u1+u3) / 2;
         f_u = a3 * ul^3 + a2 * ul^2 + a1 * ul + a0 - yev;
      elseif ( (i1_in_u04 & f0*f1>0) | (i2_in_u04 & f0*f3>0))
         if (idebug>=5), disp('  case 4.4'); end
         % three zeros, right one between u0 and u4
         ul = u4; f_u = f4;
      else
         if (idebug>=5), disp('  case 4.5'); end
         % three zeros, left one between u0 and u4
         ul = u0; f_u = f0;
      end
      df_du = 3 * a3 * ul^2 + 2 * a2 * ul + a1;
   end

   % Newton-Raphson: f(u+du) = f(u) + du * f'(u) = 0   
   %                      -->  du = - f(u) / f'(u)

   iter  = 0; tol = 1e-8; maxit = 20;
   if (~isnan(ul) & idebug>=2)
      disp(sprintf(' It %2d: ul = %8.4f in [%8.4f,%8.4f], f(u) = %12.8f, f''(u) = %8.4f', ...
                        iter, ul, u0, u4, f_u, df_du))
      if (idebug>=3)
         plot(ul, f_u, '*')
         text(ul, f_u, num2str(iter))
      end
   end

   while(~isnan(ul) & abs(f_u)>=tol & iter<maxit)
      iter  = iter + 1;
      du    = - f_u / df_du;

      % restrict step esp. in first iterations, cf. doubling of step until bracket is found
      du_max = (u4 - u0) * 2^iter;
      if (du>du_max)
         du = du_max;
      elseif (du<-du_max)
         du = -du_max;
      end

      % prevent jumping over u2 in first iteration
      if (iter<=1 & discr>0 & exist('u2'))
         if (ul<u2)
            ul = min(u2, ul+du);
         else
            ul = max(u2, ul+du);
         end
      else
         ul = ul + du;
      end

      % evaluate function and derivative at new position
      f_u   =   a3 * ul^3 +   a2 * ul^2 + a1 * ul + a0 - yev;
      df_du = 3*a3 * ul^2 + 2*a2 * ul   + a1;
      if (idebug>=2)
         disp(sprintf(' It %2d: ul = %8.4f in [%8.4f,%8.4f], f(u) = %12.8f, f''(u) = %8.4f', ...
                        iter, ul, u0, u4, f_u, df_du))
      end
      if (idebug>=3)
         plot(ul, f_u, '*')
         text(ul, f_u, num2str(iter))
      end
   end

   if (isnan(ul))
      if (idebug>=1)
         disp(sprintf('solve_newton:  no solution for jout=%3d, jseg=%3d, yev=%8.3f', jout, jseg, yev));
      end
      found = 0;
   elseif (abs(f_u)>=tol)
      disp(sprintf('solve_newton:  no convergence for jout=%3d, jseg=%3d, yev=%8.3f', jout, jseg, yev));
      found = 0;
   elseif (ul<u0 | ul>u4)
      if (idebug>=1)
         disp(sprintf('solve_newton:  ul=%8.4f outside segment [%8.3f,%8.3f] for jout=%3d, yev=%8.3f', ...
                             ul, u0, u4, jout, yev));
      end
      found = 0;
   else
      found = 1;
   end

   u_out = useg(1) + ul;

end % function solve_cubic_newton

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ u_out, found ] = solve_cubic_cardano(ppcoef, useg, yev, idebug, iout, jout, jseg)

% function [ u_out, found ] = solve_cubic_cardano(ppcoef, useg, yev, idebug, iout, jout, jseg)
%
% solve 3rd degree equation f(u) = y, returning one root in interval [useg(1),useg(2)]

   % solve f(u) = a ul^3 + b ul^2 + c ul + d = 0,  d = a3-yev,
   % with ul the local coordinate in the segment, using Cardano's method.

   d = ppcoef(1) - yev;
   c = ppcoef(2);
   b = ppcoef(3);
   a = ppcoef(4);

   if (idebug>=4)
      ul = [-25:0.01:15];
      u4 = useg(2) - useg(1);
      f_u = a * ul.^3 + b * ul.^2 + c * ul + d;

      figure(12); clf; hold on;
      plot(ul, f_u);
      v = axis;
      plot([0 0 NaN u4 u4], [v(3:4) NaN v(3:4)], '--');
      plot(v(1:2), 0*[1 1], '--');
      grid on;
      xlabel(sprintf('u_l in seg %d', jseg))
      ylabel('y_l = f(u_l)')
      title(sprintf('Cardano for (i,j)=(%d,%d)',iout,jout));
      l = legend(sprintf('cubic, jseg=%d',jseg));
      % set(l, 'autoupdate','off');
   end

   if (idebug>=4)
      disp(sprintf('  jseg = %d: a = %6.3e, b = %6.3e, c = %6.3e, d = %6.3e', jseg, a, b, c, d));
   end

   tiny_a   = 1e-10;
   tiny_b   = 1e-10;
   tiny_cfc = 1e-20;
   u0     = 0;                  % interval [u0,u4] in local coordinates
   u4     = useg(2) - useg(1);

   if (abs(b)<tiny_b & abs(a)<tiny_a)

      % segment function is linear, within roundoff precision

      if (idebug>=6), disp('  case 1'); end
      ul    = -d / c;

   elseif (abs(a)<tiny_a)

      % segment function is quadratic, within roundoff precision

      if (idebug>=6), disp('  case 2'); end

      a = b; b = c; c = d;
      discr = b^2 - 4*a*c;
      if (discr<0)
         ul = NaN;
      else
         x1 = (-b - sqrt(discr)) / (2*a);
         x2 = (-b + sqrt(discr)) / (2*a);
         if (x1>=u0 & x1<=u4)
            ul = x1;
         else
            ul = x2;
         end
      end

   else

     % Cardano's method

     delta0 = b^2 - 3*a*c;
     delta1 = 2*b^3 - 9*a*b*c + 27*a^2*d;
     discr  = (4*delta0^3 -  delta1^2) / (27*a^2);

     if (abs(discr)>tiny_cfc)
        coefc  = ( ( delta1 + sqrt(delta1^2 - 4*delta0^3) ) / 2 )^(1/3);
        if (abs(coefc)<tiny_cfc)
           disp('Impossible coefc=0?')
           coefc  = ( ( delta1 - sqrt(delta1^2 - 4*delta0^3) ) / 2 )^(1/3);
        end
        ksi  = (-1 + sqrt(-3)) / 2;
     end

     if (abs(discr)<tiny_cfc & abs(delta0)<tiny_cfc)

        % triple root

        x1 = -b / (3*a);
        x2 = x1;
        x3 = x1;
        ul = x1;
        if (idebug>=3)
           disp('cardano: triple root');
        end

     elseif (abs(discr)<tiny_cfc)

        % double root + single root
        disp(discr)

        x1 = (4*a*b*c - 9*a^2*d - b^3) / (a*delta0); % single
        x2 = (9*a*d - b*c) / (2*delta0);             % double
        x3 = x2;
        if (x1>=u0 & x1<=u4)
           ul = x1;
        else
           ul = x2;
        end
        if (idebug>=3)
           disp('cardano: double root');
        end

     elseif (discr<0)

        % one real root + two non-real complex conjugate roots
        % coefc can be complex especially when delta1 < 0.

        x1 = -( b + ksi  *coefc + delta0 / (ksi  *coefc)) / (3*a);
        x2 = -( b + ksi^2*coefc + delta0 / (ksi^2*coefc)) / (3*a);
        x3 = -( b +       coefc + delta0 /        coefc ) / (3*a);

        % select root with smallest imaginary part (theoretically zero)

        x      = [x1, x2, x3];
        [~, i] = min(abs(imag(x)));
        ul = x(i);

        if (idebug>=3)
           disp('cardano: one real root & two non-real complex conjugate roots');
        end

     else

        % three distinct real roots

        x1 = -( b + ksi  *coefc + delta0 / (ksi  *coefc)) / (3*a);
        x2 = -( b + ksi^2*coefc + delta0 / (ksi^2*coefc)) / (3*a);
        x3 = -( b +       coefc + delta0 /        coefc ) / (3*a);

        x = [x1 x2 x3];
        if (max(abs(imag(x)))>1e-12*max(abs(real(x))))
           disp(sprintf('Impossible: nonzero imag.part (%3.1e)???', max(abs(imag(x))) ))
           disp(x)
        end
        x = real(x);

        % select a root in interval [u0,u4] or closest to interval [u0,u4]

        xm = (u0 + u4) / 2;
        [~,i] = min(abs( x-xm ));
        ul = x(i);

        if (idebug>=3)
           disp('cardano: three distinct real roots');
        end

     end

     f1     = a * x1.^3 + b * x1.^2 + c * x1 + d;
     f2     = a * x2.^3 + b * x2.^2 + c * x2 + d;
     f3     = a * x3.^3 + b * x3.^2 + c * x3 + d;

     if (idebug>=4)
        disp('solutions x1, x2, x3:')
        disp([x1, x2, x3; f1, f2, f3]);
     end

   end

   if (abs(imag(ul))>1e-10)
      disp('Impossible: complex ul???')
   end
   u_out = useg(1) + real(ul);

   if (isnan(ul))
      if (idebug>=3)
         disp(sprintf('  solve_cardano: no solution for jout=%3d, jseg=%3d, yev=%8.3f', jout, jseg, yev));
      end
      found = 0;
   elseif (ul<u0 | ul>u4)
      if (idebug>=3)
         disp(sprintf('  solve_cardano: u=%8.4f outside segment [%8.3f,%8.3f] for jout=%3d, yev=%8.3f', ...
                             real(ul), u0, u4, jout, yev));
      end
      found = 0;
   else
      found = 1;
   end

end % function solve_cubic_cardano
