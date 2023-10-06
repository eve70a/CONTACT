
function [ B1, B2, B3, B4, B5, B6, B7, B8 ] = eval_bspline_basisfnc( tj, si, k, idebug )

% function [ B1, B2, B3, B4, B5, B6, B7, B8 ] = eval_bspline_basisfnc( tj, si, [k], [idebug] )
%
% evaluate B-spline basis functions B_j,k for knot-vector tj at sample locations si
%  k = spline order, default k=4 for cubic splines
%
% B1 = basis functions for k=1 using knot-vector tj  --  size (nout, nknot-1)
% B2 = basis functions for k=2 using knot-vector tj  --  size (nout, nknot-2)
% B3 = basis functions for k=3 using knot-vector tj  --  size (nout, nknot-3)
% B4 = basis functions for k=4 using knot-vector tj  --  size (nout, nknot-4)

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

if (nargin<3 | isempty(k))
   k      = 4;
end
if (nargin<4 | isempty(idebug))
   idebug = 1;
end
nknot = length(tj);
npnt  = length(si);

% determine value of each basis-function at each measurement location si

if (size(si,2)>size(si,1)), si = si'; end

% check for values si outside basic interval [tj(4),tj(nknot-k+1)]

if (idebug>=1)
   n0 = nnz(si < tj(k));
   n1 = nnz(si > tj(nknot-k+1));
   if (n0+n1>0)
      disp(sprintf('Warning: there are %d u-positions before and %d after basic interval t=[%3.1f,%3.1f]', ...
                n0, n1, tj(k), tj(nknot-k+1)));
   end
end

for ik = 1 : k

   % initialize sparse array

   Bk = sparse(npnt, nknot-ik);

   % average step sizes over ik-1 adjacent intervals

   dtj_k   = (tj(ik:end) - tj(1:nknot-ik+1)) / (ik-1);

   % inverse average step sizes, zero at zero step size

   dtj_inv = zeros(size(dtj_k)); ix  = find(dtj_k>0); dtj_inv(ix) = 1./ dtj_k(ix);

   % determine index j for start of last 'true' interval. 
   % j_last == nknot-k in case the end-knot is added k-1 times

   j_last = find(tj<tj(end), 1, 'last');

   if (ik==1) % ik=1: piecewise constant B_j,1

      for j = 1 : nknot-ik
         if (ik==1 & j==j_last)
            Bk(:,j)  = (si >= tj(j)  & si <= tj(j+1));
         else
            Bk(:,j)  = (si >= tj(j)  & si <  tj(j+1));
         end
      end

   else

      % ik=2: piecewise linear B_j,2

      % disp([ik nknot-ik length(dtj_k) length(dtj_inv)])

      for j = 1 : nknot-ik
         % disp([ik, j, length(tj) length(dtj_inv) size(Bp)])
         Bk(:,j) = ( (si - tj(j)) .* Bp(:,j) * dtj_inv(j) + ...
                   (tj(j+ik) - si) .* Bp(:,j+1) * dtj_inv(j+1) ) / (ik-1);
      end

   end

   % store result in output B1, B2, etc

   eval(sprintf('B%d = Bk;', ik));

   % copy to Bp for previous B_{ik-1} in next iteration

   Bp = Bk;

end

end % function eval_bspline_basisfnc

