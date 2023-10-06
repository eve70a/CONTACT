
function [ dif ] = diffstrs(sol1, sol2)

%
% function [ dif ] = diffstrs(sol1, sol2)
%
% Compute the difference of the results for two subsurface stress
% calculations sol1 and sol2.
%
% dif         - output struct "sol1 - sol2"
% sol1, sol2  - structs with subsurface stresses as defined by loadstrs.
%

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

dif = [];

% Check equality of the grids w.r.t. subsurface points used

if (sol1.nx~=sol2.nx)
   disp(sprintf('The two structs concern different subsurface points (nx).'));
   return
end
if (sol1.ny~=sol2.ny)
   disp(sprintf('The two structs concern different subsurface points (ny).'));
   return
end
if (sol1.nz~=sol2.nz)
   disp(sprintf('The two structs concern different subsurface points (nz).'));
   return
end

%  check x,y,z-coordinates of block

difc = max(abs(sol1.x-sol2.x));
if (difc>1e-10) 
   disp(sprintf('Warning: the x-coordinates differ by at most %f;\n         continuing with coords of first case', difc));
end
difc = max(abs(sol1.y-sol2.y));
if (difc>1e-10) 
   disp(sprintf('Warning: the y-coordinates differ by at most %f;\n         continuing with coords of first case', difc));
end
difc = max(abs(sol1.z-sol2.z));
if (difc>1e-10) 
   disp(sprintf('Warning: the z-coordinates differ by at most %f;\n         continuing with coords of first case', difc));
end

%  copy administration from sol1

dif.nx      = sol1.nx;
dif.ny      = sol1.ny;
dif.nz      = sol1.nz;
dif.npoints = sol1.npoints;
dif.x       = sol1.x;
dif.y       = sol1.y;
dif.z       = sol1.z;

%  compute differences of displacements and stress invariants

dif.ux = sol1.ux - sol2.ux;
dif.uy = sol1.uy - sol2.uy;
dif.uz = sol1.uz - sol2.uz;

dif.sighyd = sol1.sighyd - sol2.sighyd;
dif.sigvm  = sqrt(sol1.sigvm.^2 - sol2.sigvm.^2);

if (isfield(sol1,'sigxx') & isfield(sol2,'sigxx'))
   dif.sigxx = sol1.sigxx - sol2.sigxx;
   dif.sigxy = sol1.sigxy - sol2.sigxy;
   dif.sigxz = sol1.sigxz - sol2.sigxz;
   dif.sigyy = sol1.sigyy - sol2.sigyy;
   dif.sigyz = sol1.sigyz - sol2.sigyz;
   dif.sigzz = sol1.sigzz - sol2.sigzz;
end

