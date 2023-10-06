
function [ blk1, blk2, blk3, blk4, blk5, blk6, blk7, blk8 ] = ...
                                              loadstrs( expnam, icase, ipatch )
%
% [ blk1, blk2, .. ] = loadstrs( expnam, [icase], [ipatch] )
%
% Read a single <expnam>.subs-file and return structures with the contents.
%
%   expnam = the experiment name used, i.e. the input-filename without
%            extension .inp.
%   icase  = the case-number to be read
%   ipatch = contact patch number, in case of module 1, wheel/rail contact
%   blk1   = structure with results of the computation for the first
%               block (grid) of subsurface stress points
%   blk2   = structure with results of the computation for the second
%               block (grid) of subsurface stress points
%   blk3, .. = consecutive blocks
%

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% initialize optional arguments

if (nargin<2 | isempty(icase))
   icase = 1;
end
if (nargin<3 | isempty(ipatch))
   ipatch = 1;
end

% locate and read the file

fname1=[expnam,'.subs'];
if (icase<=9999)
   fname2 = [expnam, sprintf('.%04d.subs', icase)];
   fname3 = [expnam, sprintf('.%04d%c.subs', icase, 96+ipatch)];
else
   fname2 = [expnam, sprintf('.%06d.subs', icase)];
   fname3 = [expnam, sprintf('.%06d%c.subs', icase, 96+ipatch)];
end

if (exist(fname1,'file'))
   all_data = load(fname1, '-ascii');
elseif (exist(fname2,'file'))
   all_data = load(fname2, '-ascii');
elseif (exist(fname3,'file'))
   all_data = load(fname3, '-ascii');
else
   disp(['ERROR: cannot find file ',fname1,' nor ',fname2,' or ',fname3])
   sol = [];
   return
end

% decompose the file contents into the constitutive blocks

if (size(all_data,2)==14)
   full_tens= 1;
elseif (size(all_data,2)==4)
   full_tens=-1;
else
   full_tens= 0;
end

nrows = size(all_data,1);
iblk  = 0;      % shorthand for sol.nblocks

irow  = 1;
while (irow <= nrows)
   iblk = iblk + 1;
   sol = [];    % sol == the current block being processed

   % first line of the block: nx, ny, nz
   nx = all_data(irow,1); sol.nx = nx;
   ny = all_data(irow,2); sol.ny = ny;
   nz = all_data(irow,3); sol.nz = nz;
   % disp(sprintf('block %d has %d points (%d x %d x %d)', iblk, nx*ny*nz, nx, ny, nz));

   % the data lines for the block:
   rng = irow + [1:nx*ny*nz];
   blk_data = all_data(rng,:);

   % sort data such that x runs fastest, then y, then z
   [tx, iperm] = sort(blk_data(:,1));
   blk_data = blk_data(iperm,:);
   [ty, iperm] = sort(blk_data(:,2));
   blk_data = blk_data(iperm,:);
   [tz, iperm] = sort(blk_data(:,3));
   blk_data = blk_data(iperm,:);

   % the data-lines are given with iz running fastest, ix slowest.
   sol.npoints = nx * ny * nz;
   sol.x       = reshape(blk_data(:,1), nx, ny, nz);
   sol.y       = reshape(blk_data(:,2), nx, ny, nz);
   sol.z       = reshape(blk_data(:,3), nx, ny, nz);
   if (full_tens==-1)
      sol.sigvm   = reshape(blk_data(:,4), nx, ny, nz);
   else
      sol.ux      = reshape(blk_data(:,4), nx, ny, nz);
      sol.uy      = reshape(blk_data(:,5), nx, ny, nz);
      sol.uz      = reshape(blk_data(:,6), nx, ny, nz);
      sol.sighyd  = reshape(blk_data(:,7), nx, ny, nz);
      sol.sigvm   = reshape(blk_data(:,8), nx, ny, nz);
   end

   if (full_tens==1)
      sol.sigxx = reshape(blk_data(:, 9), nx, ny, nz);
      sol.sigxy = reshape(blk_data(:,10), nx, ny, nz);
      sol.sigxz = reshape(blk_data(:,11), nx, ny, nz);
      sol.sigyy = reshape(blk_data(:,12), nx, ny, nz);
      sol.sigyz = reshape(blk_data(:,13), nx, ny, nz);
      sol.sigzz = reshape(blk_data(:,14), nx, ny, nz);

      if (1==1)
         % 3D calculation of principal stresses
         princ = zeros(nx*ny*nz,3);
         sigma = blk_data(:,9:14);
         for k=1:nx*ny*nz
            princ(k,:) = sort(eig([sigma(k,1), sigma(k,2), sigma(k,3); 
                                   sigma(k,2), sigma(k,4), sigma(k,5);
                                   sigma(k,3), sigma(k,5), sigma(k,6)]),'descend');
         end
         sol.sigma1 = reshape(princ(:,1), nx, ny, nz);
         sol.sigma2 = reshape(princ(:,2), nx, ny, nz);
         sol.sigma3 = reshape(princ(:,3), nx, ny, nz);
         sol.sigtr = sol.sigma1 - sol.sigma3;
      else
         % 2D principal stress calculation for Oxz: assuming plane strain
         princ = zeros(nx*ny*nz,2);
         sigma = blk_data(:,9:14);
         for k=1:nx*ny*nz
            princ(k,:) = eig([sigma(k,1), sigma(k,3);
                              sigma(k,3), sigma(k,6)]);
         end
         sol.sigma1 = reshape(princ(:,1), nx, ny, nz);
         sol.sigma2 = reshape(princ(:,2), nx, ny, nz);
         sol.sigma3 = 0.28*(sol.sigma1+sol.sigma2);
         sol.sigtr = sol.sigma1 - sol.sigma2;
      end
   end

   sol.x       = squeeze(sol.x(:,1,1));
   sol.y       = squeeze(sol.y(1,:,1)); sol.y = sol.y';
   sol.z       = squeeze(sol.z(1,1,:));

   % advance row-number to first line of next block
   irow = irow + 1 + nx*ny*nz;
   % disp(sprintf('next block (%d) starts at row %d', iblk+1, irow));

   % copy the results to the appropriate output argument
   if (iblk==1 | iblk<=nargout)
      eval(sprintf('blk%d = sol;',iblk));
   else
      % ignore struct sol if no more output arguments are available
      % write a warning upon reading/processing of the whole file
      if (irow>=nrows)
         disp(sprintf('WARNING: this .subs-file contains data for %d blocks of points.',iblk));
         disp(sprintf('         Only %d output arguments are provided, the remaining blocks are ignored.',max(1,nargout)));
      end
   end
end

if (iblk<nargout)
   disp(sprintf('WARNING: this .subs-file contains data for %d blocks of points.',iblk));
   disp(sprintf('         %d output arguments are provided, the last %d are not filled in.', nargout, nargout-iblk));
   for jblk= iblk+1 : nargout
      eval(sprintf('blk%d = [];',jblk));
   end
end
