function [ slcs ] = read_slices(fname, mirror_y, mirror_z, scale_yz, idebug)

% function [ slcs ] = read_slices(fname, [mirror_y], [mirror_z], [scale_yz], [idebug])
%
% Lower-level routine for reading wheel/rail profiles in slices format.
% mirror_y    - -1 or 0 for no, 1 for yes
% mirror_z    - -1 for no, 0 for automatic, 1 for yes
% scale_yz    - scale-factor to convert profile data to mm, e.g. 1000. for data in meters

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

if (nargin<2 | isempty(mirror_y))
   mirror_y = 0;
end
if (nargin<3 | isempty(mirror_z))
   mirror_z = 0;
end
if (nargin<4 | isempty(scale_yz))
   scale_yz = 1;
end
if (nargin<5 | isempty(idebug))
   idebug = 1;
end
dsmax_2d = 0.5;

show_fig = [ ];
fig_ofs  = 4;
% figs 1 - 4: resampling of slices
% figs 5 - 6: slices before/after interruption

% using '.slcs' for rails and '.slcw' for wheels

[~,~,ext] = fileparts(fname);
is_wheel  = strcmp(lower(ext), '.slcw');

if (is_wheel)           % 1st surface parameter u == theta or u == s along track curve
   nam_udir = 'theta';
else
   nam_udir = 'U';
end

% Initialize output structure

slcs = struct;
slcs.slc_file = fname;
slcs.is_wheel = is_wheel;
slcs.nslc     = 0;
slcs.u_offset = 0;
slcs.u_scale  = 1;
slcs.u        = [];
slcs.fnames   = [];

% read file contents

if (~exist(fname,'file'))
   disp(['ERROR: cannot find file ',fname]);
   return
end

f = fopen(fname,'r'); iline  = 0;
[slcs, iline, ierror] = read_slices_header(slcs, nam_udir, f, iline, idebug);
if (ierror==0)
   [slcs, iline, ierror] = read_slices_fnames(slcs, nam_udir, f, iline, idebug);
end
if (ierror==0)
   [slcs, iline, ierror] = read_feature_info(slcs, nam_udir, f, iline, idebug);
end
fclose(f);
if (ierror), return; end

% apply u-offset and scaling

slcs.u = slcs.u_scale * (slcs.u + slcs.u_offset);

% read each of the profiles

[slcs, ierror] = read_slice_profiles( fname, slcs, mirror_y, mirror_z, scale_yz, idebug );

% resample all slices to the same number of points

if (ierror==0)
   [slcs, ierror] = resample_slices( slcs, dsmax_2d, idebug, show_fig, fig_ofs );
end

% add 2d spline representation 
% rails:  parameter (u,v) --> cartesian (x,y,z) with x==u.
% wheels: parameter (u,v) --> cylindrical (th,y,dr) with th==u, basic interval [-pi,pi) + wrap-around
% TODO: smoothing in x/theta-direction

if (ierror==0 & exist('make_2dspline'))
   use_approx = (slcs.u_intpol==2); use_insert = 1; use_cylindr = is_wheel; idebug = 0;
   % disp('slcs.ysurf:')
   % disp(slcs.ysurf)
   slcs.spl2d = make_2dspline(slcs.u, slcs.vj, [], slcs.ysurf, slcs.zsurf, slcs.mask_j, ...
                                                               use_approx, use_insert, use_cylindr, idebug);
   if (1==1)
      % add second spline surface for distance computation
      if (~use_approx)
         slcs.appx = make_2dspline(slcs.u, slcs.vj, [], slcs.ysurf, slcs.zsurf, slcs.mask_j, ...
                                                               1, use_insert, use_cylindr, idebug);
      else
         slcs.intp = make_2dspline(slcs.u, slcs.vj, [], slcs.ysurf, slcs.zsurf, slcs.mask_j, ...
                                                               0, use_insert, use_cylindr, idebug);
      end
   end
end

end % function read_slices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ slcs, iline, ierror ] = read_slices_header(slcs, nam_udir, f, iline, idebug)

% function [ slcs, iline ] = read_slices_header(slcs, nam_udir, f, iline, [idebug])
%
% Lower-level routine for reading contents of slcs-file 

   ierror    = 0;
   in_header = 1;

   % get next line from input

   while(~feof(f) & in_header<=4)

      s  = strtrim(fgets(f));
      iline = iline + 1;
      if (idebug>=3)
         disp(sprintf('%4d: %s',iline,s));
      end

      % remove comments, everything after %

      ix = strfind(s, '%'); if (~isempty(ix)), s=s(1:ix-1); end

      % process line

      if (~isempty(s))

         if (in_header==1)

            % 1st line: [S,TH]_OFFSET, [S,TH]_SCALE

            [tmp, nval] = sscanf(s, '%f %f');
            if (nval<2)
               disp(sprintf('ERROR: first line should have two values, %s_OFFSET, %s_SCALE, obtained %d',...
                        nam_udir, nam_udir, nval));
               ierror = 1;
            else
               slcs.u_offset = tmp(1);
               slcs.u_scale  = tmp(2);
            end
            in_header = in_header + 1;

         elseif (in_header==2)

            % 2nd line: NSLC

            slcs.nslc = sscanf(s, '%f');
            in_header = in_header + 1;

         elseif (in_header==3)

            % 3rd line: NFEAT, NKINK, NACCEL

            [tmp, nval] = sscanf(s, '%f %f %f');
            if (nval<3), tmp(3) = 0; end

            slcs.nfeat  = tmp(1);
            slcs.nkink  = tmp(2);
            slcs.naccel = tmp(3);
            in_header = in_header + 1;

         elseif (in_header==4)

            % 4th line: [U,TH]_INTPOL

            slcs.u_intpol = sscanf(s, '%f');
            in_header = in_header + 1;
            if (~any(slcs.u_intpol==[1 2]))
               disp(sprintf('ERROR: %s_INTPOL should be 1 or 2', nam_udir));
               disp(['input: ', s]);
               ierror = 4;
            end

         end

      end % ~isempty(s)
   end % while

   if (idebug>=-2)
      disp(sprintf('Slices-file "%s": %d slices, %d features, %d parts, %d kinks, %d accel', ...
                slcs.slc_file, slcs.nslc, slcs.nfeat, max(0,slcs.nfeat-1), slcs.nkink, slcs.naccel));
   end

end % function read_slices_header

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ slcs, iline, nerror ] = read_slices_fnames(slcs, nam_udir, f, iline, idebug)

% function [ slcs, iline, nerror ] = read_slices_fnames(slcs, nam_udir, f, iline, [idebug])
%
% Lower-level routine for reading filenames from slcs-file 

   nslc   = 0; % #slices read from file
   nerror = 0;

   while(~feof(f) & nslc<slcs.nslc)

      s  = strtrim(fgets(f));
      iline = iline + 1;
      if (idebug>=3)
         disp(sprintf('%4d: %s',iline,s));
      end

      % remove comments, everything after %

      ix = strfind(s, '%'); if (~isempty(ix)), s=s(1:ix-1); end

      if (~isempty(s))

         % NSLC lines: S/TH_SLC, 'R/WFNAME' or S/TH_SLC, "R/WFNAME"

         ix = findstr(s, '''');
         if (length(ix)==0)
            ix = findstr(s, '"');
         end

         if (length(ix)==2)
            nslc = nslc + 1;
            val            = sscanf(s(1:ix(1)-1), '%f');
            if (length(val)==1)
               slcs.u(nslc,1) = val;
            end
            nam            = s(ix(1)+1:ix(2)-1);
            if (~isempty(nam))
               slcs.fnames    = strvcat( slcs.fnames, nam );
            end
         end
         if (length(ix)~=2 | length(val)~=1 | isempty(nam))
            disp(sprintf('ERROR: could not interpret line %d:',iline));
            disp(sprintf('       "%s"',s));
            disp(sprintf('       line should provide %s_SLC, ''RFNAME''', nam_udir));
            nerror = nerror + 1;
         end

      end % ~isempty(s)
   end % while (~feof & islc<nslc)

   if (nslc<slcs.nslc)
      disp(sprintf('ERROR: got %d slices, expected NSLC = %d', nslc, slcs.nslc));
      nerror = nerror + 1;
   end

   % check that u-positions are in strictly increasing order, du >= 0.001 mm or >= 0.001 rad

   tiny_du = 0.001 / slcs.u_scale;
   du = diff(slcs.u);
   if (any(du<tiny_du))
      ix = find(du<tiny_du, 1, 'first');
      disp(sprintf('ERROR: slice %s-positions should be strictly increasing.', nam_udir));
      disp(sprintf('       islc= %d, %d: %s=%8.3f, %8.3f', ix, ix+1, nam_udir, slcs.u(ix), slcs.u(ix+1)));
      nerror = nerror + 1;
   end

   % for wheels, use some overlap < -pi and > pi to mitigate effect of boundary conditions
   %             require |u_i| \in [-1.5*pi, 1.5*pi], either \subset [-pi,pi) or [-pi,pi) \subset [u_i]

   if (slcs.is_wheel)
      th0 = slcs.u(1); th1 = slcs.u(end);
      if (th0<-1.5*pi | th1>1.5*pi)
         disp(sprintf('ERROR: wheel %s-positions should lie within [-1.5pi, 1.5pi].', nam_udir));
         disp(sprintf('       islc= %d, %d: %s=%8.3f, %8.3f', 1, slcs.nslc, nam_udir, th0, th1));
         nerror = nerror + 1;
      end
      eps_pi = 1e-4;
      if ((th0<-pi-eps_pi & th1<pi-eps_pi) | (th0>-pi+eps_pi & th1>pi+eps_pi))
         disp(sprintf('ERROR: wheel %s-positions should fully encompass [-pi,pi]', nam_udir));
         disp(sprintf('       or lie fully within [-pi,pi], got %s=[%8.3f,%8.3f]', nam_udir, th0, th1));
         nerror = nerror + 1;
      end
   end

end % function read_slices_fnames

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ slcs, iline, nerror ] = read_feature_info(slcs, nam_udir, f, iline, idebug)

% function [ slcs, iline, nerror ] = read_feature_info(slcs, nam_udir, f, iline, [idebug])
%
% Lower-level routine for reading feature information from a slcs-file 

   % TODO: add 2 lines: P_KINK, P_ACCEL

   nerror = 0;

   % NFEAT <= 1: using default feature information

   if (slcs.nfeat<=1)

      slcs.s_feat = ones(slcs.nslc,1) * [0, 1e6];

   else

      % next NSLC lines: u_slc, s_f, f = 1 : nfeat

      nslc_p = 0; % #slices with features read from file
      slcs.s_feat = zeros(slcs.nslc, slcs.nfeat);

      while(~feof(f) & nslc_p<slcs.nslc)

         % get next line from input

         s  = strtrim(fgets(f));
         iline = iline + 1;
         if (idebug>=3)
            disp(sprintf('%4d: %s',iline,s));
         end

         % remove comments, everything after %

         ix = strfind(s, '%'); if (~isempty(ix)), s=s(1:ix-1); end

         % process line

         if (~isempty(s))

            nslc_p  = nslc_p + 1;
            is_ok  = 1;

            [tmp, nval] = sscanf(s, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
   
            if (nval>=1)
               u_slc = tmp(1);
               i_slc = find(abs(slcs.u-u_slc)<1e-3);

               if (length(i_slc)<1)
                  is_ok = 0;
                  disp(sprintf('ERROR: could not find slice at %s_slc = %4.2f in list of slices', ...
                                                                                        nam_udir, u_slc));
                  nerror = nerror + 1;
               elseif (length(i_slc)>1)
                  is_ok = 0;
                  disp(sprintf('ERROR: multiple slices with %s_slc = %4.2f in list of slices', ...
                                                                                        nam_udir, u_slc));
                  if (~isempty(i_slc))
                     disp(sprintf('       possible slices i = %d %d %d %d %d', i_slc));
                  end
                  nerror = nerror + 1;
               elseif (any(slcs.s_feat(i_slc,:)))
                  is_ok = 0;
                  disp(sprintf('ERROR: repeated definition of features for slice with %s_slc = %4.2f', ...
                                                                                        nam_udir, u_slc));
                  nerror = nerror + 1;
               end
            end

            nbreak = nval - 1;

            if (is_ok & nbreak~=slcs.nfeat)
               disp(sprintf('ERROR: incorrect number of feature positions for slice at %s_slc = %4.2f', ...
                                                                                        nam_udir, u_slc));
               disp(sprintf('       obtained %d values, expecting n_feat = %d', nbreak, slcs.nfeat));
               is_ok = 0;
               nerror = nerror + 1;
            end

            if (is_ok)
               slcs.s_feat(i_slc, 1:nbreak) = tmp(2:nbreak+1);
            end

         end % ~isempty(s)
      end % while (~feof)

      if (nslc_p<slcs.nslc)
         disp(sprintf('ERROR: got feature information for %d slices, expected NSLC = %d', nslc_p, slcs.nslc));
         nerror = nerror + 1;
      end
   end % nfeat<=1

end % function read_feature_info

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ slcs, ierror ] = read_slice_profiles( fname, slcs, mirror_y, mirror_z, scale_yz, idebug )

% function [ slcs, ierror ] = read_slice_profiles( fname, slcs, mirror_y, mirror_z, scale_yz, idebug )
%
% read each of the per-slice profiles

ierror = 0;

% determine relative folder where slices-file resides

filepath = fileparts(fname);

nslc = slcs.nslc;

for is = 1 : nslc
   if (idebug>=1 | (idebug>=-1 & is==1) | (idebug>=-1 & is==nslc))
      disp(sprintf('Slice %3d: read file "%s"', is, strtrim(slcs.fnames(is,:)) ));
   end

   slc_file = fullfile(filepath, strtrim(slcs.fnames(is,:)));
   prf = read_profile(slc_file, slcs.is_wheel, mirror_y, mirror_z, scale_yz, [], idebug-1);

   % in case of an error (missing file), prf will be a struct with empty members

   if (isempty(prf.ProfileY))
      ierror = 1;
   else
      prf.npnt = length(prf.ProfileY);
      prf.ltot = prf.ProfileS(end) - prf.ProfileS(1);
   end

   if (is>1 & ~isempty(prf.ProfileY))
 
      % set empty for fields present in slcs.prf / missing in new prf

      slcs_fields = fieldnames(slcs.prf(1));
      for i = 1:length(slcs_fields)
         if (~isfield(prf, slcs_fields{i}))
            % disp(sprintf('prf(%d) does not have field "%s", left empty', is, slcs_fields{i}));
            prf = setfield(prf, slcs_fields{i}, []);
         end
      end

      % set empty for fields present in prf / missing in previous slcs.prf

      prf_fields = fieldnames(prf);
      for i = 1:length(prf_fields)
         if (~isfield(slcs.prf(1), prf_fields{i}))
            % disp(sprintf('prf(%d) introduces field "%s", empty for previous slices', is, prf_fields{i}));
            slcs.prf = setfield(slcs.prf, prf_fields{i}, []);
         end
      end
   end % is>1 & prf

   if (~isempty(prf.ProfileY))
      % disp(sprintf('Adding prf(%2d) = "%s"', is, prf.Fname));
      slcs.prf(is,1) = prf;
   end
end % for is

end % function read_slice_profiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

