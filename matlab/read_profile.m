function [ p ] = read_profile(fname, is_wheel, mirror_y, mirror_z, scale_yz, rgt_side, idebug, make_plot)

% function [ p ] = read_profile(fname, [is_wheel], [mirror_y], [mirror_z], [scale_yz], [rgt_side], ...
%                                                                          [idebug], [make_plot])
%
% Main routine for reading wheel/rail profiles. 
% Switches between Slices, SIMPACK, MiniProf and Vampire files using the filename extension.
% Uses modify_profiles to apply automatic corrections.
%
% is_wheel    -  0 for rail, 1 for wheel
% mirror_y    - -1 or 0 for no, 1 for yes
% mirror_z    - -1 for no, 0 for automatic, 1 for yes
% scale_yz    - scale-factor to convert to mm, e.g. 1000. for data in meters
% rgt_side    - Vampire format: select right side (1, default) or left side (0)
% make_plot   - <=0 for no, fig.number for yes

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

[~,~,ext] = fileparts(fname);

if (nargin<2 | isempty(is_wheel))
   is_wheel = (strcmp(lower(ext), '.prw') | strcmp(lower(ext), '.whe') | ...
                                                   strcmp(lower(ext), '.whl')); % default: rail
end
if (nargin<3 | isempty(mirror_y))
   mirror_y = 0;
end
if (nargin<4 | isempty(mirror_z))
   mirror_z = 0;
end
if (nargin<5 | isempty(scale_yz))
   scale_yz = 1;
end
if (nargin<6 | isempty(rgt_side))
   rgt_side = 1;
end
if (nargin<7 | isempty(idebug))
   idebug = 0;
end
if (nargin<8 | isempty(make_plot))
   make_plot = 0;
end

% switch to appropriate routine based on filename extension

is_slices = 0; 
is_simpack = 0; 
is_vampire = 0; 

if (any(strcmp(lower(ext), {'.slcs'}) ))

   p = read_slices(fname, mirror_y, mirror_z, scale_yz, idebug);
   is_slices  = 1;

elseif (any(strcmp(lower(ext), {'.prr', '.prw'}) ))

   p = read_simpack(fname, idebug, make_plot);
   is_simpack = 1;

elseif (any(strcmp(lower(ext), {'.rai', '.whe'}) ))

   p = read_vampire(fname, rgt_side, idebug, make_plot);
   is_vampire = 1;

elseif (any(strcmp(lower(ext), {'.ban', '.whl'}) ))

   p = read_miniprof(fname, idebug, make_plot);

else

   if (idebug>=1)
      disp(sprintf('Unknown file extension "%s", trying 2-column format.', ext))
   end
   p = read_miniprof(fname, idebug, make_plot);

end

if (isempty(p) | (~is_slices & ~isfield(p, 'ProfileY')) | (~is_slices & isempty(p.ProfileY)))
   if (idebug>=1)
      disp('No profile, returning');
   end
   return
end

p.is_wheel = is_wheel;

% make automatic adjustments: mirroring, swapping

if (~is_slices)

   y = p.ProfileY; ny = length(y);

   if (is_simpack)
      delete_foldback = 0;
   else
      delete_foldback = 0;
   end

   if (mirror_z==0) % automatic mode
      if (is_simpack)
         mirror_z = -1;
      else
         if (is_wheel)
            % mirror z if lowest z-value occurs in interior (flange)
            [minz,ix0]  = min(p.ProfileZ);
            [maxz,ix1]  = max(p.ProfileZ);
            if (ix0<0.05*ny | ix0>0.95*ny)
               mirror_z = -1;  % lowest z at start or end
            elseif (ix1<0.05*ny | ix1>0.95*ny)
               mirror_z =  1;  % highest z at start or end
            else
                  ix_m = round(mean([ix0 ix1]));
               z_m = mean([minz, maxz]);
               mirror_z = (p.ProfileZ(ix_m) > z_m); % mirror when profile z > straight line
            end
         else
            % mirror z if lowest z-value does not occur on the tread of the rail
            [~,ix]   = min(p.ProfileZ);
            mirror_z = (ix<0.1*ny | ix>0.9*ny);
         end
      end
   end

   if (is_wheel)
      % reverse order if first y-value < last y-value, after mirroring
      reverse_order = ( (mirror_y<=0 & y(1)<y(end)) | (mirror_y>=1 & y(1)>y(end)) );
   else
      % reverse order if last y-value < first y-value, after mirroring
      reverse_order = ( (mirror_y<=0 & y(end)<y(1)) | (mirror_y>=1 & y(end)>y(1)) );
   end

   if (scale_yz~=1)
      point_dist_min = 1e-4 / scale_yz;
   elseif (max(abs(p.ProfileY))>1) % data in [mm]
      point_dist_min = 1e-4; % mm
   else
      point_dist_min = 1e-7; % m
   end

   p = modify_profile(p, is_wheel, mirror_y, mirror_z, reverse_order, point_dist_min, ...
                                             delete_foldback, scale_yz, make_plot, idebug);

   % check for big gaps in the data

   dy = diff(p.ProfileY);
   dz = diff(p.ProfileZ);
   ds = sqrt(dy.^2 + dz.^2);

   if (idebug>=1 & max(ds) > 20*mean(ds))
      disp(sprintf('Warning: max(ds)=%6.2f >> avg(ds)=%6.2f, missing data?', max(ds), mean(ds)));
   end

   % compute s-coordinate along profile

   p.ProfileS = make_arclength(p.ProfileY, p.ProfileZ);
end % not slices

p.Fname = fname;

end % function read_profile

