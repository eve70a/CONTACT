function [ p ] = read_simpack(fname, idebug, make_plot)

% function [ p ] = read_simpack(fname, [idebug], [make_plot])
%
% Lower-level routine for reading wheel/rail profiles in SIMPACK format.
% Does not do automatic corrections like mirroring or reversing order.

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

if (nargin<2 | isempty(idebug))
   idebug = 1;
end
if (nargin<3 | isempty(make_plot))
   make_plot = 0;
end

known_fields_header = struct( ...
   'version',           'int', ...
   'type',              'int' ...
);

known_fields_spline = struct( ...
   'approx_smooth',     'real', ...
   'file',              'str', ...
   'file_mtime',        'int', ...
   'comment',           'str', ...
   'type',              'int', ...
   'point_dist_min',    'real', ...
   'shift_y',           'real', ...
   'shift_z',           'real', ...
   'rotate',            'real', ...
   'bound_y_min',       'real', ...
   'bound_y_max',       'real', ...
   'bound_z_min',       'real', ...
   'bound_z_max',       'real', ...
   'mirror_y',          'int', ...
   'mirror_z',          'int', ...
   'inversion',         'int', ...
   'units_len',         'str', ...
   'units_ang',         'str', ...
   'units_len_f',       'real', ...
   'units_ang_f',       'real' ...
);

% Initialize output structure

p = struct;
p.OrigData  = [0 0 0];
nrows_data  = 0;
p.ProfileY  = [];
p.ProfileZ  = [];
p.SplineWgt = [];

if (~exist(fname,'file'))
   disp(['ERROR: cannot find file ',fname]);
   return
end
f = fopen(fname,'r');
iline = 0;
in_header = 0;
in_spline = 0;
in_point = 0;

while(~feof(f))

   % get next line from input

   s  = strtrim(fgets(f));
   iline = iline + 1;
   if (idebug>=3)
      disp(sprintf('%4d: %s',iline,s));
   end

   % remove comments -- TODO: support comments in text strings

   ix = findstr('!', s);
   if (~isempty(ix))
      s = s(1:ix(1)-1);
   end

   % look for section start

   ix = findstr('.begin',s);
   if (~isempty(ix))
      section = strtrim(s(1:ix-1));
      if (strcmp(section,'header'))
         in_header = 1;
      elseif (strcmp(section,'spline'))
         in_spline = 1;
      elseif (strcmp(section,'point'))
         in_point = 1;
      else
         disp(['Unknown section: ',section]);
         return
      end
      continue % go read next line
   end

   % look for section end

   ix = findstr('.end',s);
   if (~isempty(ix))
      section = strtrim(s(1:ix-1));
      if (strcmp(section,'header'))
         in_header = 0;
      elseif (strcmp(section,'spline'))
         in_spline = 0;
      elseif (strcmp(section,'point'))
         in_point = 0;
      else
         disp(['Unknown section: ',section]);
         return
      end
      continue % go read next line
   end

   % look for an equals sign
   % if found: interpret keyword=value syntax

   ix = findstr('=',s);

   if (~isempty(ix))
      ix = ix(1);
      field_name = strtrim(s(1:ix-1));
      % replace dots in keyword by underscores
      jx = findstr('.',field_name);
      field_name(jx) = '_';
      field_val  = strtrim(s(ix+1:end));

      if (in_header)
         section = 'header';
         known_fields = known_fields_header;
      elseif (in_spline)
         section = 'spline';
         known_fields = known_fields_spline;
      else
         disp('ERROR: keyword outside a section');
         return
      end

      % check if the keyword is known
      if (isfield(known_fields, field_name))
         % interpret real numerical value
         if (strcmp(known_fields.(field_name), 'num') | strcmp(known_fields.(field_name), 'real'))
            p.(section).(field_name) = str2num(field_val);
         % interpret integer value
         elseif (strcmp(known_fields.(field_name), 'int'))
            p.(section).(field_name) = round(str2num(field_val));
         % store string value
         elseif (strcmp(known_fields.(field_name), 'str'))
            p.(section).(field_name) = field_val;
         % unknown type: store string value if it's non-empty
         elseif (~isempty(field_val))
            disp(sprintf('Unknown field-type %s: value="%s"',field_name,field_val));
            p.(section).(field_name) = field_val;
         end
      % unknown keyword: report & store string value
      else
         disp(sprintf('Unknown field "%s"',field_name));
         p.(section).(field_name) = field_val;
      end

   % no keyword found - interpret data line

   elseif (in_spline & in_point)

      [tmp, cnt] = sscanf(s, '%f %f %f %f %f');
      if (cnt > 0)
         nrows_data = nrows_data + 1;
         p.OrigData(nrows_data, 1:cnt) = tmp;
      end

   end % ~isempty(ix), keyword handling
end % while (~feof)
fclose(f);

% Apply modifications to original input data

opts   = p.spline;
points = p.OrigData;
np     = size(points,1);

% modification 1: remove points too close together

if (isfield(opts,'point_dist_min') & opts.point_dist_min>=1e-12)
   tmp = points;
   inew = 1;
   for ip = 2:np
      if (abs(tmp(inew,1)-points(ip,1))>opts.point_dist_min)
         inew = inew + 1;
         tmp(inew,:) = points(ip,:);
      end
   end
   points = tmp(1:inew,:);
   np     = inew;
end

% modifications 2 & 3: shift y and z

if (isfield(opts,'shift_y'))
   points(:,1) = points(:,1) + opts.shift_y;
end
if (isfield(opts,'shift_z'))
   points(:,2) = points(:,2) + opts.shift_z;
end

% modification 10: convert angle to radians

if (isfield(opts,'units_ang_f') & isfield(opts,'rotate'))
   opts.rotate = opts.rotate / opts.units_ang_f;
end

% modification 4: rotate profile about profile origin

if (isfield(opts,'rotate'))
   rot = [cos(opts.rotate), -sin(opts.rotate); sin(opts.rotate), cos(opts.rotate)];
   points(:,1:2) = points(:,1:2) * rot';
end

% modification 5: remove points outside [bound_y_min, bound_y_max]

if (isfield(opts,'bound_y_min') & isfield(opts,'bound_y_max') & opts.bound_y_min<opts.bound_y_max)
   ix = find(points(:,1)<opts.bound_y_min | points(:,1)>opts.bound_y_max);
   if (~isempty(ix))
      points(ix,:) = [];
      np = size(points,1);
   end
end

% modification 6: clipping to [bound_z_min, bound_z_max]

if (isfield(opts,'bound_z_min') & isfield(opts,'bound_z_max') & opts.bound_z_min<opts.bound_z_max)
   points(:,2) = max(opts.bound_z_min, min(opts.bound_z_max, points(:,2)));
end

% modification 7: mirror y

if (isfield(opts,'mirror_y') & opts.mirror_y==1)
   points(:,1) = - points(:,1);
end

% modification 8: mirror z

if (isfield(opts,'mirror_z') & opts.mirror_z==1)
   points(:,2) = - points(:,2);
end

% modification 9: invert data order

if (isfield(opts,'inversion') & opts.inversion==1)
   points = flipud(points);
end

% modification 10: conversion of lengths unit

if (isfield(opts,'units_len_f'))
   points(:,1:2) = 1000 * points(:,1:2) / opts.units_len_f;
end

% store end-result in structure

p.ProfileY  = points(:,1);
p.ProfileZ  = points(:,2);
p.SplineWgt = points(:,3);

if (make_plot)
   figure(make_plot); clf; hold on;
   plot(p.OrigData(:,1), p.OrigData(:,2), '-*');
   plot(p.ProfileY, p.ProfileZ, '-o');
   axis([-70 70 -5 30]);
   grid on;
   xlabel('y_{prf} [mm]'); ylabel('z_{prf} [mm]');
   set(gca,'ydir','reverse');
   legend('Original data','Modified data','location','NorthWest');
end

end % function read_simpack

