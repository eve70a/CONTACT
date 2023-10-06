function [ p ] = read_vampire(fname, right_side, idebug, make_plot)

% function [ p ] = read_vampire(fname, [right_side], [idebug], [make_plot])
%
% Lower-level routine for reading wheel/rail profiles in Vampire whe/rai format.
% Does not do automatic corrections like mirroring or reversing order.
%
%   fname      - filename of input file
%   right_side - select right side profile (1, default) or left side (0)
%   make_plot  - figure number to use for plotting the profile

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

if (nargin<2 | isempty(right_side))
   right_side = 1;
end
if (nargin<3 | isempty(idebug))
   idebug = 1;
end
if (nargin<4 | isempty(make_plot))
   make_plot = 0;
end

% Initialize output structure

p                  = struct;
p.ProfileData      = [];
nrows_data         = 0;
p.ProfileY         = [];
p.ProfileZ         = [];
inch               = 25.4; % [mm/inch]

if (~exist(fname,'file'))
   disp(['ERROR: cannot find file ',fname]);
   return
end

f = fopen(fname,'r');
iline = 0;
in_profile = 0;

% read first line, should contain 'WHEELPROFILE' or 'RAILPROFILE'

[ s, iline ] = get_next_line(f, iline, idebug);
if (strcmp(s, 'WHEELPROFILE'))
   p.is_wheel = 1;
elseif (strcmp(s, 'RAILPROFILE'))
   p.is_wheel = 0;
else
   disp(sprintf('ERROR: file ''%s'' is not in Vampire format.', fname));
   return;
end

% second line: 80-character title

[ s, iline ] = get_next_line(f, iline, idebug);
p.Title = s;

% third line: tread datum (wheel) or gauge point (rail)

if (p.is_wheel)

   % third line, wheel: tread datum 

   [ s, iline ] = get_next_line(f, iline, idebug);
   ix = length('TREADDATUM') + 1;
   p.TreadDatum = strtrim(s(ix:end));

   ix = strfind('CUSTOM', p.TreadDatum);
   if (~isempty(ix))
      p.FlangeBackPos = sscanf(p.TreadDatum(ix+6:end), ' %f');
   elseif (strfind('EUROPEAN', p.TreadDatum))
      p.FlangeBackPos = 70;
   elseif (strfind('US1', p.TreadDatum))
      p.FlangeBackPos = (2 + 27/32) * inch;
   elseif (strfind('US2', p.TreadDatum))
      p.FlangeBackPos = (3 + 1/16) * inch;
   else
      disp(['Error: unknown tread datum ',p.TreadDatum]);
   end
   p.FlangeBackPos = -p.FlangeBackPos;     % CONTACT: from datum to flange back

else

   % third line, rail: gauge point

   [ s, iline ] = get_next_line(f, iline, idebug);
   ix = length('GAUGEPOINT') + 1;
   p.GaugePoint = strtrim(s(ix:end));

   ix = strfind('CUSTOM', p.GaugePoint);
   if (~isempty(ix))
      p.GaugeHeight = sscanf(p.GaugePoint(ix+6:end), ' %f');
   elseif (strfind('EUROPEAN', p.GaugePoint))
      p.GaugeHeight = 14;
   elseif (strfind('US', p.GaugePoint))
      p.GaugeHeight = 5/8 * inch;
   else
      disp(['Error: unknown gauge point ',p.GaugePoint]);
   end
end

% fourth line: flange back distance (wheel) or gauge (rail)

[ s, iline ] = get_next_line(f, iline, idebug);

if (p.is_wheel)
   ix = length('FLANGEBACK') + 1;
   p.FlangeBackDist = sscanf(s(ix:end), '%f');
else
   ix = length('GAUGE') + 1;
   p.GaugeWidth = sscanf(s(ix:end), '%f');
end

while(~feof(f))

   % get next line from input

   [ s, iline ] = get_next_line(f, iline, idebug);

   % remove comments, everything after ", !, or %
   ix = strfind(s, '"'); if (~isempty(ix)), s=s(1:ix-1); end
   ix = strfind(s, '!'); if (~isempty(ix)), s=s(1:ix-1); end
   ix = strfind(s, '%'); if (~isempty(ix)), s=s(1:ix-1); end

   % get (yl,zl), (yr,zr) positions on left and right profiles

   ix = strfind(s, ','); if (~isempty(ix)), s(ix)=' '; end
   [tmp, cnt] = sscanf(s, '%f %f %f %f %f %f %f %f');

   if (cnt>0)
      nrows_data = nrows_data + 1;
      p.ProfileData(nrows_data, 1:cnt) = tmp;
   end

end % while (~feof)
fclose(f);

% select left or right side profile, convert to right side profile

if (right_side)
   p.ProfileY =  p.ProfileData(:,3);
   p.ProfileZ =  p.ProfileData(:,4);
else
   p.ProfileY = -p.ProfileData(:,1);
   p.ProfileZ =  p.ProfileData(:,2);
end

% renumber rails from track center to field side

if (~p.is_wheel & (p.ProfileY(end)<p.ProfileY(1)))
   p.ProfileY = flipud(p.ProfileY);
   p.ProfileZ = flipud(p.ProfileZ);
end

% shift the datum to zero

if (p.is_wheel)
   fb_dist = p.FlangeBackDist;
   fb_pos  = p.FlangeBackPos;
   p.ProfileY = p.ProfileY - fb_dist/2 + fb_pos;
else
   p.ProfileZ = -p.ProfileZ;    % Vampire rails: z positive upwards

   zmin = min(p.ProfileZ);      % interpolate to get the gauge point
   ix = find(p.ProfileZ-zmin < p.GaugeHeight, 1, 'first');
   if (~isempty(ix) & ix>1)
      ygauge = interp1( p.ProfileZ(ix+[-1,0])-zmin, p.ProfileY(ix+[-1,0]), p.GaugeHeight );
   else
      ygauge = p.ProfileY(1);
   end
   % shift rail to put the gauge point at y_r = 0
   p.ProfileY = p.ProfileY - ygauge;
   p.ProfileZ = p.ProfileZ - zmin;
   p.rail_y0  = ygauge;
   p.rail_z0  = zmin;
end

if (make_plot)
   figure(make_plot); clf; hold on;
   plot(p.ProfileY, p.ProfileZ, '-o');
   grid on;
   set(gca,'ydir','reverse');
   if (p.is_wheel)
      xlabel('y_w [mm]'); ylabel('z_w [mm]');
   else
      xlabel('y_r [mm]'); ylabel('z_r [mm]');
   end
   title(p.Title)
end

end % function read_vampire

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ s, iline ] = get_next_line(f, iline, idebug)

   s = deblank(fgets(f)); 
   iline = iline + 1;
   if (idebug>=5 & iline<10)
      disp(sprintf('%4d: %s', iline, s));
   end

end % function get_next_line
