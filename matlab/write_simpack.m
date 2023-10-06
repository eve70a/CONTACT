
function write_simpack(fname, comment, y, z)

% function write_simpack(fname, comment, y, z)
%
% Write wheel or rail profile in Simpack prw/prr format

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

if (isempty(strfind(fname, 'prw')))
   is_wheel = 0; % rail
else
   is_wheel = 1; % wheel
end

% make column vectors

if (size(y,2)>1), y=y'; end
if (size(z,2)>1), z=z'; end
if (max(size(y,2), size(z,2))>1),
   disp('ERROR: profile y, z should be 1D vectors')
   disp(size(y))
   disp(size(z))
   return
end
if (length(y)~=length(z))
   disp(sprintf('ERROR: profile y (%d), z (%d) should be the same length', ...
        length(y), length(z)))
   return
end

% make y-data decreasing for wheel

if (is_wheel & y(end)-y(1)>0)
   y = flipud(y);  z = flipud(z);
end

f = fopen(fname, 'w');
fprintf(f,'header.begin\n');
fprintf(f,'  version  = 1\n');
fprintf(f,'  type     = %d\n', is_wheel);
fprintf(f,'header.end\n');
fprintf(f,'spline.begin\n');
fprintf(f,'  approx.smooth = 0.000000\n');
fprintf(f,'  comment       = ''%s''\n', comment);
fprintf(f,'  units.len.f   = 1000.\n');
fprintf(f,'  units.ang.f   =    1.\n');
fprintf(f,'  point.begin\n');
for iy = 1 : length(y)
   fprintf(f,'   %11.6f %11.6f\n',y(iy), z(iy));
end
fprintf(f,'  point.end\n');
fprintf(f,'spline.end\n');
fclose(f);

end % write_simpack
