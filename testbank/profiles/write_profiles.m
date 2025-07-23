
function write_profiles()

write_circ_wr    = 1;
write_concave_wh = 0;
write_flat_wr    = 0;
write_groove     = 0;

% circular rail and wheel profiles

if (write_circ_wr==1)

   for R_ry = [460]
      dy = 0.1;
      if (R_ry>101)
         ny = round(100 / dy) - 1;
         y  = -50 + [1:ny]*dy;
      else
         ny = round(2*R_ry / dy) - 1;
         y  = -R_ry + [1:ny]*dy;
      end
      z  =  R_ry - sqrt(R_ry^2 - y.^2);

      % shave off minimum z-value
      % iy_mid = (ny+1) / 2;
      % z(iy_mid) = z(iy_mid-1) - 1.2e-6;
      % z  = z - min(z);

      fname   = sprintf('circ_r%d.prr', round(R_ry));
      comment = sprintf('circular rail, R_ry=%6.2f', R_ry);
      write_simpack(fname, comment, y,  z);

      fname   = sprintf('circ_r%d.prw', round(R_ry));
      comment = sprintf('circular wheel, R_wy=%6.2f', R_ry);
      write_simpack(fname, comment, y, -z);
   end

   if (1==1)
      figure(1); clf;
      plot(y, z)
      axis equal
      set(gca,'ydir','reverse');
      grid on
   end
end

% concave wheel profiles

if (write_concave_wh==1)

   for R_wy = [-10.1, -10.2, -10.5, -11.0, -12.0, -13.0 ]
      dy = 0.1;
      ny = round(-2*R_wy / dy) - 1;
      y  =  R_wy + [1:ny]*dy;
      z  = -R_wy - sqrt(R_wy^2 - y.^2);

      fname   = sprintf('concave_r%d.prw', round(-10*R_wy));
      comment = sprintf('concave wheel, R_wy=%7.2f', R_wy);
      write_simpack(fname, comment, y, z);
   end

   if (1==1)
      figure(1); clf;
      plot(y, z)
      axis equal
      set(gca,'ydir','reverse');
      grid on
   end
end

% flat (conical) wheel profiles at different cone angles

if (write_flat_wr==1)

   for ang = [0, 0.002, 0.01, 0.02, 0.1, 0.2, pi/6]
      y = [ -50 : 1 : 50 ];
      z = tan(ang) * y;

      fname   = sprintf('flat_ang%03d.prw', round(1000*ang));
      comment = sprintf('flat (conical) wheel, angle = %6.4f rad', ang);
      write_simpack(fname, comment, y, z);

      fname   = sprintf('flat_ang%03d.prr', round(1000*ang));
      comment = sprintf('flat (inclined) rail, angle = %6.4f rad', ang);
      write_simpack(fname, comment, y, z);
   end
end

% rail: groove for deep-groove rolling

if (write_groove==1)

   for R_ry = [ 10.02, 10.05, 10.1, 10.2, 10.25, 10.5, 11.0, 12.0, 15.0, 20.0, 30.0, 60.0 ]
      dy = 0.1;
      ny = round(4*R_ry / dy) - 1;
      y  = -2*R_ry + [1:ny]*dy;
      z  =   -R_ry + sqrt(max(0,R_ry^2 - y.^2));

      fname   = sprintf('groove_r%d.prr', round(100*R_ry));
      comment = sprintf('grooved surface, R_ry=%7.2f', -R_ry);
      write_simpack(fname, comment, y,  z);
   end

   if (1==1)
      figure(1); clf;
      plot(y, z)
      axis equal
      set(gca,'ydir','reverse');
      grid on
   end
end


end % write_profiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_simpack(fname, comment, y, z)

if (isempty(strfind(fname, 'prw')))
   is_wheel = 0; % rail
else
   is_wheel = 1; % wheel
end

% make y-data decreasing for wheel

if (is_wheel & any(diff(y))>0)
   y = fliplr(y);  z = fliplr(z);
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
   fprintf(f,'    %10.6f  %10.6f\n',y(iy), z(iy));
end
fprintf(f,'  point.end\n');
fprintf(f,'spline.end\n');
fclose(f);

end % write_simpack
