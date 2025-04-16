function [ p_out ] = modify_profile(p_in, is_wheel, mirror_y, mirror_z, scale_yz, y_tape, ...
                                    reverse_order, point_dist_min, delete_foldback, make_plot, idebug)

% function [ p_out ] = modify_profile(p_in, is_wheel, mirror_y, mirror_z, scale_yz, y_tape, ...
%                                     reverse_order, point_dist_min, delete_foldback, make_plot, idebug)
%
% apply a number of common filters to the input profile p_in
%   - mirror_y         - change y(i)  --> -y(i)    <=0 means no, >=1 yes
%   - mirror_z         - change z(i)  --> -z(i)    <=0 means no, >=1 yes
%   - scale_yz         - scale y(i), z(i) --> scale*y(i), scale*z(i)
%   - y_tape           - align wheel vertically with z_tape = 0 at y = y_tape
%   - reverse_order    - change [1:n] --> [n:-1:1]
%   - point_dist_min   - delete points that lie too close together
%   - delete_foldback  - identify & remove points where y(i+1)<y(i) (rail) or y(i+1)>y(i) (wheel)

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

if (nargin<2 | isempty(is_wheel))
   is_wheel = 0;
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
if (nargin<6 | isempty(y_tape))
   y_tape   = [];
end
if (nargin<7 | isempty(reverse_order))
   reverse_order = 0;
end
if (nargin<8 | isempty(point_dist_min))
   point_dist_min = -1;
end
if (nargin<9 | isempty(delete_foldback))
   delete_foldback = 0;
end
if (nargin<10 | isempty(make_plot))
   make_plot = 0;
end
if (nargin<11 | isempty(idebug))
   idebug = 1;
end

% Initialize output structure

p_out          = p_in;
if (~isfield(p_out, 'ProfileAngle') | isempty(p_out.ProfileAngle))
   p_out.ProfileAngle = zeros(size(p_out.ProfileY));
end
if (~isfield(p_out, 'ProfileCurvature') | isempty(p_out.ProfileCurvature))
   p_out.ProfileCurvature = zeros(size(p_out.ProfileY));
end
if (~isfield(p_out, 'ProfileData'))
   p_out.ProfileData = [ p_out.ProfileY, p_out.ProfileZ, ...
                         p_out.ProfileAngle, p_out.ProfileCurvature];
end
p_out.OrigData = p_out.ProfileData;

% Apply modifications to original input data

% modification 1: mirror horizontally w.r.t. y=0, changing y into -y

if (mirror_y>=1)
   if (idebug>=1)
      disp('Mirroring y-values from left-side to right-side profile...');
   end
   p_out.ProfileY     = -p_out.ProfileY;
   p_out.ProfileAngle = -p_out.ProfileAngle;
end

% modification 2: mirror vertically w.r.t. z=0, changing z into -z

if (mirror_z>0)
   if (idebug>=1)
      disp('Mirroring z-values to get z positive downwards...');
   end
   p_out.ProfileZ     = -p_out.ProfileZ;
   p_out.ProfileAngle = -p_out.ProfileAngle;
end

% modification 3: scale y,z-coordinates

if (scale_yz~=1)
   if (idebug>=1)
      disp(sprintf('Scaling y,z-values with factor %3.1f to convert to mm...',scale_yz));
   end
   p_out.ProfileY     = scale_yz * p_out.ProfileY;
   p_out.ProfileZ     = scale_yz * p_out.ProfileZ;
   p_out.ProfileCurvature = p_out.ProfileCurvature / scale_yz;
end

% modification 4: align wheel z at tape circle line

if (is_wheel & ~isempty(y_tape))
   ix = find(abs(p_out.ProfileY-y_tape)<5);
   z_tape = interp1(p_out.ProfileY(ix), p_out.ProfileZ(ix), y_tape);
   if (idebug>=1)
      disp(sprintf('Shifting wheel z-values by dz = %4.2f to have z=0 at y_tape = %3.1f...',-z_tape, y_tape));
   end
   p_out.ProfileZ     = p_out.ProfileZ - z_tape;
end

% modification 5: reverse the order of the data points

if (reverse_order)
   if (idebug>=1)
      disp('Reversing the order of the points, get inside at right hand...');
   end
   p_out.ProfileY         = flipud(p_out.ProfileY);
   p_out.ProfileZ         = flipud(p_out.ProfileZ);
   p_out.ProfileAngle     = flipud(p_out.ProfileAngle);
   p_out.ProfileCurvature = flipud(p_out.ProfileCurvature);
end

% modification 6: remove points that lie too close together

if (point_dist_min>=1e-6)
   np   = length(p_out.ProfileY);
   dy   = diff(p_out.ProfileY);
   dz   = diff(p_out.ProfileZ);
   dist = sqrt(dy.^2 + dz.^2);
   cum  = [0; cumsum(dist)];

   keep = ones(np,1);
   ik   = 1; % highest number that is kept so far
   for ip = 2 : np
      if (cum(ip)-cum(ik) < point_dist_min)
         keep(ip) = 0;
      else
         ik = ip;
      end
   end
   ix0 = find(keep==0); ix1 = find(keep==1);

   if (idebug>=1 & ~isempty(ix0))
      for i = ix0'
         disp(sprintf('Deleting (%7.2f,%7.2f) (ip=%4d), too close to adjacent', ...
                                                             p_out.ProfileY(i), p_out.ProfileZ(i), i));
      end
   end
   if ((idebug>=1 & ~isempty(ix0)) | idebug>=2)
      disp(sprintf('Point_dist_min=%6.0e: deleted %d points from profile, %d remaining', point_dist_min, ...
           length(ix0), length(ix1)));
   end

   p_out.ProfileY         = p_out.ProfileY(ix1);
   p_out.ProfileZ         = p_out.ProfileZ(ix1);
   p_out.ProfileAngle     = p_out.ProfileAngle(ix1);
   p_out.ProfileCurvature = p_out.ProfileCurvature(ix1);
end

% modification 7: remove points where the profile folds back
%                 (i.e. where the function y --> z(y) is multi-valued)

if (~delete_foldback)
   p_out.Mask             = ones(size(p_out.ProfileY));
else
   if (idebug>=2)
      disp('Checking for points where the profile folds back...');
   end

   % copy data to work variables
   points0 = [p_out.ProfileY, p_out.ProfileZ, p_out.ProfileAngle, p_out.ProfileCurvature];
   np0 = size(points0,1);

   % find "midpoint" of profile: tread of rail, bottom of flange on wheel
   % z is assumed positive downwards: tread == minimum z, flange == maximum z
   if (is_wheel)
      [~,i_mid0] = max(points0(:,2));
   else
      [~,i_mid0] = min(points0(:,2));
   end

   % set threshold "dy" below which y-values are "the same"
   dy_thresh = 1e-6 * (max(points0(:,1))-min(points0(:,1)));

   % delete points where the profile is "vertical", i.e. |y(i+1)-y(i)|<thresh
   % where (i<i_mid0) keep i, where (i>i_mid0) keep i+1.
   
   dy   = diff(points0(:,1));
   idel = find(abs(dy)<dy_thresh);
   if (isempty(idel))
      points1 = points0; np1 = np0; ikeep1 = [1:np0]; i_mid1=i_mid0;
   else
      j = find(idel<i_mid0); idel(j) = idel(j) + 1;
      keep1 = ones(np0,1); keep1(idel) = 0; ikeep1 = find(keep1);
      points1 = points0(ikeep1,:); np1 = size(points1,1);
      i_mid1 = nnz(keep1(1:i_mid0));
      if (idebug>=1)
         disp(sprintf('Identified %d points with vertical slope', length(idel)));
      end
   end

   % delete points where the profile folds back, i.e. y(i+1) < yrmx(i) (rail)
   %                                              or  y(i+1) > yrmx(i) (wheel)

   keep2 = zeros(np1,1); keep2(i_mid1) = 1;
   y = points1(:,1);
   if (is_wheel), sgn = -1; else, sgn = 1; end

   % rail: check right of mid-point (wheel: left of mid-point)
   y_runmax = sgn*y(i_mid1);
   for i = i_mid1+1 : np1
      if (sgn*y(i) >= y_runmax+dy_thresh)
         keep2(i) = 1;
      end
      y_runmax = max(y_runmax, sgn*y(i));
   end

   % rail: check left of mid-point (wheel: right) 
   y_runmin = sgn*y(i_mid1);
   for i = i_mid1-1 : -1 : 1
      if (sgn*y(i) <= y_runmin-dy_thresh)
         keep2(i) = 1;
      end
      y_runmin = min(y_runmin, sgn*y(i));
   end
   
   idel = find(keep2==0);
   if (isempty(idel))
      points2 = points1; np2 = np1; ikeep2 = [1:np1];
   else
      ikeep2 = find(keep2);
      points2 = points1(ikeep2,:); np2 = size(points2,1);
      if (idebug>=1)
         disp(sprintf('Identified %d points where profile folds back', length(idel)));
      end
   end

   % get mask that combines keep1 and keep2

   msk = zeros(np0,1);
   msk( ikeep1(ikeep2) ) = 1;
   
   if (delete_foldback<=1)
      p_out.Mask             = msk;
      p_out.ProfileY         = points0(:,1);
      p_out.ProfileZ         = points0(:,2);
      p_out.ProfileAngle     = points0(:,3);
      p_out.ProfileCurvature = points0(:,4);
   else
      p_out.Mask             = ones(np2,1);
      p_out.ProfileY         = points2(:,1);
      p_out.ProfileZ         = points2(:,2);
      p_out.ProfileAngle     = points2(:,3);
      p_out.ProfileCurvature = points2(:,4);
   end
end

% combine the columns into array ProfileData

p_out.ProfileData = [p_out.ProfileY, p_out.ProfileZ, p_out.ProfileAngle, p_out.ProfileCurvature];
p_out.XYPoints = size(p_out.ProfileData,1);

% plot the modified profile if requested

if (make_plot)
   figure(make_plot); clf; hold on;
   plot(p_out.OrigData(:,1), p_out.OrigData(:,2), '-*');
   plot(p_out.ProfileY, p_out.ProfileZ, '-o');
   % axis([-70 70 -5 30]);
   grid on;
   xlabel('y_{prf} [mm]'); ylabel('z_{prf} [mm]');
   set(gca,'ydir','reverse');
   % legend('Original data','Modified data','location','NorthWest');
end

end % function modify_profile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

