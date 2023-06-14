
% Determine contact points for MBench example including the effect of roll
% Create pictures for User Guide.

print_fig = 0;
idebug = 0;

addpath('../../../meldingen/afgesloten/m10023-example-s1002/2-contact-posit');

prw = read_simpack('../../examples/MBench_S1002_v3.prw');
prr = read_simpack('../../examples/MBench_UIC60_v3.prr');

% Note: prw-file has flange-back at -70mm, points in descending order.

ywhe = flipud(prw.ProfileY); zwhe = flipud(prw.ProfileZ);
yrai =        prr.ProfileY;  zrai =        prr.ProfileZ;
clear prr prw;

% Note: the rail local coordinates have origin at the highest point;
% find lateral position on inner side of rail (y<0) with height 14mm

ix = find(yrai<0 & zrai>14, 1, 'last');
dy = yrai(ix+1) - yrai(ix);
dz = zrai(ix+1) - zrai(ix);
ytarg = yrai(ix) + (14 - zrai(ix)) / dz * dy; % 43.031 mm
clear ix dy dz ytarg

% determine shifts for wheel and rail to (global) track-based coordinates

yshft_whe = 1360/2 - -70;   % put flange-back (-70 mm) at 1360/2 mm
yshft_rai = 1435/2 - -43.0; % put point of rail with height 14mm at 1435/2 mm

% convert the wheel and rail profiles to the (global) track-based system

ywhe = ywhe + yshft_whe;
yrai = yrai + yshft_rai;

% set the desired lateral position

for lat_pos = [ 6.0 ],

   iter      = 0;
   roll_angl = 0;      % degrees; negative phi lifts the right wheel.
   tol_dz    = 0.001;  % iteration tolerance
   dzl       = 0;
   dzr       = 10*tol_dz;

   % Perform iteration to determine the appropriate roll angle

   while (abs(dzr-dzl) > tol_dz)
      iter = iter + 1;
    
      % construct new wheel-profiles for left and right wheels, rotated clockwise
      % about the "wheel-set origin", which is at (0, -460) in the global
      % track-based coordinate system.
      % a positive phi moves the right wheel downward and lifts the left wheel.
      % [ ynew ]    [   0   ]   [ cos  -sin ]   [ y -   0   ]
      % [ znew ]  = [ -460  ] + [ sin   cos ] * [ z - -460  ]
    
      cr = cos(roll_angl*pi/180);
      sr = sin(roll_angl*pi/180);
      ywhe_rot_l =   0  + cr*(ywhe - 0) + sr*(zwhe + 460);
      zwhe_rot_l = -460 - sr*(ywhe - 0) + cr*(zwhe + 460);
      ywhe_rot_r =   0  + cr*(ywhe - 0) - sr*(zwhe + 460);
      zwhe_rot_r = -460 + sr*(ywhe - 0) + cr*(zwhe + 460);
      
      % shift the wheel profiles laterally
    
      ywhe_rot_l = ywhe_rot_l - lat_pos;
      ywhe_rot_r = ywhe_rot_r + lat_pos;
    
      % determine the contact positions for left and right wheels at current
      % shift lat_pos and roll_angle.
    
      [yr_l, zr_l, angle_l, dzl, pen_l] = ...
                             contact_point(ywhe_rot_l, zwhe_rot_l, yrai, zrai);
      [yr_r, zr_r, angle_r, dzr, pen_r] = ...
                             contact_point(ywhe_rot_r, zwhe_rot_r, yrai, zrai);
    
      % Note: pen==separation, negative means interpenetration, minimum occurs 
      % at the contact point.
    
      % Determine the location of the contact points on left/right wheels
      % in the wheel-based coordinate systems: undo lat/vert shift, roll angle
    
      yw_l_intm = yr_l + lat_pos;
      yw_r_intm = yr_r - lat_pos;
      zw_l_intm = zr_l - dzl;
      zw_r_intm = zr_r - dzr;
      yw_l =   0  + cr*(yw_l_intm - 0) - sr*(zw_l_intm + 460) - yshft_whe;
      zw_l = -460 + sr*(yw_l_intm - 0) + cr*(zw_l_intm + 460);
      yw_r =   0  + cr*(yw_r_intm - 0) + sr*(zw_r_intm + 460) - yshft_whe;
      zw_r = -460 - sr*(yw_r_intm - 0) + cr*(zw_r_intm + 460);
      clear sr cr yw_l_intm zw_l_intm yw_r_intm zw_r_intm;
      
      if (idebug>=1)
         disp(sprintf('Iteration %2d: roll_angle phi=%7.4f deg : error delta z = %6.4f mm',...
                   iter, roll_angl, dzr-dzl));
      end
    
      if (abs(dzr-dzl)>tol_dz)
      % set the desired roll angle for the next iteration
         roll_angl= roll_angl + (dzr-dzl)/1520*180/pi;
      end
    
   end % while (dzl-dzr>eps)

   if (1==1)
       % plot left wheel / left rail
      figure(3); clf; hold on
      set(gca,'ydir','reverse');
      axis equal; grid on;

      c(1) = plot(-ywhe_rot_l, zwhe_rot_l+dzl);
      % set(gca,'colororderindex',1);
      % plot(-yr_l, zr_l, '*');
      c(2) = plot(-yrai, zrai);

      siz_o = 0.12;
      fac_txt = [1 1; 1 -1; 1 1.4];
      plot_axes([-yr_l,zr_l], 12, -angle_l, -1, 'nsx', 3, [], siz_o, [], [], ...
             [], fac_txt);

      cs = cos(-angle_l*pi/180);
      sn = sin(-angle_l*pi/180);
      set(gca,'colororderindex',3);
      plot(-yr_l+[-20,25]*cs, zr_l+[-20,25]*sn, 'linewidth',1);

      xlabel('y_{tr} [mm]')
      ylabel('z_{tr} [mm]');
      legend(c(1:2), 'left wheel', 'rail profile', 'location','northeast');

      l=text(-760,18,sprintf('y_{ref}= %3.1f mm', -yr_l), ...
                                        'horizontalalignment','center');
      l=text(-760,25,sprintf('\\delta_{ref}= %5.1f deg', -angle_l), ...
                                        'horizontalalignment','center');
      axis([-800 -670 -10 30]);

      if (print_fig)
         nam=sprintf('init_cp_left_%3.1fmm',lat_pos); nam=strrep(nam,'.','_');
         eval(sprintf('print -djpeg95 %s.jpg',nam));
      end
   end
   
   if (1==1)
        % plot right wheel / right rail
      figure(4); clf; hold on
      set(gca,'ydir','reverse');
      axis equal; grid on;

      c(1) = plot(ywhe_rot_r, zwhe_rot_r+dzr);
      % set(gca,'colororderindex',1);
      % plot(yr_r,zr_r,'o');
      c(2) = plot(yrai, zrai);

      siz_o = 0.12;
      fac_txt = [1 1; 1 -1; 1 1.4];
      plot_axes([yr_r,zr_r], 12, angle_r, -1, 'nsx', 3, [], siz_o, [], [], ...
             [], fac_txt);

      cs = cos(angle_r*pi/180);
      sn = sin(angle_r*pi/180);
      set(gca,'colororderindex',3);
      plot(yr_r+[-20,25]*cs, zr_r+[-20,25]*sn, 'linewidth',1);

      xlabel('y_{tr} [mm]')
      ylabel('z_{tr} [mm]');
      legend(c(1:2), 'right wheel', 'rail profile', 'location','northwest');
      l=text(760, 18, sprintf('y_{ref}= %5.1f mm', yr_r), ...
                                        'horizontalalignment','center');
      l=text(760, 25, sprintf('\\delta_{ref}= %5.1f deg', angle_r), ...
                                        'horizontalalignment','center');
      % l=text(725,25,sprintf('right wheel \\delta z=%5.3f mm',dzr));
      axis([670 800 -10 30]);
      clear l;
   
      if (print_fig)
         nam=sprintf('init_cp_right%3.1fmm',lat_pos); nam=strrep(nam,'.','_');
         eval(sprintf('print -djpeg75 %s.jpg',nam));
      end
   end
   
   clear ywhe_rot_l zwhe_rot_l ywhe_rot_r zwhe_rot_r;
    
   disp(sprintf('Lat.pos %5.2f mm: left contact at %6.2f mm, %5.1f deg, right at %6.2f mm, %5.1f deg', ...
        lat_pos, -yr_l+yshft_rai, -angle_l, yr_r-yshft_rai,  angle_r));

   clear roll_angl dzl yr_l zr_l yw_l zw_l angle_l pen_l;
   clear           dzr yr_r zr_r yw_r zw_r angle_r pen_r;

end

clear idebug iter lat_pos tol_dz
