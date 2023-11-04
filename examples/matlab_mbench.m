
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example program using CONTACT library for Matlab. See User guide, Sections 7.5 & 5.7.
% Manchester contact benchmark - tangential contact.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1: Initialize the CONTACT library, register problems
%         "iwhe = 1" for left wheel and "iwhe = 2" for right wheel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~exist('cntc_initlibrary.m','file'))
   % set location of CONTACT installation folder
   % contactdir = 'C:\Program Files\Vtech CMCC\contact_v23.2';
   contactdir = '..';
   addpath([contactdir, '\matlab_intfc']);
   addpath([contactdir, '\matlab']);
end

if (exist('cntc_initlibrary')~=2)
   disp('ERROR: cant find CONTACT library on the Matlab search path');
   return
end
[CNTC, ifcver, ierror] = cntc_initlibrary;

if (ierror<0)
   disp(sprintf('An error occurred, ierror = %d, check output-file.',ierror));
   return
end

idebug = 1;             % 1: just a bit of information from the library
cntc_setglobalflags(CNTC.if_idebug, idebug);

for iwhe = 1 : 2 % "wheel number"
   imodul = 1; % w/r contact
   [ifcver, ierror] = cntc_initialize(iwhe, imodul);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2: Configure the main flags & control digits of the contact problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
for iwhe = 1 : 2 % "wheel number"


            % CONTACT unit convention: [mm], [mm/s], [N], acting on body (1)
            % note that w/r profiles are needed in [mm] regardless of this setting.
   clear flags values;
   flags(1) = CNTC.if_units ; values(1) = CNTC.un_cntc;

            % C1 digit (CONFIG): 0 for left side, 1 for right side
   flags(2) = CNTC.ic_config; values(2) = iwhe-1;
   flags(3) = CNTC.ic_tang  ; values(3) = 3;    % T=3: steady state rolling
   flags(4) = CNTC.ic_pvtime; values(4) = 2;    % P=2: no previous time
   flags(5) = CNTC.ic_discns; values(5) = 2;    % D=2: planar contact
   flags(6) = CNTC.if_wrtinp; values(6) = 0;    %   0: no .inp-file needed
   flags(7) = CNTC.ic_matfil; values(7) = 0;    % A=0: no .mat-file needed
   flags(8) = CNTC.ic_output; values(8) = 3;    % O=1: min. output to .out-file
   flags(9) = CNTC.ic_flow  ; values(9) = 2;    % W=2: little progress output

   cntc_setflags(iwhe, [], flags, values);

   % G=0: set maxgs, maxin, maxnr, maxout, eps

   gdigit  = 0;
   values  = [ 999, 100, 30, 1 ];
   rvalues = [ 1d-5 ];
   cntc_setsolverflags(iwhe, [], gdigit, values, rvalues);

   % material data using M = 0 (fully elastic model)

   gg    = 82000;  % [N/mm^2]
   poiss = 0.28;   % [-]

   mdigit = 0;
   rparam = [ poiss, poiss, gg, gg ];
   cntc_setmaterialparameters(iwhe, [], mdigit, rparam);

   % friction data using L = 0: Coulomb friction

   imeth = 0;      % L-digit 0: Coulomb friction
   fstat = 0.30;   % [-]
   cntc_setfrictionmethod(iwhe, [], imeth, [fstat,fstat]);

   % friction variation across rail profile:

   % imeth = 10;     % V-digit 1, L-digit 0: Coulomb friction, variable across rail
   % nvf = 2;        % 2 control points with linear interpolation, constant extrapolation
   % params = [ nvf, ...
   %            -20*pi/180, 0.20, 0.20, ...     % constant [0.2,0.2] for surface inclination <= -20deg
   %            -10*pi/180, 0.30, 0.30 ];       % linear [0.2,0.2] - [0.3,0.3] for inclin. -20 -- -10 deg
   % cntc_setfrictionmethod(iwhe, [], imeth, params);

   % total vertical force specified, setting the N-digit to 1

   fz = 10000;   % [N]

   cntc_setverticalforce(iwhe, fz);

   % grid discretization

   ipotcn = -1;     % -1 = w/r contact with fixed grid sizes & combination of patches
   dx     = 0.2;    % [mm]
   ds     = 0.2;    % [mm]
   a_sep  = pi/2;   % [rad]
   d_sep  = 8.0;    % [mm]
   d_comb = 4.0;    % [mm]
   params = [dx, ds, a_sep, d_sep, d_comb];

   cntc_setpotcontact(iwhe, [], ipotcn, params);

   dqrel = 1;
   cntc_setrollingstepsize(iwhe, [], [], dqrel);

   % track dimensions & deviations at current side (Z=3)

   ztrack = 3;
   params = [ 14, 0, 1435, 0, 0, 0, 0, 0, 0, 0, 0 ];
   cntc_settrackdimensions_new(iwhe, ztrack, params);

   % rail profile, same on both sides, scaling to [mm]

   iswheel = 0; % 0 = rail, 1 = wheel
   mirror_y = 0; iparam = [iswheel, 0, mirror_y];
   sclfac = 1; smooth = 0; rparam = [sclfac, smooth];

   cntc_setprofileinputfname(iwhe, 'MBench_UIC60_v3.prr', iparam, rparam);

   % wheelset dimensions

   ewheel = 3;
   params = [1360, -70, 460];
   cntc_setwheelsetdimensions(iwhe, ewheel, params);

   % wheel profile, provided to CONTACT as a table of values

   mirror_y = []; mirror_z = []; sclfac = []; idebug = 0;
   tmp = read_profile('MBench_S1002_v3.prw',[], mirror_y, mirror_z, sclfac, idebug);
   values = [tmp.ProfileY, tmp.ProfileZ];

   iswheel = 1; mirror_y = 0; iparam = [iswheel, 0, mirror_y];
   sclfac = 1; smooth = 0; rparam = [sclfac, smooth];

   cntc_setprofileinputvalues(iwhe, values, iparam, rparam);

   % positions for subsurface stress calculation

   iblk = 1;
   isubs = 1;
   nz = 11; zl = 1e-6; dz = 0.5;
   subs_addblock(iwhe, [], iblk, isubs, [], [], [nz, zl, dz]);

   % configure output of subsurface stress calculation

   clear flags values;
   flags(1) = CNTC.ic_sbsout; values(1) = 0;    % O_s=0: no output to out-file
   flags(2) = CNTC.ic_sbsfil; values(2) = 0;    % A_s=0: no subs-file

   cntc_setflags(iwhe, [], flags, values);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 3: perform loops for all the cases to be computed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set wheelset positions according to the example

y_ws    = [ 0 : 0.5 : 10 ];
yaw_ws  = [ 0 : 0.0012 : 0.024 ];
roll_ws = [  0.00000000, -0.00002304, -0.00005049, -0.00008103, -0.00011280, ...
            -0.00014570, -0.00018030, -0.00021680, -0.00025570, -0.00029770, ...
            -0.00035540, -0.00047770, -0.00062720, -0.00437600, -0.00639300, ...
            -0.00764200, -0.00860600, -0.00940800, -0.01010113, -0.01071386, ...
            -0.01126431 ];
vpitch  = [ -4.34811810, -4.34741340, -4.34657520, -4.34624400, -4.34591970, ...
            -4.34556270, -4.34515030, -4.34466880, -4.34409300, -4.34337150, ...
            -4.33536370, -4.33188640, -4.32937180, -4.27488340, -4.26356290, ...
            -4.25757470, -4.25348570, -4.25032450, -4.24775610, -4.24556650, ...
            -4.24363750 ];

% initialize result structures similar to parse_out1 results

for iwhe = 1 : 2
   results{iwhe} = struct('ws_pos', [], 'tot_forc', [], 'npatch', [], ...
                          'cp_pos', [], 'cp_creep', [], 'cp_force', []);
end

for iwhe = 1 : 2 % wheel number

   % loop over the lateral displacements y_ws

   for icase = 1: length(y_ws)

      ws_pos = [    0, y_ws(icase), 0, roll_ws(icase), yaw_ws(icase),   0      ];
      ws_vel = [ 2000,     0,       0,      0,            0,     vpitch(icase) ];
      ws_flx = [ 0, 0.1, 0, 0, 0.0001, 0, 0, 0, 0, 0, 0, 0.001 ];

      ewheel = 2;
      cntc_setwheelsetposition(iwhe, ewheel, ws_pos);
      cntc_setwheelsetvelocity(iwhe, ewheel, ws_vel);

      if (0==1 & iwhe==2)
         ewheel = 4;
         cntc_setwheelsetflexibility(iwhe, ewheel, ws_flx);
      end

      % compute the contact problem

      disp(sprintf('Starting case %2d for wheel %d...', icase, iwhe));

      ierror = cntc_calculate(iwhe);
      if (ierror~=0), return; end

      % compute subsurface stresses

      ierror = subs_calculate(iwhe);
      if (ierror~=0), return; end

      % get total forces on upper body (1) (global coordinates)

      values = cntc_getglobalforces(iwhe);
      results{iwhe}.tot_forc.fx_tr(icase) = values(1);
      results{iwhe}.tot_forc.fy_tr(icase) = values(2);
      results{iwhe}.tot_forc.fz_tr(icase) = values(3);
      results{iwhe}.tot_forc.fx_ws(icase) = values(7);
      results{iwhe}.tot_forc.fy_ws(icase) = values(8);
      results{iwhe}.tot_forc.fz_ws(icase) = values(9);

      % get number of contact patches

      results{iwhe}.npatch(icase) = cntc_getnumcontactpatches(iwhe);

      % get detailed results per contact patch

      for icp = 1 : results{iwhe}.npatch(icase)

         % get contact reference location

         values = cntc_getcontactlocation(iwhe, icp);

         results{iwhe}.cp_pos.xtr(icase,icp)    = values(1);
         results{iwhe}.cp_pos.ytr(icase,icp)    = values(2);
         results{iwhe}.cp_pos.ztr(icase,icp)    = values(3);
         results{iwhe}.cp_pos.delttr(icase,icp) = values(4);
         results{iwhe}.cp_pos.yr(icase,icp)     = values(6);
         results{iwhe}.cp_pos.zr(icase,icp)     = values(7);
         results{iwhe}.cp_pos.xw(icase,icp)     = values(10);
         results{iwhe}.cp_pos.yw(icase,icp)     = values(11);
         results{iwhe}.cp_pos.zw(icase,icp)     = values(12);

         % get reference velocity

         [veloc] = cntc_getreferencevelocity(iwhe, icp);
         results{iwhe}.cp_creep.veloc(icase,icp) = veloc;

         % get penetration and creepages 

         [pen] = cntc_getpenetration(iwhe, icp);
         [cksi, ceta, cphi] = cntc_getcreepages(iwhe, icp);

         results{iwhe}.cp_creep.pen(icase,icp)  = pen;
         results{iwhe}.cp_creep.cksi(icase,icp) = cksi;
         results{iwhe}.cp_creep.ceta(icase,icp) = ceta;
         results{iwhe}.cp_creep.cphi(icase,icp) = cphi;

         % get forces and moment in local coordinates

         [fn, tx, ty, mz] = cntc_getcontactforces(iwhe, icp);

         results{iwhe}.cp_force.fn(icase,icp)  = fn;
         results{iwhe}.cp_force.fx(icase,icp)  = tx;
         results{iwhe}.cp_force.fs(icase,icp)  = ty;
         results{iwhe}.cp_force.mz(icase,icp)  = mz;

         % get maximum von mises stress

         iblk = 1;
         table = subs_getresults(iwhe, icp, iblk, [1,2,3,8]);
         [vm_max, ii_max] = max(table(:,4));

         results{iwhe}.cp_force.sigvm(icase,icp) = vm_max;
         results{iwhe}.cp_force.vm_x(icase,icp)  = table(ii_max,1);
         results{iwhe}.cp_force.vm_y(icase,icp)  = table(ii_max,2);
         results{iwhe}.cp_force.vm_z(icase,icp)  = table(ii_max,3);

      end
         
      % get tractions at y_ws = 8mm (1st patch only)

      if (icase == 2 & iwhe == 1)
         [mx, my] = cntc_getnumelements;
         [pn, px, py] = cntc_gettractions(iwhe, 1);
      end

      % create pictures using plot3d

      if (icase == 11)
         clear sol;
         for icp = 1 : results{iwhe}.npatch(icase)
            sol(icp) = cntc_getcpresults(iwhe, icp);
         end

         opt = plot3d;
         if (iwhe==1)
            opt.typplot  = 'surf';
            opt.field    = 'ptabs+vec';
            opt.rw_surfc = 'prr';
         else
            opt.typplot  = 'rw_rear';
            opt.field    = 'pn';
            opt.rw_surfc = 'both';
         end

         figure(3+iwhe); clf;
         plot3d(sol, opt);
      end

   end % icase
end % iwhe

% Cleanup

for iwhe = 1 : 2
   cntc_finalize(iwhe);
end
cntc_closelibrary;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 4: plotting of the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

show_fig = [ 1 2 3 ];
make_tbl = 0;

if (~isempty(show_fig))
   if (exist('parse_out1')~=2)
      disp('ERROR: cant find Matlab script parse_out1 on search path');
      return
   end

   set(0, 'defaultaxesfontsize',14);
   set(0, 'defaulttextfontsize',14);
   set(0, 'defaultlinelinewidth',2);
   set(0, 'defaultaxeslinewidth',2);

   % read results of stand-alone program, runs for left & right wheels

   if (exist('mbench_a22_left.out'))
      lft = parse_out1('mbench_a22_left');
   elseif (exist('mbench_a22_left.ref_out'))
      lft = parse_out1('mbench_a22_left.ref_out');
   else
      disp('ERROR: no reference results for mbench_a22_left');
      return
   end
   if (exist('mbench_a22_right.out'))
      rgt = parse_out1('mbench_a22_right');
   elseif (exist('mbench_a22_right.ref_out'))
      rgt = parse_out1('mbench_a22_right.ref_out');
   else
      disp('ERROR: no reference results for mbench_a22_right');
      return
   end

   % read rail profile, convert to track coordinates

   mirror_y = []; mirror_z = []; sclfac = []; idebug = 0;
   prr = read_profile('MBench_UIC60_v3.prr',[], mirror_y, mirror_z, sclfac, idebug);
   prr.gaught = 14; prr.gaugwd = 1435;
   ix  = find(prr.ProfileZ<prr.gaught, 1, 'first');
   prr.ProfileY = prr.gaugwd/2 + prr.ProfileY - prr.ProfileY(ix);

   % 1: plot total forces Fx on left and right wheels

   if (any(show_fig==1))
      figure(1); clf; 
      subplot(2,1,1); hold on;
      plot(y_ws, results{1}.tot_forc.fx_tr);
      plot(lft.ws_pos.y, lft.tot_forc.fx_tr, '--');
      axis([0 10 -500 4000]); grid on;
      title('Left wheel');
      ylabel('Force [N]');
      legend('F_x, CONTACT library', 'F_x, stand-alone program', 'location','NorthWest');

      subplot(2,1,2); hold on;
      plot(y_ws, results{2}.tot_forc.fx_tr);
      plot(rgt.ws_pos.y, rgt.tot_forc.fx_tr, '--');
      axis([0 10 -4500 0]); grid on;
      xlabel('Wheelset displacement y_{ws} [mm]')
      title('Right wheel');
      ylabel('Force [N]');
      legend('F_x, CONTACT library', 'F_x, stand-alone program', 'location','SouthWest');
   end

   % 2: plot contact location on left and right wheels

   if (any(show_fig==2))
      figure(2); clf;
      subplot(2,1,1); hold on;
      plot(results{1}.cp_pos.ytr, y_ws);
      plot(lft.cp_pos.ytr(1,:), y_ws, '--');
      plot(-prr.ProfileY, prr.ProfileZ);
      axis([-760 -750 0 10]); grid on;
      set(gca,'ydir','reverse');
      title('Left wheel');
      ylabel('y_{ws} [mm]')
      legend('CONTACT library', 'stand-alone program', 'rail profile', 'location','SouthEast');
      text(-750.3, 2, 'User guide section 6.7', 'horizontalalignment','right');

      subplot(2,1,2); hold on;
      plot(results{2}.cp_pos.ytr(:,1), y_ws);
      plot(rgt.cp_pos.ytr(1,:), y_ws, '--');
      if (max(results{2}.npatch)>=2)
         set(gca,'colororderindex',1);
         ix = find(results{2}.npatch>=2);
         plot(results{2}.cp_pos.ytr(ix,2), y_ws(ix), '*');
      end
      if (max(rgt.npatch>=2));
         set(gca,'colororderindex',2);
         ix = find(rgt.npatch>=2);
         plot(rgt.cp_pos.ytr(2,ix), y_ws(ix), 'o', 'markersize', 12);
      end
      set(gca,'colororderindex',3);
      plot( prr.ProfileY, prr.ProfileZ);
      set(gca,'colororderindex',1);
      plot(rgt.cp_pos.ytr(1,:), rgt.cp_pos.ztr(1,:), '*', 'markersize',3);
      axis([ 715 750 0 10]); grid on;
      set(gca,'ydir','reverse');
      title('Right wheel');
      xlabel('Contact location y_{cp(tr)} [mm]');
      ylabel('y_{ws} [mm]')
      text(750, 5, 'z_{(tr)} [mm]', 'rotation',90, 'horizontalalignment','center', 'verticalalignment','top');
      legend('CONTACT library', 'stand-alone program', 'location','SouthEast');
   end

   % 3: plot maximum von Mises stresses

   if (any(show_fig==3))
      figure(3); clf; hold on;
      plot(y_ws, max(results{1}.cp_force.sigvm,[],2));
      plot(y_ws, max(results{2}.cp_force.sigvm,[],2));
      grid on;
      title('Subsurface: von Mises stress');
      legend('max(\sigma_{vm}), left wheel','max(\sigma_{vm}), right wheel', ...
             'location','NorthWest');
      xlabel('Wheelset displacement y_{ws} [mm]');
      ylabel('Stress [N/mm^2]');
   end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 5: table with overview of differences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (make_tbl)

   flds = strvcat('npatch', 'tot_forc.fx_tr', 'tot_forc.fy_tr', ...
                  'tot_forc.fz_tr', 'cp_pos.xtr', 'cp_pos.ytr', ...
                  'cp_pos.ztr', 'cp_creep.cksi', ...
                  'cp_creep.ceta', 'cp_creep.cphi', 'cp_force.fn' );

   for iwhe = 1 : 2
      if (iwhe==1)
         nam = 'Left';  ref = lft;
      else
         nam = 'Right'; ref = rgt;
      end
      disp(sprintf('%s side:', nam));

      for ifld = 1 : size(flds,1)
         fld = deblank(flds(ifld,:));
         eval(sprintf('f1 = results{iwhe}.%s;', fld));
         eval(sprintf('f2 = ref.%s;', fld));
         if (all(size(f1)==fliplr(size(f2)))), f2 = f2'; end
         if (any(size(f1)~=size(f2)))
            disp(['Different sizes for field ',fld,':'])
            disp([size(f1), size(f2)])
         else
            dif = f1 - f2;
            rel_dif = dif ./ max(abs(f1), abs(f2));
            [~,i] = max(abs(dif));
            if (length(i)<=1)
               str = [];
               if (abs(rel_dif(i))>0.0001)
                  str = sprintf(' (%8.3f %%)', 100*rel_dif(i));
               end
               disp(sprintf('    max difference %15s = %7.3f at y=%5.2f %s', ...
                                              fld, dif(i), y_ws(i), str))
            else
               disp(sprintf('    max difference %15s = %7.3f at y=%5.2f (cp1), %7.3f at y=%5.2f (cp2)', ...
                              fld, dif(i(1)), y_ws(i(1)), dif(i(2)), y_ws(i(2)) ))
            end
         end
      end
   end
end

% $Revision: 2447 $, $Date: 2023-11-04 15:03:17 +0100 (Sat, 04 Nov 2023) $
