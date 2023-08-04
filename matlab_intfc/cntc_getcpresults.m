
%------------------------------------------------------------------------------------------------------------
% function [ sol ] = cntc_getcpresults(ire, icp)
%
% get results of a single contact patch in the form used by loadcase and plot3d
%------------------------------------------------------------------------------------------------------------

% Copyright 2008-2023 by Vtech CMCC.
%
% Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

% category 6: m=*, cp     - require icp>0, default 1

function [ sol ] = cntc_getcpresults(ire, icp)

   if (nargin<1 | isempty(ire))
      ire = 1;
   end
   if (nargin<2 | isempty(icp))
      icp = 1;
   end
   if (icp<=0)
      disp(sprintf('ERROR in cntc_getcpresults: not available for icp=%d',icp));
      return
   end

   sol         = struct();

   CNTC = cntc_getmagicnumbers();

   % retrieve control digits needed for plotting

   iparam(1) = CNTC.ic_config;
   iparam(2) = CNTC.ic_tang;
   iparam(3) = CNTC.ic_frclaw;
   iparam(4) = CNTC.ic_discns;
   iparam(5) = CNTC.ic_mater;
   iparam(6) = CNTC.ic_heat;
   values = cntc_getflags(ire, icp, iparam);

   sol.config         = values(1);
   sol.kincns         = struct();
   sol.kincns.t_digit = values(2);
   % sol.fric           = struct();
   % sol.fric.frclaw    = values(3);
   sol.d_digit        = values(4);
   sol.mater          = struct();
   sol.mater.m_digit  = values(5);
   sol.h_digit        = values(6);

   % get material / kinematic parameters needed

   values = cntc_getparameters(ire, icp);
   sol.kincns.veloc = values.veloc;
   sol.kincns.chi   = values.chi;
   sol.kincns.dq    = values.dq;
   sol.mater.tau_c0 = values.tau_c0;
   use_plast = (sol.mater.m_digit==4 & (sol.mater.tau_c0>1e-10 & sol.mater.tau_c0<1e10));

   % retrieve wheel and rail profile data

   itask   = 1;
   iswheel = 0;
   isampl  = 0;
   iparam  = [iswheel, isampl];
   rparam  = [];
   val = cntc_getprofilevalues(ire, itask, iparam, rparam);
   sol.prr = struct('ProfileY',val(:,1), 'ProfileZ',val(:,2));
   itask   = 4;
   sol.prr.ProfileS = cntc_getprofilevalues(ire, itask, iparam, rparam);

   itask   = 1;
   iswheel = 1;
   iparam  = [iswheel, isampl];
   val = cntc_getprofilevalues(ire, itask, iparam, rparam);
   sol.prw = struct('ProfileY',val(:,1), 'ProfileZ',val(:,2));
   itask   = 4;
   sol.prw.ProfileS = cntc_getprofilevalues(ire, itask, iparam, rparam);

   % retrieve wheel, rail and contact position data

   ws_pos = cntc_getwheelsetposition(ire);
   cp_pos = cntc_getcontactlocation(ire, icp);
   sol.meta          = struct();
   sol.meta.s_ws     = ws_pos( 1);
   sol.meta.x_w      = cp_pos(21);
   sol.meta.y_w      = cp_pos(22);
   sol.meta.z_w      = cp_pos(23);
   sol.meta.roll_w   = ws_pos( 4);
   sol.meta.yaw_w    = ws_pos( 5);
   sol.meta.y_r      = cp_pos(26);
   sol.meta.z_r      = cp_pos(27);
   sol.meta.roll_r   = cp_pos(28);
   sol.meta.xcp_r    = cp_pos( 5);
   sol.meta.ycp_r    = cp_pos( 6);
   sol.meta.zcp_r    = cp_pos( 7);
   sol.meta.deltcp_r = cp_pos( 9);
   sol.meta.xcp_w    = cp_pos(10);
   sol.meta.ycp_w    = cp_pos(11);
   sol.meta.zcp_w    = cp_pos(12);
   sol.meta.rnom_whl =  500;
   sol.meta.rnom_rol = 1000;
   sol.meta.npatch   = cntc_getnumcontactpatches(ire);
   sol.meta.ipatch   = icp;

   % get discretization parameters

   [ mx, my, xc1, yc1, dx, dy ]  = cntc_getpotcontact(ire, icp);
   sol.mx = mx; sol.my = my;
   sol.xl = xc1 - 0.5*dx; sol.yl = yc1 - 0.5*dy;
   sol.dx = dx; sol.dy = dy;

   sol.x_offset = []; sol.y_offset = [];
   sol.x = sol.xl + ([1:mx]-0.5) * dx;
   sol.y = sol.yl + ([1:my]-0.5) * dy;

   % get grid-data

   sol.eldiv  = cntc_getelementdivision(ire, icp);
   sol.h      = cntc_getfielddata(ire, icp, CNTC.fld_h);
   sol.mu     = cntc_getfielddata(ire, icp, CNTC.fld_mu);
   sol.pn     = cntc_getfielddata(ire, icp, CNTC.fld_pn);
   sol.px     = cntc_getfielddata(ire, icp, CNTC.fld_px);
   sol.py     = cntc_getfielddata(ire, icp, CNTC.fld_py);
   sol.un     = cntc_getfielddata(ire, icp, CNTC.fld_un);
   sol.ux     = cntc_getfielddata(ire, icp, CNTC.fld_ux);
   sol.uy     = cntc_getfielddata(ire, icp, CNTC.fld_uy);
   sx         = cntc_getfielddata(ire, icp, CNTC.fld_sx);
   sy         = cntc_getfielddata(ire, icp, CNTC.fld_sy);
   sol.srel   = sqrt(sx.^2 + sy.^2);
   if (sol.kincns.t_digit>=2)
      sol.shft   = sol.srel * sol.kincns.dq;
   else
      sol.shft   = sol.srel;
   end
   if (use_plast)
      sol.taucrt     = cntc_getfielddata(ire, icp, CNTC.fld_taucrt);
      sol.uplsx      = cntc_getfielddata(ire, icp, CNTC.fld_uplsx);
      sol.uplsy      = cntc_getfielddata(ire, icp, CNTC.fld_uplsy);
      sol.trcbnd     = min(sol.mu .* sol.pn, sol.taucrt);
   else
      sol.trcbnd     = sol.mu .* sol.pn;
   end
   if (sol.h_digit>0)
      sol.temp1      = cntc_getfielddata(ire, icp, CNTC.fld_temp1);
      sol.temp2      = cntc_getfielddata(ire, icp, CNTC.fld_temp2);
   end

end % cntc_getcpresults

%------------------------------------------------------------------------------------------------------------

