!------------------------------------------------------------------------------------------------------------
! m_wr_output - print output to out-file for one case of w/r contact (module 1)
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------

module m_wr_output

use m_hierarch_data
use m_wrprof_data
use m_aijpj
use m_subsurf
implicit none
private

public  wr_output

contains

!------------------------------------------------------------------------------------------------------------

   subroutine wr_output ( wtd )
!--purpose: output routine for module 1: W/R contact
      implicit none
!--subroutine arguments:
      type(t_ws_track)         :: wtd
!--local variables:
      integer, parameter       :: idebug = 0
      logical      :: znewln, is_roll, is_ssrol, use_plast
      integer      :: icp, icpo, n_old, i, j, ii, iin, iout, iy0, iy1, ii0, ii1, lstrow, ncon, nadh,    &
                      nslip, nplast, nexter, ic_discns, is_right
      real(kind=8) :: def_turn, hmp, ptabs, ptarg, rhsx, rhsy, fxtru1, fytru1, sgn, rdum(20),           &
                      sc0, sc1, delt0, delt1, tmp1mx, tmp2mx, fac_in(nsens_in), fac_out(nsens_out)
      character(len= 9) :: str_muscal
      character(len=20) :: strng(max(20, nsens_in*nsens_out))
      character(len=80) :: str_prev
      type(t_marker)    :: whl_trk
      type(t_grid)      :: prr_glb, prw_glb
!--pointers to global data items:
      type(t_cpatch),    pointer :: cp
      type(t_probdata),  pointer :: gd
      type(t_metadata),  pointer :: meta
      character(len=5)           :: nam_side

      call timer_start(itimer_output)

      associate( ws    => wtd%ws,    trk  => wtd%trk, ic  => wtd%ic, mater => wtd%mater,                &
                 discr => wtd%discr, my_rail => wtd%trk%rai, my_wheel => wtd%ws%whl)

      is_roll      = (ic%tang.eq.2 .or. ic%tang.eq.3)
      is_ssrol     = (ic%tang.eq.3)
      use_plast    = (ic%mater.eq.4 .and. mater%tau_c0.gt.1d-10 .and. mater%tau_c0.le.1d10)

      ! set pointers to the active rail and wheel in the current configuration

      if (ic%is_left_side()) then
         nam_side = 'LEFT'
         is_right = 0
      else
         nam_side = 'RIGHT'
         is_right = 1
      endif

      ! set the sign to -1 for left and +1 for right rail/wheel combination

      sgn = 2 * is_right - 1

      ! Print the heading if "output" >= 2

      if (ic%output_surf.ge.2 .and. out_open.eq.1) then

         if (wtd%numcps.eq.1 .and. ic%config.ge.4) then
            write(lout, 1001) 'ROLLER', trim(nam_side), wtd%numcps, 'PATCH'
         elseif (wtd%numcps.eq.1) then
            write(lout, 1001) 'RAIL', trim(nam_side), wtd%numcps, 'PATCH'
         elseif (ic%config.ge.4) then
            write(lout, 1001) 'ROLLER', trim(nam_side), wtd%numcps, 'PATCHES'
         else
            write(lout, 1001) 'RAIL', trim(nam_side), wtd%numcps, 'PATCHES'
         endif
 1001    format(/,' WHEEL-',a,' CONTACT, ',a,' WHEEL,', i3,' CONTACT ',a,/)

         if (ic%gapwgt.eq.0) then
            ic_discns = ic%discns1_inp
         else
            ic_discns = ic%discns1_inp + 10*ic%gapwgt
         endif
         write(lout, 1101) ic%config, ic%pvtime, ic%bound, ic%tang, ic%norm, ic%force1, ic%stress
         write(lout, 1102) ic%varfrc, ic%frclaw_inp, ic_discns, ic%gencr_inp, ic%mater, ic%ztrack,      &
                        ic%ewheel
         write(lout, 1103) ic%heat, ic%gausei_inp, ic%iestim, ic%matfil_surf, ic%output_surf, ic%flow,  &
                        ic%return
         if (ic%xflow.ge.1) then
            write(lout, 1104) ic%x_profil, ic%x_smooth, ic%x_force, ic%x_locate, ic%x_readln,           &
                        ic%x_inflcf, ic%x_nmdbg
         endif
 1101    format (' CPBTNFS.  CONFIG',i2, ', PVTIME',i2, ', BOUND ',i2, ', TANG  ',i2, ', NORM  ',i2,    &
                        ', FORCE ',i2, ', STRESS',i2)
 1102    format (' VLDCMZE.  VARFRC',i2, ', FRCLAW',i2, ', DISCNS',i2, ', INFLCF',i2, ', MATER ',i2,    &
                        ', ZTRACK',i2, ', EWHEEL',i2)
 1103    format (' HGIAOWR.  HEAT  ',i2,', GAUSEI',i2, ', IESTIM',i2, ', MATFIL',i2, ', OUTPUT',i2,     &
                        ', FLOW  ',i2, ', RETURN',i2)
 1104    format (' PSFLRIN.  PROFIL',i2,', SMOOTH',i2, ', FORCE ',i2, ', LOCATE',i2, ', READLN',i2,     &
                        ', INFLCF',i2, ', NMDBG ',i2)

         if (ic%bound.eq.0 .and. ic%norm.eq.0) write(lout, 1201) 'POSITION',       'FULL'
         if (ic%bound.eq.0 .and. ic%norm.eq.1) write(lout, 1201) 'VERTICAL FORCE', 'FULL'
         if (ic%bound.eq.5 .and. ic%norm.eq.0) write(lout, 1201) 'POSITION',       'APPROXIMATE KPEC'
         if (ic%bound.eq.5 .and. ic%norm.eq.1) write(lout, 1201) 'VERTICAL FORCE', 'APPROXIMATE KPEC'
         if (ic%bound.eq.6 .and. ic%norm.eq.0) write(lout, 1201) 'POSITION',       'APPROXIMATE ANALYN'
         if (ic%bound.eq.6 .and. ic%norm.eq.1) write(lout, 1201) 'VERTICAL FORCE', 'APPROXIMATE ANALYN'
         if (ic%bound.eq.2) write(lout, 1203) 'ELLIPTICAL'
         if (ic%bound.eq.3) write(lout, 1203) 'PARABOLICAL'
         if (ic%bound.eq.4) write(lout, 1203) 'SDEC'
 1201    format (' VERTICAL PROBLEM WITH PRESCRIBED ',a,', ',a,' SOLUTION')
 1203    format (' VERTICAL PROBLEM SKIPPED, ',a,' TRACTION BOUND')

         if (ic%tang.ne.0) then
            if (ic%tang.eq.1) then
               strng(1) = 'SHIFT'
            else
               strng(1) = 'CREEPAGE'
            endif
               
            if (ic%force1.eq.0) then
               strng(2) = strng(1)
            elseif (ic%force1.eq.1) then
               strng(2) = 'X-FORCE'
            elseif (ic%force1.eq.2) then
               strng(2) = 'X-MOMENT'
            elseif (ic%force1.eq.3) then
               strng(2) = strng(1)
            endif

            if (ic%mater.ge.2 .and. ic%mater.le.3) then
               write(lout, 1300) trim(strng(2))
            else
               if (ic%tang.eq.1) write(lout, 1301) trim(strng(2))
               if (ic%tang.eq.2) write(lout, 1302) trim(strng(2))
               if (ic%tang.eq.3) write(lout, 1303) trim(strng(2))
            endif
 1300       format (' STEADY STATE ROLLING, "FASTSIM APPROACH", ',a,' PRESCRIBED')
 1301       format (' SHIFT TRANSIENT, ',a,' PRESCRIBED')
 1302       format (' TRANSIENT ROLLING, ',a,' PRESCRIBED')
 1303       format (' STEADY STATE ROLLING, ',a,' PRESCRIBED')

            if (wtd%mater%gencr_eff.eq.4 .and. mater%if_meth.eq.0) write(lout, 1400) mater%if_ver
            if (mater%gencr_eff.eq.4 .and. mater%if_meth.eq.1) write(lout, 1401) mater%if_ver
 1400       format (' USING CONFORMAL IF-CORRECTION: BLANCO(v',i1,'), CONSTANT CURVATURE')
 1401       format (' USING CONFORMAL IF-CORRECTION: BLANCO(v',i1,'), VARYING CURVATURE')

            if (wtd%solv%maxout.eq.1) write(lout, 1501)
            if (wtd%solv%maxout.gt.1) write(lout, 1502) wtd%solv%maxout
 1501       format (' JOHNSON"S METHOD IS USED FOR THE FRICTIONAL STRESS')
 1502       format (' PANAGIOTOPOULOS" METHOD IS USED FOR THE FRICTIONAL STRESS. MAXOUT=', I3)
         endif
      endif

      ! Print the input data for the W/R-problem if "output" >= 2.

      if (ic%output_surf.ge.2 .and. out_open.eq.1) then

         ! header per row: 2 spaces + [ 3 spaces + 9 characters ]*

         write(lout, 2100) mater%nu, fmt_gs(12,4,4,mater%ga), mater%ak, fmt_gs(12,4,4,wtd%solv%eps)
 2100    format (1x, /, ' MATERIAL CONSTANTS', /,                                                       &
            2x, 3x,'NU',7x, 3x,'G',8x, 3x,'AK',7x, 3x,'EPS',/, 2x, g12.4,a,g12.4,a,/)

         write(lout, 2101) mater%poiss, fmt_gs(12,4,4,mater%gg(1)), fmt_gs(12,4,4,mater%gg(2))
 2101    format ( 2x, 2x,'POISS(R)',2x, 2x,'POISS(W)',2x, 3x,'GG(R)',4x, 3x,'GG(W)',/, 2x, 2g12.4, 2a, /)

         if (ic%heat.ge.1 .and. ic%mater.ne.4) then
            write(lout, 2102) mater%bktemp(1), mater%heatcp(1), mater%lambda(1), mater%dens(1),         &
                              mater%bktemp(2), mater%heatcp(2), mater%lambda(2), mater%dens(2)
         elseif (ic%heat.ge.1) then
            write(lout, 2103) mater%bktemp(1), mater%heatcp(1), mater%lambda(1), mater%dens(1),         &
                mater%betapl, mater%bktemp(2), mater%heatcp(2), mater%lambda(2), mater%dens(2)
         endif
 2102    format (2x, 3x,'BKTEMP',6x, 'HEATCP',6x, 'LAMBDA',6x, 'DENSITY', /, 2(2x, 4g12.4, /))
 2103    format (2x, 3x,'BKTEMP',6x, 'HEATCP',6x, 'LAMBDA',6x, 'DENSITY',5x, 'BETAPL', /,               &
                 2x, 5g12.4, /, 2x, 4g12.4, /)

         if (ic%mater.eq.1) write(lout,2113) mater%nuv, mater%gav, mater%akv, mater%fg, &
                mater%tc, mater%vt
 2113    format (1x,'VISCOELASTIC MATERIAL CONSTANTS',/,                                                & 
            2x, 3x,'NUV',6x, 3x,'GV',7x, 3x,'AKV', /, 2x, 3g12.4,/,/,                                   &
            2x, 3x,'FG(R)',4x, 3x,'FG(W)',4x, 3x,'TC(R)',4x, 3x,'TC(W)',4x, 3x,'VT(R)',4x, 3x,'VT(W)',/, &
            2x, 6g12.4,/)

         if (ic%tang.gt.0 .and. ic%mater.eq.2) write(lout,2121) mater%flx(1)
 2121    format (1x,'SIMPLIFIED THEORY MATERIAL CONSTANTS',/, 2x, 3x,'FLX',/, 2x, g12.4,/)

         if (ic%tang.gt.0 .and. (ic%mater.eq.2 .or. ic%mater.eq.3))                                     &
                write(lout,2135) mater%k0_mf, mater%alfamf, mater%betamf
 2135    format (1x,'MODIFIED FASTSIM SLOPE REDUCTION PARAMETERS',/,                                    &
                 2x, 3x,'K0_MF',4x, 3x,'ALFAMF',3x, 3x,'BETAMF', /, 2x, 3g12.4,/)

         if (ic%tang.gt.0 .and. ic%mater.eq.4)                                                          &
                write(lout,2141) mater%gg3, mater%laythk, mater%tau_c0, mater%k_tau
 2141    format (1x,'INTERFACIAL LAYER PARAMETERS',/,                                                   &
            2x, 3x,'GG3',6x, 3x,'LAYTHK',3x, 3x,'TAU_C0',3x, 3x,'K_TAU',/, 2x, 4g12.4,/)

         if (ic%tang.gt.0) call fric_output(lout, ic%varfrc, ic%output_surf, wtd%fric)

         def_turn = 2d0 * discr%dist_sep - discr%dist_comb
         if (abs(discr%dist_turn-def_turn).lt.1d-6) then
            write(lout, 2201) discr%dx, discr%ds, discr%dqrel, discr%angl_sep, discr%dist_sep,          &
                discr%dist_comb
         else
            write(lout, 2202) discr%dx, discr%ds, discr%dqrel, discr%angl_sep, discr%dist_sep,          &
                discr%dist_comb, discr%dist_turn
         endif
 2201    format (1x, /, ' DISCRETISATION PARAMETERS', /,                                                &
            2x, 3x,'DX',7x, 3x,'DS',7x, 3x,'DQREL',4x, 3x,'A_SEP',4x, 3x,'D_SEP',4x, 3x,'D_COMB',/,     &
            2x, 6g12.4)
 2202    format (1x, /, ' DISCRETISATION PARAMETERS', /,                                                &
            2x, 3x,'DX',7x, 3x,'DS',7x, 3x,'DQREL',4x, 3x,'A_SEP',4x, 3x,'D_SEP',4x, 3x,'D_COMB',3x,    &
            3x,'D_TURN', /, 2x, 7g12.4)

         ! print information on profiles and smoothing

         if (ic%output_surf.ge.2 .and. ic%ilvout.ge.1 .and. ic%ztrack.ge.3) then

            if (ic%config.ge.4) then
               write(lout, '(/,3a)') ' ROLLER PROFILE "', trim(my_rail%prr%fname),'"'
            else
               write(lout, '(/,3a)') ' RAIL PROFILE "', trim(my_rail%prr%fname),'"'
            endif
            if (my_rail%prr%ismooth.eq.0) then
               write(lout, '(a)') '   MIRRORY     MIRRORZ     ISMOOTH      LAMBDA'
            else
               write(lout, '(a)') '   MIRRORY     MIRRORZ     ISMOOTH      L_FILT'
            endif
            write(lout, '(i7,2i12,3x,f12.4)') my_rail%prr%mirror_y, my_rail%prr%mirror_z,               &
                   my_rail%prr%ismooth, my_rail%prr%smth

         endif

         if (ic%output_surf.ge.3 .and. ic%ilvout.ge.1 .and. ic%ztrack.ge.3) then

            write(lout, '(/, a)') '   ZIGTHRS     KINKHIG     KINKLOW      KINKWID     MAXOMIT'
            write(lout, '(2x, 5g12.4)') my_rail%prr%zig_thrs, my_rail%prr%kink_high,                    &
                   my_rail%prr%kink_low, my_rail%prr%kink_wid, my_rail%prr%f_max_omit

         endif

         if (ic%output_surf.ge.2 .and. ic%ilvout.ge.1 .and. ic%ewheel.ge.3) then

            write(lout, '(/,3a)') ' WHEEL PROFILE "', trim(my_wheel%prw%fname),'"'
            if (my_wheel%prw%ismooth.eq.0) then
               write(lout, '(a)') '   MIRRORY     MIRRORZ     ISMOOTH      LAMBDA'
            else
               write(lout, '(a)') '   MIRRORY     MIRRORZ     ISMOOTH      L_FILT'
            endif
            write(lout, '(i7,2i12,3x,f12.4)') my_wheel%prw%mirror_y, my_wheel%prw%mirror_z,             &
                   my_wheel%prw%ismooth, my_wheel%prw%smth

         endif

         if (ic%output_surf.ge.3 .and. ic%ilvout.ge.1 .and. ic%ewheel.ge.3) then

            write(lout, '(/, a)') '   ZIGTHRS     KINKHIG     KINKLOW      KINKWID     MAXOMIT'
            write(lout, '(2x, 5g12.4)') my_wheel%prw%zig_thrs, my_wheel%prw%kink_high,                  &
                   my_wheel%prw%kink_low, my_wheel%prw%kink_wid, my_wheel%prw%f_max_omit

         endif

         write(lout, 4000)
 4000    format (/, 1x,'WHEEL-SET POSITION AND VELOCITY')

         if (.false.) then
            strng(1) = fmt_gs(18,10,10,ws%z)
         else
            write(strng(1), '(f12.4)') ws%z
         endif
         if (ic%config.le.1) then
            ! wheelset on track: print S_WS, VX_WS
            write(lout, 4001) 'S', fmt_gs(12,4,4,ws%s), fmt_gs(12,4,4,ws%y), trim(strng(1)),            &
                   fmt_gs(12,4,4,ws%roll), fmt_gs(12,4,4,ws%yaw), fmt_gs(12,4,4,ws%pitch)
            write(lout, 4002) fmt_gs(12,4,4,ws%vs), fmt_gs(12,4,4,ws%vy), fmt_gs(12,4,4,ws%vz),               &
                   fmt_gs(12,4,4,ws%vroll), fmt_gs(12,4,4,ws%vyaw), ws%vpitch
         else
            ! wheelset on roller rig: print X_WS, VPITCH_ROL
            write(lout, 4001) 'X', fmt_gs(12,4,4,ws%x), fmt_gs(12,4,4,ws%y), trim(strng(1)),            &
                   fmt_gs(12,4,4,ws%roll), fmt_gs(12,4,4,ws%yaw), fmt_gs(12,4,4,ws%pitch)
            write(lout, 4003) trk%vpitch_rol, fmt_gs(12,4,4,ws%vy), fmt_gs(12,4,4,ws%vz),               &
                  fmt_gs(12,4,4,ws%vroll), fmt_gs(12,4,4,ws%vyaw), ws%vpitch
         endif
 4001    format (2x, 3x,a1,'_WS',5x, 3x,'Y_WS',5x, 3x,'Z_WS',5x, 3x,'ROLL_WS',2x, 3x,'YAW_WS',3x,       &
            3x,'PITCH_WS',1x,/, 2x, 6a /)
 4002    format( 2x, 3x,'VX_WS',4x, 3x,'VY_WS',4x, 3x,'VZ_WS',4x, 3x,'VROLL_WS',1x, 3x,'VYAW_WS',2x,    &
            3x,'VPITCH_WS',/, 2x, 5a, f12.7 /)
 4003    format( 2x, 2x,'VPITCH_ROL', 3x,'VY_WS',4x, 3x,'VZ_WS',4x, 3x,'VROLL_WS',1x, 3x,'VYAW_WS',2x,  &
            3x,'VPITCH_WS',/, 2x, f12.7, 4a, f12.7 /)

         rdum(1:12) =  (/ my_wheel%dx, my_wheel%dy, my_wheel%dz, my_wheel%droll, my_wheel%dyaw,         &
                          my_wheel%dpitch, my_wheel%vx, my_wheel%vy, my_wheel%vz, my_wheel%vroll,       &
                          my_wheel%vyaw, my_wheel%vpitch /)
         if (maxval(abs(rdum(1:12))).gt.1d-10) then
            write(lout, 4101) fmt_gs(12,4,4,rdum( 1)), fmt_gs(12,4,4,rdum( 2)), fmt_gs(12,4,4,rdum( 3)), &
                              fmt_gs(12,4,4,rdum( 4)), fmt_gs(12,4,4,rdum( 5)), fmt_gs(12,4,4,rdum( 6)), &
                              fmt_gs(12,4,4,rdum( 7)), fmt_gs(12,4,4,rdum( 8)), fmt_gs(12,4,4,rdum( 9)), &
                              fmt_gs(12,4,4,rdum(10)), fmt_gs(12,4,4,rdum(11)), fmt_gs(12,4,4,rdum(12))
 4101       format (1x,'FLEXIBLE WHEEL-SET DEVIATIONS',/,                                               &
               2x, 3x,'DXWHL',4x, 3x,'DYWHL',4x, 3x,'DZWHL',4x, 3x,'DROLLW',3x, 3x,'DYAWW',4x,          &
                   3x,'DPITCHW',2x, /, 2x, 6a, /,                                                       &
               2x, 3x,'VXWHL',4x, 3x,'VYWHL',4x, 3x,'VZWHL',4x, 3x,'VROLLW',3x, 3x,'VYAWW',4x,          &
                   3x,'VPITCHW',2x, /, 2x, 6a, /)
         endif

         if (max(abs(my_rail%dy), abs(my_rail%dz), abs(my_rail%roll), abs(my_rail%vy), abs(my_rail%vz), &
                 abs(my_rail%vroll)).gt.1d-10) then
            write(lout, 4102) my_rail%dy, my_rail%dz, my_rail%roll, fmt_gs(12,4,4,my_rail%vy),          &
                              fmt_gs(12,4,4,my_rail%vz), fmt_gs(12,4,4,my_rail%vroll)
 4102       format (1x,'RAIL IRREGULARITY',/,                                                           &
               2x, 3x,'DYRAIL',3x, 3x,'DZRAIL',3x, 3x,'DROLLR',3x, 3x,'VYRAIL',3x, 3x,'VZRAIL',3x,      &
                   3x,'VROLLR',3x, /, 2x, 3f12.4, 3a,/)
         endif

         if (ic%force1.eq.3) then
            write(lout,4202)
            write(lout,4203) wtd%trk%ky_rail, wtd%trk%dy_defl, fmt_gs(12,4,4,wtd%trk%fy_rail),          &
                   wtd%trk%kz_rail, wtd%trk%dz_defl, fmt_gs(12,4,4,wtd%trk%fz_rail)
 4202       format (1x,'MASSLESS RAIL DEFLECTION')
 4203       format (2x, 3x,'KY_RAIL',2x, 3x,'DY_DEFL',2x, 3x,'FY_RAIL',2x, 3x,'KZ_RAIL',2x,             &
                   3x,'DZ_DEFL',2x, 3x,'FZ_RAIL',2x, /, 2x, 2f12.4,a, 2f12.4,a,/)
         endif

      endif

      ! Print the overall results for the W/R-problem if "output" >= 1.

      if (ic%output_surf.ge.1 .and. out_open.eq.1) then

         if (wtd%numcps.le.0) then
            write(lout, 5000)
 5000       format(' ----- NO CONTACT FOUND FOR THIS WHEEL / WHEEL-SET -----')
         endif

         if (ic%output_surf.ge.2 .and. ic%config.le.1) then
            write(lout,5010) 'RAIL'
         elseif (ic%output_surf.ge.2) then
            write(lout,5010) 'ROLLER'
         endif
         write(lout,5011)
         write(lout,5012) fmt_gs(12,4,4,wtd%ftrk%x()), fmt_gs(12,4,4,sgn*wtd%ftrk%y()),                 &
                fmt_gs(12,4,4,wtd%ftrk%z()), fmt_gs(12,4,4,wtd%fws%x()), fmt_gs(12,4,4,sgn*wtd%fws%y()), &
                fmt_gs(12,4,4,wtd%fws%z())
!               fmt_gs(16,8,8,wtd%ftrk%z()), fmt_gs(16,8,8,wtd%fws%x()), fmt_gs(12,4,4,sgn*wtd%fws%y()), &
 5010    format (1x, /, ' TOTAL FORCES AND MOMENTS ON ',a)
 5011    format (2x, 3x,'FX(TR)',3x, 3x,'FY(TR)',3x, 3x,'FZ(TR)',3x, 3x,'FX(WS)',3x, 3x,'FY(WS)',3x,    &
                     3x,'FZ(WS)')
 5012    format (6a)
      endif

      if (ic%output_surf.ge.2 .and. out_open.eq.1) then
         write(lout,5102)
         write(lout,5103) wtd%xavg%x(), sgn*wtd%xavg%y(), wtd%xavg%z(), fmt_gs(12,4,4,sgn*wtd%tavg%x()), &
                fmt_gs(12,4,4,wtd%tavg%y()), fmt_gs(12,4,4,sgn*wtd%tavg%z())
 5102    format (1x, /, ' AVERAGE CONTACT POSITION')
 5103    format (2x, 3x,'XAV(TR)',2x, 3x,'YAV(TR)',2x, 3x,'ZAV(TR)',2x, 3x,'MX(AV)',3x,                 &
                     3x,'MY(AV)',3x, 3x,'MZ(AV)', /, 3f12.4, 3a)

      endif
   
      ! Loop over all contact patches

      do icp = 1, wtd%numcps

         cp     => wtd%allcps(icp)%cp
         gd     => cp%gd
         meta   => gd%meta

         associate(cgrid  => gd%cgrid_cur,    muscal => gd%kin%muscal,    pen    => gd%kin%pen,         &
                   cksi   => gd%kin%cksi,     ceta   => gd%kin%ceta,      cphi   => gd%kin%cphi,        &
                   fntrue => gd%kin%fntrue,   fnscal => gd%kin%fnscal,    fxrel1 => gd%kin%fxrel1,      &
                   fyrel1 => gd%kin%fyrel1,   chi    => gd%kin%chi,       veloc  => gd%kin%veloc,       &
                   dq     => gd%kin%dq,       mus1   => gd%outpt1%mus,    igs1   => gd%outpt1%igs,      &
                   ps1    => gd%outpt1%ps,    shft1  => gd%outpt1%shft,   mztru1 => gd%outpt1%mztrue,   &
                   elen1  => gd%outpt1%elen,  temp1  => gd%outpt1%temp1,  temp2  => gd%outpt1%temp2,    &
                   frpow1 => gd%outpt1%frpow, pmax1  => gd%outpt1%pmax,   eps    => gd%solv%eps,        &
                   hs1    => gd%geom%hs1,     subs   => gd%subs )

         if (ic%output_surf.ge.1 .and. out_open.eq.1) then
            n_old = 0
            do icpo = 1, MAX_NUM_CPS
               if (cp%prev_icp(icpo).ge.1) n_old = n_old + 1
            enddo
            if (n_old.le.0) then
               str_prev = 'NEW PATCH'
            else
               write(str_prev, '(a,10(i3,:,'',''))') 'PREV', (cp%prev_icp(icpo), icpo=1,n_old)
            endif
            if (ic%pvtime.eq.2 .and. ic%iestim.eq.0) then
               write(lout, 6000) icp                    ! no sequence, no initial estimate
            else
               write(lout, 6001) icp, trim(str_prev)    ! sequence or initial estimate used
            endif
 6000       format(/,' ----- DATA FOR CONTACT PATCH',i3,' -----',/)
 6001       format(/,' ----- DATA FOR CONTACT PATCH',i3,' (',a,') -----',/)
         endif

         ! Print additional input data for the contact patch if "output" >= 2.

         if (ic%output_surf.ge.2 .and. out_open.eq.1) then

            ! write output for contact reference position on track / rail profile
            ! note: undo mirroring for track y but not for rail y

            write(lout, 6100)
            ! write(lout, 6101) 'RAIL'
            write(lout, 6102) meta%xcp_tr, sgn*meta%ycp_tr, meta%zcp_tr, sgn*meta%deltcp_tr,            &
                sgn*meta%ycp_r, meta%zcp_r

            ! write output for contact reference position on wheelset / wheel profile

            ! write(lout, 6101) 'WHEEL'
            write(lout, 6103) meta%xcp_w, sgn*meta%ycp_w, meta%zcp_w

 6100       format(1x,'CONTACT REFERENCE LOCATION')
 6101       format(1x,'CONTACT POINT LOCATION ON ',a)
 6102       format(2x, 3x,'XCP(TR)',2x, 3x,'YCP(TR)',2x, 3x,'ZCP(TR)',2x, 2x,'DELTCP(TR)',1x,           &
                   2x,'YCP(R)',3x, 3x,'ZCP(R)', /, 6f12.4, /)
 6103       format(2x, 3x,'XCP(W)',3x, 3x,'YCP(W)',3x, 3x,'ZCP(W)', /, 3f12.4, /)

            ! report on angle variation within the actual contact, if on-empty

            iy0 = igs1%iymin
            iy1 = igs1%iymax

            if (ic%is_conformal() .and. iy1.ge.iy0) then
               ii0 = 1 + (iy0-1) * cgrid%nx
               ii1 = 1 + (iy1-1) * cgrid%nx
               sc0 = cgrid%y( ii0 ) - 0.5d0*cgrid%dy
               sc1 = cgrid%y( ii1 ) + 0.5d0*cgrid%dy
               delt0 = atan( cp%curv_nrm%vy(iy0) / -cp%curv_nrm%vn(iy0) ) * 180d0 / pi
               delt1 = atan( cp%curv_nrm%vy(iy1) / -cp%curv_nrm%vn(iy1) ) * 180d0 / pi
               write(bufout,'(3(a,i3),2(a,f7.3),a)') ' Actual contact on columns iy = [',iy0,',',       &
                            iy1,'] (my=',cgrid%ny,'), sc = [', sc0,',', sc1,']'
               call write_log(1, bufout)
               write(bufout,'(3(a,f7.2),a)') ' Reference contact angle delttr =', cp%delttr*180d0/pi,   &
                            ' deg, range = [',delt0,',',delt1,'] deg'
               call write_log(1, bufout)

               write(bufout,'(3(a,f12.1),a)') ' Total moments Mxp =', gd%outpt1%mxtrue,', Msp =',       &
                        gd%outpt1%mytrue, ', Mnp =', gd%outpt1%mztrue, ' N.mm'
               call write_log(1, bufout)
               call write_log(' ')
            endif

            ! V = 1: report on actual friction parameters used, converted to struct with NVF = 1

            if (ic%tang.gt.0 .and. wtd%ic%varfrc.eq.1) then

               call fric_output(lout, gd%ic%varfrc, wtd%ic%output_surf, gd%fric)
               write(lout,'(1x)')

            endif

            !  M = 3: report actual flexibilities, M = 2, 3: report on actual K_EFF

            if (ic%tang.gt.0 .and. ic%mater.eq.3) then
               write(lout,6231) gd%mater%flx(1), gd%mater%flx(2), gd%mater%flx(3)
 6231          format (1x,'SIMPLIFIED THEORY MATERIAL CONSTANTS',/,                                     &
                 2x, 3x,'FLX1',5x, 3x,'FLX2',5x, 3x,'FLX3',/, 2x, 3g12.4,/)
            endif

            if (ic%tang.gt.0 .and. (ic%mater.eq.2 .or. ic%mater.eq.3)) then
               write(lout,6251) gd%mater%k_eff
 6251          format (1x,'MODIFIED FASTSIM SLOPE REDUCTION PARAMETERS',/,                              &
                 2x, 3x,'K0_EFF',4x, /, 2x, 1g12.4,/)
            endif

            ! report on creepages

            if (ic%tang.gt.0) then
               strng(1) = fmt_gs(12, 4, 4,     veloc)
               strng(2) = fmt_gs(12, 4, 4,     cksi)
               strng(3) = fmt_gs(12, 4, 4, sgn*ceta)
               strng(4) = fmt_gs(12, 4, 4, sgn*cphi)
               write(lout, 6301)
               write(lout, 6305) sgn*chi, dq, strng(1), strng(2), strng(3), strng(4)
 6301          format (' KINEMATIC CONSTANTS')
 6305          format (2x, 3x,'CHI',6x, 3x,'DQ',7x,    3x,'VELOC',4x, 3x,'CKSI',5x, 3x,'CETA',5x,       &
                    3x,'CPHI',/, 2x, 2g12.4, 4a12, /)

               if (max(abs(gd%kin%spinxo),abs(gd%kin%spinyo)).ge.tiny) then
                  write(lout, 6311) gd%kin%spinxo, gd%kin%spinyo
               endif
 6311          format (2x, 3x,'SPINXO',3x,  3x,'SPINYO',/, 2x, 2g12.4, /)
            endif
         endif

         ! Print the main results for the patch: total forces, creepages

         if (ic%output_surf.ge.1 .and. out_open.eq.1) then

            ! Filter values that are dominated by noise using filt_sml()

            fxtru1 = fxrel1 * (fntrue*muscal+tiny)
            fytru1 = fyrel1 * (fntrue*muscal+tiny)

            strng(1) = fmt_gs(12, 4, 4, fntrue)
            strng(2) = fmt_gs(12, 4, 4, filt_sml(    fxtru1,0.5d0*eps*fntrue*muscal))
            strng(3) = fmt_gs(12, 4, 4, filt_sml(sgn*fytru1,0.5d0*eps*fntrue*muscal))
            strng(4) = fmt_gs(12, 4, 4, sgn*mztru1)
            strng(5) = fmt_gs(12, 4, 4, elen1)
            strng(6) = fmt_gs(12, 4, 4, frpow1)

            if (ic%output_surf.ge.2) write(lout,6400)
            write(lout,6401)
            write(lout,6402) (strng(j), j=1,6)
 6400       format (1x, /, ' TOTAL FORCES, TORSIONAL MOMENT, ELASTIC ENERGY, FRICTIONAL POWER')
 6401       format (2x, 3x,'FN',7x, 3x,'FX',7x, 3x,'FS',7x, 3x,'MN',7x, 2x,'ELAST.EN.',1x, 2x,'FRIC.POWER')
 6402       format ( 2x, 6a12)

            strng(1) = 'FX/FSTAT/FN'
            strng(2) = 'FS/FSTAT/FN'
            if (.not.gd%kin%use_muscal) strng(1) = '  FX/FN'
            if (.not.gd%kin%use_muscal) strng(2) = ' FS/FN'

            strng(3) = fmt_gs(12, 4, 4, fnscal)
            strng(4) = fmt_gs(12, 4, 4, filt_sml(fxrel1,0.5d0*eps))
            strng(5) = fmt_gs(12, 4, 4, filt_sml(sgn*fyrel1,0.5d0*eps))
            strng(6) = fmt_gs(12, 4, 4, gd%kin%pen)

            if (ic%print_pmax) then
               strng(7)  = '    PMAX    '
               strng(9)  = fmt_gs(12, 4, 4, pmax1)
            endif
            if (ic%heat.ge.1) then
               tmp1mx = gf3_max(AllElm, temp1, ikZDIR)
               tmp2mx = gf3_max(AllElm, temp2, ikZDIR)
               if (.not.ic%print_pmax) then
                  strng(7)  = '     MAX(T1)'
                  strng(8)  = '     MAX(T2)'
                  strng(9)  = fmt_gs(12, 4, 4, tmp1mx)
                  strng(10) = fmt_gs(12, 4, 4, tmp2mx)
               else
                  strng(8)  = '  MAX(T1,T2)'
                  strng(10) = fmt_gs(12, 4, 4, max(tmp1mx,tmp2mx))
               endif
            endif

            if (ic%heat.ge.1) then
               write(lout,6411) (strng(j), j=1,2), strng(7), strng(8)
               write(lout,6412) (strng(j), j=3,6), strng(9), strng(10)
            elseif (ic%print_pmax) then
               write(lout,6411) (strng(j), j=1,2), strng(7)
               write(lout,6412) (strng(j), j=3,6), strng(9)
            else
               write(lout,6411) (strng(j), j=1,2)
               write(lout,6412) (strng(j), j=3,6)
            endif
 6411       format('     FN/G      ',a12,1x,a12,' APPROACH',:,1x,2a12)
 6412       format(2x, 6a12)

         endif

         ! Print more detailed global output when "output" >= 2

         if (ic%output_surf.ge.2 .and. ic%ilvout.ge.2) then

            ! Set mirroring of sensitivities for left-side wheel/rail

            fac_in(1:nsens_in) = 1d0
            fac_in(iin_deta1)  = sgn 
            fac_in(iin_dphi1)  = sgn

            fac_out(1:nsens_out) = 1d0
            fac_out(iout_fy1)  = sgn
            fac_out(iout_mz1)  = sgn

            ! Convert sensitivities to strings, row-wise

            do iout = 1, nsens_out
               do iin = 1, nsens_in
                  i = nsens_in*(iout-1) + iin
                  strng(i) = fmt_gs(12, 3, 3, fac_in(iin) * fac_out(iout) * gd%outpt1%sens(iout,iin))
               enddo
            enddo

            if (ic%tang.ne.0) then

               ! sensitivities for Norm and Tang, one contact-problem

               str_muscal = 'FN/FSTAT.'
               if (.not.gd%kin%use_muscal) str_muscal = 'FN.'
               write(lout, 6501) str_muscal
 6501          format(/,/, ' THE SENSITIVITIES. (FX,FY) MEANS: (FX,FY)/',a,/,                           &
                           ' A ZERO ENTRY MEANS THAT IT HAS NOT BEEN CALCULATED.',/,                    &
                       8x, ' DF/DPEN     DF/DKSI     DF/DETA     DF/DPHI',/)
               if (ic%sens.ge.2) then
                  write(lout, 6502) (strng(i),i=1,4), (strng(i),i=8,11), (strng(i),i=15,18),            &
                        (strng(i),i=22,25)
               else
                  write(lout, 6502) (strng(i),i=1,4), (strng(i),i=8,11), (strng(i),i=15,18)
               endif
 6502          format (3x,'FN ',4a12,/, 3x,'FX ',4a12,/, 3x,'FY ',4a12,:,/, 3x,'MZ ',4a12)

            else

               ! sensitivities for Norm only

               write(lout, 6511) strng(1)
 6511          format (/,/, ' THE SENSITIVITIES.', /, ' A ZERO MEANS THAT IT HAS NOT BEEN CALCULATED.', /, &
                   2x, 'DFN/DPEN', /, a12)

            endif
         endif

         ! write the statistics when "output" >= 2

         if (ic%output_surf.ge.2) then
            call eldiv_count(igs1, nadh, nslip, nplast, nexter)
            ncon = nadh + nslip + nplast
            if (use_plast) then
               write(lout, 6601) '  NPLAST', gd%potcon_cur%npot, ncon, nadh, nslip, nplast,             &
                        gd%solv%itnorm, gd%solv%ittang
            else
               write(lout, 6601) ' ', gd%potcon_cur%npot, ncon, nadh, nslip, gd%solv%itnorm, gd%solv%ittang
            endif
 6601       format (/, ' CONTACT STATISTICS', /,                                                        &
                       ' N: NUMBER OF ELEMENTS IN REGION,  I: NUMBER OF ITERATIONS,', /,                &
                       ' POT: POTENTIAL CONTACT AREA,  CON: CONTACT,  ADH: ADHESION', /,                &
                       '    NPOT   NCON   NADH  NSLIP',a,' INORM  ITANG', /, 7i7)
         endif

         ! write the picture of the contact area when "output" >= 3

         if (ic%output_surf.ge.3 .and. ic%ilvout.ge.1) then
            write(lout, 6700)
 6700       format(/, ' FORM OF THE CONTACT, ADHESION AND SLIP-AREA`S:')
            call wrigs (igs1, is_roll, gd%kin%chi)
         endif

         ! print the curved reference surface used for conformal contact when "output" >= 3 and D == 4

         if (ic%output_surf.ge.3 .and. ic%is_conformal() .and. ic%ilvout.ge.1) then

            write(lout, '(/,a,i5,a)') ' CURVED REFERENCE SURFACE, TRACK COORDINATES,', cp%curv_ref%ny, &
                ' POINTS:'
            write(lout, '(a)') '    I        Y(I)       Z(I)       A(I)'

            do ii = 1, cp%curv_ref%ny
               write(lout,6800) ii, sgn*cp%curv_ref%y(ii), cp%curv_ref%z(ii), sgn*cp%curv_incln%vy(ii)
 6800          format(i6,4(:,' ',f10.4))
            enddo

         endif

         ! Print the detailed solution for all elements, when "output" >= 5

         if (ic%output_surf.ge.5) then

            if (ic%tang.eq.0) then
               write(lout, 6900)
 6900          format (//,                                                                              &
                    ' QUANTITIES DISTRIBUTED OVER THE CONTACT AREA:', /,                                &
                    ' "H-PEN" IS THE UNDEFORMED DISTANCE WITH THE APPROACH TAKEN INTO ACCOUNT.',/,      &
                    ' PN IS THE NORMAL PRESSURE, PT IS THE TANGENTIAL TRACTION (2-VECTOR).',/,          &
                    ' X, H AND PEN ARE IN MM, PN AND PT ARE IN N/MM^2, AND ARGUMENTS ARE IN DEGREES.',//,&
                    5x,'X',10x,'H-PEN',8x, 'PN', 7x, 'ABS(PT)', 8x, 'ARG(PT)')

            elseif (.not.is_roll) then
               write(lout, 6901)
 6901          format (//,                                                                              &
                    ' QUANTITIES DISTRIBUTED OVER THE CONTACT AREA:', /,                                &
                    ' "H-PEN" IS THE UNDEFORMED DISTANCE WITH THE APPROACH TAKEN INTO ACCOUNT.',/,      &
                    ' PN IS THE NORMAL PRESSURE, PT THE TANGENTIAL TRACTION (2-VECTOR).',/,             &
                    ' S IS THE SLIP DISTANCE PER STEP, WHICH IS OPPOSITE TO THE TANGENTIAL TRACTION.',/,&
                    ' RIGID SHIFT(X, Y) ARE THE (X, Y) COMPONENTS OF THE RIGID SHIFT DISTANCE W.',/,    &
                    ' X, H, PEN, S AND W ARE IN MM, PN AND PT ARE IN N/MM^2, AND ARGUMENTS ARE IN DEGREES.'/,&
                    //,                                                                                 &
                    5x,'X',9x,'H-PEN',8x,'PN',6x,'TRCBND    ABS(PT)     ARG(PT;-S)',2x,'ABS(S)',6x,     &
                           'RIGID SHIFT(X, Y)')

            else
               write(lout, 6902)
 6902          format (//,                                                                              &
                    ' QUANTITIES DISTRIBUTED OVER THE CONTACT AREA:', /,                                &
                    ' "H-PEN" IS THE UNDEFORMED DISTANCE WITH THE APPROACH TAKEN INTO ACCOUNT.',/,      &
                    ' PN IS THE NORMAL PRESSURE, PT THE TANGENTIAL TRACTION (2-VECTOR).',/,             &
                    ' S IS THE RELATIVE SLIP VELOCITY, WHICH IS OPPOSITE TO THE TANGENTIAL TRACTION.',/ &
                    ' RIGID SLIP(X, Y) ARE THE (X, Y) COMPONENTS OF THE RELATIVE RIGID SLIP VELOCITY W.',/&
                    ' X, H AND PEN ARE IN MM, PN AND PT ARE IN N/MM^2, S AND W ARE DIMENSIONLESS, ',    &
                                                                       'ARGUMENT IS IN DEGREES.',//,    &
                    5x,'X',9x,'H-PEN',8x,'PN',6x,'TRCBND    ABS(PT)     ARG(PT;-S)',2x,'ABS(S)',6x,     &
                           'RIGID SLIP(X, Y)')

            endif

            zNewLn = .true.
            LstRow = cgrid%iy(1)

            ! for all elements inside the contact area do

            do ii = 1, gd%potcon_cur%npot
               if (cgrid%iy(ii).ne.lstrow) znewln = .true.
               lstrow = cgrid%iy(ii)
               if (igs1%el(ii).ge.Adhes .or. ic%output_surf.ge.6) then

                  ! print header when starting a new row of the potential contact:

                  if (znewln) then
                     write(lout, 6910) cgrid%y(ii), cgrid%iy(ii)
                     znewln = .false.
 6910                format (/, ' Y = ', g12.4, ' ROW', i3, ' OF THE POTENTIAL CONTACT', /)
                  endif

                  ! retrieve/determine output quantities

                  ptabs = dsqrt (ps1%vx(ii) **2 + ps1%vy(ii) **2)
                  ptarg = atan2(ps1%vy(ii), ps1%vx(ii)) * 180d0/pi
                  hmp = hs1%vn(ii) - pen
                  if (ic%tang.ne.0) then

                     ! module 1 uses F3 = 0
                     ! hs, ss  == (rigid) slip distance per (time-)step
                     ! rolling: rhs, creepage == relative slip velocity
                     ! shifts:  rhs, creepage == slip distance, dq == 1

                     rhsx = - hs1%vx(ii) / dq
                     rhsy = - hs1%vy(ii) / dq

                     write(lout, 6920) cgrid%x(ii), hmp, ps1%vn(ii), ps1%vn(ii)*mus1%vt(ii),          &
                                ptabs, ptarg, shft1%vt(ii)/dq, rhsx, rhsy
 6920                format (1x, g12.4, 4g11.3, f10.1, 2x, 3g11.3)
                  else
                     write(lout, 6930) cgrid%x(ii), hmp, ps1%vn(ii), ptabs, ptarg
 6930                format (1x, 4g12.4, f10.1)
                  endif
               endif
            enddo ! all elements
         endif ! O >= 5

         ! if S>=1, print subsurface data for patch icp

         if (ic%stress.ge.1) then

            call subsur_outfile(ic, subs, idebug)

         endif

         end associate ! reference to elements of gd on allcps(icp)

      enddo ! numcps

      ! print the profiles used in track coordinates when "output" >= 4

      if (ic%output_surf.ge.4 .and. ic%ilvout.ge.1) then

         ! get the rail profile in track coordinates

         call grid_copy(my_rail%prr%grd_data, prr_glb, with_spline=.true.)
         call cartgrid_2glob(prr_glb, my_rail%m_trk)

         write(lout, '(/,a,i5,a)') ' RAIL PROFILE USED, TRACK COORDINATES,',prr_glb%ny,' POINTS:'
         write(lout, '(a)') '    I        S(I)       Y(I)       Z(I)'

         do ii = 1, prr_glb%ny
            write(lout,8001) ii, prr_glb%s_prf(ii), sgn*prr_glb%y(ii), prr_glb%z(ii)
 8001       format(i6,4(:,' ',f10.4))
         enddo

         ! get the wheel profile in track coordinates

         whl_trk = marker_2glob( my_wheel%m_ws, ws%m_trk )
         call grid_copy(my_wheel%prw%grd_data, prw_glb, with_spline=.true.)
         call cartgrid_2glob(prw_glb, whl_trk)

         write(lout, '(/,a,i5,a)') ' PRINCIPAL WHEEL PROFILE, TRACK COORDINATES,', prw_glb%ny,' POINTS:'
         write(lout, '(a)') '    I        S(I)       X(I)       Y(I)       Z(I)'

         do ii = 1, prw_glb%ny
            if (.not.prw_glb%lies_in_oyz) then
               write(lout,8001) ii, prw_glb%s_prf(ii), prw_glb%x(ii), sgn*prw_glb%y(ii), prw_glb%z(ii)
            else
               write(lout,8001) ii, prw_glb%s_prf(ii), 0d0, sgn*prw_glb%y(ii), prw_glb%z(ii)
            endif
         enddo

         call grid_destroy(prr_glb)
         call grid_destroy(prw_glb)
      endif

      end associate

      call timer_stop(itimer_output)
   end subroutine wr_output

!------------------------------------------------------------------------------------------------------------

end module m_wr_output
