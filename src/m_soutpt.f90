!============================================================================================================
! m_soutpt - print output to out-file and/or mat-file
!
! Copyright 1979-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_soutpt

use m_hierarch_data
use m_aijpj
implicit none

private
public  calc_damping_force
public  soutpt
public  writmt

contains

!------------------------------------------------------------------------------------------------------------

   subroutine calc_damping_force( ic, mater, kin )
!--purpose: compute ad-hoc proportional damping
      implicit none
!--subroutine arguments:
      type(t_ic)                  :: ic
      type(t_material)            :: mater
      type(t_kincns)              :: kin
!--local variables :
      real(kind=8) :: df(3), dt_min
      character(len=12) :: strng(3)

      associate(cdampn => mater%cdampn,  cdampt => mater%cdampt,  dt     => kin%dt,                     &
                fcntc  => kin%fcntc,     fprev  => kin%fprev,     fdamp  => kin%fdamp )

      ! compute ad-hoc proportional damping

      if (ic%pvtime.eq.2) then
         fdamp  = (/ 0d0, 0d0, 0d0 /)
      else
         dt_min = 1d-20
         df     = (fcntc - fprev) / max(dt_min, dt)
         if (abs(df(1)).gt.mater%dftmax) df(1) = sign(mater%dftmax, df(1)) ! |dftmax| * sign(df(1))
         if (abs(df(2)).gt.mater%dftmax) df(2) = sign(mater%dftmax, df(2))
         if (abs(df(3)).gt.mater%dfnmax) df(3) = sign(mater%dfnmax, df(3))
         fdamp  = (/ cdampt * df(1), cdampt * df(2), cdampn * df(3) /)

         if (ic%x_force.ge.2) then
            ! write(bufout,'(2(a,g12.4),a,3g12.4)') ' cdampn=',cdampn,', dt=',dt, ', fdamp=',fdamp
            strng(1) = fmt_gs(12, 9, 6, fprev(3))
            strng(2) = fmt_gs(12, 9, 6, fcntc(3))
            strng(3) = fmt_gs(12, 9, 6, fdamp(3))
            write(bufout,'(2(a,g12.4),6a)') ' cdampn=',cdampn,', dt=',dt,', fprev=',strng(1),           &
                           ', fcntc=',strng(2),', fdamp=',strng(3)
            call write_log(1, bufout)
         endif
      endif

      end associate
   end subroutine calc_damping_force

!------------------------------------------------------------------------------------------------------------

   subroutine soutpt ( gd )
!--purpose: compute derived quantities, print output to out-file
      implicit none
!--subroutine arguments:
      type(t_probdata)         :: gd
!--local variables :
      integer, parameter       :: idebug = 3, nd = 4
      integer      :: ik, ii, i, j, iin, iout, lstrow, ncon, nadh, nslip, nplast, nexter
      logical      :: znewln, is_roll, is_ssrol, use_plast
      real(kind=8) :: tmp1mx, tmp2mx, hmp, trcbnd, ptabs, ptarg, rhsx, rhsy, fxtrue, fytrue, us_avg(3), &
                      frpow0
      character(len= 9) :: str_muscal
      character(len=12) :: strng(max(20, nsens_in*nsens_out))
      type(t_gridfnc3) :: dupl

      call timer_start(itimer_output)
      associate(ic     => gd%ic,            mater  => gd%mater,         pot    => gd%potcon_cur,        &
                mx     => gd%potcon_cur%mx, my     => gd%potcon_cur%my, npot   => gd%potcon_cur%npot,   &
                dx     => gd%potcon_cur%dx, dy     => gd%potcon_cur%dy, dxdy   => gd%potcon_cur%dxdy,   &
                x      => gd%cgrid_cur%x,   y      => gd%cgrid_cur%y,                                   &
                muscal => gd%kin%muscal,    pen    => gd%kin%pen,       cksi   => gd%kin%cksi,          &
                ceta   => gd%kin%ceta,      cphi   => gd%kin%cphi,      fntrue => gd%kin%fntrue,        &
                fxrel  => gd%kin%fxrel,     fyrel  => gd%kin%fyrel,     fcntc  => gd%kin%fcntc,         &
                chi    => gd%kin%chi,       dq     => gd%kin%dq,        veloc  => gd%kin%veloc,         &
                eps    => gd%solv%eps,                                                                  &
                hs1    => gd%geom%hs1,      exrhs  => gd%geom%exrhs,                                    &
                cs     => gd%influ%cs,      cv     => gd%influ%cv,                                      &
                mus1   => gd%outpt1%mus,    igs1   => gd%outpt1%igs,    igv1   => gd%outpt1%igv,        &
                ps1    => gd%outpt1%ps,     pv1    => gd%outpt1%pv,     us1    => gd%outpt1%us,         &
                uv1    => gd%outpt1%uv,     ss1    => gd%outpt1%ss,     shft1  => gd%outpt1%shft,       &
                taucs  => gd%outpt1%taucs,  upls   => gd%outpt1%upls,   uplv   => gd%outpt1%uplv,       &
                temp1  => gd%outpt1%temp1,  temp2  => gd%outpt1%temp2,                                  &
                mxtrue => gd%outpt1%mxtrue, mytrue => gd%outpt1%mytrue, mztrue => gd%outpt1%mztrue,     &
                elen   => gd%outpt1%elen,   frpow  => gd%outpt1%frpow,  pmax   => gd%outpt1%pmax)

      call gf3_new(dupl, 'output:dupl', gd%cgrid_cur, nulify=.true.)

      is_roll  = ic%tang.eq.2 .or. ic%tang.eq.3
      is_ssrol = ic%tang.eq.3
      use_plast = ic%mater.eq.4 .and. mater%tau_c0.gt.1d-10 .and. mater%tau_c0.le.1d10

      ! Print the heading if "output" >= 2

      if (ic%output_surf.ge.2 .and. out_open.eq.1) then

         write(lout, 8000)
         write(lout, 8001) ic%pvtime, ic%bound, ic%tang, ic%norm, ic%force3, ic%stress
         if (ic%varfrc.eq.0) then
            write(lout, 8002) ic%frclaw_inp, ic%discns3, ic%gencr_inp, ic%mater, ic%rznorm, ic%rztang
         else
            write(lout, 8003) ic%varfrc, ic%frclaw_inp, ic%discns3, ic%gencr_inp, ic%mater, ic%rznorm,  &
                        ic%rztang
         endif
         write(lout, 8004) ic%heat, ic%gausei_inp, ic%iestim, ic%matfil_surf, ic%output_surf,           &
                                ic%flow, ic%return
 8000    format (' ')
 8001    format ('  PBTNFS.        ',2x, '  PVTIME',i2, ', BOUND ',i2, ', TANG  ',i2, ', NORM  ',i2,    &
                        ', FORCE ',i2, ', STRESS',i2)
 8002    format ('  LDCMZE.        ',2x, '  FRCLAW',i2, ', DISCNS',i2, ', INFLCF',i2, ', MATER ',i2,    &
                        ', RZNORM',i2, ', RZTANG',i2)
 8003    format (' VLDCMZE.  VARFRC',i2, ', FRCLAW',i2, ', DISCNS',i2, ', INFLCF',i2, ', MATER ',i2,    &
                        ', RZNORM',i2, ', RZTANG',i2)
 8004    format (' HGIAOWR.  HEAT  ',i2, ', GAUSEI',i2, ', IESTIM',i2, ', MATFIL',i2, ', OUTPUT',i2,    &
                        ', FLOW  ',i2, ', RETURN',i2)

         ! z == (Norm called)

         if (ic%bound.eq.0 .and. ic%norm.eq.0) write(lout, 8011) 'APPROACH', 'FULL'
         if (ic%bound.eq.0 .and. ic%norm.eq.1) write(lout, 8011) 'NORMAL FORCE', 'FULL'
         if (ic%bound.eq.5 .and. ic%norm.eq.0) write(lout, 8011) 'APPROACH', 'APPROXIMATE KPEC'
         if (ic%bound.eq.5 .and. ic%norm.eq.1) write(lout, 8011) 'NORMAL FORCE', 'APPROXIMATE KPEC'
         if (ic%bound.eq.6 .and. ic%norm.eq.0) write(lout, 8011) 'APPROACH', 'APPROXIMATE ANALYN'
         if (ic%bound.eq.6 .and. ic%norm.eq.1) write(lout, 8011) 'NORMAL FORCE', 'APPROXIMATE ANALYN'
         if (ic%bound.eq.2) write(lout, 8013) 'ELLIPTICAL'
         if (ic%bound.eq.3) write(lout, 8013) 'PARABOLICAL'
         if (ic%bound.eq.4) write(lout, 8013) 'SDEC'
 8011    format (' NORMAL PROBLEM WITH PRESCRIBED ',a,', ',a,' SOLUTION')
 8013    format (' NORMAL PROBLEM SKIPPED, ',a,' TRACTION BOUND')

         if (ic%tang.ne.0) then
            if (ic%mater.ge.2 .and. ic%mater.le.3) then
               write(lout, 8020)
            elseif (ic%mater.eq.5) then
               write(lout, 8025)
            else
               if (ic%tang.eq.1) write(lout, 8021)
               if (ic%tang.eq.2) write(lout, 8022)
               if (ic%tang.eq.3) write(lout, 8023)
            endif
 8020       format (' STEADY STATE ROLLING,  "FASTSIM APPROACH"')
 8021       format (' SHIFT TRANSIENT')
 8022       format (' TRANSIENT ROLLING')
 8023       format (' STEADY STATE ROLLING')
 8025       format (' STEADY STATE ROLLING,  "FASTRIP APPROACH"')

            if (ic%force3.eq.0) write(lout, 8030)
            if (ic%force3.eq.1) write(lout, 8031)
            if (ic%force3.eq.2) write(lout, 8032)
 8030       format (' CREEPAGE PRESCRIBED')
 8031       format (' X-FORCE,  Y-CREEPAGE PRESCRIBED')
 8032       format (' TANGENTIAL FORCE PRESCRIBED')

            if (gd%solv%maxout.eq.1) write(lout, 8041)
            if (gd%solv%maxout.gt.1) write(lout, 8042) gd%solv%maxout
 8041       format (' JOHNSON"S METHOD IS USED FOR THE FRICTIONAL STRESS')
 8042       format (' PANAGIOTOPOULOS" METHOD IS USED FOR THE FRICTIONAL STRESS. MAXOUT=', I3)
         endif

         ! print data of Hertzian problem

         if (pot%ipotcn.lt.0) then

            if (pot%ipotcn.le.-6) then
               strng(1) = 'SDEC'
            elseif (pot%ipotcn.le.-4) then
               strng(1) = '2D HERTZIAN'
            else
               strng(1) = '3D HERTZIAN'
            endif

            if (ic%norm.eq.0) then
               strng(2) = 'APPROACH'
            else
               strng(2) = 'NORMAL FORCE'
            endif

            write(lout, 8051) trim(strng(1)), trim(strng(2))

 8051       format (/, 1x, a, ' GEOMETRY WITH PRESCRIBED ',a)

            if (pot%ipotcn.eq.-1) then
               write(lout, 8061) gd%hertz%a1, gd%hertz%b1
            elseif (pot%ipotcn.eq.-2) then
               write(lout, 8062) gd%hertz%a1, gd%hertz%aob
            elseif (pot%ipotcn.eq.-3) then
               write(lout, 8063) gd%hertz%aa, gd%hertz%bb
            elseif (pot%ipotcn.eq.-4) then
               write(lout, 8064) gd%hertz%a1, gd%hertz%bb
            elseif (pot%ipotcn.eq.-5) then
               write(lout, 8063) gd%hertz%aa, gd%hertz%bb
            elseif (pot%ipotcn.eq.-6) then
               write(lout, 8066) gd%hertz%aa, gd%hertz%bneg, gd%hertz%bpos
            endif

 8061       format (' CURVATURE PRESCRIBED, A1,B1: ', 13x, 2g12.4)
 8062       format (' CURVATURE AND ELLIPTIC. PRESCRB, A1,AA/BB:', 2g12.4)
 8063       format (' SEMIAXES PRESCRIBED, AA,BB:  ', 13x, 2g12.4)
 8064       format (' CURVATURE PRESCRIBED, A1,BB: ', 13x, 2g12.4)
 8066       format (' SEMIAXES PRESCRIBED, AA,BNEG,BPOS:', 8x, 3g12.4)

            if (pot%ipotcn.ne.-1 .and. pot%ipotcn.ne.-6) write(lout, 8071) gd%hertz%a1, gd%hertz%b1
            if (pot%ipotcn.ne.-3 .and. pot%ipotcn.ne.-6)                                                &
               write(lout, 8072) gd%hertz%aa, gd%hertz%bb, gd%hertz%aob
 8071       format (' THE CURVATURES A1,B1 ARE:    ', 13x, 2g12.4)
 8072       format (' BASIC SEMIAXES AA,BB, RATIO  ', 13x, 3g12.4)

            if (pot%ipotcn.ne.-6) write(lout, 8081) gd%hertz%rho, gd%hertz%cp
 8081       format (' EFFECTIVE RAD.CURV RHO, SEMI-AXIS CP', 6x, 2g12.4)

            write(lout, 8091)  pot%xl, pot%yl, gd%hertz%scale, mx, my, dx, dy
 8091       format (' POTENTIAL CONTACT, SCALE:    ', 13x, 3g12.4,/,                                    &
                    ' DISCRETISATION MX,MY, DX,DY: ', 13x, 2i6,2g12.4)

         endif
      endif

      ! print the input data if "output" >= 2.
      ! header per row: 2 spaces + [ 3 spaces + 9 characters ]*

      if (ic%output_surf.ge.2 .and. out_open.eq.1) then

         if (gd%solv%gausei_eff.eq.5) then
            if (gd%solv%gd_meth.le.2) then
               if (gd%solv%gd_meth.eq.1) gd%solv%kdown = max(999, mx+1)
               write(lout, 7902) gd%solv%kdown, gd%solv%kdowfb, gd%solv%betath, gd%solv%d_ifc,          &
                   gd%solv%d_lin, gd%solv%d_cns, gd%solv%d_slp, gd%solv%pow_s
            else
               write(lout, 7903) gd%solv%fdecay, gd%solv%kdowfb, gd%solv%betath, gd%solv%d_ifc,         &
                   gd%solv%d_lin, gd%solv%d_cns, gd%solv%d_slp, gd%solv%pow_s
            endif
 7902       format (/, ' PARAMETERS FOR SOLVER GDSTEADY', /,                                            &
               2x, 3x,'KDOWN',4x, 'KDOWFB  BETA', 3x,'D_IFC',4x, 3x,'D_LIN',4x, 3x,'D_CNS',4x,          &
                         3x,'D_SLP',4x, 3x,'POW_S',/, 2x, i8,4x, i5,f7.3, 6g12.4)
 7903       format (/, ' PARAMETERS FOR SOLVER GDSTEADY', /,                                            &
               2x, 3x,'FDECAY',3x, 'KDOWFB  BETA', 3x,'D_IFC',4x, 3x,'D_LIN',4x, 3x,'D_CNS',4x,         &
                         3x,'D_SLP',4x, 3x,'POW_S',/, 2x, g12.4, i5,f7.3, 5g12.4)
         endif

         write(lout, 8100) mater%nu, mater%ga, mater%ak, eps
 8100    format (/, ' MATERIAL CONSTANTS', /,                                                           &
            2x, 3x,'NU',7x, 3x,'G',8x, 3x,'AK',7x, 3x,'EPS',/, 2x, 4g12.4)

         write(lout, 8101) mater%poiss, mater%gg
 8101    format (/, 2x, 2x,'POISS(1)',2x, 2x,'POISS(2)',2x, 3x,'GG(1)',4x, 3x,'GG(2)',/, 2x, 4g12.4)

         if (ic%mater2.ge.2) then
            write(lout, 8102) fmt_gs(12,4,4,mater%cdampn), fmt_gs(12,4,4,mater%cdampt),                 &
                fmt_gs(12,4,4,mater%dfnmax), fmt_gs(12,4,4,mater%dftmax)
 8102       format (/, 2x, 2x,'CDAMPN',4x, 2x,'CDAMPT',4x, 3x,'DFNMAX',4x, 3x,'DFTMAX',/, 2x, 4a)
         endif

         if (ic%heat.ge.1) write(lout, 8108) mater%bktemp(1), mater%heatcp(1), mater%lambda(1),         &
                mater%dens(1), mater%bktemp(2), mater%heatcp(2), mater%lambda(2), mater%dens(2)
 8108    format (/, 2x, 3x,'BKTEMP',6x, 'HEATCP',6x, 'LAMBDA',6x, 'DENSITY', /, 2(2x, 4g12.4, :, /))

         if (ic%mater.eq.1) then
            strng(1) = fmt_gs(12, 4, 4, mater%tc(1))
            strng(2) = fmt_gs(12, 4, 4, mater%tc(2))
            strng(3) = fmt_gs(12, 4, 4, mater%vt(1))
            strng(4) = fmt_gs(12, 4, 4, mater%vt(2))
            write(lout,8111) mater%nuv, mater%gav, mater%akv, mater%fg, (strng(j),j=1,4)
 8111    format (/, 1x,'VISCOELASTIC MATERIAL CONSTANTS',/,                                             & 
            2x, 3x,'NUV',6x, 3x,'GV',7x, 3x,'AKV', /, 2x, 3g12.4,/,/,                                   &
            2x, 3x,'FG(1)',4x, 3x,'FG(2)',4x, 3x,'TC(1)',4x, 3x,'TC(2)',4x, 3x,'VT(1)',4x, 3x,'VT(2)',/, &
            2x, 2g12.4, 4a12)
         endif

         if (ic%mater.eq.2) write(lout,8121) mater%flx(1)
 8121    format (/, 1x,'SIMPLIFIED THEORY MATERIAL CONSTANTS',/, 2x, 3x,'FLX',/, 2x, g12.4)

         if (ic%mater.eq.3 .or. ic%mater.eq.5) write(lout,8131) mater%flx(1), mater%flx(2), mater%flx(3)
 8131    format (/, 1x,'SIMPLIFIED THEORY MATERIAL CONSTANTS',/,                                        &
                 2x, 3x,'FLX1',5x, 3x,'FLX2',5x, 3x,'FLX3',/, 2x, 3g12.4)

         if (ic%mater.eq.2 .or. ic%mater.eq.3) write(lout,8135) 'FASTSIM', mater%k0_mf, mater%alfamf,   &
                mater%betamf, mater%k_eff
         if (ic%mater.eq.5) write(lout,8135) 'FASTRIP', mater%k0_mf, mater%alfamf, mater%betamf, mater%k_eff
 8135    format (/, 1x,'MODIFIED ',a,' SLOPE REDUCTION PARAMETERS',/,                                   &
                 2x, 3x,'K0_MF',4x, 3x,'ALFAMF',3x, 3x,'BETAMF',3x, 3x,'K_EFF',/, 2x, 4g12.4)

         if (ic%mater.eq.4) write(lout,8141) mater%gg3, mater%laythk, mater%tau_c0, mater%k_tau
 8141    format (/, 1x,'INTERFACIAL LAYER PARAMETERS',/,                                                &
            2x, 3x,'GG3',6x, 3x,'LAYTHK',3x, 3x,'TAU_C0',3x, 3x,'K_TAU',/, 2x, 4g12.4)

         write(lout, *)

         call fric_output(lout, ic%output_surf, gd%fric)

         write(lout, *)

         if (gd%kin%use_muscal) then
            str_muscal = '/FSTAT/FN'
         else
            str_muscal = '/FN      '
         endif
         strng(1) = fmt_gs(12, 4, 4, veloc)
         strng(2) = fmt_gs(12, 4, 4, cksi)
         strng(3) = fmt_gs(12, 4, 4, ceta)
         strng(4) = fmt_gs(12, 4, 4, cphi)
         strng(5) = fmt_gs(12, 4, 4, fxrel)
         strng(6) = fmt_gs(12, 4, 4, fyrel)
         write(lout, 8301)
         if (.not.is_roll) then
            if (ic%force3.eq.0) write(lout, 8302)             gd%kin%dt, strng(1), strng(2), strng(3),  &
                strng(4)
            if (ic%force3.eq.1) write(lout, 8303) str_muscal, gd%kin%dt, strng(1), strng(5), strng(3),  &
                strng(4)
            if (ic%force3.eq.2) write(lout, 8304) str_muscal, str_muscal, gd%kin%dt, strng(1), strng(5), &
                strng(6), strng(4)
         else
            if (ic%force3.eq.0) write(lout, 8305) chi, dq, strng(1), strng(2), strng(3), strng(4)
            if (ic%force3.eq.1) write(lout, 8306) str_muscal, chi, dq, strng(1), strng(5), strng(3),    &
                strng(4)
            if (ic%force3.eq.2) write(lout, 8307) str_muscal, str_muscal, chi, dq, strng(1), strng(5),  &
                strng(6), strng(4)
         endif
 8301    format (' KINEMATIC CONSTANTS')
 8302    format (2x, 3x,'DT',7x,  3x,'VELOC',4x, 3x,'CKSI',5x,  3x,'CETA',5x, 3x,'CPHI',/,              &
                 2x, f8.3,4x, 4a12, /)
 8303    format (2x, 3x,'DT',7x,  3x,'VELOC',4x, 1x,'FX',a,     3x,'CETA',5x, 3x,'CPHI',/,              &
                 2x, f8.3,4x, 4a12, /)
 8304    format (2x, 3x,'DT',7x,  3x,'VELOC',4x, 1x,'FX',a,     1x,'FY',a,    3x,'CPHI',/,              &
                 2x, f8.3,4x, 4a12, /)
 8305    format (2x, 3x,'CHI',6x, 3x,'DQ',7x,    3x,'VELOC',4x, 3x,'CKSI',5x, 3x,'CETA',5x, 3x,'CPHI',/, &
                 2x, 2g12.4, 4a12, /)
 8306    format (2x, 3x,'CHI',6x, 3x,'DQ',7x,    3x,'VELOC',4x, 1x,'FX',a,    3x,'CETA',5x, 3x,'CPHI',/, &
                 2x, 2g12.4, 4a12, /)
 8307    format (2x, 3x,'CHI',6x, 3x,'DQ',7x,    3x,'VELOC',4x, 1x,'FX',a,    1x,'FY',a,    3x,'CPHI',/, &
                 2x, 2g12.4, 4a12, /)

         if (max(abs(gd%kin%spinxo),abs(gd%kin%spinyo)).ge.tiny) then
            write(lout, 8401) gd%kin%spinxo, gd%kin%spinyo
         endif
 8401    format (2x, 3x,'SPINXO',3x,  3x,'SPINYO',/, 2x, 2g12.4, /)
      endif

      ! here starts the calculating module of output.

      if (.true.) then

         ! Copy ps to pv in steady state problems, and copy corresponding el.div as well.
         ! TODO: fill all outputs, e.g. uv, sv, taucv, uplv, etc.
         ! TODO: move copying/shifting to computational routines, e.g. end of Tang or panprc.

         if (is_ssrol) then
            call gf3_copy(AllElm, ps1, pv1, ikALL)
            call eldiv_copy(igs1, igv1, ikALL)
         endif

         ! Compute displacement differences at current time, us(ii,ik).
         !  - set normal part of us to 0 when using B=2,3,4.
         !  - compute for AllElm when A>=2 for writing to mat-file
         !
         ! Note: module 1 uses A3 = -A1 (0, -1 or -2) in call to soutpt
         ! Note: ignoring reduced interaction between combined sub-contact patches

         call areas(igs1)
         if (ic%bound.ge.2 .and. ic%bound.le.4) then
            call gf3_set(AllElm, 0d0, us1, ikZDIR)
         elseif (abs(ic%matfil_surf).ge.2) then
            call VecAijPj(igs1, AllElm, us1, ikZDIR, ps1, jkALL, cs)
         else
            call VecAijPj(igs1, AllInt, us1, ikZDIR, ps1, jkALL, cs)
         endif

         !  - compute tangential part of us 
         !    Note: maintain the stored tangential part of us when M=2, 3 or 5 (Fastsim/FaStrip).
         !    TODO: exploit that tangential tractions are zero when T=0?

         if (ic%mater.ne.2 .and. ic%mater.ne.3 .and. ic%mater.ne.5) then
            if (abs(ic%matfil_surf).ge.2 .or. .false.) then
               call VecAijPj(igs1, AllElm, us1, ikTANG, ps1, jkALL, cs)
            else
               call VecAijPj(igs1, AllInt, us1, ikTANG, ps1, jkALL, cs)
            endif
         endif

         ! Compute displacement differences at previous time, uv(ii,ik).
         !  - needed for elastic energy & calculation of slip velocity
         !  - set normal part of uv to 0 when using B=2,3,4

         call areas(igv1)
         if (ic%bound.ge.2 .and. ic%bound.le.4) then
            call gf3_set(AllInt, 0d0, uv1, ikZDIR)
         else
            call VecAijPj(igs1, AllInt, uv1, ikZDIR, pv1, jkALL, cv)
         endif

         !  - maintain the stored tangential part of uv when M=2, 3 or 5 (Fastsim/FaStrip)
         !    TODO: exploit that tangential tractions are zero when T=0?

         if (ic%mater.ne.2 .and. ic%mater.ne.3 .and. ic%mater.ne.5) then
            call VecAijPj(igs1, AllInt, uv1, ikTANG, pv1, jkALL, cv)
         endif

         ! Count number of elements in contact, slip and adhesion

         call eldiv_count(igs1, nadh, nslip, nplast, nexter)
         ncon = nadh + nslip + nplast

         ! Compute total forces that are not prescribed

         if (ic%norm.eq.0) then
            fntrue = dxdy * gf3_sum(AllElm, ps1, ikZDIR)
         endif
     
         if (ic%force3.eq.0) then
            fxrel = dxdy * gf3_sum(AllElm, ps1, ikXDIR) / (fntrue*muscal + tiny)
         endif
         if (ic%force3.le.1) then
            fyrel = dxdy * gf3_sum(AllElm, ps1, ikYDIR) / (fntrue*muscal + tiny)
         endif
         fxtrue = fxrel * (fntrue*muscal+tiny)
         fytrue = fyrel * (fntrue*muscal+tiny)
         fcntc  = (/ fxtrue, fytrue, fntrue /)

         ! compute ad-hoc proportional damping

         call calc_damping_force( gd%ic, mater, gd%kin )

         ! Compute torsional moments about x, y and z-axes

         mxtrue =   dxdy * ddot(npot, ps1%vn, 1, y, 1)
         mytrue = - dxdy * ddot(npot, ps1%vn, 1, x, 1)
         mztrue = - dxdy * ddot(npot, ps1%vx, 1, y, 1)  + dxdy * ddot(npot, ps1%vy, 1, x, 1)
         if (abs(mztrue).lt.0.5d0*eps*(muscal*fntrue+tiny)) mztrue=0d0

         ! Compute Elastic Energy
         ! Note: set elen = 0 when using B=2,3,4.
         ! Note: exploiting Ps==0 in Exter, using AllInt instead of AllElm

         if (ic%bound.ge.2 .and. ic%bound.le.4) then
            elen  = 0d0
         else
            elen  = 0.5d0 * 1d-3 * dxdy * gf3_dot(Allint, us1, ps1, ikALL)
         endif

         ! Compute the deformed distance  hs - Pen + Us(.,1) and store in ss(.,1)
         ! Note: set Ss = 0 when using B=2,3,4.

         if (ic%bound.ge.2 .and. ic%bound.le.4) then
            call gf3_set(AllElm, 0d0, ss1, ikZDIR)
         else
            call gf3_set(AllElm, -pen, ss1, ikZDIR)
            call gf3_axpy(AllElm, 1d0, hs1, ss1, ikZDIR)
            call gf3_axpy(AllElm, 1d0, us1, ss1, ikZDIR)
         endif

         ! Re-allocate shft-array at appropriate size

         call gf3_new(shft1, 'output:shft1', gd%cgrid_cur)

         ! Compute the magnitude of the slip distance (shift),
         !    should only be <> 0 in the slip-area

         do ii = 1, npot
            if (igs1%el(ii).eq.Slip) then
               shft1%vt(ii) = dsqrt(ss1%vx(ii)**2 + ss1%vy(ii)**2)
            else
               ss1%vx(ii) = 0d0
               ss1%vy(ii) = 0d0
               shft1%vt(ii) = 0d0
            endif
         enddo

         ! Compute the frictional power [W] = [J/s] = [N.m/s]   (power == force times velocity)

         !    Fric = V \int_C  P_t * s_t  dxdy 

         ! Here s_t stands for the relative slip velocity of two opposing particles, multiplying by V makes
         ! this the absolute slip velocity. Array ss = S_t contains the shift distance [mm]. This is divided
         ! by 1000 dt to get the absolute velocity sa_t in [m/s].

         if (use_plast) then
            call gf3_copy(AllElm, upls, dupl, ikTANG)
            call gf3_axpy(AllElm, -1d0, uplv, dupl, ikTANG)
            frpow0 = -dxdy * gf3_dot(AllElm, ps1, dupl, ikTANG) / (1d3 * gd%kin%dt)
         else
            frpow0 = 0d0
         endif

         frpow  = frpow0 + dxdy * gf3_dot(AllElm, ps1, ss1, ikTANG) / (1d3 * gd%kin%dt)

         ! compute the maximum pressure

         pmax  = 0d0
         if (ic%print_pmax) pmax  = gf3_max(AllElm, ps1, ikZDIR)

         ! compute average displacement differences

         if (ic%print_uxavg) then
            do ik = 1, 3
               us_avg(ik) = gf3_sum(AllInt, us1, ik) / real(ncon)
            enddo
         endif

      endif ! true

      ! Here ends the calculating module of output.

      if (ic%output_surf.le.0 .or. out_open.ne.1) goto 90

      ! Arrived here, Output takes place. PRINT RESULTS. First total force
      ! and torsional moment, then elastic energy and frictional power.
      ! header per row: 2 spaces + [ 3 spaces + 9 characters ]*
      ! Filter values that are dominated by noise using filt_sml()

      strng(1) = fmt_gs(12,nd, 4, fntrue)
      strng(2) = fmt_gs(12,nd, 4, filt_sml(fxtrue,0.5d0*eps*fntrue*muscal))
      strng(3) = fmt_gs(12,nd, 4, filt_sml(fytrue,0.5d0*eps*fntrue*muscal))
      strng(4) = fmt_gs(12, 4, 4, mztrue)
      strng(5) = fmt_gs(12, 4, 4, elen)
      strng(6) = fmt_gs(12, 4, 4, frpow)

      if (ic%output_surf.ge.2) write(lout,8500)
      if (.not.is_roll) then
         write(lout,8501) 'WORK'
      else
         write(lout,8501) 'POWER'
      endif
      write(lout,8502) (strng(j), j=1,6)
 8500 format (' TOTAL FORCES, TORSIONAL MOMENT, ELASTIC ENERGY, FRICTIONAL POWER')
 8501 format (2x, 3x,'FN',7x, 3x,'FX',7x, 3x,'FY',7x, 3x,'MZ',7x, 2x,'ELAST.EN.',1x, 2x,'FRIC.',a)
 8502 format (2x, 6a12)

      strng(3) = fmt_gs(12, 4, 4, fntrue/mater%ga)
      if (ic%force3.eq.0) then
         strng(1) = 'FX/FSTAT/FN'
         if (.not.gd%kin%use_muscal) strng(1) = '  FX/FN'
         strng(4) = fmt_gs(12, 4, 4, filt_sml(fxrel, 0.5d0*eps))
      elseif (is_roll) then
         strng(1) = ' CREEP X'
         strng(4) = fmt_gs(12,nd, 4, filt_sml(cksi, ceta*eps))
      else
         strng(1) = ' SHIFT X'
         strng(4) = fmt_gs(12,nd, 4, filt_sml(cksi, ceta*eps))
      endif
      if (ic%force3.le.1) then
         strng(2) = 'FY/FSTAT/FN'
         if (.not.gd%kin%use_muscal) strng(2) = ' FY/FN'
         strng(5) = fmt_gs(12, 4, 4, filt_sml(fyrel, 0.5d0*eps))
      elseif (is_roll) then
         strng(2) = ' CREEP Y'
         strng(5) = fmt_gs(12, 4, 4, filt_sml(ceta, cksi*eps))
      else
         strng(2) = ' SHIFT Y'
         strng(5) = fmt_gs(12, 4, 4, filt_sml(ceta, cksi*eps))
      endif
      strng(6) = fmt_gs(12, 6, 4, pen)

      if (ic%print_pmax) then
         strng(7)  = '    PMAX    '
         strng(9)  = fmt_gs(12,nd, 4, pmax)
         ! write(strng(9),'(f12.8)') pmax
      endif
      if (ic%print_uxavg) then
         strng(8)  = '  UX_AVG    '
         strng(10) = fmt_gs(12,nd, 4, us_avg(1))
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

      if (ic%heat.ge.1 .or. ic%print_uxavg) then
         write(lout,8511) (strng(j), j=1,2), strng(7), strng(8)
         write(lout,8512) (strng(j), j=3,6), strng(9), strng(10)
      elseif (ic%print_pmax) then
         write(lout,8511) (strng(j), j=1,2), strng(7)
         write(lout,8512) (strng(j), j=3,6), strng(9)
      else
         write(lout,8511) (strng(j), j=1,2)
         write(lout,8512) (strng(j), j=3,6)
      endif
 8511 format('     FN/G      ',a12,1x,a12,' APPROACH',:,1x,2a12)
 8512 format(2x, 6a12)

      ! Write more detailed global output when "output" >= 2

      if (ic%output_surf.ge.2) then

         ! Convert sensitivities to strings, row-wise

         do iout = 1, nsens_out
            do iin = 1, nsens_in
               i = nsens_in*(iout-1) + iin
               strng(i) = fmt_gs(12, 3, 3, gd%outpt1%sens(iout,iin))
            enddo
         enddo

         if (ic%tang.ne.0) then

            ! sensitivities for Norm and Tang, one contact-problem

            str_muscal = 'FN/FSTAT.'
            if (.not.gd%kin%use_muscal) str_muscal = 'FN.'
            write(lout, 8621) str_muscal
 8621       format(/, ' THE SENSITIVITIES. (FX,FY) MEANS: (FX,FY)/',a,/,                                &
                      ' A ZERO ENTRY MEANS THAT IT HAS NOT BEEN CALCULATED.',/,                         &
                   8x,' DF/DPEN     DF/DKSI     DF/DETA     DF/DPHI')
            if (ic%sens.ge.2) then
               write(lout, 8622) ((strng(iin+(iout-1)*nsens_in), iin=iin_dpen,iin_dphi1),               &
                                                                              iout=iout_fn, iout_mz1)
            else
               write(lout, 8622) ((strng(iin+(iout-1)*nsens_in), iin=iin_dpen,iin_dphi1),               &
                                                                              iout=iout_fn, iout_fy1)
            endif
 8622       format (3x,'FN ',4a12,/, 3x,'FX ',4a12,/, 3x,'FY ',4a12,:,/, 3x,'MZ ',4a12)

         else

            ! sensitivities for Norm only

            write(lout, 8631) strng(1)
 8631       format (/, ' THE SENSITIVITIES.',/,                                                         &
                       ' A ZERO MEANS THAT IT HAS NOT BEEN CALCULATED.', /, &
                    2x,'DFN/DPEN', /, a12)

         endif

         ! write the Statistics

         if (use_plast) then
            write(lout, 8641) '  NPLAST', npot, ncon, nadh, nslip, nplast, gd%solv%itnorm, gd%solv%ittang
         else
            write(lout, 8641) ' ', npot, ncon, nadh, nslip, gd%solv%itnorm, gd%solv%ittang
         endif
 8641    format (/, ' CONTACT STATISTICS', /,                                                           &
                    ' N: NUMBER OF ELEMENTS IN REGION,  I: NUMBER OF ITERATIONS,', /,                   &
                    ' POT: POTENTIAL CONTACT AREA,  CON: CONTACT,  ADH: ADHESION', /,                   &
                    '    NPOT   NCON   NADH  NSLIP',a,' INORM  ITANG', /, 1x, 7i7)
 8642    format (   8x, 3i7 )

      endif  ! O-digit >= 2

      ! Show picture of contact/adhesion/slip areas when O-digit >= 3

      if (ic%output_surf.ge.3) then
         write(lout, 8650)
 8650    format(/, ' FORM OF THE CONTACT, ADHESION AND SLIP-AREA`S:')
         call wrigs (igs1, is_roll, chi)
      endif

      !------------------------------------------------------------------------------------------------------

      ! Skip detailed output for all elements when O-digit <= 4

      if (ic%output_surf.le.4) goto 90

      ! write the detailed solution,  when "output" >= 5

      if (ic%tang.eq.0) then
      write(lout, 8700)
 8700    format (/, ' QUANTITIES DISTRIBUTED OVER THE CONTACT AREA:', /,                                &
                    ' "H-PEN" IS THE UNDEFORMED DISTANCE WITH THE APPROACH TAKEN INTO ACCOUNT.',/,      &
                    ' PN IS THE NORMAL PRESSURE, PT IS THE TANGENTIAL TRACTION (2-VECTOR).',/,          &
                    ' X, H AND PEN ARE IN MM, PN AND PT ARE IN N/MM^2, AND ARGUMENTS ARE IN DEGREES.',//,&
                     5x,'X',10x,'H-PEN',8x, 'PN', 7x, 'ABS(PT)', 8x, 'ARG(PT)')

      elseif (.not.is_roll) then
      write(lout, 8701)
 8701    format (/, ' QUANTITIES DISTRIBUTED OVER THE CONTACT AREA:', /,                                &
                    ' "H-PEN" IS THE UNDEFORMED DISTANCE WITH THE APPROACH TAKEN INTO ACCOUNT.',/,      &
                    ' PN IS THE NORMAL PRESSURE, PT THE TANGENTIAL TRACTION (2-VECTOR).',/,             &
                    ' S IS THE SLIP DISTANCE PER STEP, WHICH IS OPPOSITE TO THE TANGENTIAL TRACTION.',/,&
                    ' RIGID SHIFT(X, Y) ARE THE (X, Y) COMPONENTS OF THE RIGID SHIFT DISTANCE W.',/,    &
                    ' X, H, PEN, S AND W ARE IN MM, PN AND PT ARE IN N/MM^2, AND ARGUMENTS ARE IN DEGREES.',//,&
                     /,        &
                     5x,'X',9x,'H-PEN',8x,'PN',6x,'TRCBND    ABS(PT)     ARG(PT;-S)',2x,'ABS(S)',6x,    &
                        'RIGID SHIFT(X, Y)')

      else
      write(lout, 8702)
         if (ic%heat.eq.0) then
            write(lout, 8703)
         else
            write(lout, 8704)
         endif

 8702    format (/, ' QUANTITIES DISTRIBUTED OVER THE CONTACT AREA:', /,                                &
                    ' "H-PEN" IS THE UNDEFORMED DISTANCE WITH THE APPROACH TAKEN INTO ACCOUNT.',/,      &
                    ' PN IS THE NORMAL PRESSURE, PT THE TANGENTIAL TRACTION (2-VECTOR).',/,             &
                    ' S IS THE RELATIVE SLIP VELOCITY, WHICH IS OPPOSITE TO THE TANGENTIAL TRACTION.',/ &
                    ' RIGID SLIP(X, Y) ARE THE (X, Y) COMPONENTS OF THE RELATIVE RIGID SLIP VELOCITY W.',/&
                    ' X, H AND PEN ARE IN MM, PN AND PT ARE IN N/MM^2, S AND W ARE DIMENSIONLESS, ',    &
                                                                    'ARGUMENT IS IN DEGREES.',/)
 8703    format (   5x,'X',9x,'H-PEN',8x,'PN',6x,'TRCBND    ABS(PT)     ARG(PT;-S)',2x,'ABS(S)',6x,     &
                                                                                'RIGID SLIP(X, Y)')
 8704    format (   5x,'X',9x,'H-PEN',8x,'PN',6x,'TRCBND    ABS(PT)     ARG(PT;-S)',2x,'ABS(S)',6x,     &
                                                      'RIGID SLIP(X, Y)',8x,'TEMP1',6x,'TEMP2')
      endif

      zNewLn = .true.
      LstRow = gd%cgrid_cur%iy(1)

      ! for all elements inside the contact area do

      do 871 ii = 1, npot
         if (gd%cgrid_cur%iy(ii).ne.lstrow) znewln = .true.
         lstrow = gd%cgrid_cur%iy(ii)
         if (igs1%el(ii).ge.Adhes .or. ic%output_surf.ge.6) then

            ! print header when starting a new row of the potential contact:

            if (znewln) then
               write(lout, 8710) y(ii), gd%cgrid_cur%iy(ii)
               znewln = .false.
 8710          format (/, ' Y = ', g12.4, ' ROW', i3, ' OF THE POTENTIAL CONTACT', /)
            endif

            ! retrieve/determine output quantities

            ptabs = dsqrt (ps1%vx(ii) **2 + ps1%vy(ii) **2)
            ptarg = atan2(ps1%vy(ii), ps1%vx(ii)) * 180d0/pi
            hmp = hs1%vn(ii) - pen
            if (ic%tang.ne.0) then

               ! hs, ss  == (rigid) slip distance per (time-)step
               !  - rolling: rhs, creepage == relative slip velocity
               !  - shifts:  rhs, creepage == slip distance, dq == 1

               if (ic%force3.ge.1) then
                  rhsx = filt_sml(cksi, ceta*eps) - hs1%vx(ii) / dq
               else
                  rhsx =                          - hs1%vx(ii) / dq
               endif
               if (ic%force3.eq.2) then
                  rhsy = filt_sml(ceta, cksi*eps) - hs1%vy(ii) / dq
               else
                  rhsy =                          - hs1%vy(ii) / dq
               endif
               trcbnd = ps1%vn(ii) * mus1%vt(ii)
               if (use_plast) trcbnd = min(trcbnd, taucs%vt(ii))

               if (ic%heat.eq.0) then
                  write(lout, 8720) x(ii), hmp, ps1%vn(ii), trcbnd, ptabs, ptarg, shft1%vt(ii)/dq,      &
                         rhsx, rhsy
               else
                  write(lout, 8721) x(ii), hmp, ps1%vn(ii), trcbnd, ptabs, ptarg, shft1%vt(ii)/dq,      &
                        rhsx, rhsy, temp1%vn(ii), temp2%vn(ii)
               endif
                  
 8720          format (1x, g12.4, 4g11.3, f10.1, 2x, 3g11.3)
 8721          format (1x, g12.4, 4g11.3, f10.1, 2x, 3g11.3, 2g11.4)
            else
               write(lout, 8730) x(ii), hmp, ps1%vn(ii), ptabs, ptarg
 8730          format (1x, 4g12.4, f10.1)
            endif
         endif

      ! end for (all elements inside contact area)

 871  continue

 90   continue

      ! copy penetration for new time step

      gd%kin%penv = pen

      call gf3_destroy(dupl)
      call timer_stop(itimer_output)

      end associate

   end subroutine soutpt

!------------------------------------------------------------------------------------------------------------

subroutine writmt (meta, ic, cgrid, potcon, mater, fric, kin, geom, outpt1, mirror_y)
!--purpose: writes variables into the output-file <EXPERIM>.<CASE>.MAT. This is a file that is used for
!          communication with the plot-scripts written in MatLab. See User Guide (A.6) for a specification.
   implicit none
!--subroutine parameters:
   type(t_metadata)            :: meta
   type(t_ic)                  :: ic
   type(t_grid)                :: cgrid
   type(t_potcon)              :: potcon
   type(t_material)            :: mater
   type(t_friclaw)             :: fric
   type(t_kincns)              :: kin
   type(t_geomet)              :: geom
   type(t_output)              :: outpt1
   logical,         intent(in) :: mirror_y
!--local variables:
   integer, parameter      :: maxcol = 20
   logical             :: lwrall, is_roll, use_plast, wxy_insteadof_uxy, pv_insteadof_uxy
   integer             :: lmat, ncol, itauc, iuplx, iuply, itemp1, itemp2, j, ii, ix, iy, jj, jy
   real(kind=8)        :: sgn, xl, yl, values(maxcol)
   character(len=11)   :: colnam(maxcol)
   character(len=256)  :: fname

   is_roll    = ic%tang.eq.2 .or. ic%tang.eq.3
   use_plast  = ic%mater.eq.4 .and. mater%tau_c0.gt.1d-10 .and. mater%tau_c0.le.1d10
   wxy_insteadof_uxy = .false.
   pv_insteadof_uxy  = .false.

   sgn = 1d0
   if (mirror_y) sgn = -1d0

   associate(igs  => outpt1%igs, mus  => outpt1%mus, us   => outpt1%us, ps   => outpt1%ps,              &
             pv   => outpt1%pv,  shft => outpt1%shft)

   ! check if the file unit-number is ok

   lmat = get_lunit_tmp_use()
   if (lmat.le.0) goto 996

   ! write data for whole potential contact area when requested:

   lwrall = (ic%matfil_surf.ge.2)

   ! Save mat-files using names <EXPERIM>.0001.mat, ...

   write(fname, '(a,a,i0.4)') trim(meta%expnam), '.', meta%ncase

   ! Append 'a', 'b'... to case number in case of multiple patches

   if (meta%ipatch.ge.1) then
      fname = trim(fname) // char(ichar('a')+meta%ipatch-1) // '.mat'
   else
      fname = trim(fname) // '.mat'
   endif

   call make_absolute_path(fname, meta%outdir, fname)

   ! Determine number of columns used

   ncol = 11
   if (use_plast) then
      itauc  = ncol + 1
      iuplx  = ncol + 2
      iuply  = ncol + 3
      ncol   = ncol + 3
   endif
   if (ic%heat.ge.1) then
      itemp1 = ncol + 1
      itemp2 = ncol + 2
      ncol   = ncol + 2
   endif

   if (ncol.gt.maxcol) then
      write(bufout,'(2(a,i3),a)') ' INTERNAL ERROR: ncol=',ncol,' > max=', maxcol,', aborting.'
      call write_log(1, bufout)
      call abort_run()
   endif

   ! Open file for transfer:

   open (unit=lmat, file=fname, action='readwrite', err=997)

      ! All columns are 14 positions wide. Two leading integers combined into single column.
      ! Descriptions start at 4th position in column.

      ! 1: write comment and meta-data. Last column: fmtmat, version number of mat-file

      write(lmat,101, err=998) ('  -', j=8,ncol-1), 'FMT'
      write(lmat,102) meta%tim, meta%s_ws, meta%th_ws, sgn*meta%y_r, meta%z_r, sgn*meta%roll_r,         &
                meta%rnom_rol, (0, j=8,ncol-1), fmtmat
 101  format('%  TIM',11x, 'S_WS(FC)',6x, 'TH_WS(FC)',5x, 'Y_R(TR)',7x, 'Z_R(TR)',7x, 'ROLL_R(TR)',4x,  &
                'RNOM_ROL',6x, 10(a3,:,11x))
 102  format(7g14.6, 20(i6,:,8x))

      write(lmat,103, err=998) ('  -', j=8,ncol)
      write(lmat,104) meta%x_w, sgn*meta%y_w, meta%z_w, sgn*meta%roll_w, sgn*meta%yaw_w, meta%rnom_whl, &
                (0d0, j=7,ncol)
 103  format('%  X_W(TR)',7x, 'Y_W(TR)',7x, 'Z_W(TR)',7x, 'ROLL_W(TR)',4x, 'YAW_W(TR)',5x,              &
                'RNOM_WHL',6x, 10(a3,:,11x))
 104  format(6g14.6, 20g14.6)

      ! 2: write comment and contact position

      write(lmat,105, err=998) ('-', j=10,ncol)
      write(lmat,106) meta%xcp_r, sgn*meta%ycp_r, meta%zcp_r, sgn*meta%deltcp_r, meta%xcp_w,            &
                sgn*meta%ycp_w, meta%zcp_w, meta%npatch, meta%ipatch, kin%spinxo, kin%spinyo,           &
                (0., j=12,ncol)
 105  format('%  X_CP(R)',7x, 'Y_CP(R)',7x, 'Z_CP(R)',7x, 'DELT_CP(R)',4x, 'X_CP(W)',7x, 'Y_CP(W)',7x,  &
                'Z_CP(W)',7x, 'NPATCH',8x, 'IPATCH',8x, 'SPINXO',8x, 'SPINYO',8x, 20(a3,:,11x))
 106  format(7g14.6, 2(i6,8x), 20g14.6)

      ! 3: write comment and parameters for the grid discretisation used

      xl = potcon%xl
      yl = potcon%yl
      if (mirror_y) yl = -(yl + cgrid%ny * cgrid%dy)
      write(lmat,121, err=998) ('-', j=12,ncol)
      write(lmat,122) cgrid%nx, cgrid%ny, xl, yl, cgrid%dx, cgrid%dy, kin%chi, kin%dq, ic%config,       &
                ic%discns1_eff, meta%ynom_whl, (0., j=12,ncol)
 121  format('%   MX',4x, 'MY',5x, 'XL',12x, 'YL',12x, 'DX',12x, 'DY',12x, 'CHI',11x, 'DQ',12x,          &
                'C1',12x, 'D',13x, 'YNOM_WHL',5x, 20(a3,:,11x))
 122  format(i7,1x,i4,2x, 6g14.6, 2(i6,8x), 20g14.6)

      ! 4: write comment and material parameters

      if (ic%mater.eq.0) then

         write(lmat,161) ('-', j=7,ncol)
         write(lmat,162) ic%tang, ic%mater, (mater%gg(j), j=1,2), (mater%poiss(j), j=1,2), (0., j=7,ncol)
 161     format('%    T',4x, 'M',6x, 'GG1',11x, 'GG2',11x, 'POISS1',8x, 'POISS2',8x, 20(a3,:,11x))
 162     format(i7,1x,i4,2x, 20g14.6)

      elseif (ic%mater.eq.1) then

         write(lmat,171) ('-', j=11,ncol)
         write(lmat,172) ic%tang, ic%mater, (mater%gg(j), j=1,2), (mater%poiss(j), j=1,2),              &
               (mater%fg(j), j=1,2), (mater%tc(j), j=1,2), (0., j=11,ncol)
 171     format('%    T',4x, 'M',6x, 'GG1',11x, 'GG2',11x, 'POISS1',8x, 'POISS2',8x, 'FG1',11x,         &
                'FG2',11x, 'TC1',11x, 'TC2',11x, 20(a3,:,11x))
 172     format(i7,1x,i4,2x, 20g14.6)

      elseif (ic%mater.eq.2 .or. ic%mater.eq.3 .or. ic%mater.eq.5) then

         write(lmat,181) ('-', j=12,ncol)               ! skipping poiss(2) for lack of columns
         write(lmat,182) ic%tang, ic%mater, (mater%gg(j), j=1,2), (mater%poiss(j), j=1,1),              &
               (mater%flx(j), j=1,3), mater%k0_mf, mater%alfamf, mater%betamf, (0., j=12,ncol)
 181     format('%    T',4x, 'M',6x, 'GG1',11x, 'GG2',11x, 'POISS1',8x, 'FLX1',10x, 'FLX2',10x,         &
                'FLX3',10x, 'K0_MF',9x, 'ALFAMF',8x, 'BETAMF',8x, 20(a3,:,11x))
 182     format(i7,1x,i4,2x, 20g14.6)

      elseif (ic%mater.eq.4) then

         write(lmat,191) ('-', j=11,ncol)
         write(lmat,192) ic%tang, ic%mater, (mater%gg(j), j=1,2), (mater%poiss(j), j=1,2), mater%gg3,   &
               mater%laythk, mater%tau_c0, mater%k_tau, (0., j=11,ncol)
 191     format('%    T',4x, 'M',6x, 'GG1',11x, 'GG2',11x, 'POISS1',8x, 'POISS2',8x, 'GG3',11x,         &
                'LAYTHK',8x, 'TAU_C0',8x, 'K_TAU',9x, 20(a3,:,11x))
 192     format(i7,1x,i4,2x, 20g14.6)

      endif

      ! 5: write comment and parameters for the friction law used

      if (fric%frclaw_eff.eq.0) then

         write(lmat,201) ('-', j=6,ncol)
         write(lmat,202) fric%frclaw_eff, ic%heat, kin%veloc, fric%fstat_arr(1), fric%fkin_arr(1),      &
                (0., j=6,ncol)
 201     format('%    L',4x, 'H',6x, 'VELOC',9x, 'FSTAT',9x, 'FKIN',10x, 20(a3,:,11x))
 202     format(i7,1x,i4,2x, 20g14.6)

      elseif (fric%frclaw_eff.eq.2) then

         write(lmat,221) ('-', j=11,ncol)
         write(lmat,222) fric%frclaw_eff, ic%heat, kin%veloc, fric%fkin_arr(1), fric%flin1(1),          &
               fric%sabsh1(1), fric%flin2(1), fric%sabsh2(1), fric%memdst, fric%mem_s0, (0., j=11,ncol)
 221     format('%    L',4x, 'H',6x, 'VELOC',9x, 'FKIN',10x, 'FLIN1',9x, 'SABSH1',8x, 'FLIN2',9x,       &
                'SABSH2',8x, 'MEMDST',8x, 'MEM_S0',8x, 20(a3,:,11x))
 222     format(i7,1x,i4,2x, 20g14.6)

      elseif (fric%frclaw_eff.eq.3) then

         write(lmat,231) ('-', j=11,ncol)
         write(lmat,232) fric%frclaw_eff, ic%heat, kin%veloc, fric%fkin_arr(1), fric%frat1(1),          &
               fric%sabsh1(1), fric%frat2(1), fric%sabsh2(1), fric%memdst, fric%mem_s0, (0., j=11,ncol)
 231     format('%    L',4x, 'H',6x, 'VELOC',9x, 'FKIN',10x, 'FRAT1',9x, 'SABSH1',8x, 'FRAT2',9x,       &
                'SABSH2',8x, 'MEMDST',8x, 'MEM_S0',8x, 20(a3,:,11x))
 232     format(i7,1x,i4,2x, 20g14.6)

      elseif (fric%frclaw_eff.eq.4) then

         write(lmat,241) ('-', j=11,ncol)
         write(lmat,242) fric%frclaw_eff, ic%heat, kin%veloc, fric%fkin_arr(1), fric%fexp1(1),          &
               fric%sabsh1(1), fric%fexp2(1), fric%sabsh2(1), fric%memdst, fric%mem_s0, (0., j=11,ncol)
 241     format('%    L',4x, 'H',6x, 'VELOC',9x, 'FKIN',10x, 'FEXP1',9x, 'SABSH1',8x, 'FEXP2',9x,       &
                'SABSH2',8x, 'MEMDST',8x, 'MEM_S0',8x, 20(a3,:,11x))
 242     format(i7,1x,i4,2x, 20g14.6)

      elseif (fric%frclaw_eff.eq.6) then

         write(lmat,261) ('-', j=11,ncol)
         write(lmat,262) fric%frclaw_eff, ic%heat, kin%veloc, fric%fref(1), fric%tref(1), fric%dfheat(1), &
               fric%dtheat(1), fric%memdst, fric%mem_s0, (0., j=10,ncol)
 261     format('%    L',4x, 'H',6x, 'VELOC',9x, 'FREF',10x, 'TREF',10x, 'DFHEAT',8x, 'DTHEAT',8x,      &
                'MEMDST',8x, 'MEM_S0',8x, 20(a3,:,11x))
 262     format(i7,1x,i4,2x, 20g14.6)

      endif

      ! Write comment and actual output values, 11+ columns, for interior elements in the contact area only

      colnam( 3) = 'H'
      colnam( 4) = 'MU'
      colnam( 5) = 'PN'
      colnam( 6) = 'PX'
      colnam( 7) = 'PY'
      colnam( 8) = 'UN'
      if (wxy_insteadof_uxy) then
         colnam( 9) = 'WX'
         colnam(10) = 'WY'
      elseif (pv_insteadof_uxy) then
         colnam( 9) = 'PVX'
         colnam(10) = 'PVY'
      else
         colnam( 9) = 'UX'
         colnam(10) = 'UY'
      endif
      if (is_roll) then
         colnam(11) = 'SREL'
      else
         colnam(11) = 'SHFT'
      endif
      if (use_plast) then
         colnam(itauc) = 'TAUCRT'
         colnam(iuplx) = 'UPLSX'
         colnam(iuply) = 'UPLSY'
      endif
      if (ic%heat.ge.1) then
         colnam(itemp1) = 'TEMP1'
         colnam(itemp2) = 'TEMP2'
      endif

      write(lmat,301) (colnam(j), j=3,ncol)

 301  format('%    I    IGS ', 20(3x,a,:))
 302  format(i7,1x, i4,2x, 20g14.6)

      do jy = 1, cgrid%ny
         do ix = 1, cgrid%nx
            if (.not.mirror_y) then
               iy = jy
            else
               iy = cgrid%ny + 1 - jy
            endif
            ii = ix + cgrid%nx * (iy-1) ! internal numbering ii == (ix,iy)
            jj = ix + cgrid%nx * (jy-1) ! mat-file numbering jj == (ix,jy)

            if (igs%el(ii).ge.Adhes .or. lwrall) then

               ! get values for element ii

               values( 3) =       geom%hs1%vn(ii)
               values( 4) =       mus%vt(ii)
               values( 5) =       ps%vn(ii)
               values( 6) =       ps%vx(ii)
               values( 7) = sgn * ps%vy(ii)
               values( 8) =       us%vn(ii)
               if (wxy_insteadof_uxy) then
                  values( 9) =       geom%hs1%vx(ii)/kin%dt
                  values(10) = sgn * geom%hs1%vy(ii)/kin%dt
               elseif (pv_insteadof_uxy) then
                  values( 9) =       pv%vx(ii)
                  values(10) = sgn * pv%vy(ii)
               else
                  values( 9) =       us%vx(ii)
                  values(10) = sgn * us%vy(ii)
               endif
               values(11) = shft%vt(ii)/kin%dq
               if (use_plast) then
                  values(itauc) =       outpt1%taucs%vt(ii)
                  values(iuplx) =       outpt1%upls%vx(ii)
                  values(iuply) = sgn * outpt1%upls%vy(ii)
               endif
               if (ic%heat.ge.1) then
                  values(itemp1) = outpt1%temp1%vn(ii)
                  values(itemp2) = outpt1%temp2%vn(ii)
               endif

               ! write values for element ii

               write(lmat,302) jj, igs%el(ii), (values(j), j=3,ncol)

            endif
         enddo ! ix
      enddo ! jy

      close(lmat)
      call free_lunit_tmp_use(lmat)
      return

      ! Error handling:

 996  continue
         call write_log(' ERROR: no unit-number provided for .mat-file; skipping mat-output.')
         goto 999
 997  continue
         call write_log(' ERROR: cannot open file "'//trim(fname)//'" for writing; skipping mat-output.')
         call free_lunit_tmp_use(lmat)
         goto 999
 998  continue
         call write_log(' ERROR: cannot write solution to .mat-file; continuing without write.')
         close(lmat)
         call free_lunit_tmp_use(lmat)
         goto 999
 999  continue

      end associate

end subroutine writmt

!------------------------------------------------------------------------------------------------------------

end module m_soutpt
