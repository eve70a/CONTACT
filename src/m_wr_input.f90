!------------------------------------------------------------------------------------------------------------
! m_wr_input - read input-file for one case with w/r contact (module 1)
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------

module m_wr_input

use m_hierarch_data
use m_wrprof_data
use m_subsurf       , only : rdsubs, wrsubs
implicit none
private

public  wr_input
public  wr_write_inp

contains

!------------------------------------------------------------------------------------------------------------

   subroutine wr_input (inp, ncase, linenr, wtd)
!--purpose: Input-routine for a W/R contact case. 
      implicit none
!--subroutine arguments:
      integer                   :: inp, linenr, ncase
      type(t_ws_track)          :: wtd
!--local variables:
      integer,      parameter :: mxnval = 20
      logical,      parameter :: lstop  = .true.
      integer          :: ints(mxnval), line0, nval, is_wheel, i_ftype, ldebug, ieof, ierror,           &
                          cpbtnfs, vldcmze, xhgiaowr, psflrin
      logical          :: flags(mxnval), is_roll, is_ssrol, was_roll, zerror
      real(kind=8)     :: dbles(mxnval), dx_prv, ds_prv
      character*16     :: types, namside
      character*256    :: strngs(mxnval)
      type(t_ic)       :: icold

      associate( ic    => wtd%ic,    meta  => wtd%meta, mater => wtd%mater, discr => wtd%discr,         &
                 kin   => wtd%kin,   fric  => wtd%fric, solv  => wtd%solv,  subs  => wtd%subs,          &
                 ws    => wtd%ws,    my_wheel => wtd%ws%whl,                                            &
                 trk   => wtd%trk,   my_rail  => wtd%trk%rai)

      ldebug = 1
      ieof   = -1 ! eof=error

      if (ldebug.ge.2) call write_log('--- Start subroutine wr_input ---')
      line0 = linenr

      if (inp.le.1) then
         call write_log('ERROR: screen input not supported in MODULE 1')
         call abort_run()
      endif

      ! save old ic for checks:

      icold = ic

      ! read control digits

      call readLine(linp, ncase, linenr, 'control integers cpbtnfs', 'i', ints, dbles, flags,           &
                strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
      cpbtnfs = ints(1)

      call readline(linp, ncase, linenr, 'control integers vldcmze', 'iI', ints, dbles, flags,          &
                strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
      vldcmze = ints(1)
      if (nval.ge.2) ic%gapwgt = ints(2)

      call readline(linp, ncase, linenr, 'control integers hgiaowr', 'i', ints, dbles, flags,           &
                strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
      xhgiaowr = ints(1)

      call ic_unpack (1, cpbtnfs, vldcmze, xhgiaowr, ic)

      ! set level of debug-output of input-routines. 0 = errors, 1 = warnings/info, >=2 = flow/debug

      if (ic%x_readln.le.0) then
         ldebug = ic%ilvout
      else
         ldebug = ic%x_readln
      endif

      ! Determine whether T describes rolling or not

      is_roll  = ic%tang.eq.2 .or. ic%tang.eq.3
      is_ssrol = ic%tang.eq.3
      was_roll = icold%tang.eq.2 .or. icold%tang.eq.3

      ! TEST  the input quantities (control integers)

      zerror = .false.

      zerror = zerror .or. .not.check_2rng ('Control digit C1', ic%config, 0, 1, 4, 5)
      zerror = zerror .or. .not.check_irng ('Control digit P',  ic%pvtime, 0, 2)
      zerror = zerror .or. .not.check_2rng ('Control digit B',  ic%bound , 0, 1, 5, 6)
      zerror = zerror .or. .not.check_irng ('Control digit T',  ic%tang  , 0, 3)
      zerror = zerror .or. .not.check_irng ('Control digit N1', ic%norm  , 0, 1)
      zerror = zerror .or. .not.check_2rng ('Control digit F1', ic%force1, 0, 1, 3, 3)
      zerror = zerror .or. .not.check_irng ('Control digit S',  ic%stress, 0, 3)

      zerror = zerror .or. .not.check_irng ('Control digit V',  ic%varfrc, 0, 1)
      zerror = zerror .or. .not.check_2rng ('Control digit L',  ic%frclaw_inp, 0, 4, 6, 6)
      zerror = zerror .or. .not.check_irng ('Control digit D',  ic%discns1_inp, 0, 9)
      zerror = zerror .or. .not.check_irng ('gap weighting',    ic%gapwgt, 0, 2)
      zerror = zerror .or. .not.check_irng ('Control digit C3', ic%gencr_inp , 0, 4)
      zerror = zerror .or. .not.check_irng ('Control digit M',  ic%mater , 0, 4)
      zerror = zerror .or. .not.check_irng ('Control digit Z1', ic%ztrack, 0, 3)
      zerror = zerror .or. .not.check_irng ('Control digit E1', ic%ewheel, 0, 5)

      zerror = zerror .or. .not.check_irng ('Control digit X',  ic%xflow , 0, 1)
      zerror = zerror .or. .not.check_2rng ('Control digit H',  ic%heat  , 0, 1, 3, 3)
      zerror = zerror .or. .not.check_irng ('Control digit G',  ic%gausei_inp, 0, 5)
      zerror = zerror .or. .not.check_irng ('Control digit I',  ic%iestim, 0, 3)
      zerror = zerror .or. .not.check_irng ('Control digit A',  ic%matfil_surf, 0, 2)
      zerror = zerror .or. .not.check_irng ('Control digit O',  ic%output_surf, 0, 5)
      zerror = zerror .or. .not.check_irng ('Control digit W',  ic%flow  , 0, 9)
      zerror = zerror .or. .not.check_irng ('Control digit R',  ic%return, 0, 3)

      ! Check requirements on profiles

      if (ic%config.eq.0 .or. ic%config.eq.4) then
         namside   = 'left'
      else
         namside   = 'right'
      endif

      if (my_rail%prr%fname.eq.' ' .and. (ic%ztrack.eq.0 .or. ic%ztrack.eq.2)) then
         zerror = .true.
         write(lout, 1901) trim(namside), ic%ztrack
         write(   *, 1901) trim(namside), ic%ztrack
 1901    format (' Input: ERROR. No profile has been given for the ', a, ' rail, Z1 =',i2,'.')
      endif

      if (my_wheel%prw%fname.eq.' ' .and.                                                               &
          (ic%ewheel.eq.1 .or. ic%ewheel.eq.2 .or. ic%ewheel.eq.4)) then
         zerror = .true.
         write(lout, 1902) trim(namside), ic%ewheel
         write(   *, 1902) trim(namside), ic%ewheel
 1902    format (' Input: ERROR. No profile has been given for the ', a, ' wheel, E1 =',i2,'.')
      endif

      ! Check requirements for first case

      if (ncase.eq.1 .and. ic%pvtime.le.1) then
         write(lout, 2001) ic%pvtime
         write(   *, 2001) ic%pvtime
 2001    format (' Input: WARNING. In the first case contact must be initiated.', /, &
                 16x,' Digit P=', i3,' is overruled, set to 2.')
         ic%pvtime = 2
      endif

      if (ncase.eq.1 .and. ic%frclaw_inp.eq.1) then
         zerror = .true.
         write(lout, 2002) ic%frclaw_inp
         write(   *, 2002) ic%frclaw_inp
 2002    format (' Input: ERROR. In the first case L must be.ne.1',/, ' Digit L=', i3)
      endif

      if (ncase.eq.1 .and. ic%discns1_inp.eq.0) then
         zerror = .true.
         write(lout, 2003) ic%discns1_inp
         write(   *, 2003) ic%discns1_inp
 2003    format (' Input: ERROR. In the first case D must be.ne.0',/, ' Digit D=', i3)
      endif

      if (ncase.eq.1 .and. ic%iestim.ne.0) then
         write(lout, 2004) ic%iestim
         write(   *, 2004) ic%iestim
 2004    format (' Input: WARNING. In the first case contact must be initiated.',/,     &
                 16x,' Digit I=', i3,' is overruled, set to 0.')
         ic%iestim = 0
      endif

      if (ncase.eq.1 .and. (ic%stress.eq.1 .or. ic%stress.eq.2)) then
         zerror = .true.
         write(lout, 2005) ic%stress
         write(   *, 2005) ic%stress
 2005    format (' Input: ERROR. In the first case S cannot be 1 or 2.'/,                              &
                 '               Digit S=', i3)
      elseif ((ic%stress.eq.1 .or. ic%stress.eq.2) .and. subs%nblock.le.0) then
         zerror = .true.
         write(lout, 2006) ic%stress, subs%nblock
         write(   *, 2006) ic%stress, subs%nblock
 2006    format (' Input: ERROR. No blocks are given for subsurface stresses; S cannot be 1 or 2.'/,   &
                 '               Digit S=', i3,', Nblock=',i3)
      endif

      ! Check consistency between consecutive cases

      if (ic%frclaw_inp.eq.1 .and. ic%varfrc.ne.icold%varfrc) then
         zerror = .true.
         write(lout, 2007) icold%varfrc, ic%varfrc, ic%frclaw_inp
         write(   *, 2007) icold%varfrc, ic%varfrc, ic%frclaw_inp
 2007    format (' Input: ERROR. When activating/deactivating variable friction (V=',i1,' ->',i2,')',/   &
                 '               new friction parameters are needed (L=',i1,').')
      endif

      if (ic%gencr_inp.eq.0) then
         ic%gencr_inp = 1
         write(lout, 2010) ic%gencr_inp
         write(   *, 2010) ic%gencr_inp
 2010    format (' Input: WARNING. Module 1 requires new influence functions, setting C3 :=', i3)
      endif

      if (ic%discns1_inp.ne.0 .and. ic%gencr_inp.eq.0) then
         zerror = .true.
         write(lout, 2011) ic%discns1_inp, ic%gencr_inp
         write(   *, 2011) ic%discns1_inp, ic%gencr_inp
 2011    format (' Input: ERROR. A new discretisation needs new influence functions'/, &
                ' Digits D, C3=', 2i3)
      endif

      if (is_ssrol .and. ic%pvtime.ne.2) then
         ic%pvtime = 2
         write(lout, 2021) ic%pvtime
         write(   *, 2021) ic%pvtime
 2021    format (' Input: WARNING. In steady state rolling previous Pv is ignored: digit P :=', i3)
      endif

      if (ic%ewheel.le.0 .and. (ic%norm.ne.icold%norm)) then
         zerror = .true.
         write(lout, 2031) ic%norm, icold%norm, ic%ewheel
         write(   *, 2031) ic%norm, icold%norm, ic%ewheel
 2031    format (' Input: ERROR. Switching between position and total force, new wheelset', /, &
                 '               position input is needed.  N=',i2,' /',i2,', E=',i2,'.')
      endif

      if (ic%gencr_inp.eq.0 .and. is_roll .and. .not.was_roll) then
         zerror = .true.
         write(lout, 2041) ic%tang, ic%gencr_inp
         write(   *, 2041) ic%tang, ic%gencr_inp
 2041    format (' Input: ERROR. In the first rolling case, new influence coefficients',/, &
                  15x,'are required.  T=',i3, ', C3=',i3,'.')
      endif

      if (ic%gencr_inp.eq.0 .and. .not.is_roll .and. was_roll) then
         zerror = .true.
         write(lout, 2042) ic%tang, ic%gencr_inp
         write(   *, 2042) ic%tang, ic%gencr_inp
 2042    format (' Input: ERROR. When changing from rolling to shifting, new influence',/, &
                15x, 'coefficients are required.  T=',i3,', C3=',i3,'.')
      endif

      ! Check consistency between different control integers F=1 needs N=1, F=3 needs N=1

      if ((ic%force1.eq.1 .or. ic%force1.eq.3) .and. ic%norm.ne.1) then
         zerror = .true.
         write(lout, 2051) ic%force1, ic%norm
         write(   *, 2051) ic%force1, ic%norm
 2051    format (' Input: ERROR. Total force F=',i2,' requires N=1 (',i2,').')
      endif

      ! check limitations of material model

      if (ic%mater.eq.1 .and. .not.is_ssrol) then
         zerror = .true.
         write(lout, 2061) ic%tang
         write(   *, 2061) ic%tang
 2061    format (' Input: ERROR. Visco-elastic materials require steady state rolling. T=',i3,'.')
      endif

      ! Check limitations of temperature model

      if (ic%heat.ge.1 .and. ic%tang.ne.3) then
         zerror = .true.
         write(lout, 2071) ic%heat, ic%tang
         write(   *, 2071) ic%heat, ic%tang
 2071    format (' Input: ERROR. Temperature calculation requires steady state rolling.',/,            &
                 '               H=',i3,', T='i3,'.')
      endif

      if (ic%frclaw_inp.eq.6 .and. ic%heat.eq.0) then
         zerror = .true.
         write(lout, 2073) ic%frclaw_inp, ic%heat
         write(   *, 2073) ic%frclaw_inp, ic%heat
 2073    format (' Input: ERROR. Temperature dep. friction requires that temperature is calculated.',/, &
                 '               L=',i3,', H='i3,'.')
      endif

      if (ic%frclaw_inp.eq.6 .and. (ic%mater.eq.2 .or. ic%mater.eq.3)) then
         zerror = .true.
         write(lout, 2074) ic%frclaw_inp, ic%mater
         write(   *, 2074) ic%frclaw_inp, ic%mater
 2074    format (' Input: ERROR. Temperature dep. friction not supported for FASTSIM.',/, &
                 '               L=',i3,', M='i3,'.')
      endif

      if (zerror) then
         call write_log(' Errors found. Aborting.')
         call abort_run()
      endif

      ! Store changes to ic also in packed Cpbtnfs etc.

      call ic_pack (1, cpbtnfs, vldcmze, xhgiaowr, ic)

      ! Copy the effective D-digit

      if (ic%discns1_inp.ge.2) ic%discns1_eff = ic%discns1_inp

      ! Copy the B- and M-digits to the material parameter data

      mater%bound_eff = ic%bound
      mater%mater_eff = ic%mater

      ! Copy the effective C3-digit to the material parameter data

      if (ic%gencr_inp.ge.2) mater%gencr_eff = ic%gencr_inp

      !------------------------------------------------------------------------------------------------------
      ! read & process debug parameters
      !------------------------------------------------------------------------------------------------------

      if (ic%xflow.ge.1) then
         call readline(linp, ncase, linenr, 'debug output psflrin', 'i', ints, dbles, flags, strngs,    &
                       mxnval, nval, ldebug, ieof, lstop, ierror)
         psflrin = ints(1)
      else
         psflrin = 1000000
      endif

      call ic_unpack_dbg (1, psflrin, ic)

      !------------------------------------------------------------------------------------------------------
      ! Read input for G-digit
      !------------------------------------------------------------------------------------------------------

      if (ic%gausei_inp.ne.1) then
         solv%gausei_eff = ic%gausei_inp

         call readline(linp, ncase, linenr, 'iteration constants', 'iiiidI', ints, dbles, flags,        &
                        strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
         solv%maxgs  = ints(1)  ! note: ordered alphabetically instead of inner-->outer
         solv%maxin  = ints(2)
         solv%maxnr  = ints(3)
         solv%maxout = ints(4)
         solv%eps    = dbles(1)
         if (nval.ge.6) discr%npot_max = min(1000000, max(100, ints(5)))

         if (ic%gausei_inp.eq.2 .or. ic%gausei_inp.eq.3) then
            call readline(linp, ncase, linenr, 'iteration (relaxation) parameters', 'ddid', ints,       &
                        dbles, flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
            solv%omegah = dbles(1)
            solv%omegas = dbles(2)
            solv%inislp = ints(1)
            solv%omgslp = dbles(3)
         endif

         if (ic%gausei_inp.eq.4) then
            call readline(linp, ncase, linenr, 'iteration (relaxation) parameters', 'id', ints,         &
                        dbles, flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
            solv%inislp = ints(1)
            solv%omgslp = dbles(1)
         endif

         if (ic%gausei_inp.eq.5) then
            call readline(linp, ncase, linenr, 'parameters for gdsteady', 'ddiddddd',                   &
                          ints, dbles, flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
            solv%fdecay = dbles(1)   ! <0.001 or >0.999 used for E_down / E_trl
            solv%betath = dbles(2)
            solv%kdowfb = ints(1)
            solv%d_ifc  = max(0.01d0, dbles(3))
            solv%d_lin  = dbles(4)
            solv%d_cns  = max(0.01d0, dbles(5))
            solv%d_slp  = max(0.01d0, dbles(6))
            solv%pow_s  = max(0.01d0, min(10d0, dbles(7)))

            ! check parameters diagonal scaling

            if (solv%d_lin*(solv%d_cns-solv%d_ifc).lt.0d0) then
               call write_log('ERROR: D_LIN and D_CNS-D_IFC should have same sign.')
               solv%d_lin = 0d0
               solv%d_cns = solv%d_ifc
            endif

            ! set variant used for GDsteady

            if (solv%fdecay.gt.0.999d0) then
               solv%gd_meth = 1     ! E_trl
            elseif (solv%fdecay.lt.0.001d0) then
               solv%gd_meth = 2     ! E_down(k)
               solv%kdown   = max(1, nint(-solv%fdecay))
            else
               solv%gd_meth = 3     ! E_keep(f)
            endif
         endif

         ! test the input quantities

         zerror = zerror .or. .not.check_irng ('MAXGS',  solv%maxgs , 1,   9999)
         zerror = zerror .or. .not.check_irng ('MAXIN',  solv%maxin , 1, 999999)
         zerror = zerror .or. .not.check_irng ('MAXNR',  solv%maxnr , 1,   9999)
         zerror = zerror .or. .not.check_irng ('MAXOUT', solv%maxout, 1,   9999)
      endif

      !------------------------------------------------------------------------------------------------------
      ! Read input for V- and L-digits
      !------------------------------------------------------------------------------------------------------

      if (ic%frclaw_inp.ne.1) then
         call fric_input(linp, ncase, linenr, ic%varfrc, ic%frclaw_inp, 1, fric, ldebug, ieof,          &
                lstop, zerror)
      endif

      ! Adapt the value of the friction coefficient used for scaling of tangential forces

      kin%use_muscal = ic%varfrc.eq.0
      if (kin%use_muscal) then
         kin%muscal = fric%fstat()
      else
         kin%muscal = 1d0
      endif

      !------------------------------------------------------------------------------------------------------
      ! Read input for C3 and M-digits
      !------------------------------------------------------------------------------------------------------

      ! Read elastic material constants when C3=2, 3 or 4

      if (ic%gencr_inp.eq.2 .or. ic%gencr_inp.eq.3 .or. ic%gencr_inp.eq.4) then

         if (ic%gencr_inp.eq.4) then
            call readline(linp, ncase, linenr, 'if-correction method', 'ii', ints ,dbles, flags,        &
                           strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
            mater%if_meth  = ints(1)
            mater%if_ver   = ints(2)
         endif

         call readline(linp, ncase, linenr, 'material properties', 'dddd', ints ,dbles, flags,          &
                        strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
         mater%poiss(1) = dbles(1)
         mater%poiss(2) = dbles(2)
         mater%gg(1)    = dbles(3)
         mater%gg(2)    = dbles(4)

         ! B=1: read flexibility for vertical compression of thin sheet

         if (ic%bound.eq.1) then
            call readline(linp, ncase, linenr, 'thin sheet compressibility',                            &
                          'd', ints, dbles, flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
            mater%flx_z = dbles(1)
         endif

         ! M=1: read viscoelastic material constants (Experimental; not tested)

         if (ic%mater.eq.1) then
            call readline(linp, ncase, linenr, 'viscoelastic material constants',                       &
                          'dddd', ints, dbles, flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
            mater%fg(1) = dbles(1)
            mater%fg(2) = dbles(2)
            mater%tc(1) = dbles(3)
            mater%tc(2) = dbles(4)
         endif

         ! M=2: read flexibility L1=L2=L for simplified theory and parameters of slope reduction

         if (ic%mater.eq.2) then
            call readline(linp, ncase, linenr, 'simplified theory flexibility + slope reduction',       &
                          'dddd', ints, dbles, flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
            mater%flx(1:3) = dbles(1)
            mater%k0_mf    = dbles(2)
            mater%alfamf   = dbles(3)
            mater%betamf   = dbles(4)
         endif

         ! M=3: read parameters of slope reduction for modified Fastsim algorithm

         if (ic%mater.eq.3) then
            call readline(linp, ncase, linenr, 'slope reduction parameters', 'ddd', ints, dbles,        &
                          flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
            mater%k0_mf    = dbles(1)
            mater%alfamf   = dbles(2)
            mater%betamf   = dbles(3)
         endif

         ! M=4: read parameters of elasto-plastic third body layer

         if (ic%mater.eq.4) then
            call readline(linp, ncase, linenr, 'third body layer parameters', 'dddd', ints, dbles,      &
                        flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
            mater%gg3    = dbles(1)
            mater%laythk = dbles(2)
            mater%tau_c0 = dbles(3)
            mater%k_tau  = dbles(4)
            if (mater%tau_c0.le.1d-10) mater%tau_c0 = 1d20
         endif
         if (ic%mater.eq.4) mater%flx(1:3) = mater%laythk / mater%gg3

         ! Check inputs for IF-correction

         if (ic%gencr_inp.eq.4) then
            zerror = zerror .or. .not.check_irng ('IF_METH', mater%if_meth, 0,   1)
            zerror = zerror .or. .not.check_irng ('VARIANT', mater%if_ver , 1,   4)
         endif

         ! Check material properties

         zerror = zerror .or. .not.check_range ('POISS1', mater%poiss(1), 0d0, 0.5d0)
         zerror = zerror .or. .not.check_range ('POISS2', mater%poiss(2), 0d0, 0.5d0)
         zerror = zerror .or. .not.check_range ('GG1', mater%gg(1), 1d-5, 1d20)
         zerror = zerror .or. .not.check_range ('GG2', mater%gg(2), 1d-5, 1d20)
         if (ic%bound.eq.1) then
            zerror = zerror .or. .not.check_range ('FLXZ', mater%flx_z, 1d-10, 1d20)
         endif
         if (ic%mater.eq.1) then
            zerror = zerror .or. .not.check_range ('FG1', mater%fg(1), 0d0, 1d20)
            zerror = zerror .or. .not.check_range ('FG2', mater%fg(2), 0d0, 1d20)
            zerror = zerror .or. .not.check_range ('TC1', mater%tc(1), 0d0, 1d20)
            zerror = zerror .or. .not.check_range ('TC2', mater%tc(2), 0d0, 1d20)
         endif
         if (ic%mater.eq.2) then
            zerror = zerror .or. .not.check_range ('FLX', mater%flx(1), 1d-10, 1d20)
         endif
         if (ic%mater.eq.2 .or. ic%mater.eq.3) then
            zerror = zerror .or. .not.check_range ('K0_MF', mater%k0_mf, 1d-10, 1d0)
            zerror = zerror .or. .not.check_range ('ALFAMF', mater%alfamf, 0d0, 1d0)
            zerror = zerror .or. .not.check_range ('BETAMF', mater%betamf, 0d0, 1d20)
         endif
         if (ic%mater.eq.4) then
            zerror = zerror .or. .not.check_range ('GG3',    mater%gg3, 1d-5, 1d20)
            zerror = zerror .or. .not.check_range ('LAYTHK', mater%gg3, 1d-10, 1d20)
         endif

         ! Compute combined material constants

         call combin_mater(mater)
         if (ic%mater.eq.1) then
            mater%vt(1) = mater%tc(1) * kin%veloc
            mater%vt(2) = mater%tc(2) * kin%veloc
         endif

      endif

      if (ic%tang.ne.0 .and. mater%gencr_eff.ne.4 .and. mater%gencr_eff.ne.9 .and. solv%maxout.ge.2     &
                                                                      .and. abs(mater%ak).lt.1d-4) then
         if (ic%gencr_inp.ne.0) then
            write(lout, 4011) mater%ak
            write(   *, 4011) mater%ak
         endif
 4011    format (' Input: WARNING. With equal material constants (Ak=',f7.4,'), MAXOUT is',/,           &
                 8x,'set to 1, i.e. the Panag/outer iteration process is skipped.')
         solv%maxout = 1
      endif

      if (ic%tang.ne.0 .and. solv%maxout.le.1 .and. (mater%gencr_eff.eq.4 .or. mater%gencr_eff.eq.9     &
                                                                      .or. abs(mater%ak).ge.1d-4)) then
         if (abs(mater%ak).ge.1d-4) then
            write(lout, 4013) mater%ak
            write(   *, 4013) mater%ak
 4013       format (' Input: WARNING. With dissimilar materials (Ak=',f7.4,'), MAXOUT should be > 1,',/, &
                    '        i.e. the Panag/outer iteration process should be used.')
         else
            write(lout, 4014) mater%gencr_eff
            write(   *, 4014) mater%gencr_eff
 4014       format (' Input: WARNING. With dissimilar shapes (C=',i2,'), MAXOUT should be > 1,',/,      &
                    '        i.e. the Panag/outer iteration process should be used.')
         endif
      endif

      !------------------------------------------------------------------------------------------------------
      ! Read input for H-digit: material parameters for temperature calculation
      !------------------------------------------------------------------------------------------------------

      if (ic%heat.eq.3) then
         call readline(linp, ncase, linenr, 'temperature inputs for body 1', 'dddd', ints, dbles,       &
                        flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
         mater%bktemp(1) = dbles(1)
         mater%heatcp(1) = dbles(2)
         mater%lambda(1) = dbles(3)
         mater%dens(1)   = dbles(4)

         call readline(linp, ncase, linenr, 'temperature inputs for body 2', 'dddd', ints, dbles,       &
                        flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
         mater%bktemp(2) = dbles(1)
         mater%heatcp(2) = dbles(2)
         mater%lambda(2) = dbles(3)
         mater%dens(2)   = dbles(4)

         zerror = zerror .or. .not.check_range ('HEATCP1', mater%heatcp(1), 1d-5, 1d20)
         zerror = zerror .or. .not.check_range ('LAMBDA1', mater%lambda(1), 1d-5, 1d20)
         zerror = zerror .or. .not.check_range ('DENS1',   mater%dens(1),  1d-12, 1d20)
         zerror = zerror .or. .not.check_range ('HEATCP2', mater%heatcp(2), 1d-5, 1d20)
         zerror = zerror .or. .not.check_range ('LAMBDA2', mater%lambda(2), 1d-5, 1d20)
         zerror = zerror .or. .not.check_range ('DENS2',   mater%dens(2),  1d-12, 1d20)
      endif

      if (ic%heat.eq.3 .and. ic%mater.eq.4) then
         call readline(linp, ncase, linenr, 'partitioning of plastic work', 'd', ints, dbles,           &
                        flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
         mater%betapl = dbles(1)

         zerror = zerror .or. .not.check_range ('BETAPL', mater%betapl, 0d0, 1d0)
      endif

      !------------------------------------------------------------------------------------------------------
      ! Read input for D-digit: potential contact, discretisation parameters
      !------------------------------------------------------------------------------------------------------

      if (ic%discns1_inp.ge.2 .and. ic%discns1_inp.le.9) then

         call readline(linp, ncase, linenr, 'discretisation parameters', 'dddaddD', ints, dbles,        &
                        flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
         discr%dx        = dbles(1)
         discr%ds        = dbles(2)
         discr%dqrel     = dbles(3)
         discr%angl_sep  = max(0d0, dbles(4))
         discr%dist_sep  = max(0d0, dbles(5))
         discr%dist_comb = dbles(6)
         if (nval.ge.7) then
            discr%dist_turn = dbles(7)
         else
            discr%dist_turn = 2d0 * discr%dist_sep - 1d0 * discr%dist_comb
         endif

         zerror = zerror .or. .not.check_range ('DX', discr%dx, 1d-8, 1d20)
         zerror = zerror .or. .not.check_range ('DS', discr%ds, 1d-8, 1d20)
         zerror = zerror .or. .not.check_range ('DQREL', discr%dqrel, 1d-3, 10d0)

         ! in transient cases, check element sizes of new/previous case

         if (ncase.gt.1 .and. wtd%numcps.ge.1 .and. ic%use_supergrid()) then
            associate(cp => wtd%allcps(1)%cp)
            dx_prv = cp%gd%cgrid_cur%dx
            ds_prv = cp%gd%cgrid_cur%dy
            if (abs(dx_prv-discr%dx).ge.1d-4*min(dx_prv, discr%dx) .or.                                 &
                abs(ds_prv-discr%ds).ge.1d-4*min(ds_prv, discr%ds)) then
               zerror = .true.
               write(lout,5001) ic%tang, dx_prv, ds_prv
               write(   *,5001) ic%tang, dx_prv, ds_prv
            endif
            end associate
 5001       format(' Input: ERROR. Transient contact (T=',i1,') needs constant DX, DY (',f6.3,',',f6.3, &
                   ') in subsequent cases')
         endif

         if (discr%dist_sep.lt.discr%dist_comb) then
            zerror = .true.
            write(lout, 5011) discr%dist_sep, discr%dist_comb
            write(   *, 5011) discr%dist_sep, discr%dist_comb
         endif
 5011    format (' Input: ERROR. D_SEP =',f8.2,' MUST BE >= D_COMB =',f8.2,'.')
!        if (discr%dist_sep.gt.2d0*discr%dist_comb) then
!           zerror = .true.
!           write(lout, 5012) discr%dist_sep, 2d0*discr%dist_comb
!           write(   *, 5012) discr%dist_sep, 2d0*discr%dist_comb
!        endif
!5012    format (' Input: ERROR. D_SEP =',f8.2,' MUST BE <= 2*D_COMB =',f8.2,'.')
      endif

      ! Set appropriate dq for shifts

      if (.not.is_roll) kin%dq    = 1d0
      if (.not.is_roll) kin%veloc = 1d0

      !------------------------------------------------------------------------------------------------------
      ! Read input for Z1-digit: track/roller rig design dimensions, profile(s), deviation
      !------------------------------------------------------------------------------------------------------

      if (ic%ztrack.eq.1 .or. ic%ztrack.eq.3) then

         if (ic%config.le.1) then

            ! get the track design dimensions (using 'd' for either gaugsq or raily0)

            call readline(linp, ncase, linenr, 'track design dimensions', 'ddda', ints, dbles, flags,    &
                           strngs, mxnval, nval, ldebug, ieof, lstop, ierror)

            trk%gauge_height = dbles(1)
            if (trk%gauge_height.gt.0d0) then
               trk%gauge_seqnum = nint(dbles(2))
               trk%track_gauge  = dbles(3)
            else
               trk%rail_y0      = dbles(2)
               trk%rail_z0      = dbles(3)
            endif
            trk%cant_angle   = dbles(4)

         else

            ! get the roller rig design dimensions (using 'd' for either gaugsq or raily0)

            call readline(linp, ncase, linenr, 'roller rig design dimensions', 'dddd', ints, dbles,     &
                           flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)

            trk%gauge_height = dbles(1)
            if (trk%gauge_height.gt.0d0) then
               trk%gauge_seqnum = nint(dbles(2))
               trk%track_gauge  = dbles(3)
            else
               trk%rail_y0      = dbles(2)
               trk%rail_z0      = dbles(3)
            endif
            trk%nom_radius   = dbles(4)

         endif

      endif

      if (ic%ztrack.eq.3) then

         ! Z1=3: one rail profile used for current side of track; get filename and configuration data

         call profile_read_config(trk%rai%prr, 'rail', linp, ncase, linenr, ldebug, ieof, lstop, zerror)

         ! check that variable profiles are used with absolute rail placement

         call profile_get_filetype(trk%rai%prr%fname, is_wheel, i_ftype, ldebug)

         if (i_ftype.eq.FTYPE_SLICES .and. trk%gauge_height.gt.0d0) then
            zerror = .true.
            write(lout, 5511) trk%gauge_height, trim(trk%rai%prr%fname)
            write(   *, 5511) trk%gauge_height, trim(trk%rai%prr%fname)
 5511       format (' Input: ERROR. The gauge point computation (GAUGHT =',f8.2,                        &
                    ') cannot be used when using a',/,15x,'variable profile ("',a,'").')
         endif

         ! read the rail profile

         if (.not.zerror) then
            call profile_read_file(trk%rai%prr, meta%dirnam, 0, ic%x_profil, ic%x_readln, lstop)
            zerror = zerror .or. (trk%rai%prr%ierror.ne.0)
         endif

      endif

      if (ic%ztrack.ge.2 .and. ic%ztrack.le.3) then

         ! read the track deviations for current side of the track

         call readline(linp, ncase, linenr, 'current rail deviations', 'ddadda', ints, dbles, flags,    &
                        strngs, mxnval, nval, ldebug, ieof, lstop, ierror)

         trk%rai%dy    = dbles(1)
         trk%rai%dz    = dbles(2)
         trk%rai%roll  = dbles(3)
         trk%rai%vy    = dbles(4)
         trk%rai%vz    = dbles(5)
         trk%rai%vroll = dbles(6)

      endif

      if (ic%ztrack.ge.2 .and. ic%ztrack.le.3 .and. ic%force1.eq.3) then

         ! read the track deflection parameters for current side of the track

         call readline(linp, ncase, linenr, 'current rail deflection parameters', 'dddd', ints, dbles,  &
                        flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)

         trk%ky_rail   = dbles(1)
         trk%fy_rail   = dbles(2)
         trk%kz_rail   = dbles(3)
         trk%fz_rail   = dbles(4)

      endif

      !------------------------------------------------------------------------------------------------------
      ! Read input for E-digit: wheel-set dimensions, profile, position and velocity
      !------------------------------------------------------------------------------------------------------

      ! get the wheel-set dimensions

      if (ic%ewheel.eq.3 .or. ic%ewheel.eq.5) then

         call readline(linp, ncase, linenr, 'wheel-set dimensions', 'ddd', ints, dbles, flags,          &
                        strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
         ws%flback_dist  = dbles(1)
         ws%flback_pos   = dbles(2)
         ws%nom_radius   = dbles(3)

         zerror = zerror .or. .not.check_range ('NOMRAD', ws%nom_radius, 1d-3, 1d20)

      endif

      ! get the wheel profile

      if (ic%ewheel.eq.3 .or. ic%ewheel.eq.5) then

         ! E1=3,5: get profile filename and configuration data

         call profile_read_config(ws%whl%prw, 'wheel', linp, ncase, linenr, ldebug, ieof, lstop, zerror)

         ! read the wheel profile, store in wheel data

         call profile_read_file(ws%whl%prw, meta%dirnam, 1, ic%x_profil, ic%x_readln, lstop)
         zerror = zerror .or. (ws%whl%prw%ierror.ne.0)

         if (.false. .and. ic%discns1_eff.eq.5) then
            call write_log(' D=5: converting wheel to variable profile...')
            call profile_make_varprof(ws%whl%prw, 1, ic%x_profil)
         endif

      endif

      ! get the wheel-set position and orientation

      if (ic%ewheel.ge.1 .and. ic%ewheel.le.5) then

         call readline(linp, ncase, linenr, 'wheel-set position and orientation', 'dddaaa', ints, &
                        dbles, flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
         if (ic%config.le.1) then
            ws%s      = dbles(1)
            ws%x      = 0d0
         else
            ws%s      = 0d0
            ws%x      = dbles(1)        ! C1=4,5: PITCH_ROL prescribed
         endif
         ws%y      = dbles(2)
         if (ic%norm.le.0) then
            ws%z      = dbles(3)
            ws%fz_inp = 0d0
         else
            ws%fz_inp = dbles(3)        ! N=1: FZ_TR prescribed
         endif
         ws%roll     = dbles(4)
         ws%yaw      = dbles(5)
         ws%pitch    = dbles(6)
      endif

      ! get the wheel-set velocity

      if (ic%ewheel.ge.2 .and. ic%ewheel.le.5) then

         if     (ic%config.le.1 .and. (ic%force1.le.0 .or. ic%force1.eq.3)) then
            types = 'dddaaa'    ! 1st VS_WS,      6th VPITCH_WS
         elseif (ic%config.le.1) then
            types = 'dddaad'    ! 1st VS_WS,      6th FX_WS/MY_WS
         elseif (ic%config.ge.4 .and. (ic%force1.le.0 .or. ic%force1.eq.3)) then
            types = 'addaaa'    ! 1st VPITCH_ROL, 6th VPITCH_WS
         elseif (ic%config.ge.4) then
            types = 'addaad'    ! 1st VPITCH_ROL, 6th FX_WS/MY_WS
         endif

         call readline(linp, ncase, linenr, 'wheel-set velocity and rotation', types, ints, dbles,      &
                        flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)

         if (ic%config.le.1) then
            ws%vs          = dbles(1)
            trk%vpitch_rol = 0d0
         else
            ws%vs          = 0d0
            trk%vpitch_rol = dbles(1)   ! C1=4,5: VPITCH_ROL prescribed
         endif
         ws%vy       = dbles(2)
         ws%vz       = dbles(3)
         ws%vroll    = dbles(4)
         ws%vyaw     = dbles(5)
         if (ic%force1.le.0 .or. ic%force1.eq.3) then
            ws%vpitch = dbles(6)
         elseif (ic%force1.eq.1) then
            ws%fx_inp = dbles(6)        ! F=1: FX_WS prescribed
         elseif (ic%force1.eq.2) then
            ws%my_inp = dbles(6)        ! F=2: MY_WS prescribed
         endif
      endif

      ! get flexible wheel-set deviations

      if (ic%ewheel.eq.4 .or. ic%ewheel.eq.5) then

         namside = 'current'

         ! read the wheel position deviations for current side of wheel-set

         call readline(linp, ncase, linenr, trim(namside) // ' wheel position deviations', 'dddaaa',    &
                       ints, dbles, flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)

         my_wheel%dx     = dbles(1)
         my_wheel%dy     = dbles(2)
         my_wheel%dz     = dbles(3)
         my_wheel%droll  = dbles(4)
         my_wheel%dyaw   = dbles(5)
         my_wheel%dpitch = dbles(6)

         ! read the wheel velocity deviations for current side of wheel-set

         call readline(linp, ncase, linenr, trim(namside) // ' wheel velocity deviations', 'dddaaa',    &
                       ints, dbles, flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)

         my_wheel%vx     = dbles(1)
         my_wheel%vy     = dbles(2)
         my_wheel%vz     = dbles(3)
         my_wheel%vroll  = dbles(4)
         my_wheel%vyaw   = dbles(5)
         my_wheel%vpitch = dbles(6)
      
      endif ! E=5,7

      !------------------------------------------------------------------------------------------------------
      ! read subsurface points
      !------------------------------------------------------------------------------------------------------

      if (ic%stress.ge.2) then
         call rdsubs(linp, ic, ncase, linenr, ldebug, subs)
      endif

      ! abort on errors

      if (zerror) then
         call write_log(' Errors found. Aborting.')
         call abort_run()
      endif

      if (ldebug.ge.3) call write_log('--- end subroutine wr_input ---')
      end associate

   end subroutine wr_input

!------------------------------------------------------------------------------------------------------------

   subroutine wr_write_inp(ncase, wtd)
!--purpose: Write input for the current W/R contact case to the inp-file, unit linp. 
      implicit none
!--subroutine arguments:
      integer                   :: ncase
      type(t_ws_track)          :: wtd
!--local variables:
      integer                  :: vldcmze, xhgiaowr, cpbtnfs, psflrin
      real(kind=8)             :: dflt_turn
      character(len=1)         :: namside
      type(t_wheel),   pointer :: my_wheel

      associate(ic    => wtd%ic,   mater => wtd%mater, discr => wtd%discr, kin   => wtd%kin,        &
                fric  => wtd%fric, solv  => wtd%solv,  ws    => wtd%ws,    trk   => wtd%trk  )

      if (ncase.gt.0) write(linp,'(a,i8)') '% Next case', ncase

      ! write control integers

      call ic_pack (1, cpbtnfs, vldcmze, xhgiaowr, ic)

      if (ic%xflow.le.0) then
         write(linp, 1101) cpbtnfs, vldcmze, xhgiaowr
      else
         write(linp, 1102) cpbtnfs, vldcmze, xhgiaowr
      endif

      ! write debug parameters

      if (ic%xflow.ge.1) then
         call ic_pack_dbg(psflrin, ic)
         write(linp, 1103) psflrin
      endif

      ! write parameters for the iterative solution algorithms

      if (ic%gausei_inp.ne.1) write(linp, 1201) solv%maxgs, solv%maxin, solv%maxnr, solv%maxout, solv%eps
      if (ic%gausei_inp.eq.2 .or. ic%gausei_inp.eq.3)                                                   &
         write(linp, 1211) solv%omegah, solv%omegas, solv%inislp, solv%omgslp
      if (ic%gausei_inp.eq.4) write(linp, 1212) solv%inislp, solv%omgslp
      if (ic%gausei_inp.eq.5) then
         if (solv%gd_meth.eq.1) solv%fdecay =  1d0
         if (solv%gd_meth.eq.2) solv%fdecay = -real(solv%kdown)
         write(linp, 1215) solv%fdecay, solv%betath, solv%kdowfb, solv%d_ifc, solv%d_lin, solv%d_cns,  &
            solv%d_slp, solv%pow_s
      endif

 1101 format (i8.7, 6x, 'C-P-B-T-N-F-S      CONFIG, PVTIME, BOUND,  TANG,   NORM,   FORCE,  STRESS', /,   &
              i8.7, 6x, 'V-L-D-C-M-Z-E      VARFRC, FRCLAW, DISCNS, INFLCF, MATER,  ZTRACK, EWHEEL', /,   &
              i8.7, 6x, 'H-G-I-A-O-W-R        HEAT, GAUSEI, IESTIM, MATFIL, OUTPUT, FLOW,   RETURN' )
 1102 format (i8.7, 6x, '  C-P-B-T-N-F-S         CONFIG, PVTIME, BOUND,  TANG,   NORM,   FORCE,  STRESS',/, &
              i8.7, 6x, '  V-L-D-C-M-Z-E         VARFRC, FRCLAW, DISCNS, INFLCF, MATER,  ZTRACK, EWHEEL',/, &
              i8.7, 6x, 'X-H-G-I-A-O-W-R   XFLOW,  HEAT, GAUSEI, IESTIM, MATFIL, OUTPUT, FLOW,   RETURN' )
 1103 format (i8.7, 6x, '  P-S-F-L-R-I_N         PROFIL, SMOOTH, FORCE,  LOCATE, READLN, INFLCF, NMDBG' )
 1201 format( 4i6, 3x, es8.1, 14x,   'MAXGS,  MAXIN,  MAXNR,  MAXOUT, EPS')
 1211 format( 2g12.4, i6, g12.4, 7x,  'OMEGAH, OMEGAS, INISLP, OMGSLP')
 1212 format(         i6, g12.4, 40x, 'INISLP, OMGSLP')
 1215 format( 6f8.3,  10x, 'FDECAY, D_IFC/LIN/CNS, D_SLP, POW_S')

      ! write friction description

      call fric_wrtinp(linp, ic%varfrc, ic%frclaw_inp, fric)

      ! write information for influence coefficients, esp. material constants

      if (ic%gencr_inp.eq.2 .or. ic%gencr_inp.eq.3 .or. ic%gencr_inp.eq.4) then

         ! write IF-correction method

         if (ic%gencr_inp.eq.4) then
            write(linp, 5100) mater%if_meth, mater%if_ver
 5100       format (2i6, 46x, 'IF_METH, VARIANT')
         endif

         ! write elastic material constants

         write(linp, 5200) mater%poiss, mater%gg
 5200    format (2(f8.4, 4x), 2g12.4, 10x, 'POISS 1,2,  GG 1,2')

         ! write visco-elastic material constants

         if (ic%mater.eq.1) then
            write(linp, 5301) mater%fg(1), mater%fg(2), mater%tc(1), mater%tc(2)
 5301       format (4g12.4, 10x, 'FG 1,2,  TC 1,2')
         endif

         ! write flexibility of simplified theory and slope reduction parameters

         if (ic%mater.eq.2) then
            write(linp, 5402) mater%flx(1), mater%k0_mf, mater%alfamf, mater%betamf
 5402       format (4g12.4, 10x, 'FLX, K0_MF, ALFA, BETA')
         elseif (ic%mater.eq.3) then
            write(linp, 5403) mater%k0_mf, mater%alfamf, mater%betamf
 5403       format (3g12.4, 22x, 'K0_MF, ALFAMF, BETAMF')
         endif

         ! write parameters of elasto-plastic interface layer

         if (ic%mater.eq.4) then
            write(linp, 5501) mater%gg3, mater%laythk, mater%tau_c0, mater%k_tau
 5501       format (4g12.4, 10x, 'GG3, LAYTHK, TAU_C0, K_TAU')
         endif
      endif

      ! write material parameters for temperature model

      if (ic%heat.eq.3) then
         write(linp,5611) mater%bktemp(1), mater%heatcp(1), mater%lambda(1), mater%dens(1), 1,1,1,1
         write(linp,5611) mater%bktemp(2), mater%heatcp(2), mater%lambda(2), mater%dens(2), 2,2,2,2
 5611    format (4g12.4, 10x, 'BKTEMP',i1,', HEATCP',i1,', LAMBDA',i1,', DENSITY',i1)
      endif

      if (ic%heat.eq.3 .and. ic%mater.eq.4) then
         write(linp,5612) mater%betapl
         write(linp,5612) mater%betapl
 5612    format (1g12.4, 46x, 'BETAPL')
      endif

      ! write information needed for the grid discretization

      if (ic%discns1_inp.ge.2 .and. ic%discns1_inp.le.9) then
         ! write dist_turn only if different from default value
         dflt_turn = 2d0 * discr%dist_sep - 1d0 * discr%dist_comb
         if (abs(discr%dist_turn-dflt_turn).le.1d-3) then
            write(linp, 6101) discr%dx, discr%ds, discr%dqrel, discr%angl_sep*180d0/pi, discr%dist_sep, &
                   discr%dist_comb
         else
            write(linp, 6102) discr%dx, discr%ds, discr%dqrel, discr%angl_sep*180d0/pi, discr%dist_sep, &
                   discr%dist_comb, discr%dist_turn
         endif
 6101    format( 3f9.4, f8.1,'d', 2f9.4, 4x, 'DX, DS, DQREL, A_SEP, D_SEP, D_COMB')
 6102    format( 3f9.4, f8.1,'d', 3f6.2, 4x, 'DX, DS, DQREL, A_SEP, D_SEP, D_COMB, D_TURN')
      endif

      ! write information on track geometry

      if (ic%ztrack.eq.1 .or. ic%ztrack.eq.3) then

         if     (ic%config.le.1 .and. trk%gauge_height.gt.0d0) then
            write(linp, 7101) trk%gauge_height, trk%gauge_seqnum, trk%track_gauge, trk%cant_angle, 'CANT'
         elseif (                     trk%gauge_height.gt.0d0) then
            write(linp, 7101) trk%gauge_height, trk%gauge_seqnum, trk%track_gauge, trk%nom_radius, 'NOMRAD'
         elseif (ic%config.le.1) then
            write(linp, 7102) trk%gauge_height, trk%rail_y0, trk%rail_z0, trk%cant_angle, 'CANT'
         else
            write(linp, 7102) trk%gauge_height, trk%rail_y0, trk%rail_z0, trk%nom_radius, 'NOMRAD'
         endif
 7101    format ( g14.6, i6, 2g14.6, 10x, 'GAUGHT, GAUGSQ, GAUGWD, ',a)
 7102    format (4g14.6,  2x, 'GAUGHT, RAILY0, RAILZ0, ',a)

      endif

      ! write rail profile filename(s)

      if (ic%ztrack.eq.3) call profile_write_config(wtd%meta, trk%rai%prr, 0)

      ! write rail deviations

      if (ic%ztrack.ge.2 .and. ic%ztrack.le.3) then

         write(linp, 7201) trk%rai%dy, trk%rai%dz, trk%rai%roll, trk%rai%vy, trk%rai%vz, trk%rai%vroll

 7201    format (3g14.6, 16x, '% DY, DZ, DROLLR', /, 3g14.6, 16x  '  VY, VZ, VROLL', :,i1)

      endif

      ! write massless rail deflection parameters

      if (ic%ztrack.eq.3 .and. ic%force1.eq.3) then

         write(linp, 7301) trk%ky_rail, trk%fy_rail, trk%kz_rail, trk%fz_rail

 7301    format (4g14.6, 2x, 'KYRAIL, FYRAIL, KZRAIL, FZRAIL')

      endif

      ! write information on the wheelset geometry

      if (ic%ewheel.eq.3 .or. ic%ewheel.eq.5) then

         write(linp, 8101) ws%flback_dist, ws%flback_pos, ws%nom_radius
 8101    format (3g14.6, 16x, 'FBDIST, FBPOS, NOMRAD')

      endif

      ! write wheel profile filename

      if (ic%ewheel.eq.3 .or. ic%ewheel.eq.5) then
         call profile_write_config(wtd%meta, ws%whl%prw, 1)
      endif

      ! write information on the wheelset position and velocity
      ! Using kinematic constants fz, fx, fy

      if (ic%ewheel.ge.1) then
         if (ic%config.le.1 .and. ic%norm.eq.0) then
            write(linp, 8201) ws%s, ws%y, ws%z,      'S', 'Z'
         elseif (ic%config.le.1) then
            write(linp, 8201) ws%s, ws%y, ws%fz_inp, 'S', 'FZ'
         elseif (ic%config.ge.4 .and. ic%norm.eq.0) then
            write(linp, 8201) ws%x, ws%y, ws%z,      'X', 'Z'
         elseif (ic%config.ge.4) then
            write(linp, 8201) ws%x, ws%y, ws%fz_inp, 'X', 'FZ'
         endif
         write(linp, 8203) ws%roll, ws%yaw, ws%pitch
 8201    format( 3(g14.7, 1x), 13x, '% ',a,', Y, ',a)
 8203    format( 3(g14.7, 1x), 13x, '% ROLL, YAW, PITCH')
      endif

      if (ic%ewheel.ge.2) then
         if (ic%config.le.1) then
            write(linp, 8301) ws%vs,          ws%vy,     ws%vz, 'VS'
         else
            write(linp, 8301) trk%vpitch_rol, ws%vy,     ws%vz, 'VPITCH_ROL'
         endif
 8301    format( 3(g14.7, 1x), 13x, '% ',a,', VY, VZ')

         if (ic%force1.eq.0 .or. ic%force1.eq.3) then
            write(linp, 8303) ws%vroll, ws%vyaw,   ws%vpitch, 'VPITCH'
         elseif (ic%force1.eq.1) then
            write(linp, 8303) ws%vroll, ws%vyaw,   ws%fx_inp, 'FX'
         elseif (ic%force1.eq.2) then
            write(linp, 8303) ws%vroll, ws%vyaw,   ws%my_inp, 'MY'
         endif
 8303    format( 2(g14.7, 1x), g16.9, 12x, '% VROLL, VYAW, ',a)
      endif

      ! write flexible wheelset deviations

      if (ic%ewheel.eq.4 .or. ic%ewheel.eq.5) then
         namside = ' '
         my_wheel => ws%whl
         write(linp, 8401) my_wheel%dx, my_wheel%dy, my_wheel%dz, namside,                              &
                           my_wheel%droll, my_wheel%dyaw, my_wheel%dpitch, namside,                     &
                           my_wheel%vx, my_wheel%vy, my_wheel%vz, namside,                              &
                           my_wheel%vroll, my_wheel%vyaw, my_wheel%vpitch, namside

 8401    format (3g14.6, 16x, '% DX, DY, DZ',a, /, 3g14.6, 16x, '  DROLL, DYAW, DPITCH',a, /,           &
                 3g14.6, 16x, '% VX, VY, VZ',a, /, 3g14.6, 16x, '  VROLL, VYAW, VPITCH',a)
      endif

      ! write the points in which the subsurface elastic field must be calculated. 

      if (ic%stress.ge.2) call wrsubs(linp, ic, wtd%subs)

      ! write empty line to separate from next case

      write(linp, '(1x)')

      end associate
   end subroutine wr_write_inp

!------------------------------------------------------------------------------------------------------------

end module m_wr_input
