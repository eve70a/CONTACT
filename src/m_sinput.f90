!------------------------------------------------------------------------------------------------------------
! m_sinput - read one case from the input-file and check the input
!
! Copyright 1979-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_sinput

use m_hierarch_data
use m_subsurf
use m_friclaw

implicit none
private

public input
public fric_input
public wrtinp

contains

!------------------------------------------------------------------------------------------------------------

   subroutine input (inp, linenr, meta, ic, potcon, mater, hz, cgrid, geom, fric, kin, solv, subs)
!--purpose: Input-routine for the Non-Hertzian case. Read all problem variables, unpack, and test.
      implicit none
!--subroutine parameters:
      integer          :: inp, linenr
      type(t_metadata) :: meta
      type(t_ic)       :: ic
      type(t_potcon)   :: potcon
      type(t_hertz)    :: hz
      type(t_grid)     :: cgrid
      type(t_geomet)   :: geom
      type(t_material) :: mater
      type(t_friclaw)  :: fric
      type(t_kincns)   :: kin
      type(t_solvers)  :: solv
      type(t_subsurf)  :: subs
!--local variables:
      integer, parameter :: mxnval = 20
      logical, parameter :: lstop  = .true.
      type(t_ic)     :: icold
      integer        :: ncase, vldcmze, xgiaowr, pbtnfs, nval, ieof, ierror
      integer        :: ii, iof, ip, j, k, line0, idebug, mxold, myold, nn, npot
      logical        :: z, zerror, is_roll, is_ssrol, was_roll, flags(mxnval)
      real(kind=8)   :: akeff
      integer        :: ints(mxnval), idum(1)
      real(kind=8)   :: dbles(mxnval)
      real(kind=8), dimension(:), allocatable :: tmp
      character(len=256) :: strngs(mxnval), spaces

      ncase = meta%ncase
      ieof  = -1 ! eof=error

      idebug = 0
      if (idebug.ge.1) write(*,'(a)') '--- Start subroutine input.f ---'
      line0 = linenr

      ! save old ic for checks:

      icold = ic

      !------------------------------------------------------------------------------------------------------
      ! read & check control integers
      !------------------------------------------------------------------------------------------------------

      call readLine(linp, ncase, linenr, 'control integers pbtnfs', 'i', ints, dbles, flags, strngs,    &
                    mxnval, nval, idebug, ieof, lstop, ierror)
      pbtnfs = ints(1)

      call readline(linp, ncase, linenr, 'control integers vldcmze', 'iI', ints, dbles, flags, strngs,  &
                    mxnval, nval, idebug, ieof, lstop, ierror)
      vldcmze = ints(1)
      if (nval.ge.2) ic%gapwgt = ints(2)

      call readline(linp, ncase, linenr, 'control integers giaowr', 'i', ints, dbles, flags, strngs,    &
                    mxnval, nval, idebug, ieof, lstop, ierror)
      xgiaowr = ints(1)

      call ic_unpack (3, pbtnfs, vldcmze, xgiaowr, ic)

      ! Determine whether T describes rolling or not

      is_roll  = ic%tang.eq.2 .or. ic%tang.eq.3
      is_ssrol = ic%tang.eq.3
      was_roll = icold%tang.eq.2 .or. icold%tang.eq.3

      ! Set "sens" flag, calculation of sensitivities

      ic%sens = 1
      if (.false. .and. ic%output_surf.ge.2) ic%sens = 3

      ! TEST  the input quantities (control integers)

      zerror = .false.

      zerror = zerror .or.           .not.check_irng ('Control digit P',  ic%pvtime, 0, 2)
      zerror = zerror .or.           .not.check_irng ('Control digit B',  ic%bound , 0, 6)
      zerror = zerror .or.           .not.check_irng ('Control digit T',  ic%tang  , 0, 3)
      zerror = zerror .or.           .not.check_irng ('Control digit N',  ic%norm  , 0, 1)
      zerror = zerror .or.           .not.check_irng ('Control digit F',  ic%force , 0, 2)
      zerror = zerror .or.           .not.check_irng ('Control digit S',  ic%stress, 0, 3)

      zerror = zerror .or.           .not.check_2int ('Control digit V',  ic%varfrc, 0, 2)
      zerror = zerror .or. .not.check_2rng ('Control digit L',  ic%frclaw_inp, 0, 4, 6, 6)
      zerror = zerror .or.       .not.check_irng ('Control digit D',  ic%discns_inp, 0, 2)
      zerror = zerror .or.           .not.check_irng ('gap weighting',    ic%gapwgt, 0, 2)
      zerror = zerror .or. .not.check_2rng ('Control digit C3', ic%gencr_inp , 0, 4, 9, 9)
      zerror = zerror .or.           .not.check_irng ('Control digit M',  ic%mater , 0, 4)
      zerror = zerror .or.           .not.check_irng ('Control digit Z',  ic%rznorm, 0, 2)
      zerror = zerror .or.     .not.check_2rng ('Control digit E',  ic%rztang, 0, 1, 9, 9)

      zerror = zerror .or.           .not.check_irng ('Control digit X',  ic%nmdbg , 0, 9)
      zerror = zerror .or.     .not.check_2rng ('Control digit H',  ic%heat  , 0, 1, 3, 3)
      zerror = zerror .or.       .not.check_irng ('Control digit G',  ic%gausei_inp, 0, 5)
      zerror = zerror .or.           .not.check_irng ('Control digit I',  ic%iestim, 0, 3)
      zerror = zerror .or.      .not.check_irng ('Control digit A',  ic%matfil_surf, 0, 2)
      zerror = zerror .or.      .not.check_irng ('Control digit O',  ic%output_surf, 0, 5)
      zerror = zerror .or.           .not.check_irng ('Control digit W',  ic%flow  , 0, 9)
      zerror = zerror .or.           .not.check_irng ('Control digit R',  ic%return, 0, 3)

      ! Check requirements for first case

      if (ncase.eq.1 .and. ic%pvtime.le.1) then
         write(lout, 2001) ic%pvtime
         write(   *, 2001) ic%pvtime
 2001    format (' Input: WARNING. In the first case contact must be initiated.'/,                     &
                 '                 Digit P=', i3,' is overruled, set to 2.')
         ic%pvtime = 2
      endif

      if (ncase.eq.1 .and. ic%frclaw_inp.eq.1) then
         zerror = .true.
         write(lout, 2002) ic%frclaw_inp
         write(   *, 2002) ic%frclaw_inp
 2002    format (' Input: ERROR. In the first case L must be.ne.1'/,                                   &
                 '               Digit L=', i3)
      endif

      if (ncase.eq.1 .and. ic%discns_inp.eq.0) then
         zerror = .true.
         write(lout, 2003) ic%discns_inp
         write(   *, 2003) ic%discns_inp
 2003    format (' Input: ERROR. In the first case D must be.ne.0'/,                                   &
                 '               Digit D=', i3)
      endif

      if (ncase.eq.1 .and. ic%iestim.ne.0) then
         write(lout, 2004) ic%iestim
         write(   *, 2004) ic%iestim
 2004    format (' Input: WARNING. In the first case contact must be initiated.'/,                     &
                 '                 Digit I=', i3,' is overruled, set to 0.')
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

      if (ic%discns_inp.ne.0 .and. ic%gencr_inp.eq.0) then
         zerror = .true.
         write(lout, 2011) ic%discns_inp, ic%gencr_inp
         write(   *, 2011) ic%discns_inp, ic%gencr_inp
 2011    format (' Input: ERROR. A new discretisation needs new influence functions'/,                 &
                 '               Digits D, C=', 2i3)
      endif

      if (ic%discns_inp.ne.0 .and. ic%rznorm.eq.0) then
         zerror = .true.
         write(lout, 2012) ic%discns_inp, ic%rznorm
         write(   *, 2012) ic%discns_inp, ic%rznorm
 2012    format (' Input: ERROR. A new discretisation needs new undef. distances'/,                    &
                 '               Digits D, Z=', 2i3)
      endif

      if (is_ssrol .and. ic%pvtime.ne.2) then
         ic%pvtime = 2
         write(lout, 2021) ic%pvtime
         write(   *, 2021) ic%pvtime
 2021    format (' Input: WARNING. In steady state rolling previous Pv is ignored: digit P :=', i3)
      endif

      if (ic%gencr_inp.eq.0 .and. is_roll .and. .not.was_roll) then
         zerror = .true.
         write(lout, 2041) ic%tang, ic%gencr_inp
         write(   *, 2041) ic%tang, ic%gencr_inp
 2041    format (' Input: ERROR. In the first rolling case, new influence coefficients',/,              &
                 '               are required.  T=',i3,', C=',i3,'.')
      endif

      if (ic%gencr_inp.eq.0 .and. .not.is_roll .and. was_roll) then
         zerror = .true.
         write(lout, 2042) ic%tang, ic%gencr_inp
         write(   *, 2042) ic%tang, ic%gencr_inp
 2042    format (' Input: ERROR. When changing from rolling to shifting, new influence',/,              &
                    '               coefficients are required.  T=',i3,', C=',i3,'.')
      endif

      ! check limitations of material model

      if (ic%mater.eq.1 .and. .not.is_ssrol) then
         zerror = .true.
         write(lout, 2061) ic%tang
         write(   *, 2061) ic%tang
 2061    format (' Input: ERROR. Visco-elastic materials require steady state rolling. T=',i3,'.')
      endif

      if ((ic%mater.eq.2 .or. ic%mater.eq.3) .and. ic%varfrc.ne.0) then
         zerror = .true.
         write(lout, 2063) ic%mater, ic%varfrc
         write(   *, 2063) ic%mater, ic%varfrc
 2063    format (' Input: ERROR. Friction variation is not supported for FASTSIM.',/,                   &
                 '               M=',i3,', V='i3,'.')
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

      ! Store changes to ic also in packed pbtnfs etc.

      call ic_pack (3, pbtnfs, vldcmze, xgiaowr, ic)

      ! Copy the effective D-digit

      if (ic%discns_inp.ge.2) ic%discns_eff = ic%discns_inp

      ! Copy the B- and M-digits to the material parameter data

      mater%bound_eff = ic%bound
      mater%mater_eff = ic%mater

      ! Copy the effective C-digit to the material parameter data

      if (ic%gencr_inp.ge.2) mater%gencr_eff = ic%gencr_inp

      !------------------------------------------------------------------------------------------------------
      ! read & check iteration parameters
      !------------------------------------------------------------------------------------------------------

      if (ic%gausei_inp.ne.1) then
         call readline(linp, ncase, linenr, 'iteration constants', 'iiiid',                             &
                       ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
         solv%maxgs  = ints(1)  ! note: ordered alphabetically instead of inner-->outer
         solv%maxin  = ints(2)
         solv%maxnr  = ints(3)
         solv%maxout = ints(4)
         solv%eps    = dbles(1)

         solv%gausei_eff = ic%gausei_inp
      endif

      if (ic%gausei_inp.eq.2 .or. ic%gausei_inp.eq.3) then
         call readline(linp, ncase, linenr, 'iteration (relaxation) parameters', 'ddid',                &
                       ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
         solv%omegah = dbles(1)
         solv%omegas = dbles(2)
         solv%inislp = ints(1)
         solv%omgslp = dbles(3)
      endif

      if (ic%gausei_inp.eq.4) then
         call readline(linp, ncase, linenr, 'iteration (relaxation) parameters', 'id',                  &
                       ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
         solv%inislp = ints(1)
         solv%omgslp = dbles(1)
      endif

      if (ic%gausei_inp.eq.5) then
         call readline(linp, ncase, linenr, 'parameters for gdsteady', 'ddiddddd',                      &
                       ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
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

      zerror = zerror .or. .not.check_irng ('MAXGS',  solv%maxgs , 1,   9999)
      zerror = zerror .or. .not.check_irng ('MAXIN',  solv%maxin , 1, 999999)
      zerror = zerror .or. .not.check_irng ('MAXNR',  solv%maxnr , 1,   9999)
      zerror = zerror .or. .not.check_irng ('MAXOUT', solv%maxout, 1,   9999)

      !------------------------------------------------------------------------------------------------------
      ! read & check kinematic constants
      !------------------------------------------------------------------------------------------------------

      call readline(linp, ncase, linenr, 'kinematic inputs', 'dddd',                                    &
                    ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)

      if (ic%norm.eq.0) then
         kin%pen    = dbles(1)
      else
         kin%fntrue = dbles(1)
      endif
      if (ic%force.eq.0) then
         kin%cksi   = dbles(2)
         kin%ceta   = dbles(3)
      elseif (ic%force.eq.1) then
         kin%fxrel1 = dbles(2)
         kin%ceta   = dbles(3)
      elseif (ic%force.eq.2) then
         kin%fxrel1 = dbles(2)
         kin%fyrel1 = dbles(3)
      endif
      kin%cphi = dbles(4)

      ! Check total forces Fn, Fx, Fy

      if (ic%norm.eq.1 .and. kin%fntrue.le.0d0) then
         zerror = zerror .or. .not.check_range ('the normal force', kin%fntrue, 0d0, 1d20)
      endif

      if (ic%force.ge.1 .and. abs(kin%fxrel1).gt.1d0) then
         zerror = zerror .or. .not.check_range ('Abs(Fx)', abs(kin%fxrel1), -1d0, 1d0)
      endif

      z = (kin%fxrel1 **2 + kin%fyrel1 **2).gt.1d0
      if (ic%force.eq.2 .and. z) then
         zerror = .true.
         write(lout, 2101) dsqrt (kin%fxrel1 **2 + kin%fyrel1 **2)
         write(   *, 2101) dsqrt (kin%fxrel1 **2 + kin%fyrel1 **2)
 2101    format (' Input: ERROR. Ft cannot be.gt.FStat*Fn'/,                                            &
                 ' The values prescribed by Fx and Fy give Ft = ',g12.4,'*FStat*Fn')
      endif

      !------------------------------------------------------------------------------------------------------
      ! read & check input for V- and L-digits for V = 0, single set NVF = 1
      !------------------------------------------------------------------------------------------------------

      if (ic%frclaw_inp.ne.1 .and. ic%varfrc.eq.0) then
         call fric_input(linp, ncase, linenr, ic%varfrc, ic%frclaw_inp, 1, fric, idebug, ieof,          &
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
      ! read rolling step, direction, material constants when C=2, 3, 4 or 9
      !------------------------------------------------------------------------------------------------------

      if ((ic%gencr_inp.ge.2 .and. ic%gencr_inp.le.4) .or. ic%gencr_inp.eq.9) then

         ! in case of rolling (T=2,3): read the rolling step

         if (is_roll) then
            call readline(linp, ncase, linenr, 'rolling velocity, direction and distance/time step',    &
                          'addD', ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
            kin%chi   = dbles(1)
            kin%dq    = dbles(2)
            kin%veloc = dbles(3)
            if (nval.ge.4) then
               kin%facphi = dbles(4)
            else
               kin%facphi = 1d0 / 6d0
            endif
         endif

         ! C=4: read surface inclinations needed for the Blanco-approach

         if (ic%gencr_inp.eq.4) then

            ! read number of points

            call readline(linp, ncase, linenr, 'if-correction method', 'iii', ints ,dbles, flags,       &
                           strngs, mxnval, nval, idebug, ieof, lstop, ierror)
            mater%if_meth  = ints(1)
            mater%if_ver   = ints(2)
            nn             = ints(3)

            ! check array-size, reallocate array

            zerror = zerror .or. .not.check_irng ('NN', nn, 1, 9999)
            if (allocated(mater%surf_inclin)) deallocate(mater%surf_inclin)
            allocate(mater%surf_inclin(nn,2))

            ! read surface inclinations

            allocate(tmp(nn*2))
            call read1darr(linp, ncase, linenr, 'conformal surface inclinations', 'a', nn*2, idum, tmp, &
                           idebug, lstop, ierror)

            mater%ninclin = nn
            do ii = 1, nn
               mater%surf_inclin(ii,1) = tmp(2*ii-1)
               mater%surf_inclin(ii,2) = tmp(2*ii  )
               ! write(*,*) 'i=',ii,': yi=',tmp(2*ii-1),', ai=',tmp(2*ii)
            enddo
            deallocate(tmp)
         endif

         ! C=9, numerical influence coefficients: read filename

         if (ic%gencr_inp.eq.9) then
            call readline(linp, ncase, linenr, 'file with numerical influence coefficients',            &
                          's', ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
            mater%fname_influe = trim(strngs(1))
         endif

         ! read elastic material constants

         call readline(linp, ncase, linenr, 'material properties',                                      &
                       'dddd', ints ,dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
         mater%poiss(1) = dbles(1)
         mater%poiss(2) = dbles(2)
         mater%gg(1) = dbles(3)
         mater%gg(2) = dbles(4)

         ! B=1: read flexibility for vertical compression of thin sheet

         if (ic%bound.eq.1) then
            call readline(linp, ncase, linenr, 'thin sheet compressibility',                            &
                          'd', ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
            mater%flx_z = dbles(1)
         endif

         ! M=1: read viscoelastic material constants

         if (ic%mater.eq.1) then
            call readline(linp, ncase, linenr, 'viscoelastic material constants',                       &
                          'dddd', ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
            mater%fg(1) = dbles(1)
            mater%fg(2) = dbles(2)
            mater%tc(1) = dbles(3)
            mater%tc(2) = dbles(4)
         endif

         ! M=2: read flexibility L1=L2=L for simplified theory and parameters of slope reduction

         if (ic%mater.eq.2) then
            call readline(linp, ncase, linenr, 'simplified theory flexibility + slope reduction',       &
                          'dddd', ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
            mater%flx(1:3) = dbles(1)
            mater%k0_mf    = dbles(2)
            mater%alfamf   = dbles(3)
            mater%betamf   = dbles(4)
         endif

         ! M=3: read parameters of slope reduction for modified Fastsim algorithm

         if (ic%mater.eq.3) then
            call readline(linp, ncase, linenr, 'slope reduction parameters', 'ddd', ints, dbles,        &
                          flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
            mater%k0_mf    = dbles(1)
            mater%alfamf   = dbles(2)
            mater%betamf   = dbles(3)
         endif

         ! M=4: read parameters of elasto-plastic third body layer

         if (ic%mater.eq.4) then
            call readline(linp, ncase, linenr, 'third body layer parameters',                           &
                          'dddd', ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
            mater%gg3    = dbles(1)
            mater%laythk = dbles(2)
            mater%tau_c0 = dbles(3)
            mater%k_tau  = dbles(4)
            if (mater%tau_c0.le.1d-10) mater%tau_c0 = 1d20
         endif

         if (ic%mater.eq.4) mater%flx(1:3) = mater%laythk / mater%gg3

         ! Check kinematic constants for shifts (T=0,1,5) and rolling (T=2,3)

         if (is_roll) then
            zerror = zerror .or. .not.check_range ('DQ', kin%dq, solv%eps, 1d20)
            zerror = zerror .or. .not.check_range ('VELOC', kin%veloc, solv%eps, 1d20)

            if (ic%mater.eq.1 .and. abs(kin%chi).gt.0.01d0) then
               zerror = .true.
               write(lout, 3511) kin%chi
               write(   *, 3511) kin%chi
 3511          format (' Input: ERROR. When using visco-elastic materials, the rolling direction',/     &
                       ' (Chi=',f6.2,') must be 0.')
            endif
         endif

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

         ! compute combined material constants

         call combin_mater(mater)
         if (ic%mater.eq.1) then
            mater%vt(1) = mater%tc(1) * kin%veloc
            mater%vt(2) = mater%tc(2) * kin%veloc
         endif

      ! endif (C=2, material constants)

      endif

      akeff = abs(mater%ak)

      if (ic%tang.ne.0 .and. mater%gencr_eff.ne.4 .and. mater%gencr_eff.ne.9 .and.                      &
                                                   solv%maxout.ge.2 .and. abs(akeff).lt.1d-4) then
         if (ic%gencr_inp.ne.0) then
            write(lout, 4011) akeff
            write(   *, 4011) akeff
         endif
 4011    format (' Input: WARNING. With equal material constants (Ak=',f7.4,'), MAXOUT is',/,           &
                 '        set to 1, i.e. the Panag/outer iteration process is skipped.')
         solv%maxout = 1
      endif

      if (ic%tang.ne.0 .and. solv%maxout.le.1 .and. (mater%gencr_eff.eq.4 .or. mater%gencr_eff.eq.9     &
                                                                        .or. abs(akeff).ge.1d-4)) then
         if (abs(akeff).ge.1d-4) then
            write(lout, 4013) mater%ak
            write(   *, 4013) mater%ak
 4013       format (' Input: WARNING. With dissimilar materials (Ak=',f7.4,'), MAXOUT should',/,        &
                    '        be > 1, i.e. the Panag/outer iteration process should be used.')
         else
            write(lout, 4014) mater%gencr_eff
            write(   *, 4014) mater%gencr_eff
 4014       format (' Input: WARNING. With dissimilar shapes (C=',i2,'), MAXOUT should',/,              &
                    '        be > 1, i.e. the Panag/outer iteration process should be used.')
         endif
      endif

      !------------------------------------------------------------------------------------------------------
      ! read material parameters for temperature calculation
      !------------------------------------------------------------------------------------------------------

      if (ic%heat.eq.3) then
         call readline(linp, ncase, linenr, 'temperature inputs for body 1', 'dddd', ints, dbles,       &
                        flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
         mater%bktemp(1) = dbles(1)
         mater%heatcp(1) = dbles(2)
         mater%lambda(1) = dbles(3)
         mater%dens(1)   = dbles(4)

         call readline(linp, ncase, linenr, 'temperature inputs for body 2', 'dddd', ints, dbles,       &
                        flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
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
                        flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
         mater%betapl = dbles(1)

         zerror = zerror .or. .not.check_range ('BETAPL', mater%betapl, 0d0, 1d0)
      endif

      !------------------------------------------------------------------------------------------------------
      ! read discretisation parameters
      !------------------------------------------------------------------------------------------------------

      if (ic%discns_inp.eq.2) then
         mxold = potcon%mx
         myold = potcon%my
         call readline(linp, ncase, linenr, 'specification method for potential contact area',          &
                       'i', ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
         potcon%ipotcn = ints(1)

         if (potcon%ipotcn.ge.-5 .and. potcon%ipotcn.le.-1) then

         ! Hertzian options for potential contact:

            call readline(linp, ncase, linenr, 'properties of potential contact area',                  &
                          'iiddd', ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
            potcon%mx = ints(1)
            potcon%my = ints(2)
            if (potcon%ipotcn.eq.-1) then
               hz%a1   = dbles(1)
               hz%b1   = dbles(2)
            elseif (potcon%ipotcn.eq.-2) then
               hz%a1   = dbles(1)
               hz%aob  = dbles(2)
            elseif (potcon%ipotcn.eq.-3) then
               hz%aa   = dbles(1)
               hz%bb   = dbles(2)
            elseif (potcon%ipotcn.eq.-4) then
               hz%a1   = dbles(1)
               hz%bb   = dbles(2)
            elseif (potcon%ipotcn.eq.-5) then
               hz%aa   = dbles(1)
               hz%bb   = dbles(2)
            endif
            hz%scale = dbles(3)

         elseif (potcon%ipotcn.eq.-6) then

         ! SDEC option for potential contact:

            call readline(linp, ncase, linenr, 'properties of SDEC contact area',                       &
                          'iidddd', ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
            potcon%mx = ints(1)
            potcon%my = ints(2)
            hz%aa     = dbles(1)
            hz%bpos   = dbles(2)
            hz%bneg   = dbles(3)
            hz%scale  = dbles(4)

         elseif (potcon%ipotcn.ge.1 .and. potcon%ipotcn.le.4) then

         ! Non-Hertzian options for potential contact:

            call readline(linp, ncase, linenr, 'properties of potential contact area',                  &
                          'iidddd', ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
            potcon%mx = ints(1)
            potcon%my = ints(2)
            if (potcon%ipotcn.eq.1 .or. potcon%ipotcn.eq.2) then
               potcon%xl = dbles(1)
               potcon%yl = dbles(2)
            elseif (potcon%ipotcn.eq.3 .or. potcon%ipotcn.eq.4) then
               potcon%xc1 = dbles(1)
               potcon%yc1 = dbles(2)
            endif
            if (potcon%ipotcn.eq.1 .or. potcon%ipotcn.eq.3) then
               potcon%dx = dbles(3)
               potcon%dy = dbles(4)
            elseif (potcon%ipotcn.eq.2) then
               potcon%xh = dbles(3)
               potcon%yh = dbles(4)
            elseif (potcon%ipotcn.eq.4) then
               potcon%xcm = dbles(3)
               potcon%ycm = dbles(4)
            endif
         endif

         potcon%npot = potcon%mx * potcon%my

         zerror = zerror .or. .not.check_2rng('IPOTCN', potcon%ipotcn, -6, -1, 1, 4)

         if (ic%pvtime.ne.2 .and. (potcon%mx.ne.mxold .or. potcon%my.ne.myold)) then
            zerror = .true.
            write(lout, 3001) ic%pvtime, mxold, myold, potcon%mx, potcon%my
            write(   *, 3001) ic%pvtime, mxold, myold, potcon%mx, potcon%my
 3001       format (' Input: ERROR in the input for a sequence, the number of elements MX, MY',/,      &
                    '              may not change.  P=',i3,', old:',2i5,', new:',2i5,'.')
         endif

         if (ic%iestim.ne.0 .and. (potcon%mx.ne.mxold .or. potcon%my.ne.myold)) then
            zerror = .true.
            write(lout, 3002) mxold, myold, potcon%mx, potcon%my, ic%iestim
            write(   *, 3002) mxold, myold, potcon%mx, potcon%my, ic%iestim
 3002       format (' Input: ERROR when the number of elements MX, MY changes (old:',2i5,',',/,        &
                    '              new:',2i5,'), the initial estimate cannot be used.  I=',i3,'.')
         endif

         z = potcon%mx.le.0 .or. potcon%my.le.0 .or. potcon%npot.gt.NPO
         if (z) then
            zerror = .true.
            write(lout, 3003) potcon%mx, potcon%my, NPO
            write(   *, 3003) potcon%mx, potcon%my, NPO
 3003       format (' Input: ERROR in the input for Dis, pot.con. too large',/,                        &
                    '              Mx * My > NPO : ',2i5,i8)
         endif

         if (potcon%ipotcn.eq.-6) then
            zerror = zerror .or. .not.check_range ('BPOS', hz%bpos, 1d-8, 1d20)
            zerror = zerror .or. .not.check_range ('BNEG', hz%bneg, 1d-8, 1d20)
         endif

         if (potcon%ipotcn.eq.1 .or. potcon%ipotcn.eq.3) then
            zerror = zerror .or. .not.check_range ('DX', potcon%dx, 1d-8, 1d20)
            zerror = zerror .or. .not.check_range ('DY', potcon%dy, 1d-8, 1d20)
         endif

         if (potcon%ipotcn.eq.2) then
            zerror = zerror .or. .not.check_smaller ('XL', potcon%xl , 'XH', potcon%xh)
            zerror = zerror .or. .not.check_smaller ('YL', potcon%yl , 'YH', potcon%yh)
         endif

         if (potcon%ipotcn.eq.4) then
            zerror = zerror .or. .not.check_smaller ('XC1', potcon%xc1, 'XCM', potcon%xcm)
            zerror = zerror .or. .not.check_smaller ('YC1', potcon%yc1, 'YCM', potcon%ycm)
         endif

      endif ! ic%discns_inp.eq.2

      ! combinations of ipotcn and ibound

      if (potcon%ipotcn.ge.1 .and. ic%bound.ge.2 .and. ic%bound.le.3) then
         zerror = .true.
         write(lout, 3012) ic%bound, potcon%ipotcn
         write(   *, 3012) ic%bound, potcon%ipotcn
 3012    format (' Input: ERROR. Option BOUND=',i2,' requires that a Hertzian option for the geometry', &
                 ' is used (IPOTCN=',i3,').')
      endif

      if (potcon%ipotcn.eq.-6 .and. ic%bound.ne.4) then
         zerror = .true.
         write(lout, 3013) ic%bound
         write(   *, 3013) ic%bound
 3013    format (' Input: ERROR. If IPOTCN=-6 is used, the SDEC option for the traction bound',/,       &
                 '               is required (BOUND=',i3,').')
      endif

      if (ic%bound.eq.4 .and. potcon%ipotcn.ne.-6) then
         zerror = .true.
         write(lout, 3014) potcon%ipotcn
         write(   *, 3014) potcon%ipotcn
 3014    format (' Input: ERROR. If BOUND=4 is used, the SDEC option for the geometry is ',/,           &
                 '               required (IPOTCN=',i3,').')
      endif

      ! checks on ipotcn

      if (potcon%ipotcn.ge.-5 .and. potcon%ipotcn.le.-4 .and. ic%norm.ne.1) then
         zerror = .true.
         write(lout, 3021) potcon%ipotcn
         write(   *, 3021) potcon%ipotcn
 3021    format (' Input: ERROR. Rectangular (2D) Hertzian contacts (IPOTCN=',i3,') require N=1.')
      endif

      if (potcon%ipotcn.ge.-5 .and. potcon%ipotcn.le.-4 .and. ic%norm.eq.1 .and.                        &
                                                             ic%bound.ge.2 .and. ic%bound.le.3) then
         write(lout, 3022) potcon%ipotcn, ic%bound
         write(   *, 3022) potcon%ipotcn, ic%bound
 3022    format (' Input: WARNING. Rectangular (2D) Hertzian contacts (IPOTCN=',i3,                     &
                ') do not compute the approach when B=',i1,'.')
      endif

      if (potcon%ipotcn.ge.-5 .and. potcon%ipotcn.le.-1 .and. ic%norm.eq.0 .and. kin%pen.lt.0d0) then
         zerror = .true.
         write(lout, 3023) potcon%ipotcn, kin%pen
         write(   *, 3023) potcon%ipotcn, kin%pen
 3023    format (' Input: ERROR. When using a Hertzian option for the geometry (IPOTCN=',i3,'),',/,     &
                 '               PEN must be positive (value=',g12.4,').')
      endif

      ! checks on ibound

      if (ic%bound.eq.3 .and. ic%mater.le.1) then
         write(lout, 3033) ic%bound, ic%mater
         write(   *, 3033) ic%bound, ic%mater
 3033    format (' Input: WARNING. It makes no sense to use a parabolic traction bound (BOUND=',i2,')',/, &
                 '                 when not using the simplified theory (MATER=',i2,').')
      endif

      if (ic%bound.eq.4 .and. ic%norm.ne.1) then
         zerror = .true.
         write(lout, 3035)
         write(   *, 3035)
 3035    format (' Input: ERROR. When using the SDEC approach (BOUND=4), the normal force',/,           &
                 '               must be given (needs NORM=1).')
      endif

      if (ic%bound.ge.5 .and. ic%bound.le.6 .and. ic%norm.ne.0) then
         zerror = .true.
         write(lout, 3036) ic%bound
         write(   *, 3036) ic%bound
 3036    format (' Input: ERROR. When using KPEC or ANALYN (BOUND=',i2,'), the normal force',/,         &
                 '               may not be given (needs NORM=0).')
      endif

      ! update the dimensions mx,my of the contact grid, for use in grid-functions like exrhs below

      call grid_set_dimens(cgrid, potcon%mx, potcon%my)

      !------------------------------------------------------------------------------------------------------
      ! read parameters of the undeformed distance
      !------------------------------------------------------------------------------------------------------

      if (ic%rznorm.eq.2 .and. potcon%ipotcn.lt.0) then

         ! Hertzian option: using quadratic undeformed distance,
         ! coeff. prmudf(1) and prmudf(3) determined from Hertz solution

         geom%ibase = 1
         geom%iplan = 1
         geom%prmudf(1:6) = 0d0

      elseif (ic%rznorm.eq.2) then

         ! non-Hertzian approach, undeformed distance specified by user:

         call readline(linp, ncase, linenr, 'type of geometry description, ibase and iplan', 'ii',      &
                       ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
         geom%ibase = ints(1)
         geom%iplan = ints(2)

         ! read additional parameters for undeformed distance
         !   ibase = 1: quadratic
         !   ibase = 2: circular in X
         !   ibase = 3: Hertz + Asperities in x-direction: difference
         !              of two sines.
         !   ibase = 9: pointwise given in centers of all elements

         if (geom%ibase.eq.1) then

            call read1darr(linp, ncase, linenr, 'geometry description #2', 'd', 6, idum, geom%prmudf,   &
                           idebug, lstop, ierror)

         elseif (geom%ibase.eq.2) then

            call readline(linp, ncase, linenr, 'geometry description #3', 'idddd',                      &
                          ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)

            geom%nn = ints(1)
            call reallocate_arr(geom%prmudf, 5+geom%nn)
            do k = 2, 5
               geom%prmudf(k) = dbles(k-1)
            enddo

            call read1darr(linp, ncase, linenr, 'geometry description #4', 'd', geom%nn, idum,          &
                           geom%prmudf(6:), idebug, lstop, ierror)

         elseif (geom%ibase.eq.3) then

            call read1darr(linp, ncase, linenr, 'geometry description #5', 'd', 8, idum, geom%prmudf,   &
                           idebug, lstop, ierror)

         elseif (geom%ibase.eq.9) then

            npot = potcon%mx * potcon%my
            call reallocate_arr(geom%prmudf, npot+10)

            call read1darr(linp, ncase, linenr, 'undeformed distance in all elements', 'd',             &
                           npot, idum, geom%prmudf, idebug, lstop, ierror)
         endif

         ! read additional parameters for planform
         !   iplan = 1: unrestricted

         if (geom%iplan.eq.2) then

            ! iplan = 2: elliptical planform (quadratic function)
            call read1darr(linp, ncase, linenr, 'planform description #2', 'd', 6, idum, geom%prmpln,   &
                           idebug, lstop, ierror)

         elseif (geom%iplan.eq.3) then

            ! iplan = 3: union of two rectangles

            call read1darr(linp, ncase, linenr, 'planform description #3', 'd', 8, idum, geom%prmpln,   &
                           idebug, lstop, ierror)

         elseif (geom%iplan.eq.4) then

            ! iplan = 4: weighted interaction between patches

            call readline(linp, ncase, linenr, 'planform description: npatch', 'i',                     &
                          ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)

            associate( np => geom%npatch )
            np = max(1, min(100, ints(1)))

            allocate(tmp(np*np))
            call reallocate_arr(geom%ysep, np )
            call reallocate_arr(geom%facsep, np, np)
            geom%ysep(1:np)         = 0d0
            geom%facsep(1:np, 1:np) = 0d0

            call read1darr(linp, ncase, linenr, 'planform description: ysep', 'd', np-1, idum,          &
                           geom%ysep, idebug, lstop, ierror)

            call read1darr(linp, ncase, linenr, 'planform description: fac', 'd', np*(np-1)/2, idum,    &
                           tmp, idebug, lstop, ierror)

            ! expand compressed list obtained in tmp to full matrix in facsep

            if (idebug.ge.4) then
               spaces = ' '
               iof = 0
               do ip = 1, np
                  write(bufout,'(a,f7.3,2a,7f7.3)') ' ysep=', geom%ysep(ip),', fac=', spaces(1:7*ip),   &
                           (tmp(iof+j), j=1, np-ip)
                  call write_log(1, bufout)
                  iof = iof + np-ip
               enddo
            endif

            k = 0
            do ip = 1, np
               geom%facsep(ip,ip)   = 1d0       ! diagonal entry
               do j = ip+1, np
                  k = k + 1
                  geom%facsep(ip,j) = tmp(k)    ! upper triangular part
                  geom%facsep(j,ip) = tmp(k)    ! symmetry: lower triangular part
               enddo
            enddo

            if (idebug.ge.4) then
               write(bufout,'(a,i3)') ' IPLAN=4: npatch=',np
               call write_log(1, bufout)
               do ip = 1, np
                  write(bufout,'(a,f7.3,a,7f7.3)') ' ysep=', geom%ysep(ip),', fac=',                    &
                        (geom%facsep(ip,j), j=1, np)
                  call write_log(1, bufout)
               enddo
            endif
            deallocate(tmp)
            end associate

         endif

         ! perform checks on input for undeformed distance

         zerror = zerror .or. .not.check_2rng('IBASE', geom%ibase, 1, 3, 9, 9)
         zerror = zerror .or. .not.check_irng ('IPLAN', geom%iplan, 1, 4)

         if (geom%ibase.eq.2) then
            zerror = zerror .or. .not.check_irng ('NN', geom%nn, 2, 999999)
            zerror = zerror .or. .not.check_range ('DY1', geom%prmudf(5), 1d-8, 1d20)
         endif

         if (geom%iplan.eq.3 .and. (geom%prmpln(2).le.geom%prmpln(1) .or.                              &
             geom%prmpln(4).le.geom%prmpln(3) .or. geom%prmpln(6).le.geom%prmpln(5) .or.               &
             geom%prmpln(8).le.geom%prmpln(7))) then
            zerror = .true.
            write(lout, 5004) (geom%prmpln(k), k=1,8)
            write(   *, 5004) (geom%prmpln(k), k=1,8)
 5004       format (' Input: ERROR in input for Planform. Check the relations:',/,                     &
                    '    Xl1 =',g12.4,' <= ',g12.4,' = Xh1,',/                                         &
                    '    Yl1 =',g12.4,' <= ',g12.4,' = Yh1,',/                                         &
                    '    Xl2 =',g12.4,' <= ',g12.4,' = Xh2,',/                                         &
                    '    Yl2 =',g12.4,' <= ',g12.4,' = Yh2.')
         endif

         if (geom%iplan.eq.4) then
            zerror = zerror .or. .not.check_sorted( 'YSEP', geom%npatch-1, geom%prmpln(1:), .true. )
            call write_log(nline_errmsg, errmsg)
         endif

      ! endif (Z=2, undeformed distance)

      endif

      !------------------------------------------------------------------------------------------------------
      ! read & check input for V- and L-digits for V = 2, parameters per row, NVF = MY
      !------------------------------------------------------------------------------------------------------

      if (ic%frclaw_inp.ne.1 .and. ic%varfrc.eq.2) then
         call fric_input(linp, ncase, linenr, ic%varfrc, ic%frclaw_inp, potcon%my, fric, idebug,        &
                ieof, lstop, zerror)
      endif

      ! Adapt the value of the friction coefficient used for scaling of tangential forces

      kin%use_muscal = ic%varfrc.eq.0
      if (kin%use_muscal) then
         kin%muscal = fric%fstat()
      else
         kin%muscal = 1d0
      endif

      !------------------------------------------------------------------------------------------------------
      ! read the extra term in the tangential undeformed distance
      !------------------------------------------------------------------------------------------------------

      if (ic%rztang.eq.9) then
         npot = potcon%mx * potcon%my

         ! re-allocate exrhs at the appropriate size

         call gf3_new(geom%exrhs, 'geom%exrhs', cgrid)

         allocate(tmp(2*npot))
         call read1darr(linp, ncase, linenr, 'extra term in rigid slip of elements', 'd', 2*npot,       &
                           idum, tmp, idebug, lstop, ierror)
         do ii = 1, npot
            geom%exrhs%vx(ii) = tmp(2*ii-1)
            geom%exrhs%vy(ii) = tmp(2*ii  )
         enddo
         deallocate(tmp)
      endif

      !------------------------------------------------------------------------------------------------------
      ! read subsurface points
      !------------------------------------------------------------------------------------------------------

      if (ic%stress.ge.2) then
         call rdsubs(linp, ic, ncase, linenr, idebug, subs)
      endif

      !------------------------------------------------------------------------------------------------------
      ! Stop program if an error was found, except when inp=3, when the input-file must be tested.
      !------------------------------------------------------------------------------------------------------

      if (zerror) then
         write(lout, 8002) ncase, line0
         write(   *, 8002) ncase, line0
 8002    format (' Errors found in the input.'/, ' Case-number',i8,', starting at line',i8)
         if (inp.ne.3) call abort_run()
      endif

      if (idebug.ge.1) write(*,'(a)') '--- end subroutine input.f ---'

   end subroutine input

!------------------------------------------------------------------------------------------------------------

   subroutine wrtinp(ncase, ic, potcon, mater, hz, geom, fric, kin, solv, outpt1, subs, ltight)
!--purpose: Write input for current case to INPUT-file <experim>.inp, unit linp.
      implicit none
!--subroutine parameters:
      integer          :: ncase
      type(t_ic)       :: ic
      type(t_potcon)   :: potcon
      type(t_hertz)    :: hz
      type(t_geomet)   :: geom
      type(t_material) :: mater
      type(t_friclaw)  :: fric
      type(t_kincns)   :: kin
      type(t_solvers)  :: solv
      type(t_output)   :: outpt1
      type(t_subsurf)  :: subs
      logical          :: ltight  ! flag: use solution to create tight pot.con.
!--local variables:
      integer      :: vldcmze, xgiaowr, pbtnfs
      integer      :: ii, ix, ix0, ix1, iy, iy0, iy1, iof, k, len1
      logical      :: is_roll
      real(kind=8) :: fun, fux, fuy
      character(len=80)  :: spaces

      ! Convert control integers to packed form

      call ic_pack (3, pbtnfs, vldcmze, xgiaowr, ic)
      is_roll = ic%tang.eq.2 .or. ic%tang.eq.3

      ! If tight potential contact desired: determine box around actual contact area (requires solution!)

      if (ltight) then
         call areas(outpt1%igs)
         if (outpt1%igs%ixmin.gt.outpt1%igs%ixmax) then
            ! empty contact area
            ix0 = max(        1, (potcon%mx+1)/2-2)
            ix1 = min(potcon%mx, (potcon%mx+1)/2+2)
            iy0 = max(        1, (potcon%my+1)/2-2)
            iy1 = min(potcon%my, (potcon%my+1)/2+2)
         else
            ! non-empty contact area, clipped
            ix0 = max(        1, outpt1%igs%ixmin-2)
            ix1 = min(potcon%mx, outpt1%igs%ixmax+2)
            iy0 = max(        1, outpt1%igs%iymin-2)
            iy1 = min(potcon%my, outpt1%igs%iymax+2)
         endif
      else
         ! no clipping, full potential contact area
         ix0 = 1
         ix1 = potcon%mx
         iy0 = 1
         iy1 = potcon%my
      endif

      ! Determine kinematic constants fun, fux, fuy

      if (ic%norm.eq.0) then
         fun = kin%pen
      else
         fun = kin%fntrue
      endif

      if (ic%force.eq.0) then
         fux = kin%cksi
         fuy = kin%ceta
      elseif (ic%force.eq.1) then
         fux = kin%fxrel1
         fuy = kin%ceta
      elseif (ic%force.eq.2) then
         fux = kin%fxrel1
         fuy = kin%fyrel1
      endif

      ! write flow-control variables and kinematic constants

      if (ncase.gt.0) write(linp,'(a,i8)') '% Next case', ncase
      write(linp, 1101) pbtnfs, vldcmze, xgiaowr
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
      if (ic%norm.eq.0 .and. ic%force.eq.0) write(linp, 2131) fun, fux, fuy, kin%cphi
      if (ic%norm.eq.0 .and. ic%force.eq.1) write(linp, 2132) fun, fux, fuy, kin%cphi
      if (ic%norm.eq.0 .and. ic%force.eq.2) write(linp, 2133) fun, fux, fuy, kin%cphi
      if (ic%norm.eq.1 .and. ic%force.eq.0) write(linp, 2134) fun, fux, fuy, kin%cphi
      if (ic%norm.eq.1 .and. ic%force.eq.1) write(linp, 2135) fun, fux, fuy, kin%cphi
      if (ic%norm.eq.1 .and. ic%force.eq.2) write(linp, 2136) fun, fux, fuy, kin%cphi
      if (ic%varfrc.eq.0) call fric_wrtinp(linp, ic%varfrc, ic%frclaw_inp, fric)

 1101 format (i8.6, 6x,   'P-B-T-N-F-S          PVTIME, BOUND , TANG , NORM , FORCE, STRESS', /,        &
              i8.6, 6x,   'L-D-C-M-Z-E          FRCLAW, DISCNS, INFLCF, MATER, RZNORM, EXRHS ', /,      &
              i8.7, 4x, 'H-G-I-A-O-W-R    HEAT, GAUSEI, IESTIM, MATFIL, OUTPUT, FLOW, RETURN' )
 1201 format( 2i6, 2i5, 3x, es8.1, 8x,   'MAXGS , MAXIN , MAXNR , MAXOUT, EPS')
 1211 format( 2g12.4, i6, g12.4, 7x, 'OMEGAH, OMEGAS, INISLP, OMGSLP')
 1212 format( i6, g12.4, 40x, 'INISLP, OMGSLP')
 1215 format( 6f8.3,  10x, 'FDECAY, D_IFC/LIN/CNS, D_SLP, POW_S')
 2131 format( 4g12.4, 10x, 'PEN, CKSI, CETA, CPHI')
 2132 format( 4g12.4, 10x, 'PEN, FX, CETA, CPHI')
 2133 format( 4g12.4, 10x, 'PEN, FX, FY, CPHI')
 2134 format( 4g12.4, 10x, 'FN, CKSI, CETA, CPHI')
 2135 format( 4g12.4, 10x, 'FN, FX, CETA, CPHI')
 2136 format( 4g12.4, 10x, 'FN, FX, FY, CPHI')

      ! if new influence coefficients: write rolling step, material constants when C=2, 3, 4 or 9

      if ((ic%gencr_inp.ge.2 .and. ic%gencr_inp.le.4) .or. ic%gencr_inp.eq.9) then

         ! write the rolling direction and the shift

         if (is_roll) then
            write(linp, 5001) kin%chi, kin%dq, kin%veloc
 5001       format (3g12.4, 22x, 'CHI, DQ, VELOC')
         endif

         ! write surface inclinations for the Blanco-approach

         if (ic%gencr_inp.eq.4) then
            write(linp, 5201) mater%if_meth, mater%if_ver, mater%ninclin
            do ii = 1, mater%ninclin
               write(linp, 5202) mater%surf_inclin(ii,1), mater%surf_inclin(ii,2)
            enddo
 5201       format (3i6, 40x, 'IF_METH, VARIANT, NN')
 5202       format (f8.4,4x,f8.4)
         endif

         ! filename for numerical influence coefficients

         if (ic%gencr_inp.eq.9) then
            spaces = ' '
            len1 = len(trim(mater%fname_influe))
            if (len1.le.52) then
               write(linp, 5103) '''', trim(mater%fname_influe), '''', spaces(1:53-len1)
            else
               write(linp, 5103) '''', trim(mater%fname_influe), '''', '  '
            endif
 5103       format(3x, 4a, 'CFNAME')
         endif

         ! write elastic material constants

         write(linp, 5300) mater%poiss, mater%gg
 5300    format (2(f8.4, 4x), 2g12.4, 10x, 'POISS 1,2,  GG 1,2')

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
 5501       format (4g12.4, 10x, 'GG3,LAYTHK,TAU_C0,K_TAU')
         endif
      endif

      ! write material parameters for temperature model

      if (ic%heat.eq.3) then
         write(linp,5611) mater%bktemp(1), mater%heatcp(1), mater%lambda(1), mater%dens(1), 1,1,1,1
         write(linp,5611) mater%bktemp(2), mater%heatcp(2), mater%lambda(2), mater%dens(2), 2,2,2,2
 5611    format (4g12.4, 10x, 'BKTEMP',i1,', HEATCP',i1,', LAMBDA',i1,', DENS',i1)
      endif

      if (ic%heat.eq.3 .and. ic%mater.eq.4) then
         write(linp,5612) mater%betapl
         write(linp,5612) mater%betapl
 5612    format (1g12.4, 46x, 'BETAPL')
      endif

      ! write discretisation parameters

      if (ic%discns_inp.eq.2) then
         associate(dx => potcon%dx, dy => potcon%dy)
         write(linp, 6300)  potcon%ipotcn
         if (potcon%ipotcn.eq.-1) then
            write(linp, 6291) potcon%mx, potcon%my, hz%a1, hz%b1, hz%scale
         elseif (potcon%ipotcn.eq.-2) then
            write(linp, 6292) potcon%mx, potcon%my, hz%a1, hz%aob, hz%scale
         elseif (potcon%ipotcn.eq.-3) then
            write(linp, 6293) potcon%mx, potcon%my, hz%aa, hz%bb, hz%scale
         elseif (potcon%ipotcn.eq.-4) then
            write(linp, 6294) potcon%mx, potcon%my, hz%a1, hz%bb, hz%scale
         elseif (potcon%ipotcn.eq.-5) then
            write(linp, 6293) potcon%mx, potcon%my, hz%aa, hz%bb, hz%scale
         elseif (potcon%ipotcn.eq.-6) then
            write(linp, 6296) potcon%mx, potcon%my, hz%aa, hz%bneg, hz%bpos, hz%scale
         elseif (potcon%ipotcn.eq.1) then
            write(linp, 6301) (ix1-ix0+1), (iy1-iy0+1), potcon%xl+(ix0-1)*dx, potcon%yl+(iy0-1)*dy, dx, dy
         elseif (potcon%ipotcn.eq.2) then
            write(linp, 6301) (ix1-ix0+1), (iy1-iy0+1), potcon%xl+(ix0-1)*dx, potcon%yl+(iy0-1)*dy,    &
                              potcon%xh-(potcon%mx-ix1)*dx, potcon%yh-(potcon%my-iy1)*dy
         elseif (potcon%ipotcn.eq.3) then
            write(linp, 6301) (ix1-ix0+1), (iy1-iy0+1), potcon%xc1+(ix0-1)*dx, potcon%yc1+(iy0-1)*dy, dx, dy
         elseif (potcon%ipotcn.eq.4) then
            write(linp, 6301) (ix1-ix0+1), (iy1-iy0+1), potcon%xc1+(ix0-1)*dx, potcon%yc1+(iy0-1)*dy,  &
                                           potcon%xcm-(potcon%mx-ix1)*dx, potcon%ycm-(potcon%my-iy1)*dy
         endif
         end associate
 6300    format ( i5, 53x   , 'IPOTCN')
 6291    format (2i5, 3g12.4, 12x, 'MX,MY,A1,B1,SCALE')
 6292    format (2i5, 3g12.4, 12x, 'MX,MY,A1,AOB,SCALE')
 6293    format (2i5, 3g12.4, 12x, 'MX,MY,AA,BB,SCALE')
 6294    format (2i5, 3g12.4, 12x, 'MX,MY,A1,BB,SCALE')
 6296    format (2i5, 4g12.4, 'MX/Y,AA,BNEG/POS,SCALE')
 6301    format (2i5, 4g12.4, 'MX,MY,XL,YL,DX,DY')
 6302    format (2i5, 4g12.4, 'MX,MY,XL,YL,XH,YH')
 6303    format (2i5, 4g12.4, 'MX,MY,XC1,YC1,DX,DY')
 6304    format (2i5, 4g12.4, 'MX,MY,XC1,YC1,XCM,YCM')
      endif

      ! write the parameters determining the undeformed distance

      if (ic%rznorm.eq.2 .and. potcon%ipotcn.gt.0) then
         write(linp, 7500) geom%ibase, geom%iplan
 7500    format ( 2i5, 48x, 'IBASE, IPLAN')
         if (geom%ibase.eq.1) then
            write(linp,'(a)') '%  UNDEF.DISTANCE, (1)-(2): QUADRATIC'
            write(linp, 7501) (geom%prmudf(k), k= 1, 6)
 7501       format (5g12.4,/,g12.4,46x,'B(k), k=1, 6')
         elseif (geom%ibase.eq.2) then
            write(linp,'(2a)') '%  UNDEF.DISTANCE, (1)-(2): CIRCULAR (X), POINTWISE (Y)'
            write(linp, 7502) geom%nn, (geom%prmudf(k), k=2,5)
 7502       format (i5,3x, 4g12.4,2x,'NN,XM,RM,Y1,DY1')
            write(linp, 7505) (geom%prmudf(5+k), k=1, geom%nn)
 7505       format ( 1000( 20(5g12.4,:,/),:/) )
         elseif (geom%ibase.eq.3) then
            write(linp,'(2a)') '%  UNDEF.DISTANCE, (1)-(2): QUADRATIC PLUS DIFFERENCE OF TWO SINES'
            write(linp, 7503) (geom%prmudf(k), k=1,8)
 7503       format( 5g12.4,/,3g12.4, 22x,'B(k), k=1, 8')
         elseif (geom%ibase.eq.9) then
            write(linp,'(2a)') '%  UNDEF.DISTANCE, (1)-(2): SPECIFIED PER ELEMENT'
            do iy=iy0, iy1
               iof = (iy-1)*potcon%mx
               if (iy.gt.1) write(linp,'(1x)')
               write(linp, 7509) (geom%prmudf(iof+ix), ix=ix0,ix1)
 7509          format ( 5g14.6,:,/, 9(2x,5g14.6,:,/), 10( /, 10(2x,5g14.6,:,/) ) )
            enddo
         endif

         if (geom%iplan.eq.1) then
            write(linp,'(a)') '%  UNRESTRICTED PLANFORM'
         elseif (geom%iplan.eq.2) then
            write(linp,'(a)') '%  QUADRATIC PLANFORM'
            write(linp, 7522) (geom%prmpln(k), k=1, 6)
 7522       format (5g12.4,/,g12.4,46x,'P(k), k=1, 6')
         else
            write(linp,'(a)') '%  PLANFORM = UNION OF RECTANGLES'
            write(linp, 7523) (geom%prmpln(k), k=1, 8)
 7523       format (4g12.4,/,4g12.4,10x,'P(k), k=1, 8')
         endif
      endif

      ! write the friction input

      if (ic%varfrc.eq.2) call fric_wrtinp(linp, ic%varfrc, ic%frclaw_inp, fric)

      ! write the extra term in the undeformed distance for the tangential problem

      if (ic%rztang.eq.9) then
         write(linp, '(a)') '%  EXRHS(I,2), EXRHS(I,3)'
         do iy=iy0, iy1
            do ix=ix0, ix1
               ii = ix + (iy-1)*potcon%mx
               write(linp, 8601) geom%exrhs%vx(ii), geom%exrhs%vy(ii)
            enddo
         enddo
 8601    format (100 (100 (20 (6g12.4/4g12.4/))))
      endif

      ! write the points in which the subsurface elastic field must be calculated. 

      if (ic%stress.ge.2) then
         call wrsubs(linp, ic, subs)
      endif

      ! write empty line to separate from next case

      write(linp, '(1x)')
   end subroutine wrtinp

!------------------------------------------------------------------------------------------------------------ 

end module m_sinput
