!------------------------------------------------------------------------------------------------------------
! m_scontc - driver routine for solving a basic contact problem (module 3)
!
! Copyright 1979-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_scontc

use m_hierarch_data
use m_hertz
use m_sdis
use m_visc
use m_snorm
use m_stang
use m_soutpt
use m_subsurf
use m_temperature

#ifdef _OPENMP
use omp_lib      , only : omp_get_thread_num
#endif


implicit none
private

public  contac
private check_case
public  panprc
private set_dummy_solution

contains

!------------------------------------------------------------------------------------------------------------

   subroutine contac (gd, ierror)
!--purpose: Performs calculations for one case. Calls dis, rznorm, rztang, sgencr if this is requested,
!           and if R=0,1 computes tractions ps and contact area igs.
      implicit none
!--subroutine arguments:
      type(t_probdata) :: gd
      integer          :: ierror
!--local variables :
      integer, parameter :: idebug = 0
      integer  :: ihertz, nerror
      logical  :: is_fastsm, is_roll, is_ssrol

      associate( potcon_inp => gd%potcon_inp, potcon_cur => gd%potcon_cur, cgrid      => gd%cgrid_cur,  &
                 meta   => gd%meta,   ic     => gd%ic,     mater  => gd%mater,  hertz  => gd%hertz,     &
                 geom   => gd%geom,   fric   => gd%fric,   kin    => gd%kin,    influ  => gd%influ,     &
                 solv   => gd%solv,   outpt1 => gd%outpt1, subs   => gd%subs)

      is_fastsm = ic%mater.eq.2 .or. ic%mater.eq.3
      is_roll   = ic%tang.eq.2 .or. ic%tang.eq.3
      is_ssrol  = ic%tang.eq.3

      ! check the input for the case

      nerror = 0
      ierror = 0
      if (idebug.ge.5) call write_log('calling check_case...')
      call check_case(gd, ierror, nerror)

      if (nerror.gt.0) then
         call set_dummy_solution(gd)
         return
      endif

      ! compute combined material constants

      call combin_mater(mater)

      ! (re)initialize penv when P=2

      if (ic%pvtime.eq.2) kin%penv = 0d0

      ! solve Hertzian problem when needed, complete potcon_inp

      if (potcon_inp%ipotcn.ge.-5 .and. potcon_inp%ipotcn.le.-1) then
         if (idebug.ge.5) call write_log(' contac: calling hzsol')
         call hzsol(ic, geom, potcon_inp%ipotcn, hertz, kin, mater)
         call potcon_hertz(hertz, potcon_inp)
      endif

      ! M=3: calculate flexibilities L1,L2,L3 from Hertzian contact ellipse

      if (ic%mater.eq.3 .and. potcon_inp%ipotcn.le.-1) then
         if (idebug.ge.5) call write_log(' contac: calling simpflex')
         call simpflex(ic, hertz%aa, hertz%bb, mater, fric%fstat(), kin, 0)
      endif

      ! prepare potcon and grid for current time, resize grid-functions as needed

      if (ic%discns3.ge.1) then
         if (idebug.ge.5) call write_log(' contac: calling init_curr_grid')
         call init_curr_grid(ic, potcon_inp, potcon_cur, cgrid, geom, outpt1)
      endif

      ! check/enforce requirements on rolling direction/step size

      call check_roll_stepsize(ic, solv, kin, cgrid)

      ! cycle grid-functions from current to previous time

      if (idebug.ge.5) call write_log(' contac: calling set_prev_data')
      call set_prev_data(ic, geom, mater, fric, kin, outpt1)

      ! calculate influence functions cs, csv
      !  - compute standard infl.cf when C>=1
      !  - compute Blanco-correction when C=4
      !  - load numerical influence coefficients when C=9
      !     --> overwrite cs,csv when format=0, fill cy when format=1

      if (ic%gencr_inp.ge.1) then
         if (idebug.ge.5) call write_log(' contac: calling sgencr')
         call sgencr (ic, cgrid, mater, influ, fric, kin)
      endif

      if (ic%gencr_inp.eq.4 .or. (ic%gencr_inp.eq.1 .and. mater%gencr_eff.eq.4)) then
         if (idebug.ge.5) call write_log(' contac: calling influe_blanco')
         call influe_blanco (is_roll, mater%if_meth, mater%if_ver, mater%ninclin, mater%surf_inclin,    &
                             cgrid, influ)
      endif

      if (ic%gencr_inp.eq.9 .or. (ic%gencr_inp.eq.1 .and. mater%gencr_eff.eq.9)) then
         if (idebug.ge.5) call write_log(' contac: calling influe_load')

         call influe_load(mater%fname_influe, meta%dirnam, is_roll, cgrid, influ)
         call influe_mater(influ, mater)

         if (ic%x_inflcf.eq.6) then
            write(lout,*) 'Printing influence coefficients cs (instat)'
            call inflcf_print (influ%cs, lout)
         endif
      endif

      ! form undeformed distance hs1

      if (ic%rznorm.gt.0) then
         if (idebug.ge.5) call write_log(' contac: calling set_norm_rhs')
         call set_norm_rhs(cgrid, geom)
      endif

      ! compute the constant part in the rigid slip, i.e. excluding creepages that are unknown
      ! note for Fastsim + 3 flexibilities + non-Hertzian: flexibilities not yet filled in

      if (ic%mater.ne.3 .or. potcon_cur%ipotcn.lt.0) then
         if (idebug.ge.5) call write_log(' contac: calling set_tang_rhs')
         call set_tang_rhs(ic, mater, cgrid, kin, geom)
      endif

      ! set initial traction Ps and previous traction Pv :

      ihertz = 0
      if (potcon_cur%ipotcn.lt.0) ihertz = abs(potcon_cur%ipotcn)

      if (idebug.ge.5) call write_log(' contac: calling init_curr_data')
      call init_curr_data (ic, hertz, potcon_cur, cgrid, geom, mater, kin, outpt1, ihertz)

      ! skip computing part?

      if (ic%return.le.1) then

         ! the computing part :

         ! Start the actual solution: Panagiotopoulos process

         ! call write_log('ps-z:')
         ! call print_fld(outpt1%ps, ikZDIR)

         if (idebug.ge.5) call write_log(' contac: calling panprc')
         call panprc (ic, ihertz, mater, potcon_cur, cgrid, influ, fric, kin, solv, outpt1, geom)

         ! TODO: fill outputs uv, sv, uplv, taucv, etc. in steady rolling problems

         ! Compute surface temperatures

         if (ic%heat.ge.1 .and. fric%frclaw_eff.ne.6) then
            if (idebug.ge.5) call write_log(' contac: calling calc_temp')
            call calc_temp (ic, mater, kin, cgrid, outpt1)
         endif

         ! Derive quantities like total forces and elastic energy, write output to out-file when O>=1.

         if (idebug.ge.5) call write_log(' contac: calling soutpt')
         call soutpt (gd)

         ! Output results to .mat-file if requested

         if (ic%matfil_surf.ge.1) then
            if (idebug.ge.5) call write_log(' contac: calling writmt')
            call timer_start(itimer_files)
            call writmt (meta, ic, cgrid, potcon_cur, geom%hs1, mater, fric, kin, outpt1, .false.)
            call timer_stop(itimer_files)
         endif

         ! Compute the subsurface elastic field for points specified in Xb, write output to out-file
         ! and subs-file

         if (ic%stress.ge.1) then
            if (idebug.ge.5) call write_log(' contac: calling subsur')
            call subsur(meta, ic, mater, cgrid, outpt1%igs, outpt1%ps, .false., subs)
         endif

         if (idebug.ge.5) call write_log(' contac: done, returning')

      endif ! RETURN <= 1

      end associate
   end subroutine contac

!------------------------------------------------------------------------------------------------------------

   subroutine print_fld(fld, ikarg)
      implicit none
!--subroutine arguments:
      type(t_gridfnc3) :: fld
      integer          :: ikarg
!--local variables:
      integer      :: ix, iy, ifac
      real(kind=8) :: absmax, fac

      absmax = gf3_maxabs(AllElm, fld, ikarg)
      if (absmax.gt.1d-10) then
         ifac   = int( log10( absmax )) - 1
         fac    = 10d0**ifac
      else
         ifac   = 0
         fac    = 1d0
      endif
      write(bufout,'(a,g12.4,a,i3,a,es9.1)') 'absmax=',absmax,', ifac=',ifac,', fac=',fac
      call write_log(1, bufout)

      associate(mx => fld%grid%nx, my => fld%grid%ny)
      do iy = my, 1, -1
         write(bufout,'(i3,a,20f7.3)') iy,':', (fld%val(ix+(iy-1)*mx,ikarg)/fac, ix=1,min(20,mx))
         call write_log(1, bufout)
      enddo
      end associate
   end subroutine print_fld

!------------------------------------------------------------------------------------------------------------

   subroutine check_case (gd, ierror, nerror)
!--purpose: Perform checks on the input for one case.
      implicit none
!--subroutine arguments:
      type(t_probdata)     :: gd
      integer, intent(out) :: ierror    ! first error found
      integer, intent(out) :: nerror    ! number of errors found
!--local variables :
      integer, parameter :: idebug = 0
      logical      :: zerror
      integer      :: kofs_x, kofs_y
      real(kind=8) :: delt_x, delt_y

      nerror = 0
      ierror = 0

      associate( ic => gd%ic, meta => gd%meta, geom => gd%geom, cgrid => gd%cgrid_inp, fric => gd%fric )

      if (meta%ncase.le.1 .and. ic%pvtime.ne.2) then
         write(bufout, 2001) ic%iestim
         call write_log(2, bufout)
 2001    format (' WARNING. In the first case contact must be initiated.'/,                            &
                 '          Digit P=',i3,' is overruled, set to 2.')
         ic%pvtime = 2
      endif

!     if (meta%ncase.le.1 .and. ic%iestim.ne.0) then
!        write (bufout, 2004) ic%iestim
!        call write_log(2, bufout)
!2004    format (' WARNING. In the first case contact must be initiated.'/,                            &
!                '          Digit I=', i3,' is overruled, set to 0.')
!        ic%iestim = 0
!     endif

      if (ic%pvtime.ne.2 .or. ic%iestim.ne.0) then
         associate( p0 => gd%potcon_cur, p1 => gd%potcon_inp )

         if (.false. .and. (p1%mx.ne.p0%mx .or. p1%my.ne.p0%my)) then
            write(bufout, 3005) ic%pvtime, ic%iestim, p0%mx, p0%my, p1%mx, p1%my
            call write_log(2, bufout)
 3005       format (' WARNING: P=',i3,', I=',i3,': changing the number of elements MX, MY.',/,          &
                    '                 Old:',2i5,', new:',2i5,'.')
         endif

         zerror = (abs(p0%dx-p1%dx).ge.1d-6*min(p0%dx, p1%dx) .or.                                      &
                   abs(p0%dy-p1%dy).ge.1d-6*min(p0%dy, p0%dy))

         if (zerror) then
            nerror = nerror + 1
            if (ierror.eq.0) ierror = 3006
            write(bufout, 3006) ic%pvtime, ic%iestim, p0%dx, p0%dy, p1%dx, p1%dy
            call write_log(2, bufout)
 3006       format (' ERROR: P=',i3,', I=',i3,': needs fixed grid sizes DX, DY.',/,                     &
                    '          Old:',2g12.4,', new:',2g12.4,'.')
         endif

         delt_x = p1%xc1 - p0%xc1
         delt_y = p1%yc1 - p0%yc1
         kofs_x = nint( delt_x / p1%dx )
         kofs_y = nint( delt_y / p1%dy )
         zerror = (abs(delt_x-kofs_x*p1%dx).ge.1d-6*p1%dx .or. abs(delt_y-kofs_y*p1%dy).ge.1d-6*p1%dy)

         if (zerror) then
            nerror = nerror + 1
            if (ierror.eq.0) ierror = 3006
            write(bufout, 3007) ic%pvtime, ic%iestim, delt_x, delt_y, p1%dx, p1%dy
            call write_log(2, bufout)
 3007       format (' ERROR: P=',i3,', I=',i3,': grid offsets DELTX =',g12.4,', DELTY =',g12.4, /,      &
                    '        must be multiples of DX =',g12.4,', DY =',g12.4)
         endif

         end associate
      endif ! P<>2 or I>0

      if (ic%discns3.ne.0 .and. ic%gencr_inp.eq.0 .and.                                                 &
          (.not.check_equal('new DX', gd%potcon_inp%dx, 'old DX', gd%influ%cs%dx, 1d-4, .false.) .or.   &
           .not.check_equal('new DY', gd%potcon_inp%dy, 'old DY', gd%influ%cs%dy, 1d-4, .false.))) then
         nerror = nerror + 1
         if (ierror.eq.0) ierror = 3010
         write(bufout, 3010) ic%discns3, ic%gencr_inp
         call write_log(1, bufout)
 3010    format (' ERROR. New DX, DY need new influence functions, digits D, C=', 2i3,'.')
      endif

      if (ic%discns3.ne.0 .and. ic%gencr_inp.eq.0 .and.                                                 &
          (gd%potcon_inp%mx.gt.gd%influ%cs%cf_mx .or. gd%potcon_inp%my.gt.gd%influ%cs%cf_my)) then
         nerror = nerror + 1
         if (ierror.eq.0) ierror = 3011
         write(bufout, 3011) ic%discns3, ic%gencr_inp
         call write_log(1, bufout)
 3011    format (' ERROR. A larger pot.contact needs new influence functions, digits D, C=', 2i3,'.')
      endif

      if (fric%nvf.ne.1 .and. fric%nvf.ne.cgrid%ny) then
         nerror = nerror + 1
         if (ierror.eq.0) ierror = 4001
         write(bufout, 4001) fric%nvf, cgrid%ny
         call write_log(1, bufout)
 4001    format (' Internal ERROR. The number of friction values NVF =',i6,' must be 1 or MY = ', i6,'.')
      endif

      if (geom%iplan.eq.4) then
         zerror = .not.check_sorted( 'YSEP', geom%npatch-1, geom%ysep, .true. )
         if (zerror) then
            nerror = nerror + 1
            if (ierror.eq.0) ierror = 5001
            call write_log(nline_errmsg, errmsg)
         endif
      endif

      end associate
   end subroutine check_case

!------------------------------------------------------------------------------------------------------------

      subroutine panprc(ic, ihertz, mater, potcon, cgrid, infl, fric, kin, solv, outpt1, geom)
!--purpose: Implement Panagiotopoulos-process, alternatively compute normal and tangential stresses,
!           with the other stresses fixed.
      implicit none
!--subroutine arguments:
      type(t_ic)              :: ic
      type(t_material)        :: mater
      type(t_potcon)          :: potcon
      type(t_grid)            :: cgrid
      type(t_geomet)          :: geom
      type(t_influe)          :: infl
      type(t_friclaw)         :: fric
      type(t_kincns)          :: kin
      type(t_solvers)         :: solv
      type(t_output)          :: outpt1
      integer                 :: ihertz
!--local variables:
      integer,          parameter :: idebug = 0
      type(t_gridfnc3)            :: po1
      integer                     :: ii, it, itout, j, ncon
      logical                     :: is_ssrol
      real(kind=8)                :: dif, difid, a_equiv, b_equiv
      character(len=12)           :: strng(4)

      associate(npot  => cgrid%ntot,  ps1   => outpt1%ps,    igs1  => outpt1%igs, uv1   => outpt1%uv,   &
                hs1   => geom%hs1,    hv1   => geom%hv1)
      call timer_start(itimer_panprc)

      call gf3_new(po1, 'panprc:po1', cgrid, ps1%eldiv, nulify=.true.)

      is_ssrol = ic%tang.eq.3

      solv%itnorm = 0
      solv%ittang = 0

      ! Copy initial ps1 to po1

      call gf3_copy(AllElm, ps1, po1, ikALL)

      itout = 0
      difid = 1d0
      dif   = 200d0 * difid

      ! While (dif > eps) AND (inner,outer not divergent) do

      do while (dif.gt.difid .and. itout.lt.solv%maxout .and. solv%itnorm.ge.0 .and. solv%ittang.ge.0)
         itout = itout + 1

         ! Compute new normal traction

         if (ic%bound.ge.2 .and. ic%bound.le.4 .and. ihertz.ge.1) then
            if (idebug.ge.2) call write_log('Panprc: skipping normal problem - using Hertz-solution')
            it = 0
         elseif (ic%bound.ge.5 .and. ic%bound.le.6) then
            if (idebug.ge.5) call write_log('panprc: call snorm_kpec...')
            call snorm_kpec_wrapper(ic, mater, cgrid, kin, solv, hs1, infl%cs, it, igs1, ps1)
            if (idebug.ge.5) call write_log('panprc: returned from snorm_kpec...')
            it = 0
         else
            ! Using NormCG
            if (idebug.ge.5) call write_log('panprc: calling snorm...')
            call snorm(ic, mater, cgrid, kin, geom, solv, hs1, infl%cs, infl%ms, it, igs1,              &
                       ps1, outpt1%sens)
         endif
         if (it.ge.0) then
            solv%itnorm = solv%itnorm + it
         else
            ! maxin reached in NORM
            solv%itnorm = -1
         endif

         ! Fastsim with 3 flexiblities + non-Hertzian 

         if (ic%mater.eq.3 .and. potcon%ipotcn.ge.1) then

            ! Compute semi-axes of equivalent ellipse

            call equiv_ellipse(igs1, a_equiv, b_equiv)

            ! M=3: calculate flexibilities L1,L2,L3 from Hertzian contact ellipse

            if (idebug.ge.5) call write_log('panprc: calling simpflex')
            call simpflex(ic, a_equiv, b_equiv, mater, fric%fstat(), kin, 0)

            ! fill the constant part in the rigid slip

            if (idebug.ge.5) call write_log('panprc: calling set_tang_rhs')
            call set_tang_rhs(ic, mater, cgrid, kin, geom)
         endif

         ! In steady state problems, copy element division igs to igv

         if (is_ssrol) call eldiv_copy(igs1, outpt1%igv, ikALL)

         ! Clear the exterior of traction
         ! (Note: snorm doesn't involve itself with the tangential problem)

         ncon = 0
         do ii = 1, npot
            ! if (ii.eq.1) write(*,'(i5,i3,f12.6)') ii, igs1%el(ii), ps1%vn(ii)
            if (igs1%el(ii).le.Exter) then
               ps1%vn(ii) = 0d0
               ps1%vx(ii) = 0d0
               ps1%vy(ii) = 0d0
            else
               ncon = ncon + 1
            endif
         enddo

         if (ic%tang.eq.0) then

            ! Normal problem only: Pan.process is done / "converged"

            dif = 0d0

         elseif (ncon.le.0) then

            ! Empty contact area: call dummy Tang implementation, Pan.process is done / "converged"

            if (idebug.ge.5) call write_log('panprc: call stang_empty')
            call stang_empty(ic, mater, cgrid, fric, solv, outpt1, it, hs1, infl)
            dif = 0d0

         elseif (ic%mater.eq.2 .or. ic%mater.eq.3) then

            ! Using Fastsim: Compute the new tangential traction, Pan.process is done

            if (idebug.ge.5) call write_log('panprc: call stang_fastsim...')
            call stang_fastsim (ic, mater, cgrid, fric, kin, solv, outpt1, it, hs1, infl)
            if (idebug.ge.5) call write_log('panprc: returned from stang_fastsim')
            solv%ittang = solv%ittang + it
            dif = 0d0

         else

            ! Else: Compute the new tangential traction

            if (idebug.ge.5) call write_log('panprc: call stang...')
            call stang (ic, mater, cgrid, fric, kin, solv, outpt1, it, hs1, infl)
            if (idebug.ge.5) call write_log('panprc: returned from stang')
            if (it.ge.0) then
               solv%ittang = solv%ittang + it
            else
               ! maxin reached in TANG
               solv%ittang = -1
            endif

            ! Compute difference of tractions Ps1, Ps2 with tractions of previous iter.
            ! Po1 = Po1 - Ps1;  dif = || Po1 || ;  Po1 := Ps1

            if (idebug.ge.5) call write_log('panprc: computing po-ps')
            call gf3_axpy(AllElm, -1d0, ps1, po1, ikALL)

            dif   = gf3_rms(AllInt, po1, ikALL)
            difid = 5d0 * solv%eps * gf3_rms(AllInt, ps1, ikALL)

            call gf3_copy(AllElm, ps1, po1, ikALL)

            if (solv%maxout.gt.1) then

               ! In case of non-quasiidentity: write convergence

               if (ic%flow.ge.1 .or. itout.ge.solv%maxout-4) then
                  strng(1) = fmt_gs(12,4,4, dif)
                  strng(2) = fmt_gs(12,4,4, difid)
                  write (bufout, 6000) itout, (strng(j),j=1,2)
                  call write_log(1, bufout)
                  if (ic%flow.ge.2) call write_log(' ')
 6000             format (i3, ', Panag/Out: |Pk - Pk-1|,  5 Eps |Pk| :',2a12)
               endif
            endif
         endif

      end do ! while (not converged AND it<max)

      if (dif.gt.difid .and. solv%maxout.gt.1) then
         write (bufout, 7001) solv%maxout
         call write_log(1, bufout)
 7001    format (' Panag: Outer Loop Diverges. MaxOut = ',i3)

         if (.true.) then
            strng(1) = fmt_gs(12,4,4, dif)
            strng(2) = fmt_gs(12,4,4, difid)
            write (bufout, 6000) itout, (strng(j),j=1,2)
            call write_log(1, bufout)

            write (bufout, *) 'itnorm=',solv%itnorm,', ittang=', solv%ittang
            call write_log(1, bufout)
         endif
      endif

      if (idebug.ge.5) write(*,*) 'panprc: calling gf3_destroy'
      call gf3_destroy(po1)
      if (idebug.ge.5) write(*,*) 'panprc: done, returning'

      call timer_stop(itimer_panprc)
      end associate
   end subroutine panprc

!------------------------------------------------------------------------------------------------------------

   subroutine set_dummy_solution(gd)
!--purpose: initialize output structures: grid, eldiv, pressures, etc.
      implicit none
!--subroutine arguments:
      type(t_probdata) :: gd

      ! assuming that the potcon data are ok

      call init_curr_grid(gd%ic, gd%potcon_inp, gd%potcon_cur, gd%cgrid_cur, gd%geom, gd%outpt1)
      call set_prev_data(gd%ic, gd%geom, gd%mater, gd%fric, gd%kin, gd%outpt1)

      call eldiv_exter(gd%outpt1%igs)
      call gf3_set(AllElm, 0d0, gd%outpt1%ps, ikALL)

   end subroutine set_dummy_solution

!------------------------------------------------------------------------------------------------------------

end module m_scontc
