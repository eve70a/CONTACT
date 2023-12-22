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

      is_fastsm = gd%ic%mater.eq.2 .or. gd%ic%mater.eq.3
      is_roll   = gd%ic%tang.eq.2 .or. gd%ic%tang.eq.3
      is_ssrol  = gd%ic%tang.eq.3

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

      call combin_mater(gd%mater)

      ! copy fntrue to fnscal when N=1, Fn prescribed

      if (gd%ic%norm.eq.1) gd%kin%fnscal = gd%kin%fntrue / gd%mater%ga

      ! (re)initialize penv when P=2

      if (gd%ic%pvtime.eq.2) gd%kin%penv = 0d0

      ! solve Hertzian problem when needed

      if (gd%potcon%ipotcn.ge.-5 .and. gd%potcon%ipotcn.le.-1) then
         if (idebug.ge.5) call write_log('contac: calling hzsol')
         call hzsol(gd%ic, gd%geom, gd%potcon, gd%hertz, gd%kin, gd%mater)
      endif

      ! M=3: calculate flexibilities L1,L2,L3 from Hertzian contact ellipse

      if (gd%ic%mater.eq.3 .and. gd%potcon%ipotcn.le.-1) then
         if (idebug.ge.5) call write_log('contac: calling simpflex')
         call simpflex(gd%ic, gd%hertz%aa, gd%hertz%bb, gd%mater, gd%fric%fstat(), gd%kin, 0)
      endif

      ! make new discretisation

      if (gd%ic%discns3.ge.1) then
         if (idebug.ge.5) call write_log('contac: calling sdis')
         call sdis(gd%potcon, gd%hertz, gd%cgrid)
      endif

      ! force CHI = 0 when using visco-elastic materials (T=3, M=1)

      if (is_ssrol .and. gd%ic%mater.eq.1 .and. abs(gd%kin%chi).gt.0.001d0) then
         write (bufout, 4002) gd%kin%chi, 'visco-elastic materials'
         call write_log(2, bufout)
         gd%kin%chi = 0d0
      endif

      if (is_ssrol .and. is_fastsm) then

         ! set DQ == DX when using steady rolling w. Fastsim (T=3, M=2/3)
         ! force CHI = 0 or pi when using steady rolling w. Fastsim (T=3, M=2/3)

         gd%kin%dq  = gd%cgrid%dx
         if (abs(gd%kin%chi).gt.0.01d0 .and. abs(gd%kin%chi-pi).gt.0.01d0) then
            write (bufout, 4002) gd%kin%chi, 'Fastsim'
            call write_log(2, bufout)
            gd%kin%chi = 0d0
         endif

      elseif ((gd%solv%gausei_eff.ne.2 .and. is_ssrol) .or. gd%ic%mater.eq.4) then

         ! force CHI = 0 or pi when using steady rolling w. SteadyGS (T=3, G<>2) or when M=4 (T=1,3)
         ! force DQ == DX when using steady rolling w. SteadyGS (T=3, G<>2), also when M=1 (visc.)

         if (abs(gd%kin%chi).gt.0.01d0 .and. abs(gd%kin%chi-pi).gt.0.01d0) then
            if (gd%ic%ilvout.ge.1) then
               if (gd%solv%gausei_eff.ne.2 .and. is_ssrol) then
                  write (bufout, 4002) gd%kin%chi, 'steady state rolling with solver SteadyGS'
               else
                  write (bufout, 4002) gd%kin%chi, 'solver ConvexGS with plasticity'
               endif
               call write_log(2, bufout)
            endif
            gd%kin%chi = 0d0
         endif
         if (abs(gd%kin%dq-gd%cgrid%dx).gt.0.01*gd%cgrid%dx) then
            if (gd%ic%ilvout.ge.1) then
               if (gd%solv%gausei_eff.ne.2 .and. is_ssrol) then
                  write (bufout, 4003) gd%cgrid%dx,gd%kin%dq,'steady state rolling with solver SteadyGS'
               else
                  write (bufout, 4003) gd%cgrid%dx,gd%kin%dq,'solver ConvexGS with plasticity'
               endif
               call write_log(2, bufout)
            endif
         endif
         gd%kin%dq  = gd%cgrid%dx
      endif

 4002 format (' Input: Warning: using CHI = 0 instead of',f6.1,' rad'/, 17x,                            &
              'for steady state rolling with ',a,'.')
 4003 format (' Input: Warning: Using DQ = DX =',g12.4, ' (instead of',g12.4,')', /, 17x,           &
              'for ',a,'.')

      ! set appropriate chi, dq, veloc for shifts

      if (.not.is_roll) then
         gd%kin%chi   = 0d0
         gd%kin%dq    = 1d0
         gd%kin%veloc = 1d0
      endif

      ! compute the time step size, with dq in [mm], veloc in [mm/s]

      if (is_roll .or. gd%ic%mater.ne.5) then
         gd%kin%dt = gd%kin%dq / max(1d-10, gd%kin%veloc)
      endif

      ! calculate influence functions cs, csv
      !  - compute standard infl.cf when C>=1
      !  - compute Blanco-correction when C=4
      !  - load numerical influence coefficients when C=9
      !     --> overwrite cs,csv when format=0, fill cy when format=1

      if (gd%ic%gencr_inp.ge.1) then
         if (idebug.ge.5) call write_log('contac: calling sgencr')
         call sgencr (gd%ic, gd%cgrid, gd%mater, gd%influ, gd%fric, gd%kin)
      endif

      if (gd%ic%gencr_inp.eq.4 .or. (gd%ic%gencr_inp.eq.1 .and. gd%mater%gencr_eff.eq.4)) then
         if (idebug.ge.5) call write_log('contac: calling influe_blanco')
         call influe_blanco (is_roll, gd%mater%if_meth, gd%mater%if_ver, gd%mater%ninclin,              &
                             gd%mater%surf_inclin, gd%cgrid, gd%influ)
      endif

      if (gd%ic%gencr_inp.eq.9 .or. (gd%ic%gencr_inp.eq.1 .and. gd%mater%gencr_eff.eq.9)) then
         if (idebug.ge.5) call write_log('contac: calling influe_load')

         call influe_load(gd%mater%fname_influe, gd%meta%dirnam, is_roll, gd%cgrid, gd%influ)
         call influe_mater(gd%influ, gd%mater)

         if (gd%ic%x_inflcf.eq.6) then
            write(lout,*) 'Printing influence coefficients cs (instat)'
            call inflcf_print (gd%influ%cs, lout)
         endif
      endif

      ! form undeformed distance hs1

      if (gd%ic%rznorm.gt.0) then
         if (idebug.ge.5) call write_log('contac: calling srznrm')
         call srznrm (gd%ic, gd%cgrid, gd%geom)
      endif

      ! set initial traction Ps and previous traction Pv :

      ihertz = 0
      if (gd%potcon%ipotcn.lt.0) ihertz = abs(gd%potcon%ipotcn)

      if (idebug.ge.5) call write_log('contac: calling filpvs(1)')
      call filpvs (gd%ic, gd%hertz, gd%potcon, gd%cgrid, gd%geom, gd%mater, gd%fric, gd%kin,            &
                gd%outpt1, ihertz)

      ! compute the constant part in the rigid slip, i.e. excluding creepages that are unknown
      ! note for Fastsim + 3 flexibilities + non-Hertzian: flexibilities not yet filled in

      if (gd%ic%mater.ne.3 .or. gd%potcon%ipotcn.lt.0) then
         if (idebug.ge.5) call write_log('contac: calling srztng')
         call srztng (gd%ic, gd%mater, gd%cgrid, gd%kin, gd%geom)
      endif

      ! skip computing part?

      if (gd%ic%return.ge.2) goto 90

      ! the computing part :

      ! Start the actual solution: Panagiotopoulos process

      if (idebug.ge.5) call write_log('contac: calling panprc')
      call panprc (gd%ic, ihertz, gd%mater, gd%potcon, gd%cgrid, gd%influ, gd%fric, gd%kin, gd%solv,    &
                   gd%outpt1, gd%geom)

      ! TODO: fill outputs uv, sv, uplv, taucv, etc. in steady rolling problems

      ! Compute surface temperatures

      if (gd%ic%heat.ge.1 .and. gd%fric%frclaw_eff.ne.6) then
         if (idebug.ge.5) call write_log('contac: calling calc_temp')
         call calc_temp (gd%ic, gd%mater, gd%kin, gd%cgrid, gd%outpt1)
      endif

      ! Derive quantities like total forces and elastic energy, write output to out-file when O>=1.

      if (idebug.ge.5) call write_log('contac: calling soutpt')
      call soutpt (gd)

      ! Output results to .mat-file if requested

      if (gd%ic%matfil_surf.ge.1) then
         if (idebug.ge.5) call write_log(' ...scontc: calling writmt')
         call timer_start(itimer_files)
         call writmt (gd%meta, gd%ic, gd%cgrid, gd%potcon%xl, gd%potcon%yl, gd%geom%hs1, gd%mater,      &
                      gd%fric, gd%kin, gd%outpt1, .false.)
         call timer_stop(itimer_files)
      endif

      ! Compute the subsurface elastic field for points specified in Xb, write output to out-file and subs-file

      if (gd%ic%stress.ge.1) then
         if (idebug.ge.5) call write_log('contac: calling subsur')
         call subsur(gd%meta, gd%ic, gd%mater, gd%cgrid, gd%outpt1%igs, gd%outpt1%ps, .false., gd%subs)
      endif

      if (idebug.ge.5) call write_log('contac: done, returning')
 90   continue

   end subroutine contac

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
      integer :: npot_old
      logical :: zerror

      nerror = 0
      ierror = 0

      if (gd%meta%ncase.le.1 .and. gd%ic%pvtime.le.1) then
         write(bufout, 2001) gd%ic%iestim
         call write_log(2, bufout)
 2001    format (' WARNING. In the first case contact must be initiated.'/,                            &
                 '          Digit P=',i3,' is overruled, set to 2.')
         gd%ic%pvtime = 2
      endif

      if (gd%meta%ncase.le.1 .and. gd%ic%iestim.ne.0 .and. gd%ic%iestim.ne.5) then
         write (bufout, 2004) gd%ic%iestim
         call write_log(2, bufout)
 2004    format (' WARNING. In the first case contact must be initiated.'/,                            &
                 '          Digit I=', i3,' is overruled, set to 0.')
         gd%ic%iestim = 0
      endif

      if (gd%ic%pvtime.ne.2) then
         npot_old = size(gd%outpt1%ps%val,1)
         if (npot_old.ne.gd%cgrid%ntot) then
            nerror = nerror + 1
            if (ierror.eq.0) ierror = 3001
            write(bufout, 3001) gd%ic%pvtime, npot_old, gd%cgrid%ntot
            call write_log(2, bufout)
 3001       format (' ERROR. In the input for a sequence, the number of elements MX, MY',/,            &
                    '        may not change. P=',i3,', old:',i6,', new:',i6,'.')
         endif
      endif

      if (gd%ic%iestim.ne.0) then
         npot_old = size(gd%outpt1%ps%val,1)
         if (npot_old.ne.gd%cgrid%ntot) then
            write(bufout, 3006) npot_old, gd%cgrid%ntot, gd%ic%iestim
            call write_log(2, bufout)
 3006       format (' WARNING. When the number of elements MX, MY changes (old:',i6,',',/,             &
                    '          new:',i6,'), the initial estimate cannot be used.  I=',i3,'.')
            gd%ic%iestim = 0
         endif
      endif

      if (gd%fric%nvf.ne.1 .and. gd%fric%nvf.ne.gd%cgrid%ny) then
         nerror = nerror + 1
         if (ierror.eq.0) ierror = 4001
         write(bufout, 4001) gd%fric%nvf, gd%cgrid%ny
         call write_log(1, bufout)
 4001    format (' Internal ERROR. The number of friction values NVF =',i6,' must be 1 or MY = ', i6,'.')
      endif

      if (gd%geom%iplan.eq.4) then
         zerror = .not.check_sorted( 'YSEP', gd%geom%npatch-1, gd%geom%ysep, .true. )
         if (zerror) then
            nerror = nerror + 1
            if (ierror.eq.0) ierror = 5001
            call write_log(nline_errmsg, errmsg)
         endif
      endif

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

            if (idebug.ge.5) call write_log('panprc: calling srztng')
            call srztng(ic, mater, cgrid, kin, geom)
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
                  strng(1) = fmt_gs(12,4, dif)
                  strng(2) = fmt_gs(12,4, difid)
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
            strng(1) = fmt_gs(12,4, dif)
            strng(2) = fmt_gs(12,4, difid)
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

      ! assuming that the grid data are ok

      call sdis(gd%potcon, gd%hertz, gd%cgrid)

      call eldiv_new(gd%outpt1%igs, gd%cgrid)
      call eldiv_exter(gd%outpt1%igs)

      call gf3_new(gd%outpt1%ps, 'outpt%ps', gd%cgrid, gd%outpt1%igs)
      call gf3_set(AllElm, 0d0, gd%outpt1%ps, ikALL)

   end subroutine set_dummy_solution

!------------------------------------------------------------------------------------------------------------

end module m_scontc
