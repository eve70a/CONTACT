!------------------------------------------------------------------------------------------------------------
! m_stang - set-up of tangential contact problem and TANG algorithm calling the actual solvers
!
! Copyright 1979-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_stang

use m_hierarch_data
use m_aijpj
use m_solvpt
use m_leadedge
use m_temperature

implicit none
private

public  stang
private stang_rhs
private stang_nmdbg
public  stang_fastsim
public  stang_empty

contains

!------------------------------------------------------------------------------------------------------------
   subroutine stang (ic, mater, cgrid, fric, kin, solv, outpt1, ittang, hs1, infl)
!--purpose: Tang calculates the tangential traction ps(i,ik), i=1..npot, ik=1,2, with the normal
!           pressure and contact area fixed. This normal pressure is not affected.
!           The kind of problem (shift, transient/steady rolling) is specified by ic%tang.
      implicit none
!--subroutine parameters :
      type(t_ic)               :: ic
      type(t_material)         :: mater
      type(t_grid)             :: cgrid
      type(t_friclaw)          :: fric
      type(t_kincns)           :: kin
      type(t_solvers)          :: solv
      type(t_output)           :: outpt1
      type(t_influe)           :: infl
      type(t_gridfnc3)         :: hs1
      integer                  :: ittang
!--local variables :
      type(t_gridfnc3)            :: musted, slpvel, slpprv, tmpprv, tmp, wsfix1
      type(t_eldiv)               :: igsprv
      type(t_leadedge)            :: ledg
      integer, allocatable        :: iel(:)
      logical                     :: zready, zreasl, is_roll, is_ssrol, use_plast, use_out_it
      integer                     :: i, ii, ix, ixsta, ixinc, ixend, iy, j, k, icount, info, it,        &
                                     itslp, nadh, nslip, nplast, nexter, newins,    &
                                     newadh, smller, smllww, imeth, itgs, frclaw
      real(kind=8)                :: ww, pabs, dif, difid, errpt, tol, tol1, tol2, theta, omgslp
      character(len=12)           :: strng(2)
      character(len=4), parameter :: namits(1:3) = (/ 'ItGS', 'ItCG', 'ItGD' /)

      ! Initializations:

      call timer_start(itimer_stang)
      associate( mx   => cgrid%nx,   my   => cgrid%ny,    npot => cgrid%ntot, cs    => infl%cs,         &
                 cv   => infl%cv,    igs1  => outpt1%igs, mus1  => outpt1%mus, muv1  => outpt1%muv,     &
                 ps1  => outpt1%ps,  ss1   => outpt1%ss,  temp1 => outpt1%temp1)

      is_roll   = ic%tang.eq.2 .or. ic%tang.eq.3 
      is_ssrol  = ic%tang.eq.3
      use_plast = ic%mater.eq.4 .and. mater%tau_c0.gt.1d-10 .and. mater%tau_c0.le.1d10
      frclaw    = fric%frclaw_eff
      use_out_it = ((frclaw.ge.2 .and. frclaw.le.4) .or. frclaw.eq.6)

      ! dup: nullify, copy structure and el.div from ps1, initialize at 0:

      call gf3_copy_struc(ps1, wsfix1, 'tang:wsfix1', .true.)
      call gf3_copy_struc(ps1, tmp,    'tang:tmp',    .true.)
      call gf3_copy_struc(ps1, musted, 'tang:musted', .true.)
      call gf3_copy_struc(ps1, slpvel, 'tang:slpvel', .true.)
      call gf3_copy_struc(ps1, slpprv, 'tang:slpprv', .true.)
      call gf3_copy_struc(ps1, tmpprv, 'tang:tmpprv', .true.)

      call eldiv_new(igsprv, cgrid, nulify=.true.)
      allocate(iel(npot))

      ! set start/increment/end for elements per grid row, increasing in rolling direction

      if (is_roll .and. abs(kin%chi-pi).lt.0.01d0) then
         ixsta = mx
         ixinc = -1
         ixend =  1
      else
         ixsta =  1
         ixinc =  1
         ixend = mx
      endif

      ! write type of problem to trace output

      if (ic%flow.ge.4 .or. (ic%flow.ge.2 .and. use_out_it)) then
         if (ic%mater.ge.2 .and. ic%mater.le.3) then
            call write_log(' TANG: STEADY STATE ROLLING BY THE FASTSIM APPROACH')
         elseif (ic%tang.eq.1) then
            call write_log(' TANG: SHIFT TRANSIENT')
         elseif (ic%tang.eq.2) then
            call write_log(' TANG: TRANSIENT ROLLING CONTACT')
         elseif (ic%tang.eq.3) then
            call write_log(' TANG: STEADY STATE ROLLING CONTACT')
         endif
      endif

      ! Count the number of elements "k" in the contact area, and fill
      !    indirect adressing array iel i --> ii, with the i-th active
      !    element inside the contact area being element ii of the pot.con.

      k = 0
      do ii = 1, npot
         if (igs1%el(ii).ge.Adhes) then
            k = k + 1
            iel(k) = ii
         endif
      enddo
      if (k.le.0) then
         call write_log(' INTERNAL ERROR: TANG called with empty contact')
         call abort_run()
      endif

      ! Count the number of elements at the trailing edge
      ! at the edge ix=1 or ix=mx of the potential contact area.

      if (is_roll .and. abs(kin%chi-pi).lt.0.01d0) then
         ixsta = mx
      else
         ixsta = 1
      endif
      icount = 0
      do iy = 1, my
         ii = ixsta + (iy-1)*mx
         if (igs1%el(ii).ge.Adhes) icount = icount + 1
      enddo

      ! Set the actual solver to be used
      !   M  = 2, 3 : FASTSIM
      !   G  = 0, 4 : default solver
      !                T = 3            : SteadyGS
      !                T = 1, 2, M = 4  : ConvexGS
      !                T = 1, 2, M <> 4 : TangCG,
      !   G  = 3    : same, fixed omegah, omegas
      !   G  = 2    : ConvexGS, omegah, omegas
      !   G  = 5    : GDsteady, omegas

      if (mater%mater_eff.eq.2 .or. mater%mater_eff.eq.3) then
         solv%solver_eff = isolv_fastsm
      elseif (is_ssrol) then
         if (solv%gausei_eff.eq.5 .and. ixsta.eq.1) then
            if (ixsta.ne.1) then
               write(bufout,*) 'chi=',kin%chi,', ixsta=',ixsta,', mx=',mx
               call write_log(1, bufout)
            endif
            solv%solver_eff = isolv_gdstdy      ! G=5
         elseif (solv%gausei_eff.ne.2) then
            solv%solver_eff = isolv_stdygs      ! G=0, 3 or 4
         else
            solv%solver_eff = isolv_cnvxgs      ! G=2
         endif
      else
         if (solv%gausei_eff.ne.2 .and. mater%mater_eff.ne.4) then
            solv%solver_eff = isolv_tangcg      ! G=0, 3, 4 or 5
         else
            solv%solver_eff = isolv_cnvxgs      ! G=2 or M=4
         endif
      endif

      ! Force use of ConvexGS for T=3 when guard band is too small

      if (icount.gt.0 .and. (solv%solver_eff.eq.isolv_stdygs .or. solv%solver_eff.eq.isolv_gdstdy)) then
         if (solv%solver_eff.eq.isolv_stdygs) write(bufout,1101) icount, '(faster/more robust) SteadyGS'
         if (solv%solver_eff.eq.isolv_gdstdy) write(bufout,1101) icount, '(faster) GDsteady'
         if (ic%flow.ge.3) call write_log(2, bufout)
         solv%solver_eff = isolv_cnvxgs
      endif
 1101 format(' TANG: Warning: pot.con. too small; no exterior points at trailing edge at',i4,' rows.',/, &
             '       Using solver ConvexGS instead of ',a,'.')

      ! Apply default values for relaxation parameters when G-digit=0 or 4.
      ! G=5: SteadyGS used as fall-back when GDsteady stagnates.

      if (solv%gausei_eff.eq.0 .or. solv%gausei_eff.eq.4 .or. solv%gausei_eff.eq.5) then
         if (k.le.25) then
            ! Small system of eqs
            solv%omegah = 1.0d0
            solv%omegas = 1.0d0
         elseif (is_ssrol .and. solv%solver_eff.eq.isolv_cnvxgs) then
            ! Solver ConvexGS, steady state rolling
            solv%omegah = 0.5d0
            solv%omegas = 0.5d0
         elseif (is_ssrol) then
            ! Solver SteadyGS
            if (cgrid%dx/cgrid%dy.le.5d0) then
               solv%omegah = 0.9d0
               solv%omegas = 1.0d0
            elseif (cgrid%dx/cgrid%dy.le.15d0) then
               solv%omegah = 0.8d0
               solv%omegas = 0.8d0
            else
               solv%omegah = 0.8d0
               solv%omegas = 0.6d0
            endif
         else
            ! Solver TangCG or ConvexGS, shift/transient rolling
            solv%omegah = 0.5d0
            solv%omegas = 0.5d0
         endif
      endif

      ! Administrate the leading edge positions for all rows iy

      call sxbnd(ic, solv, cgrid, igs1, kin%chi, kin%dq, ledg)

      !------------------------------------------------------------------------------------------------------
      ! Start of the outer loop for velocity-dependent friction:
      !
      !      set initial estimate for slpvel
      !      while (fixed point not achieved, |slpvel-slpprv|>tol)
      !         - copy slpprv := slpvel
      !         - compute steady state friction coefficient musted(slpvel)
      !         - compute actual friction coefficient mus1 using time evolution
      !         - solve a rolling contact problem (transient or steady)
      !           with fixed traction bound with the TANG algorithm
      !         - determine new estimate for slip velocity, slpvel
      !         - compute update |slpvel - slpprv|
      !         - determine whether slpvel is converged

      ! Set prescribed coefficients of friction

      if (frclaw.eq.0) then

         call gf3_copy_xdir(AllElm, fric%nvf, fric%fstat_arr, mus1, ikXDIR)

      endif

      ! Set initial estimate for absolute slip velocity
      !  - if IniSlp == 0: use solution S of previous case,
      !  - if IniSlp > 0: approximate final S from below,
      !  - if IniSlp < 0: approximate final S from above.
      ! note that ss == shift distance [mm], slpvel == velocity [m/s]

      do ii = 1, npot
         if (solv%inislp.eq.0) then
            slpvel%vt(ii) = sqrt(outpt1%sv%vx(ii)**2 + outpt1%sv%vy(ii)**2) / kin%dt
         elseif (solv%inislp.ge.1) then
            slpvel%vt(ii) = 0d0
         else
            slpvel%vt(ii) = kin%veloc
         endif
      enddo

      ! Set initial estimate for surface temperature

      if (frclaw.eq.6) call calc_temp(ic, mater, kin, cgrid, outpt1)

      ! initialize iteration counts for Tang and CnvxGS

      ittang = 0
      itgs   = 0

      itslp  = 0
      zreasl = .false.
      omgslp = solv%omgslp
      ! write(bufout,'(a,f8.4)') 'Initial OmgSlp=',solv%omgslp
      ! call write_log(1, bufout)

      do while (itslp.le.max(1,solv%maxin-10) .and. .not.zreasl)
         itslp  = itslp + 1
         zreasl = .true.

         ! Copy slip velocity slpvel to slpprv, temperature to tmpprv,
         ! Copy element division igs1 to igsprv.

         call gf3_copy(AllElm, slpvel, slpprv, ikXDIR)
         if (frclaw.eq.6) call gf3_copy(AllElm, temp1,  tmpprv, ikZDIR)
         call eldiv_copy(igs1, igsprv, ikALL)

         ! L=2,3,4: Set "steady state" friction coefficient \mu_s on basis of current slip velocity
         ! L=6:     Set "steady state" friction coefficient \mu_s on basis of current temperature

         do iy = 1, my
            do ix = 1, mx
               ii = ix + (iy-1) * mx
               if (frclaw.ge.2 .and. frclaw.le.4) then
                  musted%vt(ii) = fric%mu_steady(slpvel%vt(ii), iy)
               elseif (frclaw.eq.6) then
                  musted%vt(ii) = fric%mu_temp(temp1%vn(ii), iy)
               endif
            enddo
         enddo

         ! L=2,3,4, 6: Set the actual friction coefficient \mu

         if ((frclaw.ge.2 .and. frclaw.le.4) .or. frclaw.eq.6) then

            ! Falling friction: incorporate friction memory
            ! Integrate mus1 from leading ("mx" / ixend) to trailing edge ("1" / ixsta)

            do iy = 1, my
               ii = ixend + (iy-1)*mx
               if (is_ssrol) muv1%vt(ii) = fric%fstat(iy)
               do ix = ixend, ixsta, -ixinc
                  ii = ix + (iy-1)*mx
                  if (fric%memdst.le.1d-9) then
                     theta = 0d0
                  else
                     theta = fric%memdst / (fric%memdst + kin%dt * max(slpvel%vt(ii), fric%mem_s0))
                     if (.false. .and. ii.ge.750 .and. ii.le.760)                                      &
                        write(*,122) theta, fric%memdst, fric%memdst, kin%dt, max(slpvel%vt(ii), fric%mem_s0)
 122                    format('theta=',f6.3,' =',f7.4,' / (',f7.4,' +', f7.4,' *',f6.3,' )')
                  endif
                  mus1%vt(ii) = theta*muv1%vt(ii) + (1d0-theta)*musted%vt(ii)
                  if (.false. .and. ii.ge.750 .and. ii.le.760)                                         &
                     write(*,123) ix, mus1%vt(ii), theta, muv1%vt(ii), (1d0-theta), musted%vt(ii)
 123              format('ix=',i3,': mu=',f7.4,' =',f6.3,' *',f7.4,' +',f6.3,' *',f7.4)
                  if (is_ssrol .and. ix.ne.ixsta)                                                      &
                     muv1%vt(ii-ixinc) = mus1%vt(ii)
               enddo
            enddo

         endif

         !----------------------------------------------------------------------------------------------------
         ! The element division in Igs is used as initial estimate.
         !
         ! Central in the TANG algorithm stands the shift S per element I and tangential direction t (=2,3).
         ! It consists of 7 terms:
         !                                          ,      ,                ,
         !    S    = W           +  un  + ut   -  un   - ut    + upl   - upl     (1a)
         !     It     It              It    It      It     It       It      It
         !
         !                                         bnd     bnd              ,
         !    S    = W  . facdt  +  un  + ut   - un    - ut    + upl   - upl     (1b)
         !     It     It              It    It     It      It       It      It
         !
         ! (1a) is used in the interior, (1b) at the leading edge boundary.
         ! The fixed terms are collected in wsfix, the other terms are computed later on in solvpt,
         ! depending on the problem that is solved (T-digit, F-digit, K-parameter).

         call stang_rhs(ic, cgrid, igs1, k, iel, ledg, hs1, ps1, outpt1%pv, cs      , cv      , wsfix1)

         ! The fixed part wsfix = W** of the shift S_It is now complete.

         !----------------------------------------------------------------------------------------------------
         ! Start of the main loop of the TANG algorithm:
         !
         !      While (not ready .and. more its. allowed) do
         !         zready = True
         !         solve Pt (info)
         !         if (traction bound violated)
         !            zready = False
         !            expand slip area S, adjust H
         !         elseif (slip directions violated)
         !            zready = False
         !            contract slip area S, adjust H
         !         if (info > 0)
         !            zready = False
         !      end while

         zready = .false.
         do while (.not.zready .and. ittang.lt.solv%maxin)

            ! start new iteration of TANG algorithm

            ittang = ittang + 1
            zready = .true.

            call eldiv_count(igs1, nadh, nslip, nplast, nexter)
            if (ic%flow.ge.4 .or. (ic%flow.ge.3 .and. .not.use_out_it .and. ittang.ge.2)) then
               if (use_plast) then
                  write(bufout,5102) ittang, nadh, nslip, nplast
               else
                  write(bufout,5103) ittang, nadh, nslip
               endif
               call write_log(1, bufout)
            endif
 5102       format (4x,i4, ', Tang: size of Adhes, Slip, Plast area : ', 3i7)
 5103       format (4x,i4, ', Tang: size of Adhes, Slip area : ', 2i7)
            if (ic%flow.ge.9) call wrigs (igs1, is_roll, kin%chi)

            ! T4: Step T4 of TANG algorithm: solve tractions.

            ! Compute tangential tractions for current Slip, Adhesion areas. (Note that CnvxGS may change
            !     the element division as well.)

            if (ic%x_nmdbg.ge.5) then
               call stang_nmdbg ('problem to solvpt:', 5, ic, cgrid, igs1, ledg, ps1, mus1, tmp,        &
                                 wsfix1, tmp, tmp, tmp, ss1)
            endif

            call solvpt(ic, mater, cgrid, npot, k, iel, fric, kin, solv, wsfix1, infl, ledg, outpt1,    &
                        imeth, info, it, errpt)
            itgs = itgs + it

            if (ic%x_nmdbg.ge.5) then
               call stang_nmdbg ('solution of solvpt:', 6,                                              &
                                 ic, cgrid, igs1, ledg, ps1, mus1, tmp, wsfix1, tmp, tmp, tmp, ss1)
            endif

            if (ic%x_nmdbg.ge.5) then
               call stang_nmdbg ('shift distance:', 7,                                                  &
                                 ic, cgrid, igs1, ledg, ps1, mus1, tmp, wsfix1, tmp, tmp, tmp, ss1)
            endif

            if (info.ge.3) then
               write(bufout,'(a)') ' The Gauss-Seidel-based solver diverges. ',                        &
                                   ' You may re-try with a larger value of DQ. ',                      &
                                   ' In case of serious problems, please report to us.'
               call write_log(3, bufout)
               zready = .true.
               ittang = -1
            endif

            ! T5,6:  Step T5, T6 of TANG algorithm: check constraints, restore.

            ! Expand the slip area, contract the adhesion area.
            ! Search for i's in Adhes with |traction(i)| > f*pn = tract.bound.

            newins = 0
            smller = 0
            tol    = dsqrt (2d0) * errpt
            do i = 1, k
               ii = iel(i)

               if (igs1%el(ii).eq.Adhes) then
                  pabs = dsqrt (ps1%vx(ii)**2 + ps1%vy(ii)**2)
                  if (pabs.ge.mus1%vt(ii)*ps1%vn(ii)+tol) then
                     newins  = newins + 1
                     igs1%el(ii) = Slip
                     ps1%vx(ii) = ps1%vx(ii) * mus1%vt(ii) * ps1%vn(ii) / pabs
                     ps1%vy(ii) = ps1%vy(ii) * mus1%vt(ii) * ps1%vn(ii) / pabs
                  endif

                  if (abs(pabs-mus1%vt(ii)*ps1%vn(ii)).lt.tol) then
                     smller = smller + 1
                  endif
               endif
            enddo

            ! if elements were moved from Adhesion to Slip:
            !  - another iteration of TANG is required

            if (newins.ne.0) zready = .false.

            ! We may cycle as well when CnvxGS didnot converge yet.
            ! if (smller.ne.0 .and. info.gt.0) zready = .false.

            if (zready) then

               ! T7,8:  Steps T7, T8 of TANG algorithm: check Lagrange multipliers, release constraints.

               ! Expand the adhesion area, contract the slip area.
               ! Search for i's in Slip with shift(i)*traction(i) > 0,
               ! tractions not opposed to direction of shift.

               ! Compute the displacements "tol1,tol2" in the center element ii==(mx/2,my/2) due to
               ! uniform tractions of size errpt

               call areas (igs1)
               ii = max(1,mx/2) + (max(1,my/2)-1)*mx
               call gf3_set(AllElm, 1d0, tmp, ikTANG)
               tol1 = errpt * 2d0 * abs (AijPj(ii, ikXDIR, tmp, jkXDIR, cs))
               tol2 = errpt * 2d0 * abs (AijPj(ii, ikYDIR, tmp, jkYDIR, cs))

               if (ic%x_nmdbg.ge.5) then
                  call stang_nmdbg('after expanding slip-area:', 8,                                     &
                                   ic, cgrid, igs1, ledg, ps1, mus1, tmp, wsfix1, tmp, tmp, tmp, ss1)
               endif

               ! Compare the directions of tractions and the shift

               newadh = 0
               smllww = 0
               do i = 1, k
                  ii = iel(i)

                  if (igs1%el(ii).eq.Slip) then
                     ww = ss1%vx(ii) * ps1%vx(ii) + ss1%vy(ii) * ps1%vy(ii)
                     tol = tol1*abs(ps1%vx(ii)) + tol2*abs(ps1%vy(ii))                                  &
                           + errpt * (abs(ss1%vx(ii)) + abs(ss1%vy(ii)))
                     if (tol.lt.0d0) then
                        write(bufout,*) 'TANG: INTERNAL ERROR: Tol=',tol,'<0'
                        call write_log(1, bufout)
                        call abort_run()
                     endif
                     if (ww.gt.tol) then
                        igs1%el(ii) = Adhes
                        newadh = newadh + 1
                     elseif (abs(ww).le.tol) then
                        smllww = smllww + 1
                     endif
                  endif
               enddo

               ! if elements were moved from Slip to Adhesion:
               !  - another iteration of TANG is required

               if (newadh.ne.0) zready = .false.

               ! if (CnvxGS solver did not converge): abort iteration with approximate solution

               if (info.gt.0) zready = .true.

            endif ! if (zready): steps T7,8

         ! end while (Main loop of TANG algorithm)

         enddo

         call eldiv_count(igs1, nadh, nslip, nplast, nexter)

         if (.not.zready) then
            write(bufout, 7800) solv%maxin
            call write_log(2, bufout)
 7800       format (/, ' TANG: ERROR. MaxIn = ', i5, ' reached.')
            ittang = -1
         endif

         ! end of the main loop of the TANG algorithm.
         ! At this point, the pressure Ps satisfies the constraints.

         !---------------------------------------------------------------------------------------------------
         ! Outer iteration for velocity-dependent friction:

         if (.not.use_out_it) then

            ! No iteration required, friction law is not dependent on slip velocity nor temperature

            zreasl = .true.

         elseif (frclaw.ge.2 .and. frclaw.le.4) then

            !  - Count the number of changes to the element division

            icount = 0
            do ii = 1, npot
               if (igs1%el(ii).ne.igsprv%el(ii)) icount = icount + 1
            enddo
            ! if (icount.gt.0) write(*,*) 'Changed element division at',icount,' points'

            !  - Compute new estimate for absolute slip velocity sabs = shft / dt

            do ii = 1, npot
               if (igs1%el(ii).eq.Exter) then
                  slpvel%vt(ii) = 0d0
               else
                  slpvel%vt(ii) = sqrt(ss1%vx(ii)**2 + ss1%vy(ii)**2) / kin%dt
               endif
               ! if (ii.eq.156) write(*,'(a,i3,3(a,g14.6))') 'ii=',ii, ': mu_s=',mus1%vt(ii),            &
               !    ', sr_x=',ss1%vx(ii), ', sabs=',slpvel%vt(ii),', sprv=', slpprv%vt(ii)
            enddo

            if (.false. .and. ittang.ge.190) then
               do ii=1,npot
                  write(lout,1234) ii, igs1%el(ii), mus1%vt(ii), ps1%vx(ii), ps1%vy(ii), slpvel%vt(ii)
 1234             format(i4,i3,4f12.6)
               enddo
            endif

            ! Compute the update |slpvel-slpprv|

            call gf3_axpy(AllElm, -1d0, slpvel, slpprv, ikXDIR)
            dif   = gf3_rms(AllInt, slpprv, ikXDIR)
            difid = 0.001d0 * gf3_rms(AllInt, slpvel, ikXDIR)

            !  - Apply relaxation: Slpvel = Slpvel0 + (omgslp-1) * (Slpvel0 - Slpprv)
            !    Note: slpprv contains the update "slpvel0 - slpprv"
            !    Note: relaxation only in slip area where mu is increasing.

            if (itslp.gt.1) then
               ! icnt = 0
               do ii = 1, npot
                  if (igs1%el(ii).eq.Slip .and. mus1%vt(ii).gt.muv1%vt(ii)) then
                     ! icnt = icnt + 1
                     slpvel%vt(ii) = max(0d0,slpvel%vt(ii)-(omgslp-1d0)*slpprv%vt(ii))
                  endif
               enddo
            endif
            ! write(*,*) 'relaxation in',icnt,' points'

            strng(1) = fmt_gs(12, 4, 4, dif)
            strng(2) = fmt_gs(12, 4, 4, difid)
            if (ic%flow.ge.4) then
               write(bufout,8601) itslp, (strng(j), j=1,2)
            elseif (ic%flow.ge.2 .and. .not.use_plast) then
               write(bufout,8602) itslp, nadh, nslip, (strng(j), j=1,2)
            elseif (ic%flow.ge.2) then
               write(bufout,8603) itslp, nadh, nslip, nplast, (strng(j), j=1,2)
            endif
            if (ic%flow.ge.2) call write_log(1, bufout)
 8601       format (i5, ', Vel-dep: |Sk-Sk-1|, .001 |Sk| : ', 2a12)
 8602       format (i5, ', Vel-dep: A, S=',   2i6,', |Sk-Sk-1|, .001 |Sk| : ', 2a12)
 8603       format (i5, ', Vel-dep: A, S, P=',3i6,', |Sk-Sk-1|, .001 |Sk| : ', 2a12)

            !  - Set convergence flag

            if (dif.gt.difid .and. nslip.gt.0) zreasl = .false.

            !  - In case of potential alternation of elements, reduce omgslp

            if (.not.zreasl .and. icount.gt.0 .and. mod(itslp,20).eq.0) then
               omgslp = 0.8 * omgslp
               ! write(bufout,'(a,f8.4)') 'Reducing OmgSlp to',omgslp
               ! call write_log(1, bufout)
            endif

         !---------------------------------------------------------------------------------------------------
         ! Outer iteration for temperature-dependent friction:

         elseif (frclaw.eq.6) then

            !  - Count the number of changes to the element division

            icount = 0
            do ii = 1, npot
               if (igs1%el(ii).ne.igsprv%el(ii)) icount = icount + 1
            enddo

            if (ic%x_nmdbg.ge.1 .and. itslp.ge.10 .and. icount.gt.0) then
               write(bufout,'(2(a,i4),a)') '   it',itslp,': element division changed at',icount,' points'
               call write_log(1, bufout)
            endif

            !  - Compute new temperature for new slip solution

            if (ic%heat.ge.1) call calc_temp (ic, mater, kin, cgrid, outpt1)

            ! Compute the update |temp1-tmpprv|

            call gf3_axpy(AllElm, -1d0, temp1, tmpprv, ikZDIR)
            dif   = gf3_rms(AllInt, tmpprv, ikZDIR)
            difid = 0.0001d0 * gf3_rms(AllInt, temp1, ikZDIR)

            !  - Apply relaxation: Temp1 = Temp1 + (omgslp-1) * (Temp1 - Tmpprv)
            !    Note: tmpprv contains the update "temp1 - tmpprv"
            !    Note: relaxation only in slip area where mu is increasing.

            if (itslp.gt.1) then
               ! icnt = 0
               do ii = 1, npot
                  if (igs1%el(ii).eq.Slip .and. mus1%vt(ii).gt.muv1%vt(ii)) then
                     ! icnt = icnt + 1
                     temp1%vn(ii) = max(0d0, temp1%vn(ii) - (omgslp-1d0)*tmpprv%vn(ii))
                  endif
               enddo
            endif
            !  write(*,*) 'relaxation in',icnt,' points'

            strng(1) = fmt_gs(12, 4, 4, dif)
            strng(2) = fmt_gs(12, 4, 4, difid)
            if (ic%flow.ge.4) then
               write(bufout,8701) itslp, (strng(j), j=1,2)
            elseif (ic%flow.ge.2 .and. .not.use_plast) then
               write(bufout,8702) itslp, nadh, nslip, (strng(j), j=1,2)
            elseif (ic%flow.ge.2) then
               write(bufout,8703) itslp, nadh, nslip, nplast, (strng(j), j=1,2)
            endif
            if (ic%flow.ge.2) call write_log(1, bufout)
 8701       format (i5, ', Temp-dep: |Tk-Tk-1|, .001 |Tk| : ', 2a12)
 8702       format (i5, ', Temp-dep: A, S=',   2i6,', |Tk-Tk-1|, .001 |Tk| : ', 2a12)
 8703       format (i5, ', Temp-dep: A, S, P=',3i6,', |Tk-Tk-1|, .001 |Tk| : ', 2a12)

            !  - Set convergence flag

            if (dif.gt.difid .and. nslip.gt.0) zreasl = .false.

            !  - In case of potential alternation of elements, reduce omgslp

            if (.not.zreasl .and. icount.gt.0 .and. mod(itslp,15).eq.0) then
               omgslp = 0.8 * omgslp
               write(bufout,'(a,f8.4)') 'Reducing OmgSlp to',omgslp
               call write_log(1, bufout)
            endif

         endif

      ! End while (more outer iterations for slip velocity required)

      enddo

      !------------------------------------------------------------------------------------------------------
      ! Report on final solution

      if (ic%flow.ge.3) then
         if (use_plast) then
            write(bufout,8902) nadh, nslip, nplast, namits(imeth), itgs
         else
            write(bufout,8903) nadh, nslip, namits(imeth), itgs
         endif
         call write_log(1, bufout)
      endif
 8902 format (6x, 'Tang: final element division: A, S, P=', 3i7, 2x,a,'=', i6)
 8903 format (6x, 'Tang: final element division: A, S=', 2i7, 2x,a,'=', i6)
      solv%itgs = itgs

      if (smllww+smller.ge.5 .and. ic%ilvout.ge.2) then
         write(bufout, 7900) smllww, smller
         call write_log(3, bufout)
 7900    format(' TANG: WARNING. There are', i6, ' elements with small slip',/,                        &
                                 22x, 'and', i6, ' elements with tractions near the traction bound',/)
      endif

      if (ic%flow.ge.6 .and. omgslp.lt.solv%omgslp) then
         write(bufout,'(6x,a,f8.4)') 'Tang: OMGSLP was reduced to', omgslp
         call write_log(1, bufout)
      endif
      call timer_stop(itimer_stang)

      ! compute sensitivity w.r.t. creepages

      if (ic%sens.ge.3) then
         call timer_start(itimer_sens)
         call sens_tang (ic, cgrid, mater, fric, kin, solv, infl, ledg, outpt1, npot, k, iel, wsfix1)
         call timer_stop(itimer_sens)
      endif

      if (ic%flow.ge.4 .or. (ic%flow.ge.3 .and. (use_out_it .or. ittang.ge.2))) call write_log(' ')

      call gf3_destroy(wsfix1)
      call gf3_destroy(tmp)
      call gf3_destroy(musted)
      call gf3_destroy(slpvel)
      call gf3_destroy(slpprv)
      call gf3_destroy(tmpprv)
      call eldiv_destroy(igsprv)
      call leadedge_destroy(ledg)
      deallocate(iel)
      end associate

   end subroutine stang

!------------------------------------------------------------------------------------------------------------
   subroutine stang_rhs(ic, cgrid, igs, k, iel, ledg, hs, ps, pv, cs, cv, wsfix)
!--purpose: compute the right hand side of the tangential problem, i.e., the fixed part wsfix of the shift S.
      implicit none
!--subroutine parameters :
      type(t_ic)               :: ic
      type(t_grid)             :: cgrid
      type(t_eldiv)            :: igs
      type(t_leadedge)         :: ledg
      type(t_inflcf)           :: cs, cv
      type(t_gridfnc3)         :: tmp, hs, ps, pv, wsfix
      integer                  :: k
      integer, dimension(:)    :: iel(k)
!--local variables :
      type(t_gridfnc3)          :: wsrig, usn, ust, uvn, uvt
      integer                   :: i, ii, iy, ik, j, iidum
      logical                   :: is_ssrol
      character(len=100)        :: namdbg
!--functions used :
      integer ix4ii, iy4ii
!--statement functions for computing ix,iy from ii:
      ix4ii(iidum) = mod(iidum-1,cgrid%nx)+1
      iy4ii(iidum) = (iidum-1)/cgrid%nx+1

      ! Initializations:

      associate(facdt => ledg%facdt)
      call gf3_new(wsrig, 'tang:wsrig', cgrid, igs, .true.)
      call gf3_new(uvn,   'tang:uvn',   cgrid, igs, .true.)
      call gf3_new(uvt,   'tang:uvt',   cgrid, igs, .true.)
      call gf3_new(usn,   'tang:usn',   cgrid, igs, .true.)
      call gf3_new(ust,   'tang:ust',   cgrid, igs, .true.)

      is_ssrol = ic%tang.eq.3

      ! Central in the TANG algorithm stands the shift S per element I and tangential direction t (=2,3).
      ! It consists of 5 terms:
      !                                          ,      ,
      !    S    = W           +  un  + ut   -  un   - ut            (1a)
      !     It     It              It    It      It     It
      ! 
      !                                         bnd     bnd
      !    S    = W  . facdt  +  un  + ut   - un    - ut            (1b)
      !     It     It              It    It     It      It
      ! 
      ! (1a) is used in the interior, (1b) at the leading edge boundary.
      ! The fixed terms are collected in wsfix, the other terms are computed later on in solvpt,
      ! depending on the problem that is solved (T-digit, F-digit, K-parameter).
      ! 
      ! here I      stands for an element (i, ii),
      !      t      stands for tau, i.e. ik,
      !      W      is the rigid slip, prescribed by Cksi, Ceta, Cphi and ExRhs,
      !              stored in array Hs, except for unknown Cksi, Ceta,
      !      u      is the tangential displacement difference of the new time instance
      !              (consisting of usn+ust),
      !      u'     is the tangential displacement difference of the previous time instance
      !              (consisting of uvn+uvt).
      !      facdt  is the fraction of the time-step that points are in the contact area
      !      u^bnd  is the displacement difference at the leading edge
      ! 
      ! Equation (1b) is applied near the leading edge, i.e. where ii2j(ii) > 0.
      ! Equation (1a) is applied in the interior. In equation (1a) a factor facdt is used as well with facdt=1.
      ! 
      ! The displacement differences are computed using the tractions and influence coefficients:
      ! 
      !      u    =   sum   sum    A        p                        (2)
      !       It         jk    J    It Jjk   J jk
      ! 
      !       ,                     ,        ,
      !      u    =   sum   sum    A        p                        (3)
      !       It         jk    J    It Jjk   J jk
      ! 
      ! Here J ranges over all elements of the contact area at current and previous times
      ! and jk ranges over all directions 1,2,3.
      ! (Note that u is actually split into un + ut, with un corresponding to jk=3 and ut to jk=1,2.)
      ! The matrix A corresponds to influence coefficients cs.
      ! The matrix A' corresponds to influence coefficients cv.
      ! In stationary problems also csv may be used, which yields cs - cv.
      ! 
      ! Several terms in (1) - (3) do not depend on the unknown tractions p_j for jk=1,2,
      ! and are fixed during the solution process.
      ! These terms are gathered in wsfix. This is the constant W** from eq (4.55c) of [Kalker90].
      ! 
      ! wsfix is only computed for elements ii == iel(i) in the contact area.
      ! 
      ! 0) update the element divisions in row1st, rowlst for use in aijpj

      call areas(igs)

      ! 1) Initialize wsrig == -rigid shift  =  -hs . facdt
      !    Note: additional terms due to unknown creepages are added in solvpt.
      !    Note: facdt is 1 in situations/elements where no leading edge correction is required.

      do ik = 1, 2
         do i = 1, k
            ii = iel(i)
            wsrig%val(ii,ik) = - facdt%vt(ii) * hs%val(ii,ik)
         enddo
      enddo

      if (ic%x_nmdbg.ge.8) then
         namdbg = 'rigid shift):'
         call stang_nmdbg (namdbg, 1, ic, cgrid, igs, ledg, ps, tmp, wsrig, wsfix, usn, uvn, uvt, tmp)
      endif

      ! 2) Compute tangential displacements usn and uvn at times t and t' due to normal pressures ps and pv.

      ! Note: using igs to describe where uvn is required.
      ! Note: in steady state problems, ps has not yet been copied to pv.

      call VecAijPj (igs,AllInt, usn, ikTANG, ps, jkZDIR, cs)

      if (.not.is_ssrol) then
         call VecAijPj (igs,AllInt, uvn, ikTANG, pv, jkZDIR, cv)
      else
         call VecAijPj (igs,AllInt, uvn, ikTANG, ps, jkZDIR, cv)
      endif

      if (ic%x_nmdbg.ge.9) then
         namdbg = 'shift due to normal problem:'
         call stang_nmdbg (namdbg, 2, ic, cgrid, igs, ledg, ps, tmp, wsrig, wsfix, usn, uvn, uvt, tmp)
      endif

      ! 3) In transient problems, compute tangential displacements uvt of previous time t',
      !    due to tangential traction pv.

      if (.not.is_ssrol) then

         ! T=1, 2: contribution of tang.tractions in uvt.
         ! Note: using igs to describe where uvt is required.

         call VecAijPj (igs,AllInt, uvt, ikTANG, pv, jkTANG, cv)

         if (ic%x_nmdbg.ge.9) then
            namdbg = 'shift due to previous tang.traction:'
            call stang_nmdbg (namdbg, 3, ic, cgrid, igs, ledg, ps, tmp, wsrig, wsfix, usn, uvn, uvt, tmp)
         endif
      endif

      ! 4) Estimate ledg%ubnd, the displacement difference at the leading edge of the contact area
      !    Note: complete UBnd ("ubn+ubt"), due to pv in x,y,z-directions.
      !    Note: using cs, as if pv is defined on the current grid.

      call subnd(cgrid, pv, cs, ledg)

      ! 5) Add the relevant contributions of u, u' and ubnd to wsfix.

      do ik = 1, 2
         do i = 1, k
            ii = iel(i)

            ! Farther away from the leading edge:
            ! (note: leading edge correction is not used when T=1)

            if (ledg%ii2j%val(ii).le.0) then

               ! using equation (1a) as described above, w + u - u'
               ! when T=3, uvt is determined later on in solvpt

               if (.not.is_ssrol) then
                  wsfix%val(ii,ik) = wsrig%val(ii,ik) + usn%val(ii,ik) - uvn%val(ii,ik) - uvt%val(ii,ik)
               else
                  wsfix%val(ii,ik) = wsrig%val(ii,ik) + usn%val(ii,ik) - uvn%val(ii,ik)
               endif

            ! Near the leading edge:

            else

               ! using equation (1b) as described above, w + u - ubnd
               ! when T=3, ubnd is determined later on in solvpt

               if (.not.is_ssrol) then
                  iy = iy4ii(ii)
                  j  = ledg%ii2j%val(ii)
                  wsfix%val(ii,ik) = wsrig%val(ii,ik) + usn%val(ii,ik) - ledg%ubnd(j,ik)
                  if (.false. .and. iy.eq.21) then
                     write(*,1234) 'iy=',iy, ', ix=',ix4ii(ii), ', ik=',ik,                            &
                        ': applying Ubnd=',ledg%ubnd(j,ik), ' instead of Uv=',uvn%val(ii,ik)+uvt%val(ii,ik)
 1234                format(1x,3(a,i3),2(a,f9.6))
                  endif
               else
                  wsfix%val(ii,ik) = wsrig%val(ii,ik) + usn%val(ii,ik) - uvn%val(ii,ik)
               endif
            endif
         enddo
      enddo

      if (ic%x_nmdbg.ge.5) then
         namdbg = 'contributions to rhs:'
         call stang_nmdbg (namdbg, 4, ic, cgrid, igs, ledg, ps, tmp, wsrig, wsfix, usn, uvn, uvt, tmp)
      endif

      ! The fixed part wsfix = W** of the shift S_It is now complete.

      call gf3_destroy(wsrig)
      call gf3_destroy(usn)
      call gf3_destroy(ust)
      call gf3_destroy(uvn)
      call gf3_destroy(uvt)
      end associate

   end subroutine stang_rhs

!------------------------------------------------------------------------------------------------------------

   subroutine stang_nmdbg (messg, itable, ic, cgrid, igs, ledg, ps, mus, wsrig, wsfix, usn, uvn, uvt, ss)
!--purpose: Produce "nmdbg-output" for different phases in the TANG algorithm
      implicit none
!--subroutine parameters :
      type(t_ic)             :: ic
      type(t_grid)           :: cgrid
      type(t_eldiv)          :: igs
      type(t_gridfnc3)       :: ps, mus, wsrig, wsfix, usn, uvn, uvt, ss
      type(t_leadedge)       :: ledg
      integer                :: itable
      logical                :: is_ssrol
      character(len=*)       :: messg
!--local variables :
      logical, parameter    :: print_tables = .true.
      integer               :: ii, ik, j, iidum
      real(kind=8)          :: pabs, sabs, uabs, uu, ww, angl, angls, anglu
      character(len=1), parameter :: cdir(1:3) = (/ 'x', 'y', 'n' /)
!--functions used :
      integer ix4ii, iy4ii
!--statement functions for computing ix,iy from ii:
      ix4ii(iidum) = mod(iidum-1,cgrid%nx)+1
      iy4ii(iidum) = (iidum-1)/cgrid%nx+1

      is_ssrol = ic%tang.eq.3

      ! write message, e.g. describing the phase of the TANG algorithm

      write(lout,*) 'Nmdbg: ',trim(messg)

      ! if "table 1" is requested: prescribed rigid slip, lead.edge

      if (itable.eq.1 .and. print_tables) then
         write(lout,101)
 101     format('  ix,  iy,   ii,    wsrig_x,    wsrig_y,  facdt')
         do ii = 1, cgrid%ntot
            if (igs%el(ii).ge.Adhes) then
               write(lout,102) ix4ii(ii), iy4ii(ii), ii, wsrig%vx(ii), wsrig%vy(ii), ledg%facdt%vt(ii)
 102           format(2(i4,','),i5,',',2(f11.7,','),f7.3)
            endif
         enddo
      endif

      ! if "table 2" is requested: shift due to normal problem

      if (itable.eq.2 .and. print_tables) then
         write(lout,201)
 201     format('  ix,  iy,   ii,    unx        -unx'',       uny        -uny''')
         do ii = 1, cgrid%ntot
            if (igs%el(ii).ge.Adhes) then
               write(lout,202) ix4ii(ii), iy4ii(ii), ii, usn%vx(ii), uvn%vx(ii), usn%vy(ii), uvn%vy(ii)
 202           format(2(i4,','),i5,',',4(f11.7,:,','))
            endif
         enddo
      endif

      ! if "table 3" is requested: shift due to previous tang.traction

      if (itable.eq.3 .and. print_tables) then
         write(lout,301)
 301     format('  ix,  iy,   ii,   -utx'',      -uty''')
         do ii = 1, cgrid%ntot
            if (igs%el(ii).ge.Adhes) then
               write(lout,302) ix4ii(ii), iy4ii(ii), ii, uvt%vx(ii), uvt%vy(ii)
 302           format(2(i4,','),i5,',',2(f11.7,:,','))
            endif
         enddo
      endif

      ! if "table 4" is requested: break-down of right hand side

      if (itable.eq.4 .and. print_tables) then
         do ik = 1, 2
            if (is_ssrol)      write(lout,401) cdir(ik), cdir(ik), cdir(ik), cdir(ik)
            if (.not.is_ssrol) write(lout,401) cdir(ik), cdir(ik), cdir(ik), cdir(ik), cdir(ik)
 401        format('  ix,  iy,   ii,   wsfix_',a,' =   wsrig_',a,                                      &
                   '  +  usn_',a,'    -  uvn_',a,:,'    -  uvt_',a)
 402        format(2(i4,','),i5,',',5(f14.10,:,','))
 403        format(2(i4,','),i5,',',3(f14.10,:,','),4x,'0.0,',f14.10)

            do ii = 1, cgrid%ntot
               if (igs%el(ii).ge.Adhes) then
                  if (is_ssrol) then

                     ! T=3, Steady state rolling:
                     !      uvt is added later on, just as ubnd

                     write(lout,402) ix4ii(ii), iy4ii(ii), ii, wsfix%val(ii,ik), wsrig%val(ii,ik),     &
                        usn%val(ii,ik), -uvn%val(ii,ik)

                  elseif (ledg%ii2j%val(ii).le.0) then

                     ! T=1,2, Farther away from the leading edge:

                     write(lout,402) ix4ii(ii), iy4ii(ii), ii, wsfix%val(ii,ik), wsrig%val(ii,ik),     &
                        usn%val(ii,ik), -uvn%val(ii,ik), -uvt%val(ii,ik)
                  else

                     write(lout,402) ix4ii(ii), iy4ii(ii), ii, wsfix%val(ii,ik), wsrig%val(ii,ik),     &
                        usn%val(ii,ik), -uvn%val(ii,ik), -uvt%val(ii,ik)

                     ! T=1,2, Near the leading edge:

                     j  = ledg%ii2j%val(ii)
                     write(lout,403) ix4ii(ii), iy4ii(ii), ii, wsfix%val(ii,ik), wsrig%val(ii,ik),     &
                        usn%val(ii,ik), -ledg%ubnd(j,ik)
                  endif
               endif
            enddo
         enddo
      endif

      ! if "table 5" is requested: problem to solvpt:

      if (itable.eq.5 .and. print_tables) then
         write(lout,501)
 501     format('  ix,  iy,   ii,   wsfix_x,    wsfix_y,   init C, init.tract px, py')
         do ii = 1, cgrid%ntot
            if (igs%el(ii).ge.Adhes) then
               write(lout,502) ix4ii(ii), iy4ii(ii), ii, wsfix%vx(ii), wsfix%vy(ii),                   &
                  igs%el(ii), ps%vx(ii), ps%vy(ii)
 502           format(2(i4,','),i5,',',2(f11.7,','),i6,',',2(f11.4,:','))
            endif
         enddo
      endif

      ! if "table 6" is requested: solution of solvpt:

      if (itable.eq.6 .and. print_tables) then
         write(lout,601)
 601     format('  ix,  iy,   ii,   H/S,  tractions px,   py,       |pt|,       g-|pt|')
         do ii = 1, cgrid%ntot
            if (igs%el(ii).ge.Adhes) then
               pabs = dsqrt(ps%vx(ii)**2 + ps%vy(ii)**2)
               write(lout,602) ix4ii(ii), iy4ii(ii), ii, igs%el(ii), ps%vx(ii), ps%vy(ii), pabs,       &
                  mus%vt(ii)*ps%vn(ii)-pabs
 602           format(2(i4,','),i5,',',i6,',',4(f11.4,:','))
            endif
         enddo
      endif

      ! if "table 7" is requested: shift distance:

      if (itable.eq.7 .and. print_tables) then
         do ik = 1, 2
            if (is_ssrol) then
               write(lout,701) cdir(ik), ' (u-u'')', cdir(ik), cdir(ik)
            else
               write(lout,701) cdir(ik), '      u', cdir(ik), cdir(ik)
            endif
 701        format('  ix,  iy,   ii,   wsfix_',a,'  +',' ',a,'_',a,' = shift ss_',a)

            do ii = 1, cgrid%ntot
               if (igs%el(ii).ge.Adhes) then
                  write(lout,702) ix4ii(ii), iy4ii(ii), ii, wsfix%val(ii,ik),                          &
                     ss%val(ii,ik)-wsfix%val(ii,ik), ss%val(ii,ik)
 702              format(2(i4,','),i5,',',3(f11.7,:,','))
               endif
            enddo
         enddo
      endif

      ! if "table 8" is requested: after expanding slip-area:

      if (itable.eq.8 .and. print_tables) then
         write(lout,801)
 801     format('  ix,  iy,   ii, H/S,  tractions px,   py,       slip sx,     sy,   angl(s->p)')
         do ii = 1, cgrid%ntot
            if (igs%el(ii).ge.Adhes) then
               ww = (ss%vx(ii) * ps%vx(ii) + ss%vy(ii) * ps%vy(ii))
               pabs = max(1d-9, dsqrt(ps%vx(ii)**2+ps%vy(ii)**2))
               sabs = max(1d-9, dsqrt(ss%vx(ii)**2+ss%vy(ii)**2))
               angl = acos(ww/(pabs*sabs)) * 180d0/pi
               if (isnan(angl)) then
                  write(lout,*) 'NaN:', ps%vx(ii), ps%vy(ii), pabs
                  write(lout,*) '    ', ss%vx(ii), ss%vy(ii), sabs
                  write(lout,*) '    ', ww, ww/(pabs*sabs), angl
               endif
               write(lout,802) ix4ii(ii), iy4ii(ii), ii, igs%el(ii), ps%vx(ii), ps%vy(ii),             &
                  ss%vx(ii), ss%vy(ii), angl
 802           format(2(i4,','),i5,',',i4,',',2(f11.4,','),2(f11.7,','),f8.1)
            endif
         enddo
      endif

      ! if "table 9" is requested: plastic deformation:

      if (itable.eq.9 .and. print_tables) then
         write(lout,901)
 901     format('  ix,  iy,   ii, H/S/P, plast_def uplx, uply,   |upl|, angl(upl->p), angl(s->p), yield')
         do ii = 1, cgrid%ntot
            if (igs%el(ii).ge.Adhes) then
               ww = (ss%vx(ii) * ps%vx(ii) + ss%vy(ii) * ps%vy(ii))
               uu = (usn%vx(ii) * ps%vx(ii) + usn%vy(ii) * ps%vy(ii))
               uabs = max(1d-9, dsqrt(usn%vx(ii)**2+usn%vy(ii)**2))
               pabs = max(1d-9, dsqrt(ps%vx(ii)**2+ps%vy(ii)**2))
               sabs = max(1d-9, dsqrt(ss%vx(ii)**2+ss%vy(ii)**2))
               angls = acos(ww/(pabs*sabs+1.d-15)) * 180d0/pi ! TODO: test if the 1.d-15 is really nescesarry
               anglu = acos(uu/(pabs*uabs+1.d-15)) * 180d0/pi ! TODO: test if the 1.d-15 is really nescesarry
               if (isnan(angl)) then
                  write(lout,*) 'NaN:', ps%vx(ii), ps%vy(ii), pabs
                  write(lout,*) '    ', ss%vx(ii), ss%vy(ii), sabs
                  write(lout,*) '    ', ww, ww/(pabs*sabs), angl
               endif
               write(lout,902) ix4ii(ii), iy4ii(ii), ii, igs%el(ii), usn%vx(ii), usn%vy(ii),            &
                  uabs, anglu, angls, uvn%vt(ii)
 902           format(2(i4,','),i5,',',i4,',',2(f11.4,','),f9.5,',',f8.1,',   ',f8.1,',' f11.7)
            endif
         enddo
      endif

   end subroutine stang_nmdbg

!------------------------------------------------------------------------------------------------------------

   subroutine stang_fastsim (ic, mater, cgrid, fric, kin, solv, outpt1, ittang, hs1, infl)
!--purpose: Tang calculates the tangential traction ps(i,ik), i=1..npot, ik=1,2, with the normal
!           pressure and contact area fixed. This normal pressure is not affected.
!           The kind of problem (shift, transient/steady rolling) is specified by ic%tang.
      implicit none
!--subroutine parameters :
      type(t_ic)               :: ic
      type(t_material)         :: mater
      type(t_grid)             :: cgrid
      type(t_friclaw)          :: fric
      type(t_kincns)           :: kin
      type(t_solvers)          :: solv
      type(t_output)           :: outpt1
      type(t_influe)           :: infl
      type(t_gridfnc3)         :: hs1
      integer                  :: ittang
!--local variables :
      type(t_gridfnc3)            :: tmp, wsfix1
      type(t_leadedge)            :: ledg
      integer, allocatable        :: iel(:)
      logical                     :: is_roll, is_ssrol
      integer                     :: ii, k, info, it, nadh, nslip, nplast, nexter, imeth, itgs
      real(kind=8)                :: errpt
      character(len=4), parameter :: namits(1:2) = (/ 'ItGS', 'ItCG' /)

      ! Initializations:

      call timer_start(itimer_stang)

      associate(npot => cgrid%ntot, igs1 => outpt1%igs, mus1 => outpt1%mus, ps1  => outpt1%ps,  &
                ss1  => outpt1%ss)

      is_roll  = ic%tang.eq.2 .or. ic%tang.eq.3
      is_ssrol = ic%tang.eq.3

      ! dup: nullify, copy structure and el.div from ps1, initialize at 0:

      call gf3_copy_struc(ps1, wsfix1, 'tang:wsfix1', .true.)
      call gf3_copy_struc(ps1, tmp,    'tang:tmp',    .true.)

      allocate(iel(npot))

      ! check restrictions

      if (ic%mater.lt.2 .or. ic%mater.gt.3) then
         call write_log(' ERROR: tang_fastsim should only be called for Fastsim.')
         call abort_run()
      endif
      if (ic%tang.lt.2 .or. ic%tang.gt.3) then
         call write_log(' ERROR: tang_fastsim should only be called for rolling, T=2 or 3.')
         call abort_run()
      endif

      ! write type of problem to trace output

      if (ic%flow.ge.3) then
         if (ic%tang.eq.2) then
            call write_log(' TANG: TRANSIENT ROLLING CONTACT USING FASTSIM')
         elseif (ic%tang.eq.3) then
            call write_log(' TANG: STEADY STATE ROLLING BY THE FASTSIM APPROACH')
         endif
      endif

      ! Count the number of elements "k" in the contact area, and fill indirect adressing array
      !   iel i --> ii, with the i-th active element inside the contact area being element ii of the pot.con.

      k = 0
      do ii = 1, npot
         if (igs1%el(ii).ge.Adhes) then
            k = k + 1
            iel(k) = ii
         endif
      enddo
      if (k.le.0) then
         call write_log(' INTERNAL ERROR: TANG called with empty contact')
         call abort_run()
      endif

      ! Set fixed part wsfix = W** of the right hand side - unknown creepages are added later on in solvpt

      call gf3_copy(AllElm,  hs1, wsfix1, ikTANG)
      call gf3_scal(AllElm, -1d0, wsfix1, ikTANG)

      ! Set facdt to appropriate value: 1.0

      call gf3_new(ledg%facdt, 'stang:ledg%facdt', cgrid)
      call gf3_set(AllElm, 1d0, ledg%facdt, ikXDIR)

      ! Set fstat as the initial value of the coefficients of friction mu, ignoring friction variation

      call gf3_set(AllElm,  fric%fstat(), mus1, ikXDIR)

      if (ic%x_nmdbg.ge.8) then
         call stang_nmdbg ('rigid shift', 1, ic, cgrid, igs1, ledg, ps1, tmp, wsfix1, wsfix1,            &
                           tmp, tmp, tmp, tmp)
      endif

      ! Set the actual solver to be used
      !   M  = 2, 3 : FASTSIM

      solv%solver_eff = isolv_fastsm
      ittang = 1
      itgs   = 0

      ! Compute tangential tractions & corresponding Adhesion/Slip areas using modified Fastsim-algorithm.

      if (ic%x_nmdbg.ge.5) then
         call stang_nmdbg ('problem to solvpt:', 5, ic, cgrid, igs1, ledg, ps1, mus1, tmp, wsfix1,      &
                           tmp, tmp, tmp, ss1)
      endif

      call solvpt(ic, mater, cgrid, npot, k, iel, fric, kin, solv, wsfix1, infl, ledg, outpt1, imeth,   &
                  info, it, errpt)

      if (ic%x_nmdbg.ge.5) then
         call stang_nmdbg ('solution of solvpt:', 6, ic, cgrid, igs1, ledg, ps1, mus1, tmp, wsfix1,     &
                           tmp, tmp, tmp, ss1)
         call stang_nmdbg ('shift distance:', 7, ic, cgrid, igs1, ledg, ps1, mus1, tmp, wsfix1,         &
                           tmp, tmp, tmp, ss1)
      endif

      call areas (igs1)
      call eldiv_count(igs1, nadh, nslip, nplast, nexter)

      ! Report on final solution

      if (ic%flow.ge.3) then
         write(bufout, 8900) nadh, nslip, namits(imeth), itgs
         call write_log(1, bufout)
         call write_log(' ')
      endif
 8900 format (6x, 'Tang: final element division: A, S=', 2i7, 2x,a,'=', i6)
      solv%itgs = itgs

      ! compute sensitivity w.r.t. creepages

      if (ic%sens.ge.3) then
         call timer_start(itimer_sens)
         call sens_tang (ic, cgrid, mater, fric, kin, solv, infl, ledg, outpt1, npot, k, iel, wsfix1)
         call timer_stop(itimer_sens)
      endif

      call gf3_destroy(wsfix1)
      call gf3_destroy(tmp)
      call leadedge_destroy(ledg)
      deallocate(iel)

      end associate
      call timer_stop(itimer_stang)

   end subroutine stang_fastsim

!------------------------------------------------------------------------------------------------------------

   subroutine stang_empty (ic, mater, cgrid, fric, solv, outpt1, ittang, hs1, infl)
!--purpose: Defaults for Tang for empty contact area.
      implicit none
!--subroutine parameters :
      type(t_ic)               :: ic
      type(t_material)         :: mater
      type(t_grid)             :: cgrid
      type(t_friclaw)          :: fric
      type(t_solvers)          :: solv
      type(t_output)           :: outpt1
      type(t_influe)           :: infl
      type(t_gridfnc3)         :: hs1
      integer                  :: ittang
!--local variables :
      integer                   :: nadh, nslip, itgs
      real(kind=8), save        :: rdum

      associate(igs1 => outpt1%igs)

      ! avoid compiler warning wrt unused variables

      if (allocated(infl%cs%cf) .and. associated(hs1%vn) .and. associated(igs1%el)) then
         rdum = mater%ga + hs1%vn(1) + real(igs1%el(1)) + real(cgrid%nx)
      endif

      ! write type of problem to trace output

      if (ic%flow.ge.3) then
         if (ic%mater.ge.2 .and. ic%mater.le.3) then
            call write_log(' TANG: STEADY STATE ROLLING BY THE FASTSIM APPROACH')
         elseif (ic%tang.eq.1) then
            call write_log(' TANG: SHIFT TRANSIENT')
         elseif (ic%tang.eq.2) then
            call write_log(' TANG: TRANSIENT ROLLING CONTACT')
         elseif (ic%tang.eq.3) then
            call write_log(' TANG: STEADY STATE ROLLING CONTACT')
         endif
      endif

      ! initialize iteration counts for Tang and CnvxGS

      ittang = 0
      itgs   = 0

      ! Set the actual friction coefficient \mu, ignoring friction variation

      call gf3_set(AllElm, fric%fstat(), outpt1%mus, ikXDIR)

      ! Report on final solution

      nadh = 0
      nslip = 0
      if (ic%flow.ge.3) write(bufout, 8900) nadh, nslip, itgs
      if (ic%flow.ge.3) call write_log(2, bufout)
 8900 format (6x, 'Tang: final element division: A, S=', 2i7,'  ItGS=', i6, /)
      solv%itgs = itgs

      end associate
   end subroutine stang_empty

!------------------------------------------------------------------------------------------------------------

end module m_stang
