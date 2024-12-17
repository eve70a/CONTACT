!============================================================================================================
! m_solvpt - solution methods for tangential contact problem
!
! Copyright 1993-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!============================================================================================================

module m_solvpt

use m_hierarch_data
use m_aijpj
use m_leadedge
use m_hertz

! provide for compilation with "internal parallelization" (within each contact patch) enabled or disabled
#if defined OMP_INSIDE
#define _INT_OMP $omp
#else
#define _INT_OMP
#endif

implicit none
private

   public  solvpt
   private tang_solver
   private char_creep
   public  sens_tang
   private fastsim_params
   private modif_fastsim
   private modif_fastrip
   private gdsteady
   private compute_dp
   private project_searchdir
   private apply_trcbnd
   private solve_elmtrc
   private compute_diagscaling
   private perform_linesearch
   private tangcg
   private cnvxgs
   private stdygs
   private elmtrc
   private plstrc

contains

!============================================================================================================

      subroutine solvpt(ic, mater, cgrid, npot, k, iel, fric, kin, solv, wsfix1, infl, ledg, outpt1,    &
                        imeth, info1, itgs, err)
!--purpose: Compute the tangential tractions with fixed Slip- and Adhesion-Area.
!           On exit: info1 = 0  o.k.
!                    1  maxgs reached
!                    2  maxgs reached and mean rate of convergence > 0.99
!                    3  maxgs reached and mean rate of convergence > 1
      implicit none
!--subroutine parameters:
      type(t_ic)                  :: ic
      type(t_grid)                :: cgrid
      type(t_material)            :: mater
      type(t_friclaw)             :: fric
      type(t_kincns)              :: kin
      type(t_solvers)             :: solv
      type(t_gridfnc3)            :: wsfix1
      type(t_influe)              :: infl
      type(t_leadedge)            :: ledg
      type(t_output)              :: outpt1
      integer,      intent(in)    :: npot, k, iel(npot)
      integer,      intent(out)   :: info1, imeth, itgs
      real(kind=8), intent(out)   :: err
!--local variables:
      integer      :: i, ii, j, itnr, nadh, nslip, nplast, nexter, ifxy, it, it1hlf, it2hlf
      logical      :: is_roll, is_ssrol, full_sliding
      real(kind=8) :: df, dfx, dfy, dcksi, dceta, fxk, fyk, fxkp1, fykp1, dfxk, dfyk, det, conv,        &
                      dpxavg, dpyavg, sxavg, syavg
      character(len=4), parameter  :: namits(1:3) = (/ 'ItGS', 'ItCG', 'ItGD' /)
      character(len=12) :: strng(4)

      type(t_gridfnc3)          :: wstot, tmp

      ! initialize, determine number of elements in slip

      associate(dxdy  => cgrid%dxdy,  igs1  => outpt1%igs, ps1   => outpt1%ps, ss1   => outpt1%ss,      &
                sens  => outpt1%sens, facdt => ledg%facdt)

      call gf3_new(wstot, 'solvpt:wstot', cgrid, nulify=.true.)
      call gf3_copy_struc(ps1, tmp, 'solvpt:tmp', .true.)

      is_roll  = ic%tang.eq.2 .or. ic%tang.eq.3
      is_ssrol = ic%tang.eq.3

      itgs  = 0
      info1 = 0
      call areas (igs1)
      call eldiv_count(igs1, nadh, nslip, nplast, nexter)

      !-----------------------------------------------------------------------
      ! method for adhesion and slip, linear and non-linear equations:
      !          iterative solver CnvxGS (matrix-free)
      !          Newton-Raphson iteration for total forces Fx, Fy specified
      !-----------------------------------------------------------------------

      ! add contributions of cksi, ceta to vector ws if not done so before
      ! Using work-array wstot instead of input-array wsfix.

      call gf3_copy(AllElm, wsfix1, wstot, ikTANG)
      if (ic%x_nmdbg.ge.1) call gf3_check_nan(wsfix1, 'solvpt: wsfix1', AllInt, ikTANG, 2)

      if (ic%force3.ge.1) then
         do i = 1, k
            ii = iel(i)
            wstot%vx(ii) = wstot%vx(ii) + facdt%vt(ii) * kin%cksi*kin%dq
         enddo
      endif

      if (ic%force3.ge.2) then
         do i = 1, k
            ii = iel(i)
            wstot%vy(ii) = wstot%vy(ii) + facdt%vt(ii) * kin%ceta*kin%dq
         enddo
      endif

      ! solve equations with initial cksi and ceta

      call tang_solver(npot, ic, cgrid, mater, fric, kin, solv, infl, ledg, outpt1, k, iel,             &
                         wstot, .false., info1, imeth, it, err, conv)

      fxkp1 = dxdy * gf3_sum(AllElm, ps1,1) / (kin%muscal * kin%fntrue)
      fykp1 = dxdy * gf3_sum(AllElm, ps1,2) / (kin%muscal * kin%fntrue)

      itgs = itgs + it

      call eldiv_count(igs1, nadh, nslip, nplast, nexter)

      ! if total forces Fx and/or Fy are prescribed: perform Newton-Raphson iteration

      if (ic%force3.ge.1 .and. info1.le.1) then

         ! Print output on initial solution

         itnr = 0
         strng(1) = fmt_gs(12, 4, 4, kin%cksi)
         strng(2) = fmt_gs(12, 4, 4, fxkp1)
         strng(3) = fmt_gs(12, 4, 4, kin%ceta)
         strng(4) = fmt_gs(12, 4, 4, fykp1)
         if (ic%force3.eq.1) then
            if (ic%flow.ge.4) write(bufout,7000) itnr, (strng(j),j=1,2), namits(imeth), it
            if (ic%flow.ge.4) call write_log(1, bufout)
         elseif (ic%force3.eq.2) then
            if (ic%flow.ge.4) write(bufout,8000) itnr, (strng(j),j=1,4), namits(imeth), it
            if (ic%flow.ge.4) call write_log(2, bufout)
         endif
         if (ic%flow.ge.3 .and. nadh.eq.0) write(bufout,9000) nadh, nslip
         if (ic%flow.ge.3 .and. nadh.eq.0) call write_log(1, bufout)
 7000    format(4x,i6,', NR_Fx,    Cksi, Fxk: ',2g12.4, '   ',a,':', i5)
 8000    format(4x,i6,', NR_Fy,    Cksi, Fxk: ',2g12.4, /,                              &
                             22x, 'Ceta, Fyk: ',2g12.4, '   ',a,':', i5)
 9000    format(9x,'Warning: full sliding, A, S=',2i7)

         ! While (|dF| > Eps|F|) do

         df = abs(kin%fxrel - fxkp1)
         if (ic%force3.ge.2) df = df + abs(kin%fyrel - fykp1)

         do while (df.gt.solv%eps .and. itnr.lt.solv%maxnr .and. info1.le.1)
            itnr = itnr + 1

            ! print sensitivities used, when requested

            if (.false. .and. (ic%flow.ge.5 .or. ic%x_nmdbg.ge.5)) then
               if (ic%force3.eq.1) then
                  strng(1) = fmt_gs(12, 4, 4, sens(iout_fx1,iin_dksi1))
                  write(bufout,'(18x,2a)') 'Sens:',strng(1)
                  call write_log(1, bufout)
               else
                  write(bufout,'(18x,a,4g12.3)') 'Sens:',                               &
                        sens(iout_fx1,iin_dksi1), sens(iout_fy1,iin_dksi1),             &
                        sens(iout_fx1,iin_deta1), sens(iout_fy1,iin_deta1)
                  call write_log(1, bufout)
               endif
            endif

            ! determine suitable changes dcksi, dceta for cksi and ceta

            dfx = kin%fxrel - fxkp1
            dfy = kin%fyrel - fykp1

            if (ic%force3.eq.1) then
               dceta = 0d0
               if (abs(sens(iout_fx1,iin_dksi1)).gt.1e-6) then
                  dcksi = dfx / sens(iout_fx1,iin_dksi1)
               else
                  dcksi = 0.00003
               endif
            elseif (ic%force3.eq.2) then
               det =   sens(iout_fx1,iin_dksi1) * sens(iout_fy1,iin_deta1)              &
                     - sens(iout_fy1,iin_dksi1) * sens(iout_fx1,iin_deta1)
               if (det.gt.solv%eps) then
                  dcksi = ( sens(iout_fy1,iin_deta1)*dfx - sens(iout_fy1,iin_dksi1)*dfy) / det
                  dceta = (-sens(iout_fx1,iin_deta1)*dfx + sens(iout_fx1,iin_dksi1)*dfy) / det
               else
                  dcksi = 0.000003d0
                  dceta = 0.000003d0
               endif
            endif

            ! In case of full sliding, estimate the creepage needed to bring average slip back to zero

            if (nadh.le.0 .and. itnr.le.1 .and. .false.) then
               sxavg = gf3_sum(AllInt, ss1, 1) / real(k)
               syavg = gf3_sum(AllInt, ss1, 2) / real(k)
               dcksi = dcksi - sxavg / kin%dq
               dceta = dceta - syavg / kin%dq
               full_sliding = .true.
            else
               full_sliding = .false.
            endif

            ! loop over dcksi and optionally dceta

            it1hlf = 0
            it2hlf = 0

            do ifxy = 1, ic%force3

               ! store "previous" total forces fxk, fyk

               fxk = fxkp1
               fyk = fykp1

               ! set new creepage cksi or ceta, add contribution to ws

               if (ifxy.eq.1) then
                  kin%cksi = kin%cksi + dcksi
                  do i = 1, k
                     ii = iel(i)
                     wstot%vx(ii) = wstot%vx(ii) + facdt%vt(ii) * dcksi * kin%dq
                  enddo
               else
                  kin%ceta = kin%ceta + dceta
                  do i = 1, k
                     ii = iel(i)
                     wstot%vy(ii) = wstot%vy(ii) + facdt%vt(ii) * dceta * kin%dq
                  enddo
               endif

               ! solve equations for new cksi or new ceta

               call tang_solver(npot, ic, cgrid, mater, fric, kin, solv, infl, ledg, outpt1, k, iel,    &
                         wstot, .false., info1, imeth, it, err, conv)

               fxkp1 = dxdy * gf3_sum(AllElm, ps1, 1) / (kin%muscal*kin%fntrue)
               fykp1 = dxdy * gf3_sum(AllElm, ps1, 2) / (kin%muscal*kin%fntrue)

               itgs = itgs + it
               if (ifxy.eq.1) then
                  it1hlf = it
               else
                  it2hlf = it
               endif

               call eldiv_count(igs1, nadh, nslip, nplast, nexter)

               ! compute (update) sensitivities w.r.t. cksi or ceta,
               !    except when full sliding occurs.

               dfxk = fxkp1 - fxk
               dfyk = fykp1 - fyk

               if (ifxy.eq.1) then

                  ! ifxy=1: update sensitivities w.r.t. cksi

                  if (nadh.gt.0) then
                     dpxavg = dfxk * kin%muscal*kin%fntrue / ((nadh+nslip+nplast)*dxdy)
                     if (sens(iout_fx1,iin_dksi1).eq.0d0 .or. abs(dpxavg).gt.10d0*err) then
                        sens(iout_fx1,iin_dksi1) = dfxk / dcksi
                     elseif (ic%flow.ge.9) then
                        call write_log ('         skipping update of sensitivity dF/dCksi')
                     endif
                     dpyavg = dfyk * kin%muscal*kin%fntrue / ((nadh+nslip+nplast)*dxdy)
                     if (abs(dpyavg).gt.10d0*err) then
                        sens(iout_fy1,iin_dksi1) = dfyk / dcksi
                     endif
                  else
                     ! full sliding:
                     if (abs(sens(iout_fx1,iin_dksi1)).lt.1e-6) sens(iout_fx1,iin_dksi1) = fxkp1 / kin%cksi
                  endif
               endif

               if (ifxy.eq.2) then

                  ! ifxy=2: update sensitivities w.r.t. ceta

                  if (nadh.gt.0) then
                     dpxavg = dfxk * kin%muscal*kin%fntrue / ((nadh+nslip+nplast)*dxdy)
                     if (abs(dpxavg).gt.10d0*err) then
                        sens(iout_fx1,iin_deta1) = dfxk / dceta
                     elseif (ic%flow.ge.9) then
                        call write_log ('         skipping update of sensitivity dF/dCeta')
                     endif
                     dpyavg = dfyk * kin%muscal*kin%fntrue / ((nadh+nslip+nplast)*dxdy)
                     if (sens(iout_fy1,iin_deta1).eq.0d0 .or. abs(dpyavg).gt.10d0*err) then
                        sens(iout_fy1,iin_deta1) = dfyk / dceta
                     endif
                  else
                     ! full sliding:
                     if (abs(sens(iout_fy1,iin_deta1)).lt.1e-6) sens(iout_fy1,iin_deta1) = fykp1 / kin%ceta
                  endif
               endif

               ! compute error w.r.t. requested total forces Fx, Fy

               df = abs(kin%fxrel - fxkp1)
               if (ic%force3.ge.2) df = df + abs(kin%fyrel - fykp1)

               ! print progress report for half iteration

               if (ic%force3.eq.2 .and. .false.) then
                  strng(1) = fmt_gs(12, 4, 4, kin%cksi)
                  strng(2) = fmt_gs(12, 4, 4, fxkp1)
                  strng(3) = fmt_gs(12, 4, 4, kin%ceta)
                  strng(4) = fmt_gs(12, 4, 4, fykp1)
                  write(bufout,8000) itnr, (strng(j),j=1,4), namits(imeth), it
                  call write_log(2, bufout)
               endif

            ! enddo (cksi and optionally ceta)

            enddo

            ! print progress report for full iteration

            strng(1) = fmt_gs(12, 4, 4, kin%cksi)
            strng(2) = fmt_gs(12, 4, 4, fxkp1)
            strng(3) = fmt_gs(12, 4, 4, kin%ceta)
            strng(4) = fmt_gs(12, 4, 4, fykp1)
            if (ic%force3.eq.1 .and. ic%flow.ge.4) then
               write(bufout,7000) itnr, (strng(j), j=1,2), namits(imeth), it1hlf
               call write_log(1, bufout)
            elseif (ic%force3.eq.2 .and. ic%flow.ge.4) then
               write(bufout,8000) itnr, (strng(j), j=1,4), namits(imeth), it1hlf+it2hlf
               call write_log(2, bufout)
            endif
            if (ic%flow.ge.3 .and. nadh.eq.0) then
               write(bufout,9000) nadh, nslip
               call write_log(1, bufout)
            endif

         ! end while (Newton-Raphson loop)

         enddo

         ! estimate Rms-error in ps in [N/mm2]: inaccuracy from GS and from NR.
         ! The factor 2 assumes that errors due to (Fx,Fy) are spread regularly but not equally over the grid points.

         err = err + 2d0*df * kin%muscal*kin%fntrue / ((nadh+slip+nplast)*dxdy)
      else

      ! else: cksi, ceta prescribed: reset sensitivities

         sens(iout_fx1,iin_dksi1) = 0d0
         sens(iout_fy1,iin_dksi1) = 0d0
         sens(iout_fx1,iin_deta1) = 0d0
         sens(iout_fy1,iin_deta1) = 0d0

      ! endif (Fx, Fy specified)

      endif

      call gf3_destroy(wstot)
      call gf3_destroy(tmp)
      end associate
      end subroutine solvpt

!------------------------------------------------------------------------------------------------------------

      subroutine tang_solver(npot, ic, cgrid, mater, fric, kin, solv, infl, ledg, outpt1, k, iel,       &
                         wstot, is_sens, info1, imeth, it, err, conv)
!--purpose: select the appropriate solver for the case
      implicit none
!--subroutine parameters:
      type(t_ic)               :: ic
      type(t_grid)             :: cgrid
      type(t_material)         :: mater
      type(t_friclaw)          :: fric
      type(t_kincns)           :: kin
      type(t_solvers)          :: solv
      type(t_influe)           :: infl
      type(t_leadedge)         :: ledg
      type(t_output)           :: outpt1
      type(t_gridfnc3)         :: wstot
      integer                  :: npot, k, iel(npot), info1, imeth, it
      logical                  :: is_sens
      real(kind=8)             :: err, conv
!--local variables:
      integer      :: maxit_loc, ic_flow
      logical      :: is_ssrol, lstagn
      real(kind=8) :: eps_loc

      is_ssrol  = ic%tang.eq.3

      if (is_sens) then
         maxit_loc = solv%mxsens
         eps_loc   = solv%epsens
      else
         maxit_loc = solv%maxgs
         eps_loc   = solv%eps
      endif

      if (solv%solver_eff.eq.isolv_fastsm) then
         imeth = 1
         call timer_start(itimer_fastsim)
         call modif_fastsim(ic, cgrid, mater, fric, kin, wstot, outpt1)
         it  = 0
         err = 0d0
         call timer_stop(itimer_fastsim)
      elseif (solv%solver_eff.eq.isolv_fastrp) then
         imeth = 1
         call timer_start(itimer_fastrip)
         call modif_fastrip(ic, cgrid, mater, fric, kin, wstot, outpt1)
         it  = 0
         err = 0d0
         call timer_stop(itimer_fastrip)
      elseif (solv%solver_eff.eq.isolv_cnvxgs) then
         imeth = 1
         call timer_start(itimer_cnvxgs)
         call cnvxgs(ic, cgrid, mater, kin, wstot, infl, ledg, outpt1, k, iel, eps_loc,                 &
                     maxit_loc, solv%omegah, solv%omegas, is_sens, info1, it, err, conv)
         call timer_stop(itimer_cnvxgs)
      elseif (solv%solver_eff.eq.isolv_tangcg) then
         call timer_start(itimer_tangcg)
         imeth = 2
         call tangcg(ic, npot, maxit_loc, eps_loc, is_sens, wstot, infl, outpt1, it, err)
         call timer_stop(itimer_tangcg)
      elseif (solv%solver_eff.eq.isolv_stdygs) then
         imeth = 1
         call timer_start(itimer_stdygs)
         call stdygs(ic, cgrid, mater, kin, wstot, infl, outpt1, k, eps_loc, maxit_loc,                 &
                     solv%omegah, solv%omegas, is_sens, info1, it, err, conv)
         call timer_stop(itimer_stdygs)
      elseif (solv%solver_eff.eq.isolv_gdstdy) then
         if (fric%frclaw_eff.eq.0 .and. .false.) then
            ! compute initial estimate using Fastsim (note: Fastsim changes mus when L=2--4)
            imeth = 1
            call timer_start(itimer_fastsim)
            ic_flow = ic%flow
            ic%flow = 0
            call modif_fastsim(ic, cgrid, mater, fric, kin, wstot, outpt1)
            ic%flow = ic_flow
            call timer_stop(itimer_fastsim)
         endif

         ! attempt solution using GDsteady
         imeth = 3
         call timer_start(itimer_gdstdy)
         call gdsteady(ic, npot, mater, maxit_loc, eps_loc, wstot, infl, outpt1, solv, it, err, lstagn)
         call timer_stop(itimer_gdstdy)

         if (lstagn .and. .true.) then
            ! switch to SteadyGS when detecting stagnation
            imeth = 1
            if (ic%flow.ge.2) then
               call write_log('        Switching to SteadyGS, restarting solution.')
            endif
            call timer_start(itimer_stdygs)
            call stdygs(ic, cgrid, mater, kin, wstot, infl, outpt1, k, eps_loc, maxit_loc,              &
                        solv%omegah, solv%omegas, is_sens, info1, it, err, conv)
            call timer_stop(itimer_stdygs)
         endif
      else
         write(bufout,*) 'INTERNAL ERROR: incorrect value for tangential solver:', solv%solver_eff
         call write_log(1, bufout)
         call abort_run()
      endif

      if (ic%x_nmdbg.ge.1) call gf3_check_nan(outpt1%ps, 'solvpt: ps', AllInt, ikTANG, 2)

      end subroutine tang_solver

!------------------------------------------------------------------------------------------------------------

      subroutine char_creep(igs, cgrid, mater, fstat, fntrue, ch_ksi, ch_eta, ch_phi, idebug)
!--purpose: compute characteristic size of the creepages using the linear theory of rolling
      implicit none
!--subroutine arguments:
      type(t_eldiv),    intent(in)      :: igs
      type(t_grid),     intent(in)      :: cgrid
      type(t_material), intent(in)      :: mater
      real(kind=8),     intent(in)      :: fstat, fntrue
      real(kind=8),     intent(out)     :: ch_ksi, ch_eta, ch_phi
      integer,          intent(in)      :: idebug
!--local variables:
      integer      :: ncon, nadh, nslip, nexter, nplast
      real(kind=8) :: a, b, aob, c, c11, c22, c23

      ! estimate effective radius c using contact area

      call eldiv_count(igs, nadh, nslip, nplast, nexter)
      ncon = nadh + nslip + nplast
      c = sqrt(real(ncon) * cgrid%dx * cgrid%dy / pi)

      ! estimate half sizes a,b using maximum/minimum elements in contact

      call areas(igs)
      a = real(max(1, igs%ixmax - igs%ixmin + 1)) * cgrid%dx * 0.5d0
      b = real(max(1, igs%iymax - igs%iymin + 1)) * cgrid%dy * 0.5d0
      aob = a / b

      ! compute the Kalker coefficients at aob, poiss

      call linrol(mater%nu, aob, c11, c22, c23)

      ! define characteristic sizes ch_ksi, ch_eta, ch_phi

      ch_ksi = fstat * fntrue / (c**2 * mater%ga * c11)
      ch_eta = fstat * fntrue / (c**2 * mater%ga * c22)
      ch_phi = fstat * fntrue / (c**3 * mater%ga * c23)

      if (idebug.ge.2) then
         write(bufout,120) ' characteristic cksi=', ch_ksi, ', ceta=', ch_eta,', cphi=',ch_phi
 120     format(3(a,f10.6))
         call write_log(1, bufout)
      endif

      end subroutine char_creep

!------------------------------------------------------------------------------------------------------------

      subroutine sens_tang (ic, cgrid, mater, fric, kin, solv, infl, ledg, outpt1, npot, k, iel, wsfix1)
!--purpose: quickly compute senstitivity d F{x,y}/d {xi,eta,phi} with < 1% accuracy
      implicit none
!--subroutine arguments:
      type(t_ic)               :: ic
      type(t_grid)             :: cgrid
      type(t_material)         :: mater
      type(t_friclaw)          :: fric
      type(t_kincns)           :: kin
      type(t_solvers)          :: solv
      type(t_influe)           :: infl
      type(t_leadedge)         :: ledg
      type(t_output)           :: outpt1
      type(t_gridfnc3)         :: wsfix1
      integer                  :: npot, k, iel(npot)
!--local variables:
      integer, parameter :: idebug = 1
      type(t_ic)         :: icloc
      type(t_eldiv)      :: igsloc
      type(t_gridfnc3)   :: wstot, psloc, wsloc, ssloc, usloc
      integer            :: maxtmp, nadh, nslip, nplast, nexter, i, ii, icrp, itry, imeth, it, info1, itgs
      real(kind=8)       :: fxloc, fyloc, mzloc, fxrel, fyrel, mztrue, sens_fxloc, sens_fyloc,          &
                            sens_mzloc, sens_fxtru, sens_fytru, sens_mztru, errx, erry, errz,           &
                            ch_ksi, ch_eta, ch_phi, dcrp, cksi_loc, ceta_loc, cphi_loc, err, conv
      character(len= 4), parameter :: namits(1:3) = (/ 'ItGS', 'ItCG', 'ItGD' /)
      character(len= 3)  :: namcrp

      associate(igs1 => outpt1%igs, ps1 => outpt1%ps, sens => outpt1%sens, facdt => ledg%facdt)

      call gf3_copy_struc(ps1, wstot, 'wstot', .true.)
      call gf3_copy_struc(ps1, psloc, 'psloc', .true.)
      call gf3_copy_struc(ps1, wsloc, 'wsloc', .true.)
      call gf3_copy_struc(ps1, ssloc, 'ssloc', .true.)
      call gf3_copy_struc(ps1, usloc, 'usloc', .true.)
      call eldiv_nullify(igsloc)
      call eldiv_new(igsloc, ps1%grid)

      ! compute original values of total forces fx, fy and spin moment mz

      fxrel = cgrid%dxdy * gf3_sum(AllElm,ps1,1) / (kin%muscal * kin%fntrue)
      fyrel = cgrid%dxdy * gf3_sum(AllElm,ps1,2) / (kin%muscal * kin%fntrue)
      mztrue = cgrid%dxdy * ( - ddot(npot, ps1%vx, 1, cgrid%y, 1)                   &
                              + ddot(npot, ps1%vy, 1, cgrid%x, 1))
      call eldiv_count(igs1, nadh, nslip, nplast, nexter)

      ! compute characteristic sizes of the creepages
      ! Note: these apply to rolling only, T=2, 3, and are off by a factor 2-5 for shifts (?).

      call char_creep(igs1, cgrid, mater, fric%fstat(), kin%fntrue, ch_ksi, ch_eta, ch_phi, idebug)

      if (idebug.ge.2) then
         write(bufout,100) nadh, nslip, fxrel, fyrel, mztrue
         call write_log(1, bufout)
      endif

      ! add contributions of cksi, ceta to vector ws if not done so before
      ! Using work-array wstot instead of input-array wsfix.

      call gf3_copy(AllElm, wsfix1, wstot, ikTANG)

      if (ic%force3.ge.1) then
         do i = 1, k
            ii = iel(i)
            wstot%vx(ii) = wstot%vx(ii) + facdt%vt(ii) * kin%cksi * kin%dq
         enddo
      endif

      if (ic%force3.ge.2) then
         do i = 1, k
            ii = iel(i)
            wstot%vy(ii) = wstot%vy(ii) + facdt%vt(ii) * kin%ceta * kin%dq
         enddo
      endif

      ! copy control digits, set F = 0

      icloc       = ic
      icloc%force3 = 0
      icloc%flow  = ic%flow - 1
      ! icloc%flow  = 5
      maxtmp = solv%maxgs
      itgs   = 0

      ! loop over creepages cksi, ceta, cphi

      do icrp = iin_dksi1, iin_dphi1

         cksi_loc = kin%cksi
         ceta_loc = kin%ceta
         cphi_loc = kin%cphi
         call gf3_copy(AllElm, wstot, wsloc, ikTANG)

         if (icrp.eq.2) then

            ! icrp=2: perturb creepage cksi, add contribution to ws

            namcrp = 'ksi'
            dcrp  = 0.0001d0 * ch_ksi
            cksi_loc = kin%cksi + dcrp
            if (idebug.ge.2) then
               write(bufout,10) namcrp, dcrp, namcrp, cksi_loc
               call write_log(1, bufout)
            endif
   10       format(' perturbation dc',a,'=',g12.4,', c',a,'=',g12.4)

            do i = 1, k
               ii = iel(i)
               wsloc%vx(ii) = wsloc%vx(ii) + facdt%vt(ii)* dcrp * kin%dq
            enddo

         elseif (icrp.eq.3) then

            ! icrp=3: perturb creepage ceta, add contribution to ws

            namcrp = 'eta'
            dcrp  = 0.0001d0 * ch_eta
            ceta_loc = kin%ceta + dcrp
            if (idebug.ge.2) then
               write(bufout,10) namcrp, dcrp, namcrp, ceta_loc
               call write_log(1, bufout)
            endif

            do i = 1, k
               ii = iel(i)
               wsloc%vy(ii) = wsloc%vy(ii) + facdt%vt(ii)* dcrp * kin%dq
            enddo

         elseif (icrp.eq.4) then

            ! icrp=4: perturb creepage cphi, add contribution to ws

            namcrp = 'phi'
            dcrp  = 0.0001d0 * ch_phi
            cphi_loc = kin%cphi + dcrp
            if (idebug.ge.2) then
               write(bufout,10) namcrp, dcrp, namcrp, cphi_loc
               call write_log(1, bufout)
            endif

            do i = 1, k
               ii = iel(i)
               wsloc%vx(ii) = wsloc%vx(ii) - facdt%vt(ii) * dcrp * cgrid%y(ii) * kin%dq
               wsloc%vy(ii) = wsloc%vy(ii) + facdt%vt(ii) * dcrp * cgrid%x(ii) * kin%dq
            enddo

         endif

         ! copy original solution to backup

         call eldiv_copy(igs1, igsloc, ikALL)
         call gf3_copy(AllElm, ps1, psloc, ikALL)
         call gf3_copy(AllElm, outpt1%us, usloc, ikALL)
         call gf3_copy(AllElm, outpt1%ss, ssloc, ikALL)

         do itry = 1, 1

            ! solve system with original/perturbed creepages

            ! if (itry.eq.1) then
            !    solv%mxsens = 30
            ! else
            !    solv%mxsens = 12 - itry
            ! endif

            call tang_solver(npot, icloc, cgrid, mater, fric, kin, solv, infl, ledg, outpt1, k, iel,    &
                         wsloc, .true., info1, imeth, it, err, conv)

            itgs = itgs + it

            ! compute and display results

            fxloc = cgrid%dxdy * gf3_sum(AllElm,ps1,1) / (kin%muscal * kin%fntrue)
            fyloc = cgrid%dxdy * gf3_sum(AllElm,ps1,2) / (kin%muscal * kin%fntrue)
            mzloc = cgrid%dxdy*( - ddot(npot, ps1%vx,1, cgrid%y,1)                  &
                                 + ddot(npot, ps1%vy,1, cgrid%x,1))

            sens_fxloc = (fxloc - fxrel) / dcrp
            sens_fyloc = (fyloc - fyrel) / dcrp
            sens_mzloc = (mzloc - mztrue) / dcrp

            call eldiv_count(igsloc, nadh, nslip, nplast, nexter)

            if (itry.eq.1) then
               sens_fxtru = sens_fxloc
               sens_fytru = sens_fyloc
               sens_mztru = sens_mzloc
               sens(iout_fx1, icrp) = sens_fxtru
               sens(iout_fy1, icrp) = sens_fytru
               sens(iout_mz1, icrp) = sens_mztru

               if (idebug.ge.2) then
                  write(bufout,110) solv%maxgs, nadh, nslip, fxloc, fxloc-fxrel, namcrp, sens_fxloc
                  call write_log(1, bufout)
                  write(bufout,120) fyloc, fyloc-fyrel, namcrp, sens_fyloc
                  call write_log(1, bufout)
                  write(bufout,130) mzloc, mzloc-mztrue, namcrp, sens_mzloc
                  call write_log(1, bufout)
               endif

               if (ic%flow.ge.3) then
                  write(bufout,210) namcrp, nadh, nslip, namits(imeth), it
                  call write_log(1, bufout)
               endif
            else
               errx = (abs(sens_fxloc)-abs(sens_fxtru)) / max(1e-6, abs(sens_fxtru))
               erry = (abs(sens_fyloc)-abs(sens_fytru)) / max(1e-6, abs(sens_fytru))
               errz = (abs(sens_mzloc)-abs(sens_mztru)) / max(1e-6, abs(sens_mztru))

               if (idebug.ge.2) then
                  write(bufout,110) solv%maxgs, nadh, nslip, fxloc, fxloc-fxrel, namcrp,        &
                                        sens_fxloc, 100.*errx
                  call write_log(1, bufout)
                  write(bufout,120) fyloc, fyloc-fyrel, namcrp, sens_fyloc, 100.*erry
                  call write_log(1, bufout)
                  write(bufout,130) mzloc, mzloc-mztrue, namcrp, sens_mzloc, 100.*errz
                  call write_log(1, bufout)
               endif
            endif

            ! restore original solution from backup

            call eldiv_copy(igsloc, igs1, ikALL)
            call gf3_copy(AllElm, psloc, ps1, ikALL)
            call gf3_copy(AllElm, usloc, outpt1%us, ikALL)
            call gf3_copy(AllElm, ssloc, outpt1%ss, ikALL)

         enddo ! itry

 100     format(' original: H, S=',i5,',',i5,', fxtru=',f10.6,:, ', fytru =',f10.6,', mztrue=',f10.3)
 110     format(' maxgs=',i2,': H, S=',i5,',',i5,', fxrel=', f10.6,', dfxrel=',f10.6,', dfx/d',         &
                a,'=',g12.4,:, ', err=',f7.3,' %')
 120     format(29x,'fyrel=',f10.6,', dfyrel=',f10.6,', dfy/d',a,'=', g12.4,:,', err=',f7.3,' %')
 130     format(29x,'mzloc=',f10.3,', dmzloc=',f10.3,', dmz/d',a,'=', g12.4,:,', err=',f7.3,' %')
 210     format(6x,'Tang: perturbed c',a,8x,': A, S=',2i7,'  ',a,'=',i6)

      enddo ! icrp (cksi/ceta/cphi)

      ! restore original value of maxgs

      solv%maxgs = maxtmp

      call gf3_destroy(usloc)
      call gf3_destroy(ssloc)
      call gf3_destroy(wsloc)
      call gf3_destroy(psloc)
      call gf3_destroy(wstot)
      call eldiv_destroy(igsloc)
      end associate

      end subroutine sens_tang

!------------------------------------------------------------------------------------------------------------

      subroutine fastsim_params(ic, mater, fric, kin, cgrid, igs, ws, mus)
!--purpose: compute friction and slope reduction parameters for the modified fastsim algorithm
      implicit none
!--subroutine arguments:
      type(t_ic)               :: ic
      type(t_material)         :: mater
      type(t_friclaw)          :: fric
      type(t_kincns)           :: kin
      type(t_grid)             :: cgrid
      type(t_eldiv)            :: igs
      type(t_gridfnc3)         :: ws, mus
!--local variables:
      logical,       parameter :: use_spiryagin = .true.
      integer                  :: ix, iy, ii
      real(kind=8)             :: a_equiv, b_equiv, aob, area, cx, cy, cz
      real(kind=8)             :: nu_phi, wrel_tot, wabs_tot, f_eff, f_tot, epsil

      ! Compute semi-axes of equivalent ellipse

      call equiv_ellipse(igs, a_equiv, b_equiv)

      ! determine Kalker creep coefficients for the equivalent ellipse

      aob = a_equiv / b_equiv
      call linrol(mater%nu, aob, cx, cy, cz)

      if (use_spiryagin) then

         ! compute total creepage according to [Spiryagin2013], used for slope reduction

         nu_phi = kin%ceta + kin%cphi * a_equiv
         if (abs(nu_phi).lt.abs(kin%ceta)) nu_phi = kin%ceta
         wrel_tot = sqrt(kin%cksi**2 + nu_phi**2)

         ! compute total absolute rigid slip according to [Spiryagin], used for velocity dependence

         wabs_tot = kin%veloc * sqrt(kin%cksi**2 + kin%ceta**2)

      else

         ! compute total creepage and total absolute rigid slip using 2-norm of right hand side ws
         ! note that ws == -shift distance

         wrel_tot = gf3_rms(AllInt, ws, ikTANG) / kin%dq
         wabs_tot = kin%veloc * wrel_tot

      endif

      ! compute the actual coefficient of friction corresponding to wabs_tot

      if (fric%frclaw_eff.eq.0) then

         call gf3_copy_xdir(AllElm, fric%nvf, fric%fstat_arr, mus, ikXDIR)
         f_eff = sum(fric%fstat_arr(1:fric%nvf)) / real(fric%nvf)

      elseif (fric%frclaw_eff.ge.2 .and. fric%frclaw_eff.le.4) then

         f_tot = 0d0
         do iy = 1, cgrid%ny
            f_eff = fric%mu_steady(wabs_tot, iy)
            do ix = 1, cgrid%nx
               ii = ix + (iy-1) * cgrid%nx
               mus%vx(ii) = f_eff
            enddo
            f_tot = f_tot + f_eff
         enddo
         f_eff = f_tot / real(cgrid%ny)

         if (ic%flow.ge.5) then
            write(bufout,*) 'L=',fric%frclaw_eff,', wabs_tot=',wabs_tot,': f_eff=',f_eff
            call write_log(1, bufout)
         endif

      else

         write(bufout, '(a,i2,a)') ' INTERNAL ERROR in FASTSIM: no support for L=', fric%frclaw_eff, '.'
         call write_log(1, bufout)
         call abort_run()

      endif

      ! compute slope reduction parameter epsilon

      area  = pi * a_equiv * b_equiv
      epsil = mater%ga * area * mater%k0_mf * cx * wrel_tot / (4d0 * f_eff * kin%fntrue)

      if (ic%flow.ge.5) then
         write(bufout,123) ' eps=',epsil,' =',mater%ga,' *',area,' *',mater%k0_mf, ' *',cx,     &
                ' *',wrel_tot,' / (4 *', f_eff,' *',kin%fntrue,' )'
 123     format(a,f8.3,a,f7.0,a,f6.1,a,f6.3,a,f6.3,a,f6.1,a,f6.3,a,f7.0,a)
         call write_log(1, bufout)
      endif

      ! compute slope reduction factor k

      mater%k_eff = mater%k0_mf * (mater%alfamf + (1d0 - mater%alfamf) / (1d0 + mater%betamf * epsil))

      if (ic%flow.ge.5) then
         write(bufout,133) ' k_eff=', mater%k_eff,' =', mater%k0_mf,' * (',mater%alfamf,' + ',          &
                1d0-mater%alfamf,' / (1d0 + ',mater%betamf,' *',epsil,') )'
 133     format(a,f8.5,a,f8.5,a,f6.3,a,f6.3,a,f6.3,a,f8.3,a)
         call write_log(1, bufout)
      endif

      end subroutine fastsim_params

!------------------------------------------------------------------------------------------------------------

      subroutine modif_fastsim(ic, cgrid, mater, fric, kin, ws, outpt)
!--purpose: solve the tangential problem using the modified FASTSIM approach, i.e. using the Winkler
!           foundation material model, with slope reduction according to [Spiryagin2013]
      implicit none
!--subroutine arguments:
      type(t_ic)               :: ic
      type(t_grid)             :: cgrid
      type(t_material)         :: mater
      type(t_friclaw)          :: fric
      type(t_kincns)           :: kin
      type(t_gridfnc3)         :: ws
      type(t_output)           :: outpt
!--local variables :
      integer      :: ii, ix, ixsta, ixinc, ixend, iy, ixdbg, iydbg
      logical      :: is_ssrol
      real(kind=8) :: ptabs

      associate(mx    => cgrid%nx,  my    => cgrid%ny,  flx   => mater%flx, k_eff => mater%k_eff,       &
                igs   => outpt%igs, mus   => outpt%mus, ps    => outpt%ps,  us    => outpt%us,          &
                uv    => outpt%uv,  ss    => outpt%ss)

      is_ssrol = ic%tang.eq.3

      ! set row/point for which debug-output is requested

      if (.false.) then
         ixdbg = mx
         iydbg = my/2
      else
         ixdbg = -1
         iydbg = -1
      endif

      ! call fastsim_params to set the effective coefficient of friction f_eff and slope reduction k_eff

      call fastsim_params(ic, mater, fric, kin, cgrid, igs, ws, mus)

      ! For transient cases (T<>3), the grid points may be processed in any order, but for steady state
      !    rolling (T=3), they must be processed from the leading to trailing edge
      ! Set loop start/increment/end for processing rows on the basis of the rolling direction (0 or 180 deg)

      if (abs(kin%chi-pi).lt.0.01d0 .or. .not.is_ssrol) then
         ixsta = 1
         ixinc = 1
         ixend = mx
      elseif (abs(kin%chi).lt.0.01d0) then
         ixsta = mx
         ixinc = -1
         ixend = 1
      else
         write(bufout,*) 'ERROR: for steady state rolling with Fastsim',                &
                ' the rolling direction should be 0 or 180deg'
         call write_log(1, bufout)
         call abort_run()
      endif

      if (ixdbg.ge.1 .and. ixdbg.le.mx .and. iydbg.ge.1 .and. iydbg.le.my) then
         ii = ixdbg + mx*(iydbg-1)
         write(bufout,*) 'ws(iidbg)=', ws%vx(ii), ws%vy(ii)
         call write_log(1, bufout)
      endif

      ! loop over all rows of the grid -- independent calculations

      do iy = 1, my

         ! process points of row iy in the appropriate order

         do ix = ixsta, ixend, ixinc
            ii = ix + (iy-1)*mx

            ! In steady state rolling: set U' by copying U of the previous point -- assuming dq == dx

            if (is_ssrol .and. ix.eq.ixsta) then
               uv%vx(ii) = 0d0
               uv%vy(ii) = 0d0
            elseif (is_ssrol) then
               uv%vx(ii) = us%vx(ii-ixinc)
               uv%vy(ii) = us%vy(ii-ixinc)
            endif

            if (igs%el(ii).le.Exter) then

               ! Exterior elements: no traction, no slip

               ps%vx(ii) = 0d0
               ps%vy(ii) = 0d0
               us%vx(ii) = 0d0
               us%vy(ii) = 0d0
               ss%vx(ii) = 0d0
               ss%vy(ii) = 0d0
            else

               ! Interior elements: time-step equation, implicit Euler:
               !    U(i) - U'(i) = S(i) - W(i) = dq*s(i) - dq*w(i)

               ! First assume adhesion, s(i)=0

               igs%el(ii)   = Adhes

               ss%vx(ii) = 0d0
               ss%vy(ii) = 0d0

               us%vx(ii) = uv%vx(ii) - ws%vx(ii)
               us%vy(ii) = uv%vy(ii) - ws%vy(ii)

               ps%vx(ii) = us%vx(ii) * k_eff / flx(1)
               ps%vy(ii) = us%vy(ii) * k_eff / flx(2)

               ! If traction bound is exceeded

               ptabs = sqrt(ps%vx(ii)**2 + ps%vy(ii)**2)
               if (ptabs.gt.mus%vt(ii)*ps%vn(ii)) then

                  ! Solve equations for slip in element i

                  igs%el(ii)   = Slip

                  ! The direction of the tractions ( // W ) is ok, scale tractions

                  ps%vx(ii) = ps%vx(ii) * mus%vt(ii)*ps%vn(ii)/ptabs
                  ps%vy(ii) = ps%vy(ii) * mus%vt(ii)*ps%vn(ii)/ptabs

                  us%vx(ii) = flx(1) * ps%vx(ii) / k_eff
                  us%vy(ii) = flx(2) * ps%vy(ii) / k_eff

                  ! Compute slip: dq*s(i) = U(i) - U'(i) + dq*w(i)

                  ss%vx(ii) = us%vx(ii) - uv%vx(ii) + ws%vx(ii)
                  ss%vy(ii) = us%vy(ii) - uv%vy(ii) + ws%vy(ii)
               endif
            endif

            if (iy.eq.iydbg) then
               if (ix.eq.ixsta) then
                  write(bufout,122)
                  call write_log(1, bufout)
               endif
               if (igs%el(ii).ge.Adhes) then
                   write(bufout,123) ix, igs%el(ii),                                    &
                      ps%vx(ii), us%vx(ii), uv%vx(ii), ws%vx(ii), ss%vx(ii),            &
                      ps%vy(ii), us%vy(ii), uv%vy(ii), ws%vy(ii), ss%vy(ii)
                   call write_log(2, bufout)
               endif
 122           format(6x,'igs',5x,'ps',10x,'us',10x,'uv',10x,'ws',10x,'ss')
 123           format('ix=',i3,i3,5f12.6,:/,9x,5f12.6)
            endif

         ! end loop over all elements

         enddo
      enddo

      end associate
      end subroutine modif_fastsim

!------------------------------------------------------------------------------------------------------------

      subroutine modif_fastrip(ic, cgrid, mater, fric, kin, ws, outpt)
!--purpose: solve the tangential problem using the modified FASTRIP approach according to [Sichani2016b]
!           with slope reduction according to [Spiryagin2013]
      implicit none
!--subroutine arguments:
      type(t_ic)               :: ic
      type(t_grid)             :: cgrid
      type(t_material)         :: mater
      type(t_friclaw)          :: fric
      type(t_kincns)           :: kin
      type(t_gridfnc3)         :: ws
      type(t_output)           :: outpt
!--local variables :
      integer      :: ii, ix, ixsta, ixinc, ixend, iy
      real(kind=8) :: nu, mu, aa, bb, aob, cxx, cyy, cyz, cksi_p, ceta_p, cphi_p, xsta, xend, aa_j,     &
                      pmax, ptabs, y_j, dd_j, kappa, lambda, lambdap, sq_cnt, sq_adh

      associate(mx    => cgrid%nx,  my    => cgrid%ny,  flx   => mater%flx, k_eff => mater%k_eff,       &
                igs   => outpt%igs, mus   => outpt%mus, ps    => outpt%ps,  us    => outpt%us,          &
                uv    => outpt%uv,  ss    => outpt%ss)

      if (ic%tang.ne.3) then
         call write_log(' ERROR: FaStrip is supported only for steady state rolling (T=3).')
         call abort_run()
      endif
      if (abs(kin%chi-pi).ge.0.01d0 .and. abs(kin%chi).Gt.0.01d0) then
         call write_log(' ERROR: FaStrip requires the rolling direction be 0 or 180deg.')
         call abort_run()
      endif

      ! Compute semi-axes of equivalent ellipse

      call equiv_ellipse(igs, aa, bb)

      ! basic implementation for Hertzian cases

      ! 1. inputs for Hertzian case

      nu   = mater%nu
      mu   = fric%fstat()
      pmax = 1.5d0 * kin%fntrue / (pi * aa * bb)

      ! 2. preparations: Kalker coefficients, flexibilities

      aob  = aa / bb
      call linrol(nu, aob, cxx, cyy, cyz)

      ! non-dimensional creepages for strip theory

      cksi_p = -mater%ga * 4d0 * (1d0-nu) * cxx * kin%cksi / (2d0 * mu * pmax * pi**2)
      ceta_p = -mater%ga * 4d0 * (1d0-nu) * cyy * kin%ceta / (2d0 * mu * pmax * pi**2)
      cphi_p = -mater%ga * 3d0 * sqrt(aa*bb) * cyz * kin%cphi / (2d0 * mu * pmax * pi)

      ! call fastsim_params to set the effective coefficient of friction f_eff and slope reduction k_eff

      call fastsim_params(ic, mater, fric, kin, cgrid, igs, ws, mus)

      ! Points are processed from the leading to trailing edge
      ! Set loop start/increment/end on the basis of the rolling direction (0 or 180 deg)

      if (abs(kin%chi-pi).lt.0.01d0) then
         ixsta = 1
         ixinc = 1
         ixend = mx
      else
         ixsta = mx
         ixinc = -1
         ixend = 1
      endif

      ! loop over all rows of the grid -- independent calculations

      do iy = 1, my

         ii   = 1 + (iy-1)*mx
         y_j  = cgrid%y(ii)

         ! 3. determine contact length for strip y_j: [xsta(j), xend(j)]

            xsta   =  99d9
            xend   = -99d9
            do ix = 1, mx
               ii = ix + (iy-1)*mx
               if (igs%el(ii).ge.Adhes) then
               xsta  = min(xsta, cgrid%x(ii)-0.5d0*cgrid%dx)
               xend  = max(xend, cgrid%x(ii)+0.5d0*cgrid%dx)
               endif
            enddo
            aa_j = max(1d-10, 0.5d0 * (xend - xsta))

         ! 4. Fastsim, 3 flexibilities. process points ix within row iy in the appropriate order

         do ix = ixsta, ixend, ixinc
            ii = ix + (iy-1)*mx

            ! Set U' by copying U of the previous point -- assuming dq == dx

            if (ix.eq.ixsta) then
               uv%vx(ii) = 0d0
               uv%vy(ii) = 0d0
            else
               uv%vx(ii) = us%vx(ii-ixinc)
               uv%vy(ii) = us%vy(ii-ixinc)
            endif

            if (igs%el(ii).le.Exter) then

               ! Exterior elements: no traction, no slip

               ps%vx(ii) = 0d0
               ps%vy(ii) = 0d0
               us%vx(ii) = 0d0
               us%vy(ii) = 0d0
               ss%vx(ii) = 0d0
               ss%vy(ii) = 0d0
            else

               ! Interior elements: time-step equation, implicit Euler:
               !    U(i) - U'(i) = S(i) - W(i) = dq*s(i) - dq*w(i)

               ! First assume adhesion, s(i)=0

               igs%el(ii)   = Adhes

               ss%vx(ii) = 0d0
               ss%vy(ii) = 0d0

               us%vx(ii) = uv%vx(ii) - ws%vx(ii)
               us%vy(ii) = uv%vy(ii) - ws%vy(ii)

               ps%vx(ii) = us%vx(ii) * k_eff / flx(1)
               ps%vy(ii) = us%vy(ii) * k_eff / flx(2)

               ! Check traction bound for element I

               ptabs = sqrt(ps%vx(ii)**2 + ps%vy(ii)**2)

               if (ptabs.gt.mu*ps%vn(ii)) then

                  ! Solve equations for slip in element i

                  igs%el(ii)   = Slip

                  ! The direction of the tractions ( // W ) is ok, scale tractions

                  ps%vx(ii) = ps%vx(ii) * mu * ps%vn(ii) / ptabs
                  ps%vy(ii) = ps%vy(ii) * mu * ps%vn(ii) / ptabs

                  us%vx(ii) = flx(1) * ps%vx(ii) / k_eff
                  us%vy(ii) = flx(2) * ps%vy(ii) / k_eff

                  ! Compute slip: dq*s(i) = U(i) - U'(i) + dq*w(i)

                  ss%vx(ii) = us%vx(ii) - uv%vx(ii) + ws%vx(ii)
                  ss%vy(ii) = us%vy(ii) - uv%vy(ii) + ws%vy(ii)
               endif
            endif
         enddo ! ix, Fastsim 3 flexibilities

         ! 7. compute length of slip area in strip theory

         if (abs(cphi_p).ge.1d0-1d-6) then
            dd_j = aa_j
         else
            dd_j = (sqrt( ceta_p**2 + (1d0 - cphi_p**2) * (cksi_p - cphi_p * y_j / aa)**2 ) +           &
                                                ceta_p * cphi_p) / (1d0 - cphi_p**2) * aa / (1d0 - nu)
         endif

         ! 8. compute strip theory stress scaling factors

         kappa   = (cksi_p - cphi_p * y_j / aa) /                                                       &
                   sqrt( max(1d-20, ( cksi_p - cphi_p * y_j / aa)**2 + (ceta_p + cphi_p * dd_j / aa)**2 ) )

         lambda  = (ceta_p + cphi_p * dd_j / aa) /                                                      &
                   sqrt( max(1d-20, ( cksi_p - cphi_p * y_j / aa)**2 + (ceta_p + cphi_p * dd_j / aa)**2 ) )

         lambdap = lambda - cphi_p

         ! 9. compute strip theory tractions within strip-theory's adhesion area

         do ix = 1, mx
            ii = ix + (iy-1)*mx
            if (igs%el(ii).ge.Adhes .and. cgrid%x(ii).gt.-aa_j+2d0*dd_j) then
               igs%el(ii) = Adhes

               sq_cnt = sqrt( max(0d0, aa_j**2 - cgrid%x(ii)**2 ) )
               sq_adh = sqrt( max(0d0, (aa_j - dd_j)**2 - (cgrid%x(ii) - dd_j)**2 ) )

               ps%vx(ii) = mu * pmax * (kappa  * sq_cnt - kappa   * sq_adh) / aa
               ps%vy(ii) = mu * pmax * (lambda * sq_cnt - lambdap * sq_adh) / aa

            if (isnan(ps%vx(ii)) .or. isnan(ps%vy(ii))) then
                  write(bufout,'(2(a,i4),2(a,g12.4))') ' iy=',iy,', ix=',ix,': px=',ps%vx(ii),          &
                        ', py=',ps%vy(ii)
               call write_log(1, bufout)
            endif
            endif

         enddo

      enddo ! row iy

      end associate
      end subroutine modif_fastrip

!------------------------------------------------------------------------------------------------------------

      include 'gdsteady.f90'

!------------------------------------------------------------------------------------------------------------

      subroutine tangcg(ic, npot, maxcg, eps, is_sens, ws, infl, outpt, itcg, err)
!--purpose: Nonlinear Conjugate Gradient type solver for the system of equations and constraints in TANG with
!           fixed creepages. Variant 9 including inner iterations and FFT-based preconditioner. Matrix-free
!           implementation, using routine VecAijPj for evaluating "A * v", with A represented by cs.
      implicit none
!--subroutine arguments:
      type(t_ic),       intent(in)      :: ic
      integer,          intent(in)      :: npot, maxcg
      logical,          intent(in)      :: is_sens
      real(kind=8),     intent(in)      :: eps
      type(t_gridfnc3), intent(in)      :: ws
      type(t_influe),   intent(in)      :: infl
      type(t_output)                    :: outpt
      integer,          intent(out)     :: itcg
      real(kind=8),     intent(out)     :: err
!--local variables:
      integer, parameter :: num_inn = 4
      integer, parameter :: itdbg = 500
      real(kind=8), parameter :: small = 1d-6
      logical            :: use_fftprec, lchanged
      type(t_gridfnc3)   :: g, n, t, r, z, v, q, pold, ps0, tmp
      integer            :: mx, my, nadh, nslip, nplast, nexter, imvp, it_inn, ii, iidbg, iidbgx,       &
                            iidbgy, ic_nmdbg
      real(kind=8)       :: alpha, beta, rv, zq, vq, snrm, perp, ga, ptabs, facnel, dif, dif1, difid,   &
                            difinn, ptang, trsinn, conv

      associate(cs   => infl%cs,  ms   => infl%ms, igs  => outpt%igs, mu   => outpt%mus,                &
                ps   => outpt%ps, ss   => outpt%ss)

      mx = ps%grid%nx
      my = ps%grid%ny
      ga = cs%ga
      ic_nmdbg = ic%x_nmdbg
      iidbg = min(npot, 196)

#if defined WITH_MKLFFT
      use_fftprec = (ic%mater.ne.1)
#else
      use_fftprec = .false.
#endif

      if (ic%flow.ge.6 .and. use_fftprec) then
         write(bufout,'(9x,a,i1,3a)') ' TangCG(',num_inn,'):',                          &
                        ' using',' FFT-based preconditioner'
         call write_log(1, bufout)
      elseif (ic%flow.ge.6) then
         write(bufout,'(9x,a,i1,3a)') ' TangCG(',num_inn,'):',                          &
                        ' not using',' FFT-based preconditioner'
         call write_log(1, bufout)
      endif

      ! Initialize grid functions.
      ! Note: they should all be zero in the exterior area at all times.
      !       g is the traction bound, n, t normal and tangential directions,
      !       r is the residual, z the preconditioned residual,
      !       v is the search direction, q = A * v.

      call gf3_copy_struc(ps, g   , 'tangcg:g'   , .true.)
      call gf3_copy_struc(ps, n   , 'tangcg:n'   , .true.)
      call gf3_copy_struc(ps, t   , 'tangcg:t'   , .true.)
      call gf3_copy_struc(ps, r   , 'tangcg:r'   , .true.)
      call gf3_copy_struc(ps, z   , 'tangcg:z'   , .true.)
      call gf3_copy_struc(ps, v   , 'tangcg:v'   , .true.)
      call gf3_copy_struc(ps, q   , 'tangcg:q'   , .true.)
      call gf3_copy_struc(ps, pold, 'tangcg:pold', .true.)

      ! save initial estimate for modified stopping criterion

      if (is_sens) then
         call gf3_copy_struc(ps, ps0, 'tangcg:ps0', .true.)
         call gf3_copy_struc(ps, tmp, 'tangcg:tmp', .true.)
         call gf3_copy(AllElm, ps, ps0, ikTANG)
      endif

      ! Check requirement on influence coefficients (diagonal scaling)

      if (max(abs(cs%cf(0,0,1,2)), abs(cs%cf(0,0,2,1))).gt.                             &
             1e-8*min(cs%cf(0,0,1,1), cs%cf(0,0,2,2))) then
         write(bufout,111) 'Non-zero off-diagonals in A_00!',                           &
                    'A_00 = ',cs%cf(0,0,1,1), cs%cf(0,0,1,2),                           &
                    '       ',cs%cf(0,0,2,1), cs%cf(0,0,2,2)
 111     format(a,/,a,2g12.4,/,a,2g12.4)
         call write_log(3, bufout)
         call abort_run()
      endif

      ! Make preconditioner of diagonal blocks Axx and Ayy

#if defined WITH_MKLFFT
      if (use_fftprec) then
         call fft_makeprec(1, cs, 1, ms)
         call fft_makeprec(2, cs, 2, ms)
      endif
#endif

      ! count number of elements in H, S for switching off preconditioner

      call eldiv_count(igs, nadh, nslip, nplast, nexter)

      ! correction factor for rms-norm: pot.con versus actual contact

      facnel = sqrt(real(npot)/real(nadh+nslip+nplast))

      ! Compute traction-bound at each element ii - compatibility with previous implementation

      do ii = 1, npot
         g%vt(ii) = mu%vt(ii) * ps%vn(ii)
      enddo

      ! This routine uses a relative stop-criterion on the updates dif == ||pk - pk-1||_rms, which are
      !    meaningful to the user.

      itcg   = 0
      imvp   = 0
      dif    = 2d0
      difid  = 1d0
      it_inn = 0
      difinn = 0d0
      lchanged = .false.

      ! 1. Compute argument \theta, primarily for the elements in slip
      ! also: normal (n) and tangential (t) direction vectors

      do ii = 1, npot
         n%vn(ii) = atan2(ps%vy(ii), ps%vx(ii))
         n%vx(ii) = cos(n%vn(ii))
         n%vy(ii) = sin(n%vn(ii))
         t%vx(ii) = -n%vy(ii)
         t%vy(ii) =  n%vx(ii)
      enddo

      ! 1. Compute the slip ss = A * ps + ws.

      call gf3_set(AllElm, 0d0, ss, ikTANG)
      call VecAijPj(igs, AllInt, ss, ikTANG, ps, jkTANG, cs)
      call gf3_axpy(AllInt, 1.0d0, ws, ss, ikTANG)
      imvp = imvp + 1

      ! 2. Compute residual r
      !     - adhesion: r = -ss
      !     - slip:     r = t * t' * -ss

      call gf3_copy(AllElm, ss, r, ikTANG)
      call gf3_scal(AllElm, -1.0d0, r, ikTANG)

      do ii = 1, npot
         if (igs%el(ii).eq.Slip) then
            perp = t%vx(ii)*r%vx(ii) + t%vy(ii)*r%vy(ii)
            r%vx(ii) = perp * t%vx(ii)
            r%vy(ii) = perp * t%vy(ii)
         endif
      enddo

      ! WHILE ( (Changes to el.div OR dif > difid)  AND  itcg < maxcg ) do

      do while ((lchanged .or. dif.gt.difid) .and. itcg.lt.maxcg)

         ! increment iteration counter "it"

         itcg = itcg + 1
         it_inn = it_inn + 1

         ! Disable preconditioner near full sliding, >75% slip
         ! If so, restart CG-iteration
         !!! temporarily also in case of convergence problems with large dy / dx:

         if (use_fftprec) then
            if (nslip+nplast.gt.3*nadh .and. my.gt.1) use_fftprec = .false.
            if (itcg.ge.maxcg/2) use_fftprec = .false.
            if (.not.use_fftprec) lchanged = .true.
         endif

         ! 3. Compute preconditioned residual z = M * res

         if (use_fftprec) then
            call VecAijPj(igs, AllInt, z, ikXDIR, r, jkXDIR, ms)
            call VecAijPj(igs, AllInt, z, ikYDIR, r, jkYDIR, ms)
            imvp = imvp + 1
         else
            call gf3_copy(AllElm, r, z, ikTANG)
         endif

         ! 3. Preconditioning: apply diagonal scaling z = D^-1 * z using diagonal of reduced Jacobian
         !     - adhesion: d = diag(A)
         !     - slip:     d = diag(A) - (s,n) / g

         do ii = 1, npot
            if (igs%el(ii).eq.Slip) then
               snrm = -ss%vx(ii)*n%vx(ii) - ss%vy(ii)*n%vy(ii)
               z%vx(ii) = z%vx(ii) / (cs%cf(0,0,1,1) + ga*snrm/g%vt(ii))
               z%vy(ii) = z%vy(ii) / (cs%cf(0,0,2,2) + ga*snrm/g%vt(ii))
            else
               z%vx(ii) = z%vx(ii) / cs%cf(0,0,1,1)
               z%vy(ii) = z%vy(ii) / cs%cf(0,0,2,2)
            endif
         enddo

         ! 3. Preconditioning: project onto tangential direction t
         !        z := t * t' * z

         do ii = 1, npot
            if (igs%el(ii).eq.Slip) then
               perp = t%vx(ii)*z%vx(ii) + t%vy(ii)*z%vy(ii)
               z%vx(ii) = perp * t%vx(ii)
               z%vy(ii) = perp * t%vy(ii)
            endif
         enddo

         ! 4. Set (primary) search direction v

         if (itcg.le.1 .or. lchanged) then

            ! 4.2: restart CG using v = z, steepest descent direction

            call gf3_copy(AllElm, z, v, ikTANG)
         else

            ! 4.1: construct search direction conjugate to previous one
            !      beta: (r^k+2, v^k) = 0

            zq = gf3_dot(AllElm, z, q, ikTANG)
            vq = gf3_dot(AllElm, v, q, ikTANG)
            if (abs(zq).lt.1d-60 .or. abs(vq).lt.small*abs(zq)) then
               beta = 0d0
            else
               beta = - zq / vq
            endif
            if (beta.lt.0d0 .and. ic_nmdbg.ge.5) then
               ! beta<0 does not seem to harm convergence or cause any other problems
               write(bufout,'(a,i0,a,f14.9,2(a,g12.4))') ' it=',itcg,': beta= ',beta,' =',-zq,' /',vq
               call write_log(1,bufout)
            endif
            if (.false.) beta = max(0d0, beta)

               !   v := z + beta*vold

            call gf3_scal(AllElm, beta, v, ikTANG)
            call gf3_axpy(AllElm, 1.0d0, z, v, ikTANG)
         endif

         ! 6. compute change in slip q1 = A v
         !    slip area: project onto tangential direction t

         call VecAijPj(igs, AllInt, q, ikTANG, v, jkTANG, cs)
         imvp = imvp + 1

         do ii = 1, npot
            if (igs%el(ii).eq.Slip) then
               perp = t%vx(ii)*q%vx(ii) + t%vy(ii)*q%vy(ii)
               q%vx(ii) = perp * t%vx(ii)
               q%vy(ii) = perp * t%vy(ii)
            endif
         enddo

         ! 7. add effect "q2" of changing tangential direction t
         !        q2 = - v * (s,n) / g

         do ii = 1, npot
            if (igs%el(ii).eq.Slip) then
               snrm = (n%vx(ii)*ss%vx(ii) + n%vy(ii)*ss%vy(ii)) / g%vt(ii)
               q%vx(ii) = q%vx(ii) - snrm * v%vx(ii)
               q%vy(ii) = q%vy(ii) - snrm * v%vy(ii)
            endif
         enddo

         ! 8. compute step length alpha
         !   --> make new r^k+1 orthogonal to v^k
         !       residual contains two components q1 and q2
         ! Note: if A is truely pos.def. then vq is always >= 0.  vq == 0 iff v == 0.
         !       Avoid division by zero for the special case of a single element where v == 0 occurs.

         rv = gf3_dot(AllElm, r, v, ikTANG)
         vq = gf3_dot(AllElm, v, q, ikTANG)

         if (abs(rv).lt.1d-60) then
            alpha = 0d0
         elseif (abs(vq).lt.small*abs(rv)) then
            alpha = 1d0
            if (ic_nmdbg.ge.15) then
               write(bufout,*) ' it=',itcg,': alpha=',alpha,'=',rv,' / ',vq
               call write_log(1,bufout)
            endif
         else
            alpha = rv / vq
         endif
         if (isnan(alpha)) then
            write(bufout,*) ' it=',itcg,': alpha=',alpha,'=',rv,' / ',vq
            call write_log(1,bufout)
         endif

         if (min(rv,vq).lt.0d0 .and. max(rv,vq).gt.0d0 .and. ic_nmdbg.ge.3) then
            write(bufout,'(a,i0,3(a,es16.8),a)') ' it=',itcg,': alpha=',alpha,' =',rv,' /',vq,' < 0!'
            call write_log(1, bufout)
            do ii = 1, npot
               if (igs%el(ii).gt.Exter) then
                  write(bufout,'(a,i4,6(a,g12.4))') 'ii=',ii,                                           &
                     ': pn=',ps%vn(ii),', px=',ps%vx(ii),', vx=',v%vx(ii),                              &
                     ', qx=', q%vx(ii),', wx=',ws%vx(ii),', rx=',r%vx(ii)
                  call write_log(1, bufout)
                  write(bufout,'( 26x,5(a,g12.4))') 'py=',ps%vy(ii),                                    &
                     ', vy=',v%vy(ii),', qy=',q%vy(ii),', wy=',ws%vy(ii),', ry=',r%vy(ii)
                  call write_log(1, bufout)
               endif
            enddo
         endif

         if (ic_nmdbg.ge.15 .or. itcg.eq.-1) then
            ii = iidbg
            write(bufout,'(a,i4,a,g12.4)') 'iter=',itcg,', alpha=',alpha
            call write_log(1, bufout)
            write(bufout,'(a,i4,6(a,g12.4))') 'ii=',ii,                                 &
                  ': sx=',ss%vx(ii),', rx=',r%vx(ii),', zx=',z%vx(ii),                  &
                  ', vx=',v%vx(ii),', qx=',q%vx(ii),', px=',ps%vx(ii)
            call write_log(1, bufout)
            write(bufout,'(  7x,6(a,g12.4))')                                           &
                  ': sy=',ss%vy(ii),', ry=',r%vy(ii),', zy=',z%vy(ii),                  &
                  ', vy=',v%vy(ii),', qy=',q%vy(ii),', py=',ps%vy(ii)
            call write_log(1, bufout)
         endif

         ! 9. Update estimate of tangential tractions ps:
         !    \tilde{ps}(it) = ps(it-1) + alpha * v,

         call gf3_copy(AllElm, ps, pold, ikTANG)
         call gf3_axpy(AllInt, alpha, v, ps, ikTANG)

         ! 9.2 In Slip: put tractions at the traction bound.

         do ii = 1, npot
            if (igs%el(ii).eq.Slip) then
               ptabs = sqrt(ps%vx(ii)**2 + ps%vy(ii)**2)
               ps%vx(ii) = ps%vx(ii) * g%vt(ii)/ptabs
               ps%vy(ii) = ps%vy(ii) * g%vt(ii)/ptabs
            endif
         enddo

         ! 11.3 Relative stop-criterion on updates |pk - pk-1|:
         !    dif   = |pk - pk-1|_rms \approx alpha * |v|_rms
         !    difid = eps * |pk|_rms

         dif = alpha * facnel * gf3_rms(AllElm, v, ikTANG)

         if (ic_nmdbg.ge.5 .and. itcg.eq.maxcg-5) then
            iidbgx = idamax(npot, v%vx(1:),1)
            iidbgy = idamax(npot, v%vy(1:),1)
            write(bufout,*) 'max vx=',v%vx(iidbgx),' at ii=',iidbgx,                    &
                ', max vy=',v%vy(iidbgy),' at ii=',iidbgy
            call write_log(1, bufout)
            if (abs(v%vx(iidbgx)).gt.abs(v%vy(iidbgy))) then
               iidbg = iidbgx
            else
               iidbg = iidbgy
            endif
            ic_nmdbg = 15
         endif

         difinn = difinn + dif

         if (itcg.le.5 .or. mod(itcg,5).eq.1) then
            if (.not.is_sens) then
               ptang = facnel * gf3_rms(AllElm, ps, ikTANG)
            else
               call gf3_copy(AllElm, ps, tmp, ikTANG)
               call gf3_axpy(AllElm, -1d0, ps0, tmp, ikTANG)
               ptang = facnel * gf3_rms(AllElm, tmp, ikTANG)
            endif
            difid = eps * max(1d-6, ptang)
            trsinn = 0.01 * ptang
         endif

         ! 10. The inner iteration is stopped when either:
         !      a. max. number of inner iterations reached;
         !      b. overall convergence may have been reached;
         !      c. accumulated updates difinn grow larger than threshold trsinn
         !     In these cases, the outer iteration takes over.

         if (it_inn.ge.num_inn .or. dif.le.difid .or. difinn.gt.trsinn) then

            if (ic%flow.ge.5 .and. ic_nmdbg.ge.1 .and. .false.) then
               write(bufout,'(2(a,i2),2(a,g12.4))') ' it ',itcg,                        &
                         ': inn=',it_inn,': cum=',difinn,', trs=',trsinn
               call write_log(1, bufout)
            endif

            ! 11. Check traction bound and sense of slip, adjust element division accordingly

            lchanged = .false.

            ! 11.1 In Adhesion: check for tractions exceeding traction bound, move elements to the Slip area.

            do ii = 1, npot
               if (igs%el(ii).eq.Adhes) then
                  ptabs = sqrt(ps%vx(ii)**2 + ps%vy(ii)**2)
                  if (ptabs.gt.g%vt(ii)) then
                     if (ic_nmdbg.ge.5 .or. itcg.ge.maxcg-5) then
                        write(bufout,*) 'it=',itcg,': moving element ii=',ii,' to slip area'
                        call write_log(1,bufout)
                     endif
                     igs%el(ii) = Slip
                     ps%vx(ii) = ps%vx(ii) * g%vt(ii)/ptabs
                     ps%vy(ii) = ps%vy(ii) * g%vt(ii)/ptabs
                     lchanged = .true.
                     if (ic_nmdbg.ge.15 .and. ii.eq.iidbg) then
                        write(bufout,'(a,i4,2(a,g16.8),2(a,g12.4))')                    &
                           'ii=',ii, ': pmax=',g%vt(ii),', ptabs=',                     &
                           ptabs,' -> px=',ps%vx(ii),', py=',ps%vy(ii)
                        call write_log(1, bufout)
                     endif
                  endif
               endif
            enddo
            if (lchanged .and. ic_nmdbg.ge.9) call wrigs(igs, .false., 0d0)

            ! 11.2 Compute the slip ss = A * ps + ws.

            call VecAijPj(igs, AllInt, ss, ikTANG, ps, jkTANG, cs)
            call gf3_axpy(AllInt, 1.0d0, ws, ss, ikTANG)
            imvp = imvp + 1

            ! 11.2 In Slip: check direction of the slip, move elements to
            !      the Adhesion area where needed.

            do ii = 1, npot
               if (igs%el(ii).eq.Slip) then
                  snrm = ss%vx(ii)*ps%vx(ii) + ss%vy(ii)*ps%vy(ii)
                  if (snrm.gt.0d0) then
                     igs%el(ii) = Adhes
                     if (ic_nmdbg.ge.5 .or. itcg.ge.maxcg-5) then
                        write(bufout,*) 'it=',itcg,': moving element ii=',ii,' to adhesion area'
                        call write_log(1,bufout)
                     endif
                     if (ic_nmdbg.ge.15 .and. ii.eq.iidbg) then
                        write(bufout,'(a,i4,5(a,g14.6))') 'ii=',ii,                     &
                           ': wx=',ws%vx(ii),', wy=',ws%vy(ii),                         &
                           ', sx=',ss%vx(ii),', sy=',ss%vy(ii),                         &
                           ', snrm=',snrm
                        call write_log(1, bufout)
                     endif
                     lchanged = .true.
                  endif
               endif
            enddo

            ! 11.3 compute stop-criterion again for modified ps

            call gf3_axpy(AllElm, -1d0, ps, pold, ikTANG)
            dif = facnel * gf3_rms(AllElm, pold, ikTANG)

            ! 11.4 constraints checked, ss formed: restart inner iteration

            it_inn = 0
            difinn = 0d0

            if (lchanged) call eldiv_count(igs, nadh, nslip, nplast, nexter)

         endif ! check constraints

         ! 11.5 Print information of stop-criterion

         if (itcg.eq.1) dif1 = dif

         if (ic%flow.ge.5 .or. itcg.ge.maxcg-5) then
            call eldiv_count(igs, nadh, nslip, nplast, nexter)
            if (npot.lt.100000) then
               write(bufout, 5301) itcg, nadh, nslip, dif, difid
            else
               write(bufout, 5302) itcg, nadh, nslip, dif, difid
            endif
            call write_log(1, bufout)
 5301       format (6x,i6, ', size A, S =',2i6, ', |dPt(k)|, Eps|Pt|:', 2g11.3)
 5302       format (6x,i6, ', size A, S =',2i8, ', |dPt(k)|, Eps|Pt|:', 2g11.3)
         endif

         ! 12. Prepare for next iteration: set residual r

         if ((lchanged .or. dif.gt.difid) .and. itcg.lt.maxcg) then
            if (it_inn.le.0) then

               ! In outer iterations: compute full t, n, r

               ! 12.1 compute argument \theta, primarily for elements in slip
               !      also: normal (n) and tangential (t) direction vectors

               do ii = 1, npot
                  n%vn(ii) = atan2(ps%vy(ii), ps%vx(ii))
                  n%vx(ii) = cos(n%vn(ii))
                  n%vy(ii) = sin(n%vn(ii))
                  t%vx(ii) = -n%vy(ii)
                  t%vy(ii) =  n%vx(ii)
               enddo

               if (ic_nmdbg.ge.15) then
                  ii = iidbg
                  write(bufout,'(a,i4,5(a,g14.6))') 'ii=',ii,                           &
                     ': nx=',n%vx(ii),', ny=',n%vy(ii),                                 &
                     ', tx=',t%vx(ii),', ty=',t%vy(ii),                                 &
                     ', arg=',n%vn(ii)
                  call write_log(1, bufout)
               endif

               ! 12.2 Compute residual r
               !     - adhesion: r = -ss
               !     - slip:     r = t * t' * -ss

               call gf3_copy(AllElm, ss, r, ikTANG)
               call gf3_scal(AllElm, -1.0d0, r, ikTANG)

               do ii = 1, npot
                  if (igs%el(ii).eq.Slip) then
                     perp = t%vx(ii)*r%vx(ii) + t%vy(ii)*r%vy(ii)
                     r%vx(ii) = perp * t%vx(ii)
                     r%vy(ii) = perp * t%vy(ii)
                  endif
               enddo

               if (ic_nmdbg.ge.15) then
                  ii = iidbg
                  write(bufout,'(a,i4,4(a,g14.6))') 'ii=',ii,                           &
                     ': sx=',ss%vx(ii),', sy=',ss%vy(ii), ', rx=',r%vx(ii),', ry=',r%vy(ii)
                  call write_log(1, bufout)
               endif

            else

               ! 12.3 In inner iterations: update r using q
               !        r(it) = r(it-1) - alpha * q

               call gf3_axpy(AllInt, -alpha, q, r, ikTANG)
            endif
         endif

      ! end-while (not converged)

      enddo

      err = dif
      conv = 1d0
      if (dif*dif1.gt.0d0 .and. itcg.gt.1) then
         conv = exp (dlog(dif/dif1) / (itcg-1))
      endif
      ! if (conv.lt.1d0) err = dif * (conv / (1-conv))

      if (dif.gt.difid .and. conv.gt.1d0 .and. itcg.ge.maxcg) then
         write(bufout, 6000) maxcg, dif, conv
         call write_log(1, bufout)
 6000    format (' ERROR. MaxCG= ',i3, ' reached in TangCG.', /,                        &
         ' |Pk - Pk-1| = ', g12.4, '  Estim. Avg Rate of Conv =', g12.4)
         ! call abort_run()
      endif

      if (ic%flow.ge.5) then
         write(bufout, 6200) itcg, dif, conv
         call write_log(1, bufout)
 6200    format (14x, 'TangCG:', i4, ' iterations. |Pk - Pk-1|,', ' C_est:', g12.4, f8.4)
      endif

      call gf3_destroy(g)
      call gf3_destroy(n)
      call gf3_destroy(t)
      call gf3_destroy(r)
      call gf3_destroy(z)
      call gf3_destroy(v)
      call gf3_destroy(q)
      call gf3_destroy(pold)
      if (is_sens) then
         call gf3_destroy(ps0)
         call gf3_destroy(tmp)
      endif
      end associate

      end subroutine tangcg

!------------------------------------------------------------------------------------------------------------

      subroutine cnvxgs(ic, cgrid, mater, kin, ws, infl, ledg, outpt1, k, iel, eps,                     &
                        maxgs, omegah, omegas, is_sens, info, itgs, err, conv)
!--purpose: Solves the equations in TANG with a non-linear block variant of the Gauss-Seidel, that takes
!           into account the non-linear constraints. This subroutine uses the function AijPj instead of
!           using a matrix for evaluation of the shift per element.
      implicit none
!--subroutine arguments:
      type(t_ic)                :: ic
      type(t_grid)              :: cgrid
      type(t_material)          :: mater
      type(t_kincns)            :: kin
      type(t_influe)            :: infl
      type(t_leadedge)          :: ledg
      type(t_output)            :: outpt1
      type(t_gridfnc3)          :: ws
      integer,      intent(in)  :: iel(cgrid%ntot), k, maxgs
      logical,      intent(in)  :: is_sens
      integer,      intent(out) :: info, itgs
      real(kind=8), intent(in)  :: omegah, omegas, eps
      real(kind=8), intent(out) :: err, conv
!--local variables :
      integer, parameter :: itsedg = 1
      integer            :: i, ista, iinc, iend, ii, ix, iy, ixinc, iym, mx, my, npot, nadh, nslip,     &
                            nplst, nexter, elprv, idebug, difmax_i, nadh_th, nslip_th, nplst_th,        &
                            difmax_i_th
      real(kind=8)       :: coefs(2,2), coefsv(2,2), tau_c0, k_tau, dif1, difid, pr(3), s(2), dupl(2),  &
                            tauc, tauv, facnel, flx, dif, difmax, dif_th, difmax_th
      logical            :: is_ssrol, use_plast, zledge
      type(t_gridfnc3)   :: ps0, tmp
      character(len=1), parameter :: aset(0:3) = (/ 'E', 'H', 'S', 'P' /)

      associate( cs    => infl%cs,      csv  => infl%csv,    igs   => outpt1%igs,                       &
                 mus   => outpt1%mus,   ps   => outpt1%ps,   ss    => outpt1%ss,                        &
                 upls  => outpt1%upls,  uplv => outpt1%uplv, taucs => outpt1%taucs,                     &
                 taucv => outpt1%taucv, ubnd => ledg%ubnd,   ii2j  => ledg%ii2j%val)

      mx     = cgrid%nx
      my     = cgrid%ny
      npot   = cgrid%ntot
      k_tau  = mater%k_tau
      tau_c0 = mater%tau_c0

      is_ssrol  = ic%tang.eq.3
      use_plast = ic%mater.eq.4 .and. mater%tau_c0.gt.1d-10 .and. mater%tau_c0.le.1d10
      if (.not.use_plast) tau_c0 = 1d20

      ! save initial estimate for modified stopping criterion

      if (is_sens) then
         call gf3_copy_struc(ps, ps0, 'cnvxgs:ps0', .true.)
         call gf3_copy_struc(ps, tmp, 'cnvxgs:tmp', .true.)
         call gf3_copy(AllElm, ps, ps0, ikTANG)
      endif

      ! Set start/increment/end of Gauss-Seidel loop on the basis of the rolling direction:
      !     iterating from trailing to leading edge

      if (abs(kin%chi-pi).lt.0.01d0) then
         ista  =  k     ! chi = pi, rolling to the left
         iinc  = -1
         iend  =  1
         ixinc =  1     ! pointer to next point in x-direction, used for uplv(ii) = upls(ii-ixinc)
      else
         ista  =  1     ! rolling to the right or diagonally
         iinc  =  1
         iend  =  k
         ixinc = -1     ! pointer to next point in x-direction
         if (is_ssrol .and. use_plast .and. abs(kin%chi).gt.0.01d0) then
            write(bufout,'(2a,f6.1,a)') ' WARNING: ConvexGS may be inaccurate for steady rolling ',     &
                'with plasticity at chi=', kin%chi*180d0/pi,' deg.'
            call write_log(1, bufout)
         endif
      endif

      if (ic%flow.ge.5) write(bufout,100) omegah, omegas
      if (ic%flow.ge.5) call write_log(1, bufout)
 100  format(10x,'CnvxGS: Using omegah=',f7.3,', omegas=',f7.3)

      if (cs%use_3bl) then
         flx = cs%flx_3bl
         if (ic%flow.ge.5) write(bufout,110) flx
         if (ic%flow.ge.5) call write_log(1, bufout)
 110     format(10x,'CnvxGS: Using third body flexibility=',g10.3)
      else
         flx = 0d0
      endif

      if (cs%itypcf.eq.0) then
         coefs(1,1) = cs%ga_inv * cs%cf(0,0,1,1) + flx
         coefs(1,2) = cs%ga_inv * cs%cf(0,0,1,2)
         coefs(2,1) = cs%ga_inv * cs%cf(0,0,2,1)
         coefs(2,2) = cs%ga_inv * cs%cf(0,0,2,2) + flx

         coefsv(1,1) = csv%ga_inv * csv%cf(0,0,1,1) + flx
         coefsv(1,2) = csv%ga_inv * csv%cf(0,0,1,2)
         coefsv(2,1) = csv%ga_inv * csv%cf(0,0,2,1)
         coefsv(2,2) = csv%ga_inv * csv%cf(0,0,2,2) + flx
      else
         iym = (my+1) / 2
         coefs(1,1) = cs%ga_inv * cs%cy(0,0,iym,1,1) + flx
         coefs(1,2) = cs%ga_inv * cs%cy(0,0,iym,1,2)
         coefs(2,1) = cs%ga_inv * cs%cy(0,0,iym,2,1)
         coefs(2,2) = cs%ga_inv * cs%cy(0,0,iym,2,2) + flx

         coefsv(1,1) = csv%ga_inv * csv%cy(0,0,iym,1,1) + flx
         coefsv(1,2) = csv%ga_inv * csv%cy(0,0,iym,1,2)
         coefsv(2,1) = csv%ga_inv * csv%cy(0,0,iym,2,1)
         coefsv(2,2) = csv%ga_inv * csv%cy(0,0,iym,2,2) + flx
      end if
      ! write(bufout,*) 'C=',coefsv(1,1),coefsv(1,2),coefsv(2,1), coefsv(2,2)
      ! call write_log(1, bufout)
      ! write(bufout,*) 'W=',ws%vx(1), ws%vy(1)
      ! call write_log(1, bufout)

      ! correction factor for rms-norm: from npot to actual #contact

      call eldiv_count(igs, nadh, nslip, nplst, nexter)
      facnel = sqrt(real(npot)/real(nadh+nslip+nplst))

      itgs  = 0
      dif   = 2d0
      difid = 1d0

      ! Main loop: While (not converged) do

      do while (dif.ge.difid .and. itgs.lt.maxgs)

         itgs   = itgs + 1
         dif    = 0d0
         difmax = 0d0
         nadh   = 0
         nslip  = 0
         nplst  = 0

         ! Estimate the displacement difference at the leading edge of the contact area

         if (is_ssrol .and. mod(itgs,itsedg).eq.0) then
            ! write(*,*) 'Updating UBnd...'
            call subnd(cgrid, ps, cs, ledg)
         endif

         ! Loop over all elements in the contact area

         ! Using asynchronous parallelisation, i.e. the new value of element i may or may not be
         !    incorporated at element j>i.

!_INT_OMP parallel if (omp_in_parallel().eq.0 .and. mx.ge.4 .and. my.ge.4)                              &
!_INT_OMP          default(none)                                                                        &
!_INT_OMP          shared(itgs, ic, k, ista, iend, iinc, ixinc, mx, my, cgrid, iel, igs, ps, mus,       &
!_INT_OMP                 ws, ss, upls, uplv, taucs, taucv, cs, csv, coefs, coefsv, k_tau, tau_c0,      &
!_INT_OMP                 omegah, omegas, eps, ii2j, ubnd, dif, difmax, difmax_i, nadh, nslip, nplst,   &
!_INT_OMP                 aset, is_ssrol)                                                               &
!_INT_OMP          private(i, ii, ix, iy, elprv, pr, s, zledge, idebug, nadh_th, nslip_th, nplst_th,    &
!_INT_OMP                  dif_th, difmax_th, difmax_i_th)

         dif_th    = 0d0
         difmax_th = 0d0
         nadh_th   = 0
         nslip_th  = 0
         nplst_th  = 0

         ! Use static schedule in order to get large chunks, with probably less effect of the asynchronous
         !    parallelisation than when assigning consecutive elements to different threads.

!_INT_OMP do schedule(static)
         do 20 i = ista, iend, iinc

            ! Get the element number in the potential contact area

            ii = iel(i)
            ix = cgrid%ix(ii)
            iy = cgrid%iy(ii)

            ! After 1000 iterations, every other iteration the elements in the Slip-area are skipped

            if (itgs.le.1000 .or. mod(itgs,2).eq.0 .or. igs%el(ii).eq.Adhes) then

               pr(1)   = ps%vx(ii)
               pr(2)   = ps%vy(ii)
               pr(3)   = ps%vn(ii)
               tauc    = taucs%vt(ii)
               tauv    = taucv%vt(ii)
               dupl(1) = upls%vx(ii) - uplv%vx(ii)
               dupl(2) = upls%vy(ii) - uplv%vy(ii)

               ! Determine whether the leading edge equation must be applied

               zledge = (is_ssrol .and. ii2j(ii).gt.0)

               ! Compute the shift: rigid shift plus displacement difference,
               !    due to the current iterand ps (including element ii itself)

               if (.not.is_ssrol) then

                  ! T=1: no leading edge correction, T=2: included in ws

                  s(1) = ws%vx(ii) + AijPj (ii, ikXDIR, ps, jkTANG, cs)
                  s(2) = ws%vy(ii) + AijPj (ii, ikYDIR, ps, jkTANG, cs)
               elseif (.not.zledge) then

                  ! T=3, farther away from leading edge:
                  !      using equation (1a), contributions of tang. tractions to both u and u' to be
                  !      computed (csv)

                  s(1) = ws%vx(ii) + AijPj (ii, ikXDIR, ps, jkTANG, csv)
                  s(2) = ws%vy(ii) + AijPj (ii, ikYDIR, ps, jkTANG, csv)
               else

                  ! T=3, near the leading edge:
                  !      using equation (1b), u' replaced by ubnd computed above, contribution to u to be
                  !      computed (cs)

                  s(1) = ws%vx(ii) + AijPj (ii, ikXDIR, ps, jkTANG, cs) - ubnd(ii2j(ii),1)
                  s(2) = ws%vy(ii) + AijPj (ii, ikYDIR, ps, jkTANG, cs) - ubnd(ii2j(ii),2)
               endif

               idebug = 0
               ! if (ii.eq.15) idebug = 5
               elprv = igs%el(ii)

               if (.not.is_ssrol .or. zledge) then
                  call plstrc(ii, ix, iy, igs%el(ii), coefs,  eps, omegah, omegas, pr, mus%vt(ii),      &
                                tau_c0, k_tau, tauv, tauc, s, dupl, idebug)
               else
                  call plstrc(ii, ix, iy, igs%el(ii), coefsv, eps, omegah, omegas, pr, mus%vt(ii),      &
                                tau_c0, k_tau, tauv, tauc, s, dupl, idebug)
               endif

               if (igs%el(ii).ne.elprv .and. (ic%x_nmdbg.ge.4 .or. ic%flow.ge.7)) then
                  write(bufout,'(a,i4,2(a,i3),4a)') ' it=',itgs,': moving element (',cgrid%ix(ii),      &
                        ',',cgrid%iy(ii),') from ', aset(elprv),' to ', aset(igs%el(ii))
                  call write_log(1,bufout)
               endif

               ! Add update to residual (2-norm), store new tractions

               dif_th = dif_th + (ps%vx(ii)-pr(1))**2 + (ps%vy(ii)-pr(2))**2

               if (abs(ps%vx(ii)-pr(1)).gt.difmax_th .or.                                       &
                   abs(ps%vy(ii)-pr(2)).gt.difmax_th) then
                  difmax_th = max(difmax_th, abs(ps%vx(ii)-pr(1)))
                  difmax_th = max(difmax_th, abs(ps%vy(ii)-pr(2)))
                  difmax_i_th = i
               endif
               ps%vx(ii) = pr(1)
               ps%vy(ii) = pr(2)
               ss%vx(ii) = s(1)
               ss%vy(ii) = s(2)
               taucs%vt(ii) = tauc
               upls%vx(ii)  = uplv%vx(ii) + dupl(1)
               upls%vy(ii)  = uplv%vy(ii) + dupl(2)
               if (is_ssrol .and. ix+ixinc.ge.1 .and. ix+ixinc.le.mx) then
                  taucv%vx(ii+ixinc) = taucs%vx(ii)
                  uplv%vx(ii+ixinc)  = upls%vx(ii)
                  uplv%vy(ii+ixinc)  = upls%vy(ii)
               endif
            endif

            if (igs%el(ii).eq.Adhes) then
               nadh_th  = nadh_th + 1
            elseif (igs%el(ii).eq.Slip) then
               nslip_th = nslip_th + 1
            else
               nplst_th = nplst_th + 1
            endif

         ! end for (all elements in contact area)

 20      continue
!_INT_OMP end do

         ! Add contributions of current thread in parallel run

!_INT_OMP critical (cnvxgs_difmax)
         dif = dif + dif_th
         if (difmax_th.gt.difmax) then
            difmax   = difmax_th
            difmax_i = difmax_i_th
         endif
         nadh  = nadh  + nadh_th
         nslip = nslip + nslip_th
         nplst = nplst + nplst_th
!_INT_OMP end critical (cnvxgs_difmax)
!_INT_OMP end parallel

         ! Check for convergence

         dif   = dsqrt (dif/(2d0*k))
         if (.not.is_sens) then
            difid = eps * max(1d-6, facnel * gf3_rms(AllElm, ps, ikTANG))
         else
            call gf3_copy(AllElm, ps, tmp, ikTANG)
            call gf3_axpy(AllElm, -1d0, ps0, tmp, ikTANG)
            difid = eps * max(1d-6, facnel * gf3_rms(AllElm, tmp, ikTANG))
         endif

         !! Relative error w.r.t. 2-norm of solution (at least 1e-6)

         if (itgs.eq.1) dif1 = dif

         !! dif1 = |P(1) - P(0)|

         if (ic%flow.ge.5) then
            if (npot.lt.100000 .and. use_plast) then
               write(bufout, 5301) itgs, nadh, nslip, nplst, dif, difid
            elseif (use_plast) then
               write(bufout, 5302) itgs, nadh, nslip, nplst, dif, difid
            elseif (npot.lt.100000) then
               write(bufout, 5303) itgs, nadh, nslip, dif, difid
            else
               write(bufout, 5304) itgs, nadh, nslip, dif, difid
            endif
            call write_log(1, bufout)
         endif
 5301    format (6x,i6,', size A, S, P =',3i6, ', |dPt(k)|, Eps|Pt|:', 2g11.3)
 5302    format (6x,i6,', size A, S, P =',3i8, ', |dPt(k)|, Eps|Pt|:', 2g11.3)
 5303    format (6x,i6,', size A, S =',   2i6, ', |dPt(k)|, Eps|Pt|:', 2g11.3)
 5304    format (6x,i6,', size A, S =',   2i8, ', |dPt(k)|, Eps|Pt|:', 2g11.3)
         if (.false.) then
            ii = iel(difmax_i)
            write(*,'(a,g12.4,a,2i5,3a,g12.4)') '   largest dif',difmax,                &
                ' at ix,iy=',cgrid%ix(ii), cgrid%iy(ii), ' in ',aset(igs%el(ii))
         endif
         if (.false. .and. itgs.ge.241 .and. itgs.le.245) then
            do ii = 1, npot
               write(*,5400) ii, igs%el(ii), ps%vn(ii), ps%vx(ii), ps%vy(ii)
            enddo
 5400       format(i7,i3,3f16.9)
         endif

      enddo ! while(not converged)

      ! write(*,*) 'Aborting program'
      ! call abort_run()

      conv = 1d0
      if (dif*dif1.ne.0d0 .and. itgs.gt.1) then
         conv = exp (dLog(dif/dif1) / (itgs-1))
      endif
      err = dif

      !! estimate of error in ps
      ! if (conv.lt.1d0) err = dif * conv / (1-conv)

      if (dif.gt.difid .and. conv.gt.1d0 .and. itgs.ge.maxgs) then
         write(bufout, 3000) maxgs, dif, conv
         call write_log(2, bufout)
 3000    format (' ERROR. MaxGS=',i5, ' reached in CnvxGS.', /,                         &
         ' |Pk - Pk-1| = ', g12.4, '  Estim. Avg rate of conv =',g12.4)
      elseif (ic%flow.ge.5) then
         write(bufout, 4200) itgs, dif, conv
         call write_log(1, bufout)
 4200    format (10x, 'CnvxGS: ',i5, ' iterations. |Pk - Pk-1|, C_est:', g12.4,f8.4)
      endif

      info = 0
      if (itgs.ge.maxgs) info = 1
      if (itgs.ge.maxgs .and. conv.gt.0.997) info = 2
      if (itgs.ge.maxgs .and. conv.gt.1d0) info = 3

      if (is_sens) then
         call gf3_destroy(ps0)
         call gf3_destroy(tmp)
      endif

!     iy = 1
!     write(bufout,*) 'ConvexGS, iy=',iy
!     call write_log(1, bufout)
!     do ix = 5, mx
!        ii = ix + (iy-1)*mx
!        write(bufout,'(a,i3,3(a,f12.6))') 'ix=',ix,': upls=',upls%vx(ii),', uplv=',uplv%vx(ii),        &
!                       ', dupl=',upls%vx(ii)-uplv%vx(ii)
!        call write_log(1, bufout)
!     enddo
      end associate

      end subroutine cnvxgs

!------------------------------------------------------------------------------------------------------------

      subroutine stdygs(ic, cgrid, mater, kin, ws, infl, outpt1, k, eps, maxgs, omegah, omegas, &
                        is_sens, info, itgs, err, conv)
!--purpose: Solves the equations for steady state rolling with a non-linear block variant of the Gauss-
!           Seidel method, that takes into account the non-linear constraints. This subroutine uses the
!           function AijPj instead of using a matrix for evaluation of the shift per element.
      implicit none
!--subroutine arguments:
      type(t_ic)                :: ic
      type(t_grid)              :: cgrid
      type(t_material)          :: mater
      type(t_kincns)            :: kin
      type(t_influe)            :: infl
      type(t_output)            :: outpt1
      type(t_gridfnc3)          :: ws
      integer,      intent(in)  :: k, maxgs
      logical,      intent(in)  :: is_sens
      integer,      intent(out) :: info, itgs
      real(kind=8), intent(in)  :: omegah, omegas, eps
      real(kind=8), intent(out) :: err, conv
!--local variables :
      integer          :: iidbg, ii, ix, ixsta, ixinc, ixend, iy, iym, j, jj, jx, mx, my,                  &
                          elprv, npot, nadh, nslip, nplst, nexter, difmax_ii, idebug
      logical          :: use_plast
      real(kind=8)     :: coef(2,2), k_tau, tau_c0, pr(3), s(2), ptabs, ptbnd, dupl(2), tauc, tauv,     &
                          dif, dif1, difid, facnel, difmax, flx
      type(t_gridfnc3) :: dp, du, ps0, tmp
      character(len=1), parameter :: aset(0:3) = (/ 'E', 'H', 'S', 'P' /)

      associate(cs   => infl%cs,   csv  => infl%csv,  igs  => outpt1%igs,  mus  => outpt1%mus,          &
                ps   => outpt1%ps, ss   => outpt1%ss, upls => outpt1%upls, uplv => outpt1%uplv,         &
                taucs => outpt1%taucs, taucv => outpt1%taucv)

      mx     = cgrid%nx
      my     = cgrid%ny
      npot   = cgrid%ntot
      k_tau  = mater%k_tau
      tau_c0 = mater%tau_c0

      iidbg  = -1
      ! iidbg  = 91 + mx*(54-1)

      use_plast = ic%mater.eq.4 .and. mater%tau_c0.gt.1d-10 .and. mater%tau_c0.le.1d10
      if (.not.use_plast) tau_c0 = 1d20

      ! Set loop start/increment/end for processing rows on the basis of rolling direction (0 or 180 deg)

      if (abs(kin%chi).lt.0.01d0) then
         ixsta = 1      ! ixsta == trailing edge
         ixinc = 1      ! ixinc == rolling direction
         ixend = mx     ! ixend == leading edge
      elseif (abs(kin%chi-pi).lt.0.01d0) then
         ixsta = mx
         ixinc = -1
         ixend = 1
      else
         write(bufout,*) 'ERROR: in SteadyGS, the rolling direction should be 0 or 180deg'
         call write_log(1, bufout)
         call abort_run()
      endif

      ! initialise taucs := tau_c0, needed at leading edge boundaries

      call gf3_set(AllElm, tau_c0, taucs, ikXDIR)

      ! initialise grid-functions for iteration variables dp, du

      call gf3_new(dp, 'stdygs:dp', cgrid, igs, .true.)
      call gf3_new(du, 'stdygs:du', cgrid, igs, .true.)

      ! save initial estimate for modified stopping criterion

      if (is_sens) then
         call gf3_copy_struc(ps, ps0, 'stdygs:ps0', .true.)
         call gf3_copy_struc(ps, tmp, 'stdygs:tmp', .true.)
         call gf3_copy(AllElm, ps, ps0, ikTANG)
      endif

      do iy = 1, my
         do ix = 1, mx
            ii = ix + (iy-1)*mx
            if (ix.ne.ixend) then
               dp%vx(ii) = ps%vx(ii) - ps%vx(ii+ixinc)
               dp%vy(ii) = ps%vy(ii) - ps%vy(ii+ixinc)
               du%vx(ii) = upls%vx(ii) - upls%vx(ii+ixinc)
               du%vy(ii) = upls%vy(ii) - upls%vy(ii+ixinc)
            else
               dp%vx(ii) = ps%vx(ii)
               dp%vy(ii) = ps%vy(ii)
               du%vx(ii) = upls%vx(ii)
               du%vy(ii) = upls%vy(ii)
            endif
         enddo
      enddo

      if (ic%flow.ge.5) write(bufout,100) omegah, omegas
      if (ic%flow.ge.5) call write_log(1, bufout)
 100  format(10x,'StdyGS: Using Omegah=',f7.3,', Omegas=',f7.3)

      if (cs%use_3bl) then
         flx = cs%flx_3bl
         if (ic%flow.ge.5) write(bufout,110) flx
         if (ic%flow.ge.5) call write_log(1, bufout)
 110     format(10x,'StdyGS: Using third body flexibility=',g10.3)
      else
         flx = 0d0
      endif

      if (ic%x_nmdbg.ge.1 .and. iidbg.ge.1 .and. iidbg.le.npot) then
         ii = iidbg
         write(bufout,'(a,i6,4(a,g12.4))') ' stdygs: pn(',ii,')=',ps%vn(ii),', mu=', mus%vt(ii),     &
                ', wx=', ws%vx(ii),', wy=',ws%vy(ii)
         call write_log(1, bufout)
      endif

      ! correction factor for rms-norm: from npot to actual #contact

      call eldiv_count(igs, nadh, nslip, nplst, nexter)
      facnel = sqrt(real(npot)/real(nadh+nslip+nplst))

      itgs  = 0
      dif   = 2d0
      difid = 1d0

      ! Main loop: While (not converged) do

      do while (dif.ge.difid .and. itgs.lt.maxgs)

         itgs   = itgs + 1
         dif    = 0d0
         difmax = 0d0
         nadh   = 0
         nslip  = 0
         nplst  = 0

         ! Loop over all rows of the contact area
         ! Process rows in forward order (alternative: forward in odd, reversed in even passes)

         do 20 iy = 1, my

            ! Loop over all points ix within row iy, from trailing (ixsta) to leading edge (ixend)

            ix     = ixsta - ixinc
            do while (ix.ne.ixend)

               ix = ix + ixinc
               ii = ix + (iy-1)*mx

               if (igs%el(ii).ge.Adhes) then

                  ! locate previous point in Exterior or Slip area -- where dpx lands

                  jx = ix - ixinc
                  jj = jx + (iy-1)*mx
                  do while(igs%el(jj).eq.Adhes .and. jx.ne.ixsta)
                     jx = jx - ixinc
                     jj = jx + (iy-1)*mx
                  enddo

                  ! determine coefficient matrix for update -1 at jx & +1 at ix

                  if (cs%itypcf.eq.0) then
                     coef(1,1)  = cs%ga_inv * (cs%cf(0,0,1,1) - cs%cf(jx-ix,0,1,1)) + flx
                     coef(1,2)  = cs%ga_inv * (cs%cf(0,0,1,2) - cs%cf(jx-ix,0,1,2))
                     coef(2,1)  = cs%ga_inv * (cs%cf(0,0,2,1) - cs%cf(jx-ix,0,2,1))
                     coef(2,2)  = cs%ga_inv * (cs%cf(0,0,2,2) - cs%cf(jx-ix,0,2,2)) + flx
                  else
                     iym = (my+1) / 2
                     coef(1,1)  = cs%ga_inv * (cs%cy(0,0,iym,1,1) - cs%cy(jx-ix,0,iym,1,1)) + flx
                     coef(1,2)  = cs%ga_inv * (cs%cy(0,0,iym,1,2) - cs%cy(jx-ix,0,iym,1,2))
                     coef(2,1)  = cs%ga_inv * (cs%cy(0,0,iym,2,1) - cs%cy(jx-ix,0,iym,2,1))
                     coef(2,2)  = cs%ga_inv * (cs%cy(0,0,iym,2,2) - cs%cy(jx-ix,0,iym,2,2)) + flx
                  end if

                  ! Compute shift: rigid shift + displacement difference, due to the current iterand dp,
                  ! including element ii.

                  s(1) = ws%vx(ii) + AijPj (ii, ikXDIR, dp, jkTANG, cs)
                  s(2) = ws%vy(ii) + AijPj (ii, ikYDIR, dp, jkTANG, cs)

                  ! Set traction bound, copy current tractions for evaluating iteration update later on

                  pr(1) = ps%vx(ii)
                  pr(2) = ps%vy(ii)
                  pr(3) = ps%vn(ii)
                  tauc  = taucs%vt(ii)                  ! TODO: yp uses rolling in x-dir with dq=1
                  tauv  = taucs%vt(ii+ixinc)            ! TODO: ii+ixinc can be larger than the potential
                  dupl(1) = upls%vx(ii) - upls%vx(ii+ixinc)   !       contact
                  dupl(2) = upls%vy(ii) - upls%vy(ii+ixinc)
                  elprv = igs%el(ii)

                  ! Solve the equations for the current element

                  if (ic%x_nmdbg.ge.7 .and. ii.eq.iidbg) then
                     write(bufout,'(a,i4,a,7g11.3)') ' itgs',itgs,': taucs=', (taucs%vt(ii+j), j=-3,2)
                     call write_log(1,bufout)
                  endif

                  idebug = 0
                  if (ic%x_nmdbg.ge.2 .and. ii.eq.iidbg .and. itgs.ge.1 .and. itgs.le.-5) idebug = 5
                  call plstrc(ii, ix, iy, igs%el(ii), coef, eps, omegah, omegas, pr, mus%vt(ii),       &
                                 tau_c0, k_tau, tauv, tauc, s, dupl, idebug)

                  if (igs%el(ii).ne.elprv .and. (ic%x_nmdbg.ge.4 .or. ic%flow.ge.7)) then
                     write(bufout,'(a,i4,2(a,i3),4a)') ' it=',itgs,': moving element (',cgrid%ix(ii),   &
                        ',',cgrid%iy(ii),') from ', aset(elprv),' to ', aset(igs%el(ii))
                     call write_log(1,bufout)
                  endif
                  if (ic%x_nmdbg.ge.1 .and. ii.eq.iidbg) then
                     write(bufout,'(a,i4,2(a,i3),2a,2(a,f9.3),2(a,g12.4))') ' itgs',itgs,': el (',      &
                         ix,',', iy,') in ',aset(igs%el(ii)),', pr=',pr(1),',',pr(2),', s=',s(1),',',s(2)
                     call write_log(1,bufout)
                  endif
                  if (max(abs(pr(1)),abs(pr(2))).gt.1d10) then
                     write(bufout,'(a,i4,2(a,i3),2a,2(a,g12.3),2(a,g12.4))') ' itgs',itgs,': el (',      &
                         ix,',', iy,') in ',aset(igs%el(ii)),', pr=',pr(1),',',pr(2),', s=',s(1),',',s(2)
                     call write_log(1,bufout)
                     call abort_run()
                  endif

                  if (igs%el(ii).eq.Adhes) then
                     nadh  = nadh + 1
                  elseif (igs%el(ii).eq.Slip) then
                     nslip = nslip + 1
                  else
                     nplst = nplst + 1
                  endif

                  ! Add update to residual (2-norm), update location of max update

                  dif = dif + (pr(1)-ps%vx(ii))**2 + (pr(2)-ps%vy(ii))**2

                  if (abs(pr(1)-ps%vx(ii)).gt.difmax .or. abs(pr(2)-ps%vy(ii)).gt.difmax) then
                     difmax = max(difmax, abs(pr(1)-ps%vx(ii)))
                     difmax = max(difmax, abs(pr(2)-ps%vy(ii)))
                     difmax_ii = ii
                  endif

                  ! Compute new traction increments, store tractions and slip

                  dp%vx(ii) = dp%vx(ii) + pr(1) - ps%vx(ii)
                  dp%vy(ii) = dp%vy(ii) + pr(2) - ps%vy(ii)
                  du%vx(ii) = dupl(1)
                  du%vy(ii) = dupl(2)
                  ps%vx(ii) = pr(1)
                  ps%vy(ii) = pr(2)
                  ss%vx(ii) = s(1)
                  ss%vy(ii) = s(2)
                  taucs%vx(ii) = tauc
                  if (.false. .and. ii.eq.iidbg) then
                     write(bufout,*) 'new taucs=',taucs%vx(ii)
                     call write_log(1, bufout)
                  endif

                  ! Propagate change in dp to all elements with jx<ix
                  !  - jx in adhesion: keep dp, update ps
                  !  - jx in slip:     keep ps, update dp
                  !  - jx in plastic:  keep ps, update dp
                  !  - jx in exterior: update dp at trailing edge
                  !  - everywhere:     update upls, update taucs
                  ! Note: using loop over all points to the left, may be done smarter.

                  do jx = ix-ixinc, ixsta, -ixinc
                     jj = jx + (iy-1)*mx
                     if (igs%el(jj).eq.Adhes) then
                        ps%vx(jj) = ps%vx(jj+ixinc) + dp%vx(jj)
                        ps%vy(jj) = ps%vy(jj+ixinc) + dp%vy(jj)
                        ptabs = sqrt(ps%vx(jj)**2 + ps%vy(jj)**2)
                        ptbnd = min(mus%vt(jj)*ps%vn(jj), taucs%vt(jj)) 
                        if (ptabs.gt.ptbnd) then
                           ps%vx(jj) = ps%vx(jj) * ptbnd / ptabs
                           ps%vy(jj) = ps%vy(jj) * ptbnd / ptabs
                           dp%vx(jj) = ps%vx(jj) - ps%vx(jj+ixinc)
                           dp%vy(jj) = ps%vy(jj) - ps%vy(jj+ixinc)
                        endif
                     elseif (igs%el(jj).eq.Slip .or. igs%el(jj).eq.Plast) then
                        dp%vx(jj) = ps%vx(jj) - ps%vx(jj+ixinc)
                        dp%vy(jj) = ps%vy(jj) - ps%vy(jj+ixinc)
                     elseif (igs%el(jj).le.Exter .and. igs%el(jj+ixinc).ge.Adhes) then
                        dp%vx(jj) = -ps%vx(jj+ixinc)
                        dp%vy(jj) = -ps%vy(jj+ixinc)
                     endif
                     uplv%vx(jj) = upls%vx(jj+ixinc)
                     uplv%vy(jj) = upls%vy(jj+ixinc)
                     upls%vx(jj) = upls%vx(jj+ixinc) + du%vx(jj)
                     upls%vy(jj) = upls%vy(jj+ixinc) + du%vy(jj)

                     if (igs%el(jj).eq.Adhes .or. igs%el(jj).eq.Exter) then
                        taucs%vt(jj) = taucs%vt(jj+ixinc)
                     elseif (igs%el(jj).eq.Slip) then
                        ptabs = sqrt(ps%vx(jj)**2 + ps%vy(jj)**2)
                        if (taucs%vt(jj).gt.ptabs) taucs%vt(jj) = max(ptabs, taucs%vt(jj+ixinc))
                     elseif (igs%el(jj).eq.Plast) then
                        if (k_tau.ge.0d0) then
                           taucs%vt(jj) = max(taucs%vt(jj+ixinc), taucs%vt(jj))
                        else
                           taucs%vt(jj) = min(taucs%vt(jj+ixinc), taucs%vt(jj))
                        endif
                     endif
                  enddo

         ! end for (all elements in actual contact area)

               endif !          end if (interior)
            enddo    !       end while (ixend)
 20      continue    !    end for (my)

         ! Check for convergence:
         ! Compute relative error w.r.t. 2-norm of solution (at least 1e-6)

         dif   = dsqrt (dif/(2d0*k))
         if (.not.is_sens) then
            difid = eps * max(1d-6, facnel * gf3_rms(AllElm, ps, ikTANG))
         else
            call gf3_copy(AllElm, ps, tmp, ikTANG)
            call gf3_axpy(AllElm, -1d0, ps0, tmp, ikTANG)
            difid = eps * max(1d-6, facnel * gf3_rms(AllElm, tmp, ikTANG))
         endif

         ! dif1 = |P(1) - P(0)|, used for rate of convergence

         if (itgs.eq.1) dif1 = dif

         ! Write progress message & detailed info  when requested

         if (ic%flow.ge.5) then
            if (npot.lt.100000 .and. use_plast) then
               write(bufout, 5301) itgs, nadh, nslip, nplst, dif, difid
            elseif (use_plast) then
               write(bufout, 5302) itgs, nadh, nslip, nplst, dif, difid
            elseif (npot.lt.100000) then
               write(bufout, 5303) itgs, nadh, nslip, dif, difid
            else
               write(bufout, 5304) itgs, nadh, nslip, dif, difid
            endif
            call write_log(1, bufout)
         endif
 5301    format (6x,i6,', size A, S, P =',3i6, ', |dPt(k)|, Eps|Pt|:', 2g11.3)
 5302    format (6x,i6,', size A, S, P =',3i8, ', |dPt(k)|, Eps|Pt|:', 2g11.3)
 5303    format (6x,i6,', size A, S =',   2i6, ', |dPt(k)|, Eps|Pt|:', 2g11.3)
 5304    format (6x,i6,', size A, S =',   2i8, ', |dPt(k)|, Eps|Pt|:', 2g11.3)

         if (ic%x_nmdbg.ge.5 .and. .false.) then
            ii = difmax_ii
            write(bufout,'(a,g12.4,a,2i5,3a,g12.4)') '   largest dif',difmax,                &
                ' at ix,iy=',cgrid%ix(ii), cgrid%iy(ii), ' in ',aset(igs%el(ii))
            call write_log(1, bufout)
         endif
         if (ic%x_nmdbg.ge.5 .and. .false.) then
            ii = idamax(npot, ps%vy, 1)
            write(bufout,*) 'absmax(py)=',ps%vy(ii),' at ii=',ii
            call write_log(1,bufout)
            if (abs(ps%vy(ii)).gt.1d0) call abort_run()
         endif
         if (ic%x_nmdbg.ge.5 .and. itgs.eq.-1) then
            ii = difmax_ii
            write(bufout,5400) ii, igs%el(ii), ps%vn(ii), ps%vx(ii), ps%vy(ii)
            call write_log(1, bufout)
         endif
         if (ic%x_nmdbg.ge.5 .and. itgs.ge.261 .and. itgs.le.299) then
            write(bufout,'(a,i3.3,a)') 'sol', itgs, '=['
            call write_log(1, bufout)
            do ii = 1, npot
               write(bufout,5400) ii, igs%el(ii), dp%vx(ii), ps%vx(ii), dp%vy(ii), ps%vy(ii)
               call write_log(1, bufout)
            enddo
            call write_log('];')
 5400       format(i7,i3,4f16.9)
            if (itgs.eq.-1) call abort_run()
         endif
         if (ic%flow.ge.15) call wrigs (igs, .false., 0d0)

      enddo ! while(not converged)


      if (ic%x_nmdbg.ge.1 .and. iidbg.ge.1 .and. iidbg.le.npot) then
         ii = iidbg
         write(bufout,'(a,i6,4(a,g12.4))') ' stdygs: px(',ii,')=',ps%vx(ii),', py=', ps%vy(ii),     &
                ', sx=', ss%vx(ii),', sy=',ss%vy(ii)
         call write_log(1, bufout)
      endif

      ! Estimate the rate of convergence conv and the error in ps err

      conv = 1d0
      if (dif*dif1.ne.0d0 .and. itgs.gt.1) then
         conv = exp (dlog(dif/dif1) / (itgs-1))
      endif
      err = dif
      ! extrapolation: if (conv.lt.1d0) err = dif * conv / (1.-conv)

      ! Check convergence, write error message and set info-flag

      if (dif.gt.difid .and. conv.gt.1d0 .and. itgs.ge.maxgs) then
         write(bufout, 3000) maxgs, dif, conv
         call write_log(2, bufout)
 3000    format (' ERROR. MaxGS=',i5, ' reached in StdyGS.', /,                         &
         ' |Pk - Pk-1| = ', g12.4, '  Estim. Avg rate of conv =',g12.4)
      elseif (ic%flow.ge.5) then
         write(bufout, 4200) itgs, dif, conv
         call write_log(1, bufout)
 4200    format (10x, 'StdyGS: ',i5, ' iterations. |Pk - Pk-1|, C_est:', g12.4,f8.4)
      endif

      info = 0
      if (itgs.ge.maxgs) info = 1
      if (itgs.ge.maxgs .and. conv.gt.0.997) info = 2
      if (itgs.ge.maxgs .and. conv.gt.1d0) info = 3

      call gf3_destroy(dp)
      call gf3_destroy(du)
      if (is_sens) then
         call gf3_destroy(ps0)
         call gf3_destroy(tmp)
      endif

!     iy = 1
!     write(bufout,*) 'SteadyGS, iy=',iy
!     call write_log(1, bufout)
!     do ix = 1, mx
!        ii = ix + (iy-1)*mx
!        write(bufout,'(a,i3,3(a,f12.6))') 'ix=',ix,': upls=',upls%vx(ii),', uplv=',uplv%vx(ii),        &
!                       ', dupl=',upls%vx(ii)-uplv%vx(ii)
!        call write_log(1, bufout)
!     enddo
      end associate

      end subroutine stdygs

!------------------------------------------------------------------------------------------------------------

      subroutine elmtrc(ii, ix, iy, el, coef, eps, omegah, omegas, pr, mus, s, idebug)
!--Purpose: solve the tractions for one element, using plstrc, for case without plasticity
!--subroutine arguments:
      implicit none
      integer,        intent(in)    :: ii, ix, iy, idebug
      integer,        intent(inout) :: el
      real(kind=8),   intent(in)    :: eps, omegah, omegas, coef(2,2)
      real(kind=8),   intent(inout) :: pr(3), mus, s(2)
!--local variables:
      real(kind=8), parameter :: tau_c0 = 1d20, taucv = 1d20, k_tau = 0d0
      real(kind=8)            :: taucs, dupl(2)

      taucs = taucv
      call plstrc(ii, ix, iy, el, coef, eps, omegah, omegas, pr, mus, tau_c0, k_tau, taucv, taucs,      &
                  s, dupl, idebug)

      end subroutine elmtrc

!------------------------------------------------------------------------------------------------------------

      subroutine plstrc (ii, ix, iy, el, coef, eps, omegah, omegas, pr, mus, tau_c0, k_tau, taucv,      &
                        taucs, s, dupl, idebug)
!--Purpose: solve the tractions for one element ii=(ix,iy), and thereby determine whether it should
!           be in the Adhesion, Slip or Plasticity area.
!
!           The tractions for element I (ii) are given in pr. We define new tractions
!               P_I1 = pr(1) + dp(1)
!               P_I2 = pr(2) + dp(2)
!
!           The shift due to other elements and due to the initial value of pr is given in s.
!           The new shift is then described by
!               S_I1 =  C11 * dp(1) + C12 * dp(2) + s(1) + dupl(1)
!               S_I2 =  C21 * dp(1) + C22 * dp(2) + s(2) + dupl(2)
!
!           The equations for adhesion (el=Adhes) read
!               S_I1    =  0, S_I2    =  0
!               Dupl_I1 =  0, Dupl_I2 =  0, Tau_c = Tau_o
!
!           In the plastic area,
!               S_I1    =  0, S_I2    =  0
!               Tau_c   = Tau_o + ktau * |Dupl_I|
!               P_I1^2 + P_I2^2 = Tau_c^2
!               Dupl_I // P_I, or arg(Dupl_I) = arg(P_I)
!
!           For an element in the slip area (el=Slip) they are
!               Dupl_I1 =  0, Dupl_I2 =  0, Tau_c = Tau_o
!               P_I1^2 + P_I2^2 = (P_In*mus)^2
!               S_I // -P_I, or arg(S_I) = -arg(P_I)
!
!           Here S_I and P_I are the shift and the new traction in element I,
!               S    is the shift according to the rigid shift and the current iterand Ps in all
!                    elements including the current element (input value of pr)
!               Pr   is the old (input)/new (output) iterand for the traction in element I
!               dP   is the update required to the old Pr in order to satisfy the equations.
!               C    is the 2x2 matrix of influence coefficients for element I
!               mus  is the prescribed friction coefficient
!               taucs is the yield limit
!               taucv is the yield limit at the previous time
!               ktau is the rate of work hardening
!           Note: in case of a third body layer, the corresponding flexibility Flx = L = H3/G3 is
!               incorporated in C11 and C22.
!
!           On input the variable el gives the initial (i.e. preferred) state for element I;
!           on output it gives the resulting state.
!
!--method:  The equations for adhesion are linear and can be solved directly. Relaxation may be
!           applied to enlarge (over-relaxation, omegah>1) or reduce (under-relaxation, omegah<1)
!           the step dP. If the new tractions violate the traction bound, the element is moved to
!           the slip area.
!
!           The equations for the slip area are non-linear and are solved using Newton-iteration.
!           See the scientific manual (?) for a derivation of the equations. Again relaxation may be
!           applied, in this case with parameter omegas. Here the relaxation works on the direction
!           of the tractions, i.e. alpha=arg(pr), to keep the new tractions on the traction bound.
!           If the direction of the tractions does not oppose the direction of the shift, the element
!           is moved to the adhesion area.
!
!           The equations for the plastic are are 2D only, and therefore linear.
!
!           If neither of the choices appears appropriate, the initial choice is calculated again.
!--subroutine arguments:
      implicit none
      integer,        intent(in)    :: ii, ix, iy, idebug
      integer,        intent(inout) :: el
      real(kind=8),   intent(in)    :: eps, omegah, omegas, coef(2,2), k_tau, taucv, tau_c0
      real(kind=8),   intent(inout) :: pr(3), mus, s(2), taucs
      real(kind=8),   intent(out)   :: dupl(2)
!--local variables:
      integer      :: itry, violat, itnr
      real(kind=8) :: pold1, pold2, detinv, dp(2), pabs1, pabs2
      real(kind=8) :: si(2), si0(2), sabs, fk(2), gradf(2,2), det, alph0, alph1, dalph
      integer,      parameter :: maxnr = 10
      real(kind=8), parameter :: prc = 0.0001d0
      character(len=1), parameter :: aset(0:3) = (/ 'E', 'H', 'S', 'P' /)

      real(kind=8) :: trc_bound

      real(kind=8) :: theta, wstar(2), dupl_abs, s_err(2), cs, sn

      real(kind=8) :: theta_prv, utot(2), uel(2), cso(2)
      real(kind=8) :: theta_given, theta_target, theta_effective, dfdx

      ! precondition: required relative accuracy, multiplication factor on eps.
      !      Note that eps is for 2-norm over all elements, here we consider a single element only.

      if (idebug.ge.5) then
         write(bufout,'(a,i7,2(a,i3),2a)') ' Element',ii,' (',ix,',',iy,'): input state ',aset(el)
         call write_log(1, bufout)
         write(bufout,'(3(a,g12.4),a)') ' eps    =  ', eps,'; omegah =', omegah,'; omegas =',omegas,';'
         call write_log(1, bufout)
         write(bufout,'(4(a,g12.4),a)') ' tau_c0 =  ',tau_c0, '; k_tau = ',k_tau, '; taucv = ',taucv,      &
                '; taucs =', taucs,';'
         call write_log(1, bufout)
         write(bufout,'(4(a,g12.4),a)') ' coef   = [', coef(1,1),',', coef(1,2),';', coef(2,1),',',       &
                coef(2,2), '];'
         call write_log(1, bufout)
         write(bufout,'(4(a,g12.4),a)') ' pr     = [', pr(1),';', pr(2),';', pr(3), ']; mus =', mus,';'
         call write_log(1, bufout)
         write(bufout,'(4(a,g12.4),a)') ' s      = [',s(1),';', s(2),']; dupl = [', dupl(1),';',dupl(2),'];'
         call write_log(1, bufout)
      endif

      ! Store old iterand for use in relaxation later on

      pold1 = pr(1)
      pold2 = pr(2)

      ! Subtract contribution of pr from shift s, such that si0 becomes the shift due to rigid shift
      ! and to tractions of other elements

      si0(1) = s(1) - coef(1,1)*pr(1) - coef(1,2)*pr(2)
      si0(2) = s(2) - coef(2,1)*pr(1) - coef(2,2)*pr(2)

      ! Solve the equations for the current element
      !   itry = 1: current choice of element division
      !   itry = 2: constraint violated, adapted choice
      !   itry = 3: both constraints violated, adapt choice once more
      !   itry = 4: all three constraints violated, initial choice

      itry   = 0
      violat = 1
      do while (itry.le.3 .and. violat.eq.1)
         itry = itry + 1

         if (idebug.ge.5) then
            write(bufout,'(a,i7,2(a,i3),a,i2,2a)') ' Element',ii,' (',ix,',',iy,'): itry =',itry,       &
                ', state ',aset(el)
            call write_log(1, bufout)
         endif

         if (el.eq.Adhes) then

            ! Solve the equations for an element in the Adhesion area
            ! =================================================================

            ! set plasticity part of the solution: dupl = 0, taucs = unchanged

            dupl(1) = 0.d0
            dupl(2) = 0.d0
            taucs   = taucv

            ! 1) solve 2x2 system of equations for update dp

            trc_bound = mus*pr(3)

            ! AiiInv(1,1) =   coef(2,2) / det
            ! AiiInv(2,1) = - coef(2,1) / det
            ! AiiInv(1,2) = - coef(1,2) / det
            ! AiiInv(2,2) =   coef(1,1) / det

            detinv = 1d0 / (coef(1,1)*coef(2,2) - coef(1,2)*coef(2,1))

            ! si = s + coef * dp = 0 --> dp = AiiInv * -s

            dp(1) = ( - coef(2,2) * s(1) + coef(1,2) * s(2) ) * detinv
            dp(2) = (   coef(2,1) * s(1) - coef(1,1) * s(2) ) * detinv

            if (idebug.ge.5) then
               write(bufout,'(a,2g16.8)') ' elmtrc(adh):             dp=',dp(1), dp(2)
               call write_log(1,bufout)
            endif

            ! 2) apply relaxation and check the traction bound

            !  - compute pr with relaxation, check bound

            pr(1) = pold1 + omegah*dp(1)
            pr(2) = pold2 + omegah*dp(2)
            pabs2 = sqrt ( pr(1)**2 +  pr(2)**2 )

            if (idebug.ge.5) then
               write(bufout,'(a,2g16.8,a,g16.8)') ' elmtrc(adh):      pr(omega)=',pr(1), pr(2),         &
                        ', trcbnd=', min(trc_bound,taucs)
               call write_log(1,bufout)
            endif

            if (pabs2.le.min(trc_bound,taucs)) then

               ! - traction bound ok, accept solution

               violat = 0
            else

               ! - check bound for solution without relaxation

               pabs1 = sqrt ( (pold1+dp(1))**2 + (pold2+dp(2))**2 )

               if (pabs1.le.min(trc_bound,taucs)) then

               ! - reduce relaxation to stay within traction bound

                  pr(1:2) = pr(1:2) * min(trc_bound,taucs) / pabs2
                  violat = 0
               else

               ! - reducing does not help, traction bound violated

                  violat = 1

               endif
            endif

            ! compute shift according to pr

            si(1:2) = si0(1:2) + coef(1:2,1)*pr(1) + coef(1:2,2)*pr(2)

            if (idebug.ge.5) then
               write(bufout,'(a,2g16.8)') ' elmtrc(adh):             si=',si(1), si(2)
               call write_log(1,bufout)
            endif

         elseif (el.eq.Slip) then

            ! Solve the equations for an element in the Slip area
            ! =================================================================

            violat = 0

            trc_bound = mus*pr(3)

            if (idebug.ge.7) then
               write(*,'(2(a,2es13.5))') ' elmtrc(slp): si0=',si0(1),si0(2), ', mu=',mus
               write(*,'(a,4es13.5)') ' coef=',coef(1,1),coef(1,2), coef(2,1),coef(2,2)
            endif

            ! set plasticity part of solution: dupl = 0, taucs = unchanged

            dupl(1:2) = 0.d0
            taucs     = taucv

            ! Set initial estimate for pr, sabs

            sabs  = sqrt (si0(1)**2 + si0(2)**2)
            pr(1:2) = -trc_bound * si0(1:2) / sabs

            ! Compute shift according to pr

            si(1:2) = si0(1:2) + coef(1:2,1)*pr(1) + coef(1:2,2)*pr(2)

            ! Compute initial error in the equations
            !    0 = fk(1) = P(2) * S_I1 - P(1) * S_I2
            !    0 = fk(2) = P(1)^2 + P(2)^2 - (trc_bound)^2

            fk(1) = pr(2)*si(1) - pr(1)*si(2)
            fk(2) = pr(1) **2 + pr(2) **2 - (trc_bound) **2

            if (idebug.ge.5) then
               write(bufout,'(a,2g16.8)') ' elmtrc(slp):         fk-ini=',fk(1), fk(2)
               call write_log(1,bufout)
            endif

            ! Perform Newton iteration until converged,
            ! perform at least one iteration

            itnr = 0
            do while (itnr.eq.0 .or. (itnr.lt.maxnr .and. abs(fk(1))+abs(fk(2)).ge.prc*eps*trc_bound))
               itnr = itnr + 1

               ! Compute gradient of fk w.r.t. P(1) and P(2)

               gradf(1,1) = - 2d0*coef(2,1) * pr(1) - si0(2) + (coef(1,1)-coef(2,2)) * pr(2)
               gradf(1,2) =   2d0*coef(1,2) * pr(2) + si0(1) + (coef(1,1)-coef(2,2)) * pr(1)
               gradf(2,1) =   2d0*pr(1)
               gradf(2,2) =   2d0*pr(2)

               if (idebug.ge.5) then
                  write(bufout,'(a,i4,a,2g16.8,/,29x,2g16.8)') ' elmtrc(slp): it=',itnr,', gradf=',     &
                        gradf(1,1:2), gradf(2,1:2)
                  call write_log(2,bufout)
               endif

               ! Solve the Newton update  dp = -Inv(gradient) * fk

               det = gradf(1,1)*gradf(2,2) - gradf(1,2)*gradf(2,1)

               if (det.eq.0d0) then
                  write(bufout,'(2(a,i4),a,2g12.4,a,/,52x,a,2g12.4,a)')                 &
                     ' ElmTrc: ERROR. Det=0 for element (',ix,',',iy,                   &
                     '), grad=[',gradf(1,1),gradf(1,2),']',                             &
                             '[',gradf(2,1),gradf(2,2),']'
                  call write_log(2, bufout)
                  itnr = maxnr
               else
                  dp(1) = -( gradf(2,2)*fk(1) - gradf(1,2)*fk(2)) / det
                  dp(2) = -(-gradf(2,1)*fk(1) + gradf(1,1)*fk(2)) / det
                  pr(1:2) = pr(1:2) + dp(1:2)
               endif

               if (idebug.ge.5) then
                  write(bufout,'(13x,a,2g16.8)') '      update dp=',dp(1), dp(2)
                  call write_log(1,bufout)
                  write(bufout,'(13x,a,2g16.8)') '   new guess pr=',pr(1), pr(2)
                  call write_log(1,bufout)
               endif

               ! Compute new shift according to pr

               si(1:2) = si0(1:2) + coef(1:2,1)*pr(1) + coef(1:2,2)*pr(2)

               ! Compute new error/residual in the equations

               fk(1) = pr(2) * si(1) - pr(1) * si(2)
               fk(2) = pr(1) **2 + pr(2) **2 - (trc_bound) **2

               if (idebug.ge.5) then
                  write(bufout,'(13x,a,2g16.8)') '    new slip si=',si(1), si(2)
                  call write_log(1,bufout)
                  write(bufout,'(13x,a,2g16.8)') '     new res fk=',fk(1), fk(2)
                  call write_log(1,bufout)
               endif

            enddo ! end Newton iteration to find s and p

            ! write error-message when not converged

            if (abs(fk(1))+abs(fk(2)).ge.prc*eps*trc_bound .and. .false.) then
               write(bufout,3000) maxnr, mus,                                           &
                                    coef(1,1), coef(1,2), pr(1), si0(1),                &
                                    coef(2,1), coef(2,2), pr(2), si0(2)
               call write_log(3, bufout)
 3000          format(' ElmTrc: WARNING. MaxNr=',i3,'  reached. mus=', g12.4/,           &
                  ' A11, A12, ', 2g12.4, '   Pr1, S1: ', 2g12.4, /,                     &
                  ' A21, A22, ', 2g12.4, '   Pr2, S2: ', 2g12.4)
            endif

            ! Apply relaxation

            if (.true.) then

               ! method 1: relaxation on change of direction

               alph0 = atan2(pold2, pold1)
               alph1 = atan2(pr(2), pr(1))
               dalph = alph1 - alph0
               if (dalph.lt.-pi) dalph = dalph + 2d0*pi
               if (dalph.gt. pi) dalph = dalph - 2d0*pi
               alph1 = alph0 + omegas*dalph
               pr(1) = trc_bound * cos(alph1)
               pr(2) = trc_bound * sin(alph1)

               if (idebug.ge.5) then
                  write(bufout,'(13x,a,2g16.8)') '      pr(omega)=',pr(1), pr(2)
                  call write_log(1,bufout)
               endif

            endif

            ! check if there should be some plastic deformation before slip could occur 

            if (idebug.ge.6) then
               write(bufout,'(a,g12.4,a,f10.3,a,f9.1)') ' taucs=',taucs,', trcbnd=',trc_bound,          &
                        ', k_tau=',k_tau
               call write_log(1, bufout)
            endif

            if (taucs.lt.trc_bound) then

               ! if k_tau < 0 and taucs < g, we should never have been here, go to P

               if (k_tau.lt.0d0) violat = 1

               ! the solution for S+P has the tractions pr as computed above, yet part of the
               ! slip si is needed as plasticity dupl to raise taucs to the level of trc_bound

               taucs = trc_bound

               ! taucs = tauc_old + ktau * |dupl|  
               !   --> |dupl| = (trc_bound - tauc_old) / ktau
               !   -->  dupl  = direction(p) * (tauc_old - taucs) / ktau  TODO: explain sign!!

               pabs2 = sqrt (pr(1)**2 + pr(2)**2)

               dupl(1:2) = (taucs-taucv)*pr(1:2) / (pabs2*k_tau)

               si(1:2) = si(1:2) + dupl(1:2)
            endif

            ! Verify that the direction of the slip opposes the direction of tractions;
            ! if not, set violat=1 (false)

            if (abs(pr(1)).gt.abs(pr(2))) then
               if (pr(1)*si(1).gt.0d0) violat = 1
            else
               if (pr(2)*si(2).gt.0d0) violat = 1
            endif

            if (idebug.ge.5) then
               write(bufout,'(a,g12.4,a,f10.3,a,f9.1,a,i2)') ' taucs=',taucs,', trcbnd=',trc_bound,     &
                        ', k_tau=',k_tau,', violat=',violat
               call write_log(1, bufout)
            endif


         else

            ! Solve the equations for an element in the Plastic area
            ! =================================================================

            trc_bound = mus*pr(3)

            ! set the slip part of the solution, s = 0

            si(1:2) = 0d0

            ! something close to a quasi-Newton iteration with (theta, taucs) as independent variables
            ! initial estimate: tau_c = tau_old, theta = opposite to imposed slip

            taucs     = taucv
            theta     = atan2(-si0(2),-si0(1))
            theta_prv = 0d0

            itnr = 0
            do while (itnr.eq.0 .or. (itnr.lt.maxnr .and. abs(theta-theta_prv).ge.eps) )
               itnr = itnr + 1

               ! store theta for convergence checking

               theta_prv = theta

               ! initial traction in direction theta: p*

               cs = cos(theta)
               sn = sin(theta)
               cso(1) = cs
               cso(2) = sn
               pr(1:2) = taucs * cso(1:2)

               wstar(1:2) = si0(1:2) + coef(1:2,1)*pr(1) + coef(1:2,2)*pr(2)

               ! determine plastic deformation that occurs at (taucs, theta)

               dupl_abs = -(wstar(1)*cs + wstar(2)*sn) /                                                &
                 ( ((coef(1,1)*cs + coef(1,2)*sn)*cs + (coef(2,1)*cs + coef(2,2)*sn)*sn) * k_tau + 1.d0 )

               ! If dupl_abs < 0 then we have found an invalid solution
               !    the actual solution should then be dupl_abs = 0, correct for this.

               if (dupl_abs.lt.0.d0) dupl_abs = 0.d0

               dupl(1:2) = dupl_abs * cso(1:2)

               ! calculate new traction (yield limit moves due to plastic deformation)

               pr(1:2) = pr(1:2) + dupl(1:2)*k_tau

               ! determine what a better guess for theta would be
               ! calculate in what directions the current contributions are

               uel(1:2)   = coef(1:2,1)*pr(1) + coef(1:2,2)*pr(2)
               utot(1:2)  = uel(1:2) + dupl(1:2)
               s_err(1:2) = si0(1:2) + utot(1:2)

               theta_given     = theta
               theta_effective = atan2(utot(2),utot(1))
               theta_target    = atan2(-si0(2),-si0(1))

               ! use Newton: f(theta) = (theta_effective - theta_target) = 0

               det = coef(1,1)*coef(2,2) - coef(1,2)*coef(2,1)
               dfdx = det / (  (    coef(2,2)**2    +     coef(1,2)**2   )*sn**2 +                      &
                             2*(coef(2,1)*coef(2,2) + coef(1,1)*coef(1,2))*cs*sn +                      &
                               (    coef(2,1)**2    +     coef(1,1)**2   )*cs**2   )

               theta = theta_given - (theta_effective - theta_target) / dfdx

            enddo ! while itnr

            taucs = taucv + k_tau * dupl_abs
            taucs = max(0.1d0*tau_c0, taucs)

            violat = 0

            ! Check that there is some additional plasticity

            if (dupl_abs.le.tiny) violat = 1

            ! Verify that the direction of plasticity equals the direction of tractions;
            ! if not, set violat=1 (false)

            if (abs(pr(1)).gt.abs(pr(2))) then
               if (pr(1)*dupl(1).lt.0d0) violat = 1
            else
               if (pr(2)*dupl(2).lt.0d0) violat = 1
            endif

            ! Check that taucs doesnt exceed the traction bound mu*pn

            if (taucs.ge.trc_bound) violat = 1

         endif

         ! We are done solving the equations for the elements
         ! =================================================================

         if (itry.le.3 .and. violat.eq.1) then

            ! if a solution is not yet found, change the element division

            if (el.eq.Adhes) then
               if ( k_tau.gt.0d0 .or. trc_bound.le.taucs ) then
                  el = Slip
               else
                  el = Plast
               endif
            elseif (el.eq.Slip) then
               el = Plast
            else
               el = Adhes
            endif
            if (itry.le.2 .and. idebug.ge.4) then
               write(bufout,101) ii, ix, iy, aset(el)
               call write_log(1, bufout)
 101           format(' Moving element',i6,' = (',i3,',',i3,') to ',a1)
            elseif (itry.ge.3 .and. idebug.ge.2) then
               write(bufout,102) ii, ix, iy, aset(el)
               call write_log(1, bufout)
 102           format(' All three constraints violated for el.',i6,'= (',i3,',',i3,'); moving back to ',a)
            endif
         endif

      enddo

      ! Return the new slip in the element

      s(1:2) = si(1:2)

      end subroutine plstrc

!============================================================================================================

end module m_solvpt
