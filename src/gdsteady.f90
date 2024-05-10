!------------------------------------------------------------------------------------------------------------
! Steady rolling contact solver GDsteady - experimental
!
! Copyright 2020-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------

      subroutine gdsteady(ic, npot, mater, maxgd, eps, ws, infl, outpt, solv, itgd, err, lstagn)
!--purpose: Nonlinear Gradient Descent type solver for the equations and constraints for steady
!           rolling with fixed creepages. Matrix-free implementation, using routine VecAijPj for
!           evaluating "A * v", with A represented by cs.
      implicit none
!--subroutine arguments:
      type(t_ic),       intent(in)      :: ic
      integer,          intent(in)      :: npot, maxgd
      type(t_material)                  :: mater
      real(kind=8),     intent(in)      :: eps
      type(t_gridfnc3), intent(in)      :: ws
      type(t_influe)                    :: infl
      type(t_output)                    :: outpt
      type(t_solvers)                   :: solv
      integer,          intent(out)     :: itgd
      real(kind=8),     intent(out)     :: err
      logical,          intent(out)     :: lstagn
!--local variables:
      integer,          parameter :: idebug = 1
      real(kind=8),     parameter :: tiny = 1d-12, fac_v = 1000d0
      character(len=1), parameter :: aset(0:3) = (/ 'E', 'H', 'S', 'P' /)
      logical            :: lchanged, use_plast
      type(t_eldiv)      :: igsold, igsopt
      type(t_gridfnc3)   :: g, dp, dupl, dscl, n, t, r, dv, v, q, D_q, psopt, ssopt, pold, tmp
      integer            :: mx, my, nadh, nslip, nplast, nexter, newadh, newslp, imvp, it_inn, it_fb,   &
                            ii, iym, iidbg, ixdbg, iydbg, its_j, tot_j, ltmp
      real(kind=8)       :: alpha, alpha0, alpha1, beta, coefs(2,2), cosdth, facnel, flx, ga, k_tau,    &
                            tau_c0, snrm, st, ptabs, dif, dif1, difid, difinn, conv,      &
                            resmax, resrms, rii, res_dp, facdif

      associate(cs    => infl%cs,     igs   => outpt%igs,   mus  => outpt%mus,                          &
                ps    => outpt%ps,    ss    => outpt%ss,                                                &
                taucs => outpt%taucs, taucv => outpt%taucv, upls => outpt%upls, uplv => outpt%uplv )

      mx     = ps%grid%nx
      my     = ps%grid%ny
      ga     = cs%ga
      k_tau  = mater%k_tau
      tau_c0 = mater%tau_c0

      use_plast = ic%mater.eq.4 .and. mater%tau_c0.gt.1d-10 .and. mater%tau_c0.le.1d10
      if (.not.use_plast) tau_c0 = 1d20

      if (.false.) then
         ixdbg = 8
         iydbg = 6
         iidbg = ixdbg + (iydbg-1)*mx
      else
         iidbg = -1
      endif

      if (cs%use_3bl) then
         flx = cs%flx_3bl
         if (ic%flow.ge.5) write(bufout,110) flx
         if (ic%flow.ge.5) call write_log(1, bufout)
 110     format(10x,'GDstdy: Using third body flexibility=',g10.3)
      else
         flx = 0d0
      endif

      if (cs%itypcf.eq.0) then
         coefs(1,1) = cs%ga_inv * cs%cf(0,0,1,1) + flx
         coefs(1,2) = cs%ga_inv * cs%cf(0,0,1,2)
         coefs(2,1) = cs%ga_inv * cs%cf(0,0,2,1)
         coefs(2,2) = cs%ga_inv * cs%cf(0,0,2,2) + flx
      else
         iym = (my+1) / 2
         coefs(1,1) = cs%ga_inv * cs%cy(0,0,iym,1,1) + flx
         coefs(1,2) = cs%ga_inv * cs%cy(0,0,iym,1,2)
         coefs(2,1) = cs%ga_inv * cs%cy(0,0,iym,2,1)
         coefs(2,2) = cs%ga_inv * cs%cy(0,0,iym,2,2) + flx
      endif

      ! Initialize grid functions.
      ! Note: they should all be zero in the exterior area at all times.
      !       g is the traction bound, n, t normal and tangential directions,
      !       r is the residual, dscl the diagonal preconditioner
      !       dp are the traction increments, pold for stopping criterion
      !       dv is the search direction, q = A * dv, v = Sig * dv

      call eldiv_new(igsold, ps%grid, nulify=.true.)
      call eldiv_new(igsopt, ps%grid, nulify=.true.)
      call gf3_copy_struc(ps, g    , 'gdstdy:g'    , lzero=.true.)
      call gf3_copy_struc(ps, dp   , 'gdstdy:dp'   , lzero=.true.)
      call gf3_copy_struc(ps, dupl , 'gdstdy:dupl' , lzero=.true.)
      call gf3_copy_struc(ps, dscl , 'gdstdy:dscl' , lzero=.true.)
      call gf3_copy_struc(ps, n    , 'gdstdy:n'    , lzero=.true.)
      call gf3_copy_struc(ps, t    , 'gdstdy:t'    , lzero=.true.)
      call gf3_copy_struc(ps, r    , 'gdstdy:r'    , lzero=.true.)
      call gf3_copy_struc(ps, dv   , 'gdstdy:dv'   , lzero=.true.)
      call gf3_copy_struc(ps, v    , 'gdstdy:v'    , lzero=.true.)
      call gf3_copy_struc(ps, q    , 'gdstdy:q'    , lzero=.true.)
      call gf3_copy_struc(ps, D_q  , 'gdstdy:D_q'  , lzero=.true.)
      call gf3_copy_struc(ps, psopt, 'gdstdy:psopt', lzero=.true.)
      call gf3_copy_struc(ps, ssopt, 'gdstdy:ssopt', lzero=.true.)
      call gf3_copy_struc(ps, pold , 'gdstdy:pold' , lzero=.true.)

      associate(eldiv  => igs%el,  nx  => n%vx,   ny  => n%vy,   theta => n%vn, tx  => t%vx,   &
                trcbnd => g%vt,    px  => ps%vx,  py  => ps%vy,  pn  => ps%vn,  ty  => t%vy,   &
                sx     => ss%vx,   sy  => ss%vy,  rx  => r%vx,   ry  => r%vy,   dvx => dv%vx,  &
                dvy    => dv%vy,   vx  => v%vx,   vy  => v%vy,   qx  => q%vx,   qy  => q%vy,   &
                dpx    => dp%vx,   dpy => dp%vy )

      ! correction factor for rms-norm: pot.con versus actual contact

      call eldiv_count(igs, nadh, nslip, nplast, nexter)
      facnel = sqrt(real(npot)/real(nadh+nslip+nplast))

      ! Compute traction-bound at each element ii - compatibility with previous implementation

      do ii = 1, npot
         trcbnd(ii) = mus%vt(ii) * pn(ii)
      enddo

      ! set default dscl == 1

      call gf3_set(AllElm, 1d0, dscl, ikZDIR)

      ! use zero initial estimate ps == 0

      if (.false.) then
         call gf3_set(AllElm, 0d0, ps, ikXDIR)
         call gf3_set(AllElm, 0d0, ps, ikYDIR)
         do ii = 1, npot
            if (igs%el(ii).eq.Slip) igs%el(ii) = Adhes
         enddo
      endif

      ! This routine uses a relative stop-criterion on the updates dif == ||pk - pk-1||_rms, which are
      !    meaningful to the user.

      itgd   = 0
      it_fb  = -99
      imvp   = 0
      dif    = 2d0
      difid  = 1d0
      it_inn = 0
      difinn = 0d0
      tot_j  = 0
      beta   = 1d0
      lchanged = .false.
      lstagn = .false.

      ! 1. Compute argument \theta, primarily for the elements in slip
      !    also: normal (n) and tangential (t) direction vectors

      do ii = 1, npot
         theta(ii) = atan2(py(ii), px(ii))
         nx(ii)    = cos(theta(ii))
         ny(ii)    = sin(theta(ii))
         tx(ii)    = -ny(ii)
         ty(ii)    =  nx(ii)
      enddo

      call compute_dp(ps, dp)

      ! 1. Compute the slip ss = A * ps + ws.

      call gf3_set(AllElm, 0d0, ss, ikTANG)
      call VecAijPj(igs, AllInt, ss, ikTANG, dp, jkTANG, cs)
      call gf3_axpy(AllInt, 1.0d0, ws, ss, ikTANG)
      imvp = imvp + 1

      ! call gf3_print(ss, 'ss', ikTANG, 4)

      ! 1. Compute diagonal scaling factors for current tractions/slip

      call compute_diagscaling(igs, mus, tau_c0, k_tau, taucv, taucs, coefs, ps, ss, dupl, solv,        &
                dscl, idebug, iidbg)

      ! call gf3_print(dscl, 'dscl', ikZDIR, 4)

      ! 2. Compute scaled residual r
      !     - adhesion: r =          -ss * d_scl
      !     - slip:     r = t * t' * -ss * d_scl

      call gf3_set(AllElm, 0d0, r, ikTANG)
      do ii = 1, npot
         if (eldiv(ii).eq.Adhes) then
            rx(ii) =      -sx(ii) * dscl%vn(ii)
            ry(ii) =      -sy(ii) * dscl%vn(ii)
         elseif (eldiv(ii).eq.Slip) then
            st = tx(ii)*sx(ii) + ty(ii)*sy(ii)
            rx(ii) = -st * tx(ii) * dscl%vn(ii)
            ry(ii) = -st * ty(ii) * dscl%vn(ii)
         endif
      enddo

      if (.false.) then
         ltmp = get_lunit_tmp_use()
         open(ltmp, file='gdstdy0.mat')
         do ii = 1, npot
            write(ltmp,'(i4,i3,7g18.10)') ii, igs%el(ii), r%vx(ii), r%vy(ii), ss%vx(ii), ss%vy(ii), &
               ws%vx(ii), ws%vy(ii), dscl%vn(ii)
         enddo
         close(ltmp)
         call free_lunit_tmp_use(ltmp)
      endif

      ! converged if residual is all zero

      if (gf3_maxabs(AllElm, r, ikTANG).lt.tiny) dif = 0d0

      ! WHILE ( (Changes to el.div OR dif > difid)  AND  itgd < maxgd ) do

      do while ((lchanged .or. dif.gt.difid) .and. itgd.lt.maxgd)

         ! increment iteration counter "it"

         itgd   = itgd + 1
         it_inn = it_inn + 1

         ! 4. Set search direction dv, using steepest descent direction

         if (itgd.eq.-1) then
            write(bufout,'(a,4(i7,5x),2x,4(i7,5x))') 'el=', (igs%el(ii),ii=14,17), (igs%el(ii),ii=14,17)
            call write_log(1, bufout)
            write(bufout,'(a,4f12.6,2x,4f12.6)') 'r= ', (rx(ii),ii=14,17), (ry(ii),ii=14,17)
            call write_log(1, bufout)
         endif

         ! if (itgd.eq.1) call gf3_print(r, 'r', ikTANG, 4)

         call gf3_copy(AllElm, r, dv, ikTANG)

         ! 4. Project onto feasible directions "dv = Emat * r"

         call project_searchdir(igs, g, itgd, it_fb, beta, t, dv, v, solv, fac_v)

         if (itgd.eq.-1) then
            write(bufout,'(a,4f12.6,2x,4f12.6)') 'dv=', (dvx(ii),ii=14,17), (dvy(ii),ii=14,17)
            call write_log(1, bufout)
            write(bufout,'(a,4f12.6,2x,4f12.6)') 'v= ', (vx(ii),ii=14,17), (vy(ii),ii=14,17)
            call write_log(1, bufout)
         endif

         ! 6. compute effect of dv on residual: q = A dv
         !    slip area: project onto tangential direction t

         call VecAijPj(igs, AllInt, q, ikTANG, dv, jkTANG, cs)
         ! write(*,'(a,4f12.6,2x,4f12.6)') 'q= ', (qx(ii),ii=17,20), (qy(ii),ii=17,20)
         imvp = imvp + 1
         if (itgd.eq.-1) then
            call gf3_print(q, 'q', ikTANG, 4)
         endif

         ! 8. estimate order of magnitude alpha0

         if (itgd.le.1) then
            ! initialize alpha0 using fixed step-size (r'r)/(r'*D*q)
            alpha0 = 0d0
            alpha1 = 0d0
            do ii = 1, npot
               alpha0 = alpha0 + (rx(ii)**2 + ry(ii)**2)
               alpha1 = alpha1 + (rx(ii)*qx(ii) + ry(ii)*qy(ii)) * dscl%vn(ii)
            enddo
            ! write(*,*) 'a0=',a0(1), a0(2)
            alpha0 = alpha0 / alpha1
         else
            alpha0 = 0.5d0 * (alpha0 + alpha)
         endif

         ! 8. compute step length alpha
         !   --> minimize |r^k+1|

         if (itgd.eq.-3) then
            open(ltmp, file='gdstdy1.mat')
            do ii = 1, npot
               write(ltmp,'(i4,i3,9g18.10)') ii, igs%el(ii), r%vx(ii), r%vy(ii), dv%vx(ii), dv%vy(ii),  &
                  v%vx(ii), v%vy(ii), q%vx(ii), q%vy(ii), dscl%vn(ii)
            enddo
            close(ltmp)
         endif

         call perform_linesearch(igs, g, dscl, dp, ps, ss, r, dv, v, q, t, itgd, alpha0, alpha, its_j)
         tot_j = tot_j + its_j

         if (itgd.eq.-1) write(*,*) 'alpha=',alpha
         do ii = 1, npot
            d_q%vx(ii) = dscl%vn(ii) * q%vx(ii)
            d_q%vy(ii) = dscl%vn(ii) * q%vy(ii)
         enddo

         beta = alpha * gf3_nrm2(AllInt, d_q, ikTANG) / gf3_nrm2(AllInt, r, ikTANG)

         if (itgd.eq.-1) then
            write(bufout,'(4(a,g14.6))') 'alpha=',alpha,' norm(d_q)=', gf3_nrm2(AllInt, d_q, ikTANG),   &
                ', norm(r)=', gf3_nrm2(AllInt, r, ikTANG), ', beta=', beta
            call write_log(1, bufout)
            do ii = 1, npot
               write(bufout,'(i3,g13.4,g13.4,f10.6)') igs%el(ii), r%vx(ii), d_q%vx(ii), dscl%vn(ii)
               call write_log(1, bufout)
            enddo
         endif

         lchanged = .false.

         ! 9.a Update estimate of tangential tractions

         call eldiv_copy(igs, igsold, ikALL)
         call gf3_copy(AllElm, ps, pold, ikTANG)

         call gf3_axpy(AllInt, alpha, dv, dp, ikTANG)
         if (itgd.eq.-1) write(*,'(a,4f12.6,2x,4f12.6)') 'dp=', (dpx(ii),ii=16,19), (dpy(ii),ii=16,19)

         ! 9.b Clip tractions at tracion bound without changing element division

         call apply_trcbnd(igs, g, dp, ps)
         if (itgd.eq.-1) write(*,'(a,4f12.6,2x,4f12.6)') 'dp=', (dpx(ii),ii=16,19), (dpy(ii),ii=16,19)
         if (itgd.eq.-1) write(*,'(a,4f12.6,2x,4f12.6)') 'p= ', (px(ii),ii=17,20), (py(ii),ii=17,20)

         ! 1. Compute argument \theta, primarily for the elements in slip
         !    also: normal (n) and tangential (t) direction vectors

         do ii = 1, npot
            theta(ii) = atan2(py(ii), px(ii))
            nx(ii)    = cos(theta(ii))
            ny(ii)    = sin(theta(ii))
            tx(ii)    = -ny(ii)
            ty(ii)    =  nx(ii)
         enddo

         ! 1. Compute the slip ss = A * ps + ws.

         call gf3_set(AllElm, 0d0, ss, ikTANG)
         call VecAijPj(igs, AllInt, ss, ikTANG, dp, jkTANG, cs)
         call gf3_axpy(AllInt, 1.0d0, ws, ss, ikTANG)
         imvp = imvp + 1

         if (itgd.eq.-1) then
            call gf3_print(ss, 'ss', ikTANG, 4)
         endif
         if (itgd.eq.-1) then
            write(bufout,'(a,4f12.6,2x,4f12.6)') 's= ', (sx(ii),ii=17,20), (sy(ii),ii=17,20)
            call write_log(1, bufout)
         endif

         ! 9.c compute locally optimized tangential tractions using elmtrc

         call solve_elmtrc(igs, mus, tau_c0, k_tau, taucv, taucs, coefs, ps, ss, dupl, igsopt, psopt, ssopt)

         ! 10. activate constraints where traction bound is reached

         newslp = 0
         do ii = 1, npot
            if (eldiv(ii).eq.Adhes) then
               ptabs = sqrt(px(ii)**2 + py(ii)**2)
               ! if (ptabs.ge.trcbnd(ii)-tiny .and. igsopt%el(ii).eq.Slip) then
               if (ptabs.ge.trcbnd(ii)-tiny) then
                  eldiv(ii) = Slip
                  lchanged = .true.
                  newslp  = newslp + 1
               endif
            endif
         enddo

         if (idebug.ge.3 .and. newslp.gt.0) then
            write(bufout,'(a,i4,a,i4,a)') '   Act', itgd, ': moving', newslp,' elements A -> S'
            call write_log(1, bufout)
         endif

         ! 10. Check direction of slip, move elements to adhesion where needed (release constraints)

         if (.true.) then               ! previously: (.not.lchanged)
            newadh = 0
            do ii = 1, npot
               if (igsold%el(ii).eq.Slip) then
                  snrm   = nx(ii) * sx(ii) + ny(ii) * sy(ii)
                  cosdth = (nx(ii) * psopt%vx(ii) + ny(ii) * psopt%vy(ii)) / trcbnd(ii)
                  if (snrm.gt.0d0 .and.                                                                 &
                      (dscl%vn(ii).ge.0.001d0 .or. igsopt%el(ii).eq.Adhes) .or. cosdth.le.-0.5d0) then
                     eldiv(ii) = Adhes
                     lchanged = .true.
                     newadh  = newadh + 1
                  endif
               endif
            enddo

            if (idebug.ge.3 .and. newadh.gt.0) then
               write(bufout,'(a,i4,a,i4,a)') '   Act', itgd, ': moving', newadh,' elements S -> A'
               call write_log(1, bufout)
            endif
         endif

         ! 1. Compute diagonal scaling factors for current tractions/slip

         call compute_diagscaling(igs, mus, tau_c0, k_tau, taucv, taucs, coefs, ps, ss, dupl, solv,     &
                dscl, idebug, iidbg)
         if (itgd.le.-1) call gf3_print(dscl, 'dscl', ikZDIR, 4)

         ! 2. Compute scaled residual r
         !     - adhesion: r =          -ss * d_scl
         !     - slip:     r = t * t' * -ss * d_scl

         call gf3_set(AllElm, 0d0, r, ikTANG)
         do ii = 1, npot
            if (eldiv(ii).eq.Adhes) then
               rx(ii) =      -sx(ii) * dscl%vn(ii)
               ry(ii) =      -sy(ii) * dscl%vn(ii)
            elseif (eldiv(ii).eq.Slip) then
               st = tx(ii)*sx(ii) + ty(ii)*sy(ii)
               rx(ii) = -st * tx(ii) * dscl%vn(ii)
               ry(ii) = -st * ty(ii) * dscl%vn(ii)
            endif
         enddo

         if (itgd.eq.-3) then
            open(ltmp, file='gdstdy2.mat')
            do ii = 1, npot
               write(ltmp,'(i4,i3,7g18.10)') ii, igs%el(ii), r%vx(ii), r%vy(ii), ss%vx(ii), ss%vy(ii), &
                  dp%vx(ii), dp%vy(ii), dscl%vn(ii)
            enddo
            close(ltmp)
         endif

         ! 11.3 Relative stop-criterion on updates |pk - pk-1|:
         !    dif   = |pk - pk-1|_rms \approx alpha * |Sig * dv|_rms
         !    difid = eps * |pk|_rms

         call gf3_axpy(AllInt, -1.0d0, ps, pold, ikTANG)
         if (beta.ge.0.1d0) then
            facdif = 1d0
         elseif (beta.ge.0.001d0) then  ! enlarge dif to avoid premature "convergence"
            facdif = 0.1d0 / beta
         else
            facdif = 100d0
         endif
         dif   =    facnel * gf3_rms(AllElm, pold, ikTANG) * facdif
         ! dif = alpha * facnel * gf3_rms(AllElm, v, ikTANG)
         difid = eps * max(1d-6, facnel * gf3_rms(AllElm, ps, ikTANG))
         difinn = difinn + dif

         ! 11.5 Print information of stop-criterion

         if (itgd.eq.1) dif1 = dif

         if (ic%flow.ge.5 .or. itgd.ge.maxgd-5) then
            call eldiv_count(igs, nadh, nslip, nplast, nexter)
            if (npot.lt.100000 .and. use_plast) then
               write(bufout, 5301) itgd, nadh, nslip, nplast, dif, difid, beta, its_j
            elseif (use_plast) then
               write(bufout, 5302) itgd, nadh, nslip, nplast, dif, difid, beta, its_j
            elseif (npot.lt.100000) then
               write(bufout, 5303) itgd, nadh, nslip, dif, difid, beta, its_j
            else
               write(bufout, 5304) itgd, nadh, nslip, dif, difid, beta, its_j
            endif
 5301       format (6x,i6,', size A, S, P =',3i6, ', |dPt(k)|, Eps|Pt|:', 2g11.3,', beta=',f11.6,i4)
 5302       format (6x,i6,', size A, S, P =',3i8, ', |dPt(k)|, Eps|Pt|:', 2g11.3,', beta=',f11.6,i4)
 5303       format (6x,i6,', size A, S =',   2i6, ', |dPt(k)|, Eps|Pt|:', 2g11.3,', beta=',f11.6,i4)
 5304       format (6x,i6,', size A, S =',   2i8, ', |dPt(k)|, Eps|Pt|:', 2g11.3,', beta=',f11.6,i4)
            call write_log(1, bufout)
         endif

      ! end-while (not converged)

      enddo

      ! determine rms and maximum value of residual slip

      resmax = -1d0
      do ii = 1, npot
         if (eldiv(ii).eq.Adhes .or. eldiv(ii).eq.Slip) then
            rii = sqrt(rx(ii)**2 + ry(ii)**2)
            if (rii.gt.resmax) then
               resmax = rii
            endif
         endif
      enddo
      resrms = gf3_rms(AllInt, r, ikTANG)

      ! determine order-of-magnitude estimate of update dp needed for slip residual 

      res_dp = resrms / max(coefs(1,1), coefs(2,2))

      ! compute average rate of convergence

      err  = dif
      conv = 1d0
      if (dif*dif1.gt.0d0 .and. itgd.gt.1) then
         conv = exp (dlog(dif/dif1) / (itgd-1))
      endif

      if (ic%flow.ge.5) then
         write(bufout, 6100) itgd, dif, conv
         call write_log(1, bufout)
 6100    format (14x, 'GDstdy:', i4, ' iterations. |Pk - Pk-1|,', ' C_est:', g12.4, f8.4)
         write(bufout,'(19x,i6,a,f5.1,a)') tot_j, ' linesearch iterations, avg =',                      &
                        real(tot_j)/real(itgd),' per iteration'
         call write_log(1, bufout)
      endif

      ! warn in case of divergence or stagnation

      if (lchanged .and. itgd.ge.maxgd) then

         write(bufout, 6200) maxgd, dif
         call write_log(1, bufout)
 6200    format (' ERROR. MaxGD=',i4, ' reached in GDstdy.',                                           &
         ' |Pk - Pk-1| = ', g12.4, '  El.division keeps changing.')
         lstagn = .true.

      elseif (dif.gt.difid .and. conv.gt.1d0 .and. itgd.ge.maxgd) then

         write(bufout, 6250) maxgd, dif, conv
         call write_log(1, bufout)
 6250    format (' ERROR. MaxGD=',i4, ' reached in GDstdy.',                                           &
         ' |Pk - Pk-1| = ', g12.4, '  Estim. Avg Rate of Conv =', g12.4)
         ! stop
         lstagn = .true.

      elseif (conv.lt.1d0-tiny .and. res_dp.gt.5d0*difid/(1d0-conv)) then

         write(bufout, 6300) resrms, res_dp
         call write_log(1, bufout)
 6300    format (' ERROR. Stagnation in GDstdy at residual slip =', g12.4,', traction error =',g12.4)
         itgd = -itgd
         lstagn = .true.

      endif

      call eldiv_destroy(igsold)
      call eldiv_destroy(igsopt)
      call gf3_destroy(g)
      call gf3_destroy(dp)
      call gf3_destroy(dupl)
      call gf3_destroy(dscl)
      call gf3_destroy(n)
      call gf3_destroy(t)
      call gf3_destroy(r)
      call gf3_destroy(dv)
      call gf3_destroy(v)
      call gf3_destroy(q)
      call gf3_destroy(D_q)
      call gf3_destroy(psopt)
      call gf3_destroy(ssopt)
      call gf3_destroy(pold)
      call gf3_destroy(tmp)

      end associate
      end associate
      end subroutine gdsteady

!------------------------------------------------------------------------------------------------------------

   subroutine compute_dp(ps, dp)
!--purpose: compute traction increments dp from tractions ps
      implicit none
!--subroutine arguments:
      type(t_gridfnc3)  :: ps, dp
!--local variables:
      integer      :: ii, ix, iy, mx, my

      associate(px => ps%vx, py => ps%vy, dpx => dp%vx, dpy => dp%vy)

      mx = ps%grid%nx
      my = ps%grid%ny

      do iy = 1, my
         ii = mx + (iy-1) * mx
         dpx(ii) = px(ii)
         dpy(ii) = py(ii)
         do ix = mx-1, 1, -1
            ii = ix + (iy-1) * mx
            dpx(ii) = px(ii) - px(ii+1)
            dpy(ii) = py(ii) - py(ii+1)
         end do
      end do

      end associate
   end subroutine compute_dp

!------------------------------------------------------------------------------------------------------------

   subroutine project_searchdir(igs, g, itgd, it_fb, betaj, t, dv, v, solv, fac_v)
!--purpose: project search direction (dv, v) onto feasible domain, with v // p in slip area
      implicit none
!--subroutine arguments:
      type(t_eldiv),    intent(in)      :: igs
      type(t_gridfnc3)                  :: g, t, dv, v
      type(t_solvers),  intent(in)      :: solv
      integer,          intent(in)      :: itgd
      integer,          intent(inout)   :: it_fb        ! last iteration where fallback was used
      real(kind=8),     intent(in)      :: betaj, fac_v
!--local variables:
      integer      :: imeth, kdown, ii0, ii, ix, iy, k, mx, my
      real(kind=8) :: vt
      real(kind=8), dimension(:), allocatable :: dvx_in, dvy_in

      associate(eldiv  => igs%el, dvx => dv%vx, dvy => dv%vy, tx => t%vx, ty => t%vy,                   &
                trcbnd => g%vt  , vx  => v%vx,  vy  => v%vy )

      mx = dv%grid%nx
      my = dv%grid%ny

      if (betaj.gt.solv%betath .or. itgd-it_fb.le.1) then
         imeth = solv%gd_meth
         kdown = solv%kdown
      else
         imeth = 2
         kdown = solv%kdowfb
         write(bufout,'(a,i4,a,f10.6,a,f6.3,a,i2,a)') ' it=',itgd,': previous beta=',betaj,' <',        &
                solv%betath, ', switching to E_down(',kdown,')'
         call write_log(1, bufout)
         it_fb = itgd
      endif

      if (imeth.eq.1) then

         ! method 1: E_trl, compensate dv_i at trailing edge / slipping elements 

         do iy = 1, my
            ii = mx + (iy-1) * mx
            vx(ii) = 0d0
            vy(ii) = 0d0
            do ix = mx-1, 1, -1
               ii = ix + (iy-1) * mx
               if     (eldiv(ii).eq.Adhes) then   ! use dv as the prescribed value
                  vx(ii)  = vx(ii+1) + dvx(ii)
                  vy(ii)  = vy(ii+1) + dvy(ii)
               elseif (eldiv(ii).eq.Slip) then    ! use v = component of dv orthogonal to tractions p
                  vt      = tx(ii) * dvx(ii) + ty(ii) * dvy(ii)
                  ! using fac_v = O(1/alpha), assure that alpha * |vt| <= O(1) * g
                  vt      = sign(1d0, vt) * min(abs(vt), fac_v*trcbnd(ii))
                  vx(ii)  = tx(ii) * vt
                  vy(ii)  = ty(ii) * vt
                  dvx(ii) = vx(ii) - vx(ii+1)
                  dvy(ii) = vy(ii) - vy(ii+1)
               else ! (eldiv(ii).le.Exter)
                  vx(ii)  = 0d0
                  vy(ii)  = 0d0
                  dvx(ii) = -vx(ii+1)
                  dvy(ii) = -vy(ii+1)
               endif
            end do
         end do

      elseif (imeth.eq.2) then

         ! method 2: E_down(k), compensate dv_i maximally k_down elements downstream
         ! used as fall-back when previous beta_j < beta_thrs

         allocate(dvx_in(mx), dvy_in(mx))
         do iy = 1, my
            ! copy input dvx,dvy(:,iy) for current row to (1:mx)-array
            ii0 = (iy-1) * mx
            dvx_in(1:mx) = dvx(ii0+1:ii0+mx)
            dvy_in(1:mx) = dvy(ii0+1:ii0+mx)

            ! start at leading edge
            ii  = ii0 + mx
            vx(ii) = 0d0
            vy(ii) = 0d0
            do ix = mx-1, 1, -1
               ii = ii0 + ix
               if     (eldiv(ii).eq.Adhes) then   ! use dv as the prescribed value
                  if (ix+kdown.le.mx) then   ! compensate for dv of k points upstream
                     ! write(bufout,*) ' ix=',ix,': ix+k=',ix+kdown,' <= mx=',mx
                     ! call write_log(1, bufout)
                     vx(ii)  = vx(ii+1) + dvx_in(ix) - dvx_in(ix+kdown)
                     vy(ii)  = vy(ii+1) + dvy_in(ix) - dvy_in(ix+kdown)
                  else
                     vx(ii)  = vx(ii+1) + dvx_in(ix)
                     vy(ii)  = vy(ii+1) + dvy_in(ix)
                  endif
                  dvx(ii) = vx(ii) - vx(ii+1)
                  dvy(ii) = vy(ii) - vy(ii+1)
               elseif (eldiv(ii).eq.Slip) then    ! use v = component of dv orthogonal to tractions p
                  vt      = tx(ii) * dvx_in(ix) + ty(ii) * dvy_in(ix)
                  ! using fac_v = O(1/alpha), assure that alpha * |vt| <= O(1) * g
                  vt      = sign(1d0, vt) * min(abs(vt), fac_v*trcbnd(ii))
                  vx(ii)  = tx(ii) * vt
                  vy(ii)  = ty(ii) * vt
                  dvx(ii) = vx(ii) - vx(ii+1)
                  dvy(ii) = vy(ii) - vy(ii+1)
                  do k = 1, kdown            ! all upstream dv are compensated
                     if (ix+k.le.mx) then
                        dvx_in(ix+k) = 0d0
                        dvy_in(ix+k) = 0d0
                     endif
                  enddo
               else ! (eldiv(ii).le.Exter)
                  vx(ii)  = 0d0
                  vy(ii)  = 0d0
                  dvx(ii) = -vx(ii+1)
                  dvy(ii) = -vy(ii+1)
                  do k = 1, kdown            ! all upstream dv are compensated
                     if (ix+k.le.mx) then
                        dvx_in(ix+k) = 0d0
                        dvy_in(ix+k) = 0d0
                     endif
                  enddo
               endif
            end do
         end do
         deallocate(dvx_in, dvy_in)

      elseif (imeth.eq.3) then

         ! method 3: E_keep(f), keep fraction f_decay of upstream dv_i, compensate remainder

         do iy = 1, my
            ! start at leading edge
            ii  = mx + (iy-1)*mx
            vx(ii) = 0d0
            vy(ii) = 0d0
            do ix = mx-1, 1, -1
               ii = ix + (iy-1)*mx
               if     (eldiv(ii).eq.Adhes) then   ! use dv as the prescribed value
                  vx(ii)  = solv%fdecay * vx(ii+1) + dvx(ii)
                  vy(ii)  = solv%fdecay * vy(ii+1) + dvy(ii)
                  dvx(ii) = vx(ii) - vx(ii+1)
                  dvy(ii) = vy(ii) - vy(ii+1)
               elseif (eldiv(ii).eq.Slip) then    ! use v = component of dv orthogonal to tractions p
                  vt      = tx(ii) * dvx(ii) + ty(ii) * dvy(ii)
                  ! using fac_v = O(1/alpha), assure that alpha * |vt| <= O(1) * g
                  vt      = sign(1d0, vt) * min(abs(vt), fac_v*trcbnd(ii))
                  vx(ii)  = tx(ii) * vt
                  vy(ii)  = ty(ii) * vt
                  dvx(ii) = vx(ii) - vx(ii+1)
                  dvy(ii) = vy(ii) - vy(ii+1)
               else ! (eldiv(ii).le.Exter)
                  vx(ii)  = 0d0
                  vy(ii)  = 0d0
                  dvx(ii) = -vx(ii+1)
                  dvy(ii) = -vy(ii+1)
               endif
            end do
         end do

      else

         write(bufout,*) ' INTERNAL ERROR: gd_meth =', imeth,' out of range.'
         call write_log(1, bufout)

      endif

      end associate
   end subroutine project_searchdir

!------------------------------------------------------------------------------------------------------------

   subroutine apply_trcbnd(igs, g, dp, ps)
!--purpose: integrate traction increments and apply traction bound/project onto feasible area
      implicit none
!--subroutine arguments:
      type(t_eldiv),    intent(in) :: igs
      type(t_gridfnc3)             :: g, dp, ps
!--local variables:
      integer      :: ii, ix, iy, mx, my
      real(kind=8) :: ptabs

      associate(eldiv  => igs%el, dpx => dp%vx, dpy => dp%vy,                                   &
                trcbnd => g%vt  , px  => ps%vx, py  => ps%vy )

      mx = ps%grid%nx
      my = ps%grid%ny

      do iy = 1, my
         ii = mx + (iy-1) * mx
         px(ii) = 0d0
         py(ii) = 0d0
         do ix = mx-1, 1, -1
            ii = ix + (iy-1) * mx
            if (eldiv(ii).le.Exter) then
               px(ii)  = 0d0
               py(ii)  = 0d0
            else
               px(ii) = px(ii+1) + dpx(ii)
               py(ii) = py(ii+1) + dpy(ii)
               ptabs = sqrt(px(ii)**2 + py(ii)**2)
               if (eldiv(ii).eq.Slip .or. ptabs.ge.trcbnd(ii)) then
                  px(ii)  = px(ii) * trcbnd(ii) / ptabs
                  py(ii)  = py(ii) * trcbnd(ii) / ptabs
               end if
            end if
            dpx(ii) = px(ii) - px(ii+1)
            dpy(ii) = py(ii) - py(ii+1)
         end do
      end do

      end associate
   end subroutine apply_trcbnd

!------------------------------------------------------------------------------------------------------------

   subroutine solve_elmtrc(igs, mus, tau_c0, k_tau, taucv, taucs, coefs, ps, ss, dupl, igsopt, psopt, ssopt)
!--purpose: compute locally optimized tangential tractions using elmtrc
!           TODO: extension for plasticity, tauc instead of mus*pn
      implicit none
!--subroutine arguments:
      type(t_eldiv),    intent(in)      :: igs
      type(t_gridfnc3), intent(in)      :: mus, ps, ss, dupl, taucv, taucs
      real(kind=8),     intent(in)      :: tau_c0, k_tau, coefs(2,2)
      type(t_eldiv),    intent(inout)   :: igsopt
      type(t_gridfnc3), intent(inout)   :: psopt, ssopt
!--local variables:
      integer,      parameter :: dbgelm = 0
      real(kind=8), parameter :: omegah = 1d0, omegas = 1d0, epselm = 1d-6
      integer      :: ii, mx, my, npot, elnew, iidum
      real(kind=8) :: pr(3), si(2), taucvi, taucsi, dupli(2)
!--functions used :
      integer ix4ii, iy4ii
!--statement functions for computing ix,iy from ii:
      ix4ii(iidum) = mod(iidum-1,mx)+1
      iy4ii(iidum) = (iidum-1)/mx+1

      mx   = ps%grid%nx
      my   = ps%grid%ny
      npot = ps%grid%ntot

      do ii = 1, npot

         pr     = (/ ps%vx(ii), ps%vy(ii), ps%vn(ii) /)
         si     = (/ ss%vx(ii), ss%vy(ii) /)
         dupli  = (/ dupl%vx(ii), dupl%vy(ii) /)
         taucvi = taucv%vt(ii)
         taucsi = taucs%vt(ii)
            
         ! compute solution by updating element ii only

         elnew = igs%el(ii)
         call plstrc(ii, ix4ii(ii), iy4ii(ii), elnew, coefs, epselm, omegah, omegas, pr, mus%vt(ii),    &
                        tau_c0, k_tau, taucvi, taucsi, si, dupli, dbgelm)

         ! return in psopt, ssopt

         igsopt%el(ii) = elnew
         psopt%vx(ii)  = pr(1)
         psopt%vy(ii)  = pr(2)
         ssopt%vx(ii)  = si(1)
         ssopt%vy(ii)  = si(2)

      enddo

   end subroutine solve_elmtrc

!------------------------------------------------------------------------------------------------------------

   subroutine compute_diagscaling(igs, mus, tau_c0, k_tau, taucv, taucs, coefs, ps, ss, dupl, solv,     &
                dscl, idebug, iidbg)
!--purpose: determine diagonal scaling for gdsteady method
!           TODO: extension for plasticity, tauc instead of mus*pn
      implicit none
!--subroutine arguments:
      type(t_eldiv),    intent(in)      :: igs
      type(t_gridfnc3)                  :: mus, taucv, taucs, ps, ss, dupl, dscl
      type(t_solvers),  intent(in)      :: solv
      real(kind=8),     intent(in)      :: tau_c0, k_tau, coefs(2,2)
      integer,          intent(in)      :: idebug, iidbg
!--local variables:
      real(kind=8),     parameter :: tiny_err = 1d-6, epselm = 1d-6
      real(kind=8),     parameter :: omegah = 1d0, omegas = 1d0, omgscl = 1.00d0
      character(len=1), parameter :: aset(0:3) = (/ 'E', 'H', 'S', 'P' /)
      integer      :: ii, ix, iy, jx, mx, my, npot, elnew, dbgelm, iidum
      real(kind=8) :: th_p0, th_p1, th_s0, th_s1, err_0, upd_p, upd_s, fac_s, dscl_new,                 &
                      pr(3), si(2), dupli(2), taucvi, taucsi
!--functions used :
      integer ix4ii, iy4ii
!--statement functions for computing ix,iy from ii:
      ix4ii(iidum) = mod(iidum-1,mx)+1
      iy4ii(iidum) = (iidum-1)/mx+1

      associate(eldiv => igs%el, px    => ps%vx,   py    => ps%vy,   pn   => ps%vn,                   &
                mus   => mus%vt, sx    => ss%vx,   sy    => ss%vy,   dscl => dscl%vn,                 &
                                 duplx => dupl%vx, duply => dupl%vy)

      mx   = ps%grid%nx
      my   = ps%grid%ny
      npot = ps%grid%ntot

      ! for slipping elements, determine fraction resolved by A*p (fac_s) versus change of p (1-fac_s)
      ! in adhesion, determine reduction factor (fac_h) for elements close to internal boundary

      do ii = 1, npot

         if (eldiv(ii).eq.Slip) then

            pr     = (/ px(ii), py(ii), ps%vn(ii) /)
            si     = (/ sx(ii), sy(ii) /)
            dupli  = (/ duplx(ii), duply(ii) /)
            taucvi = taucv%vt(ii)
            taucsi = taucs%vt(ii)
            
            ! compute initial angle error
            th_p0 = atan2( pr(2),  pr(1))
            th_s0 = atan2(-si(2), -si(1))
            err_0 = th_s0 - th_p0
            if (abs(err_0).ge.pi) err_0 = err_0 - nint(err_0/(2d0*pi)) * 2d0*pi

            ! make perturbation of p if necessary to get some initial error
            if (abs(err_0).le.tiny_err) then
               pr(1) = mus(ii) * pr(3) * cos(th_p0+0.1d0)
               pr(2) = mus(ii) * pr(3) * sin(th_p0+0.1d0)
               si(1) = sx(ii) + coefs(1,1) * (pr(1) - px(ii)) + coefs(1,2) * (pr(2) - py(ii))
               si(2) = sy(ii) + coefs(2,1) * (pr(1) - px(ii)) + coefs(2,2) * (pr(2) - py(ii))

               th_p0 = atan2( pr(2),  pr(1))
               th_s0 = atan2(-si(2), -si(1))
               err_0 = th_s0 - th_p0
               if (abs(err_0).ge.pi) err_0 = err_0 - nint(err_0/(2d0*pi)) * 2d0*pi
            endif

            if (abs(err_0).gt.tiny_err) then
               ! compute solution by updating element ii only
               elnew = eldiv(ii)
               dbgelm = 0
               if (idebug.ge.5 .and. ii.eq.iidbg) dbgelm = 2
               call plstrc(ii, ix4ii(ii), iy4ii(ii), elnew, coefs, epselm, omegah, omegas, pr,          &
                               mus(ii), tau_c0, k_tau, taucvi, taucsi, si, dupli, dbgelm)

               if (elnew.eq.Slip) then

                  ! compute resulting angle error
                  th_p1 = atan2( pr(2),  pr(1))
                  th_s1 = atan2(-si(2), -si(1))

                  upd_p = th_p1 - th_p0
                  if (abs(upd_p).ge.pi) upd_p = upd_p - nint(upd_p/(2d0*pi)) * 2d0*pi
                  upd_s = th_s1 - th_s0
                  if (abs(upd_s).ge.pi) upd_s = upd_s - nint(upd_s/(2d0*pi)) * 2d0*pi

                  ! determine fraction of angle update obtained from updating s
                  fac_s = abs(upd_s) / abs(err_0)

                  if (ii.eq.-150) then
                     write(bufout,'(a,i4,2(a,2f8.3),2(a,f9.4))') 'ii=',ii,': th_p=',th_p0*180d0/pi,     &
                        th_p1*180d0/pi, ', th_s=',th_s0*180d0/pi, th_s1*180d0/pi, '; split',            &
                        1d0-fac_s,' vs', fac_s
                     call write_log(1, bufout)
                  endif

                  if (fac_s.lt.0d0 .or. fac_s.gt.1d0) then
                     write(bufout,'(a,i5,a,f8.4)') 'diag.scaling: at el. ii=',ii,' fac_s=', fac_s
                     call write_log(1, bufout)
                  endif
               endif
            else ! abs(err_0)<=tiny_err
               call write_log('Internal error: diag.scaling: err_0<tiny_err')
            endif
         endif ! ii in Slip

         if (eldiv(ii).eq.Slip .and. elnew.eq.Slip) then

            dscl(ii) = solv%d_slp * (fac_s**solv%pow_s)
            !write(bufout,'(a,i4,4(a,f12.4))') 'ii=',ii,': dscl =',solv%d_slp,'*',fac_s,'^',solv%pow_s, &
            !           '=',dscl(ii)
            !call write_log(1, bufout)

         elseif (eldiv(ii).eq.Adhes .or. elnew.eq.Adhes) then

            ! count number of points to Exterior/Slip downstream of ii
            ! TODO: chi=180deg

            ix = ix4ii(ii)
            jx = ix
            do while(jx.gt.1 .and. eldiv(ii-ix+jx).eq.Adhes)
               jx = jx - 1
            enddo
            if (solv%d_cns.lt.solv%d_ifc) then
               dscl_new = max(solv%d_cns, solv%d_ifc + (ix-jx-1)*min(0d0,solv%d_lin))
            else
               dscl_new = min(solv%d_cns, solv%d_ifc + (ix-jx-1)*max(0d0,solv%d_lin))
            endif
            dscl(ii) = min(50d0*dscl(ii), omgscl * dscl_new + (1d0 - omgscl) * dscl(ii))

         else
            dscl(ii) = 0d0
         endif

      enddo ! ii=1,npot

      if (idebug.ge.1 .and. iidbg.ge.1 .and. iidbg.le.npot) then
         write(bufout,'(2(a,i3),3a,f10.6)') ' el(',ix4ii(iidbg),',',iy4ii(iidbg),') in ',               &
                aset(igs%el(iidbg)), ': Dscl=',dscl(iidbg)
         call write_log(1, bufout)
      endif

      if (.false.) then
         do iy = 1, my
            write(bufout,'(20f12.4)') (dscl((iy-1)*mx+ix), ix=1,mx)
            call write_log(1, bufout)
         enddo
      endif

      end associate
   end subroutine compute_diagscaling

!------------------------------------------------------------------------------------------------------------

   subroutine perform_linesearch(igs, g, dscl, dp, ps, ss, r, dv, v, q, told, itgd, alpha0, alpha_j, j)
!--purpose: determine appropriate step-size alpha_j that minimizes new residual |r-alpha*q|^2
      implicit none
!--subroutine arguments:
      integer,          intent(in)      :: itgd
      integer,          intent(out)     :: j
      real(kind=8),     intent(in)      :: alpha0
      real(kind=8),     intent(inout)   :: alpha_j
      type(t_eldiv),    intent(in)      :: igs
      type(t_gridfnc3)                  :: g, dscl, dp, ps, ss, r, dv, v, q, told
!--local variables:
      integer,      parameter :: max_j = 99
      real(kind=8), parameter :: tiny = 1d-12, eps_j = 0.001d0
      logical      :: stop_j, ldone, has_bracket, prev_bracket, use_parab
      integer      :: ii, npot, elnew, ibrack, itb, tbl_j(max_j), idebug,                               &
                      itb_a, itb_b, itb_x, itb_w, itb_v, itb_u
      real(kind=8) :: theta, nxnew, nynew, txnew, tynew, qn, qt, rt, sn, st, vt,                        &
                      c0adh, c1adh, c0slp, c1slp, c2slp, rhonew, dtheta_j, dth_da, drho_da,             &
                      alpha_prv, dalph(2), ddrho
      real(kind=8) :: brent_p, brent_q, brent_r, dst_cur, dst_prv1, dst_prv2, xmid, tolx
      real(kind=8) :: tbl_alpha(max_j), tbl_rho(max_j), tbl_drhoda(max_j)
      type(t_gridfnc3) :: dpnew, psnew
      character(len=6) :: strpts

      associate(eldiv  => igs%el,  txold => told%vx, tyold => told%vy, vx => v%vx,    vy => v%vy,       &
                trcbnd => g%vt,    px    => ps%vx,   py    => ps%vy,   sx => ss%vx,   sy => ss%vy,      &
                dscl   => dscl%vn, rx    => r%vx,    ry    => r%vy,    qx => q%vx,    qy => q%vy )

      npot = ps%grid%ntot

      call gf3_copy_struc(ps, dpnew, 'linsrch:dpnew', .true.)
      call gf3_copy_struc(ps, psnew, 'linsrch:psnew', .true.)

      idebug    = 0
      if (itgd.eq.-1) idebug = 5

      ! perform loop to determine optimal alpha

      j           =  0
      alpha_j     =  0d0
      alpha_prv   = -1d0
      stop_j      = .false.
      has_bracket = .false.

      do while(.not.stop_j)

         ! compute new function value f(alpha) and derivative f'(alpha)

         rhonew = 0d0
         c0adh  = 0d0
         c1adh  = 0d0
         c0slp  = 0d0
         c1slp  = 0d0
         c2slp  = 0d0

         call gf3_copy(AllElm, dp, dpnew, ikTANG)
         call gf3_axpy(AllInt, alpha_j, dv, dpnew, ikTANG)
         call apply_trcbnd(igs, g, dpnew, psnew)

         do ii = 1, npot

            ! compute angle \Delta\theta and derivative d\theta/d\alpha at current alpha^j

            vt       = vx(ii) * txold(ii) + vy(ii) * tyold(ii)
            dtheta_j = atan( alpha_j * vt / max(tiny,trcbnd(ii)) )
            dth_da   = (trcbnd(ii) * vt) / (max(tiny,trcbnd(ii)**2) + alpha_j**2 * vt**2 )

            ! compute predicted n, t at alpha_j

            theta  = atan2(psnew%vy(ii), psnew%vx(ii))
            nxnew  =  cos(theta)
            nynew  =  sin(theta)
            txnew  = -nynew
            tynew  =  nxnew
            ! if (itgd.eq.2 .and. j.eq.1 .and. ii.eq.2) write(*,*) 'ii=',ii,':', psnew%vx(ii), 
            !                                                                    psnew%vy(ii), theta

            ! release constraint where slip gets positive

            sn     = (sx(ii) * nxnew + sy(ii) * nynew)
            elnew  = eldiv(ii)
            if (elnew.eq.Slip .and. sn.gt.0d0) elnew = Adhes

            ! decompose s, q into normal/tangential for new orientation

            st     = sx(ii) * txnew + sy(ii) * tynew
            rt     = -st
            qn     = qx(ii) * nxnew + qy(ii) * nynew
            qt     = qx(ii) * txnew + qy(ii) * tynew

            ! determine coefficients drho/da = c0adh + c0slp + alpha_j * (c1adh + c1slp) + alpha_j^2 * c2slp

            if (elnew.eq.Adhes) then
               c0adh  = c0adh  + 2d0 * dscl(ii)**2 * (sx(ii) * qx(ii) + sy(ii) * qy(ii))
               c1adh  = c1adh  + 2d0 * dscl(ii)**2 * (qx(ii) * qx(ii) + qy(ii) * qy(ii))
               rhonew = rhonew +       dscl(ii)**2 * ( (sx(ii) + alpha_j*qx(ii))**2 +                   &
                                                       (sy(ii) + alpha_j*qy(ii))**2 )
               if (itgd.eq.2 .and. j.eq.-1) write(*,*) 'ii=',ii,': rhonew=',rhonew,      &
                                       dscl(ii)**2 * (sx(ii) + alpha_j*qx(ii))**2,      &
                                       dscl(ii)**2 * (sy(ii) + alpha_j*qy(ii))**2 
            elseif (elnew.eq.Slip) then
               c0slp  = c0slp  - 2d0 * dscl(ii)**2 * (  rt * qt - rt * sn * dth_da )
               c1slp  = c1slp  - 2d0 * dscl(ii)**2 * ( -qt * qt + qt * sn * dth_da - rt * qn * dth_da )
               c2slp  = c2slp  - 2d0 * dscl(ii)**2 * (                               qt * qn * dth_da )
               rhonew = rhonew + dscl(ii)**2 * (st + alpha_j*qt)**2
               if (itgd.eq.2 .and. j.eq.-1) write(*,*) 'ii=',ii,': rhonew=',rhonew,      &
                                 dscl(ii)**2 * (st + alpha_j*qt)**2
               if (itgd.eq.2 .and. j.eq.-1) write(*,*) 'ii=',ii,':', dscl(ii), st, alpha_j, qt
            endif
         enddo

         if (idebug.ge.10) write(*,*) 'c0adh=',c0adh, c1adh
         drho_da = c0adh + c0slp + alpha_j * (c1adh + c1slp) + alpha_j**2 * c2slp

         ! find position to insert new point in table, keep sorted on alpha
         ! pre: j>=0, tbl_alpha(1:j) can be empty

         itb = 1
         do while(itb.le.j .and. alpha_j.gt.tbl_alpha(itb))
            itb = itb + 1
         enddo
         itb_u = itb
         ! write(*,*) 'insert alpha_j=',alpha_j,' at position itb_u=',itb_u

         ! shift values beyond insert point

         tbl_j     (itb+1:j+1) = tbl_j     (itb:j)
         tbl_alpha (itb+1:j+1) = tbl_alpha (itb:j)
         tbl_rho   (itb+1:j+1) = tbl_rho   (itb:j)
         tbl_drhoda(itb+1:j+1) = tbl_drhoda(itb:j)

         if (has_bracket) then
            if (itb_a.ge.itb_u) itb_a = itb_a + 1
            if (itb_b.ge.itb_u) itb_b = itb_b + 1
            if (itb_x.ge.itb_u) itb_x = itb_x + 1
            if (itb_w.ge.itb_u) itb_w = itb_w + 1
            if (itb_v.ge.itb_u) itb_v = itb_v + 1
         endif

         ! insert new point in table

         j = j + 1
         tbl_j     (itb) = j
         tbl_alpha (itb) = alpha_j
         tbl_rho   (itb) = rhonew
         tbl_drhoda(itb) = drho_da

         ! determine whether a bracket exists at ibrack+[0,1,2] with f(ibrack+1) < f(ibrack),f(ibrack+2)

         prev_bracket = has_bracket
         has_bracket = .false.
         ibrack = 1
         do while(ibrack.lt.j-1 .and. .not.has_bracket)
            has_bracket = (tbl_rho(ibrack+1).le.min(tbl_rho(ibrack), tbl_rho(ibrack+2)))
            if (.not.has_bracket) ibrack = ibrack + 1
         enddo

         if (has_bracket .and. .not.prev_bracket) then

            ! first time a bracket is found: initialize variables for Brent method

            itb_a = ibrack     ! low end of bracket
            itb_x = ibrack + 1 ! minimum function value
            itb_b = ibrack + 2 ! high end of bracket
            if (tbl_rho(itb_a).lt.tbl_rho(itb_b)) then
               itb_w = itb_a ! second best function value
               itb_v = itb_b ! third best function value
            else
               itb_w = itb_b ! second best function value
               itb_v = itb_a ! third best function value
            endif
            dst_prv1 = tbl_alpha(itb_x) - tbl_alpha(itb_w)  ! last step taken
            dst_prv2 = tbl_alpha(itb_w) - tbl_alpha(itb_v)  ! before last step size
            if (idebug.ge.3 .and. j.eq.-3) then
               write(bufout,'(2(a,g12.4),a,l1)') 'dst_prv1=',dst_prv1,', dst_prv2=',dst_prv2
               call write_log(1, bufout)
            endif
         
         elseif (has_bracket) then

            ! Brent method: update pointers a,b for bracket, x,w,v for 1st/2nd/3rd smallest value

            if (tbl_rho(itb_u).le.tbl_rho(itb_x)) then
               ! new overall minimum: update bracket, cycle v : = w := x := u
               if (tbl_alpha(itb_u).ge.tbl_alpha(itb_x)) then
                  itb_a = itb_x
               else
                  itb_b = itb_x
               endif
               itb_v = itb_w
               itb_w = itb_x
               itb_x = itb_u
            else
               ! no new minimum: update bracket & cycle v := w := u or v := u
               if (tbl_alpha(itb_u).lt.tbl_alpha(itb_x)) then
                  itb_a = itb_u
               else
                  itb_b = itb_u
               endif
               if (tbl_rho(itb_u).le.tbl_rho(itb_w) .or. tbl_alpha(itb_w).eq.tbl_alpha(itb_x)) then
                  itb_v = itb_w
                  itb_w = itb_u
               elseif (tbl_rho(itb_u).le.tbl_rho(itb_v) .or. tbl_alpha(itb_v).eq.tbl_alpha(itb_x)       &
                                                        .or. tbl_alpha(itb_v).eq.tbl_alpha(itb_w)) then
                  itb_v = itb_u
               endif
            endif
            dst_prv2 = dst_prv1
            dst_prv1 = dst_cur

         endif

         ! print table

         if (idebug.ge.2 .and. j.eq.8) then
            if (has_bracket) then
               write(bufout,'(2(a,i2),a)') 'A bracket exists at itb=[',ibrack,',',ibrack+2,']'
               call write_log(1, bufout)
            else
               call write_log('table does not have a bracket')
            endif
            do itb = 1, j
               strpts = ' '
               if (has_bracket) then
                  if (itb_a.eq.itb) strpts = 'A'
                  if (itb_b.eq.itb) strpts = trim(strpts) // 'B'
                  if (itb_v.eq.itb) strpts = trim(strpts) // 'V'
                  if (itb_w.eq.itb) strpts = trim(strpts) // 'W'
                  if (itb_x.eq.itb) strpts = trim(strpts) // 'X'
                  if (itb_u.eq.itb) strpts = trim(strpts) // 'U'
               endif
               write(bufout,'(a,i4,a,i2,a,f14.6,a,g18.10,a,g12.4,2a)') 'itb=',itb, ': j=',tbl_j(itb),   &
                   ', alph=',tbl_alpha(itb), ', rho=', tbl_rho(itb), ', drhoda=', tbl_drhoda(itb),      &
                   ' ', trim(strpts)
               call write_log(1, bufout)
            enddo
         endif

         alpha_prv = alpha_j

         ! actual linesearch: determine new alpha_j

         if (has_bracket) then

            ! Brent method for minimization: parabolic interpolation + golden section search

            ! check convergence; tolerance depends on size of point_x

            xmid = 0.5d0 * (tbl_alpha(itb_a) + tbl_alpha(itb_b))
            tolx = eps_j * abs(tbl_alpha(itb_x)) + tiny
            ldone = (abs(tbl_alpha(itb_x)-xmid) .le. 2d0*tolx-0.5d0*(tbl_alpha(itb_b)-tbl_alpha(itb_a)))

            if (idebug.ge.3 .and. j.eq.3) then
               write(bufout,'(2(a,g12.4),a,l1)') 'xmid=',xmid,', tolx=',tolx,', ldone=',ldone
               call write_log(1, bufout)
            endif

            ! determine update for next iteration

            if (ldone) then

               if (idebug.ge.1) call write_log('        brent method: ldone')

            else

               if (abs(dst_prv2).le.tolx) then
                  ! old step too small: use golden section step
                  use_parab = .false.
               else
                  ! prepare for parabolic fit

                  brent_r = (tbl_alpha(itb_x) - tbl_alpha(itb_w)) * (tbl_rho(itb_x) - tbl_rho(itb_v))
                  brent_q = (tbl_alpha(itb_x) - tbl_alpha(itb_v)) * (tbl_rho(itb_x) - tbl_rho(itb_w))
                  brent_p =   (tbl_alpha(itb_x) - tbl_alpha(itb_v)) * brent_q                           &
                            - (tbl_alpha(itb_x) - tbl_alpha(itb_w)) * brent_r
                  brent_q = 2d0 * (brent_q - brent_r)
                  if (brent_q.gt.0d0) brent_p = -brent_p
                  brent_q = abs(brent_q)

                  if (idebug.ge.5 .and. j.eq.8) then
                     write(bufout,'(3(a,g12.4))') 'point_x=',tbl_alpha(itb_x),', point_v=',             &
                                tbl_alpha(itb_v), ', point_w=',tbl_alpha(itb_w)
                     call write_log(1, bufout)
                     write(bufout,'(3(a,g12.4))') ' func_x=',tbl_rho(itb_x),',  func_v=',             &
                                tbl_rho(itb_v), ',  func_w=',tbl_rho(itb_w)
                     call write_log(1, bufout)
                     write(bufout,'(3(a,g12.4))') 'brent_p=',brent_p,', brent_q=',brent_q,              &
                                ', dst_prv2=',dst_prv2
                     call write_log(1, bufout)
                     write(bufout,'(3(a,g12.4))') 'point_a=',tbl_alpha(itb_a),', point_x=',             &
                                tbl_alpha(itb_x), ', point_b=',tbl_alpha(itb_b)
                     call write_log(1, bufout)
                  endif

                  ! check if parabolic fit can be used

                  if (abs(brent_p).ge.abs(0.5*brent_q*dst_prv2) .or.                                    &
                          brent_p.le.brent_q*(tbl_alpha(itb_a)-tbl_alpha(itb_x)) .or.                   &
                          brent_p.ge.brent_q*(tbl_alpha(itb_b)-tbl_alpha(itb_x))) then
                     use_parab = .false.
                  else
                     if (idebug.ge.1) call write_log('        brent method: parabolic fit')
                     use_parab = .true.
                     dst_cur = brent_p / brent_q  ! parabolic step
                     alpha_j = tbl_alpha(itb_x) + dst_cur
                     if (alpha_j-tbl_alpha(itb_a).lt.2d0*tolx .or.                                      &
                                 tbl_alpha(itb_b)-alpha_j.lt.2d0*tolx) then
                        dst_cur = tolx * sign(1d0, xmid-tbl_alpha(itb_x))
                     endif
                  endif
               endif

               if (.not.use_parab) then
                  ! golden section step

                  if (idebug.ge.1) call write_log('        brent method: golden section step')
                  if (tbl_alpha(itb_x).ge.xmid) then
                     dst_prv1 = tbl_alpha(itb_a) - tbl_alpha(itb_x)
                  else
                     dst_prv1 = tbl_alpha(itb_b) - tbl_alpha(itb_x)
                  endif
                  dst_cur = 0.381966d0 * dst_prv1
               endif

               ! Note: original Brent algorithm has safeguards to ensure that step is large enough

               alpha_j = tbl_alpha(itb_x) + dst_cur

            endif ! not ldone

            if (idebug.ge.2 .or. j.ge.max_j-1) then
               write(bufout,'(2(a,i4),a,g12.4,a,i3,a,f14.6)') 'it', itgd, ', j',j,': dalpha = ',        &
                      alpha_j - tbl_alpha(ibrack+1), ',                alpha(',j,') =', alpha_j
               call write_log(1, bufout)
            endif

         elseif (j.eq.1) then

            ! fixed step-size using (r'r)/(r'*D*q)
            if (idebug.ge.1) call write_log('        initial esimate using 0.6*alpha0')
            alpha_j = alpha_j + 0.6d0 * alpha0
            ddrho   = 10d0

            if (idebug.ge.2) then
               write(bufout,'(2(a,i4),a,g12.4,a,i3,a,f14.6)') 'it',itgd, ', j',j,': alpha0 =',alpha0,   &
                                                        ',                 alpha(',j,') =',alpha_j
               call write_log(1, bufout)
            endif

         elseif (tbl_rho(j).lt.tbl_rho(1)) then

            ! monotone rho-values, drho/da<0, smallest rho at end of table
            if (idebug.ge.1) call write_log('        extrapolation at end of table')
            itb      = j - 1
            dalph(1) =  tbl_alpha(itb+1)  - tbl_alpha(itb)
            if (j.lt.3) then
               ddrho    = (tbl_drhoda(itb+1) - tbl_drhoda(itb)) / dalph(1)
               alpha_j  = tbl_alpha(itb+1) + min(3d0*dalph(1), -tbl_drhoda(itb) / ddrho)
            else
               alpha_j  = tbl_alpha(itb+1) + 3d0 * dalph(1)
            endif

            if (idebug.ge.2 .or. j.ge.max_j-1) then
               write(bufout,'(2(a,i4),a,g12.4,a,g12.4,a,i3,a,f14.6)') 'it', itgd, ', j',j, ': dalpha = ', &
                   -tbl_drhoda(itb), ' / ', ddrho, ', alpha(',j,') =', alpha_j
               call write_log(1, bufout)
            endif

         else

            ! monotone rho-values, drho/da>0, smallest rho at start of table
            if (idebug.ge.1) call write_log('        extrapolation at start of table')
            itb      = 1
            dalph(1) =  tbl_alpha(itb+1)  - tbl_alpha(itb)
            if (j.lt.3) then
               ddrho    = (tbl_drhoda(itb+1) - tbl_drhoda(itb)) / dalph(1)
               alpha_j  = tbl_alpha(itb) + min(-3d0*dalph(1), tbl_drhoda(itb) / ddrho)
            else
               alpha_j  = tbl_alpha(itb) - 3d0 * dalph(1)
            endif

            if (idebug.ge.2 .or. j.ge.max_j-1) then
               write(bufout,'(2(a,i4),a,g12.4,a,g12.4,a,i3,a,f14.6)') 'it', itgd, ', j',j, ': dalpha = ', &
                   tbl_drhoda(itb), ' / ', ddrho, ', alpha(',j,') =', alpha_j
               call write_log(1, bufout)
            endif

         endif
         
         ! relative stopping criterion

         stop_j = (j.ge.max_j .or. abs(alpha_j-alpha_prv).lt.eps_j*max(abs(alpha_j),abs(alpha_prv)))

      enddo ! j iteration for alpha

      ! warn if no bracket was found

      if (.not.has_bracket) then
         call write_log(' INTERNAL ERROR. Line-search terminates without finding a bracket.')
         call abort_run()
      endif

      ! if no zero-crossing could be found for drho/dalpha, take minimum value from table

      if (j.ge.max_j) then
         itb = minloc( tbl_rho, 1 )
         write(bufout,'(a,i3,a,g11.3,a,g11.3)') 'no convergence, using j=', tbl_j(itb),                 &
                ' with alpha_j=', tbl_alpha(itb), ', rkr=', tbl_rho(itb)
         call write_log(1, bufout)
         alpha_j = tbl_alpha(itb)
      endif

      end associate

      call gf3_destroy(dpnew)
      call gf3_destroy(psnew)
   end subroutine perform_linesearch

