!------------------------------------------------------------------------------------------------------------
! m_solvpn - solve normal contact problem using NormCG method
!
! Copyright 1993-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_solvpn

use m_hierarch_data
use m_aijpj

implicit none
private

public  normcg
public  normcg_matvec
public  sens_norm

contains

!------------------------------------------------------------------------------------------------------------

   subroutine normcg (ic, geom, npot, dxdy, use_fftprec, maxcg, eps, is_sens, cs, ms, hstot, pen,       &
                      fntrue, igs, ps, itcg, err)
!--purpose: Conjugate Gradient type solver for the system of equations in NORM with fixed/known PENetration
!           or with prescribed force FN. Matrix-free implementation, using routine VecAijPj for evaluating
!           "A * ps" and "A * v", with A represented by Cs.
!           There are npot elements in the potential contact area, but no equations are set up and solved for
      implicit none
!--subroutine arguments:
      type(t_ic),       intent(in)      :: ic
      type(t_geomet),   intent(in)      :: geom
      type(t_inflcf),   intent(in)      :: cs, ms
      type(t_gridfnc3), intent(in)      :: hstot
      integer,          intent(in)      :: npot, maxcg
      logical,          intent(in)      :: use_fftprec, is_sens
      real(kind=8),     intent(in)      :: eps, fntrue, dxdy
      real(kind=8),     intent(inout)   :: pen
      type(t_eldiv),    intent(inout)   :: igs
      type(t_gridfnc3), intent(inout)   :: ps
      integer,          intent(out)     :: itcg
      real(kind=8),     intent(out)     :: err
!--local variables:
      type(t_gridfnc3)   :: rhs, res, z, v, q, dd, r_prv, ps0, tmp1, tmp2
      type(t_eldiv)      :: cpatch
      integer            :: ii, ncon, itinn, itchg
      logical            :: lchg_negpn, lchg_intpen, lchanged
      real(kind=8)       :: alpha, beta, rms_upd, rms_upd1, rms_xk, conv, pn, rz1, rz2, rrprv, rv,     &
                              vav, davg, fk, hsmin, hsmax, htrsh
      integer            :: numinn
      integer, parameter :: idebug = 0

      if (npot.le.150) then
         numinn = 3
      elseif (npot.le.400) then
         numinn = 2
      else
         numinn = 1
      endif
      if (idebug.ge.1) then
         if (use_fftprec) then
            write(bufout,*) 'Using BCCG+ with FFT-preconditioner, ', 'MaxCG=',maxcg,', NumInn=',numinn
         else
            write(bufout,*) 'Using BCCG+ without preconditioner, ', 'MaxCG=',maxcg,', NumInn=',numinn
         endif
         call write_log(1, bufout)
      endif

      ! Initialize element division for scaled interaction between different patches

      if (geom%iplan.eq.4) then
         call eldiv_new(cpatch, ps%grid)
         call eldiv_cpatches(cpatch, geom%npatch, geom%ysep, idebug)
      endif

      ! Initialize grid functions rhs, res, v and q.
      ! Note: rhs is computed in entire potential contact area.
      !       res should be up-to-date in the exterior area at all times.
      !       v is the search direction, which should be zero in the exterior area. q = A * v.

      call gf3_dup(rhs,   'normcg:rhs', ps, .true.)
      call gf3_dup(res,   'normcg:res', ps, .true.)
      call gf3_dup(r_prv, 'normcg:r_prv', ps, .true.)
      call gf3_dup(dd,    'normcg:dd',  ps, .true.)
      call gf3_dup(z,     'normcg:z',   ps, .true.)
      call gf3_dup(v,     'normcg:v',   ps, .true.)
      call gf3_dup(q,     'normcg:q',   ps, .true.)
      call gf3_dup(tmp1,  'normcg:tmp1', ps, .true.)
      call gf3_dup(tmp2,  'normcg:tmp2', ps, .true.)

      ! Save initial estimate for modified stopping criterion

      if (is_sens) then
         call gf3_dup(ps0, 'normcg:ps0', ps, .true.)
         call gf3_copy(AllElm, ps, ps0, ikZDIR)
      endif

      ! rhs = pen - hstot    : Right-hand side of the equations
      ! When the total force is prescribed, the rhs should not include pen
      ! The initial estimate for pen is contained in igs

      if (ic%norm.eq.1) pen = 0d0
      call gf3_set (AllElm, pen, rhs, ikZDIR)
      call gf3_axpy(AllElm, -1.0d0, hstot, rhs, ikZDIR)
      davg = 0d0

      ! Count the number of elements in the actual contact area.

      ncon = 0
      do ii = 1, npot
         if (igs%el(ii).ge.Adhes) ncon = ncon + 1
      enddo

      ! If penetration is prescribed:

      if (ic%norm.eq.0) then

         ! Check if there are points with negative separation hstot

         hsmin = gf3_min(AllElm, hstot, ikZDIR) - pen
         if (hsmin.ge.0d0) then

            ! if not: return zero solution, empty contact area

            !!! write(bufout,'(2a,g10.3,a)') ' WARNING: minimum separation', ' is',                        &
            !!!                              hsmin,' > 0, no contact at all.'
            !!! call write_log(1, bufout)
            call gf3_set(AllElm, 0d0, ps, ikZDIR)
            do ii = 1, npot
               igs%el(ii) = Exter
            enddo
            itcg = 0
            err = 0d0
            goto 9000
         endif

      else

      ! Else: normal force prescribed: 

         ! make sure that the contact area is non-empty

         if (ncon.le.0) then
            hsmin = gf3_min(AllElm, hstot, ikZDIR)
            hsmax = gf3_max(AllElm, hstot, ikZDIR)
            htrsh = hsmin + 0.1d0 * max(hsmax-hsmin, 1d-10)
            do ii = 1, npot
               if (hstot%vn(ii).lt.htrsh) then
                  igs%el(ii) = Adhes
                  ncon = ncon + 1
               endif
            enddo
            call areas(igs)
         endif

         ! scale initial estimate such that it gets the correct force

         fk = dxdy * gf3_sum(AllElm, ps, ikZDIR)
         if (abs(fk).lt.1e-3*fntrue) then
            pn = fntrue / (dxdy*real(ncon))
            call gf3_set(AllInt, pn, ps, ikZDIR)
            !!! write(*,*) 'initial fk=',fk,' small: set constant pn=',pn
            fk = dxdy * gf3_sum(AllElm, ps, ikZDIR)
         else
            call gf3_scal(AllInt, fntrue/fk, ps, ikZDIR)
            !!! write(*,*) 'scaling with fntrue/fk=',fntrue/fk
         endif
      endif

      ! res = rhs - A * ps : Residual of the linear system

      call normcg_matvec (geom, cpatch, igs, AllInt, cs, ps, res, tmp1, tmp2)
      call gf3_scal(AllInt, -1.0d0, res, ikZDIR)
      call gf3_axpy(AllInt, 1.0d0, rhs, res, ikZDIR)

      ! if normal force prescribed: take out the mean residual

      if (ic%norm.eq.1) call gf3_proj_avg(AllInt, res, ikZDIR)
      rz2 = 0d0

      ! This routine uses a relative stop-criterion on the updates 
      ! rms_upd == ||Pk - Pk-1||_rms, which are meaningful to the user. 

      itcg    = 0
      itchg   = 0
      itinn   = 0
      rms_xk  = 1d0
      rms_upd = 2d0*eps*rms_xk
      lchanged = .false.

      ! WHILE ( (Changes to el.div OR rms_upd > eps*rms_xk) ) do

      do while((lchanged .or. rms_upd.gt.eps*rms_xk) .and. itcg.lt.maxcg)

         ! increment iteration counter "it"

         itcg  = itcg + 1
         itinn = itinn + 1

         ! compute preconditioned residual z = M * res

         if (use_fftprec) then
            call VecAijPj(igs, AllInt, z, ikZDIR, res, jkZDIR, ms)
         else
            call gf3_copy(AllElm, res, z, ikZDIR)
         endif

         ! if normal force prescribed: take out the mean residual

         if (ic%norm.eq.1) call gf3_proj_avg(AllInt, z, ikZDIR)

         rz1 = rz2
         rz2 = gf3_dot(AllInt, z, res, ikZDIR)
         if (rz2.lt.0d0) then
            write(bufout,*) 'WARNING: (r,z)=',rz2,' < 0!'
            call write_log(1, bufout)
         endif

         ! Set new search direction v(it)

         if (itcg.le.1 .or. rz1.lt.tiny) then

            ! first iteration: steepest descent direction  v(1) = z(0)

            call gf3_copy(AllInt, z, v, ikZDIR)
         else

            ! use v conjugate to previous search direction(s) v(it) = z(it-1) + beta v(it-1)

            rrprv = gf3_dot(AllInt, z, r_prv, ikZDIR)
            beta = max(0d0, (rz2-rrprv) / max(tiny, rz1))
            call gf3_scal(AllInt, beta, v, ikZDIR)
            call gf3_axpy(AllInt, 1.0d0, z, v, ikZDIR)
         endif

         ! if normal force prescribed: take out the mean update of ps

         if (ic%norm.eq.1) call gf3_proj_avg(AllInt, v, ikZDIR)

         ! q = A v, interior elements only

         call normcg_matvec (geom, cpatch, igs, AllInt, cs, v, q, tmp1, tmp2)

         ! if normal force prescribed: take out the mean update du

         if (ic%norm.eq.1) call gf3_proj_avg(AllInt, q, ikZDIR)

         ! alpha = (r,v) / (v, A * v)
         ! Note: if A is truely pos.def. then vav is always >= 0.
         !       vav == 0 iff v == 0. Avoid division by zero for the
         !       special case of a single element where v == 0 occurs.

         rv    = gf3_dot(AllInt, res, v, ikZDIR)
         vav   = gf3_dot(AllInt, q, v, ikZDIR)

         if (abs(vav).gt.1d-32 .and. ncon.eq.1) then
            alpha = rv / vav
         else
            alpha = rv / max(tiny, vav)
         endif

         ! Update estimate of pressures ps:
         !    ps(it) = ps(it-1) + alpha * v

         call gf3_axpy(AllInt, alpha, v, ps, ikZDIR)

         ! Relative stop-criterion on updates |Pk - Pk-1|:
         !    rms_upd = |Pk - Pk-1| == alpha |v|;  
         !    rms_xk  = |Pk|

         rms_upd = abs(alpha) * gf3_rms(AllInt, v, ikZDIR)
         !!! write(nam, '(a,i3.3)') 'ps', itcg
         !!! call gf3_print(ps, nam, ikZDIR, 3, 12)

         if (itcg.eq.1) rms_upd1 = rms_upd
         if (itcg.le.3 .or. mod(itcg,10).eq.0) then
            if (.not.is_sens) then
               rms_xk  = gf3_rms(AllInt, ps, ikZDIR)
            else
               call gf3_copy(AllElm, ps, tmp1, ikZDIR)
               call gf3_axpy(AllElm, -1d0, ps0, tmp1, ikZDIR)
               rms_xk  = gf3_rms(AllInt, tmp1, ikZDIR)
            endif
         endif

         ! Save residual for use in Polak-Ribiere formula in next iteration

         call gf3_copy(AllElm, res, r_prv, ikZDIR)

         ! Check for negative pressures, move elements to the Exterior area.

         if (itinn.lt.numinn .and. rms_upd.ge.eps*rms_xk) then

            ! In inner iterations, we skip the checks on ps<0 and dd<0.
            ! Update the residual   res(it) = res(it-1) - alpha * q(it)

            call gf3_axpy(AllInt, -alpha, q, res, ikZDIR)

         else

            ! When the inner iteration is complete, check for negative ps

            lchg_negpn = .false.
            do ii = 1, npot
               if (igs%el(ii).ge.Adhes .and. ps%vn(ii).lt.0d0) then
                  if (idebug.ge.5) write(*,*) 'it=',itcg,': moving ','element ii=',ii,' to exterior'
                  igs%el(ii) = Exter
                  ncon = ncon - 1
                  ps%vn(ii) = 0d0
                  lchg_negpn = .true.
               endif
            enddo

            ! If there are no elements in contact anymore,
            ! select the elements with the largest interpenetration.

            if (ncon.le.0) then
               hsmin = gf3_min(AllElm, hstot, ikZDIR)
               do ii = 1, npot
                  if (hstot%vn(ii).le.hsmin+1d-5) then
                     igs%el(ii) = Adhes
                     ncon = ncon + 1
                     ps%vn(ii) = 0d0
                     lchg_negpn = .true.
                  endif
               enddo
            endif

            ! If normal force prescribed: scale pressures such that the correct force is obtained.

            if (ic%norm.eq.1 .and. lchg_negpn) then
               fk = dxdy * gf3_sum(AllElm, ps, ikZDIR)
               if (abs(fk).lt.1e-3*fntrue) then
                  call gf3_set(AllInt, 1d0, ps, ikZDIR)
                  fk = gf3_sum(AllElm, ps, ikZDIR)
               endif
               call gf3_scal(AllInt, fntrue/fk, ps, ikZDIR)
            endif

            if (lchg_negpn) call areas(igs)

            ! Compute dd in whole pot.con, for use in interpen.
            ! test (make sure that res=0 in the Exterior area)

            call normcg_matvec (geom, cpatch, igs, AllElm, cs, ps, dd, tmp1, tmp2)
            call gf3_axpy(AllElm, -1.0d0, rhs, dd, ikZDIR)

            ! if normal force prescribed: take out the mean def.dist

            if (ic%norm.eq.1) then
               davg = gf3_sum(AllInt, dd, ikZDIR) / real(ncon)
               call gf3_proj_avg(AllInt, dd, ikZDIR)
            endif

            call gf3_set(AllElm, 0d0, res, ikZDIR)
            call gf3_axpy(AllInt, -1.0d0, dd, res, ikZDIR)

            ! Check for negative deformed distances (interpenetration) in the exterior area

            lchg_intpen = .false.
            do ii = 1, npot
               if (igs%el(ii).eq.Exter) then

                  ! deformed distance < 0 ?  --> insert into contact area

                  if (dd%vn(ii)-davg.lt.0d0) then
                     igs%el(ii) = Adhes
                     ncon = ncon + 1
                     res%vn(ii) = -(dd%vn(ii)-davg)
                     if (idebug.ge.5) write(*,*) 'it=',itcg, ': moving element ii=',ii,' to interior'
                     lchg_intpen = .true.
                  endif
               endif
            enddo

            do ii = 1, npot
               if (igs%el(ii).eq.Exter) v%vn(ii) = 0d0
            enddo

            if (lchg_intpen) call areas(igs)

            itinn = 0
            lchanged = lchg_intpen .or. lchg_negpn
            if (lchanged) itchg = itchg + 1
            if (idebug.ge.5 .and. lchanged) call wrigs(igs, .false.,0d0)

         endif ! check for negative ps or interpenetration

         if (ic%flow.ge.6) then
            write(bufout,5200) itcg, ncon, npot-ncon, rms_upd, eps*rms_xk
            call write_log(1, bufout)
 5200       format (9x,i3, ', size C, E =',2i6,', |dPn(k)|, ', 'Eps|Pn|:', 2g11.3)
         endif

      ! end-while (not converged)

      enddo

      ! call gf3_print(ps, 'pn', ikZDIR, 4)
      ! call normcg_matvec (.false., geom, cpatch, igs, AllElm, cs, ps, dd, tmp1, tmp2)
      ! call gf3_print(dd, 'un', ikZDIR, 4)
      ! call normcg_matvec (.true., geom, cpatch, igs, AllElm, cs, ps, dd, tmp1, tmp2)

      ! If normal force prescribed: 
      !    obtain the penetration from the average difference h-u
      !    (assuming that pen is not included in rhs)

      if (ic%norm.eq.1) pen = davg

      conv = 1d0
      if (rms_upd*rms_upd1.gt.0d0 .and. itcg.gt.1) then
         conv = exp (dlog(rms_upd/rms_upd1) / (itcg-1))
      endif

!!!   write(*,*) 'NormCG: ItCG=',ItCG,', #act.set=',itchg

      if (rms_upd.gt.rms_xk .and. conv.gt.1d0 .and. itcg.ge.maxcg) then
         write(bufout, 6000) maxcg, rms_upd, conv
         call write_log(2, bufout)
 6000    format (' ERROR. MaxCG= ',i3, ' reached in NormCG.', /, ' |Pk - Pk-1| = ', g12.4,             &
                 '  Estim. avg rate of conv =', g12.4)
         call abort_run()
      endif

      if (ic%flow.ge.5) then
         write(bufout, 6200) itcg, rms_upd, conv
         call write_log(1, bufout)
 6200    format (14x, 'NormCG:', i4, ' iterations. |Pk - Pk-1|,', ' C_est:', g12.4, f8.4)
      endif
      if (ic%flow.ge.9) call wrigs (igs, .false., 0d0)

      err = rms_upd

      ! if (conv.lt.1d0) err = rms_upd * (conv / (1-conv))

      ! end of routine: cleanup

 9000 continue
      call gf3_destroy(rhs)
      call gf3_destroy(res)
      call gf3_destroy(r_prv)
      call gf3_destroy(dd)
      call gf3_destroy(z)
      call gf3_destroy(v)
      call gf3_destroy(q)
      call gf3_destroy(tmp1)
      call gf3_destroy(tmp2)
      if (is_sens) then
         call gf3_destroy(ps0)
      endif
      if (geom%iplan.eq.4) then
         call eldiv_destroy(cpatch)
      endif

   end subroutine normcg

!------------------------------------------------------------------------------------------------------------

   subroutine normcg_matvec (geom, cpatch, igs, iigs, cs, ps, res, tmp1, tmp2)
!--purpose: perform matrix-vector computation for NormCG
      implicit none
!--subroutine arguments:
      type(t_geomet),   intent(in)      :: geom
      type(t_eldiv),    intent(inout)   :: cpatch, igs
      integer,          intent(in)      :: iigs         ! AllInt or AllElm
      type(t_inflcf),   intent(in)      :: cs
      type(t_gridfnc3), intent(inout)   :: ps
      type(t_gridfnc3), intent(inout)   :: res, tmp1, tmp2
!--local variables:
      integer, parameter :: idebug = 0
      integer      :: ip, jp, np, mx, my, npot

      if (geom%iplan.ne.4) then

         ! traditional case, influence coefficients CS * pressures PS

         call VecAijPj(igs, iigs, res, ikZDIR, ps, jkZDIR, cs)

      else

         ! Distinguish sub-patches with reduced interaction on each other

         mx = ps%grid%nx
         my = ps%grid%ny
         npot = mx * my
         np = geom%npatch

         ! initialize result = 0

         call gf3_set(iigs, 0d0, res, ikZDIR)

         do ip = 1, np

            ! get range [yl,yh] for current sub-patch

            ! copy pressures ps to tmp1, zero outside patch ip

            call gf3_set(AllElm, 0d0, tmp1, ikZDIR)
            call gf3_msk_copy(ip, cpatch, ps, tmp1, ikZDIR)

            ! form product tmp2 = cs2 * tmp1

            call VecAijPj(igs, iigs, tmp2, ikZDIR, tmp1, jkZDIR, cs)

            ! scale contributions to other patches

            do jp = 1, np
               if (jp.ne.ip) then
                  call gf3_msk_scal(jp, cpatch, geom%facsep(ip,jp), tmp2, ikZDIR)
               endif
            enddo

            ! add result to res

            call gf3_axpy(iigs, 1d0, tmp2, res, ikZDIR)

         enddo

      endif

   end subroutine normcg_matvec

!------------------------------------------------------------------------------------------------------------

   subroutine eldiv_cpatches(eldiv, npatch, ysep, idebug)
!--function: fill element division with sub-patch numbers on the basis of y-ranges
      implicit none
!--subroutine arguments
      type(t_eldiv)             :: eldiv
      integer,      intent(in)  :: npatch, idebug
      real(kind=8), intent(in)  :: ysep(npatch-1)
!--local variables
      integer      :: mx, my, ip, ii, ix, iy, iy0, iy1
      logical      :: ldone
      real(kind=8) :: y1, yn

      mx = eldiv%grid%nx
      my = eldiv%grid%ny
      y1 = eldiv%grid%y(1)
      yn = eldiv%grid%y(mx*my)

      if (idebug.ge.2) then
         write(bufout,'(a,i4,2(a,f7.3),a,i2,a,4f7.3)') ' my=',my,', y1=',y1,', yn=', yn, ', npatch=',   &
                npatch,', ysep=',(ysep(ip), ip=1,npatch-1)
         call write_log(1, bufout)
      endif

      do ip = 1, npatch

         ! ysep-values are given in increasing order: determine range [iy0 : iy1)
         ! patch  1: y-values (-\infty   , ysep( 1))
         ! patch ip: y-values [ysep(ip-1), ysep(ip))
         ! patch np: y-values [ysep(np-1),   \infty)

         if (ip.le.1) then
            iy0 = 1
         else
            iy0 = iy1 + 1
         endif

         if (ip.ge.npatch) then
            iy1 = my            ! last segment: use all remaining iy
         elseif (iy0.gt.my) then
            iy1 = my            ! no points remaining: set empty interval
         elseif (eldiv%grid%y(iy0*mx).ge.ysep(ip)) then
            iy1 = iy0 - 1       ! first possible y already beyond upper bound ysep(ip): set empty
         else
            iy1 = iy0
            ldone = .false.
            do while (iy1.lt.my .and. .not.ldone)
               if (eldiv%grid%y((iy1+1)*mx).lt.ysep(ip)) then
                  iy1 = iy1 + 1
               else
                  ldone = .true.
               endif
            enddo
         endif

         if (iy1.lt.iy0 .and. idebug.ge.-1) then
            write(bufout,'(3(a,i3),a)') ' WARNING: sub-patch',ip,': empty range iy = [', iy0,',',iy1, '].'
            call write_log(1, bufout)
         elseif (idebug.ge.2) then
            ! if (min(iy0,iy1).ge.1 .and. max(iy0,iy1).le.my) then
            write(bufout,'(3(a,i3),2(a,f7.3),a)') ' sub-patch',ip,': selecting range [', iy0,',',       &
                   iy1,'], y=[', eldiv%grid%y(iy0*mx),',', eldiv%grid%y(iy1*mx),']'
            call write_log(1, bufout)
         endif

         ! set patch number in element division

         do iy = iy0, iy1
            do ix = 1, mx
               ii = ix + (iy-1)*mx
               eldiv%el(ii) = ip
            enddo
         enddo
      enddo

   end subroutine eldiv_cpatches

!------------------------------------------------------------------------------------------------------------

   subroutine sens_norm(ic, geom, npot, dxdy, use_fftprec, mxsens, epsens, cs, ms, hstot, pen,          &
                           fntrue, igs, ps, unn, sens, itcg)
!--purpose: quickly compute senstitivity d Fn/d pen with < 1% accuracy
      implicit none
!--subroutine arguments:
      type(t_ic),       intent(in)      :: ic
      type(t_geomet),   intent(in)      :: geom
      type(t_inflcf),   intent(in)      :: cs, ms
      type(t_gridfnc3), intent(in)      :: hstot
      integer,          intent(in)      :: npot, mxsens
      logical,          intent(in)      :: use_fftprec
      real(kind=8),     intent(in)      :: epsens, fntrue, dxdy
      real(kind=8),     intent(in)      :: pen
      type(t_eldiv),    intent(in)      :: igs
      type(t_gridfnc3), intent(in)      :: ps, unn
      integer,          intent(out)     :: itcg
      real(kind=8), dimension(:,:), intent(inout) :: sens
!--local variables:
      integer, parameter :: idebug = 0
      type(t_ic)         :: icloc
      type(t_eldiv)      :: igsloc
      type(t_gridfnc3)   :: psloc
      integer            :: ncon, nadh, nslip, nplast, nexter, max_loc
      real(kind=8)       :: pentru, dpen, penloc, fnloc, err, sens_loc, sens_true
      character(len=12)  :: strng

      call gf3_nullify(psloc)
      call gf3_dup(psloc, 'psloc', ps, .true.)
      call eldiv_nullify(igsloc)
      call eldiv_new(igsloc, ps%grid)

      ! copy control digits, set N = 0

      icloc = ic
      icloc%norm = 0
      icloc%flow = ic%flow - 1

      ! compute "true penetration" using maximum displacement unn

      pentru = gf3_max(AllInt, unn, ikZDIR)

      ! assuming Fn = C pen^1.5, dfn \approx dpen^1.5 for target dfn = 0.001, this gives dpen = 0.0067

      dpen   = 0.0067d0 * pentru
      penloc = pen + dpen

      if (idebug.ge.2) then
         write(bufout,*) 'max(unn) = ',pentru,', dpen =',dpen
         call write_log(1, bufout)
      endif

      ! set initial estimate for solver

      call gf3_copy(AllElm, ps, psloc, ikZDIR)
      call eldiv_copy(igs, igsloc, ikALL)

      ! solve system with perturbed pen

      max_loc = mxsens
      call NormCG (icloc, geom, npot, dxdy, use_fftprec, max_loc, epsens, .true., cs, ms,               &
                   hstot, penloc, fnloc, igsloc, psloc, itcg, err)

      ! compute and display results

      fnloc = dxdy * gf3_sum(AllElm, psloc, ikZDIR)
      sens_loc = (fnloc - fntrue) / dpen

      call eldiv_count(igsloc, nadh, nslip, nplast, nexter)
      ncon = nadh + nslip + nplast

      if (idebug.ge.1) then
         write(bufout,110) max_loc, ncon, fnloc, fnloc/fntrue, sens_loc
         call write_log(1, bufout)
 110     format(' maxcg=',i2,': ncon=',i5,', fn=', f10.3,', fn/fn0=',f7.4,', dfn/dpen=',f12.3,:,       &
                ', err=', f7.4,' %')
      endif
      if (icloc%flow.ge.3) then
         strng = fmt_gs(11, 3, sens_loc)
         write(bufout,120) strng, itcg
         call write_log(1, bufout)
 120     format('      Norm: estimated sensitivity: dFn/dpen=',a11, '  ItCG=',i6)
      endif

      if (.false.) then
         do max_loc = 10, 1, -1

            ! set initial estimate for solver

            call gf3_copy(AllElm, ps, psloc, ikZDIR)
            call eldiv_copy(igs, igsloc, ikALL)

            ! solve system with perturbed pen

            call NormCG(icloc, geom, npot, dxdy, use_fftprec, max_loc, epsens, .true., cs,              &
                        ms, hstot, penloc, fnloc, igsloc, psloc, itcg, err)

            ! compute and display results

            fnloc = dxdy * gf3_sum(AllElm, psloc, ikZDIR)
            sens_loc = (fnloc - fntrue) / dpen

            call eldiv_count(igsloc, nadh, nslip, nplast, nexter)
            ncon = nadh + nslip + nplast

            if (max_loc.eq.10) then
               sens_true = sens_loc
            else
               write(bufout,110) max_loc, ncon, fnloc, fnloc/fntrue, sens_loc,                         &
                                 100d0*(sens_loc-sens_true)/sens_true
               call write_log(1, bufout)
            endif

         enddo ! max_loc
      endif

      sens(iout_fn, iin_dpen) = sens_loc
      itcg = 0

      call gf3_destroy(psloc)
      call eldiv_destroy(igsloc)

   end subroutine sens_norm

!------------------------------------------------------------------------------------------------------------

end module m_solvpn

