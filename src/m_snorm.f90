!------------------------------------------------------------------------------------------------------------
! m_snorm - solve normal contact problem (basically, just delegate to m_solvpn)
!
! Copyright 1979-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_snorm

use m_hierarch_data
use m_solvpn
use m_aijpj
use m_hertz

implicit none
private

public  snorm
public  snorm_nmdbg
public  snorm_kpec_wrapper
private snorm_kpec_method
private kpec_equiv_ellipse
private snorm_analyn_method
public  test_fft

contains


!------------------------------------------------------------------------------------------------------------
   subroutine snorm (ic, mater, cgrid, kin, geom, solv, hs, cs, ms, itnorm, igs1, ps1, sens)
!--purpose: Normal problem. Norm calculates the normal pressure ps%vn(ii), for given tangential traction
!           and undeformed distance hs. If ic%norm = 1, the normal displacement is modified by a rigid
!           displacement pen in the z-direction so that the total normal force equals fn.
!           The tangential traction is not affected.
      implicit none
!--subroutine parameters :
      type(t_ic)             :: ic
      type(t_material)       :: mater
      type(t_grid),   target :: cgrid
      type(t_kincns), target :: kin
      type(t_geomet), target :: geom
      type(t_solvers)        :: solv
      type(t_inflcf)         :: cs, ms
      type(t_eldiv)          :: igs1
      type(t_gridfnc3)       :: hs, ps1
      integer                :: itnorm
      real(kind=8), dimension(:,:) :: sens
!--local variables :
      type(t_gridfnc3)      :: htang, hstot, tmp, unn
      logical               :: zready, use_vecfft, use_fftprec
      integer               :: ii, ix, iy, npot, ncon, nconprev, newext, newcon, itcg, it, smllpn, smlldd
      real(kind=8)          :: dd, errpn, tol, dxdy
      real(kind=8), pointer :: pen
      logical, parameter    :: writechg = .true.
      integer, parameter    :: iidbg = 34

      call gf3_nullify(tmp)
      npot =  cgrid%ntot
      dxdy =  cgrid%dxdy
      pen  => kin%pen

      call timer_start(itimer_snorm)

      if (ic%flow.ge.3) call write_log(' NORM')

      ! update the element divisions in row1st, rowlst for use in aijpj

      call areas (igs1)

#if defined WITH_MKLFFT
      use_vecfft = (npot.ge.700 .or. ic%x_nmdbg.ge.5)
      use_fftprec = .true.
      if (.false.) then
         call test_fft (igs1, ps1, cs)
         call abort_run()
      endif
#else
      use_vecfft = .false.
      use_fftprec = .false.
#endif

      ! No FFT's when using numerical influence coefficients (The preconditioner can still be used.)

      if (cs%itypcf.eq.1) use_vecfft = .false.

      ! No preconditioner for visco-elastic materials (M=1) or when using Blanco-correction (C=4)
      ! TODO: investigate preconditioner for case with Blanco-correction

      if (ic%mater.eq.1 .or. mater%gencr_eff.eq.4) use_fftprec = .false.

      ! Create the FFT-based approximate inverse

#if defined WITH_MKLFFT
      if (use_fftprec) then
         call fft_makeprec (3, cs, 3, ms)
      endif
#endif

      ! Create grid functions htang, hstot, tmp, unn
      ! dup: nullify, copy structure and el.div from ps1, initialize at 0:

      call gf3_dup(htang, 'snorm:htang', ps1, .true.)
      call gf3_dup(hstot, 'snorm:hstot', ps1, .true.)
      call gf3_dup(tmp,   'snorm:tmp',   ps1, .true.)
      call gf3_dup(unn,   'snorm:unn',   ps1, .true.)

      ! Initialize fixed part of deformed distance:
      !   hs    == undeformed distance, negative for interpenetration
      !   htang == normal displacement difference due to tangential tractions
      !   hstot == deformed distance, fixed part, hstot := hs + htang

      call gf3_copy(AllElm, hs, hstot, ikZDIR)

      if (cs%nt_cpl) then
         call VecAijPj(igs1, AllElm, htang, ikZDIR, ps1, jkTANG, cs)
         call gf3_axpy(AllElm, 1d0, htang, hstot, ikZDIR)
      else
         call gf3_set(AllElm, 0.d0, htang, ikZDIR) 
      endif

      if (ic%x_nmdbg.ge.5) then
         call snorm_nmdbg ('problem to NORM:', 1, cgrid, igs1, hs, pen, ps1, htang, hstot, unn)
      endif

      ! the main loop :
      ! 20  repeat
      !        zready = true
      !        set up eqs
      !        solve
      !        if (contact conds violated)
      !            zready = false
      !            change contact area
      !        if (precision not reached)
      !            zready = false
      !     until (zready)

      ! This main loop works solely on ps1 and igs1, using influence coefficients from cs to compute
      ! displacement differences.

      itnorm = 0
      itcg = 0
 20   continue
         itnorm = itnorm + 1
         zready = .true.

         ! Count the number of elements in the contact area, regularize Pn

         ncon = 0
         do ii = 1, npot
            if (igs1%el(ii).ge.Adhes) then
               ncon = ncon + 1
            else
               ps1%vn(ii) = 0d0
            endif
         enddo

         if (ic%flow.ge.3) write(bufout, 2200) itnorm, ncon, npot-ncon
         if (ic%flow.ge.3) call write_log(1, bufout)
 2200    format (4x,i4, ', Norm: size of Contact, Exterior :', 2i7)
         if (ic%flow.ge.9) call wrigs (igs1, .false., 0d0)
         call areas (igs1)

         ! Compute new normal tractions ps1 and element division igs1

         if (ic%x_nmdbg.ge.8) then
            call snorm_nmdbg ('problem to normcg:', 2, cgrid, igs1, hs, pen, ps1, htang, hstot, unn)
         endif

         call timer_start(itimer_normcg)
         call normcg (ic, geom, npot, dxdy, use_fftprec, solv%maxgs, solv%eps, .false., cs, ms, hstot,  &
                      pen, kin%fntrue, igs1, ps1, it, errpn)
         call timer_stop(itimer_normcg)

         itcg = itcg + it

         ! Count the number of elements in contact again, print output

         nconprev = ncon
         ncon = 0
         do ii = 1, npot
            if (igs1%el(ii).ge.Adhes) ncon = ncon + 1
         enddo

         if (ic%flow.ge.6 .and. ncon.ne.nconprev) then
            write(bufout, 2200) itnorm, ncon, npot-ncon
            call write_log(1, bufout)
         endif

         if (ic%x_nmdbg.ge.7) then
            call snorm_nmdbg ('solution of normcg:', 2, cgrid, igs1, hs, pen, ps1, htang, hstot, unn)
         endif

         ! Contract the contact

         newext = 0
         smllpn = 0
         do ii = 1, npot
            if (igs1%el(ii).ge.Adhes) then

               ! negative pressure? --> move to exterior, clear pressure

               if (ps1%vn(ii).lt.-errpn) then
                  if (writechg) then
                     write(bufout,*) 'Norm: move element ii=', ii,' to exterior'
                     call write_log(1, bufout)
                  endif
                  igs1%el(ii) = Exter
                  newext = newext + 1
                  ps1%vn(ii) = 0d0
               endif

               ! small pressure w.r.t. error margins? --> count

               if (abs(ps1%vn(ii)).lt.errpn) smllpn = smllpn + 1
            endif
         enddo
         ncon = ncon - newext

         ! Contact changed? --> iterate

         if (newext.gt.0) zready = .false.

         if (geom%iplan.eq.4) smlldd = 0

         if (zready .and. geom%iplan.ne.4) then
            call areas (igs1)

            ! Compute deformation differences due to normal pressures in unn

            if (use_vecfft) then
               call VecAijPj(igs1, AllElm, unn, ikZDIR, ps1, jkZDIR, cs)
            endif

            if (ic%x_nmdbg.ge.7) then
               call snorm_nmdbg ('solution of "contract" phase:', 3, cgrid, igs1, hs, pen, ps1, htang,  &
                                 hstot, unn)
            endif

            ! Compute deformation in middle element
            ! Maximal row sum, safe tolerance

            ix = max(1, cgrid%nx/2)
            iy = max(1, cgrid%ny/2)
            ii = ix + (iy-1)*cgrid%nx
            call gf3_set(AllElm, 1d0, tmp, ikZDIR)
            tol = abs(errpn * AijPj(ii, ikZDIR, tmp, jkZDIR, cs))

            ! Expand the contact area

            smlldd = 0
            newcon = 0
            do ii = 1, npot
               if (igs1%el(ii).eq.Exter .and. hstot%vn(ii)-pen.lt.0d0) then

                  ! Compute deformed distance.

                  if (use_vecfft) then
                     ! note: slower for small problems, marginally faster for 35x35 test-problem
                     dd = hstot%vn(ii) - pen + unn%vn(ii)
                  else
                     dd = hstot%vn(ii) - pen + AijPj(ii, ikZDIR, ps1, jkZDIR, cs)
                  endif

                  ! deformed distance < 0 ?  --> insert into contact area

                  if (dd.lt.-tol) then
                     if (writechg) then
                        write(bufout,*) 'Norm: move element ii=', ii,' to interior'
                        call write_log(1, bufout)
                     endif
                     igs1%el(ii) = Adhes
                     newcon = newcon + 1
                  endif

                  ! small deformed distance ? --> count
                  if (abs(dd).lt.tol) then
                     smlldd = smlldd + 1
                     if (.false.) then
                        write(bufout,'(a,i6,2(a,g12.4))') ' element ii=',ii,': dd=',dd,', tol=',tol
                        call write_log(1, bufout)
                     endif
                  endif
               endif
            enddo

            ! Contact changed? --> iterate

            ncon = ncon + newcon
            if (newcon.ne.0) zready = .false.
         endif

         ! Also iterate if the solver did not fully converge

         if (it.ge.solv%maxgs) zready = .false.

      ! main loop is repeated until contact area is not changed anymore

      if (.not.zready .and. itnorm.lt.solv%maxin) goto 20


      ! if not converged: abort execution

      if (.not.zready) then
         write(bufout, 5800) solv%maxin
         call write_log(2, bufout)
         ! stop
 5800    format (/,' NORM: ERROR. MaxIn = ', i5, ' reached.')
         itnorm = -1
      endif

      ! Regularize normal tractions, move elements to exterior

      do ii = 1, npot
         if (ps1%vn(ii).lt.0d0 .and. igs1%el(ii).ge.Adhes) then
            ps1%vn(ii) = 0d0
            igs1%el(ii)   = Exter
            if (writechg) then
               write(bufout,*) 'Norm: move element ii=', ii,' to exterior'
               call write_log(1, bufout)
            endif
            ncon = ncon - 1
         endif
      enddo

      ! Print final element division and number of iterations of CG

      if (ic%flow.ge.3) write(bufout, 6000) ncon, npot-ncon, itcg
      if (ic%flow.ge.3) call write_log(1, bufout)
 6000 format (6x, 'Norm: final element division: C, E=', 2i7, '  ItCG=', i6)
      solv%itcg = itcg

      ! Print warning for potential errors in element division

      if (smlldd+smllpn.ge.5 .and. (ic%flow.ge.3 .or. ic%x_nmdbg.ge.1)) then
         write(bufout, 6100) smlldd, smllpn
         call write_log(3, bufout)
 6100    format (/,' NORM: WARNING. There are',i6,' elements with ', 'small deformed distance', /,     &
                 22x, 'and',i6,' elements with ', 'small pressure')
      endif

      if (ic%x_nmdbg.ge.4) then
         call snorm_nmdbg ('solution of NORM:', 4, cgrid, igs1, hs, pen, ps1, htang, hstot, unn)
         write(lout,*)
      endif

      ! Ready.

      call areas (igs1)

      if (ic%norm.eq.0) then
         kin%fntrue = dxdy * gf3_sum(AllElm, ps1, ikZDIR)
         kin%fnscal = kin%fntrue / mater%ga
      endif
      call timer_stop(itimer_snorm)

      ! Compute sensitivities when requested (not yet available)

      if (ic%sens.ge.2) then
         call timer_start(itimer_sens)
         if (.not.use_vecfft) then
            call VecAijPj(igs1, AllElm, unn, ikZDIR, ps1, jkZDIR, cs)
         endif
         call sens_norm (ic, geom, npot, dxdy, use_fftprec, solv%mxsens, solv%epsens, cs, ms, hstot, pen, &
                         kin%fntrue, igs1, ps1, unn, sens, it)
         solv%itcg = solv%itcg + it
         call timer_stop(itimer_sens)
      else
         sens(iout_fn, iin_dpen) = 0d0
      endif

      if (ic%flow.ge.3) call write_log(' ')

      call gf3_destroy(unn)
      call gf3_destroy(tmp)
      call gf3_destroy(htang)
      call gf3_destroy(hstot) 
      
   end subroutine snorm

!------------------------------------------------------------------------------------------------------------

   subroutine snorm_nmdbg (messg, itable, cgrid, igs, hs, pen, ps, htang, hstot, unn)
!--purpose: Produce "nmdbg-output" for different phases in the NORM algorithm
      implicit none
!--subroutine parameters :
      type(t_grid),   target :: cgrid
      type(t_eldiv)          :: igs
      type(t_gridfnc3)       :: hs, ps, htang, hstot, unn
      real(kind=8)           :: pen
      integer                :: itable
      character(len=*)       :: messg
!--local variables :
      logical, parameter    :: print_tables = .true.
      integer               :: npot, ii, iidum
!--functions used :
      integer ix4ii, iy4ii
!--statement functions for computing ix,iy from ii:
      ix4ii(iidum) = mod(iidum-1,cgrid%nx)+1
      iy4ii(iidum) = (iidum-1)/cgrid%nx+1

      npot = cgrid%ntot

      ! write message, e.g. describing the phase of the NORM algorithm

      write(lout,*) 'Nmdbg: ',trim(messg)

      ! if "table 1" is requested: problem to NORM

      if (itable.eq.1 .and. print_tables) then
         write(lout,101)
         do ii = 1, npot
            write(lout,102) ix4ii(ii), iy4ii(ii), ii, hs%vn(ii), htang%vn(ii), hstot%vn(ii),           &
                  igs%el(ii), ps%vn(ii)
 101        format('  ix,  iy,   ii,  undef.dst + tang.part =   rhs h* , init C, init.press')
 102        format(2(i4,','),i5,',',3(f11.7,','),i6,',',f11.4)
         enddo
      endif

      ! if "table 2" is requested: problem to/solution of solvpn

      if (itable.eq.2 .and. print_tables) then
         write(lout,201)
         do ii = 1, npot
            write(lout,202) ix4ii(ii), iy4ii(ii), ii, hstot%vn(ii), igs%el(ii), ps%vn(ii)
 201        format('  ix,  iy,   ii,    rhs h* , Ck/Ek,   press.pn')
 202        format(2(i4,','),i5,',',f11.7,',',i6,',',f11.4)
         enddo
      endif

      ! if "table 3" is requested: results after contracting C

      if (itable.eq.3 .and. print_tables) then
         ! Note: htang = A * pt, H* == hstot = h + htang
         write(lout,301)
         do ii = 1, npot
            write(lout,302) ix4ii(ii), iy4ii(ii), ii, hstot%vn(ii)-pen, unn%vn(ii)-htang%vn(ii),       &
                  hs%vn(ii)+unn%vn(ii)-pen, igs%el(ii), ps%vn(ii)
 301     format('  ix,  iy,   ii, rhs h*-pen + contrib.pn = def.dist, Ck/Ek,   press.pn')
 302     format(2(i4,','),i5,',',3(f11.7,','),i6,',',f11.4)
         enddo
      endif

      ! if "table 4" is requested: results of NORM algorithm

      if (itable.eq.4 .and. print_tables) then
         ! Note: htang = A * pt, H* == hstot = h + htang
         write(lout,401)
         do ii = 1, npot
            write(lout,402) ix4ii(ii), iy4ii(ii), ii, hs%vn(ii)-pen, htang%vn(ii),                     &
                  unn%vn(ii)-htang%vn(ii), hs%vn(ii)-pen+unn%vn(ii), igs%el(ii), ps%vn(ii)
 401     format('  ix,  iy,   ii, rhs h-pen + contrib.pt + contrib.pn = def.dist, Ck/Ek, press.pn')
 402     format(2(i4,','),i5,',',4(f11.7,','),i4,',',f11.4)
         enddo
      endif

      ! if "table 5" is requested: results of NORM algorithm, compact form without el.numbers

      if (itable.eq.5 .and. print_tables) then
         ! Note: htang = A * pt, H* == hstot = h + htang
         write(lout,501)
         do ii = 1, npot
            if (igs%el(ii).gt.0) write(lout,502) hs%vn(ii)-pen, htang%vn(ii),                     &
                  unn%vn(ii)-htang%vn(ii), hs%vn(ii)-pen+unn%vn(ii), igs%el(ii), ps%vn(ii)
 501     format('  rhs h-pen + contrib.pt + contrib.pn = def.dist, Ck/Ek, press.pn')
 502     format(4(f11.7,','),i4,',',f11.4)
         enddo
      endif

   end subroutine snorm_nmdbg

!------------------------------------------------------------------------------------------------------------

   subroutine snorm_kpec_wrapper (ic, mater, cgrid, kin, geom, solv, hs, cs, itnorm, igs, ps)
!--purpose: Normal problem. Norm calculates the normal pressure ps%vn(ii), for given tangential traction
!           and undeformed distance hs. If ic%norm = 1, the normal displacement is modified by a rigid
!           displacement pen in the z-direction so that the total normal force equals fn.
!           The tangential traction is not affected.
      implicit none
!--subroutine parameters :
      type(t_ic)             :: ic
      type(t_material)       :: mater
      type(t_grid),   target :: cgrid
      type(t_kincns), target :: kin
      type(t_geomet), target :: geom
      type(t_solvers)        :: solv
      type(t_inflcf)         :: cs
      type(t_eldiv)          :: igs
      type(t_gridfnc3)       :: hs, ps
      integer                :: itnorm
!--local variables :
      integer,      parameter :: idebug = 1
      type(t_gridfnc3)      :: htang, hstot
      integer               :: ii, iimin, ixmin, iymin, ncon, my_timer
      real(kind=8)          :: hsmin, pentot

      associate(mx    => cgrid%nx,    my    => cgrid%ny,    npot  => cgrid%ntot,  dxdy  => cgrid%dxdy,  &
                dx    => cgrid%dx,    dy    => cgrid%dy,    x     => cgrid%x,     y     => cgrid%y,     &
                pen   => kin%pen)

      if (ic%bound.le.5) then
         my_timer = itimer_kpec
      else
         my_timer = itimer_analyn
      endif
      call timer_start(my_timer)

      if (ic%flow.ge.3) call write_log(' NORM-KPEC')
      if (ic%sens.ge.2  .and. ic%flow.ge.-1) call write_log(' Kpec: Sensitivities will not be computed.')
      if (ic%norm.ne.0) then
         call write_log(' ERROR: Kpec: N=1 not supported.')
         call abort_run()
      endif

      ! Create grid functions htang, hstot
      ! dup: nullify, copy structure and el.div from ps, initialize at 0:

      call gf3_dup(htang, 'snorm:htang', ps, .true.)
      call gf3_dup(hstot, 'snorm:hstot', ps, .true.)

      ! update the element divisions in row1st, rowlst for use in aijpj

      call areas(igs)

      ! Compute virtual interpenetration including elastic displ. due to tang. tractions
      ! hstot = hs - pen + htang,   gap-function: negative at interpenetration

      call gf3_set (AllElm, -pen, hstot, ikZDIR)
      call gf3_axpy(AllElm,  1d0, hs, hstot, ikZDIR)
      if (cs%nt_cpl) then
         call VecAijPj(igs, AllElm, htang, ikZDIR, ps, jkTANG, cs)
         call gf3_axpy(AllElm, 1d0, htang, hstot, ikZDIR)
      endif

      ! call gf3_print(hstot, 'hstot', ikZDIR, 5)

      ! determine element with largest (negative) interpenetration

      iimin  = idmin(npot, hstot%vn, 1)
      ixmin  = cgrid%ix(iimin)
      iymin  = cgrid%iy(iimin)
      hsmin  = hstot%vn(iimin)
      pentot = -hsmin

      if (idebug.ge.2) then
         write(bufout,'(2(a,f12.6),a,i6,2(a,i3),a)') ' pen =', pen,', min.gap = ',hsmin,' at ii=',      &
                iimin,' = (',ixmin,',',iymin,')'
         call write_log(1, bufout)
      endif

      ! Check if there are points with virtual interpenetration

      if (hsmin.ge.0d0) then

         ! no interpenetration: return zero solution, empty contact area

         call gf3_set(AllElm, 0d0, ps, ikZDIR)
         call eldiv_exter(igs)

      else 

         ! some interpenetration exists

         if (ic%bound.le.5) then
            ! B = 5: KPEC method
            call snorm_kpec_method(mater, cgrid, iimin, hstot, pentot, cs, igs, ps, idebug)
         else
            ! B = 6: ANALYN method
            call snorm_analyn_method(mater, cgrid%nx, cgrid%ny, cgrid%dx, cgrid%dy, igs%el,             &
                        cgrid%x, pentot, hstot%vn, ps%vn, idebug)
         endif

      endif

      ! Count the number of elements in the actual contact area.

      ncon = 0
      do ii = 1, npot
         if (igs%el(ii).ge.Adhes) ncon = ncon + 1
      enddo

      ! Print final element division 

      itnorm = 0
      solv%itcg = 0

      if (ic%flow.ge.3) write(bufout, 6000) ncon, npot-ncon
      if (ic%flow.ge.3) call write_log(1, bufout)
 6000 format (6x, 'Norm: final element division: C, E=', 2i7)

      ! Ready.

      call areas (igs)

      if (ic%norm.eq.0) then
         kin%fntrue = dxdy * gf3_sum(AllElm, ps, ikZDIR)
         kin%fnscal = kin%fntrue / mater%ga
      endif

      if (ic%flow.ge.3) call write_log(' ')

      call gf3_destroy(htang)
      call gf3_destroy(hstot) 

      call timer_stop(my_timer)
      end associate
      
   end subroutine snorm_kpec_wrapper

!------------------------------------------------------------------------------------------------------------

   subroutine snorm_kpec_method (mater, cgrid, iimin, hstot, pentot, cs, igs, ps, idebug)
!--purpose: Implementation of KPEC method for situation with some penetration.
      implicit none
!--subroutine parameters :
      type(t_material)       :: mater
      type(t_grid),   target :: cgrid
      integer                :: iimin
      real(kind=8)           :: pentot
      type(t_inflcf)         :: cs
      type(t_eldiv)          :: igs
      type(t_gridfnc3)       :: hstot, ps
      integer                :: idebug
!--local variables :
      real(kind=8), parameter :: eps_kp = 0.5d0
      integer               :: ii, ix, iy, iy0, iy1, &
                               nummax, numy
      real(kind=8)          :: hsmin, hsthrs, unmid, pmax, eps_kpec, xtrl, xledg, xmid, ymid,           &
                               xi, yi, leny, lenmax, a1, b1, aa, bb

      associate(mx    => cgrid%nx,    my    => cgrid%ny,    npot  => cgrid%ntot,  dxdy  => cgrid%dxdy,  &
                dx    => cgrid%dx,    dy    => cgrid%dy,    x     => cgrid%x,     y     => cgrid%y   )

      ! some interpenetration exists

      hsmin = -pentot

      ! use Hertz' theory to get shape correction factors

      call kpec_equiv_ellipse (mater, cgrid, hstot, pentot, a1, b1, aa, bb, ymid, idebug)

      ! Set contact area according to variable KPEC threshold

      do iy = 1, my
         ii       = 1 + (iy-1) * mx
         yi       = y(ii)
         eps_kpec = max(0d0, (a1 * aa**2 - (a1 * aa**2 - b1 * bb**2) * ((yi - ymid) / bb)**2) / pentot)
         hsthrs   = (1d0 - eps_kpec) * hsmin
         if (idebug.ge.3) then
            write(bufout,'(a,i4,2(a,f7.3),a,f10.6)') 'iy=',iy,', yi=',yi,': eps_kpec=',eps_kpec,        &
                ', hsthrs=',hsthrs
            call write_log(1, bufout)
         endif
         do ix = 1, mx
            ii = ix + (iy-1) * mx
            if (hstot%vn(ii).le.hsthrs) then
               igs%el(ii) = Adhes
            else
               igs%el(ii) = Exter
            endif
         enddo
      enddo

      ! determine for each row the first/last elements in C

      call areas(igs)

      ! determine the length and width of this second contact area

      iy0    = my + 1
      iy1    = 0
      nummax = 0
      do iy = 1, my
         numy   = max(0, igs%rowlst(iy) - igs%row1st(iy) + 1)
         nummax = max(nummax, numy)
         if (numy.ge.1) then
            iy0 = min(iy0, iy)
            iy1 = max(iy1, iy)
         endif
      enddo

      ! y0     = y(iy0*mx)
      ! y1     = y(iy1*mx)
      lenmax = real(nummax) * dx

      ! determine the pressure distribution with dummy pmax = 1

      do iy = 1, my
         numy = max(0, igs%rowlst(iy) - igs%row1st(iy) + 1)
         if (numy.ge.1) then
            ii       = 1 + (iy-1) * mx
            yi       = y(ii)
            eps_kpec = max(0d0, (a1 * aa**2 - (a1 * aa**2 - b1 * bb**2) * ((yi - ymid) / bb)**2) / pentot)
            hsthrs   = (1d0 - eps_kpec) * hsmin

            ! determine trailing edge position xtrl on the contour htot = hsthrs
            ix = igs%row1st(iy)
            ii = ix + (iy-1) * mx
            if (ix.le.1) then
               xtrl  = x(ii) - 0.5d0 * dx
            else
               xtrl  = x(ii-1) + dx * (hsthrs - hstot%vn(ii-1)) / (hstot%vn(ii) - hstot%vn(ii-1))
            endif
            if (idebug.ge.5 .or. xtrl.gt.x(ii)) then
               write(bufout,'(a,i4,a,2i4,a,2f7.3,a,2f10.6,a,f7.3)') ' iy=',iy,', ix_trl=',ix-1,ix,   &
                     ', x=', x(ii-1), x(ii),': h=', hstot%vn(ii-1)-hsthrs, hstot%vn(ii)-hsthrs,      &
                     ', xtrl= ', xtrl
               call write_log(1, bufout)
               xtrl = min(x(ii), xtrl)
            endif
            ! determine leading edge position xledg on the contour htot = hsthrs
            ix = igs%rowlst(iy)
            ii = ix + (iy-1) * mx
            if (ix.ge.mx) then
               xledg = x(ii) + 0.5d0 * dx
            else
               xledg = x(ii) + dx * (hsthrs - hstot%vn(ii)) / (hstot%vn(ii+1) - hstot%vn(ii))
            endif
            if (idebug.ge.5 .or. xledg.lt.x(ii)) then
               write(bufout,'(a,i4,a,2i4,a,2f7.3,a,2f10.6,a,f7.3)') ' iy=',iy,', ix_ldg=',ix,ix+1,   &
                     ', x=', x(ii), x(ii+1),': h=', hstot%vn(ii)-hsthrs, hstot%vn(ii+1)-hsthrs,      &
                     ', xledg=', xledg
               call write_log(1, bufout)
               xledg = max(x(ii), xledg)
            endif
            leny  =  xledg - xtrl

            xmid  = 0.5d0 * (xtrl + xledg)
            do ix = igs%row1st(iy), igs%rowlst(iy)
               ii = ix + (iy-1) * mx
               xi = x(ii) - xmid
               if (abs(xi).gt.abs(0.5d0*leny)) then
                  write(bufout,'(2(a,i4),2(a,f10.4))') 'iy=',iy,': ix=',ix,', xi=',xi,               &
                             ', leny/2=',0.5d0*leny
                  call write_log(1, bufout)
               endif
               ps%vn(ii) = 1d0 / (0.5d0*lenmax) * sqrt( (0.5d0*leny)**2 - xi**2 )
            enddo
         endif
      enddo

      ! call gf3_print(ps, 'ps', ikZDIR, 4)

      ! compute deformation at point with maximum penetration

      unmid = AijPj(iimin, ikZDIR, ps, jkZDIR, cs)
      pmax  = -hsmin / unmid

      if (idebug.ge.2) then
         write(bufout,'(2(a,f12.6))') ' un = ', unmid,' --> pmax = ', pmax
         call write_log(1, bufout)
      endif

      ! scale pressures w.r.t. maximum pressure

      call gf3_scal(AllElm, pmax, ps, ikZDIR)

      end associate
      
   end subroutine snorm_kpec_method

!------------------------------------------------------------------------------------------------------------

   subroutine kpec_equiv_ellipse (mater, cgrid, hstot, pentot, a1, b1, aa, bb, ymid, idebug)
!--purpose: determine equivalent ellipse as needed by the KPEC method
      implicit none
!--subroutine parameters :
      type(t_material)       :: mater
      type(t_grid),   target :: cgrid
      type(t_gridfnc3)       :: hstot
      integer                :: idebug
      real(kind=8)           :: a1, b1, aa, bb, ymid, pentot
!--local variables :
      real(kind=8), parameter :: eps_kp = 0.5d0, hz_min_curv = 1d-4
      integer               :: icnt, iofs, ii, ii0, ii1, iimin, iiext, iiint, ix, ix0, ix1, ixmin,      &
                               iy, iy0, iy1, iymin
      logical               :: lfound
      real(kind=8)          :: hsmin, hsthrs, ysta, yend, epshz, e_star, cp, rho, fn_hz

      associate(mx    => cgrid%nx,    my    => cgrid%ny,    npot  => cgrid%ntot,  dxdy  => cgrid%dxdy,  &
                dx    => cgrid%dx,    dy    => cgrid%dy,    x     => cgrid%x,     y     => cgrid%y   ,  &
                h     => hstot%vn  )

      ! determine element with largest (negative) interpenetration

      iimin  = idmin(npot, h, 1)
      ixmin  = cgrid%ix(iimin)
      iymin  = cgrid%iy(iimin)
      hsmin  = h(iimin)
      hsthrs = (1d0 - eps_kp) * hsmin
      if (idebug.ge.3) then
         write(bufout,'(a,f12.6,a,i6,a,2i4,a)') ' h  = ',hsmin,' at ii=',iimin,' = (',ixmin,iymin,')'
         call write_log(1, bufout)
      endif

      ! estimate curvature a1 = A in longitudinal direction

      !  - take three points within row iymin at equal distance from each other
      !    {x0, xmid, x1}, with {h0, h1} in the KP scaled interpenetration area htot <= hstrsh

      ix0    = ixmin - 1
      ix1    = ixmin + 1
      lfound = .false.
      do while(.not.lfound)
         ii0 = ix0 + (iymin-1)*mx
         ii1 = ix1 + (iymin-1)*mx
         lfound = (ix0.le.1 .or. ix1.ge.mx .or. max(h(ii0), h(ii1)).gt.hsthrs)
         if (.not.lfound) then
            ix0 = ix0 - 1
            ix1 = ix1 + 1
         endif
      enddo
      
      if (idebug.ge.3) then
         write(bufout,'(a,3i4,a,3f10.3)') ' selected ix=', ix0, ixmin, ix1,' with h=', h(ii0),          &
                   h(iimin), h(ii1)
         call write_log(1, bufout)
      endif

      !   - compute curvature from parabolic fit

      a1 = (h(ii1) - 2d0 * h(iimin) + h(ii0)) / (2d0 * (x(ii1) - x(iimin))**2)
      if (idebug.ge.3) then
         write(bufout,'(6(a,f9.3),a)') ' curvature a1=',a1,' = (',h(ii1),' - 2 * ',h(iimin),' +',h(ii0), &
                ') / (2 * (', x(ii1),' -',x(iimin),')**2)'
         call write_log(1, bufout)
      endif

      ! estimate curvature b1 = B in lateral direction

      !   - locate first and last rows iy0, iy1 with points in KP interpenetration area

      iy0 = my + 1
      iy1 = 0

      do iy = 1, my

         ! count number of points in KP interpenetration area htot <= hsthrs

         icnt = 0
         do ix = 1, mx
            ii = ix + (iy-1) * mx
            if (h(ii).le.hsthrs) icnt = icnt + 1
         enddo

         ! in case of contact, update iy0, iy1

         if (icnt.ge.1 .and. iy0.gt.my) iy0 = iy
         if (icnt.ge.1) iy1 = max(iy1, iy)
      enddo

      !   - determine start-position y_sta using interpolation between rows iy0-1, iy0

      if (iy0.le.1) then
         if (idebug.ge.2) then
            write(bufout,'(a,i4,a)') ' Warning: kpec_eq_ellipse: pot.contact too small for KP, iy0=',   &
                iy0, ': using extrapolation'
            call write_log(1, bufout)
         endif
         iy0 = iy0 + 1
      endif

      iofs   = (iy0-1) * mx
      ixmin  = idmin(mx, h(iofs+1:), 1)
      iiext  = ixmin + (iy0-2) * mx
      iiint  = ixmin + (iy0-1) * mx
      ysta   = y(iiint) - (hsthrs - h(iiint)) * dy / (h(iiext) - h(iiint))

      if (idebug.ge.3) then
         write(bufout,'(a,2i4,a,i4,a,2f7.3,a,2f10.6,a,f7.3)') ' iy_sta=',iy0-1,iy0,': ix=',ixmin,       &
                        ', y=', y(iiext), y(iiint),': h=', h(iiext)-hsthrs, h(iiint)-hsthrs,            &
                        ', ysta=', ysta
         call write_log(1, bufout)
      endif

      !   - determine end-position y_end using interpolation between rows iy1, iy1+1

      if (iy1.ge.my) then
         if (idebug.ge.2) then
            write(bufout,'(a,i4,a)') ' Warning: kpec_eq_ellipse: pot.contact too small for KP, iy1=',   &
                        iy1, ': using extrapolation'
            call write_log(1, bufout)
         endif
         iy1 = iy1 - 1
      endif

      iofs   = (iy1-1) * mx
      ixmin  = idmin(mx, h(iofs+1:), 1)
      iiint  = ixmin + (iy1-1) * mx
      iiext  = ixmin + (iy1  ) * mx
      yend   = y(iiint) + (hsthrs - h(iiint)) * dy / (h(iiext) - h(iiint))

      if (idebug.ge.3) then
         write(bufout,'(a,2i4,a,i4,a,2f7.3,a,2f10.6,a,f7.3)') ' iy_end=',iy1-1,iy1,': ix=',ixmin,       &
                        ', y=', y(iiint), y(iiext),': h=', h(iiint)-hsthrs, h(iiext)-hsthrs,            &
                        ', yend=', yend
         call write_log(1, bufout)
      endif

      !   - determine curvature from parabolic h(y) = hmin + B (y-ymid)^2
      !                                      hstrsh = hmin + B (yend-ymid)^2

      ymid   = 0.5d0 * (ysta + yend)
      b1     = max(hz_min_curv, (hsthrs-hsmin) / (yend-ymid)**2)

      ! use Hertz' theory to get shape correction factors

      e_star  = mater%ga / (1d0 - mater%nu)
      epshz   = 1d-6
      call hzcalc3d (e_star, epshz, -1, a1, b1, aa, bb, 0, pentot, fn_hz, cp, rho)

      if (idebug.ge.2) then
         write(bufout,'(3(a,f12.6))') ' a1 = ', a1,', b1 =', b1,', pen   =', pentot
         call write_log(1, bufout)
         write(bufout,'(4(a,f12.6))') ' aa = ', aa,', bb =', bb,', aa/bb =', aa/bb,', ymid =',ymid
         call write_log(1, bufout)
      endif

      end associate
   end subroutine kpec_equiv_ellipse

!------------------------------------------------------------------------------------------------------------

   subroutine snorm_analyn_method (mater, mx, my, dx, dy, eldiv, x, pentot, htot, pn, idebug)
!--purpose: Implementation of KPEC method for situation with some penetration.
      implicit none
!--subroutine parameters :
      type(t_material)       :: mater
      integer                :: idebug, mx, my, eldiv(mx,my)
      real(kind=8)           :: pentot, dx, dy, x(mx,my), htot(mx,my), pn(mx,my)
!--local variables :
      real(kind=8), parameter :: hz_min_curv = 1d-4
      integer      :: ic_norm, ipotcn, ixmin, ix, iy, iydbg
      real(kind=8) :: a1_y, b1_y, aa_y, bb_y, alpha1_y, beta1_y, xlc_y, xlen_y, ywid_y,                 &
                      g_y, penloc, pmax_y, fn_hz, e_star, epshz, cp, rho

      iydbg = -18

      ! loop over all rows to set A(y), B(y), a(y), pmax(y), pn(x,y)

      do iy = 1, my

         ! determine contact locus x_lc(iy) and corresponding gap g(y) = h(x_lc,y) (>=0, ex. penetration)

         ixmin  = idmin(mx, htot(1:,iy), 1)
         xlc_y  = x(ixmin,iy)
         g_y    = htot(ixmin,iy) + pentot

         ! determine curvatures A(iy), B(iy)

         if (ixmin.le.1 .or. ixmin.ge.mx) then
            write(bufout,'(a,i5)') ' ERROR (Analyn): found minimum gap at boundary, ix=',ixmin
            call write_log(1, bufout)
            ixmin = max(2, min(mx-1, ixmin))
         endif

         if (idebug.ge.2 .and. (iy.eq.1 .or. iy.eq.my) .and. htot(ixmin,iy).lt.0d0) then
            write(bufout,'(a,i5)') ' Warning (Analyn): found penetration at boundary, iy=',iy
            call write_log(1, bufout)
         endif

         a1_y   = 0.5d0 * (htot(ixmin+1,iy) - 2d0*htot(ixmin,iy) + htot(ixmin-1,iy)) / dx**2
         if (iy.le.1) then
            b1_y   = 0.5d0 * (htot(ixmin,iy+2) - 2d0*htot(ixmin,iy+1) + htot(ixmin,iy)) / dy**2
         elseif (iy.ge.my) then
            b1_y   = 0.5d0 * (htot(ixmin,iy) - 2d0*htot(ixmin,iy-1) + htot(ixmin,iy-2)) / dy**2
         else
            b1_y   = 0.5d0 * (htot(ixmin,iy+1) - 2d0*htot(ixmin,iy) + htot(ixmin,iy-1)) / dy**2
         endif

         if (idebug.ge.5 .or. iy.eq.iydbg) then
            write(bufout,'(a,i4,4(a,f9.6))') ' iy =',iy,': xlc = ',xlc_y,', g_y =',g_y,', a1_y =',     &
                        a1_y,', b1_y =',b1_y
            call write_log(1, bufout)
         endif

         ! replace negative curvature by tiny positive value

         a1_y = max(a1_y, hz_min_curv)
         b1_y = max(b1_y, hz_min_curv)

         ! use Hertz' theory with total pentot to get local semi-axis predictions aa(y), bb(y)

         e_star  = mater%ga / (1d0 - mater%nu)
         epshz   = 1d-6
         ipotcn  = -1
         ic_norm =  0
         call hzcalc3d (e_star, epshz, ipotcn, a1_y, b1_y, aa_y, bb_y, ic_norm, pentot, fn_hz, cp, rho)

         if (idebug.ge.3 .or. iy.eq.iydbg) then
            write(bufout,'(a,i4,6(a,f9.6))') ' iy =', iy,': a1 = ', a1_y,', b1 =', b1_y,', pen   =',   &
                pentot,', aa = ', aa_y,', bb =', bb_y,', aa/bb =', aa_y/bb_y
            call write_log(1, bufout)
         endif

         ! compute shape correction factors: alpha1 = 1+alpha, beta1 = 1+beta 

         alpha1_y = pentot / (a1_y * aa_y**2)
         beta1_y  = pentot / (b1_y * bb_y**2)

         ! determine if contact exists for row iy

         penloc  = pentot - beta1_y * g_y

         if (idebug.ge.3 .or. iy.eq.iydbg) then
            write(bufout,'(a,i4,3(a,f12.6))') ' iy =',iy,': alpha1 = ', alpha1_y,', beta1 =', beta1_y,  &
                        ', penloc =', penloc
            call write_log(1, bufout)
         endif

         if (penloc.le.tiny .or. min(a1_y,b1_y).le.tiny) then

            ! no contact at iy: no penetration or negative curvature

            xlen_y = 0d0
            ywid_y = 0d0

            do ix = 1, mx
               eldiv(ix,iy) = Exter
               pn(ix,iy)    = 0d0
            enddo

         else

            ! contact at iy: positive penetration and curvature

            ! compute local contact length xl(y) and local width yl(y)

            xlen_y = sqrt( penloc / (alpha1_y * a1_y) )
            ywid_y = sqrt( penloc / (beta1_y * b1_y) )

            ! use Hertz' theory to get local maximum pressure pmax(y)

            e_star  = mater%ga / (1d0 - mater%nu)
            epshz   = 1d-6
            ipotcn  = -1
            ic_norm =  0
            call hzcalc3d (e_star, epshz, ipotcn, a1_y, b1_y, aa_y, bb_y, ic_norm, penloc, fn_hz, cp, rho)

            pmax_y = 1.5d0 * fn_hz / (pi * aa_y * bb_y)

            if (idebug.ge.3 .or. iy.eq.iydbg) then
               write(bufout,'(a,i4,3(a,f9.6),a,f9.3)') ' iy =',iy,': penloc = ',penloc,', xlen_y =',    &
                   xlen_y, ', ywid_y =',ywid_y,', pmax_y =',pmax_y
               call write_log(1, bufout)
            endif

            ! set ellipsoidal pressure distribution centered at xlc_y

            do ix = 1, mx
               if (abs(x(ix,iy) - xlc_y).ge.xlen_y) then
                  eldiv(ix,iy) = Exter
                  pn(ix,iy)    = 0d0
               else
                  eldiv(ix,iy) = Adhes
                  pn(ix,iy) = pmax_y * sqrt( max(0d0, 1d0 - ((x(ix,iy) - xlc_y) / xlen_y)**2 ))
               endif
            enddo

         endif ! pen<=0 or b1<=0

      enddo ! iy

   end subroutine snorm_analyn_method

!------------------------------------------------------------------------------------------------------------

#if defined WITH_MKLFFT
   subroutine test_fft (igs, ps, cs)
!--purpose: Test routine.
      implicit none
!--subroutine parameters :
      type(t_inflcf)        :: cs
      type(t_eldiv)         :: igs
      type(t_gridfnc3)      :: ps
!--local variables :
      type(t_gridfnc3)      :: unn
      integer               :: iy

      ! Set tractions and influence grid according to Pieter's test-problem

      ps%grid%nx = 3
      ps%grid%ny = 2
      ps%vn(1:6)          = (/ 0d0, 1d0, 0d0, 0d0, 0d0, 0d0 /)

      cs%cf(-3:2,-2, 1,1) = (/ 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 /)
      cs%cf(-3:2,-1, 1,1) = (/ 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 /)
      cs%cf(-3:2, 0, 1,1) = (/ 0d0, 0d0,-1d0, 2d0, 0d0, 0d0 /)
      cs%cf(-3:2, 1, 1,1) = (/ 0d0, 0d0, 0d0, 1d0, 0d0, 0d0 /)

      do iy=1,2
         write(*,'(a,i2,a,3f12.4)') 'ps row',iy,':', ps%vn( (iy-1)*3+1 : iy*3 )
      enddo

      ! Compute product

      call gf3_dup(unn, 'unn', ps, .true.)
      call fft_VecAijPj(igs, .false., AllElm, unn, ikZDIR, ps, jkZDIR, cs)

      ! Report result

      do iy=1,2
         write(*,'(a,i2,a,3f12.4)') 'unn row',iy,':', unn%vn( (iy-1)*3+1 : iy*3 )
      enddo
      call gf3_destroy(unn)

   end subroutine test_fft
#endif

!------------------------------------------------------------------------------------------------------------

end module m_snorm
