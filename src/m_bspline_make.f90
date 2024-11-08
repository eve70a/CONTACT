!------------------------------------------------------------------------------------------------------------
! m_bspline_make - implementation of B-spline construction
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
submodule (m_bspline) m_bspline_make
   use m_globals
   use m_markers
   use m_ptrarray
   use m_interp
   implicit none

#ifdef _WIN32
   !dir$ attributes c, reference  :: dpbsv, dgbsv    ! on IA-32: use C calling convention for LAPACK routines
#endif

contains

!------------------------------------------------------------------------------------------------------------

module subroutine bspline_make_knot_vector_atmeas(nmeas, xi, use_insert, use_repl, nknot, tx, namcoor)
!--function: make basic knot vector with knots at measurement sites
!            .not.use_insert: skip 2nd and n-1th for not-a-knot b.c.
!            use_repl: replicate 1st and nth knots, .not.use_repl: replicate 1st/last intervals
   implicit none
!--subroutine arguments:
   integer,          intent(in)       :: nmeas
   real(kind=8),     intent(in)       :: xi(nmeas)
   logical,          intent(in)       :: use_insert, use_repl
   integer,          intent(out)      :: nknot
   real(kind=8),     intent(out)      :: tx(:)
   character(len=*), intent(in)       :: namcoor
!--local variables:
   integer                    :: nspl, j
   real(kind=8)    :: dx0, dx1

   nspl   = nmeas       ! number of B-splines equal to number of measured values
   nknot  = nspl + 4    ! cubic spline, order k=4: 4 additional knots
   if (use_insert) nknot = nknot + 2

   if (size(tx,1).lt.nknot) then
      call write_log(' Internal error (make_knot_atmeas): tx too small')
      call abort_run()
   endif

   tx(1:4)           = xi(1)          ! replicate x(1) and x(end) 4 times for order k=4, 
   if (use_insert) then
      tx(5:nknot-4)  = xi(2:nmeas-1)  ! keep all internal x(i)
   else
      tx(5:nknot-4)  = xi(3:nmeas-2)  ! skip x(2) and x(end-1) for not-a-knot boundaries,
   endif
   tx(nknot-3:nknot) = xi(nmeas)

   if (.not.use_repl) then
      dx0   = xi(2) - xi(1)
      tx(1) = xi(1) - 3d0 * dx0
      tx(2) = xi(1) - 2d0 * dx0
      tx(3) = xi(1) - 1d0 * dx0

      dx1   = xi(nmeas) - xi(nmeas-1)
      tx(nknot-2) = xi(nmeas) + 1d0 * dx1
      tx(nknot-1) = xi(nmeas) + 2d0 * dx1
      tx(nknot  ) = xi(nmeas) + 3d0 * dx1
   endif

   if (ldebug.ge.1) then
      write(bufout,'(3a,i4,a,i4)') ' knot vector t',namcoor(1:1),' has ',nknot,' knots, nmeas=',nmeas
      call write_log(1, bufout)
   endif
   if (ldebug.ge.4) then
      do j = 1, nknot
         write(bufout,'(a,i4,3a,f14.6)') ' j =',j,': t',namcoor(1:1),' =',tx(j)
         call write_log(1, bufout)
      enddo
   endif

   end subroutine bspline_make_knot_vector_atmeas

!------------------------------------------------------------------------------------------------------------

module subroutine bspline_make_knot_vector_simple(nprf, sprf, nkink, ikinks, naccel, iaccel, ds_out,    &
                        maxknot, nknot, tj, idebug, my_ierror)
!--function: determine knot-vector with target step ds_out and repeated knots at kinks and accelerations
!            simple version that just applies ds_out as much as possible
   implicit none
!--subroutine arguments:
   integer,      intent(in)     :: nprf, nkink, naccel, maxknot, idebug
   integer,      intent(in)     :: ikinks(nkink), iaccel(naccel)
   real(kind=8), intent(in)     :: sprf(nprf), ds_out
   real(kind=8), intent(out)    :: tj(maxknot)
   integer,      intent(out)    :: nknot, my_ierror
!--local variables:
   logical      :: kink_start, kink_end
   integer      :: i, j, ia, ik, isec, isec_sta, isec_end, nsec, nintv
   real(kind=8) :: s_sec, dt_sec

   my_ierror = 0
   if (idebug.ge.1) call write_log(' using bspline_make_knotvector_simple')

   ! using sections between kinks and acceleration points
   ! note: ikinks includes 1 and nprf

   nsec = 1 + nkink-2 + naccel

   do j = 1, 3
      tj(j) = sprf(1)
      if (idebug.ge.1) then
         write(bufout,'(a,i4,a,f12.6,a)') ' j =',j,': tj =',tj(j),' (start-point)'
         call write_log(1, bufout)
      endif
   enddo

   j    = 3
   ik   = 2 ! next kink, note: ikinks includes 1 and nprf!
   ia   = 1 ! next accel
   isec_end = 1
   kink_end = .false.

   do isec = 1, nsec

      ! section start = end of previous section

      isec_sta   = isec_end
      kink_start = kink_end

      if (isec.ge.nsec) then
         isec_end = nprf
         kink_end = .false.
      elseif (ia.gt.naccel) then
         isec_end = ikinks(ik)
         kink_end = .true.
      elseif (ik.ge.nkink) then
         isec_end = iaccel(ia)
         kink_end = .false.
      elseif (ikinks(ik).lt.iaccel(ia)) then
         isec_end = ikinks(ik)
         kink_end = .true.
      else
         isec_end = iaccel(ia)
         kink_end = .false.
      endif

      ! determine actual step dt_sec for this section \approx ds_out that fits integer multiple times

      s_sec  = sprf(isec_end) - sprf(isec_sta)
      if (isec_end-isec_sta.le.1) then
         nintv = 1      ! double kink or kink next to boundary
      else
         nintv  = max(1, nint(s_sec / ds_out))
      endif
      dt_sec = s_sec / real(nintv)

      if (idebug.ge.2) then
         write(bufout,'(3(a,i3),2(a,l2))') ' section isec=',isec,': sta=',isec_sta,', end=',isec_end,   &
                   ', kink_sta=',kink_start,', end=',kink_end
         call write_log(1, bufout)
         write(bufout,'(18x,2(a,f10.6),a,i3)') ' s_sec =',s_sec,', dt_sec=',dt_sec,', nintv=',nintv
         call write_log(1, bufout)
      endif

      ! add knots for this section including start/end-points

      do i = 0, nintv

         j     = j + 1
         if (i.lt.nintv) then
            tj(j) = sprf(isec_sta) + i * dt_sec
         else
            tj(j) = sprf(isec_end)
         endif

         if (idebug.ge.1) then
            if     (i.eq.0 .and. isec.le.1) then
               write(bufout,'(a,i4,a,f12.6,a)') ' j =',j,': tj =',tj(j),' (start-point)'
            elseif (i.eq.0 .and. kink_start) then
               write(bufout,'(a,i4,a,f12.6,a)') ' j =',j,': tj =',tj(j),' (kink)'
            elseif (i.eq.0) then
               write(bufout,'(a,i4,a,f12.6,a)') ' j =',j,': tj =',tj(j),' (acceleration)'
            elseif (i.eq.nintv .and. isec.ge.nsec) then
               write(bufout,'(a,i4,a,f12.6,a)') ' j =',j,': tj =',tj(j),' (end-point)'
            elseif (i.eq.nintv .and. kink_end) then
               write(bufout,'(a,i4,a,f12.6,a)') ' j =',j,': tj =',tj(j),' (kink)'
            elseif (i.eq.nintv) then
               write(bufout,'(a,i4,a,f12.6,a)') ' j =',j,': tj =',tj(j),' (acceleration)'
            else
               write(bufout,'(a,i4,a,f12.6  )') ' j =',j,': tj =',tj(j)
            endif
            call write_log(1, bufout)
         endif

      enddo

      ! accelerations will be repeated once, at end of isec & start of isec+1
      ! kinks need to be repeated once more

      if (kink_end) then
         j     = j + 1
         tj(j) = tj(j-1)

         if (idebug.ge.1) then
            write(bufout,'(a,i4,a,f12.6,a)') ' j =',j,': tj =',tj(j),' (kink)'
            call write_log(1, bufout)
         endif
      endif

      ! advance ik or ia for next kink or acceleration

      if (kink_end) then
         ik = ik + 1
      else
         ia = ia + 1
      endif

   enddo ! isec

   ! repeat last knot three times

   if (j+3.gt.maxknot) then
      write(bufout,'(2(a,i6))') ' INTERNAL ERROR(make_knot_simple): nknot =',j+3,' > max =',maxknot
      call write_log(1, bufout)
      call abort_run()
   else
      do i = 1, 3
         j = j + 1
         tj(j) = tj(j-1)
         if (idebug.ge.1) then
            write(bufout,'(a,i4,a,f12.6,a)') ' j =',j,': tj =',tj(j),' (end-point)'
            call write_log(1, bufout)
         endif
      enddo
   endif

   nknot = j

end subroutine bspline_make_knot_vector_simple

!------------------------------------------------------------------------------------------------------------

module subroutine bspline_make_knot_vector_advanced(nprf, sprf, nkink, ikinks, naccel, iaccel, ds_out,  &
                        maxknot, nknot, tj, idebug, my_ierror)
!--function: determine knot-vector with target step ds_out and repeated knots at kinks and accelerations
!            advanced knot-vector with at most one knot per data-interval, with measurements skipped for
!            not-a-knot conditions at boundaries, kinks and accelerations
   implicit none
!--subroutine arguments:
   integer,      intent(in)     :: nprf, nkink, naccel, maxknot, idebug
   integer,      intent(in)     :: ikinks(nkink), iaccel(naccel)
   real(kind=8), intent(in)     :: sprf(nprf), ds_out
   real(kind=8), intent(out)    :: tj(maxknot)
   integer,      intent(out)    :: nknot, my_ierror
!--local variables:
   logical      :: kink_start, kink_end, accel_end, skip_2nd, skip_2last, check_next
   integer      :: i, j, ia, ik, iprf, isec, isec_sta, isec_end, nsec, nintv
   real(kind=8) :: s_sec, dt_sec, tj_new

   my_ierror = 0
   if (idebug.ge.1) call write_log(' using bspline_make_knotvector_advanced')

   ! TODO: check unique([ikinks, iaccel])
   ! !!! NOTE: ikinks includes 1 and nprf !!!!

   nsec = 1 + nkink-2 + naccel

   do j = 1, 3
      tj(j) = sprf(1)
      if (idebug.ge.1) then
         write(bufout,'(a,i4,a,f12.6,a)') ' j =',j,': tj =',tj(j),' (start-point)'
         call write_log(1, bufout)
      endif
   enddo

   j    = 3
   ik   = 2 ! next kink
   ia   = 1 ! next accel
   isec_end = 1
   kink_end = .false.

   do isec = 1, nsec

      ! section start = end of previous section

      isec_sta   = isec_end
      kink_start = kink_end

      ! determine new section end

      if (isec.ge.nsec) then
         isec_end  = nprf
         kink_end  = .false.
         accel_end = .false.
      elseif (ia.gt.naccel) then
         isec_end  = ikinks(ik)
         kink_end  = .true.
         accel_end = .false.
      elseif (ik.ge.nkink) then
         isec_end  = iaccel(ia)
         kink_end  = .false.
         accel_end = .true.
      elseif (ikinks(ik).lt.iaccel(ia)) then
         isec_end  = ikinks(ik)
         kink_end  = .true.
         accel_end = .false.
      else
         isec_end  = iaccel(ia)
         kink_end  = .false.
         accel_end = .true.
      endif

      ! determine whether 2nd measurement and before-last measurement should be skipped

      skip_2nd   = (isec.eq.1 .or. kink_start)
      skip_2last = (isec.eq.nsec .or. kink_end .or. accel_end)

      ! determine step dt_sec for this section approximately ds_out

      s_sec  = sprf(isec_end) - sprf(isec_sta)
      nintv  = nint(s_sec / ds_out)
      dt_sec = s_sec / real(nintv)

      if (idebug.ge.1) then
         write(bufout,'(3(a,i4),2(a,f9.3),a,f6.3,a,i6)') ' section',isec,': ip=[',isec_sta,',',isec_end, &
                '], s=[', sprf(isec_sta),',', sprf(isec_end),'], dt=',dt_sec,', nintv=',nintv
         call write_log(1, bufout)
         if (isec.eq.1)    call write_log('              first section - skip 2nd measurement')
         if (kink_start)   call write_log('              kink at start - skip 2nd measurement')
         if (isec.eq.nsec) call write_log('              last section  - skip before-last measurement')
         if (kink_end)     call write_log('              kink at end   - skip before-last measurement')
         if (accel_end)    call write_log('              accel at end  - skip before-last measurement')
      endif

      ! add knots for the section

      j     = j + 1
      tj(j) = sprf(isec_sta)
      iprf  = isec_sta

      if (idebug.ge.1) then
         if (isec.le.1) then
            write(bufout,'(a,i4,a,f12.6,a)') ' j =',j,': tj =',tj(j),' (start-point)'
         elseif (kink_start) then
            write(bufout,'(a,i4,a,f12.6,a)') ' j =',j,': tj =',tj(j),' (kink)'
         else
            write(bufout,'(a,i4,a,f12.6,a)') ' j =',j,': tj =',tj(j),' (acceleration)'
         endif
         call write_log(1, bufout)
      endif

      do while(iprf.lt.isec_end)

         ! take step dt

         tj_new = tj(j) + dt_sec

         ! enlarge step if needed to go to next meas.segment

         tj_new = max(tj_new, sprf(iprf+1))

         ! enlarge step if needed to go to third meas.segment

         if (skip_2nd .and. iprf.eq.isec_sta .and. isec_end-isec_sta.ge.2) then
            tj_new = max(tj_new, sprf(iprf+2))
         endif

         ! if we are inside the last two meas.segments, jump to end of section

         if (skip_2last .and. iprf.ge.isec_end-2) then
            tj_new = sprf(isec_end)
         endif

         ! if step takes us inside the last meas.segment, jump to end of section

         if (skip_2last .and. tj_new.ge.sprf(isec_end-1)) then
            tj_new = sprf(isec_end)
         endif

         ! if step takes us inside before-last meas.segment, reduce step to beginning of this segment

         if (skip_2last .and. isec_end.ge.3) then
            if (tj_new.ge.sprf(isec_end-2) .and. tj_new.le.sprf(isec_end-1)) then
               tj_new = sprf(isec_end-2)
            endif
         endif

         ! never go beyond end of section

         tj_new = min(tj_new, sprf(isec_end))

         ! store new knot

         j     = j + 1
         tj(j) = tj_new

         ! increment iprf such that s(iprf) <= tj_new, skipping over measurements already covered

         check_next = (iprf.lt.isec_end)
         do while(check_next)
            if (sprf(iprf+1).le.tj_new+1d-10) then
               iprf = iprf + 1
               check_next = (iprf.lt.isec_end)
            else
               check_next = .false.
            endif
         enddo

         if (idebug.ge.1) then
            if (iprf.lt.isec_end) then
               write(bufout,'(a,i4,a,f12.6)') ' j =',j,': tj =',tj(j)
            elseif (isec.ge.nsec) then
               write(bufout,'(a,i4,a,f12.6,a)') ' j =',j,': tj =',tj(j),' (end-point)'
            elseif (kink_end) then
               write(bufout,'(a,i4,a,f12.6,a)') ' j =',j,': tj =',tj(j),' (kink)'
            else
               write(bufout,'(a,i4,a,f12.6,a)') ' j =',j,': tj =',tj(j),' (acceleration)'
            endif
            call write_log(1, bufout)
         endif

      enddo ! iprf<isec_end

      ! one more knot at kinks

      if (kink_end) then
         j = j + 1
         tj(j) = tj(j-1)
         if (idebug.ge.1) then
            write(bufout,'(a,i4,a,f12.6,a)') ' j =',j,': tj =',tj(j),' (kink)'
            call write_log(1, bufout)
         endif
      endif

      ! advance ik or ia for next kink or acceleration

      if (kink_end) then
         ik = ik + 1
      else
         ia = ia + 1
      endif

   enddo ! isec

   ! repeat last knot three times

   if (j+3.gt.maxknot) then
      write(bufout,'(2(a,i6))') ' INTERNAL ERROR(make_knot_advanced): nknot =',j+3,' > max =',maxknot
      call write_log(1, bufout)
      call abort_run()
   else
      do i = 1, 3
         j     = j + 1
         tj(j) = tj(j-1)
         if (idebug.ge.1) then
            write(bufout,'(a,i4,a,f12.6,a)') ' j =',j,': tj =',tj(j),' (end-point)'
            call write_log(1, bufout)
         endif
      enddo
   endif

   nknot = j

end subroutine bspline_make_knot_vector_advanced

!------------------------------------------------------------------------------------------------------------

module subroutine bspline_make_knot_vector(nprf, sprf, nkink, ikinks, naccel, iaccel, use_simple,       &
                        ds_out, maxknot, nknot, tj, keepj, nreduc, idebug, my_ierror)
!--function: determine knot-vector with target step ds_out and repeated knots at kinks and accelerations
   implicit none
!--subroutine arguments:
   integer,      intent(in)     :: nprf, nkink, naccel, maxknot, idebug
   integer,      intent(in)     :: ikinks(nkink), iaccel(naccel)
   real(kind=8), intent(in)     :: sprf(nprf), ds_out
   logical,      intent(in)     :: use_simple
   real(kind=8), intent(out)    :: tj(maxknot)
   logical,      intent(out)    :: keepj(maxknot)
   integer,      intent(out)    :: nknot, nreduc, my_ierror
!--local variables:
   integer      :: j, ik, jsta

   my_ierror = 0

   if (use_simple) then
      call bspline_make_knot_vector_simple(nprf, sprf, nkink, ikinks, naccel, iaccel, ds_out,           &
                        maxknot, nknot, tj, idebug, my_ierror)
   else
      call bspline_make_knot_vector_advanced(nprf, sprf, nkink, ikinks, naccel, iaccel, ds_out,         &
                        maxknot, nknot, tj, idebug, my_ierror)
   endif

   ! determine which segments will need linear instead of cubic approximation

   keepj(1:nknot) = .true.
   do ik = 2, nkink

      ! if 'ik' directly adjacent to previous one, find knot at end of segment

      if (ikinks(ik)-ikinks(ik-1).le.1) then

         ! j = first knot with tj(j) >= sprf(ik)

         j = 1
         do while (tj(j).lt.sprf(ikinks(ik)))
            j = j + 1
         enddo

         ! linear segment starts at jsta = j-4, knots jsta+[0:3], mark [1,2] for removal 
         ! e.g. first segment: boundary t1--t4, kink t5--t7 --> remove t2,t3, keep t1,t4,t5--t7

         jsta = j - 4
         keepj(jsta+1) = .false.
         keepj(jsta+2) = .false.

         if (idebug.ge.1) then
            write(bufout,'(2(a,i2),2(a,f12.3),a)') ' kinks at ip=',ikinks(ik-1),',',ikinks(ik),         &
                   ', lin.segment s=[',sprf(ikinks(ik-1)),',',sprf(ikinks(ik)),']'
            call write_log(1,bufout)
            write(bufout,'(a,i4,4(a,f10.3))') ' knots j=',jsta,'+[0:3], tj=',tj(jsta),',',tj(jsta+1),   &
                   ',',tj(jsta+2),',',tj(jsta+3)
            call write_log(1,bufout)
         endif

      endif ! lclose

   enddo ! ik

   ! count number of knots kept in t^*_j

   nreduc = 0
   do j = 1, nknot
      if (keepj(j)) nreduc = nreduc + 1
   enddo

   if (idebug.ge.1 .and. nreduc.ne.nknot) then
      write(bufout,'(2(a,i4),a,8l1)') ' make_knot_vector: nknot=',nknot,', nreduc=',nreduc,             &
                ', keepj(1:8)=', (keepj(j), j=1,8)
      call write_log(1, bufout)
   endif

end subroutine bspline_make_knot_vector

!------------------------------------------------------------------------------------------------------------

module subroutine bspline_make_bt_w_b(nspl, nmeas, jseg, wgt, b4, btwb, idebug)
!--function: compute matrix product B^T * W * B
   implicit none
!--subroutine arguments:
   integer,      intent(in)           :: nspl, nmeas, idebug
   integer,      intent(in)           :: jseg(nmeas)
   real(kind=8), intent(in)           :: wgt(nmeas), b4(nmeas,4)
   real(kind=8), intent(out)          :: btwb(4,nspl)
!--local variables:
   integer              :: i, j, ip, id, ii, jj, ix

   ! compute B^T * W * B; storage of lower diagonals  M(i,j) ==> arr(id,j), with id = i+1-j, i=j:j+3

   btwb(1:4,1:nspl) = 0d0

   do j = 1, nspl
      ! compute lower diagonals i>=j
      do i = j, min(j+3,nspl)
         do ip = 1, nmeas
            ii = i - jseg(ip) + 4
            jj = j - jseg(ip) + 4
            id = i + 1 - j

            if (ii.ge.1 .and. ii.le.4 .and. jj.ge.1 .and. jj.le.4) then
               btwb(id,j) = btwb(id,j) + wgt(ip) * b4(ip,ii) * b4(ip,jj)
            endif
         enddo
      enddo
   enddo

   if (idebug.ge.5) then
      call write_log('matrix B^T W B:')
      do i = 1, min(10,nspl)
         write(bufout,'(i4,1x, 15f9.4)') i, (0d0, j=1,i-4), (btwb(i+1-j,j), j=max(1,i-3),i),            &
                                                                              (0d0, j=i+1,min(10,nspl))
         ix = index(bufout(1), '0.0000')
         do while (ix.gt.0)
            bufout(1)(ix:ix+5) = '     x'
            ix = index(bufout(1), '0.0000')
         enddo
         call write_log(1, bufout)
      enddo

      if (nspl.ge.20) then
         do i = nspl-9, nspl
            write(bufout,'(i4,1x, 15f9.4)') i, (0d0, j=nspl-9,i-4), (btwb(i+1-j,j), j=max(nspl-9,i-3),i),   &
                                                                                      (0d0, j=i+1,nspl)
            ix = index(bufout(1), '0.0000')
            do while (ix.gt.0)
               bufout(1)(ix:ix+5) = '     x'
               ix = index(bufout(1), '0.0000')
            enddo
            call write_log(1, bufout)
         enddo
      endif
   endif

end subroutine bspline_make_bt_w_b

!------------------------------------------------------------------------------------------------------------

module subroutine bspline_make_bt_w_xyz(nspl_all, nspl_rdc, nmeas, jseg, wgt, b4, xi, yi, zi, has_xdata, &
                btw_xyz, idebug)
!--function: compute matrix product B^T * W * [x,y,z]
!            note: btw_xyz has space for all coefficients but is filled for reduced system
   implicit none
!--subroutine arguments:
   integer,      intent(in)           :: nspl_all, nspl_rdc, nmeas, idebug
   integer,      intent(in)           :: jseg(nmeas)
   real(kind=8), intent(in)           :: xi(nmeas), yi(nmeas), zi(nmeas), wgt(nmeas), b4(nmeas,4)
   logical,      intent(in)           :: has_xdata
   real(kind=8), intent(out)          :: btw_xyz(nspl_all,3)
!--local variables:
   integer              :: i, j, ip, jj

   ! compute B^T * W * [x,y,z]

   btw_xyz(1:nspl_all,1:3) = 0d0

   do j = 1, nspl_rdc
      do ip = 1, nmeas
         jj = j - jseg(ip) + 4

         if (jj.ge.1 .and. jj.le.4) then
            if (has_xdata) btw_xyz(j,1) = btw_xyz(j,1) + wgt(ip) * b4(ip,jj) * xi(ip)
            btw_xyz(j,2) = btw_xyz(j,2) + wgt(ip) * b4(ip,jj) * yi(ip)
            btw_xyz(j,3) = btw_xyz(j,3) + wgt(ip) * b4(ip,jj) * zi(ip)
         endif
      enddo
   enddo

   if (idebug.ge.5) then
      call write_log('right hand sides B^T W [x y z]:')
      do i = 1, nspl_rdc
         write(bufout,'(i3,1x, 3f9.4)') i, (btw_xyz(i,j), j=1, 3)
         call write_log(1, bufout)
      enddo
   endif

end subroutine bspline_make_bt_w_xyz

!------------------------------------------------------------------------------------------------------------

module subroutine bspline_make_c_and_d(nknot, tj, tiny_dt, cmat, dmat)
!--function: form the matrices C and D
   implicit none
!--subroutine arguments:
   integer,      intent(in)           :: nknot
   real(kind=8), intent(in)           :: tj(nknot), tiny_dt
   real(kind=8), intent(out)          :: cmat(nknot-1,1), dmat(nknot-1,4)
!--local variables:
   integer              :: i, j, istat
   character(len=20)    :: namvar
   character(len=256)   :: errmsg
   real(kind=8)         :: dtj_k1, dtj_k2, dtj_k3
   real(kind=8), dimension(:),   allocatable :: dtj_inv1, dtj_inv2, dtj_inv3

   namvar = 'dtj_inv1(cd)'
   allocate(dtj_inv1(nknot-1), dtj_inv2(nknot-2), dtj_inv3(nknot-3), stat=istat, errmsg=errmsg)
   if (istat.ne.0) goto 99

   ! determine inverse of average step sizes over 1/2/3 adjacent intervals, zero at zero step size

   do j = 1, nknot-1
      dtj_k1      =  tj(j+1) - tj(j)
      dtj_inv1(j) = 0d0
      if (dtj_k1.gt.tiny_dt) dtj_inv1(j) = 1d0 / dtj_k1
   enddo
   do j = 1, nknot-2
      dtj_k2      = (tj(j+2) - tj(j)) / 2d0
      dtj_inv2(j) = 0d0
      if (dtj_k2.gt.tiny_dt) dtj_inv2(j) = 1d0 / dtj_k2
   enddo
   do j = 1, nknot-3
      dtj_k3      = (tj(j+3) - tj(j)) / 3d0
      dtj_inv3(j) = 0d0
      if (dtj_k3.gt.tiny_dt) dtj_inv3(j) = 1d0 / dtj_k3
   enddo

   ! using sparse storage of C: diagonal matrix and  D = D1 * D2 * D3: diagonal + 3 lower diagonals

   ! determine matrix C: step size dt at true intervals

   cmat(1:nknot-1,1)   = 0d0

   do j = 1, nknot-1
      dtj_k1   = tj(j+1) - tj(j)
      cmat(j,1) = dtj_k1
   enddo

   ! determine matrix D: product of D1 * D2 * D3, each with diagonal + one lower diagonal
         
   dmat(1:nknot-1,1:4) = 0d0

   do i = 4, nknot-4      ! D is all zero in rows i=1:3 and nknot-3:nknot-1

      ! 3rd lower diagonal, element D_{i,i-3}

      dmat(i,1) =  -dtj_inv1(i) * -dtj_inv2(i-1) * -dtj_inv3(i-2)

      ! 2nd lower diagonal, element D_{i,i-2}

      dmat(i,2) =  -dtj_inv1(i) * -dtj_inv2(i-1) *  dtj_inv3(i-2) +                                     &
                  (-dtj_inv1(i) *  dtj_inv2(i-1) +  dtj_inv1(i)   * -dtj_inv2(i)) * -dtj_inv3(i-1)

      ! 1st lower diagonal, element D_{i,i-1}

      dmat(i,3) = (-dtj_inv1(i) *  dtj_inv2(i-1) +  dtj_inv1(i)   * -dtj_inv2(i)) *  dtj_inv3(i-1) +    &
                    dtj_inv1(i) *  dtj_inv2(i)   * -dtj_inv3(i)

      ! diagonal element D_{i,i}

      dmat(i,4) =   dtj_inv1(i) *  dtj_inv2(i)   *  dtj_inv3(i)

   enddo

   return

   ! error handling for memory allocation

99 continue

   write(bufout,'(3a,i4,2a)') ' ERROR in memory allocation (',trim(namvar),':',istat,'): ',trim(errmsg)
   call write_log(1, bufout)
   call abort_run()

end subroutine bspline_make_c_and_d

!------------------------------------------------------------------------------------------------------------

module subroutine bspline_make_dt_c_d(nknot, keepj, nspl_rdc, cmat, dmat, dtcd, idebug)
!--function: compute matrix product D^T * C * D
   implicit none
!--subroutine arguments:
   integer,      intent(in)           :: nknot, nspl_rdc, idebug
   logical,      intent(in)           :: keepj(nknot)
   real(kind=8), intent(in)           :: cmat(nknot-1,1), dmat(nknot-1,4)
   real(kind=8), intent(out)          :: dtcd(4,nspl_rdc)
!--local variables:
   integer              :: i, j, k, id, ii, jj, ix, ir, jr
   integer              :: jreduc(nknot)

   if (idebug.ge.5) then
      call write_log('matrix D:')
      do i = 1, min(10, nspl_rdc)
         write(bufout,'(i3,1x, 10f11.3)') i, (dmat(i,ii), ii=1, 4)
         call write_log(1, bufout)
      enddo
   endif

   ! determine original numbers j for reduced knots jr

   jr = 0
   do j = 1, nknot
      if (keepj(j)) jr = jr + 1
      jreduc(j) = jr
      ! jorig(jr) = j
   enddo

   ! using sparse storage of C: diagonal matrix 
   ! D = D1 * D2 * D3: diagonal + 3 lower diagonals in original numbering j
   !                   storage of lower diagonals  M(i,j) ==> arr(id,j), with id = i+1-j, i=j:j+3

   dtcd(1:4,1:nspl_rdc) = 0d0

   ! compute lower diagonals i>=j

   do k = 4, nknot-4        ! D is all zero in rows k=1:3 and nknot-3:nknot-1
      do jj = 1, 4          ! storage of D(k,j) in dmat(k,j-k+4)
         j  = k - 4 + jj    ! j \in [1, nknot-4]
         jr = jreduc(j)     ! jr \in [1, nreduc-4]
         do ii = jj, 4      ! storage of D(k,i) in dmat(k,i-k+4)
            i = k - 4 + ii  ! i \in [1, nknot-4]
            ir = jreduc(i)  ! ir \in [1, nreduc-4]
            id = ir + 1 - jr ! storage of DTCD(ir,jr) in dtcd(id,jr), id \in [1, 4]
            dtcd(id,jr) = dtcd(id,jr) + cmat(k,1) * dmat(k,ii) * dmat(k,jj)
      enddo
   enddo
   enddo

   if (idebug.ge.5) then
      call write_log('matrix D^T C D:')
      do i = 1, min(10, nspl_rdc)
         write(bufout,'(i4,1x, 10g11.3)') i, (0d0, j=1,i-4), (dtcd(i+1-j,j), j=max(1,i-3),i),           &
                                                                          (0d0, j=i+1,min(10,nspl_rdc))
         ix = index(bufout(1), ' 0.00 ')
         do while (ix.gt.0)
            bufout(1)(ix:ix+5) = '     x'
            ix = index(bufout(1), ' 0.00 ')
         enddo
         call write_log(1, bufout)
      enddo

      if (nspl_rdc.ge.20) then
         do i = nspl_rdc-9, nspl_rdc
            write(bufout,'(i4,1x, 10g11.3)') i, (0d0, j=nspl_rdc-9,i-4),                                &
                                        (dtcd(i+1-j,j), j=max(nspl_rdc-9,i-3),i), (0d0, j=i+1,nspl_rdc)
            ix = index(bufout(1), ' 0.00 ')
            do while (ix.gt.0)
               bufout(1)(ix:ix+5) = '     x'
               ix = index(bufout(1), ' 0.00 ')
            enddo
            call write_log(1, bufout)
         enddo
      endif
   endif

end subroutine bspline_make_dt_c_d

!------------------------------------------------------------------------------------------------------------

module subroutine bspline_reduce_bmat(nknot, keepj, nmeas, jseg, b4)
!--function: compute matrix product B4 * R, combining columns for segments with reduced order,
!            esp. for linear segments between double kinks
   implicit none
!--subroutine arguments:
   integer,      intent(in)           :: nknot, nmeas
   logical,      intent(in)           :: keepj(nknot)
   integer,      intent(inout)        :: jseg(nmeas)
   real(kind=8), intent(inout)        :: b4(nmeas,4)
!--local variables:
   integer      :: j, jr, ip
   integer      :: jreduc(nknot)
   real(kind=8) :: bnew0, bnew1, bnew2, bnew3, bnew4, bnew5

   ! determine new numbers jr for knots j

   jr = 0
   do j = 1, nknot
      if (keepj(j)) jr = jr + 1
      jreduc(j) = jr
   enddo

   ! loop over all measurements ip

   do ip = 1, nmeas
      j = jseg(ip)
      if     (     keepj(j-3) .and.      keepj(j-2) .and.      keepj(j-1) .and.      keepj(j)) then
         !   0: T  T  T  T
         ! all coefficients stay the same, just replace offset jseg -> jr

         jseg(ip) = jreduc(j)

      elseif (     keepj(j-3) .and.      keepj(j-2) .and. .not.keepj(j-1) .and. .not.keepj(j)) then
         !   3: T  T  F  F

         bnew1    = b4(ip,1)
         bnew2    = b4(ip,2) + 2d0/3d0 * b4(ip,3) + 1d0/3d0 * b4(ip,4)
         bnew5    = 0d0      + 1d0/3d0 * b4(ip,3) + 2d0/3d0 * b4(ip,4)
         b4(ip,1) = 0d0
         b4(ip,2) = bnew1
         b4(ip,3) = bnew2
         b4(ip,4) = bnew5
         jseg(ip) = jreduc(j+1)

      elseif (     keepj(j-3) .and. .not.keepj(j-2) .and. .not.keepj(j-1) .and.      keepj(j)) then
         !   6: T  F  F  T

         bnew1    = b4(ip,1) + 2d0/3d0 * b4(ip,2) + 1d0/3d0 * b4(ip,3)
         bnew4    = b4(ip,4) + 1d0/3d0 * b4(ip,2) + 2d0/3d0 * b4(ip,3)
         b4(ip,1) = 0d0
         b4(ip,2) = 0d0
         b4(ip,3) = bnew1
         b4(ip,4) = bnew4
         jseg(ip) = jreduc(j)

      elseif (.not.keepj(j-3) .and.      keepj(j-2) .and.      keepj(j-1) .and.      keepj(j)) then
         !   8: F  T  T  T

         bnew0    = 0d0      + 0d0                + 1d0/3d0 * b4(ip,1)
         bnew2    = b4(ip,2) + 0d0                + 2d0/3d0 * b4(ip,1)
         b4(ip,1) = bnew0
         b4(ip,2) = bnew2
         ! b4(ip,3) = bnew3
         ! b4(ip,4) = bnew4
         jseg(ip) = jreduc(j)

      elseif (.not.keepj(j-3) .and.      keepj(j-2) .and. .not.keepj(j-1) .and. .not.keepj(j)) then
         !  11: F  T  F  F

         bnew0    = 0d0      + 0d0                + 1d0/3d0 * b4(ip,1)
         bnew2    = b4(ip,2) + 0d0                + 2d0/3d0 * b4(ip,1)                                  &
                             + 2d0/3d0 * b4(ip,3) + 1d0/3d0 * b4(ip,4)
         bnew5    = 0d0      + 1d0/3d0 * b4(ip,3) + 2d0/3d0 * b4(ip,4)
         b4(ip,1) = 0d0
         b4(ip,2) = bnew0
         b4(ip,3) = bnew2
         b4(ip,4) = bnew5
         jseg(ip) = jreduc(j+1)

      elseif (.not.keepj(j-3) .and. .not.keepj(j-2) .and.      keepj(j-1) .and.      keepj(j)) then
         !  12: F  F  T  T

         bnew0    = 0d0      + 2d0/3d0 * b4(ip,1) + 1d0/3d0 * b4(ip,2)
         bnew3    = b4(ip,3) + 1d0/3d0 * b4(ip,1) + 2d0/3d0 * b4(ip,2)
         b4(ip,1) = 0d0
         b4(ip,2) = bnew0
         b4(ip,3) = bnew3
         ! b4(ip,4) = bnew4
         jseg(ip) = jreduc(j)

      else
         call write_log(' Internal error(reduce_bmat): unexpected keepj')
         write(bufout,'(a,i4,a,4l2)') ' j=',j,': keepj(j-3:j)=',keepj(j-3:j)
         call write_log(1, bufout)
         call abort_run()
      endif

   enddo

end subroutine bspline_reduce_bmat

!------------------------------------------------------------------------------------------------------------

module subroutine bspline_reduce_dmat(nknot, keepj, dmat, idebug)
!--function: compute matrix product D * R, combining columns for segments with reduced order,
!            esp. for linear segments between double kinks
   implicit none
!--subroutine arguments:
   integer,      intent(in)           :: nknot, idebug
   logical,      intent(in)           :: keepj(nknot)
   real(kind=8), intent(inout)        :: dmat(nknot-1,4)
!--local variables:
   integer    :: j

   ! matrix D = D1*D2*D3 has (nknot-1) x (nknot-4) elements: diagonal + 3 lower diagonals
   ! array dmat(i,k) stores D_{i,i-4+k}

   ! D is zero in rows j for repeated knots t_j = t_{j+1}

   ! D^* = D * R is set to zero in rows j for segments [t_j, t_{j+1}) with reduced order.
   ! for keepj(jsta+[1:4]) = [ T, F, F, T ], these segments are [t_{jsta+3}, t_{jsta+4})

   do j = 1, nknot-1
      if (keepj(j) .and. .not.keepj(j+1)) then
         dmat(j+3,1:4) = 0d0
         if (idebug.ge.5) then
            write(bufout,'(a,i3)') ' reduce_dmat: jsta=',j
            call write_log(1, bufout)
         endif
      endif
   enddo

end subroutine bspline_reduce_dmat

!------------------------------------------------------------------------------------------------------------

module subroutine bspline_expand_vector(nknot, keepj, nreduc, nspl_all, coef, idebug)
!--function: expand input coefficients for reduced basisfunctions coef^* to full coefficients coef
   implicit none
!--subroutine arguments:
   integer,      intent(in)           :: nknot, nreduc, nspl_all, idebug
   logical,      intent(in)           :: keepj(nknot)
   real(kind=8), intent(inout)        :: coef(nspl_all,3)
!--local variables:
   integer    :: jreduc, j, jsta, i

   ! shift rows down to make space for intermediate coefficient values

   jreduc = nreduc - 4

   do j = nknot-4, 1, -1
      if (keepj(j) .and. jreduc.lt.j) then
         coef(j,1:3) = coef(jreduc,1:3)
         jreduc = jreduc - 1
      endif
   enddo

   ! compute intermediate coefficient values

   do j = 1, nknot-4
      if (keepj(j) .and. .not.keepj(j+1)) then
         jsta = j
         coef(jsta+1,1:3) = 2d0/3d0 * coef(jsta,1:3) + 1d0/3d0 * coef(jsta+3,1:3)
         coef(jsta+2,1:3) = 1d0/3d0 * coef(jsta,1:3) + 2d0/3d0 * coef(jsta+3,1:3)
      endif
   enddo

   if (idebug.ge.2) then
      call write_log('Expanded coefficients coef_x, coef_y, coef_z:')
      do i = 1, nspl_all
         if (i.le.10 .or. i.ge.nspl_all-9 .or. idebug.ge.10) then
            write(bufout,'(i3,1x, 3f9.4)') i, (coef(i,j), j=1, 3)
            call write_log(1, bufout)
         elseif (i.eq.11) then
            call write_log('  ...')
         endif
      enddo
   endif

end subroutine bspline_expand_vector

!------------------------------------------------------------------------------------------------------------

module subroutine bspline_solve_coef1d_smooth(nspl_all, nspl_rdc, btwb, lambda, dtcd, btw_xyz, idebug,  &
                                my_ierror)
!--function: solve the linear systems (B^T W B + lambda D^T C D) [ coef_[xyz] ] = B^T W [x, y, z]
!            note: btw_xyz has space for all coeff but is used for reduced system
   implicit none
!--subroutine arguments:
   integer,      intent(in)           :: nspl_all, nspl_rdc, idebug
   real(kind=8), intent(in)           :: lambda
   real(kind=8), intent(in)           :: btwb(4,nspl_rdc), dtcd(4,nspl_rdc)
   real(kind=8), intent(inout)        :: btw_xyz(nspl_all,3)  ! on input:  [BtW*x, BtW*y, BtW*z]
                                                              ! on output: [coef_x, coef_y, coef_z]
   integer,      intent(out)          :: my_ierror
!--local variables:
   integer              :: kl, lda, ldb, nrhs, i, j, istat, info
   character(len=20)    :: namvar
   character(len=256)   :: errmsg
   real(kind=8), dimension(:,:), allocatable :: sysmat

   my_ierror = 0

   kl   = 3         ! number of lower diagonals in system matrix
   lda  = 4         ! leading dimension for system matrix
   nrhs = 3         ! number of right hand sides
   ldb  = nspl_all  ! leading dimension for right hand side

   namvar = 'sysmat(1d)'
   allocate(sysmat(lda,nspl_rdc), stat=istat, errmsg=errmsg)
   if (istat.ne.0) goto 99

   ! form system matrix B^T W B + lambda * D^T C D, storing the diagonal and kl=3 lower diagonals 

   do j = 1, nspl_rdc
                           sysmat(1,j) = btwb(1,j) + lambda * dtcd(1,j)
      if (j+1.le.nspl_rdc) sysmat(2,j) = btwb(2,j) + lambda * dtcd(2,j)
      if (j+2.le.nspl_rdc) sysmat(3,j) = btwb(3,j) + lambda * dtcd(3,j)
      if (j+3.le.nspl_rdc) sysmat(4,j) = btwb(4,j) + lambda * dtcd(4,j)
   enddo

   ! use LAPACK routine for (symmetric) Positive Banded system

   call dpbsv('l', nspl_rdc, kl, nrhs, sysmat, lda, btw_xyz, ldb, info)

   if (idebug.ge.3) then
      write(bufout,*) 'done dpbsv, info=',info
      call write_log(1, bufout)
   endif

   if (idebug.ge.2) then
      call write_log('B-spline coefficients coef_x, coef_y, coef_z:')
      do i = 1, nspl_rdc
         if (i.le.10 .or. i.ge.nspl_rdc-9 .or. idebug.ge.10) then
            write(bufout,'(i5,1x, 3f9.4)') i, (btw_xyz(i,j), j=1, 3)
            call write_log(1, bufout)
         endif
      enddo
   endif

   return

   ! error handling for memory allocation

99 continue

   write(bufout,'(3a,i4,2a)') ' ERROR in memory allocation (',trim(namvar),':',istat,'): ',trim(errmsg)
   call write_log(1, bufout)
   call abort_run()

end subroutine bspline_solve_coef1d_smooth

!------------------------------------------------------------------------------------------------------------

module subroutine bspline_solve_coef1d_intpol(nsplx, jseg, b4, nrow, nsply, cx, cy, cz, has_xdata,      &
                                                                                    idebug, my_ierror)
!--function: solve the linear systems B [cx2d, cy2d, cz2d] = [cx1d, cy1d, cz1d]
   implicit none
!--subroutine arguments:
   integer,      intent(in)           :: nsplx, nrow, nsply, idebug        ! nrow >= nsplx
   integer,      intent(in)           :: jseg(nsplx)
   logical,      intent(in)           :: has_xdata
   real(kind=8), intent(in)           :: b4(nsplx,4)
   real(kind=8), intent(inout)        :: cx(nrow,nsply), cy(nrow,nsply), cz(nrow,nsply)
                                                                ! on input:  c[xyz]1d, on output: c[xyz]2d
   integer,      intent(out)          :: my_ierror
!--local variables:
   integer              :: kl, ku, lda, ldb, nrhs, i, ii, j, jj, jstep, isys, istat, info, ipiv(nsplx)
   character(len=20)    :: namvar
   character(len=256)   :: errmsg
   real(kind=8), dimension(:,:), allocatable :: sysmat

   my_ierror = 0

   ku   = 2             ! number of superdiagonals in system matrix
   kl   = 2             ! number of subdiagonals in system matrix
   lda  = 2*kl + ku + 1 ! leading dimension for system matrix, including kl superdiagonals for pivoting
   nrhs = nsply         ! number of right-hand sides
   ldb  = nrow          ! leading dimension for right-hand sides

   if (idebug.ge.3) then
      write(bufout,'(5(a,i4))') ' ku=',ku,', kl=',kl,', lda=',lda,', nrhs=',nrhs,', ldb=',ldb
      call write_log(1, bufout)
   endif

   namvar = 'sysmat(2d)'
   allocate(sysmat(lda,nsplx), stat=istat, errmsg=errmsg)
   if (istat.ne.0) goto 99

   do isys = 1, 3

      ! form system matrix B, copying the coefficients from b4 to the right positions
      ! own storage b4(nspl,4)        - nonzeros shifted left, columns jj = [1:4] <--> j = jofs + jj - 4
      ! Lapack storage sysmat(7,nspl) - diagonals shifted up, rows ii = [1:5] <--> i = ii - kl - ku - 1 + j

      do j = 1, nsplx
         sysmat(1:lda,j) = 0d0
      enddo

      do i = 1, nsplx
         do jj = 1, 4
            j  = jj - 4 + jseg(i)
            ii = i + kl + ku + 1 - j
            if (ii.ge.kl+1 .and. ii.le.kl+5) sysmat(ii,j) = b4(i,jj)
         enddo
      enddo
            ! i = [-2:2] + j
      sysmat(1:4,1) = -999d0
      sysmat(1:3,2) = -999d0
      sysmat(1:2,3) = -999d0
      sysmat(1:1,4) = -999d0

      if (idebug.ge.5) then
         write(bufout,'(a,i2,a)') 'system matrix for isys=',isys,':'
         call write_log(1, bufout)

         do j = 1, nsplx
            write(bufout,'(a,i4,a,5g12.4)') 'col j=',j,': A=',(sysmat(ii,j), ii = kl+1, kl+5)
            call write_log(1, bufout)
         enddo
      endif

      ! use LAPACK routine for general banded system

      if (isys.eq.1 .and. .not.has_xdata) then
         info = 0
      elseif (isys.eq.1) then
         call dgbsv(nsplx, kl, ku, nrhs, sysmat, lda, ipiv, cx, ldb, info)
      elseif (isys.eq.2) then
         call dgbsv(nsplx, kl, ku, nrhs, sysmat, lda, ipiv, cy, ldb, info)
      else
         call dgbsv(nsplx, kl, ku, nrhs, sysmat, lda, ipiv, cz, ldb, info)
      endif

      if (info.ne.0 .or. idebug.ge.6) then
         write(bufout,'(2(a,i0))') ' done dgbsv( ',isys,' ), info= ',info
         call write_log(1, bufout)
      endif

      if (idebug.ge.2) then
         if (isys.eq.1 .and. has_xdata) call write_log('B-spline coefficients cx2d:')
         if (isys.eq.2) call write_log('B-spline coefficients cy2d:')
         if (isys.eq.3) call write_log('B-spline coefficients cz2d:')

         jstep = (nsply-1) / 10 + 1
         do i = 1, nsplx
            if (i.le.20 .or. i.ge.nsplx-19 .or. idebug.ge.10) then
               if (isys.eq.1 .and. has_xdata) write(bufout,'(i5,1x, 10g14.6)') i, (cx(i,j), j=1,nsply,jstep)
               if (isys.eq.2) write(bufout,'(i5,1x, 10g14.6)') i, (cy(i,j), j=1,nsply,jstep)
               if (isys.eq.3) write(bufout,'(i5,1x, 10g14.6)') i, (cz(i,j), j=1,nsply,jstep)
               call write_log(1, bufout)
            endif
         enddo
      endif

   enddo ! isys

   return

   ! error handling for memory allocation

99 continue

   write(bufout,'(3a,i4,2a)') ' ERROR in memory allocation (',trim(namvar),':',istat,'): ',trim(errmsg)
   call write_log(1, bufout)
   call abort_run()

end subroutine bspline_solve_coef1d_intpol

!------------------------------------------------------------------------------------------------------------

module subroutine get_next_interval(n, mask, iofs, ix0, ix1)
!--function: locate a range ix0:ix1 of true-values in mask(iofs:n), if not found, set ix1 < ix0
   implicit none
!--subroutine arguments:
   integer, intent(in)    :: n, iofs
   logical, intent(in)    :: mask(n)
   integer, intent(out)   :: ix0, ix1
!--local variables

   if (n.le.0 .or. iofs.gt.n) then
      ix0 = 1
      ix1 = 0
   else

      ! find first element >= iofs with mask=T

      ix0   = max(1, iofs)
      do while(.not.mask(ix0) .and. ix0.lt.n)
         ix0 = ix0 + 1
      enddo

      if (.not.mask(ix0)) then
         ix1 = ix0 - 1
      else
         ! find first element with mask=F

         ix1 = ix0
         do while(mask(ix1) .and. ix1.lt.n)
            ix1 = ix1 + 1
         enddo

         ! if found, step back one point
         if (.not.mask(ix1)) ix1 = ix1 - 1

      endif
   endif

end subroutine get_next_interval

!------------------------------------------------------------------------------------------------------------

module subroutine bspline_make_2d_bspline(spl2d, nmeasu, nmeasv, ui, vj, x2d, y2d, z2d, has_xdata,      &
                                                                  use_approx, my_ierror, mask)
!--function: compute parametric 2D tensor interpolating B-spline for data (x2d, y2d, z2d),
!            store in B-form (cij_x, cij_y, cij_z).
   implicit none
!--subroutine arguments:
   type(t_bspline2d)                  :: spl2d
   integer,      intent(in)           :: nmeasu, nmeasv
   logical,      intent(in)           :: has_xdata, use_approx
   real(kind=8), intent(in)           :: ui(nmeasu), vj(nmeasv)
   real(kind=8), intent(in)           :: x2d(nmeasu,nmeasv), y2d(nmeasu,nmeasv), z2d(nmeasu,nmeasv)
   integer,      intent(out)          :: my_ierror
   logical,      intent(in), optional :: mask(nmeasu,nmeasv)
!--local variables:
   logical                    :: use_mask, use_insert, use_repl, lfound, include_left, include_right
   integer                    :: istat, iofs, iu0, iu1, ni, j, j0, j1, jt, jt0, jt1, nj, njt, nknot,    &
                                 nrow, ncol, sub_ierror
   real(kind=8)               :: tnew, dt_ext, dt_int
   integer,      dimension(:),   allocatable :: jseg
   real(kind=8), dimension(:),   allocatable :: tv, tvi, tvx, tu, tui, tux
   real(kind=8), dimension(:,:), allocatable :: b1, b2, b3, b4, ci_x, ci_y, ci_z, wrk_ci_x, wrk_ci_y,   &
                                                wrk_ci_z, wrk_cij_x, wrk_cij_y, wrk_cij_z
   character(len=20)          :: namvar
   character(len=256)         :: errmsg

   my_ierror = 0
   use_mask  = present(mask)

   ! check requirements on data size

   if (nmeasu.lt.4 .or. nmeasv.lt.4) then
      write(bufout,'(a,2(a,i4))') ' ERROR: 2D spline needs at least 4 slices and 4 points per slice,', &
                ' nmeasu=', nmeasu,', nmeasv=',nmeasv
      call write_log(1, bufout)
      my_ierror = 404
   endif

   ! check requirements on 'missing parts'

   if (my_ierror.eq.0 .and. use_mask) then
      call bspline_check_mask(nmeasu, nmeasv, mask, my_ierror)
   endif

   ! initialize B-spline administration for knots ti, tj, coefficients cij_x, cij_y, cij_z

   use_insert = .true.  ! true: include 2nd and n-1th meas.position in knot-vector
   use_repl = .false.   ! false: extend knot sequence, replicating 1st and last intervals k-1 times

   call bspline2d_destroy(spl2d)
   if (use_insert) then
      spl2d%nsplu  = nmeasu + 2
      spl2d%nsplv  = nmeasv + 2
   else
      spl2d%nsplu  = nmeasu
      spl2d%nsplv  = nmeasv
   endif
   spl2d%nknotu = spl2d%nsplu + 4
   spl2d%nknotv = spl2d%nsplv + 4
   spl2d%nbrku  = spl2d%nknotu - 6
   spl2d%nbrkv  = spl2d%nknotv - 6
   spl2d%has_xdata = has_xdata

   call reallocate_arr(spl2d%tui, spl2d%nknotu)
   call reallocate_arr(spl2d%tvj, spl2d%nknotv)

   call reallocate_arr(spl2d%ubrk, spl2d%nbrku)
   call reallocate_arr(spl2d%vbrk, spl2d%nbrkv)

   call reallocate_arr(spl2d%cij_x, spl2d%nsplu, spl2d%nsplv) ! un-initialized if .not.has_xdata
   call reallocate_arr(spl2d%cij_y, spl2d%nsplu, spl2d%nsplv)
   call reallocate_arr(spl2d%cij_z, spl2d%nsplu, spl2d%nsplv)
   call reallocate_arr(spl2d%mask, spl2d%nsplu, spl2d%nsplv)

   associate(nknotu => spl2d%nknotu, nsplu => spl2d%nsplu,  tu_full => spl2d%tui,                       &
             nknotv => spl2d%nknotv, nsplv => spl2d%nsplv,  tv_full => spl2d%tvj,                       &
             cij_x  => spl2d%cij_x,  cij_y => spl2d%cij_y,  cij_z => spl2d%cij_z)

   if (ldebug.ge.-1) then
      if (     use_approx) call write_log(' make_2dspline: constructing approximating spline surface')
      if (.not.use_approx) call write_log(' make_2dspline: constructing interpolating spline surface')
   endif

   if (ldebug.ge.1) then
      write(bufout,'(8(a,i4))') ' 2dspline: nmeasu,v=',nmeasu,',',nmeasv,', nsplu,v=',nsplu,',',        &
                   nsplv,', nknotu,v=',nknotu,',',nknotv,', nbrku,v=',spl2d%nbrku,',',spl2d%nbrkv
      call write_log(1, bufout)
   endif

   ! make master knot-vectors tu_full, tv_full at all measurement positions + extension at start/end

   call bspline_make_knot_vector_atmeas(nmeasu, ui, use_insert, use_repl, nknotu, tu_full, 'u')
   call bspline_make_knot_vector_atmeas(nmeasv, vj, use_insert, use_repl, nknotv, tv_full, 'v')

   ! store the breakpoints ksi, unique knots in tu_full, tv_full excluding first k-1 and last k-1

   call bspline_make_breakpoints(nknotu, tu_full, tiny_dt, spl2d%nbrku, spl2d%ubrk, ldebug)
   call bspline_make_breakpoints(nknotv, tv_full, tiny_dt, spl2d%nbrkv, spl2d%vbrk, ldebug)

   namvar = 'jseg, b1-b4 (nmeasv)'
   allocate(jseg(nmeasv), b1(nmeasv,1), b2(nmeasv,2), b3(nmeasv,3), b4(nmeasv,4),                       &
            tv(nmeasv+4), tvi(nmeasv+5), tvx(nmeasv+6),                                                 &
            tu(nmeasu+4), tui(nmeasu+5), tux(nmeasu+6), stat=istat, errmsg=errmsg)
   if (istat.ne.0) goto 99

   !--------------------------------------------------------------------------------------------------------
   ! Phase 1: build B-splines per slice in lateral direction v -> (y,z)
   !--------------------------------------------------------------------------------------------------------

   namvar = 'ci_x, ci_y, ci_z'
   allocate(ci_x(nmeasv+2,nmeasu), wrk_ci_x(nmeasv+2,nmeasu), wrk_cij_x(nmeasu+2,nmeasv+2),             &
            ci_y(nmeasv+2,nmeasu), wrk_ci_y(nmeasv+2,nmeasu), wrk_cij_y(nmeasu+2,nmeasv+2),             &
            ci_z(nmeasv+2,nmeasu), wrk_ci_z(nmeasv+2,nmeasu), wrk_cij_z(nmeasu+2,nmeasv+2),             &
            stat=istat, errmsg=errmsg)
   if (istat.ne.0) goto 99

   ! consider 'sections' of slices [iu0:iu1] with the same mask per slice

   iu1 = 0
   do while(iu1.lt.nmeasu)

      ! select slices iu0:iu1 with the same mask per slice

      if (.not.use_mask) then
         iu0 = 1
         iu1 = nmeasu
      else
         iu0 = iu1 + 1
         iu1 = iu0
         lfound = .false.
         do while(.not.lfound .and. iu1.lt.nmeasu)
            lfound = (any(mask(iu0,1:nmeasv).ne.mask(iu1+1,1:nmeasv)))
            if (.not.lfound) iu1 = iu1 + 1
         enddo
      endif
      ni = iu1 - iu0 + 1

      ! determine the active region j0:j1

      if (.not.use_mask) then
         j0 = 1
         j1 = nmeasv
      else
         iofs = 1
         call get_next_interval(nmeasv, mask(iu0,1:nmeasv), iofs, j0, j1)
         if (j0.gt.j1) then
            write(bufout,'(3(a,i4),a)') ' Error: no active points on slice',iu0,' (j=',j0,',',j1,')'
            my_ierror = 181
         endif
      endif
      nj = j1 - j0 + 1

      if (ldebug.ge.1) then
         write(bufout,'(4(a,i4),a)') ' spline2d: compute  ci_[xyz] for i =[',iu0,',',iu1,'], active j =[', &
                   j0,',',j1,']'
         call write_log(1, bufout)
      endif

      ! determine local knot vector tv, selecting appropriate elements from master tv_full
      ! master == full knot vector; local == skipping 2nd and n-1th meas positions

      tv(   1:   4) = tv_full(3+j0-3:3+j0  )
      tv(   5:nj  ) = tv_full(3+j0+2:3+j1-2)
      tv(nj+1:nj+4) = tv_full(3+j1  :3+j1+3)
      nknot = nj + 4

      if (ldebug.ge.2) then
         write(bufout,'(a,i4,a,4f8.3,8x,2f8.3/,24x,a,2f8.3,8x,4f8.3)') ' spline2d(v): nknot=',nknot,    &
                ', tv(1:6) =',(tv(j), j=1,6), 'tv(n-5:n) =',(tv(j), j=nknot-5,nknot)
         call write_log(2, bufout)
      endif

      ! evaluate B-spline basis-functions at v-measurement locations

      if (my_ierror.eq.0) then
         call bspline_eval1d(nknot, tv, tiny_dt, nj, vj(j0:j1), jseg, b1, b2, b3, b4, ldebug, sub_ierror)
         my_ierror = sub_ierror
         if (my_ierror.ne.0) call write_log(' Error after eval1d(tv)')
      endif

      ! solve the linear systems B [ci_x, ci_y, ci_z] = [x2d, y2d, z2d]

      ! copy input (x2d,y2d,z2d) to work-arrays (ci_x,ci_y,ci_z), solve_coef1d_intpol overwrites rhs

      ! call print_2d_real(ni, nj, y2d(iu0:iu1,j0:j1), 'y2d', 10, 'g12.4')

      if (has_xdata) then
         wrk_ci_x(1:nj,1:ni) = transpose(x2d(iu0:iu1,j0:j1))
      endif
      wrk_ci_y(1:nj,1:ni) = transpose(y2d(iu0:iu1,j0:j1))
      wrk_ci_z(1:nj,1:ni) = transpose(z2d(iu0:iu1,j0:j1))

      ! call print_2d_real(nj, ni, wrk_ci_y, 'inp_ci_y', 10, 'g12.6')

      if (my_ierror.eq.0) then
         ! call write_log(' bspline_solve_coef1d_intpol(1)...')
         nrow = size(wrk_ci_y,1)
         call bspline_solve_coef1d_intpol(nj, jseg, b4, nrow, ni, wrk_ci_x, wrk_ci_y, wrk_ci_z,         &
                        has_xdata, ldebug, sub_ierror)
         my_ierror = sub_ierror
         if (my_ierror.ne.0) call write_log(' Error after solve_coef ci_xyz')
         if (ldebug.ge.-1) then
            if (has_xdata) call check_nan(wrk_ci_x, 'ci_x', 1)
            call check_nan(wrk_ci_y, 'ci_y', 1)
            call check_nan(wrk_ci_z, 'ci_z', 1)
         endif
         ! call print_2d_real(nj, ni, wrk_ci_y, 'out_ci_y', 10, 'g12.4')
      endif

      ! insert knots at not-a-knot positions, compute updated coefficient values

      tnew = tv_full(3+j0+1)
      nrow = nmeasv+2
      ncol = ni
      if (has_xdata) call bspline_insert_knot(nknot, tv, tnew, tvi, nrow, ncol, wrk_ci_x)
      call bspline_insert_knot(nknot, tv, tnew, tvi, nrow, ncol, wrk_ci_y)
      call bspline_insert_knot(nknot, tv, tnew, tvi, nrow, ncol, wrk_ci_z)

      tnew = tv_full(3+j1-1)
      nrow = nmeasv+2
      ncol = ni
      if (has_xdata) call bspline_insert_knot(nknot+1, tvi, tnew, tvx, nrow, ncol, wrk_ci_x)
      call bspline_insert_knot(nknot+1, tvi, tnew, tvx, nrow, ncol, wrk_ci_y)
      call bspline_insert_knot(nknot+1, tvi, tnew, tvx, nrow, ncol, wrk_ci_z)

      ! copy coefficients into full arrays

      if (has_xdata) ci_x(j0:j1+2, iu0:iu1) = wrk_ci_x(1:nj+2,1:ni)
      ci_y(j0:j1+2, iu0:iu1) = wrk_ci_y(1:nj+2,1:ni)
      ci_z(j0:j1+2, iu0:iu1) = wrk_ci_z(1:nj+2,1:ni)
      ! call print_2d_real(nmeasv+2, nmeasu, ci_z, 'ci_z', 10, 'g12.4')

   enddo ! while(iu1<nmeasu)

   deallocate(jseg, b1, b2, b3, b4, tv, tvi, tvx, wrk_ci_x, wrk_ci_y, wrk_ci_z)

   !--------------------------------------------------------------------------------------------------------
   ! Phase 2: build B-splines for longitudinal direction u -> (cx,cy,cz)
   !--------------------------------------------------------------------------------------------------------

   namvar = 'jseg, b1-b4 (nmeasu)'
   allocate(jseg(nmeasu), b1(nmeasu,1), b2(nmeasu,2), b3(nmeasu,3), b4(nmeasu,4), stat=istat,           &
            errmsg=errmsg)
   if (istat.ne.0) goto 99

   if (has_xdata) cij_x(1:nsplu,1:nsplv) = -999d0
   cij_y(1:nsplu,1:nsplv) = -999d0
   cij_z(1:nsplu,1:nsplv) = -999d0
   spl2d%mask(1:nsplu,1:nsplv) = 0

   ! consider 'regions' of interpolation paths [j0:j1] with the same mask per path

   j1 = 0
   do while(j1.lt.nmeasv)

      ! select interpolation paths j0:j1 with the same mask per path

      if (.not.use_mask) then
         j0 = 1
         j1 = nmeasv
      else
         j0 = j1 + 1
         j1 = j0
         lfound = .false.
         do while(.not.lfound .and. j1.lt.nmeasv)
            lfound = (any(mask(1:nmeasu,j0).ne.mask(1:nmeasu,j1+1)))
            ! write(bufout,'(2(a,i4),a,l2)') 'j0=',j0,', j1=',j1,', check j1+1: found=',lfound
            ! call write_log(1, bufout)
            if (.not.lfound) j1 = j1 + 1
         enddo
      endif
      nj = j1 - j0 + 1

      ! consider 'sections' of active points iu0:iu1 within path j0
      ! Note: for one j, mask may have multiple 'sections'   0 0 0 1 1 1 1 0 0 0 0 1 1 1 1 0 0 

      iu1 = 0
      do while(iu1.lt.nmeasu)

         if (.not.use_mask) then
            iu0 = 1
            iu1 = nmeasu
         else
            iofs = iu1 + 1
            call get_next_interval(nmeasu, mask(1:nmeasu,j0), iofs, iu0, iu1)
         endif
         ni = iu1 - iu0 + 1

         if (ni.le.0) then
            iu1 = nmeasu
         else

            ! determine whether current slices are longer or shorter than previous iu0-1 and following iu1+1
            ! longer slices: include one column from neighbour; shorter slices: give one column away

            include_left = .true.
            if (j0.gt.1 .and. iu0.gt.1) then
               include_left = include_left .and. .not.mask(iu0-1,j0-1)
            endif
            if (j0.gt.1 .and. iu1.lt.nmeasu) then
               include_left = include_left .and. .not.mask(iu1+1,j0-1)
            endif
   
            include_right = .true.
            if (j1.lt.nmeasv .and. iu0.gt.1) then
               include_right = include_right .and. .not.mask(iu0-1,j1+1)
            endif
            if (j1.lt.nmeasv .and. iu1.lt.nmeasu) then
               include_right = include_right .and. .not.mask(iu1+1,j1+1)
            endif

            if (include_left) then
               jt0 = j0     ! add guard band
            else
               jt0 = j0 + 2 ! move 1st column to neighbour
            endif
            if (include_right) then
               jt1 = j1 + 2 ! add guard band
            else
               jt1 = j1     ! move 1st column to neighbour
            endif
            njt = jt1 - jt0 + 1

            if (ldebug.ge.1) then
               write(bufout,'(4(a,i4),a)') ' spline2d: compute cij_[xyz] for jt=[',jt0,',',jt1,         &
                           '], active i =[', iu0,',',iu1,']'
               call write_log(1, bufout)
            endif

            if (use_approx) then

               ! approximating spline -- interpret (cx, cy, cz) as control points for 2D spline

               ! copy coefficients into full arrays

               if (has_xdata) cij_x(iu0+1:iu1+1,jt0:jt1) = transpose(ci_x(jt0:jt1,iu0:iu1))
               cij_y(iu0+1:iu1+1,jt0:jt1) = transpose(ci_y(jt0:jt1,iu0:iu1))
               cij_z(iu0+1:iu1+1,jt0:jt1) = transpose(ci_z(jt0:jt1,iu0:iu1))

               ! add fantom control points at row iu0, extrapolating iu0+2 -> iu0+1 -> iu0

               dt_ext = tu_full(3+iu0)   - tu_full(3+iu0-1)
               dt_int = tu_full(3+iu0+1) - tu_full(3+iu0)
               do jt = jt0, jt1
                  if (has_xdata) then
                     cij_x(iu0,jt) = (1d0+dt_ext/dt_int) * cij_x(iu0+1,jt) - dt_ext/dt_int * cij_x(iu0+2,jt)
                  endif
                  cij_y(iu0,jt) = (1d0+dt_ext/dt_int) * cij_y(iu0+1,jt) - dt_ext/dt_int * cij_y(iu0+2,jt)
                  cij_z(iu0,jt) = (1d0+dt_ext/dt_int) * cij_z(iu0+1,jt) - dt_ext/dt_int * cij_z(iu0+2,jt)
               enddo

               ! add fantom control points at row iu1+2, extrapolating iu1 -> iu1+1 -> iu1+2

               dt_ext = tu_full(3+iu1+2) - tu_full(3+iu1+1)
               dt_int = tu_full(3+iu1+1) - tu_full(3+iu1)
               do jt = jt0, jt1
                  if (has_xdata) then
                     cij_x(iu1+2,jt) = (1d0+dt_ext/dt_int) * cij_x(iu0+1,jt) - dt_ext/dt_int * cij_x(iu1,jt)
                  endif
                  cij_y(iu1+2,jt) = (1d0+dt_ext/dt_int) * cij_y(iu0+1,jt) - dt_ext/dt_int * cij_y(iu1,jt)
                  cij_z(iu1+2,jt) = (1d0+dt_ext/dt_int) * cij_z(iu0+1,jt) - dt_ext/dt_int * cij_z(iu1,jt)
               enddo

               spl2d%mask(iu0:iu1+2,jt0:jt1)  = 1

            else

               ! interpolating spline -- interpret (cx, cy, cz) as data points for 2D spline

               ! determine local knot vector tu, selecting appropriate elements from master tu_full
               ! master == full knot vector; local == skipping 2nd and n-1th meas positions

               tu(   1:   4) = tu_full(3+iu0-3:3+iu0  )
               tu(   5:ni  ) = tu_full(3+iu0+2:3+iu1-2)
               tu(ni+1:ni+4) = tu_full(3+iu1  :3+iu1+3)
               nknot = ni + 4

               if (ldebug.ge.2) then
                  write(bufout,'(a,i4,a,4g10.2,8x,2g10.2/,24x,a,2g10.2,8x,4g10.2)') ' spline2d(u): nknot=', &
                         nknot, ', tu(1:6) =',(tu(j), j=1,6), 'tu(n-5:n) =',(tu(j), j=nknot-5,nknot)
                  call write_log(2, bufout)
               endif
   
               ! evaluate B-spline basis-functions at u-measurement locations

               if (my_ierror.eq.0) then
                  call bspline_eval1d(nknot, tu, tiny_dt, ni, ui(iu0:iu1), jseg, b1, b2, b3, b4,        &
                              ldebug, sub_ierror)
                  my_ierror = sub_ierror
                  if (my_ierror.ne.0) call write_log(' Error after eval1d(tu)')
               endif

               ! copy input (ci_[xyz]) to work-arrays (cij_[xyz]), solve_coef1d_intpol overwrites rhs

               if (has_xdata) wrk_cij_x(1:ni,1:njt) = transpose(ci_x(jt0:jt1,iu0:iu1))
               wrk_cij_y(1:ni,1:njt) = transpose(ci_y(jt0:jt1,iu0:iu1))
               wrk_cij_z(1:ni,1:njt) = transpose(ci_z(jt0:jt1,iu0:iu1))

               ! solve the linear systems B [cij_x, cij_y, cij_z] = [ci_x, ci_y, ci_z]

               if (my_ierror.eq.0) then
                  ! call write_log(' bspline_solve_coef1d_intpol(1)...')
                  nrow = size(wrk_cij_y,1)
                  call bspline_solve_coef1d_intpol(ni, jseg, b4, nrow, njt, wrk_cij_x, wrk_cij_y,       &
                              wrk_cij_z, has_xdata, ldebug, sub_ierror)
                  my_ierror = sub_ierror
                  if (my_ierror.ne.0) call write_log(' Error after solve_coef cij_xyz')
                  if (ldebug.ge.-1) then
                     if (has_xdata) call check_nan(wrk_cij_x, 'cij_x', 1)
                     call check_nan(wrk_cij_y, 'cij_y', 1)
                     call check_nan(wrk_cij_z, 'cij_z', 1)
                  endif
               endif

               ! insert knots at not-a-knot positions, compute updated coefficient values

               tnew = tu_full(3+iu0+1)
               nrow = nmeasu+2
               ncol = njt
               if (has_xdata) call bspline_insert_knot(nknot, tu, tnew, tui, nrow, ncol, wrk_cij_x)
               call bspline_insert_knot(nknot, tu, tnew, tui, nrow, ncol, wrk_cij_y)
               call bspline_insert_knot(nknot, tu, tnew, tui, nrow, ncol, wrk_cij_z)

               tnew = tu_full(3+iu1-1)
               nrow = nmeasu+2
               ncol = njt
               if (has_xdata) call bspline_insert_knot(nknot+1, tui, tnew, tux, nrow, ncol, wrk_cij_x)
               call bspline_insert_knot(nknot+1, tui, tnew, tux, nrow, ncol, wrk_cij_y)
               call bspline_insert_knot(nknot+1, tui, tnew, tux, nrow, ncol, wrk_cij_z)

               ! copy coefficients into full arrays

               if (has_xdata) cij_x(iu0:iu1+2,jt0:jt1) = wrk_cij_x(1:ni+2,1:njt)
               cij_y(iu0:iu1+2,jt0:jt1) = wrk_cij_y(1:ni+2,1:njt)
               cij_z(iu0:iu1+2,jt0:jt1) = wrk_cij_z(1:ni+2,1:njt)
               spl2d%mask(iu0:iu1+2,jt0:jt1)  = 1

            endif ! use_approx
         endif ! iu1>=iu0
      enddo ! while more sections (iu1<nmeasu)
   enddo ! while more intpol.paths (j1<nmeasv)

   deallocate(jseg, b1, b2, b3, b4, ci_x, ci_y, ci_z, wrk_cij_x, wrk_cij_y, wrk_cij_z)
   end associate

   return

   ! error handling for memory allocation

99 continue

   write(bufout,'(3a,i4,2a)') ' ERROR in memory allocation (',trim(namvar),':',istat,'): ',trim(errmsg)
   call write_log(1, bufout)
   call abort_run()

end subroutine bspline_make_2d_bspline

!------------------------------------------------------------------------------------------------------------

end submodule m_bspline_make
