!------------------------------------------------------------------------------------------------------------
! m_spline_make - construction of 1D spline in PP-form (piecewise polynomial)
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_spline_make
   use m_globals
   use m_markers
   use m_ptrarray
   use m_interp_1d
   use m_spline_def
   implicit none
   private

   ! Debugging for module m_spline_make

   public  splinemake_set_debug

   integer  :: ldebug    =  0    ! local level of debugging
   integer  :: ii_debug  = -1    ! output point for which detailed info is requested (-1 = none)
   integer  :: iel_debug = -1    ! input element for which detailed info is requested (-1 = none)

   ! Additional functions defined on splines:

   public  spline_check_updates

   public  spline_add_topview

   private ppspline_make_simple_sec_intpol
   private ppspline_make_simple_sec_smoothing
   public  ppspline_make_simple_kink

   public  ppspline_make_spline
   public  ppspline_make_spline_kink
   public  ppspline_make_spline_nokink

   interface ppspline_make_spline
      module procedure ppspline_make_spline_kink
      module procedure ppspline_make_spline_nokink
   end interface ppspline_make_spline

contains

!------------------------------------------------------------------------------------------------------------

subroutine splinemake_set_debug(new_ldebug, new_ii_debug, new_iel_debug)
!--function: enable/disable debug output of spline routines
   implicit none
!--subroutine arguments:
   integer, intent(in)           :: new_ldebug       ! level of debug output required
   integer, intent(in), optional :: new_ii_debug     ! specific point of interest for debugging
   integer, intent(in), optional :: new_iel_debug    ! specific point of interest for debugging

   ldebug = new_ldebug

   if (present(new_ii_debug)) then
      ii_debug = new_ii_debug
   endif
   if (present(new_iel_debug)) then
      iel_debug = new_iel_debug
   endif

   if (ldebug.ge.3) then
      write(bufout,'(a,i3,2(a,i7))') ' spline-make: debugging level =',ldebug,', ii_debug =', ii_debug, &
                ', iel_debug =', iel_debug
      call write_log(1, bufout)
   endif

end subroutine splinemake_set_debug

!------------------------------------------------------------------------------------------------------------

subroutine spline_add_breaks(spl, nbrk, sbrk, ipbrk, my_ierror)
!--function: insert additional breaks at sbrk in PP-form of spline, except where sbrk coincides with
!            existing breaks
   implicit none
!--subroutine arguments:
   type(t_spline)                     :: spl
   integer,      intent(in)           :: nbrk
   real(kind=8), intent(in)           :: sbrk(nbrk)
   integer,      intent(out)          :: ipbrk(nbrk)
   integer,      intent(out)          :: my_ierror
!--local variables:
   real(kind=8), parameter            :: tiny_s = 1d-10
   logical      :: use_orig
   integer      :: nnew, ib, iorig, ibrk, inew
   real(kind=8) :: dist(2)
   real(kind=8), dimension(:,:), allocatable :: axbrk, aybrk, azbrk
   real(kind=8), dimension(:),   pointer     :: snew
   real(kind=8), dimension(:,:), pointer     :: axnew, aynew, aznew

   my_ierror = 0

   ! check that all sbrk lie in [s(1),s(end)]

   do ib = 1, nbrk
      if (sbrk(ib).lt.spl%s(1) .or. sbrk(ib).gt.spl%s(spl%npnt)) then
         write(bufout,'(a,i3,3(a,f12.4),a)') ' INTERNAL ERROR(add_breaks): sb(',ib,')=',sbrk(ib),       &
                ' out of range, s=[',spl%s(1),',',spl%s(spl%npnt),'].'
         call write_log(1, bufout)
         write(bufout,*) spl%s(1), sbrk(ib), spl%s(spl%npnt)
         call write_log(1, bufout)
         my_ierror = -1
      endif
   enddo

   ! check that all sbrk are in strictly increasing order

   do ib = 2, nbrk
      if (sbrk(ib).le.sbrk(ib-1)) then
         write(bufout,'(2(a,i3,a,f12.4),a)') ' INTERNAL ERROR(add_breaks): sb(',ib,')=',sbrk(ib),       &
                ' must be > sb(',ib-1,')=',sbrk(ib-1),'.'
         call write_log(1, bufout)
         my_ierror = -2
      endif
   enddo

   ! evaluate spline and derivatives at new break points

   allocate(axbrk(nbrk,4), aybrk(nbrk,4), azbrk(nbrk,4))

   if (my_ierror.eq.0) then
      if (spl%has_xdata) call spline_eval(spl, ikXDIR, nbrk, sbrk, my_ierror, 999d0, axbrk(:,1),        &
                                                                    axbrk(:,2), axbrk(:,3), axbrk(:,4))
      call spline_eval(spl, ikYDIR, nbrk, sbrk, my_ierror, 999d0, aybrk(:,1), aybrk(:,2), aybrk(:,3),   &
                                                                                            aybrk(:,4))
      call spline_eval(spl, ikZDIR, nbrk, sbrk, my_ierror, 999d0, azbrk(:,1), azbrk(:,2), azbrk(:,3),   &
                                                                                            azbrk(:,4))

      ! in breakpoints, second derivative f'' = 2*a2, third derivative f''' = 6*a3

      aybrk(:,3) = aybrk(:,3) / 2d0
      aybrk(:,4) = aybrk(:,4) / 6d0
      azbrk(:,3) = azbrk(:,3) / 2d0
      azbrk(:,4) = azbrk(:,4) / 6d0

      if (ldebug.ge.4) then
         do ib = 1, nbrk
            write(bufout,'(i4,3(a,f14.6))') ib,': s =',sbrk(ib),', y =', aybrk(ib,1), ', z =', azbrk(ib,1)
            call write_log(1, bufout)
         enddo
      endif

   endif

   ! merge original and new s-positions and corresponding data

   if (my_ierror.eq.0) then

      ! Note: some new breaks may be skipped, if coinciding with existing break points.
      !       actual nnew after merging may be lower than npnt + nbrk, some positions in arrays unused

      nnew  = spl%npnt + nbrk
      allocate(snew(nnew), aynew(nnew,4), aznew(nnew,4))
      if (spl%has_xdata) allocate(axnew(nnew,4))

      iorig = 1
      ibrk  = 1
      inew  = 1
      do while(iorig.le.spl%npnt .or. ibrk.le.nbrk)

         ! select point from original or new breakpoint data

         if (ibrk.gt.nbrk) then
            use_orig = .true.
         elseif (iorig.gt.spl%npnt) then
            use_orig = .false.
         else
            use_orig = (spl%s(iorig).le.sbrk(ibrk))
         endif

         if (use_orig) then

            ! copy point from original data

            snew(inew)      = spl%s(iorig)
            if (spl%has_xdata) axnew(inew,1:4) = spl%axspl(iorig,1:4)
            aynew(inew,1:4) = spl%ayspl(iorig,1:4)
            aznew(inew,1:4) = spl%azspl(iorig,1:4)

            ! set next point to use in output arrays

            inew = inew + 1

            ! set next point to use from original data

            iorig = iorig + 1

         else

            ! process point from new breakpoint data

            ! determine distance of s(ibrk) to nearby points s(iorig)

            dist = (/ 999d0, 999d0 /)
            if (inew.gt.1)         dist(1) = abs(sbrk(ibrk)-snew(inew-1))
            if (iorig.le.spl%npnt) dist(2) = abs(sbrk(ibrk)-spl%s(iorig))

            ! skip breaks that coincide with end-points of pre-existing segments

            if (dist(1).le.tiny_s) then

               ! record position of ibrk in output spline
 
               ipbrk(ibrk) = inew - 1

               if (ldebug.ge.2) then
                  write(bufout, '(2(a,i4))') ' add_breaks: ib=',ibrk,' coincides with ipnt=',inew-1
                  call write_log(1, bufout)
               endif

            elseif (dist(2).le.tiny_s) then

               ! record position of ibrk in output spline
 
               ipbrk(ibrk) = inew

               if (ldebug.ge.2) then
                  write(bufout, '(2(a,i4))') ' add_breaks: ib=',ibrk,' coincides with ipnt=',inew
                  call write_log(1, bufout)
               endif

            else

               ! record position of ibrk in output spline
 
               ipbrk(ibrk) = inew

               if (ldebug.ge.2) then
                  write(bufout, '(2(a,i4))') ' add_breaks: insert ib=',ibrk,' at ipnt=',inew
                  call write_log(1, bufout)
               endif

               snew(inew)      = sbrk(ibrk)
               if (spl%has_xdata) axnew(inew,1:4) = axbrk(ibrk,1:4)
               aynew(inew,1:4) = aybrk(ibrk,1:4)
               aznew(inew,1:4) = azbrk(ibrk,1:4)

               ! set next point to use in output arrays

               inew = inew + 1

            endif

            ! set next point to use from new breakpoint data

            ibrk  = ibrk + 1
         endif

      enddo

      ! deallocate existing arrays from spl-structure and replace pointers with new ones

      if (spl%has_xdata) deallocate(spl%axspl)
      deallocate(spl%s, spl%ayspl, spl%azspl)

      spl%s     => snew

      if (spl%has_xdata) then
         spl%axspl => axnew
         spl%ax0   => spl%axspl(:,1)
         spl%ax1   => spl%axspl(:,2)
         spl%ax2   => spl%axspl(:,3)
         spl%ax3   => spl%axspl(:,4)
      endif

      spl%ayspl => aynew
      spl%ay0   => spl%ayspl(:,1)
      spl%ay1   => spl%ayspl(:,2)
      spl%ay2   => spl%ayspl(:,3)
      spl%ay3   => spl%ayspl(:,4)

      spl%azspl => aznew
      spl%az0   => spl%azspl(:,1)
      spl%az1   => spl%azspl(:,2)
      spl%az2   => spl%azspl(:,3)
      spl%az3   => spl%azspl(:,4)

      spl%npnt  =  inew - 1

   endif

   deallocate(axbrk, aybrk, azbrk)

   ! clear information on uni-valued sections

   spl%nsec_uniy = 0
   if (associated(spl%ipnt_uniy)) deallocate(spl%ipnt_uniy)
   spl%ipnt_uniy => NULL()

end subroutine spline_add_breaks

!------------------------------------------------------------------------------------------------------------

subroutine spline_set_topview_sections(spl, view_minz, max_ntop, nsec_top, ysec_top, iuni_top, my_ierror)
!--function: determine the sections [ y_j, y_{j+1} ] of the top (rail) or bottom view (wheel) and
!            set pointers to corresponding uni-valued sections in spline
   implicit none
!--subroutine arguments:
   type(t_spline)                     :: spl
   logical,      intent(in)           :: view_minz   ! view surface from z=-inf (rail) or z=+inf (wheel)
   integer,      intent(in)           :: max_ntop
   integer,      intent(out)          :: nsec_top, my_ierror
   real(kind=8), intent(out)          :: ysec_top(max_ntop+1)
   integer,      intent(out)          :: iuni_top(max_ntop)
!--local variables:
   integer    :: iuni, itop, juni, ksec, ip0, ip1
   logical    :: ldone
   logical,      dimension(:), allocatable :: has_minz
   real(kind=8), dimension(:), allocatable :: sec_ymin, sec_ymax

   my_ierror = 0

   associate( npnt => spl%npnt, nuni => spl%nsec_uniy)
   allocate(sec_ymin(nuni), sec_ymax(nuni), has_minz(nuni))

   if (ldebug.ge.2) call write_log(' spline: set_topview_sections...')

   ! initialize: each section entirely visible

   do iuni = 1, nuni
      ip0  = spl%ipnt_uniy(iuni)
      ip1  = spl%ipnt_uniy(iuni+1)
      sec_ymin(iuni) = min( spl%ay0(ip0), spl%ay0(ip1) )
      sec_ymax(iuni) = max( spl%ay0(ip0), spl%ay0(ip1) )
   enddo

   if (ldebug.ge.3) then
      call write_log(' initial ranges:')
      do iuni = 1, nuni
         write(bufout,'(a,i3,2(a,f12.6),a)') '    iuni',iuni,': ymin,max = [',sec_ymin(iuni),',',       &
                sec_ymax(iuni),']'
         call write_log(1, bufout)
      enddo
   endif

   ! determine which sections lie lower (rail) or higher (wheel) than both neighbours
   !  -- assumed enitrely visible in end-result

   do iuni = 1, nuni
      ip0  = spl%ipnt_uniy(iuni)
      ip1  = spl%ipnt_uniy(iuni+1)

      if (view_minz) then       ! surface viewed from z = -inf (rail)
         has_minz(iuni) = (ip0.eq.1 .or. spl%az1(ip0).lt.0d0) .and. (ip1.eq.npnt .or. spl%az1(ip1).gt.0d0)
      else                      ! surface viewed from z = +inf (wheel)
         has_minz(iuni) = (ip0.eq.1 .or. spl%az1(ip0).gt.0d0) .and. (ip1.eq.npnt .or. spl%az1(ip1).lt.0d0)
      endif

      if (ldebug.ge.2 .and. has_minz(iuni)) then
         write(bufout,'(a,i3,4(a,f9.3),a)') '    iuni',iuni,' is fully visible: slopes dz/ds=',        &
             spl%az1(ip0),',', spl%az1(ip1),' at y=[',spl%ay0(ip0),',',spl%ay0(ip1),']'
         call write_log(1, bufout)
      endif
   enddo

   ! for each section that has_minz(iuni): flood to the left and flood to the right

   do iuni = 1, nuni
      if (has_minz(iuni)) then

         ! push coordinates to lower juni until another section is found with hard [ymin,ymax]

         juni  = iuni
         ldone = (juni.le.1)

         do while(.not.ldone)

            juni  = juni - 1
            if (juni.le.1) ldone = .true.       ! done after this one

            if (view_minz) then     ! surface viewed from z = -inf (rail); lower juni to the left of iuni
               if (.not.has_minz(juni)) then
                  ip0 = spl%ipnt_uniy(juni)
                  ip1 = spl%ipnt_uniy(juni+1)
                  sec_ymin(juni) = min(sec_ymin(juni), sec_ymin(juni+1))
                  sec_ymax(juni) = min(sec_ymax(juni), sec_ymin(juni+1))
               else
                  ldone = .true.
               endif
            else                    ! surface viewed from z = +inf (wheel); lower juni to the right of iuni
               if (.not.has_minz(juni)) then
                  ip0 = spl%ipnt_uniy(juni)
                  ip1 = spl%ipnt_uniy(juni+1)
                  sec_ymin(juni) = max(sec_ymin(juni), sec_ymax(juni+1))
                  sec_ymax(juni) = max(sec_ymax(juni), sec_ymax(juni+1))
               else
                  ldone = .true.
               endif
            endif

         enddo

         ! push coordinates to higher juni until another section is found with hard [ymin,ymax]

         juni  = iuni
         ldone = (juni.ge.nuni)

         do while(.not.ldone)

            juni  = juni + 1
            if (juni.ge.nuni) ldone = .true.       ! done after this one

            if (view_minz) then     ! surface viewed from z = -inf (rail); higher juni to the right of iuni
               if (.not.has_minz(juni)) then
                  ip0 = spl%ipnt_uniy(juni)
                  ip1 = spl%ipnt_uniy(juni+1)
                  sec_ymin(juni) = max(sec_ymin(juni), sec_ymax(juni-1))
                  sec_ymax(juni) = max(sec_ymax(juni), sec_ymax(juni-1))
               else
                  ldone = .true.
               endif
            else                    ! surface viewed from z = +inf (wheel); higher juni to the left of iuni
               if (.not.has_minz(juni)) then
                  ip0 = spl%ipnt_uniy(juni)
                  ip1 = spl%ipnt_uniy(juni+1)
                  sec_ymin(juni) = min(sec_ymin(juni), sec_ymin(juni-1))
                  sec_ymax(juni) = min(sec_ymax(juni), sec_ymin(juni-1))
               else
                  ldone = .true.
               endif
            endif
         enddo

         ! print information visibility of each section
         ! ldebug=1: final, visible only; 2: final, all; 3: after each flooding, all

      endif ! has_minz

      if (ldebug.ge.3) then
         write(bufout,'(a,i3,a)') ' after processing iuni=',iuni,':'
         call write_log(1, bufout)
      endif

      do ksec = 1, nuni
         if (sec_ymax(ksec).gt.sec_ymin(ksec)) then
            write(bufout,'(a,i3,2(a,f14.6),a)') '    iuni',ksec,': y = [',sec_ymin(ksec),',',           &
                           sec_ymax(ksec),']: visible'
            if ((ldebug.ge.2 .and. iuni.ge.nuni) .or. (ldebug.ge.3 .and. has_minz(iuni)))               &
               call write_log(1, bufout)
         else
            write(bufout,'(a,i3,2(a,f14.6),a)') '    iuni',ksec,': y = [',sec_ymin(ksec),',',           &
                           sec_ymax(ksec),']: not visible'
            if ((ldebug.ge.2 .and. iuni.ge.nuni) .or. (ldebug.ge.3 .and. has_minz(iuni)))               &
               call write_log(1, bufout)
         endif
      enddo

   enddo ! for iuni

   ! determine list of visible parts in top-view

   nsec_top = 0
   do iuni = 1, nuni
      if (sec_ymax(iuni).gt.sec_ymin(iuni)) then
         nsec_top = nsec_top + 1
         if (view_minz) then    ! rail: store y-values left-to-right
            ysec_top(nsec_top)   = sec_ymin(iuni)
            ysec_top(nsec_top+1) = sec_ymax(iuni)
         else                   ! wheel: store y-values right-to-left
            ysec_top(nsec_top)   = sec_ymax(iuni)
            ysec_top(nsec_top+1) = sec_ymin(iuni)
         endif
         iuni_top(nsec_top)   = iuni
      endif
   enddo

   ! print information on the visible parts

   if (ldebug.ge.2) then
      write(bufout,'(a,i4,a)') ' top view consists of',nsec_top,' parts [y_j, y_{j+1}]'
      call write_log(1, bufout)

      do itop = 1, nsec_top
         write(bufout,'(a,i4,2(a,f10.3),a,i4)') ' itop=',itop,': y =[', ysec_top(itop),',',            &
                        ysec_top(itop+1),'], iuni=', iuni_top(itop)
         call write_log(1, bufout)
      enddo
   endif

   deallocate(sec_ymin, sec_ymax, has_minz)
   end associate

end subroutine spline_set_topview_sections

!------------------------------------------------------------------------------------------------------------

subroutine spline_add_topview(spl, view_minz, my_ierror)
!--function: determine sections in spline with uni-valued y(s) and table { ybrk } with pointers to
!            visible parts of the surface
   implicit none
!--subroutine arguments:
   type(t_spline)                     :: spl
   logical,      intent(in)           :: view_minz   ! view surface from z=-inf (rail) or z=+inf (wheel)
   integer,      intent(out)          :: my_ierror
!--local variables:
   integer                                 :: isec, ix, nxtrm
   integer,      dimension(:), allocatable :: ip_xtrm
   real(kind=8), dimension(:), allocatable :: s_xtrm

   my_ierror = 0

   ! get the s-positions of extremal y-values, local min/max of y = y(s) with y'(s) = 0

   allocate(s_xtrm(spl%npnt))

   ! call spline_print(spl, 'rail', 5)

   ! call write_log(' locate_extremal...')
   call locate_extremal_values(spl%npnt, spl%s, spl%ay3, spl%ay2, spl%ay1, spl%ay0, nxtrm, s_xtrm, my_ierror)

   if (ldebug.ge.1) then
      write(bufout,'(a,i3,a)') ' spline has',nxtrm,' locally extremal y-values'
      call write_log(1, bufout)
   endif

   ! insert additional break points into PP-spline at extremal y-values (if any)

   if (my_ierror.eq.0 .and. nxtrm.gt.0) then
      allocate(ip_xtrm(nxtrm))
      call spline_add_breaks(spl, nxtrm, s_xtrm, ip_xtrm, my_ierror)
   endif

   if (ldebug.ge.1 .and. my_ierror.eq.0) then
      do ix = 1, nxtrm
         write(bufout,'(i4,a,f14.6,a,i5,2(a,f14.6))') ix,': s =',s_xtrm(ix),', ip =',ip_xtrm(ix),       &
                ', y =', spl%ay0(ip_xtrm(ix)),', z =', spl%az0(ip_xtrm(ix))
         call write_log(1, bufout)
      enddo
   endif

   if (my_ierror.eq.0) then

      ! add description of spline sections [y_i, y_{i+1}] between extremal values, 
      ! monotonously increasing or decreasing y(s), invertible with uni-valued s = s(y)

      spl%nsec_uniy = nxtrm + 1
      call reallocate_arr(spl%ipnt_uniy, spl%nsec_uniy+1)

      isec = 1
      spl%ipnt_uniy(isec) = 1

      do ix = 1, nxtrm
         if (ip_xtrm(ix).gt.spl%ipnt_uniy(isec)) then
            isec = isec + 1
            spl%ipnt_uniy(isec) = ip_xtrm(ix)
         endif
      enddo

      if (nxtrm.le.0) then
         isec = isec + 1
         spl%ipnt_uniy(isec) = spl%npnt
      elseif (spl%npnt.gt.ip_xtrm(nxtrm)) then
         isec = isec + 1
         spl%ipnt_uniy(isec) = spl%npnt
      endif

      spl%nsec_uniy = isec - 1

      ! determine parts per section visible in top view

      call reallocate_arr(spl%ysec_top, nxtrm+2)
      call reallocate_arr(spl%iuni_top, nxtrm+1)

      call spline_set_topview_sections(spl, view_minz, nxtrm+1, spl%nsec_top, spl%ysec_top,             &
                spl%iuni_top, my_ierror)
   endif
   
   if (nxtrm.gt.0) deallocate(ip_xtrm)
   deallocate(s_xtrm)

end subroutine spline_add_topview

!------------------------------------------------------------------------------------------------------------

subroutine spline_check_updates(spl, nprf, s_prf, x_prf, y_prf, z_prf, k_chk, dist_max, my_ierror)
!--function: determine max distance between parametric spline (x(s),y(s),z(s)) and original data
   implicit none
!--subroutine arguments:
   type(t_spline)                     :: spl
   integer,      intent(in)           :: nprf
   real(kind=8), intent(in)           :: s_prf(nprf), x_prf(nprf), y_prf(nprf), z_prf(nprf)
   integer,      intent(in)           :: k_chk      ! refinement factor for checking
   real(kind=8), intent(out)          :: dist_max
   integer,      intent(out)          :: my_ierror
!--local variables:
   integer                    :: nchk, ichk, iseg, k, imax, sub_ierror
   real(kind=8)               :: ds, dst2, dst2mx
   real(kind=8), allocatable  :: schk(:), xspl(:), yspl(:), zspl(:), xint(:), yint(:), zint(:)

   ! allocate temporary arrays for checking spline vs linear interpolant

   nchk = 1 + k_chk * (nprf-1)
   allocate(schk(nchk), xspl(nchk), yspl(nchk), zspl(nchk), xint(nchk), yint(nchk), zint(nchk))

   ! fill s-values used for checking

   ichk = 1
   schk(ichk) = s_prf(ichk)

   do iseg = 1, nprf-1
      ds = s_prf(iseg+1) - s_prf(iseg)
      do k = 1, k_chk
         ichk = ichk + 1
         schk(ichk) = s_prf(iseg) + real(k) / real(k_chk) * ds
      enddo
   enddo

   ! linear interpolation of grid x,y,z at schk-values

   if (spl%has_xdata) then
      call interp_1d( nprf, s_prf, x_prf, nchk, schk, xint, sub_ierror)
      if (my_ierror.eq.0) my_ierror = sub_ierror
   endif

   call interp_1d( nprf, s_prf, y_prf, nchk, schk, yint, sub_ierror)
   if (my_ierror.eq.0) my_ierror = sub_ierror

   call interp_1d( nprf, s_prf, z_prf, nchk, schk, zint, sub_ierror)
   if (my_ierror.eq.0) my_ierror = sub_ierror

   ! evaluate spline x,y,z at schk-values

   if (spl%has_xdata) then
      call spline_eval(spl, ikXDIR, nchk, schk, sub_ierror, f_eval=xspl)
      if (my_ierror.eq.0) my_ierror = sub_ierror
   endif

   call spline_eval(spl, ikYDIR, nchk, schk, sub_ierror, f_eval=yspl)
   if (my_ierror.eq.0) my_ierror = sub_ierror

   call spline_eval(spl, ikZDIR, nchk, schk, sub_ierror, f_eval=zspl)
   if (my_ierror.eq.0) my_ierror = sub_ierror

   ! compute max distance between spline & linear interpolant

   dst2mx = 0d0
   imax = 0

   do ichk = 1, nchk
      dst2 = (yint(ichk)-yspl(ichk))**2 + (zint(ichk)-zspl(ichk))**2
      if (spl%has_xdata) dst2 = dst2 + (xint(ichk)-xspl(ichk))**2

      if (dst2.gt.dst2mx) then
         dst2mx = dst2
         imax = ichk
      endif

      if (ldebug.ge.2) then
         write(bufout,'(i6,4f12.6)') ichk, yint(ichk), zint(ichk), yspl(ichk), zspl(ichk)
         call write_log(1, bufout)
      endif
   enddo

   dist_max = sqrt(dst2mx)

   ! report on max.distance

   iseg = (imax-1) / k_chk + 1
   k    = mod(imax-1, k_chk)

   if (dist_max.gt.0.1d0) then
      write(bufout,'(a,f7.3,a)') ' WARNING: profile spline deviates up to', dist_max,                   &
             ' mm from input profile'
      call write_log(1, bufout)
   endif

   if (ldebug.ge.2) then
      write(bufout,'(a,f7.3,a,i4,4(a,f8.3),a)') ' max.distance is', dist_max,' at seg', iseg,           &
             ', (y,z)=(',yint(imax),',',zint(imax),') vs (',yspl(imax),',',zspl(imax),')'
      call write_log(1, bufout)
   endif

   deallocate(schk, xspl, yspl, zspl, xint, yint, zint)

end subroutine spline_check_updates

!------------------------------------------------------------------------------------------------------------

subroutine ppspline_make_simple_sec_intpol(tot_pnt, ip0, ip1, s, y, a3, a2, a1, a0, ierror)
!--function: compute section [ip0:ip1] of simple cubic interpolating spline {a,b,c,d} for given data {s,y}
!            using "free end" boundary conditions / natural spline.
   implicit none
!--subroutine arguments:
   integer,      intent(in)           :: tot_pnt   ! total number of points in data
   integer,      intent(in)           :: ip0, ip1  ! range of points for current section
   real(kind=8), intent(in)           :: s(tot_pnt), y(tot_pnt)
   real(kind=8), intent(out)          :: a3(tot_pnt), a2(tot_pnt), a1(tot_pnt), a0(tot_pnt)
   integer,      intent(out)          :: ierror
!--local variables:
   integer      :: iofs, npnt, nseg, neqs, iseg, i, j, k
   real(kind=8), dimension(:),   allocatable :: h, dmat, rhs
   real(kind=8), dimension(:,:), allocatable :: qmat, amat, lmat

   ierror = 0

   ! Note: using local numbering ipnt = 1:npnt, iseg = 1:nseg,
   !       subroutine arguments addressed as (iofs+ipnt), (iofs+iseg)

   iofs = ip0 - 1
   npnt = ip1 - ip0 + 1
   nseg = npnt - 1
   neqs = npnt - 2
   allocate(h(nseg), qmat(nseg,2), amat(neqs,2), dmat(neqs), lmat(neqs,1), rhs(neqs))

   ! no smoothing, kappa = lambda = 0

   ! compute distances between points, check that s-values are strictly increasing

   do iseg = 1, nseg

      h(iseg) = s(iofs+iseg+1) - s(iofs+iseg)

      if (h(iseg).le.0d0) then
         write(bufout,'(2a,i5,2(a,g12.4))') ' Internal error: ppspline_make_intpol: s-data should be', &
                ' strictly increasing. iseg=', iseg, ': s=',s(iofs+iseg),',',s(iofs+iseg+1)
         call write_log(1, bufout)
         call abort_run()
      endif
   enddo

   if (ldebug.ge.5) then
      call write_log(' ppspline_make_intpol: h=')
      do iseg = 1, nseg
         write(bufout,'(i6,f14.8)') iseg, h(iseg)
         call write_log(1, bufout)
      enddo
   endif

   ! build Q-matrix: R * b = Q^T * d; storing upper-triangular part only, diagonal in (:,1)

   do k = 1, nseg
      if (k.gt.1) then
         qmat(k,1) = -3d0 / h(k) - 3d0 / h(k-1)
      else
         qmat(k,1) = 0d0
      endif
      qmat(k,2) = 3d0 / h(k)
   enddo

   ! form the system matrix A = R   and  the right hand side rhs = Q' * y
   ! Note: qmat(i,:) stores column i of Q == row i of Q'

   do j = 1, neqs
      amat(j,1) = 2*(h(j) + h(j+1))
      if (j.le.neqs-1) then
         amat(j,2) =        h(j+1)
      else
         amat(j,2) = h(j+1)
      endif

      rhs(j) = qmat(j  ,2) * y(iofs+j  ) + qmat(j+1,1) * y(iofs+j+1) + qmat(j+1,2) * y(iofs+j+2)
   enddo

   if (ldebug.ge.5) then
      call write_log(' ppspline_make_intpol: j, amat(j,1:2), rhs(j)')
      do j = 1, neqs
         write(bufout,'(i6,3g14.6)') j, amat(j,1), amat(j,2), rhs(j)
         call write_log(1, bufout)
      enddo
   endif

   ! form LDLt factorization of symmetric matrix amat
   ! Note: dmat could share space with amat(:,1), lmat with amat(:,2)
   !       lmat(:,1) is first subdiagonal (i+1,i).

   do j = 1, neqs
      dmat(j)   = amat(j,1)
      if (j.gt.1) dmat(j) = dmat(j) - dmat(j-1) * lmat(j-1,1) * lmat(j-1,1)

      lmat(j,1) = amat(j,2) / dmat(j)
   enddo

   if (ldebug.ge.5) then
      call write_log(' ppspline_make_intpol: j, dmat(j), lmat(j,1)')
      do j = 1, neqs
         write(bufout,'(i6,2g14.6)') j, dmat(j), lmat(j,1)
         call write_log(1, bufout)
      enddo
   endif

   ! perform forward substitution, solving L * tmp1 = rhs; tmp1 shares space with rhs

   do j = 1, neqs
      if (j.gt.1) rhs(j) = rhs(j) - lmat(j-1,1) * rhs(j-1)
   enddo

   ! perform scaling, solving D * tmp2 = tmp1; tmp2 shares space with rhs

   do j = 1, neqs
      rhs(j) = rhs(j) / dmat(j)
   enddo

   ! perform back substitution, solving Lt * b = tmp2

   a2(iofs+npnt) = 0d0
   do j = neqs, 1, -1
      i = j + 1
      a2(iofs+i) = rhs(j)
      if (j.le.neqs-1) a2(iofs+i) = a2(iofs+i) - lmat(j,1) * a2(iofs+i+1)
   enddo
   a2(iofs+1) = 0d0

   if (ldebug.ge.5) then
      call write_log(' ppspline_make_intpol: j, a2')
      do j = 1, neqs
         write(bufout,'(i6,f14.8)') j, a2(iofs+j)
         call write_log(1, bufout)
      enddo
   endif

   ! set d = y

   do i = 1, npnt
      a0(iofs+i) = y(iofs+i)
   enddo

   ! compute spline parameters a and c

   do iseg = 1, nseg
      a3(iofs+iseg) = (a2(iofs+iseg+1)-a2(iofs+iseg)) / (3d0 * h(iseg))
   enddo
   a3(iofs+nseg+1) = 0d0

   do iseg = 1, nseg
      a1(iofs+iseg) = (a0(iofs+iseg+1)-a0(iofs+iseg)) / h(iseg) - a3(iofs+iseg) * h(iseg)**2            &
                                                                        - a2(iofs+iseg) * h(iseg)
   enddo
   a1(iofs+nseg+1) = 3d0 * a3(iofs+nseg) * h(nseg)**2 + 2d0 * a2(iofs+nseg) * h(nseg) +  a1(iofs+nseg)

end subroutine ppspline_make_simple_sec_intpol

!------------------------------------------------------------------------------------------------------------

subroutine ppspline_make_simple_sec_smoothing(tot_pnt, ip0, ip1, s, y, lambda, a3, a2, a1, a0, ierror, wgt)
!--function: compute section [ip0:ip1] of simple cubic smoothing spline {a,b,c,d} for given data {s,y},
!            weights wgt and parameter lambda, using "free end" boundary conditions / natural spline.
   implicit none
!--subroutine arguments:
   integer,      intent(in)           :: tot_pnt   ! total number of points in data
   integer,      intent(in)           :: ip0, ip1  ! range of points for current section
   real(kind=8), intent(in)           :: s(tot_pnt), y(tot_pnt), lambda
   real(kind=8), intent(out)          :: a3(tot_pnt), a2(tot_pnt), a1(tot_pnt), a0(tot_pnt)
   integer,      intent(out)          :: ierror
   real(kind=8), intent(in), optional :: wgt(tot_pnt)
!--local variables:
   integer      :: iofs, npnt, nseg, neqs, iseg, ipnt, i, j, k
   real(kind=8) :: kappa
   real(kind=8), dimension(:),   allocatable :: sigma, h, dmat, rhs
   real(kind=8), dimension(:,:), allocatable :: qmat, amat, lmat

   ierror = 0

   ! Note: using local numbering ipnt = 1:npnt, iseg = 1:nseg,
   !       subroutine arguments addressed as (iofs+ipnt), (iofs+iseg)

   iofs = ip0 - 1
   npnt = ip1 - ip0 + 1
   nseg = npnt - 1
   neqs = npnt - 2
   allocate(sigma(npnt), h(nseg), qmat(nseg,2), amat(neqs,3), dmat(neqs), lmat(neqs,2), rhs(neqs))

   ! set relative importance of smoothness of the result
   ! lambda == \hat{\lambda} cf. Wikipedia, weight of 2nd derivative when weight of data is 1, kappa == \mu

   kappa  = 2d0/3d0 * lambda

   ! compute inverse weight matrix sigma ! Note: should this be squared?

   if (present(wgt)) then
      do ipnt = 1, npnt
         sigma(ipnt) = 1d0 / max(1d-10, abs(wgt(iofs+ipnt)))
      enddo
   else
      sigma(1:npnt) = 1d0
   endif

   ! compute distances between points, check that s-values are strictly increasing

   do iseg = 1, nseg

      h(iseg) = s(iofs+iseg+1) - s(iofs+iseg)

      if (h(iseg).le.0d0) then
         write(bufout,'(2a,i5,2(a,g12.4))') 'Internal error: ppspline_make_smoothing: s-data should be ', &
                'strictly increasing. iseg=', iseg, ': s=',s(iofs+iseg),',',s(iofs+iseg+1)
         call write_log(1, bufout)
         call abort_run()
      endif
   enddo

   if (ldebug.ge.5) then
      call write_log(' ppspline_make_smoothing: h=')
      do iseg = 1, nseg
         write(bufout,'(i6,f14.8)') iseg, h(iseg)
         call write_log(1, bufout)
      enddo
   endif

   ! build Q-matrix: R * b = Q^T * d; storing upper-triangular part only, diagonal in (:,1)

   do k = 1, nseg
      if (k.gt.1) then
         qmat(k,1) = -3d0 / h(k) - 3d0 / h(k-1)
      else
         qmat(k,1) = 0d0
      endif
      qmat(k,2) = 3d0 / h(k)
   enddo

   ! form the system matrix A = R + kappa * Q' * sigma * Q  and  the right hand side rhs = Q' * y
   ! Note: qmat(i,:) stores column i of Q == row i of Q'

   do j = 1, neqs
      amat(j,1) = 2*(h(j) + h(j+1)) + kappa * (   qmat(j+1,1) * sigma(j+1) * qmat(j+1,1)        &
                                                + qmat(j  ,2) * sigma(j  ) * qmat(j  ,2)        &
                                                + qmat(j+1,2) * sigma(j+2) * qmat(j+1,2) )
      if (j.le.neqs-1) then
         amat(j,2) =        h(j+1)  + kappa * (   qmat(j+1,1) * sigma(j+1) * qmat(j+1,2)        &
                                                + qmat(j+1,2) * sigma(j+2) * qmat(j+2,1) )
      else
         amat(j,2) = h(j+1)
      endif
      if (j.le.neqs-2) then
         amat(j,3) =                  kappa *     qmat(j+1,2) * sigma(j+2) * qmat(j+2,2)
      else
         amat(j,3) = 0d0
      endif

      rhs(j) = qmat(j  ,2) * y(iofs+j  ) + qmat(j+1,1) * y(iofs+j+1) + qmat(j+1,2) * y(iofs+j+2)
   enddo

   if (ldebug.ge.5) then
      call write_log(' ppspline_make_smoothing: j, amat(j,1:3), rhs(j)')
      do j = 1, neqs
         write(bufout,'(i6,4g14.6)') j, amat(j,1), amat(j,2), amat(j,3), rhs(j)
         call write_log(1, bufout)
      enddo
   endif

   ! form LDLt factorization of symmetric matrix amat
   ! Note: dmat could share space with amat(:,1), lmat with and amat(:,2:3)
   !       lmat(:,1) is second subdiagonal (i+2,i), lmat(:,2) the first one (i+1,i).

   do j = 1, neqs
      dmat(j)   = amat(j,1)
      if (j.gt.2) dmat(j) = dmat(j) - dmat(j-2) * lmat(j-2,1) * lmat(j-2,1)
      if (j.gt.1) dmat(j) = dmat(j) - dmat(j-1) * lmat(j-1,2) * lmat(j-1,2)

      lmat(j,2) = amat(j,2)
      if (j.gt.1) lmat(j,2) = lmat(j,2) - dmat(j-1) * lmat(j-1,2) * lmat(j-1,1)
      lmat(j,2) = lmat(j,2) / dmat(j)

      lmat(j,1) = amat(j,3) / dmat(j)
   enddo

   if (ldebug.ge.5) then
      call write_log(' ppspline_make_smoothing: j, dmat(j), lmat(j,1:2)')
      do j = 1, neqs
         write(bufout,'(i6,3g14.6)') j, dmat(j), lmat(j,1), lmat(j,2)
         call write_log(1, bufout)
      enddo
   endif

   ! perform forward substitution, solving L * tmp1 = rhs; tmp1 shares space with rhs

   do j = 1, neqs
      if (j.gt.1) rhs(j) = rhs(j) - lmat(j-1,2) * rhs(j-1)
      if (j.gt.2) rhs(j) = rhs(j) - lmat(j-2,1) * rhs(j-2)
   enddo

   ! perform scaling, solving D * tmp2 = tmp1; tmp2 shares space with rhs

   do j = 1, neqs
      rhs(j) = rhs(j) / dmat(j)
   enddo

   ! perform back substitution, solving Lt * b = tmp2

   a2(iofs+npnt) = 0d0
   do j = neqs, 1, -1
      i = j + 1
      a2(iofs+i) = rhs(j)
      if (j.le.neqs-1) a2(iofs+i) = a2(iofs+i) - lmat(j,2) * a2(iofs+i+1)
      if (j.le.neqs-2) a2(iofs+i) = a2(iofs+i) - lmat(j,1) * a2(iofs+i+2)
   enddo
   a2(iofs+1) = 0d0

   if (ldebug.ge.5) then
      call write_log(' ppspline_make_smoothing: j, a2')
      do j = 1, neqs
         write(bufout,'(i6,f14.8)') j, a2(iofs+j)
         call write_log(1, bufout)
      enddo
   endif

   ! compute d = y - kappa * sigma * Q * b, with column i of Q contained in qmat(i,:)

   do i = 1, npnt
      a0(iofs+i) = y(iofs+i)
      if (i.ge.3                  ) a0(iofs+i) = a0(iofs+i) - kappa * sigma(i) * qmat(i-1,2) * a2(iofs+i-1)
      if (i.ge.2 .and. i.le.neqs+1) a0(iofs+i) = a0(iofs+i) - kappa * sigma(i) * qmat(i  ,1) * a2(iofs+i  )
      if (             i.le.neqs  ) a0(iofs+i) = a0(iofs+i) - kappa * sigma(i) * qmat(i  ,2) * a2(iofs+i+1)
   enddo

   ! compute spline parameters a and c

   do iseg = 1, nseg
      a3(iofs+iseg) = (a2(iofs+iseg+1)-a2(iofs+iseg)) / (3d0 * h(iseg))
   enddo
   a3(iofs+nseg+1) = 0d0

   do iseg = 1, nseg
      a1(iofs+iseg) = (a0(iofs+iseg+1)-a0(iofs+iseg)) / h(iseg) - a3(iofs+iseg) * h(iseg)**2            &
                                                                        - a2(iofs+iseg) * h(iseg)
   enddo
   a1(iofs+nseg+1) = 3d0 * a3(iofs+nseg) * h(nseg)**2 + 2d0 * a2(iofs+nseg) * h(nseg) +  a1(iofs+nseg)

end subroutine ppspline_make_simple_sec_smoothing

!------------------------------------------------------------------------------------------------------------

subroutine ppspline_make_simple_kink(npnt, s, y, lambda, a3, a2, a1, a0, nkink, ikinks, my_ierror, wgt)
!--function: compute simple cubic smoothing spline {a,b,c,d} for given data {s,y}, weights wgt and
!            parameter lambda, using "free end" boundary conditions / natural spline.
   implicit none
!--subroutine arguments:
   integer,      intent(in)           :: npnt, nkink
   real(kind=8), intent(in)           :: s(npnt), y(npnt), lambda
   integer,      intent(in)           :: ikinks(nkink)  ! including 1, npnt+1
   real(kind=8), intent(out)          :: a3(npnt), a2(npnt), a1(npnt), a0(npnt)
   integer,      intent(out)          :: my_ierror
   real(kind=8), intent(in), optional :: wgt(npnt)
!--local variables:
   integer          :: nsec, isec, ipnt0, ipnt1, sub_ierror
   character(len=9) :: strng

   my_ierror = 0

   nsec = nkink - 1
   if (ldebug.ge.3) then
      write(bufout,'(a,i2,a)') ' ppspline_make_simple_kink:',nsec,' sections:'
      call write_log(1, bufout)
      bufout(1) = ' ['
      do isec = 1, nsec
         if (ikinks(isec+1).le.99) then
            write(strng,'(i2,a,i2)') ikinks(isec),':',ikinks(isec+1)
         elseif (ikinks(isec+1).le.999) then
            write(strng,'(i3,a,i3)') ikinks(isec),':',ikinks(isec+1)
         else
            write(strng,'(i4,a,i4)') ikinks(isec),':',ikinks(isec+1)
         endif
         if (isec.eq.1) then
            bufout(1) = trim(bufout(1)) // trim(strng)
         else
            bufout(1) = trim(bufout(1)) // ', ' // trim(strng)
         endif
      enddo
      bufout(1) = trim(bufout(1)) // ']'
      call write_log(1, bufout)
   endif

   if (ikinks(1).ne.1 .or. ikinks(nkink).ne.npnt) then
      write(bufout,'(3(a,i4))') ' Internal error(ppspline_make_simple_kink): ikinks = [',ikinks(1),      &
                '...',ikinks(nkink),'], npnt=',npnt
      call write_log(1, bufout)
      my_ierror = 372
   endif
   do isec = 1, nsec
      if (ikinks(isec).gt.ikinks(isec+1)) then
         write(bufout,'(a,i2,3(a,i4))') ' Internal error(ppspline_make_simple_kink): isec ',isec,' = [', &
                ikinks(isec), '...',ikinks(isec+1),'], npnt=',npnt
         call write_log(1, bufout)
         my_ierror = 373
      endif
   enddo

   do isec = 1, nsec

      ! determine range of points for this section, kink-points included in both sections

      ipnt0 = ikinks(isec)
      ipnt1 = ikinks(isec+1)
  
      ! make spline for current section, overwriting end-point of previous section

      if (my_ierror.eq.0) then
         if (ldebug.ge.4) then
            write(bufout,'(a,i2,2(a,i4),a)') ' make_simple for section',isec,', ip=[',ipnt0,',',ipnt1,']'
            call write_log(1, bufout)
         endif
         if (abs(lambda).le.1d-10) then
            call ppspline_make_simple_sec_intpol(npnt, ipnt0, ipnt1, s, y, a3, a2, a1, a0, sub_ierror)
         else
            call ppspline_make_simple_sec_smoothing(npnt, ipnt0, ipnt1, s, y, lambda, a3, a2, a1, a0,   &
                        sub_ierror, wgt)
         endif
         my_ierror = sub_ierror
      endif

   enddo

end subroutine ppspline_make_simple_kink

!------------------------------------------------------------------------------------------------------------

subroutine ppspline_make_spline_kink(spl, nmeas, s_prf, x_prf, y_prf, z_prf, has_xdata, lambda, use_wgt, &
                nkink, ikinks, my_ierror)
!--function: compute parametric smoothing spline for 1-d grid (x(s),y(s),z(s))
   implicit none
!--subroutine arguments:
   type(t_spline)                     :: spl
   integer                            :: nmeas
   real(kind=8), intent(inout)        :: s_prf(nmeas), x_prf(nmeas), y_prf(nmeas), z_prf(nmeas) 
   logical,      intent(in)           :: has_xdata  ! whether x_prf is a dummy or needs to be used
   real(kind=8), intent(in)           :: lambda     ! weight of 2nd derivative, relative to wgt_in
   logical,      intent(in)           :: use_wgt    ! weigh data with spacing ds
   integer,      intent(in)           :: nkink
   integer,      intent(in)           :: ikinks(nkink) ! start/end of spline sections
   integer,      intent(out)          :: my_ierror
!--local variables:
   integer                    :: ip, k, sub_ierror
   real(kind=8), allocatable  :: wgt(:)

   my_ierror = 0

   ! set weights according to grid spacing ds

   allocate(wgt(nmeas))

   if (.not.use_wgt) then
      wgt(1:nmeas) = 1d0
   else
      wgt(1) = s_prf(2) - s_prf(1)
      do ip = 2, nmeas-1
         wgt(ip) = 0.5d0 * (s_prf(ip+1) - s_prf(ip-1))
      enddo
      wgt(nmeas) = s_prf(nmeas) - s_prf(nmeas-1)
   endif

   ! increase weights  20x at begin and end points
   ! could be done for ismooth=0, breaking compatibility with earlier versions

   if (use_wgt) then
      wgt(1)     = 20d0 * wgt(1)
      wgt(nmeas) = 20d0 * wgt(nmeas)
   endif

   ! increase weights 100x at interior break points

   do k = 2, nkink-1
      wgt(ikinks(k)) = 100d0 * wgt(ikinks(k))
   enddo

                        !!!! In this routine NTOT == NP_SPL and S_PRF == S_SPL !!!!!

   ! (re-)allocate spline aspl-arrays

   call spline_allocate(spl, nmeas)

   ! copy profile s_prf to spline s_spl

   spl%s(1:spl%npnt) = s_prf(1:nmeas)

   ! compute spline for x(s)

   if (has_xdata) then
      if (.not.spl%has_xdata) call spline_initx(spl)
      call ppspline_make_simple_kink(spl%npnt, spl%s, x_prf, lambda, spl%ax3, spl%ax2, spl%ax1, spl%ax0, &
                        nkink, ikinks, sub_ierror, wgt)
      if (my_ierror.eq.0) my_ierror = sub_ierror
   endif

   ! compute spline for y(s)

   call ppspline_make_simple_kink(spl%npnt, spl%s, y_prf, lambda, spl%ay3, spl%ay2, spl%ay1, spl%ay0,   &
                        nkink, ikinks, sub_ierror, wgt)
   if (my_ierror.eq.0) my_ierror = sub_ierror

   ! compute spline for z(s)

   call ppspline_make_simple_kink(spl%npnt, spl%s, z_prf, lambda, spl%az3, spl%az2, spl%az1, spl%az0,   &
                        nkink, ikinks, sub_ierror, wgt)
   if (my_ierror.eq.0) my_ierror = sub_ierror

   if (ldebug.ge.1 .and. ii_debug.ge.1 .and. ii_debug.le.spl%npnt) then
      write(bufout,'(a,i5,3(a,g12.4))') ' make_spline: ii=',ii_debug,', s=',spl%s(ii_debug),', y=',     &
                y_prf(ii_debug), ', z=',z_prf(ii_debug)
      call write_log(1, bufout)
      write(bufout,'(a,4g12.4)') ' ay0-ay3=', spl%ay0(ii_debug), spl%ay1(ii_debug), spl%ay2(ii_debug),  &
                spl%ay3(ii_debug)
      call write_log(1, bufout)
      write(bufout,'(a,4g12.4)') ' az0-az3=', spl%az0(ii_debug), spl%az1(ii_debug), spl%az2(ii_debug),  &
                spl%az3(ii_debug)
      call write_log(1, bufout)
   endif
   deallocate(wgt)

end subroutine ppspline_make_spline_kink

!------------------------------------------------------------------------------------------------------------

subroutine ppspline_make_spline_nokink(spl, nmeas, s_prf, x_prf, y_prf, z_prf, has_xdata, lambda,       &
                        use_wgt, my_ierror)
!--function: compute parametric smoothing spline for 1-d grid (x(s),y(s),z(s))
   implicit none
!--subroutine arguments:
   type(t_spline)                     :: spl
   integer                            :: nmeas
   real(kind=8), intent(inout)        :: s_prf(nmeas), x_prf(nmeas), y_prf(nmeas), z_prf(nmeas) 
   logical,      intent(in)           :: has_xdata  ! whether x_prf is a dummy or needs to be used
   real(kind=8), intent(in)           :: lambda     ! weight of 2nd derivative, using weight 1 for the data
   logical,      intent(in)           :: use_wgt    ! weigh data with spacing ds
   integer,      intent(out)          :: my_ierror
!--local variables:
   integer    :: nkink, ikinks(2)

   nkink  = 2
   ikinks = (/ 1, nmeas /)

   call ppspline_make_spline_kink(spl, nmeas, s_prf, x_prf, y_prf, z_prf, has_xdata, lambda, use_wgt,   &
                nkink, ikinks, my_ierror)

end subroutine ppspline_make_spline_nokink

!------------------------------------------------------------------------------------------------------------

end module m_spline_make
