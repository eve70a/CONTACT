!------------------------------------------------------------------------------------------------------------
! m_bspline_get - 1D B-spline evaluation functions
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_bspline_get
   use m_globals
   use m_markers
   use m_ptrarray
   use m_interp_1d
   use m_spline_def
   use m_bspline_def
   use m_bspline_1seg
   use m_bspline_make
   implicit none
   private

   ! Debugging for module m_bspline

   integer  :: ldebug    =  0    ! local level of debugging
   integer  :: ii_debug  = -1    ! output point for which detailed info is requested (-1 = none)
   integer  :: iel_debug = -1    ! input element for which detailed info is requested (-1 = none)
   public  bspline_set_debug
   public  bsplineget_set_debug

   ! Data type for 2D B-splines:

   public t_bspline2d

   ! Functions for the creation & evaluation of 2D B-splines:

   public  bspline_eval2d_list
   public  bspline_eval2d_prod
   private bspline_eval2d_inverse
   public  bspline_eval2d_inverse_list
   public  bspline_eval2d_inverse_prod

   private bspline_make_1d_ppspline_phase1
   private bspline_make_1d_ppspline_phase2
   private bspline_make_1d_ppspline_phase3
   public  bspline_make_1d_ppspline

   ! inverse for half-parametric spline x==u, (x,v) -> (y,z)

   public  bspline_get_z_at_xy
   public  bspline_get_z_at_xy_list
   public  bspline_get_z_at_xy_prod

   interface bspline_get_z_at_xy
      module procedure bspline_get_z_at_xy_list
      module procedure bspline_get_z_at_xy_prod
   end interface bspline_get_z_at_xy

   ! inverse for full-parametric spline (u,v) -> (x,y,z)

   public  bspline_get_xz_at_uy
   public  bspline_get_xz_at_uy_list
   public  bspline_get_xz_at_uy_prod

   interface bspline_get_xz_at_uy
      module procedure bspline_get_xz_at_uy_list
      module procedure bspline_get_xz_at_uy_prod
   end interface bspline_get_xz_at_uy

contains

!------------------------------------------------------------------------------------------------------------

subroutine bspline_set_debug(new_ldebug, new_ii_debug, new_iel_debug)
!--function: enable/disable debug output of all bspline routines
   implicit none
!--subroutine arguments:
   integer, intent(in)           :: new_ldebug       ! level of debug output required
   integer, intent(in), optional :: new_ii_debug     ! specific point of interest for debugging
   integer, intent(in), optional :: new_iel_debug    ! specific point of interest for debugging

   call bsplinedef_set_debug(new_ldebug, new_ii_debug, new_iel_debug)
   call bsplinemake_set_debug(new_ldebug, new_ii_debug, new_iel_debug)
   call bsplineget_set_debug(new_ldebug, new_ii_debug, new_iel_debug)

end subroutine bspline_set_debug

!------------------------------------------------------------------------------------------------------------

subroutine bsplineget_set_debug(new_ldebug, new_ii_debug, new_iel_debug)
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
      write(bufout,'(a,i3,2(a,i7))') ' bspline_get: debugging level =',ldebug,', ii_debug =', ii_debug, &
                ', iel_debug =', iel_debug
      call write_log(1, bufout)
   endif

end subroutine bsplineget_set_debug

!------------------------------------------------------------------------------------------------------------

subroutine bspline_eval2d_list(spl2d, nout, uout, vout, xout, yout, zout, my_ierror, exterval)
!--function: evaluate 2D tensor B-spline at nout positions (uout,vout) producing (xout,yout,zout)
   implicit none
!--subroutine arguments:
   type(t_bspline2d)                           :: spl2d
   integer,      intent(in)                    :: nout
   real(kind=8), intent(in)                    :: uout(nout), vout(nout)
   real(kind=8), intent(out)                   :: xout(nout), yout(nout), zout(nout)
   integer,      intent(out)                   :: my_ierror
   real(kind=8), intent(in), optional          :: exterval
!--local variables:
   logical              :: has_exter
   integer              :: iout, i, ii, j, jj, mask, istat, sub_ierror, nsplv
   real(kind=8)         :: ci_x, ci_y, ci_z
   integer,      dimension(:),   allocatable :: isegu, jsegv
   real(kind=8), dimension(:,:), allocatable :: b1u, b2u, b3u, b4u, b1v, b2v, b3v, b4v
   character(len=20)    :: namvar
   character(len=256)   :: errmsg

   my_ierror = 0
   has_exter = present(exterval)

   namvar = 'b1u-b4v'
   allocate(isegu(nout), b1u(nout,1), b2u(nout,2), b3u(nout,3), b4u(nout,4),                            &
            jsegv(nout), b1v(nout,1), b2v(nout,2), b3v(nout,3), b4v(nout,4), stat=istat, errmsg=errmsg)
   if (istat.ne.0) goto 99

   ! evaluate u-splines at positions uout

   if (my_ierror.eq.0) then
      call bspline_eval1d(spl2d%nknotu, spl2d%tui, tiny_dt, nout, uout, isegu, b1u, b2u, b3u, b4u,      &
                   ldebug, sub_ierror)
      my_ierror = sub_ierror
      if (my_ierror.ne.0) call write_log(' Bspline_eval2d: Error after eval1d(u)')
   endif

   ! evaluate v-splines at positions vout

   nsplv = spl2d%nsplv

   if (my_ierror.eq.0) then
      call bspline_eval1d(spl2d%nknotv, spl2d%tvj, tiny_dt, nout, vout, jsegv, b1v, b2v, b3v, b4v,      &
                   ldebug, sub_ierror)
      my_ierror = sub_ierror
      if (my_ierror.ne.0) call write_log(' Bspline_eval2d: Error after eval1d(v)')
   endif

   ! evaluate at all output positions

   if (my_ierror.eq.0) then
      do iout = 1, nout

         if (spl2d%has_xdata) xout(iout) = 0d0
         yout(iout) = 0d0
         zout(iout) = 0d0
         mask = 1

         if (vout(iout).lt.spl2d%tvj(4) .or. vout(iout).gt.spl2d%tvj(nsplv+1)) then

            mask = 0
            if (ldebug.ge.1) then
               write(bufout,'(a,i4,3(a,f12.4),a,i4)') ' vout(',iout,')=',vout(iout),' out of range [',  &
                           spl2d%tvj(4),',', spl2d%tvj(nsplv+1),'], jsegv=',jsegv(iout)
               call write_log(1, bufout)
            endif

         else

            do i = isegu(iout)-3, isegu(iout)
               ii = i - isegu(iout) + 4
               if (spl2d%has_xdata) ci_x = 0d0
               ci_y = 0d0
               ci_z = 0d0
               do j = jsegv(iout)-3, jsegv(iout)
                  jj = j - jsegv(iout) + 4
                  if (spl2d%has_xdata) ci_x = ci_x + spl2d%cij_x(i,j) * b4v(iout,jj)
                  ci_y = ci_y + spl2d%cij_y(i,j) * b4v(iout,jj)
                  ci_z = ci_z + spl2d%cij_z(i,j) * b4v(iout,jj)
                  mask = mask * spl2d%mask(i,j)
               enddo
               if (spl2d%has_xdata) xout(iout) = xout(iout) + ci_x * b4u(iout,ii)
               yout(iout) = yout(iout) + ci_y * b4u(iout,ii)
               zout(iout) = zout(iout) + ci_z * b4u(iout,ii)
            enddo

         endif

         if (mask.eq.0 .and. has_exter) then
            if (spl2d%has_xdata) xout(iout) = exterval
            yout(iout) = exterval
            zout(iout) = exterval
         endif
      enddo
   endif

   deallocate(isegu, b1u, b2u, b3u, b4u, jsegv, b1v, b2v, b3v, b4v)

   return

   ! error handling for memory allocation

99 continue

   write(bufout,'(3a,i4,2a)') ' ERROR in memory allocation (',trim(namvar),':',istat,'): ',trim(errmsg)
   call write_log(1, bufout)
   call abort_run()

end subroutine bspline_eval2d_list

!------------------------------------------------------------------------------------------------------------

subroutine bspline_eval2d_prod(spl2d, nuout, nvout, uout, vout, xout, yout, zout, my_ierror, exterval)
!--function: evaluate 2D tensor B-spline at nuout x nvout positions producing (xout,yout,zout)
   implicit none
!--subroutine arguments:
   type(t_bspline2d)                  :: spl2d
   integer,      intent(in)           :: nuout, nvout
   real(kind=8), intent(in)           :: uout(nuout), vout(nvout)
   real(kind=8), intent(out)          :: xout(nuout,nvout), yout(nuout,nvout), zout(nuout,nvout)
   integer,      intent(out)          :: my_ierror
   real(kind=8), intent(in), optional :: exterval
!--local variables:
   logical              :: has_exter
   integer              :: iout, jout, i, ii, j, jj, mask, istat, sub_ierror
   real(kind=8)         :: ci_x, ci_y, ci_z
   integer,      dimension(:),   allocatable :: isegu, jsegv
   real(kind=8), dimension(:,:), allocatable :: b1u, b2u, b3u, b4u, b1v, b2v, b3v, b4v
   character(len=20)    :: namvar
   character(len=256)   :: errmsg

   my_ierror = 0

   has_exter = present(exterval)

   namvar = 'b1u-b4v'
   allocate(isegu(nuout), b1u(nuout,1), b2u(nuout,2), b3u(nuout,3), b4u(nuout,4),                       &
            jsegv(nvout), b1v(nvout,1), b2v(nvout,2), b3v(nvout,3), b4v(nvout,4), stat=istat, errmsg=errmsg)
   if (istat.ne.0) goto 99

   if (ldebug.ge.2) then
      if (nuout.ge.5) then
         write(bufout,'(a,i4,a,3g12.4,a)') ' eval2d_prod: nu=',nuout,' positions uout=',                &
                (uout(i), i=1, 3), '...'
      else
         write(bufout,'(a,i4,a,4g12.4)')   ' eval2d_prod: nu=',nuout,' positions uout=',                &
                (uout(i), i=1, nuout)
      endif
      call write_log(1, bufout)
      if (nvout.ge.5) then
         write(bufout,'(a,i4,a,3g12.4,a)') '              nv=',nvout,' positions vout=',                &
             (vout(i), i=1, 3), '...'
      else
         write(bufout,'(a,i4,a,4g12.4)')   '              nv=',nvout,' positions vout=',                &
             (vout(i), i=1, nvout)
      endif
      call write_log(1, bufout)
   endif

   ! evaluate u-splines at positions uout

   if (my_ierror.eq.0) then
      call bspline_eval1d(spl2d%nknotu, spl2d%tui, tiny_dt, nuout, uout, isegu, b1u, b2u, b3u, b4u,     &
                   ldebug, sub_ierror)
      my_ierror = sub_ierror
      if (my_ierror.ne.0) call write_log(' Bspline_eval2d: Error after eval1d(u)')
   endif

   ! evaluate v-splines at positions vout

   if (my_ierror.eq.0) then
      call bspline_eval1d(spl2d%nknotv, spl2d%tvj, tiny_dt, nvout, vout, jsegv, b1v, b2v, b3v, b4v,     &
                   ldebug, sub_ierror)
      my_ierror = sub_ierror
      if (my_ierror.ne.0) call write_log(' Bspline_eval2d: Error after eval1d(v)')
   endif

   ! evaluate at all output positions

   if (my_ierror.eq.0) then

      do iout = 1, nuout
         do jout = 1, nvout

            if (spl2d%has_xdata) xout(iout,jout) = 0d0
            yout(iout,jout) = 0d0
            zout(iout,jout) = 0d0
            mask = 1
            do i = isegu(iout)-3, isegu(iout)
               ii = i - isegu(iout) + 4

               if (spl2d%has_xdata) ci_x = 0d0
               ci_y = 0d0
               ci_z = 0d0
               do j = jsegv(jout)-3, jsegv(jout)
                  jj = j - jsegv(jout) + 4
                  if (spl2d%has_xdata) ci_x = ci_x + spl2d%cij_x(i,j) * b4v(jout,jj)
                  ci_y = ci_y + spl2d%cij_y(i,j) * b4v(jout,jj)
                  ci_z = ci_z + spl2d%cij_z(i,j) * b4v(jout,jj)
                  mask = mask * spl2d%mask(i,j)

                  if (ldebug.ge.5 .and. iout.eq.ii_debug .and. jout.eq.iel_debug) then
                     write(bufout,'(4(a,i4),3(a,g14.6))') ' (iout,jout)=(',iout,',',jout,       &
                        '): ci_y += cij_y(',i,',',j,') * b4v =', spl2d%cij_y(i,j),' *', b4v(jout,jj)
                     call write_log(1, bufout)
                  endif
               enddo

               if (spl2d%has_xdata) xout(iout,jout) = xout(iout,jout) + ci_x * b4u(iout,ii)
               yout(iout,jout) = yout(iout,jout) + ci_y * b4u(iout,ii)
               zout(iout,jout) = zout(iout,jout) + ci_z * b4u(iout,ii)

               if (ldebug.ge.5 .and. iout.eq.ii_debug .and. jout.eq.iel_debug) then
                  write(bufout,'(2(a,i4),2(a,g14.6))') ' (iout,jout)=(',iout,',',jout,       &
                     '): yout += ci_y * b4u =', ci_y,' *', b4u(iout,ii)
                  call write_log(1, bufout)
               endif
            enddo

            if (mask.eq.0) then
               if (has_exter) then
                  if (spl2d%has_xdata) xout(iout,jout) = exterval
                  yout(iout,jout) = exterval
                  zout(iout,jout) = exterval
               endif
               if (ldebug.ge.3) then
                  write(bufout,'(2(a,i4),a)') ' eval2d_prod: zero mask at (i,j) = (',iout,',',jout,')'
                  call write_log(1, bufout)
               endif
            endif
         enddo
      enddo
   endif

   if (ldebug.ge.5) then
      iout = ii_debug
      jout = iel_debug
      write(bufout,'(2(a,i4),2(a,g14.6))') ' (iout,jout)=(',iout,',',jout,'): yout =',yout(iout,jout),  &
                ', zout =',zout(iout,jout)
      call write_log(1, bufout)
   endif

   deallocate(isegu, b1u, b2u, b3u, b4u, jsegv, b1v, b2v, b3v, b4v)

   return

   ! error handling for memory allocation

99 continue

   write(bufout,'(3a,i4,2a)') ' ERROR in memory allocation (',trim(namvar),':',istat,'): ',trim(errmsg)
   call write_log(1, bufout)
   call abort_run()

end subroutine bspline_eval2d_prod

!------------------------------------------------------------------------------------------------------------

subroutine bspline_eval2d_inverse(spl2d, noutu, nouty, noutv, uout, yout, vout, my_ierror)
!--function: inverse evaluation of 2D tensor B-spline, determining vout at given (uout,yout)
   implicit none
!--subroutine arguments:
   type(t_bspline2d)                  :: spl2d
   integer,      intent(in)           :: noutu, nouty, noutv
   real(kind=8), intent(in)           :: uout(noutu), yout(nouty)
   real(kind=8), intent(out)          :: vout(noutu,noutv)
   integer,      intent(out)          :: my_ierror
!--local variables:
   integer              :: nsplu, nsplv, i, ii, j, jseg, ioutv, ioutu, iouty, iouty0, iouty1, istat,    &
                           sub_ierror  ! jj
   logical              :: found, use_list
   real(kind=8)         :: v_a, v_b, f_a, f_b, c0, c1, c2, c3, vl !  d0, d1, d2, d3
   integer,      dimension(:),   allocatable :: isegu
   real(kind=8), dimension(:),   allocatable :: cj_y, ylow, yhig
   real(kind=8), dimension(:,:), allocatable :: b1u, b2u, b3u, b4u
   character(len=20)    :: namvar
   character(len=256)   :: errmsg

   ! two cases: tensor product noutu x nouty, noutv = nouty,
   !            list noutu = nouty, noutv = 1

   use_list = (noutu.eq.nouty .and. noutv.eq.1)
   if (.not.use_list .and. noutv.ne.nouty) then
      call write_log(' Internal error: bspline_eval2d_inverse: noutv <> nouty.')
      call abort_run()
   endif
   if (noutu.le.0 .or. nouty.le.0) then
      call write_log(' Internal error: bspline_eval2d_inverse: invalid noutu or nouty <= 0.')
      call abort_run()
   endif

   ! print information on inputs used

   if (ldebug.ge.1 .and. use_list) then
      write(bufout,'(a,i4,a)') ' Inverse spline evaluation (u,y) --> (v) at nout=',noutu,' pairs (u,y)'
      call write_log(1, bufout)
   elseif (ldebug.ge.1) then
      write(bufout,'(2(a,i4),a)') ' Inverse spline evaluation (u)x(y) --> (v) at nout=',noutu,' x',     &
                nouty,' positions'
      call write_log(1, bufout)
   endif

   if (ldebug.ge.2) then
      if (noutu.ge.5) then
         write(bufout,'(a,i4,a,3g12.4,a)') ' eval2d_inv : nu=',noutu,' positions uout=',                &
                (uout(i), i=1, 3), '...'
         call write_log(1, bufout)
      else
         write(bufout,'(a,i4,a,4g12.4)')   ' eval2d_inv : nu=',noutu,' positions uout=',                &
                (uout(i), i=1, noutu)
         call write_log(1, bufout)
      endif
      if (nouty.ge.5) then
         write(bufout,'(a,i4,a,3g12.4,a)') '              ny=',nouty,' positions yout=',                &
             (yout(i), i=1, 3), '...'
         call write_log(1, bufout)
      else
         write(bufout,'(a,i4,a,4g12.4)')   '              ny=',nouty,' positions yout=',                &
             (yout(i), i=1, nouty)
         call write_log(1, bufout)
      endif
   endif

   ! determine collocation matrix for computing coefficients cj_y

   nsplu  = spl2d%nsplu
   nsplv  = spl2d%nsplv

   namvar = 'b1u-yhig'
   allocate(isegu(noutu), b1u(noutu,1), b2u(noutu,2), b3u(noutu,3), b4u(noutu,4),                       &
            cj_y(nsplv), ylow(nsplv), yhig(nsplv), stat=istat, errmsg=errmsg)
   if (istat.ne.0) goto 99

   if (my_ierror.eq.0) then
      call bspline_eval1d(spl2d%nknotu, spl2d%tui, tiny_dt, noutu, uout, isegu, b1u, b2u, b3u, b4u,     &
                   ldebug, sub_ierror)
      my_ierror = sub_ierror
      if (my_ierror.ne.0) call write_log(' Bspline_eval2d_inv: Error after eval1d(u)')
   endif

   ! loop over output positions uout(1:noutu)

   do ioutu = 1 , noutu

      if (uout(ioutu).lt.spl2d%tui(4) .or. uout(ioutu).gt.spl2d%tui(nsplu+1)) then
         if (ldebug.ge.1) then
            write(bufout,'(a,i4,3(a,f12.4),a)') ' Warning: uout(',ioutu,')=',uout(ioutu),               &
                ' is out of range tu=[',spl2d%tui(4),',', spl2d%tui(nsplu+1),']'
            call write_log(1, bufout)
         endif
      endif

      ! determine spline coefficients cj_y(:) at position uout(ioutu)

      ! cj_y: [ nv, 1 ], Bmat: [ nout, nu ], cij_y: [ nu, nv ]
      ! cj_y = (Bmatu(ioutu,:) * spl2d.cij_y)^T

      do j = 1, nsplv
         cj_y(j) = 0d0
      enddo
      do i = isegu(ioutu)-3, isegu(ioutu)
         ii = i - isegu(ioutu) + 4
         do j = 1, nsplv
            cj_y(j) = cj_y(j) + spl2d%cij_y(i,j) * b4u(ioutu,ii)
         enddo
      enddo

      ! n+k knots gives n basisfunctions and n spline coefficients.
      ! basisfunction j becomes nonzero at knot j stays nonzero until knot j+k.
      ! values of y in segment j = [t_j,t_j+1) are defined by spline coefficients j-k+1 to j
      ! n+k knots gives n+k-1 segments. The first k-1 and last k-1 segments are outside the basic interval
      ! this gives n-k+1 segments in the basic interval with numbers k to n

      ! determine interval [ylow, yhig] for segments j, j = k to n

      ylow(1:3) = 1d0
      yhig(1:3) = 0d0
      do jseg = 4, nsplv
         ylow(jseg) = cj_y(jseg)
         ylow(jseg) = min(cj_y(jseg), cj_y(jseg-1))
         ylow(jseg) = min(cj_y(jseg), cj_y(jseg-2))
         ylow(jseg) = min(cj_y(jseg), cj_y(jseg-3))
         yhig(jseg) = cj_y(jseg)
         yhig(jseg) = max(cj_y(jseg), cj_y(jseg-1))
         yhig(jseg) = max(cj_y(jseg), cj_y(jseg-2))
         yhig(jseg) = max(cj_y(jseg), cj_y(jseg-3))
      enddo

      if (ldebug.ge.5) then
         do jseg = 4, nsplv
            if (ylow(jseg).le.yout(iouty) .and. yout(iouty).le.yhig(jseg)) then
               write(bufout,'(a,i4,4(a,f12.6),a)') ' seg j=',jseg,': v = [',spl2d%tvj(jseg),',',        &
                                        spl2d%tvj(jseg+1),', y in [',ylow(jseg),',',yhig(jseg),'] *'
            else
               write(bufout,'(a,i4,4(a,f12.6),a)') ' seg j=',jseg,': v = [',spl2d%tvj(jseg),',',        &
                                        spl2d%tvj(jseg+1),', y in [',ylow(jseg),',',yhig(jseg),']'
            endif
            call write_log(1, bufout)
         enddo
      endif

      if (use_list) then
         iouty0 = ioutu
         iouty1 = ioutu
      else
         iouty0 = 1
         iouty1 = nouty
      endif

      do iouty = iouty0, iouty1

         if (use_list) then
            ioutv = 1
         else
            ioutv = iouty
         endif

         ! consider all segments jseg = k : nspl that may contain yout(iouty) at uout(ioutu)
         ! there are nknot-2k+1 segments with numbers jseg = k .. nknot-k

         jseg  = 3
         found = .false.

         do while(.not.found .and. jseg.lt.nsplv)

            jseg = jseg + 1
            if (ylow(jseg).le.yout(iouty) .and. yout(iouty).le.yhig(jseg)) then

               ! possible intersection

               if (ldebug.ge.2) then
                  write(bufout,'(a,i4)') ' checking possible segment jseg=',jseg
                  call write_log(1, bufout)
               endif

               ! spline segment j: interval [t_j, t_{j+1}]

               v_a = spl2d%tvj(jseg)
               v_b = spl2d%tvj(jseg+1)

               ! PP-spline coefficients for segment:

               call bspline_get_ppcoef_1seg_orig(spl2d%nsplv, spl2d%tvj, cj_y, jseg, c0, c1, c2, c3, ldebug)
!              call bspline_get_ppcoef_1seg_modf(spl2d%nsplv, spl2d%tvj, cj_y, jseg, d0, d1, d2, d3, ldebug)

!              if (abs(d0-c0).gt.1d-2*abs(c0) .or. abs(d1-c1).gt.1d-2*abs(c1) .or. abs(d2-c2).gt.1d-8 .or. &
!                       abs(d3-c3).gt.1d-8) then
!                 write(bufout,'(a,4g16.8)') 'modf: d0-d3=',d0, d1, d2, d3
!                 call write_log(1, bufout)
!                 write(bufout,'(a,4g16.8)') 'orig: c0-c3=',c0, c1, c2, c3
!                 call write_log(1, bufout)
!                 write(bufout,'(2(a,i6))') 'nspl=',spl2d%nsplv,', jseg=',jseg
!                 call write_log(1, bufout)
!                 write(bufout,'(a,10g16.8)') 'tvj= ',(spl2d%tvj(jj), jj=jseg-4,min(spl2d%nsplv+4,jseg+4))
!                 call write_log(1, bufout)
!                 write(bufout,'(a,10g16.8)') 'cjy= ',(cj_y(jj), jj=jseg-3,jseg)
!                 call write_log(1, bufout)

!                 call bspline_get_ppcoef_1seg_orig(spl2d%nsplv, spl2d%tvj, cj_y, jseg, d0,d1,d2,d3, ldebug+5)
!                 call abort_run()
!              endif

               ! Solve cubic equation for segment:

               f_a = c0
               f_b = c3 * (v_b-v_a)**3 + c2 * (v_b-v_a)**2 + c1 * (v_b-v_a) + c0
               call solve_cubic_segm( v_a, v_b, f_a, c1, c2, c3, f_b, yout(iouty), vl, ldebug, sub_ierror)

               found = (vl.ge.0d0 .and. vl.le.v_b-v_a)
               if (found) vout(ioutu,ioutv) = v_a + vl

               if (found .and. ldebug.ge.4) then
                  write(bufout,'(a,i4,a,f9.3)') ' jseg=',jseg,': found v=', vout(ioutu,ioutv)
                  call write_log(1, bufout)
               elseif (.not.found .and. ldebug.ge.5) then
                  write(bufout,'(a,i4,a)') ' jseg=',jseg,': no solution'
                  call write_log(1, bufout)
               endif

            endif ! y in [ylow,yhig]
         enddo ! while(~found)

         if (found .and. ldebug.ge.2) then
            write(bufout,'(2(a,i4),3(a,f10.4),a,i4)') ' iout=',ioutu,',',iouty,', (u,y)=(',             &
                   uout(ioutu),',', yout(iouty), '): found vout=', vout(ioutu,ioutv),' in jseg=',jseg
            call write_log(1, bufout)
         elseif (.not.found) then
            vout(ioutu,ioutv) = -1d0
            if (ldebug.ge.2) then
               write(bufout,'(2(a,i4),3(a,f10.4))') ' No solution for iout=',ioutu,',',iouty,           &
                   ', (u,y)=(',uout(ioutu),',', yout(iouty),'), setting vout=',vout(ioutu,ioutv)
               call write_log(1, bufout)
            endif
         endif

      enddo ! for iouty
   enddo ! for ioutu

   deallocate(isegu, b1u, b2u, b3u, b4u, cj_y, ylow, yhig)

   return

   ! error handling for memory allocation

99 continue

   write(bufout,'(3a,i4,2a)') ' ERROR in memory allocation (',trim(namvar),':',istat,'): ',trim(errmsg)
   call write_log(1, bufout)
   call abort_run()

end subroutine bspline_eval2d_inverse

!------------------------------------------------------------------------------------------------------------

subroutine bspline_eval2d_inverse_list(spl2d, nout, uout, yout, vout, my_ierror)
!--function: inverse evaluation of 2D tensor B-spline, determining vout at given (uout,yout)
   implicit none
!--subroutine arguments:
   type(t_bspline2d)                  :: spl2d
   integer,      intent(in)           :: nout
   real(kind=8), intent(in)           :: uout(nout), yout(nout)
   real(kind=8), intent(out)          :: vout(nout)
   integer,      intent(out)          :: my_ierror
!--local variables
   integer      :: noutu, nouty, noutv
   real(kind=8) :: vout2d(nout,1)

   noutu = nout
   nouty = nout
   noutv = 1

   call bspline_eval2d_inverse(spl2d, noutu, nouty, noutv, uout, yout, vout2d, my_ierror)

   vout(1:nout) = vout2d(1:nout,1)

end subroutine bspline_eval2d_inverse_list

!------------------------------------------------------------------------------------------------------------

subroutine bspline_eval2d_inverse_prod(spl2d, noutu, nouty, uout, yout, vout, my_ierror)
!--function: inverse evaluation of 2D tensor B-spline, determining vout at given (uout,yout)
   implicit none
!--subroutine arguments:
   type(t_bspline2d)                  :: spl2d
   integer,      intent(in)           :: noutu, nouty
   real(kind=8), intent(in)           :: uout(noutu), yout(nouty)
   real(kind=8), intent(out)          :: vout(noutu,nouty)
   integer,      intent(out)          :: my_ierror
!--local variables

   call bspline_eval2d_inverse(spl2d, noutu, nouty, nouty, uout, yout, vout, my_ierror)

end subroutine bspline_eval2d_inverse_prod

!------------------------------------------------------------------------------------------------------------

subroutine bspline_make_1d_ppspline_phase1(bspl, nmeas, s_prf, ds_out, lambda, nkink, ikinks, naccel, &
                iaccel, my_ierror)
!--function: compute parametric least squares smoothing B-spline for 1-d grid (x(s),y(s),z(s)),
!            stored in PP-form as used in t_spline.
!            phase 1: preparations for knot-vector
   implicit none
!--subroutine arguments:
   type(t_bspline2d)                  :: bspl
   integer,      intent(in)           :: nmeas
   real(kind=8), intent(inout)        :: s_prf(nmeas)
   real(kind=8), intent(inout)        :: ds_out     ! target spacing for B-spline
   real(kind=8), intent(inout)        :: lambda     ! weight of 2nd derivative, relative to data wgt
   integer,      intent(in)           :: nkink
   integer,      intent(in)           :: ikinks(nkink)  ! kinks = jump in 1st derivative
   integer,      intent(in)           :: naccel
   integer,      intent(in)           :: iaccel(naccel) ! accelerations = jump in 2nd derivative
   integer,      intent(out)          :: my_ierror
!--local variables:
   real(kind=8), parameter    :: lmb_min =  1d-6
   logical                    :: knot_simple
   integer                    :: maxknot, j, istat, sub_ierror
   real(kind=8)               :: ds_out_min, totlen
   character(len=20)          :: namvar
   character(len=256)         :: errmsg

   my_ierror = 0

   associate(nknot => bspl%nknotu, nreduc => bspl%nreducu, nbrk => bspl%nbrku)

   ! Check requirements on kinks and accelerations

   call bspline_check_kink_accel(nmeas, nkink, ikinks, naccel, iaccel, sub_ierror)
   my_ierror = sub_ierror
   if (my_ierror.ne.0) call write_log(' Error after check_kink_accel')

   ! Choose between simple and advanced methods for knot placement

   knot_simple = .true.         ! default
   if (lambda.lt.lmb_min) knot_simple = .false.

   ! Check requirements on lambda: need smoothing to use knot_vector_simple

   if (knot_simple .and. lambda.lt.lmb_min) then
      write(bufout,'(2(a,g12.3))') ' WARNING: overriding lambda =',lambda,                              &
                '; B-spline smoothing uses lambda >=', lmb_min
      call write_log(1, bufout)
      lambda = lmb_min
   endif

   ! Check requirements on ds_out: must be large enough to avoid poor condition number 
   ! kappa = 20 * lambda / ds^6 <= 10^6, ds >= (20 * lambda / 10^6)^(1/6)

   ds_out_min = (20d0 * lambda * 1d-6)**(1d0/6d0)

   if (knot_simple .and. ds_out.lt.ds_out_min) then
      write(bufout,'(2(a,g12.3))') ' WARNING: overriding ds_out =',ds_out,                              &
                '; must be >= 0.1 * (20 lambda)^(1/6)=', ds_out_min
      call write_log(1, bufout)
      ds_out = ds_out_min
   endif

   ! Check expected number of knots for simple knot algorithm

   totlen = s_prf(nmeas) - s_prf(1)
   nknot  = nint(totlen / ds_out) + 3*(nkink+naccel+1)
   ds_out_min = totlen / (4999 - 3*(nkink+naccel+1))

   if (knot_simple .and. nknot.gt.4999) then
      write(bufout,'(a,g12.3,a,i7,a,g12.3)') ' WARNING: overriding ds_out =',ds_out,'; nknot=',nknot,   &
                ' seems unreasonably large. Using ds_out=', ds_out_min
      call write_log(1, bufout)
      ds_out = ds_out_min
   endif

   ! Note: ikinks includes 1 and nmeas

   if (ldebug.ge.1) then
      write(bufout,'(3(a,i4),a)') ' profile has ',nmeas,' points,', nkink-2,' kinks and', naccel,       &
                ' accelerations'
      call write_log(1, bufout)
      write(bufout,'(a,20i5)') ' kinks:', (ikinks(j), j=1,nkink)
      call write_log(1, bufout)
   endif

   ! make knot-vector tj and reduced knots t^*_j, implicit via keepj

   maxknot = max( nmeas, nint(totlen/ds_out) ) + 3*(nkink+naccel+1)
   if (.not.knot_simple) then
      ! after reduction, at most nmeas+4 knots remain; reduction removes 2 knots per adjacent kinks
      maxknot = min( nmeas + 4 + 2*nkink, maxknot )
      if (ldebug.ge.4) then
         write(bufout,'(5(a,i8))') ' clipping, new maxknot=',maxknot
         call write_log(1, bufout)
      endif
   endif

   if (maxknot.ge.9999) then
      write(bufout,'(a,i9,a)') ' WARNING: maxknot=',maxknot,' seems unreasonably large.'
      call write_log(1, bufout)
   endif

   namvar = 'tj, keepj'
   allocate(bspl%tui(maxknot), bspl%keeptu(maxknot), stat=istat, errmsg=errmsg)
   if (istat.ne.0) goto 99

   if (my_ierror.eq.0) then
      call bspline_make_knot_vector(nmeas, s_prf, nkink, ikinks, naccel, iaccel, knot_simple, ds_out,   &
                            maxknot, nknot, bspl%tui, bspl%keeptu, nreduc, ldebug, sub_ierror)
      my_ierror = sub_ierror
      if (my_ierror.ne.0) call write_log(' Error after make_knot_vector')
   endif

   if (my_ierror.eq.0) then
      if (ldebug.ge.1) then
         write(bufout,'(a,i4,a,i4)') ' knot vector tj has ',nknot,' knots, nmeas=',nmeas
         call write_log(1, bufout)
      endif
      if (ldebug.ge.4) then
         do j = 1, nknot
            write(bufout,'(a,i4,a,f12.6)') ' j =',j,': tj =',bspl%tui(j)
            call write_log(1, bufout)
         enddo
      endif
   endif

   ! store the breakpoints ksi, unique knots in tj

   nbrk   = nknot - 6 - 2*(nkink-2) - naccel
   namvar = 'bspl%ubrk'
   allocate(bspl%ubrk(nbrk), stat=istat, errmsg=errmsg)
   if (istat.ne.0) goto 99

   call bspline_make_breakpoints(nknot, bspl%tui, tiny_dt, nbrk, bspl%ubrk, ldebug)

   end associate

   return

   ! error handling for memory allocation

99 continue

   write(bufout,'(3a,i4,2a)') ' ERROR in memory allocation (',trim(namvar),':',istat,'): ',trim(errmsg)
   call write_log(1, bufout)
   call abort_run()

end subroutine bspline_make_1d_ppspline_phase1

!------------------------------------------------------------------------------------------------------------

subroutine bspline_make_1d_ppspline_phase2(bspl, nspl_rdc, lambda, wgt, &
                        nmeas, s_prf, x_prf, y_prf, z_prf, has_xdata, btw_xyz, my_ierror)
!--function: compute parametric least squares smoothing B-spline for 1-d grid (x(s),y(s),z(s)),
!            stored in PP-form as used in t_spline.
!            phase 2: set up (BTWB + l*DTCD) = BTWxyz, solve B-spline coefficients
   implicit none
!--subroutine arguments:
   type(t_bspline2d)                  :: bspl
   integer,      intent(in)           :: nspl_rdc, nmeas
   real(kind=8), intent(in)           :: lambda     ! weight of 2nd derivative, relative to data wgt
   real(kind=8), intent(in)           :: wgt(nmeas) ! weights, e.g. using spacing ds
   real(kind=8), intent(in)           :: s_prf(nmeas), x_prf(nmeas), y_prf(nmeas), z_prf(nmeas) 
   logical,      intent(in)           :: has_xdata  ! whether x_prf is a dummy or needs to be used
   real(kind=8), intent(out)          :: btw_xyz(bspl%nsplu,3)
   integer,      intent(out)          :: my_ierror
!--local variables:
   integer                    :: ip, istat, sub_ierror
   integer,      dimension(:),   allocatable :: jseg
   real(kind=8), dimension(:,:), allocatable :: b1, b2, b3, b4, cmat, dmat, btwb, dtcd
   character(len=20)          :: namvar
   character(len=256)         :: errmsg

   my_ierror = 0

   associate(nknot  => bspl%nknotu,  nreduc => bspl%nreducu, nspl_all => bspl%nsplu,                    &
             tj     => bspl%tui,     keepj  => bspl%keeptu )

   ! evaluate B-spline basis-functions at measurement locations

   namvar = 'jseg, b1-b4 (1)'
   allocate(jseg(nmeas), b1(nmeas,1), b2(nmeas,2), b3(nmeas,3), b4(nmeas,4), stat=istat, errmsg=errmsg)
   if (istat.ne.0) goto 99

   if (my_ierror.eq.0) then
      call bspline_eval1d(nknot, tj, tiny_dt, nmeas, s_prf, jseg, b1, b2, b3, b4, ldebug, sub_ierror)
      my_ierror = sub_ierror
      if (my_ierror.ne.0) call write_log(' Error after eval1d')
   endif

   ! form the product B * R, i.e. combine columns of b4 in case of segments with linear order

   if (my_ierror.eq.0 .and. nreduc.lt.nknot) then
      call bspline_reduce_bmat(nknot, keepj, nmeas, jseg, b4)
   endif

   ! form the product B^T * W * B

   namvar = 'btwb'
   allocate(btwb(4,nspl_rdc), stat=istat, errmsg=errmsg)
   if (istat.ne.0) goto 99

   if (my_ierror.eq.0) then
      call bspline_make_bt_w_b(nspl_rdc, nmeas, jseg, wgt, b4, btwb, ldebug)
   endif

   ! form the product B^T * W * [x y z]

   if (my_ierror.eq.0) then
      call bspline_make_bt_w_xyz(nspl_all, nspl_rdc, nmeas, jseg, wgt, b4, x_prf, y_prf, z_prf,         &
                has_xdata, btw_xyz, ldebug)
   endif
   deallocate(jseg, b1, b2, b3, b4)

   ! form the matrices C and D

   namvar = 'cmat, dmat'
   allocate(cmat(nknot-1,1), dmat(nknot-1,4), stat=istat, errmsg=errmsg)
   if (istat.ne.0) goto 99

   if (my_ierror.eq.0) then
      call bspline_make_c_and_d(nknot, tj, tiny_dt, cmat, dmat)
   endif

   ! form the product D * R, i.e. combine columns of dmat in case of segments with linear order

   if (my_ierror.eq.0 .and. nreduc.lt.nknot) then
      call bspline_reduce_dmat(nknot, keepj, dmat, ldebug)
   endif

   ! form the product D^T * C * D

   namvar = 'dtcd'
   allocate(dtcd(4,nspl_rdc), stat=istat, errmsg=errmsg)
   if (istat.ne.0) goto 99

   if (my_ierror.eq.0) then
      call bspline_make_dt_c_d(nknot, keepj, nspl_rdc, cmat, dmat, dtcd, ldebug)
   endif
   deallocate(cmat, dmat)

   ! solve the linear systems (B^T W B + lambda D^T C D) [coef_y, coef_z] = B^T W [y, z]
   ! solution overwrites right hand side in array btw_xyz

   if (my_ierror.eq.0) then
      call bspline_solve_coef1d_smooth(nspl_all, nspl_rdc, btwb, lambda, dtcd, btw_xyz, ldebug, sub_ierror)
      my_ierror = sub_ierror
      if (my_ierror.ne.0) call write_log(' Error after solve_coef')

      if (.false.) then
         call write_log('solution solve_coef1d_smooth:')
         write(bufout,'(2(a,i4))') 'nspl=',nspl_all,', rdc=',nspl_rdc
         call write_log(1, bufout)
         do ip = 1, nspl_rdc
            write(bufout,'(a,i4,2(a,f12.4))') 'i=',ip,': cy=',btw_xyz(ip,2),', cz=',btw_xyz(ip,3)
            call write_log(1, bufout)
         enddo
      endif
   endif

   ! form product coef = R^T * coef^*, ie expand reduced coefficients coef^* to full coefficients coef

   if (my_ierror.eq.0) then
      call bspline_expand_vector(nknot, keepj, nreduc, nspl_all, btw_xyz, ldebug)
   endif

   deallocate(btwb, dtcd)

   end associate

   return

   ! error handling for memory allocation

99 continue

   write(bufout,'(3a,i4,2a)') ' ERROR in memory allocation (',trim(namvar),':',istat,'): ',trim(errmsg)
   call write_log(1, bufout)
   call abort_run()

end subroutine bspline_make_1d_ppspline_phase2

!------------------------------------------------------------------------------------------------------------

subroutine bspline_make_1d_ppspline_phase3(bspl, coef_xyz, ppspl, my_ierror)
!--function: compute parametric least squares smoothing B-spline for 1-d grid (x(s),y(s),z(s)),
!            stored in PP-form as used in t_spline.
!            phase 3: convert knots + B-spline coefficients to PP-form
   implicit none
!--subroutine arguments:
   type(t_bspline2d)                  :: bspl
   type(t_spline)                     :: ppspl
   real(kind=8), intent(in)           :: coef_xyz(bspl%nsplu,3)
   integer,      intent(out)          :: my_ierror
!--local variables:
   integer                    :: istat, sub_ierror
   integer,      dimension(:),   allocatable :: jseg
   real(kind=8), dimension(:,:), allocatable :: b1, b2, b3, b4
   character(len=20)          :: namvar
   character(len=256)         :: errmsg

   my_ierror = 0

   associate(nknot => bspl%nknotu, tj => bspl%tui, nbrk => bspl%nbrku, sbrk => bspl%ubrk )

   ! create output spline

   if (my_ierror.eq.0) then

      call spline_allocate(ppspl, nbrk)
      ppspl%s(1:nbrk) = sbrk(1:nbrk)

   endif

   ! evaluate B-spline basis-functions at breakpoints ksi

   namvar = 'jseg, b1-b4 (2)'
   allocate(jseg(nbrk), b1(nbrk,1), b2(nbrk,2), b3(nbrk,3), b4(nbrk,4), stat=istat, errmsg=errmsg)
   if (istat.ne.0) goto 99

   if (my_ierror.eq.0) then
      call bspline_eval1d(nknot, tj, tiny_dt, nbrk, ppspl%s, jseg, b1, b2, b3, b4, ldebug, sub_ierror)
      my_ierror = sub_ierror
      if (my_ierror.ne.0) call write_log(' Error after eval1d(2)')
   endif

   ! determine PP-form: value of spline & derivatives at start of each segment

   if (my_ierror.eq.0) then
      call bspline_make_ppform(nknot, tj, tiny_dt, nbrk, jseg, b1, b2, b3, b4, bspl%nsplu, coef_xyz,    &
                        ppspl, sub_ierror)
      my_ierror = sub_ierror
      if (my_ierror.ne.0) call write_log(' Error after make_ppform')
   endif

   deallocate(jseg, b1, b2, b3, b4)

   end associate

   return

   ! error handling for memory allocation

99 continue

   write(bufout,'(3a,i4,2a)') ' ERROR in memory allocation (',trim(namvar),':',istat,'): ',trim(errmsg)
   call write_log(1, bufout)
   call abort_run()

end subroutine bspline_make_1d_ppspline_phase3

!------------------------------------------------------------------------------------------------------------

subroutine bspline_make_1d_ppspline(ppspl, nmeas, s_prf, x_prf, y_prf, z_prf, has_xdata, ds_out_arg,    &
                        lambda_arg, use_wgt, nkink, ikinks, naccel, iaccel, my_ierror, wgt_arg)
!--function: compute parametric least squares smoothing B-spline for 1-d grid (x(s),y(s),z(s)),
!            stored in PP-form as used in t_spline.
!            phase 3: convert knots + B-spline coefficients to PP-form
   implicit none
!--subroutine arguments:
   type(t_spline)                     :: ppspl
   integer                            :: nmeas
   real(kind=8), intent(inout)        :: s_prf(nmeas), x_prf(nmeas), y_prf(nmeas), z_prf(nmeas) 
   logical,      intent(in)           :: has_xdata  ! whether x_prf is a dummy or needs to be used
   real(kind=8), intent(in)           :: ds_out_arg ! target spacing for B-spline
   real(kind=8), intent(in)           :: lambda_arg ! weight of 2nd derivative, relative to data wgt
   logical,      intent(in)           :: use_wgt    ! weigh data with spacing ds
   integer,      intent(in)           :: nkink
   integer,      intent(in)           :: ikinks(nkink)  ! kinks = jump in 1st derivative
   integer,      intent(in)           :: naccel
   integer,      intent(in)           :: iaccel(naccel) ! accelerations = jump in 2nd derivative
   integer,      intent(out)          :: my_ierror
   real(kind=8), intent(in), optional :: wgt_arg(nmeas)
!--local variables
   type(t_bspline2d)    :: bspl
   integer              :: nspl_rdc, ip, istat
   real(kind=8)         :: lambda, ds_out
   real(kind=8), dimension(:),   allocatable :: wgt
   real(kind=8), dimension(:,:), allocatable :: coef_xyz
   character(len=20)          :: namvar
   character(len=256)         :: errmsg

   lambda = lambda_arg
   ds_out = ds_out_arg

   call bspline2d_nullify(bspl)

   ! perform checks, adjust lambda/ds if necessary, then create a knot vector and determine reduced knots

   call bspline_make_1d_ppspline_phase1(bspl, nmeas, s_prf, ds_out, lambda, nkink, ikinks, naccel,      &
                        iaccel, my_ierror)

   ! number of spline functions

   bspl%nsplu   = bspl%nknotu - 4
   nspl_rdc     = bspl%nreducu - 4

   ! set weights according to grid spacing ds

   namvar = 'wgt'
   allocate(wgt(nmeas), stat=istat, errmsg=errmsg)
   if (istat.ne.0) goto 99

   if (.not.use_wgt) then
      wgt(1:nmeas) = 1d0
   elseif (present(wgt_arg)) then
      wgt(1:nmeas) = wgt_arg(1:nmeas)
   else
      wgt(1) = s_prf(2) - s_prf(1)
      do ip = 2, nmeas-1
         wgt(ip) = 0.5d0 * (s_prf(ip+1) - s_prf(ip-1))
      enddo
      wgt(nmeas) = s_prf(nmeas) - s_prf(nmeas-1)
   endif

   ! allocate space for B-spline coefficients

   namvar = 'coef_xyz'
   allocate(coef_xyz(bspl%nsplu,3), stat=istat, errmsg=errmsg)
   if (istat.ne.0) goto 99

   ! set up and solve equations (BtWB + lambda DtCD) * c_xyz = BtW [x y z]

   call bspline_make_1d_ppspline_phase2(bspl, nspl_rdc, lambda, wgt, nmeas, s_prf, x_prf, y_prf,        &
                        z_prf, has_xdata, coef_xyz, my_ierror)

   ! determine PP-form and store in ppspl

   call bspline_make_1d_ppspline_phase3(bspl, coef_xyz, ppspl, my_ierror)

   ! clean up bspline

   call bspline2d_destroy(bspl)

   return

   ! error handling for memory allocation

99 continue

   write(bufout,'(3a,i4,2a)') ' ERROR in memory allocation (',trim(namvar),':',istat,'): ',trim(errmsg)
   call write_log(1, bufout)
   call abort_run()

end subroutine bspline_make_1d_ppspline

!------------------------------------------------------------------------------------------------------------

subroutine bspline_get_z_at_xy_list(spl2d, nout, xout, yout, zout, my_ierror, exterval)
!--function: for nout half-parametric positions (xout,yout), determine vout in the 2D tensor B-spline
!            and produce (zout)
   implicit none
!--subroutine arguments:
   type(t_bspline2d)                    :: spl2d
   integer,                intent(in)   :: nout
   real(kind=8),           intent(in)   :: xout(nout), yout(nout)
   real(kind=8),           intent(out)  :: zout(nout)
   integer,                intent(out)  :: my_ierror
   real(kind=8), optional, intent(in)   :: exterval
!--local variables:
   integer              :: iout
   real(kind=8)         :: vout(nout), ytmp(nout), xdum2d(1,1)

   if (spl2d%has_xdata) then
      call write_log(' Internal error: bspline_get_z_at_xy: cannot be used for full-parametric splines.')
      call abort_run()
   endif

   my_ierror = 0

   call bspline_eval2d_inverse_list(spl2d, nout, xout, yout, vout, my_ierror)

   if (my_ierror.eq.0) then
      call bspline_eval2d_list(spl2d, nout, xout, vout, xdum2d, ytmp, zout, my_ierror, exterval)
   endif

   if (ldebug.ge.2) then
      do iout = 1, nout
         write(bufout,'(a,i4,5(a,f12.6))') ' bspline_get_z_at_xy: iout=',iout,': x,y=',xout(iout),',', &
                     yout(iout), ': v=',vout(iout),', z=',zout(iout)
         call write_log(1, bufout)
      enddo
   endif

end subroutine bspline_get_z_at_xy_list

!------------------------------------------------------------------------------------------------------------

subroutine bspline_get_z_at_xy_prod(spl2d, nx, ny, xout, yout, zout, my_ierror, exterval)
!--function: for nx x ny half-parametric positions (xout) x (yout), determine vij in the 2D tensor
!            B-spline and produce (zout)
   implicit none
!--subroutine arguments:
   type(t_bspline2d)                    :: spl2d
   integer,                intent(in)   :: nx, ny
   real(kind=8),           intent(in)   :: xout(nx), yout(ny)
   real(kind=8),           intent(out)  :: zout(nx,ny)
   integer,                intent(out)  :: my_ierror
   real(kind=8), optional, intent(in)   :: exterval
!--local variables:
   integer              :: iu
   real(kind=8)         :: vij(nx,ny), vj(ny), ui(1), yj(ny), zj(ny), xdum2d(1,1)

   if (spl2d%has_xdata) then
      call write_log(' Internal error: bspline_get_z_at_xy: cannot be used for full-parametric splines.')
      call abort_run()
   endif

   my_ierror = 0

   call bspline_eval2d_inverse_prod(spl2d, nx, ny, xout, yout, vij, my_ierror)

   if (my_ierror.eq.0) then
      do iu = 1, nx
         ui(1)    = xout(iu)
         vj(1:ny) = vij(iu,1:ny)
         call bspline_eval2d_prod(spl2d, 1, ny, ui, vj, xdum2d, yj, zj, my_ierror, exterval)
         zout(iu,1:ny) = zj(1:ny)
      enddo
   endif

end subroutine bspline_get_z_at_xy_prod

!------------------------------------------------------------------------------------------------------------

subroutine bspline_get_xz_at_uy_list(spl2d, nout, uout, yout, xout, zout, my_ierror, exterval)
!--function: for nout full-parametric positions (uout,yout), determine vout in the 2D tensor B-spline
!            and produce (xout,zout)
   implicit none
!--subroutine arguments:
   type(t_bspline2d)                    :: spl2d
   integer,                intent(in)   :: nout
   real(kind=8),           intent(in)   :: uout(nout), yout(nout)
   real(kind=8),           intent(out)  :: xout(nout), zout(nout)
   integer,                intent(out)  :: my_ierror
   real(kind=8), optional, intent(in)   :: exterval
!--local variables:
   integer              :: iout
   real(kind=8)         :: vout(nout), ytmp(nout)

   my_ierror = 0

   call bspline_eval2d_inverse_list(spl2d, nout, uout, yout, vout, my_ierror)

   if (my_ierror.eq.0) then
      call bspline_eval2d_list(spl2d, nout, uout, vout, xout, ytmp, zout, my_ierror, exterval)
   endif

   if (ldebug.ge.2) then
      do iout = 1, nout
         write(bufout,'(a,i4,5(a,f12.6))') ' bspline_get_xz_at_uy: iout=',iout,': u,y=',uout(iout),',', &
                     yout(iout), ': v=',vout(iout),', x=',xout(iout),', z=',zout(iout)
         call write_log(1, bufout)
      enddo
   endif

end subroutine bspline_get_xz_at_uy_list

!------------------------------------------------------------------------------------------------------------

subroutine bspline_get_xz_at_uy_prod(spl2d, nu, ny, uout, yout, xout, zout, my_ierror, exterval)
!--function: for nu x ny full-parametric positions (uout) x (yout), determine vij in the 2D tensor
!            B-spline and produce (xout,zout)
   implicit none
!--subroutine arguments:
   type(t_bspline2d)                    :: spl2d
   integer,                intent(in)   :: nu, ny
   real(kind=8),           intent(in)   :: uout(nu), yout(ny)
   real(kind=8),           intent(out)  :: xout(nu,ny), zout(nu,ny)
   integer,                intent(out)  :: my_ierror
   real(kind=8), optional, intent(in)   :: exterval
!--local variables:
   integer              :: iu
   real(kind=8)         :: vij(nu,ny), vj(ny), ui(1), xj(ny), yj(ny), zj(ny)

   my_ierror = 0

   call bspline_eval2d_inverse_prod(spl2d, nu, ny, uout, yout, vij, my_ierror)

   if (my_ierror.eq.0) then
      do iu = 1, nu
         ui(1)    = uout(iu)
         vj(1:ny) = vij(iu,1:ny)
         call bspline_eval2d_prod(spl2d, 1, ny, ui, vj, xj, yj, zj, my_ierror, exterval)
         xout(iu,1:ny) = xj(1:ny)
         zout(iu,1:ny) = zj(1:ny)
      enddo
   endif

end subroutine bspline_get_xz_at_uy_prod

!------------------------------------------------------------------------------------------------------------

subroutine test_bspline(my_ierror)
!--function: compute least squares smoothing B-spline for flat+quadratic example
   implicit none
   integer,      intent(out)          :: my_ierror
!--local variables:
   integer,      parameter    :: idebug = 5, nmeas = 11
   real(kind=8)               :: si(nmeas), xi(nmeas), yi(nmeas), zi(nmeas)
   integer                    :: ip
   integer                    :: nkink, ikinks(5), naccel, iaccel(5)
   real(kind=8)               :: ds_out, lambda
   real(kind=8), dimension(:),   allocatable :: wgt
   type(t_spline)             :: spl

   my_ierror = 0
   call write_log(' test_bspline...')

   ! measurements according to flat + quadratic

   do ip = 1, nmeas
      si(ip) = real(ip-6)
      xi(ip) = 0d0
      yi(ip) = (max(0d0, si(ip)))**2
      zi(ip) = (max(0d0, si(ip)))**2
   enddo

   ! kinks and accelerations according to flat + quadratic

   nkink  = 2
   ikinks = (/ 1, nmeas, 0, 0, 0 /)   ! Note: ikinks includes 1 and nmeas
   naccel = 1
   iaccel = (/ 6, 0, 0, 0, 0 /)

   lambda = 1d0
   ds_out = 2.5d0

   ! set weights for testing

   allocate(wgt(nmeas))
   do ip = 1, nmeas
      wgt(ip) = 1d0 + 0.1d0 * real(ip)
   enddo

   call bspline_make_1d_ppspline(spl, nmeas, si, xi, yi, zi, .true., ds_out, lambda, .true., nkink,     &
                ikinks, naccel, iaccel, my_ierror, wgt_arg=wgt)

end subroutine test_bspline

!------------------------------------------------------------------------------------------------------------

end module m_bspline_get
