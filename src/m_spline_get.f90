!------------------------------------------------------------------------------------------------------------
! m_spline_get - evaluation routines for 1D spline in PP-form (piecewise polynomial)
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_spline_get
   use m_globals
   use m_markers
   use m_ptrarray
   use m_interp_1d
   use m_spline_def
   use m_spline_make
   use m_grids
   use m_gridfunc
   implicit none
   private

   ! Debugging for module m_spline_get

   public  spline_set_debug
   public  splineget_set_debug

   integer  :: ldebug    =  0    ! local level of debugging
   integer  :: ii_debug  = -1    ! output point for which detailed info is requested (-1 = none)
   integer  :: iel_debug = -1    ! input element for which detailed info is requested (-1 = none)

   public  spline_get_s_at_minval
   public  spline_get_s_at_minval_arr
   public  spline_get_s_at_minval_spl

   interface spline_get_s_at_minval
      module procedure spline_get_s_at_minval_arr
      module procedure spline_get_s_at_minval_spl
   end interface spline_get_s_at_minval

   public  spline_get_s_at_y
   public  spline_get_s_at_y_scalar
   public  spline_get_s_at_y_arr

   interface spline_get_s_at_y
      module procedure spline_get_s_at_y_scalar
      module procedure spline_get_s_at_y_arr
   end interface spline_get_s_at_y

   public  spline_get_s_at_z
   public  spline_get_s_at_z_scalar
   public  spline_get_s_at_z_arr

   interface spline_get_s_at_z
      module procedure spline_get_s_at_z_scalar
      module procedure spline_get_s_at_z_arr
   end interface spline_get_s_at_z

   public  spline_get_xz_at_y
   public  spline_get_xz_at_y_scalar
   public  spline_get_xz_at_y_arr
   public  spline_get_xz_at_y_gout

   interface spline_get_xz_at_y
      module procedure spline_get_xz_at_y_scalar
      module procedure spline_get_xz_at_y_arr
      module procedure spline_get_xz_at_y_gout
   end interface spline_get_xz_at_y

   public  spline_get_loc_xz_at_y
   public  spline_get_loc_xz_at_y_scalar
   public  spline_get_loc_xz_at_y_arr
   public  spline_get_loc_xz_at_y_gout

   interface spline_get_loc_xz_at_y
      module procedure spline_get_loc_xz_at_y_scalar
      module procedure spline_get_loc_xz_at_y_arr
      module procedure spline_get_loc_xz_at_y_gout
   end interface spline_get_loc_xz_at_y

   public  spline_get_xy_at_z
   public  spline_get_xy_at_z_scalar
   public  spline_get_xy_at_z_arr
   public  spline_get_xy_at_z_gout

   interface spline_get_xy_at_z
      module procedure spline_get_xy_at_z_scalar
      module procedure spline_get_xy_at_z_arr
      module procedure spline_get_xy_at_z_gout
   end interface spline_get_xy_at_z

   public  spline_get_dzdy_at_s
   public  spline_get_dzdy_at_s_scalar
   public  spline_get_dzdy_at_s_arr
   public  spline_get_dzdy_at_s_gfout

   interface spline_get_dzdy_at_s
      module procedure spline_get_dzdy_at_s_scalar
      module procedure spline_get_dzdy_at_s_arr
      module procedure spline_get_dzdy_at_s_gfout
   end interface spline_get_dzdy_at_s

   public  spline_get_dzdy_at_y
   public  spline_get_dzdy_at_y_scalar
   public  spline_get_dzdy_at_y_arr
   public  spline_get_dzdy_at_y_gfout

   interface spline_get_dzdy_at_y
      module procedure spline_get_dzdy_at_y_scalar
      module procedure spline_get_dzdy_at_y_arr
      module procedure spline_get_dzdy_at_y_gfout
   end interface spline_get_dzdy_at_y

   public  spline_get_alpha_at_s
   public  spline_get_alpha_at_s_scalar
   public  spline_get_alpha_at_s_arr
   public  spline_get_alpha_at_s_gfout

   interface spline_get_alpha_at_s
      module procedure spline_get_alpha_at_s_scalar
      module procedure spline_get_alpha_at_s_arr
      module procedure spline_get_alpha_at_s_gfout
   end interface spline_get_alpha_at_s

   public  spline_get_alpha_at_y
   public  spline_get_alpha_at_y_scalar
   public  spline_get_alpha_at_y_arr
   public  spline_get_alpha_at_y_gfout

   interface spline_get_alpha_at_y
      module procedure spline_get_alpha_at_y_scalar
      module procedure spline_get_alpha_at_y_arr
      module procedure spline_get_alpha_at_y_gfout
   end interface spline_get_alpha_at_y

   public  spline_get_curv_yz_at_s
   public  spline_get_curv_yz_at_s_scalar
   public  spline_get_curv_yz_at_s_arr

   interface spline_get_curv_yz_at_s
      module procedure spline_get_curv_yz_at_s_scalar
      module procedure spline_get_curv_yz_at_s_arr
   end interface spline_get_curv_yz_at_s

   public  spline_get_xyz_at_s
   public  spline_get_xyz_at_s_scalar
   public  spline_get_xyz_at_s_arr

   interface spline_get_xyz_at_s
      module procedure spline_get_xyz_at_s_scalar
      module procedure spline_get_xyz_at_s_arr
   end interface spline_get_xyz_at_s

   public  spline_get_xyz_at_miny
   public  spline_get_xyz_at_minz

   public  spline_get_dxyz_at_s
   public  spline_get_dxyz_at_s_scalar
   public  spline_get_dxyz_at_s_arr
   public  spline_get_dxyz_at_s_gfout

   interface spline_get_dxyz_at_s
      module procedure spline_get_dxyz_at_s_scalar
      module procedure spline_get_dxyz_at_s_arr
      module procedure spline_get_dxyz_at_s_gfout
   end interface spline_get_dxyz_at_s

   public  spline_get_nvec_at_s_gfout

   public  spline_get_pt_at_line_oyz

   private spline_test_data
   public  spline_test

!------------------------------------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------------------------------------

subroutine spline_set_debug(new_ldebug, new_ii_debug, new_iel_debug)
!--function: enable/disable debug output of all spline routines
   implicit none
!--subroutine arguments:
   integer, intent(in)           :: new_ldebug       ! level of debug output required
   integer, intent(in), optional :: new_ii_debug     ! specific point of interest for debugging
   integer, intent(in), optional :: new_iel_debug    ! specific point of interest for debugging

   call splinedef_set_debug(new_ldebug, new_ii_debug, new_iel_debug)
   call splinemake_set_debug(new_ldebug, new_ii_debug, new_iel_debug)
   call splineget_set_debug(new_ldebug, new_ii_debug, new_iel_debug)

end subroutine spline_set_debug

!------------------------------------------------------------------------------------------------------------

subroutine splineget_set_debug(new_ldebug, new_ii_debug, new_iel_debug)
!--function: enable/disable debug output of spline get-routines
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
      write(bufout,'(a,i3,2(a,i7))') ' spline-get:  debugging level =',ldebug,', ii_debug =',       &
                ii_debug,', iel_debug =', iel_debug
      call write_log(1, bufout)
   endif

end subroutine splineget_set_debug

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_s_at_minval_arr(npnt, s_spl, nam_f, a3, a2, a1, a0, smin, ierror)
!--function: for a simple spline, determine the s-position with minimum function value
!            --> locate segment, solve f'(s) = 0 with quadratic f'
   implicit none
!--subroutine arguments
   integer,          intent(in)   :: npnt
   character(len=1), intent(in)   :: nam_f
   real(kind=8),     intent(in)   :: s_spl(npnt), a3(npnt), a2(npnt), a1(npnt), a0(npnt)
   real(kind=8),     intent(out)  :: smin
   integer,          intent(out)  :: ierror
!--local variables
   integer      :: imin, iseg, sub_ierror

   ierror = 0

   ! Note: using a quick search that does not guarantee overall minimum function value.

   ! determine the break point with minimum value

   imin = idmin(npnt, a0(1:), 1)

   if (ldebug.ge.5) then
      write(bufout,'(3a,i4,2(3a,g12.4))') ' s_at_minval: minimum ',nam_f(1:1),'(i) for imin=',imin,     &
                ', ',nam_f(1:1),'(imin)=',a0(imin), ', d',nam_f(1:1),'/ds(imin)=',a1(imin)
      call write_log(1, bufout)
   endif

   ! df/ds < 0 - function is decreasing - take interval to the right;
   ! df/ds > 0 - function is increasing - interval to the left

   ! Note: f(s) can wiggle more within neighbouring segments than indicated in a0(1:)

   if (a1(imin).lt.0d0) then
      iseg = imin
   else
      iseg = imin - 1
   endif

   if (iseg.le.0) then                  ! minimum at first point of spline

      smin = s_spl(1)
      if (ldebug.ge.-1) then
         write(bufout,'(3a,f8.3)') ' s_at_minval: minimum ',nam_f(1:1),                                 &
                                             ' at start of spline function, s=',smin
         call write_log(1, bufout)
      endif

   elseif (iseg.ge.npnt) then           ! minimum at last point

      smin = s_spl(npnt)
      if (ldebug.ge.-1) then
         write(bufout,'(3a,f8.3)') ' s_at_minval: minimum ',nam_f(1:1),                                 &
                                             ' at end of spline function, s=', smin
         call write_log(1, bufout)
      endif

   else                                 ! minimum in interior

      call locate_one_extremum(npnt, s_spl, a3, a2, a1, a0, iseg, -1, smin, sub_ierror)
      if (sub_ierror.ne.0) smin = s_spl(imin)

   endif

end subroutine spline_get_s_at_minval_arr

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_s_at_minval_spl(spl, ikarg, smin, ierror)
!--function: determine s-position of minimum of parametric spline 'spl' in direction ikarg
   implicit none
!--subroutine arguments
   type(t_spline)                 :: spl
   integer,          intent(in)   :: ikarg
   real(kind=8),     intent(out)  :: smin
   integer,          intent(out)  :: ierror
!--local variables:
   character(len=1)              :: nam_f
   real(kind=8), dimension(:,:), pointer :: coef

   if (ikarg.eq.ikXDIR) then
      nam_f = 'x'
      coef  => spl%axspl
   elseif (ikarg.eq.ikYDIR) then
      nam_f = 'y'
      coef  => spl%ayspl
   elseif (ikarg.eq.ikZDIR) then
      nam_f = 'z'
      coef  => spl%azspl
   else
      call write_log('INTERNAL ERROR (spline_get_s_at_minval_spl): ikarg invalid.')
      call abort_run()
   endif

   if (associated(coef)) then
      call spline_get_s_at_minval_arr( spl%npnt, spl%s, nam_f, coef(:,4), coef(:,3), coef(:,2),         &
                        coef(:,1), smin, ierror)
   else
      if (.false.) call write_log('spline_get_s_at_minval_spl: coef not associated')
      smin = 0d0
   endif

end subroutine spline_get_s_at_minval_spl

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_s_at_y_scalar( spl, yout, sout, my_ierror )
!--function: interpolate simple spline y=y(s) to the requested y-coordinate.
   implicit none
!--subroutine arguments
   type(t_spline), intent(in)  :: spl
   real(kind=8),   intent(in)  :: yout
   real(kind=8),   intent(out) :: sout
   integer,        intent(out) :: my_ierror
!--local variables
   integer      :: nout, sub_ierror
   real(kind=8) :: sarr(1), yarr(1)

   my_ierror = 0
   if (.not.associated(spl%s) .or. .not.associated(spl%ayspl)) then
      call write_log(' Internal error: get_s_at_y: spline data have not yet been filled in.')
      call abort_run()
   endif

   nout    = 1
   yarr(1) = yout

   call spline_get_s_at_f( spl, ikYDIR, nout, yarr, sarr, sub_ierror )
   if (my_ierror.eq.0) my_ierror = sub_ierror

   sout   = sarr(1)

end subroutine spline_get_s_at_y_scalar

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_s_at_y_arr( spl, nout, yout, sout, my_ierror )
!--function: interpolate simple spline y=y(s) to the requested y-coordinates.
   implicit none
!--subroutine arguments
   type(t_spline), intent(in)  :: spl
   integer,        intent(in)  :: nout
   real(kind=8),   intent(in)  :: yout(nout)
   real(kind=8),   intent(out) :: sout(nout)
   integer,        intent(out) :: my_ierror
!--local variables
   integer      :: sub_ierror

   my_ierror = 0
   if (.not.associated(spl%s) .or. .not.associated(spl%ayspl)) then
      call write_log(' Internal error: get_s_at_y: spline data have not yet been filled in.')
      call abort_run()
   endif

   call spline_get_s_at_f( spl, ikYDIR, nout, yout, sout, sub_ierror)
   if (sub_ierror.ne.0) then
      if (my_ierror.eq.0) my_ierror = sub_ierror
      write(bufout, *) 'spline_get_s_at_y_arr: error', sub_ierror,' computing s at y'
      call write_log(1, bufout)
   endif

end subroutine spline_get_s_at_y_arr

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_s_at_z_scalar( spl, zout, sout, my_ierror )
!--function: interpolate simple spline z=z(s) to the requested z-coordinate.
   implicit none
!--subroutine arguments
   type(t_spline), intent(in)  :: spl
   real(kind=8),   intent(in)  :: zout
   real(kind=8),   intent(out) :: sout
   integer,        intent(out) :: my_ierror
!--local variables
   integer      :: nout, sub_ierror
   real(kind=8) :: sarr(1), zarr(1)

   my_ierror = 0
   if (.not.associated(spl%s) .or. .not.associated(spl%azspl)) then
      call write_log(' Internal error: get_s_at_z: spline data have not yet been filled in.')
      call abort_run()
   endif

   nout    = 1
   zarr(1) = zout

   call spline_get_s_at_f( spl, ikZDIR, nout, zarr, sarr, sub_ierror )
   if (my_ierror.eq.0) my_ierror = sub_ierror

   sout   = sarr(1)

end subroutine spline_get_s_at_z_scalar

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_s_at_z_arr( spl, nout, zout, sout, my_ierror )
!--function: interpolate simple spline z=z(s) to the requested z-coordinates.
   implicit none
!--subroutine arguments
   type(t_spline), intent(in)  :: spl
   integer,        intent(in)  :: nout
   real(kind=8),   intent(in)  :: zout(nout)
   real(kind=8),   intent(out) :: sout(nout)
   integer,        intent(out) :: my_ierror
!--local variables
   integer    :: sub_ierror

   my_ierror = 0
   if (.not.associated(spl%s) .or. .not.associated(spl%azspl)) then
      call write_log(' Internal error: get_s_at_z: spline data have not yet been filled in.')
      call abort_run()
   endif

   call spline_get_s_at_f( spl, ikZDIR, nout, zout, sout, sub_ierror )
   if (my_ierror.eq.0) my_ierror = sub_ierror

end subroutine spline_get_s_at_z_arr

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_xz_at_y_arr( spl, nout, yout, my_ierror, xout, zout, exterval )
!--function: interpolate spline [x(s),y(s),z(s)] of input grid (wheel/rail profile) to output positions yout
   implicit none
!--subroutine arguments
   type(t_spline),         intent(in)  :: spl
   integer,                intent(in)  :: nout
   real(kind=8),           intent(in)  :: yout(nout)
   integer,                intent(out) :: my_ierror
   real(kind=8), optional, intent(out) :: xout(nout), zout(nout)
   real(kind=8), optional, intent(in)  :: exterval
!--local variables
   integer                   :: sub_ierror
   real(kind=8), dimension(:), allocatable :: sout

   my_ierror = 0

   ! allocate work array s(nout)

   allocate(sout(nout))

   ! compute s-positions for output y-values

   if (ldebug.ge.4) call write_log(' spline_get_xz_at_y_arr: get_s_at_y_arr')
   call spline_get_s_at_y_arr( spl, nout, yout, sout, sub_ierror )
   if (sub_ierror.ne.0) then
      if (my_ierror.eq.0) my_ierror = sub_ierror
      write(bufout,*) 'spline_xz_at_y: error',sub_ierror,' in get_s_at_y_arr'
      call write_log(1, bufout)
   endif

   ! evaluate x-positions at output s-positions

   if (present(xout) .and. .not.spl%has_xdata) then

      ! ignore when nout=1: used in xz_at_y_scalar
      if (nout.gt.1) call write_log(' ERROR: spline_xz_y: input grid has no x-data')

   elseif (present(xout)) then

      if (ldebug.ge.4) call write_log(' spline_get_xz_at_y_arr: spline_eval for xout')
      call spline_eval( spl, ikXDIR, nout, sout, sub_ierror, exterval, f_eval=xout)
      if (my_ierror.eq.0) my_ierror = sub_ierror

   endif
   
   ! evaluate z-positions at output s-positions

   if (present(zout)) then

      if (ldebug.ge.4) call write_log(' spline_get_xz_at_y_arr: spline_eval for zout')
      call spline_eval(spl, ikZDIR, nout, sout, sub_ierror, exterval, f_eval=zout)
      if (sub_ierror.ne.0) then
         if (my_ierror.eq.0) my_ierror = sub_ierror
         write(bufout,*) 'spline_xz_at_y: error',sub_ierror,' in evaluation of spline for zout'
         call write_log(1, bufout)
      endif

   endif

end subroutine spline_get_xz_at_y_arr

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_xz_at_y_scalar( spl, yout, my_ierror, xout, zout, exterval )
!--function: interpolate spline [x(s),y(s),z(s)] of input grid (wheel/rail profile) to output position yout
   implicit none
!--subroutine arguments
   type(t_spline),         intent(in)  :: spl
   real(kind=8),           intent(in)  :: yout
   integer,                intent(out) :: my_ierror
   real(kind=8), optional, intent(out) :: xout, zout
   real(kind=8), optional, intent(in)  :: exterval
!--local variables
   integer      :: nout
   real(kind=8) :: xarr(1), yarr(1), zarr(1)

   nout    = 1
   yarr(1) = yout

   call spline_get_xz_at_y_arr( spl, 1, yarr, my_ierror, xarr, zarr, exterval )

   if (present(xout)) xout = xarr(1)
   if (present(zout)) zout = zarr(1)

end subroutine spline_get_xz_at_y_scalar

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_xz_at_y_gout( spl, gout, my_ierror, exterval )
!--function: interpolate spline [x(s),y(s),z(s)] of input grid (wheel/rail profile) to output grid (surface)
!            on input, gout%ntot and  gout%y are filled in
   implicit none
!--subroutine arguments
   type(t_spline)            :: spl
   type(t_grid)              :: gout
   integer,      intent(out) :: my_ierror
   real(kind=8), optional    :: exterval

   ! TODO: performance optimization in case spl has no x-data: compute one slice + extrude?

   if (spl%has_xdata) then
      call spline_get_xz_at_y_arr( spl, gout%ntot, gout%y, my_ierror, gout%x, gout%z, exterval )
   else
      call spline_get_xz_at_y_arr( spl, gout%ntot, gout%y, my_ierror, zout=gout%z, exterval=exterval )
   endif

end subroutine spline_get_xz_at_y_gout

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_loc_xz_at_y_arr( spl, m_loc, nout, yloc_out, my_ierror, xloc_out, zloc_out, exterval )
!--function: interpolate rotated spline [x_loc(s),y_loc(s),z_loc(s)] to output positions yloc_out
!            aimed at building the undeformed distance for a contact patch. m_loc == contact reference,
!            assuming that spline arc-lengths s are largely comparable to local y-coordinates
   implicit none
!--subroutine arguments
   type(t_spline),         intent(in)  :: spl            ! spline with curve in (y,z)-coordinates
   type(t_marker),         intent(in)  :: m_loc          ! local reference within (y,z) coordinate system 
   integer,                intent(in)  :: nout
   real(kind=8),           intent(in)  :: yloc_out(nout) ! local y coordinates where output is wanted
   integer,                intent(out) :: my_ierror
   real(kind=8), optional, intent(out) :: xloc_out(nout), zloc_out(nout) ! local x,z coordinates at yloc_out
   real(kind=8), optional, intent(in)  :: exterval
!--local variables
   integer              :: ipnt, ipnt0, ipnt1, imin
   real(kind=8)         :: safety, yloc_min, yloc_max, dst, dstmin, yloc0, yloc1, sr_ref
   type(t_vec)          :: ploc
   type(t_spline)       :: spl_loc
   ! real(kind=8), dimension(:), allocatable :: sout

   my_ierror = 0

   if (ldebug.ge.2) call write_log(' starting spline_get_loc_xz_at_y_arr...')
   if (ldebug.ge.5) call spline_print( spl, 'prr', 5)
   if (ldebug.ge.4) call marker_print( m_loc, 'm_loc', 5)

   ! use restricted search range in local y coordinates
   ! spline curve should not curve too much within the search range

   yloc_min = minval(yloc_out(1:nout))
   yloc_max = maxval(yloc_out(1:nout))
   safety   = 1.001d0
   dst      = yloc_max - yloc_min
   yloc_min = yloc_min - (safety - 1d0) * dst
   yloc_max = yloc_max + (safety - 1d0) * dst

   ! determine the spline point (yi,zi) closest to marker m_loc

   imin = 1
   dstmin = 1d10

   do ipnt = 1, spl%npnt
      ! TODO: speedup by bisection, and by using 1-norm distance
      dst = sqrt( (spl%ay0(ipnt)-m_loc%y())**2 + (spl%az0(ipnt)-m_loc%z())**2 )
      if (dst.lt.dstmin) then
         dstmin = dst
         imin   = ipnt
      endif
   enddo

   if (ldebug.ge.4) then
      write(bufout,'(a,f8.4,a,i4)') ' m_loc: minimum distance',dstmin,' at spline point ipnt=',imin
      call write_log(1, bufout)
   endif

   ! determine local y-coordinate at spline point imin

   ploc = vec_2loc(vec(0d0, spl%ay0(imin), spl%az0(imin)), m_loc)
   if (ldebug.ge.4) call vec_print(ploc, 'p(imin)', 5)

   ! determine the spline sr-coordinate at the local marker (local y == 0)
   ! TODO: refine using interpolation(?)
   !    spline sr = sr_ref + local y;   local y = 0 at spline sr = sr_ref

   sr_ref = spl%s(imin) - ploc%y()

   ! determine the range of spline segments encompassing local coordinates [ylmin,ylmax]
   !    ipnt0: highest i with local y < yloc_min,  spline s = sr_ref + local y < sr_ref + yloc_min
   !    ipnt1: lowest i  with local y > yloc_max,  spline s = sr_ref + local y > sr_ref + yloc_max

   ipnt0 = imin
   do while(ipnt0.gt.1 .and. spl%s(ipnt0).ge.sr_ref+yloc_min)
      ipnt0 = ipnt0 - 1
   enddo

   ipnt1 = imin
   do while(ipnt1.lt.spl%npnt .and. spl%s(ipnt1).le.sr_ref+yloc_max)
      ipnt1 = ipnt1 + 1
   enddo

   if (ldebug.ge.4) then
      write(bufout,'(2(a,i4,a,f9.3),a,f9.3)') ' spline s: i0=',ipnt0,', s0=', spl%s(ipnt0),', i1=',     &
           ipnt1, ', s1=',spl%s(ipnt1),', sref=',sr_ref
      call write_log(1, bufout)
   endif

   ! if the spline bends a lot, the projected distance in local y direction may be considerably
   !    smaller than the arc-length. Continue shifting ipnt0, ipnt1 using the local y-coordinates
   !    of spline points

   ploc = vec_2loc(vec(0d0, spl%ay0(ipnt0), spl%az0(ipnt0)), m_loc)
   do while(ipnt0.gt.1 .and. ploc%y().ge.yloc_min)
      ipnt0 = ipnt0 - 1
      ploc = vec_2loc(vec(0d0, spl%ay0(ipnt0), spl%az0(ipnt0)), m_loc)
   enddo
   yloc0 = ploc%y()

   ploc = vec_2loc(vec(0d0, spl%ay0(ipnt1), spl%az0(ipnt1)), m_loc)
   do while(ipnt1.lt.spl%npnt .and. ploc%y().le.yloc_max)
      ipnt1 = ipnt1 + 1
      ploc = vec_2loc(vec(0d0, spl%ay0(ipnt1), spl%az0(ipnt1)), m_loc)
   enddo
   yloc1 = ploc%y()

   ! TODO: warn when ipnt0 reaches the start and/or ipnt1 reaches the end of the spline
   ! TODO: stop at vertical points in the rotated profile (conformal example, large d_comb)

   if (ldebug.ge.4) then
      write(bufout,'(2(a,i4,a,f9.3),a,f9.3)') '  local s: i0=',ipnt0,', s0=', yloc0, ', i1=',ipnt1,     &
                ', s1=', yloc1,', sref=',0d0
      call write_log(1, bufout)
   endif

   ! create trimmed spline for selected ipnt0:ipnt1

   call spline_trim(spl, spl_loc, ipnt0, ipnt1)

   ! transform spline to local coordinates

   call spline_2loc(spl_loc, m_loc)
   if (ldebug.ge.5) call spline_print(spl_loc, 'spl_loc', 5)

   ! use non-rotated subroutine to fill the outputs

   call spline_get_xz_at_y_arr(spl_loc, nout, yloc_out, my_ierror, xloc_out, zloc_out, exterval)

   ! cleanup work variable

   call spline_destroy(spl_loc)

end subroutine spline_get_loc_xz_at_y_arr

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_loc_xz_at_y_scalar( spl, m_loc, yloc_out, my_ierror, xloc_out, zloc_out, exterval )
!--function: interpolate rotated spline [x_loc(s),y_loc(s),z_loc(s)] to output positions yloc_out
   implicit none
!--subroutine arguments
   type(t_spline),         intent(in)  :: spl
   type(t_marker),         intent(in)  :: m_loc          ! local reference within (y,z) coordinate system 
   real(kind=8),           intent(in)  :: yloc_out
   integer,                intent(out) :: my_ierror
   real(kind=8), optional, intent(out) :: xloc_out, zloc_out
   real(kind=8), optional, intent(in)  :: exterval
!--local variables
   integer      :: nout
   real(kind=8) :: xarr(1), yarr(1), zarr(1)

   nout    = 1
   yarr(1) = yloc_out

   call spline_get_loc_xz_at_y_arr( spl, m_loc, 1, yarr, my_ierror, xarr, zarr, exterval )

   if (present(xloc_out)) xloc_out = xarr(1)
   if (present(zloc_out)) zloc_out = zarr(1)

end subroutine spline_get_loc_xz_at_y_scalar

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_loc_xz_at_y_gout( spl, m_loc, gout, my_ierror, exterval )
!--function: interpolate rotated spline [x(s),y(s),z(s)] (wheel/rail profile) to output grid (surface)
!            on input, gout%ntot and  gout%y are filled in
   implicit none
!--subroutine arguments
   type(t_spline)              :: spl
   type(t_marker), intent(in)  :: m_loc          ! local reference within (y,z) coordinate system 
   type(t_grid)                :: gout
   integer,        intent(out) :: my_ierror
   real(kind=8),   optional    :: exterval

   ! TODO: performance optimization in case spl has no x-data: compute one slice + extrude?

   if (spl%has_xdata) then
      call spline_get_loc_xz_at_y_arr( spl, m_loc, gout%ntot, gout%y, my_ierror, gout%x, gout%z, exterval )
   else
      call spline_get_loc_xz_at_y_arr( spl, m_loc, gout%ntot, gout%y, my_ierror, zloc_out=gout%z,       &
                exterval=exterval )
   endif

end subroutine spline_get_loc_xz_at_y_gout

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_xy_at_z_scalar( spl, zout, my_ierror, xout, yout, exterval )
!--function: interpolate spline [x(s),y(s),z(s)] of input grid (wheel/rail profile) to output position zout
   implicit none
!--subroutine arguments
   type(t_spline),         intent(in)  :: spl
   real(kind=8),           intent(in)  :: zout
   integer,                intent(out) :: my_ierror
   real(kind=8), optional, intent(out) :: xout, yout
   real(kind=8), optional, intent(in)  :: exterval
!--local variables
   integer      :: nout
   real(kind=8) :: xarr(1), yarr(1), zarr(1)

   nout    = 1
   zarr(1) = zout

   call spline_get_xy_at_z_arr( spl, nout, zarr, my_ierror, xarr, yarr, exterval )

   if (present(xout)) xout = xarr(1)
   if (present(yout)) yout = yarr(1)

end subroutine spline_get_xy_at_z_scalar

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_xy_at_z_arr( spl, nout, zout, my_ierror, xout, yout, exterval )
!--function: interpolate spline [x(s),y(s),z(s)] of input grid (wheel/rail profile) to output positions zout
   implicit none
!--subroutine arguments
   type(t_spline),         intent(in)  :: spl
   integer,                intent(in)  :: nout
   real(kind=8),           intent(in)  :: zout(nout)
   integer,                intent(out) :: my_ierror
   real(kind=8), optional, intent(out) :: xout(nout), yout(nout)
   real(kind=8), optional, intent(in)  :: exterval
!--local variables
   integer                   :: sub_ierror
   real(kind=8), dimension(:), allocatable :: sout

   my_ierror = 0

   ! allocate work array s(nout)

   allocate(sout(nout))

   ! compute s-positions for output z-values

   if (ldebug.ge.3) call write_log(' spline_get_xy_at_z_arr: get_s_at_z_arr')
   call spline_get_s_at_z_arr( spl, nout, zout, sout, sub_ierror )
   if (my_ierror.eq.0) my_ierror = sub_ierror

   ! evaluate x-positions at output s-positions

   if (present(xout) .and. .not.spl%has_xdata) then

      ! ignore xout, esp. for use by _scalar and _gout versions

   elseif (present(xout)) then

      if (ldebug.ge.3) call write_log(' spline_get_xy_at_z_arr: spline_eval for xout')
      call spline_eval( spl, ikXDIR, nout, sout, sub_ierror, exterval, f_eval=xout)
      if (my_ierror.eq.0) my_ierror = sub_ierror

   endif
   
   ! evaluate y-positions at output s-positions

   if (present(yout)) then

      if (ldebug.ge.3) call write_log(' spline_get_xy_at_z_arr: spline_eval for yout')
      call spline_eval( spl, ikYDIR, nout, sout, sub_ierror, exterval, f_eval=yout)
      if (my_ierror.eq.0) my_ierror = sub_ierror

   endif

end subroutine spline_get_xy_at_z_arr

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_xy_at_z_gout( spl, gout, my_ierror, exterval )
!--function: interpolate spline [x(s),y(s),z(s)] of input grid (wheel/rail profile) to output grid (surface)
   implicit none
!--subroutine arguments
   type(t_spline)            :: spl
   type(t_grid)              :: gout
   integer,      intent(out) :: my_ierror
   real(kind=8), optional    :: exterval
!--local variables
   integer                   :: nout

   nout = gout%ntot
   call spline_get_xy_at_z_arr( spl, nout, gout%z, my_ierror, gout%x, gout%y, exterval )

end subroutine spline_get_xy_at_z_gout

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_dzdy_at_s_scalar( spl, sout, dzdy, my_ierror )
!--function: compute profile slope dzdy of spline [y(s),z(s)] (input grid) at output position sout
   implicit none
!--subroutine arguments
   type(t_spline), intent(in)  :: spl
   real(kind=8),   intent(in)  :: sout
   real(kind=8),   intent(out) :: dzdy
   integer,        intent(out) :: my_ierror
!--local variables
   integer                   :: nout
   real(kind=8)              :: sarr(1), dzdy_arr(1)

   nout    = 1
   sarr(1) = sout

   call spline_get_dzdy_at_s_arr( spl, nout, sarr, dzdy_arr, my_ierror )

   dzdy    = dzdy_arr(1)

end subroutine spline_get_dzdy_at_s_scalar

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_dzdy_at_s_arr( spl, nout, sout, dzdy, my_ierror )
!--function: compute profile slope dzdy of spline [y(s),z(s)] (input grid) at output positions sout
   implicit none
!--subroutine arguments
   type(t_spline), intent(in)  :: spl
   integer,        intent(in)  :: nout
   real(kind=8),   intent(in)  :: sout(nout)
   real(kind=8),   intent(out) :: dzdy(nout)
   integer,        intent(out) :: my_ierror
!--local variables
   integer                   :: sub_ierror, iout
   real(kind=8), dimension(:), allocatable :: dy_ds, dz_ds

   my_ierror = 0

   ! allocate work arrays dy_ds(nout), dz_ds(nout)

   allocate(dy_ds(nout), dz_ds(nout))

   ! evaluate derivatives dy/ds and dz/ds at output s-positions

   call spline_eval(spl, ikYDIR, nout, sout, sub_ierror, f1_eval=dy_ds)
   if (my_ierror.eq.0) my_ierror = sub_ierror
   call spline_eval(spl, ikZDIR, nout, sout, sub_ierror, f1_eval=dz_ds)
   if (my_ierror.eq.0) my_ierror = sub_ierror

   do iout = 1, nout
      if (abs(dy_ds(iout)).lt.1d-6) then
         dzdy(iout) = sign(1d6, dz_ds(iout))    ! 1d6 * sign(dz_ds)
      else
         dzdy(iout) = dz_ds(iout) / dy_ds(iout)
      endif
   enddo
   deallocate(dy_ds, dz_ds)

end subroutine spline_get_dzdy_at_s_arr

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_dzdy_at_s_gfout( spl, gf_out, my_ierror )
!--function: compute profile slope dzdy of spline [y(s),z(s)] (input grid) to output grid-function gfout
   implicit none
!--subroutine arguments
   type(t_spline)               :: spl
   type(t_gridfnc3)             :: gf_out
   integer,         intent(out) :: my_ierror
!--local variables
   integer                  :: nout
   type(t_grid), pointer    :: gout

   ! get output grid gout

   gout => gf_out%grid
   nout = gout%ntot

   ! evaluate derivatives dz/dy at grid s-positions, store in column vy of grid-function gf_out

   call spline_get_dzdy_at_s_arr( spl, nout, gout%s_prf, gf_out%vy, my_ierror )

end subroutine spline_get_dzdy_at_s_gfout

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_dzdy_at_y_scalar( spl, yout, dzdy, my_ierror )
!--function: compute profile slope dzdy of spline [y(s),z(s)] (input grid) at output position yout
   implicit none
!--subroutine arguments
   type(t_spline), intent(in)  :: spl
   real(kind=8),   intent(in)  :: yout
   real(kind=8),   intent(out) :: dzdy
   integer,        intent(out) :: my_ierror
!--local variables
   integer                   :: nout
   real(kind=8)              :: yarr(1), dzdy_arr(1)

   nout    = 1
   yarr(1) = yout

   call spline_get_dzdy_at_y_arr( spl, nout, yarr, dzdy_arr, my_ierror )

   dzdy    = dzdy_arr(1)

end subroutine spline_get_dzdy_at_y_scalar

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_dzdy_at_y_arr( spl, nout, yout, dzdy, my_ierror )
!--function: compute profile slope dzdy of spline [y(s),z(s)] (input grid) at output positions yout
   implicit none
!--subroutine arguments
   type(t_spline), intent(in)  :: spl
   integer,        intent(in)  :: nout
   real(kind=8),   intent(in)  :: yout(nout)
   real(kind=8),   intent(out) :: dzdy(nout)
   integer,        intent(out) :: my_ierror
!--local variables
   integer                   :: sub_ierror
   real(kind=8), dimension(:), allocatable :: sout

   ! allocate work array s(nout)

   allocate(sout(nout))

   ! compute s-positions for output y-values

   call spline_get_s_at_y_arr( spl, nout, yout, sout, sub_ierror )
   if (my_ierror.eq.0) my_ierror = sub_ierror

   ! evaluate derivatives dy/ds and dz/ds at output s-positions, form dz/dy

   call spline_get_dzdy_at_s_arr( spl, nout, sout, dzdy, sub_ierror )
   if (my_ierror.eq.0) my_ierror = sub_ierror

   deallocate(sout)

end subroutine spline_get_dzdy_at_y_arr

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_dzdy_at_y_gfout( spl, gf_out, my_ierror )
!--function: compute profile slope dzdy of spline [y(s),z(s)] (input grid) to output grid-function gfout
   implicit none
!--subroutine arguments
   type(t_spline)               :: spl
   type(t_gridfnc3)             :: gf_out
   integer,         intent(out) :: my_ierror
!--local variables
   integer                  :: nout
   type(t_grid), pointer    :: gout

   ! get output grid gout

   gout => gf_out%grid
   nout = gout%ntot

   ! evaluate derivative dz/dy at grid s-positions, store in column vy of grid-function gf_out

   call spline_get_dzdy_at_y_arr( spl, nout, gout%y, gf_out%vy, my_ierror )

end subroutine spline_get_dzdy_at_y_gfout

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_alpha_at_s_scalar( spl, sout, alpha, my_ierror )
!--function: compute profile inclination angle alpha [rad] of spline [y(s),z(s)] (input grid) at output
!            position sout
   implicit none
!--subroutine arguments
   type(t_spline), intent(in)  :: spl
   real(kind=8),   intent(in)  :: sout
   real(kind=8),   intent(out) :: alpha
   integer,        intent(out) :: my_ierror
!--local variables
   integer                   :: nout
   real(kind=8)              :: sarr(1), alpha_arr(1)

   nout    = 1
   sarr(1) = sout

   call spline_get_alpha_at_s_arr( spl, nout, sarr, alpha_arr, my_ierror )

   alpha    = alpha_arr(1)

end subroutine spline_get_alpha_at_s_scalar

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_alpha_at_s_arr( spl, nout, sout, alpha, my_ierror )
!--function: compute profile inclination angle alpha [rad] of spline [y(s),z(s)] (input grid) at output
!            positions sout
   implicit none
!--subroutine arguments
   type(t_spline), intent(in)  :: spl
   integer,        intent(in)  :: nout
   real(kind=8),   intent(in)  :: sout(nout)
   real(kind=8),   intent(out) :: alpha(nout)
   integer,        intent(out) :: my_ierror
!--local variables
   integer                   :: sub_ierror, iout
   real(kind=8), dimension(:), allocatable :: dy_ds, dz_ds

   my_ierror = 0

   ! allocate work arrays dy_ds(nout), dz_ds(nout)

   allocate(dy_ds(nout), dz_ds(nout))

   ! evaluate derivatives dy/ds and dz/ds at output s-positions

   call spline_eval(spl, ikYDIR, nout, sout, sub_ierror, f1_eval=dy_ds)
   if (my_ierror.eq.0) my_ierror = sub_ierror
   call spline_eval(spl, ikZDIR, nout, sout, sub_ierror, f1_eval=dz_ds)
   if (my_ierror.eq.0) my_ierror = sub_ierror

   ! compute inclination angle alpha = atan2(dz, dy) at output s-positions

   do iout = 1, nout
      alpha(iout) = atan2(dz_ds(iout), dy_ds(iout))
   enddo
   deallocate(dy_ds, dz_ds)

end subroutine spline_get_alpha_at_s_arr

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_alpha_at_s_gfout( spl, gf_out, my_ierror )
!--function: compute profile inclination angle alpha [rad] of spline [y(s),z(s)] (input grid) to output
!            grid-function gfout
   implicit none
!--subroutine arguments
   type(t_spline)               :: spl
   type(t_gridfnc3)             :: gf_out
   integer,         intent(out) :: my_ierror
!--local variables
   integer                  :: nout
   type(t_grid), pointer    :: gout

   ! get output grid gout

   gout => gf_out%grid
   nout = gout%ntot

   ! compute inclination angle alpha = atan2(dz, dy) at output s-positions and store in y-column of gf_out

   call spline_get_alpha_at_s_arr( spl, nout, gout%s_prf, gf_out%vy, my_ierror )

end subroutine spline_get_alpha_at_s_gfout

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_alpha_at_y_scalar( spl, yout, alpha, my_ierror )
!--function: compute profile inclination angle alpha [rad] of spline [y(s),z(s)] (input grid) at output
!            position yout
   implicit none
!--subroutine arguments
   type(t_spline), intent(in)  :: spl
   real(kind=8),   intent(in)  :: yout
   real(kind=8),   intent(out) :: alpha
   integer,        intent(out) :: my_ierror
!--local variables
   integer                   :: nout
   real(kind=8)              :: yarr(1), alpha_arr(1)

   nout    = 1
   yarr(1) = yout

   call spline_get_alpha_at_y_arr( spl, nout, yarr, alpha_arr, my_ierror )

   alpha    = alpha_arr(1)

end subroutine spline_get_alpha_at_y_scalar

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_alpha_at_y_arr( spl, nout, yout, alpha, my_ierror )
!--function: compute profile inclination angle alpha [rad] of spline [y(s),z(s)] (input grid) at output
!            positions yout
   implicit none
!--subroutine arguments
   type(t_spline), intent(in)  :: spl
   integer,        intent(in)  :: nout
   real(kind=8),   intent(in)  :: yout(nout)
   real(kind=8),   intent(out) :: alpha(nout)
   integer,        intent(out) :: my_ierror
!--local variables
   integer                   :: sub_ierror
   real(kind=8), dimension(:), allocatable :: sout

   ! allocate work array s(nout)

   allocate(sout(nout))

   ! compute s-positions for output y-values

   call spline_get_s_at_y_arr( spl, nout, yout, sout, sub_ierror )
   if (my_ierror.eq.0) my_ierror = sub_ierror

   ! evaluate derivatives dy/ds and dz/ds at output s-positions, form alpha=atan2(dz, dy)

   call spline_get_alpha_at_s_arr( spl, nout, sout, alpha, sub_ierror )
   if (my_ierror.eq.0) my_ierror = sub_ierror

   deallocate(sout)

end subroutine spline_get_alpha_at_y_arr

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_alpha_at_y_gfout( spl, gf_out, my_ierror )
!--function: compute profile inclination angle alpha [rad] of spline [y(s),z(s)] (input grid) to output
!            grid-function gfout
   implicit none
!--subroutine arguments
   type(t_spline)               :: spl
   type(t_gridfnc3)             :: gf_out
   integer,         intent(out) :: my_ierror
!--local variables
   integer                  :: nout
   type(t_grid), pointer    :: gout

   ! get output grid gout

   gout => gf_out%grid
   nout = gout%ntot

   ! evaluate derivatives dy/ds and dz/ds at output y-positions, form alpha=atan2(dz, dy) in gf_out

   call spline_get_alpha_at_y_arr( spl, nout, gout%y, gf_out%vy, my_ierror )

end subroutine spline_get_alpha_at_y_gfout

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_curv_yz_at_s_scalar( spl, sout, kappa, ierror )
!--function: compute profile curvature kappa of spline [y(s),z(s)] (input grid) at output position sout
   implicit none
!--subroutine arguments
   type(t_spline), intent(in)  :: spl
   real(kind=8),   intent(in)  :: sout
   real(kind=8),   intent(out) :: kappa
   integer,        intent(out) :: ierror
!--local variables
   integer      :: nout
   real(kind=8) :: sarr(1), karr(1)

   nout    = 1
   sarr(1) = sout

   call spline_get_curv_yz_at_s_arr( spl, nout, sarr, karr, ierror )

   kappa  = karr(1)

end subroutine spline_get_curv_yz_at_s_scalar

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_curv_yz_at_s_arr( spl, nout, sout, kappa, my_ierror )
!--function: compute profile curvature kappa of spline [y(s),z(s)] (input grid) at output positions sout
   implicit none
!--subroutine arguments
   type(t_spline), intent(in)  :: spl
   integer,        intent(in)  :: nout
   real(kind=8),   intent(in)  :: sout(nout)
   real(kind=8),   intent(out) :: kappa(nout)
   integer,        intent(out) :: my_ierror
!--local variables
   integer                   :: sub_ierror, iout
   real(kind=8)              :: numer, denom
   real(kind=8), dimension(:), allocatable :: dy_ds, d2y_ds, dz_ds, d2z_ds

   my_ierror = 0

   ! allocate work arrays dy_ds(nout), dz_ds(nout)

   allocate(dy_ds(nout), d2y_ds(nout), dz_ds(nout), d2z_ds(nout))

   ! evaluate derivatives dy/ds, d2y/ds2, dz/ds and d2z/ds2 at output s-positions

   call spline_eval( spl, ikYDIR, nout, sout, sub_ierror, f1_eval=dy_ds, f2_eval=d2y_ds)
   if (my_ierror.eq.0) my_ierror = sub_ierror
   call spline_eval( spl, ikZDIR,  nout, sout, sub_ierror, f1_eval=dz_ds, f2_eval=d2z_ds)
   if (my_ierror.eq.0) my_ierror = sub_ierror

   ! form signed curvature k = (y'z" - z'y") / (y'^2 + z'^2)^3/2

   do iout = 1, nout
      numer = dy_ds(iout) * d2z_ds(iout) - dz_ds(iout) * d2y_ds(iout)
      denom = (dy_ds(iout)**2 + dz_ds(iout)**2)**1.5d0
      if (abs(denom).lt.1d-6) then
         kappa(iout) = sign(1d6, numer)    ! 1d6 * (y'z" - z'y")
      else
         kappa(iout) = numer / denom
      endif

      if (ldebug.ge.3) then
         write(bufout,'(a,i4,a,f9.4,5(a,g12.4))') ' i=',iout,': s=',sout,', y''=',dy_ds(iout),           &
                ', y"=', d2y_ds(iout), ', z''=',dz_ds(iout),', z"=', d2z_ds(iout), ', kappa=',kappa(iout)
         call write_log(1, bufout)
      endif
   enddo

   deallocate(dy_ds, d2y_ds, dz_ds, d2z_ds)

end subroutine spline_get_curv_yz_at_s_arr

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_xyz_at_s_scalar( spl, sout, my_ierror, xout, yout, zout)
!--function: evaluate spline [x(s),y(s),z(s)] of input grid (wheel/rail profile) to output position sout
   implicit none
!--subroutine arguments
   type(t_spline),         intent(in)  :: spl
   real(kind=8),           intent(in)  :: sout
   integer,                intent(out) :: my_ierror
   real(kind=8), optional, intent(out) :: xout, yout, zout
!--local variables
   integer      :: nout
   real(kind=8) :: sarr(1), xarr(1), yarr(1), zarr(1)

   ! evaluate x-, y- and z-positions at output s-position

   nout    = 1
   sarr(1) = sout

   call spline_get_xyz_at_s_arr( spl, nout, sarr, my_ierror, xarr, yarr, zarr )

   if (present(xout)) xout = xarr(1)
   if (present(yout)) yout = yarr(1)
   if (present(zout)) zout = zarr(1)

end subroutine spline_get_xyz_at_s_scalar

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_xyz_at_s_arr( spl, nout, sout, my_ierror, xout, yout, zout )
!--function: evaluate spline [x(s),y(s),z(s)] of input grid (wheel/rail profile) at output positions sout
   implicit none
!--subroutine arguments
   type(t_spline),         intent(in)  :: spl
   integer,                intent(in)  :: nout
   real(kind=8),           intent(in)  :: sout(nout)
   integer,                intent(out) :: my_ierror
   real(kind=8), optional, intent(out) :: xout(nout), yout(nout), zout(nout)
!--local variables
   integer                   :: sub_ierror

   my_ierror = 0

   ! evaluate x-, y- and z-positions at output s-positions

   if (present(xout) .and. spl%has_xdata) then
      call spline_eval( spl, ikXDIR, nout, sout, sub_ierror, f_eval=xout)
      if (my_ierror.eq.0) my_ierror = sub_ierror
   elseif (present(xout)) then
      xout(1:nout) = 0d0
   endif

   if (present(yout)) then
      call spline_eval( spl, ikYDIR, nout, sout, sub_ierror, f_eval=yout)
      if (my_ierror.eq.0) my_ierror = sub_ierror
   endif

   if (present(zout)) then
      call spline_eval( spl, ikZDIR, nout, sout, sub_ierror, f_eval=zout)
      if (my_ierror.eq.0) my_ierror = sub_ierror
   endif

end subroutine spline_get_xyz_at_s_arr

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_xyz_at_miny( spl, iwin0, iwin1, my_ierror, sout, xout, yout, zout)
!--function: determine s-position with minimum y-value in window in spline, evaluate [x(s),y(s),z(s)]
   implicit none
!--subroutine arguments
   type(t_spline),         intent(in)  :: spl
   integer,                intent(in)  :: iwin0, iwin1
   integer,                intent(out) :: my_ierror
   real(kind=8), optional, intent(out) :: sout, xout, yout, zout
!--local variables
   integer      :: i0, i1, nout, sub_ierror
   real(kind=8) :: sarr(1), xarr(1), yarr(1), zarr(1)

   my_ierror = 0
   if (.not.associated(spl%s)) then
      call write_log(' Internal error: xyz_miny: profile s-coordinate has not yet been filled in.')
      call abort_run()
   endif

   ! determine s-position with minimum y-value in the window [iwin0, iwin1]

   i0 = max(1, iwin0)
   i1 = min(spl%npnt, iwin1)

   if (ldebug.ge.4) then
      write(bufout,'(2(a,i4),a)') ' get_xyz_at_miny: window i=[',i0,',',i1,']'
      call write_log(1, bufout)
   endif

   if (i0.gt.i1) then           ! empty window
      sout = 0d0
   elseif (i0.eq.i1) then       ! single point
      sout = spl%s(i0)
   else                         ! regular window
      call spline_get_s_at_minval(i1-i0+1, spl%s(i0:), 'y', spl%ay3(i0:), spl%ay2(i0:),     &
                                           spl%ay1(i0:), spl%ay0(i0:), sout, sub_ierror)
      if (my_ierror.eq.0) my_ierror = sub_ierror
      if (ldebug.ge.3) then
         write(bufout,'(a,g12.4)') ' get_xyz_at_miny: minimum y at sout=', sout
         call write_log(1, bufout)
      endif
   endif

   ! evaluate x-, y- and z-positions at output s-position

   nout    = 1
   sarr(1) = sout

   call spline_get_xyz_at_s_arr( spl, nout, sarr, my_ierror, xarr, yarr, zarr )

   if (present(xout)) xout = xarr(1)
   if (present(yout)) yout = yarr(1)
   if (present(zout)) zout = zarr(1)

end subroutine spline_get_xyz_at_miny

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_xyz_at_minz( spl, my_ierror, sout, xout, yout, zout)
!--function: determine s-position with minimum z-value in spline, evaluate [x(s),y(s),z(s)]
   implicit none
!--subroutine arguments
   type(t_spline),         intent(in)  :: spl
   integer,                intent(out) :: my_ierror
   real(kind=8), optional, intent(out) :: sout, xout, yout, zout
!--local variables
   integer      :: nout, sub_ierror
   real(kind=8) :: sarr(1), xarr(1), yarr(1), zarr(1)

   my_ierror = 0
   if (.not.associated(spl%s)) then
      call write_log(' Internal error: xyz_minz: profile s-coordinate has not yet been filled in.')
      call abort_run()
   endif

   ! determine s-position with minimum z-value

   call spline_get_s_at_minval(spl, ikZDIR, sout, sub_ierror)
   if (my_ierror.eq.0) my_ierror = sub_ierror

   ! evaluate x-, y- and z-positions at output s-position

   nout    = 1
   sarr(1) = sout

   call spline_get_xyz_at_s_arr( spl, nout, sarr, my_ierror, xarr, yarr, zarr )

   if (present(xout)) xout = xarr(1)
   if (present(yout)) yout = yarr(1)
   if (present(zout)) zout = zarr(1)

end subroutine spline_get_xyz_at_minz

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_dxyz_at_s_scalar( spl, sout, my_ierror, dxout, dyout, dzout)
!--function: evaluate spline derivatives d[x(s),y(s),z(s)]/ds of input grid (wheel/rail profile) to
!            output position sout
   implicit none
!--subroutine arguments
   type(t_spline),         intent(in)  :: spl
   real(kind=8),           intent(in)  :: sout
   integer,                intent(out) :: my_ierror
   real(kind=8), optional, intent(out) :: dxout, dyout, dzout
!--local variables
   integer      :: nout
   real(kind=8) :: sarr(1), dxarr(1), dyarr(1), dzarr(1)

   ! evaluate x-, y- and z-derivatives at output s-position

   nout    = 1
   sarr(1) = sout

   call spline_get_dxyz_at_s_arr( spl, nout, sarr, my_ierror, dxarr, dyarr, dzarr )

   if (present(dxout)) dxout = dxarr(1)
   if (present(dyout)) dyout = dyarr(1)
   if (present(dzout)) dzout = dzarr(1)

end subroutine spline_get_dxyz_at_s_scalar

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_dxyz_at_s_arr( spl, nout, sout, my_ierror, dxout, dyout, dzout )
!--function: evaluate spline derivatives d[x(s),y(s),z(s)]/ds of input grid (wheel/rail profile) at
!            output positions sout
   implicit none
!--subroutine arguments
   type(t_spline),         intent(in)  :: spl
   integer,                intent(in)  :: nout
   real(kind=8),           intent(in)  :: sout(nout)
   integer,                intent(out) :: my_ierror
   real(kind=8), optional, intent(out) :: dxout(nout), dyout(nout), dzout(nout)
!--local variables
   integer                   :: sub_ierror

   my_ierror = 0

   ! evaluate x-, y- and z-derivatives at output s-positions

   if (present(dxout) .and. spl%has_xdata) then
      call spline_eval( spl, ikXDIR, nout, sout, sub_ierror, f1_eval=dxout)
      if (my_ierror.eq.0) my_ierror = sub_ierror
   elseif (present(dxout)) then
      dxout(1:nout) = 0d0
   endif

   if (present(dyout)) then
      call spline_eval( spl, ikYDIR, nout, sout, sub_ierror, f1_eval=dyout)
      if (my_ierror.eq.0) my_ierror = sub_ierror
   endif

   if (present(dzout)) then
      call spline_eval( spl, ikZDIR, nout, sout, sub_ierror, f1_eval=dzout)
      if (my_ierror.eq.0) my_ierror = sub_ierror
   endif

end subroutine spline_get_dxyz_at_s_arr

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_dxyz_at_s_gfout( spl, gf_out, my_ierror )
!--function: compute derivatives d[xyz]/ds of spline [x(s),y(s),z(s)] (input grid) to output
!            grid-function gfout
   implicit none
!--subroutine arguments
   type(t_spline)               :: spl
   type(t_gridfnc3)             :: gf_out
   integer,         intent(out) :: my_ierror
!--local variables
   integer                  :: sub_ierror, nout
   type(t_grid), pointer    :: gout

   my_ierror = 0

   ! get output grid gout

   gout => gf_out%grid
   nout = gout%ntot

   ! evaluate derivatives at output s-positions

   if (spl%has_xdata) then
      if (ldebug.ge.3) call write_log(' spline_dxyz_at_s: spline_eval(dy_ds)')
      call spline_eval( spl, ikXDIR, nout, gout%s_prf, sub_ierror, f1_eval=gf_out%vx)
      if (my_ierror.eq.0) my_ierror = sub_ierror
   endif

   if (ldebug.ge.3) call write_log(' spline_dxyz_at_s: spline_eval(dy_ds)')
   call spline_eval( spl, ikYDIR, nout, gout%s_prf, sub_ierror, f1_eval=gf_out%vy)
   if (my_ierror.eq.0) my_ierror = sub_ierror

   if (ldebug.ge.3) call write_log(' spline_dxyz_at_s: spline_eval(dz_ds)')
   call spline_eval( spl, ikZDIR, nout, gout%s_prf, sub_ierror, f1_eval=gf_out%vn)
   if (my_ierror.eq.0) my_ierror = sub_ierror

end subroutine spline_get_dxyz_at_s_gfout

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_nvec_at_s_gfout( spl, gf_out, my_ierror )
!--function: compute outward normal vector for a plane curve [0, y, z] and store in grid-function gfout
   implicit none
!--subroutine arguments
   type(t_spline)               :: spl
   type(t_gridfnc3)             :: gf_out
   integer,         intent(out) :: my_ierror
!--local variables
   integer                  :: sub_ierror, iout, nout
   real(kind=8)             :: ty, tz, tnrm
   type(t_grid), pointer    :: gout

   my_ierror = 0

   ! get output grid gout

   gout => gf_out%grid
   nout = gout%ntot

   ! evaluate derivatives at output s-positions

   if (ldebug.ge.3) call write_log(' spline_nvec_at_s: spline_eval(dy_ds)')
   call spline_eval( spl, ikYDIR, nout, gout%s_prf, sub_ierror, f1_eval=gf_out%vy)
   if (my_ierror.eq.0) my_ierror = sub_ierror

   if (ldebug.ge.3) call write_log(' spline_nvec_at_s: spline_eval(dz_ds)')
   call spline_eval( spl, ikZDIR, nout, gout%s_prf, sub_ierror, f1_eval=gf_out%vn)
   if (my_ierror.eq.0) my_ierror = sub_ierror

   ! tangent vector: [ (dx/ds), dy/ds, dz/ds ], outward normal vector [ 0, dz/ds, -dy/ds ]

   do iout = 1, nout
      ty   = gf_out%vy(iout)
      tz   = gf_out%vn(iout)
      tnrm = sqrt(ty**2 + tz**2)
      gf_out%vx(iout) = 0d0
      gf_out%vy(iout) =  tz / tnrm
      gf_out%vn(iout) = -ty / tnrm
   enddo

end subroutine spline_get_nvec_at_s_gfout

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_pt_at_line_oyz(spl, pt, dir, s, x, y, z, my_ierror)
!--purpose: find the point s, y(s), z(s) in a spline curve that lies on the straight line formed by
!           point pt and direction dir, ignoring the x-direction
   implicit none
!--subroutine arguments:
   type(t_spline), intent(in)  :: spl
   type(t_vec),    intent(in)  :: pt, dir
   real(kind=8),   intent(out) :: s, x, y, z
   integer,        intent(out) :: my_ierror
!--local variables:
   integer,      parameter   :: maxit = 10
   real(kind=8), parameter   :: eps = 1d-4
   integer      :: i0, i1, k, sub_ierror
   real(kind=8) :: pt_x, pt_y, pt_z, dir_x, dir_y, dir_z, s0, s1, s_nwt, dst, dst0, dst1, ds,           &
                   dy_ds, dz_ds, ddst_ds, lmbd

   my_ierror = 0

   ! get components of direction vector [dir_x, dir_y, dir_z]^T

   dir_x = dir%x()
   dir_y = dir%y()
   dir_z = dir%z()
   dst = sqrt(dir_y**2 + dir_z**2)

   if (dst.lt.1d-10) then
      call write_log(' INTERNAL ERROR: spline_get_pt_at_line_oyz: dir is parallel to x-direction')
      x = 0d0
      y = 0d0
      z = 0d0
      return
   endif

   ! get components of point pt [pt_x, pt_y, pt_z]^T

   pt_x = pt%x()
   pt_y = pt%y()
   pt_z = pt%z()

   if (ldebug.ge.3) then
      write(bufout,'(6(a,f10.3),a)') ' get_pt_at_line: pt=(',pt_x,',',pt_y,',',pt_z,'), dir=(',dir_x,   &
                        ',',dir_y,',',dir_z,')'
      call write_log(1, bufout)
   endif

   ! create a bracket: compute distance of spline to line at s(1) and s(end)
   ! 2D distance to line: inner product with perpendicular vector

   i0   = idmin(spl%npnt, spl%ay0, 1)
   i1   = idmax(spl%npnt, spl%ay0, 1)
   s0   = spl%s(i0)
   s1   = spl%s(i1)

   if (ldebug.ge.3) then
      write(bufout,'(2(a,f10.3))') ' get_pt_at_line: get_xyz_at_s for miny, s0=',s0,' and maxy, s1=',s1
      call write_log(1, bufout)
   endif

   call spline_get_xyz_at_s( spl, s0, sub_ierror, x, y, z )
   if (my_ierror.eq.0) my_ierror = sub_ierror
   dst0 = (y - pt_y) * dir_z - (z - pt_z) * dir_y

   if (ldebug.ge.2) then
      write(bufout,'(a,i5, 4(a,f9.4))') ' get_pt_at_line: i_sta=', i0, ' :       s=',s0,',  yz=(',y,    &
                ',',z, '), dst=',dst0
      call write_log(1, bufout)
   endif

   call spline_get_xyz_at_s( spl, s1, sub_ierror, x, y, z )
   if (my_ierror.eq.0) my_ierror = sub_ierror
   dst1 = (y - pt_y) * dir_z - (z - pt_z) * dir_y

   if (ldebug.ge.2) then
      write(bufout,'(a,i5, 4(a,f9.4))') ' get_pt_at_line: i_end=', i1, ' :       s=',s1,',  yz=(',y,    &
                ',',z, '), dst=',dst1
      call write_log(1, bufout)
   endif

   if (dst0*dst1.gt.1d-10) then
      call write_log(' WARNING: spline_get_pt_at_line_oyz: s0 and s1 lie on same side of line')
      x = 0d0
      y = 0d0
      z = 0d0
      return
   endif

   ! set initial estimate in between s0 and s1

   k   = 0
   s   = s0 - dst0 / (dst1-dst0) * (s1 - s0)
   call spline_get_xyz_at_s( spl, s, sub_ierror, x, y, z )
   if (my_ierror.eq.0) my_ierror = sub_ierror
   dst = (y - pt_y) * dir_z - (z - pt_z) * dir_y

   if (ldebug.ge.2) then
      write(bufout,'(a,i3,6(a,f9.4))') ' iter k=',k,': s=',s,' in [',s0,',',s1,'], yz=(',y,',',z,       &
                '), dst=',dst
      call write_log(1, bufout)
   endif

   ! perform Newton iteration safe-guarded with bracket

   do while (abs(dst).gt.eps .and. k.le.maxit)

      k = k + 1

      ! update bracket

      if (dst*dst0.ge.0d0) then
         s0   = s
         dst0 = dst
      else
         s1   = s
         dst1 = dst
      endif

      ! get derivative d(dst) / ds

      call spline_get_dxyz_at_s( spl, s, dyout=dy_ds, dzout=dz_ds, my_ierror=sub_ierror )
      if (my_ierror.eq.0) my_ierror = sub_ierror
      ddst_ds = dy_ds * dir_z - dz_ds * dir_y

      ! compute Newton step

      ds    = -dst / ddst_ds
      s_nwt = s + ds

      ! revert to bisection when s_nwt falls outside bracket or lies close to its sides

      if (s_nwt.gt.s0 .and. s_nwt.lt.s1) then
         s = s_nwt
      else
         s = (s0 + s1) / 2d0
         if (ldebug.ge.2) then
            write(bufout,'(2(a,f9.4))') ' rejecting s_nwt=', s_nwt,', using bisection s=',s
            call write_log(1, bufout)
         endif
      endif

      call spline_get_xyz_at_s( spl, s, sub_ierror, x, y, z )
      if (my_ierror.eq.0) my_ierror = sub_ierror
      dst = (y - pt_y) * dir_z - (z - pt_z) * dir_y

      if (ldebug.ge.2) then
         write(bufout,'(a,i3,6(a,f9.4))') ' iter k=',k,': s=',s,' in [',s0,',',s1,'], yz=(',y,',',z,    &
                '), dst=',dst
         call write_log(1, bufout)
      endif

   enddo ! while

   if (k.ge.maxit .and. abs(dst).gt.eps) then
      write(bufout,'(a,i3,a,g12.4)') ' WARNING: spline_get_pt_at_line_oyz: no convergence, k=',k,       &
                ', dst=',dst
      call write_log(1, bufout)
      x = 0d0
      y = 0d0
      z = 0d0
      return
   endif

   ! get x-position on line pt + lmbd * dir

   lmbd = (y - pt_y) * dir_y + (z - pt_z) * dir_z
   x    = pt_x + lmbd * dir_x

end subroutine spline_get_pt_at_line_oyz

!------------------------------------------------------------------------------------------------------------

subroutine spline_test_data(itest, npnt, x, y, wgt, lambda)
!--function: set values {x, y} for testing the spline computation
   implicit none
!--subroutine arguments:
   integer,      intent(in)   :: itest
   integer,      intent(out)  :: npnt
   real(kind=8), intent(out)  :: lambda
   real(kind=8), dimension(:), allocatable :: x, y, wgt
!--local variables:
   integer                   :: i
   real(kind=8)              :: dx, a_int
   real(kind=8), dimension(:), allocatable :: a0, a1, a2, a3

   if (itest.eq.1) then

      ! define natural spline, generate data accordingly

      npnt = 4
      allocate(x(npnt), y(npnt), wgt(npnt), a3(npnt), a2(npnt), a1(npnt), a0(npnt))

      x (1:npnt)  = (/  0.,  1.,   2.,   3.  /)
      a3(1:npnt)  = (/    0.,    0.,   0.,   0.  /)
      a2(1:npnt)  = (/    0.,    0.,   0.,   0.  /)
      a1(1:npnt)  = (/    1.,    1.,   1.,   1.  /)
      a0(1:npnt)  = (/   10.,   11.,  12.,  13.  /)

      ! compute data values y(i)

      y(1:npnt) = a0(1:npnt)
      wgt(1:npnt) = 1d0

      ! set smoothing parameter: no smoothing

      lambda = 0d0

   elseif (itest.eq.2) then

      ! define natural spline, generate data accordingly

      npnt = 4
      allocate(x(npnt), y(npnt), wgt(npnt), a3(npnt), a2(npnt), a1(npnt), a0(npnt))

      x (1:npnt)  = (/  0., 0.5,   2.,   3.  /)
      a3(1:npnt)  = (/    0.,    0.,   0.,   0.  /)
      a2(1:npnt)  = (/    0.,    0.,   0.,   0.  /)
      a1(1:npnt)  = (/    1.,    1.,   1.,   1.  /)
      a0(1:npnt)  = (/   10.,  10.5,  12.,  13.  /)

      ! compute data values y(i)

      y(1:npnt) = a0(1:npnt)
      wgt(1:npnt) = 1d0

      ! set smoothing parameter: no smoothing

      lambda = 0d0

   elseif (itest.eq.3) then

      ! define natural spline, generate data accordingly

      npnt = 4
      allocate(x(npnt), y(npnt), wgt(npnt), a3(npnt), a2(npnt), a1(npnt), a0(npnt))

      x (1:npnt)  = (/  0.,  1.,   2.,   3.  /)
      a3(1:npnt)  = (/    1.,   -2.,   1.,   0.  /)
      a2(1:npnt)  = (/    0.,    3.,  -3.,   0.  /)
      a1(1:npnt)  = (/    0.,    3.,   3.,   0.  /)
      a0(1:npnt)  = (/   10.,   11.,  15.,  16.  /)

      ! compute data values y(i)

      y(1:npnt) = a0(1:npnt)
      wgt(1:npnt) = 1d0

      ! set smoothing parameter: no smoothing

      lambda = 0d0

   elseif (itest.eq.4) then

      ! meat consumption in USA, from Pollock1999

      npnt = 1941 - 1919 + 1
      allocate(x(npnt), y(npnt), wgt(npnt))

      do i = 1, npnt
         x(i) = 1919 - 1 + i
      enddo
      y(1:npnt) = (/ 171, 167, 164, 169, 179, 179, 172, 170, 169, 163, 162, 161, 159, 160, 165,         &
                     163, 145, 161, 157, 157, 166, 177, 181 /)

      wgt(1:npnt) = 1d0
      lambda = 1d0/3d0

   elseif (itest.eq.5) then

      ! construct natural spline

      npnt = 7
      allocate(x(npnt), y(npnt), wgt(npnt), a3(npnt), a2(npnt), a1(npnt), a0(npnt))

      x (1:npnt)  = (/  5.,  6.5,    7.,    8.,   10.,   14.,    15.  /)
      a3(1:npnt)  = (/   0.6,    -5.,  0.3,   1.4,  -0.7,   -999.,   0.  /)
      a2(1:npnt)  = (/    0.,  -999., -999., -999., -999.,  -999.,   0.  /)
      a1(1:npnt)  = (/    1.,  -999., -999., -999., -999.,  -999., -999. /)
      a0(1:npnt)  = (/    1.,  -999., -999., -999., -999.,  -999., -999. /)

      ! a(nseg) follows from boundary condition, natural spline: b(end) = \int a dx = 0

      a_int = 0d0
      do i = 1, npnt-2
         a_int = a_int + a3(i) * (x(i+1)-x(i))
      enddo
      a3(npnt-1) = -a_int / (x(npnt)-x(npnt-1))

      ! fill spline parameters recursively, at endpoints of segments 1 .. nseg --> points 2 .. npnt

      do i = 2, npnt
         dx = x(i) - x(i-1)
         a0(i) =   a3(i-1)*dx**3 +   a2(i-1)*dx**2 + a1(i-1)*dx + a0(i-1)
         a1(i) = 3*a3(i-1)*dx**2 + 2*a2(i-1)*dx    + a1(i-1)
         a2(i) = 3*a3(i-1)*dx    +   a2(i-1)
      enddo

      ! compute data values y(i)

      y(1:npnt) = a0(1:npnt)

      ! set weights w(i)

      wgt(1:npnt) = 1d0

      ! set smoothing parameter

      lambda = 0.25d0

   else

      write(*,*) 'Test number out of range: itest=',itest
      call abort_run()

   endif

end subroutine spline_test_data

!------------------------------------------------------------------------------------------------------------

subroutine spline_test(itest)
   implicit none
!--subroutine arguments
   integer              :: itest
!--local variables:
   integer              :: npnt, neval, nkink, ikinks(99), i, jseg, my_ierror, sub_ierror
   real(kind=8)         :: lambda
   real(kind=8), dimension(:), allocatable :: x, y, wgt, a3, a2, a1, a0, xeval, yeval
!--external functions used:

   my_ierror = 0
   write(*,*) 'Obtained itest=',itest

   if (itest.le.0) then

      ! test 0: testing locate_segment

      npnt = 5
      neval = 10
      allocate(x(npnt+1), xeval(neval))
      x     = (/ 1d0, 2d0, 3d0, 4d0, 5d0, 5d0 /)
      xeval = (/ 3.5d0, 1.5d0, 1d0-1d-6, 5d0, 5d0+1d-6, 3d0, 3d0-1d-6, 3d0+1d-6, 5d0-1d-6, 1d0 /)

      jseg = 0
      do i = 1, neval
         ! find the segment iseg [ x(j), x(j+1) ] containing xeval
         call locate_segment( npnt, x, xeval(i), jseg )
         if (jseg.le.0) then
            write(*,'(a,i2,2(a,f12.6))') 'point xev(',i,')=',xeval(i),' lies before x-data, x(1)=',x(1)
         elseif (jseg.gt.npnt) then
            write(*,'(a,i2,2(a,f12.6))') 'point xev(',i,')=',xeval(i),' lies  after x-data, x(n)=',x(npnt)
         else
            write(*,'(a,i2,a,f12.6,a,i2,2(a,f12.6),a)') 'point xev(',i,')=',xeval(i),' lies in jseg=',     &
               jseg,': [',x(jseg),',',x(jseg+1),']'
         endif
      enddo

      deallocate(x, xeval)

   else

      ! Set {x, y} according to the test selected

      call spline_test_data(itest, npnt, x, y, wgt, lambda)

      write(*,*) 'npnt = ',npnt
      do i = 1, npnt
         write(*,'(a,i3,a,2f10.3,a,f10.3)') ' ip=', i,': x,y=',x(i), y(i),', w=',wgt(i)
      enddo

      ! Compute {a, b, c, d} for smoothing spline

      allocate(a3(npnt), a2(npnt), a1(npnt), a0(npnt))

      nkink = 0
      call ppspline_make_simple_kink(npnt, x, y, lambda, a3, a2, a1, a0, nkink, ikinks, sub_ierror, wgt)
      if (my_ierror.eq.0) my_ierror = sub_ierror

      write(*,*) 'nseg = ',npnt-1
      do i = 1, npnt
         write(*,'(a,i3,a,4f10.3)') ' is=', i,': a3-a0=', a3(i), a2(i), a1(i), a0(i)
      enddo

      ! evaluate spline at list of x-values

      neval = 6
      allocate(xeval(neval), yeval(neval))

      xeval = (/ -1d-10, 1d-10, 1d0, 2.5d0, 2.99d0, 3.1d0 /)
      call spline_eval(npnt, x, 'y', a3, a2, a1, a0, neval, xeval, f_eval=yeval, ierror=sub_ierror)
      if (my_ierror.eq.0) my_ierror = sub_ierror

      write(*,*) 'neval = ',neval
      do i = 1, neval
         write(*,'(a,i3,a,f10.3,a,f10.3)') ' ie=', i,': x=',xeval(i), ', y=',yeval(i)
      enddo

      xeval = (/ 3.1d0, 2.99d0, 2.5d0, 1d0, 1d-10, -1d-10 /)
      call spline_eval(npnt, x, 'y', a3, a2, a1, a0, neval, xeval, f_eval=yeval, ierror=sub_ierror)
      if (my_ierror.eq.0) my_ierror = sub_ierror

      write(*,*) 'neval = ',neval
      do i = 1, neval
         write(*,'(a,i3,a,f10.3,a,f10.3)') ' ie=', i,': x=',xeval(i), ', y=',yeval(i)
      enddo

   endif

end subroutine spline_test

!------------------------------------------------------------------------------------------------------------

end module m_spline_get
