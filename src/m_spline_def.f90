!------------------------------------------------------------------------------------------------------------
! m_spline_def - definitions and basic handling for 1D spline in PP-form (piecewise polynomial)
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_spline_def
   use m_globals
   use m_markers
   use m_ptrarray
   use m_interp
   implicit none
   private

   ! Debugging for module m_spline_def

   public  splinedef_set_debug

   integer  :: ldebug    =  0    ! local level of debugging
   integer  :: ii_debug  = -1    ! output point for which detailed info is requested (-1 = none)
   integer  :: iel_debug = -1    ! input element for which detailed info is requested (-1 = none)

   ! Data type for splines:

   public t_spline

   ! Functions defined on splines:

   public  spline_nullify
   public  spline_allocate
   public  spline_initx
   public  spline_copy
   public  spline_trim
   public  spline_print
   public  spline_destroy

   public  spline_eval
   public  spline_eval_arr
   public  spline_eval_spl

   interface spline_eval
      module procedure spline_eval_arr
      module procedure spline_eval_spl
   end interface spline_eval

   public  solve_cubic_segm
   public  solve_cubic_newton
   public  solve_cubic_cardano
   public  spline_get_s_at_f
   public  spline_get_s_at_f_spl

   interface spline_get_s_at_f
      module procedure spline_get_s_at_f_spl
   end interface spline_get_s_at_f

   public  locate_one_extremum
   public  locate_extremal_values

   public spline_shift
   public spline_mirror_y
   public spline_rotate
   public spline_roll
   public spline_2glob
   public spline_2glob_or
   public spline_2glob_m
   public spline_2loc
   public spline_2loc_m

   interface spline_2glob
      module procedure spline_2glob_m
      module procedure spline_2glob_or
   end interface spline_2glob

   interface spline_2loc
      module procedure spline_2loc_m
   end interface spline_2loc

!------------------------------------------------------------------------------------------------------------
!  data for a parametric spline (s, [x(s)], y(s), z(s)):

   type :: t_spline
      integer      :: npnt
      logical      :: has_xdata

      real(kind=8), dimension(:),   pointer  :: s     => NULL() ! (npnt)

      real(kind=8), dimension(:,:), pointer  :: axspl => NULL() ! (npnt,4)
      real(kind=8), dimension(:),   pointer  :: ax0   => NULL()
      real(kind=8), dimension(:),   pointer  :: ax1   => NULL()
      real(kind=8), dimension(:),   pointer  :: ax2   => NULL()
      real(kind=8), dimension(:),   pointer  :: ax3   => NULL()

      real(kind=8), dimension(:,:), pointer  :: ayspl => NULL() ! (npnt,4)
      real(kind=8), dimension(:),   pointer  :: ay0   => NULL()
      real(kind=8), dimension(:),   pointer  :: ay1   => NULL()
      real(kind=8), dimension(:),   pointer  :: ay2   => NULL()
      real(kind=8), dimension(:),   pointer  :: ay3   => NULL()

      real(kind=8), dimension(:,:), pointer  :: azspl => NULL() ! (npnt,4)
      real(kind=8), dimension(:),   pointer  :: az0   => NULL()
      real(kind=8), dimension(:),   pointer  :: az1   => NULL()
      real(kind=8), dimension(:),   pointer  :: az2   => NULL()
      real(kind=8), dimension(:),   pointer  :: az3   => NULL()

      integer      :: nsec_uniy
      integer,      dimension(:),   pointer  :: ipnt_uniy => NULL() ! (nsec_uni+1)

      integer      :: nsec_top
      real(kind=8), dimension(:),   pointer  :: ysec_top  => NULL() ! (nsec_top+1)
      integer,      dimension(:),   pointer  :: iuni_top  => NULL() ! (nsec_top)

      ! npnt            number of points in spline
      ! has_xdata       indicating that axspl is present, esp. for contact locus
      !
      ! s       [mm]    arc length parameter
      ! axspl, ax0-ax3  simple spline coefficients for x(s), with ax0-3 pointers to columns of axspl
      ! ayspl, ay0-ay3  simple spline coefficients for y(s), with ay0-3 pointers to columns of ayspl
      ! azspl, az0-az3  simple spline coefficients for z(s), with az0-3 pointers to columns of azspl
      !
      ! nsec_uniy       number of sections with monotonous y(s), invertible, uni-valued s(y)
      ! ipnt_uniy       start-index of each section with uni-valued y(s)
      !
      ! nsec_top        number of sections [y_j, y_{j+1}] used in 'top view'
      ! ysec_top        start-position y_j of each section used in top view
      ! iuni_top        iuni_top(j) = index of uni-valued section seen in top view section j
      !
      ! Note: uni-sections and top-sections are volatile -- discarded upon shift/rotate/trimming

   end type t_spline

contains

!------------------------------------------------------------------------------------------------------------

subroutine splinedef_set_debug(new_ldebug, new_ii_debug, new_iel_debug)
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
      write(bufout,'(a,i3,2(a,i7))') ' spline-def:  debugging level =',ldebug,', ii_debug =', ii_debug, &
                ', iel_debug =', iel_debug
      call write_log(1, bufout)
   endif

end subroutine splinedef_set_debug

!------------------------------------------------------------------------------------------------------------

subroutine spline_nullify(spl)
!--purpose: initialize spline structure, nullify pointers
   implicit none
!--subroutine parameters:
   type(t_spline)            :: spl

   spl%npnt      =  0
   spl%has_xdata = .false.
   spl%nsec_uniy =  0
   spl%nsec_top  =  0

   spl%s     => NULL()

   spl%axspl => NULL()
   spl%ax0   => NULL()
   spl%ax1   => NULL()
   spl%ax2   => NULL()
   spl%ax3   => NULL()

   spl%ayspl => NULL()
   spl%ay0   => NULL()
   spl%ay1   => NULL()
   spl%ay2   => NULL()
   spl%ay3   => NULL()

   spl%azspl => NULL()
   spl%az0   => NULL()
   spl%az1   => NULL()
   spl%az2   => NULL()
   spl%az3   => NULL()

   spl%ipnt_uniy => NULL()
   spl%ysec_top  => NULL()
   spl%iuni_top  => NULL()

end subroutine spline_nullify

!------------------------------------------------------------------------------------------------------------

subroutine spline_allocate(spl, npnt, nsec_uniy, nsec_top)
!--purpose: (re-)allocate the arrays needed for spline data
   implicit none
!--subroutine parameters:
   type(t_spline)       :: spl
   integer              :: npnt
   integer, optional    :: nsec_uniy, nsec_top

   spl%npnt = npnt
   if (present(nsec_uniy)) then
      spl%nsec_uniy = nsec_uniy
   else
      spl%nsec_uniy = 0
   endif
   if (present(nsec_top)) then
      spl%nsec_top  = nsec_top 
   else
      spl%nsec_top  = 0
   endif

   ! allocate arc-length array s(n)

   call reallocate_arr(spl%s, npnt)

   ! allocate simple spline arrays axspl(n,4), ayspl(n,4) and azspl(n,4),
   !     storing {ax0-ax4} for x(s), {ay0-ay3} for y(s) and {az0-az3} for z(s)
   ! set/update pointers

   if (spl%has_xdata) then
      call reallocate_arr(spl%axspl, spl%npnt, 4)
      spl%ax0 => spl%axspl(:,1)
      spl%ax1 => spl%axspl(:,2)
      spl%ax2 => spl%axspl(:,3)
      spl%ax3 => spl%axspl(:,4)
   endif

   call reallocate_arr(spl%ayspl, spl%npnt, 4)
   spl%ay0 => spl%ayspl(:,1)
   spl%ay1 => spl%ayspl(:,2)
   spl%ay2 => spl%ayspl(:,3)
   spl%ay3 => spl%ayspl(:,4)

   call reallocate_arr(spl%azspl, spl%npnt, 4)
   spl%az0 => spl%azspl(:,1)
   spl%az1 => spl%azspl(:,2)
   spl%az2 => spl%azspl(:,3)
   spl%az3 => spl%azspl(:,4)

   ! allocate description for uni-valued sections

   call reallocate_arr(spl%ipnt_uniy, spl%nsec_uniy+1)

   ! allocate description for top-view sections

   call reallocate_arr(spl%ysec_top, spl%nsec_top+1)
   call reallocate_arr(spl%iuni_top, spl%nsec_top)

end subroutine spline_allocate

!------------------------------------------------------------------------------------------------------------

subroutine spline_initx(spl)
!--purpose: initialize the arrays needed for x-spline data, set x==0
   implicit none
!--subroutine parameters:
   type(t_spline)            :: spl

   ! allocate simple spline array axspl(n,4)

   if (spl%npnt.le.0) then
      write(bufout,*) 'INTERNAL ERROR (spline_initx): npnt=', spl%npnt
      call write_log(1, bufout)
   endif

   ! write(bufout,*) 'add x-data to spline with npnt=',spl%npnt,' points'
   ! call write_log(1, bufout)

   call reallocate_arr(spl%axspl, spl%npnt, 4)

   ! set/update pointers

   spl%ax0 => spl%axspl(:,1)
   spl%ax1 => spl%axspl(:,2)
   spl%ax2 => spl%axspl(:,3)
   spl%ax3 => spl%axspl(:,4)

   ! set ax0--ax3 == 0 for x == 0

   spl%ax0(1:spl%npnt) = 0d0
   spl%ax1(1:spl%npnt) = 0d0
   spl%ax2(1:spl%npnt) = 0d0
   spl%ax3(1:spl%npnt) = 0d0

   ! assert that spline may go outside oyz

   spl%has_xdata = .true.

end subroutine spline_initx

!------------------------------------------------------------------------------------------------------------

subroutine spline_copy(spl_in, spl_out)
!--purpose: copy input spline spl_in to output spline spl_out
   implicit none
!--subroutine parameters:
   type(t_spline)        :: spl_in, spl_out
!--local variables:
   integer      :: ispl

   spl_out%npnt        = spl_in%npnt
   spl_out%has_xdata   = spl_in%has_xdata
   spl_out%nsec_uniy   = spl_in%nsec_uniy
   spl_out%nsec_top    = spl_in%nsec_top

   ! if s-array is present in input: allocate, copy s-coordinates, else: clear s-coordinates

   if (.not.associated(spl_in%s)) then
      call write_log('INTERNAL ERROR(spline_copy): spl_in has not been filled.')
      call abort_run()
   endif

   associate(    npnt => spl_in%npnt, nsec_uni => spl_in%nsec_uniy, nsec_top => spl_in%nsec_top)

   ! allocate output-arrays at size needed

   call spline_allocate(spl_out, npnt, nsec_uni, nsec_top)

   ! copy data

   spl_out%s(1:npnt) = spl_in%s(1:npnt)

   do ispl = 1, 4
      if (spl_in%has_xdata) spl_out%axspl(1:npnt,ispl) = spl_in%axspl(1:npnt,ispl)
      spl_out%ayspl(1:npnt,ispl) = spl_in%ayspl(1:npnt,ispl)
      spl_out%azspl(1:npnt,ispl) = spl_in%azspl(1:npnt,ispl)
   enddo

   if (nsec_uni.gt.0) then
      spl_out%ipnt_uniy(1:nsec_uni+1) = spl_in%ipnt_uniy(1:nsec_uni+1)
   endif

   if (nsec_top.gt.0) then
      spl_out%ysec_top(1:nsec_top+1) = spl_in%ysec_top(1:nsec_top+1)
      spl_out%iuni_top(1:nsec_top)   = spl_in%iuni_top(1:nsec_top)
   endif

   end associate

end subroutine spline_copy

!------------------------------------------------------------------------------------------------------------

subroutine spline_trim(spl_in, spl_out, ilow_arg, ihig_arg, s_low, s_hig)
!--purpose: create trimmed version of spline in output grid gout
   implicit none
!--subroutine parameters:
   type(t_spline)         :: spl_in, spl_out
   integer,      optional :: ilow_arg, ihig_arg
   real(kind=8), optional :: s_low, s_hig
!--local variables:
   integer      :: ilow, ihig, ispl

   if (present(ilow_arg) .and. present(ihig_arg)) then

      ilow = ilow_arg
      ihig = ihig_arg

   elseif (present(s_low) .and. present(s_hig)) then

      ! determine the indices to keep

      call locate_interval(spl_in%npnt, spl_in%s, s_low, s_hig, ilow, ihig)

      if (ldebug.ge.1) then
         write(bufout,'(2(a,f8.3),2(a,i4),a)') 'spline_trim: s=[',s_low,',',s_hig,']: keeping ip=[',  &
                ilow, ',',ihig,']'
         call write_log(1, bufout)
      endif
   else

      call write_log('INTERNAL ERROR(spline_trim): either ilow or s_low must be given.')
      if (present(ilow_arg)) call write_log('ilow_arg is present')
      if (present(ihig_arg)) call write_log('ihig_arg is present')
      if (present(s_low)) call write_log('s_low is present')
      if (present(s_hig)) call write_log('s_hig is present')
      call abort_run()

   endif

   ! set the number of points in the output spline

   spl_out%npnt = ihig - ilow + 1
   call spline_allocate(spl_out, spl_out%npnt)

   spl_out%s(1:spl_out%npnt) = spl_in%s(ilow:ihig)

   do ispl = 1, 4
      if (spl_in%has_xdata) spl_out%axspl(1:spl_out%npnt,ispl) = spl_in%axspl(ilow:ihig,ispl)
      spl_out%ayspl(1:spl_out%npnt,ispl) = spl_in%ayspl(ilow:ihig,ispl)
      spl_out%azspl(1:spl_out%npnt,ispl) = spl_in%azspl(ilow:ihig,ispl)
   enddo

   ! clear administration of uni-valued sections and top-view sections

   spl_out%nsec_uniy = 0
   spl_out%nsec_top  = 0

   if (associated(spl_out%ipnt_uniy)) deallocate(spl_out%ipnt_uniy)
   if (associated(spl_out%ysec_top))  deallocate(spl_out%ysec_top)
   if (associated(spl_out%iuni_top))  deallocate(spl_out%iuni_top)
   spl_out%ipnt_uniy => NULL()
   spl_out%ysec_top  => NULL()
   spl_out%iuni_top  => NULL()

end subroutine spline_trim

!------------------------------------------------------------------------------------------------------------

subroutine spline_print(spl, nam, idebug, ndigit)
!--function: print information on spline spl
   implicit none
!--subroutine arguments
   type(t_spline)       :: spl
   character(len=*)     :: nam
   integer              :: idebug
   integer, optional    :: ndigit       ! number of significant digits
!--local variables
   character(len= 1), parameter  :: cnams(1:3) = (/ 'x', 'y', 'z' /)
   character(len=30), parameter  :: spaces = '                              '
   integer              :: my_ndigit, my_len
   integer              :: ii, j, nline, iuni, itop
   character(len=18)    :: strng(5)

   if (present(ndigit)) then
      my_ndigit = ndigit
   else
      my_ndigit = 6
   endif
   my_ndigit = max(2, min(10, my_ndigit))
   my_len    = 8 + my_ndigit

   if (.not.associated(spl%s) .or. .not.associated(spl%ayspl)) then
      write(bufout,'(3a)') ' spline ',trim(nam),' has not yet been defined'
      call write_log(1, bufout)
      return
   endif

   ! idebug>=6: array sizes

   if (idebug.ge.6) then
      call print_array_size_1d(spl%s, 's')
      call print_array_size_2d(spl%axspl, 'axspl')
      call print_array_size_1d(spl%ax0, 'ax0')
      call print_array_size_1d(spl%ax1, 'ax1')
      call print_array_size_1d(spl%ax2, 'ax2')
      call print_array_size_1d(spl%ax3, 'ax3')
      call print_array_size_2d(spl%ayspl, 'ayspl')
      call print_array_size_1d(spl%ay0, 'ay0')
      call print_array_size_1d(spl%ay1, 'ay1')
      call print_array_size_1d(spl%ay2, 'ay2')
      call print_array_size_1d(spl%ay3, 'ay3')
      call print_array_size_2d(spl%azspl, 'azspl')
      call print_array_size_1d(spl%az0, 'az0')
      call print_array_size_1d(spl%az1, 'az1')
      call print_array_size_1d(spl%az2, 'az2')
      call print_array_size_1d(spl%az3, 'az3')
   endif

   ! idebug>=3: information on uni-valued y-sections

   if (idebug.ge.3) then
      if (spl%nsec_uniy.le.0) then
         call write_log(' spline does not have information on uni-valued y-sections')
      else
         write(bufout,'(a,i3,2a)') ' spline has', spl%nsec_uniy,' sections with monotonous y(s) / ',    &
                'uni-valued s(y)'
         call write_log(1, bufout)

         nline = spl%nsec_uniy / 10 + 1

         write(bufout, 211) '         iuni  ', (iuni, iuni=1, min(100,spl%nsec_uniy))
         call write_log(nline, bufout)
 211     format(a, 10(i9,:,','), :, /, 20(15x, 10(i9,:,','), :, /))

         write(bufout, 211) '    ipnt_uniy =', (spl%ipnt_uniy(iuni), iuni=1, min(100,spl%nsec_uniy+1))
         call write_log(nline, bufout)

         write(bufout, 212) '    ssec_uniy =', (spl%s(spl%ipnt_uniy(iuni)), iuni=1, min(100,spl%nsec_uniy+1))
 212     format(a, 10(f9.3,:,','), :, /, 20(15x, 10(f9.3,:,','), :, /))
         call write_log(nline, bufout)

         write(bufout, 212) '    ysec_uniy =',(spl%ay0(spl%ipnt_uniy(iuni)), iuni=1,min(100,spl%nsec_uniy+1))
         call write_log(nline, bufout)

         write(bufout, 212) '    zsec_uniy =',(spl%az0(spl%ipnt_uniy(iuni)), iuni=1,min(100,spl%nsec_uniy+1))
         call write_log(nline, bufout)
      endif
   endif

   ! idebug>=3: information on visible parts in top view

   if (idebug.ge.3) then
      if (spl%nsec_top.le.0) then
         call write_log(' spline does not have information on visible parts in top view')
      else
         write(bufout,'(a,i4,a)') ' spline top view consists of',spl%nsec_top,' parts [y_j, y_{j+1}]'
         call write_log(1, bufout)

         nline = spl%nsec_top / 10 + 1

         write(bufout, 211) '         itop  ', (itop, itop=1, min(100,spl%nsec_top))
         call write_log(nline, bufout)

         write(bufout, 212) '     ysec_top =', (spl%ysec_top(itop), itop=1, min(100,spl%nsec_top+1))
         call write_log(nline, bufout)

         write(bufout, 211) '     iuni_top =', (spl%iuni_top(itop), itop=1, min(100,spl%nsec_top))
         call write_log(nline, bufout)
      endif
   endif

   ! idebug>=5: full details w.r.t. spline data

   if (idebug.ge.5) then
      if (spl%has_xdata) then

         write(bufout,'(3a)') 'splx_', trim(nam), ' = [     % spline x(s):'
         call write_log(1, bufout)

         write(bufout,'(12a)') ' %       ii', spaces(1:my_len-5), 'si', spaces(1:my_len-2), 'ax0',     &
                spaces(1:my_len-2), 'ax1', spaces(1:my_len-2), 'ax2', spaces(1:my_len-2), 'ax3'
         call write_log(1, bufout)

         do ii = 1, spl%npnt
            strng(1) = fmt_gs(my_len, my_ndigit, spl%s(ii))
            do j = 1, 4
               strng(1+j) = fmt_gs(my_len, my_ndigit, spl%axspl(ii,j))
            enddo

            write(bufout,311) ii, (strng(j)(1:my_len), j=1,5)
 311        format(' ii=',i6,':',5(a,1x))
            call write_log(1, bufout)
         enddo
         call write_log('];')
      endif

      write(bufout,'(3a)') 'sply_', trim(nam), ' = [     % spline y(s):'
      call write_log(1, bufout)

      write(bufout,'(12a)') ' %       ii', spaces(1:my_len-5), 'si', spaces(1:my_len-2), 'ay0',         &
               spaces(1:my_len-2), 'ay1', spaces(1:my_len-2), 'ay2', spaces(1:my_len-2), 'ay3'
      call write_log(1, bufout)

      do ii = 1, spl%npnt
         strng(1) = fmt_gs(my_len, my_ndigit, spl%s(ii))
         do j = 1, 4
            strng(1+j) = fmt_gs(my_len, my_ndigit, spl%ayspl(ii,j))
         enddo
         write(bufout,311) ii, (strng(j)(1:my_len), j=1,5)
         call write_log(1, bufout)
      enddo
      call write_log('];')

      write(bufout,'(3a)') 'splz_', trim(nam), ' = [     % spline z(s):'
      call write_log(1, bufout)

      write(bufout,'(12a)') ' %       ii', spaces(1:my_len-5), 'si', spaces(1:my_len-2), 'az0',         &
               spaces(1:my_len-2), 'az1', spaces(1:my_len-2), 'az2', spaces(1:my_len-2), 'az3'
      call write_log(1, bufout)

      do ii = 1, spl%npnt
         strng(1) = fmt_gs(my_len, my_ndigit, spl%s(ii))
         do j = 1, 4
            strng(1+j) = fmt_gs(my_len, my_ndigit, spl%azspl(ii,j))
         enddo
         write(bufout,311) ii, (strng(j)(1:my_len), j=1,5)
         call write_log(1, bufout)
      enddo
      call write_log('];')

   endif ! idebug >= 5

end subroutine spline_print

!------------------------------------------------------------------------------------------------------------

subroutine spline_destroy(spl)
!--purpose: clean-up allocated arrays for splines, nullify pointers
   implicit none
!--subroutine parameters:
   type(t_spline)  :: spl

   if (associated(spl%s))         deallocate(spl%s)
   if (associated(spl%axspl))     deallocate(spl%axspl)
   if (associated(spl%ayspl))     deallocate(spl%ayspl)
   if (associated(spl%azspl))     deallocate(spl%azspl)
   if (associated(spl%ipnt_uniy)) deallocate(spl%ipnt_uniy)
   if (associated(spl%ysec_top))  deallocate(spl%ysec_top)
   if (associated(spl%iuni_top))  deallocate(spl%iuni_top)
   
   call spline_nullify(spl)

end subroutine spline_destroy

!------------------------------------------------------------------------------------------------------------

subroutine spline_eval_arr(npnt, s_spl, nam_f, a3, a2, a1, a0, nout, s_eval, ierror, exterval,          &
                            f_eval, f1_eval, f2_eval, f3_eval)
!--function: evaluate spline {s, a,b,c,d} at points s_eval, computing f_eval and its derivatives
!             f1=f', f2=f'', f3=f'''
   implicit none
!--subroutine arguments:
   integer,          intent(in)   :: npnt, nout
   character(len=1), intent(in)   :: nam_f
   real(kind=8),     intent(in)   :: s_spl(npnt), a3(npnt), a2(npnt), a1(npnt), a0(npnt), s_eval(nout)
   integer,          intent(out)  :: ierror
   real(kind=8),     intent(in),  optional :: exterval
   real(kind=8),     intent(out), optional :: f_eval(nout), f1_eval(nout), f2_eval(nout), f3_eval(nout)
!--local variables:
   integer              :: iseg, iout, nseg
   logical              :: has_exter, eval_f(0:3)
   real(kind=8)         :: sloc
   character(len=160)   :: strng

   ierror = 0

   ! determine presence of optional outputs

   has_exter = (present(exterval))
   eval_f(0) = (present(f_eval))
   eval_f(1) = (present(f1_eval))
   eval_f(2) = (present(f2_eval))
   eval_f(3) = (present(f3_eval))

   if (ldebug.ge.4) then
      write(bufout,'(3a,4l2)') ' eval: starting spline_eval for f="',nam_f(1:1),'" with eval_f=',eval_f(0:3)
      call write_log(1, bufout)
   endif

   nseg = npnt - 1

   if (ldebug.ge.4) then
      write(bufout,'(a,i6,3a,i6,a)') ' eval: input spline has',nseg,' segments, evaluating ',nam_f(1:1), &
                ' at',nout,' s-positions'
      call write_log(1, bufout)
   endif

   ! s_spl is assumed to be in ascending order; s_eval may be unsorted

   do iout = 1, nout

      ! locate segment iseg with s_spl(iseg) <= s_eval(iout) <= s_spl(iseg+1)
      ! using s_spl(0) = -\infty, s_spl(nseg+1) = \infty

      call locate_segment( npnt, s_spl, s_eval(iout), iseg )

      if (ldebug.ge.4) then
         write(bufout,'(a,i4,a,f10.3,a,i4,2(a,f12.6),a)') ' eval: s_eval(',iout,')=',s_eval(iout),      &
                ' lies in iseg=',iseg,', s=[',s_spl(max(1,iseg)),',',s_spl(min(npnt,iseg)),']'
         call write_log(1, bufout)
      endif

      if (iseg.le.0) then

         ! s_eval < s_spl(1) : constant extrapolation

         if (has_exter) then
            if (eval_f(0)) f_eval(iout)  = exterval
            if (eval_f(1)) f1_eval(iout) = exterval
            if (eval_f(2)) f2_eval(iout) = exterval
            if (eval_f(3)) f3_eval(iout) = exterval
         else
            if (eval_f(0)) f_eval(iout)  = a0(1)
            if (eval_f(1)) f1_eval(iout) = a1(1)
            if (eval_f(2)) f2_eval(iout) = 0d0
            if (eval_f(3)) f3_eval(iout) = 0d0
         endif

         if (ldebug.ge.3) then
            write(strng,'(a,l2)') ', has_ext=', has_exter
            if (eval_f(0)) write(strng,'(a,a,f12.6)') trim(strng), ', f=',  f_eval(iout)
            if (eval_f(1)) write(strng,'(a,a,f12.6)') trim(strng), ', f1=', f1_eval(iout)
            if (eval_f(2)) write(strng,'(a,a,f12.6)') trim(strng), ', f2=', f2_eval(iout)
            if (eval_f(3)) write(strng,'(a,a,f12.6)') trim(strng), ', f3=', f3_eval(iout)
            write(bufout,'(a,i4,a,f12.6,2a)') ' eval: iout=',iout,' has s=',s_eval(iout),               &
                ', before start', trim(strng)
            call write_log(1, bufout)
         endif

      elseif (iseg.le.nseg) then

         ! s_spl(1) <= s_eval < s_spl(nseg+1) : cubic polynomials

         sloc         = s_eval(iout) - s_spl(iseg)
         if (eval_f(0)) f_eval(iout)  =  ((    a3(iseg) * sloc  +     a2(iseg)) * sloc +                &
                                                                      a1(iseg)) * sloc + a0(iseg)
         if (eval_f(1)) f1_eval(iout) =   (3d0*a3(iseg) * sloc  + 2d0*a2(iseg)) * sloc + a1(iseg)
         if (eval_f(2)) f2_eval(iout) =    6d0*a3(iseg) * sloc  + 2d0*a2(iseg)
         if (eval_f(3)) f3_eval(iout) =    6d0*a3(iseg)

         if (ldebug.ge.3) then
            strng = ' '
            if (eval_f(0)) write(strng,'(  3a,f12.6)')              ', ',nam_f(1:1),'=',  f_eval(iout)
            if (eval_f(1)) write(strng,'(a,3a,f10.6)') trim(strng), ', ',nam_f(1:1),'1=', f1_eval(iout)
            if (eval_f(2)) write(strng,'(a,3a,f10.6)') trim(strng), ', ',nam_f(1:1),'2=', f2_eval(iout)
            if (eval_f(3)) write(strng,'(a,3a,f10.6)') trim(strng), ', ',nam_f(1:1),'3=', f3_eval(iout)
            write(bufout,'(a,i4,a,f12.6,a,i4,a)') ' eval: iout=',iout,' has s=',s_eval(iout),           &
                     ', in seg', iseg, trim(strng)
            call write_log(1, bufout)
         endif

      else

         ! s_eval >= s_spl(nseg+1): constant extrapolation

         if (has_exter) then
            if (eval_f(0)) f_eval(iout)  = exterval
            if (eval_f(1)) f1_eval(iout) = exterval
            if (eval_f(2)) f2_eval(iout) = exterval
            if (eval_f(3)) f3_eval(iout) = exterval
         else
            if (eval_f(0)) f_eval(iout)  = a0(npnt)
            if (eval_f(1)) f1_eval(iout) = a1(npnt)
            if (eval_f(2)) f2_eval(iout) = 0d0
            if (eval_f(3)) f3_eval(iout) = 0d0
         endif

         if (ldebug.ge.3) then
            strng = ' '
            if (eval_f(0)) write(strng,'(  a,f12.6)')              ', f=',  f_eval(iout)
            if (eval_f(1)) write(strng,'(a,a,f12.6)') trim(strng), ', f1=', f1_eval(iout)
            if (eval_f(2)) write(strng,'(a,a,f12.6)') trim(strng), ', f2=', f2_eval(iout)
            if (eval_f(3)) write(strng,'(a,a,f12.6)') trim(strng), ', f3=', f3_eval(iout)
            write(bufout,'(a,i4,a,f12.6,2a)') ' eval: iout=',iout,' has s=',s_eval(iout),               &
                        ', after end   ', trim(strng)
            call write_log(1, bufout)
         endif

      endif ! seval before start, in [s(1), s(npnt)], after end of spline

   enddo ! iout = 1, nout

end subroutine spline_eval_arr

!------------------------------------------------------------------------------------------------------------

subroutine spline_eval_spl(spl, ikarg, nout, s_eval, ierror, exterval, f_eval, f1_eval, f2_eval, f3_eval)
!--function: evaluate parametric spline 'spl' in direction ikarg at points s_eval, computing f_eval and
!            its derivatives f1=f', f2=f'', f3=f'''
   implicit none
!--subroutine arguments:
   type(t_spline)                 :: spl
   integer,          intent(in)   :: ikarg, nout
   real(kind=8),     intent(in)   :: s_eval(nout)
   integer,          intent(out)  :: ierror
   real(kind=8),     intent(in),  optional :: exterval
   real(kind=8),     intent(out), optional :: f_eval(nout), f1_eval(nout), f2_eval(nout), f3_eval(nout)
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
      call write_log('INTERNAL ERROR (spline_eval_spl): ikarg invalid.')
      call abort_run()
   endif

   if (associated(coef)) then
      call spline_eval_arr( spl%npnt, spl%s, nam_f, coef(:,4), coef(:,3), coef(:,2), coef(:,1),   &
                            nout, s_eval, ierror, exterval, f_eval, f1_eval, f2_eval, f3_eval)
   else
      if (.false.) call write_log('spline_eval_spl: coef not associated')
      if (present(f_eval))  f_eval(1:nout)  = 0d0
      if (present(f1_eval)) f1_eval(1:nout) = 0d0
      if (present(f2_eval)) f2_eval(1:nout) = 0d0
      if (present(f3_eval)) f3_eval(1:nout) = 0d0
   endif

end subroutine spline_eval_spl

!------------------------------------------------------------------------------------------------------------

subroutine solve_cubic_segm( s_a, s_b, f_a, c1, c2, c3, f_b, fval, sl, ldebug, ierror)
!--function: solve f(s) = c3 sl^3 + c2 sl^2 + c1 sl + f_a - fval = 0 ,
!            with sl the local coordinate in the segment.
   implicit none
!--subroutine arguments
   integer,      intent(in)  :: ldebug
   real(kind=8), intent(in)  :: s_a, s_b, f_a, f_b, c1, c2, c3, fval
   real(kind=8), intent(out) :: sl
   integer,      intent(out) :: ierror
!--local variables
   logical     :: use_newton

   use_newton = .true.

   if (use_newton) then

      call solve_cubic_newton ( s_a, s_b, f_a, c1, c2, c3, f_b, fval, sl, ldebug, ierror)

   else

      call solve_cubic_cardano( s_a, s_b, f_a, c1, c2, c3,      fval, sl, ldebug)

   endif

end subroutine solve_cubic_segm

!------------------------------------------------------------------------------------------------------------

subroutine solve_cubic_cardano( s_a, s_b, f_a, c1, c2, c3, fval, sl, ldebug)
!--function: solve f(s) = c3 sl^3 + c2 sl^2 + c1 sl + f_a - fval = 0 , using Cardano formulas,
!            with sl the local coordinate in the segment.
! solve 3rd degree equation f(s) = y, returning one root in interval [s_a,s_b]
   implicit none
!--subroutine arguments
   integer,      intent(in)  :: ldebug
   real(kind=8), intent(in)  :: s_a, s_b, f_a, c1, c2, c3, fval
   real(kind=8), intent(out) :: sl
!--local variables
   real(kind=8)    :: c0, s0, s4, delta0, delta1, discr, xr1, xr2, xr3, xm, rmax, cmax,                 &
                      tiny_a, tiny_b, tiny_cfc
   complex(kind=8) :: coefc, ksi, xc1, xc2, xc3, fc1, fc2, fc3

   c0  = f_a - fval
   s0  = 0d0                ! interval [s_a,s_b] in local coordinates
   s4  = s_b - s_a

   if (ldebug.ge.3) then
      write(bufout,'(4(a,g12.4))') ' cardano: c3 =', c3, ', c2 =', c2,', c1 =', c1, ', c0 =', c0
      call write_log(1, bufout)
   endif

   tiny_a   = 1d-10
   tiny_b   = 1d-10
   tiny_cfc = 1d-20

   if (abs(c2).lt.tiny_b .and. abs(c3).lt.tiny_a) then

      ! segment function is linear, within roundoff precision

      if (ldebug.ge.5) call write_log(' cardano: case 1')
      sl    = -c0 / c1

   elseif (abs(c3).lt.tiny_a) then

      ! segment function is quadratic, within roundoff precision

      if (ldebug.ge.5) call write_log(' cardano: case 2')

      discr = c1**2 - 4d0 * c2 * c0
      if (discr.lt.0) then
         sl = -1d9      ! outside [s0, s4]
      else
         xr1 = (-c1 - sqrt(discr)) / (2d0 * c2)
         xr2 = (-c1 + sqrt(discr)) / (2d0 * c2)
         if (xr1.ge.s0 .and. xr1.le.s4) then
            sl = xr1
         else
            sl = xr2
         endif
      endif

   else

     ! Cardano''s method

     delta0 = c2**2 - 3d0 * c3 * c1
     delta1 = 2d0 * c2**3 - 9d0 * c3 * c2 * c1 + 27d0 * c3**2 * c0
     discr  = (4d0 * delta0**3 -  delta1**2) / (27d0 * c3**2)

     if (abs(discr).gt.tiny_cfc) then
        coefc  = ( ( delta1 + sqrt(dcmplx(delta1**2 - 4d0 * delta0**3)) ) / 2d0 )**(1d0/3d0)
        if (abs(coefc).lt.tiny_cfc) then        ! delta1<0, delta0=0
           coefc  = ( ( delta1 - sqrt(dcmplx(delta1**2 - 4d0 * delta0**3)) ) / 2d0 )**(1d0/3d0)
        endif
        ksi  = (-1d0 + sqrt(dcmplx(-3d0))) / 2d0

        if (ldebug.ge.6) then
           write(bufout,'(2(a,f8.4,sp,f8.4,"i"))') ' coefc=',coefc,', ksi  =',ksi
           call write_log(1, bufout)
        endif
     endif

     if (abs(discr).lt.tiny_cfc .and. abs(delta0).lt.tiny_cfc) then

        ! triple real root

        xr1 = -c2 / (3d0 * c3)
        xr2 = xr1
        xr3 = xr1
        sl  = xr1
        if (ldebug.ge.3) call write_log(' cardano: triple root')

     elseif (abs(discr).lt.tiny_cfc) then

        ! double real root + single real root

        xr1 = (4d0 * c3 * c2 * c1 - 9d0 * c3**2 * c0 - c2**3) / (c3 * delta0) ! single
        xr2 = (9d0 * c3 * c0 - c2 * c1) / (2d0 * delta0)                      ! double
        xr3 = xr2
        if (xr1.ge.s0 .and. xr1.le.s4) then
           sl = xr1
        else
           sl = xr2
        endif
        if (ldebug.ge.3) call write_log(' cardano: double root')

     elseif (discr.lt.0d0) then

        ! one real root + two non-real complex conjugate roots
        ! coefc can be complex especially when delta1 < 0.

        xc1 = -( c2 + ksi    * coefc + delta0 / (ksi    * coefc)) / (3d0 * c3)
        xc2 = -( c2 + ksi**2 * coefc + delta0 / (ksi**2 * coefc)) / (3d0 * c3)
        xc3 = -( c2 +          coefc + delta0 /           coefc ) / (3d0 * c3)

        ! select root with smallest imaginary part (theoretically zero)

        if (abs(dimag(xc1)).lt.min(abs(dimag(xc2)),abs(dimag(xc3)))) then
           sl = dble(xc1)
        elseif (abs(dimag(xc2)).lt.abs(dimag(xc3))) then
           sl = dble(xc2)
        else
           sl = dble(xc3)
        endif
        if (ldebug.ge.3) call write_log(' cardano: one real root and two non-real complex conjugate roots')

     else

        ! three distinct real roots

        xc1 = -( c2 + ksi    * coefc + delta0 / (ksi    * coefc)) / (3d0 * c3)
        xc2 = -( c2 + ksi**2 * coefc + delta0 / (ksi**2 * coefc)) / (3d0 * c3)
        xc3 = -( c2 +          coefc + delta0 /           coefc ) / (3d0 * c3)

        rmax = max( abs(dble(xc1)),  abs(dble(xc2)),  abs(dble(xc3))  )
        cmax = max( abs(dimag(xc1)), abs(dimag(xc2)), abs(dimag(xc3)) )

        if (cmax.gt.1d-12*rmax) then
           write(bufout,'(a,g12.4,a)') ' cardano: Impossible: nonzero imag.part', cmax,' ???'
           call write_log(1, bufout)
           write(bufout,'(a,3g12.4,/,a,3g12.4)') '    real:', dble(xc1), dble(xc2), dble(xc3),          &
                ' complex:', dimag(xc1), dimag(xc2), dimag(xc3)
           call write_log(2, bufout)

           write(bufout,'(6(a,f16.10,sp,f16.10,"i",:,/))')                                              &
                '                    xc2', dble(xc2), dimag(xc2),                                       &
                '                     c2', dble(c2), 0d0,                                               &
                '         ksi**2 * coefc', dble(ksi**2 * coefc), dimag(ksi**2 * coefc),                 &
                'delta0 / ksi**2 * coefc', dble(delta0/(ksi**2*coefc)), dimag(delta0/(ksi**2*coefc)),   &
                '             1 / (3 c3)', dble(1d0/(3d0*c3)), 0d0
           call write_log(5, bufout)
        endif

        ! select a root in interval [s0,s4] or closest to interval [s0,s4]

        xm  = (s0 + s4) / 2d0
        xr1 = dble(xc1)
        xr2 = dble(xc2)
        xr3 = dble(xc3)
        if (abs(xr1-xm).lt.min(abs(xr2-xm),abs(xr3-xm))) then
           sl = xr1
        elseif (abs(xr2-xm).lt.abs(xr3-xm)) then
           sl = xr2
        else
           sl = xr3
        endif
        if (ldebug.ge.3) call write_log(' cardano: three distinct real roots')

        if (ldebug.ge.4) then
           fc1 = ((c3 * xc1 + c2) * xc1 + c1) * xc1 + c0
           fc2 = ((c3 * xc2 + c2) * xc2 + c1) * xc2 + c0
           fc3 = ((c3 * xc3 + c2) * xc3 + c1) * xc3 + c0

           write(bufout,'(a,/,6x,3(f12.4,sp,f12.4,"i"),/,6x,3(f12.4,sp,f12.4,"i"))')                    &
                        ' solutions x1, x2, x3:', dble(xc1), dimag(xc1), dble(xc2), dimag(xc2),         &
                        dble(xc3), dimag(xc3), dble(fc1), dimag(fc1), dble(fc2), dimag(fc2),            &
                        dble(fc3), dimag(fc3)
           call write_log(3, bufout)
        endif

     endif

   endif

end subroutine solve_cubic_cardano

!------------------------------------------------------------------------------------------------------------

subroutine solve_cubic_newton( s_a, s_b, f_a, c1, c2, c3, f_b, fval, sl, ldebug, ierror)
!--function: solve f(s) = c3 sl^3 + c2 sl^2 + c1 sl + f_a - fval = 0 using Newton-Raphson,
!            with sl the local coordinate in the segment.
   implicit none
!--subroutine arguments
   integer,      intent(in)  :: ldebug
   real(kind=8), intent(in)  :: s_a, s_b, f_a, f_b, c1, c2, c3, fval
   real(kind=8), intent(out) :: sl
   integer,      intent(out) :: ierror
!--local variables
   integer,      parameter  :: maxit = 50
   real(kind=8), parameter  :: tol_f = 1d-8, tiny_c3 = 1d-10, tiny_c2 = 1d-10, tiny_df = 1d-10
   integer                  :: iter
   logical                  :: s1_in_s04, s3_in_s04
   real(kind=8)             :: discr, s0, s1, s2, s3, s4, f_s, f0, f1, f2, f3, f4, df_ds, ds, ds_max,   &
                               tmp, tol_s

   ierror = 0

   ! solve f(s) = a sl^3 + b sl^2 + c sl + d - f_out = 0, with sl the local coordinate in the segment
   ! Using Newton-Raphson, cannot jump over extremal values, needs appropriate s^0

   if (ldebug.ge.5) then
      write(bufout,'(4(a,g16.8),a)') ' segment: s=[',s_a,',',s_b,'], f=[',f_a,',',f_b,']'
      call write_log(1, bufout)
      write(bufout,'(5(a,g16.8),a)') ' solve f(sl)=',c3,' sl^3 +',c2,' sl^2 +',c1,' sl +',f_a,' -',     &
                               fval,' = 0'
      call write_log(1, bufout)
   endif

   ! check for local extremal values: zeros of quadratic f'(sl) = 3a sl^2 + 2b sl + c = 0
   ! [s0, s4] act as bracket where the solution is sought

   s0 = 0d0
   s2 = 1d10
   s4 = s_b - s_a
   discr = (2d0*c2)**2 - 4d0 * (3d0*c3) * c1

   if (discr.le.0d0) then

      ! no local extremal values: use start-point of iseg as initial estimate, sl = 0

      sl    = s0
      f_s   = f_a - fval
      f0    = f_s
      df_ds = c1

      if (ldebug.ge.4) then
         write(bufout,'(2(a,g14.6))') ' solve_cubic: D=',discr,', sl=',sl
         call write_log(1, bufout)
      endif

   elseif (abs(c3).lt.tiny_c3 .and. abs(c2).lt.tiny_c2) then

      ! segment function is linear, within roundoff precision

      sl    = s0
      f_s   = f_a - fval
      f0    = f_s
      df_ds = c1

      if (ldebug.ge.4) then
         write(bufout,'(3(a,g14.6))') ' solve_cubic: c3=',c3,', c2=',c2,', sl=',sl
         call write_log(1, bufout)
      endif

   elseif (abs(c3).lt.tiny_c3) then

      ! |c2|>tiny: quadratic f, linear f' = c1 + 2*c2*s

      s1 = -0.5d0 * c1 / c2                               ! extremal value
      s1_in_s04 = ((s4-s1)*(s1-s0) .gt. 0d0)

      f0 = f_a - fval
      f1 = ((c3 * s1 + c2) * s1 + c1) * s1 + f_a - fval
      f4 = f_b - fval

      ! there can be a zero <s1 and one >s1

      if (.not.s1_in_s04) then
         ! extremum outside [s0,s4] can be ignored
         sl    = s0
         f_s   = f_a - fval
         df_ds = c1
      elseif (f0*f1.gt.0d0) then
         ! extremum on same side as s0, start from s4
         sl  = s4
         f_s = f4
         s0  = s1       ! low point of bracket
         f0  = f1
      else
         ! extremum on same side as s4, start from s0
         sl  = s0
         f_s = f0
         s4  = s1       ! high point of bracket
         f4  = f1
      endif

   else

      ! true cubic equation with D>0

      s1 = ( -(2d0*c2) - sqrt(discr) ) / (2d0 * (3d0*c3)) ! extremal value
      s2 = ( -(2d0*c2)               ) / (2d0 * (3d0*c3)) ! bending point
      s3 = ( -(2d0*c2) + sqrt(discr) ) / (2d0 * (3d0*c3)) ! extremal value
      if (s1.gt.s3) then
         tmp = s1
         s1 = s3
         s3 = tmp
      endif

      f0 = f_a - fval
      f1 = ((c3 * s1 + c2) * s1 + c1) * s1 + f_a - fval
      f2 = ((c3 * s2 + c2) * s2 + c1) * s2 + f_a - fval
      f3 = ((c3 * s3 + c2) * s3 + c1) * s3 + f_a - fval
      f4 = f_b - fval
   
      ! there can be a zero <s1, one between s1 and s3, and one >s3
      ! select initial estimate: start-point, mid-point or end-point of segment iseg

      s1_in_s04 = ((s4-s1)*(s1-s0) .gt. 0d0)
      s3_in_s04 = ((s4-s3)*(s3-s0) .gt. 0d0)

      if (f1*f3.gt.0d0 .and. c3*f1.gt.0d0) then
         ! one zero, to the left of s1
         sl  = s0
         f_s = f0
         s4  = s1       ! high point of bracket
         f4  = f1
      elseif (f1*f3.gt.0d0) then
         ! one zero, to the right of s3
         sl  = s4
         f_s = f4
         s0  = s3       ! low point of bracket
         f0  = f3
      elseif (s1_in_s04 .and. s3_in_s04) then
         ! three zeros, middle one between s0 and s4
         sl  = (s1 + s3) / 2d0
         f_s = ((c3 * sl + c2) * sl + c1) * sl + f_a - fval
         s0  = s1       ! low point of bracket
         f0  = f1
         s4  = s3       ! high point of bracket
         f4  = f3
      elseif ((s1_in_s04 .and. f0*f1.gt.0d0) .or. (s3_in_s04 .and. f0*f3.gt.0d0)) then
         ! three zeros, right one between s0 and s4
         sl  = s4
         f_s = f4
         s0  = s1       ! low point of bracket
         f0  = f1
      else
         ! three zeros, left one between s0 and s4
         sl  = s0
         f_s = f0
      endif
      df_ds = (3d0 * c3 * sl + 2d0 * c2) * sl + c1

   endif

   ! Newton-Raphson: f(s+ds) = f(s) + ds * f'(s) = 0   -->  ds = - f(s) / f'(s)

   iter  = 0
   if (ldebug.ge.5) then
      write(bufout,12) iter, sl, s0, s4, f_s, df_ds 
      call write_log(1, bufout)
 12   format(' It ',i2,': sl = ',f8.4,' in [',f8.4,',',f8.4,'], f(sl) = ',g16.8,', f''(sl) = ',g16.8)
   endif

   do while(abs(f_s).ge.tol_f .and. iter.lt.maxit)
      iter  = iter + 1

      ! use Newton step, revert to bi-section in case of small f'

      if (abs(df_ds).gt.tiny_df) then
         ds    = - f_s / df_ds
      else
         ds    = (s0 + s4) / 2d0  - sl    ! small derivative: use bisection
      endif

      ! revert to bi-section if Newton jumps out of bracket
      ! problem with rounding error w.r.t. updated s0/s4?

      ! if (sl+ds.lt.s0 .or. sl+ds.gt.s4) ds = (s0 + s4) / 2d0 - sl

      ! restrict step esp. in first iterations, cf. doubling of step until bracket is found

      ds_max = (s4 - s0) * 2d0**iter;
      if (ds.gt.ds_max) then
         ds = ds_max
      elseif (ds.lt.-ds_max) then
         ds = -ds_max
      endif

      ! first step: avoid jumping over bending point s2

      if (discr.gt.0d0 .and. iter.le.1) then
         if (sl.lt.s2) then
            sl = min(s2, sl+ds)
         else
            sl = max(s2, sl+ds)
         endif
      else
         sl = sl + ds
      endif

      ! compute new function and derivative value 

      f_s   = ((      c3 * sl +       c2) * sl  + c1) * sl + f_a - fval
      df_ds =  (3d0 * c3 * sl + 2d0 * c2) * sl  + c1

      ! update bracket [s0, s4]

      if (f_s*f0.gt.0d0) then
         s0 = sl
         f0 = f_s
      else
         s4 = sl
         f4 = f_s
      endif
   
      if (ldebug.ge.5 .or. (iter.ge.5 .and. ldebug.ge.2)) then
         write(bufout,12) iter, sl, s0, s4, f_s, df_ds 
         call write_log(1, bufout)
      endif
   enddo

   tol_s = tol_f / max(1d-3, abs(df_ds))

   if (abs(f_s).ge.tol_f) then
      ierror = maxit
   elseif (sl.lt.s0-tol_s .or. sl.gt.s4+tol_s) then
      ierror = -1
   else
      ierror = 0
   endif

end subroutine solve_cubic_newton

!------------------------------------------------------------------------------------------------------------

subroutine spline_get_s_at_f_spl( spl, ikarg, nout, f_out, s_out, my_ierror )
!--function: interpolate parametric spline 'spl' in direction ikarg, determine s_out at the requested
!            f_out-positions
!      Note: the s-values are not uniquely defined in case the inputs are non-monotonic,
!            e.g. back-of-flange, y_w oscillating around 0.
!            A unique result is obtained using the 'top view administration', if available.
!            Else, s := s(1) for all f <= f(1) and s := s(np) for f >= f(np), this is a conservative
!            approach for the range of s-values needed to cover [f0,f1].
!--subroutine arguments
   type(t_spline)                :: spl
   integer,          intent(in)  :: ikarg
   integer,          intent(in)  :: nout
   real(kind=8),     intent(in)  :: f_out(nout)
   real(kind=8),     intent(out) :: s_out(nout)
   integer,          intent(out) :: my_ierror
!--local variables
   integer                  :: sub_ierror, iseg_debug, iout, itop, iuni, ip0, ip1, np01, iseg, ip
   real(kind=8)             :: sl
   character(len=1)         :: nam_f
   real(kind=8), dimension(:),   pointer :: a, b, c, d
   real(kind=8), dimension(:,:), pointer :: coef

   my_ierror = 0

   associate(np => spl%npnt, s => spl%s)

   ! select simple-spline according to argument ikarg

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
      call write_log('INTERNAL ERROR (spline_get_s_at_f): ikarg invalid.')
      call abort_run()
   endif

   if (.not.associated(coef)) then
      if (.true.) call write_log('spline_get_s_at_f_spl: coef not associated')
      s_out(1:nout) = 0d0
      return
   endif

   ! set pointers to spline coefficients a = a3, b = a2, c = a1, d = a0

   a => coef(:,4)
   b => coef(:,3)
   c => coef(:,2)
   d => coef(:,1)

   if (ldebug.ge.2) then
      write(bufout,'(a,i6,2a,i6,1x,2a)') ' get_s_at_f: input spline has',np-1,' segments, ',            &
                'interpolating s to',nout, nam_f,'-positions'
      call write_log(1, bufout)
   endif

   if (ldebug.ge.1) then
      if (.not.is_strictly_monotonic(np, d)) then
         call write_log(' WARNING: spline_get_s_at_f: spline ' // nam_f //                              &
                '-values are NOT strictly monotonic')
         if (ikarg.eq.ikYDIR .and. spl%nsec_top.ge.1) then
            call write_log('          using top-view administration.')
         elseif (ldebug.ge.5) then
            do ip = 1, np
               write(bufout,'(a,i4,3a,2f12.6)') ' i=', ip, ', s, ',nam_f,'=', s(ip), d(ip)
               call write_log(1, bufout)
            enddo
         endif
      endif

      call check_monotone(nout, f_out, 0, iseg)
      if (iseg.gt.0) then
         call write_log(' WARNING: spline_get_s_at_f: requested ' // nam_f //                           &
                '-positions are not monotonic')
         if (ldebug.ge.5) then
            do iout = max(1,iseg-2), min(iseg+2,nout)
               write(bufout,'(a,i4,3a,f12.6)') ' i=', iout, ', ',nam_f,'_out=', f_out(iout)
               call write_log(1, bufout)
            enddo
         endif
      endif
   endif

   ! assuming that f_out is sorted, either in ascending (rail) or descending order (wheel)

   iseg = 0

   do iout = 1, nout

      ! find the segment iseg [ d(iseg), d(iseg+1) ] containing f_out

      if (ikarg.ne.ikYDIR .or. spl%nsec_top.le.0) then

         ! original case: searching whole array d for segment iseg containing f_out

         ip0  = 1
         ip1  = np
         call locate_segment( np, d, f_out(iout), iseg )

      else

         ! extended case: using subdivision into monotonous/uni-valued sections

         ! find top-section [f_j, f_{j+1}] containing f_out, j = itop

         call locate_segment( spl%nsec_top+1, spl%ysec_top, f_out(iout), itop )
         itop = max(1, min(spl%nsec_top, itop))

         ! determine the corresponding uni-valued section number iuni

         iuni = spl%iuni_top(itop)

         ! search f_out in uni-valued section iuni

         ip0  = spl%ipnt_uniy(iuni)
         ip1  = spl%ipnt_uniy(iuni+1) - 1
         np01 = ip1 - ip0 + 1

         call locate_segment( np01, d(ip0:ip1), f_out(iout), iseg )

         iseg = ip0-1 + iseg

         if (ldebug.ge.4) then
            write(bufout,'(3a,i4,a,f12.6,a,i3,3a,2(f12.6,a),3(i4,a))') ' ',nam_f,'out(',iout,')=',      &
                f_out(iout), ' lies in itop=',itop,', ',nam_f,'=[', spl%ysec_top(itop),',',             &
                spl%ysec_top(itop+1),'], iuni=',iuni,', ip=[',ip0,',',ip1,']'
            call write_log(1, bufout)
            write(bufout,'(3a,f12.6,3(a,i4),3a,2(f12.6,a))') ' found ',nam_f,'=',f_out(iout),' in seg', &
                iseg-ip0+1, ' of iuni=',iuni,', overall iseg=',iseg,', ',nam_f,'=[',d(iseg),',',d(iseg+1),']'
            call write_log(1, bufout)
         endif

      endif

      ! if f(i) is before start (ascending: f(i) <= d(1) < d(n); descending: f(i) => d(1) > d(n))

      if (iseg.lt.ip0) then

         ! set s slightly before s(1), for further processing in get_xz_at_y

         s_out(iout) = s(ip0) - 1d-9

         if (ldebug.ge.2 .or. iout.eq.ii_debug) then
            write(bufout,'(2(3a,i6),a,f12.6)') ' ',nam_f,'_out(',iout,') <= ',nam_f,'_spl(ip0), s(',    &
                iout,')=',s_out(iout)
            call write_log(1, bufout)
         endif

      ! if f(i) is after end (ascending: f(i) >= d(n) > d(n); descending: f(i) <= d(n) < d(1))

      elseif ( iseg.ge.ip1 ) then

         ! set s slightly beyond s(end), for further processing in get_xz_at_y

         s_out(iout) = s(ip1) + 1d-9

         if (ldebug.ge.2 .or. iout.eq.ii_debug) then
            write(bufout,'(2(3a,i6),a,f12.6)') ' ',nam_f,'_out(',iout,') >= ',nam_f,'_spl(ip1), s(',    &
                iout,')=',s_out(iout)
            call write_log(1, bufout)
         endif

      ! else: f(i) lies within range [d(1), d(n)], solve segment equation

      else

         if (ldebug.ge.4 .or. iout.eq.ii_debug) then
            write(bufout,'(3a,i5,a,g14.6,a,i5,3a,2(g14.6,a))') ' ',nam_f,'_out(',iout,')=',f_out(iout), &
                ' lies in segment',iseg,' with ',nam_f,' = [',d(iseg),',', d(iseg+1),']'
            call write_log(1, bufout)
         endif

         ! solve f(s) = a sl^3 + b sl^2 + c sl + d - f_out = 0, with sl the local coordinate in the segment

         iseg_debug = ldebug
         if (iout.eq.ii_debug) iseg_debug = max(5, ldebug)

         call solve_cubic_segm( s(iseg), s(iseg+1), d(iseg), c(iseg), b(iseg), a(iseg), d(iseg+1),      &
                f_out(iout), sl, iseg_debug, sub_ierror)

         if (sub_ierror.gt.0 .and. (ldebug.ge.1 .or. iout.eq.ii_debug)) then
            write(bufout,*) 'get_s_at_f: No convergence for iout=',iout
            call write_log(1, bufout)
            if (my_ierror.eq.0) my_ierror = sub_ierror
         elseif (sub_ierror.eq.-1 .and. (ldebug.ge.1 .or. iout.eq.ii_debug)) then
            write(bufout,'(a,i6,a,f8.4,a,i4,2(a,f8.4),a)') 'get_s_at_f: for iout=',iout,', s=',         &
                s(iseg)+sl, ' lies outside segment',iseg,' with s= [',s(iseg),',',s(iseg+1),']'
            call write_log(1, bufout)
            ! error not propagated upwards
         endif

         s_out(iout) = s(iseg) + sl

         if (ldebug.ge.4 .or. iout.eq.ii_debug) then
            write(bufout,'(a,i5,a,g12.4,a,i5,2(a,g12.4),a)') ' s_out(',iout,')=',s_out(iout),           &
                ' lies in segment',iseg,' with s=[',s(iseg),',', s(iseg+1),']'
            call write_log(1, bufout)
         endif
         if (ldebug.ge.2 .or. iout.eq.ii_debug) then
            write(bufout,'(3a,i5,a,f10.3, a,i5,2a,5(a,f10.3))') ' ',nam_f,'_out(',iout,')=',f_out(iout), &
                ' lies in segment',iseg,': ',nam_f,'=[',d(iseg),',', d(iseg+1),'], s=[',s(iseg),',',    &
                s(iseg+1),'], s_out=',s_out(iout)
         endif

      endif ! y<y(1) | y>y(end)

   enddo ! iout

   end associate

end subroutine spline_get_s_at_f_spl

!------------------------------------------------------------------------------------------------------------

subroutine locate_one_extremum(npnt, s_spl, a3, a2, a1, a0, iseg, typ_xtrm, s_xtrm, ierror)
!--function: for a simple spline, determine the s-position of a locally extremal function value
!            --> for segment iseg with f'(s0)*f'(s1)<0, solve f'(s) = 0 with quadratic f'
   implicit none
!--subroutine arguments
   integer,          intent(in)   :: npnt, iseg
   integer,          intent(in)   :: typ_xtrm   ! <0: minimum, 0: dont care, >0: maximum
   real(kind=8),     intent(in)   :: s_spl(npnt), a3(npnt), a2(npnt), a1(npnt), a0(npnt)
   real(kind=8),     intent(out)  :: s_xtrm
   integer,          intent(out)  :: ierror
!--local variables
   real(kind=8), parameter :: tiny_a3 = 1d-10, tiny_a2 = 1d-10
   real(kind=8) :: discr, ds, sloc, sloc1, sloc2, a1loc, floc1, floc2

   ierror = 0
   s_xtrm = 1d10

   if (ldebug.ge.3) then
      write(bufout,'(a,i5,3(a,f9.4),3g12.4,a,i3)') ' locate_xtrem: seg',iseg,': s=[', s_spl(iseg),      &
              ',', s_spl(iseg+1), '], a0-a3=',a0(iseg), a1(iseg), a2(iseg), a3(iseg),', type',typ_xtrm
      call write_log(1, bufout)
   endif

   ! segment cubic:   v(s)  = a0 + a1 * s +     a2 * s^2 +     a3 * s^3
   ! derivative:      v'(s) =      a1     + 2 * a2 * s   + 3 * a3 * s^2
   ! D = b^2 - 4ac = (2*a2)^2 - 4 * (3*a3) * a1

   ds    = s_spl(iseg+1) - s_spl(iseg)
   discr = 4d0 * a2(iseg)**2 - 12d0 * a3(iseg) * a1(iseg)

   if (abs(a3(iseg)).lt.tiny_a3 .and. abs(a2(iseg)).lt.tiny_a2) then

      ! segment function f(s) is linear within roundoff precision
      ! f'(s) is constant - changes sign only in case there's a discontinuity, kink at end of segment

      s_xtrm = s_spl(iseg+1)

      if (ldebug.ge.3) then
         write(bufout,'(a,f11.6)') '   ...linear segment, s_xtrm=', s_xtrm
         call write_log(1, bufout)
      endif

   elseif (discr.le.0d0 .and. a1(iseg)*a1(iseg+1).gt.0d0) then

      ! D<0, same sign for f'(s0) and f'(s1) -- no zeros in f'

      ierror = -1
      call write_log(' INTERNAL ERROR(spline, locate_extremum): a1(sta) * a1(end) > 0')

   elseif (discr.le.0d0) then

      ! D<0, no zeros in f' : f'(s0) * f'(s1) <= 0  indicates a discontinuity in f', kink at end of segment

      s_xtrm = s_spl(iseg+1)

      if (ldebug.ge.3) then
         write(bufout,'(a,g12.4,a,i4,a,f11.6)') '   ...D=',discr,'<0, kink at end of seg',iseg,      &
                  ', s_xtrm=', s_xtrm
         call write_log(1, bufout)
      endif

   else

      ! compute zeros sloc1 and sloc2

      if (abs(a3(iseg)).lt.tiny_a3) then          ! |a2|>tiny: quadratic f, linear f' = a1 + 2*a2*s
         sloc1 = -0.5d0 * a1(iseg) / a2(iseg)
         sloc2 =  1d10
      else                                        ! |a3|>tiny: cubic f, quadratic f'
         sloc1 = (-2d0 * a2(iseg) - sqrt(discr)) / (6d0 * a3(iseg))
         sloc2 = (-2d0 * a2(iseg) + sqrt(discr)) / (6d0 * a3(iseg))
      endif

      if (ldebug.ge.3) then
         write(bufout,'(15x,3(a,f11.4))') 'zeros at sloc1=',sloc1,', sloc2=',sloc2,', ds=',ds
         call write_log(1, bufout)
      endif

      if (min(sloc1,sloc2).ge.0d0 .and. max(sloc1,sloc2).le.ds) then

         ! two zeros within [0, ds]

         floc1 = a0(iseg) + a1(iseg) * sloc1 + a2(iseg) * sloc1**2 + a3(iseg) * sloc1**3
         floc2 = a0(iseg) + a1(iseg) * sloc2 + a2(iseg) * sloc2**2 + a3(iseg) * sloc2**3

         if (typ_xtrm.eq.0) then

            ! searching locally minimum AND locally maximum values

            s_xtrm = s_spl(iseg) + 0.5d0 * (sloc1 + sloc2)

            ierror = -2
            call write_log(' ERROR(spline, locate extremum): two zeros in segment')

         elseif ( (typ_xtrm.lt.0 .and. floc1.le.floc2) .or.                                             &
                  (typ_xtrm.gt.0 .and. floc1.ge.floc2) ) then

            ! searching (minimum & f1<=f2) or (maximum & f1>=f2)

            s_xtrm = s_spl(iseg) + sloc1

         else

            ! searching (minimum & f1> f2) or (maximum & f1< f2)

            s_xtrm = s_spl(iseg) + sloc2

         endif

         if (ldebug.ge.1 .and. typ_xtrm.ne.0)                                                           &
            call write_log(' WARNING(spline, locate extremum): solution need not be overall min/max')

      elseif (max(sloc1,sloc2).lt.0d0 .or. min(sloc1,sloc2).gt.ds .or.                                  &
              (min(sloc1,sloc2).lt.0d0 .and. max(sloc1,sloc2).gt.ds)) then

         ! no zeros within [0, ds]

         ! check for f' jumping at end of segment

         sloc  = ds
         a1loc = a1(iseg) + 2d0 * a2(iseg) * ds + 3d0 * a3(iseg) * ds**2

         if (a1loc * a1(iseg+1).lt.0d0) then      ! kink at end of segment

            s_xtrm = s_spl(iseg+1)

            if (ldebug.ge.3) then
               write(bufout,'(a,f11.6)') '    ...discontin f'', s_xtrm=', s_xtrm
               call write_log(1, bufout)
            endif

         else                                     ! continuous: should have zero within segment

            call write_log(' ERROR: no zero in segment')
            write(bufout,'(a,i4,4(a,g12.4))') ' iseg=',iseg,': sloc1,2=',sloc1,',',sloc2,         &
                          ', ds=',ds, ', D=',discr
            call write_log(1, bufout)
         endif

      elseif (sloc1.ge.0d0 .and. sloc1.le.ds) then

         ! first zero within [0, ds]

         s_xtrm = s_spl(iseg) + sloc1

         if (ldebug.ge.3) then
            write(bufout,'(a,f11.6)') '    ...using sloc1, s_xtrm=', s_xtrm
            call write_log(1, bufout)
         endif

      else

         ! second zero within [0, ds]

         s_xtrm = s_spl(iseg) + sloc2

         if (ldebug.ge.3) then
            write(bufout,'(a,f11.6)') '    ...using sloc2, s_xtrm=', s_xtrm
            call write_log(1, bufout)
         endif

      endif

   endif ! discr>0

end subroutine locate_one_extremum

!------------------------------------------------------------------------------------------------------------

subroutine locate_extremal_values(npnt, s_spl, a3, a2, a1, a0, nxtrm, s_xtrm, my_ierror)
!--function: for a simple spline, determine the s-position with locally extremal function values
!            --> for all segments segments with f'(s0)*f'(s1)<0, solve f'(s) = 0 with quadratic f'
!            NOTE: this subroutine overlaps with spline_get_s_at_minval_arr in m_spline_get.f90.
   implicit none
!--subroutine arguments
   integer,          intent(in)   :: npnt
   real(kind=8),     intent(in)   :: s_spl(npnt), a3(npnt), a2(npnt), a1(npnt), a0(npnt)
   real(kind=8),     intent(out)  :: s_xtrm(npnt)
   integer,          intent(out)  :: nxtrm, my_ierror
!--local variables
   real(kind=8), parameter :: tiny_a3 = 1d-10, tiny_a2 = 1d-10
   integer      :: iseg, sub_ierror

   my_ierror = 0
   nxtrm  = 0

   ! loop over all segments with opposite f'(s) at both ends

   do iseg = 1, npnt-1
      if (a1(iseg)*a1(iseg+1).le.0d0) then

         nxtrm = nxtrm + 1
         call locate_one_extremum(npnt, s_spl, a3, a2, a1, a0, iseg, 0, s_xtrm(nxtrm), sub_ierror)

         ! ignore segments with two zeros in f'

         if (sub_ierror.lt.0) nxtrm = nxtrm - 1

      endif ! f'*f'<0
   enddo ! iseg

end subroutine locate_extremal_values

!------------------------------------------------------------------------------------------------------------

subroutine spline_shift(spl, dx, dy, dz)
!--function: shift a spline in cartesian coordinates: spl = spl + [dx; dy; dz]
   implicit none
!--subroutine arguments
   type(t_spline)  :: spl
   real(kind=8)    :: dx, dy, dz
!--local variables
   integer         :: npnt

   ! if spline does not have x-data: initialize x-spline == 0

   if (.not.spl%has_xdata .and. abs(dx).gt.1d-9) call spline_initx(spl)

   ! apply shift to spline coefficients d

   npnt = spl%npnt
   if (spl%has_xdata) spl%ax0(1:npnt) = spl%ax0(1:npnt) + dx
   spl%ay0(1:npnt) = spl%ay0(1:npnt) + dy
   spl%az0(1:npnt) = spl%az0(1:npnt) + dz

end subroutine spline_shift

!------------------------------------------------------------------------------------------------------------

subroutine spline_mirror_y(spl)
!--function: mirror a spline in cartesian coordinates wrt plane Oxz: spl.y = -spl.y
   implicit none
!--subroutine arguments
   type(t_spline)  :: spl
!--local variables
   integer         :: npnt

   ! apply mirroring to spline coefficients ay0, ay1, ay3

   npnt   = spl%npnt
   spl%ay0(1:npnt) = -spl%ay0(1:npnt)
   spl%ay1(1:npnt) = -spl%ay1(1:npnt)
   spl%ay3(1:npnt) = -spl%ay3(1:npnt)

end subroutine spline_mirror_y

!------------------------------------------------------------------------------------------------------------

subroutine spline_rotate(spl, rot, xc, yc, zc)
!--function: rotate a spline in cartesian coordinates by rotation matrix about center-point [xc;yc;zc]
   implicit none
!--subroutine arguments
   type(t_spline)           :: spl          ! input/output spline
   type(t_rotmat)           :: rot          ! matrix with new orientation of original unit vectors
   real(kind=8), intent(in) :: xc, yc, zc   ! rotation origin
!--local variables
   integer      :: ii
   real(kind=8) :: xrel, yrel, zrel, tmpx, tmpy

   ! if curve in Oyz moves outside Oyz: initialize x-spline == 0

   if (.not.spl%has_xdata .and. (abs(rot%r(4)).gt.1d-9 .or. abs(rot%r(7)).gt.1d-9)) call spline_initx(spl)

   ! apply rotation to spline coefficients {a3, a2, a1}, without origin, and to {a0}, relative to (yc,zc)

   if (spl%has_xdata) then

      do ii = 1, spl%npnt
         tmpx        =             spl%ax3(ii)
         tmpy        =                               spl%ay3(ii)
         spl%ax3(ii) =      rot%r(1) * tmpx + rot%r(4) * tmpy + rot%r(7) * spl%az3(ii)
         spl%ay3(ii) =      rot%r(2) * tmpx + rot%r(5) * tmpy + rot%r(8) * spl%az3(ii) 
         spl%az3(ii) =      rot%r(3) * tmpx + rot%r(6) * tmpy + rot%r(9) * spl%az3(ii) 

         tmpx        =             spl%ax2(ii)
         tmpy        =                               spl%ay2(ii)
         spl%ax2(ii) =      rot%r(1) * tmpx + rot%r(4) * tmpy + rot%r(7) * spl%az2(ii)
         spl%ay2(ii) =      rot%r(2) * tmpx + rot%r(5) * tmpy + rot%r(8) * spl%az2(ii) 
         spl%az2(ii) =      rot%r(3) * tmpx + rot%r(6) * tmpy + rot%r(9) * spl%az2(ii) 

         tmpx        =             spl%ax1(ii)
         tmpy        =                               spl%ay1(ii)
         spl%ax1(ii) =      rot%r(1) * tmpx + rot%r(4) * tmpy + rot%r(7) * spl%az1(ii)
         spl%ay1(ii) =      rot%r(2) * tmpx + rot%r(5) * tmpy + rot%r(8) * spl%az1(ii) 
         spl%az1(ii) =      rot%r(3) * tmpx + rot%r(6) * tmpy + rot%r(9) * spl%az1(ii) 

         xrel        =     spl%ax0(ii) - xc
         yrel        =                       spl%ay0(ii) - yc
         zrel        =                                         spl%az0(ii) - zc
         spl%ax0(ii) = xc + rot%r(1) * xrel + rot%r(4) * yrel + rot%r(7) * zrel
         spl%ay0(ii) = yc + rot%r(2) * xrel + rot%r(5) * yrel + rot%r(8) * zrel
         spl%az0(ii) = zc + rot%r(3) * xrel + rot%r(6) * yrel + rot%r(9) * zrel
      end do

   else

      do ii = 1, spl%npnt
         tmpy        =                               spl%ay3(ii)
         spl%ay3(ii) =                        rot%r(5) * tmpy + rot%r(8) * spl%az3(ii) 
         spl%az3(ii) =                        rot%r(6) * tmpy + rot%r(9) * spl%az3(ii) 

         tmpy        =                               spl%ay2(ii)
         spl%ay2(ii) =                        rot%r(5) * tmpy + rot%r(8) * spl%az2(ii)
         spl%az2(ii) =                        rot%r(6) * tmpy + rot%r(9) * spl%az2(ii)

         tmpy        =                               spl%ay1(ii)
         spl%ay1(ii) =                        rot%r(5) * tmpy + rot%r(8) * spl%az1(ii)
         spl%az1(ii) =                        rot%r(6) * tmpy + rot%r(9) * spl%az1(ii)

         yrel        =                       spl%ay0(ii) - yc
         zrel        =                                         spl%az0(ii) - zc
         spl%ay0(ii) = yc +                   rot%r(5) * yrel + rot%r(8) * zrel
         spl%az0(ii) = zc +                   rot%r(6) * yrel + rot%r(9) * zrel
      end do

   endif

   ! clear information on uni-valued sections and top view

   spl%nsec_uniy = 0
   spl%nsec_top = 0

   if (associated(spl%ipnt_uniy)) deallocate(spl%ipnt_uniy)
   if (associated(spl%ysec_top))  deallocate(spl%ysec_top)
   if (associated(spl%iuni_top))  deallocate(spl%iuni_top)
   spl%ipnt_uniy => NULL()
   spl%ysec_top  => NULL()
   spl%iuni_top  => NULL()

end subroutine spline_rotate

!------------------------------------------------------------------------------------------------------------

subroutine spline_roll(spl, roll, yc, zc)
!--function: rotate a spline in cartesian coordinates by roll angle roll [rad] (about x-axis/point [yc;zc])
   implicit none
!--subroutine arguments
   type(t_spline)           :: spl
   real(kind=8), intent(in) :: roll     ! rotation angle [rad]
   real(kind=8), intent(in) :: yc, zc   ! rotation origin
!--local variables
   integer      :: ii
   real(kind=8) :: cs, sn, yrel, zrel, tmp

   cs = cos(roll)
   sn = sin(roll)

   ! apply rotation to spline coefficients {a3, a2, a1}, without origin, and to {a0}, relative to (yc,zc)

   do ii = 1, spl%npnt
      tmp         =           spl%ay3(ii)
      spl%ay3(ii) =      cs * spl%ay3(ii) - sn * spl%az3(ii)
      spl%az3(ii) =      sn * tmp         + cs * spl%az3(ii)

      tmp         =           spl%ay2(ii)
      spl%ay2(ii) =      cs * spl%ay2(ii) - sn * spl%az2(ii)
      spl%az2(ii) =      sn * tmp         + cs * spl%az2(ii)

      tmp         =           spl%ay1(ii)
      spl%ay1(ii) =      cs * spl%ay1(ii) - sn * spl%az1(ii)
      spl%az1(ii) =      sn * tmp         + cs * spl%az1(ii)

      yrel        =    spl%ay0(ii) - yc
      zrel        =                     spl%az0(ii) - zc
      spl%ay0(ii) = yc + cs * yrel        - sn * zrel
      spl%az0(ii) = zc + sn * yrel        + cs * zrel
   end do

   ! clear information on uni-valued sections

   spl%nsec_uniy = 0
   spl%nsec_top = 0

   if (associated(spl%ipnt_uniy)) deallocate(spl%ipnt_uniy)
   if (associated(spl%ysec_top))  deallocate(spl%ysec_top)
   if (associated(spl%iuni_top))  deallocate(spl%iuni_top)
   spl%ipnt_uniy => NULL()
   spl%ysec_top  => NULL()
   spl%iuni_top  => NULL()

end subroutine spline_roll

!------------------------------------------------------------------------------------------------------------

subroutine spline_2glob_or(spl, o, R)
!--function: compute local-to-global conversion for a spline defined with respect to local system (o, R)
!            o = origin of local system w.r.t. global system
!            R = orientation of local system w.r.t. global system
!            after the transformation, spl is defined with respect to the global system
   implicit none
!--subroutine arguments
   type(t_spline), intent(inout) :: spl
   type(t_vec),    intent(in)    :: o
   type(t_rotmat), intent(in)    :: R

   ! rotate points with respect to initial origin o

   call spline_rotate(spl, R, 0d0, 0d0, 0d0)

   ! change origin to the global system o

   call spline_shift(spl, o%v(1), o%v(2), o%v(3))

end subroutine spline_2glob_or

!------------------------------------------------------------------------------------------------------------

subroutine spline_2glob_m(spl, mref)
!--function: compute local-to-global conversion for a spline defined with respect to mref
!            o = mref%o   == origin of reference w.r.t. global system
!            R = mref%rot == orientation of reference w.r.t. global system
!            after the transformation, spl is defined with respect to the global system
   implicit none
!--subroutine arguments
   type(t_marker), intent(in)    :: mref
   type(t_spline), intent(inout) :: spl

   ! rotate points with respect to initial origin mref%o

   call spline_rotate(spl, mref%rot, 0d0, 0d0, 0d0)

   ! change origin to the global system o

   call spline_shift(spl, mref%o%v(1), mref%o%v(2), mref%o%v(3))

end subroutine spline_2glob_m

!------------------------------------------------------------------------------------------------------------

subroutine spline_2loc_m(spl, mref)
!--function: convert spline spl defined in global coordinates to local coordinates according to mref
!            o = mref%o   == origin of new reference w.r.t. global system
!            R = mref%rot == orientation of new reference w.r.t. global system
!            after the transformation, spl is defined with respect to the local system
   implicit none
!--subroutine arguments
   type(t_marker), intent(in)    :: mref
   type(t_spline), intent(inout) :: spl
!--local variables
   type(t_marker)  :: mglb_ref

   ! compute transpose of mref, i.e. the marker for the global system in terms of the local reference

   mglb_ref = marker_transpose(mref)

   ! transform spl from the new local system 'glb' to the new global system 'ref'

   call spline_2glob_m(spl, mglb_ref)

end subroutine spline_2loc_m

!------------------------------------------------------------------------------------------------------------

end module m_spline_def
