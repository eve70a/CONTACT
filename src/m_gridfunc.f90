!------------------------------------------------------------------------------------------------------------
! m_gridfunc - data-structures for grid-based data
!
! Copyright 1993-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_gridfunc
   use m_globals
   use m_markers
   use m_ptrarray
   use m_grids
   implicit none
   private

   !---------------------------------------------------------------------------------------------------------
   ! Codes for designating different types of elements

   public allelm, allint, allext, exter, adhes, slip, plast

   ! allelm   "all elements of the potential contact area"
   ! allint   "all elements in the contact area"
   ! allext   "all elements outside the contact area"
   ! (exter1   code for exterior elements that cannot enter the contact area)
   ! exter    code for exterior elements that might enter the contact area
   ! adhes    code for elements in the adhesion area
   ! slip     code for elements in the slip area
   ! plast    code for elements in the plastic area

   integer,      parameter :: allelm =     -9
   integer,      parameter :: allint =     -8
   integer,      parameter :: allext =     -7
   ! integer,      parameter :: exter1 =     -1
   integer,      parameter :: exter  =      0
   integer,      parameter :: adhes  =      1
   integer,      parameter :: slip   =      2
   integer,      parameter :: plast  =      3

   private  strlen
   integer,      parameter :: strlen =      16

   ! Data types for different kinds of grid functions:

   public t_eldiv
   public t_igrdfnc1
   public t_gridfnc3
   public p_gridfnc3

   ! Functions defined on the data types:

   public eldiv_is_defined
   public eldiv_nullify
   public eldiv_new
   public eldiv_exter
   public eldiv_copy       ! copy contents of element division
 ! public eldiv_copy_full  ! copy structure and values from one eldiv to another
   public eldiv_resize     ! resize to new matching grid, keeping data
   public eldiv_count
   public eldiv_count_atbnd
   public eldiv_cpatches   ! create mask for gf3_msk_copy
   public eldiv_print
   public areas 
   public wrigs 
   public eldiv_destroy

   public if1_nullify
   public if1_new
   public if1_destroy

   public gf3_is_defined
   public gf3_ikrange
   public gf3_nullify
   public gf3_new          ! (re-)initialize structure of grid-func
   public gf3_copy_struc   ! duplicate structure, optionally initialize at 0
   public gf3_copy_full    ! copy structure and values from one gf to another
   public gf3_eldiv
   public gf3_set
   public gf3_copy_xdir    ! copy plain val(ny) to grid.func
   public gf3_copy         ! copy values from one gf to another
   public gf3_msk_copy
   public gf3_resize       ! resize to new matching grid, keeping data
   public gf3_scal
   public gf3_msk_scal
   public gf3_axpy
   public gf3_sum
   public gf3_proj_avg     ! remove the average of a grid function
   public gf3_dot          ! inner product of two gf3-s
   public gf3_nrm2         ! vector 2-norm
   public gf3_rms          ! vector root-mean-square value 
   public gf3_min          ! smallest value: min(val[])
   public gf3_max          ! largest value:  max(val[])
   public gf3_absmax       ! returns  val[ii],  with ii = argmax( |val[]| ), in abs.sense largest value
   public gf3_maxabs       ! returns |val[ii]|, with ii = argmax( |val[]| ), largest abs.value
   public gf3_rotate       ! rotate gf-vectors by rotation matrix
   public gf3_veloc2glob   ! transform velocity-gf from local to global coordinates
   public gf3_cart2curv    ! transform gf from cartesian to curvilinear form
   public gf3_curv2cart    ! transform gf from curvilinear to cartesian form
   public gf3_dq_shift     ! in rolling, interpolate data from the current to the new grid
   public gf3_mirror_x     ! copy values with mirroring in x-direction
   public gf3_print
   public gf3_destroy
   public gf3_test

   interface gf3_dq_shift
      module procedure gf3_dq_shift_arr
      module procedure gf3_dq_shift_scalar
   end interface gf3_dq_shift

   !---------------------------------------------------------------------------------------------------------
   ! data-type for a single element division, i.e. a mask-array defined on a grid

   type :: t_eldiv
      integer, dimension(:), pointer :: el     => NULL()
      integer, dimension(:), pointer :: row1st => NULL()
      integer, dimension(:), pointer :: rowlst => NULL()
      integer                        :: ixmin
      integer                        :: ixmax
      integer                        :: iymin
      integer                        :: iymax
      type(t_grid),          pointer :: grid   => NULL()
   contains
      procedure :: is_defined => eldiv_is_defined

      ! el      element division array (1:ntot): for each element a code Exter, Adhes, Slip or Plast.
      ! row1st  for each row iy of the grid, smallest ix for which el>=Adhes,
      !         i.e. index of the first element inside the contact area.
      ! rowlst  for each row iy of the grid, largest ix for which el>=Adhes,
      !         i.e. index of the last element inside the contact area.
      ! ixmin   overall smallest ix for which elements in the contact area occur.
      ! ixmax   overall largest ix for which elements in the contact area occur.
      ! iymin   overall smallest iy for which elements in the contact area occur.
      ! iymax   overall largest iy for which elements in the contact area occur.
      ! grid    pointer to the grid on which the element division is defined, particularly nx, ny

   end type t_eldiv

   !---------------------------------------------------------------------------------------------------------
   ! data-type for a integer scalar grid function:

   type :: t_igrdfnc1
      integer,      dimension(:), pointer :: val  => NULL()
      type(t_grid),               pointer :: grid => NULL()

      ! val     values array (1:ntot): for each point an integer value
      ! grid    pointer to the grid on which the data are defined, particularly nx, ny

   end type t_igrdfnc1

   !---------------------------------------------------------------------------------------------------------
   ! data-type for a 3-column grid function:
   ! (can be used as scalar grid function as well, via pointers vn and vt)

   type :: t_gridfnc3
      character(len=STRLEN)                 :: name
      real(kind=8), dimension(:,:), pointer :: val   => NULL()
      real(kind=8), dimension(:),   pointer :: vx    => NULL()
      real(kind=8), dimension(:),   pointer :: vy    => NULL()
      real(kind=8), dimension(:),   pointer :: vn    => NULL()
      real(kind=8), dimension(:),   pointer :: vt    => NULL()
      type(t_grid),                 pointer :: grid  => NULL()
      type(t_eldiv),                pointer :: eldiv => NULL()
   contains
      procedure :: is_defined  => gf3_is_defined

      ! name    name of the grid-function, for debugging purposes
      ! val     values array (1:ntot,3): for each grid point a 3-vector with values, associated with
      !         the x,y- and n-coordinate directions
      ! vx      pointer to the values for the x-direction, that is, the first column of the val array
      ! vy      pointer to the values for the y-direction, second column of val array
      ! vn      pointer to the values for the normal direction, third column of val array
      ! vt      pointer to the values for "the tangential direction", such as magnitude of the slip.
      !         Points to same column as vx.
      ! grid    pointer to the grid on which the data are defined, particularly nx, ny
      ! eldiv   optional pointer to an element division (mask array) for subdividing the points

   end type t_gridfnc3

   type :: p_gridfnc3
      type(t_gridfnc3), pointer :: gf => NULL()      ! pointer to a gridfnc3 data-structure
   end type p_gridfnc3

   interface

      module subroutine eldiv_resize(e1, g_new)
      !--purpose: resize element division to the new (matching) grid provided
      !--subroutine arguments:
         type(t_eldiv)           :: e1
         type(t_grid),  target   :: g_new
      end subroutine eldiv_resize

      module subroutine gf3_resize(gf, g_new, defval)
      !--purpose: resize grid function to the new (matching) grid provided
      !--subroutine arguments:
         type(t_gridfnc3)          :: gf
         type(t_grid),    target   :: g_new
         real(kind=8),    optional :: defval
      end subroutine gf3_resize

      module subroutine WrIgs (igs, is_roll, chi)
      !--subroutine arguments:
         type(t_eldiv) :: igs
         logical       :: is_roll
         real(kind=8)  :: chi
      end subroutine wrigs

   end interface

contains

!------------------------------------------------------------------------------------------------------------

function eldiv_is_defined(this)
!--purpose: determine if the element division is 'defined', has memory allocated
   implicit none
!--function return value:
   logical             :: eldiv_is_defined
!--function arguments:
   class(t_eldiv)      :: this

   eldiv_is_defined = (associated(this%el))

end function eldiv_is_defined

!------------------------------------------------------------------------------------------------------------

subroutine eldiv_nullify(eldiv)
!--purpose: initialize/nullify all pointers for an element division
   implicit none
!--subroutine arguments:
   type(t_eldiv)        :: eldiv

   nullify(eldiv%el)
   nullify(eldiv%row1st)
   nullify(eldiv%rowlst)
   nullify(eldiv%grid)

end subroutine eldiv_nullify

!------------------------------------------------------------------------------------------------------------

subroutine eldiv_new(eldiv, grid, nulify)
!--purpose: (re-)allocate space for an element division on the given grid
   implicit none
!--subroutine arguments:
   type(t_grid),         target   :: grid
   type(t_eldiv)                  :: eldiv
   logical, intent(in),  optional :: nulify

   ! optionally nullify the pointers before starting

   if (present(nulify)) then
      if (nulify) call eldiv_nullify(eldiv)
   endif

   ! store pointer to the grid on which this eldiv lives

   eldiv%grid => grid

   ! allocate component-arrays if not done so before, re-allocate at correct size when necessary

   call reallocate_arr(eldiv%el, grid%ntot)
   call reallocate_arr(eldiv%row1st, grid%ny)
   call reallocate_arr(eldiv%rowlst, grid%ny)

   ! (re-)initialize row1st and rowlst pointers and bounding box

   eldiv%row1st(1:grid%ny) = grid%nx
   eldiv%rowlst(1:grid%ny) = 0

   eldiv%ixmin = grid%nx
   eldiv%ixmax = 0
   eldiv%iymin = grid%ny
   eldiv%iymax = 0

end subroutine eldiv_new

!------------------------------------------------------------------------------------------------------------

subroutine eldiv_exter(e1)
!--purpose: mark all elements as "Exterior, might enter contact area"
   implicit none
!--subroutine arguments:
   type(t_eldiv) :: e1
!--local variables:
   integer ii, iy

   if (.not.associated(e1%grid)) &
      call write_log(' eldiv_exter: Internal error: e1%grid not associated')
   if (.not.associated(e1%el)) &
      call write_log(' eldiv_exter: Internal error: e1%el-array not associated')

   do ii = 1, e1%grid%ntot
      e1%el(ii) = Exter
   enddo

   do iy = 1, e1%grid%ny
      e1%row1st(iy) = e1%grid%nx
      e1%rowlst(iy) = 0
   enddo

   e1%ixmin = e1%grid%nx
   e1%ixmax = 0
   e1%iymin = e1%grid%ny
   e1%iymax = 0

end subroutine eldiv_exter

!------------------------------------------------------------------------------------------------------------

subroutine eldiv_copy(e1, e2, imode)
!--purpose: copy contents of first element division to second element division
!          imode=1: normal part only, Slip->Adhesion
!          imode=2: full copy
   implicit none
!--subroutine arguments:
   type(t_eldiv) :: e1, e2
   integer       :: imode
!--local variables:
   integer ii, ny, ntot

   ntot = e1%grid%ntot
   if (.not.associated(e1%el)) &
      call write_log(' eldiv_copy: Internal error: e1%el-array not associated')
   if (.not.associated(e1%row1st)) &
      call write_log(' eldiv_copy: Internal error: e1%row1st-array not associated')
   if (.not.associated(e1%rowlst)) &
      call write_log(' eldiv_copy: Internal error: e1%rowlst-array not associated')
   if (.not.associated(e2%el)) &
      call write_log(' eldiv_copy: Internal error: e2%el-array not associated')
   if (.not.associated(e2%row1st)) &
      call write_log(' eldiv_copy: Internal error: e2%row1st-array not associated')
   if (.not.associated(e2%rowlst)) &
      call write_log(' eldiv_copy: Internal error: e2%rowlst-array not associated')
   if (size(e2%el,1).ne.ntot) then
      write(bufout,*) 'eldiv_copy: Internal error: size e2%el <> ntot', size(e2%el,1), ntot
      call write_log(1, bufout)
   endif

   ! copy element division itself

   if (imode.eq.ikZDIR) then

      ! mode 1: normal part only, replace Slip and Plast by Adhesion

      do ii = 1, ntot
         if (e1%el(ii).ge.Adhes) then
            e2%el(ii) = Adhes
         else
            e2%el(ii) = e1%el(ii)
         endif
      enddo

   elseif (imode.eq.ikALL) then

      ! mode 2: full copy

      do ii = 1, ntot
         e2%el(ii) = e1%el(ii)
      enddo

   else

      write(bufout,*) 'INTERNAL ERROR (eldiv_copy): invalid mode=', imode
      call write_log(1, bufout)
      call abort_run()

   endif

   ! copy row1st and rowlst arrays as well
   ! Note: may need to check size of arrays of e2?

   ny = e1%grid%ny
   e2%row1st(1:ny) = e1%row1st(1:ny)
   e2%rowlst(1:ny) = e1%rowlst(1:ny)

   e2%ixmin = e1%ixmin
   e2%ixmax = e1%ixmax
   e2%iymin = e1%iymin
   e2%iymax = e1%iymax

end subroutine eldiv_copy

!------------------------------------------------------------------------------------------------------------

subroutine eldiv_count(e1, nadh, nslip, nplast, nexter)
!--purpose: count number of elements in adhesion, slip, plasticity and exterior
   implicit none
!--subroutine arguments:
   type(t_eldiv), intent(in)  :: e1
   integer,       intent(out) :: nadh, nslip, nplast, nexter
!--local variables:
   integer ii

   nadh   = 0
   nslip  = 0
   nplast = 0
   nexter = 0
   if (.not.associated(e1%grid) .or. .not.associated(e1%el)) return

   do ii = 1, e1%grid%ntot
      if (e1%el(ii).eq.Adhes) then
         nadh   = nadh   + 1
      elseif (e1%el(ii).eq.Slip) then
         nslip  = nslip  + 1
      elseif (e1%el(ii).eq.Plast) then
         nplast = nplast + 1
      else
         nexter = nexter + 1
      endif
   enddo

end subroutine eldiv_count

!------------------------------------------------------------------------------------------------------------

function eldiv_count_atbnd(e1, rowcol, idebug)
!--Purpose: count number of elements at start/end of rows and/or columns of potential contact area
   implicit none
!--function result:
   integer      :: eldiv_count_atbnd
!--subroutine arguments:
   type(t_eldiv), intent(in)  :: e1
   integer,       intent(in)  :: rowcol         ! check start end of rows (1), columns (-1) or both (0)
   integer,       intent(in)  :: idebug
!--local variables:
   integer      :: ix, iy, ixstep, iystep, ii, nx, ny

   associate(icount => eldiv_count_atbnd)

   icount = 0
   if (.not.associated(e1%grid) .or. .not.associated(e1%el)) return

   nx = e1%grid%nx
   ny = e1%grid%ny

   if (rowcol.ge.0 .and. nx.ge.3) then
      ! include check start/end of rows
      iystep = 1
      else
      ! check start/end of columns only
      iystep = max(1, ny-1)
   endif

   ! loop over rows iy

   do iy = 1, ny, iystep

      if (rowcol.le.0 .and. ny.ge.3 .and. (iy.eq.1 .or. iy.eq.ny)) then
         ! first & last rows: include check start/end of columns
         ixstep = 1
      else
         ! interior rows: check start/end of rows only
         ixstep = max(1, nx-1)
      endif

      do ix = 1, nx, ixstep
         ii = ix + (iy-1) * nx
         if (e1%el(ii).ge.Adhes) then
            if (idebug.ge.5) then
               write(bufout,'(3(a,i3))') ' error at (',ix,',',iy,'): el=',e1%el(ii)
               call write_log(1, bufout)
            endif
            icount = icount + 1
         endif
      enddo
   enddo

   end associate
end function eldiv_count_atbnd

!------------------------------------------------------------------------------------------------------------

subroutine eldiv_cpatches(eldiv, npatch, ysep, idebug)
!--function: fill element division (mask) with sub-patch numbers 1:npatch on the basis of y-values ysep
!            used in gf3_msk_copy
   implicit none
!--subroutine arguments
   type(t_eldiv)             :: eldiv
   integer,      intent(in)  :: npatch, idebug
   real(kind=8), intent(in)  :: ysep(npatch-1)
!--local variables
   integer      :: ip, ii, ix, iy, iy0, iy1
   logical      :: ldone

   associate( mx => eldiv%grid%nx, my => eldiv%grid%ny, y1 => eldiv%grid%y(1), yn => eldiv%grid%y(mx*my))

   if (idebug.ge.2) then
      write(bufout,'(a,i4,2(a,f7.3),a,i2,a,4f7.3)') ' my=',my,', y1=',y1,', yn=', yn, ', npatch=',      &
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

   end associate
end subroutine eldiv_cpatches

!------------------------------------------------------------------------------------------------------------

subroutine eldiv_print(e1, nam)
!--purpose: print debug-output for an element division
   implicit none
!--subroutine arguments:
   type(t_eldiv),    intent(in) :: e1
   character(len=*), intent(in) :: nam
!--local variables:
   integer iy

   call write_log(' printing information on element division ' // trim(nam))
   if (.not.associated(e1%grid)) then
      call write_log(' ..eldiv does not contain a grid-reference')
   else
      write(bufout,'(2(a,i4),a,i6)') ' ..grid has size',e1%grid%nx,' x', e1%grid%ny,', ntot=',e1%grid%ntot
      call write_log(1, bufout)
   endif
   if (.not.associated(e1%el)) then
      call write_log(' ..el-array is not associated')
   else
      write(bufout,'(a,i6)') ' ..el-array has size',size(e1%el)
      call write_log(1, bufout)
   endif
   write(bufout,'(4(a,i4),a)') ' ..encompassing rectangle: [',e1%ixmin,',',e1%ixmax, '] x [',e1%iymin,   &
        ',',e1%iymax,']'
   call write_log(1, bufout)
   if (.not.associated(e1%row1st) .or. .not.associated(e1%rowlst)) then
      call write_log(' ..row1st or rowlst-array is not associated')
   else
      write(bufout,'(2(a,i6))') ' ..row1st/lst-arrays have size',size(e1%row1st), ',', size(e1%rowlst)
      call write_log(1, bufout)
      do iy = 1, min(size(e1%row1st), size(e1%rowlst))
         write(bufout,'(3(a,i4))') '   row',iy,': ix=',e1%row1st(iy),' to',e1%rowlst(iy)
         call write_log(1, bufout)
      enddo
   endif

end subroutine eldiv_print

!------------------------------------------------------------------------------------------------------------

subroutine Areas (eldiv)
!--purpose: update the arrays Row1st and RowLst, which are used in AijPj. These arrays point to the
!          first and last elements of each row in the contact area. if there are no elements in C,
!          first=nx and last=0.
   implicit none
!--subroutine arguments:
   type(t_eldiv) :: eldiv
!--local variables:
   integer nx, ny, ix, iy, i0, first, last

   nx = eldiv%grid%nx
   ny = eldiv%grid%ny

   ! initialize variables that describe the rectangle encompassing the actual contact area

   eldiv%ixmin = nx
   eldiv%ixmax = 1
   eldiv%iymin = ny
   eldiv%iymax = 1

   do iy = 1, ny

      ! compute index of 0th element on row iy

      i0 = (iy-1) * nx

      ! find first interior element of row iy

      first = 1
      do while (eldiv%el(first+i0).le.Exter .and. first.ne.nx)
         first = first + 1
      enddo

      ! find last interior element of row iy

      last = 0
      do ix = first, nx
         if (eldiv%el(ix+i0).ge.Adhes) last = ix
      enddo

      ! store the selection for row iy, note that first=nx, last=0 for empty rows

      eldiv%Row1st(iy) = first
      eldiv%RowLst(iy) = last

      ! update the extent of the encompassing rectangle

      if (first.le.last) then
         eldiv%ixmin = min(eldiv%ixmin, first)
         eldiv%ixmax = max(eldiv%ixmax, last)
         eldiv%iymin = min(eldiv%iymin, iy)
         eldiv%iymax = max(eldiv%iymax, iy)
      endif
   enddo

end subroutine areas

!------------------------------------------------------------------------------------------------------------

subroutine eldiv_destroy(eldiv)
!--purpose: cleanup, deallocate space for an element division
   implicit none
!--subroutine arguments:
   type(t_eldiv)        :: eldiv

   ! de-allocate component arrays

   if (associated(eldiv%el))     deallocate(eldiv%el)
   if (associated(eldiv%row1st)) deallocate(eldiv%row1st)
   if (associated(eldiv%rowlst)) deallocate(eldiv%rowlst)

   ! clear pointers

   nullify(eldiv%el)
   nullify(eldiv%row1st)
   nullify(eldiv%rowlst)
   nullify(eldiv%grid)

end subroutine eldiv_destroy

!------------------------------------------------------------------------------------------------------------

subroutine if1_nullify(if1)
!--purpose: initialize/nullify all pointers for an integer grid function
   implicit none
!--subroutine arguments:
   type(t_igrdfnc1)        :: if1

   nullify(if1%val)
   nullify(if1%grid)

end subroutine if1_nullify

!------------------------------------------------------------------------------------------------------------

subroutine if1_new(if1, grid)
!--purpose: (re-)allocate space for an integer grid function on the given grid
   implicit none
!--subroutine arguments:
   type(t_grid),    target :: grid
   type(t_igrdfnc1)        :: if1

   ! store pointer to grid in grid-function

   if1%grid => grid

   ! allocate values-array if not done so before, re-allocate when needed

   call reallocate_arr(if1%val, grid%ntot)

end subroutine if1_new

!------------------------------------------------------------------------------------------------------------

subroutine if1_destroy(if1)
!--purpose: cleanup, deallocate space for a grid function
   implicit none
!--subroutine arguments:
   type(t_igrdfnc1)        :: if1

   ! de-allocate values-array

   if (associated(if1%val)) deallocate(if1%val)

   ! clear pointers

   nullify(if1%val)
   nullify(if1%grid)

end subroutine if1_destroy

!------------------------------------------------------------------------------------------------------------

subroutine gf3_ikrange(ikarg, ik0, ik1)
!--purpose: Expand the possibly generic selection code "ikarg" for the coordinate directions into
!           a range of concrete coordinate directions "ik0:ik1"
   implicit none
!--subroutine parameters:
   integer, intent(in)  :: ikarg
   integer, intent(out) :: ik0, ik1

   if (ikarg.eq.ikALL) then
      ik0 = 1
      ik1 = 3
   elseif (ikarg.eq.ikTANG) then
      ik0 = ikXDIR      ! assuming XDIR < YDIR
      ik1 = ikYDIR
   elseif (ikarg.ge.1 .and. ikarg.le.3) then
      ik0 = ikarg
      ik1 = ikarg
   else
      ik0 = 1
      ik1 = 0
   endif

end subroutine gf3_ikrange

!------------------------------------------------------------------------------------------------------------

function gf3_is_defined(this)
!--purpose: determine if the gf3 is 'defined', has memory allocated
   implicit none
!--function return value:
   logical                 :: gf3_is_defined
!--function arguments:
   class(t_gridfnc3)       :: this

   gf3_is_defined = (associated(this%val) .and. associated(this%grid))

end function gf3_is_defined

!------------------------------------------------------------------------------------------------------------

subroutine gf3_nullify(gf3)
!--purpose: initialize/nullify all pointers for a grid function.
!           used to signal that there's no memory allocated for the gf.
   implicit none
!--subroutine arguments:
   type(t_gridfnc3)        :: gf3

   gf3%name = ' '
   nullify(gf3%val)
   nullify(gf3%vx)
   nullify(gf3%vy)
   nullify(gf3%vn)
   nullify(gf3%vt)
   nullify(gf3%grid)
   nullify(gf3%eldiv)

end subroutine gf3_nullify

!------------------------------------------------------------------------------------------------------------

subroutine gf3_new(gf3, gf_name, grid, eldiv, nulify, lzero)
!--purpose: (re-)allocate space for a grid function on the given grid
!           optionally nullify the pointers first
   implicit none
!--subroutine arguments:
   character(len=*)                       :: gf_name
   type(t_grid),                   target :: grid
   type(t_gridfnc3)                       :: gf3
   type(t_eldiv),        optional, target :: eldiv
   logical, intent(in),  optional         :: nulify
   logical, intent(in),  optional         :: lzero
!--local variables:
   integer, parameter      :: idebug = 0

   if (idebug.ge.1 .and. gf3%name.eq.' ') then
      write(bufout,*) 'gf3_new: initializing grid-function ',trim(gf_name)
      call write_log(1, bufout)
   endif

   ! optionally nullify the pointers before starting

   if (present(nulify)) then
      if (nulify) call gf3_nullify(gf3)
   endif

   ! store name in grid-function

   gf3%name = gf_name

   ! store pointer to grid in grid-function

   gf3%grid => grid

   ! optionally link the provided element division

   nullify(gf3%eldiv)
   if (present(eldiv)) then
      call gf3_eldiv(gf3, eldiv)
   endif

   ! allocate values-array if not done so before

   call reallocate_arr(gf3%val, grid%ntot, 3)

   ! set pointers to the separate columns of the values array

   gf3%vx => gf3%val(:,ikXDIR)
   gf3%vy => gf3%val(:,ikYDIR)
   gf3%vn => gf3%val(:,ikZDIR)
   gf3%vt => gf3%val(:,ikXDIR)

   ! optionally initialise the arrays with all zeros

   if (present(lzero)) then
      if (lzero) then
         gf3%val(1:grid%ntot,1:3) = 0d0
      endif
   endif

end subroutine gf3_new

!------------------------------------------------------------------------------------------------------------

subroutine gf3_copy_struc(old_gf, new_gf, gf_name, lzero)
!--purpose: copy the structure of an existing grid function to a new grid function
!           nullify & optionally initialize to zero
   implicit none
!--subroutine arguments:
   type(t_gridfnc3)          :: old_gf, new_gf
   character(len=*)          :: gf_name
   logical,         optional :: lzero

   ! assume that "new_gf" has not been used before, else the allocated space is lost (memory leak?)

   call gf3_nullify(new_gf)

   ! create the grid function using the grid of the existing grid function

   call gf3_new(new_gf, gf_name, old_gf%grid)

   ! link the element division of the existing grid function

   call gf3_eldiv(new_gf, old_gf%eldiv)

   ! optionally initialize the data to zero

   if (present(lzero)) then
      if (lzero) call gf3_set(AllElm, 0d0, new_gf, ikALL)
   endif

end subroutine gf3_copy_struc

!------------------------------------------------------------------------------------------------------------

subroutine gf3_copy_full(old_gf, new_gf, grid, eldiv)
!--purpose: copy the structure and values of an existing grid function to a new grid function
!           optionally provide existing copies of grid and eldiv replacing the internal pointers
   implicit none
!--subroutine arguments:
   type(t_gridfnc3)                 :: old_gf, new_gf
   type(t_grid),   optional, target :: grid
   type(t_eldiv),  optional, target :: eldiv

   if (.not.old_gf%is_defined()) then

      call gf3_nullify(new_gf)

   else

      ! copy structure of old_gf to new_gf

      call gf3_copy_struc(old_gf, new_gf, old_gf%name)

      ! replace pointers to grid, eldiv if provided

      if (present(grid)) new_gf%grid => grid

      if (present(eldiv)) new_gf%eldiv => eldiv

      ! copy values from old to new

      call gf3_copy( AllElm, old_gf, new_gf, ikALL )
   endif

end subroutine gf3_copy_full

!------------------------------------------------------------------------------------------------------------

subroutine gf3_eldiv(gf3, eldiv)
!--purpose: add reference to element division to an existing grid function
   implicit none
!--subroutine arguments:
   type(t_eldiv),   target :: eldiv
   type(t_gridfnc3)        :: gf3

   ! store pointer to grid in grid-function

   gf3%eldiv => eldiv

end subroutine gf3_eldiv

!------------------------------------------------------------------------------------------------------------

subroutine gf3_set(iigs, val, f, ikarg)
!--purpose: initialise grid function to a given value
!          ik == coordinate direction(s) of set.
!          iigs: Select elements for the computation: AllElm or AllInt.
   implicit none
!--subroutine arguments:
   type(t_gridfnc3) :: f
   integer          :: ikarg, iigs
   real(kind=8)     :: val
!--local variables:
   integer ik, ik0, ik1, ii

   ! determine range for coordinate directions ik

   call gf3_ikrange(ikarg, ik0, ik1)

   ! set grid function for coordinate directions ik

   do ik = ik0, ik1
      if (iigs.eq.AllElm) then

         call dset (f%grid%ntot, val, f%val(1:,ik), 1)

      elseif (iigs.eq.AllInt) then

         if (.not.associated(f%eldiv)) then
            call write_log(' error with el.div of gf '//trim(f%name))
            call abort_run()
         endif
         if (.not.associated(f%val)) then
            call write_log(' error with val of gf '//trim(f%name))
            call abort_run()
         endif
         do ii = 1, f%grid%ntot
            if (f%eldiv%el(ii).ge.Adhes) then
               f%val(ii,ik) = val
            endif
         enddo
      else
         call write_log(' gf3_set: Internal error: only AllElm or AllInt supported for iigs.')
         call abort_run()
      endif

   enddo

end subroutine gf3_set

!------------------------------------------------------------------------------------------------------------

subroutine gf3_copy_xdir(iigs, n_in, val, f, ikarg)
!--purpose: copy array with values iy along x-direction (same for all ikarg)
!          ik == coordinate direction(s) of set.
!          iigs: Select elements for the computation: AllElm or AllInt.
   implicit none
!--subroutine arguments:
   type(t_gridfnc3) :: f
   integer          :: iigs, n_in, ikarg
   real(kind=8)     :: val(n_in)
!--local variables:
   integer ik, ik0, ik1, nx, ny, ix, iy, ii, iv

   nx = f%grid%nx
   ny = f%grid%ny

   if (n_in.ne.1 .and. n_in.ne.ny) then
      write(bufout,'(a,i4,a)') ' gf3_copy_xdir: error: val should have 1 or ny=',ny,' values.'
      call write_log(1, bufout)
      call abort_run()
   endif
   if (.not.associated(f%val)) then
      call write_log(' error with val of gf '//trim(f%name))
      call abort_run()
   endif
   if (iigs.eq.AllInt .and. .not.associated(f%eldiv)) then
      call write_log(' error with el.div of gf '//trim(f%name))
      call abort_run()
   endif
   if (iigs.ne.AllElm .and. iigs.ne.AllInt) then
      call write_log(' gf3_copy_xdir: Internal error: only AllElm or AllInt supported for iigs.')
      call abort_run()
   endif

   ! determine range for coordinate directions ik

   call gf3_ikrange(ikarg, ik0, ik1)

   iv = 1 ! stays 1 when n_in==1

   ! set grid function for coordinate directions ik

   do ik = ik0, ik1
      do iy = 1, ny
         if (n_in.gt.1) iv = iy

         do ix = 1, nx
            ii = ix + nx * (iy-1)

            if (iigs.eq.AllElm) then

               f%val(ii,ik) = val(iv)

            elseif (iigs.eq.AllInt) then

               if (f%eldiv%el(ii).ge.Adhes) f%val(ii,ik) = val(iv)

            endif
         enddo ! ix
      enddo ! iy
   enddo ! ik

end subroutine gf3_copy_xdir

!------------------------------------------------------------------------------------------------------------

subroutine gf3_copy(iigs, f1, f2, ikarg)
!--purpose: copy values from first grid function to second grid function
!          ik == coordinate direction(s) of copy.
!          iigs: Select elements for the computation: AllElm or AllInt.
   implicit none
!--subroutine arguments:
   type(t_gridfnc3) :: f1, f2
   integer          :: ikarg, iigs
!--local variables:
   integer ik, ik0, ik1, ii

   ! check input arguments

   if (.not.f1%is_defined()) then
      call write_log(' gf3_copy: Internal error: input grid-func not initialized properly.')
      call abort_run()
   elseif (.not.f2%is_defined()) then
      call write_log(' gf3_copy: Internal error: output grid-func not initialized properly.')
      call abort_run()
   elseif (size(f1%val,1).ne.size(f2%val,1)) then
      write(bufout,*) ' gf3_copy: Internal error: incompatible sizes in grid-funcs ', trim(f1%name),    &
          ' and ',trim(f2%name), size(f1%val,1),size(f2%val,1)
      call write_log(1, bufout)
      call abort_run()
   elseif (iigs.eq.AllInt .and. .not.associated(f1%eldiv)) then
      call write_log(' gf3_copy: Internal error: no eldiv, cannot use AllInt for grid-func '//trim(f1%name))
      call abort_run()
   endif

   ! determine range for coordinate directions ik

   call gf3_ikrange(ikarg, ik0, ik1)

   ! copy grid function for coordinate directions ik

   do ik = ik0, ik1
      if (iigs.eq.AllElm) then
         call dcopy (f1%grid%ntot, f1%val(1:,ik),1, f2%val(1:,ik),1)
      elseif (iigs.eq.AllInt) then
         do ii = 1, f1%grid%ntot
            if (f1%eldiv%el(ii).ge.Adhes) then
               f2%val(ii,ik) = f1%val(ii,ik)
            endif
         enddo
      else
         call write_log('gf3_copy: error: only AllElm or AllInt supported for iigs.')
         call abort_run()
      endif
   enddo

end subroutine gf3_copy

!------------------------------------------------------------------------------------------------------------

subroutine gf3_msk_copy(imsk, mask, f1, f2, ikarg)
!--purpose: copy first grid function to second grid function at elements where mask(iel)==imsk
!          ik == coordinate direction(s) of copy.
   implicit none
!--subroutine arguments:
   type(t_gridfnc3) :: f1, f2
   type(t_eldiv)    :: mask
   integer          :: ikarg, imsk
!--local variables:
   integer ik, ik0, ik1, ii

   ! check input arguments

   if (.not.f1%is_defined()) then
      call write_log(' gf3_msk_copy: Internal error: input grid-func not initialized properly.')
      call abort_run()
   elseif (.not.f2%is_defined()) then
      call write_log(' gf3_msk_copy: Internal error: output grid-func not initialized properly.')
      call abort_run()
   elseif (.not.mask%is_defined()) then
      call write_log(' gf3_msk_copy: Internal error: mask not initialized properly.')
      call abort_run()
   elseif (size(f1%val,1).ne.size(f2%val,1)) then
      write(bufout,*) ' gf3_msk_copy: Internal error: incompatible sizes in grid-funcs ',trim(f1%name),   &
          ' and ',trim(f2%name), size(f1%val,1), size(f2%val,1)
      call write_log(1, bufout)
      call abort_run()
   elseif (size(f1%val,1).ne.size(mask%el,1)) then
      write(bufout,*) ' gf3_msk_copy: Internal error: incompatible sizes in grid-func ',trim(f1%name),    &
          ' and mask ', size(f1%val,1), size(mask%el,1)
      call write_log(1, bufout)
      call abort_run()
   endif

   ! determine range for coordinate directions ik

   call gf3_ikrange(ikarg, ik0, ik1)

   ! copy grid function for coordinate directions ik

   do ik = ik0, ik1
      do ii = 1, f1%grid%ntot
         if (mask%el(ii).eq.imsk) then
            f2%val(ii,ik) = f1%val(ii,ik)
         endif
      enddo
   enddo

end subroutine gf3_msk_copy

!------------------------------------------------------------------------------------------------------------

subroutine gf3_scal(iigs, a, f, ikarg)
!--purpose: scale grid function by factor a
!          ik == coordinate direction(s) of copy.
!          iigs: Select elements for the computation: AllElm or AllInt.
   implicit none
!--subroutine arguments:
   type(t_gridfnc3) :: f
   integer          :: ikarg, iigs
   real(kind=8)     :: a
!--local variables:
   integer ik, ik0, ik1, ii

   ! check input arguments

   if (.not.associated(f%val)) then
      call write_log(' gf3_scal: error: values-array not available in grid-func '// trim(f%name))
      call abort_run()
   endif

   ! determine range for coordinate directions ik

   call gf3_ikrange(ikarg, ik0, ik1)

   ! scale grid function for coordinate directions ik

   do ik = ik0, ik1
      if (iigs.eq.AllElm) then
         call dscal (f%grid%ntot, a, f%val(1:,ik),1)
      elseif (iigs.eq.AllInt) then
         do ii = 1, f%grid%ntot
            if (f%eldiv%el(ii).ge.Adhes) then
               f%val(ii,ik) = a * f%val(ii,ik)
            endif
         enddo
      else
         call write_log(' gf3_scal: error: only AllElm or AllInt supported for iigs.')
         call abort_run()
      endif
   enddo

end subroutine gf3_scal

!------------------------------------------------------------------------------------------------------------

subroutine gf3_msk_scal(iigs, mask, a, f, ikarg)
!--purpose: scale grid function by factor 'a' at elements where mask(iel)==iigs
!          ik == coordinate direction(s) of copy.
   implicit none
!--subroutine arguments:
   type(t_gridfnc3) :: f
   type(t_eldiv)    :: mask
   integer          :: ikarg, iigs
   real(kind=8)     :: a
!--local variables:
   integer ik, ik0, ik1, ii

   ! check input arguments

   if (.not.associated(f%val)) then
      call write_log(' gf3_scal: error: values-array not available in grid-func '// trim(f%name))
      call abort_run()
   elseif (size(f%val,1).ne.size(mask%el,1)) then
      write(bufout,*) ' gf3_scal: Internal error: incompatible sizes in grid-func ', trim(f%name),       &
          ' and mask ', size(f%val,1), size(mask%el,1)
      call write_log(1, bufout)
      call abort_run()
   endif

   ! determine range for coordinate directions ik

   call gf3_ikrange(ikarg, ik0, ik1)

   ! scale grid function for coordinate directions ik

   do ik = ik0, ik1
      do ii = 1, f%grid%ntot
         if (mask%el(ii).eq.iigs) then
            f%val(ii,ik) = a * f%val(ii,ik)
         endif
      enddo
   enddo

end subroutine gf3_msk_scal

!------------------------------------------------------------------------------------------------------------

subroutine gf3_axpy(iigs, a, f1, f2, ikarg)
!--purpose: compute axpy of grid functions: f2 := f2 + a f1
!          ik == coordinate direction(s) of copy.
!          iigs: Select elements for the computation: AllElm or AllInt.
   implicit none
!--subroutine arguments:
   type(t_gridfnc3) :: f1, f2
   integer          :: iigs, ikarg
   real(kind=8)     :: a
!--local variables:
   integer ii, ik, ik0, ik1

   ! determine range for coordinate directions ik

   call gf3_ikrange(ikarg, ik0, ik1)

   ! copy grid function for coordinate directions ik

   do ik = ik0, ik1
      if (iigs.eq.AllElm) then

         call daxpy (f1%grid%ntot, a, f1%val(1:,ik),1, f2%val(1:,ik),1)

      elseif (iigs.eq.AllInt) then

         do ii = 1, f1%grid%ntot
            if (f2%eldiv%el(ii).ge.Adhes) then
               f2%val(ii,ik) = f2%val(ii,ik) + a * f1%val(ii,ik)
            endif
         enddo
      else
         call write_log(' gf3_axpy: error: only AllElm or AllInt supported for iigs.')
         call abort_run()
      endif

   enddo

end subroutine gf3_axpy

!------------------------------------------------------------------------------------------------------------

function gf3_sum(iigs, f, ikarg)
!--purpose: compute the sum of a grid-function
!          ik == coordinate direction(s) of sum.
!          iigs: Select elements for the computation: AllElm or AllInt.
   implicit none
!--function result:
   real(kind=8) gf3_sum
!--subroutine arguments:
   type(t_gridfnc3) :: f
   integer          :: ikarg, iigs
!--local variables:
   integer ik, ik0, ik1, ii

   ! determine range for coordinate directions ik

   call gf3_ikrange(ikarg, ik0, ik1)

   ! compute sum using loop over coordinate directions ik

   gf3_sum = 0d0
   do ik = ik0, ik1
      if (iigs.eq.AllElm) then
         gf3_sum = gf3_sum + dsum(f%grid%ntot, f%val(1:,ik), 1)
      elseif (iigs.eq.AllInt) then
         do ii = 1, f%grid%ntot
            if (f%eldiv%el(ii).ge.Adhes) then
               gf3_sum = gf3_sum + f%val(ii,ik)
            endif
         enddo
      else
         call write_log(' gf3_sum: error: only AllElm or AllInt supported for iigs.')
         call abort_run()
      endif
   enddo

end function gf3_sum

!------------------------------------------------------------------------------------------------------------

subroutine gf3_proj_avg(iigs, f, ikarg)
!--purpose: compute subtract the average from a grid-function
!          (separate averages per coordinate direction)
!          ik == coordinate direction(s).
!          iigs: Select elements for the computation: AllElm or AllInt.
   implicit none
!--subroutine arguments:
   type(t_gridfnc3) :: f
   integer          :: ikarg, iigs
!--local variables:
   integer      :: nval, ik, ik0, ik1, ii, ntot
   real(kind=8) :: avg

   ! determine range for coordinate directions ik

   call gf3_ikrange(ikarg, ik0, ik1)

   ! perform operation in loop over coordinate directions ik

   ntot = f%grid%ntot
   do ik = ik0, ik1
      if (iigs.eq.AllElm) then
         avg = gf3_sum(iigs, f, ik) / real(ntot)
         do ii = 1, ntot
            f%val(ii,ik) = f%val(ii,ik) - avg
         enddo
      elseif (iigs.eq.AllInt) then
         avg  = 0d0
         nval = 0
         do ii = 1, ntot
            if (f%eldiv%el(ii).ge.Adhes) then
               avg = avg + f%val(ii,ik)
               nval = nval + 1
            endif
         enddo
         avg = avg / real(max(1,nval))
         do ii = 1, ntot
            if (f%eldiv%el(ii).ge.Adhes) f%val(ii,ik) = f%val(ii,ik) - avg
         enddo
      else
         call write_log(' gf3_proj_avg: error: only AllElm or AllInt supported for iigs.')
         call abort_run()
      endif
   enddo

end subroutine gf3_proj_avg

!------------------------------------------------------------------------------------------------------------

function gf3_dot(iigs, f1, f2, ikarg)
!--purpose: compute the inner product of two grid-functions
!          ik == coordinate direction(s) of inner product.
!          iigs: Select elements for the computation: AllElm or AllInt.
   implicit none
!--function result:
   real(kind=8) gf3_dot
!--subroutine arguments:
   type(t_gridfnc3) :: f1, f2
   integer          :: ikarg, iigs
!--local variables:
   integer ik, ik0, ik1, ii, ntot

   ! determine range for coordinate directions ik

   call gf3_ikrange(ikarg, ik0, ik1)

   ! compute inner product using loop over coordinate directions ik

   ntot = f1%grid%ntot
   gf3_dot = 0d0
   do ik = ik0, ik1
      if (iigs.eq.AllElm) then
         gf3_dot = gf3_dot + ddot(ntot, f1%val(1:,ik),1, f2%val(1:,ik),1)
      elseif (iigs.eq.AllInt) then
         do ii = 1, ntot
            if (f1%eldiv%el(ii).ge.Adhes) then
               gf3_dot = gf3_dot + f1%val(ii,ik) * f2%val(ii,ik)
            endif
         enddo
      elseif (iigs.ge.Exter .and. iigs.le.Plast) then
         do ii = 1, ntot
            if (f1%eldiv%el(ii).eq.iigs) then
               gf3_dot = gf3_dot + f1%val(ii,ik) * f2%val(ii,ik)
               ! if (ik.eq.2) then
               !     write(bufout,'(3(a,f7.3),a)') ' dot +=',f1%val(ii,ik)*f2%val(ii,ik),&
               !          ' (=', f1%val(ii,ik),' *',f2%val(ii,ik),')'
               !     call write_log(1, bufout)
               ! endif
            endif
         enddo
      else
         write(bufout,*) 'gf3_dot: error: incorrect value for iigs=', iigs
         call write_log(1, bufout)
         call abort_run()
      endif
   enddo

end function gf3_dot

!------------------------------------------------------------------------------------------------------------

function gf3_nrm2(iigs, f, ikarg)
!--purpose: compute the 2-norm of a grid-function
!          ik == coordinate direction(s) of norm.
!          iigs: Select elements for the computation: AllElm or AllInt.
   implicit none
!--function result:
   real(kind=8) gf3_nrm2
!--subroutine arguments:
   type(t_gridfnc3) :: f
   integer          :: ikarg, iigs
!--local variables:
   integer ik, ik0, ik1, ii, ntot

   ! determine range for coordinate directions ik

   call gf3_ikrange(ikarg, ik0, ik1)

   ! compute inner product using loop over coordinate directions ik

   ntot = f%grid%ntot
   gf3_nrm2 = 0d0
   do ik = ik0, ik1
      if (iigs.eq.AllElm) then
         gf3_nrm2 = gf3_nrm2 + dnrm2(ntot, f%val(1:,ik),1) **2
      elseif (iigs.eq.AllInt) then
         do ii = 1, ntot
            if (f%eldiv%el(ii).ge.Adhes) then
               gf3_nrm2 = gf3_nrm2 + f%val(ii,ik)**2
            endif
         enddo
      else
         call write_log(' gf3_nrm2: error: only AllElm or AllInt supported for iigs.')
         call abort_run()
      endif
   enddo
   gf3_nrm2 = dsqrt(gf3_nrm2)

end function gf3_nrm2

!------------------------------------------------------------------------------------------------------------

function gf3_rms(iigs, f, ikarg)
!--purpose: compute the root-mean-square value of a grid-function
!          ik == coordinate direction(s) of norm.
!          iigs: Select elements for the computation: AllElm or AllInt.
   implicit none
!--function result:
   real(kind=8) gf3_rms
!--subroutine arguments:
   type(t_gridfnc3) :: f
   integer          :: ikarg, iigs
!--local variables:
   integer ik, ik0, ik1, ii, icount, ntot

   ! determine range for coordinate directions ik

   call gf3_ikrange(ikarg, ik0, ik1)

   ! compute inner product using loop over coordinate directions ik

   ntot = f%grid%ntot
   gf3_rms = 0d0
   icount = 0
   do ik = ik0, ik1
      if (iigs.eq.AllElm) then
         gf3_rms = gf3_rms + dnrm2(ntot, f%val(1:,ik),1) **2
         icount = icount + ntot
      elseif (iigs.eq.AllInt) then
         do ii = 1, ntot
            if (f%eldiv%el(ii).ge.Adhes) then
               gf3_rms = gf3_rms + f%val(ii,ik)**2
               icount = icount + 1
            endif
         enddo
      else
         call write_log(' gf3_rms: error: only AllElm or AllInt supported for iigs.')
         call abort_run()
      endif
   enddo
   gf3_rms = dsqrt(gf3_rms/max(1,icount))

end function gf3_rms

!------------------------------------------------------------------------------------------------------------

function gf3_min(iigs, f, ikarg)
!--purpose: determine the largest element of a grid-function
!          min = x_i : { forall j : x_j <= x_i }
!          ik == coordinate direction(s) of operation.
!          iigs: Select elements for the computation: AllElm or AllInt.
   implicit none
!--function result:
   real(kind=8) gf3_min
!--subroutine arguments:
   type(t_gridfnc3) :: f
   integer          :: ikarg, iigs
!--local variables:
   integer ii, ik, ik0, ik1

   ! determine range for coordinate directions ik

   call gf3_ikrange(ikarg, ik0, ik1)

   ! update the overall minimum in a loop over coordinate directions ik,
   ! compute the minimum value per direction using BLAS function idmin.

   gf3_min = 1d20
   do ik = ik0, ik1
      if (iigs.eq.AllElm) then
         ii = idmin(f%grid%ntot, f%val(1:,ik),1)
         if (f%val(ii,ik).lt.gf3_min) gf3_min = f%val(ii,ik)
      elseif (iigs.eq.AllInt) then
         do ii = 1, f%grid%ntot
            if (f%eldiv%el(ii).ge.Adhes) then
               if (f%val(ii,ik).lt.gf3_min) gf3_min = f%val(ii,ik)
            endif
         enddo
      else
         call write_log(' gf3_min: error: only AllElm or AllInt supported for iigs.')
         call abort_run()
      endif
   enddo

end function gf3_min

!------------------------------------------------------------------------------------------------------------

function gf3_max(iigs, f, ikarg)
!--purpose: determine the largest element of a grid-function
!          max = x_i : { forall j : x_j <= x_i }
!          ik == coordinate direction(s) of operation.
!          iigs: Select elements for the computation: AllElm or AllInt.
   implicit none
!--function result:
   real(kind=8) gf3_max
!--subroutine arguments:
   type(t_gridfnc3) :: f
   integer          :: ikarg, iigs
!--local variables:
   integer ii, ik, ik0, ik1

   ! determine range for coordinate directions ik

   call gf3_ikrange(ikarg, ik0, ik1)

   ! update the overall maximum in a loop over coordinate directions ik,
   ! compute the maximum value per direction using BLAS function idmax.

   gf3_max = -1d20
   do ik = ik0, ik1
      if (iigs.eq.AllElm) then
         ii = idmax(f%grid%ntot, f%val(1:,ik),1)
         if (f%val(ii,ik).gt.gf3_max) gf3_max = f%val(ii,ik)
      elseif (iigs.eq.AllInt) then
         do ii = 1, f%grid%ntot
            if (f%eldiv%el(ii).ge.Adhes) then
               if (f%val(ii,ik).gt.gf3_max) gf3_max = f%val(ii,ik)
            endif
         enddo
      else
         call write_log(' gf3_max: error: only AllElm or AllInt supported for iigs.')
         call abort_run()
      endif
   enddo

end function gf3_max

!------------------------------------------------------------------------------------------------------------

function gf3_absmax(iigs, f, ikarg)
!--purpose: determine the in absolute sense largest element of a grid-function
!          (can be negative):  absmax = x_i : { forall j : |x_j| <= |x_i| }
!          ik == coordinate direction(s) of operation.
   implicit none
!--function result:
   real(kind=8) gf3_absmax
!--subroutine arguments:
   type(t_gridfnc3) :: f
   integer          :: ikarg, iigs
!--local variables:
   integer ii, ik, ik0, ik1

   ! determine range for coordinate directions ik

   call gf3_ikrange(ikarg, ik0, ik1)

   ! update the overall maximum in a loop over coordinate directions ik,
   ! compute the maximum value per direction using BLAS function idamax.

   gf3_absmax = 0d0
   do ik = ik0, ik1
      if (iigs.eq.AllElm) then
         ii = idamax(f%grid%ntot, f%val(1:,ik),1)
         if (abs(f%val(ii,ik)).gt.abs(gf3_absmax)) then
            gf3_absmax = f%val(ii,ik)
         endif
      elseif (iigs.eq.AllInt) then
         do ii = 1, f%grid%ntot
            if (f%eldiv%el(ii).ge.Adhes) then
               if (abs(f%val(ii,ik)).gt.gf3_absmax) gf3_absmax = f%val(ii,ik)
            endif
         enddo
      else
         call write_log(' gf3_absmax: error: only AllElm or AllInt supported for iigs.')
         call abort_run()
      endif
   enddo

end function gf3_absmax

!------------------------------------------------------------------------------------------------------------

function gf3_maxabs(iigs, f, ikarg)
!--purpose: compute the max-norm of a grid-function
!          Note: this value is always positive.
!          ik == coordinate direction(s) of norm.
!          iigs: Select elements for the computation: AllElm or AllInt.
   implicit none
!--function result:
   real(kind=8) gf3_maxabs
!--subroutine arguments:
   type(t_gridfnc3) :: f
   integer          :: ikarg, iigs

   gf3_maxabs = abs(gf3_absmax(iigs, f, ikarg))

end function gf3_maxabs

!------------------------------------------------------------------------------------------------------------

subroutine gf3_rotate(f, rot)
!--function: rotate grid function values [vx,vy,vz] by rotation matrix, fnew = R * forig
   implicit none
!--subroutine arguments
   type(t_gridfnc3)         :: f            ! input/output grid function
   type(t_rotmat)           :: rot          ! matrix with new orientation of original unit vectors
!--local variables
   integer      :: ii
   real(kind=8) :: vxtmp, vytmp, vztmp

   do ii = 1, f%grid%ntot
      vxtmp = f%vx(ii)
      vytmp = f%vy(ii)
      vztmp = f%vn(ii)
      f%vx(ii) = rot%r(1) * vxtmp + rot%r(4) * vytmp + rot%r(7) * vztmp
      f%vy(ii) = rot%r(2) * vxtmp + rot%r(5) * vytmp + rot%r(8) * vztmp
      f%vn(ii) = rot%r(3) * vxtmp + rot%r(6) * vytmp + rot%r(9) * vztmp
   end do

end subroutine gf3_rotate

!------------------------------------------------------------------------------------------------------------

subroutine gf3_veloc2glob(mloc, tvel_mloc, rvel_mloc, gsurf, veloc)
!--function: convert velocity grid function 'veloc' [vx,vy,vz] from local to global coordinates
   implicit none
!--subroutine arguments
   type(t_gridfnc3)         :: veloc        ! input/output grid function
   type(t_marker)           :: mloc         ! position and orientation of local origin in global coordinates
   type(t_vec)              :: tvel_mloc    ! velocity of local origin in global coordinates
   type(t_vec)              :: rvel_mloc    ! angular velocity vector of local origin in global coordinates
   type(t_grid)             :: gsurf        ! grid/surface in global coordinates
!--local variables
   integer                  :: ntot, ii
   type(t_vec)              :: tmpvel

   if (veloc%grid%ntot.ne.gsurf%ntot) then
      call write_log(' Internal error: gf3_veloc2glob: veloc-grid doesnt match gsurf')
   elseif (.not.associated(gsurf%coor)) then
      call write_log(' Internal error: surface doesn''t have a valid coor-array')
   else
      ntot = gsurf%ntot

      ! rotate input velocities to global coordinate system: veloc = R * veloc

      call gf3_rotate(veloc, mloc%rot)

      ! add contribution of translational velocity, veloc = veloc + velloc 

      do ii = 1, ntot
         veloc%val(ii,1:3) = veloc%val(ii,1:3) + tvel_mloc%v(1:3)
      enddo

      ! add contribution of angular velocity, veloc = veloc + angloc x (x_p - o_m)

      do ii = 1, ntot
         tmpvel = rvel_mloc .cross. (vec( gsurf%coor(ii,1:3) ) - mloc%o)
         veloc%val(ii,1:3) = veloc%val(ii,1:3) + tmpvel%v(1:3)
      enddo
   endif
end subroutine gf3_veloc2glob

!------------------------------------------------------------------------------------------------------------

subroutine gf3_cart2curv(gf, nref, cnrm, irot, idebug, l_inv_arg)
!--purpose: rotate 3-vectors in gf defined in (cartesian: constant) local reference orientation nref
!           (curvilinear: varying) global reference orientation cnrm by rotation about irot-axis
   implicit none
!--subroutine arguments:
   type(t_gridfnc3)      :: gf, cnrm
   type(t_vec)           :: nref
   integer               :: irot, idebug
   logical, optional     :: l_inv_arg
!--local variables:
   integer               :: ii_debug = 1
   logical               :: l_inverse
   integer               :: mx, my, nx, ny, ix, iy, ii, jx, jy, jj
   real(kind=8)          :: nref_x, nref_y, nref_z, fac, cs, sn, tmp

   if (idebug.ge.4) call write_log(' --- start subroutine gf3_cart2curv ---')

   ! invert rotation?

   l_inverse = .false.
   if (present(l_inv_arg)) l_inverse = l_inv_arg

   ! surface normals define orientation only if a rotation vector is given:
   ! irot=1: rvec = [1,0,0] (roll); irot=2: rvec=[0,1,0] (pitch); irot=3: rvec=[0,0,1] (yaw)

   if (irot.lt.1 .or. irot.gt.2) then
      write(bufout,'(a,i4,a)') ' INTERNAL ERROR (gf3_cart2curv): rotation direction must be 1 or 2 (', &
                irot,')'
      call write_log(1, bufout)
      call abort_run()
   endif

   mx = gf%grid%nx
   my = gf%grid%ny
   nx = cnrm%grid%nx
   ny = cnrm%grid%ny

   ! surface normals must be compatible with grid-function's grid

   if ( (nx.ne.1 .and. nx.ne.mx) .or. (ny.ne.1 .and. ny.ne.my)) then
      write(bufout,'(4(a,i4))') ' INTERNAL ERROR (gf3_cart2curv): incompatible grids',mx,' x',my,' and', &
             nx,' x',ny
      call write_log(1, bufout)
      call abort_run()
   endif

   ! loop over all elements (ix,iy) of the grid-function

   nref_x = nref%x()
   nref_y = nref%y()
   nref_z = nref%z()
   if (idebug.ge.4) then
      write(bufout,'(3(a,f8.3),a)') ' ref.normal=[',nref_x, ',', nref_y, ',', nref_z, ']'
      call write_log(1, bufout)
   endif

   do iy = 1, my
      do ix = 1, mx
         ii  = ix + (iy-1) * mx

         ! get corresponding element (jx,jy) of surface normals

         jx = ix
         if (nx.le.1) jx = 1
         jy = iy
         if (ny.le.1) jy = 1
         jj = jx + (jy-1) * nx

         ! compute angle & rotation

         if (irot.eq.1) then

            if (idebug.ge.5 .and. ii.eq.ii_debug) then
               write(bufout,'(4(a,i4),3(a,f8.3),a)') ' element',ii,' (',ix,',',iy,'): jj=',jj,          &
                        ', cnrm=[', cnrm%vx(jj), ',', cnrm%vy(jj), ',', cnrm%vn(jj), ']'
               call write_log(1, bufout)
            endif

            ! rotate gf about x-axis (roll),  using normals [*,ny,nz]

            cs = cnrm%vy(jj) * nref_y + cnrm%vn(jj) * nref_z
            sn = cnrm%vy(jj) * nref_z - cnrm%vn(jj) * nref_y
            fac = 1d0 / sqrt(cs**2 + sn**2)
            cs = cs * fac
            sn = sn * fac

            if (idebug.ge.7 .and. ii.eq.ii_debug) then
               write(bufout,'(3(a,i4),2(a,f8.3))') ' element',ii,' (',ix,',',iy,'): cs=',cs,', sn=',sn
               call write_log(1, bufout)
            endif
            if (l_inverse) sn = -sn

            tmp = gf%vy(ii)
            gf%vy(ii) = cs * gf%vy(ii) - sn * gf%vn(ii)
            gf%vn(ii) = sn * tmp       + cs * gf%vn(ii)

         elseif (irot.eq.2) then

            ! rotate gf about y-axis (pitch), using normals [nx,*,nz]

            cs = cnrm%vx(jj) * nref_x + cnrm%vn(jj) * nref_z
            sn = cnrm%vn(jj) * nref_x - cnrm%vx(jj) * nref_z
            fac = 1d0 / (cs**2 + sn**2)
            cs = cs * fac
            sn = sn * fac
            if (l_inverse) sn = -sn

            tmp = gf%vx(ii)
            gf%vx(ii) =   cs * gf%vx(ii) + sn * gf%vn(ii)
            gf%vn(ii) = - sn * tmp       + cs * gf%vn(ii)

         elseif (irot.eq.3) then

            ! rotate gf about z-axis (yaw),   using normals [nx,ny,*]
            call write_log(' irot=3: not yet implemented')
            call abort_run()

         endif

      enddo
   enddo

   if (idebug.ge.4) call write_log(' --- end subroutine gf3_cart2curv ---')
end subroutine gf3_cart2curv

!------------------------------------------------------------------------------------------------------------

subroutine gf3_curv2cart(gf, nref, cnrm, irot, idebug)
!--purpose: rotate 3-vectors in gf defined in curvilinear/conformal coordinates as indicated by surface
!           normals cnrm+irot given in global coordinates to local cartesian coordinates (mloc)
   implicit none
!--subroutine arguments:
   type(t_gridfnc3)      :: gf, cnrm
   type(t_vec)           :: nref
   integer               :: irot, idebug

   ! curv2cart == inverse of cart2curv rotation

   call gf3_cart2curv(gf, nref, cnrm, irot, idebug, l_inv_arg=.true.)

end subroutine gf3_curv2cart

!------------------------------------------------------------------------------------------------------------

subroutine gf3_dq_shift_arr(chi, dq, numext, extval, fin, fout, ikarg)
!--purpose: in rolling, interpolate data from the current to the new grid
!           ik:   select coordinate direction(s)
   implicit none
!--subroutine arguments:
   type(t_gridfnc3), intent(in)    :: fin
   integer,          intent(in)    :: numext, ikarg
   real(kind=8),     intent(in)    :: chi, dq, extval(numext)
   type(t_gridfnc3), intent(inout) :: fout
!--local variables:
   real(kind=8), parameter :: pi     = 4d0*atan(1d0)
   integer          :: ik0, ik1, ik, nx, ny, iext, ix, iy, iin, iout, kx, ky
   real(kind=8)     :: dx, dy, facx, facy, f00, f01, f10, f11

   nx = fin%grid%nx
   ny = fin%grid%ny
   dx = fin%grid%dx
   dy = fin%grid%dy

   if (fout%grid%nx.ne.nx .or. fout%grid%ny.ne.ny) then
      write(bufout,'(a,4(i4,a))') ' gf3_dq_shift: error: grids differ. In:',nx,' x',ny,          &
                ', out:', fout%grid%nx,' x',fout%grid%ny
      call write_log(1, bufout)
      call abort_run()
   endif
   if (numext.ne.1 .and. numext.ne.ny) then
      write(bufout,'(a,i4,a)') ' gf3_dq_shift: error: there should be 1 or ny=',ny,' exterior values.'
      call write_log(1, bufout)
      call abort_run()
   endif

   call gf3_ikrange(ikarg, ik0, ik1)

   ! kx = floor(dqx/dx) = number of whole elements shifted
   ! fout(ix,iy) = 
   !      (1-(dqy/dy-ky)) * { (1-(dqx/dx-kx)) * fin(ix+kx,iy  +ky) + (dqx/dx-kx) * fin(ix+kx+1,iy  +ky) }
   !    +    (dqy/dy-ky)  * { (1-(dqx/dx-kx)) * fin(ix+kx,iy+1+ky) + (dqx/dx-kx) * fin(ix+kx+1,iy+1+ky) }

   if (abs(chi).le.0.01d0 .or. abs(chi-pi).le.0.01d0) then

      ! special code for chi = 0 or 180deg : pure shift in x-direction (ky=0, facy=0)

      kx   = int(dq*cos(chi)/dx)
      facx =     dq*cos(chi)/dx - real(kx)
      iext = 1

      do ik = ik0, ik1
         do iy = 1, ny
            if (numext.gt.1) iext = iy
            do ix = 1, nx

               ! get values of input cells

               iin = ix+kx + (iy-1)*nx
               if (ix+kx  .le.0 .or. ix+kx  .gt.nx) then
                  f00 = extval(iext)
               else
                  f00 = fin%val(iin  ,ik)
               endif
               if (ix+kx+1.le.0 .or. ix+kx+1.gt.nx) then
                  f10 = extval(iext)
               else
                  f10 = fin%val(iin+1,ik)
               endif

               ! compute value for output cell

               iout = ix + (iy-1)*nx
               fout%val(iout,ik) = (1d0-facx)*f00 + facx*f10
            enddo
         enddo
      enddo

   else

      ! normal code for general chi

      kx   = int(dq*cos(chi)/dx)
      facx =     dq*cos(chi)/dx - real(kx)
      ky   = int(dq*sin(chi)/dy)
      facy =     dq*sin(chi)/dy - real(ky)
      iext = 1

      do ik = ik0, ik1
         do iy = 1, ny
            if (numext.gt.1) iext = iy
            do ix = 1, nx

               ! get values of input cells

               iin = ix+kx + (iy+ky-1)*nx
               if (ix+kx  .le.0 .or. ix+kx  .gt.nx .or. iy+ky  .le.0 .or. iy+ky  .gt.ny) then
                  f00 = extval(iext)
               else
                  f00 = fin%val(iin     ,ik)
               endif
               if (ix+kx+1.le.0 .or. ix+kx+1.gt.nx .or. iy+ky  .le.0 .or. iy+ky  .gt.ny) then
                  f10 = extval(iext)
               else
                  f10 = fin%val(iin+1   ,ik)
               endif
               if (ix+kx  .le.0 .or. ix+kx  .gt.nx .or. iy+ky  .le.0 .or. iy+ky  .gt.ny) then
                  f01 = extval(iext)
               else
                  f01 = fin%val(iin+  nx,ik)
               endif
               if (ix+kx+1.le.0 .or. ix+kx+1.gt.nx .or. iy+ky+1.le.0 .or. iy+ky+1.gt.ny) then
                  f11 = extval(iext)
               else
                  f11 = fin%val(iin+1+nx,ik)
               endif

               ! compute value for output cell

               iout = ix + (iy-1)*nx
               fout%val(iout,ik) =   (1d0-facy) * ( (1d0-facx)*f00 + facx*f10 )                     &
                                   +      facy  * ( (1d0-facx)*f01 + facx*f11 )
            enddo
         enddo
      enddo

   endif ! chi=0 or pi

end subroutine gf3_dq_shift_arr

!------------------------------------------------------------------------------------------------------------

subroutine gf3_dq_shift_scalar(chi, dq, extval, fin, fout, ikarg)
!--purpose: in rolling, interpolate data from the current to the new grid
!           ik:   select coordinate direction(s)
   implicit none
!--subroutine arguments:
   type(t_gridfnc3), intent(in)    :: fin
   integer,          intent(in)    :: ikarg
   real(kind=8),     intent(in)    :: chi, dq, extval
   type(t_gridfnc3), intent(inout) :: fout
!--local variables:
   real(kind=8)     :: extarr(1)

   extarr(1) = extval
   call gf3_dq_shift_arr(chi, dq, 1, extarr, fin, fout, ikarg)

end subroutine gf3_dq_shift_scalar

!------------------------------------------------------------------------------------------------------------

subroutine gf3_mirror_x(ikeep, fin, fout, ikarg, idebug)
!--function: copy fin to fout with mirroring in x-direction
!            ikeep<0 : fout(x<0), fout(x>0) = fin(x<0)
!            ikeep>0 : fout(x<0), fout(x>0) = fin(x>0)
   implicit none
!--subroutine arguments
   type(t_gridfnc3)             :: fin, fout
   integer,          intent(in) :: ikeep, ikarg, idebug
!--local variables
   real(kind=8), parameter :: small = 1d-9
   integer                 :: mx, my, ix, iy, ik, ik0, ik1, nval,                                       &
                              ix_negn, ix_neg1, ix_pos1, ix_posn, ix_neg, ix_pos

   if (ikeep.eq.0) then
      call write_log(' INTERNAL ERROR(mirror_x): ikeep=0')
      call abort_run()
   endif

   mx = fin%grid%nx
   my = fin%grid%ny
   call gf3_ikrange(ikarg, ik0, ik1)

   ! determine last negative ix_neg1 and first positive ix_pos1 in x-values

   ix_pos1 = 0
   do ix = 1, mx
      if (fin%grid%x(ix).lt.-small                   ) ix_neg1 = ix
      if (fin%grid%x(ix).gt. small .and. ix_pos1.eq.0) ix_pos1 = ix
   enddo

   ! determine first negative ix_negn and last positive ix_posn in x-values

   nval = min(ix_neg1, mx-ix_pos1+1)
   ix_negn = ix_neg1 - nval + 1
   ix_posn = ix_pos1 + nval - 1

   if (idebug.ge.2) then
      if (ikeep.lt.0) then
         write(bufout,'(4(a,i3),a)') ' mirror_x: copying [', ix_negn,',',ix_neg1,'] to [',ix_pos1,',',  &
                ix_posn,']'
      else
         write(bufout,'(4(a,i3),a)') ' mirror_x: copying [', ix_pos1,',',ix_posn,'] to [',ix_negn,',',  &
                ix_neg1,']'
      endif
      call write_log(1, bufout)
   endif

   if (idebug.ge.-1 .and. fin%grid%x(ix_pos1)+fin%grid%x(ix_neg1).gt.small) then
      write(bufout,'(2(a,i3,a,f12.6),a)') ' WARNING: input grid is a-symmetric.  x(',ix_neg1,')=',      &
             fin%grid%x(ix_neg1),', x(',ix_pos1,')=',fin%grid%x(ix_pos1),'.'
      call write_log(1, bufout)
   endif

   if (ikeep.lt.0) then

      ! copy values from fin to fout, mirroring negative to positive part  fout(x) = fout(-x), x>0

      do ik = ik0, ik1
         do iy = 1, my
            do ix = 1, mx
               if (ix.lt.ix_pos1) then
                  fout%val(ix,ik) = fin%val(ix,ik)
               elseif (ix.le.ix_posn) then
                  ix_neg = ix_neg1 + ix_pos1 - ix
                  fout%val(ix,ik) = fin%val(ix_neg,ik)
               else
                  fout%val(ix,ik) = 0d0
               endif
            enddo
         enddo
      enddo

!     do ival = 1, nval
!        ix_neg = ix_negn - 1 + ival
!        ix_pos = ix_posn + 1 - ival
!        if (idebug.ge.5) then
!           write(bufout,'(2(a,i3))') ' copying ix=', ix_neg,' to', ix_pos
!           call write_log(1, bufout)
!        endif

!        do iy = 1, my
!           ii_neg = ix_neg + (iy-1) * mx
!           ii_pos = ix_pos + (iy-1) * mx
!           do ik = ik0, ik1
!              fout%val(ii_neg,ik) = fin%val(ii_neg,ik)
!              fout%val(ii_pos,ik) = fin%val(ii_neg,ik)
!           enddo
!        enddo
!     enddo

   else

      ! copy values from fin to fout, mirroring positive to negative part  fout(-x) = fout(x), x>0

      do ik = ik0, ik1
         do iy = 1, my
            do ix = 1, mx
               if (ix.lt.ix_negn) then
                  fout%val(ix,ik) = 0d0
               elseif (ix.le.ix_neg1) then
                  ix_pos = ix_pos1 + ix_neg1 - ix
                  fout%val(ix,ik) = fin%val(ix_pos,ik)
               else
                  fout%val(ix,ik) = fin%val(ix,ik)
               endif
            enddo
         enddo
      enddo

!     do ival = 1, nval
!        ix_pos = ix_pos1 - 1 + ival
!        ix_neg = ix_neg1 + 1 - ival
!        if (idebug.ge.5) then
!           write(bufout,'(2(a,i3))') ' copying ix=', ix_pos,' to', ix_neg
!           call write_log(1, bufout)
!        endif

!        do iy = 1, my
!           ii_neg = ix_neg + (iy-1) * mx
!           ii_pos = ix_pos + (iy-1) * mx
!           do ik = ik0, ik1
!              fout%val(ii_neg,ik) = fin%val(ii_pos,ik)
!              fout%val(ii_pos,ik) = fin%val(ii_pos,ik)
!           enddo
!        enddo
!     enddo

   endif

end subroutine gf3_mirror_x

!------------------------------------------------------------------------------------------------------------

subroutine gf3_print(f, nam, ikarg, idebug, ndigit)
!--function: print information on grid-function f
   implicit none
!--subroutine arguments
   type(t_gridfnc3)             :: f
   character(len=*), intent(in) :: nam
   integer,          intent(in) :: ikarg
   integer,          intent(in) :: idebug
   integer, optional    :: ndigit       ! number of significant digits
!--local variables
   integer              :: nx, ny, ix, iy, ii, ik, ik0, ik1, my_ndigit, my_len
   character(len=20)    :: strng(3)

   if (present(ndigit)) then
      my_ndigit = ndigit
   else
      my_ndigit = 6
   endif
   my_ndigit = max(2, min(12, my_ndigit))
   my_len    = 8 + my_ndigit

   if (.not.associated(f%val)) then
      call write_log(' gf3_print: error: f%val-array not associated')
      return
   endif

   if (idebug.ge.1) then
      write(bufout,*) ' contents of grid-function ',trim(nam),':'
      call write_log(1, bufout)
   endif

   ! print information on the grid, depending on idebug-level

   call grid_print(f%grid, trim(nam)//'_grid', idebug)

   ! idebug>=3: full details: full list of values

   if (idebug.ge.3) then

      write(bufout,'(2a)') trim(nam), ' = ['
      call write_log(1, bufout)

      ! determine range for coordinate directions ik

      call gf3_ikrange(ikarg, ik0, ik1)

      ! write header according to ikarg

      if     (ikarg.eq.ikXDIR .or. ikarg.eq.jkXDIR) then
         call write_log('%       ii (  ix,  iy)          vx')
      elseif (ikarg.eq.ikYDIR .or. ikarg.eq.jkYDIR) then
         call write_log('%       ii (  ix,  iy)          vy')
      elseif (ikarg.eq.ikZDIR .or. ikarg.eq.jkZDIR) then
         call write_log('%       ii (  ix,  iy)          vn')
      elseif (ikarg.eq.ikTANG .or. ikarg.eq.jkTANG) then
         call write_log('%       ii (  ix,  iy)          vx         vy')
      elseif (ikarg.eq.ikALL  .or. ikarg.eq.jkALL ) then
         call write_log('%       ii (  ix,  iy)          vx         vy         vn')
      endif

      nx = f%grid%nx
      ny = f%grid%ny
      do iy = 1, ny
         do ix = 1, nx
            ii = ix + (iy-1)*nx
            do ik = ik0, ik1
               strng(ik) = fmt_gs(my_len, my_ndigit, f%val(ii,ik))
            enddo
            write(bufout, '(a,i6,a,i4,a,i4,a,3a)') ' ii=',ii,' (',ix,',',iy,'):',                       &
                                                                 (strng(ik)(1:my_len), ik=ik0, ik1)
            call write_log(1, bufout)
         enddo
      enddo

      call write_log('];')

   endif ! idebug>=4

end subroutine gf3_print

!------------------------------------------------------------------------------------------------------------

subroutine gf3_destroy(gf3)
!--purpose: cleanup, deallocate space for a grid function
   implicit none
!--subroutine arguments:
   type(t_gridfnc3)        :: gf3

   ! de-allocate values-array, nullify pointers

   if (associated(gf3%val)) deallocate(gf3%val)
   call gf3_nullify(gf3)

end subroutine gf3_destroy

!------------------------------------------------------------------------------------------------------------

subroutine gf3_test
!--purpose: test routines working on gf3s
   implicit none
!--subroutine arguments:
!--local variables:
   type(t_vec)       :: nref
   type(t_grid)      :: g
   type(t_gridfnc3)  :: gf, cnrm
   integer           :: irot, idebug, iy, my
   real(kind=8), parameter :: pi     = 4d0*atan(1d0)

   ! create uniform planar grid of 1 x 5 elements [0,-2,0] .. [0,2,0]

   my = 5
   call grid_create_uniform(g, nxarg=1, x0arg=1d0, dxarg=2d0, nyarg=my, y0arg=-2d0, dyarg=1d0, zarg=0d0)

   ! define reference orientation: 30deg w.r.t. global system

   nref = vec( 0d0, -sin(pi/6), cos(pi/6) )

   ! create input vector-field w : 15deg w.r.t. reference

   call gf3_new(gf, 'gf_test', g, nulify=.true., lzero=.true.)
   do iy = 1, my
      gf%vy(iy) = cos(pi/12)
      gf%vn(iy) = sin(pi/12)
   enddo

   ! create varying normal directions: 0, 30, 60, 90, 120deg w.r.t. global system

   call gf3_new(cnrm, 'cnrm', g, nulify=.true., lzero=.true.)
   do iy = 1, my
      cnrm%vy(iy) = -sin((iy-1)*pi/6d0)
      cnrm%vn(iy) =  cos((iy-1)*pi/6d0)
   enddo
   call gf3_print(cnrm, 'normals', ikALL, 5)

   ! rotate input v to curvilinear orientation

   irot = 1
   idebug = 1
   call gf3_cart2curv(gf, nref, cnrm, irot, idebug)

   ! print result

   call gf3_print(gf, 'rotated', ikALL, 4)

end subroutine gf3_test

!------------------------------------------------------------------------------------------------------------

end module m_gridfunc
