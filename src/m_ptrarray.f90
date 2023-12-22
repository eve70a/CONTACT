!------------------------------------------------------------------------------------------------------------
! m_ptrarray - functions for allocation/re-allocation dynamic arrays
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_ptrarray
   use m_globals
   use m_markers

   implicit none
   private

   ! Generic functions to allocate pointer-arrays:

   public  allocate_error
   public  reallocate_1d_bool
   public  reallocate_1d_int
   public  reallocate_2d_int
   public  reallocate_1d_real_ptr
   public  reallocate_1d_real_alloc
   public  reallocate_2d_real_ptr
   public  reallocate_2d_real_alloc
   public  reallocate_3d_cmplx
   public  reallocate_1d_char
   public  reallocate_1d_tvec
   public  reallocate_arr

   interface reallocate_arr
      module procedure reallocate_1d_bool
      module procedure reallocate_1d_int
      module procedure reallocate_2d_int
      module procedure reallocate_1d_real_ptr
      module procedure reallocate_1d_real_alloc
      module procedure reallocate_2d_real_ptr
      module procedure reallocate_2d_real_alloc
      module procedure reallocate_3d_cmplx
      module procedure reallocate_1d_char
      module procedure reallocate_1d_tvec
   end interface reallocate_arr

   public enlarge_1d_real
   public enlarge_2d_real
   public copy_1d_int
   public copy_1d_real

   interface copy_ptarray
      module procedure copy_1d_int
      module procedure copy_1d_real
   end interface copy_ptarray

   public destroy_1d_int
   public destroy_2d_int
   public destroy_1d_real
   public destroy_2d_real
   public destroy_1d_char
   public destroy_arr

   interface destroy_arr
      module procedure destroy_1d_int
      module procedure destroy_2d_int
      module procedure destroy_1d_real
      module procedure destroy_2d_real
      module procedure destroy_1d_char
   end interface destroy_arr

   public print_array_size_1d
   public print_array_size_2d

   public check_nan
   public check_nan_1d_real
   public check_nan_2d_real

   interface check_nan
      module procedure check_nan_1d_real
      module procedure check_nan_2d_real
   end interface check_nan

   public print_1d_int
   public print_2d_int
   public print_1d_real
   public print_2d_real

contains

!------------------------------------------------------------------------------------------------------------

subroutine allocate_error(namsub, lenarr, istat)
!--purpose: Write error message for memory allocation errors
   implicit none
!--subroutine parameters:
   integer          :: lenarr, istat
   character(len=*) :: namsub

   write(bufout,'(3a,i8,a,i4)') ' Error: memory allocation failed for ', namsub,',', lenarr,            &
                ' elements, istat=',istat
   call write_log(1, bufout)

end subroutine allocate_error

!------------------------------------------------------------------------------------------------------------

subroutine reallocate_1d_bool(larr, n1, keep)
!--purpose: Allocate array if not allocated before, else reallocate with the requested size,
!           optionally preserving the data if keep=.true.
   implicit none
!--subroutine parameters:
   integer                               :: n1
   logical,      dimension(:),   pointer :: larr
   logical,                     optional :: keep
!--local variables
   logical   :: keep_data
   integer   :: istat, ncur1, ncopy1
   logical, dimension(:), pointer  :: tmparr

   keep_data = .false.
   if (present(keep)) keep_data = keep

   if (.not.associated(larr)) then
      allocate(larr(n1), stat=istat)
      if (istat.ne.0) call allocate_error('1d_bool', n1, istat)
   endif

   if (.not.keep_data) then

      ! re-allocate without preserving data

      if (size(larr,1).ne.n1) then
         deallocate(larr)
         allocate(larr(n1), stat=istat)
         if (istat.ne.0) call allocate_error('1d_bool', n1, istat)
      endif

   else

      ! re-allocate with copying of data - enlarge or shrink existing array

      ncur1 = size(larr,1)
      if (ncur1.ne.n1) then
         allocate(tmparr(n1), stat=istat)
         if (istat.ne.0) call allocate_error('1d_bool', n1, istat)
         ncopy1 = min(n1, ncur1)
         tmparr(1:ncopy1) = larr(1:ncopy1)
         deallocate(larr)
         larr => tmparr
      endif
   endif

end subroutine reallocate_1d_bool

!------------------------------------------------------------------------------------------------------------

subroutine reallocate_1d_int(iarr, n1, keep)
!--purpose: Allocate array if not allocated before, else reallocate with the requested size,
!           optionally preserving the data if keep=.true.
   implicit none
!--subroutine parameters:
   integer                               :: n1
   integer,      dimension(:),   pointer :: iarr
   logical,                     optional :: keep
!--local variables
   logical   :: keep_data
   integer   :: istat, ncur1, ncopy1
   integer, dimension(:), pointer  :: tmparr

   keep_data = .false.
   if (present(keep)) keep_data = keep

   if (.not.associated(iarr)) then
      allocate(iarr(n1), stat=istat)
      if (istat.ne.0) call allocate_error('1d_int', n1, istat)
   endif

   if (.not.keep_data) then

      ! re-allocate without preserving data

      if (size(iarr,1).ne.n1) then
         deallocate(iarr)
         allocate(iarr(n1), stat=istat)
         if (istat.ne.0) call allocate_error('1d_int', n1, istat)
      endif

   else

      ! re-allocate with copying of data - enlarge or shrink existing array

      ncur1 = size(iarr,1)
      if (ncur1.ne.n1) then
         allocate(tmparr(n1), stat=istat)
         if (istat.ne.0) call allocate_error('1d_int', n1, istat)
         ncopy1 = min(n1, ncur1)
         tmparr(1:ncopy1) = iarr(1:ncopy1)
         deallocate(iarr)
         iarr => tmparr
      endif
   endif

end subroutine reallocate_1d_int

!------------------------------------------------------------------------------------------------------------

subroutine reallocate_2d_int(iarr, n1, n2, keep)
!--purpose: Allocate array if not allocated before, else reallocate with the requested size,
!           optionally preserving the data if keep=.true.
   implicit none
!--subroutine parameters:
   integer                               :: n1, n2
   integer,      dimension(:,:), pointer :: iarr
   logical,                     optional :: keep
!--local variables
   logical   :: keep_data
   integer   :: istat, ncur1, ncur2, ncopy1, ncopy2, icopy2
   integer, dimension(:,:), pointer  :: tmparr

   keep_data = .false.
   if (present(keep)) keep_data = keep

   if (.not.associated(iarr)) then
      allocate(iarr(n1,n2), stat=istat)
      if (istat.ne.0) call allocate_error('2d_int', n1*n2, istat)
   endif

   if (.not.keep_data) then

      ! re-allocate without preserving data

      if (size(iarr,1).ne.n1 .or. size(iarr,2).ne.n2) then
         deallocate(iarr)
         allocate(iarr(n1,n2), stat=istat)
         if (istat.ne.0) call allocate_error('2d_int', n1*n2, istat)
      endif

   else

      ! re-allocate with copying of data - enlarge or shrink existing array

      ncur1 = size(iarr,1)
      ncur2 = size(iarr,2)
      if (ncur1.ne.n1 .or. ncur2.ne.n2) then
         allocate(tmparr(n1,n2), stat=istat)
         if (istat.ne.0) call allocate_error('2d_int', n1*n2, istat)
         ncopy1 = min(n1, ncur1)
         ncopy2 = min(n2, ncur2)
         do icopy2 = 1, ncopy2
            tmparr(1:ncopy1,icopy2) = iarr(1:ncopy1,icopy2)
         enddo
         deallocate(iarr)
         iarr => tmparr
      endif
   endif

end subroutine reallocate_2d_int

!------------------------------------------------------------------------------------------------------------

subroutine reallocate_1d_real_ptr(rarr, n1, keep)
!--purpose: Allocate pointer-array if not allocated before, else reallocate with the requested size,
!           optionally preserving the data if keep=.true.
   implicit none
!--subroutine parameters:
   integer                               :: n1
   real(kind=8), dimension(:),   pointer :: rarr
   logical,                     optional :: keep
!--local variables
   logical   :: keep_data
   integer   :: istat, ncur1, ncopy1
   real(kind=8), dimension(:),  pointer  :: tmparr

   keep_data = .false.
   if (present(keep)) keep_data = keep

   if (.not.associated(rarr)) then
      allocate(rarr(n1), stat=istat)
      if (istat.ne.0) call allocate_error('1d_real', n1, istat)
   endif

   if (.not.keep_data) then

      ! re-allocate without preserving data

      if (size(rarr,1).ne.n1) then
         deallocate(rarr)
         allocate(rarr(n1), stat=istat)
         if (istat.ne.0) call allocate_error('1d_real', n1, istat)
      endif

   else

      ! re-allocate with copying of data - enlarge or shrink existing array

      ncur1 = size(rarr,1)
      if (ncur1.ne.n1) then
         allocate(tmparr(n1), stat=istat)
         if (istat.ne.0) call allocate_error('1d_real', n1, istat)
         ncopy1 = min(n1, ncur1)
         tmparr(1:ncopy1) = rarr(1:ncopy1)
         deallocate(rarr)
         rarr => tmparr
      endif
   endif

end subroutine reallocate_1d_real_ptr

!------------------------------------------------------------------------------------------------------------

subroutine reallocate_1d_real_alloc(rarr, n1)
!--purpose: Allocate allocatable array if not allocated before, else reallocate with the requested size,
   implicit none
!--subroutine parameters:
   integer                                 :: n1
   real(kind=8), dimension(:), allocatable :: rarr
!--local variables
   integer      :: istat

   if (.not.allocated(rarr)) then
      allocate(rarr(n1), stat=istat)
      if (istat.ne.0) call allocate_error('1d_real', n1, istat)
   endif

   ! re-allocate without preserving data

   if (size(rarr,1).ne.n1) then
      deallocate(rarr)
      allocate(rarr(n1), stat=istat)
      if (istat.ne.0) call allocate_error('1d_real', n1, istat)
   endif

end subroutine reallocate_1d_real_alloc

!------------------------------------------------------------------------------------------------------------

subroutine reallocate_2d_real_ptr(rarr, n1, n2, keep)
!--purpose: Allocate pointer-array if not allocated before, else reallocate with the requested size,
!           optionally preserving the data if keep=.true.
   implicit none
!--subroutine parameters:
   integer                                :: n1, n2
   real(kind=8), dimension(:,:), pointer  :: rarr
   logical,                      optional :: keep
!--local variables
   integer   :: icopy2, ncur1, ncur2, ncopy1, ncopy2, istat
   logical   :: keep_data
   real(kind=8), dimension(:,:), pointer  :: tmparr

   keep_data = .false.
   if (present(keep)) keep_data = keep

   if (.not.associated(rarr)) then
      allocate(rarr(n1,n2), stat=istat)
      if (istat.ne.0) call allocate_error('2d_real', n1*n2, istat)
   endif

   if (.not.keep_data) then

      ! re-allocate without preserving data

      if (size(rarr,1).ne.n1 .or. size(rarr,2).ne.n2) then
         deallocate(rarr)
         allocate(rarr(n1,n2), stat=istat)
         if (istat.ne.0) call allocate_error('2d_real', n1*n2, istat)
      endif

   else

      ! re-allocate with copying of data - enlarge or shrink existing array

      ncur1 = size(rarr,1)
      ncur2 = size(rarr,2)
      if (ncur1.ne.n1 .or. ncur2.ne.n2) then
         allocate(tmparr(n1,n2), stat=istat)
         if (istat.ne.0) call allocate_error('2d_real', n1*n2, istat)
         ncopy1 = min(n1, ncur1)
         ncopy2 = min(n2, ncur2)
         do icopy2 = 1, ncopy2
            tmparr(1:ncopy1,icopy2) = rarr(1:ncopy1,icopy2)
         enddo
         deallocate(rarr)
         rarr => tmparr
      endif
   endif

end subroutine reallocate_2d_real_ptr

!------------------------------------------------------------------------------------------------------------

subroutine reallocate_2d_real_alloc(rarr, n1, n2)
!--purpose: Allocate allocatable array if not allocated before, else reallocate with the requested size
   implicit none
!--subroutine parameters:
   integer                                   :: n1, n2
   real(kind=8), dimension(:,:), allocatable :: rarr
!--local variables
   integer   :: istat

   if (.not.allocated(rarr)) then
      allocate(rarr(n1,n2), stat=istat)
      if (istat.ne.0) call allocate_error('2d_real', n1*n2, istat)
   endif

   ! re-allocate without preserving data

   if (size(rarr,1).ne.n1 .or. size(rarr,2).ne.n2) then
      deallocate(rarr)
      allocate(rarr(n1,n2), stat=istat)
      if (istat.ne.0) call allocate_error('2d_real', n1*n2, istat)
   endif

end subroutine reallocate_2d_real_alloc

!------------------------------------------------------------------------------------------------------------

subroutine reallocate_3d_cmplx(carr, n1, n2, n3, keep)
!--purpose: Allocate array if not allocated before, else reallocate with the requested size,
!           optionally preserving the data if keep=.true.
   implicit none
!--subroutine parameters:
   integer                                    :: n1, n2, n3
   complex(kind=8), dimension(:,:,:), pointer :: carr
   logical,                          optional :: keep
!--local variables
   logical   :: keep_data
   integer   :: istat, ncur1, ncur2, ncur3, ncopy1, ncopy2, ncopy3, icopy2, icopy3
   complex(kind=8), dimension(:,:,:), pointer :: tmparr

   keep_data = .false.
   if (present(keep)) keep_data = keep

   if (.not.associated(carr)) then
      allocate(carr(n1,n2,n3), stat=istat)
      if (istat.ne.0) call allocate_error('3d_cmplx', n1*n2*n3, istat)
   endif

   if (.not.keep_data) then

      ! re-allocate without preserving data

      if (size(carr,1).ne.n1 .or. size(carr,2).ne.n2 .or. size(carr,3).ne.n3) then
         deallocate(carr)
         allocate(carr(n1,n2,n3), stat=istat)
         if (istat.ne.0) call allocate_error('3d_cmplx', n1*n2*n3, istat)
      endif

   else

      ! re-allocate with copying of data - enlarge or shrink existing array

      ncur1 = size(carr,1)
      ncur2 = size(carr,2)
      ncur3 = size(carr,3)
      if (ncur1.ne.n1 .or. ncur2.ne.n2 .or. ncur3.ne.n3) then
         allocate(tmparr(n1,n2,n3), stat=istat)
         if (istat.ne.0) call allocate_error('3d_cmplx', n1*n2*n3, istat)
         ncopy1 = min(n1, ncur1)
         ncopy2 = min(n2, ncur2)
         ncopy3 = min(n3, ncur3)
         do icopy3 = 1, ncopy3
            do icopy2 = 1, ncopy2
               tmparr(1:ncopy1,icopy2,icopy3) = carr(1:ncopy1,icopy2,icopy3)
            enddo
         enddo
         deallocate(carr)
         carr => tmparr
      endif
   endif

end subroutine reallocate_3d_cmplx

!------------------------------------------------------------------------------------------------------------

subroutine reallocate_1d_char(charr, n1, keep)
!--purpose: Allocate array if not allocated before, else reallocate with the requested size,
!           optionally preserving the data if keep=.true.
   implicit none
!--subroutine parameters:
   integer                                  :: n1
   character(len=*), dimension(:),  pointer :: charr
   logical,                        optional :: keep
!--local variables
   logical   :: keep_data
   integer   :: istat, ncur1, ncopy1
   character(len=:), dimension(:),  pointer :: tmparr

   keep_data = .false.
   if (present(keep)) keep_data = keep

   if (.not.associated(charr)) then
      allocate(charr(n1), stat=istat)
      if (istat.ne.0) call allocate_error('1d_char', n1, istat)
   endif

   if (.not.keep_data) then

      ! re-allocate without preserving data

      if (size(charr,1).ne.n1) then
         deallocate(charr)
         allocate(charr(n1), stat=istat)
         if (istat.ne.0) call allocate_error('1d_char', n1, istat)
      endif

   else

      ! re-allocate with copying of data - enlarge or shrink existing array

      ncur1 = size(charr,1)
      if (ncur1.ne.n1) then
         ! allocate tmparr with same len-parameter as charr (mold=)
         allocate(tmparr(n1), mold=charr, stat=istat)
         if (istat.ne.0) call allocate_error('1d_char', n1, istat)
         ncopy1 = min(n1, ncur1)
         tmparr(1:ncopy1) = charr(1:ncopy1)
         deallocate(charr)
         charr => tmparr
      endif
   endif

end subroutine reallocate_1d_char

!------------------------------------------------------------------------------------------------------------

subroutine reallocate_1d_tvec(varr, n1, keep)
!--purpose: Allocate array if not allocated before, else reallocate with the requested size,
!           optionally preserving the data if keep=.true.
   implicit none
!--subroutine parameters:
   integer                               :: n1
   type(t_vec),  dimension(:),   pointer :: varr
   logical,                     optional :: keep
!--local variables
   logical   :: keep_data
   integer   :: istat, ncur1, ncopy1
   type(t_vec), dimension(:),    pointer :: tmparr

   keep_data = .false.
   if (present(keep)) keep_data = keep

   if (.not.associated(varr)) then
      allocate(varr(n1), stat=istat)
      if (istat.ne.0) call allocate_error('1d_tvec', 3*n1, istat)
   endif

   if (.not.keep_data) then

      ! re-allocate without preserving data

      if (size(varr,1).ne.n1) then
         deallocate(varr)
         allocate(varr(n1), stat=istat)
         if (istat.ne.0) call allocate_error('1d_tvec', 3*n1, istat)
      endif

   else

      ! re-allocate with copying of data - enlarge or shrink existing array

      ncur1 = size(varr,1)
      if (ncur1.ne.n1) then
         allocate(tmparr(n1), stat=istat)
         if (istat.ne.0) call allocate_error('1d_tvec', n1, istat)
         ncopy1 = min(n1, ncur1)
         tmparr(1:ncopy1) = varr(1:ncopy1)
         deallocate(varr)
         varr => tmparr
      endif
   endif

end subroutine reallocate_1d_tvec

!------------------------------------------------------------------------------------------------------------

subroutine enlarge_1d_real(rarr, nneed, ngrow)
!--purpose: Check if array provides for n1 spaces, if not, extend to size nneed+ngrow, preserving data
   implicit none
!--subroutine parameters:
   integer                                :: nneed, ngrow
   real(kind=8), dimension(:),   pointer  :: rarr
!--local variables:
   integer                                :: ncur, istat
   real(kind=8), dimension(:),   pointer  :: tmparr

   if (.not.associated(rarr)) then
      allocate(rarr(nneed+ngrow), stat=istat)
      if (istat.ne.0) call allocate_error('+1d_real', nneed+ngrow, istat)
   endif

   ncur = size(rarr,1)
   if (ncur.lt.nneed) then
      allocate(tmparr(nneed+ngrow), stat=istat)
      if (istat.ne.0) call allocate_error('+1d_real', nneed+ngrow, istat)
      tmparr(1:ncur) = rarr(1:ncur)
      deallocate(rarr)
      rarr => tmparr
   endif

end subroutine enlarge_1d_real

!------------------------------------------------------------------------------------------------------------

subroutine enlarge_2d_real(rarr, nneed1, nneed2, ngrow1, ngrow2)
!--purpose: Check if array provides for (n1,n2) spaces, if not, extend to size nneed+ngrow, preserving data
   implicit none
!--subroutine parameters:
   integer                                :: nneed1, nneed2, ngrow1, ngrow2
   real(kind=8), dimension(:,:), pointer  :: rarr
!--local variables:
   integer                                :: ncur1, ncur2, icur2, istat
   real(kind=8), dimension(:,:), pointer  :: tmparr

   if (.not.associated(rarr)) then
      allocate(rarr(nneed1+ngrow1,nneed2+ngrow2), stat=istat)
      if (istat.ne.0) call allocate_error('+2d_real', (nneed1+ngrow1)*(nneed2+ngrow2), istat)
   endif

   ncur1 = size(rarr,1)
   ncur2 = size(rarr,2)
   if (ncur1.lt.nneed1 .or. ncur2.lt.nneed2) then
      allocate(tmparr(nneed1+ngrow1,nneed2+ngrow2), stat=istat)
      if (istat.ne.0) call allocate_error('+2d_real', (nneed1+ngrow1)*(nneed2+ngrow2), istat)
      do icur2 = 1, ncur2
         tmparr(1:ncur1,icur2) = rarr(1:ncur1,icur2)
      enddo
      deallocate(rarr)
      rarr => tmparr
   endif

end subroutine enlarge_2d_real

!------------------------------------------------------------------------------------------------------------

subroutine copy_1d_int(iarr_in, iarr_out)
!--purpose: Copy allocated array iarr_in to iarr_out, (re)allocating memory if necessary
   implicit none
!--subroutine parameters:
   integer, dimension(:),   pointer  :: iarr_in, iarr_out
!--local variables
   integer    :: n1

   if (.not.associated(iarr_in)) then
      if (associated(iarr_out)) deallocate(iarr_out)
      iarr_out => NULL()
   else
      n1 = size(iarr_in)
      call reallocate_arr(iarr_out, n1)
      iarr_out(1:n1) = iarr_in(1:n1)
   endif

end subroutine copy_1d_int

!------------------------------------------------------------------------------------------------------------

subroutine copy_1d_real(rarr_in, rarr_out)
!--purpose: Copy allocated array rarr_in to rarr_out, (re)allocating memory if necessary
   implicit none
!--subroutine parameters:
   real(kind=8), dimension(:),   pointer  :: rarr_in, rarr_out
!--local variables
   integer    :: n1

   if (.not.associated(rarr_in)) then
      if (associated(rarr_out)) deallocate(rarr_out)
      rarr_out => NULL()
   else
      n1 = size(rarr_in)
      call reallocate_arr(rarr_out, n1)
      rarr_out(1:n1) = rarr_in(1:n1)
   endif

end subroutine copy_1d_real

!------------------------------------------------------------------------------------------------------------

subroutine destroy_1d_int(iarr)
   implicit none
   integer, dimension(:), pointer :: iarr

   if (associated(iarr)) then
      deallocate(iarr)
      nullify(iarr)
   endif
end subroutine destroy_1d_int

!------------------------------------------------------------------------------------------------------------

subroutine destroy_2d_int(i2arr)
   implicit none
   integer, dimension(:,:), pointer :: i2arr

   if (associated(i2arr)) then
      deallocate(i2arr)
      nullify(i2arr)
   endif
end subroutine destroy_2d_int

!------------------------------------------------------------------------------------------------------------

subroutine destroy_1d_real(rarr)
   implicit none
   real(kind=8), dimension(:), pointer :: rarr

   if (associated(rarr)) then
      deallocate(rarr)
      nullify(rarr)
   endif
end subroutine destroy_1d_real

!------------------------------------------------------------------------------------------------------------

subroutine destroy_2d_real(r2arr)
   implicit none
   real(kind=8), dimension(:,:), pointer :: r2arr

   if (associated(r2arr)) then
      deallocate(r2arr)
      nullify(r2arr)
   endif
end subroutine destroy_2d_real

!------------------------------------------------------------------------------------------------------------

subroutine destroy_1d_char(charr)
   implicit none
   character(len=*), dimension(:),  pointer :: charr

   if (associated(charr)) then
      deallocate(charr)
      nullify(charr)
   endif
end subroutine destroy_1d_char

!------------------------------------------------------------------------------------------------------------

subroutine print_array_size_1d(arr, nam)
!--purpose: helper for grid_print, print array allocation and size
   implicit none
!--subroutine parameters:
   real(kind=8),    dimension(:), pointer  :: arr
   character(len=*)                        :: nam

   if (.not.associated(arr)) then
      call write_log('  array ' // trim(nam) // ' has not yet been allocated')
   else
      write(bufout,'(3a,i6,a)') '  array ', trim(nam), ' has', size(arr,1),' elements'
      call write_log(1, bufout)
   endif

end subroutine print_array_size_1d

!------------------------------------------------------------------------------------------------------------

subroutine print_array_size_2d(arr, nam)
!--purpose: helper for grid_print, print array allocation and size
   implicit none
!--subroutine parameters:
   real(kind=8),    dimension(:,:), pointer  :: arr
   character(len=*)                          :: nam

   if (.not.associated(arr)) then
      call write_log('  array ' // trim(nam) // ' has not yet been allocated')
   else
      write(bufout,'(2a,2(a,i6),a)') '  array ', trim(nam), ' has', size(arr,1),' x', size(arr,2),' elements'
      call write_log(1, bufout)
   endif

end subroutine print_array_size_2d

!------------------------------------------------------------------------------------------------------------

subroutine check_nan_1d_real(arr, nam, idebug)
!--purpose: helper for debugging, print number of NaN-values in array
   implicit none
!--subroutine parameters:
   integer                         :: idebug
   real(kind=8),    dimension(:)   :: arr
   character(len=*)                :: nam
!  local variables:
   integer      :: i, nnan

   nnan = 0
   do i = 1, size(arr,1)
      if (isnan(arr(i))) nnan = nnan + 1
   enddo

   if (nnan.gt.0 .and. idebug.ge.-1) then
      write(bufout,'(a,i6,3a)') ' Error: there are',nnan,' NaN-values in array "',trim(nam),'"'
      call write_log(1, bufout)
   elseif (idebug.ge.2) then
      call write_log(' no NaN-values in array "' // trim(nam) // '"')
   endif

end subroutine check_nan_1d_real

!------------------------------------------------------------------------------------------------------------

subroutine check_nan_2d_real(arr, nam, idebug)
!--purpose: helper for debugging, print number of NaN-values in array
   implicit none
!--subroutine parameters:
   integer                         :: idebug
   real(kind=8),    dimension(:,:) :: arr
   character(len=*)                :: nam
!  local variables:
   integer      :: i, j, nnan

   nnan = 0
   do j = 1, size(arr,2)
      do i = 1, size(arr,1)
         if (isnan(arr(i,j))) nnan = nnan + 1
      enddo
   enddo

   if (nnan.gt.0 .and. idebug.ge.-1) then
      write(bufout,'(a,i6,3a)') ' Error: there are',nnan,' NaN-values in array "',trim(nam),'"'
      call write_log(1, bufout)
   elseif (idebug.ge.2) then
      call write_log(' no NaN-values in array "' // trim(nam) // '"')
   endif

end subroutine check_nan_2d_real

!------------------------------------------------------------------------------------------------------------

subroutine print_1d_int(ntot, iarr, nam, n_per_line, fmtval)
!--purpose: helper for pretty printing array of numerical values
   implicit none
!--subroutine parameters:
   integer,          intent(in)  :: ntot, n_per_line
   integer,          intent(in)  :: iarr(ntot)
   character(len=*), intent(in)  :: nam, fmtval
!  local variables:
   character(len=30)  :: fmtstr
   integer            :: i, iofs, iline, nline

   ! determine number of lines

   nline = int( (ntot - 1) / n_per_line ) + 1

   call write_log(trim(nam) // ' = [')

   do iline = 1, nline
      iofs = (iline - 1) * n_per_line
      if (iline.lt.nline) then
         write(fmtstr,'(a,i2,2a)') '(',n_per_line, fmtval, ',a)'
         write(bufout, fmtstr) (iarr(i), i=iofs+1, iofs+n_per_line), ' ...'
      else
         write(fmtstr,'(a,i2,2a)') '(', ntot-iofs, fmtval, ')'
         write(bufout, fmtstr) (iarr(i), i=iofs+1, ntot)
      endif
      call write_log(1, bufout)
   enddo

   call write_log(' ];')

end subroutine print_1d_int

!------------------------------------------------------------------------------------------------------------

subroutine print_2d_int(nrow, ncol, iarr, nam, ncol_per_line, fmtval)
!--purpose: helper for pretty printing array of numerical values
   implicit none
!--subroutine parameters:
   integer,          intent(in)  :: nrow, ncol, ncol_per_line
   integer,          intent(in)  :: iarr(nrow, ncol)
   character(len=*), intent(in)  :: nam, fmtval
!  local variables:
   character(len=30)  :: fmtstr
   integer            :: irow, icol, iofs, iline, nline

   ! determine number of lines used per row

   nline = int( (ncol - 1) / ncol_per_line ) + 1

   write(bufout,'(2a,2(i4,a))') trim(nam), ' = [  % size (', nrow,',',ncol,')'
   call write_log(1, bufout)

   do irow = 1, nrow
      do iline = 1, nline
         iofs = (iline - 1) * ncol_per_line
         if (iline.lt.nline) then
            write(fmtstr,'(a,i2,2a)') '(',ncol_per_line, fmtval, ',a)'
            write(bufout, fmtstr) (iarr(irow,icol), icol=iofs+1, iofs+ncol_per_line), ' ...'
         else
            write(fmtstr,'(a,i2,2a)') '(', ncol-iofs, fmtval, ',a)'
            if (irow.lt.nrow) then
               write(bufout, fmtstr) (iarr(irow,icol), icol=iofs+1, ncol), ';'
            else
               write(bufout, fmtstr) (iarr(irow,icol), icol=iofs+1, ncol), ' '
            endif
         endif
         call write_log(1, bufout)
      enddo
   enddo

   call write_log(' ];')

end subroutine print_2d_int

!------------------------------------------------------------------------------------------------------------

subroutine print_1d_real(ntot, rarr, nam, n_per_line, fmtval)
!--purpose: helper for pretty printing array of numerical values
   implicit none
!--subroutine parameters:
   integer,          intent(in)  :: ntot, n_per_line
   real(kind=8),     intent(in)  :: rarr(ntot)
   character(len=*), intent(in)  :: nam, fmtval
!  local variables:
   character(len=30)  :: fmtstr
   integer            :: i, iofs, iline, nline

   ! determine number of lines

   nline = int( (ntot - 1) / n_per_line ) + 1

   call write_log(trim(nam) // ' = [')

   do iline = 1, nline
      iofs = (iline - 1) * n_per_line
      if (iline.lt.nline) then
         write(fmtstr,'(a,i2,2a)') '(',n_per_line, fmtval, ',a)'
         write(bufout, fmtstr) (rarr(i), i=iofs+1, iofs+n_per_line), ' ...'
      else
         write(fmtstr,'(a,i2,2a)') '(', ntot-iofs, fmtval, ')'
         write(bufout, fmtstr) (rarr(i), i=iofs+1, ntot)
      endif
      call write_log(1, bufout)
   enddo

   call write_log(' ];')

end subroutine print_1d_real

!------------------------------------------------------------------------------------------------------------

subroutine print_2d_real(nrow, ncol, rarr, nam, ncol_per_line, fmtval)
!--purpose: helper for pretty printing array of numerical values
   implicit none
!--subroutine parameters:
   integer,          intent(in)  :: nrow, ncol, ncol_per_line
   real(kind=8),     intent(in)  :: rarr(nrow, ncol)
   character(len=*), intent(in)  :: nam, fmtval
!  local variables:
   character(len=30)  :: fmtstr
   integer            :: irow, icol, iofs, iline, nline

   ! determine number of lines used per row

   nline = int( (ncol - 1) / ncol_per_line ) + 1

   write(bufout,'(2a,2(i4,a))') trim(nam), ' = [  % size (', nrow,',',ncol,')'
   call write_log(1, bufout)

   do irow = 1, nrow
      do iline = 1, nline
         iofs = (iline - 1) * ncol_per_line
         if (iline.lt.nline) then
            write(fmtstr,'(a,i2,2a)') '(',ncol_per_line, fmtval, ',a)'
            write(bufout, fmtstr) (rarr(irow,icol), icol=iofs+1, iofs+ncol_per_line), ' ...'
         else
            write(fmtstr,'(a,i2,2a)') '(', ncol-iofs, fmtval, ',a)'
            if (irow.lt.nrow) then
               write(bufout, fmtstr) (rarr(irow,icol), icol=iofs+1, ncol), ';'
            else
               write(bufout, fmtstr) (rarr(irow,icol), icol=iofs+1, ncol), ' '
            endif
         endif
         call write_log(1, bufout)
      enddo
   enddo

   call write_log(' ];')

end subroutine print_2d_real

!------------------------------------------------------------------------------------------------------------

end module m_ptrarray

!------------------------------------------------------------------------------------------------------------
