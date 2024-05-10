!------------------------------------------------------------------------------------------------------------
! m_gridfunc_sub - implementation of subroutines for m_gridfunc
!
! Copyright 1993-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
submodule (m_gridfunc) m_gridfunc_sub
   use m_globals
   use m_markers
   use m_ptrarray
   use m_grids
   implicit none

!------------------------------------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------------------------------------

module subroutine eldiv_resize(e1, g_new)
!--purpose: resize element division to the new (matching) grid provided
   implicit none
!--subroutine arguments:
   type(t_eldiv)           :: e1
   type(t_grid),  target   :: g_new
!--local variables:
   integer, parameter :: defval = Exter, idebug = 0
   logical            :: is_ok, is_equal
   integer            :: kofs_x, kofs_y, ix, ix0, ix1, iy, iy0, iy1, ii1, ii2
   integer, dimension(:), pointer :: el_new

   if (.not.associated(e1%grid)) then
      call write_log(' eldiv_resize: Internal error: e1%grid is not associated, aborting.')
      call abort_run()
   endif

   associate(g_old => e1%grid)

   ! check input grids and determine overlap region

   call grid_get_overlap(g_old, g_new, kofs_x, kofs_y, ix0, ix1, iy0, iy1, is_ok, is_equal)

   if (.not.is_ok) then
      call write_log(' eldiv_resize: Internal error: grids do not match.')
      call abort_run()
   endif

   ! check input element divisions

   if (.not.associated(e1%el)) &
      call write_log(' eldiv_resize: Internal error: e1%el-array not associated')
   if (.not.associated(e1%row1st)) &
      call write_log(' eldiv_resize: Internal error: e1%row1st-array not associated')
   if (.not.associated(e1%rowlst)) &
      call write_log(' eldiv_resize: Internal error: e1%rowlst-array not associated')

   ! reallocate/copy only if grid is changed

   if (.not.is_equal) then

      ! copy element division to new array

      allocate(el_new(g_new%ntot))

      do iy = 1, g_new%ny
         do ix = 1, g_new%nx
            ii1 = ix+kofs_x + (iy+kofs_y-1) * g_old%nx
            ii2 = ix        + (iy       -1) * g_new%nx
            if (ix.ge.ix0 .and. ix.le.ix1 .and. iy.ge.iy0 .and. iy.le.iy1) then

               el_new(ii2) = e1%el(ii1)

            else

               el_new(ii2) = defval        ! defval == Exter

            endif
         enddo ! ix
      enddo ! iy

      ! destroy original array, replace with new array

      deallocate(e1%el)
      e1%el => el_new

   endif

   ! store pointer to the grid on which this eldiv lives

   e1%grid => g_new

   ! re-allocate administration arrays, ensure consistency

   call reallocate_arr(e1%row1st, g_new%ny)
   call reallocate_arr(e1%rowlst, g_new%ny)
   call areas(e1)

   end associate
end subroutine eldiv_resize

!------------------------------------------------------------------------------------------------------------

module subroutine gf3_resize(gf, g_new, defval)
!--purpose: resize grid function to the new (matching) grid provided
   implicit none
!--subroutine arguments:
   type(t_gridfnc3)          :: gf
   type(t_grid),    target   :: g_new
   real(kind=8),    optional :: defval
!--local variables:
   integer, parameter :: idebug = 0
   logical            :: is_ok, is_equal
   integer            :: kofs_x, kofs_y, ix, ix0, ix1, iy, iy0, iy1, ii1, ii2
   real(kind=8), dimension(:,:), pointer :: val_new

   associate(g_old => gf%grid)

   ! check input grids and determine overlap region

   call grid_get_overlap(g_old, g_new, kofs_x, kofs_y, ix0, ix1, iy0, iy1, is_ok, is_equal)

   if (.not.is_ok) then
      call write_log(' gf3_resize: Internal error: cannot copy between old/new grids provided')
      call abort_run()
   endif

   ! check input grid function

   if (.not.associated(gf%val)) then
      call write_log(' gf3_resize: Internal error: values not available in input grid-func ' //   &
                       trim(gf%name))
      call abort_run()
   endif

   ! reallocate/copy only if grid is changed

   if (.not.is_equal) then

      ! allocate new values-array

      allocate(val_new(g_new%ntot,3))

      ! copy grid function data

      do iy = 1, g_new%ny
         do ix = 1, g_new%nx
            ii1 = ix+kofs_x + (iy+kofs_y-1) * g_old%nx
            ii2 = ix        + (iy       -1) * g_new%nx
            if (ix.ge.ix0 .and. ix.le.ix1 .and. iy.ge.iy0 .and. iy.le.iy1) then
               val_new(ii2,1:3) = gf%val(ii1,1:3)
            else
               if (present(defval)) val_new(ii2,1:3) = defval
            endif
         enddo ! ix
      enddo ! iy

      ! destroy original array, replace with new array

      deallocate(gf%val)
      gf%val => val_new
      gf%vx  => gf%val(:,ikXDIR)
      gf%vy  => gf%val(:,ikYDIR)
      gf%vn  => gf%val(:,ikZDIR)
      gf%vt  => gf%val(:,ikXDIR)

   endif

   ! store pointer to the grid on which this gridfunc lives

   gf%grid => g_new

   end associate

end subroutine gf3_resize

!------------------------------------------------------------------------------------------------------------

module subroutine wrigs(igs, is_roll, chi, ltight_arg)
   implicit none
!--subroutine arguments:
   type(t_eldiv), intent(in)           :: igs
   logical,       intent(in)           :: is_roll
   real(kind=8),  intent(in)           :: chi
   logical,       intent(in), optional :: ltight_arg
!--local variables:
   integer, parameter :: ScrRow = 800
   integer, parameter :: ScrCol = 800
   integer, parameter :: ArrLen = max(ScrRow, ScrCol)
        ! note that ArrLen is also hard-coded in the format strings below
   integer, parameter :: HlfCol =  36
   logical          :: ltight
   integer          :: nx, ny, nxloc, nyloc, ix, ix0, ix1, iy, iy0, iy1, i, ii
   character        :: str(ArrLen), nums(ArrLen), tens(ArrLen), hundr(ArrLen)
   real(kind=8),     parameter :: pi = 4d0*atan(1d0)
   character(len=1), parameter :: aset(Exter:Plast) = (/ '.', '*', 'S', '|' /)
   character(len=1), parameter :: dig(10) = (/ '1', '2', '3', '4', '5', '6', '7', '8', '9', '0' /)
   character(len=1000) :: bufloc(1)

   ltight = .true.
   if (present(ltight_arg)) ltight = ltight_arg

   ! initialize strings nums, tens and str
   ! TODO: only on first call, save variables

   do ix = 1, ArrLen
      i = mod(ix-1, 10) + 1
      nums(ix) = dig(i)
      if (mod(ix,10).eq.0) then
         i = mod(ix/10-1, 10) + 1
         tens(ix) = dig(i)
      else
         tens(ix) = ' '
      endif
      if (mod(ix,100).eq.0) then
         i = mod(ix/100-1, 10) + 1
         hundr(ix) = dig(i)
      else
         hundr(ix) = ' '
      endif
   enddo
   do ix = 1, ArrLen
      str(ix) = ' '
   enddo

   ! determine range to be displayed

   nx = igs%grid%nx
   ny = igs%grid%ny

   if (ltight) then
      call areas(igs)
      if (igs%ixmin.gt.igs%ixmax) then
         ! empty contact area
         ix0 = max( 1, nx/2-5)
         ix1 = min(nx, nx/2+5)
      elseif (igs%ixmin.le.5 .and. igs%ixmax.ge.nx-4) then
         ! contact close to edges of pot.con
         ix0 = 1
         ix1 = nx
      else
         ! remove white-space
         ix0 = max( 1, igs%ixmin-1)
         ix1 = min(nx, igs%ixmax+1)
      endif
      if (igs%iymin.gt.igs%iymax) then
         iy0 = max( 1, ny/2-5)
         iy1 = min(ny, ny/2+5)
      elseif (igs%iymin.le.5 .and. igs%iymax.ge.ny-4) then
         iy0 = 1
         iy1 = ny
      else
         iy0 = max( 1, igs%iymin-1)
         iy1 = min(ny, igs%iymax+1)
      endif
   else
      ix0 = 1
      ix1 = nx
      iy0 = 1
      iy1 = ny
   endif
   nxloc = ix1 - ix0 + 1
   nyloc = iy1 - iy0 + 1

   if (nxloc.le.HlfCol .and. nyloc.le.ScrRow) then

      ! If nxloc is relatively small, print the picture in the normal orientation and add spaces between
      ! the consecutive points

      do iy = iy1, iy0, -1
         do ix = ix0, ix1
            ii = ix + nx * (iy-1)
            str(ix) = aset(max(Exter,min(Plast,igs%el(ii))))
            if (igs%el(ii).le.Exter .and. (ix.eq.ix0 .or. ix.eq.ix1) &
                                    .and. (iy.eq.iy0 .or. iy.eq.iy1)) str(ix) = 'o'
         enddo
         write (bufloc, 2200) iy, (str(ix), ix=ix0,ix1)
         call write_log(1, bufloc)
 2200    format (i6, 2x, 900(' ',a,:))
      enddo

      if (nxloc.ge.200) then
         write (bufloc, 2300) (hundr(ix), ix=ix0,ix1)
         call write_log(1, bufloc)
      endif
      if (nxloc.ge.10) then
         write (bufloc, 2300) (tens(ix), ix=ix0,ix1)
         call write_log(1, bufloc)
      endif
      write (bufloc, 2300) (nums(ix), ix=ix0,ix1)
      call write_log(1, bufloc)
 2300 format (8x, 900(' ',a,:))

      if (.not.is_roll) then
         call write_log('     X  -->')
      elseif (chi.gt.-998d0) then
         write (bufloc, 2400) chi * 180d0/pi
         call write_log(1, bufloc)
 2400    format ('     X  -->        Chi=', f6.2, ' degrees')
      endif

   elseif (nxloc.le.ScrCol .and. nyloc.le.ScrRow) then

      ! Elseif nxloc is small enough w.r.t. the max number of columns, print the picture in normal
      ! orientation without the spaces

      do iy = iy1, iy0, -1
         do ix = ix0, ix1
            ii = ix + nx * (iy-1)
            str(ix) = aset(max(Exter,min(Plast,igs%el(ii))))
            if (igs%el(ii).le.Exter .and. (ix.eq.ix0 .or. ix.eq.ix1) &
                                    .and. (iy.eq.iy0 .or. iy.eq.iy1)) str(ix) = 'o'
         enddo

         ! write row number and row

         write (bufloc, 5200) iy, (str(ix), ix=ix0,ix1)
         call write_log(1, bufloc)
 5200    format (i6, 2x, 900a)
      enddo

      ! Display digits under figure

      if (nxloc.ge.200) then
         write (bufloc, 5300) (hundr(ix), ix=ix0,ix1)
         call write_log(1, bufloc)
      endif
      if (nxloc.ge.10) then
         write (bufloc, 5300) (tens(ix), ix=ix0,ix1)
         call write_log(1, bufloc)
      endif
      write (bufloc, 5300) (nums(ix), ix=ix0,ix1)
      call write_log(1, bufloc)
 5300 format (8x, 900a)

      if (.not.is_roll) then
         call write_log('     X  -->')
      elseif (chi.gt.-998d0) then
         write (bufloc, 2400) chi * 180d0/pi
         call write_log(1, bufloc)
      endif

   elseif (nyloc.gt.ScrRow .and. nyloc.le.ScrCol .and. nxloc.le.ScrRow) then

   ! If nyloc is the limiting factor and ScrCol > ScrRow, then the picture
   !    may be rotated over 90 deg counter-clockwise.
   !    Here no spaces are added in between consecutive points.

      do ix = ix1, ix0, -1
         do iy = iy1, iy0, -1
            ii = ix + nx * (iy-1)
            str(ny-iy+1) = aset(max(Exter,min(Plast,igs%el(ii))))
            if (igs%el(ii).le.Exter .and. (ix.eq.ix0 .or. ix.eq.ix1) &
                                    .and. (iy.eq.iy0 .or. iy.eq.iy1)) str(ny-iy+1) = 'o'
         enddo
         str(iy1+1) = ' '
         str(iy1+2) = ' '
         str(iy1+3) = nums(ix)
         write (bufloc, 7200) (str(iy), iy=iy0,iy1+3)
         call write_log(1, bufloc)
 7200    format (8x, 903a)
      enddo

      if (nyloc.ge.200) then
         write (bufloc, 5300) (hundr(iy), iy=iy1,iy0,-1)
         call write_log(1, bufloc)
      endif
      if (nyloc.ge.10) then
         write (bufloc, 5300) (tens(iy), iy=iy1,iy0,-1)
         call write_log(1, bufloc)
      endif
      write (bufloc, 5300) (nums(iy), iy=iy1,iy0,-1)
      call write_log(1, bufloc)
      if (.not.is_roll) then
         call write_log('     <--  Y')
      elseif (chi.gt.-998d0) then
         write (bufloc, 7300) chi * 180d0/pi
         call write_log(1, bufloc)
 7300    format ('     <--  Y        Chi=', f6.2, ' degrees')
      endif
   endif

end subroutine wrigs

!------------------------------------------------------------------------------------------------------------

end submodule m_gridfunc_sub
