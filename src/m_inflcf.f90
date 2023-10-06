!------------------------------------------------------------------------------------------------------------
! m_inflcf - data-structures for influence coefficients
!
! Copyright 1993-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_inflcf

use m_globals
use m_readline
use m_markers
use m_grids

implicit none
private

public  t_inflcf
public  t_influe

public  inflcf_new
public  inflcf_copy
public  inflcf_print
public  inflcf_destroy
public  influe_load
private influe_loadtyp0
private influe_loadtyp1
public  influe_blanco
private influe_blanco0
private influe_blanco1
private influe_cv_csv
private influe_yv_ysv

!------------------------------------------------------------------------------------------------------------
! data-type for a single influence coefficients matrix:
!------------------------------------------------------------------------------------------------------------

   type :: t_inflcf
      real(kind=8),    dimension(:,:,:,:),   allocatable :: cf
      real(kind=8),    dimension(:,:,:,:,:), allocatable :: cy
      real(kind=8)                                       :: xoffs
      real(kind=8)                                       :: yoffs
      type(t_grid),                          pointer     :: grid   => NULL()
      integer                                            :: itypcf
      logical                                            :: nt_cpl
      real(kind=8)                                       :: ga, ga_inv
      logical                                            :: use_3bl, use_flxz
      real(kind=8)                                       :: flx_3bl, flx_z
      logical,         dimension(:,:)                    :: fft_ok(3,3)
      integer                                            :: fft_mx
      integer                                            :: fft_my
      complex(kind=8), dimension(:,:,:),     pointer     :: fft_cf => NULL()

      ! itypcf   type of influence matrix.
      !            0 == standard method, using cf
      !            1 == varying with lateral (y) position, using cy
      ! nt_cpl   coupling between normal and tangential problems
      !            false == no coupling, e.g. Ak==0.
      !            true  == coupling, Ak<>0 or num.infl.cf.
      ! ga       combined modulus of rigidity - scaling factor not included in coefficients themselves
      ! ga_inv   1/ga, scaling factor not included in coefficients themselves
      ! use_3bl  flag indicating whether a 3rd body/interfacial layer is included (h3>0) or not
      ! flx_3bl  flexiblity of the 3rd body/interfacial layer, if included
      ! use_flxz flag indicating whether a compressible sheet is included (flxz>0) or not
      ! flx_z    flexiblity of the compressible sheet, if included
      ! cf       influence matrix, describing the displacements in all elements of an "observation grid"
      !          due to a unit traction in element (0,0) of a "tractions grid". The two grids are assumed
      !          to be identical except for an optional displacement (xoffs,yoffs).
      !          Element cf(ix-jx, iy-jy, ik, jk) is the displacement in ik-direction in element (ix,iy)
      !          of the observation grid due to a unit traction in jk-direction in element (jx,jy) of the
      !          tractions grid.
      !          If the two grids are the same, then element cf(ix,iy,ik,jk) is the displacement in (ix,iy)
      !          due to a unit load in the  element (0,0).
      !          Cf has size (-mx+1:mx, -my+1:my, 3, 3). Grid column (mx,:) and row (:,my) are referenced
      !          in the program but do not affect the calculations.
      !          The coefficients are scaled by GA (historical reasons(?))
      ! cy       same as cf, with separate influence matrices for all grid lines jy (3rd array dimension)
      ! xoffs    x-offset of the tractions grid with respect to the observation grid. Particularly used for
      !          computing the effect of previous tractions at the current locations of surface particles.
      !          For rolling with chi=0, the coordinate system moves forward over distance Dq per step, such
      !          that tractions pv(jx,jy) should be located at (x(jx)-dq,y(jy)) instead of (x(jx),y(jy)).
      ! yoffs    y-offset between tractions and displacements grids, see xoffs.
      ! grid     pointer to grid-data on which infl.coefficients are defined, particularly mx,my, dx,dy
                 
      ! fft_ok   flag array, describing for which coordinate direction pairs (ik,jk) the fft transform of
      !          the influence coefficients is available in fft_cf.
      ! fft_mx   grid size in x-direction used in the Fourier transform
      ! fft_my   grid size in y-direction used in the Fourier transform
      ! fft_cf   Fourier transform of influence coefficients cf. Stored in a 3D array of size (length,3,3),
      !          with length (fft_mx+1)*2*fft_my for the CCE storage format.

   end type t_inflcf

!------------------------------------------------------------------------------------------------------------
! data for all influence coefficient matrices for one problem:
!------------------------------------------------------------------------------------------------------------

   type :: t_influe
      type(t_inflcf) :: cs
      type(t_inflcf) :: cv
      type(t_inflcf) :: csv
      type(t_inflcf) :: ms

      ! cs      influence coefficients for computing displacements at current grid due to current tractions
      ! cv      influence coefficients for computing displacements at current grid due to previous tractions
      ! csv     influence coefficients for the stationary problem, cs - cv
      ! ms      influence coefficients for preconditioner to normal problem

   end type t_influe

!------------------------------------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------------------------------------

   subroutine inflcf_new(inflcf, itypcf, grid)
!--purpose: (re-)allocate space for an influence coefficients matrix on the given grid
      implicit none
!--subroutine arguments:
      integer         , intent(in)    :: itypcf
      type(t_grid),     target        :: grid
      type(t_inflcf)  , intent(inout) :: inflcf
!--local variables:
      integer    :: nmatrix
!--pointers to member variables:
      integer,                    pointer :: mx, my

      mx    => grid%nx
      my    => grid%ny

      ! Store pointer to grid in influence matrix

      inflcf%grid  => grid

      if (itypcf.eq.0) then

         ! Use standard method, same influence matrix for all elements

         inflcf%itypcf = 0

         ! Allocate component-array cf if not done so before, reallocate when size must be changed
         ! Note: the arrays can be of size (-mx+1:mx-1,-my+1:my-1). One additional grid row and column
         !       are added at the left/bottom to obtain more favourable sizes for Fast Fourier transform.
         !       Elements (-mx,:) and (:,-my) are referenced but do not affect the calculations.

         if (.not.allocated(inflcf%cf)) allocate(inflcf%cf(-mx:mx-1, -my:my-1, 3, 3))
         if (size(inflcf%cf,1).ne.2*mx .or. size(inflcf%cf,2).ne.2*my) then
            deallocate(inflcf%cf)
            allocate(inflcf%cf(-mx:mx-1, -my:my-1, 3, 3))
         endif

      elseif (itypcf.eq.1) then

         ! Separate influence matrix for each y-position

         inflcf%itypcf = 1
         nmatrix = my

         ! Allocate component-array cy if not done so before, reallocate when size must be changed

         if (.not.allocated(inflcf%cy)) allocate(inflcf%cy(-mx:mx-1, -my:my-1, nmatrix, 3, 3))
         if (size(inflcf%cy,1).ne.2*mx .or. size(inflcf%cy,2).ne.2*my .or.                              &
                                                                size(inflcf%cy,3).ne.nmatrix) then
            deallocate(inflcf%cy)
            allocate(inflcf%cy(-mx:mx-1, -my:my-1, nmatrix, 3, 3))
         endif
         inflcf%cy(0,0,1,1,1) = 27d0

      else

         write(*,*) 'invalid influence coefficients type ',itypcf
         call abort_run()

      endif

      ! Default: nonzeros in blocks A_nx, A_ny, A_xn, A_xy, i.e. coupling between normal and
      !          tangential problems

      inflcf%nt_cpl = .true.

      ! Reset Fourier transform data

      if (associated(inflcf%fft_cf)) deallocate(inflcf%fft_cf)
      inflcf%fft_cf => NULL()
      inflcf%fft_mx = 2*mx
      inflcf%fft_my = 2*my
      inflcf%fft_ok = .false.

   end subroutine inflcf_new

!------------------------------------------------------------------------------------------------------------

   subroutine inflcf_copy(if_in, if_out, grid)
!--purpose: copy data of influence coefficients matrix if_in to if_out
      implicit none
!--subroutine arguments:
      type(t_inflcf)  , intent(inout) :: if_in, if_out
      type(t_grid),     target        :: grid
!--local variables:
      integer    :: nmatrix, mx, my, ik, jk, imat, iy

      mx    = grid%nx
      my    = grid%ny

      ! copy basic parameters

      if_out%xoffs    = if_in%xoffs
      if_out%yoffs    = if_in%yoffs
      if_out%itypcf   = if_in%itypcf
      if_out%nt_cpl   = if_in%nt_cpl
      if_out%ga       = if_in%ga
      if_out%ga_inv   = if_in%ga_inv
      if_out%use_3bl  = if_in%use_3bl
      if_out%flx_3bl  = if_in%flx_3bl
      if_out%use_flxz = if_in%use_flxz
      if_out%flx_z    = if_in%flx_z

      ! copy grid pointer

      if_out%grid    = if_in%grid

      ! Reset Fourier transform data

      if (associated(if_out%fft_cf)) deallocate(if_out%fft_cf)
      if_out%fft_cf  => NULL()
      if_out%fft_mx = 2*mx
      if_out%fft_my = 2*my
      if_out%fft_ok  = .false.

      ! (re)allocate and copy cf or cy-arrays

      if (if_out%itypcf.eq.0) then

         ! Standard method, same influence matrix for all elements

         ! (re)-allocate if needed

         if (.not.allocated(if_out%cf)) allocate(if_out%cf(-mx:mx-1, -my:my-1, 3, 3))
         if (size(if_out%cf,1).ne.2*mx .or. size(if_out%cf,2).ne.2*my) then
            deallocate(if_out%cf)
            allocate(if_out%cf(-mx:mx-1, -my:my-1, 3, 3))
         endif

         ! copy data

         do ik = 1, 3
            do jk = 1, 3
               do iy = -my, my-1
                  if_out%cf(-mx:mx-1,iy,ik,jk) = if_in%cf(-mx:mx-1,iy,ik,jk)
               enddo
            enddo
         enddo

      elseif (if_out%itypcf.eq.1) then

         ! Separate influence matrix for each y-position

         nmatrix = my

         ! (re)-allocate if needed

         if (.not.allocated(if_out%cy)) allocate(if_out%cy(-mx:mx-1, -my:my-1, nmatrix, 3, 3))
         if (size(if_out%cy,1).ne.2*mx .or. size(if_out%cy,2).ne.2*my .or.                              &
                                                                size(if_out%cy,3).ne.nmatrix) then
            deallocate(if_out%cy)
            allocate(if_out%cy(-mx:mx-1, -my:my-1, nmatrix, 3, 3))
         endif

         ! copy data

         do jk = 1, 3
            do ik = 1, 3
               do imat = 1, my
                  do iy = -my, my-1
                     if_out%cy(-mx:mx-1,iy,imat,ik,jk) = if_in%cy(-mx:mx-1,iy,imat,ik,jk)
                  enddo
               enddo
            enddo
         enddo

      else

         write(*,*) 'invalid influence coefficients type ',if_out%itypcf
         call abort_run()

      endif

   end subroutine inflcf_copy

!------------------------------------------------------------------------------------------------------------

   subroutine inflcf_print (inflcf, lun)
!--purpose: Print a set of influence coefficients to the output-file
      implicit none
!--subroutine arguments:
      type(t_inflcf)  , intent(in) :: inflcf
      integer         , intent(in) :: lun
!--local variables:
      logical, parameter :: fmt_cfname = .true.
      integer            :: i, ix, iy, iy0, ik, jk, imat
      integer            :: nmatrix, iaconvex, lquasid, xyzorder, lsymmx, lsymmy, izero(10)
      real(kind=8)       :: smat, sy
      character(len=14)  :: strng(9)

      associate(mx    => inflcf%grid%nx, my    => inflcf%grid%ny, dx    => inflcf%grid%dx,              &
                dy    => inflcf%grid%dy)

      ! print influence coefficients to output-file

      ! write header to file

      izero(1:10) = 0
      if (inflcf%itypcf.eq.0) then
         nmatrix  = 1
      else
         nmatrix  = my
      endif
      if (.not.fmt_cfname) then
         iaconvex = 0
         xyzorder = 123
         lsymmx   = 0
         lsymmy   = 0
         write(lun,'(i4,i5,12x,9i2,2x,2a)') nmatrix, iaconvex, izero(1:9), '% nmatrix, iaconvex: ',     &
                             'number of infl.cf.matrices, body number of convex body'
         if (inflcf%nt_cpl) then
            lquasid = 0
            write(lun,'(i4,i5,12x,9i2,2x,2a)') lquasid, xyzorder, izero(1:9), '% lquasid, xyzordr:  ',  &
                             'lquasid=0=false, general, 9 columns per point, xyzordr = 123: x,y,z.'
         else
            lquasid = 1
            write(lun,'(i4,i5,12x,9i2,2x,2a)') lquasid, xyzorder, izero(1:9), '% lquasid, xyzordr:  ',  &
                             'lquasid=1=true, uncoupled, 5 columns per point, xyzordr = 123: x,y,z.'
         endif
         write(lun,'(i4,i5,12x,9i2,2x,2a)') lsymmx, lsymmy, izero(1:9), '% lsymmx=0, lsymmy=0: ',       &
                             'no symmetry in x, no symmetry in y'
         write(lun,'(i4,i5,2f6.2,4x,7i2,2x,a)') mx, my, dx, dy, izero(1:7), '% mx, my, dx, dy'
      else
         iaconvex = 0
         lquasid  = 0
         xyzorder = 312
         lsymmx   = 1
         lsymmy   = 0
         write(lun,'(i4,i5,12x,9i2,2x,2a)') nmatrix, iaconvex, izero(1:9), '% nmatrix, iaconvex: ',     &
                             'number of infl.cf.matrices, body number of convex body'
         write(lun,'(i4,i5,12x,9i2,2x,2a)') lquasid, xyzorder, izero(1:9), '% lquasid, xyzordr:  ',     &
                             'lquasid=0=false, general, 9 columns per point, xyzordr = 312: n,x,y.'
         write(lun,'(i4,i5,12x,9i2,2x,2a)') lsymmx, lsymmy, izero(1:9), '% lsymmx=1, lsymmy=0: ',       &
                                     'symmetry in x, no symmetry in y'
         write(lun,'(i4,i5,2f6.2,4x,7i2,2x,a)') mx, my, dx, dy, izero(1:7), '% mx, my, dx, dy'
      endif

      ! write data to file

      iy0 = -my
      if (my.eq.1) iy0 = -my+1

      if (inflcf%itypcf.eq.0) then

         if (.not.fmt_cfname) then

            ! Single matrix, printing in x-y-z order

            do iy = iy0, my-1
               write(lun,*)
               write(lun,*) 'Influence coefficients for offset iy=',iy
               if (.not.inflcf%nt_cpl) then
                  write(lun,1231)
               else
                  write(lun,1233)
               endif
               do ix = -mx, mx-1
                  if (.false.) then                ! zz
                     write(lun,1232) ix, inflcf%cf(ix,iy,3,3)
                  elseif (.not.inflcf%nt_cpl) then ! xx, xy,  yx, yy,  zz
                     write(lun,1232) ix, ((inflcf%cf(ix,iy,ik,jk), jk=1,2), ik=1,2), inflcf%cf(ix,iy,3,3)
                  else                             ! xx, xy, xz,  yx, yy, yz,  zx, zy, zz
                     write(lun,1234) ix, ((inflcf%cf(ix,iy,ik,jk), jk=1,3), ik=1,3)
                  endif
               enddo
            enddo

         else

            ! Single matrix, printing in cfname-format, lsymmx=.true., lquasid=.false.

            iy0 = max(iy0, -my+1)

            do iy = iy0, my-1
               sy = iy * inflcf%grid%dy
               write(lun,'(a,i3,a,f8.3,a)') '% offset is=',iy,' (s=', sy, ')'
               write(lun,1243)
               do ix = 0, mx-1                     ! zz, zx, zy,  xz, xx, xy,  yz, yx, yy
                  strng(1) = fmt_gs(12, 5, inflcf%cf(ix,iy,3,3))
                  strng(2) = fmt_gs(12, 5, inflcf%cf(ix,iy,3,1))
                  strng(3) = fmt_gs(12, 5, inflcf%cf(ix,iy,3,2))
                  strng(4) = fmt_gs(12, 5, inflcf%cf(ix,iy,1,3))
                  strng(5) = fmt_gs(12, 5, inflcf%cf(ix,iy,1,1))
                  strng(6) = fmt_gs(12, 5, inflcf%cf(ix,iy,1,2))
                  strng(7) = fmt_gs(12, 5, inflcf%cf(ix,iy,2,3))
                  strng(8) = fmt_gs(12, 5, inflcf%cf(ix,iy,2,1))
                  strng(9) = fmt_gs(12, 5, inflcf%cf(ix,iy,2,2))
                  write(lun,1244) ix, ix*inflcf%grid%dx, (strng(i)(1:12), i=1,9)
               enddo
            enddo

 1243       format('% ix      x  ',11x,'n-displacements',23x,'x-displacements', 23x,'s-displacements',/, &
                   '%',2x,  11x, 3(2x,'due to pn',3x,'due to px',3x,'due to ps',3x))
 1244       format(i4, f9.3, 3(3a,:,'  '))

         endif

      elseif (inflcf%itypcf.eq.1) then

         ! Separate influence matrix per y-coordinate

         if (.not.fmt_cfname) then

            do imat = 1, my
               write(lun,1239)
               write(lun,*) 'Influence matrix for load at jy=',imat
               write(lun,1239)

               do iy = iy0, my-1
                  write(lun,*)
                  write(lun,*) 'Influence coefficients for offset iy=',iy
                  write(lun,1233)
                  do ix = -mx, mx-1
                     if (.false.) then
                        write(lun,1232) ix, inflcf%cy(ix,iy,imat,3,3)
                     else
                        write(lun,1234) ix, ((inflcf%cy(ix,iy,imat,ik,jk), jk=1,3), ik=1,3)
                     endif
                  enddo
               enddo
            enddo

         else

            ! printing in cfname-format, lsymmx=.true., lquasid=.false.

            iy0 = max(iy0, -my+1)

            do imat = 1, my
               smat = real(imat-my/2-1) * inflcf%grid%dy        ! ok for odd my??
               write(lun,1249)
               write(lun,'(a,i3)') '% starting data for matrix number',imat
               write(lun,1249)
               write(lun,1246) imat, 0d0, smat
  1246         format(i4, f8.3, f9.3, '     0 0 0 0 0 0 0 0  % matrix number, loaded element (xc, sc)')

               do iy = iy0, my-1
                  sy = smat + iy * inflcf%grid%dy
                  write(lun,'(a,i3,a,f8.3,a,i4,a,2(f8.3,a))') '% offset is=',iy,' (s=', sy,             &
                                ') for matrix', imat, ' (xc,sc)=(',0d0,',',smat,'),'
                  write(lun,1243)
                  do ix = 0, mx-1                     ! zz, zx, zy,  xz, xx, xy,  yz, yx, yy
                     strng(1) = fmt_gs(12, 5, inflcf%cy(ix,iy,imat,3,3))
                     strng(2) = fmt_gs(12, 5, inflcf%cy(ix,iy,imat,3,1))
                     strng(3) = fmt_gs(12, 5, inflcf%cy(ix,iy,imat,3,2))
                     strng(4) = fmt_gs(12, 5, inflcf%cy(ix,iy,imat,1,3))
                     strng(5) = fmt_gs(12, 5, inflcf%cy(ix,iy,imat,1,1))
                     strng(6) = fmt_gs(12, 5, inflcf%cy(ix,iy,imat,1,2))
                     strng(7) = fmt_gs(12, 5, inflcf%cy(ix,iy,imat,2,3))
                     strng(8) = fmt_gs(12, 5, inflcf%cy(ix,iy,imat,2,1))
                     strng(9) = fmt_gs(12, 5, inflcf%cy(ix,iy,imat,2,2))
                     write(lun,1244) ix, ix*inflcf%grid%dx, (strng(i)(1:12), i=1,9)
                  enddo
               enddo
            enddo

         endif

      else

         write(*,*) 'invalid itypcf=',inflcf%itypcf
         call abort_run()

      endif

 1231 format(13x,'x-displacements', 18x,'y-displacements', 18x,'z-displacements')
 1232 format('ix=',i6,':',2(2f16.10,:,';'), f16.9 )
 1233 format(8x,15x,'x-displacements',29x,'y-displacements', 29x,'z-displacements',/,                   &
             8x,3(4x,'due to px',5x,'due to py',5x,'due to pn',3x))
 1234 format('ix=',i4,':',3(3f16.10,:,'; '))
 1239 format(80('-'))
 1249 format(130('%'))

      end associate

   end subroutine inflcf_print

!------------------------------------------------------------------------------------------------------------

   subroutine inflcf_destroy(inflcf)
!--purpose: de-allocate space allocated for an influence coefficients matrix
      implicit none
!--subroutine arguments:
      type(t_inflcf)  , intent(inout) :: inflcf

      if (allocated(inflcf%cf)) deallocate(inflcf%cf)
      if (allocated(inflcf%cy)) deallocate(inflcf%cy)
      if (associated(inflcf%fft_cf)) deallocate(inflcf%fft_cf)
      nullify(inflcf%grid)
      nullify(inflcf%fft_cf)

   end subroutine inflcf_destroy

!------------------------------------------------------------------------------------------------------------

   subroutine influe_load (fname_influe, dirnam, is_roll, cgrid, influ)
!--purpose: Load the influence coefficients from a file.
!           type 0: one set of coefficients, independent of y-position.
!           type 1: separate matrix of coefficients per y-position.
      implicit none
!--subroutine arguments:
      logical                  :: is_roll
      character(len=*)         :: fname_influe, dirnam
      type(t_grid),     target :: cgrid
      type(t_influe),   target :: influ
!--pointers to global data items:
      integer,                    pointer :: mxarr, myarr
      real(kind=8),               pointer :: dx, dy
      real(kind=8), dimension(:,:,:,:,:), pointer :: ys, yv, ysv
!--local variables:
      integer, parameter :: mxnval = 20
      logical, parameter :: lstop  = .true.
      integer            :: linfl, nmatrix, iaconvx, iquasid, ixyzordr, isymmx, isymmy, mxfile, myfile, &
                            itypcf, ieof, ierror, nval
      real(kind=8)       :: dxi, dyi
      logical       :: flags(mxnval)
      integer       :: ints(mxnval), ncase, linenr, idebug
      real(kind=8)  :: dbles(mxnval)
      character*256 :: strngs(mxnval), fulnam_influe

      call timer_start(itimer_infload)

      ! determine full path-name, pre-pending dirname when necessary

      if (dirnam.ne.' ' .and. index_pathsep(fname_influe).le.0) then
         fulnam_influe = trim(dirnam) // path_sep // fname_influe
      else
         fulnam_influe = fname_influe
      endif

      ! use free unit number defined in m_print_output

      linfl  = ltmp

      ! set pointers to frequently used global data items

      mxarr => cgrid%nx
      myarr => cgrid%ny
      dx    => cgrid%dx
      dy    => cgrid%dy

      ncase  = 1
      idebug = 1
      ieof   = -1 ! eof=error
      write(bufout,*) 'Opening file "',trim(fulnam_influe), '" with influence coefficients'
      call write_log(1, bufout)
      open(linfl, file=fulnam_influe, status='old', err=995)
      linenr = 0
      goto 997

      ! Error handling:

 995  continue
         write(bufout,*) 'ERROR: cannot open file "', trim(fulnam_influe),'" for reading.'
         call write_log(1, bufout)
         call abort_run()
 997  continue

      ! Read & check the type of file: nmatrix, iaconvx

      call readline(linfl, ncase, linenr, 'general configuration', 'ii', ints ,dbles, flags, strngs,    &
                        mxnval, nval, idebug, ieof, lstop, ierror)
      nmatrix = ints(1)
      iaconvx = ints(2)

      if (nmatrix.eq.1) then
         itypcf = 0
      elseif (nmatrix.eq.myarr) then
         itypcf = 1
      else
         write(bufout,'(2(a,i4),a)') 'ERROR: number of matrices (',nmatrix,') should be 1 or my (',myarr,')'
         call write_log(1, bufout)
         call abort_run()
      endif

      if (iaconvx.lt.1 .or. iaconvx.gt.2) then
         write(bufout,'(a,i4,a)') 'ERROR: convex body (',iaconvx,') should be 1 or 2'
         call write_log(1, bufout)
         call abort_run()
      endif

      ! Read & check the column definition: lquasid, xyzordr

      call readline(linfl, ncase, linenr, 'column definition', 'ii', ints ,dbles, flags, strngs,        &
                        mxnval, nval, idebug, ieof, lstop, ierror)
      iquasid  = ints(1)
      ixyzordr = ints(2)

      if (iquasid.ne.0) then
         write(bufout,'(a)') 'ERROR: non-quasiidentical input should be provided, iquasid=0.'
         call write_log(1, bufout)
         call abort_run()
      endif
      if (ixyzordr.ne.312) then
         write(bufout,'(a,i4,a)') 'ERROR: the coordinate numbering (',ixyzordr,') should be 312.'
         call write_log(1, bufout)
         call abort_run()
      endif

      ! Read & check symmetries: isymmx, isymmy

      call readline(linfl, ncase, linenr, 'flag (0/1) symmetry in x,y', 'ii', ints ,dbles, flags,       &
                        strngs, mxnval, nval, idebug, ieof, lstop, ierror)
      isymmx = ints(1)
      isymmy = ints(2)
      if (isymmx.ne.1 .and. isymmy.ne.0) then
         write(bufout,'(2(a,i4),a)') 'ERROR: symmetry should be used for x (isymmx=',isymmx,            &
                ') and not for y (isymmy=', isymmy,')'
         call write_log(1, bufout)
         call abort_run()
      endif

      ! Read & check discretization parameters: mx, my

      call readline(linfl, ncase, linenr, 'discretization parameters', 'iidd', ints ,dbles, flags,      &
                        strngs, mxnval, nval, idebug, ieof, lstop, ierror)
      mxfile = ints(1)
      myfile = ints(2)
      dxi = dbles(1)
      dyi = dbles(2)

      if (mxfile.lt.mxarr .or. myfile.lt.myarr) then
         if (.false.) then
            write(bufout,*) 'ERROR: inconsistent mx,my, CONTACT input versus influe.file'
            call write_log(1, bufout)
            call abort_run()
         else
            write(bufout,*) '!!! WARNING: file of infl.cf. has smaller grid than CONTACT, padding ',    &
                            'with zeros !!!'
            call write_log(1, bufout)
         endif
      elseif (mxfile.gt.mxarr .or. myfile.gt.myarr) then
         write(bufout,*) 'Warning: file of infl.cf. has larger mx,my than the CONTACT grid'
         call write_log(1, bufout)
      endif 
      if (abs(dx-dxi).gt.1e-5 .or. abs(dy-dyi).gt.1e-5) then
         write(bufout,*) 'ERROR: inconsistent dx,dy, CONTACT input versus influe.file'
         call write_log(1, bufout)
         call abort_run()
      endif

      if (idebug.ge.5) then
         write(bufout,'(2(a,i4))') ' Infl.cf: nmatrix=',nmatrix,', iaconvx=',iaconvx
         call write_log(1, bufout)
         write(bufout,'(4(a,i4))') '          iquasid=',iquasid,', xyzordr=',ixyzordr,', isymmx,y=',     &
                isymmx, ',', isymmy
         call write_log(1, bufout)
         write(bufout,'(2(a,i4),2(a,f7.3))') '          mx,y=',mxfile,',',myfile,', dx,y=',dxi,',',dyi
         call write_log(1, bufout)
      endif

      if (idebug.ge.1) then
         if (iaconvx.eq.1) call write_log(' Infl.cf.file: the rail is the upper body.')
         if (iaconvx.eq.2) call write_log(' Infl.cf.file: the wheel is the upper body.')
      endif

      ! (re-)allocate influence coefficient arrays at the appropriate size

      call inflcf_new(influ%cs,  itypcf, cgrid)
      call inflcf_new(influ%cv,  itypcf, cgrid)
      call inflcf_new(influ%csv, itypcf, cgrid)
      call inflcf_new(influ%ms,  0,      cgrid)
      if (idebug.ge.5) write(*,*) 'infload: allocate ok.'

      ! Coupling between normal and tangential problems (?)

      influ%cs%nt_cpl = .true.
      influ%cv%nt_cpl = .true.
      influ%csv%nt_cpl = .true.

      ys    => influ%cs%cy
      yv    => influ%cv%cy
      ysv   => influ%csv%cy

      if (itypcf.eq.0) then
         if (idebug.ge.5) write(*,*) 'starting influe_loadtyp0...'
         call influe_loadtyp0(is_roll, linfl, idebug, mxfile, myfile, mxarr, myarr, influ%cs%cf,     &
                       influ%cv%cf, influ%csv%cf)
         if (idebug.ge.5) write(*,*) 'returned influe_loadtyp0...'
      else
         if (idebug.ge.5) write(*,*) 'starting influe_loadtyp1...'
         call influe_loadtyp1(is_roll, linfl, mxfile, myfile, nmatrix, mxarr, myarr, ys, yv, ysv)
      endif
      close(linfl)

      ! print influence coefficients to output-file when requested

      if (.false.) then
         if (.false.) then
            write(lout,*) 'Printing influence coefficients csv (steady)'
            call inflcf_print (influ%csv, lout)
         else
            write(lout,*) 'Printing influence coefficients cs (instat)'
            call inflcf_print (influ%cs, lout)
         endif
      endif
      if (idebug.ge.5) call write_log('infload: done, returning')
      call timer_stop(itimer_infload)

   end subroutine influe_load

!------------------------------------------------------------------------------------------------------------

   subroutine influe_loadtyp0 (is_roll, linfl, idebug, mxfile, myfile, mxarr, myarr, cs, cv, csv)
!--purpose: Actual loading of the influence coefficients for type 0, one matrix of coefficients for all
!           y-positions.
      implicit none
!--subroutine arguments:
      logical            :: is_roll
      integer            :: linfl, idebug, mxfile, myfile, mxarr, myarr
      real(kind=8), dimension(-mxarr:,-myarr:,:,:) :: cs, cv, csv
!--local variables:
      integer            :: itmp, j, ix, iy, ik, jk
      real(kind=8)       :: tmp(10)

      ! initialize coefficients in cs, all zeros

      do jk = 1, 3
         do ik = 1, 3
            do iy = -myarr, myarr-1
               do ix = -mxarr, mxarr-1
                  cs(ix,iy,ik,jk) = 0d0
               enddo
            enddo
         enddo
      enddo

      ! read data from file: no symmetry for y, offset iy = -my+1 to my-1

      if (idebug.ge.3) write(*,*) 'load_typ0: reading'

      do iy = -myfile+1, myfile-1
         if (idebug.ge.5) write(*,*) 'starting iy=',iy
         read(linfl,*,err=999)
         read(linfl,*,err=999)
         do ix = 0, mxfile-1

            ! read ix, x, 9 columns ann,anx,any, axn,axx,axy, ayn,ayx,ayy
            ! they are total values, combined for wheel and rail

            if (idebug.ge.15) write(*,*) 'reading ix=',ix
            read(linfl,*,err=999) itmp, (tmp(j), j=1,10)

            if (ix.ge.mxarr .or. iy.lt.-myarr .or. iy.ge.myarr) then

               ! (ix,iy) outside range of arrays: ignore data

            else
               ! (ix,iy) within range: store data
               ! store coefficients for offset (ix,iy) in the actual arrays

               cs( ix,iy, ikZDIR, jkZDIR) =  tmp(2)
               cs( ix,iy, ikZDIR, jkXDIR) =  tmp(3)
               cs( ix,iy, ikZDIR, jkYDIR) =  tmp(4)
               cs( ix,iy, ikXDIR, jkZDIR) =  tmp(5)
               cs( ix,iy, ikXDIR, jkXDIR) =  tmp(6)
               cs( ix,iy, ikXDIR, jkYDIR) =  tmp(7)
               cs( ix,iy, ikYDIR, jkZDIR) =  tmp(8)
               cs( ix,iy, ikYDIR, jkXDIR) =  tmp(9)
               cs( ix,iy, ikYDIR, jkYDIR) =  tmp(10)

               ! copy to offset (-ix,iy) with mirroring for components a_nx, a_xn, a_xy and a_yx

               cs(-ix,iy, ikZDIR, jkZDIR) =  tmp(2)
               cs(-ix,iy, ikZDIR, jkXDIR) = -tmp(3)
               cs(-ix,iy, ikZDIR, jkYDIR) =  tmp(4)
               cs(-ix,iy, ikXDIR, jkZDIR) = -tmp(5)
               cs(-ix,iy, ikXDIR, jkXDIR) =  tmp(6)
               cs(-ix,iy, ikXDIR, jkYDIR) = -tmp(7)
               cs(-ix,iy, ikYDIR, jkZDIR) =  tmp(8)
               cs(-ix,iy, ikYDIR, jkXDIR) = -tmp(9)
               cs(-ix,iy, ikYDIR, jkYDIR) =  tmp(10)
            endif
         enddo
      enddo

      ! Form the array of influence numbers "cv" and "csv"

      call influe_cv_csv (is_roll, mxarr, myarr, cs, cv, csv)

      return

      ! Error handling:

 999  continue
      write(*,'(a)')       ' ERROR while reading data from infl.cf-file.'
      write(*,'(3(a,i5))') '       current position iy=',iy,', ix=',ix
      call abort_run()

   end subroutine influe_loadtyp0

!------------------------------------------------------------------------------------------------------------

   subroutine influe_loadtyp1 (is_roll, linfl, mxfile, myfile, nmatrix, mxarr, myarr, ys, yv, ysv)
!--purpose: Actual loading of the influence coefficients for type 1, separate matrix of coefficients per y-position.
      implicit none
!--subroutine arguments:
      logical       :: is_roll
      integer       :: linfl, mxarr, myarr, mxfile, myfile, nmatrix
      real(kind=8), dimension(-mxarr:mxarr-1,-myarr:myarr-1,nmatrix,3,3) :: ys, yv, ysv
!--local variables:
      integer            :: idebug = 0
      integer            :: ilin, itmp, j, imat, ix, iy, ik, jk
      real(kind=8)       :: tmp(10)

      if (idebug.ge.5) call write_log('influe_loadtyp1: starting...')
      do jk = 1, 3
         do ik = 1, 3
            do imat = 1, nmatrix
               do iy = -myarr, myarr-1
                  do ix = -mxarr, mxarr-1
                     ys(ix,iy,imat,ik,jk) = 0d0
                  enddo
               enddo
            enddo
         enddo
      enddo
      if (idebug.ge.5) call write_log('influe_loadtyp1: initialize ys ok.')

      do imat = 1, nmatrix
         ! write(*,*) 'influe_loadtyp1: starting for imat=',imat

         ! skip header (3 lines), matrix id

         iy = 0
         ix = -2
         do ilin = 1, 4
            read(linfl,*,err=999)
         enddo

         ! no symmetry for y, offset iy = -my+1 to my-1

         do iy = -myfile+1, myfile-1
            ! write(*,*) 'influe_loadtyp1: imat=',imat,': iy=',iy

            ! skip header for offset iy (2 lines)

            ix = -1
            read(linfl,*,err=999)
            read(linfl,*,err=999)

            ! read data from file: using symmetry for x, lines ix = 0 to mx-1

            do ix = 0, mxfile-1

               ! read ix, x, 9 columns ann,anx,any, axn,axx,axy, ayn,ayx,ayy
               ! they are total values, combined for wheel and rail

               read(linfl,*,err=999) itmp, (tmp(j), j=1,10)

               if (ix.ge.mxarr .or. iy.lt.-myarr .or. iy.ge.myarr) then

                  ! ignore data from file for (ix,iy) outside the size of the actual arrays

               else

                  ! store coefficients for offset (ix,iy) in the actual arrays

                  ys( ix,iy, imat, ikZDIR, jkZDIR) =  tmp(2)
                  ys( ix,iy, imat, ikZDIR, jkXDIR) =  tmp(3)
                  ys( ix,iy, imat, ikZDIR, jkYDIR) =  tmp(4)
                  ys( ix,iy, imat, ikXDIR, jkZDIR) =  tmp(5)
                  ys( ix,iy, imat, ikXDIR, jkXDIR) =  tmp(6)
                  ys( ix,iy, imat, ikXDIR, jkYDIR) =  tmp(7)
                  ys( ix,iy, imat, ikYDIR, jkZDIR) =  tmp(8)
                  ys( ix,iy, imat, ikYDIR, jkXDIR) =  tmp(9)
                  ys( ix,iy, imat, ikYDIR, jkYDIR) =  tmp(10)

                  ! copy to offset (-ix,iy) with mirroring for components a_nx, a_xn, a_xy and a_yx

                  ys(-ix,iy, imat, ikZDIR, jkZDIR) =  tmp(2)
                  ys(-ix,iy, imat, ikZDIR, jkXDIR) = -tmp(3)
                  ys(-ix,iy, imat, ikZDIR, jkYDIR) =  tmp(4)
                  ys(-ix,iy, imat, ikXDIR, jkZDIR) = -tmp(5)
                  ys(-ix,iy, imat, ikXDIR, jkXDIR) =  tmp(6)
                  ys(-ix,iy, imat, ikXDIR, jkYDIR) = -tmp(7)
                  ys(-ix,iy, imat, ikYDIR, jkZDIR) =  tmp(8)
                  ys(-ix,iy, imat, ikYDIR, jkXDIR) = -tmp(9)
                  ys(-ix,iy, imat, ikYDIR, jkYDIR) =  tmp(10)
               endif

            enddo
         enddo
      enddo
      if (idebug.ge.5) call write_log('influe_loadtyp1: done reading')

      ! Form the arrays of influence numbers "yv" and "ysv"

      call influe_yv_ysv (is_roll, nmatrix, mxarr, myarr, ys, yv, ysv)

      return

      ! Error handling:

 999  continue
      write(*,'(a)')       'ERROR while reading data from infl.cf-file.'
      write(*,'(3(a,i5))') '      current position imat=',imat, ', iy=',iy,', ix=',ix
      call abort_run()

   end subroutine influe_loadtyp1

!------------------------------------------------------------------------------------------------------------

   subroutine influe_blanco (is_roll, imethod, iversion, nincln, aincln, cgrid, influ)
!--purpose: Compute influence coefficients using the Blanco correction
      implicit none
!--subroutine arguments:
      logical                   :: is_roll
      integer,       intent(in) :: imethod, iversion, nincln
      real(kind=8),  intent(in) :: aincln(nincln,2)
      type(t_grid),      target :: cgrid
      type(t_influe),    target :: influ
!--local variables:
      integer, parameter :: idebug = 0
      integer            :: iy, ii, j, mx, my, itypcf, nmatrix
      real(kind=8)       :: s0, s1, si, a0, a1, ata(2,2), atb(2), c0, c1, det
      type(t_inflcf)     :: aij
      real(kind=8), dimension(:), allocatable :: alphai

      call timer_start(itimer_sgencr)

      mx = cgrid%nx
      my = cgrid%ny

      if (imethod.lt.0 .or. imethod.gt.1) then
         call write_log('INTERNAL ERROR(influe_blanco): needs 0 <= imethod <= 1')
         return
      endif
      if (iversion.lt.1 .or. iversion.gt.4) then
         call write_log('INTERNAL ERROR(influe_blanco): needs 1 <= iversion <= 4')
         return
      endif

      if (nincln.le.1) then
         call write_log('INTERNAL ERROR(influe_blanco): needs nincln>=2')
         return
      elseif (nincln.le.2 .or. imethod.eq.0) then
         if (idebug.ge.1) call write_log(' influe_blanco: linear alpha, using itypcf=0')
         itypcf  = 0
         nmatrix = 1
      else
         if (idebug.ge.1) call write_log(' influe_blanco: nincln>2, using itypcf=1')
         itypcf  = 1
         nmatrix = cgrid%ny
      endif

      ! get surface inclinations alphai(iy) at grid positions iy

      allocate(alphai(my))

      if (imethod.eq.0) then

         ! fast method: construct linear approximation of surface inclinations aincln

         ata(1,1) = real(nincln)
         ata(1,2) = sum( aincln(:,1) )
         ata(2,1) = ata(1,2)
         ata(2,2) = sum( aincln(:,1)**2 )
         atb(1)   = sum( aincln(:,2) )
         atb(2)   = sum( aincln(:,1)*aincln(:,2) )
         det      = ata(1,1)*ata(2,2) - ata(2,1)*ata(1,2)
         c0       = ( ata(2,2)*atb(1) - ata(1,2)*atb(2)) / det
         c1       = (-ata(2,1)*atb(1) + ata(1,1)*atb(2)) / det

         ! evaluate linear function

         do iy = 1, my
            ii = 1 + (iy-1)*mx
            alphai(iy) = c0 + c1 * cgrid%y(ii)
         enddo

      else

         ! full method: interpolate surface inclinations aincln to grid positions iy

         do iy = 1, my
            ii = 1 + (iy-1)*mx
            si = cgrid%y(ii)

            if (si.le.aincln(2,1)) then            ! <start: linear extrapolation
               j  = 1
            elseif (si.ge.aincln(nincln-1,1)) then ! >end: linear extrapolation
               j  = nincln-1
            else                                   ! interior: linear interpolation
               j  = 1
               do while(si.ge.aincln(j+1,1))
                  j = j + 1
               enddo
            endif

            s0 = aincln(j  ,1)
            s1 = aincln(j+1,1)
            a0 = aincln(j  ,2)
            a1 = aincln(j+1,2)
            alphai(iy) = a0 + (si - s0) * (a1 - a0) / (s1 - s0)

            if (idebug.ge.5) then
               write(bufout,'(a,i2,a,f6.2,a,i2,2(a,2f6.2),2(a,f6.2))') ' iy=',iy,', s=',si,': j=',j,    &
                   ', s0,a0=',s0,a0, ', s1,a1=',s1,a1,', fac=', (si-s0)/(s1-s0),', ai=',alphai(iy)
               call write_log(1,bufout)
            endif
         enddo

      endif ! imethod

      if (idebug.ge.5) then
         do iy = 1, my
            ii = 1 + (iy-1)*mx
            write(bufout,'(a,i3,2(a,f8.4))') ' iy=',iy,', s=',cgrid%y(ii),': ai=',alphai(iy)
            call write_log(1,bufout)
         enddo
      endif

      ! copy halfspace influence coefficients to temporary storage

      call inflcf_new(aij, 0, cgrid)
      call inflcf_copy(influ%cs, aij, cgrid)

      ! (re-)allocate influence coefficient arrays at the appropriate size

      call inflcf_new(influ%cs,  itypcf, cgrid)
      call inflcf_new(influ%cv,  itypcf, cgrid)
      call inflcf_new(influ%csv, itypcf, cgrid)
      call inflcf_new(influ%ms,  0,      cgrid)
      if (idebug.ge.5) write(*,*) 'infload: allocate ok.'

      ! Coupling between normal and tangential problems (?)

      influ%cs%nt_cpl = .true.
      influ%cv%nt_cpl = .true.
      influ%csv%nt_cpl = .true.

      if (itypcf.eq.0) then
         ! compute cs

         call influe_blanco0(aij%cf, iversion, alphai, mx, my, influ%cs%cf)

         ! form cv and csv

         call influe_cv_csv (is_roll, mx, my, influ%cs%cf, influ%cv%cf, influ%csv%cf)
      else
         ! compute ys

         call influe_blanco1(aij%cf, iversion, alphai, nmatrix, mx, my, influ%cs%cy)

         ! form yv and ysv

         call influe_yv_ysv (is_roll, nmatrix, mx, my, influ%cs%cy, influ%cv%cy, influ%csv%cy)
      endif

      ! print influence coefficients to output-file when requested

      if (.false.) then
         if (.false.) then
            write(lout,*) 'Printing influence coefficients csv (steady)'
            call inflcf_print (influ%csv, lout)
         else
            write(lout,*) 'Printing influence coefficients cs (instat)'
            call inflcf_print (influ%cs, lout)
         endif
      endif

      call inflcf_destroy(aij)

      if (idebug.ge.5) call write_log('infload: done, returning')
      call timer_stop(itimer_sgencr)

   end subroutine influe_blanco

!------------------------------------------------------------------------------------------------------------

   subroutine influe_blanco0 (aij, iversion, alphai, mx, my, bij)
!--purpose: Actual computation of the Blanco IF-correction for type 0, using one matrix Bij for all
!           y-positions.
      implicit none
!--subroutine arguments:
      integer            :: iversion, mx, my
      real(kind=8)       :: alphai(my)
      real(kind=8)       :: aij(-mx:mx-1,-my:my-1,3,3)
      real(kind=8)       :: bij(-mx:mx-1,-my:my-1,3,3)
!--local variables:
      integer            :: idebug = 0
      integer            :: ix, iy, ik, jk
      real(kind=8)       :: alpha, cs, sn, csb, snb, aspn, asmn

      if (idebug.ge.5) call write_log('influe_blanco0: starting, bij=0...')

      do jk = 1, 3
         do ik = 1, 3
            do iy = -my, my-1
               do ix = -mx, mx-1
                  bij(ix,iy,ik,jk) = 0d0
               enddo
            enddo
         enddo
      enddo

      if (idebug.ge.5) call write_log('influe_blanco0: initialize bij ok.')

      ! jy = 0, ky = iy
      do iy = -my+1, my-1

         if (iy.ge.0) then
            alpha = alphai( 1) - alphai( 1+iy)
         else
            alpha = alphai(my) - alphai(my+iy)
         endif
         cs    = cos(alpha)
         sn    = sin(alpha)
         csb   = cos(0.5d0*alpha)
         snb   = sin(0.5d0*alpha)

         do ix = -mx+1, mx-1

            if (iversion.eq.1) then

               ! version 1, original Blanco formulas

               ! column 1: Bxx, Bsx, Bnx
               bij(ix,iy, ikXDIR,jkXDIR) = aij(ix,iy, ikXDIR,jkXDIR)
               bij(ix,iy, ikYDIR,jkXDIR) = aij(ix,iy, ikYDIR,jkXDIR)
               bij(ix,iy, ikZDIR,jkXDIR) = aij(ix,iy, ikZDIR,jkXDIR)

               ! column 2: Bxs, Bss, Bns
               bij(ix,iy, ikXDIR,jkYDIR) = aij(ix,iy, ikXDIR,jkYDIR) * cs + aij(ix,iy, ikXDIR,jkZDIR) * sn
               bij(ix,iy, ikYDIR,jkYDIR) = aij(ix,iy, ikYDIR,jkYDIR) * cs + aij(ix,iy, ikYDIR,jkZDIR) * sn
               bij(ix,iy, ikZDIR,jkYDIR) = aij(ix,iy, ikZDIR,jkYDIR) * cs + aij(ix,iy, ikZDIR,jkZDIR) * sn

               ! column 3: Bxn, Bsn, Bnn
               bij(ix,iy, ikXDIR,jkZDIR) = aij(ix,iy, ikXDIR,jkZDIR) * cs - aij(ix,iy, ikXDIR,jkYDIR) * sn
               bij(ix,iy, ikYDIR,jkZDIR) = aij(ix,iy, ikYDIR,jkZDIR) * cs - aij(ix,iy, ikYDIR,jkYDIR) * sn
               bij(ix,iy, ikZDIR,jkZDIR) = aij(ix,iy, ikZDIR,jkZDIR) * cs - aij(ix,iy, ikZDIR,jkYDIR) * sn

            elseif (iversion.eq.2) then

               ! version 2, Blanco formulas, without xn correction

               ! column 1: Bxx, Bsx, Bnx
               bij(ix,iy, ikXDIR,jkXDIR) = aij(ix,iy, ikXDIR,jkXDIR)
               bij(ix,iy, ikYDIR,jkXDIR) = aij(ix,iy, ikYDIR,jkXDIR)
               bij(ix,iy, ikZDIR,jkXDIR) = aij(ix,iy, ikZDIR,jkXDIR)

               ! column 2: Bxs, Bss, Bns
               bij(ix,iy, ikXDIR,jkYDIR) = aij(ix,iy, ikXDIR,jkYDIR) * cs + aij(ix,iy, ikXDIR,jkZDIR) * sn
               bij(ix,iy, ikYDIR,jkYDIR) = aij(ix,iy, ikYDIR,jkYDIR) * cs + aij(ix,iy, ikYDIR,jkZDIR) * sn
               bij(ix,iy, ikZDIR,jkYDIR) = aij(ix,iy, ikZDIR,jkYDIR) * cs + aij(ix,iy, ikZDIR,jkZDIR) * sn

               ! column 3: Bxn, Bsn, Bnn
               bij(ix,iy, ikXDIR,jkZDIR) = aij(ix,iy, ikXDIR,jkZDIR)
               bij(ix,iy, ikYDIR,jkZDIR) = aij(ix,iy, ikYDIR,jkZDIR) * cs - aij(ix,iy, ikYDIR,jkYDIR) * sn
               bij(ix,iy, ikZDIR,jkZDIR) = aij(ix,iy, ikZDIR,jkZDIR) * cs - aij(ix,iy, ikZDIR,jkYDIR) * sn

            elseif (iversion.eq.3) then

               ! version 3, halfspace at J, uI = R * Aij * pJ

               ! row 1: Bxx, Bxs, Bxn
               bij(ix,iy, ikXDIR,jkXDIR) = aij(ix,iy, ikXDIR,jkXDIR)
               bij(ix,iy, ikXDIR,jkYDIR) = aij(ix,iy, ikXDIR,jkYDIR)
               bij(ix,iy, ikXDIR,jkZDIR) = aij(ix,iy, ikXDIR,jkZDIR)

               ! row 2: Bsx, Bss, Bsn
               bij(ix,iy, ikYDIR,jkXDIR) = aij(ix,iy, ikYDIR,jkXDIR) * cs - aij(ix,iy, ikZDIR,jkXDIR) * sn
               bij(ix,iy, ikYDIR,jkYDIR) = aij(ix,iy, ikYDIR,jkYDIR) * cs - aij(ix,iy, ikZDIR,jkYDIR) * sn
               bij(ix,iy, ikYDIR,jkZDIR) = aij(ix,iy, ikYDIR,jkZDIR) * cs - aij(ix,iy, ikZDIR,jkZDIR) * sn

               ! row 3: Bnx, Bns, Bnn
               bij(ix,iy, ikZDIR,jkXDIR) = aij(ix,iy, ikYDIR,jkXDIR) * sn + aij(ix,iy, ikZDIR,jkXDIR) * cs
               bij(ix,iy, ikZDIR,jkYDIR) = aij(ix,iy, ikYDIR,jkYDIR) * sn + aij(ix,iy, ikZDIR,jkYDIR) * cs
               bij(ix,iy, ikZDIR,jkZDIR) = aij(ix,iy, ikYDIR,jkZDIR) * sn + aij(ix,iy, ikZDIR,jkZDIR) * cs

            elseif (iversion.eq.4) then

               ! version 4, splitted formula, uI = R(beta) * Aij * R(beta) * pJ

               ! row 1: Bxx, Bxs, Bxn
               bij(ix,iy, ikXDIR,jkXDIR) = aij(ix,iy, ikXDIR,jkXDIR)
               bij(ix,iy, ikXDIR,jkYDIR) = aij(ix,iy, ikXDIR,jkYDIR) * csb + aij(ix,iy, ikXDIR,jkZDIR) * snb
               bij(ix,iy, ikXDIR,jkZDIR) = aij(ix,iy, ikXDIR,jkZDIR) * csb - aij(ix,iy, ikXDIR,jkYDIR) * snb

               aspn = 0.5d0 * (aij(ix,iy, ikYDIR,jkYDIR) + aij(ix,iy, ikZDIR,jkZDIR))
               asmn = 0.5d0 * (aij(ix,iy, ikYDIR,jkYDIR) - aij(ix,iy, ikZDIR,jkZDIR))

               ! row 2: Bsx, Bss, Bsn
               bij(ix,iy, ikYDIR,jkXDIR) = aij(ix,iy, ikXDIR,jkYDIR) * csb + aij(ix,iy, ikXDIR,jkZDIR) * snb
               bij(ix,iy, ikYDIR,jkYDIR) =         asmn    +    aspn * cs  + aij(ix,iy, ikYDIR,jkZDIR) * sn
               bij(ix,iy, ikYDIR,jkZDIR) = aij(ix,iy, ikYDIR,jkZDIR) * cs  -            aspn           * sn

               ! row 3: Bnx, Bns, Bnn
               bij(ix,iy, ikZDIR,jkXDIR) = aij(ix,iy, ikXDIR,jkYDIR) * snb - aij(ix,iy, ikXDIR,jkZDIR) * csb
               bij(ix,iy, ikZDIR,jkYDIR) =                      aspn * sn  - aij(ix,iy, ikYDIR,jkZDIR) * cs
               bij(ix,iy, ikZDIR,jkZDIR) =       - asmn    +    aspn * cs  + aij(ix,iy, ikYDIR,jkZDIR) * sn

            else
               call write_log('INTERNAL ERROR(blanco0): incorrect iversion')
               call abort_run()
            endif

         enddo ! ix
      enddo ! iy

      if (idebug.ge.5) call write_log('influe_blanco0: done computing')

   end subroutine influe_blanco0 

!------------------------------------------------------------------------------------------------------------

   subroutine influe_blanco1 (aij, iversion, alphai, nmatrix, mx, my, bij)
!--purpose: Actual computation of the Blanco IF-correction for type 1, using a separate matrix of
!           coefficients per y-position.
      implicit none
!--subroutine arguments:
      integer            :: iversion, mx, my, nmatrix
      real(kind=8)       :: alphai(my)
      real(kind=8)       :: aij(-mx:mx-1,-my:my-1,        3,3)
      real(kind=8)       :: bij(-mx:mx-1,-my:my-1,nmatrix,3,3)
!--local variables:
      integer            :: idebug = 0
      integer            :: imat, ix, iy, jy, ky, ik, jk
      real(kind=8)       :: alpha, cs, sn, csb, snb, aspn, asmn

      if (idebug.ge.5) call write_log('influe_blanco1: starting, bij=0...')

      do jk = 1, 3
         do ik = 1, 3
            do imat = 1, nmatrix
               do iy = -my, my-1
                  do ix = -mx, mx-1
                     bij(ix,iy,imat,ik,jk) = 0d0
                  enddo
               enddo
            enddo
         enddo
      enddo

      if (idebug.ge.5) call write_log('influe_blanco1: initialize bij ok.')

      do imat = 1, nmatrix
         jy = imat
         do iy = 1, my

            ky    = iy - jy
            alpha = alphai(jy) - alphai(iy)
            cs    = cos(alpha)
            sn    = sin(alpha)
            csb   = cos(0.5d0*alpha)
            snb   = sin(0.5d0*alpha)

            do ix = -mx+1, mx-1

               if (iversion.eq.1) then

                  ! version 1, original Blanco formulas

                  ! Bxx, Bsx, Bnx
                  bij(ix,ky, imat, ikXDIR,jkXDIR) = aij(ix,ky, ikXDIR,jkXDIR)
                  bij(ix,ky, imat, ikYDIR,jkXDIR) = aij(ix,ky, ikYDIR,jkXDIR)
                  bij(ix,ky, imat, ikZDIR,jkXDIR) = aij(ix,ky, ikZDIR,jkXDIR)

                  ! Bxs, Bss, Bns
                  bij(ix,ky, imat, ikXDIR,jkYDIR) = aij(ix,ky, ikXDIR,jkYDIR) * cs + aij(ix,ky, ikXDIR,jkZDIR) * sn
                  bij(ix,ky, imat, ikYDIR,jkYDIR) = aij(ix,ky, ikYDIR,jkYDIR) * cs + aij(ix,ky, ikYDIR,jkZDIR) * sn
                  bij(ix,ky, imat, ikZDIR,jkYDIR) = aij(ix,ky, ikZDIR,jkYDIR) * cs + aij(ix,ky, ikZDIR,jkZDIR) * sn

                  ! Bxn, Bsn, Bnn
                  bij(ix,ky, imat, ikXDIR,jkZDIR) = aij(ix,ky, ikXDIR,jkZDIR) * cs - aij(ix,ky, ikXDIR,jkYDIR) * sn
                  bij(ix,ky, imat, ikYDIR,jkZDIR) = aij(ix,ky, ikYDIR,jkZDIR) * cs - aij(ix,ky, ikYDIR,jkYDIR) * sn
                  bij(ix,ky, imat, ikZDIR,jkZDIR) = aij(ix,ky, ikZDIR,jkZDIR) * cs - aij(ix,ky, ikZDIR,jkYDIR) * sn

               elseif (iversion.eq.2) then

                  ! version 2, Blanco formulas, without xn correction

                  ! column 1: Bxx, Bsx, Bnx
                  bij(ix,ky,imat, ikXDIR,jkXDIR) = aij(ix,ky, ikXDIR,jkXDIR)
                  bij(ix,ky,imat, ikYDIR,jkXDIR) = aij(ix,ky, ikYDIR,jkXDIR)
                  bij(ix,ky,imat, ikZDIR,jkXDIR) = aij(ix,ky, ikZDIR,jkXDIR)

                  ! column 2: Bxs, Bss, Bns
                  bij(ix,ky,imat, ikXDIR,jkYDIR) = aij(ix,ky, ikXDIR,jkYDIR) * cs + aij(ix,ky, ikXDIR,jkZDIR) * sn
                  bij(ix,ky,imat, ikYDIR,jkYDIR) = aij(ix,ky, ikYDIR,jkYDIR) * cs + aij(ix,ky, ikYDIR,jkZDIR) * sn
                  bij(ix,ky,imat, ikZDIR,jkYDIR) = aij(ix,ky, ikZDIR,jkYDIR) * cs + aij(ix,ky, ikZDIR,jkZDIR) * sn

                  ! column 3: Bxn, Bsn, Bnn
                  bij(ix,ky,imat, ikXDIR,jkZDIR) = aij(ix,ky, ikXDIR,jkZDIR)
                  bij(ix,ky,imat, ikYDIR,jkZDIR) = aij(ix,ky, ikYDIR,jkZDIR) * cs - aij(ix,ky, ikYDIR,jkYDIR) * sn
                  bij(ix,ky,imat, ikZDIR,jkZDIR) = aij(ix,ky, ikZDIR,jkZDIR) * cs - aij(ix,ky, ikZDIR,jkYDIR) * sn

               elseif (iversion.eq.3) then

                  ! version 3, halfspace at J, uI = R * Aij * pJ

                  ! row 1: Bxx, Bxs, Bxn
                  bij(ix,ky,imat, ikXDIR,jkXDIR) = aij(ix,ky, ikXDIR,jkXDIR)
                  bij(ix,ky,imat, ikXDIR,jkYDIR) = aij(ix,ky, ikXDIR,jkYDIR)
                  bij(ix,ky,imat, ikXDIR,jkZDIR) = aij(ix,ky, ikXDIR,jkZDIR)

                  ! row 2: Bsx, Bss, Bsn
                  bij(ix,ky,imat, ikYDIR,jkXDIR) = aij(ix,ky, ikYDIR,jkXDIR) * cs - aij(ix,ky, ikZDIR,jkXDIR) * sn
                  bij(ix,ky,imat, ikYDIR,jkYDIR) = aij(ix,ky, ikYDIR,jkYDIR) * cs - aij(ix,ky, ikZDIR,jkYDIR) * sn
                  bij(ix,ky,imat, ikYDIR,jkZDIR) = aij(ix,ky, ikYDIR,jkZDIR) * cs - aij(ix,ky, ikZDIR,jkZDIR) * sn

                  ! row 3: Bnx, Bns, Bnn
                  bij(ix,ky,imat, ikZDIR,jkXDIR) = aij(ix,ky, ikYDIR,jkXDIR) * sn + aij(ix,ky, ikZDIR,jkXDIR) * cs
                  bij(ix,ky,imat, ikZDIR,jkYDIR) = aij(ix,ky, ikYDIR,jkYDIR) * sn + aij(ix,ky, ikZDIR,jkYDIR) * cs
                  bij(ix,ky,imat, ikZDIR,jkZDIR) = aij(ix,ky, ikYDIR,jkZDIR) * sn + aij(ix,ky, ikZDIR,jkZDIR) * cs

               elseif (iversion.eq.4) then

                  ! version 4, splitted formula, uI = R(beta) * Aij * R(beta) * pJ

                  ! row 1: Bxx, Bxs, Bxn
                  bij(ix,ky,imat, ikXDIR,jkXDIR) = aij(ix,ky, ikXDIR,jkXDIR)
                  bij(ix,ky,imat, ikXDIR,jkYDIR) = aij(ix,ky, ikXDIR,jkYDIR) * csb + aij(ix,ky, ikXDIR,jkZDIR) * snb
                  bij(ix,ky,imat, ikXDIR,jkZDIR) = aij(ix,ky, ikXDIR,jkZDIR) * csb - aij(ix,ky, ikXDIR,jkYDIR) * snb

                  aspn = 0.5d0 * (aij(ix,ky, ikYDIR,jkYDIR) + aij(ix,ky, ikZDIR,jkZDIR))
                  asmn = 0.5d0 * (aij(ix,ky, ikYDIR,jkYDIR) - aij(ix,ky, ikZDIR,jkZDIR))

                  ! row 2: Bsx, Bss, Bsn
                  bij(ix,ky,imat, ikYDIR,jkXDIR) = aij(ix,ky, ikXDIR,jkYDIR) * csb + aij(ix,ky, ikXDIR,jkZDIR) * snb
                  bij(ix,ky,imat, ikYDIR,jkYDIR) =         asmn    +    aspn * cs  + aij(ix,ky, ikYDIR,jkZDIR) * sn
                  bij(ix,ky,imat, ikYDIR,jkZDIR) = aij(ix,ky, ikYDIR,jkZDIR) * cs  -            aspn           * sn

                  ! row 3: Bnx, Bns, Bnn
                  bij(ix,ky,imat, ikZDIR,jkXDIR) = aij(ix,ky, ikXDIR,jkYDIR) * snb - aij(ix,ky, ikXDIR,jkZDIR) * csb
                  bij(ix,ky,imat, ikZDIR,jkYDIR) =                      aspn * sn  - aij(ix,ky, ikYDIR,jkZDIR) * cs
                  bij(ix,ky,imat, ikZDIR,jkZDIR) =       - asmn    +    aspn * cs  + aij(ix,ky, ikYDIR,jkZDIR) * sn

               else
                  call write_log('INTERNAL ERROR(blanco0): incorrect iversion')
                  call abort_run()
               endif

            enddo ! ix
         enddo ! iy
      enddo ! imat

      if (idebug.ge.5) call write_log('influe_blanco1: done computing')

   end subroutine influe_blanco1 

!------------------------------------------------------------------------------------------------------------

   subroutine influe_cv_csv (is_roll, mx, my, cs, cv, csv)
!--purpose: Construct influence coefficients cv, csv for the given coefficients cs
      implicit none
!--subroutine arguments:
      logical            :: is_roll
      integer            :: mx, my
      real(kind=8), dimension(-mx:,-my:,:,:) :: cs, cv, csv
!--local variables:
      integer            :: ix, iy, ik, jk

      ! 2) Form the array of influence numbers "cv".
      !    The coefficients cf(ix-jx,iy-jy,:,:) say what the deformation is in element (ix,iy) of the
      !    grid of the current time, when a unit pressure is located in the element (jx,jy) of the grid
      !    of the previous time.

      if (.not.is_roll) then

         ! in shifts, the coordinate system is world-fixed, the current grid equals the previous grid,
         ! therefore the displacements can be computed using cs.

         do jk = 1, 3
            do ik = 1, 3
               do ix = -mx, mx-1
                  do iy = -my, my-1
                     cv(ix,iy,ik,jk) = cs(ix,iy,ik,jk)
                  enddo
               enddo
            enddo
         enddo

      else

         ! in rolling, the coordinate system moves with the contact patch with direction chi and
         ! distance traversed per step dq.

         do jk = 1, 3
            do ik = 1, 3
               do ix = -mx, mx-2
                  do iy = -my, my-1
                     cv(ix,iy,ik,jk) = cs(ix+1,iy,ik,jk)
                  enddo
               enddo
            enddo
         enddo

      endif

      ! 3) Form the array of influence numbers "csv".
      !    These coefficients say what the increment of the deformation is in element (ix,iy) between
      !    the previous and current time when a unit pressure is located in the element (jx,jy) of the
      !    current grid and of the grid of the previous time.

      if (is_roll) then

         ! rolling:

         do jk = 1, 3
            do iy = -my, my-1
               do ix = -mx, mx-1
                  csv(ix,iy,ikXDIR,jk) = cs(ix,iy,ikXDIR,jk) - cv(ix,iy,ikXDIR,jk)
                  csv(ix,iy,ikYDIR,jk) = cs(ix,iy,ikYDIR,jk) - cv(ix,iy,ikYDIR,jk)
                  csv(ix,iy,ikZDIR,jk) = cs(ix,iy,ikZDIR,jk)
               enddo
            enddo
         enddo

      else

         ! shifts:

         do jk = 1, 3
            do ik = 1, 3
               do ix = -mx, mx-1
                  do iy = -my, my-1
                     csv(ix,iy,ik,jk) = cs(ix,iy,ik,jk)
                  enddo
               enddo
            enddo
         enddo

      endif

   end subroutine influe_cv_csv

!------------------------------------------------------------------------------------------------------------

   subroutine influe_yv_ysv (is_roll, nmatrix, mx, my, ys, yv, ysv)
!--purpose: Actual loading of the influence coefficients for type 1, separate matrix of coefficients per y-position.
      implicit none
!--subroutine arguments:
      logical       :: is_roll
      integer       :: nmatrix, mx, my
      real(kind=8), dimension(-mx:mx-1,-my:my-1,nmatrix,3,3) :: ys, yv, ysv
!--local variables:
      integer            :: imat, ix, iy, ik, jk

      ! Form the array of influence numbers "yv".
      !    The coefficients yv(ix-jx,iy-jy,:,:,:) say what the deformation is in element (ix,iy) of the
      !    grid of the current time, when a unit pressure is located in the element (jx,jy) of the grid
      !    of the previous time.

      if (.not.is_roll) then

         ! in shifts, the coordinate system is world-fixed, the current grid equals the previous grid,
         ! therefore the displacements can be computed using ys.

         do jk = 1, 3
            do ik = 1, 3
               do imat = 1, nmatrix
                  do ix = -mx, mx-1
                     do iy = -my, my-1
                        yv(ix,iy,imat,ik,jk) = ys(ix,iy,imat,ik,jk)
                     enddo
                  enddo
               enddo
            enddo
         enddo

      else

         ! in rolling, the coordinate system moves with the contact patch with direction chi and distance
         ! traversed per step dq.

         do jk = 1, 3
            do ik = 1, 3
               do imat = 1, nmatrix
                  do ix = -mx, mx-2
                     do iy = -my, my-1
                        yv(ix,iy,imat,ik,jk) = ys(ix+1,iy,imat,ik,jk)
                     enddo
                  enddo
               enddo
            enddo
         enddo

      endif

      ! 3) Form the array of influence numbers "ysv".
      !    These coefficients say what the increment of the deformation is in element (ix,iy) between the
      !    previous and current time when a unit pressure is located in the element (jx,jy) of the current
      !    grid and of the grid of the previous time.

      if (is_roll) then

         ! rolling:

         do jk = 1, 3
            do imat = 1, nmatrix
               do iy = -my, my-1
                  do ix = -mx, mx-1
                     ysv(ix,iy,imat,ikXDIR,jk) = ys(ix,iy,imat,ikXDIR,jk) - yv(ix,iy,imat,ikXDIR,jk)
                     ysv(ix,iy,imat,ikYDIR,jk) = ys(ix,iy,imat,ikYDIR,jk) - yv(ix,iy,imat,ikYDIR,jk)
                     ysv(ix,iy,imat,ikZDIR,jk) = ys(ix,iy,imat,ikZDIR,jk)
                  enddo
               enddo
            enddo
         enddo

      else

         ! shifts:

         do jk = 1, 3
            do ik = 1, 3
               do imat = 1, nmatrix
                  do ix = -mx, mx-1
                     do iy = -my, my-1
                        ysv(ix,iy,imat,ik,jk) = ys(ix,iy,imat,ik,jk)
                     enddo
                  enddo
               enddo
            enddo
         enddo

      endif

   end subroutine influe_yv_ysv

!------------------------------------------------------------------------------------------------------------

end module m_inflcf
