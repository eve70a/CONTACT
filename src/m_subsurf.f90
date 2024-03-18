!============================================================================================================
! m_subsurf - computation of the subsurface stresses for given surface traction
!
! Copyright 1979-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!============================================================================================================
module m_subsurf

use m_hierarch_data
use m_aijpj
#ifdef _OPENMP
   use omp_lib         , only : omp_get_num_threads, omp_get_thread_num
#endif

! provide for compilation with "internal parallelization" (within each contact patch) enabled or disabled
#if defined OMP_INSIDE
#define _INT_OMP $omp
#else
#define _INT_OMP 
#endif

implicit none
private

public  rdsubs
public  wrsubs
public  subsur
public  subsur_calc
public  subsur_matfil
public  subsur_outfile

private expand_xb
private sstres_fft
private sstres_inflcf
private sstres
private sstres_derived
private stres1_pcwcns
private stres1_bilin
private dsyevc3
private Hw
private Gw
private iFw
private Fw
public  unique_positive
public  bubsrt_int
public  bubsrt_real

type t_xyzr
   real(kind=8) :: y(-1:1,-1:1, 1:3)
   real(kind=8) :: r(-1:1,-1:1, 0:3)
   real(kind=8) :: r3(-1:1,-1:1, 0:3)
   real(kind=8) :: logx(-1:1,-1:1, 0:3)
   real(kind=8) :: logy(-1:1,-1:1, 0:3)
   real(kind=8) :: logz(-1:1,-1:1, 0:3)
   real(kind=8) :: tanKx(-1:1,-1:1, 0:3)
   real(kind=8) :: tanKy(-1:1,-1:1, 0:3)
   real(kind=8) :: tanKz(-1:1,-1:1, 0:3)
end type t_xyzr

contains

!------------------------------------------------------------------------------------------------------------

subroutine rdsubs (unitnm, ic, ncase, linenr, idebug, subs)
!--purpose: read points in which the subsurface elastic field must be computed.
!           UnitNm gives the Unit-number of the input-file.
      implicit none
!--subroutine arguments:
      integer                  :: unitnm, linenr, ncase, idebug
      type(t_ic)               :: ic
      type(t_subsurf)          :: subs
!--local variables:
      integer, parameter :: mxnval = 20
      logical, parameter :: lstop = .true.
      integer            :: ints(mxnval), idum(1), isubs, nval, ieof, ierror
      logical            :: flags(mxnval), zerror, lchanged
      real(kind=8)       :: dbles(mxnval)
      character(len=256) :: strngs(mxnval)
      real(kind=8), dimension(:), pointer :: tmp

      ieof   = -1 ! eof=error

      ! first read the A- and O-digits regarding subsurface results

      if (ic%stress.ge.2) then
         call readLine(unitnm, ncase, linenr, 'A, O-digits for subsurface stress', 'ii',                &
                          ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)

         ic%matfil_subs = max(0, min(2, ints(1)))
         ic%output_subs = max(0, min(4, ints(2)))
      endif

      if (ic%stress.ge.3) then

         ! read the input option for the first block

         call readLine(unitnm, ncase, linenr, 'subsurface input option', 'i', ints, dbles,              &
                           flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
         isubs  = ints(1)
         zerror = .not.check_3rng ('ISUBS', isubs, 1,3, 5,7, 9,9)

         ! While (more blocks) do

         subs%nblock = 0
         do while (isubs.ge.1 .and. subs%nblock.lt.MXBLCK)

            ! increment number of blocks, set reference to block data

            subs%nblock = subs%nblock + 1
            associate(b => subs%blocks(subs%nblock))
            b%isubs = isubs

            ! write(bufout,*) 'rdsubs: starting block',subs%nblock,', isubs=',isubs
            ! call write_log(1, bufout)

            ! isubs = 1 - 7: read horizontal part

            if (isubs.eq.1 .or. isubs.eq.5) then

               b%ixl_inp = 1
               b%ixinc   = 1
               b%ixh_inp = 99999
               b%iyl_inp = 1
               b%iyinc   = 1
               b%iyh_inp = 99999

            elseif (isubs.eq.2 .or. isubs.eq.6) then

               call readLine(unitnm, ncase, linenr, 'subsurf: regular selection of grid columns',       &
                          'iii', ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)

               b%ixl_inp = ints(1)
               b%ixinc   = max(1, ints(2))
               b%ixh_inp = ints(3)
               zerror  = zerror .or. .not.check_irng ('IXL', b%ixl_inp,     1,     999999) 
               zerror  = zerror .or. .not.check_irng ('IXH', b%ixh_inp, b%ixl_inp, 999999)

               call readLine(unitnm, ncase, linenr, 'subsurf: regular selection of grid rows',          &
                          'iii', ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)


               b%iyl_inp = ints(1)
               b%iyinc   = max(1, ints(2))
               b%iyh_inp = ints(3)
               zerror  = zerror .or. .not.check_irng ('IYL', b%iyl_inp,     1,     999999) 
               zerror  = zerror .or. .not.check_irng ('IYH', b%iyh_inp, b%iyl_inp, 999999)

            elseif (isubs.eq.3 .or. isubs.eq.7) then

               ! get number of points in x- and y-directions

               call readLine(unitnm, ncase, linenr, 'subsurf: irregular selection of elements, nx/ny',  &
                          'ii', ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
               b%nx_inp = ints(1)
               b%ny_inp = ints(2)
               zerror   = zerror .or. .not.check_irng ('NX', b%nx_inp, 1, MXSUBS) 
               zerror   = zerror .or. .not.check_irng ('NY', b%ny_inp, 1, MXSUBS)
               b%nx_inp = min(b%nx_inp, MXSUBS)
               b%ny_inp = min(b%ny_inp, MXSUBS)

               ! read element numbers IX as real values

               allocate(tmp(max(1,b%nx_inp,b%ny_inp)))

               if (b%nx_inp.ge.1) then

                  call reallocate_arr(b%ixlist_inp, b%nx_inp)
                  call read1darr(unitnm, ncase, linenr, 'subsurf: ix-numbers of grid elements', 'd',    &
                              b%nx_inp, idum, tmp, idebug, lstop, ierror)

                  ! sort in ascending order, copy/maintain unique positive values

                  call unique_positive('IX', b%nx_inp, tmp, b%ixlist_inp, lchanged)
                  zerror  = zerror .or. .not.check_irng ('NX', b%nx_inp, 1, MXSUBS) 
               endif

               ! read element numbers IY as real values

               if (b%ny_inp.ge.1) then

                  call reallocate_arr(b%iylist_inp, b%ny_inp)
                  call read1darr(unitnm, ncase, linenr, 'subsurf: iy-numbers of grid elements', 'd',    &
                             b%ny_inp, idum, tmp, idebug, lstop, ierror)

                  !  sort in ascending order, copy/maintain unique positive values

                  call unique_positive('IY', b%ny_inp, tmp, b%iylist_inp, lchanged)
                  zerror  = zerror .or. .not.check_irng ('NY', b%ny_inp, 1, MXSUBS)
               endif

               deallocate(tmp)

            endif

            ! isubs = 1 - 7: read vertical part

            if (isubs.ge.1 .and. isubs.le.3) then

               ! regular grid in z-direction

               call readLine(unitnm, ncase, linenr, 'subsurf: regular spacing in z',                    &
                          'idd', ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
               zerror  = zerror .or. .not.check_irng ('NZ', ints(1), 1, MXSUBS)
               b%nz    = max(1, min(ints(1), MXSUBS))
               b%zl    = dbles(1)
               b%dz    = dbles(2)
               zerror  = zerror .or. .not.check_range('DZ', b%dz, 1d-9, 1d20)

            elseif (isubs.ge.5 .and. isubs.le.7) then

               ! explicit list of z-values, sorted in ascending order

               call readLine(unitnm, ncase, linenr, 'subsurf: explicit spacing in z',                   &
                          'i', ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
               zerror  = zerror .or. .not.check_irng ('NZ', ints(1), 1, MXSUBS)
               b%nz    = max(1, min(ints(1), MXSUBS))

               call reallocate_arr(b%z, b%nz)
               call read1darr(unitnm, ncase, linenr, 'subsurf: list of z-coordinates', 'd', b%nz,       &
                              idum, b%z, idebug, lstop, ierror)
               call bubsrt_real (b%nz, b%z)

            endif

            ! isubs = 9: explicit lists x, y, z

            if (isubs.eq.9) then

               ! read number of points nx, ny, nz

               call readLine(unitnm, ncase, linenr, 'number of subsurface points in x-, y-, z-directions', &
                          'iii', ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)

               zerror  = zerror .or. .not.check_irng ('NX', ints(1), 1, MXSUBS)
               zerror  = zerror .or. .not.check_irng ('NY', ints(2), 1, MXSUBS)
               zerror  = zerror .or. .not.check_irng ('NZ', ints(3), 1, MXSUBS)
               b%nx_inp = max(1, min(MXSUBS, ints(1)))
               b%ny_inp = max(1, min(MXSUBS, ints(2)))
               b%nz     = max(1, min(MXSUBS, ints(3)))

               ! write(bufout,*) 'rdsubs: nx=',b%nx_inp,', ny=',b%ny_inp,', nz=',b%nz
               ! call write_log(1, bufout)

               ! allocate arrays x, y, z at appropriate size

               call reallocate_arr(b%x, b%nx_inp)
               call reallocate_arr(b%y, b%ny_inp)
               call reallocate_arr(b%z, b%nz)

               ! read X-coordinates, sort in ascending order

               call read1darr(unitnm, ncase, linenr, 'x-coordinates of subsurface points', 'd',         &
                              b%nx_inp, idum, b%x, idebug, lstop, ierror)
               call bubsrt_real (b%nx_inp, b%x)

               ! read Y-coordinates, sort in ascending order

               call read1darr(unitnm, ncase, linenr, 'y-coordinates of subsurface points', 'd',         &
                              b%ny_inp, idum, b%y, idebug, lstop, ierror)
               call bubsrt_real (b%ny_inp, b%y)

               ! read Z-coordinates, sort in ascending order

               call read1darr(unitnm, ncase, linenr, 'z-coordinates of subsurface points', 'd',         &
                              b%nz, idum, b%z, idebug, lstop, ierror)
               call bubsrt_real (b%nz, b%z)

            endif

            ! Stop program if an error was found

            if (zerror) then
 8000          format(' Check block',i5,' of the subsurface input, at or before line',i8,'.',/,         &
                      ' Errors found in the input for case',i8,', aborting.')
               write(bufout, 8000) subs%nblock, linenr, ncase
               call write_log(2, bufout)
               call abort_run()
            endif

            ! read whether a new block will be entered

            call readLine(unitnm, ncase, linenr, 'subsurface input option', 'i', ints, dbles,              &
                           flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)
            isubs = ints(1)

            end associate

         ! end while (more blocks)

         enddo

      endif ! ic%stress.ge.3

end subroutine rdsubs

!------------------------------------------------------------------------------------------------------------

subroutine wrsubs (unitnm, ic, subs)
!--purpose: create input for points in which the subsurface elastic field must be calculated.
      implicit none
!--subroutine arguments
      integer                      :: unitnm
      type(t_ic)                   :: ic
      type(t_subsurf)              :: subs
!--local variables
      integer                      :: iblock, i

      if (ic%stress.ge.2) then

         write(unitnm,'(a)') '%  SUBSURFACE POINTS'
         write(unitnm,'(2i6,22x,a)') ic%matfil_subs, ic%output_subs, 'MATFIL, OUTPUT'

      endif

      if (ic%stress.ge.3) then

         do iblock = 1, subs%nblock

            associate(b => subs%blocks(iblock))

            ! write comment line for 2nd and later blocks

            write(unitnm, 101) iblock
 101        format('%  BLOCK',i3)

            ! write input option for next block

            write(unitnm, 102) b%isubs
 102        format(i6,28x,'ISUBS')

            ! isubs = 1 - 7: write horizontal part

            if (b%isubs.eq.1 .or. b%isubs.eq.5) then

               ! all elements of potential contact, no input

            elseif (b%isubs.eq.2 .or. b%isubs.eq.6) then

               ! regular selection of elements

               write(unitnm, 201) b%ixl_inp, b%ixinc, b%ixh_inp, b%iyl_inp, b%iyinc, b%iyh_inp
 201           format(3i6,16x,'IXL, IXINC, IXH',/, 3i6,16x,'IYL, IYINC, IYH')

            elseif (b%isubs.eq.3 .or. b%isubs.eq.7) then

               ! irregular selection of elements

               write(unitnm, 301) b%nx_inp, b%ny_inp
 301           format(2i6, 22x, 'NX, NY')

               write(unitnm, 302) 'IX', (b%ixlist_inp(i), i=1,b%nx_inp)
               write(unitnm, 302) 'IY', (b%iylist_inp(i), i=1,b%ny_inp)
 302           format('%  elements ',a,':',/,100(10i6,:,/))

            endif

            ! isubs = 1 - 7: write vertical part

            if (b%isubs.ge.1 .and. b%isubs.le.3) then

               ! regular grid in z-direction

               write(unitnm, 401) b%nz, b%zl, b%dz
 401           format(i6, f10.3, f10.5, 8x, 'NZ, ZL, DZ')

            elseif (b%isubs.ge.5 .and. b%isubs.le.7) then

               ! explicit list of z-values

               write(unitnm, 501) b%nz
 501           format(i6, 28x, 'NZ')

               write(unitnm, 502) (b%z(i), i=1,b%nz)
 502           format('%  z-coordinates:',/,100(10f10.4,:,/))

            endif

            ! isubs = 9: write explicit lists x, y, z

            if (b%isubs.eq.9) then

               write(unitnm, 601) b%nx_inp, b%ny_inp, b%nz
 601           format(3i6, 16x, 'NX, NY, NZ')

               write(unitnm, 602) 'x', (b%x(i), i=1,b%nx_inp)
               write(unitnm, 602) 'y', (b%y(i), i=1,b%ny_inp)
               write(unitnm, 602) 'z', (b%z(i), i=1,b%nz)
 602           format('%  ',a,'-coordinates:',/,100(10f10.4,:,/))

            endif

            end associate
         enddo ! blocks

         write(unitnm, 901) 0
 901     format('%  NO MORE BLOCKS',/,i6,28x,'ISUBS')

      endif ! ic%stress.ge.3

end subroutine wrsubs

!------------------------------------------------------------------------------------------------------------

subroutine subsur (meta, ic, mater, cgrid, igs, ps, mirror_y, subs)
!--purpose: compute the subsurface elastic field in the points given in  blocks(:)
   implicit none
!--subroutine parameters:
   type(t_metadata)             :: meta
   type(t_ic)                   :: ic
   type(t_material)             :: mater
   type(t_grid)                 :: cgrid
   type(t_eldiv)                :: igs
   type(t_gridfnc3)             :: ps
   logical,          intent(in) :: mirror_y
   type(t_subsurf)              :: subs
!--local variables:
   integer,        parameter :: idebug = 0

   ! Perform actual subsurface stress calculation, store results in subs%blocks

   call subsur_calc (ic, mater, cgrid, igs, ps, mirror_y, subs, idebug)

   ! Write mat-files per block, if requested

   call subsur_matfil (meta, ic, subs, idebug)

   ! Write data to out-file, if requested

   call subsur_outfile (ic, subs, idebug)

end subroutine subsur

!------------------------------------------------------------------------------------------------------------

subroutine expand_xb (cgrid, b, ncolum, idebug)
!--purpose: Expand the compressed storage for one block of subsurface points to full list of coordinates
   implicit none
!--subroutine arguments
   type(t_grid)                 :: cgrid
   type(t_subsblk)              :: b
   integer                      :: ncolum, idebug
!--local variables
   integer      :: ix, iy, iz, i, npoint

   if (idebug.ge.5) then
      write(bufout,*) 'expand_xb: isubs=',b%isubs
      call write_log(1, bufout)
   endif

   ! ISUBS = 9: copy inputs to effective values

   if (b%isubs.eq.9) then

      b%nx_eff = b%nx_inp
      b%ny_eff = b%ny_inp

   endif

   ! ISUBS = 1, 2, 5, 6: expand regular selection of discretization elements to list
   !                     clipping element numbers outside [1,mx] x [1,my]; assuming increments >= 1

   if (b%isubs.eq.1 .or. b%isubs.eq.2 .or. b%isubs.eq.5 .or. b%isubs.eq.6) then

      b%ixl_eff = max(b%ixl_inp, 1)
      b%ixh_eff = min(b%ixh_inp, cgrid%nx)
      b%iyl_eff = max(b%iyl_inp, 1)
      b%iyh_eff = min(b%iyh_inp, cgrid%ny)

      b%nx_eff  = (b%ixh_eff - b%ixl_eff + b%ixinc) / b%ixinc
      b%ny_eff  = (b%iyh_eff - b%iyl_eff + b%iyinc) / b%iyinc

      if (idebug.ge.2) then
         write(bufout,'(4(a,i3))') ' effective ixl=',b%ixl_eff,', inc=',b%ixinc,', ixh=',b%ixh_eff,     &
                ': nx=',b%nx_eff
         call write_log(1, bufout)
         write(bufout,'(4(a,i3))') ' effective iyl=',b%iyl_eff,', inc=',b%iyinc,', iyh=',b%iyh_eff,     &
                ': ny=',b%ny_eff
         call write_log(1, bufout)
      endif

      call reallocate_arr(b%ixlist_eff, b%nx_eff)
      call reallocate_arr(b%iylist_eff, b%ny_eff)

      i = 0
      do ix = b%ixl_eff, b%ixh_eff, b%ixinc
         i = i + 1
         b%ixlist_eff(i) = ix
         if (idebug.ge.3) then
            write(bufout,*) 'expand_xb: ixlist(',i,')=',b%ixlist_eff(i)
            call write_log(1, bufout)
         endif
      enddo

      i = 0
      do iy = b%iyl_eff, b%iyh_eff, b%iyinc
         i = i + 1
         b%iylist_eff(i) = iy
         if (idebug.ge.3) then
            write(bufout,*) 'expand_xb: iylist(',i,')=',b%iylist_eff(i)
            call write_log(1, bufout)
         endif
      enddo
   endif ! isubs 1,2,5,6

   ! ISUBS = 3, 7: remove element numbers outside [1,mx] x [1,my]
   !               assuming that ixlist_inp and iylist_inp are sorted in increasing order

   if (b%isubs.eq.3 .or. b%isubs.eq.7) then

      b%nx_eff = 0
      b%ny_eff = 0
      call reallocate_arr(b%ixlist_eff, b%nx_inp)
      call reallocate_arr(b%iylist_eff, b%ny_inp)

      do i = 1, b%nx_inp
         if (b%ixlist_inp(i).ge.1 .and. b%ixlist_inp(i).le.cgrid%nx) then
            b%nx_eff = b%nx_eff + 1
            b%ixlist_eff(b%nx_eff) = b%ixlist_inp(i)
         endif
      enddo

      do i = 1, b%ny_inp
         if (b%iylist_inp(i).ge.1 .and. b%iylist_inp(i).le.cgrid%ny) then
            b%ny_eff = b%ny_eff + 1
            b%iylist_eff(b%ny_eff) = b%iylist_inp(i)
         endif
      enddo

      b%ixl_eff = b%ixlist_eff(1)
      b%iyl_eff = b%iylist_eff(1)
      b%ixh_eff = b%ixlist_eff(b%nx_eff)
      b%iyh_eff = b%iylist_eff(b%ny_eff)

   endif ! isubs 3,7

   ! ISUBS = 1,2,3, 5,6,7: set X,Y coordinates for selected element centers

   if (b%isubs.ge.1 .and. b%isubs.le.7) then

      call reallocate_arr(b%x, b%nx_eff)
      call reallocate_arr(b%y, b%ny_eff)

      do i = 1, b%nx_eff
         b%x(i) = cgrid%x(1) + real(b%ixlist_eff(i)-1) * cgrid%dx
         if (idebug.ge.3) then
            write(bufout,*) 'expand_xb: x(',i,')=',b%x(i)
            call write_log(1, bufout)
         endif
      enddo

      do i = 1, b%ny_eff
         b%y(i) = cgrid%y(1) + real(b%iylist_eff(i)-1) * cgrid%dy
         if (idebug.ge.3) then
            write(bufout,*) 'expand_xb: y(',i,')=',b%y(i)
            call write_log(1, bufout)
         endif
      enddo
   endif

   ! expand regular z-grid to list of coordinates

   if (b%isubs.ge.1 .and. b%isubs.le.3) then

      call reallocate_arr(b%z, b%nz)

      do i = 1, b%nz
         b%z(i) = b%zl + real(i-1) * b%dz
         if (idebug.ge.3) then
            write(bufout,*) 'expand_xb: z(',i,')=',b%z(i)
            call write_log(1, bufout)
         endif
      enddo
   endif

   ! (re-)allocate table for storing the results before printing
   ! note: the table must be precisely nx*ny*nz long in the 1st array dimension

   b%ncolum = ncolum
   npoint   = b%nx_eff * b%ny_eff * b%nz

   if (idebug.ge.2) then
      write(bufout,'(3(a,i4),a,i8,a,i3)') ' expand_xb: nx,ny,nz=',b%nx_eff,',', b%ny_eff,',', b%nz,      &
                ', needs table:', npoint,' x', ncolum
      call write_log(1, bufout)
   endif
   call reallocate_arr(b%table, npoint, ncolum)

   ! expand lists of x, y, z-coordinates to full list (x,y,z)-combinations

   i = 0
   do iz = 1, b%nz
      do iy = 1, b%ny_eff
         do ix = 1, b%nx_eff
            i = i + 1
            b%table(i, 1) = b%x(ix)
            b%table(i, 2) = b%y(iy)
            b%table(i, 3) = b%z(iz)
            if (idebug.ge.4) then
               write(bufout,*) 'expand_xb: table(',i,')=',b%table(i,1),',',b%table(i,2),',',b%table(i,3)
               call write_log(1, bufout)
            endif
         enddo
      enddo
   enddo

end subroutine expand_xb

!------------------------------------------------------------------------------------------------------------

subroutine subs_mirror_tractions(cgrid, igs_in, ps_in, igs_mirr, ps_mirr)
!--purpose: reorder surface tractions and element divisions for mirrored grid
   implicit none
!--subroutine parameters:
   type(t_grid)     :: cgrid
   type(t_eldiv)    :: igs_in, igs_mirr
   type(t_gridfnc3) :: ps_in, ps_mirr
!--local variables:
   integer      :: mx, my, ix, iy_in, iy_out, ii_in, ii_out

   mx = cgrid%nx
   my = cgrid%ny

   do iy_in = 1, my
      iy_out = my + 1 - iy_in
      do ix = 1, mx
         ii_in  = ix + (iy_in -1) * mx
         ii_out = ix + (iy_out-1) * mx
         igs_mirr%el(ii_out) = igs_in%el(ii_in)
         ps_mirr%vn(ii_out)  =  ps_in%vn(ii_in)
         ps_mirr%vx(ii_out)  =  ps_in%vx(ii_in)
         ps_mirr%vy(ii_out)  = -ps_in%vy(ii_in)
      enddo
   enddo

end subroutine subs_mirror_tractions

!------------------------------------------------------------------------------------------------------------

subroutine subsur_calc (ic, mater, cgrid_arg, igs_arg, ps_arg, mirror_y, subs, idebug)
!--purpose: compute the subsurface elastic field in the points given in  blocks(:)
   implicit none
!--subroutine parameters:
   type(t_ic)                   :: ic
   type(t_material)             :: mater
   type(t_grid),     target     :: cgrid_arg
   type(t_eldiv),    target     :: igs_arg
   type(t_gridfnc3), target     :: ps_arg
   logical,          intent(in) :: mirror_y
   type(t_subsurf)              :: subs
   integer                      :: idebug
!--local variables:
   integer,         parameter :: ncolum = 21
   integer            :: ii, i, j, iz, kdb, npoint
   integer            :: iblock, iofs, mythrd, numthrd
   logical            :: usefft
   real(kind=8)       :: dbstep, sighyd, sigvm, sigtr, sigmaj(3), sigma(3,3), uw(3), xw(3)
   type(t_grid),     target  :: cgrid_mirr
   type(t_eldiv),    target  :: igs_mirr
   type(t_gridfnc3), target  :: ps_mirr
   type(t_grid),     pointer :: cgrid
   type(t_eldiv),    pointer :: igs
   type(t_gridfnc3), pointer :: ps

   call timer_start(itimer_subsur)

   if (idebug.ge.10) then
      open(unit=ltmp, file='subsurf.tmp', action='write')
      call wrsubs (ltmp, ic, subs)
      close(ltmp)
   endif

   if (ic%x_nmdbg.ge.4) then
      if (mater%gencr_eff.eq.3) then
         call write_log(' subsurface stresses: using bilinear tractions profile...')
      else
         call write_log(' subsurface stresses: using piecewise constant tractions...')
      endif
   endif

   ! mirror surface tractions and element division for a left side wheel/rail combination

   if (.not.mirror_y) then
      cgrid => cgrid_arg
      igs   => igs_arg
      ps    => ps_arg
   else
      call unifgrid_mirror_y(cgrid_arg, cgrid_mirr)
      call eldiv_new(igs_mirr, cgrid_mirr, nulify=.true.)
      call gf3_new(ps_mirr, 'ps_mirr', cgrid_mirr, nulify=.true.)

      cgrid => cgrid_mirr
      igs   => igs_mirr
      ps    => ps_mirr

      call subs_mirror_tractions(cgrid_arg, igs_arg, ps_arg, igs_mirr, ps_mirr)
   endif
   
   ! For all blocks of coordinates

   do iblock = 1, subs%nblock

      associate(b => subs%blocks(iblock))
      if (idebug.ge.3) then
         write(bufout,*) 'subsur_calc: starting block',iblock,', isubs=',b%isubs
         call write_log(1, bufout)
      endif

      ! Copy/expand data of block b, create/fill first three columns of b%table

      call expand_xb(cgrid, b, ncolum, idebug-3)
      npoint = b%nx_eff * b%ny_eff * b%nz

      if (idebug.ge.5) then
         write(bufout,*) 'subsur_calc: nx=',b%nx_eff,', ny=',b%ny_eff,', nz=',b%nz
         call write_log(1, bufout)
      endif

      ! Determine if use of FFTs is possible
      !  - currently supports only the full potential contact area, isubs=1 or 5.

      usefft = (b%isubs.le.7)

      if (usefft) then

         if (idebug.ge.3) call write_log(' STRESS: using fast subsurface stress calculation with FFTs')

         ! When using FFT's: loop over all z-values of the block

!_INT_OMP parallel do if (omp_in_parallel().eq.0)                                                        &
!_INT_OMP    default(none)                                                                               &
!_INT_OMP    shared(cgrid, mater, igs, ps, b, npoint )                                                   &
!_INT_OMP    private(iofs)

         do iz = 1, b%nz

            ! compute stresses at depth z(iz)

            iofs = (iz-1) * b%nx_eff * b%ny_eff
            call sstres_fft(cgrid, mater, igs, ps, b, iz, iofs)

         enddo
!_INT_OMP end parallel do
         if (idebug.ge.3) call write_log(' STRESS: done with parallel do')

      else

         ! Not using FFT's: compute stresses for all points of the block using a parallel loop

!_INT_OMP parallel if (omp_in_parallel().eq.0)                                                          &
!_INT_OMP    default(none)                                                                              &
!_INT_OMP    shared(cgrid, mater, ps, npoint, b, idebug, lscreen )                                      &
!_INT_OMP    private(ii, i, j, sighyd, sigvm, sigtr, sigmaj, sigma, uw, xw, dbstep, kdb, numthrd, mythrd )

!_INT_OMP do schedule(static,1)

         do 21 ii = 1, npoint

#ifdef _OPENMP
            numthrd = omp_get_num_threads()
            mythrd  = omp_get_thread_num()
#else
            numthrd = 1
            mythrd  = 0
#endif

            ! copy coordinates

            do j = 1, 3
               xw(j) = b%table(ii,j)
            enddo

            ! Progress output to screen: one line per block (idebug=1)
            !    or one line per point (idebug>=2)

            if (mythrd.eq.0) then
               if (idebug.eq.1) then

                  ! write dot for every 'dbstep'-th point of the block, in total
                  ! 50 dots per block, or less if the block has less points

                  if (ii.eq.1) then
                     dbstep = max(1d0, (npoint-1)/50.)
                     if (lscreen) write(*, '(a$)') ' Processing subsurf.points '
                     kdb    = 1
                  endif
                  do while (ii.ge.kdb*dbstep) 
                     if (lscreen) write(*, '(a$)') '.'
                     kdb = kdb + 1
                  end do
               elseif (idebug.ge.2) then
                  if (lscreen) write (*, 2032) ii, xw
 2032             format (' Processing point ',i4, ' : (',g12.4, ',', g12.4, ',', g12.4, ')')
               endif
            endif

            ! compute stresses

            call sstres (cgrid, mater, ps, uw, sighyd, sigvm, sigtr, sigmaj, sigma, xw)

            ! store values in a table for writing outside the parallel loop

            do i = 1, 3
               b%table(ii, 4: 6) = uw
               b%table(ii, 7)    = sighyd
               b%table(ii, 8)    = sigvm
               b%table(ii, 9)    = sigtr
               b%table(ii,10:12) = sigmaj
               b%table(ii,13:21) = reshape(sigma, (/ 9 /) )
            enddo

 21      continue
!_INT_OMP end do
!_INT_OMP end parallel

      endif ! usefft
      if (lscreen .and. idebug.eq.1) write(*,*)

      end associate
   enddo ! all blocks of coordinates

   if (mirror_y) then
      call grid_destroy(cgrid_mirr)
      call eldiv_destroy(igs_mirr)
      call gf3_destroy(ps_mirr)
   endif
   call timer_stop(itimer_subsur)

end subroutine subsur_calc

!------------------------------------------------------------------------------------------------------------

subroutine subsur_matfil (meta, ic, subs, idebug)
!--purpose: write the subsurface elastic field to a .subs-file for importing in Matlab
   implicit none
!--subroutine parameters:
   type(t_metadata)         :: meta
   type(t_ic)               :: ic
   type(t_subsurf)          :: subs
   integer                  :: idebug
!--local variables:
   integer            :: lunmat, ncase, iblock, ii, i, j, k, npoint
   logical            :: lwrmat, full_tensor
   real(kind=8)       :: sighyd, sigvm, sigtr, sigmaj(3), sigma(3,3), uw(3), xw(3)
   character(len=256) :: fname

   ! Write to subs-file "<experim>.<case>.subs" when A=1 or 2

   lwrmat = (ic%matfil_subs.ge.1)

   ! Write full stress tensor to subs-file when A=3

   full_tensor = (ic%matfil_subs.ge.2)

   ! Open output-file for communication to Matlab
   ! using names <EXPERIM>.0001.subs, for cases >= 1M, use name  <EXPERIM>.subs.

   ncase = meta%ncase
   if (ncase.le.9999) then
      fname = trim(meta%expnam) // '.'                                                            &
                     // char(ichar('0')+    ncase/1000)                                         &
                     // char(ichar('0')+mod(ncase,1000)/100)                                    &
                     // char(ichar('0')+mod(ncase, 100)/ 10)                                    &
                     // char(ichar('0')+mod(ncase,  10)/  1)
   elseif (ncase.le.999999) then
      fname = trim(meta%expnam) // '.'                                                            &
                     // char(ichar('0')+    ncase/100000)                                       &
                     // char(ichar('0')+mod(ncase,100000)/10000)                                &
                     // char(ichar('0')+mod(ncase, 10000)/ 1000)                                &
                     // char(ichar('0')+mod(ncase,  1000)/  100)                                &
                     // char(ichar('0')+mod(ncase,   100)/   10)                                &
                     // char(ichar('0')+mod(ncase,    10)/    1)
   else
      fname = trim(meta%expnam)
   endif

   if (meta%ipatch.ge.1) then
      fname = trim(fname) // char(ichar('a')+meta%ipatch-1) // '.subs'
   else
      fname = trim(fname) // '.subs'
   endif

   if (meta%dirnam.ne.' ') then
      fname = trim(meta%dirnam) // path_sep // trim(fname)
   endif

   lunmat = lsbs
   if (lunmat.le.0) then
      write(*,*) 'ERROR: no unit-number provided for .subs-file; continuing without write.'
      lwrmat = .false.
   endif

   if (lwrmat) then
      open (unit=lunmat, file=fname, action='write', err=998)
      goto 999

      ! Error handling:

 998  continue
         write(*,*) 'ERROR: cannot write stresses to .subs-file; continuing without write.'
         lwrmat = .false.
         goto 999
 999  continue
   endif

   ! For all blocks of coordinates

   call timer_start(itimer_subsfile)
   if (lwrmat) then
      do iblock = 1, subs%nblock

         associate(b => subs%blocks(iblock))
         npoint =  b%nx_eff * b%ny_eff * b%nz

         if (idebug.ge.3) then
            write(bufout,*) 'subsur_matfil: starting block',iblock,', isubs=',b%isubs
            call write_log(1, bufout)
         endif

         ! write block-header to .subs-file

         if (.not.full_tensor) then
            write (lunmat,2121)
            write (lunmat,2122) b%nx_eff, b%ny_eff, b%nz, (0.0, k=1,5)
            write (lunmat,2123)
         else
            write (lunmat,2126)
            write (lunmat,2127) b%nx_eff, b%ny_eff, b%nz, (0.0, k=1,11)
            write (lunmat,2128)
         endif

 2121    format ('%   NX',10x,'NY',10x,'NZ',5(11x,'-'))
 2122    format (i6,2i12,6x,5g12.4)
 2123    format ('%    X',11x, 'Y',11x, 'Z',10x,'UX',10x,'UY',10x,'UZ',8x,'SIGHYD',7x,'SIGVM')

 2126    format ('%   NX',10x,'NY',10x,'NZ',11(11x,'-'))
 2127    format (i6,2i12,6x,11g12.4)
 2128    format ('%    X',11x, 'Y',11x, 'Z',10x,'UX',10x,'UY',10x,'UZ',8x,'SIGHYD',6x,'SIGVM',7x,       &
                 'SIGXX',7x,'SIGXY',7x,'SIGXZ',7x,'SIGYY',7x,'SIGYZ',7x,'SIGZZ')

         ! write data to .subs-file

         do ii = 1, npoint

            xw     = b%table(ii, 1: 3)
            uw     = b%table(ii, 4: 6)
            sighyd = b%table(ii, 7)
            sigvm  = b%table(ii, 8)
            sigtr  = b%table(ii, 9)
            sigmaj = b%table(ii,10:12)
            sigma  = reshape(b%table(ii,13:21), (/ 3, 3 /) )

            if (.not.full_tensor) then
               write (lunmat, 2251) (fmt_gs(12,6,5,xw(i)), i=1,3), (fmt_gs(12,6,5,uw(i)), i=1,3),       &
                                    fmt_gs(12,6,5,sighyd), fmt_gs(12,6,5,sigvm)
            else
               write (lunmat, 2252) (fmt_gs(12,6,5,xw(i)), i=1,3), (fmt_gs(12,6,5,uw(i)), i=1,3),       &
                                    fmt_gs(12,6,5,sighyd), fmt_gs(12,6,5,sigvm),                        &
                                    (fmt_gs(12,6,5,sigma(1,j)), j=1,3),                                 &
                                    (fmt_gs(12,6,5,sigma(2,j)), j=2,3), fmt_gs(12,6,5,sigma(3,3))
            endif

 2251       format ( 8a)
 2252       format (14a)

         enddo

         end associate
      enddo ! all blocks of coordinates

      close(lunmat)
   endif ! lwrmat
   call timer_stop(itimer_subsfile)

end subroutine subsur_matfil

!------------------------------------------------------------------------------------------------------------

subroutine subsur_outfile (ic, subs, idebug)
!--purpose: write subsurface field to out-file
   implicit none
!--subroutine parameters:
   type(t_ic)                 :: ic
   type(t_subsurf)            :: subs
   integer                    :: idebug
!--local variables:
   integer            :: iblock, ii, i, j, ii_max, npoint
   real(kind=8)       :: sighyd, sigvm, sigtr, sigmaj(3), sigma(3,3), uw(3), xw(3)

   call timer_start(itimer_subsfile)

   ! For all blocks of coordinates

   do iblock = 1, subs%nblock

      associate(b => subs%blocks(iblock))
      npoint =  b%nx_eff * b%ny_eff * b%nz

      if (idebug.ge.3) then
         write(bufout,*) 'subsur_outfile: starting block',iblock,', isubs=',b%isubs
         call write_log(1, bufout)
      endif

      ! write maximum values to .out-file when O>=1

      if (ic%output_subs.ge.1 .and. out_open.eq.1) then

         write(lout, 2101) iblock
 2101    format ( /, ' SUBSURFACE STRESSES BLOCK',i3,':' )

         write(lout, *)

         ii_max = idamax(npoint, b%table(1:, 7), 1)
         write(lout, 2061) 'ABSMAX SIGHYD', b%table(ii_max, 7), b%table(ii_max,1), b%table(ii_max,2),   &
                b%table(ii_max,3)

         ii_max = idamax(npoint, b%table(1:, 8), 1)
         write(lout, 2061) 'MAX     SIGVM', b%table(ii_max, 8), b%table(ii_max,1), b%table(ii_max,2),   &
                b%table(ii_max,3)

         ii_max = idamax(npoint, b%table(1:, 9), 1)
         write(lout, 2061) 'MAX     SIGTR', b%table(ii_max, 9), b%table(ii_max,1), b%table(ii_max,2),   &
             b%table(ii_max,3)

         ! SIGMA1: compute maximum value instead of in absolute sense maximum

         ii_max = idmax(npoint, b%table(1:,10), 1)
         write(lout, 2061) 'MAX    SIGMA1', b%table(ii_max,10), b%table(ii_max,1), b%table(ii_max,2),   &
             b%table(ii_max,3)

         ! SIGMA3: compute minimum value instead of in absolute sense maximum

         ii_max = idmin(npoint, b%table(1:,12), 1)
         write(lout, 2061) 'MIN    SIGMA3', b%table(ii_max,12), b%table(ii_max,1), b%table(ii_max,2),   &
             b%table(ii_max,3)

         ii_max = idamax(npoint, b%table(1:,13), 1)
         write(lout, 2061) 'ABSMAX  SIGXX', b%table(ii_max,13), b%table(ii_max,1), b%table(ii_max,2),   &
             b%table(ii_max,3)

         ii_max = idamax(npoint, b%table(1:,17), 1)
         write(lout, 2061) 'ABSMAX  SIGYY', b%table(ii_max,17), b%table(ii_max,1), b%table(ii_max,2),   &
             b%table(ii_max,3)

         ii_max = idamax(npoint, b%table(1:,21), 1)
         write(lout, 2061) 'ABSMAX  SIGZZ', b%table(ii_max,21), b%table(ii_max,1), b%table(ii_max,2),   &
             b%table(ii_max,3)

 2061    format(1x,a13,' =',f10.3,' AT (X,Y,Z) = (',f8.3,',',f8.3,',',f8.3,')')
      endif

      ! write point-wise output to file "<experim>.out" when O=4

      if (ic%output_subs.ge.4 .and. out_open.eq.1) then

         write(lout,*)

         do ii = 1, npoint

            ! get values from the table for writing outside the parallel loop

            xw     = b%table(ii, 1: 3)
            uw     = b%table(ii, 4: 6)
            sighyd = b%table(ii, 7)
            sigvm  = b%table(ii, 8)
            sigtr  = b%table(ii, 9)
            sigmaj = b%table(ii,10:12)
            sigma  = reshape(b%table(ii,13:21), (/ 3, 3 /) )

            ! write data to out-file

            write(lout, 2201) xw
            write(lout, 2202) (i, uw(i), (sigma(i,j), j=1,3), i=1,3)
            write(lout, 2203) sighyd, sigvm

 2201       format (' Position (X,Y,Z):',3g12.4)
 2202       format (' I,Displacement(I); Stress(I,J),J = 1,3',                                  &
                     3( / ,i3,2x,g12.4,'; ',3g12.4))
 2203       format ('  SigHyd = ',g12.4,', SigVM = ',g12.4, / )
         enddo
      endif

   ! end do (all blocks of coordinates)

      end associate
   enddo
   call timer_stop(itimer_subsfile)

end subroutine subsur_outfile

!------------------------------------------------------------------------------------------------------------

subroutine sstres_fft(cgrid, mater, igs, ps, b, iz, iofs)
!--purpose: Calculates the internal displacement, displacement gradients, strains and stresses
!           together with the first invariant of the stress and the second invariant of the
!           stress deviator at the centers of all elements at depth ZW due to the traction distribution ps.
! uses subroutine stres1, that computes the displacement and displacement gradients in XW due
!           to a unit load in one element of the contact area.
   implicit none
!--subroutine arguments:
   integer,          intent(in)    :: iz, iofs
   type(t_grid),     intent(in)    :: cgrid
   type(t_material), intent(inout) :: mater
   type(t_eldiv),    intent(in)    :: igs
   type(t_gridfnc3), intent(in)    :: ps
   type(t_subsblk)                 :: b
!--local variables :
   integer          :: ia, ix, iy, itb, ii, ii_x, ii_y, i, j, k, ik, neg
   real(kind=8)     :: sigma(3,3), sighyd, sigvm, sigtr, sigmaj(3), xw(3), uw(3), vr_ii(3,0:3)
   type(t_gridfnc3) :: tmp_vr
   type(t_inflcf)   :: ck(0:3)
   real(kind=8), dimension(:,:,:), allocatable :: vr

   ! This subroutine assumes dimensioned tractions ps [N/mm2] and position xw [mm], and
   ! produces dimensioned displacements vr(:,0) [mm], displacement gradients vr(:,1:3) [mm/mm],
   ! stress components sigma [N/mm2] and derived stresses (invariants) sighyd, sigvm, sigma1-3 [N/mm2].

   associate(table => b%table, zw => b%z(iz))

   ! Create data-structures for storing influence coefficients and other data

   allocate(vr(cgrid%ntot,3,0:3))
   call gf3_new(tmp_vr, 'tmp_vr', cgrid, igs)
   do k = 0, 3
      call inflcf_new(ck(k), 0, cgrid)
      call inflcf_mater(ck(k), mater)   ! copy material parameters from mater to ck(k)
      ck(k)%use_3bl = .false.           ! exclude 3rd body flexibility
   enddo

   ! The directions are numbered x=1, y=2, z=3.

   ! Stres1 computes v(i,j,k): the deformed distance in direction j due to a unit load (1 N/mm2) in
   ! direction i in one element ii. Index k=0 gives the deformed distance itself in [mm]. k=1-3 denote
   ! differentiation w.r.t. the three coordinate directions ([mm/mm]).
   ! This subroutine computes the aggregate deformation (vr(:,0), [mm]) and deformation gradients
   ! (vr(:,1:3), [mm/mm]) due to the contact loading in all elements jj=1,npot in all directions jk=1,3.

   ! Note that subroutine stres1 treats the positive half space only, i.e. it requires XW(3) > 0.
   ! Therefore in case of XW(3)<0 we introduce a second right handed coordinate system O(x',y',z') by
   ! x'==x, y'==-y and z'==-z.
   ! This change of variables affects XW and XP, as well as the meaning of deformations in directions
   ! i=1,2 and derivatives with respect to directions j=1,2.
   ! Further note that PS gives the traction on body 1, such that the traction on body 2 in the new
   ! coordinate system becomes
   !    Px(on 2, sys O') == -Px(on 1, sys O') == -PS(:,1) (on 1, sys O)
   !    Py(on 2, sys O') == -Py(on 1, sys O') ==  PS(:,2) (on 1, sys O)
   !    Pz(on 2, sys O') == -Pz(on 1, sys O') ==  PS(:,3) (on 1, sys O)

   ! set optional negative sign for mirroring in Z-direction (body 2, ZW<0)

   if (zw.ge.0) then
      ia  =  1
      neg =  1
   else
      ia  =  2
      neg = -1
   endif

   ! Compute influence coefficients for depth ZW for all offsets relevant to the potential contact area

   call sstres_inflcf(cgrid, mater, b%ixl_eff, b%ixh_eff, b%iyl_eff, b%iyh_eff, zw, ck)

   ! Z<0: transform to equivalent problem with Z>0 by negating px

   call timer_start(itimer_subsfft)
   if (ia.eq.2) call gf3_scal(AllElm, -1d0, ps, ikXDIR)

   ! Initialize displacements and displacement gradients for all elements

   do k = 0, 3
      do ik = 1, 3
         do ii = 1, cgrid%ntot
            vr (ii,ik,k) = 0d0
         enddo
      enddo
   enddo

   ! Compute subsurface stresses at depth ZW at all points of the potential contact area

   do k = 0, 3

      if (b%isubs.eq.1 .or. b%isubs.eq.5) then

         call VecAijPj (igs, AllElm, tmp_vr, ikALL, ps, jkALL, ck(k))
         do ik = 1, 3
            do ii = 1, cgrid%ntot
               vr(ii,ik,k) = tmp_vr%val(ii,ik)
            enddo
         enddo

      elseif (b%isubs.eq.2 .or. b%isubs.eq.3 .or. b%isubs.eq.6 .or. b%isubs.eq.7) then

         do j = 1, b%ny_eff
            iy = b%iylist_eff(j)
            do i = 1, b%nx_eff
               ix = b%ixlist_eff(i)
               ii = ix + (iy-1) * cgrid%nx
               do ik = 1, 3
                  vr(ii,ik,k) = AijPj(ii, ik, ps, jkALL, ck(k))
               enddo
            enddo
         enddo

      else

         call write_log('ERROR: sstres_fft: isubs not supported.')
         call abort_run()

      endif
   enddo
   call timer_stop(itimer_subsfft)

   ! Compute derived quantities: linearised strain, dilatation, Hookean stresses
   
   call timer_start(itimer_subsderiv)
   do ii_y = 1, b%ny_eff
      iy = b%iylist_eff(ii_y)
      do ii_x = 1, b%nx_eff
         ix  = b%ixlist_eff(ii_x)
         ii  = ix   + (iy   - 1) * cgrid%nx
         itb = ii_x + (ii_y - 1) * b%nx_eff 

         ! Note: final coordinates instead of negated ones

         xw(1) = cgrid%x(ii)
         xw(2) = cgrid%y(ii)
         xw(3) = zw

         vr_ii(1:3,0:3) = vr(ii,1:3,0:3)
         call sstres_derived(mater%gg(ia), mater%poiss(ia), neg, vr_ii, uw, sigma, sigmaj, sighyd,      &
                             sigvm, sigtr)

         ! store values in a table for writing outside sstres_fft

         table(iofs+itb, 4: 6) = uw(1:3)
         table(iofs+itb, 7)    = sighyd
         table(iofs+itb, 8)    = sigvm
         table(iofs+itb, 9)    = sigtr
         table(iofs+itb,10:12) = sigmaj(1:3)
         table(iofs+itb,13:21) = reshape(sigma, (/ 9 /) )

      enddo ! ix = ixlist
   enddo ! iy = iylist
   call timer_stop(itimer_subsderiv)

   deallocate(vr)
   call gf3_destroy(tmp_vr)
   do k = 0, 3
      call inflcf_destroy(ck(k))
   enddo

   end associate
end subroutine sstres_fft

!------------------------------------------------------------------------------------------------------------

subroutine sstres_inflcf(cgrid, mater, ix0, ix1, iy0, iy1, zw, ck)
!--purpose: Calculates the internal displacement, displacement gradients, strains and stresses
!           together with the first invariant of the stress and the second invariant of the
!           stress deviator at the centers of all elements at depth ZW due to the traction distribution ps.
! uses subroutine stres1, that computes the displacement and displacement gradients in XW due
!           to a unit load in one element of the contact area.
   implicit none
!--subroutine arguments:
   type(t_grid),     intent(in)    :: cgrid
   type(t_material), intent(inout) :: mater
   integer,          intent(in)    :: ix0, ix1, iy0, iy1
   real(kind=8),     intent(in)    :: zw
   type(t_inflcf)                  :: ck(0:3)
!--local variables :
   integer          :: ia, ix, iy, nx, ny, i, j, k, ik, jk, ifaci, ifacj, ifack, neg
   real(kind=8)     :: xp(2), xw(3), v(3,3,0:3), vnu(3,3,0:3)

   ! The directions are numbered x=1, y=2, z=3.

   ! Stres1 computes v(i,j,k): the deformed distance in direction j due to a unit load (1 N/mm2) in
   ! direction i in one element ii. Index k=0 gives the deformed distance itself in [mm]. k=1-3 denote
   ! differentiation w.r.t. the three coordinate directions ([mm/mm]).
   ! This subroutine computes the aggregate deformation (vr(:,0), [mm]) and deformation gradients
   ! (vr(:,1:3), [mm/mm]) due to the contact loading in all elements jj=1,npot in all directions jk=1,3.

   ! Note that subroutine stres1 treats the positive half space only, i.e. it requires XW(3) > 0.
   ! Therefore in case of XW(3)<0 we introduce a second right handed coordinate system O(x',y',z') by
   ! x'==x, y'==-y and z'==-z.
   ! This change of variables affects XW and XP, as well as the meaning of deformations in directions
   ! i=1,2 and derivatives with respect to directions j=1,2.
   ! Further note that PS gives the traction on body 1, such that the traction on body 2 in the new
   ! coordinate system becomes
   !    Px(on 2, sys O') == -Px(on 1, sys O') == -PS(:,1) (on 1, sys O)
   !    Py(on 2, sys O') == -Py(on 1, sys O') ==  PS(:,2) (on 1, sys O)
   !    Pz(on 2, sys O') == -Pz(on 1, sys O') ==  PS(:,3) (on 1, sys O)

   call timer_start(itimer_subsinfl)

   ! set optional negative sign for mirroring in Z-direction (body 2, ZW<0)

   if (zw.ge.0d0) then
      ia  =  1
      neg =  1
   else
      ia  =  2
      neg = -1
   endif

   ! Compute influence coefficients for depth ZW for all offsets relevant to the potential contact area

   nx = max( ix1-1, cgrid%nx-ix0 ) + 1
   ny = max( iy1-1, cgrid%ny-iy0 ) + 1

   do iy = -ny, 0
      do ix = -nx, 0

         ! put unit traction at an element centered at the origin

         xp(1) = 0d0
         xp(2) = 0d0

         ! observe displacements and displacement gradients at offset (ix*dx, iy*dy)

         xw(1) =       real(ix)*cgrid%dx
         xw(2) = neg * real(iy)*cgrid%dy
         xw(3) = neg * zw

         ! gencr_eff = 2: constant basisfunctions, 3: bilinear basisfunctions

         if (mater%gencr_eff.eq.3) then
            call stres1_bilin(cgrid%dx, cgrid%dy, mater%gg(ia), v, vnu, xw, xp)
         else
            call stres1_pcwcns(cgrid%dx, cgrid%dy, mater%gg(ia), v, vnu, xw, xp)
         endif

         ! store result in arrays of influence coefficients
         ! note that ck* will be used in gf3_AijPj and consequently must use CONTACT-indices.
         ! note that 3rd index in cf is the direction of displacement, i.e.  index j in v, vnu
         !           4th index in cf is the direction of traction, i.e. index i in v, vnu

         do k = 0, 3
            do i = 1, 3
               do j = 1, 3
                  ck(k)%cf(ix, iy, j, i) = v(i,j,k) + mater%poiss(ia) * vnu(i,j,k)
               enddo
            enddo
         enddo
      enddo
   enddo

   ! Copy values for negative ix to positive ix using appropriate minus signs

   do k = 0, 3
      ifack = 1
      if (k.eq.1) ifack = -1
      do ik = 1, 3
         ifaci = 1
         if (ik.eq.1) ifaci = -1
         do jk = 1, 3
            ifacj = 1
            if (jk.eq.1) ifacj = -1
            do iy = -ny, 0
               do ix = 1, nx-1
                  ck(k)%cf(ix, iy, jk, ik) = ifaci * ifacj * ifack * ck(k)%cf(-ix, iy, jk, ik)
               enddo
            enddo
         enddo
      enddo
   enddo

   ! Copy values for negative iy to positive iy using appropriate minus signs

   do k = 0, 3
      ifack = 1
      if (k.eq.2) ifack = -1
      do ik = 1, 3
         ifaci = 1
         if (ik.eq.2) ifaci = -1
         do jk = 1, 3
            ifacj = 1
            if (jk.eq.2) ifacj = -1
            do iy = 1, ny-1
               do ix = -nx, nx-1
                  ck(k)%cf(ix, iy, jk, ik) = ifaci * ifacj * ifack * ck(k)%cf(ix, -iy, jk, ik)
               enddo
            enddo
         enddo
      enddo
   enddo

   ! override scaling by mater%ga, stres1 delivered the full influence coefficients including ga

   do k = 0, 3
      ck(k)%ga     = 1d0
      ck(k)%ga_inv = 1d0
   enddo

   ! print influence coefficients to output-file when requested

   if (.false. .and. out_open.eq.1) then
      do k = 0, 3
         write(lout,*) 'Printing influence coefficients ck',k
         call inflcf_print (ck(k), lout)
      enddo
   endif
   call timer_stop(itimer_subsinfl)

end subroutine sstres_inflcf

!------------------------------------------------------------------------------------------------------------

subroutine sstres(cgrid, mater, ps, uw, sighyd, sigvm, sigtr, sigmaj, sigma, xw)
!--purpose: Calculates the internal displacements and stresses in the point "XW" due to the traction
!           distribution ps.
! uses subroutine stres1, that computes the displacement and displacement gradients in XW due
!           to a unit load in one element of the contact area.
      implicit none
!--subroutine parameters:
      type(t_grid)     :: cgrid
      type(t_material) :: mater
      real(kind=8)     :: xw(3)
      type(t_gridfnc3) :: ps
!--output parameters:
      real(kind=8)     :: uw(3), sighyd, sigvm, sigtr, sigmaj(3), sigma(3,3)

   ! This subroutine assumes dimensioned tractions ps [N/mm2] and position xw [mm], and produces
   ! dimensioned displacements uw(:) [mm], stress components sigma [N/mm2] and derived stresses
   ! sighyd, sigvm, sigtr, sigma1-3 [N/mm2].

!--local variables:
      integer      :: ia, ii, i, j, k, neg, ifaci
      real(kind=8) :: v(3,3,0:3), vr(3,0:3), vnu(3,3,0:3), xp(2)

   ! In subroutine stres1 the directions are numbered x=1, y=2, z=3. This subroutine largely
   ! complies with this convention. Stres1 directions 1,2,3 correspond to CONTACT directions 1,2,3

   ! Stres1 computes v(i,j,k): the deformed distance in direction j due to a unit load (1 N/mm2) in
   ! direction i in one element ii. Index k=0 gives the deformed distance itself in [mm]. k=1-3 denote
   ! differentiation w.r.t. the three coordinate directions ([mm/mm]).
   ! This subroutine computes the aggregate deformation (vr(:,0), [mm]) and deformation gradients
   ! (vr(:,1:3), [mm/mm]) due to the contact loading in all elements jj=1,npot in all directions jk=1,3.

      do j = 1, 3
         do k = 0, 3
            vr(j,k) = 0d0
         enddo
      enddo

   ! Note that subroutine stres1 treats the positive half space only, i.e. it requires XW(3) > 0.
   ! Therefore in case of XW(3)<0 we introduce a second right handed coordinate system O(x',y',z') by
   ! x'==x, y'==-y and z'==-z.
   ! This change of variables affects XW and XP, as well as the meaning of deformations in directions
   ! ik=1,2 and derivatives with respect to directions jk=1,2.
   ! Further note that PS gives the traction on body 1, such that the traction on body 2 in the new
   ! coordinate system becomes
   !    Px(on 2, sys O') == -Px(on 1, sys O') == -PS(:,1) (on 1, sys O)
   !    Py(on 2, sys O') == -Py(on 1, sys O') ==  PS(:,2) (on 1, sys O)
   !    Pz(on 2, sys O') == -Pz(on 1, sys O') ==  PS(:,3) (on 1, sys O)

      if (xw(3).ge.0) then
         ia  =  1
         neg =  1
      else
         ia  =  2
         neg = -1
      endif

      xw(3) = neg * xw(3)
      xw(2) = neg * xw(2)

      do ii = 1, cgrid%ntot

         ! Compute the influence of a unit load in element ii with position XP on the deformed distance in XW

         if (ps%vn(ii).gt.0d0) then

            xp(1) =       cgrid%x(ii)
            xp(2) = neg * cgrid%y(ii)

            ! gencr_eff = 2: constant basisfunctions, 3: bilinear basisfunctions

            if (mater%gencr_eff.eq.3) then
               call stres1_bilin(cgrid%dx, cgrid%dy, mater%gg(ia), v, vnu,  xw, xp)
            else
               call stres1_pcwcns(cgrid%dx, cgrid%dy, mater%gg(ia), v, vnu,  xw, xp)
            endif

            ! v == direct displacements & gradients; vnu == terms proportional to Poiss.ratio
            !      first index i: direction of unit load

            ! add contributions of ps(ii,1:3) to vr(j,k)

            do j = 1, 3
               do k = 0, 3
                  do i = 1, 3
                     ifaci = 1
                     if (i.eq.1) ifaci = neg
                     vr(j,k) = vr(j,k) + ifaci*ps%val(ii,i) * (v(i,j,k) + mater%poiss(ia)*vnu(i,j,k))
                  enddo
               enddo
            enddo

         endif ! pn>0
      enddo ! ii

      ! restore old value of XW(2:3)

      xw(2) = neg * xw(2)
      xw(3) = neg * xw(3)

      ! compute derived quantities sigma, sigmaj, sighyd, sigvm, etc.

      call sstres_derived(mater%gg(ia), mater%poiss(ia), neg, vr, uw, sigma, sigmaj, sighyd, sigvm, sigtr)

end subroutine sstres

!------------------------------------------------------------------------------------------------------------

subroutine sstres_derived(gg, poiss, neg, vr, uw, sigma, sigmaj, sighyd, sigvm, sigtr)
!--purpose: from displacements and gradients vr, compute stresses and other derived quantities
      implicit none
!--subroutine parameters:
      integer,      intent(in)    :: neg
      real(kind=8), intent(in)    :: gg, poiss
      real(kind=8), intent(inout) :: vr(3,0:3)
      real(kind=8), intent(out)   :: uw(3), sigma(3,3), sigmaj(3), sighyd, sigvm, sigtr
!--local variables:
      real(kind=8), parameter :: tolsml = 1d-15
      integer      :: i, j, k, ifacj, ifack, info
      real(kind=8) :: dil, er(3,3), sigii, sijsij, sigwrk(3,3), work(9)

      ! apply negation of deformations j=2,3 and derivatives k=2,3

      do j = 1, 3
         ifacj = 1
         if (j.eq.2 .or. j.eq.3) ifacj = neg
         do k = 0, 3
            ifack = 1
            if (k.eq.2 .or. k.eq.3) ifack = neg
            vr(j,k) = ifacj * ifack * vr(j,k)
         enddo
      enddo

      ! Copy displacements to output variable uw

      uw(1:3) = vr(1:3,0)

      ! Compute derived quantities: linearised strain, dilatation, Hookean stresses
      ! (cf. eq.(3.1) p.100, note G=E / 2(1 + nu))

      do i = 1, 3
         do j = 1, 3
            er(i,j) = (vr(i,j) + vr(j,i)) / 2d0
         enddo
      enddo

      dil = er(1,1) + er(2,2) + er(3,3)

      do j = 1, 3
         do i = 1, 3
            sigma(i,j) = 2d0 * gg * er(i,j)
         enddo
      enddo

      ! Note that if poiss = nu = 0.5, then dil should be 0

      do i = 1, 3
         sigma(i,i) = sigma(i,i) + 2d0 * gg * dil * poiss / max(1d-6, 1d0 - 2d0*poiss)
      enddo

      ! Compute aggregate quantities: first stress invariant sigii, hydrostatic stress sigii/3

      sigii  = sigma(1,1) + sigma(2,2) + sigma(3,3)
      sighyd = sigii / 3d0

      ! Compute second invariant sijsij

      sijsij = - (1d0 / 3d0) * sigii**2
      do i = 1, 3
         do j = 1, 3
            sijsij = sijsij + sigma(i,j)**2
         enddo
      enddo
      sijsij = 0.5d0 * sijsij

      ! Compute the stress sigvm used in the von Mises yield criterion

      sigvm  = sqrt(3d0 * sijsij)

      ! Compute principal stresses: eigenvalues of the stress tensor

      if (.true.) then
         call dsyevc3(sigma, sigmaj)
      else
         sigwrk(1:3,1:3) = sigma(1:3,1:3)
         call dsyev('N', 'L', 3, sigwrk, 3, sigmaj, work, 9, info)

         if (info.ne.0) then
            write(bufout,'(a,i4)') ' Something wrong in dsyev, info=',info
            call write_log(1, bufout)
         endif

         ! reorder eigenvalues largest to smallest

         work(1)   = sigmaj(3)
         sigmaj(3) = sigmaj(1)
         sigmaj(1) = work(1)
      endif

      ! write(bufout,'(a,3f12.6)') 'sigmaj=',sigmaj(1), sigmaj(2), sigmaj(3)
      ! call write_log(1,bufout)

      ! Compute Tresca stress: largest - smallest principal stresses

      sigtr  = sigmaj(1) - sigmaj(3)

      ! filter small values in uw and sigma

      do i = 1, 3
         if (abs(uw(i)).lt.gg*tolsml) uw(i) = 0d0
      enddo

      do i = 1, 3
         do j = 1, 3
            if (abs(sigma(i,j)).lt.gg*tolsml) sigma(i,j) = 0d0
         enddo
      enddo

end subroutine sstres_derived

!------------------------------------------------------------------------------------------------------------

subroutine stres1_pcwcns(dx, dy, gg, v, vnu, xw, xp)
!--purpose: This subroutine stres1, formerly called nstres, calculates the displacements and displacement
!           gradients in an elastic half-space due to a surface load uniformly distributed over an
!           arbitrary rectangle with sides parallel to the x,y-axes (directions 1,2) of a Cartesian
!           coordinate system, in which the half-space is given by "xw(3).ge.0".
!
!           See J.J. Kalker, "Numerical calculation of the elastic field in a half-space",
!           Comm.Appl.Num.Meth. Vol.2, 1986, pp. 401-410.
      implicit none
!--subroutine-parameters:
      real(kind=8) :: dx, dy, gg, v(3,3,0:3), vnu(3,3,0:3), xw(3), xp(2)

      ! dx       discretisation step-size in first dimension
      ! dy       discretisation step-size in second dimension
      ! gg       modulus of rigidity of half space considered [N/mm2]
      ! v(i,j,k) elements (1:3,1:3, 0 ): displacements u_j [mm] due to unit load (1 N/mm2) in direction i
      !          elements (1:3,1:3,1:3): displacement gradients u_{j,k} [mm/mm] due to unit load in direction i
      ! vnu      same as v(i,j,k), contributions of contraction (Poisson ratio effect)
      ! xp       center-coordinate of surface element with unit traction
      ! xw       subsurface location where displacements are to be calculated

!--local variables:
      real(kind=8), parameter :: epsrel = 5d-7

      ! epsrel   relative tolerance for computing length of distance vector y, to avoid ||y|| == |y_i|
      !          in any of the components i. Note: this parameter may be machine-dependent.

      integer      :: kronec(3,3), i, ip(3), j, k, sgn(3), l, jy, jx
      real(kind=8) :: epsy, pi, q, w, weps, al(3), at(3), a(3,3,0:3), t(3,3,0:3), wm(0:3), y(3), yeps(3)

      ! P o i n t   1  OF THE ALGORITHM.

      pi = 4d0 * atan(1d0)

      ! initialize sgn (sign per direction)

      sgn = (/ -1, -1, 1 /)

      ! initialize ip (next direction?)

      ip = (/ 2, 3, 1 /)

      ! initialize Kronecker delta matrix

      kronec(1,1:3) = (/ 1, 0, 0 /)
      kronec(2,1:3) = (/ 0, 1, 0 /)
      kronec(3,1:3) = (/ 0, 0, 1 /)

      ! P o i n t  4, 5, 6 OF THE ALGORITHM

      ! compute initial difference vector y and distance w

      y(1) = xp(1) - xw(1)
      y(2) = xp(2) - xw(2)
      y(3) =         xw(3)
      w    = max(1d-12, dsqrt(y(1)**2 + y(2)**2 + y(3)**2))

      ! compute perturbed difference vector yeps and distance weps, such that |yeps_i| >= epsy and
      ! keeping original sign of y_i.
      ! Note: not needed for y(3) == z (why?).
      ! This results in weps > |yeps_i|, weps-yeps_i bounded away from 0.

      epsy = epsrel*w
      do i = 1, 2
         if (y(i).ge.0d0) then
            yeps(i) = max( epsy, y(i))
         else
            yeps(i) = min(-epsy, y(i))
         endif
      enddo
      yeps(3) = y(3)

      weps = max(1d-12, dsqrt(yeps(1)**2 + yeps(2)**2 + yeps(3)**2))

      ! P o i n t  9  OF THE ALGORITHM: "NEAR" CALCULATION.


      do i = 1, 3
         do j = 1, 3
            do k = 0, 3
               vnu(i,j,k) = 0
               v(i,j,k) = 0
            enddo
         enddo
      enddo

      ! P o i n t s  10 THROUGH 15   OF THE ALGORITHM.

      ! process the four corners of the element around xp

      do jx = -1, 1, 2
         do jy = -1, 1, 2

            ! compute initial difference vector y and distance w

            y(1) = xp(1) + jx*dx/2d0 - xw(1)
            y(2) = xp(2) + jy*dy/2d0 - xw(2)
            y(3) =                     xw(3)
            w = max(1d-12, dsqrt(y(1)**2 + y(2)**2 + y(3)**2))

            ! compute perturbed difference vector yeps and distance weps, ensure that 
            ! |yeps_i| >= epsy (in each component) and keep original sign on y_i

            epsy = epsrel*w
            do i = 1, 3
               if (y(i).ge.0d0) then
                  yeps(i) = max( epsy, y(i))
               else
                  yeps(i) = min(-epsy, y(i))
               endif
            enddo

            ! compute weps for yeps, which has |weps - yeps_i| > 0 for all i

            weps = max(1d-12, dsqrt(yeps(1)**2 + yeps(2)**2 + yeps(3)**2))

            wm(0) = weps

            ! pre-compute logarithms, arctangens, distance/r

            do k = 1, 3
               al(k) = dlog(y(k) + weps)
               at(k) = datan((y(ip(k)) + y(ip(ip(k))) + w) / yeps(k))
               wm(k) = sgn(k) * y(k) / weps

               ! if (jx.eq.1 .and. jy.eq.1) write(*,'(a,i1,a,3f12.6,a,f12.6,a,f12.6)')         &
               !    'k=',k,': y=',y(ip(k)),y(ip(ip(k))),w,', yeps(k)=',yeps(k),', at=',at(k)
            enddo

            ! compute intermediate expressions a, t using al, at and various distances

            do i = 1, 3
               do j = 1, 3
                  a(i,j,0) = y(i)*al(j)
                  t(i,j,0) = y(i)*at(j)
                  do k = 1, 3
                     a(i,j,k) = sgn(k) * (                                                   &
     &                          kronec(i,k)*al(j) +                                          &
     &                          y(i)*(kronec(j,k)*w + y(k)) / (weps * (y(j) + weps))         &
     &                       )
                     if (isnan(a(i,j,k))) then
                        write(*,*) 'nan in a at',i,j,k
                        ! write(*,*) 'al,yi,yj,yk,weps=',al(j),yeps(i), yeps(j), yeps(k), weps
                        call abort_run()
                     endif
                     t(i,j,k) = sgn(k) * (                                                   &
     &                          kronec(i,k)*at(j) +                                          &
     &                          y(i)* ( y(j)*(y(k) + w) -                                    &
     &                                  kronec(j,k)*w*(y(1)+y(2)+y(3)+w) )                   &
                                    / (2*weps*(y(ip(j))+weps) * (y(ip(ip(j)))+weps))         &
     &                          )
                     if (isnan(t(i,j,k))) then
                        write(*,*) 'nan in t at',i,j,k
                        write(*,*) 'noemer=', weps,'*',y(ip(j))+weps,'*',y(ip(ip(j)))+weps
                        write(*,*) 'y=',yeps,', ip=',ip(j),ip(ip(j))
                        call abort_run()
                     endif
                  enddo
               enddo
            enddo
            ! do k=0, 0
            !    write(*,'(a,i1,a,3(3f13.6,:,/,11x))') 'a(:,:,k=',k,')=',((a(i,j,k),j=1,3),i=1,3)
            ! enddo
            ! do k=0, 3
            !    write(*,'(a,i1,a,3(3f13.6,:,3x))') 't(:,:,k=',k,')=',((t(i,j,k),j=1,3),i=1,3)
            ! enddo

            ! compute displacements and displacement gradients by adding all appropriate terms
            !    using a and t.

            do k = 0, 3
               do i = 1, 2
                  l = 3 - i
                  q = a(i,l,k) + 2*a(l,i,k) + 4*t(3,3,k) + a(i,l,k) - 2*t(3,i,k)
                  v(i,i,k) = v(i,i,k) + jx*jy*q
                  ! if (i.eq.2 .and. k.eq.0)                                                   &
                  !    write(*,'(5f14.9)') a(i,l,k), 2*a(l,i,k), 4*t(3,3,k), a(i,l,k), -2*t(3,i,k)
                  ! if (i.eq.2 .and. k.eq.0) &
                  !    write(*,'(2(a,f13.6))') 'add-0',jx*jy*q,': new v=',v(i,i,k)
                  q        = -2*(a(i,l,k) - 2*t(3,i,k))
                  vnu(i,i,k) = vnu(i,i,k) + q*jx*jy
                  q        = -a(3,3,k)
                  v(i,l,k) = v(i,l,k) + jx*jy*q
                  ! if (i.eq.3 .and. l.eq.1 .and. k.eq.3)                                      &
                  !    write(*,'(2(a,f13.6))') 'add-1',jx*jy*q,', val=',v(i,l,k)
                  q        = 2*(a(3,3,k) - wm(k))
                  vnu(i,l,k) = vnu(i,l,k) + q*jx*jy
                  q        = 2*(a(l,3,k) + a(3,l,k) + 2*t(i,i,k))
                  vnu(i,3,k) = vnu(i,3,k) + q*jx*jy
                  vnu(3,i,k) = -vnu(i,3,k)
                  q        = -a(l,3,k) - 2*t(i,i,k)
                  v(i,3,k) = v(i,3,k) + q*jx*jy
                  v(3,i,k) = v(3,i,k) + (2*a(3,l,k) - q)*jx*jy
                  ! if (i.eq.1 .and. k.eq.3)                                                   &
                  !    write(*,'(a,2f13.6,a,f13.6)') 'add-2',a(3,l,k),q,', val=',v(3,i,k)
               enddo
               q = -2*(a(1,2,k) + a(2,1,k) + 2*t(3,3,k))
               vnu(3,3,k) = vnu(3,3,k) + q*jx*jy
               q = -2*t(3,3,k) + 2*(a(1,2,k) + a(2,1,k) + 2*t(3,3,k))
               v(3,3,k) = v(3,3,k) + q*jx*jy
            enddo

      ! enddo (four corners of element)

         enddo
      enddo
      ! write(*,'(a,3f12.8)') 'near calc.:',(v(l,1,3),l=1,3)

      ! P o i n t  16, 17  OF THE ALGORITHM.

      do i = 1, 3
         do j = 1, 3
            do k = 0, 3
               v  (i,j,k) = v  (i,j,k) / (4 * pi * gg)
               vnu(i,j,k) = vnu(i,j,k) / (4 * pi * gg)
            enddo
         enddo
      enddo

end subroutine stres1_pcwcns

!------------------------------------------------------------------------------------------------------------

subroutine stres1_bilin(dx, dy, gg, v, vnu, xw, xp)
!--purpose: This subroutine stres1_bilin calculates the displacements and the displacement gradients in
!           an elastic half-space due to a surface load bilinearly distributed over an arbitrary rectangle
!           with sides parallel to the x,y-axes (directions 1,2) of a Cartesian coordinate system, in
!           which the half-space is given by "xw(3).ge.0".  See [Wekken2018a-subsurf] (submitted).
      implicit none
!--subroutine parameters:
      real(kind=8) , intent(in)  :: dx, dy, gg, xp(2), xw(3)
      real(kind=8) , intent(out) :: v(3,3,0:3), vnu(3,3,0:3)

      ! dx       discretisation step-size in first dimension
      ! dy       discretisation step-size in second dimension
      ! gg       modulus of rigidity of half space considered [N/mm2]
      ! xp       center-coordinate of surface element with unit traction
      ! xw       subsurface location where displacements are to be calculated
      ! v(i,j,k) elements (1:3,1:3, 0 ): displacements u_j [mm] due to unit load (1 N/mm2) in direction i
      !          elements (1:3,1:3,1:3): displacement gradients u_{j,k} [mm/mm] for unit load in direction i
      ! vnu      same as v(i,j,k), contributions of contraction (Poisson ratio effect)

!--local variables:
      integer      :: i, j, k, jy, jx
      real(kind=8) :: pi, eps, y(3), yeps(3), rzero, r, reps, fact
      type(t_xyzr) :: xyzr

      eps = 1d-7        ! lower bound for arguments of log's and atan's
      pi = 4d0 * datan(1d0)

      ! initialize

      v(1:3,1:3,0:3)   = 1.d99
      vnu(1:3,1:3,0:3) = 1.d99

      do jx = -1, 1
         do jy = -1, 1

            ! compute initial difference vector y and distance w

            xyzr%y(jx,jy,1) = xp(1) + jx*dx - xw(1)
            xyzr%y(jx,jy,2) = xp(2) + jy*dy - xw(2)
            xyzr%y(jx,jy,3) = xw(3)

            y(1) = xyzr%y(jx,jy,1)
            y(2) = xyzr%y(jx,jy,2)
            y(3) = xyzr%y(jx,jy,3)

            rzero = dsqrt(y(1)**2 + y(2)**2 + y(3)**2)
            r = max(eps, rzero)
            xyzr%r(jx,jy,0) = rzero
            xyzr%r3(jx,jy,0) = rzero**3

            ! Compute perturbed difference vector yeps and distance reps, such that |yeps_i| >= eps and
            ! keeping original sign of y_i.
            ! This results in reps > |yeps_i|, reps-reps_i bounded away from 0.

            do i = 1, 2
               if (y(i).ge.0.d0) then
                  yeps(i) = max( eps, y(i))
               else
                  yeps(i) = min(-eps, y(i))
               endif
            enddo
            yeps(3) = max(eps,y(3))
            reps = dsqrt(yeps(1)**2 + yeps(2)**2 + yeps(3)**2)

            ! Pre-compute radius, logarithms(x+r) and arctangens(x/y), (xy/zr), (y+z+r)/x .

            do k = 1, 3
               xyzr%r(jx,jy,k) = y(k)/r
               xyzr%r3(jx,jy,k) = 3.d0*y(k)*rzero
            enddo

            ! LOG(x + r)

            xyzr%logx(jx,jy,0) = dlog(y(1) + reps)
            xyzr%logy(jx,jy,0) = dlog(y(2) + reps)
            xyzr%logz(jx,jy,0) = dlog(y(3) + reps)

            do k = 1, 3
               if (k .eq. 1) then
                  xyzr%logx(jx,jy,1) = 1.d0/r
               else
                  xyzr%logx(jx,jy,k) = y(k)/(r*(y(1)+reps))
               endif
               if (k .eq. 2) then
                  xyzr%logy(jx,jy,2) = 1.d0/r
               else
                  xyzr%logy(jx,jy,k) = y(k)/(r*(y(2)+reps))
               endif
               if (k .eq. 3) then
                  xyzr%logz(jx,jy,3) = 1.d0/r
               else
                  xyzr%logz(jx,jy,k) = y(k)/(r*(y(3)+reps))
               endif
            enddo

            ! ATAN( (y+z+r) / x )

            xyzr%tanKx(jx,jy,0) = datan((y(2)+y(3)+rzero)/yeps(1))
            xyzr%tanKy(jx,jy,0) = datan((y(1)+y(3)+rzero)/yeps(2))
            xyzr%tanKz(jx,jy,0) = datan((y(1)+y(2)+rzero)/yeps(3))

            xyzr%tanKx(jx,jy,1) = -(y(2)*(y(2)+r) + y(3)*(y(3)+r)) / (2.d0*r*(y(2)+reps)*(y(3)+reps))
            xyzr%tanKy(jx,jy,2) = -(y(1)*(y(1)+r) + y(3)*(y(3)+r)) / (2.d0*r*(y(1)+reps)*(y(3)+reps))
            xyzr%tanKz(jx,jy,3) = -(y(1)*(y(1)+r) + y(2)*(y(2)+r)) / (2.d0*r*(y(1)+reps)*(y(2)+reps))

            xyzr%tanKx(jx,jy,2) = y(1) / (2.d0*r*(y(3)+reps))
            xyzr%tanKx(jx,jy,3) = y(1) / (2.d0*r*(y(2)+reps))

            xyzr%tanKy(jx,jy,1) = y(2) / (2.d0*r*(y(3)+reps))
            xyzr%tanKy(jx,jy,3) = y(2) / (2.d0*r*(y(1)+reps))

            xyzr%tanKz(jx,jy,1) = y(3) / (2.d0*r*(y(2)+reps))
            xyzr%tanKz(jx,jy,2) = y(3) / (2.d0*r*(y(1)+reps))
         enddo ! jy
      enddo ! jx

      ! Calculate the influence using Wekken formulae (extended Kalker)

      do k = 0, 3
        v(1,1,k) = Hw(xyzr,11,k,2)
        v(1,2,k) = Hw(xyzr,12,k,2)
        v(1,3,k) = Hw(xyzr,13,k,2)
        v(2,1,k) = v(1,2,k)
        v(2,2,k) = Hw(xyzr,22,k,2)
        v(2,3,k) = Hw(xyzr,23,k,2)
        v(3,1,k) = Hw(xyzr,31,k,2)
        v(3,2,k) = Hw(xyzr,32,k,2)
        v(3,3,k) = Hw(xyzr,33,k,2)

         vnu(1,1,k) = Hw(xyzr,11,k,3)
         vnu(1,2,k) = Hw(xyzr,12,k,3)
         vnu(1,3,k) = Hw(xyzr,13,k,3)
         vnu(2,1,k) = vnu(1,2,k)
         vnu(2,2,k) = Hw(xyzr,22,k,3)
         vnu(2,3,k) = Hw(xyzr,23,k,3)
         vnu(3,1,k) = -vnu(1,3,k)
         vnu(3,2,k) = -vnu(2,3,k)
         vnu(3,3,k) = Hw(xyzr,33,k,3)
      enddo ! k

      fact = 1.d0 / (4.d0 * pi * gg * dx*dy)
      do i = 1, 3
         do j = 1, 3
            do k = 0, 3
               v  (i,j,k) = v  (i,j,k) * fact
               vnu(i,j,k) = vnu(i,j,k) * fact
            enddo
         enddo
      enddo

end subroutine stres1_bilin

!------------------------------------------------------------------------------------------------------------

   function Hw(xyzr,ij,k,vn)
!--purpose: combine 4 Gw-functions, summing over the 4 squares that make up the full area of integration
      implicit none
!--function result:
      real(kind=8) :: Hw
!--function arguments:
      type(t_xyzr) :: xyzr
      integer      :: ij, k, vn
!--no local variables

      Hw =  Gw(-xyzr%y( 1,0,1),-xyzr%y(0, 1,2), xyzr,  0, 0, ij,k,vn)  &
          - Gw(-xyzr%y( 1,0,1),-xyzr%y(0,-1,2), xyzr,  0,-1, ij,k,vn)  &
          - Gw(-xyzr%y(-1,0,1),-xyzr%y(0, 1,2), xyzr, -1, 0, ij,k,vn)  &
          + Gw(-xyzr%y(-1,0,1),-xyzr%y(0,-1,2), xyzr, -1,-1, ij,k,vn)

   end function Hw

!------------------------------------------------------------------------------------------------------------

   function Gw(a,b,xyzr,ix,iy,ij,k,vn)
!--purpose: combines 4 iFw-functions, summing over the constant, linear in x, linear in y and bilinear
!           part of A*phi.
      implicit none
!--function result:
      real(kind=8) :: Gw
!--subroutine arguments:
      type(t_xyzr) :: xyzr
      real(kind=8) :: a, b
      integer      :: ix, iy, ij, k, vn
!--no local variables

      Gw =  a*b*iFw(xyzr,ix,iy,ij,00,k,vn)  &
          +   b*iFw(xyzr,ix,iy,ij,10,k,vn)  &
          + a*  iFw(xyzr,ix,iy,ij,01,k,vn)  &
          +     iFw(xyzr,ix,iy,ij,11,k,vn)
      if (k.eq.1) Gw = Gw + b*iFw(xyzr,ix,iy,ij,00,0,vn)  &
                          +   iFw(xyzr,ix,iy,ij,01,0,vn)
      if (k.eq.2) Gw = Gw + a*iFw(xyzr,ix,iy,ij,00,0,vn)  &
                             +   iFw(xyzr,ix,iy,ij,10,0,vn)
   end function Gw

!------------------------------------------------------------------------------------------------------------

   function iFw(xyzr,ix,iy,ij,mn,k,vn)
!--purpose: determine the integral of F_w over a square centered around (xi,yi) with width/height dx/dy
      implicit none
!--function result:
      real(kind=8) :: iFw
!--subroutine arguments:
      type(t_xyzr) :: xyzr
      integer      :: ix, iy, ij, mn, k, vn
!--no local variables
      iFw =  Fw(xyzr,ix  ,iy  ,ij,mn,k,vn)  &
           - Fw(xyzr,ix+1,iy  ,ij,mn,k,vn)  &
           - Fw(xyzr,ix  ,iy+1,ij,mn,k,vn)  &
           + Fw(xyzr,ix+1,iy+1,ij,mn,k,vn)

   end function iFw

!------------------------------------------------------------------------------------------------------------

function Fw(xyzr, ix, iy, ij, mn, k, vn)
!--purpose: determine the value of F_w, 'direct integrals / extended Kalker', at the point (x,y)
   implicit none
!--function result:
      real(kind=8) :: Fw
!--subroutine arguments:
   type(t_xyzr),               target  :: xyzr
   integer                             :: ix, iy, ij, mn, k, vn
!--local (pointer) variables:
   real(kind=8), dimension(:), pointer :: r, r3, Lx, Ly, Lz, Kx, Ky, Kz

   associate(x => xyzr%y(ix,iy,1), y => xyzr%y(ix,iy,2), z => xyzr%y(ix,iy,3))
          
   r(0:)  => xyzr%r(ix,iy,0:3)
   r3(0:) => xyzr%r3(ix,iy,0:3)

   Lx(0:) => xyzr%logx(ix,iy,0:3)
   Ly(0:) => xyzr%logy(ix,iy,0:3)
   Lz(0:) => xyzr%logz(ix,iy,0:3)

   Kx(0:) => xyzr%tanKx(ix,iy,0:3)
   Ky(0:) => xyzr%tanKy(ix,iy,0:3)
   Kz(0:) => xyzr%tanKz(ix,iy,0:3)

   Fw = 0.d0

   if (vn .eq. 2) then ! calculate the v part
      if     (ij .eq. 11) then
         if     (mn .eq. 00) then
            Fw = 2.d0*( x*Ly(k) + y*Lx(k) - z*(Kx(k) - 2.d0*Kz(k)))
            if (k.eq.1) Fw = Fw + 2.d0*Ly(0)
            if (k.eq.2) Fw = Fw + 2.d0*Lx(0)
            if (k.eq.3) Fw = Fw - 2.d0*(Kx(0) - 2.d0*Kz(0))
         elseif (mn .eq. 10) then
            Fw = (x**2 + 2.d0*z**2)*Ly(k) + y*z*Lz(k) + y*r(k)
            if (k.eq.1) Fw = Fw + 2.d0*x*Ly(0)
            if (k.eq.2) Fw = Fw + z*Lz(0) + r(0)
            if (k.eq.3) Fw = Fw + 4.d0*z*Ly(0) + y*Lz(0)
         elseif (mn .eq. 01) then
            Fw = (y**2 + z**2)*Lx(k) - x*z*Lz(k) + x*r(k)
            if (k.eq.1) Fw = Fw - z*Lz(0) + r(0)
            if (k.eq.2) Fw = Fw + 2.d0*y*Lx(0)
            if (k.eq.3) Fw = Fw + 2.d0*z*Lx(0) - x*Lz(0)
         elseif (mn .eq. 11) then
            Fw = 0.5d0*( (y**2 - x**2)*z*Lz(k) + (7.d0/3.d0)*r3(k) - (x**2 + y**2)*r(k) )
            if (k.eq.1) Fw = Fw - x*z*Lz(0) - x*r(0)
            if (k.eq.2) Fw = Fw + y*z*Lz(0) - y*r(0)
            if (k.eq.3) Fw = Fw + 0.5d0*(y**2 - x**2)*Lz(0)
         endif
      elseif (ij .eq. 12) then
         if     (mn .eq. 00) then
            Fw = - z*Lz(k)
            if (k.eq.3) Fw = Fw - Lz(0)
         elseif (mn .eq. 10) then
            Fw = z**2*Lx(k) + 2.d0*y*z*Ky(k)
            if (k.eq.2) Fw = Fw + 2.d0*z*Ky(0)
            if (k.eq.3) Fw = Fw + 2.d0*( z*Lx(0) + y*Ky(0) )
         elseif (mn .eq. 01) then
            Fw = z**2*Ly(k) + 2.d0*x*z*Kx(k)
            if (k.eq.1) Fw = Fw + 2.d0*z*Kx(0)
            if (k.eq.3) Fw = Fw + 2.d0*( z*Ly(0) + x*Kx(0) )
         elseif (mn .eq. 11) then
            Fw = z*( x**2*Kx(k) + y**2*Ky(k) - z**2*Kz(k) )
            if (k.eq.0) Fw = Fw - 0.5d0*x*y*z
            if (k.eq.1) Fw = Fw + 2.d0*x*z*Kx(0)
            if (k.eq.2) Fw = Fw + 2.d0*y*z*ky(0)
            if (k.eq.3) Fw = Fw + x**2*Kx(0) + y**2*Ky(0) - 3.d0*z**2*Kz(0) - 0.5d0*x*y
         endif
      elseif (ij .eq. 13) then
         if     (mn .eq. 00) then
            Fw = -y*Lz(k) - 2.d0*x*Kx(k)
            if (k.eq.1) Fw = Fw - 2.d0*Kx(0)
            if (k.eq.2) Fw = Fw - Lz(0)
         elseif (mn .eq. 10) then
            Fw = y**2*Ky(k) - x**2*Kx(k) - z**2*Kz(k)
            if (k.eq.0) Fw = Fw - 0.5d0*x*y
            if (k.eq.1) Fw = Fw - 2.d0*x*Kx(0)
            if (k.eq.2) Fw = Fw + 2.d0*y*Ky(0)
            if (k.eq.3) Fw = Fw - 2.d0*z*Kz(0)
         elseif (mn .eq. 01) then
            Fw = 0.5d0*( z*r(k) - (x**2 + y**2)*Lz(k) )
            if (k.eq.1) Fw = Fw - x*Lz(0)
            if (k.eq.2) Fw = Fw - y*Lz(0)
            if (k.eq.3) Fw = Fw + 0.5d0*r(0)
         elseif (mn .eq. 11) then
            Fw = ( 2.d0*y**3*Ky(k) - z**3*Lx(k) - x**3*Lz(k) + x*z*r(k) )/3.d0
            if (k.eq.0) Fw = Fw - x*y**2/3.d0
            if (k.eq.1) Fw = Fw - x**2*Lz(0) + z*r(0)/3.d0
            if (k.eq.2) Fw = Fw + 2.d0*( y**2*Ky(0) - x*y/3.d0 )
            if (k.eq.3) Fw = Fw - z**2*Lx(0) + x*r(0)/3.d0
         endif

      elseif (ij .eq. 22) then
         if     (mn .eq. 00) then
            Fw = 2.d0*( y*Lx(k) + x*Ly(k) - z*(Ky(k) - 2.d0*Kz(k)) )
            if (k.eq.1) Fw = Fw + 2.d0*Ly(0)
            if (k.eq.2) Fw = Fw + 2.d0*Lx(0)
            if (k.eq.3) Fw = Fw - 2.d0*(Ky(0) - 2.d0*Kz(0))
         elseif (mn .eq. 10) then
            Fw = (x**2 + z**2)*Ly(k) - y*z*Lz(k) + y*r(k)
            if (k.eq.1) Fw = Fw + 2.d0*x*Ly(0)
            if (k.eq.2) Fw = Fw - z*Lz(0) + r(0)
            if (k.eq.3) Fw = Fw + 2.d0*z*Ly(0) - y*Lz(0)
         elseif (mn .eq. 01) then
            Fw = (y**2 + 2.d0*z**2)*Lx(k) + x*z*Lz(k) + x*r(k)
            if (k.eq.1) Fw = Fw + z*Lz(0) + r(0)
            if (k.eq.2) Fw = Fw + 2.d0*y*Lx(0)
            if (k.eq.3) Fw = Fw + 4.d0*z*Lx(0) + x*Lz(0)
         elseif (mn .eq. 11) then
            Fw = 0.5d0*( (x**2 - y**2)*z*Lz(k) + (7.d0/3.d0)*r3(k) - (x**2 + y**2)*r(k) )
            if (k.eq.1) Fw = Fw + x*z*Lz(0) - x*r(0)
            if (k.eq.2) Fw = Fw - y*z*Lz(0) - y*r(0)
            if (k.eq.3) Fw = Fw + 0.5d0*(x**2 - y**2)*Lz(0)
         endif
      elseif (ij .eq. 23) then
         if     (mn .eq. 00) then
            Fw = - x*Lz(k) - 2.d0*y*Ky(k)
            if (k.eq.1) Fw = Fw - Lz(0)
            if (k.eq.2) Fw = Fw - 2.d0*Ky(0)
         elseif (mn .eq. 10) then
            Fw = -0.5d0*( (x**2 + y**2)*Lz(k) - z*r(k) )
            if (k.eq.1) Fw = Fw - x*Lz(0)
            if (k.eq.2) Fw = Fw - y*Lz(0)
            if (k.eq.3) Fw = Fw + 0.5d0*r(0)
         elseif (mn .eq. 01) then
            Fw = x**2*Kx(k) - y**2*Ky(k) - z**2*Kz(k)
            if (k.eq.0) Fw = Fw - 0.5d0*x*y
            if (k.eq.1) Fw = Fw + 2.d0*x*Kx(0)
            if (k.eq.2) Fw = Fw - 2.d0*y*Ky(0)
            if (k.eq.3) Fw = Fw - 2.d0*z*Kz(0)
         elseif (mn .eq. 11) then
            Fw = ( 2.d0*x**3*Kx(k) - z**3*Ly(k) - y**3*Lz(k) + y*z*r(k) )/3.d0
            if (k.eq.0) Fw = Fw - y*x**2/3.d0
            if (k.eq.1) Fw = Fw + 2.d0*( x**2*Kx(0) - x*y/3.d0 )
            if (k.eq.2) Fw = Fw - y**2*Lz(0) + z*r(0)/3.d0
            if (k.eq.3) Fw = Fw - z**2*Ly(0) + y*r(0)/3.d0
         endif

      elseif (ij .eq. 31) then
         if     (mn .eq. 00) then
            Fw = 2.d0*z*Ly(k) + y*Lz(k) + 2.d0*x*Kx(k)
            if (k.eq.1) Fw = Fw + 2.d0*Kx(0)
            if (k.eq.2) Fw = Fw + Lz(0)
            if (k.eq.3) Fw = Fw + 2.d0*Ly(0)
         elseif (mn .eq. 10) then
            Fw = x**2*Kx(k) - 2.d0*y*z*Lx(k) - y**2*Ky(k) - 3.d0*z**2*Kz(k)
            if (k.eq.0) Fw = Fw + 0.5d0*x*y
            if (k.eq.1) Fw = Fw + 2.d0*x*Kx(0)
            if (k.eq.2) Fw = Fw - 2.d0*( z*Lx(0) + y*Ky(0) )
            if (k.eq.3) Fw = Fw - 2.d0*y*Lx(0) - 6.d0*z*Kz(0)
         elseif (mn .eq. 01) then
            Fw = 0.5d0*( (x**2 + y**2)*Lz(k) + 3.d0*z*r(k) )
            if (k.eq.1) Fw = Fw + x*Lz(0)
            if (k.eq.2) Fw = Fw + y*Lz(0)
            if (k.eq.3) Fw = Fw + 1.5d0*r(0)
         elseif (mn .eq. 11) then
            Fw = -z*(y**2 + (2.d0/3.d0)*z**2)*Lx(k) + ( x**3*Lz(k) + 2.d0*( x*z*r(k) - y**3*Ky(k) ) )/3.d0
            if (k.eq.0) Fw = Fw + x*y**2/3.d0
            if (k.eq.1) Fw = Fw + x**2*Lz(0) + 2.d0*z*r(0)/3.d0
            if (k.eq.2) Fw = Fw - 2.d0*( y*z*Lx(0) + y**2*Ky(0) - x*y/3.d0 )
            if (k.eq.3) Fw = Fw - (y**2 + 2.d0*z**2)*Lx(0) + 2.d0*x*r(0)/3.d0
         endif
      elseif (ij .eq. 32) then
         if     (mn .eq. 00) then
            Fw = 2.d0*z*Lx(k) + x*Lz(k) + 2.d0*y*Ky(k)
            if (k.eq.1) Fw = Fw + Lz(0)
            if (k.eq.2) Fw = Fw + 2.d0*Ky(0)
            if (k.eq.3) Fw = Fw + 2.d0*Lx(0)
         elseif (mn .eq. 10) then
            Fw = 0.5d0*( (y**2 + x**2)*Lz(k) + 3.d0*z*r(k) )
            if (k.eq.1) Fw = Fw + x*Lz(0)
            if (k.eq.2) Fw = Fw + y*Lz(0)
            if (k.eq.3) Fw = Fw + 1.5d0*r(0)
         elseif (mn .eq. 01) then
            Fw = - 2.d0*x*z*Ly(k) + y**2*Ky(k) - x**2*Kx(k) - 3.d0*z**2*Kz(k)
            if (k.eq.0) Fw = Fw + 0.5d0*x*y
            if (k.eq.1) Fw = Fw - 2.d0*( z*Ly(0) + x*Kx(0) )
            if (k.eq.2) Fw = Fw + 2.d0*y*Ky(0)
            if (k.eq.3) Fw = Fw - 2.d0*x*Ly(0) - 6.d0*z*Kz(0)
         elseif (mn .eq. 11) then
            Fw = -z*(x**2 + (2.d0/3.d0)*z**2)*Ly(k) + ( y**3*Lz(k) + 2.d0*( y*z*r(k) - x**3*Kx(k) ) )/3.d0
            if (k.eq.0) Fw = Fw - y*x**2/3.d0
            if (k.eq.1) Fw = Fw - 2.d0*( x*z*Ly(0) + x**2*Kx(0) - x*y/3.d0 )
            if (k.eq.2) Fw = Fw + y**2*Lz(0) + 2.d0*z*r(0)/3.d0
            if (k.eq.3) Fw = Fw - (x**2 + 2.d0*z**2)*Ly(0) + 2.d0*y*r(0)/3.d0
         endif
      elseif (ij .eq. 33) then
         if     (mn .eq. 00) then
            Fw =  2.d0*( y*Lx(k) + x*Ly(k) + z*Kz(k) )
            if (k.eq.1) Fw = Fw + 2.d0*Ly(0)
            if (k.eq.2) Fw = Fw + 2.d0*Lx(0)
            if (k.eq.3) Fw = Fw + 2.d0*Kz(0)
         elseif (mn .eq. 10) then
            Fw = x**2*Ly(k) + y*r(k)
            if (k.eq.1) Fw = Fw + 2.d0*x*Ly(0)
            if (k.eq.2) Fw = Fw + r(0)
         elseif (mn .eq. 01) then
            Fw = y**2*Lx(k) + x*r(k)
            if (k.eq.1) Fw = Fw + r(0)
            if (k.eq.2) Fw = Fw + 2.d0*y*Lx(0)
         elseif (mn .eq. 11) then
            Fw = (2.d0/3.d0)*r3(k) - z**2*r(k)
            if (k.eq.3) Fw = Fw - 2.d0*z*r(0)
         endif

      else
         Fw = 1.d100
         write(*,*) 'WARNING: in the F_v function the numbers did not match any case'
      endif

   else ! calculate the vnu part
      if     (ij .eq. 11) then
         if     (mn .eq. 00) then
            Fw = -2.d0*x*Ly(k) + 4.d0*z*Kx(k)
            if (k.eq.1) Fw = Fw - 2.d0*Ly(0)
            if (k.eq.3) Fw = Fw + 4.d0*Kx(0)
         elseif (mn .eq. 10) then
            Fw = y*r(k) - 2.d0*y*z*Lz(k) - (x**2+z**2)*Ly(k)
            if (k.eq.1) Fw = Fw - 2.d0*x*Ly(0)
            if (k.eq.2) Fw = Fw + r(0) - 2.d0*z*Lz(0)
            if (k.eq.3) Fw = Fw - 2.d0*( y*Lz(0) + z*Ly(0) )
         elseif (mn .eq. 01) then
            Fw  = 2.d0*( x*z*Lz(k) - x*r(k) )
            if (k.eq.1) Fw = Fw + 2.d0*( z*Lz(0) - r(0) )
            if (k.eq.3) Fw = Fw + 2.d0*x*Lz(0)
         elseif (mn .eq. 11) then
            Fw = (x**2 - y**2)*z*Lz(k) + (2.d0*y**2 + z**2)*r(k) - (4.d0/3.d0)*r3(k)
            if (k.eq.1) Fw = Fw + 2.d0*x*z*Lz(0)
            if (k.eq.2) Fw = Fw - 2.d0*y*z*Lz(0) + 4.d0*y*r(0)
            if (k.eq.3) Fw = Fw + (x**2 - y**2)*Lz(0) + 2.d0*z*r(0)
         endif
      elseif (ij .eq. 12) then
         if     (mn .eq. 00) then
            Fw = 2.d0*( z*Lz(k) - r(k) )
            if (k.eq.3) Fw = Fw + 2.d0*Lz(0)
         elseif (mn .eq. 10) then
            Fw = (y**2 - z**2)*Lx(k) - 4.d0*y*z*Ky(k) - x*r(k)
            if (k.eq.1) Fw = Fw - r(0)
            if (k.eq.2) Fw = Fw + 2.d0*y*Lx(0) - 4.d0*z*Ky(0)
            if (k.eq.3) Fw = Fw - 2.d0*z*Lx(0) - 4.d0*y*Ky(0)
         elseif (mn .eq. 01) then
            Fw = (x**2 - z**2)*Ly(k) - 4.d0*x*z*Kx(k) - y*r(k)
            if (k.eq.1) Fw = Fw + 2.d0*x*Ly(0) - 4.d0*z*Kx(0)
            if (k.eq.2) Fw = Fw - r(0)
            if (k.eq.3) Fw = Fw - 2.d0*z*Ly(0) - 4.d0*x*Kx(0)
         elseif (mn .eq. 11) then
            Fw = (2.d0/3.d0)*( x**3*Ly(k) + y**3*Lx(k) - x*y*r(k) + z**3*Kz(k) )                       &
                - 2.d0*z*( x**2*Kx(k) + y**2*Ky(k) )
            if (k.eq.0) Fw = Fw + x*y*z
            if (k.eq.1) Fw = Fw + 2.d0*( x**2*Ly(0) - y*r(0)/3.d0 ) - 4.d0*x*z*Kx(0)
            if (k.eq.2) Fw = Fw + 2.d0*( y**2*Lx(0) - x*r(0)/3.d0 ) - 4.d0*y*z*Ky(0)
            if (k.eq.3) Fw = Fw - 2.d0*( x**2*Kx(0) + y**2*Ky(0) - z**2*Kz(0) ) + x*y
         endif
      elseif (ij .eq. 13) then
         if     (mn .eq. 00) then
            Fw = 2.d0*( y*Lz(k) + z*Ly(k) ) + 4.d0*x*Kx(k)
            if (k.eq.1) Fw = Fw + 4.d0*Kx(0)
            if (k.eq.2) Fw = Fw + 2.d0*Lz(0)
            if (k.eq.3) Fw = Fw + 2.d0*Ly(0)
         elseif (mn .eq. 10) then
            Fw = 2.d0*( x**2*Kx(k) - y**2*Ky(k) - y*z*Lx(k) - z**2*Kz(k) )
            if (k.eq.0) Fw = Fw + x*y
            if (k.eq.1) Fw = Fw + 4.d0*x*Kx(0)
            if (k.eq.2) Fw = Fw - 4.d0*y*Ky(0) - 2.d0*z*Lx(0)
            if (k.eq.3) Fw = Fw - 4.d0*z*Kz(0) - 2.d0*y*Lx(0)
         elseif (mn .eq. 01) then
            Fw  = (x**2 + y**2)*Lz(k) + z*r(k)
            if (k.eq.1) Fw = Fw + 2.d0*x*Lz(0)
            if (k.eq.2) Fw = Fw + 2.d0*y*Lz(0)
            if (k.eq.3) Fw = Fw + r(0)
         elseif (mn .eq. 11) then
            Fw = ( x*z*r(k) - (3.d0*y**2 + z**2)*z*Lx(k) + 2.d0*x**3*Lz(k) - 4.d0*y**3*Ky(k) )/3.d0
            if (k.eq.0) Fw = Fw + 2.d0*x*y**2/3.d0
            if (k.eq.1) Fw = Fw + z*r(0)/3.d0 + 2.d0*x**2*Lz(0)
            if (k.eq.2) Fw = Fw - 2.d0*y*z*Lx(0) - 4.d0*y**2*Ky(0) + 4.d0*x*y/3.d0
            if (k.eq.3) Fw = Fw + x*r(0)/3.d0 - (y**2 + z**2)*Lx(0) 
         endif

      elseif (ij .eq. 22) then
         if     (mn .eq. 00) then
            Fw = -2.d0*y*Lx(k) + 4*z*Ky(k)
            if (k.eq.2) Fw = Fw - 2.d0*Lx(0)
            if (k.eq.3) Fw = Fw + 4.d0*Ky(0)
         elseif (mn .eq. 10) then
            Fw  = 2.d0*( y*z*Lz(k) - y*r(k) )
            if (k.eq.2) Fw = Fw + 2.d0*( z*Lz(0) - r(0) )
            if (k.eq.3) Fw = Fw + 2.d0*y*Lz(0)
         elseif (mn .eq. 01) then
            Fw = x*r(k) - 2.d0*x*z*Lz(k) - (y**2+z**2)*Lx(k)
            if (k.eq.1) Fw = Fw + r(0) - 2.d0*z*Lz(0)
            if (k.eq.2) Fw = Fw - 2.d0*y*Lx(0)
            if (k.eq.3) Fw = Fw - 2.d0*( x*Lz(0) + z*Lx(0) )
         elseif (mn .eq. 11) then
            Fw = (y**2 - x**2)*z*Lz(k) + (2.d0*x**2 + z**2)*r(k) - (4.d0/3.d0)*r3(k)
            if (k.eq.1) Fw = Fw - 2.d0*x*z*Lz(0) + 4.d0*x*r(0)
            if (k.eq.2) Fw = Fw + 2.d0*y*z*Lz(0)
            if (k.eq.3) Fw = Fw + (y**2 - x**2)*Lz(0) + 2.d0*z*r(0)
         endif
      elseif (ij .eq. 23) then
         if     (mn .eq. 00) then
            Fw = 2.d0*( x*Lz(k) + z*Lx(k) ) + 4.d0*y*Ky(k)
            if (k.eq.1) Fw = Fw + 2.d0*Lz(0)
            if (k.eq.2) Fw = Fw + 4.d0*Ky(0)
            if (k.eq.3) Fw = Fw + 2.d0*Lx(0)
         elseif (mn .eq. 10) then
            Fw  = (x**2 + y**2)*Lz(k) + z*r(k)
            if (k.eq.1) Fw = Fw + 2.d0*x*Lz(0)
            if (k.eq.2) Fw = Fw + 2.d0*y*Lz(0)
            if (k.eq.3) Fw = Fw + r(0)
         elseif (mn .eq. 01) then
            Fw = 2.d0*( y**2*Ky(k) - x**2*Kx(k) - x*z*Ly(k) - z**2*Kz(k) )
            if (k.eq.0) Fw = Fw + x*y
            if (k.eq.1) Fw = Fw - 4.d0*x*Kx(0) - 2.d0*z*Ly(0)
            if (k.eq.2) Fw = Fw + 4.d0*y*Ky(0)
            if (k.eq.3) Fw = Fw - 4.d0*z*Kz(0) - 2.d0*x*Ly(0)
         elseif (mn .eq. 11) then
            Fw = ( y*z*r(k) - (3.d0*x**2 + z**2)*z*Ly(k) + 2.d0*y**3*Lz(k) - 4.d0*x**3*Kx(k) )/3.d0
            if (k.eq.0) Fw = Fw + 2.d0*y*x**2/3.d0
            if (k.eq.1) Fw = Fw - 2.d0*x*z*Ly(0) - 4.d0*x**2*Kx(0) + 4.d0*x*y/3.d0
            if (k.eq.2) Fw = Fw + z*r(0)/3.d0 + 2.d0*y**2*Lz(0)
            if (k.eq.3) Fw = Fw + y*r(0)/3.d0 - (x**2 + z**2)*Ly(0) 
         endif

      elseif (ij .eq. 33) then
         if     (mn .eq. 00) then
            Fw = -2.d0*( x*Ly(k) + y*Lx(k) ) - 4.d0*z*Kz(k)
            if (k.eq.1) Fw = Fw - 2.d0*Ly(0)
            if (k.eq.2) Fw = Fw - 2.d0*Lx(0)
            if (k.eq.3) Fw = Fw - 4.d0*Kz(0)
         elseif (mn .eq. 10) then
            Fw = -(x**2 + z**2)*Ly(k) - y*r(k)
            if (k.eq.1) Fw = Fw - 2.d0*x*Ly(0)
            if (k.eq.2) Fw = Fw - r(0)
            if (k.eq.3) Fw = Fw - 2.d0*z*Ly(0)
         elseif (mn .eq. 01) then
            Fw = -(y**2 + z**2)*Lx(k) - x*r(k)
            if (k.eq.1) Fw = Fw - r(0)
            if (k.eq.2) Fw = Fw - 2.d0*y*Lx(0)
            if (k.eq.3) Fw = Fw - 2.d0*z*Lx(0) 
         elseif (mn .eq. 11) then
            Fw = -2.d0*r3(k)/3.d0
         endif

      else
         Fw = 1.d100
         write(*,*) 'WARNING: in the Fw function the numbers did not match any case'
      endif
   endif

   end associate
end function Fw

!------------------------------------------------------------------------------------------------------------
! DSYEVC3 -- Numerical diagonalization of 3x3 real symmetric matrices
! Copyright (C) 2006  Joachim Kopp
!------------------------------------------------------------------------------------------------------------
! This subroutine "dsyevc3" is free software; you can redistribute it and/or modify it under the terms of
! the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1
! of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
! the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General
! Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License along with this library; if not,
! write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
!------------------------------------------------------------------------------------------------------------

subroutine dsyevc3(a, v)
!--purpose: Calculates the eigenvalues of a symmetric 3x3 matrix A using Cardano's analytical algorithm.
!           Only the diagonal and upper triangular parts of A are accessed. The access is read-only.
   implicit none
!--subroutine parameters:
   real(kind=8) :: a(3,3)               !   A: The symmetric input matrix
   real(kind=8) :: v(3)                 !   V: Storage buffer for eigenvalues
!--local variables
   real(kind=8), parameter :: sqrt3 = 1.73205080756887729352744634151d0 
   real(kind=8) :: m, c1, c0
   real(kind=8) :: de, dd, ee, ff
   real(kind=8) :: p, sqrtp, q, c, s, phi, tmp
  
   ! Determine coefficients of characteristic poynomial. We write
   !       | A   D   F  |
   !  A =  | D*  B   E  |
   !       | F*  E*  C  |

   de    = a(1,2) * a(2,3)
   dd    = a(1,2)**2
   ee    = a(2,3)**2
   ff    = a(1,3)**2
   m     = a(1,1) + a(2,2) + a(3,3)
   c1    = ( a(1,1)*a(2,2) + a(1,1)*a(3,3) + a(2,2)*a(3,3) ) - (dd + ee + ff)
   c0    = a(3,3)*dd + a(1,1)*ee + a(2,2)*ff - a(1,1)*a(2,2)*a(3,3) - 2.0d0 * a(1,3)*de

   p     = m**2 - 3.0d0 * c1
   q     = m*(p - (3.0d0/2.0d0)*c1) - (27.0d0/2.0d0)*c0
   sqrtp = sqrt(abs(p))

   phi   = 27.0d0 * ( 0.25d0 * c1**2 * (p - c1) + c0 * (q + (27.0d0/4.0d0)*c0) )
   phi   = (1.0d0/3.0d0) * atan2(sqrt(abs(phi)), q)

   c     = sqrtp * cos(phi)
   s     = (1.0d0/sqrt3) * sqrtp * sin(phi)

   v(2) = (1.0d0/3.0d0) * (m - c)
   v(3) = v(2) + s
   v(1) = v(2) + c
   v(2) = v(2) - s

   ! sort the eigenvalues in v such that v(3) <= v(2) <= v(1)

   if (v(1).lt.v(2)) then
      tmp  = v(2)
      v(2) = v(1)
      v(1) = tmp
   endif
   if (v(2).lt.v(3)) then
      tmp  = v(3)
      v(3) = v(2)
      v(2) = tmp
   endif
   if (v(1).lt.v(2)) then
      tmp  = v(2)
      v(2) = v(1)
      v(1) = tmp
   endif

end subroutine dsyevc3

!------------------------------------------------------------------------------------------------------------

subroutine unique_positive (nam, nval, rinput, ixlist, lchanged)
!--purpose: Filter real input values, return unique positive integer values
   implicit none
!--subroutine arguments:
   character(len=*)            :: nam
   integer,      intent(inout) :: nval
   real(kind=8), intent(in)    :: rinput(nval)
   integer,      intent(out)   :: ixlist(nval)
   logical,      intent(out)   :: lchanged
!--local variables:
   integer i, nout
   logical lwarned

   ! sort in ascending order

   call bubsrt_real(nval, rinput)

   ! check input values for non-positives

   lwarned = .false.
   do i = 1, nval
      if (nint(rinput(i)).le.0 .and. .not.lwarned) then
         call write_log(' Input (subsurf): WARNING. Non-positive values for '//nam//' are ignored.')
         lwarned = .true.
      endif
   enddo

   ! check input values for duplicates

   lwarned = .false.
   do i = 2, nval
      if (nint(rinput(i)).eq.nint(rinput(i-1)) .and. .not.lwarned) then
         call write_log(' Input (subsurf): WARNING. Duplicate values for '//nam//' are ignored.')
         lwarned = .true.
      endif
   enddo

   ! copy unique positive values to ixlist

   nout = 0
   if (nint(rinput(1)).ge.1) then
      nout = nout + 1
      ixlist(nout) = nint(rinput(1))
   endif
   do i = 2, nval
      if (nint(rinput(i)).ge.1 .and. nint(rinput(i)).gt.nint(rinput(i-1))) then
         nout = nout + 1
         ixlist(nout) = nint(rinput(i))
      endif
   enddo

   if (nout.le.0) then
      call write_log('                           There are no input values remaining.')
   elseif (nout.lt.nval) then
      write(bufout,'(a,i4,a)') '                           The number of remaining values is',nout,'.'
      call write_log(1, bufout)
   endif

   ! return number of accepted values

   lchanged = (nout.ne.nval)
   nval = nout

end subroutine unique_positive

!------------------------------------------------------------------------------------------------------------

subroutine bubsrt_int (n, iarr)
!--purpose: Sort all elements of a vector X in ascending order; Uses the Bubble sort algorithm
   implicit none
   integer      :: n, iarr(n)
   integer      :: i, i_last, i_tmp, ival

   if (n.le.1) return

   ! invariant: i_last = index of last element that can be out of place

   i_last = n
   do while(i_last.gt.1)
      i_tmp = -1

      ! i_tmp = index of last modification

      do i = 1, i_last-1
         if (iarr(i).gt.iarr(i+1)) then
            i_tmp = i

            ! swap elements:

            ival      = iarr(i+1)
            iarr(i+1) = iarr(i)
            iarr(i)   = ival
         endif
      enddo

      ! if last modification is swap (i,i+1) --> element i can be out of place.

      i_last = i_tmp
   enddo

end subroutine bubsrt_int

!------------------------------------------------------------------------------------------------------------

subroutine bubsrt_real (n, rarr)
!--purpose: Sort all elements of a vector X in ascending order; Uses the Bubble sort algorithm
   implicit none
   integer      :: n
   real(kind=8) :: rarr(n)
   integer      :: i, i_last, i_tmp
   real(kind=8) :: rval

   if (n.le.1) return

   ! Invariant: I_last = index of last element that can be out of place

   i_last = n
   do while(i_last.gt.1)

      ! i_tmp = index of last modification

      i_tmp = -1
      do i = 1, i_last-1
         if (rarr(i).gt.rarr(i+1)) then
            i_tmp = i

            ! swap elements:

            rval      = rarr(i+1)
            rarr(i+1) = rarr(i)
            rarr(i)   = rval
         endif
      enddo

      ! if last modification is swap (i,i+1) --> element i can be out of place.

      i_last = i_tmp
   enddo

end subroutine bubsrt_real

!============================================================================================================

end module m_subsurf
