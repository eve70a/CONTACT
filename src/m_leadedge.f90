!------------------------------------------------------------------------------------------------------------
! m_leadedge - administration and computation of 'leading edge correction'
!
! Copyright 2008-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------

module m_leadedge

use m_hierarch_data
use m_aijpj

implicit none
private

public sxbnd
public subnd

contains

!------------------------------------------------------------------------------------------------------------

   subroutine sxbnd(ic, solv, cgrid, igs, chi, dq, ledg)
!--purpose: estimate the positions of leading edge boundaries in all rows my of the potential contact area.
!           See m_hierarch_data.f90 for an overview of the data-structures;
!           See vollebregt2009a-cm2009.pdf for an overview of the discretization.
      implicit none
!--subroutine arguments:
      type(t_ic)               :: ic
      type(t_solvers)          :: solv
      type(t_grid),     target :: cgrid
      type(t_eldiv)            :: igs
      type(t_leadedge), target :: ledg
      real(kind=8)             :: chi, dq
!--local variables:
      real(kind=8), parameter :: pi     = 4d0*atan(1d0)
      integer          :: idebug
      integer          :: ix, ixsta, ixend, iy, ii, j, nwarn, jwarn, nposiy
      integer, pointer :: mx, my, ixinc, npos
      real(kind=8)     :: fxdfac
      logical          :: use_ledg, is_roll, is_ssrol

      mx => cgrid%nx
      my => cgrid%ny
      ixinc => ledg%ixinc
      npos  => ledg%npos

      idebug = 0
      if (ic%flow.ge.9) idebug = 2

      is_roll  = ic%tang.eq.2 .or. ic%tang.eq.3
      is_ssrol = ic%tang.eq.3

      ! set start/increment/end for elements per grid row, increasing in rolling direction

      if (is_roll .and. abs(chi-pi).lt.0.01d0) then
         ixsta = mx
         ixinc = -1
         ixend =  1
      else
         ixsta =  1
         ixinc =  1
         ixend = mx
      endif

      !------------------------------------------------------------------------------------------------------
      !     1) determine the number of leading edge positions for all rows of the potential contact
      !------------------------------------------------------------------------------------------------------

      call reallocate_arr(ledg%jbnd, my+1)

      ledg%jbnd(1) = 1
      do iy = 1, my
         npos = 0

         ! count transitions from Adhes,Slip,Plast to Exter in "1..mx-1" (0deg)

         do ix = ixsta, ixend-ixinc, ixinc
            ii = ix + (iy-1)*mx
            if (igs%el(ii).ge.Adhes .and. igs%el(ii+ixinc).le.Exter) npos = npos + 1
         enddo

         ! account for interior element "mx" (0deg)

         ii = ixend + (iy-1)*mx
         if (igs%el(ii).ge.Adhes) npos = npos + 1

         ! store number of transitions in row iy

         ledg%jbnd(iy+1) = ledg%jbnd(iy) + npos
      enddo
      if (idebug.ge.3)                                                                                 &
         write(*,'(a,i4,a,99(i4,a))') ' In total there are',ledg%jbnd(my+1)-1,' transitions C->E:',    &
            (ledg%jbnd(iy+1)-ledg%jbnd(iy),',', iy=1,my)

      !------------------------------------------------------------------------------------------------------
      !     2) allocate ixbnd at the appropriate size, fill ixbnd
      !------------------------------------------------------------------------------------------------------

      npos = ledg%jbnd(my+1)
      call reallocate_arr(ledg%ixbnd, npos)

      do iy = 1, my
         nposiy = ledg%jbnd(iy+1) - ledg%jbnd(iy)
         if (nposiy.gt.0) then
            j = ledg%jbnd(iy)
            do ix = ixsta, ixend-ixinc, ixinc
               ii = ix + (iy-1)*mx
               if (igs%el(ii).ge.Adhes .and. igs%el(ii+ixinc).le.Exter) then
                  ledg%ixbnd(j) = ix
                  j = j + 1
               endif
            enddo
            ii = ixend + (iy-1)*mx
            if (igs%el(ii).ge.Adhes) then
               ledg%ixbnd(j) = ixend
               j = j + 1
            endif
         endif
         if (idebug.ge.5) write(*,'(2(a,i4),a,5(i4,a))') ' row iy=',iy, ': there are ',nposiy,         &
            ' transitions C->E at ix=', (ledg%ixbnd(j),',', j=ledg%jbnd(iy),ledg%jbnd(iy+1)-1)
      enddo

      !------------------------------------------------------------------------------------------------------
      !     4) estimate positions of leading edge
      !------------------------------------------------------------------------------------------------------

      ! Do not use leading edge correction together with FASTSIM or SteadyGS
      ! Else, use method 1 (fixed "1.0 * dx") for leading edge position.

      use_ledg = .true.  .and. solv%solver_eff.ne.isolv_fastsm                                         &
                         .and. solv%solver_eff.ne.isolv_stdygs

      if (.not.use_ledg) then
         fxdfac = 2.0d0
      else
         fxdfac = 1.0d0
      endif

      ! reallocate arrays at the appropriate size

      npos = ledg%jbnd(my+1)
      call reallocate_arr(ledg%xbnd, npos)
      call reallocate_arr(ledg%facdx, npos)

      ! for all rows iy do

      do iy = 1, my

         ! for all leading edge positions in row iy do

         do j = ledg%jbnd(iy), ledg%jbnd(iy+1)-1
            ix = ledg%ixbnd(j)
            ii = ix + (iy-1)*mx

            ! estimate the location Xbnd for the leading edge between elements ix and "ix+1" (0deg) at row iy
            ! Note: facdx = 0 for center of ix, facdx = 1 at center of ix+1 or ix-1 (180deg)

            ! method 1: using fixed position in between of the element centers

            ledg%facdx(j) = fxdfac

            ! compute actual position of boundary and store at appropriate position in xbnd

            ledg%xbnd(j) = cgrid%x(ii) + ixinc *ledg%facdx(j) * cgrid%dx

            if (idebug.ge.5) then
               write(bufout,3011) iy, ix, ledg%xbnd(j), ledg%facdx(j)
               call write_log(1, bufout)
 3011          format(' lead.edge in row',i4,' el.',i4,': estimated x-pos',f9.5,' (',f6.3,'*dx)')
            endif

         enddo
      enddo

      !------------------------------------------------------------------------------------------------------
      !     5) for all elements, determine the fraction of the timestep in the contact area
      !------------------------------------------------------------------------------------------------------

      call if1_new(ledg%ii2j, cgrid)
      call gf3_new(ledg%facdt, 'stang:ledg%facdt', cgrid)
      nwarn = 0
      jwarn = 0

      do iy = 1, my
         j     = ledg%jbnd(iy)
         do ix = ixsta, ixend, ixinc

            ! advance j when necessary, make it point to the next leading edge element ix0 in row iy,
            ! or if there is no such element, to the first leading edge position for row iy+1.

            do while(j.lt.ledg%jbnd(iy+1) .and.                                                        &
                     ((ixinc.gt.0 .and. ix.gt.ledg%ixbnd(j)) .or. (ixinc.lt.0 .and. ix.lt.ledg%ixbnd(j))))
               j = j + 1
            enddo

            ii = ix + (iy-1)*mx

            ! determine fraction of time-step in the contact area and the corresponding leading edge point j

            if (.not.is_roll) then

               ! shift transient: no leading edge correction, facdt=1

               ledg%ii2j%val(ii) = 0
               ledg%facdt%vt(ii) = 1d0

            elseif (igs%el(ii).le.Exter) then

               ! exterior elements: no equations/no correction

               ledg%ii2j%val(ii) = 0
               ledg%facdt%vt(ii) = 0d0

            elseif ((ixinc.gt.0 .and. ix+2.gt.mx) .or. (ixinc.lt.0 .and. ix-2.le.0)) then

               ! not enough points in potential contact area:
               !    treat as whole timestep in contact area, no correction

               ledg%ii2j%val(ii) = 0
               ledg%facdt%vt(ii) = 1d0
               if (jwarn.eq.j) then
               !    already counted/warned about lead.edge position j
               else
                  nwarn = nwarn + 1
                  jwarn = j
                  if (ic%flow.ge.7) then
                     write(bufout,4001) iy
                     call write_log(1, bufout)
 4001                format(' TANG: Warning: pot.con. too small; cannot estimate Ubnd at row iy=',i4)
                  endif
               endif

            else

               ! interior elements:

               ! the leading edge is at (ix0+/-facdx)*dx, the element center is at ix*dx,
               !    the traversed distance per time-step is dq

               ledg%facdt%vt(ii) = min(1d0,ixinc*(ledg%xbnd(j) - cgrid%x(ii))/dq)

               ! record whether element ii is close to leading edge point j or not

               if (ledg%facdt%vt(ii).lt.0.9999d0) then
                  ledg%ii2j%val(ii) = j
               else
                  ledg%ii2j%val(ii) = 0
               endif

               if (idebug.ge.5 .and. ledg%ii2j%val(ii).gt.0) then
                  write(bufout,4011) ix, iy, ledg%facdt%vt(ii)
                  call write_log(1, bufout)
 4011             format(' center of element',i4,' of row',i4,' is',f6.3,'*dt in interior')
               endif

            endif
         enddo
      enddo

      if (use_ledg .and. nwarn.gt.0 .and. ic%flow.ge.3) then
         write(bufout,4021) nwarn
         call write_log(1, bufout)
 4021    format(' TANG: Warning: pot.con. too small; cannot estimate Ubnd at',i4,' rows.')
      endif

   end subroutine sxbnd

!------------------------------------------------------------------------------------------------------------

   subroutine subnd(cgrid, p, c, ledg)
!--purpose: estimate the tangential displacement difference UBnd at the leading edge of the contact area
      implicit none
!--subroutine arguments:
      type(t_grid),     target :: cgrid
      type(t_gridfnc3)         :: p
      type(t_inflcf)           :: c
      type(t_leadedge), target :: ledg
!--local variables:
      integer            :: iy, ix0, ii, ik, j, npos
      integer, pointer   :: mx, my, ixinc
      real(kind=8)       :: facdx, uii0, uii1
!--functions called:
      !x! real(kind=8) :: AijPj

      mx => cgrid%nx
      my => cgrid%ny
      ixinc => ledg%ixinc

      ! allocate ubnd at the appropriate size

      npos = ledg%jbnd(my+1)
      call reallocate_arr(ledg%ubnd, npos, 3)

      ! For each leading edge position in each row iy,
      !    estimate the displacement difference u at the edge of the contact area.

      do iy = 1, my
         do j = ledg%jbnd(iy), ledg%jbnd(iy+1)-1

            ! get index of last interior grid cell of current segment

            ix0 = ledg%ixbnd(j)

            if ((ixinc.gt.0 .and. ix0+2.gt.mx) .or. (ixinc.lt.0 .and. ix0-2.le.0)) then
               ledg%ubnd(j,1:2) = 0d0
            else
               ii = ix0 + (iy-1)*mx
               do ik = 1, 2

                  uii0 = AijPj(ii        , ik, p, jkALL, c)
                  uii1 = AijPj(ii+1*ixinc, ik, p, jkALL, c)
                  facdx = ledg%facdx(j)

                  ! method 2: Interpolating u between ix0 and ix0+1:
                  !    u(xbnd) =  u(ix0+facdx)
                  !           ~=~ u(i+0) + facdx*(u(i+1) - u(i+0))

                  ledg%ubnd(j,ik) =   (1d0-facdx) * uii0 + facdx  * uii1

                  if (.false. .and. ik.eq.1) then
                     write(bufout,123) 'ix0=',ix0,': ux0=',uii0,                                       &
                        ', ux1=',uii1,', f=',facdx,', ubnd=',ledg%ubnd(j,ik)
                     call write_log(1, bufout)
 123                 format(1x,a,i3,2(a,f10.7),a,f5.2,a,f10.7)
                  endif

               enddo
            endif
         enddo
         ! if (iy.eq.1) write(*,'(a,i4,a,2f9.5)') '  Row',iy,': UBnd=',ledg%ubnd(j,1), ledg%ubnd(j,2)
      enddo

   end subroutine subnd

!------------------------------------------------------------------------------------------------------------

end module m_leadedge
