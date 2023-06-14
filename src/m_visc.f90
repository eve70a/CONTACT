!------------------------------------------------------------------------------------------------------------
! m_visc - computation of influence coefficients for elastic and visco-elastic half-space problems
!
! Copyright 1979-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------

module m_visc

use m_hierarch_data

implicit none
private

   private mater_set_visc
   public  sgencr
   public  azz_one_element
   public  elascf_pcwcns
   public  elascf_bilin
   public  visccf

   private Hklm
   private Gklm
   private iFklm
   private Fklm
   
   private dnumin
   private fnn
   private fxx
   private fxy
   private fxn
   private fyn
   private dqagi
   private dqk15i
   private dqsort
   private dqext
   private dqmaco

   real(kind=8), parameter :: epsco = 1.0d-15

   type t_xyr
      real(kind=8), dimension(:),   allocatable :: x
      real(kind=8), dimension(:),   allocatable :: y
      real(kind=8), dimension(:,:), allocatable :: r
      real(kind=8), dimension(:,:), allocatable :: xlog
      real(kind=8), dimension(:,:), allocatable :: ylog
      real(kind=8), dimension(:,:), allocatable :: rlog
      real(kind=8), dimension(:,:), allocatable :: xlogy
      real(kind=8), dimension(:,:), allocatable :: ylogx
      real(kind=8), dimension(:,:), allocatable :: xtanyx
      real(kind=8), dimension(:,:), allocatable :: ytanxy
   end type t_xyr

contains

!------------------------------------------------------------------------------------------------------------

   subroutine mater_set_visc (ic, mater, fdg)
!--purpose: Set the combined material properties for viscoelastic computations
      implicit none
!--subroutine arguments:
      type(t_ic)               :: ic
      type(t_material), target :: mater
      real(kind=8)             :: fdg(2)
!--local variables:
      integer,      parameter  :: idebug = 0
      real(kind=8), parameter  :: epsco  = 1d-15
      integer                  :: i
      real(kind=8)             :: akt(2), dg(2), gat(2)

      associate(ak    => mater%ak,   ga    => mater%ga,    nu    => mater%nu,                           &
                akv   => mater%akv,  gav   => mater%gav,   nuv   => mater%nuv,                          &
                gg    => mater%gg,   poiss => mater%poiss, fg    => mater%fg,   vt    => mater%vt)

      ! Material parameters: elastic and viscoelastic parameters are combined

      if (ic%mater.ne.1) then
         akv = ak
         nuv = nu
         gav = ga
      else
         do i = 1, 2
            dg(i)  =   1d0 / gg(i)      ! initial glassy compliance
            fdg(i) = fg(i) / gg(i)      ! secondary compliance
            gat(i) = dg(i) + fdg(i)     ! total rubbery compliance
            akt(i) = (1d0 - 2d0*poiss(i))*dg(i)
            if (idebug.ge.2) then
 111           format(a,i2,a,f9.2,a,f5.2,a,f5.2,a,/, a,f9.6,a,f9.6,a,f9.6,a,f9.6)
               write(bufout,111) ' body i=',i,': G=',gg(i),', nu=',poiss(i),', F=',fg(i),':',         &
                    '          1/G=',dg(i),', F/G=',fdg(i),', (1+F)/G=',gat(i),', (1-2nu)/G=',akt(i)
               call write_log(2, bufout)
            endif
         enddo

         ! Note: gav is set but not used; ga is used instead.

         if (vt(1).le.epsco .and. vt(2).gt.epsco) then
            if (idebug.ge.1) call write_log(' case 1: vt(1)=0 (rubbery), vt(2)>0 (visc)')
            gav = 2d0 / (gat(1) + dg(2))
            nuv = gav/2d0 * (gat(1)*poiss(1) + dg(2)*poiss(2))
            akv = gav/4d0 * (akt(1)*(1d0+fg(1)) - akt(2))
         elseif (vt(1).gt.epsco .and. vt(2).le.epsco) then
            if (idebug.ge.1) call write_log(' case 2: vt(1)>0 (visc), vt(2)=0 (rubbery)')
            gav = 2d0 / (dg(1) + gat(2))
            nuv = gav/2d0 * (dg(1)*poiss(1) + gat(2)*poiss(2))
            akv = gav/4d0 * (akt(1) - akt(2)*(1d0+fg(2)))
         elseif (vt(1).le.epsco .and. vt(2).le.epsco) then
            if (idebug.ge.1) call write_log(' case 3: vt(1)=0 (rubbery), vt(2)=0 (rubbery)')
            gav = 2d0 / (gat(1) + gat(2))
            nuv = gav/2d0 * (gat(1)*poiss(1) + gat(2)*poiss(2))
            akv = gav/4d0 * (akt(1)*(1d0+fg(1)) - akt(2)*(1d0+fg(2)))
         else ! if (vt(1).gt.epsco .and. vt(2).gt.epsco) then
            if (idebug.ge.1) call write_log(' case 4: vt(1)>0 (visc), vt(2)>0 (visc)')
            akv = ak
            nuv = nu
            gav = ga
         endif
      endif

      end associate

   end subroutine mater_set_visc

!------------------------------------------------------------------------------------------------------------

   subroutine sgencr (ic, cgrid, mater, influ, fric, kin)
!--purpose: Construct the influence coefficients with the aid of subroutine visccf
      implicit none
!--subroutine arguments:
      type(t_ic)               :: ic
      type(t_grid),     target :: cgrid
      type(t_material), target :: mater
      type(t_friclaw),  target :: fric
      type(t_kincns),   target :: kin
      type(t_influe),   target :: influ
!--constants/parameters:
      integer,      parameter  :: idebug = 0
      real(kind=8), parameter  :: epsco  = 1d-15
      logical,      parameter  :: use_dq_scaling = .false.
!--pointers to global data items:
      real(kind=8), dimension(:,:,:,:), pointer :: cs, cv, csv
!--local variables:
      integer      :: ix, iy, ik, jk
      real(kind=8) :: cc, facdq, facdqi, sc, xshft, yshft, fdg(2)
      logical      :: is_roll

      call timer_start(itimer_sgencr)

      associate(npot  => cgrid%ntot,  mx    => cgrid%nx,  my    => cgrid%ny,   dx    => cgrid%dx,       &
                dy    => cgrid%dy,    chi   => kin%chi,   dq    => kin%dq,                              &
                ak    => mater%ak,    ga    => mater%ga,  nu    => mater%nu,   akv   => mater%akv,      &
                gav   => mater%gav,   nuv   => mater%nuv, gg    => mater%gg,   poiss => mater%poiss,    &
                fg    => mater%fg,    vt    => mater%vt)

      is_roll   = ic%tang.eq.2 .or. ic%tang.eq.3

      ! make sure that combined material parameters are up-to-date

      call combin_mater(mater)

      ! (re-)allocate influence coefficient arrays at the appropriate size
      ! Note: the arrays can be of size (-mx+1:mx-1,-my+1:my-1). One additional grid row and column are
      !       added to the left/bottom to obtain more favourable sizes for Fast Fourier transform.
      !       Elements (-mx,:) and (:,-my) are referenced but do not affect the calculations.

      call inflcf_new(influ%cs,  0, cgrid)
      call inflcf_new(influ%cv,  0, cgrid)
      call inflcf_new(influ%csv, 0, cgrid)
      call inflcf_new(influ%ms,  0, cgrid)

      ! Copy relevant material parameters to influence matrices

      call influe_mater(influ, mater)

      cs    => influ%cs%cf
      cv    => influ%cv%cf
      csv   => influ%csv%cf

      if (is_roll) then
         cc = cos(chi)
         sc = sin(chi)
      else
         cc = 0d0
         sc = 0d0
      endif

      ! Material parameters: elastic and viscoelastic parameters are combined

      call mater_set_visc (ic, mater, fdg)

      ! Enable short-cut calculation if the difference parameter is
      !    small, i.e. normal-tangential coupling is zero or negligible

      if (abs(akv).lt.1d-6) then
         influ%cs%nt_cpl = .false.
         influ%cv%nt_cpl = .false.
         influ%csv%nt_cpl = .false.
      endif

      ! Determine multiplication factor for dq for computation of time derivative (effective dq = 0.6 dx)

      if (is_roll) then
         if (use_dq_scaling .and. dq.lt.1.501d0*dx) then
            if (fric%is_veldep()) then
               facdq = 1.00d0*dx / dq
            else
               facdq = 0.60d0*dx / dq
            endif
            if (ic%flow.ge.6) then
               write(bufout,'(a,f8.3,a,f5.2,a)') ' Using facdq=',facdq,', i.e. dq* =',facdq*dq/dx,' dx'
               call write_log(1, bufout)
            endif
         else
            facdq = 1d0
         endif
      else
         facdq = 1d0
      endif
      facdqi = 1d0 / facdq

      ! 1) Form the array of influence numbers "cs".

      xshft = 0d0
      yshft = 0d0

      ! Use constant basis functions if C<>3, else use bilinear basis functions

      if (mater%gencr_eff.ne.3) then
         call elascf_pcwcns(akv, nuv, mx, my, dx, dy, xshft, yshft, influ%cs)
      else
         call elascf_bilin (akv, nuv, mx, my, dx, dy, xshft, yshft, influ%cs)
      endif
      if (ic%mater.eq.1) then
         call visccf (akv, nuv, ga, poiss, fdg, fg, vt, mx, my, dx, dy, xshft, yshft, influ%cs, idebug)
      endif

      ! 2) Form the array of influence numbers "cv".
      !    The coefficients cf(ix-jx,iy-jy,:,:) say what the deformation is in element (ix,iy) of the
      !    grid of the current time, when a unit pressure is located in the element (jx,jy) of the grid
      !    of the previous time.

      ! 3) Form the array of influence numbers "csv".
      !    These coefficients say what the increment of the deformation is in element (ix,iy) between
      !    the previous and current time when a unit pressure is located in the element (jx,jy) of the
      !    current grid and of the grid of the previous time.

      if (.not.is_roll) then

         ! in shifts, the coordinate system is world-fixed, the current grid equals the previous grid,
         ! therefore the displacements can be computed using cs.

         do jk = 1, 3
            do ik = 1, 3
               do iy = -my, my-1
                  do ix = -mx, mx-1
                     cv(ix,iy,ik,jk)  = cs(ix,iy,ik,jk)
                     csv(ix,iy,ik,jk) = cs(ix,iy,ik,jk)
                  enddo
               enddo
            enddo
         enddo

      else

         ! in rolling, the coordinate system moves with the contact patch with direction chi and
         ! distance traversed per step dq.

         ! For chi=0 the grid moves to the right. a point (xx,yy) of the current grid therefore
         ! lies at (xx+cc*dq,yy+sc*dq) in the previous grid where the tractions are thought to be
         ! centered around the origin.

         xshft = cc * facdq * dq
         yshft = sc * facdq * dq

         ! If C = 3, use bilinear basis functions to determine the influence coefficients,
         ! else use piecewise constant basis functions.

         if (mater%gencr_eff.ne.3) then
            call elascf_pcwcns(akv, nuv, mx, my, dx, dy, xshft, yshft, influ%cv)
         else
            call elascf_bilin (akv, nuv, mx, my, dx, dy, xshft, yshft, influ%cv)
         endif
         if (ic%mater.eq.1) then
            call visccf (akv, nuv, ga, poiss, fdg, fg, vt, mx, my, dx, dy, xshft, yshft, influ%cv, idebug)
         endif

         ! Copy normal part from influence matrix cs

         do jk = 1, 3
            do iy = -my, my-1
               do ix = -mx, mx-1
                  cv (ix,iy,3,jk) = cs (ix,iy,3,jk)
               enddo
            enddo
         enddo

         ! Compute "scaled derivative" in csv
         ! Adjust cv according to "scaled derivative"

         do jk = 1, 3
            do iy = -my, my-1
               do ix = -mx, mx-1
                  csv(ix,iy,1,jk) = facdqi * (cs(ix,iy,1,jk) - cv(ix,iy,1,jk))
                  csv(ix,iy,2,jk) = facdqi * (cs(ix,iy,2,jk) - cv(ix,iy,2,jk))
                  csv(ix,iy,3,jk) =           cs(ix,iy,3,jk)

                  cv(ix,iy,1,jk) = cs(ix,iy,1,jk) - csv(ix,iy,1,jk)
                  cv(ix,iy,2,jk) = cs(ix,iy,2,jk) - csv(ix,iy,2,jk)
               enddo
            enddo
         enddo
      endif

      ! print influence coefficients to output-file when requested

      if (ic%nmdbg.ge.7) then
         if (.false.) then
            write(lout,*) 'Printing influence coefficients csv (steady)'
            call inflcf_print (influ%csv, lout)
         else
            write(lout,*) 'Printing influence coefficients cs (instat)'
            call inflcf_print (influ%cs, lout)
         endif
      endif

      end associate
      call timer_stop(itimer_sgencr)

   end subroutine sgencr

!------------------------------------------------------------------------------------------------------------

   subroutine azz_one_element(ic, mater, dx, dy, azz)
!--purpose: Compute the influence coefficient Azz(0,0) with the aid of subroutine visccf
      implicit none
!--subroutine arguments:
      type(t_ic)                    :: ic
      type(t_material), target      :: mater
      real(kind=8),     intent(in)  :: dx, dy
      real(kind=8),     intent(out) :: azz
!--constants/parameters:
      integer,      parameter  :: idebug = 0
!--local variables:
      integer        :: mx, my
      real(kind=8)   :: xshft, yshft, fdg(2)
      type(t_grid)   :: cgrid
      type(t_inflcf) :: tmp_cs

      ! make sure that combined material parameters are up-to-date

      call combin_mater(mater)
      call mater_set_visc (ic, mater, fdg)

      ! compute influence coefficients using appropriate settings

      associate(gg    => mater%gg,  poiss => mater%poiss, ga    => mater%ga,                    &
                akv   => mater%akv, nuv   => mater%nuv,   fg    => mater%fg,  vt    => mater%vt)

      call grid_set_dimens(cgrid, 1, 1)
      call inflcf_new(tmp_cs,  0, cgrid)

      mx = 1
      my = 1
      xshft = 0d0
      yshft = 0d0

      if (mater%gencr_eff.ne.3) then
         call elascf_pcwcns(akv, nuv, mx, my, dx, dy, xshft, yshft, tmp_cs)
      else
         call elascf_bilin (akv, nuv, mx, my, dx, dy, xshft, yshft, tmp_cs)
      endif

      if (ic%mater.eq.1) then
         call visccf (akv, nuv, ga, poiss, fdg, fg, vt, mx, my, dx, dy, xshft, yshft, tmp_cs, idebug)
      endif
      end associate

      azz = tmp_cs%cf(0, 0, ikZDIR, ikZDIR)

   end subroutine azz_one_element

!------------------------------------------------------------------------------------------------------------

   subroutine elascf_pcwcns(akv, nuv, mx, my, dx, dy, xshft, yshft, cs)
!--purpose: Calculates the influence coefficients according to Boussinesq and Cerruti. See memo m09002.
!
!           This routine calculates the deformation at all grid points (xx,yy) of the potential contact
!           area due to a unit pressure in the element with center (0, 0). The grid may be shifted by
!           an offset (xshft,yshft).
!
!           The coefficients cf(ix-jx,iy-jy,:,:) say what the deformation is in element (ix,iy), when a
!           unit pressure is located in element (jx,jy).
      implicit none
!--subroutine arguments:
      integer        :: mx, my
      real(kind=8)   :: akv, nuv, dx, dy, xshft, yshft
      type(t_inflcf) :: cs
!--local variables:
      integer      :: ix, iy, iy0, ik, jk
      real(kind=8) :: e1, e2, e3, tolx, tolx2
      real(kind=8), dimension(:,:), allocatable :: alr, atxy, atyx, r, x, y, xly2y1, ylx2x1

      allocate(r(-mx:mx,-my:my), x(-mx:mx,-my:my), y(-mx:mx,-my:my))

      e1 = (1d0-nuv) / pi
      e2 =  1d0      / pi
      e3 =      nuv  / pi

      tolx  = dx*1d-10
      tolx2 = dx*1d-20

      ! Initialize influence coefficients

      do jk = 1, 3
         do ik = 1, 3
            do iy = -my, my-1
               do ix = -mx, mx-1
                  cs%cf(ix,iy,ik,jk) = 0d0
               enddo
            enddo
         enddo
      enddo

      ! In 2D problems, the FFT preconditioner requires cs=0 at iy=-my

      iy0 = -my
      if (my.eq.1) iy0 = -my + 1 ! = 0

      ! compute coordinates of observation point (ix,iy)+shft relative to the four corners of element (0,0)

      ! Note: x and/or y may become zero or close to zero in the determination of cv, when the traversed
      !       distance per timestep contains a half grid width, i.e. xshft=0.5*dx, 1.5*dx, etc (chi=0).
      !       This yields r=|y| or r=|x| and r+y may be 0.

      do iy = iy0, my
         do ix = -mx, mx
            x(ix,iy) = real(ix)*dx + xshft - 0.5d0*dx
            y(ix,iy) = real(iy)*dy + yshft - 0.5d0*dy
            r(ix,iy) = dsqrt(x(ix,iy)**2 + y(ix,iy)**2)
         enddo
      enddo

      ! The extended calculation (short calculation was used when mx > 1/eps)

      ! for all possible combinations of (ix-jx, iy-jy) :

      do iy = iy0, my-1
         do ix = -mx, mx-1

            ! cf(2,1) == A^vx == [[  nu j4     ]],  j4 = - r

            cs%cf(ix,iy,2,1) = cs%cf(ix,iy,2,1) - e3 * r(ix  ,iy  )
            cs%cf(ix,iy,2,1) = cs%cf(ix,iy,2,1) + e3 * r(ix  ,iy+1)
            cs%cf(ix,iy,2,1) = cs%cf(ix,iy,2,1) + e3 * r(ix+1,iy  )
            cs%cf(ix,iy,2,1) = cs%cf(ix,iy,2,1) - e3 * r(ix+1,iy+1)
         enddo
      enddo

      if (abs(akv).ge.1d-6) then
         allocate(alr(-mx:mx,-my:my), atxy(-mx:mx,-my:my), atyx(-mx:mx,-my:my))
         do iy = iy0, my
            do ix = -mx, mx
               alr (ix,iy) = dlog(r(ix,iy))
               atyx(ix,iy) = datan(y(ix,iy) / (x(ix,iy)+tolx2))
               atxy(ix,iy) = datan(x(ix,iy) / (y(ix,iy)+tolx2))
            enddo
         enddo

         do iy = iy0, my-1
            do ix = -mx, mx-1

               ! cf(1,3) == A^uz == [[  -k j5     ]],  
               !             j5 = (y-eta) log r + (x-xi) atan((y-eta)/(x-xi))

               cs%cf(ix,iy,1,3) = cs%cf(ix,iy,1,3)                                                      &
                                - e2*akv * (y(ix  ,iy  )*alr(ix  ,iy  ) + x(ix  ,iy  )*atyx(ix  ,iy  ))
               cs%cf(ix,iy,1,3) = cs%cf(ix,iy,1,3)                                                      &
                                + e2*akv * (y(ix  ,iy+1)*alr(ix  ,iy+1) + x(ix  ,iy+1)*atyx(ix  ,iy+1))
               cs%cf(ix,iy,1,3) = cs%cf(ix,iy,1,3)                                                      &
                                + e2*akv * (y(ix+1,iy  )*alr(ix+1,iy  ) + x(ix+1,iy  )*atyx(ix+1,iy  ))
               cs%cf(ix,iy,1,3) = cs%cf(ix,iy,1,3)                                                      &
                                - e2*akv * (y(ix+1,iy+1)*alr(ix+1,iy+1) + x(ix+1,iy+1)*atyx(ix+1,iy+1))

               ! cf(2,3) == A^vz == [[  -k j6     ]],  
               !             j6 = (x-xi) log r + (y-eta) atan((x-xi)/(y-eta))

               cs%cf(ix,iy,2,3) = cs%cf(ix,iy,2,3)                                                      &
                                - e2*akv * (x(ix  ,iy  )*alr(ix  ,iy  ) + y(ix  ,iy  )*atxy(ix  ,iy  ))
               cs%cf(ix,iy,2,3) = cs%cf(ix,iy,2,3)                                                      &
                                + e2*akv * (x(ix  ,iy+1)*alr(ix  ,iy+1) + y(ix  ,iy+1)*atxy(ix  ,iy+1))
               cs%cf(ix,iy,2,3) = cs%cf(ix,iy,2,3)                                                      &
                                + e2*akv * (x(ix+1,iy  )*alr(ix+1,iy  ) + y(ix+1,iy  )*atxy(ix+1,iy  ))
               cs%cf(ix,iy,2,3) = cs%cf(ix,iy,2,3)                                                      &
                                - e2*akv * (x(ix+1,iy+1)*alr(ix+1,iy+1) + y(ix+1,iy+1)*atxy(ix+1,iy+1))
            enddo
         enddo
         deallocate(alr, atxy, atyx)
      endif

      ! Determine logarithm terms, accounting for values close to zero

      allocate(xly2y1(-mx:mx,-my:my), ylx2x1(-mx:mx,-my:my))

      do iy = iy0, my-1
         do ix = -mx, mx
            if (min(abs(y(ix,iy)+r(ix,iy)), abs(y(ix,iy+1)+r(ix,iy+1))).lt.tolx) then
               xly2y1(ix,iy) = 0.d0
            else
               xly2y1(ix,iy) = x(ix,iy) * dlog( (y(ix,iy+1)+r(ix,iy+1)) / (y(ix,iy)+r(ix,iy)) )
            endif
         enddo
      enddo

      do iy = iy0, my
         do ix = -mx, mx-1
            if (min(abs(x(ix,iy)+r(ix,iy)), abs(x(ix+1,iy)+r(ix+1,iy))).lt.tolx) then
               ylx2x1(ix,iy) = 0.d0
            else
               ylx2x1(ix,iy) = y(ix,iy) * dlog( (x(ix+1,iy)+r(ix+1,iy)) / (x(ix,iy)+r(ix,iy)) )
            endif
         enddo
      enddo

      do iy = iy0, my-1
         do ix = -mx, mx-1

            ! cf(3,3) == A^wz == [[ (1-nu) j3  ]],      j3 = j1 + j2

            cs%cf(ix,iy,3,3) = cs%cf(ix,iy,3,3) -  e1*(xly2y1(ix  ,iy) +    ylx2x1(ix,iy  ))
            cs%cf(ix,iy,3,3) = cs%cf(ix,iy,3,3) +  e1*(xly2y1(ix+1,iy) +    ylx2x1(ix,iy+1))

            ! cf(1,1) == A^ux == [[ j1 + (1-nu) j2 ]],  j1 = (y-eta) log ((x-xi)+r)
            !                 == [[ j3 - nu j2 ]]

            cs%cf(ix,iy,1,1) = cs%cf(ix,iy,1,1) - (e1* xly2y1(ix  ,iy) + e2*ylx2x1(ix,iy  ))
            cs%cf(ix,iy,1,1) = cs%cf(ix,iy,1,1) + (e1* xly2y1(ix+1,iy) + e2*ylx2x1(ix,iy+1))

            ! cf(2,2) == A^vy == [[ j2 + (1-nu) j1 ]],  j2 = (x-xi) log ((y-eta)+r)
            !                 == [[ j3 - nu j1 ]]

            cs%cf(ix,iy,2,2) = cs%cf(ix,iy,2,2) - (e1* ylx2x1(ix,iy  ) + e2*xly2y1(ix  ,iy))
            cs%cf(ix,iy,2,2) = cs%cf(ix,iy,2,2) + (e1* ylx2x1(ix,iy+1) + e2*xly2y1(ix+1,iy))
         enddo
      enddo

      deallocate(x, y, r, xly2y1, ylx2x1)

      ! Impose symmetry relations

      do iy = iy0, my-1
         do ix = -mx, mx-1
            cs%cf(ix,iy,3,1) = -cs%cf(ix,iy,1,3)
            cs%cf(ix,iy,3,2) = -cs%cf(ix,iy,2,3)
            cs%cf(ix,iy,1,2) =  cs%cf(ix,iy,2,1)
         enddo
      enddo

   end subroutine elascf_pcwcns

!------------------------------------------------------------------------------------------------------------

   subroutine elascf_bilin (akv, nuv, mx, my, dx, dy, xshft, yshft, cs)
!--purpose: Calculates the influence coefficients according to Boussinesq and Cerruti with bilinear loading.
!           See memo m09002 and additions made in m17002 (CvdW)
!
!           This routine calculates the deformation at all grid points (xx,yy) of the potential contact
!           area due to a unit traction in the element with center (0,0). The grid may be shifted by an
!           offset (xshft,yshft).
!
!           The coefficients cf(ix-jx,iy-jy,:,:) say what the deformation is in element (ix,iy), when a
!           unit traction is located in element (jx,jy).
      implicit none
!--subroutine arguments:
      integer        :: mx, my
      real(kind=8)   :: akv, nuv, dx, dy, xshft, yshft
      type(t_inflcf) :: cs
!--local variables:
      integer      :: ix, iy, iy0, ik, jk
      real(kind=8) :: e1, e2, e3, tolx, toly, tolr
      type(t_xyr)  :: xyr

      allocate(xyr%x(-mx-1:mx), xyr%y(-my-1:my), xyr%r(-mx-1:mx,-my-1:my))
      allocate(xyr%xlog(-mx-1:mx,-my-1:my), xyr%ylog(-mx-1:mx,-my-1:my), xyr%rlog(-mx-1:mx,-my-1:my))
      allocate(xyr%xlogy(-mx-1:mx,-my-1:my), xyr%ylogx(-mx-1:mx,-my-1:my))
      allocate(xyr%xtanyx(-mx-1:mx,-my-1:my), xyr%ytanxy(-mx-1:mx,-my-1:my))

      e1 = (1d0-nuv) / (pi*dx*dy)
      e2 =  1d0      / (pi*dx*dy)
      e3 =      nuv  / (pi*dx*dy)

      tolx = dx*1d-10
      toly = dy*1d-10
      tolr = dx*dy*1d-10

      ! Initialize influence coefficients

      do jk = 1, 3
         do ik = 1, 3
            do iy = -my, my-1
               do ix = -mx, mx-1
                  cs%cf(ix,iy,ik,jk) = 0d0
               enddo
            enddo
         enddo
      enddo

      ! In 2D problems, the FFT preconditioner requires cs=0 at iy=-my

      iy0 = -my
      if (my.eq.1) iy0 = -my + 1 ! = 0

      ! compute coordinates of observation point (ix,iy)+shft relative to the four corners of element (0,0)

      ! Note: x and/or y may become zero or close to zero in the determination of cv, when the traversed
      !       distance per timestep contains a half grid width, i.e. xshft=0.5*dx, 1.5*dx, etc (chi=0).
      !       This yields r=|y| or r=|x| and r+y may be 0.

      do iy = iy0-1, my
         do ix = -mx-1, mx
            xyr%x(ix)    = real(ix)*dx + xshft
            xyr%y(iy)    = real(iy)*dy + yshft
            xyr%r(ix,iy) = dsqrt(xyr%x(ix)**2 + xyr%y(iy)**2)

            xyr%xlog(ix,iy) = 0.d0
            xyr%ylog(ix,iy) = 0.d0
            xyr%rlog(ix,iy) = 0.d0
            if (xyr%r(ix,iy).gt.tolr) then
               xyr%xlog(ix,iy) = xyr%x(ix)   *dlog(xyr%r(ix,iy))
               xyr%ylog(ix,iy) = xyr%y(iy)   *dlog(xyr%r(ix,iy))
               xyr%rlog(ix,iy) = xyr%r(ix,iy)*dlog(xyr%r(ix,iy))
            endif
      
            xyr%xlogy(ix,iy)  = 0.d0
            xyr%xtanyx(ix,iy) = 0.d0
            if (abs(xyr%x(ix)).gt.tolx) then
               xyr%xlogy(ix,iy)  = xyr%x(ix)*dlog(xyr%y(iy)+xyr%r(ix,iy))
               xyr%xtanyx(ix,iy) = xyr%x(ix)*datan(xyr%y(iy)/xyr%x(ix))
            endif

            xyr%ylogx(ix,iy)  = 0.d0
            xyr%ytanxy(ix,iy) = 0.d0
            if (abs(xyr%y(iy)).gt.toly) then
               xyr%ylogx(ix,iy)  = xyr%y(iy)*dlog(xyr%x(ix)+xyr%r(ix,iy))
               xyr%ytanxy(ix,iy) = xyr%y(iy)*datan(xyr%x(ix)/xyr%y(iy))
            endif
         enddo
      enddo

      ! The extended calculation (short calculation was used when mx > 1/eps)

      ! for all possible combinations of (ix-jx, iy-jy) :

      do iy = iy0, my-1
         do ix = -mx, mx-1

            ! cf(1,1) == A^ux == (1-v)*(1/r) + v*x^2/r^3

            cs%cf(ix,iy,1,1) = cs%cf(ix,iy,1,1) + e1*Hklm(xyr,ix,iy,001) + e3*Hklm(xyr,ix,iy,203)

            ! cf(2,2) == A^vy == (1-v)*(1/r) + v*y^2/r^3

            cs%cf(ix,iy,2,2) = cs%cf(ix,iy,2,2) + e1*Hklm(xyr,ix,iy,001) + e3*Hklm(xyr,ix,iy,023)

            ! cf(3,3) == A^wz == (1-v)*(1/r)

            cs%cf(ix,iy,3,3) = cs%cf(ix,iy,3,3) + e1*Hklm(xyr,ix,iy,001)

         enddo
      enddo

      do iy = iy0, my-1
         do ix = -mx, mx-1

            ! cf(2,1) == A^vx == v*x*y/r^3

            cs%cf(ix,iy,2,1) = cs%cf(ix,iy,2,1) + e3*Hklm(xyr,ix,iy,113)

         enddo
      enddo

      if (abs(akv).ge.1d-6) then
         do iy = iy0, my-1
            do ix = -mx, mx-1

               ! cf(3,1) == A^wx == K*x/r^2

               cs%cf(ix,iy,3,1) = cs%cf(ix,iy,3,1) + akv*e2*Hklm(xyr,ix,iy,102)

               ! cf(3,2) == A^wy == K*y/r^2

               cs%cf(ix,iy,3,2) = cs%cf(ix,iy,3,2) + akv*e2*Hklm(xyr,ix,iy,012)

            enddo
         enddo
      endif

      ! Impose symmetry relations

      do iy = iy0, my-1
         do ix = -mx, mx-1
            cs%cf(ix,iy,1,3) = -cs%cf(ix,iy,3,1)
            cs%cf(ix,iy,2,3) = -cs%cf(ix,iy,3,2)
            cs%cf(ix,iy,1,2) =  cs%cf(ix,iy,2,1)
         enddo
      enddo

   end subroutine elascf_bilin

!------------------------------------------------------------------------------------------------------------

   function Hklm(xyr, ix, iy, klm)
!--purpose: combines 4 Gklm-functions, summing over the 4 quadrants that make up the area of integration.
      implicit none
!--function result:
      real(kind=8) :: Hklm
!--subroutine arguments:
      type(t_xyr)  :: xyr
      integer      :: ix, iy, klm
!--local variables:

      Hklm = 0.d0
      Hklm = Hklm + Gklm(-xyr%x(ix+1),-xyr%y(iy+1), xyr, ix  ,iy  , klm)                                 &
                  - Gklm(-xyr%x(ix+1),-xyr%y(iy-1), xyr, ix  ,iy-1, klm)                                 &
                  - Gklm(-xyr%x(ix-1),-xyr%y(iy+1), xyr, ix-1,iy  , klm)                                 &
                  + Gklm(-xyr%x(ix-1),-xyr%y(iy-1), xyr, ix-1,iy-1, klm)

   end function Hklm

!------------------------------------------------------------------------------------------------------------

   function Gklm(a, b, xyr, ix, iy, klm)
!--purpose: combine 4 iFklm-functions, summing over the constant, x-linear, y-linear and bilinear parts
!           of A*phi.
      implicit none
!--function result:
      real(kind=8) :: Gklm
!--subroutine arguments:
      type(t_xyr)  :: xyr
      real(kind=8) :: a, b
      integer      :: ix, iy, klm
!--local variables:

      Gklm = 0.d0
      Gklm = Gklm + a*b*iFklm(xyr,ix,iy,klm    ) + b*iFklm(xyr,ix,iy,klm+100)
      Gklm = Gklm + a*  iFklm(xyr,ix,iy,klm+010) +   iFklm(xyr,ix,iy,klm+110)

   end function Gklm

!------------------------------------------------------------------------------------------------------------

   function iFklm(xyr, ix, iy, klm)
!--purpose: determine the integral of f_klm over a square centered around (xi,yi) with width/height Dx/Dy
      implicit none
!--function result:
      real(kind=8) :: iFklm
!--subroutine arguments:
      type(t_xyr)  :: xyr
      integer      :: ix, iy, klm
!--local variables:

      iFklm = 0.d0
      iFklm = iFklm + Fklm(xyr,ix+1,iy+1,klm) - Fklm(xyr,ix+1,iy  ,klm)
      iFklm = iFklm - Fklm(xyr,ix  ,iy+1,klm) + Fklm(xyr,ix  ,iy  ,klm)

   end function iFklm

!------------------------------------------------------------------------------------------------------------

   function Fklm(xyr,ix,iy,klm)
!--purpose: determine the value of integrating function F_klm of f_klm at the corner point (x,y)
      implicit none
!--function result:
      real(kind=8) :: Fklm
!--subroutine arguments:
      type(t_xyr),  target  :: xyr
      integer               :: ix, iy, klm
!--local variables:
      real(kind=8), pointer :: x, y, r, xlog, ylog, rlog, xlogy, ylogx, xtanyx, ytanxy

      x      => xyr%x(ix)
      y      => xyr%y(iy)
      r      => xyr%r(ix,iy)
      xlog   => xyr%xlog(ix,iy)
      ylog   => xyr%ylog(ix,iy)
      rlog   => xyr%rlog(ix,iy)
      xlogy  => xyr%xlogy(ix,iy)
      ylogx  => xyr%ylogx(ix,iy)
      xtanyx => xyr%xtanyx(ix,iy)
      ytanxy => xyr%ytanxy(ix,iy)


      if     (klm.eq.001) then
         Fklm = xlogy + ylogx
      elseif (klm.eq.101) then
         Fklm = 0.5d0*(x*xlogy + y*r)
      elseif (klm.eq.011) then
         Fklm = 0.5d0*(y*ylogx + x*r)
      elseif (klm.eq.111) then
         Fklm = r*r*r / 3.d0

      elseif (klm.eq.012) then
         Fklm = xlog + ytanxy
      elseif (klm.eq.102) then
         Fklm = ylog + xtanyx
      elseif (klm.eq.022) then
         Fklm = 0.5d0*(y*ytanxy - x*xtanyx)
      elseif (klm.eq.202) then
         Fklm = 0.5d0*(x*xtanyx - y*ytanxy)
      elseif (klm.eq.112) then
         Fklm = 0.5d0*r*rlog
      elseif (klm.eq.122) then
         Fklm = (y*y*ylog - x*x*xtanyx + y*x*x) / 3d0
      elseif (klm.eq.212) then
         Fklm = (x*x*xlog - y*y*ytanxy + x*y*y) / 3d0

      elseif (klm.eq.023) then
         Fklm = xlogy
      elseif (klm.eq.203) then
         Fklm = ylogx
      elseif (klm.eq.033) then
         Fklm = x*r
      elseif (klm.eq.303) then
         Fklm = y*r
      elseif (klm.eq.113) then
         Fklm = -r
      elseif (klm.eq.123) then
         Fklm = 0.5d0*(x*xlogy - y*r)
      elseif (klm.eq.213) then
         Fklm = 0.5d0*(y*ylogx - x*r)
      elseif (klm.eq.133) then
         Fklm = r*(2.d0*x*x - y*y) / 3d0
      elseif (klm.eq.313) then
         Fklm = r*(2.d0*y*y - x*x) / 3d0
      elseif (klm.eq.223) then
         Fklm = (x*x*xlogy + y*y*ylogx - x*y*r) / 3d0
      else
         Fklm = 1d10
      endif

   end function Fklm

!------------------------------------------------------------------------------------------------------------

   subroutine visccf (akv, nuv, ga, poiss, fdg, fg, vt, mx, my, dx, dy, xshft, yshft, cs, idebug)
!--purpose: Calculate the viscoelastic part of the influence coefficients by numerical integration.
!
!           This routine calculates the deformation at all grid points (xx,yy) of the potential contact
!           area due to a unit pressure in the element with center (0, 0). The grid may be shifted by
!           offset (xshft,yshft).
      implicit none
!--subroutine arguments:
      integer        :: mx, my, idebug
      real(kind=8)   :: akv, nuv, ga, poiss(2), fdg(2), fg(2), vt(2), dx, dy, xshft, yshft
      type(t_inflcf) :: cs
!--local variables:
      integer      :: ix, iy, iy0, ik, jk, ia, j, nvdif
      real(kind=8) :: cfv(3,3), e1, e2, e3, e1i, e3i, e4i, s1i, s31, sgn, ssnn, ssxn, ssxx, ssxy, ssyn, &
                      vti, vtm, x(2), xx, y(2), yy

      ! nvdif=1: viscoelastic contact of the same material
      ! nvdif=2: viscoelastic contact of the different material

      if (abs(akv).lt.epsco .and. abs(vt(1)-vt(2)).lt.epsco) then
         nvdif = 1
      else
         nvdif = 2
      endif

      if (idebug.ge.2) then
         write(bufout,'(2(a,f6.3),a,i2)') ' akv =',akv,', vt(1)-vt(2) =',vt(1)-vt(2),': nvdif =',nvdif
         call write_log(1,bufout)
         write(bufout,'(2(a,g12.4))') ' shift =',xshft,',', yshft
         call write_log(1,bufout)
      endif

      if (nvdif.eq.2) then
         if (vt(1).le.epsco .and. fg(1).gt.epsco) vt(1) = 10d0*epsco
         if (vt(2).le.epsco .and. fg(2).gt.epsco) vt(2) = 10d0*epsco
         ! if ((vt(1).le.epsco .and. fg(1).gt.epsco) .or. (vt(2).le.epsco .and. fg(2).gt.epsco))          &
         !    call write_log(' ERROR: different visco-elastic materials cant have TC=0.')
         ! if (vt(1).le.epsco .and. fg(1).gt.epsco)                                                       &
         !    call write_log('        the effect of FG is ignored for body 1')
         ! if (vt(2).le.epsco .and. fg(2).gt.epsco)                                                       &
         !    call write_log('        the effect of FG is ignored for body 2')
      endif

      vtm = 0.5d0*(vt(1)+vt(2))

      e1 = (1d0-nuv) / pi
      e2 =  1d0      / pi
      e3 =      nuv  / pi

      ! In 2D problems, the FFT preconditioner requires cs=0 at iy=-my

      iy0 = -my
      if (my.eq.1) iy0 = -my + 1

      ! The coefficients cf(ix-jx,iy-jy,:,:) say what the deformation is in element (ix,iy), when a
      ! unit pressure is located in element (jx,jy).

      ! for all possible combinations of (ix-jx, iy-jy) :

      do iy = iy0, my-1
         do ix = -mx, mx-1

            ! compute coordinates of observation point relative to (0,0)

            yy = iy * dy + yshft
            xx = ix * dx + xshft
            if (idebug.ge.10) then
               write(bufout,'(2(a,i4),2(a,f8.3),a)') ' starting (ix,iy)=(',ix,',',iy,'): (xx,yy)=(',     &
                        xx,',',yy,')'
               call write_log(1, bufout)
            endif

            ! compute distance vectors for the four corners of the central element

            x(1) = xx - 0.5*dx
            y(1) = yy - 0.5*dy
            x(2) = xx + 0.5*dx
            y(2) = yy + 0.5*dy

            do ik = 1, 3
               do jk = 1, 3
                  cfv(ik,jk)=0.d0
               enddo
            enddo

            ! The calculation of influence coefficients of the same material

            if (nvdif.eq.1) then

               vti = vtm
               if (vti.gt.epsco) then   ! tc>0: true visco-elastic
                  do j = 1, 2
                     sgn = 2*j - 3
                     call dnumin (fnn, 33, x(j), y, vti, ssnn)
                     call dnumin (fxx, 11, x(j), y, vti, ssxx)
                     call dnumin (fxy, 12, x(j), y, vti, ssxy)
                     cfv(3,3) = cfv(3,3) + fg(1)* e1*(ssnn+ssxx)*sgn
                     cfv(1,1) = cfv(1,1) + fg(1)*(e1*ssnn+e2*ssxx)*sgn
                     cfv(2,2) = cfv(2,2) + fg(1)*(e1*ssxx+e2*ssnn)*sgn
                     cfv(2,1) = cfv(2,1) + fg(1)*(-e3*ssxy)*sgn
                  enddo
               else                     ! rubbery state: effectively elastic
                  do ik = 1, 3
                     do jk = 1, 3
                        cs%cf (ix,iy, ik,jk) = cs%cf(ix,iy, ik,jk) * (1.d0+fg(1))
                     enddo
                  enddo
               endif

            ! The calculation of the influence coefficients of different materials

            elseif (nvdif.eq.2) then

               do ia = 1, 2
                  s1i = 3.d0-2.d0*ia
                  vti = vt(ia)
                  e1i = (1.d0-poiss(ia))/pi
                  e3i = poiss(ia)/pi
                  e4i = 0.5d0*(1.d0-2.d0*poiss(ia))/pi
                  if (vti.gt.epsco) then        ! tc>0: true visco-elastic
                     do j = 1, 2
                        sgn = 2*j-3
                        call dnumin (fnn, 33, x(j), y, vti, ssnn)
                        call dnumin (fxx, 11, x(j), y, vti, ssxx)
                        call dnumin (fxy, 12, x(j), y, vti, ssxy)
                        call dnumin (fxn, 13, x(j), y, vti, ssxn)
                        if ((dabs(y(1)+y(2))).lt.epsco) then
                           s31 = 0.d0
                        else
                           call dnumin (fyn, 23, x(j), y, vti, ssyn)
                        endif
                        cfv(3,3) = cfv(3,3) + 0.5d0*ga*fdg(ia) * sgn *  e1i*(ssnn+ssxx)
                        cfv(1,1) = cfv(1,1) + 0.5d0*ga*fdg(ia) * sgn * (e1i*ssnn+e2*ssxx)
                        cfv(2,2) = cfv(2,2) + 0.5d0*ga*fdg(ia) * sgn * (e1i*ssxx+e2*ssnn)
                        cfv(2,1) = cfv(2,1) + 0.5d0*ga*fdg(ia) * sgn * (-e3i*ssxy)
                        cfv(1,3) = cfv(1,3) + 0.5d0*ga*fdg(ia) * sgn *  e4i*ssxn*s1i
                        cfv(2,3) = cfv(2,3) + 0.5d0*ga*fdg(ia) * sgn *  e4i*ssyn*s1i
                     enddo ! j=1,2
                  else                     ! rubbery state: effectively elastic
                     ! TODO: not implemented!
                  endif
               enddo ! ia=1,2
            endif

            ! Impose symmetry relations

            cfv(3,1) = -cfv(1,3)
            cfv(3,2) = -cfv(2,3)
            cfv(1,2) =  cfv(2,1)

            ! Add to influence matrix cs at element (ix-jx,iy-jy)

            do ik = 1, 3
               do jk = 1, 3
                  cs%cf (ix,iy, ik,jk) = cs%cf(ix,iy, ik,jk) + cfv(ik,jk)
               enddo
            enddo
         enddo
      enddo

   end subroutine visccf

!------------------------------------------------------------------------------------------------------------

! The calculation of improper integrals.
! Implemented by G.Q. WANG, finished on 24 JANUARY 1992 at TU DELFT, room 825, HP-COMPUTER

   subroutine DNUMIN(FUNC, ifunc, xi, y, vti, ssk)
      implicit none
!--subroutine arguments:
      integer      :: ifunc
      real(kind=8) :: xi, y(2), vti, ssk
!--local variables:
      real(kind=8), parameter :: epsx = 1.0d-15
      integer      :: inf, ier, neval
      real(kind=8) :: epsabs, epsrel, bound, result, abserr
!--functions used:
      real(kind=8), external :: FUNC

      ssk = 0.0
      if (vti.gt.epsx) then
         inf = 1
         epsabs = 0.d0
         epsrel = 1.0d-03
         bound  = 0.d0
         call DQAGI(FUNC, xi, y, vti, bound, inf, epsabs, epsrel, result, abserr, neval, ier)
         if (ier.ne.0 .and. .not.(ier.eq.5 .and. abs(result).lt.1d-15)) then
            write(bufout,990) ifunc, xi, y(1), y(2), vti, bound, result
            call write_log(1, bufout)
            write(*,*) 'ERROR CODE =', ier

            if (ier.eq.1) then
                write (*,1001)
            elseif (ier.eq.2) then
                write (*,1002)
            elseif (ier.eq.3) then
                write (*,1003)
            elseif (ier.eq.4) then
                write (*,1004)
            elseif (ier.eq.5) then
                write (*,1005)
            endif
            write(*,*) 'INTEGRAL APPROXIMATION     =',result
            write(*,*) 'ESTIMATE OF ABSOLUTE ERROR =',abserr
            write(*,*) 'NUMBER OF FUNCTION EVALUATIONS =',neval
            write(*,*) 'PROBLEMS IN NUMERICAL INTEGRATION DQAGI DNUMIN !'
            call abort_run()
         endif
         ssk = result
      endif

 990  format( ' Problem in numerical integration for integral',i3,/                                      &
        ' xi=',g12.4,', y=[',g12.4,',', g12.4,'].',/                                                     &
        ' Vti=',f8.4,', bound=',g12.4,', ssk=',g12.4)
1001  format( ' MAXIMUM NUMBER OF ALLOWED SUBDIVISIONS HAS BEEN ACHIEVED. ONE CAN ALLOW MORE',/          &
       ' SUBDIVISIONS BY INCREASING THE DATA VALUE OF LIMIT IN DQAGI !')
1002  format( ' THE OCCURRENCE OF ROUNDOFF ERROR IS DETECTED, WHICH PREVENTS THE REQUESTED',/            &
       ' TOLERANCE FROM BEING ACHIEVED. THE ERROR MAY BE UNDER-ESTIMATED !')
1003  format( ' EXTREMELY BAD INTEGRAND BEHAVIOUR OCCURS AT SOME POINTS OF THE INTEGRATION',/            &
       ' INTERVAL.')
1004  format( ' THE ALGORITHM DOES NOT CONVERGE. ROUNDOFF ERROR IS DETECTED IN THE',/                    &
       ' EXTRAPOLATION TABLE. IT IS ASSUMED THAT THE REQUESTED TOLERANCE CANNOT BE ACHIEVED !')
1005  format(' THE INTEGRAL IS PROBABLY DIVERGENT, OR SLOWLY CONVERGENT !')

   end subroutine dnumin

!------------------------------------------------------------------------------------------------------------
! Definition of functions which will be integrated numerically
!------------------------------------------------------------------------------------------------------------

   function fnn(s,x,y,vti)
      implicit none
!--function result:
      real(kind=8) :: fnn
!--subroutine arguments:
      real(kind=8) :: s, x, y(2), vti
!--local variables:
      real(kind=8) :: u, ry1, ry2

      u = x+vti*s
      ry1=dsqrt(u**2+y(1)**2)
      ry2=dsqrt(u**2+y(2)**2)
      if (dabs(u).lt.epsco) then
         fnn = 0.0
      else
         fnn = dexp(-s)*u*dlog((y(2)+ry2)/(y(1)+ry1))
      endif

   end function fnn

!------------------------------------------------------------------------------------------------------------

   function fxx(s,x,y,vti)
      implicit none
!--function result:
      real(kind=8) :: fxx
!--subroutine arguments:
      real(kind=8) :: s, x, y(2), vti
!--local variables:
      integer      :: i
      real(kind=8) :: u, ry(2), fab(2)

      u = x+vti*s
      do 10 i = 1,2
         ry(i)=dsqrt(u**2+y(i)**2)
         if (dabs(y(i)).lt.epsco) then
             fab(i) = 0.0
         else
             fab(i) = dexp(-s)*y(i)*dlog(u+ry(i))
         endif
  10  continue
      fxx = fab(2)-fab(1)

   end function fxx

!------------------------------------------------------------------------------------------------------------

   function fxy(s, x, y, vti)
      implicit none
!--function result:
      real(kind=8) :: fxy
!--subroutine arguments:
      real(kind=8) :: s, x, y(2), vti
!--local variables:
      real(kind=8) :: u, ry1, ry2

      u   = x + vti * s
      ry1 = dsqrt(u**2 + y(1)**2)
      ry2 = dsqrt(u**2 + y(2)**2)
      fxy = dexp(-s) * (ry2 - ry1)

   end function fxy

!------------------------------------------------------------------------------------------------------------

   function fxn(s,x,y,vti)
      implicit none
!--function result:
      real(kind=8) :: fxn
!--subroutine arguments:
      real(kind=8) :: s, x, y(2), vti
!--local variables:
      integer      :: i
      real(kind=8) :: u, ry(2), f2ab(2)

      u = x+vti*s
      do 10 i = 1,2
         if (dabs(y(i)).lt.epsco .and. dabs(u).lt.epsco) then
           f2ab(i) = 0.0
         else
           ry(i) = dsqrt(u**2+y(i)**2)
           if (dabs(u).lt.epsco) then
             f2ab(i) = y(i)*dlog(ry(i))
           else
             f2ab(i) = y(i)*dlog(ry(i))+u*datan(y(i)/u)
           endif
         endif
  10  continue
      fxn = dexp(-s)*(f2ab(2)-f2ab(1))

   end function fxn

!------------------------------------------------------------------------------------------------------------

   function fyn(s,x,y,vti)
      implicit none
!--function result:
      real(kind=8) :: fyn
!--subroutine arguments:
      real(kind=8) :: s, x, y(2), vti
!--local variables:
      integer      :: i
      real(kind=8) :: u, ry(2), f3ab(2)

      u = x+vti*s
      do 10 i = 1,2
         if (dabs(y(i)).lt.epsco .and. dabs(u).lt.epsco) then
           f3ab(i) = 0.0
         else
           ry(i) = dsqrt(u**2+y(i)**2)
           if (dabs(y(i)).lt.epsco) then
             f3ab(i) = u*dlog(ry(i))
           else
             f3ab(i) = u*dlog(ry(i))+y(i)*datan(u/y(i))
           endif
         endif
   10 continue
      fyn = dexp(-s)*(f3ab(2)-f3ab(1))

   end function fyn

!------------------------------------------------------------------------------------------------------------

   subroutine DQAGI(FUNC, x, y, vti, bound, inf, epsabs, epsrel, result, abserr, neval, ier)
!--purpose: Numerical integration over infinite intervals (standard fortran subroutine) from "QUADPACK" -
!           a subroutine package for automatic integration by R.Piessens; E. De Doncker-Kapenga;
!           C.W.Ueberhuber; D.K.Kahaner, Springer-Verlag Berlin Heidelberg New York Tokyo 1983,
!           ISBN 3-540-12553-1
!           Written and tested here by G.WANG on january 11, 1992 in TU Delft (Prof. J.J.Kalker),
!           room 8.25 using a HP-workstation 700.
      implicit none
!--subroutine arguments:
      real(kind=8) x, y(2), vti, bound, epsabs, epsrel, result, abserr
      integer neval, ier, inf
      real(kind=8), external :: FUNC
!--local variables:
      real(kind=8) abseps, alist, area, area1, area12, area2, a1, a2, blist, boun, b1, b2, correc,      &
          defabs, defab1, defab2, dres, elist, epmach, erlarg, erlast, errbnd, errmax, error1,          &
          error2, erro12, errsum, ertest, oflow, resabs, reseps, res3la, rlist, rlist2, small, uflow
      integer id, ierro, iord, iroff1, iroff2, iroff3, jupbnd, k, ksgn, ktmin, last, limit, maxerr,     &
          nres, nrmax, numrl2
      logical extrap,noext

      dimension alist(500), blist(500), elist(500), iord(500), res3la(3), rlist(500), rlist2(52)
      data limit/500/

      call DQMACO(epmach,uflow,oflow)
      ier = 0
      neval = 0
      last = 0
      result = 0.0
      abserr = 0.0
      erlarg = 0d0
      alist(1) = 0.0
      blist(1) = 1.0
      rlist(1) = 0.0
      elist(1) = 0.0
      iord(1) = 0
      if (epsabs.lt.0.0 .and. epsrel.lt.0.0) ier = 6
      if (ier.eq.6) goto 999

      boun = bound
      if (inf.eq.2) boun = 0.0

      call DQK15I(FUNC,x,y,vti,boun,inf,0.0D0,1.0D0,result,abserr, defabs,resabs)

      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
      dres = dabs(result)
      errbnd = dmax1(epsabs,epsrel*dres)
      if (abserr.le.1.0d+02*epmach*defabs .and. abserr.gt.errbnd) ier=2
      if (limit.eq.1) ier = 1
      if (ier.ne.0 .or. (abserr.le.errbnd .and. abserr.ne.resabs) .or. abserr.eq.0.0) goto 130

      rlist2(1) = result
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      abserr = oflow
      nrmax = 1
      nres = 0
      ktmin = 0
      numrl2 = 2
      extrap = .false.
      noext = .false.
      ierro = 0
      iroff1 = 0
      iroff2 = 0
      iroff3 = 0
      ksgn = -1
      if (dres.ge.(1.0d+00-5.0d+01*epmach)*defabs) ksgn = 1

      do 90 last = 2,limit
         a1 = alist(maxerr)
         b1 = 5.0d-01*(alist(maxerr)+blist(maxerr))
         a2 = b1
         b2 = blist(maxerr)
         erlast = errmax

         call DQK15I(FUNC,x,y,vti,boun,inf,a1,b1,area1,error1,resabs,defab1)
         call DQK15I(FUNC,x,y,vti,boun,inf,a2,b2,area2,error2,resabs,defab2)

         area12 = area1 + area2
         erro12 = error1 + error2
         errsum = errsum+erro12-errmax
         area = area+area12-rlist(maxerr)
         if (defab1.eq.error1 .or. defab2.eq.error2) goto 15
         if (abs(rlist(maxerr)-area12).gt.1.0d-05*dabs(area12).or.erro12.lt.9.9d-01*errmax) goto 10
         if (extrap) iroff2 = iroff2+1
         if (.not.extrap) iroff1 = iroff1+1
  10     if (last.gt.10 .and. erro12.gt.errmax) iroff3 = iroff3+1
  15     rlist(maxerr) = area1
         rlist(last) = area2
         errbnd = dmax1(epsabs,epsrel*dabs(area))

         ! error in the following source programm (iroff1+iroff2) must be pp.172
         !       if (iroff1+iroff2.ge.10 .or. iroff3.ge.20) ier = 2

         if ((iroff1+iroff2).ge.10 .or. iroff3.ge.20) ier = 2
         if (iroff2.ge.5) ierro = 3
         if (last.eq.limit) ier = 1

         if (dmax1(dabs(a1),dabs(b2)).le.(1.0d+00+1.0d+03*epmach)*(dabs(a2)+1.0d+03*uflow)) ier = 4

         if (error2.gt.error1) goto 20
         alist(last) = a2
         blist(maxerr) = b1
         blist(last) = b2
         elist(maxerr) = error1
         elist(last) = error2
         goto 30
  20     alist(maxerr) = a2
         alist(last) = a1
         blist(last) = b1
         rlist(maxerr) = area2
         rlist(last) = area1
         elist(maxerr) = error2
         elist(last) = error1

  30     call DQSORT(limit,last,maxerr,errmax,elist,iord,nrmax)
         if (errsum.le.errbnd) goto 115
         if (ier.ne.0) goto 100
         if (last.eq.2) goto 80
         if (noext) goto 90
         erlarg = erlarg-erlast
         if (dabs(b1-a1).gt.small) erlarg = erlarg+erro12
         if (extrap) goto 40

         if (dabs(blist(maxerr)-alist(maxerr)).gt.small) goto 90
         extrap = .true.
         nrmax = 2
  40     if (ierro.eq.3 .or. erlarg.le.ertest) goto 60

         id = nrmax
         jupbnd = last
         if (last.gt.(2+limit/2)) jupbnd = limit+3-last
         do 50 k = id, jupbnd
            maxerr = iord(nrmax)
            errmax = elist(maxerr)
            if (abs(blist(maxerr)-alist(maxerr)).gt.small) goto 90
            nrmax = nrmax+1
  50     continue

  60     numrl2 = numrl2+1
         rlist2(numrl2) = area
         call DQEXT(numrl2,rlist2,reseps,abseps,res3la,nres)
         ktmin = ktmin+1
         if (ktmin.gt.5 .and. abserr.lt.1.0d-03*errsum) ier = 5
         if (abseps.ge.abserr) goto 70
         ktmin = 0
         abserr = abseps
         result = reseps
         correc = erlarg
         ertest = dmax1(epsabs,epsrel*dabs(reseps))
         if (abserr.le.ertest) goto 100

  70     if (numrl2.eq.1) noext = .true.
         if (ier.eq.5) goto 100
         maxerr = iord(1)
         errmax = elist(maxerr)
         nrmax = 1
         extrap = .false.
         small = small*5.0d-1
         erlarg = errsum
         goto 90
  80     small = 3.75d-01
         erlarg = errsum
         ertest = errbnd
         rlist2(2) = area
  90  continue

  100 if (abserr.eq.oflow) goto 115
      if ((ier+ierro).eq.0) goto 110
      if (ierro.eq.3) abserr = abserr+correc
      if (ier.eq.0) ier = 3
      if (result.ne.0.0d+00 .and. area.ne.0.0d+00) goto 105
      if (abserr.gt.errsum) goto 115
      if (area.eq.0.0d+00) goto 130
      goto 110
  105 if (abserr/dabs(result).gt.errsum/dabs(area)) goto 115

  110 if (ksgn.eq.(-1) .and. dmax1(dabs(result),dabs(area)).le. defabs*1.0d-02) goto 130
      if (1.0d-02.gt.(result/area) .or. (result/area).gt.1.0d+02 .or. errsum.gt.dabs(area)) ier = 6
      goto 130

  115 result = 0.0d+00
      do 120 k = 1,last
         result = result+rlist(k)
  120 continue
      abserr = errsum
  130 neval = 30*last-15
      if (inf.eq.2) neval = 2*neval
      if (ier.gt.2) ier=ier-1

  999 return

   end subroutine dqagi

!------------------------------------------------------------------------------------------------------------

   subroutine DQK15I(FUNC, x, y, vti, boun, inf, a, b, result, abserr, resabs, resasc)
      implicit none
!--subroutine arguments:
      real(kind=8) x, y(2), vti, boun, a, b, result, abserr, resabs, resasc
      integer      inf
      real(kind=8), external :: FUNC
!--local variables:
      real(kind=8) absc, absc1, absc2, centr, dinf, epmach, fc, fsum, fval1, fval2, fv1, fv2,           &
                   hlgth, oflow, resg, resk, reskh, tabsc1, tabsc2, uflow, wg, wgk, xgk
      integer j, ndinf

      dimension fv1(7), fv2(7), wgk(8), wg(8), xgk(8)

      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/                                     &
           9.914553711208126d-01,     9.491079123427585d-01,                                            &
           8.648644233597691d-01,     7.415311855993944d-01,                                            &
           5.860872354676911d-01,     4.058451513773972d-01,                                            &
           2.077849550078985d-01,     0.0/

      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/                                     &
           2.293532201052922d-02,     6.309209262997855d-02,                                            &
           1.047900103222502d-01,     1.406532597155259d-01,                                            &
           1.690047266392679d-01,     1.903505780647854d-01,                                            &
           2.044329400752989d-01,     2.094821410847278d-01/

      data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/                                             &
           0.0,                       1.294849661688697d-01,                                            &
           0.0,                       2.797053914892767d-01,                                            &
           0.0,                       3.818300505051189d-01,                                            &
           0.0,                       4.179591836734694d-01/

      call DQMACO(epmach,uflow,oflow)
      ndinf = min0(1,inf)
      dinf = dble(ndinf)

      centr = 0.5D0*(a+b)
      hlgth = 0.5D0*(b-a)
      tabsc1 = boun+dinf*(1.0d+00-centr)/centr
      fval1 = FUNC(tabsc1,x,y,vti)
      if (inf.eq.2) fval1 = fval1+FUNC(-tabsc1,x,y,vti)
      fc = (fval1/centr)/centr

      resg = wg(8)*fc
      resk = wgk(8)*fc
      resabs = dabs(resk)
      do 10 j = 1,7
        absc = hlgth*xgk(j)
        absc1 = centr-absc
        absc2 = centr+absc
        tabsc1 = boun+dinf*(1.0d+00-absc1)/absc1
        tabsc2 = boun+dinf*(1.0d+00-absc2)/absc2
        fval1 = FUNC(tabsc1,x,y,vti)
        fval2 = FUNC(tabsc2,x,y,vti)
        if (inf.eq.2) fval1 = fval1+FUNC(-tabsc1,x,y,vti)
        if (inf.eq.2) fval2 = fval2+FUNC(-tabsc2,x,y,vti)
        fval1 = (fval1/absc1)/absc1
        fval2 = (fval2/absc2)/absc2
        fv1(j) = fval1
        fv2(j) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(j)*fsum
        resabs = resabs+wgk(j)*(dabs(fval1)+dabs(fval2))
   10 continue
      reskh = resk*5.0d-01
      resasc = wgk(8)*dabs(fc-reskh)
      do 20 j = 1,7
         resasc = resasc+wgk(j)*(dabs(fv1(j)-reskh)+dabs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resasc = resasc*hlgth
      resabs = resabs*hlgth
      abserr = dabs((resk-resg)*hlgth)
      if (resasc.ne.0.0D0 .and. abserr.ne.0.0D0)                                                        &
         abserr = resasc* DMIN1(1.0D0,(2.0d+02*abserr/resasc)**1.5d+00)
      if (resabs.gt.uflow/(5.0d+01*epmach)) abserr = dmax1((epmach*5.0d+01)*resabs,abserr)

   end subroutine dqk15i

!------------------------------------------------------------------------------------------------------------

   subroutine DQSORT(limit,last,maxerr,ermax,elist,iord,nrmax)
      implicit none
!--subroutine arguments:
      real(kind=8) ermax, elist(last)
      integer limit, last, maxerr, iord(last), nrmax
!--local variables:
      real(kind=8) errmax, errmin
      integer i, ibeg, ido, isucc, j, jbnd, jupbn, k

      if (last.gt.2) goto 10
      iord(1) = 1
      iord(2) = 2
      goto 90
   10 errmax = elist(maxerr)
      if (nrmax.eq.1) goto 30
      ido = nrmax-1
      do 20 i = 1,ido
         isucc = iord(nrmax-1)
         if (errmax.le.elist(isucc)) goto 30
         iord(nrmax) = isucc
         nrmax = nrmax-1
   20 continue
   30 jupbn = last
      if (last.gt.(limit/2+2)) jupbn = limit+3-last
      errmin = elist(last)
      jbnd = jupbn-1
      ibeg = nrmax+1
      if (ibeg.gt.jbnd) goto 50
      do 40 i=ibeg,jbnd
         isucc = iord(i)
         if (errmax.ge.elist(isucc)) goto 60
         iord(i-1) = isucc
   40 continue
   50 iord(jbnd) = maxerr
      iord(jupbn) = last
      goto 90

   60 iord(i-1) = maxerr
      k = jbnd
      do 70 j=1,jbnd
         isucc = iord(k)
         if (errmin.lt.elist(isucc)) goto 80
         iord(k+1) = isucc
         k = k-1
   70 continue
      iord(i) = last
      goto 90
   80 iord(k+1) = last
   90 maxerr = iord(nrmax)
      ermax = elist(maxerr)

   end subroutine dqsort

!------------------------------------------------------------------------------------------------------------

   subroutine DQEXT(n,epstab,result,abserr,res3la,nres)
      implicit none
!--subroutine arguments:
      real(kind=8) epstab(52), result, abserr, res3la(3)
      integer n, nres
!--local variables:
      real(kind=8) delta1, delta2, delta3, epmach, epsinf, error, err1, err2, err3, e0, e1, e1abs,      &
                   e2, e3, oflow, res, ss, tol1, tol2, tol3, uflow
      integer i, ib, ib2, ie, indx, k1, k2, k3, limexp, newelm, num

      call DQMACO(epmach,uflow,oflow)
      nres = nres+1
      abserr = oflow
      result = epstab(n)
      if (n.lt.3) goto 100
      limexp = 50
      epstab(n+2) = epstab(n)
      newelm = (n-1)/2
      epstab(n) = oflow
      num = n
      k1 = n

      do 40 i = 1,newelm
         k2 = k1-1
         k3 = k1-2
         res = epstab(k1+2)
         e0 = epstab(k3)
         e1 = epstab(k2)
         e2 = res
         e1abs = dabs(e1)
         delta2 = e2-e1
         err2 = dabs(delta2)
         tol2 = dmax1(dabs(e2),e1abs)*epmach
         delta3 = e1-e0
         err3 = dabs(delta3)
         tol3 = dmax1(e1abs,dabs(e0))*epmach
         if (err2.gt.tol2 .or. err3.gt.tol3) goto 10

         result = res
         abserr = err2+err3
         goto 100
   10    e3 = epstab(k1)
         epstab(k1) = e1
         delta1 = e1-e3
         err1 = abs(delta1)
         tol1 = dmax1(e1abs,dabs(e3))*epmach
         if (err1.le.tol1 .or. err2.le.tol2 .or. err3.le.tol3) goto 20
         ss = 1.0D0/delta1+1.0D0/delta2-1.0D0/delta3
         epsinf = abs(ss*e1)
         if (epsinf.gt.1.0d-04) goto 30
   20    n = i+i-1
         goto 50
   30    res = e1+1.0D0/ss
         epstab(k1) = res
         k1 = k1-2
         error = err2+dabs(res-e2)+err3
         if (error.gt.abserr) goto 40
         abserr = error
         result = res
   40 continue

   50 if (n.eq.limexp) n = 2*(limexp/2)-1
      ib = 1
      if ((num/2)*2.eq.num) ib = 2
      ie = newelm+1
      do 60 i = 1,ie
         ib2 = ib+2
         epstab(ib) = epstab(ib2)
         ib = ib2
   60 continue
      if (num.eq.n) goto 80
      indx = num-n+1
      do 70 i = 1,n
         epstab(i) = epstab(indx)
         indx = indx+1
   70 continue
   80 if (nres.ge.4) goto 90
      res3la(nres) = result
      abserr = oflow
      goto 100

   90 abserr = dabs(result-res3la(3)) + dabs(result-res3la(2)) + dabs(result-res3la(1))
      res3la(1) = res3la(2)
      res3la(2) = res3la(3)
      res3la(3) = result
  100 abserr = dmax1(abserr,0.5D0*epmach*dabs(result))

   end subroutine dqext

!------------------------------------------------------------------------------------------------------------

   subroutine DQMACO(epmach,uflow,oflow)
! Attention: epmach, uflow and oflow are machine dependent constants, which may be quite different
!            between different computers, the following values are only valid for the HP-workstation/700.
!            For IBM 3033 they may be epmach=1.0d-7, uflow=1.0d-37 and oflow=1.0d+38 (in single precision)
!     epmach - the largest relative spacing
!     uflow  - the smallest positive magnitude
!     oflow  - the largest positive magnitude
      implicit none
!--subroutine arguments:
      real(kind=8) epmach,uflow,oflow

      epmach = 1.0d-10
      uflow = 1.0d-40
      oflow = 1.0d+30

   end subroutine dqmaco

!------------------------------------------------------------------------------------------------------------

end module m_visc
