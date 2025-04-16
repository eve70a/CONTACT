!------------------------------------------------------------------------------------------------------------
! test_mbench - This program shows the use of the CONTACT add-on interface for wheel/rail contact
!               situations ("module 1").
! It is derived from the Manchester benchmark example, as shown in Section 5.7 of the User guide.
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
program test_mbench
   use, intrinsic        :: iso_c_binding, only: C_CHAR, C_NULL_CHAR
   implicit none
   include 'caddon_flags.inc'
!--local variables:
   real(kind=8), parameter :: pi     = 4d0*atan(1d0)
   integer                 :: nunits, iunits(5)
   integer                 :: ldebug, ver, ierr, iout
   integer                 :: nwhe, iwhe, icp, imodul, iblk, isubs, gdigit, mdigit, ldigit
   integer, parameter      :: mxflgs = CNTC_FLAGS_DIM
   integer                 :: flags(mxflgs), values(mxflgs)
   real(kind=8)            :: gg, poiss, fstat, fkin, dx, ds, a_sep, d_sep, d_comb, dqrel, fz, rdum
   integer                 :: itype, iside, mirrory, ewheel, ztrack, ipotcn, istep, nstep, iparam(3)
   real(kind=8)            :: sclfac, smooth, rparam(2), rarr(1), rvalues(mxflgs)
   real(kind=8), dimension(:), allocatable :: y_ws, yaw_ws, roll_ws, vpitch
   integer                        :: len_wrkdir, len_outdir, len_expnam, len_fname
   character(len=256)             :: f_wrkdir, f_outdir, f_expnam, f_fname
   character(len=256,kind=C_CHAR) :: c_wrkdir, c_outdir, c_expnam, c_fname
!--external functions used:
#include "contact_addon.ifc"

   !---------------------------------------------------------------------------------------------------------
   ! Initialization of the CONTACT library
   !---------------------------------------------------------------------------------------------------------

   nunits = 5
   iunits = (/ 11, 12, 13, 14, 15 /)
   call cntc_setFileUnits(nunits, iunits)

   ldebug = 1

   ! Initialize using cntc_initializeFirst
   ! this may be omitted, relying on automatic initialization

   ver         = -1
   ierr        = -1
   iout        =  0     ! enable/disable (1/0) screen output
   f_wrkdir    = ' '    ! wrkdir = base for relative paths
   c_wrkdir    =     trim(f_wrkdir)  // C_NULL_CHAR
   len_wrkdir  = len(trim(f_wrkdir))
   f_outdir    = ' '    ! could be '/tmp' for instance
   c_outdir    =     trim(f_outdir)  // C_NULL_CHAR
   len_outdir  = len(trim(f_outdir))
   f_expnam    = ' '
   c_expnam    =     trim(f_expnam )  // C_NULL_CHAR
   len_expnam  = len(trim(f_expnam ))
   call cntc_initializeFirst(ver, ierr, iout, c_wrkdir, c_outdir, c_expnam, len_wrkdir, len_outdir,     &
                len_expnam)

   ! Set the amount of output of the C.library itself; enable or disable use of OpenMP

   flags(1) = CNTC_if_idebug ; values(1) = 1
   flags(2) = CNTC_if_openmp ; values(2) = 1    ! max #threads
   call cntc_setGlobalFlags(2, flags, values)

   !---------------------------------------------------------------------------------------------------------
   ! Initialize two result elements, for the left & right wheels, using module 1 for wheel/rail contact
   !---------------------------------------------------------------------------------------------------------

   nwhe = 2

   do iwhe = 1, nwhe
      imodul = 1
      call cntc_initialize(iwhe, imodul, ver, ierr, c_outdir, len_outdir)
   enddo

   ! Configure both result elements, setting permanent data

   do iwhe = 1, nwhe

      ! cp creation & numbering is automatic in module 1

      icp = -1

      ! set flags: configure control digits & output

      flags = 0 ; values = 0
      flags( 1) = CNTC_if_units  ; values( 1) = CNTC_un_cntc  ! CONTACT unit convention: mm
      flags( 2) = CNTC_ic_config ; values( 2) = iwhe-1        ! C1: 0=left wheel, 1=right
      flags( 3) = CNTC_ic_tang   ; values( 3) = 3             ! T=3: steady state rolling
      flags( 4) = CNTC_ic_pvtime ; values( 4) = 2             ! P=2: no previous time
      flags( 5) = CNTC_ic_discns ; values( 5) = 2             ! D=2: planar contact
      flags( 6) = CNTC_if_wrtinp ; values( 6) = 0             !   0: no inp-file needed
      flags( 7) = CNTC_ic_matfil ; values( 7) = 0             ! A=0: no mat-files needed
      flags( 8) = CNTC_ic_output ; values( 8) = 2             ! O=1: little output to out-file
      flags( 9) = CNTC_ic_flow   ; values( 9) = 2             ! W=2: some flow trace
      flags(10) = CNTC_ic_npomax ; values(10) = 20000         ! max #elements in pot.contact

      call cntc_setFlags(iwhe, icp, mxflgs, flags, values)

      ! G=0: set maxgs, maxin, maxnr, maxout, eps

      gdigit = 0
      values(1:4) = (/ 999, 100, 30, 1 /)
      rvalues(1)  = 1d-5
      call cntc_setSolverFlags(iwhe, icp, gdigit, 4, values, 1, rvalues)

      ! set material parameters using M = 0 (fully elastic model)

      mdigit = 0
      gg     = 82000d0  ! CONTACT units: [N/mm2]
      poiss  = 0.28d0   ! [-]
      rvalues(1:4) = (/ poiss, poiss, gg, gg /)
      call cntc_setMaterialParameters(iwhe, icp, mdigit, 4, rvalues)

      ! set friction parameters using  L = 0: Coulomb friction

      ldigit = 0
      fstat   = 0.3d0  ! [-]
      fkin    = fstat  ! [-]
      rvalues(1:2) = (/ fstat, fkin /)
      call cntc_setFrictionMethod(iwhe, icp, ldigit, 2, rvalues)

      ! friction variation across rail profile:

      ! ldigit = 10     ! V-digit 1, L-digit 0: Coulomb friction, variable across rail
      ! nvf = 2         ! 2 control points with linear interpolation, constant extrapolation
      ! params = (/ nvf,                           &
      !            -20d0*pi/180d0, 0.20d0, 0.20d0, &   ! constant [0.2,0.2] for surface inclination <= -20deg
      !            -10d0*pi/180d0, 0.30d0, 0.30d0 /)   ! linear 0.2 - 0.3 for inclination -20 -- -10 deg
      ! call cntc_setfrictionmethod(iwhe, icp, ldigit, 7, params)

      ! set grid discretization

      ipotcn = -1
      dx     = 0.2d0   ! [mm]
      ds     = 0.2d0   ! [mm]
      a_sep  = pi/2d0  ! [rad]
      d_sep  = 8.0d0   ! [mm]
      d_comb = 4.0d0   ! [mm]

      rvalues(1:5) = (/ dx, ds, a_sep, d_sep, d_comb /)
      call cntc_setPotContact(iwhe, icp, ipotcn, 5, rvalues)

      ! set rolling step size: ratio c = dq / dx

      rdum   = 0d0     ! chi: ignored
      dqrel  = 1d0     ! [-]
      call cntc_setRollingStepsize(iwhe, icp, rdum, dqrel)

      ! set track dimensions & deviations

      ztrack = 3
      rvalues(1:11) = (/ 14d0, 0d0, 1435d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 /)
      call cntc_setTrackDimensions(iwhe, ztrack, 11, rvalues)

      ! set rail profile for left and right sides

      f_fname   = 'profiles/MBench_UIC60_v3.prr'
      c_fname   =     trim(f_fname)  // C_NULL_CHAR
      len_fname = len(trim(f_fname))

      itype   = -1      ! using filename extension
      iside   =  2      ! left and right rails
      mirrory =  0      ! no mirroring
      sclfac  = 1d0     ! already in [mm], no scaling
      smooth  = 0d0     ! no smoothing
      iparam  = (/ itype, iside, mirrory /)
      rparam  = (/ sclfac, smooth /)

      call cntc_setProfileInputFname(iwhe, c_fname, len_fname, 3, iparam, 2, rparam)

      ! set wheelset dimensions

      ewheel = 3
      rvalues(1:3) = (/ 1360d0, -70d0, 460d0 /)
      call cntc_setWheelsetDimensions(iwhe, ewheel, 3, rvalues)

      ! set wheel profile. If two profiles are given, first is for left side, second for right side

      f_fname   = 'profiles/MBench_S1002_v3.prw'
      c_fname   =     trim(f_fname)  // C_NULL_CHAR
      len_fname = len(trim(f_fname))

      itype   = -1      ! using filename extension
      iside   =  2      ! left and right wheels
      mirrory =  0      ! no mirroring
      sclfac  = 1d0     ! already in [mm], no scaling
      smooth  = 0d0     ! no smoothing
      iparam  = (/ itype, iside, mirrory /)
      rparam  = (/ sclfac, smooth /)

      call cntc_setProfileInputFname(iwhe, c_fname, len_fname, 3, iparam, 2, rparam)

      ! positions used for subsurface stress calculation

      iblk  = 1
      isubs = 1                               ! X,Y: all element centers - using FFTs
      rvalues(1:3) = (/ 11d0, 1d-9, 0.5d0 /)  !   Z: NZ, ZL, DZ

      call subs_addBlock(iwhe, icp, iblk, isubs, 0, 0, 3, rarr, rarr, rvalues)

   end do

   !---------------------------------------------------------------------------------------------------------
   ! Compute results for 21 wheelset positions
   !---------------------------------------------------------------------------------------------------------

   nstep = 21
   allocate(y_ws(nstep), yaw_ws(nstep), roll_ws(nstep), vpitch(nstep))
   y_ws(1:nstep)    = (/  0.0d0, 0.5d0, 1.0d0, 1.5d0, 2.0d0, 2.5d0, 3.0d0, 3.5d0, 4.0d0, 4.5d0,         &
                          5.0d0, 5.5d0, 6.0d0, 6.5d0, 7.0d0, 7.5d0, 8.0d0, 8.5d0, 9.0d0, 9.5d0,         &
                         10.0d0 /)
   yaw_ws(1:nstep)  =     0.0024d0 * y_ws(1:nstep)
   roll_ws(1:nstep) = (/  0.00000000d0, -0.00002304d0, -0.00005049d0, -0.00008103d0, -0.00011280d0,     &
                         -0.00014570d0, -0.00018030d0, -0.00021680d0, -0.00025570d0, -0.00029770d0,     &
                         -0.00035540d0, -0.00047770d0, -0.00062720d0, -0.00437600d0, -0.00639300d0,     &
                         -0.00764200d0, -0.00860600d0, -0.00940800d0, -0.01010113d0, -0.01071386d0,     &
                         -0.01126431d0 /)
   vpitch(1:nstep)  = (/ -4.34811810d0, -4.34741340d0, -4.34657520d0, -4.34624400d0, -4.34591970d0,     &
                         -4.34556270d0, -4.34515030d0, -4.34466880d0, -4.34409300d0, -4.34337150d0,     &
                         -4.33536370d0, -4.33188640d0, -4.32937180d0, -4.27488340d0, -4.26356290d0,     &
                         -4.25757470d0, -4.25348570d0, -4.25032450d0, -4.24775610d0, -4.24556650d0,     &
                         -4.24363750d0 /)

   do istep = 1, nstep
      do iwhe = 1, nwhe

         ! set wheelset position & velocity

         rvalues(1:6) = (/ 0d0, y_ws(istep), 0d0, roll_ws(istep), yaw_ws(istep), 0d0 /)
         ewheel = 2
         call cntc_setWheelsetPosition(iwhe, ewheel, 6, rvalues)

         rvalues(1:6) = (/ 2000d0, 0d0, 0d0, 0d0, 0d0, vpitch(istep) /)
         ewheel = 2
         call cntc_setWheelsetVelocity(iwhe, ewheel, 6, rvalues)

         ! set vertical force, settign the N-digit to 1

         fz  = 10000d0    ! [N]
         call cntc_setVerticalForce(iwhe, fz)

         if (ldebug.ge.1)  &
            write(*,'(/,a,i3)') '  ...start processing of Result Element', iwhe

         ! perform the actual calculation for this wheel

         call cntc_calculate(iwhe, icp, ierr)

         ! check/report error conditions

         if (ierr.eq.CNTC_err_allow) then
            write(*,'(a,i5,a)') '      ERROR: no valid license found'
            stop
         elseif (ierr.eq.CNTC_err_profil) then
            write(*,'(a,i5,a)') '      ERROR: the rail and/or wheel profile files could not be ' //     &
                                'found or processed'
            stop
         elseif (ierr.lt.0) then
            write(*,'(a,i5,a)') '      ERROR: an error occured in the CONTACT library, ierr=',ierr,'.'
         elseif (ierr.gt.0) then
            write(*,'(a,i5,a)') '      ERROR: contact extends outside bounding box,',ierr,              &
                                ' elements at boundary'
         end if

         ! retrieve & print results

         call test_results(iwhe, ldebug)

         ! perform subsurface stress calculation when needed

         if (iwhe.eq.2 .and. istep.eq.11) call max_vonmises(iwhe, ldebug)

      end do ! iwhe
   end do ! istep

   deallocate(y_ws, yaw_ws, roll_ws, vpitch)

   ! finalize, close CONTACT library

   do iwhe = 1, nwhe
      call cntc_finalize(iwhe)
   end do

end program test_mbench

!------------------------------------------------------------------------------------------------------------

subroutine max_vonmises(iwhe, ldebug)
!--function: compute subsurface stresses, get max. von Mises and print to the screen
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: iwhe           ! result element ID
   integer,      intent(in) :: ldebug         ! level of output requested
!--local variables:
   integer               :: ierr, icp, jcp, ncp, iblk, nx, ny, nz, npnt, ncol, icol(21), ii, ii_max
   real(kind=8), dimension(:,:), allocatable :: table
!--external functions used:
#include "contact_addon.ifc"

   ! start the actual calculation of subsurface stresses, for the blocks defined at initialization

   icp = -1
   call subs_calculate(iwhe, icp, ierr)

   if (ierr.lt.0) then
      write(*,'(a,i5,a)') '      ERROR: in the calculation of subsurface stresses, ierr=',ierr,'.'
   endif

   ! loop over all contact patches, to get / display the results

   call cntc_getNumContactPatches(iwhe, ncp)
   do jcp = 1, ncp

      if (ldebug.ge.2) write(*,'(2x,a,i3,a)') 'Patch',jcp,':'

      ! get number of points used in subsurf. calculation

      iblk = 1
      call subs_getBlocksize(iwhe, jcp, iblk, nx, ny, nz)
      if (ldebug.ge.1) write(*, '(5x,3(a,i4),a)') 'subsurface stresses computed at', nx,' x',ny,' x',   &
                nz,' points'

      ! get results 1=x, 2=y, 3=z, 8=sigvm

      npnt = nx * ny * nz
      ncol = 4
      icol(1:ncol) = (/ 1, 2, 3, 8 /)
      allocate(table(npnt,ncol))

      call subs_getResults(iwhe, jcp, iblk, npnt, ncol, icol, table)

      ! locate maximum von Mises stress in table

      ii_max = 1
      do ii = 1, npnt
         if (abs(table(ii,4)).gt.abs(table(ii_max,4))) ii_max = ii
      enddo

      ! print maximum von Mises stress & corresponding location

      if (ldebug.ge.1) write(*, 2061) 'Max: Sigvm=', table(ii_max,4), table(ii_max,1), table(ii_max,2), &
                table(ii_max,3)
 2061 format(5x,a,f10.3,' [N/mm2] at (X,Y,Z) = (',f8.3,',',f8.3,',',f8.3,') [mm]')

      deallocate(table)

   enddo

end subroutine max_vonmises

!------------------------------------------------------------------------------------------------------------

subroutine test_results(iwhe, ldebug)
!--function: get the CONTACT results and print to the screen
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: iwhe           ! result element ID
   integer,      intent(in) :: ldebug         ! level of output requested
!--local variables:
   integer               :: icp, ncp, mx, my, lenarr
   integer               :: i, ixmax, iymax, izmax
   real(kind=8)          :: rvalues(20), carea, harea, sarea, pnmax, ptmax, tcpu, twall
   real(kind=8), dimension(:), allocatable :: pn, px, py, sx, sy
   real(kind=8), parameter :: pi     = 4d0*atan(1d0)
!--external functions used:
#include "contact_addon.ifc"
!--statement functions:
   integer               :: i2ix, i2iy, iidum, mxdum
   i2ix(iidum,mxdum) = mod(iidum-1,mxdum)+1
   i2iy(iidum,mxdum) = (iidum-1)/mxdum+1

   icp = -1

   ! get/print the total forces

   call cntc_getGlobalForces(iwhe, icp, 6, rvalues)
   if (ldebug.ge.1) write(*,'(2x,3(a,f8.1),a)') 'Total forces: Fx=',rvalues(1),', Fy=',rvalues(2),      &
                ', Fz=',rvalues(3),' [N]'

   ! get/print number of contact patches

   call cntc_getNumContactPatches(iwhe, ncp)
   if (ldebug.ge.1) then
      if (ncp.eq.1) then
         write(*,'(2x,a)') 'There is 1 contact patch'
      else
         write(*,'(2x,a,i3,a)') 'There are',ncp,' contact patches'
      endif
   endif

   ! get/print locations of the contact patches

   do icp = 1, ncp
      if (ldebug.ge.2) write(*,'(2x,a,i3,a)') 'Patch',icp,':'

      call cntc_getContactLocation(iwhe, icp, 20, rvalues)
      if (ldebug.ge.2) write(*,'(5x,4(a,f8.2),a)') 'X=',rvalues(1),', Y=',rvalues(2),', Z=',rvalues(3), &
                ' [mm], Delta=',rvalues(4)*180d0/pi,' [deg]'

      call cntc_getNumElements(iwhe, icp, mx, my)
      if (ldebug.ge.2) write(*,'(5x,2(a,i4))') 'Mx=',mx,', My=',my

      lenarr = mx * my

      call cntc_getContactPatchAreas(iwhe, icp, carea, harea, sarea)
      if (ldebug.ge.2) write(*,'(5x,3(a,f7.2),a)') 'C =',carea,', H =',harea,', S =',sarea,' [mm2]'

      call cntc_getMaximumPressure(iwhe, icp, pnmax)
      call cntc_getMaximumTraction(iwhe, icp, ptmax)
      if (ldebug.ge.2) write(*,'(5x,2(a,f8.1),a)') 'Pnmax=',pnmax,', Ptmax=',ptmax,' [N/mm2]'

      ! get the surface tractions, print maximum values

      allocate(pn(lenarr), px(lenarr), py(lenarr))
      call cntc_getTractions(iwhe, icp, lenarr, pn, px, py)
      ixmax = 1
      iymax = 1
      izmax = 1
      do i = 1, lenarr
         if (abs(pn(i)).gt.abs(pn(izmax))) izmax = i
         if (abs(px(i)).gt.abs(px(ixmax))) ixmax = i
         if (abs(py(i)).gt.abs(py(iymax))) iymax = i
      enddo

      if (ldebug.ge.2) write(*,'(5x, 3(2(a,i3),a,f8.1), a)')                            &
         'Max: Pn(',i2ix(izmax,mx),',',i2iy(izmax,mx),')=',pn(izmax),                   &
            ', Px(',i2ix(ixmax,mx),',',i2iy(ixmax,mx),')=',px(ixmax),                   &
            ', Py(',i2ix(iymax,mx),',',i2iy(iymax,mx),')=',py(iymax),' [N/mm2]'
      deallocate(pn, px, py)

      ! get the micro-slip, print maximum values

      allocate(sx(lenarr), sy(lenarr))
      call cntc_getMicroSlip(iwhe, icp, lenarr, sx, sy)

      ixmax = 1
      iymax = 1
      do i = 1, lenarr
         if (abs(sx(i)).gt.abs(sx(ixmax))) ixmax = i
         if (abs(sy(i)).gt.abs(sy(iymax))) iymax = i
      enddo
      if (ldebug.ge.2) write(*,'(5x, 2(2(a,i3),a,f8.5), a)')                            &
         'Max: Sx(',i2ix(ixmax,mx),',',i2iy(ixmax,mx),')=',sx(ixmax),                   &
            ', Sy(',i2ix(iymax,mx),',',i2iy(iymax,mx),')=',sy(iymax),' [-]'

      deallocate(sx, sy)

      ! get cpu/wall time and print

      call cntc_getCalculationTime(iwhe, icp, tcpu, twall)
      call cntc_resetCalculationTime(iwhe, icp)
      if (ldebug.ge.2) write(*,'(5x,2(a,f6.3))') 'CPU=',tcpu,', Wall=',twall

   enddo ! icp

end subroutine test_results

!============================================================================================================

