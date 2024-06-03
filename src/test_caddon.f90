!============================================================================================================
! Copyright 2008-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!============================================================================================================
! test_caddon - This program shows the use of the CONTACT library for Matlab, Python, Fortran and C,
!               using the basic strategy ("module 3") where the user builds the undeformed distance himself.
!
! For multiple "Result-Elements" (RE), and also multiple "Contact-Patches" (CP) per RE, 
!   1) the CONTACT input-data may be configured, 
!   2) the CONTACT calculation may be performed,
!   3) the CONTACT results may be retrieved.
!
! When using a non-Hertzian option, the input data for one contact patch consist of:
!   a) the material parameters for wheel and rail,
!   b) the static and kinetic friction coefficients,
!   c) the lower-left corner of the potential contact area
!      - pot.con == rectangular region encompassing the actual contact area
!   d) the number of discretization elements to use,
!   e) the size of the discretization elements,
!   f) the geometry of the bodies in undeformed state
!      - via the undeformed distance function, bodies touching at a point
!   g) the approach of the bodies' centers in the deformed configuration
!   h) the relative motion in tangential direction, ie. the creepages
!
! When using the Hertzian option, parts (c), (e) and (f) are not used. Instead, the user then specifies:
!   i) the principal curvatures or the semi-axes of the contact ellipse
! In this case, the program computes the geometry (f) and the total normal force from the Hertzian
! problem given by (i) and (g), and determines a suitable potential contact area (c) and the size of
! the discretization elements (e) from this Hertzian solution and the number of elements (d).
!
! For each RE the initialize function must be called once. For each "case" (timestep) per CP per RE the
! set-functions are used first, then the calculate function, and then the get-functions. The first
! set-function creates the data storage for a contact patch, and initializes it for the new case. Then
! all set-functions can be called, the calculate-function is called and the get-functions can be used.
! When next a set-function is called, this is taken as the start of a new case, where a different
! configuration can be used (approach, creepages, or also using different materials or geometry). After
! starting a new case, the results of the previous case cannot be retrieved (reliably) anymore.
!
! When all calculations for all CP's of a RE are done, the finalize function should be called. This
! cleans up the data for the RE. And when called for the last RE, it closes files and cleans up all
! other data of the CONTACT add-on as well.
!============================================================================================================
program test_caddon
   use, intrinsic        :: iso_c_binding, only: C_CHAR, C_NULL_CHAR
   implicit none
   include 'caddon_flags.inc'
!--local variables:
   integer               :: ire, icp, imodul
   integer               :: ldebug
   integer               :: ifcver, ierror, ioutput
   integer, parameter    :: mxflgs = CNTC_FLAGS_DIM
   integer               :: flags(mxflgs), values(mxflgs)
   real(kind=8)          :: g1, nu1, g3, laythk, tau_c0, k_tau, mu, a, b, rvalues(mxflgs)
   real(kind=8)          :: veloc, chi, dq, xl, yl, dx, dy, scale
   real(kind=8)          :: a1, b1, xi, yi
   real(kind=8)          :: pen, fn, vx, vy, phi
   integer               :: mx, my, lenarr, len_string
   integer               :: i, ix, iy, ipotcn, mdigit
   integer               :: nunits, iunits(5)
   real(kind=8), dimension(:), allocatable :: h
   character(len=30)             :: f_string
   character(len=30,kind=C_CHAR) :: c_string
!--external functions used (CONTACT dll):
#include "contact_addon.ifc"

   ! configure Fortran logical file units available to CONTACT

   nunits = 5
   iunits = (/ 11, 12, 13, 14, 15 /)
   call cntc_setFileUnits(nunits, iunits)

   ldebug = 2

   ! start-up, initalize internal data structures

   ifcver = -1
   ierror = -1
   if (.true.) then

      ! Initialize once using cntc_initializeFirst

      ioutput = 0       ! 0 = out-file, 1 = out-file + screen
      c_string   = ' ' // C_NULL_CHAR
      call cntc_initializeFirst_new(ifcver, ierror, ioutput, c_string, c_string, c_string, 1, 1, 1)

   elseif (.false.) then

      ! Initialize two Result Elements separately

      do ire = 1, 2
         imodul     = 3
         !!! f_string   = 'tmp/'
         f_string   = ' '
         c_string   =     trim(f_string)  // C_NULL_CHAR
         len_string = len(trim(f_string))
         call cntc_initialize(ire, imodul, ifcver, ierror, c_string, len_string)
      enddo

   else

      ! Rely on automatic initialization

   endif

   ! Set amount of output of C.library itself; enable or disable use of OpenMP

   flags(1) = CNTC_if_idebug ; values(1) = 2
   flags(2) = CNTC_if_openmp ; values(2) = 1
   call cntc_setGlobalFlags(2, flags, values)

   icp  = 2

   ! Next set all the values that are the same for both Result Elements

   do ire = 1, 2
      flags = 0 ; values = 0
      flags(1) = CNTC_if_units  ; values(1) = CNTC_un_spck      ! SI units
      flags(2) = CNTC_if_wrtinp ; values(2) = 1 !ire-1 ! off for ire=1, on for ire=2
      flags(3) = CNTC_ic_matfil ; values(3) = 0
      flags(4) = CNTC_ic_output ; values(4) = 3
      flags(5) = CNTC_ic_flow   ; values(5) = 4

      call cntc_setFlags(ire, icp, mxflgs, flags, values)

      flags = 0 ; rvalues = 0d0
      flags(1) = CNTC_mt_tim ; rvalues(1) = 0.1 + real(ire)
      flags(2) = CNTC_mt_xw  ; rvalues(2) = 12d0
      flags(3) = CNTC_mt_yw  ; rvalues(3) =  0.050d0
      flags(4) = CNTC_mt_yr  ; rvalues(4) = -0.012d0
      flags(5) = CNTC_mt_xr  ; rvalues(5) =  6d0

      call cntc_setMetadata(ire, icp, mxflgs, flags, rvalues)

      veloc = 24d0
      call cntc_setReferenceVelocity(ire, icp, veloc)

      chi   = 0d0
      dq    = 0.0005d0
      call cntc_setRollingStepsize(ire, icp, chi, dq)

      vx  = 0.001d0
      vy  = 0d0
      phi = 0.4d0
      call cntc_setCreepages(ire, icp, vx, vy, phi)
   end do

   ! Set additional values for Result Element 1
   ! First case: using Hertzian option with penetration and semi-axes prescribed

   ire = 1

   ! Elastic material model
   mdigit = 0
   g1     = 82000d6
   nu1    = 0.28d0
   rvalues(1:4) = (/ nu1, nu1, g1, g1 /)
   call cntc_setMaterialParameters(ire, icp, mdigit, 4, rvalues)

   ! Coulomb friction
   mu  = 0.3d0
   rvalues = 0d0
   rvalues(1) = mu
   rvalues(2) = mu
   call cntc_setFrictionMethod(ire, icp, 0, 2, rvalues)

   pen = 0.00005d0      ! PEN prescribed
   call cntc_setPenetration(ire, icp, pen)

   ipotcn = -3          ! Hertz, semi-axes prescribed
   mx  = 21             ! 17 elements + 4 guardband
   my  = 33
   a   = 0.00425d0      ! width 8.5 mm --> 17 elements of 0.5 mm
   b   = 0.007d0        ! height 14 mm ~~ 28 elements
   scale = real(mx) / real(max(1, mx-4))

   rvalues(1) = real(mx)
   rvalues(2) = real(my)
   rvalues(3) = a
   rvalues(4) = b
   rvalues(5) = scale
   call cntc_setHertzContact(ire, icp, ipotcn, 5, rvalues)

                        ! G=0: set maxgs, maxnr, maxin, maxout, eps
   values(1:4) = (/ 999, 30, 30, 1 /)
   rvalues(1)  = 0.0001
   call cntc_setSolverFlags(ire, icp, 0, 4, values, 1, rvalues)

   ! Set additional values for Result Element 2
   ! Second case, using non-Hertzian potential contact/undeformed distance specification

   ire = 2

   ! exponential falling friction
   rvalues = 0d0
   rvalues(1) = 0.30d0
   rvalues(2) = 0.60d0
   rvalues(3) = log(2d0)/1.25d0
   call cntc_setFrictionMethod(ire, icp, 5, 3, rvalues)

   ! Elastic material model with elastic third body layer
   mdigit = 4
   g1     = 82000d6
   nu1    = 0.28d0
   g3     = 100d6
   laythk = 0.00001
   tau_c0 = 0d0
   k_tau  = 0d0
   rvalues(1:8) = (/ nu1, nu1, g1, g1, g3, laythk, tau_c0, k_tau /)
   call cntc_setMaterialParameters(ire, icp, mdigit, 8, rvalues)

   fn  = 42083d0                 ! FN prescribed
   call cntc_setNormalForce(ire, icp, fn)

   mx  = 23
   my  = 33
   dx  =  0.0005d0
   dy  =  0.0005d0
   xl  = -0.00725d0
   yl  = -0.00725d0

   ipotcn = 1;
   rvalues(1) = mx
   rvalues(2) = my
   rvalues(3) = xl
   rvalues(4) = yl
   rvalues(5) = dx
   rvalues(6) = dy
   call cntc_setPotContact(ire, icp, ipotcn, 6, rvalues)

   a1  = 0.001214  * 1d3        ! from Hertzian case
   b1  = 0.0005731 * 1d3        ! from Hertzian case
   lenarr = mx * my
   allocate(h(lenarr))
   do iy = 1, my
      do ix = 1, mx
         i = ix + (iy-1)*mx
         xi = xl + (ix-0.5)*dx
         yi = yl + (iy-0.5)*dy
         h(i) = a1 * xi**2 + b1 * yi**2
      enddo
   enddo
   call cntc_setUndeformedDistc(ire, icp, 9, lenarr, h)
   deallocate(h)

   ! Calculate and print results for both Result Elements

   do ire = 1, 2

      if (ldebug.ge.1)  &
         write(*,'(/,a,i3)') '  ...start processing of Result Element', ire

      call cntc_calculate(ire, icp, ierror)

      if (ierror.eq.CNTC_err_allow) then
         write(*,'(a,i5,a)') '      ERROR: no valid license found'
         stop
      elseif (ierror.lt.0) then
         write(*,'(a,i5,a)') '      ERROR: code ierror=',ierror,'.'
      elseif (ierror.gt.0) then
         write(*,'(a,i5,a)') '      ERROR: contact extends outside bounding box,',ierror,               &
                                                                                ' elements at boundary'
      end if

      if (ire.eq.1) then
         call cntc_getGridDiscretization(ire, icp, dx, dy)
         if (ldebug.ge.2)       &
            write(*,'(2(a,f8.5),a)') '     Hertzian case: Dx=',dx,', Dy=',dy,' [m]'
      end if

      call test_results(ire, icp, ldebug)
   end do

   ! Calculate and print results again for Result Element 2 with different Fn, Cksi, Ceta, Cphi

   ire = 2
   if (ldebug.ge.1) write(*,'(/,a,i3,a)') '  ...start processing of Result Element', ire,               &
                ' with different Fn, Cksi etc.'

   fn  = 19563d0
   call cntc_setNormalForce(ire, icp, fn)
   vx  =  0.001d0
   vy  =  0d0
   phi = -0.3d0
   call cntc_setCreepages(ire, icp, vx, vy, phi)
   call cntc_calculate(ire, icp, ierror)
   call test_results(ire, icp, ldebug)

   call cntc_finalize(ire)

   ! Calculate and print results again for Result Element 1 with original Pen, Cksi, Ceta, Cphi

   ire = 1
   if (ldebug.ge.1) write(*,'(/,a,i3,a)') '  ...start processing of Result Element', ire,               &
                ' with original Pen, Cksi etc.'
   call cntc_calculate(ire, icp, ierror)
   call test_results(ire, icp, ldebug)
   call cntc_finalize(ire)

end program test_caddon

!------------------------------------------------------------------------------------------------------------

subroutine test_results(ire, icp, ldebug)
!--function: get the CONTACT results and print to the screen
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire           ! result element ID
   integer,      intent(in) :: icp           ! contact patch ID
   integer,      intent(in) :: ldebug         ! level of output requested
!--local variables:
   integer               :: mx, my, lenarr
   integer               :: i, ixmax, iymax, izmax
   integer               :: iidum, mxdum
   real(kind=8)          :: pen, tn, tx, ty, mz, carea, harea, sarea, pnmax, ptmax
   real(kind=8), dimension(:), allocatable :: pn, px, py, sx, sy
   real(kind=8)          :: tcpu, twall
!--external functions used:
#include "contact_addon.ifc"
!--statement functions used:
   integer               :: i2ix, i2iy
   i2ix(iidum,mxdum) = mod(iidum-1,mxdum)+1
   i2iy(iidum,mxdum) = (iidum-1)/mxdum+1

   ! Note: the output arguments of cntc-routines are actual values instead of scaled by G, Fstat

   ! get/print penetration for the contact patch

   call cntc_getPenetration(ire, icp, pen)
   if (ldebug.ge.2) write(*,'(5x,a,f9.6,a)') 'Pen=',pen,' [m]'

   ! get/print total forces for the contact patch

   call cntc_getContactForces(ire, icp, tn, tx, ty, mz)
   if (ldebug.ge.1) write(*,'(5x,4(a,f8.1),a)') 'Tn=',tn,', Fx=',tx,', Fy=',ty,' [N], Mz=',mz,' [N m]'

   ! get/print the actual number of elements used

   call cntc_getNumElements(ire, icp, mx, my)
   if (ldebug.ge.2) write(*,'(5x,2(a,i4))') 'Mx=',mx,', My=',my

   lenarr = mx * my

   call cntc_getContactPatchAreas(ire, icp, carea, harea, sarea)
   if (ldebug.ge.2) write(*,'(5x,3(a,f9.6),a)') 'C=',carea,' [m2], H=',harea,' [m2], S=',sarea,' [m2]'

   call cntc_getMaximumPressure(ire, icp, pnmax)
   call cntc_getMaximumTraction(ire, icp, ptmax)
   if (ldebug.ge.2) write(*,'(5x,2(a,f8.1),a)') 'Pnmax=',pnmax/1d6,', Ptmax=',ptmax/1d6,' [MN/m2]'

   ! get the surface tractions, print maximum values

   allocate(pn(lenarr), px(lenarr), py(lenarr))
   call cntc_getTractions(ire, icp, lenarr, pn, px, py)
   ixmax = 1
   iymax = 1
   izmax = 1
   do i = 1, lenarr
      if (abs(pn(i)).gt.abs(pn(izmax))) izmax = i
      if (abs(px(i)).gt.abs(px(ixmax))) ixmax = i
      if (abs(py(i)).gt.abs(py(iymax))) iymax = i
   enddo

   if (ldebug.ge.2) write(*,'(3(5x,2(a,i3),a,f8.1,a))')                           &
      'Max: Pn(',i2ix(izmax,mx),',',i2iy(izmax,mx),')=',pn(izmax)/1d6,' [MN/m2]', &
         ', Px(',i2ix(ixmax,mx),',',i2iy(ixmax,mx),')=',px(ixmax)/1d6,' [MN/m2]', &
         ', Py(',i2ix(iymax,mx),',',i2iy(iymax,mx),')=',py(iymax)/1d6,' [MN/m2]'
   deallocate(pn, px, py)

   ! get the micro-slip, print maximum values

   allocate(sx(lenarr), sy(lenarr))
   call cntc_getMicroSlip(ire, icp, lenarr, sx, sy)

   ixmax = 1
   iymax = 1
   do i = 1, lenarr
      if (abs(sx(i)).gt.abs(sx(ixmax))) ixmax = i
      if (abs(sy(i)).gt.abs(sy(iymax))) iymax = i
   enddo
   if (ldebug.ge.2) write(*,'(2(5x,2(a,i3),a,f8.5,a))')                         &
      'Max: Sx(',i2ix(ixmax,mx),',',i2iy(ixmax,mx),')=',sx(ixmax),' [-]',       &
         ', Sy(',i2ix(iymax,mx),',',i2iy(iymax,mx),')=',sy(iymax),' [-]'

   deallocate(sx, sy)

   ! get cpu/wall time and print

   call cntc_getCalculationTime(ire, icp, tcpu, twall)
   if (ldebug.ge.2) write(*,'(5x,2(a,f6.3))') 'CPU=',tcpu,', Wall=',twall

   if (ire.eq.1) call cntc_resetCalculationTime(ire, icp)

end subroutine test_results

!============================================================================================================

