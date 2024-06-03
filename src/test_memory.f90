!============================================================================================================
! test_memory - Test program for spotting potential memory leaks
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!============================================================================================================
program test_memory
   implicit none
   include 'caddon_flags.inc'
!--local variables:
   integer               :: ire, icp, itest, ntim, ldebug, itim, icount
   integer               :: flags(20), values(20)
   character(len=256)    :: tmpstr
!--external functions used:
#include "contact_addon.ifc"


   ! Usage:  test_memory [itest] [ntim] [ldebug]

   itest  = 1
   ntim   = 10000
   ldebug = 2

   icount = command_argument_count()
   if (icount.ge.1) then
      call get_command_argument(1, tmpstr)
      read(tmpstr, *) itest
   endif
   if (icount.ge.2) then
      call get_command_argument(2, tmpstr)
      read(tmpstr, *) ntim
   endif
   if (icount.ge.3) then
      call get_command_argument(3, tmpstr)
      read(tmpstr, *) ldebug
   endif

   write(*,'(a,i3,a,i7,a,i3)') 'Running test',itest,' for',ntim,' steps with ldebug=',ldebug

   ! Rely on automatic initialization
   ! Set amount of output of C.library itself; enable or disable use of OpenMP

   flags(1) = CNTC_if_idebug ; values(1) = 0    ! 1: little output
   flags(2) = CNTC_if_openmp ; values(2) = 1    ! 1: sequential run
   call cntc_setGlobalFlags(2, flags, values)

   ire  = 1
   icp  = 1

   if (itest.eq.1) then

      ! test 1: basic CONTACT, normal problem only, large pot.contact

      call test1_initialize(ire, icp)
      do itim = 1, ntim
         if (ldebug.ge.1) write(*,'(a,i7)') '  test 1: time step itim=', itim
         call test1_calculate(ire, icp, ldebug)
      enddo
      call cntc_finalize(ire)

   elseif (itest.eq.2) then

      ! test 2: basic CONTACT, steady rolling, large pot.contact

      call test2_initialize(ire, icp, ldebug)
      do itim = 1, ntim
         if (ldebug.ge.1) write(*,'(a,i7)') '  test 2: time step itim=', itim
         call test2_calculate(ire, icp, ldebug)
      enddo
      call cntc_finalize(ire)

   elseif (itest.eq.3) then

      ! test 3: w/r CONTACT, mbench case 1 left wheel

      call test3_initialize(ire)
      do itim = 1, ntim
         if (ldebug.ge.1) write(*,'(a,i7)') '  test 3: time step itim=', itim
         call test3_calculate(ire, ldebug)
      enddo
      call cntc_finalize(ire)

   elseif (itest.eq.4) then

      ! test 4: w/r CONTACT, Chalmers wheel flat

      call test4_initialize(ire)
      do itim = 1, ntim
         if (ldebug.ge.1) write(*,'(a,i7)') '  test 4: time step itim=', itim
         call test4_calculate(ire, ldebug)
      enddo
      call cntc_finalize(ire)

   endif

end program test_memory

!------------------------------------------------------------------------------------------------------------ 

subroutine test1_initialize(ire, icp)
!--test: basic CONTACT, normal problem only, large pot.contact
   use, intrinsic        :: iso_c_binding, only: C_CHAR, C_NULL_CHAR
   implicit none
   include 'caddon_flags.inc'
!--subroutine arguments:
   integer                       :: ire, icp
!--local variables:
   integer, parameter            :: mxflgs = 30
   integer                       :: imodul, iver, ierr, len_string, flags(mxflgs), values(mxflgs)
   integer                       :: ipotcn, mx, my, ibase, mdigit
   real(kind=8)                  :: g1, nu1, xl, yl, dx, dy, pen, rvalues(mxflgs)
   character(len=30)             :: f_string
   character(len=30,kind=C_CHAR) :: c_string
!--external functions used:
#include "contact_addon.ifc"

   ! initialize result element

   imodul     = 3
   !!! f_string   = 'tmp/'
   f_string   = ' '
   c_string   =     trim(f_string)  // C_NULL_CHAR
   len_string = len(trim(f_string))
   call cntc_initialize(ire, imodul, iver, ierr, c_string, len_string)

   ! set flags

   flags = 0 ; values = 0
   flags(1) = CNTC_if_units  ; values(1) = CNTC_un_cntc
   flags(2) = CNTC_ic_matfil ; values(2) = 0
   flags(3) = CNTC_ic_output ; values(3) = 1
   flags(4) = CNTC_ic_flow   ; values(4) = 1
   call cntc_setFlags(ire, icp, mxflgs, flags, values)

   ! set fixed problem data

   mdigit = 0
   g1     = 82000d0
   nu1    = 0.28d0
   rvalues(1:4) = (/ nu1, nu1, g1, g1 /)
   call cntc_setMaterialParameters(ire, icp, mdigit, 4, rvalues)

   rvalues = 0.3d0
   call cntc_setFrictionMethod(ire, icp, 0, 2, rvalues)

   ! using non-Hertzian potential contact/undeformed distance specification

   pen = 0.001d0
   call cntc_setPenetration(ire, icp, pen)

   mx  = 81
   my  = 41
   dx  =  0.5d0
   dy  =  0.5d0
   xl  = -20.25d0
   yl  = -10.25d0

   ipotcn = 1
   rvalues(1:6) = (/ dble(mx), dble(my), xl, yl, dx, dy /)
   call cntc_setPotContact(ire, icp, ipotcn, 6, rvalues)

   ibase = 1
   rvalues(1:6) = (/ 0.001214d0, 0d0, 0.0005731d0, 0d0, 0d0, 0d0 /)
   call cntc_setUndeformedDistc(ire, icp, ibase, 6, rvalues)

end subroutine test1_initialize

!------------------------------------------------------------------------------------------------------------ 

subroutine test1_calculate(ire, icp, ldebug)
!--test 1: basic CONTACT, normal problem only, large pot.contact
   implicit none
   include 'caddon_flags.inc'
!--subroutine arguments:
   integer                      :: ire, icp, ldebug
!--local variables:
!--local variables:
   integer                      :: ierr, mx, my
   real(kind=8)                 :: pen, tn, tx, ty, mz, carea, harea, sarea, pnmax, ptmax
!--external functions used:
#include "contact_addon.ifc"

   if (ldebug.ge.4)  &
      write(*,'(/,a,i3)') '  ...start processing of Result Element', ire

   call cntc_calculate(ire, icp, ierr)

   if (ierr.eq.CNTC_err_allow) then
      write(*,'(a,i5,a)') '      ERROR: no valid license found'
      stop
   elseif (ierr.lt.0) then
      write(*,'(a,i5,a)') '      ERROR: code ierr=',ierr,'.'
   elseif (ierr.gt.0) then
      write(*,'(a,i5,a)') '      ERROR: contact extends outside bounding box,',ierr,' elements at boundary'
   end if

   ! get/print penetration for the contact patch

   call cntc_getPenetration(ire, icp, pen)
   if (ldebug.ge.3) write(*,'(5x,a,f9.6,a)') 'Pen=',pen,' [mm]'

   ! get/print total forces for the contact patch

   call cntc_getContactForces(ire, icp, tn, tx, ty, mz)
   if (ldebug.ge.3) write(*,'(5x,4(a,f8.1),a)') 'Tn=',tn,', Fx=',tx,', Fy=',ty,' [N], Mz=',mz,' [N m]'

   ! get/print the actual number of elements used

   call cntc_getNumElements(ire, icp, mx, my)
   if (ldebug.ge.4) write(*,'(5x,2(a,i4))') 'Mx=',mx,', My=',my

   call cntc_getContactPatchAreas(ire, icp, carea, harea, sarea)
   if (ldebug.ge.3) write(*,'(5x,3(a,f8.2),a)') 'C=',carea,' [mm2], H=',harea,' [mm2], S=',sarea,' [mm2]'

   call cntc_getMaximumPressure(ire, icp, pnmax)
   call cntc_getMaximumTraction(ire, icp, ptmax)
   if (ldebug.ge.3) write(*,'(5x,2(a,f8.1),a)') 'Pnmax=',pnmax,', Ptmax=',ptmax,' [N/mm2]'

end subroutine test1_calculate

!------------------------------------------------------------------------------------------------------------ 

subroutine test2_initialize(ire, icp, ldebug)
!--test 2: basic CONTACT, steady rolling, large pot.contact
   implicit none
   include 'caddon_flags.inc'
!--subroutine arguments:
   integer                       :: ire, icp, ldebug
!--local variables:
   integer, parameter            :: mxflgs = 30
   integer                       :: mdigit
   real(kind=8)                  :: g1, nu1, pen, veloc, vx, vy, phi, rvalues(mxflgs)
!--external functions used:
#include "contact_addon.ifc"

   ! copy test 1 settings for normal problem

   call test1_initialize(ire, icp)
   if (ldebug.ge.100) write(*,*) 'test2_initialize...'

   ! test 2 additional settings

   pen = 0.001d0
   call cntc_setPenetration(ire, icp, pen)

   veloc = 24000d0
   call cntc_setReferenceVelocity(ire, icp, veloc)

   vx  = 0.001d0
   vy  = 0d0
   phi = 0.0004d0
   call cntc_setCreepages(ire, icp, vx, vy, phi)

   ! exponential falling friction cf. Polach

   rvalues(1:3) = (/ 0.30d0, 0.60d0, log(2d0)/0.00125d0 /)
   call cntc_setFrictionMethod(ire, icp, 5, 3, rvalues)

   ! interfacial layer M=4

   mdigit = 4
   g1     = 82000d0
   nu1    = 0.28d0
   rvalues(1:8) = (/ nu1, nu1, g1, g1, 100d0, 0.01d0, 0d0, 0d0 /)
   call cntc_setMaterialParameters(ire, icp, mdigit, 8, rvalues)

end subroutine test2_initialize

!------------------------------------------------------------------------------------------------------------ 

subroutine test2_calculate(ire, icp, ldebug)
!--test 2: basic CONTACT, steady rolling, large pot.contact
   implicit none
   include 'caddon_flags.inc'
!--subroutine arguments:
   integer                      :: ire, icp, ldebug
!--local variables:
!--local variables:
   real(kind=8)                 :: tn, tx, ty, mz, carea, harea, sarea, pnmax, ptmax
!--external functions used:
#include "contact_addon.ifc"

   ! copy test 1 calculation & results for normal problem

   call test1_calculate(ire, icp, ldebug)

   ! get/print total forces for the contact patch

   call cntc_getContactForces(ire, icp, tn, tx, ty, mz)
   if (ldebug.ge.1) write(*,'(5x,4(a,f8.1),a)') 'Tn=',tn,', Fx=',tx,', Fy=',ty,' [N], Mz=',mz,' [N mm]'

   call cntc_getContactPatchAreas(ire, icp, carea, harea, sarea)
   if (ldebug.ge.1) write(*,'(5x,3(a,f8.2),a)') 'C=',carea,' [mm2], H=',harea,' [mm2], S=',sarea,' [mm2]'

   call cntc_getMaximumPressure(ire, icp, pnmax)
   call cntc_getMaximumTraction(ire, icp, ptmax)
   if (ldebug.ge.1) write(*,'(5x,2(a,f8.1),a)') 'Pnmax=',pnmax,', Ptmax=',ptmax,' [N/mm2]'

end subroutine test2_calculate

!------------------------------------------------------------------------------------------------------------ 

subroutine test3_initialize(iwhe)
!--test 3: w/r CONTACT, mbench
   use, intrinsic        :: iso_c_binding, only: C_CHAR, C_NULL_CHAR
   implicit none
   include 'caddon_flags.inc'
!--subroutine arguments:
   integer                       :: iwhe
!--local variables:
   integer, parameter            :: mxflgs = 30
   integer                       :: imodul, icp, iver, ierr, mdigit
   integer                       :: len_string, flags(mxflgs), values(mxflgs)
   real(kind=8)                  :: gg, poiss, fstat, fkin, dx, ds, dqrel, rdum
   integer                       :: itype, iside, mirrory, ewheel, ztrack, iparam(3)
   real(kind=8)                  :: sclfac, smooth, rparam(2), rvalues(mxflgs)
   character(len=256)             :: f_string, f_fname
   character(len=256,kind=C_CHAR) :: c_string, c_fname
!--external functions used:
#include "contact_addon.ifc"

   ! initialize result element

   imodul     =  1
   icp        = -1
   !!! f_string   = 'tmp/'
   f_string   = ' '
   c_string   =     trim(f_string)  // C_NULL_CHAR
   len_string = len(trim(f_string))
   call cntc_initialize(iwhe, imodul, iver, ierr, c_string, len_string)

   ! set flags: configure control digits & output

   flags = 0 ; values = 0
   flags(1) = CNTC_if_units  ; values(1) = CNTC_un_cntc  ! CONTACT unit convention: mm
   flags(2) = CNTC_ic_config ; values(2) = 0             ! C1: 0=left wheel, 1=right
   flags(3) = CNTC_ic_tang   ; values(3) = 3             ! T=3: steady state rolling
   flags(4) = CNTC_ic_pvtime ; values(4) = 2             ! P=2: no previous time
   flags(5) = CNTC_ic_norm   ; values(5) = 1             ! N=1: vertical force prescribed
   flags(6) = CNTC_if_wrtinp ; values(6) = 0             !   0: no inp-file needed
   flags(7) = CNTC_ic_matfil ; values(7) = 0             !   0: no mat-files needed
   flags(8) = CNTC_ic_output ; values(8) = 0             ! O=1: little output to out-file
   flags(9) = CNTC_ic_flow   ; values(9) = 0             ! W=1: little flow trace
   call cntc_setFlags(iwhe, icp, mxflgs, flags, values)

   ! set material parameters

   mdigit = 0
   gg     = 82000d0  ! CONTACT units: [N/mm2]
   poiss  = 0.28d0   ! [-]
   rvalues(1:4) = (/ poiss, poiss, gg, gg /)
   call cntc_setMaterialParameters(iwhe, icp, mdigit, 4, rvalues)

   ! set friction parameters: L=0, Coulomb friction

   fstat   = 0.3d0  ! [-]
   fkin    = fstat  ! [-]
   rvalues(1:2) = (/ fstat, fkin /)
   call cntc_setFrictionMethod(iwhe, icp, 0, 2, rvalues)

   ! set grid discretization

   dx     = 0.5d0   ! [mm]
   ds     = 0.5d0   ! [mm]
   rvalues(1:3) = (/ dx, ds, 0d0 /)
   call cntc_setPotContact(iwhe, icp, 0, 3, rvalues)

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
   len_string = len(trim(f_fname))

   itype   = -1      ! using filename extension
   iside   =  2      ! left and right rails
   mirrory =  0      ! no mirroring
   sclfac  = 1d0     ! already in [mm], no scaling
   smooth  = 0d0     ! no smoothing
   iparam  = (/ itype, iside, mirrory /)
   rparam  = (/ sclfac, smooth /)

   call cntc_setProfileInputFname(iwhe, c_fname, len_string, 3, iparam, 2, rparam)

   ! set wheelset dimensions

   ewheel = 3
   rvalues(1:3) = (/ 1360d0, -70d0, 460d0 /)
   call cntc_setWheelsetDimensions(iwhe, ewheel, 3, rvalues)

   ! set wheel profile. If two profiles are given, first is for left side, second for right side

   f_fname   = 'profiles/MBench_S1002_v3.prw'
   c_fname   =     trim(f_fname)  // C_NULL_CHAR
   len_string = len(trim(f_fname))

   itype   = -1      ! using filename extension
   iside   =  2      ! left and right wheels
   mirrory =  0      ! no mirroring
   sclfac  = 1d0     ! already in [mm], no scaling
   smooth  = 0d0     ! no smoothing
   iparam  = (/ itype, iside, mirrory /)
   rparam  = (/ sclfac, smooth /)

   call cntc_setProfileInputFname(iwhe, c_fname, len_string, 3, iparam, 2, rparam)

end subroutine test3_initialize

!------------------------------------------------------------------------------------------------------------ 

subroutine test3_calculate(iwhe, ldebug)
!--test 3: w/r CONTACT, mbench
   implicit none
   include 'caddon_flags.inc'
!--subroutine arguments:
   integer                      :: iwhe, ldebug
!--local variables:
   integer                      :: icp, ierr
   real(kind=8)                 :: y_ws, yaw_ws, roll_ws, vpitch, fz, rvalues(30)
!--external functions used:
#include "contact_addon.ifc"

   icp     = -1

   y_ws    = 0d0
   yaw_ws  = 0d0
   roll_ws = 0d0
   vpitch  = -4.34811845d0

   ! set wheelset position & velocity

   rvalues(1:6) = (/ 0d0, y_ws, 0d0, roll_ws, yaw_ws, 0d0 /)
   call cntc_setWheelsetPosition(iwhe, 1, 6, rvalues)

   rvalues(1:6) = (/ 2000d0, 0d0, 0d0, 0d0, 0d0, vpitch /)
   call cntc_setWheelsetVelocity(iwhe, 1, 6, rvalues)

   ! set vertical force

   fz  = 10000d0    ! [N]
   call cntc_setVerticalForce(iwhe, fz)

   if (ldebug.ge.4)  &
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

end subroutine test3_calculate

!------------------------------------------------------------------------------------------------------------ 

subroutine test4_initialize(iwhe)
!--test 4: w/r CONTACT, chalmers wheel flat
   use, intrinsic        :: iso_c_binding, only: C_CHAR, C_NULL_CHAR
   implicit none
   include 'caddon_flags.inc'
!--subroutine arguments:
   integer                       :: iwhe
!--local variables:
   real(kind=8), parameter       :: pi     = 4d0*atan(1d0)
   integer, parameter            :: mxflgs = 30
   integer                       :: imodul, icp, iver, ierr, mdigit, ldigit, ipotcn
   integer                       :: len_string, flags(mxflgs), values(mxflgs)
   real(kind=8)                  :: gg, poiss, fstat, fkin, dx, ds, dqrel, rdum,                        &
                                    a_sep, d_sep, d_comb, d_turn
   integer                       :: itype, mirrory, mirrorz, ewheel, ztrack, len_fname, iparam(3)
   real(kind=8)                  :: sclfac, smooth, rparam(2), rvalues(mxflgs)
   character(len=256)             :: f_string, f_fname
   character(len=256,kind=C_CHAR) :: c_string, c_fname
!--external functions used:
#include "contact_addon.ifc"

   ! initialize result element

   imodul     =  1
   icp        = -1
   !!! f_string   = 'tmp/'
   f_string   = ' '
   c_string   =     trim(f_string)  // C_NULL_CHAR
   len_string = len(trim(f_string))
   call cntc_initialize(iwhe, imodul, iver, ierr, c_string, len_string)

   ! set flags: configure control digits & output

   flags = 0 ; values = 0
   flags(1) = CNTC_if_units  ; values(1) = CNTC_un_cntc  ! CONTACT unit convention: mm
   flags(2) = CNTC_ic_config ; values(2) = 1             ! C1: 0=left wheel, 1=right
   flags(3) = CNTC_ic_tang   ; values(3) = 3             ! T=3: steady state rolling
   flags(4) = CNTC_ic_pvtime ; values(4) = 2             ! P=2: no previous time
   flags(5) = CNTC_if_wrtinp ; values(5) = 0             !   0: no inp-file needed
   flags(6) = CNTC_ic_matfil ; values(6) = 0             !   0: no mat-files needed
   flags(7) = CNTC_ic_output ; values(7) = 0             ! O=1: little output to out-file
   flags(8) = CNTC_ic_flow   ; values(8) = 0             ! W=1: little flow trace
   call cntc_setFlags(iwhe, icp, mxflgs, flags, values)

   ! set material parameters

   mdigit = 0
   gg     = 82000d0  ! CONTACT units: [N/mm2]
   poiss  = 0.28d0   ! [-]
   rvalues(1:4) = (/ poiss, poiss, gg, gg /)
   call cntc_setMaterialParameters(iwhe, icp, mdigit, 4, rvalues)

   ! set friction parameters: L=0, Coulomb friction

   ldigit = 0
   fstat   = 0.3d0  ! [-]
   fkin    = fstat  ! [-]
   rvalues(1:2) = (/ fstat, fkin /)
   call cntc_setFrictionMethod(iwhe, icp, ldigit, 2, rvalues)

   ! set rolling step size: ratio c = dq / dx

   rdum   = 0d0     ! chi: ignored
   dqrel  = 1d0     ! [-]
   call cntc_setRollingStepsize(iwhe, icp, rdum, dqrel)

   ! set track dimensions & deviations

   ztrack = 3
   rvalues(1:11) = (/ 14d0, 0d0, 1435d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 /)
   call cntc_setTrackDimensions(iwhe, ztrack, 11, rvalues)

   ! set grid discretization

   ipotcn = -1
   dx     = 0.4d0   ! [mm]
   ds     = 0.4d0   ! [mm]
   a_sep  = pi/2d0  ! [rad]
   d_sep  =  8.0d0  ! [mm]
   d_comb =  4.0d0  ! [mm]
   d_turn = 12.0d0  ! [mm]

   rvalues(1:6) = (/ dx, ds, a_sep, d_sep, d_comb, d_turn /)
   call cntc_setPotContact(iwhe, icp, ipotcn, 6, rvalues)

   ! set constant rail profile UIC60

   f_fname   = 'profiles/uic60_true.prr'
   c_fname   =     trim(f_fname)  // C_NULL_CHAR
   len_fname = len(trim(f_fname))

   itype   = -1      ! using filename extension
   mirrory =  0      ! no mirroring
   mirrorz = -1      ! no mirroring
   sclfac  = 1d0     ! data in [mm], no scaling
   smooth  = 0d0     ! no smoothing
   iparam(1:4) = (/ itype, 0, mirrory, mirrorz /)
   rparam(1:2) = (/ sclfac, smooth /)

   call cntc_setProfileInputFname(iwhe, c_fname, len_fname, 4, iparam, 2, rparam)

   ! set wheelset dimensions

   rvalues(1:3) = (/ 1360d0, -70d0, 490d0 /)
   ewheel = 3
   call cntc_setWheelsetDimensions(iwhe, ewheel, 3, rvalues)

   ! set wheel profile

   f_fname   = 'profiles/S1002_flat.slcw'
   c_fname   =     trim(f_fname)  // C_NULL_CHAR
   len_fname = len(trim(f_fname))

   itype   = -1      ! using filename extension
   mirrory =  0      ! no y-mirroring
   mirrorz = -1      ! no z-mirroring
   sclfac  = 1d0     ! already in [mm], no scaling
   smooth  = 5d0     ! lambda smoothing
   iparam(1:4) = (/ itype, 0, mirrory, mirrorz /)
   rparam(1:2) = (/ sclfac, smooth /)

   call cntc_setProfileInputFname(iwhe, c_fname, len_fname, 4, iparam, 2, rparam)

end subroutine test4_initialize

!------------------------------------------------------------------------------------------------------------ 

subroutine test4_calculate(iwhe, ldebug)
!--test 4: w/r CONTACT, Chalmers wheel flat
   implicit none
   include 'caddon_flags.inc'
!--subroutine arguments:
   integer                      :: iwhe, ldebug
!--local variables:
   real(kind=8), parameter      :: pi    = 4d0*atan(1d0)
   integer,      parameter      :: ncase = 25
   integer                      :: icp, icase, ierr
   real(kind=8)                 :: s_ws(ncase), y_ws(ncase), fz_ws(ncase), pitch_ws(ncase),             &
                                   yaw_ws(ncase), roll_ws(ncase), vpitch(ncase), rvalues(30)
!--external functions used:
#include "contact_addon.ifc"

   icp     = -1

   s_ws(1:ncase)     =       0.0d0
   y_ws(1:ncase)     =       0.0d0
   fz_ws(1:ncase)    =     125.0d3
   pitch_ws(1:ncase) = (/ (-24.0d0 - 1.0d0 * icase, icase = 1, ncase) /) * pi/180d0
   yaw_ws(1:ncase)   =       0.0d0
   roll_ws(1:ncase)  =       0.0d0
   vpitch(1:ncase)   =      -4.08190679d0

   do icase = 1, ncase

      ! set wheelset position & velocity

      rvalues(1:6) = (/ s_ws(icase), y_ws(icase), 0d0, roll_ws(icase), yaw_ws(icase), pitch_ws(icase) /)
      call cntc_setWheelsetPosition(iwhe, 1, 6, rvalues)

      rvalues(1:6) = (/ 2000d0, 0d0, 0d0, 0d0, 0d0, vpitch(icase) /)
      call cntc_setWheelsetVelocity(iwhe, 1, 6, rvalues)

      ! set vertical force

      call cntc_setVerticalForce(iwhe, fz_ws(icase))

      if (ldebug.ge.4)  &
         write(*,'(/,2(a,i3))') '  ...start processing of Result Element', iwhe,' case',icase

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

   enddo ! icase

end subroutine test4_calculate

!------------------------------------------------------------------------------------------------------------

subroutine test_results(iwhe, ldebug)
!--function: get the CONTACT results and print to the screen
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: iwhe           ! result element ID
   integer,      intent(in) :: ldebug         ! level of output requested
!--local variables:
   integer               :: icp, ncp, mx, my
   real(kind=8)          :: rvalues(20), carea, harea, sarea, pnmax, ptmax
   real(kind=8), parameter :: pi     = 4d0*atan(1d0)
!--external functions used:
#include "contact_addon.ifc"

   icp = -1

   ! get/print the total forces

   call cntc_getGlobalForces(iwhe, icp, 6, rvalues)
   if (ldebug.ge.2) write(*,'(2x,3(a,f8.1),a)') 'Total forces: Fx=',rvalues(1),', Fy=',rvalues(2),      &
                ', Fz=',rvalues(3),' [N]'

   ! get/print number of contact patches

   call cntc_getNumContactPatches(iwhe, ncp)
   if (ldebug.ge.3) then
      if (ncp.eq.1) then
         write(*,'(2x,a)') 'There is 1 contact patch'
      else
         write(*,'(2x,a,i3,a)') 'There are',ncp,' contact patches'
      endif
   endif

   ! get/print locations of the contact patches

   do icp = 1, ncp
      if (ldebug.ge.3) write(*,'(2x,a,i3,a)') 'Patch',icp,':'

      call cntc_getContactLocation(iwhe, icp, 20, rvalues)
      if (ldebug.ge.3) write(*,'(5x,4(a,f8.2),a)') 'X=',rvalues(1),', Y=',rvalues(2),', Z=',rvalues(3), &
                ' [mm], Delta=',rvalues(4)*180d0/pi,' [deg]'

      call cntc_getNumElements(iwhe, icp, mx, my)
      if (ldebug.ge.3) write(*,'(5x,2(a,i4))') 'Mx=',mx,', My=',my

      call cntc_getContactPatchAreas(iwhe, icp, carea, harea, sarea)
      if (ldebug.ge.3) write(*,'(5x,3(a,f7.2),a)') 'C =',carea,', H =',harea,', S =',sarea,' [mm2]'

      call cntc_getMaximumPressure(iwhe, icp, pnmax)
      call cntc_getMaximumTraction(iwhe, icp, ptmax)
      if (ldebug.ge.3) write(*,'(5x,2(a,f8.1),a)') 'Pnmax=',pnmax,', Ptmax=',ptmax,' [N/mm2]'

   enddo ! icp

end subroutine test_results

!============================================================================================================

