!------------------------------------------------------------------------------------------------------------
! Test program for variable profiles: switches and crossings, wheel out-of-roundness
!
! Usage: test_varprof [expnam], with expnam = cross_brute, cross+wing, mbench_intrup, etc.
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
program test_varprof
   use, intrinsic          :: iso_c_binding, only: C_CHAR, C_NULL_CHAR
   implicit none
   include 'caddon_flags.inc'
!--local variables:
   integer,      parameter :: iwhe = 1, mxflgs = CNTC_FLAGS_DIM, mxcase = 100
   real(kind=8), parameter :: pi   = 4d0*atan(1d0)
   integer                        :: len_wrkdir, len_outdir, len_expnam, len_fname, icount, ldebug,     &
                                     ver, ierr, iout
   character(len=256)             :: f_wrkdir, f_outdir, f_expnam, f_fname, tmpstr
   character(len=256,kind=C_CHAR) :: c_wrkdir, c_outdir, c_expnam, c_fname
   integer                 :: icp, imodul, gdigit, mdigit, ldigit, itype, mirrory, mirrorz, irep, nrep
   integer                 :: flags(mxflgs), values(mxflgs)
   real(kind=8)            :: gg, poiss, fstat, fkin, dx, ds, a_sep, d_sep, d_comb, d_turn, dqrel, rdum
   integer                 :: ewheel, ztrack, ipotcn, icase, ncase, iparam(mxflgs)
   real(kind=8)            :: sclfac, smooth, rparam(2), rvalues(mxflgs)
   real(kind=8)            :: s_ws(mxcase), y_ws(mxcase), z_ws(mxcase), fz_ws(mxcase), pitch_ws(mxcase), &
                              roll_ws(mxcase), yaw_ws(mxcase), vpitch(mxcase)
!--external functions used:
#include "contact_addon.ifc"

   !---------------------------------------------------------------------------------------------------------
   ! Get experiment name
   !---------------------------------------------------------------------------------------------------------

   f_expnam = ' '
   icount = command_argument_count()
   if (icount.ge.1) then
      call get_command_argument(1, f_expnam)
   endif
   if (icount.ge.2) then
      call get_command_argument(2, tmpstr)
      read(tmpstr, *) nrep
   else
      nrep = 1
   endif

   if (trim(f_expnam).ne.'cross_brute'         .and. trim(f_expnam).ne.'cross_locus'   .and.            &
       trim(f_expnam).ne.'cross+wing'          .and. trim(f_expnam).ne.'cw_interrupt'  .and.            &
       trim(f_expnam).ne.'mbench_brute'        .and. trim(f_expnam).ne.'mbench_intrup' .and.            &
       trim(f_expnam).ne.'mbench_locus'        .and. trim(f_expnam).ne.'two_patches'   .and.            &
       trim(f_expnam).ne.'wing_brute'          .and. trim(f_expnam).ne.'wing_locus'    .and.            &
       trim(f_expnam).ne.'chalmers_flat_fz125' .and. trim(f_expnam).ne.'rounded_flat_d09') then
      write(*,*) ' Unknown experiment "',trim(f_expnam),'", using default "cross_brute".'
      f_expnam = 'cross_brute'
   endif

   !---------------------------------------------------------------------------------------------------------
   ! Initialization of the CONTACT library
   !---------------------------------------------------------------------------------------------------------

   ldebug = 1

   ! Initialize using cntc_initializeFirst, setting the experiment name

   ver         = -1
   ierr        = -1
   iout        =  0    ! enable/disable (1/0) screen output
   f_wrkdir    = ' '   ! effective working directory, starting point for relative paths
   c_wrkdir    =     trim(f_wrkdir)  // C_NULL_CHAR
   len_wrkdir  = len(trim(f_wrkdir))
   f_outdir    = ' '   ! could be '/tmp' for instance
   c_outdir    =     trim(f_outdir)  // C_NULL_CHAR
   len_outdir  = len(trim(f_outdir))
   c_expnam    =     trim(f_expnam )  // C_NULL_CHAR
   len_expnam  = len(trim(f_expnam ))
   call cntc_initializeFirst_new(ver, ierr, iout, c_wrkdir, c_outdir, c_expnam, len_wrkdir, len_outdir, &
                        len_expnam)

   ! Set the amount of output of the C.library itself; enable or disable use of OpenMP

   flags(1) = CNTC_if_idebug ; values(1) = 1
   flags(2) = CNTC_if_openmp ; values(2) = 1    ! max #threads
   call cntc_setGlobalFlags(2, flags, values)

   !---------------------------------------------------------------------------------------------------------
   ! Initialize one result element, using module 1 for wheel/rail contact
   !---------------------------------------------------------------------------------------------------------

   imodul = 1
   call cntc_initialize(iwhe, imodul, ver, ierr, c_outdir, len_outdir)

   ! Configure result element, setting permanent data

   icp = -1

   ! set flags: configure control digits & output, defaults, may be changed per experiment

   flags = 0 ; values = 0
   flags(1) = CNTC_if_units  ; values(1) = CNTC_un_cntc  ! CONTACT unit convention: mm
   flags(2) = CNTC_ic_config ; values(2) = 1             ! C1: 0=left wheel, 1=right
   flags(3) = CNTC_ic_pvtime ; values(3) = 2             ! P=2: no previous time
   flags(4) = CNTC_ic_tang   ; values(4) = 3             ! T=3: steady state rolling
   flags(5) = CNTC_if_wrtinp ; values(5) = 0             !   0: no inp-file needed
   flags(6) = CNTC_ic_matfil ; values(6) = 0             !   0: no mat-files needed
   flags(7) = CNTC_ic_output ; values(7) = 3             ! O=3: detailed output, no profiles
   flags(8) = CNTC_ic_flow   ; values(8) = 2             ! W=2: some flow trace
   flags(9) = CNTC_ic_npomax ; values(9) = 20000         ! max #elements in pot.contact

   call cntc_setFlags(iwhe, icp, mxflgs, flags, values)

                        ! G=0: set maxgs, maxin, maxnr, maxout, eps
   gdigit = 0
   values(1:4) = (/ 999, 100, 30, 1 /)
   rvalues(1)  = 1d-5
   call cntc_setSolverFlags(iwhe, icp, gdigit, 4, values, 1, rvalues)

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
   rvalues(1:11) = (/ -1d0, 750d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 /)
   call cntc_setTrackDimensions(iwhe, ztrack, 11, rvalues)

   !---------------------------------------------------------------------------------------------------------
   ! Further settings dependent on experiment name
   !---------------------------------------------------------------------------------------------------------

   if     (trim(f_expnam).eq.'cross_brute' .or. trim(f_expnam).eq.'cross_locus') then

      if (trim(f_expnam).eq.'cross_locus') then
         flags(1) = CNTC_ic_discns ; values(1) = 2             ! D=2: contact locus method (default)
      else
         flags(1) = CNTC_ic_discns ; values(1) = 5             ! D=5: brute force method
      endif
      call cntc_setFlags(iwhe, icp, 1, flags, values)

      ! set grid discretization

      ipotcn = -1
      dx     = 0.2d0   ! [mm]
      ds     = 0.05d0  ! [mm]
      a_sep  = pi/2d0  ! [rad]
      d_sep  =  8.0d0  ! [mm]
      d_comb =  4.0d0  ! [mm]
      d_turn = 12.0d0  ! [mm]

      rvalues(1:6) = (/ dx, ds, a_sep, d_sep, d_comb, d_turn /)
      call cntc_setPotContact(iwhe, icp, ipotcn, 6, rvalues)

      ! set variable rail profile for switch/crossing

      f_fname   = '../profiles/cross_nose.slcs'
      c_fname   =     trim(f_fname)  // C_NULL_CHAR
      len_fname = len(trim(f_fname))

      itype   = -1      ! using filename extension
      mirrory =  0      ! no mirroring
      mirrorz = -1      ! no mirroring
      sclfac  = 1d0     ! already in [mm], no scaling
      smooth  = 0d0     ! no smoothing
      iparam(1:4) = (/ itype, 0, mirrory, mirrorz /)
      rparam(1:2) = (/ sclfac, smooth /)

      call cntc_setProfileInputFname(iwhe, c_fname, len_fname, 4, iparam, 2, rparam)

      ncase = 1
      s_ws(1:ncase)     =     0.0d0
      y_ws(1:ncase)     =     0.0d0
      z_ws(1:ncase)     =    -2.1275d0
      fz_ws(1:ncase)    =    -1.0d0
      pitch_ws(1:ncase) =     0.0d0
      yaw_ws(1:ncase)   =     0.0d0
      roll_ws(1:ncase)  =     0.0d0
      vpitch(1:ncase)   =    -4.34811810d0

   elseif (trim(f_expnam).eq.'cross+wing') then

      flags(1) = CNTC_ic_discns ; values(1) = 2             ! D=2: contact locus method (default)
      call cntc_setFlags(iwhe, icp, 1, flags, values)

      ! set grid discretization

      ipotcn = -1
      dx     = 0.2d0   ! [mm]
      ds     = 0.05d0  ! [mm]
      a_sep  = pi/2d0  ! [rad]
      d_sep  =  8.0d0  ! [mm]
      d_comb =  4.0d0  ! [mm]
      d_turn = 12.0d0  ! [mm]

      rvalues(1:6) = (/ dx, ds, a_sep, d_sep, d_comb, d_turn /)
      call cntc_setPotContact(iwhe, icp, ipotcn, 6, rvalues)

      ! set variable rail profile for switch/crossing

      f_fname   = '../profiles/cross+wing.slcs'
      c_fname   =     trim(f_fname)  // C_NULL_CHAR
      len_fname = len(trim(f_fname))

      itype   = -1      ! using filename extension
      mirrory =  0      ! no mirroring
      mirrorz = -1      ! no mirroring
      sclfac  = 1d0     ! already in [mm], no scaling
      smooth  = 0d0     ! no smoothing
      iparam(1:4) = (/ itype, 0, mirrory, mirrorz /)
      rparam(1:2) = (/ sclfac, smooth /)

      call cntc_setProfileInputFname(iwhe, c_fname, len_fname, 4, iparam, 2, rparam)

      ncase = 1
      s_ws(1:ncase)     =   100.0d0
      y_ws(1:ncase)     =     0.0d0
      z_ws(1:ncase)     =     2.5212d0
      fz_ws(1:ncase)    =    -1.0d0
      pitch_ws(1:ncase) =     0.0d0
      yaw_ws(1:ncase)   =     0.0d0
      roll_ws(1:ncase)  =     0.0d0
      vpitch(1:ncase)   =    -4.34811810d0

   elseif (trim(f_expnam).eq.'cw_interrupt') then

      flags(1) = CNTC_ic_discns ; values(1) = 2             ! D=2: contact locus method (default)
      call cntc_setFlags(iwhe, icp, 1, flags, values)

      ! set grid discretization

      ipotcn = -1
      dx     = 0.2d0   ! [mm]
      ds     = 0.05d0  ! [mm]
      a_sep  = pi/2d0  ! [rad]
      d_sep  =  8.0d0  ! [mm]
      d_comb =  4.0d0  ! [mm]
      d_turn = 12.0d0  ! [mm]

      rvalues(1:6) = (/ dx, ds, a_sep, d_sep, d_comb, d_turn /)
      call cntc_setPotContact(iwhe, icp, ipotcn, 6, rvalues)

      ! set variable rail profile for switch/crossing

      f_fname   = '../profiles/cross+wing_extd.slcs'
      c_fname   =     trim(f_fname)  // C_NULL_CHAR
      len_fname = len(trim(f_fname))

      itype   = -1      ! using filename extension
      mirrory =  0      ! no mirroring
      mirrorz = -1      ! no mirroring
      sclfac  = 1d0     ! already in [mm], no scaling
      smooth  = 0d0     ! no smoothing
      iparam(1:4) = (/ itype, 0, mirrory, mirrorz /)
      rparam(1:2) = (/ sclfac, smooth /)

      call cntc_setProfileInputFname(iwhe, c_fname, len_fname, 4, iparam, 2, rparam)

      ncase = 2
      s_ws(1:ncase)     =       0.0d0
      y_ws(1:ncase)     = (/   40.0d0,         0.0d0   /)
      z_ws(1:ncase)     = (/   -1.6736d0,      2.5287d0  /)
      fz_ws(1:ncase)    = (/   10.0d3,        -1.0d0   /)
      pitch_ws(1:ncase) =       0.0d0
      yaw_ws(1:ncase)   =       0.0d0
      roll_ws(1:ncase)  =       0.0d0
      vpitch(1:ncase)   =      -4.34811810d0

   elseif (trim(f_expnam).eq.'mbench_brute' .or.  trim(f_expnam).eq.'mbench_locus') then

      if (trim(f_expnam).eq.'mbench_locus') then
         flags(1) = CNTC_ic_discns ; values(1) = 2          ! D=2: contact locus method (default)
      else
         flags(1) = CNTC_ic_discns ; values(1) = 5          ! D=5: brute force method
      endif
      call cntc_setFlags(iwhe, icp, 1, flags, values)

      ! set grid discretization

      ipotcn = -1
      dx     = 0.2d0   ! [mm]
      ds     = 0.2d0  ! [mm]
      a_sep  = pi/2d0  ! [rad]
      d_sep  =  8.0d0  ! [mm]
      d_comb =  4.0d0  ! [mm]
      d_turn = 12.0d0  ! [mm]

      rvalues(1:6) = (/ dx, ds, a_sep, d_sep, d_comb, d_turn /)
      call cntc_setPotContact(iwhe, icp, ipotcn, 6, rvalues)

      ! set variable rail profile for switch/crossing

      f_fname   = '../profiles/uk_crossing.slcs'
      c_fname   =     trim(f_fname)  // C_NULL_CHAR
      len_fname = len(trim(f_fname))

      itype   = -1      ! using filename extension
      mirrory =  0      ! no mirroring
      mirrorz = -1      ! no mirroring
      sclfac  = 1d3     ! data in [m], 1000x scaling
      smooth  = 0d0     ! no smoothing
      iparam(1:4) = (/ itype, 0, mirrory, mirrorz /)
      rparam(1:2) = (/ sclfac, smooth /)

      call cntc_setProfileInputFname(iwhe, c_fname, len_fname, 4, iparam, 2, rparam)

      ncase = 35
      s_ws(1:ncase)     =     239.0d0
      y_ws(1:ncase)     = (/ ( 0.80d0 + 0.01d0 * icase, icase = 1, ncase) /)
      z_ws(1:ncase)     =       0.0d0
      fz_ws(1:ncase)    =      10.0d3
      pitch_ws(1:ncase) =       0.0d0
      yaw_ws(1:ncase)   =       0.0d0
      roll_ws(1:ncase)  =       0.0d0
      vpitch(1:ncase)   =      -4.34811810d0

   elseif (trim(f_expnam).eq.'mbench_intrup') then

      flags(1) = CNTC_ic_discns ; values(1) = 2          ! D=2: contact locus method (default)
      call cntc_setFlags(iwhe, icp, 1, flags, values)

      ! set grid discretization

      ipotcn = -1
      dx     = 0.2d0   ! [mm]
      ds     = 0.2d0  ! [mm]
      a_sep  = pi/2d0  ! [rad]
      d_sep  =  8.0d0  ! [mm]
      d_comb =  4.0d0  ! [mm]
      d_turn = 12.0d0  ! [mm]

      rvalues(1:6) = (/ dx, ds, a_sep, d_sep, d_comb, d_turn /)
      call cntc_setPotContact(iwhe, icp, ipotcn, 6, rvalues)

      ! set variable rail profile for switch/crossing

      f_fname   = '../profiles/uk_interrupt_v2.slcs'
      c_fname   =     trim(f_fname)  // C_NULL_CHAR
      len_fname = len(trim(f_fname))

      itype   = -1      ! using filename extension
      mirrory =  0      ! no mirroring
      mirrorz = -1      ! no mirroring
      sclfac  = 1d3     ! data in [m], 1000x scaling
      smooth  = 0d0     ! no smoothing
      iparam(1:4) = (/ itype, 0, mirrory, mirrorz /)
      rparam(1:2) = (/ sclfac, smooth /)

      call cntc_setProfileInputFname(iwhe, c_fname, len_fname, 4, iparam, 2, rparam)

      ncase = 43
      s_ws(1:ncase)     = (/ 110.0d0, 111.0d0, 112.0d0, 113.0d0, 114.0d0, 115.0d0, 115.5d0, 116.0d0,    &
                             116.5d0, 117.0d0, 117.5d0, 118.0d0, 118.1d0, 118.2d0, 118.3d0, 118.4d0,    &
                             118.5d0, 118.6d0, 118.7d0, 118.8d0, 118.9d0, 119.0d0, 119.1d0, 119.2d0,    &
                             119.3d0, 119.4d0, 119.5d0, 119.6d0, 119.7d0, 119.8d0, 119.9d0, 120.0d0,    &
                             120.5d0, 121.0d0, 122.0d0, 123.0d0, 124.0d0, 125.0d0, 126.0d0, 127.0d0,    &
                             128.0d0, 129.0d0, 130.0d0 /)
      y_ws(1:ncase)     =     0.78d0
      z_ws(1:ncase)     =    -0.87d0
      fz_ws(1:ncase)    =    -1.0d0
      pitch_ws(1:ncase) =     0.0d0
      yaw_ws(1:ncase)   =     0.0d0
      roll_ws(1:ncase)  =     0.0d0
      vpitch(1:ncase)   =    -4.34811810d0

   elseif (trim(f_expnam).eq.'two_patches') then

      flags(1) = CNTC_ic_discns ; values(1) = 5          ! D=5: brute force method
      call cntc_setFlags(iwhe, icp, 1, flags, values)

      ! set grid discretization

      ipotcn = -1
      dx     = 0.2d0   ! [mm]
      ds     = 0.2d0  ! [mm]
      a_sep  = pi/2d0  ! [rad]
      d_sep  =  8.0d0  ! [mm]
      d_comb =  4.0d0  ! [mm]
      d_turn = 12.0d0  ! [mm]

      rvalues(1:6) = (/ dx, ds, a_sep, d_sep, d_comb, d_turn /)
      call cntc_setPotContact(iwhe, icp, ipotcn, 6, rvalues)

      ! set variable rail profile for switch/crossing

      f_fname   = '../profiles/uk_interrupt_v2.slcs'
      c_fname   =     trim(f_fname)  // C_NULL_CHAR
      len_fname = len(trim(f_fname))

      itype   = -1      ! using filename extension
      mirrory =  0      ! no mirroring
      mirrorz = -1      ! no mirroring
      sclfac  = 1d3     ! data in [m], 1000x scaling
      smooth  = 0d0     ! no smoothing
      iparam(1:4) = (/ itype, 0, mirrory, mirrorz /)
      rparam(1:2) = (/ sclfac, smooth /)

      call cntc_setProfileInputFname(iwhe, c_fname, len_fname, 4, iparam, 2, rparam)

      ncase = 1
      s_ws(1:ncase)     =    160.0d0 
      y_ws(1:ncase)     =      3.52d0
      z_ws(1:ncase)     =      0.0d0
      fz_ws(1:ncase)    =     10.0d3
      pitch_ws(1:ncase) =      0.0d0
      yaw_ws(1:ncase)   =      0.0d0
      roll_ws(1:ncase)  =      0.0d0
      vpitch(1:ncase)   =     -4.34811810d0

   elseif (trim(f_expnam).eq.'wing_brute' .or. trim(f_expnam).eq.'wing_locus') then

      if (trim(f_expnam).eq.'wing_locus') then
         flags(1) = CNTC_ic_discns ; values(1) = 2             ! D=2: contact locus method (default)
      else
         flags(1) = CNTC_ic_discns ; values(1) = 5             ! D=5: brute force method
      endif
      call cntc_setFlags(iwhe, icp, 1, flags, values)

      ! set grid discretization

      ipotcn = -1
      dx     = 0.2d0   ! [mm]
      ds     = 0.05d0  ! [mm]
      a_sep  = pi/2d0  ! [rad]
      d_sep  =  8.0d0  ! [mm]
      d_comb =  4.0d0  ! [mm]
      d_turn = 12.0d0  ! [mm]

      rvalues(1:6) = (/ dx, ds, a_sep, d_sep, d_comb, d_turn /)
      call cntc_setPotContact(iwhe, icp, ipotcn, 6, rvalues)

      ! set variable rail profile for switch/crossing

      f_fname   = '../profiles/wing_rail.slcs'
      c_fname   =     trim(f_fname)  // C_NULL_CHAR
      len_fname = len(trim(f_fname))

      itype   = -1      ! using filename extension
      mirrory =  0      ! no mirroring
      mirrorz = -1      ! no mirroring
      sclfac  = 1d0     ! already in [mm], no scaling
      smooth  = 0d0     ! no smoothing
      iparam(1:4) = (/ itype, 0, mirrory, mirrorz /)
      rparam(1:2) = (/ sclfac, smooth /)

      call cntc_setProfileInputFname(iwhe, c_fname, len_fname, 4, iparam, 2, rparam)

      ncase = 1
      s_ws(1:ncase)     =  100.0d0
      y_ws(1:ncase)     =    0.0d0
      z_ws(1:ncase)     =    0.6589d0
      fz_ws(1:ncase)    =   -1.0d0
      pitch_ws(1:ncase) =    0.0d0
      yaw_ws(1:ncase)   =    0.0d0
      roll_ws(1:ncase)  =    0.0d0
      vpitch(1:ncase)   =   -4.34811810d0

   elseif (trim(f_expnam).eq.'chalmers_flat_fz125') then

      flags(1) = CNTC_ic_discns ; values(1) = 2          ! D=2: contact locus method (default)
      flags(2) = CNTC_ic_npomax ; values(2) = 20000      ! max #elements in pot.contact
      call cntc_setFlags(iwhe, icp, 2, flags, values)

      ! set track dimensions & deviations

      ztrack = 3
      rvalues(1:11) = (/ 14d0, 0d0, 1435d0, 0.020d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 /)
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

      ! set constant rail profile

      f_fname   = '../../examples/r300_wide.prr'
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

      ncase = 26
      s_ws(1:ncase)     =       0.0d0
      y_ws(1:ncase)     =       0.0d0
      z_ws(1:ncase)     =       0.0d0
      fz_ws(1:ncase)    =     125.0d3
      pitch_ws(1:ncase) = (/ (-24.0d0 - 1.0d0 * icase, icase = 1, ncase) /) * pi/180d0
      yaw_ws(1:ncase)   =       0.0d0
      roll_ws(1:ncase)  =       0.0d0
      vpitch(1:ncase)   =      -4.08190679d0

   elseif (trim(f_expnam).eq.'rounded_flat_d09') then

      flags(1) = CNTC_ic_discns ; values(1) = 2          ! D=2: contact locus method (default)
      flags(2) = CNTC_ic_npomax ; values(2) = 20000      ! max #elements in pot.contact
      call cntc_setFlags(iwhe, icp, 2, flags, values)

      ! set track dimensions & deviations

      ztrack = 3
      rvalues(1:11) = (/ -1d0, 0d0, 0d0, 00d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0 /)
      call cntc_setTrackDimensions(iwhe, ztrack, 11, rvalues)

      ! set grid discretization

      ipotcn = -1
      dx     = 0.2d0   ! [mm]
      ds     = 0.2d0   ! [mm]
      a_sep  = pi/2d0  ! [rad]
      d_sep  =  8.0d0  ! [mm]
      d_comb =  4.0d0  ! [mm]
      d_turn = 12.0d0  ! [mm]

      rvalues(1:6) = (/ dx, ds, a_sep, d_sep, d_comb, d_turn /)
      call cntc_setPotContact(iwhe, icp, ipotcn, 6, rvalues)

      ! set constant rail profile

      f_fname   = '../profiles/circ_r300.prr'
      c_fname   =     trim(f_fname)  // C_NULL_CHAR
      len_fname = len(trim(f_fname))

      itype   = -1      ! using filename extension
      mirrory =  0      ! no mirroring
      mirrorz =  0      ! auto mirroring
      sclfac  = 1d0     ! data in [mm], no scaling
      smooth  = 0d0     ! no smoothing
      iparam(1:4) = (/ itype, 0, mirrory, mirrorz /)
      rparam(1:2) = (/ sclfac, smooth /)

      call cntc_setProfileInputFname(iwhe, c_fname, len_fname, 4, iparam, 2, rparam)

      ncase = 15
      s_ws(1:ncase)     =       0.0d0
      y_ws(1:ncase)     =       0.0d0
      z_ws(1:ncase)     =       0.0d0
      fz_ws(1:ncase)    =      10.0d3
      pitch_ws(1:ncase) = (/ ( -8.0d0 + 1.0d0 * icase, icase = 1, ncase) /) * pi/180d0
      yaw_ws(1:ncase)   =       0.0d0
      roll_ws(1:ncase)  =       0.0d0
      vpitch(1:ncase)   =      -4.44444444d0

   else

      write(*,*) 'ERROR: unknown experiment "',trim(f_expnam),'", aborting.'
      stop

   endif

   ! set wheelset dimensions

   if     (trim(f_expnam).eq.'chalmers_flat_fz125') then
      rvalues(1:3) = (/ 1360d0, -70d0, 490d0 /)
   elseif (trim(f_expnam).eq.'rounded_flat_d09') then
      rvalues(1:3) = (/    0d0,   0d0, 450d0 /)
   else
      rvalues(1:3) = (/ 1360d0, -70d0, 460d0 /)
   endif

   ewheel = 3
   call cntc_setWheelsetDimensions(iwhe, ewheel, 3, rvalues)

   ! set wheel profile

   if     (trim(f_expnam).eq.'chalmers_flat_fz125') then
      f_fname   = '../../examples/S1002_flat.slcw'
      smooth    = 5d0     ! lambda smoothing
   elseif (trim(f_expnam).eq.'rounded_flat_d09') then
      f_fname   = '../profiles/flat_d09.slcw'
      smooth    = 0d0     ! no smoothing
   else
      f_fname   = '../profiles/MBench_S1002_v3.prw'
      smooth    = 0d0     ! no smoothing
   endif
   c_fname   =     trim(f_fname)  // C_NULL_CHAR
   len_fname = len(trim(f_fname))

   itype   = -1      ! using filename extension
   mirrory =  0      ! no y-mirroring
   mirrorz = -1      ! no z-mirroring
   sclfac  = 1d0     ! already in [mm], no scaling
   iparam(1:4) = (/ itype, 0, mirrory, mirrorz /)
   rparam(1:2) = (/ sclfac, smooth /)

   call cntc_setProfileInputFname(iwhe, c_fname, len_fname, 4, iparam, 2, rparam)

   !---------------------------------------------------------------------------------------------------------
   ! Run cases, with outputs in out-file
   !---------------------------------------------------------------------------------------------------------

   do irep = 1, nrep
      if (nrep.gt.1) write(*,*) ' starting repetition irep=',irep,' of',nrep
   
      do icase = 1, ncase

         ! adjust discretisation settings for cw_interrupt

         if (icase.eq.2 .and. trim(f_expnam).eq.'cw_interrupt') then
            flags(1) = CNTC_ic_discns ; values(1) = 5          ! D=5: brute force method
            call cntc_setFlags(iwhe, icp, 1, flags, values)
         endif

         ! set wheelset position & velocity

         rvalues(1:6) = (/ s_ws(icase), y_ws(icase), z_ws(icase), roll_ws(icase), yaw_ws(icase),        &
                                                                                         pitch_ws(icase) /)
         ewheel = 2
         call cntc_setWheelsetPosition(iwhe, ewheel, 6, rvalues)

         rvalues(1:6) = (/ 2000d0, 0d0, 0d0, 0d0, 0d0, vpitch(icase) /)
         ewheel = 2
         call cntc_setWheelsetVelocity(iwhe, ewheel, 6, rvalues)

         if (fz_ws(icase).gt.0d0) then
            call cntc_setVerticalForce(iwhe, fz_ws(icase))

            flags(1) = CNTC_ic_norm; values(1) = 1             ! N=1: vert.force prescribed
            call cntc_setFlags(iwhe, icp, 1, flags, values)
         else
            flags(1) = CNTC_ic_norm; values(1) = 0             ! N=0: approach prescribed
            call cntc_setFlags(iwhe, icp, 1, flags, values)
         endif

         if (ldebug.ge.1)  &
            write(*,'(/,a,i3)') '  ...start processing case', icase

         ! perform the actual calculation for this wheel

         call cntc_calculate(iwhe, icp, ierr)

         ! check/report error conditions

         if (ierr.eq.CNTC_err_allow) then
            write(*,'(a,i5,a)') '      ERROR: no valid license found'
            stop
         elseif (ierr.eq.CNTC_err_profil) then
            write(*,'(a,i5,a)') '      ERROR: the rail and/or wheel profile files could not be ' //        &
                                'found or processed'
            stop
         elseif (ierr.lt.0) then
            write(*,'(a,i5,a)') '      ERROR: an error occured in the CONTACT library, ierr=',ierr,'.'
         elseif (ierr.gt.0) then
            write(*,'(a,i5,a)') '      ERROR: contact extends outside bounding box,',ierr,                 &
                                ' elements at boundary'
         end if

      end do ! icase
   end do ! irep

   ! finalize, close CONTACT library

   call cntc_finalize(iwhe)

end program test_varprof

!============================================================================================================

