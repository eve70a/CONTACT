!============================================================================================================
! This program tests parallel computation using the CONTACT library.
!
! Copyright 2008-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!============================================================================================================
program test_table
   use, intrinsic        :: iso_c_binding, only: C_CHAR, C_NULL_CHAR
#ifdef _OPENMP
   use omp_lib
#endif
   implicit none
   include 'caddon_flags.inc'
!--local variables:
   real(kind=8), parameter :: pi     = 4d0*atan(1d0)
   character(len=30,kind=C_CHAR) :: c_path, c_expnam
   character(len=256)      :: tmpstr, outdir, expnam
   integer, parameter      :: mxflgs = CNTC_FLAGS_DIM
   integer                 :: flags(mxflgs), values(mxflgs), iparam(mxflgs)
   real(kind=8)            :: rparam(mxflgs)
   integer                 :: ire, icp, imodul, ver, ierr, iout, l_path, l_expnam, icount, ithrd, nthrd
   integer                 :: mx, my, ipotcn, gdigit, ldigit, mdigit, maxgs, maxnr, maxin, maxout,      &
                              inislp, kdowfb
   real(kind=8)            :: g1, nu1, scale, fstat, fkin, veloc, eps, omegah, omegas, omgslp,          &
                              fdecay, betath, d_ifc, d_lin, d_cns, d_slp, pow_s, a1, aob, dx, dy
   real(kind=8)            :: cmgn, cang, cksi, ceta, cphi, fn, fx, fy, mz, carea, harea, sarea
   integer                 :: iell, ispin, ivang, ivmgn, irem, icase
   integer, parameter      :: num_ellip   = 5
   integer, parameter      :: num_spin    = 4
   integer, parameter      :: num_vecangl = 7
   integer, parameter      :: num_vecmagn = 23
   integer, parameter      :: ncase = num_ellip * num_spin * num_vecangl * num_vecmagn
   real(kind=8), parameter :: ellip(num_ellip)   = (/ 0.2      , 0.5      , 1.0      , 2.0      , 5.0       /)
   real(kind=8), parameter :: cp   (num_ellip)   = (/ 0.9063561, 0.7580279, 0.6463304, 0.5351031, 0.3978158 /)
   real(kind=8), parameter :: rho  (num_ellip)   = (/ 1.844071 , 1.479543 , 1.0      , 0.5204566, 0.1559291 /)
   real(kind=8), parameter :: spin   (num_spin)    = (/ 0.0, 0.5, 1.0, 2.0 /)
   real(kind=8), parameter :: vecangl(num_vecangl) = (/ -pi/2, -pi/3, -pi/6, 0d0, pi/6, pi/3, pi/2 /)
   real(kind=8), parameter :: vecmagn(num_vecmagn) = (/ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 1.1, 1.4, &
                                                        1.7, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, &
                                                        9.0, 10., 12.  /)
   integer                 :: ncon_list(ncase), nslp_list(ncase)
   real(kind=8)            :: fx_list(ncase), fy_list(ncase), mz_list(ncase)
!--external functions used:
#include "contact_addon.ifc"

   ! Handle optional command line arguments
   ! Usage: test_table [mx] [expnam] [gdigit] [fdecay]

   mx     = 11
   icount = command_argument_count()
   if (icount.ge.1) then
      call get_command_argument(1, tmpstr, status=ierr)
      if (ierr.ge.0) read(tmpstr,'(i)') mx
   endif
   if (icount.ge.2) then
      call get_command_argument(2, expnam, status=ierr)
   else
      write(tmpstr,'(i3.3)') mx
      expnam = 'test_table_mx' // trim(tmpstr)
   endif
   if (icount.ge.3) then
      call get_command_argument(3, tmpstr, status=ierr)
      if (ierr.ge.0) read(tmpstr,'(i)') gdigit
   else
      gdigit = 0
   endif
   if (icount.ge.4) then
      call get_command_argument(4, tmpstr, status=ierr)
      if (ierr.ge.0) read(tmpstr,'(f)') fdecay
   else
      fdecay = 1d0
   endif
   my  = mx

   ! Initialize once using cntc_initializeFirst

   outdir   = ' '
   c_path   = trim(outdir) // C_NULL_CHAR
   l_path   = len(trim(outdir))
   c_expnam = trim(expnam) // C_NULL_CHAR
   l_expnam = len(trim(expnam))
   iout = 0
   call cntc_initializeFirst_new(ver, ierr, iout, c_path, c_path, c_expnam, l_path, l_path, l_expnam)

   ! Set global flags: debug output, use of OpenMP

   flags(1) = CNTC_if_idebug ; values(1) =  0
   flags(2) = CNTC_if_openmp ; values(2) = -1
   call cntc_setGlobalFlags(2, flags, values)

   ! Initialize one result element for each thread

   nthrd = 1
#ifdef _OPENMP
   nthrd = omp_get_max_threads()
#endif
   imodul = 3
   do ire = 1, nthrd
      call cntc_initialize(ire, imodul, ver, ierr, c_path, 1)
   enddo
   icp  = 1

   ! Initialize the contact problem for each thread

   do ire = 1, nthrd
      flags = 0 ; values = 0
      flags(1) = CNTC_if_units  ; values(1) = CNTC_un_cntc
      flags(2) = CNTC_if_wrtinp ; values(2) = 0
      flags(3) = CNTC_ic_matfil ; values(3) = 0
      flags(4) = CNTC_ic_output ; values(4) = 0
      flags(5) = CNTC_ic_flow   ; values(5) = 2
      flags(6) = CNTC_ic_tang   ; values(6) = 3

      call cntc_setFlags(ire, icp, mxflgs, flags, values)

      ldigit = 0
      fstat  = 0.3d0
      fkin   = 0.3d0
      rparam(1:2) = (/ fstat, fkin /)
      call cntc_setFrictionMethod(ire, icp, ldigit, 2, rparam)

      mdigit = 0
      g1     = 1d0
      nu1    = 0.28d0
      rparam(1:4) = (/ nu1, nu1, g1, g1 /)
      call cntc_setMaterialParameters(ire, icp, mdigit, 4, rparam)

      veloc = 30000d0
      call cntc_setReferenceVelocity(ire, icp, veloc)

      fn  = 1d0
      call cntc_setNormalForce(ire, icp, fn)

      maxgs   = 299
      maxnr   =  30
      maxin   =   1
      maxout  =   1
      eps     =   1d-6
      omegah  =  0.5d0
      omegas  =  0.5d0
      omgslp  =  1.0d0
      inislp  =   1
      ! fdecay given on command line: input argument
      kdowfb  =   1
      betath  =  0.05d0
      d_ifc   =  2.0d0
      d_lin   = -1.0d0
      d_cns   =  1.0d0
      d_slp   =  2.6d0
      pow_s   =  1.0d0

      if (gdigit.eq.2) then
         iparam(1:5) = (/ maxgs, maxin, maxnr, maxout, inislp /)
         rparam(1:4) = (/ eps, omegah, omegas, omgslp /)
         call cntc_setSolverFlags(ire, icp, gdigit, 5, iparam, 4, rparam)
      elseif (gdigit.eq.5) then
         iparam(1:5) = (/ maxgs, maxin, maxnr, maxout, kdowfb /)
         rparam(1:8) = (/ eps, fdecay, betath, d_ifc, d_lin, d_cns, d_slp, pow_s /)
         call cntc_setSolverFlags(ire, icp, gdigit, 5, iparam, 8, rparam)
      else ! use gdigit = 0
         iparam(1:4) = (/ maxgs, maxin, maxnr, maxout /)
         rparam(1:1) = (/ eps /)
         call cntc_setSolverFlags(ire, icp, gdigit, 4, iparam, 1, rparam)
      endif

   end do

   ! Calculate and print results for both Result Elements

!$omp parallel &
!$omp    default(none)  &
!$omp    shared(icp, mx, my, fx_list, fy_list, mz_list, ncon_list, nslp_list, fstat)                    &
!$omp    private(ithrd, ire, ivmgn, ivang, ispin, iell, irem, ierr, ipotcn, a1, aob, scale, dx, dy,     &
!$omp            rparam,    cmgn, cang, cksi, ceta, cphi, fn, fx, fy, mz, carea, harea, sarea)
!$omp do
   do icase = 1, ncase

      ! set the result-element number - one ire per thread

      ithrd = 0
#ifdef _OPENMP
      ithrd = omp_get_thread_num()
#endif
      ire   = ithrd + 1

      ! decompose the case-number to its constituents

      ! icase =
      !     ivmgn    +
      !    (ivang-1) * num_vecmagn +
      !    (ispin-1) * num_vecmagn*num_vecangl +
      !    (iell -1) * num_vecmagn*num_vecangl*num_spin

      ivmgn = mod(icase-1, num_vecmagn) + 1
      irem  = (icase - ivmgn   ) / num_vecmagn
      ivang = mod(irem   , num_vecangl) + 1
      irem  = (irem - (ivang-1)) / num_vecangl
      ispin = mod(irem   , num_spin)    + 1
      irem  = (irem - (ispin-1)) / num_spin
      iell  = mod(irem   , num_ellip  ) + 1

      ! set the ellipticity for the case

      ipotcn = -2       ! ellipticity aob
      a1     = 1d0
      aob    = ellip(iell)
      scale  = real(mx) / real(max(1, mx-4))

      rparam(1) = real(mx)
      rparam(2) = real(my)
      rparam(3) = a1
      rparam(4) = aob
      rparam(5) = scale
      call cntc_setHertzContact(ire, icp, ipotcn, 5, rparam)

      ! solve the Hertzian problem: obtain Cp, Rho, Dx, Dy
      ! call hzsol(gd%ic, gd%geom, gd%potcon, gd%hertz, gd%kin, gd%mater)

      ! set the creepages - undo scaling of non-dimensionalised values

      cmgn = vecmagn(ivmgn)
      cang = vecangl(ivang)

      cksi = cmgn * fstat * cos(cang) * cp(iell) / rho(iell)
      ceta = cmgn * fstat * sin(cang) * cp(iell) / rho(iell)
      cphi =        fstat * spin(ispin) / rho(iell)

      call cntc_setCreepages(ire, icp, cksi, ceta, cphi)

      if (.false.) then
         write(*,'(4(a,i2),4(a,f10.4))') 'case',iell,',',ispin,',',ivang,',',ivmgn,': aob=',ellip(iell), &
                ', spin=',spin(ispin), ', angl=',vecangl(ivang),', mgn=',vecmagn(ivmgn)
      endif

#ifndef _OPENMP
      if (icase.eq.1) then
         flags = 0 ; values = 0
         flags(1) = CNTC_ic_output ; values(1) = 2
      else
         flags(1) = CNTC_ic_output ; values(1) = 0
      endif
      call cntc_setFlags(ire, icp, 1, flags, values)
#endif

      ! solve the contact problem

      call cntc_calculate(ire, icp, ierr)

      if (ierr.eq.CNTC_err_allow) then
         write(*,'(a,i5,a)') '      ERROR: no valid license found'
         stop
      elseif (ierr.lt.0) then
         write(*,'(a,i5,a)') '      ERROR: code ierr=',ierr,'.'
      elseif (ierr.gt.0) then
         write(*,'(a,i5,a)') '      ERROR: contact extends outside bounding box,',ierr,' elements at boundary'
      end if

      ! store the total forces for the contact patch

      call cntc_getContactForces(ire, icp, fn, fx, fy, mz)

      fx_list(icase) = fx / (fstat * fn)
      fy_list(icase) = fy / (fstat * fn)
      mz_list(icase) = mz / (fstat * fn * cp(iell))

      call cntc_getGridDiscretization(ire, icp, dx, dy)
      call cntc_getContactPatchAreas(ire, icp, carea, harea, sarea)
      ncon_list(icase) = nint(carea / (dx*dy))
      nslp_list(icase) = nint(sarea / (dx*dy))

   enddo
!$omp end do
!$omp end parallel

   do icase = 1, ncase
      write(*,'(3(a,f11.6),a,f6.1)') 'Fx=', fx_list(icase), ', Fy=', fy_list(icase),          &
            ', Mz=', mz_list(icase), ', %Slip=', real(100*nslp_list(icase))/real(ncon_list(icase))
   enddo

   ! Finalize all result elements used

   do ire = 1, nthrd
      call cntc_finalize(ire)
   enddo

end program test_table

!============================================================================================================
