!============================================================================================================
! This program computes USETAB tables using the CONTACT library.
!
! Copyright 2008-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!============================================================================================================
program usetab_table
   use, intrinsic        :: iso_c_binding, only: C_CHAR, C_NULL_CHAR
#ifdef _OPENMP
   use omp_lib
#endif
   implicit none
   include 'caddon_flags.inc'
!--local variables:
   integer,      parameter :: itab = 7
   logical,      parameter :: filter_small = .true.
   real(kind=8), parameter :: small_value = 5d-7
   integer,      parameter :: num_aob=itab+18, num_xi=2*(itab+1), num_eta=num_xi, num_psi=2*num_xi
   integer,      parameter :: ncase = num_aob*num_xi*num_eta*num_psi
   integer,      parameter :: mxflgs = CNTC_FLAGS_DIM
   real(kind=8), parameter :: pi = 4d0*atan(1d0)
   character(len=30,kind=C_CHAR) :: c_path, c_expnam
   character(len=256)      :: tmpstr, expnam, outfil
   integer                 :: flags(mxflgs), values(mxflgs), iparam(mxflgs)
   real(kind=8)            :: rparam(mxflgs)
   integer                 :: ire, icp, imodul, ver, ierr, iout, lout, l_expnam, icount, ithrd, nthrd
   integer                 :: mx, my, imeth, maxgs, maxnr, maxin, maxout
   real(kind=8)            :: g1, nu1, fstat, fkin, veloc, eps
   real(kind=8)            :: aob, aa, bb, cc, cx, cy, cz, scale
   real(kind=8)            :: xi, eta, psi, cksi, ceta, cphi, fn, fx, fy, mz, tmp
   integer                 :: i, iaob, ixi, ieta, ipsi, irem, icase, ipotcn
   real(kind=8)            :: aob_list(num_aob), xi_list(num_xi), eta_list(num_eta), psi_list(num_psi)
   real(kind=8)            :: fx_list(ncase), fy_list(ncase), mz_list(ncase)
   character*9             :: str_aob, str_xi, str_eta
!--external functions used:
#include "contact_addon.ifc"

   ! Handle optional command line arguments
   ! Usage: usetab_table [mx] [expnam] [tabcon]

   mx     = 12
   icount = command_argument_count()
   if (icount.ge.1) then
      call get_command_argument(1, tmpstr, status=ierr)
      if (ierr.ge.0) read(tmpstr,'(i)') mx
   endif
   if (icount.ge.2) then
      call get_command_argument(2, expnam, status=ierr)
   else
      write(tmpstr,'(i3.3)') mx
      expnam = 'contact_addon' // trim(tmpstr)
   endif
   if (icount.ge.3) then
      call get_command_argument(3, outfil, status=ierr)
   else
      if (mx.le.99) then
         write(tmpstr,'(i2.2)') mx
      else
         write(tmpstr,'(i3.3)') mx
      endif
      outfil = 'tabcon_mx' // trim(tmpstr)
   endif
   my  = mx

   ! Initialize once using cntc_initializeFirst

   c_path   = ' ' // C_NULL_CHAR
   c_expnam = trim(expnam) // C_NULL_CHAR
   l_expnam = len(trim(expnam))
   iout = 0
   call cntc_initializeFirst_new(ver, ierr, iout, c_path, c_path, c_expnam, 1, 1, l_expnam)

   ! Set global flags: debug output, use of OpenMP

   flags(1) = CNTC_if_idebug ; values(1) =  0
   flags(2) = CNTC_if_openmp ; values(2) = -1
   call cntc_setGlobalFlags(2, flags, values)

   ! Initialize one result element for each thread

   nthrd = 1
#ifdef _OPENMP
   nthrd = omp_get_max_threads()
#endif
   imodul   = 3
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
      flags(5) = CNTC_ic_flow   ; values(5) = 0
      flags(6) = CNTC_ic_tang   ; values(6) = 3

      call cntc_setFlags(ire, icp, mxflgs, flags, values)

      imeth  = 0
      fstat  = 1.0d0
      fkin   = 1.0d0
      rparam(1:2) = (/ fstat, fkin /)
      call cntc_setFrictionMethod(ire, icp, imeth, 2, rparam)

      g1  = 1d0
      nu1 = 0.28d0
      call cntc_setMaterialProperties(ire, icp, g1, nu1, g1, nu1)

      veloc = 1d0
      call cntc_setReferenceVelocity(ire, icp, veloc)

      fn  = 1d0
      call cntc_setNormalForce(ire, icp, fn)

      imeth  =   0
      maxgs  = 999
      maxnr  =  30
      maxin  =   5
      maxout =   1
      eps    = 1d-6
      iparam(1:4) = (/ maxgs, maxin, maxnr, maxout /)
      rparam(1:1) = (/ eps /)
      call cntc_setSolverFlags(ire, icp, imeth, 4, iparam, 1, rparam)

   end do

   ! Set the ellipticities and scaled creepages to use

   aob_list = (/ (real(i)/real(itab,kind=8), i=1,itab), 200d0, 130d0, 90d0, 65d0, 45d0, 32d0, 22.4d0, 16d0,  &
                  11.2d0, 8d0, 5.6d0, 4d0, 2.8d0, 2d0, 1.7d0, 1.4d0, 1.2d0, 1d0 /)
   xi_list  = (/ (real(i)/real(itab,kind=8), i=0,itab), (real(itab,kind=8)/max(1e-5,real(i)), i=0,itab) /)
   eta_list = (/ (real(i)/real(itab,kind=8), i=0,itab), (real(itab,kind=8)/max(1e-5,real(i)), i=0,itab) /)
   psi_list = (/ ( real(i)/real(itab,kind=8), i=0,itab), ( real(itab,kind=8)/max(1e-5,real(i)), i=0,itab),   &
                 (-real(i)/real(itab,kind=8), i=0,itab), (-real(itab,kind=8)/max(1e-5,real(i)), i=0,itab) /)

   write(*,17) 'aob_list =', (aob_list(i), i=1,num_aob)
   write(*,18) 'xi_list  =', (xi_list(i) , i=1,num_xi)
   write(*,18) 'eta_list =', (eta_list(i), i=1,num_eta)
   write(*,18) 'psi_list =', (psi_list(i), i=1,num_psi)
17 format(a,/,4(4x,7g12.4,:,/))
18 format(a,/,4(4x,8g12.4,:,/))

   ! Use parallel loop to calculate forces for all cases of the table

!$omp parallel &
!$omp    default(none)  &
!$omp    shared(icp, aob_list, xi_list, eta_list, psi_list, fx_list, fy_list, mz_list,                  &
!$omp           mx, my, fstat, fn, nu1, g1) &
!$omp    private(ithrd, ire, ixi, ieta, ipsi, iaob, irem, ierr, ipotcn, aob, aa, bb, cc, scale, rparam, &
!$omp            cx, cy, cz, xi, eta, psi, cksi, ceta, cphi, tmp, fx, fy, mz)
!$omp do schedule(dynamic,16)
   do icase = 1, ncase

      ! set the result-element number - one ire per thread

      ithrd = 0
#ifdef _OPENMP
      ithrd = omp_get_thread_num()
#endif
      ire   = ithrd + 1

      ! decompose the case-number to its constituents

      ! icase =
      !     ipsi     +
      !    (ieta -1) * num_psi +
      !    (ixi  -1) * num_psi * num_eta +
      !    (iaob -1) * num_psi * num_eta * num_xi

      ipsi  = mod(icase-1, num_psi) + 1
      irem  = (icase - ipsi    ) / num_psi
      ieta  = mod(irem   , num_eta) + 1
      irem  = (irem - (ieta -1)) / num_eta
      ixi   = mod(irem   , num_xi)  + 1
      irem  = (irem - (ixi  -1)) / num_xi
      iaob  = mod(irem   , num_aob) + 1
      ! write(*,*) 'icase=',icase,': iaob=',iaob,', ixi=',ixi,', ieta=',ieta,', ipsi=',ipsi

      ! Hertzian discretization with 1 free element on all sides

      scale = real(mx) / real(mx-2)

      ! semi-axes derived from ellipticity aob

      aob  = aob_list(iaob)
      cc   = 1d0
      aa   = sqrt(cc**2*aob)
      bb   = sqrt(cc**2/aob)

      ipotcn = -3  ! semi-axes
      rparam(1) = real(mx)
      rparam(2) = real(my)
      rparam(3) = aa
      rparam(4) = bb
      rparam(5) = scale
      call cntc_setHertzContact(ire, icp, ipotcn, 5, rparam)

      ! write(*,*) 'icase=',icase,': aob=',aob,', aa=',aa,', bb=',bb

      ! get Kalker coefficients, un-scale the creepages

      xi   = xi_list(ixi)
      eta  = eta_list(ieta)
      psi  = psi_list(ipsi)
      ! write(*,'(a,i6,4(a,g12.4))') 'icase=',icase,': aob=',aob,', xi=',xi,', eta=',eta,', psi=',psi

      call linrol(nu1, aob, cx, cy, cz)
      cksi = -xi  * 3d0 * fstat * fn / (g1 * cc**2 * cx)
      ceta = -eta * 3d0 * fstat * fn / (g1 * cc**2 * cy)
      cphi = -psi       * fstat * fn / (g1 * cc**3 * cz)
      call cntc_setCreepages(ire, icp, cksi, ceta, cphi)

      ! solve the contact problem

      call cntc_calculate(ire, icp, ierr)

      if (ierr.eq.CNTC_err_allow) then
         write(*,'(a,i5,a)') '      ERROR: no valid license found'
         stop
      elseif (ierr.lt.0) then
         write(*,'(a,i5,a)') '      ERROR: code ierr=',ierr,'.'
         write(*,*) 'aob=',aob,', cksi,ceta,cphi=',cksi,ceta,cphi
      elseif (ierr.gt.0) then
         write(*,'(a,i5,a)') '      ERROR: contact extends outside bounding box,',ierr,' elements at boundary'
      end if

      ! store the total forces for the contact patch

      call cntc_getContactForces(ire, icp, tmp, fx, fy, mz)

      if (filter_small) then
         if (abs(fx).lt.small_value) fx = 0d0
         if (abs(fy).lt.small_value) fy = 0d0
         if (abs(mz).lt.small_value) mz = 0d0
      endif

      fx_list(icase) = fx
      fy_list(icase) = fy
      mz_list(icase) = mz

   enddo
!$omp end do
!$omp end parallel

   ! Finalize all result elements used

   do ire = 1, nthrd
      call cntc_finalize(ire)
   enddo

   ! Write table in USETAB format

   lout = 27
   open(unit=lout, file=outfil)
   write(lout,'(a)') '! "tabcon", table of the Hertzian creep-force law (contact vXX.X, rYYYY)'
   write(lout,'(a)') '! Authors: J.J. Kalker & E.A.H. Vollebregt'
   write(lout,'(a)') '! Delft University of Technology / Vtech CMCC'
   write(lout,'(a)') ' '
   write(lout,'(5i6,a)') itab, num_aob-itab, 10, 0, 0, '  itab, naob1, imap'
   write(lout,'(20( 5(f12.6,:,'',''), /))') (aob_list(i), i=itab+1,num_aob)
   write(lout,'(a,i3,a)') '! Fx,Fy,aM computed for aob,xi,eta,psi with G=1,nu=0.28,N=1,c=Z.Z,mxeff=',mx-2,'.'
   write(lout,'(a)') ' '

   do icase = 1, ncase
      if (mod(icase, num_psi).eq.1) then
         ipsi  = mod(icase-1, num_psi) + 1
         irem  = (icase - ipsi    ) / num_psi
         ieta  = mod(irem   , num_eta) + 1
         irem  = (irem - (ieta -1)) / num_eta
         ixi   = mod(irem   , num_xi)  + 1
         irem  = (irem - (ixi  -1)) / num_xi
         iaob  = mod(irem   , num_aob) + 1
         aob  = aob_list(iaob)
         xi   = xi_list(ixi)
         eta  = eta_list(ieta)
         psi  = psi_list(ipsi)

         if (aob.lt.100d0) then
            write(str_aob,'(f9.6)') aob
         else
            write(str_aob,'(f9.3)') aob
         endif
         if (xi.lt.100d0) then
            write(str_xi,'(f9.6)') xi
         else
            write(str_xi,'(f9.1)') xi
         endif
         if (eta.lt.100d0) then
            write(str_eta,'(f9.6)') eta
         else
            write(str_eta,'(f9.1)') eta
         endif
         write(lout,'(6a)') '! aob=',str_aob,', xi=',str_xi,', eta=',str_eta
      endif
      ! write(*,'(3(a,f11.6))') 'Fx=', fx_list(icase), ', Fy=', fy_list(icase), ', Mz=', mz_list(icase)
      write(lout,'(1x,3(1x,f11.6))') fx_list(icase), fy_list(icase), mz_list(icase)
   enddo
   write(lout,*) 'End of Table'
   close(lout)

end program usetab_table

!============================================================================================================

subroutine linrol(nu, aob, cx, cy, cz)
!--purpose: determine the Hertzian creepage coefficients cx,cy,cz from Table E.3 of [4].
   implicit none
!--subroutine arguments:
   real(kind=8), intent(in)  :: nu, aob
   real(kind=8), intent(out) :: cx, cy, cz
!--local variables:
   integer      :: icase, irow, inu, icoef
   real(kind=8) :: pi, al, g, fac, c(3,3), cc(3)

   ! table E.3 of Kalker coefficients Cij of the linear theory of rolling
   !              contact for elliptical contact areas:

   real(kind=8) :: cij(2,10,3,3)
   data ((((cij(icase,irow,inu,icoef), inu=1,3), icoef=1,3), irow=1,10), icase=1,2)             &
    /  2.51, 3.31, 4.85,   2.51, 2.52, 2.53,   .334, .473, .731,                                &
       2.59, 3.37, 4.81,   2.59, 2.63, 2.66,   .483, .603, .809,                                &
       2.68, 3.44, 4.80,   2.68, 2.75, 2.81,   .607, .715, .889,                                &
       2.78, 3.53, 4.82,   2.78, 2.88, 2.98,   .720, .823, .977,                                &
       2.88, 3.62, 4.83,   2.88, 3.01, 3.14,   .827, .929, 1.07,                                &
       2.98, 3.72, 4.91,   2.98, 3.14, 3.31,   .930, 1.03, 1.18,                                &
       3.09, 3.81, 4.97,   3.09, 3.28, 3.48,   1.03, 1.14, 1.29,                                &
       3.19, 3.91, 5.05,   3.19, 3.41, 3.65,   1.13, 1.15, 1.40,                                &
       3.29, 4.01, 5.12,   3.29, 3.54, 3.82,   1.23, 1.36, 1.51,                                &
       3.40, 4.12, 5.20,   3.40, 3.67, 3.98,   1.33, 1.47, 1.63,                                &
           10.7, 11.7, 12.9,   10.7, 12.8, 16.0,   12.2, 14.6, 18.0,                            &
       6.96, 7.78, 8.82,   6.96, 8.14, 9.79,   5.72, 6.63, 7.89,                                &
       5.57, 6.34, 7.34,   5.57, 6.40, 7.51,   3.79, 4.32, 5.01,                                &
       4.84, 5.57, 6.57,   4.84, 5.48, 6.31,   2.88, 3.24, 3.70,                                &
       4.37, 5.10, 6.11,   4.37, 4.90, 5.56,   2.35, 2.62, 2.96,                                &
       4.06, 4.78, 5.80,   4.06, 4.50, 5.04,   2.01, 2.23, 2.50,                                &
       3.82, 4.54, 5.58,   3.82, 4.21, 4.67,   1.76, 1.95, 2.18,                                &
       3.65, 4.36, 5.42,   3.65, 3.99, 4.39,   1.58, 1.75, 1.94,                                &
       3.51, 4.22, 5.30,   3.51, 3.81, 4.16,   1.44, 1.59, 1.77,                                &
       3.40, 4.12, 5.20,   3.40, 3.67, 3.98,   1.33, 1.47, 1.63  /                              
                                                                                                
   pi = 4d0*datan(1d0)

   ! determine g = min(a/b, b/a), remember which case it is

   if (aob.le.1d0) then
      g = aob
      icase = 1
   else
      g = 1d0 / aob
      icase = 2
   end if

   ! determine cx, cy, cz using asymptotes or interpolation

   if (icase.eq.1 .and. g.lt.0.101) then

      ! first asymptotes, g \downarrow 0, a<b

      cx = pi**2 / 4d0 / (1d0-nu)
      cy = pi**2 / 4d0
      cz = pi*sqrt(g) / 3d0 / (1d0-nu) * (1d0 + nu * (dlog(16d0 / g) - 5d0))

   else if (icase.eq.2 .and. g.lt.0.101) then

      ! last asymptotes, g \downarrow 0, a>b

      al = dlog(16d0 / g**2)
      cx = 2d0*pi / (al-2d0*nu) / g * (1d0 + (3d0 - log(4d0)) / (al - 2d0*nu))
      cy = 2d0*pi / g * (1d0 + (1d0-nu) * (3d0-log(4d0)) /                                      &
                         (      (1d0-nu) * al + 2d0*nu)  ) /                                    &
                                            ((1d0-nu) * al + 2d0*nu) 
      cz = 2d0*pi / 3d0 / g**1.5d0 / ((1d0-nu) * al - 2d0 + 4d0*nu)

   else

      ! in range of the table, 0.101 < g <= 1.0

      ! take values from table "icase" (1: a<b, 2: a>b)
      ! take irow == floor(10*g-eps), e.g. g = 0.12 --> irow = 1
      ! interpolation weight == rem(10*g) --> fac = 0.2

      irow  = 10d0*g - 0.005
      fac = 10d0*g - real(irow)

      ! perform interpolation for all c_i ("icoef") and all nu ("inu")

      do inu = 1, 3
         do icoef = 1, 3
            c(inu,icoef) =    fac  * cij(icase, irow+1, inu, icoef) +                           &
                       (1d0 - fac) * cij(icase, irow  , inu, icoef)
         enddo
      enddo

      ! interpolate to actual nu
      ! formula?  fitting parabola?  inverse interpolation?

      do icoef = 1,3
         cc(icoef) =   (nu-0.25d0) * (nu-0.50d0) *  8d0 / c(1,icoef)                            &
                     - (nu-0.50d0) *  nu         * 16d0 / c(2,icoef)                            &
                     +  nu         * (nu-0.25d0) *  8d0 / c(3,icoef)
      enddo

      cx = 1d0 / cc(1)
      cy = 1d0 / cc(2)
      cz = 1d0 / cc(3)
   end if

end subroutine linrol

!============================================================================================================
