!------------------------------------------------------------------------------------------------------------
! m_arcfit - experimental profile smoothing using circle fitting approach
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------

module m_arcfit

use m_hierarch_data

implicit none
private

   public  make_arcfit
   public  fit_circle

#ifdef _WIN32
   !dir$ attributes c, reference  :: dgecon, dgetrf, dgetri ! on IA-32: calling convention for LAPACK routines
#endif

contains

!------------------------------------------------------------------------------------------------------------

   subroutine make_arcfit(prf_grid, prf_fname, ds_max)
!--purpose: temporary routine for testing of the arc-fit approach
!           approximate input-profile by fitting of circular arcs
      implicit none
!--subroutine parameters:
      type(t_grid)       :: prf_grid
      character(len=*)   :: prf_fname
      real(kind=8)       :: ds_max
!--local variables:
      integer            :: idebug, ipdbg, nprf, iprf, jprf, ista, iend, nmeas
      real(kind=8)       :: s_offs, s_ref
      real(kind=8), dimension(:), allocatable :: ds, ymeas, zmeas, wgt, ymid, zmid, phi, dst, rho

   write(bufout,'(3a)') ' make_arcfit: starting for "', trim(prf_fname), '"'
   call write_log(1, bufout)
   write(bufout,'(a,i6,a,f6.2)') ' make_arcfit: profile has', prf_grid%ntot,' points, ds_max=', ds_max
   call write_log(1, bufout)

   idebug = 1
   ipdbg  = -143
   nprf   = prf_grid.ntot

   allocate(ds(nprf), ymeas(nprf), zmeas(nprf), wgt(nprf), ymid(nprf), zmid(nprf), phi(nprf),           &
            dst(nprf), rho(nprf))

   ! calculate step size ds centered at profile points

   ds(1) = prf_grid%s_prf(2) - prf_grid%s_prf(1)
   do iprf = 2, nprf-1
      ds(iprf) = 0.5d0 * (prf_grid%s_prf(iprf+1) - prf_grid%s_prf(iprf-1))
   enddo
   ds(nprf) = prf_grid%s_prf(nprf) - prf_grid%s_prf(nprf-1)

   ! loop over profile points, fit circular arcs to sliding window of [-ds_max, ds_max]

   ista = 1
   iend = 1

   do iprf = 1, nprf

      ! select measured values to be fitted, assuming s in increasing order

      do while(ista.lt.iprf .and. prf_grid%s_prf(iprf)-prf_grid%s_prf(ista).gt.ds_max)
         ista = ista + 1
      enddo
      do while(iend.lt.nprf .and. prf_grid%s_prf(iend+1)-prf_grid%s_prf(iprf).lt.ds_max)
         iend = iend + 1
      enddo
      nmeas = iend - ista + 1

      if (nmeas.le.2) then
         if (idebug.ge.1) then
            write(bufout,'(a,i4,a,i2,2a,f5.2)') ' fit_profile: iprf =',iprf,': ',nmeas,' points with',  &
                ' dst <=', ds_max
            call write_log(1, bufout)
         endif

         ista  = max(1, iprf-2)
         iend  = min(nprf, iprf+2)
         nmeas = iend - ista + 1
      endif

      if (idebug.ge.3 .or. iprf.eq.ipdbg) then
         write(bufout,'(4(a,i4))') ' fit_profile: iprf =',iprf,': selected',nmeas,' points',ista,' :',iend
         call write_log(1, bufout)
         s_ref = prf_grid%s_prf(iprf)
         write(bufout,'(a,f8.3,2(a,2f8.3))') ' s_ref=',s_ref,', s=',prf_grid%s_prf(ista)-s_ref,        &
                prf_grid%s_prf(ista+1)-s_ref, '...', prf_grid%s_prf(iend-1)-s_ref,                      &
                prf_grid%s_prf(iend)-s_ref
         call write_log(1, bufout)
      endif

      do jprf = ista, iend
         ymeas(jprf)  = prf_grid%y(jprf) - prf_grid%y(iprf)
         zmeas(jprf)  = prf_grid%z(jprf) - prf_grid%z(iprf)
         s_offs = prf_grid%s_prf(jprf) - prf_grid%s_prf(iprf)
         if (abs(s_offs).gt.1d-20) then
            wgt(jprf) = ds(jprf) * sin( pi * s_offs / ds_max ) / (pi * s_offs / ds_max)
         else
            wgt(jprf) = ds(jprf)
         endif
         if (iprf.eq.ipdbg .and. jprf.eq.ista) then
            write(bufout,'(3(a,f10.4))') ' s_offs=',s_offs,', ds=',ds(jprf),', wgt=',wgt(jprf)
            call write_log(1, bufout)
         endif
      enddo

      ! fit circular arc to current data (weights set to zero for points to be excluded)

      if (iprf.eq.ipdbg) idebug = 3
      call fit_circle(nprf, ymeas, zmeas, wgt, iprf, ista, iend, ymid, zmid, phi, dst, rho, idebug)
      idebug = 1
      
      if (idebug.ge.3 .or. iprf.eq.ipdbg) then
         write(bufout,'(a,i4,a,f7.3,2(a,f6.3),a,f7.3,2(a,f6.3))') ' i=',iprf, ': R_est=',1d0/rho(iprf), &
                ', mid=(',ymid(iprf),',',zmid(iprf),'), phi=',phi(iprf),', dst=',dst(iprf),             &
                ', rho=',rho(iprf)
         call write_log(1, bufout)
      endif

   enddo ! iprf

   ! compute output positions

   do iprf = 1, nprf
      prf_grid%y(iprf) = prf_grid%y(iprf) - sin(phi(iprf)) * dst(iprf)
      prf_grid%z(iprf) = prf_grid%z(iprf) + cos(phi(iprf)) * dst(iprf)
   enddo

   ! reorder points at strong zig-zag patterns

!  if (.false.) then
!     dy     = diff(y_out)
!     dz     = diff(z_out)
!     alph   = atan2(dz, dy)
!     dalph  = diff(alph)
!     ix     = find(dalph<-pi)
!     dalph(ix) = dalph(ix) + 2*pi
!     ix     = find(dalph> pi)
!     dalph(ix) = dalph(ix) - 2*pi
!
!     zig_thrs = 135 * pi/180
!     ix_zig = find(abs(dalph(1:end-1))>zig_thrs & abs(dalph(2:end))>zig_thrs)
!     if (~isempty(ix_zig)) then
!        for ix = ix_zig
!           disp(sprintf('Zig-zag at ix=%d: alpha = %4.1f, %4.1f, %4.1f, y = %5.2f, %5.2f, %5.2f, %5.2f', ...
!                               ix, alph(ix+[0:2])*180/pi, y_out(ix+[0:3]) ))
!           ! swap middle two points
!           y_out(ix+[1,2]) = y_out(ix+[2,1])
!           z_out(ix+[1,2]) = z_out(ix+[2,1])
!        enddo
!     endif
!  endif

!  ! estimate surface inclination

!  alph = atan2( diff(z_out), diff(y_out) )
!  ix   = find(alph<-pi)
!  alph(ix) = alph(ix) + 2*pi
!  ix   = find(alph> pi)
!  alph(ix) = alph(ix) - 2*pi
!  alph = [alph(1); 0.5*(alph(1:end-1)+alph(2:end)); alph(end)]

!  ! form output structure, including Curvature

!  prf_out                  = prf_in
!  prf_out.ProfileY         = y_out
!  prf_out.ProfileZ         = z_out
!  prf_out.ProfileAngle     = alph
!  prf_out.ProfileCurvature = rho

   end subroutine make_arcfit

!------------------------------------------------------------------------------------------------------------

   subroutine fit_circle(nprf, xmeas, ymeas, wgt, ifit, ista, iend, ymid, zmid, phi, dst, rho, idebug)
!--purpose: fit circular arc to points xmeas, ymeas using Karimaki-method with weights wgt
      implicit none
!--subroutine parameters:
      integer                    :: nprf, ifit, ista, iend, idebug
      real(kind=8), dimension(:) :: xmeas(nprf), ymeas(nprf), wgt(nprf),                                &
                                    ymid(nprf), zmid(nprf), phi(nprf), dst(nprf), rho(nprf)
!--local variables:
      logical,      parameter    :: use_error_upd = .true.
      real(kind=8), parameter    :: thrs_print_err = 0.005d0
      integer                    :: nmeas, iprf, j, ipiv(3), iwork(3), lwork, info
      real(kind=8) :: r2_meas, dx, dy, Sw, Sx, Sy, Sr2, Sxx, Sxy, Syy, Sxr2, Syr2, Sr2r2,               &
                      Cxx, Cxy, Cyy, Cxr2, Cyr2, Cr2r2, q1, q2, u, phi_est, kappa, delta, sigma,        &
                      rho_est, d_est, R_est, a_est, b_est, Salpha, Sbeta, Sgamma, Sdelta, Salpha2,      &
                      err, Vm1_rho2, Vm1_rhophi, Vm1_phi2, Vm1_rhod, Vm1_phid, Vm1_d2,                  &
                      V_rho2, V_rhophi, V_rhod, V_phi2, V_phid, V_d2, DX2_drho, dX2_dd, drho, dphi, dd
      real(kind=8) :: Vm1(3,3), colsum, anorm, rcond, work(4*3)


   ! ignore positions with negative weights

   nmeas = 0
   do iprf = ista, iend
      if (wgt(iprf).le.0d0) then
         wgt(iprf) = max(0d0, wgt(iprf))
      else
         nmeas = nmeas + 1
      endif
   enddo

   if (idebug.ge.3) then
      write(bufout,'(a,i4,a,i3,3(a,f6.3))') ' ifit=',ifit,': nmeas=',nmeas,', wgt=', wgt(ista), ',',    &
                wgt(ista+1),'...', wgt(iend)
      call write_log(1, bufout)
   endif

   ! non-Iterative Karimaki solution

   if (nmeas.le.2) then
      write(bufout,'(a,i2,a,i4)') 'Warning: nmeas =',nmeas,' d for ifit =',ifit
      call write_log(1, bufout)
   endif

   Sw    = 0d0
   Sx    = 0d0
   Sy    = 0d0
   Sr2   = 0d0
   Sxx   = 0d0
   Sxy   = 0d0
   Syy   = 0d0
   Sxr2  = 0d0
   Syr2  = 0d0
   Sr2r2 = 0d0

   do iprf = ista, iend
      r2_meas = xmeas(iprf)**2 + ymeas(iprf)**2
      Sw    = Sw    + wgt(iprf)
      Sx    = Sx    + wgt(iprf) * xmeas(iprf)
      Sy    = Sy    + wgt(iprf) * ymeas(iprf)
      Sr2   = Sr2   + wgt(iprf) * r2_meas
      Sxx   = Sxx   + wgt(iprf) * xmeas(iprf) * xmeas(iprf)
      Sxy   = Sxy   + wgt(iprf) * xmeas(iprf) * ymeas(iprf)
      Syy   = Syy   + wgt(iprf) * ymeas(iprf) * ymeas(iprf)
      Sxr2  = Sxr2  + wgt(iprf) * xmeas(iprf) * r2_meas
      Syr2  = Syr2  + wgt(iprf) * ymeas(iprf) * r2_meas
      Sr2r2 = Sr2r2 + wgt(iprf) * r2_meas * r2_meas
   enddo

   Cxx   =   Sxx / Sw -    Sx**2 / Sw**2
   Cxy   =   Sxy / Sw - Sx *  Sy / Sw**2
   Cyy   =   Syy / Sw -    Sy**2 / Sw**2
   Cxr2  =  Sxr2 / Sw - Sx * Sr2 / Sw**2
   Cyr2  =  Syr2 / Sw - Sy * Sr2 / Sw**2
   Cr2r2 = Sr2r2 / Sw -   Sr2**2 / Sw**2

   if (idebug.ge.3) then
      write(bufout,'(a,i4,a,6g12.4)') ' ifit=',ifit,': Cxx--Cr2r2=', Cxx, Cxy, Cyy, Cxr2, Cyr2, Cr2r2
      call write_log(1, bufout)
   endif

   q1 = Cr2r2 * Cxy - Cxr2 * Cyr2
   q2 = Cr2r2 * (Cxx - Cyy) - Cxr2**2 + Cyr2**2

   phi_est = 0.5d0 * atan2( 2d0*q1 , q2 )
   kappa   = (sin(phi_est) * Cxr2 - cos(phi_est) * Cyr2) / Cr2r2
   delta   = -kappa * Sr2 / Sw + sin(phi_est) * Sx / Sw - cos(phi_est) * Sy / Sw

   rho_est = 2d0 * kappa / sqrt(1d0 - 4d0 * delta * kappa)
   d_est   = 2d0 * delta / (1d0 + sqrt(1d0 - 4d0 * delta * kappa))

   ! swap signs to get phi decreasing from measurement point 1 to point 2

   if (ifit.le.1) then
      dx = xmeas(ifit+1) - xmeas(ifit)
      dy = ymeas(ifit+1) - ymeas(ifit)
   else
      dx = xmeas(ifit) - xmeas(ifit-1)
      dy = ymeas(ifit) - ymeas(ifit-1)
   endif
   if (dx * cos(phi_est) + dy * sin(phi_est).gt.0d0) then
      rho_est = -rho_est
      d_est   = -d_est
      phi_est = phi_est + pi
   endif
   R_est   = 1d0 / rho_est

   a_est = (R_est + d_est) *  sin(phi_est)
   b_est = (R_est + d_est) * -cos(phi_est)

   if (idebug.ge.3) then
      write(bufout,'(a,i4,3(a,2g12.4))') ' ifit=',ifit,': a,b=',a_est,b_est,', R,rho=',R_est,rho_est,   &
                ', phi,d=',phi_est, d_est
      call write_log(1, bufout)
   endif

   ! Karimaki error estimates

   if (use_error_upd .and. nmeas.ge.4) then

      u       = 1d0 + rho_est * d_est
      Salpha  = Sx * sin(phi_est) - Sy * cos(phi_est)
      Sbeta   = Sx * cos(phi_est) + Sy * sin(phi_est)
      Sgamma  = (sin(phi_est)**2 - cos(phi_est)**2) * Sxy + sin(phi_est) * cos(phi_est) * (Sxx - Syy)
      Sdelta  = Sxr2 * sin(phi_est) - Syr2 * cos(phi_est)
      Salpha2 = Sxx * sin(phi_est)**2 - 2*sin(phi_est)*cos(phi_est)*Sxy + cos(phi_est)**2*Syy

      Vm1_rho2   = 0.25d0 * Sr2r2 - d_est * ( Sdelta - d_est * (Salpha2 + 0.5d0 * Sr2                   & 
                        - d_est * ( Salpha - 0.25d0 * d_est * Sw ) ) )
      Vm1_rhophi = -u * ( 0.5d0 * ( cos(phi_est) * Sxr2 + sin(phi_est) * Syr2 )                         &
                        - d_est * ( Sgamma - 0.5d0 * d_est * Sbeta ) )
      Vm1_phi2   = u**2 * ( cos(phi_est)**2 * Sxx + sin(2d0*phi_est) * Sxy + sin(phi_est)**2 * Syy )
      Vm1_rhod   = rho_est * (-0.5d0 * Sdelta + d_est * Salpha2) + 0.5d0 * u * Sr2                      &
                        - 0.5d0 * d_est * ( (2d0 *u  + rho_est * d_est ) * Salpha - u * d_est * Sw )
      Vm1_phid   = u * ( rho_est * Sgamma - u * Sbeta )
      Vm1_d2     = rho_est * ( rho_est * Salpha - 2d0 * u * Salpha ) - u**2 * Sw
   
   
      Vm1 = reshape( (/ Vm1_rho2, Vm1_rhophi, Vm1_rhod, Vm1_rhophi, Vm1_phi2, Vm1_phid, Vm1_rhod,       &
                        Vm1_phid, Vm1_d2 /), shape(Vm1) )

      ! Determine ||A||_1 = \max_j \sum_i |a_ij| needed by dgecon

      anorm = 0d0
      do j = 1, 3
         colsum = abs(Vm1(1,j)) + abs(Vm1(2,j)) + abs(Vm1(3,j))
         anorm = max(anorm, colsum)
      enddo

      ! compute LU factorization needed for condition number

      call dgetrf(3, 3, Vm1, 3, ipiv, info)

      if (info.ne.0) then
         write(bufout,'(a,i4)') ' Something wrong in dgetrf, info=',info
         call write_log(1, bufout)
      endif

      ! determine condition number

      call dgecon('1', 3, Vm1, 3, anorm, rcond, work, iwork, info)

      if (info.ne.0) then
         write(bufout,'(a,i4)') ' Something wrong in dgecon, info=',info
         call write_log(1, bufout)
      endif

      if (rcond.lt.1d-10) then
         write(bufout,'(a,g12.4,a,i4)') ' Error: matrix Vm1 is singular (rcond=',rcond,'), nmeas=', nmeas
         call write_log(1, bufout)
      endif

      ! compute inverse matrix V = inv(Vm1)

      lwork = 12
      call dgetri( 3, Vm1, 3, ipiv, work, lwork, info)

      if (info.ne.0) then
         write(bufout,'(a,i4)') ' Something wrong in dgetri, info=',info
         call write_log(1, bufout)
      endif

      V_rho2   = Vm1(1,1)
      V_rhophi = Vm1(1,2)
      V_rhod   = Vm1(1,3)
      V_phi2   = Vm1(2,2)
      V_phid   = Vm1(2,3)
      V_d2     = Vm1(3,3)

      ! compute improved estimates: 1 step of Newton method

      sigma    = -rho_est * Sdelta + 2d0 * u * Salpha2 - d_est * (1d0 + u) * Salpha
      dX2_drho = d_est * sigma
      dX2_dd   = rho_est * sigma

      drho = -0.5 * (V_rho2  *dX2_drho + V_rhod*dX2_dd)
      dphi = -0.5 * (V_rhophi*dX2_drho + V_phid*dX2_dd)
      dd   = -0.5 * (V_rhod  *dX2_drho + V_d2  *dX2_dd)

      err = max(abs(drho), abs(dd))
      err = max(err, abs(dphi))
      if (err.gt.thrs_print_err) then
         write(bufout,'(3(a,g12.4))') ' Error estimates: drho =',drho,', dd =',dd,', dphi =', dphi
         call write_log(1, bufout)
      endif

      rho_est = rho_est - drho
      phi_est = phi_est - dphi
      d_est   = d_est - dd

      R_est   = 1d0 / rho_est
      a_est  = (R_est + d_est) *  sin(phi_est)
      b_est  = (R_est + d_est) * -cos(phi_est)

   endif ! use_error_upd

   ! determine misfit d_i of measurement data to fitted circle

   ! th_meas = atan2( ymeas-b_est, xmeas-a_est )
   ! d_xi = xmeas - (a_est + R_est * cos(th_meas))
   ! d_yi = ymeas - (b_est + R_est * sin(th_meas))
   ! d_i  = d_xi .* cos(th_meas) + d_yi .* sin(th_meas)

   ! convert from Karimaki-convention to our own convention

   phi_est = phi_est - pi

   ymid(ifit) = a_est
   zmid(ifit) = b_est
   phi(ifit)  = phi_est
   rho(ifit)  = rho_est
   dst(ifit)  = d_est

   end subroutine fit_circle

!------------------------------------------------------------------------------------------------------------

end module m_arcfit
