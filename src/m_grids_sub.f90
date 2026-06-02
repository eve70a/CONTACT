!------------------------------------------------------------------------------------------------------------
! m_grids_sub - implementation of selected subroutines for m_grids
!
! Copyright 2016-2026 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
submodule (m_grids) m_grids_sub
   use m_globals
   use m_markers
   use m_ptrarray
   implicit none

!------------------------------------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------------------------------------

module subroutine grid_make_arclength(g, use_correction, nkink, ikinks, idebug, my_ierror)
!--purpose: compute the summed chord length s for a 1-d grid (wheel or rail profile)
   implicit none
!--subroutine arguments:
   type(t_grid)               :: g
   logical,      intent(in)   :: use_correction
   integer,      intent(in)   :: nkink, idebug
   integer,      intent(in)   :: ikinks(nkink)
   integer,      intent(out)  :: my_ierror
!--local variables:
   real(kind=8), parameter :: cos_thrs = 0.5d0, frac_ds_min = 0.1d0
   integer,      parameter :: kmax = 10
   logical            :: has_large_da
   integer            :: nseg, iseg, isg_min, ik, ip, ipmin, sub_ierror
   real(kind=8)       :: cos_da, cos_min, cos_tavg, ty_avg, tz_avg, dst, dstmin, s0
   real(kind=8), dimension(:), allocatable :: dy, dz, ds

   my_ierror = 0

   if (nkink.lt.2) then
      call write_log(' Internal error (make_arclength): need nkink>=2 section start/end-points')
      my_ierror = -2
   endif
   if (ikinks(nkink).ne.g%ntot) then
      call write_log(' Internal error (make_arclength): need last end-point at g%ntot')
      my_ierror = -3
   endif

   ! allocate work arrays

   nseg = g%ntot - 1
   allocate(dy(nseg), dz(nseg), ds(nseg))

   ! compute cumulative sum of distances between successive points

   do iseg = 1, nseg
      dy(iseg) = g%y(iseg+1) - g%y(iseg)
      dz(iseg) = g%z(iseg+1) - g%z(iseg)
      ds(iseg) = sqrt(dy(iseg)**2 + dz(iseg)**2)
   enddo

   ! check for large angles dalpha (small cosine) between consecutive segments, corrected chord length

   if (.not.use_correction) then

      has_large_da = .false.

   else

      has_large_da = .false.
      cos_min = 1d0
      isg_min = 0

      ik = 1
      iseg = 1
      do while((.not.has_large_da .or. idebug.ge.1) .and. iseg.lt.nseg)
         if (iseg+1.ge.ikinks(ik+1)) then
            ik = ik + 1         ! ignore large angle between different spline sections
         else
            cos_da = (dy(iseg) * dy(iseg+1) + dz(iseg) * dz(iseg+1)) / (ds(iseg) * ds(iseg+1))
            has_large_da = (has_large_da .or. cos_da.lt.cos_thrs)
            if (cos_da.lt.cos_min) then
               cos_min = cos_da
               isg_min = iseg
            endif
         endif
         iseg = iseg + 1
      enddo

      if (idebug.ge.2 .and. has_large_da) then
         call write_log(' make_arclength: detected large angle between segments...')
      endif
      if (idebug.ge.2) then
         write(bufout,'(a,f5.2,a,i0,a,f6.1,a)') ' min(cos) =',cos_min,' at iseg=',isg_min,', angles <=', &
                           acos(cos_min)*180d0/pi,' deg'
         call write_log(1, bufout)
      endif
   endif

   ! for wiggly signal, correct segment lengths using the average heading

   if (has_large_da) then
      do iseg = 1, nseg

         ! determine average heading (ty_avg,tz_avg) for segment iseg

         call get_average_heading( g%y, g%z, g%ntot, iseg, kmax, ty_avg, tz_avg )

         ! compute angle davg between segment and average heading

         cos_tavg = (dy(iseg) * ty_avg + dz(iseg) * tz_avg) / ds(iseg)

         ! scale segment lengths ds to the length in the average heading direction

         ds(iseg) = ds(iseg) * max(frac_ds_min, cos_tavg)
      end do
   endif

   ! re-allocate s-array at the appropriate size

   call reallocate_arr(g%s_prf, g%ntot)

   ! compute cumulative sum of distances between successive points

   g%s_prf(1) = 0d0
   do iseg = 1, nseg
      g%s_prf(iseg+1) = g%s_prf(iseg) + ds(iseg)
   enddo

   ! select reference position in profile where s:=0

   if (is_strictly_monotonic(g%ntot, g%y)) then

      ! interpolate [g%y, g%s_prf] to y = 0

      call interp_1d_to_scalar( g%ntot, g%y, g%s_prf, 0d0, s0, sub_ierror )
      if (my_ierror.eq.0) my_ierror = sub_ierror

   elseif (.false.) then

      ! shifting can be omitted cf. read_profile.m

      s0 = 0d0

   else

      ! select point (g%y, g%z) closest to origin (0, 0)

      ipmin  = -1
      dstmin = 1e12
      do ip = 1, g%ntot
         dst = g%y(ip)**2 + g%z(ip)**2
         if (dst.le.dstmin) then
            ipmin  = ip
            dstmin = dst
         endif
      enddo
      s0 = g%s_prf(ipmin)

   endif

   if (.false.) then
      write(bufout,'(a,f8.3)') ' grid_make_arclength: s0=', s0
      call write_log(1, bufout)
   endif

   ! shift s := s - s0

   do ip = 1, g%ntot
      g%s_prf(ip) = g%s_prf(ip) - s0
   enddo

   deallocate(dy, dz, ds)

end subroutine grid_make_arclength

!------------------------------------------------------------------------------------------------------------

subroutine get_average_heading( y, z, npnt, iseg, kmax, ty_avg, tz_avg )
!--purpose: compute the average heading of curve (y,z) at segment iseg
   implicit none
!--subroutine arguments:
   integer,      intent(in)   :: npnt, iseg, kmax
   real(kind=8), intent(in)   :: y(npnt), z(npnt)
   real(kind=8), intent(out)  :: ty_avg, tz_avg
!--local variables:
   integer            :: nseg, ip0, ip1, j, k, np
   real(kind=8)       :: dy, dz, ds, ty_est, tz_est, y_avg, z_avg, syy, syz, szz, alph_avg
 
   ! using orthogonal projection to fit the best line through [ y(i+k), z(i+k) ], k = [-kmax:kmax+1]
   ! with center (y_avg,z_avg) and tangent (ty_avg,tz_avg)

   nseg = npnt - 1

   if (nseg.lt.2*kmax+3) then
      write(bufout,'(2(a,i0),a)') ' Internal error: avg_heading: cannot use kmax=',kmax,' for n=',nseg, &
                        ' segments'
      call write_log(1, bufout)
      return
   endif

   ! determine initial estimate: average heading across window:  x(i+k+1) - x(i-k)

   ip0    = max(   1, iseg-kmax)
   ip1    = min(npnt, iseg+kmax+1)
   dy     = y(ip1) - y(ip0)
   dz     = z(ip1) - z(ip0)
   ds     = sqrt(dy**2 + dz**2)
   ty_est = dy / ds
   tz_est = dz / ds

   ! determine centroid of point cloud, center position (y_avg, z_avg)

   y_avg  = 0d0
   z_avg  = 0d0
   np     = 0

   do k = -kmax, kmax + 1
      j = iseg + k
      if (j.ge.1 .and. j.le.npnt) then
         np = np + 1
         y_avg = y_avg + y(j)
         z_avg = z_avg + z(j)
      endif
   enddo

   y_avg = y_avg / real(max(1,np))
   z_avg = z_avg / real(max(1,np))

   ! determine sum of squared distances syy, szz, syz

   syy  = 0d0
   syz  = 0d0
   szz  = 0d0
   np   = 0

   do k = -kmax, kmax + 1
      j = iseg + k
      if (j.ge.1 .and. j.le.npnt) then
         np = np + 1
         syy = syy + y(j)**2
         szz = szz + z(j)**2
         syz = syz + y(j) * z(j)
      endif
   enddo

   syy = syy / real(max(1,np)) - y_avg**2
   szz = szz / real(max(1,np)) - z_avg**2
   syz = syz / real(max(1,np)) - y_avg * z_avg

   ! determine average heading alph_avg based on orthogonal projection

   alph_avg = 0.5d0 * atan2( 2d0*syz, syy-szz )
   ty_avg  = cos(alph_avg)
   tz_avg  = sin(alph_avg)

   ! flip sign where t_avg goes against first estimate

   if (ty_avg*ty_est + tz_avg*tz_est .lt. 0d0 ) then
      ty_avg = -ty_avg
      tz_avg = -tz_avg
   endif

end subroutine get_average_heading

!------------------------------------------------------------------------------------------------------------

end submodule m_grids_sub
