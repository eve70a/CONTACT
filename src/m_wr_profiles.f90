!------------------------------------------------------------------------------------------------------------
! m_wr_profiles - place w/r profiles in track coordinate system (module 1)
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_wr_profiles

use m_wrprof_data
implicit none
private

public  wr_set_rail_marker
public  wr_set_wheel_marker
public  find_gauge_meas_pt

contains

!------------------------------------------------------------------------------------------------------------

   subroutine wr_set_rail_marker (trk, ic, idebug)
!--purpose: set the rail profile marker, position in the track coordinate system
      implicit none
!--subroutine arguments:
      type(t_trackdata)       :: trk
      type(t_ic)              :: ic
      integer                 :: idebug
!--local variables:
      type(t_rail),   pointer :: my_rail
      integer                 :: is_right
      real(kind=8)            :: sgn, ygauge, yvampr, zgauge, zmin, rail_y, rail_z, rail_phi
      type(t_grid),   pointer :: p_inp
      character(len=5)        :: nam_rail

      if (idebug.ge.3) then
         write(bufout,'(/a)') ' --- Start subroutine wr_set_rail_marker ---'
         call write_log(2, bufout)
      endif

      ! set pointer to the active rail in the current configuration

      my_rail     => trk%rai
      if (ic%config.eq.0 .or. ic%config.eq.4) then
         nam_rail = 'left'
         is_right = 0
      else
         nam_rail = 'right'
         is_right = 1
      endif

      ! set the sign to -1 for left and +1 for right rail/wheel combination

      sgn = 2 * is_right - 1

      ! check that the input rail profile is specified(?), if not: abort

      p_inp => my_rail%prr%grd_data

      ! print information on the input-profile

      if (idebug.ge.3) call grid_print(p_inp, nam_rail//' rail', idebug-1)

      ! if needed, determine the gauge point in the rail profile

      if (trk%gauge_height.gt.0d0) then

         ! find the left-most point in the canted profile within gauge measuring height
         ! note: ygauge,zmin are along the track y,z-directions

         call find_gauge_meas_pt(p_inp, trk%cant_angle, trk%gauge_height, ygauge, zgauge, yvampr,       &
                                 zmin, idebug)

         if (idebug.ge.2) then
            write(bufout,123) 'rail gauge point yg=',ygauge,', offset rail->track=',                    &
                   sgn*(trk%track_gauge/2 - ygauge)
            call write_log(1, bufout)
 123        format(1x, 3(a, f12.6, :))
         endif

      endif

      ! create the rail profile marker

      call marker_init(my_rail%m_trk)

      ! rotate the marker by cant and rail roll (irregularity, mirrored for left rails)

      rail_phi = -trk%cant_angle + sgn * my_rail%roll
      call marker_roll(my_rail%m_trk, rail_phi, 0d0, 0d0)

      ! shift the marker according to the rail position

      if (trk%gauge_height.le.0d0) then

         ! gauge_height <= 0: absolute rail placement + track deviations

         rail_y    =  trk%rail_y0 + sgn * my_rail%dy
         rail_z    =  trk%rail_z0 +       my_rail%dz

      else

         ! gauge_height >  0: rail placement using half gauge width + deviations, subject to mirroring

         rail_y    =  trk%track_gauge/2d0 - ygauge + sgn * my_rail%dy

         ! wheelset on track: shift vertically to put highest point on z=0 + deviations
         ! wheelset on roller rig: no vertical shifting

         if (ic%config.le.1) then
            rail_z =                      - zmin   +       my_rail%dz
         else
            rail_z =                                       my_rail%dz
         endif

      endif

      call marker_shift(my_rail%m_trk, 0d0, rail_y, rail_z)

      if (idebug.ge.2) then
         write(bufout,'(2(a,f12.6))') ' rail shift: y =',rail_y,', z =',rail_z
         call write_log(1, bufout)
      endif

      if (idebug.ge.4) call write_log('--- end subroutine wr_set_rail_marker ---')

   end subroutine wr_set_rail_marker

!------------------------------------------------------------------------------------------------------------

   subroutine wr_set_wheel_marker (ws, ic, idebug)
!--purpose: set the wheel marker in ws coordinates, wheelset marker in track coordinate system
      implicit none
!--subroutine arguments:
      type(t_wheelset)        :: ws
      type(t_ic)              :: ic
      integer                 :: idebug
!--local variables:
      type(t_wheel),  pointer :: my_wheel
      integer                 :: is_right
      character(len=5)        :: nam_wheel
      real(kind=8)            :: whl_x, whl_y, whl_z, sgn
      type(t_grid),   pointer :: p_inp

      if (idebug.ge.3) then
         write(bufout,'(/a)') ' --- Start subroutine wr_set_wheel_marker ---'
         call write_log(2, bufout)
      endif

      ! set pointer to the active wheel in the current configuration

      my_wheel    => ws%whl
      if (ic%config.eq.0 .or. ic%config.eq.4) then
         nam_wheel = 'left'
         is_right  = 0
      else
         nam_wheel = 'right'
         is_right  = 1
      endif

      p_inp => my_wheel%prw%grd_data

      ! print information on the input-profile

      if (idebug.ge.3) call grid_print(p_inp, nam_wheel//' wheel', idebug-1)

      ! set the sign to -1 for left and +1 for right rail/wheel combination

      sgn = 2 * is_right - 1

      ! define the relation from wheel profile to wheelset coordinates
      ! including flexible wheel-set deviations: axle bending, wheel deformation

      whl_x =                                            my_wheel%dx
      whl_y = ws%flback_dist/2d0 - ws%flback_pos + sgn * my_wheel%dy
      whl_z = ws%nom_radius                      +       my_wheel%dz

      call marker_init(my_wheel%m_ws)
      ! note: rotate first, then shift to desired position
      ! note: the pitch-variation dpitch is excluded: the wheel-marker is the lowest point on the wheel
      call marker_rotate(my_wheel%m_ws, sgn*my_wheel%droll, sgn*my_wheel%dyaw, 0d0, 0d0, 0d0, 0d0)
      call marker_shift(my_wheel%m_ws, whl_x, whl_y, whl_z)

      if (idebug.ge.3) then
         write(bufout,123) 'offset wheel->ws x:', whl_x, ', y:', sgn*whl_y,', z:', whl_z
         call write_log(1, bufout)
 123     format(1x, 3(a, f12.6, :))
      endif

      ! define the relation from wheelset to track coordinates, using mirroring for left wheel
      ! note: the wheelset marker excludes the pitch rotation

      call marker_init(ws%m_trk)
      call marker_rotate(ws%m_trk, sgn*ws%roll, sgn*ws%yaw, 0d0, 0d0, 0d0, 0d0)
      call marker_shift(ws%m_trk, 0d0, sgn*ws%y, ws%z-ws%nom_radius)

      if (idebug.ge.5) then
         call marker_print(my_wheel%m_ws, 'rw(ws)', 2)
         call marker_print(ws%m_trk, 'ws(tr)', 2)
         call marker_print(marker_2glob( my_wheel%m_ws, ws%m_trk ), 'rw(tr)', 2)
      endif

      if (idebug.ge.4) call write_log('--- end subroutine wr_set_wheel_marker ---')

   end subroutine wr_set_wheel_marker

!------------------------------------------------------------------------------------------------------------

   subroutine find_gauge_meas_pt (prr_orig, cant, gauge_height, ygauge, zgauge, yvampr, zmin, idebug)
!--purpose: determine the shift [ygauge,-zmin] in track [y,z] needed to position the canted rail
!           profile such that it touches the plane z=0 and has z > gauge_height for all profile points
!           y < ygauge.
!--note: assuming that the rail profile has the track center at negative y (right rail) and that cant
!        is positive if the normal is rotated towards the track center.
      implicit none
!--subroutine arguments:
      type(t_grid)              :: prr_orig
      real(kind=8), intent(in)  :: cant
      real(kind=8), intent(in)  :: gauge_height    ! gauge measuring height w.r.t. highest point on rail
      integer,      intent(in)  :: idebug
      real(kind=8), intent(out) :: ygauge          ! inner-most point on rail profile within gauge meas.range
      real(kind=8), intent(out) :: zgauge          ! rail profile z-value at gauge point ygauge
      real(kind=8), intent(out) :: yvampr          ! rail profile y-value at gauge meas.height
      real(kind=8), intent(out) :: zmin            ! lowest z-coordinate on canted rail profile
!--local variables:
      integer                   :: iminy, i1st, ilast, ir, j, ierror
      real(kind=8)              :: sout, xout, yout
      type(t_grid)              :: prr_cant, prr_trim

      if (prr_orig%y(1).ge.prr_orig%y(prr_orig%ntot)) then
         call write_log(' Internal Error: y-coordinates of rail profile should be ascending.')
         do j = 1, 6
            ir = j
            if (j.ge.4) ir = prr_orig%ntot - 6 + j
            write(bufout,'(a,i5,2(a,f9.3))') ' i=',ir,': yr(i)=',prr_orig%y(ir),', zr(i)=',prr_orig%z(ir)
            call write_log(1, bufout)
            if (j.eq.3) call write_log('     ...')
         enddo
         call abort_run()
      endif

      !  1. create canted rail profile
      ! note: for the right rail positive cant 1:40 corresponds to negative rotation alpha = -0.024995

      call grid_copy(prr_orig, prr_cant, with_spline=.true.)
      call cartgrid_roll(prr_cant, -cant, 0d0, 0d0)

      if (idebug.ge.3) then
         write(bufout,123) 'rail(1): y =',prr_orig%y(1),', z =',prr_orig%z(1)
         call write_log(1, bufout)
      endif
      if (idebug.ge.4) then
         write(bufout,123) 'cant angle=', cant
         call write_log(1, bufout)
      endif
      if (idebug.ge.3) then
         write(bufout,123) 'rail(1): y''=',prr_cant%y(1),', z''=',prr_cant%z(1)
         call write_log(1, bufout)
 123     format(1x, 3(a, f12.6, :))
      endif

      !  2. find the point with minimum z-coordinate, vertical shift zmin

      ! imin = idmin(prr_cant%ntot, prr_cant%z(1:), 1)
      ! zmin = prr_cant%z(imin)

      call spline_get_xyz_at_minz(prr_cant%spl, ierror, sout, xout, yout, zmin)

      if (idebug.ge.3) then
         write(bufout,123) 'minimum z on rail=',zmin,' at yr=', yout
         call write_log(1, bufout)
      endif

      !  3. find the first point in the profile that has z-zmin <= gauge_height
      !      and the left-most point with z-zmin in [0,gauge_height]

      !      cf. Simpack help section D.19.3.3: track gauge is minimum y-distance in z-range [0,14mm]
      !      points are considered on the inside face of the rail only, with with s <= s(zmin)

      ir    =  1
      i1st  = -1
      iminy =  1
      ygauge = 1d9
      do while( prr_cant%s_prf(ir).lt.sout .and. ir.lt.prr_cant%ntot)

         if (prr_cant%z(ir)-zmin.le.gauge_height) then

            ! first point where z-zmin <= gauge_height

            if (i1st.lt.0) i1st = ir

            ! minimum y-value over points with z in [0, gauge_height]

            if (prr_cant%y(ir).lt.ygauge) then
               ygauge = prr_cant%y(ir)
               iminy  = ir
            endif

         endif

         ir = ir + 1

      enddo

      ilast = ir

      if (idebug.ge.3) then
         write(bufout,'(3(a,i5))') ' points considered for gauge computation: i1st=',i1st,           &
                        ', iminy=', iminy,', ilast=',ilast
         call write_log(1, bufout)
         write(bufout,'(4(a,f10.6))') ' z(1)=', prr_cant%z(1)-zmin,', z(i1st)=',                     &
                     prr_cant%z(i1st)-zmin, ', z(iminy)=', prr_cant%z(iminy)-zmin,', z(ilast)=',     &
                     prr_cant%z(ilast)-zmin
         call write_log(1, bufout)
      endif

      if (iminy.eq.1) then

         ! 4.a left-most point is first point on profile

         if (idebug.ge.1) then
            call write_log(' WARNING: gauge computation failed for rail profile.')
            write(bufout,'(2(a,f9.3),a)') '          no points could be found with y <= y(1) =',     &
                        prr_cant%y(1), ' and z >=', gauge_height,' (gauge height)'
            call write_log(1, bufout)

            write(bufout,'(2(a,f9.3))') '          using 1st profile point, y(1) =', prr_cant%y(1),  &
                        ', z(1) =', prr_cant%z(1)-zmin
            call write_log(1, bufout)
         endif

         ygauge = prr_cant%y(1)
         zgauge = prr_cant%z(1)
         yvampr = ygauge

      else

         ! compute y-value at gauge-height using spline interpolation

         call grid_trim(prr_cant, prr_trim, 1, 1, 1, ilast)
         call spline_get_xy_at_z( prr_trim%spl, zmin+gauge_height, ierror, xout, yvampr, 999d0)

         if (iminy.eq.i1st) then

            ! 4.b gauge point is found at height z-zmin == gauge_height, using interpolation

            if (idebug.ge.3) call write_log(' gauge point at z=zgauge, using interpolation...')
            ygauge = yvampr
            zgauge = zmin + gauge_height

         else

            ! 4.c gauge point lies inside the interval z-zmin in [0, gauge_height],
            !     compute vertical slope dy/ds = 0

            if (idebug.ge.3) call write_log(' gauge point at zgauge<z<0, using slope dy/ds=0...')

            ! call spline_set_debug(5)
            call spline_get_xyz_at_miny(prr_cant%spl, i1st, ilast, ierror, sout, xout, ygauge, zgauge)
            ! call spline_set_debug(0)

            if (idebug.ge.1) then
               write(bufout,'(2a,f8.3,a,f7.3)') ' INFO: left-most rail point found above gauge ',    &
                        'height: zr=', zgauge-zmin, ' at yr=', ygauge
               call write_log(1, bufout)
            endif
         endif

      endif

      ! clean up local variables

      call grid_destroy(prr_cant)
      call grid_destroy(prr_trim)

   end subroutine find_gauge_meas_pt

!------------------------------------------------------------------------------------------------------------

end module m_wr_profiles
