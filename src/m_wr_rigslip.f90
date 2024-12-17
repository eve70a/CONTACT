!------------------------------------------------------------------------------------------------------------
! m_wr_rigid_slip - compute rigid slip for one contact patch in a w/r contact problem (module 1)
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_wr_rigslip

use m_hierarch_data
use m_wrprof_data
implicit none
private

public  wr_rigid_slip
private wr_creep_cref
private wr_rigslip_wrsurf
private extract_creepages
private project_midplane

contains

!------------------------------------------------------------------------------------------------------------

   subroutine wr_rigid_slip(ic, ws, trk, dqrel, icp, cp, gd, idebug)
!--purpose: compute the creepages at the contact reference and/or the rigid slip on the contact grid
      implicit none
!--subroutine arguments:
      integer               :: idebug, icp
      real(kind=8)          :: dqrel
      type(t_ic)            :: ic
      type(t_wheelset)      :: ws
      type(t_trackdata)     :: trk
      type(t_cpatch)        :: cp
      type(t_probdata)      :: gd
!--local variables:
      logical               :: use_creep_cref

      if (idebug.ge.2) call write_log(' --- Start subroutine wr_rigid_slip ---')

      use_creep_cref = (ic%discns1_eff.le.2 .or. (ic%discns1_eff.ge.5 .and. ic%discns1_eff.le.9))

      ! compute the rigid slip according to the method of choice

      if (use_creep_cref) then

         ! compute creepages at the contact reference point

         call wr_creep_cref(ic, ws, trk, dqrel, cp, gd, idebug)

      elseif (.true.) then

         ! compute rigid slip at the contact grid points, using the actual wheel & rail surfaces

         call wr_rigslip_wrsurf(ic, ws, trk, dqrel, icp, cp, gd, idebug)

      else

         gd%ic%rztang   = 0
         gd%kin%cksi    = 0d0
         gd%kin%ceta    = 0d0
         gd%kin%cphi    = 0d0

      endif
      if (idebug.ge.4) call write_log('--- end subroutine wr_rigid_slip ---')

   end subroutine wr_rigid_slip

!------------------------------------------------------------------------------------------------------------

   subroutine wr_creep_cref(ic, ws, trk, dqrel, cp, gd, idebug)
!--purpose: compute the creepages at the contact reference point
      implicit none
!--subroutine arguments:
      integer               :: idebug
      type(t_ic)            :: ic
      type(t_wheelset)      :: ws
      type(t_trackdata)     :: trk
      type(t_cpatch)        :: cp
      type(t_probdata)      :: gd
      real(kind=8)          :: dqrel
!--local variables:
      integer               :: is_right
      real(kind=8)          :: r_loc, w_loc, sgn, vs_est
      type(t_marker)        :: m_rol, mref_r, mref_rol, mw_trk, mq_w, mq_ws, mq_cp, mq_trk
      type(t_vec)           :: rol_tvel_trk, rol_rvel_trk, rail_tvel_rol, rail_rvel_rol, rail_tvel_trk, &
                               rail_rvel_trk, vp_tvel_rol, vp_rvel_rol, vp_tvel_trk, vp_rvel_trk,       &
                               vp_tvel_cp, vp_rvel_cp
      type(t_vec)           :: ws_rvel_ws, ws_tvel_trk, ws_rvel_trk, whl_tvel_ws, whl_rvel_ws,          &
                               vq_tvel_ws, vq_rvel_ws, vq_tvel_trk, vq_rvel_trk, vq_tvel_cp,            &
                               vq_rvel_cp, pitch_tvel_trk
      character(len=5)        :: nam_side 

      if (idebug.ge.3) call write_log(' --- Start subroutine wr_creep_cref ---')

      associate(my_rail   => trk%rai, my_wheel  => ws%whl)
      if (ic%is_left_side()) then
         nam_side  = 'left'
         is_right  = 0
      else
         nam_side  = 'right'
         is_right  = 1
      endif

      ! set the sign to -1 for left and +1 for right rail/wheel combination

      sgn = 2 * is_right - 1

      !------------------------------------------------------------------------------------------------------
      ! compute translation & rotation velocity of point P=mref on rail/roller

      if (.not.ic%is_roller()) then

         ! convert position of contact point P on rail from track to rail profile coordinates

         mref_r   = marker_2loc(cp%mref, my_rail%m_trk)

         ! define velocity vectors 'tvel' and 'rvel' for rail profile (origin O) w.r.t track coordinates

         rail_tvel_trk = vec( 0d0, sgn*my_rail%vy, my_rail%vz )
         rail_rvel_trk = vec( sgn*my_rail%vroll, 0d0, 0d0 )

         ! call write_log(' velocity of O_rail wrt O_trk')
         ! call vec_print(rail_tvel_trk, 'tvel_rail(trk)', 2)
         ! call vec_print(rail_rvel_trk, 'rvel_rail(trk)', 2)

         ! compute velocity vectors 'tvel' and 'rvel' of the contact point P in track coordinates

         vp_tvel_trk = vec_veloc2glob(rail_tvel_trk, rail_rvel_trk, my_rail%m_trk, mref_r, vec_zero())
         vp_rvel_trk = rail_rvel_trk

         ! call write_log(' position of O_rail and O_cp wrt O_trk')
         ! call marker_print(my_rail%m_trk, 'O_rail(trk)', 2)
         ! call marker_print(cp%mref,       'O_cp(trk)', 2)

         ! call write_log(' velocity of O_cp wrt O_trk')
         ! call vec_print(vp_tvel_trk, 'tvel_cp(trk)', 2)

         ! transform velocity vectors 'tvel' and 'rvel' from track to contact coordinates

         vp_tvel_cp = cp%mref%rot .transp. vp_tvel_trk
         vp_rvel_cp = cp%mref%rot .transp. vp_rvel_trk

         ! call write_log(' velocity of P_cp wrt O_cp')
         ! call vec_print(vp_tvel_cp, 'tvel_cp(cp)', 2)
         ! call vec_print(vp_rvel_cp, 'rvel_cp(cp)', 2)

      else

         ! create marker for roller center of mass w.r.t. track coordinates

         m_rol%o   = vec( 0d0, 0d0, trk%nom_radius )
         m_rol%rot = rotmat_identity()

         ! convert position of contact point P on roller to rail coordinates & roller coordinates 

         mref_r    = marker_2loc(cp%mref, my_rail%m_trk)
         mref_rol  = marker_2loc(cp%mref, m_rol)

         ! get local rolling radius: z-position of the contact point on the roller

         r_loc     = mref_rol%z()

         ! define velocity vectors 'tvel' and 'rvel' for rail profile origin O in roller coordinates
         !    these are equal to the velocities in track coordinates because R_{rol(trk)} == I

         rail_tvel_rol = vec( 0d0, sgn*my_rail%vy, my_rail%vz )
         rail_rvel_rol = vec( sgn*my_rail%vroll, 0d0, 0d0 )

         ! call write_log(' velocity of O_rail wrt O_rol')
         ! call vec_print(rail_tvel_rol, 'tvel_rail(rol)', 2)
         ! call vec_print(rail_rvel_rol, 'rvel_rail(rol)', 2)

         ! compute velocity vectors 'tvel' and 'rvel' of the contact reference point P in roller coordinates
         ! v_P(rol) = v_r(rol) + omega_r(rol) x ( R_r(rol) * m_P(r) ) + R_r(rol) * v_P(r)

         associate( rail_rot_rol => my_rail%m_trk%rot ) ! note: R_r(rol) == R_r(trk)
         vp_tvel_rol = vec_veloc2glob(rail_tvel_rol, rail_rvel_rol, rail_rot_rol, mref_r, vec_zero())
         vp_rvel_rol = rail_rvel_rol
         end associate

         ! call write_log(' position of P = O_cp wrt O_r')
         ! call marker_print(mref_r,        'O_cp(r)', 2)

         ! call write_log(' velocity of P = O_cp wrt O_rol')
         ! call vec_print(vp_tvel_rol, 'tvel_cp(rol)', 2)

         ! define velocity vectors 'tvel' and 'rvel' of roller in track coordinates

         rol_tvel_trk = vec_zero()
         rol_rvel_trk = m_rol%rot * vec( 0d0, trk%vpitch_rol, 0d0 )

         ! compute velocity vectors 'tvel' and 'rvel' of the contact reference point P in track coordinates
         ! v_P(trk) = v_rol(trk) + omega_rol(trk) x ( R_rol(trk) * m_P(rol) ) + R_rol(trk) * v_P(rol)

         vp_tvel_trk  = vec_veloc2glob(rol_tvel_trk, rol_rvel_trk, m_rol%rot, mref_rol, vp_tvel_rol)
         vp_rvel_trk  = rol_rvel_trk + m_rol%rot * vp_rvel_rol

         ! compute velocity vectors 'tvel' and 'rvel' of the contact reference point P in contact coordinates

         vp_tvel_cp = cp%mref%rot .transp. vp_tvel_trk
         vp_rvel_cp = cp%mref%rot .transp. vp_rvel_trk

         ! call write_log(' velocity of P_cp wrt O_cp')
         ! call vec_print(vp_tvel_cp, 'tvel_cp(cp)', 2)
         ! call vec_print(vp_rvel_cp, 'rvel_cp(cp)', 2)

      endif

      !------------------------------------------------------------------------------------------------------
      ! compute translation & rotation velocity of point Q on the wheel, at 'pen' above mref on the rail

      ! convert position of contact point Q on wheelset to wheel profile coordinates & wheelset coordinates 

      ! convert wheel-marker from wheel-set coords to track-coords

      mq_cp   = marker( 0d0, 0d0, gd%kin%pen )
      mq_trk  = marker_2glob( mq_cp, cp%mref )

      mw_trk  = marker_2glob( my_wheel%m_ws, ws%m_trk )
      mq_w    = marker_2loc( mq_trk, mw_trk )

      ! call marker_print(my_wheel%m_ws, 'm_w(ws)',  1)
      ! call marker_print(mq_trk       , 'm_cp(trk)', 1)
      ! call marker_print(mq_w        , 'm_cp(w)',  1)

      mq_ws   = marker_2loc(mq_trk, ws%m_trk)

      ! get local rolling radius: z-position of the contact point on the wheelset

      w_loc   = mq_ws%z()

      ! define velocity vectors 'tvel' and 'rvel' for flexibilities (vx,vy,vz, vroll,vyaw,vpitch:
      ! velocity of wheel profile origin O w.r.t. wheelset center

      whl_tvel_ws = vec(     my_wheel%vx,   sgn*my_wheel%vy,         my_wheel%vz )
      whl_rvel_ws = vec( sgn*my_wheel%vroll,    my_wheel%vpitch, sgn*my_wheel%vyaw )

      ! call write_log(' velocity of O_whl wrt O_ws')
      ! call vec_print(whl_tvel_ws, 'tvel_whl(ws)', 2)
      ! call vec_print(whl_rvel_ws, 'rvel_whl(ws)', 2)

      ! compute velocity vectors 'tvel' and 'rvel' of the contact reference point Q in wheelset coordinates
      ! v_Q(ws) = v_w(ws) + omega_w(ws) x ( R_w(ws) * m_Q(w) ) + R_w(ws) * v_Q(w)

      associate( whl_rot_ws => my_wheel%m_ws%rot )
      vq_tvel_ws = vec_veloc2glob(whl_tvel_ws, whl_rvel_ws, whl_rot_ws, mq_w, vec_zero())
      vq_rvel_ws = whl_rvel_ws
      end associate

      ! call write_log(' position of Q wrt O_w')
      ! call marker_print(mq_w,         'Q(w)',   2)

      ! call write_log(' velocity of Q wrt O_ws')
      ! call vec_print(vq_tvel_ws, 'tvel_cp(ws)', 2)

      ! define wheelset velocity 'tvel' and 'rvel' in terms of track coordinates
      ! angular velocity: using \omega_ws \approx \dot{\nu}_ws, \nu = [\phi,\theta,\psi] (Euler angles)

      ws_tvel_trk  = vec(     ws%vs,    sgn*ws%vy,         ws%vz   )
      ws_rvel_ws   = vec( sgn*ws%vroll,     ws%vpitch, sgn*ws%vyaw )
      ws_rvel_trk  = ws%m_trk%rot * ws_rvel_ws

      ! compute velocity vectors 'tvel' and 'rvel' of the contact reference point Q in track coordinates
      ! v_Q(trk) = v_ws(trk) + omega_ws(trk) x ( R_ws(trk) * m_Q(ws) ) + R_ws(trk) * v_Q(ws)

      if (idebug.ge.4) then
         call rotmat_print( ws%m_trk%rot,                                  '     R_ws(trk)          =', 2)
         vq_tvel_trk = ws_tvel_trk + (ws_rvel_trk .cross. (ws%m_trk%rot * mq_ws%o)) +    &
                                                             (ws%m_trk%rot * vq_tvel_ws)
         call vec_print( ws_rvel_trk,                                      'A2a: omg_ws(trk)        =', 2, 8)
         call vec_print( mq_ws%o,                                          'A2b: mQ(ws)             =', 2, 8)
         call vec_print( (ws%m_trk%rot * mq_ws%o),                         'A2c: (R*mQ)             =', 2, 8)
         call vec_print(ws_tvel_trk,                                       'A1 : v_ws(trk)          =', 2, 8)
         call vec_print( (ws_rvel_trk .cross. (ws%m_trk%rot * mq_ws%o)),   'A2 : omg_ws(trk) x (R*mQ)', 2, 8)
         call vec_print( (ws%m_trk%rot * vq_tvel_ws),                      'A3 : R*vQ(ws)           =', 2, 8)
         call vec_print(vq_tvel_trk,                                       'A  : tvel_cp(trk)       =', 2, 8)
      endif

      vq_tvel_trk = vec_veloc2glob(ws_tvel_trk, ws_rvel_trk, ws%m_trk%rot, mq_ws, vq_tvel_ws)

      if (idebug.ge.4) then
         call vec_print(vq_tvel_trk,                                       'B  : tvel_cp(trk)       =', 2, 8)
      endif

      vq_rvel_trk = ws_rvel_trk + ws%m_trk%rot * vq_rvel_ws

      ! call write_log(' velocity of Q_cp wrt O_trk')
      ! call vec_print(vq_tvel_trk, 'tvel_cp(trk)', 2)
      ! call vec_print(vq_rvel_trk, 'rvel_cp(trk)', 2)

      ! compute velocity vectors 'tvel' and 'rvel' of the contact reference point Q in contact coordinates

      vq_tvel_cp = cp%mref%rot .transp. vq_tvel_trk
      vq_rvel_cp = cp%mref%rot .transp. vq_rvel_trk

      ! call write_log(' velocity of Q_cp wrt O_cp')
      ! call vec_print(vq_tvel_cp, 'tvel_cp(cp)', 2)
      ! call vec_print(vq_rvel_cp, 'rvel_cp(cp)', 2)

      ! compute the effect of pitch rotation separately
      ! v_Qp(trk) = ( R_ws(trk) * omega_ws(ws) ) x ( R_ws(trk) * m_Q(ws) )

      pitch_tvel_trk = (ws%m_trk%rot * vec( 0d0, ws%vpitch, 0d0 )) .cross. (ws%m_trk%rot * mq_ws%o)

      if (idebug.ge.4) then
         write(bufout,'(4(a,f13.6))') ' vs=',ws%vs,', tvel_pitch=',pitch_tvel_trk%x(),              &
                ', tvel_trnsl=', vq_tvel_trk%x()-pitch_tvel_trk%x(),', tvel_tot=', vq_tvel_trk%x()
         call write_log(1, bufout)
      endif

      !------------------------------------------------------------------------------------------------------
      ! compute overall velocity & creepages

      ! set the velocity == average of forward VS + angular component OMEGA*(R+DR)
      ! set rolling direction on basis of sign(VS) or -sign(OMEGA)
      ! chi = 0 or 180: coordinate system moves along rail, even if the wheelset moves laterally

      if (ic%tang.eq.1) then
         ! in shifts, set dummy V = 1, dt = 1, dq = V*dt
         gd%kin%veloc   = 1d0
         gd%kin%chi     = 0d0
         gd%kin%dq      = 1d0
      elseif (.not.ic%is_roller()) then
         ! gd%kin%veloc =   ws%vs
         ! gd%kin%veloc =  (ws%vs - ws%vpitch * ws%nom_radius) / 2d0
         vs_est = vq_tvel_trk%x() - pitch_tvel_trk%x()
         ! gd%kin%veloc   =   abs(vs_est)
         gd%kin%veloc   =  (abs(vs_est) + abs(pitch_tvel_trk%x())) / 2d0
         if (abs(vs_est).ge.abs(pitch_tvel_trk%x()) .and. vs_est.ge.0d0) then
            gd%kin%chi  = 0d0
         elseif (abs(pitch_tvel_trk%x()).gt.abs(vs_est) .and. pitch_tvel_trk%x().le.0d0) then
            gd%kin%chi  = 0d0
         else
            gd%kin%chi  = pi
         endif
         if (.true.) then
            gd%kin%dq      =  dqrel * gd%cgrid_inp%dx
         else
            gd%kin%dq      =  gd%kin%veloc * gd%kin%dt
            write(bufout,'(6(a,g12.4))') ' dq1=',dqrel*gd%cgrid_inp%dx,' =',dqrel,' *',gd%cgrid_inp%dx, &
                                         ', dq2=',gd%kin%dq,' =',gd%kin%veloc,' *',gd%kin%dt
            call write_log(1, bufout)
         endif
      else
         vs_est = -vp_tvel_cp%x()
         gd%kin%veloc   =  abs(vp_tvel_cp%x() + vq_tvel_cp%x()) / 2d0
         if (abs(vp_tvel_cp%x()).ge.abs(vq_tvel_trk%x()) .and. vp_tvel_cp%x().le.0d0) then
            gd%kin%chi  = 0d0
         elseif (abs(vq_tvel_trk%x()).gt.abs(vp_tvel_cp%x()) .and. vq_tvel_trk%x().le.0d0) then
            gd%kin%chi  = 0d0
         else
            gd%kin%chi  = pi
         endif
         gd%kin%dq      =  dqrel * gd%cgrid_inp%dx
      endif

      ! creepage = { velocity of rail/roller (1) - velocity of wheel (2) } / reference velocity

      if (abs(gd%kin%veloc).lt.1d-6) then
         ! rolling at V=0: creepages +/- 1 for pos/neg velocity difference
         gd%kin%cksi = 1d0 * sign(1d0, vp_tvel_cp%x() - vq_tvel_cp%x())
         gd%kin%ceta = 1d0 * sign(1d0, vp_tvel_cp%y() - vq_tvel_cp%y())
         gd%kin%cphi = 1d0 * sign(1d0, vp_rvel_cp%z() - vq_rvel_cp%z())
      else
         ! rolling: compute creepage; shift: V=dt=1, cksi--cphi = shift distance
         gd%kin%cksi = (vp_tvel_cp%x() - vq_tvel_cp%x()) / gd%kin%veloc
         gd%kin%ceta = (vp_tvel_cp%y() - vq_tvel_cp%y()) / gd%kin%veloc
         gd%kin%cphi = (vp_rvel_cp%z() - vq_rvel_cp%z()) / gd%kin%veloc
      endif

      !write(bufout,'(2(a,f9.2),a,f5.1,a,f9.5)') 'vs=',vs_est,', vrot=',pitch_tvel_trk%x(),              &
      !          ': chi=',gd%kin%chi*180d0/pi, ', cksi=', gd%kin%cksi
      !call write_log(1, bufout)

      if (idebug.ge.2 .and. idebug.lt.3) then
         call write_log(' velocity of contact point on wheel w.r.t. O_trk')
         call vec_print(vq_tvel_trk,  'vq_tvel(trk)', 2, 6)
         if (ic%is_roller()) then
            call write_log(' velocity of contact point on roller w.r.t. O_trk')
            call vec_print(vp_tvel_trk,  'vp_tvel(trk)', 2, 6)
         endif
      elseif (idebug.ge.3) then
         call write_log(' wheel-set position w.r.t. O_trk')
         call marker_print(ws%m_trk,    'm_ws(trk)', 3)
         call write_log(' wheel-set ang.velocity w.r.t. O_ws')
         call vec_print(ws_rvel_ws,  ' ws_rvel(ws)', 2, 6)
         call write_log(' wheel-set velocity w.r.t. O_trk')
         call vec_print(ws_tvel_trk, 'ws_tvel(trk)', 2, 6)
         call vec_print(ws_rvel_trk, 'ws_rvel(trk)', 2, 6)
         call write_log(' contact reference position w.r.t. O_ws')
         call marker_print(mq_ws, 'm_Q(ws)', 2)
         write(bufout,'(a,f11.6)') ' actual rolling radius on wheel r(cp) =', w_loc
         call write_log(1, bufout)
         ! call write_log(' contact reference position w.r.t. O_trk')
         ! call marker_print(cp%mref, 'm_cp(trk)', 3)
         call write_log(' velocity of contact point on wheel w.r.t. O_ws (=wheelset flexibility)')
         call vec_print(vq_tvel_ws,   ' vq_tvel(ws)', 2, 6)
         call write_log(' velocity of contact point on wheel w.r.t. O_trk (eq.28: v_ws(tr) + ' //       &
                                                                                'cross(omg, R.xQ) + R.vQ)')
         call vec_print(vq_tvel_trk,  'vq_tvel(trk)', 2, 6)
         call write_log(' velocity of contact point on wheel w.r.t. O_cp (=R^T vQ(trk)))')
         call vec_print(vq_tvel_cp,   ' vq_tvel(cp)', 2, 6)
         call write_log(' rotation of contact point on wheel w.r.t. O_cp (=R^T omgQ(trk))')
         call vec_print(vq_rvel_cp,   ' vq_rvel(cp)', 2, 6)

         if (ic%is_roller()) then
            call write_log(' contact reference position w.r.t. O_rol')
            call marker_print(mref_rol, 'm_cp(rol)', 2)
            write(bufout,'(a,f11.6)') ' actual rolling radius on roller r(cp) =', r_loc
            call write_log(1, bufout)
            call write_log(' velocity of contact point on roller w.r.t. O_rol')
            call vec_print(vp_tvel_rol,  'vp_tvel(rol)', 2, 6)
            call write_log(' velocity of contact point on roller w.r.t. O_trk')
            call vec_print(vp_tvel_trk,  'vp_tvel(trk)', 2, 6)
            call write_log(' velocity of contact point on roller w.r.t. O_cp')
            call vec_print(vp_tvel_cp, ' vp_tvel(cp)', 2, 6)
            call write_log(' rotation of contact point on roller w.r.t. O_cp')
            call vec_print(vp_rvel_cp,  '  vp_rvel_cp', 2, 6)
         endif

         write(bufout, '(a,f10.1,3(a,g14.6))') ' veloc=',gd%kin%veloc,', cksi=',gd%kin%cksi,            &
                ', ceta=',gd%kin%ceta, ',  cphi=',gd%kin%cphi
         call write_log(1, bufout)
      endif

      end associate
      if (idebug.ge.4) call write_log(' --- end subroutine wr_creep_cref ---')

   end subroutine wr_creep_cref

!------------------------------------------------------------------------------------------------------------

   subroutine wr_rigslip_wrsurf(ic, ws, trk, dqrel, icp, cp, gd, idebug)
!--purpose: compute the rigid slip at each grid point in the contact reference plane
      implicit none
!--subroutine arguments:
      integer               :: idebug, icp
      type(t_ic)            :: ic
      type(t_wheelset)      :: ws
      type(t_trackdata)     :: trk
      type(t_cpatch)        :: cp
      type(t_probdata)      :: gd
      real(kind=8)          :: dqrel
!--local variables:
      logical, parameter    :: project_on_cref = .false.
      integer               :: ii, iimin, iidbg, is_right
      real(kind=8)          :: sgn, dst2, dstmin, zmin, zmax
      type(t_marker)        :: m_rol, mrol_cp, mr_rol, mr_cp, mw_arm, mw_trk, mw_cp, mws_cp
      type(t_vec)           :: rail_tvel_rol, rail_rvel_rol, rail_tvel_trk, rail_rvel_trk,              &
                               rail_tvel_cp, rail_rvel_cp,                                              &
                               rol_tvel_trk, rol_rvel_rol, rol_tvel_cp, rol_rvel_cp,                    &
                               ws_tvel_trk, ws_rvel_ws, ws_tvel_cp, ws_rvel_cp,                         &
                               whl_tvel_ws, whl_rvel_ws, whl_tvel_cp, whl_rvel_cp, nref
      type(t_gridfnc3)      :: gf_tvel_rail, gf_tvel_whl
      character(len=5)        :: nam_side 

      if (idebug.ge.3) call write_log(' --- Start subroutine wr_rigslip_wrsurf ---')
      if (idebug.ge.2) then
         if (project_on_cref) then
            call write_log(' computing rigid slip with projection on contact reference plane')
         elseif (ic%discns1_eff.le.3) then
            call write_log(' planar contact, with rigid slip on actual surfaces in planar (s,n)')
         elseif (ic%is_conformal()) then
            call write_log(' conformal contact, with rigid slip on actual surfaces, curved (s,n)')
         endif
      endif

      associate(my_rail   => trk%rai, my_wheel  => ws%whl)
      if (ic%is_left_side()) then
         nam_side  = 'left'
         is_right  = 0
      else
         nam_side  = 'right'
         is_right  = 1
      endif

      ! set the sign to -1 for left and +1 for right rail/wheel combination

      sgn = 2 * is_right - 1

      iidbg = cp%whl_srfc%nx/2 + (cp%whl_srfc%ny/2-1) * cp%whl_srfc%nx

      !------------------------------------------------------------------------------------------------------
      ! compute translational velocities of points on the rail/roller surface

      if (.not.ic%is_roller()) then

         ! convert rail profile origin to contact coordinates

         ! call marker_print(my_rail%m_trk, 'mr_trk',  2)
         ! call marker_print(cp%mpot,       'mcp_trk', 2)
         mr_cp = marker_2loc(my_rail%m_trk, cp%mpot)
         ! call marker_print(mr_cp,         'mr_cp',   2)

         ! define velocity vectors 'tvel' and 'rvel' for rail profile (origin O) w.r.t track coordinates

         rail_tvel_trk = vec( 0d0, sgn*my_rail%vy, my_rail%vz )
         rail_rvel_trk = vec( sgn*my_rail%vroll, 0d0, 0d0 )

         ! convert velocity vectors 'tvel' and 'rvel' to contact coordinates

         rail_rvel_cp = cp%mpot%rot .transp. rail_rvel_trk
         rail_tvel_cp = cp%mpot%rot .transp. rail_tvel_trk + rail_rvel_cp .cross. mr_cp%o 
         ! call vec_print(rail_tvel_cp, 'vr_cp', 2)
         ! call vec_print(rail_rvel_cp, 'wr_cp', 2)

         ! compute velocities of points P in the contact grid in contact coordinates

         if (project_on_cref) then
            ! simplfied method: measure velocities for points in the contact reference plane
            cp%rail_srfc%z(1:cp%rail_srfc%ntot) = 0d0
         endif

         call gf3_new(gf_tvel_rail, 'gf_tvel_rail', cp%rail_srfc, lzero=.true.)
         call gf3_veloc2glob(mr_cp, rail_tvel_cp, rail_rvel_cp, cp%rail_srfc, gf_tvel_rail)

         if (idebug.ge.2) then
            write(bufout,'(a,i6,a,2g12.4)') ' ii=',iidbg,': rail_tvel(x,y)=', gf_tvel_rail%vx(iidbg),    &
                gf_tvel_rail%vy(iidbg)
            call write_log(1, bufout)
         endif
      else

         !---------------------------------------------------------------------------------------------------
         ! compute translational velocities of points P_I on the roller surface in contact coordinates

         ! define marker for roller in track coords, convert to contact coordinates

         m_rol%o   = vec( 0d0, 0d0, trk%nom_radius )
         m_rol%rot = rotmat_identity()
         mr_rol    = marker_2loc( my_rail%m_trk, m_rol )
         mr_cp     = marker_2loc( my_rail%m_trk, cp%mpot )
         mrol_cp   = marker_2loc( m_rol, cp%mpot )

         if (idebug.ge.2) then
            call marker_print(m_rol,   'm_rol',   2)
            call marker_print(my_rail%m_trk, 'm_r', 2)
            call marker_print(cp%mpot, 'm_cp',    2)
            call marker_print(mr_rol,  'mr_rol',  2)
            call marker_print(mr_cp,   'mr_cp',   2)
            call marker_print(mrol_cp, 'mrol_cp', 2)
         endif

         ! define velocity vectors 'tvel' and 'rvel' of roller in track/roller coordinates
         ! angular velocity: using \omega_rol \approx \dot{\nu}_rol, \nu = [\phi,\theta,\psi] (Euler angles)

         rol_tvel_trk  = vec_zero()
         rol_rvel_rol  = vec( 0d0, trk%vpitch_rol, 0d0 )

         ! convert wheelset velocity to contact coordinates

         rol_tvel_cp   = cp%mpot%rot .transp. rol_tvel_trk
         rol_rvel_cp   = mrol_cp%rot * rol_rvel_rol

         ! define velocity vectors 'tvel' and 'rvel' for roller profile origin O w.r.t. roller coordinates
         ! (roller flexibility: vy, vz, vroll)

         rail_tvel_rol = vec( 0d0,  sgn*my_rail%vy, my_rail%vz )
         rail_rvel_rol = vec( sgn*my_rail%vroll,  0d0,  0d0 )

         ! compute velocity of the roller profile origin in contact coordinates

         rail_tvel_cp  = vec_veloc2glob(rol_tvel_cp, rol_rvel_cp, mrol_cp%rot, mr_rol, rail_tvel_rol)
         rail_rvel_cp  = rol_rvel_cp + mrol_cp%rot * rail_rvel_rol

         if (idebug.ge.2) then
            call vec_print(rol_tvel_cp,  'rol_tvel_cp',  2)
            call vec_print(rol_rvel_cp,  'rol_rvel_cp',  2)
            call vec_print(rail_tvel_cp, 'rail_tvel_cp', 2)
            call vec_print(rail_rvel_cp, 'rail_rvel_cp', 2)
         endif

         ! compute velocities of points P in the roller surface in contact coordinates

         if (idebug.ge.2) then
            call grid_get_zrange(cp%rail_srfc, zmin, zmax)
            write(bufout,'(a,i6,a,2g12.4,a,g12.4)') ' ii=',iidbg,': rail%z(1,dbg)= ', cp%rail_srfc%z(1), &
                cp%rail_srfc%z(iidbg), ', min=', zmin
            call write_log(1, bufout)
         endif

         if (project_on_cref) then
            ! simplfied method: measure velocities for points in the contact reference plane
            cp%rail_srfc%z(1:cp%rail_srfc%ntot) = 0d0
         endif

         call gf3_new(gf_tvel_rail, 'gf_tvel_rail', cp%rail_srfc, lzero=.true.)
         call gf3_veloc2glob(mr_cp, rail_tvel_cp, rail_rvel_cp, cp%rail_srfc, gf_tvel_rail)

         if (idebug.ge.2) then
            write(bufout,'(a,i6,a,2g12.4)') ' ii=',iidbg,': rail_tvel(x,y)= ', gf_tvel_rail%vx(iidbg),  &
                gf_tvel_rail%vy(iidbg)
            call write_log(1, bufout)
         endif

      endif

      !------------------------------------------------------------------------------------------------------
      ! compute translational velocities of points Q_I on the wheel surface in contact coordinates

      ! convert wheel-marker from wheel-set coords to track-coords and contact coordinates

      mw_trk  = marker_2glob( my_wheel%m_ws, ws%m_trk )
      mw_cp   = marker_2loc( mw_trk, cp%mpot )
      mws_cp  = marker_2loc( ws%m_trk, cp%mpot )

      if (idebug.ge.2) then
         call marker_print(ws%m_trk, 'm_ws',   2)
         call marker_print(my_wheel%m_ws, 'mw_ws', 2)
         call marker_print(mw_trk,   'mw_trk', 2)
         call marker_print(mw_cp,    'mw_cp',  2)
         call marker_print(mws_cp,   'mws_cp', 2)
      endif

      ! define velocity vectors 'tvel' and 'rvel' of wheelset in track/wheelset coordinates
      ! angular velocity: using \omega_ws \approx \dot{\nu}_ws, \nu = [\phi,\theta,\psi] (Euler angles)

      ws_tvel_trk  = vec(     ws%vs,    sgn*ws%vy,         ws%vz   )
      ws_rvel_ws   = vec( sgn*ws%vroll,     ws%vpitch, sgn*ws%vyaw )

      ! convert wheelset velocity to contact coordinates

      ws_tvel_cp   = cp%mpot%rot .transp. ws_tvel_trk
      ws_rvel_cp   = mws_cp%rot * ws_rvel_ws

      ! define velocity vectors 'tvel' and 'rvel' for flexibilities (vx,vy,vz, vroll,vyaw,vpitch):
      ! velocity of wheel profile origin O w.r.t. wheelset coordinates

      whl_tvel_ws  = vec(     my_wheel%vx,   sgn*my_wheel%vy,         my_wheel%vz )
      whl_rvel_ws  = vec( sgn*my_wheel%vroll,    my_wheel%vpitch, sgn*my_wheel%vyaw )

      ! convert velocity of wheel profile origin from wheelset to contact coordinates
      ! using cp as global system, ws as local system, whl_tvel_ws as flexibility

      if (idebug.ge.2) then
         call vec_print(ws_rvel_cp,  'ws_rvel_cp(1)', 2)
         call vec_print(whl_rvel_ws, 'whl_rvel_ws', 2)
      endif

      ! exclude flexibility dx in velocity whl_tvel_cp
      ! TODO: exclude flexibility dz in case of axle bending

      mw_arm = my_wheel%m_ws
      call marker_shift(mw_arm, -my_wheel%dx, 0d0, 0d0)

      if (idebug.ge.2) call marker_print(mw_arm, 'mw_ws(*)', 2)

      whl_tvel_cp = vec_veloc2glob(ws_tvel_cp, ws_rvel_cp, mws_cp%rot, mw_arm, whl_tvel_ws)
      whl_rvel_cp = ws_rvel_cp + mws_cp%rot * whl_rvel_ws

      if (idebug.ge.2) then
         call vec_print(ws_tvel_trk, 'ws_tvel_trk', 2)
         call vec_print(ws_rvel_ws,  'ws_rvel_ws',  2)
         call vec_print(ws_tvel_cp,  'ws_tvel_cp',  2)
         call vec_print(ws_rvel_cp,  'ws_rvel_cp',  2)
         call vec_print(whl_tvel_ws, 'whl_tvel_ws', 2)
         call vec_print(whl_rvel_ws, 'whl_rvel_ws', 2)
         call vec_print(whl_tvel_cp, 'whl_tvel_cp', 2)
         call vec_print(whl_rvel_cp, 'whl_rvel_cp', 2)
      endif

      ! compute velocities of points Q in the wheel surface in contact coordinates

      if (idebug.ge.2) then
         call grid_get_zrange(cp%whl_srfc, zmin, zmax)
         write(bufout,'(a,i6,a,2g12.4,a,g12.4)') ' ii=',iidbg,': whl%z(1,dbg)= ', cp%whl_srfc%z(1),     &
                cp%whl_srfc%z(iidbg), ', max=', zmax
         call write_log(1, bufout)
      endif

      if (project_on_cref) then
         ! simplfied method: measure velocities for points in the contact reference plane
         cp%whl_srfc%z(1:cp%whl_srfc%ntot) = gd%kin%pen
      endif

      call gf3_new(gf_tvel_whl, 'gf_tvel_whl', cp%whl_srfc, lzero=.true.)
      call gf3_veloc2glob(mw_cp, whl_tvel_cp, whl_rvel_cp, cp%whl_srfc, gf_tvel_whl)
      ! call gf3_print(gf_tvel_whl, 'gf_tvel_whl', ikZDIR, 5)

      ! conformal method: rotate velocities to local surface (s,n) coordinates

      if (ic%is_conformal()) then

         ! D=4: conformal rigid slip on curved potential contact

         if (idebug.ge.2) call write_log(' call gf3_cart2curv for conformal...')
         if (idebug.ge.5) call gf3_print(cp%curv_nrm, 'curv_nrm', ikALL, 5)

         nref = -cp%mref%kvec()
         call gf3_cart2curv(gf_tvel_rail, nref, cp%curv_nrm, 1, idebug)
         call gf3_cart2curv(gf_tvel_whl, nref, cp%curv_nrm, 1, idebug)

      endif

      if (idebug.ge.2) then
         write(bufout,'(a,i6,a,2g12.4)') ' ii=',iidbg,': whl_tvel= ', gf_tvel_whl%vx(iidbg),     &
                gf_tvel_whl%vy(iidbg)
         call write_log(1, bufout)
      endif

      !------------------------------------------------------------------------------------------------------
      ! compute overall velocity & rigid slip

      ! find central grid point

      iimin = 0
      dstmin = 1d10

      do ii = 1, cp%whl_srfc%ntot
         dst2 = cp%whl_srfc%x(ii)**2 + cp%whl_srfc%y(ii)**2
         if (dst2.lt.dstmin) then
            iimin = ii
            dstmin = dst2
         endif
      enddo

      ! set the velocity, rolling direction and rolling step size

      gd%kin%chi    =  0d0
      if (ic%tang.eq.1) then
         ! in shifts, set dummy V = 1, dt = 1, dq = V*dt
         gd%kin%veloc   = 1d0
      elseif (.not.ic%is_roller()) then
         ! gd%kin%veloc =  (ws%vs - ws%vpitch * ws%nom_radius) / 2d0
         gd%kin%veloc  =  ws%vs
      else
         gd%kin%veloc   =  - (gf_tvel_rail%vx(iimin) + gf_tvel_whl%vx(iimin)) / 2d0
      endif

      if (gd%kin%veloc.lt.0d0) then
         gd%kin%chi    =  pi
         gd%kin%veloc  =  abs(gd%kin%veloc)
      endif

      ! coordinate system moves along rail, even if the wheelset moves laterally

      gd%kin%dq     =  dqrel * gd%cgrid_inp%dx

      ! creepages are zero

      gd%kin%cksi   = 0d0
      gd%kin%ceta   = 0d0
      gd%kin%cphi   = 0d0
      gd%ic%rztang  = 1

      if (idebug.ge.4) then
         call write_log(' velocity of contact point on wheel w.r.t. O_trk')
         write(bufout,'(3(a,f12.5),a)') ' vq_tvel(trk) = (', gf_tvel_whl%vx(iimin), ',',                &
                gf_tvel_whl%vy(iimin), ',',gf_tvel_whl%vn(iimin), ')^T'
         call write_log(1, bufout)
         if (ic%is_roller()) then
            call write_log(' velocity of contact point on roller w.r.t. O_trk')
            write(bufout,'(3(a,f12.5),a)') ' vp_tvel(trk) = (', gf_tvel_rail%vx(iimin), ',',            &
                gf_tvel_rail%vy(iimin), ',',gf_tvel_rail%vn(iimin), ')^T'
            call write_log(1, bufout)
         endif
      endif

      ! full method: project velocities vP and vQ on inclined midplane vmid 

      if (.not.project_on_cref) then
         call project_midplane(ic%is_roller(), gf_tvel_rail, gf_tvel_whl, gd, iimin, idebug)

         if (idebug.ge.2) then
            if (.not.ic%is_roller()) then
               call write_log(' velocity of contact point on wheel w.r.t. moving O_trk')
            else
               call write_log(' velocity of contact point on wheel w.r.t. O_trk')
            endif
            write(bufout,'(3(a,f12.5),a)') ' vq_tvel(trk) = (', gf_tvel_whl%vx(iimin), ',',             &
                   gf_tvel_whl%vy(iimin), ',',gf_tvel_whl%vn(iimin), ')^T'
            call write_log(1, bufout)
            if (.not.ic%is_roller()) then
               call write_log(' velocity of contact point on rail w.r.t. moving O_trk')
            else
               call write_log(' velocity of contact point on roller w.r.t. O_trk')
            endif
            write(bufout,'(3(a,f12.5),a)') ' vp_tvel(trk) = (', gf_tvel_rail%vx(iimin), ',',            &
                   gf_tvel_rail%vy(iimin), ',',gf_tvel_rail%vn(iimin), ')^T'
            call write_log(1, bufout)
         endif
      endif

      ! store rigid slip in exrhs: { velocity of rail/roller (1) - velocity of wheel (2) } / ref. velocity

      call gf3_new(gd%geom%exrhs, 'geom%exrhs', gf_tvel_rail%grid)
      call gf3_copy(AllElm, gf_tvel_rail, gd%geom%exrhs, ikTANG)
      call gf3_axpy(AllElm, -1d0, gf_tvel_whl, gd%geom%exrhs, ikTANG)

      if (abs(gd%kin%veloc).gt.1d-6) then
         call gf3_scal(AllElm, 1d0/gd%kin%veloc, gd%geom%exrhs, ikTANG)
      endif

      ! move planar part from exrhs to creepages cksi, ceta, cphi

      if (idebug.ge.4) then
         write(bufout,'(a,i0,a)') ' patch icp= ',icp,': extracting planar creepages from exrhs'
         call write_log(1, bufout)
      endif
      call extract_creepages( gd%ic, gd%cgrid_inp, gd%kin, gd%geom%prmudf, gd%geom%exrhs )

      if (idebug.ge.4) then
         call gf3_print(gd%geom%exrhs, 'exrhs', ikTANG, 4)
      endif

      call gf3_destroy(gf_tvel_rail)
      call gf3_destroy(gf_tvel_whl)

      end associate
      if (idebug.ge.4) call write_log(' --- end subroutine wr_rigslip_wrsurf ---')

   end subroutine wr_rigslip_wrsurf

!------------------------------------------------------------------------------------------------------------

   subroutine extract_creepages( ic, cgrid, kin, prmudf, exrhs )
!--purpose: define equivalent creepages, extracting contributions from exrhs
      implicit none
!--subroutine arguments:
      type(t_ic)                 :: ic
      type(t_grid)               :: cgrid
      type(t_kincns)             :: kin
      type(t_gridfnc3)           :: exrhs
      real(kind=8), dimension(:) :: prmudf
!--local variables:
      logical               :: is_roll
      integer               :: ii, ncon
      real(kind=8)          :: gap, dst_avg, wx_avg, wy_avg, wz_avg, x_avg, y_avg, xofs, yofs,          &
                               xofs_dq, yofs_dq, denom

      associate( mx => cgrid%nx, my => cgrid%ny, npot => cgrid%ntot, x => cgrid%x, y => cgrid%y,        &
                 h  => prmudf )

      is_roll = ic%tang.eq.2 .or. ic%tang.eq.3

      ! offset 1/6 dq used in rolling problems:

      if (is_roll) then
         xofs_dq = cos(kin%chi) * kin%dq * kin%facphi
         yofs_dq = sin(kin%chi) * kin%dq * kin%facphi
      else
         xofs_dq = 0d0
         yofs_dq = 0d0
      endif

      ! solve [ cksi, ceta, cphi ] from equations [ wx = cksi - yofs * cphi, wy = ceta + xofs * cphi ]
      ! using normal equations, At * A * c = At * w
      ! At * A: a11 = a22 = sum(1), a13 = sum(-yofs), a23 = sum(xofs), a33 = sum(xofs**2 + yofs**2)
      ! At * w: d1 = sum(wx), d2 = sum(wy), d3 = sum(xofs*wy - yofs*wx)

      ! compute average position and average rigid slip over interpenetration area

      ncon    = 0
      x_avg   = 0d0
      y_avg   = 0d0
      dst_avg = 0d0
      wx_avg  = 0d0
      wy_avg  = 0d0
      wz_avg  = 0d0

      do ii = 1, npot
         gap = h(ii) - kin%pen
         if (gap.lt.0d0) then
            xofs    = x(ii) + xofs_dq - kin%spinxo
            yofs    = y(ii) + yofs_dq - kin%spinyo

            ncon    = ncon    + 1
            x_avg   = x_avg   + xofs
            y_avg   = y_avg   + yofs
            dst_avg = dst_avg + xofs**2 + yofs**2
            wx_avg  = wx_avg  + exrhs%vx(ii)
            wy_avg  = wy_avg  + exrhs%vy(ii)
            wz_avg  = wz_avg  + xofs * exrhs%vy(ii) - yofs * exrhs%vx(ii)
         endif
      enddo

      if (ncon.ge.1) then
         x_avg   = x_avg  / ncon
         y_avg   = y_avg  / ncon
         dst_avg = sqrt(dst_avg / ncon)
         wx_avg  = wx_avg / ncon
         wy_avg  = wy_avg / ncon
         wz_avg  = wz_avg / ncon
         denom   = dst_avg**2 - y_avg**2 - x_avg**2
      else
         denom   = 1d0
      endif

      if (ncon.le.0 .or. denom.lt.1d-12) then
         x_avg   = 0d0
         y_avg   = 0d0
         dst_avg = 1d0
         wx_avg  = 0d0
         wy_avg  = 0d0
         wz_avg  = 0d0
         denom   = 1d0
      endif

      ! compute least squares approximation (cksi, ceta, cphi)

      kin%cphi = (wz_avg + y_avg * wx_avg - x_avg * wy_avg) / denom
      kin%cksi =  wx_avg + y_avg * kin%cphi
      kin%ceta =  wy_avg - x_avg * kin%cphi

      ! extract (cksi, ceta, cphi) from exrhs

      do ii = 1, npot
         xofs   = x(ii) + xofs_dq - kin%spinxo
         yofs   = y(ii) + yofs_dq - kin%spinyo
         exrhs%vx(ii) = exrhs%vx(ii) - kin%cksi + yofs * kin%cphi
         exrhs%vy(ii) = exrhs%vy(ii) - kin%ceta - xofs * kin%cphi
      enddo

      if (ic%x_nmdbg.ge.1) call gf3_check_nan(exrhs, 'extract_creepage: exrhs', AllElm, ikTANG, 2)
      end associate

   end subroutine extract_creepages

!------------------------------------------------------------------------------------------------------------

   subroutine project_midplane(is_roller, gf_tvel_rail, gf_tvel_whl, gd, iimin, idebug)
!--purpose: project particle velocities vP and vQ onto the common midplane vmid
      implicit none
!--subroutine arguments:
      logical               :: is_roller
      integer               :: iimin, idebug
      type(t_gridfnc3)      :: gf_tvel_rail, gf_tvel_whl
      type(t_probdata)      :: gd
!--local variables:
      integer               :: mx, my, npot, ii, ixdbg, iydbg, iidbg
      real(kind=8)          :: fac_rail, sgn, vabsi, v1x, v1n, v2x, v2n
      type(t_gridfnc3)      :: vmid
      character(len=4)      :: namrail

      mx   = gf_tvel_rail%grid%nx
      my   = gf_tvel_rail%grid%ny
      npot = mx*my

      if (.false.) then
         iydbg = (iimin-1)/mx + 1
         ixdbg = iimin - (iydbg-1)*mx
      else
         ixdbg = -17
         iydbg =   5
      endif
      iidbg = ixdbg + (iydbg-1)*mx
      if (ixdbg.le.0 .or. iydbg.le.0 .or. iidbg.gt.npot) iidbg = -1

      ! create grid-function for mid-plane inclination

      call gf3_new(vmid, 'vmid', gf_tvel_rail%grid, lzero=.true.)

      ! set sign for forward direction

      if (abs(gd%kin%chi-pi).gt.0.01d0) then         ! chi =   0
         sgn =  1d0
      else                                           ! chi = 180 deg
         sgn = -1d0
      endif

      ! wheel-on-rail: subtract coordinate reference velocity V from x-components of vP and vQ

      if (.not.is_roller) then
         call gf3_set(AllElm, -sgn*gd%kin%veloc, vmid, ikXDIR)
         call gf3_axpy(AllElm, 1d0, vmid, gf_tvel_rail, ikXDIR)
         call gf3_axpy(AllElm, 1d0, vmid, gf_tvel_whl,  ikXDIR)
         call gf3_set(AllElm, 0d0, vmid, ikXDIR)
      endif

      ! compute vmid == tangent vectors along inclined contact reference
      ! using -veloc, chi=0 or 180, tangent pointing in negative/positive x-direction

      associate(gg => gd%mater%gg, poiss => gd%mater%poiss)
      fac_rail = (1d0-poiss(1))*gg(2) / ( (1d0-poiss(1))*gg(2) + (1d0-poiss(2))*gg(1) )
      end associate

      call gf3_axpy(AllElm,      fac_rail , gf_tvel_rail, vmid, ikALL)
      call gf3_axpy(AllElm, (1d0-fac_rail), gf_tvel_whl , vmid, ikALL)

      ! decompose vP and vQ into components along mid-plane & perpendicular to mid-plane

      do ii = 1, npot
         vabsi = -sgn / sqrt( max(1d-16, vmid%vn(ii)**2 + vmid%vx(ii)**2 ) )

         if (idebug.ge.-1 .and. ii.eq.iidbg) then
            write(bufout,'(a,i6,3(a,f11.4),a)') ' ii=',iidbg,': x=', gf_tvel_rail%grid%x(iidbg),         &
                ', y=', gf_tvel_rail%grid%y(iidbg)
            call write_log(1, bufout)
            namrail = 'rail'
            if (is_roller) namrail = 'rol'
            write(bufout,'(a,i6,2a,3(a,f11.4),a)') ' ii=',iidbg,': tvel_',trim(namrail),'=[',           &
                   gf_tvel_rail%vx(iidbg), ',', gf_tvel_rail%vy(iidbg), ',', gf_tvel_rail%vn(iidbg),']'
            call write_log(1, bufout)
            write(bufout,'(a,i6,3(a,f11.4),a)') ' ii=',iidbg,': tvel_whl =[', gf_tvel_whl%vx(iidbg), ',',&
                   gf_tvel_whl%vy(iidbg), ',', gf_tvel_whl%vn(iidbg),']'
            call write_log(1, bufout)
            write(bufout,'(a,i6,4(a,f11.4))')   ' ii=',iidbg,': vmid     =[', vmid%vx(iidbg), ',',       &
                   vmid%vy(iidbg), ',', vmid%vn(iidbg),'], vabsi=', vabsi
            call write_log(1, bufout)
         endif

         v1x  = (gf_tvel_rail%vx(ii) *  vmid%vx(ii) + gf_tvel_rail%vn(ii) * vmid%vn(ii)) * vabsi
         if (idebug.ge.-1 .and. ii.eq.iidbg) then
            write(bufout,'(a,i6,5(a,f11.4))') ' ii=',iidbg,': v1x=(', gf_tvel_rail%vx(ii),' *',          &
                 vmid%vx(ii),' +', gf_tvel_rail%vn(ii),' *', vmid%vn(ii),' ) /', 1d0/vabsi
            call write_log(1, bufout)
         endif
         v1n  = (gf_tvel_rail%vx(ii) * -vmid%vn(ii) + gf_tvel_rail%vn(ii) * vmid%vx(ii)) * vabsi
         gf_tvel_rail%vx(ii) = v1x
         gf_tvel_rail%vn(ii) = v1n

         v2x  = (gf_tvel_whl%vx(ii) *  vmid%vx(ii) + gf_tvel_whl%vn(ii) * vmid%vn(ii)) * vabsi
         v2n  = (gf_tvel_whl%vx(ii) * -vmid%vn(ii) + gf_tvel_whl%vn(ii) * vmid%vx(ii)) * vabsi
         gf_tvel_whl%vx(ii) = v2x
         gf_tvel_whl%vn(ii) = v2n

         if (idebug.ge.-1 .and. ii.eq.iidbg) then
            namrail = 'rail'
            if (is_roller) namrail = 'rol'
            write(bufout,'(a,i6,2(3a,f11.4))') ' ii=',iidbg,': new ',trim(namrail),'_x=',               &
                gf_tvel_rail%vx(iidbg), ', ',trim(namrail),'_n=', gf_tvel_rail%vn(iidbg)
            call write_log(1, bufout)
            write(bufout,'(a,i6,2(a,f11.4))') ' ii=',iidbg,': new  whl_x=', gf_tvel_whl%vx(iidbg),        &
                ', whl_n=', gf_tvel_whl%vn(iidbg)
            call write_log(1, bufout)
         endif
      enddo

      call gf3_destroy(vmid)

   end subroutine project_midplane

!------------------------------------------------------------------------------------------------------------

end module m_wr_rigslip
