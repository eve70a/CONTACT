!------------------------------------------------------------------------------------------------------------
! m_wr_undefdist - compute undeformed distance for one contact patch in a w/r contact problem (module 1)
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_wr_undefdist

use m_wrprof_data
use m_wr_locatecp,   only : find_bounding_box_2d
implicit none
private

public  wr_ud_planar
private curved_sc_to_planar_sp
private trim_rail_mesh
private trim_wheel_profile
private trim_wheel_mesh         ! unused?

contains

!------------------------------------------------------------------------------------------------------------

   subroutine wr_ud_planar (meta, ic, ws, trk, icp, cp, gd, idebug)
!--purpose: compute the undeformed distance on a potential contact area using planar contact approach.
!           TODO: propagate error conditions to calling subroutine
      implicit none
!--subroutine arguments:
      type(t_metadata)      :: meta
      type(t_ic)            :: ic
      type(t_wheelset)      :: ws
      type(t_trackdata)     :: trk
      integer               :: icp, idebug
      type(t_cpatch)        :: cp
      type(t_probdata)      :: gd
!--local variables:
      logical                 :: is_prismatic, use_bicubic
      character(len=5)        :: nam_side
      integer                 :: lud, ios, nx, nxlow, nxhig, nw, nslc, nslow, nshig, ix, iy, ii, ixdbg, &
                                 iydbg, iu, iv, icount, ierror
      real(kind=8)            :: xslc1, hmin, h_ii, xmin, xmax, ymin, ymax, zmin, zmax, smin, smax,     &
                                 swsta, swmid, swend, ds, dx_cp, zsta, zend, z_axle, xtr_sta,           &
                                 vy_ref, vz_ref, vy_iy, vz_iy, delt_iy, vlen_inv
      real(kind=8)            :: rmin, rmax, th_min, th_max, dth, v_min, v_max, dv
      real(kind=8),   dimension(:), allocatable :: u_out, v_out
      type(t_vec)             :: vtmp(4)
      type(t_marker)          :: mwhl_trk, mrai_pot, mwhl_pot, mpot_whl, mpot_rail
      type(t_grid)            :: rail_sn_full, rail_sn, bbr, bbc
      type(t_grid)            :: prw_trim, whl_srfw_full, whl_srfw, bbw
      character(len=256)      :: tmp_fname, tmp_fulnam

      if (idebug.ge.2) call write_log(' --- Start subroutine wr_ud_planar ---')
      call timer_start(itimer_udist)

      associate(my_rail   => trk%rai,      my_wheel => ws%whl,      cgrid => gd%cgrid_inp,              &
                rail_srfc => cp%rail_srfc, whl_srfc => cp%whl_srfc)

      if (ic%is_left_side()) then
         nam_side = 'left'
      else
         nam_side = 'right'
      endif

      is_prismatic  = (.not.ic%is_roller() .and. .not.my_rail%prr%is_varprof())

      if (idebug.ge.2) then
         write(bufout,'(a,i2,2(a,i4),a)') ' patch',icp,': pot.contact has',cgrid%nx,' x',cgrid%ny,      &
                ' elements'
         call write_log(1, bufout)
      endif

      ! 1. determine rail surface (x,s,n)_(pot) on potential contact grid

      if (is_prismatic .or. my_rail%prr%is_varprof()) then

         !  - mpot_rail: pot.con marker in terms of rail coordinates
    
         mpot_rail = marker_2loc(cp%mpot, my_rail%m_trk)
         if (is_prismatic) mpot_rail%o%v(1) = 0d0
         if (idebug.ge.4) call marker_print(mpot_rail, 'm_pot(r)',  2)
      
         ! define planar contact (x,sp)-positions in rail_srfc
    
         call grid_copy(cgrid, rail_srfc)

         ! conformal approach: cgrid has curved sc-coordinates, convert to sp.

         if (ic%is_conformal()) then
            call curved_sc_to_planar_sp(cp%mpot, rail_srfc, cp%curv_ref, idebug)
         endif

         if (is_prismatic) then

            ! convert and evaluate input rail profile 'prr_grd' at contact grid (x,s) coordinates
            ! evaluate local z values at local y values

            call spline_get_loc_xz_at_y( my_rail%prr%grd_data%spl, mpot_rail, rail_srfc, ierror )

            ! call grid_print(rail_srfc, 'rail_srfc', 5)

         elseif (my_rail%prr%is_varprof()) then

            ! convert and evaluate input rail profile 'vprf' at contact grid (x,s) coordinates

            ! evaluate local z values at local y values

            ! call write_log(' varprof_get_loc_z_at_xy...')
            ! call varprof_set_debug(5)
            call varprof_get_loc_z_at_xy( my_rail%prr, ws%s, mpot_rail, rail_srfc, ierror )
            call varprof_set_debug(0)

            ! call grid_print(rail_srfc, 'rail_srfc', 5)

         endif
      endif ! is_prismatic or is_varprof

      if (ic%is_roller()) then

         ! convert input rail profile 'prr_grd' to planar pot.contact coordinates
         !  1. form rail profile 'rail_sn' in planar pot.contact (s,n) coordinates

         ! determine the xtr-positions of the contact grid

         dx_cp = cgrid%dx
         nxlow = nint( (cp%mpot%x() - cp%xsta) / dx_cp ) + 1
         nxhig = nint( (cp%xend - cp%mpot%x()) / dx_cp ) + 1
         nx    = nxlow + nxhig + 1
         xtr_sta = cp%mpot%x() - nxlow * dx_cp

         if (idebug.ge.3) then
            write(bufout,'(3(a,f12.6),a)') ' cref%x(tr)=',cp%mpot%x(),', c.grid has x(tr) in [',        &
                        cp%xsta, ',', cp%xend,']'
            call write_log(1,bufout)
         endif

         ! roller rig, body of revolution: form rail surface at range [xtr_sta, xtr_end]

         ! TODO: performance optimization; restrict to relevant points in master profile

         z_axle = trk%nom_radius
         call grid_revolve_profile(my_rail%prr%grd_data, nx, xtr_sta, dx_cp, z_axle, 999d0, rail_sn_full)

         !  - mrai_pot: rail marker in terms of pot.contact coordinates, ignoring x-position
    
         mrai_pot = marker_2loc(my_rail%m_trk, cp%mpot)
    
         if (idebug.ge.4) then
            call marker_print(my_rail%m_trk, 'm_rail(trk)', 2)
            call marker_print(cp%mpot,       'm_pot(trk)',  2)
            call marker_print(mrai_pot,      'm_rail(pot)', 2)
         endif
      
         ! call grid_print(rail_sn_full,'rail_sn_full',5)
         call cartgrid_2glob(rail_sn_full, mrai_pot)
    
         if (idebug.ge.7) then
            call write_log('rail surface in pot.contact-coordinates:')
            call grid_print(rail_sn_full, 'rail_sn_full(pot)', 5)
         endif
    
         ! create trimmed grid to avoid multi-valued function at large rotation angles

         smin = cgrid%y( 1 )
         smax = cgrid%y( cgrid%ntot )
         ds   = cgrid%dy

         if (idebug.ge.4) then
            write(bufout, *) 'contact grid: s = [', smin,',', smax, '], step ds =', ds
            call write_log(1, bufout)
         endif

         call trim_rail_mesh (rail_sn_full, smin, smax, ds, rail_sn, idebug)
    
         if (idebug.ge.3) then
            call write_log('rail surface in pot.contact-coordinates:')
            call grid_print(rail_sn, 'rail_sn(pot)', 5)
         endif
    
         ! define planar contact (x,sp)-positions in rail_srfc
    
         call grid_copy(cgrid, rail_srfc)

         if (ic%is_conformal()) then
            ! conformal approach: gd%grid has curved sc-coordinates, convert to sp.
            call curved_sc_to_planar_sp(cp%mpot, rail_srfc, cp%curv_ref, idebug)
         endif
    
         ! interpolate rail height 'n' from surface rail_sn to the contact grid sp-locations in rail_srfc
         ! roller, varprof: rail_sn is 3D surface
    
         call interp_cartz2unif( rail_sn, rail_srfc, ierror, 999d0 )
    
         ! check interpolation, warn if rail grid was defined too small

         icount = 0
         do ix = 1, cgrid%nx
            do iy = 1, cgrid%ny
               ii = ix + (iy-1) * cgrid%nx
               if (rail_srfc%z(ii).ge.998d0) icount = icount + 1
            enddo
         enddo

         if (idebug.ge.2 .and. icount.gt.0) then
            write(bufout,'(a,i7,2a)') ' INFO: There are',icount, ' points in the roller with no n/z-value.'
            call write_log(1, bufout)
   
            if (idebug.ge.3) then
               call grid_get_boundbox(rail_sn, bbr)
               write(bufout,248) 'rail', bbr%x(1), bbr%x(8), bbr%y(1), bbr%y(8), bbr%z(1), bbr%z(8)
               call write_log(2, bufout)

               call grid_get_boundbox(cgrid, bbc)
               write(bufout,248) 'pot.contact', bbc%x(1), bbc%x(8), bbc%y(1), bbc%y(8)
               call write_log(2, bufout)

  248          format(' The ',a,' surface has bounding box',/,                                             &
                      '      x=[',f8.3,',',f8.3,'], s=[',f8.3,',',f8.3,']',:,', n=[',f8.3,',',f8.3,']')
            endif

            if (idebug.ge.7) then
               call grid_print(rail_sn, 'rail_sn', 5)
               call grid_print(rail_srfc, 'rail_srfc(pot)', 5)
            endif
         endif

      endif ! is_roller

      if (idebug.ge.7) then
         call write_log('interpolated to pot.contact-grid:')
         call grid_print(rail_srfc, 'rail_srfc(pot)', 5)
      endif
    
      if (idebug.ge.6) then
         call write_log('rail n-coordinate:')
         do iy = 1, cgrid%ny
            ii = 1 + (iy-1)*cgrid%nx
            write(bufout,*) 'iy=',iy,': ', rail_srfc%z(ii)
            call write_log(1, bufout)
         enddo
      endif

      !  2. compute the extent of the potential contact area in wheel-set coordinates

      ! express pot.contact reference in terms of wheel profile coordinates

      mwhl_trk  = marker_2glob( my_wheel%m_ws, ws%m_trk )
      mpot_whl   = marker_2loc( cp%mpot, mwhl_trk )

      if (idebug.ge.2) then
         call write_log(' pot.contact origin w.r.t. O_trk')
         call marker_print(cp%mpot, 'm_pot(trk)', 5)
         call write_log(' pot.contact origin w.r.t. O_whl')
         call marker_print(mpot_whl,  'm_pot(whl)',  5)
      endif

      ! the potential contact area is at [xsta,xend] x [ysta,yend] in track coordinates

      zsta = cp%mpot%z() + cp%sp_sta * sin(cp%delttr)
      zend = cp%mpot%z() + cp%sp_end * sin(cp%delttr)

      if (idebug.ge.2) then
         write(bufout,'(5(a,:,f12.6))') ' zsta=',cp%mpot%z(),' + ',cp%sp_sta, ' * sin(',cp%delttr,')'
         call write_log(1, bufout)
         call write_log(' potential contact area w.r.t. O_trk')
         write(bufout,'(2(a,f12.6),a)') '     x(tr) = [',cp%xsta,'--',cp%xend,']'
         call write_log(1, bufout)
         write(bufout,'(2(a,f12.6),a)') '     y(tr) = [',cp%ysta,'--',cp%yend,']'
         call write_log(1, bufout)
         write(bufout,'(2(a,f12.6),a)') '     z(tr) = [',   zsta,'--',   zend,']'
         call write_log(1, bufout)
      endif

      ! compute potential contact on wheel profile according to potential contact on track

      vtmp(1) = vec_2loc( vec(cp%xsta, cp%ysta, zsta), mwhl_trk )
      vtmp(2) = vec_2loc( vec(cp%xsta, cp%yend, zend), mwhl_trk )
      vtmp(3) = vec_2loc( vec(cp%xend, cp%ysta, zsta), mwhl_trk )
      vtmp(4) = vec_2loc( vec(cp%xend, cp%yend, zend), mwhl_trk )

      xmin = min(vtmp(1)%x(), vtmp(2)%x(), vtmp(3)%x(), vtmp(4)%x()) - 2*cgrid%dx
      xmax = max(vtmp(1)%x(), vtmp(2)%x(), vtmp(3)%x(), vtmp(4)%x()) + 2*cgrid%dx
      ymin = min(vtmp(1)%y(), vtmp(2)%y(), vtmp(3)%y(), vtmp(4)%y()) - 2*cgrid%dy
      ymax = max(vtmp(1)%y(), vtmp(2)%y(), vtmp(3)%y(), vtmp(4)%y()) + 2*cgrid%dy
      zmin = min(vtmp(1)%z(), vtmp(2)%z(), vtmp(3)%z(), vtmp(4)%z()) -   cgrid%dy * abs(sin(cp%delttr))
      zmax = max(vtmp(1)%z(), vtmp(2)%z(), vtmp(3)%z(), vtmp(4)%z()) +   cgrid%dy * abs(sin(cp%delttr))

      if (idebug.ge.2) then

         call write_log(' potential contact area w.r.t. O_whl')
         write(bufout,'(3(a,f12.6))')   '     x(whl) = [',xmin,'--',xmax,'] on wheel, pot at ', mpot_whl%x()
         call write_log(1, bufout)
         write(bufout,'(3(a,f12.6))')   '     y(whl) = [',ymin,'--',ymax,'] on wheel, pot at ', mpot_whl%y()
         call write_log(1, bufout)
         write(bufout,'(2(a,f12.6),a)') '     z(whl) = [',zmin,'--',zmax,']'
         call write_log(1, bufout)
      endif

      ! compute wheel surface in wheel profile coordinate system

      if (.not.my_wheel%prw%is_varprof()) then

         !  - the wheel profile prw_grd is given in wheel profile coordinates

         associate(prw_grd => my_wheel%prw%grd_data)
         nw = prw_grd%ntot

         ! find range of sw on wheel that covers [ymin, ymax]

         if (grid_has_spline(prw_grd)) then
            ! call write_log(' prw_grd has spline')
            ! call grid_print(prw_grd, 'prw_grd', 5)
            ! call spline_set_debug(4)
            call spline_get_s_at_y( prw_grd%spl, ymin,        swsta, ierror )
            call spline_get_s_at_y( prw_grd%spl, mpot_whl%y(), swmid, ierror )
            call spline_get_s_at_y( prw_grd%spl, ymax,        swend, ierror )
            ! call spline_set_debug(0)
         else
            call write_log(' Internal Error (undefdist): no spline available')
            call abort_run()
         endif

         ! swsta = swsta + 10d0 * cgrid%dy
         ! swend = swend - 10d0 * cgrid%dy

         if (idebug.ge.2) then
            write(bufout,'(3(a,f12.6))') ' pot.contact  y(whl) = [',ymin ,'--',ymax ,'] on wheel, mid=', &
                   mpot_whl%y()
            call write_log(1, bufout)
            write(bufout,'(3(a,f12.6))') ' search range s(whl) = [',swsta,'--',swend,'] on wheel, mid=',swmid
            call write_log(1, bufout)
         endif

         ! select points in [sw_sta, sw_end], refine using spline interpolation

         call trim_wheel_profile(prw_grd, swsta, swmid, swend, cgrid%dy, prw_trim, idebug)
         end associate

         ! round to grid with step of dx centered at m_pot (should still encompass pot.contact?)
         !  - define grid with xslc(i') = cp%x + i' * dx , for i' = -nslow : nshig

         dx_cp = cgrid%dx
         nslow = nint(  (mpot_whl%x() - xmin) / dx_cp )
         nshig = nint(  (xmax - mpot_whl%x()) / dx_cp )
         nslc  = nslow + nshig + 1
         xslc1 = mpot_whl%x() - nslow * dx_cp

         !  3. form wheel mesh in wheel-set coordinates using nslc slices, with resolution dx
   
         !  - initialize surface whl_srfw_full for the wheel height z_cntc above the contact plane
         !  - form wheel mesh, selecting points from profile

         z_axle = -ws%nom_radius
         call grid_revolve_profile(prw_trim, nslc, xslc1, dx_cp, z_axle, -999d0, whl_srfw_full)

      else ! my_wheel%prw%is_varprof()

         if (ic%use_oblique()) then
            call write_log(' Variable wheel (slcw) does not support oblique view direction.')
            call abort_run()
         endif

         ! convert corners vtmp(1:4) into cylindrical wcyl coordinates

         if (cp%usta.lt.cp%uend) then   ! [usta,uend] not provided: empty range [1,0]
            th_min = cp%usta
            th_max = cp%uend
         else
            rmin   = ws%nom_radius + zmin
            rmax   = ws%nom_radius + zmax
            th_min = -ws%pitch + atan(xmin/rmin)
            th_max = -ws%pitch + atan(xmax/rmin)
         endif

         if (idebug.ge.2) then
            write(bufout,'(3(a,f8.3))') ' xmin=',xmin,', th_min=',th_min,', u_min=',cp%usta
            call write_log(1, bufout)
            write(bufout,'(3(a,f8.3))') ' xmax=',xmin,', th_max=',th_max,', u_max=',cp%uend
            call write_log(1, bufout)
         endif

         ! compute nslc such that dx_ws = O(dx_cntc), dx = r * dth

         dx_cp = cgrid%dx
         nslc = nint( (th_max - th_min) * ws%nom_radius / (1d0 * dx_cp) ) + 2
         if (mod(nslc,2).eq.0) nslc = nslc + 1
         dth  = (th_max - th_min) / max(1d0, real(nslc - 1))

         ! set spline u-parameter for rolling direction

         allocate(u_out(nslc))
         do iu = 1, nslc
            u_out(iu) = wrap_around(th_min + (iu-1) * dth)
         enddo
         if (idebug.ge.2) then
            write(bufout,'(a,i4,a,3f12.6)') ' u: u_out([1,2,',nslc,'])=',u_out(1), u_out(2), u_out(nslc)
            call write_log(1, bufout)
         endif

         ! set spline v-parameter for lateral direction

         associate( spl2d => my_wheel%prw%spl2d )

         if (cp%vsta.lt.cp%vend) then   ! [vsta,vend] not provided: empty range [1,0]
            v_min = cp%vsta
            v_max = cp%vend
            nw    = nint( (cp%yend - cp%ysta) / (cgrid%dy*cos(cp%delttr)) )
            if (idebug.ge.2) then
               write(bufout,'(3(a,f10.3),2(a,f6.3),a,i4)') ' ysta=',cp%ysta,', yend=',cp%yend,' wid=', &
                        cp%yend-cp%ysta,', dy=', cgrid%dy,', cos=',cos(cp%delttr),', nw=',nw
               call write_log(1, bufout)
            endif
         else
            ! nw = my_wheel%prw%grd_data%ntot
            nw    = spl2d%nknotv - 6
            v_min = spl2d%tvj(4)
            v_max = spl2d%tvj(spl2d%nknotv - 3)
         endif
         dv  = (v_max - v_min) / real(nw - 1)

         allocate(v_out(nw))
         do iv = 1, nw
            v_out(iv) = v_min + (iv-1) * dv
         enddo
         if (idebug.ge.2) then
            write(bufout,'(2(a,f12.6),a,i4,a,f12.6)') ' v=[',v_min,',',v_max,'], nw=',nw,', dv=',dv
            call write_log(1, bufout)
            !write(bufout,'(a,i4,a,3f12.6)') ' v: v_out([1,2,',nw,'])=',v_out(1), v_out(2), v_out(nw)
            !call write_log(1, bufout)
         endif

         ! create storage for curvilinear grid with cylindrical coordinates

         call grid_create_curvil(whl_srfw_full, nslc, nw, cyl_coords=.true.)
         ! call grid_print(whl_srfw_full, 'cyl_grid', 2)

         ! evaluate 2D spline at (u,v) to get (th, y, dr)

         call bspline_eval2d_prod(spl2d, nslc, nw, u_out, v_out, whl_srfw_full%th, whl_srfw_full%y,     &
                        whl_srfw_full%r, .true., ierror, -987d0)
         deallocate(u_out, v_out)

         ! add nominal radius r = rnom + dr

         call grid_shift(whl_srfw_full, 0d0, 0d0, ws%nom_radius)

         ! call grid_print(whl_srfw_full, 'cyl', 5)

         ! convert (th, y, r) to (x, y, z) in wheel center coordinates

         call convert_cyl2cart(whl_srfw_full)

         ! bring -ws%pitch to the lowest point

         call cartgrid_pitch(whl_srfw_full, ws%pitch, 0d0, 0d0)

         ! shift to wheel profile coordinates (z=-rnom at height of axle)

         call grid_shift(whl_srfw_full, 0d0, 0d0, -ws%nom_radius)

         if (.false.) then
            call marker_print( my_wheel%m_ws, 'm_ws', 5)
            call grid_print(whl_srfw_full, 'whl_srfw_full', 5)
            call abort_run()
         endif

         end associate

      endif ! my_wheel%prw%is_varprof()

      ! express wheel profile marker in terms of pot.contact coordinates

      mwhl_trk = marker_2glob( my_wheel%m_ws, ws%m_trk )
      mwhl_pot = marker_2loc( mwhl_trk, cp%mpot )

      ! 4. convert wheel surface from wheel profile coordinates to planar pot.contact coordinates

      ! call grid_print(whl_srfw_full, 'curv', 5)
      call cartgrid_2glob( whl_srfw_full, mwhl_pot )

      ! TODO: create trimmed grid to avoid multi-valued function at large rotation angles

      call grid_copy(whl_srfw_full, whl_srfw)
      ! call grid_print(whl_srfw, 'whl', 5)

      ! copy the rail mesh containing the [x,sp]-points of the pot.contact grid

      call grid_copy(rail_srfc, whl_srfc)

      ! 5. interpolate the wheel n-coordinates to the [x,sp]-points of the pot.contact grid

      use_bicubic = (.true. .and. my_wheel%prw%is_varprof())
      ! call interp_set_debug(3, 1, 1)
      call timer_start(itimer_interp9)
      call interp_cartz2unif(whl_srfw, whl_srfc, ierror, -999d0, use_bicubic)
      call timer_stop(itimer_interp9)
      ! call interp_set_debug(0, 1, 1)

      ! check interpolation, warn if wheel grid was defined too small

      icount = 0
      do ix = 1, cgrid%nx
         do iy = 1, cgrid%ny
            ii = ix + (iy-1) * cgrid%nx
            if (abs(whl_srfc%z(ii)).ge.998d0) icount = icount + 1
         enddo
      enddo

      if (icount.gt.0) then
         if (idebug.ge.2) then
            write(bufout,'(a,i7,2a)') ' INFO: There are',icount, ' points in the wheel with no n/z-value.'
            call write_log(1, bufout)
         endif

         if (idebug.ge.2) then
            call grid_get_boundbox(whl_srfw, bbw)
            write(bufout,248) 'wheel', bbw%x(1), bbw%x(8), bbw%y(1), bbw%y(8), bbw%z(1), bbw%z(8)
            call write_log(2, bufout)

            call grid_get_boundbox(whl_srfc, bbc)
            write(bufout,248) 'pot.contact', bbc%x(1), bbc%x(8), bbc%y(1), bbc%y(8)
            call write_log(2, bufout)
         endif

         if (idebug.ge.4) then
            call write_log(' Data input to interpolation:')
            call grid_print(whl_srfw, 'whl_srfw(pot)', 5)
            call write_log(' Interpolation result:')
            call grid_print(whl_srfc, 'whl_srfc(pot)', 5)
         endif
      endif

      if (idebug.ge.7) then
         call write_log('wheel n-coordinate:')
         ix = 10
         do iy = 1, cgrid%ny
            ii = ix + (iy-1) * cgrid%nx
            write(bufout,*) 'iy=',iy,': ',whl_srfc%z(ii)
            call write_log(1, bufout)
         enddo
      endif

      ! the rail and wheel meshes are both computed on the planar pot.contact grid

      ! 6. subtract rail - wheel to get the gap, negative == interpenetration

      ixdbg = (cgrid%nx+1) / 2
      iydbg = (cgrid%ny+1) / 2

      hmin = 999d0
      vy_ref =  sin(cp%delttr)
      vz_ref = -cos(cp%delttr)

      if (idebug.ge.3 .and. ic%is_conformal()) then
         write(bufout,'(3(a,f10.6),a)') ' delttr=',cp%delttr,', nrm=[',vy_ref,',',vz_ref,']'
         call write_log(1, bufout)
      endif

      ! conformal: multiply h_ii in np-direction by cos(delt_act-delt_ref) = 
      !                                      dot(n_ref, n_actual) / (|n_ref| * |n_actual|)

      do iy = 1, cgrid%ny

         if (ic%is_conformal()) then
            ii = 1 + (iy-1) * cgrid%nx
            vy_iy = cp%curv_nrm%vy(iy)
            vz_iy = cp%curv_nrm%vn(iy)
            vlen_inv = 1d0 / max(1d-8, sqrt( vy_iy**2 + vz_iy**2 ))

            if (idebug.ge.3) then
               write(bufout,'(a,i4,3(a,f10.6))') ' iy=',iy,': nrm=[',vy_iy,',',vz_iy,'] /', 1d0/vlen_inv
               call write_log(1, bufout)
            endif
         endif

         do ix = 1, cgrid%nx
            ii = ix + (iy-1) * cgrid%nx
            h_ii = rail_srfc%z(ii) - whl_srfc%z(ii)
            if (ic%is_conformal()) h_ii = h_ii * (vy_ref*vy_iy + vz_ref*vz_iy) * vlen_inv
            if (h_ii.lt.hmin) hmin = h_ii
            gd%geom%prmudf(ii) = h_ii

            if (ix.eq.ixdbg .and. idebug.ge.7) then
               write(bufout,'(a,i5,6(a,f12.6))') ' ii=',ii,': rn=',rail_srfc%z(ii),', wn=',whl_srfc%z(ii), &
                        ' hpl=', rail_srfc%z(ii)-whl_srfc%z(ii),', vy=',vy_iy,', vz=',vz_iy,', fac=',      &
                        (vy_ref*vy_iy+vz_ref*vz_iy)*vlen_inv
               call write_log(1, bufout)
            endif

         enddo
      enddo

      ! 7. put the minimum of h in the approach pen, shift h such that minimum is 0

      do iy = 1, cgrid%ny
         do ix = 1, cgrid%nx
            ii = ix + (iy-1) * cgrid%nx
            gd%geom%prmudf(ii) = gd%geom%prmudf(ii) - hmin
         enddo
      enddo
      gd%kin%pen = -hmin

      if (idebug.ge.2) then
         write(bufout,'(a,g12.4)') ' h_min = min(ud) =', hmin
         call write_log(1, bufout)
      endif

      if (.true. .and. idebug.ge.2) then
         write(tmp_fname,'(a,i1,a)') 'dump_gap_ud',icp,'.m'
         call write_log(' Writing undef.distc to ' // trim(tmp_fname))
         call make_absolute_path(tmp_fname, meta%outdir, tmp_fulnam)
         lud = get_lunit_tmp_use()
         open(unit=lud, file=tmp_fulnam, iostat=ios, err=991)
         write(lud,'(2(a,i6),a)') 'mx=',cgrid%nx,'; my=',cgrid%ny,';'
         write(lud,*) 'tmp=['
         do iy = 1, cgrid%ny
            if (ic%is_conformal()) then
               vy_iy = cp%curv_nrm%vy(iy)
               vz_iy = cp%curv_nrm%vn(iy)
               delt_iy = atan(vz_iy/vy_iy)
            else
               delt_iy = cp%delttr
            endif
            do ix = 1, cgrid%nx
               ii = ix + (iy-1) * cgrid%nx
               write(lud,'(2i5,7g15.6)') ix, iy, cgrid%x(ii), cgrid%y(ii), whl_srfc%y(ii),       &
                        whl_srfc%z(ii), rail_srfc%z(ii), delt_iy, gd%geom%prmudf(ii)
            enddo
         enddo
         write(lud,*) '];'
         write(lud,*) 'x_ud =reshape(tmp(:,3),mx,my);'
         write(lud,*) 'sp_ud=reshape(tmp(:,4),mx,my);'
         write(lud,*) 'sc_ud=reshape(tmp(:,5),mx,my);'
         write(lud,*) 'nw_ud=reshape(tmp(:,6),mx,my);'
         write(lud,*) 'nr_ud=reshape(tmp(:,7),mx,my);'
         write(lud,*) 'delta=reshape(tmp(:,8),mx,my);'
         write(lud,*) 'h_ud =reshape(tmp(:,9),mx,my);'
         close(lud)
         call free_lunit_tmp_use(lud)
         goto 999

 991     continue
            write(bufout,'(2(a,i6))') ' Error opening lud=',lud,', ios=',ios
            call write_log(1, bufout)

 999     continue
      endif

      ! clean up local variables

      call grid_destroy(rail_sn_full)
      call grid_destroy(rail_sn)
      call grid_destroy(prw_trim)
      call grid_destroy(whl_srfw_full)
      call grid_destroy(whl_srfw)
      call grid_destroy(bbr)
      call grid_destroy(bbc)
      call grid_destroy(bbw)

      end associate

      call timer_stop(itimer_udist)
      if (idebug.ge.4) call write_log('--- end subroutine wr_ud_planar ---')

   end subroutine wr_ud_planar

!------------------------------------------------------------------------------------------------------------

   subroutine curved_sc_to_planar_sp (mpot, cgrid, curv_ref, idebug)
!--purpose: convert curved sc-coordinates in cgrid to planar sp-coordinates
      implicit none
!--subroutine arguments:
      integer                   :: idebug
      type(t_marker)            :: mpot
      type(t_grid)              :: cgrid, curv_ref
!--local variables:
      integer                   :: ix, iy, ii
      type(t_grid)              :: csrf_pot

      if (cgrid%ny.ne.curv_ref%ny) then
         write(bufout,'(a,2i5)') ' Internal ERROR: curved_to_planar: grids dont match, ny=', cgrid%ny,  &
                curv_ref%ny
         call write_log(1, bufout)
         call abort_run()
      endif

      ! input curv_ref is given in trk-coordinates. Copy and convert to planar pot.con-coords

      call grid_copy(curv_ref, csrf_pot)
      call cartgrid_2loc(csrf_pot, mpot)

      if (idebug.ge.2) call write_log(' converting sc-coordinates of contact grid to sp-values')
      if (idebug.ge.3) then
         do iy = 1, cgrid%ny
            ii = 1 + (iy-1) * cgrid%nx
            write(bufout,'(a,i3,2(a,f12.4))') ' iy=',iy,': sc=',cgrid%y(ii),' - sp=',csrf_pot%y(iy)
            call write_log(1, bufout)
         enddo
      endif

      ! overwrite curved sc-coordinates in cgrid by planar sp-coordinates

      do iy = 1, cgrid%ny
         do ix = 1, cgrid%nx
            ii = ix + (iy-1) * cgrid%nx
            cgrid%y(ii) = csrf_pot%y(iy)
         enddo
      enddo

      ! grid 'cgrid' is no longer uniform

      cgrid%is_uniform  = .false.
      cgrid%is_curvilin = .true.

      call grid_destroy(csrf_pot)

   end subroutine curved_sc_to_planar_sp

!------------------------------------------------------------------------------------------------------------

   subroutine trim_rail_mesh (rail_sn_full, smin, smax, ds, rail_sn_trim, idebug)
!--purpose: define a trimmed version of rail mesh for use in interpolation
!           the input mesh 'rail_sn_full' is defined in (x,s,n)-coordinates, with contact reference
!           at (0,0,0). The parameters [smin, smax] give the expected range of s-coordinates.
!           n-coordinates are expected to be small in the contact region.
      implicit none
!--subroutine arguments:
      integer                   :: idebug
      real(kind=8), intent(in)  :: smin, smax, ds
      type(t_grid)              :: rail_sn_full, rail_sn_trim
!--local variables:
      logical, dimension(:), allocatable :: mask
      integer               :: ii, ix, iy, ii_min, ix_min, iy_min, ix_sta, ix_end, iy_sta, iy_end,      &
                               nguard, icheck
      real(kind=8)          :: ds_tot, smin_loc, smax_loc, n_min

      associate(nx     => rail_sn_full%nx, ny     => rail_sn_full%ny,                                   &
                rail_s => rail_sn_full%y,  rail_n => rail_sn_full%z )

      ! find a centered point: minimum n-value
    
      if (idebug.ge.4) call write_log('--- trim_rail: locate overall minimum n-value ---')

      ii_min = idmin(nx*ny, rail_n, 1)
      ix_min = rail_sn_full%ix(ii_min)
      iy_min = rail_sn_full%iy(ii_min)
      n_min  = rail_sn_full%z(ii_min)

      ! make sure "centered point" is included in range

      ds_tot   = abs(smax - smin)
      smin_loc = min(smin, rail_sn_full%y(ii_min)-0.05*ds_tot)
      smax_loc = max(smax, rail_sn_full%y(ii_min)+0.05*ds_tot)

      if (idebug.ge.2) then
         write(bufout,'(2(a,i5),a,f11.6,a,i7,2(a,i5))') ' rail_sn_full: nx=',nx,', ny=',ny,             &
                ', n_min=',n_min,' at ii=',ii_min,', ix=',ix_min, ', iy=',iy_min
         call write_log(1, bufout)
         write(bufout,'(3(a,f8.3))') ' smin=',smin,', smax=',smax,', ds=',ds
         call write_log(1, bufout)
      endif

      ! make mask for points that could be included, if connected to the centered point
      ! TODO: exclude multiple points at same (x,s)-position (multi-valued function)?

      allocate(mask(nx*ny))

      ! set mask at iy = iy_min

      iy = iy_min
      do ix = 1, nx
         ii = ix + (iy-1)*nx
         mask(ii) = (rail_s(ii).ge.smin_loc .and. rail_s(ii).le.smax_loc)
      enddo

      ! expand mask to the left

      do iy = iy_min-1, 1, -1
         do ix = 1, nx
            ii = ix + (iy-1)*nx
            mask(ii) = (mask(ii+nx) .and. rail_s(ii).ge.smin_loc .and. rail_s(ii).le.smax_loc)
            ! points iy < iy_min:
            if (rail_s(ii).gt.rail_s(ii+nx)) mask(ii) = .false.
            if (idebug.ge.3 .and. ix.eq.1) then
               write(bufout,*) 'ix,iy=',ix,iy,': ii=',ii,', s=',rail_s(ii),', mask=',mask(ii)
               call write_log(1, bufout)
            endif
         enddo
      enddo

      ! expand mask to the right

      do iy = iy_min+1, ny
         do ix = 1, nx
            ii = ix + (iy-1)*nx
            mask(ii) = (mask(ii-nx) .and. rail_s(ii).ge.smin_loc .and. rail_s(ii).le.smax_loc)
            ! points iy > iy_min:
            if (rail_s(ii).le.rail_s(ii-nx)) mask(ii) = .false.
            if (idebug.ge.3 .and. ix.eq.1) then
               write(bufout,*) 'ix,iy=',ix,iy,': ii=',ii,', s=',rail_s(ii),', mask=',mask(ii)
               call write_log(1, bufout)
            endif
         enddo
      enddo

      ! determine range of indices where mask is set

      if (nx.eq.1) then

         ! rail: determine first and last indices where mask is set

         ix_sta = 1
         ix_end = 1
         iy_sta = 100000
         iy_end = 1
         do iy = 1, ny
            if (mask(iy) .and. iy_sta.gt.99999) iy_sta = iy
            if (mask(iy)) iy_end = iy
         enddo

         ! add guard band: the profile must fully encompass the contact grid

         if (iy_sta.gt.1) then
            if (rail_s(iy_sta-1).le.rail_s(iy_sta)) iy_sta = iy_sta - 1
         endif
         if (iy_end.lt.ny) then
            if (rail_s(iy_end+1).ge.rail_s(iy_end)) iy_end = iy_end + 1
         endif

      else

         ! roller: allow decreasing sp within range & at guard band

         nguard = 1
         icheck = 0
         call find_bounding_box_2d( mask, nx, ny, ix_min, iy_min, nguard, ix_sta, ix_end, iy_sta,    &
                                    iy_end, icheck, idebug)
      endif
      deallocate(mask)

      ! create trimmed copy of rail_sn

      if (idebug.ge.3) then
         write(bufout,'(4(a,i4),a)') ' trim_rail_mesh: selecting ix=[',ix_sta,'-', ix_end,           &
             '], iy=[',iy_sta,'-', iy_end,']'
         call write_log(1, bufout)
      endif
    
      call grid_trim(rail_sn_full, rail_sn_trim, ix_sta, ix_end, iy_sta, iy_end)

      end associate
   end subroutine trim_rail_mesh
    
!------------------------------------------------------------------------------------------------------------

   subroutine trim_wheel_profile(prw_grd, sw_sta, sw_mid, sw_end, ds, prw_trim, idebug)
!--purpose: define a trimmed version of the wheel profile using spline interpolation
      implicit none
!--subroutine arguments:
      integer                   :: idebug
      real(kind=8), intent(in)  :: sw_sta, sw_mid, sw_end, ds
      type(t_grid)              :: prw_grd, prw_trim
!--local variables:
      integer                   :: nx, ns, is, is_mid, ierror
      real(kind=8), dimension(:), allocatable :: s, x, y, z

      ! compute number of points for s-grid; note: s-values decreasing, sw_sta > sw_mid > sw_end
      ! add one step on each side, ns = number of steps + 1. Trimmed mesh: [sw_end-ds : ds : sw_sta+ds]

      nx     = 1
      ns     = 1 + nint( (sw_sta-sw_mid)/ds ) + 1 + nint( (sw_mid-sw_end)/ds ) + 1
      is_mid = 1 + nint( (sw_sta-sw_mid)/ds ) + 1       ! number of sw_mid

      if (idebug.ge.4) then
         call write_log(' trim_prw: input profile:')
         call grid_print( prw_grd, 'prw(whl)', 5)
      endif

      if (idebug.ge.2) then
         write(bufout,'(2(a,f11.6),a,i6)') ' target wheel s-grid: sw_sta=',sw_sta,', sw_end=',sw_end,   &
                        ', ns=',ns
         call write_log(1, bufout)
      endif

      allocate(s(ns), x(ns), y(ns), z(ns))

      do is = 1, ns
         s(is) = sw_mid - real(is-is_mid) * ds
         x(is) = 0d0
      enddo

      if (idebug.ge.2) then
         write(bufout,'(3(a,f11.6))') ' wheel s-grid: s=',s(1),' :',-ds,' :',s(ns)
         call write_log(1, bufout)
      endif

      ! compute (y,z)-positions at s-grid

      call spline_get_xyz_at_s( prw_grd%spl, ns, s, ierror, yout=y, zout=z )

      ! create grid from x, y, z-arrays

      call grid_create_curvil( prw_trim, nx, ns, x, y, z)

      if (idebug.ge.4) then
         call grid_make_arclength( prw_trim, ierror )
         call write_log(' trim_prw: output profile:')
         call grid_print( prw_trim, 'prw(whl)', 5)
      endif

      deallocate(s, x, y, z)

   end subroutine trim_wheel_profile
    
!------------------------------------------------------------------------------------------------------------

   subroutine trim_wheel_mesh (wheel_sn_full, smin, smax, wheel_sn_trim, idebug)
!--purpose: define a trimmed version of wheel mesh for use in interpolation
      implicit none
!--subroutine arguments:
      integer                   :: idebug
      real(kind=8), intent(in)  :: smin, smax
      type(t_grid)              :: wheel_sn_full, wheel_sn_trim
!--local variables:
      logical, dimension(:), allocatable :: mask
      integer               :: ii, iimin, ix_min, iy_min, ix_sta, ix_end, iy_sta, iy_end, nguard, icheck
      real(kind=8)          :: mindst2

      associate(nx      => wheel_sn_full%nx, ny      => wheel_sn_full%ny,                               &
                wheel_x => wheel_sn_full%x,  wheel_s => wheel_sn_full%y,   wheel_n => wheel_sn_full%z )

      ! mark all point in the range [smin, smax], find central point: minimum || [x,s] ||^2
      ! TODO: exclude multiple points at same (x,s)-position (multi-valued function)?
 
      allocate(mask(nx*ny))

      iimin   = 0
      mindst2 = 1d10

      do ii = 1, nx*ny
         mask(ii) = (wheel_s(ii).ge.smin .and. wheel_s(ii).le.smax)
         if (wheel_x(ii)**2 + wheel_s(ii)**2 .lt. mindst2) then
            iimin = ii
            mindst2 = wheel_x(ii)**2 + wheel_s(ii)**2
         endif
      enddo

      ix_min = mod(iimin-1,nx) + 1
      iy_min = (iimin-ix_min) / nx + 1

      ! determine range of indices where mask is set

      nguard = 1
      icheck = 0
      call find_bounding_box_2d( mask, nx, ny, ix_min, iy_min, nguard, ix_sta, ix_end, iy_sta,          &
                                 iy_end, icheck, idebug)
      deallocate(mask)

      ! create trimmed copy of wheel_sn

      if (idebug.ge.3) then
         write(bufout,'(4(a,i4),a)') ' trim_wheel_mesh: selecting ix=[',ix_sta,'-', ix_end,           &
             '], iy=[',iy_sta,'-', iy_end,']'
         call write_log(1, bufout)
      endif
 
      call grid_trim(wheel_sn_full, wheel_sn_trim, ix_sta, ix_end, iy_sta, iy_end)

      end associate
   end subroutine trim_wheel_mesh
    
!------------------------------------------------------------------------------------------------------------

end module m_wr_undefdist
