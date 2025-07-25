!------------------------------------------------------------------------------------------------------------
! m_wr_solvecp - solve one case with w/r contact (module 1) with prescribed position
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_wr_solvecp

use m_wrprof_data
use m_sdis
use m_scontc
use m_soutpt
use m_subsurf
use m_wr_profiles
use m_wr_locatecp
use m_wr_undefdist
use m_wr_rigslip
use m_wr_output
implicit none
private

public  wr_contact_pos
public  wr_subsurf
public  wr_setup_cp
private compute_sc_for_yz
public  wr_solve_cp
private aggregate_forces
private cpatch_forces_moments
public  total_forces_moments
public  set_dummy_solution
private wtd_check_nan
private merge_prev_potcon
private extend_curv_ref

contains

!------------------------------------------------------------------------------------------------------------

   subroutine wr_contact_pos(wtd, x_locate, my_ierror)
!--purpose: process the w/r problem for a case with positions & velocities prescribed
      implicit none
!--subroutine arguments:
      type(t_ws_track)              :: wtd
      integer,          intent(in)  :: x_locate
      integer,          intent(out) :: my_ierror
!--local variables:
      integer,     parameter  :: check_nan = 1
      logical                 :: is_varprof
      integer                 :: icp, ldebug, sub_ierror
      type(t_grid)            :: sf_rai

      my_ierror = 0

      if (x_locate.ge.2) call write_log(' ')
      if (x_locate.ge.2) call write_log(' --- Start subroutine wr_contact_pos ---')

      ! suppress warnings for rail/wheel profiles in subsequent total force iterations (Brent, Secant)

      wtd%meta%itforce = wtd%meta%itforce + 1

      if (wtd%ic%x_cpatch.ge.2) then
         write(bufout,'(3(a,i3))') ' itout =',wtd%meta%itforc_out,', itinn =',wtd%meta%itforc_inn,       &
                   ', ittot =',wtd%meta%itforce
         call write_log(1, bufout)
      endif

      ldebug = x_locate
      if (wtd%meta%itforc_out.ge.1 .or. wtd%meta%itforc_inn.ge.1) ldebug = ldebug - 1

      call timer_start(itimer_wrgeom)

      ! in the first outer iteration: select the rail slice when using a variable profile

      is_varprof = (wtd%trk%rai%prr%nslc.gt.0)

      if (wtd%meta%itforc_out.le.0 .and. wtd%meta%itforc_inn.le.0 .and. is_varprof) then
         ! call profile_select_slice(wtd%trk%rai%prr, wtd%ws%s, ldebug)

         if (x_locate.ge.2) call write_log(' setting up "current slice" for variable profile')
         ! call bspline2d_print(wtd%trk%rai%prr%spl2d, 'vprf.spl2d', 5)
         ! call varprof_set_debug(3)
         ! call bspline_set_debug(3)
         call varprof_intpol_xunif(wtd%trk%rai%prr, wtd%ws%s, 1, 0d0, 0d0, wtd%trk%rai%prr%grd_data, &
                        sub_ierror)
         ! call bspline_set_debug(0)
         ! call varprof_set_debug(0)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         ! call grid_print(wtd%trk%rai%prr%grd_data, 'prr', 5)
         ! call abort_run()
      endif

      if (is_varprof .and. .false.) then
         ! test code for new varprof_intpol_grid, showing wrong resuls due to compiler bug(?)
         call marker_print(wtd%trk%rai%m_trk, 'm_r(tr)', 5)
         call write_log(' ...varprof_intpol_grid')
         call grid_create_uniform(sf_rai, x0arg=0d0, dxarg=1d0, x1arg=0d0,                     &
                                         y0arg=10d0, dyarg=5d0, y1arg=20d0, zarg=0d0)
         call bspline_set_debug(1)
         call varprof_intpol_grid(wtd%trk%rai%prr, wtd%trk%rai%m_trk, wtd%ws%s, sf_rai)
         call bspline_set_debug(0)
         call grid_print(sf_rai, 'sf_rai_new', 5)

         call varprof_intpol_grid_old(wtd%trk%rai%prr, wtd%trk%rai%m_trk, wtd%ws%s, sf_rai)
         call grid_print(sf_rai, 'sf_rai_old', 5)
         call abort_run()
      endif

      ! set the rail marker, position in track coordinate system

      call wr_set_rail_marker (wtd%trk, wtd%ic, ldebug)

      ! set the wheel marker in ws coordinates, wheelset marker in track coordinate system

      call wr_set_wheel_marker (wtd%ws, wtd%ic, ldebug)

      ! locate the contact patch(es), set potential contact areas

      if (my_ierror.eq.0) then
         if (wtd%ic%x_cpatch.ge.1) then
            if (wtd%ic%force1.eq.3) then
               write(bufout,'(2(a,f9.4))') ' dy_defl=',wtd%trk%dy_defl,', z_ws=',wtd%ws%z
            else
               write(bufout,'(a,f9.4)') ' z_ws=',wtd%ws%z
            endif
            call write_log(1, bufout)
         endif
         call wr_locatecp(wtd%meta, wtd%ic, wtd%ws, wtd%trk, wtd%discr, wtd%numcps, wtd%n_miss,         &
                        wtd%numtot, wtd%allcps, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
      endif

         ! define the contact problems for all contact patches

      if (my_ierror.eq.0) then
         do icp = 1, wtd%numcps
            if (.not.associated(wtd%allcps(icp)%cp)) then
               write(bufout,'(a,i3,a)') ' Internal error: allcps(',icp,') not associated.'
               call write_log(1, bufout)
               call abort_run()
            endif
            associate( cp => wtd%allcps(icp)%cp )
            call wr_setup_cp(wtd%meta, wtd%ws, wtd%trk, wtd%numcps, icp, cp, wtd%ic, wtd%mater,         &
                             wtd%discr, wtd%fric, wtd%kin, wtd%solv, x_locate)
            end associate
         enddo
      endif
      call timer_stop(itimer_wrgeom)

      ! perform actual computation when R=0 or 1

      if (my_ierror.eq.0 .and. wtd%ic%return.le.1) then

         do icp = 1, wtd%numcps

            ! solve contact problem for contact patch icp

            associate( cp => wtd%allcps(icp)%cp )
            call wr_solve_cp(wtd%meta, wtd%ws, wtd%trk, icp, cp, wtd%ic, x_locate, sub_ierror)
            if (my_ierror.eq.0) my_ierror = sub_ierror

            if (wtd%ic%x_nmdbg.ge.3) then
               write(bufout,*) 'contact patch icp=',icp,':'
               call write_log(1, bufout)
               call wrigs(cp%gd%outpt1%igs, .true., 0d0)
            endif
            end associate

         enddo

         call total_forces_moments(wtd)

         ! if requested: check for NaNs, write diagnostic information

         if (check_nan.ge.1) then
            if (x_locate.ge.5) call write_log(' ...wr_contact: calling checknan')
            call wtd_check_nan(wtd)
         endif

      endif

   end subroutine wr_contact_pos

!------------------------------------------------------------------------------------------------------------

   subroutine wr_subsurf(wtd, idebug)
!--purpose: compute subsurface stresses for all cps in the w/r problem
      implicit none
!--subroutine arguments:
      type(t_ws_track)  :: wtd
      integer           :: idebug
!--local variables:
      integer           :: icp
!--pointers to global data items:
      type(t_probdata), pointer :: gd

      if (idebug.ge.2) call write_log(' --- Start subroutine wr_subsurf ---')

      ! loop over all contact patches, compute subsurface stresses

      do icp = 1, wtd%numcps

         associate(cp => wtd%allcps(icp)%cp)
         gd => cp%gd

         ! copy new input for subsurface stresses, except for cps that have gd%subs specified separately

         if (.not.cp%has_own_subs) call subsurf_copy(wtd%subs, gd%subs)

         ! compute subsurface stresses for contact patch, using mirroring for left rail/wheel combination

         call subsur_calc(gd%ic, gd%mater, gd%cgrid_cur, gd%outpt1%igs, gd%outpt1%ps, gd%ic%is_left_side(), &
                        gd%subs, idebug)

         ! write to subs-file if requested

         call subsur_matfil (gd%meta, gd%ic, gd%subs, idebug)

         end associate
      enddo

   end subroutine wr_subsurf

!------------------------------------------------------------------------------------------------------------

   subroutine wr_setup_cp(meta, ws, trk, numcps, icp, cp, wtd_ic, mater, discr, fric, kin, solv, x_locate)
!--purpose: define the contact problem for an initial contact point for a W/R contact case. 
      implicit none
!--subroutine arguments:
      type(t_metadata)  :: meta
      type(t_wheelset)  :: ws
      type(t_trackdata) :: trk
      integer           :: numcps, icp, x_locate
      type(t_cpatch)    :: cp
      type(t_ic)        :: wtd_ic
      type(t_material)  :: mater
      type(t_discret)   :: discr
      type(t_friclaw)   :: fric
      type(t_kincns)    :: kin
      type(t_solvers)   :: solv
!--local variables:
      character(len=5)          :: nam_side
      integer                   :: ip, ii, iy, mx, my, npot, ierror
      real(kind=8)              :: sw_ref, fac_warn, xy_ofs
      type(t_marker)            :: mref_pot, mref_rai, mref_whl, whl_trk
      type(t_potcon)            :: potcon_inp
      type(t_probdata), pointer :: gd

      if (x_locate.ge.3) then
         write(bufout,'(/a,i2,a)') ' --- Start subroutine wr_setup_cp for cp',icp,' ---'
         call write_log(2, bufout)
      endif

      ! set pointer to the active rail and wheel in the current configuration

      associate(my_rail     => trk%rai, my_wheel    => ws%whl)
      if (wtd_ic%is_left_side()) then
         nam_side = 'left'
      else
         nam_side = 'right'
      endif

      ! convert wheel-marker from wheel-set coords to track-coords

      whl_trk = marker_2glob( my_wheel%m_ws, ws%m_trk )

      ! convert contact reference from track-coords to rail coords

      mref_rai = marker_2loc( cp%mref, my_rail%m_trk )

      ! convert contact reference from track-coords to potcon-coords

      if (wtd_ic%use_supergrid()) then
         mref_pot = marker_2loc( cp%mref, cp%mpot )
      else
         call marker_init(mref_pot)
      endif

      ! convert contact reference from track-coords to wheel coords

      mref_whl = marker_2loc( cp%mref, whl_trk )

      ! get sw-position of mref on wheel profile

      call spline_get_s_at_y( my_wheel%prw%grd_data%spl, mref_whl%y(), sw_ref, ierror )

      if (x_locate.ge.3) then
         call marker_print( ws%m_trk,      'm_ws(trk)', 2 )
         call marker_print( my_wheel%m_ws, 'm_whl(ws)', 2 )
         call marker_print( whl_trk,       'm_whl(trk)', 2 )
         call marker_print( mref_whl,      'm_cp(whl)', 2 )
         call marker_print( mref_rai,      'm_cp(rai)', 2 )
         call marker_print( my_rail%m_trk, 'm_rai(trk)', 2 )
         call marker_print( cp%mref,       'm_cp(trk)', 2 )
         call marker_print( mref_pot,      'm_cp(pot)', 2 )
      endif

      ! allocate the hierarchical data-structure (TODO: already allocated?)

      if (.not.associated(cp%gd)) then
         if (wtd_ic%x_cpatch.ge.1) then
            write(bufout,'(a,i3)') ' wr_setup_cp: allocate gd for icp=',icp
            call write_log(1, bufout)
         endif
         allocate(cp%gd)
         call gd_init(cp%gd)
      endif
      gd => cp%gd

      !  - general metadata:
      gd%meta%reid      = meta%reid
      gd%meta%expnam    = meta%expnam
      gd%meta%wrkdir    = meta%wrkdir
      gd%meta%outdir    = meta%outdir
      gd%meta%npatch    = numcps
      gd%meta%ipatch    = icp

      !  - iteration-related metadata are set at the actual solution in wr_solve_cp

      ! copy wheel and rail position data to meta-data
      ! note: mirrorring for left wheel is ignored, internal data written to gd%meta

      gd%meta%tim       = meta%tim
      gd%meta%s_ws      = ws%s
      gd%meta%th_ws     = ws%pitch
      gd%meta%ynom_whl  = ws%whl%m_ws%y()
      gd%meta%rnom_whl  = ws%nom_radius
      gd%meta%rnom_rol  = 0d0
      if (wtd_ic%is_roller()) gd%meta%rnom_rol  = trk%nom_radius

      !  - position of wheel profile origin w.r.t. track origin
      gd%meta%x_w       = whl_trk%x()
      gd%meta%y_w       = whl_trk%y()
      gd%meta%z_w       = whl_trk%z()
      gd%meta%roll_w    = whl_trk%roll()
      gd%meta%yaw_w     = whl_trk%yaw()

      !  - position of rail profile origin w.r.t. track origin
      gd%meta%y_r       = my_rail%m_trk%y()
      gd%meta%z_r       = my_rail%m_trk%z()
      gd%meta%roll_r    = my_rail%m_trk%roll()

      !  - position of contact reference w.r.t. track origin
      gd%meta%xcp_tr    = cp%mref%x()
      gd%meta%ycp_tr    = cp%mref%y()
      gd%meta%zcp_tr    = cp%mref%z()
      gd%meta%deltcp_tr = cp%mref%roll()

      !  - position of contact reference w.r.t. rail origin
      gd%meta%xcp_r     = mref_rai%x()
      gd%meta%ycp_r     = mref_rai%y()
      gd%meta%zcp_r     = mref_rai%z()
      gd%meta%scp_r     = cp%sr_ref
      gd%meta%deltcp_r  = mref_rai%roll()

      !  - position of contact reference w.r.t. wheel origin
      gd%meta%xcp_w     = mref_whl%x()
      gd%meta%ycp_w     = mref_whl%y()
      gd%meta%zcp_w     = mref_whl%z()
      gd%meta%scp_w     = sw_ref
      gd%meta%deltcp_w  = mref_whl%roll()

      ! copy control digits, override some of them

      gd%ic             = wtd_ic
      gd%ic%norm        = 0
      gd%ic%force3      = 0
      gd%ic%discns3     = 1
      gd%ic%matfil_surf = -wtd_ic%matfil_surf

      ! copy solver settings

      call solv_copy( solv, gd%solv)

      ! copy material parameters

      call mater_copy( mater, gd%mater )

      ! copy kinematic parameters

      gd%kin%dt         = kin%dt
      gd%kin%facphi     = kin%facphi
      gd%kin%use_muscal = kin%use_muscal

      !  - set parameters describing the potential contact area as needed in undef.dist calculation

      gd%potcon_inp%ipotcn =  1
      gd%potcon_inp%dx     = cp%dx_eff
      gd%potcon_inp%mx     = nint((cp%xend-cp%xsta)/cp%dx_eff)
      gd%potcon_inp%xl     = cp%xsta - cp%mpot%x()

      if (gd%potcon_inp%mx.le.0 .or. gd%potcon_inp%mx.gt.discr%npot_max) then
         write(bufout,'(2(a,i0))') ' Internal error (setup_cp): mx = ',gd%potcon_inp%mx,', max = ',     &
                discr%npot_max
         call write_log(1, bufout)
         call abort_run()
      endif

      if (wtd_ic%is_conformal()) then
         ! conformal: curved surface
         gd%potcon_inp%dy  = cp%ds_eff
         gd%potcon_inp%yl  = cp%sc_sta
         gd%potcon_inp%my  = nint((cp%sc_end-cp%sc_sta)/cp%ds_eff)
         fac_warn      = 1.2d0          ! accomodate for fac_sc=1.2 in compute_curved_potcon
      else
         ! planar: tangent plane
         gd%potcon_inp%dy  = cp%ds_eff
         gd%potcon_inp%yl  = cp%sp_sta
         gd%potcon_inp%my  = nint((cp%sp_end-cp%sp_sta)/cp%ds_eff)
         fac_warn      = 1.0d0
      endif

      ! complete remaining entries in potcon_inp including npot

      call potcon_fill( gd%potcon_inp )

      ! update potcon_inp to cover the contact area of the previous time
 
      if (wtd_ic%pvtime.ne.2 .and. .not.gd%is_new) then
         potcon_inp = cp%gd%potcon_inp
         call merge_prev_potcon(icp, cp%gd, wtd_ic%x_cpatch)

         ! for conformal contact, extend curv_ref, curv_nrm, curv_incln to new potential contact size
         ! (arrays from previous time deleted in wr_update_allcps)

         if (wtd_ic%is_conformal()) then
            call extend_curv_ref(wtd_ic, potcon_inp, cp%gd%potcon_inp, cp%curv_ref, cp%curv_nrm,        &
                cp%curv_incln)
         endif
      endif

      ! check that npot_max will not be exceeded

      npot = gd%potcon_inp%npot
      if (wtd_ic%return.le.1 .and. .not.wtd_ic%use_supergrid() .and. discr%npot_max.ge.100 .and.        &
          npot.ge.fac_warn*discr%npot_max) then
         write(bufout,'(2(a,i6))') ' Internal error: npot=',npot,' >= NPOT_MAX=', discr%npot_max
         call write_log(1, bufout)
         write(bufout,'(a,2i5,a,i6,a,2f7.3,a,2f10.3)') ' mx,my=', gd%potcon_inp%mx, gd%potcon_inp%my,   &
             ', npot=',npot, ', dx,dy=', gd%potcon_inp%dx, gd%potcon_inp%dy, ', xl,yl=',                &
             gd%potcon_inp%xl, gd%potcon_inp%yl
         call write_log(1, bufout)
         call abort_run()
      endif

      !  - create the main grid "cgrid"

      call potcon_cgrid(gd%potcon_inp, gd%cgrid_inp)

      if (x_locate.ge.2) then
         mx   = gd%potcon_inp%mx
         npot = gd%potcon_inp%npot
         write(bufout,'(2(a,2f10.4))') ' x_cntc =', gd%cgrid_inp%x(1), gd%cgrid_inp%x(2),'...',         &
                        gd%cgrid_inp%x(npot-1), gd%cgrid_inp%x(npot)
         call write_log(1, bufout)
         write(bufout,'(2(a,2f10.4))') ' y_cntc =', gd%cgrid_inp%y(mx), gd%cgrid_inp%y(2*mx),'...',     &
                        gd%cgrid_inp%y(npot-mx), gd%cgrid_inp%y(npot)
         call write_log(1, bufout)
      endif

      if (wtd_ic%x_cpatch.ge.2) then
         associate(g1 => gd%cgrid_inp, g2 => gd%cgrid_cur)
         write(bufout,'(2(a,i4),2(a,f8.3))') ' cgrid_inp: mx,my=',g1%nx,',',g1%ny,', dx,dy=', g1%dx,',',g1%dy
         call write_log(1, bufout)
         write(bufout,'(2(a,i4),2(a,f8.3))') ' cgrid_cur: mx,my=',g2%nx,',',g2%ny,', dx,dy=', g2%dx,',',g2%dy
         call write_log(1, bufout)
         end associate
      endif

      ! set P-digit for this gd
      !  - zero previous tractions (P=2) for new contact patches
      !  - leave previous tractions untouched (P=3) if gd has been used already in case ncase

      if (gd%is_new .or. wtd_ic%pvtime.eq.2) then
         gd%ic%pvtime = 2
      elseif (gd%meta%ncase.ge.meta%ncase) then
         gd%ic%pvtime = 3
      else
         gd%ic%pvtime = wtd_ic%pvtime
      endif

      ! set I-digit for this gd
      !  - no initial estimate for new contact patches
      !  - no initial estimate if different grid sizes are used (dx, dy)
      !  - no initial estimate in case of reduced #elements (npot_max)
      !  - no initial estimate in case of large change in potential contact area (npot)
      !  - user setting I from wtd when a new case is started
      !  - use initial estimate if this gd was used in previous iteration
      !  - no initial estimate if this gd was inactive in previous iteration

      if (gd%is_new) then
         gd%ic%iestim = 0
         if (wtd_ic%x_cpatch.ge.1) call write_log('    no initial estimate (new patch)...')
      elseif (.not.equal_grid_sizes(gd%potcon_inp%dx, gd%potcon_cur%dx,                                 &
                                    gd%potcon_inp%dy, gd%potcon_cur%dy)) then
         gd%ic%iestim = 0
         if (wtd_ic%x_cpatch.ge.1) call write_log('    no initial estimate (change dx or dy)...')
      elseif (abs(gd%potcon_inp%npot - gd%potcon_cur%npot).gt.                                          &
                0.4d0 * min(gd%potcon_inp%npot, gd%potcon_cur%npot)) then
         gd%ic%iestim = 0
         if (wtd_ic%x_cpatch.ge.1) call write_log('    no initial estimate (large change npot)...')
      elseif (gd%meta%ncase.lt.meta%ncase) then
         gd%ic%iestim = wtd_ic%iestim
         if (gd%ic%iestim.ge.1) then
            if (wtd_ic%x_cpatch.ge.1) call write_log(' using initial estimate (new case, user setting)...')
         else
            if (wtd_ic%x_cpatch.ge.1) call write_log('    no initial estimate (new case, user setting)...')
         endif
      elseif (gd%meta%itforce.eq.meta%itforce-1) then
         gd%ic%iestim = 1
         if (wtd_ic%x_cpatch.ge.1) call write_log(' using initial estimate (prev.iteration)...')
      else
         gd%ic%iestim = 0
         if (wtd_ic%x_cpatch.ge.1) call write_log('    no initial estimate (cp skipped iteration)...')
      endif

      !  - fill in surface inclinations for conformal infl.cf (Blanco approach)

      if (mater%gencr_eff.eq.4) then
         mx   = gd%potcon_inp%mx
         my   = gd%potcon_inp%my
         gd%mater%ninclin = my

         if (allocated(gd%mater%surf_inclin)) deallocate(gd%mater%surf_inclin)
         allocate(gd%mater%surf_inclin(my,2))
         if (x_locate.ge.3) call write_log('%   iy     si      ai');
         do iy = 1, my
            gd%mater%surf_inclin(iy,1) = gd%cgrid_inp%y(iy*mx)
            gd%mater%surf_inclin(iy,2) = cp%curv_incln%vy(iy)
            if (x_locate.ge.3) then
               write(bufout,'(i6,f8.3,f11.6)') iy, gd%mater%surf_inclin(iy,1), gd%mater%surf_inclin(iy,2)
               call write_log(1, bufout)
            endif
         enddo
      endif

      ! copy/interpolate friction law data

      if (wtd_ic%varfrc.eq.0) then
         gd%ic%varfrc = 0
         call fric_copy(fric, gd%fric)
      elseif (wtd_ic%varfrc.eq.1 .and. wtd_ic%is_conformal()) then
         call fric_interp(fric, gd%potcon_inp%my, cp%curv_incln%vy, gd%fric)
         gd%ic%varfrc = 3
      elseif (wtd_ic%varfrc.eq.1) then
         gd%ic%varfrc = 0
         call fric_interp(fric, gd%meta%deltcp_w, gd%fric)
      elseif (wtd_ic%varfrc.eq.2) then
         gd%ic%varfrc = 0
         call fric_interp(fric, gd%meta%s_ws, gd%fric)
      else
         write(bufout,'(a,i0,a)') ' Internal error: incorrect V= ',wtd_ic%varfrc,', aborting.'
         call write_log(1, bufout)
         call abort_run()
      endif
      gd%fric%varfrc_eff = gd%ic%varfrc

      if (x_locate.ge.3) then
         write(bufout,'(2(a,i0),a)') ' using V = ', gd%ic%varfrc,' with NVF = ', gd%fric%nvf,' slices'
         call write_log(1, bufout)
      endif
      if (x_locate.ge.4) then
         call fric_print(fric,'wtd%fric')
         call fric_print(gd%fric,'gd%fric')
      endif

      ! constant friction parameters per patch: enable scaling output forces
      ! Note: potential conflicts for F=1 with Fx<>0.

      gd%kin%use_muscal = (gd%ic%varfrc.eq.0 .or. gd%ic%varfrc.eq.2)
      if (gd%kin%use_muscal) then
         gd%kin%muscal = gd%fric%fstat()
      else
         gd%kin%muscal = 1d0
      endif

      ! compute the undeformed distance function

      gd%ic%rznorm  = 1
      gd%geom%ibase = 9
      call reallocate_arr(gd%geom%prmudf, gd%potcon_inp%npot)

      gd%kin%pen     =  0d0

      ! perform actual computation when R=0 or 1

      if (wtd_ic%return.le.1) then

         call wr_ud_planar (meta, wtd_ic, ws, trk, icp, cp, gd, x_locate)

      endif

      ! set planform data for reduced interaction between sub-patches

      if (cp%nsub.le.1) then
         gd%geom%iplan = 1
      else
         gd%geom%iplan = 4
         associate( np => gd%geom%npatch )
         np = cp%nsub

         if (x_locate.ge.3) then
            write(bufout,'(a,i3,a)') ' npatch=',np,' sub-patches'
            call write_log(1, bufout)
            do ip = 1, np
               write(bufout,'(a,i3,6(a,f8.3),a)') ' ip=',ip,': x_tr=[',cp%xyzlim(ip,1),',',             &
                        cp%xyzlim(ip,2),'], y_tr=[', cp%xyzlim(ip,3),',',cp%xyzlim(ip,4),'], z_tr=[',   &
                        cp%xyzlim(ip,5),',',cp%xyzlim(ip,6),']'
               call write_log(1, bufout)
            enddo
            call marker_print(cp%mpot, 'mpot(trk)', 2)
         endif

         call reallocate_arr(gd%geom%xylim, np, 4)
         call reallocate_arr(gd%geom%facsep, np, np)
         gd%geom%xylim(1:np,1) = -999d0
         gd%geom%xylim(1:np,2) =  999d0
         gd%geom%xylim(1:np,3) = -999d0
         gd%geom%xylim(1:np,4) =  999d0
         gd%geom%facsep(1:np, 1:np) = 0d0

         ! copy track xsta/xend positions of sub-patches +/- safety

         xy_ofs = 0.4d0 * discr%dist_comb
         gd%geom%xylim(1:np, 1) = cp%xyzlim(1:np, 1) - cp%mpot%x() - xy_ofs
         gd%geom%xylim(1:np, 2) = cp%xyzlim(1:np, 2) - cp%mpot%x() + xy_ofs

         ! convert track ysta/yend position in rail profile to lateral s position on contact grid

         if (.not.wtd_ic%is_conformal()) then
            gd%geom%xylim(1:np, 3:4) = (cp%xyzlim(1:np, 3:4) - cp%mref%y()) / cos(cp%mref%roll()) +     &
                                                                                 cp%sr_ref - cp%sr_pot
         else
            call compute_sc_for_yz(cp, np, cp%xyzlim(1:np,3), cp%xyzlim(1:np,5), wtd_ic%x_locate,       &
                        gd%geom%xylim(1:np,3))
            call compute_sc_for_yz(cp, np, cp%xyzlim(1:np,4), cp%xyzlim(1:np,6), wtd_ic%x_locate,       &
                        gd%geom%xylim(1:np,4))
         endif

         ! add safety in s-direction

         gd%geom%xylim(1:np, 3) = gd%geom%xylim(1:np, 3) - xy_ofs
         gd%geom%xylim(1:np, 4) = gd%geom%xylim(1:np, 4) + xy_ofs

         ! copy reduction factors

         gd%geom%facsep(1:np, 1:np) = cp%f_sep2(1:np, 1:np)

         if (x_locate.ge.2) then
            write(bufout,'(a,i3)') ' IPLAN=4: npatch=',np
            call write_log(1, bufout)

            do ip = 1, np
               write(bufout,'(4(a,f8.3),a,6f7.3)') ' x_c=[',gd%geom%xylim(ip,1),',',gd%geom%xylim(ip,2), &
                        '], s_p=[', gd%geom%xylim(ip,3),',',gd%geom%xylim(ip,4),'], fac=',            &
                        (gd%geom%facsep(ip,ii), ii=1, np)
               call write_log(1, bufout)
            enddo
         endif

         end associate
      endif

      ! compute velocity and creepages

      if (wtd_ic%tang.eq.0 .or. wtd_ic%return.ge.2) then

         ! T = 0: skip tangential contact calculation
         ! R = 2 or 3: skip actual computation

         gd%ic%tang     =  0
         gd%kin%cksi    =  0d0
         gd%kin%ceta    =  0d0
         gd%kin%cphi    =  0d0
         gd%kin%spinxo  = mref_pot%x()
         gd%kin%spinyo  = mref_pot%y()

      else

         ! set spin center using offset reference marker -- pot.contact

         gd%kin%spinxo = mref_pot%x()
         gd%kin%spinyo = mref_pot%y()

         ! compute the creepages or the rigid slip function w.r.t. reference marker

         call wr_rigid_slip(wtd_ic, ws, trk, discr%dqrel, icp, cp, gd, x_locate)

      endif
 
      ! subsurface stress computation is performed separately, if needed

      gd%ic%stress = 0

      if (x_locate.ge.4) call write_log('--- end subroutine wr_setup_cp ---')
      end associate

   end subroutine wr_setup_cp

!------------------------------------------------------------------------------------------------------------

   subroutine compute_sc_for_yz(cp, ny, ylim, zlim, x_locate, sclim)
!--purpose: convert limits [yz] for sub-patches to corresponding sc on curved reference
      implicit none
!--subroutine arguments:
      type(t_cpatch)  :: cp
      integer         :: ny, x_locate
      real(kind=8)    :: ylim(ny), zlim(ny), sclim(ny)
!--local variables:
      integer         :: iy, j, jmin
      real(kind=8)    :: dstmin, dst

      ! using brute-force double loop computing ny x ntot distances
      ! TODO: optimize for large ny and large ntot

      do iy = 1, ny
         
         ! determine closest grid point on curved reference

         jmin   = 0
         dstmin = 1d9
         do j = 1, cp%curv_ref%ntot
            ! dst = abs(ylim(iy) - cp%curv_ref%y(j)) + abs(zlim(iy) - cp%curv_ref%z(j))
            dst = (ylim(iy) - cp%curv_ref%y(j))**2 + (zlim(iy) - cp%curv_ref%z(j))**2
            if (j.eq.1 .or. dst.lt.dstmin) then
               jmin   = j
               dstmin = dst
            endif
         enddo

         sclim(iy) = cp%sc_sta + (jmin-1) * cp%ds_eff

         if (x_locate.ge.3) then
            write(bufout,'(a,i3,2(a,f8.3),a,i3,3(a,f8.3),a)') ' sc_for_yz: i=',iy,' (',ylim(iy),',',    &
                zlim(iy), '): closest j=',jmin,', sc=',sclim(iy),', (',cp%curv_ref%y(jmin),',',         &
                cp%curv_ref%z(jmin),')'
            call write_log(1, bufout)
         endif
      enddo

   end subroutine compute_sc_for_yz

!------------------------------------------------------------------------------------------------------------

   subroutine wr_solve_cp(meta, ws, trk, icp, cp, ic, idebug, ierror)
!--purpose: solve the contact problem for an initial contact point for a W/R contact case. 
      implicit none
!--subroutine arguments:
      type(t_metadata)  :: meta
      type(t_wheelset)  :: ws
      type(t_trackdata) :: trk
      integer           :: icp, idebug
      type(t_cpatch)    :: cp
      type(t_ic)        :: ic
      integer           :: ierror
!--local variables:
      integer                   :: nadh, nslip, nplast, nexter, ncon, iatbnd
      type(t_probdata), pointer :: gd

      if (idebug.ge.3) then
         write(bufout,'(a,i2,a)') ' --- Start subroutine wr_solve_cp for cp',icp,' ---'
         call write_log(1, bufout)
      endif

      gd => cp%gd

      ! set iteration-related metadata -- 'last calculated time+iteration'

      gd%meta%ncase   = meta%ncase
      gd%meta%itforce = meta%itforce

      ! solve the contact problem

      gd%ic%output_surf = 0
      if (idebug.ge.5) call write_log(' ... calling contac')
      call contac(gd, ierror)

      ! Count the number of elements in contact

      call eldiv_count(gd%outpt1%igs, nadh, nslip, nplast, nexter)
      ncon = nadh + nslip + nplast

      if (ncon.le.0 .and. ic%ilvout.ge.1 .and. idebug.ge.-1) then
         write(bufout,'(2(a,i0))') ' WARNING: no actual contact in patch ',icp,': ncon = ',ncon
         call write_log(1, bufout)

         if (idebug.ge.3) call wrigs(gd%outpt1%igs, .true., 0d0)
      endif

      ! Count the number of interior elements at the boundaries of the potential contact

      iatbnd = eldiv_count_atbnd(gd%outpt1%igs, 0, idebug)
      if (iatbnd.ge.1 .and. ic%ilvout.ge.1 .and. idebug.ge.-1) then
         write(bufout,'(a,i5,a)') ' WARNING: potential contact area too small;',iatbnd,                 &
                ' elements next to boundary.'
         call write_log(1, bufout)

         if (idebug.ge.3) call wrigs(gd%outpt1%igs, .true., 0d0)
         ! ierror = -17
      endif

      ! conformal: rotate tractions, compute overall forces

      if (ic%is_conformal()) then
         call aggregate_forces(cp, gd%ic, gd%cgrid_cur, gd%kin, gd%mater, gd%outpt1, idebug)

      endif

      ! rotate forces from contact-reference coordinates to global coordinates

      call cpatch_forces_moments(ws, trk, cp, idebug)

      if (idebug.ge.4) write(*,*) '--- end subroutine wr_solve_cp ---'

   end subroutine wr_solve_cp

!------------------------------------------------------------------------------------------------------------

   subroutine aggregate_forces(cp, ic, cgrid, kin, mater, outpt1, idebug)
!--purpose: aggregate total forces from curved ref. surface to planar contact contact coordinates
      implicit none
!--subroutine arguments:
      type(t_cpatch)   :: cp
      type(t_ic)       :: ic
      type(t_grid)     :: cgrid
      type(t_kincns)   :: kin
      type(t_material) :: mater
      type(t_output)   :: outpt1
      integer          :: idebug
!--local variables:
      type(t_vec)      :: nref
      type(t_grid)     :: curv_cp
      type(t_gridfnc3) :: ps_cp
      real(kind=8)     :: dxdy, bar_fn, bar_fx, bar_fs

      associate( npot   => cgrid%ntot,    fcntc  => kin%fcntc,                                          &
                 mxtrue => outpt1%mxtrue, mytrue => outpt1%mytrue, mztrue => outpt1%mztrue )

      ! create 3-d version of curved reference surface

      call grid_extrude_profile(cp%curv_ref, cgrid%nx, cgrid%x(1)+cp%mref%x(), cgrid%dx, curv_cp)

      ! transform from global track coordinates to local planar sp coordinates

      call cartgrid_2loc(curv_cp, cp%mref)

      ! copy tractions on curved sc reference to work variable, rotate to planar sp coordinates

      if (idebug.ge.2) call grid_print(cp%curv_ref, 'curv_ref', 5)
      if (idebug.ge.4) call grid_print(curv_cp, 'curv_cp', 5)

      nref = -cp%mref%kvec()
      if (idebug.ge.3) then
         write(bufout,'(3(a,f8.3),a)') ' ref.normal=[',nref%x(), ',', nref%y(), ',', nref%z(), ']'
         call write_log(1, bufout)
      endif

      call gf3_copy_struc(outpt1%ps, ps_cp, 'ps(cp)')
      call gf3_copy(AllElm, outpt1%ps, ps_cp, ikALL)

      if (idebug.ge.4) call gf3_print(ps_cp, 'ps(curved)', ikZDIR, 4)
      if (idebug.ge.4) call gf3_print(cp%curv_nrm, 'curv_nrm', ikALL, 5)

      call gf3_curv2cart(ps_cp, nref, cp%curv_nrm, 1, idebug)

      if (idebug.ge.4) call gf3_print(ps_cp, 'ps(planar)', ikZDIR, 4)

      ! compute aggregate forces w.r.t. planar contact reference marker

      dxdy = cgrid%dxdy
      bar_fn = dxdy * gf3_sum(AllElm, ps_cp, ikZDIR)
      bar_fx = dxdy * gf3_sum(AllElm, ps_cp, ikXDIR)
      bar_fs = dxdy * gf3_sum(AllElm, ps_cp, ikYDIR)

      ! set resultant contact force at current time; note: leaving fntrue/fxrel/fyrel untouched (tilde(f))

      fcntc  = (/ bar_fx, bar_fs, bar_fn /)

      ! compute ad-hoc proportional damping using fprev, fcntc

      call calc_damping_force( ic, mater, kin )

      ! Compute torsional moments about planar contact x, sp and np-axes
      ! using (x,sp,np)-coordinates provided in curv_cp

      mxtrue = dxdy * ddot(npot, curv_cp%y, 1, ps_cp%vn, 1) - dxdy * ddot(npot, curv_cp%z, 1, ps_cp%vy, 1)
      mytrue = dxdy * ddot(npot, curv_cp%z, 1, ps_cp%vx, 1) - dxdy * ddot(npot, curv_cp%x, 1, ps_cp%vn, 1)
      mztrue = dxdy * ddot(npot, curv_cp%x, 1, ps_cp%vy, 1) - dxdy * ddot(npot, curv_cp%y, 1, ps_cp%vx, 1)

      if (idebug.ge.2) then
         write(bufout,'(3(a,g14.7))') ' Total moments Mx=',mxtrue,', Msp=',mytrue,', Mnp=',mztrue
         call write_log(1, bufout)
      endif

      end associate

      call grid_destroy(curv_cp)
      call gf3_destroy(ps_cp)

   end subroutine aggregate_forces

!------------------------------------------------------------------------------------------------------------

   subroutine cpatch_forces_moments(ws, trk, cp, idebug)
!--purpose: rotate and shift total forces and moments to track and wheelset systems
      implicit none
!--subroutine arguments:
      type(t_wheelset)  :: ws
      type(t_trackdata) :: trk
      integer           :: idebug
      type(t_cpatch)    :: cp
!--local variables:
      type(t_vec)       :: fcntc, fdamp, tcntc
      type(t_marker)    :: mwhl_trk

      associate(my_rail => trk%rai, my_wheel => ws%whl, gd => cp%gd)

      ! compute derived quantities for use in total-force iteration, based on tilde(f)

      cp%c_hz   = gd%kin%fntrue / max(1d-9, gd%kin%pen)**(1.5)
      cp%fs_rel = gd%kin%fyrel * gd%kin%muscal

      ! rotate forces from contact-reference coordinates to global coordinates

      fcntc = vec( gd%kin%fcntc )
      fdamp = vec( gd%kin%fdamp )
      tcntc = vec( gd%outpt1%mxtrue, gd%outpt1%mytrue, gd%outpt1%mztrue )

      if (idebug.ge.2) then
         call marker_print(cp%mref, 'cp(trk)', 2)
         write(bufout,800) ' rotating F_(cp)=    [',fcntc%x(),',',fcntc%y(),',',fcntc%z(),'] by delttr=', &
                cp%delttr
         call write_log(1, bufout)
         write(bufout,800) ' rotating M_@cp(cp)= [',tcntc%x(),',',tcntc%y(),',',tcntc%z(),'] by delttr=', &
                cp%delttr
         call write_log(1, bufout)
 800     format(4(a,g14.6))
      endif

      ! total force = contact force + damping

      cp%ftrk = cp%mref%rot * (fcntc + fdamp)
      tcntc   = cp%mref%rot * tcntc             ! T_@cp(cp) --> T_@cp(tr)

      ! compute moment T_@rr(tr) = (m_cp(tr) - m_r(tr)) x F_(tr)  +  T_@cp(tr)

      cp%ttrk = ((cp%mref%o - my_rail%m_trk%o) .cross. cp%ftrk) + tcntc

      if (idebug.ge.2) then
         call marker_print(my_rail%m_trk, 'r(trk)', 2)
         write(bufout,800) ' result:  F_(tr)=    [',cp%ftrk%x(),',',cp%ftrk%y(),',',cp%ftrk%z(),'].'
         call write_log(1, bufout)
         write(bufout,800) '          M_@r(tr)=  [',cp%ttrk%x(),',',cp%ttrk%y(),',',cp%ttrk%z(),'].'
         call write_log(1, bufout)
      endif

      ! rotate forces (on rail) from global coordinates to wheelset x,y,z-coordinates

      cp%fws  = ws%m_trk .transp. cp%ftrk
      ! tcntc   = ws%m_trk .transp. tcntc         ! T_@cp(tr) --> T_@cp(ws)

      ! compute moment T_@w(ws) = R_w(tr)^T * { (m_w(tr) - m_cp(tr)) x F_(tr)  +  T_@cp(tr) }

      mwhl_trk = marker_2glob( my_wheel%m_ws, ws%m_trk )

      cp%tws  = ws%m_trk .transp. (((cp%mref%o - mwhl_trk%o ) .cross. cp%ftrk) + tcntc)

      if (idebug.ge.2) then
         call marker_print(ws%m_trk, 'ws(trk)', 2)
         call marker_print(mwhl_trk, 'whl(trk)', 2)
         write(bufout,800) '          F_(ws)=    [',cp%fws%x(),',',cp%fws%y(),',',cp%fws%z(),'].'
         call write_log(1, bufout)
         write(bufout,800) '          M_@w(ws)=  [',cp%tws%x(),',',cp%tws%y(),',',cp%tws%z(),'].'
         call write_log(1, bufout)
      endif
      end associate

   end subroutine cpatch_forces_moments

!------------------------------------------------------------------------------------------------------------

   subroutine total_forces_moments(wtd)
!--purpose: compute total forces and moments
   implicit none
!--subroutine arguments:
      type(t_ws_track)          :: wtd
!--local variables:
      integer                   :: icp

      ! accumulate total forces and moments F_(trk), M_@r(trk), F_(ws) and M_@w(ws)

      wtd%ftrk = vec( 0d0, 0d0, 0d0 )
      wtd%ttrk = vec( 0d0, 0d0, 0d0 )
      wtd%fws  = vec( 0d0, 0d0, 0d0 )
      wtd%tws  = vec( 0d0, 0d0, 0d0 )

      do icp = 1, wtd%numcps
         associate(cp => wtd%allcps(icp)%cp)
         wtd%ftrk = wtd%ftrk + cp%ftrk
         wtd%ttrk = wtd%ttrk + cp%ttrk
         wtd%fws  = wtd%fws  + cp%fws 
         wtd%tws  = wtd%tws  + cp%tws 
         end associate
      enddo

   end subroutine total_forces_moments

!------------------------------------------------------------------------------------------------------------

   subroutine set_dummy_solution(cp)
!--purpose: initialize contact patch with no contact
   implicit none
!--subroutine arguments:
      type(t_cpatch)               :: cp
!--local variables:

      ! initialize the contact patch data-structure: micp = mref = 0

      call cp_init(cp)
      cp%delttr  = 0d0
      cp%gap_min = 1d0
      cp%ftrk    = vec_zero()
      cp%ttrk    = vec_zero()
      cp%fws     = vec_zero()
      cp%tws     = vec_zero()

      ! initialize gd data-structure: 1-element potcon, zero solution

      allocate(cp%gd)
      call gd_init(cp%gd)

   end subroutine set_dummy_solution

!------------------------------------------------------------------------------------------------------------

   subroutine wtd_check_nan(wtd)
!--purpose: check solution of a w/r contact problem
      implicit none
!--subroutine arguments:
      type(t_ws_track)  :: wtd
!--local variables:
      integer           :: icp

      if (any(isnan(wtd%ftrk%v)) .or. any(isnan(wtd%ttrk%v))) then

         ! write error-message

         write(bufout,'(a,i8,4(a,i4),a)') ' check_nan: found NaN-values for case', wtd%meta%ncase,      &
                ' on result element', wtd%meta%reid,', it', wtd%meta%itforce,' (',wtd%meta%itforc_out,  &
                  '/', wtd%meta%itforc_inn,')'
         call write_log(1, bufout)

         ! write out-file for the case

         call wr_output(wtd)

         ! write mat-files for each contact patch

         do icp = 1, wtd%numcps
            associate( gd => wtd%allcps(icp)%cp%gd )
            gd%ic%matfil_surf = 2
            call writmt (gd%meta, gd%ic, gd%cgrid_cur, gd%potcon_cur, gd%mater, gd%fric, gd%kin,        &
                gd%geom, gd%outpt1, wtd%ic%is_left_side())
            end associate
         enddo
      endif
         
   end subroutine wtd_check_nan

!------------------------------------------------------------------------------------------------------------

   subroutine merge_prev_potcon(icp, gd, idebug)
!--purpose: extend gd%potcon_inp to cover previous contact area
      implicit none
!--subroutine arguments:
      integer           :: icp, idebug
      type(t_probdata)  :: gd
!--local variables:
      integer           :: ierror
      type(t_potcon)    :: potcon_prv, potcon_tot

      ! determine potential contact area tight around contact area for previous time
 
      associate(igv => gd%outpt1%igv)
      call areas(igv)

      potcon_prv = gd%potcon_cur
      potcon_prv%ipotcn = 1

      if (igv%ixmax.ge.igv%ixmin) then
         potcon_prv%mx = igv%ixmax - igv%ixmin + 1
         potcon_prv%xl = potcon_prv%xl + (igv%ixmin-1) * potcon_prv%dx
      else
         potcon_prv%mx = 1
      endif

      if (igv%iymax.ge.igv%iymin) then
         potcon_prv%my = igv%iymax - igv%iymin + 1
         potcon_prv%yl = potcon_prv%yl + (igv%iymin-1) * potcon_prv%dy
      else
         potcon_prv%my = 1
      endif

      call potcon_fill( potcon_prv )
      end associate

      ! merge proposed potcon_inp with potcon_prv for previous time

      call potcon_merge( potcon_prv, gd%potcon_inp, potcon_tot, idebug, ierror )
      if (ierror.ne.0) then
         call write_log(' Internal error(merge_prev_potcon): previous/current grids cannot be merged.')
         call abort_run()
      endif

      if (idebug.ge.1) then
         write(bufout,'(a,i3,2(2(a,f8.3),a,i4))') '   prev',icp,': x=[',gd%potcon_cur%xl,',',           &
                gd%potcon_cur%xh,'], mx=',gd%potcon_cur%mx,', y=[',gd%potcon_cur%yl,',',                &
                gd%potcon_cur%yh,'], my=',gd%potcon_cur%my
         call write_log(1, bufout)
         write(bufout,'(a,i3,2(2(a,f8.3),a,i4))') '  tight',icp,': x=[',potcon_prv%xl,',',              &
                potcon_prv%xh,'], mx=',potcon_prv%mx,', y=[',potcon_prv%yl,',',                         &
                potcon_prv%yh,'], my=',potcon_prv%my
         call write_log(1, bufout)
         write(bufout,'(a,i3,2(2(a,f8.3),a,i4))') '    new',icp,': x=[',gd%potcon_inp%xl,',',           &
                gd%potcon_inp%xh,'], mx=',gd%potcon_inp%mx,', y=[',gd%potcon_inp%yl,',',                &
                gd%potcon_inp%yh,'], my=',gd%potcon_inp%my
         call write_log(1, bufout)
         write(bufout,'(a,i3,2(2(a,f8.3),a,i4))') ' merged',icp,': x=[',potcon_tot%xl,',',              &
                potcon_tot%xh,'], mx=',potcon_tot%mx,', y=[',potcon_tot%yl,',',                         &
                potcon_tot%yh,'], my=',potcon_tot%my
         call write_log(1, bufout)
      endif

      ! store extended potcon for new time instance

      gd%potcon_inp = potcon_tot

   end subroutine merge_prev_potcon

!------------------------------------------------------------------------------------------------------------

   subroutine extend_curv_ref(ic, pot_old, pot_new, curv_ref, curv_nrm, curv_incln)
!--purpose: extend curved reference grid according to new dimensions of potential contact area
      implicit none
!--subroutine arguments:
      type(t_ic)                :: ic
      type(t_potcon)            :: pot_old, pot_new        ! size of old/new curv_ref grids
      type(t_grid),    target   :: curv_ref
      type(t_gridfnc3)          :: curv_nrm, curv_incln
!--local variables:
      logical            :: is_ok, is_equal
      integer            :: ix0, ix1, iy0, iy1, kofs_x, kofs_y, iy_old, iy_new
      type(t_grid)       :: curv_new
      real(kind=8)       :: dcoor(3)

      ! check that the two grids have matching (dx,dy) and determine offset ky
      ! ky == #extra rows at start of pot_old wrt pot_new

      call potcon_get_overlap(pot_old, pot_new, kofs_x, kofs_y, ix0, ix1, iy0, iy1, is_ok, is_equal)

      if (ic%x_cpatch.ge.2) then
         write(bufout,'(2(a,f8.3,a,i4))')     ' pot_old: yl=',pot_old%yl,', my=',pot_old%my,', yh=',    &
                pot_old%yh,', offset k_y=',kofs_y
         call write_log(1, bufout)
         write(bufout,'(a,f8.3,a,i4,a,f8.3)') ' pot_new: yl=',pot_new%yl,', my=',pot_new%my,', yh=',    &
                pot_new%yh
         call write_log(1, bufout)
      endif

      ! allocate space for new grid and grid-functions

      call grid_create_curvil(curv_new, 1, pot_new%my)

      ! copy data from curv_ref to curv_new with linear extrapolation

      do iy_old = 1, pot_old%my
         iy_new = iy_old - kofs_y
         curv_new%coor(iy_new,1:3) = curv_ref%coor(iy_old,1:3)
      enddo

      ! linear extrapolation at start of old curv_ref

      dcoor(1:3) = curv_ref%coor(1,1:3) - curv_ref%coor(2,1:3)
      do iy_new = -kofs_y, 1, -1
         curv_new%coor(iy_new,1:3) = curv_new%coor(iy_new+1,1:3) + dcoor(1:3)
      enddo

      ! linear extrapolation at end of old curv_ref

      dcoor(1:3) = curv_ref%coor(pot_old%my,1:3) - curv_ref%coor(pot_old%my-1,1:3)
      do iy_new = pot_old%my-kofs_y+1, pot_new%my
         curv_new%coor(iy_new,1:3) = curv_new%coor(iy_new-1,1:3) + dcoor(1:3)
      enddo

      if (ic%x_cpatch.ge.4) then
         do iy_new = 1, pot_new%my
            iy_old = iy_new + kofs_y
            if (iy_old.ge.1 .and. iy_old.le.pot_old%my) then
               write(bufout,'(a,i3,a,3f12.3,a,i3)') ' iyn=',iy_new,': curv_new=[',curv_new%x(iy_new),   &
                   curv_new%y(iy_new), curv_new%z(iy_new),'] <-- iyo=',iy_old
               call write_log(1, bufout)
            else
               write(bufout,'(a,i3,a,3f12.3,a)') ' iyn=',iy_new,': curv_new=[',curv_new%x(iy_new),      &
                   curv_new%y(iy_new), curv_new%z(iy_new),']'
               call write_log(1, bufout)
            endif
         enddo
      endif

      ! move curv_new into subroutine argument curv_ref

      call grid_destroy(curv_ref)
      call grid_copy(curv_new, curv_ref)
      call grid_destroy(curv_new)

      ! extend grid functions curv_nrm and curv_incln

      call gf3_resize_curv(curv_nrm, curv_new, 0, kofs_y, 1, 1, iy0, iy1, 0d0)
      call gf3_resize_curv(curv_incln, curv_new, 0, kofs_y, 1, 1, iy0, iy1, 0d0)

      ! constant extrapolation at start of old curv_nrm/incln

      do iy_new = 1, -kofs_y
         curv_nrm%val(iy_new,1:3)   = curv_nrm%val(1-kofs_y,1:3)
         curv_incln%val(iy_new,1:3) = curv_incln%val(1-kofs_y,1:3)
      enddo

      ! constant extrapolation at end of old curv_nrm/incln

      do iy_new = pot_old%my-kofs_y+1, pot_new%my
         curv_nrm%val(iy_new,1:3)   = curv_nrm%val(pot_old%my-kofs_y,1:3)
         curv_incln%val(iy_new,1:3) = curv_incln%val(pot_old%my-kofs_y,1:3)
      enddo

      if (ic%x_cpatch.ge.4) call gf3_print(curv_nrm,   'curv_nrm',   ikALL, 5)
      if (ic%x_cpatch.ge.4) call gf3_print(curv_incln, 'curv_incln', ikALL, 5)

      curv_nrm%grid   => curv_ref
      curv_incln%grid => curv_ref

   end subroutine extend_curv_ref

!------------------------------------------------------------------------------------------------------------

end module m_wr_solvecp
