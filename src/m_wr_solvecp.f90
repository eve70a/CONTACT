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
public  wr_solve_cp
private aggregate_forces
private cpatch_forces_moments
public  total_forces_moments

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
            write(bufout,'(a,f9.4)') ' z_ws=',wtd%ws%z
            call write_log(1, bufout)
         endif
         call wr_locatecp(wtd%meta, wtd%ic, wtd%ws, wtd%trk, wtd%discr, wtd%numcps, wtd%n_miss,         &
                        wtd%numtot, wtd%allcps, x_locate, sub_ierror)
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

      if (wtd%ic%return.le.1 .and. my_ierror.eq.0) then

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

         call total_forces_moments(wtd, x_locate)

         ! if requested: check for NaNs, write diagnostic information

         if (check_nan.ge.1) then
            if (x_locate.ge.5) call write_log(' ...wr_contact: calling checknan')
            call wtd_check_nan(wtd)
         endif

      endif

   end subroutine wr_contact_pos

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
            call writmt (gd%meta, gd%ic, gd%cgrid_cur, gd%potcon_cur, gd%geom%hs1, gd%mater, gd%fric,   &
                gd%kin, gd%outpt1, wtd%ic%is_left_side())
            end associate
         enddo
      endif
         
   end subroutine wtd_check_nan

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

   subroutine wr_setup_cp(meta, ws, trk, numcps, icp, cp, wtd_ic, mater, discr, fric, kin, solv, idebug)
!--purpose: define the contact problem for an initial contact point for a W/R contact case. 
      implicit none
!--subroutine arguments:
      type(t_metadata)  :: meta
      type(t_wheelset)  :: ws
      type(t_trackdata) :: trk
      integer           :: numcps, icp, idebug
      type(t_cpatch)    :: cp
      type(t_ic)        :: wtd_ic
      type(t_material)  :: mater
      type(t_discret)   :: discr
      type(t_friclaw)   :: fric
      type(t_kincns)    :: kin
      type(t_solvers)   :: solv
!--local variables:
      character(len=5)          :: nam_side
      logical                   :: new_gd
      integer                   :: ip, ii, iy, mx, my, npot, ierror
      real(kind=8)              :: sw_ref, fac_warn
      type(t_marker)            :: mref_pot, mref_rai, mref_whl, whl_trk
      type(t_probdata), pointer :: gd

      if (idebug.ge.3) then
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

      if (idebug.ge.3) then
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

      new_gd = .false.
      if (.not.associated(cp%gd)) then
         if (wtd_ic%x_cpatch.ge.1) then
            write(bufout,'(a,i3)') ' wr_setup_cp: allocate gd for icp=',icp
            call write_log(1, bufout)
         endif
         new_gd = .true.
         allocate(cp%gd)
         call gd_init(cp%gd)
      endif
      gd => cp%gd

      !  - general metadata:
      gd%meta%expnam    = meta%expnam
      gd%meta%dirnam    = meta%dirnam
      gd%meta%npatch    = numcps
      gd%meta%ipatch    = icp

      !  - iteration-related metadata are set at the actual solution in wr_solve_cp

      ! copy wheel and rail position data to meta-data
      ! note: mirrorring for left wheel is ignored, internal data written to gd%meta

      gd%meta%tim       = 0d0
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
      gd%ic%matfil_surf = 0

      ! copy solver settings

      gd%solv      = solv

      ! copy material parameters

      gd%mater     = mater

      ! copy kinematic constants

      gd%kin       = kin

      !  - set parameters describing the potential contact area as needed in undef.dist calculation

      gd%potcon_inp%ipotcn =  1
      gd%potcon_inp%dx     = cp%dx_eff
      gd%potcon_inp%mx     = nint((cp%xend-cp%xsta)/cp%dx_eff)
      gd%potcon_inp%xl     = cp%xsta - cp%mpot%x()

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
 
      if (wtd_ic%pvtime.ne.2) call merge_prev_potcon(icp, cp%gd, wtd_ic%x_cpatch)

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

      if (idebug.ge.2) then
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

      if (new_gd .or. wtd_ic%pvtime.eq.2) then
         gd%ic%pvtime = 2
      elseif (gd%meta%ncase.ge.meta%ncase) then
         gd%ic%pvtime = 3
      else
         gd%ic%pvtime = wtd_ic%pvtime
      endif

      ! set I-digit for this gd
      !  - no initial estimate for new contact patches
      !  - no initial estimate in case of large change in potential contact area (npot)
      !  - user setting I from wtd when a new case is started
      !  - use initial estimate if this gd was used in previous iteration
      !  - no initial estimate if this gd was inactive in previous iteration

      if (new_gd) then
         gd%ic%iestim = 0
         if (wtd_ic%x_cpatch.ge.1) call write_log('    no initial estimate (new patch)...')
      elseif (abs(gd%potcon_inp%dx - gd%potcon_cur%dx).gt.1d-4*gd%potcon_inp%dx .or.                    &
              abs(gd%potcon_inp%dy - gd%potcon_cur%dy).gt.1d-4*gd%potcon_inp%dy) then
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
         if (idebug.ge.3) call write_log('%   iy     si      ai');
         do iy = 1, my
            gd%mater%surf_inclin(iy,1) = gd%cgrid_inp%y(iy*mx)
            gd%mater%surf_inclin(iy,2) = cp%curv_incln%vy(iy)
            if (idebug.ge.3) then
               write(bufout,'(i6,f8.3,f11.6)') iy, gd%mater%surf_inclin(iy,1), gd%mater%surf_inclin(iy,2)
               call write_log(1, bufout)
            endif
         enddo
      endif

      ! copy/interpolate friction law data

      if (wtd_ic%varfrc.eq.0) then
         gd%ic%varfrc = 0
         call fric_copy(fric, gd%fric)
      elseif (wtd_ic%is_conformal()) then
         gd%ic%varfrc = 2
         call fric_interp(fric, gd%potcon_inp%my, cp%curv_incln%vy, gd%fric)
      else
         gd%ic%varfrc = 0
         call fric_interp(fric, gd%meta%deltcp_w, gd%fric)
      endif

      if (idebug.ge.3) then
         write(bufout,*) ' Using V = ', gd%ic%varfrc,' with NVF =', gd%fric%nvf,' slices'
         call write_log(1, bufout)
      endif

      ! constant friction per patch: enable scaling output forces
      ! Note: potential conflicts for F=1 with Fx<>0.

      gd%kin%use_muscal = gd%ic%varfrc.eq.0
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

         call wr_ud_planar (wtd_ic, ws, trk, icp, cp, gd, idebug)

      endif

      ! set planform data for reduced interaction between sub-patches

      if (cp%nsub.le.1) then
         gd%geom%iplan = 1
      else
         gd%geom%iplan = 4
         associate( np => gd%geom%npatch )
         np = cp%nsub

         call reallocate_arr(gd%geom%ysep, np)
         call reallocate_arr(gd%geom%facsep, np, np)
         gd%geom%ysep(1:np) = 0d0
         gd%geom%facsep(1:np, 1:np) = 0d0

         ! convert track y_sep position in rail profile to lateral s_c position in tangent plane

         gd%geom%ysep(1:np-1) = (cp%y_sep(1:np-1) - cp%mref%y()) / cos(cp%mref%roll())

         ! fill tri-diagonal matrix fac
        
         do ip = 1, np-1                       ! nsub patches --> nsub-1 interactions
            gd%geom%facsep(ip  ,ip  ) = 1d0          ! diagonal entry
            gd%geom%facsep(ip  ,ip+1) = cp%f_sep(ip) ! upper diagonal  (i,i+1)
            gd%geom%facsep(ip+1,ip  ) = cp%f_sep(ip) ! symmetry: lower diagonal (i+1,i)
         enddo
         gd%geom%facsep(np  ,np  ) = 1d0             ! diagonal entry

         if (idebug.ge.2) then
            write(bufout,'(a,i3)') ' IPLAN=4: npatch=',np
            call write_log(1, bufout)

            do ip = 1, np
               write(bufout,'(a,f9.3,a,7f7.3)') ' sc_sep=', gd%geom%ysep(ip),', fac=',             &
                     (gd%geom%facsep(ip,ii), ii=1, np)
               call write_log(1, bufout)
            enddo

            mx   = gd%potcon_inp%mx
            my   = gd%potcon_inp%my
            ! write(bufout,*) 'my=',my,', sc1=',gd%cgrid_inp%y(1),', scn=',gd%cgrid_inp%y(mx*my)
            ! call write_log(1, bufout)
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

         ! compute the creepages or the rigid slip function w.r.t. reference marker

         call wr_rigid_slip(wtd_ic, ws, trk, discr%dqrel, cp, gd, idebug)

         ! set spin center using offset reference marker -- pot.contact

         gd%kin%spinxo = mref_pot%x()
         gd%kin%spinyo = mref_pot%y()
      endif
 
      ! subsurface stress computation is performed separately, if needed

      gd%ic%stress = 0

      if (idebug.ge.4) call write_log('--- end subroutine wr_setup_cp ---')
      end associate

   end subroutine wr_setup_cp

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
      integer                   :: iatbnd
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
         if (idebug.ge.2) then
            write(bufout,800) ' total force(cs) =  [', gd%kin%fxrel,',',gd%kin%fyrel,',',gd%kin%fntrue,']'
            call write_log(1, bufout)
 800        format(4(a,g14.6))
         endif

         call aggregate_forces(cp, gd%ic, gd%cgrid_cur, gd%kin, gd%mater, gd%outpt1, idebug)

         if (idebug.ge.2) then
            write(bufout,800) ' total force(cp) =  [', gd%kin%fxrel,',',gd%kin%fyrel,',',gd%kin%fntrue,']'
            call write_log(1, bufout)
         endif
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
      real(kind=8)     :: dxdy

      associate( npot   => cgrid%ntot,    muscal => kin%muscal,    fntrue => kin%fntrue,        &
                 fxrel  => kin%fxrel,     fyrel  => kin%fyrel,                                  &
                 mxtru1 => outpt1%mxtrue, mytru1 => outpt1%mytrue, mztru1 => outpt1%mztrue )

      ! create 3-d version of curved reference surface

      call grid_extrude_profile(cp%curv_ref, cgrid%nx, cgrid%x(1)+cp%mref%x(), cgrid%dx, curv_cp)

      ! transform from global track coordinates to local planar sp coordinates

      call cartgrid_2loc(curv_cp, cp%mref)

      ! copy tractions on curved sc reference to work variable, rotate to planar sp coordinates

      if (idebug.ge.2) call grid_print(cp%curv_ref, 'curv-srfc', 5)
      if (idebug.ge.4) call grid_print(curv_cp, 'curv-cp', 5)

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

      if (ic%norm.eq.0) then
         fntrue = dxdy * gf3_sum(AllElm, ps_cp, ikZDIR)
      endif

      if (ic%force3.eq.0) then
         fxrel = dxdy * gf3_sum(AllElm, ps_cp, ikXDIR) / (fntrue*muscal + tiny)
      endif
      if (ic%force3.le.1) then
         fyrel = dxdy * gf3_sum(AllElm, ps_cp, ikYDIR) / (fntrue*muscal + tiny)
      endif

      ! Compute torsional moments about planar contact x, sp and np-axes
      ! using (x,sp,np)-coordinates provided in curv_cp

      mxtru1 = dxdy * ddot(npot, curv_cp%y, 1, ps_cp%vn, 1) - dxdy * ddot(npot, curv_cp%z, 1, ps_cp%vy, 1)
      mytru1 = dxdy * ddot(npot, curv_cp%z, 1, ps_cp%vx, 1) - dxdy * ddot(npot, curv_cp%x, 1, ps_cp%vn, 1)
      mztru1 = dxdy * ddot(npot, curv_cp%x, 1, ps_cp%vy, 1) - dxdy * ddot(npot, curv_cp%y, 1, ps_cp%vx, 1)

      if (idebug.ge.2) then
         write(bufout,'(3(a,g14.7))') ' Total moments Mx=',mxtru1,', Msp=',mytru1,', Mnp=',mztru1
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
      type(t_vec)       :: fcntc, tcntc
      type(t_marker)    :: mwhl_trk

      associate(my_rail => trk%rai, my_wheel => ws%whl, gd => cp%gd)

      ! rotate forces from contact-reference coordinates to global coordinates

      fcntc = vec( gd%kin%fxrel * gd%kin%fntrue * gd%kin%muscal,                               &
                   gd%kin%fyrel * gd%kin%fntrue * gd%kin%muscal, gd%kin%fntrue )

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

      cp%ftrk = cp%mref%rot * fcntc
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

   subroutine total_forces_moments(wtd, idebug)
!--purpose: compute total forces and moments, location of minimum moment
   implicit none
!--subroutine arguments:
      integer,          intent(in) :: idebug
      type(t_ws_track)             :: wtd
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

end module m_wr_solvecp
