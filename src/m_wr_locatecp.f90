!------------------------------------------------------------------------------------------------------------
! m_wr_locatecp - identify contact patches for one case of w/r contact (module 1)
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_wr_locatecp

use m_hierarch_data
use m_wrprof_data
implicit none
private

public  wr_locatecp
private locate_regions
private locate_patches
private complete_patches
private compute_wr_surfc
private compute_round_whl
private compute_oor_whl
private compute_wr_locus
private locus_prismatic
private locus_iterate_varprof
private locus_iterate_roller
private make_surfc_inclination
private locate_interpen_1d
private locate_interpen_2d
private compute_contact_length_1d
private get_uv_extent
private sort_cpatches_yvalue
private combine_cpatches
private set_cpatch_reference
private compute_contact_angle
private set_planar_potcon
private compute_curved_potcon
private estimate_curvatures_v1
private estimate_curvatures_v2
private estimate_minimum
private find_gap_locmin
public  find_bounding_box_1d
public  find_bounding_box_2d
private prw_determine_foldback

contains

!------------------------------------------------------------------------------------------------------------

   subroutine wr_locatecp (ic, ws, trk, discr, numcps, allcps, idebug, my_ierror)
!--purpose: locate the initial contact point(s) for a W/R contact case.
      implicit none
!--subroutine arguments:
      type(t_ic)                    :: ic
      type(t_wheelset)              :: ws
      type(t_trackdata)             :: trk
      type(t_discret)               :: discr
      integer                       :: numcps, idebug
      type(p_cpatch)                :: allcps(MAX_NUM_CPS)
      integer,          intent(out) :: my_ierror
!--local variables:
      integer                  :: is_right, sub_ierror, irgn, num_rgn
      real(kind=8)             :: sgn
      type(t_region)           :: all_rgn(MAX_NUM_RGN)
      type(t_rail),    pointer :: my_rail
      type(t_wheel),   pointer :: my_wheel

      my_ierror = 0

      call timer_start(itimer_locatecp)
      if (idebug.ge.3) then
         write(bufout,'(/a)') ' --- Start subroutine wr_locatecp ---'
         call write_log(2, bufout)
      endif

      if (idebug.ge.3) then
         write(bufout,'(4(a,f10.6))') ' Wheelset y=',ws%y,', z=',ws%z,', yaw=',ws%yaw,', roll=',ws%roll
         call write_log(1, bufout)
      endif

      my_rail   => trk%rai
      my_wheel  => ws%whl
      if (ic%is_left_side()) then
         is_right  = 0
      else
         is_right  = 1
      endif

      ! set the sign to -1 for left and +1 for right rail/wheel combination

      sgn = 2 * is_right - 1

      if (.not.grid_has_spline(my_rail%prr%grd_data)) then
         call write_log(' Internal error: no spline in rail profile, aborting.')
         call abort_run()
      endif

      !--------------------------------------------------------------------------------------------------
      ! phase 1: perform a quick scan to find regions with possible contact
      !--------------------------------------------------------------------------------------------------

      call locate_regions(ic, my_rail, my_wheel, num_rgn, all_rgn, my_ierror)

      !--------------------------------------------------------------------------------------------------
      ! phase 2: for each region, determine contact patches
      !--------------------------------------------------------------------------------------------------

      numcps     = 0
      ws%has_overlap = .false.
      ws%gap_min = 1d10
      ws%z_cnt0  = 1d10

      do irgn = 1, num_rgn
         if (num_rgn.gt.1) then
            write(bufout,'(/,2(a,i3))') ' --- calling locate_patches for region',irgn,', #regions =',num_rgn
            call write_log(2, bufout)
         endif
         call locate_patches(ic, ws, trk, discr, all_rgn(irgn), numcps, allcps, idebug, my_ierror)
      enddo

      !--------------------------------------------------------------------------------------------------
      ! phase 3: sort and combine patches, fill in potential contact area
      !--------------------------------------------------------------------------------------------------

      if (ws%has_overlap .and. my_ierror.eq.0) then

         call complete_patches( ic, ws, trk, discr, numcps, allcps, idebug, sub_ierror )
         my_ierror = sub_ierror

      endif
      call timer_stop(itimer_locatecp)

   end subroutine wr_locatecp

!------------------------------------------------------------------------------------------------------------

   subroutine locate_regions (ic, my_rail, my_wheel, num_rgn, all_rgn, my_ierror)
!--purpose: locate the initial contact point(s) for a W/R contact case.
      implicit none
!--subroutine arguments:
      type(t_ic)                      :: ic
      type(t_rail),     target        :: my_rail
      type(t_wheel),    target        :: my_wheel
      integer,          intent(out)   :: num_rgn
      type(t_region),   target        :: all_rgn(MAX_NUM_RGN)
      integer,          intent(out)   :: my_ierror
!--local variables:
      integer                 :: nr, nw
      real(kind=8)            :: sr0, sw0, sr_len, sw_len, cos_a, sin_a
      type(t_region), pointer :: region
      type(t_grid),   pointer :: prf

      my_ierror   = 0
      num_rgn     = 0

      if (ic%use_oblique()) then

         prf    => my_rail%prr%grd_data
         nr     =  prf%ntot
         sr0    =  prf%s_prf(1)
         sr_len =  prf%s_prf(nr) - sr0
         write(bufout,'(2(a,f9.3),a)') ' Rail  profile has s=[', 0d0,',',sr_len,']'
         call write_log(1, bufout)

         prf    => my_wheel%prw%grd_data
         nw     =  prf%ntot
         sw0    =  prf%s_prf(1)
         sw_len =  prf%s_prf(nw) - sw0
         write(bufout,'(2(a,f9.3),a)') ' Wheel profile has s=[', 0d0,',',sw_len,']'
         call write_log(1, bufout)

         if (abs(sr_len-132.4d0).lt.0.1d0) then
            ! UIC60: total arc-length 132.389 mm

            call write_log(' --- locate_regions: setting up search region for mbench_right13 ---')

            ! define view point/direction for mbench_right13

            num_rgn = 1
            region => all_rgn(1)

            region%mview = marker(vec(10d0, 700d0, -5d0), rotmat_pure_roll(-0.7d0))

            ! select regions on rail and wheel profiles

            region%xsta   =       -30d0
            region%xend   =        10d0
            region%sr_sta = sr0 +  19d0
            region%sr_end = sr0 +  64d0
            region%sw_sta = sw0 +  70d0
            region%sw_end = sw0 + 119d0

         elseif (abs(sr_len-60.1d0).lt.0.1d0) then

            ! define view point/direction for guard_20_1, total arc-length 60.075 mm

            call write_log(' --- locate_regions: setting up search region for guard_20_1 ---')

            num_rgn = 1
            region => all_rgn(1)

            cos_a = 0.05d0
            sin_a = sqrt(1d0 - cos_a**2)
            region%mview = marker( vec(0d0, 50d0*sin_a, -50d0*cos_a), rotmat_pure_roll(acos(cos_a)) )

            ! select regions on rail and wheel profiles

            region%xsta   =       -30d0         ! no yaw: [-25,25], with yaw: [75,125]
            region%xend   =       130d0
            region%sr_sta = sr0 +  -5d0
            region%sr_end = sr0 +  65d0
            region%sw_sta = sw0 +  -5d0
            region%sw_end = sw0 + 120d0
            ! call marker_print( region%mview, 'm_view', 5)

         elseif (abs(sr_len-264.7d0).lt.0.1d0) then

            ! define view point/direction for willets_region, total arc-length 264.728 mm

            call write_log(' --- locate_regions: setting up search region for willets_region ---')

            num_rgn = 2

            region  => all_rgn(1)       ! 1: guard rail contact
            region%mview = marker(vec(0d0, 690d0, -25d0), rotmat_pure_roll(70d0*pi/180d0))
            region%xsta   =       -20d0
            region%xend   =       100d0
            region%sr_sta = sr0 +  74d0
            region%sr_end = sr0 +  98d0
            region%sw_sta = sw0 + 158d0
            region%sw_end = sw0 + 201d0

            region  => all_rgn(2)       ! 2: stock rail contact
            region%mview = marker(vec(0d0, 740d0, -25d0), rotmat_pure_roll( 0d0*pi/180d0))
            region%xsta   =       -60d0
            region%xend   =        60d0
            region%sr_sta = sr0 + 201d0
            region%sr_end = sr0 + 265d0
            region%sw_sta = sw0 +  25d0
            region%sw_end = sw0 + 106d0

         else

            call write_log(' ERROR: locate_regions: rail profile not recognized, using default view !!!')

         endif

      endif

      if (num_rgn.eq.0) then

         ! default: view == track direction, region == everything

         num_rgn = 1
         region  => all_rgn(1)

         call marker_init(region%mview)
         region%sr_sta = -1d3
         region%sr_end =  1d3
         region%sw_sta = -1d3
         region%sw_end =  1d3

         region%xsta   = -100d0
         region%xend   =  100d0

      endif

   end subroutine locate_regions

!------------------------------------------------------------------------------------------------------------

   subroutine locate_patches(ic, ws, trk, discr, region, numcps, allcps, idebug, my_ierror)
!--purpose: locate the initial contact point(s) for one region for a W/R contact case.
      implicit none
!--subroutine arguments:
      type(t_ic)                      :: ic
      type(t_wheelset)                :: ws
      type(t_trackdata)               :: trk
      type(t_discret)                 :: discr
      type(t_region)                  :: region
      integer,          intent(inout) :: numcps
      type(p_cpatch)                  :: allcps(MAX_NUM_CPS)
      integer,          intent(in)    :: idebug
      integer,          intent(out)   :: my_ierror
!--local variables:
      integer                  :: sub_ierror, numcp0, icp
      logical                  :: rgn_has_overlap
      real(kind=8)             :: sgn, rgn_gap_min
      type(t_vec)              :: vec_trk
      type(t_marker)           :: rai_vw
      type(t_rail),    pointer :: my_rail
      type(t_wheel),   pointer :: my_wheel
      type(t_grid)             :: prr_rgn, sf_rai, sf_whl
      type(t_gridfnc3)         :: sf_incl, uv_whl

      my_ierror = 0

      if (ic%use_oblique()) call write_log(' --- locate_patches: using oblique view direction ---')

      my_rail   => trk%rai
      my_wheel  => ws%whl

      ! set the sign to -1 for left and +1 for right rail/wheel combination

      if (ic%is_left_side()) then
         sgn = -1
      else
         sgn =  1
      endif

      !--------------------------------------------------------------------------------------------------
      ! phase 2: for each region, determine contact patches
      !--------------------------------------------------------------------------------------------------

      ! select search region from (1d) rail profile 'prr_grd'

      if (ic%use_oblique()) then
         call grid_trim(my_rail%prr%grd_data, prr_rgn, 1, 1, s_low=region%sr_sta, s_hig=region%sr_end,  &
                with_spline=.true.)
      else
         call grid_copy(my_rail%prr%grd_data, prr_rgn, with_spline=.true.)
      endif

      ! determine rail marker in terms of the view coordinate system

      rai_vw  = marker_2loc( my_rail%m_trk, region%mview )

      ! convert rail profile from rail coords to view coords

      call cartgrid_2glob(prr_rgn, rai_vw)

      if (idebug.ge.3) then
         write(bufout,123) 'rail(1): y"=',sgn*prr_rgn%y(1),', z"=',prr_rgn%z(1)
         call write_log(1, bufout)
 123     format(1x, 3(a, f12.6, :))
      endif

      ! no trimming, add top-view for consistent uni-valued interpolation
      ! TODO: avoid top-view by appropriately defined 'region' and 'view direction'

      if (idebug.ge.3) call spline_set_debug(2)
      call spline_add_topview(prr_rgn%spl, .true., sub_ierror)
      call spline_set_debug(0)
      if (idebug.ge.3) call spline_print(prr_rgn%spl, 'prr_rgn%spl', 3)

      ! define the 'gap_mesh' and compute sf_whl, sf_rai accordingly

      if (my_ierror.eq.0) then
         if (.not.use_brute(ic, ws%whl%prw)) then
            if (idebug.ge.2) call write_log(' --- calling compute_wr_locus ---')
            call compute_wr_locus (ic, ws, trk, prr_rgn, region, discr%ds, sf_whl, sf_rai, idebug,      &
                        rgn_has_overlap, sub_ierror)
            if (my_ierror.eq.0) my_ierror = sub_ierror
         else
            if (idebug.ge.2) call write_log(' --- calling compute_wr_surfc ---')
            call compute_wr_surfc (ic, ws, trk, prr_rgn, region, discr%dx, discr%ds, sf_whl, sf_rai,    &
                        uv_whl, idebug, rgn_has_overlap)
         endif
      endif

      if (rgn_has_overlap .and. my_ierror.eq.0) then

         ! determine rail surface inclinations wrt view direction, at the gap mesh of sf_rai

         call make_surfc_inclination(prr_rgn, sf_rai, sf_incl, idebug)

         ! identify initial contact points and interpenetration regions, add to allcps
         !    local minima of gap in view direction / maxima of interpenetration

         numcp0 = numcps

         if (.not.use_brute(ic, ws%whl%prw)) then
            if (idebug.ge.2) call write_log(' --- calling locate_interpen_1d ---')
            call locate_interpen_1d(ic, sf_whl, sf_rai, sf_incl, region%mview, sgn, ws%nom_radius,      &
                        trk%nom_radius, rgn_gap_min, ws%delt_min, ws%zw_min, ws%a1, ws%b1, numcps,      &
                        allcps, idebug)
         else
            if (idebug.ge.2) call write_log(' --- calling locate_interpen_2d ---')
            call timer_start(itimer_interp6)
            call locate_interpen_2d(ic, sf_whl, sf_rai, uv_whl, sf_incl, sgn, discr%dx, ws%nom_radius,  &
                        rgn_gap_min, ws%delt_min, ws%zw_min, ws%a1, ws%b1, numcps, allcps, idebug)
            call timer_stop(itimer_interp6)
         endif

         ! convert results from view direction to track orientation

         if (ic%use_oblique()) then

            do icp = numcp0+1, numcps
               associate( cp => allcps(icp)%cp, m_view => region%mview )

               ! start/end positions of interpenetration area, using micp

               vec_trk      = vec_2glob( vec( cp%xsta, cp%micp%y(), cp%micp%z() ), m_view )
               cp%xsta      = vec_trk%x()
               vec_trk      = vec_2glob( vec( cp%xend, cp%micp%y(), cp%micp%z() ), m_view )
               cp%xend      = vec_trk%x()
               vec_trk      = vec_2glob( vec( cp%micp%x(), cp%ysta, cp%zsta ), m_view )
               cp%ysta      = vec_trk%y()
               cp%zsta      = vec_trk%z()
               vec_trk      = vec_2glob( vec( cp%micp%x(), cp%yend, cp%zend ), m_view )
               cp%yend      = vec_trk%y()
               cp%zend      = vec_trk%z()

               ! position with maximum interpenetration

               cp%micp      = vec_2glob( cp%micp, m_view )

               ! weighted center position and orientation

               vec_trk      = vec_2glob( vec( cp%wgt_xgap, cp%wgt_ygap, cp%wgt_zgap ), m_view )
               cp%wgt_xgap  = vec_trk%x()
               cp%wgt_ygap  = vec_trk%y()
               cp%wgt_zgap  = vec_trk%z()
               cp%wgt_agap  = cp%wgt_agap + m_view%roll()

               end associate
            enddo

            ! convert gap-value and orientation at position with max interpenetration

            rgn_gap_min  = rgn_gap_min * cos(ws%delt_min) / cos(ws%delt_min + region%mview%roll() )
            ws%delt_min  = ws%delt_min + region%mview%roll()

         endif

         ! merge information for this region to overall information

         ws%has_overlap = .true.
         ws%gap_min = min(ws%gap_min, rgn_gap_min)
         ws%z_cnt0  = min(ws%z_cnt0, ws%z + ws%gap_min)

      endif

      ! cleanup local variables

      call grid_destroy(prr_rgn)
      call grid_destroy(sf_rai)
      call grid_destroy(sf_whl)
      call gf3_destroy(sf_incl)
      call gf3_destroy(uv_whl)

      if (my_ierror.ne.0) then
         my_ierror = CNTC_err_search
         ws%gap_min = 1d10
         ws%z_cnt0  = 1d10
      endif

   end subroutine locate_patches

!------------------------------------------------------------------------------------------------------------

   subroutine complete_patches( ic, ws, trk, discr, numcps, allcps, idebug, my_ierror )
!--purpose: sort and combine contact patches and complete the potential contact area definition
      implicit none
!--subroutine arguments:
      type(t_ic)                      :: ic
      type(t_wheelset)                :: ws
      type(t_trackdata)               :: trk
      type(t_discret)                 :: discr
      integer,          intent(inout) :: numcps
      type(p_cpatch)                  :: allcps(MAX_NUM_CPS)
      integer,          intent(in)    :: idebug
      integer,          intent(out)   :: my_ierror
!--local variables:
      integer                  :: icp, sub_ierror
      real(kind=8)             :: sgn, cref_x
      type(t_grid)             :: prr_trk
      type(t_rail),    pointer :: my_rail
      type(t_cpatch),  pointer :: cp

      my_ierror  = 0

      my_rail        => trk%rai

      ! set the sign to -1 for left and +1 for right rail/wheel combination

      if (ic%is_left_side()) then
         sgn = -1
      else
         sgn =  1
      endif

      ! sort the contact problems with |y|-position descending, i.e. from field side to track center

      if (numcps.ge.2) then
         ! call write_log(' --- calling sort_cpatches_yvalue ---')
         call sort_cpatches_yvalue( numcps, allcps, idebug )
      endif

      ! combine contact problems that are overlapping or lie too close together

      if (numcps.ge.2) then
         ! call write_log(' --- calling combine_cpatches ---')
         call combine_cpatches( ic, numcps, allcps, discr%angl_sep, discr%dist_sep, discr%dist_comb, idebug )
      endif

      ! turn contact reference angles for patches that lie close together

      if (numcps.ge.2) then
         ! call write_log(' --- calling turn_cpatch_refangle ---')
         call turn_cpatch_refangle( numcps, allcps, discr%angl_sep, discr%dist_sep, discr%dist_turn,    &
                idebug )
      endif

      ! loop over contact patches for completion of the contact patch data

      do icp = 1, numcps

         cp => allcps(icp)%cp

         if (icp.eq.1 .and. .not.my_rail%prr%is_varprof()) then

            ! initialize the profile in track coordinates; prr_grid is given in rail profile coords

            call grid_copy(my_rail%prr%grd_data, prr_trk, with_spline=.true.)
            call cartgrid_2glob(prr_trk, my_rail%m_trk)
            call spline_add_topview(prr_trk%spl, .true., sub_ierror)

         elseif (my_rail%prr%is_varprof()) then

            if (idebug.ge.2) then
               write(bufout,'(a,i3)') ' setting up "current slice" for patch',icp
               call write_log(1, bufout)
            endif

            if (ic%use_initial_cp()) then
               cref_x = cp%micp%x()
            else
               cref_x = cp%wgt_xgap
            endif

            ! call write_log(' ...varprof_intpol_xunif')
            call varprof_intpol_xunif(my_rail%prr, ws%s, 1, cref_x, 0d0, prr_trk, sub_ierror)
            my_ierror = sub_ierror

            call cartgrid_2glob(prr_trk, my_rail%m_trk)

         endif

         ! fill in the contact angle and reference marker

         ! call write_log(' --- calling set_cpatch_reference ---')
         call set_cpatch_reference( ic, icp, cp, prr_trk, my_rail%m_trk, sgn, trk%nom_radius, idebug)

         ! compute potential contact area for planar contact approach

         ! call write_log(' --- calling set_planar_potcon ---')
         call set_planar_potcon( ic, icp, cp, prr_trk, my_rail%m_trk, sgn, trk%nom_radius, discr%dx,    &
                        discr%ds, discr%npot_max, idebug )

         ! estimate curvatures a1, b1 for the contact patches

         ! call write_log(' --- calling estimate_curvatures_v2 ---')
         call estimate_curvatures_v2( ic, cp, ws, prr_trk, trk%nom_radius, idebug )

         ! compute curved contact reference surface, pot.contact for conformal contact approach

         if (ic%is_conformal()) then
            ! call write_log(' --- calling compute_curved_potcon ---')
            call compute_curved_potcon( ic, ws, prr_trk, my_rail%m_trk, trk%nom_radius, cp, idebug )
         endif

      enddo ! icp

      ! cleanup local variables

      call grid_destroy(prr_trk)

   end subroutine complete_patches

!------------------------------------------------------------------------------------------------------------

   subroutine compute_wr_surfc (ic, ws, trk, prr_rgn, region, dx_cp, ds_cp, sf_whl, sf_rai, uv_whl,     &
                                idebug, has_overlap)
!--purpose: define the 'gap_mesh' and compute the wheel and rail surfaces sf_whl and sf_rai accordingly.
      implicit none
!--subroutine arguments:
      type(t_ic)                    :: ic
      type(t_wheelset)              :: ws
      type(t_trackdata)             :: trk
      real(kind=8),     intent(in)  :: dx_cp, ds_cp  ! discretisation steps DX, DS in long/lat directions
      integer,          intent(in)  :: idebug
      type(t_grid)                  :: prr_rgn, sf_whl, sf_rai
      type(t_gridfnc3)              :: uv_whl
      type(t_region)                :: region
      logical,          intent(out) :: has_overlap
!--local variables:
      logical                 :: is_prismatic
      integer                 :: nslc, nr, nx, ny, is_right, sub_ierror
      real(kind=8)            :: fac_dyds, xmin_vw, xmax_vw, dx_vw, ymin, ymax, dy, sgn, z_axle
      type(t_marker)          :: rw_trk, rw_vw, rr_vw
      type(t_rail),   pointer :: my_rail
      type(t_wheel),  pointer :: my_wheel
      type(t_grid)            :: gap_mesh, prr_rr, prr_unif, rai_curv, whl_curv, bbr, bbw
      type(t_gridfnc3)        :: whl_curv_uv

      if (idebug.ge.3) then
         write(bufout,'(/a)') ' --- Start subroutine compute_wr_surfc ---'
         call write_log(2, bufout)
      endif

      my_rail   => trk%rai
      my_wheel  => ws%whl
      if (ic%is_left_side()) then
         is_right  = 0
      else
         is_right  = 1
      endif

      is_prismatic = (.not.ic%is_roller() .and. .not.my_rail%prr%is_varprof())

      ! set the sign to -1 for left and +1 for right rail/wheel combination

      sgn = 2 * is_right - 1

      ! form wheel surface whl_curv

      if (.not.my_wheel%prw%is_varprof()) then

         call compute_round_whl(ic, ws, region, dx_cp, nslc, whl_curv, idebug)

      else ! my_wheel%prw%is_varprof()

         call timer_start(itimer_interp1)
         call compute_oor_whl(ic, ws, dx_cp, nslc, whl_curv, whl_curv_uv, idebug)
         call timer_stop(itimer_interp1)

      endif ! prw%is_varprof()

      if (idebug.ge.2) call grid_check_nan(whl_curv, 'whl_curv(rw)', idebug)

      ! convert the wheel mesh in one step from profile to view coordinates

      call timer_start(itimer_interp5)
      rw_trk = marker_2glob( my_wheel%m_ws, ws%m_trk )
      rw_vw  = marker_2loc( rw_trk, region%mview )
      call cartgrid_2glob(whl_curv, rw_vw)

      if (.false.) then
         call grid_print(whl_curv, 'whl_curv', 5)
         call abort_run()
      endif

      ! get range of x-coordinates of wheel mesh

      call grid_get_xrange(whl_curv, xmin_vw, xmax_vw)
      dx_vw = (xmax_vw - xmin_vw) / real(nslc-1)

      if (idebug.ge.3) then
         write(bufout,'(3(a,f10.3))') ' whl_curv: x in [', xmin_vw,',', xmax_vw,'], dx=',dx_vw
         call write_log(1, bufout)
      endif

      if (idebug.ge.1) call grid_check_nan(whl_curv, 'whl_curv(tr)', idebug)

      ! the (1d) (trimmed) rail profile 'prr_rgn' is given in view coordinates

      nr   = prr_rgn%ntot

      ! define a uniform grid 'gap_mesh' in view coordinates, for computing the gap function

      ! - y-range: considering all y of the rail profile
      !            fac_dyds=0.5 gives steps of size ds along the surface where surface inclined by 60deg
      !            the resulting grid may be a little smaller at ymax
      ! TODO: increase fac_dyds based on max. surfc inclination in view direction

      fac_dyds = 0.5d0
      if (ic%use_steep_slopes()) fac_dyds = 0.2d0

      ymin = prr_rgn%y(2)
      ymax = prr_rgn%y(nr-1)
      dy   = fac_dyds * ds_cp * 1.5d0   ! increase step by 1.5 to gain some performance

      call grid_create_uniform(gap_mesh, x0arg=xmin_vw, dxarg=dx_vw, x1arg=xmax_vw,                     &
                                         y0arg=ymin, dyarg=dy, y1arg=ymax, zarg=0d0)
      call grid_get_dimens(gap_mesh, nx, ny)
      call grid_copy(gap_mesh, sf_rai)
      call timer_stop(itimer_interp5)

      if (idebug.ge.3) call grid_print(gap_mesh, 'gap_mesh', 2)

      if (my_rail%prr%is_varprof()) then

         ! variable profile - interpolate to gap_mesh and store in sf_rai

         ! call write_log(' ...varprof_intpol_grid')
         rr_vw  = marker_2loc( my_rail%m_trk, region%mview )
         call varprof_intpol_grid_old(trk%rai%prr, rr_vw, ws%s, sf_rai)

      elseif (ic%is_roller() .and. ic%use_oblique()) then

         ! wheelset on roller rig, with roller axle aligned with view coordinate y-axis

         ! interpolate rail profile to 1d intermediate 'prr_unif' 

         call grid_create_uniform(prr_unif, nxarg=1, x0arg=0d0, x1arg=0d0,                              &
                                            y0arg=ymin, dyarg=dy, y1arg=ymax, zarg=0d0)

         call spline_get_xz_at_y(prr_rgn%spl, prr_unif, sub_ierror)

         ! form the 2D rail surface in view coordinates in surface 'sf_rai'

         z_axle = trk%nom_radius
         call grid_revolve_profile(prr_unif, nx, xmin_vw, dx_vw, z_axle, 999d0, sf_rai)

      elseif (ic%is_roller()) then

         ! wheelset on roller rig, roller axle different from view coordinate y-axis

         ! profile 'prr_rgn' is given in view coordinates; get profile in rail coords 

         rr_vw  = marker_2loc( my_rail%m_trk, region%mview )
         call grid_copy(prr_rgn, prr_rr)
         call cartgrid_2loc(prr_rr, rr_vw)

         ! revolve 'prr_rr' to 2D surface 'rai_curv'

         z_axle = trk%nom_radius
         call grid_revolve_profile(prr_rr, nx, xmin_vw, dx_vw, z_axle, 999d0, rai_curv)

         ! convert rail surface from rail coords to view coords

         call cartgrid_2glob(rai_curv, rr_vw)

         ! interpolate surface rai_curv to the gap_mesh

         call grid_copy(gap_mesh, sf_rai)
         call interp_cartz2unif(rai_curv, sf_rai, sub_ierror, 999d0)

      else ! is_prismatic

         ! form the 2D rail surface in view coordinates, with prr_rgn provided in view coordinates

         ! interpolate rail profile to 1d intermediate 'prr_unif' profile

         call timer_start(itimer_interp4)
         call grid_create_uniform(prr_unif, nxarg=1, x0arg=0d0, x1arg=0d0,                              &
                                            y0arg=ymin, dyarg=dy, y1arg=ymax, zarg=0d0)

         call spline_get_xz_at_y(prr_rgn%spl, prr_unif, sub_ierror)

         ! wheelset on track: prismatic rail, extrude prr_unif along x-direction

         call grid_extrude_profile(prr_unif, nx, xmin_vw, dx_vw, sf_rai)
         call timer_stop(itimer_interp4)

      endif

      ! Check that the wheel and rail surfaces overlap

      call timer_start(itimer_interp5)
      has_overlap = .true.
      call grid_get_boundbox(whl_curv, bbw)
      call grid_get_boundbox(sf_rai, bbr)
      call timer_stop(itimer_interp5)

      if (bbw%y(1).gt.bbr%y(8) .or. bbr%y(1).gt.bbw%y(8)) then
         has_overlap = .false.
         write(bufout,'(2a,/,2(2(a,f10.3),a,:,/))') ' WARNING: the wheel and rail are LATERALLY in ',   &
                        'different places.',                                                            &
                        '          the wheel grid has y = [', bbw%y(1),',', bbw%y(8),'],',              &
                        '          the  rail grid has y = [', bbr%y(1),',', bbr%y(8),'].'
         call write_log(3, bufout)
      endif

      if (ic%norm.le.0 .and. minval(bbw%z).gt.maxval(bbr%z)) then
         has_overlap = .false.
         write(bufout,'(2a,/,2(2(a,f10.3),a,:,/))') ' WARNING: the wheel and rail are VERTICALLY in ',  &
                        'different places.',                                                            &
                        '          the wheel grid has z = [', minval(bbw%z),',', maxval(bbw%z),'],',    &
                        '          the  rail grid has z = [', minval(bbr%z),',', maxval(bbr%z),'].'
         call write_log(3, bufout)
      endif

      if (bbw%x(1).gt.bbr%x(8) .or. bbr%x(1).gt.bbw%x(8)) then
         has_overlap = .false.
         call write_log(' Internal error: x-ranges dont overlap')
      endif

      ! interpolate wheel surface to the mesh used for the gap function

      call timer_start(itimer_interp2)
      call grid_copy(gap_mesh, sf_whl)
      call interp_cartz2unif(whl_curv, sf_whl, sub_ierror, -999d0)
      call timer_stop(itimer_interp2)

      ! interpolate (u,v) positions on wheel surface to the gap mesh

      if (my_wheel%prw%is_varprof() .and. whl_curv_uv%is_defined()) then
         call timer_start(itimer_interp3)
         ! call interp_set_debug(3)
         call gf3_new(uv_whl, 'uv_whl', sf_whl, nulify=.true.)
         call interp_curvgf2unifgf(whl_curv_uv, uv_whl, sub_ierror, -888d0)
         ! call interp_set_debug(0)
         call timer_stop(itimer_interp3)
      endif

      ! for debugging: bounding-boxes of whl_curv, sf_whl and sf_rai

      if (idebug.ge.3 .or. (idebug.ge.1 .and. .not.has_overlap)) then
         write(bufout,248) 'Grid whl_curv', bbw%x(1), bbw%x(8), bbw%y(1), bbw%y(8), bbw%z(1), bbw%z(8)
         call write_log(2, bufout)

         call grid_get_boundbox(sf_whl, bbw)
         write(bufout,248) 'Grid sf_whl', bbw%x(1), bbw%x(8), bbw%y(1), bbw%y(8), bbw%z(1), bbw%z(8)
         call write_log(2, bufout)

         write(bufout,248) 'Grid sf_rai', bbr%x(1), bbr%x(8), bbr%y(1), bbr%y(8), bbr%z(1), bbr%z(8)
         call write_log(2, bufout)

  248    format(1x,a,' has bounding box',/,                                                             &
                '      x=[',f10.3,',',f10.3,'], y=[',f10.3,',',f10.3,'], z=[',f10.3,',',f10.3,']')
      endif

      ! call grid_print(whl_curv, 'whl_curv', 1)
      if (idebug.ge.5) call grid_print(sf_whl, 'sf_whl', 1)

      ! cleanup local variables

      call grid_destroy(gap_mesh)
      call grid_destroy(prr_rr)
      call grid_destroy(prr_unif)
      call grid_destroy(rai_curv)
      call grid_destroy(whl_curv)
      call grid_destroy(bbr)
      call grid_destroy(bbw)
      call gf3_destroy(whl_curv_uv)

      if (idebug.ge.4) call write_log('--- end subroutine compute_wr_surfc ---')

   end subroutine compute_wr_surfc

!------------------------------------------------------------------------------------------------------------

   subroutine compute_round_whl (ic, ws, region, dx_cp, nslc, whl_curv, idebug)
!--purpose: compute 3d mesh for out-of-round wheel surface
      implicit none
!--subroutine arguments:
      type(t_ic)                    :: ic
      type(t_wheelset)              :: ws
      real(kind=8),     intent(in)  :: dx_cp  ! discretisation step DX in long. direction
      integer,          intent(out) :: nslc
      integer,          intent(in)  :: idebug
      type(t_grid)                  :: whl_curv
      type(t_region)                :: region
!--local variables:
      integer                 :: nw
      real(kind=8)            :: xmin_ws, xmax_ws, dx_ws, z_axle
      type(t_wheel),  pointer :: my_wheel
      type(t_grid)            :: prw_rgn

      if (idebug.ge.3) then
         write(bufout,'(/a)') ' --- Start subroutine compute_round_whl ---'
         call write_log(2, bufout)
      endif

      my_wheel  => ws%whl

      ! form wheel surface whl_curv

      if (.not.my_wheel%prw%is_varprof()) then

         ! estimate appropriate range of x-coordinates, using nx = nslc slices

         if (ws%nom_radius.lt.30.5d0 .or. ic%use_oblique() .or. ic%use_steep_slopes()) then

            ! - ball-on-plane: considering an arc of [-14,14deg], sin(14d) = 0.242
            !   also use this setting for steep guard rail contacts

            xmin_ws = -ws%nom_radius * 0.242d0
            xmax_ws =  ws%nom_radius * 0.242d0
         else

            ! - wheel/rail: considering an arc of [-7,7deg], sin(7d) = 0.122 = 1/8

            xmin_ws = -ws%nom_radius / 8d0
            xmax_ws =  ws%nom_radius / 8d0
         endif

         if (ic%use_oblique()) then

            ! clip [xmin,xmax] on limits defined for region
            xmin_ws = max(xmin_ws, region%xsta)
            xmax_ws = min(xmax_ws, region%xend)

         endif

         ! compute nslc such that dx_ws = O(dx_cntc)

         nslc  = nint( (xmax_ws - xmin_ws) / (1.5d0 * dx_cp) + 2d0 )
         if (mod(nslc,2).eq.0) nslc = nslc + 1
         dx_ws = (xmax_ws - xmin_ws) / real(nslc - 1)

         if (idebug.ge.3) then
            write(bufout,'(2(a,f8.3),a,i5,a,f8.4)') ' gap function: x=[',xmin_ws,',',xmax_ws,           &
                   '], using nslc =',nslc,' slices with dx=',dx_ws
            call write_log(1, bufout)
         endif

         ! select search region from (1d) wheel profile 'grd_data'

         if (ic%use_oblique()) then
            call grid_trim(my_wheel%prw%grd_data, prw_rgn, 1, 1, s_low=region%sw_sta,                   &
                   s_hig=region%sw_end, with_spline=.true.)
         else
            call grid_copy(my_wheel%prw%grd_data, prw_rgn, with_spline=.true.)
         endif

         nw = prw_rgn%ntot

         ! form the 2D wheel surface at profile y positions (1:nw) (wheel profile coordinates)

         z_axle  = -ws%nom_radius

         call grid_revolve_profile(prw_rgn, nslc, xmin_ws, dx_ws, z_axle, -999d0, whl_curv)

      else ! my_wheel%prw%is_varprof()

         call write_log('INTERNAL ERROR(round_whl): varprof not supported.')
         call abort_run()

      endif ! prw%is_varprof()

      call grid_destroy(prw_rgn)

      if (idebug.ge.4) call write_log('--- end subroutine compute_round_whl ---')

   end subroutine compute_round_whl

!------------------------------------------------------------------------------------------------------------

   subroutine compute_oor_whl (ic, ws, dx_cp, nslc, whl_curv, whl_uv, idebug)
!--purpose: compute 3d mesh for out-of-round wheel surface
      implicit none
!--subroutine arguments:
      type(t_ic)                    :: ic
      type(t_wheelset)              :: ws
      real(kind=8),     intent(in)  :: dx_cp    ! discretisation step DX in long. direction
      integer,          intent(out) :: nslc
      integer,          intent(in)  :: idebug
      type(t_grid)                  :: whl_curv
      type(t_gridfnc3)              :: whl_uv
!--local variables:
      logical,        parameter  :: compute_uv = .false.
      integer                    :: iu, iv, ii, nw, sub_ierror
      real(kind=8)               :: th_min, th_max, dth, v_min, v_max, dv
      type(t_wheel),  pointer    :: my_wheel
      real(kind=8),   dimension(:), allocatable :: u_out, v_out

      if (idebug.ge.3) then
         write(bufout,'(/a)') ' --- Start subroutine compute_oor_whl ---'
         call write_log(2, bufout)
      endif

      my_wheel  => ws%whl

      ! form wheel surface whl_curv

      if (.not.my_wheel%prw%is_varprof()) then

         call write_log(' INTERNAL ERROR(oor_whl): needs variable profile')
         call abort_run

      else ! my_wheel%prw%is_varprof()

         if (ic%use_oblique()) then
            call write_log(' Variable wheel (slcw) does not support oblique view direction.')
            call abort_run()
         endif

         ! set search-range on wheel theta; pitch=0 in positive z-direction, counter-clockwise positive

         th_min = -ws%pitch - 0.20d0
         th_max = -ws%pitch + 0.20d0

         ! compute nslc such that dx_ws = O(dx_cntc), dx = r * dth

         nslc = nint( (th_max - th_min) * ws%nom_radius / (4.0d0 * dx_cp) + 2d0 )
         if (mod(nslc,2).eq.0) nslc = nslc + 1
         dth = (th_max - th_min) / max(1d0, real(nslc - 1))

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

         ! nw = my_wheel%prw%grd_data%ntot
         nw    = spl2d%nknotv - 6
         v_min = spl2d%tvj(4)
         v_max = spl2d%tvj(spl2d%nknotv - 3)
         dv  = (v_max - v_min) / real(nw - 1)

         allocate(v_out(nw))
         do iv = 1, nw
            v_out(iv) = v_min + (iv-1) * dv
         enddo
         if (idebug.ge.2) then
            write(bufout,'(a,i4,a,3f12.6)') ' v: v_out([1,2,',nw,'])=',v_out(1), v_out(2), v_out(nw)
            call write_log(1, bufout)
         endif

         ! create storage for curvilinear grid with cylindrical coordinates

         call grid_create_curvil(whl_curv, nslc, nw, cyl_coords=.true.)
         ! call grid_print(whl_curv, 'cyl_grid', 2)

         ! evaluate 2D spline at (u,v) to get (th, y, dr)

         call bspline_eval2d_prod(spl2d, nslc, nw, u_out, v_out, whl_curv%th, whl_curv%y, whl_curv%r,   &
                        .true., sub_ierror, -987d0)

         ! add nominal radius r = rnom + dr

         call grid_shift(whl_curv, 0d0, 0d0, ws%nom_radius)

         ! call grid_print(whl_curv, 'whl_cyl', 5)

         ! convert (th, y, r) to (x, y, z) in wheel center coordinates

         call convert_cyl2cart(whl_curv)

         ! bring -ws%pitch to the lowest point

         call cartgrid_pitch(whl_curv, ws%pitch, 0d0, 0d0)

         ! shift to wheel profile coordinates (z=-rnom at height of axle)

         call grid_shift(whl_curv, 0d0, 0d0, -ws%nom_radius)

         ! optionally create grid-function with (u,v) coordinates of wheel surface points

         if (compute_uv) then
            call gf3_new(whl_uv, 'whl_uv', whl_curv)

            do iv = 1, nw
               do iu = 1, nslc
                  ii = iu + (iv-1) * nslc
                  whl_uv%vx(ii) = u_out(iu)
                  whl_uv%vy(ii) = v_out(iv)
               enddo
            enddo
         endif

         if (.false.) then
            call marker_print( my_wheel%m_ws, 'm_ws', 5)
            call grid_print(whl_curv, 'whl_curv', 5)
            call abort_run()
         endif

         deallocate(u_out, v_out)
         end associate

      endif ! prw%is_varprof()

      if (idebug.ge.4) call write_log('--- end subroutine compute_oor_whl ---')

   end subroutine compute_oor_whl

!------------------------------------------------------------------------------------------------------------

   subroutine compute_wr_locus (ic, ws, trk, prr_rgn, region, ds_cp, sf_whl, sf_rai, idebug,            &
                        has_overlap, my_ierror)
!--purpose: compute the contact locus on the wheel surface, and define profiles sf_whl and sf_rai to
!           be used for locating the interpenetration regions
      implicit none
!--subroutine arguments:
      type(t_ic)                :: ic
      type(t_wheelset)          :: ws
      type(t_trackdata)         :: trk
      type(t_region)            :: region
      real(kind=8)              :: ds_cp
      type(t_grid)              :: prr_rgn, sf_whl, sf_rai
      integer                   :: idebug
      logical,      intent(out) :: has_overlap
      integer,      intent(out) :: my_ierror
!--local variables:
      integer,      parameter :: max_num_kinks = 99
      real(kind=8), parameter :: fac_dyds = 0.3d0
      logical                 :: is_prismatic
      integer                 :: iw, nw, iwsta, iwend, nr, nx, ny, iw_dbg, is_right, is_wheel,          &
                                 nline, nkink, ikinks(max_num_kinks), sub_ierror
      real(kind=8)            :: sgn, xmin_vw, xmax_vw, dx_vw, ymin, ymax, dy, yw_dbg,                  &
                                 kink_high, kink_low, kink_wid, scale_z, ds_thrs
      type(t_marker)          :: rw_trk, rw_vw, rr_vw
      type(t_rail),   pointer :: my_rail
      type(t_wheel),  pointer :: my_wheel
      type(t_grid)            :: prw_rgn, prw_vw, prw, prw_lc, prr_lc, bbr, bbw
      real(kind=8), dimension(:), allocatable :: dy_ds, dz_ds, dz_dy

      my_ierror = 0

      if (idebug.ge.3) then
         write(bufout,'(/a)') ' --- Start subroutine compute_wr_locus ---'
         call write_log(2, bufout)
      endif

      my_rail   => trk%rai
      my_wheel  => ws%whl
      if (ic%is_left_side()) then
         is_right  = 0
      else
         is_right  = 1
      endif

      ! set the sign to -1 for left and +1 for right rail/wheel combination

      sgn = 2 * is_right - 1

      is_prismatic = (.not.ic%is_roller() .and. .not.my_rail%prr%is_varprof())

      ! the (1d) (trimmed) rail profile 'prr_rgn' is given in view coordinates

      nr = prr_rgn%ntot

      ! select search region from (1d) wheel profile 'grd_data'

      if (ic%use_oblique()) then
         call grid_trim(my_wheel%prw%grd_data, prw_rgn, 1, 1, s_low=region%sw_sta, s_hig=region%sw_end, &
                with_spline=.true.)
      else
         call grid_copy(my_wheel%prw%grd_data, prw_rgn, with_spline=.true.)
      endif

      ! the wheel profile 'prw_rgn' is given in wheel profile coordinates
      ! get marker to convert wheel profile to view coordinates

      rw_trk = marker_2glob( my_wheel%m_ws, ws%m_trk )
      rw_vw  = marker_2loc( rw_trk, region%mview )
      if (idebug.ge.3) call marker_print(rw_vw, 'rw_vw', 5)

      ! convert prw to view-coordinates, determine fold-back there (vertical sections + roll angle)

      call grid_copy(prw_rgn, prw_vw, with_spline=.true.)
      call cartgrid_2glob(prw_vw, rw_vw)
      ! call grid_print(prw_vw, 'prw_vw', 5, 8)

      ! call write_log('prw_determine_foldback')
      call prw_determine_foldback(prw_vw, my_wheel%prw%f_max_omit, my_wheel%prw%err_hnd, iwsta, iwend,  &
                idebug, sub_ierror)
      if (my_ierror.eq.0) my_ierror = sub_ierror
      call grid_destroy(prw_vw)

      ! create trimmed version of prw_rgn (wheel coords) without fold-back when viewed in view coords

      call grid_trim(prw_rgn, prw, 1, 1, iwsta, iwend)
      nw = prw%ntot
      ! call grid_print(prw, 'prw_trim', 5, 8)

      ! copy trimmed prw without the spline data, the xw and zw-coordinates will be adapted

      call grid_copy(prw, prw_lc)
      prw_lc%lies_in_oyz = .false.

      ! select point for debugging

      iw_dbg = -143
      yw_dbg = 0d0
      if (.false.) then
         do iw = 1, nw
            if (abs(prw%y(iw)-yw_dbg).lt.0.1d0) iw_dbg = iw
         enddo
      endif
      if (iw_dbg.gt.nw) iw_dbg = -1
      if (iw_dbg.ge.1) call marker_print(rw_vw, 'rw(vw)', 5)

      ! evaluate profile slope dz/dy at s-positions of wheel profile prw
      !    (supporting vertical sections in wheel coordinates)

      allocate(dz_dy(nw))
      if (.true.) then
         allocate(dy_ds(nw), dz_ds(nw))
         call spline_eval(prw%spl, ikYDIR, nw, prw%s_prf, sub_ierror, f1_eval=dy_ds)
         call spline_eval(prw%spl, ikZDIR, nw, prw%s_prf, sub_ierror, f1_eval=dz_ds)
         do iw = 1, nw
            dz_dy(iw) = dz_ds(iw) / min(-1d-6, dy_ds(iw))
            if (iw.eq.iw_dbg) then
               write(bufout,'(a,i4,4(a,g12.4))') ' wheel iw=',iw,': yw=',prw%y(iw),', dy_ds=',dy_ds(iw), &
                     ', dz_ds=',dz_ds(iw),': dz_dy=',dz_dy(iw)
               call write_log(1, bufout)
            endif
         enddo
         deallocate(dy_ds, dz_ds)
      else
         call spline_get_dzdy_at_s( prw%spl, nw, prw%s_prf, dz_dy, sub_ierror )
      endif

      ! wheel-on-rail configurations: straight-forward calculation

      if (is_prismatic) then

         call locus_prismatic( ws, rw_vw, prw, nw, dz_dy, prw_lc, iw_dbg, idebug )

      elseif (my_rail%prr%is_varprof()) then

         ! wheel-on-variable profile: iterative calculation

         rr_vw  = marker_2loc( my_rail%m_trk, region%mview )
         call locus_iterate_varprof( my_rail%prr, ws, rr_vw, rw_vw, prw, nw, dz_dy, prr_lc, prw_lc,     &
                iw_dbg, idebug )

      elseif (ic%is_roller()) then

         ! wheel-on-roller: iterative calculation

         call locus_iterate_roller( trk, ws, rw_vw, prr_rgn, prw, nw, dz_dy, prr_lc, prw_lc, iw_dbg,    &
                idebug )

      else

         call write_log(' Internal error (wr_locus): not prismatic/varprof/roller')
         call abort_run()

      endif ! wheel-on-rail

      deallocate(dz_dy)

      ! call write_log('grid_check_nan prw_lc')
      if (idebug.ge.1) call grid_check_nan(prw_lc, 'locus(rw)', idebug)
      if (idebug.ge.4) call grid_print(prw_lc, 'locus_rw', 5)

      ! convert the wheel mesh in one step from profile to view coordinates

      ! call write_log('cartgrid_2glob prw_lc')
      call cartgrid_2glob(prw_lc, rw_vw)

      if (idebug.ge.1) call grid_check_nan(prw_lc, 'locus(vw)', idebug)
      if (idebug.ge.4) call grid_print(prw_lc, 'locus_vw', 5)

      ! detect sudden steps in contact locus, including overall start/end-points

      is_wheel  = -1
      kink_high = pi/10d0
      kink_low  = pi/10d0
      kink_wid  = 0.01d0
      scale_z   = 0.1d0
      ds_thrs   = 100d0

      call profile_find_kinks(nw, prw_lc%y, prw_lc%x, is_wheel, kink_high, kink_low, kink_wid,          &
                        scale_z, nkink, ikinks, idebug-1, sub_ierror, s=prw_lc%s_prf, ds_thrs=ds_thrs)

      if (idebug.ge.2 .and. nkink.gt.2) then
         write(bufout,'(a,i3,a,20i5,:,/,4(10x,20i5,:,/))') ' contact locus: nkink=',nkink,', ikinks= ', &
                ikinks(1:min(nkink,max_num_kinks))
         nline = int( (min(nkink,max_num_kinks)-1)/20 ) + 1
         call write_log(nline, bufout)
      endif

      ! create interpolating spline for the contact locus, preparing for interpolations

      call grid_make_ppspline(prw_lc, 0d0, .false., nkink, ikinks, my_ierror=sub_ierror)
      ! call grid_print(prw_lc, 'locus_vw', 5)

      if (.not.is_prismatic) then
         call grid_make_arclength(prr_lc, sub_ierror)

         call profile_find_kinks(nw, prr_lc%y, prr_lc%x, is_wheel, kink_high, kink_low, kink_wid,       &
                        scale_z, nkink, ikinks, idebug-1, sub_ierror, s=prr_lc%s_prf, ds_thrs=ds_thrs)

         if (idebug.ge.2 .and. nkink.gt.2) then
            write(bufout,'(a,i3,a,20i5,/,4(10x,20i5,:,/))') ' contact locus: nkink=',nkink,', ikinks=', &
                   ikinks(1:min(max_num_kinks,nkink))
            nline = int( (min(nkink,max_num_kinks)-1)/20 ) + 1
            call write_log(nline, bufout)
         endif

         call grid_make_ppspline(prr_lc, 0d0, .false., nkink, ikinks, my_ierror=sub_ierror)
         ! call grid_print(prr_lc, 'prr_lc', 5)
      endif

      ! get range of x-coordinates of contact locus on wheel

      ! call write_log('grid_get_xrange prw_lc')
      call grid_get_xrange(prw_lc, xmin_vw, xmax_vw)
      dx_vw = 0d0

      if (idebug.ge.2) then
         write(bufout,'(3(a,f10.3))') ' locus(vw): x in [', xmin_vw,',', xmax_vw,'], dx=',dx_vw
         call write_log(1, bufout)
      endif

      ! define a uniform grid 'gap_mesh' in view coordinates, for computing the gap function

      ! - y-range: considering all y of the rail profile
      !            fac_dyds=0.5 gives steps of size ds along the surface where surface inclined by 60deg
      !            the resulting grid may be a little smaller at ymax

      call grid_get_yrange(prr_rgn, ymin, ymax)
      dy   = fac_dyds * ds_cp

      call grid_create_uniform(sf_rai, nxarg=1, x0arg=0d0, x1arg=0d0,                                 &
                                         y0arg=ymin, dyarg=dy, y1arg=ymax, zarg=0d0)
      call grid_create_uniform(sf_whl, nxarg=1, x0arg=0d0, x1arg=0d0,                                 &
                                         y0arg=ymin, dyarg=dy, y1arg=ymax, zarg=0d0)

      call grid_get_dimens(sf_rai, nx, ny)

      ! interpolate rail profile to 1d uniform 'sf_rai' profile, locus prw_lc to 'sf_whl' profile

      if (is_prismatic) then

         call spline_get_xz_at_y(prr_rgn%spl, sf_rai, sub_ierror)

      else

         call spline_get_xz_at_y(prr_lc%spl, sf_rai, sub_ierror)

      endif

      call spline_get_xz_at_y(prw_lc%spl, sf_whl, sub_ierror, -888d0)

      ! Check that the wheel and rail surfaces overlap

      has_overlap = .true.
      call grid_get_boundbox(prw_lc, bbw)
      call grid_get_boundbox(sf_rai, bbr)

      ! Extend bounding box to support contact for two flat profiles

      if (bbw%z(8)-bbw%z(1).lt.0.001) then
         if (idebug.ge.2) call write_log(' wheel: flat profile')
         bbw%z(1) = bbw%z(1) - 0.05
      endif
      if (bbr%z(8)-bbr%z(1).lt.0.001) then
         if (idebug.ge.2) call write_log(' rail: flat profile')
         bbr%z(8) = bbr%z(8) + 0.05
      endif

      if (bbw%y(1).gt.bbr%y(8) .or. bbr%y(1).gt.bbw%y(8)) then
         has_overlap = .false.
         write(bufout,'(2a,/,2(2(a,f8.3),a,:,/))') ' WARNING: the wheel and rail are laterally in ',    &
                        'different places.',                                                            &
                        '          the wheel grid has y=[', bbw%y(1),',', bbw%y(8),'],',                &
                        '          the  rail grid has y=[', bbr%y(1),',', bbr%y(8),'].'
         call write_log(3, bufout)
      endif

      if (ic%norm.le.0 .and. minval(bbw%z).gt.maxval(bbr%z)) then
         has_overlap = .false.
         if (maxval(abs(bbw%z)).gt.9999. .or. maxval(abs(bbr%z)).gt.9999.) then
            write(bufout,'(2a,/,2(5a,:,/))') ' WARNING: the wheel and rail are vertically in ',         &
                        'different places.',                                                            &
                        '          the wheel grid has z=[', fmt_gs(8,3,minval(bbw%z)),',',              &
                                                            fmt_gs(8,3,maxval(bbw%z)),'],',             &
                        '          the  rail grid has z=[', fmt_gs(8,3,minval(bbr%z)),',',              &
                                                            fmt_gs(8,3,maxval(bbr%z)),'].'
         else
            write(bufout,'(2a,/,2(2(a,f8.3),a,:,/))') ' WARNING: the wheel and rail are vertically ',   &
                        'in different places.',                                                         &
                        '          the wheel grid has z=[', minval(bbw%z),',', maxval(bbw%z),'],',      &
                        '          the  rail grid has z=[', minval(bbr%z),',', maxval(bbr%z),'].'
         endif
         call write_log(3, bufout)
      endif

      ! no check needed on bbw%x versus bbr%x -- prw_lc describes the x-values for both wheel and rail

      ! for debugging: bounding-boxes of prw_lc, sf_whl and sf_rai

      if (idebug.ge.3 .or. (idebug.ge.1 .and. .not.has_overlap)) then
         write(bufout,248) 'Grid locus(vw)', bbw%x(1), bbw%x(8), bbw%y(1), bbw%y(8), bbw%z(1), bbw%z(8)
         call write_log(2, bufout)

         call grid_get_boundbox(sf_whl, bbw)
         write(bufout,248) 'Grid sf_whl', bbw%x(1), bbw%x(8), bbw%y(1), bbw%y(8), bbw%z(1), bbw%z(8)
         call write_log(2, bufout)

         write(bufout,248) 'Grid sf_rai', bbr%x(1), bbr%x(8), bbr%y(1), bbr%y(8), bbr%z(1), bbr%z(8)
         call write_log(2, bufout)

  248    format(1x,a,' has bounding box',/,                                                             &
                '      x=[',f8.3,',',f8.3,'], y=[',f8.3,',',f8.3,'], z=[',f8.3,',',f8.3,']')
      endif

      if (idebug.ge.4) call grid_print(prw_lc, 'prw_lc', 5)
      if (idebug.ge.8) call grid_print(sf_whl, 'sf_whl', 5)
      if (idebug.ge.8) call grid_print(sf_rai, 'sf_rai', 5)

      ! cleanup local variables

      call grid_destroy(prw_rgn)
      call grid_destroy(prw_vw)
      call grid_destroy(prw)
      call grid_destroy(prw_lc)
      call grid_destroy(prr_lc)
      call grid_destroy(bbr)
      call grid_destroy(bbw)

      if (idebug.ge.4) call write_log('--- end subroutine compute_wr_locus ---')

   end subroutine compute_wr_locus

!------------------------------------------------------------------------------------------------------------

   subroutine locus_prismatic( ws, rw_vw, prw, nw, dz_dy, prw_lc, iw_dbg, idebug )
!--purpose: 
      implicit none
!--subroutine arguments:
      type(t_wheelset)             :: ws
      type(t_marker)               :: rw_vw
      type(t_grid)                 :: prw, prw_lc
      integer,         intent(in)  :: nw, iw_dbg, idebug
      real(kind=8),    intent(in)  :: dz_dy(nw)
!--local variables:
      integer      :: iw
      real(kind=8) :: ry

      associate( r => rw_vw%rot )

      if (iw_dbg.ge.1) then
         write(bufout,*) 'printing data for yw=',prw_lc%y(iw_dbg)
         call write_log(1, bufout)
      endif

      ! locate x_w positions on the wheel surface with maximum z_r value

      do iw = 1, nw

         ! TODO: document using r%el instead of tan(psi)

         ry           = ws%nom_radius + prw%z(iw)
         prw_lc%x(iw) = (ry * dz_dy(iw) * r%el(1,2) - ws%nom_radius * r%el(1,3)) / r%el(1,1)

         ! avoid too large x-values at near vertical sections, creating distorted (multi-valued) y_vw

         prw_lc%x(iw) = max(-0.5d0*ws%nom_radius, min(0.5d0*ws%nom_radius, prw_lc%x(iw)))

         ! compute z-value at given x-value

         if (ry**2-prw_lc%x(iw)**2.gt.0d0) then
            prw_lc%z(iw) = sqrt( ry**2 - prw_lc%x(iw)**2 ) - ws%nom_radius
         else
            prw_lc%z(iw) = -999d0
         endif

         if (iw.eq.iw_dbg) then
            write(bufout,'(a,i4,6(a,g12.4))') ' iw=',iw,': ry=',ry, ', dz_dy=',dz_dy(iw), ', x=',    &
                 prw_lc%x(iw),' =(', ry*dz_dy(iw)*r%el(1,2), ' -', ws%nom_radius*r%el(1,3),') /', r%el(1,1)
            call write_log(1, bufout)
         endif

      enddo

      end associate

      if (idebug.ge.4) call write_log('--- end subroutine locus_prismatic ---')

   end subroutine locus_prismatic

!------------------------------------------------------------------------------------------------------------

   subroutine locus_iterate_varprof( vprf, ws, rr_vw, rw_vw, prw, nw, dz_dy, prr_lc, prw_lc, iw_dbg,    &
                idebug )
!--purpose: 
      implicit none
!--subroutine arguments:
      type(t_profile)              :: vprf
      type(t_wheelset)             :: ws
      type(t_marker)               :: rr_vw, rw_vw
      type(t_grid)                 :: prw, prr_lc, prw_lc
      integer,         intent(in)  :: nw, iw_dbg, idebug
      real(kind=8),    intent(in)  :: dz_dy(nw)
!--local variables:
      real(kind=8), parameter :: f_bisec = 0.5d0
      logical      ::    lstop
      integer            :: iter, maxit, iw, iw_max, j, nout, sub_ierror
      real(kind=8)       :: ry, rnom_z, a_y, denom, x_new, dx0, dx1, dx_max, d2zdx2, epslc, exterval
      real(kind=8), dimension(:), allocatable :: dzdx_vprf, x_lbnd, x_ubnd, dzdx_lbnd, dzdx_ubnd
      type(t_grid)       :: rsurf
      character(len=10)  :: tmpnam

      ! wheel-on-varprof: iterative calculation of contact locus

      associate( r => rw_vw%rot )      ! orientation / rotation matrix for wheel in view system
      allocate(dzdx_vprf(nw), x_lbnd(nw), x_ubnd(nw), dzdx_lbnd(nw), dzdx_ubnd(nw))

      ! evaluate the rail profile for an intpol.path in x-direction

      if (iw_dbg.ge.1 .and. iw_dbg.le.nw .and. .false.) then
      ! if (.false.) then

         nout = 1001
         if (.false.) then
            ! create grid in (s_fc,y_r) coordinates
            call grid_create_uniform(rsurf, nxarg=nout, x0arg=-100d0, x1arg=100d0, nyarg=1,             &
                        y0arg=40d0, dyarg=1d0)
         else
            ! create grid in (x_vw,y_w) coordinates
            call grid_create_uniform(rsurf, nxarg=nout, x0arg=-250d0, x1arg=250d0, nyarg=1,             &
                        y0arg=prw_lc%y(iw_dbg), dyarg=1d0)
            call cartgrid_2glob(rsurf, rw_vw)
            call cartgrid_2loc(rsurf, rr_vw)
            rsurf%x(1:nout) = rsurf%x(1:nout) + ws%s
         endif

         call bspline_get_z_at_xy(vprf%spl2d, rsurf%ntot, rsurf%x, rsurf%y, rsurf%z, sub_ierror, 888d0)

         call grid_print(rsurf, 'xslc', 5)
         if (.true.) then
            write(bufout,'(a,f8.3,a)') 's_ws =', ws%s, ';'
            call write_log(1, bufout)
            call write_log('x = xslc(:,1) - s_ws; dx = diff(x); xm = (x(1:end-1)+x(2:end))/2;')
            call write_log('y = xslc(:,2); z = xslc(:,3); dz = diff(z);')
         endif
         call abort_run()

      endif

      ! TODO: wheel-on-rail calculation can be used as initial guess:
      !       locate x_w positions on the wheel surface with maximum z_w value

      ! set bracket - lower and upper bounds on x for each grid position based on availability of varprof
      ! Note: avoid too large x-values (> 0.6 rnom), creating distorted (multi-valued) y_vw

      call grid_copy(prw_lc, prr_lc)
      call cartgrid_2glob(prr_lc, rw_vw)
      call cartgrid_2loc(prr_lc, rr_vw)
      call varprof_available_x( vprf, ws%s, nw, prr_lc%y, 0.6d0*ws%nom_radius, x_lbnd, x_ubnd)
      do iw = 1, nw
         if (x_lbnd(iw).gt.x_ubnd(iw)) then
            x_lbnd(iw) = 0d0
            x_ubnd(iw) = 0d0
         endif

         if (idebug.ge.5) then
            write(bufout,'(a,i4,3(a,f12.3),a)') ' y(',iw,')=',prr_lc%y(iw),': x in [',x_lbnd(iw),',',   &
                   x_ubnd(iw),']'
            call write_log(1, bufout)
         endif
      enddo

      dzdx_vprf(1:nw) = 0d0

      ! Perform iterations to compute the contact locus

      maxit = 150
      epslc = 1d-3
      lstop = .false.
      iter  = 0

      do while(.not.lstop .and. iter.lt.maxit)
         iter = iter + 1

         if (iw_dbg.ge.1 .and. idebug.ge.3) then
            write(bufout,'(a,i3,a,f9.4)') ' starting iter=',iter,', printing data for yw=',prw_lc%y(iw_dbg)
            call write_log(1, bufout)
         endif

         ! copy current prw_lc to prr_lc, convert prr_lc to view-coordinates to get the lc positions

         call grid_copy(prw_lc, prr_lc)
         call cartgrid_2glob(prr_lc, rw_vw)

         if (.false. .and. iw_dbg.ge.1) then
            do j = -4, 0
               write(bufout,124) 'prw_lc(',j,',rw)=[',prw_lc%x(iw_dbg+j),',', prw_lc%y(iw_dbg+j),',',   &
                        prw_lc%z(iw_dbg+j),']'
               call write_log(1, bufout)
               write(bufout,124) 'prw_lc(',j,',tr)=[',prr_lc%x(iw_dbg+j),',', prr_lc%y(iw_dbg+j),',',   &
                        prr_lc%z(iw_dbg+j),']'
            call write_log(1, bufout)
            enddo
 124        format(1x,a,i3,3(a,f12.6),a)
         endif

         ! compute varprof surface at (xlc, ylc) positions
         ! compute slope at x-positions on varprof given in prr_lc

         call cartgrid_2loc(prr_lc, rr_vw)
         ! if (iter.eq.4) call varprof_set_debug(3, iw_dbg)
         ! if (iter.eq.4) call bspline2d_print(vprf%spl2d, 'spl2d', 5)
         ! call bspline_set_debug(2, iw_dbg)
         exterval = 987d0
         call varprof_get_z_at_xy( vprf, ws%s, nw, prr_lc%x, prr_lc%y, prr_lc%z, dzdx_vprf,             &
                                   sub_ierror, exterval)
         ! call bspline_set_debug(0)
         ! call varprof_set_debug(0)
         if (sub_ierror.ne.0) then
            write(bufout,*) 'locus_iterate_varprof: error',sub_ierror,' in filling prr_lc z-data'
            call write_log(1, bufout)
         endif
         call cartgrid_2glob(prr_lc, rr_vw)

         if (iw_dbg.ge.1 .and. idebug.ge.2) then
            write(bufout,'(2(a,f12.6))') ' x_lc  =',prw_lc%x(iw_dbg),', dzdx_vprf=',dzdx_vprf(iw_dbg)
            call write_log(1, bufout)
         endif

         ! estimate new x-positions for locus in wheel-coordinates

         iw_max = 1
         dx_max = 0d0

         do iw = 1, nw

            ry     = ws%nom_radius + prw%z(iw)

            ! compute z_lc at x_lc

            if (ry**2-prw_lc%x(iw)**2.gt.0d0) then
               prw_lc%z(iw) = sqrt( ry**2 - prw_lc%x(iw)**2 ) - ws%nom_radius
            else
               prw_lc%z(iw) = -999d0
            endif

            rnom_z = ws%nom_radius + prw_lc%z(iw)

            ! compute denominator using 3rd column of R^T, 3rd row of R

            a_y   = ry * dz_dy(iw)
            denom = prw_lc%x(iw) * r%el(3,1) - a_y * r%el(3,2) + rnom_z * r%el(3,3)

            ! compute estimate for d2z_dx2

            dx0 = x_ubnd(iw) - x_lbnd(iw)
            if (dx0.le.epslc .or. dx0.gt.0.1d0*ws%nom_radius) then
               d2zdx2 = 0d0
            else
               d2zdx2 = ( dzdx_ubnd(iw) - dzdx_lbnd(iw) ) / dx0
               d2zdx2 = max(-0.5d0/denom, d2zdx2)
               d2zdx2 = min( 0.5d0/denom, d2zdx2)
               if (iw.eq.iw_dbg .and. idebug.ge.3) then
                  write(bufout,'(5(a,f12.6))') ' d2zdx2 = (', dzdx_ubnd(iw),' - ',dzdx_lbnd(iw),        &
                        ' ) / (',x_ubnd(iw),' - ',x_lbnd(iw),' ) =', d2zdx2
                  call write_log(1, bufout)
               endif
            endif

            x_new = ( a_y * r%el(1,2) - denom * dzdx_vprf(iw) + denom * d2zdx2 * prw_lc%x(iw) ) /       &
                                                                        ( r%el(1,1) + denom * d2zdx2 )

            if (iw.eq.iw_dbg .and. idebug.ge.2) then
               write(bufout,'(2(a,i4),7(a,f11.6))') ' it=',iter,', iw=',iw,': denom=  ',denom,' =  ',   &
                    prw_lc%x(iw),' * ',r%el(3,1),' - ',a_y,' * ',r%el(3,2),' + ',rnom_z,' * ',r%el(3,3)
               call write_log(1, bufout)
               write(bufout,'(2(a,i4),5(a,f11.6),a,f11.3,2(a,f11.6),a)') ' it=',iter,', iw=',iw,        &
                    ': x_new=  ', x_new,' = (', a_y,' * ', r%el(1,2),' - ', denom,' * ', dzdx_vprf(iw), &
                    ' + ', denom * d2zdx2 * prw_lc%x(iw), ') / ( ', r%el(1,1), ' + ', denom * d2zdx2, ' )'
               call write_log(1, bufout)
            endif

            dx1 = x_new - prw_lc%x(iw)
            if (dx1.lt.0d0) then
               dx0 = x_lbnd(iw) - prw_lc%x(iw)
            else
               dx0 = x_ubnd(iw) - prw_lc%x(iw)
            endif

            if (iw.eq.iw_dbg) then
               write(bufout,'(2(a,i4),4(a,f12.6))') ' it=',iter,', iw=',iw,', x_prev=', prw_lc%x(iw),   &
                        ', x_new=',x_new,', lbnd=',x_lbnd(iw),', ubnd=',x_ubnd(iw)
               call write_log(1, bufout)
               if (dx1.gt.0d0) then
                  write(bufout,'(42x,a,f12.6,a,21x,2(a,f12.6))') 'dx1=',dx1,',','dx0=',dx0,', ratio=',dx1/dx0
               else
                  write(bufout,'(42x,a,f12.6,2(a,f12.6))') 'dx1=',dx1,',  dx0=',dx0,          &
                                                                     ',                    ratio=',dx1/dx0
               endif
               call write_log(1, bufout)
            endif

            ! use bisection step when x_new jumps out of the bracket [lbnd, ubnd],
            !  or when jumping from one side of the bracket to the other

            if (x_new.lt.x_lbnd(iw)) then
               x_new = f_bisec * x_lbnd(iw) + (1d0 - f_bisec) * prw_lc%x(iw)
               if (iw.eq.iw_dbg) then
                  write(bufout,'(2(a,f12.6))') ' x_new < lbnd =',x_lbnd(iw),', bisection x_new =', x_new
                  call write_log(1, bufout)
               endif
            elseif (x_new.gt.x_ubnd(iw)) then
               x_new = f_bisec * x_ubnd(iw) + (1d0 - f_bisec) * prw_lc%x(iw)
               if (iw.eq.iw_dbg) then
                  write(bufout,'(2(a,f12.6))') ' x_new > ubnd =',x_ubnd(iw),', bisection x_new =', x_new
                  call write_log(1, bufout)
               endif
            else
               if (dx1/dx0.gt.0.7d0) then
                  x_new = prw_lc%x(iw) + 0.7d0 * dx0
                  if (iw.eq.iw_dbg .and. idebug.ge.3) then
                     write(bufout,'(30x,a,f12.6)') ' clipped x_new =', x_new
                     call write_log(1, bufout)
                  endif
               endif
            endif

            ! update the bracket [lbnd, ubnd]

            if (x_new.gt.prw_lc%x(iw)) then
               x_lbnd(iw)    = prw_lc%x(iw)
               dzdx_lbnd(iw) = dzdx_vprf(iw)
            endif
            if (x_new.lt.prw_lc%x(iw)) then
               x_ubnd(iw)    = prw_lc%x(iw)
               dzdx_ubnd(iw) = dzdx_vprf(iw)
            endif

            ! avoid too large x-values, creating distorted (multi-valued) y_vw

            x_new = max(-0.5d0*ws%nom_radius, min(0.5d0*ws%nom_radius, x_new))

            ! determine maximum update

            if (abs(x_new - prw_lc%x(iw)).gt.dx_max) then
               iw_max = iw
               dx_max = abs(x_new - prw_lc%x(iw))
            endif

            prw_lc%x(iw) = x_new

         enddo ! iw

         if (idebug.ge.2) then
            write(bufout,'(a,i3,a,f12.6,a,i4)') ' locus_varprof: it=',iter,': max(dx)=',dx_max,         &
                ' at iw=',iw_max
            call write_log(1, bufout)
         endif

         lstop = (dx_max.lt.epslc)

         if (idebug.ge.2 .and. .false.) then
            if (iter.le.9) then
               write(tmpnam,'(a,i1)') 'prr_lc', iter
            elseif (iter.le.99) then
               write(tmpnam,'(a,i2)') 'prr_lc', iter
            else
               write(tmpnam,'(a,i3)') 'prr_lc', iter
            endif
            call grid_print(prw_lc, tmpnam, 5)
         endif

      enddo ! while(not lstop)

      if (idebug.ge.2) then
         write(bufout,'(a,i3,a,f12.6,a)') ' locus_varprof: ',iter,' iterations, residual=',dx_max,' mm'
         call write_log(1, bufout)
      endif

      call grid_destroy(rsurf)
      deallocate(dzdx_vprf, x_lbnd, x_ubnd, dzdx_lbnd, dzdx_ubnd)
      end associate

      if (idebug.ge.4) call write_log('--- end subroutine locus_iterate_varprof ---')

   end subroutine locus_iterate_varprof

!------------------------------------------------------------------------------------------------------------

   subroutine locus_iterate_roller( trk, ws, rw_vw, prr, prw, nw, dz_dy, prr_lc, prw_lc, iw_dbg, idebug )
!--purpose: 
      implicit none
!--subroutine arguments:
      type(t_trackdata)            :: trk
      type(t_wheelset)             :: ws
      type(t_marker)               :: rw_vw
      type(t_grid)                 :: prr, prw, prr_lc, prw_lc
      integer,         intent(in)  :: nw, iw_dbg, idebug
      real(kind=8),    intent(in)  :: dz_dy(nw)
!--local variables:
      integer      :: iter, niter, iw, j, sub_ierror
      real(kind=8) :: sgn_rol, ry, rnom_z, a_y, denom, lhs
      real(kind=8), dimension(:), allocatable :: dzdx_rol, rrol

      ! wheel-on-roller: iterative calculation of contact locus

      associate( r => rw_vw%rot )
      allocate(dzdx_rol(nw), rrol(nw))

      ! TODO: wheel-on-rail calculation can be used as initial guess:
      !       locate x_w positions on the wheel surface with maximum z_w value

      ! Perform 3 iterations to compute the contact locus

      niter = 3
      do iter = 1, niter

         if (iw_dbg.ge.1) then
            write(bufout,'(a,i2,a,f9.4)') ' starting iter=',iter,', printing data for yw=',prw_lc%y(iw_dbg)
            call write_log(1, bufout)
         endif

         ! copy current prw_lc to prr_lc, convert prr_lc to view-coordinates to get the lc positions

         call grid_copy(prw_lc, prr_lc)
         call cartgrid_2glob(prr_lc, rw_vw)

         if (iw_dbg.ge.1) then
            do j = -4, 0
               write(bufout,124) 'prw_lc(',j,',rw)=[',prw_lc%x(iw_dbg+j),',', prw_lc%y(iw_dbg+j),',',   &
                        prw_lc%z(iw_dbg+j),']'
               call write_log(1, bufout)
               write(bufout,124) 'prw_lc(',j,',tr)=[',prr_lc%x(iw_dbg+j),',', prr_lc%y(iw_dbg+j),',',   &
                        prr_lc%z(iw_dbg+j),']'
            call write_log(1, bufout)
            enddo
 124        format(1x,a,i3,3(a,f12.6),a)
         endif

         ! interpolate roller profile z-values at locus y-positions, computing local radii rrol(y)

         if (iter.eq.-2) call spline_set_debug(3)
         call spline_get_xz_at_y( prr%spl, nw, prr_lc%y, sub_ierror, zout=rrol )
         if (sub_ierror.ne.0) then
            write(bufout,*) 'wr_contact_locus: error',sub_ierror,' in filling spline prr_lc xz-data'
            call write_log(1, bufout)
         endif
         call spline_set_debug(0)

         sgn_rol = 1d0       ! sign of m_rol z-position, negative for a raceway
         if (trk%nom_radius.lt.0d0) sgn_rol = -1d0

         do iw = 1, nw
            rrol(iw) = sgn_rol * (trk%nom_radius - rrol(iw))
         enddo

         ! compute roller surface at (xlc, ylc) positions

         do iw = 1, nw
            if (rrol(iw)**2-prr_lc%x(iw)**2.gt.0d0) then
               prr_lc%z(iw) = trk%nom_radius - sgn_rol * sqrt( rrol(iw)**2 - prr_lc%x(iw)**2 )
            else
               prr_lc%z(iw) = 888d0
            endif
         enddo

         ! compute slope at x-positions on roller given in prr_lc

         do iw = 1, nw
            dzdx_rol(iw) = sgn_rol * prr_lc%x(iw) / rrol(iw)
         enddo

         if (iw_dbg.ge.1) then
            write(bufout,'(a,f12.4,a,g12.4)') ' rrol=',rrol(iw_dbg),', dzdx_rol=',dzdx_rol(iw_dbg)
            call write_log(1, bufout)
         endif

         ! estimate new x-positions for locus in wheel-coordinates

         do iw = 1, nw

            ry     = ws%nom_radius + prw%z(iw)

            ! compute z_lc at x_lc

            if (ry**2-prw_lc%x(iw)**2.gt.0d0) then
               prw_lc%z(iw) = sqrt( ry**2 - prw_lc%x(iw)**2 ) - ws%nom_radius
            else
               prw_lc%z(iw) = -999d0
            endif

            rnom_z = ws%nom_radius + prw_lc%z(iw)

            ! compute denominator using 3rd column of R^T, 3rd row of R

            a_y   = ry * dz_dy(iw)
            denom = prw_lc%x(iw) * r%el(3,1) - a_y * r%el(3,2) + rnom_z * r%el(3,3)

            prw_lc%x(iw) = ( - denom * dzdx_rol(iw)                                                  &
                             + denom * prw_lc%x(iw) / trk%nom_radius                                 &
                             + a_y * r%el(1,2)                                                       &
                             - rnom_z * r%el(1,3)                                                    &
                           ) / (r%el(1,1) + denom / trk%nom_radius)

            if (iw.eq.iw_dbg) then
               write(bufout,'(2(a,i4),2(a,f10.3,a,f10.6))') ' it=',iter,', iw=',iw,': r(yw)=', ry,   &
                                ', zw=', prw_lc%z(iw), ', denom=',denom,', new xlc=', prw_lc%x(iw)
               call write_log(1, bufout)
               lhs = (prw_lc%x(iw)*r%el(1,1) - ry*dz_dy(iw)*r%el(1,2)) / denom
               write(bufout,'(a,f12.4,2(a,f10.6))') ' ay=',ry*dz_dy(iw),', lhs=',lhs,', rhs=',       &
                     prr_lc%x(iw)/rrol(iw)
               call write_log(1, bufout)
            endif

            ! avoid too large x-values, creating distorted (multi-valued) y_vw

            prw_lc%x(iw) = max(-0.5d0*ws%nom_radius, min(0.5d0*ws%nom_radius, prw_lc%x(iw)))

         enddo ! iw

      enddo ! iter

      deallocate(dzdx_rol, rrol)
      end associate

      if (idebug.ge.4) call write_log('--- end subroutine locus_iterate_roller ---')

   end subroutine locus_iterate_roller

!------------------------------------------------------------------------------------------------------------

   subroutine make_surfc_inclination( prr, g, sf_incl, idebug )
!--purpose: estimate the surface inclination angles of rail profile prr at y-locations of grid g (sf_rai)
      implicit none
!--subroutine arguments:
      type(t_grid),    intent(in)  :: prr, g
      integer,         intent(in)  :: idebug
      type(t_gridfnc3)             :: sf_incl
!--local variables:
      integer      :: sub_ierror, ix, iy, ii, nx, ny
      real(kind=8), dimension(:), allocatable :: y, alph

      ! allocate work-arrays and grid function

      nx = g%nx
      ny = g%ny
      allocate(y(ny), alph(ny))

      call gf3_new(sf_incl, 'sf_incl', g, nulify=.true.)

      ! get y-positions from grid g = sf_rai. 2d grid: using first column

      do iy = 1, ny
         ii = (iy-1) * nx + 1
         y(iy) = g%y(ii)
      enddo

      ! evaluate profile inclination atan2(dz, dy) at y-positions used in grid g

      call spline_get_alpha_at_y( prr%spl, ny, y, alph, sub_ierror )

      ! store surface inclination. 2d grid: copy to all columns

      do iy = 1, ny
         do ix = 1, nx
            ii = (iy-1) * nx + ix
            sf_incl%vy(ii) = alph(iy)
         enddo
      enddo

      if (idebug.ge.4) call gf3_print(sf_incl, 'sf_incl', ikYDIR, 4)

      ! clean up

      deallocate(y, alph)

   end subroutine make_surfc_inclination

!------------------------------------------------------------------------------------------------------------

   subroutine locate_interpen_1d( ic, sf_whl, sf_rai, sf_incl, m_vw, sgn, nom_radius, roller_radius,    &
                                  gap_min, delt_min, zw_min, a1, b1, numcps, allcps, idebug)
!--purpose: find the points in the contact locus on the wheel with (locally) maximum interpenetration,
!           define the corresponding interpenetration regions
      implicit none
!--subroutine arguments:
      type(t_ic)                  :: ic
      type(t_grid)                :: sf_whl, sf_rai      ! wheel contact locus & rail profile at uniform dy
                                                         ! in view coordinates
      type(t_gridfnc3)            :: sf_incl             ! rail surface inclination wrt view coordinates
                                                         ! at y-pos of sf_rai
      type(t_marker)              :: m_vw                ! view-marker in track coordinates
      real(kind=8), intent(in)    :: sgn                 ! +1/-1 for right/left w/r-combination
      real(kind=8), intent(in)    :: nom_radius          ! nominal wheel radius
      real(kind=8), intent(in)    :: roller_radius       ! nominal roller radius if >0, 0 for tracks
      real(kind=8), intent(out)   :: gap_min             ! overall minimum gap found
      real(kind=8), intent(out)   :: delt_min            ! contact angle at the overall minimum gap value
      real(kind=8), intent(out)   :: zw_min              ! wheel surface height at overall minimum gap value
      real(kind=8), intent(out)   :: a1, b1              ! curvatures at the overall minimum gap value
      integer,      intent(inout) :: numcps              ! number of contact problems
      type(p_cpatch)              :: allcps(MAX_NUM_CPS) ! data of contact problems
      integer,      intent(in)    :: idebug
!--local variables:
      integer,      parameter :: MAX_NUM_MIN = 1000
      logical        :: has_zero
      integer        :: ix, iy, iy_sta, iy_end, iy_min, nguard, numcp0
      integer        :: ii_n, ii_s, ilm, ixmin, nummin, numrem
      integer        :: locmin(2,MAX_NUM_MIN)
      real(kind=8)   :: gap_locmin, xmin, ymin, gmin, zrmin, curv_y, xmax, dgap, ysta, yend, zsta, zend
      type(t_cpatch),             pointer     :: cp
      logical,      dimension(:), pointer     :: mask => NULL()
      real(kind=8), dimension(:), allocatable :: gap

      if (idebug.ge.4) call write_log('--- start subroutine locate_interpen_1d ---')
      associate( nx => sf_whl%nx, ny => sf_whl%ny, x  => sf_whl%x,  y  => sf_whl%y,                     &
                 zw => sf_whl%z,  zr => sf_rai%z )

      if (nx.ne.1 .or. ny.le.1) then
         write(bufout,'(2(a,i4))') ' ERROR: locate_interpen_1d: grid should be 1 x ny, got',nx,' x',ny
         call write_log(1, bufout)
         call abort_run()
      endif

      ! 1. compute the gap, <0 means interpenetration.
      !    set gap to 999 outside the rail surface, where the interpolated rail height is 999

      allocate(gap(ny))
      do iy = 1, ny
         if (zr(iy).lt.888d0) then
            gap(iy) = zr(iy) - zw(iy)
         else
            gap(iy) = 999d0
         endif
      enddo

      if (.true. .and. idebug.ge.2) then
         call write_log(' Dump 1d gap-function to dump_gap1d.m ...')
         open(unit=ltmp, file='dump_gap1d.m')
         write(ltmp,'(2(a,i6),a)') 'mx=',nx,'; my=',ny,';'
         write(ltmp,'(a)') '%  iy      x_tr          y_tr          zw            zr            alph'
         write(ltmp,'(a)') 'tmp =['
         do iy = 1, ny
            write(ltmp,'(i5,5f14.6)') iy, x(iy), y(iy), zw(iy), zr(iy), sf_incl%vy(iy)
         enddo
         write(ltmp,*) '];'
         write(ltmp,*) 'xtr_gap = reshape(tmp(:,2),mx,my);'
         write(ltmp,*) 'ytr_gap = reshape(tmp(:,3),mx,my);'
         write(ltmp,*) 'zw_gap  = reshape(tmp(:,4),mx,my);'
         write(ltmp,*) 'zr_gap  = reshape(tmp(:,5),mx,my);'
         write(ltmp,*) 'a_gap   = reshape(tmp(:,6),mx,my);'
         write(ltmp,*) 'clear tmp;'
         close(ltmp)
         if (idebug.ge.4) then
            call write_log(' idebug>=4: aborting')
            call abort_run()
         endif
      endif

      if (.false.) then
         ix = (nx+1)/2
         write(bufout,'(a,i5)') ' profiles [x,y,zw,zr] at slice ix=',ix
         call write_log(1, bufout)
         do iy = 1, ny
            write(bufout,'(i5,4f10.5)') iy, x(iy), y(iy), zw(iy), zr(iy)
            call write_log(1, bufout)
         enddo
      endif

      ! 2. find overall minimum gap value

      if (idebug.ge.4) call write_log('--- locate overall minimum of gap-function ---')
      iy_min   = idmin(ny, gap, 1)
      gap_min  = gap(iy_min)
      delt_min = sf_incl%vy(iy_min)
      zw_min   = zw(iy_min)

      if (idebug.ge.4) then
         write(bufout,*) 'overall minimum is',gap_min,' at iy_min=',iy_min,', delt=', delt_min*180d0/pi, &
                ', zw=',zw_min
         call write_log(1, bufout)
      endif

      ! 3. estimate curvatures a1, b1 at the minimum gap location
      !    Note: using curve zr instead of gap function!

      if (gap_min.ge.-1d-4) then
         if (idebug.ge.4) call write_log('--- estimate curvatures ---')
         call estimate_curvatures_v1(x, y, zr, nx, ny, iy_min, nom_radius, a1, b1, idebug)
      endif

      ! 4. find local minima of gap-function, excluding those without interpenetration

      if (idebug.ge.4) call write_log('--- determine local minima of gap function ---')
      call reallocate_arr(mask, ny)
      call find_gap_locmin(x, y, zr, zw, sgn, gap, nx, ny, nummin, locmin, mask, idebug, MAX_NUM_MIN)

      ! 5. process list of local minima "(ix,iy)", find potential contact areas
      !    note: entire loop is skipped in case there's no interpenetration

      numcp0 = numcps
      numrem = nummin
      do while(numrem.gt.0)

         ! 5.a increment the number of contact patches

         numcps = numcps + 1
         if (numcps.gt.MAX_NUM_CPS) then
            call write_log('ERROR: MAX_NUM_CPS exceeded.')
            call abort_run()
         endif
         cp => allcps(numcps)%cp

         ! 5.b find the local minimum with largest interpenetration

         gap_locmin = 999d0
         ixmin  = 0
         do ilm = 1, nummin
            if (locmin(1,ilm).gt.0) then
               iy = locmin(2,ilm)
               if (gap(iy).lt.gap_locmin) then
                  ixmin  = ilm
                  gap_locmin = gap(iy)
               endif
            endif
         enddo

         ! 5.c store the initial contact point on the rail in terms of track coordinates

         iy = locmin(2,ixmin)

         if (idebug.ge.2) then
            write(bufout,'(a,i3,a,i5)') ' storing cp for local minimum',ixmin, ', iy=',iy
            call write_log(1, bufout)
         endif

         ! 3. approximate local gap function, estimate the minimum location

         if (idebug.ge.4) call write_log('--- estimate minimum ---')
         call estimate_minimum(x, y, gap, zr, ny, iy, xmin, ymin, gmin, zrmin, idebug)

         cp%micp    = vec( xmin, ymin, zrmin )
         cp%gap_min = gmin
         cp%nsub    = 1

         ! 5.d find the extent [iy_sta:iy_end] of the interpenetration region,
         !     with gap >= 0 on the outer columns, gap < 0 inside.

         do iy = 1, ny
            mask(iy) = (gap(iy).lt.0d0)
         enddo

         iy = locmin(2,ixmin)
         nguard = 1
         call find_bounding_box_1d(mask, ny, iy, nguard, iy_sta, iy_end, idebug)

         if (idebug.ge.2) then
            write(bufout,127) 'storing cp',numcps,                                                      &
                      ': (x,y,z)_icp = (',cp%micp%x(),',', sgn*cp%micp%y(),',', cp%micp%z(),')',        &
                        'range iy= [', iy_sta,':', iy_end, ']'
            call write_log(2, bufout)
 127        format(1x,a,i2,a,3(f12.6,a)/, 15x,a,2(i5,a))
         endif

         ! 5.e compute the "interpenetration center of gravity" in x and y
         !      - weighted mean \int_{gap<0} x * gap; \int_{gap<0} y * gap.
         !      - using approximation gap(x,y) = gap(xlc,y) + curv * (x-xlc)^2
         !      - integral over x evaluated analytically

         if (ic%is_roller()) then
            curv_y = 0.5d0 / nom_radius + 0.5d0 / roller_radius     ! add zw for wheel and roller
         else
            curv_y = 0.5d0 / nom_radius                             ! add zw for wheel profile
         endif

         call compute_wgt_center_1d(ic, cp, ny, x, y, zr, gap, sf_incl%vy, iy_sta, iy_end, curv_y, sgn, &
                        idebug )

         if (ic%use_initial_cp()) then
            if (idebug.ge.3) call write_log(' using the initial contact location')
            cp%wgt_xgap = cp%micp%x()
            cp%wgt_ygap = cp%micp%y()
         endif

         ! 5.h compute contact length [x0,x1] of potential contact needed for 1d gap function

         call compute_contact_length_1d( ny, x, gap, iy_sta, iy_end, curv_y, m_vw, xmin, xmax, idebug)

         ! 5.i compute improved estimates for ysta, yend using interpolation

         iy   = locmin(2,ixmin)
         ii_s = iy_sta
         ii_n = iy_end

         dgap = gap(ii_s+1) - gap(ii_s)
         has_zero = (gap(ii_s)*gap(ii_s+1).le.-1d-12 .and. abs(dgap).ge.1d-9)
         if (.not.has_zero) then
            ysta = y(ii_s)
            zsta = zr(ii_s)
         else
            ysta = y(ii_s) - gap(ii_s) * (y(ii_s+1) - y(ii_s)) / dgap
            zsta = zr(ii_s) - gap(ii_s) * (zr(ii_s+1) - zr(ii_s)) / dgap
         endif
         dgap = gap(ii_n-1) - gap(ii_n)
         has_zero = (gap(ii_n)*gap(ii_n-1).le.-1d-12 .and. abs(dgap).ge.1d-9)
         if (.not.has_zero) then
            yend = y(ii_n)
            zend = zr(ii_n)
         else
            yend = y(ii_n) - gap(ii_n) * (y(ii_n-1) - y(ii_n)) / dgap
            zend = zr(ii_n) - gap(ii_n) * (zr(ii_n-1) - zr(ii_n)) / dgap
         endif

         if (idebug.ge.4) then
            write(bufout,'(a,i4,a,2f9.3,a,2f11.6,a,2f8.3)') ' ii_s=',ii_s,', y=',y(ii_s), y(ii_s+1),    &
                   ', gap=',gap(ii_s), gap(ii_s+1), ', y/zsta=',ysta,zsta
            call write_log(1, bufout)
            write(bufout,'(a,i4,a,2f9.3,a,2f11.6,a,2f8.3)') ' ii_n=',ii_n,', y=',y(ii_n), y(ii_n-1),     &
                   ', gap=',gap(ii_n), gap(ii_n-1), ', y/zend=',yend,zend
            call write_log(1, bufout)
         endif

         ! 5.i store extent of potential contact area in terms of track coordinates

         cp%xsta = xmin
         cp%xend = xmax
         cp%ysta = ysta
         cp%yend = yend
         cp%zsta = zsta
         cp%zend = zend

         ! 5.f set dummy extent for (u,v) coordinates for interpenetration area

         cp%usta = 1d0
         cp%uend = 0d0
         cp%vsta = 1d0
         cp%vend = 0d0

         ! 5.j extend pot.contact beyond first/last points of gap function (???)

         if (iy_sta.eq. 1) cp%ysta = min(cp%ysta, 2d0*cp%wgt_ygap - y(ii_n))
         if (iy_end.eq.ny) cp%yend = max(cp%yend,                   y(ii_s))

         if (idebug.ge.3) then
            write(bufout,'(2(a,i4,a,f12.6))') '   y(',iy_sta,')=', y(ii_s),  ',   y(',iy_end,')=',y(ii_n)
            call write_log(1,bufout)
            write(bufout,'(2(a,i4,a,f12.6))') ' gap(',iy_sta,')=', gap(ii_s),', gap(',iy_end,')=',gap(ii_n)
            call write_log(1,bufout)
            write(bufout,'(4(a,f12.6))') '    y_sta =', cp%ysta, ',    y_end =',cp%yend, '    z_sta =', &
                        cp%zsta, ',    z_end =',cp%zend
            call write_log(1,bufout)
         endif

         ! 5.k remove all local minima contained in this interpenetration region

         do ilm = 1, nummin
            iy = locmin(2,ilm)
            if (iy.ge.iy_sta .and. iy.le.iy_end) then
               if (ilm.eq.ixmin) then
                  if (idebug.ge.3) then
                     write(bufout,*) 'processed local minimum',ilm
                     call write_log(1, bufout)
                  endif
               else
                  if (idebug.ge.2) then
                     write(bufout,*) 'discarding local minimum',ilm,', contained within interpen.region'
                     call write_log(1, bufout)
                  endif
               endif
               locmin(1,ilm) = -1
               numrem = numrem - 1
            endif
         enddo

      enddo ! while(numrem>0)

      if (idebug.ge.2) then
         write(bufout,'(a,i3,a,i2,a)') ' contact search: gap has',nummin,' local minima, using',        &
                numcps-numcp0, ' contact problems'
         call write_log(1, bufout)
      endif

      deallocate(gap, mask)
      end associate
      if (idebug.ge.4) call write_log('--- end subroutine locate_interpen_1d ---')

   end subroutine locate_interpen_1d

!------------------------------------------------------------------------------------------------------------

   subroutine locate_interpen_2d( ic, sf_whl, sf_rai, uv_whl, sf_incl, sgn, dx, nom_radius, gap_min,    &
                                delt_min, zw_min, a1, b1, numcps, allcps, idebug)
!--purpose: find the points in the wheel mesh with (locally) minimum gap / maximum interpenetration,
!           define the corresponding interpenetration regions
      implicit none
!--subroutine arguments:
      type(t_ic)                  :: ic
      type(t_grid)                :: sf_whl, sf_rai      ! wheel & rail surfaces at gap-mesh positions
                                                         ! in view coordinates
      type(t_gridfnc3)            :: uv_whl              ! wheel (u,v) coordinates at gap-mesh positions
      type(t_gridfnc3)            :: sf_incl             ! rail surface inclination wrt view orientation
                                                         ! at y-pos of sf_rai
      real(kind=8), intent(in)    :: sgn                 ! +1/-1 for right/left w/r-combination
      real(kind=8), intent(in)    :: dx
      real(kind=8), intent(in)    :: nom_radius          ! nominal wheel radius
      real(kind=8), intent(out)   :: gap_min             ! overall minimum gap found (view direction)
      real(kind=8), intent(out)   :: delt_min            ! contact angle at the overall minimum gap value
      real(kind=8), intent(out)   :: zw_min              ! wheel surface height at overall minimum gap value
      real(kind=8), intent(out)   :: a1, b1              ! curvatures at the overall minimum gap value
      integer,      intent(inout) :: numcps              ! number of contact problems
      type(p_cpatch)              :: allcps(MAX_NUM_CPS) ! data of contact problems
      integer,      intent(in)    :: idebug
!--local variables:
      integer,      parameter :: MAX_NUM_MIN = 1000
      logical        :: has_uv
      integer        :: ix, iy, iy_sta, iy_end, ix_sta, ix_end, ii_min, nguard, icheck
      integer        :: ii, ii_n, ii_s, ii_e, ii_w, ilm, ixmin, numcp0, nummin, numrem
      integer        :: locmin(2,MAX_NUM_MIN)
      real(kind=8)   :: gap_locmin
      type(t_cpatch),             pointer     :: cp
      logical,      dimension(:), pointer     :: mask => NULL()
      real(kind=8), dimension(:), allocatable :: gap

      if (idebug.ge.4) call write_log('--- start subroutine locate_interpen_2d ---')
      associate( nx => sf_whl%nx, ny => sf_whl%ny, x  => sf_whl%x,  y  => sf_whl%y,                     &
                 zw => sf_whl%z,  uw => uv_whl%vx, vw => uv_whl%vy, zr => sf_rai%z )

      ! (u,v) coordinates are optionally computed for variable wheel profiles

      has_uv = uv_whl%is_defined()

      ! 1. compute the gap, <0 means interpenetration.
      !    set gap to 999 outside the rail surface, where the interpolated rail height is 999

      allocate(gap(nx*ny))
      do ii = 1, nx*ny
         if (zr(ii).lt.888d0) then
            gap(ii) = zr(ii) - zw(ii)
         else
            gap(ii) = 999d0
         endif
      enddo

      if (.true. .and. idebug.ge.2) then
         call write_log(' Dump 2d gap-function to dump_gap2d.m ...')
         open(unit=ltmp, file='dump_gap2d.m')
         write(ltmp,'(2(a,i6),a)') 'nx=',nx,'; ny=',ny,';'
         if (has_uv) then
            write(ltmp,'(2a)') '%  ix   iy      x_vw         y_vw         zw           zr',             &
                                                                        '           uw           vw'
         else
            write(ltmp,'(a)') '%  ix   iy      x_vw         y_vw         zw           zr'
         endif
         write(ltmp,'(a)') 'tmp =['
         do iy = 1, ny
            do ix = 1, nx
               ii = ix + (iy-1) * nx
               if (has_uv) then
                  write(ltmp,'(2i5,6f13.6)') ix, iy, x(ii), y(ii), zw(ii), zr(ii), uw(ii), vw(ii)
               else
                  write(ltmp,'(2i5,6f13.6)') ix, iy, x(ii), y(ii), zw(ii), zr(ii)
               endif
            enddo
         enddo
         write(ltmp,*) '];'
         write(ltmp,*) 'xvw_gap = reshape(tmp(:,3),nx,ny);'
         write(ltmp,*) 'yvw_gap = reshape(tmp(:,4),nx,ny);'
         write(ltmp,*) 'zw_gap  = reshape(tmp(:,5),nx,ny);'
         write(ltmp,*) 'zr_gap  = reshape(tmp(:,6),nx,ny);'
         if (has_uv) then
            write(ltmp,*) 'uw_gap  = reshape(tmp(:,7),nx,ny);'
            write(ltmp,*) 'vw_gap  = reshape(tmp(:,8),nx,ny);'
         endif
         write(ltmp,*) 'clear tmp;'
         close(ltmp)
      endif

      ! 2. find overall minimum gap value & surface inclination in view coordinates

      if (idebug.ge.4) call write_log('--- locate overall minimum of gap-function ---')
      ii_min   = idmin(nx*ny, gap, 1)
      gap_min  = gap(ii_min)
      delt_min = sf_incl%vy(ii_min)
      zw_min   = zw(ii_min)

      if (idebug.ge.4) then
         write(bufout,*) 'overall minimum is',gap_min,' at ii_min=',ii_min,', delt=', delt_min*180d0/pi, &
                        ', zw_vw=',zw_min
         call write_log(1, bufout)
      endif

      ! 3. estimate curvatures a1, b1 at the minimum gap location

      if (gap_min.ge.-1d-4) then
         if (idebug.ge.4) call write_log('--- estimate curvatures ---')
         call estimate_curvatures_v1(x, y, gap, nx, ny, ii_min, nom_radius, a1, b1, idebug)
      endif

      ! 4. find local minima of gap-function, excluding those without interpenetration

      if (idebug.ge.4) call write_log('--- determine local minima of gap function ---')
      call reallocate_arr(mask, nx*ny)
      call find_gap_locmin(x, y, zr, zw, sgn, gap, nx, ny, nummin, locmin, mask, idebug, MAX_NUM_MIN)

      ! 5. process list of local minima "(ix,iy)", find potential contact areas
      !    note: entire loop is skipped in case there's no interpenetration

      numcp0 = numcps
      numrem = nummin
      do while(numrem.gt.0)

         ! 5.a increment the number of contact patches

         numcps = numcps + 1
         if (numcps.gt.MAX_NUM_CPS) then
            call write_log('ERROR: MAX_NUM_CPS exceeded.')
            call abort_run()
         endif
         cp => allcps(numcps)%cp

         ! 5.b find the local minimum with largest interpenetration

         gap_locmin = 999d0
         ixmin  = 0
         do ilm = 1, nummin
            if (locmin(1,ilm).gt.0) then
               ix = locmin(1,ilm)
               iy = locmin(2,ilm)
               ii = ix + (iy-1) * nx
               if (gap(ii).lt.gap_locmin) then
                  ixmin  = ilm
                  gap_locmin = gap(ii)
               endif
            endif
         enddo

         ! 5.c store the initial contact point on the rail in terms of track coordinates

         ix = locmin(1,ixmin)
         iy = locmin(2,ixmin)
         ii = ix + (iy-1) * nx

         if (idebug.ge.2) then
            write(bufout,'(a,i3,2(a,i5),a)') ' storing cp for local minimum',ixmin, ', ix=(',iy,        &
                ',',ix,')'
            call write_log(1, bufout)
         endif

         cp%micp = vec( x(ii), y(ii), zr(ii) )
         cp%nsub = 1

         ! 5.d find the extent [ix_sta:ix_end] x [iy_sta:iy_end] of the interpenetration region,
         !     with gap >= 0 on the outer rows & columns, gap < 0 inside.

         do ii = 1, nx*ny
            mask(ii) = (gap(ii).lt.0d0)
         enddo

         ix = locmin(1,ixmin)
         iy = locmin(2,ixmin)
         nguard = 1
         icheck = 1
         call find_bounding_box_2d(mask, nx, ny, ix, iy, nguard, ix_sta, ix_end, iy_sta, iy_end,        &
                                      icheck, idebug)

         if (idebug.ge.2) then
            write(bufout,127) 'storing cp',numcps,                                                      &
                      ': (x,y,z)_icp = (',cp%micp%x(),',', sgn*cp%micp%y(),',', cp%micp%z(),')',        &
                        'range iy= [', iy_sta,':', iy_end, '], ix= [', ix_sta,':', ix_end,']'
            call write_log(2, bufout)
 127        format(1x,a,i2,a,3(f12.6,a)/, 15x,a,4(i5,a))
         endif

         ! 5.e compute the "interpenetration center of gravity" in x and y
         !      - weighted mean \int_{gap<0} x * gap; \int_{gap<0} y * gap.

         call compute_wgt_center_2d( ic, cp, nx, ny, dx, x, y, zr, gap, sf_incl%vy, ix_sta, ix_end,     &
                iy_sta, iy_end, sgn, idebug )

         ! 5.f determine extent of (u,v) coordinates for interpenetration area

         if (has_uv) then
            call get_uv_extent( cp, nx, ny, uw, vw, ix_sta, ix_end, iy_sta, iy_end, idebug)
         else
            cp%usta = 1d0
            cp%uend = 0d0
            cp%vsta = 1d0
            cp%vend = 0d0
         endif

         ! 5.i store extent of potential contact area in terms of view coordinates

         ix = locmin(1,ixmin)
         iy = locmin(2,ixmin)

         ii_w = ix_sta + (iy-1) * nx
         ii_e = ix_end + (iy-1) * nx
         ii_s = ix + (iy_sta-1) * nx
         ii_n = ix + (iy_end-1) * nx

         cp%xsta = x(ii_w)
         cp%xend = x(ii_e)
         cp%ysta = y(ii_s)
         cp%yend = y(ii_n)
         cp%zsta = zr(ii_s)
         cp%zend = zr(ii_n)

         if (idebug.ge.2) then
            write(bufout,'(3(a,f8.3),a)') ' x in [',x(ii_w),',',x(ii_e),']'
            call write_log(1,bufout)
         endif

         if (idebug.ge.3) then
            write(bufout,'(2(a,i4,a,f12.6))') '   y(',iy_sta,')=', y(ii_s),  ',   y(',iy_end,')=',y(ii_n)
            call write_log(1,bufout)
            write(bufout,'(2(a,i4,a,f12.6))') ' gap(',iy_sta,')=', gap(ii_s),', gap(',iy_end,')=',gap(ii_n)
            call write_log(1,bufout)
            write(bufout,'(4(a,f12.6))') '    y_sta =', cp%ysta, ',    y_end =',cp%yend,                &
                                         ',   z_sta =', cp%zsta, ',    z_end =',cp%zend
            call write_log(1,bufout)
            !!! write(bufout,'(2(a,f12.6))') '    x_sta =', cp%xsta, ',    x_end =',cp%xend
            !!! call write_log(1,bufout)
         endif

         if (idebug.ge.6) then
            ! print values on wheel profile contained in interpenetration region
            call write_log('           ywh        zwh')
            do iy = iy_sta, iy_end
               ii = 116 + (iy-1) * nx
               write(bufout,'(a,i3,a,2f11.6)') ' iy=',iy,':',y(iy), zw(iy)
               call write_log(1, bufout)
            enddo
         endif

         ! 5.j remove all local minima contained in this interpenetration region

         do ilm = 1, nummin
            ix = locmin(1,ilm)
            iy = locmin(2,ilm)
            if (iy.ge.iy_sta .and. iy.le.iy_end .and. ix.ge.ix_sta .and. ix.le.ix_end) then
               if (ilm.eq.ixmin) then
                  if (idebug.ge.3) write(*,*) 'processed local minimum',ilm
               else
                  if (idebug.ge.2) write(*,*) 'discarding local minimum',ilm,                           &
                                ', contained within interpen.region'
               endif
               locmin(1,ilm) = -1
               numrem = numrem - 1
            endif
         enddo

      enddo ! while(numrem>0)

      if (idebug.ge.2) then
         write(bufout,'(a,i3,a,i2,a)') ' contact search: gap has',nummin,' local minima, using',        &
                numcps-numcp0,' contact problems'
         call write_log(1, bufout)
      endif

      deallocate(gap, mask)
      end associate
      if (idebug.ge.4) call write_log('--- end subroutine locate_interpen_2d ---')

   end subroutine locate_interpen_2d

!------------------------------------------------------------------------------------------------------------

   subroutine compute_contact_length_1d( ny, x, gap, iy_sta, iy_end, curv_y, m_vw, xmin, xmax, idebug)
!--purpose: estimate extent [x0,x1] of potential contact needed for 1d gap function
      implicit none
!--subroutine arguments:
      integer,                    intent(in)  :: ny, idebug, iy_sta, iy_end
      real(kind=8),               intent(in)  :: curv_y, x(ny), gap(ny)
      type(t_marker),             intent(in)  :: m_vw                ! view-marker in track coordinates
      real(kind=8),               intent(out) :: xmin, xmax
!--local variables:
      integer           :: iy
      real(kind=8)      :: cos_a, pen_y, tzero

      ! get cos( delt_vw ) - correction factor between gaps / curvatures of vertical and view directions

      cos_a = m_vw%rot%r(5)     ! == cos( m_vw%roll() )
      ! write(bufout,'(3(a,f12.4))') ' vw_roll =', m_vw%roll(),', cos_a=',cos_a,' =', m_vw%rot%r(5)
      ! call write_log(1, bufout)

      ! estimate using quadratic fit at x_lc from contact locus

      xmin =  999d0
      xmax = -999d0
      do iy = iy_sta, iy_end
         if (gap(iy).le.0d0) then
            ! gap in view direction, curv_y in vertical direction, effective curvature cos_a*curv_y
            pen_y  = -gap(iy)
            tzero  = sqrt(pen_y / (cos_a * curv_y))
            xmin   = min(xmin, x(iy)-tzero)
            xmax   = max(xmax, x(iy)+tzero)
            if (idebug.ge.4) then
               write(bufout,'(a,i4,5(a,f9.4),a)') ' iy=',iy,': xlc=',x(iy), ', pen=', pen_y,            &
                        ', curv=', curv_y, ', x in [', x(iy)-tzero,',', x(iy)+tzero, ']'
               call write_log(1,bufout)
            endif
         endif
      enddo

      if (idebug.ge.2) then
         write(bufout,'(3(a,f8.3),a)') ' x in [',xmin,',',xmax,']'
         call write_log(1,bufout)
      endif

      if (xmin.ge.xmax) then
         call write_log(' INTERNAL ERROR length_1d')
         call abort_run()
      endif

   end subroutine compute_contact_length_1d

!------------------------------------------------------------------------------------------------------------

   subroutine compute_wgt_center_1d( ic, cp, ny, x, y, z, gap, sf_incl, iy_sta, iy_end, curv_y, sgn,    &
                        idebug )
!--purpose: compute the "interpenetration center of gravity" in x and y for 1d gap function
      implicit none
!--subroutine arguments:
      type(t_ic)                              :: ic
      type(t_cpatch)                          :: cp
      integer,                    intent(in)  :: ny, iy_sta, iy_end, idebug
      real(kind=8),               intent(in)  :: curv_y, sgn
      real(kind=8), dimension(:), intent(in)  :: x(ny), y(ny), z(ny), gap(ny), sf_incl(ny)
!--local variables:
      logical      :: use_linear_wgt
      integer      :: iy
      real(kind=8) :: pen_y, tzero, wgt_xgap, wgt_ygap, wgt_zgap, wgt_agap, totgap, part_xgap,          &
                      part_ygap, part_zgap, part_agap, part_tot

      use_linear_wgt = (ic%gapwgt.eq.1)

      ! 5.e compute the "interpenetration center of gravity" in x and y, weighted average angle
      !      - weighted mean \int_{gap<0} x * gap; \int_{gap<0} y * gap; \int_{gap<0} angl * gap.
      !      - using approximation gap(x,y) = gap(xlc,y) + curv * (x-xlc)^2
      !      - integral over x evaluated analytically

      wgt_xgap = 0d0
      wgt_ygap = 0d0
      wgt_zgap = 0d0
      wgt_agap = 0d0
      totgap   = 0d0

      do iy = iy_sta, iy_end
         if (gap(iy).le.0d0) then
            pen_y  = -gap(iy)
            tzero  = sqrt(pen_y / curv_y)

            ! a == 1/2r(y) == curv_y, c == -g(x_lc,y) == pen
            if (use_linear_wgt) then
               ! linear gap:                2c t_0           -     2a/3 t0^3
               part_tot  = max(1d-10, 2d0 * pen_y * tzero - curv_y * tzero**3 * 2d0 / 3d0)
            else
               ! squared.gap:              2c^2 t_0          -     4ac/3 t0^3      +    2a^2/5 t0^5
               part_tot  = max(1d-10, 2d0 * pen_y**2 * tzero - pen_y * curv_y * tzero**3 * 4d0 / 3d0       &
                                                                  + curv_y**2 * tzero**5 * 2d0 / 5d0 )
            endif
            part_xgap = x(iy) * part_tot
            part_ygap = y(iy) * part_tot
            part_zgap = z(iy) * part_tot
            part_agap = sf_incl(iy) * part_tot

            if (idebug.ge.5 .and. abs(y(iy)-726.4d0).lt.0.05d0) then
               write(bufout,'(6g12.4)') x(iy), y(iy), gap(iy), pen_y, curv_y, tzero
               call write_log(1, bufout)
            endif
            totgap   = totgap   + part_tot
            wgt_xgap = wgt_xgap + part_xgap
            wgt_ygap = wgt_ygap + part_ygap
            wgt_zgap = wgt_zgap + part_zgap
            wgt_agap = wgt_agap + part_agap
         else
            part_tot  = 0d0
            part_xgap = 0d0
            part_ygap = 0d0
            part_zgap = 0d0
            part_agap = 0d0
         endif

         if (idebug.ge.4) then
            write(bufout, 197) iy, y(iy), part_xgap, part_ygap, part_zgap, part_agap, part_tot, gap(iy)
            call write_log(1, bufout)
 197        format(' row',i4,', y=',f6.2,': wgtx=',g12.4,', wgty=',g12.4,', wgtz=',g12.4,', wgta=',g12.4, &
                        ', tot=',g12.4,', gap=',g12.4)
         endif
      enddo

      if (idebug.ge.4) then
         write(bufout, '(14x, 5(a,g12.4))') 'total wgtx=',wgt_xgap, ', wgty=',wgt_ygap, ', wgtz=',      &
                        wgt_zgap, ', wgta=',wgt_agap, ', tot=', totgap
         call write_log(1, bufout)
      endif

      wgt_xgap = wgt_xgap / totgap
      wgt_ygap = wgt_ygap / totgap
      wgt_zgap = wgt_zgap / totgap
      wgt_agap = wgt_agap / totgap

      if (idebug.ge.2) then
         write(bufout,'(4(a,f12.6))') ' interpen. center of gravity: x_cg=',wgt_xgap,                &
             ', y_cg=', sgn*wgt_ygap, ', z_cg=',wgt_zgap, ', a_avg=', wgt_agap
         call write_log(1, bufout)
      endif

      cp%wgt_xgap = wgt_xgap
      cp%wgt_ygap = wgt_ygap
      cp%wgt_zgap = wgt_zgap
      cp%wgt_agap = wgt_agap
      cp%totgap   = totgap

   end subroutine compute_wgt_center_1d

!------------------------------------------------------------------------------------------------------------

   subroutine compute_wgt_center_2d( ic, cp, nx, ny, dx, x, y, z, gap, sf_incl, ix_sta, ix_end, iy_sta, &
                iy_end, sgn, idebug )
!--purpose: compute the "interpenetration center of gravity" in x and y for 2d gap function
      implicit none
!--subroutine arguments:
      type(t_ic)                              :: ic
      type(t_cpatch)                          :: cp
      integer,                    intent(in)  :: nx, ny, ix_sta, ix_end, iy_sta, iy_end, idebug
      real(kind=8),               intent(in)  :: dx, sgn
      real(kind=8), dimension(:), intent(in)  :: x(nx*ny), y(nx*ny), z(nx*ny), gap(nx*ny), sf_incl(nx*ny)
!--local variables:
      logical      :: use_linear_wgt
      integer      :: ii, ix, iy
      real(kind=8) :: wgt_xgap, wgt_ygap, wgt_zgap, wgt_agap, totgap, part_xgap, part_ygap, part_zgap,  &
                      part_agap, part_tot, cell_tot, part_max

      use_linear_wgt = (ic%gapwgt.eq.1)

      ! 5.e compute the "interpenetration center of gravity" in x and y, weighted avg. angle
      !      - weighted mean \int_{gap<0} x * -gap; \int_{gap<0} y * -gap; \int_{gap<0} angl * -gap.

      wgt_xgap = 0d0
      wgt_ygap = 0d0
      wgt_zgap = 0d0
      wgt_agap = 0d0
      totgap   = 0d0

      do iy = iy_sta, iy_end
         part_xgap = 0d0
         part_ygap = 0d0
         part_zgap = 0d0
         part_agap = 0d0
         part_tot  = 0d0
         part_max  = -999d0
         do ix = ix_sta, ix_end
            ii = ix + (iy-1) * nx
            if (-gap(ii).gt.part_max) part_max = -gap(ii)
            if (-gap(ii).gt.0d0) then
               if (use_linear_wgt) then
                  ! linear weighting:             -gap>0
                  cell_tot  =      dx * max(-gap(ii), 1d-10)
               else
                  ! squared weighting:            gap^2>0
                  cell_tot  =      dx * max(gap(ii)**2, 1d-10)
               endif
               part_xgap = part_xgap +       x(ii) * cell_tot
               part_ygap = part_ygap +       y(ii) * cell_tot
               part_zgap = part_zgap +       z(ii) * cell_tot
               part_agap = part_agap + sf_incl(ii) * cell_tot
               part_tot  = part_tot  +               cell_tot
               if (idebug.ge.5 .and. abs(y(ii)-13.87d0).lt.0.003d0) then
                  write(bufout,'(i6,4g14.6)') ii, x(ii), y(ii), -gap(ii), sf_incl(ii)
                  call write_log(1, bufout)
               endif
            endif
         enddo
         wgt_xgap = wgt_xgap + part_xgap
         wgt_ygap = wgt_ygap + part_ygap
         wgt_zgap = wgt_zgap + part_zgap
         wgt_agap = wgt_agap + part_agap
         totgap   = totgap   + part_tot

         if (idebug.ge.4) then
            write(bufout, 197) iy, y(ii), part_xgap, part_ygap, part_zgap, part_agap, part_tot, part_max
            call write_log(1, bufout)
 197        format(' row',i4,', y=',f7.3,': wgtx=',g12.4,', wgty=',g12.4,', wgtz=',g12.4,               &
                ', wgta=',g12.4, ', tot=',g12.4,', gap=',g12.4)
         endif
      enddo

      if (idebug.ge.4) then
         write(bufout, '(14x, 5(a,g12.4))') 'total wgtx=',wgt_xgap,', wgty=', wgt_ygap,                 &
                ', wgtz=',wgt_zgap, ', wgta=',wgt_agap, ', tot=', totgap
         call write_log(1, bufout)
      endif

      wgt_xgap = wgt_xgap / max(totgap, 1d-10)
      wgt_ygap = wgt_ygap / max(totgap, 1d-10)
      wgt_zgap = wgt_zgap / max(totgap, 1d-10)
      wgt_agap = wgt_agap / max(totgap, 1d-10)

      if (idebug.ge.2) then
         write(bufout,'(4(a,f12.6))') ' interpen. center of gravity: x_cg=',wgt_xgap,                &
             ', y_cg=', sgn*wgt_ygap, ', z_cg=',wgt_zgap, ', a_avg=', wgt_agap
         call write_log(1, bufout)
      endif

      cp%wgt_xgap = wgt_xgap
      cp%wgt_ygap = wgt_ygap
      cp%wgt_zgap = wgt_zgap
      cp%wgt_agap = wgt_agap
      cp%totgap   = totgap

   end subroutine compute_wgt_center_2d

!------------------------------------------------------------------------------------------------------------

   subroutine get_uv_extent( cp, nx, ny, uw, vw, ix_sta, ix_end, iy_sta, iy_end, idebug )
!--purpose: determine the range of surface parameters [usta,uend] x [vsta,vend] covering the whole
!           interpenetration area
      implicit none
!--subroutine arguments:
      type(t_cpatch)                          :: cp
      integer,                    intent(in)  :: nx, ny, ix_sta, ix_end, iy_sta, iy_end, idebug
      real(kind=8), dimension(:), intent(in)  :: uw(nx,ny), vw(nx,ny)
!--local variables:
      integer      :: ix, iy
      real(kind=8) :: umin, umax, vmin, vmax

      ! (u,v) can be oriented opposite to (ix,iy) (esp. v for wheel surface)

      umin =  1d10
      umax = -1d10
      vmin =  1d10
      vmax = -1d10

      ! boundary ix_sta

      ix = ix_sta
      do iy = iy_sta, iy_end
         umin = min(umin, uw(ix,iy))
         umax = max(umax, uw(ix,iy))
         vmin = min(vmin, vw(ix,iy))
         vmax = max(vmax, vw(ix,iy))
      enddo

      ! boundary ix_end

      ix = ix_end
      do iy = iy_sta, iy_end
         umin = min(umin, uw(ix,iy))
         umax = max(umax, uw(ix,iy))
         vmin = min(vmin, vw(ix,iy))
         vmax = max(vmax, vw(ix,iy))
      enddo

      ! boundary iy_sta

      iy = iy_sta
      do ix = ix_sta, ix_end
         umin = min(umin, uw(ix,iy))
         umax = max(umax, uw(ix,iy))
         vmin = min(vmin, vw(ix,iy))
         vmax = max(vmax, vw(ix,iy))
      enddo

      ! boundary iy_end

      do ix = ix_sta, ix_end
         umin = min(umin, uw(ix,iy))
         umax = max(umax, uw(ix,iy))
         vmin = min(vmin, vw(ix,iy))
         vmax = max(vmax, vw(ix,iy))
      enddo

      cp%usta = umin
      cp%uend = umax
      cp%vsta = vmin
      cp%vend = vmax

      if (idebug.ge.2) then
         write(bufout,'(4(a,f10.4),a)') ' surface parameters u = [',umin,',',umax,'], v = [',vmin,',',  &
                        vmax,']'
         call write_log(1, bufout)
      endif

   end subroutine get_uv_extent

!------------------------------------------------------------------------------------------------------------

   subroutine sort_cpatches_yvalue( numcps, allcps, idebug )
!--purpose: sort the contact problems with |y|-position descending, i.e. from field side to track center
      implicit none
!--subroutine arguments:
      type(p_cpatch)            :: allcps(MAX_NUM_CPS)  ! data of contact problems
      integer,      intent(in)  :: numcps, idebug
!--local variables:
      logical        :: lchanged
      integer        :: icp
      type(p_cpatch) :: tmp_cp

      lchanged = .true.
      do while(lchanged)

         lchanged = .false.

         do icp = 1, numcps-1

            if (idebug.ge.3) then
               write(bufout,'(2(a,i3,a,f10.6))') ' cp',icp,' has |yr|=', allcps(icp)%cp%wgt_ygap,',',   &
                               icp+1,' has |yr|=', allcps(icp+1)%cp%wgt_ygap
               call write_log(1, bufout)
            endif

            if (allcps(icp)%cp%wgt_ygap.lt.allcps(icp+1)%cp%wgt_ygap) then

               if (idebug.ge.3) then
                  write(bufout,'(2(a,i3))') '   ...swapping',icp,' and',icp+1
                  call write_log(1, bufout)
               endif

               tmp_cp        = allcps(icp)
               allcps(icp)   = allcps(icp+1)
               allcps(icp+1) = tmp_cp

               if (idebug.ge.3) then
                  write(bufout,'(2(a,i3,a,f10.6))') '   ...cp',icp,' has |yr|=',allcps(icp)%cp%wgt_ygap,   &
                        ',', icp+1,' has |yr|=',allcps(icp+1)%cp%wgt_ygap
                  call write_log(1, bufout)
               endif

               lchanged = .true.
            endif
         enddo

      enddo ! while(lchanged)

   end subroutine sort_cpatches_yvalue

!------------------------------------------------------------------------------------------------------------

   subroutine combine_cpatches( ic, numcps, allcps, angl_sep, dist_sep, dist_comb, idebug )
!--purpose: combine contact problems with overlapping sr-positions or lying too close together.
!           patches are sorted with sr-values decreasing (from field side to track center)
      implicit none
!--subroutine arguments:
      type(t_ic)                  :: ic
      type(p_cpatch)              :: allcps(MAX_NUM_CPS)  ! data of contact problems
      integer,      intent(in)    :: idebug
      integer,      intent(inout) :: numcps
      real(kind=8), intent(in)    :: angl_sep, dist_sep, dist_comb
!--local variables:
      logical, parameter         :: debug_comb = .false.
      integer                    :: icp, jcp
      real(kind=8)               :: dy, dz, dist, angl
      type(p_cpatch)             :: tmp_cp
      real(kind=8),  allocatable :: ysta(:), yend(:), zsta(:), zend(:)

      if (idebug.ge.2 .or. debug_comb) then
         write(bufout,'(a,f6.2,a,f5.1,a)') ' Separating patches with more than', angl_sep,' rad or',    &
                dist_sep,' mm separation'
         call write_log(1, bufout)
         write(bufout,'(a,13x,f5.1,a)') ' Combining patches with less than ', dist_comb,' mm separation'
         call write_log(1, bufout)
      endif

      allocate(ysta(numcps), yend(numcps), zsta(numcps), zend(numcps))
      
      icp = 1
      do while (icp.lt.numcps)

         ! get [ysta, yend] and [zsta, zend] for all contact patches

         do jcp = 1, numcps
            ysta(jcp) = allcps(jcp)%cp%ysta
            yend(jcp) = allcps(jcp)%cp%yend
            zsta(jcp) = allcps(jcp)%cp%zsta
            zend(jcp) = allcps(jcp)%cp%zend
         enddo

         ! print overview of cpatches y- and z-values (decreasing)

         if (idebug.ge.3 .or. debug_comb) then
            write(bufout,'(2(a,i3))') ' checking overlap of patches',icp,' and',icp+1
            call write_log(1, bufout)

            do jcp = 1, numcps
               write(bufout,'(a,i2,4(a,f9.4),a)') ' cp',jcp,': y=[',ysta(jcp),',', yend(jcp),           &
                        '], zr=[',zsta(jcp),',', zend(jcp),']'
               call write_log(1,bufout)
            enddo
         endif

         ! patches sorted with y-values decreasing: [ysta(i+1), yend(i+1)], dist, [ysta(i), yend(i)]

         dy    = ysta(icp) - yend(icp+1)
         dz    = zsta(icp) - zend(icp+1)
         dist  = sqrt(dy**2 + dz**2)
         angl  = allcps(icp)%cp%wgt_agap - allcps(icp+1)%cp%wgt_agap

         if (idebug.ge.4 .or. debug_comb) then
            write(bufout,'(2(a,i3),3(a,f8.3))') '   ...patches',icp,' and',icp+1,': a1=',               &
                   allcps(icp)%cp%wgt_agap, ', a2=', allcps(icp+1)%cp%wgt_agap,', diff=',angl
            call write_log(1, bufout)
            write(bufout,'(2(a,i3),3(a,f8.3))') '   ...patches',icp,' and',icp+1,': dist=', dist
            call write_log(1, bufout)
         endif

         if (abs(angl).gt.angl_sep .or. dist.gt.dist_sep) then

            ! keep separate when |a1-a2| >= angl_sep or dist >= dist_sep : move to next patch

            icp = icp + 1

         else

            ! combine patches i and i+1

            associate( cp_lft => allcps(icp+1)%cp, cp_rgt => allcps(icp)%cp )

            ! store y-value in between neighbouring patches, weighting factor

            if (dist.gt.dist_comb) then
               cp_rgt%nsub = cp_rgt%nsub + 1

               call enlarge_1d_real(cp_rgt%y_sep, cp_rgt%nsub, 5)
               call enlarge_1d_real(cp_rgt%f_sep,  cp_rgt%nsub, 5)

               cp_rgt%y_sep(cp_rgt%nsub-1) = (yend(icp+1) + ysta(icp)) / 2d0
               cp_rgt%f_sep(cp_rgt%nsub-1)  = ( (dist_sep - dist) / (dist_sep - dist_comb) )**0.8d0

            endif

            if ((idebug.ge.2 .or. debug_comb) .and. dist.gt.dist_comb) then
               write(bufout,'(2(a,i3),3(a,f8.3))') '   ...combining patches',icp,' and',icp+1,', dist=', &
                     dist,', y_sep=', cp_rgt%y_sep(cp_rgt%nsub-1),', fsep=', cp_rgt%f_sep(cp_rgt%nsub-1)
               call write_log(1, bufout)
            elseif (idebug.ge.2 .or. debug_comb) then
               write(bufout,'(2(a,i3),a,f8.3)') '   ...combining patches',icp,' and',icp+1,', dist=',dist
               call write_log(1, bufout)
            endif

            ! keep micp with largest interpenetration

            if (cp_lft%gap_min.lt.cp_rgt%gap_min) then
               cp_rgt%micp    = cp_lft%micp
               cp_rgt%gap_min = cp_lft%gap_min
            endif

            ! take union of x- and y-ranges

            cp_rgt%xsta = min(cp_rgt%xsta, cp_lft%xsta)
            cp_rgt%xend = max(cp_rgt%xend, cp_lft%xend)
            if (cp_lft%ysta.lt.cp_rgt%ysta) then
               cp_rgt%ysta = cp_lft%ysta
               cp_rgt%zsta = cp_lft%zsta
            endif
            if (cp_lft%yend.gt.cp_rgt%yend) then
               cp_rgt%yend = cp_lft%yend
               cp_rgt%zend = cp_lft%zend
            endif

            ! take union of u- and v-ranges, if not empty ([1,0])

            if (cp_lft%usta.lt.cp_lft%uend .and. cp_rgt%usta.lt.cp_rgt%uend) then
               cp_rgt%usta = min(cp_lft%usta, cp_rgt%usta)
               cp_rgt%uend = max(cp_lft%uend, cp_rgt%uend)
               cp_rgt%vsta = min(cp_lft%vsta, cp_rgt%vsta)
               cp_rgt%vend = max(cp_lft%vend, cp_rgt%vend)
            endif

            ! compute new weighted center

            if (ic%use_initial_cp()) then
               cp_rgt%wgt_xgap = cp_rgt%micp%x()
               cp_rgt%wgt_ygap = cp_rgt%micp%y()
            else
               cp_rgt%wgt_xgap = (cp_rgt%wgt_xgap * cp_rgt%totgap +                                     &
                                  cp_lft%wgt_xgap * cp_lft%totgap) / (cp_rgt%totgap + cp_lft%totgap)
               cp_rgt%wgt_ygap = (cp_rgt%wgt_ygap * cp_rgt%totgap +                                     &
                                  cp_lft%wgt_ygap * cp_lft%totgap) / (cp_rgt%totgap + cp_lft%totgap)
            endif
            cp_rgt%wgt_agap = (cp_rgt%wgt_agap * cp_rgt%totgap +                                     &
                                  cp_lft%wgt_agap * cp_lft%totgap) / (cp_rgt%totgap + cp_lft%totgap)
            cp_rgt%totgap   = cp_rgt%totgap + cp_lft%totgap
            
            end associate

            ! cycle remaining contact patches (pointers), putting icp+1 at end of list

            tmp_cp = allcps(icp+1)
            do jcp = icp+2, numcps
               allcps(jcp-1) = allcps(jcp)
            enddo
            allcps(numcps) = tmp_cp

            ! stay at icp, reduce number of patches

            numcps = numcps - 1

         endif
      enddo ! while icp<numcps

      ! reorder ysep, fsep in ascending order

      do icp = 1, numcps
         associate( cp => allcps(icp)%cp )
         if (cp%nsub.ge.3) then
            call reverse_order(cp%nsub-1, cp%y_sep)
            call reverse_order(cp%nsub-1, cp%f_sep)
         endif
         end associate
      enddo

      deallocate(ysta, yend, zsta, zend)

   end subroutine combine_cpatches

!------------------------------------------------------------------------------------------------------------

   subroutine turn_cpatch_refangle( numcps, allcps, angl_sep, dist_sep, dist_turn, idebug )
!--purpose: turn contact reference angles for patches that lie close together.
!           patches are sorted with sr-values decreasing (from field side to track center)
      implicit none
!--subroutine arguments:
      type(p_cpatch)              :: allcps(MAX_NUM_CPS)  ! data of contact problems
      integer,      intent(in)    :: idebug
      integer,      intent(inout) :: numcps
      real(kind=8), intent(in)    :: angl_sep, dist_sep, dist_turn
!--local variables:
      logical,       parameter    :: debug_turn = .false.
      real(kind=8),  parameter    :: tiny_sep = 1d-12
      logical                     :: use_wgt_angle
      integer                     :: icp, jcp
      real(kind=8)                :: dy, dz, dist, dangl, angl_comb, fac
      real(kind=8),  allocatable  :: ysta(:), yend(:), zsta(:), zend(:), angl(:), totgap(:)

      use_wgt_angle  = .true.

      ! turning not yet implemented for D=6, 8, 9 with sampled contact angle -- return

      if (.not.use_wgt_angle) return

      ! no turning in case dist_turn <= dist_sep

      if (dist_turn.le.dist_sep+tiny_sep) return

      if (idebug.ge.2 .or. debug_turn) then
         write(bufout,'(2(a,f5.1),a,f6.3,a)') ' Turning ref.angle for patches with ',dist_sep,          &
                ' <= dist <= ', dist_turn,' mm separation, |Dalpha| <=', angl_sep,' rad'
         call write_log(1, bufout)
      endif

      ! get [ysta, yend] and corresponding [zsta, zend] for all contact patches
      ! copy totgap and a-priori values of wgt_agap

      allocate(ysta(numcps), yend(numcps), zsta(numcps), zend(numcps), angl(numcps), totgap(numcps))

      do jcp = 1, numcps
         ysta(jcp)   = allcps(jcp)%cp%ysta
         yend(jcp)   = allcps(jcp)%cp%yend
         zsta(jcp)   = allcps(jcp)%cp%zsta
         zend(jcp)   = allcps(jcp)%cp%zend
         angl(jcp)   = allcps(jcp)%cp%wgt_agap
         totgap(jcp) = allcps(jcp)%cp%totgap
      enddo
      
      ! print overview of cpatches y- and sr-values (decreasing)

      if (idebug.ge.3 .or. debug_turn) then
         do jcp = 1, numcps
            write(bufout,'(a,i2,6(a,f9.4))') ' cp',jcp,': y=[',ysta(jcp),',', yend(jcp),                &
                     '], zr=[',zsta(jcp),',', zend(jcp),'], angl=',angl(jcp),', totgap=',totgap(jcp)
            call write_log(1,bufout)
         enddo
      endif

      ! for each pair (icp, icp+1): check and turn contact angles if needed

      do icp = 1, numcps-1

         ! patches sorted with y-values decreasing: [ysta(i+1), yend(i+1)], dist, [ysta(i), yend(i)]

         dy    = ysta(icp) - yend(icp+1)
         dz    = zsta(icp) - zend(icp+1)
         dist  = sqrt(dy**2 + dz**2)
         dangl = angl(icp) - angl(icp+1)

         if (abs(dangl).gt.angl_sep .or. dist.lt.dist_sep .or. dist.gt.dist_turn) then

            ! keep separate when |a1-a2| >= angl_sep or dist >= dist_turn : move to next patch

         else

            ! compute averaged angle at full combination

            angl_comb = (angl(icp) * totgap(icp) + angl(icp+1) * totgap(icp+1)) /                       &
                                                                        (totgap(icp) + totgap(icp+1))

            ! compute interpolation factor

            fac = (dist_turn - dist) / (dist_turn - dist_sep)

            ! add fac * difference --> accomodating turning from two adjacent contact patches

            allcps(icp  )%cp%wgt_agap = allcps(icp  )%cp%wgt_agap + fac * (angl_comb - angl(icp))
            allcps(icp+1)%cp%wgt_agap = allcps(icp+1)%cp%wgt_agap + fac * (angl_comb - angl(icp+1))

            if (idebug.ge.2 .or. debug_turn) then
               write(bufout,'(2(a,i2),3(a,f9.4),a)') ' cp',icp,' and',icp+1, ': dist=',dist,            &
                         ', turning with ',1d0-fac,' * a1,a2 + ',fac,' * comb'
               call write_log(1, bufout) 
               write(bufout,'(5(a,f9.4))') ' input a1=',angl(icp),', a2=',angl(icp+1),', comb=',        &
                        angl_comb,' blended=', allcps(icp)%cp%wgt_agap,',', allcps(icp+1)%cp%wgt_agap
               call write_log(1, bufout) 
            endif

         endif ! turn angles
      enddo ! for icp

      deallocate(ysta, yend, zsta, zend, angl, totgap)

   end subroutine turn_cpatch_refangle

!------------------------------------------------------------------------------------------------------------

   subroutine set_cpatch_reference( ic, icp, cp, prr, mtrk, sgn, nom_radius, idebug)
!--purpose: determine the contact angle and reference marker for one contact patch
      implicit none
!--subroutine arguments:
      type(t_ic),   intent(in)  :: ic
      integer,      intent(in)  :: icp                  ! current contact patch number
      type(t_cpatch)            :: cp
      type(t_grid)              :: prr                  ! rail profile in track coordinates
      type(t_marker)            :: mtrk                 ! roller profile marker
      real(kind=8), intent(in)  :: sgn                  ! +1/-1 for right/left w/r-combination
      real(kind=8), intent(in)  :: nom_radius           ! radius in case of a roller
      integer,      intent(in)  :: idebug
!--local variables:
      logical                   :: use_wgt_angle
      integer                   :: sub_ierror
      real(kind=8)              :: cref_x, cref_y, cref_z, zc_rol, r_y

      use_wgt_angle  = .true.

      if (ic%use_initial_cp()) then
         cref_x = cp%micp%x()
         cref_y = cp%micp%y()
      else
         cref_x = cp%wgt_xgap
         cref_y = cp%wgt_ygap
      endif

      ! 5.f interpolate rail profile to get sr and zr at reference y-position

      call spline_get_s_at_y( prr%spl, cref_y, cp%sr_ref, sub_ierror )
      call spline_get_xz_at_y( prr%spl, cref_y, sub_ierror, zout=cref_z )

      ! write(bufout,'(3(a,f12.4))') ' cref_s=',cp%sr_ref,', cref_y=',cref_y,', cref_z=',cref_z
      ! call write_log(1, bufout)

      if (ic%is_roller()) then
         ! cref_z as computed above lies in the principal profile, x==0, in track coordinates
         zc_rol = mtrk%z() + nom_radius
         ! correct cref_z for offset in x
         r_y    = zc_rol - cref_z
         cref_z = zc_rol - sqrt( r_y**2 - cref_x**2 )

         if (idebug.ge.2) then
            write(bufout,'(4(a,f12.4))') ' zc_rol=',zc_rol,', r_y=',r_y,', cref_x=',cref_x,             &
                        ', cref_z=',cref_z
            call write_log(1, bufout)
         endif
      endif

      ! 5.g compute the contact angle from track vertical to contact normal

      if (use_wgt_angle) then
         cp%delttr = cp%wgt_agap
      else
         call compute_contact_angle( prr, cref_y, cref_z, cp%delttr, idebug )
      endif

      if (idebug.ge.2) then
         write(bufout,129) ' estimated contact angle delttr=',sgn*cp%delttr,' rad =',                   &
             sgn*cp%delttr*180d0/pi,' deg'
         call write_log(1, bufout)
 129     format(a,f10.6,a,f12.6,a)
      endif

      ! 5.h create a marker for contact local coordinates w.r.t. track coordinates

      call marker_init(cp%mref)
      call marker_roll(cp%mref, cp%delttr, 0d0, 0d0)
      call marker_shift(cp%mref, cref_x, cref_y, cref_z )

      if (idebug.ge.2) then
         write(bufout,127) 'storing cp', icp,                                                        &
                   ': (x,y,z)_ref = (',cp%mref%x(),',', sgn*cp%mref%y(),',', cp%mref%z(),')',        &
                     '(x,y,z)_icp = (',cp%micp%x(),',', sgn*cp%micp%y(),',', cp%micp%z(),')'
         call write_log(2, bufout)
 127     format(1x,a,i2,a,3(f12.6,a)/, 15x,a,3(f12.6,a))
      endif

   end subroutine set_cpatch_reference

!------------------------------------------------------------------------------------------------------------

   subroutine compute_contact_angle( prr, yref, zref, delttr, idebug )
!--purpose: estimate the contact angle from track vertical to contact normal, i.e. the rail profile
!           inclination at position (yref, zref)
      implicit none
!--subroutine arguments:
      type(t_grid), intent(in)  :: prr
      real(kind=8), intent(in)  :: yref, zref
      integer,      intent(in)  :: idebug
      real(kind=8), intent(out) :: delttr
!--local variables:
      integer      :: sub_ierror
      real(kind=8) :: sr_ref(1), alph(1)

      ! compute s-position for position yref

      call spline_get_s_at_y( prr%spl, yref, sr_ref(1), sub_ierror )

      ! evaluate profile inclination atan2(dz, dy) at s-position sr_ref

      call spline_get_alpha_at_s( prr%spl, 1, sr_ref, alph, sub_ierror )

      delttr = alph(1)

      if (idebug.ge.3) then
         write(bufout,'(3(a,f12.6))') ' yref=',yref,', zref=',zref,': sr_ref=',sr_ref(1)
         call write_log(1, bufout)
         write(bufout,'(2(a,f14.8))') ' inclination on rail alph =',alph(1)
         call write_log(1, bufout)
      endif

   end subroutine compute_contact_angle

!------------------------------------------------------------------------------------------------------------

   subroutine set_planar_potcon( ic, icp, cp, prr, mtrk, sgn, nom_radius, dx_in, ds_in, npot_max, idebug )
!--purpose: define the potential contact area for planar computations
      implicit none
!--subroutine arguments:
      type(t_ic)                :: ic
      type(t_cpatch)            :: cp
      type(t_grid)              :: prr
      integer,      intent(in)  :: icp, npot_max, idebug
      type(t_marker)            :: mtrk                 ! roller profile marker
      real(kind=8), intent(in)  :: nom_radius           ! radius in case of a roller
      real(kind=8), intent(in)  :: sgn, dx_in, ds_in    ! dx, ds: user input
!--local variables:
      integer        :: sub_ierror, mx, my, ix_l, ix_h, iy_l, iy_h
      real(kind=8)   :: cref_x, zc_rol, rsta, rend, zsta, zend, xp_l, xp_h, sp_l, sp_h, fac, du, dv, tmp
      type(t_marker) :: msta, mend

      ! compute positions [ssta,send] and [zsta_tr,zend_tr] on rail profile for [ysta_tr, yend_tr],

      call spline_get_xz_at_y( prr%spl, cp%ysta, sub_ierror, zout=zsta )
      call spline_get_xz_at_y( prr%spl, cp%yend, sub_ierror, zout=zend )

      if (ic%is_roller()) then
         ! msta/end as computed above lie in the principal profile, x==0, in track coordinates
         zc_rol = mtrk%z() + nom_radius
         cref_x = cp%mref%x()

         ! correct msta%z, mend%z for offset in x
         rsta  = zc_rol - zsta
         zsta  = zc_rol - sqrt( rsta**2 - cref_x**2 )
         rend  = zc_rol - zend
         zend  = zc_rol - sqrt( rend**2 - cref_x**2 )

         if (idebug.ge.2) then
            write(bufout,'(3(a,f12.4),/,3(a,f12.4))') ' zc_rol=',zc_rol,', r_sta/end=',rsta,',',rend,   &
                           ' cref_x=',cref_x,', z_sta/end=',zsta,',',zend
            call write_log(2, bufout)
         endif
      endif

      if (idebug.ge.2) then
         write(bufout,'(a,i2,a,2f12.6)') ' patch',icp,': start of interpen.region: (y,z)_tr=', cp%ysta, zsta
         call write_log(1, bufout)
         write(bufout,'(a,2f12.6)') '            end of interpen.region: (y,z)_tr=', cp%yend, zend
         call write_log(1, bufout)
      endif

      ! transform [ysta,zsta] to contact coordinates, store sp_sta-position

      msta = marker_2loc( marker(0d0, cp%ysta, zsta), cp%mref )
      mend = marker_2loc( marker(0d0, cp%yend, zend), cp%mref )
      cp%sp_sta = msta%y()
      cp%sp_end = mend%y()

      if (idebug.ge.2) then
         call marker_print( msta, 'msta(cp)', 1)
         call marker_print( mend, 'mend(cp)', 1)
         call marker_print( cp%mref, 'mref(tr)', 1)
      endif

      ! define planar potential contact area for [xsta_tr,xend_tr] x [sp_sta,sp_end]
      !  - there is an element with center at the contact reference position, (0,0) in local coords
      !  - elements are numbered here as [ ix_l : ix_h ] x [ iy_l : iy_h ],
      !             with ix_l and iy0 <= 0, ix_h and iy_h >= 0.
      !  - we get mx = ix_h - ix_l + 1, my = iy_h - iy_l + 1
      !  - element sizes are repeatedly doubled as needed to bring mx * my below npot_max.
      !  - effective element sizes are dx_eff = dx_fac * dx_in, ds_eff = ds_fac * ds_in
      !  - the limits [xp_l, xp_h] are rounded downwards/upwards to multiples of dx_eff
      !  - one element +/-dx_eff and +/-ds_eff is added as safety, guard band
      !  - +/- 0.5 * dx_eff is added for the corners of the elements

      xp_l = cp%xsta - cp%mref%x()
      xp_h = cp%xend - cp%mref%x()

      if (xp_l.gt.0d0 .or. xp_h.lt.0d0) then
         write(bufout,'(a,2g12.4)') ' Internal error: negative size for potential contact, xl/h=',xp_l,xp_h
         call write_log(1, bufout)
         xp_l = min(0d0, xp_l)
         xp_h = max(0d0, xp_h)
      endif

      if (cp%sp_sta.gt.cp%sp_end) then
         write(bufout,'(a,2g12.4)') ' Internal error: negative size for potential contact, sl/h=',      &
                cp%sp_sta, cp%sp_end
         call write_log(1, bufout)
         tmp       = cp%sp_sta
         cp%sp_sta = cp%sp_end
         cp%sp_end = tmp
      endif

      ! determine mx, my on basis of dx_eff, ds_eff

      cp%dx_fac = 1d0
      cp%ds_fac = 1d0
      cp%dx_eff = cp%dx_fac * dx_in
      cp%ds_eff = cp%ds_fac * ds_in

      ix_l = floor  (xp_l/cp%dx_eff) - 1
      ix_h = ceiling(xp_h/cp%dx_eff) + 1
      mx   = ix_h - ix_l + 1

      iy_l = floor  (cp%sp_sta/cp%ds_eff) - 1
      iy_h = ceiling(cp%sp_end/cp%ds_eff) + 1
      my   = iy_h - iy_l + 1

      ! increase dx_eff, ds_eff as needed to bring npot below max #elements

      if (ic%return.le.1 .and. npot_max.ge.100 .and. mx*my.ge.npot_max) then
         if (idebug.ge.-1) then
            write(bufout,'(a,2i5,a,i6,a,2f7.3,4(a,i4),a)') ' setting mx,my=', mx, my, ', npot=',        &
                mx*my, ', setting  dx,dy=', cp%dx_eff, cp%ds_eff, ', ix=[',ix_l,' :',ix_h,              &
                '], iy=[',iy_l,' :', iy_h,']'
            call write_log(1, bufout)
         endif
      
         do while (npot_max.ge.100 .and. mx*my.ge.npot_max)
            if (mx.gt.10) then
               cp%dx_fac = cp%dx_fac * 2d0
               cp%dx_eff = cp%dx_fac * dx_in
               ix_l = floor  (xp_l/cp%dx_eff) - 1
               ix_h = ceiling(xp_h/cp%dx_eff) + 1
               mx   = ix_h - ix_l + 1
            endif

            if (my.gt.10) then
               cp%ds_fac = cp%ds_fac * 2d0
               cp%ds_eff = cp%ds_fac * ds_in
               iy_l = floor  (cp%sp_sta/cp%ds_eff) - 1
               iy_h = ceiling(cp%sp_end/cp%ds_eff) + 1
               my   = iy_h - iy_l + 1
            endif

            if (idebug.ge.-1) then
               write(bufout,'(a,2i5,a,i6,a,2f7.3,4(a,i4),a)') ' reduced mx,my=', mx, my, ', npot=',     &
                   mx*my, ', enlarged dx,dy=', cp%dx_eff, cp%ds_eff, ', ix=[',ix_l,' :',ix_h,           &
                   '], iy=[',iy_l, ' :',iy_h,']'
               call write_log(1, bufout)
            endif
         enddo
      endif

      ! compute extent [xl,xh] x [sl,sh] of potential contact area

      xp_l = (ix_l - 0.5d0) * cp%dx_eff
      xp_h = (ix_h + 0.5d0) * cp%dx_eff
      sp_l = (iy_l - 0.5d0) * cp%ds_eff
      sp_h = (iy_h + 0.5d0) * cp%ds_eff

      if (idebug.ge.3) then
         write(bufout,128) 'range x =',     cp%mref%x(),' + [', xp_l, ',', xp_h, '], mx=',mx
         call write_log(1, bufout)
         write(bufout,128) 'range y =', sgn*cp%mref%y(),' + [', sgn*(cp%ysta-cp%mref%y()), ',',         &
                                                                sgn*(cp%yend-cp%mref%y()), ']'
         call write_log(1, bufout)
         write(bufout,128) 'range sr=',     cp%sr_ref,  ' + [', sp_l, ',', sp_h, '], my=',my
         call write_log(1, bufout)
 128     format(1x,3(a,f11.6),a,i6)
      endif

      ! enlarge ranges [usta,uend], [vsta,vend] by same amount as [xsta,xend], [ssta,send]

      if (cp%usta.lt.cp%uend) then
         fac = (xp_h - xp_l) / (cp%xend - cp%xsta)
         du  = (fac - 1d0) * (cp%uend - cp%usta)
         cp%usta = cp%usta - 0.5d0 * du
         cp%uend = cp%uend + 0.5d0 * du

         fac = (sp_h - sp_l) / (cp%sp_end - cp%sp_sta)
         dv  = (fac - 1d0) * (cp%vend - cp%vsta)
         cp%vsta = cp%vsta - 0.5d0 * dv
         cp%vend = cp%vend + 0.5d0 * dv

         if (idebug.ge.3) then
            write(bufout,128) 'range u = [', cp%usta, ',', cp%uend, ']'
            call write_log(1, bufout)
            write(bufout,128) 'range v = [', cp%vsta, ',', cp%vend, ']'
            call write_log(1, bufout)
         endif
      endif

      ! overwrite [xsta,xend]_tr (extent of interpen.area) with [xp_l,xp_h]+xref (extent of pot.con)

      cp%xsta = cp%mref%x() + xp_l
      cp%xend = cp%mref%x() + xp_h

      ! store sp_l, sp_h in [sp_sta,sp_end]

      cp%sp_sta = sp_l
      cp%sp_end = sp_h

   end subroutine set_planar_potcon

!------------------------------------------------------------------------------------------------------------

   subroutine compute_curved_potcon( ic, ws, prr, mtrk, nom_radius, cp, idebug )
!--purpose: compute the curved reference curve for conformal contact calculation
      implicit none
!--subroutine arguments:
      type(t_ic),        intent(in)  :: ic
      type(t_wheelset),  intent(in)  :: ws
      type(t_grid),      intent(in)  :: prr         ! track coordinates
      type(t_marker),    intent(in)  :: mtrk        ! roller profile marker
      real(kind=8),      intent(in)  :: nom_radius  ! radius in case of a roller
      type(t_cpatch)                 :: cp
      integer,           intent(in)  :: idebug
!--local variables:
      logical,      parameter :: use_old_wnrm = .true.
      type(t_wheel),  pointer :: my_wheel
      type(t_marker)          :: whl_trk, mq_cp, mq_trk, mq_whl
      type(t_grid)            :: rsrf, wsrf
      type(t_gridfnc3)        :: rnrm, wnrm
      integer                 :: iter, iy_ref, iy, iy0, iy1, ky, my, ny, sub_ierror
      real(kind=8)            :: fac_sc, sr_sta, sr_end, sw_ref, z_axle, x_err, yi, dy, dz, ds_min,     &
                                 cref_x, rref, dzrol
      real(kind=8), dimension(:), allocatable :: si, xw

      ! set pointer to the active rail and wheel in the current configuration

      if (idebug.ge.3) call write_log(' compute_curved_potcon: starting')

      my_wheel  => ws%whl

      associate( prw  => my_wheel%prw%grd_data,                                                         &
                 csrf => cp%curv_ref,     cnrm => cp%curv_nrm,      calph => cp%curv_incln)

      if (idebug.ge.3) then
         write(bufout,'(3(a,f12.6),a)') ' y_sta,y_end=[', cp%ysta,',',cp%yend,']'
         call write_log(1, bufout)
      endif

      ! 1. get initial range of sc-coordinates: range [sr_sta,sr_end] - sr_ref, enlarged by 20%

      call spline_get_s_at_y( prr%spl, cp%ysta, sr_sta, sub_ierror )
      call spline_get_s_at_y( prr%spl, cp%yend, sr_end, sub_ierror )

      ! at very small gaps, the range [ysta,yend] may miss the reference yicp

      sr_sta = min(sr_sta, cp%sr_ref)
      sr_end = max(sr_end, cp%sr_ref)

      fac_sc = 1.2d0
      cp%sc_sta = fac_sc * (sr_sta - cp%sr_ref)
      cp%sc_end = fac_sc * (sr_end - cp%sr_ref)

      if (idebug.ge.3) then
         write(bufout,'(3(a,f12.6),a)') ' sr_ref=',cp%sr_ref, ', sr_sta,sr_end=[', sr_sta,',',sr_end,']'
         call write_log(1, bufout)
         write(bufout,'(3(a,f12.6),a)') ' ds_eff=',cp%ds_eff, ', sc_sta,sc_end=[', cp%sc_sta,',',       &
                cp%sc_end,']'
         call write_log(1, bufout)
      endif

      ! snap to grid points i*ds, with three extra points on each side
      ! TODO: replace own blas (idamin) by MKL functions

      cp%sc_sta = (nint(cp%sc_sta/cp%ds_eff) - 3) * cp%ds_eff
      cp%sc_end = (nint(cp%sc_end/cp%ds_eff) + 3) * cp%ds_eff
      ny        = nint((cp%sc_end-cp%sc_sta)/cp%ds_eff) + 1

      allocate( si(ny) )

      do iy = 1, ny
         si(iy) = cp%sc_sta + (iy-1) * cp%ds_eff
      enddo
      iy_ref = idamin(ny, si, 1)

      if (idebug.ge.2) then
         write(bufout,'(3(a,f12.6),a)') ' sr_ref=',cp%sr_ref, ', sr_sta,sr_end=[', sr_sta,',',sr_end,']'
         call write_log(1, bufout)
         write(bufout,'(3(a,f12.6),a,i4)') ' ds_eff=',cp%ds_eff, ', sc_sta,sc_end=[', cp%sc_sta,',',    &
                cp%sc_end, '], ny=', ny
         call write_log(1, bufout)
      endif
      if (idebug.ge.2) then
         write(bufout,'(2(a,2f12.6))') ' si=', si(1), si(2), '...', si(ny-1), si(ny)
         call write_log(1, bufout)
         write(bufout,'(a,i4,a,f12.6)') ' iy_ref=', iy_ref,', s=', si(iy_ref)
         if (idebug.ge.4) call write_log(1, bufout)
      endif

      ! 2.a fill rsrf: evaluate (yr(s), zr(s)) at positions rsrf%s == sr(i)

      call grid_create_curvil(rsrf, 1, ny)
      call reallocate_arr(rsrf%s_prf, rsrf%ntot)

      rsrf%s_prf(1:ny) = si(1:ny) + cp%sr_ref
      call spline_get_xyz_at_s( prr%spl, ny, rsrf%s_prf, sub_ierror, yout=rsrf%y, zout=rsrf%z )

      ! correction for z-offset on roller at shifted x-position

      if (ic%is_roller()) then
         cref_x = cp%mref%x()
         rref   = nom_radius + mtrk%z() - cp%mref%z()
         dzrol  = rref - sqrt(rref**2 - cref_x**2) 
         rsrf%z(1:ny) = rsrf%z(1:ny) + dzrol

         if (idebug.ge.2) then
            write(bufout,'(4(a,f12.4))') ' rref=',rref,', cref_x=',cref_x,', dz=',dzrol
            call write_log(2, bufout)
         endif
      endif

      ! 2.b fill rnrm: evaluate dydz(s) at positions sr(i)
      !     tangent vector: [ (dx/ds), dy/ds, dz/ds ], outward normal vector [ 0, dz/ds, -dy/ds ]

      call gf3_new(rnrm, 'rnrm', rsrf, lzero=.true.)
      call spline_get_nvec_at_s_gfout( prr%spl, rnrm, sub_ierror )

      ! 3. compute wsrf: points (y(s),z(s)) on wheel surface at positions wsrf%s_prf == sw(i)

      ! 3.a compute contact point on wheel in wheel coords

      if (.true.) then
         mq_trk   = cp%mref
      else
         ! TODO: kin%pen is unknown at this point
         write(bufout,*) 'pen=',cp%gd%kin%pen
         call write_log(1,bufout)

         mq_cp  = marker( 0d0, 0d0, cp%gd%kin%pen )
         mq_trk = marker_2glob( mq_cp, cp%mref )
      endif
      whl_trk = marker_2glob( my_wheel%m_ws, ws%m_trk )
      mq_whl  = marker_2loc( mq_trk, whl_trk )

      ! get sw-position of contact reference mQ on wheel profile

      call spline_get_s_at_y( prw%spl, mq_whl%y(), sw_ref, sub_ierror )

      if (idebug.ge.5) then
         call grid_print(prw, 'prw(rw)', 5)
         write(bufout,*) 'sw_ref=',sw_ref
         call write_log(1, bufout)
      endif

      ! set initial estimate for positions (sw(i), xw(i))

      allocate( xw(ny) )
      call grid_create_curvil(wsrf, 1, ny)
      call reallocate_arr(wsrf%s_prf, wsrf%ntot)

      do iter = 1, 3

         if (iter.le.1) then
            wsrf%s_prf(1:ny) = sw_ref - si(1:ny) ! note: sw is oriented opposite to sr, si
            xw(1:ny) = mq_whl%x()
         else
            ! update estimate for position xw(i). (TODO: refine estimate of sw?)
            xw(1:ny) = xw(1:ny) + ( mq_trk%x() - wsrf%x(1:ny) )
         endif

         ! evaluate (yw(s), zw(s)) at positions sw(i) below axle (wheel profile coordinates)

         wsrf%x(1:ny) = xw(1:ny)
         call spline_get_xyz_at_s( prw%spl, ny, wsrf%s_prf, sub_ierror, yout=wsrf%y, zout=wsrf%z )

         ! compute zw(s) at positions (xw(s),yw(s)) (overwriting wsrf%z permitted because nx==1.)

         z_axle = -ws%nom_radius
         call grid_revolve_profile(wsrf, 1, ny, wsrf%x, z_axle, -999d0, wsrf%z)

         ! fill wnrm: evaluate dy/ds, dz/ds at positions sw(i), define outward normal

         if (use_old_wnrm) then
            call gf3_new(wnrm, 'wnrm', wsrf, lzero=.true.)
            call spline_get_nvec_at_s_gfout( prw%spl, wnrm, sub_ierror )
         endif

         ! rotate profile from wheel to track coordinates
         ! TODO: rotation should ignore yaw angle or set nx=0 and re-normalize

         if (abs(whl_trk%yaw()).gt.0.025d0) then
            write(bufout,'(2a,g12.4,a)') ' WARNING: curved reference surface may be inaccurate for ',    &
                        'yaw angle', whl_trk%yaw(),' rad'
            call write_log(1, bufout)
         endif

         call cartgrid_2glob( wsrf, whl_trk )
         if (use_old_wnrm) call gf3_rotate( wnrm, whl_trk%rot )

         if (idebug.ge.5) then
            call grid_print(wsrf, 'wsrf(trk)', 5, 9)
         endif

         ! compute mismatch in x-positions xw(s) vs. target x_cp

         x_err = 0d0
         do iy = 1, ny
            x_err = max(x_err, abs(wsrf%x(iy) - mq_trk%x()))
         enddo

         ! compute smallest length of segments projected on Oyz

         ds_min = 999d0
         do iy = 2, ny
            dy     = wsrf%y(iy) - wsrf%y(iy-1)
            dz     = wsrf%z(iy) - wsrf%z(iy-1)
            ds_min = min(ds_min, sqrt( dy**2 + dz**2 ) )
         enddo

         if (idebug.ge.2) then
            write(bufout,'(a,i2,2(a,f12.6))') ' wsrf: iter=',iter,': dx_max=',x_err,', ds_min=',ds_min
            call write_log(1, bufout)
         endif

      enddo

      ! compute surface normals wnrm: evaluate dydz(s) at positions sw(i)
      ! outward normal: [ 0, -dy/dz, 1 ]

      if (.not.use_old_wnrm) then
         wsrf%s_prf = -wsrf%s_prf                  ! temporarily use -sw to get increasing sequence
         ! call write_log('grid_make_spline wsrf...')
         call grid_make_ppspline(wsrf, lambda=0d0, use_wgt=.false., my_ierror=sub_ierror)
         call gf3_new(wnrm, 'wnrm', wsrf, lzero=.true.)
         call spline_get_dzdy_at_s(wsrf%spl, wnrm, sub_ierror)
         call gf3_scal(AllElm, -1d0, wnrm, ikYDIR)
         call gf3_set(AllElm,   1d0, wnrm, ikZDIR)
         wsrf%s_prf = -wsrf%s_prf                  ! reset original sw
      endif

      if (idebug.ge.3) then
         call grid_print(rsrf, 'rsrf_trk', 5, 9)
         call grid_print(wsrf, 'wsrf_trk', 5, 9)
         call gf3_print(rnrm, 'rnrm_trk', ikALL, 4)
         call gf3_print(wnrm, 'wnrm_trk', ikALL, 4)
      endif

      ! 4. trim curved reference: select indices iy for which csrf%y lies in [ysta, yend]

      iy0 = iy_ref-1
      iy1 = iy_ref+1
      do iy = 1, ny
         yi = 0.5 * (rsrf%y(iy) + wsrf%y(iy))
         if (yi.ge.cp%ysta .and. yi.le.cp%yend) then
            iy0 = min(iy0, iy)
            iy1 = max(iy1, iy)
         endif
      enddo

      ! 4.b add two columns on either side

      iy0 = max( 1, iy0-2)
      iy1 = min(ny, iy1+2)
      my  = max(0, iy1 - iy0 + 1)

      if (idebug.ge.2) then
         write(bufout,'(3(a,i4))') ' curved ref: selecting iy=[',iy0,':', iy1,'] of 1:',ny
         call write_log(1,bufout)
      endif
      if (iy0.gt.ny .or. iy1.lt.1) then
         write(bufout,'(a,f6.3,a)') ' Internal ERROR: no grid points selected on curved reference.'
         call write_log(1, bufout)
         call abort_run()
      endif
      if ((iy0.le.1 .or. iy1.ge.ny) .and. idebug.ge.2) then
         write(bufout,'(a,f6.3,a)') ' WARNING: no slack in sc_sta, sc_end, fac_sc=',fac_sc,             &
                ' should be enlarged.'
         call write_log(1, bufout)
      endif

      ! overwrite sc_sta, sc_end with selected range for potential contact area

      if (my.ge.1) then
         cp%sc_sta = si(iy0) - 0.5*cp%ds_eff
         cp%sc_end = si(iy1) + 0.5*cp%ds_eff
      endif

      ! 5. create output grid with my elements: (rsrf + wsrf) / 2

      call grid_create_curvil(csrf, 1, my)
      do iy = 1, my
         ky = iy0 - 1 + iy
         csrf%x(iy) = 0d0
         csrf%y(iy) = 0.5 * (rsrf%y(ky) + wsrf%y(ky))
         csrf%z(iy) = 0.5 * (rsrf%z(ky) + wsrf%z(ky))
      enddo

      ! create normal vectors on csrf: (rnrm - wnrm) / 2

      call gf3_new(cnrm, 'cnrm', csrf)
      do iy = 1, my
         ky = iy0 - 1 + iy
         cnrm%vx(iy) = 0.5d0 * (rnrm%vx(ky) - wnrm%vx(ky))
         cnrm%vy(iy) = 0.5d0 * (rnrm%vy(ky) - wnrm%vy(ky))
         cnrm%vn(iy) = 0.5d0 * (rnrm%vn(ky) - wnrm%vn(ky))
      enddo

      if (idebug.ge.3) then
         call grid_print(csrf, 'csrf_trk', 5, 9)
         call gf3_print(cnrm, 'cnrm_trk', ikALL, 4)
      endif

      ! create surface inclination on csrf: atan2(ny, -nz)
      ! outward normal has -nz>0;  ny=0 --> angle=0;  ny>0 --> angle>0

      call gf3_new(calph, 'calph', csrf)
      do iy = 1, my
         calph%vx(iy) = 0d0
         calph%vy(iy) = atan2( cnrm%vy(iy), -cnrm%vn(iy) )
         calph%vn(iy) = 0d0
      enddo

      if (idebug.ge.3) then
         call gf3_print(calph, 'calph_trk', ikALL, 4)
      endif

      call grid_destroy(rsrf)
      call grid_destroy(wsrf)
      call gf3_destroy(rnrm)
      call gf3_destroy(wnrm)
      deallocate(si, xw)

      end associate

   end subroutine compute_curved_potcon

!------------------------------------------------------------------------------------------------------------

   subroutine estimate_curvatures_v1(x, y, gap, nx, ny, ii_min, nom_radius, a1, b1, idebug)
!--purpose: estimate the longitudinal and lateral curvatures for a quick Hertzian assessment
!           Note: gap-mesh has size 1 x ny when using the contact locus approach
      implicit none
!--subroutine arguments:
      real(kind=8), intent(in)  :: x(nx*ny), y(nx*ny)
      real(kind=8), intent(in)  :: gap(nx*ny)    ! gap function
      integer,      intent(in)  :: nx, ny        ! grid dimensions
      integer,      intent(in)  :: idebug, ii_min
      real(kind=8), intent(in)  :: nom_radius    ! nominal wheel radius
      real(kind=8), intent(out) :: a1, b1
!--local variables:
      integer      :: ix_min, iy_min, my_min, kx, ky, ii_n, ii_s, ii_e, ii_w
      real(kind=8) :: dx, dy, rx, ry

      ! decompose ii_min into (ix,iy)

      ix_min = mod(ii_min-1, nx) + 1
      iy_min = (ii_min-1) / nx + 1

      ! shift minimum by one position if found on the boundary

      if (nx.gt.1) ix_min = max(2, min(nx-1, ix_min))
      iy_min = max(2, min(ny-1, iy_min))
      my_min = ix_min + (iy_min-1) * nx

      ! take different quotient over 2*k points to reduce sensitivity to local fluctations
      ! determine number of segments used on either side of central point

      kx = min(3, ix_min-1, nx-ix_min)
      ky = min(3, iy_min-1, ny-iy_min)

      if (idebug.ge.3) then
         write(bufout,'(4(a,i4))') ' estimate_curvatures_v1: ix_min=',ix_min,', kx=',kx,           &
                ', iy_min=',iy_min,', ky=',ky
         call write_log(1, bufout)
      endif

      ! compute indices of points north (iy+k), south (iy-k), east (ix+k) and west (ix-k)

      ii_n = my_min + ky*nx
      ii_s = my_min - ky*nx
      ii_e = my_min + kx
      ii_w = my_min - kx

      ! compute grid sizes. x,y are assumed to be in increasing order.

      dx = (x(ii_e) - x(ii_w)) / 2d0
      dy = (y(ii_n) - y(ii_s)) / 2d0

      if (idebug.ge.3) then
         write(bufout,'(2(a,f10.6))') ' estimate_curvatures_v1: grid sizes dx=',dx,', dy=',dy
         call write_log(1, bufout)
      endif

      ! using quadratic expansion: gap(x,y) - gap(my_min) =~= a1 * (x-x_min)^2 + b1 * (y-y_min)^2
      ! using only direct neighbouring points for now
      ! no support for minimum at the boundaries of the gap function

      if (nx.le.1) then
         a1 = 1d0 / (2d0 * nom_radius)
         rx = nom_radius
      else
         a1 = (gap(ii_e) + gap(ii_w) - 2d0*gap(my_min)) / (2d0*dx**2)
         rx = 0.5d0 / a1
         if (idebug.ge.3) then
            write(bufout,'(a,3e14.6,a)') ' gap east/ mid/west  =', gap(ii_e), gap(my_min), gap(ii_w), &
                ' * [1 -2 1] / 2*dx**2'
            call write_log(1, bufout)
         endif
      endif

      b1 = (gap(ii_n) + gap(ii_s) - 2d0*gap(my_min)) / (2d0*dy**2)
      if (abs(b1).gt.8d-5) then
         ry = 0.5d0 / b1
      else
         ry = 6250d0
      endif

      if (idebug.ge.2) then
         write(bufout,'(a,3e14.6,a)') ' gap north/mid/south =', gap(ii_n), gap(my_min), gap(ii_s),    &
                ' * [1 -2 1] / 2*dy**2'
         call write_log(1, bufout)
         write(bufout,'(a,e12.4,a,f10.3)') ' a1=', a1, ', rx=', rx
         call write_log(1, bufout)
         write(bufout,'(a,e12.4,a,f10.3)') ' b1=', b1, ', ry=', ry
         call write_log(1, bufout)
      endif

      ! Perform sanity checks

      if (rx.lt.9.9d0) then
         if (idebug.ge.1) then
            write(bufout,'(a,g12.4)') ' v1: The effective radius rx seems unrealistically small:', rx
            call write_log(1, bufout)
         endif
      elseif (rx.gt.3001d0) then
         if (idebug.ge.1) then
            write(bufout,'(a,g12.4)') ' v1: The effective radius rx seems unrealistically large:', rx
            call write_log(1, bufout)
         endif
      endif

      if (ry.lt.0.001d0) then
         if (idebug.ge.1) then
            write(bufout,'(a,g12.4,a)') ' v1: Small or negative ry found (',ry,'), using positive value instead'
            call write_log(1, bufout)
         endif
         ry = max(ry, rx/30d0)
         b1 = 0.5d0 / ry
      elseif (ry.gt.6001d0) then        ! Rw=10, Rr=-10.02 --> 5010
         if (idebug.ge.1) then
            write(bufout,'(a,g12.4,a)') ' v1: The effective radius ry seems unrealistically large (', ry,   &
                   '), truncating'
            call write_log(1, bufout)
         endif
         ry = 5000d0
         b1 = 0.5d0 / ry
      endif

   end subroutine estimate_curvatures_v1

!------------------------------------------------------------------------------------------------------------

   subroutine estimate_curvatures_v2( ic, cp, ws, prr, roller_radius, idebug )
!--purpose: estimate curvatures for a contact patch
      implicit none
!--subroutine arguments:
      type(t_ic)                :: ic
      type(t_cpatch)            :: cp
      type(t_wheelset)          :: ws
      type(t_grid), intent(in)  :: prr
      integer,      intent(in)  :: idebug
      real(kind=8), intent(in)  :: roller_radius
!--local variables:
      integer        :: sub_ierror
      real(kind=8)   :: rrx, rwx, sw_ref, sarr(2), dzdy(2), alpha(2), kappa_r, kappa_w, a1, b1, rx, ry
      type(t_marker) :: mq_trk, whl_trk, mq_whl
      type(t_wheel),  pointer :: my_wheel

      my_wheel  => ws%whl

      ! longitudinal: using nominal radius of wheel (and roller), corrected for contact angle
      ! (possible improvement: using actual instead of nominal radius?)

      rwx  = ws%nom_radius / max(0.05d0, cos(cp%delttr)) ! cos 0.05 == 87.1 deg
      a1   = 1d0 / (2d0 * rwx)
      !write(bufout,'(a,f8.4,2(a,f8.1),a,f8.4)') ' delttr=',cp%delttr,', rnom=',ws%nom_radius,           &
      !          ', rwx=',rwx,', a1=',a1
      !call write_log(1, bufout)

      if (ic%is_roller()) then
         rrx  = roller_radius / max(0.05d0, cos(cp%delttr))
         a1   = a1 + 1d0 / (2d0 * rrx)
      endif

      ! compute contact point on wheel in wheel coords, get sw-position

      mq_trk  = cp%mref         ! TODO: kin%pen is unknown at this point
      whl_trk = marker_2glob( my_wheel%m_ws, ws%m_trk )
      mq_whl  = marker_2loc( mq_trk, whl_trk )

      if (my_wheel%prw%grd_data%spl%nsec_top.le.0) then
         if (idebug.ge.2) call write_log(' prw%grd_data: add topview')
         call spline_add_topview(my_wheel%prw%grd_data%spl, .false., sub_ierror)
      endif
      call spline_get_s_at_y( my_wheel%prw%grd_data%spl, mq_whl%y(), sw_ref, sub_ierror )

      ! lateral: compute curvatures of profiles at contact reference position

      if (.true.) then

         ! method 1: \delta s = \delta\alpha/\kappa over full width of interpenetration area

         sarr(1)  = cp%sr_ref + cp%sp_sta
         sarr(2)  = cp%sr_ref + cp%sp_end
         call spline_get_dzdy_at_s(prr%spl, 2, sarr, dzdy, sub_ierror)
         alpha(1) = atan(dzdy(1))
         alpha(2) = atan(dzdy(2))
         kappa_r  = (alpha(2) - alpha(1)) / (sarr(2) - sarr(1))

         if (idebug.ge.2) then
            write(bufout,'(a,2f9.3,a,2f8.4,a,2g12.4)') ' prr_grd: s=',sarr(1), sarr(2),', alpha=',      &
                alpha(1), alpha(2), ', kappa_r=',kappa_r
            call write_log(1, bufout)
         endif

         sarr(1)  = sw_ref - cp%sp_sta  ! note: sw is oriented opposite to sr, sp
         sarr(2)  = sw_ref - cp%sp_end
         call spline_get_dzdy_at_s(my_wheel%prw%grd_data%spl, 2, sarr, dzdy, sub_ierror)
         alpha(1) = atan(dzdy(1))
         alpha(2) = atan(dzdy(2))
         kappa_w  = (alpha(2) - alpha(1)) / (sarr(2) - sarr(1))

         if (idebug.ge.2) then
            write(bufout,'(a,2f9.3,a,2f8.4,a,2g12.4)') ' prw_grd: s=',sarr(1), sarr(2),', alpha=',      &
                alpha(1), alpha(2), ', kappa_w=',kappa_w
            call write_log(1, bufout)
         endif

      else

         ! method 2: \kappa at sr_ref, sw_ref according to spline

         call spline_get_curv_yz_at_s( prr%spl, cp%sr_ref, kappa_r, sub_ierror )
         call spline_get_curv_yz_at_s( my_wheel%prw%grd_data%spl, sw_ref, kappa_w, sub_ierror )

         if (idebug.ge.2) then
            write(bufout,'(2(a,g12.4))') 'prr_spl: kappa_r=',kappa_r, ', kappa_w=',kappa_w
            call write_log(1, bufout)
         endif

      endif

      b1 = kappa_r / 2d0 + kappa_w / 2d0

      ! Perform sanity checks

      rx = 0.5d0 / a1
      ry = 0.5d0 / b1

      if (rx.lt.9.9d0) then
         if (idebug.ge.2) then
            write(bufout,'(a,g12.4)') ' The effective radius rx seems unrealistically small:', rx
            call write_log(1, bufout)
         endif
      elseif (rx.gt.4001d0) then        ! rnom=500, cos=0.08
         if (idebug.ge.2) then
            write(bufout,'(a,g12.4)') ' The effective radius rx seems unrealistically large:', rx
            call write_log(1, bufout)
         endif
      endif

      if (ry.lt.0.001d0) then
         if (idebug.ge.2) then
            write(bufout,'(a,g12.4,a)') ' Small or negative ry found (',ry,'), using positive value instead'
            call write_log(1, bufout)
         endif
         ry = max(ry, rx/30d0)
         b1 = 0.5d0 / ry
      elseif (ry.gt.6001d0) then        ! Rw=10, Rr=-10.02 --> 5010
         if (idebug.ge.2) then
            write(bufout,'(a,g12.4,a)') ' The effective radius ry seems unrealistically large (', ry,   &
                   '), truncating'
            call write_log(1, bufout)
         endif
         ry = 5000d0
         b1 = 0.5d0 / ry
      endif

      ! store results in problem datastructure

      ws%a1 = a1
      ws%b1 = b1

      if (idebug.ge.2) then
         write(bufout,'(a,g12.4,a,f10.3)') ' new_curvatures: a1=',a1,', rx=', 0.5d0/a1
         call write_log(1, bufout)
         write(bufout,'(a,g12.4,a,f10.3)') '                 b1=',b1,', ry=', 0.5d0/b1
         call write_log(1, bufout)
      endif

   end subroutine estimate_curvatures_v2

!------------------------------------------------------------------------------------------------------------

   subroutine estimate_minimum(x, y, gap, zr, ny, iy_min, xmin, ymin, gmin, zrmin, idebug)
!--purpose: determine local minimum of gap-function close to index iy_min
      implicit none
!--subroutine arguments:
      real(kind=8), intent(in)  :: x(ny), y(ny), gap(ny), zr(ny)
      integer,      intent(in)  :: ny
      integer,      intent(in)  :: idebug, iy_min
      real(kind=8), intent(out) :: xmin, ymin, gmin, zrmin
!--local variables:
      integer      :: iy
      real(kind=8) :: dy, a, b, c

      iy = max(1, min(ny, iy_min))

      if (iy.eq.1 .or. iy.eq.ny) then

         ! use value at grid point if minimum found on the boundary

         xmin  = x(iy)
         ymin  = y(iy)
         gmin  = gap(iy)
         zrmin = zr(iy)

      elseif (gap(iy-1).ge.888d0 .or. gap(iy+1).ge.888d0) then

         ! use value at grid point if no data provided at adjacent points

         xmin  = x(iy)
         ymin  = y(iy)
         gmin  = gap(iy)
         zrmin = zr(iy)

      else

         ! use quadratic fit

         dy   = y(iy) - y(iy-1)
         a    = (gap(iy+1) - 2d0*gap(iy) + gap(iy-1)) / (2d0 * dy**2)
         b    = (gap(iy+1)               - gap(iy-1)) / (2d0 * dy)
         c    =                  gap(iy)

         ymin = -b / (2d0*a)
         gmin = a * ymin**2 + b * ymin + c

         if (idebug.ge.3) then
            write(bufout,'(2(4(a,f12.6),:,/))') ' quadratic fit: y_{i-1}=', y(iy-1), ', y_i=', y(iy),   &
                ', y_{i+1}=', y(iy+1), ', ymin=', ymin, '                g_{i-1}=', gap(iy-1),          &
                ', g_i=', gap(iy), ', g_{i+1}=', gap(iy+1), ', gmin=', gmin
            call write_log(2, bufout)
         endif

         a    = (zr(iy+1) - 2d0*zr(iy) + zr(iy-1)) / (2d0 * dy**2)
         b    = (zr(iy+1)              - zr(iy-1)) / (2d0 * dy)
         c    =                 zr(iy)
         zrmin = a * ymin**2 + b * ymin + c

         a    = (x(iy+1) - 2d0*x(iy) + x(iy-1)) / (2d0 * dy**2)
         b    = (x(iy+1)             - x(iy-1)) / (2d0 * dy)
         c    =                x(iy)
         xmin = a * ymin**2 + b * ymin + c

         ymin = ymin + y(iy)

      endif

   end subroutine estimate_minimum

!------------------------------------------------------------------------------------------------------------

   subroutine find_gap_locmin(x, y, zr, zw, sgn, gap, nx, ny, nummin, locmin, mask, idebug, len_locmin)
!--purpose: find the indices (ix,iy) where the array gap is lower than the four neighbouring points
!           excluding points that are on the same plateau +/- noise
      implicit none
!--subroutine arguments:
      real(kind=8), intent(in)  :: x(nx*ny), y(nx*ny), zr(nx*ny), zw(nx*ny)
      real(kind=8), intent(in)  :: sgn           ! +1/-1 for right/left w/r-combination
      real(kind=8), intent(in)  :: gap(nx*ny)    ! gap function
      integer,      intent(in)  :: nx, ny        ! grid dimensions
      integer,      intent(in)  :: idebug, len_locmin
      integer,      intent(out) :: nummin
      integer,      intent(out) :: locmin(2,len_locmin)
      logical                   :: mask(nx*ny)   ! work-space for plateau
!--local variables:
      real(kind=8), parameter   :: thresh = 1d-6
      integer      :: ix, iy, ii, ii_n, ii_s, ii_e, ii_w
      logical      :: is_locmin

      nummin = 0
      mask(1:nx*ny) = .false.

      if (idebug.ge.4) then
         write(bufout,*) 'starting find_gap_locmin, nx=',nx,', ny=',ny
         call write_log(1, bufout)

         do ix = 1, min(nx,10)
            write(bufout,'(a,i3,a,8f10.4)') ' ix=',ix,': ',(gap(ix+(iy-1)*nx), iy=1,min(ny,8))
            call write_log(1, bufout)
         enddo
      endif

      do ix = 1, nx
         do iy = 1, ny

            ! compute indices of points north (iy+1), south (iy-1), east (ix+1) and west (ix-1)

            ii   = ix + (iy-1) * nx
            ii_n = ii
            ii_s = ii
            ii_e = ii
            ii_w = ii
            if (iy.lt.ny) ii_n = ii + nx
            if (iy.gt. 1) ii_s = ii - nx
            if (ix.lt.nx) ii_e = ii + 1
            if (ix.gt. 1) ii_w = ii - 1

            ! exclude wheel points not on the rail surface (zw = -888.0, gap = zr + 888.0)
            ! 

            is_locmin = (gap(ii).lt.777d0)

            ! check south point iy-1

            if (is_locmin .and. iy.gt.1) then
               if (mask(ii_s) .and. abs(gap(ii_s)-gap(ii)).lt.thresh) then
                  mask(ii)  = .true.    ! expand plateau
                  is_locmin = .false.
               elseif (gap(ii_s).le.gap(ii)) then
                  is_locmin = .false.
               endif
            endif

            ! check north point iy+1

            if (is_locmin .and. iy.lt.ny) then
               if (mask(ii_n) .and. abs(gap(ii_n)-gap(ii)).lt.thresh) then
                  mask(ii)  = .true.    ! expand plateau
                  is_locmin = .false.
               elseif (gap(ii_n).lt.gap(ii)) then
                  is_locmin = .false.
               endif
            endif

            ! check west point ix-1

            if (is_locmin .and. ix.gt.1) then
               if (mask(ii_w) .and. abs(gap(ii_w)-gap(ii)).lt.thresh) then
                  mask(ii)  = .true.    ! expand plateau
                  is_locmin = .false.
               elseif (gap(ii_w).le.gap(ii)) then
                  is_locmin = .false.
               endif
            endif

            ! check east point ix+1

            if (is_locmin .and. ix.lt.nx) then
               if (mask(ii_e) .and. abs(gap(ii_e)-gap(ii)).lt.thresh) then
                  mask(ii)  = .true.    ! expand plateau
                  is_locmin = .false.
               elseif (gap(ii_e).lt.gap(ii)) then
                  is_locmin = .false.
               endif
            endif

            ! print information for all local minima

            if (is_locmin) then

               if (idebug.ge.2) then
                  write(bufout,124) ' found local minimum, gap =',gap(ii), ' at (',              &
                        x(ii), ',', sgn*y(ii), ',', zr(ii),')'
                  call write_log(1, bufout)
               endif
 124           format(a,f12.6,a,f11.6,2(a,f12.6),a)

               if (idebug.ge.3) then
                  write(bufout,*) 'ix=',ix,', iy=',iy
                  call write_log(1, bufout)
               endif

               if (idebug.ge.4) then
                   write(bufout,125) ' XS:',x(ii_s),  ', XE:',x(ii_e),  ', XN:',x(ii_n),  ', XW:',x(ii_w)
                   call write_log(1, bufout)
                   write(bufout,125) ' YS:',y(ii_s),  ', YE:',y(ii_e),  ', YN:',y(ii_n),  ', YW:',y(ii_w)
                   call write_log(1, bufout)
                   write(bufout,125) ' ZS:',zr(ii_s), ', ZE:',zr(ii_e), ', ZN:',zr(ii_n), ', ZW:',zr(ii_w)
                   call write_log(1, bufout)
                   write(bufout,125) ' ZS:',zw(ii_s), ', ZE:',zw(ii_e), ', ZN:',zw(ii_n), ', ZW:',zw(ii_w)
                   call write_log(1, bufout)
                   write(bufout,125) ' GS:',gap(ii_s),', GE:',gap(ii_e),', GN:',gap(ii_n),', GW:',gap(ii_w)
                   call write_log(1, bufout)
 125               format(4(a,f12.6))
               endif
            endif

            ! store position (ix,iy) of local minimum if gap is negative

            if (is_locmin .and. gap(ii).lt.0d0) then
               nummin = nummin + 1
               mask(ii) = .true.        ! start new plateau
               if (nummin.lt.len_locmin) then
                  locmin(1,nummin) = ix
                  locmin(2,nummin) = iy
               elseif (idebug.ge.1) then
                  write(bufout,*) 'INTERNAL ERROR: maximum number of local minima exceeded',nummin
                  call write_log(1, bufout)
               endif
            elseif (is_locmin) then
               if (idebug.ge.2) call write_log(' ...discarding local minimum, no interpenetration')
            endif

         enddo ! iy
      enddo ! ix

   end subroutine find_gap_locmin

!------------------------------------------------------------------------------------------------------------

   subroutine find_bounding_box_1d( mask, ny, iy_cntr, nguard, iy_sta, iy_end, idebug)
!--purpose: generic routine for determining the interval [iy_sta,iy_end] around the region where mask==true.
!            - [iy_cntr] >= 0 is the start point for the search
!              the smallest interval including [iy_cntr] is determined satisfying the guard condition
!            - nguard is the number of additional points in which mask==false included in the box
!              if nguard<=0 then mask==true occurs on the bounds of the interval.
!              if nguard>=1 then nguard points are added at start/end where mask==false
      implicit none
!--subroutine arguments:
      logical,      intent(in)  :: mask(ny)                   ! mask function
      integer,      intent(in)  :: ny                         ! grid dimensions
      integer,      intent(in)  :: nguard                     ! number of rows/columns in guard band
      integer,      intent(in)  :: iy_cntr                    ! start point of flooding algorithm
      integer,      intent(out) :: iy_sta, iy_end             ! extent of bounding box
      integer,      intent(in)  :: idebug
!--local variables:
      integer      :: iy
      logical      :: lfound

      if (nguard.ne.1) then
         call write_log(' ERROR: bound_box: nguard must be 1.')
         call abort_run()
      endif

      ! if no [iy_cntr] is given, we determine the maximum interval around all masked points

      if (iy_cntr.le.0) then

         ! initialize: empty box

         iy_sta = ny+1
         iy_end = 0

         ! just process all points once

         do iy = 1, ny
            if (mask(iy)) then
               iy_sta = min(iy_sta, iy)
               iy_end = max(iy_sta, iy)
            endif
         enddo

         ! extend non-empty box with guard band

         if (iy_sta.le.ny) then
            iy_sta = max( 1, iy_sta-nguard)
            iy_end = min(ny, iy_end+nguard)
         endif

      else

      ! if [iy_cntr] is given, we determine the smallest interval around this "center"

         ! initialize with interval of just one point

         iy_sta = iy_cntr
         iy_end = iy_cntr

         ! expand to the left, until no masked point at iy_sta-1

         lfound = .true.
         do while(iy_sta.gt.1 .and. lfound)
            lfound = mask(iy_sta-1)
            if (lfound) iy_sta = iy_sta - 1
         enddo

         ! expand to the right, until no masked point at iy_end+1

         lfound = .true.
         do while(iy_end.lt.ny .and. lfound)
            lfound = mask(iy_end+1)
            if (lfound) iy_end = iy_end + 1
         enddo

         ! extend the interval around the interpenetration region with one column on both sides

         if (iy_sta.gt.1)  iy_sta = iy_sta - 1
         if (iy_end.lt.ny) iy_end = iy_end + 1

         ! check: count number of masked points at first/last rows/columns

      endif ! center-point given

      if (idebug.ge.4) then
         write(bufout,'(4(a,i4),a)') ' find_bounding_box_1d: iy=[',iy_sta,',', iy_end,']'
         call write_log(1, bufout)
      endif

   end subroutine find_bounding_box_1d

!------------------------------------------------------------------------------------------------------------

   subroutine find_bounding_box_2d( mask, nx, ny, ix_cntr, iy_cntr, nguard, ix_sta, ix_end, iy_sta,     &
                                    iy_end, icheck, idebug)
!--purpose: generic routine for determining the rectangular box [ix_sta,ix_end] x [iy_sta,iy_end]
!           around the region where mask==true.
!            - [ix,iy_cntr] >= 0 is the start point for the search
!              the smallest rectangle including [ix,iy_cntr] is determined satisfying the guard condition
!            - nguard is the number of additional rows/columns in which mask==false included in the box
!              if nguard<=0 then there can be points with mask==true on the first/last rows/columns,
!              if nguard>=1 then nguard rows/columns are added at start/end where all mask==false
      implicit none
!--subroutine arguments:
      logical,      intent(in)  :: mask(nx,ny)                ! mask function
      integer,      intent(in)  :: nx, ny                     ! grid dimensions
      integer,      intent(in)  :: nguard                     ! number of rows/columns in guard band
      integer,      intent(in)  :: ix_cntr, iy_cntr           ! start point of flooding algorithm
      integer,      intent(out) :: ix_sta, ix_end, iy_sta, iy_end ! extent of bounding box
      integer,      intent(in)  :: icheck, idebug
!--local variables:
      integer      :: ix, iy, iboundx, iboundy
      logical      :: lchanged, lfound
      character    :: str(nx)

      if (nguard.ne.1) then
         call write_log(' ERROR: bound_box: nguard must be 1.')
         call abort_run()
      endif

      ! if no [ix_cntr,iy_cntr] are given, we determine the maximum bounding box around all masked points

      if (ix_cntr.le.0 .or. iy_cntr.le.0) then

         ! initialize: empty box

         ix_sta = nx+1
         ix_end = 0
         iy_sta = ny+1
         iy_end = 0

         ! just process all points once

         do iy = 1, ny
            do ix = 1, nx
               if (mask(ix,iy)) then
                  ix_sta = min(ix_sta, ix)
                  ix_end = max(ix_sta, ix)
                  iy_sta = min(iy_sta, iy)
                  iy_end = max(iy_sta, iy)
               endif
            enddo
         enddo

         ! extend non-empty box with guard band

         if (ix_sta.le.nx) then
            ix_sta = max( 1, ix_sta-nguard)
            ix_end = min(nx, ix_end+nguard)
            iy_sta = max( 1, iy_sta-nguard)
            iy_end = min(ny, iy_end+nguard)
         endif

      else

      ! if [ix_cntr,iy_cntr] are given, we determine the smallest bounding box around this "center"

         ! initialize with bounding box of just one point

         ix_sta = ix_cntr
         ix_end = ix_cntr
         iy_sta = iy_cntr
         iy_end = iy_cntr

         ! iterate until not extended anymore

         lchanged = .true.
         do while (lchanged)

            lchanged = .false.
            if (idebug.ge.6) call write_log(' bound_box: starting iteration...')

            ! expand to the west, left side of interval [ix_sta,ix_end], until no masked pts at ix_sta-1

            lfound = .true.
            do while(ix_sta.gt.1 .and. lfound)

               ! find masked pts at column ix_sta-1 including diagonal points [iy_sta-1,iy_end+1]

               lfound = .false.
               iy = max(1,iy_sta-1)
               do while(iy.le.min(ny,iy_end+1) .and. .not.lfound)
                  lfound = mask(ix_sta-1,iy)
                  iy = iy + 1
               enddo

               ! if found, move one column to the left

               if (lfound) then
                  ix_sta = ix_sta - 1
                  lchanged = .true.
               endif
            enddo

            ! expand to the east, right side of interval [ix_sta,ix_end], until no masked pts at ix_end+1

            lfound = .true.
            do while(ix_end.lt.nx .and. lfound)

               ! find masked points at column ix_end+1 including diagonal points [iysta-1, iyend+1]

               lfound = .false.
               iy = max(1,iy_sta-1)
               do while(iy.le.min(ny,iy_end+1) .and. .not.lfound)
                  lfound = mask(ix_end+1,iy)
                  iy = iy + 1
               enddo

               ! if found, move one column to the right

               if (lfound) then
                  ix_end = ix_end + 1
                  lchanged = .true.
               endif
            enddo

            ! expand to the south, bottom side of interval [iy_sta,iy_end], until no masked pts at iy_sta-1

            lfound = .true.
            do while(iy_sta.gt.1 .and. lfound)

               ! find masked points at row iy_sta-1

               lfound = .false.
               ix = max(1,ix_sta-1)
               do while(ix.le.min(nx,ix_end+1) .and. .not.lfound)
                  lfound = mask(ix,iy_sta-1)
                  ix = ix + 1
               enddo

               ! if found, move one row down

               if (lfound) then
                  iy_sta = iy_sta - 1
                  lchanged = .true.
               endif
            enddo

            ! expand to the north, top side of interval [iy_sta,iy_end], until no masked pts at iy_end+1

            lfound = .true.
            do while(iy_end.lt.ny .and. lfound)

               ! find masked points at row iy_end+1

               lfound = .false.
               ix = max(1,ix_sta-1)
               do while(ix.le.min(nx,ix_end+1) .and. .not.lfound)
                  lfound = mask(ix,iy_end+1)
                  ix = ix + 1
               enddo

               ! if found, move one row up

               if (lfound) then
                  iy_end = iy_end + 1
                  lchanged = .true.
               endif
            enddo

         enddo ! while(lchanged)

         ! extend the box around the interpenetration region by one row/column on all sides

         if (ix_sta.gt.1)  ix_sta = ix_sta - 1
         if (ix_end.lt.nx) ix_end = ix_end + 1
         if (iy_sta.gt.1)  iy_sta = iy_sta - 1
         if (iy_end.lt.ny) iy_end = iy_end + 1

         ! check: count number of masked points at first/last rows/columns

         if (icheck.ge.1 .and. idebug.ge.2) then
            iboundx = 0
            iboundy = 0

            do iy = iy_sta, iy_end
               do ix = ix_sta, ix_end
                  if (mask(ix,iy)) then
                     if (ix.eq.ix_sta .or. ix.eq.ix_end) iboundx = iboundx + 1
                     if (iy.eq.iy_sta .or. iy.eq.iy_end) iboundy = iboundy + 1
                  endif
               enddo
            enddo

            if (iboundx.gt.0) then
               write(bufout,'(3(a,i5))') ' ', iboundx,' masked points at ix',ix_sta,' or',ix_end
               call write_log(1, bufout)
            endif
            if (iboundy.gt.0) then
               write(bufout,'(3(a,i5))') ' ', iboundy,' masked points at iy',iy_sta,' or',iy_end
               call write_log(1, bufout)
            endif

            if (iboundx.gt.0 .or. iboundy.gt.0) then
               do iy = ny, 1, -1
                  str = ' '
                  do ix = 1, nx
                     if (mask(ix,iy)) then
                        str(ix) = '*'
                     else
                        str(ix) = ' '
                     endif
                     if ((ix.eq.1 .or. ix.eq.nx) .and. (iy.eq.1 .or. iy.eq.ny)) str(ix) = 'o'
                  enddo
                  write(bufout,2200) iy, (str(ix), ix=1,nx)
                  call write_log(1, bufout)
 2200             format (i6, ':', 1x, 900(a,:))
               enddo
            endif

         endif ! check

      endif ! center-point given

      if (idebug.ge.4) then
         write(bufout,'(4(a,i4),a)') ' find_bounding_box_2d: ix=[',ix_sta,',', ix_end,'], iy=[',iy_sta,',', &
                   iy_end,']'
         call write_log(1, bufout)
      endif

   end subroutine find_bounding_box_2d

!------------------------------------------------------------------------------------------------------------

   subroutine prw_determine_foldback(prw, max_omit, err_hnd, iwsta, iwend, idebug, my_ierror)
!--purpose: create trimmed version of wheel profile with uni-valued function
!           trim at vertical sections, starting at point with mean y-value
      implicit none
!--subroutine arguments:
      type(t_grid)              :: prw
      integer,      intent(in)  :: err_hnd, idebug
      real(kind=8), intent(in)  :: max_omit
      integer,      intent(out) :: iwsta, iwend, my_ierror
!--local variables:
      integer             :: iw, iymin, iymid, iymax, ikeep, inext, ndel
      real(kind=8)        :: ymid
      real(kind=8), dimension(:), allocatable :: dy_ds

      my_ierror = 0

      allocate(dy_ds(prw%ntot))
      call spline_eval(prw%spl, ikYDIR, prw%ntot, prw%s_prf, my_ierror, f1_eval=dy_ds)

      ! get mid-point 'iymid': point with middle y-value, on the tread of the wheel

      iymin = idmin(prw%ntot, prw%y, 1)
      iymax = idmax(prw%ntot, prw%y, 1)
      ymid  = (prw%y(iymin) + prw%y(iymax)) * 0.5d0

      iymid = 1
      do iw = 2, prw%ntot
         if (abs(prw%y(iw) - ymid).lt.abs(prw%y(iymid) - ymid)) iymid = iw
      enddo

      if (idebug.ge.2) then
         write(bufout,'(2(a,f8.3),a)') ' wheel profile has y in [', prw%y(iymin),',', prw%y(iymax),']'
         call write_log(1, bufout)
         write(bufout,'(a,f8.3,a,i5,a,f8.3)') ' wheel profile has y_mid=', prw%y(iymid), ' at iymid=',  &
                iymid,', dy_ds=',dy_ds(iymid)
         call write_log(1, bufout)
      endif

      ! loop over points i = iymid : 1 : ntot, to the left of iymid with decreasing y-value
      ! keep points as long as ynext < ycurr, stop at first point with ynext >= ycurr

      ikeep = iymid
      inext = min(prw%ntot, ikeep + 1)

      if (idebug.ge.4) then
         write(bufout,'(2(a,i5))') ' loop over points to the left, i=',ikeep,' to', prw%ntot
         call write_log(1, bufout)
      endif

      do while (inext.gt.ikeep)
         if (prw%y(inext).lt.prw%y(ikeep) .and. dy_ds(inext).lt.0d0) then
            ikeep = ikeep + 1
            if (inext.lt.prw%ntot) inext = inext + 1
         else
            if (idebug.ge.3) then
               write(bufout,'(2(a,i4,a,f9.4))') ' prw%y(',inext,')=',prw%y(inext),' >= prw%y(',         &
                                ikeep, ')=',prw%y(ikeep)
               call write_log(1, bufout)
            endif
            inext = -1
         endif
      enddo
      iwend = ikeep

      ! loop over points i = iymid : -1 : 1, to the right of iymid with increasing y-value
      ! keep points as long as ynext > ycurr, stop at first point with ynext <= ycurr

      ikeep = iymid
      inext = max(1, ikeep - 1)

      if (idebug.ge.4) then
         write(bufout,'(2(a,i5))') ' loop over points to the right,  i=',ikeep,' to', 1
         call write_log(1, bufout)
      endif

      do while (inext.lt.ikeep)
         if (prw%y(inext).gt.prw%y(ikeep) .and. dy_ds(inext).lt.0d0) then
            ikeep = ikeep - 1
            if (inext.gt.1) inext = inext - 1
         else
            if (idebug.ge.3) then
               write(bufout,'(2(a,i4,a,f9.4))') ' prw%y(',inext,')=',prw%y(inext),' <= prw%y(',         &
                                ikeep, ')=',prw%y(ikeep)
               call write_log(1, bufout)
            endif
            inext = 999999
         endif
      enddo
      iwsta = ikeep

      if (idebug.ge.4) then
         write(bufout,'(2(a,i5))') ' keeping points i=', iwsta,' to', iwend
         call write_log(1, bufout)
      endif
      if (idebug.ge.3 .and. (iwsta.gt.1 .or. iwend.lt.prw%ntot)) then
         write(bufout,'(a,i4,2a)') ' ignoring', prw%ntot-iwend+iwsta-1,' points from wheel profile where', &
                        ' y(i+1) >= y(i)'
         call write_log(1, bufout)
      endif

      if (idebug.ge.2) then
         if (iwsta.gt.1) then
            write(bufout,'(9x,i4,a)') iwsta-1,' points are deleted from the start of the profile'
            call write_log(1, bufout)
         endif
         if (iwend.lt.prw%ntot) then
            write(bufout,'(9x,i4,a)') prw%ntot-iwend,' points are deleted from the end of the profile'
            call write_log(1, bufout)
         endif
      endif

      ! error in case more than fraction max_omit of points is deleted from profile

      ndel = prw%ntot - (iwend - iwsta + 1)

      if (ndel.gt.nint(max_omit*prw%ntot)) then
         my_ierror = -2102
         write(bufout,'(a,i4,a,f6.2,a,i4)') 'ndel=',ndel, ', max_omit=',max_omit, ', ntot=',prw%ntot
            call write_log(1, bufout)
         if (err_hnd.ge.-1) then
            write(bufout,'(a,f6.1,a)') ' ERROR: wheel profile y-values are non-monotonic,',             &
                100d0*real(ndel)/real(prw%ntot), ' % of profile discarded'
            call write_log(1, bufout)
         endif
      endif
      deallocate(dy_ds)

   end subroutine prw_determine_foldback

!------------------------------------------------------------------------------------------------------------

end module m_wr_locatecp
