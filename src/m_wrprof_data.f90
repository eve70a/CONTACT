!------------------------------------------------------------------------------------------------------------
! m_wrprof_data - declare user-defined types for the wheel-set and track, grouping important variables
!                 into a hierarchical, tree-like data-structure.
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_wrprof_data
use m_hierarch_data
use m_profile
use m_varprof
implicit none
public

   public  t_discret

   public  t_region

   public  t_cpatch
   public  p_cpatch
   public  cp_init
   public  cp_destroy
   interface cp_destroy
      module procedure tpatch_destroy
      module procedure ppatch_destroy
   end interface cp_destroy

   public  t_wheel
   public  wheel_ini

   public  t_wheelset
   public  wheelset_ini

   public  t_rail
   public  rail_ini

   public  t_trackdata

   public  t_ws_track
   public  p_ws_track
   public  wrprof_ini
   public  wrprof_destroy

   !---------------------------------------------------------------------------------------------------------
   ! dimensions and parameters

   integer,      parameter :: MAX_NUM_RGN = 10          ! max 10 regions with possible contact per wheel
   integer,      parameter :: MAX_NUM_CPS = 10          ! max 10 contact patches per wheel

   !---------------------------------------------------------------------------------------------------------
   ! data with respect to the grid discretisation, potential contact area:

   type :: t_discret
      integer      :: npot_max
      real(kind=8) :: dx, ds
      real(kind=8) :: dqrel
      real(kind=8) :: angl_sep, dist_turn, dist_sep, dist_comb

      ! npot_max  [-]   maximum number of elements (npot=mx*my) permitted in potential contact area
      ! dx        [mm]  size of rectangular elements in rolling direction (x)
      ! ds        [mm]  size of rectangular elements in lateral direction (s/y)
      ! dqrel     [-]   rolling step size V*dt relative to dx
      ! angl_sep  [rad] angle threshold: contact patches are separated when |a1 - a2| >= angl_sep
      ! dist_turn [mm]  distance threshold: contact angles are turned when dist_sep <= dist <= dist_turn
      ! dist_sep  [mm]  upper dist. threshold: contact patches are separated when dist >= dist_sep
      ! dist_comb [mm]  lower dist. threshold: contact patches are combined when dist <= dist_comb

   end type t_discret

   !---------------------------------------------------------------------------------------------------------
   ! data with respect to a "possible contact region", i.e., search regions on wheel & rail surfaces

   type :: t_region
      type(t_marker)   :: mview
      real(kind=8)     :: xsta, xend
      real(kind=8)     :: sr_sta, sr_end, sw_sta, sw_end

      ! mview            marker w.r.t. track coordinates: local system with appropriate 'view direction'
      ! xsta/end   [mm]  search range in longitudinal direction in mview coordinates
      ! sr_sta/end [mm]  search range on the rail profile in rail arc lengths sr
      ! sw_sta/end [mm]  search range on the wheel profile in wheel arc lengths sw

   end type t_region

   !---------------------------------------------------------------------------------------------------------
   ! data with respect to a "contact patch", i.e., a separate contact problem for CONTACT:

   type :: t_cpatch
      type(t_vec)      :: micp
      real(kind=8)     :: gap_min, totgap, wgt_xgap, wgt_ygap, wgt_zgap, wgt_agap
      type(t_marker)   :: mref, mpot
      real(kind=8)     :: delttr
      real(kind=8)     :: xsta, xend, ysta, yend, zsta, zend, usta, uend, vsta, vend
      real(kind=8)     :: sr_ref, sp_sta, sp_end, sc_sta, sc_end
      real(kind=8)     :: dx_fac, ds_fac, dx_eff, ds_eff
      integer          :: nsub
      real(kind=8), dimension(:), pointer :: y_sep => NULL(), f_sep => NULL()
      integer          :: prev_icp(MAX_NUM_CPS)
      type(t_vec)      :: ftrk, ttrk, fws, tws
      type(t_grid)     :: rail_srfc, whl_srfc, curv_ref
      type(t_gridfnc3) :: curv_nrm, curv_incln
      type(t_probdata), pointer :: gd => NULL() ! pointer to a gd hierarchical data-structure
      logical          :: has_own_subs

      ! variables concerning the interpenetration region:
      ! micp             position of maximum interpenetration in contact patch wrt the track coordinate system
      ! gap_min   [mm]   minimum vertical gap for this contact patch
      ! totgap    [?]    integral of gap (squared) used for weighted center
      ! wgt_xgap  [mm]   weighted average of x-values for interpenetration area
      ! wgt_ygap  [mm]   weighted average of y-values for interpenetration area
      ! wgt_agap  [rad]  weighted average of alpha values for interpenetration area
      ! [xyz]sta  [mm]   start-position of (estimated) interpenetration area in track coordinates
      ! [xyz]end  [mm]   end-position of (estimated) interpenetration area in track coordinates
      ! [uv]sta   [-]    start-position of (estimated) interpenetration area in surface (u,v) coordinates
      ! [uv]end   [-]    end-position of (estimated) interpenetration area in surface (u,v) coordinates
      ! sr_ref    [mm]   sr-position of the contact reference in rail profile sr-coordinates
      !
      ! variables concerning the potential contact area:
      ! mref             contact reference marker w.r.t. track coordinates
      ! mpot             potential contact marker w.r.t. track coordinates (origin of local coordinates)
      ! delttr    [rad]  contact angle in track coordinates, rotation from track z-axis to contact n-axis
      ! sp_sta    [mm]   start-position of (estimated) interpenetration area in planar contact sp-coordinate
      ! sp_end    [mm]   end-position of (estimated) interpenetration area in planar contact sp-coordinate
      ! sc_sta    [mm]   start-position of (estimated) interpenetration area in conformal sc-coordinate
      ! sc_end    [mm]   end-position of (estimated) interpenetration area in conformal sc-coordinate
      ! dx_fac    [-]    multiplier for grid step size in x-direction, dx_eff = dx_fac * dx_in
      ! ds_fac    [-]    multiplier for grid step size in s-direction, ds_eff = ds_fac * ds_in
      ! dx_eff    [mm]   effective grid step size used in x-direction
      ! ds_eff    [mm]   effective grid step size used in s-direction
      !
      ! for combining contact patches with blending approach:
      ! nsub             number of sub-patches with reduced interaction
      ! y_sep     [mm]   track y-positions between neighbouring sub-patches in increasing order
      ! f_sep     [-]    weighting factors for interactions between neighbouring sub-patches
      !
      ! connection between this patch and patches at previous time instance
      ! prev_icp  [-]    contact patch numbers at previous time connected to this icp, 0 in unused positions
      !
      ! for pot.contact in conformal method:
      ! curv_ref  [mm]   curved reference [0,y,z]_c(tr) (1 x my) for conformal contact approach (trk coord)
      ! curv_nrm         non-normalized outward surface normals [0,ny,nz]_(tr) at curved reference surface
      ! curv_incln [rad] values representing surface inclination [0,alpha,0]_(tr) of curved reference surface
      !
      ! for use in wr_undefdistc and wr_rigslip:
      ! rail_srfc [mm]   rail surface [x,sp,np]_r(cp) in planar contact reference coordinates
      ! whl_srfc  [mm]   wheel surface [x,sp,np]_w(cp) in planar contact reference coordinates
      !                  note: sp_j may be non-uniform in the conformal approach where d sc=const
      !
      ! for actual computation:
      ! has_own_subs     cp has specific subs-data or use subs-data of overall problem?
      ! ftrk      [N]    total forces on rail in global track directions
      ! ttrk      [N.mm] total moments on rail about rail profile marker in global track directions
      ! fws       [N]    total forces on rail in wheelset coordinates
      ! tws       [N.mm] total moments on rail about wheel profile marker in wheelset coordinates

   end type t_cpatch

   type :: p_cpatch
      type(t_cpatch), pointer :: cp => NULL() ! pointer to a t_cpatch data-structure
   end type p_cpatch

   !---------------------------------------------------------------------------------------------------------
   ! data for one wheel:

   type :: t_wheel
      real(kind=8)              :: dx, dy, dz, droll, dyaw, dpitch, vx, vy, vz, vroll, vyaw, vpitch
      type(t_marker)            :: m_ws
      type(t_profile)           :: prw

      ! dx        [mm]   longitudinal displacement of wheel profile origin w.r.t. design layout
      ! dy        [mm]   lateral displacement of wheel profile origin w.r.t. design layout
      ! dz        [mm]   vertical displacement of wheel profile origin w.r.t. design layout
      ! droll     [rad]  rotation of wheel profile origin about x-axis w.r.t. design layout
      ! dyaw      [rad]  rotation of wheel profile origin about z-axis w.r.t. design layout
      ! dpitch    [rad]  rotation of wheel profile origin about y-axis w.r.t. design layout
      ! vx       [mm/s]  longitudinal velocity of wheel profile origin w.r.t. design layout
      ! vy       [mm/s]  lateral velocity of wheel profile origin w.r.t. design layout
      ! vz       [mm/s]  vertical velocity of wheel profile origin w.r.t. design layout
      ! vroll   [rad/s]  rotation velocity of wheel profile origin about x-axis w.r.t. design layout
      ! vyaw    [rad/s]  rotation velocity of wheel profile origin about z-axis w.r.t. design layout
      ! vpitch  [rad/s]  rotation velocity of wheel profile origin about y-axis w.r.t. design layout
      ! m_ws             position and orientation of wheel profile origin in wheelset coordinates
      !                  (excluding pitch rotation)
      ! prw              wheel profile data

   end type t_wheel

   !---------------------------------------------------------------------------------------------------------
   ! data for the wheel-set:

   type :: t_wheelset
      real(kind=8)  :: flback_dist, flback_pos, nom_radius
      real(kind=8)  :: s, x, y, z
      real(kind=8)  :: vs, vy, vz
      real(kind=8)  :: roll, yaw, pitch
      real(kind=8)  :: vroll, vyaw, vpitch
      real(kind=8)  :: fx_inp, fy_inp, fz_inp, my_inp
      real(kind=8)  :: z_cnt0, gap_min, delt_min, zw_min, a1, b1
      logical       :: has_overlap
      type(t_wheel), pointer :: whl => NULL()
      type(t_marker) :: m_trk

      ! flback_dist  [mm]   flange-back distance, e.g. 1360 mm
      ! flback_pos   [mm]   y-position of flange-back in wheel profile coordinates, e.g. -70 mm
      ! nom_radius   [mm]   wheel radius at the tape circle line, e.g. 460 mm
      ! s            [mm]   longitudinal position of wheel-set cm along the track curve - used on tracks
      ! x            [mm]   longitudinal position (shift) of wheel-set cm - used on roller rigs
      ! y            [mm]   lateral position (shift) of wheel-set cm
      ! z            [mm]   vertical position (shift) of wheel-set cm; downwards positive
      ! vs          [mm/s]  longitudinal translational velocity of wheel-set cm along track curve
      !              [mm]   when T=1, vs stores the position increment vs = shft_sws / dt* with dt* = 1
      ! vy          [mm/s]  lateral translational velocity of wheel-set cm 
      !              [mm]   when T=1, vy stores the position increment vy = shft_yws / dt* with dt* = 1
      ! vz          [mm/s]  vertical translational velocity of wheel-set cm
      !              [mm]   when T=1, vz stores the position increment vz = shft_zws / dt* with dt* = 1
      ! roll         [rad]  roll angle \phi of wheel-set
      ! yaw          [rad]  yaw angle \psi of wheel-set
      ! pitch        [rad]  pitch angle \theta of wheel-set
      ! vroll       [rad/s] roll angular velocity \dot{\phi} of wheel-set
      !              [rad]  when T=1, vroll stores the angle increment vroll = shft_rol / dt* with dt* = 1
      ! vyaw        [rad/s] yaw angular velocity \dot{\psi} of wheel-set
      !              [rad]  when T=1, vyaw stores the angle increment vyaw = shft_yaw / dt* with dt* = 1
      ! vpitch      [rad/s] pitch angular velocity \omega=\dot{\theta} of wheel-set
      !              [rad]  when T=1, vpitch stores the angle increment vpitch = shft_wpit / dt*, dt* = 1
      ! fx_inp        [N]   prescribed total force "on rail" in "wheelset" longitudinal direction (F=1)
      ! fy_inp        [N]   prescribed total force on rail in track lateral direction (n.y.a.)
      ! fz_inp        [N]   prescribed total force in track vertical direction (N>=1)
      ! my_inp      [N mm]  prescribed total moment on rail about wheelset lateral direction (F=2)
      ! z_cnt0       [mm]   z-position at which precise initial contact occurs
      ! gap_min      [mm]   overall minimum gap (max interpenetration) found, positive if no contact occurs
      ! delt_min    [rad]   contact angle at the location of the overall minimum gap value
      ! zw_min       [mm]   wheel surface height at the location of the overall minimum gap value
      ! a1, b1      [1/mm]  curvatures at the location of the overall minimum gap value
      ! has_overlap         flag indicating whether wheel and rail profiles overlap / do not overlap at all
      ! m_trk               position and orientation of wheel-set with respect to track origin,
      !                     excluding the pitch rotation
      ! whl                 data structure for wheel profile

   end type t_wheelset

   !---------------------------------------------------------------------------------------------------------
   ! data for one rail:

   type :: t_rail
      real(kind=8)               :: dy, dz, roll, vy, vz, vroll
      type(t_marker)             :: m_trk
      type(t_profile)            :: prr

      ! dy        [mm]   lateral displacement of rail w.r.t. design layout (left/right rail, non-mirrored)
      ! dz        [mm]   vertical displacement of rail w.r.t. design layout
      ! roll      [rad]  rotation of rail about x-axis w.r.t. design layout (left/right rail, non-mirrored)
      ! vy       [mm/s]  lateral velocity of rail w.r.t. design layout (track coord system, non-mirrored)
      ! vz       [mm/s]  vertical velocity of rail w.r.t. design layout
      ! vroll   [rad/s]  rotation velocity of rail about x-axis w.r.t. design layout (non-mirrored)
      ! m_trk            position and orientation of rail profile coordinates with respect to track system
      !                  (right rail, mirrored view)
      ! prr              rail profile data (right side)

   end type t_rail

   !---------------------------------------------------------------------------------------------------------
   ! data for the track or rollerset geometry:

   type :: t_trackdata
      integer               :: gauge_seqnum
      real(kind=8)          :: gauge_height, track_gauge, cant_angle
      real(kind=8)          :: rail_y0, rail_z0
      real(kind=8)          :: ky_rail, kz_rail, dy_defl, dz_defl, fy_rail, fz_rail
      real(kind=8)          :: nom_radius, vpitch_rol
      type(t_rail), pointer :: rai => NULL()

      ! gauge_height [mm]   gauge measuring height, e.g. 14mm, below track plane across top of rails.
      !                     The gauge computation is disabled when gauge_height <= 0.
      ! gauge_seqnum   -    sequence number of gauge face used: 0, 1 = inner-most face, 
      !                      2 = second inside face, particularly for situation with guard rail
      ! track_gauge  [mm]   lateral distance between inner sides of rails, e.g. 1435 mm
      ! cant_angle   [rad]  rail rotation towards track center, if not included in rail profile(s)
      !                     should not be used for roller rigs
      ! rail_y0      [mm]   lateral position of rail origin in track system, e.g. +/-755 mm
      ! rail_z0      [mm]   vertical position of rail origin in track system, e.g. +/-0.5 mm
      ! ky_rail     [N/mm]  track stiffness in lateral direction for massless rail deflection model (F=3)
      ! kz_rail     [N/mm]  track stiffness in vertical direction for massless rail deflection model (F=3)
      ! dy_defl      [mm]   lateral rail deflection (F=3, right rail, mirrored view)
      ! dz_defl      [mm]   vertical rail deflection (F=3, right rail, mirrored view)
      ! fy_rail      [N]    lateral spring force on rail at zero deflection for rail deflection model (F=3)
      !                     (left/right rail, non-mirrored)
      ! fz_rail      [N]    vertical spring force on rail at zero deflection for rail deflection model (F=3)
      ! nom_radius   [rad]  nominal radius of rollers, in case of roller rigs (config>=4)
      ! vpitch_rol  [rad/s] pitch angular velocity \omega=\dot{\theta} of rollers
      !              [rad]  when T=1, this stores the angle increment vpitch_rol = shft_rpit / dt*, dt* = 1
      ! rai                 data structure for rail

   end type t_trackdata

   !---------------------------------------------------------------------------------------------------------
   ! codes for computation of sensitivities of w/r problems

   integer, parameter :: nwsens_out = 3, nwsens_in  = 3
   integer, parameter :: iwout_fz  = 1, iwout_fx  = 2, iwout_fy  = 3
   integer, parameter :: iwin_zws  = 1, iwin_omgy = 2, iwin_vyws = 3

   ! nwsens_out    number of rows of sensitivities-matrix (1st index), output parameters,
   !               esp. forces and moments
   ! nwsens_in     number of columns of sensitivities-matrix, input parameters that can be perturbed

   !---------------------------------------------------------------------------------------------------------
   ! aggregate data for the half wheelset on track/roller combination:

   type :: t_ws_track
      type(t_metadata) :: meta    ! meta-data describing the calculation
      type(t_scaling)  :: scl     ! scaling factors for the CONTACT library
      type(t_ic)       :: ic      ! integer control digits
      type(t_material) :: mater   ! material-description of the bodies
      type(t_discret)  :: discr   ! parameters used to form potential contact areas and discretization
      type(t_subsurf)  :: subs    ! data for subsurface stresses
      type(t_friclaw)  :: fric    ! input parameters of friction law used
      type(t_kincns)   :: kin     ! kinematic description of the problem
      type(t_solvers)  :: solv    ! variables related to solution algorithms
      ! type(t_output) :: outpt   ! solution and derived quantities
      type(t_vec)      :: ftrk, ttrk, fws, tws, xavg, tavg
      real(kind=8)     :: dfz_dzws ! sensitivity of fz_tr to z_ws

      type(t_wheelset) :: ws      ! half wheel-set data, including wheel profile
      type(t_trackdata):: trk     ! half track/roller data, including rail profile
      integer          :: numcps  ! number of contact problems at current time
      integer          :: numtot  ! total number of contact problems held in allcps
      type(p_cpatch)   :: allcps(MAX_NUM_CPS) ! pointers to the data for all contact problems of ws on trk
                                  ! first numcps == current time, remainder == unconnected of prev.time
   end type t_ws_track

   type :: p_ws_track
      type(t_ws_track), pointer :: wtd => NULL() ! pointer to a t_ws_track data-structure
   end type p_ws_track

contains

!------------------------------------------------------------------------------------------------------------

   function use_brute(ic, prw)
!--function: determine from an ic-struct & profiles used whether to use the (brute force) grid-based
!            contact search
      implicit none
!--result value
      logical                     :: use_brute
!--subroutine arguments
      type(t_ic),      intent(in) :: ic
      type(t_profile), intent(in) :: prw

      use_brute = (ic%discns1_eff.eq.5 .or. ic%discns1_eff.eq.6 .or. ic%discns1_eff.eq.9)
      use_brute = use_brute .or. prw%is_varprof()

   end function use_brute

!------------------------------------------------------------------------------------------------------------

   subroutine wheel_ini(whl)
!--purpose: Initialize a wheel data-structure
      implicit none
!--subroutine parameters:
      type(t_wheel)  :: whl

      call marker_init( whl%m_ws )

      whl%dx         = 0d0
      whl%dy         = 0d0
      whl%dz         = 0d0
      whl%droll      = 0d0
      whl%dyaw       = 0d0
      whl%dpitch     = 0d0
      whl%vx         = 0d0
      whl%vy         = 0d0
      whl%vz         = 0d0
      whl%vroll      = 0d0
      whl%vyaw       = 0d0
      whl%vpitch     = 0d0

      call profile_ini(whl%prw)

   end subroutine wheel_ini

!------------------------------------------------------------------------------------------------------------

   subroutine rail_ini(rail)
!--purpose: Initialize a rail data-structure
      implicit none
!--subroutine parameters:
      type(t_rail)  :: rail

      call marker_init( rail%m_trk )

      rail%dy         = 0d0
      rail%dz         = 0d0
      rail%roll       = 0d0
      rail%vy         = 0d0
      rail%vz         = 0d0
      rail%vroll      = 0d0

      call profile_ini(rail%prr)

   end subroutine rail_ini

!------------------------------------------------------------------------------------------------------------

   subroutine wheelset_ini(ws)
!--purpose: Initialize a wheelset data-structure
      implicit none
!--subroutine parameters:
      type(t_wheelset) :: ws

      ! geometry parameters

      ws%flback_dist   = 1360d0
      ws%flback_pos    =  -70d0
      ws%nom_radius    =  460d0

      ! position & orientation
      ws%s             =    0d0
      ws%y             =    0d0
      ws%z             =    0d0
      ws%roll          =    0d0
      ws%yaw           =    0d0
      ws%pitch         =    0d0

      ! velocity parameters
      ws%vs            =   30d3
      ws%vy            =    0d0
      ws%vz            =    0d0
      ws%vroll         =    0d0
      ws%vyaw          =    0d0
      ws%vpitch        =    0d0

      ! prescribed total forces and moments
      ws%fx_inp        =    0d0
      ws%fy_inp        =    0d0
      ws%fz_inp        =    0d0
      ws%my_inp        =    0d0

      ! partial results of contact geometry analysis
      ws%has_overlap   = .false.
      ws%gap_min       =    0d0
      ws%a1            =    0d0
      ws%b1            =    0d0
 
      ! initialize the wheelset position marker

      call marker_init( ws%m_trk )

      ! initialize data for left and right wheels

      allocate(ws%whl)
      call wheel_ini(ws%whl)

   end subroutine wheelset_ini

!------------------------------------------------------------------------------------------------------------

   subroutine wrprof_ini(wtd)
!--purpose: Initialize a wheelset/track data-structure
      implicit none
!--subroutine parameters:
      type(t_ws_track)       :: wtd
!--local variables:
      integer                :: icp

      ! initialize contact patch data-structures

      wtd%numtot = 0
      wtd%numcps = 0
      do icp = 1, MAX_NUM_CPS
         wtd%allcps(icp)%cp => NULL()
      enddo

      call meta_init( wtd%meta, 2 )

      ! initialize control-digits that are not provided in the input-file
      ! see additional initializations for C.library in m_caddon_data: cntc_activate_wtd

      call ic_init( wtd%ic )

      ! initialize material data

      call mater_init( wtd%mater, wtd%ic )

      ! initialize wheel-set geometry, position

      call wheelset_ini(wtd%ws)

      ! initialize track geometry

      wtd%trk%gauge_height =   14d0
      wtd%trk%gauge_seqnum =    0
      wtd%trk%track_gauge  = 1435d0
      wtd%trk%cant_angle   =    0d0
      wtd%trk%rail_y0      =    0d0
      wtd%trk%rail_z0      =    0d0
      wtd%trk%ky_rail      =    1d3
      wtd%trk%kz_rail      =    1d3
      wtd%trk%dy_defl      =    0d0
      wtd%trk%dz_defl      =    0d0
      wtd%trk%fy_rail      =    0d0
      wtd%trk%fz_rail      =    0d0

      ! initialize data for left and right rails

      allocate(wtd%trk%rai)
      call rail_ini(wtd%trk%rai)

      ! initialize friction data

      call fric_init(wtd%fric)

      ! initialize pot.con/grid discretization data

      wtd%discr%npot_max  = 20000
      wtd%discr%dx        = 0.5d0
      wtd%discr%ds        = 0.5d0
      wtd%discr%dqrel     = 1d0
      wtd%discr%angl_sep  = pi
      wtd%discr%dist_turn = 12d0
      wtd%discr%dist_sep  = 8d0
      wtd%discr%dist_comb = 4d0

      ! initialize kinematic parameter settings

      call kincns_init( wtd%kin, wtd%ic, wtd%fric, wtd%discr%dx )

      ! initialize solver settings

      call solv_init( wtd%solv )

      ! initialize subsurface stress data

      wtd%subs%nblock =    0

      ! initialize other data in the wtd structure

   end subroutine wrprof_ini

!------------------------------------------------------------------------------------------------------------

   subroutine cp_init(cp)
!--purpose: Initialize a contact-patch data-structure
      implicit none
!--subroutine parameters:
      type(t_cpatch)    :: cp

      cp%micp = vec_zero()
      call marker_init(cp%mref)
      call marker_init(cp%mpot)
      cp%nsub = 0
      cp%y_sep => NULL()
      cp%f_sep => NULL()
      call grid_nullify(cp%rail_srfc)
      call grid_nullify(cp%whl_srfc)
      call grid_nullify(cp%curv_ref)
      call gf3_nullify(cp%curv_nrm)
      call gf3_nullify(cp%curv_incln)
      cp%gd => NULL()
      cp%prev_icp(1:MAX_NUM_CPS) = 0

!     real(kind=8)     :: gap_min, totgap, wgt_xgap, wgt_ygap, wgt_zgap, wgt_agap
!     real(kind=8)     :: delttr ...

   end subroutine cp_init

!------------------------------------------------------------------------------------------------------------

   subroutine tpatch_destroy(cp)
!--purpose: Cleanup a contact-patch data-structure
      implicit none
!--subroutine parameters:
      type(t_cpatch)  :: cp
!--local variables:

      call destroy_arr( cp%y_sep )
      call destroy_arr( cp%f_sep )
      call grid_destroy( cp%rail_srfc )
      call grid_destroy( cp%whl_srfc )
      call grid_destroy( cp%curv_ref )
      call gf3_destroy( cp%curv_nrm )
      call gf3_destroy( cp%curv_incln )

      if (associated( cp%gd )) then
         call gd_destroy( cp%gd )
         deallocate( cp%gd )
         nullify( cp%gd )
      endif

   end subroutine tpatch_destroy

!------------------------------------------------------------------------------------------------------------

   subroutine ppatch_destroy(pp)
!--purpose: Cleanup a pointer to contact-patch data-structure
      implicit none
!--subroutine parameters:
      type(p_cpatch)  :: pp

      if (associated(pp%cp)) then
         call tpatch_destroy(pp%cp)
         deallocate(pp%cp)
         pp%cp => NULL()
      endif

   end subroutine ppatch_destroy

!------------------------------------------------------------------------------------------------------------

   subroutine wrprof_destroy(wtd)
!--purpose: Cleanup a wheelset/track data-structure
      implicit none
!--subroutine parameters:
      type(t_ws_track)  :: wtd
!--local variables:
      integer                :: icp, iblk

      ! no action: meta, scl, ic, mater, discr, kin, solv

      ! cleanup subsurf data, friction data

      do iblk = 1, wtd%subs%nblock
         call subsblk_destroy( wtd%subs%blocks(iblk) )
      enddo
      call fric_destroy( wtd%fric )

      ! wheelset, trackdata: just clean up the profiles

      call profile_destroy( wtd%ws%whl%prw )
      deallocate( wtd%ws%whl )
      nullify( wtd%ws%whl )

      call profile_destroy( wtd%trk%rai%prr )
      deallocate( wtd%trk%rai )
      nullify( wtd%trk%rai )

      ! clean up contact patch data (grids, gd)

      do icp = 1, MAX_NUM_CPS
         if (associated(wtd%allcps(icp)%cp)) then
         call cp_destroy( wtd%allcps(icp)%cp )
         deallocate( wtd%allcps(icp)%cp )
         nullify( wtd%allcps(icp)%cp )
         endif
      enddo

   end subroutine wrprof_destroy

!------------------------------------------------------------------------------------------------------------

end module m_wrprof_data
