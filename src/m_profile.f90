!------------------------------------------------------------------------------------------------------------
! m_profile - data-structures for wheel/rail profiles
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_profile

use m_hierarch_data

implicit none
private

   public  t_profile
   public  profile_is_varprof
   public  profile_ini
   public  profile_copy
   public  profile_destroy
   public  varprof_destroy_slices

   public  profile_setopt
   public  profile_read_config
   public  profile_write_config

   public  profile_store_values
   public  profile_finish_grid

   public  read_simpack_profile
   public  read_miniprof_profile

   private delete_outside_boundbox
   private filter_close_points
   private remove_vertical_slopes
   private filter_double_kinks
   private check_z_pos_down
   private check_inversion
   private invert_points
   private filter_receding_order
   public  profile_find_kinks
   private segs_have_intersection

   ! codes for different file formats

   integer,      parameter :: FTYPE_SLICES   = -1
   integer,      parameter :: FTYPE_PLAIN2D  =  0
   integer,      parameter :: FTYPE_SIMPACK  =  1
   integer,      parameter :: FTYPE_MINIPROF =  2
   public FTYPE_SLICES, FTYPE_PLAIN2D, FTYPE_SIMPACK, FTYPE_MINIPROF

   ! codes for different spline variants

   integer,      parameter :: SPL2D_INTPOL   = 1
   integer,      parameter :: SPL2D_APPROX   = 2
   public SPL2D_APPROX, SPL2D_INTPOL

   ! maximum #slices in a single variable profile (avoids growing/shrinking of arrays of pointers)

   integer,      parameter :: MAX_NUM_SLCS   =  4999
   public MAX_NUM_SLCS

   !---------------------------------------------------------------------------------------------------------

   ! data with respect to a single wheel or rail profile

   type :: t_profile
      ! basic profile data + profile configuration
      character(len=256) :: fname
      type(t_grid)       :: grd_data
      type(t_gridfnc3)   :: ext_data
      logical            :: has_kyield

      integer            :: err_hnd, ierror, mirror_y, mirror_z, ismooth
      real(kind=8)       :: f_max_omit, sclfac, smth, zig_thrs, kink_high, kink_low, kink_wid

      ! variable profiles: #slices, data per slice
      integer                                     :: nslc, nfeat, nkink, naccel, u_method
      real(kind=8),       dimension(:),   pointer :: slc_u       => NULL()
      character(len=256), dimension(:),   pointer :: slc_nam     => NULL()
      type(p_grid),       dimension(:)            :: slc_grd(MAX_NUM_SLCS)
      real(kind=8),       dimension(:,:), pointer :: slc_s_f     => NULL()
      integer,            dimension(:,:), pointer :: slc_ib      => NULL()
      integer,            dimension(:),   pointer :: iseg_p      => NULL()
      integer,            dimension(:),   pointer :: kink_if     => NULL()
      integer,            dimension(:),   pointer :: accel_if    => NULL()
      type(t_bspline2d)                           :: spl2d
    contains
       procedure :: is_varprof  => profile_is_varprof

      ! fname            file-name for profile file as given by user
      ! grd_data         profile data for the rail or wheel, w.r.t. profile datum (tape circle)
      !  - y      [mm]   lateral coordinates of profile points, relative to "profile reference"
      !                  note: y-coordinates are ascending for rails, descending for wheels.
      !  - z      [mm]   vertical coordinates of profile points, relative to "profile reference"
      !  - s_prf  [mm]   arc length coordinates of profile points, relative to "profile reference"
      !  - spl    [ ]    structure with smoothing spline parameters for {y(s), z(s)}

      ! has_kyield       flag indicating whether kyield is given in ext_data
      ! ext_data         additional data at profile points: vn==kyield; n.y.a. for variable profiles

      ! err_hnd          flag for error handling: -2: continue as much as possible; -1: suppress warnings;
      !                  0: warn and continue (default); 1: signal errors and abort
      ! f_max_omit       fraction, signal error if more than maxomit of profile points are discarded
      !                  after cleanup of profile. Default 0.5. Set to 1 to disable this error.
      ! ierror           flag <0 = error codes, 0 = profile interpreted ok, 1 = not set

      ! mirror_y         flag -1 or 0=file concerns right wheel; 1=file concerns left, needs mirroring
      ! mirror_z         flag: 0=autodetect whether z needs mirroring, -1=no mirroring, 1=mirror z-values
      ! sclfac  [mm/len] profile scaling factor needed to convert input data to mm

      ! ismooth          selection of smoothing method. 
      !                    0 = original, non-weighted PP smoothing spline,
      !                        smoothing penalty on 2nd derivative, with kinks / no accelerations
      !                    1 = weighted PP smoothing spline with equal #segments as input profile
      !                        smoothing penalty on 2nd derivative, with kinks / no accelerations
      !                    2 = weighted smoothing B-spline with fewer #segments than input profile
      !                        smoothing penalty on 3rd derivative, with kinks and accelerations
      ! smth             parameter lambda = (1-p)/p for non-weighted smoothing spline (ismooth=0),
      !                  with weight p for the data and weight (1-p) for the 2nd derivative, or 
      !                  l_filt for weighted splines (ismooth=1,2)

      ! zig_thrs  [rad]  angle threshold for zig-zag pattern detection. Default pi/2. 
      !                  Set to pi to disable zig-zag detection.
      ! kink_high [rad]  angle threshold for kink detection. Default pi/6. 
      !                  Set to pi to disable kink detection.
      ! kink_low  [rad]  angle threshold for neighbouring points in kink detection. Default kink_high/5
      ! kink_wid  [mm]   half-width of window used for kink detection, default 2 mm

      ! extension for variable profiles:

      ! nslc             in case of slcs-files: number of slices (>=1);
      !                  -1 for regular profiles (prr, ban, txt)
      ! slc_u    [mm]    u-positions of slices, s_fc along the track curve or th_w on wheel, after
      !         or [rad] shifting+scaling
      ! slc_nam          file-names of profile slices as given in a slcs-file
      ! slc_grd          grid-data for each of the profile slices
      ! u_method         interpolation in longitudinal direction: SPL2D_APPROX, SPL2D_INTPOL
      !
      ! nfeat            number of feature positions in lateral direction (>= 2)
      ! nkink            number of kinks at feature positions
      ! naccel           number of accelerations at feature positions
      ! slc_s_f          nfeat lateral s_f-positions per slice at geometrical features, used as
      !                  breaks between different parts of the slice. Using s_f = 0 at start of profile.
      ! slc_ib           indices ib0, ib1 for the first and last breaks/features present in each slice
      ! iseg_p           index of first segment per part after resampling. nseg == iseg_p(npart+1)-1
      ! kink_if          feature numbers f where a kink occurs between parts f-1 and f
      ! accel_if         feature numbers f where an acceleration occurs between parts f-1 and f
      !
      ! spl2d            2D tensor B-spline representation

   end type t_profile

contains

!------------------------------------------------------------------------------------------------------------

   function profile_is_varprof(this)
!--function: tell if profile has variation in running direction
      implicit none
!--result value
      logical                      :: profile_is_varprof
!--subroutine arguments
      class(t_profile), intent(in) :: this

      profile_is_varprof = (this%nslc.gt.0)
   end function profile_is_varprof

!------------------------------------------------------------------------------------------------------------

   subroutine profile_ini(prf)
!--purpose: initialize a profile data-structure
      implicit none
!--subroutine parameters:
      type(t_profile) :: prf
!--local variables:
      real(kind=8), parameter :: pi     = 4d0*atan(1d0)

      ! clean up previous contents, deallocate memory
      call profile_destroy( prf )

      prf%fname      = ' '

      prf%has_kyield = .false.

      prf%err_hnd    = 0        ! warn and continue
      prf%ierror     = 1        ! error: not set
      prf%f_max_omit = 0.5d0    ! max fraction of profile points deleted

      prf%mirror_y   = 0        ! 0: no mirroring in y-direction
      prf%mirror_z   = 0        ! 0: autodetect mirroring in z-direction

      prf%sclfac     = 1d0      ! 1: no scaling
      prf%ismooth    = 0        ! 0: original spline method
      prf%smth       = 0d0      ! no smoothing

      prf%zig_thrs   = 150d0 * pi/180d0  ! zig-zag detection: >= +/- 150 deg

      prf%kink_high  =  30d0 * pi/180d0  ! kink detection: >= 30 deg at i,
      prf%kink_low   = prf%kink_high/5d0 !                 <=  6 deg at k=3 neighbouring points
      prf%kink_wid   = 2d0      ! mm

   end subroutine profile_ini

!------------------------------------------------------------------------------------------------------------

   subroutine profile_copy(p_in, p_out)
!--purpose: copy profile data-structure, excluding slices administration
      implicit none
!--subroutine parameters:
      type(t_profile), intent(in)  :: p_in
      type(t_profile), intent(out) :: p_out
!--local variables:

      p_out%fname       = p_in%fname

      call grid_copy(p_in%grd_data, p_out%grd_data, with_spline=.true.)

      p_out%has_kyield  = p_in%has_kyield
      if (p_in%has_kyield) then
         call gf3_new(p_out%ext_data, 'ext_data', p_out%grd_data )
         call gf3_copy(AllElm, p_in%ext_data, p_out%ext_data, ikALL )
      endif

      p_out%err_hnd     = p_in%err_hnd
      p_out%ierror      = p_in%ierror
      p_out%f_max_omit  = p_in%f_max_omit

      p_out%mirror_y    = p_in%mirror_y
      p_out%mirror_z    = p_in%mirror_z

      p_out%sclfac      = p_in%sclfac
      p_out%ismooth     = p_in%ismooth
      p_out%smth        = p_in%smth

      p_out%zig_thrs    = p_in%zig_thrs

      p_out%kink_high   = p_in%kink_high
      p_out%kink_low    = p_in%kink_low
      p_out%kink_wid    = p_in%kink_wid

      ! reset slices administration

      call varprof_destroy_slices(p_out)

   end subroutine profile_copy

!------------------------------------------------------------------------------------------------------------

   subroutine profile_destroy(prf)
!--purpose: cleanup a profile data-structure
      implicit none
!--subroutine parameters:
      type(t_profile) :: prf
!--local variables:

      call grid_destroy( prf%grd_data )
      call gf3_destroy(  prf%ext_data )
      call varprof_destroy_slices( prf )

   end subroutine profile_destroy

!------------------------------------------------------------------------------------------------------------

   subroutine varprof_destroy_slices(vprf)
!--purpose: reset the administration of slices in a variable profile
      implicit none
!--subroutine parameters:
      type(t_profile) :: vprf
!--local variables:
      integer         :: islc

      vprf%nslc     = 0
      vprf%nfeat    = 0
      vprf%nkink    = 0
      vprf%naccel   = 0
      vprf%u_method = SPL2D_INTPOL

      ! destroy arrays

      call destroy_arr(vprf%slc_u )
      call destroy_arr(vprf%slc_nam )
      call destroy_arr(vprf%slc_s_f )
      call destroy_arr(vprf%slc_ib )
      call destroy_arr(vprf%iseg_p )
      call destroy_arr(vprf%kink_if )
      call destroy_arr(vprf%accel_if )

      ! destroy slices

      do islc = 1, MAX_NUM_SLCS
         if (associated(vprf%slc_grd(islc)%g)) then
            call grid_destroy(vprf%slc_grd(islc)%g)
            deallocate(vprf%slc_grd(islc)%g )
            nullify(vprf%slc_grd(islc)%g )
         endif
      enddo

      ! destroy spline

      call bspline2d_destroy(vprf%spl2d )

   end subroutine varprof_destroy_slices

!------------------------------------------------------------------------------------------------------------

   subroutine profile_setopt(prf, nints, iparam, nreals, rparam, scl_len)
!--purpose: set options in profile data-structure
!  iparam         - integer configuration parameters
!                     1: itype     0 = rail, 1 = wheel profile, -1 = taken from file extension (default)
!                     2:  -        not used
!                     3: mirrory   0 or -1 = no mirroring (default), 1 = mirror y coordinate values
!                     4: mirrorz   0 = autodetect (default), -1 = no mirroring, 1 = mirror z values
!                     5: errhndl   configuration of error handling. 
!                                   -2 = continue as much as possible, suppress error messages; 
!                                   -1 = suppress warnings; 0: warn and continue (default);
!                                    1 = signal errors and abort
!                     6: ismooth   selection of smoothing method. 0 = original smoothing spline (default),
!                                    1 = weighted PP smoothing spline, 2 = weighted smoothing B-spline (best)
!  rparam         - real configuration parameters
!                     1: sclfac    scaling factor for conversion to [mm], e.g. 1e3 for data given in [m]
!                                  default (sclfac<=0): using the active unit convention
!                     2: smooth    smoothing parameter lambda for non-weighted spline or l_filt for
!                                  weighted spline smoothing
!                     3: maxomit   fraction: signal error if more than maxomit of profile points are
!                                  discarded after cleanup of profile. Default 0.5, use 1 to disable check.
!                     4: zigthrs   angle threshold for zig-zag detection. Default 5/6*pi, >=pi to disable.
!                     5: kinkhigh  angle threshold for kink detection. Default pi/6, >=pi to disable.
!                     6: kinklow   angle threshold for neighbouring points in kink detection. 
!                                  default kinkhigh/5.
!                     7: kinkwid   half-width of window used for kink detection, [len], default 2 mm
      implicit none
!--subroutine parameters:
      type(t_profile)          :: prf
      integer,      intent(in) :: nints, nreals
      integer,      intent(in) :: iparam(nints)
      real(kind=8), intent(in) :: rparam(nreals), scl_len
!--local variables:
      real(kind=8), parameter :: pi     = 4d0*atan(1d0)

      ! unpack option values, fill in defaults

      if (nints.ge.3) then
         prf%mirror_y = iparam(3) ! 0 or -1=no, 1=yes
      else
         prf%mirror_y = 0         ! no mirroring
      endif
      if (nints.ge.4) then
         prf%mirror_z = iparam(4) ! 0=autodetect, -1=no, 1=yes
      else
         prf%mirror_z = 0         ! autodetect
      endif
      if (nints.ge.5) then
         prf%err_hnd  = iparam(5) ! -2=continue; -1=suppress warn; 0=warn/continue; 1=abort
      else
         prf%err_hnd  = 0
      endif
      if (nints.ge.6) then
         prf%ismooth = iparam(6)  ! 0=original spline; 1=weighted PP-spline; 2=weighted B-spline
      else
         prf%ismooth = 0
      endif
   
      if (nreals.ge.1) then
         prf%sclfac  = rparam(1)  ! unit conversion, e.g. sclfac = 1000 [mm] / [user units]
         if (prf%sclfac.le.5e-4) prf%sclfac = scl_len  ! error handling
      else
         prf%sclfac  = scl_len    ! default: cf. unit convention
      endif
      if (nreals.ge.2) then
         prf%smth  = rparam(2)    ! non-weighted spline: lambda = (1-p)/p for profile smoothing
                                  ! weighted spline: wavelength l_filt
         prf%smth  = max(0d0, prf%smth)
      else
         prf%smth  = 0d0          ! no smoothing
      endif
      if (nreals.ge.3) then
         prf%f_max_omit = rparam(3) ! max. fraction of points deleted
      else
         prf%f_max_omit = 0.5d0
      endif
      if (nreals.ge.4) then
         prf%zig_thrs  = rparam(4) ! angle threshold for zig-zag detection
         prf%zig_thrs  = max(pi/4d0, prf%zig_thrs)
      else
         prf%zig_thrs  = 150d0 * pi/180d0 ! 150deg
      endif
      if (nreals.ge.5) then
         prf%kink_high = rparam(5) ! angle threshold for kink detection
      else
         prf%kink_high =  30d0 * pi/180d0  ! 30deg
      endif
      if (nreals.ge.6) then
         prf%kink_low  = rparam(6) ! angle threshold for neighbouring points in kink detection
      else
         prf%kink_low  = prf%kink_high / 5d0
      endif
      if (nreals.ge.7) then
         prf%kink_wid = rparam(7) * scl_len ! half-width window used for kink detection, default 2 mm
      else
         prf%kink_wid = 2d0
      endif

!  if (idebug.ge.2) then
!     write(bufout,'(a,a30,a,i3,3(a,i3),2(a,g10.1))') pfx,subnam,'(',ire,'): itype=',itype,             &
!               ', mirrory=',mirrory,', mirrorz=',mirrorz,', scl=',sclfac,', smooth=',smooth
!     call write_log(1, bufout)
!  endif

   end subroutine profile_setopt

!------------------------------------------------------------------------------------------------------------

   subroutine profile_read_config(prf, prfnam, linp, ncase, linenr, idebug, ieof, lstop, zerror)
!--purpose: read profile filename and configuration parameters from inp-file
      implicit none
!--subroutine parameters:
      type(t_profile)              :: prf
      integer,      intent(in)     :: linp, ncase, idebug
      logical,      intent(in)     :: lstop
      integer,      intent(inout)  :: linenr, ieof
      logical,      intent(inout)  :: zerror
      character(len=*), intent(in) :: prfnam
!--local variables:
      integer,      parameter :: mxnval = 20
      integer          :: ints(mxnval), iparam(10), nval, prfopt, sub_ierror
      logical          :: flags(mxnval)
      real(kind=8)     :: dbles(mxnval), rparam(10)
      character*256    :: strngs(mxnval), descrp
! TODO: pass sub_ierror from readline up to calling routine

      ! using profile_setopt(iparam, rparam) to store parameters and fill in default values

      iparam(1:10) = 0
      rparam(1:10) = 0d0
      associate( itype      => iparam(1),                          mirror_y   => iparam(3),             &
                 mirror_z   => iparam(4), err_hnd    => iparam(5), ismooth    => iparam(6),             &
                 sclfac     => rparam(1), smooth     => rparam(2), f_max_omit => rparam(3),             &
                 zig_thrs   => rparam(4), kink_high  => rparam(5), kink_low   => rparam(6),             &
                 kink_wid   => rparam(7) )                                                               
      itype    = 1

      ! get profile filename and basic configuration options

      descrp = 'file with ' //  trim(prfnam) // ' profile'
      call readline(linp, ncase, linenr, descrp, 'siddII', ints, dbles, flags, strngs, mxnval, nval,    &
                idebug, ieof, lstop, sub_ierror)

      ! retrieve mirror_y, mirror_z, scale and smooth

      mirror_y = ints(1)
      mirror_z = 0      ! reset to default if no value is given
      if (nval.ge.5) mirror_z = ints(2)
      prfopt   = 0
      if (nval.ge.6) prfopt   = ints(3)
      sclfac   = max(1d-6, dbles(1))
      smooth   = max(0d0,  dbles(2))

      zerror = zerror .or. .not.check_irng ('MIRROR_Y', mirror_y, -1, 1)
      zerror = zerror .or. .not.check_irng ('MIRROR_Z', mirror_z, -1, 1)
      zerror = zerror .or. .not.check_irng ('PRFOPT',   prfopt,    0, 2)

      prf%fname = trim(strngs(1))

      ! get extra configuration options, store using profile_setopt

      if (prfopt.eq.0) then

         call profile_setopt(prf, 4, iparam, 2, rparam, 1d0)

      elseif (prfopt.ge.1) then

         descrp = 'extra config for ' //  trim(prfnam) // ' profile (1)'
         call readline(linp, ncase, linenr, descrp, 'iaaad', ints, dbles, flags, strngs, mxnval,        &
                        nval, idebug, ieof, lstop, sub_ierror)

         ismooth    = max(  0, min(  2, ints(1)))
         zig_thrs   = dbles(1)
         kink_high  = dbles(2)
         kink_low   = dbles(3)
         kink_wid   = dbles(4)

         if (prfopt.ge.2) then
            descrp = 'extra config for ' //  trim(prfnam) // ' profile(2)'
            call readline(linp, ncase, linenr, descrp, 'id', ints, dbles, flags, strngs, mxnval,        &
                        nval, idebug, ieof, lstop, sub_ierror)

            err_hnd    = max( -2, min(  1, ints(1)))
            f_max_omit = max(0d0, min(1d0, dbles(1)))
         else
            err_hnd    = prf%err_hnd
            f_max_omit = prf%f_max_omit
         endif

         call profile_setopt(prf, 6, iparam, 7, rparam, 1d0)
      endif

      end associate
   end subroutine profile_read_config

!------------------------------------------------------------------------------------------------------------

   subroutine profile_write_config(prf, is_wheel)
!--purpose: write profile filename and configuration to inp-file
      implicit none
!--subroutine arguments:
      type(t_profile)            :: prf
      integer,        intent(in) :: is_wheel
!--local variables:
      integer             :: len1, prfopt
      character(len=1)    :: namprf
      character(len=6)    :: namsmth
      character(len=80)   :: spaces

      if (is_wheel.ge.1) then
         namprf = 'W'
      else
         namprf = 'R'
      endif
      if (prf%ismooth.eq.0) then
         namsmth = 'LAMBDA'
      else
         namsmth = 'L_FILT'
      endif

      spaces = ' '
      len1 = len(trim(prf%fname))

      if (prf%ismooth.eq.0) then

         ! original smoothing method - dont write detailed configuration

         if     (len1.le.27) then
            write(linp, 7100) '''', trim(prf%fname), '''', prf%mirror_y, prf%sclfac, prf%smth,          &
                   prf%mirror_z, spaces(1:28-len1), namprf, namsmth
         else
            write(linp, 7100) '''', trim(prf%fname), '''', prf%mirror_y, prf%sclfac, prf%smth,          &
                   prf%mirror_z, '  ', namprf, namsmth
         endif

      else

         ! new smoothing method - write additional configuration used
         ! ugly: repeat defaults err_hnd=0, max_omit=0.5d0

         prfopt = 1
         if (prf%err_hnd.ne.0 .or. prf%f_max_omit.ne.0.5d0) prfopt = 2

         if     (len1.le.27) then
            write(linp, 7200) '''', trim(prf%fname), '''', prf%mirror_y, prf%sclfac, prf%smth,          &
                   prf%mirror_z, prfopt, spaces(1:28-len1), namprf, namsmth
         else
            write(linp, 7200) '''', trim(prf%fname), '''', prf%mirror_y, prf%sclfac, prf%smth,          &
                   prf%mirror_z, prfopt, '  ', namprf, namsmth
         endif
         write(linp, 7300) prf%ismooth, prf%zig_thrs*180d0/pi, prf%kink_high*180d0/pi,                  &
                   prf%kink_low*180d0/pi, prf%kink_wid
         if (prfopt.ge.2) write(linp, 7400) prf%err_hnd, prf%f_max_omit

      endif

 7100 format(1x, 3a, i4, f9.3, g10.3,  i3, a, a1, 'FNAME, MIRRORY, SCALE, ', a, ', MIRRORZ')
 7200 format(1x, 3a, i3, f9.3,  f9.4, 2i3, a, a1, 'FNAME, MIRRORY, SCALE, ', a, ', MIRRORZ, PRFOPT')
 7300 format(1x, i7, 3(f8.1,'d'), f9.4, 14x, 'ISMOOTH, ZIGTHRS, KINKHIGH, KINKLOW, KINKWID')
 7400 format(1x, i7,   f9.4,       27x, 14x, 'ERRHND, MAXOMIT')

   end subroutine profile_write_config

!------------------------------------------------------------------------------------------------------------

   subroutine profile_store_values(prf, is_wheel, npoint_arg, rvalues, idebug)
!--purpose: process and store array of w/r profile values (one slice, no support for var.profiles) 
      implicit none
!--subroutine parameters:
      type(t_profile)                         :: prf
      integer,                    intent(in)  :: is_wheel, npoint_arg, idebug
      real(kind=8), dimension(:), intent(in)  :: rvalues(2*npoint_arg)  ! 2 x npoint_arg
!--local variables:
      integer,      parameter  :: ncolpnt = 2
      real(kind=8), parameter  :: point_dist_min = 1d-4
      integer                  :: npoint, ipnt, icol_kyield, my_mirror_z, flip_data
      logical                  :: fill_spline
      real(kind=8), dimension(:,:), pointer :: points => NULL()  ! npoints x 2: [y, z]

      associate(fname      => prf%fname,      prf_grd  => prf%grd_data,  mirror_y   => prf%mirror_y,   &
                mirror_z   => prf%mirror_z,   sclfac   => prf%sclfac,    ierror     => prf%ierror)

      call timer_start(itimer_profil)

      ierror = 0
      fname  = '<values>'

      if (idebug.ge.2) then
         write(bufout,*) 'profile_store_values: is_wheel=',is_wheel,', npoint=',npoint_arg
         call write_log(1, bufout)
      endif

      if (is_wheel.ne.0 .and. is_wheel.ne.1) then
         ierror = 191
         write(bufout,'(a,i4,a)') ' ERROR 191: incorrect value for itype (', is_wheel,', must be 0 or 1.'
         call write_log(1, bufout)
      endif

      ! copy data into the format used by profile_finish_grid

      npoint = npoint_arg
      call reallocate_arr(points, npoint, ncolpnt)

      do ipnt = 1, npoint
         points(ipnt,1) = rvalues(1+2*(ipnt-1))
         points(ipnt,2) = rvalues(2+2*(ipnt-1))

         if (idebug.ge.3) then
            write(bufout,*) ' i=',ipnt,': yi,zi=',points(ipnt,1), points(ipnt,2)
            call write_log(1, bufout)
         endif
      enddo

      !  - modification 1: point.dist.min - remove points too close together

      call filter_close_points(points, npoint, npoint, ncolpnt, point_dist_min, is_wheel, idebug)

      !  - modification 7: mirror y if requested by subroutine argument

      if (mirror_y.ge.1) then
         do ipnt = 1, npoint
            points(ipnt,1) = -points(ipnt,1)
         enddo
      endif

      !  - modification 8: mirror z if necessary to get z positive downwards

      my_mirror_z = mirror_z
      if (my_mirror_z.eq.0) then
         call check_z_pos_down(fname, npoint, points, is_wheel, idebug, my_mirror_z)
      endif

      if (my_mirror_z.eq.1) then
         do ipnt = 1, npoint
            points(ipnt,2) = -points(ipnt,2)
         enddo
      endif

      !  - modification 9: inversion

      if (ierror.eq.0) then
         call check_inversion(fname, npoint, points, 1, 2, is_wheel, idebug, 0, flip_data, ierror)
      endif

      if (ierror.eq.0 .and. flip_data.eq.1) then
         call invert_points(npoint, points, is_wheel, idebug)
      endif

      !  - first/last points should now be in appropriate order

      if (ierror.eq.0) then
         call check_overall_order(npoint, points, is_wheel, idebug, ierror)
      endif

      if (ierror.eq.0) then

         !  - modification 10: length scaling factor

         do ipnt = 1, npoint
            points(ipnt,1) = sclfac * points(ipnt,1)
            points(ipnt,2) = sclfac * points(ipnt,2)
         enddo

         ! store points in profile-grid, modify, and create spline representation

         icol_kyield = 0
         fill_spline = .true.
         call profile_finish_grid(fname, prf_grd, prf%ext_data, prf, npoint, points, is_wheel,          &
                        icol_kyield, fill_spline, idebug)

      endif

      call timer_stop(itimer_profil)

      end associate
   end subroutine profile_store_values

!------------------------------------------------------------------------------------------------------------

   subroutine profile_finish_grid(fname, prf_grd, prf_ext, prf_opt, npoint, points, is_wheel,           &
                icol_kyield, fill_spline, idebug)
!--purpose: copy points for one slice to a grid, make further modifications, and create smoothing spline
      implicit none
!--subroutine parameters:
      character*(*), intent(in)    :: fname
      type(t_grid)                 :: prf_grd
      type(t_gridfnc3)             :: prf_ext
      type(t_profile)              :: prf_opt
      integer,      intent(in)     :: npoint, is_wheel, icol_kyield, idebug
      logical,      intent(in)     :: fill_spline
      real(kind=8), dimension(:,:), pointer :: points
!--local variables:
      integer,      parameter :: max_num_kinks = 99, max_num_accel = 99
      integer              :: arrsiz, ncolpnt, ipnt, flip_data, nkink, naccel,                          &
                              ikinks(max_num_kinks), iaccel(max_num_accel)
      logical              :: use_wgt, use_minz
      real(kind=8)         :: lambda, ds_bspl

      associate(mirror_y   => prf_opt%mirror_y,   mirror_z   => prf_opt%mirror_z,                       &
                sclfac     => prf_opt%sclfac,     ismooth    => prf_opt%ismooth,                        &
                smooth     => prf_opt%smth,       zig_thrs   => prf_opt%zig_thrs,                       &
                kink_high  => prf_opt%kink_high,  kink_low   => prf_opt%kink_low,                       &
                kink_wid   => prf_opt%kink_wid,   has_kyield => prf_opt%has_kyield,                     &
                ierror     => prf_opt%ierror)

      ! adjust rail profile to avoid fold-back in the top surface

      if (.false. .and. is_wheel.le.0) then
         arrsiz  = size(points,1)
         ncolpnt = size(points,2)
         call remove_vertical_slopes(points, arrsiz, npoint, ncolpnt, idebug)
      endif

      !  - copy data to profil arrays
      !    rail profile: store as 1D cartesian curvilinear grid with nx=1, ny=npoint
      !    (wheel profile: store as 1D cylindrical grid with nv=npoint, nth=1)

      call grid_create_curvil(prf_grd, 1, npoint, lies_in_oyz=.true.)

      do ipnt = 1, npoint
         prf_grd%x(ipnt) = 0d0
         prf_grd%y(ipnt) = points(ipnt,1)
         prf_grd%z(ipnt) = points(ipnt,2)
      enddo

      ! check that profile does not contain any loops, segments crossing each other

      if (ierror.eq.0) then
         call profile_check_bowtie(npoint, prf_grd%y, prf_grd%z, is_wheel, idebug, ierror)
      endif

      ! check that y-values of rail are increasing / wheel-y decreasing

      if (ierror.eq.0) then
         call check_inversion(fname, prf_grd%ntot, prf_grd%coor, 2, 3, is_wheel, idebug, 1, flip_data,  &
                        ierror)
      endif

      ! create grid function for additional data: [k_yield]

      if (ierror.eq.0) then
         if (icol_kyield.ge.1) then
            ! call write_log(' Read_miniprof: profile has k_yield')
            has_kyield = .true.
            call gf3_new(prf_ext, 'ext_data', prf_grd)
   
            do ipnt = 1, npoint
               prf_ext%vn(ipnt) = points(ipnt,icol_kyield)
            enddo
            ! call gf3_print(prf_ext, 'kyield', ikZDIR, 4)
         else
            ! call write_log(' Read_miniprof: no k-yield in profile')
         endif
      endif

      deallocate(points)
      points => NULL()

      ! fill in the s-coordinate, distance along profile

      if (ierror.eq.0) then
         call grid_make_arclength(prf_grd, ierror)
      endif

      ! create smoothing spline, preparing for interpolations

      if (ierror.eq.0 .and. fill_spline) then
         if (idebug.ge.2) then
            write(bufout, '(a,i2,a,g11.3)') ' Creating spline representation using method', ismooth, &
                ', smooth=', smooth
            call write_log(1, bufout)
         endif
         if (idebug.ge.2) then
            write(bufout, '(9x,2(a,f6.3,a,f6.1),a,f6.1)') 'kink_high=', kink_high,' (',kink_high*180d0/pi, &
                ' deg), kink_low=', kink_low,' (',kink_low*180d0/pi,' deg), kink_wid=',kink_wid
            call write_log(1, bufout)
         endif

         if (ismooth.eq.0 .or. ismooth.eq.1) then

            ! ismooth = 0 or 1: using non-weighted or weighted smoothing PP-spline

            ! detect kinks in input profile, including overall start/end-points

            call profile_find_kinks(prf_grd%ntot, prf_grd%y, prf_grd%z, is_wheel, kink_high, kink_low,  &
                           kink_wid, 1d0, nkink, ikinks, idebug, ierror)

            if (nkink.gt.2 .and. ismooth.eq.1) then
               write(bufout,'(a,i3,a,10i5)') ' find_kinks: nkink=',nkink,', ikinks=',ikinks(1:nkink)
               call write_log(1, bufout)
            endif

            if (ismooth.eq.0) then
               use_wgt = .false.
               lambda  = smooth
            else
               use_wgt = .true.
               lambda  = smooth**4 / (16d0 * pi**4)
            endif

            if (ierror.eq.0) then
               use_minz = (is_wheel.eq.0)
               call spline_set_debug(1)
               call grid_make_ppspline(prf_grd, lambda, use_wgt, nkink, ikinks, ierror, k_chk=10)

               if (is_wheel.eq.1 .and. .false.) then
                  call spline_set_debug(1)
                  call spline_add_topview(prf_grd%spl, use_minz, ierror)
                  call spline_print(prf_grd%spl, 'prf_grd%spl', 3)
               endif
               call spline_set_debug(0)
            endif

            ! yout = -29.95d0
            ! call spline_get_s_at_y( prf_grd, yout, sarr(1), ierror )
            ! call spline_eval(prf_grd%spl, ikZDIR, 1, sarr, ierror, exterval=-99d0, f_eval=zarr)
            ! write(bufout,'(3(a,g14.6))') ' y_out =',yout,' found at s_out =',sarr(1),' with z_out =',zarr(1)
            ! call write_log(1, bufout)

         elseif (ismooth.eq.2) then

            if (.false.) then
               call spline_set_debug(3)
               call test_bspline(ierror)
               call spline_set_debug(0)
               call write_log(' done test_bspline...')
               call abort_run()
            endif

            ! detect kinks in input profile, including overall start/end-points

            call profile_find_kinks(prf_grd%ntot, prf_grd%y, prf_grd%z, is_wheel, kink_high, kink_low,  &
                           kink_wid, 1d0, nkink, ikinks, idebug, ierror)

            if (nkink.gt.2) then
               write(bufout,'(a,i3,a,10i5)') ' find_kinks: nkink=',nkink,', ikinks=',ikinks(1:nkink)
               call write_log(1, bufout)
            endif

            ! ismooth = 2: using weighted smoothing B-spline

            use_wgt         = .true.
            lambda          = smooth**6 / (64d0 * pi**6)
            ds_bspl         = 0.001d0
            naccel          = 0
            iaccel(1)       = 0
            use_minz        = (is_wheel.eq.0)

            ! call write_log(' call grid_make_bspline...')
            call spline_set_debug(0)
            call grid_make_bspline(prf_grd, ds_bspl, lambda, use_wgt, nkink, ikinks, naccel, iaccel,    &
                                ierror, k_chk=10)

            if (is_wheel.eq.1 .and. .false.) then
               call spline_set_debug(1)
               call spline_add_topview(prf_grd%spl, use_minz, ierror)
               call spline_print(prf_grd%spl, 'prf_grd%spl', 3)
            endif
            call spline_set_debug(0)
            ! call spline_print(prf_grd%spl, 'prf_grd%spl', 3)
            ! call write_log(' done grid_make_bspline...')

            ! yout = -30.02d0
            ! call spline_set_debug(0)
            ! call spline_get_s_at_y( prf_grd, yout, sarr(1), ierror )
            ! call spline_eval(prf_grd%spl, ikZDIR, 1, sarr, ierror, exterval=-99d0, f_eval=zarr)
            ! write(bufout,'(3(a,g14.6))') ' y_out =',yout,' found at s_out =',sarr(1),' with z_out =',zarr(1)
            ! call write_log(1, bufout)
            ! call spline_set_debug(0)

         else

            write(bufout,'(a,i4,a)') ' Internal error: smoothing method out of range (',ismooth,')'
            call write_log(1, bufout)
            call abort_run()

         endif ! case ismooth
      endif ! ierror=0

      ! print profile when needed

      if (idebug.ge.4) then
         call grid_print(prf_grd, fname, 5)
      endif

      end associate
   end subroutine profile_finish_grid

!------------------------------------------------------------------------------------------------------------

   subroutine read_simpack_profile(fname, fulnam, npoint, points, is_wheel, inp_mirror_y, inp_mirror_z, &
                        inp_sclfac, inp_zig_thrs, x_profil, x_readln, lstop, ierror)
!--purpose: Read a wheel or rail profile from a Simpack prr or prw file
      implicit none
!--subroutine arguments:
      character*(*)             :: fname, fulnam
      real(kind=8), dimension(:,:), pointer :: points ! npoint x 3:  (y,z,wgt)
      integer,      intent(out) :: npoint
      integer,      intent(in)  :: is_wheel, inp_mirror_y, inp_mirror_z, x_profil, x_readln
      real(kind=8), intent(in)  :: inp_sclfac, inp_zig_thrs
      logical,      intent(in)  :: lstop
      integer                   :: ierror
!--local variables:
      integer,     parameter :: mxitem = 20, ncolpnt = 3
      character*2, parameter :: commnt = '!%'
      integer            :: unitnm, ios, linenr, nblank, ncmtln, ieof, nitem, ipnt
      logical            :: has_key, has_str, has_match, in_header, in_spline, in_point
      character*256      :: inptxt, keywrd, valstr
      real(kind=8)       :: values(1:mxitem)
      integer            :: arrsiz, iversion, my_wheel, flip_data, file_mirror_y, file_mirror_z
      real(kind=8)       :: approx_smooth, point_dist_min, shift_y, shift_z, rotate, bound_y_min,       &
                            bound_y_max, bound_z_min, bound_z_max, units_len_fac, units_ang_fac

      if (x_profil.ge.2) then
         write(bufout,'(/,2a)') ' --- read_spck: file ',trim(fulnam)
         call write_log(2, bufout)
      endif

      ierror = 0
      unitnm = ltmp
      open(unitnm, file=fulnam, status='old', iostat=ios)
      if (ios.ne.0) then
         write(bufout,'(3a,i4,a)') ' ERROR: cannot open file "',trim(fulnam),'" (',ios,')'
         call write_log(1, bufout)
         ierror = ios
      endif
      linenr = 0

      ! allocate points-array at initial size

      npoint         =  0
      arrsiz         = 1000
      call reallocate_arr(points, arrsiz, ncolpnt)

      iversion       = -1
      my_wheel       = -1
      approx_smooth  =  0d0
      point_dist_min =  1d-4
      shift_y        =  0d0
      shift_z        =  0d0
      rotate         =  0d0
      bound_y_min    = -1d9
      bound_y_max    =  1d9
      bound_z_min    = -1d9
      bound_z_max    =  1d9
      file_mirror_y  =  0
      file_mirror_z  =  0
      flip_data      = -99      ! -99 = not specified
      units_len_fac  =  1d0
      units_ang_fac  =  1d0

      in_header      = .false.
      in_spline      = .false.
      in_point       = .false.
      ieof           = 0

      do while(ieof.le.0 .and. ierror.eq.0)

         ! Increment line number and read input, thereby skipping comments and empty lines

         call getNonemptyLine(unitnm, 0, 'line of simpack profile file', commnt, x_readln, ieof, lstop,   &
                linenr, nblank, ncmtln, inptxt)

         ! split the input at spaces and put into string array

         if (ieof.le.0) then
            call readKeywordValues(inptxt, commnt, mxitem, x_readln, has_key, keywrd, has_str, valstr,    &
                nitem, values)
         else
            ierror = 100
            if (in_point) then
               write(bufout,'(4a)') ' ERROR 100: unexpected end of file in point-section of prr/prw-',  &
                        'file "', trim(fname),'"'
            elseif (in_spline) then
               write(bufout,'(4a)') ' ERROR 100: unexpected end of file in spline-section of prr/prw-', &
                        'file "', trim(fname),'"'
            elseif (in_header) then
               write(bufout,'(4a)') ' ERROR 100: unexpected end of file in header-section of prr/prw-', &
                        'file "', trim(fname),'"'
            else
               write(bufout,'(3a)') ' ERROR 100: unexpected end of file in prr/prw-file "', trim(fname),'"'
            endif
            call write_log(1, bufout)

            has_key = .false.
            has_str = .false.
            inptxt  = ' '
            nitem   = 0
         endif

         ! process start/end of sections
         !       header
         !       spline
         !       spline - point

         if (has_key .and. .not.has_str) then
            has_match = .true.
            if (keywrd.eq.'header.begin') then
               if (in_header  .or. in_spline .or. in_point) ierror = 101
               in_header = .true.
            elseif (keywrd.eq.'header.end') then
               if (.not.in_header  .or. in_spline .or. in_point) ierror = 101
               in_header = .false.
            elseif (keywrd.eq.'spline.begin') then
               if (in_spline .or. in_header  .or. in_point) ierror = 101
               in_spline = .true.
            elseif (keywrd.eq.'spline.end') then
               if (.not.in_spline .or. in_header .or. in_point) ierror = 101
               in_spline = .false.
            elseif (keywrd.eq.'point.begin') then
               if (in_point .or. .not.in_spline .or. in_header) ierror = 101
               in_point = .true.
            elseif (keywrd.eq.'point.end') then
               if (.not.in_point .or. .not.in_spline .or. in_header) ierror = 101
               in_point = .false.
            else
               has_match = .false.
               ierror = 102
            endif
            if (ierror.eq.101) then
               write(bufout,'(3a)') ' ERROR 101: incorrect nesting of sections in prr/prw-file ("',     &
                        trim(fname),'")'
               call write_log(1, bufout)
            elseif (ierror.eq.102) then
               write(bufout,'(5a)') ' ERROR 102: unknown section name "',trim(keywrd),                  &
                        '" in prr/prw-file "', trim(fname),'")'
               call write_log(1, bufout)
            endif
            if (x_profil.ge.3 .and. has_match) then
               write(bufout,*) '   ...found section ',trim(keywrd)
               call write_log(1, bufout)
            endif
         endif

         ! process keywords of section 'header':

         if (has_key .and. has_str .and. in_header) then
            if (nitem.ne.1) then
               ierror = 111
               write(bufout,'(3a,i3,3a)') ' ERROR 111: Keyword ',trim(keywrd),                          &
                        ' requires one integer value, found',nitem,' ("',trim(fname),'")'
               call write_log(1, bufout)
            elseif (keywrd.eq.'version') then
               iversion = nint(values(1))
               if (x_profil.ge.3) then
                  write(bufout,*) '      ...found version=',iversion
                  call write_log(1, bufout)
               endif
            elseif (keywrd.eq.'type') then
               my_wheel = nint(values(1))
               if (my_wheel.ne.0 .and. my_wheel.ne.1) ierror = 121
               if (my_wheel.ne.is_wheel .and. is_wheel.eq.0) then
                  write(bufout,'(3a)') ' ERROR 112: Type must be 0 for a prr-file ("',trim(fname),'")'
                  call write_log(1, bufout)
                  ierror = 112
               elseif (my_wheel.ne.is_wheel .and. is_wheel.eq.1) then
                  write(bufout,'(3a)') ' ERROR 113: Type must be 1 for a prw-file ("',trim(fname),'")'
                  call write_log(1, bufout)
                  ierror = 113
               endif
               if (x_profil.ge.3) then
                  write(bufout,*) '      ...found is_wheel=',my_wheel
                  call write_log(1, bufout)
               endif
            endif
            if (ierror.ge.121 .and. ierror.le.129) then
               write(bufout,'(a,i4,6a)') ' ERROR',ierror,': incorrect data in header section, ',     &
                     'keyword', trim(keywrd),' ("',trim(fname),'")'
               call write_log(1, bufout)
            endif
         endif

         ! process keywords of section 'spline', excluding 'point':

         if (has_key .and. has_str .and. in_spline .and. .not.in_point) then
            if (keywrd.eq.'approx.smooth') then
               if (nitem.ne.1) then
                  ierror = 131
               else
                  approx_smooth = values(1)
               endif
            elseif (keywrd.eq.'file') then
            elseif (keywrd.eq.'file.mtime') then
            elseif (keywrd.eq.'comment') then
            elseif (keywrd.eq.'type') then
            elseif (keywrd.eq.'point.dist.min') then
               point_dist_min = max(1d-6, values(1))
            elseif (keywrd.eq.'shift.y') then
               shift_y = values(1)
            elseif (keywrd.eq.'shift.z') then
               shift_z = values(1)
            elseif (keywrd.eq.'rotate') then
               rotate  = values(1)
            elseif (keywrd.eq.'inversion') then
               flip_data = nint(values(1))
            elseif (keywrd.eq.'bound.y.min') then
               bound_y_min = values(1)
            elseif (keywrd.eq.'bound.y.max') then
               bound_y_max = values(1)
            elseif (keywrd.eq.'bound.z.min') then
               bound_z_min = values(1)
            elseif (keywrd.eq.'bound.z.max') then
               bound_z_max = values(1)
            elseif (keywrd.eq.'mirror.y') then
               file_mirror_y  = nint(values(1))
            elseif (keywrd.eq.'mirror.z') then
               file_mirror_z  = nint(values(1))
            elseif (keywrd.eq.'units.len') then
               if (x_profil.ge.1) call write_log(' The keyword "units.len" is ignored. Please use ' //    &
                        '"units.len.f" instead')
            elseif (keywrd.eq.'units.ang') then
               if (x_profil.ge.1) call write_log(' The keyword "units.ang" is ignored. Please use ' //    &
                        '"units.ang.f" instead')
            elseif (keywrd.eq.'units.len.f') then
               units_len_fac = values(1)
            elseif (keywrd.eq.'units.ang.f') then
               units_ang_fac = values(1)
            else
               write(bufout,*) 'unknown keyword "',trim(keywrd),'" in spline section, ignored.'
               call write_log(1, bufout)
            endif

            if (ierror.ge.131 .and. ierror.le.139) then
               write(bufout,'(a,i4,3a)') ' ERROR',ierror,': incorrect data in spline section, ',     &
                     'keyword', trim(keywrd)
               call write_log(1, bufout)
            endif
         endif

         ! process keywords of section 'spline' - 'point':

         if (has_key .and. has_str .and. in_spline .and. in_point) then
            write(bufout,*) 'unknown keyword "',trim(keywrd),'" in point section, ignored.'
            call write_log(1, bufout)
         endif

         ! process values of section 'spline' - 'point':

         if (.not.has_key .and. has_str .and. in_spline .and. in_point) then
            if (nitem.lt.1) then
               write(bufout,*) 'too few data values',nitem
               call write_log(1, bufout)
            elseif (nitem.gt.3) then
               write(bufout,*) 'too many data values',nitem
               call write_log(1, bufout)
            else
               npoint = npoint + 1

               ! enlarge array if more points found than current arrsiz

               if (npoint.gt.arrsiz) then
                  arrsiz = int(1.5d0 * arrsiz)
                  call reallocate_arr(points, arrsiz, ncolpnt, keep=.true.)
               endif

               ! store current data-values: y-value, z-value, [spline wgt]

               points(npoint,1) = values(1)
               points(npoint,2) = values(2)
               if (nitem.ge.3) then
                  points(npoint,3) = values(3)
               else
                  points(npoint,3) = 0d0
               endif
            endif
         endif

         ! recognize end of data

         if (index(inptxt, 'spline.end').gt.0) then
            if (x_profil.ge.3) call write_log(' found spline end.')
            ieof = 1
         endif

      ! end while (more data, error==0)

      enddo
      close(unitnm)

      ! Copy data points to either the wheel or the rail profile arrays

      if (is_wheel.ne.0 .and. is_wheel.ne.1) then
         ierror = 191
         write(bufout,'(3a)') ' ERROR 191: incorrect value for flag header - type ("',trim(fname),'")'
         call write_log(1, bufout)
      endif

      if (ierror.ne.0) then
         if (x_profil.ge.1) call write_log(' Errors found, skipping profile processing.')
      else
         if (x_profil.ge.3) then
            write(bufout,'(a,i6,a)') ' Obtained',npoint,' profile points'
            call write_log(1, bufout)
         endif

         !  - modification 1: point.dist.min - remove points too close together

         call filter_close_points(points, arrsiz, npoint, ncolpnt, point_dist_min, is_wheel, x_profil)

         !  - modification 1': remove zig-zag-patterns with double kink

         call filter_double_kinks(points, arrsiz, npoint, ncolpnt, is_wheel, inp_zig_thrs, x_profil)

         !  - modification 2: shift.y - in user units
         !  - modification 3: shift.z

         do ipnt = 1, npoint
            points(ipnt,1) = points(ipnt,1) + shift_y
            points(ipnt,2) = points(ipnt,2) + shift_z
         enddo

         !  - modification 10: units.ang.f - angle scaling factor [user] / [rad]

         rotate = rotate / units_ang_fac

         !  - modification 4: rotate  - rotate profile by prescribed rotation angle

         if (abs(rotate).gt.1d-6) then
            ierror = 182
            write(bufout,'(3a)') ' ERROR 182: profile rotation not yet implemented ("',trim(fname),'")'
            call write_log(1, bufout)
         endif

         !  - modification 5: bound.y.min, bound.y.max - removing points outside interval

         if (bound_y_min.lt.bound_y_max) then
            call delete_outside_boundbox(points, arrsiz, npoint, ncolpnt, bound_y_min, bound_y_max,     &
                         1d0, 0d0, x_profil)
         endif

         !  - modification 6: bound.z.min, bound.z.max - clipping - shifting points to interval boundary

         if (bound_z_min.lt.bound_z_max) then
            do ipnt = 1, npoint
               points(ipnt,2) = max(bound_z_min, min(bound_z_max, points(ipnt,2)))
            enddo
         endif

         !  - modification 7: mirror.y

         if (inp_mirror_y+file_mirror_y.eq.1) then
            do ipnt = 1, npoint
               points(ipnt,1) = -points(ipnt,1)
            enddo
         elseif (inp_mirror_y+file_mirror_y.ge.2) then
            call write_log(' WARNING: y-mirroring of profile file is undone by mirroring requested in input.')
         endif

         !  - modification 8: mirror.z

         if (inp_mirror_z+file_mirror_z.eq.1) then
            do ipnt = 1, npoint
               points(ipnt,2) = -points(ipnt,2)
            enddo
         elseif (file_mirror_z.eq.1 .and. abs(inp_mirror_z).eq.1) then
            call write_log(' WARNING: z-mirroring of profile file is undone by mirroring requested in input.')
         endif

         !  - modification 10: units.len.f - length scaling factor
         !    SIMPACK converts profile to [m], we add 1000 to get the end result in [mm]

         do ipnt = 1, npoint
            points(ipnt,1) = inp_sclfac * 1000d0 * points(ipnt,1) / units_len_fac
            points(ipnt,2) = inp_sclfac * 1000d0 * points(ipnt,2) / units_len_fac
         enddo

         !  - modification 9: inversion

         if (ierror.eq.0 .and. flip_data.le.-99) then
            call check_inversion(fname, npoint, points, 1, 2, is_wheel, x_profil, 0, flip_data, ierror)
         endif

         if (ierror.eq.0 .and. flip_data.eq.1) then
            call invert_points(npoint, points, is_wheel, x_profil)
         endif

         !  - first/last points should now be in appropriate order

         if (ierror.eq.0) then
            call check_overall_order(npoint, points, is_wheel, x_profil, ierror)
         endif

      endif ! ierror<>0

   end subroutine read_simpack_profile

!------------------------------------------------------------------------------------------------------------

   subroutine read_miniprof_profile(fname, fulnam, iftype, npoint, points, is_wheel, mirror_y,          &
                        inp_mirror_z, inp_sclfac, inp_zig_thrs, icol_kyield, x_profil, x_readln,        &
                        lstop, ierror)
!--purpose: Read a wheel or rail profile from a Miniprof ban or whl file
      implicit none
!--subroutine arguments:
      character*(*)             :: fname, fulnam
      integer,      intent(out) :: npoint
      real(kind=8), dimension(:,:), pointer :: points ! npoint x 6: [y, z, angle, curvature, wgt, kyield]
      integer,      intent(in)  :: iftype, is_wheel, mirror_y, inp_mirror_z, x_profil, x_readln
      real(kind=8), intent(in)  :: inp_sclfac, inp_zig_thrs
      integer,      intent(out) :: icol_kyield
      logical,      intent(in)  :: lstop
      integer                   :: ierror
!--local variables:
      integer,      parameter :: mxitem = 20, ncolpnt = 6, mxwarn_coldef = 4
      character*4,  parameter :: commnt = '!"%#'
      real(kind=8), parameter :: point_dist_min = 1d-4
      integer             :: unitnm, ios, linenr, line_prv, nblank, ncmtln, in_headers, old_style,      &
                             ieof, nwarn_coldef, nitem, i, ipnt
      logical             :: has_key, has_str
      character*256       :: inptxt, keywrd, valstr, columndef
      real(kind=8)        :: values(1:mxitem)
      integer             :: arrsiz, xypoints, mirror_z, flip_data, ncolfile, icolfile, ipnt2file(ncolpnt)
      real(kind=8)        :: xoffset, yoffset, grade, superelevation, gauge

      if (x_profil.ge.2) then
         write(bufout,'(/,2a)') '--- read_miniprof: file ',trim(fulnam)
         call write_log(2, bufout)
      endif
      icol_kyield = 0

      unitnm = ltmp
      ierror = 0
      open(unitnm, file=fulnam, status='old', iostat=ios)
      if (ios.ne.0) then
         write(bufout,'(3a,i4,a)') ' ERROR: cannot open file "',trim(fulnam),'" (',ios,')'
         call write_log(1, bufout)
         ierror = ios
      endif
      linenr = 0

      ! allocate points-array at initial size

      npoint         =  0                       ! number of points in points array
      arrsiz         = 1000                     ! length of points array
      call reallocate_arr(points, arrsiz, ncolpnt)

      xypoints  = -1                            ! number of points according to file
      columndef = 'X,Y,*'                       ! identifiers X,Y,A,C,N,K for the columns
      ncolfile  = 2                             ! number of columns in file
      ipnt2file = (/  1,  2, -1, -1, -1, -1 /)  ! for each column in points: column number in file
      xoffset        =  0d0
      yoffset        =  0d0
      grade          =  0d0
      superelevation =  0d0
      gauge          =  0d0

      old_style  = 0  ! -1 = no, 1 = yes, 0 = dont know
      in_headers = 1
      nwarn_coldef = 0
      ieof = 0
      do while(ieof.le.0 .and. ierror.eq.0)

         ! Increment line number and read input, thereby skipping comments and empty lines

         line_prv = linenr
         call getNonemptyLine(unitnm, 0, 'line of miniprof file', commnt, x_readln, ieof, lstop, linenr,  &
                        nblank, ncmtln, inptxt)

         ! Miniprof data start after an empty line

         if (nblank.ge.1) in_headers = 0

         ! Split the input at '=' sign, get keyword, string-value, numerical value(s)
         !    case 1: keyword = value (string/numeric)
         !    case 2: keyword, or keyword = <empty>
         !    case 3: numerical value(s)

         call readKeywordValues(inptxt, commnt, mxitem, x_readln, has_key, keywrd, has_str, valstr,       &
                        nitem, values)

         ! Provide support for basic 2-column X,Y values and for miniprof files with missing blank line

         if (in_headers.eq.1 .and. .not.has_key .and. old_style.le.0) then

            if (x_profil.ge.-1 .and. iftype.ne.FTYPE_PLAIN2D) then
               write(bufout,'(a,i5,a)') ' WARNING: found data at line', linenr,' not preceded by blank line'
               call write_log(1, bufout)
            endif

            in_headers = 0
         endif

         ! Handle all possible cases

         if (ieof.gt.0) then

            ! end-of-file: no more input

         elseif (in_headers.eq.0 .and. has_key) then

            ierror = 189
            write(bufout,'(3a,i5,a)') ' ERROR (189): found key "',trim(keywrd),'" at line',linenr,      &
                ' after data values on previous lines'
            call write_log(1, bufout)

         elseif (in_headers.eq.1 .and. .not.has_key) then

            if (x_profil.ge.1) then
               write(bufout,'(a,i5)') ' WARNING: found data while still in headers, ignoring line', linenr
               call write_log(1, bufout)
            endif

         elseif (in_headers.eq.1 .and. has_key .and. .not.has_str .and. nitem.le.0) then

            ! plain keyword (old-style Miniprof file?)

            if (x_profil.ge.3) then
               write(bufout,*) 'Plain keyword "',trim(keywrd),'", ignored.'
               call write_log(1, bufout)
            endif
            old_style = 1

         elseif (has_key) then

            ! has_key: process keywords (case sensitive), given here in alphabetical order:

            if (keywrd.eq.'AlarmFailureCount') then
            elseif (keywrd.eq.'AlarmWarningCount') then
            elseif (keywrd.eq.'AlignDA') then
            elseif (keywrd.eq.'AlignDX') then
            elseif (keywrd.eq.'AlignDY') then
            elseif (keywrd.eq.'AlignOX') then
            elseif (keywrd.eq.'AlignOY') then
            elseif (keywrd.eq.'AlignRX') then
            elseif (keywrd.eq.'AlignRY') then
            elseif (keywrd.eq.'AlignSX') then
            elseif (keywrd.eq.'AlignSY') then
            elseif (keywrd.eq.'AreaGain') then
            elseif (keywrd.eq.'AreaLoss') then
            elseif (keywrd.eq.'AreaRef') then
            elseif (keywrd.eq.'AxleNo') then
            elseif (keywrd.eq.'BaseOffset') then
            elseif (keywrd.eq.'BaseTilt') then
            elseif (keywrd.eq.'Battery.Capacity') then
            elseif (keywrd.eq.'Battery.Level') then
            elseif (keywrd.eq.'Bogie') then
            elseif (keywrd.eq.'CarNo') then
            elseif (keywrd.eq.'Chainage') then
            elseif (keywrd.eq.'Chars') then
            elseif (keywrd.eq.'ColumnDef') then
               columndef = valstr
               ncolfile = (len(trim(columndef)) + 1) / 2
               if (x_profil.ge.3) then
                  write(bufout,*) 'columndef = "',trim(columndef),'", ncol=',ncolfile
                  call write_log(1, bufout)
               endif
               do icolfile = 1, ncolfile
                  i = 2*icolfile-1
                  if (columndef(i:i).eq.'X') then
                     ipnt2file(1) = icolfile    ! store X-coordinates in points(:,1)
                  elseif (columndef(i:i).eq.'Y') then
                     ipnt2file(2) = icolfile    ! store Y-coordinates in points(:,2)
                  elseif (columndef(i:i).eq.'A') then
                     ipnt2file(3) = icolfile    ! store angles A      in points(:,3)
                  elseif (columndef(i:i).eq.'C') then
                     ipnt2file(4) = icolfile    ! store curvatures C  in points(:,4)
                  elseif (columndef(i:i).eq.'N') then
                     ipnt2file(5) = icolfile    ! store ??? N         in points(:,5)
                  elseif (columndef(i:i).eq.'K') then
                     ipnt2file(6) = icolfile    ! store yield K       in points(:,6)
                     icol_kyield  = 6
                  else
                     write(bufout,*) 'WARNING: Incomprehensible input in columndef="',trim(columndef),'"'
                     call write_log(1, bufout)
                  endif
               enddo
               if (x_profil.ge.3) then
                  write(bufout,*)'ipnt2file=',(ipnt2file(i), i=1,ncolpnt)
                  call write_log(1, bufout)
               endif
            elseif (keywrd.eq.'Comment') then
            elseif (keywrd.eq.'CrownR') then
            elseif (keywrd.eq.'CrownRadius') then
            elseif (keywrd.eq.'Curve') then
            elseif (keywrd.eq.'Date') then
            elseif (keywrd.eq.'DiameterFlange') then
            elseif (keywrd.eq.'DiameterFlange_AlarmStatus') then
            elseif (keywrd.eq.'DiameterTaperline') then
            elseif (keywrd.eq.'DiameterTaperline_AlarmStatus') then
            elseif (keywrd.eq.'Direction') then
            elseif (keywrd.eq.'EmplNum') then
            elseif (keywrd.eq.'Elevation') then
            elseif (keywrd.eq.'Field ID') then
            elseif (keywrd.eq.'Filename') then
            elseif (keywrd.eq.'Flag') then
            elseif (keywrd.eq.'FlangeAngleMax') then
            elseif (keywrd.eq.'FlangeAngleMax_AlarmStatus') then
            elseif (keywrd.eq.'FlangeAngleMaxPos') then
            elseif (keywrd.eq.'FlangeAngleMaxPos_AlarmStatus') then
            elseif (keywrd.eq.'FlangeWidth') then
            elseif (keywrd.eq.'FlangeWidth_AlarmStatus') then
            elseif (keywrd.eq.'Flat') then
            elseif (keywrd.eq.'Gauge') then
               gauge = values(1)
            elseif (keywrd.eq.'Gauge_AlarmStatus') then
            elseif (keywrd.eq.'GaugeRodLength') then
            elseif (keywrd.eq.'GaugeSensorDistance') then
            elseif (keywrd.eq.'GaugeSensorSamples') then
            elseif (keywrd.eq.'GaugeSensorStdDev') then
            elseif (keywrd.eq.'Grade') then
               grade = values(1)
            elseif (keywrd.eq.'GradeRatio') then
            elseif (keywrd.eq.'Hollowing') then
            elseif (keywrd.eq.'Hollowing_AlarmStatus') then
            elseif (keywrd.eq.'HollowingPos') then
            elseif (keywrd.eq.'HollowingPos_AlarmStatus') then
            elseif (keywrd.eq.'ID') then
            elseif (keywrd.eq.'Instrument.Firmware') then
            elseif (keywrd.eq.'Instrument.Connection') then
            elseif (keywrd.eq.'KP') then
            elseif (keywrd.eq.'Line') then
            elseif (keywrd.eq.'Location.Altitude') then
            elseif (keywrd.eq.'Location.HorizontalAccuracy') then
            elseif (keywrd.eq.'Location.Latitude') then
            elseif (keywrd.eq.'Location.Longitude') then
            elseif (keywrd.eq.'Location.VerticalAccuracy') then
            elseif (keywrd.eq.'Mileage') then
            elseif (keywrd.eq.'MPCalDate') then
            elseif (keywrd.eq.'MPCalTime') then
            elseif (keywrd(1:3).eq.'MPD') then
               ! MPDC, MPDN, MPD1--MPD99 - could use verify to check for known characters
            elseif (keywrd.eq.'MPSerNo') then
            elseif (keywrd.eq.'MPTypeNo') then
            elseif (keywrd.eq.'OriginalHeader') then
            elseif (keywrd.eq.'Position') then
            elseif (keywrd.eq.'ProfileAlignment') then
            elseif (keywrd.eq.'ProfileShaping') then
            elseif (keywrd.eq.'ProgramDate') then
            elseif (keywrd.eq.'ProgramName') then
            elseif (keywrd.eq.'ProgramVer') then
            elseif (keywrd.eq.'qR') then
            elseif (keywrd.eq.'qR_AlarmStatus') then
            elseif (keywrd.eq.'Rail') then
            elseif (keywrd.eq.'RailAlignType') then
            elseif (keywrd.eq.'RailAngle') then
            elseif (keywrd.eq.'RailAngle_AlarmStatus') then
            elseif (keywrd.eq.'RailGaugePoint') then
            elseif (keywrd.eq.'RailRodLength') then
            elseif (keywrd.eq.'RCF') then
            elseif (keywrd.eq.'ReferenceProfile') then
            elseif (keywrd.eq.'RefPoint1') then
            elseif (keywrd.eq.'RefPoint2') then
            elseif (keywrd.eq.'RefPoint3') then
            elseif (keywrd.eq.'RefPoint4') then
            elseif (keywrd.eq.'RefPoint5') then
            elseif (keywrd.eq.'RefPoint6') then
            elseif (keywrd.eq.'RefPoint7') then
            elseif (keywrd.eq.'RefPoint8') then
            elseif (keywrd.eq.'RefPoint9') then
            elseif (keywrd.eq.'RodLength') then
            elseif (keywrd.eq.'Sd') then
            elseif (keywrd.eq.'Sd_AlarmStatus') then
            elseif (keywrd.eq.'SectionNo') then
            elseif (keywrd.eq.'Sh') then
            elseif (keywrd.eq.'Sh_AlarmStatus') then
            elseif (keywrd.eq.'Side') then
            elseif (keywrd.eq.'Stock') then
            elseif (keywrd.eq.'SuperElevation') then
               superelevation = values(1)
            elseif (keywrd.eq.'SuperElevationHeight') then
            elseif (keywrd.eq.'Temperature') then
            elseif (keywrd.eq.'Tilt') then
            elseif (keywrd.eq.'Time') then
            elseif (keywrd.eq.'Track') then
            elseif (keywrd.eq.'TrackSection') then
            elseif (keywrd.eq.'TransformState') then
            elseif (keywrd.eq.'Transformation') then
            elseif (keywrd.eq.'UOF1') then
            elseif (keywrd.eq.'UOF2') then
            elseif (keywrd.eq.'UOF3') then
            elseif (keywrd.eq.'UOF4') then
            elseif (keywrd.eq.'UserName') then
            elseif (keywrd.eq.'W1') then
            elseif (keywrd.eq.'W1_AlarmLevel') then
            elseif (keywrd.eq.'W1_AlarmStatus') then
            elseif (keywrd.eq.'W2') then
            elseif (keywrd.eq.'W2_AlarmLevel') then
            elseif (keywrd.eq.'W2_AlarmStatus') then
            elseif (keywrd.eq.'W3') then
            elseif (keywrd.eq.'W3_AlarmLevel') then
            elseif (keywrd.eq.'W3_AlarmStatus') then
            elseif (keywrd.eq.'WheelAdjustPoint') then
            elseif (keywrd.eq.'WheelDiameterTaperline') then
            elseif (keywrd.eq.'WheelID') then
            elseif (keywrd.eq.'Windows.ComputerName') then
            elseif (keywrd.eq.'Windows.NetFramework') then
            elseif (keywrd.eq.'Windows.ServicePack') then
            elseif (keywrd.eq.'Windows.UserName') then
            elseif (keywrd.eq.'Windows.Version') then
            elseif (keywrd.eq.'Xoffset' .or. keywrd.eq.'XOffset') then
               xoffset = values(1)
            elseif (keywrd.eq.'XYPoints') then
               xypoints = nint(values(1))
            elseif (keywrd.eq.'Yoffset' .or. keywrd.eq.'YOffset') then
               yoffset = values(1)
            elseif (x_profil.ge.2) then
               write(bufout,*) 'unknown keyword "',trim(keywrd),'", ignored.'
               call write_log(1, bufout)
            endif

         elseif (nitem.lt.ncolfile) then

            ! no key: too few data values

            nwarn_coldef = nwarn_coldef + 1
            if (nwarn_coldef.le.mxwarn_coldef) then
               write(bufout,'(a,i3,3a,i5)') ' ERROR: too few data values (',nitem,') for ColumnDef ',   &
                        trim(columndef), ', ignoring line',linenr
               call write_log(1, bufout)
            endif

         elseif (nitem.gt.ncolfile .and. columndef.ne.'X,Y,*') then

            ! no key: too many data values

            nwarn_coldef = nwarn_coldef + 1
            if (nwarn_coldef.le.mxwarn_coldef) then
               write(bufout,'(a,i3,3a,i5)') ' ERROR: too many data values (',nitem,') for ColumnDef ',  &
                        trim(columndef), ', ignoring line',linenr
               call write_log(1, bufout)
            endif

         else

            ! no key: store data values:

            npoint = npoint + 1

            ! enlarge array if more points found than current arrsiz

            if (npoint.gt.arrsiz) then
               arrsiz = int(1.5d0 * arrsiz)
               call reallocate_arr(points, arrsiz, ncolpnt, keep=.true.)
            endif

            ! adjust column definition if not specified before

            if (nitem.gt.ncolfile .and. columndef.eq.'X,Y,*') then
               if (x_profil.ge.1) then
                  if (nitem.eq.3) then
                     write(bufout,'(a,i3,a,/,9x,3a)') ' WARNING: unspecified ColumnDef with',nitem,     &
                           ' values per line.',' Assuming "', trim(columndef), '", ignoring last column.'
                  else
                     write(bufout,'(a,i3,a,/,9x,3a,i2,a)') ' WARNING: unspecified ColumnDef with',nitem, &
                           ' values per line.',                                                         &
                           ' Assuming "', trim(columndef), '", ignoring last', ncolfile-2,' columns.'
                  endif
                  call write_log(2, bufout)
               endif
               ncolfile = min(ncolpnt, nitem)
            endif

            ! store current data-values: y-value, z-value, [angle], [curvature], [spline wgt], [yield]

            if (ipnt2file(1).ge.1) points(npoint,1) = values(ipnt2file(1))
            if (ipnt2file(2).ge.1) points(npoint,2) = values(ipnt2file(2))
            if (ipnt2file(3).ge.1) points(npoint,3) = values(ipnt2file(3))
            if (ipnt2file(4).ge.1) points(npoint,4) = values(ipnt2file(4))
            if (ipnt2file(5).ge.1) points(npoint,5) = values(ipnt2file(5))
            if (ipnt2file(6).ge.1) points(npoint,6) = values(ipnt2file(6))
            ! write(*,*) ' ky(',npoint,')=',points(npoint,6)
         endif

         ! put content of string array in ints, dbles and flags

         ! write(bufout,*) nitem,' words:',(' ',trim(words(i)),',', i=1,nitem)
         ! call write_log(1, bufout)

      enddo ! while (more data, error==0)
      close(unitnm)

      if (nwarn_coldef.gt.mxwarn_coldef) then
         write(bufout,'(a,i5,a)') ' ERROR: too few data values for ColumnDef, ignoring', nwarn_coldef,  &
                ' points on the profile'
         call write_log(1, bufout)
      endif

      if (x_profil.ge.3) then
         write(bufout,'(a,i6,a)') ' Obtained',npoint,' profile points'
         call write_log(1, bufout)
      endif
      if (x_profil.ge.4) then
         do ipnt = 1, npoint
            write(bufout, '(i5,2(a,f9.3))') ipnt,': y=',points(ipnt,1),', z=',points(ipnt,2)
            call write_log(1, bufout)
         enddo
      endif

      if (ierror.eq.0) then

         ! Check that any points are obtained

         if (npoint.le.1) then
            ierror = 190
            write(bufout,'(2a,i4,a)') ' ERROR 190: invalid profile. At least 2 points are needed, ',    &
                   'obtained', npoint,' points.'
            call write_log(1, bufout)
         endif

         if (is_wheel.ne.0 .and. is_wheel.ne.1) then
            ierror = 191
            write(bufout,'(3a)') ' ERROR 191: incorrect value for flag header - type ("',trim(fname),'")'
            call write_log(1, bufout)
         endif

      endif

      if (ierror.ne.0) then

         if (x_profil.ge.1) call write_log(' Errors found, skipping profile processing.')

      else

         ! Copy data points to either the wheel or the rail profile arrays

         !  - modification 1: point.dist.min - remove points too close together

         ! call write_log('...filter_close_points')
         call filter_close_points(points, arrsiz, npoint, ncolpnt, point_dist_min, is_wheel, x_profil)

         !  - modification 1': remove zig-zag-patterns with double kink

         call filter_double_kinks(points, arrsiz, npoint, ncolpnt, is_wheel, inp_zig_thrs, x_profil)

         !  - modification 2: xoffset - in [mm]
         !  - modification 3: yoffset - in [mm]

         ! call write_log('...shifting')
         do ipnt = 1, npoint
            points(ipnt,1) = points(ipnt,1) + xoffset
            points(ipnt,2) = points(ipnt,2) + yoffset
         enddo

         !  - modification 7: mirror y if requested in input-file

         ! call write_log('...mirroring')
         if (mirror_y.ge.1) then
            do ipnt = 1, npoint
               points(ipnt,1) = -points(ipnt,1)
            enddo
         endif

         !  - modification 8: mirror z if necessary to get z positive downwards

         ! call write_log('...check_z_pos_down')
         if (inp_mirror_z.ne.0) then
            mirror_z = inp_mirror_z
         else
            call check_z_pos_down(fname, npoint, points, is_wheel, x_profil, mirror_z)
         endif

         if (mirror_z.eq.1) then
            ! call write_log('...mirror_z')
            do ipnt = 1, npoint
               points(ipnt,2) = -points(ipnt,2)
            enddo
            if (x_profil.ge.4) then
               call write_log(' after mirror_z:')
               do ipnt = 1, npoint
                  write(bufout, '(i5,2(a,f9.3))') ipnt,': y=',points(ipnt,1),', z=',points(ipnt,2)
                  call write_log(1, bufout)
               enddo
            endif
         endif

         !  - modification 10: length scaling factor: Miniprof data are in [mm]
         !    apply scaling factor as specified in the inp-file

         do ipnt = 1, npoint
            points(ipnt,1) = inp_sclfac * points(ipnt,1)
            points(ipnt,2) = inp_sclfac * points(ipnt,2)
         enddo

         !  - modification 9: inversion

         call check_inversion(fname, npoint, points, 1, 2, is_wheel, x_profil, 0, flip_data, ierror)

         if (flip_data.eq.1) then
            call invert_points(npoint, points, is_wheel, x_profil)

            if (x_profil.ge.4) then
               call write_log(' after flip_data:')
               do ipnt = 1, npoint
                  write(bufout, '(i5,2(a,f9.3))') ipnt,': y=',points(ipnt,1),', z=',points(ipnt,2)
                  call write_log(1, bufout)
               enddo
            endif
         endif

         !  - first/last points should now be in appropriate order

         call check_overall_order(npoint, points, is_wheel, x_profil, ierror)

      endif ! ierror<>0

      ! call write_log('...done miniprof')

   end subroutine read_miniprof_profile

!------------------------------------------------------------------------------------------------------------

   subroutine delete_outside_boundbox(points, arrsiz, npoint, ncolpnt, bound_y_min, bound_y_max,        &
                bound_z_min, bound_z_max, ldebug)
!--purpose: modification 5: remove points from a profile that lie outside bounding box
      implicit none
!--subroutine arguments:
      integer,                      intent(in)    :: arrsiz, ncolpnt, ldebug
      integer,                      intent(inout) :: npoint
      real(kind=8),                 intent(in)    :: bound_y_min, bound_y_max, bound_z_min, bound_z_max
      real(kind=8), dimension(:,:), intent(inout) :: points(arrsiz,ncolpnt)
!--local variables:
      integer             :: ipnt, inew
      logical             :: lfilt_y, lfilt_z

      ! min >= max means no filtering

      lfilt_y = (bound_y_min.lt.bound_y_max)
      lfilt_z = (bound_z_min.lt.bound_z_max)

      if (ldebug.ge.2) then
         if (lfilt_y) then
            write(bufout,'(a,2(g11.4,a))') ' Removing points for which y not.in [',bound_y_min,         &
                        ',',bound_y_max,']'
            call write_log(1, bufout)
         endif
         if (lfilt_z) then
            write(bufout,'(a,2(g11.4,a))') ' Removing points for which z not.in [',bound_z_min,         &
                        ',',bound_z_max,']'
            call write_log(1, bufout)
         endif
      endif

      inew = 0
      do ipnt = 1, npoint
         if (lfilt_y .and. (points(ipnt,1).lt.bound_y_min .or. points(ipnt,1).gt.bound_y_max)) then
            if (ldebug.ge.3) then
               write(bufout,'(a,i5,a,g11.4,a)') ' ...discarding point',ipnt,': y=',points(ipnt,1),      &
                     ' outside [min,max]'
               call write_log(1, bufout)
            endif
         elseif (lfilt_z .and. (points(ipnt,2).lt.bound_z_min .or. points(ipnt,2).gt.bound_z_max)) then
            if (ldebug.ge.3) then
               write(bufout,'(a,i5,a,g11.4,a)') ' ...discarding point',ipnt,': z=',points(ipnt,2),      &
                     ' outside [min,max]'
               call write_log(1, bufout)
            endif
         else
            inew = inew + 1
            if (inew.lt.ipnt) points(inew,1:ncolpnt) = points(ipnt,1:ncolpnt)
         endif
      enddo
      npoint = inew
      if (ldebug.ge.2) then
         write(bufout,*) 'Boundbox: keeping',npoint,' profile points'
         call write_log(1, bufout)
      endif

   end subroutine delete_outside_boundbox

!------------------------------------------------------------------------------------------------------------

   subroutine filter_close_points(points, arrsiz, npoint, ncolpnt, point_dist_min, is_wheel, ldebug)
!--purpose: remove points from a profile that lie too close together
      implicit none
!--subroutine arguments:
      integer,                      intent(in)    :: arrsiz, ncolpnt, is_wheel, ldebug
      integer,                      intent(inout) :: npoint
      real(kind=8),                 intent(in)    :: point_dist_min
      real(kind=8), dimension(:,:), intent(inout) :: points(arrsiz,ncolpnt)
!--local variables:
      real(kind=8), dimension(:,:), allocatable   :: tmparr
      integer             :: ipnt, inew
      real(kind=8)        :: dist
      character(len=5)    :: nam_rw(0:1) = (/ 'rail ', 'wheel' /)

      !  - modification 1: point.dist.min - remove points too close together

      if (point_dist_min.ge.1d-12 .and. npoint.gt.0) then
         if (ldebug.ge.3) then
            write(bufout,'(a,g11.4,a)') ' Filtering points separated less than',point_dist_min,      &
                     ' from adjacent'
            call write_log(1, bufout)
         endif

         ! TODO: avoid copying

         allocate(tmparr(npoint,ncolpnt))
         inew = 1
         tmparr(inew,1:ncolpnt) = points(1,1:ncolpnt)

         do ipnt = 2, npoint
            dist = (tmparr(inew,1)-points(ipnt,1))**2 + (tmparr(inew,2)-points(ipnt,2))**2
            if (dist.ge.point_dist_min**2) then
               inew = inew + 1
               tmparr(inew,1:ncolpnt) = points(ipnt,1:ncolpnt)
            else
               if (ldebug.ge.3) then
                  write(bufout,'(a,i5,2(a,f10.6))') ' ...discarding point',ipnt,': y=',points(ipnt,1), &
                        ' too close to', tmparr(inew,1)
                  call write_log(1, bufout)
               endif
            endif
         enddo

         if (ldebug.ge.1 .and. npoint-inew.gt.0) then
            write(bufout,'(a,i5,3a)') ' INFO: deleting',npoint-inew,' points from ',                    &
                                        trim(nam_rw(is_wheel)), '-profile too close to adjacent...'
            call write_log(1, bufout)
         endif

         npoint = inew
         points(1:npoint,1:ncolpnt) = tmparr(1:npoint,1:ncolpnt)
         deallocate(tmparr)

         if (ldebug.ge.3) then
            write(bufout,'(a,g11.4,a,i6,a)') ' Point_dist_min: mindst=',point_dist_min,', keeping',     &
                        npoint,' profile points'
            call write_log(1, bufout)
         endif

      endif

   end subroutine filter_close_points

!------------------------------------------------------------------------------------------------------------

   subroutine remove_vertical_slopes(points, arrsiz, npoint, ncolpnt, ldebug)
!--purpose: shift profile points sideways to avoid vertical sections, to get a uni-valued function
      implicit none
!--subroutine arguments:
      integer,                      intent(in)    :: arrsiz, ncolpnt, ldebug
      integer,                      intent(in)    :: npoint
      real(kind=8), dimension(:,:), intent(inout) :: points(arrsiz,ncolpnt)
!--local variables:
      ! for dy>0 surface inclination atan2(dz,dy) is in [-120, 120] deg, negative on inside of rail
      real(kind=8), parameter :: alph_thrs    = 89.9d0 * pi/180d0
      integer             :: iymin, iymax  ! indices iy with lowest / highest y-values
      integer             :: izmin, inext
      real(kind=8)        :: dy, dz, sgn

      ! get left point 'iymin' and right point 'iymax' with lowest / highest y-values

      iymin  = idmin(npoint, points(:,1), 1)
      iymax  = idmax(npoint, points(:,1), 1)

      if (ldebug.ge.2) then
         write(bufout,'(2(a,f8.3,a,i5))') ' rail profile has min. y=', points(iymin,1),' at iymin=',    &
                iymin, ', max. y=', points(iymax,1),' at iymax=',iymax
         call write_log(1, bufout)
      endif

      ! get mid-point 'izmin': point with minimum z-value on the tread of the rail

      izmin = idmin(npoint, points(:,2), 1)

      if (ldebug.ge.2) then
         write(bufout,'(a,f8.3,a,i5)') ' rail profile has min. z=', points(izmin,2),' at izmin=',izmin
         call write_log(1, bufout)
      endif

      ! filter_dents: in range [izmin:iymax], place points in top surface with steep slopes at straight line

      do inext = izmin+1, iymax
         dy = points(inext,1) - points(inext-1,1)
         dz = points(inext,2) - points(inext-1,2)

         if (abs(atan2(dz,dy)).gt.alph_thrs) then
            sgn = sign(1d0, atan2(dz,dy))
            if (ldebug.ge.2) then
               write(bufout,'(2(a,i3),2(a,f10.5))') ' rail vertical slope at i=',inext,', y(',inext,    &
                     ')=', points(inext,1),' -->', points(inext-1,1) + sgn*dz/tan(alph_thrs)
               call write_log(1, bufout)
            endif
            points(inext,1) = points(inext-1,1) + sgn * dz / tan(alph_thrs)
         endif
      enddo

      ! filter_dents: in range [iymin:izmin], place points in top surface with steep slopes at straight line

      do inext = izmin-1, iymin, -1
         dy = points(inext+1,1) - points(inext,1)
         dz = points(inext+1,2) - points(inext,2)

       ! write(bufout,'(a,i4,7(a,f10.5))') 'inext=',inext,': lft=(',points(inext,1),',',points(inext,2), &
       !                                   '), rgt=(',points(inext+1,1),',',points(inext+1,2),'), dif=', &
       !                                   dy,',',dz,', alph=',abs(atan2(dz,dy))
       ! call write_log(1, bufout)

         if (abs(atan2(dz,dy)).gt.alph_thrs) then
            sgn = sign(1d0, atan2(dz,dy))
            if (ldebug.ge.2) then
               write(bufout,'(2(a,i3),2(a,f10.5))') ' rail vertical slope at i=',inext+1,', y(',inext,  &
                     ')=', points(inext,1),' -->', points(inext+1,1) - sgn*dz/tan(alph_thrs)
               call write_log(1, bufout)
            endif
            points(inext,1) = points(inext+1,1) - sgn*dz / tan(alph_thrs)
         endif
      enddo

   end subroutine remove_vertical_slopes

!------------------------------------------------------------------------------------------------------------

   subroutine filter_double_kinks(points, arrsiz, npoint, ncolpnt, is_wheel, zig_thrs, ldebug)
!--purpose: remove points from a profile that form a Z-pattern with double kink
      implicit none
!--subroutine arguments:
      integer,                      intent(in)    :: arrsiz, ncolpnt, is_wheel, ldebug
      integer,                      intent(inout) :: npoint
      real(kind=8),                 intent(in)    :: zig_thrs
      real(kind=8), dimension(:,:), intent(inout) :: points(arrsiz,ncolpnt)
!--local variables:
      real(kind=8), parameter    :: pi = 4d0*atan(1d0)
      real(kind=8), dimension(:),   allocatable   :: alpha
      real(kind=8), dimension(:,:), allocatable   :: tmparr
      integer             :: nseg, iseg, inew, icol, j
      real(kind=8)        :: dy, dz, dalph0, dalph1
      character(len=5)    :: nam_rw(0:1) = (/ 'rail ', 'wheel' /)

      if (npoint.le.3) return

      !  - modification 1': replace two points with Z-pattern by mean position

      if (ldebug.ge.3) then
         write(bufout,'(a,f6.3,a,f6.1,a)') ' Filtering points that make a zig-zag-pattern, thrs=',      &
                zig_thrs,' rad (', zig_thrs*180d0/pi, ' deg)'
         call write_log(1, bufout)
      endif

      nseg = npoint - 1
      allocate(alpha(nseg), tmparr(npoint,ncolpnt))

      ! compute inclinations of segments

      do iseg = 1, nseg
         dy  = points(iseg+1,1) - points(iseg,1)
         dz  = points(iseg+1,2) - points(iseg,2)
         alpha(iseg) = atan2(dz, dy)
      enddo

      inew = 1
      tmparr(1,1:ncolpnt) = points(1,1:ncolpnt)

      ! process internal segments  2..nseg-1

      iseg = 1
      do while(iseg.lt.nseg-1)

         iseg = iseg + 1

         ! compute change of angle to previous and following segments

         dalph0 = alpha(iseg) - alpha(iseg-1)
         if (dalph0.lt.-pi) dalph0 = dalph0 + 2d0 * pi
         if (dalph0.gt. pi) dalph0 = dalph0 - 2d0 * pi

         dalph1 = alpha(iseg+1) - alpha(iseg)
         if (dalph1.lt.-pi) dalph1 = dalph1 + 2d0 * pi
         if (dalph1.gt. pi) dalph1 = dalph1 - 2d0 * pi

         ! detect double kinks

         if (abs(dalph0).gt.zig_thrs .and. abs(dalph1).gt.zig_thrs) then

            if (ldebug.ge.1) then
               write(bufout,'(3a,2(f7.2,a),i5,a,3f7.1,a)') ' Zig-zag in ',trim(nam_rw(is_wheel)),       &
                        ' profile at (', points(iseg+1,1),',', points(iseg+1,2),') (ip=', iseg+1,       &
                        '): successive angles', (alpha(iseg+j)*180d0/pi, j=-1,1),' deg'
               call write_log(1, bufout)
            endif
            if (ldebug.ge.3) then
               write(bufout,'(a,i5,2(a,f8.3))') ' ...discarding point',iseg+1,': y,z=',                 &
                     points(iseg+1,1), ',', points(iseg+1,2)
               call write_log(1, bufout)
            endif

            ! replace start-point of iseg by mean(start,end)

            inew = inew + 1
            do icol = 1, ncolpnt
               tmparr(inew,icol) = (points(iseg,icol) + points(iseg+1,icol)) / 2d0
            enddo

            ! skip next segment

            iseg = iseg + 1

         else

            ! no double kink: copy start-point of iseg

            inew = inew + 1
            tmparr(inew,1:ncolpnt) = points(iseg,1:ncolpnt)

         endif

      enddo

      inew = inew + 1
      tmparr(inew,1:ncolpnt) = points(npoint-1,1:ncolpnt)
      inew = inew + 1
      tmparr(inew,1:ncolpnt) = points(npoint  ,1:ncolpnt)

      if (ldebug.ge.1 .and. npoint-inew.gt.0) then
         write(bufout,'(a,i5,3a)') ' WARNING: deleting',npoint-inew,' points from ',                    &
                                     trim(nam_rw(is_wheel)), '-profile for zig-zag pattern...'
         call write_log(1, bufout)
      endif

      npoint = inew
      points(1:npoint,1:ncolpnt) = tmparr(1:npoint,1:ncolpnt)
      deallocate(tmparr)
      deallocate(alpha)

   end subroutine filter_double_kinks

!------------------------------------------------------------------------------------------------------------

   subroutine check_z_pos_down(fname, npoint, points, is_wheel, ldebug, mirror_z)
!--purpose: check that z-values are defined positive downwards, mirror z if necessary
      implicit none
!--subroutine arguments:
      character*(*)                               :: fname
      integer,                      intent(in)    :: npoint, is_wheel, ldebug
      real(kind=8), dimension(:,:), intent(inout) :: points
      integer,                      intent(out)   :: mirror_z
!--local variables:
      integer             :: izmin, izmax, ixm, i, chng_sta, chng_end
      real(kind=8)        :: frac_min, frac_max, z_strgt, dz_tot, dz

      if (ldebug.ge.3) then
         write(bufout,'(a,i6,a)') ' Check_z_pos_down:', npoint,' points'
         call write_log(1, bufout)
      endif

      izmin = idmin(npoint, points(:,2), 1)
      izmax = idmax(npoint, points(:,2), 1)
      frac_min = (1d0*izmin) / npoint
      frac_max = (1d0*izmax) / npoint

      if (ldebug.ge.2) then
         write(bufout,'(2(a,f7.3,a,i6), a,f6.1,a)') ' Profile has zmin=',points(izmin,2),' at i=',      &
                   izmin, ', zmax=',points(izmax,2),' at i=',izmax,', fmin at', 100d0*frac_min,' %'
         call write_log(1, bufout)
      endif

      mirror_z = 0

      if (is_wheel.eq.1) then

         ! wheel: maximum z-value should occur in the interior, at the flange of the wheel
         !        minimum z-value should occur at either side, at the field side of the wheel
         !        note that minimum z-value may lie in the interior for hollow worn wheels
         ! mirror z when minimum z-value does not occur in starting or trailing 5% of the profile points

         if (frac_min.lt.0.05 .or. frac_min.gt.0.95) then
            mirror_z = 0
         elseif (frac_max.lt.0.05 .or. frac_max.gt.0.95) then
            mirror_z = 1
         else
            ixm     = nint( npoint * (frac_min+frac_max)/2d0 )
            z_strgt = (points(izmin,2)+points(izmax,2)) / 2d0
            if (points(ixm,2).gt.z_strgt) then
               mirror_z = 1  ! mirror when profile z > straight line
            else
               mirror_z = 0
            endif
            call write_log(' hollow worn wheel?')
            call write_log(trim(fname))
            write(bufout,'(2(a,f7.3,a,f5.1),a)') ' profile has zmin=',points(izmin,2),' at f=',         &
                   100d0*frac_min, ' %, zmax=',points(izmax,2),' at f=',100d0*frac_max,' %'
            call write_log(1, bufout)
         endif

         if (ldebug.ge.1 .and. mirror_z.eq.1) then
            call write_log(' INFO: mirroring z_wheel = -z_wheel to get z positive downwards...')
         endif

      elseif (is_wheel.ne.1) then

         ! rail: z-values should be decreasing from the sides towards the middle
         !       using |dz| > 0.001 * dz_tot to distinguish `essential changes' from `perturbations' 

         dz_tot = max(1d-5, points(izmax,2) - points(izmin,2))

         ! look for first i with |dz| > 0.1% of total dz,
         !    set chng_sta = sign(dz)

         chng_sta = 0   ! -1: z decreasing (ok); 0: undecided; 1: z increasing (needs mirroring)
         i = 1
         do while(chng_sta.eq.0 .and. i.lt.nint(0.1d0*npoint))
            i  = i + 1
            dz = points(i,2) - points(1,2)      ! interior - boundary
            if (abs(dz).gt.0.001d0*dz_tot) then
               chng_sta = sign(1d0, dz)
               !write(bufout,'(a,i4,2(a,f8.3),a,i2)') ' i=',i,': dz=',dz,' above 0.1% of dz_tot=',       &
               !         dz_tot,', chng_sta=',chng_sta
               !call write_log(1, bufout)
            endif
         enddo

         ! check last 10% of points: 
         !     look for highest i with change at least 0.1% of total dz,
         !       then check if interior z < boundary z

         chng_end = 0   ! -1: z decreasing; 0: undecided; 1: z increasing
         i = npoint
         do while(chng_end.eq.0 .and. i.gt.nint(0.9d0*npoint))
            i  = i - 1
            dz = points(i,2) - points(npoint,2) ! interior - boundary
            if (abs(dz).gt.0.001d0*dz_tot) then
               chng_end = sign(1d0, dz)
               !write(bufout,'(a,i4,2(a,f8.3),a,i2)') ' i=',i,': dz=',dz,' above 0.1% of dz_tot=',       &
               !         dz_tot,', chng_end=',chng_end
               !call write_log(1, bufout)
            endif
         enddo

         ! set mirror_z based on chng_sta and chng_end

         if (chng_sta.eq.chng_end .or. chng_end.eq.0) then
            mirror_z = chng_sta
         elseif (chng_sta.eq.0) then
            mirror_z = chng_end
         else
            mirror_z = -1       ! opposite tendencies
         endif

         ! fallback: mirror if minimum z-value occurs in starting or trailing 10% of the profile points

         if (mirror_z.eq.0) then
            mirror_z = -1
            if (frac_min.lt.0.1 .or. frac_min.gt.0.9) mirror_z = 1
         endif

         if (mirror_z.eq.1 .and. ldebug.ge.1) then
            call write_log(' INFO: mirroring z_rail = -z_rail to get z positive downwards...')
         endif

     endif

  end subroutine check_z_pos_down

!------------------------------------------------------------------------------------------------------------

   subroutine check_inversion(fname, npoint, points, icoly, icolz, is_wheel, ldebug, show_errors,       &
                              flip_data, ierror)
!--purpose: check that y-values of rail are in increasing order at top of rail, 
!                      wheel-y are in decreasing order at bottom of wheel profile
!                      show_errors=0: determine if inversion is needed; 1: set ierror, print message
      implicit none
!--subroutine arguments:
      character*(*)                               :: fname
      integer,                      intent(in)    :: npoint, icoly, icolz, is_wheel, show_errors, ldebug
      real(kind=8), dimension(:,:), intent(inout) :: points
      integer,                      intent(out)   :: flip_data, ierror
!--local variables:
      integer           :: i1, i2, i3, izmin, izmax
      real(kind=8)      :: dy1, dy2, dz1, dz2, alph1, alph2, dalph

      flip_data = 0
      ierror    = 0

      ! check that inside (material) is at the right hand side of the profile curve

      if (is_wheel.eq.1) then

         ! wheel: using highest z-value, typically at bottom of flange

         izmax = idmax(npoint, points(:,icolz), 1)

         ! if izmax is the first or the last point (e.g. pure conical wheel),
         !    the step to the neighbour should have the correct sign.

         if (izmax.le.1) then

            if (.not.(points(izmax+1,icoly).lt.points(izmax,icoly))) flip_data = 1

         elseif (izmax.ge.npoint) then

            if (.not.(points(izmax,icoly).lt.points(izmax-1,icoly))) flip_data = 1

         else

            dy1 = points(izmax,icoly) - points(izmax-1,icoly)
            dy2 = points(izmax+1,icoly) - points(izmax,icoly)

            if (ldebug.ge.2) then
               write(bufout,'(a,i4,2(a,f14.6))') ' wheel: max(z) at i=',izmax,', dy1=',dy1,', dy2=',dy2
               call write_log(1, bufout)
            endif

            ! if the steps left/right have the same sign, then both must be negative

            if (dy1*dy2.gt.0d0) then

               if (.not.(dy1.lt.0d0)) flip_data = 1

            else

               ! in case of a sharp corner, '<' or '>', the profile should turn to the right.

               dz1 = points(izmax,icolz) - points(izmax-1,icolz)
               dz2 = points(izmax+1,icolz) - points(izmax,icolz)
               alph1 = atan2(dz1, dy1)
               alph2 = atan2(dz2, dy2)
               dalph = alph2 - alph1
               if (dalph.le.-pi) dalph = dalph + 2d0*pi
               if (dalph.gt. pi) dalph = dalph - 2d0*pi

               if (.not.(dalph.ge.0d0)) flip_data = 1

               if (ldebug.ge.2) then
                  write(bufout,'(a,i4,3(a,f6.1))') ' wheel: max(z) at i=',izmax,', alph1=',alph1*180d0/pi, &
                                ', alph2=',alph2*180d0/pi,', dalph=',dalph*180d0/pi
                  call write_log(1, bufout)
               endif

            endif ! dy1*dy2 > 0
         endif ! izmax=1 or npoint

         ! if (mode==checking for errors): set error-code and print message

         if (flip_data.ge.1 .and. show_errors.ge.1) then
            ierror = 221
            call write_log(' ERROR 221: wheel profile y-values should be decreasing.')
            call write_log('            file "' // trim(fname) // '"')
            i1 = max(1,izmax-1)
            i2 = izmax
            i3 = min(npoint,izmax+1)
            write(bufout,221) '            y(',i1,')=',points(i1,icoly),', y(',i2,')=', points(i2,icoly), &
                        ', y(',i3,')=', points(i3,icoly)
 221        format(3(a,i4,a,f14.6))
            call write_log(1, bufout)
         endif

      else

         ! rail: using lowest z-value, typically at top of rail

         izmin = idmin(npoint, points(:,icolz), 1)

         if (izmin.ge.npoint) then
            i1 = izmin - 1
            i2 = izmin
         else
            i1 = izmin
            i2 = izmin + 1
         endif

         if (ldebug.ge.2) then
            write(bufout,'(a,i4,2(a,i4,a,f14.6))') ' rail: min(z) at i=',izmin,', y(',i1,')=',          &
                points(i1,icoly), ', y(',i2,')=',points(i2,icoly)
            call write_log(1, bufout)
         endif

         if (points(i2,icoly).lt.points(i1,icoly)) flip_data = 1

         ! if (mode==checking for errors): set error-code and print message

         if (flip_data.ge.1 .and. show_errors.ge.1) then
            ierror = 222
            call write_log(' ERROR 222: rail profile y-values should be increasing.')
            call write_log('            file "' // trim(fname) // '"')
            write(bufout,221) '            y(',i1,')=',points(i1,icoly),', y(',i2,')=', points(i2,icoly)
            call write_log(1, bufout)
         endif

      endif

   end subroutine check_inversion

!------------------------------------------------------------------------------------------------------------

   subroutine invert_points(npoint, points, is_wheel, ldebug)
!--purpose: flip points-array: [1:n] --> [n:-1:1]
      implicit none
!--subroutine arguments:
      integer,                      intent(in)    :: npoint, is_wheel, ldebug
      real(kind=8), dimension(:,:), intent(inout) :: points  ! size: (arrsiz,ncol), with arrsiz >= npoint
!--local variables:
      integer      :: ncol, icol, nhlf, ipnt, iupp
      real(kind=8) :: tmp

      ncol = size(points,2)
      nhlf = npoint / 2

      if (ldebug.ge.2) then
         if (is_wheel.eq.1) then
            call write_log(' INFO: reordering wheel-profile to make y-values decreasing...')
         else
            call write_log(' INFO: reordering rail-profile to make y-values increasing...')
         endif
      endif

      do icol = 1, ncol
         do ipnt = 1, nhlf

            ! compute companion of ipnt in upper half

            iupp = npoint + 1 - ipnt

            ! swap values

            tmp               = points(ipnt,icol)
            points(ipnt,icol) = points(iupp,icol)
            points(iupp,icol) = tmp
         enddo
      enddo

   end subroutine invert_points

!------------------------------------------------------------------------------------------------------------

   subroutine check_overall_order(npoint, points, is_wheel, ldebug, ierror)
!--purpose: check that y-values of rail start at the left, go to the right
!                      wheel-y start right, go left
      implicit none
!--subroutine arguments:
      integer,                      intent(in)    :: npoint, is_wheel, ldebug
      real(kind=8), dimension(:,:), intent(inout) :: points
      integer,                      intent(out)   :: ierror
!--local variables:
      integer              :: j, ir

      ierror = 0

      if (is_wheel.eq.1) then
         ! wheel: error if y(1) <= y(end)
         if (points(1,1).le.points(npoint,1)) ierror = 211
      else
         ! rail: error if y(1) >= y(end)
         if (points(1,1).ge.points(npoint,1)) ierror = 212
      endif

      if (ierror.ne.0) then

         if (is_wheel.eq.1 .and. ldebug.ge.1) then
            call write_log(' ERROR: y-coordinates of wheel profile should be descending.')
         elseif (ldebug.ge.1) then
            call write_log(' ERROR: y-coordinates of rail profile should be ascending.')
         endif

         if (ldebug.ge.1) then
            do j = 1, 6
               ir = j
               if (j.ge.4) ir = npoint - 6 + j
               write(bufout,'(a,i5,2(a,f9.3))') ' i=',ir,': y(i)=',points(ir,1),', z(i)=',points(ir,2)
               call write_log(1, bufout)
               if (j.eq.3) call write_log('     ...')
            enddo
         endif
      endif

   end subroutine check_overall_order

!------------------------------------------------------------------------------------------------------------

   subroutine filter_receding_order(fname, npoint, points, is_wheel, ldebug, ierror)
!--purpose: modify points-array: remove points where rail y-values are non-increasing (wheel: non-decreasing)
      implicit none
!--subroutine arguments:
      character*(*)                               :: fname
      integer,                      intent(in)    :: is_wheel, ldebug
      integer,                      intent(inout) :: npoint
      real(kind=8), dimension(:,:), intent(inout) :: points   ! size (arrsiz,ncol), arrsiz>=npoint
      integer,                      intent(out)   :: ierror
!--local variables:
      integer             :: ncol, ipnt, inew, izmin, i_run_min, i_run_max,                             &
                             idel_sta, idel_mid, idel_end, idel_tot
      integer             :: keep(npoint)
      real(kind=8)        :: dy_thresh, y_run_min, y_run_max, sgn
      character(len=5)    :: namprf

      if (ldebug.ge.3) call write_log('filter_receding_order')

      ierror = 0
      ncol = size(points,2)

      ! the difference |y(1)-y(end)| is used as order of magnitude for the data
      ! for |y| around 100 mm, points will be rejected when less than 1d-4 mm away of adjacent points

      dy_thresh = 1d-6 * abs(points(1,1) - points(npoint,1))

      if (ldebug.ge.2) then
         write(bufout,'(3a,i5,a,f9.6)') ' input profile "',trim(fname),'" has',npoint,                  &
                ' points, dy_thresh=', dy_thresh
         call write_log(1, bufout)
      endif
      if (ldebug.ge.4) then
         do ipnt = 1, npoint
            write(bufout, '(i5,2(a,f9.3))') ipnt,': y=',points(ipnt,1),', z=',points(ipnt,2)
            call write_log(1, bufout)
         enddo
      endif

      ! rail: sgn = 1, wheel: sgn = -1, applying rail algorithm on mirrored y-data

      if (is_wheel.eq.1) then
         ! wheel: get index of point with maximum z-value on the flange
         namprf = 'wheel'
         sgn    = -1d0
         izmin  = idmax(npoint, points(:,2), 1)
         if (ldebug.ge.2) then
            write(bufout,'(a,f8.3,a,i5)') ' wheel profile has max. z=', points(izmin,2), ' at izmin=',izmin
            call write_log(1, bufout)
         endif
      else
         ! rail: get mid-point 'izmin': point with minimum z-value on the tread of the rail
         namprf = 'rail'
         sgn    =  1d0
         izmin  = idmin(npoint, points(:,2), 1)
         if (ldebug.ge.2) then
            write(bufout,'(a,f8.3,a,i5)') ' rail profile has min. z=', points(izmin,2), ' at izmin=',izmin
            call write_log(1, bufout)
         endif
      endif

      ! fill the mask-array "keep" to indicate which points should be maintained

      keep(1:npoint) = 1

      ! Note: names & comments describe processing of the rail profile
      ! loop over all points to the right of izmin, track running maximum of y-values

      y_run_max = sgn * points(izmin,1)
      i_run_max = izmin

      do ipnt = izmin+1, npoint

         if (sgn*points(ipnt,1).lt.y_run_max-dy_thresh) then

            ! discard new point "i" if y-values are decreasing

            keep(ipnt) = 0

            if (ldebug.ge.3 .and. is_wheel.eq.1) then
               write(bufout,'(a,f12.4,a,i5,a)') ' wheel y(i) =', points(ipnt,1),' >  y(i-1) at i=',  &
                     ipnt,', discarding i'
               call write_log(1, bufout)
            elseif (ldebug.ge.3) then
               write(bufout,'(a,f12.4,a,i5,a)') ' rail y(i) =', points(ipnt,1),' <  y(i-1) at i=',  &
                     ipnt,', discarding i'
               call write_log(1, bufout)
            endif

         elseif (sgn*points(ipnt,1).lt.y_run_max+dy_thresh) then

            ! discard previous point "i_run_max" if y-values are (nearly) the same

            keep(i_run_max) = 0
            i_run_max = ipnt

            if (ldebug.ge.3) then
               write(bufout,'(1x,a,2(a,i5))') namprf, ' y(i) == y(i-1) at i=',ipnt,', discarding', i_run_max
               call write_log(1, bufout)
            endif

         endif

         ! update running maximum of y-values

         if (sgn*points(ipnt,1).gt.y_run_max) then
            y_run_max = sgn*points(ipnt,1)
            i_run_max = ipnt
         endif

      enddo

      ! loop over all points to the left of izmin, track running minimum of y-values

      y_run_min = sgn * points(izmin,1)
      i_run_min = izmin

      do ipnt = izmin-1, 1, -1

         if (sgn*points(ipnt,1).gt.y_run_min+dy_thresh) then

            ! discard new point "i" if y-values are increasing

            keep(ipnt) = 0

            if (ldebug.ge.3 .and. is_wheel.eq.1) then
               write(bufout,'(a,f12.4,a,i5,a)') ' wheel y(i) =', points(ipnt,1),' <  y(i+1) at i=',  &
                     ipnt,', discarding i'
               call write_log(1, bufout)
            elseif (ldebug.ge.3) then
               write(bufout,'(a,f12.4,a,i5,a)') ' rail y(i) =', points(ipnt,1),' >  y(i+1) at i=',  &
                     ipnt,', discarding i'
               call write_log(1, bufout)
            endif

         elseif (sgn*points(ipnt,1).gt.y_run_min-dy_thresh) then

            ! discard previous point "i_run_min" if y-values are (nearly) the same

            keep(i_run_min) = 0
            i_run_min = ipnt

            if (ldebug.ge.3) then
               write(bufout,'(1x,a,2(a,i5))') namprf, ' y(i) == y(i+1) at i=',ipnt,', discarding', i_run_min
               call write_log(1, bufout)
            endif

         endif

         ! update running minimum of y-values

         if (sgn*points(ipnt,1).lt.y_run_min) then
            y_run_min = sgn*points(ipnt,1)
            i_run_min = ipnt
         endif

      enddo

      ! count number of points discarded at start of profile, interior to profile and at end of profile

      idel_sta = 0
      do while(idel_sta.lt.npoint .and. keep(idel_sta+1).eq.0)
         idel_sta = idel_sta + 1
      enddo

      idel_end = 0
      do while(idel_end.lt.npoint .and. keep(npoint-idel_end).eq.0)
         idel_end = idel_end + 1
      enddo

      idel_mid = 0
      do ipnt = idel_sta+1, npoint-idel_end-1
         if (keep(ipnt).eq.0) idel_mid = idel_mid + 1
      enddo

      idel_tot = 0
      do ipnt = 1, npoint
         if (keep(ipnt).eq.0) idel_tot = idel_tot + 1
      enddo

      ! Warn the user about the numbers of points that are deleted

      if (ldebug.ge.1 .and. idel_tot.ge.1) then
         if (is_wheel.eq.1) then
            write(bufout,'(a,i4,2a)') ' WARNING:',idel_tot,' points deleted from the wheel profile ',   &
                        'where y(i+1) >= y(i)'
         else
            write(bufout,'(a,i4,2a)') ' WARNING:',idel_tot,' points deleted from rail profile where',   &
                        ' y(i+1) <= y(i)'
         endif
         call write_log(1, bufout)
      endif

      if (ldebug.ge.2) then
         if (idel_sta.ge.1) then
            write(bufout,'(9x,i4,a)') idel_sta,' points are deleted from the start of the profile'
            call write_log(1, bufout)
         endif
         if (idel_mid.ge.1) then
            write(bufout,'(9x,i4,a)') idel_mid,' points are deleted from the interior of the profile'
            call write_log(1, bufout)
         endif
         if (idel_end.ge.1) then
            write(bufout,'(9x,i4,a)') idel_end,' points are deleted from the end of the profile'
            call write_log(1, bufout)
         endif
      endif

      if (idel_sta+idel_mid+idel_end.ne.idel_tot) then
         call write_log(' Internal Error: deleted points dont add up')
      endif

      ! copy rows to be kept from original position ipnt to new position inew

      inew = 0
      do ipnt = 1, npoint
         if (keep(ipnt).ge.1) then
            inew = inew + 1
            if (ipnt.gt.inew) points(inew,1:ncol) = points(ipnt,1:ncol)

            if (ldebug.ge.5) then
               write(bufout,'(a,i4,2(a,f8.3),a)') ' new i=',inew,' = (', points(ipnt,1),',',    &
                         points(ipnt,2),')'
               call write_log(1, bufout)
            endif
         endif
      enddo

      ! return the number of points kept

      npoint = inew

   end subroutine filter_receding_order

!------------------------------------------------------------------------------------------------------------

   subroutine profile_find_kinks(npnt, y, z, is_wheel, ang_thrs_high, ang_thrs_low, dst_max, scale_z,   &
                        nkink, ikinks, ldebug, ierror, s, ds_thrs)
!--function: determine kinks in planar curve {y(s), z(s)}: discontinuous 1st derivative
      implicit none
!--subroutine arguments:
      integer,      intent(in)  :: npnt, is_wheel, ldebug
      real(kind=8), intent(in)  :: ang_thrs_high, ang_thrs_low, dst_max, scale_z
      real(kind=8), intent(in)  :: y(npnt), z(npnt)
      integer,      intent(out) :: nkink
      integer,      intent(out) :: ikinks(:)
      integer,      intent(out) :: ierror
      real(kind=8), intent(in),  optional :: ds_thrs
      real(kind=8), intent(in),  optional :: s(npnt)
!--local variables:
      real(kind=8), parameter   :: pi = 4d0*atan(1d0)
      logical      :: has_s_thrs, is_kink
      integer      :: max_kinks, nignored, ipnt, k
      real(kind=8) :: ds0, ds1, dy0, dy1, dz0, dz1, dst
      real(kind=8), dimension(:), allocatable :: dalph
      character(len=13)    :: nam_rw

      ierror = 0
      max_kinks = size(ikinks)
      ikinks(1:max_kinks) = 0

      has_s_thrs = (present(s) .and. present(ds_thrs))
   
      allocate(dalph(npnt))

      if (ldebug.ge.2) then
         write(bufout,'(3(a,f6.1),a)') ' find_kinks with thrs_high=', ang_thrs_high*180d0/pi,           &
                ' deg, thrs_low=', ang_thrs_low*180d0/pi,' deg, dst_max=',dst_max,' mm'
         call write_log(1, bufout)
      endif

      ! compute difference of surface inclinations of any two consecutive segments

      ! orientation: rail processed from left to right, moving along profile with material (inside) at
      !              right hand. -90deg vertical on inner gauge face, 0deg == horizontal to the right,
      !                           90deg vertical on outer gauge face.

      do ipnt = 2, npnt-1

         dy0 =  y(ipnt) - y(ipnt-1)
         dz0 = (z(ipnt) - z(ipnt-1)) * scale_z  ! for the contact locus, z <-- x with larger fluctuation
         dy1 =  y(ipnt+1) - y(ipnt)
         dz1 = (z(ipnt+1) - z(ipnt)) * scale_z

         dalph(ipnt) = atan2(dz1, dy1) - atan2(dz0, dy0)
         if (dalph(ipnt).lt.-pi) dalph(ipnt) = dalph(ipnt) + 2d0 * pi
         if (dalph(ipnt).gt. pi) dalph(ipnt) = dalph(ipnt) - 2d0 * pi
      enddo

      ! store start of first section

      nignored  = 0
      nkink     = 1
      ikinks(1) = 1

      do ipnt = 2, npnt-1

         is_kink = .false.

         if (has_s_thrs) then

            ! mark as kink if a segment length is larger than ds_thrs  -- contact locus, z = -999

            ds0 = s(ipnt) - s(ipnt-1)
            ds1 = s(ipnt+1) - s(ipnt)

            if (ds0.gt.ds_thrs .or. ds1.gt.ds_thrs) is_kink = .true.

            if (ldebug.ge.1 .and. is_kink) then
               if (is_wheel.lt.0) then
                  nam_rw = 'contact locus'
               elseif (is_wheel.eq.0) then
                  nam_rw = 'rail profile '
               else
                  nam_rw = 'wheel profile'
               endif
               write(bufout,'(2a,2(a,f7.2),a,i5,2(a,f7.1),a)') ' Kink in ',trim(nam_rw),' at (', &
                   y(ipnt),',',z(ipnt), ') (ip=',ipnt,'): successive steps ', ds0,',',ds1,' mm'
               call write_log(1, bufout)
            endif
         endif

         ! check for kink if |d.angle(i)| > high threshold

         if (.not.is_kink .and. abs(dalph(ipnt)).gt.ang_thrs_high) then

            ! require that |d.angle(j)| <= low threshold for j in window [-dst_max,dst_max]

            if (ldebug.ge.3) then
               write(bufout,'(a,i4,2(a,f6.2),a)') ' possible kink at ipnt=',ipnt,', (',y(ipnt),',',z(ipnt),')'
               call write_log(1, bufout)
            endif

            is_kink = .true.

            ! check neighbours ipnt+k to the left

            k = 0
            dst = 0d0
            do while(is_kink .and. ipnt+k.gt.2 .and. dst.lt.dst_max)
               k   = k - 1
               dst = sqrt( (y(ipnt)-y(ipnt+k))**2 + (z(ipnt)-z(ipnt+k))**2 )
               if (dst.lt.dst_max .and. abs(dalph(ipnt+k)).gt.ang_thrs_low) is_kink = .false.

               if (ldebug.ge.3) then
                  write(bufout,'(a,i3,3(a,f6.2),a,f7.1)') '  - ngb k=',k,': (',y(ipnt+k),',',z(ipnt+k), &
                        '), dst=', dst,', dalph=',dalph(ipnt+k)*180d0/pi
                  call write_log(1, bufout)
               endif
            enddo

            ! check neighbours ipnt+k to the right

            k = 0
            dst = 0d0
            do while(is_kink .and. ipnt+k.lt.npnt-1 .and. dst.lt.dst_max)
               k   = k + 1
               dst = sqrt( (y(ipnt+k)-y(ipnt))**2 + (z(ipnt+k)-z(ipnt))**2 )
               if (dst.lt.dst_max .and. abs(dalph(ipnt+k)).gt.ang_thrs_low) is_kink = .false.

               if (ldebug.ge.3) then
                  write(bufout,'(a,i3,3(a,f6.2),a,f7.1)') '  - ngb k=',k,': (',y(ipnt+k),',',z(ipnt+k), &
                        '), dst=', dst,', dalph=',dalph(ipnt+k)*180d0/pi
                  call write_log(1, bufout)
               endif
            enddo
         endif
   
         if (is_kink) then
            if (ldebug.ge.1) then
               if (is_wheel.lt.0) then
                  nam_rw = 'contact locus'
               elseif (is_wheel.eq.0) then
                  nam_rw = 'rail profile '
               else
                  nam_rw = 'wheel profile'
               endif
               dy0 = y(ipnt) - y(ipnt-1)
               dz0 = z(ipnt) - z(ipnt-1)
               write(bufout,'(2a,2(a,f7.2),a,i5,2(a,f7.1),a)') ' Kink in ',trim(nam_rw),' at (', &
                   y(ipnt),',',z(ipnt), ') (ip=',ipnt,'): successive angles',                           &
                   atan2(dz0,dy0)*180d0/pi,',', (atan2(dz0,dy0)+dalph(ipnt))*180d0/pi,' deg'
               call write_log(1, bufout)
            endif
            if (nkink.ge.max_kinks-2) then
               nignored = nignored + 1
               if (nignored.le.5) then
                  write(bufout,'(a,i5,a,i3,a)') ' Warning: ignoring kink at ip=',ipnt,': max.',         &
                              max_kinks, ' kinks supported.'
                  call write_log(1, bufout)
               endif
            else
               if (ipnt-ikinks(nkink).eq.2) then
                  ! avoid kinks with just one point between them, for B-spline approximation
                  nkink = nkink + 1
                  ikinks(nkink) = ipnt-1
                  if (ldebug.ge.1) then
                     write(bufout,'(2(a,i5),a)') ' This kink lies too close to previous kink at ip=',   &
                        ipnt-2, ': adding additional break at ip=',ipnt-1,'.'
                     call write_log(1, bufout)
                  endif
               endif
               nkink = nkink + 1
               ikinks(nkink) = ipnt
            endif
         endif
      enddo
   
      ! add end of last section

      if (nkink.lt.max_kinks) then
         nkink         = nkink + 1
         ikinks(nkink) = npnt
      endif

      if (nignored.gt.1) then
         write(bufout,'(2(a,i3),a)') ' Warning: max.',max_kinks,' kinks supported,',nignored,' are ignored.'
         call write_log(1, bufout)
      endif
   
      deallocate(dalph)

   end subroutine profile_find_kinks

!------------------------------------------------------------------------------------------------------------

   subroutine profile_check_bowtie(npnt, y, z, is_wheel, ldebug, ierror)
!--function: determine the number of intersections between different segments
      implicit none
!--subroutine arguments:
      integer,      intent(in)  :: npnt, is_wheel, ldebug
      real(kind=8), intent(in)  :: y(npnt), z(npnt)
      integer,      intent(out) :: ierror
!--local variables:
      integer,      parameter   :: max_nerror = 5
      real(kind=8), parameter   :: tiny = 1d-9
      logical           :: has_isec
      integer           :: num_isec, nseg, iseg, jseg
      real(kind=8)      :: y1a, y1b, y2a, y2b, z1a, z1b, z2a, z2b
      character(len= 5) :: nam_rw(0:1) = (/ 'rail ', 'wheel' /)

      ierror   = 0
      num_isec = 0
      nseg     = npnt - 1
   
      ! brute-force all-to-all comparison

      do iseg = 1, nseg-1
         y1a = y(iseg)
         z1a = z(iseg)
         y1b = y(iseg+1)
         z1b = z(iseg+1)

         ! check against all following segments except direct neighbour

         do jseg = iseg+2, nseg
            y2a = y(jseg)
            z2a = z(jseg)
            y2b = y(jseg+1)
            z2b = z(jseg+1)

            has_isec = segs_have_intersection(iseg, y1a, z1a, y1b, z1b, jseg, y2a, z2a, y2b, z2b,       &
                                                                                        tiny, ldebug)
            if (has_isec) then

               num_isec = num_isec + 1

               if (ldebug.ge.1 .and. num_isec.le.max_nerror) then
                  write(bufout,'(3a,2(i4,a))') ' ERROR: ',trim(nam_rw(is_wheel)),' profile segments',   &
                        iseg, ' and',jseg,' are intersecting.'
                  call write_log(1, bufout)
                  write(bufout,'(2(a,i4,4(a,f8.3),a,:,/))')                                             &
                        '        seg',iseg,': (',y1a,',',z1a,')--(',y1b,',',z1b,')',                    &
                        '        seg',jseg,': (',y2a,',',z2a,')--(',y2b,',',z2b,')'
                  call write_log(2, bufout)
               endif

            endif
         enddo
      enddo

      if (num_isec.ge.1) then

         ierror = 323
         if (ldebug.ge.0) then
            write(bufout,'(3a,i4,a)') ' ERROR: ',trim(nam_rw(is_wheel)),' profile has',num_isec,        &
                   ' intersections between different segments.'
            call write_log(1, bufout)
         endif

      elseif (ldebug.ge.3) then

         write(bufout,'(3a)') ' check_bowtie: ', trim(nam_rw(is_wheel)),                                &
                ' has no intersections between different segments.'
         call write_log(1, bufout)

      endif

   end subroutine profile_check_bowtie

!------------------------------------------------------------------------------------------------------------

   function segs_have_intersection(iseg, y1a, z1a, y1b, z1b, jseg, y2a, z2a, y2b, z2b, tiny, ldebug)
!--function: check whether line segments (1a)-(1b) and (2a)-(2b) are intersecting each other
!            tiny == safety margin for 'close encounters'.
!                    Positive values mean that more intersections are found, 
!                    Negative tiny could be used to exclude intersections just at the end-points.
      implicit none
!--result value
      logical      :: segs_have_intersection
!--subroutine arguments:
      integer,      intent(in)  :: iseg, jseg, ldebug
      real(kind=8), intent(in)  :: y1a, y1b, y2a, y2b, z1a, z1b, z2a, z2b, tiny
!--local variables:
      real(kind=8) :: len1, nvec1(2), tvec1(2), n2a, n2b, t2a, t2b, ti

      associate(has_isec => segs_have_intersection)

      if (max(y1a,y1b)+tiny.lt.min(y2a,y2b) .or. max(y2a,y2b)+tiny.lt.min(y1a,y1b)) then

         has_isec = .false.
         if (ldebug.ge.4) then
            write(bufout,'(2(a,i4),a)') ' segments',iseg,',',jseg,': no overlap in y'
            call write_log(1, bufout)
         endif

      elseif (max(z1a,z1b)+tiny.lt.min(z2a,z2b) .or. max(z2a,z2b)+tiny.lt.min(z1a,z1b)) then

         has_isec = .false.
         if (ldebug.ge.4) then
            write(bufout,'(2(a,i4),a)') ' segments',iseg,',',jseg,': no overlap in z'
            call write_log(1, bufout)
         endif

      else
         ! original, global coordinates

         ! tangent vector tvec1 along segment 1
         tvec1(1) = y1b - y1a
         tvec1(2) = z1b - z1a
         len1  = sqrt(tvec1(1)**2 + tvec1(2)**2)
         if (len1<1d-20) then
            len1 = 1d0
            if (ldebug.ge.2) then
               write(bufout,'(a,i4,a)') ' segs_have_intersection: segment',iseg,' has zero length!'
               call write_log(1, bufout)
            endif
         endif
         tvec1(1:2) = tvec1(1:2) / len1

         ! outward normal, on the left hand side seen from the segment
         nvec1 = (/ -tvec1(2), tvec1(1) /)

         ! new local coordinates [t1,n1]; loc = rot * (glb - g1a)
         ! rot = [  t1(1), t1(2) 
         !         -t1(2), t1(1) ]
         t2a = tvec1(1) * (y2a - y1a) + tvec1(2) * (z2a - z1a)
         n2a = nvec1(1) * (y2a - y1a) + nvec1(2) * (z2a - z1a)
         t2b = tvec1(1) * (y2b - y1a) + tvec1(2) * (z2b - z1a)
         n2b = nvec1(1) * (y2b - y1a) + nvec1(2) * (z2b - z1a)


         ! check if 2a,2b lie on same side of line 1

         if (max(n2a,n2b).lt.-tiny .or. min(n2a,n2b).gt.tiny) then

            has_isec = .false.
            if (ldebug.ge.4) then
               write(bufout,'(2(a,i4),a)') ' segments',iseg,',',jseg,': same sign for local n2a,n2b'
               call write_log(1, bufout)
            endif

         elseif (max(abs(n2a),abs(n2b)).lt.tiny) then

            has_isec = .true.
            if (ldebug.ge.4) then
               write(bufout,'(2(a,i4),2a)') ' segments',iseg,',',jseg,': (near) zero local n2a,n2b,',    &
                        ' overlapping segments'
               call write_log(1, bufout)
            endif

         else

            ! compute local ti where segment 2 crosses line 1
            ti = t2a - n2a * (t2b - t2a) / (n2b - n2a)

            if (ti.lt.-tiny .or. ti.gt.len1+tiny) then

               has_isec = .false.
               if (ldebug.ge.4) then
                  write(bufout,'(2(a,i4),a)') ' segments',iseg,',',jseg,': crossing outside segment 1'
                  call write_log(1, bufout)
               endif

            else

               has_isec = .true.
               if (ldebug.ge.4) then
                  write(bufout,'(2(a,i4),a)') ' segments',iseg,',',jseg,': crossing in segment 1'
                  call write_log(1, bufout)
               endif

            endif
         endif
      endif

      end associate

   end function segs_have_intersection

!------------------------------------------------------------------------------------------------------------

   subroutine test_intersect_one
!--function: perform unit-tests for function segs_have_intersection
      implicit none
!--subroutine arguments:
!--local variables:
      logical      :: has_isec
      integer      :: ldebug
      real(kind=8) :: tiny, y1a, y1b, y2a, y2b, z1a, z1b, z2a, z2b

      !  Testing intersection routine...
      !  segments   1,   2: crossing outside segment 1
      !  segments   3,   4: crossing in segment 1
      !  segments   5,   6: crossing outside segment 1
      !  segments   7,   8: no overlap in y
      !  segments   9,  10: (near) zero local n2a,n2b, overlapping segments

      call write_log('Testing intersection routine...')

      ldebug = 4
      tiny   = 1d-10
      y1a =  0d0; y1b =  1d0
      z1a =  0d0; z1b =  0d0
      y2a =  0d0; y2b =  2d0
      z2a = -1d0; z2b =  1d0 - 3*tiny
      has_isec = segs_have_intersection(1, y1a, z1a, y1b, z1b, 2, y2a, z2a, y2b, z2b, tiny, ldebug)

      tiny   = -1d-10
      y1a =  0d0; y1b =  1d0
      z1a = -1d0; z1b = -1d0
      y2a =  0d0; y2b =  2d0
      z2a = -2d0; z2b =  2.0001d-10
      has_isec = segs_have_intersection(3, y1a, z1a, y1b, z1b, 4, y2a, z2a, y2b, z2b, tiny, ldebug)

      tiny   = 1d-10
      y1a = -1d0; y1b =  1d0
      z1a = -1d0; z1b =  0d0
      y2a =  0d0; y2b =  1d0
      z2a = -1d0; z2b = -1.0001d-10
      has_isec = segs_have_intersection(5, y1a, z1a, y1b, z1b, 6, y2a, z2a, y2b, z2b, tiny, ldebug)

      tiny   = -1d-10
      y1a =  0d0; y1b =  0d0
      z1a =  0d0; z1b =  1d0
      y2a =  0d0; y2b =  0d0
      z2a = -2d0; z2b =  0.5d0
      has_isec = segs_have_intersection(7, y1a, z1a, y1b, z1b, 8, y2a, z2a, y2b, z2b, tiny, ldebug)

      tiny   = 1d-10
      y1a =  0d0; y1b =  1d0
      z1a =  0d0; z1b =  0d0
      y2a = -2d0; y2b = -1d-12
      z2a =  0d0; z2b =  0d0
      has_isec = segs_have_intersection(9, y1a, z1a, y1b, z1b, 10, y2a, z2a, y2b, z2b, tiny, ldebug)

   end subroutine test_intersect_one

!------------------------------------------------------------------------------------------------------------

end module m_profile
