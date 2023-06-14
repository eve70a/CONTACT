module m_caddon_data
!------------------------------------------------------------------------------------------------------------
!--PURPOSE:
!  define all the relevant state variables for the CONTACT library,
!  including a single global variable gd with hierarchical data-structure.
!
!--INTERFACE VERSIONS:
!    1      Jul-2012   all versions up to the first release
!    2   18-Jan-2013   changed setFrictionCoefficient to setFrictionMethod
!    3   01-Sep-2017   reorganized/extended for w/r contact processing
!  MMmp  01-Sep-2020   changed numbering to MAJOR (MM) minor (m) patch (p)
!
! Copyright 2008-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------------------
! macros to append underscores to exported subroutine names

#define CONCAT_(a)    a ## _
#define STRINGIFY_(a) # a
#define STRINGIFY(a)  STRINGIFY_(a)

#ifdef TARGET_nc
#define CNAME_(a)     STRINGIFY( CONCAT_(a) )
#else
#define CNAME_(a)     STRINGIFY( a )
#endif

! incorporate all global dimensions and type-definitions

use m_hierarch_data
use m_wrprof_data
use m_wrprof
#ifdef _OPENMP
   use omp_lib, only : omp_in_parallel, omp_get_max_threads, omp_set_num_threads
#endif

public

   ! declare configuration variables for the CONTACT add-on

   integer                  :: caddon_initialized = -3 ! -3: virgin state; -2: file units set;
                                                       ! -1: finalized; 1: initialized
   integer,      parameter  :: caddon_version     = 100 * my_version_major + 10 * my_version_minor      &
                                                        + my_version_update

   ! administration of "result elements" == "rail-wheel pairs" == "wheels"

   integer,      parameter  :: MAX_REid_id  = 999   ! highest REid allowed
   integer, dimension(:)    :: ix_reid(MAX_REid_id) ! position of REid in list of REids
   integer,      parameter  :: MAX_NUM_REs  = 100   ! maximum number of CONTACT result elements
   integer,      parameter  :: MAX_CPids    =   9   ! maximum number of contact patches per REid
   integer                  :: num_reids    =   0   ! counter for REid's that are in use

   ! settings w.r.t. print-output of the CONTACT add-on

   integer                  :: idebug = 0         ! level of print-output of the addon itself
                                                  ! -1 = none, 0 = timers (default), 1 = calculating/done,
                                                  !  2 = full in+output, 3+ = debug
   character(*), parameter  :: pfx = ' ** '       ! prefix for printing output of addon
   integer                  :: timer_output = 1   ! 0 = none, 1 = overview, 2 = full, 3+ = debug
                                                  ! the value is increased when idebug is increased
   character(len=256)       :: caddon_expnam = 'contact_addon' ! experiment name
   character(len=256)       :: caddon_outpath = ' '            ! full path of output directory

   ! internal data w.r.t. OpenMP and/or explicit multi-threading:

   integer                  :: using_openmp = 0     ! 1=yes, 0=unknown, -1=no

   integer                  :: orig_mxthrd          ! original number of threads available on the machine,
                                                    ! before we set omp_num_threads=1
   integer                  :: my_thread = -1       ! thread ID, -1 == uninitialized; used from m_print_output
   !$omp         threadprivate  ( my_thread )
   integer                  :: num_threads = 0      ! total number of threads
   integer                  :: num_actv_threads = 0 ! number of threads active in cntc_calculate
   logical,       parameter :: debug_openmp = .false.

   ! internal data for the timing mechanism

   integer                  :: numtmr = MAX_NUM_REs * (1+MAX_CPids) ! number of timers for REid+CPid's
   integer                  :: ioftmr             ! offset of timers for REid+CPid's

   ! define a list of module numbers for the wheel-rail pairs

   integer,          dimension(:)   :: ire_module(MAX_NUM_REs)

   ! create an array of wtd-pointers for module 1, encompassing the whole wheel-rail pair

   type(p_ws_track), dimension(:,:) :: allwtd(MAX_NUM_REs)

   ! create an array of gd-pointers for module 3, for each contact patch separately

   type(p_probdata), dimension(:,:) :: allgds(MAX_NUM_REs, MAX_CPids)

   ! declare pointers to the "active wheel-rail pair data", "active contact patch data", etc.

   type(t_ws_track), pointer        :: wtd        => NULL()
   type(t_probdata), pointer        :: gd         => NULL()
   type(t_metadata), pointer        :: my_meta    => NULL()
   type(t_scaling),  pointer        :: my_scl     => NULL()
   type(t_ic),       pointer        :: my_ic      => NULL()
   type(t_material), pointer        :: my_mater   => NULL()
   type(t_kincns),   pointer        :: my_kin     => NULL()
   type(t_solvers),  pointer        :: my_solv    => NULL()
   type(t_subsurf),  pointer        :: my_subs    => NULL()
!$omp                threadprivate   ( wtd, gd, my_meta, my_scl, my_ic, my_mater, my_kin, my_solv, my_subs )

   ! internal support functions for working with multiple contact patches and timers

   public  cntc_activate
   private cntc_activate_wtd
   private cntc_activate_gd
   private set_actv_thread
   public  cntc_destroy_gd
   public  cntc_timer_num
   public  cntc_select_units
   public  cntc_log_start

!------------------------------------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------------------------------------

subroutine cntc_activate(REid, CPid, needs_module, needs_icp, calnam, ierror)
!--function: check needs of the calling routine and activate the wheel-track data (if module=1) and/or
!               corresponding contact patch data (if module=1, 3)
!    needs_module = 0 : calling routine works for module 1 or module 3
!                   1 : calling routine works for module 1 only
!                   3 : calling routine works for module 3 only
!    needs_icp   = -1 : in module 1, calling routine works on wtd   , on overall problem
!                   0 : in module 1, calling routine works on wtd/cp, on overall problem or on one patch only
!                   1 : in module 1, calling routine works on     cp, on one contact patch only
   implicit none
!--subroutine arguments:
   integer,          intent(in)  :: REid         ! result element ID
   integer,          intent(in)  :: CPid         ! contact patch ID
   integer,          intent(in)  :: needs_module ! module(s) supported by the calling routine (0, 1, 3)
   integer,          intent(in)  :: needs_icp    ! icp(s) supported by the calling routine (-1, 0, 1)
   character(len=*), intent(in)  :: calnam       ! name of calling subroutine
   integer,          intent(out) :: ierror       ! error code
!--local variables:
   integer                  :: ixre           ! position in list of REids
   integer                  :: imodul, ifcver, len_outpath
   character(len=1)         :: c_outpath(20)
   character(len=*), parameter :: subnam = 'cntc_activate'
!--subroutines used:
   interface
      subroutine cntc_initialize(ire, imodul, ifcver, ierror, c_outpath, len_outpath) &
         bind(c,name=CNAME_(cntc_initialize))
         integer,          intent(in)    :: ire         ! result element ID
         integer,          intent(in)    :: imodul      ! module number 1=w/r contact, 3=basic contact
         integer,          intent(out)   :: ifcver      ! version of the CONTACT add-on
         integer,          intent(inout) :: ierror      ! error flag
         character,        intent(in)    :: c_outpath(*) ! C-string: full path of output directory
         integer,          intent(in)    :: len_outpath ! length of C-string
      end subroutine cntc_initialize
   end interface

   if (idebug.ge.5) call cntc_log_start(subnam, .true.)

   ierror = 0

   ! check value of REid

   if (REid.le.0 .or. REid.gt.MAX_REid_id) then
      write(bufout,'(a,i6,a,i6,a)') ' ERROR: ire=',REid,' ID too large [1,',MAX_REid_id,'], aborting.'
      call write_log(1, bufout)
      ierror = -101
      return
   endif

   ! check whether the call is made by a known or a new thread, if needed, set the thread number

   call set_actv_thread(REid, CPid)

   ! check whether result element REid has been used before, if not, initialize

!$omp critical (activate_init)
   if (ix_reid(REid).le.0) then
      if (CPid.le.0) then
         imodul = 1
      else
         imodul = 3
      endif
      ierror = 0
      ifcver = 0
      c_outpath = ' '
      len_outpath = 1
      if (idebug.ge.5) call write_log(' calling cntc_initialize for REid...')
      call cntc_initialize(REid, imodul, ifcver, ierror, c_outpath, len_outpath)
      if (idebug.ge.1 .and. ierror.ne.0) then
         write(bufout,*) 'initialize: ierror=',ierror
         call write_log(1, bufout)
      endif
   endif
!$omp end critical (activate_init)

   ixre = ix_reid(REid)

   ! check whether the calling subroutine allows for the module that's used by ire

   if (ire_module(ixre).eq.1 .and. needs_module.ne.1 .and. needs_module.ne.0) then
      write(bufout,'(4a,i3,a)') ' ERROR: subroutine ', trim(calnam),' is not available for w/r contact ', &
                '(ire=', REid, '), skipping'
      call write_log(1, bufout)
      ierror = -102
   endif
   if (ire_module(ixre).eq.3 .and. needs_module.ne.3 .and. needs_module.ne.0) then
      write(bufout,'(4a,i3,a)') ' ERROR: subroutine ', trim(calnam),' is not available for basic contact ', &
                '(ire=', REid, '), skipping'
      call write_log(1, bufout)
      ierror = -103
   endif

   ! check if CPid is valid: -1 or [1, MAX_CPids]   Value -1 means "work on w/r wtd data" (module 1)

   if (ire_module(ixre).eq.1 .and. (CPid.le.-2 .or. CPid.eq.0 .or. CPid.gt.MAX_CPids)) then
      write(bufout,'(a,i6,a,i6,a)') ' ERROR: icp=',CPid,' out of range [-1] or [1,',MAX_CPids,'], aborting.'
      call write_log(1, bufout)
      ierror = -104
   elseif (ire_module(ixre).eq.3 .and. (CPid.le.0 .or. CPid.gt.MAX_CPids)) then
      write(bufout,'(a,i6,a,i6,a)') ' ERROR: icp=',CPid,' out of range [1,',MAX_CPids,'], aborting.'
      call write_log(1, bufout)
      ierror = -104
   endif

   ! check if CPid is available in calling routine

   if (CPid.eq.-1 .and. needs_icp.ge.1) then
      write(bufout,'(4a)') ' ERROR: subroutine ', trim(calnam),' needs a valid icp number icp>0,',    &
                ' skipping calculation.'
      call write_log(1, bufout)
      ierror = -105
   elseif (ire_module(ixre).eq.1 .and. needs_icp.le.-1 .and. CPid.gt.0 ) then
      write(bufout,'(3a,i3,a)') ' ERROR: subroutine ', trim(calnam),' does not support icp>0 (', CPid,    &
                '), skipping calculation.'
      call write_log(1, bufout)
      ierror = -106
   endif
   if (ierror.ne.0) return

   ! activate wheel-track data or contact patch data depending on CPid
   ! including allocation & initialization & setting pointers

   if (ire_module(ixre).eq.1) then

      if (idebug.ge.5) call write_log(' calling cntc_activate_wtd for REid...')
      call cntc_activate_wtd(REid, CPid, ierror)

   elseif (CPid.ge.1 .and. CPid.le.MAX_CPids) then

      if (idebug.ge.5) call write_log(' calling cntc_activate_gd for REid,CPid...')
      call cntc_activate_gd(REid, CPid)

   endif

   if (idebug.ge.5) call cntc_log_start(subnam, .false.)
end subroutine cntc_activate

!------------------------------------------------------------------------------------------------------------

subroutine cntc_activate_wtd(REid, CPid, ierror)
!--function: helper to cntc_activate: activate the wheel-track data-structure for result element REid
   implicit none
!--subroutine arguments:
   integer,      intent(in)    :: REid           ! result element ID
   integer,      intent(in)    :: CPid           ! contact patch ID
   integer,      intent(inout) :: ierror         ! error code
!--local variables:
   integer                     :: ixre           ! position in list of REids
   character(len=*), parameter :: subnam = 'cntc_activate_wtd'

   if (idebug.ge.5) call cntc_log_start(subnam, .true.)

   ! check whether wtd data-structure for REid has been used (allocated) before, if not, allocate

   ixre = ix_reid(REid)
   if (.not.associated(allwtd(ixre)%wtd)) then

      if (idebug.ge.5) then
         write(bufout,'(a,i4)') ' allocating wtd for ire=',REid
         call write_log(1, bufout)
      endif

      ! allocate space for the wheel-track data-structure for REid

      allocate(allwtd(ixre)%wtd)
      wtd => allwtd(ixre)%wtd

      ! set appropriate initial values for the control/input-variables

      if (idebug.ge.5) then
         write(bufout,'(a,i4)') ' initializing wtd for ire=',REid
         call write_log(1, bufout)
      endif

      ! first get the basic settings from wrprof_ini

      call wrprof_ini(wtd)

      ! then overwrite values as needed

      wtd%meta%REid  = REid
      wtd%meta%CPid  = -1
      wtd%meta%ncase =  1

      call cntc_select_units(wtd%scl, CNTC_un_cntc) !  Default: CONTACT units

      wtd%ic%norm    = 0              ! N-digit: z_ws position prescribed i.s.o. total vertical force
      wtd%ic%matfil_surf = 0          ! A-digit: default no .mat-file
      wtd%ic%output_surf = 0          ! O-digit: default no surf-data to .out-file
      wtd%ic%matfil_subs = 0          ! A_s:     default no .subs-file
      wtd%ic%output_subs = 0          ! O_s:     default no subs-data to .out-file
      wtd%ic%flow    = 1              ! W-digit: default only overview of calculations
      wtd%ic%return  = 1              ! R-digit: return to main program (used when inp-file is created)
      wtd%ic%wrtinp  = 0              ! wrtinp:  default no input-specification to .inp-file
      wtd%ic%ilvout  = min(idebug, 1) ! ilvout:  switch off output when idebug=0

      if (caddon_outpath.eq.' ') then
         write(wtd%meta%expnam,'( a,i3.3)')                      'wraddon_re',REid
      else
         write(wtd%meta%expnam,'(2a,i3.3)') trim(caddon_outpath),'wraddon_re',REid
      endif
      if (idebug.ge.5) then
         write(bufout,*) 'experiment name=',trim(wtd%meta%expnam)
         call write_log(1, bufout)
      endif
   endif

   ! make REid the active wheel-rail pair

   wtd => allwtd(ixre)%wtd

   ! if the calling subroutine works on cp-data, activate contact patch CPid

   if (CPid.ge.1 .and. CPid.le.MAX_CPids) then
      if (.not.associated(wtd%allcps(CPid)%cp%gd)) then
         write(bufout,'(a,i2,a,i4,a)') ' ERROR: contact patch',CPid,' of wheel-rail problem',REid,      &
                ' hasnt been initialized yet.'
         call write_log(1, bufout)
         ierror = CNTC_err_icp
         return
      else
         gd       => wtd%allcps(CPid)%cp%gd
         my_scl   => wtd%scl    ! note: scl isnt copied from wtd to the cps
         my_meta  => wtd%meta   ! note: meta isnt copied from wtd to the cps
         my_ic    => gd%ic
         my_mater => gd%mater
         my_kin   => gd%kin
         my_solv  => gd%solv
         my_subs  => gd%subs
      endif

   ! else, the calling subroutine works on the wtd-problem as a whole, set pointers to its configuration data

   else
      gd       => NULL()
      my_meta  => wtd%meta
      my_scl   => wtd%scl
      my_ic    => wtd%ic
      my_mater => wtd%mater
      my_kin   => wtd%kin
      my_solv  => wtd%solv
      my_subs  => wtd%subs
   endif

   ! check if the wheel-rail pair (REid) is already in use (locked) by a different thread
   ! If so, write error message and abort

   if (my_meta%actv_thrd.ge.0) then
      write(bufout,'(a,i3,2a,i2,a)') ' ERROR: wheel-rail problem (',REid, ') is already in use (locked)', &
        ' by thread', my_meta%actv_thrd,', aborting.'
      call write_log(1, bufout)
      ierror = 305
   endif

   if (idebug.ge.5) call cntc_log_start(subnam, .false.)
end subroutine cntc_activate_wtd

!------------------------------------------------------------------------------------------------------------

subroutine cntc_activate_gd(REid, CPid)
!--function: helper to cntc_activate: activate data-structure for contact patch CPid on result element REid
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: REid           ! result element ID
   integer,      intent(in) :: CPid           ! contact patch ID
!--local variables:
   integer                  :: ixre           ! position in list of REids
   character(len=*), parameter :: subnam = 'cntc_activate_gd'

   if (idebug.ge.5) call cntc_log_start(subnam, .true.)

   ! check whether data-structure for REid,CPid has been used (allocated) before, if not, allocate

   ixre = ix_reid(REid)
   if (idebug.ge.5) then
      write(bufout,*) trim(subnam),': starting for REid,CPid=',REid,CPid,', ixre=',ixre
      call write_log(1,bufout)
   endif

   if (.not.associated(allgds(ixre,CPid)%gd)) then

      if (idebug.ge.5) then
         write(bufout,'(2(a,i4))') ' allocating gd for ire=',REid,', icp=',CPid
         call write_log(1, bufout)
      endif

      ! allocate space for the hierarchical data-structure for (REid,CPid)

      allocate(allgds(ixre,CPid)%gd)
      gd => allgds(ixre,CPid)%gd

      ! set appropriate initial values for the control/input-variables

      if (idebug.ge.5) then
         write(bufout,'(2(a,i4))') ' initializing gd for ire=',REid,', icp=',CPid
         call write_log(1, bufout)
      endif

      ! first get defaults from setini

      call setini(gd)

      ! then overwrite selected values as needed

      gd%meta%REid  = REid
      gd%meta%CPid  = CPid
      gd%meta%ncase = 1

      call cntc_select_units(gd%scl, CNTC_un_cntc) !  Default: CONTACT units

      gd%ic%norm    = 0              ! N-digit: penetration prescribed i.s.o. total normal force
      gd%ic%matfil_surf = 0          ! A-digit: default no .mat-file
      gd%ic%matfil_subs = 0          ! A-digit: default no .subs-file
      gd%ic%output_surf = 0          ! O-digit: default no surf-data to .out-file
      gd%ic%output_subs = 0          ! O-digit: default no subs-data to .out-file
      gd%ic%flow    = 1              ! W-digit: default only overview of calculations
      gd%ic%wrtinp  = 0              ! wrtinp:  default no input-specification to .inp-file
      gd%ic%ilvout  = min(idebug, 1) ! ilvout:  switch off output when idebug=0
      gd%kin%cksi   = 0d0
      gd%kin%ceta   = 0d0
      gd%kin%cphi   = 0d0            ! start with zero creepages
      gd%kin%fxrel1 = 0d0
      gd%kin%fyrel1 = 0d0

      if (caddon_outpath.eq.' ') then
         write(gd%meta%expnam,'(a,i3.3,a,i1)') 'caddon_re',REid,'_cp',CPid
      else
         write(gd%meta%expnam,'(a,a,i3.3,a,i1)') trim(caddon_outpath),'caddon_re',REid,'_cp',CPid
      endif

      if (idebug.ge.5) then
         write(bufout,*) 'experiment name=',trim(gd%meta%expnam)
         call write_log(1, bufout)
      endif
   endif

   ! set pointers to the relevant data-structures for contact patch (REid, CPid),
   ! i.e. make it the active contact patch

   gd => allgds(ixre,CPid)%gd
   my_meta    => gd%meta
   my_scl     => gd%scl
   my_ic      => gd%ic
   my_mater   => gd%mater
   my_kin     => gd%kin
   my_solv    => gd%solv
   my_subs    => gd%subs

   ! check if the contact patch (REid,CPid) is already in use (locked) by a different thread
   ! If so, write error message and abort

   if (gd%meta%actv_thrd.ge.0) then
      write(bufout,'(a,i3,a,i1,a,i2,a)') 'ERROR: contact problem (',REid,'.',CPid,                  &
                ') is already in use (locked) by thread', gd%meta%actv_thrd,', aborting.'
      call write_log(1, bufout)
      call abort_run()
   endif

   if (idebug.ge.5) call cntc_log_start(subnam, .false.)
end subroutine cntc_activate_gd

!------------------------------------------------------------------------------------------------------------

subroutine set_actv_thread(REid, CPid)
!--function: check whether the call is made by a known or a new thread, set the thread number
#ifdef _OPENMP
   use omp_lib    , only :  omp_get_thread_num, omp_in_parallel
#endif
   implicit none
!--subroutine arguments:
   integer,          intent(in)  :: REid         ! result element ID
   integer,          intent(in)  :: CPid         ! contact patch ID
!--local variables:

   ! check whether the call is made by a known or a new thread.
   ! new thread: determine OpenMP or explict threading; set the thread number

#ifdef _OPENMP
   if (my_thread.lt.0) then
      if (idebug.ge.4) then
         write(bufout,'(a,i3,a,i1,a)') ' a new thread is starting, (ire,icp)=(',REid,',',CPid,')'
         call write_log(1, bufout)
      endif

      ! first update information on using OpenMP or not
      !  - if we're in a parallel loop then OpenMP is used
      !  - if we don't know about OpenMP and this is not the first thread, then explicit threading

      if (using_openmp.eq.0) then
         if (omp_in_parallel().ne.0) then
            using_openmp = 1
         elseif (num_threads.gt.0) then
            using_openmp = -1
         endif
      endif

      ! set the thread number

      if (using_openmp.ge.1) then
         my_thread = omp_get_thread_num()
         num_threads = num_threads + 1
      elseif (using_openmp.le.-1) then
         my_thread = num_threads
         num_threads = num_threads + 1
      else ! using_openmp=0, don't know if OpenMP is used or not
         ! serial section or first explicit thread
         my_thread = 0
         num_threads = num_threads + 1
      endif

      ! tell the thread number to the timing routines and to the output routines

      call timer_thread_config(my_thread)
      call set_print_thread(my_thread)

      if (idebug.ge.4) then
         write(bufout,'(a,i2,a,i3,a,i1,a)') ' thread',my_thread,': starting for contact problem (',     &
                REid,'.',CPid,')'
         call write_log(1, bufout)
      endif
   endif
#else

   ! compiled without OpenMP, must be explict threading

   if (my_thread.lt.0) then
      if (idebug.ge.4) then
         write(bufout,'(a,i3,a,i1,a)') ' a new thread is starting, (ire,icp)=(',REid,',',CPid,')'
         call write_log(1, bufout)
      endif

      using_openmp = -1
      my_thread = num_threads
      num_threads = num_threads + 1
   endif
#endif

end subroutine set_actv_thread

!------------------------------------------------------------------------------------------------------------

subroutine cntc_destroy_gd(REid, CPid)
!--function: destroy the contact hierarchical data-structure for contact patch CPid on result element REid
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: REid           ! result element ID
   integer,      intent(in) :: CPid           ! contact patch ID
!--local variables:
   integer                  :: ixre           ! position in list of REids
   integer                  :: iblk
   character(len=*), parameter :: subnam = 'cntc_destroy_gd'

   if (idebug.ge.5) call cntc_log_start(subnam, .true.)

   ! check values of REid, CPid
   ! TODO: don't stop program execution, find something better to do.

   if (REid.le.0 .or. REid.gt.MAX_REid_id) then
      write(bufout,'(a,i6,a,i4,a)') 'ERROR: ire=',REid,' ID too large [1,',MAX_REid_id,'], aborting.'
      call write_log(1, bufout)
      call abort_run()
   endif
   if (CPid.le.0 .or. CPid.gt.MAX_CPids) then
      write(bufout,'(a,i6,a,i2,a)') 'ERROR: icp=',CPid,' number too large [1,',MAX_CPids,'], aborting.'
      call write_log(1, bufout)
      call abort_run()
   endif

   ! check whether data-structure for REid,CPid has been used (allocated), if not, return

   ixre = ix_reid(REid)
   if (.not.associated(allgds(ixre,CPid)%gd)) return

   ! clean up all space allocated for the hierarchical data-structure gd for (REid,CPid)

   gd => allgds(ixre,CPid)%gd

   ! cleanup of grid:

   call grid_destroy(gd%cgrid)

   ! allocatable & pointer arrays of t_geomet:

   call destroy_arr(gd%geom%prmudf)
   call destroy_arr(gd%geom%prmpln)
   call destroy_arr(gd%geom%ysep)
   call destroy_arr(gd%geom%facsep)
   ! TODO: cleanup subsurf.input
   !! call destroy_arr(gd%geom%xb)

   ! cleanup of influence coefficients influ:

   call inflcf_destroy(gd%influ%cs)
   call inflcf_destroy(gd%influ%cv)
   call inflcf_destroy(gd%influ%csv)
   call inflcf_destroy(gd%influ%ms)

   ! cleanup of output arrays outpt1:

   call eldiv_destroy(gd%outpt1%igs)
   call eldiv_destroy(gd%outpt1%igv)
   call destroy_arr(gd%outpt1%ledg%jbnd)
   call destroy_arr(gd%outpt1%ledg%ixbnd)
   call if1_destroy(gd%outpt1%ledg%ii2j)
   call destroy_arr(gd%outpt1%ledg%xbnd)
   call destroy_arr(gd%outpt1%ledg%facdx)
   call destroy_arr(gd%outpt1%ledg%ubnd)
   call gf3_destroy(gd%outpt1%ledg%facdt)
   call gf3_destroy(gd%outpt1%mus)
   call gf3_destroy(gd%outpt1%muv)
   call gf3_destroy(gd%outpt1%shft)
   call gf3_destroy(gd%outpt1%ps)
   call gf3_destroy(gd%outpt1%pv)
   call gf3_destroy(gd%outpt1%us)
   call gf3_destroy(gd%outpt1%uv)
   call gf3_destroy(gd%outpt1%ss)
   call gf3_destroy(gd%outpt1%sv)

   ! cleanup subsurf data

   do iblk = 1, gd%subs%nblock
      call subsblk_destroy( gd%subs%blocks(iblk) )
   enddo

   ! destroy the whole structure t_probdata

   deallocate(allgds(ixre,CPid)%gd)
   nullify(allgds(ixre,CPid)%gd)
   nullify(gd)

   if (idebug.ge.5) call cntc_log_start(subnam, .false.)
end subroutine cntc_destroy_gd

!------------------------------------------------------------------------------------------------------------

integer function cntc_timer_num(REid, CPid)
!--function: return the timer number for contact patch (REid,CPid)
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: REid           ! result element ID
   integer,      intent(in) :: CPid           ! contact patch ID
!--local variables:
   integer                  :: ixre           ! position in list of REids

   ixre = ix_reid(REid)
   if (CPid.le.0) then
      cntc_timer_num = ioftmr + 1        + (ixre-1)*(1+MAX_CPids)
   else
      cntc_timer_num = ioftmr + 1 + CPid + (ixre-1)*(1+MAX_CPids)
   endif
end function cntc_timer_num

!------------------------------------------------------------------------------------------------------------

subroutine cntc_select_units(scl, iunits)
!--function: set the scaling factors for unit scheme iunits
   implicit none
!--subroutine arguments:
   type(t_scaling), intent(out) :: scl
   integer,         intent(in)  :: iunits

   if (iunits.eq.CNTC_un_cntc) then

      ! selecting CONTACT's unit convention: length in [mm], veloc in [mm/s], force in [N], acting on body 1

      scl%units = CNTC_un_cntc
      scl%len   =  1d0                          ! [mm]   / [unit length]
      scl%area  =  1d0                          ! [mm^2] / [unit area]
      scl%forc  =  1d0                          ! [N]    / [unit force]
      scl%veloc =  1d0                          ! [mm/s] / [unit speed]
      scl%angle =  1d0                          ! [rad]  / [unit angle]
      scl%body  =  1d0                          ! [on body 1] / [on output body]

   elseif (iunits.eq.CNTC_un_spck) then

      ! selecting SIMPACK's unit convention: length in [m], veloc in [m/s], force in [N], acting on body 2

      scl%units = CNTC_un_spck
      scl%len   =  1d3                          ! [mm]   / [m]
      scl%area  =  1d6                          ! [mm^2] / [m^2]
      scl%forc  =  1d0                          ! [N]    / [N]
      scl%veloc =  1d3                          ! [mm/s] / [m/s]
      scl%angle =  1d0                          ! [rad]  / [rad]
      scl%body  = -1d0                          ! [on body 1] / [on body 2]

   elseif (iunits.eq.CNTC_un_si) then

      ! selecting SI unit convention: length in [m], veloc in [m/s], force in [N], acting on body 1

      scl%units = CNTC_un_si
      scl%len   =  1d3                          ! [mm]   / [m]
      scl%area  =  1d6                          ! [mm^2] / [m^2]
      scl%forc  =  1d0                          ! [N]    / [N]
      scl%veloc =  1d3                          ! [mm/s] / [m/s]
      scl%angle =  1d0                          ! [rad]  / [rad]
      scl%body  =  1d0                          ! [on body 1] / [on output body]

   elseif (iunits.eq.CNTC_un_imper) then

      ! selecting imperial units (n.y.a.): length in [in], force in [lbf], acting on body 1

      call write_log('imperial units not yet available')
      scl%units = CNTC_un_imper
      scl%len   =  1d0                          ! [mm]   / [unit length]
      scl%area  =  1d0                          ! [mm^2] / [unit area]
      scl%forc  =  1d0                          ! [N]    / [unit force]
      scl%veloc =  1d0                          ! [mm/s] / [unit speed]
      scl%angle =  1d0                          ! [rad]  / [rad]
      scl%body  =  1d0                          ! [on body 1] / [on output body]

   else

      write(bufout,*) 'Invalid value for flag CNTC_if_units:',iunits
      call write_log(1, bufout)

   endif

end subroutine cntc_select_units

!------------------------------------------------------------------------------------------------------------

subroutine cntc_log_start(subnam, lstart)
!--function: print log-message for starting/returning of interface function
   implicit none
!--subroutine arguments:
   character(len=*), intent(in) :: subnam     ! interface function name
   logical,          intent(in) :: lstart     ! start=T, return=F
!--local variables:

   ! in parallel runs, include the thread-number

   if (num_threads.gt.1) then
      if (num_threads.le.9) then
         if (lstart) then
            write(bufout,'(a,a31,a,i1,a)') pfx,subnam,'(tr.',my_thread,'): starting...'
         else
            write(bufout,'(a,a31,a,i1,a)') pfx,subnam,'(tr.',my_thread,'): returning...'
         endif
      else
         if (lstart) then
            write(bufout,'(a,a29,a,i3,a)') pfx,subnam,'(thr',my_thread,'): starting...'
         else
            write(bufout,'(a,a29,a,i3,a)') pfx,subnam,'(thr',my_thread,'): returning...'
         endif
      endif
   else
      if (lstart) then
         write(bufout,'(a,a37,a)') pfx,subnam,': starting...'
      else
         write(bufout,'(a,a37,a)') pfx,subnam,': returning...'
      endif
   endif
   call write_log(1, bufout)

end subroutine cntc_log_start

!------------------------------------------------------------------------------------------------------------

end module m_caddon_data

