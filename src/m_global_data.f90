module m_global_data
!------------------------------------------------------------------------------------------------------------
!--PURPOSE:
!  define the Result Element data relevant for the CONTACT library and stand-alone program.
!
! Copyright 2025 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------

! incorporate all global dimensions and type-definitions

use m_hierarch_data
use m_wrprof_data
use m_wrprof
#ifdef _OPENMP
   use omp_lib, only : omp_in_parallel, omp_get_max_threads, omp_set_num_threads
#endif

public

   ! permit different behavior dependent on type of program

   integer, parameter       :: program_cexe = 1, program_clib = 2
   integer                  :: program_id   =  -1

   ! administration of "result elements" == "rail-wheel pairs" == "wheels"

   integer,      parameter  :: MAX_REid_id  = 999   ! highest REid allowed
   integer, dimension(:)    :: ix_reid(MAX_REid_id) ! position of REid in list of REids
   integer,      parameter  :: MAX_NUM_REs  = 100   ! maximum number of CONTACT result elements
   integer,      parameter  :: MAX_CPids    =   9   ! maximum number of contact patches per REid
   integer                  :: num_reids    =   0   ! counter for REid's that are in use

   ! settings w.r.t. print-output

   integer                  :: idebug_glb  = 0     ! level of print-output of activate routines
   integer                  :: ilvout_glb  = 1     ! level of print-output of calculating routines

   ! using a single out-file and inp-file (wrtinp) for all result elements together
   character(len=256)       :: overal_expnam = 'contact_addon' ! overall experiment name, ex. folder
   character(len=256)       :: overal_wrkdir  = ' '            ! working folder for all result elements
   character(len=256)       :: overal_outdir  = ' '            ! output folder for all result elements

   ! define a list of module numbers for the wheel-rail pairs

   integer,          dimension(:)   :: ire_module(MAX_NUM_REs)

   ! module 1: create an array of wtd-pointers, encompassing the whole wheel-rail pair (including gds)
   !           ire=0: special wtd for stand-alone program

   type(p_ws_track), dimension(:)   :: allwtd(0:MAX_NUM_REs)

   ! module 3: create an array of gd-pointers, for each contact patch separately (not referenced in module 1)
   !           (ire,icp)=(0, 1): special gd for stand-alone program

   type(p_probdata), dimension(:,:) :: allgds(0:MAX_NUM_REs, MAX_CPids)

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

   ! internal support functions for working with multiple Result Elements

   public  cntc_initialize_data
   public  cntc_activate_wtd
   public  cntc_activate_gd
   public  cntc_destroy_wtd
   public  cntc_destroy_gd
   public  cntc_select_units

!------------------------------------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------------------------------------

subroutine cntc_initialize_data()
!--function: initialize the result element data-structures, esp allwtd, allgds
   implicit none
!--local variables:
   integer            :: jre, jcp

   ! Initialize the indirection table for result elements

   do jre = 1, MAX_REid_id
      ix_reid(jre) = -1
   enddo

   ! Initialize the module number per jre

   ire_module(1:MAX_NUM_REs) = -1

   ! Initialize the wheel-rail track data per jre

   do jre = 0, MAX_NUM_REs
      nullify(allwtd(jre)%wtd)
   enddo

   ! Initialize the contact data per jre+jcp

   do jre = 0, MAX_NUM_REs
      do jcp = 1, MAX_CPids
         nullify(allgds(jre,jcp)%gd)
      enddo
   enddo

   ! Initialize the active wtd and contact data pointers

   nullify(wtd)
   nullify(gd)

end subroutine cntc_initialize_data

!------------------------------------------------------------------------------------------------------------

function cntc_lookup_reid(REid, imodul, ierror) result(ixre)
!--function: helper to cntc_activate: determine the sequence number 'ixre' for result element REid,
!                                     claiming a free position if needed
   implicit none
!--function result:
   integer                     :: ixre           ! sequence number in allwtd/allgds
!--subroutine arguments:
   integer,      intent(in)    :: REid           ! result element ID
   integer,      intent(in)    :: imodul         ! module number for REid
   integer,      intent(inout) :: ierror
!--local variables:

   ! claim new position in allwtd/allgds if needed

   if (REid.ge.1) then
      if (ix_reid(REid).lt.0) then
         if (num_reids.ge.MAX_NUM_REs) then
            write(bufout,'(2(a,i4),a)') 'ERROR: ire=',REid,': exceeding the max. number of REs (',      &
                        MAX_NUM_REs,'), aborting.'
            call write_log(1, bufout)
            ierror = -1
         else
            num_reids = num_reids + 1
            ix_reid(REid) = num_reids
            ire_module(num_reids) = imodul
         endif
      endif
   endif

   ! return ixre

   if (REid.eq.0) then
      ixre = 0          ! special RE used in stand-alone program
   else
      ixre = ix_reid(REid)
   endif

end function cntc_lookup_reid

!------------------------------------------------------------------------------------------------------------

subroutine cntc_activate_wtd(REid, CPid, ierror)
!--function: helper to cntc_activate: activate the wheel-track data-structure for result element REid
   implicit none
!--subroutine arguments:
   integer,      intent(in)    :: REid           ! result element ID
   integer,      intent(in)    :: CPid           ! contact patch ID
   integer,      intent(inout) :: ierror         ! error code
!--local variables:
   integer                     :: ixre
   character(len=*), parameter :: subnam = 'cntc_activate_wtd'

   ! determine the sequence number ixre, create new one if needed

   ixre = cntc_lookup_reid(REid, 1, ierror)
   if (ierror.ne.0) return

   if (idebug_glb.ge.3) then
      write(bufout,'(1x,2a,2i5,a,i5)') trim(subnam),': starting for REid,CPid=',REid,CPid,', ixre=',ixre
      call write_log(1,bufout)
   endif

   ! check if the wtd data-structure for REid has been used (allocated) before, if not, allocate

   if (.not.associated(allwtd(ixre)%wtd)) then

      if (idebug_glb.ge.2) then
         write(bufout,'(a,i4)') ' allocating   wtd for ire=',REid
         call write_log(1, bufout)
      endif

      ! allocate space for the wheel-track data-structure for REid

      allocate(allwtd(ixre)%wtd)
      wtd => allwtd(ixre)%wtd

      ! set appropriate initial values for the control/input-variables

      if (idebug_glb.ge.2) then
         write(bufout,'(a,i4)') ' initializing wtd for ire=',REid
         call write_log(1, bufout)
      endif

      ! first get the basic settings from wrprof_init

      call wrprof_init(wtd)

      ! then overwrite values as needed

      wtd%meta%REid  = REid
      wtd%meta%CPid  = -1

      call cntc_select_units(wtd%scl, CNTC_un_cntc) !  Default: CONTACT units

      wtd%ic%norm    = 0              ! N-digit: z_ws position prescribed i.s.o. total vertical force
      wtd%ic%matfil_surf = 0          ! A-digit: default no .mat-file
      wtd%ic%output_surf = 0          ! O-digit: default no surf-data to .out-file
      wtd%ic%matfil_subs = 0          ! A_s:     default no .subs-file
      wtd%ic%output_subs = 0          ! O_s:     default no subs-data to .out-file
      wtd%ic%flow    = 1              ! W-digit: default only overview of calculations
      wtd%ic%return  = 1              ! R-digit: return to main program (used when inp-file is created)
      wtd%ic%wrtinp  = 0              ! wrtinp:  default no input-specification to .inp-file
      wtd%ic%ilvout  = ilvout_glb     ! ilvout:  switch off print-output ilvout_glb=0

      if (REid.eq.0) then
         wtd%meta%expnam = overal_expnam
      elseif (program_id.eq.program_cexe) then  ! stand-alone program
         write(wtd%meta%expnam,'(2a,i3.3)') trim(overal_expnam), '_re',REid
      else                                      ! C.library
         write(wtd%meta%expnam,'( a,i3.3)') 'caddon_re',REid
      endif
      wtd%meta%wrkdir = overal_wrkdir
      wtd%meta%outdir = overal_outdir
      if (idebug_glb.ge.3) then
         write(bufout,'(3(3a,/))') ' experiment name="',trim(wtd%meta%expnam),'",',                     &
                ' wrkdir="', trim(wtd%meta%wrkdir),'"', ' outdir="', trim(wtd%meta%outdir),'"'
         call write_log(3, bufout)
      endif
   endif

   ! make REid the active wheel-rail pair

   wtd => allwtd(ixre)%wtd

   ! if the calling subroutine works on cp-data, activate contact patch CPid

   if (CPid.ge.1 .and. CPid.le.MAX_CPids) then
      if (.not.associated(wtd%allcps(CPid)%cp) .or. .not.associated(wtd%allcps(CPid)%cp%gd)) then
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

end subroutine cntc_activate_wtd

!------------------------------------------------------------------------------------------------------------

subroutine cntc_activate_gd(REid, CPid, ierror)
!--function: helper to cntc_activate: activate data-structure for contact patch CPid on result element REid
   implicit none
!--subroutine arguments:
   integer,      intent(in)    :: REid           ! result element ID
   integer,      intent(in)    :: CPid           ! contact patch ID
   integer,      intent(inout) :: ierror         ! error code
!--local variables:
   integer                  :: ixre
   character(len=*), parameter :: subnam = 'cntc_activate_gd'

   ! determine the sequence number ixre, create new one if needed

   ixre = cntc_lookup_reid(REid, 3, ierror)
   if (ierror.ne.0) return

   if (idebug_glb.ge.3) then
      write(bufout,'(1x,2a,2i5,a,i5)') trim(subnam),': starting for REid,CPid=',REid,CPid,', ixre=',ixre
      call write_log(1,bufout)
   endif

   ! check whether data-structure for REid,CPid has been used (allocated) before, if not, allocate

   if (.not.associated(allgds(ixre,CPid)%gd)) then

      if (idebug_glb.ge.2) then
         write(bufout,'(2(a,i4))') ' allocating   gd for ire=',REid,', icp=',CPid
         call write_log(1, bufout)
      endif

      ! allocate space for the hierarchical data-structure for (REid,CPid)

      allocate(allgds(ixre,CPid)%gd)
      gd => allgds(ixre,CPid)%gd

      ! set appropriate initial values for the control/input-variables

      if (idebug_glb.ge.2) then
         write(bufout,'(2(a,i4))') ' initializing gd for ire=',REid,', icp=',CPid
         call write_log(1, bufout)
      endif

      ! first get defaults from gd_init

      call gd_init(gd)

      ! then overwrite selected values as needed

      gd%meta%REid  = REid
      gd%meta%CPid  = CPid

      call cntc_select_units(gd%scl, CNTC_un_cntc) !  Default: CONTACT units

      gd%ic%norm    = 0              ! N-digit: penetration prescribed i.s.o. total normal force
      gd%ic%matfil_surf = 0          ! A-digit: default no .mat-file
      gd%ic%matfil_subs = 0          ! A-digit: default no .subs-file
      gd%ic%output_surf = 0          ! O-digit: default no surf-data to .out-file
      gd%ic%output_subs = 0          ! O-digit: default no subs-data to .out-file
      gd%ic%flow    = 1              ! W-digit: default only overview of calculations
      gd%ic%wrtinp  = 0              ! wrtinp:  default no input-specification to .inp-file
      gd%ic%ilvout  = ilvout_glb     ! ilvout:  switch off print-output when ilvout_glb=0
      gd%kin%cksi   = 0d0
      gd%kin%ceta   = 0d0
      gd%kin%cphi   = 0d0            ! start with zero creepages
      gd%kin%fxrel  = 0d0
      gd%kin%fyrel  = 0d0

      if (REid.eq.0) then
         gd%meta%expnam = overal_expnam
      elseif (program_id.eq.program_cexe) then  ! stand-alone program
         write(gd%meta%expnam,'(2a,i3.3)') trim(overal_expnam), '_re',REid
      else                                      ! C.library
         write(gd%meta%expnam,'(a,i3.3,a,i1)') 'caddon_re',REid,'_cp',CPid
      endif
      gd%meta%wrkdir = overal_wrkdir
      gd%meta%outdir = overal_outdir

      if (idebug_glb.ge.3) then
         write(bufout,'(3(3a,/))') ' experiment name="',trim(gd%meta%expnam),'",',                      &
                ' wrkdir="', trim(gd%meta%wrkdir),'"', ' outdir="', trim(gd%meta%outdir),'"'
         call write_log(3, bufout)
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

end subroutine cntc_activate_gd

!------------------------------------------------------------------------------------------------------------

subroutine cntc_destroy_wtd(REid)
!--function: destroy the data-structure for wheel-rail contact problem on result element REid
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: REid           ! result element ID
!--local variables:
   integer                  :: ixre           ! position in list of REids
   character(len=*), parameter :: subnam = 'cntc_destroy_wtd'

   ! check value of REid
   ! TODO: don't stop program execution, find something better to do.

   if (REid.lt.0 .or. REid.gt.MAX_REid_id) then
      write(bufout,'(a,i6,a,i4,a)') 'ERROR: ire=',REid,': ID too large [0,',MAX_REid_id,'], aborting.'
      call write_log(1, bufout)
      call abort_run()
   endif

   ! check whether data-structure for REid has been used (allocated), if not, return

   ixre = ix_reid(REid)
   if (.not.associated(allwtd(ixre)%wtd)) return

   ! clean up all space allocated in the hierarchical data-structure wtd for (REid)

   call wrprof_destroy( allwtd(ixre)%wtd )

   ! destroy the wtd structure itself

   deallocate(allwtd(ixre)%wtd)
   nullify(allwtd(ixre)%wtd)

end subroutine cntc_destroy_wtd

!------------------------------------------------------------------------------------------------------------

subroutine cntc_destroy_gd(REid, CPid)
!--function: destroy the contact hierarchical data-structure for contact patch CPid on result element REid
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: REid           ! result element ID
   integer,      intent(in) :: CPid           ! contact patch ID
!--local variables:
   integer                  :: ixre           ! position in list of REids
   character(len=*), parameter :: subnam = 'cntc_destroy_gd'

   ! check values of REid, CPid
   ! TODO: don't stop program execution, find something better to do.

   if (REid.lt.0 .or. REid.gt.MAX_REid_id) then
      write(bufout,'(a,i6,a,i4,a)') 'ERROR: ire=',REid,': ID too large [0,',MAX_REid_id,'], aborting.'
      call write_log(1, bufout)
      call abort_run()
   endif
   if (CPid.le.0 .or. CPid.gt.MAX_CPids) then
      write(bufout,'(a,i6,a,i2,a)') 'ERROR: icp=',CPid,': number too large [1,',MAX_CPids,'], aborting.'
      call write_log(1, bufout)
      call abort_run()
   endif

   ! check whether data-structure for REid,CPid has been used (allocated), if not, return

   ixre = ix_reid(REid)
   if (.not.associated(allgds(ixre,CPid)%gd)) return

   ! clean up all space allocated in the hierarchical data-structure gd for (REid,CPid)

   call gd_destroy( allgds(ixre,CPid)%gd )

   ! destroy the gd structure itself

   deallocate(allgds(ixre,CPid)%gd)
   nullify(allgds(ixre,CPid)%gd)

end subroutine cntc_destroy_gd

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

end module m_global_data
