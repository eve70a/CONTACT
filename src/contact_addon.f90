!============================================================================================================
! m_contact_addon - interface functions for CONTACT library version
!
! Copyright 2008-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!============================================================================================================

module m_contact_addon

use, intrinsic        :: iso_c_binding, only: C_CHAR, C_NULL_CHAR
use m_caddon_data
use m_licensing
use m_sinput
use m_scontc
use m_aijpj
use m_subsurf
use m_wr_input
use m_wr_profiles
use m_wr_totforce
implicit none
private

private cntc_setFileUnits
private copy_c_dirnam
private cntc_initializeFirst
private cntc_initialize
private cntc_setGlobalFlags
private cntc_readInpFile
private cntc_setFlags
private cntc_setMetadata
private cntc_setSolverFlags
private cntc_setMaterialParameters
private cntc_setTemperatureData
private cntc_setTimestep
private cntc_setReferenceVelocity
private cntc_setRollingStepsize
private cntc_setFrictionMethod
private cntc_setHertzContact
private cntc_setPotContact
private cntc_setVerticalForce
private cntc_setPenetration
private cntc_setNormalForce
private cntc_setUndeformedDistc
private cntc_setCreepages
private cntc_setExtraRigidSlip
private cntc_setTangentialForces
private cntc_setProfileInputFname
private cntc_setProfileInputValues
private cntc_setTrackDimensions
private cntc_setWheelsetDimensions
private cntc_setWheelsetPosition
private cntc_setWheelsetVelocity
private cntc_setWheelsetFlexibility
private subs_addBlock
private cntc_calculate
private cntc_calculate1
private cntc_calculate3
private subs_calculate
private cntc_getFlags
private cntc_getParameters
private cntc_getProfileValues_new
private cntc_getProfileValues
private cntc_getWheelsetPosition
private cntc_getWheelsetVelocity
private cntc_getNumContactPatches
private cntc_getContactLocation
private cntc_getReferenceVelocity
private cntc_getHertzContact
private cntc_getNumElements
private cntc_getGridDiscretization
private cntc_getPotContact
private cntc_getPenetration
private cntc_getCreepages
private cntc_getContactForces
private cntc_getGlobalForces
private cntc_getContactPatchAreas
private cntc_getElementDivision
private cntc_getMaximumPressure
private cntc_getMaximumTraction
private cntc_getMaximumTemperature
private cntc_getFieldData
private cntc_getTractions
private cntc_getMicroSlip
private cntc_getDisplacements
private cntc_getSensitivities
private cntc_getCalculationTime
private cntc_resetCalculationTime
private subs_getBlocksize
private subs_getResults
private finalize_helper
private cntc_finalize
private cntc_finalizeLast

contains

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

! CONCAT_(hello)               expands to hello_
! STRINGIFY( CONCAT_(hello) )  expands to STRINGIFY_( hello_ )
! STRINGIFY_( CONCAT_(hello) ) expands to "hello_"
! CNAME_(hello)           thus expands to "hello_"

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setFileUnits(nunits, iunits) &
   bind(c,name=CNAME_(cntc_setfileunits))
!--function: used for configuring the logical file units used by the CONTACT library.
!            Five values are requested: positive, distinct and not containing 5 or 6.
!--category: 0, "m=any, glob":      not related to contact's modules 1 or 3, working on global data
   implicit none
!--subroutine arguments:
   integer, intent(in) :: nunits         ! length of array iunits, number of file units provided
   integer, intent(in) :: iunits(nunits) ! logical file units to be used by the CONTACT library
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_setFileUnits'
   integer, parameter :: num_files = 7
   integer            :: my_units(num_files)
   logical            :: lfound
   integer            :: iin, imy, jmy        ! loop counter, array indices
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setFileUnits
#endif

   if (caddon_initialized.ge.-2) then
      call write_log('ERROR: cntc_setFileUnits may not be used after cntc_initialize.')
      return
   else
#ifdef _OPENMP
   if (omp_in_parallel().ne.0) then
      call write_log('ERROR: cntc_setFileUnits may not be called from a parallel region')
      return
   endif
#endif

      ! check unit numbers provided, copy accepted ones to my_units

      imy = 0
      do iin = 1, nunits
         if (iunits(iin).le.0) then
            ! error: logical units must be positive
            ! write(*,*) 'error: file-unit (',iunits(iin),') should be positive.'
         elseif (iunits(iin).eq.5 .or. iunits(iin).eq.6) then
            ! error: reserved for stdin and stdout
            ! write(*,*) 'error: file-unit (',iunits(iin),') is reserved.'
         else
            lfound = .false.
            do jmy = 1, imy
               lfound = lfound .or. iunits(iin).eq.my_units(jmy)
            enddo
            if (lfound) then
               ! write(*,*) 'error: file-unit (',iunits(iin),') is used already.'
            elseif (imy.ge.num_files) then
               ! write(*,*) 'warning: file-unit (',iunits(iin),') ignored.'
            else
               imy = imy + 1
               my_units(imy) = iunits(iin)
               ! write(*,*) 'accepted file-unit',imy,' = ',iunits(iin)
            endif
         endif
      enddo

      if (imy.le.0) then
         call write_log('ERROR: at least one valid file-number is needed.')
         call abort_run()
      endif

      ! assign the accepted unit numbers to various files, in order of importance
      ! if there are too few unit numbers, set remaining ones to "-1"

      lout  = -1
      linp  = -1
      ltmp1 = -1
      ltmp2 = -1
      ltmp3 = -1
      if (imy.ge.1) lout  = my_units(1)
      if (imy.ge.2) linp  = my_units(2)
      if (imy.ge.3) ltmp1 = my_units(3)
      if (imy.ge.4) ltmp2 = my_units(4)
      if (imy.ge.5) ltmp3 = my_units(5)

      caddon_initialized = -2
   endif

end subroutine cntc_setFileUnits

!------------------------------------------------------------------------------------------------------------

subroutine copy_c_dirnam(c_dirnam, len_dirnam, f_dirnam, ltmpfile)
!--function: helper to convert C string to F-string and append path_sep
   implicit none
!--subroutine arguments:
   character(kind=C_CHAR), intent(in)    :: c_dirnam(*)  ! C-string: folder name
   integer,                intent(in)    :: len_dirnam   ! length of C-string
   character(len=*),       intent(out)   :: f_dirnam     ! F-string: folder name
   logical,                intent(in)    :: ltmpfile
!--local variables
   integer              :: i, ilen

   ! Copy C to Fortran

   call c_to_f_string(c_dirnam, f_dirnam, len_dirnam)

   ! append '/' on non-empty directory that doesnt end on '/' already
   ! TODO: permit use of '/' on Windows, '\' on Linux

   ilen = len(trim(f_dirnam))
   if (ilen.gt.1 .and. ilen.lt.len(f_dirnam)) then
      if (f_dirnam(ilen:ilen).ne.path_sep) then
         i = ilen + 1
         f_dirnam(i:i) = path_sep
      endif
   endif

   if (ltmpfile) then
      write(37,*) 'len=',len(f_dirnam), len_dirnam
      write(37,*) 'path=',trim(f_dirnam),'.'
   endif

end subroutine copy_c_dirnam

!------------------------------------------------------------------------------------------------------------

subroutine cntc_initializeFirst(ifcver, ierror, ioutput, c_wrkdir, c_outdir, c_expnam, len_wrkdir,      &
                                    len_outdir, len_expnam) &
   bind(c,name=CNAME_(cntc_initializefirst))
!--function: initialize the addon internal data structures and initialize output channels,
!            print version information; return the addon version number.
!--category: 0, "m=any, glob":      not related to contact's modules 1 or 3, working on global data
   implicit none
!--subroutine arguments:
   integer,                intent(out)   :: ifcver       ! version of the CONTACT add-on
   integer,                intent(inout) :: ierror       ! error flag
   integer,                intent(in)    :: ioutput      ! output channels: 0 = out-file, 1 = out-file+screen
   character(kind=C_CHAR), intent(in)    :: c_wrkdir(*)  ! C-string: effective working folder
   character(kind=C_CHAR), intent(in)    :: c_outdir(*)  ! C-string: output folder
   character(kind=C_CHAR), intent(in)    :: c_expnam(*)  ! C-string: experiment name
   integer,                intent(in)    :: len_wrkdir   ! length of C-string
   integer,                intent(in)    :: len_outdir   ! length of C-string
   integer,                intent(in)    :: len_expnam   ! length of C-string
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_initializeFirst'
   integer, parameter :: max_version = 20
   integer            :: num_version, ix
   character(len=256) :: version(max_version), fname, f_recent
   character          :: c_recent(256)
   character(len=len(bufout)-5) :: tmpbuf
   integer            :: jre, jcp, ilen, louttm
   logical            :: ltmpfile = .false.
   logical            :: l_outfil, l_screen, l_simpck
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_initializeFirst
#endif

   ifcver   = caddon_version
   if (caddon_initialized.ge.1) return

#ifdef _OPENMP
   if (omp_in_parallel().ne.0) then
      call write_log('ERROR: cntc_initializeFirst may not be called from a parallel region')
      return
   endif
#endif

   ! Upon first call: perform one-time-only initializations

   if (ltmpfile) then
      fname = 'C:\Temp\caddon_out.tmp'
      ! fname = 'caddon_out.tmp'
      open(unit=37, file=fname)
      write(37,*) 'initialize: starting...'
   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)

   ! Store folder names (note: only the value of the first call is memorized, used for all REs)

   call copy_c_dirnam(c_wrkdir, len_wrkdir, caddon_wrkdir, ltmpfile)
   call copy_c_dirnam(c_outdir, len_outdir, caddon_outdir, ltmpfile)

   ! If provided: store experiment name (note: only the value of the first call is memorized)

   if (len_expnam.ge.1 .and. c_expnam(1).ne.' ') then
      call c_to_f_string(c_expnam, caddon_expnam, len_expnam)

      if (ltmpfile) then
         write(37,*) 'len=',len(caddon_expnam), len_expnam
         write(37,*) 'expnam=',trim(caddon_expnam),'.'
      endif
   endif

   ! strip off directory name from experiment name, use as fallback for effective working folder

   ix = index_pathsep(caddon_expnam, back=.true.)
   if (ix.gt.1) then
      if (caddon_wrkdir.eq.' ') caddon_wrkdir  = caddon_expnam(1:ix-1)
      caddon_expnam  = caddon_expnam(ix+1:)

      if (ltmpfile) then
         write(37,*) '...directory separator found at ix=',ix
         write(37,*) '...removing dirname to form true experiment name'
         write(37,*) '...wrkdir= ',trim(caddon_wrkdir)
         write(37,*) '...expnam=  ',trim(caddon_expnam)
      endif
   endif

   ! Set desired configuration of output streams
   !   - default    : to out-file and Simpack
   !   - ioutput = 1: to out-file, screen and Simpack

   l_outfil = .true.
   l_screen = .false.
   l_simpck = .true.
   if (ioutput.eq.1) l_screen = .true.

   ! Activate the selected output streams, except the out-file which hasn't been opened yet

   if (ltmpfile) then
      write(37,*) 'opening streams-1=',.false., l_screen, l_simpck
   endif
   call set_print_streams(nw_outfil=.false., nw_screen=l_screen, nw_simpck=l_simpck)

   ! Try to open the .out-file in the requested output directory

   fname = trim(caddon_outdir) // trim(caddon_expnam) // '.out'
   if (ltmpfile) write(37,*) 'init=',caddon_initialized,', fname=',trim(fname)
   if (.false. .and. caddon_initialized.eq.-1) then
      ! Continuing after finalizelast (init==-1), append output after previous contents
      open(unit=lout, file=fname, action='write', access='append', err=998)
   else
      ! Overwrite in case of new run (init==-3, -2)
      open(unit=lout, file=fname, action='write', err=998)
   endif
   out_open = 1
   if (ltmpfile) write(37,*) 'open ok'

   ! If succesful: also activate output to the .out-file

   if (ltmpfile) then
      write(37,*) 'opening streams-2=',l_outfil, l_screen, l_simpck
   endif
   call set_print_streams(nw_outfil=l_outfil, nw_screen=l_screen, nw_simpck=l_simpck)
   goto 999

   ! error handling: file open failed

998  continue
      out_open = -1
      tmpbuf = fname ! avoid potential overflow of bufout
      if (ltmpfile) write(37,*) '998: could not open file for writing'
      write(bufout,'(a,/,3x,a,/,3x,3a,/)') '** ERROR in Kalker Contact Result Element:', &
             'Could not open file for writing:', '"',trim(tmpbuf),'"'
      call write_log(4, bufout)
999 continue

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   if (ltmpfile) write(37,*) 'writing to stdout...'

   ! Debugging of license check, cwd, executable path,...

   if (.false.) then
      ilen = 256
      call getMostRecentDirectory(c_recent, ilen)
      call c_to_f_string(c_recent, f_recent, ilen)
      write(bufout,*) 'mostRecent: "',trim(f_recent),'", len=',ilen
      call write_log(1, bufout)
   endif

   ! Check if there is a license present for this version of CONTACT.

   my_license%is_valid = .true.
   my_license%licensee_name = 'Open source version'

   ! Include file with version identification, display to log-streams
   ! Display information from the licensefile, if it is valid.

   ! write(bufout,*) 'LicensePath=',trim(licpath)
   ! call write_log(1, bufout)

#  include 'VERSION'
   call WriteVersion(num_version, version, my_license, show_lic_error)

   ! Initialize the indirection table for result elements

   do jre = 1, MAX_REid_id
      ix_reid(jre) = -1
   enddo

   ! Initialize the module number per jre

   ire_module(1:MAX_NUM_REs) = -1

   ! Initialize the wheel-rail track data per jre

   do jre = 1, MAX_NUM_REs
      nullify(allwtd(jre)%wtd)
   enddo

   ! Initialize the contact data per jre+jcp

   do jre = 1, MAX_NUM_REs
      do jcp = 1, MAX_CPids
         nullify(allgds(jre,jcp)%gd)
      enddo
   enddo

   ! Initialize the active wtd and contact data pointers

   nullify(wtd)
   nullify(gd)

   ! Determine the maximum number of threads available, used a.o. to define timers;
   ! Start with OpenMP disabled, using 1 thread

#ifdef _OPENMP
   orig_mxthrd = omp_get_max_threads()
   call omp_set_num_threads(1)
   if (debug_openmp) then
      write(bufout,'(a,i3)') ' According to OpenMP, max #threads =', orig_mxthrd
      call write_log(1, bufout)
      call write_log(' OpenMP is disabled initially, setting max #threads = 1')
   endif
#else
   orig_mxthrd = 1
#endif

   ! Initialize CONTACT's timing mechanism

   louttm = -1
   if (out_open) louttm = lout
   call timers_contact_init(louttm=louttm, idebug=timer_output, maxerr=-1,                              &
                            naddon=numtmr, i0addon=ioftmr, mxthrd=orig_mxthrd)

   call timer_start(itimer_main)

   caddon_initialized = 1

   if (ltmpfile) close(37)

   ! return the version number and an error flag

   if (.not.my_license%is_valid) then
      ierror = CNTC_err_allow
   else
      ierror = 0
   endif
   if (idebug.ge.3) then
      write(bufout,'(2a,i6)') trim(pfx_str(subnam,-1,-1)),' version=',ifcver
      call write_log(1, bufout)
   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)

end subroutine cntc_initializeFirst

!------------------------------------------------------------------------------------------------------------

subroutine cntc_initializeFirst_new(ifcver, ierror, ioutput, c_wrkdir, c_outdir, c_expnam, len_wrkdir,  &
                                    len_outdir, len_expnam) &
   bind(c,name=CNAME_(cntc_initializefirst_new))
!--function: initialize the addon internal data structures and initialize output channels,
!            print version information; return the addon version number.
!            Old version provided for backward compatibility. Interface will be changed in v25.1.
!--category: 0, "m=any, glob":      not related to contact's modules 1 or 3, working on global data
   implicit none
!--subroutine arguments:
   integer,                intent(out)   :: ifcver       ! version of the CONTACT add-on
   integer,                intent(inout) :: ierror       ! error flag
   integer,                intent(in)    :: ioutput      ! output channels: 0 = out-file, 1 = out-file+screen
   character(kind=C_CHAR), intent(in)    :: c_wrkdir(*)  ! C-string: effective working folder
   character(kind=C_CHAR), intent(in)    :: c_outdir(*)  ! C-string: output folder
   character(kind=C_CHAR), intent(in)    :: c_expnam(*)  ! C-string: experiment name
   integer,                intent(in)    :: len_wrkdir   ! length of C-string
   integer,                intent(in)    :: len_outdir   ! length of C-string
   integer,                intent(in)    :: len_expnam   ! length of C-string
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_initializeFirst_new
#endif

   call write_log(' initializeFirst_new: calling initializeFirst...')
   call cntc_initializeFirst(ifcver, ierror, ioutput, c_wrkdir, c_outdir, c_expnam, len_wrkdir,         &
                        len_outdir, len_expnam)

end subroutine cntc_initializeFirst_new

!------------------------------------------------------------------------------------------------------------

subroutine cntc_initialize(ire, imodul, ifcver, ierror, c_outdir, len_outdir) &
   bind(c,name=CNAME_(cntc_initialize))
!--function: upon first call: perform one-time initializations;
!            for each ire: initialize and return the addon version number.
!--category: 0, "m=any, glob":      not related to contact's modules 1 or 3, working on global data
   implicit none
!--subroutine arguments:
   integer,                intent(in)    :: ire          ! result element ID
   integer,                intent(in)    :: imodul       ! module number 1=w/r contact, 3=basic contact
   integer,                intent(out)   :: ifcver       ! version of the CONTACT add-on
   integer,                intent(inout) :: ierror       ! error flag
   character(kind=C_CHAR), intent(in)    :: c_outdir(*)  ! C-string: output folder
   integer,                intent(in)    :: len_outdir   ! length of C-string
!--local variables:
   character(len=*), parameter    :: subnam = 'cntc_initialize'
   character(len= 30,kind=C_CHAR) :: c_string
   character(len=256) :: namtim
   integer            :: jcp, itimer
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_initialize
#endif

#ifdef _OPENMP
   if (omp_in_parallel().ne.0) then
      call write_log('ERROR: cntc_initialize may not be called from a parallel region')
      return
   endif
#endif

   ! Upon first call: perform one-time-only initializations
   !    this sets the version number and an error flag

   if (caddon_initialized.le.0) then
      c_string   = ' ' // C_NULL_CHAR
      call cntc_initializeFirst(ifcver, ierror, 0, c_string, c_outdir, c_string, 1, len_outdir, 1)
   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)

   ! Check the result element and module number

   if (ire.le.0 .or. ire.gt.MAX_REid_id) then
      write(bufout,'(a,i6,a,i6,a)') 'ERROR: Initialize: ire=',ire,' ID too large [1,',MAX_REid_id,      &
                '], aborting.'
      call write_log(1, bufout)
      call abort_run()
   endif

   if (imodul.ne.1 .and. imodul.ne.3) then
      write(bufout,'(a,i6,a,i6,a)') 'ERROR: Initialize: ire=',ire,' module number (',imodul,            &
                ') must be 1 or 3, aborting.'
      call write_log(1, bufout)
      call abort_run()
   endif

   ! Increment the number of active ire's (if ire is not initialized yet)

   if (ix_reid(ire).ge.0) then
      write(bufout,'(a,i6,a,i6,a)') ' Warning: Initialize: ire=',ire,' already initialized.'
      call write_log(1, bufout)
   else
      num_reids = num_reids + 1
      ix_reid(ire) = num_reids

      ! store module number for ire

      ire_module(ire) = imodul

      ! initialize timers for ire

      itimer = cntc_timer_num(ire, -1)
      write(namtim,'(a,i3)') 'Result element',ire
      call timer_name(itimer, namtmr=namtim)

      do jcp = 1, MAX_CPids
         itimer = cntc_timer_num(ire, jcp)
         write(namtim,'(a,i2,a,i3)') 'Cp.',jcp,' on Res.el.',ire
         call timer_name(itimer, namtmr=namtim)
      enddo

      if (idebug.ge.3) then
         write(bufout,'(2a,i3)') trim(pfx_str(subnam,-1,-1)),' number of active result elements=',num_reids
         call write_log(1, bufout)
      endif
   endif

   ! return the version number and an error flag

   ifcver = caddon_version
   if (.not.my_license%is_valid) then
      ierror = CNTC_err_allow
   else
      ierror = 0
   endif

   if (idebug.ge.3) then
      write(bufout,'(2a,i6)') trim(pfx_str(subnam,ire,-1)),' version=',ifcver
      call write_log(1, bufout)
   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)

end subroutine cntc_initialize

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setGlobalFlags(lenflg, params, values) &
   bind(c,name=CNAME_(cntc_setglobalflags))
!--function: used for configuring various flags that are the same for all contact problems
!--category: 0, "m=any, glob":      not related to contact's modules 1 or 3, working on global data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: lenflg         ! length of params/values arrays
   integer,      intent(in) :: params(lenflg) ! codes of the parameters to be communicated to CONTACT
   integer,      intent(in) :: values(lenflg) ! values of the parameters to be communicated to CONTACT
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_setGlobalFlags'
   integer  :: i                              ! loop counter
#ifdef _OPENMP
   integer  :: lic_mxthrd, cur_mxthrd, numthrd ! number of threads available/requested by user
#endif
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setGlobalFlags
#endif

#ifdef _OPENMP
   if (omp_in_parallel().ne.0) then
      call write_log('ERROR: cntc_setGlobalFlags may not be called from a parallel region')
      return
   endif
#endif

   ! Note: this subroutine may be called before cntc_initializeFirst.
   !       In that case, print-output will go nowhere.

   if (caddon_initialized.ge.0 .and. idebug.ge.4) call cntc_log_start(subnam, .true.)

   do i = 1, lenflg
      if (caddon_initialized.ge.0 .and. idebug.ge.3) then
         write(bufout,'(i6,2(a,i5))') i,'-th flag: code=',params(i),', new value=',values(i)
         call write_log(1, bufout)
      endif
      if (params(i).eq.0) then

         ! do nothing, allow arrays params/values to contain unused positions

      elseif (params(i).eq.CNTC_if_idebug) then

         ! set idebug-flag in the CONTACT add-on-data

         idebug       = values(i)

      elseif (params(i).eq.CNTC_if_licdbg) then

         ! set amount of print-output on the license used

         call LicenseSetDebug( values(i) )

      elseif (params(i).eq.CNTC_if_openmp) then

         if (caddon_initialized.ge.0) then

#ifdef _OPENMP

            ! activate/deactivate use of OpenMP multi-threading in CONTACT

            !  - update orig_mxthrd if the user increased the #threads himself

            call cntc_getmaxnumthreads(lic_mxthrd)
            cur_mxthrd = omp_get_max_threads()

            if (cur_mxthrd.gt.orig_mxthrd) then
               if (debug_openmp) then
                  write(bufout,'(a,i3)') ' According to OpenMP, max #threads increased to', orig_mxthrd
                  call write_log(1, bufout)
               endif
               orig_mxthrd = cur_mxthrd
            endif

            !  - restrict the max #threads according to the license

            if (orig_mxthrd.gt.lic_mxthrd) then
               if (debug_openmp .or. idebug.ge.3) then
                  if (lic_mxthrd.le.1) then
                     call write_log(' Warning: no license for OpenMP, setting max. #threads=1')
                  else
                     write(bufout,'(2(a,i3))') ' Warning: no license for',orig_mxthrd,                  &
                                                         ' threads, setting max. #threads=',lic_mxthrd
                     call write_log(1, bufout)
                  endif
               endif
               orig_mxthrd = lic_mxthrd
            endif

            !  - set the actual number of threads according to the parameter value

            numthrd = values(i)

            if (numthrd.le.-1) then

               if (debug_openmp .or. idebug.ge.3) then
                  write(bufout,'(a,i4)') ' cntc_setGlobalFlags: set OpenMP #threads = AUTOMATIC =',     &
                                                                                        orig_mxthrd
                  call write_log(1, bufout)
               endif
               call omp_set_num_threads( orig_mxthrd )

            elseif (numthrd.le.1) then

               if (debug_openmp .or. idebug.ge.3) then
                  call write_log(' cntc_setGlobalFlags: set OpenMP #threads = 1')
               endif
               call omp_set_num_threads( 1 )

            else

               if (debug_openmp .or. idebug.ge.3) then
                  write(bufout,'(2(a,i3),a)') ' cntc_setGlobalFlags: set OpenMP #threads = min(',       &
                                                                        orig_mxthrd,',',numthrd,')'
                  call write_log(1, bufout)
               endif
               call omp_set_num_threads( min(numthrd, orig_mxthrd) )

            endif

#endif
         endif ! caddon_initialized

      else

         if (caddon_initialized.ge.0) then
            write(bufout,'(a,i8,a)') 'WARNING: unknown parameter code=',params(i),' is ignored.'
            call write_log(1, bufout)
         endif

      endif

   enddo

   if (caddon_initialized.ge.0 .and. idebug.ge.2) then
      do i = 1, lenflg
         if (params(i).ne.0) then
            write(bufout,'(2a,i2,3a,i8)') trim(pfx_str(subnam,-1,-1)),' i=',i, ', set flag ',           &
                cntc_flagName(params(i)),' to value',values(i)
            call write_log(1, bufout)
         endif
      enddo
   endif

   if (caddon_initialized.ge.0 .and. idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setGlobalFlags

!------------------------------------------------------------------------------------------------------------

subroutine cntc_readInpFile(ire, inp_type, c_fname, len_fname, ierror) &
   bind(c,name=CNAME_(cntc_readinpfile))
!--function: used for configuring multiple data items at once using the inp-file format
!--category: 2, "m=1 only, wtd":    available for module 1 only, working on wtd data
   implicit none
!--subroutine arguments:
   integer,                intent(in)  :: ire        ! result element ID
   integer,                intent(in)  :: inp_type   ! select one of available input-routines
   character(kind=C_CHAR), intent(in)  :: c_fname(*) ! C-string: full path of inp-file
   integer,                intent(in)  :: len_fname  ! length of C-string
   integer,                intent(out) :: ierror     ! error code
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_readInpFile'
   character(len=256)       :: f_fname
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_readInpFile
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, -1, 1, -1, subnam, ierror)
   if (ierror.lt.0) return

   ierror = 0

   call c_to_f_string(c_fname, f_fname, len_fname)

   ! read contents of input-file dependent on type of input

   if (inp_type.eq.CNTC_inp_modul1) then

      call wr_input_modul1(f_fname, wtd, ierror)
      if (ierror.ne.0) ierror = CNTC_err_input

   elseif (inp_type.eq.CNTC_inp_spck) then

      call wr_input_spck(f_fname, wtd, ierror)
      if (ierror.ne.0) ierror = CNTC_err_input

   else

      ierror = CNTC_err_input
      write(bufout,'(a,i3,a,i8,a)') ' ERROR(',ierror,'): unknown inp-file code=',inp_type,' is ignored.'
      call write_log(1, bufout)

   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_readInpFile

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setFlags(ire, icp, lenflg, params, values) &
   bind(c,name=CNAME_(cntc_setflags))
!--function: used for configuring various flags for a contact problem
!--category: 5, "m=any, wtd":       available for modules 1 and 3, in module 1 working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: icp            ! contact problem ID
   integer,      intent(in) :: lenflg         ! length of params/values arrays
   integer,      intent(in) :: params(lenflg) ! codes of the parameters to be communicated to CONTACT
   integer,      intent(in) :: values(lenflg) ! values of the parameters to be communicated to CONTACT
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_setFlags'
   integer  :: i, imodul, ierror, psflcin
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setFlags
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, -1, subnam, ierror)
   if (ierror.lt.0) return

   imodul = ire_module(ix_reid(ire))

   do i = 1, lenflg

      if (params(i).eq.0) then

         ! do nothing, allow arrays params/values to contain unused positions

      elseif (params(i).eq.CNTC_if_units) then          ! set the unit convention

         call cntc_select_units(my_scl, values(i))

      elseif (params(i).eq.CNTC_ic_config) then         ! set C1-digit in the RE-CP-data

         my_ic%config = max(0, min(5, values(i)))

      elseif (params(i).eq.CNTC_ic_pvtime) then         ! set P-digit in the RE-CP-data

         my_ic%pvtime = max(0, min(3, values(i)))

      elseif (params(i).eq.CNTC_ic_bound) then          ! set B-digit in the RE-CP-data

         my_ic%bound  = max(0, min(6, values(i)))

      elseif (params(i).eq.CNTC_ic_tang) then           ! set T-digit in the RE-CP-data

         my_ic%tang   = max(0, min(3, values(i)))

      elseif (params(i).eq.CNTC_ic_norm) then           ! set N-digit in the RE-CP-data

         my_ic%norm   = max(0, min(1, values(i)))

      elseif (params(i).eq.CNTC_ic_force) then          ! set F-digit in the RE-CP-data

         if (imodul.eq.3) then
            my_ic%force1  = 0
            my_ic%force3  = max(0, min(2, values(i)))
         else
            my_ic%force1  = max(0, min(3, values(i)))
            my_ic%force3  = 0
         endif

      elseif (params(i).eq.CNTC_ic_sens) then           ! set S2-digit in the RE-CP-data

         my_ic%sens   = max(0, min(3, values(i)))

      elseif (params(i).eq.CNTC_ic_discns) then 

         ! set D-digit in the RE-CP-data (only makes sense for module 1)

         if (values(i).ge.0 .and. values(i).le.9) my_ic%discns1_inp = values(i)
         if (values(i).ge.2 .and. values(i).le.9) my_ic%discns1_eff = values(i)

      elseif (params(i).eq.CNTC_ic_npomax) then         ! set max. #elements in potential contact (mod. 1)

         if (values(i).lt.100 .or. values(i).gt.1000000) then
            write(bufout,'(a,i9,a)') ' Warning: NPOMAX =',values(i),                                    &
                                                        ' is adjusted to be >= 100 and <= 1000000.'
            call write_log(1, bufout)
         endif
         wtd%discr%npot_max = min(1000000, max(100, values(i)))

      elseif (params(i).eq.CNTC_ic_inflcf) then         ! set C3-digit in the RE-CP-data

         if (values(i).ge.0 .and. values(i).le.4) my_ic%gencr_inp = values(i)
         if (values(i).ge.2 .and. values(i).le.4) my_mater%gencr_eff = values(i)

      elseif (params(i).eq.CNTC_ic_ifmeth) then         ! set IF_METH for Blanco-approach (C3=4)

         my_mater%if_meth = max(0, min(1, values(i)))

      elseif (params(i).eq.CNTC_ic_ifvari) then         ! set IF_VER for Blanco-approach (C3=4)

         my_mater%if_ver  = max(1, min(4, values(i)))

      elseif (params(i).eq.CNTC_ic_mater) then          ! set M-digit in the RE-CP-data

         call write_log(' ERROR: the M-digit cannot be set by cntc_setflags (anymore).')
         call write_log('        use cntc_setmaterialparameters instead.')
         ! my_ic%mater  = max(0, min(5, values(i)))

      elseif (params(i).eq.CNTC_ic_exrhs) then         ! set E-digit in the RE-CP-data

         my_ic%rztang = max(0, min(3, values(i)))

      elseif (params(i).eq.CNTC_ic_iestim) then         ! set I-digit in the RE-CP-data

         my_ic%iestim = max(0, min(3, values(i)))

      elseif (params(i).eq.CNTC_ic_matfil) then         ! set A-digit in the RE-CP-data

         my_ic%matfil_surf = values(i)
         if (values(i).lt.0 .or. values(i).gt.2) my_ic%matfil_surf = 0

      elseif (params(i).eq.CNTC_ic_output) then         ! set O-digit in the RE-CP-data

         my_ic%output_surf = max(0, min(5, values(i)))

      elseif (params(i).eq.CNTC_ic_sbsfil) then         ! set A_s-digit in the RE-CP-data

         my_ic%matfil_subs = values(i)
         if (values(i).lt.0 .or. values(i).gt.2) my_ic%matfil_subs = 0

      elseif (params(i).eq.CNTC_ic_sbsout) then         ! set O_s-digit in the RE-CP-data

         my_ic%output_subs = max(0, min(4, values(i)))
         if (my_ic%output_subs.eq.3) my_ic%output_subs = 2

      elseif (params(i).eq.CNTC_ic_flow) then           ! set W-digit in the RE-CP-data

         my_ic%flow   = max(0, min(9, values(i)))

      elseif (params(i).eq.CNTC_ic_xflow) then          ! set X-codes in the RE-CP-data

         ! note: usage is different from other flags; we could introduce a subroutine setdebugflags

         my_ic%xflow  = 1
         psflcin      = max(0, values(i))
         call ic_unpack_dbg(imodul, psflcin, my_ic)

      elseif (params(i).eq.CNTC_ic_return) then         ! set R-digit in the RE-CP-data

         my_ic%return = max(0, min(3, values(i)))

      elseif (params(i).eq.CNTC_if_wrtinp) then         ! set flag wrtinp in the RE-CP-data

         my_ic%wrtinp = values(i)

      elseif (params(i).eq.CNTC_if_idebug) then         ! set if_idebug using cntc_setGlobalFlags

         call cntc_setGlobalFlags(1, params(i), values(i))

      elseif (params(i).eq.CNTC_if_licdbg) then         ! set if_licdbg using cntc_setGlobalFlags

         call cntc_setGlobalFlags(1, params(i), values(i))

      elseif (params(i).eq.CNTC_if_openmp) then         ! set if_openmp using cntc_setGlobalFlags

         call cntc_setGlobalFlags(1, params(i), values(i))

      elseif (params(i).eq.CNTC_if_timers) then

         call timers_contact_outlevel(values(i))

      else

         write(bufout,'(a,i8,a)') ' WARNING: unknown parameter code=',params(i),' is ignored.'
         call write_log(1, bufout)

      endif
   enddo

   if (idebug.ge.2) then
      do i = 1, lenflg
         if (params(i).ne.0) then
            if (params(i).eq.CNTC_if_units) then
               write(bufout,'(2a,i2,4a)') trim(pfx_str(subnam,ire,icp)), ' i=',i,', set flag ',         &
                        cntc_flagName(params(i)),' to value ', trim(cntc_flagName(values(i)))
            else
               write(bufout,'(2a,i2,3a,i8)') trim(pfx_str(subnam,ire,icp)), ' i=',i,', set flag ',      &
                        cntc_flagName(params(i)),' to value',values(i)
            endif
            call write_log(1, bufout)
         endif
      enddo
   endif
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setFlags

!------------------------------------------------------------------------------------------------------------

function cntc_flagName(iflag)
!--purpose: print the name of a flag
   implicit none
!--subroutine arguments:
   integer,          intent(in)  :: iflag
!--return value:
   character(len=14)  :: cntc_flagName

   if (iflag.eq.CNTC_if_units) then
      cntc_flagName = 'CNTC_if_units'
   elseif (iflag.eq.CNTC_un_cntc) then
      cntc_flagName = 'CNTC_un_cntc'
   elseif (iflag.eq.CNTC_un_spck) then
      cntc_flagName = 'CNTC_un_spck'
   elseif (iflag.eq.CNTC_un_si) then
      cntc_flagName = 'CNTC_un_si'
   elseif (iflag.eq.CNTC_un_imper) then
      cntc_flagName = 'CNTC_un_imper'
   elseif (iflag.eq.CNTC_ic_config) then
      cntc_flagName = 'CNTC_ic_config'
   elseif (iflag.eq.CNTC_ic_pvtime) then
      cntc_flagName = 'CNTC_ic_pvtime'
   elseif (iflag.eq.CNTC_ic_bound) then
      cntc_flagName = 'CNTC_ic_bound'
   elseif (iflag.eq.CNTC_ic_tang) then
      cntc_flagName = 'CNTC_ic_tang'
   elseif (iflag.eq.CNTC_ic_norm) then
      cntc_flagName = 'CNTC_ic_norm'
   elseif (iflag.eq.CNTC_ic_force) then
      cntc_flagName = 'CNTC_ic_force'
   elseif (iflag.eq.CNTC_ic_discns) then
      cntc_flagName = 'CNTC_ic_discns'
   elseif (iflag.eq.CNTC_ic_inflcf) then
      cntc_flagName = 'CNTC_ic_inflcf'
   elseif (iflag.eq.CNTC_ic_mater) then
      cntc_flagName = 'CNTC_ic_mater'
   elseif (iflag.eq.CNTC_ic_xflow) then
      cntc_flagName = 'CNTC_ic_xflow'
   elseif (iflag.eq.CNTC_ic_iestim) then
      cntc_flagName = 'CNTC_ic_iestim'
   elseif (iflag.eq.CNTC_ic_output) then
      cntc_flagName = 'CNTC_ic_output'
   elseif (iflag.eq.CNTC_ic_flow) then
      cntc_flagName = 'CNTC_ic_flow'
   elseif (iflag.eq.CNTC_ic_return) then
      cntc_flagName = 'CNTC_ic_return'
   elseif (iflag.eq.CNTC_ic_matfil) then
      cntc_flagName = 'CNTC_ic_matfil'
   elseif (iflag.eq.CNTC_ic_sbsout) then
      cntc_flagName = 'CNTC_ic_sbsout'
   elseif (iflag.eq.CNTC_ic_sbsfil) then
      cntc_flagName = 'CNTC_ic_sbsfil'
   elseif (iflag.eq.CNTC_ic_sens) then
      cntc_flagName = 'CNTC_ic_sens'
   elseif (iflag.eq.CNTC_ic_ifmeth) then
      cntc_flagName = 'CNTC_ic_ifmeth'
   elseif (iflag.eq.CNTC_ic_ifvari) then
      cntc_flagName = 'CNTC_ic_ifvari'
   elseif (iflag.eq.CNTC_ic_npomax) then
      cntc_flagName = 'CNTC_ic_npomax'
   elseif (iflag.eq.CNTC_if_idebug) then
      cntc_flagName = 'CNTC_if_idebug'
   elseif (iflag.eq.CNTC_if_licdbg) then
      cntc_flagName = 'CNTC_if_licdbg'
   elseif (iflag.eq.CNTC_if_wrtinp) then
      cntc_flagName = 'CNTC_if_wrtinp'
   elseif (iflag.eq.CNTC_if_openmp) then
      cntc_flagName = 'CNTC_if_openmp'
   elseif (iflag.eq.CNTC_if_timers) then
      cntc_flagName = 'CNTC_if_timers'
   else
      write(cntc_flagName,'(i8)') iflag
   endif

end function cntc_FlagName

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setMetadata(ire, icp, lenmta, params, values) &
   bind(c,name=CNAME_(cntc_setmetadata))
!--function: used for providing various metadata for a contact problem
!--category: 5, "m=any, wtd":       available for modules 1 and 3, in module 1 working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: icp            ! contact problem ID
   integer,      intent(in) :: lenmta         ! length of params/values arrays
   integer,      intent(in) :: params(lenmta) ! codes of the metadata to be communicated to CONTACT
   real(kind=8), intent(in) :: values(lenmta) ! values of the metadata to be communicated to CONTACT
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_setMetadata'
   integer  :: i, ierror
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setMetadata
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, -1, subnam, ierror)
   if (ierror.lt.0) return

   do i = 1, lenmta
      if (idebug.ge.3) then
         write(bufout,'(i6,a,i5,a,g12.4)') i,'-th flag: code=',params(i),', value=',values(i)
         call write_log(1, bufout)
      endif
      if (params(i).eq.0) then
         ! do nothing, allow arrays params/values to contain unused positions
      elseif (params(i).eq.CNTC_mt_tim) then
         ! set the simulation time [s]
         my_meta%tim = values(i)
      elseif (params(i).eq.CNTC_mt_sw) then
         ! set s_ws [length], location of wheel-set CM along the track curve, converted to [mm]
         my_meta%s_ws = my_scl%len * values(i)
      elseif (params(i).eq.CNTC_mt_xr) then
         ! TODO: the SIMPACK-interface is outdated !!!
         ! [xyz]_w [mm]    location of right wheel origin m_w in terms of track coordinates
         ! [yz]_r  [mm]    location of right rail origin m_r in terms of track coordinates
         ! [xyz]_ref [mm]  location of contact reference point m_ref in terms of right rail coordinates
         my_meta%xcp_r = my_scl%len * values(i)
      elseif (params(i).eq.CNTC_mt_yr) then
         my_meta%ycp_r = my_scl%len * values(i)
      elseif (params(i).eq.CNTC_mt_xw) then
         my_meta%xcp_w = my_scl%len * values(i)
      elseif (params(i).eq.CNTC_mt_yw) then
         my_meta%ycp_w = my_scl%len * values(i)
      elseif (params(i).eq.CNTC_mt_run) then
         my_meta%irun   = nint(values(i))
      elseif (params(i).eq.CNTC_mt_axle) then
         my_meta%iax    = nint(values(i))
      elseif (params(i).eq.CNTC_mt_side) then
         my_meta%iside  = nint(values(i))
      else
         write(bufout,'(a,i8,a)') 'WARNING: unknown metadata code=',params(i),' is ignored.'
         call write_log(1, bufout)
      endif
   enddo

   if (idebug.ge.2) then
      do i = 1, lenmta
         if (params(i).ne.0) then
            write(bufout,'(2a,i2,a,i5,a,g12.4)') trim(pfx_str(subnam,ire,icp)),' i=',i,', set metadata', &
                params(i),' to value',values(i)
            call write_log(1, bufout)
         endif
      enddo
   endif
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setMetadata

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setSolverFlags(ire, icp, gdigit, nints, iparam, nreals, rparam) &
   bind(c,name=CNAME_(cntc_setsolverflags))
!--function: set parameters for the iterative solution algorithms
!    0: default solvers,    iparam = [maxgs, maxin, maxnr, maxout],         rparam = [eps]
!    1: keep parameters,    iparam = [ ],                                   rparam = [ ]
!    2: always ConvexGS,    iparam = [maxgs, maxin, maxnr, maxout, inislp], 
!                                                   rparam = [eps, omegah, omegas, omgslp]
!    3: SteadyGS if poss.,  iparam = [maxgs, maxin, maxnr, maxout, inislp], 
!                                                   rparam = [eps, omegah, omegas, omgslp]
!    4: flags veloc.dep.    iparam = [maxgs, maxin, maxnr, maxout, inislp], rparam = [eps, omgslp]
!    5: GDsteady if poss.,  iparam = [maxgs, maxin, maxnr, maxout, kdowfb], rparam = [eps, fdecay, 
!                                                       betath, d_ifc, d_lin, d_cns, d_slp, pow_s]
!    6: flags sensitivities iparam = [mxsens],                              rparam = [epsens]
!--category: 5, "m=any, wtd":       available for modules 1 and 3, in module 1 working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: icp            ! contact problem ID
   integer,      intent(in) :: gdigit         ! G-digit, selecting the iteration methods to be used
   integer,      intent(in) :: nints, nreals  ! number of integer/real parameters provided
   integer,      intent(in) :: iparam(*)      ! array of integer parameters
   real(kind=8), intent(in) :: rparam(*)      ! array of real parameters
!--local variables:
   integer, parameter  :: nints_loc(0:6)  = (/ 4, 0, 5, 5, 5, 5, 1 /)
   integer, parameter  :: nreals_loc(0:6) = (/ 1, 0, 4, 4, 2, 8, 1 /)
   integer             :: ierror
   character(len=*), parameter :: subnam = 'cntc_setSolverFlags'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setSolverFlags
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, -1, subnam, ierror)
   if (ierror.lt.0) return

   if (gdigit.lt.0 .or. gdigit.gt.6) then
      write(bufout,'(2a,i2,a)') trim(pfx_str(subnam,ire,icp)),' method G=',gdigit,' does not exist.'
      call write_log(1, bufout)
      return
   endif

   if (nints.ne.nints_loc(gdigit) .or. nreals.ne.nreals_loc(gdigit)) then
      write(bufout,'(a,3(a,i2),a)') trim(pfx_str(subnam,ire,icp)),' method G=',gdigit,' needs',        &
                nints_loc(gdigit),' integer and',nreals_loc(gdigit),' real parameters.'
      call write_log(1, bufout)
      return
   endif

   ! store input-value of G-digit

   if (gdigit.ge.0 .and. gdigit.le.5) my_ic%gausei_inp = gdigit

   if (gdigit.ne.1 .and. gdigit.ne.6) then

      ! store effective value of G-digit

      my_solv%gausei_eff = gdigit

      ! set the solver parameters except when G=1 or 6

      my_solv%eps    = max(1d-20, rparam(1))
      my_solv%maxgs  = max(1,     iparam(1))
      my_solv%maxin  = max(1,     iparam(2))
      my_solv%maxnr  = max(1,     iparam(3))
      my_solv%maxout = max(1,     iparam(4))

      if (idebug.ge.2) then
         write(bufout,'(2a,i2,4(a,i6),a)') trim(pfx_str(subnam,ire,icp)),' G=',gdigit, ', MaxGS=',      &
             my_solv%maxgs,', MaxIN=',my_solv%maxin,', MaxNR=', my_solv%maxnr,', MaxOUT=',my_solv%maxout
         call write_log(1, bufout)
      endif

      ! additional parameters for G=2,3

      if (gdigit.eq.2 .or. gdigit.eq.3) then
         my_solv%inislp = iparam(5)
         my_solv%omegah = max(1d-20, rparam(2))
         my_solv%omegas = max(1d-20, rparam(3))
         my_solv%omgslp = max(1d-20, rparam(4))

         if (idebug.ge.2) then
            write(bufout,'(2a,i2,3(a,f7.4),a,i4,a)') trim(pfx_str(subnam,ire,icp)),' G=',gdigit,        &
                ', OmegaH=',my_solv%omegah,', OmegaS=',my_solv%omegas,', OmgSlp=',my_solv%omgslp,       &
                ', IniSlp=',my_solv%inislp
            call write_log(1, bufout)
         endif
      endif

      ! additional parameters for G=4

      if (gdigit.eq.4) then
         my_solv%inislp = iparam(5)
         my_solv%omgslp = max(1d-20, rparam(2))

         if (idebug.ge.2) then
            write(bufout,'(2a,i2,a,f7.4,a,i4,a)') trim(pfx_str(subnam,ire,icp)),' G=',gdigit,           &
                ', OmgSlp=', my_solv%omgslp, ', IniSlp=',my_solv%inislp
            call write_log(1, bufout)
         endif
      endif

      ! additional parameters for G=5

      if (gdigit.eq.5) then
         my_solv%fdecay = rparam(2)   ! <0.001 or >0.999 used for E_down / E_trl
         my_solv%betath = rparam(3)
         my_solv%kdowfb = iparam(5)
         my_solv%d_ifc  = max(0.01d0, rparam(4))
         my_solv%d_lin  = rparam(5)
         my_solv%d_cns  = max(0.01d0, rparam(6))
         my_solv%d_slp  = max(0.01d0, rparam(7))
         my_solv%pow_s  = max(0.01d0, min(10d0, rparam(8)))

         ! check parameters diagonal scaling

         if (my_solv%d_lin*(my_solv%d_cns-my_solv%d_ifc).lt.0d0) then
            call write_log('ERROR: D_LIN and D_CNS-D_IFC should have same sign.')
            my_solv%d_lin = 0d0
            my_solv%d_cns = my_solv%d_ifc
         endif

         ! set variant used for GDsteady

         if (my_solv%fdecay.ge.0.999d0) then
            my_solv%gd_meth = 1     ! E_trl
         elseif (my_solv%fdecay.le.0.001d0) then
            my_solv%gd_meth = 2     ! E_down(k)
            my_solv%kdown   = max(1, nint(-my_solv%fdecay))
         else
            my_solv%gd_meth = 3     ! E_keep(f)
         endif

         if (idebug.ge.2) then
            write(bufout,'(2a,i2,3(a,f7.4))') trim(pfx_str(subnam,ire,icp)),' G=',gdigit,', f=',        &
                my_solv%fdecay, ', d_ifc=',my_solv%d_ifc,', pow_s=',my_solv%pow_s
            call write_log(1, bufout)
         endif
      endif

   elseif (gdigit.eq.6) then

      ! gdigit=6: set solver parameters related to the calculation of sensitivities

      my_solv%mxsens = max(1, iparam(1))
      my_solv%epsens = max(1d-20, rparam(1))

      if (idebug.ge.2) then
         write(bufout,'(2a,2(i0,a),es8.1)') trim(pfx_str(subnam,ire,icp)),' G=',gdigit,', set MxSENS=', &
                my_solv%mxsens, ', Epsens=', my_solv%epsens
         call write_log(1, bufout)
      endif

   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setSolverFlags

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setMaterialParameters(ire, icp, imeth, nparam, rparam) &
   bind(c,name=CNAME_(cntc_setmaterialparameters))
!--function: set the M-digit and material parameters for a contact problem
! values < 10 are used to configure the M-digit
!    0: purely elastic material,           params = [ nu1, nu2, g1, g2 ]
!    1: visco-elastic material,            params = [ nu1, nu2, g1, g2, fg1, fg2, vt1, vt2 ]
!    2: modified Fastsim, 1 flexibility    params = [ nu1, nu2, g1, g2, flx, k0_mf, alfamf, betamf ]
!    3: modified Fastsim, 3 flexibilities  params = [ nu1, nu2, g1, g2, k0_mf, alfamf, betamf ]
!    4: elastic + elasto-plastic 3rd body  params = [ nu1, nu2, g1, g2, g3, laythk, tau_c0, k_tau ]
!    5: modified FaStrip                   params = [ nu1, nu2, g1, g2, k0_mf, alfamf, betamf ]
! values >= 10 are used to configure the M2-digit, M2 = imeth - 10
!   12: force-based proportional damping   params = [ cdampn, cdampt, dfnmax, dftmax ]
!
! dimensions: nu1, nu2        [-],        g1, g2    [force/area],   
!             fg1, fg2        [-],        vt1, vt2  [time],
!             flx         [volume/force], k0_mf, alfamf, betamf [-], 
!             g3, tau_c0   [force/area],  laythk    [length],        k_tau [force/volume]
!             cdampn, cdampt  [time],     dfnmax, dftmax   [force/time]
!--category: 5, "m=any, wtd":       available for modules 1 and 3, in module 1 working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: icp            ! contact problem ID
   integer,      intent(in) :: imeth          ! M-digit
   integer,      intent(in) :: nparam         ! number of parameters in rparam
   real(kind=8), intent(in) :: rparam(nparam) ! material parameters, dependent on m-digit
!--local variables:
   integer, parameter  :: nparam_loc(0:12) = (/ 4, 8, 8, 7, 8, 7, 0, 0, 0, 0, 0, 0, 4 /)
   logical             :: set_m_digit, set_m2_digit
   integer             :: ierror
   character(len=*), parameter :: subnam = 'cntc_setMaterialParameters'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setMaterialParameters
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, -1, subnam, ierror)
   if (ierror.lt.0) return

   ! check value of M-digit. 

   set_m_digit  = (imeth.ge.0 .and. imeth.le.5)
   set_m2_digit = (imeth.ge.10 .and. imeth.le.12)
   if (.not.set_m_digit .and. .not.set_m2_digit) then
      write(bufout,'(2a,i4,a)') trim(pfx_str(subnam,ire,icp)), ' method',imeth,' does not exist.'
      call write_log(1, bufout)
      return
   endif

   ! check number of parameters provided

   if ((nparam_loc(imeth).eq.0 .and. nparam.ge.2) .or.                                                  &
       (nparam_loc(imeth).ne.0 .and. nparam.ne.nparam_loc(imeth))) then
      write(bufout,'(2a,2(i0,a))') trim(pfx_str(subnam,ire,icp)),' method ',imeth,' needs ',            &
                nparam_loc(imeth),' parameters.'
      call write_log(1, bufout)
      return
   endif

   if (set_m2_digit) then

      ! configuring M2-digit

      my_ic%mater2    = imeth - 10

      if (imeth.eq.10) then       

         ! M2 = 0 - no damping

         my_mater%cdampn = 0d0
         my_mater%cdampt = 0d0
         my_mater%dfnmax = 0d0
         my_mater%dftmax = 0d0

      elseif (imeth.eq.12) then

         ! M2 = 2 - force-based proportional damping

         my_mater%cdampn = rparam(1)                ! [s]
         my_mater%cdampt = rparam(2)                ! [s]
         my_mater%dfnmax = rparam(3) * my_scl%forc  ! [N/s]
         my_mater%dftmax = rparam(4) * my_scl%forc  ! [N/s]

      endif

   else

      ! configuring M-digit

      my_ic%mater        = imeth
      my_mater%mater_eff = my_ic%mater

      ! all methods: params(1--4) == nu1, nu2, gg1, gg2

      my_mater%poiss(1) = rparam(1) ! [-]
      my_mater%poiss(2) = rparam(2)
      my_mater%gg(1)    = rparam(3) * my_scl%forc / my_scl%area
      my_mater%gg(2)    = rparam(4) * my_scl%forc / my_scl%area

      if (imeth.eq.0) then

         ! purely elastic contact - no further data

      elseif (imeth.eq.1) then

         ! visco-elastic contact - params(5--8) == fg1, fg2, vt1, vt2

         my_mater%fg(1)    = rparam(5)  ! [-]
         my_mater%fg(2)    = rparam(6)
         my_mater%tc(1)    = rparam(7)  ! [time]
         my_mater%tc(2)    = rparam(8)

      elseif (imeth.eq.2) then

         ! Fastsim, 1 flexiblity - params(5--8) == flx, k0_mf, alfamf, betamf

         my_mater%flx(1:3) = rparam(5) * my_scl%area * my_scl%len / my_scl%forc
         my_mater%k0_mf    = max(1d-6, rparam(6)) ! [-]
         my_mater%alfamf   = max(1d-6, rparam(7))
         my_mater%betamf   = max(1d-6, rparam(8))

         if (idebug.ge.2) then
            write(bufout,'(2a,i1,4(a,f7.4))') trim(pfx_str(subnam,ire,icp)), ' M=',my_ic%mater,         &
                   ', flx=',my_mater%flx(1),', k0_mf=', my_mater%k0_mf, ', alfamf=', my_mater%alfamf,   &
                   ', betamf=', my_mater%betamf
            call write_log(1, bufout)
         endif

      elseif (imeth.eq.3 .or. imeth.eq.5) then

         ! 3: Fastsim, 3 flexiblities - params(5--7) == k0_mf, alfamf, betamf
         ! 5: FaStrip                 - params(5--7) == k0_mf, alfamf, betamf

         my_mater%k0_mf    = max(1d-6, rparam(5)) ! [-]
         my_mater%alfamf   = max(1d-6, rparam(6))
         my_mater%betamf   = max(1d-6, rparam(7))

         if (idebug.ge.2) then
            write(bufout,'(2a,i1,3(a,f7.4))') trim(pfx_str(subnam,ire,icp)), ' M=',my_ic%mater,         &
               ', k0_mf=', my_mater%k0_mf, ', alfamf=', my_mater%alfamf, ', betamf=', my_mater%betamf
            call write_log(1, bufout)
         endif

      elseif (imeth.eq.4) then

         ! elastic contact + elasto-plastic interfacial layer - params(5--8) == gg3, laythk, tau_c0, k_tau

         my_mater%gg3      = max(1d-6, rparam(5) * my_scl%forc / my_scl%area)
         my_mater%laythk   = max(1d-9, rparam(6) * my_scl%len)
         my_mater%tau_c0   = rparam(7) * my_scl%forc / my_scl%area
         my_mater%k_tau    = rparam(8) * my_scl%forc / my_scl%area / my_scl%len
         if (my_mater%tau_c0.ge.1d10) my_mater%tau_c0 = 0d0

         if (idebug.ge.2) then
            write(bufout,'(2a,i1,a,f7.1,a,f6.3,a,/, 43x,a,g9.1,a,f7.1,a)') trim(pfx_str(subnam,ire,icp)),  &
                   ' M=', my_ic%mater,', Gg3=',my_mater%gg3,' [N/mm2], Laythk=',my_mater%laythk, ' [mm]',  &
                   '  Tau_c0=', my_mater%tau_c0,' [N/mm2], K_tau=', my_mater%k_tau,' [N/mm3]'
            call write_log(2, bufout)
         endif

         ! set flexibility for third body layer

         my_mater%flx(1:3) = my_mater%laythk / my_mater%gg3

      endif

      ! compute combined material constants

      call combin_mater(my_mater)

      ! update MAXOUT for quasi-identity or non-quasi-identity

      if (my_solv%maxout.ge.2 .and. abs(my_mater%ak).lt.1d-4) my_solv%maxout =  1
      if (my_solv%maxout.le.1 .and. abs(my_mater%ak).ge.1d-4) my_solv%maxout = 10

      if (idebug.ge.2) then
         write(bufout,'(2a,f7.1,2(a,f5.2),a)') trim(pfx_str(subnam,ire,icp)),' G=',my_mater%ga,         &
                   ' [N/mm2], Nu=',my_mater%nu,', K=',my_mater%ak, ' [-]'
         call write_log(1, bufout)
      endif

   endif ! imeth>=10 (M2-digit)

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setMaterialParameters

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setTemperatureData(ire, icp, imeth, nparam, params) &
   bind(c,name=CNAME_(cntc_settemperaturedata))
!--function: set H-digit and parameters for the temperature calculation for a contact problem
!    0: no temperature calculation,     params = []
!    1: keep old parameters,            params = []
!    3: calculate temperature based on new parameters and steady rolling,
!       params = [bktemp1, heatcp1, lambda1, dens1, bktemp2, heatcp2, lambda2, dens2, betapl]
!
! dimensions:  bktemp: [C],  heatcp: [J/kg-C],  lambda: [W/length-C],  dens: [kg/length^3],  betapl: [-]
!--category: 5, "m=any, wtd":       available for modules 1 and 3, in module 1 working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: icp            ! contact problem ID
   integer,      intent(in) :: imeth          ! H-digit
   integer,      intent(in) :: nparam         ! number of parameters provided
   real(kind=8), intent(in) :: params(nparam) ! depending on method that is used
!--local variables:
   integer, parameter  :: nparam_loc(0:3) = (/ 0, 0, 0, 9 /)
   integer      :: ierror
   character(len=*), parameter :: subnam = 'cntc_setTemperatureData'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setTemperatureData
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 0, subnam, ierror)
   if (ierror.lt.0) return

   ! Check H-digit, number of values provided

   if (imeth.lt.0 .or. imeth.eq.2 .or. imeth.gt.3) then
      write(bufout,'(2a,i4,a)') trim(pfx_str(subnam,ire,icp)),' method H=',imeth,' does not exist.'
      call write_log(1, bufout)
      return
   endif

   if ((nparam_loc(imeth).eq.0 .and. nparam.ge.2) .or.                                                  &
       (nparam_loc(imeth).ne.0 .and. nparam.ne.nparam_loc(imeth))) then
      write(bufout,'(2a,2(i1,a))') trim(pfx_str(subnam,ire,icp)),' method H=',imeth,' needs ',          &
                nparam_loc(imeth),' parameters.'
      call write_log(1, bufout)
      return
   endif

   ! Store H-digit

   my_ic%heat = imeth

   if (imeth.eq.3) then

      ! Coulomb friction with static/kinematic coefficients of friction

      my_mater%bktemp(1:2) = (/ params(1), params(5) /)
      my_mater%heatcp(1:2) = (/ params(2), params(6) /)
      my_mater%lambda(1:2) = (/ params(3), params(7) /) / my_scl%len
      my_mater%dens(1:2)   = (/ params(4), params(8) /) / my_scl%len**3
      my_mater%betapl      =    params(9)

      if (idebug.ge.2) then
         write(bufout,'(2a,i1,2(a,f6.1),a,/, 3(a,39x, 3(a,g10.3),:,a,/))')                              &
            trim(pfx_str(subnam,ire,icp)), ' H=', my_ic%heat, ', BkTemp(1)=', my_mater%bktemp(1),       &
                ', BkTemp(2)=', my_mater%bktemp(2), ',',                                                &
            pfx, 'HeatCp(1)=', my_mater%heatcp(1), ', Lambda(1)=', my_mater%lambda(1), ', Dens(1)=',    &
                my_mater%dens(1), ',',                                                                  &
            pfx, 'HeatCp(2)=', my_mater%heatcp(2), ', Lambda(2)=', my_mater%lambda(2), ', Dens(2)=',    &
                my_mater%dens(2), ',',                                                                  &
            pfx, 'Betapl=', my_mater%betapl
         call write_log(3, bufout)
      endif
   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setTemperatureData

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setTimestep(ire, icp, dt) &
   bind(c,name=CNAME_(cntc_settimestep))
!--function: set the time step size dt for a contact problem, particularly for T=0 or 1
!--category: 5, "m=any, wtd":       available for modules 1 and 3, in module 1 working on wtd data
!                                   dt is currently ignored in module 1.
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: icp            ! contact problem ID
   real(kind=8), intent(in) :: dt             ! time step size [time]
!--local variables:
   integer    :: ierror
   character(len=*), parameter :: subnam = 'cntc_setTimestep'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setTimestep
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, -1, subnam, ierror)
   if (ierror.lt.0) return

   my_kin%dt  = dt

   if (idebug.ge.2) then
      write(bufout,'(2a,f9.6,a)') trim(pfx_str(subnam,ire,icp)),' dt=',my_kin%dt,' [s]'
      call write_log(1, bufout)
   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setTimestep

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setReferenceVelocity(ire, icp, veloc) &
   bind(c,name=CNAME_(cntc_setreferencevelocity))
!--function: set the rolling velocity for a contact problem
!--category: 1, "m=3 only, cp":     available for module 3 only, working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: icp            ! contact problem ID
   real(kind=8), intent(in) :: veloc          ! absolute rolling velocity [veloc]
!--local variables:
   integer    :: ierror
   character(len=*), parameter :: subnam = 'cntc_setReferenceVelocity'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setReferenceVelocity
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 3, 1, subnam, ierror)
   if (ierror.lt.0) return

   my_kin%veloc = max(1d-9, abs(veloc) * my_scl%veloc)

   if (idebug.ge.2) then
      write(bufout,'(2a,f9.1,a)') trim(pfx_str(subnam,ire,icp)),' V=',my_kin%veloc,' [mm/s]'
      call write_log(1, bufout)
   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setReferenceVelocity

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setRollingStepsize(ire, icp, chi, dq) &
   bind(c,name=CNAME_(cntc_setrollingstepsize))
!--function: set the rolling direction chi and step size dq for a contact problem
!  in w/r contact,    icp = -1
!      chi            - ignored
!      dq             - rolling step size relative to grid size dx [-]
!  in generic contact, icp > 0,
!      chi            - rolling direction [angle]
!      dq             - rolling step size [length]
!--category: 5, "m=any, wtd":       available for modules 1 and 3, in module 1 working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: icp            ! contact problem ID
   real(kind=8), intent(in) :: chi, dq        ! rolling direction [angle] and step size [length/-]
!--local variables:
   integer    :: ierror
   character(len=*), parameter :: subnam = 'cntc_setRollingStepsize'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setRollingStepsize
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, -1, subnam, ierror)
   if (ierror.lt.0) return

   if (icp.le.0) then

      ! module 1 (icp=-1): set the relative rolling step size [no unit conversion]

      wtd%discr%dqrel = dq
      wtd%ic%discns1_inp = wtd%ic%discns1_eff   ! change D=1 to original value for write inp-file

   else

      ! module 3 (icp>0): set the absolute rolling step size

      my_kin%chi = chi * my_scl%angle
      my_kin%dq  = dq  * my_scl%len

      if (idebug.ge.2) then
         write(bufout,'(2a,f6.1,a,f6.3,a)') trim(pfx_str(subnam,ire,icp)),' chi=',my_kin%chi,           &
                ' [rad], dq=',my_kin%dq,' [mm]'
         call write_log(1, bufout)
      endif
   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setRollingStepsize

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setFrictionMethod(ire, icp, imeth, nparam, params) &
   bind(c,name=CNAME_(cntc_setfrictionmethod))
!--function: set V and L-digits and parameters for the friction law for a contact problem
!  imeth          - type of friction law used: 10 * V-digit + 1 * L-digit
!  L = 0: Coulomb friction,             lparam = [fstat, fkin]
!      1: keep previous parameters,     lparam = [ ]
!      2: linear falling friction,      lparam = [fkin, flin1, sabsh1, flin2, sabsh2, memdst, mem_s0]
!      3: rational falling friction,    lparam = [fkin, frat1, sabsh1, frat2, sabsh2, memdst, mem_s0]
!      4: exponential falling friction, lparam = [fkin, fexp1, sabsh1, fexp2, sabsh2, memdst, mem_s0]
!      5: exponential falling friction, lparam = [fstat, polach_a, polach_b]
!      6: temperature dep. friction,    lparam = [fref, tref, dfheat, dtheat, memdst, mem_s0]
!
!  When V = 0, params = lparam(1:n), with n=2 for L=0, n=7 for L=2--4, n=3 for L=5, n=6 for L=6
!  When V = 1, params = [ nvf, ...
!                         alphvf(1),   lparam( 1 ,1:n) ],  ...
!                                   ...
!                         alphvf(nvf), lparam(nvf,1:n) ]
!  When V = 2, same as for V = 1, using svf instead of alphvf
!    
!  dimensions:  alphvf: [angle],  fstat, fkin, flin1,2, frat1,2, fexp1,2, polach_a, fref, dfheat: [-]
!               sabsh1,2, mem_s0: [veloc],  svf, memdst: [length],  polach_b [1/veloc],  tref, dtheat [C]
!--category: 5, "m=any, wtd":       available for modules 1 and 3, in module 1 working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: icp            ! contact problem ID
   integer,      intent(in) :: imeth          ! VL = type of friction law used
   integer,      intent(in) :: nparam         ! number of parameters provided
   real(kind=8), intent(in) :: params(nparam) ! depending on method that is used
!--local variables:
   integer, parameter  :: nparam_loc(0:6) = (/ 2, 0, 7, 7, 7, 3, 6 /)
   logical      :: use_nvf
   integer      :: ierror, ldigit, vdigit, nexpct, imodul, nvf, ivf, iofs
   real(kind=8) :: fstat, fkin, polach_a, polach_b, fexp1, shalf
   character(len=*), parameter :: subnam = 'cntc_setFrictionMethod'
   type(t_friclaw),  pointer   :: my_fric
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setFrictionMethod
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, -1, subnam, ierror)
   if (ierror.lt.0) return

   imodul = ire_module(ix_reid(ire))
   ! if (imodul.eq.1) wtd%nvf = 0   ! used to signal errors

   ! Get / check V and L-digits

   vdigit = imeth / 10
   ldigit = imeth - 10*vdigit

   if (vdigit.lt.0 .or. vdigit.gt.2 .or. ldigit.lt.0 .or. ldigit.gt.6) then
      write(bufout,'(2a,i4,a)') trim(pfx_str(subnam,ire,icp)),' method',imeth,' does not exist.'
      call write_log(1, bufout)
      return
   endif

   if (ldigit.ne.1) then
      my_ic%varfrc = vdigit
   else
      vdigit = 0
   endif
   if (ldigit.eq.5) then
      my_ic%frclaw_inp = 4
   else
      my_ic%frclaw_inp = ldigit
   endif

   use_nvf = (vdigit.eq.1 .or. vdigit.eq.2)

   if (use_nvf .and. imodul.ne.1) then
      write(bufout,'(a,2(a,i1))') trim(pfx_str(subnam,ire,icp)),' V=',vdigit,' is not allowed in module', &
                imodul
      call write_log(1, bufout)
      return
   endif

   ! Get nvf, check number of parameters

   if (use_nvf .and. nparam.ge.1) then
      nvf    = nint(params(1))
      nexpct = 1 + nvf * (1 + nparam_loc(ldigit))
   else
      nvf    = 1
      nexpct = nparam_loc(ldigit)
   endif

   if (nparam.ne.nexpct) then
      if (use_nvf) then
         write(bufout,'(2a,i2,a,2(i3,a))') trim(pfx_str(subnam,ire,icp)),' with L=',ldigit,', NVF=',    &
                nvf,' needs', nexpct,' parameters.'
      else
         write(bufout,'(2a,2(i2,a))') trim(pfx_str(subnam,ire,icp)),' method L=',imeth,' needs',        &
                nexpct,' parameters.'
      endif
      call write_log(1, bufout)
      return
   endif

   ! do iofs = 1, nparam
   !    write(bufout, *) 'params(',iofs,')=',params(iofs)
   !    call write_log(1, bufout)
   ! enddo

   if (ldigit.ne.1) then

      ! Unpack / copy the input parameters

      if (imodul.eq.1) then
         my_fric => wtd%fric
      else
         my_fric => gd%fric
      endif
      my_fric%varfrc_eff = my_ic%varfrc

      ! L = 0 -- 6: store lparam(i,:)

      if (ldigit.eq.5) then
         my_fric%frclaw_eff = 4
      else
         my_fric%frclaw_eff = ldigit
      endif

      ! resize arrays to NVF entries

      call fric_resize(my_fric, nvf)

      ! copy parameters

      iofs = 0
      if (use_nvf) iofs = 1     ! skip 1st position: nvf

      do ivf = 1, nvf

         ! V = 1 or 2: store alphvf(i) or svf(i)

         if (use_nvf) then
            if (vdigit.eq.1) then
               my_fric%paramvf(ivf) = params(iofs+1) * my_scl%angle
            else 
               my_fric%paramvf(ivf) = params(iofs+1) * my_scl%len
            endif
            iofs = iofs + 1
         endif

         if (ldigit.eq.0) then

            ! L = 0: Coulomb friction with a single coefficient of friction fstat = fkin

            my_fric%fkin_arr(ivf)   = max(1d-6, min(params(iofs+1), params(iofs+2)))
            my_fric%fstat_arr(ivf)  = my_fric%fkin_arr(ivf)

            if (idebug.ge.2) then
               write(bufout,'(2a,i2,2(a,f6.3),a)') trim(pfx_str(subnam,ire,icp)),' L=',ldigit,          &
                   ', Fstat=',my_fric%fstat_arr(ivf),', Fkin=',my_fric%fkin_arr(ivf), ' [-]'
               call write_log(1, bufout)
            endif
    
         elseif (ldigit.eq.2) then
    
            ! Linear falling friction formula according to L=2
            ! TODO: consistency/lower bounds should be moved to m_friclaw
    
            my_fric%fkin_arr(ivf)  = max(1d-6, params(iofs+1))
            my_fric%flin1(ivf)     = max(0d0,  params(iofs+2))
            my_fric%sabsh1(ivf)    = max(1d-4, params(iofs+3) * my_scl%veloc)
            my_fric%flin2(ivf)     = max(0d0,  params(iofs+4))
            my_fric%sabsh2(ivf)    = max(1d-4, params(iofs+5) * my_scl%veloc)
            my_fric%memdst         = max(0d0,  params(iofs+6) * my_scl%len)
            my_fric%mem_s0         = max(0d0,  params(iofs+7) * my_scl%veloc)
            my_fric%fstat_arr(ivf) = my_fric%fkin_arr(ivf) + my_fric%flin1(ivf) + my_fric%flin2(ivf)
    
            if (idebug.ge.2) then
               write(bufout,'(2a,i2,a,f6.3,a,f7.1,a,f6.3,a,/, 43x,2(a,f5.2,a,f7.1),a)')                    &
                   trim(pfx_str(subnam,ire,icp)), ' L=',ldigit, ', Fkin=',my_fric%fkin_arr(ivf),           &
                                ' [-], Mem_s0=',my_fric%mem_s0,' [mm/s], Memdst=',my_fric%memdst,' [mm],', &
                   'Flin1=',my_fric%flin1(ivf), ' [-], Sabsh1=',my_fric%sabsh1(ivf), ' [mm/s], Flin2=',    &
                            my_fric%flin2(ivf), ' [-], Sabsh2=',my_fric%sabsh2(ivf), ' [mm/s]'
               call write_log(2, bufout)
            endif
    
         elseif (ldigit.eq.3) then
    
            ! Rational falling friction formula according to L=3
    
            my_fric%fkin_arr(ivf)  = max(1d-6, params(iofs+1))
            my_fric%frat1(ivf)     = max(0d0,  params(iofs+2))
            my_fric%sabsh1(ivf)    = max(1d-4, params(iofs+3) * my_scl%veloc)
            my_fric%frat2(ivf)     = max(0d0,  params(iofs+4))
            my_fric%sabsh2(ivf)    = max(1d-4, params(iofs+5) * my_scl%veloc)
            my_fric%memdst         = max(0d0,  params(iofs+6) * my_scl%len)
            my_fric%mem_s0         = max(0d0,  params(iofs+7) * my_scl%veloc)
            my_fric%fstat_arr(ivf) = my_fric%fkin_arr(ivf) + my_fric%frat1(ivf) + my_fric%frat2(ivf)
    
            if (idebug.ge.2) then
               write(bufout,'(2a,i2,a,f6.3,a,f7.1,a,f6.3,a,/, 43x,2(a,f5.2,a,f7.1),a)')                    &
                   trim(pfx_str(subnam,ire,icp)), ' L=',ldigit, ', Fkin=',my_fric%fkin_arr(ivf),           &
                                ' [-], Mem_s0=',my_fric%mem_s0,' [mm/s], Memdst=',my_fric%memdst,' [mm],', &
                   'Frat1=',my_fric%frat1(ivf), ' [-], Sabsh1=',my_fric%sabsh1(ivf), ' [mm/s], Frat2=',    &
                            my_fric%frat2(ivf), ' [-], Sabsh2=',my_fric%sabsh2(ivf), ' [mm/s]'
               call write_log(2, bufout)
            endif
    
         elseif (ldigit.eq.4) then
    
            ! Exponential falling friction formula according to L=4
    
            my_fric%fkin_arr(ivf)  = max(1d-6, params(iofs+1))
            my_fric%fexp1(ivf)     = max(0d0,  params(iofs+2))
            my_fric%sabsh1(ivf)    = max(1d-4, params(iofs+3) * my_scl%veloc)
            my_fric%fexp2(ivf)     = max(0d0,  params(iofs+4))
            my_fric%sabsh2(ivf)    = max(1d-4, params(iofs+5) * my_scl%veloc)
            my_fric%memdst         = max(0d0,  params(iofs+6) * my_scl%len)
            my_fric%mem_s0         = max(0d0,  params(iofs+7) * my_scl%veloc)
            my_fric%fstat_arr(ivf) = my_fric%fkin_arr(ivf) + my_fric%fexp1(ivf) + my_fric%fexp2(ivf)    
    
            if (idebug.ge.2) then
               write(bufout,'(2a,i2,a,f6.3,a,f7.1,a,f6.3,a,/, 2(48x,a,f5.2,a,f7.1,a,:,/))')                &
                   trim(pfx_str(subnam,ire,icp)), ' L=',ldigit, ', Fkin=',my_fric%fkin_arr(ivf),           &
                                ' [-], Mem_s0=',my_fric%mem_s0,' [mm/s], Memdst=',my_fric%memdst,' [mm],', &
                   'Fexp1=',my_fric%fexp1(ivf), ' [-], Sabsh1=',my_fric%sabsh1(ivf), ' [mm/s],',           &
                   'Fexp2=',my_fric%fexp2(ivf), ' [-], Sabsh2=',my_fric%sabsh2(ivf), ' [mm/s]'
               call write_log(3, bufout)
            endif
    
         elseif (ldigit.eq.5) then
    
            ! L = 5: Polach's exponential falling friction formula, converted to L=4
            !   mu = mu0 * [ A + (1-A) * exp(-B*s) ]
            !   fstat (s=0) == mu0;  fkin (s=inf) == mu0*A;  shalf: B*s == ln(2)
    
            fstat    = max(1d-6,params(iofs+1))
            polach_a = max(0d0, min(1d0, params(iofs+2)))
            polach_b = params(iofs+3) / my_scl%veloc
    
            fkin  = polach_a * fstat
            fexp1 = fstat - fkin
            shalf = log(2d0)/polach_b
            my_fric%fstat_arr(ivf) = fstat
            my_fric%fkin_arr(ivf)  = fkin
            my_fric%fexp1(ivf)     = fexp1
            my_fric%sabsh1(ivf)    = shalf
            my_fric%fexp2(ivf)     =     0d0
            my_fric%sabsh2(ivf)    = 10000d0   ! 10 m/s
            my_fric%memdst         = 0.001d0   ! 1 micrometer
            my_fric%mem_s0         =    10d0   ! 10 mm/s
    
            if (idebug.ge.2) then
               write(bufout,'(2a,i2,2(a,f5.2),a,/, 42x,2(a,f7.1),a,f5.2,a)')                            &
                   trim(pfx_str(subnam,ire,icp)), ' L=',ldigit, ', Fkin=',my_fric%fkin_arr(ivf),        &
                           ', Fexp1=',my_fric%fexp1,' [-],', ' Sabsh1=',my_fric%sabsh1,                 &
                           ', Mem_s0=',my_fric%mem_s0, ' [mm/s], Memdst=',my_fric%memdst,' [mm]'
               call write_log(2, bufout)
            endif
    
         elseif (ldigit.eq.6) then
    
            ! L = 6: Temperature dependent friction with piecewise linear formula
    
            my_fric%fref(ivf)      = max(1d-6, params(iofs+1))
            my_fric%tref(ivf)      = params(iofs+2)
            my_fric%dfheat(ivf)    = max(1d-6-my_fric%fref(ivf), params(iofs+3))
            my_fric%dtheat(ivf)    = params(iofs+4)
            my_fric%memdst         = max(0d0,  params(iofs+5) * my_scl%len)
            my_fric%mem_s0         = max(0d0,  params(iofs+6) * my_scl%veloc)
            my_fric%fstat_arr(ivf) = my_fric%fref(ivf)
            my_fric%fkin_arr(ivf)  = my_fric%fref(ivf) + my_fric%dfheat(ivf)
    
            if (idebug.ge.2) then
               write(bufout,'(2a,i2,a,f6.3,a,f7.1,a,/, a,40x,a,f6.3,a,f7.1,a)')                            &
                   trim(pfx_str(subnam,ire,icp)), ' L=',my_fric%frclaw_eff, ', Fref=',my_fric%fref(ivf),   &
                           ' [-], Tref=',my_fric%tref(ivf),' [C],',                                        &
                   pfx,' Dfheat =',my_fric%dfheat(ivf),' [-], Dtheat =',my_fric%dtheat(ivf),' [C]'
               call write_log(2, bufout)
            endif
    
         endif

         iofs = iofs + nparam_loc(ldigit)

      enddo ! ivf

      call fric_update(my_fric)

      ! set the value of the friction coefficient used for scaling tang. forces

      my_kin%use_muscal = (my_ic%varfrc.eq.0 .or. my_ic%varfrc.eq.2)
      if (my_kin%use_muscal) then
         my_kin%muscal = my_fric%fstat()
      else
         my_kin%muscal = 1d0
      endif

   endif ! L <> 1

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setFrictionMethod

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setHertzContact(ire, icp, ipotcn, nparam, params) &
   bind(c,name=CNAME_(cntc_sethertzcontact))
!--function: set the parameters for a Hertzian contact problem
!   -6: SDEC approach, union of two half ellipses                  params = [ mx, my, aa , bneg, bpos, scale ]
!   -5: Hertzian rectangular contact, half-sizes prescribed,       params = [ mx, my, aa , bb , scale ]
!   -4: Hertzian rectangular contact, curv+half width prescribed,  params = [ mx, my, a1 , bb , scale ]
!   -3: Hertzian elliptical contact, semi-axes prescribed,         params = [ mx, my, aa , bb , scale ]
!   -2: Hertzian elliptical contact, ellipticity prescribed,       params = [ mx, my, a1 , aob, scale ]
!   -1: Hertzian elliptical contact, curvatures prescribed,        params = [ mx, my, a1 , b1 , scale ]
!
! dimensions:  mx, my, aob, scale [-],    a1, b1: [1/length],    aa, bneg, bpos, bb: [length]
!--category: 1, "m=3 only, cp":     available for module 3 only, working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: icp            ! contact problem ID
   integer,      intent(in) :: ipotcn         ! type of potcon specification
   integer,      intent(in) :: nparam         ! number of parameters provided
   real(kind=8), intent(in) :: params(nparam) ! dependent on type of specification [length]
!--local variables:
   integer, parameter  :: nparam_loc(-6:-1) = (/ 6, 5, 5, 5, 5, 5 /)
   integer      :: ierror
   character(len=*), parameter :: subnam = 'cntc_setHertzContact'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setHertzContact
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 3, 1, subnam, ierror)
   if (ierror.lt.0) return

   if (ipotcn.lt.-6 .or. ipotcn.gt.-1) then
      write(bufout,'(2a,i2,a)') trim(pfx_str(subnam,ire,icp)),' method IPOTCN=',ipotcn,' does not exist.'
      call write_log(1, bufout)
      return
   endif
   if (nparam.ne.nparam_loc(ipotcn)) then
      write(bufout,'(2a,2(i2,a))') trim(pfx_str(subnam,ire,icp)),' method IPOTCN=',ipotcn,' needs',     &
                nparam_loc(ipotcn),' parameters.'
      call write_log(1, bufout)
      return
   endif

   ! new inputs are stored in potcon_inp, will be copied to potcon_cur in cntc_calculate

   gd%potcon_inp%ipotcn = ipotcn

   if (ipotcn.eq.-6) then

      !   -6: SDEC approach, union of two half ellipses            params = [ mx, my, aa , bneg, bpos, scale ]

      gd%potcon_inp%mx   = max(1, nint(params(1)))
      gd%potcon_inp%my   = max(1, nint(params(2)))
      gd%hertz%aa    = max(1d-6,params(3)) * my_scl%len
      gd%hertz%bneg  = max(1d-6,params(4)) * my_scl%len
      gd%hertz%bpos  = max(1d-6,params(5)) * my_scl%len
      gd%hertz%scale = max(1d-6,params(6))

   elseif (ipotcn.eq.-5) then

      !   -5: Hertzian rectangular contact, half-sizes prescribed,       params = [ mx, my, aa , bb , scale ]

      gd%ic%bound    = 2        ! semi-elliptical pressure in x
      gd%potcon_inp%mx   = max(1, nint(params(1)))
      gd%potcon_inp%my   = max(1, nint(params(2)))
      gd%hertz%aa    = max(1d-6,params(3)) * my_scl%len
      gd%hertz%bb    = max(1d-6,params(4)) * my_scl%len
      gd%hertz%scale = max(1d-6,params(5))

   elseif (ipotcn.eq.-4) then

      !   -4: Hertzian rectangular contact, curv+half width prescribed,  params = [ mx, my, a1 , bb , scale ]

      gd%ic%bound    = 2        ! semi-elliptical pressure in x
      gd%potcon_inp%mx = max(1, nint(params(1)))
      gd%potcon_inp%my = max(1, nint(params(2)))
      gd%hertz%a1      = max(1d-6,params(3)) / my_scl%len
      gd%hertz%bb      = max(1d-6,params(4)) * my_scl%len
      gd%hertz%scale   = max(1d-6,params(5))

   elseif (ipotcn.eq.-3) then

      !   -3: Hertzian elliptical contact, semi-axes prescribed,         params = [ mx, my, aa , bb , scale ]

      gd%potcon_inp%mx = max(1, nint(params(1)))
      gd%potcon_inp%my = max(1, nint(params(2)))
      gd%hertz%aa      = max(1d-6,params(3)) * my_scl%len
      gd%hertz%bb      = max(1d-6,params(4)) * my_scl%len
      gd%hertz%scale   = max(1d-6,params(5))

   elseif (ipotcn.eq.-2) then

      !   -2: Hertzian elliptical contact, ellipticity prescribed,       params = [ mx, my, a1 , aob, scale ]

      gd%potcon_inp%mx = max(1, nint(params(1)))
      gd%potcon_inp%my = max(1, nint(params(2)))
      gd%hertz%a1      = max(1d-12,params(3)) / my_scl%len
      gd%hertz%aob     = max(1d-6 ,params(4))
      gd%hertz%scale   = max(1d-6 ,params(5))

   elseif (ipotcn.eq.-1) then

      !   -1: Hertzian elliptical contact, curvatures prescribed,        params = [ mx, my, a1 , b1 , scale ]

      gd%potcon_inp%mx = max(1, nint(params(1)))
      gd%potcon_inp%my = max(1, nint(params(2)))
      gd%hertz%a1      = max(1d-12,params(3)) / my_scl%len
      gd%hertz%b1      = max(1d-12,params(4)) / my_scl%len
      gd%hertz%scale   = max(1d-6 ,params(5))

   endif ! ipotcn

   gd%potcon_inp%npot  = gd%potcon_inp%mx * gd%potcon_inp%my

   if (idebug.ge.2) then
      write(bufout,'(2a,i3,a,2i5)') trim(pfx_str(subnam,ire,icp)),' ipotcn=',gd%potcon_inp%ipotcn,      &
                ', mx,my=',gd%potcon_inp%mx,gd%potcon_inp%my
      call write_log(1, bufout)
   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setHertzContact

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setPotContact(ire, icp, ipotcn, nparam, params) &
   bind(c,name=CNAME_(cntc_setpotcontact))
!--function: set the parameters of the potential contact area for a contact problem
!  for w/r contact, icp = -1:
!    0: w/r contact, fixed grid sizes,  params = [ dx, ds, n.a. ]
!   -1: w/r contact, fixed grid sizes,  params = [ dx, ds, a_sep, d_sep, d_comb, [d_turn, g_miss] ]
!
!  for generic contact, icp > 0:
!    1: lower-left + grid sizes,        params = [ mx, my, xl , yl , dx , dy  ]
!    2: lower-left + upper right,       params = [ mx, my, xl , yl , xh , yh  ]
!    3: 1st center + grid sizes,        params = [ mx, my, xc1, yc1, dx , dy  ]
!    4: 1st center + last center,       params = [ mx, my, xc1, yc1, xcm, ycm ]
!
!  dimensions: mx, my [-],   a_sep [angle],
!              dx, ds, dy, d_sep, d_comb, d_turn, g_miss, xl, yl, xh, yh, xc1, yc1, xcm, ycm [length]
!--category: 5, "m=any, wtd":       available for modules 1 and 3, in module 1 working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: icp            ! contact problem ID
   integer,      intent(in) :: ipotcn         ! type of potcon specification
   integer,      intent(in) :: nparam         ! number of parameters provided
   real(kind=8), intent(in) :: params(nparam) ! dependent on type of specification [length]
!--local variables:
   integer, parameter  :: nparam_loc(-1:4) = (/ 6, 3, 6, 6, 6, 6 /)
   integer      :: ierror
   character(len=*), parameter :: subnam = 'cntc_setPotContact'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setPotContact
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, -1, subnam, ierror)
   if (ierror.lt.0) return

   if (icp.eq.-1) then

      ! w/r contact: always using grid sizes, params = [ dx, dy, n.a. ]

      if (ipotcn.lt.-1 .or. ipotcn.gt.0) then
         write(bufout,'(2a,i4,a)') trim(pfx_str(subnam,ire,icp)),' method',ipotcn,' does not exist.'
         call write_log(1, bufout)
         return
      endif

      if (ipotcn.eq.-1) then
         if (nparam.lt.5 .or. nparam.gt.7) then
            write(bufout,'(2a,i2,a)') trim(pfx_str(subnam,ire,icp)),' method IPOTCN=',ipotcn,           &
                ' needs 5, 6 or 7 parameters.'
            call write_log(1, bufout)
            return
         endif
      else
         if (nparam.ne.nparam_loc(ipotcn)) then
            write(bufout,'(2a,2(i2,a))') trim(pfx_str(subnam,ire,icp)),' method IPOTCN=',ipotcn,        &
                ' needs', nparam_loc(ipotcn),' parameters.'
            call write_log(1, bufout)
            return
         endif
      endif

      wtd%discr%dx     = max(1d-12, params(1)) * my_scl%len
      wtd%discr%ds     = max(1d-12, params(2)) * my_scl%len

      if (ipotcn.eq.-1) then
         wtd%discr%angl_sep  = max(0d0, params(3)) * my_scl%angle
         wtd%discr%dist_sep  = max(0d0, params(4)) * my_scl%len
         wtd%discr%dist_comb = max(0d0, params(5)) * my_scl%len

         if (nparam.lt.6) then
            wtd%discr%dist_turn = 2d0 * wtd%discr%dist_sep - wtd%discr%dist_comb ! default D_TURN
         else
            wtd%discr%dist_turn = max(0d0, params(6)) * my_scl%len
         endif
         if (nparam.lt.7) then
            wtd%discr%gap_miss  = -1d0                          ! default G_MISS
         else
            wtd%discr%gap_miss  = params(7) * my_scl%len
         endif
      endif

   else

      ! module 3: various non-Hertzian options

      if (ipotcn.lt.1 .or. ipotcn.gt.4) then
         write(bufout,'(2a,i2,a)') trim(pfx_str(subnam,ire,icp)),' method',ipotcn,' does not exist.'
         call write_log(1, bufout)
         return
      endif
      if (nparam.ne.nparam_loc(ipotcn)) then
         write(bufout,'(2a,2(i2,a))') trim(pfx_str(subnam,ire,icp)),' method IPOTCN=',ipotcn,' needs',  &
                nparam_loc(ipotcn),' parameters.'
         call write_log(1, bufout)
         return
      endif

      ! new inputs are stored in potcon_inp, will be copied to potcon_cur in cntc_calculate

      gd%potcon_inp%ipotcn = ipotcn

      if (ipotcn.eq.1) then

         !    1: lower-left + grid sizes,       params = [ mx, my, xl , yl , dx , dy ]

         gd%potcon_inp%mx   = max(1, nint(params(1)))
         gd%potcon_inp%my   = max(1, nint(params(2)))
         gd%potcon_inp%xl   =            params(3)  * my_scl%len
         gd%potcon_inp%yl   =            params(4)  * my_scl%len
         gd%potcon_inp%dx   = max(1d-12, params(5)) * my_scl%len
         gd%potcon_inp%dy   = max(1d-12, params(6)) * my_scl%len

         if (idebug.ge.2) then
            write(bufout,'(2a,2f8.3,a)') trim(pfx_str(subnam,ire,icp)),' coordinates xl,yl=',           &
                      gd%potcon_inp%xl, gd%potcon_inp%yl,' [mm]'
            call write_log(1, bufout)
         endif

      elseif (ipotcn.eq.2) then

         !    2: lower-left + upper right,      params = [ mx, my, xl , yl , xh , yh ]

         gd%potcon_inp%mx   = max(1, nint(params(1)))
         gd%potcon_inp%my   = max(1, nint(params(2)))
         gd%potcon_inp%xl   = params(3) * my_scl%len
         gd%potcon_inp%yl   = params(4) * my_scl%len
         gd%potcon_inp%xh   = params(5) * my_scl%len
         gd%potcon_inp%yh   = params(6) * my_scl%len

      elseif (ipotcn.eq.3) then

         !    3: 1st center + grid sizes,       params = [ mx, my, xc1, yc1, dx , dy ]

         gd%potcon_inp%mx   = max(1, nint(params(1)))
         gd%potcon_inp%my   = max(1, nint(params(2)))
         gd%potcon_inp%xc1  =            params(3)  * my_scl%len
         gd%potcon_inp%yc1  =            params(4)  * my_scl%len
         gd%potcon_inp%dx   = max(1d-12, params(5)) * my_scl%len
         gd%potcon_inp%dy   = max(1d-12, params(6)) * my_scl%len

      elseif (ipotcn.eq.4) then

         !    4: 1st center + last center,      params = [ mx, my, xc1, yc1, xcm, ycm]

         gd%potcon_inp%mx   = max(1, nint(params(1)))
         gd%potcon_inp%my   = max(1, nint(params(2)))
         gd%potcon_inp%xc1  = params(3) * my_scl%len
         gd%potcon_inp%yc1  = params(4) * my_scl%len
         gd%potcon_inp%xcm  = params(5) * my_scl%len
         gd%potcon_inp%ycm  = params(6) * my_scl%len

      endif ! ipotcn

      ! complete remaining entries in potcon_inp including npot

      call potcon_fill( gd%potcon_inp )

   endif ! module 1

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setPotContact

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setVerticalForce(ire, fz) &
   bind(c,name=CNAME_(cntc_setverticalforce))
!--function: set the total vertical force between the contacting bodies for a w/r contact problem (module 1)
!            Note: this function sets control-digit N = 1
!--category: 2, "m=1 only, wtd":    available for module 1 only, working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   real(kind=8), intent(in) :: fz             ! total vertical force between the two bodies [force]
!--local variables:
   integer    :: ierror
   character(len=*), parameter :: subnam = 'cntc_setVerticalForce'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setVerticalForce
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, -1, 1, -1, subnam, ierror)
   if (ierror.lt.0) return

   my_ic%norm    = 1
   wtd%ws%fz_inp = fz * my_scl%forc

   if (idebug.ge.2) then
      write(bufout,'(2a,f9.2,a)') trim(pfx_str(subnam,ire,-1)),' vertical force=', wtd%ws%fz_inp, ' [N]'
      call write_log(1, bufout)
   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setVerticalForce

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setPenetration(ire, icp, pen) &
   bind(c,name=CNAME_(cntc_setpenetration))
!--function: set the approach (penetration) of the bodies as a whole for a contact problem
!            Note: this function sets control-digit N = 0
!--category: 1, "m=3 only, cp":     available for module 3 only, working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: icp            ! contact problem ID
   real(kind=8), intent(in) :: pen            ! penetration/approach of the two bodies [length]
!--local variables:
   integer    :: ierror
   character(len=*), parameter :: subnam = 'cntc_setPenetration'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setPenetration
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 3, 1, subnam, ierror)
   if (ierror.lt.0) return

   my_ic%norm = 0
   my_kin%pen = pen * my_scl%len

   if (idebug.ge.2) then
      write(bufout,'(2a,f7.4,a)') trim(pfx_str(subnam,ire,icp)),' penetration=',my_kin%pen, ' [mm]'
      call write_log(1, bufout)
   endif
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setPenetration

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setNormalForce(ire, icp, fn) &
   bind(c,name=CNAME_(cntc_setnormalforce))
!--function: set the total normal force of the bodies as a whole for a contact problem
!            Note: this function sets control-digit N = 1
!--category: 1, "m=3 only, cp":     available for module 3 only, working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: icp            ! contact problem ID
   real(kind=8), intent(in) :: fn             ! total normal force between the two bodies [force]
!--local variables:
   integer    :: ierror
   character(len=*), parameter :: subnam = 'cntc_setNormalForce'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setNormalForce
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 3, 1, subnam, ierror)
   if (ierror.lt.0) return

   my_ic%norm    = 1
   my_kin%fntrue = fn

   if (idebug.ge.2) then
      write(bufout,'(2a,f9.2,a)') trim(pfx_str(subnam,ire,icp)),' normal force=',my_kin%fntrue, ' [N]'
      call write_log(1, bufout)
   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setNormalForce

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setUndeformedDistc(ire, icp, ibase, nparam, prmudf) &
   bind(c,name=CNAME_(cntc_setundeformeddistc))
!--function: set the undeformed distance function through a formula or by element-wise specification
!    1: quadratic function            params = [b1, b2, b3, b4, b5, b6]
!    2: circular-x, piecewise-lin-y   params = [nn, xm, rm, y1, dy1], [b(k), k=1..nn]
!    3: quadratic plus two sines      params = [b1, b2, b3, b4, b5, b6, b7, b8]
!    9: element-wise specification    params = [h(i), i=1..npot] - undeformed distance per element [length]
!                                               note: positive values == separation between profiles
!--category: 1, "m=3 only, cp":     available for module 3 only, working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire            ! result element ID
   integer,      intent(in)  :: icp            ! contact problem ID
   integer,      intent(in)  :: ibase          ! type of undeformed distance specification
   integer,      intent(in)  :: nparam         ! number of parameters provided
   real(kind=8), intent(in)  :: prmudf(nparam) ! parameters of undef.dist, depending on method that is used
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_setUndeformedDistc'
   integer :: nn, nparam_loc
   integer :: ii, iimax, iimin, ierror
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setUndeformedDistc
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 3, 1, subnam, ierror)
   if (ierror.lt.0) return

   if (idebug.ge.2) then
      write(bufout,'(2a,i2,a,i6,a)') trim(pfx_str(subnam,ire,icp)),' set undeformed distance ibase=', &
                ibase,',', nparam,' parameters'
      call write_log(1, bufout)
   endif

   ! check value of ibase and check corresponding number of parameters

   if (ibase.eq.2) nn = nint(prmudf(1))

   if (ibase.eq.1) then
      nparam_loc = 6
   elseif (ibase.eq.2) then
      nparam_loc = 5 + nn
   elseif (ibase.eq.3) then
      nparam_loc = 8
   elseif (ibase.eq.9) then
      nparam_loc = gd%potcon_inp%npot
   else
      write(bufout,'(2a,i4,a)') trim(pfx_str(subnam,ire,icp)),' method',ibase,' does not exist.'
      call write_log(1, bufout)
      return
   endif

   if (nparam.ne.nparam_loc) then
      write(bufout,'(2a,i1,a,i6,a)') trim(pfx_str(subnam,ire,icp)),' method IBASE=',ibase,' needs ',    &
                nparam_loc,' parameters.'
      call write_log(1, bufout)
      return
   endif

   ! set type of specification of undeformed distance

   gd%geom%ibase = ibase

   ! copy data for the undeformed distance function

   if (ibase.eq.1) then

      ! h = b1 * x^2 + ... + b6
      ! units are [1/mm] for b1,b2,b3, [-] for b4,b5, [mm] for b6

      gd%geom%prmudf(1:3) = prmudf(1:3) / my_scl%len
      gd%geom%prmudf(4:5) = prmudf(4:5)
      gd%geom%prmudf( 6 ) = prmudf( 6 ) * my_scl%len

   elseif (ibase.eq.2) then

      ! re-allocate geom%prmudf when needed

      call reallocate_arr(gd%geom%prmudf, 5+nn)

      ! parameters are [nn, xm, rm, y1, dy1] and nn values [hi]
      ! units are [mm] for all parameters except nn

      gd%geom%nn = nn
      gd%geom%prmudf(2:5+nn) = prmudf(2:5+nn) * my_scl%len

   elseif (ibase.eq.3) then

      ! h = b1*sin(b2*(x-b3)) - b4*sin(b5*(x-b6)) - PEN + x**2/b7 + y**2/b8
      ! b1,b4 are amplitues [mm]
      ! b2,b5 are frequencies in [rad/mm]
      ! b3,b6 are phase shifts in [mm]
      ! b7,b8 are radii of curvature [mm]

      gd%geom%prmudf( 1 ) = prmudf( 1 ) * my_scl%len
      gd%geom%prmudf( 2 ) = prmudf( 2 ) / my_scl%len
      gd%geom%prmudf(3:4) = prmudf(3:4) * my_scl%len
      gd%geom%prmudf( 5 ) = prmudf( 5 ) / my_scl%len
      gd%geom%prmudf(6:8) = prmudf(6:8) * my_scl%len

   elseif (ibase.eq.9) then

      ! re-allocate prmudf at the appropriate size

      call reallocate_arr(gd%geom%prmudf, gd%potcon_inp%npot+10)

      ! copy input-data, determine minimum/maximum undef.distc as well

      iimin = 1
      iimax = 1
      do ii = 1, gd%potcon_inp%npot
         if (ii.le.nparam) then
            if (prmudf(ii).lt.prmudf(iimin)) iimin = ii
            if (prmudf(ii).gt.prmudf(iimax)) iimax = ii
            gd%geom%prmudf(ii) = prmudf(ii) * my_scl%len
         endif
      enddo
      if (idebug.ge.2) then
         write(bufout,'(2a,2(f7.4,a))') trim(pfx_str(subnam,ire,icp)),' max separation=',               &
                   min(99.9d0,gd%geom%prmudf(iimax)), ', max penetration=', gd%geom%prmudf(iimin),' [mm]'
         call write_log(1, bufout)
      endif

   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setUndeformedDistc

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setCreepages(ire, icp, vx, vy, phi) &
   bind(c,name=CNAME_(cntc_setcreepages))
!--function: set the kinematic constants (creepages) for a contact problem
!            note: vx is ignored when F=1 or 2, vy is ignored when F=2.
!--category: 1, "m=3 only, cp":     available for module 3 only, working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: icp            ! contact problem ID
   real(kind=8), intent(in) :: vx, vy, phi    ! long/lat/spin creepages [-], [angle/length]
!--local variables:
   integer    :: ierror
   character(len=*), parameter :: subnam = 'cntc_setCreepages'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setCreepages
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 3, 1, subnam, ierror)
   if (ierror.lt.0) return

   if (my_ic%tang.eq.1) then
      ! setting the rigid shift [mm], [rad]
      my_kin%cksi = vx * my_scl%len
      my_kin%ceta = vy * my_scl%len
      my_kin%cphi = phi
   else
      ! setting the creepages [-], [rad/mm]
      my_kin%cksi = vx
      my_kin%ceta = vy
      my_kin%cphi = phi / my_scl%len
   endif

   if (idebug.ge.2) then
      write(bufout,'(2a,2f7.4,a,f7.4,a)') trim(pfx_str(subnam,ire,icp)),' creepages Cksi,Ceta=',        &
                my_kin%cksi, my_kin%ceta, ' [-], Phi=',my_kin%cphi, ' [rad/mm]'
      call write_log(1, bufout)
   endif
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setCreepages

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setExtraRigidSlip(ire, icp, lenarr, wx, wy) &
   bind(c,name=CNAME_(cntc_setextrarigidslip))
!--function: set the extra term of the tangential right hand side for all elements in the potential contact
!            area for a contact problem
!--category: 1, "m=3 only, cp":     available for module 3 only, working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   integer,      intent(in)  :: lenarr        ! length of input array
   real(kind=8), intent(in)  :: wx(lenarr), wy(lenarr) ! extra term of tangential rhs for all elements of contact area [-]
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_setExtraRigidSlip'
   integer      :: ii, iimax, ierror
   real(kind=8) :: wabs, wmax
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setExtraRigidSlip
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 3, 1, subnam, ierror)
   if (ierror.lt.0) return

   if (idebug.ge.2) then
      write(bufout,'(2a)') trim(pfx_str(subnam,ire,icp)),' set extra rigid slip array exrhs [-]'
      call write_log(1, bufout)
   endif

   ! set option for adding extra term to right hand side per element

   my_ic%rztang = 1

   ! re-allocate exrhs at the appropriate size

   call grid_set_dimens(gd%cgrid_inp, gd%potcon_inp%mx, gd%potcon_inp%my)
   call gf3_new(gd%geom%exrhs, 'geom%exrhs', gd%cgrid_inp)
   call gf3_set(AllElm, 0d0, gd%geom%exrhs, ikTANG)

   ! copy input-data, determine minimum/maximum values as well

   iimax = 1
   wmax = 0d0
   do ii = 1, gd%potcon_inp%npot
      if (ii.le.lenarr) then
         wabs = sqrt(wx(ii)**2 + wy(ii)**2)
         if (wabs.gt.wmax) then
            iimax = ii
            wmax = wabs
         endif
         gd%geom%exrhs%vx(ii) = wx(ii)
         gd%geom%exrhs%vy(ii) = wy(ii)
      endif
   enddo
   if (idebug.ge.-2) then
      write(bufout,'(2a,f7.4,a,i5,a)') trim(pfx_str(subnam,ire,icp)),' max extra rigid slip=',          &
                min(9.9d0,wmax), ' at element ii=',iimax
      call write_log(1, bufout)
   endif
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setExtraRigidSlip

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setTangentialForces(ire, icp, fx, fy) &
   bind(c,name=CNAME_(cntc_settangentialforces))
!--function: set the total tangential forces for a contact problem
!            note: Fx is ignored when F=0, Fy is ignored when F=0 or 1.
!--external: category: 1, "m=3 only, cp":     available for module 3 only, working on cp data
!--internal: category: 5, "m=any, wtd":       available for modules 1 and 3, in module 1 working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: icp            ! contact problem ID
   real(kind=8), intent(in) :: fx, fy         ! total tangential forces relative to fstat*fn [-]
!--local variables:
   integer      :: ierror
   real(kind=8) :: fabs
   character(len=*), parameter :: subnam = 'cntc_setTangentialForces'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setTangentialForces
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, -1, subnam, ierror)
   if (ierror.lt.0) return

   if (icp.le.0) then

      ! module 1: absolute total force Fx on rail in wheelset x-direction, [N]

      wtd%ws%fx_inp = fx * my_scl%forc

   else

      ! module 3: scaled total forces Fx, Fy on body 1, relative to fstat*Fn [-]

      my_kin%fxrel = fx
      my_kin%fyrel = fy

      ! scale forces when necessary, such that |F| <= 1

      fabs = sqrt(fx**2 + fy**2)
      if (fabs.gt.1d0) then
         my_kin%fxrel = my_kin%fxrel / fabs
         my_kin%fyrel = my_kin%fyrel / fabs
      endif

      if (idebug.ge.2) then
         write(bufout,'(2a,2f8.4,a)') trim(pfx_str(subnam,ire,icp)),' total forces Fx,Fy=',             &
                my_kin%fxrel, my_kin%fyrel, ' [-]'
         call write_log(1, bufout)
      endif

   endif ! icp<0

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setTangentialForces

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setProfileInputFname(ire, c_fname, len_fname, nints, iparam, nreals, rparam) &
   bind(c,name=CNAME_(cntc_setprofileinputfname))
!--function: set a wheel or rail profile filename for a wheel-rail contact problem
!  fname          - string: name of profile file, absolute path or relative wrt effective working folder
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
!--category: 2, "m=1 only, wtd":    available for module 1 only, working on wtd data
   implicit none
!--subroutine arguments:
   integer,                intent(in) :: ire          ! result element ID
   character(kind=C_CHAR), intent(in) :: c_fname(*)   ! C-string: name of profile file
   integer,                intent(in) :: len_fname    ! length of filename
   integer,                intent(in) :: nints, nreals
   integer,                intent(in) :: iparam(nints)
   real(kind=8),           intent(in) :: rparam(nreals)
!--local variables:
   integer                     :: ierror, i, ilen, is_spck, itype, itimer
   type(t_profile),  pointer   :: my_prf
   character(len=256)          :: f_fname
   character(len=*), parameter :: subnam = 'cntc_setProfileInputFname'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setProfileInputFname
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, -1, 1, -1, subnam, ierror)
   if (ierror.lt.0) return

   ! initialize grid function

   ! Copy C-string (array of chars) to F-string

   f_fname = ' '
   ilen = min(len(f_fname), len_fname)
   do i = 1, ilen
      f_fname(i:i) = c_fname(i)
   enddo

   ! unpack option values, fill in defaults

   if (nints.ge.1) then
      itype   = iparam(1) ! 0=rail, 1=wheel, -1=unknown
   else
      itype   = -1        ! unknown type
   endif

   ! unknown type: taken from filename extension

   if (itype.ne.0 .and. itype.ne.1) then
      call profile_get_filetype(f_fname, itype, is_spck, idebug)
   endif

   ! set profile pointer for wheel or rail

   if (itype.eq.1) then
      my_prf => wtd%ws%whl%prw
   else
      my_prf => wtd%trk%rai%prr
   endif

   ! store configuration parameters from iparam/rparam

   call profile_setopt(my_prf, nints, iparam, nreals, rparam, my_scl%len)

   my_prf%fname      = f_fname

   ! read actual profile data

   itimer = cntc_timer_num(ire, -1)
   call timer_start(itimer)
   call profile_read_file(my_prf, my_meta%wrkdir, itype, wtd%ic%x_profil, wtd%ic%x_readln, .false.)
   call timer_stop(itimer)

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setProfileInputFname

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setProfileInputValues(ire, npoint, values, nints, iparam, nreals, rparam) &
   bind(c,name=CNAME_(cntc_setprofileinputvalues))
!--function: set a wheel or rail profile for a wheel-rail contact problem using a table of values
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
!--category: 2, "m=1 only, wtd":    available for module 1 only, working on wtd data
   implicit none
!--subroutine arguments:
   integer,                intent(in) :: ire              ! result element ID
   integer,                intent(in) :: npoint           ! number of profile points in values
   real(kind=8),           intent(in) :: values(2*npoint) ! profile values, ordered [y1,z1, y2,z2, ...]
   integer,                intent(in) :: nints, nreals
   integer,                intent(in) :: iparam(nints)
   real(kind=8),           intent(in) :: rparam(nreals)
!--local variables:
   integer                     :: ierror, itype, my_idebug
   type(t_profile),  pointer   :: my_prf
   character(len=*), parameter :: subnam = 'cntc_setProfileInputValues'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setProfileInputValues
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, -1, 1, -1, subnam, ierror)
   if (ierror.lt.0) return

   ! unpack option values, fill in defaults

   if (nints.ge.1) then
      itype   = iparam(1) ! 0=rail, 1=wheel, -1=unknown (n.a.)
      itype   = max(0, min(1, itype))
   else
      itype   = 1         ! assuming wheel!
   endif

   ! set profile pointer for wheel or rail

   if (itype.eq.1) then
      my_prf => wtd%ws%whl%prw
   else
      my_prf => wtd%trk%rai%prr
   endif

   ! store configuration parameters from iparam/rparam

   call profile_setopt(my_prf, nints, iparam, nreals, rparam, my_scl%len)

   my_prf%fname   = '<values>'

   ! store actual profile data

   my_idebug = max(idebug, wtd%ic%x_profil)
   call profile_store_values(my_prf, itype, npoint, values, my_idebug)

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setProfileInputValues

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setTrackDimensions(ire, ztrack, nparam, params) &
   bind(c,name=CNAME_(cntc_settrackdimensions))
!--function: set the track or roller-rig description for a wheel-rail contact problem
!  ztrack - control digit ZTRACK      params - depending on method that is used
!    0: maintain track dimensions     params = [ ]
!    1: new design track dimensions   params = [gaught,  dummy, gaugwd, cant, nomrad, curv],  if gaught >  0,
!                                         or   [gaught, raily0, railz0, cant, nomrad, curv],  if gaught <= 0.
!    2: new track deviations          params = [dyrail, dzrail, drollr, vyrail, vzrail, vrollr]
!    3: new dimensions & track deviations for current side of the track
!                                     params = params(1:6) cf. Z=1 followed by params(7:12) cf. Z=2;
!  ztrack >= 30 is used to configure the massless rail model, F=3
!   32/33: "F=3, Z=2/3":              params = [kyrail, fyrail, kzrail, fzrail, {dystep0} ]
!
! dimensions: gaught, gaugwd, raily0, railz0, nomrad  [length],  cant  [angle],  curv [1/length]
!             dyrail, dzrail  [length],  drollr  [angle],  vyrail, vzrail  [veloc], vrollr  [ang.veloc]
!             kyrail, kzrail  [force/length],  fyrail, fzrail  [force],  dystep0  [length]
!--category: 2, "m=1 only, wtd":    available for module 1 only, working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: ztrack         ! control digit ZTRACK
   integer,      intent(in) :: nparam         ! number of parameters provided
   real(kind=8), intent(in) :: params(nparam) ! parameters depending on method that is used
!--local variables:
   integer             :: nparam_loc(0:33, 2), iofs, ierror
   character(len=*), parameter :: subnam = 'cntc_setTrackDimensions'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setTrackDimensions
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, -1, 1, -1, subnam, ierror)
   if (ierror.lt.0) return

   ! check value of ZTRACK digit: 0--3 or 32--33

   if (.not.((ztrack.ge.0 .and. ztrack.le.3) .or. (ztrack.ge.32 .and. ztrack.le.33))) then
      write(bufout,'(2a,i4,a)') trim(pfx_str(subnam,ire,-1)),' method Z=', ztrack,' does not exist.'
      call write_log(1, bufout)
      return
   endif

   ! set minimum / maximum number of parameters for each option

   nparam_loc( 0, 1:2) = (/  0,  0 /)
   nparam_loc( 1, 1:2) = (/  5,  6 /)
   nparam_loc( 2, 1:2) = (/  6,  6 /)
   nparam_loc( 3, 1:2) = (/ 11, 12 /)
   nparam_loc(32, 1:2) = (/  4,  5 /)
   nparam_loc(33, 1:2) = (/  4,  5 /)

   ! check number of parameters supplied

   if (nparam.lt.nparam_loc(ztrack,1) .or. nparam.gt.nparam_loc(ztrack,2)) then
      if (nparam_loc(ztrack,1).eq.nparam_loc(ztrack,2)) then
         write(bufout,'(2a,3(i2,a))') trim(pfx_str(subnam,ire,-1)),' method Z=', ztrack,' needs ',      &
                nparam_loc(ztrack,1),', got ',nparam,' parameters.'
   else
         write(bufout,'(2a,4(i2,a))') trim(pfx_str(subnam,ire,-1)),' method Z=',ztrack,                 &
             ' needs ', nparam_loc(ztrack,1), '-- ', nparam_loc(ztrack,2),', got ',nparam,' parameters.'
   endif
      call write_log(1, bufout)
      return
   endif

   ! set reference to the active rail in the current configuration

   associate(my_rail  => wtd%trk%rai)

   ! store the supplied values

   if (ztrack.le.0) then
      my_ic%ztrack = 0
   else
   my_ic%ztrack = 3
   endif

   if (ztrack.eq.1) then

      wtd%trk%gauge_height    = params(1) * my_scl%len
      if (wtd%trk%gauge_height.gt.0d0) then
         wtd%trk%gauge_seqnum = nint(params(2))
         wtd%trk%track_gauge  = params(3) * my_scl%len
      else
         wtd%trk%rail_y0      = params(2) * my_scl%len
         wtd%trk%rail_z0      = params(3) * my_scl%len
      endif
      wtd%trk%cant_angle      = params(4) * my_scl%angle
      wtd%trk%nom_radius      = params(5) * my_scl%len
      if (nparam.ge.6) then
         wtd%trk%track_curv   = params(6) / my_scl%len
      else
         wtd%trk%track_curv   = 0d0     ! default curvature = 0: tangent track
      endif

   elseif (ztrack.eq.2) then

      my_rail%dy             = params(1) * my_scl%len
      my_rail%dz             = params(2) * my_scl%len
      my_rail%roll           = params(3) * my_scl%angle
      my_rail%vy             = params(4) * my_scl%veloc
      my_rail%vz             = params(5) * my_scl%veloc
      my_rail%vroll          = params(6) * my_scl%angle

   elseif (ztrack.eq.3) then

      wtd%trk%gauge_height    = params( 1) * my_scl%len
      if (wtd%trk%gauge_height.gt.0d0) then
         wtd%trk%gauge_seqnum = nint(params( 2))
         wtd%trk%track_gauge  = params( 3) * my_scl%len
      else
         wtd%trk%rail_y0      = params( 2) * my_scl%len
         wtd%trk%rail_z0      = params( 3) * my_scl%len
      endif
      wtd%trk%cant_angle      = params( 4) * my_scl%angle
      wtd%trk%nom_radius      = params( 5) * my_scl%len
      if (nparam.ge.12) then
         wtd%trk%track_curv   = params( 6) / my_scl%len
         iofs = 6
      else
         wtd%trk%track_curv   = 0d0     ! default curvature = 0: tangent track
         iofs = 5
      endif
      my_rail%dy              = params(iofs+1) * my_scl%len
      my_rail%dz              = params(iofs+2) * my_scl%len
      my_rail%roll            = params(iofs+3) * my_scl%angle
      my_rail%vy              = params(iofs+4) * my_scl%veloc
      my_rail%vz              = params(iofs+5) * my_scl%veloc
      my_rail%vroll           = params(iofs+6) * my_scl%angle

   elseif (ztrack.eq.32 .or. ztrack.eq.33) then

      wtd%trk%ky_rail         = params(1) * my_scl%forc / my_scl%len
      wtd%trk%fy_rail         = params(2) * my_scl%forc
      wtd%trk%kz_rail         = params(3) * my_scl%forc / my_scl%len
      wtd%trk%fz_rail         = params(4) * my_scl%forc
      if (nparam.ge.5) then
         wtd%trk%dystep0      = params(5) * my_scl%len
      else
         wtd%trk%dystep0      = 1d0
      endif

   endif

   if (idebug.ge.3 .and. ztrack.ge.32) then
      write(bufout,'(2a,i2,2(a,g12.4),a)') trim(pfx_str(subnam,ire,-1)),' Z=',ztrack, ', ky=',          &
                wtd%trk%ky_rail,' [N/mm], fy_rail=',wtd%trk%fy_rail, ' [N]'
      call write_log(1, bufout)
   endif

   end associate
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setTrackDimensions

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setWheelsetDimensions(ire, ewheel, nparam, params) &
   bind(c,name=CNAME_(cntc_setwheelsetdimensions))
!--function: set the wheelset description for a wheel-rail contact problem
!    E=0-2, 4: no new geometry         params = []
!    E=3, 5:   new wheelset geometry   params = [fbdist, fbpos, nomrad, {ytape}]
!
!  dimensions:  fbdist, fbpos, nomrad, ytape [length]
!--category: 2, "m=1 only, wtd":    available for module 1 only, working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: ewheel         ! control digit E of the computation
   integer,      intent(in) :: nparam         ! number of parameters provided
   real(kind=8), intent(in) :: params(nparam) ! parameters depending on method that is used
!--local variables:
   integer, parameter  :: nparam_loc(0:5) = (/ 0, 0, 0, 3, 0, 3 /)
   integer             :: ierror
   logical             :: zerror
   character(len=*), parameter :: subnam = 'cntc_setWheelsetDimensions'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setWheelsetDimensions
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, -1, 1, -1, subnam, ierror)
   if (ierror.lt.0) return

   if (ewheel.lt.0 .or. ewheel.gt.5) then
      write(bufout,'(2a,i4,a)') trim(pfx_str(subnam,ire,-1)),' method E=', ewheel,' does not exist.'
      call write_log(1, bufout)
      return
   endif

   if (ewheel.eq.3 .or. ewheel.eq.5) then       ! optional parameter {ytape}
      if (nparam.lt.3 .or. nparam.gt.4) then
         write(bufout,'(2a,i1,a)') trim(pfx_str(subnam,ire,-1)), ' method E=', ewheel,                  &
                   ' needs 3 or 4 parameters.'
         call write_log(1, bufout)
         return
      endif
   else
      if (nparam.ne.nparam_loc(ewheel)) then    ! no optional parameter, using nparam_loc
         write(bufout,'(2a,2(i1,a))') trim(pfx_str(subnam,ire,-1)), ' method E=', ewheel,' needs ',     &
                   nparam_loc(ewheel),' parameters.'
         call write_log(1, bufout)
         return
      endif
   endif

   ! store E-digit - maintaining flexible wheelset deviations

   if (my_ic%ewheel.le.3) then
      my_ic%ewheel = ewheel
   elseif (my_ic%ewheel.eq.4) then
      my_ic%ewheel = 5
   endif

   if (ewheel.eq.3 .or. ewheel.eq.5) then

      wtd%ws%flback_dist  = params(1) * my_scl%len
      wtd%ws%flback_pos   = params(2) * my_scl%len
      wtd%ws%nom_radius   = params(3) * my_scl%len
      if (nparam.ge.4) then
         wtd%ws%ytape     = params(4) * my_scl%len
      else
         wtd%ws%ytape     = 1000d0
      endif

      zerror = .not.check_range ('NOMRAD', wtd%ws%nom_radius, 1d-3, 1d20)
      wtd%ws%nom_radius   = max(1d-3, min(1d20, wtd%ws%nom_radius))

   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setWheelsetDimensions

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setWheelsetPosition(ire, ewheel, nparam, params) &
   bind(c,name=CNAME_(cntc_setwheelsetposition))
!--function: set the wheelset position data for a wheel-rail contact problem
!    1-5: new wheelset position     params = [ s_ws, y_ws, z_ws, roll, yaw, pitch ]
!       s_ws, y_ws, z_ws : [length]    roll, yaw, pitch : [angle]
!       z_ws is ignored when N=1
!--category: 2, "m=1 only, wtd":    available for module 1 only, working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: ewheel         ! control digit E of the computation
   integer,      intent(in) :: nparam         ! number of parameters provided
   real(kind=8), intent(in) :: params(nparam) ! parameters depending on method that is used
!--local variables:
   integer, parameter  :: nparam_loc(1:5) = (/ 6, 6, 6, 6, 6 /)
   integer             :: ierror
   character(len=*), parameter :: subnam = 'cntc_setWheelsetPosition'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setWheelsetPosition
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, -1, 1, -1, subnam, ierror)
   if (ierror.lt.0) return

   if (ewheel.lt.1 .or. ewheel.gt.5) then
      write(bufout,'(2a,i4,a)') trim(pfx_str(subnam,ire,-1)),' method',ewheel,' does not exist.'
      call write_log(1, bufout)
      return
   endif

   if (nparam.ne.nparam_loc(ewheel)) then
      write(bufout,'(2a,2(i1,a))') trim(pfx_str(subnam,ire,-1)),' method E=',ewheel,' needs ',          &
                nparam_loc(ewheel),' parameters.'
      call write_log(1, bufout)
      return
   endif

   ! store E-digit - maintaining velocity, geometry, flexible wheelset deviations

   if (my_ic%ewheel.le.0) my_ic%ewheel = 1

   ! store data

   if (ewheel.ge.1 .and. ewheel.le.5) then

      wtd%ws%s     = params(1) * my_scl%len
      wtd%ws%y     = params(2) * my_scl%len
      wtd%ws%z     = params(3) * my_scl%len
      wtd%ws%roll  = params(4) * my_scl%angle
      wtd%ws%yaw   = params(5) * my_scl%angle
      wtd%ws%pitch = params(6) * my_scl%angle

   endif

   if (idebug.ge.2) then
      write(bufout,'(3a,2f8.3,a)') trim(pfx_str(subnam,ire,-1)),' position=', fmt_gs(12,6,4,wtd%ws%s),  &
                wtd%ws%y, wtd%ws%z, ' [mm]'
      call write_log(1, bufout)
      write(bufout,'(2a,3f10.6,a)') trim(pfx_str(subnam,ire,-1)),' orientation=', wtd%ws%roll,          &
                wtd%ws%yaw, wtd%ws%pitch, ' [rad]'
      call write_log(1, bufout)
   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setWheelsetPosition

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setWheelsetVelocity(ire, ewheel, nparam, params) &
   bind(c,name=CNAME_(cntc_setwheelsetvelocity))
!--function: set the wheelset velocity data for a wheel-rail contact problem
!      1: no new wheelset velocity    params = [ ]
!    2-5:    new wheelset velocity    params = [ vx, vy, vz, vroll, vyaw, vpitch ]
!            for roller-rigs vx is replaced by rpitch (C1=4,5)
!            position increments v * dt are used in transient shifts (T=1)
!       vx, vy, vz    : [veloc]   vroll, vyaw, vpitch, rpitch  : [ang.veloc]
!       shft_sws--zws : [length]  shft_rol, shft_yaw, shft_pit : [angle]
!--category: 2, "m=1 only, wtd":    available for module 1 only, working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: ewheel         ! control digit E of the computation
   integer,      intent(in) :: nparam         ! number of parameters provided
   real(kind=8), intent(in) :: params(nparam) ! parameters depending on method that is used
!--local variables:
   integer, parameter  :: nparam_loc(1:5) = (/ 0, 6, 6, 6, 6 /)
   integer             :: ewheel_loc, ierror
   logical             :: is_roller
   character(len=*), parameter :: subnam = 'cntc_setWheelsetVelocity'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setWheelsetVelocity
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, -1, 1, -1, subnam, ierror)
   if (ierror.lt.0) return

   ! map previous imeth=1 to E=2

   ewheel_loc = ewheel

   if (ewheel.eq.1 .and. nparam.eq.nparam_loc(2)) then
      write(bufout,'(2a,i4,a)') trim(pfx_str(subnam,ire,-1)),' method E=', ewheel,                      &
                ' is obsolete for setting velocity data. Using E=2 instead.'
      call write_log(1, bufout)
      ewheel_loc = 2
   endif

   if (ewheel_loc.lt.2 .or. ewheel_loc.gt.5) then
      write(bufout,'(2a,i4,a)') trim(pfx_str(subnam,ire,-1)),' method E=', ewheel_loc,' does not exist.'
      call write_log(1, bufout)
      return
   endif

   if (nparam.ne.nparam_loc(ewheel_loc)) then
      write(bufout,'(2a,2(i1,a))') trim(pfx_str(subnam,ire,-1)),' method ewheel=',ewheel_loc,' needs ', &
                nparam_loc(ewheel_loc),' parameters.'
      call write_log(1, bufout)
      return
   endif

   ! store E-digit - maintaining geometry, flexible wheelset deviations

   if (my_ic%ewheel.le.1) my_ic%ewheel = ewheel_loc

   ! store data - also for E=1 for compatibility reasons

   if (ewheel_loc.ge.2 .and. ewheel_loc.le.5) then

      is_roller = (wtd%ic%config.ge.4)
      if (is_roller) then
         wtd%trk%vpitch_rol = params(1) * my_scl%angle ! angle/time
         wtd%ws%vs          = 0d0
      else
         wtd%trk%vpitch_rol = 0d0
         wtd%ws%vs          = params(1) * my_scl%veloc
      endif
      wtd%ws%vy     = params(2) * my_scl%veloc
      wtd%ws%vz     = params(3) * my_scl%veloc
      wtd%ws%vroll  = params(4) * my_scl%angle ! angle/time
      wtd%ws%vyaw   = params(5) * my_scl%angle ! angle/time
      wtd%ws%vpitch = params(6) * my_scl%angle ! angle/time

   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setWheelsetVelocity

!------------------------------------------------------------------------------------------------------------

subroutine cntc_setWheelsetFlexibility(ire, ewheel, nparam, params) &
   bind(c,name=CNAME_(cntc_setwheelsetflexibility))
!--function: set the description of wheelset flexibilities for a wheel-rail contact problem
!    1--3: no wheelset flexibility               params = []
!     4,5: new wheelset flexibility parameters   params = [dxwhl, dywhl, dzwhl, drollw, dyaww, dpitchw,
!                                                          vxwhl, vywhl, vzwhl, vrollw, vyaww, vpitchw],
!
!    dimensions:   dxwhl, dywhl, dzwhl [length],  drollw, dyaww, dpitchw [angle]
!                  vxwhl, vywhl, vzwhl [veloc],   vrollw, vyaww, vpitchw [ang.veloc]
!--category: 2, "m=1 only, wtd":    available for module 1 only, working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: ewheel         ! control digit E of the computation
   integer,      intent(in) :: nparam         ! number of parameters provided
   real(kind=8), intent(in) :: params(nparam) ! parameters depending on method that is used
!--local variables:
   integer, parameter  :: nparam_loc(1:5) = (/ 0, 0, 0, 12, 12 /)
   integer             :: ierror
   character(len=*), parameter :: subnam = 'cntc_setWheelsetFlexibility'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_setWheelsetFlexibility
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, -1, 1, -1, subnam, ierror)
   if (ierror.lt.0) return

   ! check value of E-digit

   if (ewheel.lt.1 .or. ewheel.gt.5) then
      write(bufout,'(2a,i4,a)') trim(pfx_str(subnam,ire,-1)),' method E=', ewheel,' does not exist.'
      call write_log(1, bufout)
      return
   endif

   ! check number of parameters provided

   if ((nparam_loc(ewheel).eq.0 .and. nparam.ge.2) .or.                                                 &
       (nparam_loc(ewheel).ne.0 .and. nparam.ne.nparam_loc(ewheel))) then
      write(bufout,'(2a,2(i1,a))') trim(pfx_str(subnam,ire,-1)), ' method E=', ewheel,' needs ',        &
                nparam_loc(ewheel),' parameters.'
      call write_log(1, bufout)
      return
   endif

   ! store E-digit - maintaining position/velocity/geometry specification
   !                 method: 6 x 6 table desired value for (old-ewheel, new-ewheel)
   ! e1, ewheel   0 = maintain wheel-set geometry, profile, position and velocity;
   !              1 = read new position data;
   !              2 = read new position and velocity data;
   !              3 = read new geometry, wheel profile for current side, position and velocity data;
   !              4 = as 2, including flexible wheelset deviations for current side of wheelset
   !              5 = as 3, including flexible wheelset deviations for current side of wheelset

   if (ewheel.le.3) then        ! no flexibilities, keep position, velocity, profile(s)
      if (my_ic%ewheel.le.3) then
         ! keep existing my_ic%ewheel
      elseif (my_ic%ewheel.le.5) then
         my_ic%ewheel = my_ic%ewheel - 2
      endif
   elseif (ewheel.ge.5) then    !    flexibilities + position, velocity, profile(s)
      my_ic%ewheel = ewheel
   else                         ! new ewheel 5: keep profile specification
      if (my_ic%ewheel.eq.3 .or. my_ic%ewheel.eq.5) then
         my_ic%ewheel = 5
      else
         my_ic%ewheel = ewheel
      endif
   endif

   ! store parameters provided

   associate(my_wheel => wtd%ws%whl)
   if (ewheel.ge.1 .and. ewheel.le.3) then

      my_wheel%dx     = 0d0
      my_wheel%dy     = 0d0
      my_wheel%dz     = 0d0
      my_wheel%droll  = 0d0
      my_wheel%dyaw   = 0d0
      my_wheel%dpitch = 0d0
      my_wheel%vx     = 0d0
      my_wheel%vy     = 0d0
      my_wheel%vz     = 0d0
      my_wheel%vroll  = 0d0
      my_wheel%vyaw   = 0d0
      my_wheel%vpitch = 0d0

   elseif (ewheel.eq.4 .or. ewheel.eq.5) then

      my_wheel%dx     = params( 1) * my_scl%len
      my_wheel%dy     = params( 2) * my_scl%len
      my_wheel%dz     = params( 3) * my_scl%len
      my_wheel%droll  = params( 4) * my_scl%angle
      my_wheel%dyaw   = params( 5) * my_scl%angle
      my_wheel%dpitch = params( 6) * my_scl%angle
      my_wheel%vx     = params( 7) * my_scl%veloc
      my_wheel%vy     = params( 8) * my_scl%veloc
      my_wheel%vz     = params( 9) * my_scl%veloc
      my_wheel%vroll  = params(10) * my_scl%angle
      my_wheel%vyaw   = params(11) * my_scl%angle
      my_wheel%vpitch = params(12) * my_scl%angle

   endif
   end associate

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_setWheelsetFlexibility

!------------------------------------------------------------------------------------------------------------

subroutine subs_addBlock(ire, icp, iblk, isubs, npx, npy, npz, xparam, yparam, zparam) &
   bind(c,name=CNAME_(subs_addblock))
!--function: set the parameters for a block of points for the subsurf. calculation for a contact problem
!  isubs:
!    1: xparam = [                 ]  yparam = [                 ]  zparam = [ NZ, ZL, DZ    ]
!    2: xparam = [ IXL, IXINC, IXH ]  yparam = [ IYL, IYINC, IYH ]  zparam = [ NZ, ZL, DZ    ]
!    3: xparam = [ IX(i), i=1:npx  ]  yparam = [ IY(j), j=1:npy  ]  zparam = [ NZ, ZL, DZ    ]
!    5: xparam = [                 ]  yparam = [                 ]  zparam = [ Z(k), k=1:npz ]
!    6: xparam = [ IXL, IXINC, IXH ]  yparam = [ IYL, IYINC, IYH ]  zparam = [ Z(k), k=1:npz ]
!    7: xparam = [ IX(i), i=1:npx  ]  yparam = [ IY(j), j=1:npy  ]  zparam = [ Z(k), k=1:npz ]
!    9: xparam = [ X(i),  i=1:npx  ]  yparam = [ Y(j),  j=1:npy  ]  zparam = [ Z(k), k=1:npz ]
!       ix, iy [-],    x, y, z, zl, dz [length]
!--category: 7, "m=any, wtd/cp":    available for modules 1 and 3, in module 1 working on wtd or cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in) :: ire            ! result element ID
   integer,      intent(in) :: icp            ! contact problem ID
   integer,      intent(in) :: iblk           ! subsurface stress block number
   integer,      intent(in) :: isubs          ! type of block specification
   integer,      intent(in) :: npx, npy, npz  ! number of parameters given per dimension (dependent on isubs)
   real(kind=8), intent(in) :: xparam(npx)    ! dependent on type of specification [-/length]
   real(kind=8), intent(in) :: yparam(npy)    ! dependent on type of specification [-/length]
   real(kind=8), intent(in) :: zparam(npz)    ! dependent on type of specification [-/length]
!--local variables:
   integer,          parameter :: npx_loc(1:9) = (/ 0, 3, -1, 0,  0,  3, -1, 0, -1 /),                   &
                                  npy_loc(1:9) = (/ 0, 3, -1, 0,  0,  3, -1, 0, -1 /),                   &
                                  npz_loc(1:9) = (/ 3, 3,  3, 0, -1, -1, -1, 0, -1 /)
   integer                     :: ierror, imodul, jblk
   logical                     :: lchanged
   character(len=*), parameter :: subnam = 'subs_addBlock'
#ifdef _WIN32
!dec$ attributes dllexport :: subs_addBlock
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 0, subnam, ierror)
   if (ierror.lt.0) return

   ! check 0 <= iblk <= min(max_block, nblock+1)

   if (iblk.lt.0 .or. iblk.gt.min(MXBLCK, my_subs%nblock+1)) then
      write(bufout,'(2a,i4,a)') trim(pfx_str(subnam,ire,icp)),' invalid block', iblk,                   &
                ', must be 0 <= IBLK <= NBLOCK+1.'
      call write_log(1, bufout)
      return
   endif

   ! check isubs = 1,2,3, 5,6,7 or 9

   if (isubs.lt.1 .or. isubs.eq.4 .or. isubs.eq.8 .or. isubs.gt.9) then
      write(bufout,'(2a,i4,a)') trim(pfx_str(subnam,ire,icp)),' method',isubs,' does not exist.'
      call write_log(1, bufout)
      return
   endif

   ! check number of parameters npx, npy, npz needed for isubs

   if (npx_loc(isubs).gt.0 .and. npx.ne.npx_loc(isubs)) then
      write(bufout,'(2a,2(i1,a))') trim(pfx_str(subnam,ire,icp)),' method ISUBS=',isubs,' needs NPX=', &
                npx_loc(isubs),' parameters.'
      call write_log(1, bufout)
      return
   endif
   if (npy_loc(isubs).gt.0 .and. npy.ne.npy_loc(isubs)) then
      write(bufout,'(2a,2(i1,a))') trim(pfx_str(subnam,ire,icp)),' method ISUBS=',isubs,' needs NPY=', &
                npy_loc(isubs),' parameters.'
      call write_log(1, bufout)
      return
   endif
   if (npz_loc(isubs).gt.0 .and. npz.ne.npz_loc(isubs)) then
      write(bufout,'(2a,2(i1,a))') trim(pfx_str(subnam,ire,icp)),' method ISUBS=',isubs,' needs NPZ=', &
                npz_loc(isubs),' parameters.'
      call write_log(1, bufout)
      return
   endif

   ! clear existing blocks with numbers iblk and higher

   do jblk = max(1,iblk), my_subs%nblock
      call subsblk_destroy( my_subs%blocks(jblk) )
   enddo
   my_subs%nblock = max(0, iblk - 1)

   imodul = ire_module(ix_reid(ire))

   if (iblk.le.0) then

      ! new block = 0 : just clear the list of blocks

      ! icp>0: reset flag 'has_own_subs' for contact patch icp

      if (imodul.eq.1 .and. icp.ge.1) wtd%allcps(icp)%cp%has_own_subs = .false.

   else

      ! icp>0: set flag 'has_own_subs' for contact patch icp

      if (imodul.eq.1 .and. icp.ge.1) wtd%allcps(icp)%cp%has_own_subs = .true.

      ! increment number of blocks, set reference to block data

      my_subs%nblock = my_subs%nblock + 1
      associate(b => my_subs%blocks(my_subs%nblock))

      ! store data for block iblk

      b%isubs = isubs

      ! store horizontal part for isubs = 1 - 7

      if (isubs.eq.1 .or. isubs.eq.5) then

         b%ixl_inp = 1
         b%ixinc   = 1
         b%ixh_inp = 99999
         b%iyl_inp = 1
         b%iyinc   = 1
         b%iyh_inp = 99999

      elseif (isubs.eq.2 .or. isubs.eq.6) then

         b%ixl_inp = max(nint(xparam(1)),      1)
         b%ixinc   = max(nint(xparam(2)),      1)
         b%ixh_inp = min(nint(xparam(3)), 999999)
         b%iyl_inp = max(nint(yparam(1)),      1)
         b%iyinc   = max(nint(yparam(2)),      1)
         b%iyh_inp = min(nint(yparam(3)), 999999)

      elseif (isubs.eq.3 .or. isubs.eq.7) then

         ! unique_positive: sort in ascending order, copy/maintain unique positive values

         b%nx_inp  = max(1, min(npx, MXSUBS))
         b%ny_inp  = max(1, min(npy, MXSUBS))
         call reallocate_arr(b%ixlist_inp, b%nx_inp)
         call reallocate_arr(b%iylist_inp, b%ny_inp)
         call unique_positive('IX', b%nx_inp, xparam, b%ixlist_inp, lchanged)
         call unique_positive('IY', b%ny_inp, yparam, b%iylist_inp, lchanged)

      endif

      ! store vertical part for isubs = 1 - 7

      if (isubs.ge.1 .and. isubs.le.3) then

         b%nz    = nint(zparam(1))
         b%nz    = max(1, min(b%nz, MXSUBS))
         b%zl    =     zparam(2) * my_scl%len
         b%dz    = max(zparam(3) * my_scl%len, 1d-9)

      elseif (isubs.ge.5 .and. isubs.le.7) then

         b%nz    = max(1, min(npz, MXSUBS))
         call reallocate_arr(b%z, b%nz)
         b%z(1:b%nz) = zparam(1:b%nz) * my_scl%len
         call bubsrt_real (b%nz, b%z)

      endif

      ! isubs = 9: explicit lists x, y, z

      if (isubs.eq.9) then

         b%nx_inp = max(1, min(npx, MXSUBS))
         b%ny_inp = max(1, min(npy, MXSUBS))
         b%nz     = max(1, min(npz, MXSUBS))

         call reallocate_arr(b%x, b%nx_inp)
         call reallocate_arr(b%y, b%ny_inp)
         call reallocate_arr(b%z, b%nz)

         b%x(1:b%nx_inp) = xparam(1:b%nx_inp) * my_scl%len
         b%y(1:b%ny_inp) = yparam(1:b%ny_inp) * my_scl%len
         b%z(1:b%nz)     = zparam(1:b%nz)     * my_scl%len

         call bubsrt_real (b%nx_inp, b%x)
         call bubsrt_real (b%ny_inp, b%y)
         call bubsrt_real (b%nz,     b%z)

      endif

      end associate
   endif ! iblk<=0

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine subs_addBlock

!------------------------------------------------------------------------------------------------------------

subroutine cntc_calculate(ire, icp, ierror) &
   bind(c,name=CNAME_(cntc_calculate))
!--function: perform the actual calculation for a contact problem
!--category: 5, "m=any, wtd":       available for modules 1 and 3, in module 1 working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   integer,      intent(out) :: ierror        ! error code of CONTACT calculation
!--local variables:
   character(len=*), parameter   :: subnam = 'cntc_calculate'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_calculate
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, -1, subnam, ierror)

   if (.not.my_license%is_valid) ierror = CNTC_err_allow

   call lock_contact_problem(ire, icp, ierror)

   if (ierror.eq.0) then
      if (icp.le.0) then
         call cntc_calculate1(ire, ierror)
      else
         call cntc_calculate3(ire, icp, ierror)
      endif
   endif

   call free_contact_problem(ire, icp)

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_calculate

!------------------------------------------------------------------------------------------------------------

subroutine cntc_calculate1(ire, ierror) &
   bind(c,name=CNAME_(cntc_calculate1))
!--function: perform the actual calculation for a wheel-rail contact problem
!--category: 2, "m=1 only, wtd":    available for module 1 only, working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(out) :: ierror        ! error code of CONTACT calculation
!--local variables:
   character(len=*),   parameter :: subnam = 'cntc_calculate1'
   integer                       :: is_wheel, i_ftype
   character(len=256)            :: fname
   character(len=len(bufout)-52) :: tmpbuf
   integer                       :: itimer
   logical                       :: is_ok
   real(kind=8)                  :: dx_prv, ds_prv, dx_inp, ds_inp
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_calculate1
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, -1, 1, -1, subnam, ierror)
   if (ierror.lt.0) return

   ! set reference to the active wheel and rail in the current configuration

   associate(my_wheel => wtd%ws%whl, my_rail  => wtd%trk%rai)

   if (idebug.ge.2) then
      write(bufout,'(2a)') trim(pfx_str(subnam,ire,-1)),' checking problem specification...'
      call write_log(1, bufout)
   endif

   ! We should check whether all data are provided: mx, my, etc.?
   ! We should check whether the appropriate updates are carried out: scale, ...?

   if (my_rail%prr%ierror.ne.0 .or. my_wheel%prw%ierror.ne.0) then
      if (my_rail%prr%ierror.eq.1) then
         call write_log(' ERROR. The rail profile has not been set.')
      elseif (my_rail%prr%ierror.lt.0) then
         write(bufout,'(a,i4,a)') ' ERROR. The rail profile could not be found or processed (',         &
                my_rail%prr%ierror,')'
         call write_log(1, bufout)
      endif
      if (my_wheel%prw%ierror.eq.1) then
         call write_log(' ERROR. The wheel profile has not been set.')
      elseif (my_wheel%prw%ierror.lt.0) then
         write(bufout,'(a,i4,a)') ' ERROR. The wheel profile could not be found or processed (',        &
                my_wheel%prw%ierror,')'
         call write_log(1, bufout)
      endif
      ierror = CNTC_err_profil
      return
   endif

   ! check that variable profiles are used with absolute rail placement

   call profile_get_filetype(my_rail%prr%fname, is_wheel, i_ftype, idebug)
   if (i_ftype.eq.FTYPE_SLICES .and. wtd%trk%gauge_height.gt.0d0) then
      ierror = CNTC_err_profil
      write(bufout, '(a,f8.2,a,/,15x,3a)') ' Input: ERROR. The gauge point computation (GAUGHT =',      &
                wtd%trk%gauge_height, ') cannot be used when using a', 'variable profile ("',           &
                trim(my_rail%prr%fname),'").'
      call write_log(1, bufout)
   endif

   ! TODO: signal errors in friction specification
   if (wtd%fric%nvf.le.0) then
      call write_log(' An error occurred in the specification of friction law parameters')
      ierror = CNTC_err_frclaw
      return
   endif

   ! in transient cases, check element sizes of new/previous case

   if (wtd%meta%ncase.gt.1 .and. wtd%numcps.ge.1 .and. (wtd%ic%tang.eq.1 .or. wtd%ic%tang.eq.2)) then
      associate(cp => wtd%allcps(1)%cp)
      dx_prv = cp%gd%cgrid_cur%dx
      ds_prv = cp%gd%cgrid_cur%dy
      dx_inp = wtd%discr%dx
      ds_inp = wtd%discr%ds
      is_ok = (dx_inp.ge.1d-10 .and. abs(dx_prv-dx_inp).lt.1d-4*min(dx_prv, dx_inp) .and.               &
               ds_inp.ge.1d-10 .and. abs(ds_prv-ds_inp).lt.1d-4*min(ds_prv, ds_inp))
      if (.not.is_ok) then
         ierror = CNTC_err_discr
         write(bufout,'(a,i3,a,i1,2(a,f6.3),a)') ' ERROR ',ierror,': transient contact (T=',            &
                wtd%ic%tang,') needs constant DX, DY (',dx_prv,',',ds_prv,') during the simulation'
         call write_log(1, bufout)
         return
      endif
      end associate
   endif
 
   if (idebug.ge.5) then
      write(bufout,'(a,i3,3(a,f6.3))') 'rail profile:  err_hnd=',my_rail%prr%err_hnd,', max_omit=',     &
         my_rail%prr%f_max_omit,', zig_thrs=',my_rail%prr%zig_thrs,', kink_high=',my_rail%prr%kink_high
      call write_log(1, bufout)
      write(bufout,'(a,i3,3(a,f6.3))') 'wheel profile: err_hnd=',my_wheel%prw%err_hnd,', max_omit=',    &
         my_wheel%prw%f_max_omit,', zig_thrs=',my_wheel%prw%zig_thrs,', kink_high=',my_wheel%prw%kink_high
      call write_log(1, bufout)
   endif

   if (idebug.ge.1) then
      write(bufout,'(2a)') trim(pfx_str(subnam,ire,-1)),' calculating...'
      call write_log(1, bufout)
   endif

   if ((my_ic%output_surf.ge.1 .or. my_ic%flow.ge.1) .and. out_open.eq.1) then
      if (my_meta%irun.eq.0 .and. my_meta%tim.eq.0d0) then
         write (lout,'(/,a,i0,a,i0)') ' Case ',my_meta%ncase,' for w/r contact on result element ',ire
      elseif (my_meta%irun.eq.0) then
         write (lout,'(/,a,i0,a,i0,a,f16.10)') ' Case ',my_meta%ncase,' for w/r contact on result element ', &
                ire,', t=',my_meta%tim
      else
         write (lout,'(/,5(a,i0))') ' Case ',my_meta%ncase,' for w/r contact on result element ', ire,  &
                ', run ',my_meta%irun,', axle ', my_meta%iax,', side ',my_meta%iside
      endif
   endif

   ! Open .inp-file for writing (output-folder) if not done so before

   if (my_ic%wrtinp.ge.1 .and. linp.ge.1 .and. inp_open.eq.0) then
      fname = trim(caddon_outdir) // trim(caddon_expnam) // '.inp'
      open(unit=linp, file=fname, action='write', err=998)
      inp_open = 1
      goto 999

      ! error handling: file open failed
 998  continue
         inp_open = -1
         if (idebug.ge.1) then
            tmpbuf = fname ! avoid potential overflow of bufout
            write(bufout,'(/,2a,/,36x,3a)') trim(pfx_str(subnam,ire,-1)),                               &
                ' ERROR: could not open file', '"',trim(tmpbuf), '" for writing.'
            call write_log(3, bufout)
         endif
 999  continue
   endif

   ! Optionally write description of the case to the .inp-file

   if (my_ic%wrtinp.ge.1 .and. inp_open.ge.0) then
      call timer_start(itimer_files)
      if (my_ic%return.eq.0 .or. my_ic%return.eq.2 .and. my_meta%irun.eq.0) then
         write (linp,'(a,i0,a,i3)')    '% case ',my_meta%ncase,' for result element',ire
      elseif (my_ic%return.eq.0 .or. my_ic%return.eq.2) then
         write (linp,'(a,i0,4(a,i4))') '% case ',my_meta%ncase,' for result element',ire, ': run',      &
                my_meta%irun,', axle', my_meta%iax,', side',my_meta%iside
      elseif (my_meta%irun.eq.0 .and. my_meta%tim.eq.0d0) then
         write (linp,'(a,i0,a,i3)')         ' 1  MODULE   % case ',my_meta%ncase,' for result element',ire
      elseif (my_meta%irun.eq.0) then
         write (linp,'(a,i0,a,i3,a,f16.8)') ' 1  MODULE   % case ',my_meta%ncase,' for result element', &
                ire,', t=',my_meta%tim
      else
         write (linp,'(a,i0,4(a,i4))')      ' 1  MODULE   % case ',my_meta%ncase,' for result element', &
                ire, ': run', my_meta%irun,', axle', my_meta%iax,', side',my_meta%iside
      endif

      call wr_write_inp(-1, wtd)
      call timer_stop(itimer_files)
   endif

   ! Solve the contact-problem

   if (wtd%ic%return.le.1) then
      itimer = cntc_timer_num(ire, -1)
      call timer_start(itimer)
      call timer_start(itimer_wrprof)
      call wr_contact (wtd, ierror)
      call timer_stop(itimer_wrprof)
      call timer_stop(itimer)
   endif

   ! Check for errors in the solution

   if (ierror.lt.0) then
      if (idebug.ge.0) then
         write(bufout,'(2a,i3,a)') trim(pfx_str(subnam,ire,-1)),' no solution for contact problem (', &
                ierror,').'
         call write_log(1, bufout)
      endif
   endif

   ! Increment the number of cases computed for this (ire)

   my_meta%ncase = my_meta%ncase + 1

   if (idebug.ge.1) then
      write(bufout,'(2a)')       trim(pfx_str(subnam,ire,-1)),' done...'
      call write_log(1, bufout)
   endif
   end associate
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_calculate1

!------------------------------------------------------------------------------------------------------------

subroutine cntc_calculate3(ire, icp, ierror) &
   bind(c,name=CNAME_(cntc_calculate3))
!--function: perform the actual calculation for a contact problem
!--category: 1, "m=3 only, cp":     available for module 3 only, working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   integer,      intent(out) :: ierror        ! error code of CONTACT calculation
!--local variables:
   character(len=*), parameter   :: subnam = 'cntc_calculate3'
   integer                       :: itimer
   character(len=256)            :: fname
   character(len=len(bufout)-52) :: tmpbuf
!--subroutines called:
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_calculate3
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 3, 1, subnam, ierror)
   if (ierror.lt.0) return

   ! report on the settings for the contact patch

   if (idebug.ge.2) then
      write(bufout,'(2a,f9.1,f6.3)') trim(pfx_str(subnam,ire,icp)),' g,nu1=',                           &
                my_mater%gg(1), my_mater%poiss(1)
      call write_log(1, bufout)
      write(bufout,'(2a,f9.1,f6.3)') trim(pfx_str(subnam,ire,icp)),' g,nu2=',                           &
                my_mater%gg(2), my_mater%poiss(2)
      call write_log(1, bufout)
      write(bufout,'(2a,i4,2f5.2)')  trim(pfx_str(subnam,ire,icp)),' nvf, fstat, fkin =',               &
                gd%fric%nvf, gd%fric%fstat_arr(1), gd%fric%fkin_arr(1)
      call write_log(1, bufout)
      if (gd%potcon_inp%ipotcn.lt.0) then
         write(bufout,'(2a,2f6.2)')  trim(pfx_str(subnam,ire,icp)),' semi-axes a,b=',                   &
                gd%hertz%aa, gd%hertz%bb
         call write_log(1, bufout)
      else
         write(bufout,'(2a,f6.3,f8.3,f6.3)') trim(pfx_str(subnam,ire,icp)),' step sizes=',              &
                gd%potcon_inp%dx, gd%potcon_inp%dy, my_kin%dq
         call write_log(1, bufout)
         write(bufout,'(2a,f7.2,f8.2)')  trim(pfx_str(subnam,ire,icp)),' corner xl,yl=',                &
                gd%potcon_inp%xl, gd%potcon_inp%yl
         call write_log(1, bufout)
         write(bufout,'(2a,f7.2,f8.2)')  trim(pfx_str(subnam,ire,icp)),' corner xh,yh=',              &
                gd%potcon_inp%xh, gd%potcon_inp%yh
         call write_log(1, bufout)
      endif
      if (my_ic%norm.eq.0) then
         write(bufout,'(2a,f7.4)')   trim(pfx_str(subnam,ire,icp)),' penetration=', my_kin%pen
         call write_log(1, bufout)
      else
         write(bufout,'(2a,f9.1)')   trim(pfx_str(subnam,ire,icp)),' normal force=', my_kin%fntrue
         call write_log(1, bufout)
      endif
      write(bufout,'(2a,3f9.5)')     trim(pfx_str(subnam,ire,icp)),' creepages vx,vy,phi=',           &
                my_kin%cksi, my_kin%ceta, my_kin%cphi
      call write_log(1, bufout)
   endif
   if (idebug.ge.1) then
      write(bufout,'(2a)') trim(pfx_str(subnam,ire,icp)),' calculating...'
      call write_log(1, bufout)
   endif

   ! We should check whether all data are provided: mx, my, etc.?
   ! We should check whether the appropriate updates are carried out: scale, ...?

   ! Open .inp-file for writing (output-folder) if not done so before

   if (my_ic%wrtinp.ge.1 .and. linp.ge.1 .and. inp_open.eq.0) then
      fname = trim(caddon_outdir) // trim(caddon_expnam) // '.inp'
      open(unit=linp, file=fname, action='write', err=998)
      inp_open = 1
      goto 999

      ! error handling: file open failed
 998  continue
         inp_open = -1
         if (idebug.ge.1) then
            tmpbuf = fname ! avoid potential overflow of bufout
            write(bufout,'(/,2a,/,36x,3a)') trim(pfx_str(subnam,ire,icp)),                              &
                ' ERROR: could not open file', '"',trim(tmpbuf), '" for writing.'
            call write_log(3, bufout)
         endif
 999  continue
   endif

   if (my_ic%output_surf.ge.1 .and. out_open.eq.1) then
      write (lout,'(a,i6,2(a,i3))') ' Case',my_meta%ncase,' for contact problem',icp,' on result element',ire
   endif

   ! Solve the contact-problem

   itimer = cntc_timer_num(ire, icp)
   call timer_start(itimer)
   call contac(gd, ierror)
   call timer_stop(itimer)

   ! Optionally write description of the case to the .inp-file

   if (my_ic%wrtinp.ge.1 .and. inp_open.ge.0) then
      call timer_start(itimer_files)
      write (linp,'(a,i7,2(a,i3))') ' 3  MODULE   % case', my_meta%ncase, ' for result element', ire,   &
                ', contact problem',icp
      call wrtinp (-1, my_ic, gd%potcon_cur, my_mater, gd%hertz, gd%geom, gd%fric, my_kin, my_solv,     &
                        gd%outpt1, gd%subs, .true.)
      call timer_stop(itimer_files)
   endif

   ! Increment the number of cases computed for this (ire,icp)

   my_meta%ncase = my_meta%ncase + 1

   ! Check for errors in NORM and TANG or too small pot.contact area

   if (ierror.eq.0) then
      if (my_solv%itnorm.lt.0) then
         ierror = CNTC_err_norm
         if (idebug.ge.0) then
            write(bufout,'(2a)') trim(pfx_str(subnam,ire,icp)),' no convergence in NORM algorithm.'
            call write_log(1, bufout)
         endif
      elseif (my_solv%ittang.lt.0) then
         ierror = CNTC_err_tang
         if (idebug.ge.0) then
            write(bufout,'(2a)') trim(pfx_str(subnam,ire,icp)),' no convergence in TANG algorithm.'
            call write_log(1, bufout)
         endif
      else
         ! Count the number of interior elements at the boundaries of the potential contact
         ierror = eldiv_count_atbnd(gd%outpt1%igs, 0, idebug)
         if (ierror.ge.1 .and. idebug.ge.1) then
            write(bufout,'(2a,i5,2a)') trim(pfx_str(subnam,ire,icp)),' potential contact area too small;', &
                   ierror,' elements next to boundary...'
            call write_log(1, bufout)
         endif
      endif
   endif

   if (idebug.ge.1) then
      write(bufout,'(2a)')       trim(pfx_str(subnam,ire,icp)),' done...'
      call write_log(1, bufout)
   endif
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_calculate3

!------------------------------------------------------------------------------------------------------------

subroutine subs_calculate(ire, icp, ierror) &
   bind(c,name=CNAME_(subs_calculate))
!--function: perform the actual calculation of subsurface stresses for a contact problem
!--category: 7, "m=any, wtd or cp":  available for modules 1 and 3, in module 1 working on single or all cps
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   integer,      intent(out) :: ierror        ! error code of CONTACT calculation
!--local variables:
   character(len=*), parameter   :: subnam = 'subs_calculate'
   integer                       :: imodul, itimer, jcp, jcp0, jcp1
#ifdef _WIN32
!dec$ attributes dllexport :: subs_calculate
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 0, subnam, ierror)
   if (.not.my_license%is_valid) ierror = CNTC_err_allow
   if (ierror.lt.0) return

   ! perform checks for multi-threading, lock contact problem

   call lock_contact_problem(ire, icp, ierror)

   ! skip calculation in case of errors

   if (ierror.eq.0) then

      ! TODO: check that results are available

      imodul = ire_module(ix_reid(ire))

      if (idebug.ge.2) then
         write(bufout,'(2a)') trim(pfx_str(subnam,ire,icp)),' calculating...'
         call write_log(1, bufout)
      endif

      itimer = cntc_timer_num(ire, -1)
      call timer_start(itimer)

      if (imodul.eq.1) then

         ! module 1: use loop, copy input from wtd to gd as needed

         if (icp.le.0) then
            jcp0 = 1
            jcp1 = wtd%numcps
         else
            jcp0 = icp
            jcp1 = icp
         endif

         do jcp = jcp0, jcp1

            associate(cp => wtd%allcps(jcp)%cp)

            ! copy input for wtd to gd, except for cps that have gd%subs specified separately

            if (.not.cp%has_own_subs) call subsurf_copy(wtd%subs, cp%gd%subs)

            ! calculate subsurface stresses, using mirroring for left rail/wheel combination

            cp%gd%meta%ncase = cp%gd%meta%ncase - 1
            call subsur(cp%gd%meta, cp%gd%ic, cp%gd%mater, cp%gd%cgrid_cur, cp%gd%outpt1%igs,           &
                           cp%gd%outpt1%ps, cp%gd%ic%is_left_side(), cp%gd%subs)
            cp%gd%meta%ncase = cp%gd%meta%ncase + 1

            end associate
         enddo
      else

         ! module 3: plain call to subsur_calc, gd%subs already filled in

         gd%meta%ncase = gd%meta%ncase - 1
         call subsur(gd%meta, gd%ic, gd%mater, gd%cgrid_cur, gd%outpt1%igs, gd%outpt1%ps, .false.,      &
                        gd%subs)
         gd%meta%ncase = gd%meta%ncase + 1
      endif
      call timer_stop(itimer)

   endif ! ierror==0

   ! release lock for parallel computation

   call free_contact_problem(ire, icp)

   if (idebug.ge.2) then
      write(bufout,'(2a)')       trim(pfx_str(subnam,ire,icp)),' done...'
      call write_log(1, bufout)
   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine subs_calculate

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getFlags(ire, icp, lenflg, params, values) &
   bind(c,name=CNAME_(cntc_getflags))
!--function: used for getting various configuration flags from a contact problem
!--category: 7, "m=any, wtd or cp":  available for modules 1 and 3, in module 1 working on single or all cps
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire            ! result element ID
   integer,      intent(in)  :: icp            ! contact problem ID
   integer,      intent(in)  :: lenflg         ! length of params/values arrays
   integer,      intent(in)  :: params(lenflg) ! codes of the parameters to be obtained from CONTACT
   integer,      intent(out) :: values(lenflg) ! values of the parameters obtained from CONTACT
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_getFlags'
   integer  :: i, imodul, ierror
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getFlags
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 0, subnam, ierror)
   if (ierror.lt.0) return

   imodul = ire_module(ix_reid(ire))

   do i = 1, lenflg

      ! write(bufout,'(2(a,i6),2a)') ' cntc_getflags: values(',i,') = param', params(i),': ',   &
      !         cntc_flagname(params(i))
      ! call write_log(1, bufout)

      if (params(i).eq.0) then

         ! do nothing, allow arrays params/values to contain unused positions

      elseif (params(i).eq.CNTC_ic_config) then         ! get C1-digit from the RE-CP-data

         values(i) = my_ic%config

      elseif (params(i).eq.CNTC_ic_pvtime) then         ! get P-digit from the RE-CP-data

         values(i) = my_ic%pvtime

      elseif (params(i).eq.CNTC_ic_bound) then          ! get B-digit from the RE-CP-data

         values(i) = my_ic%bound

      elseif (params(i).eq.CNTC_ic_tang) then           ! get T-digit from the RE-CP-data

         values(i) = my_ic%tang

      elseif (params(i).eq.CNTC_ic_norm) then           ! get N-digit from the RE-CP-data

         values(i) = my_ic%norm

      elseif (params(i).eq.CNTC_ic_force) then          ! get F-digit from the RE-CP-data

         if (icp.le.0) then
            values(i) = my_ic%force1
         else
            values(i) = my_ic%force3
         endif

      elseif (params(i).eq.CNTC_ic_sens) then           ! get S2-digit from the RE-CP-data

         values(i) = my_ic%sens

      elseif (params(i).eq.CNTC_ic_frclaw) then         ! get L-digit from the RE-CP-data

         values(i) = my_ic%frclaw_inp

      elseif (params(i).eq.CNTC_ic_discns) then         ! get D-digit from the RE-CP-data

         values(i) = my_ic%discns1_inp

      elseif (params(i).eq.CNTC_ic_npomax) then         ! get max. #elements in potential contact (mod. 1)

         values(i) = wtd%discr%npot_max

      elseif (params(i).eq.CNTC_ic_inflcf) then         ! get C3-digit from the RE-CP-data

         values(i) = my_ic%gencr_inp

      elseif (params(i).eq.CNTC_ic_ifmeth) then         ! get IF_METH for Blanco-approach (C3=4)

         values(i) = my_mater%if_meth

      elseif (params(i).eq.CNTC_ic_ifvari) then         ! get IF_VER for Blanco-approach (C3=4)

         values(i) = my_mater%if_ver

      elseif (params(i).eq.CNTC_ic_mater) then          ! get M-digit from the RE-CP-data

         values(i) = my_ic%mater

      elseif (params(i).eq.CNTC_ic_xflow) then          ! get X-digit from the RE-CP-data

         ! note: usage is inconsistent with setflags; we could return psflcin

         values(i) = my_ic%xflow

      elseif (params(i).eq.CNTC_ic_heat) then           ! get H-digit from the RE-CP-data

         values(i) = my_ic%heat

      elseif (params(i).eq.CNTC_ic_iestim) then         ! get I-digit from the RE-CP-data

         values(i) = my_ic%iestim

      elseif (params(i).eq.CNTC_ic_matfil) then         ! get A-digit from the RE-CP-data

         values(i) = my_ic%matfil_surf

      elseif (params(i).eq.CNTC_ic_output) then         ! get O-digit from the RE-CP-data

         values(i) = my_ic%output_surf

      elseif (params(i).eq.CNTC_ic_sbsfil) then         ! get A_s-digit from the RE-CP-data

         values(i) = my_ic%matfil_subs

      elseif (params(i).eq.CNTC_ic_sbsout) then         ! get O_s-digit from the RE-CP-data

         values(i) = my_ic%output_subs

      elseif (params(i).eq.CNTC_ic_flow) then           ! get W-digit from the RE-CP-data

         values(i) = my_ic%flow

      elseif (params(i).eq.CNTC_ic_return) then         ! get R-digit from the RE-CP-data

         values(i) = my_ic%return

      elseif (params(i).eq.CNTC_if_wrtinp) then         ! get flag wrtinp from the RE-CP-data

         values(i) = my_ic%wrtinp

      else

         write(bufout,'(a,i8,a)') ' WARNING: unknown parameter code=',params(i),' is ignored.'
         call write_log(1, bufout)

      endif
   enddo

   if (idebug.ge.2) then
      do i = 1, lenflg
         if (params(i).ne.0) then
            write(bufout,'(2a,i2,3a,i8)') trim(pfx_str(subnam,ire,icp)), ' i=',i,', got flag ',         &
                        cntc_flagName(params(i)),', value=',values(i)
            call write_log(1, bufout)
         endif
      enddo
   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getFlags

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getParameters(ire, icp, itask, lenarr, values) &
   bind(c,name=CNAME_(cntc_getparameters))
!--function: internal routine for retrieving various parameters from a contact problem
!            itask: selected group of parameters (1: cntc_getcpresults, 2: material, 3: friction)
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire            ! result element ID
   integer,      intent(in)  :: icp            ! contact problem ID
   integer,      intent(in)  :: itask          ! task: selected group of parameters
   integer,      intent(in)  :: lenarr         ! length of values array
   real(kind=8), intent(out) :: values(lenarr) ! values of the parameters obtained from CONTACT
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_getParameters'
   integer  :: imodul, ivf, iof, ierror
   type(t_friclaw), pointer :: fric
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getParameters
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)

   if (ierror.lt.0) return

   if (itask.eq.1) then

      ! 1: parameters used by plot3d

   if (lenarr.ge.1) values(1) = gd%kin%veloc    / my_scl%veloc
   if (lenarr.ge.2) values(2) = gd%kin%chi      / my_scl%angle
   if (lenarr.ge.3) values(3) = gd%kin%dq       / my_scl%len
   if (lenarr.ge.4) values(4) = gd%kin%spinxo   / my_scl%len
   if (lenarr.ge.5) values(5) = gd%kin%spinyo   / my_scl%len
   if (lenarr.ge.6) values(6) = gd%mater%tau_c0 * my_scl%area

   elseif (itask.eq.2) then

      ! 2: material parameters

      if (lenarr.ge. 1) values( 1) = gd%mater%gg(1) * my_scl%area
      if (lenarr.ge. 2) values( 2) = gd%mater%gg(2) * my_scl%area
      if (lenarr.ge. 3) values( 3) = gd%mater%ga    * my_scl%area
      if (lenarr.ge. 4) values( 4) = gd%mater%poiss(1)
      if (lenarr.ge. 5) values( 5) = gd%mater%poiss(2)
      if (lenarr.ge. 6) values( 6) = gd%mater%nu
      if (lenarr.ge. 7) values( 7) = gd%mater%ak
      if (lenarr.ge. 8) values( 8) = gd%mater%flx(1) / my_scl%len**3
      if (lenarr.ge. 9) values( 9) = gd%mater%flx(2) / my_scl%len**3
      if (lenarr.ge.10) values(10) = gd%mater%flx(3) / my_scl%len**3
      if (lenarr.ge.11) values(11) = gd%mater%k0_mf
      if (lenarr.ge.12) values(12) = gd%mater%alfamf
      if (lenarr.ge.13) values(13) = gd%mater%betamf
      if (lenarr.ge.14) values(14) = gd%mater%k_eff
      if (lenarr.ge.15) values(15) = gd%mater%gg3    * my_scl%area
      if (lenarr.ge.16) values(16) = gd%mater%laythk / my_scl%len
      if (lenarr.ge.17) values(17) = gd%mater%tau_c0 * my_scl%area
      if (lenarr.ge.18) values(18) = gd%mater%k_tau  * my_scl%area * my_scl%len
      if (lenarr.ge.19) values(19) = gd%mater%cdampn
      if (lenarr.ge.20) values(20) = gd%mater%cdampt
      if (lenarr.ge.21) values(21) = gd%mater%dfnmax
      if (lenarr.ge.22) values(22) = gd%mater%dftmax
        ! visco-elastic: fg, tc, vt, akv, gav, nuv
        ! temperature: bktemp, heatcp, lambda, dens, betapl
        ! elastic layer: flx_z
        ! numerical infl.coef: if_meth, if_ver

   elseif (itask.eq.3) then

      ! 3: friction parameters

      imodul = ire_module(ix_reid(ire))
      if (imodul.eq.1) then
         fric => wtd%fric
      else
         fric => gd%fric
      endif

      if (lenarr.ge.1) values(1) = fric%frclaw_eff
      if (lenarr.ge.2) values(2) = max(1, fric%nvf)
      if (lenarr.ge.3) values(3) = fric%memdst
      if (lenarr.ge.4) values(4) = fric%mem_s0
      do ivf = 1, max(1, fric%nvf)
         iof = 4 + 3*(ivf-1)
         if (lenarr.ge.iof+1) values(iof+1) = fric%paramvf  (ivf)
         if (lenarr.ge.iof+2) values(iof+2) = fric%fstat_arr(ivf)
         if (lenarr.ge.iof+3) values(iof+3) = fric%fkin_arr (ivf)
         ! slip-velocity dependent friction: flin1/2, frat1/2, fexp1/2, sabsh1/2
         ! temperature dependent friction: fref, tref, dfheat, dheat?
      enddo

   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getParameters

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getProfileValues(ire, itask, nints, iparam, nreals, rparam, lenarr, val) &
   bind(c,name=CNAME_(cntc_getprofilevalues))
!--function: return a wheel or rail profile for a wheel-rail contact problem as a table of values
!  itask          - select type of outputs:
!                    -1: ierr   <0: error codes, 0=profile loaded ok, 1=not set
!                     0: npnt   number of points used in requested sampling method
!                     1: r/w    get (yr,zr) values for rail or (yw,zw) for wheel profile
!                     2: trk    get (ytr,ztr) values for rail or wheel profile (principal profile)
!                     3: gaug   get left-most point within gauge height, offset, point at gauge height
!                               [ygauge1, zgauge1, yoffs, zoffs, ygauge2, zgauge2]  [length]
!                     4: arc    get arc-length parameter s along profile
!                     5: angl   get surface inclination atan2(dz, dy) [angle]
!                     6: 2d-u   evaluate varprof [xyz]_r or [xyz]_w for cross-section at constant u_i = x_i
!                     7: 2d-v   evaluate varprof [xyz]_r or [xyz]_w for interpolation path at constant v_j
!                     8: 2d-y   evaluate varprof [xyz]_r or [xyz]_w for longitudinal slice at constant y_j
!                     9: spl2d  dump 2D spline data-structure for inspection in Matlab
!  iparam         - integer configuration parameters
!                     1: itype     0 = rail, 1 = wheel profile
!                     2: isampl   -1 = sampling cf. original input data;
!                                  0 = sampling cf. spline representation (default for tasks 0--5);
!                                  1 = sampling cf. spline representation at spacing ds_out (required for
!                                      tasks 6--8)
!                            kchk>=2 = sampling cf. spline representation with integer refinement factor
!  rparam         - real configuration parameters
!                     1: ds_out  step-size ds/du/dv used with sampling method isampl=1, default 1mm
!                     2: c1_out  sample position on variable profiles (task 1, 2: x_out)
!                     3: c2_sta  start of sample range on variable profiles
!                     4: c2_end  end of sample range on variable profiles
!                   task 6: c1 = u, c2 = v; task 7: c1 = v, c2 = u; task 8: c1 = y, c2 = u.
!  units: s,x,y,z [mm], th [rad], u,v [-]
!  tasks 1--9: no unit conversion or scaling are applied for profile values
!--category: 2, "m=1 only, wtd":    available for module 1 only, working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: itask         ! integer code for requested values
   integer,      intent(in)  :: nints, nreals ! number of integer/real parameters provided
   integer,      intent(in)  :: iparam(*)     ! array of integer parameters
   real(kind=8), intent(in)  :: rparam(*)     ! array of real parameters
   integer,      intent(in)  :: lenarr        ! length of output array
   real(kind=8), intent(out) :: val(lenarr)   ! output values for selected task
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_getProfileValues'
   character(len=5)            :: namtyp(0:1) = (/ 'rail ', 'wheel' /)
   integer                     :: itype, isampl, npnt, ip, j, jp, iu, ldebug, ierror, sub_ierror
   logical                     :: is_ok
   real(kind=8)                :: s0, s1, ds, ds_out, x_out, sgn, ygauge, yvampr, zgauge, zmin,  &
                                  c2_sta, c2_end, uarr(1), varr(1), yarr(1)
   type(t_marker)              :: rw_trk
   type(t_grid),     target    :: g_wrk, g_slc
   type(t_gridfnc3)            :: dxyz
   type(t_profile),  pointer   :: my_prf
   type(t_grid),     pointer   :: g
   character(len=256)          :: fulnam
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getProfileValues
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, -1, 1, -1, subnam, ierror)
   if (ierror.lt.0) return

   itype   = 0  ! default: rail
   if (nints.ge.1) itype   = max(0, min(1, iparam(1)))     ! 0=rail, 1=wheel

   ! select rail or wheel profile data-structure

   associate(my_rail  => wtd%trk%rai, my_wheel => wtd%ws%whl)
   if (itype.eq.0) then
      my_prf   => my_rail%prr
   else
      my_prf   => my_wheel%prw
   endif

   ! check that the profile was loaded ok

   if (itask.ne.-1 .and. my_prf%ierror.ne.0) then
      write(bufout,'(a,2(a,i0))') trim(pfx_str(subnam,ire,-1)),' ERROR: task=',itask,                   &
                ': profile could not be loaded, ierror=',my_prf%ierror
      call write_log(1, bufout)
      return
   endif

   ! tasks 6--9: check that a variable profile is given

   if (itask.ge.6 .and. .not.my_prf%is_varprof()) then
      write(bufout,'(2a,2(i0,a))') trim(pfx_str(subnam,ire,-1)),' ERROR: task=',itask,', type=',itype,  &
             ': needs a variable profile'
      call write_log(1, bufout)
      return
   endif

   ! unpack option values, fill in defaults

   isampl  = 0  ! default: use s/v-positions from profile spline
   if (nints.ge.2) isampl  = max(-1, min(100, iparam(2)))
   if (itask.ge.6 .and. itask.le.8) isampl = 1          ! 2d spline: using uniform sampling

   ds_out  = 1d0
   if (nreals.ge.1) ds_out = rparam(1)

   x_out   = 0d0
   if (nreals.ge.2) x_out  = rparam(2)

   if (itask.eq.6) then                                 ! cross-section constant u
      c2_sta  = my_prf%spl2d%tvj(4)
      c2_end  = my_prf%spl2d%tvj(my_prf%spl2d%nknotv-3)
   elseif (itask.eq.7 .or. itask.eq.8) then             ! cross-section constant v or y
      c2_sta  = my_prf%spl2d%tui(4)
      c2_end  = my_prf%spl2d%tui(my_prf%spl2d%nknotu-3)
   endif
   if (nreals.ge.3) c2_sta = rparam(3)
   if (nreals.ge.4) c2_end = rparam(4)

   if (idebug.ge.2) then
      write(bufout,'(a,3(a,i2),2(a,g12.4))') trim(pfx_str(subnam,ire,-1)),' task=',itask,               &
                ', itype=',itype,', isampl=',isampl,', ds_out=',ds_out,', x_out=',x_out
      call write_log(1, bufout)
   endif

   ! reject negative ds_out; reject very small ds_out (huge npnt)

   if (isampl.eq.1 .and. ds_out.lt.1d-4 .and. itask.ge.0) then
      isampl = 0
      write(bufout,'(2a,g12.4,a)') trim(pfx_str(subnam,ire,-1)),' ds_out=',ds_out,                      &
                ' too small, using default sampling'
      call write_log(1, bufout)
   endif

   ! tasks 0--5: interpolate variable profile to requested x_out, using profile g_slc instead of varprof

   if (itask.le.5 .and. my_prf%ierror.eq.0 .and. my_prf%is_varprof()) then

      ! task 2: interpret x_out as track coordinate, add s_ws to get rail coordinate
      if (itask.eq.2 .and. itype.eq.0) x_out = x_out + wtd%ws%s

      if (idebug.ge.3 .and. itype.eq.0) then
         write(bufout,'(a,g12.4)') ' setting up "current slice" for variable rail at x_out=',x_out
         call write_log(1, bufout)
      elseif (idebug.ge.3) then
         write(bufout,'(a,g12.4)') ' setting up "current slice" for out-of-round wheel at th_out=',x_out
         call write_log(1, bufout)
      endif

      call varprof_intpol_xunif(my_prf, 0d0, 1, x_out, 0d0, g_slc, sub_ierror)
      ierror = sub_ierror
      g      => g_slc

   else

      ! set pointer g to actual data

      g     => my_prf%grd_data
      ierror = my_prf%ierror

   endif

   if (idebug.ge.3) then
      write(bufout,'(4a)') trim(pfx_str(subnam,ire,-1)),' fname="', trim(my_prf%fname),'"'
      call write_log(1, bufout)
   endif

   ! task -1: return error code

   is_ok = (ierror.eq.0)

   if (itask.eq.-1) then
      if (lenarr.ge.1 .and.      is_ok) val(1) = 0
      if (lenarr.ge.1 .and. .not.is_ok) val(1) = CNTC_err_profil
      if (lenarr.ge.2) val(2) = ierror
      return
   endif

   ! check availability of profile

   if (.not.is_ok .or. .not.associated(g%coor) .or. (itask.eq.4 .and. .not.associated(g%s_prf))) then
      write(bufout,'(3a)') trim(pfx_str(subnam,ire,-1)),' no data found for ',trim(namtyp(itype))
      call write_log(1, bufout)

      if (lenarr.ge.1) val(1) = CNTC_err_profil
      return
   endif

   ! mirror track y-values for left side

   sgn = 1d0
   if (my_ic%is_left_side()) sgn = -1d0

   ! determine number of output s-positions

   if (isampl.eq.-1) then       ! -1: s-positions of input data

      npnt = g%ntot

   elseif (isampl.eq.0) then    !  0: s-positions of spline representation

      npnt = g%spl%npnt

   elseif (isampl.eq.1) then    !  1: uniform sampling at spacing ds_out

      if (itask.ge.6 .and. itask.le.8) then
         npnt = int( (c2_end - c2_sta) / ds_out ) + 1
      else
         s0   = g%spl%s(1)
         s1   = g%spl%s(g%spl%npnt)
         npnt = int( (s1 - s0) / ds_out ) + 1
      endif
      
   elseif (isampl.ge.2 .and. isampl.le.100) then ! using 'isampl' points in each spline interval 

      npnt = isampl * (g%spl%npnt - 1) + 1

   else

      write(bufout,'(2a,i9)') trim(pfx_str(subnam,ire,-1)),' invalid isampl =',isampl
      call write_log(1, bufout)
      return

   endif
      
   if (idebug.ge.3) then
      write(bufout,'(2a,i6)') trim(pfx_str(subnam,ire,-1)),' npnt=',npnt
      call write_log(1, bufout)
   endif

   ! allocate working grid

   call grid_create_curvil(g_wrk, 1, npnt)
   call reallocate_arr(g_wrk%s_prf, npnt)

   ! fill output s-positions

   if (isampl.eq.-1) then       ! -1: s-positions of input data

      g_wrk%s_prf(1:npnt) = g%s_prf(1:npnt)

   elseif (isampl.eq.0) then    !  0: s-positions of spline representation

      g_wrk%s_prf(1:npnt) = g%spl%s(1:npnt)

   elseif (isampl.eq.1) then    !  1: uniform sampling at spacing ds_out

      if (itask.ge.6 .and. itask.le.8) then
         s0 = c2_sta
         s1 = c2_end
      else
         s0   = g%spl%s(1)
         s1   = g%spl%s(g%spl%npnt)
      endif

      do ip = 1, npnt
         g_wrk%s_prf(ip) = s0 + (ip-1) * ds_out
      enddo
      
   elseif (isampl.ge.2 .and. isampl.le.100) then ! using 'isampl' points in each spline interval 

      g_wrk%s_prf(1) = g%spl%s(1)
      do ip = 1, g%spl%npnt-1
         s0   = g%spl%s(ip)
         s1   = g%spl%s(ip+1)
         ds   = (s1 - s0) / (1d0*isampl)
         do j = 1, isampl
            jp = 1 + isampl * (ip-1) + j
            g_wrk%s_prf(jp) = s0 + j * ds
         enddo
      enddo

   endif
      
   ! evaluate spline (y,z) at requested s-positions

   if (itask.eq.1 .or. itask.eq.2) then
      call spline_eval(g%spl, ikYDIR, npnt, g_wrk%s_prf, ierror, 999d0, g_wrk%y)
      call spline_eval(g%spl, ikZDIR, npnt, g_wrk%s_prf, ierror, 999d0, g_wrk%z)
   elseif (itask.eq.6) then                     ! lateral cross-section at constant u
      uarr(1) = x_out
      call bspline_eval2d_prod(my_prf%spl2d, 1, npnt, uarr, g_wrk%s_prf, g_wrk%x, g_wrk%y, g_wrk%z,     &
                .false., sub_ierror, 99d9)
      do iu = 1, npnt
         g_wrk%x(iu) = uarr(1)
      enddo
   elseif (itask.eq.7) then                     ! longitudinal interpolation path at constant v
      varr(1) = x_out
      call bspline_eval2d_prod(my_prf%spl2d, npnt, 1, g_wrk%s_prf, varr, g_wrk%x, g_wrk%y, g_wrk%z,     &
                .false., sub_ierror, 999d0)
      do iu = 1, npnt
         g_wrk%x(iu) = g_wrk%s_prf(iu)
      enddo
   elseif (itask.eq.8) then                     ! longitudinal slice at constant y
      yarr(1) = x_out
      do iu = 1, npnt
         g_wrk%x(iu) = g_wrk%s_prf(iu)
         g_wrk%y(iu) = yarr(1)
      enddo
      if (idebug.ge.3) then
         write(bufout,'(a,f9.3,a,i0,3(a,g12.3),a)') ' evaluate spl2d at yj=',yarr(1),', npnt=',npnt,    &
                   ', ui=[',g_wrk%s_prf(1),':',g_wrk%s_prf(2)-g_wrk%s_prf(1),':',g_wrk%s_prf(npnt),']'
         call write_log(1, bufout)
      endif
      call bspline_get_z_at_xy_prod(my_prf%spl2d, npnt, 1, g_wrk%s_prf, yarr, g_wrk%z, sub_ierror, 999d0)
   endif

   ! copy output-values

   if     (itask.eq.0) then     ! 0 = npnt

      ! output number of points in profile

      if (lenarr.ge.1) val(1) = npnt

   elseif (itask.eq.1) then     ! 1 = r/w

      ! output rail or wheel profile in its own coordinates

      do ip = 1, min(lenarr     ,npnt)
         val(     ip) = g_wrk%y(ip)
      enddo
      do ip = 1, min(lenarr-npnt,npnt)
         val(npnt+ip) = g_wrk%z(ip)
      enddo

   elseif (itask.eq.2) then     ! 2 = trk

      ! convert the profile to track coordinates

      if (itype.eq.0) then      ! rail
         call cartgrid_2glob(g_wrk, my_rail%m_trk)
      else                      ! wheel
         rw_trk = marker_2glob( my_wheel%m_ws, wtd%ws%m_trk )
         ! call marker_print(my_wheel%m_ws, 'm_ws', 2)
         ! call marker_print(wtd%ws%m_trk, 'm_trk', 2)
         ! call marker_print(rw_trk, 'rw_trk', 2)
         call cartgrid_2glob(g_wrk, rw_trk)
      endif

      ! output rail or wheel profile in track coordinates

      do ip = 1, min(lenarr     ,npnt)
         val(     ip) = sgn * g_wrk%y(ip)
      enddo
      do ip = 1, min(lenarr-npnt,npnt)
         val(npnt+ip) =       g_wrk%z(ip)
      enddo

   elseif (itask.eq.3) then     ! 3 = gaug

      if (itype.eq.1) then

         ! for wheels no gauge point is available

         call write_log('ERROR: gauge data not available for wheel profiles')
         if (lenarr.ge.1) val(1) = CNTC_err_profil
      else

         ldebug = 0
         call find_gauge_meas_pt(g, wtd%trk%cant_angle, wtd%trk%gauge_height, ygauge, zgauge, yvampr,   &
                                 zmin, ldebug)

         if (lenarr.ge.1) val(1) = ygauge                                                  / my_scl%len
         if (lenarr.ge.2) val(2) = zgauge                                                  / my_scl%len
         if (lenarr.ge.3) val(3) = (sgn * (wtd%trk%track_gauge/2d0 - ygauge) + my_rail%dy) / my_scl%len
         if (lenarr.ge.4) val(4) = (my_rail%dz - zmin)                                     / my_scl%len
         if (lenarr.ge.5) val(5) = yvampr                                                  / my_scl%len
         if (lenarr.ge.6) val(6) = (wtd%trk%gauge_height + zmin)                           / my_scl%len

      endif

   elseif (itask.eq.4) then     ! 4 = arc

      ! output rail or wheel arc-length parameter

      do ip = 1, min(lenarr, npnt)
         val(ip) = g_wrk%s_prf(ip)
      enddo

   elseif (itask.eq.5) then     ! 5 = angl

      ! output surface inclination on rail or wheel

      call gf3_new(dxyz, 'dxyz', g_wrk, lzero=.true.)
      call spline_get_dxyz_at_s_gfout(g%spl, dxyz, sub_ierror)
      ! call gf3_print(dxyz, 'dxyz', ikALL, 4)
      do ip = 1, min(lenarr     ,npnt)
         val(ip) = atan2( dxyz%vn(ip), dxyz%vy(ip) ) / my_scl%angle
      enddo
      call gf3_destroy(dxyz)

   elseif (itask.ge.6 .and. itask.le.8) then     ! 6-8 = 2d-u, 2d-v, 2d-y

      if (lenarr.eq.1) then

         ! output number of points required

         val(1) = npnt

      else

         ! output rail or wheel (x,y,z) positions

         do ip = 1, npnt
            if (       ip.le.lenarr) val(       ip) = g_wrk%x(ip)
            if (  npnt+ip.le.lenarr) val(  npnt+ip) = g_wrk%y(ip)
            if (2*npnt+ip.le.lenarr) val(2*npnt+ip) = g_wrk%z(ip)
         enddo
      endif

   elseif (itask.eq.9) then     ! 9 = spl2d

      call make_absolute_path('dump_spl2d.m', wtd%meta%outdir, fulnam)

      if (itype.eq.0 .and. my_rail%prr%is_varprof()) then

         call bspline2d_dump_matlab(my_rail%prr%spl2d, itype, fulnam, 0)

      elseif (itype.eq.1 .and. my_wheel%prw%is_varprof()) then

         call bspline2d_dump_matlab(my_wheel%prw%spl2d, itype, fulnam, 0)

      else
         write(bufout,'(2a,2(i0,a))') trim(pfx_str(subnam,ire,-1)),' task=',itask,', type=',itype,      &
                ': not a variable profile'
         call write_log(1, bufout)
      endif

   endif

   call grid_destroy(g_wrk)
   if (my_prf%is_varprof()) call grid_destroy(g_slc)

   end associate
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getProfileValues

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getProfileValues_new(ire, itask, nints, iparam, nreals, rparam, lenarr, val) &
   bind(c,name=CNAME_(cntc_getprofilevalues_new))
!--function: return a wheel or rail profile for a wheel-rail contact problem as a table of values
!            see cntc_getProfileValues for documentation
!--category: 2, "m=1 only, wtd":    available for module 1 only, working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: itask         ! integer code for requested values
   integer,      intent(in)  :: nints, nreals ! number of integer/real parameters provided
   integer,      intent(in)  :: iparam(*)     ! array of integer parameters
   real(kind=8), intent(in)  :: rparam(*)     ! array of real parameters
   integer,      intent(in)  :: lenarr        ! length of output array
   real(kind=8), intent(out) :: val(lenarr)   ! output values for selected task
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getProfileValues_new
#endif

   call write_log(' getprofilevalues_new: calling getprofilevalues...')
   call cntc_getprofilevalues(ire, itask, nints, iparam, nreals, rparam, lenarr, val)

end subroutine cntc_getProfileValues_new

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getWheelsetPosition(ire, lenarr, rvalues) &
   bind(c,name=CNAME_(cntc_getwheelsetposition))
!--function: get the wheelset position data for a wheel-rail contact problem
!
!  dimensions:   s_ws, y_ws, z_ws [length],       roll, yaw, pitch [angle]
!--category: 2, "m=1 only, wtd":    available for module 1 only, working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: lenarr        ! length of output array
   real(kind=8), intent(out) :: rvalues(lenarr) ! wheelset position for current result element
!--local variables:
   integer             :: ierror
   character(len=*), parameter :: subnam = 'cntc_getWheelsetPosition'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getWheelsetPosition
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, -1, 1, -1, subnam, ierror)
   if (ierror.lt.0) return

   if (lenarr.ge.1) rvalues(1) = wtd%ws%s     / my_scl%len
   if (lenarr.ge.2) rvalues(2) = wtd%ws%y     / my_scl%len
   if (lenarr.ge.3) rvalues(3) = wtd%ws%z     / my_scl%len
   if (lenarr.ge.4) rvalues(4) = wtd%ws%roll  / my_scl%angle
   if (lenarr.ge.5) rvalues(5) = wtd%ws%yaw   / my_scl%angle
   if (lenarr.ge.6) rvalues(6) = wtd%ws%pitch / my_scl%angle

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getWheelsetPosition

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getWheelsetVelocity(ire, lenarr, rvalues) &
   bind(c,name=CNAME_(cntc_getwheelsetvelocity))
!--function: get the wheelset velocity data for a wheel-rail contact problem
!
!  dimensions:  vs_ws, vy_ws, vz_ws [veloc],     vroll, vyaw, vpitch [ang.veloc]
!--category: 2, "m=1 only, wtd":    available for module 1 only, working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: lenarr        ! length of output array
   real(kind=8), intent(out) :: rvalues(lenarr) ! wheelset velocity for current result element
!--local variables:
   integer             :: ierror
   character(len=*), parameter :: subnam = 'cntc_getWheelsetVelocity'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getWheelsetVelocity
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, -1, 1, -1, subnam, ierror)
   if (ierror.lt.0) return

   if (lenarr.ge.1) rvalues(1) = wtd%ws%vs     / my_scl%veloc
   if (lenarr.ge.2) rvalues(2) = wtd%ws%vy     / my_scl%veloc
   if (lenarr.ge.3) rvalues(3) = wtd%ws%vz     / my_scl%veloc
   if (lenarr.ge.4) rvalues(4) = wtd%ws%vroll  / my_scl%angle ! angle/time
   if (lenarr.ge.5) rvalues(5) = wtd%ws%vyaw   / my_scl%angle ! angle/time
   if (lenarr.ge.6) rvalues(6) = wtd%ws%vpitch / my_scl%angle ! angle/time

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getWheelsetVelocity

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getNumContactPatches(ire, npatch) &
   bind(c,name=CNAME_(cntc_getnumcontactpatches))
!--function: return the number of contact patches used in a w/r contact problem
!--category: 2, "m=1 only, wtd":    available for module 1 only, working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(out) :: npatch        ! number of separate contact patches
!--local variables:
   integer    :: ierror
   character(len=*), parameter :: subnam = 'cntc_getNumContactPatches'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getNumContactPatches
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, -1, 1, -1, subnam, ierror)
   if (ierror.lt.0) return

   npatch = wtd%numcps

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getNumContactPatches

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getContactLocation(ire, icp, lenarr, rvalues) &
   bind(c,name=CNAME_(cntc_getcontactlocation))
!--function: return the contact reference location for a contact patch in a w/r contact problem
!  The following values are returned, if permitted by the length of rvalues:
!   1 - XCP_TR    - x-position of the contact reference point in track coordinates
!   2 - YCP_TR    - y-position of the contact reference point in track coordinates
!   3 - ZCP_TR    - z-position of the contact reference point in track coordinates
!   4 - DELTCP_TR - contact reference angle: rotation about track x-axis from the track positive z-axis to
!                   the contact positive n-axis, with sign according the right-hand rule
!
!   5 - XCP_R     - x-position of the contact reference point in rail profile coordinates
!   6 - YCP_R     - y-position of the contact reference point in rail profile coordinates
!   7 - ZCP_R     - z-position of the contact reference point in rail profile coordinates
!   8 - SCP_R     - s-parameter of the contact reference point, measured along the rail profile
!   9 - DELTCP_R  - rotation about rail x-axis from rail positive z-axis to contact positive n-axis
!
!  10 - XCP_W     - x-position of the contact reference point in wheel profile coordinates
!  11 - YCP_W     - y-position of the contact reference point in wheel profile coordinates
!  12 - ZCP_W     - z-position of the contact reference point in wheel profile coordinates
!  13 - SCP_W     - s-parameter of the contact reference point, measured along the wheel profile
!  14 - DELTCP_W  - rotation about wheel x-axis from wheel positive z-axis to contact positive n-axis
!
!  15 - XPN_TR    - x-position of the pressure center of gravity in track coordinates
!  16 - YPN_TR    - y-position of the pressure center of gravity in track coordinates
!
!  18 - DY_DEFL   - lateral rail shift according to massless rail deflection
!  19 - DZ_DEFL   - vertical rail shift according to massless rail deflection
!
!  21 - XW_TR     - x-position of wheel profile marker in track coordinates
!  22 - YW_TR     - y-position of wheel profile marker in track coordinates
!  23 - ZR_TR     - z-position of wheel profile marker in track coordinates
!  24 - ROLLW_TR  - roll angle of wheel profile marker in track coordinates
!  25 - YAWW_TR   - yaw angle of wheel profile marker in track coordinates
!
!  26 - YR_TR     - y-position of rail profile marker in track coordinates
!  27 - ZR_TR     - z-position of rail profile marker in track coordinates
!  28 - ROLLR_TR  - roll angle of rail profile marker in track coordinates
!
!  30 - XCP_WS    - x-position of the contact reference point in wheelset coordinates
!  31 - YCP_WS    - y-position of the contact reference point in wheelset coordinates
!  32 - ZCP_WS    - z-position of the contact reference point in wheelset coordinates
!  33 - DELTCP_WS - rotation about wheelset x-axis from wheelset positive z-axis to contact positive n-axis
!
!  34 - XW_WS     - x-position of wheel profile marker in wheelset coordinates
!  35 - YW_WS     - y-position of wheel profile marker in wheelset coordinates
!  36 - ZR_WS     - z-position of wheel profile marker in wheelset coordinates
!  37 - ROLLW_WS  - roll angle of wheel profile marker in wheelset coordinates
!  38 - YAWW_WS   - yaw angle of wheel profile marker in wheelset coordinates
!
!  The "contact reference point" is the origin of the contact local coordinate system. It is determined by
!  a heuristic rule and is centered within the contact patch in a weighted sense.
!--category: 3, "m=1 only, cp":     available for module 1 only, working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   integer,      intent(in)  :: lenarr        ! length of output array
   real(kind=8), intent(out) :: rvalues(lenarr) ! contact location for current contact patch
!--local variables:
   integer                     :: ierror, ii
   real(kind=8)                :: sgn, sum_pn, sum_xpn, sum_ypn
   type(t_marker)              :: m_pn, mref_ws
   character(len=*), parameter :: subnam = 'cntc_getContactLocation'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getContactLocation
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 1, 1, subnam, ierror)
   if (ierror.lt.0) return

   if (icp.le.0 .or. icp.gt.wtd%numcps) then

      ! contact patch number icp out of range: return zero values

      rvalues(1:lenarr) = 0d0
   else

      ! set reference to the active wheel in the current configuration

      associate(my_wheel => wtd%ws%whl, cp   => wtd%allcps(icp)%cp, meta => wtd%allcps(icp)%cp%gd%meta)

      ! set the sign to -1 for left and +1 for right rail/wheel combination

      sgn = 1d0
      if (my_ic%is_left_side()) sgn = -1d0

      !  1-- 4: contact reference marker in track frame
      if (lenarr.ge.1) rvalues(1)   =       meta%xcp_tr    / my_scl%len
      if (lenarr.ge.2) rvalues(2)   = sgn * meta%ycp_tr    / my_scl%len
      if (lenarr.ge.3) rvalues(3)   =       meta%zcp_tr    / my_scl%len
      if (lenarr.ge.4) rvalues(4)   = sgn * meta%deltcp_tr / my_scl%angle

      !  5-- 9: contact reference marker in rail frame
      if (lenarr.ge.5) rvalues(5)   =       meta%xcp_r     / my_scl%len
      if (lenarr.ge.6) rvalues(6)   = sgn * meta%ycp_r     / my_scl%len
      if (lenarr.ge.7) rvalues(7)   =       meta%zcp_r     / my_scl%len
      if (lenarr.ge.8) rvalues(8)   =       meta%scp_r     / my_scl%len
      if (lenarr.ge.9) rvalues(9)   = sgn * meta%deltcp_r  / my_scl%angle

      ! 10--14: contact reference marker in wheel frame
      if (lenarr.ge.10) rvalues(10) =       meta%xcp_w     / my_scl%len
      if (lenarr.ge.11) rvalues(11) = sgn * meta%ycp_w     / my_scl%len
      if (lenarr.ge.12) rvalues(12) =       meta%zcp_w     / my_scl%len
      if (lenarr.ge.13) rvalues(13) =       meta%scp_w     / my_scl%len
      if (lenarr.ge.14) rvalues(14) = sgn * meta%deltcp_w  / my_scl%angle

      if (lenarr.ge.15) then

         ! compute the pressure center of gravity in x and y
         !      - weighted mean \int_C  x * pn / \int_C pn; \int_C y * pn / \int_C pn.

         sum_xpn = 0d0
         sum_ypn = 0d0
         sum_pn  = 1d-10

         do ii = 1, cp%gd%cgrid_cur%ntot
            sum_xpn = sum_xpn + cp%gd%cgrid_cur%x(ii) * cp%gd%outpt1%ps%vn(ii)
            sum_ypn = sum_ypn + cp%gd%cgrid_cur%y(ii) * cp%gd%outpt1%ps%vn(ii)
            sum_pn  = sum_pn  + cp%gd%outpt1%ps%vn(ii)
         enddo

         ! create marker for pressure c.g. in contact-coordinates

         m_pn = marker( sum_xpn/sum_pn, sum_ypn/sum_pn, 0d0 )
         m_pn = marker_2glob( m_pn, cp%mref )
      endif

      ! 15--16: pressure center of gravity
      if (lenarr.ge.15) rvalues(15) =       m_pn%x() / my_scl%len
      if (lenarr.ge.16) rvalues(16) = sgn * m_pn%y() / my_scl%len
      if (lenarr.ge.17) rvalues(17) = 0d0

      ! 18--19: massless rail deflections
      if (lenarr.ge.18) rvalues(18) =       wtd%trk%dy_defl / my_scl%len
      if (lenarr.ge.19) rvalues(19) =       wtd%trk%dz_defl / my_scl%len
      if (lenarr.ge.20) rvalues(20) = 0d0

      ! 21--25: wheel marker in track frame
      if (lenarr.ge.21) rvalues(21) =       meta%x_w        / my_scl%len
      if (lenarr.ge.22) rvalues(22) = sgn * meta%y_w        / my_scl%len
      if (lenarr.ge.23) rvalues(23) =       meta%z_w        / my_scl%len
      if (lenarr.ge.24) rvalues(24) = sgn * meta%roll_w     / my_scl%angle
      if (lenarr.ge.25) rvalues(25) = sgn * meta%yaw_w      / my_scl%angle

      ! 26--28: rail marker in track frame
      if (lenarr.ge.26) rvalues(26) = sgn * meta%y_r        / my_scl%len
      if (lenarr.ge.27) rvalues(27) =       meta%z_r        / my_scl%len
      if (lenarr.ge.28) rvalues(28) = sgn * meta%roll_r     / my_scl%angle

      ! 30--33: contact reference marker in wheelset frame
      mref_ws  = marker_2loc( cp%mref, wtd%ws%m_trk )
      if (lenarr.ge.30) rvalues(30) =       mref_ws%x()    / my_scl%len
      if (lenarr.ge.31) rvalues(31) = sgn * mref_ws%y()    / my_scl%len
      if (lenarr.ge.32) rvalues(32) =       mref_ws%z()    / my_scl%len
      if (lenarr.ge.33) rvalues(33) = sgn * mref_ws%roll() / my_scl%angle

      ! 34--38: wheel marker in wheelset frame
      if (lenarr.ge.34) rvalues(34) =       my_wheel%m_ws%x()    / my_scl%len
      if (lenarr.ge.35) rvalues(35) = sgn * my_wheel%m_ws%y()    / my_scl%len
      if (lenarr.ge.36) rvalues(36) =       my_wheel%m_ws%z()    / my_scl%len
      if (lenarr.ge.37) rvalues(37) = sgn * my_wheel%m_ws%roll() / my_scl%angle
      if (lenarr.ge.38) rvalues(38) = sgn * my_wheel%m_ws%yaw()  / my_scl%angle
      end associate
   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getContactLocation

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getReferenceVelocity(ire, icp, veloc) &
   bind(c,name=CNAME_(cntc_getreferencevelocity))
!--function: get the rolling velocity of a contact problem
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   real(kind=8), intent(out) :: veloc         ! absolute rolling velocity [veloc]
!--local variables:
   integer    :: ierror
   character(len=*), parameter :: subnam = 'cntc_getReferenceVelocity'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getReferenceVelocity
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)
   if (ierror.lt.0) return

   ! note: retrieving data from gd%kin instead of my_kin (=> wtd%kin)

   veloc = gd%kin%veloc / my_scl%veloc

   if (idebug.ge.2) then
      write(bufout,'(2a,f9.1,a)') trim(pfx_str(subnam,ire,icp)),' V=',gd%kin%veloc,' [mm/s]'
      call write_log(1, bufout)
   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getReferenceVelocity

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getHertzContact(ire, icp, lenarr, rvalues) &
   bind(c,name=CNAME_(cntc_gethertzcontact))
!--function: get the parameters from a Hertzian contact problem
!
! The following values are returned, if permitted by the length of rvalues:
!   1 - A1       - curvature in rolling direction [1/length]
!   2 - B1       - curvature in lateral direction [1/length]
!   3 - AA       - semi-axis in rolling direction [length]
!   4 - BB       - semi-axis in lateral direction [length]
!   5 - RHO      - effective radius of curvature, 2 / (A1 + B1) [length]
!   6 - CP       - effective semi-axis, sqrt(AA * BB) [length]
!   7 - SCALE    - potential contact scale factor [-]
!   8 - BNEG     - semi-axis of negative half-ellipse in lateral direction [length]
!   9 - BPOS     - semi-axis of positive half-ellipse in lateral direction [length]
!  10 - AOB      - ellipticity AA/BB [-]
!--category: 1, "m=3 only, cp":     available for module 3 only, working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire             ! result element ID
   integer,      intent(in)  :: icp             ! contact problem ID
   integer,      intent(in)  :: lenarr          ! length of array provided
   real(kind=8), intent(out) :: rvalues(lenarr) ! parameters of Hertzian solution
!--local variables:
   integer      :: ierror
   character(len=*), parameter :: subnam = 'cntc_getHertzContact'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getHertzContact
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 3, 1, subnam, ierror)
   if (ierror.lt.0) return

   if (lenarr.ge. 1) rvalues( 1) = gd%hertz%a1   * my_scl%len
   if (lenarr.ge. 2) rvalues( 2) = gd%hertz%b1   * my_scl%len
   if (lenarr.ge. 3) rvalues( 3) = gd%hertz%aa   / my_scl%len
   if (lenarr.ge. 4) rvalues( 4) = gd%hertz%bb   / my_scl%len
   if (lenarr.ge. 5) rvalues( 5) = gd%hertz%rho  / my_scl%len
   if (lenarr.ge. 6) rvalues( 6) = gd%hertz%cp   / my_scl%len
   if (lenarr.ge. 7) rvalues( 7) = gd%hertz%scale
   if (lenarr.ge. 8) rvalues( 8) = gd%hertz%bneg / my_scl%len
   if (lenarr.ge. 9) rvalues( 9) = gd%hertz%bpos / my_scl%len
   if (lenarr.ge.10) rvalues(10) = gd%hertz%aob

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getHertzContact

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getNumElements(ire, icp, mx, my) &
   bind(c,name=CNAME_(cntc_getnumelements))
!--function: return the number of elements in the potential contact area used for a contact problem
!            length of tractions arrays
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   integer,      intent(out) :: mx, my        ! number of discretization elements in long/lat dirs
!--local variables:
   integer    :: ierror
   character(len=*), parameter :: subnam = 'cntc_getNumElements'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getNumElements
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)
   if (ierror.lt.0) return

   mx = gd%potcon_cur%mx
   my = gd%potcon_cur%my

   if (idebug.ge.2) then
      write(bufout,'(2a,2i4)')   trim(pfx_str(subnam,ire,icp)),' #elements mx,my=', mx,my
      call write_log(1, bufout)
   endif
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getNumElements

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getGridDiscretization(ire, icp, dx, dy) &
   bind(c,name=CNAME_(cntc_getgriddiscretization))
!--function: get the grid discretization step sizes dx,dy for a contact problem
!--category: 7, "m=any, wtd or cp":  available for modules 1 and 3, in module 1 working on single or all cps
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   real(kind=8), intent(out) :: dx, dy        ! grid discretization step sizes [length]
!--local variables:
   integer    :: ierror
   character(len=*), parameter :: subnam = 'cntc_getGridDiscretization'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getGridDiscretization
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 0, subnam, ierror)
   if (ierror.lt.0) return

   if (icp.eq.-1) then

      dx = wtd%discr%dx / my_scl%len
      dy = wtd%discr%ds / my_scl%len

   else

      dx = gd%potcon_cur%dx / my_scl%len
      dy = gd%potcon_cur%dy / my_scl%len

   endif

   if (idebug.ge.2) then
      write(bufout,'(2a,2f8.5,a)') trim(pfx_str(subnam,ire,icp)),' step sizes dx,dy=',dx, dy,' [length]'
      call write_log(1, bufout)
   endif
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getGridDiscretization

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getPotContact(ire, icp, lenarr, rvalues) &
   bind(c,name=CNAME_(cntc_getpotcontact))
!--function: return the parameters of the potential contact area for a contact problem, according to ipotcn=3
!    3: 1st center + grid sizes,    params = [ mx, my, xc1, yc1, dx, dy ]
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire             ! result element ID
   integer,      intent(in)  :: icp             ! contact problem ID
   integer,      intent(in)  :: lenarr          ! length of output array
   real(kind=8), intent(out) :: rvalues(lenarr) ! num.elements [-], positions and step sizes [length]
!--local variables:
   integer    :: imodul, ierror
   character(len=*), parameter :: subnam = 'cntc_getPotContact'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getPotContact
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)
   if (ierror.lt.0) return

   ! use mirroring for left rail/wheel combination

   imodul = ire_module(ix_reid(ire))

   if (lenarr.ge.1) rvalues(1) = real(gd%potcon_cur%mx)
   if (lenarr.ge.2) rvalues(2) = real(gd%potcon_cur%my)
   if (lenarr.ge.3) rvalues(3) = gd%potcon_cur%xc1 / my_scl%len
   if (lenarr.ge.4 .and. imodul.eq.1 .and. my_ic%is_left_side()) then
      rvalues(4) = -gd%potcon_cur%ycm / my_scl%len  ! renumber [my:-1:1] and mirror y:=-y
   elseif (lenarr.ge.4) then
      rvalues(4) =  gd%potcon_cur%yc1 / my_scl%len
   endif
   if (lenarr.ge.5) rvalues(5) = gd%potcon_cur%dx / my_scl%len
   if (lenarr.ge.6) rvalues(6) = gd%potcon_cur%dy / my_scl%len
   if (lenarr.ge.7) rvalues(7:lenarr) = 0d0

   !if (idebug.ge.2) then
   !   write(bufout,'(2a,2f8.5,a)') trim(pfx_str(subnam,ire,icp)),' step sizes dx,dy=',dx, dy,' [length]'
   !   call write_log(1, bufout)
   !endif
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getPotContact

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getPenetration(ire, icp, pen) &
   bind(c,name=CNAME_(cntc_getpenetration))
!--function: return the penetration for a contact problem
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   real(kind=8), intent(out) :: pen           ! penetration [length]
!--local variables:
   integer    :: ierror
   character(len=*), parameter :: subnam = 'cntc_getPenetration'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getPenetration
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)
   if (ierror.lt.0) return

   ! note: retrieving data from gd%kin instead of my_kin (=> wtd%kin)

   pen = gd%kin%pen / my_scl%len

   if (idebug.ge.2) then
      write(bufout,'(2a,f9.6,a)') trim(pfx_str(subnam,ire,icp)),' penetration pen=',pen,' [length]'
      call write_log(1, bufout)
   endif
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getPenetration

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getCreepages(ire, icp, vx, vy, phi) &
   bind(c,name=CNAME_(cntc_getcreepages))
!--function: return the kinematic constants (creepages) for a contact problem
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire            ! result element ID
   integer,      intent(in)  :: icp            ! contact problem ID
   real(kind=8), intent(out) :: vx, vy, phi    ! long/lat/spin creepages [-], [angle/length]
!--local variables:
   integer      :: imodul, ierror
   real(kind=8) :: sgn
   character(len=*), parameter :: subnam = 'cntc_getCreepages'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getCreepages
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)
   if (ierror.lt.0) return

   ! set the sign to -1 for left and +1 for right rail/wheel combination

   sgn    = 1d0
   imodul = ire_module(ix_reid(ire))
   if (imodul.eq.1 .and. my_ic%is_left_side()) sgn = -1d0

   ! note: retrieving data from gd%kin instead of my_kin (=> wtd%kin)

   if (my_ic%tang.eq.1) then
      ! retrieving the rigid shift [mm], [rad]          ! in simpack, output body == body 2
      vx  =       my_scl%body * gd%kin%cksi / my_scl%len
      vy  = sgn * my_scl%body * gd%kin%ceta / my_scl%len
      phi = sgn * my_scl%body * gd%kin%cphi
   else
      ! retrieving the creepages [-], [rad/mm]
      vx  =       my_scl%body * gd%kin%cksi
      vy  = sgn * my_scl%body * gd%kin%ceta
      phi = sgn * my_scl%body * gd%kin%cphi * my_scl%len
   endif

   if (idebug.ge.2) then
      write(bufout,'(2a,2f7.4,a,f7.4,a)') trim(pfx_str(subnam,ire,icp)),' creepages Cksi,Ceta=',        &
                gd%kin%cksi, gd%kin%ceta, ' [-], Phi=',gd%kin%cphi, ' [rad/mm]'
      call write_log(1, bufout)
   endif
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getCreepages

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getContactForces(ire, icp, fn, tx, ty, mz) &
   bind(c,name=CNAME_(cntc_getcontactforces))
!--function: return the total forces and torsional moment for a contact problem wrt contact reference
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   real(kind=8), intent(out) :: fn            ! total normal force [force]
   real(kind=8), intent(out) :: tx, ty        ! total long/lat forces [force]
   real(kind=8), intent(out) :: mz            ! total torsional moment [force.length]
!--local variables:
   integer      :: imodul, ierror
   real(kind=8) :: sgn
   character(len=*), parameter :: subnam = 'cntc_getContactForces'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getContactForces
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)
   if (ierror.lt.0) return

   ! set the sign to -1 for left and +1 for right rail/wheel combination

   sgn    = 1d0
   imodul = ire_module(ix_reid(ire))
   if (imodul.eq.1 .and. my_ic%is_left_side()) sgn = -1d0

   fn =                     my_kin%fcntc(3)     ! note: in simpack, output body == body 2
   tx =       my_scl%body * my_kin%fcntc(1)
   ty = sgn * my_scl%body * my_kin%fcntc(2)
   mz = sgn * my_scl%body * gd%outpt1%mztrue / my_scl%len

   if (idebug.ge.2) then
      write(bufout,'(2a,2f9.1,a)') trim(pfx_str(subnam,ire,icp)),' total forces fn,mz=',fn,mz,          &
                ' [force], [force.length]'
      call write_log(1, bufout)
      write(bufout,'(2a,2f9.1,a)') trim(pfx_str(subnam,ire,icp)),' total forces tx,ty=',tx,ty,' [force]'
      call write_log(1, bufout)
   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getContactForces

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getGlobalForces(ire, icp, lenarr, rvalues) &
   bind(c,name=CNAME_(cntc_getglobalforces))
!--function: return the overall forces for a w/r contact problem (module 1 only)
!
! rvalues(lenarr)   - contact forces [force] and moments [force.length]
!
! The following values are returned, if permitted by the length of rvalues:
!   1 - FX_TR    - total force on the output body, component in track longitudinal x-direction
!   2 - FY_TR    - total force on the output body, component in track lateral y-direction
!   3 - FZ_TR    - total force on the output body, component in track vertical z-direction
!   4 - MX_R_TR  - total moment on output body about rail profile marker, component in track x-direction
!   5 - MY_R_TR  - total moment on output body about rail profile marker, component in track y-direction
!   6 - MZ_R_TR  - total moment on output body about rail profile marker, component in track z-direction
!
!   7 - FX_WS    - total force on the output body, component in wheelset longitudinal x-direction
!   8 - FY_WS    - total force on the output body, component in wheelset lateral y-direction
!   9 - FZ_WS    - total force on the output body, component in wheelset vertical z-direction
!  10 - MX_W_WS  - total moment on output body about wheel profile marker, component in wheelset x-direction
!  11 - MY_W_WS  - total moment on output body about wheel profile marker, component in wheelset y-direction
!  12 - MZ_W_WS  - total moment on output body about wheel profile marker, component in wheelset z-direction
!                  Note that the 'output body' is the rail when using CONTACTs unit convention
!
!  13 - FX_R     - total force on the output body, component in rail profile x-direction
!  14 - FY_R     - total force on the output body, component in rail profile y-direction
!  15 - FZ_R     - total force on the output body, component in rail profile z-direction
!  16 - MX_R_R   - total moment on output body about rail profile marker, component in rail x-direction
!  17 - MY_R_R   - total moment on output body about rail profile marker, component in rail y-direction
!  18 - MZ_R_R   - total moment on output body about rail profile marker, component in rail z-direction
!
!  19 - FX_W     - total force on the output body, component in wheel profile x-direction
!  20 - FY_W     - total force on the output body, component in wheel profile y-direction
!  21 - FZ_W     - total force on the output body, component in wheel profile z-direction
!  22 - MX_W_W   - total moment on output body about wheel profile marker, component in wheel x-direction
!  23 - MY_W_W   - total moment on output body about wheel profile marker, component in wheel y-direction
!  24 - MZ_W_W   - total moment on output body about wheel profile marker, component in wheel z-direction
!
!--category: 4, "m=1 only, wtd/cp": available for module 1 only, working on wtd or cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   integer,      intent(in)  :: lenarr        ! length of output array
   real(kind=8), intent(out) :: rvalues(lenarr) ! forces [force] and moments [force.length]
!--local variables:
   integer                     :: ierror
   real(kind=8)                :: sgn
   type(t_vec),      pointer   :: ftot_tr, ttot_r_tr, ftot_ws, ttot_w_ws
   type(t_vec)                 :: ftot_r, ttot_r_r, ftot_w, ttot_w_w
   type(t_vec),      target    :: fzero
   character(len=*), parameter :: subnam = 'cntc_getGlobalForces'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getGlobalForces
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 1, 0, subnam, ierror)
   if (ierror.lt.0) return

   ! set reference to the active wheel in the current configuration

   associate(my_rail  => wtd%trk%rai, my_wheel => wtd%ws%whl)
   fzero   = vec( 0d0, 0d0, 0d0 )

   ! set the sign to -1 for left and +1 for right rail/wheel combination

   sgn = 1d0
   if (my_ic%is_left_side()) sgn = -1d0

   ! total forces and moments (on rail) provided w.r.t. track/wheelset x,y,z-coordinates

   if (icp.eq.-1) then

      ftot_tr     => wtd%ftrk
      ftot_ws     => wtd%fws
      ttot_r_tr   => wtd%ttrk
      ttot_w_ws   => wtd%tws

   elseif (icp.ge.1 .and. icp.le.wtd%numcps) then

      ftot_tr     => wtd%allcps(icp)%cp%ftrk
      ftot_ws     => wtd%allcps(icp)%cp%fws
      ttot_r_tr   => wtd%allcps(icp)%cp%ttrk
      ttot_w_ws   => wtd%allcps(icp)%cp%tws

   else

      ftot_tr     => fzero
      ftot_ws     => fzero
      ttot_r_tr   => fzero
      ttot_w_ws   => fzero

   endif

   ! compute additional values in rail/wheel x,y,z-coordinates

   ftot_r     = my_rail%m_trk .transp. ftot_tr
   ttot_r_r   = my_rail%m_trk .transp. ttot_r_tr
   ftot_w     = my_wheel%m_ws .transp. ftot_ws
   ttot_w_w   = my_wheel%m_ws .transp. ttot_w_ws

   ! output of values

   if (lenarr.ge. 1) rvalues( 1) =       my_scl%body * ftot_tr%x()                 ! FX(TR)
   if (lenarr.ge. 2) rvalues( 2) = sgn * my_scl%body * ftot_tr%y()                 ! FY(TR)
   if (lenarr.ge. 3) rvalues( 3) =       my_scl%body * ftot_tr%z()                 ! FZ(TR)
   if (lenarr.ge. 4) rvalues( 4) = sgn * my_scl%body * ttot_r_tr%x() / my_scl%len  ! MX_@R(TR)
   if (lenarr.ge. 5) rvalues( 5) =       my_scl%body * ttot_r_tr%y() / my_scl%len  ! MY_@R(TR)
   if (lenarr.ge. 6) rvalues( 6) = sgn * my_scl%body * ttot_r_tr%z() / my_scl%len  ! MZ_@R(TR)

   if (lenarr.ge. 7) rvalues( 7) =       my_scl%body * ftot_ws%x()                 ! FX(WS)
   if (lenarr.ge. 8) rvalues( 8) = sgn * my_scl%body * ftot_ws%y()                 ! FY(WS)
   if (lenarr.ge. 9) rvalues( 9) =       my_scl%body * ftot_ws%z()                 ! FZ(WS)
   if (lenarr.ge.10) rvalues(10) = sgn * my_scl%body * ttot_w_ws%x() / my_scl%len  ! MX_@W(WS)
   if (lenarr.ge.11) rvalues(11) =       my_scl%body * ttot_w_ws%y() / my_scl%len  ! MY_@W(WS)
   if (lenarr.ge.12) rvalues(12) = sgn * my_scl%body * ttot_w_ws%z() / my_scl%len  ! MZ_@W(WS)

   if (lenarr.ge.13) rvalues(13) =       my_scl%body * ftot_r%x()                  ! FX(RR)
   if (lenarr.ge.14) rvalues(14) = sgn * my_scl%body * ftot_r%y()                  ! FY(RR)
   if (lenarr.ge.15) rvalues(15) =       my_scl%body * ftot_r%z()                  ! FZ(RR)
   if (lenarr.ge.16) rvalues(16) = sgn * my_scl%body * ttot_r_r%x() / my_scl%len   ! MX_@R(RR)
   if (lenarr.ge.17) rvalues(17) =       my_scl%body * ttot_r_r%y() / my_scl%len   ! MY_@R(RR)
   if (lenarr.ge.18) rvalues(18) = sgn * my_scl%body * ttot_r_r%z() / my_scl%len   ! MZ_@R(RR)

   if (lenarr.ge.19) rvalues(19) =       my_scl%body * ftot_w%x()                  ! FX(RW)
   if (lenarr.ge.20) rvalues(20) = sgn * my_scl%body * ftot_w%y()                  ! FY(RW)
   if (lenarr.ge.21) rvalues(21) =       my_scl%body * ftot_w%z()                  ! FZ(RW)
   if (lenarr.ge.22) rvalues(22) = sgn * my_scl%body * ttot_w_w%x() / my_scl%len   ! MX_@W(RW)
   if (lenarr.ge.23) rvalues(23) =       my_scl%body * ttot_w_w%y() / my_scl%len   ! MY_@W(RW)
   if (lenarr.ge.24) rvalues(24) = sgn * my_scl%body * ttot_w_w%z() / my_scl%len   ! MZ_@W(RW)

   if (lenarr.ge.25) rvalues(25:) = 0d0

   end associate
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getGlobalForces

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getContactPatchAreas(ire, icp, carea, harea, sarea) &
   bind(c,name=CNAME_(cntc_getcontactpatchareas))
!--function: return the area of contact problem
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   real(kind=8), intent(out) :: carea         ! area of contact patch [area]
   real(kind=8), intent(out) :: harea         ! area of adhesion area [area]
   real(kind=8), intent(out) :: sarea         ! area of slip area [area]
!  the size of the plasticity area follows as carea - harea - sarea
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_getContactPatchAreas'
   integer  :: ncon, nadh, nslip, nplast, nexter, ierror
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getContactPatchAreas
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)
   if (ierror.lt.0) return

   ! Count number of elements in contact, slip and adhesion

   call eldiv_count(gd%outpt1%igs, nadh, nslip, nplast, nexter)
   ncon = nadh + nslip + nplast
   carea = real(ncon)  * gd%potcon_cur%dxdy / my_scl%area
   harea = real(nadh)  * gd%potcon_cur%dxdy / my_scl%area
   sarea = real(nslip) * gd%potcon_cur%dxdy / my_scl%area

   if (idebug.ge.2) then
      write(bufout,'(2a,3(f9.6,a))') trim(pfx_str(subnam,ire,icp)),' carea=',carea,', harea=',harea,    &
                ', sarea=',sarea,' [area]'
      call write_log(1, bufout)
   endif
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getContactPatchAreas

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getElementDivision(ire, icp, lenarr, eldiv) &
   bind(c,name=CNAME_(cntc_getelementdivision))
!--function: return flags for all elements in the potential contact area for a contact problem
!            indicating whether the element is in Exerior (0), Adhesion (1), Slip (2) or Plasticity (3).
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   integer,      intent(in)  :: lenarr        ! length of output array
   integer,      intent(out) :: eldiv(lenarr) ! element division of contact area, 0=E, 1=H, 2=S, 3=P.
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_getElementDivision'
   integer    :: imodul, ix, iy, ii, jy, jj, ierror
   logical    :: mirror_y
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getElementDivision
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)
   if (ierror.lt.0) return

   ! use mirroring (renumbering) for left rail/wheel combination

   imodul   = ire_module(ix_reid(ire))
   mirror_y = (imodul.eq.1 .and. (my_ic%config.eq.0 .or. my_ic%config.eq.4))

   if (idebug.ge.3) then
      write(bufout,'(2a)') trim(pfx_str(subnam,ire,icp)),' return element division array'
      call write_log(1, bufout)
   endif

   if (.not.associated(gd%outpt1%igs%el)) then
      write(bufout,'(2a)') trim(pfx_str(subnam,ire,icp)),' element division not yet available'
      call write_log(1, bufout)
   else
      do jy = 1, gd%cgrid_cur%ny
         do ix = 1, gd%cgrid_cur%nx
            if (mirror_y) then
               iy = gd%cgrid_cur%ny + 1 - jy
            else
               iy = jy
            endif
            jj = ix + gd%cgrid_cur%nx * (jy-1)    ! external numbering jj == (ix,jy)
            ii = ix + gd%cgrid_cur%nx * (iy-1)    ! internal numbering ii == (ix,iy)
            if (jj.le.lenarr) eldiv(jj) = gd%outpt1%igs%el(ii)
         enddo
      enddo
   endif
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getElementDivision

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getMaximumPressure(ire, icp, pnmax) &
   bind(c,name=CNAME_(cntc_getmaximumpressure))
!--function: return the maximum normal pressure in contact problem
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   real(kind=8), intent(out) :: pnmax         ! maximum pressure in contact patch [force/area]
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_getMaximumPressure'
   integer    :: ii, iimax, ierror
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getMaximumPressure
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)
   if (ierror.lt.0) return

   iimax = 1
   do ii = 1, gd%cgrid_cur%ntot
      if (abs(gd%outpt1%ps%vn(ii)).gt.abs(gd%outpt1%ps%vn(iimax))) iimax = ii
   enddo
   pnmax = abs(gd%outpt1%ps%vn(iimax)) * my_scl%area

   if (idebug.ge.2) then
      write(bufout,'(2a,g10.3,a)') trim(pfx_str(subnam,ire,icp)),' maximum pressure=',pnmax,' [force/area]'
      call write_log(1, bufout)
   endif
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getMaximumPressure

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getMaximumTraction(ire, icp, ptmax) &
   bind(c,name=CNAME_(cntc_getmaximumtraction))
!--function: return the maximum tangential traction in contact problem
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   real(kind=8), intent(out) :: ptmax         ! maximum traction |pt| in contact patch [force/area]
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_getMaximumTraction'
   integer      :: ii, ierror
   real(kind=8) :: ptabs
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getMaximumTraction
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)
   if (ierror.lt.0) return

   ptmax = -1d0
   do ii = 1, gd%cgrid_cur%ntot
      ptabs = sqrt( gd%outpt1%ps%vx(ii)**2 + gd%outpt1%ps%vy(ii)**2 )
      if (ptabs.gt.ptmax) then
         ptmax = ptabs
      endif
   enddo
   ptmax = ptmax * my_scl%area

   if (idebug.ge.2) then
      write(bufout,'(2a,g10.3,a)') trim(pfx_str(subnam,ire,icp)),' maximum traction=',ptmax,' [force/area]'
      call write_log(1, bufout)
   endif
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getMaximumTraction

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getMaximumTemperature(ire, icp, t1max, t2max) &
   bind(c,name=CNAME_(cntc_getmaximumtemperature))
!--function: return the maximum temperatures in bodies 1 and 2 in a contact problem
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   real(kind=8), intent(out) :: t1max, t2max  ! maximum surface temperatures in contact patch [C]
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_getMaximumTemperature'
   integer    :: ierror
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getMaximumTemperature
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)
   if (ierror.lt.0) return

   if (my_ic%heat.eq.0 .or. .not.associated(gd%outpt1%temp1%vn)) then
      t1max = 0d0
      t2max = 0d0
      write(bufout,'(2a)') trim(pfx_str(subnam,ire,icp)),' temperature not available'
      call write_log(1, bufout)
      return
   endif

   t1max = gf3_max(AllElm, gd%outpt1%temp1, ikZDIR)
   t2max = gf3_max(AllElm, gd%outpt1%temp2, ikZDIR)

   if (idebug.ge.2) then
      write(bufout,'(2a,2(g10.3,a))') trim(pfx_str(subnam,ire,icp)),' maximum temperatures=',t1max,',', &
                t2max,' [C]'
      call write_log(1, bufout)
   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getMaximumTemperature

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getFieldData(ire, icp, ifld, lenarr, fld) &
   bind(c,name=CNAME_(cntc_getfielddata))
!--function: return an array with values for all elements in the potential contact for a contact problem
!          fields provided: CNTC_fld_h      =    1 ! for retrieving array h      [length]
!                           CNTC_fld_mu     =    2 ! for retrieving array mu     [-]
!                           CNTC_fld_px     =    3 ! for retrieving array px     [force/area]
!                           CNTC_fld_py     =    4 ! for retrieving array py     [force/area]
!                           CNTC_fld_pn     =    5 ! for retrieving array pn     [force/area]
!                           CNTC_fld_ux     =    7 ! for retrieving array ux     [length]
!                           CNTC_fld_uy     =    8 ! for retrieving array uy     [length]
!                           CNTC_fld_un     =    9 ! for retrieving array un     [length]
!                           CNTC_fld_taucrt =   11 ! for retrieving array taucrt [force/area]
!                           CNTC_fld_uplsx  =   12 ! for retrieving array uplsx  [length]
!                           CNTC_fld_uplsy  =   13 ! for retrieving array uplsy  [length]
!                           CNTC_fld_sx     =   15 ! for retrieving array sx     [-]
!                           CNTC_fld_sy     =   16 ! for retrieving array sy     [-]
!                           CNTC_fld_temp1  =   20 ! for retrieving array temp1  [C]
!                           CNTC_fld_temp2  =   21 ! for retrieving array temp2  [C]
!                           CNTC_fld_wx     =   22 ! for retrieving array wx     [-]
!                           CNTC_fld_wy     =   23 ! for retrieving array wy     [-]
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   integer,      intent(in)  :: ifld          ! integer code for requested field
   integer,      intent(in)  :: lenarr        ! length of output arrays
   real(kind=8), intent(out) :: fld(lenarr)   ! output values for all elements of contact area
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_getFieldData'
   integer           :: imodul, ix, iy, ii, jy, jj, ierror
   logical           :: mirror_y
   real(kind=8)      :: sgn, sclfac
   character(len=6)  :: fldnam
   real(kind=8), dimension(:), pointer :: val
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getFieldData
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)
   if (ierror.lt.0) return

   ! set mirroring for left rail/wheel combination

   imodul   = ire_module(ix_reid(ire))
   mirror_y = (imodul.eq.1 .and. (my_ic%config.eq.0 .or. my_ic%config.eq.4))
   sgn      = 1d0
   if (mirror_y) sgn = -1d0

   ! look up the requested field

   if (ifld.eq.CNTC_fld_h) then
      fldnam = 'h'
      sclfac =       1d0 / my_scl%len
      val    => gd%geom%hs1%vn
   elseif (ifld.eq.CNTC_fld_mu) then
      fldnam = 'mu'
      sclfac =       1d0
      val    => gd%outpt1%mus%vt
   elseif (ifld.eq.CNTC_fld_px) then
      fldnam = 'px'
      sclfac =       my_scl%area * my_scl%body
      val    => gd%outpt1%ps%vx
   elseif (ifld.eq.CNTC_fld_py) then
      fldnam = 'py'
      sclfac = sgn * my_scl%area * my_scl%body
      val    => gd%outpt1%ps%vy
   elseif (ifld.eq.CNTC_fld_pn) then
      fldnam = 'pn'
      sclfac =       my_scl%area
      val    => gd%outpt1%ps%vn
   elseif (ifld.eq.CNTC_fld_ux) then
      fldnam = 'ux'
      sclfac =       1d0 / my_scl%len
      val    => gd%outpt1%us%vx
   elseif (ifld.eq.CNTC_fld_uy) then
      fldnam = 'uy'
      sclfac = sgn * 1d0 / my_scl%len
      val    => gd%outpt1%us%vy
   elseif (ifld.eq.CNTC_fld_un) then
      fldnam = 'un'
      sclfac =       1d0 / my_scl%len
      val    => gd%outpt1%us%vn
   elseif (ifld.eq.CNTC_fld_taucrt) then
      fldnam = 'taucrt'
      sclfac =       my_scl%area
      val    => gd%outpt1%taucs%vt
   elseif (ifld.eq.CNTC_fld_uplsx) then
      fldnam = 'uplsx'
      sclfac =       1d0 / my_scl%len
      val    => gd%outpt1%upls%vx
   elseif (ifld.eq.CNTC_fld_uplsy) then
      fldnam = 'uplsy'
      sclfac = sgn * 1d0 / my_scl%len
      val    => gd%outpt1%upls%vy
   elseif (ifld.eq.CNTC_fld_sx) then
      fldnam = 'sx'
      sclfac =       my_scl%body / my_kin%dq
      val    => gd%outpt1%ss%vx
   elseif (ifld.eq.CNTC_fld_sy) then
      fldnam = 'sy'
      sclfac = sgn * my_scl%body / my_kin%dq
      val    => gd%outpt1%ss%vy
   elseif (ifld.eq.CNTC_fld_wx) then
      fldnam = 'wx'     ! ex. cksi when F>0
      sclfac =       my_scl%body / my_kin%dq
      val    => gd%geom%hs1%vx
   elseif (ifld.eq.CNTC_fld_wy) then
      fldnam = 'wy'     ! ex. ceta when F=2
      sclfac = sgn * my_scl%body / my_kin%dq
      val    => gd%geom%hs1%vy
   elseif (ifld.eq.CNTC_fld_temp1) then
      fldnam = 'temp1'
      sclfac =       1d0
      val    => gd%outpt1%temp1%vn
   elseif (ifld.eq.CNTC_fld_temp2) then
      fldnam = 'temp2'
      sclfac =       1d0
      val    => gd%outpt1%temp2%vn
   else
      write(bufout,'(2a,i4,a)') trim(pfx_str(subnam,ire,icp)),' field',ifld,' does not exist.'
      call write_log(1, bufout)
      return
   endif

   if (idebug.ge.3) then
      write(bufout,'(4a,g12.4)') trim(pfx_str(subnam,ire,icp)),' return field-array ',fldnam,         &
                ', sclfac=', sclfac
      call write_log(1, bufout)
   endif

   if (.not.associated(val)) then
      write(bufout,'(4a)') trim(pfx_str(subnam,ire,icp)),' array ', fldnam,' not available'
      call write_log(1, bufout)
   else
      do jy = 1, gd%cgrid_cur%ny
         do ix = 1, gd%cgrid_cur%nx
            if (mirror_y) then
               iy = gd%cgrid_cur%ny + 1 - jy
            else
               iy = jy
            endif
            jj = ix + gd%cgrid_cur%nx * (jy-1)    ! external numbering jj == (ix,jy)
            ii = ix + gd%cgrid_cur%nx * (iy-1)    ! internal numbering ii == (ix,iy)
            if (jj.le.lenarr) fld(jj) = sclfac * val(ii)
         enddo
      enddo
   endif

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getFieldData

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getTractions(ire, icp, lenarr, pn, px, py) &
   bind(c,name=CNAME_(cntc_gettractions))
!--function: return the tractions for all elements in the potential contact area for a contact problem
!            note the order of the arguments, with pn (z-direction) occurring before px,py.
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   integer,      intent(in)  :: lenarr        ! length of output arrays
   real(kind=8), intent(out) :: pn(lenarr), px(lenarr), py(lenarr)    ! surface tractions for all elements
                                                                      ! of contact area [force/area]
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_getTractions'
   integer :: ierror
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getTractions
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)
   if (ierror.lt.0) return

   call cntc_getFieldData(ire, icp, CNTC_fld_px, lenarr, px)
   call cntc_getFieldData(ire, icp, CNTC_fld_py, lenarr, py)
   call cntc_getFieldData(ire, icp, CNTC_fld_pn, lenarr, pn)

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getTractions

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getMicroSlip(ire, icp, lenarr, sx, sy) &
   bind(c,name=CNAME_(cntc_getmicroslip))
!--function: return the relative micro-slip velocity for all elements in the potential contact area for
!            contact problem
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   integer,      intent(in)  :: lenarr        ! length of output arrays
   real(kind=8), intent(out) :: sx(lenarr), sy(lenarr) ! relative micro-slip velocity [-]
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_getMicroSlip'
   integer :: ierror
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getMicroSlip
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)
   if (ierror.lt.0) return

   call cntc_getFieldData(ire, icp, CNTC_fld_sx, lenarr, sx)
   call cntc_getFieldData(ire, icp, CNTC_fld_sy, lenarr, sy)

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getMicroSlip

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getDisplacements(ire, icp, lenarr, un, ux, uy) &
   bind(c,name=CNAME_(cntc_getdisplacements))
!--function: return the displ.differences for all elements in the potential contact area for a contact problem
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   integer,      intent(in)  :: lenarr        ! length of output arrays
   real(kind=8), intent(out) :: un(lenarr), ux(lenarr), uy(lenarr) ! displacement difference [length]
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_getDisplacements'
   integer :: ierror
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getDisplacements
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)
   if (ierror.lt.0) return

   call cntc_getFieldData(ire, icp, CNTC_fld_ux, lenarr, ux)
   call cntc_getFieldData(ire, icp, CNTC_fld_uy, lenarr, uy)
   call cntc_getFieldData(ire, icp, CNTC_fld_un, lenarr, un)

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getDisplacements

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getSensitivities(ire, icp, lenout, lenin, sens) &
   bind(c,name=CNAME_(cntc_getsensitivities))
!--function: return the sensitivities of the total forces for a contact problem
!
!  lenout, lenin      - requested number of outputs (forces) and inputs (creepages or shifts)
!  sens(lenout,lenin) - matrix of sensitivities -- zeros if not computed
!
!  calculation of sensitivities is requested using the flag 'CNTC_ic_sens'.
!  accuracy is configured using cntc_setsolverflags with G=6.
!
!  the inputs are ordered   1: pen, 2: cksi, 3: ceta, 4: cphi
!  in rolling, T=2,3, the units are pen [length], cksi, ceta [-],      cphi [angle/length]
!  in shifts,  T=1,   the units are pen [length], cksi, ceta [length], cphi [angle]
!
!  the outputs are ordered  1: fn,  2: fx,   3: fy,   4: mz
!  the units are fn [force], fx, fy [-], mz [force.length]
!
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire                ! result element ID
   integer,      intent(in)  :: icp                ! contact problem ID
   integer,      intent(in)  :: lenout, lenin      ! requested number of output/input variables
   real(kind=8), intent(out) :: sens(lenout*lenin) ! matrix of sensitivities
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_getSensitivities'
   integer                    :: imodul, iout, iin, ii, ierror
   logical                    :: mirror_y
   real(kind=8)               :: sgn
   real(kind=8), dimension(:) :: fac_in(nsens_in), fac_out(nsens_out)
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getSensitivities
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)
   if (ierror.lt.0) return

   ! set mirroring for left rail/wheel combination

   sgn      = 1d0
   imodul   = ire_module(ix_reid(ire))
   mirror_y = (imodul.eq.1 .and. (my_ic%config.eq.0 .or. my_ic%config.eq.4))
   if (mirror_y) sgn = -1d0

   if (idebug.ge.3) then
      write(bufout,'(2a,2(i2,a))') trim(pfx_str(subnam,ire,icp)),' return sensitivities(nout=',         &
                lenout,', nin=',lenin,')'
      call write_log(1, bufout)
   endif

   ! set unit scaling factors per input- and output-variable, for scaling towards internal units
   ! e.g. pen(CONTACT) [mm] = fac_in(pen) [mm/m] * pen(external) [m]

   fac_in(1:nsens_in) = 1d0
   fac_in(iin_dpen)   = my_scl%len
   if (my_ic%tang.eq.1) then
      ! T=1, cksi,ceta are shift distances [mm], cphi is shift angle [rad]
      fac_in(iin_dksi1) =       my_scl%len
      fac_in(iin_deta1) = sgn * my_scl%len
      fac_in(iin_dphi1) = sgn * 1d0
   else
      ! T=1, cksi,ceta are relative velocities [-], cphi is angle velocity [rad/mm]
      fac_in(iin_dksi1) =       1d0
      fac_in(iin_deta1) = sgn * 1d0
      fac_in(iin_dphi1) = sgn * 1d0 / my_scl%len
   endif

   ! Fx,Fy are relative forces [-], Mz is force times distance [N.mm]
   fac_out(1:nsens_out) = 1d0
   fac_out(iout_fn)  = my_scl%forc
   fac_out(iout_fx1) =       1d0
   fac_out(iout_fy1) = sgn * 1d0
   fac_out(iout_mz1) = sgn * my_scl%forc * my_scl%len

   ! the sensitivities are d Vout / d Vin, hence the units are [Vout] / [Vin]
   ! the unit conversion is "external = internal / scl"

   do iin = 1, lenin
      do iout = 1, lenout
         ii = iout + (iin-1)*lenout
         if (iout.le.nsens_out .and. iin.le.nsens_in) then
            sens(ii) = gd%outpt1%sens(iout,iin) * fac_in(iin) / fac_out(iout)
         endif
      enddo
   enddo

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getSensitivities

!------------------------------------------------------------------------------------------------------------

subroutine cntc_getCalculationTime(ire, icp, tcpu, twall) &
   bind(c,name=CNAME_(cntc_getcalculationtime))
!--function: return accumulated cpu-time and wall-clock-time used since last timer reset for a contact problem
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   real(kind=8), intent(out) :: tcpu, twall   ! cpu- and wall-clock times used [time]
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_getCalculationTime'
   character(len=256)          :: namtmr
   integer                     :: itimer, numtms, ncontrb, ierror
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_getCalculationTime
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)
   if (ierror.lt.0) return

   itimer = cntc_timer_num(ire, icp)
   call timer_read(itimer, 1, -1, namtmr, numtms, tcpu, twall, ncontrb)

   if (idebug.ge.2) then
      write(bufout,'(2a,2(f7.1,a))') trim(pfx_str(subnam,ire,icp)),' tcpu=',tcpu,', twall=',twall,' [s]'
      call write_log(1, bufout)
   endif
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_getCalculationTime

!------------------------------------------------------------------------------------------------------------

subroutine cntc_resetCalculationTime(ire, icp) &
   bind(c,name=CNAME_(cntc_resetcalculationtime))
!--function: reset the accumulated cpu-time and wall-clock-time used for a contact problem
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_resetCalculationTime'
   integer                     :: itimer, ierror
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_resetCalculationTime
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)
   if (ierror.lt.0) return

   itimer = cntc_timer_num(ire, icp)
   call timer_reset(itimer)

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_resetCalculationTime

!------------------------------------------------------------------------------------------------------------

subroutine subs_getBlocksize(ire, icp, iblk, nx, ny, nz) &
   bind(c,name=CNAME_(subs_getblocksize))
!--function: return the number of points in a block of subsurface points used for a contact problem
!            length of results arrays
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   integer,      intent(in)  :: iblk          ! subsurface stress block number
   integer,      intent(out) :: nx, ny, nz    ! number of positions in x/y/z directions
!--local variables:
   integer    :: ierror
   character(len=*), parameter :: subnam = 'subs_getBlocksize'
#ifdef _WIN32
!dec$ attributes dllexport :: subs_getBlocksize
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)
   if (ierror.lt.0) return

   if (iblk.le.0 .or. iblk.gt.my_subs%nblock) then
      nx = -1
      ny = -1
      nz = -1
   else
      nx = my_subs%blocks(iblk)%nx_eff
      ny = my_subs%blocks(iblk)%ny_eff
      nz = my_subs%blocks(iblk)%nz
   endif

   if (idebug.ge.2) then
      write(bufout,'(2a,i3,a,3i4)') trim(pfx_str(subnam,ire,icp)),' block',iblk,' #points nx,y,z=',     &
                nx, ny, nz
      call write_log(1, bufout)
   endif
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine subs_getBlocksize

!------------------------------------------------------------------------------------------------------------

subroutine subs_getResults(ire, icp, iblk, lenarr, ncol, icol, values) &
   bind(c,name=CNAME_(subs_getresults))
!--function: return selected results of block iblk of subsurface stress calculation for a contact problem
!  icol(ncol)        - requested columns of result data
!                       1-- 3: x,y,z        - positions of points for subsurface stress calculation
!                       4-- 6: ux,uy,uz     - elastic displacements
!                       7-- 9: sighyd,vm,tr - hydrostatic, von Mises and Tresca stresses
!                      10--12: sigma1,2,3   - principal stresses
!                      13--15: sigxx,yx,zx  - components of the full stress tensor
!                      16--18: sigxy,yy,zy
!                      19--21: sigxz,yz,zz
!  values(lenarr,ncol) - result data, with min(lenarr,nx*ny*nz) entries filled per column
!--category: 6, "m=any, cp":        available for modules 1 and 3, in module 1 working on cp data
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: ire           ! result element ID
   integer,      intent(in)  :: icp           ! contact problem ID
   integer,      intent(in)  :: iblk          ! subsurface block ID
   integer,      intent(in)  :: lenarr, ncol  ! sizes of output arrays
   integer,      intent(in)  :: icol(ncol)    ! selected output variables
   real(kind=8), intent(out) :: values(lenarr,ncol)  ! output values for all points of block iblk
!--local variables:
   character(len=*), parameter :: subnam = 'subs_getResults'
   integer                  :: npnt, jcol, ii, ierror
#ifdef _WIN32
!dec$ attributes dllexport :: subs_getResults
#endif

   if (idebug.ge.4) call cntc_log_start(subnam, .true.)
   call cntc_activate(ire, icp, 0, 1, subnam, ierror)
   if (ierror.lt.0) return

   if (idebug.ge.3) then
      write(bufout,'(a,2(a,i3))') trim(pfx_str(subnam,ire,icp)),' return',ncol,                       &
                ' outputs for subs. block',iblk
      call write_log(1, bufout)
   endif

   ! set reference to requested block of subsurface data

   if (iblk.le.0 .or. iblk.gt.my_subs%nblock) then
      write(bufout,'(2a,i3)') trim(pfx_str(subnam,ire,icp)),' no data for subsurf. block',iblk
      call write_log(1, bufout)
      return
   endif

   associate(b => my_subs%blocks(iblk))
   npnt = b%nx_eff * b%ny_eff * b%nz

   if (.not.associated(b%table)) then
      write(bufout,'(2a,i2)') trim(pfx_str(subnam,ire,icp)),' no results available for block',iblk
      call write_log(1, bufout)
   else
      do jcol = 1, ncol
         if (icol(jcol).le.0 .or. icol(jcol).gt.b%ncolum) then
            values(1:lenarr,jcol) = -999d0
         else
            if (icol(jcol).le.6) then
               ! columns 1--3: x,y,z,  columns 4--6: ux,uy,uz
               do ii = 1, min(npnt, lenarr)
                  values(ii,jcol) = b%table(ii,icol(jcol)) / my_scl%len
               enddo
            elseif (icol(jcol).le.21) then
               ! columns 7--21: stresses sighyd/vm/tr, sigma1-3, sigma(i,j)
               do ii = 1, min(npnt, lenarr)
                  values(ii,jcol) = b%table(ii,icol(jcol)) * my_scl%area
               enddo
            endif
         endif
      enddo
   endif

   end associate
   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine subs_getResults

!------------------------------------------------------------------------------------------------------------

subroutine finalize_helper(ire)
!--function: cleanup for one ire, common code for cntc_finalize and cntc_finalizeLast
   implicit none
!--subroutine arguments:
   integer,      intent(in)    :: ire           ! result element ID
!--local variables:
   integer                     :: icp, jre
   character(len=*), parameter :: subnam = 'finalize_helper'

   ! Clean-up ire and decrement the number of active ire's remaining

   if (ix_reid(ire).le.-1) then

      ! skip, should not happen

   else

      ! clean up the gd-structures for all icps of ire

      do icp = 1, MAX_CPids
         if (idebug.ge.4) then
            write(bufout,'(a,i6,a,i6,a)') 'cntc_finalize: destroying gd for ire=',ire,', icp=',icp,'.'
            call write_log(1, bufout)
         endif
         call cntc_destroy_gd(ire, icp)
      enddo

      ! deactivate ire, decrement number of ire's remaining

      ix_reid(ire) = -1
      num_reids = num_reids - 1

      if (idebug.ge.3) then
         write(bufout,'(2a,i3)') trim(pfx_str(subnam,ire,-1)),' number of active result elements=',   &
                num_reids
         call write_log(1, bufout)
      endif

      if (idebug.ge.4) then
         do jre = 1, MAX_REid_id
            if (ix_reid(jre).ge.0) then
               write(bufout,*) 'REid',jre,' is active'
               call write_log(1, bufout)
            endif
         enddo
      endif
   endif

end subroutine finalize_helper

!------------------------------------------------------------------------------------------------------------

subroutine cntc_finalize(ire) &
   bind(c,name=CNAME_(cntc_finalize))
!--function: cleanup for one ire, and in case of last one: call overall finalization
!--category: 5, "m=any, wtd":       available for modules 1 and 3, in module 1 working on wtd data
   implicit none
!--subroutine arguments:
   integer,      intent(in)    :: ire           ! result element ID
!--local variables:
   character(len=*), parameter :: subnam = 'cntc_finalize'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_finalize
#endif

   ! abort if the addon wasn't initialized appropriately

   if (caddon_initialized.le.0) return

#ifdef _OPENMP
   if (omp_in_parallel().ne.0) then
      call write_log('ERROR: cntc_finalize may not be called from a parallel region')
      return
   endif
#endif

   if (idebug.ge.3) then
      write(bufout,'(2a)') trim(pfx_str(subnam,ire,-1)),' finalizing...'
      call write_log(1, bufout)
   endif

   ! Decrement the number of active ire's

   if (ix_reid(ire).le.-1) then
      write(bufout,'(a,i6,a,i6,a)') 'WARNING: Finalize: result element ire=',ire,' is not active (anymore).'
      call write_log(1, bufout)
   else

      ! clean up the gd-structures for ire, decrement number of ire's remaining

      call finalize_helper(ire)

   endif

   ! If this was the last active ire:

   if (num_reids.le.0) call cntc_finalizeLast

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_finalize

!------------------------------------------------------------------------------------------------------------

subroutine cntc_finalizeLast() &
   bind(c,name=CNAME_(cntc_finalizelast))
!--function: clean-up all remaining ire's, print timings, close files, cleanup
!--category: 0, "m=any, glob":      not related to contact's modules 1 or 3, working on global data
   implicit none
!--subroutine arguments:
!--local variables:
   integer                     :: ire
#ifdef _OPENMP
   integer                     :: cur_mxthrd
#endif
   character(len=*), parameter :: subnam = 'cntc_finalizeLast'
#ifdef _WIN32
!dec$ attributes dllexport :: cntc_finalizeLast
#endif

   ! abort if the addon wasn't initialized appropriately

   if (caddon_initialized.le.0) return

#ifdef _OPENMP
   if (omp_in_parallel().ne.0) then
      call write_log('ERROR: cntc_finalizeLast may not be called from a parallel region')
      return
   endif
#endif

   ! disable the thread prefix

   call set_print_thread(-1)

   if (idebug.ge.3) then
      write(bufout,'(2a)') trim(pfx_str(subnam,-1,-1)),': finalizing entirely...'
      call write_log(1, bufout)
   endif

   ! clean up for all ire's that are still active

   do ire = 1, MAX_REid_id
      if (ix_reid(ire).ge.0) then

         ! clean up the gd-structures for ire, decrement number of ire's remaining

         call finalize_helper(ire)

      endif
   enddo

   ! when using parallel processing via OpenMP: print num.threads used

#ifdef _OPENMP
   cur_mxthrd = omp_get_max_threads()
   if (cur_mxthrd.gt.1 .and. idebug.ge.0) then
      write(bufout,'(a,i3,a)')'Parallel run using',cur_mxthrd,' threads'
      call write_log(1, bufout)
   endif
#endif

   ! print timings

   call timer_stop(itimer_main)
   if (idebug.ge.0) call timers_contact_print

   ! close files

   if (inp_open.eq.1) then
      write (linp,*) ' 0  MODULE'
      close (unit=linp)
      inp_open = 0
   endif

   if (out_open.eq.1) then
      close (unit=lout)
      out_open = 0
   endif

   ! allow for a complete new simulation to be started

   caddon_initialized = -1

   if (idebug.ge.4) call cntc_log_start(subnam, .false.)
end subroutine cntc_finalizeLast

!------------------------------------------------------------------------------------------------------------

end module m_contact_addon
