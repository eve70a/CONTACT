!============================================================================================================
! THE FORTRAN PROGRAM "CONTACT" FOR ROLLING AND SLIDING CONTACT PROBLEMS WITH FRICTION.
!
! Copyright 1979-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!============================================================================================================

program contact
   use m_global_data
   use m_hierarch_data
   use m_wrprof_data
   use m_nonhz
   use m_wrprof
#ifdef _OPENMP
   use omp_lib      , only : omp_get_max_threads, omp_set_num_threads
#endif
   use m_aijpj
   implicit none

   ! Work variables for writing the CONTACT version information, see VERSION

   integer, parameter :: max_version = 20
   integer            :: num_version, ix
   logical            :: show_error
   character*256      :: version(max_version)

   ! Work variables for reading from stdin or the input-file

   logical, parameter :: lstop = .true.
   integer, parameter :: mxnval = 20
   integer            :: ints(mxnval)
   logical            :: flags(mxnval), zerror
   real(kind=8)       :: dbles(mxnval)
   character*256      :: strngs(mxnval)
   integer            :: linenr, nval, idebug, ieof, ierror

   ! Other 'local' variables :

   integer, parameter :: timer_output = 2   !  0 = none, 1 = overview, 2 = full, 3+ = debug
   integer            :: imode, imodul, REid, CPid, mxthrd


   program_id = program_cexe

   ! Obtain mode of operation and experiment name, open file-units linp, lout

   call ReadModeExp(imode, linenr)

   call set_print_streams(nw_outfil=.true., nw_screen=.true., nw_simpck=.false.)

   ! For testing whether screen-output is redirected ok:

   if (.false.) then
      call set_print_streams(nw_outfil=.false., nw_screen=.false., nw_simpck=.false.)
   endif

   ! Initialise performance timers

   call timers_contact_init(louttm=lout, idebug=timer_output)
   call timer_start(itimer_main)

   ! Initialize the result element data-structures, esp allwtd, allgds, wtd, gd

   ierror = 0
   call cntc_initialize_data()
   call cntc_activate_wtd(0, -1, ierror)
   call cntc_activate_gd(0, 1, ierror)

   ! Check if there is a license present for this version of CONTACT.
   ! Display information from the licensefile, if it is valid.

   my_license%is_valid = .true.
   my_license%licensee_name = 'Open source version'

   ! Include file with version identification, display to the screen

#  include "VERSION"
   show_error = (imode.ge.2)
   call WriteVersion(num_version, version, my_license, show_error)

   ! IMODE = 1: License Management actions

   if (imode.le.1) then

      call write_log(' Mode 1 not available in this version.')

   ! IMODE = 2 & 3: Require a valid license

   elseif (.not.my_license%is_valid) then

   ! IMODE = 2 & 3: Actual processing

   else

      ! Perform a loop for executing all "cases" of the simulation

      imodul = 1

      do while (imodul.ne.0)

         ! Get the module-number for the next case or sequence of cases

         idebug = 0
         ieof   = -1 ! eof=error
         call readLine(linp, my_meta%ncase+1, linenr, 'module number', 'iII', ints, dbles, flags,       &
                       strngs, mxnval, nval, idebug, ieof, lstop, ierror)
         zerror = (ierror.ne.0)

         imodul = ints(1)
         if (imodul.ne.0 .and. imodul.ne.1 .and. imodul.ne.3) then
            zerror = .true.
            write(bufout,7200) imodul, linenr
            call write_log(2, bufout)
 7200       format (/' ERROR: invalid module number (',i6,') at line',i6,'.')
         endif

         ! get optional REid, CPid

         if (nval.ge.2) then    ! if specified, require REid >= 1
            REid = ints(2)
            zerror = zerror .or. .not.check_irng('Result Element number', REid, 1, MAX_REid_id)
         else
            REid = 0
         endif
         if (imodul.eq.1) then
            CPid = -1
         elseif (nval.ge.3) then
            CPid = ints(3)
            zerror = zerror .or. .not.check_irng('Contact patch number', CPid, 1, MAX_CPids)
         else
            CPid = 1            ! default Result element: (REid, CPid) = (0, 1)
         endif

         if (zerror) then

            call write_log(' Errors found, aborting execution')
            imodul = 0

         elseif (imodul.ne.0) then

            ! Activate data-structure for current result element

            if (imodul.eq.1 .or. REid.eq.0) then
               call cntc_activate_wtd(REid, -1, ierror)
            endif
            if (imodul.eq.3 .or. REid.eq.0) then
               call cntc_activate_gd(REid, max(1,CPid), ierror)
            endif
            if (ierror.ne.0) call abort_run()

            ! Start processing the next case
            ! Note: each call to WrProf or NonHz may handle more than one case.

            if (imodul.eq.1) then
               call wrprof (imode, wtd, linenr)
               if (REid.eq.0) gd%meta%ncase = wtd%meta%ncase
            elseif (imodul.eq.3) then
               call nonhz (imode, gd, linenr)
               if (REid.eq.0) wtd%meta%ncase = gd%meta%ncase
            endif
         endif ! zerror

      ! end while (imodul.ne.0)

      enddo

      ! Write timing information to .out-file

      call timer_stop(itimer_main)
      call save_print_streams
      call set_print_streams(nw_outfil=.true.,  nw_screen=.false., nw_simpck=.false.)
#ifdef _OPENMP
      mxthrd = omp_get_max_threads()
#else
      mxthrd = 1
#endif
      if (mxthrd.gt.1) then
         write(bufout,'(/,a,i3,a)') 'Parallel run using', mxthrd, ' threads'
         call write_log(2,bufout)
      endif
      call timers_contact_print
      call restore_print_streams

   endif ! IMODE = 1 / 2 & 3

   ! clean up global datastructures

   call cntc_destroy_wtd(0)
   call cntc_destroy_gd(0, 1)

   do REid = 1, MAX_REid_id
      if (ix_reid(REid).gt.0) then
         if (ire_module(ix_reid(REid)).eq.1) then
            call cntc_destroy_wtd(REid)
         else
            do CPid = 1, MAX_CPids
               call cntc_destroy_gd(REid, CPid)
            enddo
         endif
      endif
   enddo

#if   defined WITH_MKLFFT
   call fft_cleanup()
#endif

   call write_log(' The CONTACT program ended succesfully.')
   stop
   end program contact

!------------------------------------------------------------------------------------------------------------

   subroutine ReadModeExp(imode, linenr)
!--purpose: read mode of operation for entire run, read experiment name, check existence of files
!           and open files.
   use m_global_data
   use m_hierarch_data
   implicit none
!--subroutine arguments:
   integer, intent(out)       :: imode, linenr
!--local variables :
   integer       :: icount, ierror, ln, ix
   character*256 :: fname, tmpstr, progname
   logical       :: z

   ! Handle optional command line arguments

#ifdef PLATF_ppc
   icount = IARGC()
#else
   icount = command_argument_count()
#endif

   !  - first optional argument, the mode-option:

   if (icount.ge.1) then
#ifdef PLATF_ppc
      call GETARG(1, tmpstr)
#else
      call get_command_argument(1, tmpstr, status=ierror)
      if (ierror.lt.0) goto 981
#endif
      read(tmpstr, *, err=981) imode
   endif

   ! repeatedly inquire for the mode of operation when not yet provided

   if (icount.le.0) then
      z = .false.
      do while (.not.z)
         write (*, 100)
 100     format (/,                                                                             &
                ' Mode of operation of this program.'//                                         &
                ' IMODE=1: license management.'/                                                &
                ' IMODE=2: start from input file <EXPERIM>.inp.'/                               &
                ' IMODE=3: check input file <EXPERIM>.inp.')
         read (*,*) imode
         z = (imode.ge.1 .and. imode.le.3)
      enddo
   endif

   if (imode.le.1) then

      ! mode=1: arguments handled in ManageLicense

   else 

      ! mode=2, 3: read the experiment name 

      if (icount.ge.2) then
#ifdef PLATF_ppc
         call GETARG(2, overal_expnam)
#else
         call get_command_argument(2, overal_expnam, status=ierror)
         if (ierror.lt.0) goto 982
#endif
      else
         write (*, 250)
 250     format(//, ' Give base-name of experiment, default "contact": ',$)
         read (*,'(a)') overal_expnam
         if (overal_expnam.eq.' ') overal_expnam = 'contact'
      endif

      ! strip off the extension .inp from experiment name

      ln = len(trim(overal_expnam))
      if (ln.gt.4 .and. overal_expnam(ln-3:ln).eq.'.inp') then
         write(*,*) '...removing filename extension ".inp" to form experiment name'
         overal_expnam = overal_expnam(1:ln-4)
         ! write(*,*) 'new overal_expnam="',trim(overal_expnam),'"'
      endif

      ! strip off the directory name from experiment name, store as input+output folder

      ix = index_pathsep(overal_expnam, back=.true.)
      if (ix.gt.1) then
         ! write(*,*) '...directory separator found at ix=',ix
         ! write(*,*) '...removing wrkdir to form true experiment name'
         overal_wrkdir = overal_expnam(1:ix-1)
         overal_outdir = overal_wrkdir
         overal_expnam = overal_expnam(ix+1:)
         ! write(*,*) '...wrkdir=',trim(overal_wrkdir)
         ! write(*,*) '...expnam= ',trim(overal_expnam)
      endif

      ! Open input-file <EXPERIM>.inp or <STDIN>, unit linp
      !  - mode 2-3: batch usage, input from <EXPERIM>.inp

      if (imode.ge.2 .and. imode.le.3) then
         fname = trim(overal_expnam) // '.inp'
         if (overal_wrkdir.ne.' ') then
            fname = trim(overal_wrkdir) // path_sep // trim(fname)
         endif

         inp_open = -1
         open (linp, file=fname, status='old', err=985)
         inp_open =  2
         linenr = 0
      endif

      ! Open output-file <EXPERIM>.out, unit lout
      !  - used for flow-output and results in human-readable form

      fname = trim(overal_expnam) // '.out'
      if (overal_outdir.ne.' ') then
         fname = trim(overal_outdir) // path_sep // trim(fname)
      endif

      out_open = -1
      open (lout, file=fname, action='write', err=986)
      out_open = 1

   endif

   return

   ! Error handling:

 981  continue
         write(*,*) 'ERROR reading 1st command argument.'
         goto 999

 982  continue
         write(*,*) 'ERROR reading 2nd command argument.'
         goto 999

 983  continue
         write(*,*) 'ERROR reading 3rd command argument.'
         goto 999

 985  continue
         write(*,*) 'ERROR: cannot open input-file: ', trim(fname)
         goto 999

 986  continue
         write(*,*) 'ERROR: cannot open file "', trim(fname), '" for writing'
         goto 999

 999  continue
#ifdef PLATF_ppc
      call GETARG(0, progname)
#else
      call get_command_argument(0, progname, status=ierror)
#endif
      write(*,*) 'USAGE: ',trim(progname), ' [IMODE [EXPNAM]]'
      write(*,*) 'Aborting program execution.'
      call abort_run()

   end subroutine ReadModeExp

!------------------------------------------------------------------------------------------------------------
