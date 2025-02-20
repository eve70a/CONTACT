!============================================================================================================
! THE FORTRAN PROGRAM "CONTACT" FOR ROLLING AND SLIDING CONTACT PROBLEMS WITH FRICTION.
!
! Copyright 1979-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!============================================================================================================

program contact
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
   logical            :: flags(mxnval)
   real(kind=8)       :: dbles(mxnval)
   character*256      :: strngs(mxnval)
   integer            :: linenr, nval, idebug, ieof, ierror

   ! Other 'local' variables :

   integer, parameter :: timer_output = 2   !  0 = none, 1 = overview, 2 = full, 3+ = debug
   integer            :: imode, modul, mxthrd

   ! Pointer to the main problem data (-structure)

   type(t_probdata), pointer :: gd
   type(t_ws_track), pointer :: wtd

   ! Create a single global instance "gd" of the hierarchical data-structure,
   ! initialize with a sensible default problem

   allocate(gd)
   call gd_init(gd)

   allocate(wtd)
   call wrprof_init(wtd)

   ! Obtain mode of operation and experiment name, open file-units linp, lout

   call ReadModeExp(imode, gd, linenr)
   wtd%meta%wrkdir  = gd%meta%wrkdir
   wtd%meta%outdir  = gd%meta%outdir
   wtd%meta%expnam  = gd%meta%expnam

   call set_print_streams(nw_outfil=.true., nw_screen=.true., nw_simpck=.false.)

   ! For testing whether screen-output is redirected ok:

   if (.false.) then
      call set_print_streams(nw_outfil=.false., nw_screen=.false., nw_simpck=.false.)
   endif

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

      ! Initialise performance timers

      call timers_contact_init(louttm=lout, idebug=timer_output)
      call timer_start(itimer_main)

      ! Perform a loop for executing all "cases" of the simulation

      modul = 1

      do while (modul.ne.0)

         ! Get the module-number for the next case or sequence of cases

         idebug = 0
         ieof   = -1 ! eof=error
         call readLine(linp, gd%meta%ncase+1, linenr, 'module number', 'i', ints, dbles, flags,         &
                       strngs, mxnval, nval, idebug, ieof, lstop, ierror)
         modul = ints(1)
         if (modul.ne.0 .and. modul.ne.1 .and. modul.ne.3) then
            write(bufout,7200) modul, linenr
            call write_log(2, bufout)
 7200       format (/' ERROR: invalid module number (',i6,') at line',i6,'.')
            if (modul.eq.1) modul = 0
         endif

         ! Start processing the next case
         ! Note: each call to WrProf or NonHz may handle more than one case.

         if (modul.eq.1) then
            call wrprof (imode, wtd, linenr)
            gd%meta%ncase = wtd%meta%ncase
         elseif (modul.eq.3) then
            call nonhz (imode, gd, linenr)
            wtd%meta%ncase = gd%meta%ncase
         endif

      ! end while (modul.ne.0)

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

   ! clean up global datastructure

   call gd_destroy(gd)
   deallocate(gd)
   nullify(gd)
   call wrprof_destroy(wtd)
   deallocate(wtd)
   nullify(wtd)

#if   defined WITH_MKLFFT
   call fft_cleanup()
#endif

   call write_log(' The CONTACT program ended succesfully.')
   stop
   end program contact

!------------------------------------------------------------------------------------------------------------

   subroutine ReadModeExp(imode, gd, linenr)
!--purpose: read mode of operation for entire run, read experiment name, check existence of files
!           and open files.
   use m_hierarch_data
   implicit none
!--subroutine arguments:
   type(t_probdata)           :: gd
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
         call GETARG(2, gd%meta%expnam)
#else
         call get_command_argument(2, gd%meta%expnam, status=ierror)
         if (ierror.lt.0) goto 982
#endif
      else
         write (*, 250)
 250     format(//, ' Give base-name of experiment, default "contact": ',$)
         read (*,'(a)') gd%meta%expnam
         if (gd%meta%expnam.eq.' ') gd%meta%expnam = 'contact'
      endif

      ! strip off the extension .inp from experiment name

      ln = len(trim(gd%meta%expnam))
      if (ln.gt.4 .and. gd%meta%expnam(ln-3:ln).eq.'.inp') then
         write(*,*) '...removing filename extension ".inp" to form experiment name'
         gd%meta%expnam = gd%meta%expnam(1:ln-4)
         ! write(*,*) 'new gd%meta%expnam="',trim(gd%meta%expnam),'"'
      endif

      ! strip off the directory name from experiment name, store as input+output folder

      ix = index_pathsep(gd%meta%expnam, back=.true.)
      if (ix.gt.1) then
         ! write(*,*) '...directory separator found at ix=',ix
         ! write(*,*) '...removing wrkdir to form true experiment name'
         gd%meta%wrkdir = gd%meta%expnam(1:ix-1)
         gd%meta%outdir = gd%meta%wrkdir
         gd%meta%expnam = gd%meta%expnam(ix+1:)
         ! write(*,*) '...wrkdir=',trim(gd%meta%wrkdir)
         ! write(*,*) '...expnam= ',trim(gd%meta%expnam)
      endif

      ! Open input-file <EXPERIM>.inp or <STDIN>, unit linp
      !  - mode 2-3: batch usage, input from <EXPERIM>.inp

      if (imode.ge.2 .and. imode.le.3) then
         fname = trim(gd%meta%expnam) // '.inp'
         if (gd%meta%wrkdir.ne.' ') then
            fname = trim(gd%meta%wrkdir) // path_sep // trim(fname)
         endif

         inp_open = -1
         open (linp, file=fname, status='old', err=985)
         inp_open =  2
         linenr = 0
      endif

      ! Open output-file <EXPERIM>.out, unit lout
      !  - used for flow-output and results in human-readable form

      fname = trim(gd%meta%expnam) // '.out'
      if (gd%meta%outdir.ne.' ') then
         fname = trim(gd%meta%outdir) // path_sep // trim(fname)
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
