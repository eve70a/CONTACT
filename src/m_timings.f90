!-- MODULE m_timings ----------------------------------------------------------------------------------------
!-- DESCRIPTION ---------------------------------------------------------------------------------------------
!
!!  Purpose:
!!  Subroutines that allow detailed measurement of cpu- and wallclock-time used by different portions of 
!!  code of a simulation program.
!
!-- Copyright 2008-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!
!------------------------------------------------------------------------------------------------------------
!!
!! A timer can be compared to a stop-watch.
!!
!! You can reset it, start, stop, start again, stop again, read it, and so on. When you read it you get
!! the accumulated time over all periods that it was started ("active").
!!
!! - our timers measure both cpu-time (cpu-usage) and wallclock-time (elapsed time). Also the number of
!!   times that it was started is reported.
!! - each timer records the contributions of different threads to its timing. The total as well as the
!!   contributions may be queried.
!! - a name/description of up to 20 characters can be assigned to each timer. This name/description is
!!   reported again when the timer is read.
!! - start and stop are implemented using different functions, instead of by a single button. It is a
!!   (non-fatal) error to call start again before stop has been called. A warning will be printed and the
!!   second start will be ignored.
!! - it is possible to read a timer while it is running. This reports the accumulated time of all previous
!!   elapsed periods plus the time of the running period.
!!
!! This module is built around a "large" array of timers: timer_table.
!!
!! The user may define the size of the table, and may decide for himself for which purpose each timer is
!! used. A possible layout is to use timer 1 for the total time for the application, timers 2-9 for the
!! main chunks, 20-29 for further subdivision of timer 2, and so on.
!!
!! The table is two-dimensional: for each timer (first dimension) there may be one or more "instances"
!! (second dimension). The instance-argument of the various routines is optional. These instances can be
!! used when a segment of code is traversed multiple times. For instance to measure times for first and
!! second half time-steps separately. Or to time all iterations of a loop separately and calculate
!! statistics of the timings afterwards.
!!
!! In the implementation a third dimension is used for distinguishing different threads. This is handled
!! transparently.
!
module m_timings

#ifdef _OPENMP
   use omp_lib    , only :  omp_in_parallel, omp_get_thread_num, omp_get_wtime
#endif
implicit none

public  timer_table_size
public  timer_output_config
public  timer_error_config
public  timer_thread_config
private timer_wrong_instn
private timer_not_allocated
public  timer_name
public  timer_start
public  timer_stop
public  timer_read
public  timer_reset_all
public  timer_reset

!! Logical unit-number of stdout, used for error messages and debug info
!! (LOUTTM may be redirected to a file, * is not); see timer_output_config.
!! (Output may be disabled by using a negative value, LOUTTM = -1.)
integer            :: LOUTTM = 6

! Level of debug-output of these routines, 0==none; see timer_output_config.
integer            :: idebug = 2

! Stop program when more than maxerr errors occurred; see timer_error_config
integer            :: maxerr = 0
integer            :: numerr = 0

private

! definition of one timer:
type, private :: t_one_timer
   character(len=20) :: namtmr  !! name of the timer, used for output
   real(kind=8)      :: totcpu  !! accumulated cpu-time for this timer
   real(kind=8)      :: totwal  !! accumulated wall-clock time for this timer
   integer           :: numtms  !! number of times timer has been called and times have been added to
                                !! totcpu and totwal
   integer           :: isactv  !! flag indicating that a new timing is active (has been started) (1)
                                !! or not (0)
   real(kind=8)      :: stacpu  !! cpu-time at start of current (active) timing
   real(kind=8)      :: stawal  !! walltime at start of current (active) timing
end type t_one_timer

! maximum extent of 3D table (particularly for print-statements)
integer, parameter :: MAX_NTIMER =99999
integer, parameter :: MAX_NINSTN =  999
integer, parameter :: MAX_NTHREAD=  999

! definition of the entire 3D table of timers:
type(t_one_timer), dimension(:,:,:), pointer :: timer_table => NULL()

! size of the allocated table
integer            :: ntimer_tbl  = -1
integer            :: ninstn_tbl  = -1
integer            :: nthread_tbl = -1

! explicit thread number set by calling environment; -1 == automatic, using omp_get_thread_num()
integer            :: tmr_thread   = -1
!$omp threadprivate ( tmr_thread )

contains

!-- SUBROUTINE timer_table_size -----------------------------------------------------------------------------
!-- DESCRIPTION ---------------------------------------------------------------------------------------------
!!
!!  Set the size of the timer_table, allocate the table, and initialize all timers.
!!
!!> On entry:
!!  ntimer       number of timers to be used in the table
!!  ninstn_arg   optional: number instances per timer to be used (default 1)
!!  nthread_arg  optional: maximal number of threads to be supported (default 1)
!!
!!  On return:
!!< -            the timer_table memory has been allocated and initialized
!!
subroutine timer_table_size(ntimer, ninstn_arg, nthread_arg)

!-- HEADER VARIABLES/ARGUMENTS
   implicit none
   integer, intent(in)           :: ntimer    
   integer, intent(in), optional :: ninstn_arg
   integer, intent(in), optional :: nthread_arg

!-- LOCAL VARIABLES
   ! loop counters, actual number of instances
   integer :: itimer, iinstn, ithread, ninstn, nthread
   type(t_one_timer), dimension(:,:,:), pointer :: tmp_table

   ! This routine changes global data, and may not be called from parallel regions

#ifdef _OPENMP
   if (omp_in_parallel().ne.0) return
#endif

   ! Handle optional arguments

   if (present(ninstn_arg)) then
      ninstn = ninstn_arg
   else
      ninstn = 1
   endif
   if (present(nthread_arg)) then
      nthread = nthread_arg
   else
      nthread = 1
   endif

   ! Check that requested size of the table is valid

   if (ntimer.le.0 .or. ninstn.le.0 .or. nthread.le.0) then
      write(LOUTTM,*) 'timer_table_size: Error: dimensions must be >0!'
      write(LOUTTM,*) '   ntimer=', ntimer,', ninstn=',ninstn,', nthread=',nthread
      numerr = numerr + 1
      if (maxerr.ge.0 .and. numerr.gt.maxerr) stop
   endif

   if (ntimer.gt.MAX_NTIMER .or. ninstn.gt.MAX_NINSTN .or. nthread.gt.MAX_NTHREAD) then
      write(LOUTTM,*) 'timer_table_size: Error: requested dimensions too large!'
      write(LOUTTM,*) '   ntimer= ', ntimer, ', max=',MAX_NTIMER
      write(LOUTTM,*) '   ninstn= ', ninstn, ', max=',MAX_NINSTN
      write(LOUTTM,*) '   nthread=', nthread,', max=',MAX_NTHREAD
      numerr = numerr + 1
      if (maxerr.ge.0 .and. numerr.gt.maxerr) stop
   endif

   ! write(LOUTTM,*) 'timer_table_size: creating table with ntimer=', ntimer,', ninstn=',ninstn,          &
   !                  ', nthread=',nthread

   ! No action required if the table was allocated before with the right size.

   if (.not.associated(timer_table) .or. ntimer.ne.ntimer_tbl .or. ninstn.ne.ninstn_tbl .or.            &
                                                                      nthread.ne.nthread_tbl) then

   ! If the new size is different than the previous size:
      ! - allocate temporary table with new size

      ! write(LOUTTM,*) 'timer_table_size: allocate tmp_table'
      allocate(tmp_table(ntimer,ninstn,0:nthread-1))

      ! - initialize all timers and instances

      do ithread = 0, nthread-1
         do iinstn = 1, ninstn
            do itimer = 1, ntimer
               tmp_table(itimer,iinstn,ithread)%namtmr = ' '
               tmp_table(itimer,iinstn,ithread)%numtms = 0
               tmp_table(itimer,iinstn,ithread)%totcpu = 0d0
               tmp_table(itimer,iinstn,ithread)%totwal = 0d0
               tmp_table(itimer,iinstn,ithread)%isactv = 0
            enddo
         enddo
      enddo

      ! - copy data from previous table, if available

      if (associated(timer_table)) then
         ! write(LOUTTM,*) 'timer_table_size: timer_table already allocated, copying data'
         do ithread = 0, min(nthread, nthread_tbl)-1
            do iinstn = 1, min(ninstn, ninstn_tbl)
               do itimer = 1, min(ntimer, ntimer_tbl)
                  tmp_table(itimer,iinstn,ithread) = timer_table(itimer,iinstn,ithread)
               enddo
            enddo
         enddo
      endif

      ! - de-allocate previous table, if available

      if (associated(timer_table)) then
         ! write(LOUTTM,*) 'timer_table_size: timer_table already allocated, de-allocating'
         deallocate(timer_table)
      endif

      ! - set pointer to new table

      ! write(LOUTTM,*) 'timer_table_size: set pointer to tmp_table'
      timer_table => tmp_table
      ntimer_tbl  = ntimer
      ninstn_tbl  = ninstn
      nthread_tbl = nthread

   else

      ! write(LOUTTM,*) 'timer_table_size: previous table is of the correct size.'

   endif

end subroutine timer_table_size
!------------------------------------------------------------------------------------------------------------


!-- SUBROUTINE timer_output_config --------------------------------------------------------------------------
!-- DESCRIPTION ---------------------------------------------------------------------------------------------
!!
!!  Configure the unit number -output and amount of print output.
!!
!!> On entry:
!!  louttm      optional: logical unit number for print output
!!  idebug      optional: level of print output
!!
!!  On return:
!!< -           the configuration has been stored
!!
subroutine timer_output_config(louttm_arg, idebug_arg)

!-- HEADER VARIABLES/ARGUMENTS
   implicit none
   integer, intent(in), optional :: louttm_arg
   integer, intent(in), optional :: idebug_arg

   ! This routine changes global data, and may not be called from parallel regions

#ifdef _OPENMP
   if (omp_in_parallel().ne.0) return
#endif

   ! Handle optional unit number for print-output

   if (present(louttm_arg)) LOUTTM = louttm_arg

   ! Handle optional level of print-output

   if (present(idebug_arg)) idebug = idebug_arg

end subroutine timer_output_config
!------------------------------------------------------------------------------------------------------------


!-- SUBROUTINE timer_error_config --------------------------------------------------------------------------
!-- DESCRIPTION ---------------------------------------------------------------------------------------------
!!
!!  Configure the behavior on errors
!!
!!> On entry:
!!  maxerr      maximum number of errors allowed, when exceeded program execution is aborted
!!              default 0: first error stops program
!!              set to -1: program will never be stopped.
!!
!!  On return:
!!< -           the configuration has been stored
!!
subroutine timer_error_config(maxerr_arg)

!-- HEADER VARIABLES/ARGUMENTS
   implicit none
   integer, intent(in), optional :: maxerr_arg

   ! This routine changes global data, and may not be called from parallel regions

#ifdef _OPENMP
   if (omp_in_parallel().ne.0) return
#endif

   maxerr = maxerr_arg

end subroutine timer_error_config
!------------------------------------------------------------------------------------------------------------


!-- SUBROUTINE timer_thread_config --------------------------------------------------------------------------
!-- DESCRIPTION ---------------------------------------------------------------------------------------------
!!
!!  Set the thread number for the calling thread
!!
!!> On entry:
!!  tmr_thread   thread number
!!
!!  On return:
!!< -           the configuration has been stored
!!
subroutine timer_thread_config(tmr_thread_arg)

!-- HEADER VARIABLES/ARGUMENTS
   implicit none
   integer, intent(in) :: tmr_thread_arg

   tmr_thread = tmr_thread_arg

end subroutine timer_thread_config
!------------------------------------------------------------------------------------------------------------


!-- SUBROUTINE timer_wrong_instn ----------------------------------------------------------------------------
!-- DESCRIPTION ---------------------------------------------------------------------------------------------
!! 
!!  Handle an incorrect timer- or instance number
!!
!!> On entry:
!!  subnam      name of calling subroutine
!!  itimer      number of requested timer
!!  iinstn      number of requested instance
!!
!!  On return:
!!< -           appropriate error handling is done
!!
subroutine timer_wrong_instn(subnam, itimer, iinstn)

!-- HEADER VARIABLES/ARGUMENTS
   implicit none
   character(len=*)              :: subnam
   integer, intent(in)           :: itimer    
   integer, intent(in)           :: iinstn

   ! This routine mainly provides a critical section for accessing I/O and global data

!$omp critical (timer_wrong_instn_print)
   if (LOUTTM.ge.1) write(LOUTTM,91) trim(subnam), ': Error: requested timer (',itimer,',',iinstn,       &
                  ') is out of bounds; ntimer=',ntimer_tbl,', ninstn=',ninstn_tbl,'.'
91 format(2a,i5,a,i3,a,i10,a,i5)
   numerr = numerr + 1
   if (maxerr.ge.0 .and. numerr.gt.maxerr) stop
!$omp end critical (timer_wrong_instn_print)

end subroutine timer_wrong_instn
!------------------------------------------------------------------------------------------------------------


!-- SUBROUTINE timer_not_allocated --------------------------------------------------------------------------
!-- DESCRIPTION ---------------------------------------------------------------------------------------------
!!
!!  Error handling for incorrect initialization
!!
!!> On entry:
!!  subnam      name of calling subroutine
!!
!!  On return:
!!< -           appropriate error handling is done
!!
subroutine timer_not_allocated(subnam)

!-- HEADER VARIABLES/ARGUMENTS
   implicit none
   character(len=*)              :: subnam

   ! This routine mainly provides a critical section for accessing I/O and global data

!$omp critical (timer_not_allocated_print)
   if (LOUTTM.ge.1) write(LOUTTM,*) trim(subnam),': Error: timer_table has not been allocated yet!'
   numerr = numerr + 1
   if (maxerr.ge.0 .and. numerr.gt.maxerr) stop
!$omp end critical (timer_not_allocated_print)

end subroutine timer_not_allocated
!------------------------------------------------------------------------------------------------------------


!-- SUBROUTINE timer_name -----------------------------------------------------------------------------------
!-- DESCRIPTION ---------------------------------------------------------------------------------------------
!!
!!  Set the name of timer (stop-watch) with coordinates (itimer, iinstn_arg)
!!
!!> On entry:
!!  itimer      number of timer to be started
!!  iinstn_arg  optional: instance of timer to be started, default 1
!!  namtmr      name/description of the timer
!!
!!  On return:
!!< -           the name has been registered in the timer_table
!!
subroutine timer_name(itimer, iinstn_arg, namtmr)

!-- HEADER VARIABLES/ARGUMENTS
   implicit none
   integer,          intent(in)           :: itimer    
   integer,          intent(in), optional :: iinstn_arg 
   character(len=*), intent(in)           :: namtmr

!-- LOCAL VARIABLES
   integer      :: iinstn        ! actual instance
   integer      :: jthread       ! loop variable for threads

   ! This routine changes global data, and may hence not be called from parallel regions

#ifdef _OPENMP
   if (omp_in_parallel().ne.0) return
#endif

   ! Handle optional iinstn-argument

   if (present(iinstn_arg)) then
      iinstn = iinstn_arg
   else
      iinstn = 1
   endif

   ! Check that table is allocated and that requested timer is within bounds

   if (.not.associated(timer_table)) call timer_not_allocated('timer_name')
   if (itimer.le.0 .or. itimer.gt.ntimer_tbl .or. iinstn.le.0 .or. iinstn.gt.ninstn_tbl)                &
      call timer_wrong_instn('timer_name', itimer, iinstn)

   ! register the name/description of the timer in timer_table

   do jthread = 0, nthread_tbl-1
      timer_table(itimer,iinstn,jthread)%namtmr = namtmr
   enddo

end subroutine timer_name
!------------------------------------------------------------------------------------------------------------


!-- SUBROUTINE timer_start ----------------------------------------------------------------------------------
!-- DESCRIPTION ---------------------------------------------------------------------------------------------
!!
!!  Start the timer (stop-watch) with coordinates (itimer, iinstn_arg)
!!
!!> On entry:
!!  itimer      number of timer to be started
!!  iinstn_arg  optional: instance of timer to be started, default 1
!!
!!  On return:
!!< -           the start-time has been registered in the timer_table
!!
subroutine timer_start(itimer, iinstn_arg)

!-- HEADER VARIABLES/ARGUMENTS
   implicit none
   integer, intent(in)           :: itimer    
   integer, intent(in), optional :: iinstn_arg 

!-- LOCAL VARIABLES
   integer                    :: iinstn         !  actual instance
   integer                    :: ithread        !  actual thread number
   real(kind=8)               :: curcpu, curwal !  current cpu/wallclock-time

   ! Handle optional iinstn-argument, get current thread number

   if (present(iinstn_arg)) then
      iinstn = iinstn_arg
   else
      iinstn = 1
   endif
#ifdef _OPENMP
   if (tmr_thread.ge.0) then
      ithread = tmr_thread
   else
      ithread = omp_get_thread_num()
   endif
#else
   ithread = 0
#endif

   ! Check that table is allocated and that requested timer is within bounds

   if (.not.associated(timer_table)) call timer_not_allocated('timer_start')
   if (itimer.le.0 .or. itimer.gt.ntimer_tbl .or. iinstn.le.0 .or. iinstn.gt.ninstn_tbl)                &
      call timer_wrong_instn('timer_start', itimer, iinstn)

   ! Ignore timing if the actual thread number is higher than the max.threads of the table

   if (ithread.ge.nthread_tbl) return

   ! Check that no timing is active yet for the requested timer

   if (timer_table(itimer,iinstn,ithread)%isactv.gt.0) then

!$omp critical (timer_start_print)
      if (LOUTTM.ge.1) write(LOUTTM,92) 'timer_start: Error: requested timer (',itimer,',',             &
                     iinstn,',',ithread,') has already been started; current start will be ignored.'
92    format(a,i5,a,i3,a,i3,a)
!$omp end critical (timer_start_print)

   else

      ! get current cpu/wallclock time and store in timer_table

      call my_cpu_time(curcpu)
      call wall_time(curwal)
      timer_table(itimer,iinstn,ithread)%isactv = 1
      timer_table(itimer,iinstn,ithread)%stacpu = curcpu
      timer_table(itimer,iinstn,ithread)%stawal = curwal

   endif
      
end subroutine timer_start
!------------------------------------------------------------------------------------------------------------


!-- SUBROUTINE timer_stop -----------------------------------------------------------------------------------
!-- DESCRIPTION ---------------------------------------------------------------------------------------------
!!
!!  Stop the timer (stop-watch) with coordinates (itimer, iinstn_arg)
!!
!!> On entry:
!!  itimer      number of timer to be stopped
!!  iinstn_arg  optional: instance of timer to be stopped, default 1
!!
!!  On return:
!!  -           the elapsed (cpu and wallclock) time since the start of the timing has been added to the
!!<             timer in the timer_table
!!
subroutine timer_stop(itimer, iinstn_arg)

!-- HEADER VARIABLES/ARGUMENTS
   implicit none
   integer, intent(in)           :: itimer    
   integer, intent(in), optional :: iinstn_arg  

!-- LOCAL VARIABLES
   integer                    :: iinstn         !  actual instance
   integer                    :: ithread        !  actual thread number
   real(kind=8)               :: curcpu, curwal !  current cpu/wallclock-time
   type(t_one_timer), pointer :: timer          !  pointer to current timer

   ! Handle optional iinstn-argument, get current thread number

   if (present(iinstn_arg)) then
      iinstn = iinstn_arg
   else
      iinstn = 1
   endif
#ifdef _OPENMP
   if (tmr_thread.ge.0) then
      ithread = tmr_thread
   else
      ithread = omp_get_thread_num()
   endif
#else
   ithread = 0
#endif

   ! Check that table is allocated and that requested timer is within bounds

   if (.not.associated(timer_table)) call timer_not_allocated('timer_stop')
   if (itimer.le.0 .or. itimer.gt.ntimer_tbl .or. iinstn.le.0 .or. iinstn.gt.ninstn_tbl)                &
      call timer_wrong_instn('timer_stop', itimer, iinstn)

   ! Ignore timing if the actual thread number is higher than the max.threads of the table

   if (ithread.ge.nthread_tbl) return

   ! Check that a timing is active for the requested timer

   if (timer_table(itimer,iinstn,ithread)%isactv.le.0) then

!$omp critical (timer_stop_print)
      if (LOUTTM.ge.1) write(LOUTTM,92) 'timer_stop: Error: requested timer (',itimer,',',iinstn,       &
                  ',',ithread,') has not been started; current stop will be ignored.'
92    format(a,i5,a,i3,a,i3,a)
!$omp end critical (timer_stop_print)

   else

      ! get current cpu/wallclock time; compute elapsed time since start and
      ! add to totals in timer_table

      call my_cpu_time(curcpu)
      call wall_time(curwal)

      timer => timer_table(itimer,iinstn,ithread)
      timer%isactv = 0
      timer%numtms = timer%numtms + 1
      timer%totcpu = timer%totcpu + curcpu - timer%stacpu
      timer%totwal = timer%totwal + curwal - timer%stawal 

   endif
      
end subroutine timer_stop
!------------------------------------------------------------------------------------------------------------


!-- SUBROUTINE timer_read -----------------------------------------------------------------------------------
!-- DESCRIPTION ---------------------------------------------------------------------------------------------
!!
!!  Return the name and times for timer (stop-watch) with coordinates (itimer, iinstn_arg, ithread_arg)
!!
!!> On entry:
!!  itimer       number of timer to be read
!!  iinstn_arg   optional: instance of timer to be read, default 1
!!  ithread_arg  optional: thread number of timer to be read, default -1 = aggregate of all threads
!!
!!  On return:
!!  namtmr      name/description of timer
!!  numtms      number of timings performed with this timer; number of times that the timer was started
!!              and stopped since the timer was initialized/reset.
!!  cputim      total amount of CPU-time used in the timings with this timer
!!  waltim      total amount of wallclock-time elapsed in the timings with this timer
!!  ncontrb     number of threads contributing to the results.
!!
subroutine timer_read(itimer, iinstn_arg, ithread_arg, namtmr, numtms, cputim, waltim, ncontrb)

!-- HEADER VARIABLES/ARGUMENTS
   implicit none
   integer,          intent(in)            :: itimer    
   integer,          intent(in),  optional :: iinstn_arg  
   integer,          intent(in),  optional :: ithread_arg  
   character(len=*), intent(out), optional :: namtmr
   integer,          intent(out), optional :: numtms
   real(kind=8),     intent(out), optional :: cputim
   real(kind=8),     intent(out), optional :: waltim
   integer,          intent(out), optional :: ncontrb

!-- LOCAL VARIABLES
   integer      :: iinstn                ! actual instance
   integer      :: ithread               ! actual thread number
   real(kind=8) :: curcpu, curwal        ! current cpu/wallclock-time
   integer      :: jthread, jthread0, jthread1 ! variables for loop over threads
   integer      :: my_numtms, my_ncontrb
   real(kind=8) :: my_cputim, my_waltim

   ! Handle optional arguments

   if (present(iinstn_arg)) then
      iinstn = iinstn_arg
   else
      iinstn = 1
   endif
   if (present(ithread_arg)) then
      ithread = ithread_arg
   else
      ithread = -1
   endif

   ! Initialize the outputs

   my_numtms = 0
   my_cputim = 0d0
   my_waltim = 0d0
   my_ncontrb = 0

   ! If ithread<0 this routine changes global data, and may then not be called from parallel regions

#ifdef _OPENMP
   if (ithread.lt.0 .and. omp_in_parallel().ne.0) return
#else
   ithread = 0
#endif

   ! Ignore requests if the actual thread number is higher than the max.threads of the table

   if (ithread.ge.nthread_tbl) return

   ! Check that table is allocated and that requested timer is within bounds

   if (.not.associated(timer_table)) call timer_not_allocated('timer_read')
   if (itimer.le.0 .or. itimer.gt.ntimer_tbl .or. iinstn.le.0 .or. iinstn.gt.ninstn_tbl)                &
      call timer_wrong_instn('timer_read', itimer, iinstn)

   ! Get range of threads for which data are to be aggregated

   if (ithread.lt.0) then
      jthread0 = 0
      jthread1 = nthread_tbl-1
   else
      jthread0 = ithread
      jthread1 = ithread
   endif

   ! Get current accumulated timings from timer_table

   if (present(namtmr)) namtmr = timer_table(itimer,iinstn,jthread0)%namtmr 

   do jthread = jthread0, jthread1

      my_numtms = my_numtms + timer_table(itimer,iinstn,jthread)%numtms 
      my_cputim = my_cputim + timer_table(itimer,iinstn,jthread)%totcpu 
      my_waltim = my_waltim + timer_table(itimer,iinstn,jthread)%totwal 

      ! If a timing is active for the requested timer, get time and add to accumulated time

      if (timer_table(itimer,iinstn,jthread)%isactv.gt.0) then
         call my_cpu_time(curcpu)
         call wall_time(curwal)

         my_numtms = my_numtms + 1
         my_cputim = my_cputim + curcpu - timer_table(itimer,iinstn,jthread)%stacpu
         my_waltim = my_waltim + curwal - timer_table(itimer,iinstn,jthread)%stawal 
      endif

      ! Increment number of contributing threads when appropriate

      if (timer_table(itimer,iinstn,jthread)%numtms.gt.0 .or.                                           &
          timer_table(itimer,iinstn,jthread)%isactv.gt.0) my_ncontrb = my_ncontrb + 1

   enddo

   if (present(numtms)) numtms = my_numtms
   if (present(cputim)) cputim = my_cputim
   if (present(waltim)) waltim = my_waltim
   if (present(ncontrb)) ncontrb = my_ncontrb

end subroutine timer_read
!------------------------------------------------------------------------------------------------------------


!-- SUBROUTINE timer_reset_all ------------------------------------------------------------------------------
!-- DESCRIPTION ---------------------------------------------------------------------------------------------
!!
!!  Reset all timers (stop-watches) of the table
!!
!!> On entry:
!!  -
!!
!!  On return:
!!< -           the timers have been reset
!!
subroutine timer_reset_all()

!-- HEADER VARIABLES/ARGUMENTS
   implicit none
!-- LOCAL VARIABLES
   integer      :: itimer, iinstn        ! actual timer, instance

   ! This routine changes global data, and may hence not be called from parallel regions

#ifdef _OPENMP
   if (omp_in_parallel().gt.0) return
#endif

   ! Loop over all timers/instances, delegate actual work to timer_reset

   do itimer = 1, ntimer_tbl
      do iinstn = 1, ninstn_tbl
         call timer_reset(itimer, iinstn)
      enddo
   enddo

end subroutine timer_reset_all
!------------------------------------------------------------------------------------------------------------


!-- SUBROUTINE timer_reset ----------------------------------------------------------------------------------
!-- DESCRIPTION ---------------------------------------------------------------------------------------------
!!
!!  Reset the timer (stop-watch) with coordinates (itimer, iinstn, ithread)
!!
!!> On entry:
!!  itimer      number of timer to be reset
!!  iinstn      optional: instance of timer to be reset
!!  ithread     optional: instance of timer to be reset, default = -1 = all threads
!!
!!  On return:
!!< -           the timer has been reset
!!
subroutine timer_reset(itimer, iinstn_arg, ithread_arg)

!-- HEADER VARIABLES/ARGUMENTS
   implicit none
   integer, intent(in)           :: itimer    
   integer, intent(in), optional :: iinstn_arg  
   integer, intent(in), optional :: ithread_arg  

!-- LOCAL VARIABLES
   integer      :: iinstn        ! actual instance
   integer      :: ithread       ! actual thread
   integer      :: jthread, jthread0, jthread1 ! variables for loop over threads

   ! Handle optional arguments

   if (present(iinstn_arg)) then
      iinstn = iinstn_arg
   else
      iinstn = 1
   endif
   if (present(ithread_arg)) then
      ithread = ithread_arg
   else
      ithread = -1
   endif

   ! If ithread<0 this routine changes global data, and may then not be called from parallel regions

#ifdef _OPENMP
   if (ithread.lt.0 .and. omp_in_parallel().ne.0) return
#else
   ithread = 0
#endif

   ! Ignore request if the actual thread number is higher than the max.threads of the table

   if (ithread.ge.nthread_tbl) return

   ! Check that table is allocated and that requested timer is within bounds

   if (.not.associated(timer_table)) call timer_not_allocated('timer_reset')
   if (itimer.le.0 .or. itimer.gt.ntimer_tbl .or. iinstn.le.0 .or. iinstn.gt.ninstn_tbl)                &
      call timer_wrong_instn('timer_reset', itimer, iinstn)

   ! Get range of threads for which timers are to be reset

   if (ithread.lt.0) then
      jthread0 = 0
      jthread1 = nthread_tbl-1
   else
      jthread0 = ithread
      jthread1 = ithread
   endif

   do jthread = jthread0, jthread1

      ! Check that no timing is active for the requested timer

      if (timer_table(itimer,iinstn,jthread)%isactv.gt.0) then
!$omp    critical (timer_reset_print)
         if (LOUTTM.ge.1) write(LOUTTM,92) 'timer_reset: Warning: requested timer (',itimer,',',           &
            iinstn,',',jthread, ') is active; timer will be reset anyway.'
92       format(a,i5,a,i3,a,i3,a)
!$omp    end critical (timer_reset_print)
      endif

      ! Reset requested timer in timer_table

      timer_table(itimer,iinstn,jthread)%numtms = 0
      timer_table(itimer,iinstn,jthread)%totcpu = 0d0
      timer_table(itimer,iinstn,jthread)%totwal = 0d0
      timer_table(itimer,iinstn,jthread)%isactv = 0

   enddo

end subroutine timer_reset
!------------------------------------------------------------------------------------------------------------


!-- SUBROUTINE my_cpu_time ----------------------------------------------------------------------------------
!-- DESCRIPTION ---------------------------------------------------------------------------------------------
!!+
!!  Deliver the estimated total cpu time used by the current thread
!!
!!> On entry:
!!  -
!!
!!  On return:
!!< curcpu      current cpu time used
!!
subroutine my_cpu_time(curcpu)

!-- HEADER VARIABLES/ARGUMENTS
   implicit none
   real(kind=8), intent(out) :: curcpu

   ! The system function cpu_time() delivers the total cpu-time used by a process, i.e. the aggregate of
   ! cpu-times for all threads that are used. This is ok for our timers that are used in serial regions.
   ! However, in parallel regions it means that the cpu-time spent by one thread is compounded in the
   ! timers of all other threads. In that case it's better to return the wall-clock time instead.

#ifdef _OPENMP
   if (omp_in_parallel().ne.0 .or. tmr_thread.gt.0) then

      curcpu = omp_get_wtime()

   else
#endif

      call cpu_time(curcpu)

#ifdef _OPENMP
   endif
#endif

end subroutine my_cpu_time


!-- SUBROUTINE wall_time ------------------------------------------------------------------------------------
!-- DESCRIPTION ---------------------------------------------------------------------------------------------
!!+
!!  Deliver the elapsed wallclock time in seconds with respect to an unknown
!!  but fixed reference
!!
!!> On entry:
!!  -
!!
!!  On return:
!!< curwal      current wallclock time
!!
subroutine wall_time(curwal)

!-- HEADER VARIABLES/ARGUMENTS
   implicit none
   real(kind=8), intent(out) :: curwal

!-- LOCAL VARIABLES
   logical,      save :: is_initialized = .false.
   real(kind=8), save :: mulfac_to_secs  ! conversion-factor from clock-ticks to seconds
   integer,      save :: num_periods = 0 ! number of times clock has been reset during program execution
   real(kind=8), save :: secs_per_period ! the amount of time per period in which the clock is reset once
   integer,      save :: last_count = 0  ! value of the system_clock in the previous call to this routine,
                                         ! for detecting that the system_clock was reset
!$omp threadprivate(is_initialized, mulfac_to_secs, num_periods, secs_per_period, last_count)
   integer            :: curr_ticks, tick_rate, tick_max ! auxiliary variables for calling the system_clock
   real(kind=8)       :: rem_secs                        ! auxiliary variables for debug-output

   ! Initialization of this routine:
   !  - determine the measuring-unit of subroutine system_clock (e.g. 10000 ticks per second)
   !  - derive the multiplication factor for conversion to seconds (e.g. 0.00001)
   !  - determine the maximum value of this clock (e.g. 2147483647)
   !  - derive the number of seconds after which the timer tops over (e.g. 214748, which is about 59 hours)

   if (.not.is_initialized) then
      call system_clock(curr_ticks, tick_rate, tick_max)
      rem_secs = (real(tick_max,8) - real(curr_ticks,8)) / real(tick_rate,8)
      mulfac_to_secs  = 1d0 / real(tick_rate,8)
      secs_per_period = tick_max * mulfac_to_secs
      is_initialized  = .true.
      if (LOUTTM.ge. 1 .and. idebug.ge.1) then
!$omp    critical (wall_time_init)
         if (idebug.ge.2) write(LOUTTM,*) 'curr_ticks=',curr_ticks,', tick_max=', tick_max,             &
            'tick_rate=',tick_rate
         if (idebug.ge.1) write(LOUTTM,*) 'time to timer-roll-over is approx.', int(rem_secs/3600.),    &
            ' hr', nint( (rem_secs-3600.*int(rem_secs/3600.))/60. ), ' min'
!$omp    end critical (wall_time_init)
      endif
   else
      call system_clock(curr_ticks)
   endif

   ! Detect whether the timer has been reset since the last time that this routine was called. If this
   ! routine is called often enough (e.g. at least once in every 24 hours) then all resets will be detected.

   if (curr_ticks.lt.last_count) then
      num_periods = num_periods + 1
      if (LOUTTM.ge.1 .and. idebug.ge.2) then
!$omp    critical (wall_time_rollover)
         write(LOUTTM,*) 'wall_time: timer roll-over. curr_ticks=', curr_ticks, ', last_count=',last_count
         write(LOUTTM,*) 'num_periods=',num_periods
!$omp    end critical (wall_time_rollover)
      endif
   endif
   last_count = curr_ticks

   ! Compute the number of seconds since the last time the system clock was
   ! reset before the application was started: the current number of ticks plus
   ! the number of times the clock was reset during program execution

   curwal = real(curr_ticks,8) * mulfac_to_secs + num_periods * secs_per_period

end subroutine wall_time
!------------------------------------------------------------------------------------------------------------

end module m_timings
!------------------------------------------------------------------------------------------------------------
