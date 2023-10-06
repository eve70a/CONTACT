!------------------------------------------------------------------------------------------------------------
! m_timers_contact - timing of compute-intensive parts of the CONTACT computations
!
! Copyright 2008-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_timers_contact
!!--declarations---------------------------------------------------------------------------------------------
use m_timings
use m_print_output
implicit none

public timers_contact_init
public timers_contact_print
public timer_start
public timer_stop

integer, parameter :: itimer_main               = 1
integer, parameter :: itimer_wrprof             = itimer_main      + 1
integer, parameter :: itimer_nonhz              = itimer_wrprof    + 1
integer, parameter :: itimer_subsur             = itimer_nonhz     + 1

integer, parameter :: itimer_profil             = itimer_subsur    + 1
integer, parameter :: itimer_wrgeom             = itimer_profil    + 1
integer, parameter :: itimer_panprc             = itimer_wrgeom    + 1
integer, parameter :: itimer_snorm              = itimer_panprc    + 1
integer, parameter :: itimer_stang              = itimer_snorm     + 1
integer, parameter :: itimer_sens               = itimer_stang     + 1
integer, parameter :: itimer_output             = itimer_sens      + 1

integer, parameter :: itimer_normcg             = itimer_output    + 1
integer, parameter :: itimer_kpec               = itimer_normcg    + 1
integer, parameter :: itimer_analyn             = itimer_kpec      + 1
integer, parameter :: itimer_tangcg             = itimer_analyn    + 1
integer, parameter :: itimer_cnvxgs             = itimer_tangcg    + 1
integer, parameter :: itimer_stdygs             = itimer_cnvxgs    + 1
integer, parameter :: itimer_gdstdy             = itimer_stdygs    + 1
integer, parameter :: itimer_fastsim            = itimer_gdstdy    + 1
integer, parameter :: itimer_temper             = itimer_fastsim   + 1

integer, parameter :: itimer_initcf             = itimer_temper    + 1
integer, parameter :: itimer_input              = itimer_initcf    + 1
integer, parameter :: itimer_locatecp           = itimer_input     + 1
integer, parameter :: itimer_udist              = itimer_locatecp  + 1
integer, parameter :: itimer_sdis               = itimer_udist     + 1
integer, parameter :: itimer_sgencr             = itimer_sdis      + 1
integer, parameter :: itimer_infload            = itimer_sgencr    + 1
integer, parameter :: itimer_filpvs             = itimer_infload   + 1
integer, parameter :: itimer_srznrm             = itimer_filpvs    + 1
integer, parameter :: itimer_srztng             = itimer_srznrm    + 1
integer, parameter :: itimer_files              = itimer_srztng    + 1

integer, parameter :: itimer_fftdsc             = itimer_files     + 1
integer, parameter :: itimer_ffttot             = itimer_fftdsc    + 1

integer, parameter :: itimer_subsinfl           = itimer_ffttot    + 1
integer, parameter :: itimer_subsfft            = itimer_subsinfl  + 1
integer, parameter :: itimer_subsderiv          = itimer_subsfft   + 1
integer, parameter :: itimer_subsfile           = itimer_subsderiv + 1

integer, parameter :: itimer_last               = itimer_subsfile      ! last one of regular timers

! Note: when adding timers, don't forget to add timer names below.             

! Total number of timers
integer            :: num_timers                = itimer_last

! Offset and number of timers for CONTACT add-on
integer, parameter :: my_i0addon                = itimer_last
integer            :: my_naddon

! Level of output w.r.t. timers:
integer, private   :: timer_output              = 1

contains

!------------------------------------------------------------------------------------------------------------

subroutine timers_contact_init(louttm, idebug, maxerr, naddon, i0addon, mxthrd)
!-- HEADER VARIABLES/ARGUMENTS
implicit none
integer, intent(in),  optional :: louttm
integer, intent(in),  optional :: idebug
integer, intent(in),  optional :: maxerr
integer, intent(in),  optional :: naddon
integer, intent(out), optional :: i0addon
integer, intent(in),  optional :: mxthrd

!! executable statements -------------------------------------------------------

   ! Handle optional unit number for print-output of m_timings

   if (present(louttm)) then
      call timer_output_config(louttm_arg=louttm)
   endif

   ! Handle optional level of print-output

   if (present(idebug)) then
      timer_output = idebug
      call timer_output_config(idebug_arg=idebug-2)
   endif

   ! Handle optional error handling for m_timings

   if (present(maxerr)) then
      call timer_error_config(maxerr_arg=maxerr)
   endif

   ! Handle optional flag for CONTACT add-on

   if (present(naddon) .and. present(i0addon)) then
      ! return offset to 1st timer for addon
      i0addon    = my_i0addon
      ! increment number of timers with number requested by addon
      num_timers = my_i0addon + naddon
      my_naddon  = naddon
   endif

   ! set appropriate dimension of table of timers

   if (present(mxthrd)) then
      call timer_table_size(num_timers, 1, mxthrd)
   else
      call timer_table_size(num_timers)
   endif

   ! register all timers for CONTACT

   call timer_name(itimer_main      , namtmr='Total               ')

   call timer_name(itimer_wrprof    , namtmr='Wheel-rail cases    ')
   call timer_name(itimer_nonhz     , namtmr='Non-Hertzian cases  ')
   call timer_name(itimer_subsur    , namtmr='Subsurface points   ')

   call timer_name(itimer_profil    , namtmr='W/r profiles        ')
   call timer_name(itimer_wrgeom    , namtmr='Geometric analysis  ')
   call timer_name(itimer_panprc    , namtmr='Panag.process       ')
   call timer_name(itimer_snorm     , namtmr='Algorithm Norm      ')
   call timer_name(itimer_stang     , namtmr='Algorithm Tang      ')
   call timer_name(itimer_sens      , namtmr='Sensitivities       ')
   call timer_name(itimer_output    , namtmr='Output quantities   ')

   call timer_name(itimer_normcg    , namtmr='Norm: ConjGrd solver')
   call timer_name(itimer_kpec      , namtmr='Norm: KPEC method   ')
   call timer_name(itimer_analyn    , namtmr='Norm: ANALYN method ')
   call timer_name(itimer_tangcg    , namtmr='Tang: TangCG solver ')
   call timer_name(itimer_cnvxgs    , namtmr='Tang: ConvxGS solver')
   call timer_name(itimer_stdygs    , namtmr='Tang: StedyGS solver')
   call timer_name(itimer_gdstdy    , namtmr='Tang: GDstedy solver')
   call timer_name(itimer_fastsim   , namtmr='Tang: Fastsim solver')
   call timer_name(itimer_temper    , namtmr='Temperature calc.   ')

   call timer_name(itimer_initcf    , namtmr='Initial time-step   ')
   call timer_name(itimer_input     , namtmr='Reading input-file  ')
   call timer_name(itimer_locatecp  , namtmr='Contact location    ')
   call timer_name(itimer_udist     , namtmr='Undeformed distance ')
   call timer_name(itimer_sdis      , namtmr='Discretisation      ')
   call timer_name(itimer_sgencr    , namtmr='Influence coeffic.  ')
   call timer_name(itimer_infload   , namtmr='Load infl.coeffic.  ')
   call timer_name(itimer_filpvs    , namtmr='Initial state/estim.')
   call timer_name(itimer_srznrm    , namtmr='Rhs normal problem  ')
   call timer_name(itimer_srztng    , namtmr='Rhs tang. problem   ')
   call timer_name(itimer_files     , namtmr='Writing output-files')

   call timer_name(itimer_fftdsc    , namtmr='FFT descriptors     ')
   call timer_name(itimer_ffttot    , namtmr='FFT calculations    ')

   call timer_name(itimer_subsinfl  , namtmr='Subsurf infl.cf.    ')
   call timer_name(itimer_subsfft   , namtmr='Subsurf FFTs        ')
   call timer_name(itimer_subsderiv , namtmr='Subsurf derived qnt ')
   call timer_name(itimer_subsfile  , namtmr='Subsurf writing     ')

end subroutine timers_contact_init

!------------------------------------------------------------------------------------------------------------

subroutine timers_contact_print
   implicit none
!--local variables
   integer                                :: i
   real(kind=8), dimension(2)             :: cum_time
   real(kind=8)                           :: total_cpu, total_wall
   character(len=20)                      :: namtmr
   integer                                :: numtms, ncontrb
   real(kind=8)                           :: cputim, waltim
   integer                                :: ilines

   ! executable statements -------------------------------------------------------

   ! obtain total cpu-time and wall-clock time for computing percentages

   call timer_read(itimer_main, 1, -1, namtmr, numtms, total_cpu, total_wall, ncontrb)
   total_cpu = max(1e-6, total_cpu)
   total_wall = max(1e-6, total_wall)

   call write_log('Performance timers:')
   call write_log('|------------------------------------------------------------------------|')
   call write_log('|Timer name                     |   cpu time        |    wall clock      |')
   call write_log('|                               |-------------------|--------------------|')
   call write_log('|                     | #times  |   sec     |  %    |    sec     |  %    |')
   call write_log('|------------------------------------------------------------------------|')

   cum_time = 0
   ilines = 0
   do i = 1, num_timers

      if (timer_output.ge.2 .or. i.lt.itimer_initcf .or. i.gt.my_i0addon) then

         ! print separator when reaching first timer of new block

         if (i.eq.itimer_normcg .or. i.eq.itimer_initcf .or. i.eq.my_i0addon+1) then
            if (ilines.gt.0) call write_log('|------------------------------------------------------------------------|')
            ilines = 0
         endif

         ! read data of timer

         call timer_read(i, 1, -1, namtmr, numtms, cputim, waltim, ncontrb)

         ! if timer was used, print its results

         if (numtms.gt.0) then
            if (ncontrb.le.1) then
               write(bufout,111) namtmr, numtms, cputim, 100.*cputim/total_cpu,                         &
                                                 waltim, 100.*waltim/total_wall
            else
               ! write(bufout,112) namtmr, numtms, cputim, 100.*cputim/total_cpu,                         &
               !                                   waltim, 100.*waltim/total_wall, ncontrb
               ! print average per thread, print * to indicate parallel region
               waltim = waltim / real(ncontrb)
               write(bufout,113) namtmr, numtms, cputim, 100.*cputim/total_cpu,                         &
                                                 waltim, 100.*waltim/total_wall
            endif
            call write_log(1, bufout)
            cum_time(1) = cum_time(1) + cputim
            cum_time(2) = cum_time(2) + waltim
            ilines = ilines + 1
         endif

      endif
   enddo
   if (ilines.gt.0) call write_log('|------------------------------------------------------------------------|')

 111 format('|',a,' |',i8,' |',2(f10.1,' | ',f5.1,' | '))
 112 format('|',a,' |',i8,' |',2(f10.1,' | ',f5.1,' | '),'(',i3,')')
 113 format('|',a,' |',i8,' |',2(f10.1,' | ',f5.1,' | '),'(*)')

end subroutine timers_contact_print

end module m_timers_contact
