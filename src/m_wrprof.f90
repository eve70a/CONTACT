!------------------------------------------------------------------------------------------------------------
! m_wrprof - driver routine for w/r contact (module 1) in stand-alone CONTACT program
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_wrprof

use m_wrprof_data
use m_wr_input
use m_wr_totforce
implicit none
private

public  wrprof

contains

!------------------------------------------------------------------------------------------------------------

   subroutine wrprof (inp, wtd, linenr)
!--purpose: driver routine for module 1. Performs input and controls Contac in the case of W/R contact.
      implicit none
!--subroutine parameters:
      type(t_ws_track)         :: wtd
      integer   :: inp, linenr
!--local variables:
      integer   :: ierror
      logical   :: lfirst

      lfirst = .true.

      ! Repeat executing cases until R=1 or R=3

      do while(lfirst .or. wtd%ic%return.eq.0 .or. wtd%ic%return.eq.2)
         call timer_start(itimer_wrprof)

         ! Increment the case number

         wtd%meta%ncase = wtd%meta%ncase + 1
         lfirst = .false.

         if (wtd%meta%REid.ge.1 .and. wtd%meta%ncase.le.99999) then
            write(bufout,'(/,a,i6,a,i4)') ' Case',wtd%meta%ncase,', result element',wtd%meta%REid
         elseif (wtd%meta%REid.ge.1) then
            write(bufout,'(/,a,i8,a,i4)') ' Case',wtd%meta%ncase,', result element',wtd%meta%REid
         elseif (                     wtd%meta%ncase.le.99999) then
            write(bufout,'(/,a,i6)')      ' Case',wtd%meta%ncase
         else
            write(bufout,'(/,a,i8)')      ' Case',wtd%meta%ncase
         endif
         call write_log(2, bufout)

         ! read and test new input quantities

         call timer_start(itimer_input)
         call wr_input(wtd%meta%wrkdir, linp, inp, wtd%meta%ncase, linenr, wtd, .true., ierror)
         call timer_stop(itimer_input)

         ! determine the initial contact points, setup & solve the contact problems, output results

         if (inp.eq.2 .and. wtd%ic%return.le.1) call wr_contact(wtd, ierror)

         call timer_stop(itimer_wrprof)
      end do ! end while

   end subroutine wrprof

!------------------------------------------------------------------------------------------------------------

end module m_wrprof
