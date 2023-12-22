!------------------------------------------------------------------------------------------------------------
! m_nonhz - driver routine for module 3, basic (non-Hertzian) contact
!
! Copyright 1979-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------

module m_nonhz

use m_hierarch_data
use m_scontc
use m_sinput
use m_soutpt

implicit none

private
public  nonhz

contains

!------------------------------------------------------------------------------------------------------------

subroutine nonhz (inp, gd, linenr)
!--purpose: Performs input and controls Contac in the case of Non-Hertzian contact.
   implicit none
!--subroutine parameters
   type(t_probdata)         :: gd
   integer                  :: inp, linenr
!--local variables:
   logical       :: lfirst
   integer       :: ierror

   lfirst = .true.

   ! Repeat executing cases until R=1 or R=3

   do while(lfirst .or. (gd%ic%return.eq.0 .or. gd%ic%return.eq.2))
      call timer_start(itimer_nonhz)

      ! Increment the case number

      gd%meta%ncase = gd%meta%ncase + 1
      lfirst = .false.

      write(bufout,1001) gd%meta%ncase
      call write_log(2, bufout)
 1001 format (/' Case',i6)

      ! read and test new input quantities

      call timer_start(itimer_input)
      call input (inp, linenr, gd%meta, gd%ic, gd%potcon_inp, gd%potcon_cur, gd%mater, gd%hertz,        &
                  gd%cgrid_inp, gd%geom, gd%influ, gd%fric, gd%kin, gd%solv, gd%subs)
      call timer_stop(itimer_input)

      ! call contac to solve this case

      if (inp.eq.2) call contac (gd, ierror)

   ! end while(lfirst or R=0 or R=2)

      call timer_stop(itimer_nonhz)
   end do

end subroutine nonhz

!------------------------------------------------------------------------------------------------------------

end module m_nonhz
