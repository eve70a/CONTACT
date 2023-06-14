!------------------------------------------------------------------------------------------------------------
! m_globals - provide access to global parameters and functions
!
! Copyright 2008-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_globals

! using "public", all generic modules used are exported to all using modules

use m_print_output
use m_licensing
use m_errormsg
use m_timers_contact
use m_blas

implicit none

public

   !---------------------------------------------------------------------------------------------------------

   ! flags and error codes

   include 'caddon_flags.inc'

   !---------------------------------------------------------------------------------------------------------

   ! dimensions and parameters

   integer,      parameter :: MXBLCK  =      10
   integer,      parameter :: MXSUBS  =    9999
   integer,      parameter :: npo     = 2000000
   integer,      parameter :: strlen  =      16
   real(kind=8), parameter :: tiny    =      1d-20
   real(kind=8), parameter :: pi      = 4d0*atan(1d0)

   !---------------------------------------------------------------------------------------------------------
   ! codes for the coordinate directions

   integer, parameter :: ikXDIR = 1, ikYDIR = 2, ikZDIR = 3, ikTANG = -2, ikALL = -3
   integer, parameter :: jkXDIR = 1, jkYDIR = 2, jkZDIR = 3, jkTANG = -2, jkALL = -3

   !---------------------------------------------------------------------------------------------------------

end module m_globals
