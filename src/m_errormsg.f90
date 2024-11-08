!------------------------------------------------------------------------------------------------------------
! m_errormsg - This module provides facilities for error handling.
!
! Copyright 2008-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_errormsg
use m_print_output
implicit none
public

   ! routine(s) used for printing information to the output streams

   public write_errmsg

   ! last_error, last_nline: last error message

   character(len=240)   :: last_error(10)
   integer              :: last_nline

   ! bufout   convenience: provide a generic buffer for preparing output to log-streams, instead of
   !          having to declare this everywhere.

   character(len=240)   :: errmsg(10)
   integer              :: nline_errmsg
   public                  errmsg, nline_errmsg
!$omp threadprivate      ( errmsg, nline_errmsg )

   ! ints, reals, chars: provide generic arrays for error reporting

   integer              :: err_ints(10)
   real(kind=8)         :: err_reals(10)
   character(len=240)   :: err_chars(10)
contains

!------------------------------------------------------------------------------------------------------------

subroutine write_errmsg(nline_arg, errmsg_arg)
!--purpose: Store error message and display to the error output streams
   implicit none
!--subroutine arguments
   integer,          intent(in) :: nline_arg
   character(len=*), intent(in) :: errmsg_arg(nline_arg)
!--local variables
   integer        :: iline

   ! copy error message to last error

   do iline = 1, nline_arg
      last_error(iline) = errmsg_arg(iline)
   enddo
   last_nline = nline_arg

   ! write error message to output stream

   call write_log(nline_arg, errmsg_arg)

end subroutine write_errmsg

!------------------------------------------------------------------------------------------------------------

end module m_errormsg

! -------------- file m_profile:
! prf%ierror:              1  == ierror not set
! profile_values_store:  191  == incorrect itype (is_wheel)
! profile_read_simpack:  ios  == file-system error opening file
! profile_read_simpack:  100  == unexpected end-of-file
! profile_read_simpack:  101  == incorrect nesting of sections
! profile_read_simpack:  102  == unknown section name
! profile_read_simpack:  111  == keyword requires one integer value
! profile_read_simpack:  112  == type must be 0 for a prr-file
! profile_read_simpack:  113  == type must be 1 for a prw-file
! profile_read_simpack:  121  == incorrect data in header section
! profile_read_simpack:  131  == incorrect data in spline section
! profile_read_simpack:  191  == incorrect data for flag header - type
! profile_read_simpack:  182  == profile rotation not yet implemented
! profile_read_miniprof: ios  == file-system error opening file
! profile_read_miniprof: 189  == found keyword after data values on previous lines
! profile_read_miniprof: 190  == invalid profile, at least 2 points needed
! profile_read_miniprof: 191  == incorrect data for flag header - type
! check_inversion:       221  == wheel profile y-values should be decreasing
! check_inversion:       222  == rail profile y-values should be increasing
! check_overall_order:   211  == wheel profile y-values should be descending
! check_overall_order:   212  == rail profile y-values should be ascending
! profile_check_bowtie:  323  == profile has intersections between different segments
! profile_finish_grid:   241  == length of profile curve less than minimum
! profile_finish_grid:   242  == length of profile curve more than maximum
!
! -------------- file m_varprof:
! profile_read_varprf:   151  == incorrect parameters nslc, nintrup, or npart
! profile_read_varprf:   152  == s-positions must be strictly increasing
! profile_read_varprf:   153  == incorrect s_slc in parts-information
! profile_read_varprf:   154  == kink or acceleration out of range
! profile_read_varprf:   155  == incorrect s_slc in parts-information
! profile_read_slice:    301  == incorrect file extension for profile file
!
! -------------- file m_readline:
! readline:               -1  == reached end of file
! readline:               -6  == incorrect format string, i, d, a, l, s
! readline:               -7  == expected element is wrong or missing
! readdbles:              -6  == incorrect format string, d, a
! readdbles:              -2  == expected element is wrong or missing
! readdbles:              -3  == an element of line <> is wrong or missing
! check_sorted:           iv  == values out of order at position iv
!
! -------------- file m_wr_brentmeth:
! wr_contact:       known_err == (/ CNTC_err_allow, CNTC_err_search, CNTC_err_ftot, CNTC_err_tol, 
!                                   CNTC_err_norm, CNTC_err_tang, CNTC_err_icp, CNTC_err_profil, 
!                                   CNTC_err_frclaw /)
!                   my_ierror == CNTC_err_other
! wr_contact_fx_secant:    -1 == NaN in update dfk
! wr_contact_fx_secant:    -2 == residual in fk > tolerance
! brent_its_add_iterate:   -1 == NaN in residual rk
! wr_contact_fxz_brent:    -2 == residual in fk > tolerance
! wr_contact_fxz_broyden:  -1 == NaN in update dfk
! wr_contact_fxz_broyden:  -2 == residual in fk > tolerance
! broyden_solve_dxk:       -3 == zero determinant
!
! -------------- file m_spline_def:
! solve_cubic_segm:        -1 == solution outside segment
! solve_cubic_segm:        maxit == no convergence
! spline_get_s_at_f_spl:   -1 == error -1 in solve_cubic_segm
! locate_one_extremum:     -1 == called for segment with no sign-change
! locate_one_extremum:     -2 == two zeros within segment
