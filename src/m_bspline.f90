!------------------------------------------------------------------------------------------------------------
! m_bspline - interfaces regarding B-splines
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_bspline
   use m_globals
   use m_spline
   implicit none
   private

   ! Debugging for module m_bspline

   integer  :: ldebug    =  0    ! local level of debugging
   integer  :: ii_debug  = -1    ! output point for which detailed info is requested (-1 = none)
   integer  :: iel_debug = -1    ! input element for which detailed info is requested (-1 = none)
   public  bspline_set_debug

   ! global parameter for B-spline knots

   real(kind=8), parameter    :: tiny_dt = 1d-10
   public tiny_dt

   ! Data type for 2D B-splines:

   public t_bspline2d

   !---------------------------------------------------------------------------------------------------------
   ! m_bspline_def: Basic functions for B-splines:
   !---------------------------------------------------------------------------------------------------------

   public  bspline2d_nullify
   public  bspline2d_destroy
   public  bspline2d_print
   public  bspline_make_breakpoints
   public  bspline_insert_knot
   public  bspline_eval1d
   public  bspline_make_ppform
   public  bspline_check_kink_accel
   public  bspline_check_mask

   !---------------------------------------------------------------------------------------------------------
   ! m_bspline_make: Functions for the creation of 1D B-splines:
   !---------------------------------------------------------------------------------------------------------

   public  bspline_make_knot_vector_atmeas
   public  bspline_make_knot_vector_simple
   public  bspline_make_knot_vector_advanced
   public  bspline_make_knot_vector

   public  bspline_make_bt_w_b
   public  bspline_make_bt_w_xyz
   public  bspline_make_c_and_d
   public  bspline_make_dt_c_d
   public  bspline_reduce_bmat
   public  bspline_reduce_dmat
   public  bspline_expand_vector

   public  bspline_solve_coef1d_smooth
   public  bspline_solve_coef1d_intpol

   public  bspline_make_2d_bspline

   !---------------------------------------------------------------------------------------------------------
   ! m_bspline_get: Functions for the creation & evaluation of 2D B-splines:
   !---------------------------------------------------------------------------------------------------------

   public  bspline_eval2d_list
   public  bspline_eval2d_prod
   private bspline_eval2d_inverse
   public  bspline_eval2d_inverse_list
   public  bspline_eval2d_inverse_prod

   private bspline_make_1d_ppspline_phase1
   private bspline_make_1d_ppspline_phase2
   private bspline_make_1d_ppspline_phase3
   public  bspline_make_1d_ppspline

   public  bspline_get_ppcoef_1seg_orig
   public  bspline_get_ppcoef_1seg_modf

   ! inverse for half-parametric spline x==u, (x,v) -> (y,z)

   public  bspline_get_z_at_xy
   public  bspline_get_z_at_xy_list
   public  bspline_get_z_at_xy_prod

   interface bspline_get_z_at_xy
      module procedure bspline_get_z_at_xy_list
      module procedure bspline_get_z_at_xy_prod
   end interface bspline_get_z_at_xy

   ! inverse for full-parametric spline (u,v) -> (x,y,z)

   public  bspline_get_xz_at_uy
   public  bspline_get_xz_at_uy_list
   public  bspline_get_xz_at_uy_prod

   interface bspline_get_xz_at_uy
      module procedure bspline_get_xz_at_uy_list
      module procedure bspline_get_xz_at_uy_prod
   end interface bspline_get_xz_at_uy

   public  test_bspline

!------------------------------------------------------------------------------------------------------------
!  data for a parametric tensor B-spline (u, v, [x(u,v)], y(u,v), z(u,v)):

   type :: t_bspline2d
      integer      :: nknotu, nknotv, nreducu, nreducv, nsplu, nsplv, nbrku, nbrkv
      logical      :: has_xdata

      real(kind=8), dimension(:),   pointer  :: tui    => NULL() ! (nknotu)
      real(kind=8), dimension(:),   pointer  :: tvj    => NULL() ! (nknotv)
      logical,      dimension(:),   pointer  :: keeptu => NULL() ! (nknotu)
      logical,      dimension(:),   pointer  :: keeptv => NULL() ! (nknotv)

      real(kind=8), dimension(:),   pointer  :: ubrk   => NULL() ! (nbrku)
      real(kind=8), dimension(:),   pointer  :: vbrk   => NULL() ! (nbrkv)

      real(kind=8), dimension(:,:), pointer  :: cij_x  => NULL() ! (nsplu,nsplv)
      real(kind=8), dimension(:,:), pointer  :: cij_y  => NULL() ! (nsplu,nsplv)
      real(kind=8), dimension(:,:), pointer  :: cij_z  => NULL() ! (nsplu,nsplv)
      integer,      dimension(:,:), pointer  :: mask   => NULL() ! (nsplu,nsplv)

      ! nknotu   number of knots used in u-direction
      ! nknotv   number of knots used in v-direction
      ! tui      knots in u-direction
      ! tvj      knots in v-direction
      ! keeptu   flags indicating kept knots in u-direction
      ! keeptv   flags indicating kept knots in v-direction
      ! nreducu  number of knots kept in u-direction after reduction
      ! nreducv  number of knots kept in v-direction after reduction

      ! has_xdata  flag indicating that cij_x is present (3d space curve or 3d surface)
      ! nsplu    number of splines used in u-direction
      ! nsplv    number of splines used in v-direction
      ! cij_x    tensor B-spline coefficients for x(u,v)
      ! cij_y    tensor B-spline coefficients for y(u,v)
      ! cij_z    tensor B-spline coefficients for z(u,v)
      ! mask     flags indicating which cij are available (1) or missing (0)

      ! nbrku    number of break points in u-direction
      ! nbrkv    number of break points in v-direction
      ! ubrk     break points in u-direction
      ! vbrk     break points in v-direction

   end type t_bspline2d

!------------------------------------------------------------------------------------------------------------

   interface

      !------------------------------------------------------------------------------------------------------
      ! m_bspline_def: Basic functions for B-splines:
      !------------------------------------------------------------------------------------------------------

      module subroutine bspline_set_debug(new_ldebug, new_ii_debug, new_iel_debug)
      !--function: enable/disable debug output of spline routines
      !--subroutine arguments:
         integer, intent(in)           :: new_ldebug       ! level of debug output required
         integer, intent(in), optional :: new_ii_debug     ! specific point of interest for debugging
         integer, intent(in), optional :: new_iel_debug    ! specific point of interest for debugging
      end subroutine bspline_set_debug

      !------------------------------------------------------------------------------------------------------

      module subroutine bspline2d_nullify(bspl)
      !--function: initialize B-spline structure, nullify pointers
      !--subroutine parameters:
         type(t_bspline2d) :: bspl
      end subroutine bspline2d_nullify

      !------------------------------------------------------------------------------------------------------

      module subroutine bspline2d_destroy(bspl)
      !--function: clean-up allocated arrays for B-spline, nullify alias-pointers
      !--subroutine parameters:
         type(t_bspline2d)  :: bspl
      end subroutine bspline2d_destroy

      !------------------------------------------------------------------------------------------------------

      module subroutine bspline2d_print(bspl, nam, idebug, ndigit)
      !--function: print information on spline bspl
      !--subroutine arguments
         type(t_bspline2d)    :: bspl
         character(len=*)     :: nam
         integer              :: idebug
         integer, optional    :: ndigit       ! number of significant digits
      end subroutine bspline2d_print

      !------------------------------------------------------------------------------------------------------

      module subroutine bspline_make_breakpoints(nknot, tj, tiny_dt, nbrk_arg, spnt, idebug)
      !--function: determine breakpoints in knot vector tj: unique knots, excluding 1:k-1 and end-[0:k-2]
      !--subroutine arguments:
         integer,      intent(in)     :: nknot, nbrk_arg, idebug
         real(kind=8), intent(in)     :: tj(nknot), tiny_dt
         real(kind=8), intent(out)    :: spnt(nbrk_arg)
      end subroutine bspline_make_breakpoints

      !------------------------------------------------------------------------------------------------------

      module subroutine bspline_insert_knot(nknot, tj, tnew, tjx, nrow, ncol, coef)
      !--function: shift and update B-spline coefficients for new knot tnew
      !--subroutine arguments:
         integer,      intent(in)           :: nknot, nrow, ncol
         real(kind=8), intent(in)           :: tj(nknot), tnew
         real(kind=8), intent(inout)        :: coef(nrow,ncol)
         real(kind=8), intent(out)          :: tjx(nknot+1)
      end subroutine bspline_insert_knot

      !------------------------------------------------------------------------------------------------------

      module subroutine bspline_eval1d(nknot, tj, tiny_dt, npnt, si, jseg, b1, b2, b3, b4, idebug,      &
                                my_ierror)
      !--function: evaluate B-spline basis functions B_{j,k} for knot-vector tj at sample locations si
      !--subroutine arguments:
         integer,      intent(in)    :: nknot, npnt, idebug
         real(kind=8), intent(in)    :: tj(nknot), tiny_dt, si(npnt)
         integer,      intent(out)   :: jseg(npnt)
         real(kind=8), intent(out)   :: b1(npnt,1), b2(npnt,2), b3(npnt,3), b4(npnt,4)
         integer,      intent(out)   :: my_ierror
      end subroutine bspline_eval1d

      !------------------------------------------------------------------------------------------------------

      module subroutine bspline_make_ppform(nknot, tj, tiny_dt, nbrk, jseg, b1, b2, b3, b4, nspl,       &
                                coef_xyz, spl, my_ierror)
      !--function: create PP-form of B-spline, evaluating {ay0-ay3}, {az0-az3} at start of each segment
      !--subroutine arguments:
         integer,      intent(in)           :: nknot, nbrk, nspl
         integer,      intent(in)           :: jseg(nbrk)
         real(kind=8), intent(in)           :: tj(nknot), tiny_dt,                                      &
                                               b1(nbrk,1), b2(nbrk,2), b3(nbrk,3), b4(nbrk,4),          &
                                               coef_xyz(nspl,3)
         type(t_spline)                     :: spl
         integer,      intent(out)          :: my_ierror
      end subroutine bspline_make_ppform

      !------------------------------------------------------------------------------------------------------

      module subroutine bspline_check_kink_accel(nprf, nkink, ikinks, naccel, iaccel, my_ierror)
      !--function: check requirements of B-spline method on kinks and accelerations
      !--subroutine arguments:
         integer,      intent(in)     :: nprf, nkink, naccel
         integer,      intent(in)     :: ikinks(nkink), iaccel(naccel)
         integer,      intent(out)    :: my_ierror
      end subroutine bspline_check_kink_accel

      !------------------------------------------------------------------------------------------------------

      module subroutine bspline_check_mask(nmeasu, nmeasv, mask, ierror)
      !--function: check requirements on `missing parts' encoded in mask array
      !--subroutine arguments:
         integer                            :: nmeasu, nmeasv
         logical,      intent(in)           :: mask(nmeasu,nmeasv)
         integer,      intent(out)          :: ierror
      end subroutine bspline_check_mask

      !------------------------------------------------------------------------------------------------------
      ! m_bspline_make: Functions for the creation of 1D B-splines:
      !------------------------------------------------------------------------------------------------------

      module subroutine bspline_make_knot_vector_atmeas(nmeas, xi, use_insert, use_repl, nknot, tx,     &
                                namcoor)
      !--function: make basic knot vector with knots at measurement sites
      !            .not.use_insert: skip 2nd and n-1th for not-a-knot b.c.
      !            use_repl: replicate 1st and nth knots, .not.use_repl: replicate 1st/last intervals
      !--subroutine arguments:
         integer,          intent(in)       :: nmeas
         real(kind=8),     intent(in)       :: xi(nmeas)
         logical,          intent(in)       :: use_insert, use_repl
         integer,          intent(out)      :: nknot
         real(kind=8),     intent(out)      :: tx(:)
         character(len=*), intent(in)       :: namcoor
         end subroutine bspline_make_knot_vector_atmeas
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_make_knot_vector_simple(nprf, sprf, nkink, ikinks, naccel, iaccel,      &
                              ds_out, maxknot, nknot, tj, idebug, my_ierror)
      !--function: determine knot-vector with target step ds_out and repeated knots at kinks and
      !            accelerations; simple version that just applies ds_out as much as possible
      !--subroutine arguments:
         integer,      intent(in)     :: nprf, nkink, naccel, maxknot, idebug
         integer,      intent(in)     :: ikinks(nkink), iaccel(naccel)
         real(kind=8), intent(in)     :: sprf(nprf), ds_out
         real(kind=8), intent(out)    :: tj(maxknot)
         integer,      intent(out)    :: nknot, my_ierror
      end subroutine bspline_make_knot_vector_simple
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_make_knot_vector_advanced(nprf, sprf, nkink, ikinks, naccel, iaccel,    &
                              ds_out, maxknot, nknot, tj, idebug, my_ierror)
      !--function: determine knot-vector with target step ds_out and repeated knots at kinks and 
      !            accelerations; advanced knot-vector with at most one knot per data-interval, with
      !            measurements skipped for not-a-knot conditions at boundaries, kinks and accelerations
      !--subroutine arguments:
         integer,      intent(in)     :: nprf, nkink, naccel, maxknot, idebug
         integer,      intent(in)     :: ikinks(nkink), iaccel(naccel)
         real(kind=8), intent(in)     :: sprf(nprf), ds_out
         real(kind=8), intent(out)    :: tj(maxknot)
         integer,      intent(out)    :: nknot, my_ierror
      end subroutine bspline_make_knot_vector_advanced
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_make_knot_vector(nprf, sprf, nkink, ikinks, naccel, iaccel, use_simple, &
                              ds_out, maxknot, nknot, tj, keepj, nreduc, idebug, my_ierror)
      !--function: determine knot-vector with target step ds_out and repeated knots at kinks and
      !            accelerations
      !--subroutine arguments:
         integer,      intent(in)     :: nprf, nkink, naccel, maxknot, idebug
         integer,      intent(in)     :: ikinks(nkink), iaccel(naccel)
         real(kind=8), intent(in)     :: sprf(nprf), ds_out
         logical,      intent(in)     :: use_simple
         real(kind=8), intent(out)    :: tj(maxknot)
         logical,      intent(out)    :: keepj(maxknot)
         integer,      intent(out)    :: nknot, nreduc, my_ierror
      end subroutine bspline_make_knot_vector
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_make_bt_w_b(nspl, nmeas, jseg, wgt, b4, btwb, idebug)
      !--function: compute matrix product B^T * W * B
      !--subroutine arguments:
         integer,      intent(in)           :: nspl, nmeas, idebug
         integer,      intent(in)           :: jseg(nmeas)
         real(kind=8), intent(in)           :: wgt(nmeas), b4(nmeas,4)
         real(kind=8), intent(out)          :: btwb(4,nspl)
      end subroutine bspline_make_bt_w_b
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_make_bt_w_xyz(nspl_all, nspl_rdc, nmeas, jseg, wgt, b4, xi, yi, zi,     &
                      has_xdata, btw_xyz, idebug)
      !--function: compute matrix product B^T * W * [x,y,z]
      !            note: btw_xyz has space for all coefficients but is filled for reduced system
      !--subroutine arguments:
         integer,      intent(in)           :: nspl_all, nspl_rdc, nmeas, idebug
         integer,      intent(in)           :: jseg(nmeas)
         real(kind=8), intent(in)           :: xi(nmeas), yi(nmeas), zi(nmeas), wgt(nmeas), b4(nmeas,4)
         logical,      intent(in)           :: has_xdata
         real(kind=8), intent(out)          :: btw_xyz(nspl_all,3)
      end subroutine bspline_make_bt_w_xyz
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_make_c_and_d(nknot, tj, tiny_dt, cmat, dmat)
      !--function: form the matrices C and D
      !--subroutine arguments:
         integer,      intent(in)           :: nknot
         real(kind=8), intent(in)           :: tj(nknot), tiny_dt
         real(kind=8), intent(out)          :: cmat(nknot-1,1), dmat(nknot-1,4)
      end subroutine bspline_make_c_and_d
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_make_dt_c_d(nknot, keepj, nspl_rdc, cmat, dmat, dtcd, idebug)
      !--function: compute matrix product D^T * C * D
      !--subroutine arguments:
         integer,      intent(in)           :: nknot, nspl_rdc, idebug
         logical,      intent(in)           :: keepj(nknot)
         real(kind=8), intent(in)           :: cmat(nknot-1,1), dmat(nknot-1,4)
         real(kind=8), intent(out)          :: dtcd(4,nspl_rdc)
      end subroutine bspline_make_dt_c_d
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_reduce_bmat(nknot, keepj, nmeas, jseg, b4)
      !--function: compute matrix product B4 * R, combining columns for segments with reduced order,
      !            esp. for linear segments between double kinks
      !--subroutine arguments:
         integer,      intent(in)           :: nknot, nmeas
         logical,      intent(in)           :: keepj(nknot)
         integer,      intent(inout)        :: jseg(nmeas)
         real(kind=8), intent(inout)        :: b4(nmeas,4)
      end subroutine bspline_reduce_bmat
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_reduce_dmat(nknot, keepj, dmat, idebug)
      !--function: compute matrix product D * R, combining columns for segments with reduced order,
      !            esp. for linear segments between double kinks
      !--subroutine arguments:
         integer,      intent(in)           :: nknot, idebug
         logical,      intent(in)           :: keepj(nknot)
         real(kind=8), intent(inout)        :: dmat(nknot-1,4)
      end subroutine bspline_reduce_dmat
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_expand_vector(nknot, keepj, nreduc, nspl_all, coef, idebug)
      !--function: expand input coefficients for reduced basisfunctions coef^* to full coefficients coef
      !--subroutine arguments:
         integer,      intent(in)           :: nknot, nreduc, nspl_all, idebug
         logical,      intent(in)           :: keepj(nknot)
         real(kind=8), intent(inout)        :: coef(nspl_all,3)
      end subroutine bspline_expand_vector
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_solve_coef1d_smooth(nspl_all, nspl_rdc, btwb, lambda, dtcd, btw_xyz, &
                                                                                      idebug, my_ierror)
      !--function: solve the linear systems (B^T W B + lambda D^T C D) [ coef_[xyz] ] = B^T W [x, y, z]
      !            note: btw_xyz has space for all coeff but is used for reduced system
      !--subroutine arguments:
         integer,      intent(in)           :: nspl_all, nspl_rdc, idebug
         real(kind=8), intent(in)           :: lambda
         real(kind=8), intent(in)           :: btwb(4,nspl_rdc), dtcd(4,nspl_rdc)
         real(kind=8), intent(inout)        :: btw_xyz(nspl_all,3)  ! on input:  [BtW*x, BtW*y, BtW*z]
                                                                    ! on output: [coef_x, coef_y, coef_z]
         integer,      intent(out)          :: my_ierror
      end subroutine bspline_solve_coef1d_smooth
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_solve_coef1d_intpol(nsplx, jseg, b4, nrow, nsply, cx, cy, cz,           &
                                                                          has_xdata, idebug, my_ierror)
      !--function: solve the linear systems B [cx2d, cy2d, cz2d] = [cx1d, cy1d, cz1d]
      !--subroutine arguments:
         integer,      intent(in)           :: nsplx, nrow, nsply, idebug        ! nrow >= nsplx
         integer,      intent(in)           :: jseg(nsplx)
         logical,      intent(in)           :: has_xdata
         real(kind=8), intent(in)           :: b4(nsplx,4)
         real(kind=8), intent(inout)        :: cx(nrow,nsply), cy(nrow,nsply), cz(nrow,nsply)
                                                              ! on input:  c[xyz]1d, on output: c[xyz]2d
         integer,      intent(out)          :: my_ierror
      end subroutine bspline_solve_coef1d_intpol
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine get_next_interval(n, mask, iofs, ix0, ix1)
      !--function: locate a range ix0:ix1 of true-values in mask(iofs:n), if not found, set ix1 < ix0
      !--subroutine arguments:
         integer, intent(in)    :: n, iofs
         logical, intent(in)    :: mask(n)
         integer, intent(out)   :: ix0, ix1
      end subroutine get_next_interval
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_make_2d_bspline(spl2d, nmeasu, nmeasv, ui, vj, x2d, y2d, z2d,           &
                                                               has_xdata, use_approx, my_ierror, mask)
      !--function: compute parametric 2D tensor interpolating B-spline for data (x2d, y2d, z2d),
      !            store in B-form (cij_x, cij_y, cij_z).
      !--subroutine arguments:
         type(t_bspline2d)                  :: spl2d
         integer,      intent(in)           :: nmeasu, nmeasv
         logical,      intent(in)           :: has_xdata, use_approx
         real(kind=8), intent(in)           :: ui(nmeasu), vj(nmeasv)
         real(kind=8), intent(in)           :: x2d(nmeasu,nmeasv), y2d(nmeasu,nmeasv), z2d(nmeasu,nmeasv)
         integer,      intent(out)          :: my_ierror
         logical,      intent(in), optional :: mask(nmeasu,nmeasv)
      end subroutine bspline_make_2d_bspline

      !------------------------------------------------------------------------------------------------------
      ! m_bspline_get: Functions for the creation & evaluation of 2D B-splines:
      !------------------------------------------------------------------------------------------------------

      module subroutine bspline_eval2d_list(spl2d, nout, uout, vout, xout, yout, zout, has_xout,        &
                                my_ierror, exterval)
      !--function: evaluate 2D tensor B-spline at nout positions (uout,vout) producing (xout,yout,zout)
      !--subroutine arguments:
         type(t_bspline2d)                           :: spl2d
         integer,      intent(in)                    :: nout
         real(kind=8), intent(in)                    :: uout(nout), vout(nout)
         real(kind=8), intent(out)                   :: xout(nout), yout(nout), zout(nout)
         logical,      intent(in)                    :: has_xout         ! xout may be a dummy array
         integer,      intent(out)                   :: my_ierror
         real(kind=8), intent(in), optional          :: exterval
      end subroutine bspline_eval2d_list
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_eval2d_prod(spl2d, nuout, nvout, uout, vout, xout, yout, zout, has_xout, &
                              my_ierror, exterval)
      !--function: evaluate 2D tensor B-spline at nuout x nvout positions producing (xout,yout,zout)
      !--subroutine arguments:
         type(t_bspline2d)                  :: spl2d
         integer,      intent(in)           :: nuout, nvout
         real(kind=8), intent(in)           :: uout(nuout), vout(nvout)
         real(kind=8), intent(out)          :: xout(nuout,nvout), yout(nuout,nvout), zout(nuout,nvout)
         logical,      intent(in)           :: has_xout         ! xout may be a dummy array
         integer,      intent(out)          :: my_ierror
         real(kind=8), intent(in), optional :: exterval
      end subroutine bspline_eval2d_prod
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_eval2d_inverse(spl2d, noutu, nouty, noutv, uout, yout, vout, my_ierror)
      !--function: inverse evaluation of 2D tensor B-spline, determining vout at given (uout,yout)
      !--subroutine arguments:
         type(t_bspline2d)                  :: spl2d
         integer,      intent(in)           :: noutu, nouty, noutv
         real(kind=8), intent(in)           :: uout(noutu), yout(nouty)
         real(kind=8), intent(out)          :: vout(noutu,noutv)
         integer,      intent(out)          :: my_ierror
      end subroutine bspline_eval2d_inverse
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_eval2d_inverse_list(spl2d, nout, uout, yout, vout, my_ierror)
      !--function: inverse evaluation of 2D tensor B-spline, determining vout at given (uout,yout)
      !--subroutine arguments:
         type(t_bspline2d)                  :: spl2d
         integer,      intent(in)           :: nout
         real(kind=8), intent(in)           :: uout(nout), yout(nout)
         real(kind=8), intent(out)          :: vout(nout)
         integer,      intent(out)          :: my_ierror
      end subroutine bspline_eval2d_inverse_list
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_eval2d_inverse_prod(spl2d, noutu, nouty, uout, yout, vout, my_ierror)
      !--function: inverse evaluation of 2D tensor B-spline, determining vout at given (uout,yout)
      !--subroutine arguments:
         type(t_bspline2d)                  :: spl2d
         integer,      intent(in)           :: noutu, nouty
         real(kind=8), intent(in)           :: uout(noutu), yout(nouty)
         real(kind=8), intent(out)          :: vout(noutu,nouty)
         integer,      intent(out)          :: my_ierror
      end subroutine bspline_eval2d_inverse_prod
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_make_1d_ppspline_phase1(bspl, nmeas, s_prf, ds_out, lambda, nkink,      &
                                ikinks, naccel, iaccel, my_ierror)
      !--function: compute parametric least squares smoothing B-spline for 1-d grid (x(s),y(s),z(s)),
      !            stored in PP-form as used in t_spline.
      !            phase 1: preparations for knot-vector
      !--subroutine arguments:
         type(t_bspline2d)                  :: bspl
         integer,      intent(in)           :: nmeas
         real(kind=8), intent(inout)        :: s_prf(nmeas)
         real(kind=8), intent(inout)        :: ds_out     ! target spacing for B-spline
         real(kind=8), intent(inout)        :: lambda     ! weight of 2nd derivative, relative to data wgt
         integer,      intent(in)           :: nkink
         integer,      intent(in)           :: ikinks(nkink)  ! kinks = jump in 1st derivative
         integer,      intent(in)           :: naccel
         integer,      intent(in)           :: iaccel(naccel) ! accelerations = jump in 2nd derivative
         integer,      intent(out)          :: my_ierror
      end subroutine bspline_make_1d_ppspline_phase1
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_make_1d_ppspline_phase2(bspl, nspl_rdc, lambda, wgt,                    &
                              nmeas, s_prf, x_prf, y_prf, z_prf, has_xdata, btw_xyz, my_ierror)
      !--function: compute parametric least squares smoothing B-spline for 1-d grid (x(s),y(s),z(s)),
      !            stored in PP-form as used in t_spline.
      !            phase 2: set up (BTWB + l*DTCD) = BTWxyz, solve B-spline coefficients
      !--subroutine arguments:
         type(t_bspline2d)                  :: bspl
         integer,      intent(in)           :: nspl_rdc, nmeas
         real(kind=8), intent(in)           :: lambda     ! weight of 2nd derivative, relative to data wgt
         real(kind=8), intent(in)           :: wgt(nmeas) ! weights, e.g. using spacing ds
         real(kind=8), intent(in)           :: s_prf(nmeas), x_prf(nmeas), y_prf(nmeas), z_prf(nmeas) 
         logical,      intent(in)           :: has_xdata  ! whether x_prf is a dummy or needs to be used
         real(kind=8), intent(out)          :: btw_xyz(bspl%nsplu,3)
         integer,      intent(out)          :: my_ierror
      end subroutine bspline_make_1d_ppspline_phase2
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_make_1d_ppspline_phase3(bspl, coef_xyz, ppspl, my_ierror)
      !--function: compute parametric least squares smoothing B-spline for 1-d grid (x(s),y(s),z(s)),
      !            stored in PP-form as used in t_spline.
      !            phase 3: convert knots + B-spline coefficients to PP-form
      !--subroutine arguments:
         type(t_bspline2d)                  :: bspl
         type(t_spline)                     :: ppspl
         real(kind=8), intent(in)           :: coef_xyz(bspl%nsplu,3)
         integer,      intent(out)          :: my_ierror
      end subroutine bspline_make_1d_ppspline_phase3
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_make_1d_ppspline(ppspl, nmeas, s_prf, x_prf, y_prf, z_prf, has_xdata,   &
                                ds_out_arg, lambda_arg, use_wgt, nkink, ikinks, naccel, iaccel,         &
                                my_ierror, wgt_arg)
      !--function: compute parametric least squares smoothing B-spline for 1-d grid (x(s),y(s),z(s)),
      !            stored in PP-form as used in t_spline.
      !            phase 3: convert knots + B-spline coefficients to PP-form
      !--subroutine arguments:
         type(t_spline)                     :: ppspl
         integer                            :: nmeas
         real(kind=8), intent(inout)        :: s_prf(nmeas), x_prf(nmeas), y_prf(nmeas), z_prf(nmeas) 
         logical,      intent(in)           :: has_xdata  ! whether x_prf is a dummy or needs to be used
         real(kind=8), intent(in)           :: ds_out_arg ! target spacing for B-spline
         real(kind=8), intent(in)           :: lambda_arg ! weight of 2nd derivative, relative to data wgt
         logical,      intent(in)           :: use_wgt    ! weigh data with spacing ds
         integer,      intent(in)           :: nkink
         integer,      intent(in)           :: ikinks(nkink)  ! kinks = jump in 1st derivative
         integer,      intent(in)           :: naccel
         integer,      intent(in)           :: iaccel(naccel) ! accelerations = jump in 2nd derivative
         integer,      intent(out)          :: my_ierror
         real(kind=8), intent(in), optional :: wgt_arg(nmeas)
      end subroutine bspline_make_1d_ppspline
      
      !------------------------------------------------------------------------------------------------------

      module subroutine bspline_get_ppcoef_1seg_orig(nspl, tj, cj_y, jseg, c0, c1, c2, c3, idebug)
      !--function: determine PP-spline coefficients for one segment jseg
      !--subroutine arguments:
         integer,      intent(in)           :: nspl, jseg, idebug
         real(kind=8), intent(in)           :: cj_y(nspl), tj(nspl+4)
         real(kind=8), intent(out)          :: c0, c1, c2, c3
      end subroutine bspline_get_ppcoef_1seg_orig

      !------------------------------------------------------------------------------------------------------

      module subroutine bspline_get_ppcoef_1seg_modf(nspl, tj, cj_y, jseg, c0, c1, c2, c3, idebug)
      !--function: determine PP-spline coefficients for one segment jseg
      !--subroutine arguments:
         integer,      intent(in)           :: nspl, jseg, idebug
         real(kind=8), intent(in)           :: cj_y(nspl), tj(nspl+4)
         real(kind=8), intent(out)          :: c0, c1, c2, c3
      end subroutine bspline_get_ppcoef_1seg_modf

      !------------------------------------------------------------------------------------------------------

      module subroutine bspline_get_z_at_xy_list(spl2d, nout, xout, yout, zout, my_ierror, exterval)
      !--function: for nout half-parametric positions (xout,yout), determine vout in the 2D tensor B-spline
      !            and produce (zout)
      !--subroutine arguments:
         type(t_bspline2d)                    :: spl2d
         integer,                intent(in)   :: nout
         real(kind=8),           intent(in)   :: xout(nout), yout(nout)
         real(kind=8),           intent(out)  :: zout(nout)
         integer,                intent(out)  :: my_ierror
         real(kind=8), optional, intent(in)   :: exterval
      end subroutine bspline_get_z_at_xy_list
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_get_z_at_xy_prod(spl2d, nx, ny, xout, yout, zout, my_ierror, exterval)
      !--function: for nx x ny half-parametric positions (xout) x (yout), determine vij in the 2D tensor
      !            B-spline and produce (zout)
      !--subroutine arguments:
         type(t_bspline2d)                    :: spl2d
         integer,                intent(in)   :: nx, ny
         real(kind=8),           intent(in)   :: xout(nx), yout(ny)
         real(kind=8),           intent(out)  :: zout(nx,ny)
         integer,                intent(out)  :: my_ierror
         real(kind=8), optional, intent(in)   :: exterval
      end subroutine bspline_get_z_at_xy_prod
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_get_xz_at_uy_list(spl2d, nout, uout, yout, xout, zout, my_ierror, exterval)
      !--function: for nout full-parametric positions (uout,yout), determine vout in the 2D tensor B-spline
      !            and produce (xout,zout)
      !--subroutine arguments:
         type(t_bspline2d)                    :: spl2d
         integer,                intent(in)   :: nout
         real(kind=8),           intent(in)   :: uout(nout), yout(nout)
         real(kind=8),           intent(out)  :: xout(nout), zout(nout)
         integer,                intent(out)  :: my_ierror
         real(kind=8), optional, intent(in)   :: exterval
      end subroutine bspline_get_xz_at_uy_list
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine bspline_get_xz_at_uy_prod(spl2d, nu, ny, uout, yout, xout, zout, my_ierror,     &
                                exterval)
      !--function: for nu x ny full-parametric positions (uout) x (yout), determine vij in the 2D tensor
      !            B-spline and produce (xout,zout)
      !--subroutine arguments:
         type(t_bspline2d)                    :: spl2d
         integer,                intent(in)   :: nu, ny
         real(kind=8),           intent(in)   :: uout(nu), yout(ny)
         real(kind=8),           intent(out)  :: xout(nu,ny), zout(nu,ny)
         integer,                intent(out)  :: my_ierror
         real(kind=8), optional, intent(in)   :: exterval
      end subroutine bspline_get_xz_at_uy_prod
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine test_bspline(my_ierror)
      !--function: compute least squares smoothing B-spline for flat+quadratic example
         integer,      intent(out)          :: my_ierror
      end subroutine test_bspline

   end interface

!------------------------------------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------------------------------------

end module m_bspline
