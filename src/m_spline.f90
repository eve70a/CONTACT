!------------------------------------------------------------------------------------------------------------
! m_spline - interface definitions for 1D spline in PP-form (piecewise polynomial)
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_spline
   use m_globals
   use m_markers
   use m_ptrarray
   use m_interp
   implicit none
   private

   ! Debugging for modules m_spline_def and m_spline_make

   public  splinedef_set_debug

   integer  :: ldebug    =  0    ! local level of debugging
   integer  :: ii_debug  = -1    ! output point for which detailed info is requested (-1 = none)
   integer  :: iel_debug = -1    ! input element for which detailed info is requested (-1 = none)

   ! Data type for splines:

   public t_spline

   !---------------------------------------------------------------------------------------------------------
   ! m_spline_def: Basic functions for splines:
   !---------------------------------------------------------------------------------------------------------

   public  spline_nullify
   public  spline_allocate
   public  spline_initx
   public  spline_copy
   public  spline_trim
   public  spline_print
   public  spline_destroy

   public  spline_eval
   public  spline_eval_arr
   public  spline_eval_spl
   public  spline_eval_scalar

   interface spline_eval
      module procedure spline_eval_arr
      module procedure spline_eval_spl
      module procedure spline_eval_scalar
   end interface spline_eval

   public  solve_cubic_segm
   public  solve_cubic_newton
   public  solve_cubic_cardano
   public  spline_get_s_at_f
   public  spline_get_s_at_f_spl

   interface spline_get_s_at_f
      module procedure spline_get_s_at_f_spl
   end interface spline_get_s_at_f

   public  locate_one_extremum
   public  locate_extremal_values

   public spline_shift
   public spline_mirror_y
   public spline_rotate
   public spline_roll
   public spline_2glob
   public spline_2glob_or
   public spline_2glob_m
   public spline_2loc
   public spline_2loc_m

   interface spline_2glob
      module procedure spline_2glob_m
      module procedure spline_2glob_or
   end interface spline_2glob

   interface spline_2loc
      module procedure spline_2loc_m
   end interface spline_2loc

   !---------------------------------------------------------------------------------------------------------
   ! m_spline_make: Functions for the creation of 1D PP-splines:
   !---------------------------------------------------------------------------------------------------------

   public  spline_check_updates

   public  spline_add_topview

   private ppspline_make_simple_sec_intpol
   private ppspline_make_simple_sec_smoothing
   public  ppspline_make_simple_kink

   public  ppspline_make_spline
   public  ppspline_make_spline_kink
   public  ppspline_make_spline_nokink

   interface ppspline_make_spline
      module procedure ppspline_make_spline_kink
      module procedure ppspline_make_spline_nokink
   end interface ppspline_make_spline

!------------------------------------------------------------------------------------------------------------
!  data for a parametric spline (s, [x(s)], y(s), z(s)):

   type :: t_spline
      integer      :: npnt
      logical      :: has_xdata

      real(kind=8), dimension(:),   pointer  :: s     => NULL() ! (npnt)

      real(kind=8), dimension(:,:), pointer  :: axspl => NULL() ! (npnt,4)
      real(kind=8), dimension(:),   pointer  :: ax0   => NULL()
      real(kind=8), dimension(:),   pointer  :: ax1   => NULL()
      real(kind=8), dimension(:),   pointer  :: ax2   => NULL()
      real(kind=8), dimension(:),   pointer  :: ax3   => NULL()

      real(kind=8), dimension(:,:), pointer  :: ayspl => NULL() ! (npnt,4)
      real(kind=8), dimension(:),   pointer  :: ay0   => NULL()
      real(kind=8), dimension(:),   pointer  :: ay1   => NULL()
      real(kind=8), dimension(:),   pointer  :: ay2   => NULL()
      real(kind=8), dimension(:),   pointer  :: ay3   => NULL()

      real(kind=8), dimension(:,:), pointer  :: azspl => NULL() ! (npnt,4)
      real(kind=8), dimension(:),   pointer  :: az0   => NULL()
      real(kind=8), dimension(:),   pointer  :: az1   => NULL()
      real(kind=8), dimension(:),   pointer  :: az2   => NULL()
      real(kind=8), dimension(:),   pointer  :: az3   => NULL()

      integer      :: nsec_uniy
      integer,      dimension(:),   pointer  :: ipnt_uniy => NULL() ! (nsec_uni+1)

      integer      :: nsec_top
      real(kind=8), dimension(:),   pointer  :: ysec_top  => NULL() ! (nsec_top+1)
      integer,      dimension(:),   pointer  :: iuni_top  => NULL() ! (nsec_top)

      ! npnt            number of points in spline
      ! has_xdata       indicating that axspl is present, esp. for contact locus
      !
      ! s       [mm]    arc length parameter
      ! axspl, ax0-ax3  simple spline coefficients for x(s), with ax0-3 pointers to columns of axspl
      ! ayspl, ay0-ay3  simple spline coefficients for y(s), with ay0-3 pointers to columns of ayspl
      ! azspl, az0-az3  simple spline coefficients for z(s), with az0-3 pointers to columns of azspl
      !
      ! nsec_uniy       number of sections with monotonous y(s), invertible, uni-valued s(y)
      ! ipnt_uniy       start-index of each section with uni-valued y(s)
      !
      ! nsec_top        number of sections [y_j, y_{j+1}] used in 'top view'
      ! ysec_top        start-position y_j of each section used in top view
      ! iuni_top        iuni_top(j) = index of uni-valued section seen in top view section j
      !
      ! Note: uni-sections and top-sections are volatile -- discarded upon shift/rotate/trimming

   end type t_spline

!------------------------------------------------------------------------------------------------------------

   interface

      !------------------------------------------------------------------------------------------------------
      ! m_spline_def: Basic functions for splines:
      !------------------------------------------------------------------------------------------------------

      module subroutine splinedef_set_debug(new_ldebug, new_ii_debug, new_iel_debug)
      !--function: enable/disable debug output of spline routines
      !--subroutine arguments:
         integer, intent(in)           :: new_ldebug       ! level of debug output required
         integer, intent(in), optional :: new_ii_debug     ! specific point of interest for debugging
         integer, intent(in), optional :: new_iel_debug    ! specific point of interest for debugging
      end subroutine splinedef_set_debug
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_nullify(spl)
      !--purpose: initialize spline structure, nullify pointers
      !--subroutine parameters:
         type(t_spline)            :: spl
      end subroutine spline_nullify
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_allocate(spl, npnt, nsec_uniy, nsec_top)
      !--purpose: (re-)allocate the arrays needed for spline data
      !--subroutine parameters:
         type(t_spline)       :: spl
         integer              :: npnt
         integer, optional    :: nsec_uniy, nsec_top
      end subroutine spline_allocate
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_initx(spl)
      !--purpose: initialize the arrays needed for x-spline data, set x==0
      !--subroutine parameters:
         type(t_spline)            :: spl
      end subroutine spline_initx
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_copy(spl_in, spl_out)
      !--purpose: copy input spline spl_in to output spline spl_out
      !--subroutine parameters:
         type(t_spline)        :: spl_in, spl_out
      end subroutine spline_copy
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_trim(spl_in, spl_out, ilow_arg, ihig_arg, s_low, s_hig)
      !--purpose: create trimmed version of spline in output grid gout
      !--subroutine parameters:
         type(t_spline)         :: spl_in, spl_out
         integer,      optional :: ilow_arg, ihig_arg
         real(kind=8), optional :: s_low, s_hig
      end subroutine spline_trim
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_print(spl, nam, idebug, ndigit)
      !--function: print information on spline spl
      !--subroutine arguments
         type(t_spline)       :: spl
         character(len=*)     :: nam
         integer              :: idebug
         integer, optional    :: ndigit       ! number of significant digits
      end subroutine spline_print
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_destroy(spl)
      !--purpose: clean-up allocated arrays for splines, nullify pointers
      !--subroutine parameters:
         type(t_spline)  :: spl
      end subroutine spline_destroy
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_eval_arr(npnt, s_spl, nam_f, a3, a2, a1, a0, nout, s_eval, ierror,       &
                                  exterval, f_eval, f1_eval, f2_eval, f3_eval)
      !--function: evaluate spline {s, a,b,c,d} at points s_eval, computing f_eval and its derivatives
      !             f1=f', f2=f'', f3=f'''
      !--subroutine arguments:
         integer,          intent(in)   :: npnt, nout
         character(len=1), intent(in)   :: nam_f
         real(kind=8),     intent(in)   :: s_spl(npnt), a3(npnt), a2(npnt), a1(npnt), a0(npnt), s_eval(nout)
         integer,          intent(out)  :: ierror
         real(kind=8),     intent(in),  optional :: exterval
         real(kind=8),     intent(out), optional :: f_eval(nout), f1_eval(nout), f2_eval(nout), f3_eval(nout)
      end subroutine spline_eval_arr
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_eval_spl(spl, ikarg, nout, s_eval, ierror, exterval, f_eval, f1_eval,    &
                                f2_eval, f3_eval)
      !--function: evaluate parametric spline 'spl' in direction ikarg at points s_eval, computing f_eval
      !            and its derivatives f1=f', f2=f'', f3=f'''
      !--subroutine arguments:
         type(t_spline)                 :: spl
         integer,          intent(in)   :: ikarg, nout
         real(kind=8),     intent(in)   :: s_eval(nout)
         integer,          intent(out)  :: ierror
         real(kind=8),     intent(in),  optional :: exterval
         real(kind=8),     intent(out), optional :: f_eval(nout), f1_eval(nout), f2_eval(nout), f3_eval(nout)
      end subroutine spline_eval_spl
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_eval_scalar(spl, ikarg, s_eval, ierror, exterval, f_eval, f1_eval,       &
                                f2_eval, f3_eval)
      !--function: evaluate parametric spline 'spl' in direction ikarg at point s_eval, computing f_eval
      !            and its derivatives f1=f', f2=f'', f3=f'''
      !--subroutine arguments:
         type(t_spline)                 :: spl
         integer,          intent(in)   :: ikarg
         real(kind=8),     intent(in)   :: s_eval
         integer,          intent(out)  :: ierror
         real(kind=8),     intent(in),  optional :: exterval
         real(kind=8),     intent(out), optional :: f_eval, f1_eval, f2_eval, f3_eval
      end subroutine spline_eval_scalar
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine solve_cubic_segm( s_a, s_b, f_a, c1, c2, c3, f_b, fval, sl, ldebug, ierror)
      !--function: solve f(s) = c3 sl^3 + c2 sl^2 + c1 sl + f_a - fval = 0 ,
      !            with sl the local coordinate in the segment.
      !--subroutine arguments
         integer,      intent(in)  :: ldebug
         real(kind=8), intent(in)  :: s_a, s_b, f_a, f_b, c1, c2, c3, fval
         real(kind=8), intent(out) :: sl
         integer,      intent(out) :: ierror
      end subroutine solve_cubic_segm
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine solve_cubic_cardano( s_a, s_b, f_a, c1, c2, c3, fval, sl, ldebug)
      !--function: solve f(s) = c3 sl^3 + c2 sl^2 + c1 sl + f_a - fval = 0 , using Cardano formulas,
      !            with sl the local coordinate in the segment.
      ! solve 3rd degree equation f(s) = y, returning one root in interval [s_a,s_b]
      !--subroutine arguments
         integer,      intent(in)  :: ldebug
         real(kind=8), intent(in)  :: s_a, s_b, f_a, c1, c2, c3, fval
         real(kind=8), intent(out) :: sl
      end subroutine solve_cubic_cardano
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine solve_cubic_newton( s_a, s_b, f_a, c1, c2, c3, f_b, fval, sl, ldebug, ierror)
      !--function: solve f(s) = c3 sl^3 + c2 sl^2 + c1 sl + f_a - fval = 0 using Newton-Raphson,
      !            with sl the local coordinate in the segment.
      !--subroutine arguments
         integer,      intent(in)  :: ldebug
         real(kind=8), intent(in)  :: s_a, s_b, f_a, f_b, c1, c2, c3, fval
         real(kind=8), intent(out) :: sl
         integer,      intent(out) :: ierror
      end subroutine solve_cubic_newton
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_get_s_at_f_spl( spl, ikarg, nout, f_out, s_out, my_ierror )
      !--function: interpolate parametric spline 'spl' in direction ikarg, determine s_out at the requested
      !            f_out-positions
      !      Note: the s-values are not uniquely defined in case the inputs are non-monotonic,
      !            e.g. back-of-flange, y_w oscillating around 0.
      !            A unique result is obtained using the 'top view administration', if available.
      !            Else, s := s(1) for all f <= f(1) and s := s(np) for f >= f(np), this is a conservative
      !            approach for the range of s-values needed to cover [f0,f1].
      !--subroutine arguments
         type(t_spline)                :: spl
         integer,          intent(in)  :: ikarg
         integer,          intent(in)  :: nout
         real(kind=8),     intent(in)  :: f_out(nout)
         real(kind=8),     intent(out) :: s_out(nout)
         integer,          intent(out) :: my_ierror
      end subroutine spline_get_s_at_f_spl
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine locate_one_extremum(npnt, s_spl, a3, a2, a1, a0, iseg, typ_xtrm, s_xtrm, ierror)
      !--function: for a simple spline, determine the s-position of a locally extremal function value
      !            --> for segment iseg with f'(s0)*f'(s1)<0, solve f'(s) = 0 with quadratic f'
      !--subroutine arguments
         integer,          intent(in)   :: npnt, iseg
         integer,          intent(in)   :: typ_xtrm   ! <0: minimum, 0: dont care, >0: maximum
         real(kind=8),     intent(in)   :: s_spl(npnt), a3(npnt), a2(npnt), a1(npnt), a0(npnt)
         real(kind=8),     intent(out)  :: s_xtrm
         integer,          intent(out)  :: ierror
      end subroutine locate_one_extremum
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine locate_extremal_values(npnt, s_spl, a3, a2, a1, a0, nxtrm, s_xtrm, my_ierror)
      !--function: for a simple spline, determine the s-position with locally extremal function values
      !            --> for all segments segments with f'(s0)*f'(s1)<0, solve f'(s) = 0 with quadratic f'
      !            NOTE: this subroutine overlaps with spline_get_s_at_minval_arr in m_spline_get.f90.
      !--subroutine arguments
         integer,          intent(in)   :: npnt
         real(kind=8),     intent(in)   :: s_spl(npnt), a3(npnt), a2(npnt), a1(npnt), a0(npnt)
         real(kind=8),     intent(out)  :: s_xtrm(npnt)
         integer,          intent(out)  :: nxtrm, my_ierror
      end subroutine locate_extremal_values
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_shift(spl, dx, dy, dz)
      !--function: shift a spline in cartesian coordinates: spl = spl + [dx; dy; dz]
      !--subroutine arguments
         type(t_spline)  :: spl
         real(kind=8)    :: dx, dy, dz
      end subroutine spline_shift
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_mirror_y(spl)
      !--function: mirror a spline in cartesian coordinates wrt plane Oxz: spl.y = -spl.y
      !--subroutine arguments
         type(t_spline)  :: spl
      end subroutine spline_mirror_y
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_rotate(spl, rot, xc, yc, zc)
      !--function: rotate a spline in cartesian coordinates by rotation matrix about center-point [xc;yc;zc]
      !--subroutine arguments
         type(t_spline)           :: spl          ! input/output spline
         type(t_rotmat)           :: rot          ! matrix with new orientation of original unit vectors
         real(kind=8), intent(in) :: xc, yc, zc   ! rotation origin
      end subroutine spline_rotate
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_roll(spl, roll, yc, zc)
      !--function: rotate a spline in cartesian coordinates by roll angle roll [rad] (about x-axis/
      !            point [yc;zc])
      !--subroutine arguments
         type(t_spline)           :: spl
         real(kind=8), intent(in) :: roll     ! rotation angle [rad]
         real(kind=8), intent(in) :: yc, zc   ! rotation origin
      end subroutine spline_roll
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_2glob_or(spl, o, R)
      !--function: compute local-to-global conversion for a spline defined with respect to local system (o, R)
      !            o = origin of local system w.r.t. global system
      !            R = orientation of local system w.r.t. global system
      !            after the transformation, spl is defined with respect to the global system
      !--subroutine arguments
         type(t_spline), intent(inout) :: spl
         type(t_vec),    intent(in)    :: o
         type(t_rotmat), intent(in)    :: R
      end subroutine spline_2glob_or
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_2glob_m(spl, mref)
      !--function: compute local-to-global conversion for a spline defined with respect to mref
      !            o = mref%o   == origin of reference w.r.t. global system
      !            R = mref%rot == orientation of reference w.r.t. global system
      !            after the transformation, spl is defined with respect to the global system
      !--subroutine arguments
         type(t_marker), intent(in)    :: mref
         type(t_spline), intent(inout) :: spl
      end subroutine spline_2glob_m
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_2loc_m(spl, mref)
      !--function: convert spline spl defined in global coordinates to local coordinates according to mref
      !            o = mref%o   == origin of new reference w.r.t. global system
      !            R = mref%rot == orientation of new reference w.r.t. global system
      !            after the transformation, spl is defined with respect to the local system
      !--subroutine arguments
         type(t_marker), intent(in)    :: mref
         type(t_spline), intent(inout) :: spl
      end subroutine spline_2loc_m
      
      !------------------------------------------------------------------------------------------------------
      ! m_spline_make: Functions for the creation of 1D PP-splines:
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_add_breaks(spl, nbrk, sbrk, ipbrk, my_ierror)
      !--function: insert additional breaks at sbrk in PP-form of spline, except where sbrk coincides with
      !            existing breaks
      !--subroutine arguments:
         type(t_spline)                     :: spl
         integer,      intent(in)           :: nbrk
         real(kind=8), intent(in)           :: sbrk(nbrk)
         integer,      intent(out)          :: ipbrk(nbrk)
         integer,      intent(out)          :: my_ierror
      end subroutine spline_add_breaks
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_set_topview_sections(spl, view_minz, max_ntop, nsec_top, ysec_top,       &
                                iuni_top, my_ierror)
      !--function: determine the sections [ y_j, y_{j+1} ] of the top (rail) or bottom view (wheel) and
      !            set pointers to corresponding uni-valued sections in spline
      !--subroutine arguments:
         type(t_spline)                     :: spl
         logical,      intent(in)           :: view_minz   ! view surface from z=-inf (rail) or z=+inf (wheel)
         integer,      intent(in)           :: max_ntop
         integer,      intent(out)          :: nsec_top, my_ierror
         real(kind=8), intent(out)          :: ysec_top(max_ntop+1)
         integer,      intent(out)          :: iuni_top(max_ntop)
      end subroutine spline_set_topview_sections
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_add_topview(spl, view_minz, my_ierror)
      !--function: determine sections in spline with uni-valued y(s) and table { ybrk } with pointers to
      !            visible parts of the surface
      !--subroutine arguments:
         type(t_spline)                     :: spl
         logical,      intent(in)           :: view_minz   ! view surface from z=-inf (rail) or z=+inf (wheel)
         integer,      intent(out)          :: my_ierror
      end subroutine spline_add_topview
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine spline_check_updates(spl, nprf, s_prf, x_prf, y_prf, z_prf, k_chk, dist_max,    &
                                my_ierror)
      !--function: determine max distance between parametric spline (x(s),y(s),z(s)) and original data
      !--subroutine arguments:
         type(t_spline)                     :: spl
         integer,      intent(in)           :: nprf
         real(kind=8), intent(in)           :: s_prf(nprf), x_prf(nprf), y_prf(nprf), z_prf(nprf)
         integer,      intent(in)           :: k_chk      ! refinement factor for checking
         real(kind=8), intent(out)          :: dist_max
         integer,      intent(out)          :: my_ierror
      end subroutine spline_check_updates
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine ppspline_make_simple_sec_intpol(tot_pnt, ip0, ip1, s, y, a3, a2, a1, a0, ierror)
      !--function: compute section [ip0:ip1] of simple cubic interpolating spline {a,b,c,d} for given
      !            data {s,y} using "free end" boundary conditions / natural spline.
      !--subroutine arguments:
         integer,      intent(in)           :: tot_pnt   ! total number of points in data
         integer,      intent(in)           :: ip0, ip1  ! range of points for current section
         real(kind=8), intent(in)           :: s(tot_pnt), y(tot_pnt)
         real(kind=8), intent(out)          :: a3(tot_pnt), a2(tot_pnt), a1(tot_pnt), a0(tot_pnt)
         integer,      intent(out)          :: ierror
      end subroutine ppspline_make_simple_sec_intpol
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine ppspline_make_simple_sec_smoothing(tot_pnt, ip0, ip1, s, y, lambda, a3, a2, a1, &
                                a0, ierror, wgt)
      !--function: compute section [ip0:ip1] of simple cubic smoothing spline {a,b,c,d} for given data {s,y},
      !            weights wgt and parameter lambda, using "free end" boundary conditions / natural spline.
      !--subroutine arguments:
         integer,      intent(in)           :: tot_pnt   ! total number of points in data
         integer,      intent(in)           :: ip0, ip1  ! range of points for current section
         real(kind=8), intent(in)           :: s(tot_pnt), y(tot_pnt), lambda
         real(kind=8), intent(out)          :: a3(tot_pnt), a2(tot_pnt), a1(tot_pnt), a0(tot_pnt)
         integer,      intent(out)          :: ierror
         real(kind=8), intent(in), optional :: wgt(tot_pnt)
      end subroutine ppspline_make_simple_sec_smoothing
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine ppspline_make_simple_kink(npnt, s, y, lambda, a3, a2, a1, a0, nkink, ikinks,    &
                                my_ierror, wgt)
      !--function: compute simple cubic smoothing spline {a,b,c,d} for given data {s,y}, weights wgt and
      !            parameter lambda, using "free end" boundary conditions / natural spline.
      !--subroutine arguments:
         integer,      intent(in)           :: npnt, nkink
         real(kind=8), intent(in)           :: s(npnt), y(npnt), lambda
         integer,      intent(in)           :: ikinks(nkink)  ! including 1, npnt+1
         real(kind=8), intent(out)          :: a3(npnt), a2(npnt), a1(npnt), a0(npnt)
         integer,      intent(out)          :: my_ierror
         real(kind=8), intent(in), optional :: wgt(npnt)
      end subroutine ppspline_make_simple_kink
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine ppspline_make_spline_kink(spl, nmeas, s_prf, x_prf, y_prf, z_prf, has_xdata,    &
                                lambda, use_wgt, nkink, ikinks, my_ierror)
      !--function: compute parametric smoothing spline for 1-d grid (x(s),y(s),z(s))
      !--subroutine arguments:
         type(t_spline)                     :: spl
         integer                            :: nmeas
         real(kind=8), intent(inout)        :: s_prf(nmeas), x_prf(nmeas), y_prf(nmeas), z_prf(nmeas) 
         logical,      intent(in)           :: has_xdata  ! whether x_prf is a dummy or needs to be used
         real(kind=8), intent(in)           :: lambda     ! weight of 2nd derivative, relative to wgt_in
         logical,      intent(in)           :: use_wgt    ! weigh data with spacing ds
         integer,      intent(in)           :: nkink
         integer,      intent(in)           :: ikinks(nkink) ! start/end of spline sections
         integer,      intent(out)          :: my_ierror
      end subroutine ppspline_make_spline_kink
      
      !------------------------------------------------------------------------------------------------------
      
      module subroutine ppspline_make_spline_nokink(spl, nmeas, s_prf, x_prf, y_prf, z_prf, has_xdata,  &
                                lambda, use_wgt, my_ierror)
      !--function: compute parametric smoothing spline for 1-d grid (x(s),y(s),z(s))
      !--subroutine arguments:
         type(t_spline)                     :: spl
         integer                            :: nmeas
         real(kind=8), intent(inout)        :: s_prf(nmeas), x_prf(nmeas), y_prf(nmeas), z_prf(nmeas) 
         logical,      intent(in)           :: has_xdata  ! whether x_prf is a dummy or needs to be used
         real(kind=8), intent(in)           :: lambda     ! weight of 2nd derivative, using weight 1 for data
         logical,      intent(in)           :: use_wgt    ! weigh data with spacing ds
         integer,      intent(out)          :: my_ierror
      end subroutine ppspline_make_spline_nokink
      
      !------------------------------------------------------------------------------------------------------

   end interface

!------------------------------------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------------------------------------

end module m_spline
