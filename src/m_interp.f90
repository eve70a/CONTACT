!------------------------------------------------------------------------------------------------------------
! m_interp - linear 1D interpolation, bilinear 2D, bicubic 2D interpolation
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_interp
   use m_globals
   use m_markers
   implicit none
   private

   ! Debugging for module m_interp

   integer  :: ldebug    =  0    ! local level of debugging
   integer  :: ii_debug  = -1    ! output point for which detailed info is requested (-1 = none)
   integer  :: iel_debug = -1    ! input element for which detailed info is requested (-1 = none)

   public  interp1d_set_debug

   ! Generic helper routines

   public  is_strictly_monotonic
   public  check_monotone
   public  locate_segment
   public  locate_interval

   ! Interpolation routines defined on 1d input data:

   public  interp_1d
   public  interp_1d_to_2d
   public  interp_1d_to_1d
   public  interp_1d_to_scalar

   interface interp_1d
      module procedure interp_1d_to_scalar
      module procedure interp_1d_to_1d
      module procedure interp_1d_to_2d
   end interface interp_1d

   ! Interpolation routines defined on 2d logically rectangular data:

   public  interp_wgt_surf2unif
   private compute_logical_uv
   public  interp_aply_surf2unif_1d
   public  interp_aply_surf2unif_1d_bicubic
   public  interp_aply_surf2unif_3d

   interface

      module subroutine interp1d_set_debug(new_ldebug, new_ii_debug, new_iel_debug)
      !--function: enable/disable debug output of interpolation routines
      !--subroutine arguments:
         integer, intent(in)           :: new_ldebug       ! level of debug output required
         integer, intent(in), optional :: new_ii_debug     ! specific point of interest for debugging
         integer, intent(in), optional :: new_iel_debug    ! specific point of interest for debugging
      end subroutine interp1d_set_debug

      !------------------------------------------------------------------------------------------------------

      module function is_strictly_monotonic( nin, xin )
      !--function: determine whether array values are strictly increasing or decreasing
      !--result value
         logical  :: is_strictly_monotonic
      !--subroutine arguments:
         integer,      intent(in)  :: nin                   ! number of points
         real(kind=8), intent(in)  :: xin(nin)              ! input array
      end function is_strictly_monotonic

      !------------------------------------------------------------------------------------------------------

      module subroutine check_monotone(npnt, xpnt, istrict, ierror)
      !--function: check whether array xpnt is monotonically increasing or decreasing
      !--subroutine arguments:
         integer,      intent(in)   :: npnt, istrict
         real(kind=8), intent(in)   :: xpnt(npnt)
         integer,      intent(out)  :: ierror
      end subroutine check_monotone

      !------------------------------------------------------------------------------------------------------

      module subroutine locate_segment( np, vp, vi, iseg )
      !--function: find segment 'iseg' in list of segments { vp(ip), ip=1,np } that contains position vi,
      !            starting with input iseg as initial estimate. 
      !            On output: iseg in { 0, .., np },  iseg=0 when vi < min(vp), iseg=np when vi >= max(vp).
      !            cf. hunt-algorithm of Numerical Recipes, sec. 3.4.
      !--subroutine arguments
         integer,      intent(in)    :: np
         real(kind=8), intent(in)    :: vp(np), vi
         integer,      intent(inout) :: iseg
      end subroutine locate_segment
      
      !------------------------------------------------------------------------------------------------------

      module subroutine locate_interval(np, vp, v_low, v_hig, ilow, ihig)
      !--purpose: determine ilow, ihig such that vp( [ilow:ihig] ) just encompasses [v_low, v_hig]
      !           return ilow=1 in case v_low < vp(1) and ihig=np in case v_hig > vp(np)
      !--subroutine parameters:
         integer,      intent(in)  :: np
         real(kind=8), intent(in)  :: vp(np), v_low, v_hig
         integer,      intent(out) :: ilow, ihig
      end subroutine locate_interval

      !------------------------------------------------------------------------------------------------------

      module subroutine interp_1d_to_2d( nin, xin, fin, nout1, nout2, xintp, fintp, ierror, exterval)
      !--purpose: interpolate (non equidistant) 1d input data [xin, fin] to the points listed in xintp,
      !           which can be 2d. Using linear interpolation, constant extrapolation of first/last values.
      !--subroutine arguments:
         integer,      intent(in)  :: nin, nout1, nout2     ! number of points in, out
         real(kind=8), intent(in)  :: xin(nin), fin(nin)    ! input data [x,fx]
         real(kind=8), intent(in)  :: xintp(nout1, nout2)   ! output locations [x]
         real(kind=8), intent(out) :: fintp(nout1, nout2)   ! output values [fx]
         integer,      intent(out) :: ierror
         real(kind=8), optional    :: exterval
      end subroutine interp_1d_to_2d

      !------------------------------------------------------------------------------------------------------

      module subroutine interp_1d_to_1d( nin, xin, fin, nout1, xintp, fintp, my_ierror, exterval)
      !--purpose: interpolate (non equidistant) 1d input data [xin, fin] to the points listed in xintp(:).
      !           Using linear interpolation, constant extrapolation of first/last values.
      !--subroutine arguments:
         integer,      intent(in)  :: nin, nout1            ! number of points in, out
         real(kind=8), intent(in)  :: xin(nin), fin(nin)    ! input data [x,fx]
         real(kind=8), intent(in)  :: xintp(nout1)          ! output locations [x]
         real(kind=8), intent(out) :: fintp(nout1)          ! output values [fx]
         integer,      intent(out) :: my_ierror
         real(kind=8), optional    :: exterval
      end subroutine interp_1d_to_1d

      !------------------------------------------------------------------------------------------------------

      module subroutine interp_1d_to_scalar( nin, xin, fin, xintp, fintp, my_ierror, exterval)
      !--purpose: interpolate (non equidistant) 1d input data [xin, fin] to the point xintp.
      !           Using linear interpolation, constant extrapolation of first/last values.
      !--subroutine arguments:
         integer,      intent(in)  :: nin                   ! number of points in
         real(kind=8), intent(in)  :: xin(nin), fin(nin)    ! input data [x,fx]
         real(kind=8), intent(in)  :: xintp                 ! output location [x]
         real(kind=8), intent(out) :: fintp                 ! output value [fx]
         integer,      intent(out) :: my_ierror
         real(kind=8), optional    :: exterval
      end subroutine interp_1d_to_scalar

      !------------------------------------------------------------------------------------------------------

      module subroutine interp_wgt_surf2unif(nnode_x, nnode_y, nnode, x_node, y_node, z_node, z_thrs,   &
                                      mx, my, nout, dx_arg, dy_arg, x_out, y_out, ii2iel, ii2nod,       &
                                      wii2nod, fac_uv, my_ierror)
      !--function: construct interpolation tables between a curvi-linear input grid (surface) and a uniform
      !            output-grid.
      !            Returns weighting factors for bilinear and relative positions for bicubic interpolation
      !--subroutine arguments:
         integer,      intent(in)  :: nnode_x, nnode_y, nnode ! number of points in the input grid
         real(kind=8), intent(in)  :: x_node(nnode), y_node(nnode), z_node(nnode)
                                                        ! (x,y,z) coordinates of input grid points
         real(kind=8), intent(in)  :: z_thrs            ! typical dimension in z-direction, threshold for bounding box
         integer,      intent(in)  :: mx, my, nout      ! number of points in the output grid
         real(kind=8), intent(in)  :: dx_arg, dy_arg    ! step-sizes used in the output grid, <0 for non-uniform
         real(kind=8), intent(in)  :: x_out(nout), y_out(nout) ! coordinates of the output grid points
         integer,      intent(out) :: ii2iel(nout)      ! input element number for output grid points
         integer,      intent(out) :: ii2nod(4,nout)    ! input node numbers for output grid points
         real(kind=8), intent(out) :: wii2nod(4,nout)   ! interpolation weights per input node per output grid point
         real(kind=8), intent(out) :: fac_uv(2,nout)    ! relative position of output grid point in input element
                                                        ! used for bicubic interpolation
         integer,      intent(out) :: my_ierror
      end subroutine interp_wgt_surf2unif

      !------------------------------------------------------------------------------------------------------

      module subroutine compute_logical_uv(nnode, x_node, y_node, i_ll, i_lr, i_ur, i_ul, x_out, y_out, &
                                    ii, u_out, v_out, loc_debug, ierror)
      !--function: compute relative [u,v] position within quadrilateral element for given physical [x,y]
      !            coordinates for use in bilinear interpolation
      !--subroutine arguments:
         integer,      intent(in)  :: nnode                         ! number of points in the input grid
         real(kind=8), intent(in)  :: x_node(nnode), y_node(nnode)  ! input grid [x,y] coordinates
         integer,      intent(in)  :: i_ll, i_lr, i_ur, i_ul        ! corners of quadrilateral A-B-C-D
         real(kind=8), intent(in)  :: x_out, y_out                  ! physical position of output point P
         integer,      intent(in)  :: ii                            ! point number - for print output only
         real(kind=8), intent(out) :: u_out, v_out                  ! logical position of output point P
         integer,      intent(in)  :: loc_debug                     ! level of debug output for current point
         integer,      intent(out) :: ierror
      end subroutine compute_logical_uv

      !------------------------------------------------------------------------------------------------------

      module subroutine interp_aply_surf2unif_1d(nnode, arr_node, nout, arr_out, ii2iel, ii2nod,        &
                                wii2nod, ierror, defval)
      !--function: interpolate scalar (1d) data given at the nodes of an input surface to the points of the
      !            output grid
      !--subroutine arguments:
         integer,      intent(in)  :: nnode             ! number of nodes in the input surface
         real(kind=8), intent(in)  :: arr_node(nnode)   ! function values of surface nodes 
                                                        ! (nval=3: associated with (x,y,z) dirs)
         integer,      intent(in)  :: nout              ! number of points in the output grid
         integer,      intent(in)  :: ii2iel(nout)      ! surface element number for output grid points
         integer,      intent(in)  :: ii2nod(4,nout)    ! surface node numbers for output grid points
         real(kind=8), intent(in)  :: wii2nod(4,nout)   ! interpolation weights per surface node per output point
         real(kind=8), intent(out) :: arr_out(nout)     ! function values interpolated to the output grid
         integer,      intent(out) :: ierror
         real(kind=8), optional    :: defval            ! default value
      end subroutine interp_aply_surf2unif_1d

      !------------------------------------------------------------------------------------------------------

      module subroutine interp_aply_surf2unif_1d_bicubic(nnode_x, nnode_y, nnode, f_node, nout, ii2nod, &
                                                  fac_uv, f_out, defval)
      !--function: perform bicubic interpolation based on relative positions
      !--subroutine arguments:
         integer,      intent(in)  :: nnode_x, nnode_y, nnode ! number of points in the input grid
         real(kind=8), intent(in)  :: f_node(nnode)           ! function values of input grid points
         integer,      intent(in)  :: nout                    ! number of points in the output grid
         integer,      intent(in)  :: ii2nod(4,nout)          ! input node numbers for output grid points
         real(kind=8), intent(in)  :: fac_uv(2,nout)          ! relative pos. of output points in input elements
         real(kind=8), intent(out) :: f_out(nout)             ! interpolated function values
         real(kind=8), optional    :: defval                  ! default value
      end subroutine interp_aply_surf2unif_1d_bicubic

      !------------------------------------------------------------------------------------------------------

      module subroutine interp_aply_surf2unif_3d(nnode, arr_node, nout, arr_out, ii2iel, ii2nod,        &
                        wii2nod, ierror, defval)
      !--function: interpolate 3d data given at the nodes of an input surface to the points of the 
      !            output grid
      !--subroutine arguments:
         integer,      intent(in)  :: nnode             ! number of nodes in the input surface
         real(kind=8), intent(in)  :: arr_node(nnode,3) ! function values of surface nodes 
                                                        ! (nval=3: associated with (x,y,z) dirs)
         integer,      intent(in)  :: nout              ! number of points in the output grid
         integer,      intent(in)  :: ii2iel(nout)      ! surface element number for output grid points
         integer,      intent(in)  :: ii2nod(4,nout)    ! surface node numbers for output grid points
         real(kind=8), intent(in)  :: wii2nod(4,nout)   ! interpolation weights per surface node per output point
         real(kind=8)              :: arr_out(nout,3)   ! function values interpolated to the output grid
         integer,      intent(out) :: ierror
         real(kind=8), optional    :: defval            ! default value
      end subroutine interp_aply_surf2unif_3d

      !------------------------------------------------------------------------------------------------------

      module subroutine interp_test
      !--function: test interpolation routines
         implicit none
      end subroutine interp_test

      !------------------------------------------------------------------------------------------------------

   end interface

!------------------------------------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------------------------------------

end module m_interp
