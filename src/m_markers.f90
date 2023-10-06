!------------------------------------------------------------------------------------------------------------
! m_markers - data-structures for 'markers' defining local coordinate systems
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_markers
use m_globals
implicit none
private

public coordsys_create
public coordsys_id
public coordsys_name

public vec
public vec_zero
public vec_norm
public vec_dot
public vec_cross
public vec_2glob
public vec_2loc
public vec_veloc2glob
public vec_print

interface vec
   module procedure vec111
   module procedure vec3
end interface vec

interface vec_2glob
   module procedure vec_2glob_m
   module procedure vec_2glob_or
end interface vec_2glob

interface vec_2loc
   module procedure vec_2loc_m
   module procedure vec_2loc_or
end interface vec_2loc

interface vec_veloc2glob
   module procedure vec_veloc2glob_r_o
   module procedure vec_veloc2glob_r_m
   module procedure vec_veloc2glob_m_o
   module procedure vec_veloc2glob_m_m
end interface vec_veloc2glob

public rotmat
public rotmat_identity
public rotmat_pure_roll
public rotmat_pure_yaw
public rotmat_pure_pitch
public rotmat_roll_yaw
public rotmat_transpose
public rotmat_mirror_y
public rotmat_print

public marker
public marker_init
public marker_shift
public marker_roll
public marker_yaw
public marker_pitch
public marker_rotate
public marker_2glob
public marker_2loc
public marker_transpose
public marker_mirror_y
public marker_print

public markers_test

interface marker
   module procedure marker111
   module procedure marker3
   module procedure markerv
end interface marker

interface marker_rotate
   module procedure marker_rotate_eulang
   module procedure marker_rotate_rotmat
end interface marker_rotate

public operator (+)
interface operator (+)
   procedure vec_add
end interface operator (+)

public operator (-)
interface operator (-)
   procedure vec_neg
   procedure vec_diff
end interface operator (-)

public operator (*)
interface operator (*)
   procedure dble_times_vec
   procedure vec_times_dble
   procedure rotmat_vec
   procedure rotmat_product
end interface operator (*)

public operator (/)
interface operator (/)
   procedure vec_div_dble
end interface operator (/)

public operator (.dot.)
interface operator (.dot.)
   procedure vec_dot
end interface operator (.dot.)

public operator (.cross.)
interface operator (.cross.)
   procedure vec_cross
end interface operator (.cross.)

public operator (.transp.)
interface operator (.transp.)
   procedure rotmat_transp_vec
   procedure marker_transp_vec
   procedure rotmat_transp_product
end interface operator (.transp.)

!------------------------------------------------------------------------------------------------------------
!  data-type for a coordinate system

   integer, parameter             :: CSYS_LEN = 32

   type, public :: t_coordsys
      character(len=CSYS_LEN)     :: name
   end type t_coordsys

!------------------------------------------------------------------------------------------------------------
!  data-type for a point or vector in 3D space

   type, public :: t_vec
      real(kind=8) :: v(3)   ! (x,y,z) coordinates of the point/vector
   contains
      procedure :: x    => vec_get_x
      procedure :: y    => vec_get_y
      procedure :: z    => vec_get_z
      procedure :: norm => vec_norm
   end type t_vec

!------------------------------------------------------------------------------------------------------------
!  data-type for an orientation matrix in 3D space

   type, public :: t_rotmat
      real(kind=8) :: r(9)      ! 3x3 matrix, stored in (Fortran-like) column-major numbering
   contains
      procedure :: el     => rotmat_get_element
      procedure :: ivec   => rotmat_get_ivec
      procedure :: jvec   => rotmat_get_jvec
      procedure :: kvec   => rotmat_get_kvec
      procedure :: transp => rotmat_transpose
      procedure :: roll   => rotmat_get_rollangle
      procedure :: yaw    => rotmat_get_yawangle
      procedure :: pitch  => rotmat_get_pitchangle

      ! r          orientation of a marker with respect to a global system. This consists of the local
      !            vectors e1, e2, e3 expressed as i, j, k in terms of coordinates of the global system
      ! i,j,k      pointers to the local unit vectors in terms of global coordinates

   end type t_rotmat

!------------------------------------------------------------------------------------------------------------
!  data-type for a marker: a local coordinate system within a global system

   type, public :: t_marker
      type(t_vec)     :: o
      type(t_rotmat)  :: rot
      ! integer         :: csys
   contains
      procedure :: x      => marker_get_x
      procedure :: y      => marker_get_y
      procedure :: z      => marker_get_z
      procedure :: ivec   => marker_get_ivec
      procedure :: jvec   => marker_get_jvec
      procedure :: kvec   => marker_get_kvec
      procedure :: roll   => marker_get_rollangle
      procedure :: yaw    => marker_get_yawangle
      procedure :: pitch  => marker_get_pitchangle

      ! csys       integer handle to (type of) parent coordinate system: trk, ws, rw, rr, cp or unkn
      ! o          (x,y,z) coordinates of the local origin within the global system
      ! rot        orientation of the marker with respect to the global system. This consists of the local
      !            vectors e1, e2, e3 expressed as i, j, k in terms of coordinates of the global system

   end type t_marker

!------------------------------------------------------------------------------------------------------------

   integer, parameter      :: CSYS_UNKN    =  -1
   integer, parameter      :: CSYS_TRK     =   1
   integer, parameter      :: CSYS_WS      =   2
   integer, parameter      :: CSYS_LW      =   3
   integer, parameter      :: CSYS_RW      =   4
   integer, parameter      :: CSYS_LR      =   5
   integer, parameter      :: CSYS_RR      =   6
   integer, parameter      :: CSYS_CP      =   7

   integer, parameter      :: MAX_NUM_CSYS = 999      ! maximum number of coordinate systems
   integer                 :: num_csys = 0            ! number of coordinate-systems in use
   type(t_coordsys)        :: csys_list(MAX_NUM_CSYS) ! administration of coordinate systems
!  character(len=CSYS_LEN) :: csys_names(MAX_NUM_CSYS)

contains

!------------------------------------------------------------------------------------------------------------

function coordsys_create(name)
!--function: define a new coordinate system
   implicit none
!--result value
   integer :: coordsys_create
!--subroutine arguments
   character(len=*), intent(in) :: name

   ! TODO: check if name does not yet exist

   ! - increment number of coordinate systems defined
   num_csys = num_csys + 1
   ! - store information on new coordinate system
   csys_list(num_csys)%name = name
   ! - return the ID
   coordsys_create = num_csys

end function coordsys_create

!------------------------------------------------------------------------------------------------------------

function coordsys_id(name)
!--function: get ID of an existing coordinate system, or return -1 as error code
   implicit none
!--result value
   integer :: coordsys_id
!--subroutine arguments
   character(len=*), intent(in) :: name
!  local variables
   integer :: isys
   logical :: lfound

   isys = 0
   lfound = .false.
   do while(.not.lfound .and. isys.lt.num_csys)
      isys = isys + 1
      lfound = (trim(name).eq.trim(csys_list(isys)%name))
   enddo
   if (.not.lfound) isys = -1

   coordsys_id = isys

end function coordsys_id

!------------------------------------------------------------------------------------------------------------

function coordsys_name(isys)
!--function: get name of an existing coordinate system, or return ' ' as error code
   implicit none
!--result value
   character(len=CSYS_LEN) :: coordsys_name
!--subroutine arguments
   integer, intent(in)     :: isys

   ! TODO: check if name does not yet exist

   if (isys.ge.1 .and. isys.le.num_csys) then
      coordsys_name = csys_list(isys)%name
   else
      coordsys_name = ' '
   endif

end function coordsys_name

!------------------------------------------------------------------------------------------------------------

function vec_get_x(this)
!--function: get vector element v_x
   implicit none
!--result value
   real(kind=8)             :: vec_get_x
!--subroutine arguments
   class(t_vec), intent(in) :: this

   vec_get_x = this%v(1)
end function vec_get_x

!------------------------------------------------------------------------------------------------------------

function vec_get_y(this)
!--function: get vector element v_y
   implicit none
!--result value
   real(kind=8)             :: vec_get_y
!--subroutine arguments
   class(t_vec), intent(in) :: this

   vec_get_y = this%v(2)
end function vec_get_y

!------------------------------------------------------------------------------------------------------------

function vec_get_z(this)
!--function: get vector element v_z
   implicit none
!--result value
   real(kind=8)             :: vec_get_z
!--subroutine arguments
   class(t_vec), intent(in) :: this

   vec_get_z = this%v(3)
end function vec_get_z

!------------------------------------------------------------------------------------------------------------

function dble_times_vec(a, v)
!--function: compute scaled vector a * v
   implicit none
!--result value
   type(t_vec)                  :: dble_times_vec
!--subroutine arguments
   real(kind=8), intent(in) :: a
   type(t_vec),  intent(in) :: v

   dble_times_vec%v(1:3) = a * v%v(1:3)
end function dble_times_vec

!------------------------------------------------------------------------------------------------------------

function vec_times_dble(this, a)
!--function: compute scaled vector a * v
   implicit none
!--result value
   type(t_vec)              :: vec_times_dble
!--subroutine arguments
   class(t_vec), intent(in) :: this
   real(kind=8), intent(in) :: a

   vec_times_dble%v(1:3) = a * this%v(1:3)
end function vec_times_dble

!------------------------------------------------------------------------------------------------------------

function vec_div_dble(this, a)
!--function: compute scaled vector v / a
   implicit none
!--result value
   type(t_vec)              :: vec_div_dble
!--subroutine arguments
   class(t_vec), intent(in) :: this
   real(kind=8), intent(in) :: a

   vec_div_dble%v(1:3) = this%v(1:3) / a
end function vec_div_dble

!------------------------------------------------------------------------------------------------------------

function vec111(x, y, z)
!--function: create vector v from its elements
   implicit none
!--result value
   type(t_vec)              :: vec111
!--subroutine arguments
   real(kind=8), intent(in) :: x, y, z

   vec111%v(1) = x
   vec111%v(2) = y
   vec111%v(3) = z
end function vec111

!------------------------------------------------------------------------------------------------------------

function vec3(xyz)
!--function: create vector v from its elements
   implicit none
!--result value
   type(t_vec)              :: vec3
!--subroutine arguments
   real(kind=8), intent(in) :: xyz(3)

   vec3%v(1:3) = xyz(1:3)
end function vec3

!------------------------------------------------------------------------------------------------------------

function vec_zero()
!--function: return zero vector
   implicit none
!--result value
   type(t_vec)              :: vec_zero

   vec_zero%v = (/ 0d0, 0d0, 0d0 /)
end function vec_zero

!------------------------------------------------------------------------------------------------------------

function vec_neg(v)
!--function: compute negative of vector vout = -vin
   implicit none
!--result value
   type(t_vec)              :: vec_neg
!--subroutine arguments
   type(t_vec),  intent(in) :: v

   vec_neg%v(1:3) = - v%v(1:3)
end function vec_neg

!------------------------------------------------------------------------------------------------------------

function vec_add(v1, v2)
!--function: compute vector sum  v = v1 + v2
   implicit none
!--result value
   type(t_vec)              :: vec_add
!--subroutine arguments
   type(t_vec),  intent(in) :: v1, v2

   vec_add%v(1:3) = v1%v(1:3) + v2%v(1:3)
end function vec_add

!------------------------------------------------------------------------------------------------------------

function vec_diff(v1, v2)
!--function: compute vector sum  v = v1 - v2
   implicit none
!--result value
   type(t_vec)              :: vec_diff
!--subroutine arguments
   type(t_vec),  intent(in) :: v1, v2

   vec_diff%v(1:3) = v1%v(1:3) - v2%v(1:3)
end function vec_diff

!------------------------------------------------------------------------------------------------------------

function vec_norm(v)
!--function: compute length |v|
   implicit none
!--result value
   real(kind=8) :: vec_norm
!--subroutine arguments
   class(t_vec), intent(in) :: v

   vec_norm = sqrt(v%v(1)**2 + v%v(2)**2 + v%v(3)**2)
end function vec_norm

!------------------------------------------------------------------------------------------------------------

function vec_dot(v1, v2)
!--function: compute dot product (v1, v2)
   implicit none
!--result value
   real(kind=8) :: vec_dot
!--subroutine arguments
   type(t_vec), intent(in) :: v1, v2

   vec_dot = v1%v(1)*v2%v(1) + v1%v(2)*v2%v(2) + v1%v(3)*v2%v(3)
end function vec_dot

!------------------------------------------------------------------------------------------------------------

function vec_cross(v1, v2)
!--function: compute cross product v3 = v1 x v2
   implicit none
!--result value
   type(t_vec)              :: vec_cross
!--subroutine arguments
   type(t_vec), intent(in)  :: v1, v2

   vec_cross%v(1) = v1%v(2) * v2%v(3) - v1%v(3) * v2%v(2)
   vec_cross%v(2) = v1%v(3) * v2%v(1) - v1%v(1) * v2%v(3)
   vec_cross%v(3) = v1%v(1) * v2%v(2) - v1%v(2) * v2%v(1)
end function vec_cross

!------------------------------------------------------------------------------------------------------------

function vec_2glob_or(vl, o, R)
!--function: compute local-to-global conversion vg = o + R vl
!            o == origin of local system w.r.t. global reference
!            R == orientation of local system w.r.t. global reference
   implicit none
!--result value
   type(t_vec)                 :: vec_2glob_or
!--subroutine arguments
   type(t_vec),    intent(in)  :: vl, o
   type(t_rotmat), intent(in)  :: R

   vec_2glob_or%v(1:3) = o%v(1:3) + R%r(1:3) * vl%v(1) + R%r(4:6) * vl%v(2) + R%r(7:9) * vl%v(3)
end function vec_2glob_or

!------------------------------------------------------------------------------------------------------------

function vec_2glob_m(vl, m)
!--function: compute local-to-global conversion vg = o + R vl
!            o = m%o   == origin of local system w.r.t. global reference
!            R = m%rot == orientation of local system w.r.t. global reference
   implicit none
!--result value
   type(t_vec)                 :: vec_2glob_m
!--subroutine arguments
   type(t_vec),    intent(in)  :: vl
   type(t_marker), intent(in)  :: m

   vec_2glob_m%v(1:3) = m%o%v(1:3) + m%rot%r(1:3) * vl%v(1) + m%rot%r(4:6) * vl%v(2) + m%rot%r(7:9) * vl%v(3)
end function vec_2glob_m

!------------------------------------------------------------------------------------------------------------

function vec_2loc_or(vg, o, R)
!--function: compute global-to-local conversion vl = R^t ( vg - o )
!            o == origin of local system w.r.t. global reference
!            R == orientation of local system w.r.t. global reference
   implicit none
!--result value
   type(t_vec)                 :: vec_2loc_or
!--subroutine arguments
   type(t_vec),    intent(in)  :: vg, o
   type(t_rotmat), intent(in)  :: R
!--local variables
   type(t_vec)  :: tmp
   associate( vl => vec_2loc_or )

   tmp%v(1:3) = vg%v(1:3) - o%v(1:3)
   vl%v(1) = vec_dot(R%ivec(), tmp)
   vl%v(2) = vec_dot(R%jvec(), tmp)
   vl%v(3) = vec_dot(R%kvec(), tmp)
   end associate
end function vec_2loc_or

!------------------------------------------------------------------------------------------------------------

function vec_2loc_m(vg, m)
!--function: compute global-to-local conversion vl = R^t ( vg - o )
!            o = m%o   == origin of local system w.r.t. global reference
!            R = m%rot == orientation of local system w.r.t. global reference
   implicit none
!--result value
   type(t_vec)                 :: vec_2loc_m
!--subroutine arguments
   type(t_vec),    intent(in)  :: vg
   type(t_marker), intent(in)  :: m
!--local variables
   type(t_vec)  :: tmp

   tmp%v(1:3) = vg%v(1:3) - m%o%v(1:3)
   vec_2loc_m%v(1) = vec_dot(m%rot%ivec(), tmp)
   vec_2loc_m%v(2) = vec_dot(m%rot%jvec(), tmp)
   vec_2loc_m%v(3) = vec_dot(m%rot%kvec(), tmp)
end function vec_2loc_m

!------------------------------------------------------------------------------------------------------------

function vec_veloc2glob_r_o(tvel_o_g, rvel_o_g, rot_o, p_l, tvel_p_l)
!--function: compute global velocity of a point from local velocity & velocity of coordinate origin
!            cf. x = o + R * bar(x) -->   dot(x) = dot(o) + dot(R) * bar(x) + R * dot(bar(x))
   implicit none
!--subroutine arguments
   type(t_vec)          :: tvel_o_g  ! translational velocity of local origin O in global coordinates
   type(t_vec)          :: rvel_o_g  ! angular velocity vector of local origin O in global coordinates
   type(t_rotmat)       :: rot_o     ! orientation of local origin O in global coordinates
   type(t_vec)          :: p_l       ! position of P in local coordinates
   type(t_vec)          :: tvel_p_l  ! translational velocity of point P w.r.t. O in local coordinates
!--function result:
   type(t_vec)          :: vec_veloc2glob_r_o

   ! translation \dot{p} w.r.t. global origin has three components:
   !   1. translation of \dot{o} of local origin within global system
   !   2. \omega \times (R \bar{p}), i.e. rotation of local origin within global system
   !   3. R * \dot{\bar{p}}, translation of point within local system rotated to global system

   vec_veloc2glob_r_o = tvel_o_g + (rvel_o_g .cross. (rot_o * p_l)) + (rot_o * tvel_p_l)

end function vec_veloc2glob_r_o

!------------------------------------------------------------------------------------------------------------

function vec_veloc2glob_r_m(tvel_o_g, rvel_o_g, rot_o, m_p, tvel_p_l)
!--function: compute global velocity of a point from local velocity & velocity of coordinate origin
!            cf. x = o + R * bar(x) -->   dot(x) = dot(o) + dot(R) * bar(x) + R * dot(bar(x))
   implicit none
!--subroutine arguments
   type(t_vec)          :: tvel_o_g  ! translational velocity of local origin O in global coordinates
   type(t_vec)          :: rvel_o_g  ! angular velocity vector of local origin O in global coordinates
   type(t_rotmat)       :: rot_o     ! orientation of local origin O in global coordinates
   type(t_marker)       :: m_p       ! position of P in local coordinates
   type(t_vec)          :: tvel_p_l  ! translational velocity of point P w.r.t. O in local coordinates
!--function result:
   type(t_vec)          :: vec_veloc2glob_r_m

   ! translation \dot{p} w.r.t. global origin has three components:
   !   1. translation of \dot{o} of local origin within global system
   !   2. \omega \times (R \bar{p}), i.e. rotation of local origin within global system
   !   3. R * \dot{\bar{p}}, translation of point within local system rotated to global system

   vec_veloc2glob_r_m = tvel_o_g + (rvel_o_g .cross. (rot_o * m_p%o)) + (rot_o * tvel_p_l)

end function vec_veloc2glob_r_m

!------------------------------------------------------------------------------------------------------------

function vec_veloc2glob_m_o(tvel_o_g, rvel_o_g, m_l, p_l, tvel_p_l)
!--function: compute global velocity of a point from local velocity & velocity of coordinate origin
!            cf. x = o + R * bar(x) -->   dot(x) = dot(o) + dot(R) * bar(x) + R * dot(bar(x))
   implicit none
!--subroutine arguments
   type(t_vec)          :: tvel_o_g  ! translational velocity of local origin O in global coordinates
   type(t_vec)          :: rvel_o_g  ! angular velocity vector of local origin O in global coordinates
   type(t_marker)       :: m_l       ! marker for local origin O in global coordinates
   type(t_vec)          :: p_l       ! position of P in local coordinates
   type(t_vec)          :: tvel_p_l  ! translational velocity of point P w.r.t. O in local coordinates
!--function result:
   type(t_vec)          :: vec_veloc2glob_m_o

   ! translation \dot{p} w.r.t. global origin has three components:
   !   1. translation of \dot{o} of local origin within global system
   !   2. \omega \times (R \bar{p}), i.e. rotation of local origin within global system
   !   3. R * \dot{\bar{p}}, translation of point within local system rotated to global system

   vec_veloc2glob_m_o = tvel_o_g + (rvel_o_g .cross. (m_l%rot * p_l)) + (m_l%rot * tvel_p_l)

end function vec_veloc2glob_m_o

!------------------------------------------------------------------------------------------------------------

function vec_veloc2glob_m_m(tvel_o_g, rvel_o_g, m_l, m_p, tvel_p_l)
!--function: compute global velocity of a point from local velocity & velocity of coordinate origin
!            cf. x = o + R * bar(x) -->   dot(x) = dot(o) + dot(R) * bar(x) + R * dot(bar(x))
   implicit none
!--subroutine arguments
   type(t_vec)          :: tvel_o_g  ! translational velocity of local origin O in global coordinates
   type(t_vec)          :: rvel_o_g  ! angular velocity vector of local origin O in global coordinates
   type(t_marker)       :: m_l       ! marker for local origin O in global coordinates
   type(t_marker)       :: m_p       ! marker for point P in local coordinates
   type(t_vec)          :: tvel_p_l  ! translational velocity of point P w.r.t. O in local coordinates
!--function result:
   type(t_vec)          :: vec_veloc2glob_m_m

   ! translation \dot{p} w.r.t. global origin has three components:
   !   1. translation of \dot{o} of local origin within global system
   !   2. \omega \times (R \bar{p}), i.e. rotation of local origin within global system
   !   3. R * \dot{\bar{p}}, translation of point within local system rotated to global system

   vec_veloc2glob_m_m = tvel_o_g + (rvel_o_g .cross. (m_l%rot * m_p%o)) + (m_l%rot * tvel_p_l)

end function vec_veloc2glob_m_m

!------------------------------------------------------------------------------------------------------------

subroutine vec_print(v, nam, idebug, ndigit)
!--function: print information on a vec
   implicit none
!--subroutine arguments
   type(t_vec)          :: v
   character(len=*)     :: nam
   integer              :: idebug
   integer, optional    :: ndigit       ! number of significant digits
!--local variables
   integer              :: my_ndigit, my_len
   character(len=16)    :: strng(3)

   if (present(ndigit)) then
      my_ndigit = ndigit
   else
      my_ndigit = 4
   endif
   my_ndigit = max(2, min(8, my_ndigit))
   my_len    = 8 + my_ndigit

   strng(1) = fmt_gs(my_len, my_ndigit, v%v(1))
   strng(2) = fmt_gs(my_len, my_ndigit, v%v(2))
   strng(3) = fmt_gs(my_len, my_ndigit, v%v(3))

   if (idebug.eq.3) then
      ! RecurDyn format:
      write(bufout,*) trim(nam),' = [', v%v(1),',', v%v(2),',', v%v(3),']^T.'
      call write_log(1, bufout)
   elseif (idebug.ge.1) then
      write(bufout,*) trim(nam), ' = (', strng(1)(1:my_len),',', strng(2)(1:my_len),',',                &
                strng(3)(1:my_len),')^T'
      call write_log(1, bufout)
   endif

end subroutine vec_print

!------------------------------------------------------------------------------------------------------------

function rotmat_get_ivec(this)
!--function: get first column i of rotation matrix
   implicit none
!--result value
   type(t_vec)                 :: rotmat_get_ivec
!--subroutine arguments
   class(t_rotmat), intent(in) :: this

   rotmat_get_ivec = t_vec( this%r(1:3) )
end function rotmat_get_ivec

!------------------------------------------------------------------------------------------------------------

function rotmat_get_jvec(this)
!--function: get first column i of rotation matrix
   implicit none
!--result value
   type(t_vec)                 :: rotmat_get_jvec
!--subroutine arguments
   class(t_rotmat), intent(in) :: this

   rotmat_get_jvec = t_vec( this%r(4:6) )
end function rotmat_get_jvec

!------------------------------------------------------------------------------------------------------------

function rotmat_get_kvec(this)
!--function: get first column i of rotation matrix
   implicit none
!--result value
   type(t_vec)                 :: rotmat_get_kvec
!--subroutine arguments
   class(t_rotmat), intent(in) :: this

   rotmat_get_kvec = t_vec( this%r(7:9) )
end function rotmat_get_kvec

!------------------------------------------------------------------------------------------------------------

function rotmat_get_element(this, i, j)
!--function: get element i,j of rotation matrix
   implicit none
!--result value
   real(kind=8)                :: rotmat_get_element
!--subroutine arguments
   class(t_rotmat), intent(in) :: this
   integer,         intent(in) :: i, j
!--local variables
   integer                     :: ij

   ij = i + 3*(j-1)
   rotmat_get_element = this%r(ij)
end function rotmat_get_element

!------------------------------------------------------------------------------------------------------------

function rotmat(ivec, jvec, kvec)
!--function: create rotation matrix r from direction vectors i, j, k
   implicit none
!--result value
   type(t_rotmat)           :: rotmat
!--subroutine arguments
   type(t_vec)              :: ivec, jvec, kvec

   rotmat%r = (/ ivec%v(1:3), jvec%v(1:3), kvec%v(1:3) /)
end function rotmat

!------------------------------------------------------------------------------------------------------------

function rotmat_identity()
!--function: return rotation matrix R = I
   implicit none
!--result value
   type(t_rotmat) :: rotmat_identity

!  note: this code shows the transpose R^T, with rows [ i^T; j^T; k^T ]:

   rotmat_identity%r = (/ 1d0,  0d0, 0d0,     0d0,  1d0, 0d0,    0d0,  0d0, 1d0 /)

end function rotmat_identity

!------------------------------------------------------------------------------------------------------------

function rotmat_pure_roll(roll)
!--function: initialize rotation matrix at R = R_x(roll)
   implicit none
!--result value
   type(t_rotmat)              :: rotmat_pure_roll
!--subroutine arguments
   real(kind=8),   intent(in)  :: roll
!--local variables
   real(kind=8)                :: cs, sn

   cs = cos(roll)
   sn = sin(roll)

!  note: this code shows the transpose R^T, with rows [ i^T; j^T; k^T ]:

   rotmat_pure_roll%r(1:3) = (/  1d0,  0d0, 0d0 /)
   rotmat_pure_roll%r(4:6) = (/  0d0,   cs,  sn /)
   rotmat_pure_roll%r(7:9) = (/  0d0,  -sn,  cs /)

end function rotmat_pure_roll

!------------------------------------------------------------------------------------------------------------

function rotmat_pure_yaw(yaw)
!--function: initialize rotation matrix at R = R_z(yaw)
   implicit none
!--result value
   type(t_rotmat)              :: rotmat_pure_yaw
!--subroutine arguments
   real(kind=8),   intent(in)  :: yaw
!--local variables
   real(kind=8)                :: cs, sn

   cs = cos(yaw)
   sn = sin(yaw)

!  note: this code shows the transpose R^T, with rows [ i^T; j^T; k^T ]:

   rotmat_pure_yaw%r(1:3) = (/   cs,   sn, 0d0 /)
   rotmat_pure_yaw%r(4:6) = (/  -sn,   cs, 0d0 /)
   rotmat_pure_yaw%r(7:9) = (/  0d0,  0d0, 1d0 /)

end function rotmat_pure_yaw

!------------------------------------------------------------------------------------------------------------

function rotmat_pure_pitch(pitch)
!--function: initialize rotation matrix at R = R_y(pitch)
   implicit none
!--result value
   type(t_rotmat)              :: rotmat_pure_pitch
!--subroutine arguments
   real(kind=8),   intent(in)  :: pitch
!--local variables
   real(kind=8)                :: cs, sn

   cs = cos(pitch)
   sn = sin(pitch)

!  note: this code shows the transpose R^T, with rows [ i^T; j^T; k^T ]:

   rotmat_pure_pitch%r(1:3) = (/   cs,  0d0, -sn /)
   rotmat_pure_pitch%r(4:6) = (/  0d0,  1d0, 0d0 /)
   rotmat_pure_pitch%r(7:9) = (/   sn,  0d0,  cs /)

end function rotmat_pure_pitch

!------------------------------------------------------------------------------------------------------------

function rotmat_roll_yaw(roll, yaw)
!--function: compute rotation matrix for roll-and-yaw angles
!            roll is about original x-axis, yaw about modified z-axis
   implicit none
!--result value
   type(t_rotmat)              :: rotmat_roll_yaw
!--subroutine arguments
   real(kind=8),   intent(in)  :: roll, yaw    ! rotation angle [rad]

   ! note: this code shows the transpose R^T, with rows [ i^T; j^T; k^T ]:

   rotmat_roll_yaw%r(1:3) = (/  cos(yaw),  cos(roll) * sin(yaw), sin(roll) * sin(yaw) /)
   rotmat_roll_yaw%r(4:6) = (/ -sin(yaw),  cos(roll) * cos(yaw), sin(roll) * cos(yaw) /)
   rotmat_roll_yaw%r(7:9) = (/      0d0 , -sin(roll)           , cos(roll)            /)

   ! equivalent to: rotmat_roll_yaw = rotmat_pure_roll(roll) * rotmat_pure_yaw(yaw)

end function rotmat_roll_yaw

!------------------------------------------------------------------------------------------------------------

function rotmat_transpose(R)
!--function: compute transpose of rotation matrix R, Rout = R^T
   implicit none
!--result value
   type(t_rotmat)              :: rotmat_transpose
!--subroutine arguments
   class(t_rotmat), intent(in) :: R

!  note: this code shows the transpose Rout^T, with rows [ iout^T; jout^T; kout^T ]:

   rotmat_transpose%r(1:3) = (/ R%r(1),  R%r(4), R%r(7) /)
   rotmat_transpose%r(4:6) = (/ R%r(2),  R%r(5), R%r(8) /)
   rotmat_transpose%r(7:9) = (/ R%r(3),  R%r(6), R%r(9) /)

end function rotmat_transpose

!------------------------------------------------------------------------------------------------------------

function rotmat_mirror_y(R)
!--function: compute mirrored version of rotation matrix R, y' = -y
   implicit none
!--result value
   type(t_rotmat)              :: rotmat_mirror_y
!--subroutine arguments
   class(t_rotmat), intent(in) :: R

!  note: this code shows the transpose Rout^T, with rows [ iout^T; jout^T; kout^T ]:

   rotmat_mirror_y%r(1:3) = (/  R%r(1), -R%r(4),  R%r(7) /)
   rotmat_mirror_y%r(4:6) = (/ -R%r(2),  R%r(5), -R%r(8) /)
   rotmat_mirror_y%r(7:9) = (/  R%r(3), -R%r(6),  R%r(9) /)

end function rotmat_mirror_y

!------------------------------------------------------------------------------------------------------------

function rotmat_vec(R, vin)
!--function: compute matrix vector product vout = R * vin
   implicit none
!--result value
   type(t_vec)                 :: rotmat_vec
!--subroutine arguments
   type(t_rotmat), intent(in)  :: R
   type(t_vec),    intent(in)  :: vin

   rotmat_vec%v(1:3) = R%r(1:3) * vin%v(1) + R%r(4:6) * vin%v(2) + R%r(7:9) * vin%v(3)

end function rotmat_vec

!------------------------------------------------------------------------------------------------------------

function rotmat_transp_vec(R, vin)
!--function: compute matrix vector product vout = R^T * vin
   implicit none
!--result value
   type(t_vec)                 :: rotmat_transp_vec
!--subroutine arguments
   type(t_rotmat), intent(in)  :: R
   type(t_vec),    intent(in)  :: vin

   rotmat_transp_vec%v(1) = R%r(1) * vin%v(1) + R%r(2) * vin%v(2) + R%r(3) * vin%v(3)
   rotmat_transp_vec%v(2) = R%r(4) * vin%v(1) + R%r(5) * vin%v(2) + R%r(6) * vin%v(3)
   rotmat_transp_vec%v(3) = R%r(7) * vin%v(1) + R%r(8) * vin%v(2) + R%r(9) * vin%v(3)

end function rotmat_transp_vec

!------------------------------------------------------------------------------------------------------------

function rotmat_product(R1, R2)
!--function: multiply rotation matrices Rprod = R1 * R2
   implicit none
!--result value
   type(t_rotmat)              :: rotmat_product
!--subroutine arguments
   type(t_rotmat), intent(in)  :: R1, R2
!--local variables
   integer                     :: i, j, k, ij, ik, kj
   associate( Rprod => rotmat_product )

   do j = 1, 3
      do i = 1, 3
         ! r(i,j) = sum_{k=1:3} r1(i,k) * r2(k,j)
         ij = i + 3*(j-1)
         Rprod%r(ij) = 0d0
         do k = 1, 3
            ik = i + 3*(k-1)
            kj = k + 3*(j-1)
            Rprod%r(ij) = Rprod%r(ij) + R1%r(ik) * R2%r(kj)
         enddo
      enddo
   enddo

   end associate
end function rotmat_product

!------------------------------------------------------------------------------------------------------------

function rotmat_transp_product(R1, R2)
!--function: multiply rotation matrices Rprod = R1^T * R2
   implicit none
!--result value
   type(t_rotmat)              :: rotmat_transp_product
!--subroutine arguments
   type(t_rotmat), intent(in)  :: R1, R2
!--local variables
   integer                     :: i, j, k, ij, ki, kj
   associate( Rprod => rotmat_transp_product )

   do j = 1, 3
      do i = 1, 3
         ! r(i,j) = sum_{k=1:3} r1(k,i) * r2(k,j)
         ij = i + 3*(j-1)
         Rprod%r(ij) = 0d0
         do k = 1, 3
            ki = k + 3*(i-1)
            kj = k + 3*(j-1)
            Rprod%r(ij) = Rprod%r(ij) + R1%r(ki) * R2%r(kj)
         enddo
      enddo
   enddo

   end associate
end function rotmat_transp_product

!------------------------------------------------------------------------------------------------------------

function rotmat_get_rollangle(r)
!--function: compute roll of rotation matrix assuming roll-yaw-pitch convention
   implicit none
!--function result:
   real(kind=8) rotmat_get_rollangle
!--subroutine arguments
   class(t_rotmat)     :: r
!--local variables
   real(kind=8)        :: roll, yaw

   yaw  = asin( min(1d0, max(-1d0, -r%el(1,2) )))
   roll = asin( min(1d0, max(-1d0,  r%el(3,2) / cos(yaw) )))

   rotmat_get_rollangle = roll

end function rotmat_get_rollangle

!------------------------------------------------------------------------------------------------------------

function rotmat_get_yawangle(r)
!--function: compute yaw of rotation matrix assuming roll-yaw-pitch convention
   implicit none
!--function result:
   real(kind=8) rotmat_get_yawangle
!--subroutine arguments
   class(t_rotmat)     :: r
!--local variables
   real(kind=8)        :: yaw

   yaw  = asin( min(1d0, max(-1d0, -r%el(1,2) )))

   rotmat_get_yawangle = yaw

end function rotmat_get_yawangle

!------------------------------------------------------------------------------------------------------------

function rotmat_get_pitchangle(r)
!--function: compute pitch of rotation matrix assuming roll-yaw-pitch convention
   implicit none
!--function result:
   real(kind=8) rotmat_get_pitchangle
!--subroutine arguments
   class(t_rotmat)     :: r
!--local variables
   real(kind=8)        :: yaw, pitch

   yaw   = asin( min(1d0, max(-1d0, -r%el(1,2) )))
   pitch = atan2( r%el(1,3) , r%el(1,1) )

   rotmat_get_pitchangle = pitch

end function rotmat_get_pitchangle

!------------------------------------------------------------------------------------------------------------

subroutine rotmat_print(rot, nam, idebug)
!--function: print information on a rotation matrix
   implicit none
!--subroutine arguments
   class(t_rotmat)      :: rot
   character(len=*)     :: nam
   integer              :: idebug
!--local variables
   character(len=80)    :: spaces = ' '

   if (idebug.eq.3) then
      ! RecurDyn format:
      write(bufout,*)              nam,   ' = [', rot%r(1),',', rot%r(4),',', rot%r(7),']'
      call write_log(1, bufout)
      write(bufout,*) spaces(1:len(nam)), '   |', rot%r(2),',', rot%r(5),',', rot%r(8),'|'
      call write_log(1, bufout)
      write(bufout,*) spaces(1:len(nam)), '   [', rot%r(3),',', rot%r(6),',', rot%r(9),'].'
      call write_log(1, bufout)

   elseif (idebug.eq.4) then
      ! RecurDyn format, cntc_x/y/z:
      write(bufout,*) nam,'_i = [', rot%r(1),'], ',nam,'_j=[',rot%r(4),'], ', nam,'_k=[',rot%r(7),']'
      call write_log(1, bufout)
      write(bufout,*) '         [', rot%r(2),'],        [',rot%r(5),'],        [',rot%r(8),']'
      call write_log(1, bufout)
      write(bufout,*) '         [', rot%r(3),'],        [',rot%r(6),'],        [',rot%r(9),']'
      call write_log(1, bufout)

   elseif (idebug.ge.1) then
      write(bufout,*) 'rotmat "',trim(nam),'":'
      call write_log(1, bufout)
      write(bufout,123) '    (', rot%r(1),',', rot%r(4),',', rot%r(7),')'
      call write_log(1, bufout)
      write(bufout,123) '    (', rot%r(2),',', rot%r(5),',', rot%r(8),')'
      call write_log(1, bufout)
      write(bufout,123) '    (', rot%r(3),',', rot%r(6),',', rot%r(9),')'
      call write_log(1, bufout)
 123  format(1x,3(a,f11.6),a)
   endif

   if (idebug.eq.2) then
      write(bufout,123) ' roll=',rot%roll(),', yaw=',rot%yaw(),', pitch=',rot%pitch(),' [rad]'
      call write_log(1, bufout)
   endif

end subroutine rotmat_print

!------------------------------------------------------------------------------------------------------------

function marker111(x, y, z, rot)
!--function: create vector v from position-values and orientation matrix
   implicit none
!--result value
   type(t_marker)                       :: marker111
!--subroutine arguments
   real(kind=8),   intent(in)           :: x, y, z
   type(t_rotmat), intent(in), optional :: rot

   marker111%o = t_vec( (/ x, y, z /) )
   if (present(rot)) then
      marker111%rot = rot
   else
      marker111%rot = rotmat_identity()
   endif
end function marker111

!------------------------------------------------------------------------------------------------------------

function marker3(xyz, rot)
!--function: create marker m from position-array and orientation matrix
   implicit none
!--result value
   type(t_marker)                       :: marker3
!--subroutine arguments
   real(kind=8),   intent(in)           :: xyz(3)
   type(t_rotmat), intent(in), optional :: rot

   marker3%o = t_vec( xyz(1:3) )
   if (present(rot)) then
      marker3%rot = rot
   else
      marker3%rot = rotmat_identity()
   endif
end function marker3

!------------------------------------------------------------------------------------------------------------

function markerv(o, rot)
!--function: create marker m from position-vector and orientation matrix
   implicit none
!--result value
   type(t_marker)                       :: markerv
!--subroutine arguments
   type(t_vec),    intent(in)           :: o
   type(t_rotmat), intent(in), optional :: rot

   markerv%o = o
   if (present(rot)) then
      markerv%rot = rot
   else
      markerv%rot = rotmat_identity()
   endif
end function markerv

!------------------------------------------------------------------------------------------------------------

subroutine marker_init(m)
!--function: initialize a marker: o == 0, R == I
   implicit none
!--subroutine arguments
   type(t_marker)      :: m

   m%o   = t_vec( (/ 0d0, 0d0, 0d0 /) )
   m%rot = rotmat_identity()

end subroutine marker_init

!------------------------------------------------------------------------------------------------------------

subroutine marker_shift(m, dx, dy, dz)
!--function: shift a marker: o := o + [dx; dy; dz]
   implicit none
!--subroutine arguments
   type(t_marker)  :: m
   real(kind=8)    :: dx, dy, dz

   m%o%v(1) = m%o%v(1) + dx
   m%o%v(2) = m%o%v(2) + dy
   m%o%v(3) = m%o%v(3) + dz

end subroutine marker_shift

!------------------------------------------------------------------------------------------------------------

subroutine marker_roll(m, roll, yc_arg, zc_arg)
!--function: rotate a marker by roll angle roll [rad] (about x-axis/point [yc;zc])
   implicit none
!--subroutine arguments
   type(t_marker)                     :: m
   real(kind=8), intent(in)           :: roll           ! rotation angle [rad]
   real(kind=8), intent(in), optional :: yc_arg, zc_arg ! rotation origin, default (0,0)
!--local variables
   real(kind=8)             :: yc, zc, cs, sn, yrel, zrel

   yc = 0d0
   zc = 0d0
   if (present(yc_arg)) yc = yc_arg
   if (present(zc_arg)) zc = zc_arg

   ! apply rotation of origin m%o

   cs = cos(roll)
   sn = sin(roll)

   yrel = m%o%v(2) - yc
   zrel = m%o%v(3) - zc
   m%o%v(2) = yc + cs * yrel - sn * zrel
   m%o%v(3) = zc + sn * yrel + cs * zrel

   ! apply effect on orientation m%rot

   m%rot = m%rot * rotmat_pure_roll(roll)

end subroutine marker_roll

!------------------------------------------------------------------------------------------------------------

subroutine marker_yaw(m, yaw, xc_arg, yc_arg)
!--function: rotate a marker by pure yaw angle yaw [rad] (about z-axis/point [xc;yc])
   implicit none
!--subroutine arguments
   type(t_marker)                     :: m
   real(kind=8), intent(in)           :: yaw            ! rotation angle [rad]
   real(kind=8), intent(in), optional :: xc_arg, yc_arg ! rotation origin, default (0,0)
!--local variables
   real(kind=8)             :: xc, yc, cs, sn, xrel, yrel

   xc = 0d0
   yc = 0d0
   if (present(xc_arg)) xc = xc_arg
   if (present(yc_arg)) yc = yc_arg

!  apply rotation of origin m%o

   cs = cos(yaw)
   sn = sin(yaw)

   xrel = m%o%v(1) - xc
   yrel = m%o%v(2) - yc
   m%o%v(1) = xc + cs * xrel - sn * yrel
   m%o%v(2) = yc + sn * xrel + cs * yrel

!  apply rotation on orientation m%rot

   m%rot = m%rot * rotmat_pure_yaw(yaw)

end subroutine marker_yaw

!------------------------------------------------------------------------------------------------------------

subroutine marker_pitch(m, pitch, xc, zc)
!--function: rotate a marker by pure pitch angle pitch [rad] (about y-axis/point [xc;zc])
   implicit none
!--subroutine arguments
   type(t_marker)           :: m
   real(kind=8), intent(in) :: pitch    ! rotation angle [rad]
   real(kind=8), intent(in) :: xc, zc   ! rotation origin
!--local variables
   real(kind=8)             :: cs, sn, xrel, zrel

!  apply rotation of origin m%o

   cs = cos(pitch)
   sn = sin(pitch)

   xrel = m%o%v(1) - xc
   zrel = m%o%v(3) - zc
   m%o%v(1) = xc + cs * xrel + sn * zrel
   m%o%v(3) = zc - sn * xrel + cs * zrel

!  apply rotation on orientation m%rot

   m%rot = m%rot * rotmat_pure_pitch(pitch)

end subroutine marker_pitch

!------------------------------------------------------------------------------------------------------------

subroutine marker_rotate_rotmat(m, r, xc, yc, zc)
!--function: rotate a marker about rotation center [xc;yc;zc] by multiplying with rotation matrix [rad]
   implicit none
!--subroutine arguments
   type(t_marker)           :: m
   type(t_rotmat)           :: r
   real(kind=8), intent(in) :: xc, yc, zc         ! rotation origin
!--local variables
   type(t_vec)              :: vc, vrel, vrot

!  o_new = vc + R * ( o_old - vc )

   vc    = vec( xc, yc, zc )
   vrel  = m%o - vc
   vrot  = r * vrel
   m%o   = vc + vrot

!  R_new = R * R_old

   m%rot = r * m%rot

end subroutine marker_rotate_rotmat

!------------------------------------------------------------------------------------------------------------

subroutine marker_rotate_eulang(m, roll, yaw, pitch, xc, yc, zc)
!--function: rotate a marker by roll-yaw-pitch euler angles [rad] about rotation center [xc;yc;zc]
   implicit none
!--subroutine arguments
   type(t_marker)           :: m
   real(kind=8), intent(in) :: roll, yaw, pitch   ! rotation angles [rad]
   real(kind=8), intent(in) :: xc, yc, zc         ! rotation origin
!--local variables
   type(t_rotmat)           :: rot

!  compute rotation-matrix for roll-yaw-pitch angles

   rot = rotmat_pure_roll(roll) * rotmat_pure_yaw(yaw) * rotmat_pure_pitch(pitch)

!  delegate actual rotation to marker_rotate_rotmat

   call marker_rotate_rotmat(m, rot, xc, yc, zc)

end subroutine marker_rotate_eulang

!------------------------------------------------------------------------------------------------------------

function marker_transpose(mref_glb)
!--function: compute "transpose of marker", i.e. the conversion from global to local coordinates
!            mref_glb:  xglb = oref + Rref * xref
!            mglb_ref:  xref = oglb + Rglb * xglb  = -Rref^T * oref + Rref^T * xglb
   implicit none
!--function result:
   type(t_marker)      :: marker_transpose
!--subroutine arguments
   type(t_marker)      :: mref_glb
!--local variables
   type(t_vec)         :: oneg
   associate( mglb_ref => marker_transpose )

!  the "inverse" rotation matrix is \bar{R} = R^T

   mglb_ref%rot = rotmat_transpose( mref_glb%rot )

!  the "inverse" origin is \bar{o} = - R^T o

   oneg = dble(-1d0) * mref_glb%o
   mglb_ref%o = rotmat_vec( mglb_ref%rot, oneg )

   end associate
end function marker_transpose

!------------------------------------------------------------------------------------------------------------

function marker_mirror_y(m)
!--function: compute "mirrored version of marker", with y' = -y
   implicit none
!--function result:
   type(t_marker)      :: marker_mirror_y
!--subroutine arguments
   type(t_marker)      :: m

   marker_mirror_y%rot = rotmat_mirror_y( m%rot )
   marker_mirror_y%o   = vec( m%o%v(1), -m%o%v(2), m%o%v(3) )

end function marker_mirror_y

!------------------------------------------------------------------------------------------------------------

function marker_2glob(mloc_ref, mref_glb)
!--function: compute local-to-global conversion for a marker defined with respect to mref
!            mref_glb:  xglb = oref + Rref * xref
!            mloc_ref:                       xref = oloc + Rloc * xloc
!            mloc_glb:  xglb = oref + Rref *      ( oloc + Rloc * xloc )
!            the output corresponds to mloc_glb defined with respect to global system
   implicit none
!--function result:
   type(t_marker)      :: marker_2glob
!--subroutine arguments
   type(t_marker), intent(in)  :: mloc_ref, mref_glb

!  transform origin o = oref + Rref * oloc

   marker_2glob%o = vec_2glob_m(mloc_ref%o, mref_glb)

!  transform rotation matrix R = Rref * Rloc

   marker_2glob%rot = mref_glb%rot * mloc_ref%rot

end function marker_2glob

!------------------------------------------------------------------------------------------------------------

function marker_2loc(mloc_glb, mref_glb)
!--function: compute global-to-local conversion for a marker defined with respect to global system
!            for given markers A = mloc_glb, B = mref_glb defined as local systems with respect to
!            a global reference C, determine the marker A* with respect to B.
!
!         B: mref_glb:                       xglb   =   oref + Rref * xref
!         A: mloc_glb:  oloc + Rloc * xloc = xglb
!        A*: mloc_ref:                xloc = Rloc^T * ( oref + Rref * xref - oloc )
!                                          = Rloc^T * ( oref - oloc ) + Rloc^T * Rref * xref
!            the output mloc_ref corresponds to mloc_glb defined with respect to the reference
   implicit none
!--function result:
   type(t_marker)              :: marker_2loc
!--subroutine arguments
   type(t_marker), intent(in)  :: mloc_glb, mref_glb
!--local variables
   type(t_marker)              :: mglb_ref

   ! compute transpose of B = mref_glb, i.e. the marker for the global system C in terms of the reference B

   mglb_ref = marker_transpose(mref_glb)

   ! transform A = mloc_glb from the new local system C = 'glb' to the new global system B = 'ref'

   marker_2loc = marker_2glob( mloc_glb, mglb_ref )

end function marker_2loc

!------------------------------------------------------------------------------------------------------------

function marker_transp_vec(m, vin)
!--function: rotate a vector by the rotmat of a marker: compute vector product vout = m%R^T * vin
   implicit none
!--result value
   type(t_vec)                 :: marker_transp_vec
!--subroutine arguments
   type(t_marker), intent(in)  :: m
   type(t_vec),    intent(in)  :: vin

   associate( R => m%rot )
   marker_transp_vec%v(1) = R%r(1) * vin%v(1) + R%r(2) * vin%v(2) + R%r(3) * vin%v(3)
   marker_transp_vec%v(2) = R%r(4) * vin%v(1) + R%r(5) * vin%v(2) + R%r(6) * vin%v(3)
   marker_transp_vec%v(3) = R%r(7) * vin%v(1) + R%r(8) * vin%v(2) + R%r(9) * vin%v(3)
   end associate

end function marker_transp_vec

!------------------------------------------------------------------------------------------------------------

function marker_get_x(this)
!--function: get marker x-position o_x
   implicit none
!--result value
   real(kind=8)                :: marker_get_x
!--subroutine arguments
   class(t_marker), intent(in) :: this

   marker_get_x = this%o%v(1)
end function marker_get_x

!------------------------------------------------------------------------------------------------------------

function marker_get_y(this)
!--function: get marker y-position o_y
   implicit none
!--result value
   real(kind=8)                :: marker_get_y
!--subroutine arguments
   class(t_marker), intent(in) :: this

   marker_get_y = this%o%v(2)
end function marker_get_y

!------------------------------------------------------------------------------------------------------------

function marker_get_z(this)
!--function: get marker z-position o_z
   implicit none
!--result value
   real(kind=8)                :: marker_get_z
!--subroutine arguments
   class(t_marker), intent(in) :: this

   marker_get_z = this%o%v(3)
end function marker_get_z

!------------------------------------------------------------------------------------------------------------

function marker_get_ivec(this)
!--function: get first column i of rotation matrix
   implicit none
!--result value
   type(t_vec)                 :: marker_get_ivec
!--subroutine arguments
   class(t_marker), intent(in) :: this

   marker_get_ivec = t_vec( this%rot%r(1:3) )
end function marker_get_ivec

!------------------------------------------------------------------------------------------------------------

function marker_get_jvec(this)
!--function: get first column i of rotation matrix
   implicit none
!--result value
   type(t_vec)                 :: marker_get_jvec
!--subroutine arguments
   class(t_marker), intent(in) :: this

   marker_get_jvec = t_vec( this%rot%r(4:6) )
end function marker_get_jvec

!------------------------------------------------------------------------------------------------------------

function marker_get_kvec(this)
!--function: get first column i of rotation matrix
   implicit none
!--result value
   type(t_vec)                 :: marker_get_kvec
!--subroutine arguments
   class(t_marker), intent(in) :: this

   marker_get_kvec = t_vec( this%rot%r(7:9) )
end function marker_get_kvec

!------------------------------------------------------------------------------------------------------------

function marker_get_rollangle(m)
!--function: compute roll of rotation matrix assuming roll-yaw-pitch convention
   implicit none
!--function result:
   real(kind=8) marker_get_rollangle
!--subroutine arguments
   class(t_marker)     :: m

   marker_get_rollangle = m%rot%roll()

end function marker_get_rollangle

!------------------------------------------------------------------------------------------------------------

function marker_get_yawangle(m)
!--function: compute yaw of rotation matrix assuming roll-yaw-pitch convention
   implicit none
!--function result:
   real(kind=8) marker_get_yawangle
!--subroutine arguments
   class(t_marker)     :: m

   marker_get_yawangle = m%rot%yaw()

end function marker_get_yawangle

!------------------------------------------------------------------------------------------------------------

function marker_get_pitchangle(m)
!--function: compute pitch of rotation matrix assuming roll-yaw-pitch convention
   implicit none
!--function result:
   real(kind=8) marker_get_pitchangle
!--subroutine arguments
   class(t_marker)     :: m

   marker_get_pitchangle = m%rot%pitch()

end function marker_get_pitchangle

!------------------------------------------------------------------------------------------------------------

subroutine marker_print(m, nam, idebug)
!--function: print information on a marker
   implicit none
!--subroutine arguments
   type(t_marker)       :: m
   character(len=*)     :: nam
   integer              :: idebug
!--local variables
   real(kind=8)         :: roll, yaw, pitch

   if (idebug.ge.1) then
      write(bufout,*) 'marker "',trim(nam),'":'
      call write_log(1, bufout)
      write(bufout,123) '    o = (', m%o%v(1),',', m%o%v(2),',', m%o%v(3),')^T'
      call write_log(1, bufout)
   endif

   if (idebug.ge.3) then
      write(bufout,123) '    R = (', m%rot%el(1,1),',', m%rot%el(1,2),',', m%rot%el(1,3),')'
      call write_log(1, bufout)
      write(bufout,123) '        (', m%rot%el(2,1),',', m%rot%el(2,2),',', m%rot%el(2,3),')'
      call write_log(1, bufout)
      write(bufout,123) '        (', m%rot%el(3,1),',', m%rot%el(3,2),',', m%rot%el(3,3),')'
      call write_log(1, bufout)
 123  format(1x,3(a,f11.6),a)
   endif

   if (idebug.ge.2) then
      yaw   = asin( min(1d0, max(-1d0, -m%rot%el(1,2) )))
      roll  = asin( min(1d0, max(-1d0,  m%rot%el(3,2) / cos(yaw) )))
      pitch = atan2( m%rot%el(1,3) , m%rot%el(1,1) )
      write(bufout,123) '    roll=',roll,', yaw=',yaw,', pitch=',pitch,' [rad]'
      call write_log(1, bufout)
   endif

end subroutine marker_print

!------------------------------------------------------------------------------------------------------------

subroutine markers_test
!--function: test global-to-local and local-to-global conversions
   implicit none
!--local variables:
   real(kind=8)   :: pi = 4d0*atan(1d0)
   type(t_marker) :: m_fixed, m_cntc, mf_cntc, m_glob, mf_glob, m_rdyn
   type(t_vec)    :: o_fixed, o_cntc, oc_fixed, o_glob, ivec, jvec, kvec, p_fixed, p_rdyn, p_cntc, p_glob
   type(t_rotmat) :: r_fixed, r_cntc, r_glob

   call write_log(' ----------------------------- Markers_test: ------------------------------')
   !  - construct marker m_fixed in recurdyn coordinates: of_r = [ 20, -15, 0 ], 
   !                                                      Rf_r = [0.87,-0.5,0; 0.5,0.87,0; 0,0,1]

   o_fixed = vec( 20d0, -15d0, 0d0 )
   r_fixed = rotmat_pure_yaw( pi/6 )
   m_fixed = t_marker(o_fixed, r_fixed)
   call marker_print(m_fixed, 'm_fixed', 2)

   ! - create point P in fixed-roller coordinates: p_f = [ 10, 5, 0 ]

   p_fixed = vec( 10d0, 5d0, 0d0 )
   call vec_print(p_fixed, 'p_fixed', 2)

   ! - convert P to recurdyn coordinates: p_r = of_r + Rf_r * p_f = [ 26.2, -5.7, 0 ]

   p_rdyn = vec_2glob(p_fixed, m_fixed)
   call vec_print(p_rdyn, 'p_rdyn', 2)

   !  - construct marker m_cntc in recurdyn coordinates: oc_r = [ 20, 15, 0 ], 
   !                                                     Rc_r = [ -1,0,0; 0,0,1; 0,1,0 ]

   o_cntc  = vec( 20d0,  15d0, 0d0 )
   ivec    = vec( -1d0, 0d0, 0d0 )
   jvec    = vec(  0d0, 0d0, 1d0 )
   kvec    = vec(  0d0, 1d0, 0d0 )
   r_cntc  = rotmat( ivec, jvec, kvec )
   m_cntc  = t_marker(o_cntc, r_cntc)
   call marker_print(m_cntc, 'm_cntc', 2)

   ! - convert O_c to fixed-roller coordinates: oc_f = Rf_r^T * ( oc_r - of_r ) = [ 15.0, 26.0, 0 ]

   oc_fixed = vec_2loc(o_cntc, m_fixed)
   call vec_print(oc_fixed, 'o_c wrt m_fixed', 2)

   ! - convert marker m_fixed to contact coordinates: of_c = Rc_r^T * (of_r - oc_r) = [ 0, 0, -30 ]
   !                                      Rf_c = Rc_r^T * Rf_r = [ -0.87,0,0.5; 0.5,0,0.87; 0,1,0 ]

   mf_cntc = marker_2loc(m_fixed, m_cntc)
   call marker_print(mf_cntc, 'mf_cntc', 2)

   ! - convert point P to contact coordinates:  p_c = of_c + Rf_c * p_f = [ -6.2, 0, -20.7 ]
   !                                            p_c = Rc_r^T * ( p_r - oc_r ) = [ -6.2,  0, -20.7 ]

   p_cntc = vec_2glob(p_fixed, mf_cntc)
   call vec_print(p_cntc, 'p_cntc1, via mf_cntc', 2)

   p_cntc = vec_2loc(p_rdyn, m_cntc)
   call vec_print(p_cntc, 'p_cntc2, via m_rdyn', 2)

   ! - construct marker m_glob at og_r = [ -15, 0, 0 ] at angle -15deg (-pi/12 rad) w.r.t. m_rdyn
   !                                            Rg_r = [ 0.97,0.26,0; -0.26,0.97,0; 0,0,1 ]

   o_glob = vec( -15d0, 0d0, 0d0 )
   r_glob = rotmat_pure_yaw( -pi/12 )
   m_glob = t_marker(o_glob, r_glob)
   call marker_print(m_glob, 'm_glob', 2)

   ! - compute transpose: position of m_rdyn within global coordinates: 
   !                                            or_g = -Rg_r^T * og_r = [ 14.5, 3.9, 0 ]
   !                                            Rr_g =  Rg_r^T = [ 0.97,-0.26,0; 0.26,0.97,0; 0,0,1 ]

   m_rdyn = marker_transpose(m_glob)
   call marker_print(m_rdyn, 'm_rdyn', 2)

   ! - convert marker m_fixed to global coordinates: of_g = or_g + Rr_g * of_r = [ 37.7, -5.4, 0 ]
   !                                    Rf_g = Rr_g * Rf_r = [ 0.71,-0.71,0; 0.71,0.71,0; 0,0,1 ]

   mf_glob = marker_2glob(m_fixed, m_rdyn)
   call marker_print(mf_glob, 'mf_glob', 2)

   ! - convert point P to global coordinates: p_g = of_g + Rf_g * p_f = [ 41.2, 5.2, 0 ]
   !                                          p_g = Rg_r^T * ( p_r - og_r ) = [ 41.2,  5.2, 0 ]

   p_glob = vec_2glob(p_fixed, mf_glob)
   call vec_print(p_glob, 'p_glob1, via mf_glob', 2)

   p_glob = vec_2loc(p_rdyn, m_glob)
   call vec_print(p_glob, 'p_glob2, via m_glob', 2)
   call write_log(' ---------------------------- end Markers_test -----------------------------')

end subroutine markers_test

!------------------------------------------------------------------------------------------------------------

end module m_markers
