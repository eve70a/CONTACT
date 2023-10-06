!------------------------------------------------------------------------------------------------------------
! m_grids - data-structures for grids - uniform / non-uniform
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_grids
   use m_globals
   use m_markers
   use m_ptrarray
   use m_interp_1d
   use m_spline_def
   use m_spline_make
   use m_bspline_get
   implicit none
   private

   ! Data type for different kinds of grids:

   public t_grid
   public p_grid

   ! Functions defined on grids:

   public  grid_nullify
   public  grid_set_dimens
   public  grid_get_ix4ii
   public  grid_get_iy4ii
   public  grid_get_dimens
   public  grid_get_xrange
   public  grid_get_yrange
   public  grid_get_zrange
   public  grid_get_boundbox
   public  grid_check_nan

   public  grid_create_uniform
   public  grid_create_reg_cyl
   public  grid_create_cylindr
   public  grid_create_curvil

   public  grid_has_spline
   public  grid_alloc_splines
   public  grid_xspline_init

   public  grid_make_arclength
   public  grid_make_ppspline
   public  grid_make_ppspline_kink
   public  grid_make_ppspline_nokink
   public  grid_make_bspline

   interface grid_make_ppspline
      module procedure grid_make_ppspline_kink
      module procedure grid_make_ppspline_nokink
   end interface grid_make_ppspline

   public  grid_extrude_profile
   public  grid_revolve_profile
   public  grid_spline_copy
   public  grid_copy
   public  grid_trim
   public  grid_print
   public  grid_destroy
   public  grid_destroy_spline

   interface grid_extrude_profile
      module procedure grid_extrude_profile_xarr
      module procedure grid_extrude_profile_xunif
   end interface grid_extrude_profile

   interface grid_revolve_profile
      module procedure grid_revolve_profile_xarr
      module procedure grid_revolve_profile_xarr_zarr
      module procedure grid_revolve_profile_xunif
   end interface grid_revolve_profile

   public  unifgrid_mirror_y

   public  cartgrid_shift
   public  cartgrid_mirror_y
   public  cartgrid_rotate
   public  cartgrid_roll
   public  cartgrid_yaw
   public  cartgrid_roll_yaw
   public  cartgrid_2glob
   public  cartgrid_2loc

   interface cartgrid_2glob
      module procedure cartgrid_2glob_m
      module procedure cartgrid_2glob_or
   end interface cartgrid_2glob

   interface cartgrid_2loc
      module procedure cartgrid_2loc_m
   end interface cartgrid_2loc

   public convert_cyl2curv

!------------------------------------------------------------------------------------------------------------
!  data with respect to the discretisation grid:

   type :: t_grid
      logical      :: is_uniform, is_curvilin, is_cylindr, lies_in_oyz
      integer      :: nx, ny, ntot, ivdir
      real(kind=8) :: dx, dy, dxdy
      real(kind=8), dimension(:,:), pointer  :: coor  => NULL() ! ntot x 3
      real(kind=8), dimension(:),   pointer  :: x     => NULL()
      real(kind=8), dimension(:),   pointer  :: y     => NULL()
      real(kind=8), dimension(:),   pointer  :: z     => NULL()
      real(kind=8), dimension(:),   pointer  :: r     => NULL()
      real(kind=8), dimension(:),   pointer  :: th    => NULL()
      real(kind=8), dimension(:),   pointer  :: v     => NULL()

      real(kind=8), dimension(:),   pointer  :: s_prf => NULL()

      type(t_spline)                         :: spl

   contains
      procedure :: ix   => grid_get_ix4ii
      procedure :: iy   => grid_get_iy4ii

      ! is_uniform     indicating original uniform CONTACT grid, with fixed dx, dy
      ! is_curvilin    indicating curvilinear surface grid, with [x(i,j), y(i,j), z(i,j)]
      !                e.g. rail surface
      ! is_cylindr     indicating cylindrical product grid [v(nv)] x [th(nth)] -> [r(nv,nth)]
      !                stored as [r,theta,v], [r,v,theta] or [v,r,theta] depending on ivdir
      ! lies_in_oyz    indicating that x_ij == 0, esp. for wheel/rail profiles
      ! ivdir          coordinate axis x=1, y=2 or z=3 to which the cylinder axis is aligned
      !
      ! nx             number of points in grid in x-direction; cylindrical: v running fastest, nx == nv
      ! ny             number of points in grid in y-direction; cylindrical: ny == nth
      ! ntot           total number of points nx * ny
      ! coor   [mm]    coordinate array (1:ntot,3): for each grid point a 3-vector with coordinates
      !                uniform & curvilinear: [x,y,z]-coordinates,
      !                cylindrical: [v,r,th], [r,v,th] or [r,th,v], depending on ivdir
      !
      ! - for uniform grids:
      !   dx, dy [mm]  constant grid step-sizes
      !   dxdy         area of each grid cell for uniform grids
      !
      ! - for uniform & curvilinear grids:
      !   x    [mm]    pointer to the x-coordinates, that is, the first column of the coor array
      !   y    [mm]    pointer to the y-coordinates, that is, the second column of the coor array
      !   z    [mm]    pointer to the z-coordinates, that is, the third column of the coor array
      ! - for cylindrical coordinate grids:
      !   r    [mm]    pointer to the r-coordinates, stored in second (ivdir=1) or first column (2,3) of coor
      !   th   [deg]   pointer to theta-coordinates, stored in second (ivdir=3) or third column (1,2) of coor
      !   v    [mm]    pointer to the v-coordinates, stored in column 'ivdir' of the coor array
      !
      ! - for 1d grids (curves):
      !   s_prf [mm]    arc length parameter in profile points
      !
      ! - data for splines:
      !   spl           structure with spline data

   end type t_grid

   type :: p_grid
      type(t_grid), pointer :: g => NULL()      ! pointer to a grid data-structure
   end type p_grid

contains

!------------------------------------------------------------------------------------------------------------

subroutine grid_nullify(g)
!--purpose: initialize grid structure, nullify pointers
   implicit none
!--subroutine parameters:
   type(t_grid)              :: g

   g%coor   => NULL()
   g%x      => NULL()
   g%y      => NULL()
   g%z      => NULL()
   g%r      => NULL()
   g%th     => NULL()
   g%v      => NULL()

   g%s_prf  => NULL()

   call spline_nullify(g%spl)

end subroutine grid_nullify

!------------------------------------------------------------------------------------------------------------

subroutine grid_set_dimens(g, nx, ny)
!--purpose: Set/update number of points in grid
   implicit none
!--subroutine parameters:
   type(t_grid)              :: g
   integer,      intent(in)  :: nx, ny

   g%nx = nx
   g%ny = ny
   g%ntot = nx * ny
end subroutine grid_set_dimens

!------------------------------------------------------------------------------------------------------------

function grid_get_ix4ii(this, ii)
!--function: compute first 2-d index ix for given 1-d point number ii
   implicit none
!--result value
   integer                   :: grid_get_ix4ii
!--subroutine arguments
   class(t_grid), intent(in) :: this
   integer,       intent(in) :: ii

   grid_get_ix4ii = mod(ii-1, this%nx) + 1
end function grid_get_ix4ii

!------------------------------------------------------------------------------------------------------------

function grid_get_iy4ii(this, ii)
!--function: compute second 2-d index iy for given 1-d point number ii
   implicit none
!--result value
   integer                   :: grid_get_iy4ii
!--subroutine arguments
   class(t_grid), intent(in) :: this
   integer,       intent(in) :: ii

   grid_get_iy4ii = (ii-1) / this%nx + 1
end function grid_get_iy4ii

!------------------------------------------------------------------------------------------------------------

subroutine grid_get_dimens(g, nx, ny)
!--purpose: Get the number of points in grid
   implicit none
!--subroutine parameters:
   type(t_grid)              :: g
   integer,      intent(out) :: nx, ny

   nx = g%nx
   ny = g%ny
end subroutine grid_get_dimens

!------------------------------------------------------------------------------------------------------------

subroutine grid_get_xrange(g, xmin, xmax)
!--purpose: Get the extent of the grid in x-direction
   implicit none
!--subroutine parameters:
   type(t_grid)              :: g
   real(kind=8), intent(out) :: xmin, xmax
!--local variables:
   integer    :: ii

   if (g%lies_in_oyz) then
      xmin =  0d0
      xmax =  0d0
   elseif (g%is_cylindr) then
      xmin =  1d10
      xmax = -1d10
   elseif (g%is_uniform) then
      xmin = g%x( 1 )
      xmax = g%x( g%ntot )
   else
      xmin =  1d10
      xmax = -1d10
      do ii = 1, g%ntot
         if (g%x(ii).lt.xmin) xmin = g%x(ii)
         if (g%x(ii).gt.xmax) xmax = g%x(ii)
      enddo
   endif
end subroutine grid_get_xrange

!------------------------------------------------------------------------------------------------------------

subroutine grid_get_yrange(g, ymin, ymax)
!--purpose: Get the extent of the grid in y-direction
   implicit none
!--subroutine parameters:
   type(t_grid)              :: g
   real(kind=8), intent(out) :: ymin, ymax
!--local variables:
   integer    :: ii

   if (g%is_cylindr) then
      ymin =  1d10
      ymax = -1d10
   elseif (g%is_uniform) then
      ymin = g%y( 1 )
      ymax = g%y( g%ntot )
   else
      ymin =  1d10
      ymax = -1d10
      do ii = 1, g%ntot
         if (g%y(ii).lt.ymin) ymin = g%y(ii)
         if (g%y(ii).gt.ymax) ymax = g%y(ii)
      enddo
   endif
end subroutine grid_get_yrange

!------------------------------------------------------------------------------------------------------------

subroutine grid_get_zrange(g, zmin, zmax)
!--purpose: Get the extent of the grid in z-direction
   implicit none
!--subroutine parameters:
   type(t_grid)              :: g
   real(kind=8), intent(out) :: zmin, zmax
!--local variables:
   integer    :: ii

   if (g%is_cylindr) then
      zmin =  1d10
      zmax = -1d10
   else
      zmin =  1d10
      zmax = -1d10
      do ii = 1, g%ntot
         if (g%z(ii).lt.zmin) zmin = g%z(ii)
         if (g%z(ii).gt.zmax) zmax = g%z(ii)
      enddo
   endif
end subroutine grid_get_zrange

!------------------------------------------------------------------------------------------------------------

subroutine grid_get_boundbox(g, bb)
!--purpose: Get a grid that describes the bounding box of the input grid
!           curvilinear: 8x1 grid with (1:nx,1:ny,1:nz), i.e. x running fastest
!           cylindrical: 4x2 grid with ((1:nr,1:nv),1:nth), i.e. r running fastest
   implicit none
!--subroutine parameters:
   type(t_grid)              :: g, bb
!--local variables:
   integer      :: ii, ik, iv, ir, ith, nx, ny, nv, nth
   real(kind=8) :: vmin(3), vmax(3)
   real(kind=8) :: bbx(8), bby(8), bbz(8), bbv(4), bbth(2), bbr(4,2)

   ! get min/max values for all three coordinate directions

   vmin = (/  1d10,  1d10,  1d10 /)
   vmax = (/ -1d10, -1d10, -1d10 /)
   do ik = 1, 3
      do ii = 1, g%ntot
         if (g%coor(ii,ik).lt.vmin(ik)) vmin(ik) = g%coor(ii,ik)
         if (g%coor(ii,ik).gt.vmax(ik)) vmax(ik) = g%coor(ii,ik)
      enddo
   enddo

   if (.not.g%is_cylindr) then

      ! create curvilinear grid with 8 x 1 points, 1D version of 3D grid
      !    bb%x(nx,ny,nz) -- ix running fastest, iz slowest

      nx  = 8
      ny  = 1

      !        (1,1,1)  (2,1,1)  (1,2,1)  (2,2,1)  (1,1,2)  (2,1,2)  (1,2,2)  (2,2,2)
      bbx = (/ vmin(1), vmax(1), vmin(1), vmax(1), vmin(1), vmax(1), vmin(1), vmax(1) /)
      bby = (/ vmin(2), vmin(2), vmax(2), vmax(2), vmin(2), vmin(2), vmax(2), vmax(2) /)
      bbz = (/ vmin(3), vmin(3), vmin(3), vmin(3), vmax(3), vmax(3), vmax(3), vmax(3) /)
      call grid_create_curvil(bb, nx, ny, bbx, bby, bbz)

   else

      ! create cylindrical grid with 4 points in v-direction, 2 in theta
      !              th0    th1
      !            ---------------
      !         v0 |  r0  |  r0  |
      !         v0 |  r1  |  r1  |
      !         v1 |  r0  |  r0  |
      !         v1 |  r1  |  r1  |
      !            ---------------

      nv  = 4
      nth = 2

      ! r-coordinates are stored in second column if (ivdir=1) or first column if (ivdir=2,3)
      ! th-coordinates are stored in second column if (ivdir=3) or third column if (ivdir=1,2)

      iv  = g%ivdir
      if (iv.eq.1) then
         ir  = 2
         ith = 3
      elseif (iv.eq.2) then
         ir  = 1
         ith = 3
      else
         ir  = 1
         ith = 2
      endif

      bbv  = (/ vmin(iv),  vmin(iv),  vmax(iv),  vmax(iv)  /)
      bbth = (/ vmin(ith), vmax(ith) /) 
      bbr  = reshape( (/ vmin(ir),  vmin(ir),  vmax(ir),  vmax(ir),                                      &
                         vmin(ir),  vmin(ir),  vmax(ir),  vmax(ir)  /),  (/ nv, nth /) )

      call grid_create_cylindr(bb, iv, nv, nth, bbv, bbth, bbr)

   endif

end subroutine grid_get_boundbox

!------------------------------------------------------------------------------------------------------------

subroutine grid_check_nan(g, nam, idebug, nnan_out)
!--purpose: Count number of NaN-values in grid g
   implicit none
!--subroutine parameters:
   type(t_grid), intent(in)            :: g
   character(len=*)                    :: nam
   integer,      intent(in)            :: idebug
   integer,      intent(out), optional :: nnan_out(3)
!--local variables:
   integer      :: ii, ii0(3), ii1(3), ik, nnan(3), nnan_tot

   ! count number of nan values for all three coordinate directions, determine first/last positions ii

   nnan = (/  0,  0,  0 /)
   ii0  = (/ -1, -1, -1 /)
   ii1  = (/ -1, -1, -1 /)
   do ik = 1, 3
      do ii = 1, g%ntot
         if (isnan(g%coor(ii,ik))) then
            nnan(ik) = nnan(ik) + 1
            if (ii0(ik).lt.0) ii0(ik) = ii
            ii1(ik) = ii
         endif
      enddo
   enddo
   nnan_tot = nnan(1) + nnan(2) + nnan(3)

   ! print information as requested by idebug

   if (nnan_tot.le.0 .and. idebug.ge.3) then
      write(bufout,'(3a)') ' Grid ',trim(nam),' has no NaN-values'
      call write_log(1, bufout)
   endif

   if (nnan_tot.ge.1 .and. idebug.ge.1) then
      write(bufout,'(3a,i7,3(a,i7),a)') ' Grid ',trim(nam),' has',nnan_tot,' NaN-values (',     &
                nnan(1),'/',nnan(2),'/',nnan(3),')'
      call write_log(1, bufout)
   endif

   if (idebug.ge.2) then
      do ik = 1, 3
         if (nnan(ik).gt.0) then
            write(bufout,'(a,i1,2(a,i7))') ' Coord.dir.',ik,': first NaN at ii=',ii0(ik),', last at',ii1(ik)
            call write_log(1,bufout)
         endif
      enddo
   endif

   ! copy nnan to nnan_out if necessary

   if (present(nnan_out)) nnan_out(1:3) = nnan(1:3)

end subroutine grid_check_nan

!------------------------------------------------------------------------------------------------------------

subroutine grid_create_uniform(g, nxarg, x0arg, dxarg, x1arg, nyarg, y0arg, dyarg, y1arg, zarg)
!--purpose: Create a regular grid with cartesian coordinates from the specified parameters
   implicit none
!--subroutine parameters:
   type(t_grid)                        :: g
   integer,      intent(in), optional  :: nxarg, nyarg
   real(kind=8), intent(in), optional  :: x0arg, dxarg, x1arg, y0arg, dyarg, y1arg, zarg
!--local variables:
   integer      :: nx, ny, ii, ix, iy
   real(kind=8) :: x0, dx, x1, y0, dy, y1, yy, zz

   ! get input-data for x-direction. Three out of four must be provided: nx, x0, dx, x1

   if (present(nxarg) .and. present(x0arg) .and. present(dxarg)) then
      nx = nxarg
      x0 = x0arg
      dx = dxarg
      x1 = x0 + (nx-1) * dx
   elseif (present(nxarg) .and. present(x0arg) .and. present(x1arg)) then
      nx = nxarg
      x0 = x0arg
      x1 = x1arg
      if (nx.le.1) then
         dx = 1d0
         x1 = x0
      else
         dx = (x1 - x0) / real(nx-1)
      endif
   elseif (present(nxarg) .and. present(dxarg) .and. present(x1arg)) then
      nx = nxarg
      dx = dxarg
      x1 = x1arg
      x0 = x1 - (nx-1) * dx
   elseif (present(x0arg) .and. present(dxarg) .and. present(x1arg)) then
      x0 = x0arg
      dx = dxarg
      x1 = x1arg
      ! Note: using int, the actual x1 will be <= the input value
      nx = int( (x1-x0) / dx + 1d-6 ) + 1
      x1 = x0 + (nx-1) * dx
   else
      call write_log('grid_create_uniform: insufficient parameters for x-direction')
      return
   endif

   ! get input-data for y-direction. Three out of four must be provided: ny, y0, dy, y1

   if (present(nyarg) .and. present(y0arg) .and. present(dyarg)) then
      ny = nyarg
      y0 = y0arg
      dy = dyarg
      y1 = y0 + (ny-1) * dy
   elseif (present(nyarg) .and. present(y0arg) .and. present(y1arg)) then
      ny = nyarg
      y0 = y0arg
      y1 = y1arg
      if (ny.le.1) then
         dy = 1d0
         y1 = y0
      else
         dy = (y1 - y0) / real(ny-1)
      endif
   elseif (present(nyarg) .and. present(dyarg) .and. present(y1arg)) then
      ny = nyarg
      dy = dyarg
      y1 = y1arg
      y0 = y1 - (ny-1) * dy
   elseif (present(y0arg) .and. present(dyarg) .and. present(y1arg)) then
      y0 = y0arg
      dy = dyarg
      y1 = y1arg
      ! Note: using int, the actual y1 will be <= the input value
      ny = int( (y1-y0) / dy + 1d-6 ) + 1
      y1 = y0 + (ny-1) * dy
   else
      call write_log('grid_create_uniform: insufficient parameters for y-direction')
      return
   endif

   ! get input-data for z-direction / set default

   if (present(zarg)) then
      zz = zarg
   else
      zz = 0d0
   endif

   ! store the type of grid

   g%is_uniform  = .true.
   g%is_curvilin = .false.
   g%is_cylindr  = .false.
   g%lies_in_oyz = .false.
   g%ivdir       = 0

   ! store the number of points in the grid

   g%nx   = nx
   g%ny   = ny
   g%ntot = nx * ny

   ! store the grid sizes

   g%dx   = dx
   g%dy   = dy
   g%dxdy = dx * dy

   ! re-allocate coordinate array at the appropriate size

   call reallocate_arr(g%coor, g%ntot, 3)

   ! set pointers to x,y,z-coordinates

   g%x  => g%coor(:,1)
   g%y  => g%coor(:,2)
   g%z  => g%coor(:,3)
   g%r  => NULL()
   g%th => NULL()
   g%v  => NULL()

   ! compute the (x,y,z)-coordinates of all points of the grid

   do iy = 1, ny
      yy = y0 + (iy-1) * dy
      do ix = 1, nx
         ii = ix + (iy-1) * nx
         g%x(ii) = x0 + (ix-1) * dx
         g%y(ii) = yy
         g%z(ii) = zz
      enddo
   enddo

   ! clear previous arc-length s and spline data

   if (associated(g%s_prf)) deallocate(g%s_prf)
   g%s_prf  => NULL()

   call grid_destroy_spline(g)

end subroutine grid_create_uniform

!------------------------------------------------------------------------------------------------------------

subroutine grid_create_reg_cyl(g, ivdir, nvarg, v0arg, dvarg, v1arg, ntharg, th0arg, dtharg, th1arg, rarg)
!--purpose: Create a regular grid with cylindrical coordinates from the specified parameters
   implicit none
!--subroutine parameters:
   type(t_grid)                        :: g
   integer,      intent(in)            :: ivdir
   integer,      intent(in), optional  :: nvarg, ntharg
   real(kind=8), intent(in), optional  :: v0arg, dvarg, v1arg, th0arg, dtharg, th1arg, rarg
!--local variables:
   integer      :: nv, nth, ii, ith, iv
   real(kind=8) :: v0, dv, v1, th0, dth, th1, r, th

   ! check v-direction

   if (ivdir.ne.ikXDIR .and. ivdir.ne.ikYDIR .and. ivdir.ne.ikZDIR) then
      write(bufout,*) 'Internal error: create_reg_cyl: ivdir=',ivdir,' must be X,Y or Z'
      call write_log(1, bufout)
      call abort_run()
   endif

   ! get input-data for v-direction. Three out of four must be provided: nv, v0, dv, v1

   if (present(nvarg) .and. present(v0arg) .and. present(dvarg)) then
      nv = nvarg
      v0 = v0arg
      dv = dvarg
      v1 = v0 + (nv-1) * dv
   elseif (present(nvarg) .and. present(v0arg) .and. present(v1arg)) then
      nv = nvarg
      v0 = v0arg
      v1 = v1arg
      if (nv.le.1) then
         dv = 1d0
         v1 = v0
      else
         dv = (v1 - v0) / real(nv-1)
      endif
   elseif (present(nvarg) .and. present(dvarg) .and. present(v1arg)) then
      nv = nvarg
      dv = dvarg
      v1 = v1arg
      v0 = v1 - (nv-1) * dv
   elseif (present(v0arg) .and. present(dvarg) .and. present(v1arg)) then
      v0 = v0arg
      dv = dvarg
      v1 = v1arg
      nv = nint( (v1-v0) / dv ) + 1
      v1 = v0 + (nv-1) * dv
   else
      call write_log('cylgrid_create_reg_cyl: insufficient parameters for v-direction')
      return
   endif

   ! get input-data for th-direction. Three out of four must be provided: nth, th0, dth, th1

   if (present(ntharg) .and. present(th0arg) .and. present(dtharg)) then
      nth = ntharg
      th0 = th0arg
      dth = dtharg
      th1 = th0 + (nth-1) * dth
   elseif (present(ntharg) .and. present(th0arg) .and. present(th1arg)) then
      nth = ntharg
      th0 = th0arg
      th1 = th1arg
      if (nth.le.1) then
         dth = 1d0
         th1 = th0
      else
         dth = (th1 - th0) / real(nth-1)
      endif
   elseif (present(ntharg) .and. present(dtharg) .and. present(th1arg)) then
      nth = ntharg
      dth = dtharg
      th1 = th1arg
      th0 = th1 - (nth-1) * dth
   elseif (present(th0arg) .and. present(dtharg) .and. present(th1arg)) then
      th0 = th0arg
      dth = dtharg
      th1 = th1arg
      nth = nint( (th1-th0) / dth ) + 1
      th1 = th0 + (nth-1) * dth
   else
      call write_log('cylgrid_create_reg_cyl: insufficient parameters for theta-direction')
      return
   endif

   ! get input-data for r-direction / set default

   if (present(rarg)) then
      r = rarg
   else
      r = 1d0
   endif

   ! store the type of grid

   g%is_uniform  = .false.
   g%is_curvilin = .false.
   g%is_cylindr  = .true.
   g%lies_in_oyz = .false.
   g%ivdir       = ivdir

   ! store the number of points in the grid

   g%nx   = nv
   g%ny   = nth
   g%ntot = nv * nth

   ! store the grid sizes

   g%dx   = dv
   g%dy   = dth
   g%dxdy = dv * dth

   ! re-allocate coordinate array at the appropriate size

   call reallocate_arr(g%coor, g%ntot, 3)

   ! set pointers to r,th,v-coordinates

   if (ivdir.eq.ikXDIR) then
      g%v  => g%coor(:,1)
      g%r  => g%coor(:,2)
      g%th => g%coor(:,3)
   elseif (ivdir.eq.ikYDIR) then
      g%r  => g%coor(:,1)
      g%v  => g%coor(:,2)
      g%th => g%coor(:,3)
   elseif (ivdir.eq.ikZDIR) then
      g%r  => g%coor(:,1)
      g%th => g%coor(:,2)
      g%v  => g%coor(:,3)
   else
      call write_log(' Internal error: create_reg_cyl: incorrect ivdir')
   endif
   g%x  => NULL()
   g%y  => NULL()
   g%z  => NULL()

   ! compute the (r,th,v)-coordinates of all points of the grid

   do ith = 1, nth
      th = th0 + (ith-1) * dth
      do iv = 1, nv
         ii = iv + (ith-1) * nv
         g%r(ii)  = r
         g%th(ii) = th
         g%v(ii)  = v0 + (iv-1) * dv
      enddo
   enddo

   ! clear previous arc-length s and spline data

   if (associated(g%s_prf)) deallocate(g%s_prf)
   g%s_prf  => NULL()

   call grid_destroy_spline(g)

end subroutine grid_create_reg_cyl

!------------------------------------------------------------------------------------------------------------

subroutine grid_create_cylindr(g, ivdir, nv, nth, v, th, rvth)
!--purpose: Create a cylindr grid from the specified parameters
   implicit none
!--subroutine parameters:
   type(t_grid)              :: g
   integer,      intent(in)  :: ivdir, nv, nth
   real(kind=8), intent(in)  :: v(nv), th(nth), rvth(nv,nth)
!--local variables:
   integer      :: iv, it, ii

   ! check v-direction

   if (ivdir.ne.ikXDIR .and. ivdir.ne.ikYDIR .and. ivdir.ne.ikZDIR) then
      write(bufout,*) 'Internal error: create_reg_cyl: ivdir=',ivdir,' must be X,Y or Z'
      call write_log(1, bufout)
      return
   endif

   ! store the type of grid

   g%is_uniform  = .false.
   g%is_curvilin = .false.
   g%is_cylindr  = .true.
   g%lies_in_oyz = .false.
   g%ivdir       = ivdir

   ! store the number of points in the grid

   g%nx   = nv
   g%ny   = nth
   g%ntot = nv * nth

   ! store dummy grid sizes

   g%dx   = -1d0
   g%dy   = -1d0
   g%dxdy =  1d0

   ! re-allocate coordinate array at the appropriate size

   call reallocate_arr(g%coor, g%ntot, 3)

   ! set pointers to r,th,v-coordinates

   if (ivdir.eq.ikXDIR) then
      g%v  => g%coor(:,1)
      g%r  => g%coor(:,2)
      g%th => g%coor(:,3)
   elseif (ivdir.eq.ikYDIR) then
      g%r  => g%coor(:,1)
      g%v  => g%coor(:,2)
      g%th => g%coor(:,3)
   elseif (ivdir.eq.ikZDIR) then
      g%r  => g%coor(:,1)
      g%th => g%coor(:,2)
      g%v  => g%coor(:,3)
   else
      call write_log(' Internal error: create_cylindr: incorrect ivdir')
   endif
   g%x  => NULL()
   g%y  => NULL()
   g%z  => NULL()

   ! store the (r,th,v)-coordinates of all points of the grid

   do it = 1, nth
      do iv = 1, nv
         ii = iv + (it-1) * nv
         g%r(ii)  = rvth(iv,it)
         g%th(ii) = th(it)
         g%v(ii)  = v(iv)
      enddo
   enddo

   ! clear previous arc-length s and spline data

   if (associated(g%s_prf)) deallocate(g%s_prf)
   g%s_prf  => NULL()

   call grid_destroy_spline(g)

end subroutine grid_create_cylindr

!------------------------------------------------------------------------------------------------------------

subroutine grid_create_curvil(g, nx, ny, x, y, z, lies_in_oyz)
!--purpose: Create a curvilinear grid from the specified parameters
!           note: x may be nx or nx*ny long; y may be ny or nx*ny elements; z may be scalar
   implicit none
!--subroutine parameters:
   type(t_grid)                        :: g
   integer,      intent(in)            :: nx, ny
   real(kind=8), intent(in), optional  :: x(:), y(:), z(:)
   logical,      intent(in), optional  :: lies_in_oyz
!--local variables:
   integer      :: ntot, iy, ii0

   ! store the type of grid

   g%is_uniform  = .false.
   g%is_curvilin = .true.
   g%is_cylindr  = .false.
   g%ivdir       = 0

   if (.not.present(lies_in_oyz)) then
      g%lies_in_oyz = .false.
   else
      g%lies_in_oyz = lies_in_oyz
   endif

   ! store the number of points in the grid

   g%nx   = nx
   g%ny   = ny
   g%ntot = nx * ny

   ! store dummy grid sizes

   g%dx   = -1d0
   g%dy   = -1d0
   g%dxdy =  1d0

   ! re-allocate coordinate array at the appropriate size

   call reallocate_arr(g%coor, g%ntot, 3)

   ! set pointers to x,y,z-coordinates

   g%x  => g%coor(:,1)
   g%y  => g%coor(:,2)
   g%z  => g%coor(:,3)
   g%r  => NULL()
   g%th => NULL()
   g%v  => NULL()

   ! optionally: store the (x,y,z)-coordinates of the grid points

   ntot = g%ntot

   if (.not.g%lies_in_oyz .and. present(x)) then
      if (size(x).eq.1) then
         g%x(1:ntot) = x(1)
      elseif (size(x).eq.nx) then
         do iy = 1, ny
            ii0 = (iy-1) * nx
            g%x(ii0+1:ii0+nx) = x(1:nx)
         enddo
      elseif (size(x).eq.ntot) then
         g%x(1:ntot) = x(1:ntot)
      else
         write(bufout,*) 'grid_create_curvil: Internal error: incorrect length of x=',size(x),', nx=',nx
         call write_log(1, bufout)
      endif
   else
      g%x(1:ntot) = 0d0
   endif

   if (present(y)) then
      if (size(y).eq.1) then
         g%y(1:ntot) = y(1)
      elseif (size(y).eq.ny) then
         do iy = 1, ny
            ii0 = (iy-1) * nx
            g%y(ii0+1:ii0+nx) = y(iy)
         enddo
      elseif (size(y).eq.ntot) then
         g%y(1:ntot) = y(1:ntot)
      else
         write(bufout,*) 'grid_create_curvil: Internal error: incorrect length of y=',size(y),', ny=',ny
         call write_log(1, bufout)
      endif
   else
      g%y(1:ntot) = 0d0
   endif

   if (present(z)) then
      if (size(z).eq.1) then
         g%z(1:ntot) = z(1)
      elseif (size(z).eq.ntot) then
         g%z(1:ntot) = z(1:ntot)
      else
         write(bufout,*) 'grid_create_curvil: Internal error: incorrect length of z=',size(z),', ntot=',ntot
         call write_log(1, bufout)
      endif
   else
      g%z(1:ntot) = 0d0
   endif

   ! clear previous arc-length s and spline data

   if (associated(g%s_prf)) deallocate(g%s_prf)
   g%s_prf  => NULL()

   call grid_destroy_spline(g)

end subroutine grid_create_curvil

!------------------------------------------------------------------------------------------------------------

function grid_has_spline(g)
!--purpose: Tell whether the spline data are filled in or not
   implicit none
!--function result:
   logical grid_has_spline
!--subroutine parameters:
   type(t_grid)              :: g

   grid_has_spline = (associated(g%spl%ayspl))

end function grid_has_spline

!------------------------------------------------------------------------------------------------------------

subroutine grid_alloc_splines(g, np_spl)
!--purpose: (re-)allocate the arrays needed for spline data
   implicit none
!--subroutine parameters:
   type(t_grid)  :: g
   integer       :: np_spl

   ! allocate simple spline arrays, set/update pointers

   call spline_allocate(g%spl, np_spl)

end subroutine grid_alloc_splines

!------------------------------------------------------------------------------------------------------------

subroutine grid_xspline_init(g)
!--purpose: initialize the arrays needed for x-spline data, set x==0
   implicit none
!--subroutine parameters:
   type(t_grid)              :: g

   ! allocate simple spline array axspl(n,4)

   call spline_initx(g%spl)

   ! assert that grid/spline may go outside oyz

   g%lies_in_oyz = .false.

end subroutine grid_xspline_init

!------------------------------------------------------------------------------------------------------------

subroutine grid_make_arclength(g, my_ierror)
!--purpose: compute the arc-length s for a 1-d grid (wheel or rail profile)
!  TODO: add options for length in (x,y), (x,z), (y,z)-planes or (x,y,z)-space
   implicit none
!--subroutine arguments:
   type(t_grid)               :: g
   integer,      intent(out)  :: my_ierror
!--local variables:
   integer            :: ip, ipmin, sub_ierror
   real(kind=8)       :: dy, dz, ds, dst, dstmin, s0

   ! re-allocate s-array at the appropriate size

   my_ierror = 0
   call reallocate_arr(g%s_prf, g%ntot)

   ! compute cumulative sum of distances between successive points

   g%s_prf(1) = 0d0
   do ip = 2, g%ntot
      dy = g%y(ip) - g%y(ip-1)
      dz = g%z(ip) - g%z(ip-1)
      ds = sqrt(dy**2 + dz**2)
      ! write(*,*) 'ip=',ip,': ds=',ds
      g%s_prf(ip) = g%s_prf(ip-1) + ds
   enddo

   ! select reference position in profile where s:=0

   if (is_strictly_monotonic(g%ntot, g%y)) then

      ! interpolate [g%y, g%s_prf] to y = 0

      call interp_1d_to_scalar( g%ntot, g%y, g%s_prf, 0d0, s0, sub_ierror )
      if (my_ierror.eq.0) my_ierror = sub_ierror

   elseif (.false.) then

      ! shifting can be omitted cf. read_profile.m

      s0 = 0d0

   else

      ! select point (g%y, g%z) closest to origin (0, 0)

      ipmin  = -1
      dstmin = 1e12
      do ip = 1, g%ntot
         dst = g%y(ip)**2 + g%z(ip)**2
         if (dst.le.dstmin) then
            ipmin  = ip
            dstmin = dst
         endif
      enddo
      s0 = g%s_prf(ipmin)

   endif

   if (.false.) then
      write(bufout,'(a,f8.3)') ' grid_make_arclength: s0=', s0
      call write_log(1, bufout)
   endif

   ! shift s := s - s0

   do ip = 1, g%ntot
      g%s_prf(ip) = g%s_prf(ip) - s0
   enddo

end subroutine grid_make_arclength

!------------------------------------------------------------------------------------------------------------

subroutine grid_make_ppspline_kink(g, lambda, use_wgt, nkink, ikinks, my_ierror, k_chk)
!--function: compute parametric smoothing spline for 1-d grid (x(s),y(s),z(s))
   implicit none
!--subroutine arguments:
   type(t_grid)                       :: g
   real(kind=8), intent(in)           :: lambda     ! weight of 2nd derivative, relative to wgt_in
   logical,      intent(in)           :: use_wgt    ! weigh data with spacing ds
   integer,      intent(in)           :: nkink
   integer,      intent(in)           :: ikinks(nkink) ! start/end of spline sections
   integer,      intent(out)          :: my_ierror
   integer,      intent(in), optional :: k_chk      ! refinement factor for checking
!--local variables:
   integer      :: sub_ierror
   logical      :: has_xdata
   real(kind=8) :: dist_max

   my_ierror = 0

   ! check requirements: 1-d grid, s-data available

   if (g%nx.ne.1 .and. g%ny.ne.1) then
      write(bufout,*) 'Internal error: make_ppspline works on 1-d grids only: ',g%nx,' x',g%ny
      call write_log(1, bufout)
      call abort_run()
   endif
   if (.not.associated(g%s_prf)) then
      call write_log(' Internal error: make_ppspline: profile s-coordinate has not yet been filled in.')
      call abort_run()
   endif

   ! determine the PP-spline coefficients ax0--ax3, ay0--ay3, az0--az3

   has_xdata = (.not.g%lies_in_oyz)
   call ppspline_make_spline_kink(g%spl, g%ntot, g%s_prf, g%x, g%y, g%z, has_xdata, lambda, use_wgt,    &
                nkink, ikinks, sub_ierror)
   if (sub_ierror.ne.0) my_ierror = sub_ierror

   ! determine max.distance between spline and input data

   if (my_ierror.eq.0 .and. present(k_chk)) then
      call spline_check_updates(g%spl, g%ntot, g%s_prf, g%x, g%y, g%z, k_chk, dist_max, sub_ierror)
      if (sub_ierror.ne.0) my_ierror = sub_ierror
   endif

   ! replace input-data by result of smoothing

   if (my_ierror.eq.0) then
      if (has_xdata) g%x(1:g%ntot) = g%spl%ax0(1:g%spl%npnt)
      g%y(1:g%ntot) = g%spl%ay0(1:g%spl%npnt)
      g%z(1:g%ntot) = g%spl%az0(1:g%spl%npnt)
   endif

end subroutine grid_make_ppspline_kink

!------------------------------------------------------------------------------------------------------------

subroutine grid_make_ppspline_nokink(g, lambda, use_wgt, my_ierror, k_chk)
!--function: compute parametric smoothing spline for 1-d grid (x(s),y(s),z(s))
   implicit none
!--subroutine arguments:
   type(t_grid)                       :: g
   real(kind=8), intent(in)           :: lambda     ! weight of 2nd derivative, using weight 1 for the data
   logical,      intent(in)           :: use_wgt    ! weigh data with spacing ds
   integer,      intent(out)          :: my_ierror
   integer,      intent(in), optional :: k_chk      ! refinement factor for checking
!--local variables:
   integer    :: nkink, ikinks(2)

   nkink  = 2
   ikinks = (/ 1, g%ntot /)

   call grid_make_ppspline_kink(g, lambda, use_wgt, nkink, ikinks, my_ierror, k_chk)

end subroutine grid_make_ppspline_nokink

!------------------------------------------------------------------------------------------------------------

subroutine grid_make_bspline(g, ds_bspl, lambda, use_wgt, nkink, ikinks, naccel, iaccel, my_ierror, k_chk)
!--function: compute parametric smoothing spline for 1-d grid (x(s),y(s),z(s))
   implicit none
!--subroutine arguments:
   type(t_grid)                       :: g
   real(kind=8), intent(in)           :: ds_bspl    ! target stepsize
   real(kind=8), intent(in)           :: lambda     ! weight of 2nd derivative, relative to wgt_in
   logical,      intent(in)           :: use_wgt    ! weigh data with spacing ds
   integer,      intent(in)           :: nkink
   integer,      intent(in)           :: ikinks(nkink)  ! start/end, jumps in 1st derivative
   integer,      intent(in)           :: naccel
   integer,      intent(in)           :: iaccel(naccel) ! jumps in radius of curvature
   integer,      intent(out)          :: my_ierror
   integer,      intent(in), optional :: k_chk      ! refinement factor for checking
!--local variables:
   integer      :: sub_ierror
   logical      :: has_xdata
   real(kind=8) :: dist_max

   my_ierror = 0

   ! check requirements: 1-d grid, s-data available

   if (g%nx.ne.1 .and. g%ny.ne.1) then
      write(bufout,*) 'Internal error: make_bspline works on 1-d grids only: ',g%nx,' x',g%ny
      call write_log(1, bufout)
      call abort_run()
   endif
   if (.not.associated(g%s_prf)) then
      call write_log(' Internal error: make_bspline: profile s-coordinate has not yet been filled in.')
      call abort_run()
   endif

   ! determine the PP-spline coefficients ax0--ax3, ay0--ay3, az0--az3

   has_xdata = (.not.g%lies_in_oyz)
   call bspline_make_1d_ppspline(g%spl, g%ntot, g%s_prf, g%x, g%y, g%z, has_xdata, ds_bspl, lambda,     &
                use_wgt, nkink, ikinks, naccel, iaccel, sub_ierror)
   if (sub_ierror.ne.0) my_ierror = sub_ierror

   ! determine max.distance between spline and input data

   if (my_ierror.eq.0 .and. present(k_chk)) then
      call spline_check_updates(g%spl, g%ntot, g%s_prf, g%x, g%y, g%z, k_chk, dist_max, sub_ierror)
      if (sub_ierror.ne.0) my_ierror = sub_ierror
   endif

   ! replace input-data by spline evaluated at profile s-positions

   if (my_ierror.eq.0) then
      if (g%spl%has_xdata) then
         call spline_eval(g%spl, ikXDIR, g%ntot, g%s_prf, sub_ierror, f_eval=g%x)
         if (my_ierror.eq.0) my_ierror = sub_ierror
      endif

      call spline_eval(g%spl, ikYDIR, g%ntot, g%s_prf, sub_ierror, f_eval=g%y)
      if (my_ierror.eq.0) my_ierror = sub_ierror

      call spline_eval(g%spl, ikZDIR, g%ntot, g%s_prf, sub_ierror, f_eval=g%z)
      if (my_ierror.eq.0) my_ierror = sub_ierror
   endif

end subroutine grid_make_bspline

!------------------------------------------------------------------------------------------------------------

subroutine grid_extrude_profile_xarr(gprof, nx, x_in, gsurf)
!--purpose: Construct a prismatic grid specified by list of {xi} and lateral profile gprof = {yj,zj}
   implicit none
!--subroutine parameters:
   type(t_grid), intent(in)    :: gprof
   integer,      intent(in)    :: nx
   real(kind=8), intent(in)    :: x_in(nx)
   type(t_grid), intent(inout) :: gsurf
!--local variables:
   integer      :: ix, iy, ii

   if (gprof%nx.ne.1) then
      write(bufout,*) 'Internal error: grid_extrude requires nx=1 (',nx,')'
      call write_log(1, bufout)
      call abort_run()
   endif

   ! create output grid with nx slices, ny profile points

   call grid_create_curvil(gsurf, nx, gprof%ny, lies_in_oyz=.false.)

   ! fill grid 

   do iy = 1, gprof%ny
      do ix = 1, nx
         ii    = ix + (iy-1) * nx
         gsurf%x(ii) = x_in(ix)
         gsurf%y(ii) = gprof%y(iy)
         gsurf%z(ii) = gprof%z(iy)
      enddo
   enddo

end subroutine grid_extrude_profile_xarr

!------------------------------------------------------------------------------------------------------------

subroutine grid_extrude_profile_xunif(gprof, nx, x0, dx, gsurf)
!--purpose: Construct a prismatic grid specified by uniform {x0:dx:x1} and lateral profile gprof = {yj,zj}
   implicit none
!--subroutine parameters:
   type(t_grid), intent(in)    :: gprof
   integer,      intent(in)    :: nx
   real(kind=8), intent(in)    :: x0, dx
   type(t_grid), intent(inout) :: gsurf
!--local variables:
   integer    :: ix
   real(kind=8), dimension(:), allocatable :: xarr

   allocate(xarr(nx))
   do ix = 1, nx
      xarr(ix) = x0 + (ix-1) * dx
   enddo

   call grid_extrude_profile_xarr(gprof, nx, xarr, gsurf)
   deallocate(xarr)

end subroutine grid_extrude_profile_xunif

!------------------------------------------------------------------------------------------------------------

subroutine grid_revolve_profile_xarr(gprof, nx, ny, x_in, z_axle, defval, gsurf)
!--purpose: Consider a body of revolution specified by gprof = {yj,zj} and z_axle, 
!           create curvilinear grid gsurf sampling the surface at x-positions {xi} x {yj}.
!           {xi} can be 1d list or 2d grid.
   implicit none
!--subroutine parameters:
   type(t_grid), intent(in)    :: gprof
   integer,      intent(in)    :: nx, ny
   real(kind=8), intent(in)    :: x_in(nx*ny)
   real(kind=8), intent(in)    :: z_axle, defval
   type(t_grid), intent(inout) :: gsurf
!--local variables:
   integer      :: ix, iy, ii
   real(kind=8) :: nom_radius, ry, sgn

   if (gprof%nx.ne.1) then
      write(bufout,*) 'ERROR: grid_revolve requires nx=1 (',nx,')'
      call write_log(1, bufout)
      call abort_run()
   endif
   if (ny.gt.1 .and. ny.ne.gprof%ny) then
      write(bufout,*) 'ERROR: grid_revolve requires ny=1 or gin%ny=', gprof%ny,' (',ny,')'
      call write_log(1, bufout)
      call abort_run()
   endif

   ! create output grid with nx slices, ny profile points

   call grid_create_curvil(gsurf, nx, gprof%ny, lies_in_oyz=.false.)

   ! fill grid 

   nom_radius = abs(z_axle)

   if (z_axle.lt.0d0) then
      sgn =  1d0
   else
      sgn = -1d0
   endif

   ! take positive sqrt for the wheel (or concave rail),
   !      negative sqrt for roller (or inside of hollow wheel)

   associate(x => gsurf%x, y => gsurf%y, z => gsurf%z)

   do iy = 1, gprof%ny
      do ix = 1, nx
         ii    = ix + (iy-1) * nx
         if (ny.le.1) then
            x(ii) = x_in(ix)
         else
            x(ii) = x_in(ii)
         endif
         y(ii) =     gprof%y(iy)
         ry    = abs(gprof%z(iy) - z_axle)

         if (ry**2-x(ii)**2.gt.0d0) then
            z(ii) =  z_axle + sgn*sqrt( ry**2 - x(ii)**2 )
         else
            z(ii) = defval
         endif
      enddo
   enddo

   end associate

end subroutine grid_revolve_profile_xarr

!------------------------------------------------------------------------------------------------------------

subroutine grid_revolve_profile_xarr_zarr(gprof, nx, ny, x_in, z_axle, defval, z_out)
!--purpose: Consider a body of revolution specified by gprof = {yj,zj} and z_axle, 
!           create heights 'z_out' sampling the surface at x-positions {xi} x {yj}.
!           {xi} can be 1d list or 2d grid.
   implicit none
!--subroutine parameters:
   type(t_grid)               :: gprof
   integer,      intent(in)   :: nx, ny
   real(kind=8), intent(in)   :: x_in(nx*ny)
   real(kind=8), intent(in)   :: z_axle, defval
   real(kind=8)               :: z_out(:)               ! note: z_out may be gprof%z if nx==1
!--local variables:
   integer      :: ix, iy, ii
   real(kind=8) :: nom_radius, ry, sgn, x, y

   if (gprof%nx.ne.1) then
      write(bufout,*) 'ERROR: grid_revolve requires nx=1 (',nx,')'
      call write_log(1, bufout)
      call abort_run()
   endif
   if (ny.gt.1 .and. ny.ne.gprof%ny) then
      write(bufout,*) 'ERROR: grid_revolve requires ny=1 or gin%ny=', gprof%ny,' (',ny,')'
      call write_log(1, bufout)
      call abort_run()
   endif

   ! check size of output-array: nx slices, ny profile points

   if (size(z_out,1).lt.nx*gprof%ny) then
      write(bufout,*) 'ERROR: grid_revolve requires output-array of nx*ny=', nx, gprof%ny,      &
                ' (',size(z_out,1),')'
      call write_log(1, bufout)
      call abort_run()
   endif

   ! fill grid 

   nom_radius = abs(z_axle)

   if (z_axle.lt.0d0) then
      sgn =  1d0
   else
      sgn = -1d0
   endif

   ! take positive sqrt for the wheel (or concave rail),
   !      negative sqrt for roller (or inside of hollow wheel)

   do iy = 1, gprof%ny
      do ix = 1, nx
         ii    = ix + (iy-1) * nx
         if (ny.le.1) then
            x = x_in(ix)
         else
            x = x_in(ii)
         endif
         y  =     gprof%y(iy)
         ry = abs(gprof%z(iy) - z_axle)

         if (ry**2-x**2.gt.0d0) then
            z_out(ii) =  z_axle + sgn*sqrt( ry**2 - x**2 )
         else
            z_out(ii) = defval
         endif
      enddo
   enddo

end subroutine grid_revolve_profile_xarr_zarr

!------------------------------------------------------------------------------------------------------------

subroutine grid_revolve_profile_xunif(gprof, nx, x0, dx, z_axle, defval, gsurf)
!--purpose: Consider a body of revolution specified by gprof = {yj,zj} and nom_radius, 
!           create curvilinear grid gsurf sampling the surface at x-positions {xi} x {yj}.
   implicit none
!--subroutine parameters:
   type(t_grid), intent(in)    :: gprof
   integer,      intent(in)    :: nx
   real(kind=8), intent(in)    :: x0, dx
   real(kind=8), intent(in)    :: z_axle, defval
   type(t_grid), intent(inout) :: gsurf
!--local variables:
   integer    :: ix
   real(kind=8), dimension(:), allocatable :: xarr

   allocate(xarr(nx))
   do ix = 1, nx
      xarr(ix) = x0 + (ix-1) * dx
   enddo

   call grid_revolve_profile_xarr(gprof, nx, 1, xarr, z_axle, defval, gsurf)
   deallocate(xarr)

end subroutine grid_revolve_profile_xunif

!------------------------------------------------------------------------------------------------------------

subroutine grid_spline_copy(gin, gout)
!--purpose: copy spline data from input grid gin to output grid gout
   implicit none
!--subroutine parameters:
   type(t_grid)          :: gin, gout

   ! copy the spline data

   call spline_copy(gin%spl, gout%spl)

end subroutine grid_spline_copy

!------------------------------------------------------------------------------------------------------------

subroutine grid_copy(gin, gout, with_spline)
!--purpose: copy input grid gin to output grid gout
   implicit none
!--subroutine parameters:
   type(t_grid)          :: gin, gout
   logical,     optional :: with_spline
!--local variables:
   logical      :: copy_spline
   integer      :: it

   copy_spline = .false.
   if (present(with_spline)) copy_spline = with_spline

   ! copy the type of grid

   gout%is_uniform  = gin%is_uniform
   gout%is_curvilin = gin%is_curvilin
   gout%is_cylindr  = gin%is_cylindr
   gout%lies_in_oyz = gin%lies_in_oyz
   gout%ivdir       = gin%ivdir

   ! copy the number of points in the grid

   gout%nx   = gin%nx
   gout%ny   = gin%ny
   gout%ntot = gin%ntot

   ! copy grid sizes

   gout%dx   = gin%dx
   gout%dy   = gin%dy
   gout%dxdy = gin%dxdy

   ! allocate coordinate array at the appropriate size

   call reallocate_arr(gout%coor, gout%ntot, 3)

   ! set pointers to x,y,z- or r,th,z-coordinates

   if (gout%is_cylindr) then
      gout%x  => NULL()
      gout%y  => NULL()
      gout%z  => NULL()
      if (gout%ivdir.eq.ikXDIR) then
         gout%v  => gout%coor(:,1)
         gout%r  => gout%coor(:,2)
         gout%th => gout%coor(:,3)
      elseif (gout%ivdir.eq.ikYDIR) then
         gout%r  => gout%coor(:,1)
         gout%v  => gout%coor(:,2)
         gout%th => gout%coor(:,3)
      elseif (gout%ivdir.eq.ikZDIR) then
         gout%r  => gout%coor(:,1)
         gout%th => gout%coor(:,2)
         gout%v  => gout%coor(:,3)
      else
         call write_log(' Internal error: grid_copy: incorrect ivdir')
      endif
   else
      gout%x  => gout%coor(:,1)
      gout%y  => gout%coor(:,2)
      gout%z  => gout%coor(:,3)
      gout%r  => NULL()
      gout%th => NULL()
      gout%v  => NULL()
   endif

   ! copy the coordinates of the grid points

   associate(ntot => gout%ntot)
   do it = 1, 3
      gout%coor(1:ntot,it) = gin%coor(1:ntot,it)
   enddo

   ! if s-array is present in input grid: allocate in gout, copy s-coordinates, else: clear s-coordinates

   if (associated(gin%s_prf)) then
      call reallocate_arr(gout%s_prf, ntot)
      gout%s_prf(1:ntot) = gin%s_prf(1:ntot)
   elseif (associated(gout%s_prf)) then
      deallocate(gout%s_prf)
      gout%s_prf  => NULL()
   endif
   end associate

   ! if spline-data are present in input grid: copy spline-data, else: clear spline-data

   if (copy_spline .and. grid_has_spline(gin)) then
      call grid_spline_copy(gin, gout)
   else
      call grid_destroy_spline(gout)
   endif

end subroutine grid_copy

!------------------------------------------------------------------------------------------------------------

subroutine grid_trim(gin, gout, ix_low_arg, ix_hig_arg, iy_low_arg, iy_hig_arg, s_low, s_hig,           &
                with_spline, idebug_arg)
!--purpose: create trimmed version of input grid gin in output grid gout
   implicit none
!--subroutine parameters:
   type(t_grid)           :: gin, gout
   integer                :: ix_low_arg, ix_hig_arg
   integer,      optional :: iy_low_arg, iy_hig_arg, idebug_arg
   real(kind=8), optional :: s_low, s_hig
   logical,      optional :: with_spline
!--local variables:
   logical      :: copy_s, copy_spline
   integer      :: ix_low, ix_hig, iy_low, iy_hig, idebug
   integer      :: ix_of, ix, iy, ii_in, ii_out

   idebug = 0
   if (present(idebug_arg)) idebug = idebug_arg

   copy_spline = .true. ! default
   if (present(with_spline)) copy_spline = with_spline

   copy_s      = associated(gin%s_prf)

   if (present(iy_low_arg) .and. present(iy_low_arg)) then

      iy_low = iy_low_arg
      iy_hig = iy_hig_arg

   elseif (present(s_low) .and. present(s_hig) .and. copy_s) then

      ! determine the indices to keep

      call locate_interval(gin%ntot, gin%s_prf, s_low, s_hig, iy_low, iy_hig)

      if (idebug.ge.1) then
         write(bufout,'(2(a,f8.3),2(a,i4),a)') ' grid_trim: s=[',s_low,',',s_hig,']: keeping iy=[',     &
                iy_low, ',',iy_hig,']'
         call write_log(1, bufout)
      endif

   else

      call write_log(' Internal error(grid_trim): either iy_low or s_low must be given.')
      if (present(iy_low_arg)) call write_log('iy_low_arg is present')
      if (present(iy_hig_arg)) call write_log('iy_hig_arg is present')
      if (present(s_low)) call write_log('s_low is present')
      if (present(s_hig)) call write_log('s_hig is present')
      if (copy_s) call write_log('array s_prf is present')
      call abort_run()

   endif

   ix_low = max(     1, ix_low_arg)
   ix_hig = min(gin%nx, ix_hig_arg)
   iy_low = max(     1, iy_low)
   iy_hig = min(gin%ny, iy_hig)

   ! copy the type of grid

   gout%is_uniform  = gin%is_uniform
   gout%is_curvilin = gin%is_curvilin
   gout%is_cylindr  = gin%is_cylindr
   gout%lies_in_oyz = gin%lies_in_oyz
   gout%ivdir       = gin%ivdir

   ! set the number of points in the output grid

   gout%nx   = ix_hig - ix_low + 1
   gout%ny   = iy_hig - iy_low + 1
   gout%ntot = gout%nx * gout%ny

   ! copy grid sizes

   gout%dx   = gin%dx
   gout%dy   = gin%dy
   gout%dxdy = gin%dxdy

   ! allocate coordinate array at the appropriate size

   call reallocate_arr(gout%coor, gout%ntot, 3)

   if (idebug.ge.3) then
      write(bufout,'(4(a,i6))') ' grid_trim: output grid has',gout%nx,' x',gout%ny,' points, coor=',   &
           size(gout%coor,1),' x', size(gout%coor,2)
      call write_log(1, bufout)
   endif

   ! set pointers to x,y,z- or r,th,z-coordinates

   if (gout%is_cylindr) then
      gout%x  => NULL()
      gout%y  => NULL()
      gout%z  => NULL()
      if (gout%ivdir.eq.ikXDIR) then
         gout%v  => gout%coor(:,1)
         gout%r  => gout%coor(:,2)
         gout%th => gout%coor(:,3)
      elseif (gout%ivdir.eq.ikYDIR) then
         gout%r  => gout%coor(:,1)
         gout%v  => gout%coor(:,2)
         gout%th => gout%coor(:,3)
      elseif (gout%ivdir.eq.ikZDIR) then
         gout%r  => gout%coor(:,1)
         gout%th => gout%coor(:,2)
         gout%v  => gout%coor(:,3)
      else
         call write_log(' Internal error: grid_trim: incorrect ivdir')
      endif
   else
      gout%x  => gout%coor(:,1)
      gout%y  => gout%coor(:,2)
      gout%z  => gout%coor(:,3)
      gout%r  => NULL()
      gout%th => NULL()
      gout%v  => NULL()
   endif

   ! allocate s-array if present in input

   if (copy_s) then
      call reallocate_arr(gout%s_prf, gout%ntot)
   else
      if (associated(gout%s_prf)) deallocate(gout%s_prf)
      gout%s_prf => NULL()
   endif

   ! copy the coordinates of the selected grid points

   ix_of = ix_low - 1

   do iy = iy_low, iy_hig
      do ix = ix_low, ix_hig
         ii_in  =  ix        + (iy-1)     *gin%nx
         ii_out = (ix-ix_of) + (iy-iy_low)*gout%nx
         gout%coor(ii_out,1:3) = gin%coor(ii_in,1:3)
         if (copy_s) gout%s_prf(ii_out) = gin%s_prf(ii_in)
      enddo
   enddo

   ! if spline-data are present in input grid: copy spline-data, else: clear spline-data

   if (copy_spline .and. grid_has_spline(gin)) then
      if     (gin%ntot.eq.gin%spl%npnt .and. gout%nx.eq.1) then
         call spline_trim(gin%spl, gout%spl, iy_low, iy_hig)
      elseif (gin%ntot.eq.gin%spl%npnt .and. gout%ny.eq.1) then
         call spline_trim(gin%spl, gout%spl, ix_low, ix_hig)
      else
         call spline_trim(gin%spl, gout%spl, s_low=gout%s_prf(1), s_hig=gout%s_prf(gout%ntot))
      endif
   else
      call grid_destroy_spline(gout)
   endif

end subroutine grid_trim

!------------------------------------------------------------------------------------------------------------

subroutine grid_print(g, nam, idebug, ndigit)
!--function: print information on grid g
   implicit none
!--subroutine arguments
   type(t_grid)         :: g
   character(len=*)     :: nam
   integer              :: idebug
   integer, optional    :: ndigit       ! number of significant digits
!--local variables
   character(len=1), parameter  :: cnams(1:3) = (/ 'x', 'y', 'z' /)
   integer              :: my_ndigit, my_len
   integer              :: ix, iy, ii, j, nval, nline
   character(len=18)    :: strng(5)

   if (present(ndigit)) then
      my_ndigit = ndigit
   else
      my_ndigit = 6
   endif
   my_ndigit = max(2, min(10, my_ndigit))
   my_len    = 8 + my_ndigit

   if (.not.associated(g%coor)) then
      write(bufout,'(3a)') ' grid ',trim(nam),' has not yet been defined'
      call write_log(1, bufout)
      return
   endif

   ! idebug>=1: overview: type of grid, number of points

   if (idebug.ge.1) then
      if (g%is_uniform) then
         write(bufout,'(3a,2(i4,a))') ' grid ',trim(nam),' is uniform with', g%nx,' x',g%ny,' points'
      elseif (g%is_curvilin) then
         if (g%nx.eq.1) then
            write(bufout,'(3a,i6,a)') ' profile "',trim(nam),'" has', g%ny,' points'
         else
            write(bufout,'(3a,2(i4,a))') ' grid ',trim(nam),' is curvilinear with', g%nx,' x',g%ny,' points'
         endif
      elseif (g%is_cylindr) then
         write(bufout,'(3a,2(i4,a))') ' grid ',trim(nam),' is cylindrical with', g%nx,' x',g%ny,' points'
         call write_log(1, bufout)
         write(bufout,'(3a)') ' the cylinder axis is aligned with the ',cnams(g%ivdir),'-axis'
      endif
      call write_log(1, bufout)
      if (g%lies_in_oyz) then
         write(bufout,'(3a)') ' grid ',trim(nam),' lies in the oyz plane, x==0'
         call write_log(1, bufout)
      endif
   endif

   ! idebug>=2: overall details: structured grid coordinates

   if (idebug.ge.2) then
      if (g%is_uniform) then
         write(bufout,121) '  x = [', g%x(1), ' :', g%dx, ' :', g%x(g%ntot), ']'
         call write_log(1, bufout)
         write(bufout,121) '  y = [', g%y(1), ' :', g%dy, ' :', g%y(g%ntot), ']'
         call write_log(1, bufout)
 121     format(a, 3(f11.6, a))
      elseif (g%is_cylindr) then
         nline = min(10, int( (g%nx-1)/10 ) + 1)
         write(bufout,122) ' v  = ', (g%v(ii), ii=1, min(100,g%nx))
         call write_log(nline, bufout)

         nline = min(10, int( (g%ny-1)/10 ) + 1)
         write(bufout,122) ' th = ', (g%th(ii), ii=1, min(100*g%nx,g%ntot), g%nx)
         call write_log(nline, bufout)
 122     format(a, 10(f9.3,:,','), /, 19( 6x, 10(f9.3,:,','), /) )
      endif
   endif

   ! idebug>=6: array sizes

   if (idebug.ge.6) then
      call print_array_size_2d(g%coor, 'coor')
      call print_array_size_1d(g%x, 'x')
      call print_array_size_1d(g%y, 'y')
      call print_array_size_1d(g%z, 'z')
      call print_array_size_1d(g%r, 'r')
      call print_array_size_1d(g%th, 'th')
      call print_array_size_1d(g%v, 'v')
      call print_array_size_1d(g%s_prf, 's_prf')
   endif

   ! idebug>=4: full details: full list of coordinates

   if (idebug.ge.4) then
 
      nval = 3
      if (associated(g%s_prf)) nval = 4

      write(bufout,'(2a)') trim(nam), ' = ['
      call write_log(1, bufout)

      if (g%is_cylindr) then
         strng(1) = '(  iv, ith)'
         strng(2) = '         r'
         strng(3) = '        th'
         strng(4) = '         r'
      else
         strng(1) = '(  ix,  iy)'
         strng(2) = '         x'
         strng(3) = '         y'
         strng(4) = '         z'
      endif
      strng(5) = '         s'
      write(bufout, 210) strng(1)(1:11), (strng(j+1)(1:my_len), j=1,nval)
      call write_log(1, bufout)
 210  format(' %      ii ',5(a,1x))

      do iy = 1, g%ny
         do ix = 1, g%nx
            ii = ix + (iy-1)*g%nx

            if (g%is_cylindr) then
               strng(1) = fmt_gs(my_len, my_ndigit, g%r(ii))
               strng(2) = fmt_gs(my_len, my_ndigit, g%th(ii))
               strng(3) = fmt_gs(my_len, my_ndigit, g%v(ii))
            else
               strng(1) = fmt_gs(my_len, my_ndigit, g%x(ii))
               strng(2) = fmt_gs(my_len, my_ndigit, g%y(ii))
               strng(3) = fmt_gs(my_len, my_ndigit, g%z(ii))
            endif
            if (nval.ge.4) strng(4) = fmt_gs(my_len, my_ndigit, g%s_prf(ii))

            write(bufout,211) ii, ix, iy, (strng(j)(1:my_len), j=1,nval)
 211        format(' ii=',i6,' (',i4,',',i4,'):',4(a,1x))
            call write_log(1, bufout)

         enddo
      enddo

      call write_log('];')
   endif

   ! idebug>=1: overview w.r.t. spline

   if (idebug.ge.1 .and. grid_has_spline(g)) then
      write(bufout,'(3a,i6,a)') ' grid "',trim(nam),'" has spline data,', g%spl%npnt,' points'
      call write_log(1, bufout)
   endif

   ! idebug>=3: spline top-view, idebug>=5: print spline data

   call spline_print(g%spl, nam, idebug, ndigit)

end subroutine grid_print

!------------------------------------------------------------------------------------------------------------

subroutine grid_destroy(g)
!--purpose: clean-up allocated array, nullify alias-pointers
   implicit none
!--subroutine parameters:
   type(t_grid)              :: g

   if (associated(g%coor)) deallocate(g%coor)
   g%coor => NULL()
   g%x    => NULL()
   g%y    => NULL()
   g%z    => NULL()
   g%r    => NULL()
   g%th   => NULL()
   g%v    => NULL()

   if (associated(g%s_prf))    deallocate(g%s_prf)
   g%s_prf    => NULL()

   call grid_destroy_spline(g)

end subroutine grid_destroy

!------------------------------------------------------------------------------------------------------------

subroutine grid_destroy_spline(g)
!--purpose: clean-up allocated arrays for splines, nullify alias-pointers
   implicit none
!--subroutine parameters:
   type(t_grid)              :: g

   call spline_destroy(g%spl)

end subroutine grid_destroy_spline

!------------------------------------------------------------------------------------------------------------

subroutine unifgrid_mirror_y(g_in, g_mirr)
!--purpose: create grid with y-coordinates mirrored and in reversed order
   implicit none
!--subroutine parameters:
   type(t_grid)              :: g_in, g_mirr
!--local variables:
   integer      :: nx, ny, ix, iy, ii
   real(kind=8) :: y0_in, y0_mirr, yy

   if (.not.g_in%is_uniform) then
      call write_log(' Internal error: unifgrid_mirror_y implemented for uniform grids only.')
      call abort_run()
   endif

   ! copy type of grid, #points, grid sizes; copy original coordinate values; set pointers

   call grid_copy(g_in, g_mirr)

   ! mirror the y coordinates of all points of the grid and store in increasing order

   nx      = g_in%nx
   ny      = g_in%ny
   y0_in   = g_in%y(1)
   y0_mirr = -(y0_in + (ny-1)*g_in%dy)

   do iy = 1, ny
      yy = y0_mirr + (iy-1) * g_in%dy
      do ix = 1, nx
         ii = ix + (iy-1) * nx
         g_mirr%y(ii) = yy
      enddo
   enddo

   ! clear previous arc-length s and spline data

   if (associated(g_mirr%s_prf)) deallocate(g_mirr%s_prf)
   g_mirr%s_prf  => NULL()

   call grid_destroy_spline(g_mirr)

end subroutine unifgrid_mirror_y

!------------------------------------------------------------------------------------------------------------

subroutine cartgrid_shift(g, dx, dy, dz)
!--function: shift a 2d surface in cartesian coordinates: g = g + [dx; dy; dz]
   implicit none
!--subroutine arguments
   type(t_grid)    :: g
   real(kind=8)    :: dx, dy, dz
!--local variables
   integer         :: ntot

   if (g%is_cylindr) then
      call write_log(' Internal error: cartgrid_shift: g is in cylindrical coordinates')
   else
      ntot   = g%ntot
      g%x(1:ntot) = g%x(1:ntot) + dx
      g%y(1:ntot) = g%y(1:ntot) + dy
      g%z(1:ntot) = g%z(1:ntot) + dz
      if (g%lies_in_oyz .and. abs(dx).gt.1d-9) g%lies_in_oyz = .false.
   endif

   ! apply shift to spline coefficients d

   if (grid_has_spline(g)) call spline_shift(g%spl, dx, dy, dz)

end subroutine cartgrid_shift

!------------------------------------------------------------------------------------------------------------

subroutine cartgrid_mirror_y(g)
!--function: mirror a 2d surface in cartesian coordinates wrt plane Oxz: g.y = -g.y
!            Note: no reordering of grid points here.
   implicit none
!--subroutine arguments
   type(t_grid)    :: g
!--local variables
   integer         :: ntot

   if (g%is_cylindr) then
      call write_log(' Internal error: cartgrid_mirror_y: g is in cylindrical coordinates')
   else
      ntot   = g%ntot
      g%y(1:ntot) = - g%y(1:ntot)
   endif

   ! apply mirroring to spline coefficients ay0, ay1, ay3

   if (grid_has_spline(g)) call spline_mirror_y(g%spl)

end subroutine cartgrid_mirror_y

!------------------------------------------------------------------------------------------------------------

subroutine cartgrid_rotate(g, rot, xc, yc, zc)
!--function: rotate a 2d surface in cartesian coordinates by rotation matrix about center-point [xc;yc;zc]
   implicit none
!--subroutine arguments
   type(t_grid)             :: g            ! input/output surface
   type(t_rotmat)           :: rot          ! matrix with new orientation of original unit vectors
   real(kind=8), intent(in) :: xc, yc, zc   ! rotation origin
!--local variables
   integer      :: ii
   real(kind=8) :: xrel, yrel, zrel

   if (g%is_cylindr) then
      call write_log(' Internal error: cartgrid_rotate: g is in cylindrical coordinates')
   else
      do ii = 1, g%ntot
         xrel = g%x(ii) - xc
         yrel = g%y(ii) - yc
         zrel = g%z(ii) - zc
         g%x(ii) = xc + rot%r(1) * xrel + rot%r(4) * yrel + rot%r(7) * zrel
         g%y(ii) = yc + rot%r(2) * xrel + rot%r(5) * yrel + rot%r(8) * zrel
         g%z(ii) = zc + rot%r(3) * xrel + rot%r(6) * yrel + rot%r(9) * zrel
      end do
      if (g%lies_in_oyz .and. (abs(rot%r(4)).gt.1d-9 .or. abs(rot%r(7)).gt.1d-9)) g%lies_in_oyz = .false.
   endif

   ! apply rotation to spline coefficients {a3, a2, a1}, without origin, and to {a0}, relative to (yc,zc)

   if (grid_has_spline(g)) call spline_rotate(g%spl, rot, xc, yc, zc)

end subroutine cartgrid_rotate

!------------------------------------------------------------------------------------------------------------

subroutine cartgrid_roll(g, roll, yc, zc)
!--function: rotate a 2d surface in cartesian coordinates by roll angle roll [rad] (about x-axis/point [yc;zc])
   implicit none
!--subroutine arguments
   type(t_grid)             :: g
   real(kind=8), intent(in) :: roll     ! rotation angle [rad]
   real(kind=8), intent(in) :: yc, zc   ! rotation origin
!--local variables
   integer      :: ii
   real(kind=8) :: cs, sn, yrel, zrel

   if (g%is_cylindr) then
      call write_log(' Internal error: cartgrid_roll: g is in cylindrical coordinates')
   else
      cs = cos(roll)
      sn = sin(roll)

      do ii = 1, g%ntot
         yrel = g%y(ii) - yc
         zrel = g%z(ii) - zc
         g%y(ii) = yc + cs * yrel - sn * zrel
         g%z(ii) = zc + sn * yrel + cs * zrel
      end do
   endif

   ! apply rotation to spline coefficients {a3, a2, a1}, without origin, and to {a0}, relative to (yc,zc)

   if (grid_has_spline(g)) call spline_roll(g%spl, roll, yc, zc)

end subroutine cartgrid_roll

!------------------------------------------------------------------------------------------------------------

subroutine cartgrid_yaw(g, yaw, xc, yc)
!--function: rotate a 2d surface in cartesian coordinates by pure yaw angle yaw [rad] (about z-axis,
!            centered at point [xc;yc])
   implicit none
!--subroutine arguments
   type(t_grid)             :: g
   real(kind=8), intent(in) :: yaw      ! rotation angle [rad]
   real(kind=8), intent(in) :: xc, yc   ! rotation origin
!--local variables
   integer      :: ii
   real(kind=8) :: cs, sn, xrel, yrel

   if (g%is_cylindr) then
      call write_log(' Internal error: cartgrid_yaw: g is in cylindrical coordinates')
   else
      cs = cos(yaw)
      sn = sin(yaw)

      do ii = 1, g%ntot
         xrel = g%x(ii) - xc
         yrel = g%y(ii) - yc
         g%x(ii) = xc + cs * xrel - sn * yrel
         g%y(ii) = yc + sn * xrel + cs * yrel
      end do
   endif

   ! note: no action for spline coefficients (?)

   if (grid_has_spline(g) .and. g%lies_in_oyz) then
      call write_log(' Internal warning: cartgrid_yaw: spline coefficients not rotated')
   elseif (grid_has_spline(g)) then
      call write_log(' Internal error: cartgrid_yaw: spline with x-values not supported')
   endif

end subroutine cartgrid_yaw

!------------------------------------------------------------------------------------------------------------

subroutine cartgrid_roll_yaw(g, roll, yaw, xc, yc, zc)
!--function: rotate a 2d surface in cartesian coordinates by roll-and-yaw [rad] about center-point [xc;yc;zc]
!            roll is about original x-axis, yaw about modified z-axis
   implicit none
!--subroutine arguments
   type(t_grid)             :: g
   real(kind=8), intent(in) :: roll, yaw    ! rotation angle [rad]
   real(kind=8), intent(in) :: xc, yc, zc   ! rotation origin
!--local variables
   type(t_rotmat) :: rot

   rot = rotmat_roll_yaw(roll, yaw)
   call cartgrid_rotate(g, rot, xc, yc, zc)

end subroutine cartgrid_roll_yaw

!------------------------------------------------------------------------------------------------------------

subroutine cartgrid_2glob_or(g, o, R)
!--function: compute local-to-global conversion for a surface defined with respect to local system (o, R)
!            o = origin of local system w.r.t. global system
!            R = orientation of local system w.r.t. global system
!            after the transformation, g is defined with respect to the global system
   implicit none
!--subroutine arguments
   type(t_grid),   intent(inout) :: g
   type(t_vec),    intent(in)    :: o
   type(t_rotmat), intent(in)    :: R

   ! rotate points with respect to initial origin o

   call cartgrid_rotate(g, R, 0d0, 0d0, 0d0)

   ! change origin to the global system o

   call cartgrid_shift(g, o%v(1), o%v(2), o%v(3))

end subroutine cartgrid_2glob_or

!------------------------------------------------------------------------------------------------------------

subroutine cartgrid_2glob_m(g, mref)
!--function: compute local-to-global conversion for a surface defined with respect to mref
!            o = mref%o   == origin of reference w.r.t. global system
!            R = mref%rot == orientation of reference w.r.t. global system
!            after the transformation, g is defined with respect to the global system
   implicit none
!--subroutine arguments
   type(t_marker), intent(in)    :: mref
   type(t_grid),   intent(inout) :: g

   ! rotate points with respect to initial origin mref%o

   call cartgrid_rotate(g, mref%rot, 0d0, 0d0, 0d0)

   ! change origin to the global system o

   call cartgrid_shift(g, mref%o%v(1), mref%o%v(2), mref%o%v(3))

end subroutine cartgrid_2glob_m

!------------------------------------------------------------------------------------------------------------

subroutine cartgrid_2loc_m(g, mref)
!--function: convert surface g defined in global coordinates to local coordinates according to mref
!            o = mref%o   == origin of new reference w.r.t. global system
!            R = mref%rot == orientation of new reference w.r.t. global system
!            after the transformation, g is defined with respect to the local system
   implicit none
!--subroutine arguments
   type(t_marker), intent(in)    :: mref
   type(t_grid),   intent(inout) :: g
!--local variables
   type(t_marker)              :: mglb_ref

   ! compute transpose of mref, i.e. the marker for the global system in terms of the local reference

   mglb_ref = marker_transpose(mref)

   ! transform g from the new local system 'glb' to the new global system 'ref'

   call cartgrid_2glob_m(g, mglb_ref)

end subroutine cartgrid_2loc_m

!------------------------------------------------------------------------------------------------------------

subroutine convert_cyl2curv(g)
!--function: compute curvilinear cartesian coordinates (x,y,z) for cylindrical coordinates (r,th,v)
   implicit none
!--subroutine arguments:
   type(t_grid)              :: g                 ! in: grid with cylindrical (r,th,v) filled in
                                                  ! out: grid with curvilinear (x,y,z) filled in
!--local variables:
!  character(len=*), parameter  :: subnam = 'convert_cyl2curv'
   integer      :: ii
   real(kind=8) :: x, y, z, r, th
   real(kind=8), parameter :: pi     = 4d0*atan(1d0)

   if (g%ivdir.eq.ikXDIR) then

      ! cylindrical grid with [x,r,th]-coordinates

      do ii = 1, g%ntot
         r  = g%coor(ii,2)
         th = g%coor(ii,3)
         y  = r * cos(th*pi/180d0)
         z  = r * sin(th*pi/180d0)
         g%coor(ii,2) = y
         g%coor(ii,3) = z
      enddo

   elseif (g%ivdir.eq.ikYDIR) then

      ! cylindrical grid with [r,y,th]-coordinates

      do ii = 1, g%ntot
         r  = g%coor(ii,1)
         th = g%coor(ii,3)
         x  = r * cos(th*pi/180d0)
         z  = r * sin(th*pi/180d0)
         g%coor(ii,1) = x
         g%coor(ii,3) = z
      enddo

   elseif (g%ivdir.eq.ikZDIR) then

      ! cylindrical grid with [r,th,z]-coordinates

      do ii = 1, g%ntot
         r  = g%coor(ii,1)
         th = g%coor(ii,2)
         x  = r * cos(th*pi/180d0)
         y  = r * sin(th*pi/180d0)
         g%coor(ii,1) = x
         g%coor(ii,2) = y
      enddo

   else
      call write_log(' Internal error: cyl2curv: incorrect ivdir')
   endif

   g%is_cylindr  = .false.
   g%is_curvilin = .true.
   g%lies_in_oyz = .false.
   g%ivdir       = 0
      
   g%x  => g%coor(:,1)
   g%y  => g%coor(:,2)
   g%z  => g%coor(:,3)
   g%r  => NULL()
   g%th => NULL()
   g%v  => NULL()

end subroutine convert_cyl2curv

!------------------------------------------------------------------------------------------------------------

end module m_grids
