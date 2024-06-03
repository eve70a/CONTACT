!------------------------------------------------------------------------------------------------------------
! m_interp_gf - interpolating routines for grid-based data
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_interp_gf
   use m_globals
   use m_markers
   use m_interp
   use m_grids
   use m_gridfunc
   implicit none
   private

   ! Debugging for module m_interp_gf

   public  interp_set_debug
   public  interpgf_set_debug

   integer  :: ldebug    =  0    ! local level of debugging
   integer  :: ii_debug  = -1    ! output point for which detailed info is requested (-1 = none)
   integer  :: iel_debug = -1    ! input element for which detailed info is requested (-1 = none)

   ! Interpolation routines defined on grids and grid functions:

   public  interp_aply_surf2unif_gf3
   public  interp_aply_unif2surf
   public  interp_surf2unif_gf3
   public  interp_curvgf2unifgf
   public  interp_splgf2unifgf

   public  interp_test

   public  interp_cartz2unif

   interface interp_cartz2unif
      module procedure interp_cartz2unif_grid
      module procedure interp_cartz2unif_gf3
   end interface interp_cartz2unif

!------------------------------------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------------------------------------

subroutine interp_set_debug(new_ldebug, new_ii_debug, new_iel_debug)
!--function: enable/disable debug output of all interpolation routines
   implicit none
!--subroutine arguments:
   integer, intent(in)           :: new_ldebug       ! level of debug output required
   integer, intent(in), optional :: new_ii_debug     ! specific point of interest for debugging
   integer, intent(in), optional :: new_iel_debug    ! specific point of interest for debugging

   call interp1d_set_debug(new_ldebug, new_ii_debug, new_iel_debug)
   call interpgf_set_debug(new_ldebug, new_ii_debug, new_iel_debug)

end subroutine interp_set_debug

!------------------------------------------------------------------------------------------------------------

subroutine interpgf_set_debug(new_ldebug, new_ii_debug, new_iel_debug)
!--function: enable/disable debug output of interpolation routines
   implicit none
!--subroutine arguments:
   integer, intent(in)           :: new_ldebug       ! level of debug output required
   integer, intent(in), optional :: new_ii_debug     ! specific point of interest for debugging
   integer, intent(in), optional :: new_iel_debug    ! specific point of interest for debugging

   ldebug = new_ldebug

   if (present(new_ii_debug)) then
      ii_debug = new_ii_debug
   endif
   if (present(new_iel_debug)) then
      iel_debug = new_iel_debug
   endif

   if (ldebug.ge.3) then
      write(bufout,'(a,i3,2(a,i7))') ' interpolations: debugging level =',ldebug,', ii_debug =',        &
                ii_debug,', iel_debug =', iel_debug
      call write_log(1, bufout)
   endif

end subroutine interpgf_set_debug

!------------------------------------------------------------------------------------------------------------

subroutine interp_aply_surf2unif_gf3(nnode, gf_node, nout, gf_out, ii2iel, ii2nod, wii2nod, ikarg,      &
                                     ierror, defval)
!--function: interpolate data given at the nodes of an input surface to the points of the output grid
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: nnode             ! number of nodes in the input surface
   real(kind=8), intent(in)  :: gf_node(nnode,3)  ! function values of surface nodes 
                                                  ! (nval=3: associated with (x,y,z) dirs)
   integer,      intent(in)  :: nout              ! number of points in the output grid
   integer,      intent(in)  :: ii2iel(nout)      ! surface element number for output grid points
   integer,      intent(in)  :: ii2nod(4,nout)    ! surface node numbers for output grid points
   real(kind=8), intent(in)  :: wii2nod(4,nout)   ! interpolation weights per surface node per output point
   type(t_gridfnc3)          :: gf_out            ! function values interpolated to the output grid
   integer,      intent(in)  :: ikarg             ! coordinate direction(s) of gf3 to be interpolated
   integer,      intent(out) :: ierror
   real(kind=8), optional    :: defval            ! default value
!--local variables:
!  character(len=*), parameter  :: subnam = 'interp_aply_surf2unif_gf3'
   integer          :: ii, j, ik0, ik1, ik
   logical          :: set_defval

   if (.not.gf_out%is_defined()) then
      call write_log('Internal error (interp_aply_surf2unif_gf3): gf_out not initialized properly.')
      call abort_run()
   endif

   ierror = 0
   set_defval = .false.
   if (present(defval)) set_defval = .true.

   ! expand range of coordinate directions

   call gf3_ikrange(ikarg, ik0, ik1)

   ! compute values of function gf_node at locations of the output grid

   do ik = ik0, ik1
      do ii = 1, nout

         ! if a surface element iel is found for output grid point ii:

         if (ii2iel(ii).ne.0) then

            ! sum up contributions of corners i_ll, i_lr, i_ur, i_ul (counter-clockwise)

            gf_out%val(ii,ik) = 0d0

            do j = 1, 4

               if (ldebug.ge.80 .and. ii.eq.26) then
                  write(bufout,*) 'ii_out=',ii,', cornr j=',j,': surf.nod=',ii2nod(j,ii)
                  call write_log(1, bufout)
               endif

               gf_out%val(ii,ik) =  gf_out%val(ii,ik) + wii2nod(j,ii) * gf_node(ii2nod(j,ii),ik)

               if (ldebug.ge.80 .and. ii.eq.26) then
                  write(bufout,*) 'comp. ik=',ik,': wgt=',wii2nod(j,ii),', gf_node=', gf_node(ii2nod(j,ii),ik)
                  call write_log(1, bufout)
               endif

            enddo

         else
            if (set_defval) gf_out%val(ii,ik) = defval
            if (ldebug.ge.20 .and. ii.ge.-10) then
               write(bufout, '(2(a,i3),a)') '  no value assigned to output grid point (',            &
                              gf_out%grid%ix(ii), ',',  gf_out%grid%iy(ii),')'
               call write_log(1, bufout)
            endif
         endif
      enddo ! ii
   enddo ! ik

end subroutine interp_aply_surf2unif_gf3

!------------------------------------------------------------------------------------------------------------

subroutine interp_aply_unif2surf(nin, gf_unif, nnode, f_node, ii2iel, ii2nod, wii2nod, ierror)
!--function: distribute grid-function defined on uniform input-grid to the nodal points of the output surface
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: nin              ! number of input grid points
   type(t_gridfnc3)          :: gf_unif          ! values defined on input grid
   integer,      intent(in)  :: nnode            ! number of nodes in the output surface
   type(t_vec),  intent(out) :: f_node(nnode)    ! value assigned to each node of the surface ((x,y,z) dirs)
   integer,      intent(in)  :: ii2iel(nin)      ! surface element number for input grid points
   integer,      intent(in)  :: ii2nod(4,nin)    ! surface node numbers for input grid points
   real(kind=8), intent(in)  :: wii2nod(4,nin)   ! interpolation weights per surface node per input grid point
   integer,      intent(out) :: ierror
!--local variables:
!  character(len=*), parameter  :: subnam = 'interp_aply_unif2surf'
   integer                      :: ii, inod, j

   ierror = 0

   ! initialize forces for all surface nodes

   if (nnode.le.0) then
      write(bufout,*) 'distribute: nnode=',nnode
      call write_log(1, bufout)
   endif
   if (size(f_node,1).le.0) then
      write(bufout,*) 'distribute: size f_node=',size(f_node,1)
      call write_log(1, bufout)
   endif

   do inod = 1, nnode
      f_node(inod)%v(1:3) = 0d0
   enddo

   ! loop over input grid points, distribute values to corresponding surface nodes

   do ii = 1, nin

      ! if a surface element iel is found for input grid point ii:

      if (ii2iel(ii).ne.0) then

         ! add contributions to corners i_ll, i_lr, i_ur, i_ul (counter-clockwise)

         do j = 1, 4
            if (ii2nod(j,ii).le.0 .or. ii2nod(j,ii).gt.nnode) then
               write(bufout,*) 'distribute: Internal error: for ii=',ii,', j=',j,' ii2nod=',ii2nod(j,ii)
               call write_log(1,bufout)
            else
               associate( fx => f_node(ii2nod(j,ii))%v(1), fy => f_node(ii2nod(j,ii))%v(2) )
               fx = fx + wii2nod(j,ii) * gf_unif%vx(ii)
               fy = fy + wii2nod(j,ii) * gf_unif%vy(ii)
!              fz = fz + wii2nod(j,ii) * gf_unif%vn(ii)
               end associate
            endif
         enddo

      endif
   enddo

end subroutine interp_aply_unif2surf

!------------------------------------------------------------------------------------------------------------

subroutine interp_cartz2unif_grid(g_surf, g_unif, my_ierror, defval, bicubic)
!--function: perform interpolation of z(x,y)-coordinates of a cartesian grid to a uniform (x,y)-grid
   implicit none
!--subroutine arguments:
   type(t_grid)              :: g_surf            ! input-grid with z(x,y)-function filled in
   type(t_grid)              :: g_unif            ! regular output-grid with (x,y) filled in
   integer,      intent(out) :: my_ierror
   real(kind=8), optional    :: defval            ! default value
   logical,      optional    :: bicubic           ! switch, default false
!--local variables:
!  character(len=*), parameter  :: subnam = 'interp_cartz2unif_grid'
   integer      :: nin, nout, sub_ierror
   real(kind=8) :: z_thrs, dx, dy
   integer,      dimension(:),   allocatable :: ii2iel  ! input element number for output grid points
   integer,      dimension(:,:), allocatable :: ii2nod  ! input node numbers for output grid points
   real(kind=8), dimension(:,:), allocatable :: wii2nod ! interpolation weights per node per grid point
   real(kind=8), dimension(:,:), allocatable :: fac_uv  ! relative (u,v) positions per grid point
   real(kind=8), dimension(:),   allocatable :: zdum
   logical                                   :: use_bicubic

   if (.not.g_unif%is_defined()) then
      call write_log('Internal error (interp_cartz2unif_grid): g_unif not initialized properly.')
      call abort_run()
   endif

   use_bicubic = .false.
   if (present(bicubic)) use_bicubic = bicubic

   my_ierror = 0
   if (min(g_surf%nx, g_surf%ny).eq.1) then
      write(bufout,'(2(a,i5),a)') ' Internal error: cartz2unif_grid: 1-d grids (', g_surf%nx,' x',      &
             g_surf%ny,') are not supported.'
      call write_log(1, bufout)
   endif
   if (ldebug.ge.2 .and. .not.g_unif%is_uniform) then
      call write_log(' interp_cartz2unif_grid: assuming the output grid is a cartesian product')
   endif
   if (ldebug.ge.1 .and. use_bicubic) then
      call write_log('interp_cartz2unif_grid: using bicubic interpolation')
   endif

   nin  = g_surf%ntot
   nout = g_unif%ntot

   if (g_unif%is_uniform) then
      dx = g_unif%dx
      dy = g_unif%dy
   else
      dx = -1d0
      dy = -1d0
   endif

   allocate(zdum(nin), ii2iel(nout), ii2nod(4,nout), wii2nod(4,nout), fac_uv(2,nout))

   ! calculate interpolation weights

   z_thrs      = 1d10
   zdum(1:nin) =  0d0
   call interp_wgt_surf2unif(g_surf%nx, g_surf%ny, nin, g_surf%x, g_surf%y, zdum, z_thrs,               &
                             g_unif%nx, g_unif%ny, nout, dx, dy, g_unif%x, g_unif%y,                    &
                             ii2iel, ii2nod, wii2nod, fac_uv, sub_ierror)
   if (my_ierror.eq.0) my_ierror = sub_ierror

   ! calculate interpolation itself

   if (present(defval) .and. use_bicubic) then
      call interp_aply_surf2unif_1d_bicubic(g_surf%nx, g_surf%ny, nin, g_surf%z, nout, ii2nod, fac_uv,  &
                            g_unif%z, defval)
   elseif (present(defval)) then
      call interp_aply_surf2unif_1d(nin, g_surf%z, nout, g_unif%z, ii2iel, ii2nod, wii2nod, sub_ierror, &
                                defval)
   elseif (use_bicubic) then
      call interp_aply_surf2unif_1d_bicubic(g_surf%nx, g_surf%ny, nin, g_surf%z, nout, ii2nod, fac_uv,  &
                            g_unif%z)
   else
      call interp_aply_surf2unif_1d(nin, g_surf%z, nout, g_unif%z, ii2iel, ii2nod, wii2nod, sub_ierror)
   endif
   if (my_ierror.eq.0) my_ierror = sub_ierror

   deallocate(zdum, ii2iel, ii2nod, wii2nod, fac_uv)
end subroutine interp_cartz2unif_grid

!------------------------------------------------------------------------------------------------------------

subroutine interp_cartz2unif_gf3(g_surf, gf_unif, my_ierror, defval)
!--function: perform interpolation of z(x,y)-coordinates of a cartesian grid to a gf3 defined on a
!            uniform (x,y)-grid
   implicit none
!--subroutine arguments:
   type(t_grid)                  :: g_surf            ! input-grid with z(x,y)-function filled in
   type(t_gridfnc3)              :: gf_unif           ! regular output-grid with (x,y) filled in
   integer,          intent(out) :: my_ierror
   real(kind=8),     optional    :: defval            ! default value
!--local variables:
!  character(len=*), parameter  :: subnam = 'interp_cartz2unif_gf3'
   integer      :: nin, nout, sub_ierror
   real(kind=8) :: z_thrs, dx, dy
   integer,      dimension(:),   allocatable :: ii2iel  ! input element number for output grid points
   integer,      dimension(:,:), allocatable :: ii2nod  ! input node numbers for output grid points
   real(kind=8), dimension(:,:), allocatable :: wii2nod ! interpolation weights per node per grid point
   real(kind=8), dimension(:,:), allocatable :: fac_uv  ! relative (u,v) positions per grid point
   real(kind=8), dimension(:),   allocatable :: zdum

   my_ierror = 0
   associate(g_unif => gf_unif%grid)

   if (.not.gf_unif%is_defined()) then
      call write_log('Internal error (interp_cartz2unif_gf3): gf_unif not initialized properly.')
      call abort_run()
   endif
   if (min(g_surf%nx, g_surf%ny).eq.1) then
      write(bufout,'(2(a,i5),a)') ' Internal error: cartz2unif_gf3: 1-d grids (', g_surf%nx,' x',       &
             g_surf%ny,') are not supported.'
      call write_log(1, bufout)
   endif
   if (ldebug.ge.2 .and. .not.g_unif%is_uniform) then
      call write_log(' interp_cartz2unif_gf3: assuming the output grid is a cartesian product')
   endif

   nin  = g_surf%ntot
   nout = g_unif%ntot

   if (g_unif%is_uniform) then
      dx = g_unif%dx
      dy = g_unif%dy
   else
      dx = -1d0
      dy = -1d0
   endif

   allocate(zdum(nin), ii2iel(nout), ii2nod(4,nout), wii2nod(4,nout), fac_uv(2,nout))

   ! calculate interpolation weights

   z_thrs      = 1d10
   zdum(1:nin) =  0d0
   call interp_wgt_surf2unif(g_surf%nx, g_surf%ny, nin, g_surf%x, g_surf%y, zdum, z_thrs,               &
                             g_unif%nx, g_unif%ny, nout, dx, dy, g_unif%x, g_unif%y,                    &
                             ii2iel, ii2nod, wii2nod, fac_uv, sub_ierror)
   if (my_ierror.eq.0) my_ierror = sub_ierror

   ! calculate interpolation itself

   if (present(defval)) then
      call interp_aply_surf2unif_1d(nin, g_surf%z, nout, gf_unif%vn, ii2iel, ii2nod, wii2nod,           &
                             sub_ierror, defval)
   else
      call interp_aply_surf2unif_1d(nin, g_surf%z, nout, gf_unif%vn, ii2iel, ii2nod, wii2nod, sub_ierror)
   endif
   if (my_ierror.eq.0) my_ierror = sub_ierror

   deallocate(zdum, ii2iel, ii2nod, wii2nod, fac_uv)
   end associate
end subroutine interp_cartz2unif_gf3

!------------------------------------------------------------------------------------------------------------

subroutine interp_curvgf2unifgf(gf_curv, gf_unif, my_ierror, defval)
!--function: perform interpolation of a grid function on a curvilinear grid to a grid function defined on a
!            uniform (x,y)-grid
   implicit none
!--subroutine arguments:
   type(t_gridfnc3)              :: gf_curv           ! gf3 on curvilinear input-grid
   type(t_gridfnc3)              :: gf_unif           ! gf3 on uniform output-grid
   integer,          intent(out) :: my_ierror
   real(kind=8),     optional    :: defval            ! default value
!--local variables:
!  character(len=*), parameter  :: subnam = 'interp_cartz2unif_gf3'
   integer      :: nin, nout, sub_ierror
   real(kind=8) :: z_thrs, dx, dy
   integer,      dimension(:),   allocatable :: ii2iel  ! input element number for output grid points
   integer,      dimension(:,:), allocatable :: ii2nod  ! input node numbers for output grid points
   real(kind=8), dimension(:,:), allocatable :: wii2nod ! interpolation weights per node per grid point
   real(kind=8), dimension(:,:), allocatable :: fac_uv  ! relative (u,v) positions per grid point
   real(kind=8), dimension(:),   allocatable :: zdum

   my_ierror = 0
   associate(g_curv => gf_curv%grid, g_unif => gf_unif%grid)

   if (.not.gf_unif%is_defined()) then
      call write_log('Internal error (interp_curvgf2unifgf): gf_unif not initialized properly.')
      call abort_run()
   endif
   if (min(g_curv%nx, g_curv%ny).eq.1) then
      write(bufout,'(2(a,i5),a)') ' Internal error: curvgf2unifgf: 1-d grids (', g_curv%nx,' x',        &
             g_curv%ny,') are not supported.'
      call write_log(1, bufout)
   endif
   if (ldebug.ge.2 .and. .not.g_unif%is_uniform) then
      call write_log(' interp_curvgf2unifgf: assuming the output grid is a cartesian product')
   endif

   nin  = g_curv%ntot
   nout = g_unif%ntot

   if (g_unif%is_uniform) then
      dx = g_unif%dx
      dy = g_unif%dy
   else
      dx = -1d0
      dy = -1d0
   endif

   allocate(zdum(nin), ii2iel(nout), ii2nod(4,nout), wii2nod(4,nout), fac_uv(2,nout))

   ! calculate interpolation weights

   z_thrs      = 1d10
   zdum(1:nin) =  0d0
   call interp_wgt_surf2unif(g_curv%nx, g_curv%ny, nin, g_curv%x, g_curv%y, zdum, z_thrs,               &
                             g_unif%nx, g_unif%ny, nout, dx, dy, g_unif%x, g_unif%y,                    &
                             ii2iel, ii2nod, wii2nod, fac_uv, sub_ierror)
   if (my_ierror.eq.0) my_ierror = sub_ierror

   ! calculate interpolation itself

   if (present(defval)) then
      call interp_aply_surf2unif_gf3(nin, gf_curv%val, nout, gf_unif, ii2iel, ii2nod, wii2nod, ikALL,   &
                             sub_ierror, defval)
   else
      call interp_aply_surf2unif_gf3(nin, gf_curv%val, nout, gf_unif, ii2iel, ii2nod, wii2nod, ikALL,   &
                             sub_ierror)
   endif
   if (my_ierror.eq.0) my_ierror = sub_ierror

   deallocate(zdum, ii2iel, ii2nod, wii2nod, fac_uv)
   end associate
end subroutine interp_curvgf2unifgf

!------------------------------------------------------------------------------------------------------------

subroutine interp_splgf2unifgf(gf_spl, gf_unif, sg_unif, ikarg, my_ierror, defval)
!--function: perform interpolation of a grid function on a grid with spline (arc-length) parametrization
!            to a grid function defined with uniform (x, s_l) parametrization
   implicit none
!--subroutine arguments:
   type(t_gridfnc3)              :: gf_spl            ! gf3 on input-grid with spline representation
   type(t_gridfnc3)              :: gf_unif           ! gf3 on uniform output-grid
   real(kind=8)                  :: sg_unif           ! s_g-value for the output grid ref. marker
   integer,          intent(in)  :: ikarg             ! coordinate direction(s) of gf3 to be interpolated
   integer,          intent(out) :: my_ierror
   real(kind=8),     optional    :: defval            ! default value
!--local variables:
   integer      :: nin, nout, ik0, ik1, ik, sub_ierror
   real(kind=8) :: z_thrs, dx, dy
   integer,      dimension(:),   allocatable :: ii2iel  ! input element number for output grid points
   integer,      dimension(:,:), allocatable :: ii2nod  ! input node numbers for output grid points
   real(kind=8), dimension(:,:), allocatable :: wii2nod ! interpolation weights per node per grid point
   real(kind=8), dimension(:,:), allocatable :: fac_uv  ! relative (u,v) positions per grid point
   real(kind=8), dimension(:),   allocatable :: zdum, s_unif

   my_ierror = 0
   associate(g_spl  => gf_spl%grid, g_unif => gf_unif%grid)

   if (.not.gf_unif%is_defined()) then
      call write_log('Internal error (interp_splgf2unifgf): gf_unif not initialized properly.')
      call abort_run()
   endif

   nin  = g_spl%ntot
   nout = g_unif%ntot

   ! expand range of coordinate directions

   call gf3_ikrange(ikarg, ik0, ik1)

   ! map output y/s_l-coordinates to input s_g arc-length coordinates

   allocate(s_unif(nout))
   s_unif(1:nout) = sg_unif + g_unif%y(1:nout)

   if (min(g_spl%nx, g_spl%ny).eq.1) then

      if (ldebug.ge.2) call write_log(' interp_splgf2unifgf: using 1-d interpolation')

      do ik = ik0, ik1
         call interp_1d_to_2d( nin, g_spl%s_prf, gf_spl%val(:,ik), g_unif%nx, g_unif%ny, s_unif,    &
                gf_unif%val(:,ik), sub_ierror, defval)
         if (my_ierror.eq.0) my_ierror = sub_ierror
      enddo

   else

      if (ldebug.ge.2) call write_log(' interp_splgf2unifgf: using 1-d interpolation')

      if (ldebug.ge.2 .and. .not.g_unif%is_uniform) then
         call write_log(' interp_splgf2unifgf: assuming the output grid is a cartesian product')
      endif

      if (g_unif%is_uniform) then
         dx = g_unif%dx
         dy = g_unif%dy
      else
         dx = -1d0
         dy = -1d0
      endif

      allocate(zdum(nin), ii2iel(nout), ii2nod(4,nout), wii2nod(4,nout), fac_uv(2,nout))

      ! calculate interpolation weights using g_spl%s_prf <--> s_unif

      z_thrs      = 1d10
      zdum(1:nin) =  0d0
      call interp_wgt_surf2unif(g_spl%nx,  g_spl%ny,  nin,  g_spl%x, g_spl%s_prf, zdum, z_thrs,         &
                                g_unif%nx, g_unif%ny, nout, dx, dy, g_unif%x, s_unif,                   &
                                ii2iel, ii2nod, wii2nod, fac_uv, sub_ierror)
      if (my_ierror.eq.0) my_ierror = sub_ierror

      ! calculate interpolation itself

      if (present(defval)) then
         call interp_aply_surf2unif_gf3(nin, gf_spl%val, nout, gf_unif, ii2iel, ii2nod, wii2nod, ikarg, & 
                        sub_ierror, defval)
      else
         call interp_aply_surf2unif_gf3(nin, gf_spl%val, nout, gf_unif, ii2iel, ii2nod, wii2nod, ikarg, &
                        sub_ierror)
      endif
      if (my_ierror.eq.0) my_ierror = sub_ierror

      deallocate(zdum, ii2iel, ii2nod, wii2nod, fac_uv)
   endif
   deallocate(s_unif)
   end associate

end subroutine interp_splgf2unifgf

!------------------------------------------------------------------------------------------------------------

subroutine interp_surf2unif_gf3(nnode_x, nnode_y, nnode, pos_node, gf_node, z_thrs, mx, my, nout,       &
                                dx, dy, pos_out, gf_out, my_ierror, defval)
!--function: perform interpolation from a curvi-linear input grid (surface) to a uniform output-grid.
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: nnode_x, nnode_y, nnode ! number of points in the input grid
   real(kind=8), intent(in)  :: pos_node(nnode,3) ! (x,y,z) coordinates of input grid points
   real(kind=8), intent(in)  :: gf_node(nnode,3)  ! function values of surface nodes 
                                                  ! (nval=3: associated with (x,y,z) dirs)
   real(kind=8), intent(in)  :: z_thrs            ! typical dimension in z-direction, threshold for bounding box
   integer,      intent(in)  :: mx, my, nout      ! number of points in the output grid
   real(kind=8), intent(in)  :: dx, dy            ! step-sizes used in the output grid
   real(kind=8), intent(in)  :: pos_out(nout,3)   ! coordinates of the output grid points
   type(t_gridfnc3)          :: gf_out            ! function values interpolated to the output grid
   integer,      intent(out) :: my_ierror
   real(kind=8), optional    :: defval            ! default value
!--local variables:
!  character(len=*), parameter  :: subnam = 'interp_surf2unif_gf3'
   integer                   :: sub_ierror
   integer,      dimension(:),   allocatable :: ii2iel  ! input element number for output grid points
   integer,      dimension(:,:), allocatable :: ii2nod  ! input node numbers for output grid points
   real(kind=8), dimension(:,:), allocatable :: wii2nod ! interpolation weights per node per grid point
   real(kind=8), dimension(:,:), allocatable :: fac_uv  ! relative (u,v) positions per grid point

   if (.not.gf_out%is_defined()) then
      call write_log('Internal error (interp_surf2unif_gf3): gf_out not initialized properly.')
      call abort_run()
   endif

   my_ierror = 0
   allocate(ii2iel(nout), ii2nod(4,nout), wii2nod(4,nout), fac_uv(2,nout))

   ! calculate interpolation weights

   call interp_wgt_surf2unif(nnode_x, nnode_y, nnode, pos_node(:,1), pos_node(:,2), pos_node(:,3),      &
                             z_thrs, mx, my, nout, dx, dy, pos_out(:,1), pos_out(:,2), ii2iel, ii2nod,  &
                             wii2nod, fac_uv, sub_ierror)
   if (my_ierror.eq.0) my_ierror = sub_ierror

   ! calculate interpolation itself

   if (present(defval)) then
      call interp_aply_surf2unif_gf3(nnode, gf_node, nout, gf_out, ii2iel, ii2nod, wii2nod, ikALL,      &
                        sub_ierror, defval)
   else
      call interp_aply_surf2unif_gf3(nnode, gf_node, nout, gf_out, ii2iel, ii2nod, wii2nod, ikALL,      &
                        sub_ierror)
   endif
   if (my_ierror.eq.0) my_ierror = sub_ierror

   deallocate(ii2iel, ii2nod, wii2nod, fac_uv)
end subroutine interp_surf2unif_gf3

!------------------------------------------------------------------------------------------------------------

subroutine interp_test
!--function: test interpolation routines
   implicit none
!--local variables:
   integer       :: nx_in, ny_in, nx_out, ny_out, nin, nout1, nout2, iin, iout, my_ierror, sub_ierror
   real(kind=8)  :: exterval, x1_out, y1_out, dx_out, dy_out
   real(kind=8), dimension(:),   allocatable :: xin, fin, x_in, y_in, z_in
   real(kind=8), dimension(:,:), allocatable :: xout, fout
   type(t_grid)  :: g_curv, g_unif

   call write_log(' ----------------------------- Interp_test: -------------------------------')
   my_ierror = 0

   ! linear test function, x-values increasing

   if (.false.) then

      nin = 5
      allocate(xin(nin), fin(nin))

      do iin = 1, nin
         xin(iin) = real(iin)
         fin(iin) = 2d0 * xin(iin)
      enddo

      write(bufout,'(a,i3,a)') ' input: linear function with', nin,' points:'
      call write_log(1, bufout)
      write(bufout,'(a,10f8.3)') '    xi = ',(xin(iin), iin=1,nin)
      call write_log(1, bufout)
      write(bufout,'(a,10f8.3)') ' f(xi) = ',(fin(iin), iin=1,nin)
      call write_log(1, bufout)

      ! interpolation to 2d array with exterval

      nout1 = 7
      nout2 = 2
      allocate(xout(nout1,nout2), fout(nout1,nout2))

      do iout = 1, nout1
         xout(iout,1) = real(iout) - 0.3d0
         xout(iout,2) = 1.25d0 * real(nout1-iout) - 2d0
      enddo
      fout(1:nout1,1:nout2) = -999d0

      exterval = 99d0

      call interp_1d(nin, xin, fin, nout1, nout2, xout, fout, sub_ierror, exterval)
      if (my_ierror.eq.0) my_ierror = sub_ierror

      call write_log(' output of interpolation with default value:')
      write(bufout,'(a,10f8.3)') '    xj = ',(xout(iout,1), iout=1,nout1)
      call write_log(1, bufout)
      write(bufout,'(a,10f8.3)') ' f(xj) = ',(fout(iout,1), iout=1,nout1)
      call write_log(1, bufout)
      write(bufout,'(a,10f8.3)') '    xj = ',(xout(iout,2), iout=1,nout1)
      call write_log(1, bufout)
      write(bufout,'(a,10f8.3)') ' f(xj) = ',(fout(iout,2), iout=1,nout1)
      call write_log(1, bufout)

      ! decreasing input-function

      do iin = 1, nin
         xin(iin) = real(nin - iin)
         fin(iin) = 2d0 * xin(iin)
      enddo
      xin(3) = real(5-2)

      write(bufout,'(/,a,i3,a)') ' input: function with jump,', nin,' points with x decreasing:'
      call write_log(2, bufout)
      write(bufout,'(a,10f8.3)') '    xi = ',(xin(iin), iin=1,nin)
      call write_log(1, bufout)
      write(bufout,'(a,10f8.3)') ' f(xi) = ',(fin(iin), iin=1,nin)
      call write_log(1, bufout)

      ! interpolation to 2d array without exterval

      do iout = 1, nout1
         xout(iout,1) = 0.9d0 * real(iout) - 2d0
         xout(iout,2) = 1.1d0 * real(nout1-iout) - 2d0
      enddo
      fout(1:nout1,1:nout2) = -999d0

      call interp_1d(nin, xin, fin, nout1, nout2, xout, fout, sub_ierror)
      if (my_ierror.eq.0) my_ierror = sub_ierror

      call write_log(' output of interpolation with no default (constant extrapolation):')
      write(bufout,'(a,10f8.3)') '    xj = ',(xout(iout,1), iout=1,nout1)
      call write_log(1, bufout)
      write(bufout,'(a,10f8.3)') ' f(xj) = ',(fout(iout,1), iout=1,nout1)
      call write_log(1, bufout)
      write(bufout,'(a,10f8.3)') '    xj = ',(xout(iout,2), iout=1,nout1)
      call write_log(1, bufout)
      write(bufout,'(a,10f8.3)') ' f(xj) = ',(fout(iout,2), iout=1,nout1)
      call write_log(1, bufout)
   endif

   ! test interpolation of z(x,y)-coordinates from curvilinear cartesian grid to a uniform (x,y)-grid

   nx_in = 3
   ny_in = 2
   nin = nx_in * ny_in
   allocate(x_in(nin), y_in(nin), z_in(nin))

   x_in = (/ 0d0, 1d0, 2d0,             &
             0d0, 1d0, 2d0 /)
   y_in = (/ 0d0, 0d0, 0d0,             &
             1d0, 1d0, 1d0 /)
   z_in = (/ 0d0, 2d0, 4d0,             &
             6d0, 8d0, 10d0 /)
   call grid_create_curvil(g_curv, nx_in, ny_in, x_in, y_in, z_in)
   call grid_print(g_curv, 'input: curv.grid', 5)

   nx_out = 6
   ny_out = 1
   x1_out = -0.333d0
   y1_out = 0.5d0
   dx_out = 0.334d0
   dy_out = 1d0
   call grid_create_uniform(g_unif, nxarg=nx_out, x0arg=x1_out, dxarg=dx_out,                           &
                                    nyarg=ny_out, y0arg=y1_out, dyarg=dy_out)
   call grid_print(g_unif, 'input: unif.grid', 5)

   exterval = -777d0
   call interp_set_debug(10, 1, 1)
   call interp_cartz2unif_grid(g_curv, g_unif, sub_ierror, exterval)
   if (my_ierror.eq.0) my_ierror = sub_ierror
   call interp_set_debug(0)
   call grid_print(g_unif, 'output: unif.grid', 5)

   call grid_destroy(g_curv)
   call grid_destroy(g_unif)

   call write_log(' ')
   call write_log(' ---------------------------- end Interp_test ------------------------------')

end subroutine interp_test

!------------------------------------------------------------------------------------------------------------

end module m_interp_gf
