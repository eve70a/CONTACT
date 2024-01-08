!------------------------------------------------------------------------------------------------------------
! m_wr_totforce - iterative solution of w/r contact (module 1) with total forces prescribed
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_wr_totforce

use m_wrprof_data
use m_hertz
use m_visc
use m_soutpt
use m_wr_solvecp
use m_wr_output
use m_wr_brentmeth

implicit none
private

public  wr_contact

private wr_contact_init_hertz
private wr_contact_fz_brent

private wr_contact_fx_secant
private secant_fx_print

private wr_contact_fy_brent

contains

!------------------------------------------------------------------------------------------------------------

   subroutine wr_contact(wtd, my_ierror)
!--purpose: process the whole w/r problem for a case, including NR-loops for prescribed forces
      implicit none
!--subroutine arguments:
      type(t_ws_track)              :: wtd
      integer,          intent(out) :: my_ierror
!--local variables:
      logical            :: lfound
      integer            :: iestim_br, x_locate, icp, icpo, i, sub_ierror, known_err(10)

      my_ierror = 0
      wtd%meta%itforce    = 0
      wtd%meta%itforc_out = 0
      wtd%meta%itforc_inn = 0

      ! idebug: 0 = errors/silent, 1 = warnings/info, >=2 = flow/debug

      x_locate = wtd%ic%x_locate

      ! set rail deflections if not used in this case

      if (wtd%ic%force1.ne.3) then
         wtd%trk%dy_defl = 0d0
         wtd%trk%dz_defl = 0d0
      endif

      ! adjust contact grid when switching from T=0 or 3 to T=1 or 2

      if (wtd%meta%ncase.gt.1) call shift_grid_super(wtd)

      ! switch between solvers dependent on prescribed total forces

      if (wtd%ic%norm.le.0) then

         ! N=0, F=0: solve w/r problem with positions & velocities prescribed

         if (x_locate.ge.3) call write_log(' wr_contact: solving for given position...')
         call wr_contact_pos(wtd, x_locate, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

      elseif (wtd%ic%force1.le.0) then

         ! N=1, F=0: solve w/r problem with prescribed vertical force

         if (x_locate.ge.3) call write_log(' wr_contact: solving for given vertical force...')

         iestim_br = 0
         if (wtd%ic%iestim.ge.1) iestim_br = 1

         call wr_contact_fz_brent(wtd, iestim_br, 0d0, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

      elseif (wtd%ic%force1.eq.1 .or. wtd%ic%force1.eq.2) then

         ! N=1, F=1,2: solve w/r problem with prescribed vertical and longitudinal forces

         if (x_locate.ge.3) call write_log(' wr_contact: solving for given vert./long. forces (secant)...')
         call wr_contact_fx_secant(wtd, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

      elseif (wtd%ic%force1.eq.3) then

         ! N=1, F=3: solve w/r problem with prescribed vertical and lateral forces - Brent method

         if (x_locate.ge.3) call write_log(' wr_contact: solving for given vert./lat. forces (brent)...')
         call wr_contact_fy_brent(wtd, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

      else

         call write_log(' Internal Error: invalid N/F-digits')
         call abort_run()

      endif

      ! translate errors

      if (my_ierror.eq.-1 .or. my_ierror.eq.-3) then ! -1 = no solution: NaN in residual
         my_ierror = CNTC_err_ftot                   ! -3 = zero determinant
      elseif (my_ierror.eq.-2) then                  ! -2 = residual > tolerance
         my_ierror = CNTC_err_tol
      endif

      if (my_ierror.ne.0) then
         known_err(1:9) = (/ CNTC_err_allow, CNTC_err_search, CNTC_err_ftot, CNTC_err_tol,              &
                             CNTC_err_norm,  CNTC_err_tang, CNTC_err_icp, CNTC_err_profil,              &
                             CNTC_err_frclaw /)
         lfound = .false.
         do i = 1, 9
            if (my_ierror.eq.known_err(i)) lfound = .true.
         enddo
         if (.not.lfound) my_ierror = CNTC_err_other
      endif

      ! if requested: write mat-files for final solution

      if (wtd%ic%matfil_surf.ge.1) then
         if (x_locate.ge.5) call write_log(' ...wr_contact: calling writmt')
         do icp = 1, wtd%numcps
            associate( gd => wtd%allcps(icp)%cp%gd )
            gd%ic%matfil_surf = wtd%ic%matfil_surf
            call writmt (gd%meta, gd%ic, gd%cgrid_cur, gd%potcon_cur, gd%geom%hs1, gd%mater, gd%fric,   &
                gd%kin, gd%outpt1, wtd%ic%is_left_side())
            end associate
         enddo
      endif

      ! compute subsurface stresses when requested

      if (wtd%ic%stress.ge.1) then
         call wr_subsurf(wtd, x_locate)
      endif

      ! write output when R=0 or 1

      if (wtd%ic%return.le.1) call wr_output(wtd)

      ! destroy remaining contact patches (gds) of the previous time

      do icpo = wtd%numcps+1, wtd%numtot
         if (x_locate.ge.1) then
            write(bufout,'(2(a,i3),a)') ' There are',wtd%numcps,' contact patches and',         &
                wtd%numtot-wtd%numcps,' unconnected patches for previous time'
            call write_log(1, bufout)
         endif
         if (associated(wtd%allcps(icpo)%cp)) then
            if (x_locate.ge.1 .and. associated(wtd%allcps(icpo)%cp%gd)) then
               write(bufout,'(a,i3)') ' wr_contact: destroy gd for icpo=',icpo
               call write_log(1, bufout)
            elseif (x_locate.ge.3) then
               write(bufout,'(a,i3)') ' wr_contact: destroy allcps cp=',icpo
               call write_log(1, bufout)
            endif
            call cp_destroy(wtd%allcps(icpo))
         endif
      enddo
      wtd%numtot = wtd%numcps

      ! set ic_prv for next case

      wtd%ic_prv = wtd%ic

   end subroutine wr_contact

!------------------------------------------------------------------------------------------------------------

   subroutine shift_grid_super(wtd)
!--purpose: add supergrid offset to contact grid when switching from T=0 or 3 to T=1 or 2
      implicit none
!--subroutine arguments:
      type(t_ws_track)                 :: wtd
!--local variables:
      logical           :: use_super_prv, use_super_new
      integer           :: icp
      real(kind=8)      :: s_ws, cref_x, sr_ref, xref_pot, yref_pot

      use_super_prv = wtd%ic_prv%tang.eq.1 .or. wtd%ic_prv%tang.eq.2
      use_super_new = wtd%ic%tang.eq.1 .or. wtd%ic%tang.eq.2

      if (.not.use_super_prv .and. use_super_new) then

         ! no supergrid used in previous case / using supergrid in new case

         if (wtd%ic%x_force.ge.1) call write_log(' shift_grid_super: adding supergrid offset')

         do icp = 1, wtd%numcps
            associate( gd => wtd%allcps(icp)%cp%gd )

            ! get s_ws, cref_x, sr_ref from previous case

            s_ws   = gd%meta%s_ws
            cref_x = gd%meta%xcp_tr
            sr_ref = gd%meta%scp_r

            ! set position of contact reference in terms of super-grid coordinates

            if (wtd%ic%tang.eq.1) then     ! transient shift: world/material-fixed super-grid
               xref_pot  = s_ws + cref_x
               yref_pot  = sr_ref
            elseif (wtd%ic%tang.eq.2) then ! transient rolling: super-grid fixed in lateral direction
               xref_pot  = 0d0
               yref_pot  = sr_ref
            else                           ! steady rolling: contact grid centered at contact reference marker
               xref_pot  = 0d0
               yref_pot  = 0d0
            endif

            ! snap contact reference to nearest supergrid location

            if (wtd%ic%x_force.ge.1) then
               write(bufout,'(2(a,f12.6),a)') ' yref =', yref_pot,' =', yref_pot/gd%cgrid_cur%dy,' * dy'
               call write_log(1, bufout)
            endif

            xref_pot = nint(xref_pot / gd%cgrid_cur%dx) * gd%cgrid_cur%dx
            yref_pot = nint(yref_pot / gd%cgrid_cur%dy) * gd%cgrid_cur%dy

            ! shift grid from contact reference to supergrid coordinate reference

            gd%potcon_cur%xl = gd%potcon_cur%xl + xref_pot
            gd%potcon_cur%yl = gd%potcon_cur%yl + yref_pot
            call potcon_fill(  gd%potcon_cur )
            call potcon_cgrid( gd%potcon_cur, gd%cgrid_cur)

            end associate
         enddo ! icp

      endif ! add offset

   end subroutine shift_grid_super

!------------------------------------------------------------------------------------------------------------

   subroutine wr_contact_init_hertz(wtd, its, dfz_dzws, aa, bb, imeth, x_force, x_locate, my_ierror)
!--purpose: determine an initial bracket for z_ws for starting the Brent algorithm
      implicit none
!--subroutine arguments:
      type(t_ws_track)                 :: wtd
      integer,           intent(in)    :: imeth
      integer,           intent(in)    :: x_force, x_locate
      integer,           intent(out)   :: my_ierror
      real(kind=8),      intent(out)   :: dfz_dzws, aa, bb
      type(t_brent_its), intent(inout) :: its
!--local variables:
      integer                  :: k, ic_norm, ipotcn, ic_return, sub_ierror
      real(kind=8)             :: cp, rho, pen, e_star, epshz, azz, x_0, x_new, r_new

      if (x_force.ge.3) call write_log(' --- Start subroutine wr_contact_init_hertz ---')

      my_ierror = 0

      ! initial estimate for z-position

      k   = 0
      wtd%meta%itforc_inn = k
      wtd%ws%z = 0d0

      ! analyze geometrical problem for the initial z-position, without actual solving

      ic_return = wtd%ic%return
      wtd%ic%return = 3

      call wr_contact_pos(wtd, x_locate, sub_ierror)
      if (my_ierror.eq.0) my_ierror = sub_ierror

      wtd%ic%return = ic_return

      if (max(x_force,x_locate).ge.3) then
         write(bufout,'(3(a,g14.6))') ' z_ws =',wtd%ws%z
         call write_log(1, bufout)
         write(bufout,'(3(a,g14.6))') ' Overall minimum gap=', wtd%ws%gap_min,', a1=',wtd%ws%a1,     &
                ', b1=',wtd%ws%b1
         call write_log(1, bufout)
      endif

      ! Abort in case of an error

      if (my_ierror.ne.0) then
         wtd%numcps = 0
         if (max(x_force,x_locate).ge.1) then
            write(bufout,'(a,i6,a)') ' An error occurred (',my_ierror,'), skipping computation.'
            call write_log(1, bufout)
         endif
         return
      endif

      ! Abort if the profiles have no overlap at all

      if (.not.wtd%ws%has_overlap) then
         wtd%numcps = 0
         if (max(x_force,x_locate).ge.1) then
            call write_log(' The profiles have no overlap at all, skipping computation.')
         endif
         return
      endif
         
      ! store initial point { x_0, r_0 }, print output

      x_0 = wtd%ws%z + wtd%ws%gap_min

      if (imeth.eq.IMETH_BRENT) then

         call brent_its_add_iterate(its, k, x_0, -wtd%ws%fz_inp, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, x_force)

      endif

      ! Use the Hertzian approach to get a sensible initial z-position

      ! solve the Hertzian problem for estimated curvatures, approximating fn ~=~ fz

      ic_norm =  1
      ipotcn  = -1
      e_star  = wtd%mater%ga / (1d0 - wtd%mater%nu)
      epshz   = wtd%solv%eps
         
      call hzcalc3d (e_star, epshz, ipotcn, wtd%ws%a1, wtd%ws%b1, aa, bb, ic_norm, pen,                 &
                     wtd%ws%fz_inp, cp, rho)

      ! estimate derivative dFz/dzws ~=~ dFn/dpen

      if (pen.le.1d-6) then
         ! no contact: using influence coefficient Azz(0,0) for 1 element in contact

         call azz_one_element(wtd%ic, wtd%mater, wtd%discr%dx, wtd%discr%ds, azz)
         dfz_dzws = wtd%discr%dx * wtd%discr%ds * wtd%mater%ga / azz 
         write(bufout,*) 'gap_min=',wtd%ws%gap_min,', pen=',pen,', using dfz_dzws=',dfz_dzws
         call write_log(1, bufout)
      else
         dfz_dzws = (1.5d0 * wtd%ws%fz_inp / pen) * cos(wtd%ws%delt_min)

         if (x_force.ge.3) then
            write(bufout,'(a,f12.6,a,f8.4,a,f12.1)') ' pen=',pen,', delt=',wtd%ws%delt_min,             &
                                ': est. dfz_dzws=',dfz_dzws
            call write_log(1, bufout)
            write(bufout,'(a,f12.3,a,f12.6,a,f8.4,a,f12.1)') '    1.5 * ',wtd%ws%fz_inp,' /',pen,' *',  &
                                cos(wtd%ws%delt_min),' =', dfz_dzws
            call write_log(1, bufout)
         endif
      endif

      ! lower the wheel-set by gap_min and pen to get near the desired approach

      ! write(bufout,*) ' x_1=', x_0, ' +', pen, ' /', cos(wtd%ws%delt_min)
      ! call write_log(1, bufout)

      k     = 1
      wtd%meta%itforc_inn = k
      x_new = x_0 + pen / max(0.5d0, cos(wtd%ws%delt_min))

      wtd%ws%z = x_new

      ! Brent: solve contact problem for the current estimate

      if (imeth.eq.IMETH_BRENT) then

         call wr_contact_pos(wtd, x_locate, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

         r_new    = wtd%ftrk%z() - wtd%ws%fz_inp
         dfz_dzws = 1.5d0 * wtd%ftrk%z() / (x_new - x_0)

         ! write(bufout,*) 'dfz_dzws= =', dfz_dzws
         ! call write_log(1, bufout)

         call brent_its_add_iterate(its, k, x_new, r_new, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, x_force)

      endif
         
   end subroutine wr_contact_init_hertz

!------------------------------------------------------------------------------------------------------------

   subroutine wr_contact_fz_brent(wtd, iestim_br, z_offset, my_ierror)
!--purpose: process the w/r problem for a case with prescribed vertical force with Brents algorithm
      implicit none
!--subroutine arguments:
      type(t_ws_track)              :: wtd
      integer,          intent(in)  :: iestim_br
      real(kind=8),     intent(in)  :: z_offset
      integer,          intent(out) :: my_ierror
!--local variables:
      integer, parameter       :: maxit = 100
      logical                  :: used_bisec, ldone
      integer                  :: k, x_force, x_locate, ic_return, sub_ierror
      real(kind=8)             :: ftarg, res_fk, tol_fk, dfk, tol_xk, x_0, x_new, r_new, dfz_dzws, aa, bb
      type(t_brent_its)        :: its

      my_ierror = 0
      x_force  = wtd%ic%x_force
      x_locate = wtd%ic%x_locate
      ! if (wtd%meta%itforc_out.eq.2) x_force = 3

      if (max(x_force,x_locate).ge.2) call write_log(' --- Start subroutine wr_contact_fz_brent ---')

      ! Abort if a zero or negative total force is prescribed

      if (wtd%ws%fz_inp.lt.1d-9) then
         wtd%numcps = 0
         if (x_locate.ge.1) then  
            call write_log(' Prescribed total force is zero or negative, skipping computation.')
         endif
         return
      endif

      ! prepare structure to hold iterates

      ftarg = wtd%ws%fz_inp
      call brent_its_init(its, ikZDIR, maxit, ftarg)

      ! initialize iteration; set zero point and initial estimate

      !   0 = full initialization - using geometrical analysis + Hertz approximation
      !   1 = continuation, using estimate z_ws from wtd%ws
      !   2 = continuation, using geometrical analysis + offset z_offset provided

      if (iestim_br.eq.0) then
         wtd%ws%z  = 0d0
      endif

      if (iestim_br.eq.2) then

         ! analyze geometrical problem for the initial z-position, without actual solving

         wtd%ws%z  = 0d0
         ic_return = wtd%ic%return
         wtd%ic%return = 3

         call wr_contact_pos(wtd, x_locate, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         wtd%ic%return = ic_return

         wtd%ws%z  = wtd%ws%z_cnt0 + z_offset
      endif

      if (iestim_br.ge.1) then

         ! methods 1, 2: use initial estimate from ws%z, wtd%dfz_dzws

         if (x_force.ge.3) then
            write(bufout,'(a,i3)')   ' wr_contact_fz_brent: starting iteration k=',1
            call write_log(1, bufout)
         endif

         wtd%meta%itforc_inn = 1
         call wr_contact_pos(wtd, x_locate, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

         ! store initial point { x_0, r_0 } at which initial contact occurs

         k   = 0
         x_0 = wtd%ws%z_cnt0
         call brent_its_add_iterate(its, k, x_0, -wtd%ws%fz_inp, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, x_force)

         ! store first point { x_1, r_1 } for initial estimate ws%z

         k   = 1
         wtd%meta%itforc_inn = k
         x_new    = wtd%ws%z
         r_new    = wtd%ftrk%z() - wtd%ws%fz_inp
         dfz_dzws = wtd%dfz_dzws

         call brent_its_add_iterate(its, k, x_new, r_new, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, x_force)

      else

         ! analyze geometry for overall minimum gap and curvatures, set initial estimates k=0, k=1

         call wr_contact_init_hertz(wtd, its, dfz_dzws, aa, bb, IMETH_BRENT, x_force, x_locate,     &
                        sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

      endif

      ! Abort in case of error or if the profiles have no overlap at all

      if (my_ierror.ne.0 .or. .not.wtd%ws%has_overlap) then
         call write_log(' no overlap...')
         call total_forces_moments(wtd, x_locate)
         return
      endif

      ! Set tolerances on |fk| and |xk(2)-xk(1)|

      dfk     = abs(its%r_b - its%r_a) / max(1d-6, abs(its%x_b - its%x_a))
      tol_fk  = max(1d-12*dfk, wtd%solv%eps * abs(ftarg))
      tol_xk  = wtd%solv%eps

      ! initialize the Brent loop for obtaining the prescribed vertical force

      k  = its%numit - 1
      used_bisec = .true.

      ! check convergence

      res_fk = its%r_b
      ldone = (abs(res_fk).lt.tol_fk .or. abs(its%x_b-its%x_a).lt.tol_xk .or. k.ge.wtd%solv%maxnr .or.  &
               wtd%ic%return.ge.2)

      ! while not "done" do

      do while (.not.ldone .and. my_ierror.eq.0)

         ! increment iteration number

         k = k + 1
         wtd%meta%itforc_inn = k

         if (x_force.ge.3) then
            write(bufout,'(a,i3)')   ' wr_contact_fz_brent: starting iteration k=',k
            call write_log(1, bufout)
         endif

         ! determine new point x_new

         call brent_set_xnew(k, its, dfz_dzws, used_bisec, tol_xk, x_new, x_force)

         ! compute problem with new xk, get new fk

         wtd%ws%z = x_new
         call wr_contact_pos(wtd, x_locate, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

         r_new = wtd%ftrk%z() - ftarg

         call brent_its_add_iterate(its, k, x_new, r_new, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, x_force)

         ! check convergence

         res_fk = its%r_b
         ldone = (abs(res_fk).lt.tol_fk .or. abs(its%x_b-its%x_a).lt.tol_xk .or. k.ge.wtd%solv%maxnr    &
                        .or. my_ierror.ne.0 .or. brent_its_has_jump(its, x_force) )

         if (my_ierror.lt.0) then
            write(bufout,'(a,i3,a)') ' An error occurred the in Brent algorithm (',my_ierror,           &
                '), aborting contact solution'
            call write_log(1, bufout)
         endif
      enddo  ! while not done: Brent loop

      ! compute sensitivity at final solution

      wtd%dfz_dzws = brent_sensitivity_k(k, its, x_force)

      if (my_ierror.eq.0 .and. abs(res_fk).ge.tol_fk .and. abs(its%x_b-its%x_a).ge.tol_xk) then
         write(bufout,'(a,i3,a)') ' Brent algorithm stopped with residual > tolerance.'
         call write_log(1, bufout)
         my_ierror = -2
      endif

   end subroutine wr_contact_fz_brent

!------------------------------------------------------------------------------------------------------------

   subroutine wr_contact_fx_secant(wtd, my_ierror)
!--purpose: process the w/r problem for a case with prescribed forces with Secant for Fx, Brent for Fz
      implicit none
!--subroutine arguments:
      type(t_ws_track)              :: wtd
      integer,          intent(out) :: my_ierror
!--local variables:
      integer,       parameter :: neq = 2
      logical                  :: has_fz, has_fx, ldone, laccept
      integer                  :: k, iestim_br, x_force, x_locate, itry, maxtry, sub_ierror
      real(kind=8)             :: aa, aob, bb, c11, c22, c23, dfx_domg, dfz_dzws, veloc, r_act, solv_eps
      real(kind=8)             :: xk(neq), xkm1(neq), xchr(neq), dxk(neq), fk(neq), fkm1(neq), dfk(neq), &
                                  ftarg(neq), res_fk(neq), tol_fk(neq), dfk_max(neq), bk(neq,neq),       &
                                  bkm1(neq,neq)
      type(t_brent_its)        :: its

      my_ierror = 0

      ! x_force = 0
      x_force   = wtd%ic%x_force
      x_locate  = wtd%ic%x_locate
      if (max(x_force,x_locate).ge.2) call write_log(' --- Start subroutine wr_contact_fx_secant ---')

      has_fz    = (wtd%ic%norm.ge.1)
      has_fx    = (wtd%ic%tang.ge.1 .and. (wtd%ic%force1.eq.1 .or. wtd%ic%force1.eq.2))

      if (.not.has_fz .or. .not.has_fx) then
         call write_log(' INTERNAL ERROR: secant method needs N=1, F=1')
         call abort_run()
      endif

      k = 0

      ! preparations for initial estimate for z-position

      if (has_fz) then

         k = k + 1
         if (x_force.ge.3) then
            write(bufout,'(a,i3)')   ' wr_contact_fx_secant: starting iteration k=',k
            call write_log(1, bufout)
         endif

         ! analyze geometry for overall minimum gap and curvatures, set initial estimate k=0 for Fz

         call wr_contact_init_hertz(wtd, its, dfz_dzws, aa, bb, IMETH_SECANT, x_force, x_locate,    &
                        sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

         ! Abort if the profiles have no overlap at all

         if (.not.wtd%ws%has_overlap) return
         
      endif

      if (x_force.ge.2) then
         write(bufout,'(3(a,g16.8))') ' fz_inp=', wtd%ws%fz_inp, ', dfz_dzws=', dfz_dzws
         call write_log(1, bufout)
      endif

      ! prepare initial estimate for pitch velocity omg_ws

      if (has_fx) then

         if (wtd%ic%is_roller()) then
            veloc = wtd%trk%vpitch_rol * wtd%trk%nom_radius
         else
            veloc = wtd%ws%vs
         endif

         ! using actual radius at the initial contact position, could use weighted average over all patches

         r_act = wtd%ws%nom_radius + wtd%ws%zw_min
         wtd%ws%vpitch = -veloc / r_act

         ! compute the Kalker coefficients at aob, poiss

         aob = aa / bb
         call linrol(wtd%mater%nu, aob, c11, c22, c23)

         ! compute stiffness according to linear theory
         ! fx = - gg * cc**2 * c11 * cksi

         dfx_domg  = wtd%mater%ga * aa * bb * c11 * wtd%ws%nom_radius / abs(veloc)

         ! set characteristic size of x-values: perturbation = xchr / 1000.

         xchr = 0.1d0 * veloc / wtd%ws%nom_radius

      endif

      ! Solve equations Fz(z_ws, omg_ws) / Fzinp       = 1              (1)
      !                 Fx(z_ws, omg_ws) / Fz(zws,omg) = Fxinp/Fzinp    (2)

      ! Using Secant-method for (2), with Brent-method as inner iteration to solve (1)

      ! target forces

      ftarg(1) = wtd%ws%fz_inp / wtd%ws%fz_inp
      ftarg(2) = wtd%ws%fx_inp / wtd%ws%fz_inp

      ! tolerances. max. resolution for Fz: z_ws \pm 1e-14, Fz \pm 1e-14 * dfz_dzws

      tol_fk(1) = max(1d-12*dfz_dzws/wtd%ws%fz_inp, wtd%solv%eps)
      tol_fk(2) = wtd%solv%eps * wtd%fric%fstat_min()

      ! maximum change in forces per iteration

      dfk_max(1) =  99d0
      dfk_max(2) = 1.0d0 * wtd%fric%fstat_min()

      if (x_force.ge.2) then
         write(bufout,'(2(a,g16.6))') ' tol_fk=',tol_fk(2),', dfk_max=',dfk_max(2)
         call write_log(1, bufout)
      endif

      ! initialize the Secant loop for obtaining the prescribed longitudinal force

      ! initial estimate [ z_ws, omg_ws ]

      xk(1) = wtd%ws%z
      xk(2) = wtd%ws%vpitch
      solv_eps  = wtd%solv%eps

      ! solve contact problem for the initial estimate

      iestim_br = 0
      wtd%solv%eps = max(1d-2, solv_eps)
      wtd%meta%itforc_out = k
      call wr_contact_fz_brent(wtd, iestim_br, 0d0, sub_ierror)
      if (my_ierror.eq.0 .and. sub_ierror.ne.-2) my_ierror = sub_ierror    ! ignore -2: |res|>tol
      wtd%solv%eps = solv_eps

      fk(1) = wtd%ftrk%z() / wtd%ws%fz_inp
      fk(2) = wtd%fws%x()  / wtd%ftrk%z()       ! fx_rel @ current fz
      xk(1) = wtd%ws%z

      ! set approximate Jacobian matrix B_0

      bk      = 0d0
      bk(2,2) = dfx_domg / wtd%ws%fz_inp

      ! check convergence

      res_fk = fk - ftarg

      ldone = (all(abs(res_fk).lt.tol_fk) .or. k.ge.wtd%solv%maxnr .or. wtd%ic%return.ge.2)

      if (x_force.ge.2) then
         write(bufout, '(3(a,g16.8))') ' fk=',fk(2),', res=',res_fk(2), ', tol=',tol_fk(2)
         call write_log(1, bufout)
      endif

      ! print output on Secant process

      call secant_fx_print(neq, k, xk, fk, ftarg, bk, wtd%ic, x_force)

      ! while not "done" do

      do while (.not.ldone .and. my_ierror.eq.0)

         ! increment iteration number

         k = k + 1

         if (k.eq.28) x_force = 3

         if (x_force.ge.3) then
            write(bufout,'(a,i3)')   ' wr_contact_fx_secant: starting iteration k=',k
            call write_log(1, bufout)
         endif

         ! if (x_force.ge.2) stop

         ! cycle previous values

         xkm1 = xk
         fkm1 = fk
         bkm1 = bk

         ! set new pitch velocity

         dxk(1) = 0d0
         dxk(2) = -res_fk(2) / bkm1(2,2)
         xk     = xkm1 + dxk

         ! Use back-tracking when dF becomes too large

         itry    =  0
         maxtry  = 10
         laccept = .false.

         do while(.not.laccept .and. my_ierror.eq.0)

            itry = itry + 1

            ! reduce step-size dx on consecutive inner iterations

            if (itry.gt.1) then
               dxk = 0.5d0 * dxk
               if (x_force.ge.2) then
                  write(bufout,'(2(a,f9.4))') ' ...halved step dom_ws=',dxk(2)
                  call write_log(1, bufout)
               endif
            endif

            xk   = xkm1 + dxk

            ! compute problem with new xk, get new fk

            wtd%ws%vpitch = xk(2)

            iestim_br = 1
            wtd%solv%eps = max(res_fk(2)/5d0, solv_eps)
            wtd%meta%itforc_out = k

            call wr_contact_fz_brent(wtd, iestim_br, 0d0, sub_ierror)
            if (my_ierror.eq.0 .and. sub_ierror.ne.-2) my_ierror = sub_ierror    ! ignore -2: |res|>tol

            wtd%solv%eps = solv_eps

            fk(1) = wtd%ftrk%z() / wtd%ws%fz_inp
            fk(2) = wtd%fws%x()  / wtd%ftrk%z()       ! fx_rel @ current fz
            xk(1) = wtd%ws%z
            dfk   = fk - fkm1

            ! check solution for fx

            if (any(isnan(dfk))) my_ierror = -1

            laccept = (itry.ge.maxtry .or. all(abs(dfk).lt.dfk_max))

            if (any(isnan(dfk))) then
               write(bufout,'(2(a,f11.3),a)') ' ...NaN-values found (', dfk(1), ',', dfk(2),            &
                        '), aborting'
               call write_log(1, bufout)
            elseif (.not.laccept) then
               write(bufout,'(2(a,f11.3),a)') ' ...update too large (', dfk(1), ',', dfk(2),            &
                     '), rejecting step'
               call write_log(1, bufout)
            endif
         enddo

         ! compute/update approximate Jacobian matrix

         if (abs(fk(2)-fkm1(2)).gt.1d3*tol_fk(2)) then
            bk(2,2) = (fk(2) - fkm1(2)) / dxk(2)
         endif

         ! check convergence

         res_fk = fk - ftarg
         ldone = (all(abs(res_fk).lt.tol_fk) .or. k.ge.wtd%solv%maxnr)

         ! print output on Secant process

         call secant_fx_print(neq, k, xk, fk, ftarg, bk, wtd%ic, x_force)

         if (my_ierror.lt.0) then
            write(bufout,'(a,i3,a)') ' An error occurred the in Secant algorithm (',my_ierror,          &
                '), aborting contact solution'
            call write_log(1, bufout)
         endif
      enddo  ! while not done: Secant loop

      if (my_ierror.eq.0 .and. any(abs(res_fk).ge.tol_fk)) my_ierror = -2

   end subroutine wr_contact_fx_secant

!------------------------------------------------------------------------------------------------------------

   subroutine secant_fx_print(neq, k, xk, fk, ftarg, bk, ic, x_force)
!--purpose: print information on Secant process for Fx
      implicit none
!--subroutine arguments:
      type(t_ic)               :: ic
      integer                  :: neq, k, x_force
      real(kind=8)             :: xk(neq), fk(neq), ftarg(neq), bk(neq,neq)
!--local variables:
      ! character(len=16)        :: str16((2+neq)*neq)

      ! horizontal and vertical forces prescribed: {F_z,F_x} = F( z_ws, omg_ws )

      if (x_force.ge.2 .or. ic%flow.ge.1) then
         ! str16(1) = fmt_gs(12,6, xk(1))
         ! str16(2) = fmt_gs(13,8, xk(2))
         ! write(str16(3), '(f13.8)') fk(1) - ftarg(1)
         ! write(str16(4), '(f13.8)') fk(2) - ftarg(2)
         ! write(str16(5), '(f11.6)') bk(2,2)
         write(bufout,8000) k, xk(1), xk(2), fk(2)-ftarg(2), bk(2,2)
         call write_log(1, bufout)
      endif

 8000 format(2x, i4,', SC, z, omg_ws:',f12.6,f13.8,', res: ',f13.8,', dFx/domg: ',f11.6)

   end subroutine secant_fx_print

!------------------------------------------------------------------------------------------------------------

   subroutine wr_contact_fy_brent(wtd, my_ierror)
!--purpose: process the w/r problem for a case with prescribed forces with Brent for Fy, Brent for Fz
      implicit none
!--subroutine arguments:
      type(t_ws_track)              :: wtd
      integer,          intent(out) :: my_ierror
!--local variables:
      integer, parameter       :: maxit = 100
      logical                  :: has_fz, has_fy, used_bisec, ldone
      integer                  :: k, x_force, x_locate, iestim_br, sub_ierror
      real(kind=8)             :: ftarg, fy_tot, dfy_dy, dz_ws, x_new, r_new, tol_xk, tol_fk, res_fk,   &
                                  eps_out, eps_inn
      type(t_brent_its)        :: its

      my_ierror = 0
      x_force  = wtd%ic%x_force
      x_locate = wtd%ic%x_locate

      if (x_force.ge.2) call write_log(' --- Start subroutine wr_contact_fy_brent ---')

      has_fz    = (wtd%ic%norm.ge.1)
      has_fy    = (wtd%ic%tang.ge.1 .and. wtd%ic%force1.eq.3)

      if (.not.has_fz .or. .not.has_fy) then
         call write_log(' INTERNAL ERROR: Brent method for Fy needs N=1, F=3')
         call abort_run()
      endif

      ! Solve equations F(xk) = fk = ftarg,  Fk : contact force on rail
      !    (1)   Fy(y_ws+dy, z_ws+dz) - ky * dy = -Fy_rail   -- outer iteration solving dy (this sub)
      !    (2)   Fz(y_ws+dy, z_ws+dz) - kz * dz = -Fz_rail   -- inner iteration solving dz (fz_brent)
      ! y_ws, z_ws considered fixed parameters, unknowns are the rail deflections dy, dz
      ! F[yz]_rail: spring force on rail

      ! if N=1: override vertical inputs

      if (wtd%ic%norm.eq.1) then
         wtd%trk%fz_rail = -wtd%ws%fz_inp
         wtd%trk%kz_rail =  0d0
      endif

      ! prepare structure to hold iterates with target force -fy_rail

      ftarg = -wtd%trk%fy_rail
      if (wtd%ic%is_left_side()) ftarg = -ftarg      ! mirroring left-side input to right-side internal

      call brent_its_init(its, ikYDIR, maxit, ftarg)

      ! set initial estimate dy = 0 and solve for z_ws using inner iteration

      k               = 0
      wtd%meta%itforc_out = k
      iestim_br       = 0
      wtd%trk%dy_defl = 0d0
      eps_out         = wtd%solv%eps
      eps_inn         = wtd%solv%eps
      wtd%solv%eps    = eps_inn

      call wr_contact_fz_brent(wtd, iestim_br, 0d0, sub_ierror)
      if (my_ierror.eq.0 .and. sub_ierror.ne.-2) my_ierror = sub_ierror    ! ignore -2: |res|>tol

      x_new  = wtd%trk%dy_defl
      fy_tot = wtd%ftrk%y() - wtd%trk%ky_rail * wtd%trk%dy_defl
      r_new  = fy_tot - ftarg
      dz_ws  = wtd%ws%z - wtd%ws%z_cnt0

      call brent_its_add_iterate(its, k, x_new, r_new, sub_ierror)
      call brent_its_print(k, its, wtd%ic, x_force)

      ! set second guess dy = 1 and solve for z_ws using inner iteration

      if (r_new.gt.0d0) then
         wtd%trk%dy_defl =  1d0
      else
         wtd%trk%dy_defl = -1d0
      endif

      k         = 1
      wtd%meta%itforc_out = k
      iestim_br = 0
      call wr_contact_fz_brent(wtd, iestim_br, dz_ws, sub_ierror)
      if (my_ierror.eq.0 .and. sub_ierror.ne.-2) my_ierror = sub_ierror    ! ignore -2: |res|>tol

      x_new  = wtd%trk%dy_defl
      fy_tot = wtd%ftrk%y() - wtd%trk%ky_rail * wtd%trk%dy_defl
      r_new  = fy_tot - ftarg
      dz_ws  = wtd%ws%z - wtd%ws%z_cnt0

      call brent_its_add_iterate(its, k, x_new, r_new, sub_ierror)
      call brent_its_print(k, its, wtd%ic, x_force)

      ! while not "done" do

      ldone   = .false.
      tol_xk  = eps_out
      tol_fk  = max(1d-6, eps_out * abs(ftarg))
      used_bisec = .false.

      do while (.not.ldone .and. my_ierror.eq.0)

         ! increment iteration number

         k = k + 1
         wtd%meta%itforc_out = k

         if (x_force.ge.3) then
            write(bufout,'(a,i3)')   ' wr_contact_fy_brent: starting iteration k=',k
            call write_log(1, bufout)
         endif

         ! determine new point x_new

         dfy_dy = brent_sensitivity_k(k-1, its, x_force)

         call brent_set_xnew(k, its, dfy_dy, used_bisec, tol_xk, x_new, x_force)

         ! solve for z_ws using inner iteration, get new fk

         iestim_br       = 0
         if (abs(x_new-wtd%trk%dy_defl).lt.0.1d0) iestim_br = 2
         wtd%trk%dy_defl = x_new
         eps_inn         = eps_out
         wtd%solv%eps    = eps_inn

         call wr_contact_fz_brent(wtd, iestim_br, dz_ws, sub_ierror)
         if (my_ierror.eq.0 .and. sub_ierror.ne.-2) my_ierror = sub_ierror    ! ignore -2: |res|>tol

         fy_tot = wtd%ftrk%y() - wtd%trk%ky_rail * wtd%trk%dy_defl
         r_new  = fy_tot - ftarg
         dz_ws  = wtd%ws%z - wtd%ws%z_cnt0

         call brent_its_add_iterate(its, k, x_new, r_new, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, x_force)
         ! write(bufout,'(a,f12.4)') ' dz_ws=',dz_ws
         ! call write_log(1, bufout)

         ! check convergence

         res_fk = its%r_b
         ldone = (abs(res_fk).lt.tol_fk .or. abs(its%x_b-its%x_a).lt.tol_xk .or. k.ge.wtd%solv%maxnr    &
                        .or. my_ierror.ne.0 .or. brent_its_has_jump(its, x_force) )

         if (my_ierror.lt.0) then
            write(bufout,'(a,i3,a)') ' An error occurred the in Brent algorithm (',my_ierror,           &
                '), aborting contact solution'
            call write_log(1, bufout)
         endif
      enddo  ! while not done: Brent loop

      if (my_ierror.eq.0 .and. abs(res_fk).ge.tol_fk .and. abs(its%x_b-its%x_a).ge.tol_xk) then
         write(bufout,'(a,i3,a)') ' Brent algorithm (fy) stopped with residual > tolerance.'
         call write_log(1, bufout)
         my_ierror = -2
      endif

   end subroutine wr_contact_fy_brent

!------------------------------------------------------------------------------------------------------------

end module m_wr_totforce

