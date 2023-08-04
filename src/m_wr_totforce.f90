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

private wr_contact_fyz_broyden

contains

!------------------------------------------------------------------------------------------------------------

   subroutine wr_contact(wtd, my_ierror)
!--purpose: process the whole w/r problem for a case, including NR-loops for prescribed forces
      implicit none
!--subroutine arguments:
      type(t_ws_track), target      :: wtd
      integer,          intent(out) :: my_ierror
!--local variables:
      logical            :: lfound, is_left
      integer            :: iestim_br, idebug_cnt, icp, i, sub_ierror, known_err(10)

      my_ierror = 0

      ! idebug: 0 = errors/silent, 1 = warnings/info, >=2 = flow/debug

      if (wtd%ic%x_locate.le.0) then
         idebug_cnt = wtd%ic%ilvout
      else
         idebug_cnt = wtd%ic%x_locate
      endif

      ! set rail deflections if not used in this case

      if (wtd%ic%force1.ne.3 .and. wtd%ic%force1.ne.4) then
         wtd%trk%dy_defl = 0d0
         wtd%trk%dz_defl = 0d0
      endif

      ! switch between solvers dependent on prescribed total forces

      if (wtd%ic%norm.le.0) then

         ! N=0, F=0: solve w/r problem with positions & velocities prescribed

         if (idebug_cnt.ge.3) call write_log(' wr_contact: solving for given position...')
         wtd%meta%itforc_out = 0
         wtd%meta%itforc_inn = 0
         call wr_contact_pos(wtd, idebug_cnt, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

      elseif (wtd%ic%force1.le.0) then

         ! N=1, F=0: solve w/r problem with prescribed vertical force

         if (idebug_cnt.ge.3) call write_log(' wr_contact: solving for given vertical force...')
         iestim_br = 0
         call wr_contact_fz_brent(wtd, iestim_br, 0d0, idebug_cnt, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

      elseif (wtd%ic%force1.eq.1 .or. wtd%ic%force1.eq.2) then

         ! N=1, F=1,2: solve w/r problem with prescribed vertical and longitudinal forces

         if (idebug_cnt.ge.3) call write_log(' wr_contact: solving for given vert./long. forces (secant)...')
         call wr_contact_fx_secant(wtd, idebug_cnt, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

      elseif (wtd%ic%force1.eq.4) then

         ! N=1, F=3: solve w/r problem with prescribed vertical and lateral forces - Broyden method

         if (idebug_cnt.ge.3) call write_log(' wr_contact: solving for given vert./lat. forces (Broyden)...')
         call wr_contact_fyz_broyden(wtd, idebug_cnt, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

      elseif (wtd%ic%force1.eq.3) then

         ! N=1, F=3: solve w/r problem with prescribed vertical and lateral forces - Brent method

         if (idebug_cnt.ge.3) call write_log(' wr_contact: solving for given vert./lat. forces (brent)...')
         call wr_contact_fy_brent(wtd, idebug_cnt, sub_ierror)
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
         if (idebug_cnt.ge.5) call write_log(' ...wr_contact: calling writmt')
         do icp = 1, wtd%numcps
            associate( gd => wtd%allcps(icp)%cp%gd )
            gd%ic%matfil_surf = wtd%ic%matfil_surf
            is_left = (wtd%ic%config.eq.0 .or. wtd%ic%config.eq.4)
            call writmt (gd%meta, gd%ic, gd%cgrid, gd%potcon%xl, gd%potcon%yl, gd%geom%hs1, gd%mater,   &
                         gd%fric, gd%kin, gd%outpt1, is_left)
            end associate
         enddo
      endif

      ! compute subsurface stresses when requested

      if (wtd%ic%stress.ge.1) then
         call wr_subsurf(wtd, idebug_cnt)
      endif

      ! write output when R=0 or 1

      if (wtd%ic%return.le.1) call wr_output(wtd)

   end subroutine wr_contact

!------------------------------------------------------------------------------------------------------------

   subroutine wr_contact_init_hertz(wtd, its, dfz_dzws, aa, bb, imeth, idebug_br, idebug_cnt, my_ierror)
!--purpose: determine an initial bracket for z_ws for starting the Brent algorithm
      implicit none
!--subroutine arguments:
      type(t_ws_track),  target        :: wtd
      integer,           intent(in)    :: imeth
      integer,           intent(in)    :: idebug_br, idebug_cnt
      integer,           intent(out)   :: my_ierror
      real(kind=8),      intent(out)   :: dfz_dzws, aa, bb
      type(t_brent_its), intent(inout) :: its
!--local variables:
      integer                  :: k, ic_norm, ipotcn, ic_return, sub_ierror
      real(kind=8)             :: cp, rho, pen, e_star, epshz, azz, x_0, x_new, r_new

      if (idebug_br.ge.3) call write_log(' --- Start subroutine wr_contact_init_hertz ---')

      my_ierror = 0

      ! initial estimate for z-position

      k   = 0
      wtd%meta%itforc_inn = k
      wtd%ws%z = 0d0

      ! analyze geometrical problem for the initial z-position, without actual solving

      ic_return = wtd%ic%return
      wtd%ic%return = 3

      call wr_contact_pos(wtd, idebug_cnt, sub_ierror)
      if (my_ierror.eq.0) my_ierror = sub_ierror

      wtd%ic%return = ic_return

      if (max(idebug_br,idebug_cnt).ge.3) then
         write(bufout,'(3(a,g14.6))') ' z_ws =',wtd%ws%z
         call write_log(1, bufout)
         write(bufout,'(3(a,g14.6))') ' Overall minimum gap=', wtd%ws%gap_min,', a1=',wtd%ws%a1,     &
                ', b1=',wtd%ws%b1
         call write_log(1, bufout)
      endif

      ! Abort in case of an error

      if (my_ierror.ne.0) then
         wtd%numcps = 0
         if (max(idebug_br,idebug_cnt).ge.1) then
            write(bufout,'(a,i6,a)') ' An error occurred (',my_ierror,'), skipping computation.'
            call write_log(1, bufout)
         endif
         return
      endif

      ! Abort if the profiles have no overlap at all

      if (.not.wtd%ws%has_overlap) then
         wtd%numcps = 0
         if (max(idebug_br,idebug_cnt).ge.1) then
            call write_log(' The profiles have no overlap at all, skipping computation.')
         endif
         return
      endif
         
      ! store initial point { x_0, r_0 }, print output

      x_0 = wtd%ws%z + wtd%ws%gap_min

      if (imeth.eq.IMETH_BRENT) then

         call brent_its_add_iterate(its, k, x_0, -wtd%ws%fz_inp, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, idebug_br)

      elseif (imeth.eq.IMETH_BROYDN) then

         if (idebug_br.ge.2 .or. wtd%ic%flow.ge.1) then
            write(bufout,7020) -1, x_0, -1d0
            call write_log(1, bufout)
         endif
 7020    format(4x,i6,', NR,  z_ws, Fz: ',f16.6, f16.8)

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

         if (idebug_br.ge.3) then
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

         call wr_contact_pos(wtd, idebug_cnt, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

         r_new    = wtd%ftrk%z() - wtd%ws%fz_inp
         dfz_dzws = 1.5d0 * wtd%ftrk%z() / (x_new - x_0)

         ! write(bufout,*) 'dfz_dzws= =', dfz_dzws
         ! call write_log(1, bufout)

         call brent_its_add_iterate(its, k, x_new, r_new, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, idebug_br)

      endif
         
   end subroutine wr_contact_init_hertz

!------------------------------------------------------------------------------------------------------------

   subroutine wr_contact_fz_brent(wtd, iestim_br, z_offset, idebug_cnt, my_ierror)
!--purpose: process the w/r problem for a case with prescribed vertical force with Brents algorithm
      implicit none
!--subroutine arguments:
      type(t_ws_track), target      :: wtd
      integer,          intent(in)  :: iestim_br, idebug_cnt
      real(kind=8),     intent(in)  :: z_offset
      integer,          intent(out) :: my_ierror
!--local variables:
      integer, parameter       :: maxit = 100
      logical                  :: used_bisec, ldone
      integer                  :: k, idebug_br, ic_return, sub_ierror
      real(kind=8)             :: ftarg, res_fk, tol_fk, dfk, tol_xk, x_0, x_new, r_new, dfz_dzws, aa, bb
      type(t_brent_its)        :: its

      my_ierror = 0
      idebug_br = 2
      ! if (wtd%meta%itforc_out.eq.2) idebug_br = 3

      if (idebug_cnt.ge.2) call write_log(' --- Start subroutine wr_contact_fz_brent ---')

      ! Abort if a zero or negative total force is prescribed

      if (wtd%ws%fz_inp.lt.1d-9) then
         wtd%numcps = 0
         if (idebug_cnt.ge.1) then  
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

      if (iestim_br.eq.2) then

         ! analyze geometrical problem for the initial z-position, without actual solving

         wtd%ws%z  = 0d0
         ic_return = wtd%ic%return
         wtd%ic%return = 3

         call wr_contact_pos(wtd, idebug_cnt, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         wtd%ic%return = ic_return

         wtd%ws%z  = wtd%ws%z_cnt0 + z_offset
      endif

      if (iestim_br.ge.1) then

         ! methods 1, 2: use initial estimate from ws%z, wtd%dfz_dzws

         wtd%meta%itforc_inn = 1
         call wr_contact_pos(wtd, idebug_cnt, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

         ! store initial point { x_0, r_0 } at which initial contact occurs

         k   = 0
         x_0 = wtd%ws%z_cnt0
         call brent_its_add_iterate(its, k, x_0, -wtd%ws%fz_inp, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, idebug_br)

         ! store first point { x_1, r_1 } for initial estimate ws%z

         k   = 1
         wtd%meta%itforc_inn = k
         x_new    = wtd%ws%z
         r_new    = wtd%ftrk%z() - wtd%ws%fz_inp
         dfz_dzws = wtd%dfz_dzws

         call brent_its_add_iterate(its, k, x_new, r_new, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, idebug_br)

      else

         ! analyze geometry for overall minimum gap and curvatures, set initial estimates k=0, k=1

         call wr_contact_init_hertz(wtd, its, dfz_dzws, aa, bb, IMETH_BRENT, idebug_br, idebug_cnt,     &
                        sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

      endif

      ! Abort in case of error or if the profiles have no overlap at all

      if (my_ierror.ne.0 .or. .not.wtd%ws%has_overlap) then
         call write_log(' no overlap...')
         call total_forces_moments(wtd, idebug_cnt)
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

         if (idebug_br.ge.3) then
            write(bufout,'(a,i3)')   ' wr_contact_fz_brent: starting iteration k=',k
            call write_log(1, bufout)
         endif

         ! determine new point x_new

         call brent_set_xnew(k, its, dfz_dzws, used_bisec, tol_xk, x_new, idebug_br)

         ! compute problem with new xk, get new fk

         wtd%ws%z = x_new
         call wr_contact_pos(wtd, idebug_cnt, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

         r_new = wtd%ftrk%z() - ftarg

         call brent_its_add_iterate(its, k, x_new, r_new, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, idebug_br)

         ! check convergence

         res_fk = its%r_b
         ldone = (abs(res_fk).lt.tol_fk .or. abs(its%x_b-its%x_a).lt.tol_xk .or. k.ge.wtd%solv%maxnr    &
                        .or. my_ierror.ne.0 .or. brent_its_has_jump(its, idebug_br) )

         if (my_ierror.lt.0) then
            write(bufout,'(a,i3,a)') ' An error occurred the in Brent algorithm (',my_ierror,           &
                '), aborting contact solution'
            call write_log(1, bufout)
         endif
      enddo  ! while not done: Brent loop

      ! compute sensitivity at final solution

      wtd%dfz_dzws = brent_sensitivity_k(k, its, idebug_br)

      if (my_ierror.eq.0 .and. abs(res_fk).ge.tol_fk .and. abs(its%x_b-its%x_a).ge.tol_xk) then
         write(bufout,'(a,i3,a)') ' Brent algorithm stopped with residual > tolerance.'
         call write_log(1, bufout)
         my_ierror = -2
      endif

   end subroutine wr_contact_fz_brent

!------------------------------------------------------------------------------------------------------------

   subroutine wr_contact_fx_secant(wtd, idebug_cnt, my_ierror)
!--purpose: process the w/r problem for a case with prescribed forces with Secant for Fx, Brent for Fz
      implicit none
!--subroutine arguments:
      type(t_ws_track), target      :: wtd
      integer,          intent(in)  :: idebug_cnt
      integer,          intent(out) :: my_ierror
!--local variables:
      integer,       parameter :: neq = 2
      logical                  :: has_fz, has_fx, is_roller, ldone, laccept
      integer                  :: k, iestim_br, idebug_br, itry, maxtry, sub_ierror
      real(kind=8)             :: aa, aob, bb, c11, c22, c23, dfx_domg, dfz_dzws, veloc, r_act, solv_eps
      real(kind=8)             :: xk(neq), xkm1(neq), xchr(neq), dxk(neq), fk(neq), fkm1(neq), dfk(neq), &
                                  ftarg(neq), res_fk(neq), tol_fk(neq), dfk_max(neq), bk(neq,neq),       &
                                  bkm1(neq,neq)
      type(t_brent_its)        :: its

      my_ierror = 0

      idebug_br = idebug_cnt
      idebug_br = 0
      if (idebug_br.ge.2) call write_log(' --- Start subroutine wr_contact_fx_secant ---')

      has_fz    = (wtd%ic%norm.ge.1)
      has_fx    = (wtd%ic%tang.ge.1 .and. (wtd%ic%force1.eq.1 .or. wtd%ic%force1.eq.2))
      is_roller = (wtd%ic%config.ge.4)

      if (.not.has_fz .or. .not.has_fx) then
         call write_log(' INTERNAL ERROR: secant method needs N=1, F=1')
         call abort_run()
      endif

      ! preparations for initial estimate for z-position

      if (has_fz) then

         ! analyze geometry for overall minimum gap and curvatures, set initial estimate k=0 for Fz

         call wr_contact_init_hertz(wtd, its, dfz_dzws, aa, bb, IMETH_SECANT, idebug_br, idebug_cnt,    &
                        sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

         ! Abort if the profiles have no overlap at all

         if (.not.wtd%ws%has_overlap) return
         
      endif

      if (idebug_br.ge.2) then
         write(bufout,'(3(a,g14.6))') ' fz_inp=', wtd%ws%fz_inp, ', dfz_dzws=', dfz_dzws
         call write_log(1, bufout)
      endif

      ! prepare initial estimate for pitch velocity omg_ws

      if (has_fx) then

         if (is_roller) then
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

      if (idebug_br.ge.2) then
         write(bufout,'(2(a,g14.4))') ' tol_fk=',tol_fk(2),', dfk_max=',dfk_max(2)
         call write_log(1, bufout)
      endif

      ! initialize the Secant loop for obtaining the prescribed longitudinal force

      ! initial estimate [ z_ws, omg_ws ]

      k     = 0
      xk(1) = wtd%ws%z
      xk(2) = wtd%ws%vpitch
      solv_eps  = wtd%solv%eps

      ! solve contact problem for the initial estimate

      iestim_br = 0
      wtd%solv%eps = max(1d-2, solv_eps)
      call wr_contact_fz_brent(wtd, iestim_br, 0d0, idebug_cnt, sub_ierror)
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

      if (idebug_br.ge.2) then
         write(bufout, '(4(a,f11.3))') ' fk=',fk(2),', res=',res_fk(2), ', tol=',tol_fk(2)
         call write_log(1, bufout)
      endif

      ! print output on Secant process

      call secant_fx_print(neq, k, xk, fk, ftarg, bk, wtd%ic, idebug_br)

      ! while not "done" do

      do while (.not.ldone .and. my_ierror.eq.0)

         ! increment iteration number

         k = k + 1

         if (k.eq.28) idebug_br = 3

         if (idebug_br.ge.3) then
            write(bufout,'(a,i3)')   ' wr_contact_fx_secant: starting iteration k=',k
            call write_log(1, bufout)
         endif

         ! if (idebug_br.ge.2) stop

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
               if (idebug_br.ge.2) then
                  write(bufout,'(2(a,f9.4))') ' ...halved step dom_ws=',dxk(2)
                  call write_log(1, bufout)
               endif
            endif

            xk   = xkm1 + dxk

            ! compute problem with new xk, get new fk

            wtd%ws%vpitch = xk(2)

            iestim_br = 1
            wtd%solv%eps = max(res_fk(2)/5d0, solv_eps)

            call wr_contact_fz_brent(wtd, iestim_br, 0d0, idebug_cnt, sub_ierror)
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

         bk(2,2) = (fk(2) - fkm1(2)) / dxk(2)

         ! check convergence

         res_fk = fk - ftarg
         ldone = (all(abs(res_fk).lt.tol_fk) .or. k.ge.wtd%solv%maxnr)

         ! print output on Secant process

         call secant_fx_print(neq, k, xk, fk, ftarg, bk, wtd%ic, idebug_br)

         if (my_ierror.lt.0) then
            write(bufout,'(a,i3,a)') ' An error occurred the in Secant algorithm (',my_ierror,          &
                '), aborting contact solution'
            call write_log(1, bufout)
         endif
      enddo  ! while not done: Secant loop

      if (my_ierror.eq.0 .and. any(abs(res_fk).ge.tol_fk)) my_ierror = -2

   end subroutine wr_contact_fx_secant

!------------------------------------------------------------------------------------------------------------

   subroutine secant_fx_print(neq, k, xk, fk, ftarg, bk, ic, idebug_br)
!--purpose: print information on Secant process for Fx
      implicit none
!--subroutine arguments:
      type(t_ic)               :: ic
      integer                  :: neq, k, idebug_br
      real(kind=8)             :: xk(neq), fk(neq), ftarg(neq), bk(neq,neq)
!--local variables:
      ! character(len=16)        :: str16((2+neq)*neq)

      ! horizontal and vertical forces prescribed: {F_z,F_x} = F( z_ws, omg_ws )

      if (idebug_br.ge.2 .or. ic%flow.ge.1) then
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

   subroutine wr_contact_fy_brent(wtd, idebug_cnt, my_ierror)
!--purpose: process the w/r problem for a case with prescribed forces with Brent for Fy, Brent for Fz
      implicit none
!--subroutine arguments:
      type(t_ws_track), target      :: wtd
      integer,          intent(in)  :: idebug_cnt
      integer,          intent(out) :: my_ierror
!--local variables:
      integer, parameter       :: maxit = 100
      logical                  :: has_fz, has_fy, is_left_side, used_bisec, ldone
      integer                  :: k, idebug_br, iestim_br, sub_ierror
      real(kind=8)             :: ftarg, fy_tot, dfy_dy, dz_ws, x_new, r_new, tol_xk, tol_fk, res_fk,   &
                                  eps_out, eps_inn
      type(t_brent_its)        :: its

      my_ierror = 0
      idebug_br = idebug_cnt
      idebug_br = 1

      if (idebug_br.ge.2) call write_log(' --- Start subroutine wr_contact_fy_brent ---')

      has_fz    = (wtd%ic%norm.ge.1)
      has_fy    = (wtd%ic%tang.ge.1 .and. (wtd%ic%force1.eq.3 .or. wtd%ic%force1.eq.4))
      is_left_side = (wtd%ic%config.eq.0 .or. wtd%ic%config.eq.4)

      if (.not.has_fz .or. .not.has_fy) then
         call write_log(' INTERNAL ERROR: Brent method for Fy needs N=1, F=4')
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
      if (is_left_side) ftarg = -ftarg          ! mirroring left-side input to right-side internal

      call brent_its_init(its, ikYDIR, maxit, ftarg)

      ! set initial estimate dy = 0 and solve for z_ws using inner iteration

      k               = 0
      wtd%meta%itforc_out = k
      iestim_br       = 0
      wtd%trk%dy_defl = 0d0
      eps_out         = wtd%solv%eps
      eps_inn         = wtd%solv%eps
      wtd%solv%eps    = eps_inn

      call wr_contact_fz_brent(wtd, iestim_br, 0d0, idebug_cnt, sub_ierror)
      if (my_ierror.eq.0 .and. sub_ierror.ne.-2) my_ierror = sub_ierror    ! ignore -2: |res|>tol

      x_new  = wtd%trk%dy_defl
      fy_tot = wtd%ftrk%y() - wtd%trk%ky_rail * wtd%trk%dy_defl
      r_new  = fy_tot - ftarg
      dz_ws  = wtd%ws%z - wtd%ws%z_cnt0

      call brent_its_add_iterate(its, k, x_new, r_new, sub_ierror)
      call brent_its_print(k, its, wtd%ic, idebug_br)

      ! set second guess dy = 1 and solve for z_ws using inner iteration

      if (r_new.gt.0d0) then
         wtd%trk%dy_defl =  1d0
      else
         wtd%trk%dy_defl = -1d0
      endif

      k         = 1
      wtd%meta%itforc_out = k
      iestim_br = 0
      call wr_contact_fz_brent(wtd, iestim_br, dz_ws, idebug_cnt, sub_ierror)
      if (my_ierror.eq.0 .and. sub_ierror.ne.-2) my_ierror = sub_ierror    ! ignore -2: |res|>tol

      x_new  = wtd%trk%dy_defl
      fy_tot = wtd%ftrk%y() - wtd%trk%ky_rail * wtd%trk%dy_defl
      r_new  = fy_tot - ftarg
      dz_ws  = wtd%ws%z - wtd%ws%z_cnt0

      call brent_its_add_iterate(its, k, x_new, r_new, sub_ierror)
      call brent_its_print(k, its, wtd%ic, idebug_br)

      ! while not "done" do

      ldone   = .false.
      tol_xk  = eps_out
      tol_fk  = max(1d-6, eps_out * abs(ftarg))

      do while (.not.ldone .and. my_ierror.eq.0)

         ! increment iteration number

         k = k + 1
         wtd%meta%itforc_out = k

         if (idebug_br.ge.3) then
            write(bufout,'(a,i3)')   ' wr_contact_fy_brent: starting iteration k=',k
            call write_log(1, bufout)
         endif

         ! determine new point x_new

         dfy_dy = brent_sensitivity_k(k-1, its, idebug_br)

         call brent_set_xnew(k, its, dfy_dy, used_bisec, tol_xk, x_new, idebug_br)

         ! solve for z_ws using inner iteration, get new fk

         iestim_br       = 0
         if (abs(x_new-wtd%trk%dy_defl).lt.0.1d0) iestim_br = 2
         wtd%trk%dy_defl = x_new
         eps_inn         = eps_out
         wtd%solv%eps    = eps_inn

         call wr_contact_fz_brent(wtd, iestim_br, dz_ws, idebug_cnt, sub_ierror)
         if (my_ierror.eq.0 .and. sub_ierror.ne.-2) my_ierror = sub_ierror    ! ignore -2: |res|>tol

         fy_tot = wtd%ftrk%y() - wtd%trk%ky_rail * wtd%trk%dy_defl
         r_new  = fy_tot - ftarg
         dz_ws  = wtd%ws%z - wtd%ws%z_cnt0

         call brent_its_add_iterate(its, k, x_new, r_new, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, idebug_br)
         ! write(bufout,'(a,f12.4)') ' dz_ws=',dz_ws
         ! call write_log(1, bufout)

         ! check convergence

         res_fk = its%r_b
         ldone = (abs(res_fk).lt.tol_fk .or. abs(its%x_b-its%x_a).lt.tol_xk .or. k.ge.wtd%solv%maxnr    &
                        .or. my_ierror.ne.0 .or. brent_its_has_jump(its, idebug_br) )

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

   subroutine wr_contact_init_broyden(wtd, dfn_dpen, idebug_br, idebug_cnt, my_ierror)
!--purpose: determine initial position and Jacobian for the Broyden algorithm for massless rail deflection
      implicit none
!--subroutine arguments:
      type(t_ws_track),  target        :: wtd
      integer,           intent(in)    :: idebug_br, idebug_cnt
      integer,           intent(out)   :: my_ierror
      real(kind=8),      intent(out)   :: dfn_dpen
!--local variables:
      integer                  :: k, ic_norm, ipotcn, ic_return, sub_ierror
      real(kind=8)             :: aa, bb, cp, rho, pen, e_star, epshz, azz, fn_inp

      if (idebug_br.ge.3) call write_log(' --- Start subroutine wr_contact_init_broyden ---')

      my_ierror = 0

      ! initial estimate for z-position

      k        = 0
      wtd%meta%itforc_inn = k
      wtd%ws%z = 0d0
      wtd%ftrk = vec_zero()

      ! analyze geometrical problem for the initial z-position, without actual solving

      ic_return = wtd%ic%return
      wtd%ic%return = 3

      call wr_contact_pos(wtd, idebug_cnt, sub_ierror)
      if (my_ierror.eq.0) my_ierror = sub_ierror

      wtd%ic%return = ic_return

      if (max(idebug_br,idebug_cnt).ge.3) then
         write(bufout,'(4(a,g14.6))') ' Overall minimum gap=', wtd%ws%gap_min,', delta=',wtd%ws%delt_min, &
                ', a1=',wtd%ws%a1, ', b1=',wtd%ws%b1
         call write_log(1, bufout)
      endif

      ! Abort in case of an error

      if (my_ierror.ne.0) then
         wtd%numcps = 0
         if (max(idebug_br,idebug_cnt).ge.1) then
            write(bufout,'(a,i6,a)') ' An error occurred (',my_ierror,'), skipping computation.'
            call write_log(1, bufout)
         endif
         return
      endif

      ! Abort if the profiles have no overlap at all

      if (.not.wtd%ws%has_overlap) then
         wtd%numcps = 0
         if (max(idebug_br,idebug_cnt).ge.1) then
            call write_log(' The profiles have no overlap at all, skipping computation.')
         endif
         return
      endif
         
      ! Use the Hertzian approach to get a sensible initial z-position,
      ! solve the Hertzian problem for estimated curvatures and initial contact angle

      ic_norm =  1
      ipotcn  = -1
      e_star  = wtd%mater%ga / (1d0 - wtd%mater%nu)
      epshz   = wtd%solv%eps
         
      if (wtd%ic%norm.eq.1) wtd%trk%fz_rail = -wtd%ws%fz_inp
      fn_inp  = -wtd%trk%fz_rail / max(0.5d0,cos(wtd%ws%delt_min))

      call hzcalc3d (e_star, epshz, ipotcn, wtd%ws%a1, wtd%ws%b1, aa, bb, ic_norm, pen,                 &
                     fn_inp, cp, rho)

      ! estimate derivative dFn/dpen

      if (pen.gt.1d-6) then
         dfn_dpen = 1.5d0 * wtd%ws%fz_inp / pen
      else
         ! no contact: using influence coefficient Azz(0,0) for 1 element in contact

         call azz_one_element(wtd%ic, wtd%mater, wtd%discr%dx, wtd%discr%ds, azz)
         dfn_dpen = wtd%discr%dx * wtd%discr%ds * wtd%mater%ga / azz 
         write(bufout,*) 'gap_min=',wtd%ws%gap_min,', pen=',pen,', using dfn_dpen=',dfn_dpen
         call write_log(1, bufout)
      endif

   end subroutine wr_contact_init_broyden

!------------------------------------------------------------------------------------------------------------

   subroutine wr_contact_fyz_broyden(wtd, idebug_cnt, my_ierror)
!--purpose: process the w/r problem for a case with prescribed lat/vert forces (massless rail model)
!           with Broyden algorithm, using estimated Jacobian Bk
      implicit none
!--subroutine arguments:
      type(t_ws_track), target      :: wtd
      integer,          intent(in)  :: idebug_cnt
      integer,          intent(out) :: my_ierror
!--local variables:
      integer,       parameter :: neq = 2
      logical                  :: use_findiff_1st, use_findiff, is_left_side, ldone, laccept
      integer                  :: k, idebug_br, ntry, maxtry, sub_ierror
      real(kind=8)             :: dfn_dpen
      real(kind=8)             :: xk(neq), xkm1(neq), xchr(neq), dxk(neq), fk(neq), fkm1(neq), dfk(neq), &
                                  ftarg(neq), res_fk(neq), tol_fk(neq), dfk_max(neq), bk(neq,neq),       &
                                  bkm1(neq,neq)

      my_ierror = 0
      idebug_br = idebug_cnt
      ! idebug_br = 2
      if (idebug_br.ge.2) call write_log(' --- Start subroutine wr_contact_fyz ---')

      use_findiff_1st = .true.
      use_findiff     = .true.
      is_left_side    = (wtd%ic%config.eq.0 .or. wtd%ic%config.eq.4)

      ! if N=1: override vertical inputs

      if (wtd%ic%norm.eq.1) then
         wtd%trk%fz_rail = -wtd%ws%fz_inp
         wtd%trk%kz_rail =  0d0
      endif

      ! set target forces -fy_rail, -fz_rail

      ftarg(1) = -wtd%trk%fy_rail
      ftarg(2) = -wtd%trk%fz_rail
      if (is_left_side) ftarg(1) = -ftarg(1)     ! mirroring left-side input to right-side internal

      ! relative tolerances. max. resolution for dy,dz \pm 1e-12, Fy,Fz \pm Jac * 1e-12

      tol_fk(1:2) = wtd%solv%eps * max(ftarg(1), ftarg(2))

      ! maximum change in forces per iteration

      dfk_max(1:2) = 1.5d0 * max(abs(ftarg(1)),abs(ftarg(2)))

      ! initial estimates for outputs dy, dz

      wtd%trk%dy_defl = 0d0
      wtd%trk%dz_defl = 0d0

      ! characteristic sizes of outputs dy, dz

      xchr(1:2) = 0.05d0

      ! analyze geometry for overall minimum gap and curvatures, set initial estimates k=0

      call wr_contact_init_broyden(wtd, dfn_dpen, idebug_br, idebug_cnt, sub_ierror)
      if (my_ierror.eq.0) my_ierror = sub_ierror

      ! Abort if the profiles have no overlap at all

      if (.not.wtd%ws%has_overlap) return
         
      ! overwrite z_ws such that there's precise contact at dy = dz = 0

      if (wtd%ic%norm.eq.1) then
         wtd%ws%z = wtd%ws%z + wtd%ws%gap_min
      endif

      ! Broyden: solve equations F(xk) = fk = ftarg,  Fk : contact force on rail
      !     vector xk    == [ dy ; dz ]; 
      !     vector fk    == [ Fy(y_ws+dy, z_ws+dz) - ky*dy ; Fz(y_ws+dy, z_ws+dz) - kz*dz ]
      !     vector ftarg == [ -Fy_rail ; -Fz_rail ],    F[yz]_rail: spring force on rail

      ! initial estimate [ dy, dz ] = [ 0, 0 ]

      k  = 0
      xk = (/ 0d0 , 0d0 /)
      fk = (/ 0d0 , 0d0 /)

      ! solve contact problem for the initial estimate
      !if (wtd%ic%norm.eq.0) then
      !   call broyden_set_state(neq, xk, wtd, idebug_br)
      !   call wr_contact_pos(wtd, idebug_cnt, sub_ierror)
      !   if (my_ierror.eq.0) my_ierror = sub_ierror
      !   call broyden_get_state(wtd, neq, fk, idebug_br)
      !endif

      ! set approximate Jacobian matrix B_0

      bk      =  0d0
      bk(1,1) =   sin(wtd%ws%delt_min) * dfn_dpen - wtd%trk%ky_rail
      bk(2,2) =  -cos(wtd%ws%delt_min) * dfn_dpen - wtd%trk%ky_rail

      ! lower the wheel-set by pen to get near the desired approach

      ! k     = 1
      ! wtd%meta%itforc_inn = k
      ! wtd%trk%dy_defl = pen * sin(wtd%ws%delt_min)
      ! wtd%trk%dz_defl = pen * cos(wtd%ws%delt_min)

      ! check convergence

      res_fk = fk - ftarg

      ldone = (all(abs(res_fk).lt.tol_fk) .or. k.ge.wtd%solv%maxnr .or. wtd%ic%return.ge.2)

      ! print output on Broyden process

      ntry = 1
      call broyden_print(neq, k, xk, fk, ftarg, bk, wtd%ic, ntry, idebug_br)

      ! while not "done" do

      do while (.not.ldone .and. my_ierror.eq.0)

         ! increment iteration number

         k = k + 1

         if (k.eq.28) idebug_br = 3

         if (idebug_br.ge.3) then
            write(bufout,'(a,i3)')   ' wr_contact_fyz: starting iteration k=',k
            call write_log(1, bufout)
         endif

         ! if (idebug_br.ge.2) stop

         ! cycle previous values

         xkm1 = xk
         fkm1 = fk
         bkm1 = bk

         ! compute increment dxk and new iterand xk

         call broyden_solve_dxk(neq, bkm1, ftarg, fk, res_fk, dxk, idebug_br, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

         ! Solve contact problem, using backtracking for dy,dz if dF is too large or contact is lost

         ntry =  0
         maxtry = 10
         laccept = .false.

         do while(.not.laccept .and. my_ierror.eq.0)

            ntry = ntry + 1

            ! reduce step-size dx on consecutive inner iterations

            if (ntry.gt.1) then
               dxk = 0.5d0 * dxk
               if (idebug_br.ge.2) then
                  write(bufout,'(2(a,f8.4))') ' ...halved step dy_defl=',dxk(1),', dz_defl=',dxk(2)
                  call write_log(1, bufout)
               endif
            endif

            xk   = xkm1 + dxk

            if (idebug_br.ge.3) then
               write(bufout, '(2(a,f12.6),a)') ' Broyden: xk     = [', xk(1),    ',', xk(2),    ']^T'
               call write_log(1, bufout)
            endif

            ! compute problem with new xk, get new fk

            call broyden_set_state(neq, xk, wtd, idebug_br)

            call wr_contact_pos(wtd, idebug_cnt, sub_ierror)
            if (my_ierror.eq.0) my_ierror = sub_ierror

            call broyden_get_state(wtd, neq, fk, idebug_br)

            dfk   = fk - fkm1

            ! check solution for fx

            if (any(isnan(dfk))) my_ierror = -1
            ! write(bufout,*) 'abs(dfk)=',abs(dfk),', dfk_max=',dfk_max,', test', (abs(dfk).lt.dfk_max)
            ! call write_log(1, bufout)

            laccept = ntry.ge.maxtry .or. (wtd%numcps.ge.1 .and. all(abs(dfk).lt.dfk_max))

            if (any(isnan(dfk))) then

               write(bufout,'(2(a,g12.4),a)') ' ...NaN-values found (', dfk(1), ',', dfk(2),            &
                        '), aborting'
               call write_log(1, bufout)

            elseif (.not.laccept .and. wtd%numcps.le.0) then

               call write_log(' ...update too large (loss of contact), rejecting step')
               write(bufout,'(2(a,g12.4))') ' Fy=',wtd%ftrk%y(),', Fz=',wtd%ftrk%z()
               call write_log(1, bufout)

            elseif (.not.laccept) then

               write(bufout,'(2(a,g12.4),a)') ' ...update too large (', dfk(1), ',', dfk(2),            &
                     '), rejecting step'
               call write_log(1, bufout)

            endif
         enddo

         ! check convergence

         res_fk = fk - ftarg
         ldone  = (all(abs(res_fk).lt.tol_fk) .or. k.ge.wtd%solv%maxnr)

         if (.not.ldone) then

            ! compute/update approximate Jacobian matrix
            ! Note: this overwrites the total forces, esp. wtd%ftrk

            if (use_findiff) then

               ! compute Jacobian using finite differences

               call broyden_findiff_jac(neq, wtd, xchr, xk, fk, bk, idebug_br, idebug_cnt, sub_ierror)
               if (my_ierror.eq.0) my_ierror = sub_ierror

            elseif (any(abs(dfk).gt.10d0*tol_fk)) then

               ! Broyden update: Bk = Bkm1 + (dfk - Bkm1 * dxk) * dxk^T / (dxk^T * dxk)

               call broyden_jac_update(neq, dxk, dfk, bkm1, bk)

            endif
         endif

         ! print output on Broyden process

         call broyden_print(neq, k, xk, fk, ftarg, bk, wtd%ic, ntry, idebug_br)

         if (my_ierror.lt.0) then
            write(bufout,'(a,i3,a)') ' An error occurred the in Broyden algorithm (',my_ierror,         &
                '), aborting contact solution'
            call write_log(1, bufout)
         endif

      enddo  ! while not done: Broyden loop

      ! N=1: subtract rail deflection from wheelset position z_ws

      if (wtd%ic%norm.eq.1) then
         wtd%ws%z = wtd%ws%z - wtd%trk%dz_defl
         wtd%trk%dz_defl = 0d0
      endif

      ! no convergence is treated here as -1 'no solution' (error) rather than -2 'inaccurate' (warning)

      if (my_ierror.eq.0 .and. any(abs(res_fk).ge.tol_fk)) my_ierror = -1

   end subroutine wr_contact_fyz_broyden

!------------------------------------------------------------------------------------------------------------

end module m_wr_totforce

