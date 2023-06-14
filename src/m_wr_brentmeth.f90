!------------------------------------------------------------------------------------------------------------
! m_wr_brentmeth - iterative solution of w/r contact (module 1) with total forces prescribed
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_wr_brentmeth

use m_wrprof_data
use m_hertz
use m_visc
use m_soutpt
use m_wr_solvecp
use m_wr_output
implicit none
private

public  wr_contact
private wr_contact_fx_secant
private secant_print

private brent_its_init
private brent_its_find_k
private brent_its_update_ptr
private brent_its_add_iterate
private brent_sensitivity_k
private brent_its_print
private brent_its_has_bracket
private brent_its_has_jump
private brent_its_has_chng_stiff
private wr_contact_init_hertz
private brent_set_xnew
private wr_contact_fxz_brent

private wr_contact_fxz_broyden
private broyden_set_state
private broyden_get_state
private broyden_get_target
private broyden_solve_dxk
private broyden_findiff_jac
private broyden_jac_update
private broyden_invjac
private broyden_print

integer, parameter :: IMETH_BRENT = 1, IMETH_BROYDN = 2, IMETH_SECANT = 3

!------------------------------------------------------------------------------------------------------------
!--data to keep track of iterates in Brent's algorithm
!       Brent's method keeps four iterates [a,b,c,d] with function values [fa,fb,fc,fd]
!       b is the best guess, c == b^{k-1}, d == b^{k-2}

   type :: t_brent_its
      integer                                 :: maxit, numit
      real(kind=8)                            :: ftarg
      integer,      dimension(:), allocatable :: k_all
      real(kind=8), dimension(:), allocatable :: xk_all, rk_all, dr_sec
      integer                                 :: k_0, k_1, k_a, k_b, k_c, k_d
      real(kind=8), pointer     :: x_a => NULL(), x_b => NULL(), x_c => NULL(), x_d => NULL(),          &
                                   r_a => NULL(), r_b => NULL(), r_c => NULL(), r_d => NULL()
      real(kind=8), pointer     :: x_0 => NULL(), r_0 => NULL(), x_1 => NULL(), r_1 => NULL()
      real(kind=8), pointer     :: x_k => NULL(), x_km1 => NULL(), x_km2 => NULL()

      ! maxit       size of arrays used, max. #iterates
      ! numit       actual #iterates stored
      ! ftarg       target force value
      ! xk_all      x-values used sorted in ascending order
      ! rk_all      residual values corresponding to x-values used
      ! dr_sec      derivatives dr/dx according to secant method
      ! k_all       iteration numbers corresponding to x-values used
      ! k_0, k_1    iteration numbers of current bracket points
      ! x_0, r_0    lower side of current bracket
      ! k_a -k_d    iteration numbers of current/previous best guesses
      ! x_1, r_1    upper side of current bracket
      ! x_a, r_a    contra-point to current best guess (other end of bracket)
      ! x_b, r_b    current best guess
      ! x_c, r_c    previous best guess
      ! x_d, r_d    previous previous best guess
      ! x_k, x_km1, x_km2   most recent iterates x^k, x^{k-1}, x^{k-2}
   end type t_brent_its

contains

!------------------------------------------------------------------------------------------------------------

   subroutine wr_contact(wtd, my_ierror)
!--purpose: process the whole w/r problem for a case, including NR-loops for prescribed forces
      implicit none
!--subroutine arguments:
      type(t_ws_track), target      :: wtd
      integer,          intent(out) :: my_ierror
!--local variables:
      logical, parameter :: use_brent = .true.
      logical            :: lfound, is_left
      integer            :: iestim_br, idebug_cnt, icp, i, sub_ierror, known_err(10)

      my_ierror = 0

      ! idebug: 0 = errors/silent, 1 = warnings/info, >=2 = flow/debug

      if (wtd%ic%nmdbg.le.0) then
         idebug_cnt = wtd%ic%ilvout
      else
         idebug_cnt = wtd%ic%nmdbg
      endif

      if (wtd%ic%norm.le.0) then

         ! solve w/r problem with positions & velocities prescribed

         if (idebug_cnt.ge.3) call write_log(' wr_contact: solving for given position...')
         wtd%meta%itbrent = 0
         call wr_contact_pos(wtd, idebug_cnt, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

      elseif (wtd%ic%force.le.0 .and. use_brent) then

         ! solve w/r problem with prescribed vertical force

         if (idebug_cnt.ge.3) call write_log(' wr_contact: solving for given vertical force...')
         iestim_br = 0
         call wr_contact_fxz_brent(wtd, iestim_br, idebug_cnt, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

      elseif (wtd%ic%force.le.0) then

         ! solve w/r problem with prescribed vertical and longitudinal forces

         if (idebug_cnt.ge.3) call write_log(' wr_contact: solving for given vert./long. forces (Broyden)...')
         call wr_contact_fxz_broyden(wtd, idebug_cnt, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

      else

         ! solve w/r problem with prescribed vertical and longitudinal forces

         if (idebug_cnt.ge.3) call write_log(' wr_contact: solving for given vert./long. forces (secant)...')
         call wr_contact_fx_secant(wtd, idebug_cnt, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

      endif

      ! translate errors

      if (my_ierror.eq.-1) then         ! -1 = NaN in residual/...
         my_ierror = CNTC_err_ftot
      elseif (my_ierror.eq.-2) then     ! -2 = residual > tolerance
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
      has_fx    = (wtd%ic%tang.ge.1 .and. wtd%ic%force.ge.1)
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
      call wr_contact_fxz_brent(wtd, iestim_br, idebug_cnt, sub_ierror)
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

      call secant_print(neq, k, xk, fk, ftarg, bk, wtd%ic, idebug_br)

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
                  write(bufout,'(2(a,f8.4))') ' ...halved step dom_ws=',dxk(2)
                  call write_log(1, bufout)
               endif
            endif

            xk   = xkm1 + dxk

            ! compute problem with new xk, get new fk

            wtd%ws%vpitch = xk(2)

            iestim_br = 1
            wtd%solv%eps = max(res_fk(2)/5d0, solv_eps)

            call wr_contact_fxz_brent(wtd, iestim_br, idebug_cnt, sub_ierror)
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

         call secant_print(neq, k, xk, fk, ftarg, bk, wtd%ic, idebug_br)

         if (my_ierror.lt.0) then
            write(bufout,'(a,i3,a)') ' An error occurred the in Secant algorithm (',my_ierror,          &
                '), aborting contact solution'
            call write_log(1, bufout)
         endif
      enddo  ! while not done: Secant loop

      if (my_ierror.eq.0 .and. any(abs(res_fk).ge.tol_fk)) my_ierror = -2

   end subroutine wr_contact_fx_secant

!------------------------------------------------------------------------------------------------------------

   subroutine secant_print(neq, k, xk, fk, ftarg, bk, ic, idebug_br)
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

   end subroutine secant_print

!------------------------------------------------------------------------------------------------------------

   subroutine brent_its_init(its, maxit, ftarg)
!--purpose: initialize the iterates structure for Brent's algorithm
      implicit none
!--subroutine arguments:
      type(t_brent_its), target :: its
      integer                   :: maxit
      real(kind=8)              :: ftarg

      ! call write_log(' brent_its_init...')
      its%maxit = max(1, min(10000, maxit))
      its%numit = 0
      its%ftarg = ftarg
      allocate(its%k_all(its%maxit))
      allocate(its%xk_all(its%maxit))
      allocate(its%rk_all(its%maxit))
      allocate(its%dr_sec(its%maxit))

      call brent_its_update_ptr(its)
      ! call write_log(' brent_its_init ok...')
   end subroutine brent_its_init

!------------------------------------------------------------------------------------------------------------

   function brent_its_find_k(its, k)
!--purpose: locate sequence number it for iteration k
      implicit none
!--function result:
      integer                   :: brent_its_find_k
!--subroutine arguments:
      type(t_brent_its), target :: its
      integer                   :: k
!--local variables:
      integer                   :: it

      brent_its_find_k = -1
      if (k.ge.0 .and. k.le.its%numit-1) then
         do it = 1, its%numit
            if (its%k_all(it).eq.k) brent_its_find_k = it
         enddo
      endif

      ! write(bufout,*) 'found k=',k,' at it=', brent_its_find_k
      ! call write_log(1, bufout)

   end function brent_its_find_k

!------------------------------------------------------------------------------------------------------------

   subroutine brent_its_update_ptr(its)
!--purpose: update the pointers x_a--x_d for current and previous best points
      implicit none
!--subroutine arguments:
      type(t_brent_its), target :: its
!--local variables
      integer      :: it, it_br0, imin

      ! call write_log(' brent_its_update_ptr...')

      if (its%numit.le.0) then

         its%k_0 = -1           ! current bracket
         its%k_1 = -1
         its%x_0 => NULL()
         its%r_0 => NULL()
         its%x_1 => NULL()
         its%r_1 => NULL()

         its%k_a = -1           ! points a, b, c, d
         its%k_b = -1
         its%k_c = -1
         its%k_d = -1
         its%x_a => NULL()
         its%x_b => NULL()
         its%x_c => NULL()
         its%x_d => NULL()
         its%r_a => NULL()
         its%r_b => NULL()
         its%r_c => NULL()
         its%r_d => NULL()

         its%x_k   => NULL()    ! most recent iterates
         its%x_km1 => NULL()
         its%x_km2 => NULL()

      elseif (its%numit.le.1) then

         its%k_0 =  its%k_all(1)
         its%k_1 =  its%k_all(1)
         its%x_0 => its%xk_all(1)
         its%x_1 => its%xk_all(1)
         its%r_0 => its%rk_all(1)
         its%r_1 => its%rk_all(1)

         its%k_a =  its%k_all(1)
         its%k_b =  its%k_all(1)
         its%k_c =  its%k_all(1)
         its%k_d =  its%k_all(1)
         its%x_a => its%xk_all(1)
         its%x_b => its%xk_all(1)
         its%x_c => its%xk_all(1)
         its%x_d => its%xk_all(1)
         its%r_a => its%rk_all(1)
         its%r_b => its%rk_all(1)
         its%r_c => its%rk_all(1)
         its%r_d => its%rk_all(1)

         its%x_k   => its%xk_all(1)
         its%x_km1 => its%xk_all(1)
         its%x_km2 => its%xk_all(1)
      else

         ! cycle previous best iterates

         its%k_d =  its%k_c
         its%k_c =  its%k_b

         it = brent_its_find_k(its, its%k_c)
         its%x_c => its%xk_all(it)
         its%r_c => its%rk_all(it)

         it = brent_its_find_k(its, its%k_d)
         its%x_d => its%xk_all(it)
         its%r_d => its%rk_all(it)

         ! determine current best iterate

         if (brent_its_has_bracket(its, it_br0)) then
            if (abs(its%rk_all(it_br0)).lt.abs(its%rk_all(it_br0+1))) then
               imin = it_br0
            else
               imin = it_br0 + 1
            endif
         else
            imin = idamin(its%numit, its%rk_all, 1)
         endif

         its%k_b =  its%k_all(imin)
         its%x_b => its%xk_all(imin)
         its%r_b => its%rk_all(imin)

         ! determine contra-point of current bracket (no bracket: sensible adjacent point)

         if (imin.le.1) then
            its%k_a =  its%k_all(imin+1)
            its%x_a => its%xk_all(imin+1)
            its%r_a => its%rk_all(imin+1)
         elseif (imin.ge.its%numit) then
            its%k_a =  its%k_all(imin-1)
            its%x_a => its%xk_all(imin-1)
            its%r_a => its%rk_all(imin-1)
         elseif (its%rk_all(imin-1)*its%r_b.lt.0d0) then
            its%k_a =  its%k_all(imin-1)
            its%x_a => its%xk_all(imin-1)
            its%r_a => its%rk_all(imin-1)
         else
            its%k_a =  its%k_all(imin+1)
            its%x_a => its%xk_all(imin+1)
            its%r_a => its%rk_all(imin+1)
         endif

         ! set pointers to bracket

         if (its%x_a.lt.its%x_b) then
            its%k_0 =  its%k_a
            its%k_1 =  its%k_b
            its%x_0 => its%xk_all(imin-1)
            its%x_1 => its%xk_all(imin)
            its%r_0 => its%rk_all(imin-1)
            its%r_1 => its%rk_all(imin)
         else
            its%k_0 =  its%k_b
            its%k_1 =  its%k_a
            its%x_0 => its%xk_all(imin)
            its%x_1 => its%xk_all(imin+1)
            its%r_0 => its%rk_all(imin)
            its%r_1 => its%rk_all(imin+1)
         endif

         ! set pointers to most recent iterates

         ! write(bufout,*) 'k=',its%numit,', it=', brent_its_find_k(its, its%numit)
         ! call write_log(1, bufout)

         it = brent_its_find_k(its,       its%numit-1)
         its%x_k   => its%xk_all(it)
         it = brent_its_find_k(its,       its%numit-2)
         its%x_km1 => its%xk_all(it)
         it = brent_its_find_k(its, max(0,its%numit-3))
         its%x_km2 => its%xk_all(it)

      endif
      ! call write_log(' brent_its_update_ptr ok...')

   end subroutine brent_its_update_ptr

!------------------------------------------------------------------------------------------------------------

   subroutine brent_its_add_iterate(its, k, xk, rk, drdx, my_ierror)
!--purpose: insert an iterate in the sorted structure for Brent's algorithm
      implicit none
!--subroutine arguments:
      type(t_brent_its) :: its
      integer           :: k, my_ierror
      real(kind=8)      :: xk, rk, drdx
!--local variables
      integer           :: i, j
      logical           :: ldone

      my_ierror = 0
      if (isnan(rk)) then
         my_ierror = -1
         write(bufout,'(a,f11.3,a)') ' ...NaN-values found (', rk, '), aborting'
         call write_log(1, bufout)
      endif

      ! determine first position i with xk < xk_all(i)

      ! call write_log(' brent_its_add_iterate...')

      i     = 0
      ldone = (i.ge.its%numit)
      do while (.not.ldone)
         i     = i + 1
         ldone = .true.
         if (i.le.its%numit) ldone = (xk.lt.its%xk_all(i))
      enddo
      ! write(bufout,*) 'i=',i,', num=',its%numit
      ! call write_log(1, bufout)

      ! shift iterates [i--end] one position

      do j = its%numit, max(i,1), -1
         ! write(bufout,*) 'j=',j
         ! call write_log(1, bufout)
         its%k_all(j+1)  = its%k_all(j)
         its%xk_all(j+1) = its%xk_all(j)
         its%rk_all(j+1) = its%rk_all(j)
         its%dr_sec(j+1) = its%dr_sec(j)
      enddo

      ! insert iterate at position i

      if (i.le.0) i = 1

      its%numit     = its%numit + 1
      its%k_all(i)  = k
      its%xk_all(i) = xk
      its%rk_all(i) = rk
      its%dr_sec(i) = drdx

      ! update pointers to current and previous best estimates

      call brent_its_update_ptr(its)
      ! call write_log(' brent_its_add_iterate ok...')

   end subroutine brent_its_add_iterate

!------------------------------------------------------------------------------------------------------------

   function brent_sensitivity_k(k, its, idebug_br)
!--purpose: estimate derivative df/dx at Brent iterate k
      implicit none
!--function result:
      real(kind=8)             :: brent_sensitivity_k
!--subroutine arguments:
      integer                  :: k, idebug_br
      type(t_brent_its)        :: its
!--local variables:
      integer                  :: it, j
      real(kind=8)             :: dfdx

      it = -1
      do j = 1, its%numit
         if (its%k_all(j).eq.k) it = j
      enddo

      if (it.eq.-1) then
         dfdx = -99d9
      elseif (its%numit.le.1) then
         dfdx = 0d0
      elseif (it.gt.1 .and. it.lt.its%numit) then
         ! central difference
         dfdx = (its%rk_all(it+1) - its%rk_all(it-1)) / max(tiny, its%xk_all(it+1) - its%xk_all(it-1))
      elseif (it.lt.its%numit) then
         ! forward difference
         dfdx = (its%rk_all(it+1) - its%rk_all(it)) / max(tiny, its%xk_all(it+1) - its%xk_all(it))
      else
         ! backward difference
         dfdx = (its%rk_all(it) - its%rk_all(it-1)) / max(tiny, its%xk_all(it) - its%xk_all(it-1))
      endif
      brent_sensitivity_k = dfdx

      if (idebug_br.ge.5) then
         write(bufout,'(2(a,i3),a,f12.4)') ' brent_sens: it=',it,', numit=',its%numit,', dfdx=',dfdx
         call write_log(1, bufout)
      endif

   end function brent_sensitivity_k

!------------------------------------------------------------------------------------------------------------

   subroutine brent_its_print(k, its, ic, idebug_br)
!--purpose: print information on Brent process
      implicit none
!--subroutine arguments:
      integer                  :: k, idebug_br
      type(t_ic)               :: ic
      type(t_brent_its)        :: its
!--local variables:
      integer                  :: it, j
      real(kind=8)             :: dfdx
      character(len=12)        :: str12(6)
      character(len=16)        :: str16(6)

      ! call write_log(' brent_its_print...')
      ! write(bufout,*) 'x_a:x_d point to k=',its%k_a, its%k_b, its%k_c, its%k_d
      ! call write_log(1, bufout)
      do it = 1, its%numit
         if (idebug_br.ge.3 .or. (ic%flow.ge.1 .and. its%k_all(it).eq.k)) then
            dfdx = brent_sensitivity_k(k, its, idebug_br)

            if (idebug_br.ge.2 .or. .false.) then
               str16(1) = fmt_gs(16,8, its%xk_all(it))
               str16(3) = fmt_gs(16,8, its%dr_sec(it))
               str16(2) = fmt_gs(16,8, its%rk_all(it))
               write(bufout,7010) its%k_all(it), (str16(j),j=1,3)
            else
               write(str12(1), '(f12.4)') its%xk_all(it)
               str12(2) = fmt_gs(12,4, its%rk_all(it)+its%ftarg)
               ! str12(3) = fmt_gs(12,4, dfdx)
               str12(3) = fmt_gs(12,4, its%dr_sec(it))
               write(bufout,7020) its%k_all(it), (str12(j),j=1,3)
            endif
            call write_log(1, bufout)
 7010       format(4x, i6,', NR,  z_ws, res: ',2a,', dFz/dz:',a)
 7020       format(4x, i6,', NR,  z_ws, Fz: ',2a,', dFz/dz:',a)
         endif
      enddo

      if (idebug_br.ge.3 .or. ic%flow.ge.6) then

         ! Process oriented: [xa,xb,xc], residues [fa,fb,fc]

         str12(1) = fmt_gs(12,4, its%x_a)
         str12(2) = fmt_gs(12,4, its%x_b)
         str12(3) = fmt_gs(12,4, its%x_c)
         str12(4) = fmt_gs(12,4, its%r_a)
         str12(5) = fmt_gs(12,4, its%r_b)
         str12(6) = fmt_gs(12,4, its%r_c)
         write(bufout,8006) k, (str12(j),j=1,6)
         call write_log(1, bufout)
 8006    format(4x, i6,', BR,  z_abc=[',a,',',a,',',a,'], res_abc=[',a,',',a,',',a,']')
      endif

   end subroutine brent_its_print

!------------------------------------------------------------------------------------------------------------

   function brent_its_has_bracket(its, it0)
!--purpose: determine whether solution is bracketed, return lower index of bracket
      implicit none
!--function result:
      logical           :: brent_its_has_bracket
!--subroutine arguments:
      type(t_brent_its)                       :: its
      integer,          intent(out), optional :: it0
!--local variables:
      integer           :: it

      ! call write_log(' brent_its_has_bracket...')
      if (its%numit.le.1) then
         brent_its_has_bracket = .false.
         it = 0
      else
         it = 1
         do while(it.lt.its%numit .and.  its%rk_all(it)*its%rk_all(it+1).gt.0d0)
            it = it + 1
         enddo
         brent_its_has_bracket = (it.lt.its%numit)
      endif
      if (present(it0)) it0 = it
      ! call write_log(' brent_its_has_bracket ok...')
   end function brent_its_has_bracket

!------------------------------------------------------------------------------------------------------------

   function brent_its_has_jump(its, idebug_br)
!--purpose: determine whether solution jumps across zero in the current bracket
      implicit none
!--function result:
      logical           :: brent_its_has_jump
!--subroutine arguments:
      type(t_brent_its) :: its
      integer           :: idebug_br
!--local variables:
      real(kind=8)      :: jump_thrs_df = 10d0, jump_thrs_dx = 0.01d0
      integer           :: it0, it1
      real(kind=8)      :: dx_brack, dfdx_brack, dfdx_left, dfdx_right

      if (its%numit.le.3 .or. .not.brent_its_has_bracket(its)) then
         brent_its_has_jump = .false.
      else
         it0 = brent_its_find_k(its, its%k_0)
         it1 = brent_its_find_k(its, its%k_1)

         dx_brack   = (its%x_1 - its%x_0)
         dfdx_brack = (its%r_1 - its%r_0) / dx_brack

         if (it0.gt.1) then
            ! backward difference
            dfdx_left = (its%r_0 - its%rk_all(it0-1)) / (its%x_0 - its%xk_all(it0-1))
         else
            dfdx_left = 1d20
         endif
         if (it1.lt.its%numit) then
            ! forward difference
            dfdx_right = (its%rk_all(it1+1) - its%r_1) / (its%xk_all(it1+1) - its%x_1)
         else
            dfdx_right = 1d20
         endif

         ! jump == dfdx in bracket >> dfdx on both sides

         brent_its_has_jump =     dx_brack .lt. jump_thrs_dx .and.                                      &
                              abs(dfdx_brack) .gt. jump_thrs_df*max(abs(dfdx_left), abs(dfdx_right))

         if (brent_its_has_jump .and. idebug_br.ge.-1) then
            write(bufout,'(2(a,f12.6),a)') '  ...detected jump in x=[', its%x_0,',', its%x_1,']'
            call write_log(1, bufout)
            write(bufout,'(3(a,g12.4))') '     df/dx=', dfdx_left,',', dfdx_brack,',', dfdx_right
            call write_log(1, bufout)
         endif
      endif

   end function brent_its_has_jump

!------------------------------------------------------------------------------------------------------------

   function brent_its_has_chng_stiff(its, idebug_br, x0_new)
!--purpose: determine whether solution has stiffness increasing across the current bracket
      implicit none
!--function result:
      logical           :: brent_its_has_chng_stiff
!--subroutine arguments:
      type(t_brent_its) :: its
      integer           :: idebug_br
      real(kind=8)      :: x0_new
!--local variables:
      real(kind=8)      :: stiff_thrs = 2d0, frac_last = 0.8d0
      integer           :: it0, it1
      real(kind=8)      :: dfdx_brack, dfdx_left, dfdx_right, x0_left, x0_right

      if (its%numit.le.2 .or. .not.brent_its_has_bracket(its)) then

         ! not enough points to compare stiffnesses

         brent_its_has_chng_stiff = .false.

         ! dummy output for x_new: bisection

         if (its%numit.le.1) then
            x0_new  = its%xk_all(1)
         else
            x0_new  = 0.5d0 * (its%xk_all(1) + its%xk_all(2))
         endif

      else

         ! determine left/right sides of current bracket

         it0 = brent_its_find_k(its, its%k_0)
         it1 = brent_its_find_k(its, its%k_1)

         if (it0.le.1) then

            ! bracket is first segment of list: report no stiffness change

            brent_its_has_chng_stiff = .false.

            ! dummy output for x_new: bisection

            x0_new    = 0.5d0 * (its%x_0 + its%x_1)

         elseif (it1.ge.its%numit) then

            ! bracket is last segment of list: compare left to current stiffness

            dfdx_left  = (its%r_0 - its%rk_all(it0-1)) / (its%x_0 - its%xk_all(it0-1))
            dfdx_brack = (its%r_1 - its%r_0)           / (its%x_1 - its%x_0)
            dfdx_right = 0d0

            brent_its_has_chng_stiff = (dfdx_brack .gt. stiff_thrs*dfdx_left)

            if (brent_its_has_chng_stiff) then
               x0_new = (1d0-frac_last) * its%x_0 + frac_last * its%x_1
            else
               x0_new = 0.5d0 * (its%x_0 + its%x_1)
            endif

         else

            ! bracket has left and right segments: compare left to right stiffness

            dfdx_left  = (its%r_0 - its%rk_all(it0-1)) / (its%x_0 - its%xk_all(it0-1))
            dfdx_brack = (its%r_1 - its%r_0)           / (its%x_1 - its%x_0)
            dfdx_right = (its%rk_all(it1+1) - its%r_1) / (its%xk_all(it1+1) - its%x_1)

            brent_its_has_chng_stiff = (dfdx_right .gt. stiff_thrs*dfdx_left)

            ! in case stiffness changes: use extrapolation from both sides

            if (brent_its_has_chng_stiff) then
               x0_left  = its%x_0 - its%r_0 / dfdx_left
               x0_right = its%x_1 - its%r_1 / dfdx_right

               x0_new   = min(x0_left, x0_right)
               x0_new   = max(x0_new, its%x_0+0.0001d0*(its%x_1-its%x_0))
            else
               x0_new   = 0.5d0 * (its%x_0 + its%x_1)
            endif

         endif

         if (brent_its_has_chng_stiff .and. idebug_br.ge.1) then
            write(bufout,'(2(a,f12.6),3(a,g14.6))') '  ...stiffness changing in x=[', its%x_0, ',',     &
                its%x_1,'], df/dx=', dfdx_left,',', dfdx_brack,',', dfdx_right
            call write_log(1, bufout)
         endif
      endif

   end function brent_its_has_chng_stiff

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
      wtd%meta%itbrent = k
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

         call brent_its_add_iterate(its, k, x_0, -wtd%ws%fz_inp, 0d0, sub_ierror)
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
         ! write(bufout,*) 'pen=',pen,', delt=',wtd%ws%delt_min,': est. dfz_dzws=',dfz_dzws
         ! call write_log(1, bufout)
         ! write(bufout,*) '1.5 * ',wtd%ws%fz_inp,' /',pen,' *',cos(wtd%ws%delt_min),' =', dfz_dzws
         ! call write_log(1, bufout)
      endif

      ! lower the wheel-set by gap_min and pen to get near the desired approach

      ! write(bufout,*) ' x_1=', x_0, ' +', pen, ' /', cos(wtd%ws%delt_min)
      ! call write_log(1, bufout)

      k     = 1
      wtd%meta%itbrent = k
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

         call brent_its_add_iterate(its, k, x_new, r_new, dfz_dzws, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, idebug_br)

      endif
         
   end subroutine wr_contact_init_hertz

!------------------------------------------------------------------------------------------------------------

   subroutine brent_set_xnew(k, its, dfz_dzws, used_bisec, tol_xk, x_new, idebug_br)
!--purpose: determine the next trial point using Brent's algorithm
      implicit none
!--subroutine arguments:
      integer,           intent(in)    :: k, idebug_br
      logical,           intent(inout) :: used_bisec
      real(kind=8),      intent(in)    :: dfz_dzws, tol_xk
      type(t_brent_its), target        :: its
      real(kind=8),      intent(out)   :: x_new
!--local variables:
      logical         :: use_inv_interp = .false.
      logical         :: ztest(5)
      real(kind=8)    :: dx, dr, dx_max, dx_new

      dr = its%rk_all(its%numit) - its%rk_all(1)

      if (.not.brent_its_has_bracket(its) .and. dr.lt.0d0 .and. k.le.5) then

         ! double the step

         x_new = its%xk_all(1) + 2d0 * (its%xk_all(its%numit) - its%xk_all(1))
         if (idebug_br.ge.2) call write_log(' negative force Fz, double step x(0)->x(k)')

      elseif (.not.brent_its_has_bracket(its) .and. k.le.3) then

         ! extrapolate using current estimate dfz_dzws

         x_new = its%x_b - its%r_b / dfz_dzws

         ! restrict maximum step: allow pen * 1.5

         dx_max = 0.5d0 * (its%x_b - its%xk_all(1))
         if (x_new.gt.its%x_b+dx_max .and. idebug_br.ge.3) then
            write(bufout,'(4(a,f9.4))') ' brent_set_xnew: x_0=', its%xk_all(1),', x_k=', its%x_b,       &
                ', x_new=',x_new,' >', its%x_b+dx_max
            call write_log(1, bufout)
         endif
         x_new  = min(x_new, its%x_b+dx_max)

      elseif (.not.brent_its_has_bracket(its)) then

         ! extrapolate using [x_0, r_0] and [x_n, r_n]

         dx = its%xk_all(its%numit) - its%xk_all(1)
         dr = its%rk_all(its%numit) - its%rk_all(1)

         ! write(bufout,*) 'x0,xn=',its%xk_all(1), its%xk_all(its%numit),', dx=', dx
         ! call write_log(1, bufout)
         ! write(bufout,*) 'r0,rn=',its%rk_all(1), its%rk_all(its%numit),', dr=', dr
         ! call write_log(1, bufout)

         ! restrict maximum step: in range [1.05, 1.2] * current interval

         dx_new = - its%rk_all(its%numit) * dx / dr
         dx_new = max(0.05d0*dx, min(0.2d0*dx, dx_new))
         x_new  = its%xk_all(its%numit) + dx_new

         ! write(bufout,*) 'dx_new=',dx_new,', x_new=',x_new
         ! call write_log(1, bufout)

      elseif (.false.) then

         ! stiffness increases over the current bracket: x_new set by chng_stiff
         ! code deleted in r.1932

      else ! has_bracket

         if (use_inv_interp .and. abs(its%r_a-its%r_c).gt.tiny .and. abs(its%r_b-its%r_c).gt.tiny) then

            ! attempt inverse quadratic interpolation

            x_new = its%x_a * its%r_b * its%r_c / ((its%r_a-its%r_b) * (its%r_a-its%r_c)) +             &
                    its%x_b * its%r_a * its%r_c / ((its%r_b-its%r_a) * (its%r_b-its%r_c)) +             &
                    its%x_c * its%r_a * its%r_b / ((its%r_c-its%r_a) * (its%r_c-its%r_b))
            if (idebug_br.ge.3) then
               write(bufout,'(a,2g12.4,a,f12.4)') ' dfk=', abs(its%r_a-its%r_c), abs(its%r_b-its%r_c),  &
                                      ', inverse interpolation: x_new =', x_new
               call write_log(1, bufout)
            endif

         elseif (k.le.2) then

            ! attempt the secant method using dfz_dzws from Hertz solution

            x_new = its%x_b - its%r_b / dfz_dzws
            if (idebug_br.ge.3) then
               write(bufout,'(a,g12.4,10x,a,f12.4)') ' dfz_dzws=', dfz_dzws,                            &
                                      ', secant method: x_new =', x_new
               call write_log(1, bufout)
            endif

         elseif (abs(its%r_b-its%r_c).gt.tiny) then

            ! attempt the secant method using current and previous iterates

            x_new = its%x_b - its%r_b * (its%x_b-its%x_c) / (its%r_b-its%r_c)
            if (idebug_br.ge.3) then
               write(bufout,'(a,2g12.4,a,f12.4)') ' dfk=', abs(its%r_a-its%r_c), abs(its%r_b-its%r_c),  &
                                      ', secant method: x_new =', x_new
               call write_log(1, bufout)
            endif
   
         else
   
            ! attempt linear interpolation (equals secant when c==a)
   
            x_new = its%x_b - its%r_b * (its%x_b-its%x_a) / (its%r_b-its%r_a)
            if (idebug_br.ge.3) then
               write(bufout,'(a,2g12.4,a,f12.4)') ' dfk=', abs(its%r_a-its%r_c), abs(its%r_b-its%r_c),  &
                                      ', linear interpolation: x_new =', x_new
               call write_log(1, bufout)
            endif
         endif

         ztest(1) = (x_new-(3d0*its%x_a+its%x_b)/4d0) * (x_new-its%x_b) .gt. 0d0
         ztest(2) =      used_bisec .and. abs(x_new-its%x_b) .ge. 0.5d0*abs(its%x_k  -its%x_km1)
         ztest(3) = .not.used_bisec .and. abs(x_new-its%x_b) .ge. 0.5d0*abs(its%x_km1-its%x_km2)
         ztest(4) =      used_bisec .and. abs(its%x_k  -its%x_km1) .lt. tol_xk
         ztest(5) = .not.used_bisec .and. abs(its%x_km1-its%x_km2) .lt. tol_xk

         if (idebug_br.ge.4) then
            if (ztest(1)) then
               write(bufout,*) ' 1: x_new=',x_new,' not in [', (3d0*its%x_a+its%x_b)/4d0, ',', its%x_b, ']'
               call write_log(1, bufout)
            elseif (ztest(2)) then
               ! call brent_its_print(k, its, ic, idebug_br+2)
               write(bufout,*) ' 2: x_new-its%x_b=',x_new-its%x_b,' >= ', 0.5d0*abs(its%x_k-its%x_km1)
               call write_log(1, bufout)
            elseif (ztest(3)) then
               write(bufout,*) ' 3: x_new-its%x_b=',x_new-its%x_b,' >= ', 0.5d0*abs(its%x_km1-its%x_km2)
               call write_log(1, bufout)
            elseif (ztest(4)) then
               write(bufout,*) ' 4: x_k  -its%x_km1=',its%x_k  -its%x_km1,' <= tol=', tol_xk
               call write_log(1, bufout)
            elseif (ztest(5)) then
               write(bufout,*) ' 5: x_km1-its%x_km2=',its%x_km1-its%x_km2,' <= tol=', tol_xk
               call write_log(1, bufout)
            endif
         endif

         ! reject interpolation when any condition fails, use bisection instead

         if ( ztest(1) .or. ztest(2) .or. ztest(3) .or. ztest(4) .or. ztest(5) ) then
            if (idebug_br.ge.3) call write_log('  ...reject estimate, using bisection')
            x_new = 0.5d0 * (its%x_a + its%x_b)
            used_bisec = .true.
         else
            used_bisec = .false.
         endif

      endif ! not has_bracket

   end subroutine brent_set_xnew

!------------------------------------------------------------------------------------------------------------

   subroutine wr_contact_fxz_brent(wtd, iestim_br, idebug_cnt, my_ierror)
!--purpose: process the w/r problem for a case with prescribed vertical force with Brents algorithm
      implicit none
!--subroutine arguments:
      type(t_ws_track), target      :: wtd
      integer,          intent(in)  :: iestim_br, idebug_cnt
      integer,          intent(out) :: my_ierror
!--local variables:
      integer, parameter       :: maxit = 100
      logical                  :: used_bisec, ldone
      integer                  :: k, idebug_br, sub_ierror
      real(kind=8)             :: ftarg, res_fk, tol_fk, dfk, tol_xk, x_0, x_new, r_new, drdx,           &
                                  dfz_dzws, aa, bb
      type(t_brent_its)        :: its

      my_ierror = 0
      idebug_br = 1

      if (idebug_cnt.ge.2) call write_log(' --- Start subroutine wr_contact_fxz_brent ---')

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
      call brent_its_init(its, maxit, ftarg)

      ! set zero point and initial estimate

      if (iestim_br.ge.1) then

         ! use initial estimate from ws%z, wtd%dfz_dzws; store initial point { x_0, r_0 }, print output

         k   = 0
         wtd%meta%itbrent = k
         call wr_contact_pos(wtd, idebug_cnt, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

         x_0 = wtd%ws%z_cnt0
         call brent_its_add_iterate(its, k, x_0, -wtd%ws%fz_inp, 0d0, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, idebug_br)

         k   = 1
         wtd%meta%itbrent = k
         x_new    = wtd%ws%z
         r_new    = wtd%ftrk%z() - wtd%ws%fz_inp
         dfz_dzws = wtd%dfz_dzws

         call brent_its_add_iterate(its, k, x_new, r_new, dfz_dzws, sub_ierror)
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
         call total_forces_moments(wtd, idebug_cnt)
         return
      endif

      ! Set tolerances on |fk| and |xk(2)-xk(1)|

      dfk     = abs(its%r_b - its%r_a) / max(1d-6, abs(its%x_b - its%x_a))
      tol_fk  = max(1d-12*dfk, wtd%solv%eps * abs(ftarg))
      tol_xk  = 1d-6

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
         wtd%meta%itbrent = k

         if (idebug_br.ge.3) then
            write(bufout,'(a,i3)')   ' wr_contact_fxz_brent: starting iteration k=',k
            call write_log(1, bufout)
         endif

         ! determine new point x_new

         call brent_set_xnew(k, its, dfz_dzws, used_bisec, tol_xk, x_new, idebug_br)

         ! compute problem with new xk, get new fk

         wtd%ws%z = x_new
         call wr_contact_pos(wtd, idebug_cnt, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

         r_new = wtd%ftrk%z() - ftarg
         if (abs(x_new - its%x_b).lt.tiny) then
            drdx = 1d20
         else
            drdx = (r_new - its%r_b) / (x_new - its%x_b)
         endif

         call brent_its_add_iterate(its, k, x_new, r_new, drdx, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, idebug_br)

         ! if (k.eq.15) call brent_its_print(k, its, wtd%ic, 2)

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

   end subroutine wr_contact_fxz_brent

!------------------------------------------------------------------------------------------------------------

   subroutine wr_contact_fxz_broyden(wtd, idebug_cnt, my_ierror)
!--purpose: process the w/r problem for a case with prescribed force(s) with Broyden algorithm,
!           using estimated Jacobian Bk
      implicit none
!--subroutine arguments:
      type(t_ws_track), target      :: wtd
      integer,          intent(in)  :: idebug_cnt
      integer,          intent(out) :: my_ierror
!--local variables:
      integer,       parameter :: neq = 2
      logical                  :: use_findiff_1st, use_findiff, has_fz, has_fx, is_roller, ldone, laccept
      integer                  :: k, idebug_br, itry, maxtry, sub_ierror
      real(kind=8)             :: aa, aob, bb, c11, c22, c23, dfx_domg,     &
                                  dfz_dzws, veloc, r_act
      real(kind=8)             :: xk(neq), xkm1(neq), xchr(neq), dxk(neq), fk(neq), fkm1(neq), dfk(neq), &
                                  ftarg(neq), res_fk(neq), tol_fk(neq), dfk_max(neq), bk(neq,neq),       &
                                  bkm1(neq,neq)
      type(t_brent_its)        :: its

      my_ierror = 0
      idebug_br = idebug_cnt
      idebug_br = 0
      if (idebug_br.ge.2) call write_log(' --- Start subroutine wr_contact_fxz ---')

      use_findiff_1st = .true.
      use_findiff     = .false.
      has_fz    = (wtd%ic%norm.ge.1)
      has_fx    = (wtd%ic%tang.ge.1 .and. wtd%ic%force.ge.1)
      is_roller = (wtd%ic%config.ge.4)

      ! preparations for initial estimate for z-position

      if (has_fz) then

         ! analyze geometry for overall minimum gap and curvatures, set initial estimates k=0, k=1

         call wr_contact_init_hertz(wtd, its, dfz_dzws, aa, bb, IMETH_BROYDN, idebug_br, idebug_cnt,    &
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

      endif

      ! Broyden: solve equations Fz(z_ws, omg_ws) / Fzinp       = 1
      !                          Fx(z_ws, omg_ws) / Fz(zws,omg) = Fxinp/Fzinp
      ! vector xk == [ z_ws; omg_ws ]; fk == [ Fz/Fzinp ; Fx/Fz ]

      ! target values for equations f(x) = ftarg

      call broyden_get_target(wtd, neq, ftarg, tol_fk, dfk_max, dfz_dzws, idebug_br)

      ! initialize the Broyden loop for obtaining the prescribed vertical force

      ! initial estimate [ z_ws, omg_ws ]

      k  = 0
      xk = (/ wtd%ws%z , wtd%ws%vpitch /)

      ! solve contact problem for the initial estimate

      call broyden_set_state(neq, xk, wtd, idebug_br)

      call wr_contact_pos(wtd, idebug_cnt, sub_ierror)
      if (my_ierror.eq.0) my_ierror = sub_ierror

      call broyden_get_state(wtd, neq, fk, idebug_br)

      ! set approximate Jacobian matrix B_0

      bk      =  0d0
      bk(1,1) =  -1.5d0 * fk(1) / wtd%ws%gap_min
      bk(2,2) =  dfx_domg / wtd%ws%fz_inp
      if (.not.has_fx) bk(2,2) = 1d0

      ! set characteristic sizes of x-values: perturbation = xchr / 1000.

      xchr(1) = ftarg(1) / bk(1,1)
      xchr(2) = 0.1d0 * veloc / wtd%ws%nom_radius

      ! compute Jacobian using finite differences

      if (has_fx .and. use_findiff_1st) then

         call broyden_findiff_jac(neq, wtd, xchr, xk, fk, bk, idebug_br, idebug_cnt, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

      endif

      ! check convergence

      res_fk = fk - ftarg

      ldone = (all(abs(res_fk).lt.tol_fk) .or. k.ge.wtd%solv%maxnr .or. wtd%ic%return.ge.2)

      if (idebug_br.ge.2) then
         write(bufout, '(4(a,f11.3))') ' fz_inp=',ftarg(1),', fk=',fk(1),', res=',res_fk(1),            &
                                                                                ', tol=',tol_fk(1)
         call write_log(1, bufout)
      endif

      ! print output on Broyden process

      call broyden_print(neq, k, xk, fk, ftarg, bk, wtd%ic, idebug_br)

      ! while not "done" do

      do while (.not.ldone .and. my_ierror.eq.0)

         ! increment iteration number

         k = k + 1

         if (k.eq.28) idebug_br = 3

         if (idebug_br.ge.2) then
            write(bufout,'(a,i3)')   ' wr_contact_fxz: starting iteration k=',k
            call write_log(1, bufout)
         endif

         ! if (idebug_br.ge.2) stop

         ! cycle previous values

         xkm1 = xk
         fkm1 = fk
         bkm1 = bk

         ! compute increment dxk and new iterand xk

         call broyden_solve_dxk(neq, bkm1, res_fk, dxk, idebug_br, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

         ! no change on omega_ws in first (few) iterations

         if (has_fx .and. k.le.1) dxk(2) = 0d0

         ! reduce step-size to avoid loss of contact

         do while(xkm1(1)+dxk(1).lt.wtd%ws%z_cnt0)
            if (idebug_br.ge.-1) then
               write(bufout,'(3(a,f8.4))') ' ...z_ws=',xkm1(1)+dxk(1),', z_ws0=',wtd%ws%z_cnt0,         &
                        ', halved step dz_ws=',0.5d0*dxk(1)
               call write_log(1, bufout)
            endif
            dxk = 0.5d0 * dxk
         enddo

         ! Use back-tracking for dz, domega_ws when dF becomes too large

         itry =  0
         maxtry = 10
         laccept = .false.

         do while(.not.laccept .and. my_ierror.eq.0)

            itry = itry + 1

            ! reduce step-size dx on consecutive inner iterations

            if (itry.gt.1) then
               dxk = 0.5d0 * dxk
               write(bufout,'(2(a,f8.4))') ' ...halved step dz_ws=',dxk(1),', dom_ws=',dxk(2)
               call write_log(1, bufout)
            endif

            xk   = xkm1 + dxk

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

            laccept = (itry.ge.maxtry .or. all(abs(dfk).lt.dfk_max))

            if (any(isnan(dfk))) then
               write(bufout,'(2(a,f11.3),a)') ' ...NaN-values found (', dfk(1), ',', dfk(2),            &
                        '), aborting'
               call write_log(1, bufout)
            elseif (.not.laccept) then
               write(bufout,'(2(a,f11.3),a)') ' ...update too large (', dfk(1), ',', dfk(2),         &
                     '), rejecting step'
               call write_log(1, bufout)
            endif
         enddo

         ! compute/update approximate Jacobian matrix

         if (has_fx .and. use_findiff) then

            ! compute Jacobian using finite differences

            call broyden_findiff_jac(neq, wtd, xchr, xk, fk, bk, idebug_br, idebug_cnt, sub_ierror)
            if (my_ierror.eq.0) my_ierror = sub_ierror

         elseif (any(abs(dfk).gt.10d0*tol_fk)) then

            ! Broyden update: Bk = Bkm1 + (dfk - Bkm1 * dxk) * dxk^T / (dxk^T * dxk)

            call broyden_jac_update(neq, dxk, dfk, bkm1, bk, idebug_br)

         endif

         ! check convergence

         res_fk = fk - ftarg
         ldone = (all(abs(res_fk).lt.tol_fk) .or. k.ge.wtd%solv%maxnr)

         ! print output on Broyden process

         call broyden_print(neq, k, xk, fk, ftarg, bk, wtd%ic, idebug_br)

         if (my_ierror.lt.0) then
            write(bufout,'(a,i3,a)') ' An error occurred the in Broyden algorithm (',my_ierror,         &
                '), aborting contact solution'
            call write_log(1, bufout)
         endif
      enddo  ! while not done: Broyden loop

      if (my_ierror.eq.0 .and. any(abs(res_fk).ge.tol_fk)) my_ierror = -2

   end subroutine wr_contact_fxz_broyden

!------------------------------------------------------------------------------------------------------------

   subroutine broyden_set_state(neq, xk, wtd, idebug_br)
!--purpose: set values in wtd corresponding to Broyden vector xk 
      implicit none
!--subroutine arguments:
      integer,          intent(in)  :: neq, idebug_br
      real(kind=8),     intent(in)  :: xk(neq)
      type(t_ws_track), target      :: wtd
!--local variables:

      wtd%ws%z      = xk(1)
      wtd%ws%vpitch = xk(2)

      if (idebug_br.ge.2) then
         write(bufout, '(a,g14.6,a,g15.8)') ' Broyden: set state z_ws=', xk(1),', vpitch=', xk(2)
         call write_log(1, bufout)
      endif

   end subroutine broyden_set_state

!------------------------------------------------------------------------------------------------------------

   subroutine broyden_get_state(wtd, neq, fk, idebug_br)
!--purpose: get values from wtd corresponding to Broyden vector fk 
      implicit none
!--subroutine arguments:
      integer,          intent(in)  :: neq, idebug_br
      real(kind=8),     intent(out) :: fk(neq)
      type(t_ws_track), target      :: wtd
!--local variables:
      logical           :: has_fx

      has_fx = (wtd%ic%tang.ge.1 .and. wtd%ic%force.ge.1)

      fk(1) = wtd%ftrk%z() / wtd%ws%fz_inp
      fk(2) = wtd%fws%x() / wtd%ftrk%z()       ! fx_rel @ current fz
      if (.not.has_fx) fk(2) = wtd%ws%vpitch

      if (idebug_br.ge.2) then
         write(bufout, '(2(a,f12.6,a,f12.3),a)') ' Broyden: get state Fz=', fk(1),' (',wtd%ftrk%z(),    &
                '), Fx=', fk(2),' (',wtd%fws%x(),')'
         call write_log(1, bufout)
      endif

   end subroutine broyden_get_state

!------------------------------------------------------------------------------------------------------------

   subroutine broyden_get_target(wtd, neq, ftarg, tol_fk, dfk_max, dfz_dzws, idebug_br)
!--purpose: get target values from wtd corresponding to Broyden vector ftarg 
      implicit none
!--subroutine arguments:
      type(t_ws_track), target      :: wtd
      integer,          intent(in)  :: neq, idebug_br
      real(kind=8),     intent(in)  :: dfz_dzws
      real(kind=8),     intent(out) :: ftarg(neq), tol_fk(neq), dfk_max(neq)
!--local variables:
      logical           :: has_fx

      has_fx = (wtd%ic%tang.ge.1 .and. wtd%ic%force.ge.1)

      ! target forces

      ftarg(1) = wtd%ws%fz_inp / wtd%ws%fz_inp
      ftarg(2) = wtd%ws%fx_inp / wtd%ws%fz_inp
      if (.not.has_fx) ftarg(2) = wtd%ws%vpitch

      ! tolerances. max. resolution for Fz: z_ws \pm 1e-14, Fz \pm 1e-14 * dfz_dzws

      tol_fk(1) = max(1d-12*dfz_dzws/wtd%ws%fz_inp, wtd%solv%eps)
      tol_fk(2) = wtd%solv%eps * wtd%fric%fstat_min()
      if (.not.has_fx) tol_fk(2) = wtd%solv%eps

      ! maximum change in forces per iteration

      dfk_max(1) =  99d0
      dfk_max(2) = 1.0d0 * wtd%fric%fstat_min()

      if (idebug_br.ge.2) then
         write(bufout,'(2(a,2g12.4))') ' tol_fk=',tol_fk(1), tol_fk(2),', dfk_max=',dfk_max(1), dfk_max(2)
         call write_log(1, bufout)
      endif

   end subroutine broyden_get_target

!------------------------------------------------------------------------------------------------------------

   subroutine broyden_solve_dxk(neq, bmat, res_fk, dxk, idebug_br, my_ierror)
!--purpose: solve update: B_k-1 * dxk = dfk
      implicit none
!--subroutine arguments:
      integer,     intent(in)  :: neq, idebug_br
      integer,     intent(out) :: my_ierror
      real(kind=8)             :: bmat(neq,neq), res_fk(neq), dxk(neq)
!--local variables:
      integer                    :: j
      real(kind=8)               :: det, binv(neq,neq)

      my_ierror = 0

      ! compute inverse of 2x2 matrix

      det = bmat(1,1) * bmat(2,2) - bmat(1,2) * bmat(2,1)

      ! abort if zero determinant is found

      if (abs(det).le.1e-10) then
         if (idebug_br.ge.1) call write_log(' Internal error in Broyden algorithm: zero determinant')
         my_ierror = -3
         return
      endif

      binv(1,1) =  bmat(2,2) / det
      binv(1,2) = -bmat(1,2) / det
      binv(2,1) = -bmat(2,1) / det
      binv(2,2) =  bmat(1,1) / det

      dxk  = 0d0
      do j = 1, neq
         dxk  = dxk - res_fk(j) * binv(:,j)
      enddo

      if (idebug_br.ge.2) then
         write(bufout, '(2(a,g14.6),a)')        ' Broyden: res_fk = [', res_fk(1), ',', res_fk(2), ']^T'
         call write_log(1, bufout)
         write(bufout, '(2(2(a,g14.6),a,:,/))') ' Broyden: B_inv  = [', binv(1,1), ',', binv(1,2), ']', &
                                                '                   [', binv(2,1), ',', binv(2,2), '].'
         call write_log(2, bufout)
         write(bufout, '(2(a,g14.6),a)')        ' Broyden: d_xk   = [', dxk(1),    ',', dxk(2),    ']^T'
         call write_log(1, bufout)
      endif

   end subroutine broyden_solve_dxk

!------------------------------------------------------------------------------------------------------------

   subroutine broyden_findiff_jac(neq, wtd, xchr, xk, fk, bk, idebug_br, idebug_cnt, my_ierror)
!--purpose: compute approximate Jacobian for Broyden process using finite differences
      implicit none
!--subroutine arguments:
      integer,          intent(in)  :: neq, idebug_br, idebug_cnt
      integer,          intent(out) :: my_ierror
      type(t_ws_track), target      :: wtd
      real(kind=8)                  :: xchr(neq), xk(neq), fk(neq), bk(neq,neq)
!--local variables:
      integer           :: j, sub_ierror
      real(kind=8)      :: xtmp(neq), ftmp(neq)

      my_ierror = 0

      do j = 1, neq

         ! solve contact problem with one input perturbed

         xtmp = xk
         xtmp(j) = xtmp(j) + 0.001d0 * xchr(j)

         call broyden_set_state(neq, xtmp, wtd, idebug_br)

         call wr_contact_pos(wtd, idebug_cnt, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

         call broyden_get_state(wtd, neq, ftmp, idebug_br)

         bk(1,j) = (ftmp(1) - fk(1)) / (xtmp(j) - xk(j))
         bk(2,j) = (ftmp(2) - fk(2)) / (xtmp(j) - xk(j))

      enddo

   end subroutine broyden_findiff_jac

!------------------------------------------------------------------------------------------------------------

   subroutine broyden_jac_update(neq, dxk, dfk, bkm1, bk, idebug_br)
!--purpose: update approximate Jacobian for Broyden process
      implicit none
!--subroutine arguments:
      integer,     intent(in)  :: neq, idebug_br
      real(kind=8)             :: dxk(neq), dfk(neq), bkm1(neq,neq), bk(neq,neq)
!--local variables:
      integer                  :: j
      real(kind=8)             :: fac, vec1(neq)

      ! Broyden: Bk = Bkm1 + (dfk - Bkm1 * dxk) * dxk^T / (dxk^T * dxk)

      ! numerator: column-vector v1 = dfk - Bkm1 * dxk

      vec1 = dfk
      do j = 1, neq
         vec1 = vec1 - dxk(j) * bkm1(:,j)
      enddo

      ! denominator: scalar fac = dxk .dot. dxk

      fac = sum( dxk(:) * dxk(:) )

      ! update approximate Jacobian: bk = bkm1 + column vec1 * row dxk / fac

      bk  = bkm1
      do j = 1, neq
         bk(:,j) = bk(:,j) + vec1(:) * dxk(j) / fac
      enddo

      ! check sign on dfz/dpen, reject bk when negative

      if (bk(1,1).lt.0d0) then
         if (idebug_br.ge.1) call write_log('   ...negative derivative dfz/dpen, rejecting Jacobian')
         bk = bkm1
      endif

   end subroutine broyden_jac_update

!------------------------------------------------------------------------------------------------------------

   subroutine broyden_invjac(neq, dxk, dfk, ckm1, ck)
!--purpose: update inverse Jacobian for Broyden process
      implicit none
!--subroutine arguments:
      integer                  :: neq
      real(kind=8)             :: dxk(neq), dfk(neq), ckm1(neq,neq), ck(neq,neq)
!--local variables:
      integer                  :: i, j
      real(kind=8)             :: fac, vec1(neq), vec2(neq)

      ! Broyden: Ck = Ckm1 + (dxk - Ckm1 * dfk) / (dxk^T * Ckm1 * dfk) * dxk^T * Ckm1

      ! numerator: column-vector v1 = dxk - Ckm1 * dfk

      vec1 = dxk
      do j = 1, neq
         vec1 = vec1 - dfk(j) * ckm1(:,j)
      enddo

      ! final part: row-vector v2 = dxk^T * Ckm1

      vec2 = 0d0
      do i = 1, neq
         vec2(i) = sum( dxk(:) * ckm1(:,i) )
      enddo

      ! denominator: scalar fac = (dxk^T * Ckm1 * dfk) = v2 .dot. dfk

      fac = sum( vec2(:) * dfk(:) )

      ! update inverse Jacobian: ck = ckm1 + column vec1 * row vec2 / fac

      ck  = ckm1
      do j = 1, neq
         ck(:,j) = ck(:,j) + vec1(:) * vec2(j) / fac
      enddo

   end subroutine broyden_invjac

!------------------------------------------------------------------------------------------------------------

   subroutine broyden_print(neq, k, xk, fk, ftarg, bk, ic, idebug_br)
!--purpose: print information on Broyden process
      implicit none
!--subroutine arguments:
      type(t_ic)               :: ic
      integer                  :: neq, k, idebug_br
      real(kind=8)             :: xk(neq), fk(neq), ftarg(neq), bk(neq,neq)
!--local variables:
      logical                  :: has_fz, has_fx
      integer                  :: j
      character(len=12)        :: str12((2+neq)*neq)
      character(len=16)        :: str16((2+neq)*neq)

      has_fz = (ic%norm.ge.1)
      has_fx = (ic%tang.ge.1 .and. ic%force.ge.1)

      if (idebug_br.ge.4) then
         write(bufout,*) 'Bk = [', bk(1,1),',', bk(1,2),']'
         call write_log(1, bufout)
         write(bufout,*) '     [', bk(2,1),',', bk(2,2),']'
         call write_log(1, bufout)
      endif

      if (.not.has_fx) then

         ! vertical force prescribed: F_z = F_z( z_ws )

         if (idebug_br.ge.2 .or. .false.) then
            str16(1) = fmt_gs(16,10, xk(1))
            str16(2) = fmt_gs(16,10, fk(1))
            str16(3) = fmt_gs(12,4, bk(1,1))
            write(bufout,7000) k, (str16(j),j=1,3)
            call write_log(1, bufout)
         elseif (idebug_br.ge.2 .or. ic%flow.ge.1) then
            str12(1) = fmt_gs(12,4, xk(1))
            str12(2) = fmt_gs(12,4, fk(1))
            str12(3) = fmt_gs(12,4, bk(1,1))
            write(bufout,7000) k, (str12(j),j=1,3)
            call write_log(1, bufout)
         endif

      else

         ! horizontal and vertical forces prescribed: {F_z,F_x} = F( z_ws, omg_ws )

         if (idebug_br.ge.2 .or. ic%flow.ge.1) then
            str16(1) = fmt_gs(16,8, xk(1))
            write(str16(2), '(f16.8)') fk(1) - ftarg(1)
            write(str16(3), '(f12.6)') bk(1,1)
            write(str16(4), '(f12.6)') bk(1,2)

            str16(5) = fmt_gs(16,8, xk(2))
            write(str16(6), '(f16.8)') fk(2) - ftarg(2)
            write(str16(7), '(f12.6)') bk(2,1)
            write(str16(8), '(f12.6)') bk(2,2)
            write(bufout,8000) k, (str16(j),j=1,8)
            call write_log(2, bufout)
         endif

      endif

 7000 format(4x, i6,', NR,  z_ws, Fz: ',2a,', dFz/dz:',a)
 8000 format(4x, i6,', NR,  z_ws, Fz: ',2a,', dFz: ',2a, /,                                          &
             4x, 6x, 4x,'  om_ws, Fx: ',2a,', dFx: ',2a)

   end subroutine broyden_print

!------------------------------------------------------------------------------------------------------------

end module m_wr_brentmeth

