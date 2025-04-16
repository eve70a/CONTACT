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

      if (wtd%ic%norm.le.0 .and. wtd%ic%force1.le.0) then

         ! N=0, F=0: solve w/r problem with positions & velocities prescribed

         if (x_locate.ge.3) call write_log(' wr_contact: solving for given position...')
         call wr_contact_pos(wtd, x_locate, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

      elseif (wtd%ic%force1.le.0) then

         ! N=1, F=0: solve w/r problem with prescribed vertical force

         if (x_locate.ge.3) call write_log(' wr_contact: solving for given vertical force...')

         iestim_br = 0
         if (wtd%ic%iestim.ge.1) iestim_br = 1

         call wr_contact_fz_brent(wtd, iestim_br, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

      elseif (wtd%ic%force1.eq.1 .or. wtd%ic%force1.eq.2) then

         ! N=1, F=1,2: solve w/r problem with prescribed vertical and longitudinal forces

         if (x_locate.ge.3) call write_log(' wr_contact: solving for given vert./long. forces (secant)...')
         call wr_contact_fx_secant(wtd, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

      elseif (wtd%ic%force1.eq.3) then

         ! N=0/1, F=3: solve w/r problem with prescribed vertical and lateral forces - Brent method

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
            call writmt (gd%meta, gd%ic, gd%cgrid_cur, gd%potcon_cur, gd%mater, gd%fric, gd%kin,        &
                gd%geom, gd%outpt1, wtd%ic%is_left_side())
            end associate
         enddo
      endif

      ! compute subsurface stresses when requested

      if (wtd%ic%stress.ge.1) then
         call wr_subsurf(wtd, x_locate)
      endif

      ! write output when R=0 or 1

      if (wtd%ic%return.le.1) call wr_output(wtd)

      ! destroy remaining (unconnected) contact patches (gds) of the previous time

      if (x_locate.ge.1) then
         write(bufout,'(3(a,i3),a)') ' There are',wtd%numcps,' contact patches,',wtd%n_miss,            &
                ' near misses, and', wtd%numtot-wtd%n_miss-wtd%numcps,' patches from previous time'
         call write_log(1, bufout)
      endif

      do icpo = wtd%numcps+wtd%n_miss+1, wtd%numtot
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
      wtd%numtot = wtd%numcps + wtd%n_miss

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

         if (wtd%ic%x_cpatch.ge.1) call write_log(' shift_grid_super: adding supergrid offset')

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
               xref_pot  = cref_x          !                    moving longitudinally at x_tr==0
             ! xref_pot  = 0d0             ! HACK for old behavior!!!
               yref_pot  = sr_ref
            else                           ! steady rolling: contact grid centered at contact reference marker
               xref_pot  = 0d0
               yref_pot  = 0d0
            endif

            ! snap contact reference to nearest supergrid location

            if (wtd%ic%x_cpatch.ge.1) then
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

   subroutine wr_contact_init_hertz(wtd, its, ftarg, dftot_dz, aa, bb, imeth, x_force, x_locate, my_ierror)
!--purpose: determine an initial bracket for vertical position for starting the Brent algorithm
      implicit none
!--subroutine arguments:
      type(t_ws_track)                 :: wtd
      real(kind=8),      intent(in)    :: ftarg
      integer,           intent(in)    :: imeth
      integer,           intent(in)    :: x_force, x_locate
      integer,           intent(out)   :: my_ierror
      real(kind=8),      intent(out)   :: dftot_dz, aa, bb
      type(t_brent_its), intent(inout) :: its
!--local variables:
      integer          :: k, ic_norm, ipotcn, ic_return, sub_ierror
      real(kind=8)     :: z_shift, z_shift0, fz_cntc, ftot, dftot_dpen, r_new,                          &
                          cp, rho, pen, e_star, epshz, fz_hz, azz

      if (x_force.ge.3) call write_log(' --- Start subroutine wr_contact_init_hertz ---')

      my_ierror = 0

      ! initial estimate for z-position: z_shift = z_ws - dz_defl = 0 

      k       = 0
      z_shift = 0d0
      call set_z_shift(wtd, z_shift)

      ! analyze geometrical problem for the initial z-position without actual solving

      ic_return = wtd%ic%return
      wtd%ic%return = 3
      wtd%meta%itforc_inn = k

      call wr_contact_pos(wtd, x_locate, sub_ierror)
      if (my_ierror.eq.0) my_ierror = sub_ierror

      wtd%ic%return = ic_return

      if (x_force.ge.3) then
         write(bufout,'(3(a,g14.6))') ' z_shift =',z_shift,', z_ws =',wtd%ws%z,', dz_defl =',wtd%trk%dz_defl
         call write_log(1, bufout)
         write(bufout,'(3(a,g14.6))') ' Overall minimum gap=', wtd%ws%gap_min,', a1=',wtd%ws%a1,        &
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
         
      ! store initial point { z_shift0, res0 }, print output
      ! TODO: check r_new < 0, support loss of contact with massless rail model

      z_shift0 = wtd%ws%gap_min
      call set_z_shift(wtd, z_shift0)
      fz_cntc  = 0d0
      ftot     = get_f_totz(wtd, fz_cntc)
      r_new    = ftot - ftarg

      ! write(bufout,'(4(a,g12.4))') ' z_shift0=',z_shift0,', ftot=',ftot,', ftarg=',ftarg,', res=',r_new
      ! call write_log(1, bufout)

      if (imeth.eq.IMETH_BRENT) then

         call brent_its_add_iterate(its, k, wtd%numcps, z_shift0, r_new, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, x_force)

      endif

      ! Use the Hertzian approach to get a sensible initial z-position

      ! solve the Hertzian problem for estimated curvatures, approximating fn ~=~ fz
      ! Note: for cylinder on plane configurations it may useful to switch to 2D Hertz

      ic_norm =  1
      ipotcn  = -1
      e_star  = wtd%mater%ga / (1d0 - wtd%mater%nu)
      epshz   = wtd%solv%eps
         
      fz_hz   = -r_new  ! Note: does not account for kz * dz_defl
      call hzcalc3d (e_star, epshz, ipotcn, wtd%ws%a1, wtd%ws%b1, aa, bb, ic_norm, pen, fz_hz, cp, rho)
      ! write(bufout,'(2(a,g12.4),a,i0,2(a,g12.4))') ' e_star=',e_star,', eps=',epshz,', ipotcn=',     &
      !        ipotcn,', a1=',wtd%ws%a1, ', b1=',wtd%ws%b1
      ! call write_log(1, bufout)
      ! write(bufout,'(2(a,f8.3),a,i0,4(a,g12.4))') ' aa=',aa,', bb=',bb,', ic_norm=',ic_norm,', pen=',&
      !        pen,', fz_hz=',fz_hz, ', cp=',cp,', rho=',rho
      ! call write_log(1, bufout)

      ! estimate derivative dF_tot/dz_shift ~=~ dFn/dpen

      if (pen.le.1d-6) then
         ! no contact: using influence coefficient Azz(0,0) for 1 element in contact

         call azz_one_element(wtd%ic, wtd%mater, wtd%discr%dx, wtd%discr%ds, azz)
         dftot_dpen = wtd%discr%dx * wtd%discr%ds * wtd%mater%ga / azz 
         write(bufout,*) 'gap_min=',wtd%ws%gap_min,', pen=',pen,', using dftot_dpen=',dftot_dpen
         call write_log(1, bufout)
      else
         ! when F=0,N=1, Fz_hz == wtd%ws%fz_inp
         dftot_dpen = 1.5d0 * fz_hz / pen

         if (x_force.ge.3) then
            write(bufout,'(3(a,g14.6))') ' fz_hz=',fz_hz,', pen=',pen,': est. dftot_dpen=',dftot_dpen
            call write_log(1, bufout)
         endif
      endif

      dftot_dz = dftot_dpen * cos(wtd%ws%delt_min) - wtd%trk%kz_rail

      ! shift the wheel or rail by gap_min and pen to get near the desired approach

      k     = 1
      wtd%meta%itforc_inn = k
      z_shift = z_shift0 + pen / max(0.5d0, cos(wtd%ws%delt_min))
      call set_z_shift(wtd, z_shift)

      ! Brent: solve contact problem for the current estimate

      if (imeth.eq.IMETH_BRENT) then

         call wr_contact_pos(wtd, x_locate, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

         ftot   = get_f_totz(wtd)
         r_new  = ftot - ftarg
         dftot_dz = 1.5d0 * wtd%ftrk%z() / (z_shift - z_shift0)

         ! write(bufout,'(4(a,g12.4))') ' x_1=',x_new,', ftot=',ftot,', ftarg=',ftarg,', res=',r_new
         ! call write_log(1, bufout)

         ! write(bufout,*) 'dftot_dz= =', dftot_dz
         ! call write_log(1, bufout)

         call brent_its_add_iterate(its, k, wtd%numcps, z_shift, r_new, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, x_force)

      endif
         
   end subroutine wr_contact_init_hertz

!------------------------------------------------------------------------------------------------------------

   subroutine solve_zshift_powerlaw(wtd, k, its, ftarg, tol_zk, zshift_new, dftot_dz, x_force,          &
                                 was_complete, use_pow)
!--purpose: determine new z_shift using approximate power law for each contact patch: F_n = C_i pen_i^1.5 
      implicit none
!--subroutine arguments:
      type(t_ws_track)                 :: wtd
      integer,           intent(in)    :: k, x_force
      type(t_brent_its), intent(in)    :: its
      real(kind=8),      intent(in)    :: ftarg, tol_zk
      real(kind=8),      intent(out)   :: zshift_new, dftot_dz
      logical,           intent(inout) :: was_complete
      logical,           intent(out)   :: use_pow
!--local variables:
      integer,      parameter :: max_nr = 20
      real(kind=8), parameter :: eps_fz = 1d-3
      logical        :: is_complete, new_c_hz(MAX_NUM_CPS)
      integer        :: icp, it_nr, ic_norm, ipotcn
      real(kind=8)   :: z_shift, z_shift0, dzshift, zshift_prv, dz_defl, pen, res, sgn,                 &
                        alph, c_hz, fs_rel, fz_est, e_star, epshz, fn_hz, aa, bb, c_p, rho

      if (wtd%numcps+wtd%n_miss.le.0) then
         zshift_new = 0d0
         dftot_dz   = 0d0
         use_pow    = .false.
         return
      endif

      if (wtd%ic%is_left_side()) then
         sgn = -1d0
      else
         sgn =  1d0
      endif

      ! Determine sensible C_hz and fs_rel for new contact patches & near miss patches

      new_c_hz(1:MAX_NUM_CPS) = .false.

      do icp = 1, wtd%numcps+wtd%n_miss
         associate(cp => wtd%allcps(icp)%cp)
         if (cp%c_hz.le.0d0 .or. icp.gt.wtd%numcps) then

            ! use curvatures a1, b1 to get sensible C_hz, assume fs_rel proportional to contact angle

            ic_norm =  1
            ipotcn  = -1
            e_star  = wtd%mater%ga / (1d0 - wtd%mater%nu)
            epshz   = wtd%solv%eps
            fn_hz   = 10d3
         
            call hzcalc3d (e_star, epshz, ipotcn, cp%a1, cp%b1, aa, bb, ic_norm, pen, fn_hz, c_p, rho)
            cp%c_hz = fn_hz / max(1d-9,pen)**(1.5)
            new_c_hz(icp) = .true.
            cp%fs_rel = wtd%fric%fstat_min() * sin(cp%delttr)
            if (wtd%ic%tang.eq.0) cp%fs_rel = 0d0

            if (x_force.ge.4) then
               write(bufout,'(2(a,g12.4),a,i0,2(a,g12.4))') ' e_star=',e_star,', eps=',epshz,           &
                      ', ipotcn=', ipotcn,', a1=',cp%a1, ', b1=',cp%b1
               call write_log(1, bufout)
               write(bufout,'(2(a,f8.3),a,i0,4(a,g12.4))') ' aa=',aa,', bb=',bb,', ic_norm=',ic_norm,   &
                      ', pen=',pen,', fn_hz=',fn_hz, ', cp=',c_p,', rho=',rho
               call write_log(1, bufout)
            endif
         endif
         end associate
      enddo

      ! get current estimate z_shift and value z_shift0 at initial contact

      z_shift  = get_z_shift(wtd)
      z_shift0 = its%itdata(1)%xk

      ! print offset z_shift0 per actual patch or z_shift1 for near miss patches

      if (x_force.ge.2 .or. k.le.-3) then
         write(bufout,'(a,i3,3(a,f10.6))') ' k=',k,': z_shift=',z_shift,', z_ws=',wtd%ws%z,             &
                ', dz_defl=',wtd%trk%dz_defl
         call write_log(1, bufout)

         do icp = 1, wtd%numcps+wtd%n_miss
            associate(cp => wtd%allcps(icp)%cp)
            if (icp.gt.wtd%numcps) then
               write(bufout,'(a,i2,a,f9.5,a,es12.5,2(a,f9.6))') ' mis',icp-wtd%numcps,': alph_i=',      &
                    sgn*cp%delttr,', C_i=',cp%c_hz, ', F_s,i=',cp%fs_rel, ', zshft1=',z_shift+cp%gap_min
            else
               write(bufout,'(a,i2,a,f9.5,a,es12.5,2(a,f9.6))') ' icp',icp,': alph_i=',sgn*cp%delttr,   &
                        ', C_i=', cp%c_hz, ', F_s,i=',cp%fs_rel, ', zshft0=',z_shift+cp%gap_min
            endif
            call write_log(1, bufout)
            end associate
         enddo
      endif

      ! approximate Fn = c_hz * pen^1.5, Fs = fs_rel * Fn
      ! rotate [Fy;Fz] = R * [Fs;Fn], [shft_s;shft_n] = R^T * [shft_y;shft_z]

      zshift_new = max(z_shift, z_shift0+1d-3)  ! initial estimate for Newton iteration

      ! Newton iteration for zshift_new

      res    = 1d0
      it_nr  =  0
      do while(abs(res).gt.eps_fz .and. it_nr.lt.max_nr)

         it_nr    = it_nr + 1
         is_complete = .true.

         ! sum ftot and dftot_dz over contact patches

         if (wtd%ic%norm.eq.1) then
            fz_est   = -wtd%trk%kz_rail * wtd%trk%dz_defl       ! fixed dz_defl
            dftot_dz = 0d0
         else
            dz_defl  =  wtd%ws%z - zshift_new                   ! variable dz_defl
            fz_est   = -wtd%trk%kz_rail * dz_defl
            dftot_dz =  wtd%trk%kz_rail
            if (x_force.ge.4) then
               write(bufout,'(4(a,f12.4))') ' dz_defl=(',wtd%ws%z,'-',zshift_new,'), fz_est=',fz_est,   &
                          ', dftot_dz=',dftot_dz
               call write_log(1, bufout)
            endif
         endif

         do icp = 1, wtd%numcps+wtd%n_miss
            c_hz     = wtd%allcps(icp)%cp%c_hz
            fs_rel   = wtd%allcps(icp)%cp%fs_rel
            alph     = wtd%allcps(icp)%cp%delttr
            pen      = (zshift_new - z_shift - wtd%allcps(icp)%cp%gap_min) * cos(alph)
            if (x_force.ge.4) then
               write(bufout,'(5(a,f12.4))') ' pen=(',zshift_new,'-',z_shift,'-',                        &
                          wtd%allcps(icp)%cp%gap_min,') *',cos(alph),' =',pen
               call write_log(1, bufout)
            endif

            ! determine whether the new estimate is based on true patches with complete C_hz, fs_rel
            if (pen.ge.1d-9) is_complete = is_complete .and. icp.le.wtd%numcps .and. .not.new_c_hz(icp)

            pen      = max(1d-9, pen)
            fz_est   = fz_est +           (fs_rel*sin(alph) + cos(alph)) * c_hz * pen**(1.5d0)
            dftot_dz = dftot_dz + 1.5d0 * (fs_rel*sin(alph) + cos(alph)) * c_hz * pen**(0.5d0) * cos(alph)
            if (x_force.ge.4) then
               write(bufout,'(a,i2,a,f9.6,2(a,f13.2))') ' icp',icp,': pen=',pen,', fz_est=',fz_est,     &
                        ', dfz=',dftot_dz
               call write_log(1, bufout)
            endif
         enddo

         ! determine residual, update, new estimate

         res     = ftarg - fz_est
         dzshift = res / dftot_dz
         zshift_prv = zshift_new
         zshift_new = zshift_new + dzshift

         if (x_force.ge.3) then
            write(bufout,'(a,i2,a,f9.4,2(a,f13.2),2(a,f9.4))') ' it ',it_nr,': zshift=',zshift_prv,     &
                ', fz_est=',fz_est,', res=',res, ', dzshift=',dzshift,', zshift_new=',zshift_new
            call write_log(1, bufout)
         endif

      enddo ! while not converged

      ! tell the iteration method whether to use zshift_new or not

      use_pow      = (.not.brent_its_has_bracket(its) .and. (k.le.1 .or. zshift_new.gt.z_shift) .and. &
                                                           (.not.was_complete .or. .not.is_complete))
      was_complete = is_complete

      ! when using powerlaw, make sure iteration isnt aborted for a small step size

      if (use_pow .and. abs(zshift_new-z_shift).lt.tol_zk) then
         if (zshift_new.gt.z_shift) then
            zshift_new = z_shift + 1.1d0 * tol_zk
         else
            zshift_new = z_shift - 1.1d0 * tol_zk
         endif
      endif

      if (use_pow .and. isnan(zshift_new)) then
         call write_log(' NaN in solve_zshift_powerlaw')
         call abort_run()
      endif

   end subroutine solve_zshift_powerlaw

!------------------------------------------------------------------------------------------------------------

   function get_z_shift(wtd)
!--purpose: get the unknown z_shift for the vertical problem from the contact datastructure
      implicit none
!--function result:
      real(kind=8)      :: get_z_shift
!--subroutine arguments:
      type(t_ws_track)  :: wtd

      get_z_shift = wtd%ws%z - wtd%trk%dz_defl

   end function get_z_shift

!------------------------------------------------------------------------------------------------------------

   subroutine set_z_shift(wtd, z_shift)
!--purpose: copy the unknown z_shift = z_ws - dz_defl for the vertical problem into the contact problem
      implicit none
!--subroutine arguments:
      type(t_ws_track)              :: wtd
      real(kind=8),     intent(in)  :: z_shift

      if (wtd%ic%norm.eq.1) then
         wtd%ws%z        = z_shift + wtd%trk%dz_defl    ! cases 0--2: unknown z_ws
      else
         wtd%trk%dz_defl = wtd%ws%z - z_shift           ! case  3:    unknown dz_defl
      endif

      if (wtd%ic%x_force.ge.5) then
         write(bufout,'(3(a,f12.6))') ' set_z_shift: z_shift=',z_shift,': z_ws=',wtd%ws%z,', dz_defl=', &
                           wtd%trk%dz_defl
         call write_log(1, bufout)
      endif

   end subroutine set_z_shift

!------------------------------------------------------------------------------------------------------------

   function get_f_totz(wtd, fz_cntc_arg)
!--purpose: get the vertical force on the rail for the vertical problem: contact force + spring force
      implicit none
!--function result:
      real(kind=8)           :: get_f_totz
!--subroutine arguments:
      type(t_ws_track)       :: wtd
      real(kind=8), optional :: fz_cntc_arg
!--local variables:
      real(kind=8)           :: ftotz

      ! get contact force F_z(tr) on rail

      if (present(fz_cntc_arg)) then
         ftotz = fz_cntc_arg     ! F_z(tr) may be specified without going through the contact calculation
      else
         ftotz = wtd%ftrk%z()
      endif

      ! add spring force from ground on rail

      ftotz = ftotz - wtd%trk%kz_rail * wtd%trk%dz_defl

      get_f_totz = ftotz

   end function get_f_totz

!------------------------------------------------------------------------------------------------------------

   function get_y_shift(wtd)
!--purpose: get the unknown y_shift for the lateral problem from the contact datastructure
      implicit none
!--function result:
      real(kind=8)      :: get_y_shift
!--subroutine arguments:
      type(t_ws_track)  :: wtd

      get_y_shift = -wtd%trk%dy_defl  ! dy_defl: non-mirrored

   end function get_y_shift

!------------------------------------------------------------------------------------------------------------

   subroutine set_y_shift(wtd, y_shift)
!--purpose: copy the unknown y_shift = -dy_defl for the lateral problem into the contact problem
      implicit none
!--subroutine arguments:
      type(t_ws_track)              :: wtd
      real(kind=8),     intent(in)  :: y_shift

      wtd%trk%dy_defl = -y_shift

      if (wtd%ic%x_force.ge.5) then
         write(bufout,'(3(a,f12.6))') ' set_y_shift: y_shift=',y_shift,': y_ws=',wtd%ws%y,', dy_defl=', &
                           wtd%trk%dy_defl
         call write_log(1, bufout)
      endif

   end subroutine set_y_shift

!------------------------------------------------------------------------------------------------------------

   function get_f_toty(wtd)
!--purpose: get the lateral force on the rail: contact force + spring force
      implicit none
!--function result:
      real(kind=8)           :: get_f_toty
!--subroutine arguments:
      type(t_ws_track)       :: wtd
!--local variables:
      real(kind=8)           :: sgn

      ! get contact force F_y(tr) on rail + spring force from ground on rail

      if (wtd%ic%is_left_side()) then
         sgn = -1d0
      else
         sgn =  1d0
      endif

      get_f_toty = sgn * wtd%ftrk%y() - wtd%trk%ky_rail * wtd%trk%dy_defl

   end function get_f_toty

!------------------------------------------------------------------------------------------------------------

   subroutine wr_contact_fz_brent(wtd, iestim_br, my_ierror)
!--purpose: process the w/r problem for a case with prescribed vertical force with Brents algorithm
      implicit none
!--subroutine arguments:
      type(t_ws_track)              :: wtd
      integer,          intent(in)  :: iestim_br
      integer,          intent(out) :: my_ierror
!--local variables:
      integer, parameter       :: maxit = 100
      logical, parameter       :: use_powerlaw = .true.
      logical                  :: used_bisec, was_complete, use_pow, ldone
      integer                  :: k, x_force, x_locate, ic_return, sub_ierror
      real(kind=8)             :: ftot, ftarg, res_fk, tol_fk, dfk, tol_zk, z_shift0, z_shift, fz_cntc, &
                                  dftot_dz, r_new, aa, bb
      type(t_brent_its)        :: its

      my_ierror = 0
      x_force  = wtd%ic%x_force
      x_locate = wtd%ic%x_locate
      ! if (wtd%meta%itforc_out.eq.2) x_force = 3

      if (max(x_force,x_locate).ge.2) call write_log(' --- Start subroutine wr_contact_fz_brent ---')

      ! case 0--2, N=1: Abort if a zero or negative total force is prescribed

      if (wtd%ic%norm.eq.1 .and. wtd%ws%fz_inp.lt.1d-9) then
         wtd%numcps = 0
         if (x_force.ge.1) then  
            call write_log(' Prescribed total force is zero or negative, skipping computation.')
         endif
         return
      endif

      ! case 1--2, F=3, N=1: massless rail with known forces: solve vertical deflection directly

      if (wtd%ic%force1.eq.3 .and. wtd%ic%norm.eq.1) then
         if (abs(wtd%trk%kz_rail).le.1d-10) then
            wtd%trk%dz_defl = 0d0                                                  ! case 1
         else
            wtd%trk%dz_defl = (wtd%ws%fz_inp + wtd%trk%fz_rail) / wtd%trk%kz_rail  ! case 2
         endif
      endif

      ! case 0--3: Iteration works on equation F_tot,z := F_z(tr) - k_z dz_defl
      !                          using unknown z_shift := z_ws - dz_defl

      ! case 0, F=0: ensure k_z = 0

      if (wtd%ic%force1.ne.3) wtd%trk%kz_rail = 0d0

      ! case 0--1: N=1, k_z=0: ensure F_defl,z = -F_z,inp

      if (wtd%ic%norm.eq.1 .and. abs(wtd%trk%kz_rail).le.1d-10) wtd%trk%fz_rail = -wtd%ws%fz_inp

      ! prepare structure to hold iterates

      ftarg = -wtd%trk%fz_rail
      if (wtd%ic%force1.eq.0 .or. wtd%ic%force1.eq.1) then
         call brent_its_init(its,  ikZDIR, maxit, ftarg)
      else
         call brent_its_init(its, -ikZDIR, maxit, ftarg)
      endif

      ! initialize iteration; set zero point and initial estimate

      !   0 = full initialization, using geometrical analysis + Hertz approximation
      !   1 = continuation, using estimates z_ws from wtd%ws, dz_defl from wtd%trk
      !   2 = continuation, using geometrical analysis + offset z_offset provided

      was_complete = .false.

      if (iestim_br.le.0 .and. .not.use_powerlaw) then

         ! method 0: analyze geometry for overall minimum gap and curvatures, set initial estimates k=0, k=1

         call write_log(' use_powerlaw=F: contact_init_hertz')
         call wr_contact_init_hertz(wtd, its, ftarg, dftot_dz, aa, bb, IMETH_BRENT, x_force, x_locate,  &
                        sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

      elseif (iestim_br.eq.2 .and. .not.use_powerlaw) then

         call write_log(' Internal error: iestim_br=2 based on powerlaw method')

      elseif (iestim_br.eq.0 .or. iestim_br.eq.2) then

         ! method 0 & 2: use initial estimate from powerlaw approach

         if (iestim_br.eq.0) then
            z_shift = 0d0
         else
            z_shift = get_z_shift(wtd)
         endif

         ! analyze geometrical problem for the initial z-position, without actual solving

         if (x_force.ge.3) then
            write(bufout,'(a,i3)')   ' wr_contact_fz_brent: starting iteration k=', 0
            call write_log(1, bufout)
         endif

         call set_z_shift(wtd, z_shift)
         ic_return = wtd%ic%return
         wtd%ic%return = 3

         call wr_contact_pos(wtd, x_locate, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         wtd%ic%return = ic_return

         ! store initial point { z_shift0, r_0 } at which initial contact occurs

         k        = 0
         z_shift0 = z_shift + wtd%ws%gap_min    ! N=0: rail moving up; N=1: wheel moving down
         fz_cntc  = 0d0
         ftot     = get_f_totz(wtd, fz_cntc)
         if (wtd%ic%norm.eq.0) ftot = ftot + wtd%trk%kz_rail * wtd%ws%gap_min
         r_new    = ftot - ftarg
         if (x_force.ge.3) then
            write(bufout,'(3(a,f12.6))') ' z_shift=',z_shift,': gap_min=',wtd%ws%gap_min,', fz_cntc=',  &
                        fz_cntc
            call write_log(1, bufout)
            write(bufout,'(a,f12.6,4(a,f12.3))') ' z_shift=',z_shift,': ftot=',get_f_totz(wtd),         &
                        ', ftot0=',ftot,', ftarg=',ftarg,', res=',r_new
            call write_log(1, bufout)
         endif
         call brent_its_add_iterate(its, k, wtd%numcps, z_shift0, r_new, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, x_force)

         if (x_force.ge.3) then
            write(bufout,'(a,i3)')   ' wr_contact_fz_brent: starting iteration k=', 1
            call write_log(1, bufout)
         endif

         ! estimate z_shift using powerlaw method

         k   = 1
         wtd%meta%itforc_inn = k

         tol_zk  = wtd%solv%eps
         call solve_zshift_powerlaw(wtd, k, its, ftarg, tol_zk, z_shift, dftot_dz, x_force,                &
                        was_complete, use_pow)

         if (use_pow) then
            if (x_force.ge.2) then
               write(bufout,'(a,i2,a,f12.6,a,l2)') ' iteration k=',k,': using z_shift=',z_shift,            &
                        ' from powerlaw method, is_complete=',was_complete
               call write_log(1, bufout)
            endif
         else
            call write_log(' Error: use_pow=.false., no initial estimate?')
         endif

         ! compute problem with new z_shift, get new fk

         call set_z_shift(wtd, z_shift)
         call wr_contact_pos(wtd, x_locate, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

         ! store first point { z_shift1, r_1 } for initial estimate z_shift

         ftot  = get_f_totz(wtd)
         r_new = ftot - ftarg
         call brent_its_add_iterate(its, k, wtd%numcps, z_shift, r_new, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, x_force)

      elseif (iestim_br.eq.1) then

         ! method 1: use initial estimate for z_shift from existing ws%z, trk%dz_defl

         if (x_force.ge.3) then
            write(bufout,'(a,i3)')   ' wr_contact_fz_brent: starting iteration k=',1
            call write_log(1, bufout)
         endif

         ! compute contact problem for z_shift

         wtd%meta%itforc_inn = 1
         call wr_contact_pos(wtd, x_locate, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

         ! store initial point { z_shift0, r_0 } at which initial contact occurs

         k        = 0
         z_shift0 = get_z_shift(wtd) + wtd%ws%gap_min
         ftot     = get_f_totz(wtd, 0d0)
         r_new    = ftot - ftarg
         call brent_its_add_iterate(its, k, wtd%numcps, z_shift0, r_new, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, x_force)

         ! store first point { z_shift1, r_1 } for initial estimate ws%z, trk%dz_defl

         k        = 1
         wtd%meta%itforc_inn = k
         z_shift  = get_z_shift(wtd)
         ftot     = get_f_totz(wtd)
         r_new    = ftot - ftarg
         dftot_dz = wtd%dfz_dzws
         call brent_its_add_iterate(its, k, wtd%numcps, z_shift, r_new, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, x_force)

         was_complete = .true.

      endif

      ! Abort in case of error or if the profiles have no overlap at all

      if (my_ierror.ne.0 .or. .not.wtd%ws%has_overlap) then
         call write_log(' no overlap...')
         call total_forces_moments(wtd)
         ! call brent_its_destroy(its)
         return
      endif

      ! Set tolerances on |fk| and |xk(2)-xk(1)|

      dfk     = abs(its%it_b%rk - its%it_a%rk) / max(1d-6, abs(its%it_b%xk - its%it_a%xk))
      tol_fk  = max(1d-12*dfk, wtd%solv%eps * abs(ftarg))
      tol_zk  = wtd%solv%eps

      ! initialize the Brent loop for obtaining the prescribed vertical force

      k  = its%numit - 1
      used_bisec = .true.

      ! check convergence

      res_fk = its%it_b%rk
      ldone = (abs(res_fk).lt.tol_fk .or. abs(its%it_b%xk-its%it_a%xk).lt.tol_zk .or.                   &
                                                  k.ge.wtd%solv%maxnr .or.  wtd%ic%return.ge.2)

      ! while not "done" do

      do while (.not.ldone .and. my_ierror.eq.0)

         ! increment iteration number

         k = k + 1
         wtd%meta%itforc_inn = k

         if (x_force.ge.3) then
            write(bufout,'(a,i3)')   ' wr_contact_fz_brent: starting iteration k=',k
            call write_log(1, bufout)
         endif

         ! determine new z_shift using power law for each patch: F_n = C_i pen_i^1.5 

         if (use_powerlaw) then
            call solve_zshift_powerlaw(wtd, k, its, ftarg, tol_zk, z_shift, dftot_dz, x_force,          &
                        was_complete, use_pow)
         else
            use_pow = .false.
         endif

         ! fall-back: determine new point z_shift using Brent update

         if (use_pow) then
            if (x_force.ge.2) then
               write(bufout,'(a,i2,a,f11.7,a,l2)') ' iteration k=',k,': using z_shift=',z_shift,        &
                        ' from powerlaw method, is_complete=',was_complete
               call write_log(1, bufout)
            endif
         else
            call brent_set_xnew(k, its, dftot_dz, used_bisec, tol_zk, z_shift, x_force)
            if (x_force.ge.2) then
               write(bufout,'(a,i2,a,f11.7,a)') ' iteration k=',k,': using z_shift=',z_shift,           &
                        ' from Brents method'
               call write_log(1, bufout)
            endif
         endif

         ! compute problem with new z_shift^k, get new fk

         call set_z_shift(wtd, z_shift)
         call wr_contact_pos(wtd, x_locate, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror

         ftot = get_f_totz(wtd)
         r_new = ftot - ftarg

         call brent_its_add_iterate(its, k, wtd%numcps, z_shift, r_new, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, x_force)

         ! check convergence

         res_fk = its%it_b%rk
         ldone = (abs(res_fk).lt.tol_fk .or. abs(its%it_b%xk-its%it_a%xk).lt.tol_zk .or.                &
                        k.ge.wtd%solv%maxnr .or. my_ierror.ne.0 .or. brent_its_has_jump(its, x_force) )

         if (my_ierror.lt.0) then
            write(bufout,'(a,i3,a)') ' An error occurred the in Brent algorithm (',my_ierror,           &
                '), aborting contact solution'
            call write_log(1, bufout)
         endif
      enddo  ! while not done: Brent loop

      ! compute sensitivity at final solution

      wtd%dfz_dzws = brent_sensitivity_k(k, its, x_force)

      if (my_ierror.eq.0 .and. abs(res_fk).ge.tol_fk .and. abs(its%it_b%xk-its%it_a%xk).ge.tol_zk) then
         write(bufout,'(a,i3,a)') ' Brent algorithm stopped with residual > tolerance.'
         call write_log(1, bufout)
         my_ierror = -2
      endif
      ! call brent_its_destroy(its)

   end subroutine wr_contact_fz_brent

!------------------------------------------------------------------------------------------------------------

   subroutine wr_contact_fy_brent(wtd, my_ierror)
!--purpose: process the w/r problem for a case with prescribed forces with Brent for Fy, Brent for Fz
      implicit none
!--subroutine arguments:
      type(t_ws_track)              :: wtd
      integer,          intent(out) :: my_ierror
!--local variables:
      integer, parameter       :: maxit = 100
      type(t_brent_its)        :: its
      logical                  :: has_fz, has_fy, used_bisec, use_pow2, ldone
      integer                  :: k, x_force, x_locate, iestim_br, sub_ierror
      real(kind=8)             :: sgn, ftarg(2), fy_tot, dfy_dy, r_new, tol_xk, tol_fk, res_fk, eps_out, &
                                  eps_inn, alph, y_shift, z_shift, yshift_prv, yshift_pow, zshift_pow

      my_ierror = 0
      x_force  = wtd%ic%x_force
      x_locate = wtd%ic%x_locate
      wtd%ic%x_force = max(0, x_force-4)   ! used in inner loop

      if (x_force.ge.2) call write_log(' --- Start subroutine wr_contact_fy_brent ---')

      ! set the sign to -1 for left and +1 for right rail/wheel combination

      if (wtd%ic%is_left_side()) then
         sgn = -1
      else
         sgn =  1
      endif

      has_fz = (wtd%ic%norm.ge.1 .or. abs(wtd%trk%kz_rail).ge.1d-10)
      has_fy = (wtd%ic%tang.ge.1 .and. wtd%ic%force1.eq.3)

      if (.not.has_fy) then
         call write_log(' INTERNAL ERROR: Brent method for Fy needs T>=1, F=3')
         call abort_run()
      endif
      if (.not.has_fz) then
         call write_log(' ERROR: Brent method for Fy not yet supported for N=0, kz_rail=0')
         call abort_run()
      endif

      ! Solve equations Ftot(xk) = Ftarg,
      !    xk = [ y_shift ] = [  0   - dy_defl ], Ftarg = [ -Fy_rail ]
      !         [ z_shift ] = [ z_ws - dz_defl ]          [ -Fz_rail ]
      !    F_tot = [ F_{y(tr)} - ky * dy_defl ] -- contact force on rail + spring force on rail
      !            [ F_{z(tr)} - kz * dz_defl ]

      ! case 1--2, F=3, N=1: massless rail with known forces: solve vertical deflection directly

      if (wtd%ic%force1.eq.3 .and. wtd%ic%norm.eq.1) then
         if (abs(wtd%trk%kz_rail).le.1d-10) then
            wtd%trk%dz_defl = 0d0                                                  ! case 1
         else
            wtd%trk%dz_defl = (wtd%ws%fz_inp + wtd%trk%fz_rail) / wtd%trk%kz_rail  ! case 2
         endif
      endif

      ! N=1: kz=0 needs compatible Fz_rail = -Fz_cntc, override vertical inputs

      if (wtd%ic%norm.eq.1 .and. abs(wtd%trk%kz_rail).le.1d-10) then
         wtd%trk%fz_rail = -wtd%ws%fz_inp
      endif

      ftarg(1) = -wtd%trk%fy_rail
      ftarg(2) = -wtd%trk%fz_rail

      ! prepare structure to hold iterates

      call brent_its_init(its, ikYDIR, maxit, ftarg(1))

      ! set initial estimate y_shift = 0 and solve for z_shift using inner iteration

      k               = 0
      wtd%meta%itforc_out = k
      eps_out         = wtd%solv%eps
      eps_inn         = max(1d-4, 10d0 * wtd%solv%eps)

      tol_xk          = eps_out
      tol_fk          = max(1d-6, eps_out * abs(ftarg(1)))

      if (wtd%ic%iestim.eq.1) then
         iestim_br    = 1
         y_shift      = get_y_shift(wtd)
      else
         iestim_br    = 0
         y_shift      = 0d0
         call set_y_shift(wtd, y_shift)
      endif

      wtd%solv%eps    = eps_inn
      call wr_contact_fz_brent(wtd, iestim_br, sub_ierror)
      wtd%solv%eps    = eps_out
      if (my_ierror.eq.0 .and. sub_ierror.ne.-2) my_ierror = sub_ierror    ! ignore -2: |res|>tol

      ! store first point { x_1, r_1 } for initial estimate y_shift

      fy_tot = get_f_toty(wtd)
      r_new  = fy_tot - ftarg(1)

      call brent_its_add_iterate(its, k, wtd%numcps, y_shift, r_new, sub_ierror)
      call brent_its_print(k, its, wtd%ic, x_force)

      ldone  = (abs(r_new).lt.tol_fk .or. my_ierror.ne.0)

      ! determine new point y_shift using 2x2 powerlaw method

      if (.not.ldone) then

         call solve_yshift_powerlaw(wtd, k, its, ftarg, yshift_pow, z_shift, x_force, use_pow2)

         ! fall-back: set second guess dy_defl = +/- 1, solve for z_shift using inner iteration

         if (use_pow2) then
            y_shift = yshift_pow
            if (x_force.ge.2) then
               write(bufout,'(a,i2,2(a,f10.6),a)') ' fy, iteration k=',k,': using y_shift=',y_shift,    &
                        ', z_shift=',z_shift,' from powerlaw method'
               call write_log(1, bufout)
            endif
         else
            if (r_new.gt.0d0) then
               y_shift = y_shift - 1d0
            else
               y_shift = y_shift + 1d0
            endif
            if (x_force.ge.2) then
               write(bufout,'(a,i2,a,f10.6)') ' fy, iteration k=',k,': using fall-back y_shift=', y_shift
               call write_log(1, bufout)
            endif
         endif
         call set_y_shift(wtd, y_shift)

         k         = 1
         wtd%meta%itforc_out = k
         iestim_br = 0
         wtd%solv%eps = eps_inn
         call wr_contact_fz_brent(wtd, iestim_br, sub_ierror)
         wtd%solv%eps = eps_out
         if (my_ierror.eq.0 .and. sub_ierror.ne.-2) my_ierror = sub_ierror    ! ignore -2: |res|>tol

         ! store second point { x_2, r_2 } for second guess y_shift

         fy_tot = get_f_toty(wtd)
         r_new  = fy_tot - ftarg(1)

         call brent_its_add_iterate(its, k, wtd%numcps, y_shift, r_new, sub_ierror)
         call brent_its_print(k, its, wtd%ic, x_force)

         ldone  = (abs(r_new).lt.tol_fk .or. my_ierror.ne.0)

      endif

      ! while not "done" do

      used_bisec = .false.
      eps_inn    = eps_out

      do while (.not.ldone .and. my_ierror.eq.0)

         ! increment iteration number

         k = k + 1
         wtd%meta%itforc_out = k
         yshift_prv = y_shift

         if (x_force.ge.3) then
            write(bufout,'(a,i3)')   ' wr_contact_fy_brent: starting iteration k=',k
            call write_log(1, bufout)
         endif

         ! determine new point y_shift using 2x2 powerlaw method

         z_shift  = get_z_shift(wtd)
         call solve_yshift_powerlaw(wtd, k, its, ftarg, y_shift, zshift_pow, x_force, use_pow2)

         ! fall-back: determine new point y_shift using Brent update

         if (use_pow2) then
            z_shift = zshift_pow
            if (x_force.ge.2) then
               write(bufout,'(a,i2,2(a,f11.7),a)') ' fy, iteration k=',k,': using y_shift=',y_shift,    &
                        ', z_shift=',z_shift,' from powerlaw method'
               call write_log(1, bufout)
            endif
         else
            dfy_dy = brent_sensitivity_k(k-1, its, x_force)
            call brent_set_xnew(k, its, dfy_dy, used_bisec, tol_xk, y_shift, x_force)

            if (x_force.ge.2) then
               write(bufout,'(a,i2,2(a,f11.7))') ' fy, iteration k=',k,': using y_shift=',y_shift,      &
                        ' from Brents method, keep z_shift=',z_shift
               call write_log(1, bufout)
            endif
         endif

         ! solve for z_shift using inner iteration, get new fk

         call set_y_shift(wtd, y_shift)
         iestim_br       = 2
         if (wtd%numcps.ge.2 .and. abs(y_shift-yshift_prv).lt.0.1d0) then
            iestim_br = 1
            call set_z_shift(wtd, z_shift)
         elseif (wtd%numcps.eq.1 .and. abs(y_shift-yshift_prv).lt.0.1d0) then
            iestim_br = 1
            alph     = wtd%allcps(1)%cp%delttr
            wtd%ws%z = wtd%ws%z + sgn * sin(alph) * (y_shift - yshift_prv)
         endif

         wtd%solv%eps = eps_inn
         call wr_contact_fz_brent(wtd, iestim_br, sub_ierror)
         wtd%solv%eps = eps_out
         if (my_ierror.eq.0 .and. sub_ierror.ne.-2) my_ierror = sub_ierror    ! ignore -2: |res|>tol

         ! store new point { x_k, r_k } for y_shift

         fy_tot = get_f_toty(wtd)
         r_new  = fy_tot - ftarg(1)

         call brent_its_add_iterate(its, k, wtd%numcps, y_shift, r_new, sub_ierror)
         if (my_ierror.eq.0) my_ierror = sub_ierror
         call brent_its_print(k, its, wtd%ic, x_force)

         ! check convergence

         res_fk = its%it_b%rk
         ldone  = (abs(res_fk).lt.tol_fk .or. abs(its%it_b%xk-its%it_a%xk).lt.tol_xk .or.               &
                        k.ge.wtd%solv%maxnr .or. my_ierror.ne.0 .or. brent_its_has_jump(its, x_force) )

         if (my_ierror.lt.0) then
            write(bufout,'(a,i3,a)') ' An error occurred the in Brent algorithm (',my_ierror,           &
                '), aborting contact solution'
            call write_log(1, bufout)
         endif
      enddo  ! while not done: Brent loop

      if (my_ierror.eq.0 .and. abs(res_fk).ge.tol_fk .and. abs(its%it_b%xk-its%it_a%xk).ge.tol_xk) then
         write(bufout,'(a,i3,a)') ' Brent algorithm (fy) stopped with residual > tolerance.'
         call write_log(1, bufout)
         my_ierror = -2
      endif
      ! call brent_its_destroy(its)

      wtd%ic%x_force = x_force

   end subroutine wr_contact_fy_brent

!------------------------------------------------------------------------------------------------------------

   subroutine solve_yshift_powerlaw(wtd, k, its, ftarg, yshift_new, zshift_new, x_force, use_pow2)
!--purpose: determine new (y_shift,z_shift) using approximate 2x2 power law for each contact patch:
!           F_{n,i} = C_i pen_i^1.5 , F_{s,i} = f_{s,i} * F_{n,i}
      implicit none
!--subroutine arguments:
      type(t_ws_track)                 :: wtd
      integer,           intent(in)    :: k, x_force
      type(t_brent_its), intent(in)    :: its
      real(kind=8),      intent(in)    :: ftarg(2)
      real(kind=8),      intent(out)   :: yshift_new, zshift_new
      logical,           intent(out)   :: use_pow2
!--local variables:
      integer,      parameter :: max_nr = 20
      real(kind=8), parameter :: eps_fz = 1d-3
      logical        :: new_c_hz(MAX_NUM_CPS), ldone
      integer        :: icp, ic_norm, ipotcn, iter, ncontrb
      real(kind=8)   :: y_shift, z_shift, dz_defl, c_hz, fs_rel, alph, pen, sgn, fac,                   &
                        e_star, epshz, fn_hz, aa, bb, c_p, rho, fac_max, dx_n
      real(kind=8)   :: rot(2,2), jac_loc(2,2), jac_rt(2,2), rjrt(2,2), jac_glb(2,2), j_inv(2,2),       &
                        det, res(2), dx(2), fn_est, fs_est, fy_est, fz_est, f_glb(2), tol_brack

      if (wtd%ic%is_left_side()) then
         sgn = -1d0
      else
         sgn =  1d0
      endif

      ! Determine sensible c_hz and fs_rel for new contact patches & near miss patches

      new_c_hz(1:MAX_NUM_CPS) = .false.

      do icp = 1, wtd%numcps+wtd%n_miss
         associate(cp => wtd%allcps(icp)%cp)
         if (cp%c_hz.le.0d0 .or. icp.gt.wtd%numcps) then

            ! use curvatures a1, b1 to get sensible C_hz, assume fs_rel proportional to contact angle

            ic_norm =  1
            ipotcn  = -1
            e_star  = wtd%mater%ga / (1d0 - wtd%mater%nu)
            epshz   = wtd%solv%eps
            fn_hz   = 10d3
         
            call hzcalc3d(e_star, epshz, ipotcn, cp%a1, cp%b1, aa, bb, ic_norm, pen, fn_hz, c_p, rho)
            cp%c_hz = fn_hz / max(1d-9,pen)**(1.5)
            new_c_hz(icp) = .true.
            cp%fs_rel = wtd%fric%fstat_min() * sin(cp%delttr)

            if (x_force.ge.4) then
               write(bufout,'(2(a,g12.4),a,i0,2(a,g12.4))') ' e_star=',e_star,', eps=',epshz,           &
                      ', ipotcn=', ipotcn,', a1=',cp%a1, ', b1=',cp%b1
               call write_log(1, bufout)
               write(bufout,'(2(a,f8.3),a,i0,4(a,g12.4))') ' aa=',aa,', bb=',bb,', ic_norm=',ic_norm,   &
                      ', pen=',pen,', fn_hz=',fn_hz, ', cp=',c_p,', rho=',rho
               call write_log(1, bufout)
            endif
         endif
         end associate
      enddo

      ! get current estimates y_shift, z_shift

      y_shift  = get_y_shift(wtd)
      z_shift  = get_z_shift(wtd)

      if (x_force.ge.2 .or. k.le.-3) then
         write(bufout,'(a,i3,5(a,f11.6))') ' k=',k,': y_shift=',y_shift,', dy_defl=',wtd%trk%dy_defl,   &
                ', z_shift=',z_shift, ', z_ws=',wtd%ws%z, ', dz_defl=',wtd%trk%dz_defl
         call write_log(1, bufout)
      endif

      ! print offset z_shift0 per actual patch or z_shift1 for near miss patches

      if (x_force.ge.2 .or. k.le.-3) then
         do icp = 1, wtd%numcps+wtd%n_miss
            associate(cp => wtd%allcps(icp)%cp)
            if (icp.gt.wtd%numcps) then
               write(bufout,'(a,i2,a,f9.5,a,es12.5,2(a,f9.6))') ' mis',icp-wtd%numcps,': alph_i=',      &
                    sgn*cp%delttr,', C_i=',cp%c_hz, ', F_s,i=',cp%fs_rel, ', zshft1=',z_shift+cp%gap_min
            else
               write(bufout,'(a,i2,a,f9.5,a,es12.5,2(a,f9.6))') ' icp',icp,': alph_i=',sgn*cp%delttr,   &
                        ', C_i=', cp%c_hz, ', F_s,i=',cp%fs_rel, ', zshft0=',z_shift+cp%gap_min
            endif
            call write_log(1, bufout)
            end associate
         enddo
      endif

      ! determine estimate for y_shift, z_shift from 2x2 powerlaw method using Newton iteration

      ! - initial estimate: current values

      yshift_new = y_shift
      zshift_new = z_shift

      if (x_force.ge.2) then
         write(bufout,'(a,i3,a,71x,a,2f14.7)') ' solve_yshift: k=',iter,':','[yz]_new=', yshift_new,    &
                zshift_new
         call write_log(1, bufout)
      endif

      ldone = .false.
      iter  = 0
      do while(iter.lt.max_nr .and. .not.ldone)

         iter    = iter + 1
         ncontrb = 0

         ! - form forces and Jacobian per patch & add to overall forces and Jacobian

         f_glb(1:2)       = 0d0
         jac_glb(1:2,1:2) = 0d0

         do icp = 1, wtd%numcps + wtd%n_miss

            ! c_hz, fs_rel, alph, pen: local results for mirrored configuration, yz_shift: non-mirrored

            c_hz     = wtd%allcps(icp)%cp%c_hz
            fs_rel   = wtd%allcps(icp)%cp%fs_rel
            alph     = wtd%allcps(icp)%cp%delttr
            if (icp.le.wtd%numcps) then
               pen   = wtd%allcps(icp)%cp%gd%kin%pen + sgn * (yshift_new - y_shift) * -sin(alph)        &
                                                     +       (zshift_new - z_shift) *  cos(alph)
            else
               pen   = -wtd%allcps(icp)%cp%gap_min * cos(alph)                                          &
                             + sgn * (yshift_new - y_shift) * -sin(alph)                                &
                             +       (zshift_new - z_shift) *  cos(alph)
            endif
            if (x_force.ge.4) then
               write(bufout,'(9(a,f10.4))') ' pen=',pen,' =',-wtd%allcps(icp)%cp%gap_min*cos(alph),     &
                     ' + (',zshift_new,' -',z_shift,') *',cos(alph),' - ',sgn,' * (',yshift_new,        &
                     ' -',y_shift,') *', sin(alph)
               call write_log(1, bufout)
            endif
            if (pen.gt.1d-6) ncontrb = ncontrb + 1
            pen      = max(1d-6, pen)

            ! forces and Jacobian in local coordinates, mirrored configuration

            fn_est   = c_hz * pen**1.5d0
            fs_est   = fs_rel * fn_est

!           write(bufout,'(a,i3,4(a,f14.6))') ' icp=',icp,': alph=',sgn*alph,', pen=',pen,' fn=',fn_est, &
!                    ', fs=',fs_est
!           call write_log(1, bufout)

            jac_loc(1,1) = 0d0  ! jac_loc within mirrored configuration
            jac_loc(2,1) = 0d0
            jac_loc(1,2) = 1.5d0 * c_hz * sqrt(max(0d0, pen)) * fs_rel
            jac_loc(2,2) = 1.5d0 * c_hz * sqrt(max(0d0, pen))

            ! R = [ ca  -sa ],   R^T = [  ca  sa ]
            !     [ sa   ca ]          [ -sa  ca ]

            rot(1,1)    =      cos(alph)   ! rot within global, non-mirrored configuration
            rot(2,1)    =  sgn*sin(alph)
            rot(1,2)    = -sgn*sin(alph)
            rot(2,2)    =      cos(alph)

            fy_est      = rot(1,1) * sgn*fs_est + rot(1,2) * fn_est  ! f_glb: non-mirrored configuration
            fz_est      = rot(2,1) * sgn*fs_est + rot(2,2) * fn_est

            ! Jac * R^T = [ 0  J12 ] * [  ca  sa ],  global, non-mirrored configuration
            !             [ 0  J22 ] * [ -sa  ca ]

            jac_rt(1,1) = sgn*jac_loc(1,2) * rot(1,2)
            jac_rt(2,1) =     jac_loc(2,2) * rot(1,2)
            jac_rt(1,2) = sgn*jac_loc(1,2) * rot(2,2)
            jac_rt(2,2) =     jac_loc(2,2) * rot(2,2)

            ! R * (Jac * R^T) = [ R11  R12 ] * [ JR11  JR12 ],  global, non-mirrored configuration
            !                   [ R21  R22 ] * [ JR21  JR22 ]

            rjrt(1,1)   = rot(1,1) * jac_rt(1,1) + rot(1,2) * jac_rt(2,1)
            rjrt(2,1)   = rot(2,1) * jac_rt(1,1) + rot(2,2) * jac_rt(2,1)
            rjrt(1,2)   = rot(1,1) * jac_rt(1,2) + rot(1,2) * jac_rt(2,2)
            rjrt(2,2)   = rot(2,1) * jac_rt(1,2) + rot(2,2) * jac_rt(2,2)

            ! add to total forces and total Jacobian

            f_glb(1)    = f_glb(1) + fy_est
            f_glb(2)    = f_glb(2) + fz_est

            jac_glb(1:2,1:2) = jac_glb(1:2,1:2) + rjrt(1:2,1:2)

            if (x_force.ge.5) then
               write(bufout,'(a,i3,8(a,f12.1),a,/,18x,8(a,f12.1),a)') ' icp=',icp,                      &
                ': F_loc = [', fs_est, '], F_glb = [', fy_est, '], J_loc = [', jac_loc(1,1),',',        &
                                       jac_loc(1,2), '], J*R^T = [', jac_rt(1,1), ',',jac_rt(1,2),      &
                                                     '], J_glb = [', rjrt(1,1), ',',rjrt(1,2), ']',     &
                          '[', fn_est, ']          [', fz_est, '],         [', jac_loc(2,1),',',        &
                                       jac_loc(2,2), ']          [', jac_rt(2,1), ',',jac_rt(2,2),      &
                                                     ']          [', rjrt(2,1), ',',rjrt(2,2), ']'
               call write_log(2, bufout)
            endif

         enddo ! icp

         ! add spring force (dy_defl = -yshift_new) using overall track orientation

         f_glb(1)     = f_glb(1) + wtd%trk%ky_rail * yshift_new
         jac_glb(1,1) = jac_glb(1,1) + wtd%trk%ky_rail
         if (ncontrb.le.0) jac_glb(2,2) = jac_glb(2,2) + 1d2 * wtd%trk%ky_rail

         if (wtd%ic%norm.eq.1) then
            dz_defl = wtd%trk%dz_defl
         else
            dz_defl = wtd%ws%z - zshift_new
         endif
         f_glb(2)     = f_glb(2) - wtd%trk%kz_rail * dz_defl
         if (wtd%ic%norm.eq.0) jac_glb(2,2) = jac_glb(2,2) + wtd%trk%kz_rail

!        write(bufout,'(2(a,2f14.6))') ' ftrue=',get_f_toty(wtd),get_f_totz(wtd),', rtrue=',            &
!            get_f_toty(wtd)-ftarg(1), get_f_totz(wtd)-ftarg(2)
!        call write_log(1, bufout)
!        write(bufout,'(2(a,2f14.6))') ' fest= ', f_glb(1),f_glb(2),', rest= ',f_glb(1)-ftarg(1),       &
!            f_glb(2)-ftarg(2)
!        call write_log(1, bufout)

         ! determine residual estimate

         res(1) = f_glb(1) - ftarg(1)
         res(2) = f_glb(2) - ftarg(2)
         ldone  = (max(abs(res(1)), abs(res(2))).lt.eps_fz)

         ! solve J * dx = -res

         det = jac_glb(1,1) * jac_glb(2,2) - jac_glb(2,1) * jac_glb(1,2)
         j_inv(1,1) =  jac_glb(2,2) / det
         j_inv(2,1) = -jac_glb(2,1) / det
         j_inv(1,2) = -jac_glb(1,2) / det
         j_inv(2,2) =  jac_glb(1,1) / det

         if (x_force.ge.5) then
            write(bufout,'(a,i3,a,g12.4)') '   k=',iter,': det=',det
            call write_log(1, bufout)
            write(bufout,'(a,i3,a,25x,5(a,f12.1),a,/,42x,5(a,f12.1),a)') '   k=',iter,':',              &
               'F_tot = [', f_glb(1),'], J_tot = [', jac_glb(1,1), ',', jac_glb(1,2),                   &
                                     '], J_inv = [', j_inv(1,1)*1d9, ',', j_inv(1,2)*1d9,'] * 1d-9',    &
                       '[', f_glb(2),']          [', jac_glb(2,1), ',', jac_glb(2,2),                   &
                                     ']          [', j_inv(2,1)*1d9, ',', j_inv(2,2)*1d9,']'
            call write_log(2, bufout)
         endif

         dx(1) = -( j_inv(1,1) * res(1) + j_inv(1,2) * res(2) )
         dx(2) = -( j_inv(2,1) * res(1) + j_inv(2,2) * res(2) )

         ! apply back-tracking to ensure pen>=0

         fac_max = 1d0

         do icp = 1, wtd%numcps
            alph    = wtd%allcps(icp)%cp%delttr
            dx_n    = sgn * dx(1) * -sin(alph) + dx(2) * cos(alph)
            pen     = wtd%allcps(icp)%cp%gd%kin%pen + sgn * (yshift_new - y_shift) * -sin(alph)         &
                                                    +       (zshift_new - z_shift) *  cos(alph)
            if (.true. .and. dx_n.lt.0d0) then
               fac_max = min( fac_max, pen / abs(dx_n) )
!              write(bufout,'(a,i2,3(a,f12.6))') ' icp',icp,': dx_n=',dx_n,' < 0, pen=',pen,            &
!                    ', fac_max=', pen/abs(dx_n)
!              call write_log(1, bufout)
            endif
         enddo

         fac = fac_max
         yshift_new = yshift_new + fac * dx(1)
         zshift_new = zshift_new + fac * dx(2)

         if (x_force.ge.2) then
            write(bufout,'(a,i3,a,2f14.3,2(a,2f14.7))') ' solve_yshift: k=',iter,': res=',res(1:2),     &
                        ', d[yz]=',dx(1:2),', [yz]_new=',yshift_new,zshift_new
            call write_log(1, bufout)
         endif

      enddo ! while(.not.done)

      ! tell the iteration method whether to use [yz]shift_new or not

      use_pow2 = (ldone .and. k.le.10)        !   (.not.was_complete .or. .not.is_complete))

      ! for one-patch cases, reject large steps in the initial iterations

      if (ncontrb.le.1 .and. k.le.3 .and. (abs(zshift_new-z_shift).gt.1d0 .or.                          &
                                           abs(yshift_new-y_shift).gt.5d0)) use_pow2 = .false.

      ! for one-patch cases, switch to Brent as soon as a bracket is found

      if (ncontrb.le.1 .and. brent_its_has_bracket(its)) use_pow2 = .false.

      ! for two-patch cases, require new yshift in interior of bracket, away from the sides

      if (use_pow2 .and. brent_its_has_bracket(its)) then
         tol_brack = 0.1d0 * (its%it_1%xk - its%it_0%xk)
         use_pow2 = (yshift_new.ge.its%it_0%xk+tol_brack .and. yshift_new.le.its%it_1%xk-tol_brack)
      endif
            
      if (isnan(yshift_new) .or. isnan(zshift_new)) then
         call write_log(' NaN in solve_yshift_powerlaw')
         call abort_run()
      endif

   end subroutine solve_yshift_powerlaw

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
      !!! call brent_its_init(its, ikXDIR, maxit, ftarg) ??

      ! preparations for initial estimate for z-position

      if (has_fz) then

         k = k + 1
         if (x_force.ge.3) then
            write(bufout,'(a,i3)')   ' wr_contact_fx_secant: starting iteration k=',k
            call write_log(1, bufout)
         endif

         ! analyze geometry for overall minimum gap and curvatures, set initial estimate k=0 for Fz

         call wr_contact_init_hertz(wtd, its, wtd%ws%fz_inp, dfz_dzws, aa, bb, IMETH_SECANT, x_force, &
                        x_locate, sub_ierror)
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
      call wr_contact_fz_brent(wtd, iestim_br, sub_ierror)
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

            call wr_contact_fz_brent(wtd, iestim_br, sub_ierror)
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
         ! str16(1) = fmt_gs(12,6,6, xk(1))
         ! str16(2) = fmt_gs(13,8,8, xk(2))
         ! write(str16(3), '(f13.8)') fk(1) - ftarg(1)
         ! write(str16(4), '(f13.8)') fk(2) - ftarg(2)
         ! write(str16(5), '(f11.6)') bk(2,2)
         write(bufout,8000) k, xk(1), xk(2), fk(2)-ftarg(2), bk(2,2)
         call write_log(1, bufout)
      endif

 8000 format(2x, i4,', SC, z, omg_ws:',f12.6,f13.8,', res: ',f13.8,', dFx/domg: ',f11.6)

   end subroutine secant_fx_print

!------------------------------------------------------------------------------------------------------------

end module m_wr_totforce

