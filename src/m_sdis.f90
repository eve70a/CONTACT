!------------------------------------------------------------------------------------------------------------
! m_sdis - basic tasks (preparations) for a basic contact problem (module 3)
!
! Copyright 1979-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------

module m_sdis

   use m_hierarch_data

   implicit none
   private

   public  init_curr_grid
   private create_gd_gridfunc
   public  check_roll_stepsize
   public  set_prev_data
   public  set_norm_rhs
   public  set_tang_rhs
   public  init_curr_data
   private sdec_pn
   private eldiv0

contains

!------------------------------------------------------------------------------------------------------------

   subroutine init_curr_grid(ic, potcon_inp, potcon_cur, cgrid, geom, outpt)
!--purpose: prepare potcon and grid for current time, resize grid-functions correspondingly
      implicit none
!--subroutine arguments:
      type(t_ic)       :: ic
      type(t_potcon)   :: potcon_inp, potcon_cur
      type(t_grid)   :: cgrid
      type(t_geomet)   :: geom
      type(t_output)   :: outpt
!--local variables :
      type(t_grid)     :: cgrid_new

      ! if "previous time" won't be needed:

      if (ic%pvtime.eq.2 .and. ic%iestim.eq.0) then

         ! copy inputs, construct contact grid for current time (uniform grid for element centers)

         ! call write_log(' init_curr_grid: create new grid, gridfnc...')
         call potcon_copy(potcon_inp, potcon_cur)
         call potcon_cgrid(potcon_cur, cgrid)

         ! create (reallocate) grid-functions used in hierarchical data-structure

         call create_gd_gridfunc(cgrid, geom, outpt)

      else

         ! copy inputs, construct contact grid for current time (uniform grid for element centers)

         ! call write_log(' init_curr_grid: keep previous time, resize...')
         call potcon_copy(potcon_inp, potcon_cur)
         call potcon_cgrid(potcon_cur, cgrid_new)

         ! resize grid-functions used in hierarchical data-structure, keeping existing data

         call gd_resize_gridfunc(cgrid, cgrid_new, geom, outpt, .true.)

         ! copy new grid into gd

         call grid_copy(cgrid_new, cgrid)

      endif

   end subroutine init_curr_grid

!------------------------------------------------------------------------------------------------------------

   subroutine create_gd_gridfunc(cgrid, geom, outpt)
!--purpose: create or reallocate grid-functions used in hierarchical data-structure
      implicit none
!--subroutine arguments:
      type(t_grid),  target :: cgrid
      type(t_geomet)        :: geom
      type(t_output)        :: outpt

      ! (re-)allocate array for element division

      call eldiv_new(outpt%igs, cgrid)

      ! (re-)allocate hs1 at the appropriate size

      call gf3_new(geom%hs1, 'geom%hs1', cgrid)

      ! (re-)allocate arrays for friction coefficients, tractions, displacements and shift

      call gf3_new(outpt%mus, 'outpt%mus', cgrid, outpt%igs)
      call gf3_new(outpt%ps,  'outpt%ps',  cgrid, outpt%igs)
      call gf3_new(outpt%us,  'outpt%us',  cgrid, outpt%igs)
      call gf3_new(outpt%ss,  'outpt%ss',  cgrid, outpt%igs)

      ! (re-)allocate arrays for plastic deformation, yield point, temperature

      call gf3_new(outpt%taucs, 'outpt%taucs', cgrid, outpt%igs)
      call gf3_new(outpt%upls,  'outpt%upls',  cgrid, outpt%igs)
      call gf3_new(outpt%temp1, 'outpt%temp1', cgrid, outpt%igs)
      call gf3_new(outpt%temp2, 'outpt%temp2', cgrid, outpt%igs)

      ! (re-)allocate grid-functions for previous time

      call eldiv_new(outpt%igv, cgrid)
      call gf3_new(geom%hv1,    'geom%hv1',    cgrid)
      call gf3_new(outpt%muv,   'outpt%muv',   cgrid, outpt%igv)
      call gf3_new(outpt%pv,    'outpt%pv',    cgrid, outpt%igv)
      call gf3_new(outpt%uv,    'outpt%uv',    cgrid, outpt%igv)
      call gf3_new(outpt%sv,    'outpt%sv',    cgrid, outpt%igv)
      call gf3_new(outpt%taucv, 'outpt%taucv', cgrid, outpt%igv)
      call gf3_new(outpt%uplv,  'outpt%uplv',  cgrid, outpt%igv)

   end subroutine create_gd_gridfunc

!------------------------------------------------------------------------------------------------------------

   subroutine check_roll_stepsize(ic, solv, kin, cgrid)
!--purpose: check and enforce restrictions the rolling direction chi and step size dq
      implicit none
!--subroutine arguments:
      type(t_ic)       :: ic
      type(t_solvers)  :: solv
      type(t_kincns)   :: kin
      type(t_grid)     :: cgrid
!--local variables :
      logical  :: is_fastsm, is_roll, is_ssrol

      is_fastsm = ic%mater.eq.2 .or. ic%mater.eq.3
      is_roll   = ic%tang.eq.2 .or. ic%tang.eq.3
      is_ssrol  = ic%tang.eq.3

      ! force CHI = 0 when using visco-elastic materials (T=3, M=1)

      if (is_ssrol .and. ic%mater.eq.1 .and. abs(kin%chi).gt.0.001d0) then
         write (bufout, 4002) kin%chi, 'visco-elastic materials'
         call write_log(2, bufout)
         kin%chi = 0d0
      endif

      if (is_ssrol .and. is_fastsm) then

         ! set DQ == DX when using steady rolling w. Fastsim (T=3, M=2/3)
         ! force CHI = 0 or pi when using steady rolling w. Fastsim (T=3, M=2/3)

         kin%dq  = cgrid%dx
         if (abs(kin%chi).gt.0.01d0 .and. abs(kin%chi-pi).gt.0.01d0) then
            write (bufout, 4002) kin%chi, 'Fastsim'
            call write_log(2, bufout)
            kin%chi = 0d0
      endif

      elseif ((solv%gausei_eff.ne.2 .and. is_ssrol) .or. ic%mater.eq.4) then

         ! force CHI = 0 or pi when using steady rolling w. SteadyGS (T=3, G<>2) or when M=4 (T=1,3)
         ! force DQ == DX when using steady rolling w. SteadyGS (T=3, G<>2), also when M=1 (visc.)

         if (abs(kin%chi).gt.0.01d0 .and. abs(kin%chi-pi).gt.0.01d0) then
            if (ic%ilvout.ge.1) then
               if (solv%gausei_eff.ne.2 .and. is_ssrol) then
                  write (bufout, 4002) kin%chi, 'steady state rolling with solver SteadyGS'
               else
                  write (bufout, 4002) kin%chi, 'solver ConvexGS with plasticity'
               endif
               call write_log(2, bufout)
            endif
            kin%chi = 0d0
         endif
         if (abs(kin%dq-cgrid%dx).gt.0.01*cgrid%dx) then
            if (ic%ilvout.ge.1) then
               if (solv%gausei_eff.ne.2 .and. is_ssrol) then
                  write (bufout, 4003) cgrid%dx, kin%dq, 'steady state rolling with solver SteadyGS'
      else
                  write (bufout, 4003) cgrid%dx, kin%dq, 'solver ConvexGS with plasticity'
      endif
               call write_log(2, bufout)
            endif
         endif
         kin%dq  = cgrid%dx
      endif

 4002 format (' Warning: using CHI = 0 instead of',f6.1,' rad'/, 17x,                                   &
              'for steady state rolling with ',a,'.')
 4003 format (' Warning: Using DQ = DX =',g12.4, ' (instead of',g12.4,')', /, 17x, 'for ',a,'.')

      ! set appropriate chi, dq, veloc for shifts

      if (.not.is_roll) then
         kin%chi   = 0d0
         kin%dq    = 1d0
         kin%veloc = 1d0
      endif

      ! compute the time step size, with dq in [mm], veloc in [mm/s]

      if (is_roll .or. ic%mater.ne.5) then
         kin%dt = kin%dq / max(1d-10, kin%veloc)
      endif

   end subroutine check_roll_stepsize

!------------------------------------------------------------------------------------------------------------

   subroutine set_prev_data(ic, geom, mater, fric, kin, outpt)
!--purpose: copy/set variables concerning the previous time instance
      implicit none
!--subroutine arguments:
      type(t_ic)       :: ic
      type(t_geomet)   :: geom
      type(t_material) :: mater
      type(t_friclaw)  :: fric
      type(t_kincns)   :: kin
      type(t_output)   :: outpt
!--local variables :
      logical          :: is_roll

      is_roll = (ic%tang.eq.2 .or. ic%tang.eq.3)

      associate(hs1   => geom%hs1,    hv1   => geom%hv1,    igs   => outpt%igs,   igv   => outpt%igv,   &
                mus   => outpt%mus,   muv   => outpt%muv,   ps    => outpt%ps,    pv    => outpt%pv,    &
                us    => outpt%us,    uv    => outpt%uv,    ss    => outpt%ss,    sv    => outpt%sv,    &
                upls  => outpt%upls,  uplv  => outpt%uplv,  taucs => outpt%taucs, taucv => outpt%taucv)

      ! shift data from current time to previous time
      !  P = 0: pure transient sequence, pv = ps
      !  P = 1: sequence for Norm only, pv_n = ps_n, pv_t = 0
      !  P = 2: initiation of contact, pv = 0
      !  P = 3: obsolete? pv untouched

      ! Element division igv:

      if (ic%pvtime.eq.0) then
         call eldiv_copy(igs, igv, ikALL)
      elseif (ic%pvtime.eq.1) then
         call eldiv_copy(igs, igv, ikZDIR)
      elseif (ic%pvtime.eq.2) then
         call eldiv_exter(igv)
      endif
      call areas(igv)

      ! right-hand-side hv:

      if (ic%pvtime.eq.0) then
         call gf3_copy(AllElm, hs1, hv1, ikALL)
      elseif (ic%pvtime.eq.1) then
         call gf3_copy(AllElm, hs1, hv1, ikZDIR)
      elseif (ic%pvtime.eq.2) then
         call gf3_set (AllElm, 0d0, hv1, ikALL)
      endif

      ! Normal part of pv, fn:

      if (ic%pvtime.le.1) then
         call gf3_copy(AllElm, ps, pv, ikZDIR)
         kin%fprev(ikZDIR) = kin%fcntc(ikZDIR)
      elseif (ic%pvtime.eq.2) then
         call gf3_set(AllElm, 0d0, pv, ikZDIR)
         kin%fprev(ikZDIR) = 0d0
      endif

      ! Normal part of uv:

      if (ic%pvtime.le.1 .and. .not.is_roll) then
         call gf3_copy(AllElm, us, uv, ikZDIR)
      elseif (ic%pvtime.le.1) then
         call gf3_dq_shift(kin%chi, kin%dq, 0d0, us, uv, ikZDIR)
      elseif (ic%pvtime.eq.2) then
         call gf3_set(AllElm, 0d0, uv, ikZDIR)
      endif

      ! Tangential part of pv: copy from current to new grid:

      if (ic%pvtime.eq.0) then
         call gf3_copy(AllElm, ps, pv, ikTANG)
         kin%fprev(ikXDIR) = kin%fcntc(ikXDIR)
         kin%fprev(ikYDIR) = kin%fcntc(ikYDIR)
      elseif (ic%pvtime.le.2) then
         call gf3_set(AllElm, 0d0, pv, ikTANG)
         kin%fprev(ikXDIR) = 0d0
         kin%fprev(ikYDIR) = 0d0
      endif

      ! Friction coefficients muv of previous time instance, Tangential parts of sv, uv:
      !    PvTime=0: pure transient sequence, muv=mus (interpolated).
      !    PvTime=1: sequence for Norm only, muv=fstat
      !    PvTime=2: new initiation of contact, muv=fstat.
      !    PvTime=3: muv untouched.

      if (ic%pvtime.eq.0 .and. .not.is_roll) then
         ! Shift: copying from the current to the new grid
         call gf3_copy(AllElm, mus,   muv,   ikXDIR)   ! scalar
         call gf3_copy(AllElm, taucs, taucv, ikXDIR)   ! scalar
         call gf3_copy(AllElm, ss,    sv,    ikTANG)
         call gf3_copy(AllElm, us,    uv,    ikTANG)
         call gf3_copy(AllElm, upls,  uplv,  ikTANG)
      elseif (ic%pvtime.eq.0) then
         ! Rolling: Interpolation from the current to the new grid, with b.c. at sides of pot.con
         ! TODO: in T=3, uv, sv, etc. should be considered outputs of the previous step, filled in Panag.
         call gf3_dq_shift(kin%chi, kin%dq, fric%nvf, fric%fstat_arr, mus, muv, ikXDIR)
         call gf3_dq_shift(kin%chi, kin%dq, mater%tau_c0, taucs, taucv, ikXDIR)
         call gf3_dq_shift(kin%chi, kin%dq, 0d0, ss,    sv,    ikTANG)
         call gf3_dq_shift(kin%chi, kin%dq, 0d0, us,    uv,    ikTANG)
         call gf3_dq_shift(kin%chi, kin%dq, 0d0, upls,  uplv,  ikTANG)
      elseif (ic%pvtime.le.2) then
         call gf3_copy_xdir(AllElm, fric%nvf, fric%fstat_arr, muv, ikTANG)
         call gf3_set(AllElm, mater%tau_c0, taucv, ikTANG)
         call gf3_set(AllElm, 0d0, sv,    ikTANG)
         call gf3_set(AllElm, 0d0, uv,    ikTANG)
         call gf3_set(AllElm, 0d0, uplv,  ikTANG)
      endif

      end associate
   end subroutine set_prev_data

!------------------------------------------------------------------------------------------------------------

   subroutine set_norm_rhs (cgrid, geom)
!--purpose: construct the Undeformed Distances hs1(i,1) from ibase, iplan, prmudf, prmpln
!           ibase = 1: Undeformed Distance Quadratic in x and y.
!                   2: Undeformed Distance Circular in x, pointwise given in y.
!                   3: Undeformed Distance Difference of two sines, circ. in x.
!                   9: Undeformed Distance Specified per Element.
!           iplan = 1: Unrestricted planform.
!                   2: Elliptic planform.
!                   3: Planform: Union of two rectangles.
!                   4: Weighted interaction between patches.
!           prmudf describes the gap between the profiles (1) and (2)
      implicit none
!--subroutine arguments:
      type(t_grid)           :: cgrid
      type(t_geomet)         :: geom
!--local variables:
      integer, parameter :: idebug = 0
      integer      :: ii, mleft, nplan
      logical      :: z1, z2
      real(kind=8) :: a, dy1, rm, xm, y1, yleft, yn, rc

      associate(prmudf => geom%prmudf, prmpln => geom%prmpln, hs1    => geom%hs1)

      if (geom%ibase.eq.1) then

         ! ibase==1: Quadratic undeformed distance:

         do ii = 1, cgrid%ntot
            hs1%vn(ii) = prmudf(1) * cgrid%x(ii)**2 +                                               &
                       prmudf(2) * cgrid%x(ii) * cgrid%y(ii) +                                      &
                       prmudf(3) * cgrid%y(ii)**2 +                                                 &
                       prmudf(4) * cgrid%x(ii) +                                                    &
                       prmudf(5) * cgrid%y(ii) +                                                    &
                       prmudf(6)
         enddo

      elseif (geom%ibase.eq.2) then

         ! ibase==2: circular in X, pointwise given in Y in nn points

         ! nn  == number of points in Y-direction, read by input
         ! xm  == x-coordinate of the center line of the circles
         ! rm  == radius of curvature in X-direction
         ! y1  == initial Y-Coordinate
         ! dy1 == increment of Y-Coordinate between successive points
         ! yn  == Y-Coordinate of last given value

         xm  = prmudf(2)
         rm  = prmudf(3)
         y1  = prmudf(4)
         dy1 = prmudf(5)
         yn  = y1 + (geom%nn-1)*dy1

         do 22 ii = 1, cgrid%ntot

            ! Ymin <= Y <= Ymax: interpolate between py(mleft) and py(mleft+1)
            ! Y < Ymin or Y > Ymax: extrapolate

            if (cgrid%y(ii).lt.y1) then
               yleft = y1
               mleft = 1
            elseif (cgrid%y(ii).ge.yn) then
               yleft = yn - dy1
               mleft = geom%nn - 1
            else
               mleft = int ((cgrid%y(ii) - y1) / dy1) + 1
               yleft = y1 + (mleft-1) * dy1
            endif

            ! mleft = index of left point in prmudf

            mleft = mleft + 5

            ! Compute the hs1-value at (xm, Y(ii))
            !    basic equation of a straight line: y0 + r.c. * (y-y0)

            rc = (prmudf(mleft+1) - prmudf(mleft)) / dy1
            hs1%vn(ii) = prmudf(mleft) + rc * (cgrid%y(ii) - yleft)

            ! Compute hs1-values at (x(ii), y(ii)), linearized

            hs1%vn(ii) = hs1%vn(ii) + (cgrid%x(ii)-xm)**2 / (2d0 * rm)

 22      continue

      elseif (geom%ibase.eq.3) then

         ! ibase==3: Undeformed distance is the difference of two sines:
         !    B1*sin(B2*(x-B3)) - B4*sin(B5*(x-B6)) - PEN + x**2/B7 + y**2/B8

         do ii = 1, cgrid%ntot
            hs1%vn(ii) =                                                                        &
                    prmudf(1) * sin( prmudf(2) * (cgrid%x(ii) - prmudf(3)) )                        &
                  - prmudf(4) * sin( prmudf(5) * (cgrid%x(ii) - prmudf(6)) )                        &
                  + cgrid%x(ii)**2 / prmudf(7) + cgrid%y(ii)**2 / prmudf(8)
         enddo

      elseif (geom%ibase.eq.9) then

         ! ibase==9: Undeformed distance specified per element in prmudf

         ! copy the values from the input-array to the problem array hs1

         do ii = 1, cgrid%ntot
            hs1%vn(ii) = prmudf(ii)
         enddo

      endif

      ! Restrict contact to planform - if provided

      if (geom%iplan.eq.1) then

         ! iplan == 1: Unrestricted planform

      elseif (geom%iplan.eq.2) then

         ! iplan == 2: Quadratic planform

         do ii = 1, cgrid%ntot
            a = prmpln(1) * cgrid%x(ii)**2 + prmpln(2) * cgrid%x(ii)*cgrid%y(ii) +                      &
                prmpln(3) * cgrid%y(ii)**2 + prmpln(4) * cgrid%x(ii) +                                  &
                prmpln(5) * cgrid%y(ii)    + prmpln(6)
            if (a.ge.0) hs1%vn(ii) = 1d30
         enddo

      elseif (geom%iplan.eq.3) then

         ! iplan == 3: Planform is union of two rectangles

         do ii = 1, cgrid%ntot

            ! z1 = ( (x,y) in rectangle 1 )

            z1 = (prmpln(1).le.cgrid%x(ii) .and. cgrid%x(ii).le.prmpln(2)) .and.                        &
                 (prmpln(3).le.cgrid%y(ii) .and. cgrid%y(ii).le.prmpln(4))

            ! z2 = ( (x,y) in rectangle 2 )

            z2 = (prmpln(5).le.cgrid%x(ii) .and. cgrid%x(ii).le.prmpln(6)) .and.                        &
                 (prmpln(7).le.cgrid%y(ii) .and. cgrid%y(ii).le.prmpln(8))
            if (.not.(z1 .or. z2)) hs1%vn(ii) = 1d30
         enddo

      endif

      ! check that the planform does not exclude everything

      if (geom%iplan.eq.2 .or. geom%iplan.eq.3) then
         nplan = 0
         do ii = 1, cgrid%ntot
            if (hs1%vn(ii).lt.1d20) nplan = nplan + 1
         enddo
         if (nplan.le.0) then
            write(bufout,*) 'ERROR: all elements are excluded by the planform.'
            call write_log(1, bufout)
            call abort_run()
         endif
      endif

      end associate
   end subroutine set_norm_rhs

!------------------------------------------------------------------------------------------------------------

   subroutine set_tang_rhs(ic, mater, cgrid, kin, geom)
      implicit none
!--purpose: set fixed part of the Right-hand-side of tangential problem.
!           exrhs is the extra term in the right-hand-side, input by the user.
!--subroutine arguments :
      type(t_ic)              :: ic
      type(t_material)        :: mater
      type(t_grid)            :: cgrid
      type(t_kincns)          :: kin
      type(t_geomet)          :: geom
!--local variables:
      integer      :: ii
      logical                   :: is_roll
      real(kind=8)              :: cc, sc, facx, facy

      associate(chi    => kin%chi,    dq     => kin%dq,     facphi => kin%facphi,                       &
                cksi   => kin%cksi,   ceta   => kin%ceta,   cphi   => kin%cphi,                         &
                spinxo => kin%spinxo, spinyo => kin%spinyo, hs1    => geom%hs1, exrhs  => geom%exrhs)

      is_roll   = ic%tang.eq.2 .or. ic%tang.eq.3

      if (ic%rztang.eq.0) then
         call gf3_new(exrhs, 'geom%exrhs', cgrid)
         call gf3_set(AllElm, 0d0, exrhs, ikTANG)
      endif

      ! M=3, using FASTSIM approach with 3 flexibilities: scale spin creepage by factors L1/L3 and L2/L3

      if (ic%mater.eq.3) then
         facx = mater%flx(1) / mater%flx(3)
         facy = mater%flx(2) / mater%flx(3)
      else
         facx = 1d0
         facy = 1d0
      endif

      ! compute the rigid slip (wx,wy) of body (1) w.r.t. body (2) and place it in the tangential
      !    right-hand side as hs1 = -dq w.

      if (is_roll) then

         ! rolling:

         cc = cos(chi)
         sc = sin(chi)
         ! write(bufout,'(a,f6.3)') ' using f=',f
         ! call write_log(1, bufout)

         do ii = 1, cgrid%ntot
            hs1%vx(ii) = facx * cphi *(-cgrid%y(ii) + spinyo - sc*dq*facphi) + exrhs%vx(ii) ! below: * -dq 
            hs1%vy(ii) = facy * cphi *( cgrid%x(ii) - spinxo + cc*dq*facphi) + exrhs%vy(ii) !        * -dq
         enddo
      else

         ! shift:

         do ii = 1, cgrid%ntot
            hs1%vx(ii) = facx * cphi *(-cgrid%y(ii) + spinyo) + exrhs%vx(ii)
            hs1%vy(ii) = facy * cphi *( cgrid%x(ii) - spinxo) + exrhs%vy(ii)
         enddo
      endif

      ! if x-creepage prescribed:

      if (ic%force3.eq.0) then
         do ii = 1, cgrid%ntot
            hs1%vx(ii) = hs1%vx(ii) + cksi
         enddo
      endif

      ! if y-creepage prescribed:

      if (ic%force3.le.1) then
         do ii = 1, cgrid%ntot
            hs1%vy(ii) = hs1%vy(ii) + ceta
         enddo
      endif

      ! change sign and multiply with dq (=1 in case of shifts)

      call gf3_scal(AllElm, -dq, hs1, ikTANG)

      end associate

   end subroutine set_tang_rhs

!------------------------------------------------------------------------------------------------------------

   subroutine init_curr_data(ic, hz, potcon, cgrid, geom, mater, kin, outpt, ihertz)
!--purpose: initialize contact area, tractions and other arrays for the new case
!            2) all cases: initialise s-arrays for new case
      implicit none
!--subroutine parameters:
      type(t_ic)              :: ic
      type(t_hertz)           :: hz
      type(t_potcon)          :: potcon
      type(t_grid)            :: cgrid
      type(t_geomet)          :: geom
      type(t_material)        :: mater
      type(t_kincns)          :: kin
      type(t_output)          :: outpt
      integer                 :: ihertz
!--local variables:
      integer      :: ii
      real(kind=8) :: pabs, pnmax, pnfac
      logical      :: is_roll, use_iestim_hertz

      associate(igs   => outpt%igs,   mus   => outpt%mus,  ps    => outpt%ps,    us    => outpt%us,     &
                ss    => outpt%ss,    upls  => outpt%upls, taucs => outpt%taucs, temp1 => outpt%temp1,  &
                temp2 => outpt%temp2, sens  => outpt%sens, aa    => hz%aa,       bb    => hz%bb)
      call timer_start(itimer_filpvs)

      is_roll   = ic%tang.eq.2 .or. ic%tang.eq.3

      !------------------------------------------------------------------------------------------------------
      ! Phase 2) in all cases: initialise s-arrays for new case
      !------------------------------------------------------------------------------------------------------
      ! Set the initial traction ps, according to iestim :
      !   iestim=0 : initiate contact, init.est. according to undeformed distance Hs
      !   iestim=1 : init.est. according to igs, regularize ps
      !   iestim=2 : init.est. according to C,E,  set ps t = 0, i in Adhes
      !   iestim=3 : init.est. according to H,S,P,E  leave igs, ps untouched

      ! clear/initialize all displacement differences

      call gf3_set(AllElm, 0d0, us, ikALL)
      call gf3_set(AllElm, 0d0, upls, ikTANG)

      ! reset normal tractions when I=0

      if (ic%iestim.eq.0) then
         call gf3_set(AllElm, 0d0, ps, ikZDIR)
      endif

      ! reset tangential tractions when I=0 or I=2, or when T=0

      if (ic%iestim.eq.0 .or. ic%iestim.eq.2 .or. ic%tang.eq.0) then
         call gf3_set(AllElm, 0d0, ps, ikTANG)
         call gf3_set(AllElm, mater%tau_c0, taucs, ikTANG)
      endif

      ! set elliptical traction bound when requested (B=2) or
      !    when Hertzian solution may be used as initial estimate (I=0)

      use_iestim_hertz = (ic%iestim.eq.0 .and. ic%bound.eq.0 .and. ihertz.ge.1)

      if (ic%bound.eq.2 .or. use_iestim_hertz) then

         if (ihertz.le.3) then

            ! 3D elliptical traction bound

            pnmax = 3d0 * kin%fntrue / (2d0 * pi * aa * bb)
            do ii = 1, cgrid%ntot
               pnfac = 1d0 - (cgrid%x(ii)/aa)**2 - (cgrid%y(ii)/bb)**2
               ps%vn(ii) = pnmax * sqrt(max(0d0, pnfac))
            enddo
         else

            ! 2D elliptical traction bound, uniform in y

            pnmax = kin%fntrue / (pi * aa * bb)
            do ii = 1, cgrid%ntot
               pnfac = 1d0 - (cgrid%x(ii)/aa)**2
               if (abs(cgrid%y(ii)).gt.bb) pnfac = 0d0
               ps%vn(ii) = pnmax * sqrt(max(0d0, pnfac))
            enddo
         endif
      endif

      ! set parabolic traction bound when requested (B=3)

      if (ic%bound.eq.3 .and. ihertz.le.3) then

         ! 3D parabolical traction bound

         pnmax = 2d0 * kin%fntrue / (pi * aa * bb)
         do ii = 1, cgrid%ntot
            pnfac = 1d0 - (cgrid%x(ii)/aa)**2 - (cgrid%y(ii)/bb)**2
            ps%vn(ii) = pnmax * max(0d0, pnfac)
         enddo

      elseif (ic%bound.eq.3) then

         ! 2D parabolical traction bound, uniform in y

         pnmax = 3d0 * kin%fntrue / (8d0 * aa * bb)
         do ii = 1, cgrid%ntot
            pnfac = 1d0 - (cgrid%x(ii)/aa)**2
            if (abs(cgrid%y(ii)).gt.bb) pnfac = 0d0
            ps%vn(ii) = pnmax * max(0d0, pnfac)
         enddo

      endif

      ! set SDEC pressures when requested (B=4)

      if (ic%bound.eq.4) call sdec_pn(aa, hz%bpos, hz%bneg, kin%fntrue, cgrid, ps)

      ! set normal part of element division: C <--> E

      if ((ic%bound.ge.2 .and. ic%bound.le.4) .or. use_iestim_hertz) then

         ! if Hertzian (elliptical, parabolic) or SDEC solution is used:
         !    adjust element division C,E accordingly

         do ii = 1, cgrid%ntot
            if (ps%vn(ii).gt.1e-20) then
               igs%el(ii) = Adhes
            else
               igs%el(ii) = Exter
            endif
         enddo

      elseif (ic%iestim.eq.0) then

         ! I=0, Otherwise: use Eldiv0 to set initial element division
         !      Note: for N=1 this also estimates Pen and d Fntrue/d Pen.

         call eldiv0(ic, potcon, cgrid, geom, mater, kin, outpt)

      endif

      ! regularize tractions Pn, Pt when requested (I=1)

      if (ic%iestim.eq.1) then
         do ii = 1, cgrid%ntot
            if (igs%el(ii).le.Exter) then
               ps%vn(ii) = 0d0
               ps%vx(ii) = 0d0
               ps%vy(ii) = 0d0
            elseif (igs%el(ii).eq.Slip) then
               pabs= max(1d-10, sqrt(ps%vx(ii)**2 + ps%vy(ii)**2))
               ps%vx(ii) = mus%vt(ii) * ps%vn(ii) * ps%vx(ii) / pabs
               ps%vy(ii) = mus%vt(ii) * ps%vn(ii) * ps%vy(ii) / pabs
            endif
         enddo
      endif

      ! adjust tangential element division H <--> S, P when required (I=2)

      if (ic%iestim.eq.2) then
         do ii = 1, cgrid%ntot
            if (igs%el(ii).ge.Adhes) igs%el(ii) = Adhes
         enddo
      endif

      ! filter elements outside the planform:

      do ii = 1, cgrid%ntot
         if (geom%hs1%vn(ii).gt.1e29) igs%el(ii) = Exter
      enddo
      call areas(igs)

      ! initial estimate of sensitivity to penetration when I=0:

      if (ic%iestim.eq.0) sens(iout_fx1:iout_mz1, iin_dpen) = 0d0
      if (ic%iestim.eq.0 .and. (ic%norm.ne.1 .or. ihertz.ne.0)) sens(iout_fn, iin_dpen) = 0d0

      ! initial estimates of creepages and tangential sensitivities:

      if (ic%iestim.eq.0 .or. ic%iestim.eq.2) then
         if (ic%force3.ge.1) kin%cksi = 1d-6
         if (ic%force3.eq.2) kin%ceta = 0d0
         sens(iout_fn:iout_mz1, iin_dksi1:iin_dphi1) = 0d0
      endif

      call timer_stop(itimer_filpvs)
      end associate
   end subroutine init_curr_data

!------------------------------------------------------------------------------------------------------------

   subroutine sdec_pn(aa, bpos, bneg, fntrue, cgrid, ps)
!--purpose: set the pressure pn according the SDEC parameters
      implicit none
!--subroutine arguments :
      real(kind=8)            :: aa, bpos, bneg, fntrue
      type(t_grid)            :: cgrid
      type(t_gridfnc3)        :: ps
!--local variables :
      integer                   :: ix, iy, ii
      real(kind=8)              :: bb, pnmax, psi, xl, xx, y0, yy

      ! bb = (bpos + bneg) / 2,  y0 = -(bpos - bneg) / 2,  psi = -y0 / b
      ! psi > 0 when bpos > bneg, in which case y0 < 0

      bb    = (bpos + bneg) / 2d0
      psi   = (bpos - bneg) / (2d0*bb)
      y0    = -bb * psi
      pnmax = 1.5d0 * fntrue / (pi * aa * bb)

      do iy = 1, cgrid%ny
         ii = 1 + (iy-1) * cgrid%nx
         yy = cgrid%y(ii)
        
         if (yy.gt.y0) then
            ! upper half ellipse: x_l,pos
            xl = aa * sqrt( 1d0 - (yy-y0)**2 / (bb*(1d0+psi))**2 )
         else
            ! lower half ellipse: x_l,neg
            xl = aa * sqrt( 1d0 - (yy-y0)**2 / (bb*(1d0-psi))**2 )
         endif

         do ix = 1, cgrid%nx
            ii = ix + (iy-1) * cgrid%nx
            xx = cgrid%x(ii)
            if (abs(xx).lt.xl) then
               ps%vn(ii) = pnmax/aa * sqrt(xl**2 - xx**2)
            else
               ps%vn(ii) = 0d0
            endif
         enddo
      enddo

   end subroutine sdec_pn

!------------------------------------------------------------------------------------------------------------

   subroutine eldiv0(ic, potcon, cgrid, geom, mater, kin, outpt)
!--purpose: Guess the element division C <--> E for the normal problem.
!           Note: for N=1 this also estimates Pen and d Fntrue/d Pen.
      implicit none
!--subroutine arguments :
      type(t_ic)              :: ic
      type(t_potcon)          :: potcon
      type(t_grid)            :: cgrid
      type(t_geomet)          :: geom
      type(t_material)        :: mater
      type(t_kincns)          :: kin
      type(t_output)          :: outpt
!--local variables :
      integer,      parameter   :: maxit = 50, idebug = 0
      real(kind=8), parameter   :: facpen = 0.60, reltol = 0.01
      integer                   :: ii, iter, ncnmin, ncnmid, ncnmax
      real(kind=8)              :: hsmin, pentru, rm, cdy, rmn, fac, penmin, penmid, penmax,            &
                                   fnmin, fnmid, fnmax, fnscal

      associate(igs  => outpt%igs, hs   => geom%hs1)

      fnscal = kin%fntrue / mater%ga

      ! 1) Determine the minimum of the undeformed distance in Hs.

      !    The "Undeformed profile height" becomes Hs - min(Hs) >= 0.
      !    The "True penetration" becomes Pen - min(Hs) > 0.
      !    The "Undeformed distance" stays Hs - Pen = Hs - hsmin - pentru

      ii     = idmin (cgrid%ntot, hs%vn, 1)
      hsmin  = hs%vn(ii)
      pentru = kin%pen - hsmin
      if (idebug.ge.3) then
         write(bufout,*) 'Minimum Hs=',hsmin
         call write_log(1, bufout)
      endif

      ! 2) N=1, total force Fn prescribed: Determine estimate for pentru.
      !    See explanation in report-m09036.

      if (ic%norm.eq.1 .and. cgrid%ny.eq.1) then

      ! Case A: 2D problem, MY=1.

         ! The relation between true penetration and force per unit width was fitted as:
         !      C = 0.35 - 0.05 * (log(dy) - log(200))
         !      Fn/G/dy = C * 0.6/(1-nu) * pen^1.08 * R^0.1
         !    <-->
         !      pen = { Fn/G/dy * (1-nu)/0.6 / C / R^0.1 }^0.926

         if (geom%ibase.eq.1) then
            rm = 0.5d0 / max(1d-6, geom%prmudf(1))
         elseif (geom%ibase.eq.2) then
            rm = geom%prmudf(3)
         elseif (geom%ibase.eq.3) then
            rm = 0.5d0 * geom%prmudf(7)
         else
            rm = 1d0
         endif
         cdy = 0.35d0 - 0.05d0 * (log(potcon%dy)-log(200d0))

         pentru = ( fnscal/potcon%dy * (1d0-mater%nu)/0.6d0 / cdy / rm**0.1 )**0.926
         kin%pen = pentru + hsmin

         ! Set estimate of sensitivity dfntrue/dpen

         outpt%sens(iout_fn, iin_dpen) = potcon%dy * cdy * mater%ga *                             &
                  0.6d0/(1d0-mater%nu) * rm**0.1 * 1.08d0 * pentru**0.08
         if (idebug.ge.2) then
            write(*,'(4(a,f7.3),a)') ' A: pen=(',fnscal/potcon%dy,                                &
              ' *',(1d0-mater%nu)/0.6d0,' /',cdy,' /',rm**0.1,')**0.926'
            write(*,*) 'A: pen=',kin%pen,' (hsmin=',hsmin,')'
            write(*,*) 'A: df/dpen=',outpt%sens(iout_fn, iin_dpen)
         endif

      elseif (ic%norm.eq.1) then

      ! Case B: 3D problem, MY>1.

         ! Assume that the deformation is <= facpen times the true penetration
         ! Assume that pentru*Area = 1 * 0.75*sqrt(pi) * (1-nu) * Fn/G
         ! Use a bi-section algorithm to determine pentru that delivers
         ! the requested Fn/G.

         ! set constants, rmn == r * sqrt(m*n) of Hertzian formulae

         rmn    = 1d0
         fac    = rmn * 0.75 * sqrt(pi) * (1d0 - mater%nu)

         ! set lower and upper bound for true penetration, determine
         ! corresponding Ncon

         penmin = 0d0
         ncnmin = 1
         fnmin  = 0d0

         penmax = fac * fnscal / sqrt(potcon%dxdy)
         ncnmax = 0
         do ii = 1, cgrid%ntot
            if (hs%vn(ii)-hsmin.lt.facpen*penmax) ncnmax = ncnmax + 1
         enddo
         fnmax  = penmax * sqrt(potcon%dxdy * ncnmax) / fac

         ! perform bi-section algorithm until Fn is approximated within
         ! 1% accuracy

         iter = 0
         do while ( iter.lt.maxit .and. abs(fnmax-fnmin).gt.reltol*fnscal )
            iter = iter + 1
            if (idebug.ge.2) then
               write(bufout,101) iter, penmin, penmax, fnmin, fnmax
               call write_log(1, bufout)
 101           format(' iter=',i3,': Pen in [',f8.5,',',f8.5, ' ], Fn in [',f7.4,',',f7.4,' ]')
            endif

            ! set intermediate point, determine corresponding Ncon

            penmid = 0.5d0 * (penmax + penmin)
            ncnmid = 0
            do ii = 1, cgrid%ntot
               if (hs%vn(ii)-hsmin.lt.facpen*penmid) ncnmid=ncnmid+1
            enddo

            ! estimate corresponding total force

            fnmid = penmid * sqrt(potcon%dxdy * ncnmid) / fac
            if (idebug.ge.2) then
               write(bufout,102) penmid, ncnmid, fnmid
               call write_log(1, bufout)
 102           format(18x,'mid= ',f8.5,' : ncon=',i7,', fn=', f7.4)
            endif

            ! update either lower or upper bound for pentru, depending on
            ! whether fnmid is too low or too high

            if (fnmid.lt.fnscal) then

               ! fnmid too low, pentru must be >= penmid

               penmin = penmid
               ncnmin = ncnmid
               fnmin  = fnmid
            else

               ! fnmid too high, pentru must be <= penmid

               penmax = penmid
               ncnmax = ncnmid
               fnmax  = fnmid
            endif
         enddo ! end-while (bisection algorithm)

         ! Set estimate of true penetration

         pentru  = penmax
         kin%pen = pentru + hsmin


         ! Set estimate of sensitivity dfntrue / dpen

         if (ncnmin.lt.ncnmax) then
            outpt%sens(iout_fn, iin_dpen) = mater%ga * (fnmax - fnmin) / max(1d-6,penmax-penmin)
         else
            outpt%sens(iout_fn, iin_dpen) = mater%ga * sqrt(potcon%dxdy * ncnmax) / fac
         endif

         if (idebug.ge.1) then
            write(bufout,111) iter, kin%pen, ncnmax, fnmax, penmax-penmin, fnmax-fnmin,         &
                outpt%sens(iout_fn,iin_dpen)
            call write_log(2, bufout)
 111        format(' Bi-section: niter=',i4,', Pen=',f8.5,', Ncon=',i7, ', Fn=',f7.4,/,         &
                        13x,'dPen=',f8.5,', dFn=',f7.4, ', sens=',f9.3)
         endif
      endif

      ! 3) Set initial guess for contact area, penetration prescribed:

      !    Assume that the deformation is <= facpen times the true penetration
      !    Contact area == all elements with undeformed profile height less
      !    than facpen times the true penetration

      do ii = 1, cgrid%ntot
         if (hs%vn(ii)-hsmin.lt.facpen*pentru) then
            igs%el(ii) = Adhes
         else
            igs%el(ii) = Exter
         endif
      enddo
      end associate
   end subroutine eldiv0

!------------------------------------------------------------------------------------------------------------

end module m_sdis

