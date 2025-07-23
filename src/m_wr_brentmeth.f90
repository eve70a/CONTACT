!------------------------------------------------------------------------------------------------------------
! m_wr_brentmeth - building blocks for Brent / Broyden methods for cases with total forces prescribed
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_wr_brentmeth

use m_wrprof_data
use m_wr_solvecp

implicit none
private

   ! Data type for Brent iterates:

   public t_brent_its

   ! Building blocks for Brent method:

   public  brent_its_init
   public  brent_its_find_k
   public  brent_its_update_ptr
   public  brent_its_add_iterate
   public  brent_sensitivity_k
   public  brent_its_print
   public  brent_its_destroy
   public  brent_its_has_bracket
   public  brent_its_has_jump
   public  brent_its_has_chng_stiff
   public  brent_set_xnew

   public   IMETH_BRENT, IMETH_BROYDN, IMETH_SECANT
   integer, parameter :: IMETH_BRENT = 1, IMETH_BROYDN = 2, IMETH_SECANT = 3

!------------------------------------------------------------------------------------------------------------
!--data to keep track of iterates in Brent's algorithm
!       Brent's method keeps four iterates [a, b, c, d] with function values [fa, fb, fc, fd]
!           b is the best guess, a is neighbour of b with opposite residual, c == b^{k-1}, d == b^{k-2}
!       we additionally use [x0, x1] for the bracket, and we store all iterates in itdata (sorted)

   type :: t_brent_one_it
      integer      :: k, np
      real(kind=8) :: xk, rk

      ! k           iteration number
      ! xk          x-value of iteration
      ! rk          residual value of iteration
      ! np          number of contact patches
   end type t_brent_one_it

   type :: t_brent_its
      integer                          :: ikarg, maxit, numit
      real(kind=8)                     :: ftarg
      type(t_brent_one_it), pointer    :: it_0 => NULL(), it_1 => NULL(), it_a => NULL(),               &
                                          it_b => NULL(), it_c => NULL(), it_d => NULL(),               &
                                          it_k => NULL(), it_km1 => NULL(), it_km2 => NULL()
      type(t_brent_one_it), dimension(:), allocatable :: itdata

      ! ikarg       meta-information on type of problem: 
      !                -ikYDIR for solving Fy(y_shift), left side,
      !                 ikYDIR for solving Fy(y_shift), right side,
      !                 ikZDIR when solving Fz(z_ws),
      !                -ikZDIR when solving Ftot(z_shift)
      ! maxit       size of arrays used, max. #iterates
      ! numit       actual #iterates stored
      ! ftarg       target force value
      ! itdata      per-iteration data sorted on x in ascending order
      ! it_0        lower side of current bracket, NULL as long as no bracket is found
      ! it_1        upper side of current bracket, NULL as long as no bracket is found
      ! it_a        contra-point to current best guess (other end of bracket)
      ! it_b        current best guess
      ! it_c        previous best guess
      ! it_d        previous previous best guess
      ! it_k        most recent iterate x^k
      ! it_km1      previous iterate x^{k-1}
      ! it_km2      previous iterate x^{k-2}
   end type t_brent_its

contains

!------------------------------------------------------------------------------------------------------------

   subroutine brent_its_init(its, ikarg, maxit, ftarg)
!--purpose: initialize the iterates structure for Brent's algorithm
      implicit none
!--subroutine arguments:
      type(t_brent_its), target :: its
      integer                   :: ikarg, maxit
      real(kind=8)              :: ftarg

      its%ikarg = ikarg
      its%maxit = max(1, min(10000, maxit))
      its%numit = 0
      its%ftarg = ftarg
      allocate(its%itdata(its%maxit))

      call brent_its_update_ptr(its)
   end subroutine brent_its_init

!------------------------------------------------------------------------------------------------------------

   subroutine brent_its_destroy(its)
!--purpose: cleanup iterates structure for Brent's algorithm
      implicit none
!--subroutine arguments:
      type(t_brent_its), target :: its

      deallocate(its%itdata)
   end subroutine brent_its_destroy

!------------------------------------------------------------------------------------------------------------

   function brent_its_find_k(its, k)
!--purpose: locate sequence number it for iteration k
      implicit none
!--function result:
      integer                   :: brent_its_find_k
!--subroutine arguments:
      type(t_brent_its)         :: its
      integer                   :: k
!--local variables:
      integer                   :: it

      brent_its_find_k = -1
      if (k.ge.0 .and. k.le.its%numit-1) then
         do it = 1, its%numit
            if (its%itdata(it)%k.eq.k) brent_its_find_k = it
         enddo
      endif

   end function brent_its_find_k

!------------------------------------------------------------------------------------------------------------

   subroutine brent_its_update_ptr(its)
!--purpose: update the pointers it_0--it_1 (bracket), it_a--it_d for current and previous best points,
!                               it_k, it_km1, it_km2 (current/previous iterates)
      implicit none
!--subroutine arguments:
      type(t_brent_its), target :: its
!--local variables
      integer      :: it, it_br0, imin
      real(kind=8) :: rmin

      if (its%numit.le.0) then

         its%it_0   => NULL()   ! current bracket
         its%it_1   => NULL()

         its%it_a   => NULL()   ! points a, b, c, d
         its%it_b   => NULL()
         its%it_c   => NULL()
         its%it_d   => NULL()

         its%it_k   => NULL()   ! most recent iterates
         its%it_km1 => NULL()
         its%it_km2 => NULL()

      elseif (its%numit.le.1) then

         its%it_0   => NULL()
         its%it_1   => NULL()

         its%it_a   => its%itdata(1)
         its%it_b   => its%itdata(1)
         its%it_c   => its%itdata(1)
         its%it_d   => its%itdata(1)

         its%it_k   => its%itdata(1)
         its%it_km1 => its%itdata(1)
         its%it_km2 => its%itdata(1)

      else

         ! cycle previous best iterates

         it = brent_its_find_k(its, its%it_c%k)         ! Note: using it_c
         its%it_d   => its%itdata(it)                   ! Note: point to itdata() instead of it_c
         it = brent_its_find_k(its, its%it_b%k)
         its%it_c   => its%itdata(it)

         ! determine current best iterate

         if (brent_its_has_bracket(its, it_br0)) then
            if (abs(its%itdata(it_br0)%rk).lt.abs(its%itdata(it_br0+1)%rk)) then
               imin = it_br0
            else
               imin = it_br0 + 1
            endif
         else
            imin = 0
            rmin = 1d99
            do it = 1, its%numit
               if (abs(its%itdata(it)%rk).lt.rmin) then
                  imin = it
                  rmin = abs(its%itdata(it)%rk)
               endif
            enddo
         endif

         its%it_b   => its%itdata(imin)

         ! determine contra-point of current bracket (no bracket: sensible adjacent point)

         if (imin.le.1) then
            its%it_a   => its%itdata(imin+1)
         elseif (imin.ge.its%numit) then
            its%it_a   => its%itdata(imin-1)
         elseif (its%itdata(imin-1)%rk * its%it_b%rk.lt.0d0) then
            its%it_a   => its%itdata(imin-1)
         else
            its%it_a   => its%itdata(imin+1)
         endif

         ! set pointers to bracket

         if (.not.brent_its_has_bracket(its, it_br0)) then
            its%it_0 => NULL()
            its%it_1 => NULL()
         elseif (its%it_a%xk.lt.its%it_b%xk) then
            its%it_0 => its%it_a
            its%it_1 => its%it_b
         else
            its%it_0 => its%it_b
            its%it_1 => its%it_a
         endif

         ! set pointers to most recent iterates

         it = brent_its_find_k(its,       its%numit-1)
         its%it_k   => its%itdata(it)
         it = brent_its_find_k(its,       its%numit-2)
         its%it_km1 => its%itdata(it)
         it = brent_its_find_k(its, max(0,its%numit-3))
         its%it_km2 => its%itdata(it)

      endif
      ! call write_log(' brent_its_update_ptr ok...')

   end subroutine brent_its_update_ptr

!------------------------------------------------------------------------------------------------------------

   subroutine brent_its_add_iterate(its, k, npatch, xk, rk, my_ierror)
!--purpose: insert an iterate in the sorted structure for Brent's algorithm
      implicit none
!--subroutine arguments:
      type(t_brent_its), target :: its
      integer                   :: k, npatch, my_ierror
      real(kind=8)              :: xk, rk
!--local variables
      integer           :: i, j
      logical           :: ldone

      my_ierror = 0
      if (isnan(rk)) then
         my_ierror = -1
         write(bufout,'(a,f11.3,a)') ' ...NaN-values found (', rk, '), aborting'
         call write_log(1, bufout)
      endif

      ! determine first position i with xk < its%xk(i)

      i     = 0
      ldone = (i.ge.its%numit)
      do while (.not.ldone)
         i     = i + 1
         ldone = .true.
         if (i.le.its%numit) ldone = (xk.lt.its%itdata(i)%xk)
      enddo

      ! shift iterates [i--end] one position

      do j = its%numit, max(i,1), -1
         its%itdata(j+1) = its%itdata(j)
         ! update pointers it_b, it_c, used in update_ptr
         if (its%it_b%k.eq.its%itdata(j)%k) its%it_b => its%itdata(j+1)
         if (its%it_c%k.eq.its%itdata(j)%k) its%it_c => its%itdata(j+1)
      enddo

      ! insert iterate at position i

      if (i.le.0) i = 1

      its%numit        = its%numit + 1
      its%itdata(i)%k  = k
      its%itdata(i)%np = npatch
      its%itdata(i)%xk = xk
      its%itdata(i)%rk = rk

      ! update pointers to current and previous best estimates

      call brent_its_update_ptr(its)

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
         if (its%itdata(j)%k.eq.k) it = j
      enddo

      if (it.eq.-1) then
         dfdx = -99d9
      elseif (its%numit.le.1) then
         dfdx = 0d0
      elseif (it.le.1) then
         ! forward difference
         dfdx = (its%itdata(it+1)%rk - its%itdata(it)%rk) /                                             &
                                                max(tiny, its%itdata(it+1)%xk - its%itdata(it)%xk)
      elseif (it.ge.its%numit) then
         ! backward difference
         dfdx = (its%itdata(it)%rk - its%itdata(it-1)%rk) /                                             &
                                                max(tiny, its%itdata(it)%xk - its%itdata(it-1)%xk)
         ! initial estimate for dFz/dpen: Fz ~ pen^1.5
         if (its%ikarg.eq.ikZDIR .and. it.eq.2) dfdx = 1.5d0 * dfdx
      elseif (.false.) then
         ! central difference
         dfdx = (its%itdata(it+1)%rk - its%itdata(it-1)%rk) /                                           &
                                                max(tiny, its%itdata(it+1)%xk - its%itdata(it-1)%xk)
      elseif (abs(its%itdata(it+1)%rk).lt.abs(its%itdata(it-1)%rk)) then
         ! forward difference
         dfdx = (its%itdata(it+1)%rk - its%itdata(it)%rk) /                                             &
                                                max(tiny, its%itdata(it+1)%xk - its%itdata(it)%xk)
      else
         ! backward difference
         dfdx = (its%itdata(it)%rk - its%itdata(it-1)%rk) /                                             &
                                                max(tiny, its%itdata(it)%xk - its%itdata(it-1)%xk)
      endif
      brent_sensitivity_k = dfdx

      if (idebug_br.ge.6) then
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
      integer                  :: it, j, ilen
      real(kind=8)             :: dfdx
      character(len=12)        :: str12(6)
      character(len=16)        :: str16(6)
      character(len=24)        :: strptr

      if (ic%flow.ge.1 .and. idebug_br.le.1) then

         ! print one line with current iterate / residual

         it   = brent_its_find_k(its, k)
         dfdx = brent_sensitivity_k(k, its, idebug_br)

         write(str12(1), '(f12.4)') its%itdata(it)%xk
         str12(2) = fmt_gs(12, 4, 4, its%itdata(it)%rk+its%ftarg)
         str12(3) = fmt_gs(12, 4, 4, dfdx)
         if (abs(its%ikarg).eq.ikYDIR) then
            write(bufout,1020) its%itdata(it)%k, (str12(j),j=1,3), its%itdata(it)%np
         elseif (its%ikarg.eq. ikZDIR) then
            write(bufout,1030) its%itdata(it)%k, (str12(j),j=1,3)
         elseif (its%ikarg.eq.-ikZDIR) then
            write(bufout,1040) its%itdata(it)%k, (str12(j),j=1,3)
         else
            call write_log(' Internal error (Brent): invalid ikarg')
            call abort_run()
         endif
         call write_log(1, bufout)

 1020    format(2x, i4,', BR,  yshift, Fytot: ',2a,', dFy/dy:',a,',',i3,' patches')
 1030    format(4x, i6,', NR,  z_ws, Fz: ',2a,', dFz/dz:',a)
 1040    format(4x, i6,', NR,  zshift, Fztot: ',2a,', dFz/dz:',a)

      endif

      if (ic%flow.ge.1 .and. idebug_br.eq.2) then

         ! print one line with current iterate / residual

         it   = brent_its_find_k(its, k)
         dfdx = brent_sensitivity_k(k, its, idebug_br)

         write(str12(1), '(f12.4)') its%itdata(it)%xk
         str12(2) = fmt_gs(12, 4, 4, its%itdata(it)%rk+its%ftarg)
         str12(3) = fmt_gs(12, 4, 4, its%itdata(it)%rk)
         str12(4) = fmt_gs(12, 4, 4, dfdx)

         if (abs(its%ikarg).eq.ikYDIR) then
            write(bufout,2020) its%itdata(it)%k, (str12(j),j=1,4), its%itdata(it)%np
         elseif (its%ikarg.eq. ikZDIR) then
            write(bufout,2030) its%itdata(it)%k, (str12(j),j=1,4)
         elseif (its%ikarg.eq.-ikZDIR) then
            write(bufout,2040) its%itdata(it)%k, (str12(j),j=1,4)
         endif
         call write_log(1, bufout)
      endif

      if (idebug_br.ge.3) then

         ! print table with all iterates and residuals

         write(bufout,'(a,i3,2a)') ' --- brent table for',its%numit,' iterations. B=best, A=contra, ',  &
                'C,D=prev best, bracket [x0,x1] --------------------'
         call write_log(1, bufout)

         do it = 1, its%numit
            dfdx = brent_sensitivity_k(its%itdata(it)%k, its, idebug_br)

            strptr = ' <--'
            if (associated(its%it_0)) then
               if (its%itdata(it)%k.eq.its%it_0%k) strptr = trim(strptr) // ' x0,'
            endif
            if (associated(its%it_1)) then
               if (its%itdata(it)%k.eq.its%it_1%k) strptr = trim(strptr) // ' x1,'
            endif
            if (its%itdata(it)%k.eq.its%it_a%k) strptr = trim(strptr) // ' A,'
            if (its%itdata(it)%k.eq.its%it_b%k) strptr = trim(strptr) // ' B,'
            if (its%itdata(it)%k.eq.its%it_c%k) strptr = trim(strptr) // ' C,'
            if (its%itdata(it)%k.eq.its%it_d%k) strptr = trim(strptr) // ' D,'
            if (idebug_br.le.1 .or. len(trim(strptr)).le.4) strptr = ' '
            ilen = max(1, len(trim(strptr)))

            if (idebug_br.ge.6 .or. .true.) then
               str16(1) = fmt_gs(16, 8, 8, its%itdata(it)%xk)
               str16(2) = fmt_gs(16, 8, 8, its%itdata(it)%rk+its%ftarg)
               str16(3) = fmt_gs(16, 8, 8, its%itdata(it)%rk)
               str16(4) = fmt_gs(16, 8, 8, dfdx)
               if (abs(its%ikarg).eq.ikYDIR) then
                  write(bufout,2020) its%itdata(it)%k, (str16(j),j=1,4), its%itdata(it)%np, strptr(1:ilen-1)
               elseif (its%ikarg.eq. ikZDIR) then
                  write(bufout,2030) its%itdata(it)%k, (str16(j),j=1,4), its%itdata(it)%np, strptr(1:ilen-1)
               elseif (its%ikarg.eq.-ikZDIR) then
                  write(bufout,2040) its%itdata(it)%k, (str16(j),j=1,4), strptr(1:ilen-1)
               endif
            else
               write(str12(1), '(f12.4)') its%itdata(it)%xk
               str12(2) = fmt_gs(12, 4, 4, its%itdata(it)%rk+its%ftarg)
               str12(3) = fmt_gs(12, 4, 4, its%itdata(it)%rk)
               str12(4) = fmt_gs(12, 4, 4, dfdx)
               if (abs(its%ikarg).eq.ikYDIR) then
                  write(bufout,2020) its%itdata(it)%k, (str12(j),j=1,4), its%itdata(it)%np, strptr(1:ilen-1)
               elseif (its%ikarg.eq. ikZDIR) then
                  write(bufout,2030) its%itdata(it)%k, (str12(j),j=1,4), its%itdata(it)%np, strptr(1:ilen-1)
               elseif (its%ikarg.eq.-ikZDIR) then
                  write(bufout,2040) its%itdata(it)%k, (str12(j),j=1,4), strptr(1:ilen-1)
               endif
            endif
            call write_log(1, bufout)
         enddo

 2020    format(2x, i4,', BR,  yshift, Fytot, res:',3a,', dFy/dy:',a,',',i3,' patches',a)
 2030    format(4x, i6,', NR,  z_ws, Fz:',2a,', res: ',a,', dFz/dz:',a,',',i3,' patches',a)
 2040    format(4x, i6,', NR,  zshift, Fztot:',2a,', res: ',a,', dFz/dz:',2a)

         call write_log(' --- end of brent iteration table -----------------------' //                  &
                                                   '------------------------------------------------')
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

      if (its%numit.le.1) then
         brent_its_has_bracket = .false.
         it = 0
      else
         it = 1
         do while(it.lt.its%numit .and. its%itdata(it)%rk*its%itdata(it+1)%rk.ge.0d0)
            it = it + 1
         enddo
         brent_its_has_bracket = (it.lt.its%numit)
      endif
      if (present(it0)) it0 = it
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
         it0 = brent_its_find_k(its, its%it_0%k)
         it1 = brent_its_find_k(its, its%it_1%k)

         dx_brack   = (its%it_1%xk - its%it_0%xk)
         dfdx_brack = (its%it_1%rk - its%it_0%rk) / dx_brack

         if (it0.gt.1) then
            ! backward difference
            dfdx_left = (its%it_0%rk - its%itdata(it0-1)%rk) / (its%it_0%xk - its%itdata(it0-1)%xk)
         else
            dfdx_left = 1d20
         endif
         if (it1.lt.its%numit) then
            ! forward difference
            dfdx_right = (its%itdata(it1+1)%rk - its%it_1%rk) / (its%itdata(it1+1)%xk - its%it_1%xk)
         else
            dfdx_right = 1d20
         endif

         ! jump == dfdx in bracket >> dfdx on both sides

         brent_its_has_jump =     dx_brack .lt. jump_thrs_dx .and.                                      &
                              abs(dfdx_brack) .gt. jump_thrs_df*max(abs(dfdx_left), abs(dfdx_right))

         if (brent_its_has_jump .and. idebug_br.ge.-1) then
            write(bufout,'(2(a,f12.6),a)') '  ...detected jump in x=[', its%it_0%xk,',', its%it_1%xk,']'
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
            x0_new  = its%itdata(1)%xk
         else
            x0_new  = 0.5d0 * (its%itdata(1)%xk + its%itdata(2)%xk)
         endif

      else

         ! determine left/right sides of current bracket

         it0 = brent_its_find_k(its, its%it_0%k)
         it1 = brent_its_find_k(its, its%it_1%k)

         if (it0.le.1) then

            ! bracket is first segment of list: report no stiffness change

            brent_its_has_chng_stiff = .false.

            ! dummy output for x_new: bisection

            x0_new    = 0.5d0 * (its%it_0%xk + its%it_1%xk)

         elseif (it1.ge.its%numit) then

            ! bracket is last segment of list: compare left to current stiffness

            dfdx_left  = (its%it_0%rk - its%itdata(it0-1)%rk) / (its%it_0%xk - its%itdata(it0-1)%xk)
            dfdx_brack = (its%it_1%rk - its%it_0%rk)          / (its%it_1%xk - its%it_0%xk)
            dfdx_right = 0d0

            brent_its_has_chng_stiff = (dfdx_brack .gt. stiff_thrs*dfdx_left)

            if (brent_its_has_chng_stiff) then
               x0_new = (1d0-frac_last) * its%it_0%xk + frac_last * its%it_1%xk
            else
               x0_new = 0.5d0 * (its%it_0%xk + its%it_1%xk)
            endif

         else

            ! bracket has left and right segments: compare left to right stiffness

            dfdx_left  = (its%it_0%rk - its%itdata(it0-1)%rk) / (its%it_0%xk - its%itdata(it0-1)%xk)
            dfdx_brack = (its%it_1%rk - its%it_0%rk)          / (its%it_1%xk - its%it_0%xk)
            dfdx_right = (its%itdata(it1+1)%rk - its%it_1%rk) / (its%itdata(it1+1)%xk - its%it_1%xk)

            brent_its_has_chng_stiff = (dfdx_right .gt. stiff_thrs*dfdx_left)

            ! in case stiffness changes: use extrapolation from both sides

            if (brent_its_has_chng_stiff) then
               x0_left  = its%it_0%xk - its%it_0%rk / dfdx_left
               x0_right = its%it_1%xk - its%it_1%rk / dfdx_right

               x0_new   = min(x0_left, x0_right)
               x0_new   = max(x0_new, its%it_0%xk+0.0001d0*(its%it_1%xk-its%it_0%xk))
            else
               x0_new   = 0.5d0 * (its%it_0%xk + its%it_1%xk)
            endif

         endif

         if (brent_its_has_chng_stiff .and. idebug_br.ge.1) then
            write(bufout,'(2(a,f12.6),3(a,g14.6))') '  ...stiffness changing in x=[', its%it_0%xk, ',', &
                its%it_1%xk,'], df/dx=', dfdx_left,',', dfdx_brack,',', dfdx_right
            call write_log(1, bufout)
         endif
      endif

   end function brent_its_has_chng_stiff

!------------------------------------------------------------------------------------------------------------

   subroutine brent_set_xnew(k, its, dfx_dx, used_bisec, tol_xk, x_new, idebug_br)
!--purpose: determine the next trial point using Brent's algorithm
      implicit none
!--subroutine arguments:
      integer,           intent(in)    :: k, idebug_br
      logical,           intent(inout) :: used_bisec
      real(kind=8),      intent(in)    :: dfx_dx, tol_xk
      type(t_brent_its)                :: its
      real(kind=8),      intent(out)   :: x_new
!--local variables:
      logical         :: use_inv_interp = .true.
      logical         :: ztest(5)
      real(kind=8)    :: dx, dr, dx_est, dx_max, dx_new

      associate(x_a => its%it_a%xk, x_b   => its%it_b%xk,   x_c   => its%it_c%xk,                       &
                x_k => its%it_k%xk, x_km1 => its%it_km1%xk, x_km2 => its%it_km2%xk,                     &
                r_a => its%it_a%rk, r_b   => its%it_b%rk,   r_c   => its%it_c%rk)

      !------------------------------------------------------------------------------------------------
      ! case 1: searching a bracket for solving Fz(z_ws) or Ftot(z_shift)
      !------------------------------------------------------------------------------------------------

      if (.not.brent_its_has_bracket(its) .and. abs(its%ikarg).eq.ikZDIR) then

         ! search bracket -- solving F_z or F_{tot,z}

         dr = its%itdata(its%numit)%rk - its%itdata(1)%rk

         if (dr.le.0d0 .and. k.le.5) then

            ! case 1.a: Fz decreasing across brent-table: 
            !                               double the step to reach positive dFz/dz

            x_new = its%itdata(1)%xk + 2d0 * (its%itdata(its%numit)%xk - its%itdata(1)%xk)
            if (idebug_br.ge.2) call write_log(' negative slope Fz, double step x(0)->x(k)')

         elseif (k.le.3) then

            ! case 1.b: start of iteration: extrapolate using current estimate dFz/dz (at end of table)
            !                               we could add a factor 1.05 to favour overshoot, to get a bracket
            !                               maximum step: dpen <= 0.5 * pen, new pen <= 1.5 * pen

            x_new  = x_b - r_b / dfx_dx
            dx_max = 0.5d0 * (x_b - its%itdata(1)%xk)

            if (idebug_br.ge.3 .and. x_new.gt.x_b+dx_max) then
               write(bufout,'(4(a,f9.4))') ' brent_set_xnew: x_0=', its%itdata(1)%xk,', x_k=', x_b,     &
                   ', x_new=',x_new,' >', x_b+dx_max
               call write_log(1, bufout)
            endif
            x_new  = min(x_new, x_b+dx_max)

         else

            ! after some iterations: extrapolate at end of table, using full table [x_0, x_n]
            !                        maximum step: new x in range [1.05, 1.2] * current interval

            dx     = its%itdata(its%numit)%xk - its%itdata(1)%xk
            dr     = its%itdata(its%numit)%rk - its%itdata(1)%rk

            dx_new = - its%itdata(its%numit)%rk * dx / dr
            dx_new = max(0.05d0*dx, min(0.2d0*dx, dx_new))
            x_new  = its%itdata(its%numit)%xk + dx_new

         endif

      !------------------------------------------------------------------------------------------------
      ! case 2: searching a bracket for solving Ftoty(y_shift)
      !------------------------------------------------------------------------------------------------

      elseif (.not.brent_its_has_bracket(its) .and. abs(its%ikarg).eq.ikYDIR) then

         ! y_shift <--> Fy: extrapolate at start or end of table, with point B == either x_0 or x_n
         !                  max. step == fac_dx * | x_n - x_0 |

         dx = its%itdata(its%numit)%xk - its%itdata(1)%xk
         dr = its%itdata(its%numit)%rk - its%itdata(1)%rk

         ! left side: expecting flange contact at y_shift < 0, Fy on rail < 0
         ! right side: expecting flange contact at y_shift > 0, Fy on rail > 0
         ! to increase Fy, use y_shift > 0

         if (its%itdata(1)%rk.gt.0d0) then

            ! expecting Fy(-50) < 0, extrapolate at start of table + 5% to shoot past zero
            !                        keep moving in direction dx even if dr suggests otherwise
            !                        require step dx_new \in [ -1.5, -0.2 ] * (x_n - x_0)

            dx_est = -r_b * dx / dr * 1.05d0
            dx_new = max(-1.5d0*dx, min(-0.2d0*dx, dx_est))
            x_new  =  its%itdata(1)%xk + dx_new

         else

            ! expecting Fy(50) > 0, extrapolate at end of table + 5% to shoot past zero
            !                       keep moving in direction dx even if dr suggests otherwise
            !                       require step dx_new \in [ 0.2, 1.5 ] * (x_n - x_0)

            dx = its%itdata(its%numit)%xk - its%itdata(1)%xk
            dr = its%itdata(its%numit)%rk - its%itdata(1)%rk

            dx_est = -r_b * dx / dr * 1.05d0
            dx_new = max(0.2d0*dx, min(1.5d0*dx, dx_est))
            x_new  =  its%itdata(its%numit)%xk + dx_new

         endif

         if (idebug_br.ge.3) then
            write(bufout,'(3(a,g12.4))') ' extrapolated dx_new=',dx_est,', clipped=',dx_new,            &
                     ', x_new=',x_new
            call write_log(1, bufout)
         endif

      !------------------------------------------------------------------------------------------------
      ! case 3: a bracket exists in the table
      !------------------------------------------------------------------------------------------------

      else

         if (use_inv_interp .and. abs(r_a-r_c).gt.tiny .and. abs(r_b-r_c).gt.tiny) then

            ! attempt inverse quadratic interpolation

            x_new = x_a * r_b * r_c / ((r_a-r_b) * (r_a-r_c)) +                                         &
                    x_b * r_a * r_c / ((r_b-r_a) * (r_b-r_c)) +                                         &
                    x_c * r_a * r_b / ((r_c-r_a) * (r_c-r_b))
            if (idebug_br.ge.3) then
               write(bufout,'(a,2g12.4,a,f12.4)') ' dfk=', abs(r_a-r_c), abs(r_b-r_c),                  &
                                      ', inverse interpolation: x_new =', x_new
               call write_log(1, bufout)
            endif

         elseif (k.le.2) then

            ! attempt the secant method using dfx_dx obtained from calling subroutine

            x_new = x_b - r_b / dfx_dx
            if (idebug_br.ge.3) then
               write(bufout,'(a,g12.4,10x,a,f12.6)') ' dfx_dx=', dfx_dx,                                &
                                      ', secant method: x_new =', x_new
               call write_log(1, bufout)
            endif

         elseif (abs(r_b-r_c).gt.tiny) then

            ! attempt the secant method using current and previous iterates

            x_new = x_b - r_b * (x_b-x_c) / (r_b-r_c)
            if (idebug_br.ge.3) then
               write(bufout,'(a,2g12.4,a,f12.6)') ' dfk=', abs(r_a-r_c), abs(r_b-r_c),                  &
                                      ', secant method: x_new =', x_new
               call write_log(1, bufout)
            endif
   
         else
   
            ! attempt linear interpolation (equals secant when c==a)
   
            x_new = x_b - r_b * (x_b-x_a) / (r_b-r_a)
            if (idebug_br.ge.3) then
               write(bufout,'(a,2g12.4,a,f12.4)') ' dfk=', abs(r_a-r_c), abs(r_b-r_c),  &
                                      ', linear interpolation: x_new =', x_new
               call write_log(1, bufout)
            endif
         endif

         ztest(1) = (x_new-(3d0*x_a+x_b)/4d0) * (x_new-x_b) .gt. 0d0
         ztest(2) =      used_bisec .and. abs(x_new-x_b) .ge. 0.5d0*abs(x_k -x_km1)
         ztest(3) = .not.used_bisec .and. abs(x_new-x_b) .ge. 0.5d0*abs(x_km1-x_km2) .and. k.gt.2
         ztest(4) =      used_bisec .and. abs(x_k  -x_km1) .lt. tol_xk
         ztest(5) = .not.used_bisec .and. abs(x_km1-x_km2) .lt. tol_xk .and. k.gt.2

         if (idebug_br.ge.2) then
            if (ztest(1)) then
               write(bufout,*) ' 1: x_new=',x_new,' not in [', (3d0*x_a+x_b)/4d0, ',', x_b, ']'
               call write_log(1, bufout)
            elseif (ztest(2)) then
               ! call brent_its_print(k, its, ic, idebug_br+2)
               write(bufout,*) ' 2: |x_new-its%x_b|=',abs(x_new-x_b),' >= ',0.5d0*abs(x_k-x_km1)
               call write_log(1, bufout)
            elseif (ztest(3)) then
               write(bufout,*) ' 3: |x_new-its%x_b|=',abs(x_new-x_b),' >= ',0.5d0*abs(x_km1-x_km2)
               call write_log(1, bufout)
            elseif (ztest(4)) then
               write(bufout,*) ' 4: |x_k  -its%x_km1|=',abs(x_k  -x_km1),' <= tol=', tol_xk
               call write_log(1, bufout)
            elseif (ztest(5)) then
               write(bufout,*) ' 5: |x_km1-its%x_km2|=',abs(x_km1-x_km2),' <= tol=', tol_xk
               call write_log(1, bufout)
            endif
         endif

         ! reject interpolation when any condition fails, use bisection instead

         if ( ztest(1) .or. ztest(2) .or. ztest(3) .or. ztest(4) .or. ztest(5) ) then
            if (idebug_br.ge.3) call write_log('  ...reject estimate, using bisection')
            x_new = 0.5d0 * (x_a + x_b)
            used_bisec = .true.
         else
            used_bisec = .false.
         endif

      endif ! not has_bracket
      end associate

   end subroutine brent_set_xnew

!------------------------------------------------------------------------------------------------------------

end module m_wr_brentmeth

