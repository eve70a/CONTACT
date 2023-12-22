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
!       we additionally use [x0, x1] for the bracket, and we store all iterates in xk_all (sorted)

   type :: t_brent_its
      integer                                 :: ikarg, maxit, numit
      real(kind=8)                            :: ftarg
      integer,      dimension(:), allocatable :: k_all
      real(kind=8), dimension(:), allocatable :: xk_all, rk_all
      integer                                 :: k_0, k_1, k_a, k_b, k_c, k_d
      real(kind=8), pointer     :: x_a => NULL(), x_b => NULL(), x_c => NULL(), x_d => NULL(),          &
                                   r_a => NULL(), r_b => NULL(), r_c => NULL(), r_d => NULL()
      real(kind=8), pointer     :: x_0 => NULL(), r_0 => NULL(), x_1 => NULL(), r_1 => NULL()
      real(kind=8), pointer     :: x_k => NULL(), x_km1 => NULL(), x_km2 => NULL()

      ! ikarg       meta-information for printing: ikarg==ikYDIR for solving Fy, ikZDIR when solving Fz
      ! maxit       size of arrays used, max. #iterates
      ! numit       actual #iterates stored
      ! ftarg       target force value
      ! xk_all      x-values used sorted in ascending order
      ! rk_all      residual values corresponding to x-values used
      ! k_all       iteration numbers corresponding to x-values used (starting at k=0?)
      ! k_0, k_1    iteration numbers of current bracket points
      ! x_0, r_0    lower side of current bracket
      ! k_a--k_d    iteration numbers of current/previous best guesses
      ! x_1, r_1    upper side of current bracket
      ! x_a, r_a    contra-point to current best guess (other end of bracket)
      ! x_b, r_b    current best guess
      ! x_c, r_c    previous best guess
      ! x_d, r_d    previous previous best guess
      ! x_k, x_km1, x_km2   most recent iterates x^k, x^{k-1}, x^{k-2}
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

      ! call write_log(' brent_its_init...')
      its%ikarg = ikarg
      its%maxit = max(1, min(10000, maxit))
      its%numit = 0
      its%ftarg = ftarg
      allocate(its%k_all(its%maxit))
      allocate(its%xk_all(its%maxit))
      allocate(its%rk_all(its%maxit))

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
      type(t_brent_its)         :: its
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

   subroutine brent_its_add_iterate(its, k, xk, rk, my_ierror)
!--purpose: insert an iterate in the sorted structure for Brent's algorithm
      implicit none
!--subroutine arguments:
      type(t_brent_its)      :: its
      integer                :: k, my_ierror
      real(kind=8)           :: xk, rk
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

      i     = 0
      ldone = (i.ge.its%numit)
      do while (.not.ldone)
         i     = i + 1
         ldone = .true.
         if (i.le.its%numit) ldone = (xk.lt.its%xk_all(i))
      enddo

      ! shift iterates [i--end] one position

      do j = its%numit, max(i,1), -1
         its%k_all(j+1)  = its%k_all(j)
         its%xk_all(j+1) = its%xk_all(j)
         its%rk_all(j+1) = its%rk_all(j)
      enddo

      ! insert iterate at position i

      if (i.le.0) i = 1

      its%numit     = its%numit + 1
      its%k_all(i)  = k
      its%xk_all(i) = xk
      its%rk_all(i) = rk

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
         if (its%k_all(j).eq.k) it = j
      enddo

      if (it.eq.-1) then
         dfdx = -99d9
      elseif (its%numit.le.1) then
         dfdx = 0d0
      elseif (it.le.1) then
         ! forward difference
         dfdx = (its%rk_all(it+1) - its%rk_all(it)) / max(tiny, its%xk_all(it+1) - its%xk_all(it))
      elseif (it.ge.its%numit) then
         ! backward difference
         dfdx = (its%rk_all(it) - its%rk_all(it-1)) / max(tiny, its%xk_all(it) - its%xk_all(it-1))
         ! initial estimate for dFz/dpen: Fz ~ pen^1.5
         if (its%ikarg.eq.ikZDIR .and. it.eq.2) dfdx = 1.5d0 * dfdx
      elseif (.false.) then
         ! central difference
         dfdx = (its%rk_all(it+1) - its%rk_all(it-1)) / max(tiny, its%xk_all(it+1) - its%xk_all(it-1))
      elseif (abs(its%rk_all(it+1)).lt.abs(its%rk_all(it-1))) then
         ! forward difference
         dfdx = (its%rk_all(it+1) - its%rk_all(it)) / max(tiny, its%xk_all(it+1) - its%xk_all(it))
      else
         ! backward difference
         dfdx = (its%rk_all(it) - its%rk_all(it-1)) / max(tiny, its%xk_all(it) - its%xk_all(it-1))
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

         write(str12(1), '(f12.4)') its%xk_all(it)
         str12(2) = fmt_gs(12,4, its%rk_all(it)+its%ftarg)
         str12(3) = fmt_gs(12,4, dfdx)
         if (its%ikarg.eq.ikYDIR) then
            write(bufout,1020) its%k_all(it), (str12(j),j=1,3)
         else
            write(bufout,1030) its%k_all(it), (str12(j),j=1,3)
         endif
         call write_log(1, bufout)

 1020    format(2x, i4,', BR,  dyrail, Fy: ',2a,', dFy/dy:',a)
 1030    format(4x, i6,', NR,  z_ws, Fz: ',2a,', dFz/dz:',a)

      endif

      if (ic%flow.ge.1 .and. idebug_br.eq.2) then

         ! print one line with current iterate / residual

         it   = brent_its_find_k(its, k)
         dfdx = brent_sensitivity_k(k, its, idebug_br)

         write(str12(1), '(f12.4)') its%xk_all(it)
         str12(2) = fmt_gs(12,4, its%rk_all(it)+its%ftarg)
         str12(3) = fmt_gs(12,4, its%rk_all(it))
         str12(4) = fmt_gs(12,4, dfdx)

         if (its%ikarg.eq.ikYDIR) then
            write(bufout,2020) its%k_all(it), (str12(j),j=1,4)
         else
            write(bufout,2030) its%k_all(it), (str12(j),j=1,4)
         endif
         call write_log(1, bufout)
      endif

      if (idebug_br.ge.3) then

         ! print table with all iterates and residuals

         write(bufout,'(a,i3,2a)') ' --- brent table for',its%numit,' iterations. B=best, A=contra, ',  &
                'C,D=prev best, bracket [x0,x1] --------------------'
         call write_log(1, bufout)

         do it = 1, its%numit
            dfdx = brent_sensitivity_k(its%k_all(it), its, idebug_br)

            strptr = ' <--'
            if (its%k_all(it).eq.its%k_0) strptr = trim(strptr) // ' x0,'
            if (its%k_all(it).eq.its%k_1) strptr = trim(strptr) // ' x1,'
            if (its%k_all(it).eq.its%k_a) strptr = trim(strptr) // ' A,'
            if (its%k_all(it).eq.its%k_b) strptr = trim(strptr) // ' B,'
            if (its%k_all(it).eq.its%k_c) strptr = trim(strptr) // ' C,'
            if (its%k_all(it).eq.its%k_d) strptr = trim(strptr) // ' D,'
            if (idebug_br.le.1 .or. len(trim(strptr)).le.4) strptr = ' '
            ilen = max(1, len(trim(strptr)))

            if (idebug_br.ge.6 .or. .true.) then
               str16(1) = fmt_gs(16,8, its%xk_all(it))
               str16(2) = fmt_gs(16,8, its%rk_all(it)+its%ftarg)
               str16(3) = fmt_gs(16,8, its%rk_all(it))
               str16(4) = fmt_gs(16,8, dfdx)
               if (its%ikarg.eq.ikYDIR) then
                  write(bufout,2020) its%k_all(it), (str16(j),j=1,4), strptr(1:ilen-1)
               else
                  write(bufout,2030) its%k_all(it), (str16(j),j=1,4), strptr(1:ilen-1)
               endif
            else
               write(str12(1), '(f12.4)') its%xk_all(it)
               str12(2) = fmt_gs(12,4, its%rk_all(it)+its%ftarg)
               str12(3) = fmt_gs(12,4, its%rk_all(it))
               str12(4) = fmt_gs(12,4, dfdx)
               if (its%ikarg.eq.ikYDIR) then
                  write(bufout,2020) its%k_all(it), (str12(j),j=1,4), strptr(1:ilen-1)
               else
                  write(bufout,2030) its%k_all(it), (str12(j),j=1,4), strptr(1:ilen-1)
               endif
            endif
            call write_log(1, bufout)
         enddo

 2020    format(2x, i4,', BR,  dy, Fy, res:',3a,', dFy/dy:',2a)
 2030    format(4x, i6,', NR,  z_ws, Fz:',2a,', res: ',a,', dFz/dz:',2a)

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
      logical         :: use_inv_interp = .false.
      logical         :: ztest(5)
      real(kind=8)    :: dx, dr, dx_est, dx_max, dx_new

      !------------------------------------------------------------------------------------------------
      ! case 1: searching a bracket for solving z_ws - Fz
      !------------------------------------------------------------------------------------------------

      if (.not.brent_its_has_bracket(its) .and. its%ikarg.eq.ikZDIR) then

         ! search bracket -- solving Fz

         dr = its%rk_all(its%numit) - its%rk_all(1)

         if (dr.lt.0d0 .and. k.le.5) then

            ! case 1.a: Fz decreasing across brent-table: double the step to reach positive dFz/dz

            x_new = its%xk_all(1) + 2d0 * (its%xk_all(its%numit) - its%xk_all(1))
            if (idebug_br.ge.2) call write_log(' negative force Fz, double step x(0)->x(k)')

         elseif (k.le.3) then

            ! case 1.b: start of iteration: extrapolate using current estimate dfx_dx (at end of table)
            !                               maximum step: dpen <= 0.5 * pen, new pen <= 1.5 * pen

            x_new  = its%x_b - its%r_b / dfx_dx
            dx_max = 0.5d0 * (its%x_b - its%xk_all(1))

            if (idebug_br.ge.3 .and. x_new.gt.its%x_b+dx_max) then
               write(bufout,'(4(a,f9.4))') ' brent_set_xnew: x_0=', its%xk_all(1),', x_k=', its%x_b,       &
                   ', x_new=',x_new,' >', its%x_b+dx_max
               call write_log(1, bufout)
            endif
            x_new  = min(x_new, its%x_b+dx_max)

         else

            ! after some iterations: extrapolate at end of table, using full table [x_0, x_n]
            !                        maximum step: new x in range [1.05, 1.2] * current interval

            dx     = its%xk_all(its%numit) - its%xk_all(1)
            dr     = its%rk_all(its%numit) - its%rk_all(1)

            dx_new = - its%rk_all(its%numit) * dx / dr
            dx_new = max(0.05d0*dx, min(0.2d0*dx, dx_new))
            x_new  = its%xk_all(its%numit) + dx_new

            ! write(bufout,*) 'dx_new=',dx_new,', x_new=',x_new
            ! call write_log(1, bufout)

         endif

      !------------------------------------------------------------------------------------------------
      ! case 2: searching a bracket for solving dy_rail - Fy
      !------------------------------------------------------------------------------------------------

      elseif (.not.brent_its_has_bracket(its) .and. its%ikarg.eq.ikYDIR) then

         ! dy_rail <--> Fy: extrapolate at start or end of table, with point B == either x_0 or x_n
         !                  max. step == fac_dx * | x_n - x_0 |

         dx = its%xk_all(its%numit) - its%xk_all(1)
         dr = its%rk_all(its%numit) - its%rk_all(1)

         if (its%rk_all(1).lt.0d0) then

            ! expecting Fy(-50) > 0, extrapolate at start of table + 5% to shoot past zero
            !                        require step dx_new \in [ -1.5, -0.1 ] * (x_n - x_0)

            dx_est = -its%r_b * dx / dr * 1.05d0
            dx_new = max(-1.5d0*dx, min(-0.1d0*dx, dx_est))
            x_new  =  its%xk_all(1) + dx_new

         else

            ! expecting Fy(50) < 0, extrapolate at end of table + 5% to shoot past zero
            !                       require step dx_new \in [ 0.1, 1.5 ] * (x_n - x_0)

            dx = its%xk_all(its%numit) - its%xk_all(1)
            dr = its%rk_all(its%numit) - its%rk_all(1)

            dx_est = -its%r_b * dx / dr * 1.05d0
            dx_new = max(0.1d0*dx, min(1.5d0*dx, dx_est))
            x_new  =  its%xk_all(its%numit) + dx_new

         endif

         if (idebug_br.ge.3) then
            write(bufout,'(3(a,g12.4))') ' extrapolated dx_new=',dx_est,', clipped=',dx_new,         &
                     ', x_new=',x_new
            call write_log(1, bufout)
         endif

      !------------------------------------------------------------------------------------------------
      ! case 3: a bracket exists in the table
      !------------------------------------------------------------------------------------------------

      else

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

            ! attempt the secant method using dfx_dx obtained from calling subroutine

            x_new = its%x_b - its%r_b / dfx_dx
            if (idebug_br.ge.3) then
               write(bufout,'(a,g12.4,10x,a,f12.6)') ' dfx_dx=', dfx_dx,                                &
                                      ', secant method: x_new =', x_new
               call write_log(1, bufout)
            endif

         elseif (abs(its%r_b-its%r_c).gt.tiny) then

            ! attempt the secant method using current and previous iterates

            x_new = its%x_b - its%r_b * (its%x_b-its%x_c) / (its%r_b-its%r_c)
            if (idebug_br.ge.3) then
               write(bufout,'(a,2g12.4,a,f12.6)') ' dfk=', abs(its%r_a-its%r_c), abs(its%r_b-its%r_c),  &
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
               write(bufout,*) ' 1: |x_new=',x_new,' not in [', (3d0*its%x_a+its%x_b)/4d0, ',', its%x_b, ']'
               call write_log(1, bufout)
            elseif (ztest(2)) then
               ! call brent_its_print(k, its, ic, idebug_br+2)
               write(bufout,*) ' 2: |x_new-its%x_b|=',abs(x_new-its%x_b),' >= ',0.5d0*abs(its%x_k-its%x_km1)
               call write_log(1, bufout)
            elseif (ztest(3)) then
               write(bufout,*) ' 3: |x_new-its%x_b|=',abs(x_new-its%x_b),' >= ',0.5d0*abs(its%x_km1-its%x_km2)
               call write_log(1, bufout)
            elseif (ztest(4)) then
               write(bufout,*) ' 4: |x_k  -its%x_km1|=',abs(its%x_k  -its%x_km1),' <= tol=', tol_xk
               call write_log(1, bufout)
            elseif (ztest(5)) then
               write(bufout,*) ' 5: |x_km1-its%x_km2|=',abs(its%x_km1-its%x_km2),' <= tol=', tol_xk
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

end module m_wr_brentmeth

