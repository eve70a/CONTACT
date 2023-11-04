!------------------------------------------------------------------------------------------------------------
! m_bspline_def - definitions regarding 1D B-splines
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_bspline_def
   use m_globals
   use m_markers
   use m_ptrarray
   use m_interp
   use m_spline_def
   implicit none
   private

   ! Debugging for module m_bspline_def

   integer  :: ldebug    =  0    ! local level of debugging
   integer  :: ii_debug  = -1    ! output point for which detailed info is requested (-1 = none)
   integer  :: iel_debug = -1    ! input element for which detailed info is requested (-1 = none)
   public  bsplinedef_set_debug

   ! global parameter for B-spline knots

   real(kind=8), parameter    :: tiny_dt = 1d-10
   public tiny_dt

   ! Data type for 2D B-splines:

   public t_bspline2d

   ! Functions for the creation & evaluation of 2D B-splines:

   public  bspline2d_nullify
   public  bspline2d_destroy
   public  bspline2d_print
   public  bspline_make_breakpoints
   public  bspline_insert_knot
   public  bspline_eval1d
   public  bspline_make_ppform
   public  bspline_check_kink_accel
   public  bspline_check_mask

!------------------------------------------------------------------------------------------------------------
!  data for a parametric tensor B-spline (u, v, [x(u,v)], y(u,v), z(u,v)):

   type :: t_bspline2d
      integer      :: nknotu, nknotv, nreducu, nreducv, nsplu, nsplv, nbrku, nbrkv
      logical      :: has_xdata

      real(kind=8), dimension(:),   pointer  :: tui    => NULL() ! (nknotu)
      real(kind=8), dimension(:),   pointer  :: tvj    => NULL() ! (nknotv)
      logical,      dimension(:),   pointer  :: keeptu => NULL() ! (nknotu)
      logical,      dimension(:),   pointer  :: keeptv => NULL() ! (nknotv)

      real(kind=8), dimension(:),   pointer  :: ubrk   => NULL() ! (nbrku)
      real(kind=8), dimension(:),   pointer  :: vbrk   => NULL() ! (nbrkv)

      real(kind=8), dimension(:,:), pointer  :: cij_x  => NULL() ! (nsplu,nsplv)
      real(kind=8), dimension(:,:), pointer  :: cij_y  => NULL() ! (nsplu,nsplv)
      real(kind=8), dimension(:,:), pointer  :: cij_z  => NULL() ! (nsplu,nsplv)
      integer,      dimension(:,:), pointer  :: mask   => NULL() ! (nsplu,nsplv)

      ! nknotu   number of knots used in u-direction
      ! nknotv   number of knots used in v-direction
      ! tui      knots in u-direction
      ! tvj      knots in v-direction
      ! keeptu   flags indicating kept knots in u-direction
      ! keeptv   flags indicating kept knots in v-direction
      ! nreducu  number of knots kept in u-direction after reduction
      ! nreducv  number of knots kept in v-direction after reduction

      ! has_xdata  flag indicating that cij_x is present (3d space curve or 3d surface)
      ! nsplu    number of splines used in u-direction
      ! nsplv    number of splines used in v-direction
      ! cij_x    tensor B-spline coefficients for x(u,v)
      ! cij_y    tensor B-spline coefficients for y(u,v)
      ! cij_z    tensor B-spline coefficients for z(u,v)
      ! mask     flags indicating which cij are available (1) or missing (0)

      ! nbrku    number of break points in u-direction
      ! nbrkv    number of break points in v-direction
      ! ubrk     break points in u-direction
      ! vbrk     break points in v-direction

   end type t_bspline2d

contains

!------------------------------------------------------------------------------------------------------------

subroutine bsplinedef_set_debug(new_ldebug, new_ii_debug, new_iel_debug)
!--function: enable/disable debug output of spline routines
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
      write(bufout,'(a,i3,2(a,i7))') ' bspline_def:  debugging level =',ldebug,', ii_debug =', ii_debug, &
                ', iel_debug =', iel_debug
      call write_log(1, bufout)
   endif

end subroutine bsplinedef_set_debug

!------------------------------------------------------------------------------------------------------------

subroutine bspline2d_nullify(bspl)
!--purpose: initialize B-spline structure, nullify pointers
   implicit none
!--subroutine parameters:
   type(t_bspline2d) :: bspl

   bspl%nknotu  =  0
   bspl%nknotv  =  0
   bspl%nreducu =  0
   bspl%nreducv =  0
   bspl%nsplu   =  0
   bspl%nsplv   =  0
   bspl%nbrku   =  0
   bspl%nbrkv   =  0
   bspl%has_xdata = .false.

   bspl%tui    => NULL()
   bspl%tvj    => NULL()
   bspl%keeptu => NULL()
   bspl%keeptv => NULL()

   bspl%nbrku   =  0
   bspl%nbrkv   =  0
   bspl%ubrk   => NULL()
   bspl%vbrk   => NULL()

   bspl%nsplu   =  0
   bspl%nsplv   =  0
   bspl%cij_x  => NULL()
   bspl%cij_y  => NULL()
   bspl%cij_z  => NULL()
   bspl%mask   => NULL()

end subroutine bspline2d_nullify

!------------------------------------------------------------------------------------------------------------

subroutine bspline2d_destroy(bspl)
!--purpose: clean-up allocated arrays for B-spline, nullify alias-pointers
   implicit none
!--subroutine parameters:
   type(t_bspline2d)  :: bspl

   if (associated(bspl%tui))    deallocate(bspl%tui)
   if (associated(bspl%tvj))    deallocate(bspl%tvj)
   if (associated(bspl%keeptu)) deallocate(bspl%keeptu)
   if (associated(bspl%keeptv)) deallocate(bspl%keeptv)
   if (associated(bspl%ubrk))   deallocate(bspl%ubrk)
   if (associated(bspl%vbrk))   deallocate(bspl%vbrk)
   if (associated(bspl%cij_x))  deallocate(bspl%cij_x)
   if (associated(bspl%cij_y))  deallocate(bspl%cij_y)
   if (associated(bspl%cij_z))  deallocate(bspl%cij_z)
   if (associated(bspl%mask))   deallocate(bspl%mask)

   call bspline2d_nullify(bspl)

end subroutine bspline2d_destroy

!------------------------------------------------------------------------------------------------------------

subroutine bspline2d_print(bspl, nam, idebug, ndigit)
!--function: print information on spline bspl
   implicit none
!--subroutine arguments
   type(t_bspline2d)    :: bspl
   character(len=*)     :: nam
   integer              :: idebug
   integer, optional    :: ndigit       ! number of significant digits
!--local variables
   ! character(len= 1), parameter  :: cnams(1:3) = (/ 'x', 'y', 'z' /)
   ! character(len=30), parameter  :: spaces = '                              '
   integer              :: my_ndigit, my_len
   ! character(len=18)    :: strng(5)

   if (present(ndigit)) then
      my_ndigit = ndigit
   else
      my_ndigit = 6
   endif
   my_ndigit = max(2, min(10, my_ndigit))
   my_len    = 8 + my_ndigit

   if (.not.associated(bspl%tui) .or. .not.associated(bspl%cij_y)) then
      write(bufout,'(3a)') ' spline ',trim(nam),' has not yet been defined'
      call write_log(1, bufout)
      return
   endif

   ! idebug>=6: array sizes

   if (idebug.ge.6) then
      call print_array_size_1d(bspl%tui,   'tui')
      call print_array_size_1d(bspl%tvj,   'tvj')
      ! call print_array_size_1d(bspl%keeptu,'keeptu')
      ! call print_array_size_1d(bspl%keeptv,'keeptv')
      call print_array_size_1d(bspl%ubrk,  'ubrk')
      call print_array_size_1d(bspl%vbrk,  'vbrk')
      call print_array_size_2d(bspl%cij_x, 'cij_x')
      call print_array_size_2d(bspl%cij_y, 'cij_y')
      call print_array_size_2d(bspl%cij_z, 'cij_z')
   endif

   ! idebug>=5: full details w.r.t. spline data

   if (idebug.ge.5) then

      call write_log(trim(nam) // ' = struct;')

      ! knots ti, breaks in x-direction

      write(bufout,'(2a,i5,a)') trim(nam), '.nknotu =', bspl%nknotu,';'
      call write_log(1, bufout);

      call print_1d_real(bspl%nknotu, bspl%tui, trim(nam)//'.tui', 10, 'f13.4')

      write(bufout,'(2a,i5,a)') trim(nam), '.nbrku = ', bspl%nbrku,';'
      call write_log(1, bufout);

      call print_1d_real(bspl%nbrku, bspl%ubrk, trim(nam)//'.ubrk', 10, 'f13.4')

      ! knots tj, breaks in u-direction

      write(bufout,'(2a,i5,a)') trim(nam), '.nknotv =', bspl%nknotv,';'
      call write_log(1, bufout);

      call print_1d_real(bspl%nknotv, bspl%tvj, trim(nam)//'.tvj', 10, 'f12.4')

      write(bufout,'(2a,i5,a)') trim(nam), '.nbrkv = ', bspl%nbrkv,';'
      call write_log(1, bufout);

      call print_1d_real(bspl%nbrkv, bspl%vbrk, trim(nam)//'.ubrk', 10, 'f12.4')

      ! coefficients cij_[xyz]      ! TODO: use fmt_gs for requested #digits

      if (bspl%has_xdata) then
         call print_2d_real(bspl%nsplu, bspl%nsplv, bspl%cij_x, trim(nam)//'.cij_x', 10, 'g12.4')
      endif
      call print_2d_real(bspl%nsplu, bspl%nsplv, bspl%cij_y, trim(nam)//'.cij_y', 10, 'g12.4')
      call print_2d_real(bspl%nsplu, bspl%nsplv, bspl%cij_z, trim(nam)//'.cij_z', 10, 'g12.4')
      call print_2d_int(bspl%nsplu, bspl%nsplv, bspl%mask, trim(nam)//'.mask_j', 40, 'i3')

   endif ! idebug >= 5

end subroutine bspline2d_print

!------------------------------------------------------------------------------------------------------------

subroutine bspline_make_breakpoints(nknot, tj, tiny_dt, nbrk_arg, spnt, idebug)
!--function: determine breakpoints in knot vector tj: unique knots, excluding 1:k-1 and end-[0:k-2]
   implicit none
!--subroutine arguments:
   integer,      intent(in)     :: nknot, nbrk_arg, idebug
   real(kind=8), intent(in)     :: tj(nknot), tiny_dt
   real(kind=8), intent(out)    :: spnt(nbrk_arg)
!--local variables:
   integer      :: nbrk, ipnt, j

   nbrk = 1
   spnt(nbrk) = tj(4)
   do j = 5, nknot-3
      if (tj(j)-tj(j-1).gt.tiny_dt) then
         nbrk = nbrk + 1
         if (nbrk.le.nbrk_arg) spnt(nbrk) = tj(j)
      endif
   enddo

   if (idebug.ge.1) then
      write(bufout,'(a,i4,a)') ' PP-spline has ',nbrk,' breakpoints ksi'
      call write_log(1, bufout)
   endif

   if (nbrk.ne.nbrk_arg) then
      call write_log(' Internal error(make_breakpoints): incorrect nbrk.')
      write(bufout,'(2(a,i4))') ' expected ',nbrk_arg,' breakpoints, found', nbrk
      call write_log(1, bufout)
      call abort_run()
   endif

   if (idebug.eq.2) then
      do ipnt = 1, nbrk
         write(bufout,'(a,i4,a,f12.6)') 'i =',ipnt,': ksi =',spnt(ipnt)
         call write_log(1, bufout)
      enddo
   endif

end subroutine bspline_make_breakpoints

!------------------------------------------------------------------------------------------------------------

subroutine bspline_insert_knot(nknot, tj, tnew, tjx, nrow, ncol, coef)
!--function: shift and update B-spline coefficients for new knot tnew
   implicit none
!--subroutine arguments:
   integer,      intent(in)           :: nknot, nrow, ncol
   real(kind=8), intent(in)           :: tj(nknot), tnew
   real(kind=8), intent(inout)        :: coef(nrow,ncol)
   real(kind=8), intent(out)          :: tjx(nknot+1)
!--local variables:
   integer, parameter :: k = 4  ! spline order
   integer            :: nspl, jnew, j
   real(kind=8)       :: fj

   ! on input, coef holds nspl rows of data

   nspl = nknot - k

   ! array coef must have space for an additional row, i.e. nrow >= nknot+1

   if (nrow.le.nspl) then
      call write_log(' Internal error(bspline_insert_knot): nrow <= nknot-k')
      write(bufout,'(4(a,i4))') ' nknot=',nknot,', k=',k,', nspl=',nspl,', nrow=',nrow
      call write_log(1, bufout)
      call abort_run()
   endif

   ! determine position jnew for tnew in the extended knot vector

   jnew = 1
   do while(jnew.lt.nknot .and. tnew.gt.tj(jnew))
      jnew = jnew + 1
   enddo

   if (jnew.le.1 .or. jnew.ge.nknot) then
      write(bufout,'(a,i4)') ' nknot=',nknot
      call write_log(1, bufout)
      call write_log(' Internal error(bspline_insert_knot): jnew<=1 or >=nknot.')
      if (jnew.le.1) then
         write(bufout,'(a,i4,a,f8.3,a,4f8.3)') ' jnew=',jnew,', tnew=',tnew,', tj=',(tj(j), j=1,4)
      else
         write(bufout,'(2(a,i4),a,f8.3,a,4f8.3)') ' jnew=',jnew,', nknot=',nknot,', tnew=',tnew,        &
                ', tj=',(tj(j), j=nknot-3,nknot)
      endif
      call write_log(1, bufout)
      call abort_run()
   endif

   ! determine extended knot vector

   tjx(     1:jnew-1)  = tj(1:jnew-1)
   tjx(jnew  )         = tnew
   tjx(jnew+1:nknot+1) = tj(jnew:nknot)
   
   ! compute new coefficients last to first, overwriting input values

   do j = nspl+1, 1, -1

      if (j.ge.jnew) then

         ! j = jnew : nspl_out

         coef(j,1:ncol) = coef(j-1,1:ncol)

      elseif (j.ge.jnew-k+1) then

         ! j = jnew-k+1 : jnew-1

         fj = (tnew - tj(j)) / (tj(j+3) - tj(j))
         coef(j,1:ncol) = (1d0 - fj) * coef(j-1,1:ncol) + fj * coef(j,1:ncol)

      else

         ! j = 1 : jnew-k:  a_out(j) = a_in(j)

      endif
   enddo

end subroutine bspline_insert_knot

!------------------------------------------------------------------------------------------------------------

subroutine bspline_eval1d(nknot, tj, tiny_dt, npnt, si, jseg, b1, b2, b3, b4, idebug, my_ierror)
!--function: evaluate B-spline basis functions B_{j,k} for knot-vector tj at sample locations si
   implicit none
!--subroutine arguments:
   integer,      intent(in)    :: nknot, npnt, idebug
   real(kind=8), intent(in)    :: tj(nknot), tiny_dt, si(npnt)
   integer,      intent(out)   :: jseg(npnt)
   real(kind=8), intent(out)   :: b1(npnt,1), b2(npnt,2), b3(npnt,3), b4(npnt,4)
   integer,      intent(out)   :: my_ierror
!--local variables:
   integer              :: ip, j, ix, jj1, jj2, jj3, jj4
   real(kind=8)         :: dtj_k1l, dtj_k1r, dtj_k2l, dtj_k2r, dtj_k3l, dtj_k3r, s_ip,                  &
                           dtj_inv1l, dtj_inv1r, dtj_inv2l, dtj_inv2r, dtj_inv3l, dtj_inv3r

   my_ierror = 0

   ! check all s_i in range [t(4), t(n-3)] (cubic spline, k=4)

   if (ldebug.ge.1) then
      do ip = 1, npnt
         if (si(ip).lt.tj(4)-tiny_dt .or. si(ip).gt.tj(nknot-3)+tiny_dt) then
            write(bufout,'(a,i4,3(a,f8.3),a)') ' Warning(bspline_eval1d): point si(',ip,')=',si(ip),    &
                ' lies outside tj(4:n-3) = [', tj(4),',',tj(nknot-3),']'
            call write_log(1, bufout)
         endif
      enddo
   endif

   ! for each ip: determine first true interval where s_i \in [tj-1,tj], store j-1
   ! TODO: this subroutine should signal s_i outside the basic interval [tj(4),tj(n-k+1)],
   !       e.g. return jseg<4 for points that lie before tj(4), jseg>n-k for points after tj(n-k+1)

   do ip = 1, npnt
      j = 5
      do while (j.lt.nknot-3 .and. si(ip).ge.tj(j))
         j = j + 1
      enddo
      jseg(ip) = j - 1
      if (jseg(ip).gt.nknot-4) then
         call write_log('INTERNAL ERROR: jseg > nknot-4')
         call abort_run()
      endif
   enddo

   if (ldebug.ge.3 .and. ii_debug.ge.1 .and. ii_debug.le.npnt) then
      ip = ii_debug
      write(bufout,'(a,i4,a,f12.4,a,i4,2(a,f12.4),a)') ' eval1d: point iout=',ip,', si=',si(ip),        &
                ' lies in true segment', jseg(ip),', tj=[',tj(jseg(ip)),',', tj(jseg(ip)+1),')'
      call write_log(1, bufout)
   endif

   ! determine value of each basis-function at each measurement location si

   b1(1:npnt,1:1) = 0d0
   b2(1:npnt,1:2) = 0d0
   b3(1:npnt,1:3) = 0d0
   b4(1:npnt,1:4) = 0d0

   do ip = 1, npnt

      ! clip si(ip) \in [t(4), t(n-3)] -- constant extrapolation

      s_ip = min( max(si(ip), tj(4)), tj(nknot-3) )

      ! piecewise constant B_{j,1}

      do j = jseg(ip), jseg(ip)
         if (j.eq.nknot-4) then
            if (s_ip.ge.tj(j) .and. s_ip.le.tj(j+1)) b1(ip,j-jseg(ip)+1) = 1d0
         else
            if (s_ip.ge.tj(j) .and. s_ip.lt.tj(j+1)) b1(ip,j-jseg(ip)+1) = 1d0
         endif
      enddo

      ! piecewise linear B_{j,2}

      do j = jseg(ip)-1, jseg(ip)

         ! step length over 1 interval, inverse 0 for length 0

         dtj_k1l   = tj(j+1) - tj(j)
         dtj_k1r   = tj(j+2) - tj(j+1)
         dtj_inv1l = 0d0
         dtj_inv1r = 0d0
         if (dtj_k1l.gt.tiny_dt) dtj_inv1l = 1d0 / dtj_k1l
         if (dtj_k1r.gt.tiny_dt) dtj_inv1r = 1d0 / dtj_k1r
         
         ! b1 is stored only for j = jseg(ip), b2 is stored for jseg(ip) + [-1,0]

         jj1 = j - jseg(ip) + 1
         jj2 = j - jseg(ip) + 2

         b2(ip,jj2) = 0d0
         if (jj1.eq.1)   b2(ip,jj2) = b2(ip,jj2) + (s_ip  - tj(j))  * b1(ip,jj1)   * dtj_inv1l 
         if (jj1+1.eq.1) b2(ip,jj2) = b2(ip,jj2) + (tj(j+2) - s_ip) * b1(ip,jj1+1) * dtj_inv1r
      enddo

      ! piecewise quadratic B_{j,3}

      do j = jseg(ip)-2, jseg(ip)

         ! average step length over 2 intervals, inverse 0 for length 0

         dtj_k2l   = (tj(j+2) - tj(j)  ) / 2d0
         dtj_k2r   = (tj(j+3) - tj(j+1)) / 2d0
         dtj_inv2l = 0d0
         dtj_inv2r = 0d0
         if (dtj_k2l.gt.tiny_dt) dtj_inv2l = 1d0 / dtj_k2l
         if (dtj_k2r.gt.tiny_dt) dtj_inv2r = 1d0 / dtj_k2r
         
         ! b2 is stored only for j = jseg(ip) + [-1,0], b3 is stored for jseg(ip) + [-2:0]

         jj2 = j - jseg(ip) + 2
         jj3 = j - jseg(ip) + 3

         b3(ip,jj3) = 0d0
         if (jj2.ge.1 .and. jj2.le.2)                                                                   &
            b3(ip,jj3) = b3(ip,jj3) + (s_ip  - tj(j))  * b2(ip,jj2)   * dtj_inv2l / 2d0
         if (jj2+1.ge.1 .and. jj2+1.le.2)                                                               &
            b3(ip,jj3) = b3(ip,jj3) + (tj(j+3) - s_ip) * b2(ip,jj2+1) * dtj_inv2r / 2d0
      enddo

      ! piecewise cubic B_{j,4}

      ! do j = 1, nknot-4
      do j = jseg(ip)-3, jseg(ip)

         ! average step length over 3 intervals, inverse 0 for length 0

         dtj_k3l   = (tj(j+3) - tj(j)  ) / 3d0
         dtj_k3r   = (tj(j+4) - tj(j+1)) / 3d0
         dtj_inv3l = 0d0
         dtj_inv3r = 0d0
         if (dtj_k3l.gt.tiny_dt) dtj_inv3l = 1d0 / dtj_k3l
         if (dtj_k3r.gt.tiny_dt) dtj_inv3r = 1d0 / dtj_k3r

         ! b3 is stored only for j = jseg(ip) + [-2:0], b4 is stored for jseg(ip) + [-3:0]

         jj3 = j - jseg(ip) + 3
         jj4 = j - jseg(ip) + 4

         b4(ip,jj4) = 0d0
         if (jj3.ge.1 .and. jj3.le.3)                                                                   &
            b4(ip,jj4) = b4(ip,jj4) + (s_ip  - tj(j))  * b3(ip,jj3)   * dtj_inv3l / 3d0
         if (jj3+1.ge.1 .and. jj3+1.le.3)                                                               &
            b4(ip,jj4) = b4(ip,jj4) + (tj(j+4) - s_ip) * b3(ip,jj3+1) * dtj_inv3r / 3d0
      enddo

   enddo  ! ip

   if (idebug.ge.5 .and. nknot.le.16) then
      call write_log('basisfunctions B_{j,1}:')
      write(bufout,'(a, 15i9)') '  s_i  jseg', (j, j=1, min(nknot-1,15) )
      call write_log(1, bufout)
      do ip = 1, npnt
                                                                ! TODO: j=1:15
         write(bufout,'(f5.1, i5,1x, 15f9.4)') si(ip), jseg(ip), (-1d0, j=1,jseg(ip)-1),                &
                                                (b1(ip,j), j=1,1), (-1d0, j=jseg(ip)+1,nknot-1)
         ix = index(bufout(1), '-1.0000')
         do while (ix.gt.0)
            bufout(1)(ix:ix+6) = '      x'
            ix = index(bufout(1), '-1.0000')
         enddo
         call write_log(1, bufout)
      enddo

      call write_log('basisfunctions B_{j,2}:')
      do ip = 1, npnt
                                                                ! TODO: j=1:15
         write(bufout,'(f5.1, i5,1x, 15f9.4)') si(ip), jseg(ip), (-1d0, j=1,jseg(ip)-2),                &
                                                (b2(ip,j), j=1,2), (-1d0, j=jseg(ip)+1,nknot-2)
         ix = index(bufout(1), '-1.0000')
         do while (ix.gt.0)
            bufout(1)(ix:ix+6) = '      x'
            ix = index(bufout(1), '-1.0000')
         enddo
         call write_log(1, bufout)
      enddo

      call write_log('basisfunctions B_{j,3}:')
      do ip = 1, npnt
                                                                ! TODO: j=1:15
         write(bufout,'(f5.1, i5,1x, 15f9.4)') si(ip), jseg(ip), (-1d0, j=1,jseg(ip)-3),                &
                                                (b3(ip,j), j=1,3), (-1d0, j=jseg(ip)+1,nknot-3)
         ix = index(bufout(1), '-1.0000')
         do while (ix.gt.0)
            bufout(1)(ix:ix+6) = '      x'
            ix = index(bufout(1), '-1.0000')
         enddo
         call write_log(1, bufout)
      enddo

      call write_log('basisfunctions B_{j,4}:')
      do ip = 1, npnt
                                                                ! TODO: j=1:15
         write(bufout,'(f5.1, i5,1x, 15f9.4)') si(ip), jseg(ip), (-1d0, j=1,jseg(ip)-4),                &
                                               (b4(ip,j), j=1,4), (-1d0, j=jseg(ip)+1,nknot-4)
         ix = index(bufout(1), '-1.0000')
         do while (ix.gt.0)
            bufout(1)(ix:ix+6) = '      x'
            ix = index(bufout(1), '-1.0000')
         enddo
         call write_log(1, bufout)
      enddo
   endif

end subroutine bspline_eval1d

!------------------------------------------------------------------------------------------------------------

subroutine bspline_make_ppform(nknot, tj, tiny_dt, nbrk, jseg, b1, b2, b3, b4, nspl, coef_xyz, spl,     &
                        my_ierror)
!--function: create PP-form of B-spline, evaluating {ay0-ay3}, {az0-az3} at start of each segment
   implicit none
!--subroutine arguments:
   integer,      intent(in)           :: nknot, nbrk, nspl
   integer,      intent(in)           :: jseg(nbrk)
   real(kind=8), intent(in)           :: tj(nknot), tiny_dt,                                            &
                                         b1(nbrk,1), b2(nbrk,2), b3(nbrk,3), b4(nbrk,4),                &
                                         coef_xyz(nspl,3)
   type(t_spline)                     :: spl
   integer,      intent(out)          :: my_ierror
!--local variables:
   integer              :: j, jj, ip, istat
   real(kind=8)         :: dtj_k1, dtj_k2, dtj_k3
   real(kind=8), dimension(:),   allocatable :: dtj_inv1, dtj_inv2, dtj_inv3
   real(kind=8), dimension(:,:), allocatable :: d3coef, d23coef, d123coef
   character(len=20)    :: namvar
   character(len=256)   :: errmsg

   my_ierror = 0

   namvar = 'dtj_inv1(pp)'
   allocate(dtj_inv1(nknot-1), dtj_inv2(nknot-2), dtj_inv3(nknot-3), stat=istat, errmsg=errmsg)
   if (istat.ne.0) goto 99

   ! determine inverse of average step sizes over 1/2/3 adjacent intervals, zero at zero step size

   do j = 1, nknot-1
      dtj_k1      =  tj(j+1) - tj(j)
      dtj_inv1(j) = 0d0
      if (dtj_k1.gt.tiny_dt) dtj_inv1(j) = 1d0 / dtj_k1
   enddo
   do j = 1, nknot-2
      dtj_k2      = (tj(j+2) - tj(j)) / 2d0
      dtj_inv2(j) = 0d0
      if (dtj_k2.gt.tiny_dt) dtj_inv2(j) = 1d0 / dtj_k2
   enddo
   do j = 1, nknot-3
      dtj_k3      = (tj(j+3) - tj(j)) / 3d0
      dtj_inv3(j) = 0d0
      if (dtj_k3.gt.tiny_dt) dtj_inv3(j) = 1d0 / dtj_k3
   enddo

   ! compute PP-spline coefficients [ ay0, az0 ] = B4 * [ coef_x, coef_y, coef_z ]

   if (spl%has_xdata) spl%ax0(1:nbrk) = 0d0
   spl%ay0(1:nbrk) = 0d0
   spl%az0(1:nbrk) = 0d0

   do ip = 1, nbrk
      do jj = 1, 4
         j = jseg(ip) + jj - 4
         if (spl%has_xdata) spl%ax0(ip) = spl%ax0(ip) + b4(ip,jj) * coef_xyz(j,1)
         spl%ay0(ip) = spl%ay0(ip) + b4(ip,jj) * coef_xyz(j,2)
         spl%az0(ip) = spl%az0(ip) + b4(ip,jj) * coef_xyz(j,3)
      enddo
   enddo

   ! compute B-spline coefficients [ d3coef ] := D3 * [ coef_x, coef_y, coef_z ]

   namvar = 'd3coef'
   allocate(d3coef(nspl+1,3), d23coef(nspl+2,3), d123coef(nspl+3,3), stat=istat, errmsg=errmsg)
   if (istat.ne.0) goto 99

   do j = 2, nspl
      d3coef(j,1) = dtj_inv3(j) * (coef_xyz(j,1) - coef_xyz(j-1,1))
      d3coef(j,2) = dtj_inv3(j) * (coef_xyz(j,2) - coef_xyz(j-1,2))
      d3coef(j,3) = dtj_inv3(j) * (coef_xyz(j,3) - coef_xyz(j-1,3))
      ! if (j.ge.218 .and. j.le.220) then
      !    write(bufout,'(a,i4,4(a,g12.4))') '  d3coef(',j,') = ',dtj_inv3(j),' * (', coef_xyz(j,2),      &
      !                   ' - ',coef_xyz(j-1,2),' ) =', d3coef(j,2)
      !    call write_log(1, bufout)
      ! endif
   enddo
   j = nspl+1
   d3coef(j,1) = dtj_inv3(j) * (0d0 - coef_xyz(j-1,1))
   d3coef(j,2) = dtj_inv3(j) * (0d0 - coef_xyz(j-1,2))
   d3coef(j,3) = dtj_inv3(j) * (0d0 - coef_xyz(j-1,3))

   ! compute PP-spline coefficients [ ay1, az1 ] = B3 * D3 * [ coef_x, coef_y, coef_z ]

   if (spl%has_xdata) spl%ax1(1:nbrk) = 0d0
   spl%ay1(1:nbrk) = 0d0
   spl%az1(1:nbrk) = 0d0

   do ip = 1, nbrk
      do jj = 1, 3
         j = jseg(ip) + jj - 3
         if (spl%has_xdata) spl%ax1(ip) = spl%ax1(ip) + b3(ip,jj) * d3coef(j,1)
         spl%ay1(ip) = spl%ay1(ip) + b3(ip,jj) * d3coef(j,2)
         spl%az1(ip) = spl%az1(ip) + b3(ip,jj) * d3coef(j,3)
      enddo
   enddo

   ! compute B-spline coefficients [ d23coef ] := D2 * D3 * [ coef_x, coef_y, coef_z ]

   do j = 2, nspl+1
      d23coef(j,1) = dtj_inv2(j) * (d3coef(j,1) - d3coef(j-1,1))
      d23coef(j,2) = dtj_inv2(j) * (d3coef(j,2) - d3coef(j-1,2))
      d23coef(j,3) = dtj_inv2(j) * (d3coef(j,3) - d3coef(j-1,3))
      ! if (j.eq.219 .or. j.eq.220) then
      !    write(bufout,'(a,i4,4(a,g12.4))') ' d23coef(',j,') = ',dtj_inv2(j),' * (', d3coef(j,2),      &
      !                   ' - ',d3coef(j-1,2),' ) =', d23coef(j,2)
      !    call write_log(1, bufout)
      ! endif
   enddo
   j = nspl+2
   d23coef(j,1) = dtj_inv2(j) * (0d0 - d3coef(j-1,1))
   d23coef(j,2) = dtj_inv2(j) * (0d0 - d3coef(j-1,2))
   d23coef(j,3) = dtj_inv2(j) * (0d0 - d3coef(j-1,3))

   ! compute PP-spline coefficients [ ay2, az2 ] = B2 * D2 * D3 * [ coef_x, coef_y, coef_z ] / 2

   if (spl%has_xdata) spl%ax2(1:nbrk) = 0d0
   spl%ay2(1:nbrk) = 0d0
   spl%az2(1:nbrk) = 0d0

   do ip = 1, nbrk
      do jj = 1, 2
         j = jseg(ip) + jj - 2
         if (spl%has_xdata) spl%ax2(ip) = spl%ax2(ip) + b2(ip,jj) * d23coef(j,1) / 2d0
         spl%ay2(ip) = spl%ay2(ip) + b2(ip,jj) * d23coef(j,2) / 2d0
         spl%az2(ip) = spl%az2(ip) + b2(ip,jj) * d23coef(j,3) / 2d0
      enddo
   enddo

   ! compute B-spline coefficients [ d123coef ] := D1 * D2 * D3 * [ coef_x, coef_y, coef_z ]

   do j = 2, nspl+2
      d123coef(j,1) = dtj_inv1(j) * (d23coef(j,1) - d23coef(j-1,1))
      d123coef(j,2) = dtj_inv1(j) * (d23coef(j,2) - d23coef(j-1,2))
      d123coef(j,3) = dtj_inv1(j) * (d23coef(j,3) - d23coef(j-1,3))
      ! if (j.eq.220) then
      !    write(bufout,'(a,i4,4(a,g12.4))') 'd123coef(',j,') = ',dtj_inv1(j),' * (', d23coef(j,2),      &
      !                   ' - ',d23coef(j-1,2),' ) =', d123coef(j,2)
      !    call write_log(1, bufout)
      ! endif
   enddo
   j = nspl+3
   d123coef(j,1) = dtj_inv1(j) * (0d0 - d23coef(j-1,1))
   d123coef(j,2) = dtj_inv1(j) * (0d0 - d23coef(j-1,2))
   d123coef(j,3) = dtj_inv1(j) * (0d0 - d23coef(j-1,3))

   ! compute PP-spline coefficients [ ay3, az3 ] = B1 * D1 * D2 * D3 * [ coef_x, coef_y, coef_z ] / 6

   if (spl%has_xdata) spl%ax3(1:nbrk) = 0d0
   spl%ay3(1:nbrk) = 0d0
   spl%az3(1:nbrk) = 0d0

   do ip = 1, nbrk
      do jj = 1, 1
         j = jseg(ip) + jj - 1
         if (spl%has_xdata) spl%ax3(ip) = spl%ax3(ip) + b1(ip,jj) * d123coef(j,1) / 6d0
         spl%ay3(ip) = spl%ay3(ip) + b1(ip,jj) * d123coef(j,2) / 6d0
         spl%az3(ip) = spl%az3(ip) + b1(ip,jj) * d123coef(j,3) / 6d0
         ! if (ip.eq.205) then
         !    write(bufout,'(2(a,i4),3(a,g12.4))') 'j=',j,', ip=',ip,': ay3 =',b1(ip,jj),' * ',           &
         !                d123coef(j,2),' / 6d0 = ',spl%ay3(ip)
         !    call write_log(1, bufout)
         ! endif
      enddo
   enddo

   if (ldebug.ge.6) call spline_print(spl, 'pp_form', 5)

   return

   ! error handling for memory allocation

99 continue

   write(bufout,'(3a,i4,2a)') ' ERROR in memory allocation (',trim(namvar),':',istat,'): ',trim(errmsg)
   call write_log(1, bufout)
   call abort_run()

end subroutine bspline_make_ppform

!------------------------------------------------------------------------------------------------------------

subroutine bspline_check_kink_accel(nprf, nkink, ikinks, naccel, iaccel, my_ierror)
!--function: check requirements of B-spline method on kinks and accelerations
   implicit none
!--subroutine arguments:
   integer,      intent(in)     :: nprf, nkink, naccel
   integer,      intent(in)     :: ikinks(nkink), iaccel(naccel)
   integer,      intent(out)    :: my_ierror
!--local variables:
   integer                      :: ia, ik, itot, idist, ntot, istat
   character(len=20)            :: namvar
   character(len=256)           :: errmsg
   integer, dimension(:,:), allocatable :: kink_acc

   my_ierror = 0

   ! require that kinks and accelerations are sorted in strictly increasing order

   do ik = 2, nkink
      if (ikinks(ik).le.ikinks(ik-1)) then
         my_ierror = 2001
         write(bufout,'(a)') ' ERROR: B-spline: array ikinks must be strictly increasing'
         call write_log(1, bufout)
         write(bufout,'(3(a,i4))') '        pos',ik,': ikinks=',ikinks(ik-1),' >= ',ikinks(ik)
         call write_log(1, bufout)
      endif
   enddo

   do ia = 2, naccel
      if (iaccel(ia).le.iaccel(ia-1)) then
         my_ierror = 2002
         write(bufout,'(a)') ' ERROR: B-spline: array iaccel must be strictly increasing'
         call write_log(1, bufout)
         write(bufout,'(3(a,i4))') '        pos',ia,': iaccel=',iaccel(ia-1),' >= ',iaccel(ia)
         call write_log(1, bufout)
      endif
   enddo

   ! require that kinks lie within measurement range
   ! Note that ikinks includes boundaries [1, nprf]

   if (nkink.ge.3) then
      if (ikinks(2).le.1 .or. ikinks(2).eq.3 .or. ikinks(nkink-1).eq.nprf-2 .or.                        &
                                                                        ikinks(nkink-1).ge.nprf) then
         my_ierror = 2003
         write(bufout,'(2(a,i4),a)') ' ERROR: B-spline: kinks must be either 2 or nprf-1 =',nprf-1,     &
                   ', or 4 <= ikink <= nprf-3 =',nprf-3
         call write_log(1, bufout)
         write(bufout,'(2(a,i4))') '        1st =',ikinks(2),', last =',ikinks(nkink-1)
         call write_log(1, bufout)
      endif
   endif

   ! require that accelerations lie within measurement range, away from boundaries

   if (naccel.ge.1) then
      if (iaccel(1).le.3 .or. iaccel(naccel).ge.nprf-2) then
         my_ierror = 2004
         write(bufout,'(a,i4)') ' ERROR: B-spline: accelerations must be 4 <= iaccel <= nprf-3 =',nprf-3
         call write_log(1, bufout)
         write(bufout,'(2(a,i4))') '        1st =',iaccel(1),', last =',iaccel(naccel)
         call write_log(1, bufout)
      endif
   endif

   ! create sorted list of boundaries, kinks and accelerations
   ! note that ikinks includes boundaries [1, nprf]

   ntot = nkink + naccel
   namvar = 'kink_acc'
   allocate(kink_acc(ntot,2), stat=istat, errmsg=errmsg)
   if (istat.ne.0) goto 99

   kink_acc(1,1) = 1    ! position
   kink_acc(1,2) = 0    ! type: 0=boundary, 1=kink, 2=accel
   itot = 1  ! last used position
   ik = 2    ! next unused position
   ia = 1    ! next unused position

   do while(ik.le.nkink-1 .and. ia.le.naccel)
      itot = itot + 1
      if (ikinks(ik).le.iaccel(ia)) then
         kink_acc(itot,1) = ikinks(ik)
         kink_acc(itot,2) = 1
         ik = ik + 1
      else
         kink_acc(itot,1) = iaccel(ia)
         kink_acc(itot,2) = 2
         ia = ia + 1
      endif
   enddo

   do while(ik.le.nkink-1)
      itot = itot + 1
      kink_acc(itot,1) = ikinks(ik)
      kink_acc(itot,2) = 1
      ik = ik + 1
   enddo
         
   do while(ia.le.naccel)
      itot = itot + 1
      kink_acc(itot,1) = iaccel(ia)
      kink_acc(itot,2) = 2
      ia = ia + 1
   enddo

   itot = itot + 1
   kink_acc(itot,1) = nprf
   kink_acc(itot,2) = 0

   ! require that kinks have 0 or >= 2 points between them, kinks and accel have >= 2 points between them

   do itot = 2, ntot

      idist = kink_acc(itot,1) - kink_acc(itot-1,1)

      if (idist.eq.2 .and. kink_acc(itot-1,2).eq.1 .and. kink_acc(itot,2).eq.1) then
         my_ierror = 2005
         write(bufout,'(a)') ' ERROR: B-spline: kinks must be adjacent or have >= 2 points between them'
         call write_log(1, bufout)
         write(bufout,'(2(a,i4,a,i1),a)') '       1st =',kink_acc(itot-1,1),' (',kink_acc(itot-1,2),    &
                        '), 2nd =',kink_acc(itot,1),' (',kink_acc(itot,2),')'
         call write_log(1, bufout)
      endif

      if (idist.le.2 .and. (kink_acc(itot-1,2).eq.2 .or. kink_acc(itot,2).eq.2)) then
         my_ierror = 2006
         write(bufout,'(2a,/,a)') ' ERROR: B-spline: accelerations must have at least 2 points between', &
                ' them ', '                 and kinks, boundaries and other accelerations'
         call write_log(2, bufout)
         write(bufout,'(2(a,i4,a,i1),a)') '       1st =',kink_acc(itot-1,1),' (',kink_acc(itot-1,2),    &
                        '), 2nd =',kink_acc(itot,1),' (',kink_acc(itot,2),')'
         call write_log(1, bufout)
      endif

   enddo ! itot

   return

   ! error handling for memory allocation

99 continue

   write(bufout,'(3a,i4,2a)') ' ERROR in memory allocation (',trim(namvar),':',istat,'): ',trim(errmsg)
   call write_log(1, bufout)
   call abort_run()

end subroutine bspline_check_kink_accel

!------------------------------------------------------------------------------------------------------------

subroutine bspline_check_mask(nmeasu, nmeasv, mask, ierror)
!--function: check requirements on `missing parts' encoded in mask array
   implicit none
!--subroutine arguments:
   integer                            :: nmeasu, nmeasv
   logical,      intent(in)           :: mask(nmeasu,nmeasv)
   integer,      intent(out)          :: ierror
!--local variables:
   integer, parameter :: minsep = 4, minpts = 4, maxreg = 9
   logical     :: in_region, is_ok
   integer     :: iu, jv, nregion, npts

   do iu = 1, nmeasu

      ! determine the number of regions of active points
      !  - each region must hold at least minpts points
      !  - there must be at least minsep points between regions

      in_region = .false.
      nregion = 0
      npts    = 0

      do jv = 1, nmeasv
         if (.not.in_region .and. mask(iu,jv)) then

            ! entering a new active region. Check previous inactive part

            is_ok = (npts.ge.minsep .or. (nregion.eq.0 .and. npts.eq.0))
            if (.not.is_ok) then
               ierror = 171
               write(bufout,'(a,i4,a,i2,2a,i2,a)') ' slice',iu,' has',npts,' inactive points at start ', &
                        'or between active regions, need separation >=',minsep,' points'
               call write_log(1, bufout)
            endif

            ! set up new active region

            in_region = .true.
            nregion   = nregion + 1
            npts      = 1

         elseif (in_region .and. .not.mask(iu,jv)) then

            ! entering a new inactive part. Check previous active region

            if (npts.lt.minpts) then
               ierror = 172
               write(bufout,'(a,i4,a,i2,2a,i2,a)') ' slice',iu,' has a region of',npts,' active points,', &
                        ' must be >=',minpts,' points'
               call write_log(1, bufout)
            endif

            ! set up new inactive part

            in_region = .false.
            npts      = 1

         else

            ! continue with existing inactive part or active region

            npts      = npts + 1

         endif
      enddo ! points jv within slice

      if (nregion.gt.maxreg) then
         ierror = 173
         write(bufout,'(a,i4,a,i2,2a,i2,a)') ' slice',iu,' has',nregion,' regions of active points,',    &
                  ' whereas max.',maxreg,' regions are permitted.'
         call write_log(1, bufout)
      endif

   enddo ! slices iu

   ! check the active regions per column ('interpolation path')

   do jv = 1, nmeasv

      ! determine the number of regions of active points
      !  - each region must hold at least minpts points
      !  - there must be at least minsep points between regions

      in_region = .false.
      nregion = 0
      npts    = 0

      do iu = 1, nmeasu
         if (.not.in_region .and. mask(iu,jv)) then

            ! entering a new active region. Check previous inactive part

            is_ok = (npts.ge.minsep .or. (nregion.eq.0 .and. npts.eq.0))
            if (.not.is_ok) then
               ierror = 171
               write(bufout,'(a,i4,a,i2,2a,i2,a)') ' intpol.path',jv,' has',npts,' inactive points ', &
                        'at start or between active regions, need separation >=',minsep,' points'
               call write_log(1, bufout)
            endif

            ! set up new active region

            in_region = .true.
            nregion   = nregion + 1
            npts      = 1

         elseif (in_region .and. .not.mask(iu,jv)) then

            ! entering a new inactive part. Check previous active region

            if (npts.lt.minpts) then
               ierror = 172
               write(bufout,'(a,i4,a,i2,2a,i2,a)') ' intpol.path',jv,' has a region of',npts,' active', &
                        ' points, must be >=',minpts,' points'
               call write_log(1, bufout)
            endif

            ! set up new inactive part

            in_region = .false.
            npts      = 1

         else

            ! continue with existing inactive part or active region

            npts      = npts + 1

         endif
      enddo ! points jv within slice

      if (nregion.gt.maxreg) then
         ierror = 173
         write(bufout,'(a,i4,a,i2,2a,i2,a)') ' intpol.path',jv,' has',nregion,' regions of active',    &
                  ' points, whereas max.',maxreg,' regions are permitted.'
         call write_log(1, bufout)
      endif

   enddo ! paths jv

end subroutine bspline_check_mask

!------------------------------------------------------------------------------------------------------------

end module m_bspline_def
