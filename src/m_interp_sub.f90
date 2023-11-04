!------------------------------------------------------------------------------------------------------------
! m_interp_sub - implementation of linear 1D interpolation, bilinear 2D, bicubic 2D interpolation
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
submodule (m_interp) m_interp_sub
   use m_globals
   use m_markers
   implicit none

!------------------------------------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------------------------------------

module subroutine interp1d_set_debug(new_ldebug, new_ii_debug, new_iel_debug)
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

end subroutine interp1d_set_debug

!------------------------------------------------------------------------------------------------------------

module function is_strictly_monotonic( nin, xin )
!--function: determine whether array values are strictly increasing or decreasing
   implicit none
!--result value
   logical  :: is_strictly_monotonic
!--subroutine arguments:
   integer,      intent(in)  :: nin                   ! number of points
   real(kind=8), intent(in)  :: xin(nin)              ! input array
!--local variables:
   logical      :: is_mono, has_increasing, has_equal, has_decreasing
   integer      :: i

   if (nin.le.1) then

      is_mono = .true.

   else

      has_increasing = .false.
      has_equal      = .false.
      has_decreasing = .false.
      is_mono        = .not.has_equal .and. (.not.has_increasing .or. .not.has_decreasing)

      i = 1
      do while(i.lt.nin .and. is_mono)
         i = i + 1
         if (xin(i).gt.xin(i-1)) then
            has_increasing = .true.
         elseif (xin(i).lt.xin(i-1)) then
            has_decreasing = .true.
         else
            has_equal      = .true.
         endif
         is_mono = .not.has_equal .and. (.not.has_increasing .or. .not.has_decreasing)
      enddo

   endif

   is_strictly_monotonic = is_mono

end function is_strictly_monotonic

!------------------------------------------------------------------------------------------------------------

module subroutine check_monotone(npnt, xpnt, istrict, ierror)
!--function: check whether array xpnt is monotonically increasing or decreasing
   implicit none
!--subroutine arguments:
   integer,      intent(in)   :: npnt, istrict
   real(kind=8), intent(in)   :: xpnt(npnt)
   integer,      intent(out)  :: ierror
!--local variables:
   integer              :: ipnt
   logical              :: is_ascending
   real(kind=8)         :: tiny = 1d-10

   ! check that xpnt is sorted in increasing or decreasing order

   ierror = 0
   is_ascending = (xpnt(npnt).gt.xpnt(1))

   if (is_ascending) then

      ipnt = 1
      do while(ipnt.lt.npnt .and. ierror.eq.0)
         ipnt = ipnt + 1
         if (istrict.ge.1) then
            if (xpnt(ipnt).le.xpnt(ipnt-1)) ierror = ipnt
         else
            if (xpnt(ipnt).lt.xpnt(ipnt-1)-tiny) ierror = ipnt
         endif
      enddo

   else

      ipnt = 1
      do while(ipnt.lt.npnt .and. ierror.eq.0)
         ipnt = ipnt + 1
         if (istrict.ge.1) then
            if (xpnt(ipnt).ge.xpnt(ipnt-1)) ierror = ipnt
         else
            if (xpnt(ipnt).gt.xpnt(ipnt-1)+tiny) ierror = ipnt
         endif
      enddo

   endif

end subroutine check_monotone

!------------------------------------------------------------------------------------------------------------

module subroutine locate_segment( np, vp, vi, iseg )
!--function: find segment 'iseg' in list of segments { vp(ip), ip=1,np } that contains position vi,
!            starting with input iseg as initial estimate. 
!            On output: iseg in { 0, .., np },  iseg=0 when vi < min(vp), iseg=np when vi >= max(vp).
!            cf. hunt-algorithm of Numerical Recipes, sec. 3.4.
   implicit none
!--subroutine arguments
   integer,      intent(in)    :: np
   real(kind=8), intent(in)    :: vp(np), vi
   integer,      intent(inout) :: iseg
!--local variables
   logical      :: ascending
   integer      :: il, iu, inc

   ! assuming np > 0, sorted vp(i) strictly increasing or strictly decreasing

   ascending = (vp(np).ge.vp(1))

   if (iel_debug.ge.1) then
      write(bufout,*) 'np=',np,', v(1)=',vp(1),', v(end)=',vp(np),', ascending=',ascending
      call write_log(1, bufout)
   endif

   ! if search position lies before start of data: done

   if ( (vi.lt.vp(1)) .eqv. ascending ) then

      iseg = 0
      if (iel_debug.ge.1) then
         if (ascending) then
            write(bufout,*) 'ascending, vi=',vi,'< v(1) ==> iseg=0'
         else
            write(bufout,*) 'descending, vi=',vi,'> v(1) ==> iseg=0'
         endif
         call write_log(1, bufout)
      endif

   ! elseif search position lies after end of data: done

   elseif ( (vi.gt.vp(np)) .eqv. ascending ) then

      iseg = np
      if (iel_debug.ge.1) then
         if (ascending) then
            write(bufout,*) 'ascending, v(end) < vi=',vi,' ==> iseg=', iseg
         else
            write(bufout,*) 'descending, v(end) > vi=',vi,' ==> iseg=', iseg
         endif
         call write_log(1, bufout)
      endif

   ! elseif search position equals the end-point: done

   elseif ( vi.eq.vp(np) ) then

      iseg = np - 1
      if (iel_debug.ge.1) then
         write(bufout,*) 'vi=',vi,'== v(end) ==> iseg=', iseg
         call write_log(1, bufout)
      endif

   else

      ! hunting: find a bracket [il, iu], starting with iseg as initial estimate

      if (iseg.le.0 .or. iseg.gt.np) then

         ! initial estimate is out of range, use whole range as bracket

         il = 0
         iu = np + 1

      else

         il = iseg
         iu = iseg

         ! vp(il) <= vi: hunting upwards: push vp(iu) to the right, beyond vi

         if ( (vi.ge.vp(il)) .eqv. ascending ) then

            inc = 1
            do while( iu.lt.np .and. (vi.ge.vp(iu) .eqv. ascending) )
               iu = min(iu+inc, np)
               inc = inc * 2
            enddo

         ! vi <  vp(iu): hunting downwards: push vp(il) to the left, beyond vi

         else

            inc = 1
            do while( il.gt.0 .and. (vi.lt.vp(il) .eqv. ascending) )
               il = max(il-inc, 1)
               inc = inc * 2
            enddo
         endif

      endif

      ! locating: shrink the bracket until just two points remain

      do while(iu-il.gt.1)

         iseg = (iu + il) / 2
         if (iel_debug.ge.1) then
            write(bufout,*) 'il=',il,', iu=',iu,', iseg=',iseg,', v=',vp(iseg)
            call write_log(1, bufout)
         endif
         if ( (vi.ge.vp(iseg)) .eqv. ascending ) then
            il = iseg
         else
            iu = iseg
         endif
      enddo

      ! select final value

      iseg = il
      if (iel_debug.ge.1) then
         write(bufout,*) 'il=',il,', iu=',iu,', iseg=',iseg,', v=',vp(iseg)
         call write_log(1, bufout)
      endif

   endif

end subroutine locate_segment

!------------------------------------------------------------------------------------------------------------

module subroutine locate_interval(np, vp, v_low, v_hig, ilow, ihig)
!--purpose: determine ilow, ihig such that vp( [ilow:ihig] ) just encompasses [v_low, v_hig]
!           return ilow=1 in case v_low < vp(1) and ihig=np in case v_hig > vp(np)
   implicit none
!--subroutine parameters:
   integer,      intent(in)  :: np
   real(kind=8), intent(in)  :: vp(np), v_low, v_hig
   integer,      intent(out) :: ilow, ihig
!--local variables:
   logical, parameter :: use_linear = .false.
   integer            :: jlow, jhig
   logical            :: check_next

   ! assuming np > 0, vp(i) sorted and strictly increasing

   ilow = nint( (v_low - vp(1)) * (np - 1) / (vp(np) - vp(1)) ) ! initial estimate assuming equal dv
   call locate_segment( np, vp, v_low, ilow )                   ! returns ilow=0 for v_low < vp(1)
   ilow = max(1, ilow)

   ihig = nint( (v_hig - vp(1)) * (np - 1) / (vp(np) - vp(1)) ) ! initial estimate assuming equal dv
   call locate_segment( np, vp, v_hig, ihig )                   ! returns ihig<=np : vp(ihig) <= v_hig
   ihig = min(np, ihig+1)

   if (use_linear) then

      ! check computation using brute-force linear search

      ! determine the first index to keep, jlow = find(vp <= v_low, 1, 'last')

      jlow = 1
      check_next = (jlow.lt.np)
      do while(check_next)
         if (vp(jlow+1).le.v_low) then
            jlow = jlow + 1
            check_next = (jlow.lt.np)
         else
            check_next = .false.
         endif
      enddo

      ! determine the last index to keep, jhig = find(vp >= v_hig, 1, 'first')

      jhig = jlow
      check_next = (jhig.lt.np)
      do while(check_next)
         if (vp(jhig).lt.v_hig) then
            jhig = jhig + 1
            check_next = (jhig.lt.np)
         else
            check_next = .false.
         endif
      enddo

      if (ilow.ne.jlow .or. ihig.ne.jhig .or. ldebug.ge.2) then
         write(bufout,'(2(a,f8.3),4(a,i4),a)') ' locate_interval: v=[',v_low,',',v_hig,']: range i=[',  &
                ilow, ',',ihig,'], range j=[',jlow,',',jhig,']'
         call write_log(1, bufout)
      endif

   endif

end subroutine locate_interval

!------------------------------------------------------------------------------------------------------------

module subroutine interp_1d_to_2d( nin, xin, fin, nout1, nout2, xintp, fintp, ierror, exterval)
!--purpose: interpolate (non equidistant) 1d input data [xin, fin] to the points listed in xintp,
!           which can be 2d. Using linear interpolation, constant extrapolation of first/last values.
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: nin, nout1, nout2     ! number of points in, out
   real(kind=8), intent(in)  :: xin(nin), fin(nin)    ! input data [x,fx]
   real(kind=8), intent(in)  :: xintp(nout1, nout2)   ! output locations [x]
   real(kind=8), intent(out) :: fintp(nout1, nout2)   ! output values [fx]
   integer,      intent(out) :: ierror
   real(kind=8), optional    :: exterval
!--local variables:
   logical      :: has_exter, ascending
   integer      :: iin, iout, jout
   real(kind=8) :: dx, df, fac

   ierror    = 0
   has_exter = (present(exterval))

   ! facilitate decreasing/increasing x-data

   ascending = xin(nin).gt.xin(1)

   if (ldebug.ge.10) then
      write(bufout,'(3(a,i4))') ' interp_1d: starting with nin=',nin,', nout=',nout1,' x',nout2
      call write_log(1, bufout)
      if (     ascending) call write_log('            input values in increasing order')
      if (.not.ascending) call write_log('            input values in decreasing order')

      if (     has_exter) then
         write(bufout,'(a,f10.1)') '            using default value',exterval
         call write_log(1, bufout)
      endif
      if (.not.has_exter) call write_log('            using constant extrapolation')
   endif
   if (ldebug.ge.1) then
      if (.not.is_strictly_monotonic(nin, xin)) then
         call write_log(' WARNING: interp_1d: input values are NOT strictly monotonic')
      endif
   endif

   iin = 0

   do jout = 1, nout2
      do iout = 1, nout1

         ! find interval iin = [ xin(iin), xin(iin+1) ] containing xintp(i,j)
         ! using previous iin as initial estimate

         call locate_segment( nin, xin, xintp(iout,jout), iin )

         if (iin.le.0) then

            ! point xintp lies before the input-range
            !  - ascending:  xintp(i,j) <= xin(1) < xin(end)
            !  - descending: xintp(i,j) >= xin(1) > xin(end)

            if (has_exter) then
               fintp(iout,jout) = exterval      ! prescribed default value
            else
               fintp(iout,jout) = fin(1)        ! constant extrapolation
            endif

            if (ldebug.ge.15) then
               write(bufout,'(2(a,i5),a)') ' iout=',iout,',',jout,': below xin(1)'
               call write_log(1, bufout)
            endif

         elseif (iin.ge.nin) then

            ! point xintp lies after the input range
            !  - ascending:   xin(1) < xin(end) <= xintp(i,j)
            !  - descending:  xin(1) > xin(end) >= xintp(i,j)

            if (has_exter) then
               fintp(iout,jout) = exterval      ! prescribed default value
            else
               fintp(iout,jout) = fin(nin)      ! constant extrapolation
            endif

            if (ldebug.ge.15) then
               write(bufout,'(2(a,i5),a)') ' iout=',iout,',',jout,': above xin(end)'
               call write_log(1, bufout)
            endif

         else

            ! perform linear interpolation between iin and iin+1

            if (ascending) then
               dx = max( 1d-9, xin(iin+1) - xin(iin))
            else
               dx = min(-1d-9, xin(iin+1) - xin(iin))
            endif

            df = fin(iin+1) - fin(iin)
            fac = (xintp(iout,jout) - xin(iin)) / dx
            fintp(iout,jout) = (1d0-fac) * fin(iin) + fac * fin(iin+1)

         endif

      enddo ! iout
   enddo ! jout
         
end subroutine interp_1d_to_2d

!------------------------------------------------------------------------------------------------------------

module subroutine interp_1d_to_1d( nin, xin, fin, nout1, xintp, fintp, my_ierror, exterval)
!--purpose: interpolate (non equidistant) 1d input data [xin, fin] to the points listed in xintp(:).
!           Using linear interpolation, constant extrapolation of first/last values.
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: nin, nout1            ! number of points in, out
   real(kind=8), intent(in)  :: xin(nin), fin(nin)    ! input data [x,fx]
   real(kind=8), intent(in)  :: xintp(nout1)          ! output locations [x]
   real(kind=8), intent(out) :: fintp(nout1)          ! output values [fx]
   integer,      intent(out) :: my_ierror
   real(kind=8), optional    :: exterval
!--local variables:
   logical      :: has_exter
   integer      :: nout2, sub_ierror

   my_ierror = 0

   nout2 = 1
   has_exter = (present(exterval))

   if (has_exter) then
      call interp_1d_to_2d(nin, xin, fin, nout1, nout2, xintp, fintp, sub_ierror, exterval)
   else
      call interp_1d_to_2d(nin, xin, fin, nout1, nout2, xintp, fintp, sub_ierror)
   endif
   my_ierror = sub_ierror

end subroutine interp_1d_to_1d

!------------------------------------------------------------------------------------------------------------

module subroutine interp_1d_to_scalar( nin, xin, fin, xintp, fintp, my_ierror, exterval)
!--purpose: interpolate (non equidistant) 1d input data [xin, fin] to the point xintp.
!           Using linear interpolation, constant extrapolation of first/last values.
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: nin                   ! number of points in
   real(kind=8), intent(in)  :: xin(nin), fin(nin)    ! input data [x,fx]
   real(kind=8), intent(in)  :: xintp                 ! output location [x]
   real(kind=8), intent(out) :: fintp                 ! output value [fx]
   integer,      intent(out) :: my_ierror
   real(kind=8), optional    :: exterval
!--local variables:
   logical      :: has_exter
   integer      :: nout1, nout2, sub_ierror
   real(kind=8) :: xarr(1), farr(1)

   my_ierror = 0
   nout1   = 1
   nout2   = 1
   xarr(1) = xintp
   has_exter = (present(exterval))

   if (has_exter) then
      call interp_1d_to_2d(nin, xin, fin, nout1, nout2, xarr, farr, sub_ierror, exterval)
   else
      call interp_1d_to_2d(nin, xin, fin, nout1, nout2, xarr, farr, sub_ierror)
   endif
   my_ierror = sub_ierror

   fintp = farr(1)

end subroutine interp_1d_to_scalar

!------------------------------------------------------------------------------------------------------------

module subroutine interp_wgt_surf2unif(nnode_x, nnode_y, nnode, x_node, y_node, z_node, z_thrs, mx, my,        &
                                nout, dx_arg, dy_arg, x_out, y_out, ii2iel, ii2nod, wii2nod,            &
                                fac_uv, my_ierror)
!--function: construct interpolation tables between a curvi-linear input grid (surface) and a uniform
!            output-grid.
!            Returns weighting factors for bilinear and relative positions for bicubic interpolation
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: nnode_x, nnode_y, nnode ! number of points in the input grid
   real(kind=8), intent(in)  :: x_node(nnode), y_node(nnode), z_node(nnode)
                                                  ! (x,y,z) coordinates of input grid points
   real(kind=8), intent(in)  :: z_thrs            ! typical dimension in z-direction, threshold for bounding box
   integer,      intent(in)  :: mx, my, nout      ! number of points in the output grid
   real(kind=8), intent(in)  :: dx_arg, dy_arg    ! step-sizes used in the output grid, <0 for non-uniform
   real(kind=8), intent(in)  :: x_out(nout), y_out(nout) ! coordinates of the output grid points
   integer,      intent(out) :: ii2iel(nout)      ! input element number for output grid points
   integer,      intent(out) :: ii2nod(4,nout)    ! input node numbers for output grid points
   real(kind=8), intent(out) :: wii2nod(4,nout)   ! interpolation weights per input node per output grid point
   real(kind=8), intent(out) :: fac_uv(2,nout)    ! relative position of output grid point in input element
                                                  ! used for bicubic interpolation
   integer,      intent(out) :: my_ierror
!--local variables:
!  character(len=*), parameter  :: subnam = 'interp_wgt_surf2unif'
   logical          :: l_nan, l_far, in_cell, is_uniform_x, is_uniform_y
   integer          :: nel_x, nel_y, nel, iel_x, iel_y, iel, i, i_ll, i_lr, i_ur, i_ul, inod,           &
                       ix, iy, ii, ix0, ix1, iy0, iy1, num_error_msg, my_ldebug, sub_ierror
   real(kind=8)     :: dx, dy, xout_low, yout_low, xout_hig, yout_hig, xin_min, xin_max,        &
                       yin_min, yin_max, zin_min, zin_max
   real(kind=8)     :: tlow_x, tlow_y, trgt_x, trgt_y, tupp_x, tupp_y, tlft_x, tlft_y,          &
                       nlow_x, nlow_y, nrgt_x, nrgt_y, nupp_x, nupp_y, nlft_x, nlft_y,          &
                       ofs_low, ofs_rgt, ofs_upp, ofs_lft, fac_long, fac_lat, rgt_sense
   type(t_vec)      :: vec_x, vec_y, vec_z
   integer,      dimension(:,:), allocatable :: iel2node
   real(kind=8), dimension(:),   allocatable :: xo_list, yo_list

   my_ierror = 0
   num_error_msg = 0

   ! determine whether the output grid is uniform/non-uniform in x and y

   is_uniform_x = (dx_arg.gt.0d0)
   is_uniform_y = (dy_arg.gt.0d0)

   if (.not.is_uniform_x) then
      allocate(xo_list(mx))
      do ix = 1, mx
         ii = ix
         xo_list(ix) = x_out(ii)
      enddo
      dx = (xo_list(mx) - xo_list(1)) / real(max(1,mx-1))
   else
      dx = dx_arg
   endif

   if (.not.is_uniform_y) then
      allocate(yo_list(my))
      do iy = 1, my
         ii = 1 + (iy-1)*mx
         yo_list(iy) = y_out(ii)
      enddo
      dy = (yo_list(my) - yo_list(1)) / real(max(1,my-1))
   else
      dy = dy_arg
   endif

   if (.not.is_uniform_x .and. mx.le.1 .and. my.gt.1) dx = dy
   if (.not.is_uniform_y .and. my.le.1 .and. mx.gt.1) dy = dx

   ! get extent of the uniform output grid considered as elements

   xout_low = x_out(1) - 0.5*dx
   xout_hig = xout_low + mx*dx
   yout_low = y_out(1) - 0.5*dy
   yout_hig = yout_low + my*dy

   ! administration of elements of the input grid
   ! note: we assume that node numbers are running fastest in the direction of nnode_x

   nel_x = nnode_x - 1
   nel_y = nnode_y - 1
   nel   = nel_x * nel_y

   if (nel.le.0) then
      if (ldebug.ge.-1 .and. num_error_msg.le.5) then
         write(bufout,'(2(a,i5),a)') ' Internal error: interp_wgt: input grid has',nel_x,' x',          &
                nel_y,' elements.'
         call write_log(1, bufout)
         num_error_msg = num_error_msg + 1
      endif
   endif

   if (ldebug.ge.3) then
      write(bufout,'(a,2(i5,a))') ' interp_wgt_surf2unif: input grid has ', nel_x,' x', nel_y,' elements'
      call write_log(1, bufout)
      write(bufout,'(a,2(i5,a))') '                       output grid has', mx,' x', my,' points'
      call write_log(1, bufout)

      if (ii_debug.ge.1 .and. ii_debug.ge.1 .and. ii_debug.le.nout) then
         write(bufout,'(a,i7,2(a,f10.3),a)') ' printing debug for  ii_debug=', ii_debug,' (x,y) = (',    &
                x_out(ii_debug),',', y_out(ii_debug),')'
         call write_log(1, bufout)
      endif
   endif

   if (ldebug.ge.5) then
      write(bufout,'(2(a,f8.4))') '                     output grid sizes',dx,' x',dy
      call write_log(1, bufout)
      write(bufout,'(4(a,f10.4))') '                 x=', xout_low,' to ', xout_hig,', y=',             &
                yout_low, ' to ', yout_hig
      call write_log(1, bufout)

      zin_min = z_node(1)
      zin_max = zin_min
      do inod = 1, nnode
         zin_min = min(zin_min, z_node(inod))
         zin_max = max(zin_max, z_node(inod))
      enddo
      write(bufout,'(17x,2(a,f10.4),a,g12.4)') 'z=', zin_min,' to ', zin_max,', z_thrs=', z_thrs
      call write_log(1, bufout)

      if (.not.is_uniform_x) then
         write(bufout,'(a,i4,a)') ' non-uniform x, mx=',mx,':'
         call write_log(1, bufout)
         do ix = 1, mx
            write(bufout,'(a,i4,a,f12.3)') ' x(',ix,') =',xo_list(ix)
            call write_log(1, bufout)
         enddo
      endif

      if (.not.is_uniform_y) then
         write(bufout,'(a,i4,a)') ' non-uniform y, my=',my,':'
         call write_log(1, bufout)
         do iy = 1, my
            write(bufout,'(a,i4,a,f12.3)') ' y(',iy,') =',yo_list(iy)
            call write_log(1, bufout)
         enddo
      endif
   endif


   if (ldebug.ge.3 .and. iel_debug.ge.1 .and. iel_debug.le.nel) then
      iel_x = mod(iel_debug-1, nel_x) + 1
      iel_y = (iel_debug-iel_x) / nel_x + 1
      write(bufout,'(a,i7,2(a,i4),a)') ' printing debug for iel_debug=', iel_debug,' (',iel_x,',',iel_y,')'
      call write_log(1, bufout)
   endif

   allocate(iel2node(4,nel))

   do iel_y = 1, nel_y
      do iel_x = 1, nel_x
         ! note: this is our own element numbering; doesn't have to be the same as the element numbering
         iel = iel_x + (iel_y-1)*nel_x
         ! note: "lower-left" refers to a "logical grid layout"; doesn't refer to CONTACT-coordinates
         iel2node(1,iel) = iel_x   + (iel_y-1)*nnode_x   ! node at lower-left corner
         iel2node(2,iel) = iel_x+1 + (iel_y-1)*nnode_x   ! node at lower-right corner
         iel2node(3,iel) = iel_x+1 + (iel_y  )*nnode_x   ! node at upper-right corner
         iel2node(4,iel) = iel_x   + (iel_y  )*nnode_x   ! node at upper-left corner

         if (ldebug.ge.3 .and. iel.eq.iel_debug) then
            write(bufout,'(2(a,i4),4(a,i7))') ' element (',iel_x,',',iel_y,'): nodes ll=',iel2node(1,iel), &
                        ', lr=',iel2node(2,iel), ', ur=',iel2node(3,iel), ', ul=',iel2node(4,iel)
            call write_log(1,bufout)

            i_ll = iel2node(1,iel)
            i_lr = iel2node(2,iel)
            i_ur = iel2node(3,iel)
            i_ul = iel2node(4,iel)
            write(bufout, 980) 'lower-left: ',i_ll, x_node(i_ll), y_node(i_ll)
            call write_log(1,bufout)
            write(bufout, 980) 'lower-right:',i_lr, x_node(i_lr), y_node(i_lr)
            call write_log(1,bufout)
            write(bufout, 980) 'upper-right:',i_ur, x_node(i_ur), y_node(i_ur)
            call write_log(1,bufout)
            write(bufout, 980) 'upper-left: ',i_ul, x_node(i_ul), y_node(i_ul)
            call write_log(1,bufout)
         endif
      enddo
   enddo

   ! determine the normal direction to input element 1, in order to define the "upper side" of the surface
   ! element 1 x-direction: lower-left --> lower-right
   ! element 1 y-direction: lower-left --> upper-left
   ! element 1 z-direction: cross-product, vec_z = vec_x x vec_y

   iel_x = max(1, nel_x/2)
   iel_y = max(1, nel_y/2)
   iel   = iel_x + (iel_y-1)*nel_x
   i_ll = iel2node(1,iel)
   i_lr = iel2node(2,iel)
   i_ul = iel2node(4,iel)

   vec_x%v(1:3) = (/ x_node(i_lr), y_node(i_lr), z_node(i_lr) /) -                                      &
                                                        (/ x_node(i_ll), y_node(i_ll), z_node(i_ll) /)
   vec_y%v(1:3) = (/ x_node(i_ul), y_node(i_ul), z_node(i_ul) /) -                                      &
                                                        (/ x_node(i_ll), y_node(i_ll), z_node(i_ll) /)
   vec_z = vec_x .cross. vec_y
   if (vec_z%z().ge.0d0) then
      rgt_sense =  1d0
   else
      rgt_sense = -1d0
   endif

   if (ldebug.ge.10) then
      call vec_print(vec_x, 'vec_x', 2)
      call vec_print(vec_y, 'vec_y', 2)
      call vec_print(vec_z, 'vec_z', 2)
      write(bufout,*) 'rgt_sense=',rgt_sense
      call write_log(1, bufout)
   endif

   ! initialize output-variables: element, node numbers and weights per output grid point

   do ii = 1, nout
      ii2iel(ii)      = 0
      ii2nod(1:4,ii)  = 0
      wii2nod(1:4,ii) = 0d0
   enddo

   ! loop over all input surface elements, determine which output grid points are covered

   do iel_y = 1, nel_y
      do iel_x = 1, nel_x
         iel = iel_x + (iel_y-1)*nel_x

         ! compute bounding box for input surface element iel

         xin_min = x_node(iel2node(1,iel))
         xin_max = xin_min
         yin_min = y_node(iel2node(1,iel))
         yin_max = yin_min
         zin_min = z_node(iel2node(1,iel))
         zin_max = zin_min
         do i = 2, 4
            xin_min = min(xin_min, x_node(iel2node(i,iel)))
            xin_max = max(xin_max, x_node(iel2node(i,iel)))
            yin_min = min(yin_min, y_node(iel2node(i,iel)))
            yin_max = max(yin_max, y_node(iel2node(i,iel)))
            zin_min = min(zin_min, z_node(iel2node(i,iel)))
            zin_max = max(zin_max, z_node(iel2node(i,iel)))
         enddo

         ! detect invalid data for input surface: NaN positions

         l_nan = (isnan(xin_min) .or. isnan(xin_max) .or. isnan(yin_min) .or. isnan(yin_max) .or.       &
                  isnan(zin_min) .or. isnan(zin_max))

         ! ignore surface elements that don't overlap with the extent of the output grid

         l_far = (xin_min.gt.xout_hig .or. yin_min.gt.yout_hig .or. zin_min.gt. z_thrs .or.             &
                  xin_max.lt.xout_low .or. yin_max.lt.yout_low .or. zin_max.lt.-z_thrs)

         if (l_nan) then
            if (ldebug.ge.-1 .and. num_error_msg.le.5) then
               write(bufout,'(a,i7,2(a,i4),a)') '  surface element',iel,' (',iel_x,',',iel_y,           &
                        '): NaN coordinates'
               call write_log(1, bufout)
               if (ldebug.ge.-2) then
                  i_ll = iel2node(1,iel)
                  i_lr = iel2node(2,iel)
                  i_ur = iel2node(3,iel)
                  i_ul = iel2node(4,iel)
                  write(bufout, 980) 'lower-left: ',i_ll, x_node(i_ll), y_node(i_ll)
                  call write_log(1,bufout)
                  write(bufout, 980) 'lower-right:',i_lr, x_node(i_lr), y_node(i_lr)
                  call write_log(1,bufout)
                  write(bufout, 980) 'upper-right:',i_ur, x_node(i_ur), y_node(i_ur)
                  call write_log(1,bufout)
                  write(bufout, 980) 'upper-left: ',i_ul, x_node(i_ul), y_node(i_ul)
                  call write_log(1,bufout)
               endif
               num_error_msg = num_error_msg + 1
            endif
         elseif (l_far) then
            if (ldebug.ge.15 .or. iel.eq.iel_debug) then
               write(bufout,'(a,i7,a)') '  surface element',iel,': no overlap with output grid'
               call write_log(1, bufout)
            endif
         else
            if (ldebug.ge.14 .or. iel.eq.iel_debug) then
               write(bufout,'(a,i7,2(a,i4),a)') '  Processing surface element',iel,' (',iel_x,',',iel_y,')'
               call write_log(1, bufout)
            endif

            ! node numbers of four corners of surface element iel (counter-clockwise, i_ll, i_lr, i_ur, i_ul)

            i_ll = iel2node(1,iel)
            i_lr = iel2node(2,iel)
            i_ur = iel2node(3,iel)
            i_ul = iel2node(4,iel)

            if (ldebug.ge.13) then
               write(bufout, 980) 'lower-left: ',i_ll, x_node(i_ll), y_node(i_ll)
               call write_log(1,bufout)
               write(bufout, 980) 'lower-right:',i_lr, x_node(i_lr), y_node(i_lr)
               call write_log(1,bufout)
               write(bufout, 980) 'upper-right:',i_ur, x_node(i_ur), y_node(i_ur)
               call write_log(1,bufout)
               write(bufout, 980) 'upper-left: ',i_ul, x_node(i_ul), y_node(i_ul)
               call write_log(1,bufout)
 980           format(2x,a,' node',i7,' at (',f9.3,',',f9.3,')')
            endif

            ! vectors for the four sides of surface element iel, going round counter-clockwise
            !   lower-left -> lower-right -> upper-right -> upper-left -> lower-left

            tlow_x = x_node(i_lr) - x_node(i_ll)
            tlow_y = y_node(i_lr) - y_node(i_ll)
            trgt_x = x_node(i_ur) - x_node(i_lr)
            trgt_y = y_node(i_ur) - y_node(i_lr)
            tupp_x = x_node(i_ul) - x_node(i_ur)
            tupp_y = y_node(i_ul) - y_node(i_ur)
            tlft_x = x_node(i_ll) - x_node(i_ul)
            tlft_y = y_node(i_ll) - y_node(i_ul)

            ! inward normal vectors for the four sides of surface element iel

            nlow_x = -tlow_y
            nlow_y =  tlow_x
            nrgt_x = -trgt_y
            nrgt_y =  trgt_x
            nupp_x = -tupp_y
            nupp_y =  tupp_x
            nlft_x = -tlft_y
            nlft_y =  tlft_x

            ! inner product value at the four sides of surface element iel

            ofs_low = rgt_sense * (x_node(i_ll)*nlow_x + y_node(i_ll)*nlow_y)
            ofs_rgt = rgt_sense * (x_node(i_lr)*nrgt_x + y_node(i_lr)*nrgt_y)
            ofs_upp = rgt_sense * (x_node(i_ur)*nupp_x + y_node(i_ur)*nupp_y)
            ofs_lft = rgt_sense * (x_node(i_ul)*nlft_x + y_node(i_ul)*nlft_y)

            ! find output grid points that are within the bounding box of input element iel

            if (is_uniform_x) then
               ix0 = max( 1, nint((xin_min-0.1*dx - xout_low)/dx)+1 )
               ix1 = min(mx, nint((xin_max+0.1*dx - xout_low)/dx)   )
            else
               ix0 = -1
               ix1 = -1
               call locate_segment( mx, xo_list, xin_min-0.1*dx, ix0 )
               call locate_segment( mx, xo_list, xin_max+0.1*dx, ix1 )
               ix0 = min(mx, max( 1, ix0+1))
               ix1 = min(mx, max( 1, ix1))
            endif

            if (is_uniform_y) then
               iy0 = max( 1, nint((yin_min-0.1*dy - yout_low)/dy)+1 )
               iy1 = min(my, nint((yin_max+0.1*dy - yout_low)/dy)   )
            else
               iy0 = -1
               iy1 = -1
               call locate_segment( my, yo_list, yin_min-0.1*dy, iy0 )
               call locate_segment( my, yo_list, yin_max+0.1*dy, iy1 )
               iy0 = min(my, max( 1, iy0+1))
               iy1 = min(my, max( 1, iy1))
            endif

            if (ldebug.ge.12 .or. iel.eq.iel_debug) then
               write(bufout,'(4(a,i4),a)') '  possible overlap: ix_out = [',ix0,',',ix1,'], iy = [',    &
                        iy0,',',iy1,']'
               call write_log(1, bufout)
            endif

            ! for all output grid points within bounding box,
            !     compute inner products with four sides of element iel

            do iy = iy0, iy1
               do ix = ix0, ix1
                  ii = ix + (iy-1) * mx

                  if (ldebug.ge.20 .or. (ldebug.ge.5 .and. iel.eq.iel_debug)) then
                     write(bufout,'(2(a,i4),2(a,f9.3),a)') '  checking (ix,iy) = (',ix,',',iy,          &
                                '), (x,y) = (',x_out(ii),',',y_out(ii),')'
                     call write_log(1, bufout)
                  endif

                  ! check if a value was computed before for another surface element

                  in_cell = (ii2iel(ii).eq.0)
                  if (ldebug.ge.9 .and. iel.eq.iel_debug) then
                     write(bufout,*) 'ii=',ii,', ii2iel=',ii2iel(ii),', in_cell=',in_cell
                     call write_log(1,bufout)
                  endif

                  ! check if the output grid point is above the lower edge

                  if (in_cell .and. ldebug.ge.9 .and. iel.eq.iel_debug) then
                     write(bufout,1970) rgt_sense, x_out(ii), nlow_x, y_out(ii), nlow_y,                &
                                        rgt_sense * (x_out(ii)*nlow_x + y_out(ii)*nlow_y),              &
                                        ofs_low, 'low'
                     call write_log(1, bufout)
  1970               format(' val=',f4.1,' *(',f8.3,'*',f8.3,' +', f8.3,'*',f8.3,') =', f10.6,          &
                            '   ?>=?', f10.6,' = ofs_',a)
                  endif
                  if (in_cell) in_cell = rgt_sense * (x_out(ii)*nlow_x + y_out(ii)*nlow_y) .ge. ofs_low

                  ! check if the output grid point is to the left of the right edge

                  if (in_cell .and. ldebug.ge.9 .and. iel.eq.iel_debug) then
                     write(bufout,1970) rgt_sense, x_out(ii), nrgt_x, y_out(ii), nrgt_y,                &
                                        rgt_sense * (x_out(ii)*nrgt_x + y_out(ii)*nrgt_y),              &
                                        ofs_rgt, 'rgt'
                     call write_log(1, bufout)
                  endif
                  if (in_cell) in_cell = rgt_sense * (x_out(ii)*nrgt_x + y_out(ii)*nrgt_y) .ge. ofs_rgt

                  ! check if the output grid point is below the upper edge

                  if (in_cell .and. ldebug.ge.9 .and. iel.eq.iel_debug) then
                     write(bufout,1970) rgt_sense, x_out(ii), nupp_x, y_out(ii), nupp_y,                &
                                        rgt_sense * (x_out(ii)*nupp_x + y_out(ii)*nupp_y),              &
                                        ofs_upp, 'upp'
                     call write_log(1, bufout)
                  endif
                  if (in_cell) in_cell = rgt_sense * (x_out(ii)*nupp_x + y_out(ii)*nupp_y) .ge. ofs_upp

                  ! check if the output grid point is to the right of the left edge

                  if (in_cell .and. ldebug.ge.9 .and. iel.eq.iel_debug) then
                     write(bufout,1970) rgt_sense, x_out(ii), nlft_x, y_out(ii), nlft_y,                &
                                        rgt_sense * (x_out(ii)*nlft_x + y_out(ii)*nlft_y),              &
                                        ofs_lft, 'lft'
                     call write_log(1, bufout)
                  endif
                  if (in_cell) in_cell = rgt_sense * (x_out(ii)*nlft_x + y_out(ii)*nlft_y) .ge. ofs_lft

                  ! remember the input surface element number for this output grid point

                  if (in_cell) ii2iel(ii) = iel

                  ! compute weights for bi-linear interpolation (is it truely bilinear or just a good
                  !    approximation?)

                  if (in_cell) then

                     ! store surface node numbers for output grid point ii (counter-clockwise, ll, lr, ur, ul)

                     ii2nod(1,ii) = i_ll
                     ii2nod(2,ii) = i_lr
                     ii2nod(3,ii) = i_ur
                     ii2nod(4,ii) = i_ul

                     ! store corresponding interpolation weights for output grid point ii

                     my_ldebug = 0
                     if (ldebug.ge.9 .and. ii.eq.ii_debug) my_ldebug = 9
                     call compute_logical_uv(nnode, x_node, y_node, i_ll, i_lr, i_ur, i_ul,             &
                                             x_out(ii), y_out(ii), ii, fac_long, fac_lat, my_ldebug,    &
                                             sub_ierror)
                     if (my_ierror.eq.0) my_ierror = sub_ierror

                     ! relative positions needed for bicubic interpolation

                     fac_uv(1,ii)  = fac_long
                     fac_uv(2,ii)  = fac_lat

                     ! weighting factors for bilinear interpolation

                     wii2nod(1,ii) = (1d0-fac_long) * (1d0-fac_lat) ! * f_node(i_ll)
                     wii2nod(2,ii) =      fac_long  * (1d0-fac_lat) ! * f_node(i_lr)
                     wii2nod(3,ii) =      fac_long  *      fac_lat  ! * f_node(i_ur)
                     wii2nod(4,ii) = (1d0-fac_long) *      fac_lat  ! * f_node(i_ul)

                     if ((ldebug.ge.9 .and. ii.eq.ii_debug) .or. ldebug.ge.18 .and. ii.ge.-10) then
                        write(bufout, 990) ix, iy, x_out(ii), y_out(ii), fac_long, fac_lat
                        call write_log(1, bufout)
 990                    format('  interpolation for grid point (',i3,',',i3,') at (',f7.3,',',f9.3,'),', &
                               ' position in cell (fx,fy) = (',f6.3,',',f6.3,')')
                     endif

                  elseif (ldebug.ge.20 .or. (ldebug.ge.9 .and. iel.eq.iel_debug)) then

                     write(bufout, '(2(a,i4),a,i7)') '  output grid point (',ix,',',iy,                &
                               ') does not lie in element',iel
                     call write_log(1, bufout)

                  endif
               enddo
            enddo

         endif ! .not.l_far
      enddo ! iel_x
   enddo ! iel_y

   deallocate(iel2node)
end subroutine interp_wgt_surf2unif

!------------------------------------------------------------------------------------------------------------

module subroutine compute_logical_uv(nnode, x_node, y_node, i_ll, i_lr, i_ur, i_ul, x_out, y_out, ii,          &
                              u_out, v_out, loc_debug, ierror)
!--function: compute relative [u,v] position within quadrilateral element for given physical [x,y]
!            coordinates for use in bilinear interpolation
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: nnode                         ! number of points in the input grid
   real(kind=8), intent(in)  :: x_node(nnode), y_node(nnode)  ! input grid [x,y] coordinates
   integer,      intent(in)  :: i_ll, i_lr, i_ur, i_ul        ! corners of quadrilateral A-B-C-D
   real(kind=8), intent(in)  :: x_out, y_out                  ! physical position of output point P
   integer,      intent(in)  :: ii                            ! point number - for print output only
   real(kind=8), intent(out) :: u_out, v_out                  ! logical position of output point P
   integer,      intent(in)  :: loc_debug                     ! level of debug output for current point
   integer,      intent(out) :: ierror
!--local variables:
   integer,      parameter :: maxit = 20    ! max.iterations for Newton algorithm
   real(kind=8), parameter :: eps = 1d-8    ! absolute tolerance on (u,v)
   integer,      save :: num_errors = 0
   real(kind=8), save :: max_dist   = 0d0
   integer       :: k
   real(kind=8)  :: pk(2), fk(2), jac(2,2), dpk(2), det, upd
   real(kind=8)  :: tlow_x, tlow_y, trgt_x, trgt_y, tupp_x, tupp_y, tlft_x, tlft_y
   real(kind=8)  :: nlong_x, nlong_y, sq_long, ofs_long, fac_long,                      &
                    nlat_x, nlat_y, sq_lat, ofs_lat, fac_lat, dist

   ierror = 0
   associate( xa => x_node(i_ll), xb => x_node(i_lr), xc => x_node(i_ur), xd => x_node(i_ul),        &
              ya => y_node(i_ll), yb => y_node(i_lr), yc => y_node(i_ur), yd => y_node(i_ul),        &
              xp => x_out,        yp => y_out,        uk => pk(1),        vk => pk(2))

   ! Original ad-hoc linear calculation

   tlow_x = xb - xa
   tlow_y = yb - ya
   trgt_x = xc - xb
   trgt_y = yc - yb
   tupp_x = xd - xc
   tupp_y = yd - yc
   tlft_x = xa - xd
   tlft_y = ya - yd
   
   ! normals for computing local coordinates for bi-linear interpolation

   nlong_x = 0.5*(tlow_x - tupp_x)
   nlong_y = 0.5*(tlow_y - tupp_y)
   sq_long = nlong_x**2 + nlong_y**2
   ofs_long = (xa*nlong_x + ya*nlong_y) / sq_long

   nlat_x  = 0.5*(trgt_x - tlft_x)
   nlat_y  = 0.5*(trgt_y - tlft_y)
   sq_lat  = nlat_x**2 + nlat_y**2
   ofs_lat  = (xa*nlat_x  + ya*nlat_y ) / sq_lat

   fac_long = (x_out*nlong_x + y_out*nlong_y) / sq_long - ofs_long
   fac_lat  = (x_out*nlat_x  + y_out*nlat_y ) / sq_lat  - ofs_lat

! using counter-clockwise ordering of corners:
!                 D         C
!    ^             o <---- o
!    |             |    p  ^
!  v |             v       |
!    ----->        o ----> o
!      u          A         B
!
! solving the equations
!       xp  = xa + (xb - xa) * u + (xd - xa) * v + ((xc - xb) - (xd - xa)) * u * v
!       yp  = ya + (yb - ya) * u + (yd - ya) * v + ((yc - yb) - (yd - ya)) * u * v
! 
!   fx(u,v) = xa + (xb - xa) * u + (xd - xa) * v + ((xc - xb) - (xd - xa)) * u * v - xp = 0
!   fy(u,v) = ya + (yb - ya) * u + (yd - ya) * v + ((yc - yb) - (yd - ya)) * u * v - yp = 0
!
! using Newton,  fk = f( pk ) = (uk,vk), Jk = [ df / dp ] at pk,  
!     fk1 \approx fk + Jk (pk1 - pk) = 0  
!       -->  dpk = - Jk \ fk 
!       -->  pk1 = pk - Jk \ fk

   ! New implementation using bilinear model

   k     = 0
   pk    = (/ 0.5d0, 0.5d0 /)
   upd   = 1d0

   do while(upd.gt.eps .and. k.le.maxit)

      k = k + 1

      ! compute residual fk

      fk(1) = xa + (xb - xa) * uk + (xd - xa) * vk + ((xc - xb) - (xd - xa)) * uk * vk - xp
      fk(2) = ya + (yb - ya) * uk + (yd - ya) * vk + ((yc - yb) - (yd - ya)) * uk * vk - yp

      ! compute Jacobian Jk

      jac(1,1) =   (xb - xa)                       + ((xc - xb) - (xd - xa))      * vk
      jac(1,2) =                    (xd - xa)      + ((xc - xb) - (xd - xa)) * uk     
      jac(2,1) =   (yb - ya)                       + ((yc - yb) - (yd - ya))      * vk
      jac(2,2) =                    (yd - ya)      + ((yc - yb) - (yd - ya)) * uk     

      ! solve update dpk

      det = jac(1,1) * jac(2,2) - jac(2,1) * jac(1,2)
      if (abs(det).lt.1d-8 .and. num_errors.lt.10) then
         num_errors = num_errors + 1
         write(bufout, '(a,/, 2(a,f14.8),a,/, 3(a,f14.8))')                                             &
                ' Internal error (u,v): matrix is near singular',                                       &
                ' Jac= [', jac(1,1),',', jac(1,2),']',                                                  &
                '      [', jac(2,1),',', jac(2,2),'], det=', det
         call write_log(3, bufout)
         write(bufout,'(5(a,i5,2(a,f12.6),a,:/))')                                                      &
               '   A(', i_ll,') = (', xa,',', ya,')', '   B(', i_lr,') = (', xb,',', yb,')',            &
               '   C(', i_ur,') = (', xc,',', yc,')', '   D(', i_ul,') = (', xd,',', yd,')',            &
               '   P(', ii  ,') = (', xp,',', yp,')'
         call write_log(5, bufout)
      endif

      !      1  [  j22  -j12 ]    [ f1 ]
      !   - --- |            |  * |    |
      !     det [ -j21   j11 ]    [ f2 ]

      dpk(1) = - ( jac(2,2) * fk(1) - jac(1,2) * fk(2)) / det
      dpk(2) = - (-jac(2,1) * fk(1) + jac(1,1) * fk(2)) / det

      ! compute max-norm of update

      upd = max(abs(dpk(1)), abs(dpk(2)))

      ! update vector pk = [uk; vk]

      pk(1:2) = pk(1:2) + dpk(1:2)

      if (loc_debug.ge.9 .or. (abs(det).lt.1d-8 .and. num_errors.le.10)) then
         write(bufout, '(a,i3,3(a,f12.8))') ' k=', k,': u=', pk(1),', v=', pk(2),' dv=', upd
         call write_log(1, bufout)
      endif

   enddo

   if (num_errors.lt.20 .and. upd.gt.eps) then
      num_errors = num_errors + 1
      write(bufout,*) 'ERROR in logical [u,v], no convergence for ii=',ii
      call write_log(1, bufout)
   endif

   if (.true.) then
      u_out = pk(1)
      v_out = pk(2)
   else
      u_out = fac_long
      v_out = fac_lat
   endif

   dist = sqrt( (u_out-fac_long)**2 + (v_out-fac_lat)**2 )
   if (.false. .and. dist.gt.max_dist+0.01) then
      max_dist = dist

      write(bufout,'(a,i5,a,f6.3,/, 4(a,f7.3),a)') ' Wrong results for [u,v] for ii=',ii,', dist=',dist,   &
                                          ' Orig: (', fac_long,',', fac_lat,'), new: (', u_out,',', v_out,')'
      call write_log(2, bufout)
      write(bufout,'(5(a,i5,2(a,f12.6),a,:/))')                                                         &
               '   A(', i_ll,') = (', xa,',', ya,')', '   B(', i_lr,') = (', xb,',', yb,')',            &
               '   C(', i_ur,') = (', xc,',', yc,')', '   D(', i_ul,') = (', xd,',', yd,')',            &
               '   P(', ii  ,') = (', xp,',', yp,')'
      call write_log(5, bufout)
   endif

   end associate
end subroutine compute_logical_uv

!------------------------------------------------------------------------------------------------------------

module subroutine interp_aply_surf2unif_1d(nnode, arr_node, nout, arr_out, ii2iel, ii2nod, wii2nod, ierror, defval)
!--function: interpolate scalar (1d) data given at the nodes of an input surface to the points of the
!            output grid
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: nnode             ! number of nodes in the input surface
   real(kind=8), intent(in)  :: arr_node(nnode)   ! function values of surface nodes 
                                                  ! (nval=3: associated with (x,y,z) dirs)
   integer,      intent(in)  :: nout              ! number of points in the output grid
   integer,      intent(in)  :: ii2iel(nout)      ! surface element number for output grid points
   integer,      intent(in)  :: ii2nod(4,nout)    ! surface node numbers for output grid points
   real(kind=8), intent(in)  :: wii2nod(4,nout)   ! interpolation weights per surface node per output point
   real(kind=8), intent(out) :: arr_out(nout)     ! function values interpolated to the output grid
   integer,      intent(out) :: ierror
   real(kind=8), optional    :: defval            ! default value
!--local variables:
!  character(len=*), parameter  :: subnam = 'interp_aply_surf2unif_1d'
   integer          :: ii, j
   logical          :: set_defval

   ierror = 0
   set_defval = .false.
   if (present(defval)) set_defval = .true.

   ! compute values of function arr_node at locations of the output grid

   do ii = 1, nout

      ! if a surface element iel is found for output grid point ii:

      if (ii2iel(ii).ne.0) then

         ! sum up contributions of corners i_ll, i_lr, i_ur, i_ul (counter-clockwise)

         arr_out(ii) = 0d0

         do j = 1, 4

            if (ldebug.ge.9 .and. ii.eq.ii_debug) then
               write(bufout,*) 'ii_out=',ii,', cornr j=',j,': surf.nod=',ii2nod(j,ii)
               call write_log(1, bufout)
            endif

            arr_out(ii) =  arr_out(ii) + wii2nod(j,ii) * arr_node(ii2nod(j,ii))

            if (ldebug.ge.9 .and. ii.eq.ii_debug) then
               write(bufout,*) ' wgt=',wii2nod(j,ii),', arr_node=', arr_node(ii2nod(j,ii))
               call write_log(1, bufout)
            endif

         enddo

      else
         if (set_defval) arr_out(ii) = defval
         if (ldebug.ge.20 .and. ii.ge.-10) then
            write(bufout, '(a,i7,a)') '  no value assigned to output grid point iout=', ii,'.'
            call write_log(1, bufout)
         endif
      endif
   enddo

end subroutine interp_aply_surf2unif_1d

!------------------------------------------------------------------------------------------------------------

module subroutine interp_aply_surf2unif_1d_bicubic(nnode_x, nnode_y, nnode, f_node, nout, ii2nod, &
                                            fac_uv, f_out, defval)
!--function: perform bicubic interpolation based on relative positions
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: nnode_x, nnode_y, nnode ! number of points in the input grid
   real(kind=8), intent(in)  :: f_node(nnode)           ! function values of input grid points
   integer,      intent(in)  :: nout                    ! number of points in the output grid
   integer,      intent(in)  :: ii2nod(4,nout)          ! input node numbers for output grid points
   real(kind=8), intent(in)  :: fac_uv(2,nout)          ! relative pos. of output points in input elements
   real(kind=8), intent(out) :: f_out(nout)             ! interpolated function values
   real(kind=8), optional    :: defval                  ! default value
!--local variables
   integer,      parameter :: ll_ngb = 1, lr_ngb = 2, ul_ngb = 3, ur_ngb = 4 ! rows in deriv_ngb
   integer,      parameter :: i_ngb  = 1, j_ngb  = 2, ij_ngb = 3             ! columns in deriv_ngb
   real(kind=8), parameter :: bicubic_mat(4,4) = reshape( &
                                             (/ 1d0, 0d0, -3d0,  2d0, &
                                                0d0, 0d0,  3d0, -2d0, &
                                                0d0, 1d0, -2d0,  1d0, &
                                                0d0, 0d0, -1d0,  1d0 /), (/ 4, 4 /))
   real(kind=8), parameter :: bicubic_mat_t(4,4) = transpose(bicubic_mat)
   integer                 :: ii, i_ll, i_lr, i_ur, i_ul
   integer                 :: deriv_ngb(4,3)
   real(kind=8)            :: deriv_dist(4,2)
   real(kind=8)            :: deriv_dist_ij(4,2)
   real(kind=8)            :: fac_long, fac_lat
   real(kind=8)            :: deriv_low, deriv_hgh
   real(kind=8)            :: pow_u(4), pow_v(4)
   real(kind=8)            :: bicubic_f_arr(4,4), alpha(4,4)
   logical                 :: set_defval

   set_defval = .false.
   if (present(defval)) set_defval = .true.

   ! initialize interpolation result

   f_out(1:nout) = 0d0

   ! loop over output positions

   do ii = 1, nout

      ! skip this node if it is not in an input element
      if (ii2nod(1,ii).eq.0) then
         if (set_defval) f_out(ii) = defval
         if (ldebug.ge.3 .and. ii.eq.ii_debug) then
            write(bufout,'(a,i6,a,i4,a,g12.4)') ' point ii=',ii,': ii2nod=',ii2nod(1,ii),', f_out=',f_out(ii)
            call write_log(1, bufout)
         endif
         cycle
      endif

      ! get corners of element surrounding point ii, get relative position within element

      i_ll = ii2nod(1,ii)
      i_lr = ii2nod(2,ii)
      i_ur = ii2nod(3,ii)
      i_ul = ii2nod(4,ii)
      fac_long = fac_uv(1,ii)
      fac_lat  = fac_uv(2,ii)

      ! determine derivative neighbors
      ! TODO: These neighbors only have to be determined once per source element,
      !       code should be added that stores the neighbors for re-use

      ! determine neighbors below element

      if (mod(i_ll-1, nnode_x).eq.0) then
         ! nodes ll and ul have no neighbors with lower i index
         deriv_ngb(ll_ngb, i_ngb) = i_ll
         deriv_ngb(ul_ngb, i_ngb) = i_ul
         deriv_dist(ll_ngb, i_ngb) = 1d0
         deriv_dist(ul_ngb, i_ngb) = 1d0
      else
         deriv_ngb(ll_ngb, i_ngb) = i_ll - 1
         deriv_ngb(ul_ngb, i_ngb) = i_ul - 1
         deriv_dist(ll_ngb, i_ngb) = 2d0
         deriv_dist(ul_ngb, i_ngb) = 2d0
      endif

      ! determine neighbors above element

      if (mod(i_lr, nnode_x).eq.0) then
         ! nodes lr and ur have no neighbors with higher i index
         deriv_ngb(lr_ngb, i_ngb) = i_lr
         deriv_ngb(ur_ngb, i_ngb) = i_ur
         deriv_dist(lr_ngb, i_ngb) = 1d0
         deriv_dist(ur_ngb, i_ngb) = 1d0
      else
         deriv_ngb(lr_ngb, i_ngb) = i_lr + 1
         deriv_ngb(ur_ngb, i_ngb) = i_ur + 1
         deriv_dist(lr_ngb, i_ngb) = 2d0
         deriv_dist(ur_ngb, i_ngb) = 2d0
      endif

      ! determine neighbors left of element

      if (i_ll.le.nnode_x) then
         ! nodes ll and lr have no neighbors with lower j index
         deriv_ngb(ll_ngb, j_ngb) = i_ll
         deriv_ngb(lr_ngb, j_ngb) = i_lr
         deriv_dist(ll_ngb, j_ngb) = 1d0
         deriv_dist(lr_ngb, j_ngb) = 1d0
      else
         deriv_ngb(ll_ngb, j_ngb) = i_ll - nnode_x
         deriv_ngb(lr_ngb, j_ngb) = i_lr - nnode_x
         deriv_dist(ll_ngb, j_ngb) = 2d0
         deriv_dist(lr_ngb, j_ngb) = 2d0
      endif

      ! determine neighbors right of element

      if (i_ul.gt.nnode_x*(nnode_y-1)) then
         ! nodes ul and ur have no neighbors with higher j index
         deriv_ngb(ul_ngb, j_ngb) = i_ul
         deriv_ngb(ur_ngb, j_ngb) = i_ur
         deriv_dist(ul_ngb, j_ngb) = 1d0
         deriv_dist(ur_ngb, j_ngb) = 1d0
      else
         deriv_ngb(ul_ngb, j_ngb) = i_ul + nnode_x
         deriv_ngb(ur_ngb, j_ngb) = i_ur + nnode_x
         deriv_dist(ul_ngb, j_ngb) = 2d0
         deriv_dist(ur_ngb, j_ngb) = 2d0
      endif

      ! determine diagonal neighbor below left of element

      if (mod(i_ll-1, nnode_x).eq.0) then
         ! node ll has no neighbors with lower i index
         if (i_ll.le.nnode_x) then
            ! node ll has no neighbors with lower j index
            deriv_ngb(ll_ngb, ij_ngb) = i_ll
            deriv_dist_ij(ll_ngb, i_ngb) = 1d0
            deriv_dist_ij(ll_ngb, j_ngb) = 1d0
         else
            ! node ll has neighbors with lower j index
            deriv_ngb(ll_ngb, ij_ngb) = i_ll - nnode_x
            deriv_dist_ij(ll_ngb, i_ngb) = 1d0
            deriv_dist_ij(ll_ngb, j_ngb) = 2d0
         endif
      else
         ! node ll has neighbors with lower i index
         if (i_ll.le.nnode_x) then
            ! node ll has no neighbors with lower j index
            deriv_ngb(ll_ngb, ij_ngb) = i_ll - 1
            deriv_dist_ij(ll_ngb, i_ngb) = 2d0
            deriv_dist_ij(ll_ngb, j_ngb) = 1d0
         else
            ! node ll has neighbors with lower j index
            deriv_ngb(ll_ngb, ij_ngb) = i_ll - 1 - nnode_x
            deriv_dist_ij(ll_ngb, i_ngb) = 2d0
            deriv_dist_ij(ll_ngb, j_ngb) = 2d0
         endif
      endif

      ! determine diagonal neighbor below right of element

      if (mod(i_lr, nnode_x).eq.0) then
         ! node lr has no neighbors with higher i index
         if (i_ll.le.nnode_x) then
            ! node lr has no neighbors with lower j index
            deriv_ngb(lr_ngb, ij_ngb) = i_lr
            deriv_dist_ij(lr_ngb, i_ngb) = 1d0
            deriv_dist_ij(lr_ngb, j_ngb) = 1d0
         else
            ! node lr has neighbors with lower j index
            deriv_ngb(lr_ngb, ij_ngb) = i_lr - nnode_x
            deriv_dist_ij(lr_ngb, i_ngb) = 1d0
            deriv_dist_ij(lr_ngb, j_ngb) = 2d0
         endif
      else
         ! node lr has neighbors with higher i index
         if (i_ll.le.nnode_x) then
            ! node lr has no neighbors with lower j index
            deriv_ngb(lr_ngb, ij_ngb) = i_lr + 1
            deriv_dist_ij(lr_ngb, i_ngb) = 2d0
            deriv_dist_ij(lr_ngb, j_ngb) = 1d0
         else
            ! node lr has neighbors with lower j index
            deriv_ngb(lr_ngb, ij_ngb) = i_lr + 1 - nnode_x
            deriv_dist_ij(lr_ngb, i_ngb) = 2d0
            deriv_dist_ij(lr_ngb, j_ngb) = 2d0
         endif
      endif

      ! determine diagonal neighbor above left of element

      if (mod(i_ll-1, nnode_x).eq.0) then
         ! node ul has no neighbors with lower i index
         if (i_ul.gt.nnode_x*(nnode_y-1)) then
            ! node ul has no neighbors with higher j index
            deriv_ngb(ul_ngb, ij_ngb) = i_ul
            deriv_dist_ij(ul_ngb, i_ngb) = 1d0
            deriv_dist_ij(ul_ngb, j_ngb) = 1d0
         else
            ! node ul has neighbors with higher j index
            deriv_ngb(ul_ngb, ij_ngb) = i_ul + nnode_x
            deriv_dist_ij(ul_ngb, i_ngb) = 1d0
            deriv_dist_ij(ul_ngb, j_ngb) = 2d0
         endif
      else
         ! node ul has neighbors with lower i index
         if (i_ul.gt.nnode_x*(nnode_y-1)) then
            ! node ul has no neighbors with higher j index
            deriv_ngb(ul_ngb, ij_ngb) = i_ul - 1
            deriv_dist_ij(ul_ngb, i_ngb) = 2d0
            deriv_dist_ij(ul_ngb, j_ngb) = 1d0
         else
            ! node ul has neighbors with higher j index
            deriv_ngb(ul_ngb, ij_ngb) = i_ul - 1 + nnode_x
            deriv_dist_ij(ul_ngb, i_ngb) = 2d0
            deriv_dist_ij(ul_ngb, j_ngb) = 2d0
         endif
      endif

      ! determine diagonal neighbor above right of element

      if (mod(i_lr, nnode_x).eq.0) then
         ! node ur has no neighbors with higher i index
         if (i_ur.gt.nnode_x*(nnode_y-1)) then
            ! node ur has no neighbors with higher j index
            deriv_ngb(ur_ngb, ij_ngb) = i_ur
            deriv_dist_ij(ur_ngb, i_ngb) = 1d0
            deriv_dist_ij(ur_ngb, j_ngb) = 1d0
         else
            ! node ur has neighbors with higher j index
            deriv_ngb(ur_ngb, ij_ngb) = i_ur + nnode_x
            deriv_dist_ij(ur_ngb, i_ngb) = 1d0
            deriv_dist_ij(ur_ngb, j_ngb) = 2d0
         endif
      else
         ! node ur has neighbors with higher i index
         if (i_ur.gt.nnode_x*(nnode_y-1)) then
            ! node ur has no neighbors with higher j index
            deriv_ngb(ur_ngb, ij_ngb) = i_ur + 1
            deriv_dist_ij(ur_ngb, i_ngb) = 2d0
            deriv_dist_ij(ur_ngb, j_ngb) = 1d0
         else
            ! node ur has neighbors with higher j index
            deriv_ngb(ur_ngb, ij_ngb) = i_ur + 1 + nnode_x
            deriv_dist_ij(ur_ngb, i_ngb) = 2d0
            deriv_dist_ij(ur_ngb, j_ngb) = 2d0
         endif
      endif

      ! get function values at four nodes surrounding output position

      bicubic_f_arr(1,1) = f_node(i_ll)
      bicubic_f_arr(2,1) = f_node(i_lr)
      bicubic_f_arr(1,2) = f_node(i_ul)
      bicubic_f_arr(2,2) = f_node(i_ur)

      ! get function derivatives in i direction

      ! f_node derivative in i direction at ll node

      bicubic_f_arr(3,1) = (f_node(i_lr) - f_node(deriv_ngb(ll_ngb, i_ngb))) / deriv_dist(ll_ngb, i_ngb)

      ! f_node derivative in i direction at lr node

      bicubic_f_arr(4,1) = (f_node(deriv_ngb(lr_ngb, i_ngb)) - f_node(i_ll)) / deriv_dist(lr_ngb, i_ngb)

      ! f_node derivative in i direction at ul node

      bicubic_f_arr(3,2) = (f_node(i_ur) - f_node(deriv_ngb(ul_ngb, i_ngb))) / deriv_dist(ul_ngb, i_ngb)

      ! f_node derivative in i direction at ur node

      bicubic_f_arr(4,2) = (f_node(deriv_ngb(ur_ngb, i_ngb)) - f_node(i_ul)) / deriv_dist(ur_ngb, i_ngb)

      ! get function derivatives in j direction

      ! f_node derivative in j direction at ll node

      bicubic_f_arr(1,3) = (f_node(i_ul) - f_node(deriv_ngb(ll_ngb, j_ngb))) / deriv_dist(ll_ngb, j_ngb)

      ! f_node derivative in j direction at lr node

      bicubic_f_arr(2,3) = (f_node(i_ur) - f_node(deriv_ngb(lr_ngb, j_ngb))) / deriv_dist(lr_ngb, j_ngb)

      ! f_node derivative in j direction at ul node

      bicubic_f_arr(1,4) = (f_node(deriv_ngb(ul_ngb, j_ngb)) - f_node(i_ll)) / deriv_dist(ul_ngb, j_ngb)

      ! f_node derivative in j direction at ur node

      bicubic_f_arr(2,4) = (f_node(deriv_ngb(ur_ngb, j_ngb)) - f_node(i_lr)) / deriv_dist(ur_ngb, j_ngb)

      ! get function cross derivatives in ij

      ! f_node cross derivative in ij at ll node

      deriv_low = (f_node(deriv_ngb(lr_ngb, j_ngb)) - f_node(deriv_ngb(ll_ngb, ij_ngb)))               &
                                                                             / deriv_dist_ij(ll_ngb, i_ngb)
      deriv_hgh = bicubic_f_arr(3,2)

      bicubic_f_arr(3,3) = (deriv_hgh - deriv_low) / deriv_dist_ij(ll_ngb, j_ngb)

      ! f_node cross derivative in ij at lr node

      deriv_low = (f_node(deriv_ngb(lr_ngb, ij_ngb)) - f_node(deriv_ngb(ll_ngb, j_ngb)))               &
                                                                             / deriv_dist_ij(lr_ngb, i_ngb)
      deriv_hgh = bicubic_f_arr(4,2)

      bicubic_f_arr(4,3) = (deriv_hgh - deriv_low) / deriv_dist_ij(lr_ngb, j_ngb)

      ! f_node cross derivative in ij at ul node

      deriv_low = bicubic_f_arr(3,1)
      deriv_hgh = (f_node(deriv_ngb(ur_ngb, j_ngb)) - f_node(deriv_ngb(ul_ngb, ij_ngb)))                &
                                                                             / deriv_dist_ij(ul_ngb, i_ngb)
      bicubic_f_arr(3,4) = (deriv_hgh - deriv_low) / deriv_dist_ij(ul_ngb, j_ngb)

      ! f_node cross derivative in ij at ur node

      deriv_low = bicubic_f_arr(4,1)
      deriv_hgh = (f_node(deriv_ngb(ur_ngb, ij_ngb)) - f_node(deriv_ngb(ul_ngb, j_ngb)))                &
                                                                             / deriv_dist_ij(ur_ngb, i_ngb)
      bicubic_f_arr(4,4) = (deriv_hgh - deriv_low) / deriv_dist_ij(ur_ngb, j_ngb)

      ! compute coefficients as explained on wikipedia

      alpha = matmul(matmul(bicubic_mat, bicubic_f_arr), bicubic_mat_t)

      ! compute powers of u, v needed for evaluation

      pow_u = (/ 1d0, fac_long, fac_long**2, fac_long**3 /)
      pow_v = (/ 1d0, fac_lat,  fac_lat**2,  fac_lat**3 /)
      f_out(ii) = dot_product(matmul(pow_u, alpha), pow_v)
   end do

end subroutine interp_aply_surf2unif_1d_bicubic

!------------------------------------------------------------------------------------------------------------

module subroutine interp_aply_surf2unif_3d(nnode, arr_node, nout, arr_out, ii2iel, ii2nod, wii2nod, ierror, defval)
!--function: interpolate 3d data given at the nodes of an input surface to the points of the output grid
   implicit none
!--subroutine arguments:
   integer,      intent(in)  :: nnode             ! number of nodes in the input surface
   real(kind=8), intent(in)  :: arr_node(nnode,3) ! function values of surface nodes 
                                                  ! (nval=3: associated with (x,y,z) dirs)
   integer,      intent(in)  :: nout              ! number of points in the output grid
   integer,      intent(in)  :: ii2iel(nout)      ! surface element number for output grid points
   integer,      intent(in)  :: ii2nod(4,nout)    ! surface node numbers for output grid points
   real(kind=8), intent(in)  :: wii2nod(4,nout)   ! interpolation weights per surface node per output point
   real(kind=8)              :: arr_out(nout,3)   ! function values interpolated to the output grid
   integer,      intent(out) :: ierror
   real(kind=8), optional    :: defval            ! default value
!--local variables:
!  character(len=*), parameter  :: subnam = 'interp_aply_surf2unif_3d'
   integer          :: ii, j, k
   logical          :: set_defval

   ierror = 0
   set_defval = .false.
   if (present(defval)) set_defval = .true.

   ! compute values of function arr_node at locations of the output grid

   do ii = 1, nout

      ! if a surface element iel is found for output grid point ii:

      if (ii2iel(ii).ne.0) then

         ! sum up contributions of corners i_ll, i_lr, i_ur, i_ul (counter-clockwise)

         arr_out(ii,1:3) = 0d0

         do j = 1, 4

            if (ldebug.ge.80 .and. ii.eq.26) then
               write(bufout,*) 'ii_out=',ii,', cornr j=',j,': surf.nod=',ii2nod(j,ii)
               call write_log(1, bufout)
            endif

            do k = 1, 3
               arr_out(ii,k) =  arr_out(ii,k) + wii2nod(j,ii) * arr_node(ii2nod(j,ii),k)

               if (ldebug.ge.80 .and. ii.eq.26) then
                  write(bufout,*) 'comp. k=',k,': wgt=',wii2nod(j,ii),', arr_node=', arr_node(ii2nod(j,ii),k)
                  call write_log(1, bufout)
               endif
            enddo
         enddo

      else
         if (set_defval) arr_out(ii,1:3) = defval
         if (ldebug.ge.20 .and. ii.ge.-10) then
            write(bufout, '(a,i7,a)') '  no value assigned to output grid point iout=',ii,'.'
            call write_log(1, bufout)
         endif
      endif
   enddo

end subroutine interp_aply_surf2unif_3d

!------------------------------------------------------------------------------------------------------------

end submodule m_interp_sub
