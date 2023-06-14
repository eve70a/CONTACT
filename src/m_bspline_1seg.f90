!------------------------------------------------------------------------------------------------------------
! m_bspline_1seg - get cubic function for one spline segment
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_bspline_1seg
   use m_globals
   use m_bspline_def
   implicit none
   private

   public  bspline_get_ppcoef_1seg_orig
   public  bspline_get_ppcoef_1seg_modf

contains

!------------------------------------------------------------------------------------------------------------

subroutine bspline_get_ppcoef_1seg_orig(nspl, tj, cj_y, jseg, c0, c1, c2, c3, idebug)
!--function: determine PP-spline coefficients for one segment jseg
   implicit none
!--subroutine arguments:
   integer,      intent(in)           :: nspl, jseg, idebug
   real(kind=8), intent(in)           :: cj_y(nspl), tj(nspl+4)
   real(kind=8), intent(out)          :: c0, c1, c2, c3
!--local variables:
   real(kind=8), parameter    :: small_dt = 1d-7
   integer              :: jj
   real(kind=8)         :: dtj, dtj_inv1(4), dtj_inv2(4), dtj_inv3(4), b1(4), b2(4), b3(4), b4(4),      &
                           v, d1_cy(4), d2_cy(4), d3_cy(4)

   ! average step sizes over 1/2/3 adjacent intervals

   do jj = 3, 4         ! dtj_k1: used for jseg:jseg
      dtj = tj(jseg-4+jj+1) - tj(jseg-4+jj)
      dtj_inv1(jj) = 0d0
      if (dtj.gt.small_dt) dtj_inv1(jj) = 1d0 / dtj
   enddo

   do jj = 2, 4         ! dtj_k2: used for jseg-1:jseg
      dtj = (tj(jseg-4+jj+2) - tj(jseg-4+jj)) / 2d0
      dtj_inv2(jj) = 0d0
      if (dtj.gt.small_dt) dtj_inv2(jj) = 1d0 / dtj
   enddo

   do jj = 1, 4         ! dtj_k3: used for jseg-2:jseg
      dtj = (tj(jseg-4+jj+3) - tj(jseg-4+jj)) / 3d0
      dtj_inv3(jj) = 0d0
      if (dtj.gt.small_dt) dtj_inv3(jj) = 1d0 / dtj
   enddo

   ! evaluate function and derivative values at position tj(jseg)

   v = tj(jseg) + 1d-9
   b1(3)   = 0d0
   b1(4)   = 1d0
   b2(2)   = 0d0
   do jj = 3 , 4
      b2(jj) =                       (v-tj(jseg-4+jj  )) * b1(jj  ) * dtj_inv1(jj  ) / 1d0
      if (jj.lt.4) b2(jj) = b2(jj) + (tj(jseg-4+jj+2)-v) * b1(jj+1) * dtj_inv1(jj+1) / 1d0
   enddo
   b3(1)   = 0d0
   do jj = 2 , 4
      b3(jj) =                       (v-tj(jseg-4+jj  )) * b2(jj  ) * dtj_inv2(jj  ) / 2d0
      if (jj.lt.4) b3(jj) = b3(jj) + (tj(jseg-4+jj+3)-v) * b2(jj+1) * dtj_inv2(jj+1) / 2d0
   enddo
   do jj = 1 , 4
      b4(jj) =                       (v-tj(jseg-4+jj  )) * b3(jj  ) * dtj_inv3(jj  ) / 3d0
      if (jj.lt.4) b4(jj) = b4(jj) + (tj(jseg-4+jj+4)-v) * b3(jj+1) * dtj_inv3(jj+1) / 3d0
   enddo

   if (idebug.ge.5) then
      write(bufout,'(a,36x, f12.4)') 'b1=', b1(4)
      call write_log(1, bufout)
      write(bufout,'(a,24x,2f12.4)') 'b2=', b2(3), b2(4)
      call write_log(1, bufout)
      write(bufout,'(a,12x,3f12.4)') 'b3=', b3(2), b3(3), b3(4)
      call write_log(1, bufout)
      write(bufout,'(a,    4f12.4)') 'b4=', b4(1), b4(2), b4(3), b4(4)
      call write_log(1, bufout)
   endif

   do jj = 2, 4
      d1_cy(jj) = (cj_y(jseg-4+jj) - cj_y(jseg-4+jj-1)) * dtj_inv3(jj)
   enddo
   do jj = 3, 4
      d2_cy(jj) = (d1_cy(jj) - d1_cy(jj-1)) * dtj_inv2(jj)
   enddo
   do jj = 4, 4
      d3_cy(jj) = (d2_cy(jj) - d2_cy(jj-1)) * dtj_inv1(jj)
   enddo

   if (idebug.ge.5) then
      write(bufout,'(a,12x,3f12.4)') 'd1_cy=', d1_cy(2), d1_cy(3), d1_cy(4)
      call write_log(1, bufout)
      write(bufout,'(a,24x,2f12.4)') 'd2_cy=', d2_cy(3), d2_cy(4)
      call write_log(1, bufout)
      write(bufout,'(a,36x, f12.4)') 'd3_cy=', d3_cy(4)
      call write_log(1, bufout)
   endif

   c0  = 0d0
   c1  = 0d0
   c2  = 0d0
   c3  = 0d0
   do jj = 1, 4
                   c0 = c0 + b4(jj) * cj_y(jseg-4+jj)
      if (jj.ge.2) c1 = c1 + b3(jj) * d1_cy(jj)
      if (jj.ge.3) c2 = c2 + b2(jj) * d2_cy(jj) / 2d0
      if (jj.ge.4) c3 = c3 + b1(jj) * d3_cy(jj) / 6d0
   enddo

   if (idebug.ge.5) then
      write(bufout,'(a,7f12.4)') 'c0-c3=',c0, c1, c2, c3
      call write_log(1, bufout)
   endif

   if (.false. .and. abs(c0-cj_y(jseg)).gt.1d0) then
      call write_log('orig:')
      write(bufout,'(a,24x,4f12.4)') 'dt1  =', dtj_inv1(3),dtj_inv1(4)
      call write_log(1, bufout)
      write(bufout,'(a,12x,4f12.4)') 'dt2  =', dtj_inv2(2),dtj_inv2(3),dtj_inv2(4)
      call write_log(1, bufout)
      write(bufout,'(a,    4f12.4)') 'dt3  =', dtj_inv3(1),dtj_inv3(2),dtj_inv3(3),dtj_inv3(4)
      call write_log(1, bufout)
      write(bufout,'(a,36x, f12.4)') 'b1   =', b1(4)
      call write_log(1, bufout)
      write(bufout,'(a,24x,2f12.4)') 'b2   =', b2(3),b2(4)
      call write_log(1, bufout)
      write(bufout,'(a,12x,3f12.4)') 'b3   =', b3(2),b3(3),b3(4)
      call write_log(1, bufout)
      write(bufout,'(a,    4f12.4)') 'b4   =', b4(1),b4(2),b4(3),b4(4)
      call write_log(1, bufout)
      write(bufout,'(a,12x,3f12.4)') 'd1_cy=', d1_cy(2),d1_cy(3),d1_cy(4)
      call write_log(1, bufout)
      write(bufout,'(a,24x,2f12.4)') 'd2_cy=', d2_cy(3),d2_cy(4)
      call write_log(1, bufout)
      write(bufout,'(a,36x, f12.4)') 'd3_cy=', d3_cy(4)
      call write_log(1, bufout)
   endif

end subroutine bspline_get_ppcoef_1seg_orig

!------------------------------------------------------------------------------------------------------------

subroutine bspline_get_ppcoef_1seg_modf(nspl, tj, cj_y, jseg, c0, c1, c2, c3, idebug)
!--function: determine PP-spline coefficients for one segment jseg
   implicit none
!--subroutine arguments:
   integer,      intent(in)           :: nspl, jseg, idebug
   real(kind=8), intent(in)           :: cj_y(nspl), tj(nspl+4)
   real(kind=8), intent(out)          :: c0, c1, c2, c3
!--local variables:
   real(kind=8), parameter    :: small_dt = 1d-7
   integer              :: jj
   real(kind=8)         :: dtj, dtj_inv1(4), dtj_inv2(4), dtj_inv3(4), b1(4), b2(4), b3(4), b4(4),      &
                           v, d1_cy(4), d2_cy(4), d3_cy(4)

   ! average step sizes over 1/2/3 adjacent intervals

   jj = 3
   dtj = tj(jseg-4+jj+1) - tj(jseg-4+jj)
   dtj_inv1(jj) = 0d0
   if (dtj.gt.small_dt) dtj_inv1(jj) = 1d0 / dtj
   jj = 4
   dtj = tj(jseg-4+jj+1) - tj(jseg-4+jj)
   dtj_inv1(jj) = 0d0
   if (dtj.gt.small_dt) dtj_inv1(jj) = 1d0 / dtj

   jj = 2
   dtj = (tj(jseg-4+jj+2) - tj(jseg-4+jj)) / 2d0
   dtj_inv2(jj) = 0d0
   if (dtj.gt.small_dt) dtj_inv2(jj) = 1d0 / dtj
   jj = 3
   dtj = (tj(jseg-4+jj+2) - tj(jseg-4+jj)) / 2d0
   dtj_inv2(jj) = 0d0
   if (dtj.gt.small_dt) dtj_inv2(jj) = 1d0 / dtj
   jj = 4
   dtj = (tj(jseg-4+jj+2) - tj(jseg-4+jj)) / 2d0
   dtj_inv2(jj) = 0d0
   if (dtj.gt.small_dt) dtj_inv2(jj) = 1d0 / dtj

   jj = 1
   dtj = (tj(jseg-4+jj+3) - tj(jseg-4+jj)) / 3d0
   dtj_inv3(jj) = 0d0
   if (dtj.gt.small_dt) dtj_inv3(jj) = 1d0 / dtj
   jj = 2
   dtj = (tj(jseg-4+jj+3) - tj(jseg-4+jj)) / 3d0
   dtj_inv3(jj) = 0d0
   if (dtj.gt.small_dt) dtj_inv3(jj) = 1d0 / dtj
   jj = 3
   dtj = (tj(jseg-4+jj+3) - tj(jseg-4+jj)) / 3d0
   dtj_inv3(jj) = 0d0
   if (dtj.gt.small_dt) dtj_inv3(jj) = 1d0 / dtj
   jj = 4
   dtj = (tj(jseg-4+jj+3) - tj(jseg-4+jj)) / 3d0
   dtj_inv3(jj) = 0d0
   if (dtj.gt.small_dt) dtj_inv3(jj) = 1d0 / dtj

   ! evaluate function and derivative values at position tj(jseg)

   v = tj(jseg) + 1d-9
   b1(3)  = 0d0
   b1(4)  = 1d0

   b2(2)  = 0d0
   jj = 3
   b2(jj) =          (v-tj(jseg-4+jj  )) * b1(jj  ) * dtj_inv1(jj  ) / 1d0
   b2(jj) = b2(jj) + (tj(jseg-4+jj+2)-v) * b1(jj+1) * dtj_inv1(jj+1) / 1d0
   jj = 4
   b2(jj) =          (v-tj(jseg-4+jj  )) * b1(jj  ) * dtj_inv1(jj  ) / 1d0

   b3(1)  = 0d0
   jj = 2
   b3(jj) =          (v-tj(jseg-4+jj  )) * b2(jj  ) * dtj_inv2(jj  ) / 2d0
   b3(jj) = b3(jj) + (tj(jseg-4+jj+3)-v) * b2(jj+1) * dtj_inv2(jj+1) / 2d0
   jj = 3
   b3(jj) =          (v-tj(jseg-4+jj  )) * b2(jj  ) * dtj_inv2(jj  ) / 2d0
   b3(jj) = b3(jj) + (tj(jseg-4+jj+3)-v) * b2(jj+1) * dtj_inv2(jj+1) / 2d0
   jj = 4
   b3(jj) =          (v-tj(jseg-4+jj  )) * b2(jj  ) * dtj_inv2(jj  ) / 2d0

   jj = 1
   b4(jj) =          (v-tj(jseg-4+jj  )) * b3(jj  ) * dtj_inv3(jj  ) / 3d0
   b4(jj) = b4(jj) + (tj(jseg-4+jj+4)-v) * b3(jj+1) * dtj_inv3(jj+1) / 3d0
   jj = 2
   b4(jj) =          (v-tj(jseg-4+jj  )) * b3(jj  ) * dtj_inv3(jj  ) / 3d0
   b4(jj) = b4(jj) + (tj(jseg-4+jj+4)-v) * b3(jj+1) * dtj_inv3(jj+1) / 3d0
   jj = 3
   b4(jj) =          (v-tj(jseg-4+jj  )) * b3(jj  ) * dtj_inv3(jj  ) / 3d0
   b4(jj) = b4(jj) + (tj(jseg-4+jj+4)-v) * b3(jj+1) * dtj_inv3(jj+1) / 3d0
   jj = 4
   b4(jj) =          (v-tj(jseg-4+jj  )) * b3(jj  ) * dtj_inv3(jj  ) / 3d0

   jj = 2
   d1_cy(jj) = (cj_y(jseg-4+jj) - cj_y(jseg-4+jj-1)) * dtj_inv3(jj)
   jj = 3
   d1_cy(jj) = (cj_y(jseg-4+jj) - cj_y(jseg-4+jj-1)) * dtj_inv3(jj)
   jj = 4
   d1_cy(jj) = (cj_y(jseg-4+jj) - cj_y(jseg-4+jj-1)) * dtj_inv3(jj)

   jj = 3
   d2_cy(jj) = (d1_cy(jj) - d1_cy(jj-1)) * dtj_inv2(jj)
   jj = 4
   d2_cy(jj) = (d1_cy(jj) - d1_cy(jj-1)) * dtj_inv2(jj)

   jj = 4
   d3_cy(jj) = (d2_cy(jj) - d2_cy(jj-1)) * dtj_inv1(jj)

   c0  = 0d0
   c1  = 0d0
   c2  = 0d0
   c3  = 0d0

   jj = 1
   c0 = c0 + b4(jj) * cj_y(jseg-4+jj)
   jj = 2
   c0 = c0 + b4(jj) * cj_y(jseg-4+jj)
   c1 = c1 + b3(jj) * d1_cy(jj)
   jj = 3
   c0 = c0 + b4(jj) * cj_y(jseg-4+jj)
   c1 = c1 + b3(jj) * d1_cy(jj)
   c2 = c2 + b2(jj) * d2_cy(jj) / 2d0
   jj = 4
   c0 = c0 + b4(jj) * cj_y(jseg-4+jj)
   c1 = c1 + b3(jj) * d1_cy(jj)
   c2 = c2 + b2(jj) * d2_cy(jj) / 2d0
   c3 = c3 + b1(jj) * d3_cy(jj) / 6d0

   if (.false. .and. idebug.ge.0 .and. abs(c0-cj_y(jseg)).gt.1d0) then
      call write_log('modf:')
      write(bufout,'(a,24x,4f12.4)') 'dt1  =',dtj_inv1(3),dtj_inv1(4)
      call write_log(1, bufout)
      write(bufout,'(a,12x,4f12.4)') 'dt2  =',dtj_inv2(2),dtj_inv2(3),dtj_inv2(4)
      call write_log(1, bufout)
      write(bufout,'(a,    4f12.4)') 'dt3  =',dtj_inv3(1),dtj_inv3(2),dtj_inv3(3),dtj_inv3(4)
      call write_log(1, bufout)
      write(bufout,'(a,36x, f12.4)') 'b1   =', b1(4)
      call write_log(1, bufout)
      write(bufout,'(a,24x,2f12.4)') 'b2   =', b2(3),b2(4)
      call write_log(1, bufout)
      write(bufout,'(a,12x,3f12.4)') 'b3   =', b3(2),b3(3),b3(4)
      call write_log(1, bufout)
      write(bufout,'(a,    4f12.4)') 'b4   =',b4(1),b4(2),b4(3),b4(4)
      call write_log(1, bufout)
      write(bufout,'(a,12x,3f12.4)') 'd1_cy=',d1_cy(2),d1_cy(3),d1_cy(4)
      call write_log(1, bufout)
      write(bufout,'(a,24x,2f12.4)') 'd2_cy=',d2_cy(3),d2_cy(4)
      call write_log(1, bufout)
      write(bufout,'(a,36x, f12.4)') 'd3_cy=',d3_cy(4)
      call write_log(1, bufout)
   endif

end subroutine bspline_get_ppcoef_1seg_modf

!------------------------------------------------------------------------------------------------------------

end module m_bspline_1seg
