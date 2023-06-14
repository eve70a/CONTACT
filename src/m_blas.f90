!------------------------------------------------------------------------------------------------------------
! m_blas - BLAS and BLAS-like functionality
!
! Copyright 1993-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------

module m_blas

implicit none
private

public swap
public reverse_order
public dset
public iset
public dsum
public idmin
public idamin
public idmax
public idamax
public ddot
public dasum
public dnrm2
public daxpy
public dscal
public dcopy

contains

!------------------------------------------------------------------------------------------------------------

      subroutine swap (a, b)
      implicit none
      real(kind=8), intent(inout) :: a, b
      real(kind=8)                :: tmp

      tmp = a
      a   = b
      b   = tmp
      end subroutine swap

!------------------------------------------------------------------------------------------------------------

      subroutine reverse_order (n, a)
      implicit none
      integer,      intent(in)    :: n
      real(kind=8), intent(inout) :: a(n)
      integer                     :: i

      if (n.le.1) return
      do i = 1, n/2
         call swap(a(i), a(n+1-i))
      enddo
      end subroutine reverse_order

!------------------------------------------------------------------------------------------------------------

      subroutine DSet (N, A, Dx, IncX)
!--purpose: Sets all elements of a vector X to the value A
      implicit none
      integer N, IncX, Ix, i
      real(kind=8) A, Dx(*)

      if (N.le.0) Return
      if (IncX.eq.1) then

      ! Code for vector-elements in contiguous memory locations

         do 10 i = 1, N
            Dx (i) = A
 10      continue
      else

      ! Code for  increment not equal to 1

         ix = 1
         if (IncX.lt.0)  ix = (-N+1)*IncX + 1
         do 20 i = 1, N
            Dx (ix) = A
            ix = ix + IncX
 20      continue
      endif
      end subroutine dset

!------------------------------------------------------------------------------------------------------------

      subroutine ISet (N, iA, iX, IncX)
!--purpose: Sets all elements of a vector iX to the value iA

      implicit none
      integer N, IncX, idx, i, iA, iX(1)

      if (N.le.0) Return
      if (IncX.eq.1) then

      ! Code for vector-elements in contiguous memory locations

         do 10 i = 1, N
            iX (i) = iA
 10      continue
      else

      ! Code for  increment not equal to 1

         idx = 1
         if (IncX.lt.0)  idx = (-N+1)*IncX + 1
         do 20 i = 1, N
            iX (idx) = iA
            idx = idx + IncX
 20      continue
      endif
      end subroutine iset

!------------------------------------------------------------------------------------------------------------

      real(kind=8) function dsum(n,dx,incx)
!--purpose: takes the sum of the values in vector dx
!     uses unrolled loops for increment equal to one.
!     modified dasum, by edwin vollebregt 30/7/92.
!     dasum by jack dongarra, linpack, 3/11/78.
!     modified to correct problem with negative increment, 8/21/90.

      real(kind=8) dx(*),dtemp
      integer i,incx,ix,m,mp1,n

      dsum = 0.0d0
      dtemp = 0.0d0
      if (n.le.0) return
      if (incx.eq.1) goto 20

         ! code for increment not equal to 1

      ix = 1
      if (incx.lt.0) ix = (-n+1)*incx + 1
      do 10 i = 1, n
        dtemp = dtemp + dx(ix)
        ix = ix + incx
   10 continue
      dsum = dtemp
      return

         ! code for increment equal to 1

         ! clean-up loop

   20 m = mod(n,6)
      if (m.eq.0) goto 40
      do 30 i = 1, m
        dtemp = dtemp + dx(i)
   30 continue
      if (n.lt.6) goto 60
   40 mp1 = m + 1
      do 50 i = mp1, n, 6
        dtemp = dtemp + dx(i) + dx(i+1) + dx(i+2) + dx(i+3) + dx(i+4) + &
     &                       dx(i+5)
   50 continue
   60 dsum = dtemp
      end function dsum

!------------------------------------------------------------------------------------------------------------

      integer function IdMin (N, dX, IncX)
!--purpose: find the index i with minimum x(i) (not in abs.sense)
      implicit none
      integer IncX, i, ix, N
      real(kind=8) dMin, dX(1)

      IdMin = 0
      if (N.le.0) Return
      if (IncX.eq.1) then
         dMin = dX(1)
         IdMin = 1
         do 10 i = 2, N
            if (dX(i).lt.dMin) then
               dMin = dX(i)
               IdMin = i
            endif
 10      continue
      else
         ix = 1
         if (IncX.lt.0) ix = (-N+1) * IncX + 1
         dMin = dX(ix)
         do 20 i = 2, N
            ix = ix + IncX
            if (dX(ix).lt.dMin) then
               dMin = dX(ix)
               IdMin = ix
            endif
 20      continue
      endif
      end function idmin

!------------------------------------------------------------------------------------------------------------

      integer function IdaMin (N, dX, IncX)
!--purpose: find the in absolute value smallest x(i).
      implicit none
      integer IncX, i, ix, N
      real(kind=8) dMin, dX(1)

      IdaMin = 0
      if (N.le.0) Return
      if (IncX.eq.1) then
         dMin = Abs (dX(1))
         IdaMin = 1
         do 10 i = 2, N
            if (Abs (dX(i)).lt.dMin) then
               dMin = Abs (dX(i))
               IdaMin = i
            endif
 10      continue
      else
         ix = 1
         if (IncX.lt.0) ix = (-N+1) * IncX + 1
         dMin = Abs (dX(ix))
         do 20 i = 2, N
            ix = ix + IncX
            if (Abs (dX(ix)).lt.dMin) then
               dMin = Abs (dX(ix))
               IdaMin = ix
            endif
 20      continue
      endif
      end function idamin

!------------------------------------------------------------------------------------------------------------

      subroutine daxpy(n,da,dx,incx,dy,incy)
!--purpose: constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.

      real(kind=8) dx(*), dy(*), da
      integer i, incx, incy, ix, iy, n

      if (n.le.0) return
      if (da.eq.0.0d0) return
      if (incx.ne.1 .or. incy.ne.1) then

         ! code for unequal increments or equal increments not equal to 1

         ix = 1
         iy = 1
         if (incx.lt.0) ix = (-n+1)*incx + 1
         if (incy.lt.0) iy = (-n+1)*incy + 1
         do 10 i = 1, n
            dy(iy) = dy(iy) + da*dx(ix)
            ix = ix + incx
            iy = iy + incy
   10    continue

      else

         ! code for both increments equal to 1

         do 50 i = 1, n
           dy(i) = dy(i) + da*dx(i)
   50    continue
      endif
      end subroutine daxpy

!------------------------------------------------------------------------------------------------------------

      subroutine  dscal(n,da,dx,incx)
!--purpose: scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified to correct problem with negative increment, 8/21/90.

      real(kind=8) da,dx(1)
      integer i,incx,ix,m,mp1,n

      if (n.le.0)return
      if (incx.eq.1)go to 20

         ! code for increment not equal to 1

      ix = 1
      if (incx.lt.0)ix = (-n+1)*incx + 1
      do 10 i = 1,n
        dx(ix) = da*dx(ix)
        ix = ix + incx
   10 continue
      return

         ! code for increment equal to 1

         ! clean-up loop

   20 m = mod(n,5)
      if (m.eq.0 ) goto 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if (n.lt.5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      end subroutine dscal

!------------------------------------------------------------------------------------------------------------

      integer function IdMax (N, dX, IncX)
!--purpose: find the index i with maximum x(i) (not in abs.sense)
      implicit none
      integer IncX, i, ix, N
      real(kind=8) dMax, dX(1)

      IdMax = 0
      if (N.le.0) Return
      if (IncX.eq.1) then
         dMax = dX(1)
         IdMax = 1
         do 10 i = 2, N
            if (dX(i).gt.dMax) then
               dMax = dX(i)
               IdMax = i
            endif
 10      continue
      else
         ix = 1
         if (IncX.lt.0) ix = (-N+1) * IncX + 1
         dMax = dX(ix)
         do 20 i = 2, N
            ix = ix + IncX
            if (dX(ix).gt.dMax) then
               dMax = dX(ix)
               IdMax = ix
            endif
 20      continue
      endif
      end function idmax

!------------------------------------------------------------------------------------------------------------

      integer function idamax(n,dx,incx)
!--purpose: finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified to correct problem with negative increment, 8/21/90.

      real(kind=8) dx(1),dmax
      integer i,incx,ix,n

      idamax = 0
      if (n.lt.1 ) return
      idamax = 1
      if (n.eq.1)return
      if (incx.eq.1)go to 20

         ! code for increment not equal to 1

      ix = 1
      if (incx.lt.0)ix = (-n+1)*incx + 1
      dmax = dabs(dx(ix))
      ix = ix + incx
      do 10 i = 2,n
         if (dabs(dx(ix)).le.dmax) goto 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return

         ! code for increment equal to 1

   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if (dabs(dx(i)).le.dmax) goto 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      end function idamax

!------------------------------------------------------------------------------------------------------------

      real(kind=8) function ddot(n,dx,incx,dy,incy)
!--purpose: forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.

      real(kind=8) dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n

      ddot = 0.0d0
      dtemp = 0.0d0
      if (n.le.0)return
      if (incx.eq.1 .and. incy.eq.1)go to 20

         ! code for unequal increments or equal increments not equal to 1

      ix = 1
      iy = 1
      if (incx.lt.0)ix = (-n+1)*incx + 1
      if (incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return

         ! code for both increments equal to 1

         ! clean-up loop

   20 m = mod(n,5)
      if (m.eq.0 ) goto 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if (n.lt.5 ) goto 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +             &
     &   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      end function ddot

!------------------------------------------------------------------------------------------------------------

      real(kind=8) function dasum(n,dx,incx)
!--purpose: takes the sum of the absolute values.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified to correct problem with negative increment, 8/21/90.

      real(kind=8) dx(1),dtemp
      integer i,incx,ix,m,mp1,n

      dasum = 0.0d0
      dtemp = 0.0d0
      if (n.le.0)return
      if (incx.eq.1)go to 20

         ! code for increment not equal to 1

      ix = 1
      if (incx.lt.0)ix = (-n+1)*incx + 1
      do 10 i = 1,n
        dtemp = dtemp + dabs(dx(ix))
        ix = ix + incx
   10 continue
      dasum = dtemp
      return

         ! code for increment equal to 1

         ! clean-up loop

   20 m = mod(n,6)
      if (m.eq.0 ) goto 40
      do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
   30 continue
      if (n.lt.6 ) goto 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2)) &
     &  + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 continue
   60 dasum = dtemp
      end function dasum

!------------------------------------------------------------------------------------------------------------

      real(kind=8) function dnrm2 ( n, dx, incx)
      integer i, incx, ix, j, n, next
      real(kind=8)   dx(1), cutlo, cuthi, hitest, sum, xmax,zero,one
      data   zero, one /0.0d0, 1.0d0/

      ! euclidean norm of the n-vector stored in dx() with storage
      ! increment incx .
      ! if    n.le.0 return with result = 0.
      ! if n.ge.1 then incx must be.ge.1
      !
      !       c.l.lawson, 1978 jan 08
      ! modified to correct problem with negative increment, 8/21/90.
      ! modified to correct failure to update ix, 1/25/92.
      !
      ! four phase method     using two built-in constants that are
      ! hopefully applicable to all machines.
      !     cutlo = maximum of  dsqrt(u/eps)  over all known machines.
      !     cuthi = minimum of  dsqrt(v)      over all known machines.
      ! where
      !     eps = smallest no. such that eps + 1..gt.1.
      !     u   = smallest positive no.   (underflow limit)
      !     v   = largest  no.            (overflow  limit)
      !
      ! brief outline of algorithm..
      !
      ! phase 1    scans zero components.
      ! move to phase 2 when a component is nonzero and.le.cutlo
      ! move to phase 3 when a component is.gt.cutlo
      ! move to phase 4 when a component is.ge.cuthi/m
      ! where m = n for x() real and m = 2*n for complex.
      !
      ! values for cutlo and cuthi..
      ! from the environmental parameters listed in the imsl converter
      ! document the limiting values are as follows..
      ! cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
      !               univac and dec at 2**(-103)
      !               thus cutlo = 2**(-51) = 4.44089e-16
      ! cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
      !               thus cuthi = 2**(63.5) = 1.30438e19
      ! cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
      !               thus cutlo = 2**(-33.5) = 8.23181d-11
      ! cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
      ! data cutlo, cuthi / 8.232d-11,  1.304d19 /
      ! data cutlo, cuthi / 4.441e-16,  1.304e19 /

      data cutlo, cuthi / 8.232d-11,  1.304d19 /

      if (n.gt.0) goto 10
         dnrm2  = zero
         goto 300

   10 assign 30 to next
      sum = zero
      i = 1
      if (incx.lt.0 )i = (-n+1)*incx + 1
      ix = 1

                      ! begin main loop

   20    goto next,(30, 50, 70, 110)
   30 if (dabs(dx(i)).gt.cutlo) goto 85
      assign 50 to next
      xmax = zero

                      ! phase 1.  sum is zero

   50 if (dx(i).eq.zero) goto 200
      if (dabs(dx(i)).gt.cutlo) goto 85

                      ! prepare for phase 2.

      assign 70 to next
      goto 105

                      ! phase 2.  sum is small.
                      !         scale to avoid destructive underflow.

   70 if (dabs(dx(i)).gt.cutlo ) goto 75

                      ! common code for phases 2 and 4.
                      ! in phase 4 sum is large.  scale to avoid overflow.

  110 if (dabs(dx(i)).le.xmax ) goto 115
         sum = one + sum * (xmax / dx(i))**2
         xmax = dabs(dx(i))
         goto 200

  115 sum = sum + (dx(i)/xmax)**2
      goto 200

                      ! prepare for phase 3.

   75 sum = (sum * xmax) * xmax


      ! for real or d.p. set hitest = cuthi/n
      ! for complex      set hitest = cuthi/(2*n)

   85 hitest = cuthi/dble( n )

                      ! phase 3.  sum is mid-range.  no scaling.

      do 95 j = ix,n
      if (dabs(dx(i)).ge.hitest) goto 100
         sum = sum + dx(i)**2
         i = i + incx
   95 continue
      dnrm2 = dsqrt( sum )
      goto 300

                      ! prepare for phase 4.

  100 continue
      ix = j
      assign 110 to next
      sum = (sum / dx(i)) / dx(i)
  105 xmax = dabs(dx(i))
      goto 115

  200 continue
      ix = ix + 1
      i = i + incx
      if (ix.le.n ) goto 20

                      ! end of main loop.

      ! compute square root and adjust for scaling.

      dnrm2 = xmax * dsqrt(sum)
  300 continue
      end function dnrm2

!------------------------------------------------------------------------------------------------------------

      subroutine dcopy(n,dx,incx,dy,incy)
!--purpose: copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.

      real(kind=8) dx(1),dy(1)
      integer i,incx,incy,ix,iy,m,mp1,n

      if (n.le.0)return
      if (incx.eq.1 .and. incy.eq.1)go to 20

         ! code for unequal increments or equal increments not equal to 1

      ix = 1
      iy = 1
      if (incx.lt.0)ix = (-n+1)*incx + 1
      if (incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return

         ! code for both increments equal to 1

         ! clean-up loop

   20 m = mod(n,7)
      if (m.eq.0 ) goto 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if (n.lt.7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      end subroutine dcopy

!------------------------------------------------------------------------------------------------------------

end module m_blas

