!------------------------------------------------------------------------------------------------------------
! m_hertz - computation of 2D and 3D Hertz solution
!
! Copyright 1979-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_hertz

use m_hierarch_data
implicit none
private

public  hzsol
private elliptic_to_curv
public  hzcalc2d
public  hzcalc3d
private bisec
public  eli
private ellip_bd
public  equiv_ellipse
public  simpflex
public  linrol

contains

!------------------------------------------------------------------------------------------------------------

subroutine hzsol(ic, geom, pot, hz, kin, mater)
!--purpose: Driver-routine to compute the Hertz-solution, using hzcalc3d directly or through elliptic_to_curv
   implicit none
!--subroutine arguments:
   type(t_ic)       :: ic
   type(t_geomet)   :: geom
   type(t_potcon)   :: pot
   type(t_hertz)    :: hz
   type(t_kincns)   :: kin
   type(t_material) :: mater
!--local variables:
   integer          :: idebug
   real(kind=8)     :: e_star
   real(kind=8), parameter :: epshz = 1e-9

   e_star = mater%ga / (1d0 - mater%nu)

   if (pot%ipotcn.eq.-4 .or. pot%ipotcn.eq.-5) then

      ! Finite line contact: 2d Hertz calculation

      if (ic%norm.ne.1) then
         call write_log('Internal error (Hzsol): 2d Hertz requires N=1')
         call abort_run()
      endif

      call hzcalc2d (pot%ipotcn, hz%a1, hz%aa, hz%b1, hz%bb, hz%cp, hz%rho, kin%fntrue, e_star)

   elseif (pot%ipotcn.eq.-2) then

      ! Curvature+ellipticity prescribed: use an iterative root-finding algorithm

      idebug = 0
      call elliptic_to_curv(ic, hz%a1, hz%aa, hz%aob, hz%b1, hz%bb, hz%cp, hz%rho, kin%fntrue,       &
                     kin%pen, e_star, epshz, idebug)

   else

      ! Point contact, curvatures or semi-axes prescribed: 3d Hertz calculation

      call hzcalc3d (e_star, epshz, pot%ipotcn, hz%a1, hz%b1, hz%aa, hz%bb, ic%norm, kin%pen,           &
                     kin%fntrue, hz%cp, hz%rho)

   endif

   ! Compute aa/bb, allowing for aa, bb >= 0

   if (pot%ipotcn.ne.-2) then
      if (hz%aa.le.tiny .and. hz%bb.le.tiny) then
         hz%aob = 1d0
      else
         hz%aob = hz%aa / hz%bb
      endif
   endif

   ! Store solution variables

   geom%ibase       = 1
   geom%prmudf(1:6) = 0d0
   geom%prmudf(1)   = hz%a1
   geom%prmudf(3)   = hz%b1

   if (ic%norm.eq.0) kin%fnscal = kin%fntrue / mater%ga

end subroutine hzsol

!------------------------------------------------------------------------------------------------------------

subroutine elliptic_to_curv(ic, a1, aa, aob, b1, bb, cp, rho, fntrue, pen, e_star, epshz, idebug)
!--purpose: For given curvature a1 and ellipticity ratio aob determine the corresponding curvature b1
!           using an iterative root-finding algorithm.
   implicit none
!--subroutine arguments:
   type(t_ic)   :: ic
   integer      :: idebug
   real(kind=8) :: a1, aa, aob, b1, bb, cp, rho, fntrue, pen, e_star, epshz
!--local variables:
   integer iter, ipotcn
   real(kind=8) aoblow, aobhigh, daob, b1low, b1high, epsaob, rico

   epsaob = 10d0*epshz
   if (idebug.ge.1) then
      write(bufout,*) 'Special calculation of curvature b1 from given A/B =', aob
      call write_log(1, bufout)
   endif

   ! initialize variables for call of hzcalc3d

   ipotcn = -1

   ! Set initial estimate b1 := a1
   ! Solve the Hertz-problem for aa and bb, and compute initial aa/bb

   iter  = 0
   b1 = a1
   call hzcalc3d (e_star, epshz, ipotcn, a1, b1, aa, bb, ic%norm, pen, fntrue, cp, rho)

   if (aa.lt.tiny .and. bb.lt.tiny) then
      if (idebug.ge.1) call write_log('Hzcalc3d: empty contact area')
      return
   endif

   if (idebug.ge.2) then
      write(bufout, 19) 'Iter0', iter, 0., 99.99d0, b1, aa/bb
      call write_log(1, bufout)
19    format(a,i4,': B in [',f12.9,',',f12.9,'], b1=',f12.9, ', A/B=',f12.9)
   endif

   ! determine a second value for b1 such that desired A/B is included in an interval [b1low, b1high]

   if (aa/bb.gt.aob) then

      ! initial ellipticity too large: curvature b1 too large:
      !    use successive halvation of b1 until ellipticity is too small

      b1high  = b1
      aobhigh = aa/bb
 20   if (iter.le.100 .and. aa/bb.gt.aob) then
         iter    = iter + 1
         b1high  = b1
         aobhigh = aa/bb
         b1      = b1*0.5d0
         call hzcalc3d (e_star, epshz, ipotcn, a1, b1, aa, bb, ic%norm, pen, fntrue, cp, rho)
         if (idebug.ge.2) then
            write(bufout, 19) 'Iter1', iter, 0., b1high, b1, aa/bb
            call write_log(1, bufout)
         endif
         goto 20
      endif
      b1low  = b1
      aoblow = aa/bb
   else

      ! initial ellipticity too small: curvature b1 too small:
      !    use successive doubling of b1 until ellipticity is too small

      b1low  = b1
      aoblow = aa/bb
 25   if (iter.le.100 .and. aa/bb.lt.aob) then
         iter   = iter + 1
         b1low  = b1
         aoblow = aa/bb
         b1     = b1*2.0d0
         call hzcalc3d (e_star, epshz, ipotcn, a1, b1, aa, bb, ic%norm, pen, fntrue, cp, rho)
         if (idebug.ge.2) then
            write(bufout, 19) 'Iter2', iter, b1low, 99.99d0, b1, aa/bb
            call write_log(1, bufout)
         endif
         goto 25
      endif
      b1high  = b1
      aobhigh = aa/bb
   endif

   ! iteratively refine the interval using "Regula Falsi":
   !    select a value b1 in the interval using linear interpolation,
   !   compute A/B, set new interval to either [b1low,b1] or [b1,b1high]

 30 if (iter.le.100 .and. dabs(aa/bb-aob).gt.epsaob*aob) then
      iter = iter + 1

      daob = aob - aoblow
      rico = (aobhigh - aoblow) / (b1high - b1low)
      b1   = b1low + daob / rico

      call hzcalc3d (e_star, epshz, ipotcn, a1, b1, aa, bb, ic%norm, pen, fntrue, cp, rho)
      if (idebug.ge.2) then
         write(bufout, 19) 'Iter3', iter, b1low, b1high, b1, aa/bb
         call write_log(1, bufout)
      endif
      if (aa/bb.gt.aob) then

         ! ellipticity at b1 too large: interval [b1low,b1]

         b1high = b1
         aobhigh = aa/bb
      else

         ! ellipticity at b1 too small: interval [b1,b1high]

         b1low = b1
         aoblow = aa/bb
      endif
      goto 30
   endif

end subroutine elliptic_to_curv

!------------------------------------------------------------------------------------------------------------

subroutine hzcalc2d (ipotcn, a1, aa, b1, bb, cp, rho, fntrue, e_star)
!--purpose: Solve the Hertz-problem for 2D line contacts.
   implicit none
!--subroutine parameters:
   integer      :: ipotcn
   real(kind=8) :: a1, b1, aa, bb, cp, rho, fntrue, e_star

   ! ipotcn    type of calculation, must be -4 (a1,bb specified) or -5 (aa,bb) specified
   ! a1        combined curvature in rolling direction
   ! b1        combined curvature in transverse direction
   ! aa        semi-length of contact rectangle in rolling direction
   ! bb        semi-width of contact rectangle in transverse direction
   ! cp        effective contact radius, area 4*aa*bb == pi*cp^2
   ! e_star    ga / (1 - nu)
   ! rho       effective radius of curvature 2 / (a1 + b1)
   ! fntrue    total normal force
!--local variables:
   real(kind=8) :: fn_2d, rx

   fn_2d  = fntrue / (2d0 * bb)
   b1     = 0d0

   if (ipotcn.eq.-4) then

      ! curvature a1 prescribed, normal force fntrue

      rx     = 1d0 / (2d0 * a1)
      aa     = sqrt( 4d0 * rx * fn_2d / (pi * e_star) )

   elseif (ipotcn.eq.-5) then

      ! semi-length aa prescribed, normal force fntrue

      rx     = pi * aa**2 * e_star / (4d0 * fn_2d)
      a1     = 1d0 / (2d0 * rx)

   else

      call write_log('Internal error (hertz2d): invalid ipotcn')
      call abort_run()

   endif

   rho = 2d0 / (a1 + b1)
   cp  = sqrt(4d0 * aa * bb / pi)

end subroutine hzcalc2d

!------------------------------------------------------------------------------------------------------------

subroutine hzcalc3d (e_star, epshz, ipotcn, a1, b1, aa, bb, ic_norm, pen, fntrue, cp, rho)
!--purpose: Solve the Hertz-problem for 3D point contacts.
!           Note: the case ipotcn=-2, given a/b, is handled outside this routine.
   implicit none
!--subroutine parameters:
   integer,      intent(in)    :: ic_norm, ipotcn
   real(kind=8), intent(in)    :: e_star, epshz
   real(kind=8), intent(out)   :: cp, rho
   real(kind=8), intent(inout) :: a1, b1, aa, bb, fntrue, pen

   ! ipotcn    type of calculation, -1 = a1,b1 given, aa,bb output; -3 = aa,bb given, a1,b1 output
   ! ic_norm   type of calculation, 0 = pen given, fntrue output; 1 = fntrue given, pen output
   ! a1        combined curvature in rolling direction
   ! b1        combined curvature in transverse direction
   ! aa        semi-axis of contact ellipse in rolling direction
   ! bb        semi-axis of contact ellipse in transverse direction
   ! cp        sqrt(aa*bb), area of contact ellipse: pi*cp^2=pi*aa*bb
   ! e_star    ga / (1 - nu)
   ! rho       effective radius of curvature 2 / (a1 + b1)
   ! fntrue    total normal force
!--local variables:
   real(kind=8), parameter :: penmin = 1d-9
   real(kind=8) ama, ami, apb, elmodk, ellipk, ellipe, g, sqrtg, x, y
   logical zbla

   ! g         ratio of minimum over maximum semi-axes
   ! elmodk    elliptic modulus k of equation (1.57a), k = sqrt(1 - g^2),  g = min(a,b) / max(a,b)

   if (ipotcn.eq.-1) then

      ! curvatures prescribed in a1, b1

      if (min(a1,b1).lt.tiny) then
         write(bufout,'(2(a,g12.4),a)') ' ERROR (Hertz3D): curvatures A=',a1,' and B=',b1,              &
                ' must be positive.'
         call write_log(1, bufout)
      endif

      zbla = (b1.le.a1)
      if (zbla) then
         y = b1/a1
      else
         y = a1/b1
      endif

      ! find point x with eli(x) = y in the interval [0, 1]

      call bisec (ellipk, ellipe, epshz, x, 0.d0, 1.d0, y)

      elmodk = x
      g = dsqrt(max(1d-40, 1d0-elmodk**2))
      sqrtg = dsqrt(g)
      rho = 2d0 / (a1 + b1)

      ! compute cp,pen when fntrue is known, or cp,fntrue when pen is known
      ! note that aa^2*bb == cp^3*sqrt(g), aa^2 = (cp*sqrt(g))^2

      if (ic_norm.eq.1) then
         cp  = ((3d0 * max(0d0,fntrue) * rho * ellipe) / (4d0 * pi * e_star * sqrtg)) **(1d0/3d0)
         pen = 2d0 * (cp * sqrtg)**2 * ellipk / (rho * ellipe)
      else
         cp  = sqrt((max(0d0,pen) * rho * ellipe) / (2d0 * ellipk * sqrtg**2))
         fntrue = 4d0 * pi * cp**3 * e_star * sqrtg / (3d0 * rho * ellipe)
      endif

      ! determine corresponding semi-axes aa, bb

      if (zbla) then
         aa = cp * sqrtg
         bb = cp / sqrtg
      else
         aa = cp / sqrtg
         bb = cp * sqrtg
      endif

   else

      ! semiaxes prescribed in aa, bb

      zbla = (bb.le.aa)
      if (zbla) then
         g = bb/aa
      else
         g = aa/bb
      endif
      sqrtg  = sqrt (g)
      elmodk = sqrt (1d0 - g*g)
      y = eli (elmodk, ellipk, ellipe)
      cp = sqrt (aa*bb)

      ! compute rho,pen when fntrue is known, or rho,fntrue when pen is known

      if (ic_norm.eq.1) then
         rho    = 4d0 * pi * cp**3 * e_star * sqrtg / (3d0 * fntrue * ellipe)
         pen    = 2d0 * (cp * sqrtg)**2 * ellipk / (rho * ellipe)
      else
         pen    = max(pen, penmin)
         rho    = 2d0 * (cp * sqrtg)**2 * ellipk / (pen * ellipe)
         fntrue = 4d0 * pi * cp**3 * e_star * sqrtg / (3d0 * rho * ellipe)
      endif

      ! compute major and minor curvatures, copy to a1, b1

      apb = 2d0 / rho
      ama = apb / (y + 1d0)
      ami = y * ama
      if (zbla) then
         a1 = ami
         b1 = ama
      else
         a1 = ama
         b1 = ami
      endif
   endif

end subroutine hzcalc3d

!------------------------------------------------------------------------------------------------------------

subroutine bisec(ellipk, ellipe, epshz, x, xlft, xrght, y)
!--purpose: Simple bisection routine to find the value of x (elliptic modulus k) between xlft and xrght
!           such that the ratio of curvatures D1/D2 == eli(x) equals y. 
!           The values of ellipk and ellipe are a side-product of function eli, and are passed
!           through to the calling routine.
   implicit none
!--subroutine parameters:
   real(kind=8), intent(in)  :: epshz, xlft, xrght, y
   real(kind=8), intent(out) :: ellipk, ellipe, x
!--local variables:
   real(kind=8) :: ell, elr, elx, xl, xr

   ! initialize left/right points

   xl = xlft
   xr = xrght
   elr = eli (xr, ellipk, ellipe)
   ell = eli (xl, ellipk, ellipe)

   ! iterate until requested precision is reached

   do while (abs(xr-xl).gt.epshz)

      ! verify that y is between elr and ell

      if ((y-elr)*(y-ell).gt.0d0) then
         write (*,10) elr,y,ell
 10      format (' Internal error (Bisec): Y should be between ELR and ELL.', ' ELR,Y,ELL:',3g12.4)
         call abort_run()
      endif

      ! compute point half-way

      x = (xl + xr) / 2d0
      elx = eli(x, ellipk, ellipe)

      ! update either lower or upper boundary

      if ((y-elr)*(y-elx).le.0d0) then

         ! root is situated between xr and x

         xl = x
         ell = elx
      else

         ! root is situated between xl and x

         xr = x
         elr = elx
      endif
   enddo

   ! obtain final values for ellipk, ellipe

   elx = eli (x, ellipk, ellipe)
end subroutine bisec

!------------------------------------------------------------------------------------------------------------

function eli (elmodk, ell_k, ell_e)
!--purpose: routine that determines the complete elliptic integrals "K" and "E" (see [Kalker1990] App.E).
!           The function result is a combination of these elliptic integrals such that when e is
!           the eccentricity of the contact ellipse, eli(e)= min(A1/B1,B1/A1), where A1 and B1 are the
!           curvatures of the undeformed distance function h(x,y)=A1*x**2 + B1*y**2 - pen.
   implicit none
!--function result:
   real(kind=8) :: eli
!--subroutine arguments:
   real(kind=8), intent(in)  :: elmodk
   real(kind=8), intent(out) :: ell_k, ell_e

   ! elmodk    : elliptic modulus k of equation (1.57a), k = sqrt(1 - g^2),  g = min(a,b) / max(a,b)
   ! ell_k     : complete elliptic integral of the first kind K
   ! ell_e     : complete elliptic integral of the second kind E
   ! eli       : ratio of curvatures D_rel = D_min/D_max, D_min=D1=min(a1,b1), D_max=D2=max(a1,b1)
!--local variables:
   real(kind=8) mc, ell_b, ell_d 
   real(kind=8) e2, fln, g2, u, elldmc
   ! mc        : "complementary modulus" cf. Fukushima, square of g, g = sqrt(1-k^2)

   if (.true.) then
      mc    = 1d0 - elmodk**2
      call ellip_bd(mc, ell_b, ell_d)
      ell_k = ell_b + ell_d
      ell_e = ell_b + mc * ell_d
      eli   = (ell_e - ell_b) / ell_b
   else

      e2 = elmodk**2
      g2 = 1.d0 - e2
      if (e2.lt.0.05) fln = g2 * (1/2d0 + e2/3d0 + e2**2/4d0 + e2**3/5d0 + e2**4/6d0)
      if (e2.ge.0.05) fln = g2 * (dlog (1/(g2+1d-20)) - e2) / e2**2

      ! ell_k IS THE COMPLETE EllIPTIC INTEGRAL OF THE 1ST KIND "K".
      ! K = \int_0^{\pi/2} 1 / sqrt(1 - k^2 sin^2(\psi)) d\psi

      ell_k = 1.3862944 + 0.1119723*g2 + 0.0725296*g2**2 - ( 0.5 + 0.1213478*g2 + 0.0288729*g2**2 )    &
                                                                                      * dlog (g2+1d-20)

      ! ell_e IS THE COMPLETE EllIPTIC INTEGRAL OF THE 2ND KIND "E".
      ! E = \int_0^{\pi/2} sqrt(1 - k^2 sin^2(\psi)) d\psi

      ell_e = 1 + 0.4630151*g2 + 0.1077812*g2**2 + ( 0.2452727*g2 + 0.0412496*g2**2) * dlog (1/(g2+1d-20))

      ! elldmc IS THE COMPLETE EllIPTICAL INTEGRAL "D-C", SEE /1/

      elldmc = 1 - 0.064301*g2 + 0.0168037*g2**2 + (-0.7547273 + 0.0438267*g2 - 0.0164962*g2**2) * fln

      u = e2 * elldmc / ell_e
      eli = (1d0 - u) / (1d0 + u)
   endif

end function eli

!------------------------------------------------------------------------------------------------------------

subroutine ellip_bd(mc, elb, eld)
!
!  Double precision associate complete elliptic integral of the second kind
!
!  Reference: T. Fukushima, (2011), Math. Comp., 80, 1725-1743
!     "Precise and Fast Computation of General Complete Elliptic Integral of Second Kind"
!
!  Author: T. Fukushima Toshio.Fukushima@nao.ac.jp
!
!  Inputs: mc = complementary parameter, 0 < mc <= 1
!
!  Outputs: elb = B(m), eld = D(m)
!
   implicit none
   real*8 mc,elb,eld
   real*8 m,nome,dkkc,dddc,mx,kkc,logq2,elk,dele,elk1,delb

   real*8 PIQ,PIHALF,PI,PIINV
   parameter (PIQ=0.78539816339744830961566084581988d0)
   parameter (PIHALF=1.5707963267948966192313216916398d0)
   parameter (PI=3.1415926535897932384626433832795d0)
   parameter (PIINV=0.31830988618379067153776752674503d0)

   real*8 mcold,elbold,eldold
   save mcold,elbold,eldold

   real*8 Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,Q16
   parameter (Q1=1.d0/16.d0,Q2=1.d0/32.d0,Q3=21.d0/1024.d0)
   parameter (Q4=31.d0/2048.d0,Q5=6257.d0/524288.d0)
   parameter (Q6=10293.d0/1048576.d0,Q7=279025.d0/33554432.d0)
   parameter (Q8=483127.d0/67108864.d0)
   parameter (Q9=435506703.d0/68719476736.d0)
   parameter (Q10=776957575.d0/137438953472.d0)
   parameter (Q11=22417045555.d0/4398046511104.d0)
   parameter (Q12=40784671953.d0/8796093022208.d0)
   parameter (Q13=9569130097211.d0/2251799813685248.d0)
   parameter (Q14=17652604545791.d0/4503599627370496.d0)
   parameter (Q15=523910972020563.d0/144115188075855872.d0)
   parameter (Q16=976501268709949.d0/288230376151711744.d0)

   real*8 K1,K2,K3,K4,K5,K6,K7
   parameter (K1=1.d0/4.d0)
   parameter (K2=9.d0/64.d0)
   parameter (K3=25.d0/256.d0)
   parameter (K4=1225.d0/16384.d0)
   parameter (K5=3969.d0/65536.d0)
   parameter (K6=53361.d0/1048576.d0)
   parameter (K7=184041.d0/4194304.d0)

   real*8 B1,B2,B3,B4,B5,B6,B7,B8
   parameter (B1=1.d0/2.d0)
   parameter (B2=1.d0/16.d0)
   parameter (B3=3.d0/128.d0)
   parameter (B4=25.d0/2048.d0)
   parameter (B5=245.d0/32768.d0)
   parameter (B6=1323.d0/262144.d0)
   parameter (B7=7623.d0/2097152.d0)
   parameter (B8=184041.d0/67108864.d0)

   real*8 D1,D2,D3,D4,D5,D6,D7,D8
   parameter (D1=1.d0/2.d0)
   parameter (D2=3.d0/16.d0)
   parameter (D3=15.d0/128.d0)
   parameter (D4=175.d0/2048.d0)
   parameter (D5=2205.d0/32768.d0)
   parameter (D6=14553.d0/262144.d0)
   parameter (D7=99099.d0/2097152.d0)
   parameter (D8=2760615.d0/67108864.d0)

   logical first/.TRUE./

   if(first) then
      first=.FALSE.
      mcold=1.d0
      elbold=PIQ
      eldold=PIQ
   endif

   m=1.d0-mc

   if(abs(mc-mcold).lt.1.11d-16*mc) then
      elb=elbold
      eld=eldold
   elseif(m.lt.1.11d-16) then
      elb=PIQ
      eld=PIQ
   elseif(mc.lt.1.11d-32) then       ! EV: add guard against Inf
      elb=1.d0
      eld=0.3862943611198906188344642429164d0-0.5d0*log(1.11d-32)
   elseif(mc.lt.1.11d-16) then
      elb=1.d0
      eld=0.3862943611198906188344642429164d0-0.5d0*log(mc)
   elseif(mc.lt.0.1d0) then
      nome=mc*(Q1+mc*(Q2+mc*(Q3+mc*(Q4+mc*(Q5+mc*(Q6 &
          +mc*(Q7+mc*(Q8+mc*(Q9+mc*(Q10+mc*(Q11+mc*(Q12 &
          +mc*(Q13+mc*(Q14+mc*(Q15+mc*Q16))))))))))))))) 
      if(mc.lt.0.01d0) then
         dkkc=mc*(K1+mc*(K2+mc*(K3+mc*(K4+mc*(K5+mc*(K6+mc*K7))))))
         dddc=mc*(D1+mc*(D2+mc*(D3+mc*(D4+mc*(D5+mc*(D6+mc*D7))))))
      else
         mx=mc-0.05d0

         ! (K'-1)/(pi/2)

         dkkc=    0.01286425658832983978282698630501405107893d0 &
            +mx*(0.26483429894479586582278131697637750604652d0 &
            +mx*(0.15647573786069663900214275050014481397750d0 &
            +mx*(0.11426146079748350067910196981167739749361d0 &
            +mx*(0.09202724415743445309239690377424239940545d0 &
            +mx*(0.07843218831801764082998285878311322932444d0 &
            +mx*(0.06935260142642158347117402021639363379689d0 &
            +mx*(0.06293203529021269706312943517695310879457d0 &
            +mx*(0.05821227592779397036582491084172892108196d0 &
            +mx*(0.05464909112091564816652510649708377642504d0 &
            +mx*(0.05191068843704411873477650167894906357568d0 &
            +mx*(0.04978344771840508342564702588639140680363d0 &
            +mx*(0.04812375496807025605361215168677991360500d0 &
            ))))))))))))

         ! (K'-E')/(pi/2)

         dddc=    0.02548395442966088473597712420249483947953d0 &
            +mx*(0.51967384324140471318255255900132590084179d0 &
            +mx*(0.20644951110163173131719312525729037023377d0 &
            +mx*(0.13610952125712137420240739057403788152260d0 &
            +mx*(0.10458014040566978574883406877392984277718d0 &
            +mx*(0.08674612915759188982465635633597382093113d0 &
            +mx*(0.07536380269617058326770965489534014190391d0 &
            +mx*(0.06754544594618781950496091910264174396541d0 &
            +mx*(0.06190939688096410201497509102047998554900d0 &
            +mx*(0.05771071515451786553160533778648705873199d0 &
            +mx*(0.05451217098672207169493767625617704078257d0 &
            +mx*(0.05204028407582600387265992107877094920787d0 &
            +mx*(0.05011532514520838441892567405879742720039d0 &
            ))))))))))))
      endif
      kkc=1.d0+dkkc
      logq2=-0.5d0*log(nome)
      elk=kkc*logq2
      dele=-dkkc/kkc+logq2*dddc
      elk1=elk-1.d0
      delb=(dele-mc*elk1)/m
      elb=1.d0+delb
      eld=elk1-delb
   elseif(m.le.0.01d0) then
      elb=PIHALF*(B1+m*(B2+m*(B3+m*(B4+m*(B5+m*(B6+m*(B7+m*B8)))))))
      eld=PIHALF*(D1+m*(D2+m*(D3+m*(D4+m*(D5+m*(D6+m*(D7+m*D8)))))))
   elseif(m.le.0.1d0) then
      mx=0.95d0-mc
      elb=     0.790401413584395132310045630540381158921005d0 &
         +mx*(0.102006266220019154892513446364386528537788d0 &
         +mx*(0.039878395558551460860377468871167215878458d0 &
         +mx*(0.021737136375982167333478696987134316809322d0 &
         +mx*(0.013960979767622057852185340153691548520857d0 &
         +mx*(0.009892518822669142478846083436285145400444d0 &
         +mx*(0.007484612400663335676130416571517444936951d0 &
         +mx*(0.005934625664295473695080715589652011420808d0 &
         +mx*(0.004874249053581664096949448689997843978535d0 &
         +mx*(0.004114606930310886136960940893002069423559d0 &
         +mx*(0.003550452989196176932747744728766021440856d0 &
         +mx*(0.003119229959988474753291950759202798352266d0 &
         )))))))))))
      eld=     0.800602040206397047799296975176819811774784d0 &
         +mx*(0.313994477771767756849615832867393028789057d0 &
         +mx*(0.205913118705551954501930953451976374435088d0 &
         +mx*(0.157744346538923994475225014971416837073598d0 &
         +mx*(0.130595077319933091909091103101366509387938d0 &
         +mx*(0.113308474489758568672985167742047066367053d0 &
         +mx*(0.101454199173630195376251916342483192174927d0 &
         +mx*(0.0929187842072974367037702927967784464949434d0 &
         +mx*(0.0865653801481680871714054745336652101162894d0 &
         +mx*(0.0817279846651030135350056216958053404884715d0 &
         +mx*(0.0779906657291070378163237851392095284454654d0 &
         +mx*(0.075080426851268007156477347905308063808848d0 &
         )))))))))))
   elseif(m.le.0.2d0) then
      mx=0.85d0-mc
      elb=     0.80102406445284489393880821604009991524037d0 &
         +mx*(0.11069534452963401497502459778015097487115d0 &
         +mx*(0.047348746716993717753569559936346358937777d0 &
         +mx*(0.028484367255041422845322166419447281776162d0 &
         +mx*(0.020277811444003597057721308432225505126013d0 &
         +mx*(0.015965005853099119442287313909177068173564d0 &
         +mx*(0.013441320273553634762716845175446390822633d0 &
         +mx*(0.011871565736951439501853534319081030547931d0 &
         +mx*(0.010868363672485520630005005782151743785248d0 &
         +mx*(0.010231587232710564565903812652581252337697d0 &
         +mx*(0.009849585546666211201566987057592610884309d0 &
         +mx*(0.009656606347153765129943681090056980586986d0 &
         )))))))))))
      eld=     0.834232667811735098431315595374145207701720d0 &
         +mx*(0.360495281619098275577215529302260739976126d0 &
         +mx*(0.262379664114505869328637749459234348287432d0 &
         +mx*(0.223723944518094276386520735054801578584350d0 &
         +mx*(0.206447811775681052682922746753795148394463d0 &
         +mx*(0.199809440876486856438050774316751253389944d0 &
         +mx*(0.199667451603795274869211409350873244844882d0 &
         +mx*(0.204157558868236842039815028663379643303565d0 &
         +mx*(0.212387467960572375038025392458549025660994d0 &
         +mx*(0.223948914061499360356873401571821627069173d0 &
         +mx*(0.238708097425597860161720875806632864507536d0 &
         +mx*(0.256707203545463755643710021815937785120030d0 &
         )))))))))))
   elseif(m.le.0.3d0) then
      mx=0.75d0-mc
      elb=     0.81259777291992049322557009456643357559904d0 &
         +mx*(0.12110961794551011284012693733241967660542d0 &
         +mx*(0.057293376831239877456538980381277010644332d0 &
         +mx*(0.038509451602167328057004166642521093142114d0 &
         +mx*(0.030783430301775232744816612250699163538318d0 &
         +mx*(0.027290564934732526869467118496664914274956d0 &
         +mx*(0.025916369289445198731886546557337255438083d0 &
         +mx*(0.025847203343361799141092472018796130324244d0 &
         +mx*(0.026740923539348854616932735567182946385269d0 &
         +mx*(0.028464314554825704963640157657034405579849d0 &
         +mx*(0.030995446237278954096189768338119395563447d0 &
         +mx*(0.034384369179940975864103666824736551261799d0 &
         +mx*(0.038738002072493935952384233588242422046537d0 &
         ))))))))))))
      eld=     0.873152581892675549645633563232643413901757d0 &
         +mx*(0.420622230667770215976919792378536040460605d0 &
         +mx*(0.344231061559450379368201151870166692934830d0 &
         +mx*(0.331133021818721761888662390999045979071436d0 &
         +mx*(0.345277285052808411877017306497954757532251d0 &
         +mx*(0.377945322150393391759797943135325823338761d0 &
         +mx*(0.427378012464553880508348757311170776829930d0 &
         +mx*(0.494671744307822405584118022550673740404732d0 &
         +mx*(0.582685115665646200824237214098764913658889d0 &
         +mx*(0.695799207728083164790111837174250683834359d0 &
         +mx*(0.840018401472533403272555302636558338772258d0 &
         +mx*(1.023268503573606060588689738498395211300483d0 &
         +mx*(1.255859085136282496149035687741403985044122d0 &
         ))))))))))))
   elseif(m.le.0.4d0) then
      mx=0.65d0-mc
      elb=     0.8253235579835158949845697805395190063745d0 &
         +mx*(0.1338621160836877898575391383950840569989d0 &
         +mx*(0.0710112935979886745743770664203746758134d0 &
         +mx*(0.0541784774173873762208472152701393154906d0 &
         +mx*(0.0494517449481029932714386586401273353617d0 &
         +mx*(0.0502221962241074764652127892365024315554d0 &
         +mx*(0.0547429131718303528104722303305931350375d0 &
         +mx*(0.0627462579270016992000788492778894700075d0 &
         +mx*(0.0746698810434768864678760362745179321956d0 &
         +mx*(0.0914808451777334717996463421986810092918d0 &
         +mx*(0.1147050921109978235104185800057554574708d0 &
         +mx*(0.1465711325814398757043492181099197917984d0 &
         +mx*(0.1902571373338462844225085057953823854177d0 &
         ))))))))))))
      eld=     0.9190270392420973478848471774160778462738d0 &
         +mx*(0.5010021592882475139767453081737767171354d0 &
         +mx*(0.4688312705664568629356644841691659415972d0 &
         +mx*(0.5177142277764000147059587510833317474467d0 &
         +mx*(0.6208433913173031070711926900889045286988d0 &
         +mx*(0.7823643937868697229213240489900179142670d0 &
         +mx*(1.0191145350761029126165253557593691585239d0 &
         +mx*(1.3593452027484960522212885423056424704073d0 &
         +mx*(1.8457173023588279422916645725184952058635d0 &
         +mx*(2.5410717031539207287662105618152273788399d0 &
         +mx*(3.5374046552080413366422791595672470037341d0 &
         +mx*(4.9692960029774259303491034652093672488707d0 &
         +mx*(7.0338228700300311264031522795337352226926d0 &
         +mx*(10.020043225034471401553194050933390974016d0 &
         )))))))))))))
   elseif(m.le.0.5d0) then
      mx=0.55d0-mc
      elb=     0.8394795702706129706783934654948360410325d0 &
         +mx*(0.1499164403063963359478614453083470750543d0 &
         +mx*(0.0908319358194288345999005586556105610069d0 &
         +mx*(0.0803470334833417864262134081954987019902d0 &
         +mx*(0.0856384405004704542717663971835424473169d0 &
         +mx*(0.1019547259329903716766105911448528069506d0 &
         +mx*(0.1305748115336160150072309911623351523284d0 &
         +mx*(0.1761050763588499277679704537732929242811d0 &
         +mx*(0.2468351644029554468698889593583314853486d0 &
         +mx*(0.3564244768677188553323196975301769697977d0 &
         +mx*(0.5270025622301027434418321205779314762241d0 &
         +mx*(0.7943896342593047502260866957039427731776d0 &
         +mx*(1.2167625324297180206378753787253096783993d0 &
         ))))))))))))
      eld=     0.9744043665463696730314687662723484085813d0 &
         +mx*(0.6132468053941609101234053415051402349752d0 &
         +mx*(0.6710966695021669963502789954058993004082d0 &
         +mx*(0.8707276201850861403618528872292437242726d0 &
         +mx*(1.2295422312026907609906452348037196571302d0 &
         +mx*(1.8266059675444205694817638548699906990301d0 &
         +mx*(2.8069345309977627400322167438821024032409d0 &
         +mx*(4.4187893290840281339827573139793805587268d0 &
         +mx*(7.0832360574787653249799018590860687062869d0 &
         +mx*(11.515088120557582942290563338274745712174d0 &
         +mx*(18.931511185999274639516732819605594455165d0 &
         +mx*(31.411996938204963878089048091424028309798d0 &
         +mx*(52.520729454575828537934780076768577185134d0 &
         +mx*(88.384854735065298062125622417251073520996d0 &
         +mx*(149.56637449398047835236703116483092644714d0 &
         +mx*(254.31790843104117434615624121937495622372d0 &
         )))))))))))))))
   elseif(m.le.0.6d0) then
      mx=0.45d0-mc
      elb=     0.8554696151564199914087224774321783838373d0 &
         +mx*(0.1708960726897395844132234165994754905373d0 &
         +mx*(0.1213352290269482260207667564010437464156d0 &
         +mx*(0.1282018835749474096272901529341076494573d0 &
         +mx*(0.1646872814515275597348427294090563472179d0 &
         +mx*(0.2374189087493817423375114793658754489958d0 &
         +mx*(0.3692081047164954516884561039890508294508d0 &
         +mx*(0.6056587338479277173311618664015401963868d0 &
         +mx*(1.0337055615578127436826717513776452476106d0 &
         +mx*(1.8189884893632678849599091011718520567105d0 &
         +mx*(3.2793776512738509375806561547016925831128d0 &
         +mx*(6.0298883807175363312261449542978750456611d0 &
         +mx*(11.269796855577941715109155203721740735793d0 &
         +mx*(21.354577850382834496786315532111529462693d0 &
         )))))))))))))
      eld=     1.04345529511513353426326823569160142342838d0 &
         +mx*(0.77962572192850485048535711388072271372632d0 &
         +mx*(1.02974236093206758187389128668777397528702d0 &
         +mx*(1.62203722341135313022433907993860147395972d0 &
         +mx*(2.78798953118534762046989770119382209443756d0 &
         +mx*(5.04838148737206914685643655935236541332892d0 &
         +mx*(9.46327761194348429539987572314952029503864d0 &
         +mx*(18.1814899494276679043749394081463811247757d0 &
         +mx*(35.5809805911791687037085198750213045708148d0 &
         +mx*(70.6339354619144501276254906239838074917358d0 &
         +mx*(141.828580083433059297030133195739832297859d0 &
         +mx*(287.448751250132166257642182637978103762677d0 &
         +mx*(587.115384649923076181773192202238389711345d0 &
         +mx*(1207.06543522548061603657141890778290399603d0 &
         +mx*(2495.58872724866422273012188618178997342537d0 &
         +mx*(5184.69242939480644062471334944523925163600d0 &
         +mx*(10817.2133369041327524988910635205356016939d0 &
         ))))))))))))))))
   elseif(m.le.0.7d0) then
      mx=0.35d0-mc
      elb=     0.8739200618486431359820482173294324246058d0 &
         +mx*(0.1998140574823769459497418213885348159654d0 &
         +mx*(0.1727696158780152128147094051876565603862d0 &
         +mx*(0.2281069132842021671319791750725846795701d0 &
         +mx*(0.3704681411180712197627619157146806221767d0 &
         +mx*(0.6792712528848205545443855883980014994450d0 &
         +mx*(1.3480084966817573020596179874311042267679d0 &
         +mx*(2.8276709768538207038046918250872679902352d0 &
         +mx*(6.1794682501239140840906583219887062092430d0 &
         +mx*(13.935686010342811497608625663457407447757d0 &
         +mx*(32.218929281059722026322932181420383764028d0 &
         +mx*(76.006962959226101026399085304912635262362d0 &
         +mx*(182.32144908775406957609058046006949657416d0 &
         +mx*(443.51507644112648158679360783118806161062d0 &
         +mx*(1091.8547229028388292980623647414961662223d0 &
         +mx*(2715.7658664038195881056269799613407111521d0 &
         )))))))))))))))
      eld=     1.13367833657573316571671258513452768536080d0 &
         +mx*(1.04864317372997039116746991765351150490010d0 &
         +mx*(1.75346504119846451588826580872136305225406d0 &
         +mx*(3.52318272680338551269021618722443199230946d0 &
         +mx*(7.74947641381397458240336356601913534598302d0 &
         +mx*(17.9864500558507330560532617743406294626849d0 &
         +mx*(43.2559163462326133313977294448984936591235d0 &
         +mx*(106.681534454096017031613223924991564429656d0 &
         +mx*(268.098486573117433951562111736259672695883d0 &
         +mx*(683.624114850289804796762005964155730439745d0 &
         +mx*(1763.49708521918740723028849567007874329637d0 &
         +mx*(4592.37475383116380899419201719007475759114d0 &
         +mx*(12053.4410190488892782190764838488156555734d0 &
         +mx*(31846.6630207420816960681624497373078887317d0 &
         +mx*(84621.2213590568080177035346867495326879117d0 &
         +mx*(225956.423182907889987641304430180593010940d0 &
         +mx*(605941.517281758859958050194535269219533685d0 &
         +mx*(1.63108259953926832083633544697688841456604d6 &
         )))))))))))))))))
   elseif(m.le.0.8d0) then
      mx=0.25d0-mc
      elb=     0.895902820924731621258525533131864225704d0 &
         +mx*(0.243140003766786661947749288357729051637d0 &
         +mx*(0.273081875594105531575351304277604081620d0 &
         +mx*(0.486280007533573323895498576715458103274d0 &
         +mx*(1.082747437228230914750752674136983406683d0 &
         +mx*(2.743445290986452500459431536349945437824d0 &
         +mx*(7.555817828670234627048618342026400847824d0 &
         +mx*(22.05194082493752427472777448620986154515d0 &
         +mx*(67.15640644740229407624192175802742979626d0 &
         +mx*(211.2722537881770961691291434845898538537d0 &
         +mx*(681.9037843053270682273212958093073895805d0 &
         +mx*(2246.956231592536516768812462150619631201d0 &
         +mx*(7531.483865999711792004783423815426725079d0 &
         +mx*(25608.51260130241579018675054866136922157d0 &
         +mx*(88140.74740089604971425934283371277143256d0 &
         +mx*(306564.4242098446591430938434419151070722d0 &
         +mx*(1.076036077811072193752770590363885180738d6 &
         +mx*(3.807218502573632648224286313875985190526d6 &
         +mx*(1.356638224422139551020110323739879481042d7 &
         ))))))))))))))))))
      eld=     1.26061282657491161418014946566845780315983d0 &
         +mx*(1.54866563808267658056930177790599939977154d0 &
         +mx*(3.55366941187160761540650011660758187283401d0 &
         +mx*(9.90044467610439875577300608183010716301714d0 &
         +mx*(30.3205666174524719862025105898574414438275d0 &
         +mx*(98.1802586588830891484913119780870074464833d0 &
         +mx*(329.771010434557055036273670551546757245808d0 &
         +mx*(1136.65598974289039303581967838947708073239d0 &
         +mx*(3993.83433574622979757935610692842933356144d0 &
         +mx*(14242.7295865552708506232731633468180669284d0 &
         +mx*(51394.7572916887209594591528374806790960057d0 &
         +mx*(187246.702914623152141768788258141932569037d0 &
         +mx*(687653.092375389902708761221294282367947659d0 &
         +mx*(2.54238553565398227033448846432182516906624d6 &
         +mx*(9.45378121934749027243313241962076028066811d6 &
         +mx*(3.53283630179709170835024033154326126569613d7 &
         +mx*(1.32593262383393014923560730485845833322771d8 &
         +mx*(4.99544968184054821463279808395426941549833d8 &
         +mx*(1.88840934729443872364972817525484292678543d9 &
         +mx*(7.16026753447893719179055010636502508063102d9 &
         +mx*(2.72233079469633962247554894093591262281929d10 &
      ))))))))))))))))))))
   elseif(m.le.0.85d0) then
      mx=0.175d0-mc
      elb=     0.915922052601931494319853880201442948834592d0 &
         +mx*(0.294714252429483394379515488141632749820347d0 &
         +mx*(0.435776709264636140422971598963772380161131d0 &
         +mx*(1.067328246493644238508159085364429570207744d0 &
         +mx*(3.327844118563268085074646976514979307993733d0 &
         +mx*(11.90406004445092906188837729711173326621810d0 &
         +mx*(46.47838820224626393512400481776284680677096d0 &
         +mx*(192.7556002578809476962739389101964074608802d0 &
         +mx*(835.3356299261900063712302517586717381557137d0 &
         +mx*(3743.124548343029102644419963712353854902019d0 &
         +mx*(17219.07731004063094108708549153310467326395d0 &
         +mx*(80904.60401669850158353080543152212152282878d0 &
         +mx*(386808.3292751742460123683674607895217760313d0 &
         +mx*(1.876487670110449342170327796786290400635732d6 &
         +mx*(9.216559908641567755240142886998737950775910d6 &
         ))))))))))))))
      eld=     1.402200569110579095046054435635136986038164d0 &
         +mx*(2.322205897861749446534352741005347103992773d0 &
         +mx*(7.462158366466719682730245467372788273333992d0 &
         +mx*(29.43506890797307903104978364254987042421285d0 &
         +mx*(128.1590924337895775262509354898066132182429d0 &
         +mx*(591.0807036911982326384997979640812493154316d0 &
         +mx*(2830.546229607726377048576057043685514661188d0 &
         +mx*(13917.76431889392229954434840686741305556862d0 &
         +mx*(69786.10525163921228258055074102587429394212d0 &
         +mx*(355234.1420341879634781808899208309503519936d0 &
         +mx*(1.830019186413931053503912913904321703777885d6 &
         +mx*(9.519610812032515607466102200648641326190483d6 &
         +mx*(4.992086875574849453986274042758566713803723d7 &
         +mx*(2.635677009826023473846461512029006874800883d8 &
         +mx*(1.399645765120061118824228996253541612110338d9 &
         +mx*(7.469935792837635004663183580452618726280406d9 &
         +mx*(4.004155595835610574316003488168804738481448d10 &
         +mx*(2.154630668144966654449602981243932210695662d11 &
         )))))))))))))))))
   else
      mx=0.125d0-mc
      elb=     0.931906061029524827613331428871579482766771d0 &
         +mx*(0.348448029538453860999386797137074571589376d0 &
         +mx*(0.666809178846938247558793864839434184202736d0 &
         +mx*(2.210769135708128662563678717558631573758222d0 &
         +mx*(9.491765048913406881414290930355300611703187d0 &
         +mx*(47.09304791027740853381457907791343619298913d0 &
         +mx*(255.9200460211233087050940506395442544885608d0 &
         +mx*(1480.029532675805407554800779436693505109703d0 &
         +mx*(8954.040904734313578374783155553041065984547d0 &
         +mx*(56052.48220982686949967604699243627330816542d0 &
         +mx*(360395.7241626000916973524840479780937869149d0 &
         +mx*(2.367539415273216077520928806581689330885103d6 &
         +mx*(1.582994957277684102454906900025484391190264d7 &
         +mx*(1.074158093278511100137056972128875270067228d8 &
         +mx*(7.380585460239595691878086073095523043390649d8 &
         +mx*(5.126022002555101496684687154904781856830296d9 &
         +mx*(3.593534065502416588712409180013118409428367d10 &
         +mx*(2.539881257612812212004146637239987308133582d11 &
         +mx*(1.808180007145359569674767150594344316702507d12 &
         ))))))))))))))))))
      eld=     1.541690112721819084362258323861459983048179d0 &
         +mx*(3.379176214579645449453938918349243359477706d0 &
         +mx*(14.94058385670236671625328259137998668324435d0 &
         +mx*(81.91773929235074880784578753539752529822986d0 &
         +mx*(497.4900546551479866036061853049402721939835d0 &
         +mx*(3205.184010234846235275447901572262470252768d0 &
         +mx*(21457.32237355321925571253220641357074594515d0 &
         +mx*(147557.0156564174712105689758692929775004292d0 &
         +mx*(1.035045290185256525452269053775273002725343d6 &
         +mx*(7.371922334832212125197513363695905834126154d6 &
         +mx*(5.314344395142401141792228169170505958906345d7 &
         +mx*(3.868823475795976312985118115567305767603128d8 &
         +mx*(2.839458401528033778425531336599562337200510d9 &
         +mx*(2.098266122943898941547136470383199468548861d10 &
         +mx*(1.559617754017662417944194874282275405676282d11 &
         +mx*(1.165096220419884791236699872205721392201682d12 &
         +mx*(8.742012983013913804987431275193291316808766d12 &
         +mx*(6.584725462672366918676967847406180155459650d13 &
         +mx*(4.976798737062434393396993620379481464465749d14 &
         +mx*(3.773018634056605404718444239040628892506293d15 &
         +mx*(2.868263194837819660109735981973458220407767d16 &
         ))))))))))))))))))))
   endif
 
   mcold=mc
   elbold=elb
   eldold=eld

end subroutine ellip_bd

!------------------------------------------------------------------------------------------------------------

subroutine equiv_ellipse(igs, a_equiv, b_equiv)
!--purpose: compute equivalent ellipse for contact area given in igs
   implicit none
!--subroutine parameters:
   type(t_eldiv)             :: igs
   real(kind=8), intent(out) :: a_equiv, b_equiv
!--local variables:
   logical, parameter :: use_kp_method = .true.
   integer            :: nadh, nslip, nplast, nexter
   real(kind=8)       :: tot_area, tot_len, tot_wid

   ! Compute area of the actual contact area

   call eldiv_count(igs, nadh, nslip, nplast, nexter)
   tot_area = igs%grid%dxdy * (nadh+nslip+nplast)

   ! Compute the length and width of an equivalent ellipse approximating the actual contact area

   call areas(igs)
   tot_len  = igs%grid%dx * max(1, igs%ixmax - igs%ixmin + 1)
   tot_wid  = igs%grid%dy * max(1, igs%iymax - igs%iymin + 1)

   if (use_kp_method) then
      ! method of Kik-Piotrowski (2008): maintain area A and ellipticity a/b
      ! call write_log(' Equivalent ellipse by Kik-Piotrowski')
      a_equiv  = sqrt( tot_area/pi ) * sqrt( tot_len/tot_wid )
      b_equiv  = sqrt( tot_area/pi ) * sqrt( tot_wid/tot_len )
   else
      ! our old method (v17.1): maintain semi-axis a and area A
      ! call write_log(' Equivalent ellipse by our old method')
      a_equiv  = 0.5d0 * tot_len
      b_equiv  = tot_area / (pi * a_equiv)
   endif

end subroutine equiv_ellipse

!------------------------------------------------------------------------------------------------------------

subroutine simpflex(ic, aa, bb, mater, fstat, kin, idebug)
!--purpose: compute flexibilities for the simplified theory for given semi-axes aa, bb
   implicit none
!--subroutine parameters:
   type(t_ic)                    :: ic
   type(t_material)              :: mater
   type(t_kincns)                :: kin
   real(kind=8),     intent(in)  :: aa, bb, fstat
   integer,          intent(in)  :: idebug
!--local variables:
   real(kind=8)   :: aob, cc, cxx, cyy, cyz, ch_ksi, ch_eta, ch_phi, f3, f4, f_crp, f_wgt, l_crss

   associate(ga => mater%ga, flx => mater%flx, cksi => kin%cksi, ceta => kin%ceta, cphi => kin%cphi)

   ! ellipticity a/b, effective radius c
   ! Note: aa, bb can be zero for an equivalent ellipse with no contact

   if (bb.gt.1d-10) then
      aob = aa / bb
   elseif (aa.gt.1d-10) then
      aob = 1d6
   else
      aob = 1d0
   endif
   cc  = sqrt(aa * bb)

   if (idebug.ge.4) then
      write(bufout,'(3(a,f9.5),a,f9.2)') ' Hertz-data:   aa= ',aa, ', bb= ',bb,', cc= ',cc,', ga=',ga
      call write_log(1, bufout)
   endif

   ! determine appropriate creepage coefficients

   call linrol(mater%nu, aob, cxx, cyy, cyz)
   if (idebug.ge.2) then
      write(bufout,'(3(a,f9.5))') ' creep-coeff:  Cxx=',cxx, ', Cyy=',cyy,', Cyz=',cyz
      call write_log(1, bufout)
   endif

   ! define characteristic sizes ch_ksi, ch_eta, ch_phi

   ch_ksi = fstat * kin%fntrue / (cc**2 * mater%ga * cxx)
   ch_eta = fstat * kin%fntrue / (cc**2 * mater%ga * cyy)
   ch_phi = fstat * kin%fntrue / (cc**3 * mater%ga * cyz)

   if (idebug.ge.2) then
      write(bufout,'(3(a,f10.6))') ' characteristic cksi=', ch_ksi, ', ceta=', ch_eta,', cphi=',ch_phi
      call write_log(1, bufout)
   endif

   ! compute three flexibilities according to Kalker

   flx(1) = 8d0*aa / (3d0*cxx*ga)
   flx(2) = 8d0*aa / (3d0*cyy*ga)
   flx(3) =  pi*aa*sqrt(aob) / (4d0*cyz*ga)

   if (idebug.ge.1) then
      write(bufout,'(3(a,g12.5))') ' flexbilities: L1=', flx(1),', L2=',flx(2), ', L3=',flx(3)
      call write_log(1, bufout)
   endif

   ! compute combined flexibility L_crss

   if (ic%force.le.0) then
      f3 = 8d0 * aa / 3d0
      f4 = pi * aa**2 / 4d0
      l_crss = (cksi**2 + ceta**2 + (cphi*f4/f3)**2) /                                                  &
                           (cksi**2/flx(1) + ceta**2/flx(2) + (cphi*f4/f3)**2/flx(3))

      ! apply blending approach from 3 to 1 flexibilities

      f_crp  = sqrt( (cksi**2 + ceta**2 + cphi**2) / (ch_ksi**2 + ch_eta**2 + ch_phi**2) )
      f_wgt  = min(1d0, max(0d0, (f_crp - 0d0)/(3d0 - 0d0) ))
      if (idebug.ge.1) then
         write(bufout,'(3(a,g12.5))') '      f_crp=',f_crp,', f_wgt=',f_wgt,', L_crss=',l_crss
         call write_log(1, bufout)
      endif

      flx(1) = (1d0 - f_wgt) * flx(1) + f_wgt * l_crss
      flx(2) = (1d0 - f_wgt) * flx(2) + f_wgt * l_crss
      flx(3) = (1d0 - f_wgt) * flx(3) + f_wgt * l_crss

      if (idebug.ge.1) then
         write(bufout,'(3(a,g12.5))') '      blended: L1=', flx(1),', L2=',flx(2), ', L3=',flx(3)
         call write_log(1, bufout)
      endif
   endif

   end associate

end subroutine simpflex

!------------------------------------------------------------------------------------------------------------

subroutine linrol(nu, aob, cxx, cyy, cyz)
!--purpose: determine the Hertzian creepage coefficients cxx,cyy,cyz from Table E.3 of [4].
   implicit none
!--subroutine arguments:
   real(kind=8), intent(in)  :: nu, aob
   real(kind=8), intent(out) :: cxx, cyy, cyz
!--local variables:
   integer      :: icase, irow, inu, icoef
   real(kind=8) :: pi, al, g, fac, c(3,3), cc(3)
!--table E.3 of Kalker coefficients Cij of the linear theory of rolling contact for elliptical contact areas:
   real(kind=8) :: cij(2,10,3,3)
   data ((((cij(icase,irow,inu,icoef), inu=1,3), icoef=1,3), irow=1,10), icase=1,2)          &
       /  2.51, 3.31, 4.85,   2.51, 2.52, 2.53,   .334, .473, .731,                             &
          2.59, 3.37, 4.81,   2.59, 2.63, 2.66,   .483, .603, .809,                             &
          2.68, 3.44, 4.80,   2.68, 2.75, 2.81,   .607, .715, .889,                             &
          2.78, 3.53, 4.82,   2.78, 2.88, 2.98,   .720, .823, .977,                             &
          2.88, 3.62, 4.83,   2.88, 3.01, 3.14,   .827, .929, 1.07,                             &
          2.98, 3.72, 4.91,   2.98, 3.14, 3.31,   .930, 1.03, 1.18,                             &
          3.09, 3.81, 4.97,   3.09, 3.28, 3.48,   1.03, 1.14, 1.29,                             &
          3.19, 3.91, 5.05,   3.19, 3.41, 3.65,   1.13, 1.25, 1.40,                             &
          3.29, 4.01, 5.12,   3.29, 3.54, 3.82,   1.23, 1.36, 1.51,                             &
          3.40, 4.12, 5.20,   3.40, 3.67, 3.98,   1.33, 1.47, 1.63,                             &

          10.7, 11.7, 12.9,   10.7, 12.8, 16.0,   12.2, 14.6, 18.0,                             &
          6.96, 7.78, 8.82,   6.96, 8.14, 9.79,   5.72, 6.63, 7.89,                             &
          5.57, 6.34, 7.34,   5.57, 6.40, 7.51,   3.79, 4.32, 5.01,                             &
          4.84, 5.57, 6.57,   4.84, 5.48, 6.31,   2.88, 3.24, 3.70,                             &
          4.37, 5.10, 6.11,   4.37, 4.90, 5.56,   2.35, 2.62, 2.96,                             &
          4.06, 4.78, 5.80,   4.06, 4.50, 5.04,   2.01, 2.23, 2.50,                             &
          3.82, 4.54, 5.58,   3.82, 4.21, 4.67,   1.76, 1.95, 2.18,                             &
          3.65, 4.36, 5.42,   3.65, 3.99, 4.39,   1.58, 1.75, 1.94,                             &
          3.51, 4.22, 5.30,   3.51, 3.81, 4.16,   1.44, 1.59, 1.77,                             &
          3.40, 4.12, 5.20,   3.40, 3.67, 3.98,   1.33, 1.47, 1.63  /

   pi = 4d0*datan(1d0)

   ! determine g = min(a/b, b/a), remember which case it is

   if (aob.le.1d0) then
      g = aob
      icase = 1
   else
      g = 1d0 / aob
      icase = 2
   end if

   ! determine cxx, cyy, cyz using asymptotes or interpolation

   if (icase.eq.1 .and. g.lt.0.101d0) then

      ! first asymptotes, g \downarrow 0, a<b

      cxx = pi**2 / 4d0 / (1d0-nu)
      cyy = pi**2 / 4d0
      cyz = pi*sqrt(g) / 3d0 / (1d0-nu) * (1d0 + nu * (dlog(16d0 / g) - 5d0))

   else if (icase.eq.2 .and. g.lt.0.101d0) then

      ! last asymptotes, g \downarrow 0, a>b

      al = dlog(16d0 / g**2)
      cxx = 2d0*pi / (al-2d0*nu) / g * (1d0 + (3d0 - log(4d0)) / (al - 2d0*nu))
      cyy = 2d0*pi / g * (1d0 + (1d0-nu) * (3d0-log(4d0)) /                                           &
                            (      (1d0-nu) * al + 2d0*nu)  ) / ((1d0-nu) * al + 2d0*nu) 
      cyz = 2d0*pi / 3d0 / g**1.5d0 / ((1d0-nu) * al - 2d0 + 4d0*nu)

   else

      ! in range of the table, 0.101 < g <= 1.0

      ! take values from table "icase" (1: a<b, 2: a>b)
      ! take irow == floor(10*g-eps), e.g. g = 0.12 --> irow = 1
      ! interpolation weight == rem(10*g) --> fac = 0.2

      irow  = 10d0*g - 0.005d0
      fac = 10d0*g - real(irow)

      ! perform interpolation for all c_i ("icoef") and all nu ("inu")

      do inu = 1, 3
         do icoef = 1, 3
            c(inu,icoef) =    fac  * cij(icase, irow+1, inu, icoef) +                                &
                       (1d0 - fac) * cij(icase, irow  , inu, icoef)
         enddo
      enddo

      ! interpolate to actual nu
      ! formula?  fitting parabola?  inverse interpolation?

      do icoef = 1, 3
         cc(icoef) =   (nu-0.25d0) * (nu-0.50d0) *  8d0 / c(1,icoef)                                 &
                     - (nu-0.50d0) *  nu         * 16d0 / c(2,icoef)                                 &
                     +  nu         * (nu-0.25d0) *  8d0 / c(3,icoef)
      enddo

      cxx = 1d0 / cc(1)
      cyy = 1d0 / cc(2)
      cyz = 1d0 / cc(3)
   end if

end subroutine linrol

!------------------------------------------------------------------------------------------------------------

end module m_hertz
