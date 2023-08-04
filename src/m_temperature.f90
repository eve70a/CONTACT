!------------------------------------------------------------------------------------------------------------
! m_temperature - computation of surface temperatures due to frictional heating
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_temperature

use m_hierarch_data

implicit none
private

public  calc_temp
private calc_temp_steady
! private calc_temp_transient

contains

!------------------------------------------------------------------------------------------------------------

   subroutine calc_temp(ic, mater, kin, cgrid, outpt)
!--purpose: switch between different versions of temperature calculation
      implicit none
!--subroutine-parameters:
      type(t_ic),        intent(in)    :: ic
      type(t_material),  intent(in)    :: mater
      type(t_kincns),    intent(in)    :: kin
      type(t_grid),      intent(in)    :: cgrid
      type(t_output),    intent(inout) :: outpt
!--local variables:

      call timer_start(itimer_temper)

      if (ic%tang.eq.3) then
         call calc_temp_steady(ic, mater, kin, cgrid, outpt)
      else
         call write_log('INTERNAL ERROR: Temperature calculation while T<>3.')
         call abort_run()

         ! call calc_temp_transient(ic, mater, kin, outpt)
      endif

      call timer_stop(itimer_temper)
   end subroutine calc_temp

!------------------------------------------------------------------------------------------------------------

   subroutine calc_temp_steady(ic, mater, kin, cgrid, outpt)
!--purpose: Compute the surface temperatures of wheel and rail in steady rolling scenarios. Using analytic
!           formulas for convection (rail chill) and heat input (friction) according to (Ertz, 2002).
      implicit none
!--subroutine-parameters:
      type(t_ic),       intent(in)          :: ic
      type(t_material), intent(in)          :: mater
      type(t_kincns),   intent(in)          :: kin
      type(t_grid),     intent(in),  target :: cgrid
      type(t_output),   intent(inout)       :: outpt
!--local variables:
      real(kind=8), parameter :: pi = 4d0 * atan(1d0)
      real(kind=8) :: beta1, beta2    ! effective themal penetration coefficient for bodies 1, 2
      real(kind=8) :: eps_b1          ! heat partitioning factor into body 1
      real(kind=8) :: fact, sumflx, ptabs, sabs, dupl, tmean, alen, xlen, hdx
      integer      :: ix, iy, ii, jx, jj, ixsta, ixend, ixinc
      logical      :: use_plast

      associate(mx     => cgrid%nx,     my     => cgrid%ny,    npot   => cgrid%ntot,                    &
                dx     => cgrid%dx,     x      => cgrid%x,     y      => cgrid%y,                       &
                chi    => kin%chi,      veloc  => kin%veloc,   dt     => kin%dt,                        &
                igs    => outpt%igs,    ps     => outpt%ps,    ss     => outpt%ss,                      &
                upls   => outpt%upls,   uplv   => outpt%uplv,  bta_pl => mater%betapl,                  &
                row1st => igs%row1st,   rowlst => igs%rowlst,                                           &
                bktemp => mater%bktemp, t1     => outpt%temp1, t2     => outpt%temp2)

      hdx = dx / 2d0
      use_plast = ic%mater.eq.4 .and. mater%tau_c0.gt.1d-10 .and. mater%tau_c0.le.1d10

      ! compute thermal penetration coefficient for bodies 1 and 2
      ! with lambda in [W/mm-C], dens in [kg/mm^3], heatcp in [J/kg-C], [W]=[N-m/s], [J]=[W-s]
      ! beta is in mm/m * sqrt[ N-m/s-mm-C * kg/mm^3 * N-m/kg-C ] = [N/mm-sqrt(s)-C]

      beta1   = 1d3 * sqrt(mater%lambda(1) * mater%dens(1) * mater%heatcp(1))
      beta2   = 1d3 * sqrt(mater%lambda(2) * mater%dens(2) * mater%heatcp(2))

      ! write(bufout,*) ' beta=', beta1, beta2, ' [1e3 W s05 / mm2 C]'
      ! call write_log(1, bufout)

      ! compute heat partitioning epsilon_w according to Ertz, eq.25, ignoring v0 and vw

      eps_b1    = beta1 / (beta1 + beta2)

      !------------------------------------------------------------------------------------------------------
      ! Part 1: Set temperatures due to convection and conduction
      !------------------------------------------------------------------------------------------------------

      ! Initialize the rail and wheel temperatures at the background / initial values

      call gf3_set(AllElm, bktemp(1), t1, ikZDIR)
      call gf3_set(AllElm, bktemp(2), t2, ikZDIR)

      ! Set temperatures due to convection: constant T_mean in the contact area

      tmean = bktemp(1) + eps_b1 * (bktemp(2) - bktemp(1))
      call gf3_set(AllInt, tmean, t1, ikZDIR)
      call gf3_set(AllInt, tmean, t2, ikZDIR)

      ! make sure that row1st and rowlst are up to date

      call areas(igs)

      ! After leaving the contact area, the temperatures regress to the background value, cf. Ertz eq.(47)

      do iy = 1, my
         alen = max(0, rowlst(iy)-row1st(iy)+1) * hdx              ! half length of contact strip at row iy

         ! set range of elements at trailing edge of contact

         if (abs(chi-pi).lt.1d0) then
            ixsta = max(1,rowlst(iy))  ! int. element at trailing edge, 0 if no contact occurs on row iy
            ixend = mx                 ! last element of row
            ixinc = 1                  ! direction of increasing time
         else
            ixsta = min(mx,row1st(iy)) ! int. element at trailing edge, mx+1 if no contact occurs on row iy
            ixend =  1
            ixinc = -1
         endif

         ! write(*,*) 'iy=',iy,': row',row1st(iy),'-',rowlst(iy),', ix=',ixsta,':',ixinc,':',ixend

         ! loop over exterior elements from trailing edge onwards (excluding interior element ixsta)

         do ix = ixsta+ixinc, ixend, ixinc
            ii = ix + (iy-1)*mx
            xlen = alen + dx * (ix-ixsta-0.5d0) * ixinc      ! distance from leading edge of contact
            ! if (ix.eq.ixsta) write(*,*) 'iy=',iy,': a,x=',alen,xlen,', asin(',sqrt(2d0*alen/xlen),')=', &
            !        asin(sqrt(2d0*alen/xlen))
            !   write(*,*)
            !   write(*,*) 'Hier gaan we t2%vn(ii) opbouwen voor ii = ', ii
            !   write(*,*)
            !   write(*,*) 'bktemp(2)                ', bktemp(2)
            !   write(*,*) 'bktemp(1)                ', bktemp(1)
            !   write(*,*) '(bktemp(2) - bktemp(1))  ', (bktemp(2) - bktemp(1))
            !   write(*,*) '(2d0/pi) * (1d0-eps_b1)  ', (2d0/pi) * (1d0-eps_b1)
            !   write(*,*) '    alen                 ',     alen
            !   write(*,*) '         xlen            ',          xlen
            !   write(*,*) 'sqrt(2d0*alen/xlen)      ', sqrt(2d0*alen/xlen)
            !   write(*,*) 'asin(sqrt(2d0*alen/xlen))', asin(sqrt(2d0*alen/xlen))
            !   write(*,*)
            t1%vn(ii) = bktemp(1) - (bktemp(1) - bktemp(2)) * (2d0/pi) *      eps_b1  *                 &
                           asin(sqrt( 2d0*alen / (alen+xlen) ))
            t2%vn(ii) = bktemp(2) - (bktemp(2) - bktemp(1)) * (2d0/pi) * (1d0-eps_b1) *                 &
                           asin(sqrt( 2d0*alen / (alen+xlen) ))
            !   write(*,*) 'Dus t2%vn(ii) = ', t2%vn(ii)
            !   write(*,*)
         enddo
      enddo

      !------------------------------------------------------------------------------------------------------
      ! Part 2: Add temperature increase due to frictional heat input
      !------------------------------------------------------------------------------------------------------

      ! the temperature increase is obtained by summing the contributions of all previous heat inputs,
      ! cf. Ertz eq.(30).

      if (.false.) then
         write(bufout,*) 'use_plast=', use_plast,', beta_pl=',bta_pl
         call write_log(1, bufout)
      endif

      do iy = 1, my
         do ix = 1, mx
            ii = ix + (iy-1)*mx

            ! sum contributions of all elements jx from the leading edge to point 'ix-1'

            sumflx = 0d0
            dupl   = 0d0

            if (abs(chi-pi).lt.1d0) then

               ! rolling to the left, leading edge at '1', trailing at 'mx', step '1'
               ! note that outpt%ss contains the shift distance [mm]

               do jx = 1, ix-1
                  jj = jx + (iy-1)*mx
                  if (igs%el(jj).eq.Slip .or. igs%el(jj).eq.Plast) then
                     if (use_plast) then
                        dupl = sqrt((upls%vx(jj)-uplv%vx(jj))**2 + (upls%vy(jj)-uplv%vy(jj))**2) / dt
                     endif
                     sabs   = sqrt(ss%vx(jj)**2 + ss%vy(jj)**2) / dt
                     ptabs  = sqrt(ps%vx(jj)**2 + ps%vy(jj)**2)
                     sumflx = sumflx + (sabs + bta_pl*dupl) * ptabs *                                   &
                                                (sqrt(x(ii)-x(jj)+hdx) - sqrt(x(ii)-x(jj)-hdx))
                  endif
               enddo

            else

               ! rolling to the right, leading edge at 'mx', trailing at '1', step '-1'

               do jx = mx, ix+1, -1
                  jj = jx + (iy-1)*mx
                  if (igs%el(jj).eq.Slip .or. igs%el(jj).eq.Plast) then
                     if (use_plast) then
                        dupl = sqrt((upls%vx(jj)-uplv%vx(jj))**2 + (upls%vy(jj)-uplv%vy(jj))**2) / dt
                     endif
                     sabs   = sqrt(ss%vx(jj)**2 + ss%vy(jj)**2) / dt
                     ptabs  = sqrt(ps%vx(jj)**2 + ps%vy(jj)**2)
                     sumflx = sumflx + (sabs + bta_pl*dupl) * ptabs *                                   &
                                                (sqrt(x(jj)-x(ii)+hdx) - sqrt(x(jj)-x(ii)-hdx))
                  endif
               enddo
            endif

            ! add contribution of the point ix itself
            ! With ss in [mm], sabs is in [mm/s], ptabs is in [N/mm^2], hdx in [mm]
            ! Therefore, sumflx is in [ mm/s * N/mm^2 * sqrt(mm) ] = [ N / s*sqrt(mm) ]

            ii = ix + (iy-1)*mx
            if (igs%el(ii).eq.Slip .or. igs%el(ii).eq.Plast) then
               dupl   = sqrt((upls%vx(ii)-uplv%vx(ii))**2 + (upls%vy(ii)-uplv%vy(ii))**2) / dt
               sabs   = sqrt(ss%vx(ii)**2 + ss%vy(ii)**2) / dt
               ptabs  = sqrt(ps%vx(ii)**2 + ps%vy(ii)**2)
               sumflx = sumflx + (sabs + bta_pl*dupl) * ptabs * 2d0 * sqrt(hdx)

               if (.false. .and. ix.eq.5 .and. iy.eq.1) then
                  write(bufout,'(a,i4,5(a,f12.6))') 'ii=',ii,': sx=', ss%vx(ii),', sabs=',sabs,           &
                          ', uplsx=', upls%vx(ii), ', uplvx=',uplv%vx(ii), ', dupl=',dupl
                  call write_log(1, bufout)
               endif
            endif

            ! fact is the multiplication factor for the integrated heat production.
            ! With beta in [N/mm-sqrt(s)-C], veloc in [mm/s], fact is in [sqrt(mm)-s-C / N].
            ! Therefore, t1, t2 are in [sqrt(mm)-s-C / N] * [ N / s*sqrt(mm) ] = [C]

            fact = (2.d0/beta1) / sqrt(pi*veloc)
            t1%vn(ii) = t1%vn(ii) + (1d0-eps_b1) * fact * sumflx
            t2%vn(ii) = t2%vn(ii) +      eps_b1  * fact * sumflx

         enddo  ! ix = 1, mx

         ! BEGIN FORMATTING FOR PUTTING RESULTS IN MATLAB ==================================

         if (ic%x_nmdbg.ge.5) then
            ii = 1 + (iy-1)*mx
            if (abs(y(ii)) < 1d-10) then
               call write_log(' ')
               write(bufout,*) 'Y positie = ', y(ii)
               call write_log(1, bufout)
               call write_log(' ')
               call write_log('A = [ ...')
               alen = max(0, rowlst(iy)-row1st(iy)+1) * hdx
               do ix = 1, mx
                  ii = ix + (iy-1)*mx
                  write(bufout,*) x(ii)/alen, t1%vn(ii), t2%vn(ii)
                  call write_log(1, bufout)
               enddo
               call write_log('];')
               call write_log('figure(6); plot(flipud(A(:,1)),A(:,2),''r*-'');')
               ! call write_log('figure(6); plot(flipud(A(:,1)),A(:,2),''''')
               ! call write_log('hold on; plot(flipud(A(:,1)),A(:,3),''r*-'');')
               call write_log('xlabel(''x/a'');'                       )
               call write_log('ylabel(''Temperature [C]'')')
               call write_log('axis([.6 1.4001 245 275]);')
               call write_log(' ')
            endif
         endif

         ! END FORMATTING FOR PUTTING RESULTS IN MATLAB ====================================

      enddo  ! iy = 1, my

      end associate
   end subroutine calc_temp_steady

!------------------------------------------------------------------------------------------------------------

!  subroutine calc_temp_transient(ic, mater, kin, outpt)
!--purpose: ...
!     implicit none
!--subroutine-parameters:
!     type(t_ic),        intent(in)    :: ic
!     type(t_material),  intent(in)    :: mater
!     type(t_kincns),    intent(in)    :: kin
!     type(t_output),    intent(inout) :: outpt
!--local variables:
!
!  end subroutine calc_temp_transient

!------------------------------------------------------------------------------------------------------------

end module m_temperature
