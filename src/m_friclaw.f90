!------------------------------------------------------------------------------------------------------------
! m_friclaw - data-structures for friction laws and friction variation
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_friclaw

use m_globals
use m_readline
use m_ptrarray
use m_grids
use m_interp

implicit none

public

   public  t_friclaw
   public  fric_resize
   public  fric_init
   public  fric_copy
   public  fric_destroy
   public  fric_print
   public  fric_update
   public  fric_is_veldep
   private fric_get_fstat
   private fric_get_fstat_min
   private fric_get_fstat_avg
   public  fric_veldep_mu_steady
   public  fric_tempdep_mu_temp
   private fric_check_param
   public  fric_interp
   public  fric_input
   public  fric_wrtinp
   public  fric_output

   interface fric_interp
      module procedure fric_interp_arr
      module procedure fric_interp_scalar
   end interface fric_interp

   !---------------------------------------------------------------------------------------------------------

   ! dimensions and parameters

   integer,      parameter :: MAX_NVF =      10     ! #control points for variable friction in user input

   !---------------------------------------------------------------------------------------------------------

   ! data with respect to the friction law used -- see report-m20002
   ! In module 1, the user may use V=1 or 2 and define nvf arbitrarily
   ! In module 3, the program uses either V=0 with nvf=1 or V=2 with nvf=my

   type :: t_friclaw
      integer      :: frclaw_eff
      integer      :: varfrc_eff
      integer      :: nvf
      real(kind=8), dimension(:),   pointer   :: paramvf => NULL()
      real(kind=8), dimension(:),   pointer   :: fstat_arr => NULL()
      real(kind=8), dimension(:),   pointer   :: fkin_arr  => NULL()
      real(kind=8), dimension(:),   pointer   :: flin1   => NULL(), flin2 => NULL()
      real(kind=8), dimension(:),   pointer   :: frat1   => NULL(), frat2 => NULL()
      real(kind=8), dimension(:),   pointer   :: fexp1   => NULL(), fexp2 => NULL()
      real(kind=8), dimension(:),   pointer   :: sabsh1  => NULL(), sabsh2 => NULL()
      real(kind=8), dimension(:),   pointer   :: fref    => NULL(), dfheat => NULL()
      real(kind=8), dimension(:),   pointer   :: tref    => NULL(), dtheat => NULL()
      real(kind=8) :: memdst
      real(kind=8) :: mem_s0
   contains
      procedure :: is_veldep  => fric_is_veldep
      procedure :: mu_steady  => fric_veldep_mu_steady
      procedure :: mu_temp    => fric_tempdep_mu_temp
      procedure :: fstat      => fric_get_fstat
      procedure :: fstat_min  => fric_get_fstat_min
      procedure :: fstat_avg  => fric_get_fstat_avg

      ! frclaw_eff     effective friction law used in current case, see frclaw_inp in t_ic
      ! varfrc_eff     type of friction variation cf varfrc in t_ic
      ! nvf            number of sets of friction parameters
      ! paramvf        independent parameter: ALPHVF for V=1, SVF for V=2
      ! alphvf   [rad] surface inclination alpha for which friction parameters are given
      ! svf      [mm]  track position s which friction parameters are given

      ! fkin_arr  [-]  kinetic coefficients of friction (L=0),
      !                limit values for coefficient of friction at large slip velocities (L=2-4,6).
      ! fstat_arr [-]  static coefficients of friction (L=0-6), derived from other input when L=2-6.
      ! flin1     [-]  size of first linear term in velocity dependent friction law (L=2)
      ! flin2     [-]  size of second linear term in veloc.dep. friction law (L=2)
      ! frat1     [-]  size of first rational term in vel.dep. friction law (L=3)
      ! frat2     [-]  size of second rational term in vel.dep. friction law (L=3)
      ! fexp1     [-]  size of first exponential term in vel.dep. frict.law (L=4)
      ! fexp2     [-]  size of second exponential term in vel.dep. frict.law (L=4)
      ! sabsh1 [mm/s]  absolute slip velocity at which first term is halved (L=2-4)
      ! sabsh2 [mm/s]  absolute slip velocity at which second term is halved (L=2-4)
      ! fref      [-]  coefficient of friction at low surface temperature in case L=6
      ! tref      [C]  lower temperature at which coefficient of friction starts changing (L=6)
      ! dfheat    [-]  change in friction from low to high surface temperatures (L=6)
      ! dtheat    [C]  change temperature after which friction stops changing (L=6)
      ! memdst   [mm]  characteristic distance d_c for the friction memory (L=2-4,6)
      ! mem_s0 [mm/s]  background velocity s_0 for the friction memory (L=2-4,6)

   end type t_friclaw

contains

!------------------------------------------------------------------------------------------------------------

   subroutine fric_resize(fric, nvf)
!--purpose: Resize the arrays in a friclaw data-structure
      implicit none
!--subroutine parameters:
      type(t_friclaw) :: fric
      integer         :: nvf

      fric%nvf  =  nvf
      call reallocate_arr(fric%paramvf, nvf)
      call reallocate_arr(fric%fstat_arr, nvf)
      call reallocate_arr(fric%fkin_arr,  nvf)
      call reallocate_arr(fric%flin1,   nvf)
      call reallocate_arr(fric%flin2,   nvf)
      call reallocate_arr(fric%frat1,   nvf)
      call reallocate_arr(fric%frat2,   nvf)
      call reallocate_arr(fric%fexp1,   nvf)
      call reallocate_arr(fric%fexp2,   nvf)
      call reallocate_arr(fric%sabsh1,  nvf)
      call reallocate_arr(fric%sabsh2,  nvf)
      call reallocate_arr(fric%fref,    nvf)
      call reallocate_arr(fric%dfheat,  nvf)
      call reallocate_arr(fric%tref,    nvf)
      call reallocate_arr(fric%dtheat,  nvf)

   end subroutine fric_resize

!------------------------------------------------------------------------------------------------------------

   subroutine fric_init(fric)
!--purpose: Initialize a friclaw data-structure
      implicit none
!--subroutine parameters:
      type(t_friclaw) :: fric
!--local variables:
      integer         :: nvf

      ! initialize arrays at length 1

      nvf = 1
      call fric_resize(fric, nvf)

      ! store sensible default values

      fric%frclaw_eff =  0
      fric%varfrc_eff =  0

      fric%paramvf   = (/  0.0d0  /)
      fric%fkin_arr  = (/  0.30d0 /)
      fric%fstat_arr = (/  0.30d0 /)
      fric%flin1     = (/  0.0d0  /)
      fric%flin2     = (/  0.0d0  /)
      fric%frat1     = (/  0.0d0  /)
      fric%frat2     = (/  0.0d0  /)
      fric%fexp1     = (/  0.0d0  /)
      fric%fexp2     = (/  0.0d0  /)
      fric%sabsh1    = (/ 10000d0 /)       ! 10 m/s
      fric%sabsh2    = (/ 10000d0 /)       ! 10 m/s
      fric%fref      = (/  0.30d0 /)
      fric%tref      = (/  0.0d0  /)
      fric%dfheat    = (/ -0.10d0 /)
      fric%dtheat    = (/   800d0 /)
      fric%memdst    =     0.001d0         ! 1 micrometer
      fric%mem_s0    =     10d0            ! 10 mm/s

   end subroutine fric_init

!------------------------------------------------------------------------------------------------------------

   subroutine fric_copy(f_in, f_out)
!--purpose: copy friction law data
      implicit none
!--subroutine arguments:
      type(t_friclaw), intent(in)    :: f_in
      type(t_friclaw), intent(inout) :: f_out
!--local variables:
      integer   :: nvf

      ! resize arrays of f_out to same size as used in f_in

      nvf = f_in%nvf
      call fric_resize(f_out, nvf)

      ! copy contents of f_in

      f_out%frclaw_eff     = f_in%frclaw_eff
      f_out%varfrc_eff     = f_in%varfrc_eff

      f_out%paramvf(1:nvf) = f_in%paramvf(1:nvf)
      f_out%fstat_arr(1:nvf) = f_in%fstat_arr(1:nvf)
      f_out%fkin_arr (1:nvf) = f_in%fkin_arr (1:nvf)
      f_out%flin1 (1:nvf)  = f_in%flin1 (1:nvf)
      f_out%flin2 (1:nvf)  = f_in%flin2 (1:nvf)
      f_out%frat1 (1:nvf)  = f_in%frat1 (1:nvf)
      f_out%frat2 (1:nvf)  = f_in%frat2 (1:nvf)
      f_out%fexp1 (1:nvf)  = f_in%fexp1 (1:nvf)
      f_out%fexp2 (1:nvf)  = f_in%fexp2 (1:nvf)
      f_out%sabsh1(1:nvf)  = f_in%sabsh1(1:nvf)
      f_out%sabsh2(1:nvf)  = f_in%sabsh2(1:nvf)
      f_out%fref  (1:nvf)  = f_in%fref  (1:nvf)
      f_out%tref  (1:nvf)  = f_in%tref  (1:nvf)
      f_out%dfheat(1:nvf)  = f_in%dfheat(1:nvf)
      f_out%dtheat(1:nvf)  = f_in%dtheat(1:nvf)
      f_out%memdst         = f_in%memdst  
      f_out%mem_s0         = f_in%mem_s0  

   end subroutine fric_copy

!------------------------------------------------------------------------------------------------------------

   subroutine fric_destroy(fric)
!--purpose: Destroy the arrays in a friclaw data-structure
      implicit none
!--subroutine parameters:
      type(t_friclaw) :: fric

      call destroy_arr(fric%paramvf)
      call destroy_arr(fric%fstat_arr)
      call destroy_arr(fric%fkin_arr)
      call destroy_arr(fric%flin1)
      call destroy_arr(fric%flin2)
      call destroy_arr(fric%frat1)
      call destroy_arr(fric%frat2)
      call destroy_arr(fric%fexp1)
      call destroy_arr(fric%fexp2)
      call destroy_arr(fric%sabsh1)
      call destroy_arr(fric%sabsh2)
      call destroy_arr(fric%fref)
      call destroy_arr(fric%dfheat)
      call destroy_arr(fric%tref)
      call destroy_arr(fric%dtheat)

   end subroutine fric_destroy

!------------------------------------------------------------------------------------------------------------

   subroutine fric_print( f, nam )
!--purpose: Print information on friction data f
      implicit none
!--subroutine arguments:
      type(t_friclaw)  :: f
      character(len=*) :: nam
!--local variables:
      integer          :: ivf
      character(len=5) :: namprm

      namprm = '     '
      if (f%varfrc_eff.eq.1) namprm = 'ALPHA'
      if (f%varfrc_eff.eq.2) namprm = 'SVF  '
      if (f%varfrc_eff.eq.3) namprm = 'ALPHA'

      write(bufout,'(3a,2(i2,a),i4,a)') ' friction law ',trim(nam),' uses L=', f%frclaw_eff,' with V=', &
                f%varfrc_eff, ',', f%nvf,' slices'
      call write_log(1, bufout)

      if (f%frclaw_eff.eq.0) write(bufout, 3200) namprm
      if (f%frclaw_eff.eq.2) write(bufout, 3202) namprm, 'FLIN1', 'FLIN2'
      if (f%frclaw_eff.eq.3) write(bufout, 3202) namprm, 'FRAT1', 'FRAT2'
      if (f%frclaw_eff.eq.4) write(bufout, 3202) namprm, 'FEXP1', 'FEXP2'
      if (f%frclaw_eff.eq.6) write(bufout, 3206) namprm
      call write_log(1, bufout)
 3200 format ( 2x, 3x,a,4x, 3x,'FSTAT',4x, 3x,'FKIN')
 3202 format ( 2x, 3x,a,4x, 3x,'FKIN',5x, 3x,a5,4x, 3x,'SABSH1',3x, 3x,a5,4x, 3x,'SABSH2')
 3206 format ( 2x, 3x,a,4x, 3x,'FREF',5x, 3x,'TREF',5x, 3x,'DFHEAT',3x, 3x,'DTHEAT')

      do ivf = 1, f%nvf
         if (f%frclaw_eff.eq.0) write(bufout, 3300) f%paramvf(ivf), f%fstat_arr(ivf), f%fkin_arr(ivf)
         if (f%frclaw_eff.eq.2) write(bufout, 3300) f%paramvf(ivf), f%fkin_arr(ivf), f%flin1(ivf),      &
                f%sabsh1(ivf), f%flin2(ivf), f%sabsh2(ivf)
         if (f%frclaw_eff.eq.3) write(bufout, 3300) f%paramvf(ivf), f%fkin_arr(ivf), f%frat1(ivf),      &
                f%sabsh1(ivf), f%frat2(ivf), f%sabsh2(ivf)
         if (f%frclaw_eff.eq.4) write(bufout, 3300) f%paramvf(ivf), f%fkin_arr(ivf), f%fexp1(ivf),      &
                f%sabsh1(ivf), f%fexp2(ivf), f%sabsh2(ivf)
         if (f%frclaw_eff.eq.6) write(bufout, 3300) f%paramvf(ivf), f%fref(ivf), f%tref(ivf),           &
                f%dfheat(ivf), f%dtheat(ivf)
         call write_log(1, bufout)
 3300    format ( 2x, 6g12.4)
      enddo

   end subroutine fric_print

!------------------------------------------------------------------------------------------------------------

   subroutine fric_update( f )
!--purpose: Update derived data: static/kinetic coefficients of friction
      implicit none
!--subroutine arguments:
      type(t_friclaw) :: f
!--local variables:
      integer         :: ivf

      do ivf = 1, f%nvf
         if (f%frclaw_eff.eq.0) f%fkin_arr(ivf)  = f%fstat_arr(ivf)

         if (f%frclaw_eff.eq.2) f%fstat_arr(ivf) = f%fkin_arr(ivf) + f%flin1(ivf) + f%flin2(ivf)
         if (f%frclaw_eff.eq.3) f%fstat_arr(ivf) = f%fkin_arr(ivf) + f%frat1(ivf) + f%frat2(ivf)
         if (f%frclaw_eff.eq.4) f%fstat_arr(ivf) = f%fkin_arr(ivf) + f%fexp1(ivf) + f%fexp2(ivf)

         if (f%frclaw_eff.eq.6) f%fstat_arr(ivf) = f%fref(ivf)
         if (f%frclaw_eff.eq.6) f%fkin_arr(ivf)  = f%fref(ivf) + f%dfheat(ivf)
      enddo

   end subroutine fric_update

!------------------------------------------------------------------------------------------------------------

   function fric_is_veldep(f)
   !--function: determine whether the friction law is velocity dependent
      implicit none
   !--result value
      logical           :: fric_is_veldep
   !--subroutine arguments
      class(t_friclaw), intent(in) :: f

      fric_is_veldep  = (f%frclaw_eff.ge.2 .and. f%frclaw_eff.le.4)
   end function fric_is_veldep

!------------------------------------------------------------------------------------------------------------

   function fric_get_fstat(f, iy)
   !--function: obtain overall static coefficient of friction from friclaw data-structure
      implicit none
   !--result value
      real(kind=8)      :: fric_get_fstat
   !--subroutine arguments
      class(t_friclaw),           intent(in) :: f
      integer,          optional, intent(in) :: iy
   !--local variables:
      integer           :: ivf

      ! select slice ivf = 1 or iy

      if (present(iy)) then
         ivf = max(1, min(f%nvf, iy))
      else
         ivf = 1
      endif

      fric_get_fstat = f%fstat_arr(ivf)

   end function fric_get_fstat

!------------------------------------------------------------------------------------------------------------

   function fric_get_fstat_min(f)
   !--function: obtain minimum static coefficient of friction from friclaw data-structure
      implicit none
   !--result value
      real(kind=8)      :: fric_get_fstat_min
   !--subroutine arguments
      class(t_friclaw),           intent(in) :: f
   !--local variables:
      integer           :: ivf
      real(kind=8)      :: fmin

      fmin = f%fstat_arr(1)
      do ivf = 2, f%nvf
         if (f%fstat_arr(ivf).lt.fmin) fmin = f%fstat_arr(ivf)
      enddo

      fric_get_fstat_min = fmin

   end function fric_get_fstat_min

!------------------------------------------------------------------------------------------------------------

   function fric_get_fstat_avg(f)
   !--function: obtain average static coefficient of friction from friclaw data-structure
      implicit none
   !--result value
      real(kind=8)      :: fric_get_fstat_avg
   !--subroutine arguments
      class(t_friclaw),           intent(in) :: f
   !--local variables:
      integer           :: ivf
      real(kind=8)      :: fsum

      fsum = f%fstat_arr(1)
      do ivf = 2, f%nvf
         fsum = fsum + f%fstat_arr(ivf)
      enddo

      fric_get_fstat_avg = fsum / real(f%nvf)

   end function fric_get_fstat_avg

!------------------------------------------------------------------------------------------------------------

   function fric_veldep_mu_steady(f, slpvel, iy)
   !--function: compute the effective coefficient of friction at steady sliding at velocity slpvel
      implicit none
   !--result value
      real(kind=8)      :: fric_veldep_mu_steady
   !--subroutine arguments
      class(t_friclaw),          intent(in) :: f
      real(kind=8),              intent(in) :: slpvel
      integer,         optional, intent(in) :: iy
   !--local variables:
      integer           :: ivf
      real(kind=8)      :: f_eff

      ! select slice ivf = 1 or iy

      if (present(iy)) then
         ivf = max(1, min(f%nvf, iy))
      else
         ivf = 1
      endif

      if (f%frclaw_eff.eq.0) then

         f_eff = f%fstat_arr(ivf)

      elseif (f%frclaw_eff.eq.2) then

         f_eff = f%fkin_arr(ivf)                                                                        &
            + f%flin1(ivf) * max(0d0, 1d0-0.5d0*slpvel/f%sabsh1(ivf))                                   &
            + f%flin2(ivf) * max(0d0, 1d0-0.5d0*slpvel/f%sabsh2(ivf))

      elseif (f%frclaw_eff.eq.3) then

         f_eff = f%fkin_arr(ivf)                                                                        &
            + f%frat1(ivf) / (1d0 +  slpvel/f%sabsh1(ivf))                                              &
            + f%frat2(ivf) / (1d0 + (slpvel/f%sabsh2(ivf))**2)

      elseif (f%frclaw_eff.eq.4) then

         f_eff = f%fkin_arr(ivf)                                                                        &
            + f%fexp1(ivf) * exp(-log(2d0)*slpvel/f%sabsh1(ivf))                                        &
            + f%fexp2(ivf) * exp(-log(2d0)*slpvel/f%sabsh2(ivf))

      elseif (f%frclaw_eff.eq.6) then

         f_eff = f%fref(ivf)

      else
         f_eff = 0d0
      endif

      fric_veldep_mu_steady = f_eff

   end function fric_veldep_mu_steady

!------------------------------------------------------------------------------------------------------------

   function fric_tempdep_mu_temp(f, temp, iy)
   !--function: compute the effective coefficient of friction at current temperature
      implicit none
   !--result value
      real(kind=8)      :: fric_tempdep_mu_temp
   !--subroutine arguments
      class(t_friclaw),          intent(in) :: f
      real(kind=8),              intent(in) :: temp
      integer,         optional, intent(in) :: iy
   !--local variables:
      integer           :: ivf
      real(kind=8)      :: dtemp, f_eff

      ! select slice ivf = 1 or iy

      if (present(iy)) then
         ivf = max(1, min(f%nvf, iy))
      else
         ivf = 1
      endif

      if (f%frclaw_eff.eq.6) then

         dtemp = max(0d0, min(f%dtheat(ivf), temp - f%tref(ivf)))
         f_eff = f%fref(ivf) + dtemp * f%dfheat(ivf) / f%dtheat(ivf)

      else

         f_eff = 0d0

      endif

      fric_tempdep_mu_temp = f_eff

   end function fric_tempdep_mu_temp

!------------------------------------------------------------------------------------------------------------

   subroutine fric_check_param( nvf, paramvf, ierror )
!--purpose: Check param-values of variable friction: must be strictly increasing
   implicit none
!--subroutine arguments:
   integer,         intent(in)  :: nvf
   real(kind=8),    intent(in)  :: paramvf(nvf)
   integer,         intent(out) :: ierror
!--local variables:
   real(kind=8), parameter :: thresh = 1d-6
   integer                 :: ivf
 
   ierror = 0

   ! check that variable friction values are given with param's in increasing order
   ! repeated values are prohibited to avoid division by zero when using interpolation

   ivf  = 1
   do while(ierror.le.0 .and. ivf.lt.nvf)
      ivf = ivf + 1
      if (paramvf(ivf)-paramvf(ivf-1).lt.thresh) ierror = ivf
   enddo

   end subroutine fric_check_param

!------------------------------------------------------------------------------------------------------------

   subroutine fric_interp_arr( f_in, n_out, param_out, f_out )
!--purpose: Interpolate friction law parameters from f_in to f_out
   implicit none
!--subroutine arguments:
   integer,         intent(in)  :: n_out
   real(kind=8),    intent(in)  :: param_out(n_out)
   type(t_friclaw), intent(in)  :: f_in
   type(t_friclaw)              :: f_out
!--local variables:
   integer       :: frclaw, my_error, sub_error

   my_error = 0

   ! resize arrays in f_out to appropriate size

   call fric_resize(f_out, n_out)

   ! copy/interpolate friction parameters & enforce consistency rules

   frclaw           = f_in%frclaw_eff
   f_out%frclaw_eff = frclaw

   ! store parameter positions used for interpolation

   f_out%paramvf(1:n_out) = param_out(1:n_out)

   ! L=0: Coulomb friction

   if (frclaw.eq.0) then
      call interp_1d( f_in%nvf, f_in%paramvf, f_in%fstat_arr, n_out, param_out, f_out%fstat_arr, sub_error )
      if (my_error.eq.0) my_error = sub_error
   endif

   ! L=2-4: Velocity-dependent friction

   if (frclaw.ge.2 .and. frclaw.le.4) then
      call interp_1d( f_in%nvf, f_in%paramvf, f_in%fkin_arr, n_out, param_out, f_out%fkin_arr, sub_error )
      if (my_error.eq.0) my_error = sub_error
      call interp_1d( f_in%nvf, f_in%paramvf, f_in%sabsh1, n_out, param_out, f_out%sabsh1, sub_error )
      if (my_error.eq.0) my_error = sub_error
      call interp_1d( f_in%nvf, f_in%paramvf, f_in%sabsh2, n_out, param_out, f_out%sabsh2, sub_error )
      if (my_error.eq.0) my_error = sub_error
      f_out%memdst = f_in%memdst
      f_out%mem_s0 = f_in%mem_s0
   endif
   if (frclaw.eq.2) then
      call interp_1d( f_in%nvf, f_in%paramvf, f_in%flin1, n_out, param_out, f_out%flin1, sub_error )
      if (my_error.eq.0) my_error = sub_error
      call interp_1d( f_in%nvf, f_in%paramvf, f_in%flin2, n_out, param_out, f_out%flin2, sub_error )
      if (my_error.eq.0) my_error = sub_error
   endif
   if (frclaw.eq.3) then
      call interp_1d( f_in%nvf, f_in%paramvf, f_in%frat1, n_out, param_out, f_out%frat1, sub_error )
      if (my_error.eq.0) my_error = sub_error
      call interp_1d( f_in%nvf, f_in%paramvf, f_in%frat2, n_out, param_out, f_out%frat2, sub_error )
      if (my_error.eq.0) my_error = sub_error
   endif
   if (frclaw.eq.4) then
      call interp_1d( f_in%nvf, f_in%paramvf, f_in%fexp1, n_out, param_out, f_out%fexp1, sub_error )
      if (my_error.eq.0) my_error = sub_error
      call interp_1d( f_in%nvf, f_in%paramvf, f_in%fexp2, n_out, param_out, f_out%fexp2, sub_error )
      if (my_error.eq.0) my_error = sub_error
   endif

   ! L=6: Temperature-dependent friction

   if (frclaw.eq.6) then
      call interp_1d( f_in%nvf, f_in%paramvf, f_in%fref,   n_out, param_out, f_out%fref, sub_error )
      if (my_error.eq.0) my_error = sub_error
      call interp_1d( f_in%nvf, f_in%paramvf, f_in%tref,   n_out, param_out, f_out%tref, sub_error )
      if (my_error.eq.0) my_error = sub_error
      call interp_1d( f_in%nvf, f_in%paramvf, f_in%dfheat, n_out, param_out, f_out%dfheat, sub_error )
      if (my_error.eq.0) my_error = sub_error
      call interp_1d( f_in%nvf, f_in%paramvf, f_in%dtheat, n_out, param_out, f_out%dtheat, sub_error )
      if (my_error.eq.0) my_error = sub_error
      f_out%memdst = f_in%memdst
      f_out%mem_s0 = f_in%mem_s0
   endif

   ! update derived data: static/kinetic coefficients of friction

   call fric_update(f_out)

   end subroutine fric_interp_arr

!------------------------------------------------------------------------------------------------------------

   subroutine fric_interp_scalar( f_in, param_out, f_out )
!--purpose: Interpolate friction law parameters from f_in to f_out
   implicit none
!--subroutine arguments:
   real(kind=8),    intent(in)  :: param_out
   type(t_friclaw), intent(in)  :: f_in
   type(t_friclaw)              :: f_out
!--local variables:
   real(kind=8)                 :: rarr(1)

   rarr(1) = param_out
   call fric_interp_arr( f_in, 1, rarr, f_out )

   end subroutine fric_interp_scalar

!------------------------------------------------------------------------------------------------------------

   subroutine fric_input(linp, ncase, linenr, ic_varfrc, ic_frclaw_inp, my, fric, idebug, ieof,         &
                lstop, zerror)
!--purpose: read input for V- and L-digits: friction and friction variation
      implicit none
!--subroutine parameters:
      integer                  :: linp, ncase, ic_varfrc, ic_frclaw_inp, my, linenr, ieof, idebug
      logical                  :: lstop, zerror
      type(t_friclaw)          :: fric
!--local variables:
      integer, parameter :: mxnval = 20
      integer            :: ints(mxnval), ivf, iofs, nval, ierror
      logical            :: lparam, flags(mxnval)
      real(kind=8)       :: dbles(mxnval)
      character(len=256) :: strngs(mxnval)
      character(len=1)   :: tparam
      character(len=5)   :: namprm
      character(len=6)   :: types

      ! Copy effective friction law used

      fric%frclaw_eff = ic_frclaw_inp
      fric%varfrc_eff = ic_varfrc

      ! Determine NVF, number of sets of parameters used for friction variation

      if (ic_varfrc.eq.0) then                          ! V = 0:  NVF = 1

         fric%nvf = 1

      elseif (ic_varfrc.eq.1 .or. ic_varfrc.eq.2) then  ! V = 1, 2:  NVF from input

         call readline(linp, ncase, linenr, 'number of points for var.friction', 'i', ints, dbles,      &
                        flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)

         fric%nvf = ints(1)
         zerror = zerror .or. .not.check_irng ('NVF', fric%nvf, 1, MAX_NVF)

      elseif (ic_varfrc.eq.3) then                      ! V = 3:  NVF == MY

         fric%nvf = my

      endif

      ! Resize arrays for friction parameters

      call fric_resize(fric, fric%nvf)

      ! Read input depending on V-digit:  ALPHVF values for V = 1, SVF values for V = 2

      if (ic_varfrc.eq.0 .or. ic_varfrc.eq.3) then
         lparam  = .false.
         tparam  = ' '
         iofs    = 0
      elseif (ic_varfrc.eq.1) then
         lparam  = .true.
         tparam  = 'a'
         iofs    = 1
      else
         lparam  = .true.
         tparam  = 'd'
         iofs    = 1
      endif

      ! read each set of parameters in turn

      do ivf = 1, fric%nvf

         if (.not.lparam) fric%paramvf(ivf) = real(ivf)
         
         if (ic_frclaw_inp.eq.0) then           ! L = 0: [ALPHVF/SVF], FSTAT, FKIN

            types = trim(tparam) // 'dd'
            call readline(linp, ncase, linenr, 'coefficients of friction', types,                       &
                          ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)

            if (lparam) fric%paramvf(ivf) = dbles(1)

            fric%fstat_arr(ivf)  = dbles(iofs+1)
            fric%fkin_arr(ivf)   = dbles(iofs+2)

         elseif (ic_frclaw_inp.eq.2) then       ! L = 2: [ALPHVF/SVF], FKIN, FLIN1, SABSH1, FLIN2, SABSH2

            types = trim(tparam) // 'ddddd'
            call readline(linp, ncase, linenr, 'velocity dep. friction, linear formula', types,         &
                          ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)

            if (lparam) fric%paramvf(ivf) = dbles(1)

            fric%fkin_arr(ivf) = dbles(iofs+1)
            fric%flin1(ivf)    = dbles(iofs+2)
            fric%sabsh1(ivf)   = dbles(iofs+3)
            fric%flin2(ivf)    = dbles(iofs+4)
            fric%sabsh2(ivf)   = dbles(iofs+5)

         elseif (ic_frclaw_inp.eq.3) then       ! L = 3: [ALPHVF/SVF], FKIN, FRAT1, SABSH1, FRAT2, SABSH2

            types = trim(tparam) // 'ddddd'
            call readline(linp, ncase, linenr, 'velocity dep. friction, rational formula', types,       &
                          ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)

            if (lparam) fric%paramvf(ivf) = dbles(1)

            fric%fkin_arr(ivf) = dbles(iofs+1)
            fric%frat1(ivf)    = dbles(iofs+2)
            fric%sabsh1(ivf)   = dbles(iofs+3)
            fric%frat2(ivf)    = dbles(iofs+4)
            fric%sabsh2(ivf)   = dbles(iofs+5)

         elseif (ic_frclaw_inp.eq.4) then       ! L = 4: [ALPHVF/SVF], FKIN, FEXP1, SABSH1, FEXP2, SABSH2

            types = trim(tparam) // 'ddddd'
            call readline(linp, ncase, linenr, 'velocity dep. friction, exponential formula', types,    &
                          ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)

            if (lparam) fric%paramvf(ivf) = dbles(1)

            fric%fkin_arr(ivf) = dbles(iofs+1)
            fric%fexp1(ivf)    = dbles(iofs+2)
            fric%sabsh1(ivf)   = dbles(iofs+3)
            fric%fexp2(ivf)    = dbles(iofs+4)
            fric%sabsh2(ivf)   = dbles(iofs+5)

         elseif (ic_frclaw_inp.eq.6) then       ! L = 6: [ALPHVF/SVF], FREF, TREF, DFHEAT, DTHEAT

            types = trim(tparam) // 'dddd'
            call readline(linp, ncase, linenr, 'temperature dependent friction', types,                 &
                          ints, dbles, flags, strngs, mxnval, nval, idebug, ieof, lstop, ierror)

            if (lparam) fric%paramvf(ivf) = dbles(1)

            fric%fref(ivf)     = dbles(iofs+1)
            fric%tref(ivf)     = dbles(iofs+2)
            fric%dfheat(ivf)   = dbles(iofs+3)
            fric%dtheat(ivf)   = dbles(iofs+4)

         endif
      enddo ! nvf

      ! Read MEMDST, same for all positions across the rail profile

      if ((ic_frclaw_inp.ge.2 .and. ic_frclaw_inp.le.4) .or. ic_frclaw_inp.eq.6) then

         call readline(linp, ncase, linenr, 'friction memory', 'dd', ints, dbles, flags, strngs,        &
                       mxnval, nval, idebug, ieof, lstop, ierror)

         fric%memdst = dbles(1)
         fric%mem_s0 = dbles(2)
      endif

      ! Fill derived data: Update static/kinetic coefficients of friction

      call fric_update( fric )

      if (idebug.ge.3) call fric_print( fric, 'fric' )

      ! Check parameters of friction law

      if (ic_varfrc.eq.1 .or. ic_varfrc.eq.2) then
         if (ic_varfrc.eq.1) namprm = 'ALPHA'
         if (ic_varfrc.eq.2) namprm = 'SVF  '

         call fric_check_param( fric%nvf, fric%paramvf, ierror )

         if (ierror.gt.0) then
            zerror = .true.
            write(lout, 2101) namprm, namprm, ierror-1, fric%paramvf(ierror-1),                         &
                                      namprm, ierror, fric%paramvf(ierror)
            write(   *, 2101) namprm, namprm, ierror-1, fric%paramvf(ierror-1),                         &
                                      namprm, ierror, fric%paramvf(ierror)
 2101       format (' Input: ERROR. For variable friction, ',a,' must be strictly increasing.', /,      &
                    '               ',a,'(',i2,')=',f10.6,', ',a,'(',i2,')=',f10.6,'.')
         endif
      endif

      do ivf = 1, fric%nvf

         if (ic_frclaw_inp.eq.0) then

         ! - when L=0, fstat must be equal to fkin, >= 1d-4
         !             Note: fric_update copies fkin := fstat

            zerror = zerror .or. .not.check_range('FKIN', fric%fkin_arr(ivf), 1d-4, 1d10)

            if (fric%fstat_arr(ivf).lt.0.999*fric%fkin_arr(ivf) .or.                                    &
                fric%fstat_arr(ivf).gt.1.001*fric%fkin_arr(ivf)) then
               zerror = .true.
               write(lout, 2111) fric%fstat_arr(ivf), fric%fkin_arr(ivf)
               write(   *, 2111) fric%fstat_arr(ivf), fric%fkin_arr(ivf)
 2111          format (' Input: ERROR. FSTAT and FKIN must be the same:', f8.4,',',f8.4 )
            endif

         elseif (ic_frclaw_inp.ge.2 .and. ic_frclaw_inp.le.4) then

         !  - when L=2,3 or 4, a lower bound of 0.0 applies for fkin

            zerror = zerror .or. .not.check_range ('FKIN', fric%fkin_arr(ivf), 0d0, 1d10)

         elseif (ic_frclaw_inp.eq.6) then

         !  - when L=6, a lower bound of 0.0 applies for fref and fref+dfheat

            zerror = zerror .or. .not.check_range ('FREF', fric%fref(ivf), 0d0, 1d10)
            zerror = zerror .or. .not.check_range ('FREF+DFHEAT', fric%fkin_arr(ivf), 0d0, 1d10)

         endif

         !  - for L=2-4, the coefficients flin, frat, etc must be positive

         if (ic_frclaw_inp.eq.2) then
            zerror = zerror .or. .not.check_range ('FLIN1', fric%flin1(ivf) , 0d0, 1d20)
            zerror = zerror .or. .not.check_range ('FLIN2', fric%flin2(ivf) , 0d0, 1d20)
            if (fric%flin1(ivf).lt.1d-8) fric%sabsh1(ivf) = max(1d-4, fric%sabsh1(ivf))
            if (fric%flin2(ivf).lt.1d-8) fric%sabsh2(ivf) = max(1d-4, fric%sabsh2(ivf))
         endif
         if (ic_frclaw_inp.eq.3) then
            zerror = zerror .or. .not.check_range ('FRAT1', fric%frat1(ivf) , 0d0, 1d20)
            zerror = zerror .or. .not.check_range ('FRAT2', fric%frat2(ivf) , 0d0, 1d20)
            if (fric%frat1(ivf).lt.1d-8) fric%sabsh1(ivf) = max(1d-4, fric%sabsh1(ivf))
            if (fric%frat2(ivf).lt.1d-8) fric%sabsh2(ivf) = max(1d-4, fric%sabsh2(ivf))
         endif
         if (ic_frclaw_inp.eq.4) then
            zerror = zerror .or. .not.check_range ('FEXP1', fric%fexp1(ivf) , 0d0, 1d20)
            zerror = zerror .or. .not.check_range ('FEXP2', fric%fexp2(ivf) , 0d0, 1d20)
            if (fric%fexp1(ivf).lt.1d-8) fric%sabsh1(ivf) = max(1d-4, fric%sabsh1(ivf))
            if (fric%fexp2(ivf).lt.1d-8) fric%sabsh2(ivf) = max(1d-4, fric%sabsh2(ivf))
         endif

         ! - for L=2-4, the coefficients sabsh1,2 must be > 0
         ! - for L=2-4,6, the memory coefficients dc and s0 must be >= 0

         if (ic_frclaw_inp.ge.2 .and. ic_frclaw_inp.le.4) then
            zerror = zerror .or. .not.check_range ('SABSH1', fric%sabsh1(ivf) , 1d-5, 1d20)
            zerror = zerror .or. .not.check_range ('SABSH2', fric%sabsh2(ivf) , 1d-5, 1d20)
         endif
         if ((ic_frclaw_inp.ge.2 .and. ic_frclaw_inp.le.4) .or. ic_frclaw_inp.eq.6) then
            if (ivf.eq.1) then
               zerror = zerror .or. .not.check_range ('MEMDST', fric%memdst, 0d0, 1d20)
               zerror = zerror .or. .not.check_range ('MEM_S0', fric%mem_s0, 0d0, 1d20)
            endif
         endif

      enddo ! nvf

   end subroutine fric_input

!------------------------------------------------------------------------------------------------------------

   subroutine fric_wrtinp(linp, ic_frclaw_inp, fric)
!--purpose: Write friction input for current case to INPUT-file <experim>.inp, unit linp.
      implicit none
!--subroutine parameters:
      integer          :: linp, ic_frclaw_inp
      type(t_friclaw)  :: fric
!--local variables:
      integer      :: ivf
      character(len=5) :: namprm

      ! keep previous settings when L=1 

      if (ic_frclaw_inp.eq.1) return

      ! V = 0: constant coefficient; V = 3: parameter for each row w.o. parameter value

      if (fric%varfrc_eff.eq.0 .or. fric%varfrc_eff.eq.3) then

         do ivf = 1, fric%nvf

            if (ic_frclaw_inp.eq.0) write(linp, 3101) fric%fstat_arr(ivf), fric%fkin_arr(ivf)

            if (ic_frclaw_inp.eq.2) write(linp, 3102) fric%fkin_arr(ivf), fric%flin1(ivf),              &
                fric%sabsh1(ivf), fric%flin2(ivf), fric%sabsh2(ivf), 'FLIN1'
            if (ic_frclaw_inp.eq.3) write(linp, 3102) fric%fkin_arr(ivf), fric%frat1(ivf),              &
                fric%sabsh1(ivf), fric%frat2(ivf), fric%sabsh2(ivf), 'FRAT1'
            if (ic_frclaw_inp.eq.4) write(linp, 3102) fric%fkin_arr(ivf), fric%fexp1(ivf),              &
                fric%sabsh1(ivf), fric%fexp2(ivf), fric%sabsh2(ivf), 'FEXP1'

            if (ic_frclaw_inp.eq.6) write(linp, 3106) fric%fref(ivf), fric%tref(ivf),                   &
                fric%dfheat(ivf), fric%dtheat(ivf)

         enddo ! ivf

 3101    format( 2g12.4, 34x, 'FSTAT, FKIN')
 3102    format( f8.4, 4x, 2(f8.4, g12.4), 6x, 'FKIN,',a5,',SH1,F2,SH2')
 3106    format( 4g12.4, 10x, 'FREF, TREF, DFHEAT, DTHEAT')

      ! V = 1: variable friction with parameter alpha; V = 2: variable friction with parameter sfc

      else ! varfrc = 1 or 2

         if (fric%varfrc_eff.eq.1) namprm = 'ALPHA'
         if (fric%varfrc_eff.eq.2) namprm = 'SVF  '

         write(linp, 3200) fric%nvf
         do ivf = 1, fric%nvf

            if (ic_frclaw_inp.eq.0) write(linp, 3201) fric%paramvf(ivf), fric%fstat_arr(ivf),            &
                fric%fkin_arr(ivf), namprm

            if (ic_frclaw_inp.eq.2) write(linp, 3202) fric%paramvf(ivf), fric%fkin_arr(ivf),             &
                fric%flin1(ivf), fric%sabsh1(ivf), fric%flin2(ivf), fric%sabsh2(ivf), namprm, 'FLIN1'
            if (ic_frclaw_inp.eq.3) write(linp, 3202) fric%paramvf(ivf), fric%fkin_arr(ivf),             &
                fric%frat1(ivf), fric%sabsh1(ivf), fric%frat2(ivf), fric%sabsh2(ivf), namprm, 'FRAT1'
            if (ic_frclaw_inp.eq.4) write(linp, 3202) fric%paramvf(ivf), fric%fkin_arr(ivf),             &
                fric%fexp1(ivf), fric%sabsh1(ivf), fric%fexp2(ivf), fric%sabsh2(ivf), namprm, 'FEXP1'

            if (ic_frclaw_inp.eq.6) write(linp, 3206) fric%paramvf(ivf), fric%fref(ivf), fric%tref(ivf), &
                fric%dfheat(ivf), fric%dtheat(ivf), namprm

         enddo ! ivf

 3200    format(  i8,    50x, 'NVF')
 3201    format( 3g12.4, 22x, a,', FSTAT, FKIN')
 3202    format( 2f8.4, 2(f8.4, g12.4), 2x, a,', FKIN,',a5,',SHLF')
 3206    format( 2f8.4, g12.4, f8.4, g12.4, 10x, a,', FREF, TREF, DF,DTHEAT')

      endif ! v-digit
     
      if ((ic_frclaw_inp.ge.2 .and. ic_frclaw_inp.le.4) .or. ic_frclaw_inp.eq.6)                        &
         write(linp, 3302) fric%memdst, fric%mem_s0

 3302 format( 2g12.4, 34x, 'MEMDST, MEM_S0')

   end subroutine fric_wrtinp

!------------------------------------------------------------------------------------------------------------

   subroutine fric_output(lout, ic_output, fric)
!--purpose: Write friction input for current case to output-file <experim>.out, unit lout.
      implicit none
!--subroutine parameters:
      integer          :: lout, ic_output
      type(t_friclaw)  :: fric
!--local variables:
      integer      :: ivf
      character(len=5) :: namprm

      write(lout, 3001)
 3001 format (/, ' FRICTION LAW PARAMETERS')

      if (fric%varfrc_eff.le.0) then

         if (fric%frclaw_eff.eq.0) write(lout, 3002)
         if (fric%frclaw_eff.eq.2) write(lout, 3003) 'FLIN1', 'FLIN2'
         if (fric%frclaw_eff.eq.3) write(lout, 3003) 'FRAT1', 'FRAT2'
         if (fric%frclaw_eff.eq.4) write(lout, 3003) 'FEXP1', 'FEXP2'
         if (fric%frclaw_eff.eq.6) write(lout, 3006) 
 3002    format (2x, 3x,'FSTAT',4x, 3x,'FKIN')
 3003    format (2x, 3x,'FKIN',5x, 3x,a5,4x, 3x,'SABSH1',3x, 3x,a5,4x, 3x,'SABSH2')
 3006    format (2x, 3x,'FREF',5x, 3x,'TREF',5x, 3x,'DFHEAT',3x, 3x,'DTHEAT')

         ivf = 1
         if (fric%frclaw_eff.eq.0) write(lout, 3100) fric%fstat_arr(ivf), fric%fkin_arr(ivf)

         if (fric%frclaw_eff.eq.2) write(lout, 3100) fric%fkin_arr(ivf), fric%flin1(ivf),               &
             fric%sabsh1(ivf), fric%flin2(ivf), fric%sabsh2(ivf)
         if (fric%frclaw_eff.eq.3) write(lout, 3100) fric%fkin_arr(ivf), fric%frat1(ivf),               &
             fric%sabsh1(ivf), fric%frat2(ivf), fric%sabsh2(ivf)
         if (fric%frclaw_eff.eq.4) write(lout, 3100) fric%fkin_arr(ivf), fric%fexp1(ivf),               &
             fric%sabsh1(ivf), fric%fexp2(ivf), fric%sabsh2(ivf)

         if (fric%frclaw_eff.eq.6) write(lout, 3100) fric%fref(ivf), fric%tref(ivf),                    &
             fric%dfheat(ivf), fric%dtheat(ivf)

 3100    format (2x, 6g12.4)

      elseif (fric%varfrc_eff.eq.1 .or. fric%varfrc_eff.eq.2) then

         if (fric%varfrc_eff.eq.1) namprm = 'ALPHA'
         if (fric%varfrc_eff.eq.2) namprm = 'SVF  '

         if (fric%frclaw_eff.eq.0) write(lout, 4000) namprm
         if (fric%frclaw_eff.eq.2) write(lout, 4002) namprm, 'FLIN1', 'FLIN2'
         if (fric%frclaw_eff.eq.3) write(lout, 4002) namprm, 'FRAT1', 'FRAT2'
         if (fric%frclaw_eff.eq.4) write(lout, 4002) namprm, 'FEXP1', 'FEXP2'
         if (fric%frclaw_eff.eq.6) write(lout, 4006) namprm
 4000    format ( 2x, 3x,a,4x, 3x,'FSTAT',4x, 3x,'FKIN')
 4002    format ( 2x, 3x,a,4x, 3x,'FKIN',5x, 3x,a5,4x, 3x,'SABSH1',3x, 3x,a5,4x, 3x,'SABSH2')
 4006    format ( 2x, 3x,a,4x, 3x,'FREF',5x, 3x,'TREF',5x, 3x,'DFHEAT',3x, 3x,'DTHEAT')

         do ivf = 1, fric%nvf
            if (fric%frclaw_eff.eq.0) write(lout, 3100) fric%paramvf(ivf), fric%fstat_arr(ivf),          &
                fric%fkin_arr(ivf)
            if (fric%frclaw_eff.eq.2) write(lout, 3100) fric%paramvf(ivf), fric%fkin_arr(ivf),           &
                fric%flin1(ivf), fric%sabsh1(ivf), fric%flin2(ivf), fric%sabsh2(ivf)
            if (fric%frclaw_eff.eq.3) write(lout, 3100) fric%paramvf(ivf), fric%fkin_arr(ivf),           &
                fric%frat1(ivf), fric%sabsh1(ivf), fric%frat2(ivf), fric%sabsh2(ivf)
            if (fric%frclaw_eff.eq.4) write(lout, 3100) fric%paramvf(ivf), fric%fkin_arr(ivf),           &
                fric%fexp1(ivf), fric%sabsh1(ivf), fric%fexp2(ivf), fric%sabsh2(ivf)
            if (fric%frclaw_eff.eq.6) write(lout, 3100) fric%paramvf(ivf), fric%fref(ivf),               &
                fric%tref(ivf), fric%dfheat(ivf), fric%dtheat(ivf)
         enddo

      elseif (fric%varfrc_eff.ge.3) then

         if (fric%frclaw_eff.eq.0) write(lout, 5002)
         if (fric%frclaw_eff.eq.2) write(lout, 5003) 'FLIN1', 'FLIN2'
         if (fric%frclaw_eff.eq.3) write(lout, 5003) 'FRAT1', 'FRAT2'
         if (fric%frclaw_eff.eq.4) write(lout, 5003) 'FEXP1', 'FEXP2'
         if (fric%frclaw_eff.eq.6) write(lout, 5006) 

 5002    format (2x, 4x,'IY',6x, 3x,'FSTAT',4x, 3x,'FKIN')
 5003    format (2x, 4x,'IY',6x, 3x,'FKIN',5x, 3x,a5,4x, 3x,'SABSH1',3x, 3x,a5,4x, 3x,'SABSH2')
 5006    format (2x, 4x,'IY',6x, 3x,'FREF',5x, 3x,'TREF',5x, 3x,'DFHEAT',3x, 3x,'DTHEAT')

         do ivf = 1, fric%nvf
            if (ic_output.ge.4 .or. ivf.le.2 .or. ivf.ge.fric%nvf-1) then
               if (fric%frclaw_eff.eq.0) write(lout, 5100) ivf, fric%fstat_arr(ivf), fric%fkin_arr(ivf)

               if (fric%frclaw_eff.eq.2) write(lout, 5100) ivf, fric%fkin_arr(ivf), fric%flin1(ivf),    &
                   fric%sabsh1(ivf), fric%flin2(ivf), fric%sabsh2(ivf)
               if (fric%frclaw_eff.eq.3) write(lout, 5100) ivf, fric%fkin_arr(ivf), fric%frat1(ivf),    &
                   fric%sabsh1(ivf), fric%frat2(ivf), fric%sabsh2(ivf)
               if (fric%frclaw_eff.eq.4) write(lout, 5100) ivf, fric%fkin_arr(ivf), fric%fexp1(ivf),    &
                   fric%sabsh1(ivf), fric%fexp2(ivf), fric%sabsh2(ivf)

               if (fric%frclaw_eff.eq.6) write(lout, 5100) ivf, fric%fref(ivf), fric%tref(ivf),         &
                   fric%dfheat(ivf), fric%dtheat(ivf)
            elseif (ivf.eq.3) then
               write(lout, 5150)
            endif
         enddo

 5100    format (2x, i8,4x, 6g12.4)
 5150    format (2x,   12x, 8x, '...')

      endif ! V-digit

   end subroutine fric_output

!------------------------------------------------------------------------------------------------------------

end module m_friclaw
