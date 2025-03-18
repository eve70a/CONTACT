!------------------------------------------------------------------------------------------------------------
! m_wr_input - read input-file for one case with w/r contact (module 1)
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------

module m_wr_input

use m_hierarch_data
use m_sinput
use m_wrprof_data
use m_subsurf       , only : subsurf_input, subsurf_wrtinp
implicit none
private

   ! module numbers used for different purposes

   integer, parameter :: modul_spck       = 11

   public  wr_input
   public  wr_input_modul1
   public  wr_input_spck
   private warn_ic_changed
   public  discret_input
   public  trackdata_input
   public  wheelset_input
   public  wr_write_inp

contains

!------------------------------------------------------------------------------------------------------------

   subroutine wr_input (inpdir, lunit, inp, ncase, linenr, wtd, lstop, ierror)
!--purpose: Input-routine for a W/R contact case. 
      implicit none
!--subroutine arguments:
      character*(*), intent(in) :: inpdir
      integer                   :: lunit, inp, linenr, ncase, ierror
      logical                   :: lstop
      type(t_ws_track)          :: wtd
!--local variables:
      integer,      parameter :: mxnval = 20, modul = 1
      integer            :: line0, ldebug, ieof
      logical            :: zerror
      character(len=16)  :: namside

      associate( ic    => wtd%ic,    meta  => wtd%meta, mater => wtd%mater, discr => wtd%discr,         &
                 kin   => wtd%kin,   fric  => wtd%fric, solv  => wtd%solv,  subs  => wtd%subs,          &
                 ws    => wtd%ws,    trk   => wtd%trk,  spck  => wtd%spck,  my_wheel => wtd%ws%whl,     &
                 my_rail => wtd%trk%rai)

      ierror = 0
      ldebug = 1
      ieof   = -1 ! eof=error

      if (ldebug.ge.2) call write_log('--- Start subroutine wr_input ---')
      line0 = linenr

      if (inp.le.1) then
         call write_log('ERROR: screen input not supported in MODULE 1')
         if (lstop) then
             call abort_run()
         else
            ierror = CNTC_err_input
            return
         endif
      endif

      !------------------------------------------------------------------------------------------------------
      ! read & check control integers
      !------------------------------------------------------------------------------------------------------

      zerror = .false.
      call ic_input(lunit, ncase, linenr, 1, ic, subs, ldebug, ieof, lstop, zerror)

      ! Copy the effective D-digit

      if (ic%discns1_inp.ge.2) ic%discns1_eff = ic%discns1_inp

      ! Copy the B- and M-digits to the material parameter data

      mater%bound_eff = ic%bound
      mater%mater_eff = ic%mater

      ! Copy the effective C3-digit to the material parameter data

      if (ic%gencr_inp.ge.2) mater%gencr_eff = ic%gencr_inp

      if (ic%config.eq.0 .or. ic%config.eq.4) then
         namside   = 'left'
      else
         namside   = 'right'
      endif

      if (zerror) then
         if (lstop) then
            call write_log(' Errors found. Aborting.')
            call abort_run()
         else
            call write_log(' Errors found. Returning.')
            ierror = CNTC_err_input
            return
         endif
      endif

      ! set level of debug-output of input-routines. 0 = errors, 1 = warnings/info, >=2 = flow/debug

      if (ic%x_readln.le.0) then
         ldebug = ic%ilvout
      else
         ldebug = ic%x_readln
      endif

      !------------------------------------------------------------------------------------------------------
      ! read & check input for G-digit: iteration parameters
      !------------------------------------------------------------------------------------------------------

      if (ic%gausei_inp.ne.1) then
         call solv_input(lunit, ncase, linenr, ic%gausei_inp, solv, discr%npot_max, ldebug, ieof,       &
                         lstop, zerror)
      endif

      !------------------------------------------------------------------------------------------------------
      ! Read input for V- and L-digits
      !------------------------------------------------------------------------------------------------------

      if (ic%frclaw_inp.ne.1) then
         call fric_input(lunit, ncase, linenr, ic%varfrc, ic%frclaw_inp, 1, fric, ldebug, ieof,         &
                lstop, zerror)
      endif

      ! Adapt the value of the friction coefficient used for scaling of tangential forces

      kin%use_muscal = (ic%varfrc.eq.0 .or. ic%varfrc.eq.2)
      if (kin%use_muscal) then
         kin%muscal = fric%fstat()
      else
         kin%muscal = 1d0
      endif

      !------------------------------------------------------------------------------------------------------
      ! Read input for C3 and M-digits
      !------------------------------------------------------------------------------------------------------

      ! Read elastic material constants when C3=2, 3 or 4

      if (ic%gencr_inp.eq.2 .or. ic%gencr_inp.eq.3 .or. ic%gencr_inp.eq.4) then
         call mater_input(lunit, ncase, linenr, 1, ic, kin, mater, solv, ldebug, ieof, lstop, zerror)
      endif

      !------------------------------------------------------------------------------------------------------
      ! Read input for M2-digit: material parameters for damping
      !------------------------------------------------------------------------------------------------------

      call mater2_input(lunit, ncase, linenr, ic, mater, ldebug, ieof, lstop, zerror)

      !------------------------------------------------------------------------------------------------------
      ! Read input for H-digit: material parameters for temperature calculation
      !------------------------------------------------------------------------------------------------------

      if (ic%heat.eq.3) then
         call heat_input(lunit, ncase, linenr, ic, mater, ldebug, ieof, lstop, zerror)
      endif

      !------------------------------------------------------------------------------------------------------
      ! Read input for D-digit: potential contact, discretisation parameters
      !------------------------------------------------------------------------------------------------------

      if (ic%discns1_inp.ge.2 .and. ic%discns1_inp.le.9) then
         call discret_input(lunit, ncase, linenr, wtd, ic, discr, ldebug, ieof, lstop, zerror)
      endif

      !------------------------------------------------------------------------------------------------------
      ! Read input for Z1-digit: track/roller rig design dimensions, profile(s), deviation
      !------------------------------------------------------------------------------------------------------

      call trackdata_input(inpdir, lunit, ncase, linenr, modul, ic, trk, spck, ldebug, ieof, lstop, zerror)

      if (my_rail%prr%fname.eq.' ' .and. (ic%ztrack.eq.0 .or. ic%ztrack.eq.2)) then
         zerror = .true.
         write(lout, 1901) trim(namside), ic%ztrack
         write(   *, 1901) trim(namside), ic%ztrack
 1901    format (' Input: ERROR. No profile has been given for the ', a, ' rail, Z1 =',i2,'.')
      endif

      !------------------------------------------------------------------------------------------------------
      ! Read input for E-digit: wheel-set dimensions, profile, position and velocity
      !------------------------------------------------------------------------------------------------------

      call wheelset_input(inpdir, lunit, ncase, linenr, modul, ic, ws, trk, spck, ldebug, ieof, lstop,  &
                zerror)

      if (my_wheel%prw%fname.eq.' ' .and.                                                               &
          (ic%ewheel.eq.1 .or. ic%ewheel.eq.2 .or. ic%ewheel.eq.4)) then
         zerror = .true.
         write(lout, 1902) trim(namside), ic%ewheel
         write(   *, 1902) trim(namside), ic%ewheel
 1902    format (' Input: ERROR. No profile has been given for the ', a, ' wheel, E1 =',i2,'.')
      endif

      !------------------------------------------------------------------------------------------------------
      ! read subsurface points
      !------------------------------------------------------------------------------------------------------

      if (ic%stress.ge.2) then
         call subsurf_input(lunit, ic, ncase, linenr, subs, ldebug, zerror)
      endif

      ! abort on errors

      if (zerror) then
         if (lstop) then
            call write_log(' Errors found. Aborting.')
            call abort_run()
         else
            call write_log(' Errors found. Returning.')
            ierror = CNTC_err_input
            return
         endif
      endif

      if (ldebug.ge.3) call write_log('--- end subroutine wr_input ---')
      end associate

   end subroutine wr_input

!------------------------------------------------------------------------------------------------------------

   subroutine wr_input_modul1 (fname, wtd, ierror)
!--purpose: Input-routine for initialization of a W/R configuration in the CONTACT library
      implicit none
!--subroutine arguments:
      character(len=*)          :: fname
      type(t_ws_track)          :: wtd
      integer,      intent(out) :: ierror
!--local variables:
      integer,      parameter :: mxnval = 20
      logical,      parameter :: lstop  = .false.
      integer          :: ints(mxnval), linp, ldebug, inp, ncase, linenr, ix, nval, ieof, modul
      logical          :: flags(mxnval)
      real(kind=8)     :: dbles(mxnval)
      character*256    :: strngs(mxnval), inpdir, fulnam

      ierror = 0

      ldebug = 0
      if (ldebug.ge.2) call write_log('--- Start subroutine wr_input_modul1 ---')

      ! determine full path-name, pre-pending wrkdir when necessary

      call make_absolute_path(fname, wtd%meta%wrkdir, fulnam)

      linp  = get_lunit_tmp_use()
      open(linp, file=fulnam, status='old', err=985)

      ! determine the folder name from the input filename

      inpdir = ' '
      ix = index_pathsep(fulnam, back=.true.)
      if (ix.gt.0) then
         inpdir = fulnam(1:ix-1) // path_sep
      endif

      ! Get the module-number for first case

      ldebug = 1
      ieof   = -1 ! eof=error
      call readLine(linp, ncase, linenr, 'module number', 'i', ints, dbles, flags, strngs, mxnval,      &
                    nval, ldebug, ieof, lstop, ierror)
      modul = 0
      if (ierror.eq.0) modul = ints(1)

      if (modul.ne.1) then

         write(bufout,7200) modul, linenr
         call write_log(2, bufout)
 7200    format (/' ERROR: invalid module number (',i6,') at line',i6,', skip reading inp-file')
         ierror = -1
      
      else

         inp    =  2
         ncase  =  0
         linenr =  0
         call wr_input (inpdir, linp, inp, ncase, linenr, wtd, .false., ierror)

      endif

      close(linp)
      call free_lunit_tmp_use(linp)

      return

 985  continue
         ierror = -2
         write(bufout,'(3a)') ' ERROR: cannot open input-file: "', trim(fulnam),'"'
         call write_log(1, bufout)
         call free_lunit_tmp_use(linp)
         return

   end subroutine wr_input_modul1

!------------------------------------------------------------------------------------------------------------

   subroutine wr_input_spck (fname, wtd, ierror)
!--purpose: Input-routine for a rail-wheel pair in the user-subroutine to Simpack
      implicit none
!--subroutine arguments:
      character(len=*)          :: fname
      type(t_ws_track)          :: wtd
      integer,      intent(out) :: ierror
!--local variables:
      integer,      parameter :: mxnval = 20
      logical,      parameter :: lstop  = .false.
      integer          :: ints(mxnval), lspck, nval, ldebug, ieof, my_ierror, ncase, linenr, modul, ix
      logical          :: flags(mxnval), zerror
      real(kind=8)     :: dbles(mxnval)
      character*256    :: strngs(mxnval), inpdir, fulnam
      type(t_ic)       :: ic0

      ierror = 0
      zerror = .false.
      ldebug = 0
      if (ldebug.ge.2) call write_log('--- Start subroutine wr_input_spck ---')

      associate( ic    => wtd%ic,    meta  => wtd%meta, mater => wtd%mater, discr => wtd%discr,         &
                 kin   => wtd%kin,   fric  => wtd%fric, solv  => wtd%solv,  subs  => wtd%subs,          &
                 ws    => wtd%ws,    trk   => wtd%trk,  spck  => wtd%spck,  my_wheel => wtd%ws%whl,     &
                 my_rail => wtd%trk%rai)

      ieof   = -1 ! eof is considered an error
      linenr =  0
      ncase  =  0

      if (.false.) wtd%meta%ncase = 0

      ! determine full path-name, pre-pending wrkdir when necessary

      call make_absolute_path(fname, meta%wrkdir, fulnam)

      lspck = get_lunit_tmp_use()
      open(lspck, file=fulnam, status='old', err=985)

      ! determine the folder name from the input filename

      inpdir = ' '
      ix = index_pathsep(fulnam, back=.true.)
      if (ix.gt.0) then
         inpdir = fulnam(1:ix-1) // path_sep
      endif

      ! Get the module-number for first case

      ldebug = 1
      ieof   = -1 ! eof=error
      call readLine(lspck, ncase, linenr, 'module number', 'i', ints, dbles, flags, strngs, mxnval,      &
                    nval, ldebug, ieof, lstop, my_ierror)
      zerror = zerror .or. (my_ierror.ne.0)

      modul = ints(1)
      if (modul.ne.0 .and. modul.ne.modul_spck) then
         write(bufout,7200) modul, linenr
         call write_log(2, bufout)
 7200    format (/' ERROR: invalid module number (',i6,') at line',i6,', skip reading inp-file')
         ierror = -1
      endif

      if (ierror.eq.0) then

         ! read & check control integers

         ic0 = ic
         call ic_input(lspck, ncase, linenr, 1, ic, subs, ldebug, ieof, lstop, zerror)

         ! revert control digits that cannot be changed through the spck.inp-file

         if (warn_ic_changed('C1 (CONFIG)', ic0%config, ic%config, .false.)) ic%config = ic0%config
         if (warn_ic_changed('P (PVTIME)',  ic0%pvtime, ic%pvtime, .true.))  ic%pvtime = ic0%pvtime
         if (warn_ic_changed('N (NORM)',    ic0%norm,   ic%norm,   .true.))  ic%norm   = ic0%norm
         if (warn_ic_changed('F1 (FORCE)',  ic0%force1, ic%force1, .true.))  ic%force1 = ic0%force1
         if (warn_ic_changed('H (HEAT)',    ic0%heat,   ic%heat,   .true.))  ic%heat   = ic0%heat
         if (warn_ic_changed('S (STRESS)',  ic0%stress, ic%stress, .true.))  ic%stress = ic0%stress
         if (warn_ic_changed('I (IESTIM)',  ic0%iestim, ic%iestim, .false.)) ic%iestim = ic0%iestim

         ! reject features not supported in Simpack input

         zerror = zerror .or. .not.check_irng ('Control digit T',  ic%tang,    0, 3)
         zerror = zerror .or. .not.check_irng ('Control digit S',  ic%stress,  0, 0)
         zerror = zerror .or. .not.check_2rng ('Control digit Z1', ic%ztrack,  0, 0, 3, 3)
         zerror = zerror .or. .not.check_2rng ('Control digit E1', ic%ewheel,  0, 0, 3, 3)

         ! Copy the effective D-digit

         if (ic%discns1_inp.ge.2) ic%discns1_eff = ic%discns1_inp

         ! Copy the B- and M-digits to the material parameter data

         mater%bound_eff = ic%bound
         mater%mater_eff = ic%mater

         ! Copy the effective C3-digit to the material parameter data

         if (ic%gencr_inp.ge.2) mater%gencr_eff = ic%gencr_inp

         if (zerror) then
            call write_log(' Errors found. Returning.')
            ierror = -1
         endif
      endif ! ierror==0

      if (ierror.eq.0) then

         ! read & check input for G-digit: iteration parameters

         if (ic%gausei_inp.ne.1) then
            call solv_input(lspck, ncase, linenr, ic%gausei_inp, solv, discr%npot_max, ldebug, ieof,    &
                            lstop, zerror)
         endif

         ! Read input for V- and L-digits

         if (ic%frclaw_inp.ne.1) then
            call fric_input(lspck, ncase, linenr, ic%varfrc, ic%frclaw_inp, 1, fric, ldebug, ieof,      &
                   lstop, zerror)
         endif

         ! Adapt the value of the friction coefficient used for scaling of tangential forces

         kin%use_muscal = (ic%varfrc.eq.0 .or. ic%varfrc.eq.2)
         if (kin%use_muscal) then
            kin%muscal = fric%fstat()
         else
            kin%muscal = 1d0
         endif

         ! Read input for C3 and M-digits when C3=2, 3 or 4

         if (ic%gencr_inp.eq.2 .or. ic%gencr_inp.eq.3 .or. ic%gencr_inp.eq.4) then
            call mater_input(lspck, ncase, linenr, 11, ic, kin, mater, solv, ldebug, ieof, lstop, zerror)
         endif

         ! Read input for M2-digit: material parameters for damping

         call mater2_input(lspck, ncase, linenr, ic, mater, ldebug, ieof, lstop, zerror)

         ! Read input for H-digit: material parameters for temperature calculation

         if (ic%heat.eq.3) then
            call heat_input(lspck, ncase, linenr, ic, mater, ldebug, ieof, lstop, zerror)
         endif

         ! Read input for D-digit: potential contact, discretisation parameters

         if (ic%discns1_inp.ge.2 .and. ic%discns1_inp.le.9) then
            call discret_input(lspck, ncase, linenr, wtd, ic, discr, ldebug, ieof, lstop, zerror)
         endif

         ! Read input for Z1-digit: track/roller rig design dimensions, profile(s), deviation

         call trackdata_input(inpdir, lspck, ncase, linenr, modul, ic, trk, spck, ldebug, ieof, lstop,  &
                zerror)

         ! Read input for E-digit: wheel-set dimensions, profile, position and velocity

         call wheelset_input(inpdir, lspck, ncase, linenr, modul, ic, ws, trk, spck, ldebug, ieof,      &
                lstop, zerror)

         ! read subsurface points

         if (ic%stress.ge.2) then
            call subsurf_input(lspck, ic, ncase, linenr, subs, ldebug, zerror)
         endif

         if (zerror) then
            call write_log(' Errors found. Returning.')
            ierror = -1
         endif

      endif ! ierror==0

      close(lspck)
      call free_lunit_tmp_use(lspck)

      return

 985  continue
         ierror = -2
         write(bufout,'(3a)') ' ERROR: cannot open input-file: "', trim(fulnam),'"'
         call write_log(1, bufout)
         call free_lunit_tmp_use(lspck)
         return

      end associate
   end subroutine wr_input_spck

!------------------------------------------------------------------------------------------------------------

   logical function warn_ic_changed (descrp, ic1, ic2, lprint)
!--purpose: perform check on ic-values ic1 == ic2
      implicit none
      character(len=*)  :: descrp
      integer           :: ic1, ic2
      logical, optional :: lprint
      logical           :: zwarn, my_lprint

      my_lprint = .false.
      if (present(lprint)) my_lprint = lprint

      zwarn = .false.
      if (ic1.ne.ic2) then
         zwarn = .true.
         if (my_lprint) then
            write(bufout, 1000) descrp, ic1
            call write_log(1, bufout)
 1000       format (' Warning: Control digit ',a,' =',i2,' cannot be changed through the Simpack inp-file.')
         endif
      endif
      warn_ic_changed = zwarn
   end function warn_ic_changed

!------------------------------------------------------------------------------------------------------------

   subroutine discret_input(lunit, ncase, linenr, wtd, ic, discr, ldebug, ieof, lstop, zerror)
!--purpose: read input for D1-digit: discretisation parameters
      implicit none
!--subroutine arguments:
      integer, intent(in)      :: lunit, ncase, ldebug
      integer, intent(inout)   :: linenr, ieof
      type(t_ws_track)         :: wtd
      type(t_ic)               :: ic
      type(t_discret)          :: discr
      logical, intent(in)      :: lstop
      logical, intent(inout)   :: zerror
!--local variables:
      integer, parameter :: mxnval = 10
      integer            :: ints(mxnval), nval, ierror
      logical            :: flags(mxnval)
      real(kind=8)       :: dbles(mxnval), dx_prv, ds_prv
      character(len=256) :: strngs(mxnval)

      call readline(lunit, ncase, linenr, 'discretisation parameters', 'dddaddDD', ints, dbles,         &
                     flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
      zerror = zerror .or. (ierror.ne.0)
      discr%dx        = dbles(1)
      discr%ds        = dbles(2)
      discr%dqrel     = dbles(3)
      discr%angl_sep  = max(0d0, dbles(4))
      discr%dist_sep  = max(0d0, dbles(5))
      discr%dist_comb = dbles(6)
      if (nval.ge.7) then
         discr%dist_turn = dbles(7)
      else
         discr%dist_turn = 2d0 * discr%dist_sep - 1d0 * discr%dist_comb
      endif
      if (nval.ge.8) discr%gap_miss  = dbles(8)

      zerror = zerror .or. .not.check_range ('DX', discr%dx, 1d-8, 1d20)
      zerror = zerror .or. .not.check_range ('DS', discr%ds, 1d-8, 1d20)
      zerror = zerror .or. .not.check_range ('DQREL', discr%dqrel, 5d-7, 10d0)

      ! in transient cases, check element sizes of new/previous case

      if (ncase.gt.1 .and. wtd%numcps.ge.1 .and. ic%use_supergrid()) then
         associate(cp => wtd%allcps(1)%cp)
         dx_prv = cp%gd%cgrid_cur%dx
         ds_prv = cp%gd%cgrid_cur%dy
         if (abs(dx_prv-discr%dx).ge.1d-4*min(dx_prv, discr%dx) .or.                                    &
             abs(ds_prv-discr%ds).ge.1d-4*min(ds_prv, discr%ds)) then
            zerror = .true.
            write(lout,5001) ic%tang, dx_prv, ds_prv
            write(   *,5001) ic%tang, dx_prv, ds_prv
         endif
         end associate
 5001    format(' Input: ERROR. Transient contact (T=',i1,') needs constant DX, DY (',f6.3,',',f6.3,    &
                ') in subsequent cases')
      endif

      if (discr%dist_sep.lt.discr%dist_comb) then
         zerror = .true.
         write(lout, 5011) discr%dist_sep, discr%dist_comb
         write(   *, 5011) discr%dist_sep, discr%dist_comb
      endif
 5011 format (' Input: ERROR. D_SEP =',f8.2,' MUST BE >= D_COMB =',f8.2,'.')

!     if (discr%dist_sep.gt.2d0*discr%dist_comb) then
!        zerror = .true.
!        write(lout, 5012) discr%dist_sep, 2d0*discr%dist_comb
!        write(   *, 5012) discr%dist_sep, 2d0*discr%dist_comb
!     endif
!5012 format (' Input: ERROR. D_SEP =',f8.2,' MUST BE <= 2*D_COMB =',f8.2,'.')

   end subroutine discret_input

!------------------------------------------------------------------------------------------------------------

   subroutine trackdata_input(inpdir, lunit, ncase, linenr, modul, ic, trk, spck, ldebug, ieof, lstop,  &
                        zerror)
!--purpose: read material parameters for temperature calculation
      implicit none
!--subroutine arguments:
      character*(*), intent(in) :: inpdir
      integer, intent(in)       :: lunit, ncase, modul, ldebug
      integer, intent(inout)    :: linenr, ieof
      type(t_ic)                :: ic
      type(t_trackdata)         :: trk
      type(t_simpack)           :: spck
      logical, intent(in)       :: lstop
      logical, intent(inout)    :: zerror
!--local variables:
      integer, parameter :: mxnval = 10
      integer            :: ints(mxnval), nval, ierror, is_wheel, i_ftype
      logical            :: flags(mxnval)
      real(kind=8)       :: dbles(mxnval)
      character(len=256) :: strngs(mxnval)

      if (ic%ztrack.eq.1 .or. ic%ztrack.eq.3) then

         if (ic%config.le.1) then

            ! get the track design dimensions (using 'd' for either gaugsq or raily0)

            call readline(lunit, ncase, linenr, 'track design dimensions', 'ddda', ints, dbles, flags,  &
                           strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
            zerror = zerror .or. (ierror.ne.0)

            trk%gauge_height = dbles(1)
            if (trk%gauge_height.gt.0d0) then
               trk%gauge_seqnum = nint(dbles(2))
               trk%track_gauge  = dbles(3)
            else
               trk%rail_y0      = dbles(2)
               trk%rail_z0      = dbles(3)
            endif
            trk%cant_angle   = dbles(4)

         else

            ! get the roller rig design dimensions (using 'd' for either gaugsq or raily0)

            call readline(lunit, ncase, linenr, 'roller rig design dimensions', 'dddd', ints, dbles,    &
                           flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
            zerror = zerror .or. (ierror.ne.0)

            trk%gauge_height = dbles(1)
            if (trk%gauge_height.gt.0d0) then
               trk%gauge_seqnum = nint(dbles(2))
               trk%track_gauge  = dbles(3)
            else
               trk%rail_y0      = dbles(2)
               trk%rail_z0      = dbles(3)
            endif
            trk%nom_radius   = dbles(4)

         endif

      endif

      if (ic%ztrack.eq.3) then

         ! Z1=3: one rail profile used for current side of track; get filename and configuration data

         call profile_read_config(lunit, trk%rai%prr, 'rail', ncase, linenr, ldebug, ieof, lstop, zerror)
         zerror = zerror .or. (ierror.ne.0)

         ! check that variable profiles are used with absolute rail placement

         call profile_get_filetype(trk%rai%prr%fname, is_wheel, i_ftype, ldebug)

         if (i_ftype.eq.FTYPE_SLICES .and. trk%gauge_height.gt.0d0) then
            zerror = .true.
            write(lout, 5511) trk%gauge_height, trim(trk%rai%prr%fname)
            write(   *, 5511) trk%gauge_height, trim(trk%rai%prr%fname)
 5511       format (' Input: ERROR. The gauge point computation (GAUGHT =',f8.2,                        &
                    ') cannot be used when using a',/,15x,'variable profile ("',a,'").')
         endif

         ! read the rail profile

         if (.not.zerror) then
            call profile_read_file(trk%rai%prr, inpdir, 0, ic%x_profil, ic%x_readln, lstop)
            zerror = zerror .or. (trk%rai%prr%ierror.ne.0)
         endif

      endif

      if (modul.ne.modul_spck .and. ic%ztrack.ge.2 .and. ic%ztrack.le.3) then

         ! read the track deviations for current side of the track

         call readline(lunit, ncase, linenr, 'current rail deviations', 'ddadda', ints, dbles, flags,   &
                        strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
         zerror = zerror .or. (ierror.ne.0)

         trk%rai%dy    = dbles(1)
         trk%rai%dz    = dbles(2)
         trk%rai%roll  = dbles(3)
         trk%rai%vy    = dbles(4)
         trk%rai%vz    = dbles(5)
         trk%rai%vroll = dbles(6)

      endif

      if (modul.ne.modul_spck .and. ic%ztrack.ge.2 .and. ic%ztrack.le.3 .and. ic%force1.eq.3) then

         ! read the track deflection parameters for current side of the track

         call readline(lunit, ncase, linenr, 'current rail deflection parameters', 'dddd', ints, dbles, &
                        flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
         zerror = zerror .or. (ierror.ne.0)

         trk%ky_rail   = dbles(1)
         trk%fy_rail   = dbles(2)
         trk%kz_rail   = dbles(3)
         trk%fz_rail   = dbles(4)

      endif

   end subroutine trackdata_input

!------------------------------------------------------------------------------------------------------------

   subroutine wheelset_input(inpdir, lunit, ncase, linenr, modul, ic, ws, trk, spck, ldebug, ieof,      &
                lstop, zerror)
!--purpose: read input for E-digit: wheel-set dimensions, profile, position and velocity
      implicit none
!--subroutine arguments:
      character*(*), intent(in) :: inpdir
      integer, intent(in)       :: lunit, ncase, modul, ldebug
      integer, intent(inout)    :: linenr, ieof
      type(t_ic)                :: ic
      type(t_wheelset)          :: ws
      type(t_trackdata)         :: trk
      type(t_simpack)           :: spck
      logical, intent(in)       :: lstop
      logical, intent(inout)    :: zerror
!--local variables:
      integer, parameter :: mxnval = 10
      integer            :: ints(mxnval), nval, ierror
      logical            :: flags(mxnval)
      real(kind=8)       :: dbles(mxnval)
      character(len=16)  :: types, namside
      character(len=256) :: strngs(mxnval)

      ! get the wheel-set dimensions

      if (ic%ewheel.eq.3 .or. ic%ewheel.eq.5) then

         call readline(lunit, ncase, linenr, 'wheel-set dimensions', 'ddd', ints, dbles, flags,          &
                        strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
         zerror = zerror .or. (ierror.ne.0)
         ws%flback_dist  = dbles(1)
         ws%flback_pos   = dbles(2)
         ws%nom_radius   = dbles(3)

         zerror = zerror .or. .not.check_range ('NOMRAD', ws%nom_radius, 1d-3, 1d20)

      endif

      ! get the wheel profile

      if (ic%ewheel.eq.3 .or. ic%ewheel.eq.5) then

         ! E1=3,5: get profile filename and configuration data

         call profile_read_config(lunit, ws%whl%prw, 'wheel', ncase, linenr, ldebug, ieof, lstop, zerror)

         ! read the wheel profile, store in wheel data

         call profile_read_file(ws%whl%prw, inpdir, 1, ic%x_profil, ic%x_readln, lstop)
         zerror = zerror .or. (ws%whl%prw%ierror.ne.0)

         if (.false. .and. ic%discns1_eff.eq.5) then
            call write_log(' D=5: converting wheel to variable profile...')
            call profile_make_varprof(ws%whl%prw, 1, ic%x_profil)
         endif

      endif

      ! get the wheel-set position and orientation

      if (modul.ne.modul_spck .and. ic%ewheel.ge.1 .and. ic%ewheel.le.5) then

         call readline(lunit, ncase, linenr, 'wheel-set position and orientation', 'dddaaa', ints,      &
                        dbles, flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
         zerror = zerror .or. (ierror.ne.0)
         if (ic%config.le.1) then
            ws%s      = dbles(1)
            ws%x      = 0d0
         else
            ws%s      = 0d0
            ws%x      = dbles(1)        ! C1=4,5: PITCH_ROL prescribed
         endif
         ws%y      = dbles(2)
         if (ic%norm.le.0) then
            ws%z      = dbles(3)
            ws%fz_inp = 0d0
         else
            ws%fz_inp = dbles(3)        ! N=1: FZ_TR prescribed
         endif
         ws%roll     = dbles(4)
         ws%yaw      = dbles(5)
         ws%pitch    = dbles(6)
      endif

      ! get the wheel-set velocity

      if (modul.ne.modul_spck .and. ic%ewheel.ge.2 .and. ic%ewheel.le.5) then

         if     (ic%config.le.1 .and. (ic%force1.le.0 .or. ic%force1.eq.3)) then
            types = 'dddaaa'    ! 1st VS_WS,      6th VPITCH_WS
         elseif (ic%config.le.1) then
            types = 'dddaad'    ! 1st VS_WS,      6th FX_WS/MY_WS
         elseif (ic%config.ge.4 .and. (ic%force1.le.0 .or. ic%force1.eq.3)) then
            types = 'addaaa'    ! 1st VPITCH_ROL, 6th VPITCH_WS
         elseif (ic%config.ge.4) then
            types = 'addaad'    ! 1st VPITCH_ROL, 6th FX_WS/MY_WS
         endif

         call readline(lunit, ncase, linenr, 'wheel-set velocity and rotation', types, ints, dbles,     &
                        flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
         zerror = zerror .or. (ierror.ne.0)

         if (ic%config.le.1) then
            ws%vs          = dbles(1)
            trk%vpitch_rol = 0d0
         else
            ws%vs          = 0d0
            trk%vpitch_rol = dbles(1)   ! C1=4,5: VPITCH_ROL prescribed
         endif
         ws%vy       = dbles(2)
         ws%vz       = dbles(3)
         ws%vroll    = dbles(4)
         ws%vyaw     = dbles(5)
         if (ic%force1.le.0 .or. ic%force1.eq.3) then
            ws%vpitch = dbles(6)
         elseif (ic%force1.eq.1) then
            ws%fx_inp = dbles(6)        ! F=1: FX_WS prescribed
         elseif (ic%force1.eq.2) then
            ws%my_inp = dbles(6)        ! F=2: MY_WS prescribed
         endif
      endif

      ! get flexible wheel-set deviations

      if (ic%ewheel.eq.4 .or. ic%ewheel.eq.5) then

         namside = 'current'

         ! read the wheel position deviations for current side of wheel-set

         call readline(lunit, ncase, linenr, trim(namside) // ' wheel position deviations', 'dddaaa',   &
                       ints, dbles, flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
         zerror = zerror .or. (ierror.ne.0)

         ws%whl%dx     = dbles(1)
         ws%whl%dy     = dbles(2)
         ws%whl%dz     = dbles(3)
         ws%whl%droll  = dbles(4)
         ws%whl%dyaw   = dbles(5)
         ws%whl%dpitch = dbles(6)

         ! read the wheel velocity deviations for current side of wheel-set

         call readline(lunit, ncase, linenr, trim(namside) // ' wheel velocity deviations', 'dddaaa',   &
                       ints, dbles, flags, strngs, mxnval, nval, ldebug, ieof, lstop, ierror)
         zerror = zerror .or. (ierror.ne.0)

         ws%whl%vx     = dbles(1)
         ws%whl%vy     = dbles(2)
         ws%whl%vz     = dbles(3)
         ws%whl%vroll  = dbles(4)
         ws%whl%vyaw   = dbles(5)
         ws%whl%vpitch = dbles(6)
      
      endif ! E=5,7

   end subroutine wheelset_input

!------------------------------------------------------------------------------------------------------------

   subroutine wr_write_inp(ncase, wtd)
!--purpose: Write input for the current W/R contact case to the inp-file, unit linp. 
      implicit none
!--subroutine arguments:
      integer                   :: ncase
      type(t_ws_track)          :: wtd
!--local variables:
      integer                  :: vldcmze, xhgiaowr, cpbtnfs, psflcin
      real(kind=8)             :: dflt_turn
      character(len=1)         :: namside
      character(len=10)        :: strng
      type(t_ic)               :: tmp_ic
      type(t_wheel),   pointer :: my_wheel

      call ic_copy(wtd%ic, tmp_ic)      ! create local copy of configuration

      associate(ic    => tmp_ic,   mater => wtd%mater, discr => wtd%discr, kin   => wtd%kin,            &
                fric  => wtd%fric, solv  => wtd%solv,  ws    => wtd%ws,    trk   => wtd%trk)

      if (ic%wrtinp.ge.2) then
         ic%frclaw_inp  = wtd%fric%frclaw_eff     ! write friction to inp-file
         ic%discns1_inp = wtd%ic%discns1_eff      ! write discretisation parameters to inp-file
         ic%gencr_inp   = wtd%mater%gencr_eff     ! write material parameters to inp-file
         ic%mater2      = 0
         if (max(mater%cdampn,mater%cdampt).ge.1d-15) ic%mater2 = 2
         ic%ewheel      = max(3, ic%ewheel)       ! TODO: check writing of wheelset flexibility
         ic%ztrack      = 3
      endif

      if (ncase.gt.0) write(linp,'(a,i0)') '% Next case ', ncase

      ! write control integers

      call ic_pack (1, cpbtnfs, vldcmze, xhgiaowr, ic)

      write(linp, 1101) cpbtnfs
      if (ic%mater2.le.0) then
         write(linp, 1102) vldcmze
      else
         write(linp, 1103) vldcmze, ic%mater2
      endif
      if (ic%xflow.le.0) then
         write(linp, 1104) xhgiaowr
      else
         write(linp, 1105) xhgiaowr
      endif

      ! write debug parameters

      if (ic%xflow.ge.1) then
         call ic_pack_dbg(psflcin, ic)
         if (ic%x_readln.le.0) then
            write(linp, 1106) psflcin
         else
            write(linp, 1107) psflcin, ic%x_readln
         endif
      endif

 1101 format(i8.7, 6x, '  C-P-B-T-N-F-S       CONFIG, PVTIME, BOUND,  TANG,   NORM,   FORCE,  STRESS')
 1102 format(i8.7, 6x, '  V-L-D-C-M-Z-E       VARFRC, FRCLAW, DISCNS, INFLCF, MATER,  ZTRACK, EWHEEL')
 1103 format(i8.7,i3,4x,' V-L-D-C-M-Z-E  M    VARFRC, FRCLAW, DISCNS, INFLCF, MATER,  ZTRACK, EWHEEL, MATER2')
 1104 format(i8.7, 6x, '  H-G-I-A-O-W-R         HEAT, GAUSEI, IESTIM, MATFIL, OUTPUT, FLOW,   RETURN')
 1105 format(i8.7, 6x, 'X-H-G-I-A-O-W-R  XFLOW, HEAT, GAUSEI, IESTIM, MATFIL, OUTPUT, FLOW,   RETURN')
 1106 format(i8.7, 6x, '  P-S-F-L-C-I-N       PROFIL, SMOOTH, FORCE,  LOCATE, CPATCH, INFLCF, NMDBG' )
 1107 format(i8.7,i3,4x,' P-S-F-L-C-I-N  R    PROFIL, SMOOTH, FORCE,  LOCATE, CPATCH, INFLCF, NMDBG,  READLN')

      ! write parameters for the iterative solution algorithms

      if (ic%gausei_inp.ne.1) write(linp, 1201) solv%maxgs, solv%maxin, solv%maxnr, solv%maxout, solv%eps
      if (ic%gausei_inp.eq.2 .or. ic%gausei_inp.eq.3)                                                   &
         write(linp, 1211) solv%omegah, solv%omegas, solv%inislp, solv%omgslp
      if (ic%gausei_inp.eq.4) write(linp, 1212) solv%inislp, solv%omgslp
      if (ic%gausei_inp.eq.5) then
         if (solv%gd_meth.eq.1) solv%fdecay =  1d0
         if (solv%gd_meth.eq.2) solv%fdecay = -real(solv%kdown)
         write(linp, 1215) solv%fdecay, solv%betath, solv%kdowfb, solv%d_ifc, solv%d_lin, solv%d_cns,  &
            solv%d_slp, solv%pow_s
      endif

 1201 format( 4i6,    3x, es8.1, 14x, 'MAXGS,  MAXIN,  MAXNR,  MAXOUT, EPS')
 1211 format( 2g12.4, i6, g12.4,  7x, 'OMEGAH, OMEGAS, INISLP, OMGSLP')
 1212 format(         i6, g12.4, 40x, 'INISLP, OMGSLP')
 1215 format( 6f8.3,  10x, 'FDECAY, D_IFC/LIN/CNS, D_SLP, POW_S')

      ! write friction description

      call fric_wrtinp(linp, ic%frclaw_inp, fric)

      ! write information for influence coefficients, esp. material constants

      if (ic%gencr_inp.eq.2 .or. ic%gencr_inp.eq.3 .or. ic%gencr_inp.eq.4) then

         ! write IF-correction method

         if (ic%gencr_inp.eq.4) then
            write(linp, 5100) mater%if_meth, mater%if_ver
 5100       format (2i6, 46x, 'IF_METH, VARIANT')
         endif

         ! write elastic material constants

         write(linp, 5200) mater%poiss, mater%gg
 5200    format (2(f8.4, 4x), 2g12.4, 10x, 'POISS 1,2,  GG 1,2')

         ! write visco-elastic material constants

         if (ic%mater.eq.1) then
            write(linp, 5301) mater%fg(1), mater%fg(2), mater%tc(1), mater%tc(2)
 5301       format (4g12.4, 10x, 'FG 1,2,  TC 1,2')
         endif

         ! write flexibility of simplified theory and slope reduction parameters

         if (ic%mater.eq.2) then
            write(linp, 5402) mater%flx(1), mater%k0_mf, mater%alfamf, mater%betamf
 5402       format (4g12.4, 10x, 'FLX, K0_MF, ALFA, BETA')
         elseif (ic%mater.eq.3 .or. ic%mater.eq.5) then
            write(linp, 5403) mater%k0_mf, mater%alfamf, mater%betamf
 5403       format (3g12.4, 22x, 'K0_MF, ALFAMF, BETAMF')
         endif

         ! write parameters of elasto-plastic interface layer

         if (ic%mater.eq.4) then
            write(linp, 5501) mater%gg3, mater%laythk, mater%tau_c0, mater%k_tau
 5501       format (4g12.4, 10x, 'GG3, LAYTHK, TAU_C0, K_TAU')
         endif
      endif

      ! write material parameters for damping model

      if (ic%mater2.eq.2) then
         write(linp,5611) mater%cdampn, mater%cdampt, mater%dfnmax, mater%dftmax
 5611    format (4g12.4, 10x, 'CDAMPN, CDAMPT, DFNMAX, DFTMAX')
      endif

      ! write material parameters for temperature model

      if (ic%heat.eq.3) then
         write(linp,5711) mater%bktemp(1), mater%heatcp(1), mater%lambda(1), mater%dens(1), 1,1,1,1
         write(linp,5711) mater%bktemp(2), mater%heatcp(2), mater%lambda(2), mater%dens(2), 2,2,2,2
 5711    format (4g12.4, 10x, 'BKTEMP',i1,', HEATCP',i1,', LAMBDA',i1,', DENSITY',i1)
      endif

      if (ic%heat.eq.3 .and. ic%mater.eq.4) then
         write(linp,5712) mater%betapl
 5712    format (1g12.4, 46x, 'BETAPL')
      endif

      ! write information needed for the grid discretization

      if (ic%discns1_inp.ge.2 .and. ic%discns1_inp.le.9) then
         strng     = fmt_gs(10, 3, 3, discr%dqrel)
         ! write dist_turn only if different from default value
         dflt_turn = 2d0 * discr%dist_sep - 1d0 * discr%dist_comb
         if (abs(discr%dist_turn-dflt_turn).le.1d-3) then
            write(linp, 6101) discr%dx, discr%ds, strng, discr%angl_sep*180d0/pi, discr%dist_sep, &
                   discr%dist_comb
         else
            write(linp, 6102) discr%dx, discr%ds, strng, discr%angl_sep*180d0/pi, discr%dist_sep, &
                   discr%dist_comb, discr%dist_turn
         endif
 6101    format( 2f9.4, a, f8.1,'d', 2f9.4, 3x, 'DX, DS, DQREL, A_SEP, D_SEP, D_COMB')
 6102    format( 2f9.4, a, f8.1,'d', 3f6.2, 3x, 'DX, DS, DQREL, A_SEP, D_SEP, D_COMB, D_TURN')
      endif

      ! write information on track geometry

      if (ic%ztrack.eq.1 .or. ic%ztrack.eq.3) then

         if     (ic%config.le.1 .and. trk%gauge_height.gt.0d0) then
            write(linp, 7101) trk%gauge_height, trk%gauge_seqnum, trk%track_gauge, trk%cant_angle, 'CANT'
         elseif (                     trk%gauge_height.gt.0d0) then
            write(linp, 7101) trk%gauge_height, trk%gauge_seqnum, trk%track_gauge, trk%nom_radius, 'NOMRAD'
         elseif (ic%config.le.1) then
            write(linp, 7102) trk%gauge_height, trk%rail_y0, trk%rail_z0, trk%cant_angle, 'CANT'
         else
            write(linp, 7102) trk%gauge_height, trk%rail_y0, trk%rail_z0, trk%nom_radius, 'NOMRAD'
         endif
 7101    format ( g14.6, i6, 2g14.6, 10x, 'GAUGHT, GAUGSQ, GAUGWD, ',a)
 7102    format (4g14.6,  2x, 'GAUGHT, RAILY0, RAILZ0, ',a)

      endif

      ! write rail profile filename(s)

      if (ic%ztrack.eq.3) call profile_write_config(linp, wtd%meta, trk%rai%prr, 0)

      ! write rail deviations

      if (ic%ztrack.ge.2 .and. ic%ztrack.le.3) then

         write(linp, 7201) trk%rai%dy, trk%rai%dz, trk%rai%roll, trk%rai%vy, trk%rai%vz, trk%rai%vroll

 7201    format (3g14.6, 16x, '% DY, DZ, DROLLR', /, 3g14.6, 16x  '  VY, VZ, VROLL', :,i1)

      endif

      ! write massless rail deflection parameters

      if (ic%ztrack.eq.3 .and. ic%force1.eq.3) then

         write(linp, 7301) trk%ky_rail, trk%fy_rail, trk%kz_rail, trk%fz_rail

 7301    format (4g14.6, 2x, 'KYRAIL, FYRAIL, KZRAIL, FZRAIL')

      endif

      ! write information on the wheelset geometry

      if (ic%ewheel.eq.3 .or. ic%ewheel.eq.5) then

         write(linp, 8101) ws%flback_dist, ws%flback_pos, ws%nom_radius
 8101    format (3g14.6, 16x, 'FBDIST, FBPOS, NOMRAD')

      endif

      ! write wheel profile filename

      if (ic%ewheel.eq.3 .or. ic%ewheel.eq.5) then
         call profile_write_config(linp, wtd%meta, ws%whl%prw, 1)
      endif

      ! write information on the wheelset position and velocity
      ! Using kinematic constants fz, fx, fy

      if (ic%ewheel.ge.1) then
         if (ic%config.le.1 .and. ic%norm.eq.0) then
            write(linp, 8201) ws%s, ws%y, ws%z,      'S', 'Z'
         elseif (ic%config.le.1) then
            write(linp, 8201) ws%s, ws%y, ws%fz_inp, 'S', 'FZ'
         elseif (ic%config.ge.4 .and. ic%norm.eq.0) then
            write(linp, 8201) ws%x, ws%y, ws%z,      'X', 'Z'
         elseif (ic%config.ge.4) then
            write(linp, 8201) ws%x, ws%y, ws%fz_inp, 'X', 'FZ'
         endif
         write(linp, 8203) ws%roll, ws%yaw, ws%pitch
 8201    format( 3(g15.8, 1x), 13x, '% ',a,', Y, ',a)
 8203    format( 3(g15.8, 1x), 13x, '% ROLL, YAW, PITCH')
      endif

      if (ic%ewheel.ge.2) then
         if (ic%config.le.1) then
            write(linp, 8301) ws%vs,          ws%vy,     ws%vz, 'VS'
         else
            write(linp, 8301) trk%vpitch_rol, ws%vy,     ws%vz, 'VPITCH_ROL'
         endif
 8301    format( 3(g15.8, 1x), 13x, '% ',a,', VY, VZ')

         if (ic%force1.eq.0 .or. ic%force1.eq.3) then
            write(linp, 8303) ws%vroll, ws%vyaw,   ws%vpitch, 'VPITCH'
         elseif (ic%force1.eq.1) then
            write(linp, 8303) ws%vroll, ws%vyaw,   ws%fx_inp, 'FX'
         elseif (ic%force1.eq.2) then
            write(linp, 8303) ws%vroll, ws%vyaw,   ws%my_inp, 'MY'
         endif
 8303    format( 2(g15.8, 1x), g16.9, 12x, '% VROLL, VYAW, ',a)
      endif

      ! write flexible wheelset deviations

      if (ic%ewheel.eq.4 .or. ic%ewheel.eq.5) then
         namside = ' '
         my_wheel => ws%whl
         write(linp, 8401) my_wheel%dx, my_wheel%dy, my_wheel%dz, namside,                              &
                           my_wheel%droll, my_wheel%dyaw, my_wheel%dpitch, namside,                     &
                           my_wheel%vx, my_wheel%vy, my_wheel%vz, namside,                              &
                           my_wheel%vroll, my_wheel%vyaw, my_wheel%vpitch, namside

 8401    format (3g14.6, 16x, '% DX, DY, DZ',a, /, 3g14.6, 16x, '  DROLL, DYAW, DPITCH',a, /,           &
                 3g14.6, 16x, '% VX, VY, VZ',a, /, 3g14.6, 16x, '  VROLL, VYAW, VPITCH',a)
      endif

      ! write the points in which the subsurface elastic field must be calculated. 

      if (ic%stress.ge.2) call subsurf_wrtinp(linp, ic, wtd%subs)

      ! write empty line to separate from next case

      write(linp, '(1x)')

      end associate
   end subroutine wr_write_inp

!------------------------------------------------------------------------------------------------------------

end module m_wr_input
