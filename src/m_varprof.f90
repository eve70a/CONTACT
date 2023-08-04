!------------------------------------------------------------------------------------------------------------
! m_varprof - working with 'variable profiles' consisting of multiple 'profile slices'
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_varprof

use m_hierarch_data
use m_profile

implicit none
private

   ! Debugging for module m_varprof

   public  varprof_set_debug

   integer  :: ldebug    =  0    ! local level of debugging
   integer  :: ii_debug  = -1    ! output point for which detailed info is requested (-1 = none)
   integer  :: iel_debug = -1    ! input element for which detailed info is requested (-1 = none)

   ! extension of type t_profile for variation in longitudinal direction, using multiple slices

   public  profile_get_filetype
   public  profile_read_file
   public  varprof_read_file
   public  profile_read_slice

   private varprof_make_vj
   public  varprof_make_bspline
   public  varprof_get_z_at_xy
   public  varprof_get_loc_z_at_xy
   public  varprof_intpol_grid_old
   public  varprof_intpol_grid
   public  varprof_intpol_xunif

   private varprof_find_slices_xrange
   public  varprof_available_x

contains

!------------------------------------------------------------------------------------------------------------

subroutine varprof_set_debug(new_ldebug, new_ii_debug, new_iel_debug)
!--function: enable/disable debug output of varprof routines
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
      write(bufout,'(a,i3,2(a,i7))') ' varprof:  debugging level =',ldebug,', ii_debug =', ii_debug,    &
                ', iel_debug =', iel_debug
      call write_log(1, bufout)
   endif

end subroutine varprof_set_debug

!------------------------------------------------------------------------------------------------------------

   subroutine profile_get_filetype(fname, is_wheel, i_ftype, idebug)
!--purpose: Determine from profile filename extension if its a rail, wheel or unknown (0 / 1 / -1),
!           in Miniprof or Simpack format (0 / 1).
      implicit none
!--subroutine arguments:
      character*(*)             :: fname
      integer,      intent(in)  :: idebug
      integer,      intent(out) :: is_wheel, i_ftype
!--local variables:
      integer              :: ix
      character(len=20)    :: file_ext

      ! find index of last '.' in filename

      ix = index(fname, '.', back=.true.) + 1
      file_ext = to_lower( trim(fname(ix:)) )

      if (idebug.ge.4) then
         write(bufout,*) 'File extension is "', trim(file_ext), '"'
         call write_log(1, bufout)
      endif

      ! Determine file format: Simpack uses extensions prr or prw, Miniprof uses ban or whl

      if     ( file_ext(1:4).eq.'slcs' ) then
         i_ftype = FTYPE_SLICES
      elseif ( file_ext(1:3).eq.'prr' .or. file_ext(1:3).eq.'prw' ) then
         i_ftype = FTYPE_SIMPACK
      elseif ( file_ext(1:3).eq.'ban' .or. file_ext(1:3).eq.'whl' ) then
         i_ftype = FTYPE_MINIPROF
      else
         i_ftype = FTYPE_PLAIN2D
      endif

      ! Determine whether file is rail or wheel

      if     ( file_ext(1:3).eq.'prr' .or. file_ext(1:3).eq.'ban' .or. file_ext(1:4).eq.'slcs' ) then
         is_wheel = 0
      elseif ( file_ext(1:3).eq.'prw' .or. file_ext(1:3).eq.'whl' ) then
         is_wheel = 1
      else
         is_wheel = -1
         if (idebug.ge.5) then
            write(bufout,'(3a)') ' INFO: unknown file extension for profile file "', trim(fname), '"'
            call write_log(1, bufout)
         endif
      endif

   end subroutine profile_get_filetype

!------------------------------------------------------------------------------------------------------------

   subroutine profile_read_file(prf, dirnam, is_wheel_arg, x_profil, x_readln, lstop)
!--purpose: reads a single w/r profile (Simpack, Miniprof, Ascii) or a variable profile (slcs)
      implicit none
!--subroutine parameters:
      type(t_profile)              :: prf
      integer,       intent(in)    :: is_wheel_arg, x_profil, x_readln
      logical,       intent(in)    :: lstop
      character*(*), intent(in)    :: dirnam
!--local variables:
      integer              :: is_wheel, i_ftype
      logical              :: fill_spline

      call timer_start(itimer_profil)

      ! determine the type of file

      call profile_get_filetype(prf%fname, is_wheel, i_ftype, x_profil)

      ! switch on basis of the type of file

      if (i_ftype.eq.FTYPE_SLICES) then

         call varprof_read_file(prf, dirnam, is_wheel_arg, x_profil, x_readln, lstop)

      else

         fill_spline = .true.
         call profile_read_slice(prf%fname, dirnam, prf%grd_data, prf%ext_data, prf, is_wheel_arg,      &
                        fill_spline, x_profil, x_readln, lstop)

      endif

      call timer_stop(itimer_profil)

   end subroutine profile_read_file

!------------------------------------------------------------------------------------------------------------

   subroutine varprof_read_file(vprf, dirnam, is_wheel_arg, x_profil, x_readln, lstop)
!--purpose: read, process and store contents of a so-called slices file
      implicit none
!--subroutine parameters:
      type(t_profile)              :: vprf
      integer,       intent(in)    :: is_wheel_arg, x_profil, x_readln
      logical,       intent(in)    :: lstop
      character*(*), intent(in)    :: dirnam
!--local variables:
      integer,      parameter :: mxnval = 20
      character*3,  parameter :: commnt = '!%"'
      real(kind=8), parameter :: tiny_ds = 1d-3  ! cf. tiny_dt in m_bspline
      integer          :: ints(mxnval), idum(1), lslcs, ncase, linenr, nval, ieof, my_profil
      integer          :: ia, ik, ip, islc, ix, nout
      integer          :: iout, sub_ierror
      logical          :: flags(mxnval), any_kyield, fill_spline, zerror
      real(kind=8)     :: dbles(mxnval), rdum(1), s_offset, s_scale, dst,                               &
                          xout(5001), yout(5001), zout(5001)
      character*256    :: strngs(mxnval), fmtstr, descrp, fulnam, slcdir
      real(kind=8), dimension(:), allocatable :: tmp

      associate(fname      => vprf%fname,    nslc       => vprf%nslc,                                   &
                nfeat      => vprf%nfeat,    nkink      => vprf%nkink,    naccel    => vprf%naccel,     &
                s_method   => vprf%s_method, ierror     => vprf%ierror)

      ! determine the relative folder from the slices filename

      slcdir = ' '
      ix = index_pathsep(fname, back=.true.)
      if (ix.gt.0) then
         slcdir = fname(1:ix-1)
      endif

      ! determine full path-name, pre-pending dirnam when necessary

      if (dirnam.ne.' ' .and. index_pathsep(fname).le.0) then
         fulnam = trim(dirnam) // path_sep // fname
         slcdir = trim(dirnam) // path_sep // trim(slcdir)
      else
         fulnam = fname
      endif

      if (x_profil.ge.1) then
         write(bufout,*) 'Reading file "',trim(fulnam), '" with profile slices'
         call write_log(1, bufout)
      endif

      ! use free unit number defined in m_print_output

      lslcs = ltmp

      open(lslcs, file=fulnam, status='old', err=995)
      linenr = 0
      goto 997

      ! Error handling:

 995  continue
         write(bufout,*) 'ERROR: cannot open file "', trim(fulnam),'" for reading.'
         call write_log(1, bufout)
         call abort_run()
 997  continue

      ierror =  0
      ncase  =  0
      ieof   = -1  ! eof considered an error
      nslc   =  0

      ! Read offset and scaling factor

      call readline(lslcs, ncase, linenr, 'offset and scaling factor', 'dd', ints, dbles, flags,        &
                strngs, mxnval, nval, x_readln, ieof, lstop, ierror)

      s_offset = dbles(1)
      s_scale  = dbles(2)

      if (x_profil.ge.1) then
         write(bufout,'(a,g12.4,a,g12.4)') ' s_offset =', s_offset,', s_scale =',s_scale
         call write_log(1, bufout)
      endif

      ! Read number of slices, number of parts, kinks, accelerations, spline method

      call readline(lslcs, ncase, linenr, 'number of slices', 'i', ints, dbles, flags, strngs, mxnval,  &
                nval, x_readln, ieof, lstop, ierror)

      nslc     = ints(1)

      call readline(lslcs, ncase, linenr, 'number of parts, kinks, accel', 'iii', ints, dbles, flags,   &
                strngs, mxnval, nval, x_readln, ieof, lstop, ierror)

      nfeat    = ints(1)
      if (nfeat.le.1) then      ! NFEAT <= 1 means 1 part with default S_F = [ 0, 1e6 ]
         nkink  = 0
         naccel = 0
      else
         nkink    = ints(2)
         naccel   = ints(3)
      endif

      call readline(lslcs, ncase, linenr, 'longitudinal spline method', 'i', ints, dbles, flags,        &
                strngs, mxnval, nval, x_readln, ieof, lstop, ierror)

      s_method = ints(1)

      ! Check dimensions

      zerror =             .not.check_irng ('NSLC',    nslc,    1, MAX_NUM_SLCS)
      zerror = zerror .or. .not.check_irng ('NFEAT',   nfeat,   0, 999)
      if (nfeat.ge.2) zerror = zerror .or. .not.check_irng ('NKINK',   nkink,   0, nfeat-2)
      if (nfeat.ge.2) zerror = zerror .or. .not.check_irng ('NACCEL',  naccel,  0, nfeat-2)

      if (s_method.ne.SPL2D_INTPOL .and. s_method.ne.SPL2D_APPROX) then
         write(bufout,'(3(a,i4),a)') ' ERROR: invalid S_METHOD =',s_method,' must be', SPL2D_INTPOL,    &
                ' (intpol) or',SPL2D_APPROX,' (approx)'
         call write_log(1, bufout)
         zerror = .true.
      endif

      if (zerror) ierror = 151

      if (ierror.eq.0) then

         ! allocate arrays slc_s and slc_nam for nslc slices

         call reallocate_arr(vprf%slc_s, nslc)
         call reallocate_arr(vprf%slc_nam, nslc)

         ! Read s-position and profile filename per slice

         do islc = 1, nslc

            write(descrp,'(a,i4)') 'profile slice',islc

            call readline(lslcs, ncase, linenr, descrp, 'ds', ints, dbles, flags, strngs, mxnval, nval,    &
                   x_readln, ieof, lstop, ierror)

            vprf%slc_s(islc)   = s_scale * (s_offset + dbles(1))
            vprf%slc_nam(islc) = trim(strngs(1))

            if (x_profil.ge.3) then
               write(bufout,'(a,i4,a,f11.3,3a)') ' slice ',islc,': slc_s =',vprf%slc_s(islc),              &
                   ', fname = "', trim(vprf%slc_nam(islc)),'"'
               call write_log(1, bufout)
            endif

         enddo ! islc

         if (x_profil.ge.1) then
            write(bufout,'(a,i3,a)') ' obtained ',nslc,' slices'
            call write_log(1, bufout)
         endif

      endif

      ! check that s-positions are in strictly increasing order

      if (ierror.eq.0) then
         do islc = 1, nslc-1
            dst = vprf%slc_s(islc+1) - vprf%slc_s(islc)
            if (dst.le.tiny_ds) then
               write(bufout,'(a,2(a,i4),2(a,g12.4))') ' ERROR: slice s-positions should be strictly ',  &
                   'increasing. islc=',islc,',', islc+1,': s=', vprf%slc_s(islc),',', vprf%slc_s(islc+1)
               call write_log(1, bufout)
               ierror = 152
            endif
         enddo
      endif

      ! Read information on parts in lateral direction: kinks and accelerations

      if (ierror.eq.0 .and. nkink.ge.1) then
         call reallocate_arr(vprf%kink_if, nkink)
         call read1darr(lslcs, ncase, linenr, 'kinks between parts', 'i', nkink, vprf%kink_if,           &
                        rdum, x_readln, lstop, ierror)

         if (x_profil.ge.-1) then
            if (nkink.eq.1) then
               write(bufout,'(a,i5)') ' there is 1 kink after part p =', vprf%kink_if(1)
            else
               write(bufout,'(a,i4,a,20i4)') ' there are',nkink,' kinks after parts p =',               &
                                                                       (vprf%kink_if(ik), ik=1,nkink)
            endif
            call write_log(1, bufout)
         endif
      endif

      if (ierror.eq.0 .and. naccel.ge.1) then
         call reallocate_arr(vprf%accel_if, naccel)
         call read1darr(lslcs, ncase, linenr, 'accel between parts', 'i', naccel, vprf%accel_if,         &
                        rdum, x_readln, lstop, ierror)

         if (x_profil.ge.-1) then
            if (naccel.eq.1) then
               write(bufout,'(a,i5)') ' there is 1 acceleration after part p =', vprf%accel_if(1)
            else
               write(bufout,'(a,i4,a,20i4)') ' there are',naccel,' accelerations after parts p =',      &
                                                                     (vprf%accel_if(ia), ia=1,naccel)
            endif
            call write_log(1, bufout)
         endif
      endif

      ! check kink and accel positions

      if (ierror.eq.0) then
         fmtstr = '(a,i3,a)'
         if (max(nkink,naccel).le.99) fmtstr = '(a,i2,a)'
         if (max(nkink,naccel).le. 9) fmtstr = '(a,i1,a)'
         zerror = .false.

         do ik = 1, nkink
            write(descrp, fmtstr) 'KINK(',ik,')'
            zerror = zerror .or. .not.check_irng(trim(descrp), vprf%kink_if(ik), 1, nfeat-2)
         enddo
         do ia = 1, naccel
            write(descrp, fmtstr) 'ACCEL(',ia,')'
            zerror = zerror .or. .not.check_irng(trim(descrp), vprf%accel_if(ia), 1, nfeat-2)
         enddo

         if (zerror) ierror = 153
      endif

      ! Read information on parts per slice: s_f feature positions

      if (ierror.eq.0) then
         if (nfeat.le.1) then   ! NFEAT<=1: default S_F at start/after end of profile

            nfeat = 2
            call reallocate_arr(vprf%slc_s_f, nslc, nfeat)
            do islc = 1, nslc
               vprf%slc_s_f(islc,1:2) = (/ 0d0, 1d6 /)
            enddo

         else

            call reallocate_arr(vprf%slc_s_f, nslc, nfeat)
            allocate(tmp(nfeat+1))
            do islc = 1, nslc

               write(descrp,'(a,i4)') 'features s_f for slice', islc
               call read1darr(lslcs, ncase, linenr, descrp, 'd', nfeat+1, idum, tmp, x_readln, lstop,   &
                        ierror)

               tmp(1) = s_scale * (s_offset + tmp(1))

               if (abs(tmp(1)-vprf%slc_s(islc)).ge.tiny_ds) then
                  ierror = 154
                  write(bufout,'(2(a,f12.3))') ' ERROR: incorrect s_slc =', tmp(1),                     &
                        ' in feature information, expecting s_slc =', vprf%slc_s(islc)
                  call write_log(1, bufout)
               endif
               vprf%slc_s_f(islc,1:nfeat) = tmp(2:nfeat+1)

            enddo ! islc
            deallocate(tmp)

            if (x_profil.ge.-1) then
               write(bufout,'(2(a,i4),a)') ' feature information for',nslc,' slices with',nfeat,        &
                              ' features per slice'
               call write_log(1, bufout)
               if (nfeat.le.21) then
                  do islc = 1, nslc
                     write(bufout,'(i4,21f9.2)') islc, (vprf%slc_s_f(islc,ip), ip=1, nfeat)
                     call write_log(1, bufout)
                  enddo
               else
                  do islc = 1, nslc
                     write(bufout,'(i4,10f9.2,a,10f9.2)') islc, (vprf%slc_s_f(islc,ip), ip=1,10),       &
                                                    ' ... ', (vprf%slc_s_f(islc,ip), ip=nfeat-9, nfeat)
                     call write_log(1, bufout)
                  enddo
               endif
            endif
         endif ! nfeat<=1

         ! TODO: check feature positions

      endif

      close(lslcs)

      ! read files for each profile slice

      any_kyield = .false.

      do islc = 1, nslc

         if (ierror.eq.0) then

            if (x_profil.ge.3) then
               write(bufout,'(a,i4,3a)') ' reading slice',islc,', file="',trim(vprf%slc_nam(islc)),'"'
               call write_log(1, bufout)
            endif

            ! destroy grid/grid-function if allocated before

            if (associated(vprf%slc_grd(islc)%g)) then
               call grid_destroy(vprf%slc_grd(islc)%g)
               deallocate(vprf%slc_grd(islc)%g)
            endif

            ! allocate space for grid/gridfunc datastructures

            allocate(vprf%slc_grd(islc)%g)

            ! read file for slice islc

            my_profil = x_profil
            if (islc.gt.1) my_profil = my_profil - 1
            fill_spline = .true.

            call profile_read_slice(vprf%slc_nam(islc), slcdir, vprf%slc_grd(islc)%g, vprf%ext_data,       &
                        vprf, is_wheel_arg, fill_spline, my_profil, x_readln, lstop)
            ! ierror returned via vprf%ierror

            if (vprf%has_kyield) any_kyield = .true.

            if (ierror.ne.0) then
               write(bufout,'(a,i4,3a,i4,a)') ' An error occurred for slice islc=',islc,' file="',      &
                        trim(vprf%slc_nam(islc)),'" (',ierror,'), aborting.'
               call write_log(1, bufout)
            endif

         endif

      enddo

      ! build the 2d spline representation

      if (ierror.eq.0) then
         ! call write_log(' ...varprof_make_bspline')
         ! call varprof_set_debug(1)
         call varprof_make_bspline(vprf)
         ! call varprof_set_debug(0)
         ! call bspline2d_print(vprf%spl2d, 'spl2d', 5)

         if (.false.) then
            ! testing get_z_at_xy...
            nout    = 5001
            do iout = 1, nout
               xout(iout) = (iout-1) * 0.01d0
               yout(iout) = 24.90d0
               zout(iout) = -1234d0
            enddo
            call write_log(' ...bspline_get_z_at_xy')
            call bspline_set_debug(1)
            call spline_set_debug(1)
            call bspline_get_z_at_xy(vprf%spl2d, nout, xout, yout, zout, sub_ierror)
            call spline_set_debug(0)
            call bspline_set_debug(0)

            do iout = 1, nout
               write(bufout,'(a,i4,3(a,f14.8))') ' i=',iout,', (x,y)=(',xout(iout),',',yout(iout),         &
                           '): z=',zout(iout)
               call write_log(1, bufout)
            enddo
            call abort_run()
         endif
      endif

      ! copy data for first slice to the output profile

      if (ierror.eq.0) then
         call grid_copy(vprf%slc_grd(1)%g, vprf%grd_data, with_spline=.true.)
      endif

      ! warn in case kyield was provided for any of the profile slices

      if (ierror.eq.0 .and. any_kyield) then
         call write_log(' Variable profiles do not support kyield hardness data.')
         vprf%has_kyield = .false.
         call gf3_destroy(vprf%ext_data)
      endif

      if (ierror.ne.0 .and. lstop) then
         call write_log(' Errors found. Aborting.')
         call abort_run()
      endif
      end associate

   end subroutine varprof_read_file

!------------------------------------------------------------------------------------------------------------

   subroutine profile_read_slice(fname, slcdir, prf_grd, prf_ext, prf_opt, is_wheel_arg, fill_spline,   &
                        x_profil, x_readln, lstop)
!--purpose: read, process and store contents of a single w/r profile file (Simpack, Miniprof, Ascii)
      implicit none
!--subroutine parameters:
      character*(*), intent(in)    :: fname, slcdir
      type(t_grid)                 :: prf_grd   ! output: profile data
      type(t_gridfnc3)             :: prf_ext   ! output: extra data
      type(t_profile)              :: prf_opt   ! input: configuration options
      integer,       intent(in)    :: is_wheel_arg, x_profil, x_readln
      logical,       intent(in)    :: fill_spline, lstop
!--local variables:
      real(kind=8), parameter :: pi     = 4d0*atan(1d0)
      integer              :: is_wheel, i_ftype, npoint, icol_kyield
      character(len=5)     :: nam_rw
      character(len=256)   :: fulnam
      real(kind=8), dimension(:,:), pointer :: points => NULL()

      associate(mirror_y   => prf_opt%mirror_y,   mirror_z   => prf_opt%mirror_z,                       &
                sclfac     => prf_opt%sclfac,     zig_thrs   => prf_opt%zig_thrs,                       &
                has_kyield => prf_opt%has_kyield, ierror     => prf_opt%ierror)

      ierror     = 0
      has_kyield = .false.

      if (is_wheel_arg.eq.0) then
         nam_rw = 'rail '
      else
         nam_rw = 'wheel'
      endif

      ! determine full path-name, pre-pending slcdir when necessary

      if (slcdir.ne.' ' .and. index_pathsep(fname).ne.1) then
         fulnam = trim(slcdir) // path_sep // fname
      else
         fulnam = fname
      endif
      call set_platform_filesep(fulnam)

      if (x_profil.ge.1) call write_log(' reading ' // nam_rw // ' profile "' // trim(fulnam) // '"')

      ! get type of file from filename extension

      call profile_get_filetype(fulnam, is_wheel, i_ftype, x_profil)
      if (is_wheel.eq.-1) is_wheel = is_wheel_arg

      ! read file with appropriate routine 

      if     ( i_ftype.eq.FTYPE_SIMPACK .and. is_wheel.ne.is_wheel_arg ) then

         ! incorrect file extension: not supported

         write(bufout,'(5a)') ' ERROR: incorrect file extension for ',trim(nam_rw), ' profile file "',  &
                                trim(fname),'"'
         call write_log(1, bufout)
         ierror = 301

      elseif ( i_ftype.eq.FTYPE_SIMPACK ) then

         if (x_profil.ge.3) then
            write(bufout,*) 'Opening Simpack ',trim(nam_rw),' profile "',trim(fname),'"...'
            call write_log(1, bufout)
         endif
         call read_simpack_profile(fname, fulnam, npoint, points, is_wheel, mirror_y, mirror_z, sclfac, &
                        zig_thrs, x_profil, x_readln, lstop, ierror)
         icol_kyield = 0

      elseif ( i_ftype.eq.FTYPE_MINIPROF .and. is_wheel.eq.is_wheel_arg ) then

         if (x_profil.ge.3) then
            write(bufout,*) 'Opening Miniprof ',trim(nam_rw),' profile "',trim(fname),'"...'
            call write_log(1, bufout)
         endif
         call read_miniprof_profile(fname, fulnam, i_ftype, npoint, points, is_wheel, mirror_y, mirror_z, &
                        sclfac, zig_thrs, icol_kyield, x_profil, x_readln, lstop, ierror)

      elseif ( i_ftype.eq.FTYPE_PLAIN2D ) then

         ! unknown file extension: assume plain X,Y columns, process via miniprof routine 
         ! using is_wheel_arg as provided by user instead of is_wheel obtained from file extension

         if (x_profil.ge.3) then
            write(bufout,'(5a,/,a)') ' INFO: unknown file extension for ',trim(nam_rw),                 &
                        ' profile file "',     trim(fname),'",', '       assuming plain X,Y-values'
            call write_log(2, bufout)
         endif

         if (x_profil.ge.3) then
            write(bufout,*) 'Opening plain ',trim(nam_rw),' profile "',trim(fname),'"...'
            call write_log(1, bufout)
         endif
         call read_miniprof_profile(fname, fulnam, i_ftype, npoint, points, is_wheel_arg, mirror_y,     &
                        mirror_z, sclfac, zig_thrs, icol_kyield, x_profil, x_readln, lstop, ierror)

      else

         ! unknown file extension: not supported

         write(bufout,'(5a)') ' ERROR: unknown file extension for ',trim(nam_rw), ' profile file "',    &
                                trim(fname),'"'
         call write_log(1, bufout)
         ierror = 301

      endif

      if (ierror.eq.0) then

         ! store points in profile-grid, modify, and create spline representation

         call profile_finish_grid(fname, prf_grd, prf_ext, prf_opt, npoint, points, is_wheel,           &
                icol_kyield, fill_spline, x_profil)

      endif

      call destroy_arr(points)
      end associate
   end subroutine profile_read_slice

!------------------------------------------------------------------------------------------------------------

   subroutine varprof_make_vj(vprf, ds_max2d, vj)
!--purpose: for a variable profile, determine the knot-vector vj for resampling of slices
      implicit none
!--subroutine arguments:
      type(t_profile),            target      :: vprf
      real(kind=8),               intent(in)  :: ds_max2d
      real(kind=8), dimension(:), allocatable :: vj
!--local variables:
      real(kind=8), parameter   :: tiny = 1d-6
      integer                   :: islc, ib0, ib1, ib, ip, i0, i1, j, npart, nseg_p, nseg
      real(kind=8)              :: ds, dv, len_part, len_part_max(vprf%nfeat-1)
      type(t_grid), pointer     :: prr

      call reallocate_arr(vprf%iseg_p, vprf%nfeat)
      call reallocate_arr(vprf%slc_ib, vprf%nslc, 2)

      associate( nfeat  => vprf%nfeat,  s_f    => vprf%slc_s_f, iseg_p => vprf%iseg_p,          &
                 ierror => vprf%ierror )
      npart = nfeat - 1

      ! check and adapt the feature information
      ! determine the maximum length per part over all of the slices

      len_part_max(1:npart) = -1d0

      do islc = 1, vprf%nslc

         prr => vprf%slc_grd(islc)%g

         ! determine the active features ('breaks') [ib0:ib1] for each slice

         ib0 = 1
         do while(ib0.lt.vprf%nfeat .and. s_f(islc,ib0).lt.-tiny)
            ib0 = ib0 + 1
         enddo

         ib1 = vprf%nfeat
         do while(ib1.gt.1 .and. s_f(islc,ib1).lt.-tiny)
            ib1 = ib1 - 1
         enddo

         if (ib1.le.ib0) then
            ierror = 161
            write(bufout,'(a,i4)') ' ERROR: no features defined for slice', islc
            call write_log(1, bufout)
         endif

         ! store (ib0,ib1) in vprf for later use

         vprf%slc_ib(islc,1) = ib0
         vprf%slc_ib(islc,2) = ib1

         ! shift the active features to the s-coordinates as used in the profile

         if (ldebug.ge.3) then
            write(bufout,'(a,i4,4(a,f8.2))') ' Slice',islc,': s_f=', s_f(islc,ib0),',', s_f(islc,ib1),  &
                                ', prr_s=', prr%s_prf(1),',', prr%s_prf(prr%ntot)
            call write_log(1, bufout)
         endif

         do ib = ib0, ib1
            s_f(islc,ib) = s_f(islc,ib) + prr%s_prf(1)
         enddo

         ! shift the last feature to the end of the profile if needed

         s_f(islc,ib1) = min(s_f(islc,ib1), prr%s_prf(prr%ntot))

         ! check that s_f are in increasing order and in range of profile
         ! note that s_f(1)>=s_prf(1) and s_f(end)<=s_prf(end)

         do ib = ib0+1, ib1
            ds = s_f(islc,ib) - s_f(islc,ib-1)
            if (ds.le.tiny) then
               ierror = 162
               write(bufout,'(2a,i4,2(a,i3),2(a,f8.2),a)') ' ERROR: feature positions s_f must be ',    &
                        'strictly increasing. Slice',islc,', features',ib-1, ' and',ib,', s_f=',        &
                        s_f(islc,ib-1),',', s_f(islc,ib),'.'
               call write_log(1, bufout)
            endif
         enddo

         ! determine the length of each part in this slice, update the maximum over all slices

         do ip = ib0, ib1-1
            len_part = s_f(islc,ip+1) - s_f(islc,ip)
            len_part_max(ip) = max(len_part_max(ip), len_part)
         enddo

      enddo ! islc: check/update feature info

      ! determine the number of segments nseg per part and the first index in the vector vj

      iseg_p(1) = 1
      do ip = 1, npart
         nseg_p       = 1 + int( (len_part_max(ip)+tiny) / ds_max2d )
         iseg_p(ip+1) = iseg_p(ip) + nseg_p
      enddo

      ! get the total number of points used for resampling

      nseg = iseg_p(npart+1) - 1

      ! determine the vj-positions to be used for resampling

      allocate(vj(nseg+1))

      do ip = 1, npart
         i0 = iseg_p(ip)
         i1 = iseg_p(ip+1)  ! npnt = nseg+1
         dv = 1d0
         do j = i0, i1
            vj(j) = real(j)
         enddo

         if (ldebug.ge.1) then
            write(bufout,'(a,i3,3(a,i4),3(a,f6.1),a)') ' part',ip,':',i1-i0,' segments, points [',i0,   &
                        ',',i1,'], vj = [', vj(i0), ',', vj(i0+1),' ...', vj(i1),' ]'
            call write_log(1, bufout)
         endif
      enddo

      end associate
   end subroutine varprof_make_vj

!------------------------------------------------------------------------------------------------------------

   subroutine varprof_make_bspline(vprf)
!--purpose: for a variable profile, build a 2D tensor B-spline representation in struct vprf%spl2d
!           using approximation of slices (resampling) with nsplv points in the lateral direction
      implicit none
!--subroutine arguments:
      type(t_profile)           :: vprf
!--local variables:
      real(kind=8), parameter   :: ds_max2d = 0.5d0
      logical            :: has_xdata, use_approx
      integer            :: nsplx, nsplv, islc, jy, ib0, ib1, ip, i0, i1, npnt, sub_ierror
      real(kind=8)       :: len_s, len_v
      logical,      dimension(:,:), allocatable :: mask
      real(kind=8), dimension(:),   allocatable :: vj, sj, y1d, z1d
      real(kind=8), dimension(:,:), allocatable :: x2d, y2d, z2d

      associate(ierror => vprf%ierror)

      if (vprf%nslc.le.0) then

         if (ldebug.ge.2) call write_log(' varprof_make_bspline: no slices, nothing to do')

      else

         if (ldebug.ge.1) call write_log(' setting up 2D B-spline for variable profile')

         ! set up the sampling for the lateral direction, based on feature information

         if (ierror.eq.0) then
            call varprof_make_vj(vprf, ds_max2d, vj)
         endif

         ! allocate work arrays for resampling of slices

         nsplx = vprf%nslc
         nsplv = size(vj, 1)
         allocate(sj(nsplv), y1d(nsplv), z1d(nsplv))
         allocate(x2d(1,1), y2d(nsplx,nsplv), z2d(nsplx,nsplv), mask(nsplx,nsplv))

         x2d(1,1) = -1234d0

         do jy = 1, nsplv
            y2d(1:nsplx,jy) = 0d0
            z2d(1:nsplx,jy) = 0d0
         enddo

         ! fill mask array for active points

         do islc = 1, vprf%nslc
            mask(islc,1:nsplv) = .false.
            ib0 = vprf%slc_ib(islc,1)   ! features [ib0:ib1]
            ib1 = vprf%slc_ib(islc,2)
            i0  = vprf%iseg_p(ib0)      ! corresponding points [i0:i1]
            i1  = vprf%iseg_p(ib1)
            mask(islc,i0:i1) = .true.
         enddo

         ! loop over all slices, evaluate 1d B-spline at uniform s-sampling

         if (ierror.eq.0) then

            do islc = 1, vprf%nslc

               ! get the active features ('breaks') [ib0:ib1] for the slice

               ib0 = vprf%slc_ib(islc,1)
               ib1 = vprf%slc_ib(islc,2)

               ! for each part that is used, determine the scaling

               do ip = ib0, ib1-1

                  ! length of part within profile

                  len_s = vprf%slc_s_f(islc,ip+1) - vprf%slc_s_f(islc,ip)

                  ! point numbers after resampling

                  i0    = vprf%iseg_p(ip)
                  i1    = vprf%iseg_p(ip+1)
                  len_v = vj(i1) - vj(i0)          ! == 1d0

                  ! stretching of s-values

                  sj(i0:i1) = vprf%slc_s_f(islc,ip) + len_s * (vj(i0:i1) - vj(i0)) / len_v

                  if (ldebug.ge.3) then
                     write(bufout,'(a,i4,a,i3,2(a,i4),4(a,f8.3),a)') ' slice',islc,', part',ip,         &
                           ': points [', i0,',',i1,'], vj = [',vj(i0),',',vj(i1),', sj = [', sj(i0),    &
                           ',',sj(i1),']'
                     call write_log(1, bufout)
                  endif

               enddo

               ! evaluate spline for slice islc at sj-positions

               associate(gslc => vprf%slc_grd(islc)%g, spl => vprf%slc_grd(islc)%g%spl)

               i0    = vprf%iseg_p(ib0)
               i1    = vprf%iseg_p(ib1)
               npnt  = i1 - i0 + 1

               ! note: using constant extrapolation to permit sj > spl%s by round-off error

               call spline_eval(gslc%spl, ikYDIR, npnt, sj(i0:i1), sub_ierror, f_eval=y1d)
               y2d(islc,i0:i1) = y1d(1:npnt)

               call spline_eval(gslc%spl, ikZDIR, npnt, sj(i0:i1), sub_ierror, f_eval=z1d)
               z2d(islc,i0:i1) = z1d(1:npnt)

               end associate ! gslc

            enddo ! islc

            if (ldebug.ge.5) then
               call print_2d_real(vprf%nslc, nsplv, y2d, 'y2d', 10, 'g12.4')
               call print_2d_real(vprf%nslc, nsplv, z2d, 'z2d', 10, 'g12.4')
            endif
         endif ! ierror.eq.0

         ! create 2D tensor spline for data x2d, y2d, z2d

         if (ierror.eq.0) then
            ! call write_log('bspline_make_2d_bspline...')
            ! call bspline_set_debug(3)
            has_xdata  = .false.
            use_approx = (vprf%s_method.eq.SPL2D_APPROX)
            call bspline_make_2d_bspline(vprf%spl2d, vprf%nslc, nsplv, vprf%slc_s, vj, x2d, y2d, z2d,   &
                                        has_xdata, use_approx, ierror, mask)
            ! call bspline_set_debug(0)
         endif

         deallocate(vj, sj, y1d, z1d, x2d, y2d, z2d)

      endif ! no slices, nothing to do
      end associate

   end subroutine varprof_make_bspline

!------------------------------------------------------------------------------------------------------------

   subroutine varprof_get_z_at_xy(vprf, s_ws, nout, xout, yout, zout, dzout, my_ierror, exterval)
!--purpose: interpolate variable profile at n positions (x,y), producing both z and dz/dx
      implicit none
!--subroutine arguments:
      type(t_profile)           :: vprf
      integer,      intent(in)  :: nout
      real(kind=8), intent(in)  :: s_ws, xout(nout), yout(nout)
      real(kind=8), intent(out) :: zout(nout), dzout(nout)
      integer,      intent(out) :: my_ierror
      real(kind=8), optional    :: exterval
!--local variables:
      integer           :: iout, i0, i1, sub_ierror
      real(kind=8)      :: dx, x_s(nout), zout1(nout)

      my_ierror = 0

      if (ldebug.ge.3) then
         write(bufout,'(a,i5,a)') ' varprof_get_z_at_xy: get z at', nout,' xy-positions'
         call write_log(1, bufout)
      endif

      x_s(1:nout) = s_ws + xout(1:nout)

      ! if (ldebug.ge.3) call splineget_set_debug(3, ii_debug)
      call bspline_get_z_at_xy(vprf%spl2d, nout, x_s, yout, zout, sub_ierror, exterval)
      ! call splineget_set_debug(0)

      ! use finite difference to get dz/dx

      dx = 0.01d0
      x_s(1:nout) = x_s(1:nout) + dx

      call bspline_get_z_at_xy(vprf%spl2d, nout, x_s, yout, zout1, sub_ierror, exterval)

      dzout(1:nout) = (zout1(1:nout) - zout(1:nout)) / dx

      if (ldebug.ge.3) then
         if (ldebug.eq.3) then
            i0 = max(1, ii_debug)
            i1 = min(nout, ii_debug)
         else
            i0 = 1
            i1 = nout
         endif

         do iout = i0, i1
            write(bufout,'(a,i4,5(a,f12.4))') '  iout=',iout,': (x,y)=(',xout(iout),',',yout(iout),     &
                   '), z=',zout(iout),', z1=',zout1(iout),', dz/dx=',dzout(iout)
            call write_log(1, bufout)
         enddo
      endif

   end subroutine varprof_get_z_at_xy

!------------------------------------------------------------------------------------------------------------

   subroutine varprof_get_loc_z_at_xy( vprf, s_ws, m_loc, g_out, my_ierror )
!--purpose: for a variable profile, compute intermediate slices, rotate to local coordinates, and
!           evaluate local z at the positions (x,s) of the output grid
      implicit none
!--subroutine arguments:
      type(t_profile)             :: vprf
      real(kind=8),   intent(in)  :: s_ws
      type(t_marker)              :: m_loc
      type(t_grid)                :: g_out
      integer,        intent(out) :: my_ierror
!--local variables:
      real(kind=8),   parameter :: defval = 999d0, tiny = 1d-6
      type(t_grid)              :: g_curv, g_slc, g_trim, g_out1d
      integer                   :: ix, iy, iy0, iy1, ii, ii0, nx, ns, nline, nbefor, nafter, sub_ierror
      real(kind=8), dimension(:), allocatable :: xi
      real(kind=8)              :: xdum2d(1,1)

      my_ierror = 0
      nx = g_out%nx
      allocate(xi(nx))

      if (vprf%nslc.le.0) then

         call write_log(' Internal Error: (varprof_get_loc_z): not a variable profile.')
         call abort_run()

      endif

      ! Note: g_out has cp-coordinates, transform to rail/track and then to track-curve 

      xi(1:nx) = s_ws + m_loc%x() + g_out%x(1:nx)

      if (ldebug.ge.2) then
         call write_log(' output positions s_ws+xi along track curve:')
         write(bufout,'(10( 20(f10.3),:,/ ))') (xi(ix), ix=1, min(200,nx))
         nline = int( (min(200,nx)-1) / 20 ) + 1
         call write_log(nline, bufout)
      endif

      ! check for xi-positions before first slice or after last slice

      nbefor = 0
      nafter = 0
      do ix = 1, nx
         if (xi(ix).lt.vprf%slc_s(1)) nbefor = nbefor + 1
         if (xi(ix).gt.vprf%slc_s(vprf%nslc)) nafter = nafter + 1
      enddo

      if (ldebug.ge.1) then
         if (nbefor+nafter.gt.0) then
            write(bufout,'(a,i4,2(a,f10.3),a)') ' varprof_get_loc_z: requested profile at', nx,         &
                   ' positions s_ws+x in [', xi(1),',',xi(nx),']'
            call write_log(1, bufout)
         endif
         if (nbefor.gt.0) then
            write(bufout,'(a,i4,a,f10.3,a)') ' WARNING:',nbefor,' positions lie before 1st slice at s_1=', &
                   vprf%slc_s(1),', using constant extrapolation'
            call write_log(1, bufout)
         endif
         if (nafter.gt.0) then
            write(bufout,'(a,i4,a,f10.3,a)') ' WARNING:',nafter,' positions lie after last slice at s_n=', &
                   vprf%slc_s(vprf%nslc),', using constant extrapolation'
            call write_log(1, bufout)
         endif
      endif

      ! clip xi-values outside [s1,sn]

      if (nbefor+nafter.gt.0) then
         do ix = 1, nx
            xi(ix) = max(vprf%slc_s(1), min(vprf%slc_s(vprf%nslc), xi(ix)))
         enddo
      endif

      ! create temporary grid with nx slices, ns profile points

      ns = vprf%spl2d%nbrkv
      call grid_create_curvil(g_curv, nx, ns, lies_in_oyz=.false.)

      ! extrude x-positions from first row of g_out

      do iy = 1, ns
         ii0 = (iy-1) * nx
         g_curv%x(ii0+1:ii0+nx) = g_out%x(1:nx)
      enddo

      ! evaluate 2D spline at given x-values and s-breakpoints

      if (my_ierror.eq.0) then
         if (ldebug.ge.1) then
            write(bufout,'(2(a,i4),a)') ' varprof_get_loc_z: evaluate 2d-spline at', nx,' x', ns,       &
                        ' positions'
            call write_log(1, bufout)
         endif

         call bspline_eval2d_prod(vprf%spl2d, nx, ns, xi, vprf%spl2d%vbrk, xdum2d, g_curv%y, g_curv%z,  &
                        sub_ierror, defval)
         my_ierror = sub_ierror
         if (my_ierror.ne.0) call write_log(' Varprof_get_loc_z: Error after eval2d')
         ! call grid_print(g_curv, 'g_curv', 5)
      endif

      ! create temporary grid for one column of output grid

      call grid_create_curvil(g_out1d, 1, g_out%ny, lies_in_oyz=.true.)
      do iy = 1, g_out%ny
         ii = 1 + (iy-1)*nx
         g_out1d%y(iy) = g_out%y(ii)
      enddo

      ! create temporary grid with 1 slice, ns profile points

      call grid_create_curvil(g_slc, 1, ns, lies_in_oyz=.true.)

      ! call marker_print(m_loc, 'mcp_r', 5)

      ! loop over slices, compute spline, evaluate using 1D spline_get_loc_xz

      do ix = 1, nx

         if (ldebug.ge.5) then
            write(bufout,'(a,i4,a,f8.3)') ' varprof_get_loc_z: compute slice ix=',ix,', xtr=',g_out%x(ix)
            call write_log(1, bufout)
         endif

         ! copy slice ix from 2d g_curv into 1d g_slc

         do iy = 1, ns
            ii = ix + (iy-1)*nx
            g_slc%y(iy) = g_curv%y(ii)   ! y(ix,1:ns)
            g_slc%z(iy) = g_curv%z(ii)   ! z(ix,1:ns)
         enddo

         ! compute spline (fast hack, no smoothing, no kinks)

         iy0 = ns
         iy1 =  1
         do iy = 1, ns
            if (g_slc%y(iy).lt.defval-tiny) then
               iy0 = min(iy, iy0)       ! first non-default
               iy1 = max(iy, iy1)       ! last non-default
            endif
         enddo

         if (ldebug.ge.2) then
            write(bufout,'(3(a,i4),a)') 'ix=', ix,': active region iy=[',iy0,',',iy1,'], trimming...'
            call write_log(1, bufout)
         endif

         call grid_trim(g_slc, g_trim, 1, 1, iy0, iy1, with_spline=.false.)
         ! call grid_print(g_trim, 'g_trim', 5)

         call grid_make_arclength(g_trim, sub_ierror)
         my_ierror = sub_ierror

         if (ldebug.ge.1) then
            call spline_set_debug(1)
            call write_log(' loc_z_at_xy: kchk...')
            call grid_make_ppspline(g_trim, 0d0, .false., sub_ierror, k_chk=10)
            call spline_set_debug(0)
         else
            call grid_make_ppspline(g_trim, 0d0, .false., sub_ierror)
         endif
         my_ierror = sub_ierror

         ! if (ix.eq.31) call spline_print(g_trim%spl, 'g_trim31', 5)

         ! evaluate spline at local y-positions of output grid

         call spline_get_loc_xz_at_y( g_trim%spl, m_loc, g_out%ny, g_out1d%y, sub_ierror,               &
                        zloc_out=g_out1d%z )

         ! store results in g_out

         do iy = 1, g_out%ny
            ii = ix + (iy-1)*nx
            g_out%z(ii) = g_out1d%z(iy)
            if (ix.eq.-31) then
               write(bufout,'(2(a,i3),2(a,f10.4))') 'ix=',ix,', iy=',iy,': y=',g_out1d%y(iy),', z=',    &
                        g_out1d%z(iy)
               call write_log(1, bufout)
            endif
         enddo

      enddo ! slice ix

      call grid_destroy(g_curv)
      call grid_destroy(g_slc)
      call grid_destroy(g_trim)
      call grid_destroy(g_out1d)
      deallocate(xi)

   end subroutine varprof_get_loc_z_at_xy

!------------------------------------------------------------------------------------------------------------

   subroutine varprof_intpol_grid_old(vprf, m_trk, s_ws, g_out)
!--purpose: interpolate between the slices of a variable profile as needed to fill the grid g_out.
      implicit none
!--subroutine arguments:
      type(t_profile), target      :: vprf
      type(t_marker),  intent(in)  :: m_trk
      real(kind=8),    intent(in)  :: s_ws      ! offset between slice s and track x_tr
      type(t_grid)                 :: g_out
!--local variables:
      real(kind=8),   parameter :: defval = 999d0, tiny = 1d-6
      type(t_grid)              :: g_curv, g_slc, g_trim, g_out1d
      integer                   :: ix, iy, iy0, iy1, ii, ii0, nx, ns, nline, nbefor, nafter,            &
                                   my_ierror, sub_ierror
      real(kind=8)              :: xdum2d(1,1)
      real(kind=8), dimension(:), allocatable :: xi

      my_ierror = 0
      if (vprf%nslc.le.0) then

         call write_log(' Internal Error: (varprof_intpol_grid_old): not a variable profile.')
         call abort_run()

      endif

      if (ldebug.ge.1) then
         write(bufout,'(a,i4,a,i5,a)') ' interpolating variable profile to grid of',g_out%nx,' x',      &
                     g_out%ny, ' points'
         call write_log(1, bufout)
      endif

      nx = g_out%nx
      allocate(xi(nx))

      xi(1:nx) = s_ws + g_out%x(1:nx)

      if (ldebug.ge.2) then
         call write_log(' output positions s_ws+xi along track curve:')
         write(bufout,'(10( 20(f10.3),:,/ ))') (xi(ix), ix=1, min(200,nx))
         nline = int( (min(200,nx)-1) / 20 ) + 1
         call write_log(nline, bufout)
      endif

      ! check for xi-positions before first slice or after last slice

      nbefor = 0
      nafter = 0
      do ix = 1, nx
         if (xi(ix).lt.vprf%slc_s(1)) nbefor = nbefor + 1
         if (xi(ix).gt.vprf%slc_s(vprf%nslc)) nafter = nafter + 1
      enddo

      if (nbefor+nafter.gt.0) then
         write(bufout,'(a,i4,2(a,f10.3),a)') ' varprof_intpol_grid_old: requested profile at', nx,      &
                ' positions s_ws+x in [', xi(1),',',xi(nx),']'
         call write_log(1, bufout)
      endif
      if (nbefor.gt.0) then
         write(bufout,'(a,i4,a,f10.3,a)') ' WARNING:',nbefor,' positions lie before 1st slice at s_1=', &
                vprf%slc_s(1),', using constant extrapolation'
         call write_log(1, bufout)
      endif
      if (nafter.gt.0) then
         write(bufout,'(a,i4,a,f10.3,a)') ' WARNING:',nafter,' positions lie after last slice at s_n=', &
                vprf%slc_s(vprf%nslc),', using constant extrapolation'
         call write_log(1, bufout)
      endif

      ! clip xi-values outside [s1,sn]

      if (nbefor+nafter.gt.0) then
         do ix = 1, nx
            xi(ix) = max(vprf%slc_s(1), min(vprf%slc_s(vprf%nslc), xi(ix)))
         enddo
      endif

      ! create temporary grid with nx slices, ns profile points

      ns = vprf%spl2d%nbrkv
      call grid_create_curvil(g_curv, nx, ns, lies_in_oyz=.false.)

      ! extrude x-positions from first row of g_out

      do iy = 1, ns
         ii0 = (iy-1) * nx
         g_curv%x(ii0+1:ii0+nx) = g_out%x(1:nx)
      enddo

      ! evaluate 2D spline at given x-values and s-breakpoints

      if (my_ierror.eq.0) then
         if (ldebug.ge.1) then
            write(bufout,'(2(a,i4),a)') ' varprof_intpol_grid_old: evaluate 2d-spline at', nx,' x',     &
                        ns, ' positions'
            call write_log(1, bufout)
         endif

         call bspline_eval2d_prod(vprf%spl2d, nx, ns, xi, vprf%spl2d%vbrk, xdum2d, g_curv%y, g_curv%z,  &
                        sub_ierror, defval)
         my_ierror = sub_ierror
         if (my_ierror.ne.0) call write_log(' Varprof_intpol_grid_old: Error after eval2d')
         ! call grid_print(g_curv, 'g_curv', 5)
      endif

      ! create temporary grid for one column of output grid

      call grid_create_curvil(g_out1d, 1, g_out%ny, lies_in_oyz=.true.)
      do iy = 1, g_out%ny
         ii = 1 + (iy-1)*nx
         g_out1d%y(iy) = g_out%y(ii)
      enddo

      ! create temporary grid with 1 slice, ns profile points

      call grid_create_curvil(g_slc, 1, ns, lies_in_oyz=.true.)

      ! loop over slices, compute spline, evaluate using 1D spline_get_xz_at_y

      do ix = 1, nx

         if (ldebug.ge.5) then
            write(bufout,'(a,i4,a,f8.3)') ' varprof_intpol_grid_old: compute slice ix=',ix,', xtr=',    &
                        g_out%x(ix)
            call write_log(1, bufout)
         endif

         ! copy slice ix from 2d g_curv into 1d g_slc

         do iy = 1, ns
            ii = ix + (iy-1)*nx
            g_slc%y(iy) = g_curv%y(ii)   ! y(ix,1:ns)
            g_slc%z(iy) = g_curv%z(ii)   ! z(ix,1:ns)
         enddo

         ! compute spline (fast hack, no smoothing, no kinks)

         iy0 = ns
         iy1 =  1
         do iy = 1, ns
            if (g_slc%y(iy).lt.defval-tiny) then
               iy0 = min(iy, iy0)       ! first non-default
               iy1 = max(iy, iy1)       ! last non-default
            endif
         enddo

         if (ldebug.ge.2) then
            write(bufout,'(3(a,i4),a)') 'ix=', ix,': active region iy=[',iy0,',',iy1,'], trimming...'
            call write_log(1, bufout)
         endif

         call grid_trim(g_slc, g_trim, 1, 1, iy0, iy1, with_spline=.false.)
         ! call grid_print(g_trim, 'g_trim', 5)

         ! convert interpolated slice to track coordinate system

         call cartgrid_2glob(g_trim, m_trk)

         call grid_make_arclength(g_trim, sub_ierror)
         my_ierror = sub_ierror

         if (ldebug.ge.1) then
            call spline_set_debug(1)
            call write_log(' intpol_grid_old: kchk...')
            call grid_make_ppspline(g_trim, 0d0, .false., sub_ierror, k_chk=10)
            call spline_set_debug(0)
         else
            call grid_make_ppspline(g_trim, 0d0, .false., sub_ierror)
         endif
         my_ierror = sub_ierror

         ! if (ix.eq.31) call spline_print(g_trim%spl, 'g_trim31', 5)

         ! evaluate spline at local y-positions of output grid

         call spline_get_xz_at_y( g_trim%spl, g_out%ny, g_out1d%y, sub_ierror, zout=g_out1d%z,          &
                exterval=defval )

         ! store results in g_out

         do iy = 1, g_out%ny
            ii = ix + (iy-1)*nx
            g_out%z(ii) = g_out1d%z(iy)
            if (ix.eq.-31) then
               write(bufout,'(2(a,i3),2(a,f10.4))') 'ix=',ix,', iy=',iy,': y=',g_out1d%y(iy),', z=',    &
                        g_out1d%z(iy)
               call write_log(1, bufout)
            endif
         enddo

      enddo ! slice ix

      call grid_destroy(g_curv)
      call grid_destroy(g_slc)
      call grid_destroy(g_trim)
      call grid_destroy(g_out1d)
      deallocate(xi)

   end subroutine varprof_intpol_grid_old

!------------------------------------------------------------------------------------------------------------

   subroutine varprof_intpol_grid(vprf, m_trk, s_ws, g_out)
!--purpose: interpolate between the slices of a variable profile as needed to fill the grid g_out.
      implicit none
!--subroutine arguments:
      type(t_profile), target      :: vprf
      type(t_marker),  intent(in)  :: m_trk
      real(kind=8),    intent(in)  :: s_ws      ! offset between slice s and track x_tr
      type(t_grid)                 :: g_out
!--local variables:
      real(kind=8),   parameter :: defval = 999d0
      integer                   :: ix, iy, ii, nx, ny, nline, nbefor, nafter,            &
                                   my_ierror
      real(kind=8), dimension(:),   allocatable :: xi, yj, xtmp, ztmp
      real(kind=8), dimension(:,:), allocatable :: zij

      my_ierror = 0
      if (vprf%nslc.le.0) then
         call write_log(' Internal Error: (varprof_intpol_grid): not a variable profile.')
         call abort_run()
      endif

      if (ldebug.ge.1) then
         write(bufout,'(a,i4,a,i5,a)') ' interpolating variable profile to grid of',g_out%nx,' x',      &
                     g_out%ny, ' points'
         call write_log(1, bufout)
      endif

      ! convert output grid temporarily to rail coordinate system

      call cartgrid_2loc(g_out, m_trk)

      nx = g_out%nx
      ny = g_out%ny
      allocate(xi(nx), yj(ny), zij(nx,ny), xtmp(ny), ztmp(ny))

      xi(1:nx) = s_ws + g_out%x(1:nx)
      do iy = 1, ny
         yj(iy) = g_out%y(iy*nx)
      enddo

      if (ldebug.ge.2) then
         call write_log(' output positions s_ws+xi along track curve:')
         write(bufout,'(10( 20(f10.3),:,/ ))') (xi(ix), ix=1, min(200,nx))
         nline = int( (min(200,nx)-1) / 20 ) + 1
         call write_log(nline, bufout)
      endif

      ! check for xi-positions before first slice or after last slice

      nbefor = 0
      nafter = 0
      do ix = 1, nx
         if (xi(ix).lt.vprf%slc_s(1)) nbefor = nbefor + 1
         if (xi(ix).gt.vprf%slc_s(vprf%nslc)) nafter = nafter + 1
      enddo

      if (nbefor+nafter.gt.0) then
         write(bufout,'(a,i4,2(a,f10.3),a)') ' varprof_intpol_grid: requested profile at', nx,         &
                ' positions s_ws+x in [', xi(1),',',xi(nx),']'
         call write_log(1, bufout)
      endif
      if (nbefor.gt.0) then
         write(bufout,'(a,i4,a,f10.3,a)') ' WARNING:',nbefor,' positions lie before 1st slice at s_1=', &
                vprf%slc_s(1),', using constant extrapolation'
         call write_log(1, bufout)
      endif
      if (nafter.gt.0) then
         write(bufout,'(a,i4,a,f10.3,a)') ' WARNING:',nafter,' positions lie after last slice at s_n=', &
                vprf%slc_s(vprf%nslc),', using constant extrapolation'
         call write_log(1, bufout)
      endif

      ! clip xi-values outside [s1,sn]

      if (nbefor+nafter.gt.0) then
         do ix = 1, nx
            xi(ix) = max(vprf%slc_s(1), min(vprf%slc_s(vprf%nslc), xi(ix)))
         enddo
      endif

      if (.false.) then
         call write_log(' varprof_intpol_grid: using list version')
         do ix = 1, nx
            xtmp(1:ny) = xi(ix)
            call bspline_get_z_at_xy_list(vprf%spl2d, ny, xtmp, yj, ztmp, my_ierror, defval)
   
            do iy = 1, ny
               ii = ix + (iy-1)*nx
               g_out%z(ii) = ztmp(iy)
            enddo
         enddo
      else
         call write_log(' varprof_intpol_grid: using prod version')
         call bspline_get_z_at_xy_prod(vprf%spl2d, nx, ny, xi, yj, zij, my_ierror, defval)
         g_out%z(1:nx*ny) = reshape(zij, (/ nx*ny /) )
      endif

      ! convert output grid back to track coordinate system

      call cartgrid_2glob(g_out, m_trk)

      deallocate(xi, yj, zij)

   end subroutine varprof_intpol_grid

!------------------------------------------------------------------------------------------------------------

   subroutine varprof_intpol_xunif(vprf, s_ws, nx, x0, dx, g_out, my_ierror)
!--purpose: for a variable profile, compute one or more intermediate slices at requested x_tr positions
!           method 2: using evaluation of 2D tensor B-spline representation
      implicit none
!--subroutine arguments:
      type(t_profile)             :: vprf
      integer,        intent(in)  :: nx
      real(kind=8),   intent(in)  :: x0, dx, s_ws
      type(t_grid)                :: g_out
      integer,        intent(out) :: my_ierror
!--local variables:
      real(kind=8), parameter :: defval = 999d0, tiny = 1d-6
      integer                 :: ix, iy, iy0, iy1, is, ii, ns, nline, nbefor, nafter, sub_ierror
      real(kind=8)            :: xi(nx), xdum2d(1,1)
      type(t_grid)            :: g_trim

      my_ierror = 0

      if (vprf%nslc.le.0) then

         call write_log(' Internal Error: (intpol_xunif2): not a variable profile.')
         call abort_run()

      endif

      ! create output grid with nx slices, ns profile points

      ns = vprf%spl2d%nbrkv
      call grid_create_curvil(g_out, nx, ns, lies_in_oyz=(nx.le.1))

      if (ldebug.ge.1) then
         write(bufout,'(2(a,i4),a)') ' xunif2: output grid has',nx,' x',ns,' points'
         call write_log(1, bufout)
      endif

      ! set uniform x-positions

      do iy = 1, g_out%ny
         do ix = 1, nx
            ii = ix + (iy-1) * g_out%nx
            g_out%x(ii) = x0 + (ix-1) * dx
         enddo
      enddo

      xi(1:nx) = s_ws + g_out%x(1:nx)

      if (ldebug.ge.2) then
         call write_log(' output positions s_ws+xi along track curve:')
         write(bufout,'(10( 20(f10.3),:,/ ))') (xi(ix), ix=1, min(200,nx))
         nline = int( (min(200,nx)-1) / 20 ) + 1
         call write_log(nline, bufout)
      endif
      if (ldebug.ge.2) then
         call write_log(' output positions sj across rail profile:')
         write(bufout,'(10( 20(f10.3),:,/ ))') (vprf%spl2d%vbrk(is), is=1, min(200,ns))
         nline = int( (min(200,ns)-1) / 20 ) + 1
         call write_log(nline, bufout)
      endif

      ! check for xi-positions before first slice or after last slice

      nbefor = 0
      nafter = 0
      do ix = 1, nx
         if (xi(ix).lt.vprf%slc_s(1)) nbefor = nbefor + 1
         if (xi(ix).gt.vprf%slc_s(vprf%nslc)) nafter = nafter + 1
      enddo

      if (ldebug.ge.1) then
         if (nbefor.gt.0) then
            write(bufout,'(a,i4,a,f10.3,a)') ' WARNING:',nbefor,' positions lie before 1st slice at s_1=', &
                   vprf%slc_s(1),', using constant extrapolation'
            call write_log(1, bufout)
         endif
         if (nafter.gt.0) then
            write(bufout,'(a,i4,a,f10.3,a)') ' WARNING:',nafter,' positions lie after last slice at s_n=', &
                   vprf%slc_s(vprf%nslc),', using constant extrapolation'
            call write_log(1, bufout)
         endif
         if (nbefor+nafter.gt.0) then
            write(bufout,'(a,i4,2(a,f10.3),a)') ' varprof_intpol_xunif: requested profile at', nx,         &
                   ' positions s_ws+x in [', xi(1),',',xi(nx),']'
            call write_log(1, bufout)
         endif
      endif

      ! clip xi-values outside [s1,sn]

      if (nbefor+nafter.gt.0) then
         do ix = 1, nx
            xi(ix) = max(vprf%slc_s(1), min(vprf%slc_s(vprf%nslc), xi(ix)))
         enddo
      endif

      ! evaluate 2D spline at given x-values and s-breakpoints

      if (my_ierror.eq.0) then
         ! call write_log(' ...bspline_eval2d_prod')
         ! call bspline_set_debug(5, 1, 1)
         call bspline_eval2d_prod(vprf%spl2d, nx, ns, xi, vprf%spl2d%vbrk, xdum2d, g_out%y, g_out%z,    &
                sub_ierror, defval)
         call bspline_set_debug(0)

         if (ldebug.ge.3) then
            write(bufout,'(2(a,g14.6))') 'g_out%y(1) =', g_out%y(1),', z(1) =',g_out%z(1)
            call write_log(1, bufout)
         endif

         my_ierror = sub_ierror
         if (my_ierror.ne.0) call write_log(' Varprof_intpol_xunif: Error after eval2d')
      endif

      ! re-compute spline for interpolated profile (fast hack, ismooth=0, no kinks)

      if (my_ierror.eq.0 .and. nx.eq.1) then
         iy0 = ns
         iy1 =  1
         do iy = 1, ns
            if (g_out%y(iy).lt.defval-tiny) then
               iy0 = min(iy, iy0)       ! first non-default
               iy1 = max(iy, iy1)       ! last non-default
            endif
         enddo

         call grid_trim(g_out, g_trim, 1, 1, iy0, iy1, with_spline=.false.)
         call grid_copy(g_trim, g_out)
         call grid_destroy(g_trim)

         ! call write_log(' ...grid_make_arclength')
         call grid_make_arclength(g_out, sub_ierror)
         my_ierror = sub_ierror

         ! call grid_print(g_out, 'g_slc', 5)

         ! call write_log(' ...grid_make_ppspline')
         if (ldebug.ge.1) then
            call write_log(' xunif: kchk...')
            call grid_make_ppspline(g_out, vprf%smth, .false., sub_ierror, k_chk=10)
         else
            call grid_make_ppspline(g_out, vprf%smth, .false., sub_ierror)
         endif
         my_ierror = sub_ierror
      endif

   end subroutine varprof_intpol_xunif

!------------------------------------------------------------------------------------------------------------

   subroutine varprof_find_slices_xrange(vprf, s_ws, xmin, xmax, islc0, islc1)
!--purpose: determine the slices [islc0, islc1] needed for interpolation to s_ws+[xmin, xmax]
      implicit none
!--subroutine arguments:
      type(t_profile)           :: vprf
      real(kind=8), intent(in)  :: s_ws, xmin, xmax
      integer,      intent(out) :: islc0, islc1
!--local variables:
      real(kind=8)       :: smin, smax

      smin = s_ws + xmin
      smax = s_ws + xmax

      ! no interpolation is needed if all si lie before start or after end of variable profile

      if (smax.le.vprf%slc_s(1)) then

         ! whole grid before start of variable profile

         if (ldebug.ge.1) then
            write(bufout,'(3(a,f12.3),a)') ' grid has s_tr \in [',smin,',',smax,                        &
                           '], before start of varprof, s1=', vprf%slc_s(1)
            call write_log(1, bufout)
         endif

         islc0 = 1
         islc1 = 1

      elseif (smin.ge.vprf%slc_s(vprf%nslc)) then

         ! whole grid after end of variable profile

         if (ldebug.ge.1) then
            write(bufout,'(3(a,f12.3),a)') ' grid has s_tr \in [',smin,',',smax,                        &
                        '], after end of varprof, sn=', vprf%slc_s(vprf%nslc)
            call write_log(1, bufout)
         endif

         islc0 = vprf%nslc
         islc1 = vprf%nslc

      else

         ! determine the range of 'relevant slices' islc0 to islc1

         islc0 = vprf%nslc
         do while( vprf%slc_s(islc0).gt.smin .and. islc0.gt.1 )
            islc0 = islc0 - 1
         enddo
         islc1 = 1
         do while( vprf%slc_s(islc1).le.smax .and. islc1.lt.vprf%nslc )
            islc1 = islc1 + 1
         enddo

      endif ! [smin,smax] before start / after end

      if (ldebug.ge.2) then
         write(bufout,'(2(a,i5),4(a,f11.3),a)') ' using slices',islc0,':',islc1,' with s_slc = [',      &
                  vprf%slc_s(islc0),',',vprf%slc_s(islc1),'] for grid s_tr \in [',smin,',',smax,']'
         call write_log(1, bufout)
      endif

   end subroutine varprof_find_slices_xrange

!------------------------------------------------------------------------------------------------------------

   subroutine varprof_available_x(vprf, s_ws, ny, yr, xmax, x_lbnd, x_ubnd)
!--purpose: for each position yr, set range of x-values where profile is available
      implicit none
!--subroutine arguments:
      type(t_profile)           :: vprf
      integer,      intent(in)  :: ny
      real(kind=8), intent(in)  :: s_ws, xmax, yr(ny)
      real(kind=8), intent(out) :: x_lbnd(ny), x_ubnd(ny)
!--local variables:
      integer            :: islc, islc0, islc1, iy
      real(kind=8)       :: xslc, ymin, ymax

      ! determine the slices [islc0, islc1] needed for interpolation to s_ws+[-xmax, xmax]

      call varprof_find_slices_xrange(vprf, s_ws, -xmax, xmax, islc0, islc1)

      ! initialize bracket to 'no data at all'

      do iy = 1, ny
         x_lbnd(iy) =  1d9
         x_ubnd(iy) = -1d9
      enddo

      ! apply availability for available slices
      !  - assuming no vertical slopes within each slice, using y-range = [y(1), y(n)]
      !  - assuming continuous availability from x(i0) to x(i1) if available at i0 and i1

      do islc = islc0, islc1
         associate(gslc => vprf%slc_grd(islc)%g)
         xslc = vprf%slc_s(islc) - s_ws
         ymin = gslc%y(1)
         ymax = gslc%y(gslc%ntot)
         do iy = 1, ny
            if (yr(iy).ge.ymin .and. yr(iy).le.ymax) then
               x_lbnd(iy) = min(x_lbnd(iy), xslc)
               x_ubnd(iy) = max(x_ubnd(iy), xslc)
            endif
         enddo

         if (ldebug.ge.5) then
            write(bufout,'(a,i3,4(a,f12.3),a)') ' slice',islc,': s=',vprf%slc_s(islc),', xslc=',xslc,   &
                   ': y = [',ymin,',',ymax,']'
            call write_log(1, bufout)
         endif
         end associate
      enddo

      ! clip at input [-xmax, xmax]

      do iy = 1, ny
         x_lbnd(iy) = max(-xmax, x_lbnd(iy))
         x_ubnd(iy) = min( xmax, x_ubnd(iy))
      enddo

   end subroutine varprof_available_x

!------------------------------------------------------------------------------------------------------------

end module m_varprof
