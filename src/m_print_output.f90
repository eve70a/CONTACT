!------------------------------------------------------------------------------------------------------------
! m_print_output - This module governs the output of CONTACT to the screen or log-file(s).
!
! The idea is that eventually the program does not contain print-statements anymore,
! only uses the functions "write_log" provided here.
!
! By design there are multiple "streams" for the print-output that can be activated and turned off.
! Currently two streams are implemented:
!    "lout" == the .out-file, note that this is used for logging and for "true" output as well.
!    "*"    == the screen
! Further, the module declares the unit numbers for in- and output-files, such as the .out-file
! where the results are stored.
!
! Copyright 2008-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_print_output
use, intrinsic   :: iso_c_binding, only: C_CHAR, C_NULL_CHAR
#ifdef _OPENMP
   use omp_lib    , only :  omp_in_parallel
#endif
implicit none
public

   ! fmtmat   version number of (ASCII) matlab-file (100*major + 10*minor + 1*patch)
   ! linp     logical unit number of .inp file
   ! lout     logical unit number of .out file
   ! lsbs     logical unit number of .subs file for subsurface stresses
   ! ltmp*    logical units for temporary use, see get_lunit_tmp_use

   integer,      parameter :: fmtmat =   2410
   integer                 :: linp   =     11
   integer                 :: lout   =     12
   integer                 :: ltmp1  =     17
   integer                 :: ltmp2  =     18
   integer                 :: ltmp3  =     19
   public  fmtmat, linp, lout, ltmp1, ltmp2, ltmp3

   ! inp_open  flag indicating the status of the .inp-file
   ! out_open  flag indicating the status of the .out-file
   !              -1 == could not open;
   !               0 == uninitialized;
   !               1 == open for output (write_inp);
   !               2 == open for input (read_inp);
   ! ltmp_open flag indicating the status of logical units for temporary use
   !               0 == free
   !               3 == in use (get_lunit_tmp_use, free_lunit_tmp_use).

   integer                 :: inp_open = 0
   integer                 :: out_open = 0
   integer                 :: ltmp1_open = 0
   integer                 :: ltmp2_open = 0
   integer                 :: ltmp3_open = 0
   public  inp_open, out_open, ltmp1_open, ltmp2_open, ltmp3_open

   public get_lunit_tmp_use
   public free_lunit_tmp_use

   ! MAX_CHAR_INP  length of character variables for reading input, max. length of input lines

   integer,      parameter :: MAX_CHAR_INP = 256
   public  MAX_CHAR_INP

   ! platform  integer indicating the platform used
   ! path_sep  direction/file separation character

   integer,          parameter :: plat_lnx = 1991
   integer,          parameter :: plat_win = 1983

   character(len=1), parameter :: path_sep_win = '\\'
   character(len=1), parameter :: path_sep_lnx = '/'
#if defined _WIN32 || defined _WIN64
   integer,          parameter :: platform = plat_win
   character(len=1), parameter :: path_sep = path_sep_win
#else
   integer,          parameter :: platform = plat_lnx
   character(len=1), parameter :: path_sep = path_sep_lnx
#endif

   ! helper routines to find first or last '/' or '\' in a filename, or to set for current platform

   public  index_pathsep
   public  set_platform_filesep
   public  is_absolute_path
   public  make_absolute_path

   ! routines for activating and de-activating the different output streams

   public set_print_streams
   public get_print_streams
   public save_print_streams
   public restore_print_streams

   ! routine(s) used for printing information to the output streams

   public write_log
   interface write_log
      module procedure write_log_str
      module procedure write_log_msg
   end interface write_log

   ! bufout   convenience: provide a generic buffer for preparing output to log-streams, instead of
   !          having to declare this everywhere.

   character(len=240)   :: bufout(10)
   public                  bufout
!$omp threadprivate      ( bufout )

   ! prn_thread thread number to be printed along with the output " 1 #", <0: suppress prefix

   integer              :: prn_thread  = -1
   character(len=4)     :: pfx_thread = ' '
!$omp threadprivate      ( prn_thread, pfx_thread )

   ! routine for configuring the thread number (prefix)

   public set_print_thread

   ! Internal variables:
   ! loutfil   flag indicating that writing to the out-file is active
   ! lscreen   flag indicating that writing to the screen is active
   ! lsimpck   flag indicating that writing to the SIMPACK log-file is active

   logical              :: loutfil = .false.
   logical              :: lscreen = .false.
   logical              :: lsimpck = .false.

   ! savd_outfil, savd_screen, savd_simpck, ...: allow memorizing the flags loutfil, lscreen, lsimpck
   !           for temporary re-routing the output

   logical              :: savd_outfil = .false.
   logical              :: savd_screen = .false.
   logical              :: savd_simpck = .false.

#ifdef _WIN32
   !dir$ attributes stdcall :: cntc_print
#endif

   ! generic helper routines: filter small values, format numbers without leading zero

   public  fmt_gs
   public  filt_sml

   ! generic helper routines: convert strings to upper or lower case, length/trim C-strings

   public  to_upper
   public  to_lower
   public  c_len
   public  c_trim
   public  c_to_f_string

contains

!------------------------------------------------------------------------------------------------------------

pure function index_pathsep(fname, back)
!--purpose: determine the index of the first or last path_sep character in a filename
   implicit none
!--function argument:
   character(len=*),           intent(in) :: fname
   logical,          optional, intent(in) :: back
!--local variables:
   integer        :: ix
!--return value:
   integer        :: index_pathsep

   ! first look for the path_sep ('/' or '\') for the current platform

   ix = index(fname, path_sep, back)

   ! if not found, look for the path_sep of the alternative platform

   if (ix.le.0) then
      if (platform.eq.plat_win) then
         ix = index(fname, path_sep_lnx, back)
      else
         ix = index(fname, path_sep_win, back)
      endif
   endif
   index_pathsep = ix

end function index_pathsep
 
!------------------------------------------------------------------------------------------------------------

subroutine set_platform_filesep(fname)
!--purpose: replace path_sep characters '/' and '\' in a filename by path_sep for current platform
   implicit none
!--subroutine argument:
   character(len=*), intent(inout) :: fname
!--local variables:
   integer          :: ix
   character(len=1) :: my_sep, alt_sep

   if (platform.eq.plat_win) then
      my_sep  = path_sep_win
      alt_sep = path_sep_lnx
   else
      my_sep  = path_sep_lnx
      alt_sep = path_sep_win
   endif

   ! look for the path_sep ('/' or '\') for the other platform

   ix = index(fname, alt_sep)

   do while(ix.ge.1)

      ! if found, replace and look again

      fname(ix:ix) = my_sep
      ix = index(fname, alt_sep)
   enddo

end subroutine set_platform_filesep
 
!------------------------------------------------------------------------------------------------------------

pure function is_absolute_path(fname)
!--purpose: determine whether a filename is absolute or relative
   implicit none
!--function argument:
   character(len=*),           intent(in) :: fname
!--return value:
   logical        :: is_absolute_path

   ! absolute path if fname starts with '/' or '\' 

   if (index_pathsep(fname).eq.1) then

      is_absolute_path = .true.

   elseif (platform.eq.plat_win .and. len(fname).ge.3) then

      ! Windows: absolute path if fname starts with drive e.g. 'C:\'

      is_absolute_path = (fname(2:2).eq.':' .and. index_pathsep(fname).eq.3)

   else

      is_absolute_path = .false.

   endif

end function is_absolute_path

!------------------------------------------------------------------------------------------------------------

subroutine make_absolute_path(fname, dirnam, fulnam)
!--purpose: in case of relative path-name, form absolute path-name
!           Note: fulnam may refer to fname for in-place updating
   implicit none
!--subroutine arguments
   character(len=*)  :: fname, dirnam
   character(len=*)  :: fulnam          ! Note: fulnam can be equal to fname
!--local variables

   if (dirnam.ne.' ' .and. .not.is_absolute_path(fname)) then
      fulnam = trim(dirnam) // path_sep // fname
   else
      fulnam = fname
   endif
   call set_platform_filesep(fulnam)

end subroutine make_absolute_path

!------------------------------------------------------------------------------------------------------------

subroutine set_print_streams(nw_outfil, nw_screen, nw_simpck)
!--purpose: set activity of the output streams
   implicit none
!--subroutine arguments
   logical, intent(in) :: nw_outfil, nw_screen, nw_simpck
!--local variables

   loutfil = nw_outfil
   lscreen = nw_screen
   lsimpck = nw_simpck
end subroutine set_print_streams

!------------------------------------------------------------------------------------------------------------

subroutine get_print_streams(cur_outfil, cur_screen, cur_simpck)
!--purpose: get activity status of the output streams
   implicit none
!--subroutine arguments
   logical, intent(out) :: cur_outfil, cur_screen, cur_simpck
!--local variables

   cur_outfil = loutfil
   cur_screen = lscreen 
   cur_simpck = lsimpck
end subroutine get_print_streams

!------------------------------------------------------------------------------------------------------------

subroutine save_print_streams
!--purpose: save the current activity status of the output streams
   implicit none
!--subroutine arguments
!--local variables

   savd_outfil = loutfil
   savd_screen = lscreen 
   savd_simpck = lsimpck
end subroutine save_print_streams

!------------------------------------------------------------------------------------------------------------

subroutine restore_print_streams
!--purpose: restore the last saved activity status of the output streams
   implicit none
!--subroutine arguments
!--local variables

   loutfil = savd_outfil
   lscreen = savd_screen
   lsimpck = savd_simpck
end subroutine restore_print_streams

!------------------------------------------------------------------------------------------------------------

subroutine set_print_thread(prn_thread_arg)
!--purpose: set the thread number to be printed as prefix to the output " 1 #"
   implicit none
!--subroutine arguments
   integer, intent(in) :: prn_thread_arg
!--local variables

   prn_thread = prn_thread_arg
   write(pfx_thread,'(i2,a)') prn_thread, ' #'
end subroutine set_print_thread

!------------------------------------------------------------------------------------------------------------

subroutine write_log_str(str)
!--purpose: Write a single string to the log output streams (screen/.out-file/SIMPACK).
   implicit none
!--subroutine arguments
   character(len=*), intent(in) :: str
!--local variables
   integer :: iparll

#ifdef _OPENMP
   iparll = omp_in_parallel()
#else
   iparll = 0
#endif

   ! write message to .out-file, if this is active:

   if (loutfil .and. out_open.eq.1) then
!$omp critical (writelog_str_outfil)
      if (prn_thread.ge.1 .or. iparll.ne.0) then
         write(lout, '(2a)') pfx_thread, trim(str)
      else
         write(lout, '(a)') trim(str)
      endif
!$omp end critical (writelog_str_outfil)
   endif

   ! write message to the screen, if this is active:

   if (lscreen) then
!$omp critical (writelog_str_screen) 
      if (prn_thread.ge.1 .or. iparll.ne.0) then
         write(lout, '(2a)') pfx_thread, trim(str)
      else
         write(*, '(a)') trim(str)
      endif
!$omp end critical (writelog_str_screen) 
   endif

end subroutine write_log_str

!------------------------------------------------------------------------------------------------------------

subroutine write_log_msg(nlines, msg)
!--purpose: Write an array of strings to the log output streams (screen/.out-file/SIMPACK)
   implicit none
!--subroutine arguments
   integer,          intent(in) :: nlines
   character(len=*), intent(in) :: msg(nlines)
!--local variables
   integer   :: j, iparll

#ifdef _OPENMP
   iparll = omp_in_parallel()
#else
   iparll = 0
#endif

   ! write message to .out-file, if this is active:

   if (loutfil .and. out_open.eq.1) then
!$omp critical (writelog_msg_outfil)
      if (prn_thread.ge.1 .or. iparll.ne.0) then
         write(lout, '(2a)') (pfx_thread, trim(msg(j)), j=1,nlines)
      else
         write(lout, '(a)') (trim(msg(j)), j=1,nlines)
      endif
!$omp end critical (writelog_msg_outfil)
   endif

   ! write message to the screen, if this is active:

   if (lscreen) then
!$omp critical (writelog_msg_screen)
      if (prn_thread.ge.1 .or. iparll.ne.0) then
         write(*, '(2a)') (pfx_thread, trim(msg(j)), j=1,nlines)
      else
         write(*, '(a)') (trim(msg(j)), j=1,nlines)
      endif
!$omp end critical (writelog_msg_screen)
   endif

end subroutine write_log_msg

!------------------------------------------------------------------------------------------------------------

function get_lunit_tmp_use()
!--purpose: claim one of the logical unit numbers available for temporary use
   implicit none
!--return value:
   integer             :: get_lunit_tmp_use

   if     (ltmp1.gt.0 .and. ltmp1_open.eq.0) then
      ltmp1_open = 3
      get_lunit_tmp_use = ltmp1
   elseif (ltmp2.gt.0 .and. ltmp2_open.eq.0) then
      ltmp2_open = 3
      get_lunit_tmp_use = ltmp2
   elseif (ltmp3.gt.0 .and. ltmp3_open.eq.0) then
      ltmp3_open = 3
      get_lunit_tmp_use = ltmp3
   else
      get_lunit_tmp_use = -1
      ! call write_log(' ERROR: no free lunit')
      ! stop
   endif

end function get_lunit_tmp_use

!------------------------------------------------------------------------------------------------------------

subroutine free_lunit_tmp_use(ltmp)
!--purpose: release logical unit number claimed previously for temporary use
   implicit none
!--subroutine arguments:
   integer             :: ltmp

   if (ltmp.eq.ltmp1) then
      ltmp1_open = 0
   elseif (ltmp.eq.ltmp2) then
      ltmp2_open = 0
   elseif (ltmp.eq.ltmp3) then
      ltmp3_open = 0
   else
      ! call write_log(' ERROR: incorrect ltmp')
      ! stop
   endif

end subroutine free_lunit_tmp_use

!------------------------------------------------------------------------------------------------------------

function fmt_gs(iw, idf, ides, value)
!--purpose: Print a value with "GS"-format descriptor: switch between Fiw.idf and ESiw.ides depending
!           on size of value. Output string of iw characters with idf/ides significant digits. 
!           examples: fmt_gs(12, 6, 6, -12345678.9)   --> '-1.23457E+07'
!                     fmt_gs(12, 8, 6, -12345678.9)   --> '  -12345679.'
!                     fmt_gs(12, 6, 6, -123.456789)   --> '    -123.457'
!                     fmt_gs(12, 6, 6, -0.123456789)  --> '   -0.123457'
!                     fmt_gs(12, 6, 6, -0.0123456789) --> '-1.23457E-02'
!           F iw.idf: 3 characters for minus + leading zero + dot, so idef <= iw-3
!           ES iw.ides: 6 characters for minus + dot + exponent, so ides <= iw-6.
   implicit none
!--subroutine arguments:
   integer,          intent(in)  :: iw, idf, ides
   real(kind=8),     intent(in)  :: value
!--return value:
   character(len=iw)  :: fmt_gs
!--local variables :
   integer, parameter :: ie = 2  ! #digits used for exponent e.g. E+07, E-02
   integer            :: lg      ! 10-log of input value
   character(len=12)  :: fmtstr
   character(len=iw)  :: strng

   if (abs(value).lt.1d-20) then
      lg = 0
   else
      lg = floor( log(abs(value)) / log(10d0) + 1d-7) 
   endif
   if (lg.ge.idf .or. lg.le.-2) then
      write(fmtstr,'(3(a,i2),a)') '(ES',iw,'.',ides-1,'E',ie,')'
      write(strng,fmtstr) value
   else
      write(fmtstr,'(2(a,i2),a)') '(F',iw,'.',max(0,idf-lg-1),')' ! idf-lg-1 digits for fraction
      write(strng,fmtstr) value
   endif
   fmt_gs = strng
end function fmt_gs

!------------------------------------------------------------------------------------------------------------

function filt_sml(rvalue, thresh)
!--function: filter values that are close to zero
   implicit none
!--result value
   real(kind=8) :: filt_sml
!--subroutine arguments:
   real(kind=8),     intent(in)  :: rvalue, thresh

   if (abs(rvalue).lt.thresh) then
      filt_sml = 0d0
   else
      filt_sml = rvalue
   endif

end function filt_sml

!------------------------------------------------------------------------------------------------------------

pure function to_upper (in_str) result (out_str)
!--purpose: change a string to upper case
   implicit none
!--function argument:
   character(*), intent(in) :: in_str
!--return value:
   character(len(in_str))   :: out_str
!--local variables:
   integer :: il, is
   character(len=26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   character(len=26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

   ! Capitalize each letter if it is lowercase

   out_str = in_str
   do is = 1, len_trim(in_str)
       il = index(low, in_str(is:is))
       if (il.gt.0) out_str(is:is) = cap(il:il)
   end do

end function to_upper

!------------------------------------------------------------------------------------------------------------

pure function to_lower (in_str) result (out_str)
!--purpose: change a string to lower case
   implicit none
!--function argument:
   character(*), intent(in) :: in_str
!--return value:
   character(len(in_str))   :: out_str
!--local variables:
   integer :: ic, is
   character(len=26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   character(len=26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

   ! De-capitalize each letter if it is upper case

   out_str = in_str
   do is = 1, len_trim(in_str)
       ic = index(cap, in_str(is:is))
       if (ic.gt.0) out_str(is:is) = low(ic:ic)
   end do

end function to_lower

!------------------------------------------------------------------------------------------------------------

pure function c_len(f_str)
!--purpose: get length of string that may be C-terminated or tab/blank-terminated
   implicit none
!--function argument:
   character(len=*), intent(in) :: f_str
!--return value:
   integer                      :: c_len

   c_len = scan(f_str, C_NULL_CHAR, back=.false.) - 1
   if (c_len.ne.-1) then
      c_len = len_trim(f_str(:c_len)) ! C string: ignore NULL and everything thereafter
   else
      c_len = len_trim(f_str)         ! else: return length of trimmed string
   end if
end function c_len
 
!------------------------------------------------------------------------------------------------------------

pure function c_trim(in_str) result (out_str)
!--purpose: trim a string that may be C-terminated or tab/blank-terminated, removing any null characters
   implicit none
!--function argument:
   character(len=*), intent(in) :: in_str
!--return value:
   character(len= c_len(in_str//"") ) :: out_str
!--local variables:
 
   out_str = in_str

   ! to keep a null-character if present in input:
   ! character,    parameter  :: null_ch = char(0)
   ! character(len= c_len(in_str//"") + min(1, scan(in_str,null_ch)) ) :: out_str
   ! if (scan(in_str,null_ch).ne.0) out_str = out_str(:c_len(out_str)) // ""C
   ! ""C looks like old-fashioned way to say C_NULL_CHAR?
end function c_trim

!------------------------------------------------------------------------------------------------------------

subroutine c_to_f_string(c_str, f_str, len_c_arg)
!--purpose: copy a C-string argument to a Fortran character variable
!           c_str may be NULL-terminated or may have length len_c_arg
   implicit none
!--subroutine arguments:
   character(len=1,kind=C_CHAR), intent(in)  :: c_str(*)
   character(len=*),             intent(out) :: f_str
   integer,            optional, intent(in)  :: len_c_arg
!--local variables:
   integer      :: i, len_f, len_c

   len_c = 999999
   if (present(len_c_arg)) len_c = len_c_arg

   len_f = min(len(f_str), len_c)
   f_str = ' '

   i = 1
   do while(i.le.len_f .and. c_str(i).ne.C_NULL_CHAR)
      f_str(i:i) = c_str(i)
      i = i + 1
   enddo

end subroutine c_to_f_string

!------------------------------------------------------------------------------------------------------------

end module m_print_output
