!------------------------------------------------------------------------------------------------------------
! m_licensing -  implementation of license checking functionality - dummy version
!
! Copyright 2016-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_licensing
use, intrinsic        :: iso_c_binding, only: C_INT, C_SHORT, C_CHAR, C_DOUBLE, C_BOOL, C_NULL_CHAR
use m_print_output
implicit none
public

   integer            :: lic_debug  = 1
                                ! 0 = silent; 1 = user feedback; 2 = extended user feedback; >=3 = internal

   !---------------------------------------------------------------------------------------------------------

   ! - Version number for current version.  For the trunk: version number for next release
   integer, parameter :: my_version_major  = 23
   integer, parameter :: my_version_minor  =  1
   integer, parameter :: my_version_update =  1

   !---------------------------------------------------------------------------------------------------------

   ! information contained in a CONTACT license

   type :: t_license
      ! licensee information:
      character(len=64) :: licensee_name
      ! licensed product:
      integer           :: max_ver_major, max_ver_minor
      ! license status:
      character(len=64) :: license_error
      logical           :: is_valid
   end type t_license

   ! global instance of license structure

   type(t_license)      :: my_license

   ! routine to display the startup-message and license information

   public  WriteVersion
   public  LicenseGetDebug
   public  LicenseSetDebug

   ! generic helper routines

   public  abort_run            ! clean-up and stop program execution
   public  GetEnvironmentVar    ! get environment variable with error reporting
   public  DirIsWritable        ! check that folder exists and is writable
   public  ScanFolder           ! find files using a pattern

contains

!------------------------------------------------------------------------------------------------------------

subroutine WriteVersion(num_version, version, license)
!--purpose: write version strings and license info to out-file/screen/..
   implicit none
!--subroutine arguments:
   integer            :: num_version
   character(len=*)   :: version(num_version)
   type(t_license)    :: license
!--local variables:
   integer            :: iv, len1
   character(len=80)  :: spaces

   ! Display the contents of string array "version"

   write(bufout,123)
   call write_log(1, bufout)
   do iv = 1, num_version
      ! skip empty lines version(iv), except version(3)
      if (iv.eq.3 .or. version(iv).ne.' ') then
         write(bufout,124) version(iv)
         call write_log(1, bufout)
      endif
   enddo

   if (license%is_valid) then
      spaces = ' '
      len1 = len(c_trim(license%licensee_name))

      if (license%licensee_name(1:12).eq.'Open source version') then

         write(bufout,125) c_trim(license%licensee_name), ' ',' ',' ', spaces(1:65-len1)
         call write_log(1, bufout)

      endif
   endif

   ! Display error message if no valid license was found

   if (.not.license%is_valid) then
      write(bufout,124)  c_trim(my_license%license_error)
      call write_log(1, bufout)
   endif

   write(bufout,123)
   call write_log(1, bufout)

123 format(3x,'------------------------------------------------------------------------')
124 format(3x,'|  ',a66,'  |')
125 format(3x,'|  ',5a,   '|')

end subroutine WriteVersion

!------------------------------------------------------------------------------------------------------------

function LicenseGetDebug
!--function: get level of debug output of license management routines
   implicit none
!--return value:
   integer         :: LicenseGetDebug

   LicenseGetDebug = lic_debug

end function LicenseGetDebug

!------------------------------------------------------------------------------------------------------------

subroutine LicenseSetDebug(new_lic_debug)
!--function: enable/disable debug output of license management routines
   implicit none
!--subroutine arguments:
   integer, intent(in)           :: new_lic_debug       ! level of debug output required

   lic_debug = new_lic_debug

   if (lic_debug.ge.3) then
      write(bufout,'(a,i3,2(a,i7))') ' license: debugging level =',lic_debug
      call write_log(1, bufout)
   endif
end subroutine LicenseSetDebug

!------------------------------------------------------------------------------------------------------------

subroutine abort_run
!--purpose: clean up and abort
   implicit none

   ! End program execution

   stop
end subroutine abort_run

!------------------------------------------------------------------------------------------------------------

subroutine GetEnvironmentVar(envnam, envvar, istat, idebug)
!--purpose: helper that gets environment variable and prints error messages
   implicit none
!--subroutine arguments:
   integer,          intent(in)  :: idebug
   character(len=*), intent(in)  :: envnam
   character(len=*), intent(out) :: envvar
   integer,          intent(out) :: istat
!--local variables:
   integer                 :: lenenv

   call get_environment_variable(envnam, envvar, length=lenenv, status=istat, trim_name=.true.)

   ! in case of failure, write error message

   if (istat.eq.-1) then

      envvar = ' '
      if (idebug.ge.1) call write_log(' Environment variable ' // trim(envnam) // ' too long for lic.path')

   elseif (istat.ge.2) then

      envvar = ' '
      if (idebug.ge.1) then
         write(bufout,'(3a,i5,a,i5)') ' Environment variable ', trim(envnam),': obtained status=',      &
                istat,', length=',lenenv
         call write_log(1, bufout)
      endif

   elseif (istat.eq.1) then

      if (idebug.ge.2) then
         write(bufout,'(3a)') ' Environment variable ', trim(envnam),' does not exist'
         call write_log(1, bufout)
      endif

   endif

end subroutine GetEnvironmentVar

!------------------------------------------------------------------------------------------------------------

function DirIsWritable(folder)
!--purpose: check that folder exists and is writable
   implicit none
!--subroutine arguments:
   character(len=*)     :: folder
!--local variables:
   integer, parameter   :: idebug = 0
   integer              :: irnd, ios
   real(kind=8)         :: rnd
   character(len=256)   :: test_file
!--return value:
   logical              :: DirIsWritable

   call random_seed()
   call random_number(rnd)
   irnd = nint(100000d0 * rnd)

   write(test_file,'(a,i6.6,a)') 'this_is_a_test_', irnd, '.txt'
   test_file = trim(folder) // path_sep // trim(test_file)

   open (ltmp, file=test_file, action='write', status='new', dispose='delete', iostat=ios)

   if (idebug.ge.1) then
      if (ios.eq.9) then
         write(bufout,'(3a,i3,a)') ' ERROR: permission denied on file "', trim(test_file), '" (ios=', ios,')'
         call write_log(1, bufout)
      elseif (ios.eq.10) then
         write(bufout,'(3a,i3,a)') ' ERROR: file exists "', trim(test_file), '" (ios=', ios,')'
         call write_log(1, bufout)
      elseif (ios.eq.43) then
         write(bufout,'(3a,i3,a)') ' ERROR: invalid filename "', trim(test_file), '" (ios=', ios,')'
         call write_log(1, bufout)
      elseif (ios.ne.0) then
         write(bufout,'(3a,i6,a)') ' ERROR: cannot open file "', trim(test_file), '" for writing (ios=', &
                ios,')'
         call write_log(1, bufout)
      endif
   endif

   if (idebug.ge.2 .and. ios.eq.0) then
      write(bufout,'(3a,i6,a)') ' Test file "', trim(test_file), '" opened ok.'
      call write_log(1, bufout)
   endif

   ! close file & discard (dispose=delete)

   if (ios.eq.0) close(ltmp)

   ! return success or failure

   DirIsWritable = (ios.eq.0)

end function DirIsWritable

!------------------------------------------------------------------------------------------------------------

subroutine ScanFolder(folder, pattern, nfound, fullname, idebug)
!--purpose: scan one folder for filename matching the pattern, return last one in fullname
   use ifport
   implicit none
!--subroutine arguments:
   character(len=*), intent(in)  :: folder, pattern     ! pattern: C-string
   character(len=*), intent(out) :: fullname
   integer,          intent(in)  :: idebug
   integer,          intent(out) :: nfound
!--local variables:
   character(240)                :: fullpatt
   type(file$info)               :: match
   integer(kind=int_ptr_kind())  :: handle
   integer(kind=4)               :: result
 
   if (idebug.ge.2) then
      call write_log(' Searching files matching "'//c_trim(pattern)//'" in folder "'//trim(folder)//'".')
   endif

   fullpatt = c_trim(folder) // path_sep // c_trim(pattern)
   fullname = ' '
   nfound   = 0
   result   = 1
  
   handle   = file$first
   do while(result.gt.0 .and. handle.ne.file$last .and. handle.ne.file$error)
      ! (Linux:) returns file$last if no files match at all, file$error if no more files could be found

      result = getfileinfoqq(fullpatt, match, handle)

      ! write(bufout,*) 'result=',result,', handle=',handle,' (first=',file$first,', last=',file$last,', error=',file$error,')'
      ! call write_log(1, bufout)

      if (result.gt.0) then
         nfound   = nfound + 1
         fullname = trim(folder) // path_sep // trim(match%name)
         ! fullname = trim(match%name)
         if (idebug.ge.3) call write_log(' ScanFolder: found file "' // trim(match%name) // '".')
      endif

      if (handle.eq.file$last .or. handle.eq.file$error) then
         if (idebug.ge.3 .and. nfound.ge.1) call write_log(' ScanFolder: no more files.')
      endif
   enddo

   if (nfound.le.0) then
      if (idebug.ge.2) call write_log(' ScanFolder: no files matching "' // trim(fullpatt) // '".')
   endif

end subroutine ScanFolder

!------------------------------------------------------------------------------------------------------------

end module m_licensing
