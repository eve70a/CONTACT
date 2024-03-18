!------------------------------------------------------------------------------------------------------------
! m_readline - functions for managing the input-file
!
! Copyright 2008-2023 by Vtech CMCC.
!
! Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.
!------------------------------------------------------------------------------------------------------------
module m_readline

use m_globals

implicit none
private

public  readline
public  read1darr
public  readkeywordvalues
public  readplainvalues
public  getnonemptyline
private splitlinetowords

! functions:
public check_irng
public check_2rng
public check_2int
public check_3rng
public check_range
public check_equal
public check_smaller
public check_sorted

contains

!------------------------------------------------------------------------------------------------------------

   subroutine readline(unitnm, ncase, linenr, descrip, types, ints, dbles, flags, strngs, mxnval,       &
                       nval, idebug, ieof_arg, lstop, ierror)
!--purpose: Read a line from unit unitnm, retrieve values as indicated by types (i/I == integer,
!           d/D == real(kind=8), a/A == angle(kind=8), l/L == logical, s/S == string; lower case
!           mandatory, upper case: optional) and return in ints, dbles, flags, strngs.
!           Write error message with linenr and descrip in case of errors in the user input.
!           Return ieof=1 for end-of-file, treated as error when ieof=-1 on input.
      implicit none
!--subroutine arguments:
      integer,       intent(in)    :: unitnm, ncase, mxnval, idebug
      integer,       intent(inout) :: linenr, ieof_arg
      integer,       intent(out)   :: ints(mxnval), nval, ierror
      logical,       intent(in)    :: lstop
      logical,       intent(out)   :: flags(mxnval)
      real(kind=8),  intent(out)   :: dbles(mxnval)
      character(len=*), intent(in)    :: descrip, types
      character(len=*), intent(out)   :: strngs(mxnval)
!--local variables:
      integer,      parameter :: mxitem = 20
      real(kind=8), parameter :: pi     = 4d0*atan(1d0)
      integer            :: nblank, ncmtln, nval_tot, nval_req, nval_opt, ieof, iint, idbl, iflg,       &
                            istr, ilen, nitem, ii, iostat
      logical            :: is_deg
      character(len=1)   :: squote = ''''
      character(len=MAX_CHAR_INP) :: inptxt, words(1:mxitem)

      if (idebug.ge.5) then
         write(bufout,'(/,2a)') '--- readLine: ',trim(descrip)
         call write_log(2, bufout)
      endif

      ieof   = ieof_arg  ! -1: end-of-file is considered an error, 0: no error for eof
      ierror = 0

      ! initialise required and optional number of values to be read

      nval_tot = len(trim(types))
      nval_req = 0
      nval_opt = 0
      do ii = 1, nval_tot
         if (types(ii:ii).eq.'i' .or. types(ii:ii).eq.'d' .or. types(ii:ii).eq.'a' .or.                 &
             types(ii:ii).eq.'l' .or. types(ii:ii).eq.'s') then
            nval_req = nval_req + 1
            if (nval_opt.gt.0) then
               call write_log(' Internal ERROR: mandatory value after optional value(s): ' // trim(types))
               call abort_run()
            endif
         elseif (types(ii:ii).eq.'I' .or. types(ii:ii).eq.'D' .or. types(ii:ii).eq.'A' .or.             &
             types(ii:ii).eq.'L' .or. types(ii:ii).eq.'S') then
            nval_opt = nval_opt + 1
         else
            call write_log(' Internal ERROR: invalid data-type ' // types(ii:ii))
            call abort_run()
         endif
      enddo

      if (nval_req.le.0) then
         call write_log(' Internal ERROR: no mandatory value(s): ' // trim(types))
         call abort_run()
      endif

      ! initialise counters

      nval  = 0
      iint  = 0
      idbl  = 0
      iflg  = 0
      istr  = 0

      ! read new lines as long as not all required values have been filled

      do while (nval.lt.nval_req .and. ieof.le.0)

         ! increment line number and read input, thereby skipping comments and empty lines

         call getnonemptyline(unitnm, ncase, descrip, '%', idebug, ieof, lstop, linenr, nblank,         &
                ncmtln, inptxt)

         if (ieof.gt.0) then

            ! reached end of file; considered an error when input ieof_arg < 0

            if (ieof_arg.lt.0) then
               if (lstop) call abort_run()
               ierror = -1
            endif

         else

            ! split the input at spaces and put into string array

            call splitlinetowords(inptxt, '%', mxitem, words, nitem)

            ! put content of string array in ints, dbles and flags

            ! while (more (req/opt) values needed & more words available) do

            ii     = 0
            iostat = 0
            do while (nval.lt.nval_tot .and. ii.lt.nitem .and. iostat.eq.0)

               ii   = ii + 1
               nval = nval + 1

               if (to_lower(types(nval:nval)).eq.'i') then

                  ! obtain an integer value from string words(ii)

                  iint = iint + 1
                  read(words(ii), *, IOSTAT=iostat) ints(iint)
                  if (idebug.ge.9) then
                     write(bufout,'(2(a,i3),2(a,i7),a,i3,3a)') 'item', ii,                              &
                        ', type 1 -> ints(', iint, ') =',ints(iint),' ; iostat =', iostat,              &
                        ' ; words(', ii, ') = "', trim(words(ii)), '"'
                     call write_log(1, bufout)
                  endif
                  if (iostat.ne.0) iint = iint - 1

               elseif (to_lower(types(nval:nval)).eq.'d' .or. to_lower(types(nval:nval)).eq.'a') then

                  ! type 'a': angles in degrees given by '20.5d' are converted to radians

                  is_deg = .false.
                  if (to_lower(types(nval:nval)).eq.'a') then
                     ilen = len(trim(words(ii)))
                     ! check whether last character is 'd' -> angle in degrees
                     if (to_lower(words(ii)(ilen:ilen)).eq.'d') is_deg = .true.
                     ! strip last character 'd' or 'r'
                     if (to_lower(words(ii)(ilen:ilen)).eq.'d' .or.                                     &
                         to_lower(words(ii)(ilen:ilen)).eq.'r') then
                        words(ii) = words(ii)(1:ilen-1)
                     endif
                  endif
   
                  ! obtain a real(kind=8) value from string words(ii)
                  ! Note: according to C. Page, "Prof.Prog's Guide", p.93, the format "(Fw.0)" will cope
                  !       with all common floating point formats, including integer forms.
   
                  idbl = idbl + 1
                  read(words(ii), fmt='(f40.0)', IOSTAT=iostat) dbles(idbl)
                  if (idebug.ge.9) then
                     write(bufout,'(a,i3,a,i2,a,e15.4,a,i7,a,i3,3a)') 'item', ii,                       &
                        ', type 2 -> dbles(', idbl, ') = ',dbles(idbl),' ; iostat = ', iostat,          &
                        ' ; words(',ii,') "' , trim(words(ii)), '"'
                     call write_log(1, bufout)
                  endif
   
                  if (is_deg) dbles(idbl) = dbles(idbl) * pi/180d0
   
                  if (is_deg .and. idebug.ge.9) then
                     write(bufout,'(a,i2,a,e15.4)') '       in radians: dbles(', idbl, ') = ', dbles(idbl)
                     call write_log(1, bufout)
                  endif
                  if (iostat.ne.0) idbl = idbl - 1
   
               elseif (to_lower(types(nval:nval)).eq.'l') then
   
                  ! obtain a logical value from string words(ii)
   
                  iflg = iflg + 1
                  read(words(ii), fmt='(l)', IOSTAT=iostat) flags(iflg)
                  if (idebug.ge.9) then
                     write(bufout,'(2(a,i3),a,l5,a,i7,a,i3,3a)') 'item', ii,                            &
                        ', type 3 -> flags(', iflg, ') = ',flags(iflg),' - iostat = ', iostat,          &
                        ' - words(', ii,' ) = "', trim(words(ii)), '"'
                     call write_log(1, bufout)
                  endif
                  if (iostat.ne.0) iflg = iflg - 1
   
               elseif (to_lower(types(nval:nval)).eq.'s') then
   
                  ! obtain a string value from string words(ii)
   
                  istr = istr + 1
                  ilen = len(trim(words(ii)))
                  if (words(ii)(1:1).ne.squote .or. words(ii)(ilen:ilen).ne.squote) then
                     iostat = 1
                     strngs(istr) = ' '
                  else
                     iostat = 0
                     strngs(istr) = words(ii)(2:ilen-1)
                  endif
                  if (idebug.ge.9) then
                     write(bufout,'(2(a,i3),3a,i7,a,i3,3a)') 'item',ii,                                 &
                        ', type 4 -> strngs(', istr, ') = "', trim(strngs(istr)),'", iostat = ',        &
                        iostat, ' - words(',ii,' ) = "', trim(words(ii)), '"'
                     call write_log(3, bufout)
                  endif
                  if (iostat.ne.0) istr = istr - 1
   
               else

                  write(bufout,*) 'ERROR reading "', trim(descrip),'": ',                               &
                     'the string types must consist only of i, d, a, l and s'
                  call write_log(1, bufout)
                  if (lstop) call abort_run()
                  ierror = -6

               endif

               ! in case of an error in parsing: write error message

               if (iostat.ne.0 .and. nval.gt.nval_req) then

                  ! no success parsing optional value

                  nval = nval - 1

               elseif (iostat.ne.0) then

                  ! no success parsing mandatory value

                  write(bufout,'(2(a,i7),a,2(/,3a),i3,a)') ' ERROR processing input for case', ncase,   &
                     ', line', linenr,':',                                                              &
                     '  "',trim(inptxt),'"',                                                            &
                     ' expecting "',trim(descrip),'"; element',ii,' is wrong or missing'
                  call write_log(3, bufout)
                  if (lstop) call abort_run()
                  ierror = -7
               endif

            end do ! end do ii = 1, nitem (one line)

         endif ! reached end-of-file

      enddo ! end while (nval < nval_req): more lines

      ! debug output

      if (idebug.ge.7) then
         write(bufout, '(a,a,/, a,a,/, a,i6,/, a,i6)') 'descrip: ', descrip, 'types:   ', types,       &
            'linenr:  ', linenr,  'iostat:  ', iostat
         call write_log(4, bufout)
         write(bufout,*) iint,' ints:    ',(ints(ii), ii=1,min(4,iint))
         call write_log(1, bufout)
         bufout(2) = ' '
         write(bufout,*) idbl,' dbles:   ',(dbles(ii), ii=1,min(4,idbl))
         ii = max(1, min(4,idbl+1)/2 )
         call write_log(ii, bufout)
         write(bufout,*) iflg,' flags:   ',(flags(ii), ii=1,min(6,iflg))
         call write_log(1, bufout)
         write(bufout,*) istr,' strings: '
         call write_log(1, bufout)
         do ii = 1, istr
            call write_log(trim(strngs(ii)))
         enddo
      endif

      ! return current end-of-file status

      ieof_arg = ieof

   end subroutine readline

!------------------------------------------------------------------------------------------------------------

   subroutine read1darr(unitnm, ncase, linenr, descrip, types, nval, ints, dbles, idebug, lstop, ierror)
!--purpose: obtain nval ints (types = 'i') or double values (types = 'd' or 'a') from input file "unitnm",
!           or write an error message (descrip, linenr) in case of errors in the user input.
      implicit none
!--subroutine arguments:
      integer,      intent(in)    :: unitnm, ncase, nval, idebug
      integer,      intent(inout) :: linenr
      integer,      intent(out)   :: ierror
      logical,      intent(in)    :: lstop
      integer                     :: ints(nval)   ! output when types = 'i'
      real(kind=8)                :: dbles(nval)  ! output when types = 'd' or 'a'
      character(len=*), intent(in)    :: descrip, types
!--local variables:
      integer,      parameter :: mxitem = 50         ! max #values per input-line
      real(kind=8), parameter :: pi     = 4d0*atan(1d0)
      integer            :: nblank, ncmtln, iostat, ieof, ilen, ii, ival, nitem
      logical            :: use_ints, use_angles, is_deg
      character(len= 1 ) :: ext
      character(len=MAX_CHAR_INP) :: inptxt, words(mxitem)

      if (idebug.ge.5) then
         write(bufout,'(/,3a,i5,3a)') '--- read1darr: ', trim(descrip),',', nval,' values of type "', &
                 types, '"'
         call write_log(2, bufout)
      endif

      if (to_lower(types(1:1)).eq.'i') then
         use_ints   = .true.
         use_angles = .false.
      elseif (to_lower(types(1:1)).eq.'d') then
         use_ints   = .false.
         use_angles = .false.
      elseif (to_lower(types(1:1)).eq.'a') then
         use_ints   = .false.
         use_angles = .true.
      else
         write(bufout,*) 'ERROR reading "', trim(descrip),'": the string types may only contain i, d or a'
         call write_log(1, bufout)
         ierror = -6
         if (lstop) stop
      endif

      ieof   = -1  ! end-of-file is considered an error
      ival   =  0
      ierror =  0

      ! while (ival < nval) do

      do while (ival.lt.nval .and. ierror.ge.0)

         ! Increment line number and read input, thereby skipping comments and empty lines

         call getnonemptyline(unitnm, ncase, descrip, '%', idebug, ieof, lstop, linenr, nblank,         &
                ncmtln, inptxt)

         ! split the line into separate words

         call splitlinetowords(inptxt, '%', mxitem, words, nitem)

         ! process the strings, read values to ints or dbles array

         do ii = 1, nitem
            if (ival.lt.nval) then

               ival = ival + 1

               if (use_ints) then

                  ! obtain an integer value from string words(ii)

                  read(words(ii), fmt='(i)', IOSTAT=iostat) ints(ival)
                  if (idebug.ge.8) then
                     write(bufout,'(a,i5,a,i7,a,i7,a,i3,3a)') ' ints(', ival, ') = ',ints(ival),        &
                           ' ; iostat = ', iostat, ' ; words(',ii,') "' , trim(words(ii)), '"'
                     call write_log(1, bufout)
                  endif

               else

                  ! type 'a': angles in degrees given by '20.5d' are converted to radians

                  is_deg = .false.
                  if (use_angles) then
                     ilen = len(trim(words(ii)))
                     ! check whether last character is 'd' -> angle in degrees
                     ext = to_lower(words(ii)(ilen:ilen))
                     if (ext.eq.'d') is_deg = .true.
                     ! strip last character 'd' or 'r'
                     if (ext.eq.'d' .or. ext.eq.'r') words(ii) = words(ii)(1:ilen-1)
                  endif

                  ! obtain a real(kind=8) value from string words(ii)
                  ! Note: according to C. Page, "Prof.Prog's Guide", p.93, the format "(Fw.0)" will cope
                  !       with all common floating point formats, including integer forms.

                  read(words(ii), fmt='(f40.0)', IOSTAT=iostat) dbles(ival)
                  if (idebug.ge.8) then
                     write(bufout,'(a,i5,a,e15.4,a,i7,a,i3,3a)') ' dbles(', ival, ') = ',dbles(ival),   &
                           ' ; iostat = ', iostat, ' ; words(',ii,') "' , trim(words(ii)), '"'
                     call write_log(1, bufout)
                  endif

                  if (is_deg) dbles(ival) = dbles(ival) * pi/180d0
                  if (is_deg .and. idebug.ge.9) then
                     write(bufout,'(a,i5,a,e15.4)') '        in radians: dbles(',ival,') = ',dbles(ival)
                     call write_log(1, bufout)
                  endif

               endif ! use_ints

               ! in case of an error in parsing: write error message

               if (iostat.ne.0) then
                  write(bufout,'(a,i7,a,i7,a,/,3a,/,3a,i7,a)') ' ERROR processing input for case',      &
                     ncase,', line', linenr,':',                                                        &
                     '  "',trim(inptxt),'"',                                                            &
                     ' expecting "',trim(descrip),'"; element',ival,' is wrong or missing'
                  call write_log(3, bufout)
                  ierror = -2
                  if (lstop) stop
               endif

            endif ! ival < nval
         enddo ! ii=1,nitem
      enddo ! end while (ival < nval)

      ! write debug output when requested

      if (idebug.ge.7) then
         write(bufout, '(a,a,/, a,i6,/, a,i6,/, a,i6)') 'descrip: ',descrip, 'linenr:  ',linenr,       &
            'iostat:  ',iostat, 'nval:  ',nval
         call write_log(4, bufout)

         if (use_ints) then
            write(bufout, '(a,8i10)') 'ints: ', ints(1:min(8,nval))
            call write_log(1, bufout)
         else
            bufout(2) = ' '
            write(bufout, *) 'dbles: ', dbles(1:min(4,nval))
            ii = max(1, min(4,nval+1)/2 )
            call write_log(ii, bufout)
         endif
      endif

      ! stop in case of error

      if (iostat.ne.0) then
         write(bufout,'(a,a,a,i7,a)') 'ERROR reading ', descrip, ': an element of line ', linenr,      &
            ' of the input file is wrong or missing'
         call write_log(1, bufout)
         ierror = -3
         if (lstop) stop
      endif

   end subroutine read1darr

!------------------------------------------------------------------------------------------------------------

      subroutine readkeywordvalues(inptxt, commnt, mxitem, idebug, has_key, keywrd, has_str, valstr,    &
                        nitem, values)
!--purpose: Split a line of input-file: 1) keyword='string', keyword=real1, or keyword=<empty>, 
!                                       2) plain keyword w.o. '=',
!                                       3) numeric values real1 [real2...]
!           Detect keywords and values, ignoring comma's and comments but honouring single quotes
      implicit none
!--subroutine arguments:
      character(len=*), intent(in)  :: inptxt              ! a complete line of input, not a comment line
      character(len=*), intent(in)  :: commnt              ! comment character(s)
      integer,       intent(in)  :: mxitem              ! length of output-array
      integer,       intent(in)  :: idebug              ! level of debug output
      logical,       intent(out) :: has_key             ! T/F there's a keyword present
      character(len=*), intent(out) :: keywrd              ! keyword (if found)
      logical,       intent(out) :: has_str             ! T/F there's a value string present
      character(len=*), intent(out) :: valstr              ! value (if not numeric)
      integer,       intent(out) :: nitem               ! number of numeric values
      real(kind=8),  intent(out) :: values(1:mxitem)    ! numeric output values
!--local variables:
      integer, parameter :: itabvl = ichar('\t')
      character(len=1)   :: squote = '''', dquote = '"'
      integer            :: ncmt_char, ic, ipos, ios, ix, k, nval
      logical            :: in_str, in_cmt, in_val, has_equals
      real(kind=8)       :: tmpval(mxitem)
      character(len=len(inptxt)) :: mytext

      if (idebug.ge.8) call write_log('   ...starting keyVal for str=' // trim(inptxt) // '.')

      ncmt_char = len(commnt)

      ! Locate comments, everything after first non-quoted comment-sign

      ipos = 0
      in_str = .false.
      in_cmt = .false.
      do while (ipos.lt.len(inptxt) .and. .not.in_cmt)
         ipos = ipos + 1
         if (inptxt(ipos:ipos).eq.squote .or. inptxt(ipos:ipos).eq.dquote) in_str = .not.in_str
         do ic = 1, ncmt_char
            if (.not.in_str .and. inptxt(ipos:ipos).eq.commnt(ic:ic)) in_cmt = .true.
         enddo
      enddo

      ! Copy non-comment part of line

      if (ipos.le.0) then
         mytext = inptxt
      elseif (ipos.eq.1) then
         mytext = ' '
      else
         mytext = inptxt(1:ipos-1)
      endif
      if (idebug.ge.7) call write_log('   ...removed comment, mytext=' // trim(mytext) // '.')

      ! Locate key/value split, split at first non-quoted equals-sign

      ipos = 0
      in_str = .false.
      in_val = .false.
      do while (ipos.lt.len(mytext) .and. .not.in_val)
         ipos = ipos + 1
         if (mytext(ipos:ipos).eq.squote .or. inptxt(ipos:ipos).eq.dquote) in_str = .not.in_str
         if (.not.in_str .and. mytext(ipos:ipos).eq.'=') in_val = .true.
      enddo
      has_equals = in_val
      if (has_equals) then
         if (idebug.ge.7) then
            write(bufout,'(a,i6)') '   ...found equals at ipos=',ipos
            call write_log(1, bufout)
         endif
      else
         if (idebug.ge.7) call write_log('   ...no equals found in string mytext')
      endif

      ! Three cases: 1) keyword=value, keyword=<empty>, 2) plain keyword w.o. '=', 3) only numeric values, 

      if (has_equals) then

         ! 1) keyword=value, keyword=<empty>

         ! copy keyword without leading blanks

         has_key = .true.
         if (ipos.eq.1) then
            keywrd = ' '
         else
            keywrd = adjustl(mytext(1:ipos-1))
         endif

         ! get value-string

         if (ipos.lt.len(mytext)) then
            has_str = .true.
            valstr = adjustl(mytext(ipos+1:))
            if (len(trim(valstr)).le.0) has_str = .false.
         else
            has_str = .false.
            valstr = ' '
         endif
         if (idebug.ge.7) then
            write(bufout,'(5a)') '   ...obtained key/value, str=',trim(keywrd),', ',trim(valstr),'.'
            call write_log(1, bufout)
         endif

         ! check for '/' in value-string - used in Miniprof dates

         ix = index(valstr, '/')
         if (ix.gt.0) then

            nval = 0

         else

            ! count number of values in value-string

            nval = 0
            ios = 0
            do while (nval.lt.mxitem .and. ios.eq.0)
               ! attempt to read nval+1 values
               read(valstr,*,iostat=ios,err=899) (tmpval(k), k=1,nval+1)
               ! if ok, increment nval, else jump out of loop (err=)
               if (ios.eq.0) nval = nval + 1
            enddo
899         continue
         endif
         if (idebug.ge.7) then
            write(bufout,'(3a,i4,a)') '      valstr="',trim(valstr),'" has',nval,' values'
            call write_log(1, bufout)
         endif

         ! in case of numeric values, copy to values array

         if (nval.gt.0) then
            read(valstr,*) (values(k), k=1,nval)
            nitem = nval
         else
            nitem = 0
         endif

         if (idebug.ge.6 .and. nitem.gt.0) then
            write(bufout,'(3a,20(10g12.4,:,/,20x))') '   ...key "',trim(keywrd),'": value=',            &
                                        (values(k),k=1,nitem)
            call write_log(1, bufout)
         elseif (idebug.ge.6) then
            write(bufout,'(5a)') '   ...key "',trim(keywrd),'": value="',trim(valstr),'"'
            call write_log(1, bufout)
         endif

      else ! (.not.has_equals)

         ! 2), 3): No '='-sign, can be plain keyword or list of plain values

         has_str = .false.
         valstr  = ' '

         ! Count the number of numeric values in mytext

         nval = 0
         ios = 0
         do while (nval.lt.mxitem .and. ios.eq.0)
            ! attempt to read nval+1 values
            read(mytext,*,iostat=ios,err=999) (tmpval(k), k=1,nval+1)
            ! if ok, increment nval, else jump out of loop (err=)
            if (ios.eq.0) nval = nval + 1
         enddo
999      continue
         if (idebug.ge.7) then
            write(bufout,'(3a,i4,a)') '   ...mytext="',trim(mytext),'" has',nval,' values'
            call write_log(1, bufout)
         endif

         ! in case of numeric values, copy to values array and reset keywrd string
         ! else no values and set keywrd string

         if (nval.gt.0) then
            read(mytext,*) (values(k), k=1,nval)
            nitem = nval
            has_key = .false.
            keywrd = ' '
            if (idebug.ge.6) then
               write(bufout,'(a,20(10g12.4,:,/,10x))') '   ...values=',(values(k),k=1,nitem)
               call write_log(1, bufout)
            endif
         else
            keywrd  = adjustl(mytext)
            has_key = .true.
            nitem   = 0
            if (idebug.ge.6) call write_log('   ...key "' // trim(keywrd) // '"')
         endif

      endif ! has_equals

      end subroutine readkeywordvalues

!------------------------------------------------------------------------------------------------------------

      subroutine readplainvalues(inptxt, commnt, nexpct, idebug, nitem, values)
!--purpose: Read nexpct numeric values from inptxt, with trailing comments indicated e.g. by commnt='!%' 
      implicit none
!--subroutine arguments:
      character(len=*), intent(in)  :: inptxt              ! a complete line of input, not a comment line
      character(len=*), intent(in)  :: commnt              ! comment character(s)
      integer,          intent(in)  :: nexpct              ! #expected values, length of output-array
      integer,          intent(in)  :: idebug              ! level of debug output
      integer,          intent(out) :: nitem               ! number of numeric values
      real(kind=8),     intent(out) :: values(1:nexpct)    ! numeric output values
!--local variables:
      integer            :: ncmt_char, ic, ipos, ios, iend, k, nval
      real(kind=8)       :: tmpval(nexpct)
      character(len=len(inptxt)) :: mytext

      ! Try to read the expected number of numeric values from inptxt

      read(inptxt,*,iostat=ios) (values(k), k=1,nexpct)

      ! If successful, return

      if (ios.eq.0) then

         nitem = nexpct

      ! Else start extended processing

      else

         ! Locate comments, everything after first comment-sign (disregarding quotes in inptxt)

         ncmt_char = len(commnt)
         iend = len(inptxt)

         do ic = 1, ncmt_char
            ipos = index(inptxt(1:iend), commnt(ic:ic))
            if (ipos.gt.0) iend = ipos-1
         enddo

         ! Copy non-comment part of line

         if (iend.le.0) then
            mytext = inptxt
         elseif (iend.eq.1) then
            mytext = ' '
         else
            mytext = inptxt(1:iend)
         endif
         if (idebug.ge.7) call write_log('   ...removed comment, mytext=' // trim(mytext) // '.')

         ! Determine number of values available

         nval = 0
         ios = 0
         do while (nval.lt.nexpct .and. ios.eq.0)
            ! attempt to read nval+1 values
            read(mytext,*,iostat=ios,err=999) (tmpval(k), k=1,nval+1)
            ! if ok, increment nval, else jump out of loop (err=)
            if (ios.eq.0) nval = nval + 1
         enddo
 999     continue

         if (idebug.ge.7) then
            write(bufout,'(3a,i4,a)') '   ...mytext="',trim(mytext),'" has',nval,' values'
            call write_log(1, bufout)
         endif

         ! Store available values

         nitem = nval
         read(mytext,*) (values(k), k=1,nval)

         if (idebug.ge.6) then
            write(bufout,'(a,20(10g12.4,:,/,10x))') '   ...values=',(values(k),k=1,nitem)
            call write_log(1, bufout)
         endif

      endif ! ios=0, shortcut processing

      end subroutine readplainvalues

!------------------------------------------------------------------------------------------------------------

   subroutine getnonemptyline(unitnm, ncase, descrip, commnt, idebug, ieof, lstop, linenr, nblank,      &
                ncmtln, line)
!--purpose: Read lines from unit unitnm until a non-empty line is found.
!           return the whole line, including trailing comments
      implicit none
!--subroutine arguments:
      integer,       intent(in)    :: unitnm, ncase, idebug
      logical,       intent(in)    :: lstop
      character(len=*), intent(in)    :: descrip
      character(len=*), intent(in)    :: commnt
      integer,       intent(inout) :: ieof, linenr
                                        ! on input,  ieof<0  means end-of-file is considered an error,
                                        !            ieof>=0 means end-of-file is ok
                                        ! on output, ieof=1 means that end-of-file is encountered
      integer,       intent(out)   :: nblank, ncmtln
                                        ! count number of blank lines & comment lines skipped over
      character(len=*), intent(out)   :: line
!--local variables:
      integer iostat, ipos, ncmt_char, ic, ix
      logical lcomment

      ncmt_char = len(commnt)
      ! write(bufout,*) trim(descrip),': ncmt_char=',ncmt_char,', commnt=', (commnt(ic:ic),',',         &
      !         ic=1, ncmt_char)
      ! call write_log(1, bufout)

      ! count number of blank lines and comment lines skipped (Miniprof: detect empty line)

      nblank = 0
      ncmtln = 0

      ! increment line number and read input;
      ! check whether line is empty if comments are removed, until a non-empty line is found.

      lcomment = .true.
      do while (lcomment .and. ieof.le.0)

         ! read line from file, increment line-number

         read(unitnm,'(a)',iostat=iostat) line
         linenr = linenr + 1

         ! check for I/O errors
         ! Note: end-of-file is considered an error when ieof<0 is specified

         if (iostat.lt.0 .and. ieof.lt.0) then

            ! Error: reached end of file

            write(bufout,'(a,i7,a,i7,a,/,3a)') ' ERROR processing input for case',ncase,               &
               ', line',linenr,':','       reached end of file while expecting "',trim(descrip),'"'
            call write_log(2, bufout)
            ieof = 1
            if (lstop) stop

         elseif (iostat.lt.0) then

            ! Return: reached end of file
            ! Note: end-of-file is considered ok when ieof=0 is specified

            ieof = 1

         elseif (iostat.gt.0) then

            ! Error in reading line from file

            write(bufout,'(a,i7,a,i7,a,/,a,i6,/,3a)') ' ERROR processing input for case',ncase,      &
               ', line',linenr,':','      low-level file-reading-error, IOSTAT=',iostat,                  &
               ' while expecting      "',trim(descrip),'"'
            call write_log(3, bufout)
            if (lstop) stop
         endif

         if (ieof.gt.0) then

            ! reached end of file

            line = ' '

         else

            if (idebug.ge.4) then
               write(bufout,'(i5,2a)') linenr, 'g: ', trim(line)
               call write_log(1, bufout)
            endif

            ! Locate comments, everything after first 'commnt' sign
            ! index=0 if comment-character not found

            ipos = 999
            do ic = 1, ncmt_char
               ix = index(line, commnt(ic:ic))
               if (ix.gt.0) ipos = min(ipos, ix)
            enddo
            if (ipos.eq.999) ipos = 0

            ! Determine whether the line is non-empty if comments are ignored

            if (ipos.le.0) then         ! may be empty or non-empty with no comments
               lcomment = (line.eq.' ')
               if (lcomment) nblank = nblank + 1
            elseif (ipos.eq.1) then     ! entirely comment line
               lcomment = .true.
               ncmtln = ncmtln + 1
            else                        ! may be empty or non-empty + comments
               lcomment = (line(1:ipos-1).eq.' ')
               if (lcomment) ncmtln = ncmtln + 1
            endif

            ! write(bufout,*) 'ipos=',ipos,', lcomment=',lcomment
            ! call write_log(1,bufout)
         endif

      enddo

   end subroutine getnonemptyline

!------------------------------------------------------------------------------------------------------------

   subroutine splitlinetowords(inptxt, commnt, mxitem, words, nitem)
!--purpose: Split a line of input at spaces into separate words, thereby ignoring comma's and comments
!           but honouring single quotes, return in array of character strings.
      implicit none
!--subroutine arguments:
      character(len=*), intent(in)  :: inptxt
      character(len=*), intent(in)  :: commnt
      integer,       intent(in)  :: mxitem
      character(len=*), intent(out) :: words(1:mxitem)
      integer,       intent(out) :: nitem
!--local variables:
      integer, parameter :: itabvl = ichar('\t')
      character(len=1)   :: squote = ''''
      integer            :: ncmt_char, ipos, ii, ic
      logical            :: in_str, in_cmt
      character(len=len(inptxt)) :: mytext

      ncmt_char = len(commnt)

      ! Locate comments, everything after first non-quoted comment-sign

      ipos = 0
      in_str = .false.
      in_cmt = .false.
      do while (ipos.lt.len(inptxt) .and. .not.in_cmt)
         ipos = ipos + 1
         if (inptxt(ipos:ipos).eq.squote) in_str = .not.in_str
         do ic = 1, ncmt_char
            if (.not.in_str .and. inptxt(ipos:ipos).eq.commnt(ic:ic)) in_cmt = .true.
         enddo
      enddo

      ! Copy non-comment part of line

      if (ipos.le.0) then
         mytext = inptxt
      elseif (ipos.eq.1) then
         mytext = ' '
      else
         mytext = inptxt(1:ipos-1)
      endif

      ! Replace commas in input line by spaces

      ii = index(mytext,',')
      do while (ii.gt.0)
         mytext(ii:ii) = ' '
         ii = index(mytext,',')
      enddo

      ! initialise the output string array

      do ii = 1, mxitem
         words(ii) = ' '
      enddo

      ! split the input at spaces and '='-signs and put into string array

      ipos = 0
      nitem = 0
      in_str = .false.

      ! TODO: more than mxitem items in a line, "ii-1" for ii=1.

      ! loop over all characters
      !    if (in string parameter or not a space/equals character)
      !       set the correct position in the "words" array
      !       copy the character to the correct position in the "words" array
      !    endif

      do ii = 1, len(mytext)

         ! if (in string parameter or not a space/equals character)

         if (in_str .or. (ichar(mytext(ii:ii)).ne.itabvl .and.                                         &
                          mytext(ii:ii).ne.' ' .and. mytext(ii:ii).ne.'=')) then

            ! very first character:

            if (ii.eq.1) then

               ! start new element of array, nitem:=1, ipos=0

               nitem = nitem + 1
               ipos  = 0
               in_str = (mytext(ii:ii).eq.squote)

            ! previous character was a space/equals and does not belong to a string:

            elseif (.not.in_str .and. (ichar(mytext(ii-1:ii-1)).eq.itabvl .or.                         &
                                       mytext(ii-1:ii-1).eq.' ' .or. mytext(ii-1:ii-1).eq.'=')) then

               ! start new element of array, nitem+=1, ipos=0

               nitem = nitem + 1
               ipos  = 0
               in_str = (mytext(ii:ii).eq.squote)

            ! terminating quote of a string is found

            elseif (in_str .and. mytext(ii:ii).eq.squote) then

               in_str = .false.

            endif
            ipos = ipos + 1
            words(nitem)(ipos:ipos) = mytext(ii:ii)
         endif
      end do

   end subroutine splitlinetowords

!------------------------------------------------------------------------------------------------------------

   logical function check_irng (descrp, ival, imin, imax)
      implicit none
      character(len=*) descrp
      integer ival, imin, imax
      logical zerror

      zerror = .false.
      if (ival.lt.imin .or. ival.gt.imax) then
         zerror = .true.
         write(lout, 1000) descrp, ival, imin, imax
         write(   *, 1000) descrp, ival, imin, imax
 1000    format (' Input: ERROR. ',a,' (=',i5,') is out of range: ', i4,' <= value <=',i6,'.')
      endif
      check_irng = .not.zerror
   end function check_irng

!------------------------------------------------------------------------------------------------------------

   logical function check_2rng (descrp, ival, imin1, imax1, imin2, imax2)
      implicit none
      character(len=*) descrp
      integer ival, imin1, imax1, imin2, imax2
      logical zerror

      zerror = .false.
      if (.not.( (ival.ge.imin1 .and. ival.le.imax1) .or. (ival.ge.imin2 .and. ival.le.imax2) ) ) then
         zerror = .true.
         write(lout, 1000) descrp, ival, imin1, imax1, imin2, imax2
         write(   *, 1000) descrp, ival, imin1, imax1, imin2, imax2
 1000    format (' Input: ERROR. ',a,' (=',i6,') is out of range: ',/,                                 &
                 7x,i6,' <= value <=',i6,' or ',i6,' <= value <=',i6,'.')
      endif
      check_2rng = .not.zerror
   end function check_2rng

!------------------------------------------------------------------------------------------------------------

   logical function check_2int (descrp, ival, ival1, ival2)
      implicit none
      character(len=*) descrp
      integer ival, ival1, ival2
      logical zerror

      zerror = .false.
      if (.not.( ival.eq.ival1 .or. ival.ge.ival2 ) ) then
         zerror = .true.
         write(lout, 1000) descrp, ival, ival1, ival2
         write(   *, 1000) descrp, ival, ival1, ival2
 1000    format (' Input: ERROR. ',a,' (=',i6,') is invalid, should be',i6,' or ',i6,'.')
      endif
      check_2int = .not.zerror
   end function check_2int

!------------------------------------------------------------------------------------------------------------

   logical function check_3rng (descrp, ival, imin1, imax1, imin2, imax2, imin3, imax3)
      implicit none
      character(len=*) descrp
      integer ival, imin1, imax1, imin2, imax2, imin3, imax3
      logical zerror

      zerror = .false.
      if (.not.( (ival.ge.imin1 .and. ival.le.imax1) .or. (ival.ge.imin2 .and. ival.le.imax2) .or.     &
                 (ival.ge.imin3 .and. ival.le.imax3) ) ) then
         zerror = .true.
         write(lout, 1000) descrp, ival, imin1, imax1, imin2, imax2, imin3, imax3
         write(   *, 1000) descrp, ival, imin1, imax1, imin2, imax2, imin3, imax3
 1000    format (' Input: ERROR. ',a,' (=',i6,') is out of range: ',/,                                 &
                 7x,i6,'--',i6,' or ',i6,'--',i6,' or ',i6,'--',i6,'.')
      endif
      check_3rng = .not.zerror
      end function check_3rng

!------------------------------------------------------------------------------------------------------------

   logical function check_range (descrp, dval, dmin, dmax)
!--purpose: perform check on input-value dval: dval \in [dmin, dmax]
      implicit none
      character(len=*) :: descrp
      real(kind=8)     :: dval, dmin, dmax
      logical          :: zerror

      zerror = .false.
      if (dval.lt.dmin) then
         zerror = .true.
         write(lout, 1000) descrp, dval, dmin
         write(   *, 1000) descrp, dval, dmin
 1000    format (' Input: ERROR. ',a,'=',g12.4,' may not be less than', g12.4,'.')
      elseif (dval.gt.dmax) then
         zerror = .true.
         write(lout, 2000) descrp, dval, dmax
         write(   *, 2000) descrp, dval, dmax
 2000    format (' Input: ERROR. ',a,'=',g12.4,' may not be larger than',g12.4,'.')
      endif
      check_range = .not.zerror
   end function check_range

!------------------------------------------------------------------------------------------------------------

   logical function check_equal (descrp1, dval1, descrp2, dval2, reltol, lprint)
!--purpose: perform check on input-values dval1 == dval2
      implicit none
      character(len=*)  :: descrp1, descrp2
      real(kind=8)      :: dval1, dval2, reltol
      logical, optional :: lprint
      logical           :: zerror, my_lprint

      my_lprint = .false.
      if (present(lprint)) my_lprint = lprint

      zerror = .false.
      if (abs(dval1-dval2).gt.reltol*min(abs(dval1),abs(dval2))) then
         zerror = .true.
         if (my_lprint) then
            write(lout, 1000) descrp1, dval1, descrp2, dval2
            write(   *, 1000) descrp1, dval1, descrp2, dval2
 1000       format (' Input: ERROR. ',a,'=',g12.4,' must be equal to ',a,'=',g12.4,'.')
         endif
      endif
      check_equal = .not.zerror
   end function check_equal

!------------------------------------------------------------------------------------------------------------

   logical function check_smaller (descrp1, dval1, descrp2, dval2)
!--purpose: perform check on input-values dval1 <= dval2
      implicit none
      character(len=*) :: descrp1, descrp2
      real(kind=8)     :: dval1, dval2
      logical          :: zerror

      zerror = .false.
      if (dval1.gt.dval2) then
         zerror = .true.
         write(lout, 1000) descrp1, dval1, descrp2, dval2
         write(   *, 1000) descrp1, dval1, descrp2, dval2
 1000    format (' Input: ERROR. ',a,'=',g12.4,' must be less than ',a,'=',g12.4,'.')
      endif
      check_smaller = .not.zerror
   end function check_smaller

!------------------------------------------------------------------------------------------------------------

   logical function check_sorted( namval, nval, val, ascend )
!--purpose: check whether array of values is strictly increasing (decreasing)
      implicit none
!--subroutine arguments:
      character(len=*), intent(in)  :: namval
      integer,          intent(in)  :: nval
      real(kind=8),     intent(in)  :: val(nval)
      logical,          intent(in)  :: ascend
!--local variables:
      real(kind=8), parameter :: thresh = 1d-9
      integer                 :: iv, ierror
 
      ! check that values are given in strictly increasing (decreasing) order

      ierror = 0
      iv     = 1
      do while(ierror.le.0 .and. iv.lt.nval)
         iv = iv + 1
         if (     ascend .and. val(iv)-val(iv-1).lt.thresh) ierror = iv
         if (.not.ascend .and. val(iv-1)-val(iv).lt.thresh) ierror = iv
      enddo

      if (ierror.ne.0) then
         if (ascend) then
            write(errmsg, 1000) namval, 'in', namval, ierror-1, val(ierror-1), namval, ierror, val(ierror)
         else                                                                                
            write(errmsg, 1000) namval, 'de', namval, ierror-1, val(ierror-1), namval, ierror, val(ierror)
         endif
         nline_errmsg = 2
 1000    format (' Input: ERROR. ',a,' must be strictly ',a,'creasing.', /,                             &
                 '               ',a,'(',i2,')=',f10.6,', ',a,'(',i2,')=',f10.6,'.')
      endif

      check_sorted = (ierror.eq.0)
   end function check_sorted

end module m_readline
