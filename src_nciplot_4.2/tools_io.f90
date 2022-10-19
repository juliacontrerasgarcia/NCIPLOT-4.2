! Copyright (c) 2013 Alberto Otero de la Roza <aoterodelaroza@ucmerced.edu>,
! Julia Conteras-Garcia <julia.contreras.garcia@gmail.com>,
! Erin R. Johnson <ejohnson29@ucmerced.edu>, and Weitao Yang
! <weitao.yang@duke.edu>
!
! nciplot is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! tools_io: Tools for managing input/output
module tools_io
   implicit none

   private
   public :: tictac
   public :: getargs
   public :: error
   public :: header
   public :: zatguess, nameguess
   public :: upper
   public :: realloc
   public :: getword, isinteger, isreal
   public :: DoRangeWarning, DoRangeWarning2
contains

   !> Outputs a tictac: date and time.
   function ffdate()
      use param
      character*(mline) :: ffdate

      integer :: values(8)

      call date_and_time(values=values)
      write (ffdate, '(I4.4,".",I2.2,".",I2.2,", ",&
         &I2.2,":",I2.2,":",I2.2,".",I3.3)') values(1:3), values(5:8)

   end function ffdate

   !> Interface to the timer routines of different computers. Returns
   !> machine seconds at the calling time.
   subroutine tictac(mesg)
      use param
      character*(*), intent(in) :: mesg !< Prefix message

      character*(mline) :: sdate

      sdate = ffdate()
      write (uout, '(A," -- ",A/)') trim(mesg), trim(adjustl(sdate))

   end subroutine tictac

   !> Get command line arguments and argument count
   !> on return, argc will contain the number of arguments and
   !> argv will be a character array containing the arguments.
   subroutine getargs(argc, argv)
      use param

      integer, intent(out) :: argc
      character*(*), intent(out) :: argv(2)

      character*(mline) :: line  !local string for parsing
      integer :: length         !total length of command line

      argc = 0
      length = 1
      do while (length > 0)
         call getarg(argc + 1, line)
         length = len(trim(adjustl(line)))
         if (length .gt. 0) then
            argc = argc + 1
            if (argc > 2) call error('getargs', 'Usage: nciplot [file.in [file.out]]', faterr)
            argv(argc) = line(1:length)
         endif
      end do

   end subroutine getargs

   !> Send an error message 'message' to stdout, coming from routine
   !> 'routine'. errortype is the error code (see mod_param.f90).
   subroutine error(routine, message, errortype)
      use param

      character*(*), intent(in) :: routine !< routine calling the error
      character*(*), intent(in) :: message !< the message
      integer, intent(in) :: errortype !< fatal, warning or info

      character*(20) chtype

      ! message styles shamelessly copied from abinit.
      if (errortype .eq. faterr) then
         chtype = 'ERROR'
      else if (errortype .eq. warning) then
         chtype = 'WARNING'
      else if (errortype .eq. noerr) then
         chtype = 'COMMENT'
      else
         chtype = 'UNKNOWN'
      endif
      write (uout, 100) trim(adjustl(chtype)), trim(adjustl(routine)), &
         trim(adjustl(message))
      if (errortype .eq. faterr) then
         write (stderr, 100) trim(adjustl(chtype)), trim(adjustl(routine)), &
            trim(adjustl(message))
         if (uin /= stdin) close (uin)
         if (uout /= stdout) close (uout)
         stop 1
      else if (errortype .eq. warning) then
         nwarns = nwarns + 1
      else if (errortype .eq. noerr) then
         ncomms = ncomms + 1
      endif

100   format(A, "(", A, "): ", A)

   end subroutine error

   subroutine header()
      use param

      write (uout, 1000)
1000  format(' # ----------------- NCIPLOT ------------------------', &
             /, ' # --- PLOTTING NON COVALENT INTERACTION REGIONS ----', /, &
             ' # ---             E.R. Johnson                  ----', /, &
             ' # ---          J. Contreras-Garcia              ----', /, &
             ' # ----------    Duke University         ------------', /, &
             ' #                                                   ', /, &
             ' # ---             A. de la Roza                  ---', /, &
             ' # --------- University of California Merced --------', /, &
             ' #                                                   ', /, &
             ' # ---               R. A. Boto                   ---', /, &
             ' # ---                 C. Quan                     --', /, &
             ' # --------  Université Pierre et Marie Curie -------', / &
             ' # --------------------------------------------------', /, &
             ' # ---              Please cite                  ----', /, &
             ' # --J. Am. Chem. Soc., 2010, 132 (18), pp 6498–6506-', /, &
             ' # --------------------------------------------------', / &
             ' # --------------------------------------------------', /, &
             ' # ---     Contributions for the wfn properties  ----', /, &
             ' # ---      from H. L. Schmider are acknowledged  ---', /, &
             ' # --------------------------------------------------', /, &
             ' # --------------------------------------------------', /, &
             ' # ---     Contributions for the wfx reader      ----', /, &
             ' # ---      from Dave Arias are acknowledged      ---', /, &
             ' # --------------------------------------------------', /, &
             ' # --------------------------------------------------', /, &
             ' # ---     Contributions for the integration --------', /, &
             ' # ---      algorithms from Erna Wieduwilt    --------',/, &
             ' # ---             are acknowledged           --------',/, &
             ' # ---------------------------------------------------',/, &   
             ' #')

   end subroutine header

   !> Guess the atomic number from the atomic name. A value
   !> of -1 is returned if the atom is not recognized.
   !> A valid atomic name is made of:
   !>   1) one or two capital or small case letters corresponding to
   !>      a legal atomic symbol
   !>   2) any other additional characters with the only excepcion that
   !>      the first character following (1) should not be a letter
   function zatguess(atname)

      integer :: zatguess !< Ouptut atomic number
      character*(*), intent(in) :: atname !< Input atomic symbol (case insensitive)

      character*(2) :: atsymbol
      integer       :: i, idx
      character*(2) :: an(1:103)

      data(an(i), i=1, 103) &
           & /'H ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O ', 'F ', 'NE',&
           &  'NA', 'MG', 'AL', 'SI', 'P ', 'S ', 'CL', 'AR', 'K ', 'CA',&
           &  'SC', 'TI', 'V ', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN',&
           &  'GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y ', 'ZR',&
           &  'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN',&
           &  'SB', 'TE', 'I ', 'XE', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND',&
           &  'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB',&
           &  'LU', 'HF', 'TA', 'W ', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG',&
           &  'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH',&
           &  'PA', 'U ', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM',&
           &  'MD', 'NO', 'LW'/

      idx = 1
      do while (atname(idx:idx) == ' ' .and. idx < len(atname) - 1)
         idx = idx + 1
      end do
      if (len(trim(adjustl(atname))) > 1) then
         atsymbol = atname(idx:idx + 1)
      else
         atsymbol = atname(idx:idx)
      end if
      call upper(atsymbol)
      if (atsymbol(2:2) .lt. 'A' .or. atsymbol(2:2) .gt. 'Z') atsymbol(2:2) = ' '
      do i = 1, 103
         if (atsymbol .eq. an(i)) then
            zatguess = i
            return
         endif
      enddo
      zatguess = -1
   end function zatguess

   !> The reverse of zatguess. Given the atomic number, return the
   !> chemical symbol of the atom, or XX if it is not recognized.
   function nameguess(zat)

      integer, intent(in) :: zat !< Input atomic number
      character*(2) :: nameguess !< Output atomic symbol

      character*(2) :: an(1:103)
      integer :: i

      data(an(i), i=1, 103) &
           & /'h_', 'he', 'li', 'be', 'b_', 'c_', 'n_', 'o_', 'f_', 'ne',&
           &  'na', 'mg', 'al', 'si', 'p_', 's_', 'cl', 'ar', 'k_', 'ca',&
           &  'sc', 'ti', 'v_', 'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn',&
           &  'ga', 'ge', 'as', 'se', 'br', 'kr', 'rb', 'sr', 'y_', 'zr',&
           &  'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in', 'sn',&
           &  'sb', 'te', 'i_', 'xe', 'cs', 'ba', 'la', 'ce', 'pr', 'nd',&
           &  'pm', 'sm', 'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm', 'yb',&
           &  'lu', 'hf', 'ta', 'w_', 're', 'os', 'ir', 'pt', 'au', 'hg',&
           &  'tl', 'pb', 'bi', 'po', 'at', 'rn', 'fr', 'ra', 'ac', 'th',&
           &  'pa', 'u_', 'np', 'pu', 'am', 'cm', 'bk', 'cf', 'es', 'fm',&
           &  'md', 'no', 'lw'/

      if (zat > 0 .and. zat < 104) then
         nameguess = an(zat)
      else
         nameguess = 'xx'
      end if

   end function nameguess

   !> In-place convert string to uppercase.
   subroutine upper(string)

      character*(*), intent(inout) :: string !< Input string, same as the function result on output

      integer :: i

      do i = 1, len(string)
         if (string(i:i) .ge. 'a' .and. string(i:i) .le. 'z') &
            string(i:i) = char(ichar(string(i:i)) + ichar('A') - ichar('a'))
      enddo

   end subroutine upper

   !> Adapt the size of an allocatable 1D real*8 array
   subroutine realloc(a, nnew)
      use param

      real*8, intent(inout), allocatable :: a(:) !< Input array, real*8, 1D
      integer, intent(in) :: nnew !< new dimension

      real*8, allocatable :: temp(:)
      integer :: nold

      if (.not. allocated(a)) &
         call error('realloc', 'array not allocated', faterr)
      nold = size(a)
      allocate (temp(nnew))

      temp(1:min(nnew, nold)) = a(1:min(nnew, nold))
      call move_alloc(temp, a)

   end subroutine realloc

   !> Get next word from line at lp and increase lp
   !>   - a word is defined as any sequence on nonblanks
   !>   - word and getword will return the same
   function getword(word, line, lp)

      character*(*), intent(out) :: word !< Word read, same as function result
      character*(*), intent(in) :: line !< Input line
      integer, intent(inout) :: lp !< Pointer to position on input line, updated after reading.
      character*(len(word)) :: getword

      integer wp, l, i

      l = len(word)
      do while (line(lp:lp) .eq. ' ' .or. line(lp:lp) .eq. char(9))   !skip blanks
         lp = lp + 1
      enddo
      if (line(lp:lp) .eq. char(10)) then
         !        word = char(10)//char(0)
         word = char(0)
      else
         wp = 1
         do while (wp .lt. l .and. line(lp:lp) .ne. ' ' .and. line(lp:lp) .ne. char(9) .and. &
                   line(lp:lp) .ne. char(0) .and. line(lp:lp) .ne. char(10))
            word(wp:wp) = line(lp:lp)
            wp = wp + 1
            lp = lp + 1
         enddo
         word(wp:wp) = char(0)

         !.clean rest of the word:
         !
         do i = wp + 1, l
            word(i:i) = ' '
         enddo
      endif
      getword = word

   end function getword

   !> Get a real number from line and sets rval to it.
   !> If a valid real number is not found, isreal returns .false.
   logical function isreal(rval, line, lp)
      use param

      character*(*), intent(in) :: line !< Input string
      integer, intent(inout) :: lp !< Pointer to current position on string
      real*8, intent(out) :: rval !< Real value read

      character*(mline) dumchar
      integer tp, i
      character*(1) ch
      logical matched, isdigit
      isdigit(ch) = ch .ge. '0' .and. ch .le. '9'

      do while (line(lp:lp) .eq. ' ')
         lp = lp + 1
      end do

      i = lp
      if (line(i:i) .eq. '+' .or. line(i:i) .eq. '-') i = i + 1
      if (isdigit(line(i:i))) then
         do while (isdigit(line(i:i)))
            i = i + 1
         enddo
         if (line(i:i) .eq. '.') then
            i = i + 1
            do while (isdigit(line(i:i)))
               i = i + 1
            enddo
         endif
         matched = .true.
      else if (line(i:i) .eq. '.') then
         i = i + 1
         if (isdigit(line(i:i))) then
            do while (isdigit(line(i:i)))
               i = i + 1
            enddo
            matched = .true.
         else
            matched = .false.
         endif
      else
         matched = .false.
      endif

      !.....get optional exponent
      tp = i - 1
      if (matched) then
         if (line(i:i) == 'e' .or. line(i:i) == 'E' .or. line(i:i) == 'd' .or. line(i:i) == 'D' .or. &
             line(i:i) == '-' .or. line(i:i) == '+') then
            i = i + 1
            if (line(i:i) .eq. '+' .or. line(i:i) .eq. '-') i = i + 1
            if (isdigit(line(i:i))) then
               do while (isdigit(line(i:i)))
                  i = i + 1
               enddo
               if (index(' '//','//char(0)//char(10), line(i:i)) .gt. 0) then
                  dumchar = line(lp:i - 1)//char(0)
                  rval = atof(dumchar)
                  lp = i
               else
                  matched = .false.
                  rval = 0d0
               endif
            else
               matched = .false.
            endif
         else
            if (index(' '//','//char(0)//char(10), line(i:i)) .gt. 0) then
               dumchar = line(lp:tp)//char(0)
               rval = atof(dumchar)
               lp = i
            else
               matched = .false.
               rval = 0d0
            endif
         endif
      else
         rval = 0d0
      endif
      !
      isreal = matched
      return
   end function isreal

   !> Get integer value from input text. If a valid integer is not
   !> found, then return .false.
   logical function isinteger(ival, line, lp)
      use param

      character*(mline), intent(in) :: line !< Input string
      integer, intent(inout) :: lp !< Pointer to current position on string
      integer, intent(out) :: ival !< Integer value read, same as function result

      integer :: i

      ival = 0
      isinteger = .false.
      do while (line(lp:lp) .eq. " ")
         lp = lp + 1
         if (lp > len(line)) return
      enddo
      i = lp
      if (line(i:i) .eq. '+' .or. line(i:i) .eq. '-') i = i + 1
      if (i > len(line)) return
      if (isdigit(line(i:i))) then
         do while (isdigit(line(i:i)))
            i = i + 1
            if (i >= len(line)) exit
         enddo
         if (line(i:i) .eq. " " .or. line(i:i) .eq. char(0) .or. line(i:i) .eq. char(10)) then
            ival = atoi(line(lp:i - 1))
            lp = i
            isinteger = .true.
         elseif (i == len(line)) then
            ival = atoi(line(lp:i))
            lp = i + 1
            isinteger = .true.
         else
            return
         endif
      else
         return
      endif

   end function isinteger

   !> Return true if c is a digit (private).
   logical function isdigit(c)

      character*(1), intent(in) :: c !< Is the character c a digit?

      isdigit = c .ge. '0' .and. c .le. '9'

   end function isdigit

   !> Convert ascii string to integer (private).
   integer function atoi(string)

      character*(*), intent(in) :: string !< Input string

      integer i, sign

      atoi = 0
      i = 1
      do while (string(i:i) .eq. " ")
         i = i + 1
      end do
      sign = 1
      if (string(i:i) .eq. '+' .or. string(i:i) .eq. '-') then
         if (string(i:i) .eq. '-') sign = -1
         i = i + 1
      endif
      do while (isdigit(string(i:i)))
         atoi = 10*atoi + ichar(string(i:i)) - ichar('0')
         i = i + 1
         if (i > len(string)) then
            exit
         end if
      end do
      atoi = atoi*sign
      return

   end function atoi

   !> Convert ascii string to real (private).
   real*8 function atof(str)

      character*(*), intent(in) :: str !< Input string

      real*8, parameter :: ten = 10d0

      real*8 val, power
      integer exponent, sign, esign, i

      sign = 1
      val = 0d0
      power = 1d0
      exponent = 0
      esign = 1
      i = 1
      do while (str(i:i) .eq. ' ')
         i = i + 1
      end do
      if (str(i:i) .eq. '+' .or. str(i:i) .eq. '-') then
         if (str(i:i) .eq. '-') sign = -1
         i = i + 1
      endif
      do while (isdigit(str(i:i)))
         val = ten*val + ichar(str(i:i)) - ichar('0')
         i = i + 1
      enddo
      if (str(i:i) .eq. '.') then
         i = i + 1
         do while (isdigit(str(i:i)))
            val = ten*val + ichar(str(i:i)) - ichar('0')
            i = i + 1
            power = power*ten
         enddo
      endif
      if (str(i:i) .eq. 'e' .or. str(i:i) .eq. 'E' .or. str(i:i) .eq. 'd' .or. str(i:i) .eq. 'D') then
         i = i + 1
         if (str(i:i) .eq. '+' .or. str(i:i) .eq. '-') then
            if (str(i:i) .eq. '-') esign = -1
            i = i + 1
         endif
         do while (isdigit(str(i:i)))
            exponent = int(ten)*exponent + ichar(str(i:i)) - ichar('0')
            i = i + 1
         end do
      elseif (str(i:i) .eq. '-' .or. str(i:i) .eq. '+') then
         esign = 1
         if (str(i:i) .eq. '-') esign = -1
         i = i + 1
         do while (isdigit(str(i:i)))
            exponent = int(ten)*exponent + ichar(str(i:i)) - ichar('0')
            i = i + 1
         end do
      endif

      atof = (sign*val/power)*ten**(esign*exponent)
      return
   end function atof

   subroutine DoRangeWarning()
      use param
      write (uout, 5000)

5000  format('#                     Warning!                           #', /, &
             '#       Criterion for active box has changed.            #', /, &
             '# A box is defined as active if ALL its vertices are     #', /, &
             '# within the density and reduced density gradient range. #', /, &
             '#                                                        # ', /, &
             '#                                                        # ', /, &
             '                                                         ', /, &
             '                                                         ', /, &
             '                                                         ')

   end subroutine DoRangeWarning

   subroutine DoRangeWarning2()
      use param
      write (uout, 5001)

5001  format('#                     Warning!                           #', /, &
             '# Some range limit exceeds rho cutoff.                   #', /, &
             '# Points with sign(lambda2)rho value above electron      #', /, &
             '# density cutoff and below minus electron density cutoff #', /, &
             '# are not considered.                                    #', /, &
             '#                                                        # ', /, &
             '                                                           ', /, &
             '                                                           ', /, &
             '                                                           ')

   end subroutine DoRangeWarning2

end module tools_io
