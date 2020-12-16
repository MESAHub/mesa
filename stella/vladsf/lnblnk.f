      integer function lnblnk(string)
      character*(*) string
      integer strlen, i, j, len, null
c      integer strlen, i, j
      character*1 nullchar
      equivalence (nullchar, null)
      data null/0/
c      nullchar=char(null) 
      length = len(string)

      do i = 1, length
         if (string(i:i) .eq. nullchar) go to 10
      end do

 10   length = i - 1
c
      do i = 1, length
         j = length - i + 1
         if (string(j:j) .ne. ' ') then
            lnblnk = j
            return
         end if
      end do
c
      lnblnk = 0
      return
      end
