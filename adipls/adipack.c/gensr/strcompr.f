      character*(*) function strcompr(s)
c
c  removes blanks from string s and returns result
c  Original version 12/12/04
c
      character*(*) s
      character*80 s1
      l=len(s)
c
      i1=0
      do i=1,l
	if(s(i:i).ne.' ') then
	  i1=i1+1
	  s1(i1:i1) = s(i:i)
        end if
      end do
c
      if(i1.gt.0) then
        strcompr=s1(1:i1)
      else
	strcompr=''
      end if
      return
      end
