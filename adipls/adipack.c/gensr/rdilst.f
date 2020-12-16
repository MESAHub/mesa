      subroutine rdilst(ids,ia,n,ierr)
c
c  reads integer array ia(i), i = 1,..., n from unit ids, 
c  using list-directed read 
c  does not need place-holding commas at the end of input string
c  ierr is returned as 0, 1 or 2 for normal read, end of file
c  or error
c
c  Note: use of scratch file is extremely clumsy. However
c  list-directed input is not allowed in internal read
c
      character sa*80, comma*20, acomma*200
      dimension ia(n)
      data comma/',,,,,,,,,,,,,,,,,,,,'/
      data iopen /0/
c
c  test for opening scratch file
c
      if(iopen.eq.0) then
	open(99,status='scratch')
	iopen=1
      else
	rewind 99
      end if
c
      ierr=0
c
      read(ids,'(a)',end=80,err=90) sa
c
      acomma = sa//comma//comma//comma
c
      write(99,'(a)') acomma
      rewind 99
      read(99,*,end=80,err=90) (ia(i),i=1,n)
c
      return
c
   80 ierr=1
      return
c
   90 ierr=2
      return
      end
