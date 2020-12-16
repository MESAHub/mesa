      subroutine iprsgn(ids,ia,nn,s)
c
c  prints the array ia(n), n = 1,...,nn, excluding possible
c  list of zeros at the end, and followed by the string s
c
c   
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 15/7/90.
c
c  Note: this double precision version of the routine has same name
c  as single precision version, but is distinguished by its file name
c
      character*(*) s
      dimension ia(nn)
c
      n1=nn+1
      do 10 n=1,nn
      n1=n1-1
      if(ia(n1).ne.0) go to 15
   10 continue
c
   15 write(ids,*) (ia(n),n=1,n1),s
c
      return
      end
      subroutine prtsgn(ids,a,nn,s)
c
c  prints the array a(n), n = 1,...,nn, excluding possible
c  list of zeros at the end, and followed by the string s
c
      implicit double precision(a-h,o-z)
      character*(*) s
      dimension a(nn)
c
      n1=nn+1
      do 10 n=1,nn
      n1=n1-1
      if(a(n1).ne.0) go to 15
   10 continue
c
   15 write(ids,*) (a(n),n=1,n1),s
c
      return
      end
