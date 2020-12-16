      subroutine store(a,b,n)
c
c     stores first n elements of single precision a 
c        into single precision b
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90.
c
c  Note: this double precision version of the routine has same name
c  as single precision version, but is distinguished by its file name
c
      implicit double precision(a-h,o-z)
      dimension a(n),b(n)
c
   10 do 11 i=1,n
   11 b(i)=a(i)
      return
c
      end
