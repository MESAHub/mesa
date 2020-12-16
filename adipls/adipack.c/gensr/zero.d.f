      subroutine  zero(a,nn)
c
c
c     sets a(n)=0., n=1,nn
c
c
c
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
      dimension a(nn)
      do 1 n=1,nn
    1 a(n)=0.d0
      return
c
      end
