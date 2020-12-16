      double precision function danorm(a,n)
c  sets l2-norm of array a(n)
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      dimension a(n)
      common/cprcns/ epsprc,epsufl,epsofl,eprufl
c
      save
c
c  find maximum absolute element of a
      amx=0
      do 10 i=1,n
   10 amx=max(amx,abs(a(i)))
c
      if(amx.gt.0) go to 15
      danorm=0
      return
c
   15 sum=0
      do 20 i=1,n
      b=abs(a(i)/amx)
      if(b.lt.eprufl) b=0
   20 sum=sum+b*b
c
      danorm=sqrt(sum)*amx
      return
      end
c
c
