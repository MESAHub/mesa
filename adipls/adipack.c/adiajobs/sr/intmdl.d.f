      subroutine intmdl(x,nn,x1,aa,nn1,ia,w)
c
c  interpolates model given in (x1(n),aa(i,n),n=1,nn1), to mesh
c  (x(n),n=1,nn)
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 14/5/90
c
      implicit double precision (a-h, o-z)
      dimension x(nn),x1(nn),aa(ia,nn),w(5,nn)
      int=0
      n1=1
      do 20 n=1,nn
      if(x(n).eq.x1(n)) go to 20
      call lir(x(n),x1,w(1,n),aa,5,ia,nn1,n1,inter)
      n1=2
      int=1
   20 continue
c  test for interpolation
      if(int.eq.0) return
c  reset x1,aa
      do 30 n=1,nn
      x1(n)=x(n)
c
c  note that x = 0 is already ok
c
      if(x(n).gt.0) then
        do 25 i=1,5
   25   aa(i,n)=w(i,n)
      end if
c
   30 continue
      return
      end
