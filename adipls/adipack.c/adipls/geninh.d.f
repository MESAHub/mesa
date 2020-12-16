      subroutine geninh(x,y,rhs,ii,iy,ig,isn,nd1,nn,mdintg,iclcn)
c  driving subroutine for initial value integration of ordinary
c  differential equations. method depends on value of mdintg:
c   mdintg = 1: second-order finite differences
c   mdintg = 2: constant coefficient method (only for ii = 2, ig = 1)
c   mdintg = 5: fourth-order finite differences
c  for mdintg = 2 iclcn is set to contribution to 'eckart'
c  classification of solution from interval considered, otherwise
c  iclcn is returned as 0.
c
c  Modified 29/3/05 adding fourth-order integration method
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      dimension x(nn),y(iy,nn)
      common/clscon/ icl
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      external rhs
c
      go to (10,20,30,40,50), mdintg
c
   10 call lininh(x,y,rhs,ii,iy,ig,isn,nd1,nn)
      iclcn=0
      return
c
   20 call eiginh(x,y,rhs,ii,iy,ig,isn,nd1,nn)
      iclcn=icl
      return
c
c
   30 write(istdou,100) 
      stop
   40 call eigin4(x,y,rhs,ii,iy,ig,isn,nd1,nn)
      iclcn=0
      return
c
   50 call lininh4(x,y,rhs,ii,iy,ig,isn,nd1,nn)
      iclcn=0
      return
  100 format(//' ***** error. geninh called with mdintg = 3')
      end
      subroutine eigin4(x,y,rhs,ii,iy,ig,isn,nd1,nn)
c
c  dummy routine
c
      implicit double precision (a-h,o-z)
      external rhs
      return
      end
