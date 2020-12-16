      subroutine eiginh(x,y,rhs,ii,iy,ig,isn,nd1,nn)
c
c  Initial value integrator.
c  =========================
c  Integrates a set of two first-order ordinary
c  linear homogeneous differential equations, by means of 
c  constant-coefficient technique (Gabriel & Noels, A and A, vol. 53,
c  p. 149, 1976).
c
c  The equations are on the form
c
c       d y(i; x)/dx = sum ( a(i,j; x) * y(j; x))
c                       j
c
c  The coefficients are set up in the routine rhs, which must be
c  supplied by the user as external.
c  
c  The routine finds the solution
c  at every isn-th point in the input mesh. For ordinary use 
c  isn is set to 1.
c  
c  nd1 specifies the initial value of the mesh index passed into the
c  routine rhs (see below).
c
c  iy is the first dimension of y in the calling programme.
c
c  The arguments ii and ig are included for consistency with
c  subroutine linint.
c  
c  On input x(1+isn*(n-1)), n = 1,nn, must contain the mesh in the 
c  independent variable. The initial conditions for set k must be
c  set into y(i,1), i = 1,ii.
c  
c  The routine returns the solution for initial value set k in
c  y(i,1+isn*(n-1)), i=1,ii, n=1,nn.
c  
c  The right hand side subroutine rhs:
c  -----------------------------------
c  
c  This is called by linint and should be defined as
c  
c        subroutine rhs(x,aa,finh,iaa,nd)
c        dimension aa(iaa,1)
c           .
c           .
c           .
c  
c  A single call must set the coefficient matrix at a single point.
c  On input x (a scalar) gives the value of the independent variable
c  at the given point; iaa is the first dimension of aa;
c  nd is passed from linint and may be used to address a data
c  array in rhs. 
c
c  finh may later be used to set inhomogeneous term
c  
c  The routine should return
c  
c  aa(i,j) = a(i,j; x)
c  
c  where a(i,j; x) as defined above is the coefficient matrix in
c  the equations.
c  
c  The definition of nd is slightly convoluted. It is assumed that the
c  use might need certain variables to define the right hand side
c  of the equations. Assume, for example, that these are passed into 
c  the routine in
c  
c        common/rhsdat/ data(10,200)
c  
c  When rhs is called at the meshpoint x(1+isn*(n-1)), nd is
c  set to nd1+isn*(n-1). The meshpoint should therefore correspond
c  to the data in data(i,nd). In most cases presumably x and the
c  data would be given on the same mesh, and nd1 would be 1.
c   
c  Note: this version is consistent with s/r lininh to integrate
c  inhomogeneous equations, in the call of rhs.
c  However, a non-zero inhomogeneous terms has not yet been
c  implemented.
c
c  Original version: 9/7/95
c
      implicit double precision(a-h,o-z)
      common/clscon/ icls
c
      dimension x(nn),y(iy,nn),f(2,4),finh(2),cc(2,20),al(2),a(2,2)
      external rhs
      ifd=2
      icls=0
c  right hand sides at first point
    5 n=1
      nc=1
      nd=nd1
      if1=1
      if2=3
      x2=x(1)
      call rhs(x2,f(1,if1),finh,ifd,nd)
c  right hand sides at next point
   10 nd=nd+isn
      nc=nc+1
      n1=n
      n=n+isn
      x1=x2
      x2=x(n)
      dx=x2-x1
c
      call rhs(x2,f(1,if2),finh,ifd,nd)
c  set mean coefficient matrix
      j1=if1-1
      j2=if2-1
      do 15 i=1,2
      do 15 j=1,2
   15 a(i,j)=0.5d0*(f(i,j1+j)+f(i,j2+j))
c
      do 17 i=1,2
   17 y(i,n)=y(i,n1)
c  integrate constant coefficient equation
      call ieigst(a,y(1,n),dx,ig,al,cc,ie)
c  increment icls by contribution from this interval
      icls=icls+iclcon(x2,a,dx,al,cc,ie)
c..      write(72,72090) n, x(n), y(1,n), y(2,n), icls
c..72090 format(i5,f12.7,1p2e13.5,i5)
c
      if(nc.eq.nn) return
      i=if1
      if1=if2
      if2=i
c
      go to 10
c
      end
      subroutine ieigst(a,yy,tt,ig,al,cc,ie)
c  integrates second order constant coefficient equation from 0
c  to tt
      implicit double precision(a-h,o-z)
      dimension a(2,2),yy(1),al(2),cc(2,1),x(2,2)
c  eigenvalues of a
      do 12 i=1,2
      al(i)=0
      do 12 j=1,2
   12 x(i,j)=0
      call egenv2(a,x,al,ie)
c  set constants in analytical solution and solution at tt
      dd=x(1,1)*x(2,2)-x(1,2)*x(2,1)
c
      i1=-1
      do 30 k=1,ig
      i1=i1+2
      i2=i1+1
      i0=i1-1
      y1=yy(i1)
      y2=yy(i2)
      ca=(x(2,2)*y1-x(2,1)*y2)/dd
      cb=(x(1,1)*y2-x(1,2)*y1)/dd
c
      if(ie-1) 14,17,22
   14 do 15 i=1,2
      cc(i,i1)=ca*x(1,i)+cb*x(2,i)
   15 cc(i,i2)=cb*x(1,i)-ca*x(2,i)
      amt=al(2)*tt
      elt=exp(al(1)*tt)
      c=cos(amt)*elt
      s=sin(amt)*elt
      do 16 i=1,2
   16 yy(i+i0)=cc(i,i1)*c+cc(i,i2)*s
      go to 30
c
   17 all=0.5d0*(a(1,1)-a(2,2))
      if(x(2,1).eq.0) go to 18
      cc(1,i2)=all
      cc(2,i2)=a(2,1)
      go to 19
   18 cc(2,i2)=-all
      cc(1,i2)=a(1,2)
   19 do 20 i=1,2
      cc(i,i1)=ca*x(1,i)+cb*x(2,i)
   20 cc(i,i2)=cb*cc(i,i2)
      elt=exp(al(1)*tt)
      do 21 i=1,2
   21 yy(i+i0)=(cc(i,i1)+tt*cc(i,i2))*elt
      go to 30
c
   22 do 23 i=1,2
      cc(i,i1)=ca*x(1,i)
   23 cc(i,i2)=cb*x(2,i)
      el1=exp(al(1)*tt)
      el2=exp(al(2)*tt)
      do 24 i=1,2
   24 yy(i+i0)=cc(i,i1)*el1+cc(i,i2)*el2
   30 continue
      return
      end
      subroutine egenv2(a,x,al,i)
c  finds eigenvectors and eigenvalues of real 2*2 matrix a.
c   output:
c   ======
c  for i.ge.1: i (1 or 2) real eigenvalues corresponding to different
c  eigenvectors (the eigenvalues may be the same). the j-th eigenvalue
c  is in al(j) and the j-th eigenvector in x(j,1),x(j,2), for j=1,i.
c
c  for i = -1: two complex conjugated eigenvalues. the eigenvalues are
c   al(1) +- i*al(2), and the eigenvectors have the k-th component
c   x(1,k) +- i*x(2,k), k=1,2.
      implicit double precision(a-h,o-z)
      dimension a(2,2),x(2,2),al(2)
      eps=1.d-10
      eps2=eps*eps
      alp=(a(1,1)+a(2,2))/2
      y=a(2,2)-a(1,1)
      x1=a(1,2)*a(2,1)
      if(abs(x1).lt.eps2)goto 2
      dis=y*y+4*x1
      if(abs(dis).gt.eps2) goto 6
c  dis = 0, one real eigenvalue
      i=1
      al(1)=alp
      x(1,1)=2*a(1,2)/y
      x(1,2)=1.d0
c  set x(2,.)
      if(abs(a(1,2)).gt.abs(a(2,1))) go to 9
      x(2,1)=1.d0
      x(2,2)=0.d0
      return
    9 x(2,1)=0.d0
      x(2,2)=1.d0
      return
c  dis.ne.0
    6 if(dis.lt.0.0d0) goto 1
c  two real eigenvalues
      dis=sqrt(dis)/2
      d=a(1,1)*a(2,2)-x1
      if(alp) 11,11,12
   11 al(2)=alp-dis
      al(1)=d/al(2)
      go to 14
   12 al(1)=alp+dis
      al(2)=d/al(1)
   14 do 16 j=1,2
      x(j,2)=1
   16 x(j,1)=a(1,2)/(al(j)-a(1,1))
      i=2
      return
c  two complex eigenvalues
    1 i=-1
      al(1)=alp
      al(2)=sqrt(-dis)/2
      xx=alp-a(1,1)
      xx1=xx*xx-dis/4
      x(1,1)=a(1,2)*xx/xx1
      x(1,2)=1.0d0
      x(2,1)=-a(1,2)*al(2)/xx1
      x(2,2)=0.0d0
      return
c  zero off-diagonal element
    2 if(abs(y).lt.eps) goto 3
c  two different real eigenvalues
      al(1)=a(1,1)
      al(2)=a(2,2)
      x(1,1)=1.d0
      x(2,2)=1.d0
      x(1,2)=-a(2,1)/y
      x(2,1)=a(1,2)/y
      i=2
      return
c  diagonal elements equal. one or two eigenvectors
    3 i=0
      if(abs(a(2,1)).gt.eps) goto 4
      i=1
      al(1)=a(1,1)
      x(1,1)=1.d0
      x(1,2)=0.d0
      if(abs(a(1,2)).lt.eps) goto 5
c  a(1,2).ne.0, set x(2,.)
      x(2,1)=0.d0
      x(2,2)=1.d0
      return
c  a(2,1).ne.0., set x(2,.)
    4 x(2,1)=1.d0
      x(2,2)=0.d0
c
    5 i=i+1
      al(i)=a(2,2)
      x(i,1)=0.d0
      x(i,2)=1.d0
      return
      end
      integer function iclcon(x,a,tt,al,cc,ie)
c  finds contribution to classification index from interval (0,tt),
c  from solution of constant coefficient equation.
      implicit double precision(a-h,o-z)
      dimension a(2,2),al(2),cc(2,*)
      data pi/3.14159265358979d0/
c  find number of zeros of y1
      c1=cc(1,1)
      c2=cc(1,2)
      icl=0
      if(c1.eq.0.and.c2.eq.0) go to 50
c
   10 if(ie-1) 20,30,40
c  oscillatory solution
   20 if(al(2).eq.0) go to 50
      delta=atan2(c2,c1)
      zt0=-delta/pi-0.5d0
      ztt=al(2)*tt/pi+zt0
c  test for direction of integration
      if(tt) 22,50,24
   22 zt1=ztt
      zt2=zt0
      go to 26
c
   24 zt1=zt0
      zt2=ztt
c  now zt2 corresponds to the larger x, zt1 to the smaller
   26 izt1=intgpt(-zt1)
      izt2=intgpt(zt2)
      icl=1+izt1+izt2
c  test for zero at smallest x
      if(izt1.eq.-zt1) icl=icl-1
      go to 50
c  degenerate and exponential cases
c
c  find zero
   30 if(c2.eq.0) go to 50
      tz=-c1/c2
      go to 45
c
   40 if(c1*c2.ge.0.or.al(1).eq.al(2)) go to 50
      tz=log(-c2/c1)/(al(1)-al(2))
c  test for direction of integration
   45 if(tt) 46,50,48
c  test for position of zero
   46 if(tt.lt.tz.and.tz.le.0) icl=1
      go to 50
   48 if(0.lt.tz.and.tz.le.tt) icl=1
c  find sign of a(1,2)
   50 isg=0
      if(a(1,2).ne.0) isg=sign(1.d0,a(1,2))
c  contribution to index
      iclcon=-isg*icl
c..      write(73,73090) x,zt1,zt2,izt1,izt2,icl
c..73090 format(f12.7,1p2e13.5,3i10)
      return
      end
