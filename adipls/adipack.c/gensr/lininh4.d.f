      subroutine lininh4(x,y,rhs,ii,iy,ig,isn,nd1,nn)
c
c
c  Initial value integrator.
c  =========================
c  Integrates a set of first-order ordinary
c  linear inhomogeneous differential equations, by means of second-order
c  centred difference approximation, from given initial values.
c
c  The equations are on the form
c
c       d y(i; x)/dx = sum ( a(i,j; x) * y(j; x)) + f(i; x)
c                       j
c
c  The order of the equations is ii.
c  
c  The coefficients and the inhomogeneous terms
c  are set up in the routine rhs, which must be
c  supplied by the user as external.
c  
c  The routine has the option for determining in a single call the
c  solutions to a given set of equations for several initial conditions;
c  ig specifies the number of such sets.
c  
c  For use with Richardson extrapolation the routine finds the solution
c  at every isn-th point in the input mesh. For ordinary use 
c  isn is set to 1.
c  
c  nd1 specifies the initial value of the mesh index passed into the
c  routine rhs (see below).
c
c  iy is the first dimension of y in the calling programme.
c  
c  On input x(1+isn*(n-1)), n = 1,nn, must contain the mesh in the 
c  independent variable. The initial conditions for set k must be
c  set into y(i + ii*(k-1),1), i = 1,ii, for k = 1,ig.
c  
c  The routine returns the solution for initial value set k in
c  y(i+ii*(k-1),1+isn*(n-1)), i=1,ii, n=1,nn, k=1,ig.
c  
c  The right hand side subroutine rhs:
c  -----------------------------------
c  
c  This is called by lininh and should be defined as
c  
c        subroutine rhs(x,aa,finh,iaa,nd)
c        dimension aa(iaa,1),finh(1)
c           .
c           .
c           .
c  
c  A single call must set the coefficient matrix and inhomogeneous
c  term at a single point.
c  On input x (a scalar) gives the value of the independent variable
c  at the given point; iaa is the first dimension of aa;
c  nd is passed from lininh and may be used to address a data
c  array in rhs. 
c  
c  The routine should return
c  
c  aa(i,j) = a(i,j; x)
c  finh(i + ii*(k-1)) = f(i; x) for the k-th set, k = 1,ig.
c  
c  where a(i,j; x) and f(i; x) as defined above are the 
c  coefficient matrix and inhomogeneous term in the equations.
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
c  4th order version, using algorithm of Cash & Moore (1980; BIT 20:1, 44 - 52)
c   
c  original version: 11/3/05
c
      implicit double precision(a-h,o-z)
      dimension x(nn),y(iy,nn),fd(20,20),fdp(20,20),fdh(20,20),finh(20),
     *  fp(20),w(20,20),y1(20)
      common/modfac/ r21,ntld
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      external rhs
c
c..      write(6,*) 'Enter lininh with ii, iy, ig =',ii, ig, iy
c..      write(6,*) 'y(1-2,1) =',(y(i,1),i=1,2)
c
c  explicitly disable check for linear dependence
c
      ifd=20
c
      iw=20
c  right hand sides at first point
    5 n=1
      nc=1
      nd=nd1
      ntst=nn-2
      iig=ii*ig
      x2=x(1)
      call rhs(x2,fd,finh,ifd,nd)
c
c  right hand sides at next point
c
   10 do i=1,ii
	do j=1,ii
	  fdp(i,j)=fd(i,j)
        end do
      end do
c
c  set up right hand side of linear equation at next point
c
      n1=n
      n=n+isn
      nc=nc+1
      x1=x2
      x2=x(n)
c
      nd=nd+isn
c
      call rhs(x2,fd,finh,ifd,nd)
c
      do i=1,ii
	do j=1,ii
	  fdh(i,j)=0.5d0*(fd(i,j)+fdp(i,j))
        end do
      end do
      dx=x2-x1
      kg=-ii
      do 25 l=1,ig
      kg=kg+ii
      do 23 i=1,ii
      k=kg+i
      sum=0
      do 22 j=1,ii
      sum=sum+fdp(i,j)*y(kg+j,n1)
   22 continue
      fp(i)=sum
   23 continue
c
      do i=1,ii
	sum1=0
	sum2=0
	do j=1,ii
	  sum1=sum1+fdh(i,j)*y(kg+j,n1)
	  sum2=sum2+fdh(i,j)*fp(j)
        end do
	y(kg+i,n)=y(kg+i,n1)+dx*fp(i)/6.d0+dx*sum1/3.d0
     *            +dx*dx*sum2/12.d0
      end do
c
   25 continue
c  set coefficient matrix for linear equations
      do i=1,ii
        do j=1,ii
	  sum=0
	  do k=1,ii
	    sum=sum+fdh(i,k)*fd(k,j)
          end do
	  w(i,j)=-(fd(i,j)+2.d0*fdh(i,j))*dx/6.d0+sum*dx*dx/12.d0
        end do
        w(i,i)=1.d0+w(i,i)
      end do
c
c  include inhomogeneous contribution from next point (needs to be added)
c
c..      do 15 i=1,iig
c..   15 y(i,n)=y(i,n)-dx*finh(i)
c
c..      write(6,15091) n,x(n),((w(i,j),j=1,2),y(i,n),i=1,2)
c..15091 format(' Equations at n, x =',i5,f10.5/(1p3e13.5))
c
c  solve linear equations by gaussian elimination without pivoting
c
      ii1=ii-1
      if(ii1.gt.0) then
c
c  triangularize matrix
c
        do 27 i=1,ii1
        r=w(i,i)
c
c  test for non-zero diagonal element
c
        if(r.eq.0) then
          if(istdpr.gt.0) write(istdpr,110) i,n
          return
        end if
c
        r=-1.d0/r
c
        i1=i+1
c
        do 27 j=i1,ii
        rj=r*w(j,i)
        do 26 k=i1,ii
   26   w(j,k)=w(j,k)+rj*w(i,k)
c
        js=j
        is=i
c
        do 27 l=1,ig
        y(js,n)=y(js,n)+rj*y(is,n)
        js=js+ii
   27   is=is+ii
c
      end if
c
c  now matrix is on triangular form
c
c  start solution
c
      i=ii+1
      do 28 icnt=1,ii
      i=i-1
c
      r=w(i,i)
c
c  test for non-zero diagonal element
c
      if(r.eq.0) then
        if(istdpr.gt.0) write(istdpr,110) i,n
        return
      end if
c
      r=1.d0/r
      is=i
      i1=i+1
c
      do 28 l=1,ig
      sum=y(is,n)
c
      if(i.lt.ii) then
c
        j1=is
        do j=i1,ii
          j1=j1+1
          sum=sum-w(i,j)*y(j1,n)
        end do
c
      end if
c
      y(is,n)=r*sum
c
   28 is=is+ii
      if(nc.eq.nn) return
c
      go to 10
c
  110 format(//'  ***** in s/r lininh zero diagonal element at i =',
     *  i5,'   n =',i5)
      end
