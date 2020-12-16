C SUBROUTINE TO CALCULATE THE SECOND ORDER DISTORTION TERMS WHICH ARE 
C INDEPENDENT OF m, (AS IN KJELDSEN ET AL, 1998)
C AUTHOR KARA BURKE, JULY 2006

      subroutine delta1(x,y,del1,del1as,aa,omgrtp,sig,iy,ia,nn,el)

      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER*4 v
      include 'adipls.c.d.incl'
      parameter(iwork=10*nnmax)
      REAL mass, rad, gconst, oms2, oms

      DIMENSION aa(ia,1)
      DIMENSION omgrtp(*),dom(nnmax), omega(nnmax)
      DIMENSION x(1), y(iy,1)
      DIMENSION rho(nnmax)
      DIMENSION xi(nnmax), eta(nnmax)
      DIMENSION f1(nnmax), f2(nnmax), f1int(nnmax),f2int(nnmax)
      DIMENSION ci(nnmax), cint(nnmax)
      DIMENSION del1(nnmax), del1as(nnmax)

      COMMON/work/work(iwork)
      COMMON/csumma/ cs(50)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder

      pi = 4.d0*atan(1.d0)

      ell = el*(el+1.)

      freq = sqrt(cs(20))

      mass = cs(2)/1000
      rad = cs(3)/100
      gconst = 6.67e-11
      
      oms2 = gconst*mass/(rad*rad*rad)
      oms = sqrt(oms2)

      do n=1,nn
      omega(n) = omgrtp(n)/sqrt(oms2)
      end do

      call derive(x,omega,dom,nn,1,1,1,1)


c
c set rho, xi, eta
c
      do 500 n=1,nn
         rho(n) = aa(1,n)*aa(5,n)/(4.*pi)         
         xi(n) = y(1,n)
 500  continue

c
c set integrand for I
c

      if(el.gt.1.e-6) then

      do 600 n=1,nn
           eta(n) = y(2,n)/ell
           ci(n)=rho(n)*(xi(n)**2.+eta(n)**2.*ell)*x(n)**2.
c
c 1) Coefficient of m^0 from xi0.N0xi0
c
      f1(n)=rho(n)*omega(n)*dom(n)*x(n)*x(n)*
     *           (xi(n)*xi(n)*x(n)*(ell-1.)/(4.*ell-1.)
     2           -xi(n)*eta(n)*x(n)*ell/(4.*ell-3.))

c 2) Coefficient of m^0 from xi0.M0xi1
c
      f2(n)=2.*rho(n)*x(n)*x(n)*omega(n)*omega(n)/(2.*el+1.)*
     *     ((el+1.)*(el+2.)/(2.*el+3.)*(xi(n)-el*eta(n))**2.
     2   +el*(el-1.)/(2.*el-1.)*(xi(n)+(el+1.)*eta(n))**2.)

 600  continue
      
      else 
c equations for l=0
      do 800 n=1,nn
         ci(n)= rho(n)*xi(n)**2.*x(n)**2.
         f1(n) = 0.
         f2(n) = 4./3.*rho(n)*xi(n)**2.*x(n)**2.*omega(n)*omega(n)
c..         write(6,*) n, x(n), rho(n), xi(n), omega(n)
 800     continue
       end if


      call vinta(x,ci,cint,nn,1,1)
      call vinta(x,f1,f1int,nn,1,1)
      call vinta(x,f2,f2int,nn,1,1)


      write(istdpr,*) '#D# cint(nn) etc.',cint(nn), f1int(nn),f2int(nn)
      do 900 n=1,nn
         f1int(n)=f1int(n)/cint(n)/omega(n)/omega(n)
         f2int(n)=f2int(n)/cint(n)/omega(n)/omega(n)
 900  continue
 
c total coefficient of m^0

      do 950 n=1,nn
         del1(n) = f1int(n)+f2int(n)
 950  continue

 
c asymptotic approximation - for large values of n

      do n =1,nn
      del1as(n) = 2./(2.*el + 1.)*((el+1.)*(el+2.)/(2.*el+3.)
     +     + el*(el-1.)/(2.*el-1.))
      end do

      return
      end
