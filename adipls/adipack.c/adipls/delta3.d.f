C A SUBROUTINE TO CALCULATE DELTA3, AS DEFINED IN KJELDSEN 98 
C (WITH OMEGA not CONSTANT INSIDE THE STAR)


C THE EQUATIONS USED HERE NEED TO BE CHECKED THOROUGHLY AND REPROGRAMMED 
C CORRECTLY. CURRENTLY I AM USING EQUATIONS LIFTED FROM MJT. -- KDB 11/04/07
      subroutine delta3(x,y,del3,del3as,u,aa,omgrtp,sig,iy,ia,nn,el)


      implicit real*8 (a-h,o-z)
      integer*4 v
      include 'adipls.c.d.incl'
      parameter(iwork=10*nnmax)

      real mass, rad, gconst, oms2


      dimension aa(ia,1)
      dimension omgrtp(*)
      dimension x(1), y(iy,1)
      dimension u(2,nnmax)
      dimension rho(nnmax),dpr(nnmax)
      dimension xi(nnmax),eta(nnmax)
      dimension q(nnmax)
      dimension del3(nnmax), del3as(nnmax)
      dimension ci(nnmax),cint(nnmax)

      dimension uu(nnmax), sl(nnmax)

      dimension t1(nnmax), t2(nnmax), t3(nnmax), t4(nnmax)
      dimension t1int(nnmax), t2int(nnmax), t3int(nnmax), t4int(nnmax)

      common/work/work(iwork)
      common/csumma/ cs(50)

      iw=1

      pi = 4.d0*atan(1.d0)

      ell = el*(el+1.)

      freq = sqrt(sig)

      mass = cs(2)/1000
      rad = cs(3)/100
      gconst = 6.67e-11
c
c Set conversion factor
c      
      oms2 = gconst*mass/(rad*rad*rad)

      do 500 n=1,nn
c
c set pressure, rho, xi, eta
c
         dpr(n) = -aa(1,n)*aa(1,n)*aa(5,n)*x(n)/(4.*pi)
         rho(n) = aa(1,n)*aa(5,n)/(4.*pi) 
         xi(n) = y(1,n)
         eta(n) = y(2,n)/ell
c
c set q, U
c
         q(n) = aa(1,n)*x(n)**3.

         uu(n) = aa(4,n)/x(n)

 500  continue

      uu(1) = 0.

c
c set integrand for I
c
      do 650 n=1,nn     
         ci(n)=rho(n)*(xi(n)**2.+eta(n)**2.*ell)*x(n)**2.
 650  continue

C ALTERNATIVE USING EQUATIONS FROM MJT CODE


      do n=1,nn

      t1(n) = u(2,n)*(dpr(n)*uu(n) - rho(n) *freq*freq)*x(n)*x(n) 
     +     + 2.*u(1,n)*(dpr(n)*uu(n) + rho(n) *freq*freq)*x(n)

      t1(n) = t1(n) *xi(n) *xi(n) 

      t2(n) = u(2,n)*(ell*rho(n) - freq*freq*rho(n)*aa(2,n)/aa(1,n))
     8        *freq*freq
     +        - 2.*u(1,n)*freq**4.*rho(n)*aa(2,n)/aa(1,n)/x(n)

      t2(n) = t2(n)*eta(n)*eta(n)*x(n)*x(n)

      t3(n) = 2.*u(1,n)*(2.*dpr(n)*uu(n) + rho(n)*freq*freq)

      t3(n) = t3(n)*xi(n)*eta(n)*3.*x(n)

      t4(n) = u(2,n)*rho(n)*freq*freq + 
     +     u(1,n)*(2.*(aa(4,n)-aa(2,n))/x(n))*rho(n)*freq*freq

      t4(n) = t4(n) *eta(n)*eta(n)*x(n)*x(n)*3.

      end do

      t1(1) = 0.
      t2(1) = 0.
      t3(1) = 0.
      t4(1) = 0.

      call vinta(x,t1,t1int,nn,1,1)
      call vinta(x,t2,t2int,nn,1,1)
      call vinta(x,t3,t3int,nn,1,1)
      call vinta(x,t4,t4int,nn,1,1)

      call vinta(x,ci,cint,nn,1,1)

      open(76,file="del3test")
      do n=1,nn
         write(76,*) el, x(n), u(1,n), u(2,n)
      end do



      do n=1,nn
         del3(n) = t1int(n) + t2int(n) + t3int(n) + t4int(n)
         del3(n) =-del3(n)/freq/freq/2./cint(n)/omgrtp(n)/omgrtp(n)*oms2
      end do

c asymptotic approximation - for large values of n
      do n=1,nn
         sl(n) = aa(1,n)*x(n)**2./aa(2,n) 
         t1(n) = 4./3.*x(n)*x(n)*x(n)/q(n)*omgrtp(n)*omgrtp(n)/oms2
         t1(n) = t1(n)/sqrt(sl(n))
         t2(n) = 1./sqrt(sl(n))
      end do

      t1(1) = 0.
      t2(1) = 0.

      call vinta(x,t1,t1int,nn,1,1)
      call vinta(x,t2,t2int,nn,1,1)

      t1int(nn) = t1int(nn)/omgrtp(nn)/omgrtp(nn)*oms2
      del3as(nn) = t1int(nn)/t2int(nn)

      print*, "del3",el,cs(19), x(nn), del3as(nn), del3(nn)
       
      return
      end
