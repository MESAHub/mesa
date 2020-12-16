C SUBROUTINE TO CALCULATE THE SECOND ORDER DISTORTION TERMS WHICH ARE 
C DEPENDENT ON m^2, (AS IN KJELDSEN ET AL, 1998)
C AUTHOR KARA BURKE, JULY 2006

      subroutine delta2(x,y,del2,del2as,xis,aa,omgrtp,sig,iy,ia,nn,el)
       
      implicit real*8 (a-h,o-z)
      integer*4 v
      include 'adipls.c.d.incl'
      parameter(iwork=10*nnmax)
      real mass, rad, gconst, oms2, f4, f4int, oms

      dimension aa(ia,1)
      dimension omgrtp(*),dom(nnmax)
      dimension x(1), y(iy,1)
      dimension rho(nnmax)
      dimension xi(nnmax), eta(nnmax)
      dimension xis(2,nnmax)
      dimension f1(nnmax), f2(nnmax), f3(nnmax)
      dimension f1int(nnmax), f2int(nnmax), f3int(nnmax)
      dimension ci(nnmax), cint(nnmax)
      dimension del2(nnmax), del2as(nnmax)

      dimension f3a(nnmax), f3b(nnmax),f3aint(nnmax), f3bint(nnmax)

      common/crot_split/ beta, split(1)
      common/work/work(iwork)
      common/csumma/ cs(50)
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

c
c Set conversion factor      
c
      oms2 = gconst*mass/(rad*rad*rad)
      oms=sqrt(oms2)

      call derive(x,omgrtp,dom,nn,1,1,1,1)
c
c set rho, xi, eta
c
      do 500 n=1,nn
         rho(n) = aa(1,n)*aa(5,n)/(4.*pi)         
         xi(n) = y(1,n)
         eta(n) = y(2,n)/ell

 500  continue
c
c set integrand for I
c
      do 600 n=1,nn     
         ci(n)=rho(n)*(xi(n)**2.+eta(n)**2.*ell)*x(n)**2.
 600  continue

      if(dom(nn).eq.0.) then
      do 700 n=1,nn
c
c 1) Coefficient of m^2 from xi0.N0xi0
c
      f1(n)=1./2.*(-xi(n)*xi(n)+eta(n)*eta(n)*(2.-ell)+4.*eta(n)*xi(n))
     r        *rho(n)*x(n)*x(n)
c
c 2) Coefficient of m^2 from xi
c
      f2(n)=-(xi(n)*xis(1,n)+ell*eta(n)*xis(2,n))*x(n)*x(n)*rho(n)

c
c 3) Coefficient of m^2 from xi0.M0xi1
c
c in this term xis and etas have a factor of Omega in them so the second half 
c of the term is therefore also proportional to Omega^2 and this factor is 
c taken out.

      f3a(n) = -4.*rho(n)*((el+2.)*(xi(n)-el*eta(n))**2./(2.*el+1.)/
     1 (2.*el+3.)/(el+1.) + (el-1.)*(xi(n)+(el+1.)*eta(n))**2./el/
     2 (2.*el+1.)/(2.*el-1.))*x(n)*x(n)/freq
     
      f3b(n) = 2. *rho(n)*(xi(n)*xis(1,n)+ell*eta(n)*xis(2,n) - 
     4  eta(n)*xis(1,n) - xis(2,n)*xi(n)-eta(n)*xis(2,n))*x(n)*x(n)

c
c 4) Coeffiecient of m^2 from omega1 terms
c
      f4=beta*beta/2.

 700  continue

      else
      do 800 n=1,nn
c
c 1) Coefficient of m^2 from xi0.N0xi0
c
      f1(n)=(xi(n)*xi(n)*(4.*omgrtp(n)*dom(n)*x(n)/(4.*ell-3.)-omgrtp(n)
     1* omgrtp(n))
     1     + eta(n)*eta(n)*(2.-ell)*omgrtp(n)*omgrtp(n)
     2    +eta(n)*xi(n)*(4.*omgrtp(n)*omgrtp(n)+x(n)*omgrtp(n)*dom(n)*6.
     2        /(4.*ell-3.)))*rho(n)*x(n)*x(n)
c
c 2) Coefficient of m^2 from xi
c
      f2(n)=(xi(n)*xis(1,n)+ell*eta(n)*xis(2,n))*x(n)*x(n)*rho(n)
c
c 3) Coefficient of m^2 from xi0.M0xi1
c
      f3(n)=(2.*rho(n)*omgrtp(n)*(xi(n)*xis(1,n)+ell*eta(n)*xis(2,n)
     r     -eta(n)*xis(1,n)-xis(2,n)*xi(n)-eta(n)*xis(2,n))
     1     -4.*rho(n)*omgrtp(n)*omgrtp(n)*((el+2.)/((2.*el+1.)*(2.*el+3.
     1     )*(el+1.))*(xi(n)-el*eta(n))*(xi(n)-el*eta(n))
     2     +(el-1.)/(el*(2.*el+1.)*(2.*el-1.))*(xi(n)+(el+1.)*eta(n))*
     2     (xi(n)+(el+1.)*eta(n))))*x(n)*x(n)
c
c 4) Coefficient of m^2 from omega1 terms
c
      f4=beta*beta/2.
 800  continue
      end if
      
      call vinta(x,ci,cint,nn,1,1)
      call vinta(x,f1,f1int,nn,1,1)     
      call vinta(x,f2,f2int,nn,1,1)     
      call vinta(x,f3a,f3aint,nn,1,1)
      call vinta(x,f3b,f3bint,nn,1,1)



      if(dom(nn).eq.0.) then
      do 850 n=1,nn
         f1int(n)=f1int(n)/cint(n)
         f2int(n)=f2int(n)*beta/cint(n)*freq
         f3aint(n)=f3aint(n)/cint(n)/2.*freq
         f3bint(n)=f3bint(n)/cint(n)/2.*freq
         f4int=f4

         f3int(n) = f3aint(n) + f3bint(n)

 850     continue
      else
c
c need to check the omega not constant equations thoroughly
c

      do 900 n=1,nn  
         f1int(n)=f1int(n)/(2.*cint(n))*oms2/omgrtp(n)/omgrtp(n)
         f2int(n)=f2int(n)*beta/(2.*cint(n))*oms2/omgrtp(n)/omgrtp(n)
         f3int(n)=f3int(n)/(2.*cint(n))*oms2/omgrtp(n)/omgrtp(n)
         f4int=f4
 900  continue
      end if

c total coefficient of m^2
      do 950 n=1,nn
         del2(n) = f1int(n)+f2int(n)+f3int(n)+f4int
 950  continue
      if(abs(del2(nn)).ge.10) then
        write(istder,*) 'Excessive del2(nn) = ',del2(nn)
        write(91,*) el,sig,f4int,del2(nn)
        write(91,'(i5,1p6e12.4)') (n, x(n), f1int(n),f2int(n),
     *    f3int(n),cint(n),del2(n),n=1,nn)
        write(92,*) el, sig
        write(92,'(i5,1p5e12.4)') (n,x(n),xi(n),eta(n),xis(1,n),
     *    xis(2,n),n=1,nn)
      end if
c asymptotic approxiamtion

      do n=1,nn
         del2as(n) = 4./(2.*el-1.)/(2.*el+3.)
      end do
      
      return
      end
