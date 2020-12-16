C PROGRAM TO FIND xi_s and eta_s, USING GOUGH AND THOMPSON 1990  - 
C CORRECTED EQUATIONS. 
C AUTHOR KARA BURKE, JULY 2006

      subroutine sph(x,y,xis,aa,omgrtp,el,sig,iy,ia,nn)

      implicit real*8 (a-h,o-z)
      integer*4 v
      include 'adipls.c.d.incl'
      parameter (iw=10*nnmax)

      real mass, rad, gconst, oms, oms2

      dimension aa(ia,nn),omgrtp(*)
      dimension xis(2,nnmax)
      dimension dom(nnmax),d2om(nnmax), beta(nnmax)
      dimension aq(20,nnmax), om1(1)
      dimension x(nn),y(iy,nnmax)
      dimension ea(2,3),v(2), zk(1), ap(1)
      dimension h(nnmax)
      dimension q(nnmax),rho(nnmax), rhoc(nnmax)
      dimension xi(nnmax),eta(nnmax),deta(nnmax)
      dimension dpr(nnmax),dlnpr(nnmax),drho(nnmax),sl(nnmax),uu(nnmax)
      dimension etaxi(nnmax),detaxi(nnmax), dxi(nnmax)

      common/nnpts/ nn1
      common/work/work(iw)      
      common/crot_split/ beta, split(4)

      external rhskdbsph,bckdbsph

      nn1 = nn-1

c calculate domega/dx, d2omega/dx2


      call derive(x,omgrtp,dom,nn,1,1,1,1)
      call derive(x,dom,d2om,nn,1,1,1,1)

      pi = 4.d0*atan(1.d0)
      ell = el*(el+1.)
      freq = sqrt(sig)
      gconst = 6.67e-11

      do n=1,nn

         q(n) = aa(1,n)*x(n)**3.
         rho(n) = aa(1,n)*aa(5,n)/(4.*pi)

         xi(n) = y(1,n)
         eta(n) = y(2,n)/ell
	 if(n.le.2.or.n.eq.nn) 
     *     write(6,*) 'in sph n, ell, y(2,n), eta(n) =',
     *     n, ell, y(2,n), eta(n)

         dpr(n) = -aa(1,n)*aa(1,n)*aa(5,n)*x(n)/(4.*pi)
         drho(n) = -aa(1,n)*aa(5,n)*(aa(2,n)+aa(4,n))/(4.*pi*x(n))
         sl(n) = aa(1,n)*x(n)**2./aa(2,n)
         dlnpr(n) = - aa(2,n)/x(n)
         uu(n) = aa(4,n)/x(n)

         aq(1,n) = q(n)
         aq(2,n) = rho(n)
         aq(3,n) = dom(n)
         aq(4,n) = d2om(n)
         aq(5,n) = omgrtp(n)
         aq(6,n) = xi(n)
         aq(7,n) = eta(n)
         aq(8,n) = dpr(n)
         aq(9,n) = drho(n)
         aq(10,n) = sl(n)
         aq(11,n) = dlnpr(n)
         aq(12,n) = uu(n)
         aq(13,n) = freq
         aq(14,n) = el
c         aq(16,n) = split(n)


         etaxi(n) = eta(n) + xi(n)

      end do

      call derive(x,etaxi,detaxi,nn,1,1,1,1)
      call derive(x,eta,deta,nn,1,1,1,1)
      call derive(x,xi,dxi,nn,1,1,1,1)

      do n=1,nn
         aq(15,n) = detaxi(n)
         aq(17,n) = deta(n)
         aq(18,n) = dxi(n)
      end do

      ii=2
      id=2
      kk=1
      ka=1
      kb=2
      ki=0

      v(1) = 1
      v(2) = 2
      ucy = 1.
 

      iter = 0

      do 130 l1=1,2
      do 130 l2 = 1,nn
 130     xis(l1,l2)= 0.0

      zk(1)= 1.

  150 continue

      write(*,*) 'Before call to nrk in sph'

      call nrk(x(2),xis(1,2),zk,ap,aq(1,2),rhskdbsph,bckdbsph,ii,kk,ka
     *,kb,ki,nn-1,id,ucy,ea,det,v)

      write(*,*) 'After call to nrk in sph'
      write(*,*) 'ea:'
      write(*,*) ea
      iter = iter+1
      if (iter.lt.3) goto 150

      xis(1,1) = xis(1,2)
      xis(2,1) = xis(2,2)

c..      write(0,*) '#D# returning from sph'
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c RHS
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      subroutine rhskdbsph(x,y,zk,ap,aq,f,fd,h,hd,ifd,n)

      implicit real*8 (a-h,o-z)

      include 'adipls.c.d.incl'
      dimension y(*),f(*),fd(ifd,*),h(*),hd(ifd,*),zk(1),ap(*)
      dimension aq(20,*)

      common/csumma/ cs(50)
      
C RHS FOR CALCULATING xi_s and eta_s      

      pi = 4.d0*atan(1.d0)

      q = aq(1,n)
      rho = aq(2,n)
      dom = aq(3,n)
      d2om = aq(4,n)
      omgrtp = 1.
      xi = aq(6,n)
      eta = aq(7,n)
      dpr = aq(8,n) 
      drho = aq(9,n) 
      sl = aq(10,n)
      dlnpr = aq(11,n)
      uu = aq(12,n)
      freq = aq(13,n)
      el = aq(14,n)
      detaxi = aq(15,n)
      deta = aq(17,n)
      dxi = aq(18,n)
      
      ell = el*(el+1.)

c defining p

      p1 = 2.*omgrtp*x*freq/sl *(xi+eta)/ell
      p2 = 2.*x*freq/sl *omgrtp*eta

      p = (-p1 + p2)

c defining s

      s1 = -2.*omgrtp/freq*(xi/x + (uu-1./x)*eta - deta)
      s2 = 2.*omgrtp/(freq*ell)*(ell*eta/x +(uu - 1./x)*(eta+xi)-deta-
     *     dxi)
      s3 = 2.*dom/freq*(eta-(eta+xi)/ell)
      
      s = (s1 + s2 + s3)

c      define rhs

      fd(1,1) = -(2./x + drho/rho + uu)
      fd(1,2) = ell/x - freq*freq * x/sl
c
c for comparisons sake have changed from freq to freq^2
c (see MJT code)
c ALSO CHECK THIS EQUATION!
c
      fd(1,3) = -2.*x*freq*freq/sl*eta
      f(1) = fd(1,1)*y(1) + fd(1,2)*y(2) + fd(1,3)*zk(1) + p
      fd(2,1) = 1./x + dpr*uu/(rho*freq**2.*x)
      fd(2,2) = -(1./x - uu)
c
c again for comparisons sake have changed to freq^2
c

      fd(2,3) = -2.*uu*dpr/rho/freq/freq/x*xi

      f(2) = fd(2,1)*y(1) + fd(2,2)*y(2) + fd(2,3)*zk(1) + s


      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c BOUNDARY CONDITIONS
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bckdbsph(x1,x2,y1,y2,zk,ap,aq,g,gd,ig,id,n)
      
      implicit real*8(a-h,o-z)

      include 'adipls.c.d.incl'
      dimension y1(*),y2(*),g(*),gd(ig,*),zk(1),ap(*),aq(20,*)
      common/nnpts/ nn1
      common/csumma/ cs(50)

C BOUNDARY CONDITIONS FOR CALCULATING xi_s and eta_s

      nn = nn1

      pi = 4.d0*atan(1.d0)

      l=2

      q1 = aq(1,l)
      q2 = aq(1,nn)
      rho1 = aq(2,l)
      rho2 = aq(2,nn)
      dom1 = aq(3,l)
      dom2 = aq(3,nn)
      d2om1 = aq(4,l)
      d2om2 = aq(4,nn)

      omgrtp1 = 1.
      omgrtp2 = 1.
   
      xi1 = aq(6,l)
      xi2 = aq(6,nn)
      eta1 = aq(7,l)
      eta2 = aq(7,nn)
      dpr1 = aq(8,l)
      dpr2 = aq(8,nn)
      drho1 = aq(9,l)
      drho2 = aq(9,nn)
      sl1 =  aq(10,l)
      sl2 =  aq(10,nn)
      dlnpr1 =  aq(11,l)
      dlnpr2 =  aq(11,nn)
      uu1 =  aq(12,l)
      uu2 =  aq(12,nn)

      freq = aq(13,1)
      el = aq(14,1)

      ell = el*(el+1.)
     
      
      t = 2.*omgrtp1*freq*x1/ell/sl1*(ell*eta1-xi1-eta1)
      p = -2.*omgrtp2*x2*freq/sl2/ell*(xi2+eta2-ell*eta2)

      gd(1,1) = -(uu1 + 2./x1 + drho1/rho1) - (el-1.)/x1
      gd(1,2) = -x1*freq*freq/sl1 + ell/x1
      gd(1,3) = -2.*x1*freq*freq*eta1/sl1
      g(1) = gd(1,1)*y1(1) + gd(1,2)*y1(2) + gd(1,3)*zk(1) + t
      gd(2,1) = 1.
      gd(2,2) = 0.
      gd(2,3) = 0.
      g(2) = gd(2,1)*y2(1) + gd(2,2)*y2(2) + gd(2,3)*zk(1)
      gd(3,1) = -uu2 - drho2/rho2
      gd(3,2) = -x2*freq*freq/sl2
      gd(3,3) = -2.*x2*freq*freq/sl2*eta2
      g(3) = gd(3,1)*y2(1) + gd(3,2)*y2(2) + gd(3,3)*zk(1) + p

      write(6,*) 'in bckdbsph x1,freq,eta1,sl1,sl2,eta2 =',
     *  x1,freq,eta1,sl1,sl2,eta2
      write(6,*) 'gd in bckdbsph:'
      write(6,'(1p3e15.7)') ((gd(i,j),i=1,3),j=1,3)

      return
      end
