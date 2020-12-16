C PROGRAM TO FIND u(x)=h(x)x, USING GOUGH AND THOMPSON 1990 (3.25)-(3.27) -
C CORRECTED EQUATIONS.
C AUTHOR KARA BURKE, MAY 2006

      subroutine uhx(x,u,aa,omgrtp,ia,nn)

      implicit real*8 (a-h,o-z)
      integer*4 v
      include 'adipls.c.d.incl'
      parameter(iwork=10*nnmax)

      dimension aa(ia,nn),omgrtp(*)
      dimension dom(nnmax),d2om(nnmax)
      dimension aq(10,nnmax)
      dimension x(nn),u(2,nnmax)
      dimension ea(2,3),v(2)
      dimension h(nnmax)
      dimension q(nnmax),rho(nnmax)

      common/nnpts/ nn1
      common/work/work(iwork)

      external rhskdbh,bckdbh

      nn1 = nn-1

c calculate domega/dx, d2omega/dx2

      call derive(x,omgrtp,dom,nn,1,1,1,1)
      call derive(x,dom,d2om,nn,1,1,1,1)

      pi = 4.d0*atan(1.d0)

      do n=1,nn

         q(n) = aa(1,n)*x(n)**3.
         rho(n) = aa(1,n)*aa(5,n)/(4.*pi)

         aq(1,n) = q(n)
         aq(2,n) = rho(n)
         aq(3,n) = dom(n)
         aq(4,n) = d2om(n)
         aq(5,n) = omgrtp(n)


      end do

      ii=2
      id=2
      kk=0
      ka=1
      kb=1
      ki=0

      v(1) = 1
      v(2) = 2
      ucy = 1.


      iter = 0
  150 continue

      write(*,*) 'Before call to nrk in uhx'

      call nrk(x(2),u(1,2),zk,ap,aq(1,2),rhskdbh,bckdbh,ii,kk,ka,kb,
     * ki,nn-1,id,ucy,ea,det,v)
      u(1,1)=0
      u(2,1)=u(2,2)
      write(*,*) 'After call to nrk in uhx'
      write(*,*) 'ea:'
      write(*,*) ea
      iter = iter+1
      if (iter.lt.3) goto 150

      print*, "xuhx", x(nn), u(1,nn), u(2,nn)
 701  continue
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c RHS
c
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhskdbh(x,y,zk,ap,aq,f,fd,h,hd,ifd,n)
      implicit real*8 (a-h,o-z)
      real mass, rad, gconst, oms2

      include 'adipls.c.d.incl'
      dimension y(*),f(*),fd(ifd,*),h(*),hd(ifd,*),zk(*),ap(*),aq(10,*)

      common/csumma/ cs(50)

C RHS FOR CALCULATING U(X)=H(X).X (g&t 90 eq.3.25)

C DEFINING CONVERSION FACTOR

      mass = cs(2)/1000
      rad = cs(3)/100
      gconst = 6.67e-11

      oms2 = gconst*mass/(rad*rad*rad)

C DEFINING PI, Q, RHO, OMEGA

      pi = 4.d0*atan(1.d0)
      q = aq(1,n)
      rho = aq(2,n)
      dom = aq(3,n)
      d2om = aq(4,n)
      omgrtp = aq(5,n)

c DEFINING B(x)

      B1 = 8.*pi*x**6. * rho * omgrtp * dom / (3.*q**2.)
      B2 = x**3. * (4.*omgrtp*dom + 2.*x*(dom**2 + omgrtp*d2om)/3.)/q
      B = -(B1+B2)/oms2

C RHS

      fd(1,1) = 0.0
      fd(1,2) = 1.0
      f(1) = fd(1,1)*y(1) + fd(1,2)*y(2)
      fd(2,1) = 4./(x*x)
      fd(2,2) = (2./x) - 8.*pi*x*x*rho/q
      f(2) = fd(2,1)*y(1) + fd(2,2)*y(2) + B

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c BOUNDARY CONDITIONS
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bckdbh(x1,x2,y1,y2,zk,ap,aq,g,gd,ig,id,n)
      implicit real*8(a-h,o-z)
      real mass, rad, gconst, oms2

      include 'adipls.c.d.incl'
      dimension y1(*),y2(*),g(*),gd(ig,*),zk(*),ap(*),aq(10,*)
      common/nnpts/ nn
      common/csumma/ cs(50)

C BOUNDARY CONDITIONS FOR CALCULATING U(X)=H(X).X (g&t90 eq.3.25)

C DEFINING CONVERSION FACTOR

      mass = cs(2)/1000
      rad = cs(3)/100
      gconst = 6.67e-11

      oms2 = gconst*mass/(rad*rad*rad)

C DEFINING PI, Q, RHO, OMEGA

      pi = 4.d0*atan(1.d0)
      q = aq(1,nn)
      rho = aq(2,nn)
      dom = aq(3,nn)
      d2om = aq(4,nn)
      omgrtp = aq(5,nn)

c DEFINING S(x)

      S = x2**4. *omgrtp*(5.*omgrtp + 2.*x2*dom)/(3.*q)/oms2

C BC

      gd(1,1) = -1.0
      gd(1,2) = x1
      g(1) = gd(1,1)*y1(1) + gd(1,2)*y1(2)
      gd(2,1) = 1.0
      gd(2,2) = x2
      g(2) = gd(2,1)*y2(1) + gd(2,2)*y2(2) - S

      return
      end


