      subroutine squint(x,f,g,y,nn,if,ig,iy)
c     sets integral 
c     int(f(x)(g(x))**(-1/2)dx, from x(1) to x(n),
c     where f(x(n)) = f(1,n) and g(x(n)) = g(1,n), into
c     real y(1,n),  n=1,nn   
c
c     Note: the integrand is set to zero where undefined (i.e. where
c     g .lt. 0)
c
c     the independent variable x, which need not be uniformly divided
c     or increasing with n (but must be monotonic), must be supplied 
c     by the calling programme.
c
      implicit double precision (a-h,o-z)
      dimension x(1),f(if,1),g(ig,1),y(iy,1)
c
      y(1,1)=0
c
c  start loop in x
c
      do 50 n=2,nn
      n1=n-1
      g1=g(1,n1)
      g2=g(1,n)
      delg=g2-g1
c
c  test for cases on g
c
      if(g1.le.0.and.g2.le.0) then
        aint0=0
        aint1=0
      else if(g1.le.0) then
        sqr2=sqrt(g2)
        aint0=2*sqr2/delg
        aint1=-4*(2*g1-delg)*sqr2/(3*delg*delg)
      else if(g2.le.0) then
        sqr1=sqrt(g1)
        aint0=-2*sqr1/delg
        aint1=4*sqr1*g1/(3*delg*delg)
      else
c
c  general case. test for using expansion
c
        eps=delg/g1
        sqri1=1./sqrt(g1)
        if(abs(eps).lt.0.01) then
          aint0=sqri1*(1+eps*(0.125*eps-0.25))
          aint1=0.5*sqri1*(1+eps*(0.1875*eps-0.333333333))
        else
          sqreps=sqrt(1+eps)
          aint0=2*sqri1*(sqreps-1)/eps
          aint1=4*sqri1*(1-(1-eps/2)*sqreps)/(3*eps*eps)
        end if
      end if
c
c  now g integrals are set. set contribution to final integral
c
   50 y(1,n)=y(1,n1)+
     *       (x(n)-x(n1))*(f(1,n1)*aint0+(f(1,n)-f(1,n1))*aint1)
c
      return
      end
