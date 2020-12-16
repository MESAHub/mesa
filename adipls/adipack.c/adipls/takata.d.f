      subroutine takata(x,y,nn,iasn,data,aa,iy,iaa,iord,sig)
!
!  this subroutine intends(!) to define the Takata's variables given
!  in Takata,M. 2006, ESASP, 624, 26 (section3) and to identify the 
!  modes according to the scheme therein.
!
!  Dated Feb2008
!
!  Modified 5/5/10 setting nnmax with include file and using consistent
!  value of G
!
!  Modified 29/7/10, accounting for using every second mesh point with
!  Richardson extrapolation.
!  Also fixing dividing by zero gravity at core, in setting T2.
!
!  Modified 15/9/10, fixing dividing by zero gravity at centre, in
!  setting yy(2,.)
!
!  Modified 11/10/10, correcting the use of aa(i,.), replacing 
!  aa(i,.) by aa(i-1,.) and aa(1,.) by x(.).
!
      IMPLICIT NONE
      integer iaa,iord,iordg,iordp,iy,k,n,nn,iasn,nnmax,ns,
     *        istdin,istdou,istdpr,istder,icount
c..      parameter (nnmax=10000)
      include 'adipls.c.d.incl'
      real*8 aa(iaa,nn),data(*),dfyy,G,J(nnmax),phx1,phx2,
     *       phy1,phy2,phy2m,pi,rho(nnmax),sig,T2(nnmax),x(nn),
     *       y(iy,nn),yy(2,nnmax),gr(nnmax),dphipdr(nnmax), cgrav
      common/ccgrav/ cgrav
      common /cstdio/ istdin, istdou, istdpr, istder
      data icount /0/
      save

      icount=icount+1
!
      pi=3.141593d0
      G=cgrav
!
      do n=1,nn
        if(n.eq.1.or.mod(n-2,iasn).eq.0) then
          rho(n)=data(1)*aa(1,n)*aa(5,n)/(4*pi*(data(2))**3)
          J(n)=1-aa(5,n)/3
          gr(n)=G*aa(1,n)*x(n)*data(1)/data(2)**2
!       gr: gravitational acceleration
          if(x(n).gt.0) then
            dphipdr(n)=-(gr(n)/x(n))*(y(4,n)-y(3,n))-(y(3,n)*G*data(2)*
     *        4*pi*rho(n))
!       dphipdr is d(phi')/dr 
            T2(n)=-y(3,n)/(3*x(n))-dphipdr(n)/(3*gr(n))
!       T2 is the 2nd term in both of the new variables
            yy(1,n)=(J(n)*y(1,n)/x(n))+T2(n)
            yy(2,n)=J(n)*y(2,n)*G*data(1)*sig/(2*gr(n)*data(2)**2)+
     *        J(n)*y(3,n)/x(n)+T2(n)   
          end if
        end if
      end do
      if(x(1).eq.0) then
        dphipdr(1)=dphipdr(2)
        T2(1)=T2(2)
        yy(1,1)=yy(1,2)
        yy(2,1)=yy(2,2)
      end if
!
      iordp=0
      iordg=0
      do 100 k=1,nn-1
        if(k.eq.1.or.mod(k-2,iasn).eq.0) then
          phx1=yy(1,k)
          phy1=yy(2,k)
          phx2=yy(1,k+iasn)
          phy2=yy(2,k+iasn)
!
!     dfyy was originally evaluated using these two lines:
!          phy2m=(phx1*phy2-phx2*phy1)/(phx1-phx2)
!          dfyy=phy2m*(phx2-phx1)
!     where the first introduces a potential FPE.  But they appear to
!     simplify to this expression (WHB, 2020-01-21)
          dfyy = phx2*phy1-phx1*phy2
           if (phx1*phx2.LE.0.and.dfyy.LT.0) then
             iordp=iordp+1
           else if (phx1*phx2.LE.0.and.dfyy.gt.0) then
             iordg=iordg+1
           end if
         end if
  100 CONTINUE
      if (iordp.GE.iordg) then 
        iord=iordp-iordg+1
      else 
        iord=iordp-iordg
      end if
c..      if(icount.le.4) write(99,'(i2,1p10e13.5)')
c..     *  (icount,x(n),y(1,n),y(2,n),rho(n),J(n),gr(n),dphipdr(n),T2(n),
c..     *  yy(1,n),yy(2,n),n=1,nn)
!
      end
