      subroutine vinta(x,a,y,nn,ia,id)   
c     vinta  
c     sets integral of the real variable a(1,n)  ( =a(1,x(n)) into   
c     real y(1,n),  n=1,nn   
c
c
c     the independent variable x, which need not be uniformly divided
c     or increasing with n (but must be monotonic), must be supplied 
c     by the calling programme   
c
c
c
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90.
c
c  Note: this double precision version of the routine has same name
c  as single precision version, but is distinguished by its file name
c
      implicit double precision(a-h,o-z)
      dimension a(nn),y(nn),x(nn)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      if(nn.eq.1) then
	if(istdpr.gt.0) write(istdpr,100)
	y(1)=0
        return
      else if(nn.eq.2) then
c
c  use trapezoidal integration
c
	y(1)=0
	y(1+id)=0.5d0*(x(2)-x(1))*(a(1)+a(1+ia))
	return
      end if
c
    1 wc=a(1)
      ir=1   
    3 is=ir*ia   
      ird=ir*id  
      n=nn   
      m=n/2  
      ie=n-m*2   
      if(ie.eq.0) m=m-1  
      w=1.0d0
      r=x(n)-x(1)
      if(r.eq.0.d0) go to 50   
c
      j3=1   
      k3=1   
      y(1)=0.d0
      ah2=0.d0 
      vsimf=0.d0   
c
    4 do 7 i=1,m 
      i2=i*2 
      i1=i2-1
      i3=i2+1
      hn=(x(i2)-x(i1))/6.d0  
      hn1=(x(i3)-x(i2))/6.d0 
      if(r*hn.le.0.d0.or.r*hn1.le.0.d0) go to 51 
      wa=hn+hn1  
      wb=hn-hn1  
      wd=wa/hn   
      ah=wd*(hn+wb)  
      ah2=wa/hn1 
      ah1=wa*wd*ah2  
      ah2=ah2*(hn1-wb)   
      bh1=hn+3.d0*hn1
      bh2=hn/hn1 
      bh=(bh1+hn)*hn/wa  
      bh1=bh1*bh2
      bh2=-hn*hn*bh2/wa  
      wa=wc  
      k2=k3+is   
      k3=k2+is   
      wb=a(k2)   
      wc=a(k3)   
    6 j2=j3+ird  
      j3=j2+ird  
      y(j2)=vsimf+bh*wa+bh1*wb+bh2*wc
      vsimf=vsimf+ah*wa+ah1*wb+ah2*wc
      y(j3)=vsimf
    7 continue   
c
      if(ie.ne.0) return 
c
      wa=x(n-1)  
      hn=(wa-x(n-2))/6.d0
      hn1=(x(n)-wa)/6.d0 
      if(r*hn.le.0.d0.or.r*hn1.le.0.d0) go to 52 
      wa=hn+hn1  
      wd=hn1/hn  
      ah=-hn1*hn1*wd/wa  
      vsimf=vsimf+ah*wb 
      wb=3.0d0*hn+hn1 
      ah1=wd*wb 
      ah2=hn1*(wb+hn1)/wa   
      wa=wc 
      wb=a(k3+is)   
    9 vsimf=vsimf+ah1*wa+ah2*wb 
      y(j3+ird)=vsimf   
c   
      return
c   
c   
c     diagnostics   
   50 if(istdpr.gt.0) write(istdpr,60) n,x(1)
      return
   51 i=i2+1
      if(istdpr.gt.0) 
     *  write(istdpr,61) i1,x(i1),i2,x(i2),i,x(i),n,x(1),x(n)  
      return
  52  if(istdpr.gt.0) 
     *  write(istdpr,61) i2,x(i2),i1,x(i1),n,x(n),n,x(1),x(n)  
      return
c   
c   
c   
   60 format(//1x,10('*'),5x,'null range of independent variable in ',  
     *       'intf/inta',5x,10('*')//21x,'x(1) = x(',i4,') =',1pe14.6/) 
   61 format(//1x,10('*'),5x,'independent variable not monotonic in',   
     *       ' intf/inta',5x,10('*')//16x,'x(',i4,') =',1pe13.5,
     *       ',    x(',i4,') =',1pe13.5,',    x(',i4,') =',1pe13.5/ 
     *       16x,'number of mesh points =',i4, 3x,'x(1) =',1pe13.5, 
     *       8x,'x(n) =',1pe13.5/)  
c   
  100 format(/' **** Error in vinta: nn = 1')
      end   
