      subroutine nrk(x,y,zk,ap,aq,rhs,bc,ii,kk,ka,kb,ki,nn,id,ucy,ea,
     .det,v)
c
c  ***********************************************************************
c
c  package for solving boundary value problems for ordinary differential
c  equations with optional eigenvalues. accepts  one or two point boundary
c  conditions, and integral constraints.
c
c  method: replaces differential equations by centred, second order
c          accuracy difference equations, and solves these by
c          newton-raphson-kantorovich iteration, using forward
c          elimination and back-substitution. see baker, moore and
c          spiegel, q. j. mech. appl. math., vol. 24, p. 391
c          (1971).
c
c  package written by d. o. gough, institute of astronomy and
c  department of applied mathematics and theoretical phyisics,
c  university of cambridge, england.
c
c  **********************************************************************
c
c  newton-raphson-kantorovich setting up routine
c
c   space required by common/work/ is
c
c  (2*((ii+kk+1)*(6*ii+2*kk+1)+ii*(ii-ka)+ki*(ii+kk+2))
c      +  (ii+kk+1-ka)*ii*nn)*4  bytes
c
c   the programmer has the option of defining common/nrkchk/nrchk
c   and setting it to the number of single precision words that has been
c   reserved in labelled common/work/
c   if nrkwsp.gt.0,  nrk checks that it is great enough to accomodate th
c   problem.
c   if the work space is too small or if nrkwsp.le.0, appropriate
c   diagnostics are written
c
c   if an error is detected  v(1) is set to zero before return
c
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 16/4/90.
c
c  Note: this double precision version of the routine has same name
c  as single precision version, but is distinguished by its file name
c
      implicit double precision (a-h,o-z)
      integer v
      dimension x(nn),y(id,nn),zk(kk),ea(id,3),ap(1),aq(1),v(ii)
      common/nrkchk/nrkwsp
      common/work/ w(2000)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      external rhs,bc
      data ichk/0/
c
c
c   line printer  dsrn
c
c      set dimension paraeters
c
   20 ik=ii+kk
      ik1=ik+1
      ikka1=ik-ka+1
      ip=ii*ikka1
      ipn=ip*nn
      iiik=ii*ik
      iiik1=ii+ik1
      kid=ki
      if(ki.le.0)kid=1
c
c      location of variables
c  note: a,d,p must be stored consecutively and in order
c
   30 lq=1
      lf=lq+ip
      lh=lf+ii+ii
      lfd=lh+kid
      lhd=lfd+iiik+iiik
      lg=lhd+ki*ik1
      lgp=lg+ik1
      la=lgp+2*ik1*ik1
      ld=la+ii*iiik1
      lp=ld+iiik+ii
c
c   size of common work/
c
      iwork=(lp-lq)+ ipn
c
c   if nrkwsp has been set, check that enough space has been reserved
c   in  /work/
c
   40 if(nrkwsp.gt.0)go to 42
   41 if(ichk.ne.iwork.and.istdpr.gt.0) write(istdpr,1001)iwork
      ichk=iwork
      go to 50
   42 if(nrkwsp.ge.iwork)go to 50
   43 write(istdou,1002)nrkwsp,iwork
      if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,1002)
     *  nrkwsp,iwork
      v(1)=0
      return
c
c   call newton-raphson-kantorovich routine
c  note:logically gd and gp are equivalent
c   g,gst,eb and gs,p share store to save space, consequently ik1*8 byte
c   have been reserved for g in /work/ to accomodate gst
c
   50 call nrke(x,y,zk,ap,aq,rhs,bc,ii,kk,ka,kb,ki,nn,id,ucy,ea,det,v
     .,w(lf),w(lfd),w(lg),w(lg),w(lgp),w(lgp),w(lp),w(la),w(ld),w(lh)
     .,w(lhd),w(lq),w(lp),w(lg),ik,ik1,ikka1,iiik1,ip,ipn,iw,kid)
      return
c
c      diagnostic messages
c
 1001 format(1x,5('%'),3x,'size of common/work/ not checked by nrk',
     .i9,'*4bytes required ... execution continuing')
 1002 format(1x,5('%'),i9,'*4 bytes in common/work/ insufficient for
     . nrk',i9,'*4 bytes required ...  nrk execution terminated')
c
      end
      subroutine nrke(x,y,zk,ap,aq,rhs,bc,ii,kk,ka,kb,ki,nn,id,ucy,ea,
     .               det,v,f,fd,g,gst,gd,gp,gs,a,d,h,hd,q,p,eb,ik,
     .               ik1,ikka1,iiik1,ip,ipn,iw,kid)
c
c
c                  newton-raphson-kantorovich programme
c                  ************************************
c
c     if an error is detected, v(1) is set to zero before return
c
c
      implicit double precision (a-h,o-z)
      integer v
      dimension f(ii,2),fd(ii,ik,2),g(ik),gst(ik1),gd(ik1,ik1),
     .          gp(ik1,ik1,2),gs(ik,ik1,2),a(ii,iiik1),d(ii,ik1),
     .          h(kid),hd(kid,ik1),q(ii,ikka1),p(ipn)
c
      dimension x(nn),y(id,nn),zk(kk),ea(id,3),
     .          eb(ii),v(ii),ap(1),aq(1)
c
c
      external rhs,bc
c      note: rhs sets derivatives of h into d and not hd. therefore
c            it must assume that the first dimension of hd is the
c            same as that of fd
c
c
c
c
c     line printer dsrn
c
c     set counting limits
      kab=ka+kb
      kap=ik-ki
      ika=ii-ka
      ikka=ik-ka
      kp=kap-kab
      kip=ki+kp
c
    3 n1=nn-1
      i1=ii+1
      k1=kk+1
      kab1=kab+1
      kap1=kap+1
      ka1=ka+1
      ika1=ika+1
      kaik=ka*ik1
      ikka1k=ikka1*ik1
      ii2=ii*ii
      ikka1i=ikka1*ii
      ikkaik=ikka*ik1
c
c
c     compatibility test
      if(ka.le.ii.and.kp.ge.0.and.ka.ge.0.and.kb.ge.0.and.ki.ge.0.and.
     .   kk.ge.0) go to 10
      write(istdou,1000) ii,kk,ka,kb,ki
      if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,1000) 
     *  ii,kk,ka,kb,ki
      v(1)=0
      return
   10 continue
c
c
c     set v if necessary
   20 iv=0
      do 21 i=1,ii
      if(v(i).gt.0) go to 21
      iv=1
      v(i)=i
      if(istdpr.gt.0) write(istdpr,1100) i,i
   21 continue
      if(iv.eq.1.and.istdpr.gt.0) write(istdpr,1101) (v(i), i=1,ii)
c
c     set range of independent variable
   25 r=x(nn)-x(1)
      if(r.eq.0.) go to 407
c
c
c     **********
c     conditions solely at first boundary
c
c     empty derivative matrices
   30 do 34 j=1,ik
      do 31 i=1,ii
   31 fd(i,j,1)=0.0
      do 32 k=1,2
      do 32 i=1,kap
   32 gs(i,j,k)=0.0
      if(ki.eq.0) go to 34
      do 33 ig=1,ki
   33 d(ig,j)=0.0
   34 continue
c
c
c     boundary conditions and equations at first boundary
   35 ia=ik
      ib=ik1
      in=ii
      i=1
   36 call bc(x(1),x(nn),y,y(1,nn),zk,ap,aq,g,gs,ia,ib,nn)
      call rhs(x(1),y,zk,ap,aq,f,fd,h,d,in,i)
c     note that although the dimensions of d might not be adequate to
c     store hd, in that event there is always sufficient unused space
c     at the beginning of p to accomodate the overflow
c
c     set boundary derivative matrix
   40 do 43 i=1,ik
      iv=i
      if(i.le.ii) iv=v(i)
      do 41 j=1,kap
      do 41 k=1,2
   41 gp(j,i,k)=gs(j,iv,k)
c     (ensure that first and second b.c. matrices are loaded into gd)
      if(kab.eq.0) go to 43
      do 42 j=1,kab
   42 gd(j,i)=gd(j,i)+gp(j,i,2)
   43 gd(i,ik1)=g(i)
c     (ensure that eigenvalue dependence is loaded into gd)
      if(kp.eq.0.or.kk.eq.0) go to 45
      do 44 is=kab1,kap
      do 44 j=i1,ik
   44 gd(is,j)=gd(is,j)+gp(is,j,2)
c
c     first inversion (for equation 5)
   45 detsgn=1.
      if(ka.ne.0) go to 46
      det=0.
      go to 48
   46 call leqsd(gd,gd(1,ka1),ka,ikka1,ik1,ik1,err,kaik,ikka1k)
      if(err)  47,400,48
   47 detsgn=-detsgn
   48 det=log10( abs(err))
c     note that now gd(ial,nu)=-h(ial,nu,1) in equation 5
c
c     initial integration coefficients
      hn=(x(2)-x(1))/6.e0
      hn1=(x(3)-x(2))/6.e0
      if(ki.eq.0) go to 50
      wa=hn+hn1
      wb=hn-hn1
      wd=wa/hn
      ah=wd*(hn+wb)
      ah1=wa*wa*wd/hn1
      ah2=wa*(hn1-wb)/hn1
c     set integral constraint matrix
      do 49 ig=1,ki
      hd(ig,ik1)=h(ig)
      do 49 i=1,ik
      iv=i
      if(i.le.ii) iv=v(i)
   49 hd(ig,i)=d(ig,iv)
c
c
c     compute p and gamma (gd) at first mesh point
c     and store p in real*8 q for computation of e in loops 103 and 105
   50 if(ka.eq.0) go to 55
      do 54 nu=ka1,ik1
      j=(nu-ka-1)*ii+ika
      do 51 ial=1,ka
      wd=-gd(ial,nu)
      q(ika+ial,nu-ka)=wd
   51 p(ial+j)=wd
      if(kp.eq.0) go to 54
      do 53 is=kab1,kap
      wd=0.0
      do 52 ial=1,ka
   52 wd=wd-gd(is,ial)*gd(ial,nu)
   53 gd(is,nu)=wd+gd(is,nu)
   54 continue
c     contribution at first boundary to integral constraints
   55 if(ki.eq.0) go to 58
      do 57 nu=ka1,ik1
      do 57 ig=1,ki
      wd=0.0
      if(ka.eq.0) go to 57
      do 56 ial=1,ka
   56 wd=wd-hd(ig,ial)*gd(ial,nu)
   57 gd(ig+kap,nu)=ah*(wd+hd(ig,nu))
   58 continue
c
c
c
c     **********
c     beginning of preliminary outer loop
c
c     set storage indices
   70 in=1
      in1=2
c
   80 do 130 n=1,n1
      np1=n+1
      dx=3.e0*hn
c
c     check that independent variable is strictly monotonic
      if(r*dx.le.0.) go to 408
c
c     compute integration coefficients
      hn=hn1
      if(n+3.gt.nn) go to 81
      hn1=(x(n+3)-x(n+2))/6.e0
   81 if(ki.eq.0) go to 90
      if(n.eq.n1) go to 85
      if(in.eq.1) go to 84
c     check whether np1=nn-1 when nn is even
      if(np1.eq.n1) go to 84
      wa=hn+hn1
      wb=hn-hn1
      wd=wa/hn
      ah=wd*(hn+wb)+ah2
      ah2=wa/hn1
      ah1=wa*wd*ah2
      ah2=ah2*(hn1-wb)
      go to 90
c     integration coefficients at np1 when np1 is even
c     and nn-1 when nn is even
   84 ah=ah1
c     check whether np1=nn-2
      if(n.ne.nn-3) go to 90
c     integration coefficients at nn-2,nn-1 and nn when nn is even
      wa=hn+hn1
      wb=3.e0*hn+hn1
      ah=ah-hn1*hn1*hn1/(hn*wa)
      ah1=hn1*wb/hn+ah2
      ah2=hn1*(wb+hn1)/wa
      go to 90
c     integration coefficient at nn
   85 ah=ah2
c
c     equations at mesh point np1
   90 do 92 i=1,ik
      do 91 j=1,ii
   91 fd(j,i,in1)=0.e0
      do 92 ig=1,ki
   92 d(ig,i)=0.e0
      ia=ii
      i=np1
      call rhs(x(np1),y(1,np1),zk,ap,aq,f(1,in1),fd(1,1,in1),h,d,ia,i)
      do 93 i=1,ii
   93 gst(i)=y(v(i),n)-y(v(i),np1)
   96 if(ki.eq.0) go to 100
      do 97 i=1,ik
      iv=i
      if(i.le.ii) iv=v(i)
      do 97 ig=1,ki
   97 hd(ig,i)=d(ig,iv)
c
c     beginning of loop to construct e matrix
  100 do 109 i=1,ii
c
c     construct a and b matrices (b is loaded into d)
      do 101 j=1,ii
      a(i,j)=dx*fd(v(i),v(j),in)
  101 d(i,j)=dx*fd(v(i),v(j),in1)
      a(i,i)=a(i,i)+1.0e0
      d(i,i)=d(i,i)-1.0e0
      do 102 k= i1,ik
  102 d(i,k)=dx*(fd(v(i),k,in)+fd(v(i),k,in1))
      d(i,ik1)=dx*(f(v(i),in)+f(v(i),in1))+gst(i)
c     construct d matrix (b block is already loaded)
      if(ka.eq.0) go to 107
      do 104 m=i1,ik1
      wd=0.0
      do 103 ial=1,ka
  103 wd=wd+a(i,ial)*q(ika+ial,m-ka)
  104 d(i,m)=d(i,m)+wd
c
c     construct beginning of e matrix and set into end of a
      if(ika.eq.0) go to 107
      do 106 j=1,ika
      wd=0.0
      do 105 ial=1,ka
  105 wd=wd+a(i,ial)*q(ika+ial,j)
  106 a(i,j+ka)=wd+a(i,j+ka)
  107 continue
c
c      douglas is playing around here and we need to dort it out ]]
c
       do 108 m=1,ik1
  108 a(i,ii+m)=d(i,m)
c
c     in 110 and 111,121 loops a(i,ii+m) will be used in place of d(i,m)
  109 continue
c     end of loop to construct e matrix
c
c
c     compute p matrix
c     and store current value in real*8 q for computation of e
c     in loops 103 and 105
c
  110 call leqsd(a(1,ka1),d(1,ka1),ii,ikka1,ii,ii,err,ii2,ikka1i)
      if(err) 111,403,112
  111 detsgn=-detsgn
  112 k=n*ip
      det=det+log10( abs(err))
      do 115 nu=ka1,ik1
      j=(nu-ka-1)*ii+k
      do 115 i=1,ii
      wd=-d(i,nu)
      q(i,nu-ka)=wd
  115 p(i+j)=wd
c
c     compute gamma matrix
  120 if(kip.eq.0.or.ika.eq.0) go to 125
      do 124 is=kab1,ik
      do 123 nu=ka1,ik1
      wd=0.0
      do 121 ib=ka1,ii
  121 wd=wd-gd(is,ib)*d(ib-ka,nu)
      if(nu.gt.ii) go to 122
      gst(nu-ka)=wd
      go to 123
  122 gst(nu-ka)=wd+gd(is,nu)
  123 continue
      do 124 nu=ka1,ik1
  124 gd(is,nu)=gst(nu-ka)
c     integral constraints at point np1
  125 if(ki.eq.0) go to 129
      do 127 ig=1,ki
      hd(ig,ik1)=h(ig)
      do 127 nu=ka1,ik1
      wd=0.0
      if(ka.eq.0) go to 127
      do 126 ial=1,ka
  126 wd=wd-hd(ig,ial)*d(ika+ial,nu)
  127 gd(ig+kap,nu)=gd(ig+kap,nu)+ah*(wd+hd(ig,nu))
c
c     reset storage indices
  129 i=in
      in=in1
      in1=i
  130 continue
c
c     end of preliminary outer loop
c     **********
c
c
c     **********
c     remaining boundary conditions
c
  200 if(ka.eq.ik) go to 233
c     compute coefficients for equation 12a and load into end of gd
      if(kp.eq.0) go to 210
      do 205 is=kab1,kap
      do 202 nu=ka1,ik1
      wd=0.0
      do 201 ial=1,ka1
  201 wd=wd-gp(is,ial,2)*d(ika+ial,nu)
  202 gd(is,nu)=gd(is,nu)+wd
      do 203 ia=ka1,ii
  203 gd(is,ia)=gd(is,ia)+gp(is,ia,2)
  205 continue
  210 continue
c
c
c     compute coefficients for equation 13 and load into middle of gd
  220 if(ka.eq.0.or.kb.eq.0) go to 230
      do 223 m=ka1,kab
      do 222 nu=ka1,ik1
      wd=0.0
      do 221 ial=1,ka
  221 wd=wd-gd(m,ial)*d(ika+ial,nu)
  222 gd(m,nu)=wd+gd(m,nu)
  223 continue
c
c     solve equations 12a and 13 (answer has wrong sign)
  230 if(ikka.eq.0) go to 233
  231 call leqsd(gd(ka1,ka1),gd(ka1,ik1),ikka,1,ik1,ikka,err,ikkaik,
     .  ikka)
      if(err) 233,404,234
  233 detsgn=-detsgn
  234 det=det+log10( abs(err))
      gd(ik1,ik1)=-1.0e0
      gd(ik1,1)=-1.0e0
c     (note that the first element of the second half of gp may now
c     have been overwritten)
c
c
c
c     **********
c     iterated solution
c
c     empty ea,eb
  300 do 301 i=1,ii
      eb(i)=0.0
      do 301 j=1,3
  301 ea(i,j)=0.0
c
c     solution at second boundary
  310 if(kk.eq.0) go to 313
      do 312 k=1,kk
      wd=gd(ii+k,ik1)
      gd(ii+k,1)=wd
      zk(k)=zk(k)-ucy*wd
  312 continue
  313 if(ka.ge.ii) go to 320
      do 318 ia=ka1,ii
      wa= abs(gd(ia,ik1))
      wd=y(v(ia),nn)
  315 wd=wd-ucy*gd(ia,ik1)
      wb= abs(wd)
      ea(v(ia),1)=wa
      eb(ia)=wb
      if(wb.lt.1.0e-10) go to 316
      ea(v(ia),2)=wa/wb
      ea(v(ia),3)= float(nn)
  316 continue
      y(v(ia),nn)=wd
  318 continue
c
c     remainder of solution   (beginning of second outer loop)
c
c     set storage indices
  320 in=1
      in1=ik1
c
      do 330 m=1,nn
      n=nn+1-m
      np1=ip*(n-1)
      do 328 i=1,ii
      if(i.gt.ika) go to 321
      k=n-1
      if(k.lt.1) go to 328
      j=i+ka
      go to 322
  321 j=i-ika
      k=n
  322 wd=0.0
      ia=i+np1
      do 323 nu=ka1,ik1
  323 wd=wd+p(ia+(nu-ka-1)*ii)*gd(nu,in1)
      wa= abs(wd)
      wb=y(v(j),k)-ucy*wd
      y(v(j),k)=wb
  325 wb= abs(wb)
      ea(v(j),1)=ea(v(j),1)+wa
      eb(j)=eb(j)+wb
      if(wb.lt.1.0e-10) go to 326
      wb=wa/wb
  326 if(wb.lt.ea(v(j),2)) go to 327
      ea(v(j),2)=wb
      ea(v(j),3)= float(k)
  327 continue
      if(i.le.ika) gd(j,in)=wd
  328 continue
c
c     reset storage indices
      i=in
      in=in1
      in1=i
  330 continue
c
c     end of second outer loop
c
c     iteration complete
c     **********
c
c
c     set mean relative corrections
  340 do 341 i=1,ii
  341 ea(v(i),1)=ea(v(i),1)/eb(i)
c     set det
      if(detsgn) 345,350,350
  345 det=det*1.e10
c
c
  350 return
c
c
c     diagnostics
c     **********
  400 write(istdou,1001)
      if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdou,1001)
  401 do 402 k=1,2
      do 402 j=1,ik
      write(istdou,1002) (gs(i,j,k), i=1,kap)
      if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,1002)
     *  (gs(i,j,k), i=1,kap)
  402 continue
      write(istdou,1101) (v(i), i=1,ii)
      if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,1101)
     *  (v(i), i=1,ii)
      go to 500
  403 write(istdou,1003)n
      if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,1003)n
      go to 500
  404 write(istdou,1004)
      if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,1004)
      do 405 i=1,kap
      do 405 j=1,ik
      do 405 k=1,2
  405 gs(i,j,k)=0.0
      call bc(x(1),x(nn),y,y(1,nn),zk,ap,aq,g,gs,ik,ik1,nn)
      go to 401
  407 write(istdou,1006) nn,x(1)
      if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdou,1006) nn,x(1)
      go to 500
  408 dx=2.0*dx
      write(istdou,1007) n,x(n),np1,x(np1),dx,r,nn,(x(i), i=1,nn)
      if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,1007) 
     *  n,x(n),np1,x(np1),dx,r,nn,(x(i), i=1,nn)
c
  500 v(1)=0
      return
c
c
 1000 format(//1x,10('*'),10x,'improper formulation detected by nrk',
     . 10x,10('*')/17x,'ii =',i3,4x,'kk =',i3,4x,'ka =',i3,4x,'kb =',i3,
     . 4x,'ki =',i3/)
 1001 format(//1x,10('*'),5x,'first boundary condition matrix singular i
     .n nrk',5x,10('*')//' gd ='/)
 1002 format(1x,1p14e9.1/)
 1003 format(//1x,10('*'),5x,'e matrix singular at n =',i4,' in nrk',
     . 5x,10('*'))
 1004 format(//1x,10('*'),5x,'second boundary matrix singular in nrk',
     .  5x,10('*')//' gd ='/)
 1006 format(//1x,10('*'),5x,'null range of independent variable in ',
     .       'nrk',5x,10('*')//21x,'x(1) = x(',i4,') =',1pe14.6/)
 1007 format(//1x,10('*'),5x,'independent variable not monotonic in',
     .       ' nrk',5x,10('*')//16x, 'x(',i4,') =',1pe12.4,
     .       ',    x(',i4,') =',1pe12.4/16x,'difference =',1pe12.4,
     .       ',    range =',1pe12.4//1x,'x(n), n=1,',i4,/
     .       (1x,1p10e13.5))
c
c     warning message
 1100 format(1x,5('*'),3x,'v set in nrk        v(',i2,') =',i3)
 1101 format(/9x,'v =',7(1x,5i3))
c
c
      end
c..      subroutine leqsd(a,b,n,m,ia,ib,det,isa,isb)
c..      implicit double precision (a-h,o-z)
c..      dimension a(ia,1),b(ib,1)
c..c
c..      call leq(a,b,n,m,ia,ib,det)
c..      return
c..      end
