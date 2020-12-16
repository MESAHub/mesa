      subroutine nrkm(x,y,zk,ap,aq,rhs,bc,ii,kk,nr,ki,nn,id,ucy,ea,
     .det,v)
c
c  boundary value problem with matching
c  ************************************
c
c  the following note applies only on an ibm computer
c
c  (note that only double precision version has been implemented. thus x
c   real*8 and y and zk are complex*16. for historical reasons ucy,ea and
c   det are real*4)
c
c  newton-raphson-kantorovich setting up routine
c
c   space required by common/work/ is certainly not
c
c  (2*((ii+kk+1)*(6*ii+2*kk+1)+ii*(ii-ka)+ki*(ii+kk+2))
c      +  (ii+kk+1-ka)*ii*nn)*8 bytes
c  this formula should be fixed up. however, with present logic, where
c  ii, kk, etc are set as the regions are entered, the exact value of
c  the space required can clearly not be known beforehand; but its
c  maximum value could, given the values of ii and kk passed to the
c  routine.
c
c   the programmer has the option of defining common/nrmchk/nrkwsp
c   and setting it to the number of 8 byte words that has been
c   reserved in labelled common/work/
c   if nrkwsp.gt.0, nrkm checks that it is great enough to accomodate th
c   problem.
c   if the work space is too small or if nrkwsp.le.0, appropriate
c   diagnostics are written
c   if ncfac (also in common/nrmchk/) is .ne.1 centered differences are
c   used to represent derivatives.
c   if ncfac.eq.1 centralization factors phi(i), i=1,ii, must have been
c   set by the user in common/cntfct/phi
c   if irsdif.eq.1 difference equations may be reset by call of
c   user-supplied subroutine resdif
c
c   if an error is detected  v(1) is set to zero before return
c
c  nr is number of regions
c  ii and kk are maximum order and maximum number of eigenvalues
c
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
c  Modified 26/4/99, adding initializations of kab1 and kap
c  (Rather trivial, as it happens).
c
      implicit double precision (a-h,o-z)
      integer v
      dimension x(nn),y(id,nn),zk(kk),ea(id,nr,3),ap(1),aq(1),v(ii)
c
c  to avoid problems with equivalencing real and integer arrays,
c  introduce an array to store counting indices for regions
c  (note: number of regions must be less than 142)
c
      dimension iwlstor(1000)
      common/nrmchk/nrkwsp,ncfac,irsdif
      common/work/ w(1000)
      common/stps81/ lq,lf,lh,lfd,lhd,lg,lgp,lgs,listor,la,ld,lp,lpmax,
     *  nrkwsa,nvrtot
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      external rhs,bc
      data ichk/0/
c
c   line printer  dsrn
c
   10 iw=istdou
c
c      set dimension parameters
c
   20 ik=ii+kk
      ik1=ik+1
      ip=ii*ik1
      ik2=ik+ik
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
      lgp=lg+ik2
      lgs=lgp+ik2*ik2+ik2
      listor=lgs+ik2*ik2+ik2
      la=listor+7*nr
      ld=la+ii*iiik1
      lp=ld+iiik+ii
c
c   call newton-raphson-kantorovich routine
c  note:logically gd and gp are equivalent
c   g,gst,eb and gs,p share store to save space, consequently ik1*16
c   bytes have been reserved for g in /work/ to accomodate gst
c
   50 continue
c
      call nrkme(x,y,zk,ap,aq,rhs,bc,ii,kk,nr   ,ki,nn,id,ucy,ea,det,
     .  v,w(lf),w(lfd),w(lg),w(lg),w(lgp),w(lgp),w(lgs),w(la),
     .  w(ld),w(lh),w(lhd),w(lq),w(lp),w(lg),iwlstor,ik,ik2,ik1,
     .  iiik1,lp,iw,kid)
      return
c
c      diagnostic messages
c
c
      end
      subroutine nrkme(x,y,zk,ap,aq,rhs,bc,iimax,kkmax,nr,kimx,nnmx,
     .  id,ucy,ea,det,v,f,fd,g,gst,gd,gp,gs,a,d,h,hd,q,p,eb,iorst,
     .  ikmax,ikmax2,ik1max,iiik1m,lp,iw,kid)
c
c
c                  newton-raphson-kantorovich programme
c                  ************************************
c
c     if an error is detected, v(1) is set to zero before return
c
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
c
      implicit double precision (a-h,o-z)
      integer v
      logical cfct
      dimension f(iimax,2),fd(iimax,ikmax,2),g(ikmax2),gst(ikmax2),
     .  gd(ikmax2,ikmax2),gp(ikmax2,ikmax2),gs(ikmax2,ikmax2),
     .   a(iimax,iiik1m),d(iimax,ik1max),h(kid),hd(kid,ik1max),
     .  q(iimax,ik1max),p(1),iorst(7,nr)
c
      dimension x(nnmx),y(id,nnmx),zk(kkmax),ea(id,nr,3),
     .  eb(iimax),v(iimax),ap(1),aq(1)
c
c
      common/nrmchk/ nrkwsp,ncfac,irsdif
      common/cntfct/ phi(10)
      common/stps81/ ist81(12),lpmax,nrkwsa,nvrtot
      external rhs,bc
c
c  ilngwk is 1 or 2, depending on whether first part of work is
c  single or double precision
c
      data ilngwk /1/
c
c      note: rhs sets derivatives of h into d and not hd. therefore
c            it must assume that the first dimension of hd is the
c            same as that of fd
c
c
c
c
c
c  start step in regions
c  $$$$$$$$$$$$$$$$$$$$$
c
      ipn=0
      iqs=0
      nqs=1
      nqf=1
      kqs=0
      kqs1=1
      kaprv=0
      iiprv=0
      kkprv=0
      ikprv=0
      do 150 nnqq=1,nr
c  set left hand conditions
      kqspr1=kqs1
      kqs1=kqs+1
      nu=nnqq
      ia=ikmax2
c  empty derivative matrix
      do 2 i=1,ia
      do 2 j=1,ia
    2 gs(i,j)=0
c
      i=nqs
c
      call bc(x(nqf),x(nqs),y(1,nqf),y(1,nqs),zk(kqspr1),zk(kqs1),
     .  ap,aq,g,gs,ia,i,km,ii,kk,nn,nu)
c  force no integral constraints or two point boundary conditions
      ki=0
      kp=0
c  force kip to zero (meaning of kip is unclear, but it is used in
c  statement no 120)
c
c  this should be checked
c  **********************
c
      kip=0
c
      kmapr=km+kaprv
      ka=kmapr-ikprv
c
c     set counting limits
      kab=ka+kb
      ik=ii+kk
      ik1=ik+1
      ikka1=ik1-ka
      ika=ii-ka
      ikka=ik-ka
      kap=ik-ki
c
    3 n1=nn-1
      i1=ii+1
      k1=kk+1
      kab1=kab+1
      km1=km+1
      kaprv1=kaprv+1
      ikprv1=ikprv+1
      kmapr1=km1+kaprv
      ikkapr=ikprv-kaprv
      ikikp1=ik1+ikprv
      ka1=ka+1
      ika1=ika+1
      kaik=ka*ik1
      ikka1k=ikka1*ik1
      ii2=ii*ii
      ikka1i=ikka1*ii
      ikkaik=ikka*ik1
c
c     compatibility test
      if(ka.ge.0.and.ka.le.ii) go to 10
      write(iw,1000) ii,kk,ka,km,iiprv,kkprv,kaprv
      v(1)=0
      return
   10 continue
c
c
c     set v if necessary
   20 iv=0
      do 21 i=1,ii
      j=i+iqs
      if(v(j).gt.0 .and.v(j).le.ii) go to 21
      iv=1
      v(j)=i
      write(iw,1100) j,i
   21 continue
      if(iv.eq.1) write(iw,1101) (v(i+iqs), i=1,ii)
c
c  set centralization factors to 0.5 if ncfac.ne.1
      if(ncfac.eq.1) go to 25
      do 22 i=1,ii
   22 phi(i+iqs)=0.5
c
c     set range of independent variable
   25 nqf=nqs+nn-1
      r=x(nqf)-x(nqs)
      if(r.eq.0.) go to 407
c  centralization factors
      cfct=ncfac.ne.1
      if(.not.cfct) go to 28
      pii=0.5
      qii=0.5
   28 continue
c
c
c     **********
c     conditions solely at first boundary
c
c     empty derivative matrices
   30 do 34 j=1,ik
      do 31 i=1,ii
   31 fd(i,j,1)=0.
      if(ki.eq.0) go to 34
      do 33 ig=1,ki
   33 d(ig,j)=0.
   34 continue
c
c
c    equations at first boundary
   35 in=iimax
      i=nqs
      nu=nnqq
      call rhs(x(nqs),y(1,nqs),zk(kqs1),ap,aq,f,fd,h,d,in,i,nu)
c     note that although the dimensions of d might not be adequate to
c     store hd, in that event there is always sufficient unused space
c     at the beginning of p to accomodate the overflow
c
c     set boundary derivative matrix
c
   40 do 4002 iu=1,km
 4002 gp(iu,ikikp1)=g(iu)
      if(iiprv.eq.0) go to 42
      j=iqsprv
      do 41 i=1,ikprv
      iv=i
      if(i.le.iiprv) iv=v(j+i)
      do 41 iu=1,km
   41 gp(iu,i)=gs(iu,iv)
c
      if(kaprv.le.0) go to 42
      do 415 iu=1,km
      do 415 it=kaprv1,ikprv1
      wd=0.
      do 413 ial=1,kaprv
  413 wd=wd+gp(iu,ial)*q(ikaprv+ial,it-kaprv)
      j=it
      if(it.eq.ikprv1) j=ikikp1
  415 gp(iu,j)=gp(iu,j)+wd
c
   42 do 43 i=1,ik
      iv=i
      if(i.le.ii) iv=v(i+iqs)
      is=i+ikprv
      ivs=iv+ikprv
      do 43 iu=1,km
   43 gp(iu,is)=gs(iu,ivs)
c
c     first inversion (for equation 5)
   45 detsgn=1.
      if(km.ne.0) go to 46
      det=0.
      go to 48
   46 call leqsd(gd(1,kaprv1),gd(1,kmapr1),km,ikka1,ikmax2,
     .  ikmax2,cerr,kaik,ikka1k)
      err=abs(cerr)
      if(err)  47,400,47
   47 det=log10(err)
c     note that now gd(ial,nu)=-h(ial,nu,1) in equation 5
c
c     initial integration coefficients
   48 n2=nqs+1
      n3=nqs+2
      hn=(x(n2)-x(nqs))/6.
      hn1=(x(n3)-x(n2))/6.
      if(ki.eq.0) go to 50
      wa=hn+hn1
      wb=hn-hn1
      we=wa/hn
      ah=we*(hn+wb)
      ah1=wa*wa*we/hn1
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
c     and store p in complex q for computation of e in loops 103 and
c     105
   50 if(ikkapr.le.0) go to 508
      do 502 it=kmapr1,ikikp1
      do 501 i=1,ikkapr
  501 p(i+ipn)=-gd(i,it)
  502 ipn=ipn+ikkapr
c  test size of ipn
      if(nrkwsp.gt.0.and.ipn+ilngwk*lp.gt.nrkwsp) go to 410
c
  508 if(ka.gt.0) go to 509
      ipn=ipn+ii*ikka1
      go to 55
  509 ipn=ipn+ika
      do 54 it=kmapr1,ikikp1
      do 51 ial=1,ka
      wd=-gd(ial+ ikkapr,it)
      q(ika+ial,it-kmapr)=wd
   51 p(ial+ipn)=wd
      ipn=ipn+ii
      if(kp.eq.0) go to 54
      do 53 is=kab1,kap
      wd=0.
      do 52 ial=1,ka
   52 wd=wd-gd(is,ial)*gd(ial,nu)
   53 gd(is,nu)=wd+gd(is,nu)
   54 continue
      ipn=ipn-ika
c     contribution at first boundary to integral constraints
   55 if(ki.eq.0) go to 58
      do 57 nu=ka1,ik1
      do 57 ig=1,ki
      wd=0.
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
   80 nqf1=nqf-1
      do 130 n=nqs,nqf1
      np1=n+1
      dx=6.*hn
c
c     check that independent variable is strictly monotonic
      if(r*dx.le.0.) go to 408
c
c     compute integration coefficients
      hn=hn1
      if(n+3.gt.nqf) go to 81
      hn1=(x(n+3)-x(n+2))/6.
   81 if(ki.eq.0) go to 90
      if(n.eq.n1) go to 85
      if(in.eq.1) go to 84
c     check whether np1=nn-1 when nn is even
      if(np1.eq.n1) go to 84
      wa=hn+hn1
      wb=hn-hn1
      we=wa/hn
      ah=we*(hn+wb)+ah2
      ah2=wa/hn1
      ah1=wa*we*ah2
      ah2=ah2*(hn1-wb)
      go to 90
c     integration coefficients at np1 when np1 is even
c     and nn-1 when nn is even
   84 ah=ah1
c     check whether np1=nn-2
      if(n.ne.nn-3) go to 90
c     integration coefficients at nn-2,nn-1 and nn when nn is even
      wa=hn+hn1
      wb=3.*hn+hn1
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
   91 fd(j,i,in1)=0.
      do 92 ig=1,ki
   92 d(ig,i)=0.
      ia=iimax
      i=np1
      nu=nnqq
      call rhs(x(np1),y(1,np1),zk(kqs1),ap,aq,f(1,in1),fd(1,1,in1),h,
     .  d,ia,i,nu)
      do 93 i=1,ii
      iv=v(i+iqs)
   93 gst(i)=y(iv,n)-y(iv,np1)
   96 if(ki.eq.0) go to  98
      do 97 i=1,ik
      iv=i
      if(i.le.ii) iv=v(i+iqs)
      do 97 ig=1,ki
   97 hd(ig,i)=d(ig,iv)
c
c     beginning of loop to construct e matrix
   98 do 10201 i=1,ii
      iv=v(i+iqs)
c  centralization factors
      if(cfct) go to 100
      pii=phi(iv+iqs)
      qii=1.-pii
c
c     construct a and b matrices (b is loaded into d)
  100 do 101 j=1,ii
      jv=v(j+iqs)
      a(i,j)=dx*qii*fd(iv,jv,in)
  101 d(i,j)=dx*pii*fd(iv,jv,in1)
      a(i,i)=a(i,i)+1
      d(i,i)=d(i,i)-1
      if(i1.gt.ik) go to 10201
      do 102 k= i1,ik
  102 d(i,k)=dx*(qii*fd(iv,k,in)+pii*fd(iv,k,in1))
10201 d(i,ik1)=dx*(qii*f(iv,in)+pii*f(iv,in1))+gst(i)
c  test for resetting of difference equations
      if(irsdif.ne.1) go to 10205
      ia=iimax
      ib=ikmax
      i=np1
      nu=nnqq
      ial=ik1
      j=in1
      k=in
      call resdif(x(np1),dx,y(1,n),y(1,np1),zk(kqs1),ap,aq,f,fd,
     .  v(iqs+1),phi(iqs+1),a,d,ia,ib,i,ial,nu,j,k)
c     construct d matrix (b block is already loaded)
10205 do 109 i=1,ii
      iv=v(i+iqs)
      if(ka.eq.0) go to 107
      do 104 m=i1,ik1
      wd=0.
      do 103 ial=1,ka
  103 wd=wd+a(i,ial)*q(ika+ial,m-ka)
  104 d(i,m)=d(i,m)+wd
c
c     construct beginning of e matrix and set into end of a
      if(ika.eq.0) go to 107
      do 106 j=1,ika
      wd=0.
      do 105 ial=1,ka
  105 wd=wd+a(i,ial)*q(ika+ial,j)
  106 a(i,j+ka)=wd+a(i,j+ka)
  107 continue
c
c      douglas is playing around here and we need to sort it out ]]
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
c     and store current value in complex q  for computation of e
c     in loops 103 and 105
c
  110 call leqsd(a(1,ka1),d(1,ka1),ii,ikka1,iimax,iimax,cerr,ii2,
     .  ikka1i)
      err=abs(cerr)
      if(err) 112,403,112
  112 det=det+log10(err)
      do 116 nu=ka1,ik1
      do 115 i=1,ii
      wd=-d(i,nu)
      q(i,nu-ka)=wd
  115 p(i+ipn)=wd
  116 ipn=ipn+ii
c  test size of ipn
      nrkwsa=ipn+ilngwk*lp
c*      write(6,11601) np1,ipn,nrkwsa
11601 format(' n, ipn, nrkwsa =',3i10)
      if(nrkwsp.gt.0.and.ipn+ilngwk*lp.gt.nrkwsp) go to 410
c
c     compute gamma matrix
  120 if(kip.eq.0.or.ika.eq.0) go to 125
      do 124 is=kab1,ik
      do 123 nu=ka1,ik1
      wd=0.
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
      wd=0.
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
c     *****************************
c
      if(nnqq.eq.nr) go to 150
      iorst(1,nnqq)=ii
      iorst(2,nnqq)=kk
      iorst(3,nnqq)=ka
      iorst(4,nnqq)=nn
      iorst(5,nnqq)=iqs
      iorst(6,nnqq)=kqs
      iorst(7,nnqq)=nqs
      iqsprv=iqs
      iqs=iqs+ii
      kqs=kqs+kk
      nqs=nqs+nn
      iiprv=ii
      kkprv=kk
      kaprv=ka
      ikprv=ik
      ikaprv=ika
c
c  end step in regions
c  *******************
c
  150 continue
c
      lpmax=ipn
      nvrtot=iqs+ii
c
c  set, and possibly print, actual workspace needed
      nrkwsa=lpmax+ilngwk*lp
c
      if(nrkwsp.gt.0) go to 200
      write(iw,1200) nrkwsa
      nrkwsp=nrkwsa
c
c     **********
c     remaining boundary conditions
c
  200 if(ka.eq.ik) go to 234
c  set boundary conditions at final boundary
      ia=ikmax2
      nu=nr+1
c
      do 2001 i=1,ia
      do 2001 j=1,ia
 2001 gs(i,j)=0.
c
      i=nqf
c
      call bc(x(nqf),x(nqf),y(1,nqf),y(1,nqf),zk(kqs1),zk(kqs1),ap,aq,
     .  g,gs,ia,i,km,i1,i2,i3,nu)
      do 2003 i=1,ik
      iv=i
      if(i.le.ii) iv=v(iqs+i)
      do 2003 j=1,ikka
 2003 gd(j+ka,i)=gs(j,iv)
      do 2005 j=1,ikka
 2005 gd(j+ka,ik1)=g(j)
c     compute coefficients for equation 12a and load into end of gd
      if(kp.eq.0) go to 210
      do 205 is=kab1,kap
      do 202 nu=ka1,ik1
      wd=0.
      do 201 ial=1,ka1
  201 wd=wd-gp(is,ial)*d(ika+ial,nu)
  202 gd(is,nu)=gd(is,nu)+wd
      do 203 ia=ka1,ii
  203 gd(is,ia)=gd(is,ia)+gp(is,ia)
  205 continue
  210 continue
c
c
c     compute coefficients for equation 13 and load into middle of gd
  220 if(ka.eq.0.or.ka1.gt.ik) go to 230
      do 223 m=ka1,ik
      do 222 nu=ka1,ik1
      wd=0.
      do 221 ial=1,ka
  221 wd=wd-gd(m,ial)*d(ika+ial,nu)
  222 gd(m,nu)=wd+gd(m,nu)
  223 continue
c
c     solve equations 12a and 13 (answer has wrong sign)
  230 if(ikka.eq.0) go to 240
      call leqsd(gd(ka1,ka1),gd(ka1,ik1),ikka,1,ikmax2,ikmax2,
     .  cerr,ikkaik,ikka)
      err=abs(cerr)
      if(err) 234,404,234
  234 det=det+log10(err)
  240 gd(ik1,ik1)=-1
      gd(ik1,1)=-1
c     (note that the first element of the second half of gp may now
c     have been overwritten)
c
c
c  final step through regions
c  **************************
c
      nr1=nr+1
      do 360 nnqq1=1,nr
      nnqq=nr1-nnqq1
      if(nnqq1.eq.1) go to 300
c  set order etc from array
      ii =iorst(1,nnqq)
      kk =iorst(2,nnqq)
      ka =iorst(3,nnqq)
      nn =iorst(4,nnqq)
      iqs=iorst(5,nnqq)
      kqs=iorst(6,nnqq)
      nqs=iorst(7,nnqq)
c
      nqf=nqs+nn-1
      ikt1=ik1
      kat1=ka1
      ik=ii+kk
      ik1=ik+1
      ikkat1=ikka1
      ka1=ka+1
      ikka1=ik1-ka
      ika=ii-ka
      ikka=ik-ka
c  set correction to solution and eigenvalues at last point in new
c  region
      do 252 nut=kat1,ikt1
  252 gd(nut,1)=gd(nut,in)
      ipn=ipn-ikka*ikkat1
      do 255 nu=ka1,ik
      wd=0.
      ia=ipn+nu-ka
      do 253 nut=kat1,ikt1
  253 wd=wd+p(ia+(nut-kat1)*ikka)*gd(nut,1)
  255 gd(nu,ik1)=wd
c
      gd(ik1,ik1)=-1
      gd(ik1,1)=-1
c
c     **********
c     iterated solution
c
c     empty ea,eb
  300 do 301 i=1,ii
      eb(i)=0.
      do 301 j=1,3
  301 ea(i,nnqq,j)=0.0
c
c     solution at second boundary
  310 if(kk.eq.0) go to 313
      do 312 k=1,kk
      j=k+kqs
      wd=gd(ii+k,ik1)
      gd(ii+k,1)=wd
      zk(j)=zk(j)-ucy*wd
  312 continue
  313 if(ka.ge.ii) go to 320
      do 318 ia=ka1,ii
      iva=v(ia+iqs)
      wa=abs(gd(ia,ik1))
      wd=y(iva,nqf)
  315 wd=wd-ucy*gd(ia,ik1)
      wb=abs(wd)
      ea(iva,nnqq,1)=wa
      eb(ia)=wb
      if(wb.lt.1.0e-10) go to 316
      ea(iva,nnqq,2)=wa/wb
      ea(iva,nnqq,3)= float(nqf)
  316 continue
      y(iva,nqf)=wd
  318 continue
c
c     remainder of solution   (beginning of second outer loop)
c
c     set storage indices
  320 in=1
      in1=ik1
c
      do 330 m=1,nn
      n=nqf+1-m
      ipn=ipn-ii*ikka1
      do 328 i=1,ii
      if(i.gt.ika) go to 321
      k=n-1
      if(k.lt.nqs) go to 328
      j=i+ka
      go to 322
  321 j=i-ika
      k=n
  322 wd=0.
      ia=i+ipn
      do 323 nu=ka1,ik1
  323 wd=wd+p(ia+(nu-ka-1)*ii)*gd(nu,in1)
      wa=abs(wd)
      jv=v(j+iqs)
      cwb=y(jv,k)-ucy*wd
      y(jv,k)=cwb
  325 wb=abs(cwb)
      ea(jv,nnqq,1)=ea(jv,nnqq,1)+wa
      eb(j)=eb(j)+wb
      if(wb.lt.1.0e-10) go to 326
      wb=wa/wb
  326 if(wb.lt.ea(jv,nnqq,2)) go to 327
      ea(jv,nnqq,2)=wb
      ea(jv,nnqq,3)= float(k)
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
c  test for zero solutions
      ebmax=0
      do 335 i=1,ii
      ebmax=max(ebmax,eb(i))
      if(eb(i).eq.0) eb(i)=1.e-35
  335 continue
c
      if(ebmax.eq.0) write(iw,1102)
c
c     set mean relative corrections
  340 do 341 i=1,ii
      iv=v(i+iqs)
  341 ea(iv,nnqq,1)=ea(iv,nnqq,1)/eb(i)
c
  360 continue
c
  380 return
c
c
c     diagnostics
c     **********
  400 write(iw,1001) nnqq
  401 do 402 k=1,2
      do 402 j=1,km
  402 write(iw,1002) (gs(j,i), i=1,ikikp1)
      write(iw,1101) (v(i+iqs), i=1,ii)
      go to 500
  403 write(iw,1003)n
      go to 500
  404 write(iw,1004)
c
c  reset boundary conditions at final boundary
c
      ia=ikmax2
      nu=nr+1
c
      do 405 i=1,ia
      do 405 j=1,ia
  405 gs(i,j)=0.
c
      i=nqf
c
      call bc(x(nqf),x(nqf),y(1,nqf),y(1,nqf),zk(kqs1),zk(kqs1),ap,aq,
     .  g,gs,ia,i,km,i1,i2,i3,nu)
      go to 401
  407 write(iw,1006) nnqq,nqf,nqs,x(nqf)
      go to 500
  408 write(iw,1007) nnqq,n,x(n),np1,x(np1),dx,r,nqs,nqf,
     *  (x(i), i=nqs,nqf)
      go to 500
c
  410 nrkwsa=ipn+ilngwk*lp
      write(iw,1008) nrkwsa,nrkwsp,np1,nnqq
c
  500 v(1)=0
      lpmax=ipn
      return
c
c
 1000 format(//1x,10('*'),10x,'improper formulation detected by nrkm',
     . 10x,10('*')/17x,'ii =',i3,4x,'kk =',i3,4x,'ka =',i3,4x,'km =',i3,
     .  '  iiprv =',i3,'  kkprv =',i3,'  kaprv =',i3/)
 1001 format(//1x,10('*'),5x,'first boundary condition matrix singular i
     .n nrkm, region no',i4,5x,10('*')//' gd ='/)
 1002 format(/1x,1p14e9.1)
 1003 format(//1x,10('*'),5x,'e matrix singular at n =',i4,' in nrk',
     . 5x,10('*'))
 1004 format(//1x,10('*'),5x,'second boundary matrix singular in nrk',
     .  5x,10('*')//' gd ='/)
 1005 format(' ')
 1006 format(//1x,10('*'),5x,'null range of independent variable in ',
     .       ' nrkm, region no.',i3,5x,10('*')//
     .       21x,'x(',i5,') = x(',i5,') =',1pe14.6/)
 1007 format(//1x,10('*'),5x,'independent variable not monotonic in',
     .       ' nrkm, region no.',i3,5x,10('*')//
     .       16x, 'x(',i4,') =',1pe12.4,
     .       ',    x(',i4,') =',1pe12.4/16x,'difference =',1pe12.4,
     .       ',    range =',1pe12.4//1x,'x(n), n=',i5,',',i5/
     .       (1x,1p10e13.5))
c
 1008 format(//1x,10(1h*),'  workspace needed in nrkm is more than',
     .  i10,' 8 byte words and exceeds workspace provided =',
     .  i10,' 8 byte words'
     .  //11x,' this happens at point no.',i4,
     .  ' in x, in region no.',i4)
c     warning message
 1100 format(1x,5('*'),3x,'v set in nrk        v(',i2,') =',i3)
 1101 format(/9x,'v =',7(1x,5i3))
 1102 format(//1x,5(1h*),' solution is identically zero in nrkm')
c
c
 1200 format(///' size of common/work/ not checked in s/r nrkm.',
     .  i10,' 8 byte words required.')
c
      end
c..      subroutine leqsd(a,b,n,m,ia,ib,det,isa,isb)
c..      implicit double precision (a-h,o-z)
c..      dimension a(ia,1),b(ib,1)
c..c
c..      call leq(a,b,n,m,ia,ib,det)
c..      return
c..      end
