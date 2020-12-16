      subroutine mchcff(ds,c,is,js,isetsl,err)
c
c  computs matching coefficients c from almost singular
c  matrix ds.
c  if isetsl .ne. 1 equation is and variable js are eliminated.
c  otherwise all possibilities are tested, and the one minimizing
c  the error is selected.
c
c  original version: 2/4/1985
c
c  modified 11/7/1985 to be single precision throughout.
c
c  modified 13/8/87 to standardize output
c
c  ..............................................................
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      dimension ds(4,4),c1(4),c(4)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      if(isetsl.eq.1) go to 20
c
      call mchcf1(ds,c,is,js,err)
      return
c
c  step through is, js
c
   20 ismin=0
c
      do 30 is1=1,4
      do 30 js1=1,4
c
      call mchcf1(ds,c1,is1,js1,err)
c
c  test for degeneracy
c
      if(err.lt.0) go to 30
c
      if(err.gt.errmin.and.ismin.gt.0) go to 30
      ismin=is1
      jsmin=js1
      errmin=err
      do 25 i=1,4
   25 c(i)=c1(i)
c
   30 continue
c
c  test for no solution found
c
      if(ismin.eq.0) return
c
c  set is, js to values giving minimum, print diagnostics
c
      is=ismin
      js=jsmin
      err=errmin
      if(istdpr.gt.0) write(istdpr,130) is,js
      return
  120 format(/' is, js, err =',2i4,1pe13.5)
  130 format(//' error minimizing is, js =',2i5)
      end
      subroutine mchcf1(ds,c,is,js,err)
c
c  finds solution to almost singular matrix ds, with
c  maximum element normalized to 1.
c  err returns rms error
c
c  original version: 2/4/1985
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      dimension ds(4,4),ds1(4,4),c(4),ss(4)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      do 20 i=1,4
      do 20 j=1,4
   20 ds1(i,j)=ds(i,j)
c
      do 25 i=1,4
      ds1(is,i)=0
      ds1(i,js)=0
   25 c(i)=-ds(i,js)
c
      ds1(is,js)=1
      c(js)=0
c
      call leq(ds1,c,4,1,4,4,det)
c
c  test for no solution
c
      if(det.eq.0) go to 80
c
      c(js)=1
c
c  set maximum to 1
c
      cmax=0
      do 30 i=1,4
      ca=abs(c(i))
      if(ca.lt.cmax) go to 30
        imax=i
        cmax=ca
   30 continue
c
      cimax=1./c(imax)
      do 35 i=1,4
   35 c(i)=cimax*c(i)
c
c  test solution
c
      err=0
      do 45 i=1,4
      ss(i)=0
      do 42 j=1,4
   42 ss(i)=ss(i)+ds(i,j)*c(j)
c
   45 err=err+ss(i)*ss(i)
c
      err=sqrt(err)
      return
c
c  diagnostics for zero determinant
c
   80 write(istdou,100) is,js
      if(istdpr.ne.istdou) write(istdpr,100) is,js
      err=-1
      return
  100 format(//' ***** in s/r mchcf1 determinant is zero',
     *  ' for is, js =',2i5)
      end
c
c
