      subroutine rnmean(a,b,nn,ia,ib,imean)
c
c  sets gaussian weighted running mean over imean points
c  (incremented by 1 if even on input)
c  Weight is defined to be 0.05 at outermost point in mean
c
c  Original version: 30/5/92
c
      implicit double precision(a-h, o-z)
      parameter(iwmax=101)
      dimension a(ia,1), b(ib,1), w(iwmax)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  test for odd imean
c
      if(mod(imean,2).eq.0) then
	imean = imean+1
	if(istdpr.gt.0) write(istdpr,120) imean
      end if
c
      imh=(imean-1)/2
      do 15 i=1,imean
      xx=float(i-imh-1)/float(imh)
   15 w(i)=exp(-3.*xx*xx)
c
c  now set running mean
c
      do 30 n=1,nn
      irange=min(imh,n-1,nn-n)
      n1=n-irange
      n2=n+irange
      sum=0
      wsum=0
      i=imh-irange
      do 25 k=n1,n2
      i=i+1
      sum=sum+w(i)*a(1,k)
   25 wsum=wsum+w(i)
   30 b(1,n)=sum/wsum
      return
  120 format(/' **** Warning in s/r rnmean. imean increased to',i4)
      end
