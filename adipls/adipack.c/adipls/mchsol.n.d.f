      subroutine mchsol(x,y,iy,aa,iaa,nw1,nibc,nn,iasn,ii,ig,
     *    nfit,npout,imissl,imjssl,imstsl,det,ds1,icry)
c
c  find matched solution, after shooting method integration
c
c  original version: 14/2/89
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Modified 28/6/95 to remove nfit1 as an argument, and include
c  moddet in common/cincnt/
c
      implicit double precision (a-h,o-z)
      dimension x(1),y(iy,1),aa(iaa,1),ds1(4,5),c(4)
      common/cincnt/ xmnevn,xfit, fcnorm, eps, epssol, dsigmx, fsig, 
     *  dsigre,icow,iturpr,iplneq,iriche,nnwwin, moddet, itmax, irsevn, 
     *  nftmax,irsord,iekinr,itsord
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c..      write(6,*) 'Enter mchsol with iasn, npout =',iasn,npout
c
c  find coefficients
c  *****************
c
      if(ig.eq.1) then
c
c  cowling approximation or radial case
c
        if=1
        if(y(1,nfit).eq.0) if=2
        if(y(if,nfit).eq.0) then
c  zero solution at fitting point. write diagnostic
          nwd1=nfit-20
          nwd2=nfit+20
          write(istdou,132) (n,x(n),(y(i,n),i=1,2),n=nwd1,nwd2)
          if(istdpr.ne.istdou.and.istdpr.gt.0)
     *      write(istdpr,132) (n,x(n),(y(i,n),i=1,2),n=nwd1,nwd2)
          go to 80
        end if
c
        c(1)=y(if,nfit+1)/y(if,nfit)
c
c  output of individual solutions
c
        if(npout.gt.0.and.istdpr.gt.0) then
          nnw=nn-nw1+1
          ndpout=iasn*max0(1,nnw/npout)
          write(istdpr,135) nw1,x(nw1),(y(i,nw1),i=1,2)
          if(nw1.lt.nibc) then
            nw11=nibc
          else
            nw11=nw1+ndpout
          end if
          write(istdpr,132) (n,x(n),(y(i,n),i=1,2),
     *      n=nw11,nfit,ndpout)
          write(istdpr,132) (n,x(n),(y(i,n+1),i=1,2),
     *      n=nfit,nn,ndpout)
        end if
c
c  find combined solution
c
        do 20 n=nw1,nfit
        do 20 i=1,2
   20   y(i,n)=c(1)*y(i,n)
        do 25 n=nfit+1,nn
        n1=n+1
        do 25 i=1,2
   25   y(i,n)=y(i,n1)
c
      else
c
c  full set
c
        call mchcff(ds1,c,imissl,imjssl,imstsl,det)
c
        if(det.lt.0) then
c
c  degenerate solution at fitting point. write diagnostic
c
          nwd1=nfit-20
          nwd2=nfit+20
          write(istdou,120) (n,x(n),(y(i,n),i=1,8),n=nwd1,nwd2)
          if(istdpr.ne.istdou.and.istdpr.gt.0)
     *      write(istdpr,120) (n,x(n),(y(i,n),i=1,8),n=nwd1,nwd2)
          go to 80
        end if
c
c  analyze determinant?
c
        if(moddet.eq.1) call anldet(ds1,c,x,y,iy,aa,iaa,nn,nfit)
c
        if(istdpr.gt.0) write(istdpr,130) c
c
c  output of individual solutions
c
        if(npout.gt.0) then
          nnw=nn-nw1+1
          ndpout=iasn*max0(1,nnw/npout)
          if(istdpr.gt.0) write(istdpr,140) nw1,x(nw1),(y(i,nw1),i=1,8)
          if(nw1.lt.nibc) then
            nw11=nibc
          else
            nw11=nw1+ndpout
          end if
          if(istdpr.gt.0) write(istdpr,145) (n,x(n),(y(i,n),i=1,8),
     *      n=nw11,nfit,ndpout)
          if(istdpr.gt.0) write(istdpr,145) (n,x(n),(y(i,n+1),i=1,8),
     *      n=nfit,nn,ndpout)
        end if
c
c  find combined solution
c
        fccan=0
        do 30 n=nw1,nfit
        do 30 i=1,4
        yy1=c(1)*y(i,n)
        yy2=c(2)*y(i+4,n)
        if(n.gt.1) fccan=fccan+abs(yy1+yy2)/(abs(yy1)+abs(yy2)+1.d-10)
   30   y(i,n)=yy1+yy2
c
        do 35 n=nfit+1,nn
        n1=n+1
        do 35 i=1,4
        yy1=-c(3)*y(i,n1)
        yy2=-c(4)*y(i+4,n1)
        y(i,n)=yy1+yy2
   35   fccan=fccan+abs(yy1+yy2)/(abs(yy1)+abs(yy2)+1.d-10)
c
        fccan=fccan/(4*nn)
c
c  print fccan
c
        if(fccan.le.1.d-5) then
          write(istdou,150) fccan
          if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,150) fccan
        else
          if(istdpr.gt.0) write(istdpr,155) fccan
        end if
c
      end if
c
      icry = 0
c
      return
c
c  exit with errors
c
   80 icry = -1
      return
  110 format(///1x,10(1h*),' degeneracy at fitting point. solution:'//
     *  (i5,0pf12.7,1p2e13.5))
  120 format(///1x,10(1h*),' degeneracy at fitting point. solution:'//
     *  (i5,0pf12.7,1p8e13.5))
  130 format(//'  c:',1p4e13.5)
  132 format(/(i5,0pf10.7,1p2e11.3))
  135 format(//'  n,x,y(1-2):'//(i5,0pf10.7,1p2e11.3))
  140 format(//'  n,x,y(1-8):'//(i5,0pf10.7,1p8e11.3))
  145 format(/(i5,0pf10.7,1p8e11.3))
  150 format(//1x,10(1h*),' *** warning. cancellation factor =',1pe13.5,
     *  '  is dangerously small')
  155 format(/' cancellation factor =',1pe13.5)
      end
