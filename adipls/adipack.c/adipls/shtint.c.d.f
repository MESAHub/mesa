      subroutine shtint(x,y,iy,nw1,nibc,nn,nev1,nd,ny,iasn,ii,ig,
     *  nfit,mdintg,det,ddsol,ds1,iord,icry)
c
c  uses linint or eigint integration from centre and surface
c  iasn gives absolute value of step in model. When doing
c  Richardson extrapolation, call initially with iasn = 2 and
c  subsequently with iasn = 1.
c
c  sets matching determinant det, scaled determinant ddsol and,
c  if mdintg = 2, order of mode iord
c
c  original version 13/2/89
c
c  Modified 28/6/95, removing nfit1 as argument. Moving iprdet to
c  common /cdiagn/
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      dimension x(1),y(iy,1),ds1(4,5),ds(4,5)
      common/rhsdat/ el,ell,alb,els,el1,sig,anres,perfac,data(8)
      common/cdiagn/ idgrhs, idgrh1, iprdet
      common/cincnt/ xmnevn,xfit, fcnorm, eps, epssol, dsigmx, fsig, 
     *  dsigre,icow,iturpr,iplneq,iriche,nnwwin, moddet, itmax, irsevn, 
     *  nftmax,irsord,iekinr,itsord,imissl,imjssl,imstsl
      common/csexpc/ y11, y13, y21, y23, y31, y33, y41, y43,
     *  yt11, yt13, yt21, yt23, yt31, yt33, yt41, yt43, 
     *  yt241, yt243, ytt41, ytt43 
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      external adirhs
c
      data iopen/0/
c
c  initial points and numbers of integration steps (only used when
c  mdintg = 1 or 2)
c
c  hardwired flag for diagnostics output
c
      idiag = 0
c
      if(idiag.gt.0.and.iopen.eq.0) then
	iopen=1
	open(51,file='ttt.shtint.out',status='unknown')
      end if
      if(idiag.gt.0) write(51,69091) el,sig,ell,data(7),
     *  y11, y13, y21, y23, y31, y33, y41, y43,
     *  yt11, yt13, yt21, yt23, yt31, yt33, yt41, yt43, 
     *  yt241, yt243, ytt41, ytt43 
69091 format(1p4e15.7)
c
      icry=0
c
      if(iasn.gt.1) then
        if(mod(nev1-nibc,iasn).ne.0) then
          write(istdou,110) nev1, nibc, iasn
          if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,110) 
     *      nev1, nibc, iasn
          icry = -1
          return
        else if(mod(nfit-nev1,iasn).ne.0) then
          write(istdou,115) nfit, nev1, iasn
          if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,110) 
     *      nfit, nev1, iasn
          icry = -1
          return
        else if(mod(nd-nfit,iasn).ne.0) then
          write(istdou,120) nd, nfit, iasn
          if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,120) 
     *      nd, nfit, iasn
          icry = -1
          return
        end if
      end if
c
      nd1=nibc
      nni1=(nev1-nd1)/iasn+1
      nd2=nev1
      nni2=(nfit-nd2)/iasn+1
      nnis=(nd-nfit)/iasn+1
c
      iord=0
c
c  integrate from centre
c
      isn=iasn
      anres=0
c
      if(nni1.gt.1) then
        if(iplneq.eq.1) then
          el1=el
        else
          el1=el-1
        end if
c
        call geninh(x(nibc),y(1,nibc),adirhs,ii,iy,ig,isn,nd1,nni1,
     *    mdintg,iclcn)
c..	write(6,*) 'end first integration'
c  contribution to order from inner region (only applicable for
c  mdintg = 2)
        iord=iclcn
c
      end if
c
c  integrate in intermediate region
c
      if(nni2.gt.1) then
        el1=0
        call geninh(x(nev1),y(1,nev1),adirhs,ii,iy,ig,isn,nd2,nni2,
     *    mdintg,iclcn)
c..	write(6,*) 'end second integration'
c  contribution to order from this region
        iord=iord+iclcn
c
      end if
c
c  integrate from surface
c
      isn=-iasn
c
      el1=0
c
c  test for taking out t**n from y4 for n .lt. 1
c
      if(data(7).lt.1.and.data(7).gt.0.and.ii.eq.4
     *  .and.mdintg.eq.1) then
	anres=data(7)
      else
	anres=0
      end if
c
c  the integration
c
      call geninh(x(nd),y(1,ny),adirhs,ii,iy,ig,isn,nd,nnis,
     *  mdintg,iclcn)
c
c  for simplicity, reset y4 here, if term in t**n has been taken out
c
      if(anres.ne.0) then
	do 25 n=1,nnis
	n1=nd+1-n
	n2=ny+1-n
	t=1-x(n1)
	tn=t**anres
c..        write(6,*) t,y(4,n2),yt41*tn,y(8,n2),yt43*tn
        y(4,n2)=yt41*tn+y(4,n2)
   25   y(8,n2)=yt43*tn+y(8,n2)
c
      end if
c
c  contribution from outer region to order
c
      iig=ii*ig
      if(idiag.gt.0) then
	write(51,*) iig
        do 27 n=1,nnis
        t=1.d0-x(nd+2-n)
   27   write(51,69092) t,(y(i,ny+2-n),i=1,iig)
69092   format(1p9e15.7)
c
      end if
c
      iord=iord+iclcn
c
c  for mdintg=2, test for correcting order for contribution
c  at innermost meshpoint
c
      if(mdintg.eq.2.and.y(1,nibc)*y(2,nibc).lt.0) iord=iord+1
c
c  calculate matching determinant
c  ******************************
c
      if(ig.eq.1) then
c
c  cowling approximation or radial case
c
        n=nfit
        n1=nfit+1
        det=y(1,n)*y(2,n1)-y(2,n)*y(1,n1)
c  print determinant?
        if(iprdet.eq.1.and.istdpr.gt.0) write(istdpr,130) 
     *    ((y(i,j),i=1,2),j=n,n1)
c  norm of determinant
        detnrm=danorm(y(1,n),2)*danorm(y(1,n1),2)
c
      else
c
c  full set
c
        j=0
        do 30 n=nfit,nfit+1
        l=0
        do 30 k=1,2
        j=j+1
        do 30 i=1,4
        l=l+1
        ds(i,j)=y(l,n)
   30   ds1(i,j)=ds(i,j)
c  print determinant?
        if(iprdet.eq.1.and.istdpr.gt.0) write(istdpr,140) 
     *    ((ds1(i,j),j=1,4),i=1,4)
c  norm of determinant
        detnrm=1
        do 35 j=1,4
   35   detnrm=detnrm*danorm(ds(1,j),4)
c
        call leq(ds,ds(1,5),4,0,4,4,det)
c
      end if
c
      if(detnrm.eq.0) then
        write(istdou,150) 
        if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,150) 
        icry=-1
      else
        ddsol=abs(det)/detnrm
      end if
      return
c
  110 format(//' **** error in s/r shtint. nev1, nibc =',2i5,
     *  '  iasn =',i3)
  115 format(//' **** error in s/r shtint. nfit, nev1 =',2i5,
     *  '  iasn =',i3)
  120 format(//' **** error in s/r shtint. nd, nfit =',2i5,
     *  '  iasn =',i3)
  130 format(//'  matching determinant:'//(1p2e16.8))
  140 format(//'  matching determinant:'//(1p4e16.8))
  150 format(//' **** error in s/r shtint. detnrm = 0')
      end
