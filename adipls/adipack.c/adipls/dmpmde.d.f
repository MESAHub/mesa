      subroutine dmpmde(ids,iw,k)
c
c  dump summary from adiabatic mode file on d/s ids in standard form
c  (for k = 1) or completely (for k = 2), on d/s iw
c  if ids.lt.0 s/r dumps contents of common/csumma/, according to
c  the value of k
c
c  ....................................................................
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      logical aeq,smdat
      dimension gs(50),mlnamp(4),dat(11),datp(8)
      common/csumma/ xmod1,data(8),datmd1(2),xtrnct,frfit,xfit,
     *  fctsbc,fcttbc,albsum,
     *  elsum,ordsum,sigsum,sigc,dmx,xmx,ekin,per,perv,frqv,
     *  ddsig,ddsol,ysum(4),dosdum(5),
     *  in,nn,mdintg,ivarf,icase,iorign,ifilec,nfmodc,mlname(4)
      common/work/ xy(7,1)
      equivalence (gs(1),xmod1)
c
      save
c
      data icasep,mlnamp,xtrncp /-1,4*8h        ,0./
      data datp /8*0./
c
      if(ids.ge.0) go to 5
      go to (1,2), k
c
    1 nm=2
      go to 12
c
    2 write(iw,117)
      nm=0
      go to 30
c
    5 rewind ids
      write(iw,100) ids
      nm=0
c
c
   10 nm=nm+1
      read(ids,end=90,err=90) gs,nnmd,(xy(1,n),(xy(i,n),i=2,7),
     *   n=1,nnmd)
      go to (12,30), k
c  standard summary
   12 call cmparr(mlname,mlnamp,4,aeq,1)
c  test for change in data
      smdat=.true.
      do 13 i=1,8
      smdat=smdat.and.(data(i).eq.datp(i))
   13 datp(i)=data(i)
      if(nm.gt.1.and.smdat) go to 15
c
      if(nm.gt.1) write(iw,104)
      write(iw,105) (data(i),i=1,4)
c
      ih=1
      nwmod=1
      write(iw,110) mlname
      go to 20
c  test for output of name or heading
   15 if(aeq) go to 17
      write(iw,110) mlname
      nwmod=1
      go to 20
c
   17 if(icase.ne.icasep) ih=1
      if(xtrnct.ne.xtrncp) nwmod=1
c  set summary
   20 isolcv=1
      if(ddsol.gt.0.05) isolcv=-1
      call setsum(isolcv,dat)
c  output summary
      call sumout(iw,ih,nwmod,data(7),xtrnct,icase,dat,1,1,nm)
   25 nwmod=0
      ih=0
      icasep=icase
      xtrncp=xtrnct
c
      if(ids.ge.0) go to 10
c
      return
c
c  complete summary
c
   30 write(iw,120) nm,mlname,(gs(i),i=1,11),in,nn,(gs(i),i=12,17),
     *  mdintg,ivarf,icase,(gs(i),i=18,38),iorign,ifilec,
     *  nfmodc
      if(ids.ge.0) go to 10
c
   90 continue
      return
  100 format(1h1,' modes  on d/s',i4//)
  104 format(//40x,50(1h*))
  105 format(/' new model. mass =',1pe13.5,'  radius =',e13.5,
     *  '  pc =',e13.5,'  rhoc =',e13.5)
  110 format(/' name of model: ',4a8)
  117 format(/' dump of common/csumma/:')
  120 format(/' mode no',i5/1x,4a8/1p11e11.3/2i11,6e11.3,3i11/
     *  10e11.3/11e11.3/3i11)
      end
c
c
