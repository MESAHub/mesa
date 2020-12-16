      subroutine setsum(isolcv,dat)
c  sets summary in dat from data in common/csumma/
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      dimension dat(*)
      common/csumma/ xmod1,data(8),datmd1(2),xtrnct,frfit,xfit,
     *  fctsbc,fcttbc,albsum,
     *  elsum,ordsum,sigsum,sigc,dmx,xmx,ekin,per,perv,frqv,
     *  ddsig,ddsol,ysum(4),dosdum(5),
     *  in,nn,mdintg,ivarf,icase,iorign,ifilec,nfmodc,mlname(4)
      common/cicode/ icow,iper,intmod,ialb,istsbc,icc(4)
c
      save
c  decode icase
      call decicd(icase,iel,iekin)
c
      i=0
      if(albsum.eq.1) go to 10
      i=1
      dat(i)=albsum
      go to 15
c
   10 if(fctsbc.eq.0) go to 15
      i=1
      dat(2)=fctsbc
c
   15 i=i+1
      dat(i)=elsum
      dat(i+1)=ordsum
      i=i+2
      dat(i)=sigsum
c  cowling approximation?
      if(icow.eq.0) go to 20
      i=i+1
      dat(i)=sigc
c
   20 dat(i+1)=dmx
      dat(i+2)=xmx
      ekn=ekin
      if(isolcv.eq.-1) ekn=-ekn
      dat(i+3)=ekn
      dat(i+4)=per
c  periods?
      if(iper.ne.1) return
      dat(i+5)=perv
      dat(i+6)=frqv
      return
      end
c
c
