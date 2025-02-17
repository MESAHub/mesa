      subroutinelineexpop(linelist,readdump,nwave,mxwave,hfreq3c2,dlnfreq,wlgrd,scatop,opacity,emis,plnkfnc,init,mkshtlst,taumin,lon
     *glist,nradii,nrdim,eden,temp,dvdravg,iondenz,ionstg1,nr1,nr2)
      implicitreal*8(a-h,o-z)
      PARAMETER(MAXIONM1=(6-1))
      PARAMETER(NSBINTVL=3*10)
      PARAMETER(lindim=300000)
      datapi/3.141592653589793d+00/,fourpi/12.5637061d+00/
      databc/1.380626d-16/,h/6.626205d-27/,hev/4.1357d-15/
      datac/2.997925d+10/,elecxsec/6.6524d-25/
      dataevtoerg/1.6022d-12/,a/7.56464d-15/
      datastefbltz/5.66956d-05/,hmass/1.67352d-24/
      dataesu/4.80298d-10/,emass/9.1091d-28/
      datasrpi/1.772453851d+00/,bcev/8.617064d-05/
      character*80linelist,longlist
      logicalreaddump,mkshtlst
      integerionstg1(nrdim,99)
      logicalwrtdump,fexist
      integernwave,mxwave,nradii,nrdim,init,nr1,nr2
      real*8wlgrd(mxwave+1),taumin
      real*8dlnfreq(mxwave),hfreq3c2(mxwave)
      real*8scatop(mxwave,nradii),opacity(mxwave,nradii)
      real*8emis(mxwave,nradii)
      real*8plnkfnc(mxwave,nradii)
      real*8eden(nradii),temp(nradii),dvdravg(nradii)
      real*8iondenz(nrdim,6,99)
      real*8expvecx(75*200)
      real*8lscat(10000),lopac(10000)
      real*8lemis(10000)
      real*8idenz(0:31*99*75)
      integerlinecountr(4950)
      integerlinecountw(4950)
      real*8enrvec(200)
      real*4wavelen0(lindim),gf0(lindim)
      real*4elower0(lindim),eupper0(lindim)
      real*4element0(lindim)
      integerznucline0(lindim),ionstage0(lindim)
      real*4wavelen(lindim),gf(lindim)
      real*4elower(lindim),eupper(lindim)
      real*4enrwt(lindim)
      real*4dumvec(lindim),element(lindim)
      real*8dvdrlog(75),hckt(75)
      logicallabund(99)
      integeridumvec(lindim),levindex(lindim)
      integerznucline(lindim),ionstage(lindim)
      integerindexv(lindim),lteflgx(lindim)
      integerenrindex(lindim),wlindex(lindim)
      real*4jupper,jlower,tenlog
      character*4ref
      character*8lnamelw,lnameup
      integerptrinit
      common/ltelnz/ptrwavln,pelement,pidmvec,pdumvec,ptznucln,ptenrwt,ptenrvec,ptgf,ptelow,ptenrndx,ptionstg,ptwlindx,ptrlvndx,nlin
     *es
      real*8quench_fac(2)
      integerdiscard1,discard2,type1,type2,type3,type4
      saveptrinit
      realtm1,tm2
      datac18/2.997925d+18/,ptrinit/0/,smallnum/1.d-70/
      dataediffmax/1.d-01/
      dataquench_fac/2.9d+17,3.4d+16/
      iiii=1
      saha=0.5d+00*((h/emass)*(h/bc)/(2.d+00*pi))**1.5
      cross0=pi*(esu/emass)*(esu/c)*1.d-08
      doiz=1,99
      donr=1,nradii
      idenz(nr+(1-1)*nradii+(iz-1)*31*nradii)=0.d+00
      enddo
      doi=2,31
      calldcopy(nradii,idenz(1+(i-1-1)*nradii+(iz-1)*31*nradii),1,idenz(1+(i-1)*nradii+(iz-1)*31*nradii),1)
      enddo
      donr=1,nradii
      doii=1,6
      i=ionstg1(nr,iz)+ii-1
      idenz(nr+(i-1)*nradii+(iz-1)*31*nradii)=iondenz(nr,ii,iz)
      enddo
      enddo
      enddo
      donr=1,nradii
      hckt(nr)=h*c/(bc*temp(nr))
      enddo
      wrtdump=(.not.readdump)
      if(mkshtlst)then
      inquire(file=longlist,exist=fexist)
      print*,'LongList=',longlist
      print*,'linelist=',linelist
      if(.not.fexist)then
      write(7,'(3a)')' lineexpop: unable to',' find file ',linelist(1:lnblnk(linelist))
      stop
      endif
      open(file=longlist,unit=20,status='old')
      open(file=linelist,unit=21,status='unknown')
      tsminln=dlog(cross0*10.d+00/taumin)
      do3000iz=1,99
      do2620i=1,min(iz,6)
      do2610nr=1,nradii
      idenz(nr+(i-1)*nradii+(iz-1)*31*nradii)=dlog(max(idenz(NR+(I-1)*nradii+(IZ-1)*31*nradii),1.d-99))
2610  continue
2620  continue
3000  continue
      do3010nr=1,nradii
      dvdrlog(nr)=dlog(dabs(dvdravg(nr)))
3010  continue
      dol=1,4950
      linecountr(l)=0
      enddo
      dol=1,4950
      linecountw(l)=0
      enddo
      nlcount=0
      do3130il=1,9999999
      read(20,30,end=3140)wl,ref,gflog,jlower,el,jupper,eu,code,lnamelw,lnameup
      iz=int(code+0.001)
      i=int(100.d+00*(code-dble(iz))+0.1d+00)+1
      linecountr(i+((iz-1)*iz)/2)=linecountr(i+((iz-1)*iz)/2)+1
      el=dabs(el)
      eu=dabs(eu)
      if(el.gt.eu)then
      oldeu=eu
      eu=el
      el=oldeu
      oldju=jupper
      jupper=jlower
      jlower=oldju
      endif
      wllog=dlog(wl)
      denmaxln=idenz(1+(i-1)*nradii+(iz-1)*31*nradii)-el*hckt(1)-dvdrlog(1)
      donr=2,nradii
      denmaxln=max(denmaxln,idenz(nr+(i-1)*nradii+(iz-1)*31*nradii)-el*hckt(nr)-dvdrlog(nr))
      enddo
      if((tsminln+denmaxln+gflog+wllog).gt.0.d+00)then
      nlcount=nlcount+1
      write(21,30)wl,ref,gflog,jlower,el,jupper,eu,code,lnamelw,lnameup
      linecountw(i+((iz-1)*iz)/2)=linecountw(i+((iz-1)*iz)/2)+1
      endif
      wlo=wl
      donl=1,300
      read(20,30,end=3140)wl,ref,gflog,jlower,el,jupper,eu,code,lnamelw,lnameup
      dlnwl=2.d+00*(wl-wlo)/(wl+wlo)
      wllog=wllog+dlnwl
      wlo=wl
      iz=int(code+0.001)
      i=int(100.d+00*(code-dble(iz))+0.1d+00)+1
      denmaxln=idenz(1+(i-1)*nradii+(iz-1)*31*nradii)-el*hckt(1)-dvdrlog(1)
      donr=2,nradii
      denmaxln=max(denmaxln,idenz(nr+(i-1)*nradii+(iz-1)*31*nradii)-el*hckt(nr)-dvdrlog(nr))
      enddo
      if((tsminln+denmaxln+gflog+wllog).gt.0.d+00)then
      nlcount=nlcount+1
      write(21,30)wl,ref,gflog,jlower,el,jupper,eu,code,lnamelw,lnameup
      endif
      enddo
3130  continue
3140  continue
      write(7,'(a,a,a,a,a,i8)')' Number of lines written from ',longlist(1:lnblnk(longlist)),' to ',linelist(1:lnblnk(linelist)),' e
     *quals ',nlcount
      write(7,*)
      write(7,'(a)')' Lines per ionization stage read and written:'
      write(7,'(a)')'    z    i      read   written'
      doiz=1,99
      doi=1,iz
      if(linecountr(i+((iz-1)*iz)/2).gt.0)then
      write(7,'(1x,i4,1x,i4,1x,i9,1x,i9)')iz,i,linecountr(i+((iz-1)*iz)/2),linecountw(i+((iz-1)*iz)/2)
      endif
      enddo
      enddo
      write(7,*)
      stop
      endif
      donr=1,nradii
      lscat(nr+(1-1)*nradii)=0.d+00
      lopac(nr+(1-1)*nradii)=0.d+00
      enddo
      donf=2,nwave
      calldcopy(nradii,lscat(1+(nf-1-1)*nradii),1,lscat(1+(nf-1)*nradii),1)
      calldcopy(nradii,lopac(1+(nf-1-1)*nradii),1,lopac(1+(nf-1)*nradii),1)
      enddo
      if(init.eq.0)then
      if(readdump)then
      inquire(file=linelist,exist=fexist)
      print*,' Linelist ',linelist
      if(.not.fexist)then
      write(7,'(3a)')' lineexpop: readdump=true: unable to',' find file: ',linelist(1:lnblnk(linelist))
      print*,' program stoped'
      stop
      endif
      print*,' program continues , file exists ',linelist
      open(file=linelist,form='unformatted',unit=35,status='old',action='read')
      read(35)nlines0
      callreadbinr(35,elower0,nlines0)
      callreadbinr(35,eupper0,nlines0)
      callreadbinr(35,gf0,nlines0)
      callreadbini(35,ionstage0,nlines0)
      callreadbinr(35,wavelen0,nlines0)
      callreadbini(35,znucline0,nlines0)
      callreadbinr(35,element0,nlines0)
      close(unit=35)
      donl=1,nlines0
      elower0(nl)=abs(elower0(nl))
      eupper0(nl)=abs(eupper0(nl))
      enddo
      else
      inquire(file=linelist,exist=fexist)
      if(.not.fexist)then
      write(7,'(3a)')' lineexpop: unable to',' find file ',linelist(1:lnblnk(linelist))
      stop
      endif
      open(file=linelist,unit=20,status='old')
19    format(i10)
      read(20,19)nlines0
      do20nl=1,nlines0
      read(20,30)wavelen0(nl),ref,gf0(nl),jlower,elower0(nl),jupper,eupper0(nl),element0(nl),lnamelw,lnameup
30    format(f10.4,1x,a4,1x,f7.3,f4.1,f11.3,f4.1,1x,f11.3,f7.2,a8,2x,a8)
      znucline0(nl)=int(element0(nl)+0.001)
      ionstage0(nl)=int(100.d+00*(element0(nl)-dble(znucline0(nl)))+0.1d+00)+1
      eupper0(nl)=abs(eupper0(nl))
      elower0(nl)=abs(elower0(nl))
20    continue
      close(20)
      if(wrtdump)then
      open(file='linedata.dump',form='unformatted',unit=35,status='unknown')
      write(35)nlines0
      callwritbinr(35,elower0,nlines0)
      callwritbinr(35,eupper0,nlines0)
      callwritbinr(35,gf0,nlines0)
      callwritbini(35,ionstage0,nlines0)
      callwritbinr(35,wavelen0,nlines0)
      callwritbini(35,znucline0,nlines0)
      callwritbinr(35,element0,nlines0)
      close(unit=35)
      endif
      endif
      write(7,*)' lineexpop: number of lines read in = ',nlines0
      write(7,*)' lineexpop: uncomment write(7,*) for debugging '
      enrmax=0.d+00
      donl=1,nlines0
      if(elower0(nl).ge.enrmax)enrmax=elower0(nl)+1.d+00
      enddo
      tenlog=log(10.e+00)
      donl=1,nlines0
      gf0(nl)=exp(tenlog*gf0(nl))
      wavelen0(nl)=10.e+00*wavelen0(nl)
      enddo
      endif
      nlines=nlines0
      donl=1,nlines
      elower(nl)=elower0(nl)
      eupper(nl)=eupper0(nl)
      gf(nl)=gf0(nl)
      ionstage(nl)=ionstage0(nl)
      wavelen(nl)=wavelen0(nl)
      znucline(nl)=znucline0(nl)
      element(nl)=element0(nl)
      enddo
      wlmin=wlgrd(nwave+1)
      wlmax=wlgrd(1)
      if(dble(wavelen(1)).gt.wlmin)then
      nl1=1
      else
      nl1=ndexs(sngl(wlmin),wavelen,nlines)+1
      endif
      if(dble(wavelen(nlines)).lt.wlmax)then
      nl2=nlines
      else
      nl2=ndexs(sngl(wlmax),wavelen,nlines)
      endif
      if(nl1.gt.1)then
      do45nnl=nl1,nl2
      nl=nnl-nl1+1
      elower(nl)=elower(nnl)
      eupper(nl)=eupper(nnl)
      gf(nl)=gf(nnl)
      ionstage(nl)=ionstage(nnl)
      wavelen(nl)=wavelen(nnl)
      znucline(nl)=znucline(nnl)
      element(nl)=element(nnl)
45    continue
      nlines=nl2-nl1+1
      else
      if(nl2.lt.nlines)then
      nlines=nl2-nl1+1
      endif
      endif
      doiz=1,99
      labund(iz)=.false.
      doii=1,6
      labund(iz)=(labund(iz).or.(dasum(nradii,iondenz(1,ii,iz),1).gt.1.d-10))
      enddo
      enddo
      donl=1,nlines
      iz=znucline(nl)
      if(labund(iz))then
      lteflgx(nl)=0
      else
      lteflgx(nl)=10
      endif
      enddo
      dol=1,4950
      linecountr(l)=0
      enddo
      donl=1,nlines
      if(lteflgx(nl).lt.10)then
      iz=znucline(nl)
      i=ionstage(nl)
      linecountr(i+((iz-1)*iz)/2)=linecountr(i+((iz-1)*iz)/2)+1
      endif
      enddo
      callitablsrt(nlines,lteflgx,indexv,1)
      callsgthr(nlines,wavelen,dumvec,indexv)
      callscopy(nlines,dumvec,1,wavelen,1)
      callsgthr(nlines,elower,dumvec,indexv)
      callscopy(nlines,dumvec,1,elower,1)
      calligather(nlines,idumvec,znucline,indexv)
      callicopy(nlines,idumvec,1,znucline,1)
      calligather(nlines,idumvec,ionstage,indexv)
      callicopy(nlines,idumvec,1,ionstage,1)
      calligather(nlines,idumvec,levindex,indexv)
      callicopy(nlines,idumvec,1,levindex,1)
      calligather(nlines,idumvec,lteflgx,indexv)
      callicopy(nlines,idumvec,1,lteflgx,1)
      callsgthr(nlines,gf,dumvec,indexv)
      callscopy(nlines,dumvec,1,gf,1)
      callsgthr(nlines,element,dumvec,indexv)
      callscopy(nlines,dumvec,1,element,1)
      donl=1,nlines
      if(lteflgx(nl).gt.0)goto600
      enddo
600   continue
      nlines=nl-1
      callstblsort(nlines,elower,indexv,1)
      callsgthr(nlines,wavelen,dumvec,indexv)
      callscopy(nlines,dumvec,1,wavelen,1)
      callsgthr(nlines,elower,dumvec,indexv)
      callscopy(nlines,dumvec,1,elower,1)
      calligather(nlines,idumvec,znucline,indexv)
      callicopy(nlines,idumvec,1,znucline,1)
      calligather(nlines,idumvec,ionstage,indexv)
      callicopy(nlines,idumvec,1,ionstage,1)
      callsgthr(nlines,gf,dumvec,indexv)
      callscopy(nlines,dumvec,1,gf,1)
      callsgthr(nlines,element,dumvec,indexv)
      callscopy(nlines,dumvec,1,element,1)
      denergy=enrmax/dble(200-1)
      enrvec(1)=0.d+00
      doi=1,200
      enrvec(i)=dble(i-1)*denergy
      enddo
      enrindex(1)=max(ndexd(dble(elower(1)),enrvec,200),1)
      do650nl=2,nlines
      if(elower(nl).eq.elower(nl-1))then
      enrindex(nl)=enrindex(nl-1)
      else
      enrindex(nl)=max(ndexd(dble(elower(nl)),enrvec,200),1)
      endif
      enrwt(nl)=(elower(nl)-enrvec(enrindex(nl)))/denergy
650   continue
      do420i=1,200
      donr=nr1,nr2
      expvecx(nr+(i-1)*nradii)=dexp(-enrvec(i)*hckt(nr))
      enddo
420   continue
      do660nl=1,nlines
      wlindex(nl)=ndexd(dble(wavelen(nl)),wlgrd,nwave+1)
660   continue
      ndtfac=128
      dtfac=1.d+00/float(ndtfac)
      do920nl=1,nlines
      ne=enrindex(nl)
      nf=wlindex(nl)
      scatx=cross0*(gf(nl)*wavelen(nl))
      iz=znucline(nl)
      io=ionstage(nl)
      do910nr=nr1,nr2
      expterm=enrwt(nl)*expvecx(nr+(ne+1-1)*nradii)+(1.d+00-enrwt(nl))*expvecx(nr+(ne-1)*nradii)
      qfactor=1.d0
      tausob=scatx*expterm*idenz(nr+(io-1)*nradii+(iz-1)*31*nradii)/dvdravg(nr)
      opx=dvdravg(nr)*(1.d+00-(1.d+00/(1.d+00+tausob*dtfac))**ndtfac)/(dlnfreq(nf)*c)
      lscat(nr+(nf-1)*nradii)=lscat(nr+(nf-1)*nradii)+(1.d+00-qfactor)*opx
      lopac(nr+(nf-1)*nradii)=lopac(nr+(nf-1)*nradii)+qfactor*opx
910   continue
920   continue
      donr=nr1,nr2
      donw=1,nwave
      scatop(nw,nr)=scatop(nw,nr)+lscat(nr+(nw-1)*nradii)
      opacity(nw,nr)=opacity(nw,nr)+lopac(nr+(nw-1)*nradii)
      enddo
      enddo
      ptrinit=1
      init=1
      return
      end
      subroutinereadbini(iunit,x,lenx)
      integerx(lenx)
      read(iunit)x
      return
      end
      subroutinereadbinr(iunit,x,lenx)
      real*4x(lenx)
      read(iunit)x
      return
      end
      subroutinewritbini(iunit,x,lenx)
      integerx(lenx)
      write(iunit)x
      return
      end
      subroutinewritbinr(iunit,x,lenx)
      real*4x(lenx)
      write(iunit)x
      return
      end
