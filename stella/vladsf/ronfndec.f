      Programronfict
      implicitreal*8(a-h,o-z)
      character*132tablestr
      PARAMETER(NVARS=3)
      include '../obj/nfreq_and_mzone.inc'
      PARAMETER(NYDIM=(NVARS+2*NFREQ)*Mzon,MAXDER=4)
      Parameter(Is=5)
      PARAMETER(NZ=3000000)
      Parameter(Nstage=28,Natom=15)
      PARAMETER(KOMAX=80)
      LogicalLSYSTEM
      Parameter(LSystem=.FALSE.)
      Parameter(Pi=3.1415926535897932d+00,hPlanc=1.0545716280D-27,Cs=2.9979245800D+10,Boltzk=1.3806504000D-16,Avogar=6.0221417900D+2
     *3,AMbrun=1.6605387832D-24,AMelec=9.1093821500D-28,echarg=4.8032042700D-10,CG=6.6742800000D-08,CMS=1.9884000000D+33,RSol=6.9551
     *000000D+10,ULGR=1.4000000000D+01,UPURS=1.0000000000D+07,ULGPU=7.0000000000D+00,ULGEU=1.3000000000D+01,UPC=3.0856776000D+18,UTP
     *=1.0000000000D+05,URHO=1.0000000000D-06,CARAD=7.5657680191D-15,CSIGM=5.6704004778D-05,ERGEV=1.6021764864D-12,GRADeV=1.16045052
     *85D+04,RADC=7.5657680191D-02,CTOMP=4.0062048575D-01,CCAPS=2.6901213726D+01,CCAPZ=9.8964034725D+00)
      IntegerZn(Natom),ZnCo(Natom+1)
      DimensionAZ(Natom)
      Common/AZZn/AZ,Zn,ZnCo
      Common/NiAdap/tday,t_eve,XNifor(Mzon),AMeveNi,KNadap
      LOGICALFRST
      Parameter(Mfreq=130)
      Common/Kmzon/km,kmhap,Jac,FRST
      COMMON/STCOM1/t,H,HMIN,HMAX,EPS,N,METH,KFLAG,JSTART
      COMMON/YMAX/YMAX(NYDIM)
      COMMON/YSTIF/Y(NYDIM,MAXDER+1)
      COMMON/HNUSED/HUSED,NQUSED,NFUN,NJAC,NITER,NFAIL
      COMMON/HNT/HNT(7)
      PARAMETER(DELTA=1.d-05)
      PARAMETER(LICN=4*NZ,LIRN=2*NZ)
      LogicalNEEDBR
      COMMON/STJAC/THRMAT,HL,AJAC(LICN),IRN(LIRN),ICN(LICN),WJAC(NYDIM),FSAVE(NYDIM*2),IKEEP(5*NYDIM),IW(8*NYDIM),IDISP(11),NZMOD,NE
     *EDBR
      LOGICALCONV,CHNCND,SCAT,SEP
      COMMON/CUTOFF/FLOOR(NVARS+1),Wacc(NVARS+1),FitTau,TauTol,Rvis,CONV,CHNCND,SCAT,SEP
      LogicalLTHICK
      COMMON/THICK/LTHICK(Nfreq,Mzon)
      COMMON/CONVEC/UC(Mzon),YAINV(Mzon)
      COMMON/RAD/EDDJ(Mzon,Nfreq),EDDH(Mzon),HEDD(Nfreq),HEDRAD,CLIGHT,CKRAD,UFREQ,CFLUX,CCL,CLUM,CLUMF,CIMP,NTHICK(NFREQ),NTHNEW(NF
     *REQ),NCND,KRAD,NFRUS
      LOGICALEDTM
      COMMON/RADOLD/HEDOLD,HINEDD,EDTM
      Common/newedd/EddN(Mzon,Nfreq),HEdN(Nfreq),tfeau
      Common/oldedd/EddO(Mzon,Nfreq),HEdo(Nfreq),trlx
      Common/cnlast/Cnlast
      Common/Dhap/DHaphR(Mzon,Nfreq)
      COMMON/BAND/FREQ(NFREQ+1),FREQMN(NFREQ),WEIGHT(130),HAPPAL(NFREQ),HAPABSRON(NFREQ),HAPABS(NFREQ),DLOGNU(NFREQ)
      PARAMETER(NFRMIN=Nfreq/2)
      IntegerdLfrMax
      Common/observer/wH(Mfreq),cH(Mfreq),zerfr,Hcom(Mfreq),Hobs(Mfreq),freqob(Mfreq),dLfrMax
      Parameter(NP=15+15-1)
      Common/famu/fstatic(0:NP+1,Nfreq),fobs_corr(0:NP+1,Mfreq),fcom(0:NP+1,Mfreq),amustatic(0:NP+1)
      Common/rays/Pray(0:Np+1),fout(0:NP+1,Mfreq),abermu(0:NP+1),NmuNzon
      COMMON/LIM/Uplim,Haplim
      COMMON/AMM/DMIN,DM(Mzon),DMOUT,AMINI,AM(Mzon),AMOUT
      COMMON/Centr/RCE,Nzon
      Common/InEn/AMHT,EBurst,tBurst,tbeght
      COMMON/RADPUM/AMNI,XMNi,XNi,KmNick
      COMMON/RADGAM/FJgam(Mzon,2),toldg,tnewg
      COMMON/RADGlg/FJglog(Mzon,2)
      COMMON/CHEM/CHEM0(Mzon),RTphi(0:Mzon+1),EpsUq
      COMMON/REGIME/NREG(Mzon)
      doubleprecisionNRT
      COMMON/AQ/AQ,BQ,DRT,NRT
      COMMON/AZNUC/ACARB,ZCARB,ASI,ZSI,ANI,ZNI,QCSI,QSINI
      COMMON/QNRGYE/QNUC,RGASA,YELECT
      COMMON/CKN1/CK1,CK2,CFR,CRAP,CRAOLD
      LOGICALEVALJA,OLDJAC,BADSTE
      COMMON/EVAL/EVALJA,BADSTE,OLDJAC
      LogicalRadP
      COMMON/RadP/RadP
      COMMON/ARG/TP,PL,CHEM,LST,KENTR,JURS
      COMMON/RESULT/P,Egas,Sgas,ENG,CAPPA,PT,ET,ST,ENGT,CAPT,NZR
      COMMON/ABUND/XYZA,Yat(Natom)
      COMMON/AZ/AS,ZS,SCN
      COMMON/STR/PPL,EPL,SPL,ENGPL,CAPPL,CP,GAM,DA,DPE,DSE,DSP,BETgas
      COMMON/XELECT/XE,XET,XEPL,PE,Ycomp
      COMMON/URScap/Tpsqrt,Psicap,Scap,ScapT,ScapPl,ZMean,YZMean,ZMT,ZMPl,YZMT,YZMPl
      COMMON/BURNCC/CC,CCTP,CCPL,YDOT
      COMMON/ABarr/YABUN(Natom,Mzon)
      COMMON/UNSTL/UL,UPRESS,UE,UEPS,UCAP,UTIME,UV,UFLUX,UP
      COMMON/TAIL/KTAIL
      COMMON/UNINV/UPI,UEI
      COMMON/UNBSTL/UR,UM,UEPRI,ULGP,ULGE,ULGV,ULGTM,ULGEST,ULGFL,ULGCAP,ULGEPS
      Parameter(epsNi=3.9656224881D-02,epsCo=7.0232694190D-03,epspst=2.3346047930D-04,tNi=8.8000000000D+00,tCo=1.1130000000D+02)
      Character*80Opafile
      real*4hptabab,hptababron,hptabsc,hpsavsc
      Common/Opsave/hpsavsc(Nfreq,14,14,Mzon/(Mzon/50)+1,6)
      Common/OpBand/TpTab(14),RhoTab(14),STab(6),Wavel(Nfreq),YATab(Natom),hptabab(Nfreq,14,14,Mzon/(Mzon/50)+1,2),hptababron(Nfreq,
     *14,14,Mzon/(Mzon/50)+1,2),hptabsc(Nfreq,14,14,Mzon/(Mzon/50)+1,2),Msta,Nrho,Ntp,im
      real*4hpbanab1(Nfreq,14,14,Mzon/(Mzon/50)+1),hpbanabron1(Nfreq,14,14,Mzon/(Mzon/50)+1),hpbansc1(Nfreq,14,14,Mzon/(Mzon/50)+1),
     *hpbanab2(Nfreq,14,14,Mzon/(Mzon/50)+1),hpbanabron2(Nfreq,14,14,Mzon/(Mzon/50)+1),hpbansc2(Nfreq,14,14,Mzon/(Mzon/50)+1)
      Equivalence(hptabab(1,1,1,1,1),hpbanab1(1,1,1,1)),(hptababron(1,1,1,1,1),hpbanabron1(1,1,1,1)),(hptabsc(1,1,1,1,1),hpbansc1(1,
     *1,1,1)),(hptabab(1,1,1,1,2),hpbanab2(1,1,1,1)),(hptababron(1,1,1,1,2),hpbanabron2(1,1,1,1)),(hptabsc(1,1,1,1,2),hpbansc2(1,1,1
     *,1))
      Common/tintrp/stmlog(6),tdlog,thaplog1,thaplog2,istold,Opafile
      Common/dumfreq/dumFreq(Nfreq+1),dumFreqmn(Nfreq),dumwavel(Nfreq)
      Parameter(NatomMax=99)
      logicalfexist
      Character*80Nfile,model,NiDist,Rfile
      Character*1app
      datadaybegl/0.d0/,dayendl/2.d0/
      PARAMETER(MAXIONM1=(6-1))
      PARAMETER(NSBINTVL=3*10)
      real*8scatop(Nfreq),opac(Nfreq)
      real*8eden,temp,nucden
      real*8abund(99)
      character*80linelist,longlist,xsecdatadir
      logicalrdlndump
      logicaltstwritten
      callGetArg(1,tablestr)
      read(tablestr,*)itable
      write(*,*)'itable',itable
      if(itable.LT.1.or.itable.GT.6)stop'bad value for itable: must be 1..6'
      Open(unit=1,file='ronfict.1',status='old',action='read')
      READ(1,'(A)')
      read(1,*)Model,NiDist,Rfile,Opafile
      Open(unit=9,file=Model,form='unformatted',action='read')
      Open(unit=10,file=NiDist,form='unformatted',action='read')
      app=char(ichar('0')+itable)
      Rfile=Rfile(1:Length(Rfile)-5)//app//Rfile(Length(Rfile)-3:Length(Rfile))
      inquire(file=Rfile,exist=fexist)
      nures=4
      write(*,*)'result file '//trim(Rfile)
      if(fexist)then
      open(file=Rfile,unit=nures,status='old')
      else
      open(file=Rfile,unit=nures,status='new')
      endif
      inquire(file='lines'//app//'.out',exist=fexist)
      if(fexist)then
      open(file='lines'//app//'.out ',unit=7,status='old')
      else
      open(file='lines'//app//'.out ',unit=7,status='new')
      endif
      open(15,file='testabun'//app//'.res')
      DO09999nf=itable,itable
      nu=70+nf
      app=char(ichar('0')+nf)
      Opafile=Opafile(1:Length(Opafile)-1)//app
      Open(unit=nu,file=Opafile,form='unformatted')
      write(*,*)'Opafile '//trim(Opafile)
      if(nf.EQ.6)Open(unit=69,file=Opafile(1:Length(Opafile)-1)//'ab',form='unformatted')
09999 CONTINUE
      INQUIRE(FILE='linedata.dump',EXIST=rdlndump)
      if(rdlndump)then
      linelist='linedata.dump'
      else
      linelist='lineatom.dat '
      endif
      xsecdatadir='./yakovlev '
      Lunit=9
      NFILE=Model
      callStradIO('rm',Lunit,NFILE)
      read(10)(XNifor(k),k=1,NZon)
      wlmax=5.d+04
      wlmin=1.d0
      freq(1)=Cs*1.d+08/wlmax
      freq(nfreq+1)=Cs*1.d+08/wlmin
      Basis=(WLMAX/WLMIN)**(1.D0/(DBLE(NFREQ)))
      DO09996i=2,NFREQ
      FREQ(i)=FREQ(I-1)*BASIS
09996 CONTINUE
      DO09993L=1,NFREQ
      FREQMN(L)=SQRT(FREQ(L)*FREQ(L+1))
09993 CONTINUE
      dolf=1,nfreq
      wavel(lf)=Cs*1.d+08/freqmn(lf)
      enddo
      tstwritten=.false.
      DO09990nw=itable,itable
      nu=70+nw
      time=1.d+01**(daybegl+(dayendl-daybegl)/dble(6-1)*(nw-1))
      dvdr=1.0d0/(time*86400.d0)
      dt=time-t*Utime/86400.d0
      write(*,*)' time in tables and T_eve:',time,t*Utime/86400.d0
      if(dt.LT.0.d0)then
      write(*,*)' time in tables is less than T_eve:',time,t*Utime/86400.d0
      stop
      endif
      DO09987km=1,Nzon,(Mzon/50)
      doiz=1,99
      abund(iz)=0.d0
      enddo
      DO09984ii=1,Natom-2
      abund(Zn(ii))=Yabun(ii,km)*Az(ii)
09984 CONTINUE
      write(nures,'(a,1p,e10.3)')' Fe abund in eve model=',Yabun(Natom-1,km)*AZ(Natom-1)
      write(nures,'(a,1p,e10.3)')' Ni stable abund in eve model=',Yabun(Natom,km)*AZ(Natom)
      write(nures,'(a,1p,e10.3)')' Ni56 abund in eve model=',Xnifor(Km)
      abund(28)=Xnifor(Km)*exp(-dt/tNi)
      abund(27)=Xnifor(Km)*(tCo/(tCo-tNi))*(exp(-dt/tCo)-exp(-dt/tNi))
      abund(26)=Xnifor(Km)-abund(28)-abund(27)
      write(nures,'(a,1p,3e10.3)')' Fe Co Ni from Ni56=',abund(26),abund(27),abund(28)
      abund(28)=abund(28)+Yabun(Natom,km)*AZ(Natom)
      abund(26)=Yabun(Natom-1,km)*AZ(Natom-1)-Xnifor(Km)+abund(26)
      write(nures,'(a,1p,3e10.3)')' Fe Co Ni abund here=',abund(26),abund(27),abund(28)
      if(.NOT.tstwritten)then
      write(15,'(i5,99e10.3)')km,abund
      endif
      lastTp=1
      DO09981ir=1,14
      rho=1.d+01**(-18+ir)
      rhotab(ir)=log(rho)
      if(ir.EQ.1)then
      rho=1.d+01*rho
      elseif(ir.EQ.14)then
      rho=rho/1.d+01
      endif
      nucden=abund(1)*rho/Ambrun
      doiz=2,99
      nucden=nucden+abund(iz)*rho/(dble(2*iz)*Ambrun)
      enddo
      If(lastTp.EQ.1)then
      ktpb=1
      lastTp=14
      kstep=1
      else
      ktpb=14
      lastTp=1
      kstep=-1
      endif
      DO09978itp=ktpb,lastTp,kstep
      temp=1.d+01**(3.4d0+0.2d0*(itp-1))
      Tptab(itp)=log(temp)
      if(itp.EQ.1)then
      temp=1.d+01**(3.4d0+0.2d0*itp)
      elseif(itp.EQ.14)then
      temp=1.d+01**(3.4d0+0.2d0*(14-2))
      endif
      callopacity(rho,temp,abund,nfreq,freq,freqmn,dvdr,linelist,rdlndump,.false.,longlist,0.0,xsecdatadir,opac,scatop,eden)
15    format(' rho: ',1pe10.3,' temp: ',1pe10.3,' time(days): ',1pg12.3,/,' eden: ',1pe10.3,' eden/nucden: ',1pe10.3,/)
      dolf=1,nfreq
      if(nw.EQ.6)hpbanab1(lf,itp,ir,(km-1)/(Mzon/50)+1)=alog(sngl(opac(lf)*1.d+01**UlgR/(rho/Urho)))
      hpbansc1(lf,itp,ir,(km-1)/(Mzon/50)+1)=alog(sngl((scatop(lf)+opac(lf))*1.d+01**UlgR/(rho/Urho)))
      hbab=hpbanab1(lf,itp,ir,(km-1)/(Mzon/50)+1)
      hbal=hpbansc1(lf,itp,ir,(km-1)/(Mzon/50)+1)
      if(hbal.GT.1.d50.or.hbab.GT.1.d50)then
      write(nures,'(4(a,i4)/,1p,2(a,e12.4))')' L=',lf,' itp=',itp,' iro=',ir,' im=',(km-1)/(Mzon/50)+1,'  hbal=',hbal,'  hbab=',hbab
      write(nures,'(1p,2(a,e12.4))')'  opac=',opac(lf),'  scatop=',scatop(lf)
      stop' Ronfict: bad happa !!!'
      endif
      enddo
09978 CONTINUE
      write(nures,15)rho,temp,time,eden,eden/nucden
09981 CONTINUE
      write(nures,'(A,F8.2,A,I4)')' done time=',time,'  zone=',km
09987 CONTINUE
      tstwritten=.true.
      Msta=6
      Nrho=14
      NTp=14
      Stime=time
      write(nu)nw,Stime,Nfreq,Msta,Nrho,NTp,Wavel,Freq,Freqmn,TpTab,RhoTab,hpbansc1
      close(unit=nu)
      if(nw.EQ.6)then
      write(69)hpbanab1
      endif
09990 CONTINUE
      close(nures)
      close(9)
      close(10)
      close(15)
      end
