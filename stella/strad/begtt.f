      SUBROUTINEBEGINTT
      IMPLICITREAL*8(A-H,O-Z)
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
     *REQ),bolM,NCND,KRAD,NFRUS
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
      COMMON/CONUR/EIT,DST,BBRCNR(5)
      COMMON/BAL/EL(MAXDER+1),YENTOT(MAXDER+1),ETOT0,ELVOL,ELSURF,ELTOT,TPSURF,HOLDBL,ELOST,EKO,RADBEG
      common/NSTEP/NSTEP,NDebug,MAXER,IOUT,NOUT
      common/CREAD/TAUOLD,NSTMAX,MBATCH,MAXORD
      common/debug/LfrDebug,Nperturb,Kbad
      REAL*8TPMAX(MAXDER+1),TQ(4)
      COMMON/TAU/TAU(Mzon+1),FLUX(Mzon)
      common/tauubvri/tauU(Mzon),tauB(Mzon),tauV(Mzon),tauR(Mzon),tauI(Mzon)
      COMMON/PHOT/XJPH,DMPH,RPH,TPH,PLPH,VPH,CHEMPH,GRVPH,HP,JPH
      PARAMETER(NFUNC=6)
      REAL*4WORK(Mzon+2,NFREQ),WRK(Mzon,4)
      REAL*8WRKX(Mzon),WORKX(Mzon+2)
      COMMON/STEPD/WRKX,WORKX,TPHOT,TEFF,WORK,WRK,NPHOT,NZM
      PARAMETER(TMCRIT=1.D-6,TPNSE=5.D0,EPGROW=0.02D0)
      Common/RUTP/Ry(Mzon),Uy(Mzon),Ty(Mzon),Press(Mzon),Rho(Mzon)
      COMMON/TOO/TOO,KO,KNTO,TO(KOMAX),STO(KOMAX),NTO(KOMAX)
      Parameter(Lcurdm=1000)
      RealTcurv
      IntegerNFRUSED
      REAL*8Flsave
      Common/Curve/tcurv(8,Lcurdm),Depos(Lcurdm),Flsave(MFREQ+1,Lcurdm),NFRUSED(Lcurdm),Lsaved
      LOGICALBEGRUN
      Common/BEGR/BEGRUN
      CHARACTER*80Model,Sumprf,Sumcur,Depfile,Flxfile
      COMMON/Files/Model,Sumprf,Sumcur,Depfile,Flxfile
      CHARACTER*1app
      LogicalGivdtl
      Common/ABGrap/NSTA,NSTB,TcurA,TcurB,Givdtl
      REAL*8MBOL,MU,MB,MV,MR,MI,MBOL1
      COMMON/COLOR/MBOL,MU,MB,MV,MR,MI,UMB,BMV,MBOL1,LubvU,LubvB,LubvV,LubvR,LubvI,Lyman
      COMMON/DETAIL/QRTarr(Mzon),UUarr(Mzon),ArrLum(Mzon),Acc(Mzon)
      Common/XYZ/XA,YA,URM
      CHARACTER*80STRING
      Dimensionindfr(6),indop(6)
      CHARACTER*3status,NFILE*80
      Logicalpetread
      Datapetread/.false./
      Dataindfr/1,2,2,3,3,3/,indop/1,1,2,3,3,3/
      DATAWLMAX/5.0D+04/,WLMIN/1.D+00/
      Nc=8
      READ(Nc,'(A)')
      Read(Nc,*)EPS,HMIN,HMAX
      READ(Nc,'(A)')
      Read(Nc,*)METHN,JSTART,MAXORD,KNadap
      READ(Nc,'(A)')
      Read(Nc,*)NSTA,NSTB,TcurA,TcurB
      READ(Nc,'(A)')
      Read(Nc,*)NSTMAX,NDebug,NOUT,IOUT,MBATCH
      Mbatch=Min(Mbatch,Lcurdm)
      READ(Nc,'(A)')
      Read(Nc,*)AMNI,XMNI
      READ(Nc,'(A)')
      Read(Nc,*)AMHT,EBurst,tBurst,tbeght
      READ(Nc,'(A)')
      Read(Nc,*)EKOB,AI1,AI2,AI3,US
      READ(Nc,'(A)')
      Read(Nc,*)THRMAT,CRAP,CONV,EDTM,CHNCND,Givdtl
      READ(Nc,'(A)')
      Read(Nc,*)FLOOR(1),FLOOR(2),FLOOR(3),FLOOR(4)
      READ(Nc,'(A)')
      Read(Nc,*)Wacc(1),Wacc(2),Wacc(3),Wacc(4)
      READ(Nc,'(A)')
      Read(Nc,*)FitTau,TauTol,Rvis,AQ,BQ,DRT,NRT,SCAT
      READ(Nc,'(A)')
      Read(Nc,*)NnTO
      DO09997ito=1,NnTO
      Read(Nc,*)TO(ito)
09997 CONTINUE
      If(LSystem)then
      READ(5,*)IRC
      If(IRC.EQ.0)then
      BEGRUN=.FALSE.
      else
      BEGRUN=.TRUE.
      endif
      endif
      IF(BEGRUN.AND.IABS(NSTB).EQ.1)THEN
      Lunit=9
      NFILE=Model
      write(*,*)'>>',Nfile,'<<'
      callStradIO('rm',Lunit,NFILE)
      write(*,*)' Eko, Ekob: ',Eko,Ekob
      EKO=EKO+EKOB
      If(EKO.GT.0.)then
      i=0
09994 CONTINUE
      i=i+1
      IF(.NOT.((am(i)-amini)/amout.GE.ai1))GOTO09994
      i1=i
09991 CONTINUE
      i=i+1
      IF(.NOT.((am(i)-amini)/amout.GE.ai2))GOTO09991
      i2=i
09988 CONTINUE
      i=i+1
      IF(.NOT.((am(i)-amini)/amout.GE.ai3))GOTO09988
      i3=i
      Z1=US/(AM(I2)-AM(I1))
      IF(I2.NE.I3)Z2=US/(AM(I2)-AM(I3))
      DO09985K=I1,I3
      If(K.LE.I2)then
      Uy(K)=Z1*(AM(K)-AM(I1))
      else
      Uy(K)=Z2*(AM(K)-AM(I3))
      endif
09985 CONTINUE
      Z1=0.D0
      K1=I1+1
      DO09982K=K1,I3
      Z1=Z1+Uy(K)**2*(DM(K)+DM(K+1))
09982 CONTINUE
      DO09979K=I1,I3
      Y(Nzon+K,1)=Uy(K)*2.D0*SQRT((EKO/UEPRI)/Z1)
09979 CONTINUE
      endif
      NCND=NZON
      NFRUS=NFREQ
      KRAD=(NZON-NCND)*NFRUS
      TAUOLD=0.D0
      ELSE
      if(NSTA.LT.0)then
      Lunit=10
      NFILE=Sumprf
      else
      Lunit=12
      NFILE='test.prf'
      endif
      callStradIO('cm',Lunit,NFILE)
      if(NSTA.LT.0)then
      Lunit=11
      NFILE=Sumcur
      else
      Lunit=13
      NFILE='test.crv'
      endif
      callStradIO('rc',Lunit,NFILE)
      ENDIF
      Z1=0.d0
      Z2=0.d0
      DO09976K=1,Nzon
      Ry(K)=Y(K,1)
      Uy(K)=Y(Nzon+K,1)
      If(abs(KNadap).EQ.5)then
      If(K.LT.NZON)Then
      DM2=DM(K+1)
      else
      DM2=DMOUT
      endif
      Z1=Z1+Uy(K)**2*(DM(K)+DM2)
      Z2=Z2+Ry(K)**2*(DM(K)+DM2)
      endif
      Ty(K)=Y(2*Nzon+K,1)
09976 CONTINUE
      If(abs(KNadap).EQ.5)then
      URout=sqrt(Z1/Z2)
      DO09973K=1,Nzon
      Uy(K)=URout*Ry(K)
      Y(Nzon+K,1)=Uy(K)
09973 CONTINUE
      endif
      km=1
09970 IF(.NOT.((AM(km)-AMini)*UM.LE.AMNi*1.0000000001d0.AND.km.LT.nzon))GOTO09969
      km=km+1
      GOTO09970
09969 CONTINUE
      If(km.GE.nzon)then
      write(*,*)' Begrad: AMNI in error! km=',km,'>=nzon =',nzon
      stop
      Endif
      kmnick=km-1
      AMNi=(AM(kmnick)-AMini)*UM
      XNI=XMNI/AMNi
      If(XNI.GE.1.d0.AND.KNadap.GT.0)then
      write(*,*)' Begrad: AMNI too small ! XNI=',XNI
      stop
      Endif
      If(KNadap.GT.0)then
      iferrum=1
09967 IF(.NOT.(iferrum.LE.Natom.AND.Zn(Iferrum).NE.26))GOTO09966
      iferrum=iferrum+1
      GOTO09967
09966 CONTINUE
      If(iferrum.GT.Natom)Then
      write(*,*)' Begrad: Ferrum not found !'
      stop
      Endif
      DO09964km=1,kmnick
      YABUN(iferrum,km)=XNI/AZ(iferrum)
      sum=0.d0
      DO09961j=1,Natom
      if(j.NE.iferrum)sum=sum+YABUN(j,km)*AZ(j)
09961 CONTINUE
      If(sum.LT.1.d-5)Then
      write(*,*)' Begrad: Sum too small=',sum
      stop
      Endif
      DO09958j=1,Natom
      if(j.NE.iferrum)YABUN(j,km)=YABUN(j,km)*(1.d0-XNI)/sum
09958 CONTINUE
09964 CONTINUE
      else
      read(28)XNIfor
      Xmni=0.d0
      DO09955km=1,Nzon
      Xmni=Xmni+Xnifor(km)*dm(km)*UM
09955 CONTINUE
      endif
      IF(CONV)THEN
      DO09952I=1,NZON
      Y(I+(NVARS-2)*NZON,1)=UC(I)
      Y(I+(NVARS-1)*NZON,1)=YAINV(I)
09952 CONTINUE
      ENDIF
      NSTMAX=NSTEP+NSTMAX
      IF(MOD(NSTMAX,MBATCH).NE.0)THEN
      NSTMAX=(NSTMAX/MBATCH+1)*MBATCH
      WRITE(4,'( '' NSTMAX CHANGED:'',I6)')NSTMAX
      ENDIF
      JSTART=0
      UR=10.D0**ULGR
      CLIGHT=CS*UTIME/UR
      IF(CRAOLD.NE.CRAP)THEN
      T=T*(CRAOLD/CRAP)
      H=H*(CRAOLD/CRAP)
      CK1=CK1*(CRAOLD/CRAP)
      CK2=CK2*(CRAOLD/CRAP)
      ENDIF
      UFREQ=BOLTZK*UTP/(2.d0*PI*HPLANC)
      CKRAD=6.D+01/PI**4*CSIGM*UTP**4*UTIME**3/(URHO*UR**3)
      CCL=CS*1.D+08/UFREQ
      CFLUX=60.D0*CSIGM*(UTP/PI)**4
      CLUM=32.D0*PI/3.D0*(CSIGM*UR*UTP**4/URHO)
      CLUMF=4.D0*PI*UR**2*CFLUX
      CIMP=CFLUX*UTIME**2/(CS*UR**2*URHO)
      AS=ACARB
      ZS=ZCARB
      N=NZON*NVARS+2*KRAD
      HMIN=HMIN/UTIME
      HMAX=HMAX/UTIME
      METH=METHN
      Haplim=1.D0/(3.D0*FitTau)
      Uplim=1.D0+Haplim
      WRITE(4,'(''%RUN:'')')
      WRITE(4,'(//30X,''<===== HYDRODYNAMIC RUN OF MODEL '',A,                    ''=====>'',/)')Sumprf
      WRITE(4,'(30X,''MASS(SOLAR)='',F6.3,7X,''RADIUS(SOLAR)='',F9.3/,              30X,''EXPLOSION ENERGY(10**50 ERG)='',1P,E12.5,/
     *)')AMOUT*UM,RADBEG,EKO+Eburst
      WRITE(4,'(30X,''<====='',33X,''=====>'',//)')
      WRITE(4,'(''  INPUT PARAMETERS     '')')
      WRITE(4,'('' EPS   = '',F15.5,9X,'' Rce   = '',1P,E15.5,A)')EPS,RCE*UR/RSOL,' SOL.Rad.'
      WRITE(4,'('' HMIN  = '',1P,E15.5,A,7X,'' AMht  = '',1P,E15.5,A)')HMIN*UTIME,' S',AMht,' SOL.MASS'
      WRITE(4,'('' HMAX  = '',1P,E15.5,A,7X,'' Tburst= '',1P,E15.5,A)')HMAX*UTIME,' S',TBurst,' S'
      WRITE(4,'('' THRMAT= '',1P,E15.5,9X,'' Ebstht= '',1P,E15.5,A)')THRMAT,EBurst,' 1e50 ergs'
      WRITE(4,'('' METH  = '',I15,9X,'' CONV  = '',L15)')METH,CONV
      WRITE(4,'('' JSTART= '',I15,9X,'' EDTM  = '',L15)')JSTART,EDTM
      WRITE(4,'('' MAXORD= '',I15,9X,'' CHNCND= '',L15)')MAXORD,CHNCND
      WRITE(4,'('' NSTA  = '',I15,9X,'' FitTau= '',1P,E15.5)')NSTA,FitTau
      WRITE(4,'('' NSTB  = '',I15,9X,'' TauTol= '',1P,E15.5)')NSTB,TauTol
      WRITE(4,'('' NOUT  = '',I15,9X,'' IOUT  = '',I15)')NOUT,IOUT
      WRITE(4,'('' TcurA = '',F15.5,9X,'' Rvis   ='',F15.5)')TcurA,Rvis
      WRITE(4,'('' TcurB = '',F15.5,9X,'' BQ    = '',1P,E15.5)')TcurB,BQ
      WRITE(4,'('' NSTMAX= '',I15,9X,'' DRT   = '',1P,E15.5)')NSTMAX,DRT
      WRITE(4,'('' XMNI  = '',1P,E15.5,A,'' NRT   = '',G15.3)')XMNI,' SOL.MASS',NRT
      If(KNadap.GT.0)then
      WRITE(4,'('' XNI   = '',1P,E15.5)')XNI
      WRITE(4,'('' CONTM.= '',1P,E15.5,A,'' SCAT  = '',L15)')AMNI,' SOL.MASS',SCAT
      else
      WRITE(4,'('' XNifor= '',1P,E15.5)')XNifor(1)
      WRITE(4,'('' MNicor= '',1P,E15.5,A,'' SCAT  = '',L15)')AMNI,' SOL.MASS',SCAT
      endif
09998 FORMAT(' CK1  =',1P,E12.5,'  CK2=',E12.5,'   CFR=',E12.5,'  CKRAD=',E12.5/' CFLUX=',E12.5,'  CLUM=',E12.5,'  CLUMF=',E12.5)
09999 FORMAT(' UTP=',1P,E12.5,' UTIME=',E12.5,' URHO=',E12.5,' UFREQ=',E12.5)
      WRITE(4,'('' FLOOR :''/1X,1P,10E10.2)')FLOOR
      IF(NSTEP.EQ.0)THEN
      H=1.D+05*HMIN
      HOLDBL=H
      YENTOT(1)=ELOST
      YENTOT(2)=H*ELTOT
      GOTO(09949,09948,09947),indfr(abs(KNadap))
      GOTO09946
09949 CONTINUE
      FREQ(1)=CS*1.D+08/(WLMAX*UFREQ)
      FREQ(NFREQ+1)=CS*1.D+08/(WLMIN*UFREQ)
      Basis=(WLMAX/WLMIN)**(1.D0/(DBLE(NFREQ)))
      DO09945i=2,NFREQ
      FREQ(i)=FREQ(I-1)*BASIS
09945 CONTINUE
      DO09942L=1,NFREQ
      FREQMN(L)=SQRT(FREQ(L)*FREQ(L+1))
09942 CONTINUE
      GOTO09946
09948 CONTINUE
      If(.NOT.petread)then
      read(2)Nfreq0,Msta,Nrho,Ntp,YATab,(Wavel(iif),iif=1,36),TpTab,RhoTab,STab,EpsBand,EppBand
      petread=.true.
      endif
      If(Nfreq.GT.36)then
      DO09939L=36+1,Nfreq
      Wavel(L)=Wavel(L-1)/2.d0
09939 CONTINUE
      endif
      DO09936L=1,NFREQ
      FREQMN(L)=CS*1.D+08/(WaveL(L)*UFREQ)
09936 CONTINUE
      FREQ(1)=0.5d0*Freqmn(1)
      FREQ(Nfreq+1)=2.d0*Freqmn(Nfreq)
      DO09933L=2,NFREQ
      FREQ(L)=0.5d0*(Freqmn(L-1)+Freqmn(L))
09933 CONTINUE
      GOTO09946
09947 CONTINUE
      DO09930i=1,NFREQ+1
      FREQ(i)=dumFREQ(i)/Ufreq
09930 CONTINUE
      DO09927L=1,NFREQ
      FREQMN(L)=dumFREQMN(L)/Ufreq
09927 CONTINUE
09946 CONTINUE
      DO09924L=1,NFREQ
      WEIGHT(L)=(FREQ(L+1)-FREQ(L))*FREQMN(L)**3
      IF(L.LT.NFREQ)THEN
      DLOGNU(L)=1.D0/LOG(FREQMN(L)/FREQMN(L+1))
      ELSE
      DLOGNU(L)=1.D0/LOG(FREQMN(L)/FREQ(L+1))
      ENDIF
09924 CONTINUE
      NFRUS=NFREQ
      DO09921L=1,NFRUS
      NTHICK(L)=ncnd
09921 CONTINUE
      DO09918KM=NCND+1,NZON
      DO09915L=1,NFRUS
      LTHICK(L,KM)=.FALSE.
09915 CONTINUE
09918 CONTINUE
      ENDIF
      doL=1,Nfreq
      freqob(L)=freqmn(L)
      enddo
      If(Mfreq.GT.Nfreq)then
      doL=Nfreq+1,Mfreq
      freqob(L)=freqmn(Nfreq)*exp(-dble(L-Nfreq)/dlognu(1))
      enddo
      endif
      toldg=-1.d0
      WRITE(4,'(''   FREQ:'',1P,10E12.5)')FREQ
      WRITE(4,'(''   FREQMN:'',1P,10E12.5)')FREQMN
      WRITE(4,'(''  WAVE BOUNDS:'',1P,11E10.3)')(CCL/FREQ(LW),LW=1,NFREQ+1)
      WRITE(4,'(''  WAVES:'',1P,10E12.5)')(CCL/FREQMN(LW),LW=1,NFREQ)
      WRITE(4,'('' WEIGHT:'',1P,10E12.5)')WEIGHT
      LubvU=NFRUS
09912 IF(.NOT.(LubvU.GT.1.AND.FREQMN(LubvU).GT.(CCL/3.65D+03)))GOTO09911
      LubvU=LubvU-1
      GOTO09912
09911 CONTINUE
      LubvB=NFRUS
09909 IF(.NOT.(LubvB.GT.1.AND.FREQMN(LubvB).GT.(CCL/4.4D+03)))GOTO09908
      LubvB=LubvB-1
      GOTO09909
09908 CONTINUE
      LubvV=NFRUS
09906 IF(.NOT.(LubvV.GT.1.AND.FREQMN(LubvV).GT.(CCL/5.5D+03)))GOTO09905
      LubvV=LubvV-1
      GOTO09906
09905 CONTINUE
      LubvR=NFRUS
09903 IF(.NOT.(LubvR.GT.1.AND.FREQMN(LubvR).GT.(CCL/7.D+03)))GOTO09902
      LubvR=LubvR-1
      GOTO09903
09902 CONTINUE
      LubvI=NFRUS
09900 IF(.NOT.(LubvI.GT.1.AND.FREQMN(LubvI).GT.(CCL/9.D+03)))GOTO09899
      LubvI=LubvI-1
      GOTO09900
09899 CONTINUE
      Lyman=NFRUS
09897 IF(.NOT.(Lyman.GT.1.AND.FREQMN(Lyman).GT.(CCL/912.D0)))GOTO09896
      Lyman=Lyman-1
      GOTO09897
09896 CONTINUE
      WRITE(4,*)' L UBVRI Lyman:',LubvU,LubvB,LubvV,LubvR,LubvI,Lyman
      LST=2
      KENTR=0
      RHO(1)=3.D0*DM(1)/(Ry(1)**3-Rce**3)
      DO09894I=2,NZON
      RHO(I)=3.D0*DM(I)/(Ry(I)**3-Ry(I-1)**3)
09894 CONTINUE
      IF(NSTEP.EQ.0)THEN
      CALLLOSSEN
      HOLDBL=H
      YENTOT(1)=ELOST
      YENTOT(2)=H*ELTOT
      ENDIF
      chem=0.5d0
      RETURN
9     WRITE(4,*)' BEGIN: error in read file ',Lunit
      stop
      END
