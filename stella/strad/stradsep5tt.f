      PROGRAMSTELLA
      IMPLICITREAL*8(A-H,O-Z)
      REALLTIM,MTIM
      REALB0001,TIMEND
      CHARACTER*10D0001
      CHARACTER*8CT001
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
      Common/Volanm/Ryzer,Kmcor
      PARAMETER(NzMIN=Mzon*5/6)
      PARAMETER(NCNMIN=5)
      Parameter(TREZON=3.D18)
      Parameter(tsmoth=4.D6)
      Parameter(tret=1.d0/3.d0)
      LogicalHOLDFR
      LogicalLEXIST1,LEXIST2
      Character*80Dfile,Rfile,bmfile,Nidist,swdet,balafile
      character*160wordm(10),runname,filestr
      integerlwordm(10),lrun,nwordm,errnum
      Character*8STRBAT
      DimensionYsave(NYdim)
      saveYsave
      Common/FRAD/FRADJ(Mzon,Nfreq),FRADH(Mzon,Nfreq)
      Common/out/tretard,tout,Yout(NYdim)
      BLACK(Lbl,Tpbl)=(exp(-(FREQMN(Lbl)/Tpbl)))/(1.d0-(exp(-(FREQMN(Lbl)/Tpbl))))
      BLACKD(Lbl,Tpbl)=(FREQMN(Lbl)/Tpbl)*(exp(-(FREQMN(Lbl)/Tpbl)))/(1.d0-(exp(-(FREQMN(Lbl)/Tpbl))))**2
      Open(unit=1,file='strad.1',status='old')
      READ(1,'(A)')
      read(1,'(A)')filestr
      close(1)
      callwords(filestr,wordm,lwordm,nwordm)
      if(nwordm.LT.5)then
      write(*,*)' check your start file, something is missing !'
      stop12
      endif
      lrun=lwordm(1)
      runname(1:)=wordm(1)(1:lrun)
      rfile(1:lwordm(2))=wordm(2)(1:lwordm(2))
      Narg=Iargc()
      GOTO(09999,09998),Narg+1
      write(*,*)' Extra arguments.'
      Stop16
      GOTO09997
09999 CONTINUE
      GOTO09997
09998 CONTINUE
      callGetArg(1,wordm(1))
      runname=wordm(1)(1:len_trim(wordm(1)))
      rfile=runname(1:len_trim(runname))//'.res'
      lrun=len_trim(runname)
      lwordm(3)=8
      wordm(3)=runname(1:lwordm(3))//'_2000h'
      lwordm(3)=lwordm(3)+6
      lwordm(4)=lwordm(3)
      wordm(4)=wordm(3)
09997 CONTINUE
      dfile(1:)=runname(1:lrun)//'.dat'
      bmfile(1:)=runname(1:lrun)//'.bm'
      model(1:)='../../eve/run/'//wordm(3)(1:lwordm(3))//'.mod'
      Nidist(1:)='../../eve/run/'//wordm(4)(1:lwordm(4))//'.xni'
      Opafile(1:)='../../vladsf/'//wordm(5)(1:lwordm(5))
      Sumprf(1:)=runname(1:lrun)//'.prf'
      SumCur(1:)=runname(1:lrun)//'.crv'
      Depfile(1:)=runname(1:lrun)//'.dep'
      Flxfile(1:)=runname(1:lrun)//'.flx'
      swdet(1:)=runname(1:lrun)//'.swd'
      balafile(1:)=runname(1:lrun)//'.balance'
      INQUIRE(FILE=Sumprf,EXIST=LEXIST1)
      INQUIRE(FILE=sumcur,EXIST=LEXIST2)
      IF(LEXIST1.AND.LEXIST2)THEN
      BEGRUN=.FALSE.
      ELSEIF(.NOT.LEXIST1.AND..NOT.LEXIST2)THEN
      BEGRUN=.TRUE.
      ELSE
      Write(*,'(A)')' File Sumprf or SUMCUR is missing'
      stop28
      ENDIF
      Open(unit=8,file=Dfile,status='old')
      Open(unit=4,file=Rfile,status='new')
      Open(unit=60,file=bmfile,status='new')
      Open(unit=81,file=swdet,status='new')
      Open(unit=28,file=Nidist,form='unformatted',status='unknown')
      Open(unit=76,file=balafile,status='unknown')
      callVtime(B0001)
      If(Lsystem)then
      else
      CallTremain(LTIM,B0001)
      endif
      EIT=1.D-5
      DST=1.D-4
      CALLBEGIN
      close(8)
      Close(unit=28)
      FLnois=Floor(NVARS+1)
      If(Lsystem)CALLBASTAT('BEGIN RD')
      NEEDBR=.TRUE.
      NSTEP0=NSTEP
      doK=1,Nzon
      Ry(K)=Y(K,1)
      Uy(K)=Y(Nzon+K,1)
      Ty(K)=Y(2*Nzon+K,1)
      enddo
      Rho(1)=3.D0*DM(1)/(Ry(1)**3-Rce**3)
      doI=2,Nzon
      Rho(I)=3.D0*DM(I)/(Ry(I)**3-Ry(I-1)**3)
      enddo
      CallLOSSEN
      tday=t*UTIME/8.64d+04
      if(NSTEP.EQ.0)t_eve=tday
      I=NZON
      TAU(I)=0.D0
09996 IF(.NOT.(I.GT.1))GOTO09995
      TP=Ty(I)
      Pl=3.D0*DM(I)/(Ry(I)**3-Ry(I-1)**3)
      DO09993ii=1,Natom
      Yat(ii)=YABUN(ii,I)
09993 CONTINUE
      RADP=.FALSE.
      CALLURSOS
      Rho(I)=Pl
      Press(I)=P*UPI
      kmhap=I
      CALLHAPPA
      Hapmin=Happal(1)
      Ltau=1
      DO09990L=2,Nfrus
      If(Happal(L).LT.Hapmin)then
      Hapmin=Happal(L)
      Ltau=L
      endif
09990 CONTINUE
      TAUste=(Ry(I)-Ry(I-1))*HAPPAL(LTAU)*PL
      If(Scat)then
      TAUabs=(Ry(I)-Ry(I-1))*HAPabs(LTAU)*PL
      TAUsca=Tauste-Tauabs
      If(Tausca.GT.1.d0)then
      Tauste=sqrt(Tauabs*Tausca)
      else
      Tauste=Tauabs
      endif
      endif
      TAU(I-1)=TAU(I)+Tauste
      I=I-1
      GOTO09996
09995 CONTINUE
      Ncnew=Nthick(Ltau)
      DO09987K=Ncnew-1,1,-1
      tau(K)=tau(K)-tau(Ncnew)
09987 CONTINUE
09984 IF(.NOT.(TAU(max(Ncnew,1)).LE.FitTau.AND.Ncnew.GE.NCNMIN))GOTO09983
      Ncnew=Ncnew-1
      GOTO09984
09983 CONTINUE
      IF(Ncnew.LT.NCNMIN)Ncnew=0
      Ncnew=MIN(Ncnew,NZON-3)
      HOLDFR=.TRUE.
      IF(.NOT.((TAU(max(Ncnd,1)).GT.TauTol*FitTau.or.TAU(max(Ncnd,1))*TauTol.LT.FitTau).AND.Ncnew.NE.Ncnd.AND.(CHNCND.or.Nstep.EQ.0)
     *))GOTO09981
      IF(.NOT.(Ncnew.LT.Ncnd))GOTO09978
      DO09974L=1,Nfrus
      DO09971IK=Ncnd+1,NZON
      EddO(IK,L)=EddJ(IK,L)
09971 CONTINUE
      HEdO(L)=HEdd(L)
09974 CONTINUE
      CALLFeau(Ncnew,Ncnd)
      LFR=NFRUS
      Tpcomp=Ty(NZON)
      DO09968I=Ncnew+1,NZON
      If(Ty(I).GT.Tpcomp)Tpcomp=Ty(I)
09968 CONTINUE
09965 IF(.NOT.(BLACK(MIN(LFR+1,NFREQ),Tpcomp)*FREQMN(MIN(LFR+1,NFREQ))**3.GT.FLNOIS))GOTO09964
      LFR=LFR+1
      IF(.NOT.(LFR.GE.NFREQ))GOTO09965
09964 CONTINUE
      LFR=min(LFR,NFREQ)
      IF(LFR.NE.NFRUS)THEN
      WRITE(4,*)'Tpcomp=',Tpcomp
      WRITE(4,*)' NEW NFRUS -->',LFR,' AT STEP=',NSTEP
      ENDIF
      KRAD1=(NZON-Ncnew)*LFR
      DO09962L=NFRUS,1,-1
      IX=NVARS*NZON+(NZON-Ncnd)*(L-1)
      DO09959K=NZON-1,Ncnd+1,-1
      Y(IX+(K-Ncnd)+(Ncnd-Ncnew)*L+KRAD1,1)=Y(IX+(K-Ncnd)+KRAD,1)
09959 CONTINUE
09962 CONTINUE
      DO09956L=NFRUS,1,-1
      IX=NVARS*NZON+(NZON-Ncnd)*(L-1)
      DO09953K=NZON,Ncnd+1,-1
      Y(IX+(K-Ncnd)+(Ncnd-Ncnew)*L,1)=Y(IX+(K-Ncnd),1)
09953 CONTINUE
09956 CONTINUE
      DO09950K=Ncnd+1,NZON
      DO09947L=NFRUS+1,LFR
      IX=NVARS*NZON+(NZON-Ncnew)*(L-1)
      Y(IX+(K-Ncnew),1)=0.D0
      Y(IX+(K-Ncnew)+KRAD1,1)=0.D0
09947 CONTINUE
09950 CONTINUE
      NFRUSO=NFRUS
      NFRUS=LFR
      CallNTH(Ncnew)
      DO09944L=NFRUSO+1,NFRUS
      DO09941IK=Ncnd+1,NZON
      EddO(IK,L)=Eddn(IK,L)
      EDDJ(IK,L)=Eddn(IK,L)
09941 CONTINUE
      HEdO(L)=HEDD(L)
09944 CONTINUE
      trlx=25.d0*h
      CALLGDEPOS
      KTAIL=NSTEP
      Tp=Ty(Ncnew+1)
      If(Ncnew.GT.0)Then
      PL=3.D0*DM(Ncnew+1)/(Ry(Ncnew+1)**3-Ry(Ncnew)**3)
      Else
      PL=3.D0*DM(Ncnew+1)/(Ry(Ncnew+1)**3-Rce**3)
      Endif
      RADP=.FALSE.
      DO09938ii=1,Natom
      Yat(ii)=YABUN(ii,Ncnew+1)
09938 CONTINUE
      CALLURSOS
      kmhap=Ncnew+1
      CALLOPACIT
      CAP2=CAPPA
      DO09935K=Ncnew+1,Ncnd
      If(K.LT.Nzon)Then
      TP1=Tp
      CAP1=CAP2
      Tp=Ty(K+1)
      PL=3.D0*DM(K+1)/(Ry(K+1)**3-Ry(K)**3)
      DO09932ii=1,Natom
      Yat(ii)=YABUN(ii,K+1)
09932 CONTINUE
      CALLURSOS
      kmhap=K+1
      CALLOPACIT
      callHAPPA
      CAP2=CAPPA
      FLLF1=((Ry(K)*TP1)**4-(Ry(K)*TP)**4)*TP1**4/(CAP1*(DM(K)+DM(K+1))*(TP1**4+TP**4))
      FLRT1=((Ry(K)*TP1)**4-(Ry(K)*TP)**4)*TP**4/(CAP2*(DM(K)+DM(K+1))*(TP1**4+TP**4))
      FLUM=CLUM*(FLLF1+FLRT1)
      Endif
      If(NSTEP.EQ.0.AND.EddJ(K,LubvB).GT.0.25D0)Then
      dillm=0.75d0-0.25d0*sqrt(12.d0*EddJ(K,LubvB)-3.D0)
      Tpm=Freqmn(LubvB)/LOG(1.D0+dillm/FradJ(K,LubvB))
      else
      dillm=0.75d0
      TPM=MAX(Ty(K),2.D+3/UTP)
      endif
      FLUN=0.D0
      DO09929L=1,NFRUS
      IX=NVARS*NZON+(NZON-Ncnew)*(L-1)
      If(K.LE.Nthick(L))then
      Y(IX+(K-Ncnew),1)=Black(L,Ty(k))
      else
      Y(IX+(K-Ncnew),1)=FradJ(K,L)
      endif
      If(K.LT.Nzon)Then
      If(Nstep.GT.1.AND.K.GT.Nthick(L))then
      Y(IX+(K-Ncnew)+KRAD1,1)=FradH(K,L)
      elseif(Nstep.GT.1)then
      Y(IX+(K-Ncnew)+KRAD1,1)=(FLLF1+FLRT1)*Blackd(L,0.5D0*(Tp+Tp1))*CAP2/(6.D0*(0.5D0*(Tp+Tp1))**4*Ry(K)**2*HAPPAL(L))
      else
      Y(IX+(K-Ncnew)+KRAD1,1)=0.1d0*FradH(K,L)
      endif
      FLUN=FLUN+Y(IX+(K-Ncnew)+KRAD1,1)*WEIGHT(L)
      Endif
09929 CONTINUE
      If(K.LT.Nzon)Then
      FLUN=CLUMF*FLUN*Ry(K)**2
      WRITE(4,'(2(A,1P,E12.5),A,I3)')' Cond Lum=',Flum,'      New  Lum=',Flun,'   Km=',k
      Endif
09935 CONTINUE
      Ncnd=Ncnew
      KRAD=KRAD1
      N=NZON*NVARS+2*KRAD
      WRITE(4,*)' NEW Ncnd KRAD N:   ',Ncnd,KRAD,N,' AT STEP=',NSTEP
      HOLDFR=.FALSE.
      GOTO09975
09978 CONTINUE
      IF(.NOT.(Ncnew.GT.Ncnd.AND.Ncnew.LE.NZON-3))GOTO09977
      DO09926L=1,Nfrus
      DO09923IK=Ncnd+1,NZON
      EddO(IK,L)=EddJ(IK,L)
09923 CONTINUE
      HEdO(L)=HEdd(L)
09926 CONTINUE
      CALLFeau(Ncnew,Ncnd)
      KRAD1=(NZON-Ncnew)*NFRUS
      DO09920L=1,NFRUS
      IX=NVARS*NZON+(NZON-Ncnd)*(L-1)
      DO09917K=Ncnew+1,NZON
      Y(NVARS*NZON+(NZON-Ncnew)*(L-1)+(K-Ncnew),1)=Y(IX+(K-Ncnd),1)
09917 CONTINUE
09920 CONTINUE
      DO09914L=1,NFRUS
      IX=NVARS*NZON+(NZON-Ncnd)*(L-1)
      DO09911K=Ncnew+1,NZON-1
      Y(NVARS*NZON+(NZON-Ncnew)*(L-1)+(K-Ncnew)+KRAD1,1)=Y(IX+(K-Ncnd)+KRAD,1)
09911 CONTINUE
09914 CONTINUE
      trlx=25.d0*h
      CALLGDEPOS
      KTAIL=NSTEP
      Ncnd=Ncnew
      KRAD=KRAD1
      N=NZON*NVARS+2*KRAD
      WRITE(4,*)' NEW Ncnd KRAD N:   ',Ncnd,KRAD,N,' AT STEP=',NSTEP
      WRITE(4,*)'Ty(Ncnd+1)=',Ty(Ncnd+1)
      GOTO09975
09977 CONTINUE
      IF(.NOT.(Ncnew.GT.NZON-3))GOTO09976
      KFLAG=-3
      GOTO09975
09976 CONTINUE
09975 CONTINUE
      If(KFLAG.EQ.0)Then
      JSTART=0
      NEEDBR=.TRUE.
      endif
09981 CONTINUE
      RTphi(0)=0.
      RTphi(Nzon+1)=0.
      Uscale=1.d9*Utime/UR
      Ryzer=Ry(Kmcor)
      RTphi(Nzon+1)=0.d0
      EpsUq=0.d0
      DO09908K=1,NZON
      EpsUq=max(EpsUq,DRT*Eps*abs(Uy(K)))
      TauRT=(Amout-Am(K))*(0.4d0*Urho*UR)/Ry(K)**2
      RTphi(K)=BQ*dM(K)*Ry(K)**(NRT-1.d0)/(1.d0+0.033d0*TauRT)
09908 CONTINUE
      Tpcomp=Ty(NZON)
      DO09905I=Ncnd+1,NZON
      If(Ty(I).GT.Tpcomp)Tpcomp=Ty(I)
09905 CONTINUE
      Ycomp=(Boltzk*UTp*Tpcomp)/(Amelec*Cs**2)
      LFR=NFRUS
09902 IF(.NOT.(BLACK(LFR,max(Tpcomp,floor(3)))*FREQMN(LFR)**3.LT.FLNOIS/1.6D0.AND..NOT.HOLDFR.AND.Lfr.GT.NFRMIN))GOTO09901
      LFR=LFR-1
      GOTO09902
09901 CONTINUE
      IF(LFR.LT.NFRMIN)LFR=NFRMIN
      IF(LFR.NE.NFRUS)THEN
      NFRUS=LFR
      KRAD1=(NZON-Ncnd)*NFRUS
      DO09899L=1,NFRUS
      IX=NVARS*NZON+(NZON-Ncnd)*(L-1)
      DO09896K=Ncnd+1,NZON-1
      Y(IX+(K-Ncnd)+KRAD1,1)=Y(IX+(K-Ncnd)+KRAD,1)
09896 CONTINUE
09899 CONTINUE
      KRAD=KRAD1
      N=NZON*NVARS+2*KRAD
      JSTART=0
      Hmax=max(Hmax,t*1.d-3)
      NEEDBR=.TRUE.
      WRITE(4,*)' NEW NFRUS KRAD N-->',NFRUS,KRAD,N,' AT STEP=',NSTEP
      CallNTH(Ncnd)
      KTAIL=-50
      ELSEIF(LFR.LT.NFREQ)THEN
09893 IF(.NOT.(BLACK(LFR+1,max(Tpcomp,floor(3)))*FREQMN(LFR+1)**3.GT.FLNOIS))GOTO09892
      LFR=LFR+1
      IF(.NOT.(LFR.EQ.NFREQ))GOTO09893
09892 CONTINUE
      IF(LFR.NE.NFRUS)THEN
      WRITE(4,*)'Tpcomp=',Tpcomp
      KRAD1=(NZON-Ncnd)*LFR
      DO09890L=NFRUS,1,-1
      IX=NVARS*NZON+(NZON-Ncnd)*(L-1)
      DO09887K=NZON-1,Ncnd+1,-1
      Y(IX+(K-Ncnd)+KRAD1,1)=Y(IX+(K-Ncnd)+KRAD,1)
09887 CONTINUE
09890 CONTINUE
      DO09884K=Ncnd+1,NZON
      DO09881L=NFRUS+1,LFR
      IX=NVARS*NZON+(NZON-Ncnd)*(L-1)
      Y(IX+(K-Ncnd),1)=Y(IX+(K-Ncnd)-(NZON-Ncnd),1)*EXP((FREQMN(L-1)-FREQMN(L))/Tpcomp)
      Y(IX+(K-Ncnd)+KRAD1,1)=0.D0
09881 CONTINUE
09884 CONTINUE
      NFRUSO=NFRUS
      NFRUS=LFR
      KRAD=KRAD1
      N=NZON*NVARS+2*KRAD
      JSTART=0
      Hmax=max(Hmax,t*1.d-3)
      NEEDBR=.TRUE.
      WRITE(4,*)' NEW NFRUS KRAD N-->',NFRUS,KRAD,N,' AT STEP=',NSTEP
      CallNTH(Ncnd)
      KTAIL=-50
      ENDIF
      ENDIF
      DO09878L=1,NFREQ
      DO09875IK=1,NZON
      Eddn(IK,L)=tret
      EddO(IK,L)=Eddn(IK,L)
      EDDJ(IK,L)=Eddn(IK,L)
09875 CONTINUE
      HEdN(L)=0.5d0
      HEdO(L)=HEdN(L)
      HEDD(L)=HEDN(L)
09878 CONTINUE
      If(Givdtl)Then
      Open(unit=18,file='details.res')
      Endif
      WRITE(4,'(1X,''   ====> STEP:'',I6,2X,A,2X,'' ON  '',A)')NSTEP,CT001,D0001
      CALLFeau(Ncnd,Ncnd)
      trlx=25.d0*h
      DO09872L=1,Nfrus
      DO09869IK=Ncnd+1,NZON
      EddO(IK,L)=EddN(IK,L)
      EddJ(IK,L)=EddN(IK,L)
09869 CONTINUE
      HEdO(L)=HEdN(L)
      HEdd(L)=HEdN(L)
09872 CONTINUE
      DO09866L=Nfrus+1,Nfreq
      HEdO(L)=HEdN(Nfrus)
      HEdN(L)=HEdN(Nfrus)
      HEdd(L)=HEdN(Nfrus)
09866 CONTINUE
      CALLGDEPOS
      NEEDBR=.TRUE.
      IF(NSTEP.EQ.0.or.NOUT.GE.200)THEN
      ELOST=YENTOT(1)
      doK=1,Nzon
      Ry(K)=Y(K,1)
      Uy(K)=Y(Nzon+K,1)
      Ty(K)=Y(2*Nzon+K,1)
      enddo
      Rho(1)=3.D0*DM(1)/(Ry(1)**3-Rce**3)
      doI=2,Nzon
      Rho(I)=3.D0*DM(I)/(Ry(I)**3-Ry(I-1)**3)
      enddo
      tout=t
      doI=1,N
      Yout(I)=Y(I,1)
      enddo
      tretard=Ry(NZON-1)/CLIGHT
      t_ob=(t-tretard)*UTIME
      CALLPRIBAL(IOUT)
      ENDIF
      I=NZON
      TAU(I)=0.D0
09863 IF(.NOT.(I.GT.1.AND.I.GE.Ncnd))GOTO09862
      TP=Ty(I)
      PL=3.D0*DM(I)/(Ry(I)**3-Ry(I-1)**3)
      DO09860ii=1,Natom
      Yat(ii)=YABUN(ii,I)
09860 CONTINUE
      RADP=.FALSE.
      CALLURSOS
      kmhap=I
      CALLHAPPA
      TAUste=(Ry(I)-Ry(I-1))*HAPPAL(LubvB)*PL
      If(Scat)then
      TAUabs=(Ry(I)-Ry(I-1))*HAPabs(LubvB)*PL
      TAUsca=Tauste-Tauabs
      If(Tausca.GT.1.d0)then
      Tauste=sqrt(Tauabs*Tausca)
      else
      Tauste=Tauabs
      endif
      endif
      TAU(I-1)=TAU(I)+Tauste
      I=I-1
      GOTO09863
09862 CONTINUE
      DO09857kko=1,komax
      TO(kko)=TO(kko)*8.64d+04
09857 CONTINUE
      KFLAG=0
      KNTO=1
      nmom=1
      TOO=TO(1)
      KO=nmom
09854 IF(.NOT.(KFLAG.EQ.0.and.NSTEP.LT.NSTMAX))GOTO09853
      doI=1,N
      Ysave(I)=Y(I,1)
      enddo
      doIFL=1,NVARS
      doI=1,NZON
      YMAX(I+(IFL-1)*NZON)=MAX(ABS(Y(I+(IFL-1)*NZON,1)),FLOOR(IFL))*Wacc(IFL)
      enddo
      enddo
      DO09851IFL=1,Nfrus
      DO09848I=1,NZON-Ncnd
      If(Freqmn(IFL).LT.1.d0)Then
      YMAX(I+Nvars*Nzon+(IFL-1)*(NZON-Ncnd))=MAX(ABS(Y(I+Nvars*Nzon+(IFL-1)*(NZON-Ncnd),1)),FLOOR(Nvars+1))
      else
      YMAX(I+Nvars*Nzon+(IFL-1)*(NZON-Ncnd))=MAX(ABS(Y(I+Nvars*Nzon+(IFL-1)*(NZON-Ncnd),1)),FLOOR(Nvars+1)/(5.d0*Freqmn(IFL)**3))
      endif
09848 CONTINUE
09851 CONTINUE
      DO09845I=NZON*NVARS+Krad+1,N
      YMAX(I)=MAX(ABS(Y(I,1)),FLOOR(NVARS+1))
09845 CONTINUE
      tday=t*UTIME/8.64d+04
      if(NSTEP.EQ.0)t_eve=tday
      IF(.NOT.(abs(KNadap).GE.4))GOTO09842
      tdlog=log(max(tday,hmin))
      IF(.NOT.(tdlog.LE.stmlog(1)-(stmlog(2)-stmlog(1))/2.d0.or.tdlog.GT.stmlog(6).or.abs(Knadap).EQ.6))GOTO09839
      istim=0
      GOTO09837
09839 CONTINUE
      IF(.NOT.(tdlog.GT.stmlog(1)-(stmlog(2)-stmlog(1))/2.d0.AND.tdlog.LE.stmlog(1)))GOTO09838
      istim=1
      GOTO09837
09838 CONTINUE
      istim=2
09836 IF(.NOT.(stmlog(istim).LT.tdlog))GOTO09835
      istim=istim+1
      GOTO09836
09835 CONTINUE
09837 CONTINUE
      if(istim.NE.istold)then
      istold=istim
      GOTO(09833,09832),istim+1
      if(istim.GT.1.and.istim.LE.6)then
      ihp=istim
      thaplog1=stmlog(istim-1)
      thaplog2=stmlog(istim)
      else
      write(*,*)' in strad istim=',istim
      stop' wrong istim in strad'
      endif
      GOTO09831
09833 CONTINUE
      ihp=6
      thaplog1=stmlog(1)
      thaplog2=stmlog(6)
      GOTO09831
09832 CONTINUE
      ihp=1
      thaplog1=stmlog(1)-(stmlog(2)-stmlog(1))/2.d0
      thaplog2=stmlog(1)
09831 CONTINUE
      DO09830im=1,Nzon/(Mzon/50)
      DO09827iro=1,14
      DO09824itp=1,14
      DO09821L=1,Nfreq
      hpbansc1(L,itp,iro,im)=hpbansc2(L,itp,iro,im)
      hpbansc2(L,itp,iro,im)=hpsavsc(L,itp,iro,im,ihp)
09821 CONTINUE
09824 CONTINUE
09827 CONTINUE
09830 CONTINUE
      endif
09842 CONTINUE
      CALLVTSTIF
      Badste=.true.
09818 IF(.NOT.(Badste))GOTO09817
      Badste=.false.
      K=1
09815 IF(.NOT.(.not.Badste.and.K.LE.Nzon))GOTO09814
      Ry(K)=Y(K,1)
      Uy(K)=Y(Nzon+K,1)
      if(K.GT.1)then
      if(Ry(K).LE.Ry(K-1))then
      write(*,*)' BADSTE stradsep5tt Km=',K-1,'  R0=',Ry(K-1),'  u0=',Uy(K-1)
      write(*,*)' BADSTE stradsep5tt Km=',K,'  R1=',Ry(K),'  u1=',Uy(K)
      Badste=.true.
      endif
      endif
      Ty(K)=Y(2*Nzon+K,1)
      if(Ty(K).LE.0.d0)then
      write(*,*)' stradsep5tt: bad step - T, Km=',K,'  R1=',Ry(K),'  u1=',Uy(K),'  Ty=',Ty(K)
      Badste=.true.
      endif
      K=K+1
      GOTO09815
09814 CONTINUE
      if(Badste)then
      t=t-Hused
      doI=1,N
      Y(I,1)=Ysave(I)
      enddo
      H=Hused*1.d-1
      Jstart=0
      CALLVTSTIF
      Badste=.true.
      endif
      GOTO09818
09817 CONTINUE
      NSTEP=NSTEP+1
      NZM=NZMOD
      if(tday.GT.2.d0)CALLGDEPOS
      Rho(1)=3.D0*DM(1)/(Ry(1)**3-Rce**3)
      DO09812I=2,NZON
      Rho(I)=3.D0*DM(I)/(Ry(I)**3-Ry(I-1)**3)
09812 CONTINUE
      If(Givdtl)Then
      write(18,'(a)')'  nstep, Ncnd, Ty(max(Ncnd,1)), Ty(Ncnd+1)'
      write(18,'(2I5,1P,4E12.4)')nstep,Ncnd,Ty(max(Ncnd,1)),Ty(Ncnd+1)
      write(18,'(a,1P,4E12.4)')' k=max(Ncnd,1)    L=40 FJ Bb: ',Y(NVARS*NZON+(NZON-Ncnd)*(40-1)-Ncnd+max(Ncnd,1),1),Black(40,Ty(max(
     *Ncnd,1)))
      write(18,'(a,1P,4E12.4)')' k=Ncnd+1 L=40 FJ Bb: ',Y(NVARS*NZON+(NZON-Ncnd)*(40-1)+1,1),Black(40,Ty(Ncnd+1))
      Endif
      RX=1.D0
      IF(HUSED.NE.HOLDBL)THEN
      DO09809J=2,NQused+1
      RX=RX*(HUSED/HOLDBL)
      YENTOT(J)=YENTOT(J)*RX
09809 CONTINUE
      endif
      DO09806J1=1,NQUSED
      DO09803J2=J1,NQUSED
      J=NQUSED-J2+J1
      YENTOT(J)=YENTOT(J)+YENTOT(J+1)
09803 CONTINUE
09806 CONTINUE
      doK=1,Nzon
      Ry(K)=Y(K,1)
      Uy(K)=Y(Nzon+K,1)
      Ty(K)=Y(2*Nzon+K,1)
      enddo
      Rho(1)=3.D0*DM(1)/(Ry(1)**3-Rce**3)
      doI=2,Nzon
      Rho(I)=3.D0*DM(I)/(Ry(I)**3-Ry(I-1)**3)
      enddo
      CALLLOSSEN
      ERLUM=HUSED*ELTOT-YENTOT(2)
      CALLCOSET(METH,max(NQUSED,1),EL,TQ,MAXDER,IDOUB)
      DO09800J=1,NQused+1
      YENTOT(J)=YENTOT(J)+EL(J)*ERLUM
09800 CONTINUE
      IF(IABS(JSTART).GT.NQUSED)YENTOT(NQused+1+1)=ERLUM*EL(NQused+1)/FLOAT(NQused+1)
      HOLDBL=HUSED
      CallNTH(Ncnd)
      I=NZON
      TAU(I)=0.D0
09797 IF(.NOT.(I.GT.1))GOTO09796
      TP=Ty(I)
      Pl=3.D0*DM(I)/(Ry(I)**3-Ry(I-1)**3)
      DO09794ii=1,Natom
      Yat(ii)=YABUN(ii,I)
09794 CONTINUE
      RADP=.FALSE.
      CALLURSOS
      Rho(I)=Pl
      Press(I)=P*UPI
      kmhap=I
      CALLHAPPA
      Hapmin=Happal(1)
      Ltau=1
      DO09791L=2,Nfrus
      If(Happal(L).LT.Hapmin)then
      Hapmin=Happal(L)
      Ltau=L
      endif
09791 CONTINUE
      TAUste=(Ry(I)-Ry(I-1))*HAPPAL(LTAU)*PL
      If(Scat)then
      TAUabs=(Ry(I)-Ry(I-1))*HAPabs(LTAU)*PL
      TAUsca=Tauste-Tauabs
      If(Tausca.GT.1.d0)then
      Tauste=sqrt(Tauabs*Tausca)
      else
      Tauste=Tauabs
      endif
      endif
      TAU(I-1)=TAU(I)+Tauste
      I=I-1
      GOTO09797
09796 CONTINUE
      Ncnew=Nthick(Ltau)
      DO09788K=Ncnew-1,1,-1
      tau(K)=tau(K)-tau(Ncnew)
09788 CONTINUE
09785 IF(.NOT.(TAU(max(Ncnew,1)).LE.FitTau.AND.Ncnew.GE.NCNMIN))GOTO09784
      Ncnew=Ncnew-1
      GOTO09785
09784 CONTINUE
      if(mod(NSTEP,50).EQ.0)then
      WRITE(*,'(a,I6,99(2x,a,1P,G15.5))',advance='no')'step',NSTEP,'day',tday,'v kms',Uy(Nzon)*1d3,'R e14',Ry(Nzon)
      endif
      IF(Ncnew.LT.NCNMIN)Ncnew=0
      Ncnew=MIN(Ncnew,NZON-3)
      HOLDFR=.TRUE.
      IF(.NOT.((TAU(max(Ncnd,1)).GT.TauTol*FitTau.OR.TAU(max(Ncnd,1))*TauTol.LT.FitTau).AND.Ncnew.NE.Ncnd.AND.(CHNCND.OR.Nstep.EQ.0)
     *))GOTO09782
      IF(.NOT.(Ncnew.LT.Ncnd))GOTO09779
      DO09775L=1,Nfrus
      DO09772IK=Ncnd+1,NZON
      EddO(IK,L)=EddJ(IK,L)
09772 CONTINUE
      HEdO(L)=HEdd(L)
09775 CONTINUE
      CALLFeau(Ncnew,Ncnd)
      LFR=NFRUS
      Tpcomp=Ty(NZON)
      DO09769I=Ncnew+1,NZON
      If(Ty(I).GT.Tpcomp)Tpcomp=Ty(I)
09769 CONTINUE
09766 IF(.NOT.(BLACK(MIN(LFR+1,NFREQ),Tpcomp)*FREQMN(MIN(LFR+1,NFREQ))**3.GT.FLNOIS))GOTO09765
      LFR=LFR+1
      IF(.NOT.(LFR.GE.NFREQ))GOTO09766
09765 CONTINUE
      LFR=min(LFR,NFREQ)
      IF(LFR.NE.NFRUS)THEN
      WRITE(4,*)'Tpcomp=',Tpcomp
      WRITE(4,*)' NEW NFRUS -->',LFR,' AT STEP=',NSTEP
      ENDIF
      KRAD1=(NZON-Ncnew)*LFR
      DO09763L=NFRUS,1,-1
      IX=NVARS*NZON+(NZON-Ncnd)*(L-1)
      DO09760K=NZON-1,Ncnd+1,-1
      Y(IX+(K-Ncnd)+(Ncnd-Ncnew)*L+KRAD1,1)=Y(IX+(K-Ncnd)+KRAD,1)
09760 CONTINUE
09763 CONTINUE
      DO09757L=NFRUS,1,-1
      IX=NVARS*NZON+(NZON-Ncnd)*(L-1)
      DO09754K=NZON,Ncnd+1,-1
      Y(IX+(K-Ncnd)+(Ncnd-Ncnew)*L,1)=Y(IX+(K-Ncnd),1)
09754 CONTINUE
09757 CONTINUE
      DO09751K=Ncnd+1,NZON
      DO09748L=NFRUS+1,LFR
      IX=NVARS*NZON+(NZON-Ncnew)*(L-1)
      Y(IX+(K-Ncnew),1)=0.D0
      Y(IX+(K-Ncnew)+KRAD1,1)=0.D0
09748 CONTINUE
09751 CONTINUE
      NFRUSO=NFRUS
      NFRUS=LFR
      CallNTH(Ncnew)
      DO09745L=NFRUSO+1,NFRUS
      DO09742IK=Ncnd+1,NZON
      EddO(IK,L)=Eddn(IK,L)
      EDDJ(IK,L)=Eddn(IK,L)
09742 CONTINUE
      HEdO(L)=HEDD(L)
09745 CONTINUE
      trlx=25.d0*h
      CALLGDEPOS
      KTAIL=NSTEP
      Tp=Ty(Ncnew+1)
      If(Ncnew.GT.0)Then
      PL=3.D0*DM(Ncnew+1)/(Ry(Ncnew+1)**3-Ry(Ncnew)**3)
      Else
      PL=3.D0*DM(Ncnew+1)/(Ry(Ncnew+1)**3-Rce**3)
      Endif
      RADP=.FALSE.
      DO09739ii=1,Natom
      Yat(ii)=YABUN(ii,Ncnew+1)
09739 CONTINUE
      CALLURSOS
      kmhap=Ncnew+1
      CALLOPACIT
      CAP2=CAPPA
      DO09736K=Ncnew+1,Ncnd
      If(K.LT.Nzon)Then
      TP1=Tp
      CAP1=CAP2
      Tp=Ty(K+1)
      PL=3.D0*DM(K+1)/(Ry(K+1)**3-Ry(K)**3)
      DO09733ii=1,Natom
      Yat(ii)=YABUN(ii,K+1)
09733 CONTINUE
      CALLURSOS
      kmhap=K+1
      CALLOPACIT
      callHAPPA
      CAP2=CAPPA
      FLLF1=((Ry(K)*TP1)**4-(Ry(K)*TP)**4)*TP1**4/(CAP1*(DM(K)+DM(K+1))*(TP1**4+TP**4))
      FLRT1=((Ry(K)*TP1)**4-(Ry(K)*TP)**4)*TP**4/(CAP2*(DM(K)+DM(K+1))*(TP1**4+TP**4))
      FLUM=CLUM*(FLLF1+FLRT1)
      Endif
      If(NSTEP.EQ.0.AND.EddJ(K,LubvB).GT.0.25D0)Then
      dillm=0.75d0-0.25d0*sqrt(12.d0*EddJ(K,LubvB)-3.D0)
      Tpm=Freqmn(LubvB)/LOG(1.D0+dillm/FradJ(K,LubvB))
      else
      dillm=0.75d0
      TPM=MAX(Ty(K),2.D+3/UTP)
      endif
      FLUN=0.D0
      DO09730L=1,NFRUS
      IX=NVARS*NZON+(NZON-Ncnew)*(L-1)
      If(K.LE.Nthick(L))then
      Y(IX+(K-Ncnew),1)=Black(L,Ty(k))
      else
      Y(IX+(K-Ncnew),1)=FradJ(K,L)
      endif
      If(K.LT.Nzon)Then
      If(Nstep.GT.1.AND.K.GT.Nthick(L))then
      Y(IX+(K-Ncnew)+KRAD1,1)=FradH(K,L)
      elseif(Nstep.GT.1)then
      Y(IX+(K-Ncnew)+KRAD1,1)=(FLLF1+FLRT1)*Blackd(L,0.5D0*(Tp+Tp1))*CAP2/(6.D0*(0.5D0*(Tp+Tp1))**4*Ry(K)**2*HAPPAL(L))
      else
      Y(IX+(K-Ncnew)+KRAD1,1)=0.1d0*FradH(K,L)
      endif
      FLUN=FLUN+Y(IX+(K-Ncnew)+KRAD1,1)*WEIGHT(L)
      Endif
09730 CONTINUE
      If(K.LT.Nzon)Then
      FLUN=CLUMF*FLUN*Ry(K)**2
      WRITE(4,'(2(A,1P,E12.5),A,I3)')' Cond Lum=',Flum,'      New  Lum=',Flun,'   Km=',k
      Endif
09736 CONTINUE
      Ncnd=Ncnew
      KRAD=KRAD1
      N=NZON*NVARS+2*KRAD
      WRITE(4,*)' NEW Ncnd KRAD N:   ',Ncnd,KRAD,N,' AT STEP=',NSTEP
      HOLDFR=.FALSE.
      GOTO09776
09779 CONTINUE
      IF(.NOT.(Ncnew.GT.Ncnd.AND.Ncnew.LE.NZON-3))GOTO09778
      DO09727L=1,Nfrus
      DO09724IK=Ncnd+1,NZON
      EddO(IK,L)=EddJ(IK,L)
09724 CONTINUE
      HEdO(L)=HEdd(L)
09727 CONTINUE
      CALLFeau(Ncnew,Ncnd)
      KRAD1=(NZON-Ncnew)*NFRUS
      DO09721L=1,NFRUS
      IX=NVARS*NZON+(NZON-Ncnd)*(L-1)
      DO09718K=Ncnew+1,NZON
      Y(NVARS*NZON+(NZON-Ncnew)*(L-1)+(K-Ncnew),1)=Y(IX+(K-Ncnd),1)
09718 CONTINUE
09721 CONTINUE
      DO09715L=1,NFRUS
      IX=NVARS*NZON+(NZON-Ncnd)*(L-1)
      DO09712K=Ncnew+1,NZON-1
      Y(NVARS*NZON+(NZON-Ncnew)*(L-1)+(K-Ncnew)+KRAD1,1)=Y(IX+(K-Ncnd)+KRAD,1)
09712 CONTINUE
09715 CONTINUE
      trlx=25.d0*h
      CALLGDEPOS
      KTAIL=NSTEP
      Ncnd=Ncnew
      KRAD=KRAD1
      N=NZON*NVARS+2*KRAD
      WRITE(4,*)' NEW Ncnd KRAD N:   ',Ncnd,KRAD,N,' AT STEP=',NSTEP
      WRITE(4,*)'Ty(Ncnd+1)=',Ty(Ncnd+1)
      GOTO09776
09778 CONTINUE
      IF(.NOT.(Ncnew.GT.NZON-3))GOTO09777
      KFLAG=-3
      GOTO09776
09777 CONTINUE
09776 CONTINUE
      If(KFLAG.EQ.0)Then
      JSTART=0
      NEEDBR=.TRUE.
      endif
09782 CONTINUE
      If(KFLAG.EQ.0)then
      NTRANS=NZON
09709 IF(.NOT.(TAU(NTRANS).LT.3.D-4.AND.NTRANS.GT.Ncnd+1))GOTO09708
      NTRANS=NTRANS-1
      GOTO09709
09708 CONTINUE
      NTRANS=NZON-NTRANS
      IF(.NOT.(Nzon.GT.Nzmin.AND.NTRANS.GE.20.AND.T*Utime.GT.TREZON))GOTO09706
      NDELZN=NTrans
      WRITE(4,'('' ==> Dropping '',I3,'' zones at STEP:'',I5)')Ndelzn,Nstep
      WRITE(*,'('' ==> Dropping '',I3,'' zones at STEP:'',I5)')Ndelzn,Nstep
      Nznew=max(Nzmin,Nzon-Ndelzn)
      DO09703K=1,Nznew
      Y(Nznew+K,1)=Uy(K)
      Y(2*Nznew+K,1)=Ty(K)
09703 CONTINUE
      KRAD1=(Nznew-Ncnd)*NFRUS
      DO09700L=2,NFRUS
      IX=NVARS*NZON+(NZON-Ncnd)*(L-1)
      DO09697K=Ncnd+1,Nznew
      Y(NVARS*Nznew+(Nznew-Ncnd)*(L-1)+(K-Ncnd),1)=Y(IX+(K-Ncnd),1)
09697 CONTINUE
09700 CONTINUE
      DO09694L=2,NFRUS
      IX=NVARS*NZON+(NZON-Ncnd)*(L-1)
      DO09691K=Ncnd+1,Nznew-1
      Y(NVARS*Nznew+(Nznew-Ncnd)*(L-1)+(K-Ncnd)+KRAD1,1)=Y(IX+(K-Ncnd)+KRAD,1)
09691 CONTINUE
09694 CONTINUE
      Nzon=Nznew
      KRAD=KRAD1
      N=NZON*NVARS+2*KRAD
      CallNTH(Ncnd)
      PRINT*,' @@@>>> New Nzon Krad N:   ',Nzon,KRAD,N
      JSTART=0
      Hmax=max(Hmax,t*1.d-3)
      DO09688L=1,Nfrus
      DO09685IK=Ncnd+1,NZON
      EddO(IK,L)=EddJ(IK,L)
09685 CONTINUE
      HEdO(L)=HEdd(L)
09688 CONTINUE
      CALLFeau(Ncnd,Ncnd)
      trlx=25.d0*hCALLGDEPOS
      NEEDBR=.TRUE.
09706 CONTINUE
      endif
      Tpcomp=Ty(NZON)
      DO09682I=Ncnd+1,NZON
      If(Ty(I).GT.Tpcomp)Tpcomp=Ty(I)
09682 CONTINUE
      LFR=NFRUS
09679 IF(.NOT.(BLACK(LFR,max(Tpcomp,floor(3)))*FREQMN(LFR)**3.LT.FLNOIS/1.6D0.AND..NOT.HOLDFR.AND.Lfr.GT.NFRMIN))GOTO09678
      LFR=LFR-1
      GOTO09679
09678 CONTINUE
      IF(LFR.LT.NFRMIN)LFR=NFRMIN
      IF(LFR.NE.NFRUS)THEN
      NFRUS=LFR
      KRAD1=(NZON-Ncnd)*NFRUS
      DO09676L=1,NFRUS
      IX=NVARS*NZON+(NZON-Ncnd)*(L-1)
      DO09673K=Ncnd+1,NZON-1
      Y(IX+(K-Ncnd)+KRAD1,1)=Y(IX+(K-Ncnd)+KRAD,1)
09673 CONTINUE
09676 CONTINUE
      KRAD=KRAD1
      N=NZON*NVARS+2*KRAD
      JSTART=0
      Hmax=max(Hmax,t*1.d-3)
      NEEDBR=.TRUE.
      WRITE(4,*)' NEW NFRUS KRAD N-->',NFRUS,KRAD,N,' AT STEP=',NSTEP
      CallNTH(Ncnd)
      KTAIL=-50
      ELSEIF(LFR.LT.NFREQ)THEN
09670 IF(.NOT.(BLACK(LFR+1,max(Tpcomp,floor(3)))*FREQMN(LFR+1)**3.GT.FLNOIS))GOTO09669
      LFR=LFR+1
      IF(.NOT.(LFR.EQ.NFREQ))GOTO09670
09669 CONTINUE
      IF(LFR.NE.NFRUS)THEN
      WRITE(4,*)'Tpcomp=',Tpcomp
      KRAD1=(NZON-Ncnd)*LFR
      DO09667L=NFRUS,1,-1
      IX=NVARS*NZON+(NZON-Ncnd)*(L-1)
      DO09664K=NZON-1,Ncnd+1,-1
      Y(IX+(K-Ncnd)+KRAD1,1)=Y(IX+(K-Ncnd)+KRAD,1)
09664 CONTINUE
09667 CONTINUE
      DO09661K=Ncnd+1,NZON
      DO09658L=NFRUS+1,LFR
      IX=NVARS*NZON+(NZON-Ncnd)*(L-1)
      Y(IX+(K-Ncnd),1)=Y(IX+(K-Ncnd)-(NZON-Ncnd),1)*EXP((FREQMN(L-1)-FREQMN(L))/Tpcomp)
      Y(IX+(K-Ncnd)+KRAD1,1)=0.D0
09658 CONTINUE
09661 CONTINUE
      NFRUSO=NFRUS
      NFRUS=LFR
      KRAD=KRAD1
      N=NZON*NVARS+2*KRAD
      JSTART=0
      Hmax=max(Hmax,t*1.d-3)
      NEEDBR=.TRUE.
      WRITE(4,*)' NEW NFRUS KRAD N-->',NFRUS,KRAD,N,' AT STEP=',NSTEP
      CallNTH(Ncnd)
      KTAIL=-50
      ENDIF
      ENDIF
      If(N.GT.NYDIM)then
      WRITE(4,*)'N=',N,' > NYDIM'
      STOP390
      endif
      DO09655K=NZON*NVARS+1,NZON*NVARS+KRAD
      Y(K,1)=max(0.D0,Y(K,1))
09655 CONTINUE
      If((NSTEP-KTAIL.GE.50.OR.MOD(NSTEP,MBATCH).EQ.0.OR.MOD(NSTEP,50).EQ.0))then
      Ryzer=Ry(Kmcor)
      RTphi(Nzon+1)=0.d0
      EpsUq=0.d0
      DO09652K=1,NZON
      EpsUq=max(EpsUq,DRT*Eps*abs(Uy(K)))
      TauRT=(Amout-Am(K))*(0.4d0*Urho*UR)/Ry(K)**2
      RTphi(K)=BQ*dM(K)*Ry(K)**(NRT-1.d0)/(1.d0+0.033d0*TauRT)
09652 CONTINUE
      if(NSTEP.GT.0)then
      DO09649L=1,NFRUSO
      DO09646IK=Ncnd+1,NZON
      EddO(IK,L)=EddJ(IK,L)
09646 CONTINUE
      HEdO(L)=HEdd(L)
09649 CONTINUE
      endif
      CALLFeau(Ncnd,Ncnd)
      DO09643L=NFRUSO+1,NFRUS
      DO09640IK=Ncnd+1,NZON
      EddO(IK,L)=Eddn(IK,L)
      EDDJ(IK,L)=Eddn(IK,L)
09640 CONTINUE
      HEdO(L)=HEdn(L)
      HEDD(L)=HEDN(L)
09643 CONTINUE
      trlx=25.d0*h
      CALLGDEPOS
      KTAIL=NSTEP
      Ycomp=(Boltzk*UTp*Tpcomp)/(Amelec*Cs**2)
      NEEDBR=.TRUE.
      Evalja=.TRUE.
      endif
      If(JSTART.EQ.0)then
      doK=1,Nzon
      Ry(K)=Y(K,1)
      Uy(K)=Y(Nzon+K,1)
      Ty(K)=Y(2*Nzon+K,1)
      enddo
      Rho(1)=3.D0*DM(1)/(Ry(1)**3-Rce**3)
      doI=2,Nzon
      Rho(I)=3.D0*DM(I)/(Ry(I)**3-Ry(I-1)**3)
      enddo
      CALLLOSSEN
      YENTOT(2)=H*ELTOT
      HOLDBL=H
      endif
      tretard=Ry(NZON-1)/CLIGHT
      t_ob=(t-tretard)*UTIME
09637 IF(.NOT.(t_ob.GE.TOO.AND.KO.LE.KOMAX.AND.Nstep.GT.10))GOTO09636
      tout=TOO/UTIME+tretard
      DO09634I=1,N
      Yout(I)=0.D0
      RX=1.D0
      DO09631K=1,JSTART+1
      Yout(I)=Y(I,K)*RX+Yout(I)
      RX=RX*((TOO/UTIME+tretard-t)/H)
09631 CONTINUE
09634 CONTINUE
      DO09628I=1,NZON
      Ry(i)=Yout(i)
      Uy(I)=Yout(NZON+I)
      Ty(I)=Yout(2*NZON+I)
      if(Ty(I).LT.0.d0)then
      write(*,*)' I, Ty =',I,Ty(I)
      stop' strad %ALOF: Ty<0'
      endif
09628 CONTINUE
      RX=1.D0
      ELOST=0.D0
      DO09625K=1,JSTART+1
      ELOST=YENTOT(K)*RX+ELOST
      RX=RX*((TOO/UTIME-T)/H)
09625 CONTINUE
      Rho(1)=3.D0*DM(1)/(Ry(1)**3-Rce**3)
      DOI=2,NZON
      Rho(I)=3.D0*DM(I)/(Ry(I)**3-Ry(I-1)**3)
      ENDDO
      CALLPRIBAL(IOUT)
      CALLPRIBAL(+100)
      if(NTO(KO).NE.0)then
      if(KNTO.LE.NTO(KO))then
      TOO=TOO+STO(KO)
      KNTO=KNTO+1
      else
      KNTO=1
      KO=KO+1
      IF(KO.LE.KOMAX)TOO=TO(KO)
      endif
      else
      KO=KO+1
      IF(KO.LE.KOMAX)TOO=TO(KO)
      endif
      GOTO09637
09636 CONTINUE
      ELOST=YENTOT(1)
      if(MOD(NSTEP,NOUT).EQ.0)then
      WRITE(4,'(1X,''   ====> STEP:'',I6,2X,A,2X,'' ON  '',A)')NSTEP,CT001,D0001
      tout=t
      doI=1,N
      Yout(I)=Y(I,1)
      enddo
      doK=1,Nzon
      Ry(K)=Y(K,1)
      Uy(K)=Y(Nzon+K,1)
      Ty(K)=Y(2*Nzon+K,1)
      enddo
      Rho(1)=3.D0*DM(1)/(Ry(1)**3-Rce**3)
      doI=2,Nzon
      Rho(I)=3.D0*DM(I)/(Ry(I)**3-Ry(I-1)**3)
      enddo
      CALLPRIBAL(IOUT)
      CALLPRIMAG
      else
      if(Nstep.GT.1)CALLPRIMAG
      endif
      IF(KFLAG.EQ.0)THEN
      If(t*UTIME.GE.8.64D4*tcurB)NSTMAX=NSTEP
      IF(MOD(NSTEP,10).EQ.0)THEN
      ierr=0
      open(unit=77,file='stop',status='old',action='read',iostat=ierr)
      if(ierr.EQ.0)then
      NSTMAX=NSTEP
      write(*,*)'found file named stop, so terminating run now'
      endif
      ENDIF
      If((MOD(NSTEP,MBATCH).EQ.0.OR.NSTEP.EQ.NSTMAX).AND..NOT.Givdtl)then
      CALLSAVEM
      WRITE(4,'(1X,''   ====> SAVED:'',I6,2X,A,2X,'' ON  '',A)')NSTEP,CT001,D0001
      callTremain(MTIM,B0001)
      LTIM=MTIM
      endif
      ENDIF
      GOTO09854
09853 CONTINUE
      GOTO(09622,09621,09620,09619),-KFLAG
C      Write(6,*)' finish with KFLAG=',KFLAG
      GOTO09618
09622 CONTINUE
      Write(6,*)' the error test fails for H=HMIN: KFLAG=',KFLAG
      GOTO09618
09621 CONTINUE
      Write(6,*)' corrector iterations fail to converge for',' H=HMIN: KFLAG=',KFLAG
      Write(*,*)' corrector iterations fail to converge for',' H=HMIN: KFLAG=',KFLAG
      GOTO09618
09620 CONTINUE
      WRITE(4,*)' Ncnew=',Ncnew,'>=@NCOND_MAX'
      GOTO09618
09619 CONTINUE
      WRITE(4,*)' time limit: KFLAG=',KFLAG
09618 CONTINUE
      tretard=Ry(NZON-1)/CLIGHT
      t_ob=(t-tretard)*UTIME
      tout=t
      doK=1,Nzon
      Ry(K)=Y(K,1)
      Uy(K)=Y(Nzon+K,1)
      Ty(K)=Y(2*Nzon+K,1)
      enddo
      Rho(1)=3.D0*DM(1)/(Ry(1)**3-Rce**3)
      doI=2,Nzon
      Rho(I)=3.D0*DM(I)/(Ry(I)**3-Ry(I-1)**3)
      enddo
      doI=1,N
      Yout(I)=Y(I,1)
      enddo
      IF(KFLAG.LT.0)THEN
      CALLPRIBAL(0)
      ELSE
      CALLPRIBAL(Iout)
      ENDIF
      CALLSAVEM
      WRITE(4,'(1X,''   ====> SAVED:'',I6,2X,A,2X,'' ON  '',A)')NSTEP,CT001,D0001
      CALLVTIME(TIMEND)
      TIMEND=(TIMEND-B0001)
      WRITE(4,'('' RUN_TIME(s):'',1P,G12.3)')TIMEND
      write(*,'(a,f8.1)')'elaspsed time (minutes)',TIMEND/60.0
      WRITE(4,'(1X,''   ====> STEP:'',I6,2X,A,2X,'' ON  '',A)')NSTEP,CT001,D0001
      Close(unit=2)
      Close(unit=8)
      Close(unit=4)
      Close(unit=9)
      Close(unit=10)
      Close(unit=11)
      Close(unit=22)
      close(unit=23)
      close(unit=60)
      close(unit=81)
      close(unit=76)
      If(Givdtl)Then
      Close(unit=18)
      Endif
      if(KFLAG.EQ.0)then
      calltt4strad(runname)
      callLbol(runname)
      endif
      Isys=SYSTEM('mv  '//runname(1:lrun)//'.ph'//' ../../res/')
      Isys=SYSTEM('mv  '//runname(1:lrun)//'.tt'//' ../../res/')
      Isys=SYSTEM('mv  '//runname(1:lrun)//'.swd'//' ../../res/')
      Isys=SYSTEM('mv  '//runname(1:lrun)//'.lbol'//' ../../res/')
      Isys=SYSTEM('mv  '//runname(1:lrun)//'.res'//' ../../res/')
      If(Isys.EQ.-1)then
      errnum=ierrno()
      write(*,'(a,i2)')' Error mv to res= ',errnum
      endif
      Isys=SYSTEM('rm -f '//runname(1:lrun)//'.bm')
      Isys=SYSTEM('rm -f '//runname(1:lrun)//'.avg')
      If(Isys.EQ.-1)then
      errnum=ierrno()
      write(*,'(a,i2)')' Error  rm bm,avg = ',errnum
      endif
      STOP' Normal stop of stella'
      END
