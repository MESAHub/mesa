      SUBROUTINEBEGIN
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
      callStradIO('rm',Lunit,NFILE)
      EKO=EKO+EKOB
      If(EKO.GT.0.)then
      i=0
09994 CONTINUE
      i=i+1
      IF(.NOT.((am(i)-amini)/(amout-amini).GE.ai1.or.i.EQ.nzon))GOTO09994
      i1=i
09991 IF(.NOT.(i.LT.nzon))GOTO09990
      i=i+1
      IF(.NOT.((am(i)-amini)/(amout-amini).GE.ai2.or.i.EQ.nzon))GOTO09991
09990 CONTINUE
      i2=i
09988 IF(.NOT.(i.LT.nzon))GOTO09987
      i=i+1
      IF(.NOT.((am(i)-amini)/(amout-amini).GE.ai3.or.i.EQ.nzon))GOTO09988
09987 CONTINUE
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
      if(K.LT.Nzon)then
      Z1=Z1+Uy(K)**2*(DM(K)+DM(K+1))
      else
      Z1=Z1+Uy(K)**2*(DM(K)+DMout)
      endif
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
      km=1
      Endif
      kmnick=km-1
      AMNi=(AM(max(kmnick,1))-AMini)*UM
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
      If(ABS(Knadap).LE.3)then
      Open(unit=2,file=Opafile,form='unformatted')
      Else
      DO09949nf=1,6
      nu=30+nf
      app=char(ichar('0')+nf)
      Opafile=Opafile(1:Length(Opafile)-1)//app
      Open(unit=nu,file=Opafile,form='unformatted')
09949 CONTINUE
      Open(unit=29,file=Opafile(1:Length(Opafile)-1)//'ab',form='unformatted')
      endif
      GOTO(09946,09945,09944),indop(abs(KNadap))
      GOTO09943
09946 CONTINUE
      GOTO09943
09945 CONTINUE
      If(.NOT.petread)then
      read(2)Nfreq0,Msta,Nrho,Ntp,YATab,(Wavel(iif),iif=1,36),TpTab,RhoTab,STab,EpsBand,EppBand
      petread=.true.
      endif
      GOTO09943
09944 CONTINUE
      DO09942ihp=1,6
      read(30+ihp)nw,Stime,Nfreq0,Msta,Nrho,NTp,dumWavel,dumFreq,dumFreqmn,TpTab,RhoTab,hpbansc2
      DO09939im=1,Nzon/(Mzon/50)
      DO09936iro=1,Nrho
      DO09933itp=1,NTp
      DO09930L=1,Nfreq0
      hpsavsc(L,itp,iro,im,ihp)=hpbansc2(L,itp,iro,im)
09930 CONTINUE
09933 CONTINUE
09936 CONTINUE
09939 CONTINUE
      close(30+ihp)
      If(ihp.EQ.nw.OR.ABS(knadap).EQ.6.OR.ihp.EQ.1)then
      stmlog(ihp)=log(Stime)
      else
      write(*,*)' error reading haptab, ihp=',ihp,' nw=',nw
      stop324
      endif
09942 CONTINUE
      read(29)hpbanab2
      close(29)
09943 CONTINUE
      tdlog=log(max(t*Utime/8.6d+04,hmin))
      IF(.NOT.(tdlog.LE.stmlog(1)-(stmlog(2)-stmlog(1))/2.d0.or.tdlog.GT.stmlog(6).or.abs(Knadap).EQ.6))GOTO09927
      istim=0
      GOTO09925
09927 CONTINUE
      IF(.NOT.(tdlog.GT.stmlog(1)-(stmlog(2)-stmlog(1))/2.d0.AND.tdlog.LE.stmlog(1)))GOTO09926
      istim=1
      GOTO09925
09926 CONTINUE
      istim=2
09924 IF(.NOT.(stmlog(istim).LT.tdlog))GOTO09923
      istim=istim+1
      GOTO09924
09923 CONTINUE
09925 CONTINUE
      istold=istim
      GOTO(09921,09920),istim+1
      if(istim.GT.1.and.istim.LE.6)then
      thaplog1=stmlog(istim-1)
      thaplog2=stmlog(istim)
      else
      write(*,*)' in begrad istim=',istim
      stop' wrong istim in begrad'
      endif
      GOTO09919
09921 CONTINUE
      DO09918im=1,Nzon/(Mzon/50)
      DO09915iro=1,14
      DO09912itp=1,14
      DO09909L=1,Nfreq
      hpbanab1(L,itp,iro,im)=hpbanab2(L,itp,iro,im)
      hpbansc1(L,itp,iro,im)=hpbansc2(L,itp,iro,im)
09909 CONTINUE
09912 CONTINUE
09915 CONTINUE
09918 CONTINUE
      thaplog1=stmlog(1)
      thaplog2=stmlog(6)
      GOTO09919
09920 CONTINUE
      DO09906im=1,Nzon/(Mzon/50)
      DO09903iro=1,14
      DO09900itp=1,14
      DO09897L=1,Nfreq
      hpbanab1(L,itp,iro,im)=hpbanab2(L,itp,iro,im)
      hpbansc1(L,itp,iro,im)=hpbansc2(L,itp,iro,im)
09897 CONTINUE
09900 CONTINUE
09903 CONTINUE
09906 CONTINUE
      thaplog1=stmlog(1)-(stmlog(2)-stmlog(1))/2.d0
      thaplog2=stmlog(1)
09919 CONTINUE
      If(istim.NE.0)then
      app=char(ichar('0')+istim)
      Opafile=Opafile(1:Length(Opafile)-1)//app
      Open(unit=30,file=Opafile,form='unformatted')
      read(30)nw,Stime,Nfreq0,Msta,Nrho,NTp,dumWavel,dumFreq,dumFreqmn,TpTab,RhoTab,hpbansc2
      close(30)
      endif
      If(istim.NE.0.AND.istim.NE.1)then
      app=char(ichar('0')+istim-1)
      Opafile=Opafile(1:Length(Opafile)-1)//app
      Open(unit=30,file=Opafile,form='unformatted')
      read(30)nw,Stime,Nfreq0,Msta,Nrho,NTp,dumWavel,dumFreq,dumFreqmn,TpTab,RhoTab,hpbansc1
      close(30)
      close(29)
      Open(unit=29,file=Opafile(1:Length(Opafile)-1)//'ab',form='unformatted')
      read(29)hpbanab1
      close(29)
      endif
      DO09894im=1,Nzon/(Mzon/50)
      DO09891iro=1,14
      DO09888itp=1,14
      DO09885L=1,Nfreq
      hbab=hpbanab2(L,itp,iro,im)
      hbal=hpbansc2(L,itp,iro,im)
      if(hbal.GT.1.d50.or.hbab.GT.1.d50)then
      write(*,'(4(a,i4),1p,2(a,e12.4))')' L=',L,' itp=',itp,' iro=',iro,' im=',im,'  hbal=',hbal,'  hbab=',hbab
      endif
09885 CONTINUE
09888 CONTINUE
09891 CONTINUE
09894 CONTINUE
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
      WRITE(4,'('' XMNI  = '',1P,E15.5,A,'' NRT   = '',G10.3)')XMNI,' SOL.MASS',NRT
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
      YENTOT(2)=0.d0
      GOTO(09882,09881,09880),indfr(abs(KNadap))
      GOTO09879
09882 CONTINUE
      FREQ(1)=CS*1.D+08/(WLMAX*UFREQ)
      FREQ(NFREQ+1)=CS*1.D+08/(WLMIN*UFREQ)
      Basis=(WLMAX/WLMIN)**(1.D0/(DBLE(NFREQ)))
      DO09878i=2,NFREQ
      FREQ(i)=FREQ(i-1)*BASIS
09878 CONTINUE
      DO09875L=1,NFREQ
      FREQMN(L)=SQRT(FREQ(L)*FREQ(L+1))
09875 CONTINUE
      GOTO09879
09881 CONTINUE
      If(.NOT.petread)then
      read(2)Nfreq0,Msta,Nrho,Ntp,YATab,(Wavel(iif),iif=1,36),TpTab,RhoTab,STab,EpsBand,EppBand
      petread=.true.
      endif
      If(Nfreq.GT.36)then
      DO09872L=36+1,Nfreq
      Wavel(L)=Wavel(L-1)/2.d0
09872 CONTINUE
      endif
      DO09869L=1,NFREQ
      FREQMN(L)=CS*1.D+08/(WaveL(L)*UFREQ)
09869 CONTINUE
      FREQ(1)=0.5d0*Freqmn(1)
      FREQ(Nfreq+1)=2.d0*Freqmn(Nfreq)
      DO09866L=2,NFREQ
      FREQ(L)=0.5d0*(Freqmn(L-1)+Freqmn(L))
09866 CONTINUE
      GOTO09879
09880 CONTINUE
      DO09863i=1,NFREQ+1
      FREQ(i)=dumFREQ(i)/Ufreq
09863 CONTINUE
      DO09860L=1,NFREQ
      FREQMN(L)=dumFREQMN(L)/Ufreq
09860 CONTINUE
09879 CONTINUE
      DO09857L=1,NFREQ
      WEIGHT(L)=(FREQ(L+1)-FREQ(L))*FREQMN(L)**3
      IF(L.LT.NFREQ)THEN
      DLOGNU(L)=1.D0/LOG(FREQMN(L)/FREQMN(L+1))
      ELSE
      DLOGNU(L)=1.D0/LOG(FREQMN(L)/FREQ(L+1))
      ENDIF
09857 CONTINUE
      NFRUS=NFREQ
      DO09854L=1,NFRUS
      NTHICK(L)=ncnd
09854 CONTINUE
      DO09851KM=NCND+1,NZON
      DO09848L=1,NFRUS
      LTHICK(L,KM)=.FALSE.
09848 CONTINUE
09851 CONTINUE
      ENDIF
      doL=1,Nfreq
      freqob(L)=freqmn(L)
      enddo
      BASIS=FREQ(2)/FREQ(1)
      If(Mfreq.GT.Nfreq)then
      FREQprev=FREQ(Nfreq+1)
      doL=Nfreq+1,Mfreq
      freqob(L)=freqmn(Nfreq)*exp(-dble(L-Nfreq)/dlognu(1))
      FREQnext=FREQprev*BASIS
      WEIGHT(L)=(FREQnext-FREQprev)*freqob(L)**3
      testFreqob=sqrt(FREQnext*FREQprev)
      FREQprev=FREQnext
      enddo
      endif
      toldg=-1.d0
      WRITE(4,'(''   FREQ:'',1P,10E12.5)')FREQ
      WRITE(4,'(''   FREQMN:'',1P,10E12.5)')FREQMN
      WRITE(4,'(''  WAVE BOUNDS:'',1P,11E10.3)')(CCL/FREQ(LW),LW=1,NFREQ+1)
      WRITE(4,'(''  WAVES:'',1P,10E12.5)')(CCL/FREQMN(LW),LW=1,NFREQ)
      WRITE(4,'('' WEIGHT:'',1P,10E12.5)')WEIGHT
      LubvU=NFRUS
09845 IF(.NOT.(LubvU.GT.1.AND.FREQMN(LubvU).GT.(CCL/3.65D+03)))GOTO09844
      LubvU=LubvU-1
      GOTO09845
09844 CONTINUE
      LubvB=NFRUS
09842 IF(.NOT.(LubvB.GT.1.AND.FREQMN(LubvB).GT.(CCL/4.4D+03)))GOTO09841
      LubvB=LubvB-1
      GOTO09842
09841 CONTINUE
      LubvV=NFRUS
09839 IF(.NOT.(LubvV.GT.1.AND.FREQMN(LubvV).GT.(CCL/5.5D+03)))GOTO09838
      LubvV=LubvV-1
      GOTO09839
09838 CONTINUE
      LubvR=NFRUS
09836 IF(.NOT.(LubvR.GT.1.AND.FREQMN(LubvR).GT.(CCL/7.D+03)))GOTO09835
      LubvR=LubvR-1
      GOTO09836
09835 CONTINUE
      LubvI=NFRUS
09833 IF(.NOT.(LubvI.GT.1.AND.FREQMN(LubvI).GT.(CCL/9.D+03)))GOTO09832
      LubvI=LubvI-1
      GOTO09833
09832 CONTINUE
      Lyman=NFRUS
09830 IF(.NOT.(Lyman.GT.1.AND.FREQMN(Lyman).GT.(CCL/912.D0)))GOTO09829
      Lyman=Lyman-1
      GOTO09830
09829 CONTINUE
      WRITE(4,*)' L UBVRI Lyman:',LubvU,LubvB,LubvV,LubvR,LubvI,Lyman
      LST=2
      KENTR=0
      RHO(1)=3.D0*DM(1)/(Ry(1)**3-Rce**3)
      DO09827I=2,NZON
      RHO(I)=3.D0*DM(I)/(Ry(I)**3-Ry(I-1)**3)
09827 CONTINUE
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
      SUBROUTINESAVEM
      IMPLICITREAL*8(A-H,O-Z)
      INTEGERLUNIT
      Character*80NFILE*80
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
      if(NSTB.LT.0)then
      Lunit=10
      NFILE=Sumprf
      else
      Lunit=12
      NFILE='test.prf'
      endif
      callStradIO('wm',Lunit,NFILE)
      if(NSTB.LT.0)then
      Lunit=11
      NFILE=Sumcur
      else
      Lunit=13
      NFILE='test.crv'
      endif
      callStradIO('wc',Lunit,NFILE)
      RETURN
      END
      BLOCKDATASTDATA
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
      COMMON/QNRGYE/QNUC,RGASA,YELECT
      COMMON/HNUSED/HUSED,NQUSED,NFUN,NJAC,NITER,NFAIL
      COMMON/AZNUC/ACARB,ZCARB,ASI,ZSI,ANI,ZNI,QCSI,QSINI
      COMMON/TOO/TOO,KO,KNTO,TO(KOMAX),STO(KOMAX),NTO(KOMAX)
      DATATO/komax*1.D+19/,STO/komax*1.D-5/,NTO/komax*0/
      DATAQNUC,YELECT/8.8861D+6,0.5D0/
      DATAACARB,ZCARB,ASI,ZSI,ANI,ZNI,QCSI,QSINI/12.D0,6.D0,28.D0,14.D0,56.D0,28.D0,88.861D0,1.884D0/
      DATANFUN,NJAC,NITER,NFAIL/4*0/
      END
