      SUBROUTINEDIFMAT
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
      common/NSTEP/NSTEP,NDebug,MAXER,IOUT,NOUT
      common/CREAD/TAUOLD,NSTMAX,MBATCH,MAXORD
      common/debug/LfrDebug,Nperturb,Kbad
      PARAMETER(ALFCON=1.D0)
      Common/Hydro/dm0,dm1,dm2,Rodm,RSodm,R0,R1,R2,RSm1,RS0,RS1,RS2,RC0,RC1,RC2,Rcinv0,RCinv1,RCinv2,U0,U1,U2,Q1,Q2,QU1,QU2,UDIVR,UU
     *1,UU2,UQ1,UQ2,AQ1,AQ2,QUHALF
      Common/Termo/PL0,PL1,PL2,ET1,ET2,Tp0,Tp1,Tp2,PQ1,PQ2,PQpl1,PQpl2,PT1,PT2,TPTQ1,TPTQ2,TpT1,TpT2,Tp0up4,Tp1up4,Tp2up4,PTPL1,PTPL
     *2,PTT1,PTT2,ETPL1,ETPL2,ETT1,ETT2,DMETIN,CDMETI,CK2ET1,HEATRA,FLUM,FIMP,DTPOLD,ENG0,ENG1,ENG2,ENGPL1,ENGPL2,ENGT1,ENGT2
      Common/Tratst/TOTJ,TOTW,TOTWPR
      Common/Convo/PP0,PP1,PP2,CP0,CP1,CP2,GRA0,GRA1,GRA2,VCSQ0,VCSQ1,FC0,FC1,ELMIX0,ELMIX1,CFC0,CFC1,CPPL0,CPPL1,CPPL2,CPT0,CPT1,CP
     *T2,GRAPL0,GRAPL1,GRAPL2,GRAT0,GRAT1,GRAT2,PRPL0,PRPL1,PRPL2,PRT0,PRT1,PRT2
      Common/Fluxo/tau1(Nfreq),tau2(Nfreq),FL0,FL1,Flcor1,Flcor2,Fllf0,Flrt0,Fllf1,Flrt1,FLcore,CAPcor,CAP0,CAP1,CAP2,CAP1D,CAP1T,CA
     *PT0,CAPT1,CAPT2,CAPPL0,CAPPL1,CAPPL2,CM0INV,CM1INV
      Common/Transo/DJNU(NFREQ),B24(NFREQ),DFH(NFREQ),DFJ(NFREQ),DUR(NFREQ),DFJDPL(NFREQ),DFHDPL(NFREQ),HAPL(NFREQ),HAPL1(NFREQ),HAP
     *L2(NFREQ),HAPL1T(NFREQ),HAPL2T(NFREQ),HAPL1D(NFREQ),HAPL2D(NFREQ),HAPH(NFREQ),HAPH1(NFREQ),HAAB(NFREQ),HAAB1(NFREQ),HAAB2(NFRE
     *Q),HAAB1T(NFREQ),HAAB2T(NFREQ),HAAB1D(NFREQ),HAAB2D(NFREQ),FHcore(NFREQ),BLCKD(NFREQ)
      Common/Burno/YDOT1,YDOT2,CCPL1,CCPL2,CCTP1,CCTP2,AS1,AS2,YCARB1,YCARB2
      Parameter(RADC3=RADC/3.D0)
      COMMON/PHOT/XJPH,DMPH,RPH,TPH,PLPH,VPH,CHEMPH,GRVPH,HP,JPH
      REAL*8Ry(Mzon),Uy(Mzon),Ty(Mzon)
      BLACK(Lbl,Tpbl)=(exp(-(FREQMN(Lbl)/Tpbl)))/(1.d0-(exp(-(FREQMN(Lbl)/Tpbl))))
      BLACKD(Lbl,Tpbl)=(FREQMN(Lbl)/Tpbl)*(exp(-(FREQMN(Lbl)/Tpbl)))/(1.d0-(exp(-(FREQMN(Lbl)/Tpbl))))**2
      Theta(Z)=0.5D0*(SIGN(1.D0,Z)+1.D0)
      GOTO09989
09999 CONTINUE
      Tprad4=0.D0
      TRD4=0.D0
      DO09986L=1,NFRUS
      If(LTHICK(L,KMth))Then
      Tprad4=Tprad4+BLACK(L,Tp)*WEIGHT(L)
      TRD4=TRD4+BLACKD(L,Tp)*WEIGHT(L)
      endif
09986 CONTINUE
      If(Tprad4.NE.0.D0)Then
      Tprad4=Tprad4*(15.D0/PI**4)
      TRD4=TRD4*(15.D0/PI**4)
      EPL=EPL-RADC*(Tprad4/PL)
      ET=ET+RADC*(TRD4/Tp)/PL
      PTth=PTth+RADC3*TRD4
      endif
      GOTO(09998,09997,09996,09995,09994,09993,09992,09991,09990),I6R001
09989 CONTINUE
      NFUN=NFUN+1
      Kbad=0
      DO09983K=1,Nzon
      Ry(K)=Fsave(K)
      If(K.GT.1)then
      if(Ry(K).LE.Ry(K-1))then
      Badste=.true.
      Kbad=K
      Return
      endif
      endif
      Uy(K)=Fsave(Nzon+K)
      Ty(K)=Fsave(2*Nzon+K)
      If(Ty(K).LE.3.d1/UTp)then
      Badste=.true.
      Return
      endif
09983 CONTINUE
      If(EVALJA)then
      CallDHapH
      endif
      IF(EDTM)THEN
      tw=min((t-tfeau)/trlx,1.d0)
      ww=tw**2*(3.d0-2.d0*tw)
      DO09980L=1,Nfrus
      DO09977IK=Ncnd+1,NZON
      EddJ(IK,L)=(1.d0-ww)*EddO(IK,L)+ww*EddN(IK,L)
09977 CONTINUE
      HEdd(L)=(1.d0-ww)*HEdO(L)+ww*HEdN(L)
09980 CONTINUE
      ELSE
      DO09974L=1,Nfrus
      DO09971IK=Ncnd+1,NZON
      EddJ(IK,L)=EddN(IK,L)
09971 CONTINUE
      HEdd(L)=HEdN(L)
09974 CONTINUE
      ENDIF
      Km=0
      JAC=0
      R1=RCE
      U1=0.D0
      FL1=0.D0
      FC1=0.D0
      RS1=R1**2
      RC1=RS1*R1
      DM1=0.D0
      Tp1=0.D0
      Tp1up4=0.D0
      RCinv1=0.D0
      Fllf1=0.D0
      Flrt1=0.D0
      CM1INV=0.D0
      DO09968ii=1,Natom
      Yat(ii)=YABUN(ii,1)
09968 CONTINUE
      R2=Ry(Km+1)
      RS2=R2**2
      RC2=RS2*R2
      U2=Uy(Km+1)
      DM2=DM(Km+1)
      RCinv2=3.D0/(RC2-RC1)
      PL2=DM2*RCinv2
      Tp2=Ty(Km+1)
      PL=PL2
      Tp=Tp2
      Tp2up4=Tp2**4
      If(Km.LE.Ncnd.AND.Ncnd.NE.0)then
      RADP=.TRUE.
      else
      RADP=.FALSE.
      endif
      CHEM=CHEM0(Km+1)
      CALLURSOS
      CALLVOLEN(Km+1)
      CHEM0(Km+1)=CHEM
      PTth=Tp*PT
      kmhap=1
      IF(.NOT.(Ncnd.GE.1))GOTO09965
      CALLOPACIT
      GOTO09964
09965 CONTINUE
      CALLOPACIT
      CALLHAPPA
      DO09962L=1,NFRUS
      hapL2(L)=hapPAL(L)
      HAAB2(L)=hapABS(L)
09962 CONTINUE
      KMth=1
      I6R001=7
      GOTO09999
09992 CONTINUE
09964 CONTINUE
      CHEM0(1)=CHEM
      UU2=U2-U1*(RS1/RS2)
      UQ2=UU2-EpsUq
      AQ2=AQ
      IF(UU2.LT.0.D0)THEN
      QUHALF=UU2*(AQ2*PL2)
      ELSE
      QUHALF=0.D0
      ENDIF
      Q2=QUHALF*UU2
      PQ2=P*UPI+Q2
      ENG2=ENG
      CAP2=CAPPA
      TPTQ2=PTth*UPI+Q2
      ET2=ET*UEI
      PT2=PT*UPI
      CP2=CP*UEI/PL2
      GRA2=DA
      IF(.NOT.(EVALJA))GOTO09959
      DO09956L=1,NFRUS
      B24(L)=0.D0
09956 CONTINUE
      PT2=PT*UPI
      ENGPL2=ENGPL
      ENGT2=ENGT
      PQPL2=PL2*((PPL*UPI)+UU2*MIN(UU2,0.D0)*AQ2)
      CAPPL1=0.D0
      CAPT1=0.D0
      QU2=2.D0*QUHALF
      CAPT2=CAPT
      CAPPL2=CAPPL
      PL=PL2*(1.D0+DELTA)
      Tp=Tp2
      CHEM=CHEM0(Km+1)
      CALLURSOS
      PTth=Tp*PT
      IF(.NOT.(Km+1.GT.Ncnd))GOTO09953
      kmhap=km+1
      CALLHAPPA
      DO09950L=1,NFRUS
      hapL2D(L)=hapPAL(L)
      HAAB2D(L)=hapABS(L)
09950 CONTINUE
      KMth=Km+1
      I6R001=8
      GOTO09999
09991 CONTINUE
09953 CONTINUE
      PTPL2=PTth*UPI
      ETPL2=ET*UEI
      PL=PL2
      Tp=Tp2*(1.D0+DELTA)
      CHEM=CHEM0(Km+1)
      CALLURSOS
      PTth=Tp*PT
      IF(.NOT.(Km+1.GT.Ncnd))GOTO09947
      kmhap=km+1
      CALLHAPPA
      DO09944L=1,NFRUS
      hapL2T(L)=hapPAL(L)
      HAAB2T(L)=hapABS(L)
09944 CONTINUE
      I6R001=9
      GOTO09999
09990 CONTINUE
09947 CONTINUE
      PTT2=PTth*UPI
      ETT2=ET*UEI
      Tp=Tp2
09959 CONTINUE
09941 IF(.NOT.(Km.LT.NZON))GOTO09940
      Km=Km+1
      If(Km.NE.1)then
      RSm1=RS0
      else
      RSm1=0.
      endif
      R0=R1
      R1=R2
      RS0=RS1
      RS1=RS2
      RC0=RC1
      RC1=RC2
      U0=U1
      U1=U2
      DM0=DM1
      DM1=DM2
      FL0=FL1
      UU1=UU2
      UQ1=UQ2
      Tp0=Tp1
      Tp1=Tp2
      Q1=Q2
      Tp0up4=Tp1up4
      Tp1up4=Tp2up4
      ENG0=ENG1
      ENG1=ENG2
      PL0=PL1
      PL1=PL2
      RCinv0=RCinv1
      RCinv1=RCinv2
      IF(.NOT.(Km.LE.Ncnd))GOTO09938
      ET1=ET2
      TPTQ1=TPTQ2
      CAP0=CAP1
      CAP1=CAP2
      PT1=PT2
      PQ1=PQ2
      PQPL1=PQPL2
      GOTO09936
09938 CONTINUE
      IF(.NOT.(Km.EQ.Ncnd+1))GOTO09937
      RADP=.FALSE.
      PL=PL1
      Tp=Tp1
      DO09935i=1,Natom
      Yat(i)=YABUN(I,Km)
09935 CONTINUE
      CHEM=CHEM0(Km)
      CALLURSOS
      kmhap=km
      CALLHAPPA
      PTth=Tp*PT
      KMth=Km
      I6R001=1
      GOTO09999
09998 CONTINUE
      PQ1=P*UPI+Q1
      ET1=ET*UEI
      TPTQ1=PTth*UPI+Q1
      CAP0=CAP1
      CAP1=CAP2
      UDIVR=(RS0*U0+RS1*U1)/(RC0+RC1)
      DO09932L=1,NFRUS
      DUR(L)=(U1-U0)/(R1-R0)*EDDJ(Km,L)+UDIVR*(1.D0-EDDJ(Km,L))
      hapL1(L)=hapPAL(L)
      HAAB1(L)=hapABS(L)
09932 CONTINUE
      GOTO09936
09937 CONTINUE
      ET1=ET2
      TPTQ1=TPTQ2
      PT1=PT2
      PQ1=PQ2
      PQPL1=PQPL2
      UDIVR=(RS0*U0+RS1*U1)/(RC0+RC1)
      DO09929L=1,NFRUS
      DUR(L)=(U1-U0)/(R1-R0)*EDDJ(Km,L)+UDIVR*(1.D0-EDDJ(Km,L))
      hapL1(L)=hapL2(L)
      HAAB1(L)=HAAB2(L)
09929 CONTINUE
09936 CONTINUE
      IF(.NOT.(EVALJA))GOTO09926
      CAPT0=CAPT1
      CAPPL0=CAPPL1
      QU1=QU2
      CM0INV=CM1INV
      AQ1=AQ2
      CAPT1=CAPT2
      CAPPL1=CAPPL2
      FLLF0=FLLF1
      FLRT0=FLRT1
      ENGPL1=ENGPL2
      ENGT1=ENGT2
      IF(.NOT.(Km.LE.Ncnd))GOTO09923
      PTPL1=PTPL2
      PTT1=PTT2
      ETPL1=ETPL2
      ETT1=ETT2
      GOTO09921
09923 CONTINUE
      IF(.NOT.(Km.EQ.Ncnd+1))GOTO09922
      PT1=PT*UPI
      PQPL1=PL*((PPL*UPI)+UU1*MIN(UU1,0.D0)*AQ1)
      RADP=.FALSE.
      PL=PL1*(1.D0+DELTA)
      Tp=Tp1
      CHEM=CHEM0(Km)
      CALLURSOS
      kmhap=km
      CALLOPACIT
      CALLHAPPA
      DO09920L=1,NFRUS
      hapL1D(L)=hapPAL(L)
      HAAB1D(L)=hapABS(L)
09920 CONTINUE
      PTth=Tp*PT
      KMth=Km
      I6R001=2
      GOTO09999
09997 CONTINUE
      PTPL1=PTth*UPI
      ETPL1=ET*UEI
      CAP1D=CAPPA
      PL=PL1
      CHEM=CHEM0(Km)
      Tp=Tp1*(1.D0+DELTA)
      CALLURSOS
      kmhap=km
      CALLOPACIT
      CALLHAPPA
      DO09917L=1,NFRUS
      hapL1T(L)=hapPAL(L)
      HAAB1T(L)=hapABS(L)
09917 CONTINUE
      PTth=Tp*PT
      I6R001=3
      GOTO09999
09996 CONTINUE
      PTT1=PTth*UPI
      ETT1=ET*UEI
      CAP1T=CAPPA
      GOTO09921
09922 CONTINUE
      PTPL1=PTPL2
      PTT1=PTT2
      ETPL1=ETPL2
      ETT1=ETT2
      DO09914L=1,NFRUS
      hapL1T(L)=hapL2T(L)
      HAAB1T(L)=HAAB2T(L)
      hapL1D(L)=hapL2D(L)
      HAAB1D(L)=HAAB2D(L)
09914 CONTINUE
09921 CONTINUE
09926 CONTINUE
      IF(.NOT.(Km.NE.Nzon))GOTO09911
      DO09908i=1,Natom
      Yat(i)=YABUN(I,Km+1)
09908 CONTINUE
      R2=Ry(Km+1)
      RS2=R2**2
      RC2=RS2*R2
      U2=Uy(Km+1)
      DM2=DM(Km+1)
      RCinv2=3.D0/(RC2-RC1)
      PL2=DM2*RCinv2
      Tp2=Ty(Km+1)
      PL=PL2
      Tp=Tp2
      Tp2up4=Tp2**4
      If(Km.LE.Ncnd.AND.Ncnd.NE.0)then
      RADP=.TRUE.
      else
      RADP=.FALSE.
      endif
      CHEM=CHEM0(Km+1)
      CALLURSOS
      CALLVOLEN(Km+1)
      CHEM0(Km+1)=CHEM
      PTth=Tp*PT
      kmhap=km+1
      IF(.NOT.(Km+1.LE.Ncnd))GOTO09905
      CALLOPACIT
      GOTO09904
09905 CONTINUE
      CALLOPACIT
      CALLHAPPA
      DO09902L=1,NFRUS
      hapL2(L)=hapPAL(L)
      HAAB2(L)=hapABS(L)
09902 CONTINUE
      KMth=Km+1
      I6R001=4
      GOTO09999
09995 CONTINUE
09904 CONTINUE
      UU2=U2-U1*(RS1/RS2)
      UQ2=UU2-EpsUq
      AQ2=AQ
      IF(UU2.LT.0.D0)THEN
      QUHALF=UU2*(AQ2*PL2)
      ELSE
      QUHALF=0.D0
      ENDIF
      Q2=QUHALF*UU2
      PQ2=P*UPI+Q2
      ENG2=ENG
      CAP2=CAPPA
      TPTQ2=PTth*UPI+Q2
      ET2=ET*UEI
      PT2=PT*UPI
      IF(.NOT.(Km.LT.Ncnd))GOTO09899
      FLLF1=((R1*Tp1)**4-(R1*Tp2)**4)*Tp1up4/(CAP1*(DM1+DM2)*(Tp1up4+Tp2up4))
      FLRT1=((R1*Tp1)**4-(R1*Tp2)**4)*Tp2up4/(CAP2*(DM1+DM2)*(Tp1up4+Tp2up4))
      CM1INV=(Tp1up4/(CAP1*(DM1+DM2))+Tp2up4/(CAP2*(DM1+DM2)))/(Tp1up4+Tp2up4)
      FL1=FLLF1+FLRT1
      GOTO09897
09899 CONTINUE
      IF(.NOT.(Km.EQ.Ncnd))GOTO09898
      FLcor1=((R1*Tp1)**4-(R1*Tp2)**4)*Tp1up4/(CAP1*(DM1+DM2)*(Tp1up4+Tp2up4))
      FLcor2=((R1*Tp1)**4-(R1*Tp2)**4)*Tp2up4/((CAP2*(DM1+DM2)+0.d0*RS1)*(Tp1up4+Tp2up4))
      CM1INV=(Tp1up4/(CAP1*(DM1+DM2))+Tp2up4/(CAP2*(DM1+DM2)+0.d0*RS1))/(Tp1up4+Tp2up4)
      FL1=FLcor1+FLcor2
      FLcore=FL1
      CAPcor=CAP2
      GOTO09897
09898 CONTINUE
      FL1=0.D0
09897 CONTINUE
      IF(.NOT.(EVALJA))GOTO09896
      PT2=PT*UPI
      ENGPL2=ENGPL
      ENGT2=ENGT
      PQPL2=PL2*((PPL*UPI)+UU2*MIN(UU2,0.D0)*AQ2)
      QU2=2.D0*QUHALF
      CAPT2=CAPT
      CAPPL2=CAPPL
      PL=PL2*(1.D0+DELTA)
      CHEM=CHEM0(Km+1)
      Tp=Tp2
      CALLURSOS
      PTth=Tp*PT
      IF(.NOT.(Km+1.GT.Ncnd))GOTO09893
      kmhap=km+1
      CALLHAPPA
      DO09890L=1,NFRUS
      hapL2D(L)=hapPAL(L)
      HAAB2D(L)=hapABS(L)
09890 CONTINUE
      KMth=Km+1
      I6R001=5
      GOTO09999
09994 CONTINUE
09893 CONTINUE
      PTPL2=PTth*UPI
      ETPL2=ET*UEI
      PL=PL2
      CHEM=CHEM0(Km+1)
      Tp=Tp2*(1.D0+DELTA)
      CALLURSOS
      PTth=Tp*PT
      IF(.NOT.(Km+1.GT.Ncnd))GOTO09887
      kmhap=km+1
      CALLHAPPA
      DO09884L=1,NFRUS
      hapL2T(L)=hapPAL(L)
      HAAB2T(L)=hapABS(L)
09884 CONTINUE
      I6R001=6
      GOTO09999
09993 CONTINUE
09887 CONTINUE
      PTT2=PTth*UPI
      ETT2=ET*UEI
      Tp=Tp2
09896 CONTINUE
      GOTO09910
09911 CONTINUE
      DM2=DMOUT
      R2=R1+(R1-R0)
      Q2=0.D0
      CM1INV=1.D0/(CAP1*(DM1+DM2)+RS1*CFR)
      PQ2=0.D0
      FC1=0.D0
      CFC1=0.D0
      If(EVALJA)then
      PQPL2=0.D0
      QU2=0.D0
      PT2=0.D0
      CAPT2=0.D0
      CAPPL2=0.D0
      PL2=0.D0
      Tp2=0.D0
      ENGPL2=0.D0
      ENGT2=0.D0
      endif
09910 CONTINUE
      FSAVE(NYDIM+Km)=U1
      RODM=2.D0/(DM1+DM2)
      RSODM=RS1*RODM
      FSAVE(NYDIM+Nzon+Km)=RSODM*(PQ1-PQ2)-AM(Km)/RS1
      DMETIN=1.D0/(DM1*ET1)
      CDMETI=CK1*DMETIN
      CK2ET1=CK2/ET1
      IF(.NOT.(Km.LE.Ncnd))GOTO09881
      HEATRA=(FL0-FL1)*CDMETI
      FLUM=CLUM*FL1
      FLX0=FL0
      FLX=FL1
      FIMP=0.D0
      GOTO09880
09881 CONTINUE
      FSAVE(NYDIM+Nzon+Km)=FSAVE(NYDIM+Nzon+Km)+RODM*(RTphi(Km)*U0*(-min(0.d0,UQ1))/(R1-R0)**NRT-RTphi(Km+1)*U2*(-min(0.d0,UQ2))/(R2
     *-R1)**NRT)
      Tp=Tp1
      PL=PL1
      DO09878L=1,NFRUS
      hapL(L)=hapL1(L)
      HAAB(L)=HAAB1(L)
      tau1(L)=HapL1(L)*(dm1+dm2)/rs1
      tau2(L)=HapL2(L)*(dm1+dm2)/rs1
      hapH(L)=(Tp1up4+Tp2up4)/((Tp2up4/hapL2(L))*(1.d0+1.5d0*(exp(-tau2(L))+tau2(L)-1.d0)*exp(-tau1(L)))+(Tp1up4/hapL1(L))*(1.d0+1.5
     *d0*(exp(-tau1(L))+tau1(L)-1.d0)*exp(-tau2(L))))
      hapH1(L)=hapH(L)
09878 CONTINUE
      FRST=.TRUE.
      calltraneq
      If(Badste)return
      FRST=.FALSE.
      FIMP=0.D0
      FLUM=0.D0
      IF(Km.LT.NZON)THEN
      DO09875L=1,NFRUS
      FIMP=FIMP+FSAVE(NVARS*NZON-Ncnd+KRAD+(NZON-Ncnd)*(L-1)+Km)*WEIGHT(L)*hapH1(L)
      FLUM=FLUM+FSAVE(NVARS*NZON-Ncnd+KRAD+(NZON-Ncnd)*(L-1)+Km)*WEIGHT(L)
      If(L.EQ.LfrDebug.AND.Km.LT.Ncnd+5)WRITE(4,'(A,2I5,1X,6(A,1P,E20.12))')' km, Nperturb =',km,Nperturb,' @FH1=',FSAVE(NVARS*NZON-
     *Ncnd+KRAD+(NZON-Ncnd)*(L-1)+Km)
09875 CONTINUE
      ELSE
      DO09872L=1,NFRUS
      FIMP=FIMP+FSAVE(NVARS*NZON-Ncnd+(NZON-Ncnd)*(L-1)+Km)*HEDD(L)*WEIGHT(L)*hapL1(L)
      FLUM=FLUM+FSAVE(NVARS*NZON-Ncnd+KRAD+(NZON-Ncnd)*(L-1)+Km)*WEIGHT(L)*HEDD(L)
09872 CONTINUE
      ENDIF
      FIMP=CIMP*FIMP
      FLUM=CLUMF*FLUM*RS1
      FLX0=FLX
      FLX=FLUM/CLUM
      FSAVE(NYDIM+Nzon+Km)=FSAVE(NYDIM+Nzon+Km)+FIMP
      DO09869L=1,NFRUS
      FSAVE(NVARS*NZON-Ncnd+(NZON-Ncnd)*(L-1)+Km+NYDIM)=DFJ(L)
      FSAVE(NVARS*NZON-Ncnd+KRAD+(NZON-Ncnd)*(L-1)+Km+NYDIM)=DFH(L)
09869 CONTINUE
      HEATRA=HEATRA/ET1
09880 CONTINUE
      If(abs(KNadap).EQ.5)FSAVE(NYDIM+Nzon+Km)=0.d0
      If(Tp1.GT.1.D+1/UTP)then
      DTPOLD=ENG1*(CK2ET1)+(RS0*U0-RS1*U1)*(TPTQ1*DMETIN)+HEATRA
      else
      DTPOLD=0.D0
      endif
      FSAVE(NYDIM+2*Nzon+Km)=DTPOLD
      if(NSTEP.GE.NDebug)then
      LfrDebug=56
      IF(Km.GT.Ncnd.AND.Km.LT.Ncnd+4)then
      WRITE(4,'(3(A,I10))')' LfrDebug=',LfrDebug,' Nperturb=',Nperturb
      WRITE(4,'(4(A,I10),2(A,L5))')' NSTEP=',NSTEP,' Km=',Km,' Ncnd=',Ncnd,' Nfrus=',Nfrus,'  EVJA=',EVALJA,'  LTH(40)=',LTHICK(40,K
     *m)
      WRITE(4,*)' Tdot=',FSAVE(NYDIM+2*Nzon+Km)
      WRITE(4,'(1X,6(A,1P,E10.2))')'   R1=',R1,'   U1=',U1,'  Tp=',Tp1,'  PL=',PL1,'  Fl0=',fl0,'   FJ1=',FSAVE(NVARS*NZON-Ncnd+(NZO
     *N-Ncnd)*(Lfrdebug-1)+Km),'   FH1=',FSAVE(NVARS*NZON-Ncnd+(NZON-Ncnd)*(Lfrdebug-1)+Km+Krad),'   EDDJ(Km,LfrDebug)=',EDDJ(Km,Lfr
     *Debug),'   EDDJ(Km+1,LfrDebug)=',EDDJ(Km+1,LfrDebug),'   HEDD(LfrDebug)=',HEDD(LfrDebug)
      endif
      endif
      IF(.NOT.(EVALJA))GOTO09866
      callDIFJAC
09866 CONTINUE
      GOTO09941
09940 CONTINUE
      RETURN
      END
