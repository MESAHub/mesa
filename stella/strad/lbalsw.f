C-_TRACE "@wterm' depos(1), depos(20)',depos(1),depos(20),"
C- _TRACE "@wterm' LubvV=',LubvV,' Nfrus=',Nfrus,"
C- _TRACE "@wterm' Natom=',Natom,"
C- _TRACE "@wterm' Irph=',Irph,"
C- _TRACE "@wterm' Lb Hedd(63)=',Hedd(63),"
C- _TRACE "@wterm' Tp, Pl: =',Tp,Pl,"
C: LOSSEN - LOSSES AND GAINS OF ENERGY *
      SUBROUTINELOSSEN
      IMPLICITREAL*8(A-H,O-Z)
C- INPUT  RHO,Ty,NREG,Ry,DM,CK1,2,- ARGUMENTS IN COMMONS
C- OUTPUT  ELTOT,ELVOL,ELSURF,TPSURF IN COM/BAL/
C-      @wterm  "print*,";
C-NVARS - number of independent variables
      PARAMETER(NVARS=3)
C- 3 - NOCONV,4 - CONV
      include '../obj/nfreq_and_mzone.inc'
C-  PARAMETER(Mzon=90); -- h toy model
C- PARAMETER(Mzon=171); -- zoning for W7
C- PARAMETER(Mzon=340); -- zoning for W7+wind
C-PARAMETER(Mzon=86); -- zoning for W7fhh+wind
C-PARAMETER(Mzon=43); -- zoning for W7fhh
C-   PARAMETER (Mzon=521); -- zoning for ntomi Crab model
      PARAMETER(NYDIM=(NVARS+2*NFREQ)*Mzon,MAXDER=4)
C- PARAMETER(NYDIM=2*NFREQ*Mzon,MAXDER=4);
C- Nydim must be  2*(@Npe+@Npc)*Mzon  to use in Feau
C- PARAMETER(NYDIM=2*(@Npe+@Npc)*Mzon,MAXDER=4); -- to use in Feau
      Parameter(Is=5)
C- for test -- chego???
C-   PARAMETER (NZ=1200000); --  for Nfreq=200, Mzon=200
      PARAMETER(NZ=3000000)
C-  for Nfreq=100, Mzon=600
      Parameter(Nstage=28,Natom=15)
      PARAMETER(KOMAX=80)
C-MAX NUMBER OF STEPS IN PRESCRIBED MOMENTS
      LogicalLSYSTEM
C- /.TRUE./; -- for IBM
      Parameter(LSystem=.FALSE.)
      Parameter(Pi=3.1415926535897932d+00,hPlanc=1.0545716280D-27,Cs=2.9979245800D+10,Boltzk=1.3806504000D-16,Avogar=6.0221417900D+2
     *3,AMbrun=1.6605387832D-24,AMelec=9.1093821500D-28,echarg=4.8032042700D-10,CG=6.6742800000D-08,CMS=1.9884000000D+33,RSol=6.9551
     *000000D+10,ULGR=1.4000000000D+01,UPURS=1.0000000000D+07,ULGPU=7.0000000000D+00,ULGEU=1.3000000000D+01,UPC=3.0856776000D+18,UTP
     *=1.0000000000D+05,URHO=1.0000000000D-06,CARAD=7.5657680191D-15,CSIGM=5.6704004778D-05,ERGEV=1.6021764864D-12,GRADeV=1.16045052
     *85D+04,RADC=7.5657680191D-02,CTOMP=4.0062048575D-01,CCAPS=2.6901213726D+01,CCAPZ=9.8964034725D+00)
      IntegerZn(Natom),ZnCo(Natom+1)
      DimensionAZ(Natom)
      Common/AZZn/AZ,Zn,ZnCo
C
C      include for eve and strad
C*
      Common/NiAdap/tday,t_eve,XNifor(Mzon),AMeveNi,KNadap
      LOGICALFRST
      Parameter(Mfreq=130)
C- think about max Nfreq !
      Common/Kmzon/km,kmhap,Jac,FRST
C-   Common/NiAdap/tday,t_eve,XNifor(Mzon),AMeveNi,KNadap; -- must go to commonEve.inc
C- since it enters also bgcon*trf and is forgotten there
      COMMON/STCOM1/t,H,HMIN,HMAX,EPS,N,METH,KFLAG,JSTART
      COMMON/YMAX/YMAX(NYDIM)
      COMMON/YSTIF/Y(NYDIM,MAXDER+1)
      COMMON/HNUSED/HUSED,NQUSED,NFUN,NJAC,NITER,NFAIL
      COMMON/HNT/HNT(7)
C- FOR COSETBGH
      PARAMETER(DELTA=1.d-05)
C- DISTURBANCE CONSTANT standard was 1.d-04
C- cannot be less ~1.d-6 since opacity is saved in single precision
C-   PARAMETER (LICN=4*NZ,LIRN=8*NZ/5); -- FOR F01BRF, M28Y12
      PARAMETER(LICN=4*NZ,LIRN=2*NZ)
C- FOR F01BRF, M28Y12
C- maximum LIRN=2*NZ
      LogicalNEEDBR
C- F01BRF or M28Y12 needed
      COMMON/STJAC/THRMAT,HL,AJAC(LICN),IRN(LIRN),ICN(LICN),WJAC(NYDIM),FSAVE(NYDIM*2),IKEEP(5*NYDIM),IW(8*NYDIM),IDISP(11),NZMOD,NE
     *EDBR
      LOGICALCONV,
C- .TRUE. IF CONVECTION IS INCLUDED
     *CHNCND,
C- .TRUE. IF NCND IS CHANGED
     *SCAT,
C- .TRUE. IF Scattering included (Hapabs^=Happa)
     *SEP
C- .TRUE. IF absorption and scattering are in separate files
C-           (i.e. absorption without expansion effect)
      COMMON/CUTOFF/FLOOR(NVARS+1),Wacc(NVARS+1),FitTau,TauTol,Rvis,CONV,CHNCND,SCAT,SEP
      LogicalLTHICK
      COMMON/THICK/LTHICK(Nfreq,Mzon)
      COMMON/CONVEC/UC(Mzon),YAINV(Mzon)
      COMMON/RAD/EDDJ(Mzon,Nfreq),EDDH(Mzon),HEDD(Nfreq),HEDRAD,CLIGHT,CKRAD,UFREQ,CFLUX,CCL,CLUM,CLUMF,CIMP,NTHICK(NFREQ),NTHNEW(NF
     *REQ),bolM,NCND,KRAD,NFRUS
C- KRAD=(Mzon-NCND)*NFRUS IN STELLA
      LOGICALEDTM
C- Eddington factors time-dependent==.true.
      COMMON/RADOLD/HEDOLD,HINEDD,EDTM
      Common/newedd/EddN(Mzon,Nfreq),HEdN(Nfreq),tfeau
C- for EDTM==T
      Common/oldedd/EddO(Mzon,Nfreq),HEdo(Nfreq),trlx
C- for EDTM==T
      Common/cnlast/Cnlast
      Common/Dhap/DHaphR(Mzon,Nfreq)
      COMMON/BAND/FREQ(NFREQ+1),
C- frequency bins boundaries,
     *
C- frequency in units of h/(k*Tpunit)
     *FREQMN(NFREQ),
C- frequency mean (middle of bin)
     *WEIGHT(130),
C- bandwidth*freq(mean)**3
     *HAPPAL(NFREQ),
C- dimens-less kappa for Ron's absorp.+scatt.
     *HAPABSRON(NFREQ),
C- dimens-less kappa for Ron's absorption
     *HAPABS(NFREQ),
C- dimens-less kappa for S-B absorption
     *DLOGNU(NFREQ)
C- (log step in frequency)**-1
C- PARAMETER(NFRMIN=Nfreq/5); -- MINIMUM NFRUS
      PARAMETER(NFRMIN=Nfreq/2)
C- MINIMUM NFRUS
      IntegerdLfrMax
      Common/observer/wH(Mfreq),cH(Mfreq),zerfr,Hcom(Mfreq),Hobs(Mfreq),freqob(Mfreq),dLfrMax
      Parameter(NP=15+15-1)
C- Total number of tangent rays
      Common/famu/fstatic(0:NP+1,Nfreq),fobs_corr(0:NP+1,Mfreq),fcom(0:NP+1,Mfreq),amustatic(0:NP+1)
      Common/rays/Pray(0:Np+1),fout(0:NP+1,Mfreq),abermu(0:NP+1),NmuNzon
C- fout probably not needed
      COMMON/LIM/Uplim,Haplim
C- to BEGRAD
      COMMON/AMM/DMIN,DM(Mzon),DMOUT,AMINI,AM(Mzon),AMOUT
C- exactly as in VELTEM
      COMMON/Centr/RCE,Nzon
C- central radius & current Number of zones
      Common/InEn/AMHT,EBurst,tBurst,tbeght
C- Mass of Heated Core, Energy & time
      COMMON/RADPUM/AMNI,XMNi,XNi,KmNick
C- MASS OF Ni CORE, Ni ABUNDANCE
      COMMON/RADGAM/FJgam(Mzon,2),toldg,tnewg
C- zero-moment for gamma-radiation
C- 1 for old time (toldg)      2 for present time (tnewg)
      COMMON/RADGlg/FJglog(Mzon,2)
C- log FJgam
      COMMON/CHEM/CHEM0(Mzon),RTphi(0:Mzon+1),EpsUq
      COMMON/REGIME/NREG(Mzon)
      doubleprecisionNRT
      COMMON/AQ/AQ,BQ,DRT,NRT
C- BQ pseudo-viscosity for R-T
C- DRT is the mass distance or optical thickness
C- and used for weight in artificial viscosity
C- NRT is the power in pseudo-viscosity (may be noninteger)
      COMMON/AZNUC/ACARB,ZCARB,ASI,ZSI,ANI,ZNI,QCSI,QSINI
      COMMON/QNRGYE/QNUC,RGASA,YELECT
C-   COMMON/CKN1/CK1,CK2,CFR,CRAP;
      COMMON/CKN1/CK1,CK2,CFR,CRAP,CRAOLD
      LOGICALEVALJA,
C- EVALUATE JACOBIAN
     *OLDJAC,
C- JACOBIAN MAY NEED UPDATING
     *BADSTE
C- BAD STEP: RETURNED BY DIFMAT - TRY SMALLER
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
C- UNITS:
      COMMON/UNSTL/UL,UPRESS,UE,UEPS,UCAP,UTIME,UV,UFLUX,UP
      COMMON/TAIL/KTAIL
      COMMON/UNINV/UPI,UEI
C- 1./UP, 1./UE
      COMMON/UNBSTL/UR,UM,UEPRI,ULGP,ULGE,ULGV,ULGTM,ULGEST,ULGFL,ULGCAP,ULGEPS
      COMMON/CONUR/EIT,DST,BBRCNR(5)
      COMMON/BAL/EL(MAXDER+1),YENTOT(MAXDER+1),ETOT0,ELVOL,ELSURF,ELTOT,TPSURF,HOLDBL,ELOST,EKO,RADBEG
      common/NSTEP/NSTEP,NDebug,MAXER,IOUT,
C- IOUT      LINE PRINT INTERVAL
     *NOUT
C- NOUT      PRINT STEP
      common/CREAD/TAUOLD,NSTMAX,MBATCH,MAXORD
      common/debug/LfrDebug,Nperturb,Kbad
      REAL*8TPMAX(MAXDER+1),TQ(4)
      COMMON/TAU/TAU(Mzon+1),FLUX(Mzon)
      common/tauubvri/tauU(Mzon),tauB(Mzon),tauV(Mzon),tauR(Mzon),tauI(Mzon)
      COMMON/PHOT/XJPH,DMPH,RPH,TPH,PLPH,VPH,CHEMPH,GRVPH,HP,JPH
      PARAMETER(NFUNC=6)
      REAL*4WORK(Mzon+2,NFREQ)
CNFUNC IF NFREQ < NFUNC*
     *,WRK(Mzon,4)
      REAL*8WRKX(Mzon),WORKX(Mzon+2)
      COMMON/STEPD/WRKX,WORKX,TPHOT,TEFF,WORK,WRK,NPHOT,NZM
C       1 - LG(T), 2 - LG(PL), 3 - LG(P), 4 - LG(S)      *
      PARAMETER(TMCRIT=1.D-6,TPNSE=5.D0,EPGROW=0.02D0)
      Common/RUTP/Ry(Mzon),Uy(Mzon),Ty(Mzon),Press(Mzon),Rho(Mzon)
      COMMON/TOO/TOO,KO,KNTO,TO(KOMAX),STO(KOMAX),NTO(KOMAX)
      Parameter(Lcurdm=1000)
C- Dimension of the Tcurv array
      RealTcurv
      IntegerNFRUSED
C- an integer array which store exact number of used freqs
      REAL*8Flsave
C- remove this to read old flx files!
      Common/Curve/tcurv(8,Lcurdm),Depos(Lcurdm),Flsave(MFREQ+1,Lcurdm),NFRUSED(Lcurdm),Lsaved
      LOGICALBEGRUN
      Common/BEGR/BEGRUN
      CHARACTER*80Model,Sumprf,Sumcur,Depfile,Flxfile
      COMMON/Files/Model,Sumprf,Sumcur,Depfile,Flxfile
      CHARACTER*1app
C- dummy var for constructing of Opafile
      LogicalGivdtl
      Common/ABGrap/NSTA,NSTB,TcurA,TcurB,Givdtl
C-No. of steps & Time in days
      REAL*8MBOL,MU,MB,MV,MR,MI,MBOL1
      COMMON/COLOR/MBOL,MU,MB,MV,MR,MI,UMB,BMV,MBOL1,LubvU,LubvB,LubvV,LubvR,LubvI,Lyman
      COMMON/DETAIL/QRTarr(Mzon),UUarr(Mzon),ArrLum(Mzon),Acc(Mzon)
      Common/XYZ/XA,YA,URM
C-_TRACE "@wterm' LubvV=',LubvV,' Nfrus=',Nfrus,"
      ELVOL=0.D0
      DO09999I=1,NZON
      PL=RHO(I)
      Tp=Ty(I)
      CALLVOLEN(I)
      If(NREG(I).EQ.1)Then
C: FIND YDOT & HEATING ENGNUC *
      YCARB=1.D0/Y((NVARS-1)*NZON+I,1)
C-YCINV
      CALLBURNC
      ENGNUC=QNUC*YDOT*YCARB**2
      else
      ENGNUC=0.D0
      endif
      ENG=ENG+ENGNUC
      ELVOL=ELVOL+ENG*DM(I)
09999 CONTINUE
      ELVOL=CK2*ELVOL
C-  @wterm ' eLvol=',eLvol;
      Flum=0.D0
      DO09996L=1,NFRUS
      Flum=Flum+Y(NVARS*NZON+(NZON-NCND)*L-1+KRAD,1)*WEIGHT(L)
C- FOR I=NZON-1
     *
09996 CONTINUE
C-      @wterm 'Flum 1=',Flum;
      Flum=(CLUMF*1.D-50)*Flum*Ry(NZON-1)**2
C-      @wterm 'Flum 2=',Flum;
      ELSURF=Flum/(UEPRI/(UTIME*CRAP))
C-      @wterm 'Flum=',Flum,'   Ls=',Elsurf;
C- CALL URSOS; CALL OPACIT;
C- ELSURF=(Ry(NZON))**2/(CAPPA*(DM(NZON)+DMOUT)+CFR*Ry(NZON)**2);
C- TPSURF=Ty(NZON)*(ELSURF*CFR)**.25D0;
C- ELSURF=CK1*ELSURF*(Ry(NZON)*Ty(NZON))**2*Ty(NZON)**2;
C- TEFF=TPSURF*1.189207D0*UTP;
      ELTOT=ELVOL-ELSURF
      RETURN
      END
C: PRIBAL - PRINT MODEL & BALANCES *
C-  DEBUG SUBCHK,UNIT(6);ENDDEBUG;
      SUBROUTINEPRIBAL(IPOUT
Cline print interval*
     *)
      IMPLICITREAL*8(A-H,O-Z)
C-      @wterm  "print*,";
C-NVARS - number of independent variables
      PARAMETER(NVARS=3)
C- 3 - NOCONV,4 - CONV
      include '../obj/nfreq_and_mzone.inc'
C-  PARAMETER(Mzon=90); -- h toy model
C- PARAMETER(Mzon=171); -- zoning for W7
C- PARAMETER(Mzon=340); -- zoning for W7+wind
C-PARAMETER(Mzon=86); -- zoning for W7fhh+wind
C-PARAMETER(Mzon=43); -- zoning for W7fhh
C-   PARAMETER (Mzon=521); -- zoning for ntomi Crab model
      PARAMETER(NYDIM=(NVARS+2*NFREQ)*Mzon,MAXDER=4)
C- PARAMETER(NYDIM=2*NFREQ*Mzon,MAXDER=4);
C- Nydim must be  2*(@Npe+@Npc)*Mzon  to use in Feau
C- PARAMETER(NYDIM=2*(@Npe+@Npc)*Mzon,MAXDER=4); -- to use in Feau
      Parameter(Is=5)
C- for test -- chego???
C-   PARAMETER (NZ=1200000); --  for Nfreq=200, Mzon=200
      PARAMETER(NZ=3000000)
C-  for Nfreq=100, Mzon=600
      Parameter(Nstage=28,Natom=15)
      PARAMETER(KOMAX=80)
C-MAX NUMBER OF STEPS IN PRESCRIBED MOMENTS
      LogicalLSYSTEM
C- /.TRUE./; -- for IBM
      Parameter(LSystem=.FALSE.)
      Parameter(Pi=3.1415926535897932d+00,hPlanc=1.0545716280D-27,Cs=2.9979245800D+10,Boltzk=1.3806504000D-16,Avogar=6.0221417900D+2
     *3,AMbrun=1.6605387832D-24,AMelec=9.1093821500D-28,echarg=4.8032042700D-10,CG=6.6742800000D-08,CMS=1.9884000000D+33,RSol=6.9551
     *000000D+10,ULGR=1.4000000000D+01,UPURS=1.0000000000D+07,ULGPU=7.0000000000D+00,ULGEU=1.3000000000D+01,UPC=3.0856776000D+18,UTP
     *=1.0000000000D+05,URHO=1.0000000000D-06,CARAD=7.5657680191D-15,CSIGM=5.6704004778D-05,ERGEV=1.6021764864D-12,GRADeV=1.16045052
     *85D+04,RADC=7.5657680191D-02,CTOMP=4.0062048575D-01,CCAPS=2.6901213726D+01,CCAPZ=9.8964034725D+00)
      IntegerZn(Natom),ZnCo(Natom+1)
      DimensionAZ(Natom)
      Common/AZZn/AZ,Zn,ZnCo
C
C      include for eve and strad
C*
      Common/NiAdap/tday,t_eve,XNifor(Mzon),AMeveNi,KNadap
      LOGICALFRST
      Parameter(Mfreq=130)
C- think about max Nfreq !
      Common/Kmzon/km,kmhap,Jac,FRST
C-   Common/NiAdap/tday,t_eve,XNifor(Mzon),AMeveNi,KNadap; -- must go to commonEve.inc
C- since it enters also bgcon*trf and is forgotten there
      COMMON/STCOM1/t,H,HMIN,HMAX,EPS,N,METH,KFLAG,JSTART
      COMMON/YMAX/YMAX(NYDIM)
      COMMON/YSTIF/Y(NYDIM,MAXDER+1)
      COMMON/HNUSED/HUSED,NQUSED,NFUN,NJAC,NITER,NFAIL
      COMMON/HNT/HNT(7)
C- FOR COSETBGH
      PARAMETER(DELTA=1.d-05)
C- DISTURBANCE CONSTANT standard was 1.d-04
C- cannot be less ~1.d-6 since opacity is saved in single precision
C-   PARAMETER (LICN=4*NZ,LIRN=8*NZ/5); -- FOR F01BRF, M28Y12
      PARAMETER(LICN=4*NZ,LIRN=2*NZ)
C- FOR F01BRF, M28Y12
C- maximum LIRN=2*NZ
      LogicalNEEDBR
C- F01BRF or M28Y12 needed
      COMMON/STJAC/THRMAT,HL,AJAC(LICN),IRN(LIRN),ICN(LICN),WJAC(NYDIM),FSAVE(NYDIM*2),IKEEP(5*NYDIM),IW(8*NYDIM),IDISP(11),NZMOD,NE
     *EDBR
      LOGICALCONV,
C- .TRUE. IF CONVECTION IS INCLUDED
     *CHNCND,
C- .TRUE. IF NCND IS CHANGED
     *SCAT,
C- .TRUE. IF Scattering included (Hapabs^=Happa)
     *SEP
C- .TRUE. IF absorption and scattering are in separate files
C-           (i.e. absorption without expansion effect)
      COMMON/CUTOFF/FLOOR(NVARS+1),Wacc(NVARS+1),FitTau,TauTol,Rvis,CONV,CHNCND,SCAT,SEP
      LogicalLTHICK
      COMMON/THICK/LTHICK(Nfreq,Mzon)
      COMMON/CONVEC/UC(Mzon),YAINV(Mzon)
      COMMON/RAD/EDDJ(Mzon,Nfreq),EDDH(Mzon),HEDD(Nfreq),HEDRAD,CLIGHT,CKRAD,UFREQ,CFLUX,CCL,CLUM,CLUMF,CIMP,NTHICK(NFREQ),NTHNEW(NF
     *REQ),bolM,NCND,KRAD,NFRUS
C- KRAD=(Mzon-NCND)*NFRUS IN STELLA
      LOGICALEDTM
C- Eddington factors time-dependent==.true.
      COMMON/RADOLD/HEDOLD,HINEDD,EDTM
      Common/newedd/EddN(Mzon,Nfreq),HEdN(Nfreq),tfeau
C- for EDTM==T
      Common/oldedd/EddO(Mzon,Nfreq),HEdo(Nfreq),trlx
C- for EDTM==T
      Common/cnlast/Cnlast
      Common/Dhap/DHaphR(Mzon,Nfreq)
      COMMON/BAND/FREQ(NFREQ+1),
C- frequency bins boundaries,
     *
C- frequency in units of h/(k*Tpunit)
     *FREQMN(NFREQ),
C- frequency mean (middle of bin)
     *WEIGHT(130),
C- bandwidth*freq(mean)**3
     *HAPPAL(NFREQ),
C- dimens-less kappa for Ron's absorp.+scatt.
     *HAPABSRON(NFREQ),
C- dimens-less kappa for Ron's absorption
     *HAPABS(NFREQ),
C- dimens-less kappa for S-B absorption
     *DLOGNU(NFREQ)
C- (log step in frequency)**-1
C- PARAMETER(NFRMIN=Nfreq/5); -- MINIMUM NFRUS
      PARAMETER(NFRMIN=Nfreq/2)
C- MINIMUM NFRUS
      IntegerdLfrMax
      Common/observer/wH(Mfreq),cH(Mfreq),zerfr,Hcom(Mfreq),Hobs(Mfreq),freqob(Mfreq),dLfrMax
      Parameter(NP=15+15-1)
C- Total number of tangent rays
      Common/famu/fstatic(0:NP+1,Nfreq),fobs_corr(0:NP+1,Mfreq),fcom(0:NP+1,Mfreq),amustatic(0:NP+1)
      Common/rays/Pray(0:Np+1),fout(0:NP+1,Mfreq),abermu(0:NP+1),NmuNzon
C- fout probably not needed
      COMMON/LIM/Uplim,Haplim
C- to BEGRAD
      COMMON/AMM/DMIN,DM(Mzon),DMOUT,AMINI,AM(Mzon),AMOUT
C- exactly as in VELTEM
      COMMON/Centr/RCE,Nzon
C- central radius & current Number of zones
      Common/InEn/AMHT,EBurst,tBurst,tbeght
C- Mass of Heated Core, Energy & time
      COMMON/RADPUM/AMNI,XMNi,XNi,KmNick
C- MASS OF Ni CORE, Ni ABUNDANCE
      COMMON/RADGAM/FJgam(Mzon,2),toldg,tnewg
C- zero-moment for gamma-radiation
C- 1 for old time (toldg)      2 for present time (tnewg)
      COMMON/RADGlg/FJglog(Mzon,2)
C- log FJgam
      COMMON/CHEM/CHEM0(Mzon),RTphi(0:Mzon+1),EpsUq
      COMMON/REGIME/NREG(Mzon)
      doubleprecisionNRT
      COMMON/AQ/AQ,BQ,DRT,NRT
C- BQ pseudo-viscosity for R-T
C- DRT is the mass distance or optical thickness
C- and used for weight in artificial viscosity
C- NRT is the power in pseudo-viscosity (may be noninteger)
      COMMON/AZNUC/ACARB,ZCARB,ASI,ZSI,ANI,ZNI,QCSI,QSINI
      COMMON/QNRGYE/QNUC,RGASA,YELECT
C-   COMMON/CKN1/CK1,CK2,CFR,CRAP;
      COMMON/CKN1/CK1,CK2,CFR,CRAP,CRAOLD
      LOGICALEVALJA,
C- EVALUATE JACOBIAN
     *OLDJAC,
C- JACOBIAN MAY NEED UPDATING
     *BADSTE
C- BAD STEP: RETURNED BY DIFMAT - TRY SMALLER
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
C- UNITS:
      COMMON/UNSTL/UL,UPRESS,UE,UEPS,UCAP,UTIME,UV,UFLUX,UP
      COMMON/TAIL/KTAIL
      COMMON/UNINV/UPI,UEI
C- 1./UP, 1./UE
      COMMON/UNBSTL/UR,UM,UEPRI,ULGP,ULGE,ULGV,ULGTM,ULGEST,ULGFL,ULGCAP,ULGEPS
      COMMON/CONUR/EIT,DST,BBRCNR(5)
      COMMON/BAL/EL(MAXDER+1),YENTOT(MAXDER+1),ETOT0,ELVOL,ELSURF,ELTOT,TPSURF,HOLDBL,ELOST,EKO,RADBEG
      common/NSTEP/NSTEP,NDebug,MAXER,IOUT,
C- IOUT      LINE PRINT INTERVAL
     *NOUT
C- NOUT      PRINT STEP
      common/CREAD/TAUOLD,NSTMAX,MBATCH,MAXORD
      common/debug/LfrDebug,Nperturb,Kbad
      REAL*8TPMAX(MAXDER+1),TQ(4)
      COMMON/TAU/TAU(Mzon+1),FLUX(Mzon)
      common/tauubvri/tauU(Mzon),tauB(Mzon),tauV(Mzon),tauR(Mzon),tauI(Mzon)
      COMMON/PHOT/XJPH,DMPH,RPH,TPH,PLPH,VPH,CHEMPH,GRVPH,HP,JPH
      PARAMETER(NFUNC=6)
      REAL*4WORK(Mzon+2,NFREQ)
CNFUNC IF NFREQ < NFUNC*
     *,WRK(Mzon,4)
      REAL*8WRKX(Mzon),WORKX(Mzon+2)
      COMMON/STEPD/WRKX,WORKX,TPHOT,TEFF,WORK,WRK,NPHOT,NZM
C       1 - LG(T), 2 - LG(PL), 3 - LG(P), 4 - LG(S)      *
      PARAMETER(TMCRIT=1.D-6,TPNSE=5.D0,EPGROW=0.02D0)
      Common/RUTP/Ry(Mzon),Uy(Mzon),Ty(Mzon),Press(Mzon),Rho(Mzon)
      COMMON/TOO/TOO,KO,KNTO,TO(KOMAX),STO(KOMAX),NTO(KOMAX)
      Parameter(Lcurdm=1000)
C- Dimension of the Tcurv array
      RealTcurv
      IntegerNFRUSED
C- an integer array which store exact number of used freqs
      REAL*8Flsave
C- remove this to read old flx files!
      Common/Curve/tcurv(8,Lcurdm),Depos(Lcurdm),Flsave(MFREQ+1,Lcurdm),NFRUSED(Lcurdm),Lsaved
      LOGICALBEGRUN
      Common/BEGR/BEGRUN
      CHARACTER*80Model,Sumprf,Sumcur,Depfile,Flxfile
      COMMON/Files/Model,Sumprf,Sumcur,Depfile,Flxfile
      CHARACTER*1app
C- dummy var for constructing of Opafile
      LogicalGivdtl
      Common/ABGrap/NSTA,NSTB,TcurA,TcurB,Givdtl
C-No. of steps & Time in days
      REAL*8MBOL,MU,MB,MV,MR,MI,MBOL1
      COMMON/COLOR/MBOL,MU,MB,MV,MR,MI,UMB,BMV,MBOL1,LubvU,LubvB,LubvV,LubvR,LubvI,Lyman
      COMMON/DETAIL/QRTarr(Mzon),UUarr(Mzon),ArrLum(Mzon),Acc(Mzon)
      Common/XYZ/XA,YA,URM
C:  PARAMETERS & COMMONS *
C-  COMMON/IONI/XH,XHEI,XHEII,XP,XCR;
C-    Dimension Xion(Nstage,Natom);
C-    Dimension tauU(Mzon),tauB(Mzon),tauV(Mzon),tauR(Mzon),tauI(Mzon);
      Dimensiontauzon(Nfreq)
      LogicalLOW,HIGH
      COMMON/NIT/Xelow,Xion(0:Nstage,Natom),Nit,LOW,HIGH
      REAL*8BL(NFREQ),FJ(NFREQ),FH(NFREQ)
      REAL*8FJL(NFREQ)
C-Dimension workm(Mzon),Xlow(2),Xup(2);
      COMMON/CONV/ELMX(Mzon),FLCON(Mzon)
C-    Common/observer/wH(Nfreq),cH(Nfreq),zerfr;
C-    Parameter(NP=@Npe+@Npc-1);   -- Total number of tangent rays
C-    Common/rays/Pray(0:NP+1),fout(0:NP+1,Nfreq);
      common/phofit/Tphfit,Rphfit
      common/corrttt/ttt(nfreq)
      Common/out/tretard,tout,Yout(NYdim)
C- for interpolation
      CHARACTER*120STR2
      character*75symb
      Parameter(RADC3=RADC/3.D0)
      Parameter(FJnois=1.d-30)
      Parameter(tolIni=1.d-2)
C- how large is tolerable initial condition disbalance
C-  Real*4  RTpsi(Mzon),RTchi(Mzon);
C-  Equivalence (RTpsi,WRK),(RTchi,WRK(1,3));
      Logicalfullp,fullgp
C: Vars for GPLOT    *
      CHARACTERTHD(4)*60,THX(4)*30,THY(4)*30
      DIMENSIONMX(6),XMN(6),YMX(6),YMN(6),ITYP(6)
      DATATHD/' Tp  Rho  P  S       ','           Spectra   ','               FIG. 3','               FIG. 4'/
      DATATHX/'          Zone Number','        Freq. group Number','     Lagr. mass (Solar units)','        Radius (1.E14 cm)'/
      DATATHY/'              LG(Tp D P) ','              log( J) ','               LOG(RPH/RSUN)','                   -MV'/
      DATANHD/40/,NHX/30/,NHY/30/
C-    DATA ITYP/4*1/;
C-    DATA THDB/' '/,THXB/' '/,THYB/' '/;
      DATAHTEXT/0.35D0/
      BLACK(Lbl,Tpbl)=(exp(-(FREQMN(Lbl)/Tpbl)))/(1.d0-(exp(-(FREQMN(Lbl)/Tpbl))))
      BLACKD(Lbl,Tpbl)=(FREQMN(Lbl)/Tpbl)*(exp(-(FREQMN(Lbl)/Tpbl)))/(1.d0-(exp(-(FREQMN(Lbl)/Tpbl))))**2
      Theta(Z)=0.5D0*(SIGN(1.D0,Z)+1.D0)
C-  DQRGAS=CTOMP*(UR/Utime)**2/(Boltzk*UTp*Avogar)*UR*URHO;;
      Fullgp=Ipout.EQ.0
C- full print with graphics
      Fullp=Ipout.LE.0
C- full print
      Lpout=max(abs(ipout),1)
C: BODY OF PRIBAL *
C-    @wterm ' do not forget ubvnew !!!';
C-  pause;
C: print head_line : stiff & head of zone output *
      WRITE(4,'(''%H:'')')
      WRITE(4,'(''   NSTEP   KFLAG  JSTART  NQUSED    NFUN    NJAC'',              ''       N   NFrus    Ncnd'',              ''   N
     *ITER   NFAIL   NZMOD'')')
      WRITE(4,'(12I8)')NSTEP,KFLAG,JSTART,NQUSED,NFUN,NJAC,N,NFrus,Ncnd,NITER,NFAIL,NZM
      WRITE(4,'(''  OBS.TIME='',F15.5,'' D   PROPER T='',1P,E15.8,'' S'',              ''  STEP USED='',E11.5,              ''  STEP
     * TRIED='',E11.5,'' S'')')(tout-tretard)*(UTIME*CRAP)/86400.D0,tout*(UTIME*CRAP),HUSED*(UTIME*CRAP),H*(UTIME*CRAP)
      If(CONV)Then
      WRITE(4,95)ULGR,ULGR-6,LOG10(UTP),LOG10(URHO),ULGP-LOG10(UP)
      else
      WRITE(4,94)ULGR,ULGR-6,LOG10(UTP),LOG10(URHO),ULGP-LOG10(UP)
      endif
      KPRINT=Lpout
C- step of print
C: prepare
C           balances, arrays and line for zone print & graphics*
C-    @wterm ' Ty(Nzon)=',Ty(Nzon),'  rho=',rho(Nzon);
      KENTR=1
      ETERM=0.D0
      EKIN=0.D0
      EGRAV=0.D0
      EVIR=0.D0
      VISVIR=0.D0
      ACVIR=0.D0
      ERADD=0.D0
      If(NCND.GT.1)then
      RNcnd2=Ry(min(NCND+1,Mzon))
      RNCND=Ry(min(NCND,Mzon-1))
      else
      RNcnd2=0.D0
      RNCND=0.D0
      endif
C: Find Teff,U,B,V & store HAPPAL in WORK(.,.), cappa in WRKX()
C    <*Tpcomp: find in Tpcomp the highest Tp(i) for i>ncnd
C              for use in approximate treatment of Compton *
C: FIND PARAMETERS OF PHOTOSPHERE ( %MA ) *
      Km=NZON
      tau(Km)=0.D0
      tauTp=0.D0
      IRph=0
09993 CONTINUE
      PL=RHO(Km)
C-         If(Ty(km)<0) @wterm ' Ty(km)=',Ty(km),'  in km=',km;;
C-         Tp=abs(Ty(Km));
      Tp=Ty(Km)
C-         CHEM=CHEM0(Km);
      DO09990ii=1,Natom
      Yat(ii)=YABUN(ii,Km)
09990 CONTINUE
      If(Km.GT.Ncnd)then
      RADP=.False.
      else
      RADP=.True.
      endif
      CALLURSOS
      CHEM0(Km)=CHEM
      kmhap=Km
      CALLHAPPA
      CallOpacit
      Press(Km)=P
      WRK(Km,1)=Xion(0,1)
C- HI
      WRK(Km,2)=Egas
      WRKX(Km)=Cappa
      DO09987L=1,NFRUS
      WORK(Km,L)=SNGL(HAPPAL(L))
09987 CONTINUE
C-    Km=Km-1; -- wrong place???
      If(Km.GT.1)Then
      tau(Km-1)=tau(Km)+(Ry(Km)-Ry(Km-1))*HAPPAL(LubvB)*PL
      tauTp=0.5D0*(tau(Km)+tau(min(Km+1,Nzon)))
      Endif
      If(tau(Km).GE.0.64D0.AND.IRph.EQ.0)IRph=Km
C-  @wterm ' km pl @freqtau hap ',km,pl,@freqtau,happal(@freqtau);
      Km=Km-1
C- new place
      IF(.NOT.(Km.EQ.0.OR.(tau(Km+1).GE.0.64D0)))GOTO09993
C-   _UNTIL Km==0 ! (tau(Km) >= @tauPH ! tauTp >= @tauPH);
      IRph=min(IRph,Nzon-2)
      JPH=IRph
      If(NCND.GT.0)Then
C: find Rph, Vph & Tph interpolating in tau *
      RPH=Ry(IRph+1)+(Ry(IRph)-Ry(IRph+1))/(tau(IRph)-tau(IRph+1))*(0.64D0-tau(IRph+1))
      RPH=max(RPH,Ry(1))
      VPH=Uy(IRph+1)+(Uy(IRph)-Uy(IRph+1))/(tau(IRph)-tau(IRph+1))*(0.64D0-tau(IRph+1))
      XPH=YABUN(1,IRph+1)+(YABUN(1,IRph)-YABUN(1,IRph+1))/(tau(IRph)-tau(IRph+1))*(0.64D0-tau(IRph+1))
C-   0.5D0*(tau(IRph)+tau(IRph+1)) -- tau IN THE MIDDLE OF ZONE IRph+1
      TPH4=Ty(IRph+2)**4+(Ty(IRph+1)**4-Ty(IRph+2)**4)/(0.5D0*(tau(IRph)-tau(IRph+2)))*(0.64D0-0.5D0*(tau(IRph+1)+tau(IRph+2)))
      TPH=ABS(TPH4)**.25D0
C-   TPH=SIGN(TPH,TPH4);
      else
C: fitting Tph for transparent core, then Rph & Vph*
      TPH=Ty(1)
C- CHANGE!!! FOR FTEFF
      RPH=SQRT(ABS(ELSURF)/(TPH**4*CSIGM*UTp**4*UTIME**3/(URHO*UR**3)))
C-   @wterm' Lbalrad RPH=',rph;
      I=1
09984 IF(.NOT.(I.LT.NZON-1.AND.Ry(I+1).LT.RPH))GOTO09983
      I=I+1
      GOTO09984
09983 CONTINUE
C- Ry(I+1)>=RPH
      VPH=Uy(I+1)+(Uy(I)-Uy(I+1))/(Ry(I)-Ry(I+1))*(RPH-Ry(I+1))
      endif
      FLOUT=CFR*ELSURF/CK1*2.D0
      TPH=TPH*UTp
      TEFF=TPH
      CALLUBV
      TPH=TPH/UTp
      Km=JPH
      Km=Nzon
C- test
      tauU(Km)=0.d0
      tauB(Km)=0.d0
      tauV(Km)=0.d0
      tauR(Km)=0.d0
      tauI(Km)=0.d0
      open(40,file='hapint.tst')
      write(40,*)'nstep=',nstep,'   tday=',tout*UTIME/86400.D0
C- tday
09981 IF(.NOT.(Km.GT.0))GOTO09980
C-    _WHILE Km>NCND _DO
C-      @wterm ' rho km ', km, rho(km),rho(km-1),rho(km-2),dm(km),dm(km-1);
      PL=RHO(Km)
      Tp=Ty(Km)
C- CHEM=CHEM0(Km);
      DO09978ii=1,Natom
      Yat(ii)=YABUN(ii,Km)
09978 CONTINUE
      RADP=.FALSE.
C-      @wterm' before URSOS';
C-      @wterm' km rho pl: ',Km,Rho(KM), Pl;
      CALLURSOS
      CHEM0(Km)=CHEM
      kmhap=Km
      CALLHAPPA
      write(40,'(1x,i4,1p,2e10.2,0p,(20f6.2))')km,pl*urho,Tp*uTp,(log10(max(1.d-50,hapabs(ll)/(ur*urho))),ll=1,nfrus)
      If(Scat)write(40,'('' scat'',20x,20f6.2)')(log10(max(1.d-50,(happal(ll)-hapabs(ll))/(ur*urho))),ll=1,nfrus)
      CallOpacit
      Press(Km)=P
C-         WRK(Km,1)=Xion(1,1); -- HII
      WRK(Km,1)=Xion(0,1)
C- HI
      WRK(Km,2)=Egas
      WRKX(Km)=Cappa
      DO09975L=1,NFRUS
      WORK(Km,L)=SNGL(HAPPAL(L))
09975 CONTINUE
      If(Km.GT.1)then
      tau(Km-1)=tau(Km)+(Ry(Km)-Ry(Km-1))*HAPPAL(LubvB)*PL
      tauU(Km-1)=tauU(Km)+(Ry(Km)-Ry(Km-1))*HAPPAL(LubvU)*PL
      tauB(Km-1)=tauB(Km)+(Ry(Km)-Ry(Km-1))*HAPPAL(LubvB)*PL
      tauV(Km-1)=tauV(Km)+(Ry(Km)-Ry(Km-1))*HAPPAL(LubvV)*PL
      tauR(Km-1)=tauR(Km)+(Ry(Km)-Ry(Km-1))*HAPPAL(LubvR)*PL
      tauI(Km-1)=tauI(Km)+(Ry(Km)-Ry(Km-1))*HAPPAL(LubvI)*PL
      endif
      Km=Km-1
      GOTO09981
09980 CONTINUE
      close(40)
C-  Ycomp=(Boltzk*UTp*Tpcomp)/(Amelec*Cs**2); -- Comptonization par.
C-  WRITE(@Wres,*)' Tpcomp   Ycomp:  ',Tpcomp,Ycomp;
C: put  pressure into  Press() for Km<=Ncnd
C                   Cappa    into  WRKX()  for Km<=Ncnd
C                   UU       into  UUarr()
C              put  QRT      into  QRTarr()        *
      DO09972km=1,Nzon
      WORKX(Km)=FLOAT(Km)
C- for BARKUK
      If(Km.GT.1)Then
      U0=Uy(Km-1)
      R0=Ry(Km-1)
      else
      U0=0.D0
      R0=Rce
      endif
      UUarr(Km)=Uy(Km)-U0*(R0/Ry(Km))**2
      If(Km.LE.Ncnd)then
C- P for other Km already in %PA_Pre.T
      RadP=.true.
C-   @wterm ' km Rho(km) dm ',km,Rho(km),dm(km),km+1,rho(km+1),dm(km+1);
      Pl=Rho(Km)
      Tp=Ty(Km)
      DO09969ii=1,Natom
      Yat(ii)=YABUN(ii,Km)
09969 CONTINUE
      CallURSOS
      kmhap=Km
      CallOpacit
      Press(Km)=P
      WRK(Km,1)=Xion(0,1)
      WRK(Km,2)=Egas
      WRKX(Km)=Cappa
C-    QRTarr(Km)=0.;   -- now non-zero everywehere
C- else;
C-  QRTarr(Km)= RTphi(Km)*(-min(0.D0,UUarr(Km)))/(Ry(Km)-R0)**NRT;
C- to multiply by U0 or U2 for Qcold
C- to multiply by abs(U1)*RS1   for Q hot
      endif
      QRTarr(Km)=RTphi(Km)*(-min(0.D0,UUarr(Km)-EpsUq))/(Ry(Km)-R0)**NRT
C-
C If(Km>Ncnd & Km>2 & Km<NZON)then;  -- crudely:
C   -- QRT1=QRT(Km);
C   hot   QRT1=RTphi(Km)*abs(Uy(Km))*(-min(0.D0,UU1))/(Ry(Km+1)-Ry(Km))**NRT;
C   cold  QRT1=RTphi(Km)*Uy(Km)*(-min(0.D0,UU1))/(Ry(Km+1)-Ry(Km))**NRT;
C   -- QRT1=RTphi(Km)*((Ry(Km)-Ry(Km-1))/(Ry(Km+1)-Ry(Km))-1.);
C   -- QRT1=RTpsi(Km)*( 1./(Ry(Km+1)-Ry(Km)) - 1./(Ry(Km)-Ry(Km-1)))*
C     --          (Ry(Km)-Ry(Km-1))*Theta(RTchi(Km));
C   else;
C      QRT1=0.;
C   endif;
C*
      If(Km.EQ.1)Prce=P*UPI
      If(Km.EQ.ncnd+1)then
      Pcnd2=P*UPI
      Tpcnd2=Tp
      endif
09972 CONTINUE
      CAP2=WRKX(1)
C-  CM2 =WRKX(1)*DM(1);
      DO09966Km=1,NZON
      If(Km.LT.NZON)Then
      DM2=DM(Km+1)
      else
      DM2=DMOUT
      endif
      RX=(DM(Km)+DM2)/2.D0
C: find ENG, EGAS, ERAD, Q1, Accrad & put
C                   udot       into  Acc()
C                   luminosity into  ArrLum()         *
CBurn*     *
      If(Km.GT.1)Then
      U0=Uy(Km-1)
      else
      U0=0.D0
      endif
      If(Km.LT.NZON)Then
      PQ2=Press(Km+1)*UPI+AQ*Rho(Km+1)*UUarr(Km+1)*(-min(0.D0,UUarr(Km+1)))
      U2=Uy(Km+1)
      dm2=dm(Km+1)
      else
      PQ2=0.
      U2=0.
      dm2=dmout
      endif
C- for enhanced Q1 it was used AQ1 instead of AQ:
C-      If(Km>NCND)Then;
C-        WXX=DRT*(AM(Km)-AM(max(NCND,1)))**2;
C-        AQ1=AQ+BQ*WXX/(1.D0+WXX);
C-      else;
C-        AQ1=AQ;
C-      endif;
      Q1=AQ*Rho(Km)*UUarr(Km)*min(0.D0,UUarr(Km))
      If(Km.EQ.1)Qrce=Q1
      PQ1=Press(Km)*UPI+Q1
      Pl=rho(Km)
      Tp=Ty(Km)
      CALLVOLEN(Km)
      RSODM=Ry(Km)**2/RX
      Accel=
C- dU/dt
     *RSODM*(PQ1-PQ2)-AM(Km)/Ry(Km)**2
C- dU/dt correction
     *+(QRTarr(Km)*U0-QRTarr(min(Km+1,Nzon))*U2)/RX
C- Qcold
C-+RSODM*( QRTarr(Km)*abs(Uy(Km))-QRTarr(min(Km+1,Nzon))*abs(U2));
C-hot
      If(Km.EQ.Ncnd)Accel=Accel-RSODM*RadC3*UPI*Tpcnd2**4
      QRTarr(Km)=QRTarr(Km)*Uy(Km)
C- crude, only for output (cold)
C- QRTarr(Km)=QRTarr(Km)*abs(Uy(Km)) ; -- hot
C- Q1=Q1+QRTarr(Km);
C:  calculate output flux in Flux(Km) *
C-      workm(Km)=Am(Km)*Um;
C-      CM1=CM2;
      Tp1=Ty(Km)
      CAP1=CAP2
      If(Km.LT.NZON)Then
C-         CM2=CAPPA*DM(Km+1); CMINV=1.D0/(CM1+CM2);
      CAP2=WRKX(Km+1)
      Tp2=Ty(Km+1)
      ELSE
C-         CM2=CM1*(DMOUT/DM(NZON));
      Tp2=0.D0
C- Cap2 the same as previous zone
C-         CMINV=1.D0/(CM1+CM2+CFR*Ry(NZON)**2);
      ENDIF
C-        print *,Tp,Pl,CAP1,CAP2;
C-      _SELECT
C-        _ Km<=NCND
C-             [
      FLLF1=((Ry(Km)*Tp1)**4-(Ry(Km)*Tp2)**4)*Tp1**4/(CAP1*(DM(Km)+DM2)*(Tp1**4+Tp2**4))
      FLRT1=((Ry(Km)*Tp1)**4-(Ry(Km)*Tp2)**4)*Tp2**4/(CAP2*(DM(Km)+DM2)*(Tp1**4+Tp2**4))
      FL1=FLLF1+FLRT1
C- ]
C-             [ FL1=((Ry(Km)*Tp1)**4-(Ry(Km)*Tp2)**4)*CMINV ]
C        _ Km==NCND [ SUMJ=0.D0;
C                 _DO L=1,NFRUS; SUMJ=SUMJ+@FJ2*WEIGHT(L)_OD;
C                 SUMJ=(1.5D1/PI**4)*SUMJ;
C                 FLCOR1=((Ry(Km)*Tp1)**4-SUMJ*Ry(Km)**4)*Tp1**4
C                        /(CAP1*(DM(Km)+DM(Km+1))*(Tp1**4+Tp2**4));
C                 FLCOR2=((Ry(Km)*Tp1)**4-SUMJ*Ry(Km)**4)*Tp2**4
C                    /((CAP2*(DM(Km)+DM(Km+1))+@tauLIM*Ry(Km)**2)
C                    *(Tp1**4+Tp2**4));
C                 FL1=FLCOR1+FLCOR2 ]
C          _ Km==NCND+1
C                [ SUMJ=0.D0;
C                  _DO L=1,NFRUS; SUMJ=SUMJ+@FJ1*WEIGHT(L)_OD;
C                  SUMJ=(1.5D1/PI**4)*SUMJ;
C                  --TPRAD=SUMJ**.25D0;
C                   FL1=0.D0 ]
C          _OTHER [ FL1=0.D0 ]
C        _END;                        *
      FLUX(Km)=FL1
C- now only for conductivity
      IF(Km.LE.NCND)THEN
      Flum=CLUM*FLUX(Km)
      Accrad=0.D0
      ERAD=0.D0
      TpRAD=0.D0
      ELSE
      Flum=0.D0
      Accrad=0.D0
      SUMJ=0.D0
      SUMB=0.D0
      DO09963L=1,NFRUS
      SUMJ=SUMJ+Yout(NVARS*NZON+(NZON-NCND)*(L-1)-NCND+Km)*WEIGHT(L)
C-       If(Km==1) SUMB=SUMB+Black(L,Tp)*WEIGHT(L);
      If(Km.LT.NZON)Then
C-           Flum=Flum+Y(@NFHL+Km,1)*WEIGHT(L);
      Flum=Flum+Yout(NVARS*NZON+(NZON-NCND)*(L-1)-NCND+KRAD+Km)*WEIGHT(L)
C-         HAPHL=(WORK(Km+1,L)*DM(Km+1)+WORK(Km,L)*DM(Km))
C-                /(DM(Km)+DM(Km+1));
      HAPHL=(Tp**4+Ty(Km+1)**4)/(Ty(Km+1)**4*Uplim/(WORK(Km+1,L)+Haplim*WORK(Km,L))+Tp**4/WORK(Km,L))
C-           Accrad=Accrad+Y(@NFHL+Km,1)*WEIGHT(L)*HAPHL;
      Accrad=Accrad+Yout(NVARS*NZON+(NZON-NCND)*(L-1)-NCND+KRAD+Km)*WEIGHT(L)*HAPHL
      else
C-           Flum=Flum+Y(@NFHL+Km-KRAD,1)*WEIGHT(L)*HEDD(L);
      Flum=Flum+Yout(NVARS*NZON+(NZON-NCND)*(L-1)-NCND+KRAD+Km-KRAD)*WEIGHT(L)*HEDD(L)
      HAPHL=WORK(NZON,L)
C-           Accrad=Accrad+Y(@NFHL+Km-KRAD,1)*WEIGHT(L)*HAPHL*HEDD(L);
      Accrad=Accrad+Yout(NVARS*NZON+(NZON-NCND)*(L-1)-NCND+KRAD+Km-KRAD)*WEIGHT(L)*HAPHL*HEDD(L)
      endif
09963 CONTINUE
CScattst* testing scattering Tp balance *
      SUMJ=(1.5D1/PI**4)*SUMJ
      if(SUMJ.LE.0.d0)then
C-           @wterm ' SUMJ',SUMJ;
C-         pause;
C-       else;
C-         TPRAD=ABS(SUMJ)**.25D0; -- fails on some compilers for zero
      endif
      TPRAD=sqrt(sqrt(ABS(SUMJ)))
      TPRAD=SIGN(TPRAD,SUMJ)
      ERAD=RadC*SUMJ/rho(Km)
      Flum=CLUMF*Flum*Ry(Km)**2
      Accrad=CIMP*Accrad
CTHICK* *
C-     _If fullgp _Then
C-       <*THICK: *>
C-     _fi;
      ENDIF
      Acc(Km)=Accel+Accrad
      ArrLum(Km)=Flum
      EKIN=EKIN+RX*Uy(Km)**2
      EGRAV=EGRAV-RX*AM(Km)/Ry(Km)
C- ETERM=ETERM+EGAS*DM(Km)*UEI;
      ETERM=ETERM+WRK(Km,2)*DM(Km)*UEI
      ERADD=ERADD+ERAD*DM(Km)*UEI
      ACVIR=ACVIR+Ry(Km)*(Accel-Accrad)*RX
      IF(Km.NE.1.AND.Km.NE.NCND+1)EVIR=EVIR+Press(km)*UPI*(DM(Km)/Rho(km))
      IF(Km.NE.1)VISVIR=VISVIR+Q1*(DM(Km)/Rho(km))
C- VISVIR=VISVIR+QRT1*(DM(Km)/Rho(km)); -- R-T
      If(Km.EQ.IBURNT.or.Km.LE.20.or.Km.GE.NCND-10)KPRINT=Lpout
      If(fullgp.or.fullp)Then
      IF(.NOT.(KPRINT.EQ.Lpout))GOTO09960
      IF(.NOT.(CONV))GOTO09957
CC* PRINT LINE WITH CONVECTION *
      GOTO09956
09957 CONTINUE
C: PRINT LINE WITHOUT CONVECTION *
C-  XCARB=ACARB*YCARB;
      PLLOG=LOG10(max(1.D-50,Rho(km)))
      PLOG=LOG10(max(1.D-50,Press(km)))
      QVLOG=LOG10(max(1.D-50,Q1*UP))
      QRTLG=LOG10(max(1.D-50,QRTarr(km)*UP))
      AMPR=AM(Km)*UM
      If(AMPR.LT.1.01*(AMini*UM).AND.Km.GT.1)AMPR=AMPR-(AMini*UM)
      PL=RHO(Km)
      Tp=Ty(Km)
C- CHEM=CHEM0(Km);
      DO09954ii=1,Natom
      Yat(ii)=YABUN(ii,Km)
09954 CONTINUE
      RADP=.FALSE.
      CALLURSOS
      CHEM0(Km)=CHEM
      barionn=rho(km)*Urho*Avogar
      If(AMPR.GT.0.99D0*(AMOUT*UM).AND.Km.GT.Nzon/3)then
C-      AMPR=AMPR-(AMOUT*UM);
C: calculate AMPR - mass in external layers via dm() *
      AMPR=DMOUT
      DO09951ikm=Nzon,km,-1
      AMPR=AMPR+DM(ikm)
09951 CONTINUE
      AMPR=AMPR*UM
      WRITE(4,193)Km,-AMPR,min(999.99d0,Ry(Km)),min(99.9999d0,Uy(Km)*1.D+6/(UTIME*CRAP)),
     * Ty(Km),TpRAD,PLLOG,PLOG,QVLOG,QRTLG,WRK(Km,1),ENG,Flum,WRKX(Km),
     * Km,barionn,Xe*barionn,
C- barion n cm^{-3}, n_e
     *Yat(14)*barionn,xion(1,14),xion(2,14)
C- Fe, FeII, FeIII
C         Yat(7)*barionn,xion(1,7),xion(2,7), -- NaII, NaII xion(2,7)*barionn, -- NaII, NaIII
C         Yat(13)*barionn,xion(1,13),xion(2,13), -- CaII, CaIII
C*
      else
      WRITE(4,93)Km,AMPR,min(999.99d0,Ry(Km)),min(99.99,Uy(Km)*1.D+6/(UTIME*CRAP)),
C-      Ty(Km),PLLOG,PLOG,NREG(Km),ENGNUC,XCARB,ENG,
     *
C-      Ty(Km),PLLOG,PLOG,NREG(Km),XH,Accel,Accrad,
     *
C-      Ty(Km),PLLOG,PLOG,NREG(Km),Xion(1,1),Xion(1,2),DPE,
     *
C-      Ty(Km),PLLOG,Press(km)  ,Q1*UP ,Xion(1,1),Xion(1,2),TpRAD,
     *
C-      Ty(Km),(UTp/1.d3)*TpRAD,PLLOG,PLOG,QVLOG,QRTLG,
     *Ty(Km),TpRAD,PLLOG,PLOG,QVLOG,QRTLG,
C-      WRK(Km,1),CLUM*FLUX(Km),Flum,WRKX(Km),Km;
     *WRK(Km,1),ENG,Flum,WRKX(Km),Km,barionn,Xe*barionn,
C- barion n cm^{-3}, n_e
     *Yat(14)*barionn,xion(1,14),xion(2,14)
C- Fe, FeII, FeIII
C         Yat(7)*barionn,xion(1,7),xion(2,7), -- NaII, NaIIxion(2,7)*barionn, -- NaII, NaIII
C         Yat(13)*barionn,xion(1,13),xion(2,13), -- CaII, CaIII
C         Yat(14)*barionn,xion(1,14),xion(2,14); -- FeII, FeIII
C*
      endif
09956 CONTINUE
      KPRINT=0
09960 CONTINUE
      else
C: shock wave details in separate file *
      PLLOG=LOG10(max(1.D-50,Rho(km)))
      PLOG=LOG10(max(1.D-50,Press(km)))
      QVLOG=LOG10(max(1.D-50,Q1*UP))
C-    QRTLG=LOG10(max(1.D-50, QRTarr(km)*UP)); -- in place of eng
C: calculate AMPR - mass in external layers via dm() *
      AMPR=DMOUT
      DO09948ikm=Nzon,km,-1
      AMPR=AMPR+DM(ikm)
09948 CONTINUE
      AMPR=AMPR*UM
      PL=RHO(Km)
      Tp=Ty(Km)
C- CHEM=CHEM0(Km);
      If(Km.EQ.1.or.Km.EQ.Nzon)then
      WRITE(81,181)(tout-tretard)*(UTIME*CRAP)/86400.D0,
C- t in days
     *
C- tout*Utime/8.64d4,
     *Km,log10(AMPR),log10(UR*Ry(Km)),Uy(Km)*1.D+6/(UTIME*CRAP),log10(max(UTP*Ty(Km),1.d0)),log10(max(UTP*TpRAD,1.d0)),PLLOG,PLOG,QV
     *LOG,log10(max(eng,1.d-50)),Flum*1.d-40,WRKX(Km)
      else
      WRITE(81,281)0,Km,log10(AMPR),log10(UR*Ry(Km)),Uy(Km)*1.D+6/(UTIME*CRAP),log10(max(UTP*Ty(Km),1.d0)),log10(max(UTP*TpRAD,1.d0)
     *),PLLOG,PLOG,QVLOG,log10(max(eng,1.d-50)),Flum*1.d-40,WRKX(Km)
      endif
C- 181:FORMAT(1X,F10.2,I5,F12.6,F14.7,F9.4,6F9.3,1P,2E11.3); -- for secs
181   FORMAT(1X,F10.6,I5,F12.6,F14.7,F9.4,6F9.3,1P,2E11.3)
C- for days
281   FORMAT(10X,i1,I5,F12.6,F14.7,F9.4,6F9.3,1P,2E11.3)
      endif
      If(fullgp)Then
C:  PREPARE ARRAYS FOR BARKUK *
      WORK(Km,1)=LOG10(Ty(Km))
      WORK(Km,2)=LOG10(Rho(Km))
      WORK(Km,3)=LOG10(Press(km))
C-  WORK(Km,4)=SGAS;
      endif
      KPRINT=KPRINT+1
09966 CONTINUE
      IBURNT=0
      KENTR=0
      EKIN=EKIN/2.D0
      EVIR=3.D0*EVIR+Prce*Ry(1)**3-(Pcnd2+RadC3*UPI*Tpcnd2**4)*RNCND**3+Pcnd2*RNcnd2**3
C- EVIR=3.D0*EVIR-Pout*Rout**3+Prce*Rce**3;
      VISVIR=3.D0*VISVIR+Qrce*Ry(1)**3
      ETOT=ETERM+ERADD+EGRAV+EKIN
      If(NSTEP.EQ.0)ETOT0=ETOT
      BALENG=ETOT-ELOST-ETOT0
C-    write(*,'(a,i5,1p,4g15.3)')' nstep, Etot, Elost, Etot0, Baleng=',
C-          nstep, Etot, Elost, Etot0, Baleng;
C-pause;
      BALVIR=EVIR+VISVIR+EGRAV-ACVIR
C-TEFF=TPSURF*1.189207D0*UTP;
C: write integral values *
      WRITE(4,'(''%B:'')')
      write(4,'(4x,8(8x,a2))')'10','20','30','40','50','60','70'
      DO09945km=1,Nzon
      DO09942lfr=1,75
      if(lthick(lfr,km))then
      symb(lfr:lfr)='*'
      else
      symb(lfr:lfr)='o'
      endif
09942 CONTINUE
      write(4,'(i3,1x,a)')km,symb
09945 CONTINUE
      WRITE(4,'(A,(10I6))')'NTHICK: ',NTHICK
      WRITE(4,'(''  K/    tau:'',a)')' LubvB' 
      WRITE(4,'(8(1X,I3,''->'',1P,E9.2))')(K,tau(K),K=NZON,NCND+1,-1)
      WRITE(4,'(''  K/    tauU:'',a)')
      WRITE(4,'(8(1X,I3,''->'',1P,E9.2))')(K,tauU(K),K=NZON,NCND+1,-1)
      WRITE(4,'(''  K/    tauB:'',a)')
      WRITE(4,'(8(1X,I3,''->'',1P,E9.2))')(K,tauB(K),K=NZON,NCND+1,-1)
      WRITE(4,'(''  K/    tauV:'',a)')
      WRITE(4,'(8(1X,I3,''->'',1P,E9.2))')(K,tauV(K),K=NZON,NCND+1,-1)
      WRITE(4,'(''  K/    tauR:'',a)')
      WRITE(4,'(8(1X,I3,''->'',1P,E9.2))')(K,tauR(K),K=NZON,NCND+1,-1)
      WRITE(4,'(''  K/    tauI:'',a)')
      WRITE(4,'(8(1X,I3,''->'',1P,E9.2))')(K,tauI(K),K=NZON,NCND+1,-1)
      WRITE(4,'('' EFFECTIVE TEMPERATURE(KELVINS) & LOG10 ='',1P,2E15.5)')TEFF,DLOG10(TEFF)
      WRITE(4,'('' VOLUME GAINS POWER (1E+50 ERG/S) ='',1P,E15.5)')ELVOL*(UEPRI/(UTIME*CRAP))
      WRITE(4,'('' SURFACE LUMINOSITY (1E+50 ERG/S) ='',1P,E15.5)')ELSURF*(UEPRI/(UTIME*CRAP))
      WRITE(4,'('' TOTAL  GAINS POWER (1E+50 ERG/S) ='',1P,E15.5)')ELTOT*(UEPRI/(UTIME*CRAP))
      WRITE(4,'('' RADIAT.  ENERGY  (1E+50 ERG) ='',1P,E15.5,10X,            '' TOTAL    ENERGY  (1E+50 ERG) ='',E15.5)')ERADD*UEPRI
     *,ETOT*UEPRI
      WRITE(4,'('' KINETIC  ENERGY  (1E+50 ERG) ='',1P,E15.5,10X,            '' GAINED   ENERGY  (1E+50 ERG) ='',E15.5)')EKIN*UEPRI,
     *ELOST*UEPRI
      WRITE(4,'('' GRAVIT.  ENERGY  (1E+50 ERG) ='',1P,E15.5,10X,            '' VISCOUS  VIRIAL  (1E+50 ERG) ='',E15.5)')EGRAV*UEPRI
     *,VISVIR*UEPRI
      WRITE(4,'('' VIRIAL   ENERGY  (1E+50 ERG) ='',1P,E15.5,10X,            '' VIRIAL  BALANCE  (1E+50 ERG) ='',E15.5)')EVIR*UEPRI,
     *BALVIR*UEPRI
      WRITE(4,'('' THERMAL  ENERGY  (1E+50 ERG) ='',1P,E15.5,10X,            '' TOTAL   BALANCE  (1E+50 ERG) ='',E15.5)')ETERM*UEPRI
     *,BALENG*UEPRI
      WRITE(76,'(''ENERGY BALANCE AT OBS. TIME ='',1P,F15.5,2x,''DAYS : '',               E15.5)')(tout-tretard)*(UTIME*CRAP)/86400.
     *D0/10.,BALENG*UEPRI/(ETOT*UEPRI)
      WRITE(4,'('' <======== PARAMETERS OF PHOTOSPHERE ========> '')')
      WRITE(4,'(''   OB.T(D)         TEFF(10**3) RPH         VPH'',            ''         Xph     HUSED     MBOL      -U-       -B-'
     *',            ''       -V-       -U-B-     -B-V-'')')
      WRITE(4,'(1X,1P,G12.5,'' E'',0P,F9.3,3F12.3,1P,E10.2,0P,6F10.3)')(tout-tretard)*(UTIME*CRAP)/8.64D+04,TEFF*1.D-03,RPH,VPH*1.D+
     *06/(UTIME*CRAP),XPH,HUSED*UTIME*CRAP,MBOL,MU,MB,MV,MU-MB,MB-MV
C-  CALL PRIMAG;
C- WRITE(@Wres,*)' fullgp=',fullgp,'   Ipout=',Ipout;
      IF(.NOT.(fullgp))GOTO09939
C: draw graphic *
C- WRITE(@Wres,'(''%GR:'')');
CBar* Barkuk for plot T,P,Rho of the model *
C: transform outgoing flux (Nzon-1) to rest observer frame *
      nnx=0
      Nrec=0
      CLUMNU=LOG10(8.D0*UR**2*(BOLTZK*UTP)**3/(CS*HPLANC)**2)
C: obtain in ttt(NFREQ) corrected output flux
C                      na R**2  umnozhitj ! *
      Hnorm=0.d0
      DO09936Lfr=1,NFRUS+dLfrMax
      fr0=freqob(Lfr)*(1.d0-uy(Nzon)/clight)
      fr0log=log(fr0)
      xnu=(zerfr-fr0log)*dlognu(1)
C- dlognu < 0 !
      Lfr0=min(Nfreq-1,Nfrus,max(int(xnu),1))
      L=Lfr0
C-        Hobsg=cH(Lfr)*exp((1.d0-wH(Lfr))*log(max(Y(@L1),1.d-100))+
C-                       wH(Lfr)*log(max(Y(@L),1.d-100)));
      if(Lfr0.LT.Nfrus)then
      Hobsg=cH(Lfr)*exp((1.d0-wH(Lfr))*log(max(Yout(NVARS*NZON+(NZON-NCND)*(L+1)-1+KRAD),1.d-100))+wH(Lfr)*log(max(Yout(NVARS*NZON+(
     *NZON-NCND)*L-1+KRAD),1.d-100)))
      else
      Hobsg=cH(Lfr)*Yout(NVARS*NZON+(NZON-NCND)*L-1+KRAD)
      endif
      ttt(Lfr)=Hobsg*RY(NZON-1)**2
C-comoving  L=Lfr;     Yout(@LL)*RY(NZON-1)**2;
      Hnorm=Hnorm+ttt(Lfr)*WEIGHT(Lfr)
09936 CONTINUE
      Obflum=CLUMF*Hnorm
C- obs.lum.
      SUM=0.
      DO09933L=1,NFRUS+dLfrMax
      ttt(L)=SNGL(LOG10(MAX(ABS(ttt(L)*freqob(L)**3),1.d-100)))
      ttt(L)=ttt(L)+CLUMNU
C- erg/(sec*Hz)
      SUM=SUM+10.**ttt(L)*(freq(L+1)-freq(L))*(BoltzK*UTp/(2.d0*pi*Hplanc))
09933 CONTINUE
      WRITE(4,*)' Rest observer Lum=',Obflum,SUM
C-    pause;
      WRITE(4,*)'   L        lg nu         comov        rest',' observer lg(L_nu)  erg/(s Hz)'
      DO09930L=1,Nfrus
      write(4,'(i4,3f14.3)')L,log10(freqmn(L)*UFREQ),
C-               SNGL(LOG10(MAX(ABS(Y(@L)*RY(NZON-1)**2*FREQMN(L)**3),
     *SNGL(LOG10(MAX(ABS(Yout(NVARS*NZON+(NZON-NCND)*L-1+KRAD)*RY(NZON-1)**2*FREQMN(L)**3),1.d-100)))+Clumnu,ttt(L)
09930 CONTINUE
C:  PRINT FJ & FH *
      WRITE(STR2,'(10I9)')(I,I=1,MIN(NFRUS,10))
C- HEADER FOR FREQ
      WRITE(4,'(''  L:'',A)')STR2
      WRITE(4,'(''  K/  lg[ J_nu erg/(cm**2 s Hz strad)]  '')')
      Lfin=MIN(NFRUS,10)
      CJnu=(BOLTZK*UTP)**3/(2.d0*(pi*CS*HPLANC)**2)
C- Hplanck is hbar !!!
      DO09927K=NCND+1,NZON
      WRITE(4,'(1X,I4,2X,10F9.3)')K,(log10(max(CJnu*freqmn(L)**3
C-           *Y(NVARS*NZON+(K-NCND)+(NZON-NCND)*(L-1),1),1.d-299)),
     **Yout(NVARS*NZON+(K-NCND)+(NZON-NCND)*(L-1)),1.d-299)),L=1,Lfin)
C- FJ
      If(NFRUS.GT.10)WRITE(4,'(7X,10F9.3)')(log10(max(CJnu*freqmn(L)**3
C-           *Y(NVARS*NZON+(K-NCND)+(NZON-NCND)*(L-1),1),1.d-299)),
     **Yout(NVARS*NZON+(K-NCND)+(NZON-NCND)*(L-1)),1.d-299)),L=11,NFRUS)
C- FJ
09927 CONTINUE
C- H:
      WRITE(STR2,'(10I9)')(I,I=1,MIN(NFRUS,10))
C- HEADER FOR FREQ
      WRITE(4,'(''  L:'',A)')STR2
      WRITE(4,'(''  K/  lg[ H_nu erg/(cm**2 s Hz strad)]  '')')
      DO09924K=NCND+1,NZON
      WRITE(4,'(1X,I4,2X,10F9.3)')K,(log10(max(CJnu*freqmn(L)**3*Yout(NVARS*NZON+KRAD+(K-NCND)+(NZON-NCND)*(L-1)),1.d-299)),L=1,Lfin
     *)
C- FH
      If(NFRUS.GT.10)WRITE(4,'(7X,10F9.3)')(log10(max(CJnu*freqmn(L)**3*Yout(NVARS*NZON+KRAD+(K-NCND)+(NZON-NCND)*(L-1)),1.d-299)),L
     *=11,NFRUS)
C- FH
09924 CONTINUE
C
C    WRITE(STR2,'(10I9)') (I,I=1,MIN(NFRUS,10)); -- HEADER FOR FREQ
C    WRITE(@Wres,'(''  L:'',A)') STR2;
C    WRITE(@Wres,'(''  K/   FJ  '')');
C    Lfin=MIN(NFRUS,10);
C    _Do K=NCND+1,NZON;
C      WRITE(@Wres,45) K,(
C--         Y(NVARS*NZON+(K-NCND)+(NZON-NCND)*(L-1),1),LTHICK(L,K),
C         Yout(NVARS*NZON+(K-NCND)+(NZON-NCND)*(L-1)),LTHICK(L,K),
C          L=1,Lfin);  -- FJ
C      If(NFRUS>10)
C        WRITE(@Wres,46)((
C--           Y(NVARS*NZON+(K-NCND)+(NZON-NCND)*(L-1),1)),LTHICK(L,K),
C           Yout(NVARS*NZON+(K-NCND)+(NZON-NCND)*(L-1))),LTHICK(L,K),
C           L=11,NFRUS);  -- FJ
C    _od;
C    WRITE(@Wres,'(''  L:'',A)') STR2;
C    WRITE(@Wres,'(''    K/ FH  '')');
C    _Do K=NCND+1,NZON;
C      WRITE(@Wres,45)K,(
C--         Y(NVARS*NZON+KRAD+(K-NCND)+(NZON-NCND)*(L-1),1),LTHICK(L,K),
C         Yout(NVARS*NZON+KRAD+(K-NCND)+(NZON-NCND)*(L-1)),LTHICK(L,K),
C            L=1,Lfin);  -- FH
C      If(NFRUS>10)
C         WRITE(@Wres,46)((
C--         Y(NVARS*NZON+KRAD+(K-NCND)+(NZON-NCND)*(L-1),1)),LTHICK(L,K),
C         Yout(NVARS*NZON+KRAD+(K-NCND)+(NZON-NCND)*(L-1))),LTHICK(L,K),
C                      L=11,NFRUS);  -- FH
C    _od;
C*
45    FORMAT(1X,I4,2X,1P,10(E9.2,L2))
46    FORMAT(7X,1P,10(E9.2,L2))
      Lbeg=1
      Lfin=MIN(NFRUS,10)
09921 IF(.NOT.(Lbeg.LT.Nfrus))GOTO09920
      WRITE(4,'(//)')
      WRITE(STR2,'(10I9)')(I,I=Lbeg,Lfin)
C- HEADER FOR FREQ
      WRITE(4,'(''        L:'',A)')STR2
      WRITE(4,'(''  ray(cm) /  lg[ I_nu erg/(cm**2 s Hz strad)]  '')')
      DO09918ip=0,NP
C- rays
      WRITE(4,'(1X,1p,e9.3,2X,0p,10F9.3)')pray(ip)*UR,(log10(max(CJnu*freqmn(L)**3*fout(ip,L),1.d-299)),L=Lbeg,Lfin)
09918 CONTINUE
      Lbeg=Lfin+1
      Lfin=MIN(NFRUS,Lfin+10)
      GOTO09921
09920 CONTINUE
C:  plots  for spectra *
      DO09915K=NCND+1,NZON
C- RADIATION
      If(NZON-NCND.GT.5)Then
      IF(.NOT.(K.EQ.NCND+1))GOTO09912
      KK=2
      GOTO09906
09912 CONTINUE
      IF(.NOT.(K.EQ.NCND+(NZON-NCND)/4))GOTO09911
      KK=3
      GOTO09906
09911 CONTINUE
      IF(.NOT.(K.EQ.NCND+(NZON-NCND)/2))GOTO09910
      KK=4
      GOTO09906
09910 CONTINUE
      IF(.NOT.(K.EQ.NCND+3*(NZON-NCND)/4))GOTO09909
      KK=5
      GOTO09906
09909 CONTINUE
      IF(.NOT.(K.EQ.NZON-1))GOTO09908
      KK=6
      GOTO09906
09908 CONTINUE
      IF(.NOT.(K.EQ.42.or.K.EQ.43))GOTO09907
      kk=7
      GOTO09906
09907 CONTINUE
      GOTO09915
09906 CONTINUE
      else
      KK=K-NCND+1
      endif
      DO09905L=1,NFRUS
      BL(L)=MAX(BLACK(L,Ty(K))*FREQMN(L)**3,FJNOIS)
      FJ(L)=Y(NVARS*NZON+K-NCND+(NZON-NCND)*(L-1),1)*FREQMN(L)**3
      FH(L)=Y(NVARS*NZON+KRAD+K-NCND+(NZON-NCND)*(L-1),1)*FREQMN(L)**3
C- If(KK==2)WORK(L,1)=SNGL(LOG10(BL(L))); -- was used for Nfreq<Mzon
C- WORK(L,KK)=SNGL(LOG10(MAX(ABS(FJ(L)),FJNOIS)));
C- WORKX(L)=FLOAT(L);
09905 CONTINUE
      WRITE(4,'(A,I3,A,1P,E12.3,A,(1P,8E12.3))')' K=',K,'   tau=',tau(K),' FJ(L,K)=',(FJ(ill),ill=1,NFRUS)
09915 CONTINUE
CBAR*  use  BARKUK *
C: tau in all frequencies for selected zone *
      Nztau=43
      write(4,*)' L / lambda / tau in zon / tauabs / hapabs / hapsum /',' hapscat / cooling / heating ',' for Km = ',Nztau
      DO09902LL=1,Nfrus
      tauzon(LL)=0.d0
09902 CONTINUE
      DO09899Km=Nzon,Nztau,-1
      PL=RHO(Km)
      Tp=Ty(Km)
C- CHEM=CHEM0(Km);
      DO09896ii=1,Natom
      Yat(ii)=YABUN(ii,Km)
09896 CONTINUE
      RADP=.FALSE.
      CALLURSOS
      kmhap=Km
      CALLHAPPA
      CallOpacit
      DO09893LL=1,Nfrus
      tauzon(LL)=tauzon(LL)+(Ry(Km)-Ry(Km-1))*HAPPAL(LL)*PL
      if(Km.EQ.Nztau)then
      write(4,'(1x,i4,1p,2e13.4,6e11.3)')ll,Cs*1.d8/(freqmn(ll)*Ufreq),(Ry(Km)-Ry(Km-1))*HAPPAL(LL)*PL,(Ry(Km)-Ry(Km-1))*HAPAbs(LL)*
     *PL,log10(max(1.d-50,hapabs(ll)/(ur*urho))),log10(max(1.d-50,happal(ll)/(ur*urho))),log10(max(1.d-50,(happal(ll)-hapabs(ll))/(u
     *r*urho))),HAPAbs(LL)*black(LL,Tp)*weight(LL),HAPAbs(LL)*Y(NVARS*NZON+Nztau-NCND+(NZON-NCND)*(LL-1),1)*weight(LL)
      endif
09893 CONTINUE
09899 CONTINUE
      write(4,*)' L / lambda / tau for Km = ',Nztau
      write(4,'(1x,i5,2f15.5)')(ll,Cs*1.d8/(freqmn(ll)*Ufreq),tauzon(ll),ll=1,NFRUS)
09939 CONTINUE
      RETURN
C-  ENTRy PRIBUR(IPOUT);  -- PRINT BURNT ZONE NUMBER IPOUT
C-  <*E* PRINT ZONES IPOUT-1, IPOUT, IPOUT+1 *>;
C-  RETURN;
C: FORMATS*
C-93:FORMAT(1X,I3,1P,4E13.5,0P,2F8.4,I5,1P,5E10.2,I5);
93    FORMAT(1X,I3,F9.5,F12.7,F8.4,F10.3,F8.3,4F7.2,1P,4E10.2,I5,11e10.2)
193   FORMAT(1X,I3,1P,E9.2,0P,F12.7,F8.4,F10.3,F8.3,4F7.2,1P,4E10.2,I5,11e10.2)
94    FORMAT(' ZON',3X,'AM/SOL',6X,'R',F3.0,6X,'V',F3.0,5X,'T',F3.0,4X,'Trad5 lgD',F3.0,' lgP',F3.0,
C- 4X,'Trad3 lgD',F3.0,' lgP',F3.0,
     *
C- ' NREG   XHII      XHEII     DPE       LUM      CAPPA    ZON');
     *
C- ' NREG   XH        Accel     Accrad    LUM      CAPPA    ZON');
     *
C- ' NREG   XHII      XHEII     TpRAD     LUM      CAPPA    ZON');
     *
C- ' Qvis       XHII      XHEII     LUM      CAPPA    ZON');
     *
C- '  lgQv  lgQRT      XHII     LUMcnd     LUM     CAPPA   ZON');
     *'  lgQv  lgQRT      XHI      ENG        LUM     CAPPA   ZON','    n_bar     n_e          Fe        II        III')
C-     ,'    n_bar     n_e       Na        II        III       '
C-       ,  'Ca        II        III       Fe        II        III');
C- ' NREG   ENGNUC    XCARB     ENG       LUM      CAPPA    ZON');
95    FORMAT(' ZON',5X,'AM/SOL',6X,'R',F3.0,10X,'V',F3.0,10X,'T',F3.0,5X,'LG D',F3.0,' LG P',F3.0,' NREG   ENGNUC    XCARB     UCONV
     *     ENTROP   FLCONV   ZON')
      END
C: PRIMAG - PRINT MAGNITUDES *
      SUBROUTINEPRIMAG
      IMPLICITREAL*8(A-H,O-Z)
      Parameter(MAXIT=15)
      Parameter(EPSTPH=1.D-3)
C-      @wterm  "print*,";
C-NVARS - number of independent variables
      PARAMETER(NVARS=3)
C- 3 - NOCONV,4 - CONV
      include '../obj/nfreq_and_mzone.inc'
C-  PARAMETER(Mzon=90); -- h toy model
C- PARAMETER(Mzon=171); -- zoning for W7
C- PARAMETER(Mzon=340); -- zoning for W7+wind
C-PARAMETER(Mzon=86); -- zoning for W7fhh+wind
C-PARAMETER(Mzon=43); -- zoning for W7fhh
C-   PARAMETER (Mzon=521); -- zoning for ntomi Crab model
      PARAMETER(NYDIM=(NVARS+2*NFREQ)*Mzon,MAXDER=4)
C- PARAMETER(NYDIM=2*NFREQ*Mzon,MAXDER=4);
C- Nydim must be  2*(@Npe+@Npc)*Mzon  to use in Feau
C- PARAMETER(NYDIM=2*(@Npe+@Npc)*Mzon,MAXDER=4); -- to use in Feau
      Parameter(Is=5)
C- for test -- chego???
C-   PARAMETER (NZ=1200000); --  for Nfreq=200, Mzon=200
      PARAMETER(NZ=3000000)
C-  for Nfreq=100, Mzon=600
      Parameter(Nstage=28,Natom=15)
      PARAMETER(KOMAX=80)
C-MAX NUMBER OF STEPS IN PRESCRIBED MOMENTS
      LogicalLSYSTEM
C- /.TRUE./; -- for IBM
      Parameter(LSystem=.FALSE.)
      Parameter(Pi=3.1415926535897932d+00,hPlanc=1.0545716280D-27,Cs=2.9979245800D+10,Boltzk=1.3806504000D-16,Avogar=6.0221417900D+2
     *3,AMbrun=1.6605387832D-24,AMelec=9.1093821500D-28,echarg=4.8032042700D-10,CG=6.6742800000D-08,CMS=1.9884000000D+33,RSol=6.9551
     *000000D+10,ULGR=1.4000000000D+01,UPURS=1.0000000000D+07,ULGPU=7.0000000000D+00,ULGEU=1.3000000000D+01,UPC=3.0856776000D+18,UTP
     *=1.0000000000D+05,URHO=1.0000000000D-06,CARAD=7.5657680191D-15,CSIGM=5.6704004778D-05,ERGEV=1.6021764864D-12,GRADeV=1.16045052
     *85D+04,RADC=7.5657680191D-02,CTOMP=4.0062048575D-01,CCAPS=2.6901213726D+01,CCAPZ=9.8964034725D+00)
      IntegerZn(Natom),ZnCo(Natom+1)
      DimensionAZ(Natom)
      Common/AZZn/AZ,Zn,ZnCo
C
C      include for eve and strad
C*
      Common/NiAdap/tday,t_eve,XNifor(Mzon),AMeveNi,KNadap
      LOGICALFRST
      Parameter(Mfreq=130)
C- think about max Nfreq !
      Common/Kmzon/km,kmhap,Jac,FRST
C-   Common/NiAdap/tday,t_eve,XNifor(Mzon),AMeveNi,KNadap; -- must go to commonEve.inc
C- since it enters also bgcon*trf and is forgotten there
      COMMON/STCOM1/t,H,HMIN,HMAX,EPS,N,METH,KFLAG,JSTART
      COMMON/YMAX/YMAX(NYDIM)
      COMMON/YSTIF/Y(NYDIM,MAXDER+1)
      COMMON/HNUSED/HUSED,NQUSED,NFUN,NJAC,NITER,NFAIL
      COMMON/HNT/HNT(7)
C- FOR COSETBGH
      PARAMETER(DELTA=1.d-05)
C- DISTURBANCE CONSTANT standard was 1.d-04
C- cannot be less ~1.d-6 since opacity is saved in single precision
C-   PARAMETER (LICN=4*NZ,LIRN=8*NZ/5); -- FOR F01BRF, M28Y12
      PARAMETER(LICN=4*NZ,LIRN=2*NZ)
C- FOR F01BRF, M28Y12
C- maximum LIRN=2*NZ
      LogicalNEEDBR
C- F01BRF or M28Y12 needed
      COMMON/STJAC/THRMAT,HL,AJAC(LICN),IRN(LIRN),ICN(LICN),WJAC(NYDIM),FSAVE(NYDIM*2),IKEEP(5*NYDIM),IW(8*NYDIM),IDISP(11),NZMOD,NE
     *EDBR
      LOGICALCONV,
C- .TRUE. IF CONVECTION IS INCLUDED
     *CHNCND,
C- .TRUE. IF NCND IS CHANGED
     *SCAT,
C- .TRUE. IF Scattering included (Hapabs^=Happa)
     *SEP
C- .TRUE. IF absorption and scattering are in separate files
C-           (i.e. absorption without expansion effect)
      COMMON/CUTOFF/FLOOR(NVARS+1),Wacc(NVARS+1),FitTau,TauTol,Rvis,CONV,CHNCND,SCAT,SEP
      LogicalLTHICK
      COMMON/THICK/LTHICK(Nfreq,Mzon)
      COMMON/CONVEC/UC(Mzon),YAINV(Mzon)
      COMMON/RAD/EDDJ(Mzon,Nfreq),EDDH(Mzon),HEDD(Nfreq),HEDRAD,CLIGHT,CKRAD,UFREQ,CFLUX,CCL,CLUM,CLUMF,CIMP,NTHICK(NFREQ),NTHNEW(NF
     *REQ),bolM,NCND,KRAD,NFRUS
C- KRAD=(Mzon-NCND)*NFRUS IN STELLA
      LOGICALEDTM
C- Eddington factors time-dependent==.true.
      COMMON/RADOLD/HEDOLD,HINEDD,EDTM
      Common/newedd/EddN(Mzon,Nfreq),HEdN(Nfreq),tfeau
C- for EDTM==T
      Common/oldedd/EddO(Mzon,Nfreq),HEdo(Nfreq),trlx
C- for EDTM==T
      Common/cnlast/Cnlast
      Common/Dhap/DHaphR(Mzon,Nfreq)
      COMMON/BAND/FREQ(NFREQ+1),
C- frequency bins boundaries,
     *
C- frequency in units of h/(k*Tpunit)
     *FREQMN(NFREQ),
C- frequency mean (middle of bin)
     *WEIGHT(130),
C- bandwidth*freq(mean)**3
     *HAPPAL(NFREQ),
C- dimens-less kappa for Ron's absorp.+scatt.
     *HAPABSRON(NFREQ),
C- dimens-less kappa for Ron's absorption
     *HAPABS(NFREQ),
C- dimens-less kappa for S-B absorption
     *DLOGNU(NFREQ)
C- (log step in frequency)**-1
C- PARAMETER(NFRMIN=Nfreq/5); -- MINIMUM NFRUS
      PARAMETER(NFRMIN=Nfreq/2)
C- MINIMUM NFRUS
      IntegerdLfrMax
      Common/observer/wH(Mfreq),cH(Mfreq),zerfr,Hcom(Mfreq),Hobs(Mfreq),freqob(Mfreq),dLfrMax
      Parameter(NP=15+15-1)
C- Total number of tangent rays
      Common/famu/fstatic(0:NP+1,Nfreq),fobs_corr(0:NP+1,Mfreq),fcom(0:NP+1,Mfreq),amustatic(0:NP+1)
      Common/rays/Pray(0:Np+1),fout(0:NP+1,Mfreq),abermu(0:NP+1),NmuNzon
C- fout probably not needed
      COMMON/LIM/Uplim,Haplim
C- to BEGRAD
      COMMON/AMM/DMIN,DM(Mzon),DMOUT,AMINI,AM(Mzon),AMOUT
C- exactly as in VELTEM
      COMMON/Centr/RCE,Nzon
C- central radius & current Number of zones
      Common/InEn/AMHT,EBurst,tBurst,tbeght
C- Mass of Heated Core, Energy & time
      COMMON/RADPUM/AMNI,XMNi,XNi,KmNick
C- MASS OF Ni CORE, Ni ABUNDANCE
      COMMON/RADGAM/FJgam(Mzon,2),toldg,tnewg
C- zero-moment for gamma-radiation
C- 1 for old time (toldg)      2 for present time (tnewg)
      COMMON/RADGlg/FJglog(Mzon,2)
C- log FJgam
      COMMON/CHEM/CHEM0(Mzon),RTphi(0:Mzon+1),EpsUq
      COMMON/REGIME/NREG(Mzon)
      doubleprecisionNRT
      COMMON/AQ/AQ,BQ,DRT,NRT
C- BQ pseudo-viscosity for R-T
C- DRT is the mass distance or optical thickness
C- and used for weight in artificial viscosity
C- NRT is the power in pseudo-viscosity (may be noninteger)
      COMMON/AZNUC/ACARB,ZCARB,ASI,ZSI,ANI,ZNI,QCSI,QSINI
      COMMON/QNRGYE/QNUC,RGASA,YELECT
C-   COMMON/CKN1/CK1,CK2,CFR,CRAP;
      COMMON/CKN1/CK1,CK2,CFR,CRAP,CRAOLD
      LOGICALEVALJA,
C- EVALUATE JACOBIAN
     *OLDJAC,
C- JACOBIAN MAY NEED UPDATING
     *BADSTE
C- BAD STEP: RETURNED BY DIFMAT - TRY SMALLER
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
C- UNITS:
      COMMON/UNSTL/UL,UPRESS,UE,UEPS,UCAP,UTIME,UV,UFLUX,UP
      COMMON/TAIL/KTAIL
      COMMON/UNINV/UPI,UEI
C- 1./UP, 1./UE
      COMMON/UNBSTL/UR,UM,UEPRI,ULGP,ULGE,ULGV,ULGTM,ULGEST,ULGFL,ULGCAP,ULGEPS
      COMMON/CONUR/EIT,DST,BBRCNR(5)
      COMMON/BAL/EL(MAXDER+1),YENTOT(MAXDER+1),ETOT0,ELVOL,ELSURF,ELTOT,TPSURF,HOLDBL,ELOST,EKO,RADBEG
      common/NSTEP/NSTEP,NDebug,MAXER,IOUT,
C- IOUT      LINE PRINT INTERVAL
     *NOUT
C- NOUT      PRINT STEP
      common/CREAD/TAUOLD,NSTMAX,MBATCH,MAXORD
      common/debug/LfrDebug,Nperturb,Kbad
      REAL*8TPMAX(MAXDER+1),TQ(4)
      COMMON/TAU/TAU(Mzon+1),FLUX(Mzon)
      common/tauubvri/tauU(Mzon),tauB(Mzon),tauV(Mzon),tauR(Mzon),tauI(Mzon)
      COMMON/PHOT/XJPH,DMPH,RPH,TPH,PLPH,VPH,CHEMPH,GRVPH,HP,JPH
      PARAMETER(NFUNC=6)
      REAL*4WORK(Mzon+2,NFREQ)
CNFUNC IF NFREQ < NFUNC*
     *,WRK(Mzon,4)
      REAL*8WRKX(Mzon),WORKX(Mzon+2)
      COMMON/STEPD/WRKX,WORKX,TPHOT,TEFF,WORK,WRK,NPHOT,NZM
C       1 - LG(T), 2 - LG(PL), 3 - LG(P), 4 - LG(S)      *
      PARAMETER(TMCRIT=1.D-6,TPNSE=5.D0,EPGROW=0.02D0)
      Common/RUTP/Ry(Mzon),Uy(Mzon),Ty(Mzon),Press(Mzon),Rho(Mzon)
      COMMON/TOO/TOO,KO,KNTO,TO(KOMAX),STO(KOMAX),NTO(KOMAX)
      Parameter(Lcurdm=1000)
C- Dimension of the Tcurv array
      RealTcurv
      IntegerNFRUSED
C- an integer array which store exact number of used freqs
      REAL*8Flsave
C- remove this to read old flx files!
      Common/Curve/tcurv(8,Lcurdm),Depos(Lcurdm),Flsave(MFREQ+1,Lcurdm),NFRUSED(Lcurdm),Lsaved
      LOGICALBEGRUN
      Common/BEGR/BEGRUN
      CHARACTER*80Model,Sumprf,Sumcur,Depfile,Flxfile
      COMMON/Files/Model,Sumprf,Sumcur,Depfile,Flxfile
      CHARACTER*1app
C- dummy var for constructing of Opafile
      LogicalGivdtl
      Common/ABGrap/NSTA,NSTB,TcurA,TcurB,Givdtl
C-No. of steps & Time in days
      REAL*8MBOL,MU,MB,MV,MR,MI,MBOL1
      COMMON/COLOR/MBOL,MU,MB,MV,MR,MI,UMB,BMV,MBOL1,LubvU,LubvB,LubvV,LubvR,LubvI,Lyman
      COMMON/DETAIL/QRTarr(Mzon),UUarr(Mzon),ArrLum(Mzon),Acc(Mzon)
      Common/XYZ/XA,YA,URM
C- _TRACE "WRITE(0,'(A,1P,2E12.3,A)')'YL YL1 ',Y(@L),Y(@L1),"
C-_TRACE "@wterm' Jph=',Jph,' SumH=',SumH,' Teff=',Teff,"
C-_TRACE "@wterm' depos(1), depos(20)',depos(1),depos(20),"
C:  PARAMETERS & COMMONS *
C-    Common/observer/wH(Nfreq),cH(Nfreq),zerfr;
      common/phofit/Tphfit,Rphfit
C-  REAL*8 BL(NFREQ),FJ(NFREQ),FH(NFREQ); -- FOR %MG:
      REAL*8Ryp(Mzon),Uyp(Mzon),Typ(Mzon),Rhop(Mzon)
C- for primag
      Character*80Opafile
C- Common/opabs100/opacabs100(Nfreq); -- bad if here; now in argument for subr. 'opacity"
C- real*4 hptabab,hptabsc,hpsavsc; -- new!
      real*4hptabab,hptababron,hptabsc,hpsavsc
C- new!
      Common/Opsave/hpsavsc(Nfreq,14,14,Mzon/(Mzon/50)+1,6)
      Common/OpBand/TpTab(14),RhoTab(14),STab(6),Wavel(Nfreq),YATab(Natom),hptabab(Nfreq,14,14,Mzon/(Mzon/50)+1,2),hptababron(Nfreq,
     *14,14,Mzon/(Mzon/50)+1,2),hptabsc(Nfreq,14,14,Mzon/(Mzon/50)+1,2),
C-Peter Hoeflich:       EpsBand(@Nfreq,@Ntab,@Ntab,@Ns),
     *
C-        EppBand(@Nfreq,@Ntab,@Ntab,@Ns),
     *Msta,Nrho,Ntp,im
C- im number of mixture
C- old Dimension hpbanab1(Nfreq,@Ntab,@Ntab,Mzon/@skip),
      real*4hpbanab1(Nfreq,14,14,Mzon/(Mzon/50)+1),hpbanabron1(Nfreq,14,14,Mzon/(Mzon/50)+1),hpbansc1(Nfreq,14,14,Mzon/(Mzon/50)+1),
     *hpbanab2(Nfreq,14,14,Mzon/(Mzon/50)+1),hpbanabron2(Nfreq,14,14,Mzon/(Mzon/50)+1),hpbansc2(Nfreq,14,14,Mzon/(Mzon/50)+1)
      Equivalence(hptabab(1,1,1,1,1),hpbanab1(1,1,1,1)),(hptababron(1,1,1,1,1),hpbanabron1(1,1,1,1)),(hptabsc(1,1,1,1,1),hpbansc1(1,
     *1,1,1)),(hptabab(1,1,1,1,2),hpbanab2(1,1,1,1)),(hptababron(1,1,1,1,2),hpbanabron2(1,1,1,1)),(hptabsc(1,1,1,1,2),hpbansc2(1,1,1
     *,1))
      Common/tintrp/stmlog(6),tdlog,thaplog1,thaplog2,istold,Opafile
      Common/dumfreq/dumFreq(Nfreq+1),dumFreqmn(Nfreq),dumwavel(Nfreq)
      BLACK(Lbl,Tpbl)=(exp(-(FREQMN(Lbl)/Tpbl)))/(1.d0-(exp(-(FREQMN(Lbl)/Tpbl))))
      BLACKD(Lbl,Tpbl)=(FREQMN(Lbl)/Tpbl)*(exp(-(FREQMN(Lbl)/Tpbl)))/(1.d0-(exp(-(FREQMN(Lbl)/Tpbl))))**2
      DO09890K=1,Nzon
      Ryp(K)=Y(K,1)
      Uyp(K)=Y(Nzon+K,1)
C- NR
C-      Uyp(K)=Y(Nzon+K,1)/sqrt(1.d0+(Y(Nzon+K,1)/clight)**2); -- Rel
      Typ(K)=Y(2*Nzon+K,1)
09890 CONTINUE
      Rhop(1)=3.D0*DM(1)/(Ryp(1)**3-Rce**3)
      DO09887I=1,NZON
      IF(I.GT.1)Rhop(I)=3.D0*DM(I)/(Ryp(I)**3-Ryp(I-1)**3)
09887 CONTINUE
C: FIND PARAMETERS OF PHOTOSPHERE *
      I=NZON
      tau(I)=0.D0
      tauTp=0.D0
      IRph=0
      kmrep=0
09884 IF(.NOT.((I.GT.1.AND.(tau(I).LT.0.64D0.or.tauTp.LT.0.64D0)).or.I.GT.Nzon-5))GOTO09883
      PL=rhop(I)
      Tp=Typ(I)
C-       write(*,'(a,i5,1p,2e12.3)')' in primag I rhop Typ:',I,rhop(I),Typ(I);
      Doii=1,Natom
      Yat(ii)=YABUN(ii,I)
      enddo
      RADP=.FALSE.
      CALLURSOS
      kmhap=I
      CALLHAPPA
      tau(I-1)=tau(I)+(Ryp(I)-Ryp(I-1))*HAPPAL(LubvB)*PL
      I=I-1
      If(tau(I).GE.0.1d0.AND.kmRep.EQ.0)kmRep=max(I,kmnick)
      If(tau(I).GE.0.64D0.AND.IRph.EQ.0)IRph=I
      tauTp=0.5D0*(tau(I)+tau(I+1))
      GOTO09884
09883 CONTINUE
C- tau(I)>=@tauPH & tauTp>=@tauPH or I==1
      IRph=max(IRph,1)
      JPH=IRph
      RPHtau=Ryp(IRph+1)+(Ryp(IRph)-Ryp(IRph+1))/(tau(IRph)-tau(IRph+1))*(0.64D0-tau(IRph+1))
C If(NCND>0)Then;
C     <*COND: FIND RPH, VPH & TPH INTERPOLATING IN tau *;
C   ELSE;
C     <*TRAN: FITTING TPH FOR TRANSPARENT CORE, Then RPH & VPH*;
C   ENDIF;  *
C: FITTING TPH FOR TRANSPARENT CORE, Then RPH & VPH*
      SUMH=0.D0
      DO09881L=1,NFRUS
C-       write(*,'(a,i5,1p,e12.4)')' primag L Hedd:', L,Hedd(L);
      SUMH=SUMH+Y(NVARS*NZON+(NZON-NCND)*L-1+KRAD,1)*WEIGHT(L)/(HEDD(L)*(1.D0-HEDD(L)))
09881 CONTINUE
      SUMH=(1.5D1/PI**4)*SUMH
      TPHo=ABS(SUMH)**.25D0
      TPHo=max(TPHo,1.d-2)
      Tph=TPHo
      RPH=SQRT(ABS(ELSURF)/(TPHo**4*CSIGM*UTp**4*UTIME**3/(URHO*UR**3)))
C- old Rph now via f_Edd at tau \sim 0.1 instead of h_Edd:
C    SUMH=0.D0;
C    _DO L=1,NFRUS;
C       SUMH=SUMH+Y(@NFHL+Kmrep,1)*WEIGHT(L)*
C           8.d0/(3.d0-6.d0*EDDJ(Kmrep,L)
C                   + sqrt(12.d0*EDDJ(Kmrep,L)-3.d0))
C    _OD;
C    SUMH=(1.5D1/PI**4)*SUMH;
C    TPHo=ABS(SUMH)**.25D0;
C    TPH=TPHo;
C    RPH=SQRT(ABS(ELSURF)/(TPHo**4*CSIGM*UTp**4*UTIME**3/(URHO*UR**3))); *
C-  If( MOD( NSTEP,NOUT) ==0)  then;
CRPHOT* least squares method *
C-      Tph=Tphfit; -- comment this for real run
C-      RPH=SQRT(ABS(ELSURF)/(TPH**4*CSIGM*UTp**4*UTIME**3/(URHO*UR**3)));
C-      @wterm ' new Tphfit Rph ', Tphfit, Rph;
C-  Endif;
      I=1
09878 IF(.NOT.(I.LT.NZON-1.AND.Ryp(I+1).LT.RPH))GOTO09877
      I=I+1
      GOTO09878
09877 CONTINUE
C- Ryp(I+1)>=RPH
      VPH=uyp(I+1)+(uyp(I)-uyp(I+1))/(Ryp(I)-Ryp(I+1))*(RPH-Ryp(I+1))
      XPH=YABUN(1,I+1)+(YABUN(1,I)-YABUN(1,I+1))/(Ryp(I)-Ryp(I+1))*(RPH-Ryp(I+1))
C-    If( MOD( NSTEP, MAX0(1,NOUT/10) )==0)
C-           @wterm ' primag Tph Rph Iph=',I, Tph, Rph;
C%MG:
C  _DO I=1,NZON;
C    WRKX(I)=FLOAT(I);
C    WRK(I,1)=SNGL(LOG10(Typ(I)));
C    WRK(I,2)=SNGL(SIGN(LOG10(ABS(uyp(I))+1.D-08),uyp(I)));
C  _OD;
C/_DO K=NCND+1,NZON;
C      _SELECT
C        _ K==NCND+1 [KK=2]
C        _ K==NCND+(NZON-NCND)/4   [KK=3]
C        _ K==NCND+(NZON-NCND)/2   [KK=4]
C        _ K==NCND+3*(NZON-NCND)/4   [KK=5]
C        _ K==NZON-1 [KK=6]
C        _OTHER  [_ITERATE K]
C      _END;
C    _DO L=1,NFRUS;
C      BL(L)=BLACK(L,Typ(K))*FREQMN(L)**3;
C      FJ(L)=Y(NVARS*NZON+K-NCND+(NZON-NCND)*(L-1),1)*FREQMN(L)**3;
C      FH(L)=Y(NVARS*NZON+KRAD+K-NCND+(NZON-NCND)*(L-1),1)*FREQMN(L)**3;
C      If(KK==2)WORK(L,KK-1)=SNGL(LOG10(ABS(BL(L))));
C      WORK(L,KK)=SNGL(LOG10(ABS(FJ(L))));
C      WORKX(L)=FLOAT(L);
C    _OD;
C  _OD;*/ CALL BARKUK (WORKX, WORK,Mzon+2,NFREQ,6,WORKX(1), WORKX(NFREQ),
C                  0.5,0., 0.,0.,60, 20, NFREQ*2, 0, 6 ); *
C: FIND TEFF,U,B,V *
      TEFF=TPH*UTp
C: FIND FH FOR U,B,V BY INTERPOLATION IN FHU,FHB,FHV *
      L=LubvU
C- FOR U
C-   print *, L,Y(@L),Y(@L1);
      If(L.LT.NFRUS.AND.abs(Y(NVARS*NZON+(NZON-NCND)*(L+1)-1+KRAD,1)).GT.1.d-100)Then
      FHU=EXP(LOG(ABS(Y(NVARS*NZON+(NZON-NCND)*(L+1)-1+KRAD,1)))+LOG(ABS(Y(NVARS*NZON+(NZON-NCND)*L-1+KRAD,1)/Y(NVARS*NZON+(NZON-NCN
     *D)*(L+1)-1+KRAD,1)))*DLOGNU(L)*LOG((CCL/3.65D+03)/FREQMN(L+1)))
C-           LOG(ABS(Y(@L)/Y(@L1)))/(FREQMN(L)-FREQMN(L+1))
C-       *(@U-FREQMN(L+1)));
      else
      FHU=ABS(Y(NVARS*NZON+(NZON-NCND)*L-1+KRAD,1))
      endif
      L=LubvB
C- FOR B
C-   print *, L, Y(@L),Y(@L1);
      If(L.LT.NFRUS.AND.abs(Y(NVARS*NZON+(NZON-NCND)*(L+1)-1+KRAD,1)).GT.1.d-100)Then
      FHB=EXP(LOG(ABS(Y(NVARS*NZON+(NZON-NCND)*(L+1)-1+KRAD,1)))+LOG(ABS(Y(NVARS*NZON+(NZON-NCND)*L-1+KRAD,1)/Y(NVARS*NZON+(NZON-NCN
     *D)*(L+1)-1+KRAD,1)))*DLOGNU(L)*LOG((CCL/4.4D+03)/FREQMN(L+1)))
C-           LOG(ABS(Y(@L)/Y(@L1)))/(FREQMN(L)-FREQMN(L+1))
C-       *(@B-FREQMN(L+1)));
      else
      FHB=ABS(Y(NVARS*NZON+(NZON-NCND)*L-1+KRAD,1))
      endif
      L=LubvV
C- FOR V
C-   print *, L,Y(@L),Y(@L1);
      If(L.LT.NFRUS.AND.abs(Y(NVARS*NZON+(NZON-NCND)*(L+1)-1+KRAD,1)).GT.1.d-100)Then
      FHV=EXP(LOG(max(1.d-100,ABS(Y(NVARS*NZON+(NZON-NCND)*(L+1)-1+KRAD,1))))+LOG(max(1.d-100,ABS(Y(NVARS*NZON+(NZON-NCND)*L-1+KRAD,
     *1)/Y(NVARS*NZON+(NZON-NCND)*(L+1)-1+KRAD,1))))*DLOGNU(L)*LOG((CCL/5.5D+03)/FREQMN(L+1)))
C-           LOG(ABS(Y(@L)/Y(@L1)))/(FREQMN(L)-FREQMN(L+1))
C-       *(@V-FREQMN(L+1)));
      else
      FHV=ABS(Y(NVARS*NZON+(NZON-NCND)*L-1+KRAD,1))
      endif
C-   PRINT*,@U,@B,@V,UFREQ;
C-   PRINT*,FHU,FHB,FHV
C- CFLUX=60/PI**4 * SIGMA * UTp**4
      MU=-2.5D0*(LOG10(max(CFLUX*FHU*(CCL/3.65D+03)**5/CCL*(Ryp(NZON-1)*UR/3.0857D19)**2,1.d-10))+8.37D0)
      MB=-2.5D0*(LOG10(max(CFLUX*FHB*(CCL/4.4D+03)**5/CCL*(Ryp(NZON-1)*UR/3.0857D19)**2,1.d-10))+8.18D0)
      MV=-2.5D0*(LOG10(max(CFLUX*FHV*(CCL/5.5D+03)**5/CCL*(Ryp(NZON-1)*UR/3.0857D19)**2,1.d-10))+8.42D0)
CG* PREPARE WRK,WRKX & CALL BARKUK TO CONTROL GAS DYNAMICS *
      ObsLum=0.d0
      DO09875L=1,NFRUS
      Obslum=Obslum+Y(NVARS*NZON+(NZON-NCND)*L-1+KRAD,1)*WEIGHT(L)*(1.d0+(uyp(Nzon-1)/clight)*(1.d0+EddJ(Nzon-1,L))/HEdd(L))
C- observer
C- in first order u/c -- need be corrected for all orders
09875 CONTINUE
      Obslum=max(CLUMF*Obslum*Ryp(NZON-1)**2,3.86d+31)
      MBOL=4.75D0-2.5D0*LOG10(ABS(ObsLum)/3.86d+33)
      if(mod(NSTEP,50).EQ.0)write(*,'(a,1p,2g15.5)')'logL',log10(Obslum)
C-write(@term,'(a,1p,2e12.3)')
C-  ' Obslum  Elsurf in lbalsw ',Obslum, Elsurf;
C-       ELSURF=Obslum/(UEPRI/(UTIME*CRAP));
C-comov  MBOL=4.75D0-2.5D0*LOG10(ABS(ELSURF)*(UEPRI/(UTIME*CRAP))/3.86D-17);
      time=(t-Ryp(NZON-1)/CLIGHT)*(UTIME*CRAP)/8.64D+04
C- Below we save outgoing flux and exact number of frequencies we used:
      Lsaved=Lsaved+1
      NFRUSED(Lsaved)=NFRUS+dLfrMax
C- for flx files and plots in tt* routines
      Hnorm=0.d0
      DO09872Lfr=1,NFRUS+dLfrMax
      fr0=freqob(Lfr)*(1.d0-uyp(Nzon)/clight)
      fr0log=log(fr0)
      xnu=(zerfr-fr0log)*dlognu(1)
C- dlognu < 0 !
C-        Lfr0=min(Nfreq-1,Nfrus,max(int(xnu),1));
      Lfr0=int(xnu)
      if(Lfr0.GE.1)then
C-        wfr0=(freqmn(Lfr0+1)-fr0)/(freqmn(Lfr0+1)-freqmn(Lfr0));
C- in linear scale
      L=Lfr0
C-	if(mod(Lfr,20)==0 .and. mod(Nstep,100)==0)then;
C-        if(Nstep==1000 .or. Nstep==1050 .or. Nstep==1100)then;
C-        if(Nstep==350 .or. Nstep==400 .or. Nstep==450)then;
C-          write(*,'(a,2I5)')' in lbalsw Lfr, Lfr0:',Lfr,Lfr0;
C-        pause;
C-        endif;
      if(Lfr0.LT.Nfrus)then
      Hobsg=cH(Lfr)*exp((1.d0-wH(Lfr))*log(max(Y(NVARS*NZON+(NZON-NCND)*(L+1)-1+KRAD,1),1.d-100))+wH(Lfr)*log(max(Y(NVARS*NZON+(NZON
     *-NCND)*L-1+KRAD,1),1.d-100)))
      else
      Hobsg=cH(Lfr)*Y(NVARS*NZON+(NZON-NCND)*L-1+KRAD,1)
      endif
      Flsave(Lfr,Lsaved)=Hobsg*Ryp(NZON-1)**2
C-comoving       Flsave(Lfr,Lsaved)=Y(@L)*Ryp(NZON-1)**2;
      Hnorm=Hnorm+Flsave(Lfr,Lsaved)*WEIGHT(Lfr)
C-         if(mod(Lfr,20)==0 .and. mod(Nstep,100)==0)then;
C-        if(Nstep==1000 .or. Nstep==1050 .or. Nstep==1100)then;
C-        if(Nstep>=350 .and. Nstep<=550)then;
C-          write(*,'(a,2i5,1p,3e12.3)')
C-            ' lbalsw:  wH(Lfr),cH(Lfr):',Lfr,Lfr0,wH(Lfr),cH(Lfr);
C-          write(*,'(a,1p,3e12.3)')
C-            ' Hobsg, Y(@L), Hnorm:     ',Hobsg, Y(@L), Hnorm;
C-        endif;
      else
C- Lfr0<1
      Flsave(Lfr,Lsaved)=Hobs(Lfr)*Ryp(NZON-1)**2
      endif
09872 CONTINUE
      Obflum=CLUMF*Hnorm
C- obs.lum. according Flsave
      if(mod(Nstep,100).EQ.0)then
C-       write(*,'(a,1p,4e12.3)')'  lbalsw: Obslum Obflum CLUMF Hnorm: ',
C-       Obslum, Obflum, CLUMF, Hnorm;
C-       pause;
      endif
C-    Obflum=max(CLUMF*Hnorm,1.d33); -- obs.lum. according Flsave
      doL=1,NFRUS
      if(Obflum.GT.1.d-100)then
      Flsave(L,Lsaved)=Flsave(L,Lsaved)*Obslum/Obflum
      else
      Flsave(L,Lsaved)=FJnois
      endif
      enddo
C-    Flsave(Mfreq,Lsaved)=TPHfit; -- save Catchpole Tphotosphere
      Flsave(Mfreq,Lsaved)=Rphtau
C- save Rph from tau=0.64
      callubvnew
      If(MOD(NSTEP,MAX0(1,NOUT/10)).EQ.0)then
      WRITE(4,'(1X,1P,G12.5,0P,F11.3,3F12.3,1P,E10.2,0P,6F10.3)')Time,TEFF*1.D-03,RPH,VPH*1.D+06/(UTIME*CRAP),XPH,HUSED*UTIME*CRAP,M
     *BOL,MU,MB,MV,MU-MB,MB-MV
      endif
C-   Write(11,'(8A4)') SNGL(Time),SNGL(TEFF*1.D-03),
C-      SNGL(RPH),SNGL(VPH*1.D+06/(UTIME*CRAP)),SNGL(MBOL),
C-    SNGL(MU),SNGL(MB),SNGL(MV);
      Tcurv(1,Lsaved)=Time
      Tcurv(2,Lsaved)=TEFF*1.D-03
C- Our Tphotosphere
      Tcurv(3,Lsaved)=RPH
      Tcurv(4,Lsaved)=VPH*1.D+06/(UTIME*CRAP)
      Tcurv(5,Lsaved)=MBOL
      Tcurv(6,Lsaved)=MU
      Tcurv(7,Lsaved)=MB
      Tcurv(8,Lsaved)=MV
      depos(Lsaved)=ELVOL*UEPRI/UTIME
C-    depos(Lsaved)=TPHfit;  -- for experiments
      RETURN
      END
