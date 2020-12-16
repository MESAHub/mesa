C- GEAR integrator for ODE with sparse Jacobian and option METH=3
C- for BRAYTON, GUSTAVSON & HACHTEL method
      SUBROUTINEVTSTIF
      IMPLICITREAL*8(A-H,O-Z)
C VTSTIF performs one step of the integration of an initial
C   value problem for a system of ordinary differential equations.
C   Communication with VTSTIF is done with the following variables:
C   Y       an NYDIM by LMAX array containing the dependent variables
C             and their scaled derivatives. LMAX is currently 13 for
C             the Adams methods and 6 for the Gear methods (though
C             LMAX=5 is recommended due to stability considerations).
C             LMAX-1 is MAXORD, the maximum order used. See subroutine
C             COSET.
C               Y(I,J+1) contains the J-th derivative of Y(I), scaled
C             by H**J/Factorial(J). Only Y(I), 1 <= I <= N, need be set
C             by the calling program on the first entry.
C               If it is desired to interpolate to non-mesh points,
C             the Y array can be used. If the current step size is H
C             and the value at t+E is needed, form  S = E/H, and then
C             compute
C                            NQ
C               Y(I)(t+E) = SUM   Y(I,J+1)*S**J .
C                           J=0
C             The Y array should not be altered by the calling program.
C             When referencing Y as a 2-dimensional array, use a column
C             length of NYDIM, as this is the value used in STIFF.
C   N       The number of first order differential equations. N may be
C             decreased on later calls if the number of active equations
C             reduces, but it must not be increased without calling with
C             JSTART=0.
C   NYDIM   A constant integer >= N, used for dimensioning purposes.
C             NYDIM must not be changed without setting JSTART=0.
C   t       The independent variable. t is updated on each step taken.
C   H       The step size to be attempted on the next step. H may be
C             adjusted up or down by the routine in order to achieve
C             an economical integration. However, if the H provided by
C             the user does not cause a larger error than requested, it
C             will be used. To save computer time, the user is advised
C             to use a fairly small step for the first call. It will be
C             automatically increased later. H can be either positive
C             or negative, but its sign must remain constant throughout
C             the problem.
C   HMIN    The minimum absolute value of the step size that will be
C             used for the problem. On starting this must be much smaller
C             than the average Abs(H) expected, since a first order
C             method is used initially.
C   HMAX    The maximum absolute value of the step size that will be
C             used for the problem.
C   EPSJ    The relative error test constant. Single step error estimates
C             divided by YMAX(I) must be less than this in the Euclidean
C             norm. The step and/or order is adjusted to achieve this.
C   METH    The method flag.
C             METH=1  means the implicit ADAMS methods
C             METH=2  means the GEAR method for stiff problems
C             METH=3  means the GEAR method, corrected by
C             R.K.Brayton, F.G.Gustavson & G.D.Hachtel in the
C             Proceedings of the IEEE v.60, No.1, January 1972, p.98-108
C   YMAX    An array of N locations which contains the maximum absolute
C             value of each Y seen so far.
C   ERROR   An array of N elements proportional to the estimated one step
C             error in each component.
C   KFLAG   A completion code with the following meanings:
C             0 the step was succesful
C            -1 the requested error could not be achieved with
C               abs(H)=HMIN
C            -2 corrector convergence could not be achieved for
C               abs(H)>HMIN.
C             On a return with KFLAG negative, the values of t and the Y
C             array are as of the beginning of the last step, and H is
C             the last step size attempted.
C   JSTART  An integer used on input and output.
C             On input it has the following values and meanings:
C              0 perform the first step. This value enables the
C                subroutine to initialize itself.
C             >0 take a new step continuing from the last.
C                Assumes the last step was succesful and user has not
C                changed any parameters.
C             <0 repeat the last step with a new value of H and/or
C                EPS and/or METH. This may be either in redoing a step
C                that failed, or in continuing from a succesful step.
C             On exit, JSTART is set to NQ, the current order of the
C                method. This is also the order of the maximum derivative
C                available in the Y array. After a succesful step, JSTART
C                need not be reset for the next call.
C   AJAC    A block used for partial derivatives. It keeps the
C             matrix AJAC=1-EL(1)*H*Jacobian & is computed by DIFMAT.
C   HL      = H*EL(1) is communicated to DIFMAT to compute AJAC.
C   FSAVE   A block of at least 2*NYDIM locations for temporary storage.
C
C The parameters which must be input by the user are:
C   N, NYDIM, t, Y, H, HMIN, HMAX, EPS, METH, JSTART, MAXORD.
C
C Additional Subroutines required are:
C  DIFMAT computes DY/Dt given t and Y, and if EVALJA computes AJAC,
C         stored in the form appropriate for M28Y12,
C  COSET,
C  RESCAL,
C  M28Y12, M30Y12, M28BYS - sparse matrix solvers.
C
C The calling program must contain the common declarations given in
C node %C:
C                                                                      *
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
C: COMMONS & VAR STIFF *
      common/NSTEP/NSTEP,NDebug,MAXER,IOUT,
C- IOUT      LINE PRINT INTERVAL
     *NOUT
C- NOUT      PRINT STEP
      common/CREAD/TAUOLD,NSTMAX,MBATCH,MAXORD
      common/debug/LfrDebug,Nperturb,Kbad
C- COMMON/NSTEP/NSTEP,NDebug,MAXER,IOUT,NOUT;
C- COMMON/CREAD/TAUOLD,NSTMAX,MBATCH,MAXORD;
      COMMON/STSAVE/RMAX,TREND,OLDL0,RC,HOLD,EDN,E,EUP,BND,EPSOLD,TOLD,MEO,NOLD,NQ,LNQ,IDOUB
C- COMMON/STSAVE/RMAX,TREND,OLDL0,RC,HOLD,EDN,E,EUP,BND,EPSOLD,TOLD,
C-           MEO,NOLD,NQ,LNQ,IDOUB;
      REAL*8ERROR(NYDIM),EL(13),TQ(4)
      REAL*8XMOD,ERREST,ERROLD
      PARAMETER(NCHANL=6,ETA=1.D-04)
C- PIVOT THRESHOLD IN M28BYS
      PARAMETER(DeltaDeb=1.d-10)
C- FOR DEBUGGER
      CHARACTERSTR020*20,LIN130*130
C- FOR DEBUGGER
      CHARACTER*20CharVarName(5)
C- FOR DEBUGGER
      LogicalFOUNDI
C- FOR DEBUGGER
      REAL*8DNDVAR(NYDIM),DSAVE(NYDIM)
C- FOR DEBUGGER
      IntegerRETCOD/0/,Ifail/0/
C- RETURN CODE
      RealTm1,Tm2
C- Timer for CMS -- CRAY REAL
      LOGICALCORRCO
C- CORRECTOR ITERATIONS HAVE CONVERGED
      LOGICALNUMSING
C- Numerical singularity in decomposition
C: COMMON & VARIABLES FOR LINEAR EQS SOLVER F01BRF or M28Y12 *
      Parameter(MAXIT=15)
      LOGICALLBLOCK,GROW,ABORT(4)
      COMMON/JSAVE/ASAVE(NZ),IRS(NZ),ICS(NZ)
      DimensionDB(NYDIM),XSAVE(NYDIM)
C: DATA FOR STIFF *
      DATAANOISE/1.D-12/
C-should be set to the noise level of the machine
      DATADFLTZR/1.D-75/,MAXITE/3/,MAXFAI/3/
      DATARMXINI/1.D+4/,RMXNOR/10.D0/,RMXFAI/2.0D0/,IDELAY/10/
      DATARHCORR/.25D0/,RHERR3/.1D0/,RCTEST/.3D0/
      DATABIAS1/1.3D0/,BIAS2/1.2D0/,BIAS3/1.4D0/
      DATACharVarName/'radius','velocity','temperature','FJ','FH'/
C: DATA FOR SPARSE MATRIX SubroutineS M28Y12 *
      DATALBLOCK/.True./,GROW/.TRUE./,ABORT/.FALSE.,.TRUE.,.FALSE.,.TRUE./,PIVOT/1.D-01/
C: PREPARATIONS DEPENDING ON JSTART *
C------ '-->Entering Node %F:'
      KFLAG=0
      tOLD=t
      NUMSING=.FALSE.
C-  WRITE(@term,*)'Entering Stiff NSTEP=', NSTEP;
C  On the first call, the order is set to 1 and the initial
C       derivatives are calculated. RMAX is the maximum ratio
C       by which H can be increased in a single step. It is
C       initially RMXINI to compensate for the small initial H,
C       but then is normally equal to RMXNOR. If a failure occurs
C       (in corrector convergence or error test), RMAX is set
C       at RMXFAI for the next increase. *
      IF(.NOT.(JSTART.EQ.0))GOTO09999
C: INITIAL *
C------ '-->Entering Node %FB:'
      HUSED=0.D0
      NQUSED=0
      NOLD=N
      EVALJA=.FALSE.
      NQ=1
      DO09996I=1,N
      FSAVE(I)=Y(I,1)
09996 CONTINUE
      CALLCOSET(METH,NQ,EL,TQ,MAXORD,IDOUB)
C- needed to define EL(1)
      HL=H*EL(1)
      CALLDIFMAT
      DO09993I=1,N
      Y(I,2)=FSAVE(NYDIM+I)*H
09993 CONTINUE
      LNQ=2
      IDOUB=LNQ+1
      RMAX=RMXINI
      EPSJ=EPS*SQRT(Real(N))
      TREND=1.D0
      OLDL0=1.D0
      RC=0.D0
C- RC is the ratio of new to old values
C- of the coefficient H*EL(1).
C- When RC differs from 1 by more than RCTEST,
C- EVALJA is set to .TRUE. to force
C- Jacobian to be calculated in DIFMAT.
      HOLD=H
      MEO=METH
      EVALJA=.TRUE.
      OLDJAC=.FALSE.
      NEEDBR=.TRUE.
      LMAX=MAXORD+1
      EPSOLD=EPSJ
C------ '<--Leaving  Node %FB:'
      GOTO09997
09999 CONTINUE
      IF(.NOT.(JSTART.LT.0))GOTO09998
C: reset Y & other variables for changed
C                            METH, EPS, N, H *
C------ '-->Entering Node %FC:'
C- If the caller change EPSJ or METH , the constants E,EUP,EDN,BND
C- must be reseted.
      EPSJ=EPS*SQRT(Real(N))
      IF(METH.EQ.MEO)THEN
      IF(EPSJ.NE.EPSOLD)THEN
      EPSOLD=EPSJ
      ENDIF
      ELSE
      IDOUB=LNQ+1
      ENDIF
C: TEST N==NOLD, H==HOLD ? *
C------ '-->Entering Node %FCC:'
      IF(N.NE.NOLD)THEN
      IDOUB=LNQ+1
      EVALJA=.TRUE.
      OLDJAC=.FALSE.
      NOLD=N
      ENDIF
      IF(H.NE.HOLD)THEN
      RH=H/HOLD
      H=HOLD
      CALLRESCAL(Y,NYDIM,RH,RMAX,RC,LNQ)
      IDOUB=LNQ+1
      ENDIF
C------ '<--Leaving  Node %FCC:'
      EPSOLD=EPSJ
C------ '<--Leaving  Node %FC:'
      GOTO09997
09998 CONTINUE
09997 CONTINUE
C------ '<--Leaving  Node %F:'
      HNT(1)=H
      CALLCOSET(METH,NQ,EL,TQ,MAXORD,IDOUB)
C- needed for continuation
C-runs
      RC=RC*EL(1)/OLDL0
      OLDL0=EL(1)
C- RC is the ratio of new to old values
C- of the coefficient H*EL(1).
C- When RC differs from 1 by more than RCTEST,
C- EVALJA is set to .TRUE. to force
C- jacobian to be calculated in DIFMAT.
C: the constants E,EUP,EDN,BND *
C------ '-->Entering Node %A:'
      EDN=(TQ(1)*EPSJ)**2
C- is to test for decreasing the order,
      E=(TQ(2)*EPSJ)**2
C- comparison for errors of current order NQ,
      EUP=(TQ(3)*EPSJ)**2
C- is to test for increasing the order,
      BND=(TQ(4)*EPSJ)**2
C- is used to test for convergence of the
C- correction iterates
C------ '<--Leaving  Node %A:'
      IF(MOD(NSTEP,MBATCH).EQ.0)NEEDBR=.TRUE.
      RETCOD=0
C: GENERAL STEP OF STIFF WHILE RETCOD==0 *
C------ '-->Entering Node %G:'
09990 IF(.NOT.(RETCOD.EQ.0))GOTO09989
      IF(ABS(RC-1.D0).GT.RCTEST.or.MOD(NSTEP,MBATCH).EQ.0)EVALJA=.TRUE.
      If(ifail.EQ.13)Then
      Needbr=.true.
      Evalja=.true.
      Endif
      t=t+h
C: Predictor.
C        This section computes the predicted values by effectively
C        multiplying the Y array by the Pascal triangle matrix *
C------ '-->Entering Node %GE:'
      DO09987J1=1,NQ
      DO09984J2=J1,NQ
      J=NQ-J2+J1
      DO09981I=1,N
      Y(I,J)=Y(I,J)+Y(I,J+1)
09981 CONTINUE
09984 CONTINUE
09987 CONTINUE
C------ '<--Leaving  Node %GE:'
09978 CONTINUE
C- loop until .not.EVALJA
      ITER=0
      DO09975I=1,N
      ERROR(I)=0.D0
09975 CONTINUE
      DO09972I=1,N
      FSAVE(I)=Y(I,1)
09972 CONTINUE
      CORRCO=.FALSE.
C- Corrector convergence. It is tested
C- by requiring changes relative to YMAX(I)
C- to be less, in Euclidian norm, than BND,
C- which is dependent on EPSJ.
09969 IF(.NOT.(.NOT.CORRCO.AND.ITER.LT.MAXITE))GOTO09968
C: Up to MAXITE corrector iterations are taken.
C                 The sum of the corrections is accumulated in the
C                 vector ERROR(i). It is approximately equal to the
C                 Lnq-th derivative of Y multiplied by
C                    H**Lnq/(factorial(Lnq-1)*EL(Lnq)),
C                 and is thus proportional to the actual errors
C                 to the lowest power of H present (H**Lnq).
C                 The Y array is not altered in the correction
C                 loop. The updated Y vector is stored temporarily
C                 in FSAVE. The norm of the iterate difference
C                 is stored in D. *
C------ '-->Entering Node %G_Cor:'
      BADSTE=.FALSE.
C- Logical variable BADSTE may be changed by DIFMAT
      HL=H*EL(1)
C- if(Nstep>1750 & EVALJA)then;
C- if(Nstep>1586 & EVALJA)then;
      if(Nstep.GE.NDebug.AND.EVALJA)then
C: full debugger of Jacobian *
C------ '-->Entering Node %G_COR_DEB:'
      If(.NOT.BADSTE)Then
      t=tOLD
      DO09966J1=1,NQ
      DO09963J2=J1,NQ
      J=NQ-J2+J1
      DO09960I=1,N
      Y(I,J)=Y(I,J)-Y(I,J+1)
09960 CONTINUE
09963 CONTINUE
09966 CONTINUE
      DO09957Ich=1,Mzon
      Chem0(Ich)=0.1
C- check if we need this
09957 CONTINUE
      DO09954I=1,N
      FSAVE(I)=Y(I,1)
09954 CONTINUE
C-************** ************** ************** **************
C-   EVALJA=.FALSE.; -- never put it here: some errors may be lost
C-************** ************** ************** **************
      Nperturb=0
C- index of perturbation, now NO perturbation
      CALLDIFMAT
      DO09951KTST=1,N
C- save old YDOT
      DSAVE(KTST)=FSAVE(NYDIM+KTST)
09951 CONTINUE
      DO09948JTST=1,N
C-@VAR
C-     _DO JTST=Ncnd+1,Ncnd+2; -- only d over dr is tested
C-      @wterm' JTST=',JTST;
      EVALJA=.FALSE.
      OLD=FSAVE(JTST)
      If(JTST.LT.Nzon)then
C- radius
      If(JTST.EQ.1)then
      DeltY=DeltaDeb*min(FSAVE(JTST+1)-OLD,OLD-Rce)
      else
      DeltY=DeltaDeb*min(FSAVE(JTST+1)-OLD,OLD-FSAVE(JTST-1))
      endif
      elseIf(Nzon.LT.JTST.AND.JTST.LT.2*Nzon)then
C- velocity
      If(JTST.EQ.Nzon+1)then
      DeltY=DeltaDeb*min(abs(FSAVE(JTST+1)-OLD),abs(OLD))
      else
      DeltY=DeltaDeb*min(abs(FSAVE(JTST+1)-OLD),abs(OLD-FSAVE(JTST-1)))
      endif
      else
      DeltY=OLD*DeltaDeb
      endif
      FSAVE(JTST)=OLD+DeltY
      DO09945Ich=1,Mzon
      Chem0(Ich)=0.1
09945 CONTINUE
      Nperturb=1
C- index of perturbation, now  IS perturbation
      CALLDIFMAT
      FSAVE(JTST)=OLD
C: put in IRN(NZ+1:ILAST) all indeces I such that
C                ICN(I)==JTST *
C------ '-->Entering Node %G_COR_DEB_ROW:'
      ILAST=0
      DO09942I=1,NZMOD
      IF(ICS(I).EQ.JTST)THEN
      ILAST=ILAST+1
      IRN(NZ+ILAST)=I
      ENDIF
09942 CONTINUE
C-  WRITE(5,*)'ILAST=',ILAST;
C------ '<--Leaving  Node %G_COR_DEB_ROW:'
      DO09939KTST=1,N
C-@NAMELEFT,@NAMERIGHT;
      DNDVAR(KTST)=(FSAVE(NYDIM+KTST)-DSAVE(KTST))
C: IF in IRN(NZ+1:ILAST) there is I such that
C                     IRN(I)==KTST THEN
C                       IF AJAC(I) differs strongly
C                          from -HL*@DNAMEDVAR for KTST^=JTST or
C                          from 1.-HL*@DNAMEDVAR for KTST==JTST
C                       THEN output relevant vars
C                     ELSE
C                       IF @DNAMEDVAR is HIGHER than @NOISE
C                       THEN output relevant vars          *
C------ '-->Entering Node %G_COR_DEB_COMPR:'
      FOUNDI=.FALSE.
      IL=1
09936 IF(.NOT.(IL.LE.ILAST.AND..NOT.FOUNDI))GOTO09935
      IF(IRS(IRN(NZ+IL)).EQ.KTST)THEN
      FOUNDI=.TRUE.
      ELSE
      IL=IL+1
      ENDIF
      GOTO09936
09935 CONTINUE
      IF(FOUNDI)THEN
      If(KTST.NE.JTST)THEN
      If(ABS(AJAC(IRN(NZ+IL))*DeltY+HL*DNDVAR(KTST)).GT.1.D-1*ABS(AJAC(IRN(NZ+IL))*DeltY).AND.ABS(AJAC(IRN(NZ+IL))/HL).GT.1.D-4)Then
      If(abs(DeltY).GT.DFLTZR*abs(DNDVAR(KTST)))then
      WRITE(4,'(/A/2I7,1P,3G12.3)')'  K,J,AJAC,Stella df/dY,Num DN(K)/DV(J):   ',KTST,JTST,AJAC(IRN(NZ+IL)),-AJAC(IRN(NZ+IL))/HL,DND
     *VAR(KTST)/DeltY
      else
      WRITE(4,'(/A/2I7,1P,3G12.3)')'  K,J,AJAC,(df/dY,DN/DV)*DeltaDeb:  ',KTST,JTST,AJAC(IRN(NZ+IL)),-AJAC(IRN(NZ+IL))*DeltY/HL,DNDV
     *AR(KTST)
      endif
C: give info on variable number KTST *
C------ '-->Entering Node %G_COR_DEB_COMPR_PhysInfo:'
      	doI=1,NZMOD
      	if(IRS(I).EQ.KTST.and.ICS(I).EQ.JTST)then
      	WRITE(4,'(3(A,1P,E12.4))')	'  Y(K)=',Y(KTST,1),'    Fsave(K)=',Fsave(KTST),	'    Ydot(K)=',FSAVE(NYDIM+KTST)
      WRITE(4,'(2(A,I5),A,1P,E12.4)')' Control IR=',KTST,	'  IC=',JTST,'   AJAC=',AJAC(I)
      endif
      	enddo
C- this is for control of found AJAC;
      ICW=JTST
      If(ICW.LE.NZON*NVARS)then
      Izon=MOD(ICW-1,NZON)+1
      Ivar=(ICW-1)/NZON+1
      WRITE(4,'(3A,I3)')' Derivative over ',CharVarName(Ivar)(1:length(CharVarName(Ivar))),' in  zon=',Izon
      else
      ICW=ICW-NZON*NVARS
      If(ICW.LE.KRAD)Then
      L=(ICW-1)/(NZON-NCND)+1
      Izon=NCND+MOD(ICW-1,(NZON-NCND))+1
      WRITE(4,'(2(A,I3))')'  over  FJ in L=',L,'    Izon=',Izon
      else
      ICW=ICW-KRAD
      L=(ICW-1)/(NZON-NCND)+1
      Izon=NCND+MOD(ICW-1,(NZON-NCND))+1
      WRITE(4,'(2(A,I3))')'  over  FH in L=',L,'    Izon=',Izon
      endif
      endif
      If(KTST.LE.NZON*NVARS)then
      Izon=MOD(KTST-1,NZON)+1
      Ivar=(KTST-1)/NZON+1
      WRITE(4,'(3A,I3)')' for Ydot of ',CharVarName(Ivar)(1:length(CharVarName(Ivar))),' in  zon=',Izon
      else
      KTST1=KTST-NZON*NVARS
      If(KTST1.LE.KRAD)Then
      L=(KTST1-1)/(NZON-NCND)+1
      Izon=NCND+MOD(KTST1-1,(NZON-NCND))+1
      WRITE(4,'(2(A,I5))')'  Der. of dot(FJ) in L=',L,'    Izon=',Izon
      else
      KTST1=KTST1-KRAD
      L=(KTST1-1)/(NZON-NCND)+1
      Izon=NCND+MOD(KTST1-1,(NZON-NCND))+1
      WRITE(4,'(2(A,I5))')'  Der. of dot(FH) in L=',L,'    Izon=',Izon
      endif
      endif
C------ '<--Leaving  Node %G_COR_DEB_COMPR_PhysInfo:'
      endif
      else
      If(ABS((AJAC(IRN(NZ+IL))-1.D0)*DeltY+HL*DNDVAR(KTST)).GT.1.D-1*ABS((AJAC(IRN(NZ+IL))-1.D0)*DeltY).AND.ABS((AJAC(IRN(NZ+IL))-1.
     *D0)/HL).GT.1.D-4)Then
      If(abs(DeltY).GT.DFLTZR*abs(DNDVAR(KTST)))then
      WRITE(4,'(/A/2I7,1P,3G12.3)')'  K,J,AJAC,Stella df/dY,Num DN(K)/DV(J):  ',KTST,JTST,AJAC(IRN(NZ+IL)),(1.D0-AJAC(IRN(NZ+IL)))/H
     *L,DNDVAR(KTST)/DeltY
      else
      WRITE(4,'(/A/2I7,1P,3G12.3)')'  K,J,AJAC,(df/dY,DN/DV)*DeltaDeb:  ',KTST,JTST,AJAC(IRN(NZ+IL)),(1.D0-AJAC(IRN(NZ+IL)))*DeltY/H
     *L,DNDVAR(KTST)
      endif
C: give info on variable number KTST *
C------ '-->Entering Node %G_COR_DEB_COMPR_PhysInfo:'
      	doI=1,NZMOD
      	if(IRS(I).EQ.KTST.and.ICS(I).EQ.JTST)then
      	WRITE(4,'(3(A,1P,E12.4))')	'  Y(K)=',Y(KTST,1),'    Fsave(K)=',Fsave(KTST),	'    Ydot(K)=',FSAVE(NYDIM+KTST)
      WRITE(4,'(2(A,I5),A,1P,E12.4)')' Control IR=',KTST,	'  IC=',JTST,'   AJAC=',AJAC(I)
      endif
      	enddo
C- this is for control of found AJAC;
      ICW=JTST
      If(ICW.LE.NZON*NVARS)then
      Izon=MOD(ICW-1,NZON)+1
      Ivar=(ICW-1)/NZON+1
      WRITE(4,'(3A,I3)')' Derivative over ',CharVarName(Ivar)(1:length(CharVarName(Ivar))),' in  zon=',Izon
      else
      ICW=ICW-NZON*NVARS
      If(ICW.LE.KRAD)Then
      L=(ICW-1)/(NZON-NCND)+1
      Izon=NCND+MOD(ICW-1,(NZON-NCND))+1
      WRITE(4,'(2(A,I3))')'  over  FJ in L=',L,'    Izon=',Izon
      else
      ICW=ICW-KRAD
      L=(ICW-1)/(NZON-NCND)+1
      Izon=NCND+MOD(ICW-1,(NZON-NCND))+1
      WRITE(4,'(2(A,I3))')'  over  FH in L=',L,'    Izon=',Izon
      endif
      endif
      If(KTST.LE.NZON*NVARS)then
      Izon=MOD(KTST-1,NZON)+1
      Ivar=(KTST-1)/NZON+1
      WRITE(4,'(3A,I3)')' for Ydot of ',CharVarName(Ivar)(1:length(CharVarName(Ivar))),' in  zon=',Izon
      else
      KTST1=KTST-NZON*NVARS
      If(KTST1.LE.KRAD)Then
      L=(KTST1-1)/(NZON-NCND)+1
      Izon=NCND+MOD(KTST1-1,(NZON-NCND))+1
      WRITE(4,'(2(A,I5))')'  Der. of dot(FJ) in L=',L,'    Izon=',Izon
      else
      KTST1=KTST1-KRAD
      L=(KTST1-1)/(NZON-NCND)+1
      Izon=NCND+MOD(KTST1-1,(NZON-NCND))+1
      WRITE(4,'(2(A,I5))')'  Der. of dot(FH) in L=',L,'    Izon=',Izon
      endif
      endif
C------ '<--Leaving  Node %G_COR_DEB_COMPR_PhysInfo:'
      endif
      endif
      ELSE
      IF(ABS(DNDVAR(KTST)).GT.1.D-4)then
      If(abs(DeltY).GT.DFLTZR*abs(DNDVAR(KTST)))then
      WRITE(4,'(/A,2I7,1P,2G12.3)')' NO JAC: K,J,DN(K)/DV(J)  ',KTST,JTST,DNDVAR(KTST)/DeltY
      else
      WRITE(4,'(/A,2I7,1P,2G12.3)')' NO JAC: K,J,DN(K)/DV(J)*DeltaDeb ',KTST,JTST,DNDVAR(KTST)
      endif
      endif
      ENDIF
C------ '<--Leaving  Node %G_COR_DEB_COMPR:'
09939 CONTINUE
09948 CONTINUE
      STOP
      else
      WRITE(4,*)' Testing not possible: BADSTE=',BADSTE
      endif
C------ '<--Leaving  Node %G_COR_DEB:'
      endif
      CALLDIFMAT
C- if necessary (i.e. EVALJA=.true.) the partials
C- are reevaluated by DIFMAT
CMatrOut : output Jacobian matrix *
C- stop ' Jacobian is output';
      IF(.NOT.(BADSTE))GOTO09933
C- Something is wrong in DIFMAT (i.e. Y(I) is out of physical
C- range for some I)
C-       WRITE(*,*)' G_Cor: BADSTE=',BADSTE;
      ITER=MAXITE
      EVALJA=.FALSE.
      OLDJAC=.FALSE.
      Ifail=13
C- Variables are given the values such as to guarantee that
C- control passes to the node %GY - to try a new smaller step
      GOTO09932
09933 CONTINUE
      DO09930Isw=1,N
      XSAVE(Isw)=0.D0
09930 CONTINUE
      IF(EVALJA)THEN
C: L-U Decomposition of matrix AJAC=1-El*H*Jacobian *
C------ '-->Entering Node %G_Cor.L:'
C_Do I=1,NZMOD;
C      If(IRS(I)==201 & Nstep>=160 & Nstep<=180) then;
C        WRITE(NCHANL,'(2(A,I5),1P,A,E15.8,A,L3)')
C              ' ROW=',IRS(I),'   COL=',ICS(I),'   AJAC=',AJAC(I),
C              '   Needbr=',Needbr;
C      endif
C  _od; *
C-  write(Nchanl,*)' Nzmod=',Nzmod;
      DO09927Isw=1,NZMOD
      ASAVE(Isw)=AJAC(Isw)
09927 CONTINUE
C-_repeat
      If(Needbr)then
C:   *
C------ '-->Entering Node %G_Cor.L_SOLY12:'
C-    IFAIL=110; -- hard failure
      IFAIL=111
C- soft failure
C     if(NUMSING)then;
C       NZMOD=N;
C       do isw=1,N;
C         IRN(isw)=isw;
C         ICN(isw)=isw;
C         AJAC(isw)=1.d0;
C       enddo;
C     else;
C *
      DO09924Isw=1,NZMOD
      IRN(Isw)=IRS(Isw)
      ICN(Isw)=ICS(Isw)
09924 CONTINUE
C    endif; *
C-    Call VTIME(Tm1);
C-    CALL F01BRF(N , NZMOD, AJAC, LICN, IRN, LIRN, ICN, PIVOT,
C-                IKEEP, IW, WJAC, LBLOCK, GROW, ABORT, IDISP, IFAIL);
C-      write(*,*) ' In stiffbgh.trf  NZMOD=',NZMOD;
C-      pause;
      CALLM28Y12(N,NZMOD,AJAC,LICN,IRN,LIRN,ICN,PIVOT,IKEEP,IW,WJAC,LBLOCK,GROW,ABORT,IDISP,IFAIL,1.D-12)
      NEEDBR=.FALSE.
C- F01BRF or M28Y12 NOT NEEDED
C-    Call VTIME(Tm2);
C-      IF (GROW) WRITE (NCHANL,89996) WJAC(1);
89996 FORMAT(' ON EXIT FROM M28Y12:  W(1)  = ',1P,E12.4)
C-      WRITE (NCHANL,99996)IDISP(2), IDISP(6),
C-             IDISP(7), IDISP(3),IDISP(4);
C-      IF (LBLOCK) WRITE (NCHANL,99994) (IDISP(I),I=8,10);
99999 FORMAT(6A4,A3)
99998 FORMAT(4(1X/),1X,5A4,A3,'RESULTS'/1X)
99997 FORMAT(5(2X,G8.0,2X,I1,2X,I1))
99996 FORMAT(' NUMBER OF NON-ZEROS IN DECOMPOSITION',9X,'=',I9/' MINIMUM SIZE OF ARRAY IRN',20X,'=',I9/' MINIMUM SIZE OF ARRAYS A AN
     *D ICN',13X,'=',I9/' NUMBER OF COMPRESSES ON IRN (IDISP(3))',7X,'=',I7/' NUMBER OF COMPRESSES ON A AND ICN (IDISP(4)) =',I7)
99994 FORMAT(' STRUCTURAL RANK',16X,'=',I7/' NUMBER OF DIAGONAL BLOCKS',6X,'=',I7/' SIZE OF LARGEST DIAGONAL BLOCK =',I7)
C- WRITE(NCHANL,
C- '('' M28Y12-1 done at step :'',I6,'' Time(mS):'',I6)')
C- CRAY  '('' M28Y12 done at step :'',I6,'' Time(S):'',1P,G11.5)')
C-       NSTEP,Tm2-Tm1;
C------ '<--Leaving  Node %G_Cor.L_SOLY12:'
      IRW=Idisp(11)
      IF(IFAIL.EQ.2)then
      NUMSING=.true.
      NEEDBR=.true.
      ITER=MAXITE
      EVALJA=.FALSE.
      OLDJAC=.FALSE.
      Ifail=13
C- Variables are given the values such as to guarantee that
C- control passes to the node %GY - to try a new smaller step
      else
      NUMSING=.false.
      endif
      else
C:   *
C------ '-->Entering Node %G_Cor.L_SOLBSF:'
      IFAIL=11
C-   Call VTIME(Tm1);
      CALLM28BYS(N,NZMOD,AJAC,LICN,IRS,ICS,ICN,IKEEP,IW,WJAC,GROW,ETA,RPMIN,ABORT(4),IDISP,IFAIL)
C-   Call VTIME(Tm2);
C-  WRITE(NCHANL,'('' M28BYS done at step :'',I6,'' Time(ms):'',I6)')
C-      NSTEP,Tm2-Tm1;
C   IF (GROW) WRITE (NCHANL,89995) WJAC(1);
C      WRITE (NCHANL,89994) RPMIN;
C89995:FORMAT (' ON EXIT FROM M28BYS:  W(1)  = ',1P,E12.4);
C89994:FORMAT (' VALUE OF RPMIN = ', G12.4);                *
      IF(IFAIL.NE.0.OR.WJAC(1).GT.1.D50)THEN
      IF(IFAIL.EQ.6)THEN
C: process the numerical singularity in row IR *
C------ '-->Entering Node %G_Cor.L_SOLBSF_NUMSING:'
C------ '<--Leaving  Node %G_Cor.L_SOLBSF_NUMSING:'
      ENDIF
C-    IFAIL=110; -- hard failure
      IFAIL=111
C- soft failure
      DO09921Isw=1,NZMOD
      IRN(Isw)=IRS(Isw)
      ICN(Isw)=ICS(Isw)
      AJAC(Isw)=ASAVE(Isw)
09921 CONTINUE
C-      Call VTIME(Tm1);
      CALLM28Y12(N,NZMOD,AJAC,LICN,IRN,LIRN,ICN,PIVOT,IKEEP,IW,WJAC,LBLOCK,GROW,ABORT,IDISP,IFAIL,1.D-12)
      NEEDBR=.FALSE.
C- F01BRF or M28Y12 NOT NEEDED
C-      Call VTIME(Tm2);
C-      WRITE(NCHANL,
C-        '('' M28Y12 done again at step :'',I6,'' Time(ms):'',I6)')
C- CRAY   '('' M28Y12 done at step :'',I6,'' Time(S):'',1P,G11.5)')
C-          NSTEP,Tm2-Tm1;
      ENDIF
C------ '<--Leaving  Node %G_Cor.L_SOLBSF:'
      IRW=IW(1)
      IF(IFAIL.EQ.6)then
      NUMSING=.true.
      NEEDBR=.true.
      ITER=MAXITE
      EVALJA=.FALSE.
      OLDJAC=.FALSE.
      Ifail=13
C- Variables are given the values such as to guarantee that
C- control passes to the node %GY - to try a new smaller step
      else
      NUMSING=.false.
      NEEDBR=.false.
      endif
      endif
      IF(NUMSING)THEN
C: process the numerical singularity in row IR *
C------ '-->Entering Node %G_Cor.L_NUMSING:'
      DO09918I=1,NZMOD
      IF(IRS(I).EQ.IRW)THEN
      WRITE(*,'(2(A,I5),A,1P,E12.4)')'IR=',IRS(I),'  IC=',ICS(I),'   AJAC=',ASAVE(I)
      ICW=ICS(I)
      If(ICW.LE.NZON*NVARS)then
      Izon=MOD(ICW-1,NZON)+1
      Ivar=(ICW-1)/NZON+1
      WRITE(*,'(2(A,I5))')'Ivar=',Ivar,'    Izon=',Izon
      else
      ICW=ICW-NZON*NVARS
      If(ICW.LE.KRAD)Then
      L=(ICW-1)/(NZON-NCND)+1
      Izon=NCND+MOD(ICW-1,(NZON-NCND))+1
      WRITE(*,'(2(A,I5))')'FJ in L=',L,'    Izon=',Izon
      else
      ICW=ICW-KRAD
      L=(ICW-1)/(NZON-NCND)+1
      Izon=NCND+MOD(ICW-1,(NZON-NCND))+1
      WRITE(*,'(2(A,I5))')'FH in L=',L,'    Izon=',Izon
      endif
      endif
      ENDIF
09918 CONTINUE
      WRITE(*,'(A)')'For Ydot of: '
      If(IRW.LE.NZON*NVARS)then
      Izon=MOD(IRW-1,NZON)+1
      Ivar=(IRW-1)/NZON+1
      WRITE(*,'(2(A,I5))')'Ivar=',Ivar,'    Izon=',Izon
      else
      IRW=IRW-NZON*NVARS
      If(IRW.LE.KRAD)Then
      L=(IRW-1)/(NZON-NCND)+1
      Izon=NCND+MOD(IRW-1,(NZON-NCND))+1
      WRITE(*,'(2(A,I5))')'FJ in L=',L,'    Izon=',Izon
      else
      IRW=IRW-KRAD
      L=(IRW-1)/(NZON-NCND)+1
      Izon=NCND+MOD(IRW-1,(NZON-NCND))+1
      WRITE(*,'(2(A,I5))')'FH in L=',L,'    Izon=',Izon
      endif
      endif
C------ '<--Leaving  Node %G_Cor.L_NUMSING:'
C-      stop 26;
      ENDIF
C-_until ^NUMSING ;
C------ '<--Leaving  Node %G_Cor.L:'
      EVALJA=.FALSE.
      OLDJAC=.FALSE.
      RC=1.D0
      NJAC=NJAC+1
      ENDIF
C: Find in FSAVE(1:N) the solution  X of equation
C
C                 AJAC!X>=H*F(Y)-Y(I,2)-ERROR(I),
C
C            where AJAC is the matrix AJAC=1-El*H*Jacobian,
C            then obtain in ERROR the full correction,
C            in D - the norm of the iterate difference &
C            in FSAVE a new approximation to Y *
C------ '-->Entering Node %G_Cor.U:'
C: Find solution by M28CYN  *
C------ '-->Entering Node %G_Cor.UK:'
      DO09915I=1,N
      FSAVE(I+NYDIM)=FSAVE(I+NYDIM)*H-Y(I,2)-ERROR(I)
C- R-H SIDE IN FSAVE(NYDIM+1:NYDIM+N)
      DB(I)=FSAVE(I+NYDIM)
C- TO SAVE R-H SIDE
09915 CONTINUE
C- ITERATIVE REFINEMENT:
      ITQ=0
      If(Ifail.NE.13)Then
09912 CONTINUE
      CALLM28CYN(N,AJAC,LICN,ICN,IKEEP,FSAVE(NYDIM+1),WJAC,IW,1,IDISP,RESID,Ifail)
C- SOLUTION IN FSAVE(NYDIM+1:NYDIM+N)
      If(Ifail.NE.13)Then
      ERROLD=ERREST
      ERREST=0.D0
      XMOD=0.D0
      DO09909ID=1,N
      XSAVE(ID)=XSAVE(ID)+FSAVE(ID+NYDIM)
      XMOD=XMOD+ABS(XSAVE(ID))
C- NORMA
      ERREST=ERREST+ABS(FSAVE(ID+NYDIM))
      FSAVE(ID)=0.D0
09909 CONTINUE
      DO09906K=1,NZMOD
      FSAVE(IRS(K))=FSAVE(IRS(K))
C- Polish strategy
     *+XSAVE(ICS(K))*ASAVE(K)
09906 CONTINUE
      DO09903IL=1,N
      FSAVE(IL+NYDIM)=DB(IL)-FSAVE(IL)
09903 CONTINUE
C-  WRITE(NCHANL,'(2(A,1P,G12.3))')' ERROLD=',ERROLD,' ERREST=',ERREST;
      ITQ=ITQ+1
      Else
      DO09900Isw=1,NZMOD
      IRN(Isw)=IRS(Isw)
      ICN(Isw)=ICS(Isw)
      AJAC(Isw)=ASAVE(Isw)
09900 CONTINUE
      NEEDBR=.TRUE.
      EVALJA=.TRUE.
C-        @wterm ' G_Cor.UK: ifail 13';
      Endif
C- _Until ERREST<=ANOISE*XMOD ! ITQ > MAXIT;
      IF(.NOT.(ERREST.LE.ANOISE*XMOD.or.Ifail.EQ.13.or.(ITQ.GT.2.AND.ERREST.GT.ERROLD).or.ITQ.GT.MAXIT))GOTO09912
      IF(ERREST.GT.XMOD*.1d0*EPS)NEEDBR=.TRUE.
C- Esaulov correction
C-  IF(XMOD>EPS**3)WRITE(NCHANL,*)' RELEST=',ERREST/XMOD;
      Endif
C------ '<--Leaving  Node %G_Cor.UK:'
      If(Ifail.NE.13)Then
C: Find ERROR, D and put new Y in FSAVE *
C------ '-->Entering Node %G_Cor.UV:'
      D=0.D0
      DO09897I=1,N
      ERROR(I)=ERROR(I)+XSAVE(I)
      D=D+(XSAVE(I)/YMAX(I))**2
      FSAVE(I)=Y(I,1)+EL(1)*ERROR(I)
09897 CONTINUE
C------ '<--Leaving  Node %G_Cor.UV:'
      Endif
C------ '<--Leaving  Node %G_Cor.U:'
      If(Ifail.NE.13)Then
      IF(ITER.NE.0)TREND=MAX(.9D0*TREND,D/D1)
C- If ITER > 0 an estimate of the convergence rate constant
C- is stored in TREND, and this is used in the convergence test
      CORRCO=D*MIN(1.D0,2.D0*TREND).LE.BND
      D1=D
      ITER=ITER+1
      Else
C-WRITE(@Wres,*)'Ifail=',Ifail;
      ITER=MAXITE
      EVALJA=.FALSE.
      OLDJAC=.FALSE.
C- Variables are given the values such as to guarantee that
C- control passes to the node %GY - to try a new smaller step
      Endif
09932 CONTINUE
C- FOR BADSTE
C------ '<--Leaving  Node %G_Cor:'
      GOTO09969
09968 CONTINUE
C- CORRCO ! ITER==MAXITE
      IF(.NOT.(CORRCO))GOTO09894
C: The corrector has converged. OLDJAC is set to .true.
C               if partial derivatives were used, to signal that they
C               may need updating on subsequent steps. The error test
C               is made and control passes to node %GXA if it
C               succeeds. *
C------ '-->Entering Node %GX:'
      D=0.D0
      ERMAX=0.D0
      MAXER=0
      DO09891I=1,N
      IF((ERROR(I)/YMAX(I))**2.GT.ERMAX)THEN
      ERMAX=(ERROR(I)/YMAX(I))**2
      MAXER=I
      ENDIF
      D=D+(ERROR(I)/YMAX(I))**2
09891 CONTINUE
      OLDJAC=.TRUE.
      IF(.NOT.(D.LE.E))GOTO09888
C- After a succesful step, update the Y array and YMAX.
C- Consider changing H if IDOUB==1. Otherwise decrease
C- IDOUB by 1.
C- If a change in H is considered, an increase
C- or decrease in order by one is considered also.
      HUSED=H
      NQUSED=NQ
      KFLAG=0
      IF(METH.EQ.3)THEN
C- Save time intervals in HNT
      DO09885J=7,3,-1
      HNT(J)=HNT(J-1)+H
C- HNT(LNQ+1) saves old H for possible order increase
C- Initially HNT(2:6) is undefined here!
09885 CONTINUE
      HNT(2)=H
      ENDIF
      DO09882J=1,LNQ
      DO09879I=1,N
      Y(I,J)=Y(I,J)+EL(J)*ERROR(I)
09879 CONTINUE
09882 CONTINUE
      IF(.NOT.(IDOUB.EQ.1))GOTO09876
C: Factors PR1, PR2 and PR3 are computed, by which H
C                 could be multiplied at order NQ-1, order NQ,
C                 or order NQ+1, respectively. The largest of these
C                 is determined and the new order is chosen accordingly.
C                 If the order is to be increased, we compute one
C                 additional scaled derivative. *
C------ '-->Entering Node %GXA:'
      PR3=DFLTZR
      IF(NQ.NE.MAXORD)THEN
C: RESET PR3 *
C------ '-->Entering Node %GXA3:'
      SUM=0.D0
      IF(METH.NE.3)THEN
      DO09873I=1,N
      SUM=SUM+((ERROR(I)-Y(I,LMAX))/YMAX(I))**2
09873 CONTINUE
      PR3=1.D0/(((SUM/EUP)**(.5D0/DBLE(LNQ+1)))*BIAS3+DFLTZR)
      ELSE
      DO09870I=1,N
      SUM=SUM+((ERROR(I)*EL(1)/DBLE(NQ+2)-Y(I,LMAX)/TQ(3))/YMAX(I))**2
09870 CONTINUE
      PR3=1.D0/(((SUM/EPSJ**2)**(.5D0/DBLE(LNQ+1)))*BIAS3+DFLTZR)
      ENDIF
C------ '<--Leaving  Node %GXA3:'
      ENDIF
C: RESET PR2 *
C------ '-->Entering Node %GX2:'
      PR2=1.D0/(((D/E)**(.5D0/DBLE(LNQ)))*BIAS2+DFLTZR)
C------ '<--Leaving  Node %GX2:'
      PR1=DFLTZR
      IF(NQ.NE.1)THEN
C: RESET PR1 *
C------ '-->Entering Node %GX1:'
      SUM=0.D0
      DO09867I=1,N
      SUM=SUM+(Y(I,LNQ)/YMAX(I))**2
09867 CONTINUE
      PR1=1.D0/(((SUM/EDN)**(.5D0/DBLE(NQ)))*BIAS1+DFLTZR)
C------ '<--Leaving  Node %GX1:'
      ENDIF
      IF((PR2.GT.PR1).AND.(PR2.GT.PR3))THEN
      RH=PR2
      ELSE
C
C        IF(PR1>PR3) THEN;
C           NEWQ=NQ-1;RH=PR1;
C        ELSE;
C           NEWQ=LNQ;RH=PR3; -- PR3 IS GREATER
C           _DO I=1,N;Y(I,NEWQ+1)=ERROR(I)*EL(LNQ)/DBLE(LNQ)_OD;
C        ENDIF;
C*
C- bug  (for METH==3) in above lines corrected:
      IF(PR1.GT.PR3)THEN
      NEWQ=NQ-1
      RH=PR1
      ELSEIF(NQ.NE.MAXORD)THEN
      NEWQ=LNQ
      RH=PR3
C- PR3 IS GREATER
      DO09864I=1,N
      Y(I,NEWQ+1)=ERROR(I)*EL(LNQ)/DBLE(LNQ)
09864 CONTINUE
      ELSE
C- NQ==MAXORD
      NEWQ=NQ
      RH=PR3
C- PR3 IS GREATER
      ENDIF
C- If there is a change of order, reset NQ and LNQ.
C- In any case H is reset according to
C- RH and the Y array and the coefficients EL()
C- are rescaled. Then exit.
      NQ=NEWQ
      LNQ=NQ+1
      ENDIF
      CALLRESCAL(Y,NYDIM,RH,RMAX,RC,LNQ)
      IDOUB=LNQ+1
C------ '<--Leaving  Node %GXA:'
      RETCOD=3
C- EXIT
      GOTO09875
09876 CONTINUE
      IDOUB=IDOUB-1
      IF((IDOUB.EQ.1).AND.(NQ.NE.MAXORD))THEN
C-  Error is saved for use in a possible
C-  order increase on the next step.
      DO09861I=1,N
      Y(I,LMAX)=ERROR(I)
09861 CONTINUE
      ENDIF
      IF(METH.EQ.3)THEN
C: RESET PR2 *
C------ '-->Entering Node %GX2:'
      PR2=1.D0/(((D/E)**(.5D0/DBLE(LNQ)))*BIAS2+DFLTZR)
C------ '<--Leaving  Node %GX2:'
      RH=PR2
      CALLRESCAL(Y,NYDIM,RH,RMAX,RC,LNQ)
C- NO IDOUB HERE!!!
      ENDIF
      RETCOD=4
C- EXIT
09875 CONTINUE
      GOTO09887
09888 CONTINUE
C- The error test failed. KFLAG keeps track
C- of multiple failures. Restore t and the Y array to their
C- previous values, and prepare to try the step again.
C- Compute the optimum step size for this or one lower order.
      NFAIL=NFAIL+1
      KFLAG=KFLAG-1
      t=tOLD
      DO09858J1=1,NQ
      DO09855J2=J1,NQ
      J=NQ-J2+J1
      DO09852I=1,N
      Y(I,J)=Y(I,J)-Y(I,J+1)
09852 CONTINUE
09855 CONTINUE
09858 CONTINUE
      RMAX=RMXFAI
      IF(.NOT.(ABS(H).LE.(HMIN*1.00001D0)))GOTO09849
      RETCOD=1
      GOTO09848
09849 CONTINUE
      IF(.NOT.(KFLAG.LE.-MAXFAI))GOTO09846
C: Control reaches this section if MAXFAI or more
C                failures have occured. It is assumed that the
C                derivatives that have accumulated in the Y array
C                have errors of the wrong order.
C                Hence the first derivative is recomputed, and the order
C                is set to one. Then H is reduced by factor of RHERR3,
C                and the step is retried. *
C------ '-->Entering Node %GXZ:'
      RH=RHERR3
      RH=MAX(HMIN/ABS(H),RH)
      H=H*RH
      HL=H*EL(1)
      DO09843I=1,N
      FSAVE(I)=Y(I,1)
09843 CONTINUE
      CALLDIFMAT
      DO09840I=1,N
      Y(I,2)=FSAVE(NYDIM+I)*H
09840 CONTINUE
C-      EVALJA=.TRUE.;
      EVALJA=.False.
C- to get out the _repeat loop
      RC=0.d0
C- to force a new Jacobian
      IDOUB=IDELAY
      IF(NQ.NE.1)THEN
      NQ=1
      LNQ=2
      CALLCOSET(METH,NQ,EL,TQ,MAXORD,IDOUB)
      OLDL0=EL(1)
C: EDN,E,EUP,BND *
C------ '-->Entering Node %A:'
      EDN=(TQ(1)*EPSJ)**2
C- is to test for decreasing the order,
      E=(TQ(2)*EPSJ)**2
C- comparison for errors of current order NQ,
      EUP=(TQ(3)*EPSJ)**2
C- is to test for increasing the order,
      BND=(TQ(4)*EPSJ)**2
C- is used to test for convergence of the
C- correction iterates
C------ '<--Leaving  Node %A:'
      ENDIF
C------ '<--Leaving  Node %GXZ:'
      GOTO09845
09846 CONTINUE
C: RESET PR2 *
C------ '-->Entering Node %GX2:'
      PR2=1.D0/(((D/E)**(.5D0/DBLE(LNQ)))*BIAS2+DFLTZR)
C------ '<--Leaving  Node %GX2:'
      IF(NQ.EQ.1)THEN
      RH=PR2
      CALLRESCAL(Y,NYDIM,RH,RMAX,RC,LNQ)
      IDOUB=LNQ+1
      ELSE
C: RESET PR1 *
C------ '-->Entering Node %GX1:'
      SUM=0.D0
      DO09837I=1,N
      SUM=SUM+(Y(I,LNQ)/YMAX(I))**2
09837 CONTINUE
      PR1=1.D0/(((SUM/EDN)**(.5D0/DBLE(NQ)))*BIAS1+DFLTZR)
C------ '<--Leaving  Node %GX1:'
      IF(PR1.LE.PR2)THEN
      RH=PR2
      IF(METH.EQ.3)THEN
      HNT(1)=H*RH
      CALLCOSET(METH,NQ,EL,TQ,MAXORD,IDOUB)
      RC=RC*EL(1)/OLDL0
      OLDL0=EL(1)
C: the constants E,EUP,EDN,BND *
C------ '-->Entering Node %A:'
      EDN=(TQ(1)*EPSJ)**2
C- is to test for decreasing the order,
      E=(TQ(2)*EPSJ)**2
C- comparison for errors of current order NQ,
      EUP=(TQ(3)*EPSJ)**2
C- is to test for increasing the order,
      BND=(TQ(4)*EPSJ)**2
C- is used to test for convergence of the
C- correction iterates
C------ '<--Leaving  Node %A:'
      CALLRESCAL(Y,NYDIM,RH,RMAX,RC,LNQ)
      ELSE
      CALLRESCAL(Y,NYDIM,RH,RMAX,RC,LNQ)
      IDOUB=LNQ+1
      ENDIF
      ELSE
C: CHANGE NQ,RESET CONST.  *
C------ '-->Entering Node %GXR:'
      LNQ=NQ
      NQ=NQ-1
      RH=PR1
      HNT(1)=H*RH
      CALLCOSET(METH,NQ,EL,TQ,MAXORD,IDOUB)
      RC=RC*EL(1)/OLDL0
      OLDL0=EL(1)
C: EDN,E,EUP,BND *
C------ '-->Entering Node %A:'
      EDN=(TQ(1)*EPSJ)**2
C- is to test for decreasing the order,
      E=(TQ(2)*EPSJ)**2
C- comparison for errors of current order NQ,
      EUP=(TQ(3)*EPSJ)**2
C- is to test for increasing the order,
      BND=(TQ(4)*EPSJ)**2
C- is used to test for convergence of the
C- correction iterates
C------ '<--Leaving  Node %A:'
      CALLRESCAL(Y,NYDIM,RH,RMAX,RC,LNQ)
      IDOUB=LNQ+1
C------ '<--Leaving  Node %GXR:'
      ENDIF
      ENDIF
09845 CONTINUE
09848 CONTINUE
09887 CONTINUE
C------ '<--Leaving  Node %GX:'
      GOTO09893
09894 CONTINUE
C- ITER==MAXITE
C- The corrector iterations failed to converge in MAXITE
C- tries.
      NITER=NITER+1
      IF(OLDJAC)THEN
C- If partials are not up to date , they are
C- reevaluated for the next try.
      EVALJA=.TRUE.
C- Refresh JACOBIAN
      ELSE
C: Otherwise, the Y array is retracted to its values
C                 before prediction, and H is reduced, if possible.
C                 If not, a no convergence exit is taken. *
C------ '-->Entering Node %GY:'
      t=tOLD
      RMAX=RMXFAI
      DO09834J1=1,NQ
      DO09831J2=J1,NQ
      J=NQ-J2+J1
      DO09828I=1,N
      Y(I,J)=Y(I,J)-Y(I,J+1)
09828 CONTINUE
09831 CONTINUE
09834 CONTINUE
      IF(ABS(H).LE.(HMIN*1.00001D0))THEN
      RETCOD=2
      ELSE
C-        write(*,'(a,i7,2i3,1p,2e12.3)') ' GY: maxer, NQ, kflag, H, Hu:',
C-                                             maxer, NQ, kflag, H, Hused;
      RH=RHCORR
      if(ifail.EQ.13.and.Kbad.GT.0)then
C-          @wterm ' GY: ifail 13';
C- NQ=1; LNQ=2; -- if we want 1st order
C-           @wterm' Initial from Y  R0=',Y(Kbad-1,1),'  u0=',Y(Nzon+Kbad-1,1);
C-           @wterm' Initial from Y  R1=',Y(Kbad,1),'  u1=',Y(Nzon+Kbad,1);
C-           @wterm' Initial dr=',Y(Kbad,1)-Y(Kbad-1,1),
C-                '  predicted dr=',
C-                 Y(Kbad,1)-Y(Kbad-1,1)+(Y(Nzon+Kbad,1)-Y(Nzon+Kbad-1,1))*h;
      endif
      IF(METH.EQ.3)THEN
      HNT(1)=H*RH
      CALLCOSET(METH,NQ,EL,TQ,MAXORD,IDOUB)
C: the constants E,EUP,EDN,BND *
C------ '-->Entering Node %A:'
      EDN=(TQ(1)*EPSJ)**2
C- is to test for decreasing the order,
      E=(TQ(2)*EPSJ)**2
C- comparison for errors of current order NQ,
      EUP=(TQ(3)*EPSJ)**2
C- is to test for increasing the order,
      BND=(TQ(4)*EPSJ)**2
C- is used to test for convergence of the
C- correction iterates
C------ '<--Leaving  Node %A:'
      RC=RC*EL(1)/OLDL0
      OLDL0=EL(1)
      CALLRESCAL(Y,NYDIM,RH,RMAX,RC,LNQ)
      ELSE
      CALLRESCAL(Y,NYDIM,RH,RMAX,RC,LNQ)
      IDOUB=LNQ+1
      ENDIF
      ENDIF
C------ '<--Leaving  Node %GY:'
      ENDIF
09893 CONTINUE
      IF(.NOT.(.NOT.EVALJA.or.RETCOD.NE.0))GOTO09978
      GOTO09990
09989 CONTINUE
C- while RETCOD==0
C------ '<--Leaving  Node %G:'
C: ALL RETURNS ARE MADE THROUGH THIS SECTION DEPENDING ON RETCOD.
C        H IS USED IN HOLD TO ALLOW THE CALLER TO CHANGE H ON THE NEXT STEP. *
C------ '-->Entering Node %X:'
      GOTO(09825,09824,09823),RETCOD
C- Control passes to _Esac for RETCOD==4
      GOTO09822
09825 CONTINUE
      KFLAG=-1
      GOTO09822
09824 CONTINUE
      KFLAG=-2
      GOTO09822
09823 CONTINUE
      RMAX=RMXNOR
09822 CONTINUE
      HOLD=H
      JSTART=NQ
C------ '<--Leaving  Node %X:'
      RETURN
      END
