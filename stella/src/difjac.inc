_Define
@hapHlim  1.5d0 -- 4*.375
@FJBEG  FSAVE(NVARS*NZON-Ncnd+(NZON-Ncnd)
@FJ2    @FJBEG*(L-1)+Km+1)       -- FJ(Km+1,L)
@FJ1    @FJBEG*(L-1)+Km)         -- FJ(Km,L)
@FJ0    @FJBEG*(L-1)+Km-1)       -- FJ(Km-1,L)
@DFJ1   @FJBEG*(L-1)+Km+NYDIM)   -- DFJ(Km,L)
@FHBEG  FSAVE(NVARS*NZON-Ncnd+KRAD+(NZON-Ncnd)
@FH2    @FHBEG*(L-1)+Km+1)       -- FH(Km+1,L)
@FH1    @FHBEG*(L-1)+Km)         -- FH(Km,L)
@FH0    @FHBEG*(L-1)+Km-1)       -- FH(Km-1,L)
@DFH1   @FHBEG*(L-1)+Km+NYDIM)   -- DFH(Km,L)
@FJ1C   @FJBEG*L+Km)             -- FJ(L+1)
@FJ1A   @FJBEG*(L-2)+Km)         -- FJ(L-1)
@FH1C   @FHBEG*L+Km+KRAD)        -- FH(L+1)
@FH1A   @FHBEG*(L-2)+Km+KRAD)    -- FH(L-1)
@Tp_CUT  1.D+1/UTP
@TAULIM 0.d0  --- CFR   -- == (4/3)/(UR*URHO==(Cappa unit)**-1)
--@EddH1  EddH(Km)    -- for EDTM==T & EDJRAD(Km)    -- for EDTM==T
@AJAC   'JAC=JAC+1;IRS(JAC)=KROW;AJAC(JAC)'
@LEFT       6.D0 ;    --  COEFF. FOR LEFT BOUNDARY AT R=Rcore
--  @LEFT==6 AT core WHERE FH==.5*INTEGRAL(F*MU*DMU)
  ;
   PARAMETER (ALFCON=1.D0);
   Common/Hydro/dm0,dm1,dm2,Rodm,RSodm,R0,R1,R2,RSm1,RS0,RS1,RS2,RC0,RC1,RC2,
                Rcinv0,RCinv1,RCinv2,U0,U1,U2,Q1,Q2,QU1,QU2,UDIVR,
                UU1,UU2,UQ1,UQ2,AQ1,AQ2,QUHALF;
   Common/Termo/PL0,PL1,PL2,ET1,ET2,
                Tp0,Tp1,Tp2,PQ1,PQ2,PQpl1,PQpl2,PT1,PT2,TPTQ1,TPTQ2,
                TpT1,TpT2,Tp0up4,Tp1up4,Tp2up4,
                PTPL1,PTPL2,PTT1,PTT2,ETPL1,ETPL2,ETT1,ETT2,
                DMETIN,CDMETI,CK2ET1,
                HEATRA,FLUM,FIMP,DTPOLD,
                ENG0,ENG1,ENG2,
                ENGPL1,ENGPL2,
                ENGT1,ENGT2;
   Common/Tratst/TOTJ,TOTW,TOTWPR;
   Common/Convo/PP0,PP1,PP2,CP0,CP1,CP2,GRA0,GRA1,GRA2,VCSQ0,VCSQ1,
                FC0,FC1,ELMIX0,ELMIX1,CFC0,CFC1,
                CPPL0,CPPL1,CPPL2,
                CPT0,CPT1,CPT2,
                GRAPL0,GRAPL1,GRAPL2,
                GRAT0,GRAT1,GRAT2,
                PRPL0,PRPL1,PRPL2,
                PRT0,PRT1,PRT2;
   Common/Fluxo/tau1(Nfreq),tau2(Nfreq), -- for hapH
                FL0,FL1,Flcor1,Flcor2,Fllf0,Flrt0,Fllf1,Flrt1,
                FLcore, CAPcor,
                CAP0,CAP1,CAP2,
                CAP1D,CAP1T,
                CAPT0,CAPT1,CAPT2,
                CAPPL0,CAPPL1,CAPPL2,
                CM0INV,CM1INV;
   Common/Transo/DJNU(NFREQ),B24(NFREQ),
                 DFH(NFREQ),DFJ(NFREQ),DUR(NFREQ),
                 DFJDPL(NFREQ),DFHDPL(NFREQ),
                 HAPL(NFREQ),
                 HAPL1(NFREQ),HAPL2(NFREQ),
                 HAPL1T(NFREQ),HAPL2T(NFREQ),
                 HAPL1D(NFREQ),HAPL2D(NFREQ),
                 HAPH(NFREQ),HAPH1(NFREQ),
                 HAAB(NFREQ),HAAB1(NFREQ),HAAB2(NFREQ),
                 HAAB1T(NFREQ),HAAB2T(NFREQ),
                 HAAB1D(NFREQ),HAAB2D(NFREQ),
                 FHcore(NFREQ),BLCKD(NFREQ);
   Common/Burno/YDOT1,YDOT2,CCPL1,CCPL2,CCTP1,CCTP2,AS1,AS2,YCARB1,YCARB2;
