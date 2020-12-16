  _DEFINE
          @TOLF   2.0D0*EDDH(KM) -- TOLERANCE FOR FH FJ -- try EddJ!!!
          @TP_FLOOR       2.D+3/UTP  -- FOR FJ
  _INCLUDE zone;  -- Mzon
  _INCLUDE fundrad;  -- Fundamental constants
  _INCLUDE commonEve;  
   LOGICAL FRST;
 Parameter(Mfreq=@Mfreq);
     -- think about max Nfreq !
   Common/Kmzon/km,kmhap,Jac,FRST;
--   Common/NiAdap/tday,t_eve,XNifor(Mzon),AMeveNi,KNadap; -- must go to commonEve.inc
    -- since it enters also bgcon*trf and is forgotten there
   COMMON/STCOM1/t,H,HMIN,HMAX,EPS,N,METH,KFLAG,JSTART;
   COMMON/YMAX/YMAX(NYDIM);
   COMMON/YSTIF/Y(NYDIM,MAXDER+1);
   COMMON/HNUSED/HUSED,NQUSED,NFUN,NJAC,NITER,NFAIL;
   COMMON/HNT/HNT(7); -- FOR COSETBGH
   PARAMETER (DELTA=1.d-05);  -- DISTURBANCE CONSTANT standard was 1.d-04
-- cannot be less ~1.d-6 since opacity is saved in single precision
--   PARAMETER (LICN=4*NZ,LIRN=8*NZ/5); -- FOR F01BRF, M28Y12
   PARAMETER (LICN=4*NZ,LIRN=2*NZ); -- FOR F01BRF, M28Y12
             -- maximum LIRN=2*NZ
   Logical NEEDBR; -- F01BRF or M28Y12 needed
   COMMON/STJAC/THRMAT,HL,AJAC(LICN),IRN(LIRN),ICN(LICN),WJAC(NYDIM),
                FSAVE(NYDIM*2),IKEEP(5*NYDIM),IW(8*NYDIM),IDISP(11),
                NZMOD,NEEDBR;
   LOGICAL CONV,   -- .TRUE. IF CONVECTION IS INCLUDED
           CHNCND, -- .TRUE. IF NCND IS CHANGED
           SCAT,   -- .TRUE. IF Scattering included (Hapabs^=Happa)
           SEP;    -- .TRUE. IF absorption and scattering are in separate files
                   --           (i.e. absorption without expansion effect)
   COMMON/CUTOFF/FLOOR(NVARS+1),Wacc(NVARS+1),FitTau,TauTol,
                  Rvis,CONV,CHNCND,SCAT,SEP;
   Logical LTHICK;
   COMMON/THICK/ LTHICK(Nfreq,Mzon);
   COMMON/CONVEC/UC(Mzon),YAINV(Mzon);
   COMMON/RAD/EDDJ(Mzon,Nfreq),
              EDDH(Mzon),
              HEDD(Nfreq),HEDRAD,
              CLIGHT,CKRAD,UFREQ,CFLUX,CCL,CLUM,CLUMF,CIMP,
              NTHICK(NFREQ),NTHNEW(NFREQ),
              NCND,KRAD,NFRUS;-- KRAD=(Mzon-NCND)*NFRUS IN STELLA
   LOGICAL EDTM; -- Eddington factors time-dependent==.true.
   COMMON/RADOLD/HEDOLD,HINEDD,EDTM;
   Common/newedd/EddN(Mzon,Nfreq),HEdN(Nfreq),tfeau; -- for EDTM==T
   Common/oldedd/EddO(Mzon,Nfreq),HEdo(Nfreq),trlx; -- for EDTM==T
   Common/cnlast/Cnlast;
   Common/Dhap/DHaphR(Mzon,Nfreq);

   COMMON/BAND/FREQ(NFREQ+1), -- frequency bins boundaries,
                              -- frequency in units of h/(k*Tpunit)
               FREQMN(NFREQ), -- frequency mean (middle of bin)
               WEIGHT(@MFREQ), -- bandwidth*freq(mean)**3
               HAPPAL(NFREQ), -- dimens-less kappa for Ron's absorp.+scatt.
               HAPABSRON(NFREQ), -- dimens-less kappa for Ron's absorption 
               HAPABS(NFREQ), -- dimens-less kappa for S-B absorption
               DLOGNU(NFREQ); -- (log step in frequency)**-1
-- PARAMETER(NFRMIN=Nfreq/5); -- MINIMUM NFRUS
 PARAMETER(NFRMIN=Nfreq/2); -- MINIMUM NFRUS
 Integer dLfrMax;
 Common/observer/wH(Mfreq),cH(Mfreq),zerfr,Hcom(Mfreq),Hobs(Mfreq),
                               freqob(Mfreq),dLfrMax;
 Parameter(NP=@Npe+@Npc-1);   -- Total number of tangent rays
 Common/famu/fstatic(0:NP+1,Nfreq),
    fobs_corr(0:NP+1,Mfreq),fcom(0:NP+1,Mfreq),amustatic(0:NP+1);
 Common/rays/Pray(0:Np+1),fout(0:NP+1,Mfreq),abermu(0:NP+1),NmuNzon;
    -- fout probably not needed
   COMMON/LIM/ Uplim,Haplim;  -- to BEGRAD
   COMMON/AMM/DMIN,DM(Mzon),DMOUT,AMINI,AM(Mzon),AMOUT;
  -- exactly as in VELTEM
   COMMON/Centr/RCE,Nzon; -- central radius & current Number of zones
   Common/InEn/AMHT,EBurst,tBurst,tbeght; -- Mass of Heated Core, Energy & time
   COMMON/RADPUM/AMNI,XMNi,XNi,KmNick;-- MASS OF Ni CORE, Ni ABUNDANCE
   COMMON/RADGAM/FJgam(Mzon,2),toldg,tnewg; -- zero-moment for gamma-radiation
      -- 1 for old time (toldg)      2 for present time (tnewg)
   COMMON/RADGlg/FJglog(Mzon,2); -- log FJgam
   COMMON/CHEM/CHEM0(Mzon),RTphi(0:Mzon+1),EpsUq;
   COMMON/REGIME/NREG(Mzon);
   double precision NRT;
   COMMON/AQ/AQ,BQ,DRT,NRT; -- BQ pseudo-viscosity for R-T
             -- DRT is the mass distance or optical thickness
             -- and used for weight in artificial viscosity
             -- NRT is the power in pseudo-viscosity (may be noninteger)
   COMMON/AZNUC/ACARB,ZCARB,ASI,ZSI,ANI,ZNI,QCSI,QSINI;
   COMMON/QNRGYE/QNUC,RGASA,YELECT;
--   COMMON/CKN1/CK1,CK2,CFR,CRAP;
   COMMON/CKN1/CK1,CK2,CFR,CRAP,CRAOLD;
   LOGICAL EVALJA, -- EVALUATE JACOBIAN
           OLDJAC, -- JACOBIAN MAY NEED UPDATING
           BADSTE; -- BAD STEP: RETURNED BY DIFMAT - TRY SMALLER
   COMMON/EVAL/EVALJA,BADSTE,OLDJAC;
   _Include stateq;
   COMMON/BURNCC/CC,CCTP,CCPL,YDOT;
   COMMON/ABarr/ YABUN(Natom,Mzon);
 -- UNITS:
   COMMON/UNSTL/UL,UPRESS,UE,UEPS,UCAP,UTIME,UV,UFLUX,UP;
   COMMON/TAIL/KTAIL;
   COMMON/UNINV/UPI,UEI;  -- 1./UP, 1./UE
   COMMON/UNBSTL/UR,UM,UEPRI,ULGP,ULGE,ULGV,ULGTM,
                 ULGEST,ULGFL,ULGCAP,ULGEPS;
