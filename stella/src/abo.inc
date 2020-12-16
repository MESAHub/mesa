_Define
   @LM   LubvB -- freq. bin for typical Feautrier band
   @Feautrier  Feau
   @HOLD_BAL   HOLDBL   --   OLD STEP FOR BALANCE
   @ffloor     1.d-50   -- f (occ.n) floor - machine dependent

-- Parameter(Mfreq=@Mfreq);
     -- think about max Nfreq !
 COMMON/CONUR/EIT,DST,BBRCNR(5);
 COMMON/BAL/EL(MAXDER+1),YENTOT(MAXDER+1),ETOT0,ELVOL,ELSURF,ELTOT,
            TPSURF,@HOLD_BAL,ELOST,EKO,RADBEG;
 _include nstep;
 REAL*8 TPMAX(MAXDER+1),TQ(4);
 COMMON/TAU/TAU(Mzon+1),FLUX(Mzon);
 common/tauubvri/tauU(Mzon),tauB(Mzon),tauV(Mzon),tauR(Mzon),tauI(Mzon);
 COMMON/PHOT/XJPH,DMPH,RPH,TPH,PLPH,VPH,CHEMPH,GRVPH,HP,JPH;
 PARAMETER (NFUNC=6);
 REAL*4 WORK(Mzon+2,NFREQ)/*NFUNC IF NFREQ < NFUNC*/,WRK(Mzon,4);
 REAL*8 WRKX(Mzon),WORKX(Mzon+2);
 COMMON/STEPD/WRKX,WORKX,TPHOT,TEFF,WORK,WRK,NPHOT,NZM;
       /*       1 - LG(T), 2 - LG(PL), 3 - LG(P), 4 - LG(S)      */
 PARAMETER(TMCRIT=1.D-6,TPNSE=5.D0,EPGROW=0.02D0);
 Common/RUTP/Ry(Mzon),Uy(Mzon),Ty(Mzon),Press(Mzon),Rho(Mzon);
 COMMON/TOO/TOO,KO,KNTO,TO(KOMAX),STO(KOMAX),NTO(KOMAX);
 Parameter(Lcurdm=1000); -- Dimension of the Tcurv array
 Real Tcurv;
 Integer NFRUSED; -- an integer array which store exact number of used freqs
 REAL*8 Flsave; -- remove this to read old flx files!
 Common/Curve/tcurv(8,Lcurdm),Depos(Lcurdm),Flsave(MFREQ+1,Lcurdm),
              NFRUSED(Lcurdm),Lsaved;
 LOGICAL BEGRUN;
 Common/BEGR/BEGRUN;
 CHARACTER*80 Model,Sumprf,Sumcur, Depfile, Flxfile;
 COMMON/Files/Model,Sumprf,Sumcur, Depfile, Flxfile;
 CHARACTER*1 app; -- dummy var for constructing of Opafile
 Logical Givdtl;
 Common/ABGrap/NSTA,NSTB,TcurA,TcurB,Givdtl;
         --No. of steps & Time in days
 REAL*8 MBOL,MU,MB,MV,MR,MI,MBOL1;
 COMMON/COLOR/MBOL,MU,MB,MV,MR,MI,UMB,BMV,MBOL1,LubvU,LubvB,LubvV,LubvR,LubvI,
    Lyman;
 COMMON/DETAIL/QRTarr(Mzon),UUarr(Mzon),ArrLum(Mzon),Acc(Mzon);
 Common/XYZ/XA,YA,URM;
