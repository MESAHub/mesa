  _Define @skip  4    -- step in mass zones for tables in composition
          @skip (Mzon/2) -- test and uniform runs
--          @skip (Mzon/50) -- standard
--          @skip (Mzon/45) -- for Mzon=45
          ;
Character*80 Opafile;
-- Common/opabs100/opacabs100(Nfreq); -- bad if here; now in argument for subr. 'opacity"
-- real*4 hptabab,hptabsc,hpsavsc; -- new!
real*4 hptabab,hptababron,hptabsc,hpsavsc; -- new!
Common/Opsave/
        hpsavsc(Nfreq,@Ntab,@Ntab,Mzon/@skip+1,@Ns);
Common/OpBand/
        TpTab(@Ntab),RhoTab(@Ntab),STab(@Ns),
        Wavel(Nfreq),YATab(Natom),
        hptabab(Nfreq,@Ntab,@Ntab,Mzon/@skip+1,2), 
        hptababron(Nfreq,@Ntab,@Ntab,Mzon/@skip+1,2),
        hptabsc(Nfreq,@Ntab,@Ntab,Mzon/@skip+1,2), 
--Peter Hoeflich:       EpsBand(@Nfreq,@Ntab,@Ntab,@Ns),
--        EppBand(@Nfreq,@Ntab,@Ntab,@Ns),
        Msta,Nrho,Ntp,im; -- im number of mixture
-- old Dimension hpbanab1(Nfreq,@Ntab,@Ntab,Mzon/@skip), 
   real*4 hpbanab1(Nfreq,@Ntab,@Ntab,Mzon/@skip+1), 
          hpbanabron1(Nfreq,@Ntab,@Ntab,Mzon/@skip+1),
          hpbansc1(Nfreq,@Ntab,@Ntab,Mzon/@skip+1), 
          hpbanab2(Nfreq,@Ntab,@Ntab,Mzon/@skip+1), 
          hpbanabron2(Nfreq,@Ntab,@Ntab,Mzon/@skip+1),
          hpbansc2(Nfreq,@Ntab,@Ntab,Mzon/@skip+1);
Equivalence (hptabab(1,1,1,1,1),hpbanab1(1,1,1,1)),
            (hptababron(1,1,1,1,1),hpbanabron1(1,1,1,1)),
            (hptabsc(1,1,1,1,1),hpbansc1(1,1,1,1)),
            (hptabab(1,1,1,1,2),hpbanab2(1,1,1,1)),
            (hptababron(1,1,1,1,2),hpbanabron2(1,1,1,1)),   
            (hptabsc(1,1,1,1,2),hpbansc2(1,1,1,1));
Common/tintrp/ 
          stmlog(@ns),tdlog,thaplog1,thaplog2,istold,Opafile;
Common/dumfreq/dumFreq(Nfreq+1),dumFreqmn(Nfreq),dumwavel(Nfreq);
