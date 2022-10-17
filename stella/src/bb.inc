-- _define @Mfreq  150
--_define @Mfreq  126
-- parameter (Tbba = 2000.d0, epsbb = 1.d-3);
parameter (Tbba = 1.d04, epsbb = 1.d-3);
--Parameter(Mfreq = @Mfreq);
-- think about max Nfreq !
real*8 alumnu(Mfreq), Blc(Mfreq), Ab(Mfreq), nu(Mfreq), numn(Mfreq), Xpl(Mfreq),
wbb(Mfreq);
common/frq/ alumnu, Blc, Ab, nu, numn, R14, frbba;
