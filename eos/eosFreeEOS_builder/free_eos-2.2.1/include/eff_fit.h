      integer*4 maxnf, maxng, nf, ng
      parameter (maxnf = 200)
      parameter (maxng = 200)
      real*8 fbasis(maxnf), gbasis(maxng)
      integer*4 mf, mg, mflg, mglg, mflg2, mglg2, mflf, mglf,
     &  maxmf, maxmg
      parameter (maxmf = 10)
      parameter (maxmg = 10)
      common /eff_fit_block/fbasis, gbasis, nf, ng,
     &  mf, mg, mflg, mglg, mflg2, mglg2, mflf, mglf
