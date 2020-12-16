      integer*4 maxionstage_nuvar, nelements_nuvar
      parameter(maxionstage_nuvar = 28)
      parameter(nelements_nuvar = 20)
      integer*4 nuvar_index_element(nelements_nuvar),
     &  nuvar_atomic_number(nelements_nuvar), nuvar_nelements
      real*8
     &  nuvar(maxionstage_nuvar+1,nelements_nuvar),
     &  nuvarf(maxionstage_nuvar,nelements_nuvar),
     &  nuvart(maxionstage_nuvar,nelements_nuvar),
     &  nuvar_dv(maxionstage_nuvar,maxionstage_nuvar,nelements_nuvar)
      common/nuvarblk/nuvar, nuvarf, nuvart, nuvar_dv,
     &  nuvar_index_element, nuvar_atomic_number, nuvar_nelements
