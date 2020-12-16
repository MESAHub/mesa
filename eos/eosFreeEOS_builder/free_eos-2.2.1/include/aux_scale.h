      integer*4 maxnextrasum_scale, nxextrasum_scale,
     &  iextraoff_scale, naux_scale
      parameter (maxnextrasum_scale = 9)
      parameter (nxextrasum_scale = 4)
      parameter (iextraoff_scale = 7)
      parameter (naux_scale = iextraoff_scale +
     &  maxnextrasum_scale + nxextrasum_scale + 1)
      real*8 aux_underflow
      parameter(aux_underflow = 1.d-250)
      real*8 aux_scale(naux_scale),
     &  sum0_scale, sum2_scale,
     &  extrasum_scale(maxnextrasum_scale),
     &  xextrasum_scale(nxextrasum_scale)
      equivalence
     &  (aux_scale(iextraoff_scale-1), sum0_scale),
     &  (aux_scale(iextraoff_scale  ), sum2_scale),
     &  (aux_scale(iextraoff_scale+1), extrasum_scale),
     &  (aux_scale(iextraoff_scale+maxnextrasum_scale+1),
     &  xextrasum_scale)
      common/aux_scale_block/aux_scale
