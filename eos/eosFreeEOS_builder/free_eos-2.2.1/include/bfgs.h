      integer*4 nmax_bfgs
      parameter(nmax_bfgs=300)
      logical fletcher_estimate
      real*8 alpha_previous, deltaf_previous
      common/block_bfgs/alpha_previous, deltaf_previous,
     &  fletcher_estimate
