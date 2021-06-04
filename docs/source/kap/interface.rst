====================
kap module Interface
====================

The primary entry point to the kap module is through the routine ``kap_get``.  Broadly, one provides the density, temperature and full composition information, as well as information about the free electron fraction and electron chemical potential (required for calculation of the Compton scattering opacity). Then the module returns the opacity and its temperature and density derivatives.


.. literalinclude:: ../../../kap/public/kap_lib.f90
   :language: fortran
   :start-at: subroutine kap_get
   :end-at: end subroutine kap_get
   :linenos:
