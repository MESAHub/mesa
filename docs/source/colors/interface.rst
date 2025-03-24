=======================
colors module Interface
=======================

XXX TODO XXX
The primary entry point to the colors module is through the routine ``colors_get``.  Broadly, one provides the density, temperature and full composition information, as well as information about the free electron fraction and electron chemical potential (required for calculation of the Compton scattering opacity). Then the module returns the opacity and its temperature and density derivatives.


.. literalinclude:: ../../../colors/public/colors_lib.f90
   :language: fortran
   :start-at: subroutine colors_get
   :end-at: end subroutine colors_get
   :linenos:
