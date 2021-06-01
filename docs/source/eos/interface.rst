====================
eos module Interface
====================

The primary entry point to the eos module is through the routine ``eosDT_get``.  Broadly, one provides the density, temperature and full composition information.  Then the EOS returns its main set of results, their temperature and density derivatives, and the composition derivatives of a small subset of the results.


.. literalinclude:: ../../../eos/public/eos_lib.f90
   :language: fortran
   :start-at: subroutine eosDT_get
   :end-at: end subroutine eosDT_get
   :linenos:


The underlying EOS is in a density-temperature basis, but if one has only density or temperature, there are search interfaces (``eosDT_get_Rho`` and ``eosDT_get_T``) for searching for the other, given some another known EOS quantity (e.g., ``lnE``).

.. literalinclude:: ../../../eos/public/eos_lib.f90
   :language: fortran
   :start-at: subroutine eosDT_get_T
   :end-at: end subroutine eosDT_get_T
   :linenos:


For legacy reasons, there is also an ``eosPT_get`` entry point.  (The same result could ultimately be achieved via ``eosDT_get_Rho``.)  Internally, this is also implemented as a root-find using the standard density-temperature EOS.

.. literalinclude:: ../../../eos/public/eos_lib.f90
   :language: fortran
   :start-at: subroutine eosPT_get
   :end-at: end subroutine eosPT_get
   :linenos:
