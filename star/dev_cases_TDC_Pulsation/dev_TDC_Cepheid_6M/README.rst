.. dev_TDC_Cepheid_6M:

******************
dev_TDC_Cepheid_6M
******************

This test case evolves a 6 M star into the core Helium burning phase, and stops the evolution when the star reach log Teff = 3.75 during the blue loop.
The model then has the core removed, is remeshed and delivered a kick in the fundamental radial mode using GYRE in MESA. The model should display a
growing kinetic energy given by a positive kinetic energy growth rate per cycle, until finite ampltiude pulsations are achieved, similar to RSP. The
metallicity is chosen such that this test case can be directly compared to the RSP_6M_Cepheid test case.

Initialization of the model
===========================
The initial mass and metallicity of the star is set in ``inlist_extra``

.. literalinclude:: ../../../star/dev_cases_TDC_Pulsation/dev_TDC_Cepheid_6M/inlist_extra

Last-Updated: 2025-10-02 (mesa r25+) by Ebraheem Farag

