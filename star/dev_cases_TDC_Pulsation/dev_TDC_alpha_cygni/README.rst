.. _ppisn:

*****
alpha_cygni
*****


Properties:
(See F. Schiller and N.Przybilla 2008 and Guzik et al. 2024)
[Fe/H] ~ solar - ish
MZAMS ~ 23 Msun
teff ~ 8525 =-75
log_g ~ 1.1 +- 0.05

This test case evolves a very massive main sequence star from the pms to the alpha-cygni instability strip
up to the ocurrence of a pulsational pair-instability event (see |Marchant2019|).

.. |Marchant2019| replace:: `Marchant et al. 2019 <https://ui.adsabs.harvard.edu/abs/2019ApJ...882...36M/abstract>`__

Initialization of the model
===========================
The initial mass of the helium star is set in ``inlist_extra``

.. literalinclude:: ../../../star/test_suite/ppisn/inlist_extra

In this case we use a :math:`72 M_\odot`

Last-Updated: 2019-11-12 (mesa r12413) by Pablo Marchant

