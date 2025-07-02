.. _ppisn:

*****
dev_TDC_Mira
*****

! M-type Mira variable stars
! Ireland, Scholz, and Wood, MNRAS 418, 114-128 (2011)
! C50 5050 L  series
! M = 1.35 Msun         when reaches AGB.  start with 1.5 at ZAMS
! at phase = 0
! L = 6722 Lsun      logL = 3.8
! T = 3271           logT = 3.5
! period = 427 days
! mlt_alpha = 2.0
! age ~ 3e9 years
! see Fig 4 for history plots of logL and logTeff and R

This test case evolves a very massive helium star from the He-ZAMS
up to the ocurrence of a pulsational pair-instability event (see |Marchant2019|).

.. |Marchant2019| replace:: `Marchant et al. 2019 <https://ui.adsabs.harvard.edu/abs/2019ApJ...882...36M/abstract>`__

Initialization of the model
===========================
The initial mass of the helium star is set in ``inlist_extra``

.. literalinclude:: ../../../star/test_suite/ppisn/inlist_extra

In this case we use a :math:`72 M_\odot`

Last-Updated: 2019-11-12 (mesa r12413) by Pablo Marchant

