.. _wind_fed_bhhmxb:

****************
wind_fed_bhhmxb
****************

.. tags:: binary, binary-evolution, high-mass-x-ray-binary, massive-star, black-hole, black-hole-accretion, wind-accretion, stellar-winds, mass-loss, eddington, mass-transfer, roche-lobe-overflow

This test evolves a wind-fed high-mass X-ray binary with a massive
donor and a black-hole companion.

The case checks wind accretion and Roche-lobe overflow while limiting
black-hole accretion at the Eddington rate.  The binary extras verify
both sub-Eddington and super-Eddington phases and check the associated
mass-loss bookkeeping.

How this test passes
====================

The active test-suite entry expects the output string ``Properly
tested sub and super Eddington phases`` and checks the final saved
model ``final.mod``.
