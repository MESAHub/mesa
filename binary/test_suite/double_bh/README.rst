.. _double_bh:

*********
double_bh
*********

.. tags:: binary, binary-evolution, massive-star, black-hole, black-hole-binary, chemically-homogeneous-evolution, rotation, mass-transfer, roche-lobe-overflow, tidal-synchronization, contact-binary, point-mass, helium-depletion

This test evolves a close binary with initially massive stars toward a
double black-hole configuration through chemically homogeneous
evolution.

The case evolves both stars, enables tidal synchronization and mass
transfer, and turns a component into a point mass after helium
depletion.  The test terminates after the other component also reaches
helium depletion.

How this test passes
====================

The active test-suite entry expects the output string ``Terminate due
to helium depletion`` and checks the final saved model ``final.mod``.
