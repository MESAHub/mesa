.. _high_rot_darkening:

******************
high_rot_darkening
******************

.. tags:: star, rotation, critical-rotation, gravity-darkening, mass-loss, intermediate-mass, main-sequence

This test evolves a rapidly rotating 3 |Msun|, Z=0.02 model with
gravity darkening enabled.

The case relaxes the model toward high rotation, keeps the model near
a specified fraction of critical rotation early in the run, and
exercises rotationally enhanced mass loss near critical rotation.

How this test passes
====================

The active test-suite entry expects the termination string
``termination code: xa_central_lower_limit`` and checks the final saved
model ``final.mod``.
