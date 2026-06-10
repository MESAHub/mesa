.. _check_redo:

**********
check_redo
**********

.. tags:: star, redo, forced-redo, restart, photo-restart, checksums, prebuilt-model, hydrogen-depletion

This developer test checks that a run with a forced redo gives the
same answer when restarted from an earlier photo.

The test loads ``zams.mod`` and uses ``run_star_extras`` to request a
redo after model 10 in ``./rn`` but not during ``./re``.  Restarting
from an earlier photo should reproduce the same final answer.

How this test passes
====================

The active test-suite entry expects the usual hydrogen-depletion
termination string, checks ``final.mod``, and restarts from an
automatically selected photo.
