.. highlight:: console
.. _known_bugs:

**********
Known bugs
**********

This page lists a number of known bugs or issues in released versions of MESA. Where possible 
we will also list work arounds, but for some bugs the only option will be to update to
a newer version of MESA. Note this list is NOT comprehensive, users should check this first if they have an 
issue but it may not be complete.


r21.12.1
========

Al26 isomers
------------

After running a model with the al26 isomers in your net, when you run the model again, it may
immediately crashes and prints a backtrace containing:

.. code-block:: shell

    create initial model
    create rate data for r_al26-1_to_al26-2
    create rate data for r_al26-2_to_al26-1

    Program received signal SIGSEGV: Segmentation fault - invalid memory reference.

    Backtrace for this error:
    #0  0x7fd69c02931f in ???
        at /usr/src/debug/glibc-2.33-20.fc34.x86_64/signal/../sysdeps/unix/sysv/linux/x86_64/sigaction.c:670
    #1  0x99a523 in __interp_1d_misc_MOD_do_interp_values

The solution for now is to remove all files in ``$MESA_DIR/data/rates_data/cache/`` before
each MESA run, you may also find that changing the number of OMP threads also fixes the problem.

See `gh-360 <https://github.com/MESAHub/mesa/issues/360>`_

  
Atmosphere in pulse data
------------------------

The control ``add_atmosphere_to_pulse_data`` does not work properly with an Eddington atmosphere (the default), and also crashes if ``atm_T_tau_opacity = 'varying'`` is set. 

See `gh-375 <https://github.com/MESAHub/mesa/issues/375>`_


Colors: bad filter name
-----------------------

If you get an error:

.. code-block:: shell

    bad filter name: 

First check that the name matches in your history_coloumns.list file and your color file. Next check for non-printing characters history_coloumns.list in the filter name. This can bee checked with:

.. code-block:: shell

    cat -A history_columns.list | grep "abs_mag"

Finally, there is a bug if you name any column with ``/`` in it (for instance ``[Fe/H]``). The solution is to rename the column to remove the forward slash.

See `gh-379 <https://github.com/MESAHub/mesa/issues/379>`_

RSP
---

An experimental RSP solver feature was turned on by default, leading to convergence issues in nonlinear model integration. Users should include RSP_do_check_omega = .true. in the &controls section of their inlists to get rid of this issue.



r15140
======


r12778
======


r12115
======


