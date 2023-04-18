.. highlight:: console
.. _known_bugs:

**********
Known bugs
**********

This page lists a number of known bugs or issues in released versions of MESA. Where possible 
we will also list work arounds, but for some bugs the only option will be to update to
a newer version of MESA. Note this list is NOT comprehensive, users should check this first if they have an 
issue but it may not be complete.


r22.11.1
========

Rates
=====

There is a bug in the rate selection code that certain endothermic weak reactions are not added to the nuclear network. These are 
r_be10_wk-minus_b10, r_ni66_wk-minus_cu66, and r_h3_wk-minus_he3. Other weak reactions with heavier parents may also be affected.

A separate issue also meant we are missing the rate r_he4_ap_li7 as the reverse rate of r_li7_pa_he4.

Both issues will effect previous versions of MESA as well.

Both issues have been fixed in the git main branch.

See `gh-491 <https://github.com/MESAHub/mesa/issues/491>`_ and `gh-497 <https://github.com/MESAHub/mesa/issues/497>`_

RTI
===

A bug has existed since shortly after r15140 where RTI mixing will be effectively zero in a model even with the ``RTI_flag=.true.``

This has now been fixed in the git main.

See `gh-503 <https://github.com/MESAHub/mesa/issues/503>`_



r22.05.1
========

Convective Premixing
--------------------

Convective premixing (CPM) has not worked properly since release r15140. CPM was broken by the
removal of the ``lnPgas_flag``, which caused some of the necessary EOS updates to be missed after
CPM updates the abundances in mixed cells. CPM does not need ``lnPgas_flag``, but it does require
EOS updates at constant pressure. This will be fixed in future releases.

See `gh-425 <https://github.com/MESAHub/mesa/issues/425>`_


Invalid location for overshoot boundary
---------------------------------------

Sometimes MESA will crash with an error similar to this:

.. code-block:: shell

    s%top_conv_bdy(i)= F
    D(k)   0.0000000000000000    
    s%D_mix(k-1)   1.1101956346180402    
    s%overshoot_D_min   100.00000000000000    
    Invalid location for overshoot boundary: cz_bdy_dq, dq= -0.13040604669743103        1.4532774141478022E-003
            0 terminate reason: nonzero_ierr

This bug effects many previous versions of MESA as well. This has been fixed in `gh-400 <https://github.com/MESAHub/mesa/issues/400>`_ .
The solution is to update to a newer MESA version.



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


