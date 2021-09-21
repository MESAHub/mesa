===================
A Tour of MESA/star
===================

This document gives a mid-level overview of the structure and control
flow of MESA/star.

The star program
----------------

The ``star`` executable is generated in the context of a work
directory, the template of which lives at ``star/work``.  This
executable is linked against the various MESA modules, including star
itself.

The source file ``src/run.f90`` contains the ``program`` statement.
This file invokes routines that read the ``star_job`` namelist and
then start a run.

The invoked routines live in ``star/job`` in the files
``run_star.f90`` and ``run_star_support.f90``.  (This directory also
contains the code associated with the standard ``run_star_extras``
hooks and the special ``test_suite_extras`` routines used to produce
output for the TestHub.)

The routines in ``run_star_support.f90`` contain key parts of the
control flow.  The routine ``run1_star`` knows how to initialize the
various things needed for a star run and contains the ``evolve_loop``.
This is a simple do loop, each iteration of which is a stellar
evolution step.

.. literalinclude:: ../../../star/job/run_star_support.f90
   :language: fortran
   :start-at: evolve_loop:
   :end-at: end do evolve_loop
   :linenos:                  


Within that loop are calls the routine ``do_evolve_one_step``, which
contains the ``step_loop``.  This loop consists of attempting a
stellar evolution step, and then either accepting the step or
rejecting it and repeating the attempt (via a redo or retry).

.. literalinclude:: ../../../star/job/run_star_support.f90
   :language: fortran
   :start-at: step_loop:
   :end-at: end do step_loop
   :linenos:                  


The stellar evolution step itself is triggered via a call to
``star_evolve_step``, which is a function defined in the
public ``star_lib`` interface.

In ``star/public``, in the file ``star_lib.f90``, this function is

.. literalinclude:: ../../../star/public/star_lib.f90
   :language: fortran
   :start-at: integer function star_evolve_step
   :end-at: end function star_evolve_step
   :emphasize-lines: 9, 18-20
   :linenos:                  
            

You can see that this splits the step into two parts (part1 and
part2).  Those are implemented in ``star/private``.

Before we go deeper, we should become familiar with the ``star_info``
type, which is conventionally assigned the variable name ``s``.

The star_info structure
-----------------------

The ``star_info`` type is defined in the ``star_data`` module.  (It is
split off in order to allow for star to depend on other modules that
need to know the definition of the ``star_info`` type, without
introducing circular dependencies.)

The definition of this derived type occurs in
``public/star_data_def.f90``, but most of the details are deferred to
various include files (``*.inc``).

First, is ``star_data.inc``, which has information about the structure
of the star itself.

.. literalinclude:: ../../../star_data/public/star_data.inc
   :language: fortran

This defers to another set of files that divide things up based on
whether they are inputs for a step or things that are calculated from
these inputs during the process of taking a step.

Some of the most important variables defined in
``star_data_step_input.inc`` are the arrays ``xh(:,:)`` and
``xa(:,:)``, which hold the structure variables and abundance
variables for the model.  These are the "truth" about the stellar
structure.  Other quantities are derived from these through simple
transformations (e.g., ``Rho``) or by the interaction with other
modules (e.g., ``Cp`` via the ``eos``).

The various namelist options for ``star_job``, ``controls``, and
``pgstar`` are also part of the ``star_info`` structure.  They are
also added through ``*.inc`` files, but these live in
``star_data/private``.

.. note:: If you are adding a new inlist option, your best bet is to
   pick an existing option that it is similar to and grep in ``star``
   and ``star_data``.  That will generally show you all the places you
   will need to add code.

There is far too much contained in ``star_info`` to go through it
piece by piece, but you should spend some time perusing theses files
to get a feel for what sorts of things go where.
   
star/private
------------

The guts of MESA live in ``star/private``.

evolve.f90
^^^^^^^^^^

There are many files, but the starting point is ``evolve.f90``.  This
file contains the definitions of the ``part1`` and ``part2`` functions
we saw previously.

Part 1 is mainly preparing for the step (see the functions
``do_evolve_step_part1`` and ``do_step_part1``).  It ends with making
an initial estimate for the mass transfer rate.  (The necessity of the
part1/part2 split is primarily driven by binary, where the two stars
interact in a way that means their mass change rates are related.)

Part 2 contains the bulk of the work associate with a step (see the
functions ``do_evolve_step_part2`` and ``do_step_part2``).

It contains a loop (``implicit_mdot_loop``) that allows this part to
be repeated multiple times in cases where the mass change rate may
depend on the calculated stellar structure.

Within this loop is what one generally thinks of as taking a stellar
evolution step.  First, some of the operator-split actions are taken
(i.e., CPM, diffusion).  Then, the fully-coupled parts are solved in
the call to ``do_struct_burn_mix``.


do_struct_burn_mix
^^^^^^^^^^^^^^^^^^

The function ``do_struct_burn_mix`` lives in ``struct_burn_mix.f90``.

.. warning::
   Forthcoming
