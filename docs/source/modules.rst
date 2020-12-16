********************
Module documentation
********************

This page lists the MESA modules alphabetically by name and briefly
summarizes their purpose.

Each MESA module has its own directory with the same general
structure, including a standard set of subdirectories and scripts. The
standard subdirectories for each module are ``make``, ``private``,
``public``, and ``test``. The test directory has ``make`` and ``src``
directories for the program that tests the module when it is
created. The make directory has the makefile for the library and will
hold the object files and ``.mod`` files that are created by the
compiler. The public directory has the sources for the interface to
the library, while the private directory has sources for the parts of
the implementation that are meant for internal use only. For example,
if you want to see what is available in the eos module, look in
``eos/public/eos_lib.f90`` for the routines and
``eos/public/eos_def.f90`` for the data.


.. _atm:

Atmospheres (``atm``)
=====================

MESA uses the atmosphere (``atm``) module to obtain a surface
temperature (:math:`T_{\rm surf}`) and surface pressure (:math:`P_{\rm
surf}`), representing the conditions at the base of the stellar
atmosphere. These values are applied as boundary conditions when
evolving the interior model.

.. note::

   MESA treats the atmosphere separately, via these boundary
   conditions, because the physics governing the atmosphere is often
   quite different than in the interior.

Historically, MESA decided how :math:`T_{\rm surf}` and :math:`P_{\rm
surf}` are calculated using the ``which_atm_option`` control. To
improve consistency between the atmosphere and interior calculations,
and to offer more flexibility to users, this control was replaced in
revision 11869 by the ``atm_option`` control, plus a number of other
subsidiary controls.

.. toctree::

   atm/overview
   atm/mapping
   atm/structure


.. _chem:

Element data (``chem``)
=======================

The ``chem`` module provides data on the properties of elements and
isotopes (e.g., atomic masses).  It also defines solar abundance
patterns as reported in various references.


.. _const:

Constants (``const``)
=====================

The ``const`` module defines a range of mathematical constants (e.g.,
pi), physical constants (e.g., hbar), astronomical constants (e.g.,
Msun), and other fixed values (e.g., version number).


.. _eos:

Equation of state (``eos``)
===========================

The ``eos`` module provides the equation of state.

.. toctree::
   :maxdepth: 1

   eos/overview
   eos/defaults


.. _kap:

Opacities (``kap``)
===================

The ``kap`` module provides radiative opacities combined with
conductive opacities.

.. toctree::
   :maxdepth: 1

   kap/overview
   kap/defaults

.. _net:

Nuclear reaction networks (``net``)
===================================

The ``net`` module implements nuclear reaction networks.  


.. _neu:

Thermal neutrinos (``neu``)
===========================

The ``neu`` module provides the specific rates of energy loss via
various thermal neutrino processes.  (Nuclear neutrinos are handed in
``rates`` and ``net``.)


.. _rates:

Nuclear reaction rates (``rates``)
==================================

The ``rates`` module collects thermonuclear reaction rates and weak
reaction rates from a range of sources.


