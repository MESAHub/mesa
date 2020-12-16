Overview of Atmosphere Choices
==============================

.. toctree::
   :maxdepth: 2
   :hidden:

   t-tau
   table
   irradiated
   paczynski
   fixed-teff
   fixed-tsurf
   fixed-psurf
   fixed-psurf-tsurf

The ``atm_option`` control determines MESA's overall approach to
evaluating :math:`T_{\rm surf}` and :math:`P_{\rm surf}`, as follows:

.. list-table:: The ``atm_option`` control
   :widths: 25 75

   * - Value
     - Meaning

   * - ``'T_tau'``
     - :doc:`t-tau`

   * - ``'table'``
     - :doc:`table`

   * - ``'irradiated_grey'``
     - :doc:`irradiated`

   * - ``'paczynski_grey'``
     - :doc:`paczynski`

   * - ``'fixed_Teff'``
     - :doc:`fixed-teff`

   * - ``'fixed_Tsurf'``
     - :doc:`fixed-tsurf`

   * - ``'fixed_Psurf'``
     - :doc:`fixed-psurf`

   * - ``'fixed_Psurf_and_Tsurf'``
     - :doc:`fixed-psurf-tsurf`
