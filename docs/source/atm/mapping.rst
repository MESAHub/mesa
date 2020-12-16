Mapping from Old to New
=======================

For those wanting to convert inlists from the old (pre-11869) approach
to the new one described here, the mapping is as follows:

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Old setting
     - New settings

   * - ``which_atm_option == 'simple_photosphere'``
     - | ``atm_option == 'T_tau'``
       | ``atm_T_tau_relation == 'Eddington'``
       | ``atm_T_tau_opacity = 'fixed'``

   * - ``which_atm_option == 'grey_and_kap'``
     - | ``atm_option == 'T_tau'``
       | ``atm_T_tau_relation == 'Eddington'``
       | ``atm_T_tau_opacity = 'iterated'``

   * - ``which_atm_option == 'Eddington_grey'``
     - | ``atm_option == 'T_tau'``
       | ``atm_T_tau_relation == 'Eddington'``
       | ``atm_T_tau_opacity = 'varying'``

   * - ``which_atm_option == 'XXXXXX_tables'``
     - | ``atm_option == 'table'``
       | ``atm_table = 'XXXXXX'``

   * - ``which_atm_option == 'grey_irradiated'``
     - | ``atm_option == 'irradiated_grey'``
       | ``atm_irradiated_opacity == 'fixed | varying'`` 

.. note::
   If MESA is run using one of the old options, it will print out the
   appropriate mapping from the table above, and then halt execution.
