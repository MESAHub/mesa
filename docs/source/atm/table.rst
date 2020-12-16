Table Atmospheres
=================

When the ``atm_option`` control is set to ``'table'``, MESA uses
pre-computed tables from separate model atmosphere calculations to
interpolate :math:`T_{\rm surf}` and :math:`P_{\rm surf}` values for
the boundary conditions. The interpolations are peformed using
:math:`T_{\rm eff}` and :math:`\log g` as independent variables. The
``atm_table`` control determines the choice of table, as follows:

.. list-table:: The ``atm_table`` control
   :widths: 25 75

   * - Value
     - Meaning

   * - ``'photosphere'``
     - Use model atmosphere tables computed for a range of
       metallicites, via a combination of PHOENIX (Hauschildt et
       al. 1999a,b) and Castelli & Kurucz (2003) treatments. Surface
       optical depth is defined so that :math:`T(\tau_{\rm surf}) =
       T_{\rm eff}`.

   * - ``'tau_100'``
     - Use model atmosphere tables computed for solar metallicites,
       via a combination of COND (Allard et al. 2001) and Castelli &
       Kurucz (2003) treatments. Surface optical depth is
       :math:`\tau_{\rm surf} = 100`.

   * - ``'tau_10'``
     - Same as ``'tau_100'``, except that surface optical depth is
       :math:`\tau_{\rm surf} = 10`.
  
   * - ``'tau_1'``
     - Same as ``'tau_100'``, except that surface optical depth is
       :math:`\tau_{\rm surf} = 1`.
  
   * - ``'tau_1m1'``
     - Same as ``'tau_100'``, except that surface optical depth is
       :math:`\tau_{\rm surf} = 0.1`.
  
   * - ``'WD_tau_25'``
     - Use model atmospheres for cool white dwarfs, via Rohrmann et
       al. (2011, MNRAS, 411, 781) treatment. Surface optical depth
       is :math:`\tau_{\rm surf} = 25.12`.

   * - ``'DB_tau_25'``
     - Use model atmospheres for hot(ter) white dwarrfs, provided by Odette
       Toloza based on Detlev Koester’s atmosphere code. Surface
       optical depth is :math:`\tau_{\rm surf} = 25.12`.

When any of these ``'table'`` choices is used, it's possible that the
requested :math:`T_{\rm eff}` and :math:`\log g` lie outside the table
range. When this happens, MESA can optionally switch to an alternative
approach for evaluating :math:`T_{\rm surf}` and :math:`P_{\rm
surf}`. The ``atm_off_table_option`` control determines the
alternative approach adopted, as follows:

.. list-table:: The ``atm_off_table_option`` control
   :widths: 25 75

   * - Value
     - Meaning

   * - ``'T_tau'``
     - Use the same :math:`T(\tau)` approach as described in the
       :doc:`t-tau` section. The ``atm_T_tau_relation`` and
       ``atm_T_tau_opacity`` controls (and others) are used in this case
       for the off-table evaluation.

If the ``atm_off_table_option`` control is left blank, an attempt to interpolate
off the table will cause MESA to halt execution.
