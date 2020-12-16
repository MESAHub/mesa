Atmosphere Structure Building
=============================

MESA's ``atm`` module is primarily used to evaluate the :math:`T_{\rm
surf}` and :math:`P_{\rm surf}` values for the boundary
conditions. However, when writing output files for pulsation codes
(using the ``save_pulse_data...``  controls), it is sometimes desirable
to append the (spatial) atmosphere structure data to the interior
structure data; that way, the extension of modes into the stellar
atmosphere can naturally be modeled by the pulsation code.

To append atmosphere structure data in this way, set the
``add_atmosphere_to_pulse_data`` control to ``.true.``. A couple of
caveats apply:

- This capability is only available when ``atm_option`` has the values
  ``'T_tau'`` or ``'Paczynski_grey'``. For other values (e.g.,
  ``'table'``), there are insufficient data to determine what the
  spatial structure of the atmosphere should be.

- It's important to realize that calculation of atmosphere structure
  data has NO effect on the evolution of the star, OR on MESA's
  standard profile, history and photo outputs. Only the pulsation
  output files are affected.

When building atmosphere structure data, a number of controls influence
the calculation, as follows:

.. list-table:: Controls for building atmosphere structure data
   :widths: 25 75

   * - Control
     - Meaning

   * - ``atm_build_errtol``
     - Error tolerance for integrating the :math:`r(\tau)` [and
       possibly also :math:`T(\tau)`] equations

   * - ``atm_build_dlogtau``
     - Maximum logarithmic spacing :math:`\Delta\log_{10}(\tau)` for atmosphere structure data

   * - ``atm_build_tau_outer``
     - Outermost optical depth :math:`\tau_{\rm outer}` for atmosphere
       structure data. The structure will extend from :math:`\tau_{\rm
       surf}` to this value, and contain approximately
       :math:`[\log_{10}(\tau_{\rm surf}) - \log_{10}(\tau_{\rm
       outer})]/\Delta\log_{10}(\tau)` points
