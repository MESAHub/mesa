T(tau) Atmospheres
==================

When the ``atm_option`` control is set to ``'T_tau'``, MESA defines
the temperature structure of the atmosphere using a :math:`T(\tau)`
relation. The ``atm_T_tau_relation`` control determines the choice of
relation, as follows:

.. list-table:: The ``atm_T_tau_relation`` control
   :widths: 25 75
	    
   * - Value
     - Meaning

   * - ``'Eddington'``
     - Use the Eddington grey relation :math:`T^{4}(\tau) = 3/4 T_{\rm eff}^{4} (\tau + 2/3)`.

   * - ``'solar_Hopf'``
     - Use an approximate Hopf function calibrated against the Sun
       (see section A.5 of Paxton et al. 2013, MESA II). This is
       equivalent to the fit given by Sonoi et al. (2019, A&A, 621,
       84) to the Vernazza et al. (1981) VAL-C model.

   * - ``'Krishna_Swamy'``
     - Use the approximate Hopf function given by Krishna-Swamy
       (1966, ApJ, 145, 174).
    
For the selected :math:`T(\tau)` relation, the pressure structure of
the atmosphere is obtained by integrating the equation of hydrostatic
equilibrium

.. math::

   \frac{{\rm d}P}{{\rm d}\tau} = \frac{g}{\kappa}

under the assumption that the gravity :math:`g` is spatially
constant. The ``atm_T_tau_opacity`` control determines how the opacity
:math:`\kappa` is evaluated, as follows:

.. list-table:: The ``atm_T_tau_opacity`` control
   :widths: 25 75

   * - Value
     - Meaning

   * - ``'varying'``
     - At each optical depth, evaluate the opacity consistent with the
       local thermodynamic state (typically, density and temperature);
       this is done via calls to MESA's equation-of-state and opacity
       modules. The equation of hydrostatic equilibrium is integrated
       numerically. Controls that affect this integration include
       ``atm_T_tau_errtol`` (sets the error tolerance) and
       ``atm_T_tau_max_steps`` (sets the maximum number of steps).

   * - ``'fixed'``
     - Assume the same opacity at each optical depth (i.e. uniform
       opacity), equal to the the opacity in the outermost cell of
       the interior model. The equation of hydrostatic
       equilibrium is integrated analytically.

   * - ``'iterated'``
     - Assume the same opacity at each optical depth (again,
       uniform opacity), but use iteration to find an opacity
       consistent with the thermodynamic state at the base of the
       atmosphere. The equation of hydrostatic
       equilibrium is integrated analytically. Controls that affect the
       opacity iteration include ``atm_T_tau_errtol`` (sets the error
       tolerance) and ``atm_T_tau_max_iters`` (sets the maximum number of
       iterations).

Once the temperature and pressure structure of the atmosphere is
known, MESA evaluates :math:`T(\tau)` and :math:`P(\tau)` at
:math:`\tau=\tau_{\rm surf}`, providing the necessary :math:`T_{\rm
surf}` and :math:`P_{\rm surf}` values for the boundary
conditions.

.. _tau-surf:

The nominal surface optical depth :math:`\tau_{\rm surf}` is
calculated via

.. math::

   \tau_{\rm surf} = \tau_{\rm base} \times {\tt tau\_factor}

Here, :math:`\tau_{\rm base}` is the optical depth where
:math:`T(\tau) = T_{\rm eff}` (the stellar effective temperature) for
the chosen :math:`T(\tau)` relation, while ``tau_factor`` is a parameter
determined from various controls in the ``&star_job`` namelist (e.g.,
``set_tau_factor``).
