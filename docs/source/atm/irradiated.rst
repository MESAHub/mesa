Irradiated Atmosphere
=====================

When the ``atm_option`` control has the value ``'irradiated_grey'``,
MESA defines the temperature structure of the atmosphere using the
:math:`T(\tau)` relation for the irradiated grey-atmosphere model of
Guillot & Havel (2011, A&A 527, A20). A number of controls influence
the evaluation of this model, as follows:

.. list-table:: Controls for irradiated grey-atmosphere
   :widths: 25 75

   * - Control
     - Meaning

   * - ``atm_irradiated_T_eq``
     - Equilibrium irradiation temperature :math:`T_{\rm eq}`

   * - ``atm_irradiated_kap_v``
     - Mean visible opacity :math:`\kappa_{\rm v}`

   * - ``atm_irradiated_kap_v_div_kap_th``
     - Ratio of mean visible opacity :math:`\kappa_{\rm v}` to thermal
       opacity :math:`\kappa_{\rm th}`

The pressure structure of the atmosphere is obtained by integrating
the equation of hydrostatic equilibrium

.. math::

   \frac{{\rm d}P}{{\rm d}\tau} = \frac{g}{\kappa_{\rm th}}

under the assumption that the gravity :math:`g` is spatially
constant. The mean visible opacity, appearing in the :math:`T(\tau)`
relation, is evaluated via

.. math::

   \kappa_{\rm v} =
   \begin{cases}
   \kappa_{\rm th} \times {\tt atm\_irradiated\_kap\_v\_div\_v\_th} &  {\tt atm\_irradiated\_kap\_v\_div\_v\_th} > 0 \\
   {\tt atm\_irradiated\_kap\_v} & \text{otherwise}
   \end{cases}

The ``atm_irradiated_opacity`` control determines how the thermal opacity :math:`\kappa_{\rm th}` is
evaluated, as follows:

.. list-table:: The ``atm_irradiated_opacity`` control
   :widths: 25 75

   * - Value
     - Meaning

   * - ``'fixed'``
     - Assume the same opacity at each optical depth (i.e. uniform
       opacity), equal to the the opacity in the outermost cell of the
       interior model. The equation of hydrostatic equilibrium is
       integrated analytically.

   * - ``'iterated'``
     - Assume the same opacity at each optical depth (again, uniform
       opacity), but use iteration to find an opacity consistent with
       the thermodynamic state at the base of the atmosphere. The
       equation of hydrostatic equilibrium is integrated
       analytically. Controls that affect the opacity iteration
       include ``atm_irradiated_errtol`` (sets the error tolerance)
       and ``atm_irradiated_max_iters`` (sets the maximum number of
       iterations).

.. warning::

  MESA can sometimes run into convergence problems when the
  ``'iterated'`` option is used.

Once the temperature and pressure structure of the atmosphere is
known, MESA evaluates :math:`T(\tau)` and :math:`P(\tau)` at
:math:`\tau=\tau_{\rm surf}`, providing the necessary :math:`T_{\rm
surf}` and :math:`P_{\rm surf}` values for the boundary
conditions. The nominal surface optical depth :math:`\tau_{\rm surf}`
is determined so that

.. math::

    P(\tau_{\rm surf}) = {\tt atm\_irradiated\_P\_surf}
