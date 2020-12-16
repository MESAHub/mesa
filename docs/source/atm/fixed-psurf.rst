Fixed Surface Pressure
======================

When the ``atm_option`` control has the value ``'fixed_Psurf'``, MESA
sets :math:`P_{\rm surf}` to the value specified by the
``atm_fixed_Psurf`` control.

Then, MESA sets :math:`T_{\rm surf}` from the Eddington-grey
:math:`T(\tau)` relation

.. math::

   T_{\rm surf}^{4} = \frac{3}{4} T_{\rm eff}^{4} \left( \tau_{\rm surf} + 2/3 \right).

The effective temperature :math:`T_{\rm eff}` is determined from the
stellar radius and luminosity; see :ref:`here <tau-surf>` for a
discussion of how the surface optical depth :math:`\tau_{\rm surf}` is
evaluated.
