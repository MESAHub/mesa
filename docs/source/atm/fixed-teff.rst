Fixed Effective Temperature
===========================

When the ``atm_option`` control has the value ``'fixed_Teff'``, MESA
sets :math:`T_{\rm surf}` from the Eddington-grey :math:`T(\tau)`
relation

.. math::

   T_{\rm surf}^{4} = \frac{3}{4} T_{\rm eff}^{4} \left( \tau_{\rm surf} + 2/3 \right).

The effective temperature :math:`T_{\rm eff}` is specified by the
``atm_fixed_Teff`` control; see :ref:`here <tau-surf>` for a
discussion of how the surface optical depth :math:`\tau_{\rm surf}` is
evaluated.

MESA then sets :math:`P_{\rm surf}` to equal the radiation pressure
associated with :math:`T_{\rm surf}`.
