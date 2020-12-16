Fixed Surface Temperature
=========================

When the ``atm_option`` control has the value ``'fixed_Tsurf'``, MESA
sets :math:`T_{\rm surf}` to the value specified by the
``atm_fixed_Tsurf`` control. The effective temperature is determined
from the Eddington-grey :math:`T(\tau)` relation

.. math::

   T_{\rm surf}^{4} = \frac{3}{4} T_{\rm eff}^{4} \left( \tau_{\rm surf} + 2/3 \right).

See :ref:`here <tau-surf>` for a discussion of how the surface
optical depth :math:`\tau_{\rm surf}` is evaluated.

Then, MESA sets :math:`P_{\rm surf}` to equal the radiation pressure
associated with :math:`T_{\rm surf}`.
