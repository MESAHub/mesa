Overview of auto_diff module
============================

.. toctree::
   :maxdepth: 1

The ``auto_diff`` module provides Fortran derived types that support automatic
differentiation via operator overloading. Users will not generally
need to interact with this module, but it can be used within
run_star_extras to make derivatives easier to calculate (e.g. in the
implicit hooks like ``other_surface``).

Usage is by writing ``use auto_diff`` at the top of a module or routine.
This imports types such as ``auto_diff_real_4var_order1``, which supports first-order derivatives
with respect to up to four independent variables.
A variable of this type could be declared via::

    type(auto_diff_real_4var_order1) :: x

This variable then holds five fields: ``x%val`` stores the value of ``x``.
``x%d1val1`` stores the derivative of `x` with respect to the first independent
variable. ``x%d1val2`` is the same for the second independent variable, and so on.
All ``d1val_`` fields are initialized to zero when the variable is first set.
Once an `auto_diff` variable is initialized, all mathematical operations can be performed
as they would be on a ``real(dp)`` variable. `auto_diff` variables also interoperate with
``real(dp)`` and ``integer`` types.
So for instance in the following ``f%d1val1`` stores df/dx and ``f%d1val2`` stores df/dy.::

    x = 3d0
    x%d1val1 = 1d0
    y = 2d0
    y%d1val2 = 1d0
    f = exp(x) * y + x + 4

Similar types are included supporting higher-order and mixed-partial
derivatives.  These derivatives are accessed via e.g. ``d2val1``
(:math:`\partial^2 f/\partial x^2`), ``d1val1_d2val2`` (:math:`\partial^3 f/\partial x \partial y^2`).

An additional special type ``auto_diff_real_star_order1`` provides support
for first-order derivatives accessed using arrays.
This type contains a value (``x%val``) and an array of first partial derivatives
with respect to at least as many variables as the solver in ``MESA/star``.
This type is meant to make it easy to write equations and then, after the fact,
change the basis of independent variables or re-index them.
The current indices are defined in ``star_data/public/star_data_def.inc``.
E.g., if ``my_var`` is of type ``auto_diff_real_star_order1``,
then ``my_var% d1Array(i_lnT_00)`` contains the derivative of ``my_var``
with respect to ``lnT`` at the same mesh point.
