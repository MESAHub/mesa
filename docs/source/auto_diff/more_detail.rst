===================================
More detail on the auto_diff module
===================================

Adapted from a `post <https://adamjermyn.com/posts/auto_diff_mesa/>`_ on Adam Jermyn's website describing the implementation of ``auto_diff`` in MESA.

What is auto_diff?
==================

Forward-mode automatic differentiation via operator overloading:

-  Forward-mode means we calculate the chain rule as we go.
-  Each variable in the calculation needs to be able to track derivative
   information.
-  Variables need to know how the chain rule applies to each operation.
-  Fortran source files are generated automatically by a python program.

   -  This allows robust support for many different
      functions/operators/derivative configurations.

What does it look like in Fortran?
==================================

.. code:: fortran

      type :: auto_diff_real_1var_order1
         real(dp) :: val
         real(dp) :: d1val1
      end type auto_diff_real_1var_order1

The types themselves are simple! Here’s the ``1var_order1`` type, which
supports 1 independent variable through 1 (first order) derivative.
``val`` stores the value, ``d1val1`` stores the derivative with respect
to the independent variable.

Concretely, we might set up a variable like this:

.. code:: fortran

   type(auto_diff_real_1var_order1) :: x
   x = 1d0 ! Sets val to 1, zero's out d1val1
   x%d1val1 = 1d0 ! Says dx/d(val1) = 1.
                  ! Often used as a shorthand for saying 'x' is the independent variable.

And we might perform operations with these variables:

.. code:: fortran

   f = sin(x) ! Now f%val = sin(1), f%d1val1 = cos(1)
   f = exp(x) ! Now f%val = e, f%d1val1 = e
   f = pow3(x + 1) ! Now f%val = 8, f%d1val1 = 6
   f = f + x ! Now f%val = 9, f%d1val1 = 7

Note that we **do not** support assignments like
``(real) = (auto_diff)``. Why? Because we don’t want to accidentally
lose the derivative information, and a ``real`` type doesn’t have
anywhere to put it!

So if you want to get the value you have to use ``f%val``, and if you
want the derivative info that’s in ``f%d1val1``.

Other types support more derivatives and more variables. The general
pattern is ``Nvar_orderM`` will support **all** derivatives through
:math:`m`-th total order in all combinations of :math:`n` variables. So
for instance ``auto_diff_real_2var_order2`` supports ``d1val1``,
``d1val2``, ``d1val1_d1val2``, ``d2val1``, ``d2val2``, which are the
five combinations of (mixed) partial derivatives of two variables up to
total order 2. Note that the mixed ones are always ordered by the
``val`` index, not the ``d`` index, e.g. \ ``d2val1_d1val2`` is how
you’d write one of the third order mixed partials.

How does it work in Fortran?
============================

Behind the scenes are ludicrously large Fortran files, which begin like:

.. code:: fortran

   module auto_diff_real_1var_order1_module
         use const_def, only: dp, ln10, pi
         use utils_lib
         use support_functions
         use math_lib
      
         implicit none
         private
      public :: auto_diff_real_1var_order1, &
         assignment(=), &
         operator(.eq.), &
         operator(.ne.), &
         operator(.gt.), &
         operator(.lt.), &
         operator(.le.), &
         operator(.ge.), &
         make_unop, &
         make_binop, &
         sign, &
         safe_sqrt, &
         operator(-), &
         exp, &
         expm1, &
         ...

It goes on for a while. This is just exporting all the many (many)
operators that get overloaded. Scrolling down we find the
implementations of these overloaded routines, like

.. code:: fortran

      function expm1_self(x) result(unary)
         type(auto_diff_real_1var_order1), intent(in) :: x
         type(auto_diff_real_1var_order1) :: unary
         unary%val = expm1(x%val)
         unary%d1val1 = x%d1val1*exp(x%val)
      end function expm1_self

and

.. code:: fortran

      function add_self(x, y) result(binary)
         type(auto_diff_real_1var_order1), intent(in) :: x
         type(auto_diff_real_1var_order1), intent(in) :: y
         type(auto_diff_real_1var_order1) :: binary
         binary%val = x%val + y%val
         binary%d1val1 = x%d1val1 + y%d1val1
      end function add_self

The operators are all labelled as either unary or binary. Binary
operators generally are named by the types they work with
(e.g. ``add_self`` adds two ``auto_diff`` types, ``add_self_int`` adds
an ``auto_diff`` type and an ``integer``, etc.).

Sometimes the operators are a little inscrutable:

.. code:: fortran

      function dim_self(x, y) result(binary)
         type(auto_diff_real_1var_order1), intent(in) :: x
         type(auto_diff_real_1var_order1), intent(in) :: y
         type(auto_diff_real_1var_order1) :: binary
         real(dp) :: q0
         q0 = x%val - y%val
         binary%val = -0.5_dp*y%val + 0.5_dp*x%val + 0.5_dp*Abs(q0)
         binary%d1val1 = -0.5_dp*y%d1val1 + 0.5_dp*x%d1val1 + 0.5_dp*(x%d1val1 - y%d1val1)*sgn(q0)
      end function dim_self

The reason for this is that they’re all auto-generated by python
scripts, in a way that optimizes for (Fortran) runtime speed at all
costs.

Derivative Operators
--------------------

Each ``auto_diff`` type additionally has some number of derivative
operators, one per independent variable. These work like:

.. code:: fortran

   df_dx = differentiate_1(f)
   df_dy = differentiate_2(f)

The idea here is that you might want an ``auto_diff`` type which itself
represents the derivatives of another ``auto_diff`` variable (so you can
propagate higher order derivatives through later operations). This is
what let’s Skye do things like writing the pressure as

.. code:: fortran

   p = pow2(den) * differentiate_2(F)

and have ``p`` still contain derivative information.

These methods can’t fill in higher order derivatives than exist. In the
above example ``F`` has a third derivative with respect to ``rho``.
``p`` is a derivative of ``F`` with respect to ``rho``, so we don’t know
enough to construct the third derivative of ``p`` with respect to
``rho``. This is handled by just zeroing out the derivatives we don’t
know.

We considered using NaN’s instead of zeros, following a philosophy that
you should know very clearly when you’ve mistakenly read a missing entry
(silent failure is bad). The problem with using NaN’s here is that we
want to be able to run MESA with floating point error (FPE) checking
turned on as a way to catch numerical problems, and if we assign NaN to
variables routinely that becomes impossible.

Custom Operators
----------------

Two functions I want to highlight are ``make_unop`` and ``make_binop``:

.. code:: fortran

      function make_unary_operator(x, z_val, z_d1x) result(unary)
         type(auto_diff_real_1var_order1), intent(in) :: x
         real(dp), intent(in) :: z_val
         real(dp), intent(in) :: z_d1x
         type(auto_diff_real_1var_order1) :: unary
         unary%val = z_val
         unary%d1val1 = x%d1val1*z_d1x
      end function make_unary_operator
      
      function make_binary_operator(x, y, z_val, z_d1x, z_d1y) result(binary)
         type(auto_diff_real_1var_order1), intent(in) :: x
         type(auto_diff_real_1var_order1), intent(in) :: y
         real(dp), intent(in) :: z_val
         real(dp), intent(in) :: z_d1x
         real(dp), intent(in) :: z_d1y
         type(auto_diff_real_1var_order1) :: binary
         binary%val = z_val
         binary%d1val1 = x%d1val1*z_d1x + y%d1val1*z_d1y
      end function make_binary_operator

Let’s focus on ``make_unop``. It takes as input an ``auto_diff``
variable and ``z_val`` and ``z_d1x``. The latter two are the
specification of a function and its derivative with respect to an
arbitrary independent variable, evaluated at that value of ``x``.
``make_unop`` then propagates through the chain rule to apply that
function to ``x`` and give a result which inherits derivatives from
``x``. These helper routines are there to perform **variable
substitutions**. The idea is you might know the function
``z(independent_variable)`` but want to instead have ``z(x)`` (which has
different derivatives because ``x`` may itself be a complicated function
of independent variables). ``make_binop`` does the same but for binary
operators.

As far as I’m aware these functions only get used in the Skye equation
of state in MESA, where we play some tricks with custom operators, but
they’re there if you ever need to do a variable substitution.

Array Types
-----------

There are two special types that break the mold:

.. code:: fortran

      type :: auto_diff_real_star_order1
         real(dp) :: val
         real(dp) :: d1Array(33)
      end type auto_diff_real_star_order1

This is the basic ``auto_diff`` type used everywhere in ``MESA/star``.
Instead of 33 different independent variable entries it puts them all in
an array. The meaning of these is set in
``MESA/star_data/public/star_data_def.inc``, where you’ll find

.. code:: fortran

         ! auto_diff constants for solver variables
         ! used to access auto_diff_real_star_order1 d1Array
         integer, parameter :: i_lnd_m1 = 1
         integer, parameter :: i_lnd_00 = 2
         integer, parameter :: i_lnd_p1 = 3
         integer, parameter :: i_lnT_m1 = 4
         integer, parameter :: i_lnT_00 = 5
         integer, parameter :: i_lnT_p1 = 6
         integer, parameter :: i_w_m1 = 7
         integer, parameter :: i_w_00 = 8
         integer, parameter :: i_w_p1 = 9
         integer, parameter :: i_lnR_m1 = 10
         integer, parameter :: i_lnR_00 = 11
         integer, parameter :: i_lnR_p1 = 12
   ...

which tells the solver which indices correspond to which variables in
the array. Hence ``d1Array(5)`` corresponds to the derivative with
respect to ``lnT`` in the current cell, ``d1Array(6)`` with respect to
``lnT`` in the next cell, and so on.

If you need to change the number of independent variables, you can do
that by updating (1) the entry in the auto_diff config file (both for
the star and tdc types), (2) adding new indexing parameters to
``star_data_def.inc``, and (3) adding new helper routines to
``MESA/star/private/auto_diff_support.f90`` to handle your new
independent variables.

There are also lots of helper routines in
``MESA/star/private/auto_diff_support.f90`` for manipulating these
objects, including ways to shift the indexing so ``p1 -> 00`` (and
vice-versa), ways to generate e.g. \ ``lnT(k)`` with the appropriate
derivative setup (``d1Array(1:4,6:33)==0``, ``d1Array(5)==1``), etc.

The other special one is

.. code:: fortran

      type :: auto_diff_real_tdc
         real(dp) :: val
         real(dp) :: d1val1
         real(dp) :: d1Array(33)
         real(dp) :: d1val1_d1Array(33)
      end type auto_diff_real_tdc

This type is only used in the time-dependent convection (TDC) code, and
exists because we needed a type that has a derivative with respect to
one additional variable (the superadiabaticity on a face) and needed all
mixed partial derivatives with all of the star solver variables.

How does it work in Python?
===========================

So how does the python side generate these files?

Config Files
------------

In ``MESA/auto_diff/config`` there are a bunch of files, one per
``auto_diff`` type. These are yaml files, and look like:

::

   name: auto_diff_real_2var_order1
   orders: [[1,0],[0,1]]
   array: False

This says:

   Make a type named ``auto_diff_real_2var_order1``. It has to have all
   partial derivatives up to and including the first derivative with
   respect to the first variable and the first derivative with respect
   to the second variable.It does not store derivatives as an array.

Another example:

::

   name: auto_diff_real_2var_order3
   orders: [[3,0],[2,1],[1,2],[0,3]]
   array: False

which says

   Make a type named ``auto_diff_real_2var_order3``. It has to have all
   partial derivatives up to and including the third derivative with
   respect to the first variable, the (2,1) mixed partial, the (1,2)
   mixed partial, and the third derivative with respect to the third
   variable. It does not store derivatives as an array.

Finally, the star example:

::

   name: auto_diff_real_star_order1
   orders: [[1]]
   array: True
   fixed_length: True
   array_length: 33

which says

   Make a type named ``auto_diff_real_star_order1``. It stores
   derivatives as arrays of fixed length 33 and has to have all partial
   derivatives up to and including the first derivative with respect to
   each component of the array.

Parser
------

You can regenerate the ``auto_diff`` Fortran source by calling
``python parser.py`` in the ``python`` directory. The parser is
reasonably straightforward. It begins by getting the list of config
files:

.. code:: python

   # Get config files
   config_path = '../config'
   config_files = [f for f in listdir(config_path) if isfile(join(config_path, f)) and '.config' in f]
   config_files = [join(config_path, f) for f in config_files]

It then makes two lists of files. The ``compilation_list`` are all the
files that ``make`` will need to act on, and the ``use_list`` is all the
modules that need to be shared by the ``auto_diff`` public interface.

.. code:: python

   # compilation_list stores a list of all the fortran files that will need compiling.
   # This is used in the makefile.
   compilation_list = []
   compilation_list.append('support_functions.f90')

   # use_list stores a list of all private auto_diff modules that need importing into the public auto_diff module.
   use_list = []
   use_list.append(tab + 'use support_functions')

We then loop over all config files. For each, we read out the relevant
info:

.. code:: python

   data = load(fi, Loader=Loader)

   # gfortran does not (as of September 2021) support variable-length
   # arrays in parameterized-derived-types. So stick with fixed-length
   # arrays. If this changes in the future you can set fixed_length
   # to False and use variable-length arrays as desired.
   if data['array'] and data['fixed_length']:
     array_length = data['array_length']
   else:
     array_length = None

Note that we can’t do variable length arrays. The Python side can
generate parameterized derived ``auto_diff`` types supporting variable
length arrays, but gfortran doesn’t actually implement the F2003 spec
and so won’t compile it. Some versions of ifort worked with this
functionality but I can’t remember which. The gfortran bug report is
`here <https://gcc.gnu.org/bugzilla/show_bug.cgi?id=97818>`__.

Then construct the list of all partial derivatives required:

.. code:: python

   # Read desired highest-order partial derivatives
   partials = list(Partial(orders, data['array']) for orders in data['orders'])

This fills in all lower-order derivatives needed to fulfill the
requested list of higher-order ones (e.g. if you request a third order
derivative, this adds in a second and a first as well). The ``Partial``
data type is defined in ``partial.py`` and just has some helper methods
for helping implement the chain rule.

That done, we build the types and write them out to files, appending
them to the ``compilation_list`` and ``use_list``:

.. code:: python

   # Build auto_diff type with those and all lower-order derivatives.
   adr = AutoDiffType(data['name'], data['array'], array_length, partials)

   out_fi = open('../private/' + data['name'] + '_module.f90', 'w+')
   out_fi.write(py_to_fort(make_auto_diff_type(adr, unary_operators, binary_operators, comparison_operators, intrinsics)))
   out_fi.close()
   compilation_list.append(data['name'] + '.f90')
   use_list.append(tab + 'use ' + data['name'] + '_module')

AutoDiffType
------------

The ``AutoDiffType`` class lives in ``auto_diff_type.py``. This type is
the internal representation of an ``auto_diff`` Fortran type on the
Python side. It’s initialized as

.. code:: python

   class AutoDiffType:
       def __init__(self, name, array, array_length, partials):
           '''
           Stores a list of partials that is complete, in the sense that there is enough information
           to compute the chain rule within that set of partials, and sorted by total order.
           '''

So you pass the partials you want, the variable name, and some
information about arrays.

Now that I look again, it seems that we don’t need the array
information, because the ``Partial`` type already has that. So the
``array`` and ``array_length`` entries could be safely removed here.

The initialization has a few important pieces. First, we work out the
complete set of partial derivatives we need:

.. code:: python

   # Complete the partials list
   partials = set(partials)
   complete = False
   while not complete:
     ps = list(partials)
     for p in ps:
       partials.update(p.completion_partials())
     if len(partials) == len(ps):
       complete = True

The routine ``completion_partials`` returns any additional partial
derivatives that a given partial needs to be able to propagate in the
chain rule. For instance :math:`\partial_x\partial_y^2` in the chain
rule needs access to :math:`\partial_x\partial_y` and
:math:`\partial_y^2`, so its completion will return those two. We just
keep calling ``completion_partials`` till it stops returning new
derivatives.

Next, we put these in a sorted order so we can refer to them
consistently:

.. code:: python

   self.partials = sorted(list(partials), key=lambda p: [p.net_order, tuple(-o for o in p.orders)])

Finally, we construct the sets of partials of unary operators and binary
operators out to the maximum order represented. These, too, are needed
by the chain rule:

.. code:: python

   self.unary_partials = sorted(list(Partial((i,), False) for i in range(self.max_order+1)), key=lambda p: [p.net_order, tuple(-o for o in p.orders)])
   self.binary_partials = sorted(list(Partial((i,j), False) for i in range(self.max_order+1) for j in range(self.max_order+1) if i+j <= self.max_order), key=lambda p: [p.net_order, tuple(-o for o in p.orders)])

You can think of these as :math:`\partial_x f(x,y)` and
:math:`\partial_y f(x,y)`, which you need to compute the chain rule for
:math:`\partial_u (f(x(u),y(u)))`.

The rest of the class specification is full of functions that construct
the various operators that appear on the Fortran side. For instance

.. code:: python

       def specific_unary_operator_function(self, operator_name, operator):
           '''
           Returns a function which implements the specified unary operator.
           '''

           function_name = operator_name + '_self'
           function_arguments = [('x', self.declare_name(ref='*'), 'in')]
           function_result = ('unary', self.declare_name(ref='x'))
           function_body, function_declarations = unary_specific_chain_rule(self, operator, fixed_length=self.array_length)
           function_body = function_declarations + function_body

           # Special case handling for safe_log
           if 'safe' in operator_name:
               for i in range(len(function_body)):
                   function_body[i] = function_body[i].replace('log', 'safe_log')

           return FortranFunction(function_name, function_arguments, function_result, function_body)

takes as input an operator’s name and the operator itself (as a
``sympy`` function) and returns a valid Fortran function (as a string)
implementing the derivative propagation logic. Most of this is wrapper
logic: all the magic and complicated stuff that goes in the body gets
constructed in
``unary_specific_chain_rule(self, operator, fixed_length=self.array_length)``
(and there are equivalent functions for binary operators).

chain_rule
----------

The real magic on the Python side all happens in ``chain_rule.py``.
That’s where functions like
``unary_specific_chain_rule(self, operator, fixed_length=self.array_length``
are defined.

There are four functions in this file. They are each labelled ``unary``
or ``binary``, after the kind of operator they represent, and
``specific`` or ``generic``. The ``generic`` ones are used to write the
`Custom Operator <#custom-operators>`__ routines and the ``specific``
ones are used to implement actual specific operators like ``exp`` and
``+`` and so on.

How do they work?

specific
~~~~~~~~

This is a bit complicated.

Everything here uses ``sympy`` for calculus and algebra, which means
most of what we’re doing is setting up lots of ``sympy`` variables and
manipulating them.

We start by making symbols for the independent variables:

.. code:: python

   # Construct sympy variables corresponding to the various independent variables.
   # These never appear on the Fortran side, but we keep the naming consistent to correspond to the
   # names in partial_orders.
   # So these are called val1, val2, ..., valN.
   indep_symbol_str = ' '.join(auto_diff_type.partials[0].val_name(i) for i in range(auto_diff_type.num_independent))
   indep_syms = wrap_element(symbols(indep_symbol_str, real=True))

Then we make symbols for the places we’ll store the derivatives (this is
where ``d1val1`` comes from!):

.. code:: python

   # Construct sympy variables corresponding to the various derivatives dx/d(...).
   # Note that these variable names correspond to the names we'll use on the Fortran side, so
   # we can just directly map sympy expressions to strings and get valid Fortran :-)
   # Hence these are called x%d1val1, x%d2val1, ..., x%d1val2, x%d2val2, ..., x%d1val1_d1val2, ...
   # The first integer in each 'd_val_' block is the number of derivatives,
   # the second is the independent variable those derivatives are with respect to.
   x_symbol_str = ' '.join(auto_diff_type.partial_str_in_instance('x', p).replace(':','colon') for p in partials)
   x_syms = wrap_element(symbols(x_symbol_str, real=True))

We then represent ``x`` as a power series in terms of its partial
derivative symbols:

.. code:: python

   # Construct x as a power series in terms of its partial derivatives (sym) with respect to the independent
   # variables (indep).
   x = 0
   for p,sym in zip(*(partials, x_syms)):
     term = sym
     for order, indep in zip(*(p.orders, indep_syms)):
       term = term * indep ** order / factorial(order)
     x = x + term

And then call our operator on ``x`` to get ``z(x)``:

.. code:: python

   z = operator(x)

The reason we play around with the power series is so that ``z`` has an
explicit representation in terms of the partial derivatives of ``x``,
which in turn are explicitly represented as individual ``sympy``
symbols.

With all that done, we actually extract derivatives. This starts with a
few lists:

.. code:: python

   expressions = []
   left_hand_names = []
   derivatives = []

Here ``expressions`` is the list of derivative expressions we’ll build,
``left_hand_names`` is the corresponding list of e.g. \ ``d1val1``
(which appear on the left-hand side in the Fortran code), and
``derivatives`` appears to be an unused list that I forgot to delete.

We then iterate over all required partials:

.. code:: python

       for p in partials:

For each, we construct the left-hand side of the expression:

.. code:: python

   unary_symbol_str = auto_diff_type.partial_str_in_instance('unary', p).replace(':','colon')

The ``replace`` business here is just to make sure we only use valid
``sympy`` symbols. There’s a lot of that all over this code (string
replacements to avoid invalid or reserved ``sympy`` symbols, followed by
back-replacement at the very end right before we write the file).

If life were simple, we’d then just ask ``sympy`` for the derivative at
the right order. But some use cases require ``auto_diff`` to support
non-differentiable functions like ``abs`` and ``>`` and ``min`` and so
on. Those spawn Dirac Delta’s when you try differentiate them. Which is
awful because (1) we can’t do anything with those in any numerical
methods and (2) they’re zero everywhere but a set of measure zero, so we
don’t care about them. So we get a bunch of logic that special cases
Dirac Delta and a few related objects and zero’s them out:

.. code:: python

   d = z
   for order, indep in zip(*(p.orders, indep_syms)):
     d = diff(d, indep, order)
     d = d.replace(DiracDelta, zero_function) # Diract Delta is only non-zero in a set of measure zero, so we set it to zero.
     d = d.replace(sign, sgn) # simplify can do weird things with sign, turning it into the piecewise operator. We want to avoid that so we redirect to a custom sgn function.
     d = d.replace(Derivative, zero_function) # Eliminates derivatives of the Dirac Delta and Sign, which are non-zero only at sets of measure zero.
     d = d.subs(indep, 0)

This is taking the derivatives one at a time, clearing out garbage as it
arises.

Life would be nice if this were all we had to do, but we want the
resulting Fortran code to be fast, so we do some simplifications:

.. code:: python

   d = simplify(d, measure=weighted_count_ops, force=True, ratio=1)

More on that later.

The rest of the routine is just a bunch of string manipulation to get
everything into the right format for Fortran. You can print the results
as they accumulate if you’re interested to see how the substitutions
gradually turn a pile of algebra into valid Fortran.

generic
~~~~~~~

The generic ones cheat by calling the specific ones on a dummy function.
They start by constructing symbols corresponding to the partial
derivatives of some general function ``z(x)`` out through the highest
order we care about:

.. code:: python

   # Construct the symbols corresponding to the partial derivatives of z.
   # These are d1z, d2z, ..., dNz, giving dz/dx, d^2 / dx^2, and so on.
   z_symbol_strs = ['z_' + str(p).replace('val1','x') for p in auto_diff_type.unary_partials]
   z_symbol_str = ' '.join(z_symbol_strs)
   z_syms = wrap_element(symbols(z_symbol_str, real=True))

We then construct a Taylor series out of these symbols:

.. code:: python

   def operator(x):
     # Construct z as a power series in terms of its partial derivatives (z_syms) with
     # respect to the x.
     z = sum(sym * x**p.orders[0] / factorial(p.orders[0]) for sym,p in zip(*(z_syms, auto_diff_type.unary_partials)))
     return z

Then we call ``unary_specific_chain_rule`` to give us the chain rule
code for this dummy operator, and that gets everything in terms of the
partial derivatives of ``z(x)``, which we can then supply as inputs to
the custom operator builders.

make_auto_diff_type
-------------------

This file puts it all together, going over all the functions and all the
Fortran boiler plate and doing a bunch of accounting to make sure every
``function blah`` gets closed by ``end function blah`` and so on. It’s
super boring.

measure
-------

This is where all performance optimizations happen. We use the built-in
``sympy`` function ``simplify``, but with a twist. We don’t care how
complicated the functions are, we care how fast they are. And moreover
speed is actually set by how many divides and special function calls we
have. So we have to tell ``simplify`` about all of that. In
``measure.py`` specify crudely how much each function call costs:

.. code:: python

   # 'basic' here means roughly a one-cycle op.
   # 'div' is division, which takes ~30 cycles.
   # 'special' is a special function, which takes ~1000 cycles.
   # DIRACDELTA and DERIVATIVE get eliminated in post-processing and so are free.
   special = 1000
   div = 30
   basic = 1
   weights = {
       'SIN': special,
       'COS': special,
       'TAN': special,
       'TANH': special,
       'COSH': special,
       'SINH': special,
       'ASIN': special,
       'ACOS': special,
       'ATAN': special,
       'ATANH': special,
       'ACOSH': special,
       'ASINH': special,
       'EXP': special,
       'LOG': special,
       'POW': special,
       'ADD': basic,
       'MUL': basic,
       'NEG': basic,
       'SUB': basic,
       'HEAVISIDE': basic,
       'ABS': basic,
       'DIV': div,
       'SGN': basic,
       'POWM1': div,
       'SSQRT': special,
       'DIRACDELTA': 0,
       'DERIVATIVE': 0
   }

Then we have a function that goes through a ``sympy`` expression
counting function calls and tallying them up:

.. code:: python

   def weighted_count_ops(expr_original, verbose=False):

Understanding this code requires a decent amount of knowledge of
``sympy``\ ’s API, but suffice it to say that we’re crawling an abstract
syntax tree and counting instances of functions as we encounter them.

functions
---------

The ``functions.py`` file defines all the supported ``auto_diff``
functions in ``sympy`` language. Not much more to say there.

Helper Methods
--------------

There are a bunch of boring helper methods in ``routine.py`` (for
spitting out valid Fortran routines), and ``utils.py`` (for string
manipulation and a few performance optimizations like
``pow(x,N) -> powN(x)``).
