Coding style
============

In free form Fortran files, indentation is 3 spaces.

Indentation/whitespace rules are codified in an `EditorConfig`_ file located in the root of the repository (``.editorconfig``).

.. _EditorConfig: https://editorconfig.org/


Exponents
---------

When raising a number to a floating-point power, use the ``pow()`` routine provided by ``math_lib``:

.. code-block:: fortran

    pow(x,1.5d0)

instead of

.. code-block:: fortran

    x**1.5d0

Use ``powN(x)`` to raise ``x`` to an integer power ``N``.
For ``real`` types you should *always* use these for integer powers greater than 2.
``**2`` has so far not been a problem and is tolerated though not recommended.
For ``auto_diff`` types you *must* use ``powN`` for integer powers *including* 2.
Doing otherwise will result in a compiler error.


Errors
------

In general a subroutine should return as its last argument a integer ``ierr`` which denotes whether the 
subroutine succeeded or not. A value of ``ierr=0`` is to be used to show the subroutine succeeded, 
any other value of ierr denotes an error. Any meaning in different non-zero values of ``ierr`` 
can only be determined by inspecting the subroutine that set the ``ierr``. 

The typical way to handle things is:

.. code-block:: fortran

    integer, intent(inout) :: ierr
    ierr = 0

    call my_subroutine(x,ierr)
    if (ierr /= 0) return

for non-critical errors.


Stop
----

When signalling a critical error which can not be recovered from you should use:

.. code-block:: fortran

    use utils_lib, only: mesa_error

    call mesa_error(__FILE__,__LINE__)

Which will generate a backtrace where the error occurs, which aids in debugging. The 
subroutine ``mesa_error`` can also have a optional error message:


.. code-block:: fortran

    use utils_lib, only: mesa_error

    call mesa_error(__FILE__,__LINE__,"An error occurred")


Do not use:

.. code-block:: fortran

    stop 1

As it can become difficult to distinguish where in the code called ``stop`` if multiple functions use ``stop 1``


Doubles
-------

The preferred manner to declare a double precision variable is:

.. code-block:: fortran

    real(dp) :: x

instead of 

.. code-block:: fortran

    double precision :: x

When using a numerical value in an expression you must make sure it is evaluated as a double.
Thus use:

.. code-block:: fortran

    y1 = 1.1d0 * x
    ! or
    y2 = 1.1_dp * x

Do not leave values as a bare float:

.. code-block:: fortran

    y3 = 1.1 * x

As the ``1.1`` gets interpreted as a single precision value, and will lead ``y3`` to have a different value 
to ``y1`` or ``y2``.


OMP critical blocks
-------------------

OMP critical blocks allow the programmer to specify that a section of code should only be executed by one thread at a time.
They can also be given a name:

.. code-block:: fortran

    !$omp critical my_block

and this name should differ from any other code entities (e.g. subroutines).

Each named critical block will be executed by one thread at a time. Different named critical blocks can be executed
at the same time. However, all unnamed critical blocks act like one block and thus can not be executed in parallel.
Therefore you should always named your OMP critical blocks to ensure the best performance.  

Do not name your OMP critical block with a name that has already been used for a variable, procedure, module or any other object.


Formatting
----------

Use explicit formats for any ``write`` statements.  Different compilers use different default formats, which can lead to spurious
failures when strings are compared. e.g. when printing some floating point number ``x``, instead of ::

  write(*,*) x

use ::

  write(*, '(1pd26.16)') x

Unformatted statements are likely to cause unit tests to fail.  They also make it difficult to compare output from runs with
different compilers.

Some helpful formats are provided in ``include/formats``.


Constants
---------

The ``const`` module defines many commonly used mathematical
(e.g. ``pi``) and physical constants (e.g. ``hbar``), which should be
used for consistency across the code.  This includes simple fractions
(e.g. ``one_third``) and simple functions of mathematical constants
(e.g. ``sqrt2``, ``pi4 = 4*pi``).


Environment variables
---------------------

If making a new environment variable then the variable should be prefixed with ``MESA_`` to ensure we don’t collide with other variables.
