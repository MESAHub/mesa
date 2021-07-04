.. _test_memory:

***********
test_memory
***********

This test case program checks MESA's memory management.
It is designed primarily to be run inside the valgrind leak-checking tool,
and is based on code provided originally by Warrick Ball.

The program ``test_memory`` builds a pre-main sequence model, and then deletes it,
twice in succession - see ``src/test_memory.f90``.
For speed reasons, the ./rn script does NOT use valgrind. 
Rather, valgrind must be launched by hand:

.. code-block:: console

  valgrind --leak-check=full --error-limit=no  ./test_memory

To hide a number of false positives, due to internal gfortran issues, one can also use the supplied supression file ``mesa.supp``:

.. code-block:: console

 valgrind --leak-check=full --error-limit=no --suppressions=mesa.supp ./test_memory

Other usefull valgrind options:

.. code-block:: console

 --show-leak-kinds=all --track-origins=yes  : To add additonal output
 --log-file=valgrind.out : Redirect valgrind's output to the file ``valgrind.out''
 --gen-suppressions=all : Generates supression rules as valgrind runs, helpful when you want to edit mesa.supp


Last-Updated: 03Jul2021 (MESA 094ff71) by fxt.
