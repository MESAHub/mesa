This test-suite program is designed as a simple test of MESA's memory
management, and is based on code provided originally by Warrick Ball.

The program builds a pre-main sequence model, and then deletes it,
twice in succession. It is designed primarily to be run inside the
valgrind leak-checking tool.

For speed reasons, the ./rn script does NOT use valgrind. Rather,
valgrind must be launched by hand, for instance via

valgrind --leak-check=full --error-limit=no  ./test_memory

To hide a number of false positives, due to internal gfortran issues, we can also use a supression file:

valgrind --leak-check=full --error-limit=no --suppressions=mesa.supp ./test_memory

Other usefull valgrind options:
--show-leak-kinds=all --track-origins=yes  : To add additonal output
--log-file=valgrind.out : Redirect valgrind's output to the file ``valgrind.out''
--gen-suppressions=all : Generates supression rules as valgrind runs, helpful when you want to edit mesa.supp
