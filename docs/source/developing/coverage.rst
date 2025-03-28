Coverage
========

MESA uses `gcov <https://gcc.gnu.org/onlinedocs/gcc/Gcov-Intro.html>`_
to measure code coverage on the units test run by the ``.\install`` script.
``gcov`` comes as a standard utility with the GNU Compiler Collection (gcc).
The measurement is automatically done with a GitHub action,
which generates html content with `lcov <https://github.com/linux-test-project/lcov>`_
and posts the `coverage report <https://mesastar.org/mesa/>`_ to GitHub Pages.
