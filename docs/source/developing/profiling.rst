Profiling
=========

The MESA SDK contains the `valgrind <https://www.valgrind.org/>`_
instrumentation framework.  A full introduction to valgrind is beyond
the scope of this documentation, but here are some sample applications.


Memory Usage
------------

To profile MESA's memory usage, run inside valgrind's massif tool.

.. note::  Massif is a heap profiler. It performs detailed heap profiling by taking regular snapshots of a program's heap. It produces a graph showing heap usage over time, including information about which parts of the program are responsible for the most memory allocations. The graph is supplemented by a text or HTML file that includes more information for determining where the most memory is being allocated. Massif runs programs about 20x slower than normal.  [vgtools]_

valgrind must be launched by hand, for instance via

.. code-block:: console

    valgrind --tool=massif ./star

Other useful options:

.. code-block:: bash
                
    --stacks=yes # Makes it run slower but will include information about the stack usage

Once finished the tool will write a file massif.out.$PID where $PID is the PID of the 
./star run. If you ctrl-c the run you will still get the output up to that point, but
otherwise you will need to wait till the run finishes before being able to view the output.

To visualize the data use:

.. code-block:: console

    massif-visualizer massif.out.$PID

which is avilable from `<https://github.com/KDE/massif-visualizer>`_ or through your local package manager.


CPU Usage
---------

To profile MESA's CPU usage, run inside valgrind's callgrind tool.

.. note:: Cachegrind is a cache profiler. It performs detailed simulation of the I1, D1 and L2 caches in your CPU and so can accurately pinpoint the sources of cache misses in your code. It identifies the number of cache misses, memory references and instructions executed for each line of source code, with per-function, per-module and whole-program summaries. It is useful with programs written in any language. Cachegrind runs programs about 20--100x slower than normal.  Callgrind is an extension to Cachegrind.  It provides all the information that Cachegrind does, plus extra information about callgraphs. [vgtools]_

valgrind must be launched by hand, for instance via

.. code-block:: console

    valgrind --tool=callgrind ./star


Other useful options:

.. code-block:: bash
                
    --dump-instr=yes --collect-jumps=yes #  Dump extra information about the run into the output
    --dump-every-bb=1000000000 # How often to output a file with the profile information, this will output maybe once every few steps. This lets you start profiling while the run is still going.

Once finished the tool will write a file callgrind.out.$PID.$NUM where $PID is the PID of the 
./star run. $NUM is a seqeuential output number for the each output file.

To visualize the data use:

.. code-block:: console

    kcachegrind callgrind.out.*

which is avilable from `<https://kcachegrind.github.io/html/Home.html>`_ or through your local package manager.


An alternative non-interactive visualization can be produced via `gprof2dot.py <https://github.com/jrfonseca/gprof2dot>`_ and `graphviz <https://graphviz.org/>`_.  

.. code-block:: console

    gprof2dot ./callgrind.out.619277 -f callgrind > callgrind.dot
    dot -Tpng callgrind.dot -o callgrind.png


.. [vgtools] `<https://www.valgrind.org/info/tools.html>`_
