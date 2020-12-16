These plotting programs and scripts aim to make it simple to plot EOS
quantities.

To use them, first compile the eos plotter program::

  ./clean
  ./mk

The options that control what data is output by the program are
documented in the file ``inlist_plotter``.  This inlist contains two
namelists: the ``eos`` namelist that controls the MESA ``eos`` module
and the ``plotter`` namelist that controls the plotter program (see
source in ``src/eos_plotter.f90``).  Edit these namelists so that the
plotter will output the desired quantities.

Then, run the plotter::

  ./plotter

This will create an output data file ``eos_plotter.dat``.

A python script that knows how to read this file and plot it using
matplotlib is provided.  You can invoke it via::

  ./plotter.py

This will produce a plot file ``eos_plotter.png`` that you can then
view.  You may need to edit the python file to manually adjust various
aspects of the plotting (e.g., colorbar limits).




  


