These plotting programs and scripts aim to make it simple to plot thermal neutrino energy loss rates.

To use them, first compile the neu plotter program::

  ./clean
  ./mk

The options that control what data is output by the program are
documented in the file ``inlist_plotter``.  This inlist contains two
namelists: the ``neu`` namelist that controls the MESA ``neu`` module
and the ``plotter`` namelist that controls the plotter program (see
source in ``src/neu_plotter.f90``).  Edit these namelists so that the
plotter will output the desired quantities.

Then, run the plotter::

  ./plotter

This will create an output data file ``neu_plotter.dat``.

A python script that knows how to read this file and plot it using
matplotlib is provided.  You can invoke it via::

  ./plotter.py

This will produce a plot file ``neu_plotter.png`` that you can then
view.  You may need to edit the python file to manually adjust various
aspects of the plotting (e.g., colorbar limits).




  


