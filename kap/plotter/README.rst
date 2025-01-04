These plotting programs and scripts aim to make it simple to plot opacities.

The options that control what data is output by the program are documented in
the file ``inlist_plotter``.  This inlist contains two namelists: the ``kap``
namelist that controls the MESA ``kap`` module and the ``plotter`` namelist
that controls the plotter program (see source in ``src/kap_plotter.f90``). Edit
these namelists so that the plotter will output the desired quantities.

Then, run the plotter::

  make plot

This will produce a plot file ``kap_plotter.png`` that you can then view. You
may need to edit the python file to manually adjust various aspects of the
plotting (e.g., colorbar limits).
