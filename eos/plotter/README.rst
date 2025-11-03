These plotting programs and scripts aim to make it simple to plot EOS
quantities.

The options that control what data is output by the program are
documented in the file ``inlist_plotter``. This inlist contains two
namelists: the ``eos`` namelist that controls the MESA ``eos`` module
and the ``plotter`` namelist that controls the plotter program (see
source in ``src/eos_plotter.f90``). Edit these namelists so that the
plotter will output the desired quantities.

Then, run the plotter::

  make plot
