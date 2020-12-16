====================
mesa/star Quickstart
====================

..
  If you can't stand reading anything that isn't on the web, skip
  this and go directly to https://docs.mesastar.org.  Even if you do
  read this file, when you are done you should still go to that site!

These directions assume you have already installed MESA.

Copy this work directory to somewhere outside the mesa directory tree
and name it anything you like.

Compile by executing the ``mk`` script::

    ./mk

This compiles the files in ``src/``, links them against MESA, and
produces the ``star`` executable file.

By convention, run the program using the  ``rn`` script::

    ./rn

When MESA runs, it first reads the ``inlist`` file.  This file can
point to other inlist files.  Here, it points to ``inlist_project``
and ``inlist_pgstar``.  You might want to change the file name of
``inlist_project`` to something more appropriate for your work, such
as ``inlist_hot_jupiter`` -- if you do, then change the names in
``inlist`` to match.

You can control MESA by editing the options in the various sectinons
of the inlist.  The full set of parameters and their default values
can be found in the defaults files listed in the inlists.

To restart MESA from a saved photo file, run the program using the
``re`` script::

    ./re [photo]

where ``[photo]`` is one of the saved snapshot files in ``photos``.
If no file is specified, MESA restarts from the most recent photo.

You can remove the compiled files and the star program by executing
the ``clean`` script::

    ./clean


