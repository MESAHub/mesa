.. _make_brown_dwarf:

****************
make_brown_dwarf
****************

This test case checks the creation of a 1.05 Mjup, Z=1e-4 metallicity model and its subsequent evolution for 20 billion years.

This test case has 1 part. Click to see a larger version of a plot.

* Part 1 (``inlist_make_brown_dwarf``) creates a 1.05 Mjup, Z=1e-4 metallicity model, using the less frequently deployed ``create_initial_model`` capability, and evolves the brown drawf model for 20 billion years.

.. image:: ../../../star/test_suite/make_brown_dwarf/docs/grid1_000237.svg
   :scale: 100%

pgstar commands used for the plots above:


.. code-block:: console

 &pgstar

  file_white_on_black_flag = .true. ! white_on_black flags -- true means white foreground color on black background
  !file_device = 'png'            ! png
  !file_extension = 'png'

  file_device = 'vcps'          ! postscript
  file_extension = 'ps'

   pgstar_interval = 10

  Grid_win_flag(1) = .true.
  Grid_win_width(1) = 14
  Grid_win_aspect_ratio(1) = 0.75
  Grid_title(1) = 'inlist_make_brown_dwarf'
      
  TRho_Profile_xmin = -8
  TRho_Profile_xmax = 4
  TRho_Profile_ymin = 2
  TRho_Profile_ymax = 6

  Grid_file_flag(1) = .true.
  Grid_file_dir(1) = 'png'
  Grid_file_prefix(1) = 'grid1_'
  Grid_file_interval(1) = 10000 ! output when mod(model_number,file_interval)==0
  Grid_file_width(1) = 14        ! (inches) negative means use same value as for window
  Grid_file_aspect_ratio(1) = -1 ! negative means use same value as for window

 / ! end of pgstar namelist


Last-Updated: 17Jun2021 (MESA e2acbc2) by fxt.
