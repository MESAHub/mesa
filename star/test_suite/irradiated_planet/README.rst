.. _irradiated_planet:

*****************
irradiated_planet
*****************

This test case checks the evolution of an ~1 Mjup model after the surface has been irradiated.

This test case has 2 parts. Click to see a larger version of a plot.

* Part 1 (``inlist_create``) loads a 1 Mjup model, relaxes the surface to an optical depth of 10, irradiates the surface with a flux of 1e7 erg cm :sup:`-2` s :sup:`-1` to a column depth of 1 g cm :sup:`-2`, and terminates when the inflated envelope reaches a surface pressure of log10 (P/(erg cm :sup:`-3`)) = 6.0  (about 1 bar):

.. image:: ../../../star/test_suite/irradiated_planet/docs/pres_000084.svg
   :scale: 100%

* Part 2 (``inlist_evolves``) continues the evolution, without irradiation, for 15 Gyr. The irradiated planet model cools and shrinks with time:

.. image:: ../../../star/test_suite/irradiated_planet/docs/track1_000127.svg
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

    History_Track_win_flag(1) = .true.
    History_Track_win_width(1) = 12
    History_Track_win_aspect_ratio(1) = 0.75
    History_Track_title(1) = 'irradiated_planet'

    History_Track_xname(1) = 'log_star_age'
    History_Track_yname(1) = 'log_surf_cell_pressure'
    History_Track_xaxis_label(1) = 'log10 (Time/year) '
    History_Track_yaxis_label(1) = 'log10 (Surface Pressure/(erg/cm\u3\u))'
    History_Track_reverse_xaxis(1) = .false.
    History_Track_reverse_yaxis(1) = .false.
    History_Track_log_xaxis(1) = .false.
    History_Track_log_yaxis(1) = .false.

    History_Track_xmin(1) = 3.0
    History_Track_xmax(1) = 6.0
    History_Track_ymin(1) = 5.7
    History_Track_ymax(1) = 6.1

    History_Track_file_flag(1) = .true.
    History_Track_file_dir(1) = 'png'
    History_Track_file_prefix(1) = 'pres_'
    History_Track_file_interval(1) = 10000
    History_Track_file_width(1) = -1
    History_Track_file_aspect_ratio(1) = -1


    History_Track_win_flag(2) = .true.
    History_Track_win_width(2) = 12
    History_Track_win_aspect_ratio(2) = 0.75
    History_Track_title(2) = 'irradiated_planet'

    History_Track_xname(2) = 'log_star_age'
    History_Track_yname(2) = 'radius_cm'
    History_Track_xaxis_label(2) = 'Time (log years)'
    History_Track_yaxis_label(2) = 'Radius (cm)'
    History_Track_reverse_xaxis(2) = .false.
    History_Track_reverse_yaxis(2) = .false.
    History_Track_log_xaxis(2) = .false.
    History_Track_log_yaxis(2) = .false.

    History_Track_xmin(2) = 5.0
    History_Track_xmax(2) = 10.3
    History_Track_ymin(2) = 0.7e10
    History_Track_ymax(2) = 1.2e10

    History_Track_file_flag(2) = .true.
    History_Track_file_dir(2) = 'png'
    History_Track_file_prefix(2) = 'track1_'
    History_Track_file_interval(2) = 1000
    History_Track_file_width(2) = -1
    History_Track_file_aspect_ratio(2) = -1

 / ! end of pgstar namelist


Last-Updated: 14Jun2021 (MESA 5be9e57) by fxt.
