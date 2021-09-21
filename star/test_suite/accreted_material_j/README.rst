.. _accreted_material_j:

*******************
accreted_material_j
*******************

This test suite example checks the accretion of material and angular momentum onto a 20 Msun model.


This test case has two parts. Click to see a larger view of a plot.

* Part 1 (``inlist_zams``) creates a 20 Msun, Z=0.02 metallicity, main-sequence model.

* Part 2 (``inlist_accreted_material_j``) continues the evolution by first applying uniform rotation Omega/Omega_crit = 0.1, and then accreting material with the same composition as the outermost cell at a rate of 0.001 Msun/year with an angular momentum of 0.1 Keplerian. The model terminates when the mass reaches 25 Msun:

.. image:: ../../../star/test_suite/accreted_material_j/docs/track1_000452.svg
   :scale: 100%


pgstar commands used for the plots above:

.. code-block:: console

 &pgstar

   file_white_on_black_flag = .true. ! white_on_black flags -- true means white foreground color on black background
   file_device = 'png'            ! png
   file_extension = 'png'

   !file_device = 'vcps'          ! postscript
   !file_extension = 'ps'

    pgstar_interval = 10


    History_Track1_win_flag = .true.
    History_Track1_win_width = 12
    History_Track1_title = 'accreted_material_j'                                                                                                                                                                                   
    History_Track1_xname = 'star_mass'
    History_Track1_yname = 'log_total_angular_momentum'
    History_Track1_yaxis_label = 'log J'
    History_Track1_xaxis_label = 'M\d\(2281)'
    History_Track1_reverse_xaxis = .false.
    History_Track1_reverse_yaxis = .false.

    History_Track1_xmin = 20.0
    History_Track1_xmax = 25.0
    History_Track1_ymin = 51.9
    History_Track1_ymax = 52.7

 ! file output
    History_Track1_file_flag = .true.
    History_Track1_file_dir = 'png'
    History_Track1_file_prefix = 'track1_'
    History_Track1_file_interval = 10
    History_Track1_file_width =12
    History_Track1_file_aspect_ratio = -1

 / ! end of pgstar namelist



Last-Updated: 30May2021 (MESA e37f76f) by fxt

