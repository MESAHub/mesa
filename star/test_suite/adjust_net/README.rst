.. _adjust_net:

**********
adjust_net
**********

This test suite example checks the functionality of the adaptive nuclear reaction network.

Physical checks
===============

This test case tracks the number of isotopes it has at the end of the run.

Inlists
=======

This test case has two parts. Click to see a larger view of a plot.

* Part 1 (``inlist_zams``) creates a 15 Msun, Z=0.02 metallicity, main-sequence model using the default 8 isotope ``basic.net``.

* Part 2 (``inlist_adjust_net_header``) continues the evolution, activates the adaptive nuclear reaction network, and terminates at about 2 million years with 62 isotopes in the reaction network:

.. image:: ../../../star/test_suite/adjust_net/docs/track1_000234.svg
   :width: 100%

.. image:: ../../../star/test_suite/adjust_net/docs/Network_000205.svg
   :width: 100%

.. image:: ../../../star/test_suite/adjust_net/docs/Network_000234.svg
   :width: 100%


pgstar commands used for the plots above:

.. code-block:: console

 &pgstar

   file_white_on_black_flag = .true. ! white_on_black flags -- true means white foreground color on black background
   file_device = 'png'            ! png
   file_extension = 'png'

   !file_device = 'vcps'          ! postscript
   !file_extension = 'ps'

    pgstar_interval = 10


    Network_win_flag = .true.
    Network_win_width = 12
    Network_win_aspect_ratio = 0.75
    Network_title = 'adjust_net'

    Network_nmin = -101d0
    Network_nmax = 20.0
    Network_zmin = -101d0
    Network_zmax = 20.0

    Network_show_mass_fraction = .true.
    Network_show_element_names = .true.
    Network_show_colorbar = .true.

    Network_log_mass_frac_min = -5.0d0
    Network_log_mass_frac_max = 0.0d0

    Network_file_flag = .true.
    Network_file_dir = 'png'
    Network_file_prefix = 'Network_'
    Network_file_interval = 10
    Network_file_width = 12
    Network_file_aspect_ratio = -1


    History_Track_win_flag(1) = .true.
    History_Track_win_width(1) = 12
    History_Track_title(1) = 'adjust_net'
    History_Track_xname(1) = 'model_number'
    History_Track_xaxis_label(1) = 'model number'

    History_Track_yname(1) = 'species'
    History_Track_yaxis_label(1) = 'species in network'

    History_Track_reverse_xaxis(1) = .false.
    History_Track_reverse_yaxis(1) = .false.

    History_Track_xmin(1) = 180.0
    History_Track_xmax(1) = 240.0
    History_Track_ymin(1) = 0.0
    History_Track_ymax(1) = 70.0

    History_Track_file_flag(1) = .true.
    History_Track_file_dir(1) = 'png'
    History_Track_file_prefix(1) = 'track1_'
    History_Track_file_interval(1) = 10
    History_Track_file_width(1) =12
    History_Track_file_aspect_ratio(1) = -1

 / ! end of pgstar namelist



Last-Updated: 31May2021 (MESA e37f76f) by fxt

