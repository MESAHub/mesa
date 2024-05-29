.. _make_env:

********
make_env
********

This test case checks the creation and stability of a pure iron neutron star envelope.

This test case has 1 part. Click to see a larger version of a plot.

* Part 1 (``inlist_env_header``) first creates an initial neutron star envelope through its significant ``run_star_extras.f90`` and saves the result as ``start.mod``.  This initial model is then loaded and the evolution begins. The envelope should remain stable - no changes in the thermodynamic or structure profiles - over the 10 million year evolution.

.. image:: ../../../star/test_suite/make_env/docs/profile1_000172.svg
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

        Profile_Panels_win_flag(1) = .true.
        Profile_Panels_win_width(1) = 12
        Profile_Panels_win_aspect_ratio(1) = 1.0
        Profile_Panels_txt_scale(1) = 0.8

        Profile_Panels_xaxis_name(1) = 'logxm' ! 'log_q' ! 'logR_cm'
        Profile_Panels_xmin(1) = -101d0
        Profile_Panels_xmax(1) = -101d0
        Profile_Panels_xaxis_reversed(1) = .true.

        Profile_Panels_num_panels(1) = 3

        Profile_Panels_yaxis_name(1, 1) = 'logRho'
           Profile_Panels_ymin(1, 1) = -101d0
           Profile_Panels_ymax(1, 1) = -101d0
           Profile_Panels_ymargin(1, 1) = 0.1

        Profile_Panels_other_yaxis_name(1, 1) = 'entropy'
           Profile_Panels_other_ymin(1, 1) = -101d0
           Profile_Panels_other_ymax(1, 1) = -101d0
           Profile_Panels_other_ymargin(1, 1) = 0.1

        Profile_Panels_yaxis_name(1, 2) = 'logT'
           Profile_Panels_ymin(1, 2) = -101d0
           Profile_Panels_ymax(1, 2) = -101d0
           Profile_Panels_ymargin(1, 2) = 0.1

        Profile_Panels_yaxis_name(1, 3) = 'v_div_csound'
           Profile_Panels_ymin(1, 3) = -0.01
           Profile_Panels_ymax(1, 3) = 0.01
           Profile_Panels_ymargin(1, 3) = 0.1

   Profile_Panels_file_flag(1) = .true.
   Profile_Panels_file_dir(1) = 'pgstar_out'
   Profile_Panels_file_prefix(1) = 'profile1_'
   Profile_Panels_file_interval(1) = 1000
   Profile_Panels_file_width(1) = -1
   Profile_Panels_file_aspect_ratio(1) = -1

 / ! end of pgstar namelist


Last-Updated: 17Jun2021 (MESA e2acbc2) by fxt.
