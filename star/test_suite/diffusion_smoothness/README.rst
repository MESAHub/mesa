.. _diffusion_smoothness:

********************
diffusion_smoothness
********************

This test suite case checks that element diffusion produces a sufficiently smooth Brunt profile.

This test case has 2 part2. Click to see a larger view of a plot.

* Part 1 (``inlist_zams``) builds a 1.0 Msun, Z=0.02 metallicity, main-sequence model.

* Part 2 (``inlist_diffusion_smoothness``) activates element diffusion and continues the evolution for 3e9 years. The ``run-star_extras.f90`` then checks that the Brunt profile is sufficiently smooth in the central portion of the model.

.. image:: ../../../star/test_suite/diffusion_smoothness/docs/brunt_000036.svg
   :scale: 100%

pgstar commands used for the plot:

.. code-block:: console

 &pgstar

   file_white_on_black_flag = .true. ! white_on_black flags -- true means white foreground color on black background
   file_device = 'png'            ! png
   file_extension = 'png'

   !file_device = 'vcps'          ! postscript
   !file_extension = 'ps'

    pgstar_interval = 1

 ! Brunt B profile
    Profile_Panels_win_flag(1) = .true.
    Profile_Panels_win_width(1) = 14
    Profile_Panels_show_grid(1) = .false
    Profile_Panels_title(1) = 'diffusion_smoothness'

    Profile_Panels_xaxis_name(1) = 'zone'
    Profile_Panels_xmin(1) = -101d0
    Profile_Panels_xmax(1) = -101d0
    Profile_Panels_xaxis_reversed(1) = .true.
    Profile_Panels_show_mix_regions_on_xaxis(1) = .true.

    Profile_Panels_num_panels(1) = 1
    Profile_Panels_yaxis_name(1, 1) = 'brunt_B'
    Profile_Panels1_ymin(:) = -101d0
    Profile_Panels1_ymax(:) = -101d0
    Profile_Panels_other_yaxis_name(1, 1) = ''

    Profile_Panels_file_flag(1) = .true.
    Profile_Panels_file_dir(1) = 'png'
    Profile_Panels_file_prefix(1) = 'brunt_'
    Profile_Panels_file_interval(1) = 10000
    Profile_Panels_file_width(1) = 14

 / ! end of pgstar namelist

Last-Updated: 06Jun2021 (MESA 5be9e57) by fxt

