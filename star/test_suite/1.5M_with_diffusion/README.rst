.. _1.5M_with_diffusion:

*******************
1.5M_with_diffusion
*******************

The test checks the functionality of element diffusion. 
The test vehicle is a 1.5 Msun solar metallicity model.


This test case has two parts. Click to see a larger view of a plot.

* Part 1 (``inlist_to_ZAMS``) creates the pre-main-sequence model and evolves it to the main sequence. In the mass fraction plot, note the roughness of the 12C mass fraction profile (orange curve) from convection and that the surface hydrogen profile is flat (yellow curve). 

.. image:: ../../../star/test_suite/1.5M_with_diffusion/docs/zams.svg
   :scale: 100%

* Part 2 (``1.5M_with_diffusion``) continues the evolution until the central hydrogen mass fraction drops below 0.01. In the mass fraction plot, note the 12C mass fraction profile (orange) has been smoothed by element diffusion and that the surface is nearly pure hydrogen (yellow upwards spike) from heavier elements diffusing inwards (e.g., orange downwards spike in 16O).

.. image:: ../../../star/test_suite/1.5M_with_diffusion/docs/h_depletion.svg
   :scale: 100%


pgstar commands used for the plots above:

.. code-block:: console

 &pgstar

     pgstar_interval = 10

  ! device

   file_white_on_black_flag = .true. ! white_on_black flags -- true means white foreground color on black background
   !file_device = 'png'            ! png
   !file_extension = 'png'           

   file_device = 'vcps'          ! postscript
   file_extension = 'ps'           


      Grid2_win_width = 15
      Grid2_win_flag = .true.
      Grid2_file_flag = .true.
      file_digits = 7
      Grid2_file_dir = 'png2'
      Grid2_file_prefix = 'grid'
      Grid2_file_interval = 5
      Grid2_file_width = 15

 / ! end of pgstar namelist




Last-Updated: 27May2021 (MESA ebecc10) by fxt

