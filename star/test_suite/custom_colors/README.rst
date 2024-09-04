.. _custom_colors:

*************
custom_colors
*************

This test suite example shows how to use user-defined color filter and extinction files.

This test case has 1 part. Click to see a larger view of a plot.

* Part 1 (``inlist_1.0``) builds a 1.0 Msun, Z=0.02 metallicity, pre-main sequence model and evolves until core hydrogen depletion (mass fraction h1 < 0.1). This example loads the default |LCB98| color filter ''lcb98cor.dat'', a custom color filter ``data/blackbody_bc_v.txt`` which in this case is blackbody V band filter, and a custom extinction color correction file ``data/fake_av_v.txt``. Example color-color, color-magnitude, magnitude-color and magnitude-magnitude plots:

.. image:: ../../../star/test_suite/custom_colors/docs/Color_magnitude1_000241.svg
   :width: 100%

.. image:: ../../../star/test_suite/custom_colors/docs/Color_magnitude2_000241.svg
   :width: 100%

.. image:: ../../../star/test_suite/custom_colors/docs/Color_magnitude3_000241.svg
   :width: 100%


pgstar commands used for the first 7 plots:

.. code-block:: console

 &pgstar

   file_white_on_black_flag = .true. ! white_on_black flags -- true means white foreground color on black background
   file_device = 'png'            ! png
   file_extension = 'png'

   !file_device = 'vcps'          ! postscript
   !file_extension = 'ps'

    pgstar_interval = 1


 !# Color Magnitude Panels
   ! Plots either color-color, color-magnitude, magnitude-color or magnitude-magnitude

      !### Color_magnitude1

         Color_Magnitude_win_flag(1) = .true.

         Color_Magnitude_win_width(1) = 15
         Color_Magnitude_win_aspect_ratio(1) = 0.75 ! aspect_ratio = height/width

         Color_Magnitude_xleft(1) = 0.15
         Color_Magnitude_xright(1) = 0.85
         Color_Magnitude_ybot(1) = 0.15
         Color_Magnitude_ytop(1) = 0.85
         Color_Magnitude_txt_scale(1) = 1.0
         Color_Magnitude_title(1) = 'Color_magnitude1'

         ! setup default
         Color_Magnitude_num_panels(1) = 2

         ! Plots xaxis1-xaxis2 leave xaxis2 blank if you only want to plot xaxis1.
         Color_Magnitude_xaxis1_name(1) = 'model_number'
         Color_Magnitude_xaxis2_name(1) = ''


         ! Plots yaxis1-yaxis2 leave yaxis2 blank if you only want to plot yaxis1.
         Color_Magnitude_yaxis1_name(1)(1) = 'bc_B'
         Color_Magnitude_yaxis2_name(1)(1) = 'bc_U'
         Color_Magnitude_yaxis_reversed(1)(1) = .false.
         
         ! Plots `other_yaxis1-other_yaxis2` leave `other_yaxis2` blank if you only want to plot `other_yaxis1`.
         Color_Magnitude_other_yaxis1_name(1)(1) = 'abs_mag_V'
         Color_Magnitude_other_yaxis2_name(1)(1) = ''
         Color_Magnitude_other_yaxis_reversed(1)(1) = .true.


         Color_Magnitude_yaxis1_name(2)(1) = 'bc_B'
         Color_Magnitude_other_yaxis1_name(2)(1) = 'bc_U'
         
         ! Enables calling a subroutine to add extra information to a plot
         ! see `$MESA_DIR/star/other/pgstar_decorator.f90`
         Color_Magnitude_use_decorator(1) = .true.

         ! file output
         Color_Magnitude_file_flag(1) = .true.
         Color_Magnitude_file_dir(1) = 'png'
         Color_Magnitude_file_prefix(1) = 'Color_magnitude1_'
         Color_Magnitude_file_interval(1) = 5 ! output when `mod(model_number,Color_magnitude1_file_interval)==0`
         Color_Magnitude_file_width(1) = -1 ! (inches) negative means use same value as for window
         Color_Magnitude_file_aspect_ratio(1) = -1 ! negative means use same value as for window


      !### Color_magnitude2

         Color_Magnitude_win_flag(2) = .true.

         Color_Magnitude_win_width(2) = 15
         Color_Magnitude_win_aspect_ratio(2) = 0.75 ! aspect_ratio = height/width

         Color_Magnitude_xleft(2) = 0.15
         Color_Magnitude_xright(2) = 0.85
         Color_Magnitude_ybot(2) = 0.15
         Color_Magnitude_ytop(2) = 0.85
         Color_Magnitude_txt_scale(2) = 1.0
         Color_Magnitude_title(2) = 'Color_magnitude2'

         ! Plots xaxis1-xaxis2 leave xaxis2 blank if you only want to plot xaxis1.
         Color_Magnitude_xaxis1_name(2) = 'abs_mag_B'
         Color_Magnitude_xaxis2_name(2) = 'abs_mag_U'

         ! Plots yaxis1-yaxis2 leave yaxis2 blank if you only want to plot yaxis1.
         Color_Magnitude_yaxis1_name(1)(2) = 'abs_mag_R'
         Color_Magnitude_yaxis2_name(1)(2) = 'abs_mag_J'

         ! setup default
         Color_Magnitude_num_panels(2) = 1
         ! file output
         Color_Magnitude_file_flag(2) = .true.
         Color_Magnitude_file_dir(2) = 'png'
         Color_Magnitude_file_prefix(2) = 'Color_magnitude2_'
         Color_Magnitude_file_interval(2) = 5 ! output when `mod(model_number,Color_magnitude2_file_interval)==0`
         Color_Magnitude_file_width(2) = -1 ! (inches) negative means use same value as for window
         Color_Magnitude_file_aspect_ratio(2) = -1 ! negative means use same value as for window


      !### Color_magnitude3

         Color_Magnitude_win_flag(3) = .true.

         Color_Magnitude_win_width(3) = 15
         Color_Magnitude_win_aspect_ratio(3) = 0.75 ! aspect_ratio = height/width

         Color_Magnitude_xleft(3) = 0.15
         Color_Magnitude_xright(3) = 0.85
         Color_Magnitude_ybot(3) = 0.15
         Color_Magnitude_ytop(3) = 0.85
         Color_Magnitude_txt_scale(3) = 1.0
         Color_Magnitude_title(3) = 'Color_magnitude3'

         ! Plots xaxis1-xaxis2 leave xaxis2 blank if you only want to plot xaxis1.
         Color_Magnitude_xaxis1_name(3) = 'model_number'
         Color_Magnitude_xaxis2_name(3) = ''

         ! Plots yaxis1-yaxis2 leave yaxis2 blank if you only want to plot yaxis1.
         Color_Magnitude_yaxis1_name(1)(3) = 'bc_v_bb'
         
         Color_Magnitude_other_yaxis1_name(1)(3) = 'av_v'
         
         ! setup default
         Color_Magnitude_num_panels(3) = 1
         ! file output
         Color_Magnitude_file_flag(3) = .true.
         Color_Magnitude_file_dir(3) = 'png'
         Color_Magnitude_file_prefix(3) = 'Color_magnitude3_'
         Color_Magnitude_file_interval(3) = 5 ! output when `mod(model_number,Color_magnitude3_file_interval)==0`
         Color_Magnitude_file_width(3) = -1 ! (inches) negative means use same value as for window
         Color_Magnitude_file_aspect_ratio(3) = -1 ! negative means use same value as for window


 / ! end of pgstar namelist

.. |LCB98| replace:: `Lejeune, Cuisinier, & Buser (1998) <https://ui.adsabs.harvard.edu/abs/1998A%26AS..130...65L/abstract>`__

Last-Updated: 05Jun2021 (MESA 5be9e57) by fxt

