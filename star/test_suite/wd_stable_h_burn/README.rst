.. _wd_stable_h_burn:

****************
wd_stable_h_burn
****************

This test case checks the evolution stable hydrogen burning on a white dwarf.

This test case has 1 parts. Click to see a larger version of a plot.

* Part 1 (``inlist_wd_stable_h_burn``) loads ``1.1M_lgTc_7.7.mod``, a prebuilt 1.1 Msun carbon oxygen white dwarf from the :ref:`make_co_wd` test suite in r13738 with an an ``initial_mass`` of 6.4 Msun. The mass is relaxed to 1.0 Msun, the optical depth is relaxed to 300, and a hydrogen-rich composition is accreted at 2.5e7 Msun/yr. After about 100 year of evolution there is a hydrogen burning outburst followed by a period of stable hydrogen burning:


.. image:: ../../../star/test_suite/wd_stable_h_burn/docs/grid_000829.svg
   :width: 100%


|br|
pgstar commands used for the plots above:

.. code-block:: console

 &pgstar

  file_white_on_black_flag = .true. ! white_on_black flags -- true means white foreground color on black background
  !file_device = 'png'            ! png
  !file_extension = 'png'

  file_device = 'vcps'          ! postscript
  file_extension = 'ps'

          pgstar_interval = 1

   pgstar_left_yaxis_label_disp = 4.0

   History_Track_win_flag(1) = .false.
   History_Track_win_width(1) = 12
   History_Track_win_aspect_ratio(1) = 0.75
   History_Track_txt_scale(1) = 0.8
   History_Track_title(1) = ' '

   History_Track_xname(1) = 'star_age'
   History_Track_yname(1) = 'total_mass_h1'
   History_Track_xaxis_label(1) = ' '
   History_Track_yaxis_label(1) = 'H1 Mass (M\d\(2281)\u)'
   History_Track_reverse_xaxis(1) = .false.
   History_Track_reverse_yaxis(1) = .false.

   History_Track_xmin(1) = 0.0
   History_Track_xmax(1) = 420.0
   History_Track_ymin(1) = 0.0
   History_Track_ymax(1) = 2.0e-5

   History_Track_file_flag(1) = .false.
   History_Track_file_dir(1) = 'pgstar_out'
   History_Track_file_prefix(1) = 'track1_'
   History_Track_file_interval(1) = 10000
   History_Track_file_width(1) = -1
   History_Track_file_aspect_ratio(1) = -1


   History_Track_win_flag(2) = .false.
   History_Track_win_width(2) = 12
   History_Track_win_aspect_ratio(2) = 0.75
   History_Track_txt_scale(2) = 0.8
   History_Track_title(2) = ' '

   History_Track_xname(2) = 'star_age'
   History_Track_yname(2) = 'log_LH'
   History_Track_xaxis_label(2) = ' '
   History_Track_yaxis_label(2) = 'log10 (LH/L\d\(2281)\u)'
   History_Track_reverse_xaxis(2) = .false.
   History_Track_reverse_yaxis(2) = .false.

   History_Track_xmin(2) = 0
   History_Track_xmax(2) = 420.0
   History_Track_ymin(2) = 0
   History_Track_ymax(2) = 8.0

   History_Track_file_flag(2) = .false.
   History_Track_file_dir(2) = 'pgstar_out'
   History_Track_file_prefix(2) = 'track2_'
   History_Track_file_interval(2) = 10000
   History_Track_file_width(2) = -1
   History_Track_file_aspect_ratio(2) = -1


   History_Track_win_flag(3) = .false.
   History_Track_win_width(3) = 12
   History_Track_win_aspect_ratio(3) = 0.75
   History_Track_txt_scale(3) = 0.8
   History_Track_title(3) = ' '

   History_Track_xname(3) = 'star_age'
   History_Track_yname(3) = 'log_R'
   History_Track_xaxis_label(3) = 'Time (yr)'
   History_Track_yaxis_label(3) = 'log10 (R/R\d\(2281)\u)'
   History_Track_reverse_xaxis(3) = .false.
   History_Track_reverse_yaxis(3) = .false.

   History_Track_xmin(3) = 0
   History_Track_xmax(3) = 420.0
   History_Track_ymin(3) = -2.5
   History_Track_ymax(3) = 0.2

   History_Track_file_flag(3) = .false.
   History_Track_file_dir(3) = 'pgstar_out'
   History_Track_file_prefix(3) = 'track3_'
   History_Track_file_interval(3) = 10000
   History_Track_file_width(3) = -1
   History_Track_file_aspect_ratio(3) = -1

   Grid_win_flag(1) = .true.
   Grid_win_width(1) = 10
   Grid_win_aspect_ratio(1) = 1.2

   Grid_plot_name(1, :) = ''
   Grid_txt_scale_factor(1, :) = 1.0 ! multiply txt_scale for subplot by this

   Grid_title(1) = 'wd_stable_h_burn'

   Grid_num_cols(1) = 1 ! divide plotting region into this many equal width cols
   Grid_num_rows(1) = 3 ! divide plotting region into this many equal height rows
   Grid_num_plots(1) = 3 ! <= 10

   Grid_plot_name(1, 1) = 'History_Track1'
   Grid_plot_row(1, 1) = 1           ! number from 1 at top
   Grid_plot_rowspan(1, 1) = 1       ! plot spans this number of rows
   Grid_plot_col(1, 1) =  1          ! number from 1 at left
   Grid_plot_colspan(1, 1) = 1       ! plot spans this number of columns

   Grid_plot_pad_left(1, 1) = 0.00    ! fraction of full window width for padding on left
   Grid_plot_pad_right(1, 1) = -0.02   ! fraction of full window width for padding on right
   Grid_plot_pad_top(1, 1) = -0.02     ! fraction of full window height for padding at top
   Grid_plot_pad_bot(1, 1) = 0.05     ! fraction of full window height for padding at bottom
   Grid_txt_scale_factor(1, 1) = 0.8 ! multiply txt_scale for subplot by this


   Grid_plot_name(1, 2) = 'History_Track2'
   Grid_plot_row(1, 2) = 2           ! number from 1 at top
   Grid_plot_rowspan(1, 2) = 1       ! plot spans this number of rows
   Grid_plot_col(1, 2) =  1          ! number from 1 at left
   Grid_plot_colspan(1, 2) = 1       ! plot spans this number of columns

   Grid_plot_pad_left(1, 2) = 0.0    ! fraction of full window width for padding on left
   Grid_plot_pad_right(1, 2) = -0.02   ! fraction of full window width for padding on right
   Grid_plot_pad_top(1, 2) = 0.0     ! fraction of full window height for padding at top
   Grid_plot_pad_bot(1, 2) = 0.00     ! fraction of full window height for padding at bottom
   Grid_txt_scale_factor(1, 2) = 0.8 ! multiply txt_scale for subplot by this


   Grid_plot_name(1, 3) = 'History_Track3'
   Grid_plot_row(1, 3) = 3           ! number from 1 at top
   Grid_plot_rowspan(1, 3) = 1       ! plot spans this number of rows
   Grid_plot_col(1, 3) =  1          ! number from 1 at left
   Grid_plot_colspan(1, 3) = 1       ! plot spans this number of columns

   Grid_plot_pad_left(1, 3) = 0.0    ! fraction of full window width for padding on left
   Grid_plot_pad_right(1, 3) = -0.02   ! fraction of full window width for padding on right
   Grid_plot_pad_top(1, 3) = 0.05     ! fraction of full window height for padding at top
   Grid_plot_pad_bot(1, 3) = -0.03     ! fraction of full window height for padding at bottom
   Grid_txt_scale_factor(1, 3) = 0.8 ! multiply txt_scale for subplot by this


   Grid_file_flag(1) = .true.
   Grid_file_dir(1) = 'pgstar_out'
   Grid_file_prefix(1) = 'grid_'
   Grid_file_interval(1) = 10000
   Grid_file_width(1) = -1
   Grid_file_aspect_ratio(1) = -1

 / ! end of pgstar namelist



Last-Updated: 08Jul2021 (MESA 094ff71) by fxt.


.. # define a hard line break for HTML
.. |br| raw:: html

      <br>
