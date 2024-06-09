.. _wd_diffusion:

************
wd_diffusion
************

This test case the checks element diffusion in a 0.6 Msun carbon-oxygen white dwarf.

This test case has 1 part. Click to see a larger version of a plot.

* Part 1 (``inlist_wd_diff``) loads ``co_wd_settled.mod``, a 0.596 Msun carbon-oxygen white dwarf model with hydrogen-rich outer layers, produced from ``final.mod`` by :ref:`make_co_wd` as of version 14731. The run terminates when the effective temperatures drops below 9000 K. The ``run_star_extras.f90`` then checks that the electrostatic/gravitational force ratio at the 0.5 Msun point are between 1.95 and 2.05:

.. code-block:: console

 Core eE/mg =    2.0005029376686494     
 passed test for electric field in the core


.. image:: ../../../star/test_suite/wd_diffusion/docs/profile_000527.svg
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

        Profile_Panels_win_flag(2) = .true.
        Profile_Panels_win_width(2) = 12
        Profile_Panels_title(2) = 'wd_diffusion'

        Profile_Panels_yaxis_name(2, 2) = 'eE_div_mg_element_diffusion'

        Profile_Panels_xaxis_name(2) = 'logxm'
        Profile_Panels_xaxis_reversed(2) = .true.
        Profile_Panels_xmin(2) = -8
        Profile_Panels_xmax(2) = 0
        Profile_Panels_show_mix_regions_on_xaxis(2) = .true.

        Profile_Panels_xright(2) = 0.92
        Profile_Panels_ytop(2) = 0.92

        num_abundance_line_labels = 5
        Abundance_legend_max_cnt = 0

        Profile_Panels_file_flag(2) = .true.
        Profile_Panels_file_dir(2) = 'pgstar_out'
        Profile_Panels_file_prefix(2) = 'profile_'
        Profile_Panels_file_interval(2) = 100000
        Profile_Panels_file_width(2) = -1
        Profile_Panels_file_aspect_ratio(2) = -1

 / ! end of pgstar namelist



Last-Updated: 06Jul2021 (MESA 094ff71) by fxt.


.. # define a hard line break for HTML
.. |br| raw:: html

      <br>
