.. _wd_he_shell_ignition:

********************
wd_he_shell_ignition
********************

This test case checks the ignition of a helium layer in an accreting in a 0.96 Msun carbon-oxygen white dwarf model.

This test case has 1 part. Click to see a larger version of a plot.

* Part 1 (``inlist_he_shell_ignition``) loads 'co_wd_0.96M.mod', a 0.96 Msun carbon-oxygen white dwarf with a helium layer and then a hydrogen layer. The model is from a 6 Msun ZAMS progenitor with :ref:`make_co_wd` in r13718.  Pure helum then acceretes onto the hydrogen surface at a rate of 1e-9 Msun/yr. The helium layer ignites as teh total mass approches 1.36 Msun, and the run terminates when the helium burning luminosity first exceeds 1e11 Lsun:


.. image:: ../../../star/test_suite/wd_he_shell_ignition/docs/track1_000922.svg
   :width: 100%

|br|
The initial hydrogen burns as its pushed deeper into the model by accreretion. The initial helium layer and the accreted helium 
thus merge to form a single helium layer. This single helium layer ignites at it is pushed to higher densities and temperatures
by the continued accretion:

.. image:: ../../../star/test_suite/wd_he_shell_ignition/docs/profile_000922.svg
   :width: 100%

|br|
Temperature and desnity profile at ignition:

.. image:: ../../../star/test_suite/wd_he_shell_ignition/docs/trho_000922.svg
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
    Profile_Panels_win_width(2) = 10
    Profile_Panels_title(2) = 'wd_he_shell_ignition'

    Profile_Panels_xaxis_name(2) = 'mass'
    Profile_Panels_xaxis_reversed(2) = .false.
    Profile_Panels_xmin(2) = 0.90
    Profile_Panels_xmax(2) = 1.2
    Profile_Panels_show_mix_regions_on_xaxis(2) = .false.

    Profile_Panels_xright(2) = 0.92
    Profile_Panels_ytop(2) = 0.92

    num_abundance_line_labels = 5
    Abundance_legend_max_cnt = 0

    Profile_Panels_yaxis_name(2, 2) = 'Power'
    Profile_Panels_ymin(2, 2) = 5.0
    Profile_Panels_ymax(2, 2) = 15.0

    Profile_Panels_file_flag(2) = .true.
    Profile_Panels_file_dir(2) = 'pgstar_out'
    Profile_Panels_file_prefix(2) = 'profile_'
    Profile_Panels_file_interval(2) = 100000
    Profile_Panels_file_width(2) = -1
    Profile_Panels_file_aspect_ratio(2) = -1


    History_Track_win_flag(1) = .true.
    History_Track_win_width(1) = 12
    History_Track_win_aspect_ratio(1) = 0.75
    History_Track_title(1) = 'wd_he_shell_ignition'

    History_Track_xname(1) = 'log_star_age'
    History_Track_yname(1) = 'log_Lnuc'
    History_Track_xaxis_label(1) = 'log10(star_age/yr)'
    History_Track_yaxis_label(1) = 'log Lnuc/L\d\(2281)'
    History_Track_reverse_xaxis(1) = .false.
    History_Track_reverse_yaxis(1) = .false.

    History_Track_xmin(1) = 6.0
    History_Track_xmax(1) = 9.0
    History_Track_ymin(1) = -11.0
    History_Track_ymax(1) = 12.0

    History_Track_file_flag(1) = .true.
    History_Track_file_dir(1) = 'pgstar_out'
    History_Track_file_prefix(1) = 'track1_'
    History_Track_file_interval(1) = 10000
    History_Track_file_width(1) = -1
    History_Track_file_aspect_ratio(1) = -1


    TRho_Profile_win_flag = .true.
    TRho_Profile_win_width = 10
    TRho_Profile_win_aspect_ratio = 0.75 ! aspect_ratio = height/width
    TRho_Profile_title = 'wd_he_shell_ignition'      

    TRho_Profile_xmin = -6.0
    TRho_Profile_xmax = 10.0
    TRho_Profile_ymin = 5.0
    TRho_Profile_ymax = 9.0        
         
    TRho_Profile_xleft = 0.10
    TRho_Profile_xright = 0.93
    TRho_Profile_ybot = 0.10
    TRho_Profile_ytop = 0.90
    TRho_Profile_txt_scale = 0.9
         
    show_TRho_Profile_legend = .true.
    TRho_Profile_legend_coord = 0.07
    TRho_Profile_legend_fjust = 0.0
    TRho_Profile_legend_disp1 = -2.0
    TRho_Profile_legend_del_disp = -1.3
    TRho_Profile_legend_txt_scale = 0.9

    show_TRho_Profile_eos_regions = .true.
    show_TRho_Profile_degeneracy_line = .true.
    show_TRho_Profile_Pgas_Prad_line = .true.
    show_TRho_Profile_burn_lines = .true.
    show_TRho_Profile_burn_labels = .true.
      
    show_TRho_Profile_mass_locs = .true.
    num_profile_mass_points = 2 

    profile_mass_point_q(1) = 0.5
    profile_mass_point_color_index(1) = 1
    profile_mass_point_symbol(1) = -6
    profile_mass_point_symbol_scale(1) = 1.0
    profile_mass_point_str(1) = '  q=0.5'
    profile_mass_point_str_clr(1) = 1
    profile_mass_point_str_scale(1) = 0.8
         
    profile_mass_point_q(2) = 0.99
    profile_mass_point_color_index(2) = 1
    profile_mass_point_symbol(2) = -6
    profile_mass_point_symbol_scale(2) = 1.0
    profile_mass_point_str(2) = '  q=0.99'
    profile_mass_point_str_clr(2) = 1
    profile_mass_point_str_scale(2) = 0.8
         
    TRho_Profile_file_flag = .true.
    TRho_Profile_file_dir = 'pgstar_out'
    TRho_Profile_file_prefix = 'trho_'
    TRho_Profile_file_interval = 10000
    TRho_Profile_file_width = -1 
    TRho_Profile_file_aspect_ratio = -1 

 / ! end of pgstar namelist



Last-Updated: 07Jul2021 (MESA 094ff71) by fxt.


.. # define a hard line break for HTML
.. |br| raw:: html

      <br>
