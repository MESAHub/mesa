.. _hse_riemann:

***********
hse_riemann
***********

This test case checks Riemann HLLC solver can hold an envelope model in hydrostatic equilibrium.

This test case has 1 parts. Click to see a larger version of a plot.

* Part 1 (``inlist_finish``) loads a model containing the envelope from a 12 msun model with its core excised. The model is then evolved for 1e5 s, and the ``run_star_extras.f90`` then checks that the peak mach number is less than 0.1. This is an example of a null test; if the Riemann HLLC solver is working, then the profiles should remain static and the peak mach number should be small:

.. raw:: html

  <video width="640" height="280" controls>
    <source src="../../../../star/test_suite/hse_riemann/docs/profile1.mp4" type="video/mp4">
    Your browser does not support the video tag.
  </video>


pgstar commands used for the plots that ake the movie above:


.. code-block:: console

 &pgstar

  file_white_on_black_flag = .true. ! white_on_black flags -- true means white foreground color on black background
  !file_device = 'png'            ! png
  !file_extension = 'png'

  file_device = 'vcps'          ! postscript
  file_extension = 'ps'

   pgstar_interval = 2

         Profile_Panels1_win_flag = .true.
         Profile_Panels1_win_width = 12
         Profile_Panels1_win_aspect_ratio = 1.0
         
         Profile_Panels1_xaxis_name = 'logR_cm'
         Profile_Panels1_xmin = -101d0 
         Profile_Panels1_xmax = -101d0 

         Profile_Panels1_num_panels = 3

         Profile_Panels1_yaxis_name(1) = 'logRho'         
            Profile_Panels1_ymin(1) = -101d0 
            Profile_Panels1_ymax(1) = -101d0 
            Profile_Panels1_ymargin(1) = 0.1
      
         Profile_Panels1_other_yaxis_name(1) = 'entropy' 
            Profile_Panels1_other_ymin(1) = -101d0 
            Profile_Panels1_other_ymax(1) = -101d0 
            Profile_Panels1_other_ymargin(1) = 0.1
          
         Profile_Panels1_yaxis_name(2) = 'logT'
            Profile_Panels1_ymin(2) = -101d0 
            Profile_Panels1_ymax(2) = -101d0 
            Profile_Panels1_ymargin(2) = 0.1
         
         Profile_Panels1_yaxis_name(3) = 'v_div_csound'
            Profile_Panels1_ymin(3) = -0.06  
            Profile_Panels1_ymax(3) = 0.06 
            Profile_Panels1_ymargin(3) = 0.1

  Profile_Panels1_file_flag = .true.
  Profile_Panels1_file_dir = 'pgstar_out'
  Profile_Panels1_file_prefix = 'profile1_'
  Profile_Panels1_file_interval = 5     ! output when mod(model_number,file_interval)==0
  Profile_Panels1_file_width = -1        ! (inches) negative means use same value as for window
  Profile_Panels1_file_aspect_ratio = -1 ! negative means use same value as for window

 / ! end of pgstar namelist


Last-Updated: 14Jun2021 (MESA 5be9e57) by fxt.
