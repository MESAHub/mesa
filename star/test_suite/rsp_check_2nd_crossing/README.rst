.. _rsp_check_2nd_crossing:

**********************
rsp_check_2nd_crossing
**********************

This test case exercises the RSP model building and linear nonadiabatic stability analysis 
to find the instability strip edges, and effective temperatures offset from the blue edge of the instability strip.

This test case has 1 part. 

* Part 1 (``inlist_rsp_check_2nd_crossing``) opens the MESA history file ``history_logs/history_7M.data`` for a 7 Msun, Z=0.008 constructed from the :ref:`5M_cepheid_blue_loop` test case. Using a significant ``run_star_extras.f90``, the RSP initial model building and linear nonadiabatic stability analysis modules (see Section 2.2 of |MESA V| for details) are exercised to find the edges of the instability strip. The run then continues by interpolating the HR diagram evolutionary track to find points at selected effective temperatures Teff offset from the blue edge of the instabilty strip. The results are written to the terminal and a user-specifed file, ``7B_2nd_crossing.data`` in this case:

.. code-block:: console

       Teff_red_edge     5046.6949704733
      Teff_blue_edge     5997.8180161307

    offset                Teff                   L
         0     5997.8180161307     5576.7666901507
       100     5897.8180161307     5557.6216073044
       200     5797.8180161307     5539.0539422822
       300     5697.8180161307     5520.5578751096
       400     5597.8180161307     5502.2654214400
       500     5497.8180161307     5483.2779156791
       600     5397.8180161307     5463.2566431324
       700     5297.8180161307     5441.8094957682
       800     5197.8180161307     5418.5159995785
       900     5097.8180161307     5392.7973651162
       951     5046.6949704733     5378.3439119466

 done rsp_check_2nd_crossing
 7B_2nd_crossing.data

Last-Updated: 30Jun2021 (MESA 094ff71) by fxt.
