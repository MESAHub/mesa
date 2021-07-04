.. _timing:

******
timing
******

This test checks the counter and timing routines with a 1.5 Msun, Z=0.02 metallicity model.

This test case has 2 part2. 

* Part 1 (``inlist_zams``) builds a 1.5 Msun, Z=0.02 metallicity and evolves it to the zero-age main sequence.

* Part 2 (``inlist_timing``) continues the evolutioon until the central mass fraction of hydrogen drops below 0.5. To get a counter and timing breakdown the  star_job namelist parameter ``first_model_for_timing=2`` is set. At the end of the run, various counters and timings are written to the terminal:

.. code-block:: console

                                                nz               853
                                        nvar_total                12
                                 basic.net species                 8
                       total_num_solver_iterations               298
                          timing_num_get_eos_calls            719175
                        timing_num_solve_eos_calls              1299
                          timing_num_get_kap_calls            386434

                                           threads                12


                                             total    9.571    1.000

                                            matrix    2.568    0.268
                                               net    1.973    0.206
                                             solve    1.947    0.203
                                               eos    1.033    0.108
                                        hydro_vars    0.666    0.070
                                       neu_and_kap    0.460    0.048
                                               mlt    0.429    0.045
                                       evolve_step    0.240    0.025
                                            remesh    0.168    0.018
                                         run1_star    0.056    0.006
                                       adjust_mass    0.019    0.002
                                       mixing_info    0.012    0.001



Last-Updated: 03Jul2021 (MESA 094ff71) by fxt.
