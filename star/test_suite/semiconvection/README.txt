This test is to show the capability of MESA to correctly place the convective and semiconvective boundaries using the Ledoux criterion and predictive mixing (see MESAIV paper), and allow mixing from semiconvection according to the Langer prescription.

It evolves a 1.5M model from zero age main sequence to when the central hydrogen mass fraction drops below 0.4. At that time the convective core has almost (but not quite) reached its maximum mass extension, and there is a semiconvective layer above the convective core.

To verify that the test ran successfully, MESA checks the mixing types at three points (0.12, 0.135, and 0.15Msun), where the convective types should be convection, semiconvection, and no mixing, and the average temperature and density between two points in the star. Target values and ranges are listed in src/run_star_extras.f. If the mixing types match and the temperature and density values fall within the given range, the terminal output at the end of the run should read ‘‘all values are within tolerances’’.

The inlist used to produce the starting zams model (inlist_to_ZAMS) is provided.

Note that the choice made for the inital mixture (a09 in this test case) and for the nuclear reaction network (here: pp_and_cno_extras.net) influence the growth of the convective core. The values used in this test case for the predictive_superad_thresh parameters work well here but need to be adjusted when using other mixtures or nuclear reaction networks.

Also note that for speed purposes the values used here for the maximum allowed timestep and the mesh size do not produce completely converged models. A convergence study should be done when using this inlist for science purposes.





 