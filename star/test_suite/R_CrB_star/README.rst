.. _R_CrB_star:

**********
R_CrB_star
**********

This test case creates and evolves a simple model of an R Corona
Borealis star.  The model is constructed from an 0.875 |Msun|
homogeneous He star and evolved until it first reaches a low effective
temperature.  The calculation is based on the |Schwab2019|
reproduction of the work of |Weiss1987|.

This test suite provides an example of how to use AESOPUS opacity
tables in MESA.  See $MESA_DIR/kap/preprocessor/AESOPUS/README for
more information.

The (xz compressed) directory RCrB_GS98 contains the opacity files
obtained from the webform.  The preprocessing steps

.. code-block:: sh
  
    tar xf RCrB_GS98.tar.xz
    $MESA_DIR/kap/preprocessor/AESOPUS/aesopus.py RCrB_GS98.yaml

generated the AESOPUS_GS98_RCrB.h5 file.  This is the opacity file
that is read by MESA.  The .dat files in the RCrB_GS98 tarball are
included only to demonstrate the workflow.

The commands that configure the opacity tables are in
``inlist_common``.

This test case has 3 parts, but by default, saved models are used to
skip to part 3.  Options shared between all parts are contained in
``inlist_common``.

* Part 1 (``inlist_He_star``) creates a hydrogen-free initial model
  and evolves it from the pre-main-sequence through He core burning
  until there is less than 0.45 |Msun| of He remaining in the
  envelope.

* Part 2 (``inlist_change_abundances``) manually alters the envelope
  composition so that it has a carbon-enhanced composition,
  corresponding to mixture R2 in |Weiss1987|.  It then evolves the
  model until the envelope begins to expand and the effective
  temperature falls to 20,000 K.

* Part 3 (``inlist_R_CrB_star``) evolves the models from effective
  temperature 20,000 K to 6,000 K so that the outer layers cross the
  opacity blend from the OPAL opacities to the AESOPUS opacities.


.. |Weiss1987| replace:: `Weiss (1987) <https://ui.adsabs.harvard.edu/abs/1987A%26A...185..165W/abstract>`__           
.. |Schwab2019| replace:: `Schwab (2019) <https://ui.adsabs.harvard.edu/abs/2019ApJ...885...27S/abstract>`__


Last-Updated: 2020-11-16 (mesa r14909) by Josiah Schwab

