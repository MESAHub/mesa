.. _starspots:

******************
starspots
******************

This test suite case evolves a 1 solar mass, solar metallicity (Z=0.02) model to an age of 0.8 Gyr using MESA's implementation of the YREC (Yale Rotating Evolution Code) SPOTS formalism introduced by Somers et al. (2020; ApJ).

The impedence of the surface flux due to magnetic pressure from starspots is parameterized in the style of an atmospheric boundary modification.
As first described by Somers et al. (2015; ApJ), the degree of "spottiness" on the stellar surface is characterized using two parameters:
 * SPOTF (hereafter fspot): a coverage fraction, or "spot filling factor" (in the notation of the YREC documentation); and

* SPOTX (hereafter xspot): the temperature contrast between the spotted and unspotted regions: xspot = T_spot/T_photosphere.

Currently, x_ctrl(1) corresponds to fspot, and x_ctrl(2) corresponds to xspot. The coverage fraction is set to fspot = 0.34 (for consistency with observations of low-mass stars: Cao et al., 2022) and the temperature contrast is set to xspot = 0.85 (also from fits to observations).
 
Detailed discussion of this functionality can be found in MESA Instrument Paper VI: Starspots

Last-Updated: 2022-07-28 by Meridith Joyce

