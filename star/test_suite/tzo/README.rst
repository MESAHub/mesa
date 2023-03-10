.. _tzo:

***
TZO
***

This test case builds and evolves a Thorne-Zytkow object (TZO) (Farmer et al in prep)

The model stats by building a 20 |Msun| star mid way through the main sequence (inlist_initial). The precise physics do not matter here,
but having a large convective burning core is necessary for the the TZO formation to proceed smoothly.

inlist_initial_make then moves the inner boundary to a desired inner mass and radius, assumed for the NS. It is helpful to have some
assumed energy injection at this point in time (core_avg_eps).

inlist_initial_post_relax then alters the total mass of the TZO and the composition to match the required starting values. While relax_Z does change the total
metal mass fraction, it does not change the metal distribution. Thus you may want to use relax_composition_filename instead to have complete control of the composition.

inlist_evolve_tzo then evolves the TZO for as long as possible. For the test suite the end condition if the NS sufficient mass via Eddington limited accretion.
With this check removed the end condition is when the surface begins undergoing large amplitude pulsations and we can no longer follow the evolution
(this usually occurs around v_surf/v_escape ~ 40%). Note that as the NS will accrete at Eddington limited rates (10^-8 Msun/yr) it grows very little over the typical
lifetime of a TZO, thus you should not use the test suite ending condition for actual science.
