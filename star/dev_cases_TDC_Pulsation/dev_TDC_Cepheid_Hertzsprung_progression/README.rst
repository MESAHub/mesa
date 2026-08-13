.. TDC_Cepheid:

***********
TDC_Cepheid
***********

This test case loads an evolutionary Cepheid model from MESA-star. The model then has the core removed, is remeshed and delivered a kick in the fundamental radial mode using GYRE in MESA. The model should display a
growing kinetic energy given by a positive kinetic energy growth rate per cycle, until finite amplitude pulsations are achieved, similar to RSP. The
metallicity is chosen such that this test case can be directly compared to the Lab 1 and Lab 2 LMC Cepheid models.

Set ``load_model_filename`` in ``inlist_pulses`` to a saved evolutionary
Cepheid model before running this example. Unlike the mass-specific dev
cases, this directory is a template for models spanning the Hertzsprung
progression and does not include a starting model.


Last-Updated: 2026-05-14 (mesa r26.04.1) by Ebraheem Farag
