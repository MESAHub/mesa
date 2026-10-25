.. TDC_Cepheid:

***********
TDC_Cepheid
***********

This example loads an evolutionary Cepheid model, removes its core, remeshes
the envelope, and uses GYRE to excite the fundamental radial mode. Positive
kinetic energy growth per cycle indicates growing pulsations. The run can be
continued to finite amplitude, as in RSP. The metallicity matches the LMC
Cepheid models in Labs 1 and 2.

Set ``load_model_filename`` in ``inlist_pulses`` to a saved evolutionary
Cepheid model. This directory does not include a starting model and can be
used for Cepheid models spanning the Hertzsprung progression.


Last-Updated: 2026-05-14 (mesa r26.04.1) by Ebraheem Farag
