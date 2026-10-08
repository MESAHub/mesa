dev_test_star_lna
================

Load one saved model, run radial ``star_LNA``, export ``gyre.data``, and
exit. The driver has no evolution loop. It does not kick, relax or remesh
the model. No batch runner or external model catalog is required.

From this directory, after installing this MESA checkout::

   export MESA_DIR="$(cd ../../.. && pwd)"
   make
   ./rn
   python plotter/plot_modes.py

``rn`` runs the existing executable without compiling it. The plotter
requires NumPy and Matplotlib and reads only the saved output. It writes
``plots/spectrum.png`` and ``plots/mode_0.png``. Select another zero-based
terminal-table index with ``python plotter/plot_modes.py --mode 3``.
The corresponding MESA eigenfunction filename is numbered one higher.

Models
------

``inlist_model`` selects the tracked RGB-tip model from
``../../test_suite/1M_pre_ms_to_wd/standard_start_he_core_flash.mod``.
The model has about 0.804 solar masses. Its source case supplies the
GS98 opacity family, Zbase and mixing-length setting used here.

To select the 5 solar mass example::

   cp inlist_5M inlist_model
   ./rn
   python plotter/plot_modes.py

This loads ``5M_cepheid_blue_loop/standard_start.mod``, the red starting
model before its Cepheid crossing. It is not a saved model inside the
instability strip. A later saved model from that case can be selected by
changing ``load_model_filename``. Match the opacity, composition and
mixing-length controls to any replacement model.

The common analysis controls are in ``inlist_lna``. The default
``star_LNA_T_inner = 1d7`` restricts the eigenproblem to the envelope.
Set it to zero to analyze the full model. This does not truncate the
loaded structure or the GYRE export. The supplied atmosphere and TDC
controls define the analysis background; loading a model does not
reconverge it under changed controls.
``include_mlt_corr_to_TDC = .true.`` retains the MLT efficiency correction
for the supplied RGB model.

Output
------

``LNA/`` contains periods, growth rates, eigenvector residuals,
eigenfunctions, diagnostic work terms and matrix/background diagnostics.
The spectrum plot shows the selected modes. A selection index is not
radial order, and mode zero is not automatically the fundamental.
Positive ``logKE_per_cycle`` denotes growth. Work curves are diagnostics;
their sum is not required to reproduce the eigenvalue growth rate.

The default residual limit is ``1d-6``. Inspect rejected-mode messages
and background residuals if no modes pass. This is a development example,
not a calibrated period or growth-rate regression. Loading the supplied
RGB model does not guarantee a steady TDC background.
Check ``star_LNA_tdc_face_audit.data`` before interpreting growth rates.

``gyre.data`` uses experimental schema 130. A compatible GYRE reader is
required. For this single snapshot the driver sets ``dt = 0`` only during
export, then restores it. The exporter writes zero background
``eps_grav`` without changing the stored heating or the LNA thermal
perturbation. A ``.mod`` file does not store ``eps_grav`` or restore the
previous thermal state needed for its finite-step calculation. Its saved
timestep is insufficient, also with ``energy_eqn_option = 'dedt'``.
Zero is a static assumption, not an estimate of evolutionary heating.
An estimate from the luminosity balance would require accounting for the
energy sources and other storage and work terms; this driver does not
make that estimate.

``add_atmosphere_to_pulse_data = .false.`` keeps the GYRE export at the
MESA surface. Setting it true appends a ``T(tau)`` atmosphere, by default
to optical depth ``1d-3``. It does not change the star LNA domain or its
atmospheric boundary condition. The appended points have zero convective
flux and velocity; the extra schema-130 thermodynamic fields, including
``Cp_face`` and ``Hp_face``, are also left at zero. They are not an
extension of the TDC solution.

Each run uses the same output filenames. Archive ``LNA``, ``gyre.data``
and ``plots`` before changing the model if both results are needed.
