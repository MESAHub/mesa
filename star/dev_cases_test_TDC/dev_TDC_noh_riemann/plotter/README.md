# Noh energy-balance plotter

Run from the Noh work directory after `noh_energy_balance.txt` has been
created:

```sh
conda run -n base python plotter/plot_energy_balance.py
```

The parser streams the diagnostic file, ignores a final model that is still
being written, and keeps complete zone data only for the requested profile
models. Output is written to `plotter/plots`:

- `energy_balance_history.png` compares the integrated energy terms and their
  closure.
- `wall_heating_history.png` compares the innermost postshock state with the
  analytic Noh solution and locates the density deficit and shock-local
  internal-energy update relative to the analytic shock radius.
- `wall_heating_profiles.png` follows density, pressure, specific energy, and
  the entropy proxy `P/rho^gamma` through the wall region and shock.
- `energy_balance_model_*.png` shows the cell terms near the origin for the
  selected models.
- `energy_balance_summary.csv` records the integrated and maximum local
  diagnostics for every accepted model.

Select different snapshots or a fixed radial range with, for example,

```sh
conda run -n base python plotter/plot_energy_balance.py \
   --models 100 400 800 --radius-max 0.05
```
