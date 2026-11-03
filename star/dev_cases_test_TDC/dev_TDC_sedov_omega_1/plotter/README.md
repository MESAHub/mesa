# Sedov comparison plot

`plot_sedov_comparison.py` compares saved MESA profiles with the public-domain
Kamm-Timmes Sedov solver distributed by Cococubed. The analytic modules are a
Python port of the solver described by Kamm & Timmes (2007), *On Efficient
Generation of Numerically Robust Sedov Solutions*.

Source: <https://cococubed.com/research_pages/sedov.shtml>

From this directory, run:

```console
conda run -n base python plot_sedov_comparison.py
```

The defaults compare `LOGS_high_res` and `LOGS_sedov` near `t = 0.36 s`.

To plot one named run at a selected time and show its AMR thresholds, use:

```console
conda run -n base python plot_sedov_comparison.py \
   --logs LOGS_sedov --labels test --time 0.36 \
   --max-long 1.25 --max-short 2.5 \
   --output sedov_comparison_test_t036.pdf
```
