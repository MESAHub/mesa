# Turbulent energy during mass adjustment

This regression calls the private `adjust_mass` implementation in the current
source build. It checks identity, zero turbulence, cold accretion, wind energy
export, the limit of vanishing mass change, repeated calls from an unchanged
checkpoint, and the face shared by changed and unchanged cells. Independent cell
averaging checks the total turbulent-energy budget. Additional intersection
checks cover incoming and outgoing mass, including a single cell.

After building the star library with the same SDK and build settings, run:

```sh
export MESA_DIR=/path/to/mesa
make -C "$MESA_DIR/star/test/tdc_mass_remap" run
```

The four-cell fixtures keep the innermost cell mass fixed. Incoming turbulent
energy is zero. The single-cell case checks the helper's algebraic boundary
behavior. These checks exercise the native remap rather than a stellar solve;
stellar remeshing, retries and photo restart require integration tests.
