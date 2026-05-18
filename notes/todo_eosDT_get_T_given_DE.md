# TODO: public density-energy EOS T wrapper

Status: note only. The source edits described here were reverted after capture.

## Goal

Add a public EOS-library convenience wrapper equivalent to
`star/private/eos_support.f90::solve_eos_given_DE`, but without requiring a
`star_info` pointer.

Candidate public routine:

```fortran
subroutine eosDT_get_T_given_DE( &
      handle, species, chem_id, net_iso, xa, &
      logRho, logE, logT_guess, logT_tol, logE_tol, &
      logT_result, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
      d_dxa_const_TRho, eos_calls, ierr)
```

## Intended implementation

Place the wrapper after `eosDT_get_T` in `eos/public/eos_lib.f90`.

The wrapper should call `eosDT_get_T` with `which_other = i_lnE`:

```fortran
call eosDT_get_T( &
      handle, &
      species, chem_id, net_iso, xa, &
      logRho, i_lnE, logE*ln10, &
      logT_tol, logE_tol*ln10, max_iter_for_solve, logT_guess, &
      arg_not_provided, arg_not_provided, arg_not_provided, arg_not_provided, &
      logT_result, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
      d_dxa_const_TRho, eos_calls, ierr)
```

This mirrors `solve_eos_given_DE` in `star/private/eos_support.f90`.

## Unit conversion check

`eosDT_get_T` expects EOS result variables in their stored form. For energy,
`res(i_lnE)` stores

```text
res(i_lnE) = ln(E)
```

The star wrapper input is base-10:

```text
logE = log10(E)
```

Therefore the value passed to `eosDT_get_T` should be

```text
ln(E) = log10(E) * ln(10)
```

and the tolerance should be converted consistently:

```text
other_tol = logE_tol * ln(10)
```

The root solve is therefore

```text
f(logT) = res(i_lnE; logRho, logT, xa) - logE*ln(10) = 0
```

## Code references

- `eos/public/eos_lib.f90`: `eosDT_get_T` is the public lower-level root finder.
- `star/private/eos_support.f90`: `solve_eos_given_DE` is the existing star-side wrapper.
- `eos/public/eos_def.f90`: `i_lnE` indexes the natural-log EOS energy result.

## TODO

- [ ] Add `eosDT_get_T_given_DE` to `eos/public/eos_lib.f90`.
- [ ] Use `i_lnE`, `logE*ln10`, and `logE_tol*ln10`.
- [ ] Use `arg_not_provided` for unknown bounds, matching `solve_eos_given_DE`.
- [ ] Keep `max_iter_for_solve = 100` unless the API should expose `max_iter`.
- [ ] Optionally rewire `solve_eos_given_DE` to call the public wrapper.
- [ ] Compile or run tests only after explicit permission.
