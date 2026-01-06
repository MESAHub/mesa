# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MESA (Modules for Experiments in Stellar Astrophysics) is an open-source software suite for 1D stellar evolution calculations. Written primarily in Fortran 90/2008, it simulates stellar evolution from pre-main sequence through advanced stages including binary interactions, pulsations, and compact objects.

## Build Commands

```bash
# Full build (requires MESA SDK and MESA_DIR set)
export MESA_DIR=$(pwd)
./install

# Build single module (from module directory)
make

# Build with options
make PROFILE=debug              # debug|release|release-with-dbg-info (default)
make WITH_OPENMP=no            # Disable OpenMP
make WITH_PGSTAR=no            # Disable PGPLOT visualization

# Clean
./clean                        # Remove cached data
make clean                     # Remove build directory
```

## Running Tests

```bash
# Module unit tests (from module directory)
make check

# Test suite (from star/test_suite, binary/test_suite, or astero/test_suite)
./count_tests                  # Count available tests
./list_tests                   # List all tests
./list_tests 13                # Show test 13 name
./each_test_run                # Run all tests
./each_test_run 13             # Run single test by number

# Test output in out.txt and err.txt after run
```

## Linting

```bash
# Fortran
fortitude check

# Python
ruff format --check
ruff check

# RST documentation
sphinx-lint --ignore=gyre/gyre

# Custom MESA linters (from linters/ directory)
python check_columns.py
python check_defaults.py
python check_pgstar.py
python check_photos.py
```

## Architecture

### Module Build Order (dependencies flow left to right)
```
const → utils → math → mtx → auto_diff → num → interp_1d/interp_2d
                                                      ↓
chem, eos, forum, colors, rates, neu, net, kap, ionization, atm, turb
                                                      ↓
                              gyre ←   star_data  →  star
                                                      ↓
                                              astero, binary
```

### Module Structure
Each module directory contains:
- `Makefile` - Build configuration (MODULE_NAME, SRCS, BINTYPE, dependencies)
- `private/` - Internal implementation
- `public/` - Public API (module interfaces)
- `defaults/` - Default control parameters
- `test/` - Unit tests

### Key Directories
- `make/` - Build system (subdirs.mk defines build order)
- `data/` - Physics data caches (EOS, opacity, rates tables)
- `star/` - Main stellar evolution engine (~130k LOC)
- `star/test_suite/` - ~107 test cases
- `linters/` - Custom Python linting tools

## Fortran Code Style

- 3-space indentation (see `.editorconfig`)
- Use `real(dp)` not `double precision`
- Use `1.1d0` or `1.1_dp` for double literals, never bare `1.1`
- Use `pow(x, 1.5d0)` and `powN(x)` from `math_lib` instead of `**` operator
- For `auto_diff` types, must use `powN` even for squaring
- Return `ierr` (0=success) from subroutines for error handling
- Use `call mesa_error(__FILE__, __LINE__)` for critical errors, not `stop`
- Use explicit formats for `write` statements (different compilers use different defaults)
- Name OMP critical blocks uniquely to avoid serialization

## Test Suite Guidelines

- Tests must not write to stderr (causes automatic failure)
- Multi-part tests use `rn` script to run sequential `inlist_*_header` files
- Use `test_suite_helpers` for TestHub-compatible output
- Include `test_suite_extras` calls in `run_star_extras`
- Set `required_termination_code_string` for stopping condition validation

## CI Keywords in Commit Messages

- `[ci skip]` - Compile only, skip test suite
- `[ci optional]` - Run all test parts (no skipping)
- `[ci fpe]` - Enable floating-point exception checks
- `[ci split]` - Split test suite between machines
