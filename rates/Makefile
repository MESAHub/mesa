include ../make/defaults-module.mk

# Build

MODULE_NAME := rates
SRCS := public/rates_def.f90 \
        public/rates_lib.f90 \
        private/rates_names.f90 \
        private/rates_reverses.f90 \
        private/weak_support.f90 \
        private/weaklib_tables.f90 \
        private/suzuki_tables.f90 \
        private/load_weak.f90 \
        private/eval_weak.f90 \
        private/coulomb.f90 \
        private/load_ecapture.f90 \
        private/eval_psi.f90 \
        private/eval_coulomb.f90 \
        private/eval_ecapture.f90 \
        private/reaclib_eval.f90 \
        private/reaclib_support.f90 \
        private/reaclib_print.f90 \
        private/reaclib_input.f90 \
        private/ratelib.f90 \
        private/pycno.f90 \
        private/raw_rates.f90 \
        private/rates_initialize.f90 \
        private/rates_support.f90 \
        private/screen5.f90 \
        private/screen_chugunov.f90 \
        private/screen.f90
SRCS_CHECK := test/src/test_rates.f90
INTERNAL_DEPENDS_ON := const utils math auto_diff interp_1d interp_2d chem eos forum
EXTERNAL_DEPENDS_ON :=
BINTYPE := static-lib

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := rates_def.mod rates_lib.mod
INSTALL_INCLUDES :=

include $(MAKE_DIR)/Makefile
