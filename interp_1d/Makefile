include ../make/defaults-module.mk

# Build

MODULE_NAME := interp_1d
SRCS := public/interp_1d_def.f90 \
        public/interp_1d_lib_sg.f90 \
        public/interp_1d_lib.f90 \
        private/interp_1d_misc_sg.f90 \
        private/interp_1d_mp.f90 \
        private/interp_1d_mp_sg.f90 \
        private/interp_1d_pm.f90 \
        private/interp_1d_pm_sg.f90 \
        private/interp_1d_misc.f90 \
        private/interp_1d_mp_autodiff.f90 \
        private/interp_1d_pm_autodiff.f90
SRCS_CHECK = test/src/interp_1d_support.f \
	     test/src/interp_1d_support_sg.f \
	     test/src/test_interp_1d.f
INTERNAL_DEPENDS_ON := const math utils auto_diff
EXTERNAL_DEPENDS_ON :=
BINTYPE := static-lib

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := interp_1d_lib.mod interp_1d_def.mod interp_1d_lib_sg.mod
INSTALL_INCLUDES :=

include $(MAKE_DIR)/Makefile
