include ../make/defaults-module.mk

# Build

MODULE_NAME := net
SRCS := public/net_def.f90 \
        public/net_lib.f90 \
        private/net_screen.f90 \
        private/net_derivs_support.f90 \
        private/net_derivs.f90 \
        private/net_approx21.f90 \
        private/net_burn_support.f90 \
        private/net_initialize.f90 \
        private/net_eval.f90 \
        private/net_burn.f90 \
        private/net_burn_const_density.f90 \
        private/net_burn_const_p.f90
SRCS_CHECK := test/src/test_net.f90
INTERNAL_DEPENDS_ON := const utils math mtx num chem rates interp_1d eos
EXTERNAL_DEPENDS_ON :=
BINTYPE := static-lib

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := net_def.mod net_lib.mod
INSTALL_INCLUDES := public/burner_const_density_get_eos_info.inc \
                    public/burner_finish_substep.inc

include $(MAKE_DIR)/Makefile
