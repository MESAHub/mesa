include ../make/defaults-module.mk

# Build

MODULE_NAME := auto_diff
SRCS := public/auto_diff.f90 \
        private/auto_diff_real_15var_order1_module.f90 \
        private/auto_diff_real_1var_order1_module.f90 \
        private/auto_diff_real_2var_order1_module.f90 \
        private/auto_diff_real_2var_order2_module.f90 \
        private/auto_diff_real_2var_order3_module.f90 \
        private/auto_diff_real_4var_order1_module.f90 \
        private/auto_diff_real_star_order1_module.f90 \
        private/auto_diff_real_tdc_module.f90 \
        private/support_functions.f90
SRCS_CHECK := test/src/test_auto_diff.f90
INTERNAL_DEPENDS_ON := const utils math
EXTERNAL_DEPENDS_ON :=
BINTYPE := static-lib

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := auto_diff.mod
INSTALL_INCLUDES :=

include $(MAKE_DIR)/Makefile
