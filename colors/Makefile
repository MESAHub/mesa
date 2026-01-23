include ../make/defaults-module.mk

# Build

MODULE_NAME := colors
SRCS := public/colors_def.f90 \
        public/colors_lib.f90 \
        private/bolometric.f90 \
        private/colors_iteration.f90 \
        private/colors_ctrls_io.f90 \
        private/colors_history.f90 \
        private/colors_utils.f90 \
        private/hermite_interp.f90 \
        private/knn_interp.f90 \
        private/linear_interp.f90 \
        private/synthetic.f90
SRCS_CHECK := test/src/test_colors.f90
INTERNAL_DEPENDS_ON := const utils math
EXTERNAL_DEPENDS_ON :=
INCLUDE_DIRS := -Idefaults -Ipublic -Iprivate
BINTYPE := static-lib

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := colors_def.mod colors_lib.mod
INSTALL_INCLUDES :=
EXTRACT_DATA_COMMANDS += install-colors-data

include $(MAKE_DIR)/Makefile

install-colors-data:
	./build_data_and_export
