include ../make/defaults-module.mk

# Build

MODULE_NAME := star_data
SRCS := public/star_data_def.f90 \
	public/star_data_lib.f90 \
	public/star_pgstar.f90

SRCS_CHECK = test/src/test_star_data.f
INTERNAL_DEPENDS_ON := const chem net kap colors neu
EXTERNAL_DEPENDS_ON :=
BINTYPE := static-lib

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := star_data_lib.mod star_data_def.mod star_pgstar.mod
INSTALL_INCLUDES := private/pgstar_controls.inc \
                    private/star_controls_dev.inc \
                    private/star_controls.inc \
                    private/star_job_controls_dev.inc \
                    private/star_job_controls.inc \
                    private/star_job_controls_params.inc

include $(MAKE_DIR)/Makefile
