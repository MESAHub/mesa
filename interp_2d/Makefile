include ../make/defaults-module.mk

# Build

MODULE_NAME := interp_2d
SRCS := public/interp_2d_lib_db.f90 \
        public/interp_2d_lib_sg.f90 \
        private/akima760_db.f \
        private/akima760_sg.f \
        private/bicub_db.f \
        private/bicub_sg.f \
        private/bipm_db.f \
        private/bipm_sg.f \
        private/renka790_db.f \
        private/renka790_sg.f
SRCS_CHECK = test/src/interp_2d_support.f \
             test/src/test_akima_db.f \
             test/src/test_akima_sg.f \
             test/src/test_interp_2d.f \
             test/src/test_renka790_db.f \
             test/src/test_renka790_sg.f
INTERNAL_DEPENDS_ON := const utils math num interp_1d
EXTERNAL_DEPENDS_ON :=
BINTYPE := static-lib

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := interp_2d_lib_db.mod interp_2d_lib_sg.mod
INSTALL_INCLUDES :=

include $(MAKE_DIR)/Makefile
