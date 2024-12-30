include ../make/defaults-module.mk

# Build

MODULE_NAME := atm
SRCS := public/atm_def.f90 \
        public/atm_lib.f90 \
        private/atm_t_tau_uniform.f90 \
        private/atm_t_tau_varying.f90 \
        private/atm_t_tau_relations.f90 \
        private/atm_irradiated.f90 \
        private/atm_table.f90 \
        private/atm_utils.f90 \
        private/table_atm.f90
SRCS_CHECK = test/src/test_atm.f
INTERNAL_DEPENDS_ON := const math utils chem num auto_diff interp_1d interp_2d eos
EXTERNAL_DEPENDS_ON :=
BINTYPE := static-lib

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := atm_def.mod atm_lib.mod
INSTALL_INCLUDES :=
INSTALL_COMMANDS += install-atm-data

include $(MAKE_DIR)/Makefile

install-atm-data:
	mkdir -p ../data/atm_data && cp atm_data/* ../data/atm_data
