include ../make/defaults-module.mk

# Build

MODULE_NAME := eos
SRCS := public/eos_def.f90 \
        public/eos_lib.f90 \
        private/create_EXCOR7_table.f90 \
        private/create_FSCRliq8_table.f90 \
        private/eos_blend.f90 \
        private/eoscms_eval.f90 \
        private/eos_ctrls_io.f90 \
        private/eosde_eval.f90 \
        private/eosdt_eval.f90 \
        private/eosdt_load_tables.f90 \
        private/eosdt_support.f90 \
        private/eos_helm_eval.f90 \
        private/eos_initialize.f90 \
        private/eospc_eval.f90 \
        private/eospt_eval.f90 \
        private/gauss_fermi.f90 \
        private/helm_alloc.f90 \
        private/helm.f90 \
        private/helm_polynomials.f90 \
        private/ideal.f90 \
        private/ion_offset.f90 \
        private/pc_eos.f90 \
        private/pc_support.f90 \
        private/skye_coulomb.f90 \
        private/skye_coulomb_liquid.f90 \
        private/skye_coulomb_solid.f90 \
        private/skye.f90 \
        private/skye_ideal.f90 \
        private/skye_thermodynamics.f90 \
        other/other_eos.f90
SRCS_CHECK = test/src/test_eos.f90 \
             test/src/eos_support.f90 \
             test/src/test_eos_blend.f90 \
             test/src/test_eos_support.f90
INTERNAL_DEPENDS_ON := const math utils chem num auto_diff interp_1d interp_2d
EXTERNAL_DEPENDS_ON :=
BINTYPE := lib
INCLUDE_DIRS := -Iprivate -Ipublic -Idefaults

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := eos_def.mod eos_lib.mod
INSTALL_INCLUDES := public/helm_def.dek
EXTRACT_DATA_COMMANDS += ../data/eosFreeEOS_data/.extract-completed ../data/eosCMS_data/.extract-completed ../data/eosDT_data/.extract-completed ../data/eosPC_support_data/.extract-completed

include $(MAKE_DIR)/Makefile

../data/eosFreeEOS_data/.extract-completed: eosFreeEOS_data.tar.xz
	rm -rf ../data/eosFreeEOS_data/
	rm -rf ../data/eosDT_data/cache/
	tar -xJf eosFreeEOS_data.tar.xz -C ../data
	touch $@

../data/eosCMS_data/.extract-completed: eosCMS_data.tar.xz
	rm -rf ../data/eosCMS_data/
	tar -xJf eosCMS_data.tar.xz -C ../data
	touch $@

../data/eosDT_data/.extract-completed: eosDT_data.tar.xz helm_table.dat.xz
	rm -rf ../data/eosDT_data/
	tar -xJf eosDT_data.tar.xz -C ../data
	xz --decompress --keep --stdout helm_table.dat.xz > ../data/eosDT_data/helm_table.dat
	touch $@

../data/eosPC_support_data/.extract-completed:
	rm -rf eosPC_support_data
	mkdir -p ../data/eosPC_support_data
	mkdir -p ../data/eosPC_support_data/cache
	touch $@
