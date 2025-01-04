include ../make/defaults-module.mk

# Build

MODULE_NAME := kap
SRCS = \
   private/cubic_interpolator.f90 \
   private/op_def.f90 \
   private/op_load.f \
   private/op_load_master.f \
   private/op_common.f \
   private/op_ev.f \
   private/op_osc.f \
   private/op_radacc.f \
   private/op_eval.f \
   private/op_eval_mombarg.f90 \
   private/kap_eval_support.f90 \
   private/condint.f90 \
   private/load_co_kap.f90 \
   private/load_kap.f90 \
   private/kap_eval_fixed.f90 \
   private/kap_eval_co.f90 \
   private/kap_eval.f90 \
   private/kap_ctrls_io.f90 \
   private/kapcn.f90 \
   private/kap_aesopus.f90 \
   public/kap_def.f90 \
   public/kap_lib.f90 \
   other/other_elect_cond_opacity.f90 \
   other/other_compton_opacity.f90 \
   other/other_radiative_opacity.f90
SRCS_CHECK = \
   test/src/test_kap.f90 \
   test/src/test_kap_support.f90
INTERNAL_DEPENDS_ON := const utils math chem interp_2d forum eos
EXTERNAL_DEPENDS_ON :=
BINTYPE := static-lib
INCLUDE_DIRS := -Iprivate -Ipublic -Idefaults

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := kap_def.mod kap_lib.mod
INSTALL_COMMANDS += install-kap

include $(MAKE_DIR)/Makefile

install-kap:
	if [ ! -r ../data/kap_data ]; then \
                mkdir -p ../data/kap_data; \
		tar -xJf kap_data.tar.xz -C ../data; \
		tar -xJf kapcn_data.txz -C ../data; \
		cp AESOPUS_AGSS09.h5 ../data/kap_data; \
	fi
