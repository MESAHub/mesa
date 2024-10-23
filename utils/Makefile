include ../make/defaults-module.mk

# Build

MODULE_NAME := utils
SRCS := public/utils_def.f90 \
	public/utils_lib.f90 \
	private/utils_c_system.c \
	private/utils_dict.f90 \
	private/utils_idict.f90 \
	private/utils_nan_dp.f90 \
	private/utils_nan.f90 \
	private/utils_nan_qp.f90 \
	private/utils_nan_sp.f90 \
	private/utils_system.f90

ifeq ($(WITH_OPENMP),yes)
  SRCS += private/utils_openmp.f90
else
  SRCS += private/utils_openmp_stub.f90
endif

SRCS_CHECK = test/src/test_utils.f90
INTERNAL_DEPENDS_ON := const
EXTERNAL_DEPENDS_ON :=
BINTYPE := static-lib

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := utils_def.mod utils_lib.mod
INSTALL_INCLUDES := formats

include $(MAKE_DIR)/Makefile
