include ../make/defaults-module.mk

# Build

MODULE_NAME := gyre

ifeq ($(WITH_GYRE),yes)
SRCS :=
SRCS_CHECK := test/src/test_gyre.f90
INTERNAL_DEPENDS_ON := forum const
EXTERNAL_DEPENDS_ON := hdf5_fortran lapack95 lapack crlibm-fortran
CHECK_RESULTS_GOLDEN := test/test_output
CHECK_DIFF_PROG := ndiff -quiet -relerr 1.0e-6
else
SRCS := public/gyre_mesa_m_stub.f90
SRCS_CHECK :=
INTERNAL_DEPENDS_ON := forum
EXTERNAL_DEPENDS_ON :=
CHECK_RESULTS_GOLDEN :=
endif
BINTYPE := static-lib
INCLUDE_DIRS :=

# Testing


# Install

MODULES := gyre_mesa_m.mod
INSTALL_INCLUDES :=

include $(MAKE_DIR)/Makefile

ifeq ($(WITH_GYRE),yes)

ifeq ($(PROFILE),debug)
  GYRE_DEBUG = yes
else
  GYRE_DEBUG = false
endif

ifeq ($(MESASDK_ROOT),)
  FPP=./fypp_deps $(FLAGS_DEPS)
else
  FPP= PYTHONPATH=$(MESASDK_ROOT)/lib/python ./fypp_deps $(FLAGS_DEPS)
endif

BUILD_DIR_GYRE_BUILD := $(BUILD_DIR_MODULE)/src

$(BUILD_DIR_GYRE_BUILD): gyre-7.2.1.tar.gz | $(BUILD_DIR_GYRE_BUILD)/ $(BUILD_DIR_MODULE)/include/
	tar --strip-components=1 -xf $< -C $(BUILD_DIR_GYRE_BUILD)

wrap-gyre: $(BUILD_DIR_GYRE_BUILD)
	$(MAKE) -C $(BUILD_DIR_GYRE_BUILD) SHARED=no DEBUG=$(GYRE_DEBUG) OMP=$(WITH_OPENMP) CRMATH=$(WITH_CRLIBM) FPE=yes VERSION_CHECK=passed FPP="$(FPP)" LDFLAGS="$(LIB_DEP_ARGS)" FFLAGS="$(FLAGS_DEPS)" FRONTENDS=no TOOLS=no FORUM=no IFACES=yes


$(OBJ_OUT): wrap-gyre | $(BUILD_DIR_MODULE)/lib/
	echo -e "create $@\naddlib $(BUILD_DIR_GYRE_BUILD)/build/libgyre.a\naddlib $(BUILD_DIR_GYRE_BUILD)/build/libgyre_mesa.a\nsave\nend" | ar -M

$(BUILD_DIR_MODULE)/modules/gyre_mesa_m.mod: wrap-gyre | $(BUILD_DIR_MODULE)/modules/
	cp $(BUILD_DIR_GYRE_BUILD)/build/$(notdir $@) $@

.PHONY: wrap-gyre

else

check: $(OBJ_OUT);

endif
