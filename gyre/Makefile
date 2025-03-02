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
BINTYPE := manual
LIB_NAMES = gyre_mesa gyre
INSTALL_COMMANDS := install-lib

else

SRCS := public/gyre_mesa_m_stub.f90
SRCS_CHECK :=
INTERNAL_DEPENDS_ON := forum
EXTERNAL_DEPENDS_ON :=
CHECK_RESULTS_GOLDEN :=
BINTYPE := static-lib

endif

INCLUDE_DIRS :=
INSTALL_INCLUDES :=
MODULES := gyre_mesa_m.mod

include $(MAKE_DIR)/Makefile

ifeq ($(WITH_GYRE),yes)

include $(MAKE_DIR)/pkg-config.mk

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
OBJ_OUT := $(BUILD_DIR_MODULE)/lib/libgyre.a $(BUILD_DIR_MODULE)/lib/libgyre_mesa.a

$(BUILD_DIR_GYRE_BUILD): gyre-7.2.1.tar.gz | $(BUILD_DIR_GYRE_BUILD)/ $(BUILD_DIR_MODULE)/include/
	tar --strip-components=1 -xf $< -C $(BUILD_DIR_GYRE_BUILD)

$(BUILD_DIR_GYRE_BUILD)/build-complete: $(BUILD_DIR_GYRE_BUILD)
	$(MAKE) -C $(BUILD_DIR_GYRE_BUILD) SHARED=no DEBUG=$(GYRE_DEBUG) OMP=$(WITH_OPENMP) CRMATH=$(WITH_CRLIBM) FPE=yes VERSION_CHECK=passed FPP="$(FPP)" LDFLAGS="$(LIB_DEP_ARGS)" FFLAGS="$(FLAGS_DEPS)" FRONTENDS=no TOOLS=no FORUM=no IFACES=yes
	touch $@


$(BUILD_DIR_MODULE)/lib/libgyre.a $(BUILD_DIR_MODULE)/lib/libgyre_mesa.a: $(BUILD_DIR_GYRE_BUILD)/build-complete | $(BUILD_DIR_MODULE)/lib/
	cp $(BUILD_DIR_GYRE_BUILD)/build/$(notdir $@) $@

$(BUILD_DIR_MODULE)/modules/gyre_mesa_m.mod: $(BUILD_DIR_GYRE_BUILD)/build-complete | $(BUILD_DIR_MODULE)/modules/
	cp $(BUILD_DIR_GYRE_BUILD)/build/$(notdir $@) $@

all: $(PKG_CONFIG) $(BUILD_DIR_MODULE)/lib/libgyre_mesa.a $(BUILD_DIR_MODULE)/lib/libgyre.a

else

check: $(OBJ_OUT);

endif
