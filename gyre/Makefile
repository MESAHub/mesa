include ../make/defaults-module.mk

# Build

MODULE_NAME := gyre

ifeq ($(WITH_GYRE),yes)
SRCS :=
SRCS_CHECK := test/src/test_gyre.f90
INTERNAL_DEPENDS_ON := forum const
EXTERNAL_DEPENDS_ON := hdf5_fortran lapack95 lapack blas
CHECK_RESULTS_GOLDEN := test/test_output
CHECK_DIFF_PROG := ndiff -quiet -relerr 1.0e-6
BINTYPE := lib
LIB_NAMES = gyre_mesa gyre

ifeq ($(WITH_CRLIBM),yes)
  EXTERNAL_DEPENDS_ON += crmath
endif

else

SRCS := public/gyre_mesa_m_stub.f90
SRCS_CHECK :=
INTERNAL_DEPENDS_ON := forum
EXTERNAL_DEPENDS_ON :=
CHECK_RESULTS_GOLDEN :=
BINTYPE := lib

endif

INCLUDE_DIRS :=
INSTALL_INCLUDES :=
MODULES := gyre_mesa_m.mod

include $(MAKE_DIR)/Makefile

ifeq ($(WITH_GYRE),yes)

ifeq ($(PROFILE),debug)
  GYRE_DEBUG = yes
else
  GYRE_DEBUG = false
endif

BUILD_DIR_GYRE_BUILD := $(BUILD_DIR_MODULE)/src

$(BUILD_DIR_GYRE_BUILD): gyre-HEAD.tar.gz | $(BUILD_DIR_GYRE_BUILD)/. $(BUILD_DIR_MODULE)/include/.
	tar --strip-components=1 -xf $< -C $(BUILD_DIR_GYRE_BUILD)

$(BUILD_DIR_GYRE_BUILD)/build-complete: $(BUILD_DIR_GYRE_BUILD)
	$(MAKE) -C $(BUILD_DIR_GYRE_BUILD) SHARED=$(DYNAMIC) DEBUG=$(GYRE_DEBUG) OMP=$(WITH_OPENMP) \
            CRMATH=$(WITH_CRLIBM) FRONTENDS=no TOOLS=no EXTERNAL_FORUM=yes EXTERNAL_FORUM_PKG=mesa-forum IFACES=yes
	touch $@


$(OBJ_OUT): $(BUILD_DIR_GYRE_BUILD)/build-complete | $(BUILD_DIR_MODULE)/lib/.
	cp $(BUILD_DIR_GYRE_BUILD)/build/$(notdir $@) $@

$(BUILD_DIR_MODULE)/modules/gyre_mesa_m.mod: $(BUILD_DIR_GYRE_BUILD)/build-complete | $(BUILD_DIR_MODULE)/modules/.
	cp $(BUILD_DIR_GYRE_BUILD)/build/$(notdir $@) $@

all: $(PKG_CONFIG) $(OBJ_OUT)

else

check: $(OBJ_OUT);

endif
