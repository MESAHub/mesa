include ../make/defaults-module.mk

# Build

MODULE_NAME := forum
SRCS :=
SRCS_CHECK :=
INTERNAL_DEPENDS_ON :=
EXTERNAL_DEPENDS_ON := hdf5_fortran
BINTYPE := static-lib
INCLUDE_DIRS :=

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := forum_m.mod
INSTALL_INCLUDES :=

include $(MAKE_DIR)/Makefile

ifeq ($(PROFILE),debug)
  FORUM_DEBUG = yes
else
  FORUM_DEBUG = false
endif

ifeq ($(MESASDK_ROOT),)
  FPP=./fypp_deps
else
  FPP= PYTHONPATH=$(MESASDK_ROOT)/lib/python ./fypp_deps
endif

BUILD_DIR_FORUM_BUILD := $(BUILD_DIR_MODULE)/src
INSTALL_INCLUDES_BUILD += $(BUILD_DIR_MODULE)/include/forum.inc

$(BUILD_DIR_FORUM_BUILD): forum-1.0.3.tar.gz | $(BUILD_DIR_FORUM_BUILD)/ $(BUILD_DIR_MODULE)/include/
	tar --strip-components=1 -xf $< -C $(BUILD_DIR_FORUM_BUILD)
	cp $(BUILD_DIR_FORUM_BUILD)/src/include/forum.inc $(BUILD_DIR_MODULE)/include/forum.inc

$(BUILD_DIR_FORUM_BUILD)/build-complete: $(BUILD_DIR_FORUM_BUILD)
	$(MAKE) -C $(BUILD_DIR_FORUM_BUILD) SHARED=no DEBUG=$(FORUM_DEBUG) OMP=$(WITH_OPENMP) FPE=yes VERSION_CHECK=passed FPP="$(FPP)" LDFLAGS="$(LIB_DEP_ARGS)" FFLAGS="$(FLAGS_DEPS)"
	touch $@

$(OBJ_OUT): $(BUILD_DIR_FORUM_BUILD)/build-complete | $(BUILD_DIR_MODULE)/lib/
	cp $(BUILD_DIR_FORUM_BUILD)/build/$(notdir $@) $@

$(BUILD_DIR_MODULE)/modules/forum_m.mod: $(BUILD_DIR_FORUM_BUILD)/build-complete | $(BUILD_DIR_MODULE)/modules/
	cp $(BUILD_DIR_FORUM_BUILD)/build/$(notdir $@) $@

check:
	pushd test > /dev/null ; ../$(BUILD_DIR_FORUM_BUILD)/build/test_order > ../$(CHECK_RESULTS); popd > /dev/null
	diff -b $(CHECK_RESULTS) $(CHECK_RESULTS_GOLDEN)
