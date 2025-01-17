# Makefile sanity
MAKEFLAGS += --no-builtin-rules
MAKEFLAGS += --no-builtin-variables

.DELETE_ON_ERROR:
.SHELLFLAGS := -eu -c
.DEFAULT_GOAL := all
.PHONY: clean all

NOCOMPILE := $(filter clean,$(MAKECMDGOALS))

include $(MAKE_DIR)/subdirs.mk
include $(MAKE_DIR)/helpers.mk

include $(MAKE_DIR)/setup-builddir.mk

BUILD_DIR_MODULE := $(call clean-comp,$(shell echo $(dir $(firstword $(MAKEFILE_LIST))) | sed "s|/$$||")/$(BUILD_DIR))

$(BUILD_DIR_MODULE)/:
	mkdir -p $@

$(BUILD_DIR_MODULE)%/:
	mkdir -p $@

INTERNAL_DEPENDS_ON += $(SUBDIRS)
INCLUDE_DIRS ?=

ifeq ($(NOCOMPILE),)
  # Only the dependencies need to come from the original build dir
  include $(MAKE_DIR)/setup-depends.mk

  ifeq ($(MODULE_NAME),)
    $(error Missing module name)
  endif

  ifeq ($(COMPILER),gfortran)
    include $(MAKE_DIR)/compile-settings-gnu.mk
  else
    $(error Unknown or unset COMPILER)
  endif

  include $(MAKE_DIR)/compile.mk
  include $(MAKE_DIR)/link-exec.mk

  all: $(OBJ_OUT)
endif

clean:
	@rm -rf $(BUILD_DIR_MODULE)
