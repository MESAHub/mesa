# Don't forget to update the documentation when you change something here!
# Makefile sanity
MAKEFLAGS += --no-builtin-rules
MAKEFLAGS += --no-builtin-variables

.DELETE_ON_ERROR:
.SHELLFLAGS := -eu -c
.DEFAULT_GOAL := all
.PHONY: clean print-mesa-deps install all check extract-data

NOCOMPILE := $(or \
  $(filter clean,$(MAKECMDGOALS)), \
  $(filter print-mesa-deps,$(MAKECMDGOALS)), \
  $(filter extract-data,$(MAKECMDGOALS)), \
)

include $(MAKE_DIR)/helpers.mk
include $(MAKE_DIR)/setup-builddir.mk
BUILD_DIR_MODULE := $(BUILD_DIR_)/$(MODULE_NAME)

ifeq ($(NOCOMPILE),)
  include $(MAKE_DIR)/setup-depends.mk

  ifeq ($(MODULE_NAME),)
    $(error Missing module name)
  endif

  INCLUDE_DIRS ?= -Iprivate -Ipublic
  ifeq ($(COMPILER),gfortran)
    include $(MAKE_DIR)/compile-settings-gnu.mk
  else
    $(error Unknown or unset COMPILER)
  endif

  include $(MAKE_DIR)/compile.mk
  include $(MAKE_DIR)/install-commands.mk

  ifeq ($(BINTYPE), lib)
    include $(MAKE_DIR)/link-library.mk
  else ifeq ($(BINTYPE), executable)
    include $(MAKE_DIR)/link-exec.mk
  else ifeq ($(BINTYPE), manual)
  else
    $(error Unknown or unset BINTYPE)
  endif

  ALL_DEPS += $(OBJ_OUT) $(MODULES) $(INSTALL_INCLUDES_BUILD)

  all: $(ALL_DEPS)

  install: $(INSTALL_COMMANDS)

  include $(MAKE_DIR)/check.mk
endif

extract-data: $(EXTRACT_DATA_COMMANDS)

clean:
	@rm -rf $(BUILD_DIR_MODULE)

print-mesa-deps:
	@echo $(INTERNAL_DEPENDS_ON)
