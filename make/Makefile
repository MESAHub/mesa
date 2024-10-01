# Don't forget to update the documentation when you change something here!
# Makefile sanity
MAKEFLAGS += --no-builtin-rules
MAKEFLAGS += --no-builtin-variables

.DELETE_ON_ERROR:
.SHELLFLAGS := -eu -c
.DEFAULT_GOAL := all
.PHONY: clean print-mesa-deps install all check

# Position independent

MAKEFLAGS += -I$(MAKE_DIR)

NOCOMPILE := $(or \
  $(filter clean,$(MAKECMDGOALS)), \
  $(filter print-mesa-deps,$(MAKECMDGOALS)), \
)

include helpers.mk
include setup-builddir.mk
BUILD_DIR_MODULE := $(BUILD_DIR_)/$(MODULE_NAME)

ifeq ($(NOCOMPILE),)
  include setup-depends.mk

  ifeq ($(MODULE_NAME),)
    $(error Missing module name)
  endif

  INCLUDE_DIRS ?= -Iprivate -Ipublic
  ifeq ($(COMPILER),gfortran)
    include compile-settings-gnu.mk
  else
    $(error Unknown or unset COMPILER)
  endif

  include compile.mk
  include install-commands.mk

  ifeq ($(BINTYPE), static-lib)
    include link-static.mk
  else ifeq ($(BINTYPE), dynamic-lib)
    include link-dynamic.mk
  else ifeq ($(BINTYPE), executable)
    include link-exec.mk
  else
    $(error Unknown or unset BINTYPE)
  endif

  ALL_DEPS += $(OBJ_OUT) $(MODULES) $(INSTALL_INCLUDES_BUILD)

  all: $(ALL_DEPS)

  install: $(INSTALL_COMMANDS)

  include check.mk
endif

clean:
	@rm -rf $(BUILD_DIR_MODULE)

print-mesa-deps:
	@echo $(INTERNAL_DEPENDS_ON)
