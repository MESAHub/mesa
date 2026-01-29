.DELETE_ON_ERROR:
.DEFAULT_GOAL := all
.SHELLFLAGS := -eu -c

BUILD_SUBDIR := @$(MAKE) -C
CHECK_SUBDIR := @$(MAKE) check -C
MAKE_DIR := make

BUILD_DIR ?= build

include $(MAKE_DIR)/helpers.mk
include $(MAKE_DIR)/subdirs.mk
include $(MAKE_DIR)/setup-builddir.mk
include $(MAKE_DIR)/version.mk
include $(MAKE_DIR)/subdir-deps.mk

clean:
	@rm -rf $(BUILD_DIR_)

all-build: $(SUBDIRS)
all-extract-data: $(addsuffix -extract-data,$(SUBDIRS))
all: $(addsuffix -check,$(SUBDIRS))

.PHONY: all all-build all-extract-data $(SUBDIRS) $(addsuffix -extract-data,$(SUBDIRS)) $(addsuffix -check,$(SUBDIRS)) $(addsuffix -install,$(SUBDIRS)) clean
