include ../../make/defaults-module.mk

MODULE_NAME := chem-preprocessor
SRCS := src/chem_support.f90 \
        src/create_table.f90
SRCS_CHECK :=
INTERNAL_DEPENDS_ON := const utils math chem
EXTERNAL_DEPENDS_ON :=
CHECK_RESULTS_GOLDEN :=
INCLUDE_DIRS := -Isrc

# Install

MODULES :=
INSTALL_INCLUDES :=
BINTYPE := executable

include ../../make/Makefile

$(BUILD_DIR_MODULE)/chem-out/isotopes.data: $(OBJ_OUT) inlist chem_input_data.tar.xz | $(BUILD_DIR_MODULE)/chem-out/ $(BUILD_DIR_MODULE)/chem-input/
	rm -rf $(BUILD_DIR_MODULE)/chem-input/*
	tar -xf chem_input_data.tar.xz --strip-components=1 -C $(BUILD_DIR_MODULE)/chem-input/
	$(OBJ_OUT) inlist

process-data: $(BUILD_DIR_MODULE)/chem-out/isotopes.data
	cp $^ ../data
