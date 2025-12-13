include ../../make/defaults-module.mk

# Build

MODULE_NAME := eos-dt-builder
SRCS := src/azbar.f \
        src/create_eos_files.f \
        src/eos5_xtrin_h-he.f \
        src/helm_opal_scvh_driver.f \
        src/opal_core.f \
        src/opal_scvh_driver.f \
        src/scvh_core.f
SRCS_GENERATED := \
        private/helm_alloc.f90 \
        private/helm_polynomials.f90 \
        private/gauss_fermi.f90 \
        private/helm.f90
SRCS_CHECK :=
INTERNAL_DEPENDS_ON := const utils math chem num interp_2d eos
EXTERNAL_DEPENDS_ON :=
CHECK_RESULTS_GOLDEN :=
INCLUDE_DIRS := -Isrc -I../private

# Install

MODULES :=
INSTALL_INCLUDES :=
BINTYPE := executable

include ../../make/Makefile

# Generated sources
# This just copies over private implementations from eos

$(BUILD_DIR_MODULE)/private/helm_alloc.f90: ../private/helm_alloc.f90 | $(BUILD_DIR_MODULE)/private/
	@cp $^ $@

$(BUILD_DIR_MODULE)/private/helm_polynomials.f90: ../private/helm_polynomials.f90 | $(BUILD_DIR_MODULE)/private/
	@cp $^ $@

$(BUILD_DIR_MODULE)/private/gauss_fermi.f90: ../private/gauss_fermi.f90 | $(BUILD_DIR_MODULE)/private/
	@cp $^ $@

$(BUILD_DIR_MODULE)/private/helm.f90: ../private/helm.f90 | $(BUILD_DIR_MODULE)/private/
	@cp $^ $@
