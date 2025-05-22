include ../../make/defaults-module.mk

MODULE_NAME := ionization-preprocessor
SRCS := src/def_args.f90 \
        src/def_const.f90 \
        src/def_type.f90 \
        src/lib_alert.f90 \
        src/lib_rkqsstep.f90 \
        src/lib_util.f90 \
        src/mod_bi.f90 \
        src/mod_density.f90 \
        src/mod_fast_tfdh.f90 \
        src/mod_gfdi.f90 \
        src/mod_run_tfdh.f90 \
        src/mod_tfdh.f90 \
        src/root_chi.f90 \
        src/root_vex.f90 \
	src/main.f90
SRCS_CHECK :=
INTERNAL_DEPENDS_ON := ionization
EXTERNAL_DEPENDS_ON :=
CHECK_RESULTS_GOLDEN :=
INCLUDE_DIRS := -Isrc

# Install

MODULES :=
INSTALL_INCLUDES :=
BINTYPE := executable

include ../../make/Makefile

build-data:

