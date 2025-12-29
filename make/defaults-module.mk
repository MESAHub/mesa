MAKE_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
MAKE_DIR := $(MAKE_DIR:/=)

ifeq ($(MAKE_DIR),)
  $(error extracting `make` directory based on defaults-module.mk location failed)
endif

BUILD_DIR ?= build
COMPILER ?= gfortran
VERSION ?= $(file <$(MAKE_DIR)../data/version_number)
PROFILE ?= release-with-dbg-info
WITH_OPENMP ?= yes
WITH_CRLIBM ?= yes
ifeq ($(MESA_FPE_CHECKS_ON),1)
WITH_FPE_CHECKS ?= yes
else
WITH_FPE_CHECKS ?= no
endif
WITH_GYRE ?= yes
WITH_ADIPLS ?= yes
WITH_PGSTAR ?= yes
WITH_COVERAGE ?= no
PREFIX ?= $(MAKE_DIR)/..
