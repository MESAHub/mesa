# A bit of sed fun to ensure MAKE_DIR does not end with a '/'
# This is not a problem, but it makes the output nicer
MAKE_DIR := $(shell echo $(dir $(lastword $(MAKEFILE_LIST))) | sed "s|/$$||")

ifneq ($(.SHELLSTATUS),0)
  $(error extracting `make` directory based on defaults-module.mk location failed)
endif

BUILD_DIR ?= build
COMPILER ?= gfortran
VERSION ?= $(file <$(MAKE_DIR)../data/version_number)
PROFILE ?= release-with-dbg-info
WITH_OPENMP ?= yes
WITH_CRLIBM ?= yes
WITH_GYRE ?= yes
WITH_PGSTAR ?= yes
PREFIX ?= $(MAKE_DIR)/..
