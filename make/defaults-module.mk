# A bit of sed fun to ensure MAKE_DIR does not end with a '/'
# This is not a problem, but it makes the output nicer
MAKE_DIR := $(shell echo $(dir $(lastword $(MAKEFILE_LIST))) | sed "s|/$$||")

BUILD_DIR ?= build
COMPILER ?= gfortran
VERSION ?= $(file <$(MAKE_DIR)../data/version_number)
PROFILE ?= release-with-dbg-info
WITH_OPENMP ?= yes
WITH_CRLIBM ?= yes
PREFIX ?= $(MAKE_DIR)/..
