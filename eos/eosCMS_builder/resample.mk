include ../../make/defaults-module.mk

MODULE_NAME := eos-cms-resample
SRCS := src/cms_resample.f90
INTERNAL_DEPENDS_ON := interp_1d
INCLUDE_DIRS := -Isrc -I../private
BINTYPE := executable

include ../../make/Makefile
