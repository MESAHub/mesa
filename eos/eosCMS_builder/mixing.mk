include ../../make/defaults-module.mk

MODULE_NAME := eos-cms-mixing
SRCS := src/cms_mixing.f90
INTERNAL_DEPENDS_ON := const math
INCLUDE_DIRS := -Isrc -I../private
BINTYPE := executable

include ../../make/Makefile
