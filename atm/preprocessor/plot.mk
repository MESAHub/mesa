include ../../make/defaults-module.mk

MODULE_NAME := atm-plot
SRCS := src/plot.f90
INTERNAL_DEPENDS_ON := const utils
INCLUDE_DIRS := -Isrc
BINTYPE := executable

include ../../make/Makefile
