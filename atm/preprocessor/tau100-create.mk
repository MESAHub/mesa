include ../../make/defaults-module.mk

MODULE_NAME := atm-tau100-create
SRCS := src/create_tau100.f90
INTERNAL_DEPENDS_ON := const utils math kap eos chem atm
INCLUDE_DIRS := -Isrc
BINTYPE := executable

include ../../make/Makefile
