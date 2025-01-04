include ../../make/defaults-module.mk

MODULE_NAME := atm-table-create
SRCS := src/create_table_atm.f90
INTERNAL_DEPENDS_ON := const utils math kap eos chem atm
INCLUDE_DIRS := -Isrc
BINTYPE := executable

include ../../make/Makefile
