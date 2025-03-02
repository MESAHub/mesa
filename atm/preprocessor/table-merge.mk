include ../../make/defaults-module.mk

MODULE_NAME := atm-table-merge
SRCS := src/table_merge.f90
INTERNAL_DEPENDS_ON := const
INCLUDE_DIRS := -Isrc
BINTYPE := executable

include ../../make/Makefile
