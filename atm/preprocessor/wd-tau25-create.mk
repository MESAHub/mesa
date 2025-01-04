include ../../make/defaults-module.mk

MODULE_NAME := atm-wd-tau25-create
SRCS := src/create_wd_tau25.f90
INTERNAL_DEPENDS_ON := const utils
INCLUDE_DIRS := -Isrc
BINTYPE := executable

include ../../make/Makefile
