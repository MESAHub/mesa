include ../../make/defaults-module.mk

MODULE_NAME := atm-ck-read
SRCS := src/ckread.f90
INTERNAL_DEPENDS_ON := const
INCLUDE_DIRS := -Isrc
BINTYPE := executable

include ../../make/Makefile
