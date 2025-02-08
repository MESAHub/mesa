include ../../make/defaults-module.mk

MODULE_NAME := atm-ck-read
SRCS := src/ckread.f
INTERNAL_DEPENDS_ON := const
INCLUDE_DIRS := -Isrc
BINTYPE := executable

include ../../make/Makefile
