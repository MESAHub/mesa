include ../../make/defaults-module.mk

MODULE_NAME := atm-bt-settl-create
SRCS := src/create_BT_Settl.f90
INTERNAL_DEPENDS_ON := const math chem num
INCLUDE_DIRS := -Isrc
BINTYPE := executable

include ../../make/Makefile
