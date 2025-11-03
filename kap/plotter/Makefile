include ../../make/defaults-module.mk

MODULE_NAME := kap-plotter
SRCS := src/kap_plotter.f90
SRCS_CHECK :=
INTERNAL_DEPENDS_ON := kap
EXTERNAL_DEPENDS_ON :=
CHECK_RESULTS_GOLDEN :=
INCLUDE_DIRS := -Isrc

# Install

MODULES :=
INSTALL_INCLUDES :=
BINTYPE := executable

include ../../make/Makefile

kap_plotter.dat: $(OBJ_OUT)
	$(OBJ_OUT)

plot: kap_plotter.dat
	./plotter.py

.PHONY: plot
