include ../../make/defaults-module.mk

MODULE_NAME := eos-plotter
SRCS := src/eos_plotter.f90
INTERNAL_DEPENDS_ON := const utils math chem eos
INCLUDE_DIRS :=
BINTYPE := executable

include ../../make/Makefile

plot: $(OBJ_OUT)
	$(OBJ_OUT)
	./regions.py
	./plotter.py
