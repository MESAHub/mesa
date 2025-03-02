include ../../make/defaults-module.mk

MODULE_NAME := neu-plotter
SRCS := src/neu_plotter.f90
INTERNAL_DEPENDS_ON := utils neu num
INCLUDE_DIRS :=
BINTYPE := executable

include ../../make/Makefile

plot: $(OBJ_OUT)
	$(OBJ_OUT)
	./plotter.py
