include ../../../make/defaults-module.mk

MODULE_NAME := eos-composition-partials-check
SRCS := src/eos_composition_partials_check.f90
INTERNAL_DEPENDS_ON := const utils math chem eos
INCLUDE_DIRS :=
BINTYPE := executable

include ../../../make/Makefile

plot: $(OBJ_OUT)
	mkdir -p data figures data/.mplconfig
	$(OBJ_OUT)
	MPLCONFIGDIR=data/.mplconfig python plot_composition_partials.py
