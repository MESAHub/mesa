include ../../make/defaults-module.mk

MODULE_NAME := stella-extras
SRCS := stella_extras.f90
SRCS_CHECK :=
INTERNAL_DEPENDS_ON := const math utils colors interp_1d interp_2d
EXTERNAL_DEPENDS_ON :=
CHECK_RESULTS_GOLDEN :=
INCLUDE_DIRS :=

# Install

MODULES :=
INSTALL_INCLUDES :=
BINTYPE := executable

include $(MAKE_DIR)/Makefile

run: $(OBJ_OUT)
	$(OBJ_OUT)
