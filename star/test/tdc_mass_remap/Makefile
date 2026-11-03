ifeq ($(MESA_DIR),)
  $(error MESA_DIR environment variable is not set)
endif

include $(MESA_DIR)/make/defaults-module.mk

MODULE_NAME := star
SRCS := test_tdc_mass_remap.f90

# The regression intentionally exercises the private implementation in this
# source checkout, including the module produced when the star library is built.
INCLUDE_DIRS = -I$(BUILD_DIR_)/star/modules

# This driver needs the star dependency closure, not astero or binary.
override SUBDIRS := const utils math mtx auto_diff num interp_1d interp_2d chem eos forum colors rates neu net kap ionization atm turb star_data star

include $(MAKE_DIR)/work.mk

run: $(OBJ_OUT)
	$(OBJ_OUT)

.PHONY: run
