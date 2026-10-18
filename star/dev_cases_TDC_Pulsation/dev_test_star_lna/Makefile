ifeq ($(MESA_DIR),)
  $(error MESA_DIR environment variable is not set)
endif

include $(MESA_DIR)/make/defaults-module.mk

MODULE_NAME := star
SRCS := src/run.f90 src/run_star_extras.f90

include $(MAKE_DIR)/work.mk
