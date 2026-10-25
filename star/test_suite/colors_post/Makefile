ifeq ($(MESA_DIR),)
  $(error MESA_DIR environment variable is not set)
endif

include $(MESA_DIR)/make/defaults-module.mk

MODULE_NAME := colors_post
SRCS := src/run.f90

include $(MAKE_DIR)/work.mk

run: $(OBJ_OUT)
	@$(OBJ_OUT)

.PHONY: run
