ifeq ($(MESA_DIR),)
  $(error MESA_DIR enviroment variable is not set)
endif

include $(MESA_DIR)/make/defaults-module.mk

MODULE_NAME := star
SRCS := \
   src/run.f90 \
   src/run_star_extras.f90 \

include $(MAKE_DIR)/work.mk

run: $(OBJ_OUT)
	@date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"
	@$(OBJ_OUT)
	@date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"

.PHONY: run
