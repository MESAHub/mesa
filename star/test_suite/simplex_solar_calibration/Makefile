ifeq ($(MESA_DIR),)
  $(error MESA_DIR environment variable is not set)
endif

include $(MESA_DIR)/make/defaults-module.mk

MODULE_NAME := star
SRCS := \
   src/run.f90 \
   src/run_star_extras.f90 \
   src/simplex_search_data.f90 \
   src/simplex_search_run_support.f90

include $(MAKE_DIR)/work.mk

run: $(OBJ_OUT)
	@date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"
	@$(OBJ_OUT)
	@date "+DATE: %Y-%m-%d%nTIME: %H:%M:%S"

.PHONY: run
