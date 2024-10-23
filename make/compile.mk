MODULE_SUBDIR := modules
MODULE_DIR := $(BUILD_DIR_MODULE)/$(MODULE_SUBDIR)
SCRATCH_DIR := $(BUILD_DIR_MODULE)/scratch

escape = '$(subst ','\'',$(1))'

define to_obj_ext
  $(patsubst %.f90,%.o, \
  $(patsubst %.F90,%.o, \
  $(patsubst %.f,%.o,	\
  $(patsubst %.c,%.o,	\
  $(addprefix $(BUILD_DIR_MODULE)/, $(1) \
  )))))
endef

OBJS = $(call to_obj_ext,$(SRCS) $(SRCS_GENERATED))
ifneq ($(SRCS_CHECK),)
  OBJS_CHECK = $(call to_obj_ext,$(SRCS_CHECK) $(SRCS_CHECK_GENERATED))
endif
MODULES := $(addprefix $(MODULE_DIR)/, $(MODULES))

FORTRAN_SOURCES := $(filter-out %.c, $(SRCS) $(SRCS_CHECK) $(addprefix $(BUILD_DIR_MODULE)/,$(SRCS_GENERATED) $(SRCS_CHECK_GENERATED)))

$(BUILD_DIR_MODULE)/depend : $(FORTRAN_SOURCES) | $(BUILD_DIR_MODULE)/
	INSTALL_INCLUDES=$(call escape,$(INSTALL_INCLUDES)) \
	MODULES=$(call escape,$(MODULES)) \
	BUILD_DIR_MODULE=$(call escape,$(BUILD_DIR_MODULE)) \
	$(MAKE_DIR)/gen-compile-tree  $(INCLUDE_DIRS) $(FLAGS_DEPS) $(FORTRAN_SOURCES) > $(BUILD_DIR_MODULE)/depend

include $(BUILD_DIR_MODULE)/depend

# Rule specific variables get filled too late for variable expansion
# in dependencies. The following option enables a second round of
# expansions when these variables are set.
.SECONDEXPANSION:

$(BUILD_DIR_MODULE)/%.o: %.c | $$(dir $$@0)/
	$(CCOMPILE) $< -o $@
