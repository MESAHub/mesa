ifeq ($(DYNAMIC),yes)
  UNAME := $(shell uname)
  LIB_TOOL := $(LIB_TOOL_DYNAMIC)
  LIB_TOOL_ARGS := $(LIB_DEP_ARGS)
  ifeq ($(UNAME),Linux)
    LIB_EXT := so
  else ifeq ($(UNAME),Darwin)
    LIB_EXT := dylib
    LIB_TOOL_ARGS += -install_name @rpath/lib$(MODULE_NAME).$(LIB_EXT)
  else
    $(warning unknown OS type, library extension could be wrong)
    LIB_EXT := so
  endif
else
  LIB_EXT := a
  LIB_TOOL := $(LIB_TOOL_STATIC)
  LIB_TOOL_ARGS :=
endif

LIB_NAMES ?= $(MODULE_NAME)
OBJ_OUT := $(patsubst %,$(BUILD_DIR_MODULE)/lib/lib%.$(LIB_EXT),$(LIB_NAMES))
INSTALL_COMMANDS += install-lib

$(OBJ_OUT): $(OBJS) | $(BUILD_DIR_MODULE)/lib/.
	$(LIB_TOOL) $(OBJ_OUT) $(OBJS) $(LIB_TOOL_ARGS)

include $(MAKE_DIR)/pkg-config.mk

