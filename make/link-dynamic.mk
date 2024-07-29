OBJ_OUT := $(BUILD_DIR_MODULE)/lib/lib$(MODULE_NAME).so
INSTALL_COMMANDS += install-lib

$(OBJ_OUT): $(OBJS) | $(BUILD_DIR_MODULE)/lib/
	$(LIB_TOOL_DYNAMIC) -o $(OBJ_OUT) $(OBJS) $(LIB_DEP_ARGS)

include $(MAKE_DIR)/pkg-config.mk
