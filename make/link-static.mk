OBJ_OUT := $(BUILD_DIR_MODULE)/lib/lib$(MODULE_NAME).a
INSTALL_COMMANDS += install-lib

$(OBJ_OUT): $(OBJS) | $(BUILD_DIR_MODULE)/lib/
	$(LIB_TOOL_STATIC) -o $(OBJ_OUT) $(OBJS)

include $(MAKE_DIR)/pkg-config.mk
