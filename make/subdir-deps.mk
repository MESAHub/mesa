NODEPS := $(filter clean,$(MAKECMDGOALS))

ifeq ($(NODEPS),)
$(BUILD_DIR_)/depend: $(addsuffix /Makefile,$(SUBDIRS)) make/gen-folder-deps | $(BUILD_DIR_)
	make/gen-folder-deps "$(MAKE)" $(SUBDIRS) > $(BUILD_DIR_)/depend
include $(BUILD_DIR_)/depend
endif
