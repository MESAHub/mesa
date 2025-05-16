BUILD_DIR_ := $(call clean-comp,$(MAKE_DIR)/../$(BUILD_DIR))

$(BUILD_DIR_):
	mkdir -p $@

$(BUILD_DIR_)%/:
	mkdir -p $@
