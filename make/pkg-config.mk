PKG_CONFIG := $(BUILD_DIR_MODULE)/lib/pkgconfig/mesa-$(MODULE_NAME).pc
ALL_DEPS += $(PKG_CONFIG) $(foreach dir,$(patsubst mesa-%,%,$(INTERNAL_DEPENDS_ON)),$(subst %,$(dir),$(BUILD_DIR_)/%/lib/pkgconfig/mesa-%.pc))
INSTALL_COMMANDS += install-pkgconfig

define PKG_CONFIG_CONTENTS
prefix=$${pcfiledir}/../..

Version: $(VERSION)
Name: mesa-$(MODULE_NAME)
Requires.private: $(INTERNAL_DEPENDS_ON)
Description: MESA $(MODULE_NAME) module
Cflags: -I$${prefix}/include
Libs: -L$${prefix}/lib $(addprefix -l,$(LIB_NAMES))
Libs.private: $(call pkg-config, --libs,$(EXTERNAL_DEPENDS_ON))
endef

export PKG_CONFIG_CONTENTS

$(PKG_CONFIG): $(OBJ_OUT) | $(BUILD_DIR_MODULE)/lib/pkgconfig/. $(BUILD_DIR_MODULE)/include/.
	echo "$${PKG_CONFIG_CONTENTS}" > $@
