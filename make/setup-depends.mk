DEPENDS_ON := $(EXTERNAL_DEPENDS_ON)

ifneq ($(INTERNAL_DEPENDS_ON),)
  DEPENDS_ON += $(addprefix mesa-,$(INTERNAL_DEPENDS_ON))
endif

include $(MAKE_DIR)/subdirs.mk

PKG_CONFIG_PATH := $(shell BUILD_DIR=$(BUILD_DIR_) $(MAKE_DIR)/gen-pkgconfig-path $(SUBDIRS)):$(PKG_CONFIG_PATH)
export PKG_CONFIG_PATH

define pkg-config-inner
ifneq ($(2),)
  TMP_VAR := $$(shell pkg-config $(1) $(2))
  ifneq ($$(.SHELLSTATUS),0)
    $$(warning PKG_CONFIG_PATH is $$(PKG_CONFIG_PATH))
    $$(error pkg-config failed to find some of $(2), check PKG_CONFIG_PATH is correct)
  endif
endif
endef

define pkg-config
  $(eval $(call pkg-config-inner,$(1),$(2)))$(TMP_VAR)
endef
