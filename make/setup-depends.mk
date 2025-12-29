INTERNAL_DEPENDS_ON := $(addprefix mesa-,$(INTERNAL_DEPENDS_ON))

include $(MAKE_DIR)/subdirs.mk

PKG_CONFIG_FLAGS =

ifneq ($(MESASDK_ROOT),)
  ifeq ($(WITH_CRLIBM),yes)
    PKG_CONFIG_FLAGS += --define-variable=MATH_SLOT=crmath
    MESASDK_MATH_SLOT = crmath
  else
    PKG_CONFIG_FLAGS += --define-variable=MATH_SLOT=default
    MESASDK_MATH_SLOT = default
  endif
  export MESASDK_MATH_SLOT
endif

PKG_CONFIG_PATH := $(shell BUILD_DIR=$(BUILD_DIR_) $(MAKE_DIR)/gen-pkgconfig-path $(SUBDIRS)):$(PKG_CONFIG_PATH)
export PKG_CONFIG_PATH

define pkg-config-inner
ifneq ($(2),)
  TMP_VAR := $$(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config $(PKG_CONFIG_FLAGS) $(1) $(2) || echo "failed")
  ifeq ($$(TMP_VAR),failed)
    $$(warning PKG_CONFIG_PATH is $$(PKG_CONFIG_PATH))
    $$(error pkg-config failed to find some of $(2), check PKG_CONFIG_PATH is correct)
  endif
else
  TMP_VAR :=
endif
endef

define pkg-config
  $(eval $(call pkg-config-inner,$(1),$(2)))$(TMP_VAR)
endef
