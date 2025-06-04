install-pkgconfig:
	install -pd $(PREFIX)/lib/pkgconfig
	install -pm 644 $(PKG_CONFIG) $(PREFIX)/lib/pkgconfig

install-lib:
	install -pd $(PREFIX)/lib
	install -p $(OBJ_OUT) $(PREFIX)/lib

install-exec:
	install -pd $(PREFIX)/bin
	install -p $(OBJ_OUT) $(PREFIX)/bin

INSTALL_INCLUDES_BUILD := $(addprefix $(BUILD_DIR_MODULE)/include/,$(notdir $(INSTALL_INCLUDES) $(MODULES)))

install-includes:
	install -pd $(PREFIX)/include
	if [ "$(INSTALL_INCLUDES_BUILD)" ]; then \
	  for f in "$(INSTALL_INCLUDES_BUILD)"; do \
	    install -pm 644 $$f $(PREFIX)/include; \
	  done; \
	fi

INSTALL_COMMANDS += install-includes
