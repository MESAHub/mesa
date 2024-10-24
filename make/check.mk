CHECKER := $(BUILD_DIR_MODULE)/bin/checker
CHECK_RESULT_DIR := $(BUILD_DIR_MODULE)
CHECK_RESULTS := $(CHECK_RESULT_DIR)/check-results
CHECK_DIFF_PROG ?= diff -b

$(CHECKER): $(ALL_DEPS) $(OBJS_CHECK) | $(BUILD_DIR_MODULE)/bin/
	$(EXECUTABLE) -o $(CHECKER) $(OBJS_CHECK) $(call pkg-config,--libs,mesa-$(MODULE_NAME))

ifneq ($(OBJS_CHECK),)
  $(CHECK_RESULTS) : $(CHECKER) | $(CHECK_RESULT_DIR)/
	pushd test > /dev/null ; ../$(CHECKER) > ../$(CHECK_RESULTS); popd > /dev/null
	$(CHECK_DIFF_PROG) $(CHECK_RESULTS) $(CHECK_RESULTS_GOLDEN)

  check: $(CHECK_RESULTS)
else
  check:
	$(error No test sources specified)
endif
