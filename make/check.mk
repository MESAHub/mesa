CHECKER := $(BUILD_DIR_MODULE)/bin/checker
CHECK_RESULT_DIR := $(BUILD_DIR_MODULE)
CHECK_RESULTS := $(CHECK_RESULT_DIR)/check-results
CHECK_DIFF_PROG ?= diff -b
CHECK_FILTER_CACHE_WRITE ?= yes

$(CHECKER): $(ALL_DEPS) $(OBJS_CHECK) | $(BUILD_DIR_MODULE)/bin/
	$(EXECUTABLE) -o $(CHECKER) $(OBJS_CHECK) $(call pkg-config,--libs,mesa-$(MODULE_NAME))

ifeq ($(CHECK_FILTER_CACHE_WRITE),yes)
CHECK_FILTER_PROG ?= grep -Ev "^(write)|(create rate data for)|( read )|( write )|(                           number not already in cache:)"
else
CHECK_FILTER_PROG ?= cat
endif

ifneq ($(OBJS_CHECK),)
  $(CHECK_RESULTS) : $(CHECKER) | $(CHECK_RESULT_DIR)/
	cd test; ../$(CHECKER) | $(CHECK_FILTER_PROG) > ../$(CHECK_RESULTS)
	$(CHECK_DIFF_PROG) $(CHECK_RESULTS) $(CHECK_RESULTS_GOLDEN)

  check: $(CHECK_RESULTS)
else
  check:
endif
