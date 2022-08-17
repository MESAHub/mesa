# QUIET makes the output from the compilation process much less verbose -- basically,
# show which files are being compiled with which generic compilation command. The
# goal of QUIET is to make it easier to spot warnings and/or circular dependency
# errors.

include $(MESA_DIR)/utils/makefile_header

QUIET = yes

export QUIET


ifeq ($(USE_ADIPLS),YES)
all install:
	@rm -f ../test/build_skip_test
	@make -f makefile_adipls $@
	@make -f makefile_adiajobs $@
else
all install:
	@touch ../test/build_skip_test
endif


clean:
	@make -f makefile_adipls $@
	@make -f makefile_adiajobs $@
	@rm -f ../test/build_skip_test
