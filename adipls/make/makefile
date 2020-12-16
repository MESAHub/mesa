# QUIET makes the output from the compilation process much less verbose -- basically,
# show which files are being compiled with which generic compilation command. The
# goal of QUIET is to make it easier to spot warnings and/or circular dependency
# errors.

QUIET = yes

export QUIET

all clean install:
	@make -f makefile_adipls $@
	@make -f makefile_adiajobs $@
