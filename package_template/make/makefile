# This is the makefile for the package template

#################################################################
#
# PREAMBLE

include $(MESA_DIR)/utils/makefile_header

#################################################################
#
# SOURCE FILES

SRCS = pkg_mod.f90 \
   pkg_def.f90 \
   pkg_lib.f90

#################################################################
#
# TARGETS

LIB = libpkg.$(LIB_SUFFIX)
ifeq ($(USE_SHARED), YES)
    LIB_DEPENDS = 
endif
LIB_OBJS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SRCS)))
$(LIB) : $(LIB_OBJS)
	$(LIB_TOOL) $(LIB) $(LIB_OBJS) $(LIB_DEPENDS)

#################################################################
#
# COMPILATION RULES

COMPILE = $(COMPILE_TO_DEPLOY) $(FCfree)
#COMPILE = $(COMPILE_TO_TEST) $(FCfree)

COMPILE_FIXED = $(COMPILE_TO_DEPLOY) $(FCfixed)
#COMPILE_FIXED = $(COMPILE_TO_TEST) $(FCfixed)

COMPILE_CMD = $(COMPILE)

%.o : %.f
	$(COMPILE_CMD) $<
%.o : %.f90
	$(COMPILE_CMD) $<
%.mod : %.o
	@true
%.o : %.mod

clean:
	-@rm -f *.o *.mod *genmod.f90 *.so *.a .depend

install:
	@cp xxx_lib.mod $(MESA_DIR)/include
	@cp xxx_def.mod $(MESA_DIR)/include
	@cp libxxx.$(LIB_SUFFIX) $(MESA_DIR)/lib


#################################################################
#
# DEPENDENCIES

# Paths to search for source & include files

VPATH = $(MOD_PUBLIC_DIR):$(MOD_PRIVATE_DIR)

# Automatic dependency generation using makedepf90

nodeps : $(.DEFAULT_GOAL)

NODEPS = $(or $(filter nodeps,$(MAKECMDGOALS)),$(filter clean,$(MAKECMDGOALS)))

ifeq ($(NODEPS),)

  DEPEND_DIRS = $(MOD_PUBLIC_DIR):$(MOD_PRIVATE_DIR)

  .depend :
	makedepf90 -m %m.mod -I$(DEPEND_DIRS) $(SRCS) > .depend
  -include .depend

endif

### END ###
