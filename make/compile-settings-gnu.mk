FLAGS_CODE_SANITY := \
  -Wall \
  -Wextra \
  -Wno-unused-label \
  -Wno-unused-parameter \
  -fstack-protector-all \
  -fstack-clash-protection \
  -fcheck=bounds \
  -D_FORTIFY_SOURCE=2
FFLAGS_FP_SANITY := -finit-real=snan -finit-derived -ffpe-trap=invalid,overflow,zero
FFLAGS_FORTRAN_SANITY := -std=f2008 -ffree-line-length-none -ffixed-line-length-none -Wno-unused-dummy-argument -Wno-compare-reals -Wno-do-subscript
FLAGS_REPRO := -ffp-contract=off
FFLAGS_PREPROCESSOR := -cpp

ifeq ($(PROFILE),release)
  FLAGS_OPT := -O2
  FLAGS_DEBUG :=
else ifeq ($(PROFILE),release-with-dbg-info)
  FLAGS_OPT := -O2
  FLAGS_DEBUG := -ggdb
else ifeq ($(PROFILE),debug)
  FLAGS_OPT := -Og
  FLAGS_DEBUG := -ggdb
else
  $(error Unknown or unset PROFILE)
endif

ifeq ($(WITH_OPENMP),yes)
  FLAGS_OPENMP := -fopenmp
else
  FLAGS_OPENMP :=
endif

# no-sign-zero only affects output formatting
FFLAGS_COMPAT := -fno-sign-zero
FFLAGS_FREE :=  -ffree-form -x f95-cpp-input -fimplicit-none
FFLAGS_FIXED := -ffixed-form -x f77-cpp-input
FLAGS_DEPS := $(call pkg-config,--cflags,$(DEPENDS_ON)) $(INCLUDE_DIRS)

FFLAGS_SHARED  := \
  $(FLAGS_CODE_SANITY) \
  $(FFLAGS_FP_SANITY) \
  $(FFLAGS_FORTRAN_SANITY) \
  $(FLAGS_REPRO) \
  $(FFLAGS_PREPROCESSOR) \
  $(FLAGS_OPT) \
  $(FLAGS_DEBUG) \
  $(FLAGS_OPENMP) \
  $(FFLAGS_COMPAT) \
  $(FLAGS_DEPS)
FFLAGS_FIXED := $(FFLAGS_SHARED) $(FFLAGS_FIXED) $(FFLAGS)
_FFLAGS := $(FFLAGS_SHARED) $(FFLAGS_FREE) $(FFLAGS)
_CFLAGS := $(FLAGS_CODE_SANITY) $(FLAGS_REPRO) $(FLAGS_OPT) $(FLAGS_DEBUG) $(FLAGS_DEPS) $(FLAGS_OPENMP) $(CFLAGS)

PREPROCESS := gfortran -cpp -E
FCOMPILE := gfortran $(_FFLAGS) -c
FCOMPILE_MODULE := gfortran $(_FFLAGS) -c -fsyntax-only
FCOMPILE_FIXED := gfortran $(FFLAGS_FIXED) -c
FCOMPILE_MODULE_FIXED:= gfortran $(FFLAGS_FIXED) -c -fsyntax-only
CCOMPILE := gcc $(_CFLAGS) -c
LIB_DEP_ARGS := $(call pkg-config, --libs,$(DEPENDS_ON))
LIB_TOOL_STATIC := ar rcs
LIB_TOOL_DYNAMIC := gfortran -shared $(FLAGS_OPENMP)
EXECUTABLE := gfortran $(FLAGS_OPENMP)
