FLAGS_CODE_SANITY := \
  -Wall \
  -Wextra \
  -Wno-unused-label \
  -Wno-unused-parameter \
  -fstack-protector-all \
  -fstack-clash-protection \
  -fcheck=bounds \
  -D_FORTIFY_SOURCE=2 \
  -fPIC
FFLAGS_FP_SANITY := -finit-derived
ifeq ($(WITH_FPE_CHECKS),yes)
FFLAGS_FP_SANITY += -ffpe-trap=invalid,overflow,zero -finit-real=snan
endif
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
FFLAGS_FREE := -ffree-form -fimplicit-none
FFLAGS_FIXED := -ffixed-form
FLAGS_DEPS := $(call pkg-config,--cflags,$(EXTERNAL_DEPENDS_ON) $(INTERNAL_DEPENDS_ON)) $(INCLUDE_DIRS)

ifeq ($(WITH_COVERAGE),yes)
  FLAGS_COVERAGE := --coverage
else
  FLAGS_COVERAGE :=
endif

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
  $(FLAGS_DEPS) \
  $(FLAGS_COVERAGE)
FFLAGS_FIXED := $(FFLAGS_SHARED) $(FFLAGS_FIXED) $(FFLAGS)
_FFLAGS := $(FFLAGS_SHARED) $(FFLAGS_FREE) $(FFLAGS)
_CFLAGS := $(FLAGS_CODE_SANITY) $(FLAGS_REPRO) $(FLAGS_OPT) $(FLAGS_DEBUG) $(FLAGS_DEPS) $(FLAGS_OPENMP) $(FLAGS_COVERAGE) $(CFLAGS)

PREPROCESS := gfortran -cpp -E
FCOMPILE := gfortran $(_FFLAGS) -c
FCOMPILE_MODULE := gfortran $(_FFLAGS) -w -c -fsyntax-only
FCOMPILE_FIXED := gfortran $(FFLAGS_FIXED) -c
FCOMPILE_MODULE_FIXED:= gfortran $(FFLAGS_FIXED) -w -c -fsyntax-only
CCOMPILE := gcc $(_CFLAGS) -c
LIB_DEP_ARGS := $(call pkg-config, --libs,$(EXTERNAL_DEPENDS_ON)) $(call pkg-config, --libs --static,$(INTERNAL_DEPENDS_ON))
LIB_TOOL_STATIC := ar rcs
LIB_TOOL_DYNAMIC := gfortran -shared $(FLAGS_OPENMP) $(FLAGS_COVERAGE) -o
EXECUTABLE := gfortran $(FLAGS_OPENMP) $(FLAGS_COVERAGE)
