# Timing tests for star bcyclic solver
# Run from star/time/ directory

# Check that MESA_DIR is set
ifndef MESA_DIR
  $(error MESA_DIR environment variable is not set)
endif

FC = gfortran
FFLAGS = -O2 -fopenmp

# MESA include and library paths
# Include both public and private module paths for star_bcyclic
MESA_INCLUDE = -I$(MESA_DIR)/include -I$(MESA_DIR)/star/make
MESA_LIB = -L$(MESA_DIR)/lib

# MESASDK library paths
MESASDK_MATH_LIB = -L$(MESASDK_ROOT)/math-slots/default/lib
MESASDK_LIB = -L$(MESASDK_ROOT)/lib

# Libraries to link (star and its dependencies)
# Note: star library includes star_bcyclic module
MESA_LIBS = $(MESA_LIB) $(MESASDK_MATH_LIB) $(MESASDK_LIB) \
	-lstar -lstar_data -lgyre -latm -lturb -lionization -lkap -lnet -lneu -lrates -lcolors -lchem \
	-leos -lforum -linterp_2d -linterp_1d -lnum -lauto_diff -lmtx -lmath -lutils -lconst \
	-llapack -lblas -lcrmath

# Source files
SRCS = src/time_star_bcyclic.f90 src/run_timing.f90
OBJS = make/time_star_bcyclic.o make/run_timing.o

# Output executable
TARGET = timer

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJS) $(MESA_LIBS)

make/time_star_bcyclic.o: src/time_star_bcyclic.f90 | make/.
	$(FC) $(FFLAGS) $(MESA_INCLUDE) -c $< -o $@ -Jmake

make/run_timing.o: src/run_timing.f90 make/time_star_bcyclic.o | make/.
	$(FC) $(FFLAGS) $(MESA_INCLUDE) -c $< -o $@ -Imake -Jmake

make/.:
	mkdir -p make

clean:
	rm -rf make $(TARGET)

run: $(TARGET)
	./$(TARGET)
