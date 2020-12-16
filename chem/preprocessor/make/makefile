# This is the makefile for creating the program that builds the chem table

PACKAGE_DIR = ../..
MESA_DIR = ../../..

#################################################################

# STEP 1: get the standard compiler setup

include $(MESA_DIR)/utils/makefile_header

#################################################################

# STEP 2: info about the libraries

OTHER_LIBS_DIR = $(MESA_DIR)/lib
INCLUDES = -I$(MESA_DIR)/include
LOAD_OTHER_LIBS = -L$(OTHER_LIBS_DIR) $(LOAD_MESA_NUMERICS)

#################################################################

# STEP 3: build it

OBJS = chem_support.o

CCHEM = cchem

all : $(CCHEM)

$(CCHEM) : $(OBJS) create_table.o
	$(LOADER) $(FCopenmp) -o ../$(CCHEM) $(OBJS) create_table.o \
	$(LOAD_OTHER_LIBS)

#################################################################

SRC_DIR = ../src

#COMPILE = $(COMPILE_NO_CHECKS)
COMPILE = $(COMPILE_TO_TEST)

%.o: $(SRC_DIR)/%.f90
	$(COMPILE) $(FCfree) $<

clean:
	-@rm -f *.o *.mod ../$(CCHEM)

