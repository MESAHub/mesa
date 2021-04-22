# This is the makefile for creating the program that builds the FreeEOS tables

MESA_DIR = ../../..

#################################################################

# STEP 1: get the standard compiler setup

include $(MESA_DIR)/utils/makefile_header

#################################################################

# STEP 2: info about the libraries

LOAD_OTHER_LIBS = -L$(MESA_LIB_DIR) $(LOAD_MESA_MICRO)

#################################################################

# STEP 3: build it

OBJS = cms_mixing.o cms_resample.o

MIX = cms_mixing
RES = cms_resample

all : $(MIX) $(RES)

$(MIX) : cms_mixing.o
	$(FC) $(FCopenmp) -o ../$@ $^ $(LOAD_OTHER_LIBS)

$(RES) : cms_resample.o
	$(FC) $(FCopenmp) -o ../$@ $^ $(LOAD_OTHER_LIBS)


#################################################################

SRC_DIR = ../src

COMPILE_SIMPLE = $(FC) -pedantic -Wall -fcheck=all -fbacktrace $(FCdebug) -I$(MESA_INCLUDE_DIR) -c

#COMPILE = $(FC) $(FCbasic) $(FCopenmp) $(FCopt) $(FCstatic) $(FCstandard) $(FCdebug) -I$(MESA_INCLUDE_DIR)  -c

cms_mixing.o: $(SRC_DIR)/cms_mixing.f90
	$(COMPILE_SIMPLE) $(FCfree) $^

cms_resample.o: $(SRC_DIR)/cms_resample.f90
	$(COMPILE_SIMPLE) $(FCfree) $^

clean:
	-@rm -f *.o *mod ../$(MIX) ../$(RES)
