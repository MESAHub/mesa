# This is the makefile for creating the program that builds the FreeEOS tables

MESA_DIR = ../../..
FREE_EOS_DIR = $(MESA_DIR)/eos/eosFreeEOS_builder/free_eos-2.2.1

#################################################################

# STEP 1: get the standard compiler setup

include $(MESA_DIR)/utils/makefile_header

#################################################################

# STEP 2: info about the libraries

LOAD_OTHER_LIBS = -L$(FREE_EOS_DIR)/lib -lfree_eos -L$(MESA_LIB_DIR) $(LOAD_MESA_MICRO)

#################################################################

# STEP 3: build it

OBJS = free_eos_table.o

TABLE = free_eos_table

all : $(TABLE)

$(TABLE) : $(OBJS)
	$(COMPILE_BASIC) -o ../$(TABLE) $(OBJS) $(LOAD_OTHER_LIBS)


#################################################################

SRC_DIR = ../src

COMPILE = $(FC) $(FCbasic) $(FCopenmp) $(FCopt) $(FCstatic) $(FCstandard) $(FCdebug) -I$(MESA_INCLUDE_DIR)  -c

free_eos_table.o: $(SRC_DIR)/free_eos_table.f90
	$(COMPILE) $(FCfree) $^

clean:
	-@rm -f *.o *mod ../$(TABLE)
