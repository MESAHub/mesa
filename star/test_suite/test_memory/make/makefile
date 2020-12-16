
MESA_DIR = ../../../..

# change this as necessary.  see utils/makefile_header for definitions.

include $(MESA_DIR)/utils/makefile_header

LOAD_MESA = $(LOAD_EXTRAS1) $(LOAD_MESA_STAR) $(LOAD_EXTRAS)

#################################################################

EXEC = test_memory
OBJS = test_memory.o

WORK_DIR = ..
WORK_SRC_DIR = $(WORK_DIR)/src
STAR_JOB_DIR = $(MESA_DIR)/star/job

$(EXEC) : $(OBJS)
	$(LOADER) $(FCopenmp) -o $(WORK_DIR)/$(EXEC) $(OBJS) $(LOAD_MESA)

#################################################################

# WORK_COMPILE defined in makefile_header

ifeq ($(SKIP_NAN_TRAPS),YES)
test_memory.o: $(WORK_SRC_DIR)/test_memory.f90
	$(WORK_COMPILE) $<
else
test_memory.o: $(WORK_SRC_DIR)/test_memory.f90
	$(WORK_COMPILE) $(FCtrapNANs) $<
endif

%.o: $(STAR_JOB_DIR)/%.f
	$(WORK_COMPILE) $<

clean:
	-@rm -f *.o *.mod $(WORK_DIR)/$(EXEC)

remk:
	-@rm -f $(OBJS) $(WORK_DIR)/$(EXEC)

