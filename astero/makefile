# astero/makefile

include $(MESA_DIR)/utils/makefile_header

ifeq ($(USE_GYRE),YES)
LOAD_GYRE = -lgyre $(LOAD_LAPACK) $(LOAD_BLAS)
else
LOAD_GYRE = 
endif

ifeq ($(USE_ADIPLS),YES)
ADIPLS_LIB = -ladipls
else
ADIPLS_LIB = 
endif

LOAD_MESA = $(LOAD_EXTRAS1) -L$(MESA_LIB_DIR) -lastero $(ADIPLS_LIB) $(LOAD_MESA_STAR) -lastero $(LOAD_GYRE) $(LOAD_EXTRAS)
# the replication of -lastero is required to deal with the lack of "hooks" in adipls.  
# it needs astero, and astero needs it and order of loading matters on linux.  so need to list it twice.
# check this again someday to see if the problem has gone away.

#################################################################

ifndef STAR
STAR = star
endif

# STAR_WORK_OBJS defined in makefile_header
OBJS = $(STAR_WORK_OBJS) run.o

WORK_DIR = ..
WORK_SRC_DIR = $(WORK_DIR)/src
STAR_JOB_DIR = $(MESA_DIR)/star/job

$(STAR) : $(OBJS)
	$(LOADER) $(FCopenmp) -o $(WORK_DIR)/$(STAR) $(OBJS) $(LOAD_MESA)

#################################################################

# WORK_COMPILE defined in makefile_header

ifeq ($(SKIP_NAN_TRAPS),YES)
run.o: $(WORK_SRC_DIR)/run.f90
	$(WORK_COMPILE) $<
else
run.o: $(WORK_SRC_DIR)/run.f90
	$(WORK_COMPILE) $(FCtrapNANs) $<
endif

ifeq ($(MESA_CONTRIB_DIR),)
run_star_extras.o: $(WORK_SRC_DIR)/run_star_extras.f90
	$(WORK_COMPILE) $<
else
run_star_extras.o: $(WORK_SRC_DIR)/run_star_extras.f90
	$(WORK_COMPILE) -I $(MESA_CONTRIB_DIR)/hooks $<
endif

%.o: $(STAR_JOB_DIR)/%.f90
	$(WORK_COMPILE) $<

clean:
	-@rm -f *.o *.mod $(WORK_DIR)/$(STAR)

remk:
	-@rm -f run.o $(WORK_DIR)/$(STAR)
