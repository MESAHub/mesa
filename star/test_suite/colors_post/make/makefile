include $(MESA_DIR)/utils/makefile_header

LOAD_MESA = -L$(MESA_DIR)/lib $(LOAD_MESA_MACRO)

WORK_DIR = ..
WORK_SRC_DIR = $(WORK_DIR)/src

EXE = colors_post_proc

OBJS = run.o

$(EXE): $(OBJS)
	$(LOADER) $(FCopenmp) -o $(WORK_DIR)/$(EXE) $(OBJS) $(LOAD_MESA)

run.o: $(WORK_SRC_DIR)/run.f90
	$(WORK_COMPILE) $<

clean:
	-@rm -f *.o *.mod *.smod $(WORK_DIR)/$(EXE)
