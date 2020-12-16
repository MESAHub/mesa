
#################################################################


MESA_DIR = ../..

include $(MESA_DIR)/utils/makefile_header


#################################################################

PROG = stella_extras
PROG_OBJS = stella_extras.o
PROG_DIR = .


$(PROG) : $(PROG_OBJS)
	$(LOADER) $(FCopenmp) -o $(PROG_DIR)/$(PROG) $(PROG_OBJS) \
	-L$(MESA_LIB_DIR) $(LOAD_MESA_MACRO)

# for more options, see utils/makefile_header

#################################################################

MY_FC_FLAGS = $(FCfree)
SRC_DIR = .

COMPILE = $(FC) $(FCbasic) $(FCopenmp) $(FCstatic)  $(OTHER_INCLUDES) $(FCwarn) $(FCimpno) $(FCchecks) $(FCopt) $(FCdebug) -c

%.o: $(SRC_DIR)/%.f90
	$(COMPILE) $(MY_FC_FLAGS) $<

clean:
	-@rm -f *.o *.mod $(PROG_DIR)/$(PROG)
