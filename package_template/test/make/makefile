
MESA_DIR = ../../..

#################################################################

# STEP 1: set up the standard info

include $(MESA_DIR)/utils/makefile_header

#################################################################

# STEP 2: info about the libraries

LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -lpkg
LOAD_OTHER = -L$(MESA_LIB_DIR) $(LOAD_MESA_MICRO)
#### change LOAD_OTHER as needed.  e.g., $(LOAD_MESA_STAR)

#################################################################

# STEP 3: build the test program

TEST_OBJS = test_pkg.o
TEST = tester

$(TEST) : $(TEST_OBJS)
	$(LOADER) $(FCopenmp) -o $(TEST_DIR)/$(TEST) $(TEST_OBJS) \
	$(LOAD_LOCAL) $(LOAD_OTHER) $(LD_FLAGS)

#################################################################

%.o: $(TEST_SRC_DIR)/%.f90
	$(TEST_COMPILE) $(FCfree) $<

clean:
	-@rm -f *.o *.mod $(TEST_LIB) $(TEST_DIR)/$(TEST)
