#############COMPILER STUFF###############
FC=gfortran -std=legacy

WARN=-Wall
CHECK= #-pedantic -fcheck=all
OPT=-O2

DEBUG=-g -fbacktrace

COMPILE=$(FC) $(OPT) $(WARN) $(CHECK) $(DEBUG)

INCLUDE = ../include

LOAD_LIBS = -L../lib -lfree_eos -L$(MESASDK_ROOT)/lib -llapack -lblas

############SOURCE FILES##################

SRCS = abund_read.f abund_process.f free_eos_test_main.f

OBJS = $(patsubst %.f, %.o, $(SRCS))

TEST = free_eos_test
FD = fermi_dirac_test
EX = exchange_test
ROSE = test_rosenbrock

all: $(TEST) $(FD) $(EX) $(ROSE)

$(FD): fermi_dirac_test_main.f
	$(COMPILE) -I$(INCLUDE) -o $(FD) $< $(LOAD_LIBS)

$(EX): exchange_test_main.f
	$(COMPILE) -I$(INCLUDE) -o $(EX) $< $(LOAD_LIBS)

$(ROSE): test_rosenbrock.f
	$(COMPILE) -I$(INCLUDE) -o $(ROSE) $< $(LOAD_LIBS)

$(TEST): $(OBJS)
	$(COMPILE) -o $(TEST) $(OBJS) $(LOAD_LIBS)


%.o : %.f
	$(COMPILE) -I$(INCLUDE) -c $<


clean:
	rm -f $(TEST) $(FD) $(EX) $(ROSE) $(OBJS)
