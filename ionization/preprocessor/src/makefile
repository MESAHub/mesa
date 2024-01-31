FC=gfortran
FCFLAGS=-ffree-form -O2 #-std=f2003
LDFLAGS=

SOURCES=def_args.f\
def_const.f\
def_type.f\
lib_alert.f\
lib_rkqsstep.f\
lib_util.f\
mod_bi.f\
mod_density.f\
mod_fast_tfdh.f\
mod_gfdi.f\
mod_run_tfdh.f\
mod_tfdh.f\
root_chi.f\
root_vex.f
OBJECTS=$(SOURCES:.f=.o)

SRCMAIN=main.f
OBJMAIN=$(SRCMAIN:.f=.o)
EXECUTABLE=tfdh

SRCTEST=maintest.f
OBJTEST=$(SRCTEST:.f=.o)
TESTEXECUTABLE=tfdh_test

#################################################

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) $(OBJMAIN)
	$(FC) $(FCFLAGS) $(OBJECTS) $(OBJMAIN) -o $(EXECUTABLE)


test: $(TESTEXECUTABLE)

$(TESTEXECUTABLE): $(OBJECTS) $(OBJTEST)
	$(FC) $(FCFLAGS) $(OBJECTS) $(OBJTEST) -o $(TESTEXECUTABLE)


.f.o:
	$(FC) $(FCFLAGS) -c $<


clean:
	rm -f *.o *.mod $(EXECUTABLE) $(TESTEXECUTABLE)


#################################################
# dependencies in the code
main.o: def_args.o def_const.o def_type.o lib_util.o mod_run_tfdh.o
maintest.o: def_args.o def_const.o def_type.o lib_util.o mod_run_tfdh.o
mod_run_tfdh.o: def_args.o def_const.o def_type.o lib_util.o mod_tfdh.o root_chi.o root_vex.o
mod_tfdh.o: def_args.o def_const.o def_type.o lib_alert.o lib_rkqsstep.o lib_util.o mod_density.o
root_vex.o: def_args.o def_type.o mod_fast_tfdh.o
mod_fast_tfdh.o: def_args.o def_const.o def_type.o lib_rkqsstep.o mod_density.o
root_chi.o: def_type.o mod_density.o
mod_density.o: lib_alert.o def_args.o def_type.o mod_gfdi.o mod_bi.o
mod_bi.o: def_args.o def_type.o lib_alert.o lib_rkqsstep.o
lib_rkqsstep.o: def_args.o def_type.o lib_alert.o
mod_gfdi.o: def_type.o def_const.o lib_alert.o
lib_util.o: def_type.o
def_const.o: def_type.o 
def_args.o: def_type.o
def_type.o:
lib_alert.o:
#################################################


