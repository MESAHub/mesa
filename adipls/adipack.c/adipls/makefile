
ROOT = ..
include $(ROOT)/../../utils/makefile_header

FFLAGS = -w

GDIR = $(ROOT)/gensr
ADIAJOBS = $(ROOT)/adiajobs
ADISR = $(ADIAJOBS)/sr

LOFILES = \
$(GDIR)/aramax.d.o $(GDIR)/derive.d.o $(GDIR)/derivk.new.d.o $(GDIR)/leq.d.o \
$(GDIR)/lininh.new.d.o $(GDIR)/lininh4.d.o $(GDIR)/ofiles.c.o \
$(GDIR)/nrk.d.o $(GDIR)/vinta.d.o $(GDIR)/lir.d.o $(GDIR)/rnmean.d.o \
$(GDIR)/vintk.new.d.o $(GDIR)/store.d.o $(GDIR)/notwin.d.o $(GDIR)/notwni.o \
$(GDIR)/skpcom.o $(GDIR)/zero.d.o $(GDIR)/eiginh.d.o $(GDIR)/intgpt.d.o \
$(GDIR)/spline.d.o $(GDIR)/length.o $(GDIR)/strcompr.o $(GDIR)/lirt.d.o $(GDIR)/blstio.o

OFILES = \
adipls.c.d.o anldet.d.o cmparr.d.o danorm.d.o decicd.d.o \
dmpgsm.d.o dmpmde.d.o dmpsum.d.o geninh.d.o gm1ker.n.d.o gravpo.n.d.o \
kiner.n.d.o mchcff.d.o order.d.o takata.d.o readml.n.d.o res_adimod.c.d.o \
readsm.d.o adirhs.c.d.o rotker.c.d.o \
trnmod.n.d.o setbcs.n.d.o setssm.n.d.o signew.n.d.o rdfrqr.d.o varfrq.n.d.o \
xlmult.d.o nrkint.c.d.o nrtssl.d.o stevft.n.d.o shtint.c.d.o mchsol.n.d.o \
sclasl.c.d.o nrkm.d.o setsum.d.o sumout.d.o \
sigint.n.d.o sigitr.n.d.o sigitc.c.d.o sigsol.c.d.o sigout.c.d.o sigscn.c.d.o \
rsnfit.n.d.o res_adipar.c.d.o setups_adi.c.d.o set_rotation.c.d.o \
setobs_st.c.d.o dump_obs.c.d.o \
delta1.d.o delta2.d.o delta3.d.o \
sph.d.o uhx.d.o


REDIST_OBS = $(ADISR)/srdist.c.d.o $(ADISR)/rseta4.c.d.o





LIB = libadipls.$(LIB_SUFFIX)
LIB_OBJS = $(OFILES)
LIB_TOOL = ar crs

default: adipls.c.d.x

$(LIB) : $(LIB_OBJS) $(REDIST_OBS) lofiles
	$(LIB_TOOL) $(LIB) $(LIB_OBJS) $(REDIST_OBS) $(LOFILES)

#On windows the blstio common block isnt linked in properly 
ifeq ($(BUILD_ON_WINDOWS),YES)
   EXTRAS=$(GDIR)/blstio.o
endif
	
adipls.c.d.x: main.c.d.o $(LIB) lofiles user_adi.c.d.o
	$(FC) $(FFLAGS) -o adipls.c.d.x main.c.d.o $(LIB) user_adi.c.d.o $(EXTRAS)

lofiles:
	cd $(GDIR); make $(LOFILES)
	
#%.o: %.f
#	 $(FC) $(FFLAGS) -fbounds-check -Wuninitialized -Warray-bounds -g -c $<

# to use libadipls, you must provide your own user_adi.c.d.o
# the routine spcout_adi is called for each solution


clean:
	-@rm -f *.o $(LIB) *.so
