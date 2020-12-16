
ROOT = ..
include $(ROOT)/../../utils/makefile_header

FFLAGS = -w

ADIR = $(ROOT)/adipls
GDIR = $(ROOT)/gensr
ADIAJOBS = .
ADISR = $(ADIAJOBS)/sr

LOFILES = \
$(GDIR)/ofiles.c.o \
$(GDIR)/lirt.d.o \
$(GDIR)/derive.d.o \
$(GDIR)/craysub.d.o \
$(GDIR)/iprsgn.d.o \
$(GDIR)/length.o \
$(GDIR)/leq.d.o \
$(GDIR)/lir.d.o \
$(GDIR)/lsqpol.d.o \
$(GDIR)/lstsqr.o \
$(GDIR)/blstio.o \
$(GDIR)/cpycom.o \
$(GDIR)/ofiles.o \
$(GDIR)/rdilst.o \
$(GDIR)/rnmean.d.o \
$(GDIR)/skpcom.o \
$(GDIR)/sort.d.o \
$(GDIR)/spline.d.o \
$(GDIR)/vinta.d.o \
$(GDIR)/squint.d.o \
$(GDIR)/getcas.o \
$(GDIR)/zero.d.o \
$(GDIR)/izero.o \
$(GDIR)/blstio.o

ALL = \
asymp-scale.d.x compamod.n.d.x diff-fgong.d.x diff-sum.d.x \
fdif.sin.d.x fgong-amdl.d.x \
filt-gsm.d.x fit-dnl.x form-agsm.d.x form-amdl.d.x freqdif.d.x \
lsqfreq.d.x mer-amde.d.x mer-sum.d.x prt-amdl.d.x redistrb.c.d.x \
redistrb.d.x res-amde.d.x scan-agsm.d.x scan-amde.d.x scan-amdl.d.x \
scan-ssm.d.x selamde.d.x sel-amdl.d.x selsum.d.x set-asscal.d.x set-dnl.d.x \
setexec.d.x set-obs.d.x sortsum.d.x 
 
default: all

all: $(ALL)

asymp-scale.d.x: asymp-scale.d.o srfiles lofiles
	$(FC) $(FFLAGS) -o asymp-scale.d.x asymp-scale.d.o \
	$(GDIR)/lir.d.o $(GDIR)/getcas.o $(ADISR)/rdamdl.n.d.o \
	$(GDIR)/skpcom.o $(GDIR)/cpycom.o $(GDIR)/length.o \
	$(GDIR)/ofiles.o $(GDIR)/blstio.o

compamod.n.d.x: compamod.n.d.o lofiles
	$(FC) $(FFLAGS) -o compamod.n.d.x compamod.n.d.o \
	$(GDIR)/ofiles.o $(GDIR)/lir.d.o \
	$(GDIR)/vinta.d.o $(GDIR)/length.o $(GDIR)/blstio.o

diff-fgong.d.x: diff-fgong.d.o lofiles srfiles
	$(FC) $(FFLAGS) -o diff-fgong.d.x diff-fgong.d.o $(GDIR)/lir.d.o $(GDIR)/blstio.o

diff-sum.d.x: diff-sum.d.o lofiles srfiles
	$(FC) $(FFLAGS) -o diff-sum.d.x diff-sum.d.o \
	$(ADIAJOBS)/sr/rdfreq.d.o $(ADIAJOBS)/sr/wrfreq.d.o \
	$(GDIR)/skpcom.o $(GDIR)/blstio.o

fdif.sin.d.x: fdif.sin.d.o lofiles
	$(FC) $(FFLAGS) -o fdif.sin.d.x fdif.sin.d.o \
	$(GDIR)/ofiles.o $(GDIR)/zero.d.o $(GDIR)/lir.d.o \
	$(GDIR)/skpcom.o $(GDIR)/length.o $(GDIR)/blstio.o

fgong-amdl.d.x: fgong-amdl.d.o
	$(FC) $(FFLAGS) -o fgong-amdl.d.x fgong-amdl.d.o $(GDIR)/blstio.o

filt-gsm.d.x: filt-gsm.d.o
	$(FC) $(FFLAGS) -o filt-gsm.d.x filt-gsm.d.o $(GDIR)/blstio.o

fit-dnl.x: fit-dnl.o lofiles
	$(FC) $(FFLAGS) -o fit-dnl.x fit-dnl.o \
	$(GDIR)/lstsqr.o $(GDIR)/skpcom.o $(GDIR)/length.o $(GDIR)/blstio.o

form-amdl.d.x: form-amdl.d.o
	$(FC) $(FFLAGS) -o form-amdl.d.x form-amdl.d.o $(GDIR)/blstio.o

form-agsm.d.x: form-agsm.d.o
	$(FC) $(FFLAGS) -o form-agsm.d.x form-agsm.d.o $(GDIR)/blstio.o

freqdif.d.x: freqdif.d.o lofiles srfiles
	$(FC) $(FFLAGS) -o freqdif.d.x freqdif.d.o \
	$(GDIR)/ofiles.o $(GDIR)/zero.d.o $(GDIR)/izero.o $(GDIR)/lir.d.o \
	$(ADIAJOBS)/sr/rdfreq.d.o $(ADIAJOBS)/sr/rdfrqe.d.o \
	$(GDIR)/skpcom.o $(GDIR)/length.o $(GDIR)/blstio.o

lsqfreq.d.x: lsqfreq.d.o lofiles srfiles
	$(FC) $(FFLAGS) -o lsqfreq.d.x lsqfreq.d.o \
	$(ADIAJOBS)/sr/rdfrqm.d.o $(GDIR)/ofiles.o \
	$(GDIR)/lsqpol.d.o $(GDIR)/leq.d.o $(GDIR)/skpcom.o \
	$(GDIR)/craysub.d.o $(GDIR)/length.o $(GDIR)/blstio.o

mer-amde.d.x: mer-amde.d.o
	$(FC) $(FFLAGS) -o mer-amde.d.x mer-amde.d.o $(GDIR)/blstio.o

mer-sum.d.x: mer-sum.d.o lofiles srfiles
	$(FC) $(FFLAGS) -o mer-sum.d.x mer-sum.d.o \
	$(ADIAJOBS)/sr/rdfreq.d.o $(ADIAJOBS)/sr/wrfreq.d.o \
	$(GDIR)/skpcom.o $(GDIR)/blstio.o

prt-amdl.d.x: prt-amdl.d.o
	$(FC) $(FFLAGS) -o prt-amdl.d.x prt-amdl.d.o $(GDIR)/blstio.o

redistrb.d.x: redistrb.d.o lofiles
	$(FC) $(FFLAGS) -o redistrb.d.x redistrb.d.o \
	$(GDIR)/ofiles.o $(GDIR)/length.o $(GDIR)/lir.d.o $(GDIR)/vinta.d.o \
	$(GDIR)/rnmean.d.o $(GDIR)/spline.d.o $(GDIR)/blstio.o

redistrb.c.d.x: redistrb.c.d.o srfiles lofiles
	$(FC) $(FFLAGS) -o redistrb.c.d.x redistrb.c.d.o \
	$(ADISR)/srdist.c.d.o $(ADISR)/rseta4.c.d.o \
	$(GDIR)/ofiles.c.o $(GDIR)/length.o \
	$(GDIR)/lir.d.o $(GDIR)/lirt.d.o $(GDIR)/vinta.d.o $(GDIR)/derive.d.o \
	$(GDIR)/rnmean.d.o $(GDIR)/blstio.o

res-amde.d.x: res-amde.d.o
	$(FC) $(FFLAGS) -o res-amde.d.x res-amde.d.o $(GDIR)/blstio.o

scan-agsm.d.x: scan-agsm.d.o
	$(FC) $(FFLAGS) -o scan-agsm.d.x scan-agsm.d.o $(GDIR)/blstio.o

scan-amde.d.x: scan-amde.d.o
	$(FC) $(FFLAGS) -o scan-amde.d.x scan-amde.d.o $(GDIR)/blstio.o

scan-amdl.d.x: scan-amdl.d.o
	$(FC) $(FFLAGS) -o scan-amdl.d.x scan-amdl.d.o $(GDIR)/blstio.o

scan-ssm.d.x: scan-ssm.d.o
	$(FC) $(FFLAGS) -o scan-ssm.d.x scan-ssm.d.o $(GDIR)/blstio.o

sel-amdl.d.x: sel-amdl.d.o
	$(FC) $(FFLAGS) -o sel-amdl.d.x sel-amdl.d.o $(GDIR)/blstio.o

selamde.d.x: selamde.d.o lofiles srfiles
	$(FC) $(FFLAGS) -o selamde.d.x selamde.d.o \
	$(GDIR)/ofiles.o $(GDIR)/length.o $(ADIAJOBS)/sr/inrnge.o \
	$(ADIAJOBS)/sr/rdfreq.d.o $(ADIAJOBS)/sr/setssm.d.o \
	$(GDIR)/skpcom.o $(GDIR)/blstio.o

selsum.d.x: selsum.d.o lofiles srfiles
	$(FC) $(FFLAGS) -o selsum.d.x selsum.d.o \
	$(GDIR)/ofiles.o $(GDIR)/length.o $(ADIAJOBS)/sr/inrnge.o \
	$(ADIAJOBS)/sr/rdfreq.d.o $(ADIAJOBS)/sr/setssm.d.o \
	$(GDIR)/skpcom.o $(GDIR)/blstio.o

set-asscal.d.x: set-asscal.d.o lofiles 
	$(FC) $(FFLAGS) -o set-asscal.d.x set-asscal.d.o \
	$(ADIAJOBS)/sr/rdamdl.n.d.o \
	$(GDIR)/ofiles.o $(GDIR)/length.o \
	$(GDIR)/squint.d.o $(GDIR)/vinta.d.o \
	$(GDIR)/lir.d.o $(GDIR)/blstio.o

set-dnl.d.x: set-dnl.d.o lofiles
	$(FC) $(FFLAGS) -o set-dnl.d.x set-dnl.d.o \
	sr/rdfreq.d.o $(GDIR)/skpcom.o $(GDIR)/blstio.o

set-obs.d.x: set-obs.d.o lofiles srfiles
	$(FC) $(FFLAGS) -o set-obs.d.x set-obs.d.o \
	$(ADIAJOBS)/sr/rdfreq.d.o $(GDIR)/skpcom.o $(GDIR)/blstio.o

setexec.d.x: setexec.d.o lofiles
	$(FC) $(FFLAGS) -o setexec.d.x setexec.d.o sr/rdfreq.d.o \
	$(GDIR)/ofiles.o $(GDIR)/length.o $(GDIR)/skpcom.o $(GDIR)/blstio.o

sortsum.d.x: sortsum.d.o lofiles
	$(FC) $(FFLAGS) -o sortsum.d.x sortsum.d.o sr/rdfreq.d.o \
	$(GDIR)/rdilst.o $(GDIR)/iprsgn.d.o $(GDIR)/sort.d.o \
	$(GDIR)/skpcom.o $(GDIR)/ofiles.o $(GDIR)/length.o $(GDIR)/blstio.o

lofiles:
	cd $(GDIR); make $(LOFILES)

srfiles:
	cd $(ADISR); make

clean:
	-@rm -f *.o *.x sr/*.o
