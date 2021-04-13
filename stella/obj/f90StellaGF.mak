##############################
# MAKE FOR STELLA
##############################

# for standalone use define HOMEStella:
# tcsh:
# setenv HOMEStella `pwd`
# in ../ or
# bash:
# export HOMEStella=`pwd`
# in ../

HOMEStella := ../


#sys := $(shell uname)
#ifneq (,$(findstring CYGWIN,$(strip $(sys))))
#   SYSTYPE := "cygwin_f90"
#else
#  SYSTYPE := "ifort"
#  SYSTYPE="pgf"
#  SYSTYPE := "MPA_f95i"
#endif
#SYSTYPE="ifort"
# SYSTYPE := "lf"
SYSTYPE=gfortran
# SYSTYPE="cygwin_ifort"
# SYSTYPE="cygwin_f90"
# SYSTYPE="MPA_f95i"
# SYSTYPE="MPA_f95n"
UnixOrWin = "unix" # default value

########################################
## For Intel fortran
ifeq ($(SYSTYPE),"ifort")
  FC = ifort  # Intel fortran
#   FFLAGS = -c -save -zero -O3 -mcmodel=medium -shared-intel
#   FFLAGS = -c -save -zero -O3 -fp-model precise
#   FFLAGS = -c -save -zero -O3 -fp-model strict -CB -g -traceback # -- for debug
#   FFLAGS = -c -save -zero -fpe0 -no-ftz  -check bounds -g -traceback # -- for debug
#   FFLAGS = -c -save -zero -O3 -fp-model strict
#   FFLAGS = -c -O3 -save -zero -tpp7   -ip # -xW # -- optimize for Pentium
#   FFLAGS =   -c -save -static -zero -O2 -fpe0 -traceback # -CB -traceback
 FFLAGS_FIX = -c -132 -save -zero -O3 -fp-model strict
#  FFLAGS_FIX = -c -132 -save -zero -fp-model strict -CB -g -traceback -debug # inline-debug-info # -- for debug
  FFLAGS := -free
  FFLAGS := $(FFLAGS_FIX) $(FFLAGS)
#define
  FFLAGS := $(FFLAGS) -D__INTEL

#  FFLAGS = -c -save -zero -fpe0 -no-ftz  -check bounds -g -traceback # -- for debug
#  FFLAGS = -c -save -zero -check bounds -g -inline_debug_info -traceback # -- for debug
  #FFLAGS = -c -O3 -save -zero -tpp7   -ip # -xW # -- optimize for Pentium
  #FFLAGS =   -c -save -static -zero -O2 -fpe0 -traceback # -CB -traceback
  #LDFLAGS = -i-static
#  LDFLAGS = -mcmodel=medium -shared-intel # for 64b
#  FFLAGS = -c -save -fpe0 -g -CB -traceback -inline_debug_info  -DD
  # LIBS = /afs/mpa/HOMEStella/seb/lib/libsparsifc.a /afs/mpa/HOMEStella/seb/lib/liblnagifc.a
  #LIBS = /HOMEStella/blinn/lib/libsparseifc.a /HOMEStella/blinn/lib/liblnagifc.a
  LIBS= $(HOMEStella)/lib/libsparse.a $(HOMEStella)/lib/liblnag.a
  #LIBS= $(HOME)/lib/libsparse.a $(HOME)/lib/liblnag.a
  LIBSM = /usr/local/lib/libplotsub.a \
        /usr/local/lib/libdevices.a \
        /usr/local/lib/libutils.a -lX11
  UnixOrWin := "unix"
endif

########################################
## For Windows Intel fortran
ifeq ($(SYSTYPE),cygwin_ifort)
  FC = ifl  # Intel fortran
  FFLAGS = -c -save -zero  -tpp7   -ip  -static -fast
#  FFLAGS = -c -O3 -save -zero -tpp7   -ip
#  FFLAGS = -compile_only -static   -architecture:k7 -optimize:5 -fast -tune:k7 # for optimization Pentium
#  FFLAGS = -c -static -debug:full -check:bounds,power,overflow -traceback -list -show:all
  LDFLAGS =
#  LDFLAGS = /libs:static
  LIBS = ..\\lib\\liblnag.a	..\\lib\\libspars.a
  LIBDIRSYS = /I"/cygdrive/c/Program Files/Intel/Compiler70/IA32/Lib/"
#  LIBDIRSYS = /I"/cygdrive/c/Program Files/Intel/Compiler70/IA32/Lib/"
  UnixOrWin := "win"
endif

########################################
## For Visual Compaq Fortran
ifeq ($(SYSTYPE),cygwin_f90)
  FC = f90  # Compaq fortran
  LINK = link.exe
  # FFLAGS = -compile_only -static  -fast -tune:k7
  LDFLAGS = -architecture:k7 -optimize:5 -fast -tune:k7 # /ignore:505
  FFLAGS = -compile_only -static   -architecture:k7 -optimize:5 -fast -tune:k7 # for optimization Pentium
#  FFLAGS = -c -static -debug:full -check:bounds,power,overflow -traceback -list -show:all
  LIBS = d:\\prg\\lib\\libsparscvf.a d:\\prg\\lib\\liblnagcvf.a
#  LIBS = ..\\lib\\liblnagcvf.a	..\\lib\\libsparscvf.a
  LIBDIRSYS =
#  LIBDIRSYS = /libpath:"c:\\Program Files\\Microsoft Visual Studio\\DF98\\LIB\\" /libpath:"c:\\Program Files\\Microsoft Visual Studio\\VC98\\LIB\\"

#  LIBDIRSYS = "c:/Program Files/Microsoft Visual Studio/DF98/LIB/"
  UnixOrWin := "win"
endif

########################################
## For Intel fortran
ifeq ($(SYSTYPE),MPA_f95i)
  FC = f95i  # Intel fortran
  FFLAGS = -c -static -save -zero -fpe0  -check bounds -g -traceback # -- for debug
  FFLAGS = -c -O3 -save -zero -tpp7   -ip # -xW # -- optimize for Pentium
  LIBS= $(HOMEStella)/lib/libsparse.a $(HOMEStella)/lib/liblnag.a
  # LIBS = /afs/ipp/HOMEStella/s/seb/lib/libsparsifc.a  /afs/ipp/HOMEStella/s/seb/lib/liblnagifc.a
  # LIBS = /afs/mpa/HOMEStella/seb/lib/libsparsifc.a /afs/mpa/HOMEStella/seb/lib/liblnagifc.a
  LDFLAGS =
  UnixOrWin := "unix"
endif

########################################
## For NAG fortran
ifeq ($(SYSTYPE),MPA_f95n)
  FC = f95n  # NAG fortran
  FFLAGS = -c  -save
  FFLAGS = -c -O3 -save -f77
  LDFLAGS =
  LIBS = ../lib/libsparsNag.a ../lib/liblnagNag.a	# for NAG
  UnixOrWin := "unix"
endif
########################################
########################################

########################################
## For Portland Group fortran
ifeq ($(SYSTYPE),pgf)
  FC = pgf95  # pgf fortran
  FFLAGS = -c -fast -Minform,warn -Msave
#  FFLAGS = -c -C -g -Ktrap=fp -Mbounds -Minform,warn -Msave
  LIBS= $(HOMEStella)/lib/libsparse.a $(HOMEStella)/lib/liblnag.a
  UnixOrWin := "unix"
endif
########################################
## For Lahey-Fujutsu fortran
ifeq ($(SYSTYPE),"lf")
  FC = lf95  # Lahey-fujutsu fortran
#  FFLAGS = -c  -sav -g  --trap  --ap   --chkglobal   --pca # debug
  FFLAGS = -c  -sav --ap --tp4 --sse2 --zfm --o2
# --tp4 --sse2 --zfm --o2
# -x name.f  for inlining code name.f
#  FFLAGS = -c  -sav --ap -O --ml cdecl
#  FFLAGS = -c  -sav --trap  --tpp   --o2   --ntrace   --f95   --info
#  FFLAGS = -c  -sav  --tpp   --o2   --ntrace   --f95   --info
  LIBS =  $(HOME)/lib/libsparseLF.a $(HOME)/lib/liblnagLF.a
#   LIBSM = /usr/local/lib/libplotsub.a \
#         /usr/local/lib/libdevices.a \
#         /usr/local/lib/libutils.a -L/usr/lib -lX11
  UnixOrWin := "unix"
endif

########################################

########################################
## For g77 fortran
ifeq ($(SYSTYPE),g77)
  FC = g77  # g77 fortran
  FFLAGS = -c -g -fbounds-check -finit-local-zero -fno-automatic
  FFLAGS = -c -g -finit-local-zero -fno-automatic
  FFLAGS = -c -O2 -finit-local-zero -fno-automatic -Wall
  LDFLAGS =
  LIBS= $(HOMEStella)/lib/libsparse77.a $(HOMEStella)/lib/liblnag77.a
  LIBSM = /usr/local/lib/libplotsub.a \
        /usr/local/lib/libdevices.a \
        /usr/local/lib/libutils.a -L/usr/lib -lX11
  UnixOrWin := "unix"
endif
########################################


########################################
## For gfortran fortran
ifeq ($(SYSTYPE),gfortran)
  FC = gfortran
  FFLAGS_FIX = -c -g -fbounds-check
  FFLAGS_FIX = -c -g  
  FFLAGS_FIX = -c -O0  -ffixed-line-length-132  -Wall
  FFLAGS_FIX = -c -O2  -ffixed-line-length-132  -w
#  FFLAGS_FIX = -c -O2 -g -fbounds-check  -ffixed-line-length-132 -fno-automatic -Wall
  FFLAGS := -ffree-form
  FFLAGS := $(FFLAGS_FIX) $(FFLAGS)
  FFLAGS_FIX += -fno-automatic
  LDFLAGS =
  #LIBS= $(HOME)/lib/libsparseGF.a $(HOME)/lib/liblnagGF.a
  #LIBS= $(HOMEStella)/lib/libsparseGF.a $(HOMEStella)/lib/liblnagGF.a
  LIBS= 
  LIBSM = /usr/local/lib/libplotsub.a \
        /usr/local/lib/libdevices.a \
        /usr/local/lib/libutils.a -L/usr/lib -lX11
  UnixOrWin := "unix"
endif
########################################



########################################
## For Windows
ifeq ($(UnixOrWin),"win")
 FFLAGS += -object:$@
 LDFLAGS += /exe:
 slash_l := \\
endif
########################################
## For Linux
ifeq ($(UnixOrWin),"unix")
 LDFLAGS += -o
 slash_l := /
endif
#-------------------------------------------------------
#---
#-------------------------------------------------------

DEL = rm -f

STLHOME  := ../
DIR_SRC := $(STLHOME)

vpath %.trf ../eve:../strad:../vladsf:../src
vpath %.f ../eve:../strad:../vladsf:../src
vpath %.f90 ../src
vpath %.o ./
vpath %.inc ../src


VPATH := $(STLHOME)src $(STLHOME)src/cp2k $(STLHOME)src/stl $(STLHOME)src/util $(STLHOME)vladsf
INCL_DIR := -I$(STLHOME)src/ -I$(STLHOME)src/cp2k -I$(STLHOME)src/stl -I$(STLHOME)src/util -I$(STLHOME)vladsf

#--------------------------
FILES_base =  kinds.F90  math_constants.f90  phys_constants.f90
TEMP := $(FILES_base:.trf=.F)
TEMP := $(patsubst %.F90,%.f90, $(TEMP:.f=.f90))
TEMP := $(patsubst %.F90,%.f90, $(TEMP:.F=.f90))
OBJS_base := $(patsubst %.f90,%.o, $(TEMP))

FILES_cross =  rad_photo_cross_section.f90
TEMP := $(FILES_cross:.trf=.F)
TEMP := $(patsubst %.F90,%.f90, $(TEMP:.F=.f90))
TEMP := $(TEMP:.f=.o)
OBJScross := $(patsubst %.f90,%.o, $(TEMP))

#--------------------------


INCL = opacity.inc stateq.inc sahaandd.inc\
       abo.inc nstep.inc stsave.inc fundrad.inc \
       black.inc snrad.inc formod.inc \
       zone.inc commonEve.inc azzn.inc

INCLOUTUNI = opacityOutUni.inc stateq.inc sahaandd.inc \
       abo.inc nstep.inc stsave.inc fundrad.inc \
       black.inc snrad.inc formod.inc \
       zone.inc commonEve.inc azzn.inc

OBJDIR = obj

PROGEVE1A = ../eve/run/eve1a.exe

PROGEVEPAB = ../eve/run/evenewpab.exe
PROGEVEPAB := $(subst /,$(slash_l),$(PROGEVEPAB))


PROGEVE = ../eve/run/evenew.exe
PROGEVE := $(subst /,$(slash_l),$(PROGEVE))

PROGEVENT = ../eve/run/event.exe
PROGEVENT := $(subst /,$(slash_l),$(PROGEVENT))

PROGEVE1A := $(subst /,$(slash_l),$(PROGEVE1A))

PROGEVE2 = ../eve/run/eve2.exe
PROGEVE2     := $(subst /,$(slash_l),$(PROGEVE2))

PROGEVEW7 = ../eve/run/evew7.exe
PROGEVEW7     := $(subst /,$(slash_l),$(PROGEVEW7))

PROGEVE2NONI = ../eve/run/eve2.exe
PROGEVE2NONI     := $(subst /,$(slash_l),$(PROGEVE2NONI))

PROGEVE2c = ../eve/run/eve2c.exe
PROGEVE2c     := $(subst /,$(slash_l),$(PROGEVE2c))

PROGEVEWS = ../eve/run/evewindsi.exe
PROGEVEWS := $(subst /,$(slash_l),$(PROGEVEWS))

PROGSTNR = ../strad/run/stella.exe # non-relativ. as royac1

PROGSTR = ../strad/run/stellarel.exe
PROGSTR := $(subst /,$(slash_l),$(PROGSTR))

PROGSTNR5 = ../strad/run/xstella5.exe # non-relativ. with gdepos5 & volennoint
PROGSTNR5   := $(subst /,$(slash_l),$(PROGSTNR5))

PROGSTNR6 = ../strad/run/xstella6.exe # non-relativ. with gdepos6 & volennoint
PROGSTNR6   := $(subst /,$(slash_l),$(PROGSTNR6))

PROGSTNR6_MESA = ../strad/run/xstella6_mesa.exe # non-relativ. with gdepos6 & volennoint
PROGSTNR6_MESA   := $(subst /,$(slash_l),$(PROGSTNR6_MESA))

PROGSTNR6Y12M = ../strad/run/xstella6y12m.exe # non-relativ. with gdepos6 & volennoint
PROGSTNR6Y12M   := $(subst /,$(slash_l),$(PROGSTNR6Y12M))

PROGSTNR6magn = ../strad/run/xstella6magn.exe # non-relativ. with gdepos6 & volennointmag for magnetar, C Takashi Moriya
PROGSTNR6magn   := $(subst /,$(slash_l),$(PROGSTNR6magn))

PROGSTCON = ../strad/run/xstellacon.exe # non-relativ. with constant opacity
PROGSTCON   := $(subst /,$(slash_l),$(PROGSTCON))

PROGSTNR6WC = ../strad/run/xstellaWC.exe # non-relativ. with gdepos6 & volennoint and compton
PROGSTNR6WC   := $(subst /,$(slash_l),$(PROGSTNR6WC))

PROGSTNR6WEAVER = ../strad/run/xstellaWCWeaver.exe # non-relativ. with gdepos6 & volennoint and compton for weaver
PROGSTNR6WEAVER   := $(subst /,$(slash_l),$(PROGSTNR6WEAVER))

PROGSTNR6WEAVNC = ../strad/run/xstellaNCWeaver.exe # non-relativ. with gdepos6 & volennoint and NO compton for weaver
PROGSTNR6WEAVNC   := $(subst /,$(slash_l),$(PROGSTNR6WEAVNC))

PROGSTFEAU = ../strad/run/xstellaFeau.exe # non-relativ. with gdepos6 & volennoint, feau
PROGSTFEAU   := $(subst /,$(slash_l),$(PROGSTFEAU))

PROGSTFEAOU = ../strad/run/xstellaFeaou.exe # non-relativ. with gdepos6 & volennoint, feau Out Uniform
PROGSTFEAOU   := $(subst /,$(slash_l),$(PROGSTFEAOU))

PROGSTLFRQ = ../strad/run/xstellaLfrq.exe # non-relativ. with gdepos6 & volennoint
PROGSTLFRQ   := $(subst /,$(slash_l),$(PROGSTLFRQ))

PROGSTRZN = ../strad/run/xstellarzn.exe # non-relativ. with rezone, gdepos6 & volennoint
PROGSTRZN   := $(subst /,$(slash_l),$(PROGSTRZN))

PROGSTOUTUNI = ../strad/run/xstellaOutUniFH0.exe # non-relativ. with gdepos6 & volennoint Out Uniform
PROGSTOUTUNI  := $(subst /,$(slash_l),$(PROGSTOUTUNI))

PROGSTNRTST = ../strad/run/xstellatst.exe # as NR6 with test EddN=1/3
PROGSTNRTST   := $(subst /,$(slash_l),$(PROGSTNRTST))

PROGSTArnt = ../strad/run/xstellaArnett.exe #  relativ. with gdepos6 & volennoint
PROGSTArnt   := $(subst /,$(slash_l),$(PROGSTArnt))

PROGSTNRRADA = ../strad/run/stellarada.exe # integration with radiation transport
PROGSTNRRADA := $(subst /,$(slash_l),$(PROGSTNRRADA))

PROGSTRON = ../strad/run/stellaron.exe # non-relativ. direct Ron opacity
PROGSTRON := $(subst /,$(slash_l),$(PROGSTRON))

PROGSTRONALL = ../strad/run/stellaronall.exe # non-relativ. direct Ron opacity
PROGSTRONALL := $(subst /,$(slash_l),$(PROGSTRONALL))

PROGSTNR7BQEPS = ../strad/run/xstella7BQeps.exe # with stiffBGHeps no interp for TOO non-relativ. with gdepos6 & volennoint
PROGSTNR7BQEPS   := $(subst /,$(slash_l),$(PROGSTNR7BQEPS))


PROGTT = ../strad/run/xttfit.exe
PROGTT := $(subst /,$(slash_l),$(PROGTT))

PROGTTOU = ../strad/run/xttfitou.exe
PROGTTOU := $(subst /,$(slash_l),$(PROGTTOU))

PROGTTS = ../strad/run/xttfits.exe
PROGTTS := $(subst /,$(slash_l),$(PROGTTS))

PROGTTSE = ../strad/run/ettfits.exe
PROGTTSE := $(subst /,$(slash_l),$(PROGTTSE))

PROGTTPS = ../strad/run/xttfitps.exe
PROGTTPS := $(subst /,$(slash_l),$(PROGTTPS))

PROGTTNR = ../strad/run/xttfitnr.exe
PROGTTNR := $(subst /,$(slash_l),$(PROGTTNR))

PROGTTR = ../strad/run/xttfitrel.exe
PROGTTR := $(subst /,$(slash_l),$(PROGTTR))

PROGRONFICTD = ../vladsf/xronfictd.exe
PROGRONFICTD := $(subst /,$(slash_l),$(PROGRONFICTD))

#PROGRONFICT = ../run/vladsf/xvladsf.exe
PROGRONFICT = ../vladsf/xronfict.exe
PROGRONFICT := $(subst /,$(slash_l),$(PROGRONFICT))

# PROGRONFTF = ../run/vladsf/xronftf.exe
PROGRONFTF = ../vladsf/xronftf.exe
PROGRONFTF := $(subst /,$(slash_l),$(PROGRONFTF))


PROGRONLFRQ = ../vladsf/xronLfrq.exe
PROGRONLFRQ := $(subst /,$(slash_l),$(PROGRONLFRQ))

PROGRONFICTOUTUNI = ../vladsf/xronfictOutUni.exe
PROGRONFICTOUTUNI := $(subst /,$(slash_l),$(PROGRONFICTOUTUNI))

PROGINSH = ../vladsf/xinsh.exe
PROGINSH:= $(subst /,$(slash_l),$(PROGINSH))

PROGINSHSHB = ../vladsf/xinshb.exe
PROGINSHSHB:= $(subst /,$(slash_l),$(PROGINSHSHB))

PROGBFTOT = ../vladsf/xronfbftot.exe
PROGBFTOT:= $(subst /,$(slash_l),$(PROGBFTOT))

PROGBFTOTSHB = ../vladsf/xronfbftotshb.exe
PROGBFTOTSHB:= $(subst /,$(slash_l),$(PROGBFTOTSHB))

PROGRONFSHBVR = ../vladsf/xronfshbvR.exe
PROGRONFSHBVR := $(subst /,$(slash_l),$(PROGRONFSHBVR))

PROGRONFSHB = ../vladsf/xronfshb.exe
PROGRONFSHB := $(subst /,$(slash_l),$(PROGRONFSHB))

PROGOPATOLD = ../vladsf/opatold.exe
PROGOPATOLD := $(subst /,$(slash_l),$(PROGOPATOLD))

PROGOPATINSH = ../vladsf/opatinsh.exe
PROGOPATINSH := $(subst /,$(slash_l),$(PROGOPATINSH))

PROGOPATBFTOT = ../vladsf/opatbftot.exe
PROGOPATBFTOT := $(subst /,$(slash_l),$(PROGOPATBFTOT))

PROGTTR = ../strad/run/xttfit.exe
PROGTTR := $(subst /,$(slash_l),$(PROGTTR))

FILES_EVEPAB = evenewpab.trf evebeg.trf bgconhepb.trf resdia.trf evefun.trf \
         evesub.trf forback.trf sahaz.trf \
         sahaandd.trf ubv.trf volenpum.trf opazr.trf stradio.trf \
         azdat.trf
OBJEVEPAB := $(FILES_EVEPAB:.trf=.o)

FILES_EVE =evenew.trf  evebeg.trf bgconhe.trf resdia.trf evefun.trf  \
         evesub.trf forback.trf sahaz.trf \
         sahaandd.trf ubv.trf volenpum.trf opazr.trf stradio.trf \
         azdat.trf
OBJEVE := $(FILES_EVE:.trf=.o)

FILES_EVENT =event.trf  evebeg.trf bgconnt.trf resdia.trf evefun.trf  \
         evesub.trf forback.trf sahaz.trf \
         sahaandd.trf ubv.trf volenpum.trf opazr.trf stradio.trf \
         azdat.trf
OBJEVENT := $(FILES_EVENT:.trf=.o)

FILES_EVE1A = evesn1a.trf sahaz.trf sahaandd.trf opazr.trf stradio.trf \
         azdat.trf length.trf
OBJEVE1A := $(FILES_EVE1A:.trf=.o)

FILES_EVE2 = evesn2.trf sahaz.trf sahaandd.trf opazr.trf stradio.trf \
         azdat.trf length.trf
OBJEVE2 := $(FILES_EVE2:.trf=.o)

FILES_EVEW7 = evesn2w7.trf sahaz.trf sahaandd.trf opazr.trf stradio.trf \
         azdat.trf length.trf
OBJEVEW7 := $(FILES_EVEW7:.trf=.o)

FILES_EVE2NONI =evesn2NoNi.trf sahaz.trf sahaandd.trf opazr.trf stradio.trf \
         azdat.trf length.trf
OBJEVE2NONI := $(FILES_EVE2NONI:.trf=.o)

FILES_EVEWS = evewindSi.trf evebeg.trf bgconhepb.trf resdia.trf evefun.trf \
         evesub.trf forback.trf sahaz.trf \
         sahaandd.trf ubv.trf volenpum.trf opazr.trf stradio.trf \
         azdat.trf
OBJEVEWS := $(FILES_EVEWS:.trf=.o)

FILES_STOUTUNI = stradOutUni.trf begradOutUni.trf cosetbgh.trf hcdfjrad.trf \
        hcdfnrad.trf traneq.trf eddi.trf  gdepos6.trf nthnew.trf \
        stiffbgh.trf lbalsw.trf stradio.trf \
        vtimef90.trf sahaandd.trf ubv.trf obsubvri.trf \
        tt4strad.trf begttOutUni.trf lbol.trf \
        burnc.trf volenpumnoint.trf hapOutUni.trf hcdhaph.trf \
        oparonOutUni.trf length.trf words.trf azdat.trf
OBJSTOUTUNI := $(FILES_STOUTUNI:.trf=.o)


FILES_RONFD = ronfndecd.trf bessi0.f bessk0ex.f \
            blas.f dmach.trf edensol.trf \
            gffcalc.f gshfdxsec.f hydxsecl.f \
            hydxsecn.f hypho.f lineexpab_cor.trf \
            lnblnk.f ndex.f \
            pfsaha.f sahaeqn.f setnucms.f \
            sparseblas.f tablsort.f valence_nl.f \
            azdat.trf length.trf stradio.trf \
            vtimef90.trf hMinusAbsorp.trf opacitytd.trf
OBJS_RONFD := $(OBJS_base) $(OBJScross) $(patsubst %.f,%.o, $(FILES_RONFD:.trf=.f))

FILES_RONF = ronfndec.trf bessi0.f bessk0ex.f \
            blas.f dmach.trf edensol.trf \
            gffcalc.f gshfdxsec.f hydxsecl.f \
            hydxsecn.f hypho.f lineexpab_cor.trf \
            lnblnk.f ndex.f \
            pfsaha.f sahaeqn.f setnucms.f \
            sparseblas.f tablsort.f valence_nl.f \
            azdat.trf length.trf stradio.trf \
            vtimef90.trf hMinusAbsorp.trf opacityt.trf
OBJS_RONF := $(OBJS_base) $(OBJScross) $(patsubst %.f,%.o, $(FILES_RONF:.trf=.f))


FILES_RONFTF = ronftfixed.f bessi0.f bessk0ex.f \
            blas.f dmach.trf edensol.trf \
            gffcalc.f gshfdxsec.f hydxsecl.f \
            hydxsecn.f hypho.f lineexpab_cor.trf \
            lnblnk.f ndex.f opacityt.trf \
            pfsaha.f sahaeqn.f setnucms.f \
            sparseblas.f tablsort.f valence_nl.f \
            azdat.trf length.trf stradio.trf \
            vtimef90.trf hMinusAbsorp.trf
OBJS_RONFTF := $(patsubst %.f,%.o, $(FILES_RONFTF:.trf=.f)) \
        $(OBJS_base) $(OBJScross)

FILES_RONLFRQ = ronfndecLFRQ.trf bessi0.f bessk0ex.f \
            blas.f dmach.trf edensol.trf \
            gffcalc.f gshfdxsec.f hydxsecl.f \
            hydxsecn.f hypho.f lineexpab_cor.trf \
            lnblnk.f ndex.f opacity.trf \
            pfsaha.f sahaeqn.f setnucms.f \
            sparseblas.f tablsort.f valence_nl.f \
            azdat.trf length.trf stradio.trf \
            vtimef90.trf hMinusAbsorp.trf
OBJS_RONLFRQ := $(patsubst %.f,%.o, $(FILES_RONLFRQ:.trf=.f)) \
        $(OBJS_base) $(OBJScross)

FILES_RONFOUTUNI = ronfndecOutUni.trf bessi0.f bessk0ex.f \
            blas.f dmach.trf edensol.trf \
            gffcalc.f gshfdxsec.f hydxsecl.f \
            hydxsecn.f hypho.f lineexpab_cor.trf \
            lnblnk.f ndex.f opacityt.trf \
            pfsaha.f sahaeqn.f setnucms.f \
            sparseblas.f tablsort.f valence_nl.f \
            azdat.trf length.trf stradio.trf \
            vtimef90.trf hMinusAbsorp.trf
OBJS_RONFOUTUNI := $(patsubst %.f,%.o, $(FILES_RONFOUTUNI:.trf=.f)) \
        $(OBJS_base) $(OBJScross)


FILES_RONFSHB = ronfshb.trf bessi0.f bessk0ex.f \
            blas.f dmach.trf edensol.trf \
            gffcalc.f gshfdxsec.f hydxsecl.f \
            hydxsecn.f hypho.f lineexpab_cor.trf \
            lnblnk.f ndex.f opacityt.trf \
            pfsaha.f sahaeqn.f setnucms.f \
            sparseblas.f tablsort.f valence_nl.f \
            azdat.trf length.trf stradio.trf \
            vtimef90.trf hMinusAbsorp.trf
OBJS_RONFSHB := $(patsubst %.f,%.o, $(FILES_RONFSHB:.trf=.f)) \
        $(OBJS_base) $(OBJScross)

FILES_RONFSHBVR = ronfshb.trf bessi0.f bessk0ex.f \
            blas.f dmach.trf edensol.trf \
            gffcalc.f gshfdxsec.f hydxsecl.f \
            hydxsecn.f hypho.f lineexpVRegem.trf \
            lnblnk.f ndex.f opacityt.trf \
            pfsaha.f sahaeqn.f setnucms.f \
            sparseblas.f tablsort.f valence_nl.f \
            azdat.trf length.trf stradio.trf \
            vtimef90.trf hMinusAbsorp.trf
OBJS_RONFSHBVR := $(patsubst %.f,%.o, $(FILES_RONFSHBVR:.trf=.f)) \
        $(OBJS_base) $(OBJScross)

FILES_INSH = ronfndec.trf azdat.trf bessi0.f bessk0ex.f \
        blas.f dmach.trf edensol.trf  gffcalc.f gshfdxsec.f \
        hydxsecl.f hydxsecn.f hypho.f length.trf lineexpab_cor.trf \
        lnblnk.f ndex.f opacityInSh.trf pfsaha.f \
        sahaeqn.f setnucms.f sparseblas.f tablsort.f stradio.trf \
        valence_nl.f vtimef90.trf hMinusAbsorp.trf
OBJS_INSH := $(patsubst %.f,%.o, $(FILES_INSH:.trf=.f)) \
        $(OBJS_base) $(OBJScross)

FILES_INSHSHB = ronfshb.trf azdat.trf bessi0.f bessk0ex.f \
        blas.f dmach.trf edensol.trf  gffcalc.f gshfdxsec.f \
        hydxsecl.f hydxsecn.f hypho.f length.trf lineexpab_cor.trf \
        lnblnk.f ndex.f opacityInSh.trf pfsaha.f \
        sahaeqn.f setnucms.f sparseblas.f tablsort.f stradio.trf \
        valence_nl.f vtimef90.trf hMinusAbsorp.trf
OBJS_INSHSHB := $(patsubst %.f,%.o, $(FILES_INSHSHB:.trf=.f)) \
        $(OBJS_base) $(OBJScross)

# OBJS_INSHSHB = ronfshb.o azdat.o bessi0.o bessk0ex.o \
#         blas.o dmach.o edensol.o gffcalc.o gshfdxsec.o \
#         hydxsecl.o hydxsecn.o hypho.o length.o lineexpab_cor.o \
#         lnblnk.o ndex.o opacityInSh.o pfsaha.o \
#         sahaeqn.o setnucms.o sparseblas.o tablsort.o  stradio.o \
#         valence_nl.o vtimef90.o  hMinusAbsorp.o\
#         $(OBJS_base) $(OBJScross)

FILES_BFTOT = ronfndecbftot.trf bessi0.f bessk0ex.f \
            blas.f dmach.trf edensol.trf \
            gffcalc.f gshfdxsec.f hydxsecl.f \
            hydxsecn.f hypho.f lineexpab_cor.trf \
            lnblnk.f ndex.f opacitybftot.trf \
            pfsaha.f sahaeqn.f setnucms.f \
            sparseblas.f tablsort.f valence_nl.f \
            azdat.trf length.trf stradio.trf \
            vtimef90.trf hMinusAbsorp.trf
OBJS_BFTOT := $(patsubst %.f,%.o, $(FILES_BFTOT:.trf=.f)) \
         $(OBJS_base) $(OBJScross)

FILES_BFTOTSHB = ronfshbbftot.trf bessi0.f bessk0ex.f \
            blas.f dmach.trf edensol.trf \
            gffcalc.f gshfdxsec.f hydxsecl.f \
            hydxsecn.f hypho.f lineexpab_cor.trf \
            lnblnk.f ndex.f opacitybftot.trf \
            pfsaha.f sahaeqn.f setnucms.f \
            sparseblas.f tablsort.f valence_nl.f \
            azdat.trf length.trf stradio.trf \
            vtimef90.trf hMinusAbsorp.trf
OBJS_BFTOTSHB := $(patsubst %.f,%.o, $(FILES_BFTOTSHB:.trf=.f)) \
         $(OBJS_base) $(OBJScross)

FILES_TNR6 = stradsep5tt.trf begradsep.trf cosetbgh.trf hcdfjrad.trf \
        hcdfnrad.trf traneq.trf eddi.trf  gdepos6.trf nthnew.trf \
        stiffbgh.trf lbalsw.trf stradio.trf \
        vtimef90.trf sahaandd.trf ubv.trf obsubvri.trf \
        tt4strad.trf begtt.trf lbol.trf \
        burnc.trf volenpumnoint.trf hapsepnc.trf hcdhaph.trf \
        oparon.trf length.trf words.trf azdat.trf
OBJSTNR6 := $(FILES_TNR6:.trf=.o)

FILES_TNR6_MESA = stradsep5tt.trf begradsep.trf cosetbgh.trf hcdfjrad.trf \
        hcdfnrad.trf traneq.trf eddi.trf  gdepos6.trf nthnew.trf \
        stiffbgh.trf lbalsw.trf stradio.trf \
        vtimef90.trf sahaandd.trf ubv.trf obsubvri.trf \
        tt4strad.trf begtt.trf lbol.trf \
        burnc.trf volenpumnoint.trf hapsepnc.trf hcdhaph.trf \
        oparon.trf length.trf words.trf azdat.trf \
        f01brqt.f f01brrt.f f01brst.f f01brtt.f \
        f01brut.f f01brvt.f f01brwt.f f01brxt.f f01bryt.f f01brzt.f \
        p01abf.f p01abz.f x04aaf.f x04abf.f x04baf.f \
        m28aux.f m28y12.f
OBJSTNR6_MESA := $(patsubst %.f,%.o, $(FILES_TNR6_MESA:.trf=.f))
# m28aux.f contains m30y12.f && m28bys.f && m28cys.f && m28cyn.f

FILES_TNR6Y12M = stradsep5tt.trf begradsep.trf cosetbgh.trf hcdfjrad.trf \
        hcdfnrad.trf traneq.trf eddi.trf  gdepos6.trf nthnew.trf \
        stiffbghY12m.trf lbalsw.trf stradio.trf \
        vtimef90.trf sahaandd.trf ubv.trf obsubvri.trf \
        tt4strad.trf begtt.trf lbol.trf \
        burnc.trf volenpumnoint.trf hapsepnc.trf hcdhaph.trf \
        oparon.trf length.trf words.trf azdat.trf y12m.f
OBJSTNR6Y12M := $(patsubst %.f,%.o, $(FILES_TNR6Y12M:.trf=.f))

FILES_TNR6magn = stradsep5tt.trf begradsep.trf cosetbgh.trf hcdfjrad.trf \
        hcdfnrad.trf traneq.trf eddi.trf  gdepos6.trf nthnew.trf \
        stiffbgh.trf lbalsw.trf stradio.trf \
        vtimef90.trf sahaandd.trf ubv.trf obsubvri.trf \
        tt4strad.trf begtt.trf lbol.trf \
        burnc.trf volenpumnointmag.trf hapsepnc.trf hcdhaph.trf \
        oparon.trf length.trf words.trf azdat.trf
OBJSTNR6magn := $(FILES_TNR6magn:.trf=.o)

#hapsepnc.trf No Compton happa !

FILES_STCON = stradsep5tt.trf begradsep.trf cosetbgh.trf hcdfjrad.trf \
        hcdfnrad.trf traneq.trf eddi.trf  gdepos6.trf nthnew.trf \
        stiffbgh.trf lbalsw.trf stradio.trf findTradByMNKLum.trf \
        vtimef90.trf sahaandd.trf ubv.trf obsubvri.trf \
        tt4strad.trf begtt.trf lbol.trf \
        burnc.trf volenpumnoint.trf hapconst.trf hcdhaph.trf \
        opaconst.trf length.trf words.trf azdat.trf
OBJSTCON := $(FILES_STCON:.trf=.o)

FILES_OPATOLD = opatest.trf azdat.trf bessi0.f bessk0ex.f \
        blas.f dmach.trf edensol.trf gffcalc.f gshfdxsec.f \
        hydxsecl.f hydxsecn.f hypho.f lineexpab_cor.trf \
        lnblnk.f ndex.f opacityt.trf pfsaha.f \
        sahaeqn.f setnucms.f sparseblas.f tablsort.f \
        valence_nl.f vtimef90.trf  hMinusAbsorp.trf
OBJS_OPATOLD := $(patsubst %.f,%.o, $(FILES_OPATOLD:.trf=.f)) \
        $(OBJS_base) $(OBJScross)

FILES_OPATINSH = opatest.trf azdat.trf bessi0.f bessk0ex.f \
        blas.f dmach.trf edensol.trf  gffcalc.f gshfdxsec.f \
        hydxsecl.f hydxsecn.f hypho.f lineexpab_cor.trf \
        lnblnk.f ndex.f opacityInSh.trf pfsaha.f \
        sahaeqn.f setnucms.f sparseblas.f tablsort.f stradio.trf \
        valence_nl.f vtimef90.trf hMinusAbsorp.trf
OBJS_OPATINSH := $(patsubst %.f,%.o, $(FILES_OPATINSH:.trf=.f)) \
        $(OBJS_base) $(OBJScross)

FILES_OPATBFTOT = opatestbftot.trf azdat.trf bessi0.f bessk0ex.f \
        blas.f dmach.trf edensol.trf  gffcalc.f gshfdxsec.f \
        hydxsecl.f hydxsecn.f hypho.f lineexpab_cor.trf \
        lnblnk.f ndex.f opacitybftot.trf pfsaha.f \
        sahaeqn.f setnucms.f sparseblas.f tablsort.f stradio.trf \
        valence_nl.f vtimef90.trf hMinusAbsorp.trf
OBJS_OPATBFTOT := $(patsubst %.f,%.o, $(FILES_OPATBFTOT:.trf=.f)) \
        $(OBJS_base) $(OBJScross)

FILES_STR =  stradsep5ttrel.trf begradsep.trf cosetbgh.trf dfnradrel.trf \
        dfjradrel.trf traneq.trf eddirel.trf gdepos6.trf nthnew.trf \
        stiffbgh.trf lbalswrel.trf stradio.trf findTradByMNKLum.trf \
        vtimef90.trf sahaandd.trf ubv.trf obsubvri.trf \
        tt4strad.trf begtt.trf lbol.trf \
        burnc.trf volenpumnoint.trf hapsepnc.trf hcdhaph.trf \
        oparon.trf length.trf words.trf azdat.trf
OBJSTR := $(FILES_STR:.trf=.o)

FILES_TNR6WC = stradsep5tt.trf begradsep.trf cosetbgh.trf hcdfjrad.trf \
        hcdfnrad.trf traneq.trf eddi.trf  gdepos6.trf nthnew.trf \
        stiffbgh.trf lbalsw.trf stradio.trf \
        vtimef90.trf sahaandd.trf ubv.trf obsubvri.trf \
        tt4strad.trf begtt.trf lbol.trf \
        burnc.trf volenpumnoint.trf hapsepWithComp.trf hcdhaph.trf \
        oparon.trf length.trf words.trf azdat.trf
OBJSTNR6WC := $(FILES_TNR6WC:.trf=.o)

FILES_TNR6WEAVER = stradweaver.trf begradsep.trf cosetbgh.trf hcdfjrad.trf \
        hcdfnrad.trf traneq.trf eddi.trf  gdepos6.trf nthnew.trf \
        stiffbgh.trf lbalsw.trf stradio.trf \
        vtimef90.trf sahaandd.trf ubv.trf obsubvri.trf \
        tt4strad.trf begtt.trf lbol.trf \
        burnc.trf volenpumnoint.trf hapsepWithComp.trf hcdhaph.trf \
        oparon.trf length.trf words.trf azdat.trf
OBJSTNR6WEAVER := $(FILES_TNR6WEAVER:.trf=.o)

FILES_TNR6WEAVNC = stradweaver.trf begradsep.trf cosetbgh.trf hcdfjrad.trf \
        hcdfnrad.trf traneq.trf eddi.trf  gdepos6.trf nthnew.trf \
        stiffbgh.trf lbalsw.trf stradio.trf \
        vtimef90.trf sahaandd.trf ubv.trf obsubvri.trf \
        tt4strad.trf begtt.trf lbol.trf \
        burnc.trf volenpumnoint.trf hapsepnc.trf hcdhaph.trf \
        oparon.trf length.trf words.trf azdat.trf
OBJSTNR6WEAVNC := $(FILES_TNR6WEAVNC:.trf=.o)

FILES_FEAU = stradsep5tt.trf begradsep.trf cosetbgh.trf hcdfjrad.trf \
        hcdfnrad.trf traneq.trf obfeaw.trf  gdepos6.trf nthnew.trf \
        stiffbgh.trf lbalsw.trf stradio.trf \
        vtimef90.trf sahaandd.trf ubv.trf obsubvri.trf \
        tt4strad.trf begtt.trf lbol.trf \
        burnc.trf volenpumnoint.trf hapsepnc.trf hcdhaph.trf \
        oparon.trf length.trf words.trf azdat.trf
OBJSTFEAU := $(FILES_FEAU:.trf=.o)

FILES_FEAOU = stradOutUni.trf begradOutUni.trf cosetbgh.trf hcdfjrad.trf \
        hcdfnrad.trf traneq.trf obfeaw.trf  gdepos6.trf nthnew.trf \
        stiffbgh.trf lbalsw.trf stradio.trf \
        vtime.trf sahaandd.trf ubv.trf obsubvri.trf \
        tt4strad.trf begttOutUni.trf lbol.trf \
        burnc.trf volenpumnoint.trf hapOutUni.trf hcdhaph.trf \
        oparonOutUni.trf length.trf words.trf azdat.trf
OBJSTFEAOU := $(FILES_FEAOU:.trf=.o)

#hapsepnc.trf No Compton happa !

FILES_LFRQ = stradsep5tt.trf begradLowFreq.trf cosetbgh.trf hcdfjrad.trf \
        hcdfnrad.trf traneq.trf eddi.trf  gdepos6.trf nthnew.trf \
        stiffbgh.trf lbalsw.trf stradio.trf \
        vtimef90.trf sahaandd.trf ubv.trf obsubvri.trf \
        tt4strad.trf begtt.trf lbol.trf \
        burnc.trf volenpumnoint.trf hapsepnc.trf hcdhaph.trf \
        oparon.trf length.trf words.trf azdat.trf
OBJSTLFRQ := $(FILES_LFRQ:.trf=.o)

FILES_STNR7BQEPS = stradsep7ttBQ.trf begradsep.trf cosetbgh.trf hcdfjrad.trf \
        hcdfnrad.trf traneq.trf eddi.trf  gdepos6.trf nthnew.trf \
        stiffbghEps.trf lbalswBUG.trf stradio.trf findTradByMNKLum.trf \
        vtimef90.trf sahaandd.trf ubv.trf obsubvri.trf \
        tt4strad.trf begtt.trf lbol.trf \
        burnc.trf volenpumnoint.trf hapsepnc.trf hcdhaph.trf \
        oparon.trf length.trf words.trf azdat.trf
OBJSTNR7BQEPS := $(FILES_STNR7BQEPS:.trf=.o)

FILES_TTFIT = ttfit5.trf begradsep.trf vtimef90.trf \
        nthnew.trf  lbalsw.trf    words.trf \
        stradio.trf  sahaandd.trf  ubv.trf \
        obsubvri.trf  burnc.trf  volenpumnoint.trf \
        hapsepnc.trf  oparon.trf  length.trf \
        azdat.trf lbol.trf
OBJTT := $(FILES_TTFIT:.trf=.o)

FILES_TTFITOU = ttfit5OutUni.trf begradOutUni.trf vtimef90.trf \
        nthnew.trf  lbalsw.trf    words.trf \
        stradio.trf  sahaandd.trf  ubv.trf \
        obsubvri.trf  burnc.trf  volenpumnoint.trf \
        hapsepnc.trf  oparon.trf  length.trf \
        azdat.trf lbol.trf
OBJTTOU := $(FILES_TTFITOU:.trf=.o)

FILES_TTFITS = ttfitsimple.trf vecsubs.trf begradsep.trf eddi.trf vtimef90.trf \
               nthnew.trf  lbalsw.trf  lbol.trf   words.trf  \
               stradio.trf    sahaandd.trf  ubv.trf \
               obsubvri.trf   burnc.trf  volenpumnoint.trf \
               hapsepnc.trf  oparon.trf  length.trf \
               azdat.trf
OBJTTS := $(FILES_TTFITS:.trf=.o)

FILES_TTFITSE = ttfitsimple.trf vecsubs.trf begradsep.trf eddiplot.trf vtimef90.trf \
        nthnew.trf  lbalsw.trf  lbol.trf   words.trf \
        stradio.trf    sahaandd.trf  ubv.trf \
        obsubvri.trf   burnc.trf  volenpumnoint.trf \
        hapsepnc.trf  oparon.trf  length.trf \
        azdat.trf
OBJTTSE := $(FILES_TTFITSE:.trf=.o)


FILES_TTFITPS = ttfitps.trf vecsubs.trf begradsep.trf eddiplot.trf vtimef90.trf \
        nthnew.trf  lbalsw.trf  lbol.trf   words.trf \
        stradio.trf    sahaandd.trf  ubv.trf \
        obsubvri.trf   burnc.trf  volenpumnoint.trf \
        hapsepnc.trf  oparon.trf  length.trf \
        azdat.trf ps_write.f90
OBJTTPS := $(patsubst %.f90,%.o, $(FILES_TTFITPS:.trf=.o))



TREFOR = trf

# for ctrf:
TRFFLAGS = -nfs



.IGNORE:

# use cleanf to remove the .f files before doing make for trf's
#%.o : %.trf
#	$(TREFOR) $(TRFFLAGS) $<
#	$(FC) $(FFLAGS_FIX)  $(patsubst %.trf,%.f,$<) $(INCL_DIR)

# this uses existing .f files, bypassing trf
%.o : %.trf
	$(FC) $(FFLAGS_FIX)  $(patsubst %.trf,%.f,$<) $(INCL_DIR)



%.o: %.f
	$(FC) $(FFLAGS_FIX) $< $(INCL_DIR)

%.o: %.f90
	$(FC) $(FFLAGS) $< $(INCL_DIR)
%.o: %.F90
	$(FC) $(FFLAGS)  $< $(INCL_DIR)
%.o: %.F
	$(FC) $(FFLAGS)  $<  $(INCL_DIR)

%.mod : %.o
	@if [! -f $@ ]; then \
	rm $< \
	$(MAKE) $< \
	fi


# %.mod :
# 	@if [! -f $@ ]; then \
#           rm $(*F).o; \
#         fi
# 	$(MAKE) $<
#
# %.o : %.f90
# 	$(FC) -c -o $(*F).o $<
#
# %.o : %.F90
# 	$(FC) -c -o $(*F).o $<
#
# %.o : %.F
# 	$(FC) -c -o $(*F).o $<
#

all: help

evenewpab: $(OBJEVEPAB)
	$(FC) $(LDFLAGS) $(PROGEVEPAB) $(OBJEVEPAB) $(LIBS) $(LIBDIRSYS)


evenew: $(OBJEVE)
	$(FC) $(LDFLAGS) $(PROGEVE) $(OBJEVE) $(LIBS) $(LIBDIRSYS)

event: $(OBJEVENT)
	$(FC) $(LDFLAGS) $(PROGEVENT) $(OBJEVENT) $(LIBS) $(LIBDIRSYS)

eve1a: $(OBJEVE1A)
	$(FC) $(LDFLAGS) $(PROGEVE1A) $(OBJEVE1A)

eve2: $(OBJEVE2)
	$(FC) $(LDFLAGS) $(PROGEVE2) $(OBJEVE2) $(LIBS) $(LIBDIRSYS)

evew7: $(OBJEVEW7)
	$(FC) $(LDFLAGS) $(PROGEVEW7) $(OBJEVEW7) $(LIBS) $(LIBDIRSYS)

eve2NoNi: $(OBJEVE2NONI)
	$(FC) $(LDFLAGS) $(PROGEVE2NONI) $(OBJEVE2NONI) $(LIBS) $(LIBDIRSYS)

evewindsi: $(OBJEVEWS)
	$(FC) $(LDFLAGS) $(PROGEVEWS) $(OBJEVEWS) $(LIBS) $(LIBDIRSYS)

stellarel: $(OBJSTR)   #   quasi-relativistic
	$(FC) $(LDFLAGS) $(PROGSTR) $(OBJSTR) $(LIBS)

stella6: $(OBJSTNR6)
	$(FC) $(LDFLAGS) $(PROGSTNR6) $(OBJSTNR6) $(LIBS)

stella6_mesa: $(OBJSTNR6_MESA)
	$(FC) $(LDFLAGS) $(PROGSTNR6_MESA) $(OBJSTNR6_MESA) $(LIBS)

stella6y12m: $(OBJSTNR6Y12M)
	$(FC) $(LDFLAGS) $(PROGSTNR6Y12M) $(OBJSTNR6Y12M) $(LIBS)

stella6magn: $(OBJSTNR6magn)
	$(FC) $(LDFLAGS) $(PROGSTNR6magn) $(OBJSTNR6magn) $(LIBS)

stellacon: $(OBJSTCON)
	$(FC) $(LDFLAGS) $(PROGSTCON) $(OBJSTCON) $(LIBS)

stellaWCWeaver: $(OBJSTNR6WEAVER)
	$(FC) $(LDFLAGS) $(PROGSTNR6WEAVER) $(OBJSTNR6WEAVER) $(LIBS)

stellaNCWeaver: $(OBJSTNR6WEAVNC)
	$(FC) $(LDFLAGS) $(PROGSTNR6WEAVNC) $(OBJSTNR6WEAVNC) $(LIBS)

stellafeau: $(OBJSTFEAU)
	$(FC) $(LDFLAGS) $(PROGSTFEAU) $(OBJSTFEAU) $(LIBS)

stellafeaou: $(OBJSTFEAOU)
	$(FC) $(LDFLAGS) $(PROGSTFEAOU) $(OBJSTFEAOU) $(LIBS)

stellaLowFreq: $(OBJSTLFRQ)
	$(FC) $(LDFLAGS) $(PROGSTLFRQ) $(OBJSTLFRQ) $(LIBS)

stellaOutUni: $(OBJSTOUTUNI)
	$(FC) $(LDFLAGS) $(PROGSTOUTUNI) $(OBJSTOUTUNI) $(LIBS)

stellatst: $(OBJSTNRTST)
	$(FC) $(LDFLAGS) $(PROGSTNRTST) $(OBJSTNRTST) $(LIBS)

stella7BQeps: $(OBJSTNR7BQEPS)
	$(FC) $(LDFLAGS) $(PROGSTNR7BQEPS) $(OBJSTNR7BQEPS) $(LIBS)

ttfit:  $(OBJTT)
	$(FC) $(LDFLAGS) $(PROGTT) $(OBJTT) $(LIBS)

ttfitou:  $(OBJTTOU)
	$(FC) $(LDFLAGS) $(PROGTTOU) $(OBJTTOU) $(LIBS) $(LIBSM)

ttfits:  $(OBJTTS)
	$(FC) $(LDFLAGS) $(PROGTTS) $(OBJTTS) $(LIBS) $(LIBSM)

ettfits:  $(OBJTTSE)
	$(FC) $(LDFLAGS) $(PROGTTSE) $(OBJTTSE) $(LIBS) $(LIBSM)

ttfitps:  $(OBJTTPS)
	$(FC) $(LDFLAGS) $(PROGTTPS) $(OBJTTPS) $(LIBS)

ronfictd:  $(OBJS_RONFD)
	$(FC) $(LDFLAGS) $(PROGRONFICTD) $(OBJS_RONFD) $(LIBS)

ronfict:  $(OBJS_RONF)
	$(FC) $(LDFLAGS) $(PROGRONFICT) $(OBJS_RONF) $(LIBS)

ronficttf:  $(OBJS_RONFTF)
	$(FC) $(LDFLAGS) $(PROGRONFTF) $(OBJS_RONFTF) $(LIBS)

ronfshb:  $(OBJS_RONFSHB)
	$(FC) $(LDFLAGS) $(PROGRONFSHB) $(OBJS_RONFSHB) $(LIBS)

ronfshbtf:  $(OBJS_RONFSHBTF)
	$(FC) $(LDFLAGS) $(PROGRONFSHBTF) $(OBJS_RONFSHBTF) $(LIBS)

ronfshbvR:  $(OBJS_RONFSHBVR)
	$(FC) $(LDFLAGS) $(PROGRONFSHBVR) $(OBJS_RONFSHBVR) $(LIBS)

bftot: $(OBJS_BFTOT)
	$(FC) $(LDFLAGS) $(PROGBFTOT) $(OBJS_BFTOT) $(LIBS)

bftotshb: $(OBJS_BFTOTSHB)
	$(FC) $(LDFLAGS) $(PROGBFTOTSHB) $(OBJS_BFTOTSHB) $(LIBS)

inner: $(OBJS_INSH)
	$(FC) $(LDFLAGS) $(PROGINSH)  $(OBJS_INSH) $(LIBS)

innersh: $(OBJS_INSHSHB)
	$(FC) $(LDFLAGS) $(PROGINSHSHB)  $(OBJS_INSHSHB) $(LIBS)

ronfLfrq:  $(OBJS_RONLFRQ)
	$(FC) $(LDFLAGS) $(PROGRONLFRQ) $(OBJS_RONLFRQ) $(LIBS)

ronfictOutUni:  $(OBJS_RONFOUTUNI)
	$(FC) $(LDFLAGS) $(PROGRONFICTOUTUNI) $(OBJS_RONFOUTUNI) $(LIBS)

opatold:  $(OBJS_OPATOLD)
	$(FC) $(LDFLAGS) $(PROGOPATOLD) $(OBJS_OPATOLD)  $(LIBS)

opatinsh:  $(OBJS_OPATINSH)
	$(FC) $(LDFLAGS) $(PROGOPATINSH) $(OBJS_OPATINSH)  $(LIBS)

opatbftot:  $(OBJS_OPATBFTOT)
	$(FC) $(LDFLAGS) $(PROGOPATBFTOT) $(OBJS_OPATBFTOT)  $(LIBS)


help:
	@echo "!!! DO NOT FORGET to define HOMEStella!!! "
	@echo "You can do: "
	@echo " evenewpab      --  compile evenewpab"
	@echo " evenew         --  compile evenew smoothing kepler models"
	@echo " event          --  compile event -- smoothing NTominaga models"
	@echo " eve1a          --  compile eve1a"
	@echo " eve2           --  compile eve2 for models from modmake"
	@echo " evew7          --  compile evew7 for models W7* from modmake"
	@echo " eve2NoNi       --  compile eve2NoNi" no Ni floor value
	@echo " evewindsi      --  compile evewindsi"
	@echo " ronfict        --  compile opacity tables in vladsf (use better rparlnx.mak)"
	@echo " ronfictd       --  compile opacity tables in vladsf for earlier time"
	@echo " ronficttf      --  compile opacity tables t_fixed"
	@echo " ronfshb        --  compile opacity tables for shocks (uniform composition)"
	@echo " ronfshbvR      --  compile opacity tables for shocks with vRegemorter lines"
	@echo " bftot          --  compile opacity tables with bf from excited levels"
	@echo " bftotshb       --  opacity tables  bf from excited levels -- shocks uniform comp."
	@echo " inner          --  inner-shell opacity tables (use better rparlnx.mak)"
	@echo " innersh        --  inner-shell opacity tables for shocks (uniform composition)"
	@echo " ronfLfrq       --  compile opacity tables in vladsf for Low Freq"
	@echo " ronfictOutUni  --  compile ronfictOutUni in vladsf: outer zones have uniform composition"
	@echo " opatold        --  build opatold.exe with opacityt.trf in vladsf"
	@echo " opatinsh       --  build opatinsh.exe with opacityInSh.trf in vladsf"
	@echo " opatbftot      --  build opatbftot.exe with opacitybftot.trf in vladsf"
	@echo " stellarel      --  compile strad with quasi-relativism"
	@echo " stellaOutUni   --  compile NR stella for hapOutUni: outer zones have uniform composition"
	@echo " stella6        --  compile strad without relativism and with introduced opacity"
	@echo " stella6y12m    --  compile strad with y12m without relativism and with introduced opacity"
	@echo " stella6magn    --  compile strad without relativism and with introduced opacity and magnetar triviality"
	@echo " stellaWithComp --  compile strad like stella6 but with approx. compton"
	@echo " stellacon      --  compile strad without relativism and with 'constant' opacity"
	@echo " stellaWCWeaver --  compile strad like stella6 but with approx. compton for Weaver's problem"
	@echo " stellaNCWeaver --  compile strad like stella6 (no Compton) for Weaver's problem"
	@echo " stellafeau     --  as stella6, but feauw used in place of eddi"
	@echo " stellafeaou    --  as stellafeau,  but outer zones have uniform composition"
	@echo " stellaLowFreq  --  compile strad with Low Freq opacity"
	@echo " stella7BQeps   --  with new stiffBGHeps large BQ in deep via RTphi"
	@echo " ttfit          --  compile ttfit"
	@echo " ttfits         --  compile ttfit simple with sm"
	@echo " ttfitou        --  compile ttfitou: outer zones have uniform composition "
	@echo " ettfits        --  compile ttfit with eddiplot"
	@echo " ttfitps        --  compile ttfitps with ps_write"
	@echo " clean          --  rm -f *.o ../strad/*.f ../src/*.f ../eve/*.f"
	@echo " cleandata      --  rm -f *.o core*"
	@echo " help           --  print this help"

#
# Here are all the dependencies:
#

# '../eve', '../src' etc  not needed if VPATH works:

######
# EVE
######
evenewpab.o  :  evenewpab.trf $(INCL)
evenewpa.o   :  evenewpa.trf $(INCL)
evenew.o     :  evenew.trf $(INCL)
event.o      :  event.trf $(INCL)
eveflexwindSi.o : eveflexwindSi.trf $(INCL)
evewindSi.o  :  evewindSi.trf $(INCL)
evewind.o    :  evewind.trf $(INCL)
evesn1a.o    :  evesn1a.trf $(INCL)
evesn2.o     :  evesn2.trf $(INCL)
evesn2w7.o   :  evesn2w7.trf $(INCL)
evesn2NoNi.o:    evesn2NoNi.trf $(INCL)

evesn2crab.o :  evesn2crab.trf $(INCL)
evesn2kep.o  :  evesn2kep.trf $(INCL)
evebeg.o     :  evebeg.trf $(INCL)
bgconnt.o    :  bgconnt.trf $(INCL)
bgconhep.o   :  bgconhep.trf $(INCL)
bgconhepb.o  :  bgconhepb.trf $(INCL)
resdia.o     :  resdia.trf $(INCL)
evefun.o     :  evefun.trf $(INCL)
evesub.o     :  evesub.trf $(INCL)
forback.o    :  forback.trf $(INCL)
sahaz.o      :  sahaz.trf $(INCL)
############
# STELLA
############
stradsep5tt.o:    stradsep5tt.trf $(INCL)
stradweaver.o:    stradweaver.trf $(INCL)
stradrezon5tt.o:  stradrezon5tt.trf $(INCL)
stradsep5ttrel.o: stradsep5ttrel.trf $(INCL)
stradsep5.o:      stradsep5.trf $(INCL)
stradrada.o:      stradrada.trf $(INCL)
stradOutUni.o:    stradOutUni.trf $(INCLOUTUNI)
stradsep7ttBQ.o:  stradsep7ttBQ.trf $(INCL)
begradOutUni.o:   begradOutUni.trf $(INCLOUTUNI)
begradsep.o:      begradsep.trf $(INCL)
begradLowFreq.o:  begradLowFreq.trf $(INCL)
cosetbgh.o:       cosetbgh.trf
dfnradrel.o:      dfnradrel.trf $(INCL)
dfjradrel.o:      dfjradrel.trf $(INCL)
hcdfjrad.o:       hcdfjrad.trf $(INCL)
hcdfnrad.o:       hcdfnrad.trf $(INCL)
traneq.o:         traneq.trf $(INCL)
eddi.o:           eddi.trf $(INCL)
obfeaw.o:         obfeaw.trf $(INCL)
obfearel.o:       obfearel.trf $(INCL)
feauj2.o:         feauj2.trf $(INCL)
gdepos6.o:        gdepos6.trf $(INCL)
nthnew.o:         nthnew.trf $(INCL)
stiffbgh.o:       stiffbgh.trf $(INCL)
stiffbghY12m.o:   stiffbghY12m.trf $(INCL)
stiffbghEps.o:    stiffbghEps.trf $(INCL)
y12m.o        :   y12m.f
lbalsw.o:         lbalsw.trf $(INCL)
lbalswrel.o:      lbalswrel.trf $(INCL)
lbalswBUG.o:      lbalswBUG.trf $(INCL)
hcdhaph.o:        hcdhaph.trf $(INCL)
obsubvri.o:       obsubvri.trf $(INCL)
#findTradByMNKLum.o:  findTradByMNKLum.trf $(INCL)
begtt.o:          begtt.trf $(INCL)
begttOutUni.o:    begttOutUni.trf $(INCLOUTUNI)
lbol.o:           lbol.trf $(INCL)

###########
# TTFIT
###########
tt4strad.o  : tt4strad.trf $(INCL)
ttfit5.o    :  ttfit5.trf $(INCL)
eddiplot.o  :  eddiplot.trf $(INCL)
ps_write.o  :  ps_write.f90

#########
# VLADSF
#########
ronfndec.o:        ronfndec.trf  $(INCL) $(SRC_DIR_FORT)nico.inc
ronfndecd.o:       ronfndecd.trf  $(INCL) $(SRC_DIR_FORT)nico.inc
ronftfixed.o:      ronftfixed.trf  $(INCL) $(SRC_DIR_FORT)nico.inc
ronfndecLFRQ.o:    ronfndecLFRQ.trf $(INCL) $(SRC_DIR_FORT)nico.inc
ronfshb.o:         ronfshb.trf  $(INCL) $(SRC_DIR_FORT)nico.inc
bessi0.o:          bessi0.f
bessk0ex.o:        bessk0ex.f
blas.o:            blas.f
dmach.o:           dmach.trf
edensol.o:         edensol.trf
gffcalc.o:         gffcalc.f
gshfdxsec.o:       gshfdxsec.f
hydxsecl.o:        hydxsecl.f
hydxsecn.o:        hydxsecn.f
hypho.o:           hypho.f
kinds.o         :  kinds.F
lineexpab_cor.o:   lineexpab_cor.trf
lineexpVRegem.o:   lineexpVRegem.trf
lnblnk.o:          lnblnk.f
ndex.o:            ndex.f
opacity.o:         opacity.trf $(SRC_DIR_FORT)zone.inc
opacityInSh.o:     opacityInSh.trf ../src/zone.inc # ../src/opacity.inc
opacitybftot.o:    opacitybftot.trf ../src/zone.inc ../vladsf/wmbasron.inc # ../src/opacity.inc
pfsaha.o:          pfsaha.f
sahaeqn.o:         sahaeqn.f
setnucms.o:        setnucms.f
sparseblas.o:      sparseblas.f
tablsort.o:        tablsort.f
valence_nl.o:      valence_nl.f
hMinusAbsorp.o:    hMinusAbsorp.trf
opatest.o :        opatest.trf $(INCL)
opatestbftot.o :   opatestbftot.trf $(INCL) ../src/nlte_param1.inc ../vladsf/mblock.inc

######
# SRC
######
burnc.o:           burnc.trf
sahaandd.o:        sahaandd.trf $(INCL)
ubv.o:             ubv.trf $(INCL)
volenpum.o:        volenpum.trf $(INCL)
volenpumnoint.o:   volenpumnoint.trf $(INCL)
volenpumnointmagn.o:volenpumnointmagn.trf $(INCL)
opazr.o:           opazr.trf $(INCL)
stradio.o:         stradio.trf $(INCL)
azdat.o:           azdat.trf $(INCL)
vtime.o:           vtime.trf
vtimeGF.o:         vtimeGF.trf
vtimef90.o:        vtimef90.trf
vecsubs.o:         vecsubs.trf
hapsepnc.o:        hapsepnc.trf $(INCL)
hapsepWithComp.o:  hapsepWithComp.trf $(INCL)
hapOutUni.o:       hapOutUni.trf $(INCL)
hap3.o:            hap3.trf $(INCL)
oparon.o:          oparon.trf $(INCL)
oparonOutUni.o:    oparonOutUni.trf $(INCLOUTUNI)
hapconst.o:        hapconst.trf $(INCL)
opaconst.o:        opaconst.trf $(INCL)
length.o:          length.trf
words.o:           words.trf

###################################################


cleandata:
	$(DEL) ../strad/run/nohup.out ../strad/run/core* ../strad/run/fort.*

cleantt:
	$(DEL) $(OBJTT)

cleanf:
	$(DEL) *.o  *.il ../strad/*.f ../src/*.f ../eve/*.f
	$(DEL) *.mod ../strad/*.lst ../src/*.lst ../eve/*.lst
	( cd ../vladsf ; \
        	pwd ; \
	$(DEL) *.o *.mod  mach.f dmach.f edensol.f linexpabav.f lineexpab_cor.f lineexpfast.f \
               lineexpfastKur.f lineexpfastKur_taudist.f linterpol.f \
               ronfshb.f ronfndec.f ronftfixed.f ronfshbtfixed.f ronfndecbftot.f ronfshbbftot.f \
               opatest.f opatestgplot.f opatestswd.f opatestbftot.f opatestbf.f opaabtest.f \
               opatestnltab.f opacity100.f opacityt.f opacita3.f opacitybftot.f opacitybf.f \
               opacity.f opacitynltab.f vtime.f opacityInSh.f hMinusAbsorp.f ; )

clean:
	$(DEL) *.o  *.il
	$(DEL) *.mod ../strad/*.lst ../src/*.lst ../eve/*.lst
	( cd ../vladsf ; \
        	pwd ; \
	$(DEL) *.o *.mod ; )

