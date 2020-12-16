##############################
# MAKE FOR STELLA
##############################

#sys := $(shell uname)
#ifneq (,$(findstring CYGWIN,$(strip $(sys))))
#   SYSTYPE := "cygwin_f90"
#else
#  SYSTYPE := "ifort"
#  SYSTYPE="pgf"
#  SYSTYPE := "MPA_f95i"
#endif
# SYSTYPE="ifort"
# SYSTYPE := "lf"
# SYSTYPE=gfortran
# SYSTYPE="cygwin_ifort"
# SYSTYPE="cygwin_f90"
# SYSTYPE="MPA_f95i"
# SYSTYPE="MPA_f95n"
UnixOrWin = "unix" # default value

########################################
## For Intel fortran
ifeq ($(SYSTYPE),ifort)
  FC = ifort  # Intel fortran
#   FFLAGS = -c -save -zero -O3 -mcmodel=medium -shared-intel
#   FFLAGS = -c -save -zero -O3 -fp-model precise
#   FFLAGS = -c -save -zero -O3 -fp-model strict -CB -g -traceback # -- for debug
#   FFLAGS = -c -save -zero -fpe0 -no-ftz  -check bounds -g -traceback # -- for debug
#   FFLAGS = -c -save -zero -O3 -fp-model strict
#   FFLAGS = -c -O3 -save -zero -tpp7   -ip # -xW # -- optimize for Pentium
#   FFLAGS =   -c -save -static -zero -O2 -fpe0 -traceback # -CB -traceback
 FFLAGS_FIX = -c -132 -save -zero -O3 -fp-model strict
#   FFLAGS_FIX = -c -save -zero -check bounds -g -traceback -debug inline-debug-info # -- for debug
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
  LIBS=
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
  LIBS =
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
  LIBS=
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
  LIBS =
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
  LIBS=
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
  LIBS =
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
  LIBS=
  UnixOrWin := "unix"
endif
########################################


########################################
## For gfortran fortran
ifeq ($(SYSTYPE),gfortran)
  FC = gfortran
  FFLAGS_FIX = -c -g -fbounds-check  -fno-automatic
  FFLAGS_FIX = -c -g  -fno-automatic
  FFLAGS_FIX = -c -O2  -ffixed-line-length-132 -fno-automatic -Wall
  FFLAGS := -ffree-form
  FFLAGS := $(FFLAGS_FIX) $(FFLAGS)
  LDFLAGS =
  LIBS=
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
vpath %.inc ../src:../vladsf



VPATH := $(STLHOME)src $(STLHOME)src/cp2k $(STLHOME)src/stl $(STLHOME)src/util $(STLHOME)vladsf $(STLHOME)sparse
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


INCL = opacity.inc stateq.inc \
       abo.inc nstep.inc stsave.inc fundrad.inc \
       black.inc snrad.inc \
       zone.inc commonEve.inc azzn.inc

INCLOUTUNI = opacityOutUni.inc stateq.inc \
       abo.inc nstep.inc stsave.inc fundrad.inc \
       black.inc snrad.inc \
       zone.inc commonEve.inc azzn.inc

OBJDIR = obj

PROGEVE1A = ../eve/run/eve1a.exe
PROGEVE1A := $(subst /,$(slash_l),$(PROGEVE1A))

PROGEVE2 = ../eve/run/eve2.exe
PROGEVE2     := $(subst /,$(slash_l),$(PROGEVE2))


PROGSTNR6Y12M = ../strad/run/xstella6y12m.exe # non-relativ. with gdepos6 & volennoint
PROGSTNR6Y12M   := $(subst /,$(slash_l),$(PROGSTNR6Y12M))


PROGTT = ../strad/run/xttfit.exe
PROGTT := $(subst /,$(slash_l),$(PROGTT))


PROGRONFICT = ../vladsf/xronfict.exe
PROGRONFICT := $(subst /,$(slash_l),$(PROGRONFICT))



FILES_EVE1A =evesn1a.trf sahaz.trf sahaandd.trf opazr.trf stradio.trf \
         azdat.trf length.trf
OBJEVE1A := $(FILES_EVE1A:.trf=.o)


FILES_EVE2 =evesn2.trf sahaz.trf sahaandd.trf opazr.trf stradio.trf \
         azdat.trf length.trf
OBJEVE2 := $(FILES_EVE2:.trf=.o)


FILES_RONF = ronfndec.trf bessi0.f bessk0ex.f \
            blas.f dmach.trf edensol.trf \
            gffcalc.f gshfdxsec.f hydxsecl.f \
            hydxsecn.f hypho.f lineexpab_cor.trf \
            lnblnk.f ndex.f opacityt.trf \
            pfsaha.f sahaeqn.f setnucms.f \
            sparseblas.f tablsort.f valence_nl.f \
            azdat.trf length.trf stradio.trf \
            vtimef90.trf hMinusAbsorp.trf
OBJS_RONF := $(patsubst %.f,%.o, $(FILES_RONF:.trf=.f)) \
        $(OBJS_base) $(OBJScross)

FILES_TNR6Y12M = stradsep5tt.trf begradsep.trf cosetbgh.trf hcdfjrad.trf \
        hcdfnrad.trf traneq.trf eddi.trf  gdepos6.trf nthnew.trf \
        stiffbghY12m.trf lbalsw.trf stradio.trf \
        vtimef90.trf sahaandd.trf ubv.trf obsubvri.trf \
        tt4strad.trf begtt.trf lbol.trf \
        burnc.trf volenpumnoint.trf hapsepnc.trf hcdhaph.trf \
        oparon.trf length.trf words.trf azdat.trf y12m.f
OBJSTNR6Y12M := $(patsubst %.f,%.o, $(FILES_TNR6Y12M:.trf=.f))

FILES_TTFIT = ttfit5.trf begradsep.trf vtimef90.trf \
        nthnew.trf  lbalsw.trf    words.trf \
        stradio.trf  sahaandd.trf  ubv.trf \
        obsubvri.trf  burnc.trf  volenpumnoint.trf \
        hapsepnc.trf  oparon.trf  length.trf \
        azdat.trf lbol.trf
OBJTT := $(FILES_TTFIT:.trf=.o)


TREFOR = trf

# for ctrf:
TRFFLAGS = -nfs



.IGNORE:

%.o : %.trf
	$(TREFOR) $(TRFFLAGS) $<
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


eve1a: $(OBJEVE1A)
	$(FC) $(LDFLAGS) $(PROGEVE1A) $(OBJEVE1A)

eve2: $(OBJEVE2)
	$(FC) $(LDFLAGS) $(PROGEVE2) $(OBJEVE2) $(LIBS) $(LIBDIRSYS)

stella6y12m: $(OBJSTNR6Y12M)
	$(FC) $(LDFLAGS) $(PROGSTNR6Y12M) $(OBJSTNR6Y12M) $(LIBS)

ttfit:  $(OBJTT)
	$(FC) $(LDFLAGS) $(PROGTT) $(OBJTT) $(LIBS) $(LIBSM)

ronfict:  $(OBJS_RONF)
	$(FC) $(LDFLAGS) $(PROGRONFICT) $(OBJS_RONF) $(LIBS)

help:
	@echo "You can do: "
	@echo " eve1a          --  compile eve1a"
	@echo " eve2           --  compile eve2"
	@echo " ronfict        --  compile opacity tables in vladsf"
	@echo " stella6y12m    --  compile strad with y12m without relativism and with introduced opacity"
	@echo " ttfit          --  compile ttfit"
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
evesn1a.o:   evesn1a.trf $(INCL)
evesn2.o:    evesn2.trf $(INCL)
evesn2NoNi.o:    evesn2NoNi.trf $(INCL)

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
begradOutUni.o:   begradOutUni.trf $(INCLOUTUNI)
begradsep.o:      begradsep.trf $(INCL)
begradLowFreq.o:  begradLowFreq.trf $(INCL)
cosetbgh.o:       cosetbgh.trf
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
y12m.o        :   y12m.f
lbalsw.o:         lbalsw.trf $(INCL)
lbalswrel.o:      lbalswrel.trf $(INCL)
hcdhaph.o:        hcdhaph.trf $(INCL)
obsubvri.o:       obsubvri.trf $(INCL)
tt4strad.o:       tt4strad.trf $(INCL)
findTradByMNKLum.o:  findTradByMNKLum.trf $(INCL)
begtt.o:          begtt.trf $(INCL)
begttOutUni.o:    begttOutUni.trf $(INCLOUTUNI)
lbol.o:           lbol.trf $(INCL)

###########
# TTFIT
###########
ttfit5.o:  ttfit5.trf $(INCL)
eddiplot.o :  eddiplot.trf $(INCL)
ps_write.o :  ps_write.f90

#########
# VLADSF
#########
ronfndec.o:        ronfndec.trf  $(INCL) $(SRC_DIR_FORT)nico.inc
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
length.o:          length.trf
words.o:           words.trf

###################################################


cleandata:
	$(DEL) ../strad/run/nohup.out ../strad/run/core* ../strad/run/fort.*

cleantt:
	$(DEL) $(OBJTT)

clean:
	$(DEL) *.o  *.il ../strad/*.f ../src/*.f ../eve/*.f
	$(DEL) *.mod ../strad/*.lst ../src/*.lst ../eve/*.lst
	( cd ../vladsf ; \
        	pwd ; \
	$(DEL) mach.f dmach.f edensol.f linexpabav.f lineexpab_cor.f lineexpfast.f \
               lineexpfastKur.f lineexpfastKur_taudist.f linterpol.f \
               ronfshb.f ronfndec.f ronfndecbftot.f ronfshbbftot.f \
               opatest.f opatestgplot.f opatestswd.f opatestbftot.f opatestbf.f opaabtest.f \
               opatestnltab.f opacity100.f opacityt.f opacita3.f opacitybftot.f opacitybf.f \
               opacity.f opacitynltab.f vtime.f opacityInSh.f hMinusAbsorp.f ; )

