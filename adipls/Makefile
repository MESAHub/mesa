include ../make/defaults-module.mk

FFLAGS += -fno-check=bounds

# Build

MODULE_NAME := adipls
SRCS := adipack.c/adipls/adipls.c.d.f \
        adipack.c/adipls/adirhs.c.d.f \
        adipack.c/adipls/anldet.d.f \
        adipack.c/adipls/cmparr.d.f \
        adipack.c/adipls/danorm.d.f \
        adipack.c/adipls/decicd.d.f \
        adipack.c/adipls/delta1.d.f \
        adipack.c/adipls/delta2.d.f \
        adipack.c/adipls/delta3.d.f \
        adipack.c/adipls/dmpgsm.d.f \
        adipack.c/adipls/dmpmde.d.f \
        adipack.c/adipls/dmpsum.d.f \
        adipack.c/adipls/dump_obs.c.d.f \
        adipack.c/adipls/geninh.d.f \
        adipack.c/adipls/gm1ker.n.d.f \
        adipack.c/adipls/gravpo.n.d.f \
        adipack.c/adipls/kiner.n.d.f \
        adipack.c/adipls/mchcff.d.f \
        adipack.c/adipls/mchsol.n.d.f \
        adipack.c/adipls/nrkint.c.d.f \
        adipack.c/adipls/nrkm.d.f \
        adipack.c/adipls/nrtssl.d.f \
        adipack.c/adipls/order.d.f \
        adipack.c/adipls/rdfrqr.d.f \
        adipack.c/adipls/readml.n.d.f \
        adipack.c/adipls/readsm.d.f \
        adipack.c/adipls/res_adimod.c.d.f \
        adipack.c/adipls/res_adipar.c.d.f \
        adipack.c/adipls/rotker.c.d.f \
        adipack.c/adipls/rsnfit.n.d.f \
        adipack.c/adipls/sclasl.c.d.f \
        adipack.c/adipls/setbcs.n.d.f \
        adipack.c/adipls/setobs_st.c.d.f \
        adipack.c/adipls/set_rotation.c.d.f \
        adipack.c/adipls/setssm.n.d.f \
        adipack.c/adipls/setsum.d.f \
        adipack.c/adipls/setups_adi.c.d.f \
        adipack.c/adipls/shtint.c.d.f \
        adipack.c/adipls/sigint.n.d.f \
        adipack.c/adipls/sigitc.c.d.f \
        adipack.c/adipls/sigitr.n.d.f \
        adipack.c/adipls/signew.n.d.f \
        adipack.c/adipls/sigout.c.d.f \
        adipack.c/adipls/sigscn.c.d.f \
        adipack.c/adipls/sigsol.c.d.f \
        adipack.c/adipls/sph.d.f \
        adipack.c/adipls/stevft.n.d.f \
        adipack.c/adipls/sumout.d.f \
        adipack.c/adipls/takata.d.f \
        adipack.c/adipls/trnmod.n.d.f \
        adipack.c/adipls/uhx.d.f \
        adipack.c/adipls/user_adi.c.d.f \
        adipack.c/adipls/varfrq.n.d.f \
        adipack.c/adipls/xlmult.d.f \
        adipack.c/gensr/aramax.d.f \
        adipack.c/gensr/blstio.f \
        adipack.c/gensr/cpycom.f \
        adipack.c/gensr/craysub.d.f \
        adipack.c/gensr/derive.d.f \
        adipack.c/gensr/derivk.new.d.f \
        adipack.c/gensr/eiginh.d.f \
        adipack.c/gensr/getcas.f \
        adipack.c/gensr/intgpt.d.f \
        adipack.c/gensr/iprsgn.d.f \
        adipack.c/gensr/izero.f \
        adipack.c/gensr/length.f \
        adipack.c/gensr/leq.d.f \
        adipack.c/gensr/lininh4.d.f \
        adipack.c/gensr/lininh.new.d.f \
        adipack.c/gensr/lir.d.f \
        adipack.c/gensr/lirt.d.f \
        adipack.c/gensr/lsqpol.d.f \
        adipack.c/gensr/lstsqr.f \
        adipack.c/gensr/notwin.d.f \
        adipack.c/gensr/notwni.f \
        adipack.c/gensr/nrk.d.f \
        adipack.c/gensr/ofiles.c.f \
        adipack.c/gensr/ofiles.f \
        adipack.c/gensr/rdilst.f \
        adipack.c/gensr/rnmean.d.f \
        adipack.c/gensr/skpcom.f \
        adipack.c/gensr/sort.d.f \
        adipack.c/gensr/spline.d.f \
        adipack.c/gensr/squint.d.f \
        adipack.c/gensr/store.d.f \
        adipack.c/gensr/strcompr.f \
        adipack.c/gensr/vinta.d.f \
        adipack.c/gensr/vintk.new.d.f \
        adipack.c/gensr/zero.d.f \
        adipack.c/adiajobs/sr/inrnge.f \
        adipack.c/adiajobs/sr/intmdl.d.f \
        adipack.c/adiajobs/sr/rdamdl.n.d.f \
        adipack.c/adiajobs/sr/rdfreq.d.f \
        adipack.c/adiajobs/sr/rdfrqe.d.f \
        adipack.c/adiajobs/sr/rdfrqm.d.f \
        adipack.c/adiajobs/sr/rseta4.c.d.f \
        adipack.c/adiajobs/sr/setssm.d.f \
        adipack.c/adiajobs/sr/srdist.c.d.f \
        adipack.c/adiajobs/sr/wrfreq.d.f \

SRCS_CHECK :=
INTERNAL_DEPENDS_ON := const
EXTERNAL_DEPENDS_ON :=
BINTYPE = static-lib
INCLUDE_DIRS := -Iadipack.c/adipls

CHECK_RESULTS_GOLDEN :=

include $(MAKE_DIR)/Makefile
