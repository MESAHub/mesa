include ../../make/defaults-module.mk

MODULE_NAME := kap-preprocessor
SRCS := src/co_enhanced.f90 \
	src/create_tables.f90 \
	src/ferg_logP.f90 \
	src/ferguson.f90 \
	src/fixed_metal.f90 \
	src/freedman.f90 \
	src/kap_support.f90 \
	src/opal_type1.f \
	src/opal_type2.f
SRCS_CHECK :=
INTERNAL_DEPENDS_ON := kap
EXTERNAL_DEPENDS_ON :=
CHECK_RESULTS_GOLDEN :=
INCLUDE_DIRS := -Isrc

# Install

MODULES :=
INSTALL_INCLUDES :=
BINTYPE := executable

include ../../make/Makefile

$(BUILD_DIR_MODULE)/kap_input_data/: kap_input_data.tar.xz
	tar -xf kap_input_data.tar.xz -C $(BUILD_DIR_MODULE)

$(BUILD_DIR_MODULE)/ags09/: AESOPUS/aesopus_ags09.tar.xz
	mkdir -p $@
	tar -xf AESOPUS/aesopus_ags09.tar.xz -C $@

$(BUILD_DIR_MODULE)/logT_points_Fe $(BUILD_DIR_MODULE)/logT_points_highT_8.0 $(BUILD_DIR_MODULE)/logT_points_highT_8.7 $(BUILD_DIR_MODULE)/logT_points_lowT_2.7 $(BUILD_DIR_MODULE)/logT_points_lowT_3.0 $(BUILD_DIR_MODULE)/logT_points_lowT_Freedman11:
	cp $(notdir $@) $@

logT_points: $(BUILD_DIR_MODULE)/logT_points_Fe $(BUILD_DIR_MODULE)/logT_points_highT_8.0 $(BUILD_DIR_MODULE)/logT_points_highT_8.7 $(BUILD_DIR_MODULE)/logT_points_lowT_2.7 $(BUILD_DIR_MODULE)/logT_points_lowT_3.0 $(BUILD_DIR_MODULE)/logT_points_lowT_Freedman11


define kap
@echo kap-preprocessor $1 $2 $3
@pushd $(BUILD_DIR_MODULE) > /dev/null; $(OBJ_OUT) $(CURDIR)/$1 $2 $3; popd > /dev/null
endef

define kap-fixed-metal
$(call kap,$1,0.0000,-1)
$(call kap,$1,0.0001,-1)
$(call kap,$1,0.0003,-1)
$(call kap,$1,0.0010,-1)
$(call kap,$1,0.0020,-1)
$(call kap,$1,0.0040,-1)
$(call kap,$1,0.0100,-1)
$(call kap,$1,0.0200,-1)
$(call kap,$1,0.0300,-1)
$(call kap,$1,0.0400,-1)
$(call kap,$1,0.0600,-1)
$(call kap,$1,0.0800,-1)
$(call kap,$1,0.1000,-1)
endef

define kap-co-enhanced
$(call kap,$1,0.001,-1)
$(call kap,$1,0.004,-1)
$(call kap,$1,0.010,-1)
$(call kap,$1,0.020,-1)
$(call kap,$1,0.030,-1)
$(call kap,$1,0.050,-1)
$(call kap,$1,0.100,-1)
endef

define kap-freedman
$(call kap,$1,0.01,0)
$(call kap,$1,0.02,0)
$(call kap,$1,0.04,0)
$(call kap,$1,0.10,0)
$(call kap,$1,0.20,0)
$(call kap,$1,0.63,0)
$(call kap,$1,1.00,0)
endef

define kap-oplib
python $(CURDIR)/spline_OPLIB_tables.py $(BUILD_DIR_MODULE)/kap_input_data/oplib/$1/ $(BUILD_DIR_MODULE)/data/kap_data/
endef

build-data: $(BUILD_DIR_MODULE)/kap_input_data/ logT_points $(OBJ_OUT) | $(BUILD_DIR_MODULE)/data/kap_data/
	# As of r8118 file is empty kap_input_data/opal/A05+Ne
	# do_one inlist_a05+Ne
# This is broken (SIGFPE due to negative opacity), 1.0 is way outside the table
# boundaries (also doesn't seem to be in use)
#	$(call kap,inlist_Fe,1.0,0.0)

	$(call kap-fixed-metal,inlist_a09)
	$(call kap-co-enhanced,inlist_a09_CO)

	$(call kap-fixed-metal,inlist_OP_a09_nans_removed_by_hand)
	$(call kap-fixed-metal,inlist_OP_gs98)

	$(call kap-freedman,inlist_lowT_Freedman11)

	$(call kap-fixed-metal,inlist_lowT_fa05_gs98)
	$(call kap-fixed-metal,inlist_lowT_fa05_gs98_aFe_m2)
	$(call kap-fixed-metal,inlist_lowT_fa05_gs98_aFe_p2)
	$(call kap-fixed-metal,inlist_lowT_fa05_gs98_aFe_p4)
	$(call kap-fixed-metal,inlist_lowT_fa05_gs98_aFe_p6)
	$(call kap-fixed-metal,inlist_lowT_fa05_gs98_aFe_p8)

	$(call kap-fixed-metal,inlist_lowT_fa05_gn93)
	$(call kap-fixed-metal,inlist_lowT_fa05_ags04)
	$(call kap-fixed-metal,inlist_lowT_fa05_a09)

	$(call kap-fixed-metal,inlist_lowT_af94_gn93)

	$(call kap-fixed-metal,inlist_gn93)
	$(call kap-co-enhanced,inlist_gn93_CO)

	$(call kap-fixed-metal,inlist_gs98)
	$(call kap-co-enhanced,inlist_gs98_CO)
	$(call kap-fixed-metal,inlist_gs98_aFe_m2)
	$(call kap-fixed-metal,inlist_gs98_aFe_p2)
	$(call kap-fixed-metal,inlist_gs98_aFe_p3)
	$(call kap-fixed-metal,inlist_gs98_aFe_p4)
	$(call kap-co-enhanced,inlist_gs98_aFe_p4_CO)
	$(call kap-fixed-metal,inlist_gs98_aFe_p6)
	$(call kap-fixed-metal,inlist_gs98_aFe_p8)

	$(call kap-oplib,OPLIB_GS98_LOG_1194_25E)
	$(call kap-oplib,OPLIB_AGSS09_LOG_1194_25E)
	$(call kap-oplib,OPLIB_AAG21_LOG_1194_25E)
	$(call kap-oplib,OPLIB_MB22_LOG_1194_25E)

	$(call kap-fixed-metal,inlist_lowT_fa05_mb22)

	cp -r $(BUILD_DIR_MODULE)/kap_input_data/wichita/AAG21_lowT_MESA_tables/*  $(BUILD_DIR_MODULE)/data/kap_data/
	cp condtable.data $(BUILD_DIR_MODULE)/data/kap_data/

$(BUILD_DIR_MODULE)/AESOPUS_AGSS09.h5: $(BUILD_DIR_MODULE)/ags09/
	pushd $(BUILD_DIR_MODULE) > /dev/null; python $(CURDIR)/AESOPUS/aesopus.py $(CURDIR)/AESOPUS/AGSS09.yaml ; popd > /dev/null

build-aesopus: $(BUILD_DIR_MODULE)/AESOPUS_AGSS09.h5

export-data:
	rm -f ../kap_data.tar ../kap_data.tar.xz ../AESOPUS_AGSS09.h5
	tar -cRf ../kap_data.tar -C $(BUILD_DIR_MODULE)/data kap_data
	xz ../kap_data.tar
	cp $(BUILD_DIR_MODULE)/AESOPUS_AGSS09.h5 ..

.PHONY: build-data export-data
