include ../make/defaults-module.mk

# Build

MODULE_NAME := binary
SRCS := \
   public/binary_pgbinary.f90 \
   public/binary_def.f90 \
   private/binary_private_def.f90 \
   private/binary_utils.f90 \
   private/binary_photos.f90 \
   private/binary_jdot.f90 \
   private/binary_wind.f90 \
   private/binary_mdot.f90 \
   private/binary_tides.f90 \
   private/binary_edot.f90 \
   private/binary_ctrls_io.f90 \
   private/binary_irradiation.f90 \
   private/binary_history_specs.f90 \
   private/binary_history.f90 \
   private/binary_evolve.f90 \
   private/binary_timestep.f90\
   private/binary_ce.f90\
   private/binary_job_ctrls_io.f90\
   private/binary_do_one_utils.f90\
   private/run_binary_support.f90 \
   public/binary_lib.f90 \
   other/mod_other_rlo_mdot.f90 \
   other/mod_other_tsync.f90 \
   other/mod_other_sync_spin_to_orbit.f90 \
   other/mod_other_mdot_edd.f90 \
   other/mod_other_adjust_mdots.f90 \
   other/mod_other_accreted_material_j.f90 \
   other/mod_other_binary_jdot.f90 \
   other/mod_other_binary_wind_transfer.f90 \
   other/mod_other_binary_edot.f90 \
   other/mod_other_binary_ce.f90 \
   other/mod_other_implicit_rlo.f90 \
   other/mod_other_binary_extras.f90 \
   other/mod_other_binary_photo_read.f90 \
   other/mod_other_binary_photo_write.f90 \
   other/mod_other_e2.f90 \
   other/mod_other_pgbinary_plots.f90 \
   other/pgbinary_decorator.f90

ifeq ($(WITH_PGSTAR),yes)
SRCS += \
   private/pgbinary_support.f90 \
   private/pgbinary_ctrls_io.f90 \
   private/pgbinary_summary.f90 \
   private/pgbinary_summary_history.f90 \
   private/pgbinary_hist_track.f90 \
   private/pgbinary_history_panels.f90 \
   private/pgbinary_orbit.f90 \
   private/pgbinary_star.f90 \
   private/pgbinary_grid.f90 \
   private/pgbinary_lib.f90 \
   private/pgbinary_full.f90
EXTERNAL_DEPENDS_ON := pgplot
else
SRCS += private/pgstar_stub.f90
EXTERNAL_DEPENDS_ON :=
endif

SRCS_CHECK = test/src/test_binary.f90
INTERNAL_DEPENDS_ON := utils star
BINTYPE := static-lib
INCLUDE_DIRS := -Iprivate -Ipublic -Idefaults

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := binary_lib.mod binary_def.mod binary_pgbinary.mod
INSTALL_INCLUDES := job/binary_test_suite_extras.inc job/binary_test_suite_extras_def.inc

include $(MAKE_DIR)/Makefile
