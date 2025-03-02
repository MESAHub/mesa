include ../make/defaults-module.mk

# Build

MODULE_NAME := num
SRCS := public/accurate_sum.f90 \
        public/accurate_sum_auto_diff_star_order1.f90 \
        public/num_def.f \
        public/num_lib.f \
        private/mod_dc_decsol.f \
        private/mod_dop853.f \
        private/mod_dopri5.f \
        private/mod_integrate.f90 \
        private/mod_qsort.f90 \
        private/mod_random.f90 \
        private/mod_root.f \
        private/mod_rosenbrock.f \
        private/mod_simplex.f90 \
        private/mod_solout_root.f \
        private/mod_bobyqa.f \
        private/mod_brent.f90 \
        private/mod_newton.f \
        private/mod_newuoa.f
SRCS_CHECK := test/src/bari_beam.f \
              test/src/bari_chemakzo.f \
              test/src/bari_medakzo.f \
              test/src/bari_vdpol.f \
              test/src/test_beam.f90 \
              test/src/test_bobyqa.f90 \
              test/src/test_brent.f90 \
              test/src/test_chemakzo.f90 \
              test/src/test_diffusion.f90 \
              test/src/test_int_support.f90 \
              test/src/test_integrate.f90 \
              test/src/test_medakzo.f90 \
              test/src/test_newton.f90 \
              test/src/test_newuoa.f90 \
              test/src/test_num.f90 \
              test/src/test_simplex.f90 \
              test/src/test_support.f90 \
              test/src/test_vdpol.f90
INTERNAL_DEPENDS_ON := const utils math mtx auto_diff
EXTERNAL_DEPENDS_ON :=
BINTYPE := static-lib

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := num_def.mod num_lib.mod accurate_sum.mod accurate_sum_auto_diff_star_order1.mod
INSTALL_INCLUDES := public/newton_procs.dek \
                    public/num_binary_search.dek \
                    public/num_bobyqa.dek \
                    public/num_bobyqa_proc.dek \
                    public/num_brent.dek \
                    public/num_dfridr.dek \
                    public/num_dfx.dek \
                    public/num_dop853.dek \
                    public/num_dopri5.dek \
                    public/num_fcn.dek \
                    public/num_fcn_blk_dble.dek \
                    public/num_interp_y.dek \
                    public/num_isolve.dek \
                    public/num_jac.dek \
                    public/num_jac_blk_dble.dek \
                    public/num_mas.dek \
                    public/num_newton.dek \
                    public/num_newuoa.dek \
                    public/num_newuoa_proc.dek \
                    public/num_radau5_debug_routine.dek \
                    public/num_rk2.dek \
                    public/num_rk4.dek \
                    public/num_rodas3.dek \
                    public/num_rodas4.dek \
                    public/num_rodasp.dek \
                    public/num_root_fcn.dek \
                    public/num_ros2.dek \
                    public/num_ros3p.dek \
                    public/num_ros3pl.dek \
                    public/num_rose2.dek \
                    public/num_safe_root.dek \
                    public/num_safe_root_x.dek \
                    public/num_simplex.dek \
                    public/num_simplex_fcn.dek \
                    public/num_sjac.dek \
                    public/num_solout.dek \
                    public/num_solout_root.dek \
                    public/num_solout_root_fcn.dek \
                    public/num_solvers_options.dek

include $(MAKE_DIR)/Makefile
