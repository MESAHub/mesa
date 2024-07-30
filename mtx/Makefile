include ../make/defaults-module.mk

# Build

MODULE_NAME := mtx
SRCS := public/mtx_def.f \
        public/mtx_lib.f90 \
        private/DGBSVX_wrapper.f90 \
        private/bcyclic.f90 \
        private/mtx_support.f90 \
        private/pre_conditioners.f90 \
        private/my_lapack95.f90
SRCS_CHECK := test/src/test_mtx.f90 \
              test/src/test_mtx_support.f90 \
              test/src/test_square.f90 \
              test/src/test_block_tridiagonal.f90
# TODO: Also contains test dependencies sadly, I should fix this at some point
INTERNAL_DEPENDS_ON := const utils math
EXTERNAL_DEPENDS_ON := lapack
BINTYPE := static-lib

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := mtx_def.mod mtx_lib.mod
INSTALL_INCLUDES := public/mtx_bcyclic_dble_decsol.dek \
                    public/mtx_dble_decsol.dek \
                    public/mtx_debug_decsol.dek \
                    public/mtx_decsol.dek \
                    public/mtx_decsolblk.dek \
                    public/mtx_decsolblk_dble.dek \
                    public/mtx_decsolc.dek \
                    public/mtx_decsolcs.dek \
                    public/mtx_decsols.dek \
                    public/mtx_jac.dek \
                    public/mtx_lapack95.dek \
                    public/mtx_null_decsol.dek \
                    public/mtx_rcond.dek \
                    public/mtx_tridiag_decsol.dek \
                    public/mtx_solve_routines.inc

include $(MAKE_DIR)/Makefile
