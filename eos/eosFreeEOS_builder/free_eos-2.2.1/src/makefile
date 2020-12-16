#############COMPILER STUFF###############
FC=gfortran -std=legacy -ffixed-line-length-132

WARN=#-Wall
CHECK= #-pedantic -fcheck=all
OPT=-O2

DEBUG=-g -fbacktrace

COMPILE=$(FC) $(OPT) $(WARN) $(CHECK) $(DEBUG)

INCLUDE = ../include

#############SOURCE FILES#################
SRCS = aux_to_traditional.f bfgs.f coulomb_adjust.f \
	coulomb.f deriv_product.f deriv_ratio.f \
	effective_radius.f effsum_calc.f eos_bfgs.f \
	eos_calc.f eos_cold_start.f eos_free_calc.f \
	eos_jacobian.f eos_sum_calc.f eos_tqft.f \
	eos_tqsum_calc.f eos_warm_step.f exchange_coeff.f \
	exchange_gcpf.f excitation_pi_end.f excitation_pi.f \
	excitation_sum.f exct_calc.f fermi_dirac_coeff.f \
	fermi_dirac_ct.f fermi_dirac_direct.f fermi_dirac_exchange.f \
	fermi_dirac.f fermi_dirac_minusnr_coeff.f fermi_dirac_original_coeff.f \
	fermi_dirac_recursion.f fjs_pi.f fp12_calc.f \
	f_psi.f free_eos_detailed.f free_eos.f \
	free_eos_interpolated.f free_non_ideal_calc.f ionize.f \
	lnblnk.f lnx_calc.f master_coulomb.f mdh_pi.f molecular_hydrogen.f \
	pl_prepare.f plsum_approx.f plsum.f poly_sum.f pteh_pi.f pteh_theta.f \
	qmhd_calc.f qryd_approx.f qryd_calc.f qstar_calc.f quad4b.f \
	round_ln.f tau_calc.f version.f

##########################################
LIB=libfree_eos.a
OBJS = $(patsubst %.f,%.o,$(SRCS))

$(LIB) : $(OBJS)
	ar crs $(LIB) $(OBJS)

##########################################
%.o : %.f
	$(COMPILE) -I$(INCLUDE) -c $<

install:
	mkdir -p ../lib
	cp $(LIB) ../lib
#ifdef MESA_DIR
#	cp $(LIB) $(MESA_DIR)/lib/
#	cp ../include/*.h $(MESA_DIR)/include/
#endif


clean:
	-@rm -f *.o $(LIB) ../lib/$(LIB)
