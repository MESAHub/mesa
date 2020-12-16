      subroutine get_eos_info_for_burn_at_const_density( &
            z, xh, abar, zbar, xa, rho, logRho, T, logT, &
            Cv, d_Cv_dlnT, eta, d_eta_dlnT, ierr)
         double precision, intent(in) :: &
            z, xh, abar, zbar, xa(:), rho, logRho, T, logT
         double precision, intent(out) :: &
            Cv, d_Cv_dlnT, eta, d_eta_dlnT
         integer, intent(out) :: ierr
      end subroutine get_eos_info_for_burn_at_const_density
