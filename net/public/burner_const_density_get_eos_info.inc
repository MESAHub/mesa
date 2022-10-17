subroutine get_eos_info_for_burn_at_const_density(&
      eos_handle, species, chem_id, net_iso, xa, &
      Rho, logRho, T, logT, &
      Cv, d_Cv_dlnT, eta, d_eta_dlnT, ierr)
   use const_def, only : dp
   integer, intent(in) :: eos_handle, species
   integer, pointer :: chem_id(:) ! maps species to chem id
   integer, pointer :: net_iso(:) ! maps chem id to species number
   real(dp), intent(in) :: &
         xa(:), rho, logRho, T, logT
   real(dp), intent(out) :: &
         Cv, d_Cv_dlnT, eta, d_eta_dlnT
   integer, intent(out) :: ierr
end subroutine get_eos_info_for_burn_at_const_density
