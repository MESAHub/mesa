module ideal

   use const_def, only: dp
   use math_lib
   use auto_diff
   use eos_def

   implicit none

   private
   public :: get_ideal_eos_results, get_ideal_alfa, get_ideal_for_eosdt

   contains

   subroutine get_ideal_alfa( & 
            rq, logRho, logT, Z, abar, zbar, &
            alfa, d_alfa_dlogT, d_alfa_dlogRho, &
            ierr)
      use const_def
      use eos_blend
      type (EoS_General_Info), pointer :: rq
      real(dp), intent(in) :: logRho, logT, Z, abar, zbar
      real(dp), intent(out) :: alfa, d_alfa_dlogT, d_alfa_dlogRho
      integer, intent(out) :: ierr

      alfa = 0d0
      d_alfa_dlogRho = 0d0
      d_alfa_dlogT = 0d0
   end subroutine get_ideal_alfa

   subroutine get_ideal_for_eosdt(handle, dbg, Z, X, abar, zbar, species, chem_id, net_iso, xa, &
                              rho, logRho, T, logT, remaining_fraction, res, d_dlnd, d_dlnT, d_dxa, skip, ierr)
   integer, intent(in) :: handle
   logical, intent(in) :: dbg
   real(dp), intent(in) :: &
      Z, X, abar, zbar, remaining_fraction
   integer, intent(in) :: species
   integer, pointer :: chem_id(:), net_iso(:)
   real(dp), intent(in) :: xa(:)
   real(dp), intent(in) :: rho, logRho, T, logT
   real(dp), intent(inout), dimension(nv) :: res, d_dlnd, d_dlnT
   real(dp), intent(inout), dimension(nv, species) :: d_dxa
   logical, intent(out) :: skip
   integer, intent(out) :: ierr
   type (EoS_General_Info), pointer :: rq

   rq => eos_handles(handle)

   call get_ideal_eos_results(rq, Z, X, abar, zbar, rho, logRho, T, logT, species, chem_id, xa, &
                              res, d_dlnd, d_dlnT, d_dxa, ierr)
   skip = .false.

   ! zero all components
   res(i_frac:i_frac+num_eos_frac_results-1) = 0.0
   d_dlnd(i_frac:i_frac+num_eos_frac_results-1) = 0.0
   d_dlnT(i_frac:i_frac+num_eos_frac_results-1) = 0.0

   ! mark this one
   res(i_frac_ideal) = 1.0

   end subroutine get_ideal_for_eosdt

   subroutine get_ideal_eos_results( &
         rq, Z, X, abar, zbar, Rho, logRho, T, logT, &
         species, chem_id, xa, res, d_dlnd, d_dlnT, d_dxa, ierr)   
   type (EoS_General_Info), pointer :: rq
   real(dp), intent(in) :: Z, X, abar, zbar
   real(dp), intent(in) :: Rho, logRho, T, logT
   integer, intent(in) :: species
   integer, pointer :: chem_id(:)
   real(dp), intent(in) :: xa(:)
   integer, intent(out) :: ierr
   real(dp), intent(out), dimension(nv) :: res, d_dlnd, d_dlnT
   real(dp), intent(out), dimension(nv, species) :: d_dxa
   
   real(dp) :: logT_ion, logT_neutral

   ierr = 0

   call ideal_eos( &
      T, Rho, X, abar, zbar, &
      species, chem_id, xa, &
      res, d_dlnd, d_dlnT, d_dxa, ierr)

   ! composition derivatives not provided
   d_dxa = 0

   end subroutine get_ideal_eos_results

   subroutine ideal_eos( &
         temp_in, den_in, Xfrac, abar, zbar,  &
         species, chem_id, xa, &
         res, d_dlnd, d_dlnT, d_dxa, ierr)

      use eos_def
      use const_def, only: dp
      use utils_lib, only: is_bad
      use chem_def, only: chem_isos
      use ion_offset, only: compute_ion_offset
      use skye_ideal
      use skye_thermodynamics
      use auto_diff

      implicit none
      integer :: j
      integer, intent(in) :: species
      integer, pointer :: chem_id(:)
      real(dp), intent(in) :: xa(:)
      real(dp), intent(in) :: temp_in, den_in
      real(dp), intent(in) :: Xfrac, abar, zbar
      integer, intent(out) :: ierr
      real(dp), intent(out), dimension(nv) :: res, d_dlnd, d_dlnT
      real(dp), intent(out), dimension(nv, species) :: d_dxa
      real(dp), parameter :: mass_fraction_limit = 1d-10
      
      integer :: relevant_species
      type(auto_diff_real_2var_order3) :: temp, den
      real(dp) :: ACMI(species), A(species), ya(species), select_xa(species), norm
      type(auto_diff_real_2var_order3) :: etaele, xnefer, phase, latent_ddlnT, latent_ddlnRho
      type(auto_diff_real_2var_order3) :: F_rad, F_ideal_ion, F_ele, F_coul

      ! Set up
      ierr = 0
      F_ele = 0d0
      F_coul = 0d0
      
      ! No electrons, so extreme negative chemical potential
      etaele = -1d42
      xnefer = 1d-99

      ! no latent heat
      latent_ddlnT = 0d0
      latent_ddlnRho = 0d0

      ! Construct rho,T and partials
      temp = temp_in
      temp%d1val1 = 1d0
      den = den_in
      den%d1val2 = 1d0      

      ! Count and pack relevant species for Coulomb corrections. Relevant means mass fraction above limit.
      relevant_species = 0
      norm = 0d0
      do j=1,species
         if (xa(j) > mass_fraction_limit) then
            relevant_species = relevant_species + 1
            ACMI(relevant_species) = chem_isos% W(chem_id(j))
            A(relevant_species) = chem_isos% Z_plus_N(chem_id(j))
            select_xa(relevant_species) = xa(j)
            norm = norm + xa(j)
         end if
      end do

      ! Normalize
      do j=1,relevant_species
         select_xa(j) = select_xa(j) / norm
      end do

      ! Compute number fractions
      norm = 0d0
      do j=1,relevant_species
         ya(j) = select_xa(j) / A(j)
         norm = norm + ya(j)
      end do
      do j=1,relevant_species
         ya(j) = ya(j) / norm
      end do

      ! Ideal ion free energy, only depends on abar
      F_ideal_ion = compute_F_ideal_ion(temp, den, abar, relevant_species, ACMI, ya)
      F_ideal_ion = F_ideal_ion + compute_ion_offset(relevant_species, select_xa, chem_id) ! Offset so ion ground state energy is zero.

      ! Radiation free energy, independent of composition
      F_rad = compute_F_rad(temp, den)


      call  pack_for_export(F_ideal_ion, F_coul, F_rad, F_ele, temp, den, xnefer, etaele, abar, zbar, &
                        phase, latent_ddlnT, latent_ddlnRho, res, d_dlnd, d_dlnT)

   end subroutine ideal_eos

end module ideal
