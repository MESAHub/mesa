module skye
      use const_def, only: dp
      use math_lib
      use auto_diff
      use eos_def

      implicit none

      logical, parameter :: dbg = .false.
      !logical, parameter :: dbg = .true.
      
      
      private
      public :: Get_Skye_EOS_Results, Get_Skye_alfa, get_Skye_for_eosdt

      contains

      subroutine Get_Skye_alfa( & 
            rq, logRho, logT, Z, abar, zbar, &
            alfa, d_alfa_dlogT, d_alfa_dlogRho, &
            ierr)
         use const_def
         use eos_blend
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: logRho, logT, Z, abar, zbar
         real(dp), intent(out) :: alfa, d_alfa_dlogT, d_alfa_dlogRho
         integer, intent(out) :: ierr

         logical :: contained
         type(auto_diff_real_2var_order1) :: p(2), blend, dist

         ! Blend parameters
         real(dp) :: big
         real(dp) :: skye_blend_width
         integer, parameter :: num_points = 8
         real(dp) :: bounds(8,2)
         type (Helm_Table), pointer :: ht

         ierr = 0

         ht => eos_ht

         big = 12d0
         skye_blend_width = 0.1d0
         bounds(1,1) = ht% logdlo
         bounds(1,2) = 7.5d0

         bounds(2,1) = 4d0
         bounds(2,2) = 7.5d0

         bounds(3,1) = 0.6d0
         bounds(3,2) = 6.2d0

         bounds(4,1) = 4d0
         bounds(4,2) = 6.2d0

         bounds(5,1) = 4d0
         bounds(5,2) = ht% logtlo

         bounds(6,1) = ht% logdhi
         bounds(6,2) = ht% logtlo

         bounds(7,1) = ht% logdhi
         bounds(7,2) = ht% logthi

         bounds(8,1) = ht% logdlo
         bounds(8,2) = ht% logthi

         ! Set up auto_diff point
         p(1) = logRho
         p(1)%d1val1 = 1d0
         p(2) = logT
         p(2)%d1val2 = 1d0

         contained = is_contained(num_points, bounds, p)
         dist = min_distance_to_polygon(num_points, bounds, p)

         if (contained) then ! Make distance negative for points inside the polygon
            dist = -dist
         end if

         dist = dist / skye_blend_width
         blend = max(dist, 0d0)
         blend = min(blend, 1d0)

         alfa = blend%val
         d_alfa_dlogRho = blend%d1val1
         d_alfa_dlogT = blend%d1val2

      end subroutine Get_Skye_alfa


      subroutine get_Skye_for_eosdt(handle, dbg, Z, X, abar, zbar, species, chem_id, net_iso, xa, &
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

         call Get_Skye_EOS_Results(rq, Z, X, abar, zbar, rho, logRho, T, logT, species, chem_id, xa, &
                                    res, d_dlnd, d_dlnT, d_dxa, ierr)
         skip = .false.

         ! zero all components
         res(i_frac:i_frac+num_eos_frac_results-1) = 0.0
         d_dlnd(i_frac:i_frac+num_eos_frac_results-1) = 0.0
         d_dlnT(i_frac:i_frac+num_eos_frac_results-1) = 0.0

         ! mark this one
         res(i_frac_Skye) = 1.0

      end subroutine get_Skye_for_eosdt            

      subroutine Get_Skye_EOS_Results( &
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
         
         include 'formats'

         ierr = 0

         call skye_eos( &
            T, Rho, X, abar, zbar, &
            rq%Skye_min_gamma_for_solid, rq%Skye_max_gamma_for_liquid, &
            rq%Skye_solid_mixing_rule, rq%mass_fraction_limit_for_Skye, &
            species, chem_id, xa, &
            res, d_dlnd, d_dlnT, d_dxa, ierr)

         ! composition derivatives not provided
         d_dxa = 0

         if (ierr /= 0) then
            if (dbg) then
               write(*,*) 'failed in Get_Skye_EOS_Resultskye_EOS'
               write(*,1) 'T', T
               write(*,1) 'logT', logT
               write(*,1) 'Rho', Rho
               write(*,1) 'logRho', logRho
               write(*,1) 'abar', abar
               write(*,1) 'zbar', zbar
               write(*,1) 'X', X
               stop 'Get_Skye_EOS_Results'
            end if
            return
         end if     

      end subroutine Get_Skye_EOS_Results


      !>..given a temperature temp [K], density den [g/cm**3], and a composition 
      !!..this routine returns most of the other 
      !!..thermodynamic quantities. of prime interest is the pressure [erg/cm**3], 
      !!..specific thermal energy [erg/gr], the entropy [erg/g/K], along with 
      !!..their derivatives with respect to temperature, density, abar, and zbar.
      !!..other quantites such the normalized chemical potential eta (plus its
      !!..derivatives), number density of electrons and positron pair (along 
      !!..with their derivatives), adiabatic indices, specific heats, and 
      !!..relativistically correct sound speed are also returned.
      !!..
      !!..this routine assumes planckian photons, an ideal gas of ions, 
      !!..and an electron-positron gas with an arbitrary degree of relativity
      !!..and degeneracy. interpolation in a table of the helmholtz free energy
      !!..is used to return the electron-positron thermodynamic quantities.
      !!..all other derivatives are analytic.
      !!..
      !!..references: cox & giuli chapter 24 ; timmes & swesty apj 1999

      !!..this routine assumes a call to subroutine read_helm_table has
      !!..been performed prior to calling this routine.
      subroutine skye_eos( &
            temp_in, den_in, Xfrac, abar, zbar,  &
            Skye_min_gamma_for_solid, Skye_max_gamma_for_liquid, &
            Skye_solid_mixing_rule, &
            mass_fraction_limit, species, chem_id, xa, &
            res, d_dlnd, d_dlnT, d_dxa, ierr)

         use eos_def
         use const_def, only: dp
         use utils_lib, only: is_bad
         use chem_def, only: chem_isos
         use skye_ideal
         use skye_coulomb
         use skye_thermodynamics
         use auto_diff

         implicit none
         integer :: j
         integer, intent(in) :: species
         integer, pointer :: chem_id(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(in) :: temp_in, den_in, mass_fraction_limit, Skye_min_gamma_for_solid, Skye_max_gamma_for_liquid
         real(dp), intent(in) :: Xfrac, abar, zbar
         character(len=128), intent(in) :: Skye_solid_mixing_rule
         integer, intent(out) :: ierr
         real(dp), intent(out), dimension(nv) :: res, d_dlnd, d_dlnT
         real(dp), intent(out), dimension(nv, species) :: d_dxa
         
         integer :: relevant_species, lookup(species)
         type(auto_diff_real_2var_order3) :: temp, logtemp, den, logden, din
         real(dp) :: AZION(species), ACMI(species), select_xa(species), ya(species)
         type (Helm_Table), pointer :: ht
         real(dp) :: ytot1, ye, norm
         type(auto_diff_real_2var_order3) :: etaele, xnefer, phase, latent_ddlnT, latent_ddlnRho
         type(auto_diff_real_2var_order3) :: F_ion_gas, F_rad, F_ideal_ion, F_coul
         type(auto_diff_real_2var_order3) :: pele, eele, eep, sele

         ht => eos_ht

         ierr = 0

         temp = temp_in
         temp%d1val1 = 1d0
         logtemp = log10(temp)

         den = den_in
         den%d1val2 = 1d0
         logden = log10(den)

         ! HELM table lookup uses din rather than den
         ytot1 = 1.0d0 / abar
         ye = ytot1 * zbar
         din = ye*den

         F_rad = 0d0
         F_ion_gas = 0d0
         F_ideal_ion = 0d0
         F_coul = 0d0
         sele = 0d0
         pele = 0d0
         eele = 0d0

         ! Radiation free energy, independent of composition
         F_rad = compute_F_rad(temp, den)

         ! Count and pack relevant species for Coulomb corrections. Relevant means mass fraction above limit.
         relevant_species = 0
         norm = 0d0
         do j=1,species
            if (xa(j) > mass_fraction_limit) then
               relevant_species = relevant_species + 1
               AZION(relevant_species) = chem_isos% Z(chem_id(j))
               ACMI(relevant_species) = chem_isos% W(chem_id(j))
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
            ya(j) = select_xa(j) / ACMI(j)
            norm = norm + ya(j)
         end do
         do j=1,relevant_species
            ya(j) = ya(j) / norm
         end do

         ! Ideal ion free energy, only depends on abar
         F_ideal_ion = compute_F_ideal_ion(temp, den, abar, relevant_species, ACMI, ya)

         ! Ideal electron-positron thermodynamics (s, e, p)
         ! Derivatives are handled by HELM code, so we don't pass *in* any auto_diff types (just get them as return values).
         call compute_ideal_ele(temp%val, den%val, din%val, logtemp%val, logden%val, zbar, ytot1, ye, ht, &
                               sele, eele, pele, etaele, xnefer, ierr)

         xnefer = compute_xne(den, ytot1, zbar)

         ! Normalize mass fractions
         do j=1,relevant_species
            select_xa(j) = select_xa(j) / norm
         end do

         ! Compute non-ideal corrections
         call nonideal_corrections(relevant_species, ya(1:relevant_species), &
                                     AZION(1:relevant_species), ACMI(1:relevant_species), &
                                     Skye_min_gamma_for_solid, Skye_max_gamma_for_liquid, &
                                     Skye_solid_mixing_rule, den, temp, xnefer, abar, &
                                     F_coul, latent_ddlnT, latent_ddlnRho, phase)


         call  pack_for_export(F_ideal_ion, F_coul, F_rad, temp, den, xnefer, etaele, abar, zbar, &
                                 pele, eele, sele, phase, latent_ddlnT, latent_ddlnRho, &
                                 res, d_dlnd, d_dlnT)

      end subroutine skye_eos


end module skye
