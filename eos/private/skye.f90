module skye
      use const_def, only: dp
      use math_lib
      use auto_diff
      use eos_def

      implicit none

      logical, parameter :: dbg = .false.
      !logical, parameter :: dbg = .true.
      
      
      private
      public :: Get_Skye_EOS_Results, Get_Skye_alfa, Get_Skye_alfa_simple, get_Skye_for_eosdt

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

         ! Avoid catastrophic loss of precision in HELM tables
         bounds(1,1) = ht% logdlo
         bounds(1,2) = 8.3d0 

         ! Rough ionization temperature from Jermyn+2021 Equation 52 (treating denominator as ~1).
         ! We put a lower bound of logT=7.3 to ensure that solar models never use Skye.
         ! This is because the blend even in regions that are 99+% ionized produces noticeable
         ! kinks in the sound speed profile on a scale testable by the observations.
         bounds(2,1) = ht% logdlo
         bounds(2,2) = max(7.3d0,log10(1d5 * pow2(zbar))) + skye_blend_width

         ! Rough ionization density from Jermyn+2021 Equation 53, dividing by 3 so we get closer to Dragons.
         bounds(3,1) = max(2d0,log10(abar * pow3(zbar))) + skye_blend_width
         bounds(3,2) = max(7.3d0,log10(1d5 * pow2(zbar))) + skye_blend_width

         ! HELM low-T bound
         bounds(4,1) = max(2d0,log10(abar * pow3(zbar))) + skye_blend_width
         bounds(4,2) = ht% logtlo

         ! Lower-right of (rho,T) plane
         bounds(5,1) = ht% logdhi
         bounds(5,2) = ht% logtlo

         ! Upper-right of (rho,T) plane
         bounds(6,1) = ht% logdhi
         bounds(6,2) = ht% logthi

         ! Avoid catastrophic loss of precision in HELM tables
         bounds(7,1) = 3d0 * ht% logthi + log10(abar * mp * crad / (3d0 * kerg * (zbar + 1d0))) - 6d0
         bounds(7,2) =  ht% logthi

         ! Avoid catastrophic loss of precision in HELM tables
         bounds(8,1) = 3d0 * 8.3d0 + log10(abar * mp * crad / (3d0 * kerg * (zbar + 1d0))) - 6d0
         bounds(8,2) = 8.3d0

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


      subroutine Get_Skye_alfa_simple( &
            rq, logRho, logT, Z, abar, zbar, &
            alfa, d_alfa_dlogT, d_alfa_dlogRho, &
            ierr)
         use const_def
         use eos_blend
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: logRho, logT, Z, abar, zbar
         real(dp), intent(out) :: alfa, d_alfa_dlogT, d_alfa_dlogRho
         integer, intent(out) :: ierr

         type(auto_diff_real_2var_order1) :: logT_auto, logRho_auto
         type(auto_diff_real_2var_order1) :: blend, blend_logT, blend_logRho

         include 'formats'

         ierr = 0

         ! logRho is val1
         logRho_auto% val = logRho
         logRho_auto% d1val1 = 1d0
         logRho_auto% d1val2 = 0d0

         ! logT is val2
         logT_auto% val = logT
         logT_auto% d1val1 = 0d0
         logT_auto% d1val2 = 1d0

         ! logT blend
         if (logT_auto < rq% logT_min_for_any_Skye) then
            blend_logT = 0d0
         else if (logT_auto <= rq% logT_min_for_all_Skye) then
            blend_logT = (logT_auto - rQ% logT_min_for_any_Skye) / (rq% logT_min_for_all_Skye - rq% logT_min_for_any_Skye)
         else if (logT_auto > rq% logT_min_for_all_Skye) then
            blend_logT = 1d0
         end if


         ! logRho blend
         if (logRho_auto < rq% logRho_min_for_any_Skye) then
            blend_logRho = 0d0
         else if (logRho_auto <= rq% logRho_min_for_all_Skye) then
            blend_logRho = (logRho_auto - rQ% logRho_min_for_any_Skye) / (rq% logRho_min_for_all_Skye - rq% logRho_min_for_any_Skye)
         else if (logRho_auto > rq% logRho_min_for_all_Skye) then
            blend_logRho = 1d0
         end if

         ! combine blends
         blend = (1d0 - blend_logRho) * (1d0 - blend_logT)

         alfa = blend% val
         d_alfa_dlogRho = blend% d1val1
         d_alfa_dlogT = blend% d1val2

      end subroutine get_Skye_alfa_simple


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
               write(*,*) 'failed in Get_Skye_EOS_Results'
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
         use ion_offset, only: compute_ion_offset
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
         real(dp) :: AZION(species), ACMI(species), A(species), select_xa(species), ya(species)
         type (Helm_Table), pointer :: ht
         real(dp) :: ytot1, ye, norm
         type(auto_diff_real_2var_order3) :: etaele, xnefer, phase, latent_ddlnT, latent_ddlnRho
         type(auto_diff_real_2var_order3) :: F_ion_gas, F_rad, F_ideal_ion, F_coul
         type(auto_diff_real_2var_order3) :: F_ele

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
         F_ele = 0d0

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

         ! Ideal electron-positron thermodynamics (s, e, p)
         ! Derivatives are handled by HELM code, so we don't pass *in* any auto_diff types (just get them as return values).
         call compute_ideal_ele(temp%val, den%val, din%val, logtemp%val, logden%val, zbar, ytot1, ye, ht, &
                               F_ele, etaele, xnefer, ierr)

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

         call  pack_for_export(F_ideal_ion, F_coul, F_rad, F_ele, temp, den, xnefer, etaele, abar, zbar, &
                                 phase, latent_ddlnT, latent_ddlnRho, res, d_dlnd, d_dlnT)

      end subroutine skye_eos


end module skye
