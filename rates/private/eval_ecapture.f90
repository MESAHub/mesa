! ***********************************************************************
!
!   Copyright (C) 2013-2019  Josiah Schwab, Bill Paxton & The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

module eval_ecapture

  use rates_def
  use const_def, only: dp, ln2
  use math_lib

  implicit none


contains

  subroutine do_eval_ecapture_reaction_info( &
       n, ids, cc, T9, YeRho, &
       etak, d_etak_dlnT, d_etak_dlnRho, &  ! these are kinetic chemical potentials
       lambda, dlambda_dlnT, dlambda_dlnRho, &
       Q, dQ_dlnT, dQ_dlnRho, &
       Qneu, dQneu_dlnT, dQneu_dlnRho, &
       ierr)
    use rates_def, only: Coulomb_Info
    use eval_psi
    use eval_coulomb
    use const_def, only: ln10, kerg, mev_to_ergs, keV, me, clight, ev2erg, fine
    use const_def, only: const_pi => pi
    use chem_def
    use utils_lib, only: is_bad, &
      integer_dict_lookup, integer_dict_size

    integer, intent(in) :: n, ids(:)
    type(Coulomb_Info), intent(in) :: cc
    real(dp), intent(in) :: T9, YeRho, etak, d_etak_dlnT, d_etak_dlnRho
    real(dp), dimension(:), intent(inout), pointer :: &
         lambda, dlambda_dlnT, dlambda_dlnRho, &
         Q, dQ_dlnT, dQ_dlnRho, &
         Qneu, dQneu_dlnT, dQneu_dlnRho
    integer, intent(out) :: ierr

    logical, parameter :: dbg = .false.

    integer :: i, ir, in, out, j, lhs, rhs
    integer :: offset, ntrans, lo, hi
    integer :: offset_lhs, nstates_lhs, lo_lhs, hi_lhs
    integer :: offset_rhs, nstates_rhs, lo_rhs, hi_rhs

    real(dp) :: mue, dmue_dlnRho, dmue_dlnT ! chem pot and derivatives

    character(len=iso_name_length) :: ecapture_lhs, ecapture_rhs

    real(dp) :: Qx, T, beta, kT, mec2, Z, Eavg, Ei, Ef, G_GM
    real(dp) :: eta, deta_dlnT, deta_dlnRho

    real(dp), dimension(max_num_ecapture_states) :: E_lhs, E_rhs, J_lhs, J_rhs
    real(dp), dimension(max_num_ecapture_states) :: bf, Pi, dPi_dlnT

    integer, dimension(max_num_ecapture_transitions) :: Si, Sf
    real(dp), dimension(max_num_ecapture_transitions) :: logft

    real(dp), dimension(max_num_ecapture_transitions) :: zeta
    real(dp), dimension(max_num_ecapture_transitions) :: Ie, dIe_dlnT, dIe_dlnRho
    real(dp), dimension(max_num_ecapture_transitions) :: Je, dJe_dlnT, dJe_dlnRho
    real(dp), dimension(max_num_ecapture_transitions) :: Pj, dPj_dlnT

    real(dp), dimension(max_num_ecapture_transitions) :: lambdaj, lambdaj_neu

    ! for coulomb corrections
    real(dp) :: Z_lhs, Z_rhs, muI_lhs, muI_rhs, delta_muI, delta_Q, Vs, f_muE

    real(dp) :: neutrino, dneutrino_dlnT, dneutrino_dlnRho, &
         gammaray, dgammaray_dlnT, dgammaray_dlnRho, lambdam1
    character(len=2*iso_name_length+1) :: key

    integer :: rxn_type, ii
    integer, parameter :: rxn_ecapture = 1, rxn_betadecay = -1

    include 'formats'

    ierr = 0

    ! first, translate the density, temperature, etc into appropriate units

    kT = 1d3 * keV * T9 ! in MeV
    mec2 = me * clight*clight / mev_to_ergs ! in MeV
    beta = mec2/kT ! dimesionless

    ! the chemical potentials from the equation of state are kinetic
    ! so add in the rest mass terms

    eta = etak + beta
    deta_dlnT = d_etak_dlnT - beta
    deta_dlnRho = d_etak_dlnRho

    ! only evaluate this for really degenerate stuff
    if (eta .lt. 2*beta) return

    ! also need chemical potential in MeV
    mue = eta * kT
    dmue_dlnRho = deta_dlnRho * kT
    dmue_dlnT = deta_dlnT * kT + mue

    do i = 1, n ! loop over reactions

       ! if there's not a weak reaction index, don't bother
       ir = ids(i)
       if (ir <= 0) then
          if (dbg) write(*,'(a,i3)') "No weak reaction for ", ir
          cycle
       endif

       ! get reactant names & ids
       ecapture_lhs = weak_lhs_nuclide_name(ir)
       ecapture_rhs = weak_rhs_nuclide_name(ir)

       ! generate reaction dictionary key (e.g. "mg24 na24")
       call create_ecapture_dict_key(ecapture_lhs, ecapture_rhs, key)

       ! get the transition data
       call integer_dict_lookup(ecapture_transitions_number_dict, key, ntrans, ierr)
       if (ierr /=0) then
          if (dbg) write(*,*) key, "is not a reaction included in ecapture module"
          ierr = 0
          cycle
       endif
       if (dbg) write(*,*) key, "is a reaction included in ecapture module"

       call integer_dict_lookup(ecapture_transitions_offset_dict, key, offset, ierr)
       if (ierr /=0) stop "ERROR: ecapture (transitions)"
       if (dbg) write(*,*) ntrans, offset

       ! get nuclide indicies
       lhs = weak_lhs_nuclide_id(ir)
       rhs = weak_rhs_nuclide_id(ir)

       ! get species charges (which will be one different)
       Z_lhs = chem_isos% Z(lhs)
       Z_rhs = chem_isos% Z(rhs)

       ! determine whether this is a beta-decay or an electron capture
       rxn_type = int(Z_lhs - Z_rhs)

       Qx = 0
       G_GM = 1
       select case (rxn_type)
       case(rxn_ecapture)
          ! use the *atomic* isotope data to get the *nuclear* mass excess
          Qx = chem_isos% mass_excess(lhs) - chem_isos% mass_excess(rhs) - mec2
          G_GM = exp(const_pi* fine * Z_lhs)
          if (dbg) write(*,*) "ecapture ", key
       case(rxn_betadecay)
          ! use the *atomic* isotope data to get the *nuclear* mass excess
          Qx = chem_isos% mass_excess(lhs) - chem_isos% mass_excess(rhs) + mec2
          G_GM = exp(const_pi * fine * Z_rhs)
          if (dbg) write(*,*) "betadecay ", key
       end select

       lo = offset + 1
       hi = offset + ntrans

       Si(1:ntrans) = ecapture_transitions_data(lo:hi, i_Si)
       Sf(1:ntrans) = ecapture_transitions_data(lo:hi, i_Sf)
       logft(1:ntrans) = ecapture_logft_data(lo:hi)

       ! get the left state info (energies and spins)

       call integer_dict_lookup(ecapture_states_number_dict, ecapture_lhs, nstates_lhs, ierr)
       if (ierr /=0) then
          write(*,*) 'ecapture_lhs ' // trim(ecapture_lhs)
          write(*,*) 'ecapture_rhs ' // trim(ecapture_rhs)
          write(*,*) 'size of dict', integer_dict_size(ecapture_states_number_dict)
          stop "ERROR: ecapture_states_number_dict (lhs states)"
       end if
       call integer_dict_lookup(ecapture_states_offset_dict, ecapture_lhs, offset_lhs, ierr)
       if (ierr /=0) stop "ERROR: ecapture_states_offset_dict (lhs states)"

       lo_lhs = offset_lhs + 1
       hi_lhs = offset_lhs + nstates_lhs

       E_lhs(1:nstates_lhs) = ecapture_states_data(lo_lhs:hi_lhs, i_E)
       J_lhs(1:nstates_lhs) = ecapture_states_data(lo_lhs:hi_lhs, i_J)

       ! get the right state info (energies and spins)

       call integer_dict_lookup(ecapture_states_number_dict, ecapture_rhs, nstates_rhs, ierr)
       if (ierr /=0) stop "ERROR: ecapture_states_number_dict (rhs states)"
       call integer_dict_lookup(ecapture_states_offset_dict, ecapture_rhs, offset_rhs, ierr)
       if (ierr /=0) stop "ERROR: ecapture_states_offset_dict (rhs states)"

       lo_rhs = offset_rhs + 1
       hi_rhs = offset_rhs + nstates_rhs

       E_rhs(1:nstates_rhs) = ecapture_states_data(lo_rhs:hi_rhs, i_E)
       J_rhs(1:nstates_rhs) = ecapture_states_data(lo_rhs:hi_rhs, i_J)

       ! we assume the left hand side states have thermal occupation fractions

       ! calculate boltzmann factor (bf), partition function (Z), and <E>
       do ii=1,nstates_lhs
         bf(ii) = (2 * J_lhs(ii) + 1) * exp(-E_lhs(ii)/ kT)
       end do
       Z = sum(bf(1:nstates_lhs))
       Eavg = sum(bf(1:nstates_lhs)* E_lhs(1:nstates_lhs))/ Z

       ! occupation probability (& derivative)
       Pi(1:nstates_lhs) = bf(1:nstates_lhs)/Z
       dPi_dlnT(1:nstates_lhs) = Pi(1:nstates_lhs) * (E_lhs(1:nstates_lhs) - Eavg) / kT

       ! now rearrange these rates to they apply to the transitions
       do j = 1, ntrans
          Pj(j) = Pi(Si(j))
          dPj_dlnT(j) = dPi_dlnT(Si(j))
       end do

       if (dbg) then
          do j = 1, ntrans
             write(*,"(2(F5.2, I2, F5.2))")  E_lhs(Si(j)), int(J_lhs(Si(j))), Pj(j), &
             E_rhs(Sf(j)), int(J_rhs(Sf(j))), logft(j)
          end do
       end if

       ! get chemical potential difference, related to coulomb correction

       ! get ion chemical potentials (already in units of kT)
       muI_lhs = do_muI_coulomb(cc, Z_lhs)
       muI_rhs = do_muI_coulomb(cc, Z_rhs)

       ! shift should be negative (for e-capture) given our definitions
       delta_muI = muI_lhs - muI_rhs
       delta_Q = delta_muI * kT

       ! get screening potential (in fraction of fermi energy)
       Vs = do_Vs_coulomb(cc, Z_lhs)
       f_muE = 1d0 - Vs

       if (dbg) write(*,*) eta, muI_lhs, muI_rhs, delta_muI, Vs

       ! do the phase space integrals
       do j = 1, ntrans

          Ei = E_lhs(Si(j))
          Ef = E_rhs(Sf(j))

          ! calculate energy difference, including electron rest mass
          zeta(j) = (Qx + delta_Q - Ef + Ei)/ kT

          select case(rxn_type)
          case(rxn_ecapture)
             call do_psi_Iec_and_Jec(beta, zeta(j), &
                  f_muE * eta, f_muE * deta_dlnT, f_muE * deta_dlnRho, &
                  Ie(j), dIe_dlnT(j), dIe_dlnRho(j), Je(j), dJe_dlnT(j), dJe_dlnRho(j))
          case(rxn_betadecay)
             call do_psi_Iee_and_Jee(beta, zeta(j), &
                  f_muE * eta, f_muE * deta_dlnT, f_muE * deta_dlnRho, &
                  Ie(j), dIe_dlnT(j), dIe_dlnRho(j), Je(j), dJe_dlnT(j), dJe_dlnRho(j))
          end select

       end do

       ! apply fermi-function correction factor
       Ie(1:ntrans) = G_GM * Ie(1:ntrans)
       dIe_dlnT(1:ntrans) = G_GM * dIe_dlnT(1:ntrans)
       dIe_dlnRho(1:ntrans) = G_GM * dIe_dlnRho(1:ntrans)

       Je(1:ntrans) = G_GM * Je(1:ntrans)
       dJe_dlnT(1:ntrans) = G_GM * dJe_dlnT(1:ntrans)
       dJe_dlnRho(1:ntrans) = G_GM * dJe_dlnRho(1:ntrans)

       ! add approximate fermi-function corrections

       do j=1,ntrans
         lambdaj(j) =  Ie(j) * exp10(-logft(j)) * ln2
         lambdaj_neu(j) = Je(j) * exp10(-logft(j)) * ln2
       end do

       ! masked arrays prevent dividing by zero in derivatives
       ! if Ie is 0, we also expect dIe_dlnT, dIe_dlnRho is 0

       lambda(i) = sum(lambdaj(1:ntrans) * Pj(1:ntrans))
       dlambda_dlnT(i) = sum(lambdaj(1:ntrans) *  &
            ( dPj_dlnT(1:ntrans) + Pj(1:ntrans) * dIe_dlnT(1:ntrans)/Ie(1:ntrans)), &
            MASK=(Ie(1:ntrans).ne.0))
       dlambda_dlnRho(i) = sum(lambdaj(1:ntrans) * &
            Pj(1:ntrans) * dIe_dlnRho(1:ntrans)/Ie(1:ntrans), &
            MASK=(Ie(1:ntrans).ne.0))

       neutrino = mec2 * sum(lambdaj_neu(1:ntrans) * Pj(1:ntrans))
       dneutrino_dlnT = mec2 * sum(lambdaj_neu(1:ntrans) * &
            ( dPj_dlnT(1:ntrans) + Pj(1:ntrans) * dJe_dlnT(1:ntrans)/Je(1:ntrans)), &
            MASK=(Je(1:ntrans).ne.0))
       dneutrino_dlnRho = mec2 * sum(lambdaj_neu(1:ntrans) * &
            Pj(1:ntrans) * dJe_dlnRho(1:ntrans)/Je(1:ntrans), &
            MASK=(Je(1:ntrans).ne.0))

       if (lambda(i) .gt. 1d-30) then
          lambdam1 = 1.0_dp / lambda(i)
       else
          lambdam1 = 0
       end if

       Qneu(i) = neutrino * lambdam1
       dQneu_dlnT(i) = (dneutrino_dlnT - dlambda_dlnT(i) * Qneu(i)) * lambdam1
       dQneu_dlnRho(i) = (dneutrino_dlnRho - dlambda_dlnRho(i) * Qneu(i)) * lambdam1

       ! this is the *total* energy per decay (neu losses are subtracted later)

       ! in the past, these Q values used to include terms associated
       ! with the electron and ion chemical potentials.
       ! these terms are now handled elsewhere, so Q is just the change in rest mass.
       ! note that Qx is defined differently here than in eval_weak,
       ! which is why we explicitly add/subtract mec2.

       select case(rxn_type)
       case(rxn_ecapture)

          Q(i) = Qx + mec2
          dQ_dlnT(i) = 0
          dQ_dlnRho(i) = 0

       case(rxn_betadecay)

          Q(i) = Qx - mec2
          dQ_dlnT(i) = 0
          dQ_dlnRho(i) = 0

       end select

    end do

    ierr = 0

  end subroutine do_eval_ecapture_reaction_info

end module eval_ecapture
