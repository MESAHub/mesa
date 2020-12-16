! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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
!
! ***********************************************************************

      module eos_HELM_eval
      use eos_def
      use const_def, only: avo, crad, ln10, arg_not_provided, mp, kerg, dp, qp
      use utils_lib, only: is_bad, mesa_error
      use math_lib
      use eospc_eval
      use helm

      implicit none
      
      logical, parameter :: stop_for_is_bad = .false.
      logical, parameter :: dbg = .false.

      contains


      subroutine get_helm_for_eosdt( &
            handle, dbg, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            rho, logRho, T, logT, remaining_fraction, &
            res, d_dlnd, d_dlnT, d_dxa, &
            skip, ierr)
         use chem_lib, only: composition_info
         integer, intent(in) :: handle
         logical, intent(in) :: dbg
         real(dp), intent(in) :: &
            Z, X, abar, zbar, remaining_fraction
         integer, intent(in) :: species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:)
         real(dp), intent(in) :: rho, logRho, T, logT
         real(dp), intent(inout), dimension(nv) :: &
            res, d_dlnd, d_dlnT
         real(dp), intent(inout), dimension(nv, species) :: d_dxa
         logical, intent(out) :: skip
         integer, intent(out) :: ierr
         type (EoS_General_Info), pointer :: rq
         real(dp), dimension(nv) :: d_dabar, d_dzbar
         real(dp) :: helm_res(num_helm_results)

         real(dp) :: xh, xhe, zz, abar_ci, zbar_ci, z2bar_ci, z53bar_ci, ye_ci, mass_correction, sumx
         real(dp), dimension(species) :: dabar_dx, dzbar_dx, dmc_dx
         integer :: i

         rq => eos_handles(handle)
         call Get_HELMEOS_Results( &
            rq, Z, X, abar, zbar, Rho, logRho, T, logT, &
            res, d_dlnd, d_dlnT, d_dabar, d_dzbar, &
            helm_res, skip, ierr)

         ! need this call to get dabar_dx, dzbar_dx
         ! might want to pass these in eventually
         call composition_info( &
            species, chem_id, xa, xh, xhe, zz, &
            abar_ci, zbar_ci, z2bar_ci, z53bar_ci, ye_ci, mass_correction, &
            sumx, dabar_dx, dzbar_dx, dmc_dx)

         ! zero these for now
         d_dabar(i_phase:i_latent_ddlnRho) = 0d0
         d_dzbar(i_phase:i_latent_ddlnRho) = 0d0
         d_dabar(i_frac:i_frac+num_eos_frac_results-1) = 0d0
         d_dzbar(i_frac:i_frac+num_eos_frac_results-1) = 0d0

         do i=1, species
            d_dxa(:,i) = d_dabar(:)*dabar_dx(i) + d_dzbar(:)*dzbar_dx(i)
         end do

         ! zero phase information
         res(i_phase:i_latent_ddlnRho) = 0d0
         d_dlnT(i_phase:i_latent_ddlnRho) = 0d0
         d_dlnd(i_phase:i_latent_ddlnRho) = 0d0

         ! zero all components
         res(i_frac:i_frac+num_eos_frac_results-1) = 0.0
         d_dlnd(i_frac:i_frac+num_eos_frac_results-1) = 0.0
         d_dlnT(i_frac:i_frac+num_eos_frac_results-1) = 0.0

         ! mark this one
         res(i_frac_HELM) = 1.0

      end subroutine get_helm_for_eosdt

         
      subroutine Get_HELMEOS_Results( &
               rq, Z, X, abar, zbar, Rho, logRho, T, logT, &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar, &
               helm_res, off_table, ierr)   
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: Z, X, abar, zbar
         real(dp), intent(in) :: Rho, logRho, T, logT
         real(dp), intent(inout), dimension(nv) :: &
            res, d_dlnd, d_dlnT, d_dabar, d_dzbar
         real(dp), intent(inout) :: helm_res(num_helm_results)
         logical, intent(out) :: off_table
         integer, intent(out) :: ierr

         real(dp) :: dx_0, err, dvardx
         logical :: doing_d_dlnd
         integer :: i_var         
         
         logical, parameter :: clip_to_table_boundaries = .true.
         
         real(dp), dimension(nv) :: &
            res_1, d_dlnd_1, d_dlnT_1, d_dabar_1, d_dzbar_1, &
            res_2, d_dlnd_2, d_dlnT_2, d_dabar_2, d_dzbar_2 
         real(dp), dimension(num_helm_results) :: helm_res_1, helm_res_2
         real(dp) :: logT_ion, logT_neutral, alfa, beta, d_alfa_dlogT, &
            alfa0, d_alfa0_dlnT, A, d_alfa_dlnT, d_beta_dlnT, AA, dA_dlnT
         integer :: j
         logical :: always_skip_elec_pos, always_include_elec_pos, &
            include_radiation
         
         include 'formats.dek'

         ierr = 0
         off_table = .false.
         
         always_skip_elec_pos = rq% always_skip_elec_pos
         always_include_elec_pos = rq% always_include_elec_pos
         include_radiation = rq% include_radiation
         logT_ion = rq% logT_ion_HELM        
         logT_neutral = rq% logT_neutral_HELM
         
         call helmeos2( &
            T, logT, Rho, logRho, X, abar, zbar, &
            rq% coulomb_temp_cut_HELM, rq% coulomb_den_cut_HELM, &
            helm_res, clip_to_table_boundaries, include_radiation, &
            always_skip_elec_pos, always_include_elec_pos, &
            logT_ion, logT_neutral, off_table, ierr)
         if (off_table) return
         if (ierr /= 0) then
            if (dbg) then
               write(*,*) 'failed in helmeos2'
               write(*,1) 'logT_ion', logT_ion
               write(*,1) 'logT_neutral', logT_neutral
               write(*,1) 'T', T
               write(*,1) 'logT', logT
               write(*,1) 'Rho', Rho
               write(*,1) 'logRho', logRho
               write(*,1) 'abar', abar
               write(*,1) 'zbar', zbar
               write(*,1) 'X', X
               write(*,*) 'clip_to_table_boundaries', clip_to_table_boundaries
               write(*,*) 'include_radiation', include_radiation
               write(*,*) 'always_skip_elec_pos', always_skip_elec_pos
               write(*,*) 'always_include_elec_pos', always_include_elec_pos
               stop 'Get1_HELMEOS_Results'
            end if
            !write(*,*) 'failed in helmeos2'
            return
         end if
         call do_convert_helm_results( &
               helm_res, Z, X, abar, zbar, Rho, T, &
               res, d_dlnd, d_dlnT, d_dabar, d_dzbar, ierr)
         if (ierr /= 0) then
            if (dbg) write(*,*) 'failed in do_convert_helm_results'
            return
         end if         

      end subroutine Get_HELMEOS_Results


      subroutine do_convert_helm_results( &
               helm_res, Z, X, abar, zbar, Rho, T, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, d_dabar_c_TRho, d_dzbar_c_TRho, ierr)
         use helm
         use chem_def
         real(dp), intent(in) :: helm_res(num_helm_results)
         real(dp), intent(in) :: Z, X, abar, zbar, Rho, T
         real(dp), intent(inout) :: res(nv)
         real(dp), intent(inout) :: d_dlnRho_c_T(nv)
         real(dp), intent(inout) :: d_dlnT_c_Rho(nv)
         real(dp), intent(inout) :: d_dabar_c_TRho(nv)
         real(dp), intent(inout) :: d_dzbar_c_TRho(nv)
         integer, intent(out) :: ierr

         real(dp) :: mu, P, Pgas, energy, entropy, free_e, dse, dpe, dsp
         integer :: j, k, ci
         
         include 'formats.dek'
         
         ierr = 0
         
         if (.false. .and. eos_test_partials) then   
            eos_test_partials_val = helm_res(h_etot)
            eos_test_partials_dval_dx = helm_res(h_dea)
            write(*,1) 'logRho', log10(Rho)
            write(*,1) 'logT', log10(T)
            write(*,1) 'Rho', Rho
            write(*,1) 'T', T
            write(*,1) 'X', X
            write(*,1) 'Z', Z
            write(*,1) 'abar', abar
            write(*,1) 'zbar', zbar
            write(*,1) 'etot', helm_res(h_etot)
            write(*,1) 'detot_dT', helm_res(h_det)
            write(*,1) 'detot_dRho', helm_res(h_ded)
            write(*,1) 'detot_dabar', helm_res(h_dea)
            write(*,1) 'detot_dzbar', helm_res(h_dez)
            stop 'do_convert_helm_results'
         end if
         
         energy = helm_res(h_etot)
         entropy = helm_res(h_stot)
         P = helm_res(h_ptot)
         Pgas = helm_res(h_pgas)
         
         res(i_lnE) = log(energy)
         res(i_lnS) = log(entropy)
         res(i_lnPgas) = log(Pgas)

         res(i_grad_ad) = helm_res(h_nabad)
         res(i_chiRho) = helm_res(h_chid)
         res(i_chiT) = helm_res(h_chit)
         res(i_Cp) = helm_res(h_cp)
         res(i_Cv) = helm_res(h_cv)
         res(i_dE_dRho) = helm_res(h_ded)
         res(i_dS_dT) = helm_res(h_dst)
         res(i_dS_dRho) = helm_res(h_dsd)
         mu = abar/(1 + zbar)
         res(i_mu) = mu
         free_e = max(1d-99, helm_res(h_xne))/(avo*Rho) ! assuming complete ionization
         res(i_lnfree_e) = log(free_e)
         res(i_gamma1) = helm_res(h_gam1)
         res(i_gamma3) = helm_res(h_gam3)
         res(i_eta) = helm_res(h_etaele)
         
         d_dlnRho_c_T(i_lnS) = helm_res(h_dsd)*Rho/entropy
         d_dlnRho_c_T(i_lnPgas) = helm_res(h_dpgasd)*Rho/Pgas
         d_dlnRho_c_T(i_lnE) = helm_res(h_ded)*Rho/energy
         
         d_dlnRho_c_T(i_grad_ad) = helm_res(h_dnabdd)*Rho
         d_dlnRho_c_T(i_chiRho) = helm_res(h_dchiddd)*Rho
         d_dlnRho_c_T(i_chiT) = helm_res(h_dchitdd)*Rho
         d_dlnRho_c_T(i_Cp) = helm_res(h_dcpdd)*Rho
         d_dlnRho_c_T(i_Cv) = helm_res(h_dcvdd)*Rho
         d_dlnRho_c_T(i_dE_dRho) = helm_res(h_dedd)*Rho
         d_dlnRho_c_T(i_dS_dT) = helm_res(h_dsdt)*Rho
         d_dlnRho_c_T(i_dS_dRho) = helm_res(h_dsdd)*Rho
         d_dlnRho_c_T(i_mu) = 0
         d_dlnRho_c_T(i_lnfree_e) = helm_res(h_dxned)/(avo*free_e) - 1d0
         d_dlnRho_c_T(i_gamma1) = helm_res(h_dgam1dd)*Rho
         d_dlnRho_c_T(i_gamma3) = helm_res(h_dgam3dd)*Rho
         d_dlnRho_c_T(i_eta) = helm_res(h_detad)*Rho
           
         d_dlnT_c_Rho(i_lnS) = helm_res(h_dst)*T/entropy
         d_dlnT_c_Rho(i_lnPgas) = helm_res(h_dpgast)*T/Pgas
         d_dlnT_c_Rho(i_lnE) = helm_res(h_det)*T/energy
         
         d_dlnT_c_Rho(i_grad_ad) = helm_res(h_dnabdt)*T
         d_dlnT_c_Rho(i_chiRho) = helm_res(h_dchiddt)*T
         d_dlnT_c_Rho(i_chiT) = helm_res(h_dchitdt)*T
         d_dlnT_c_Rho(i_Cp) = helm_res(h_dcpdt)*T
         d_dlnT_c_Rho(i_Cv) = helm_res(h_dcvdt)*T
         d_dlnT_c_Rho(i_dE_dRho) = helm_res(h_dedt)*T
         d_dlnT_c_Rho(i_dS_dT) = helm_res(h_dstt)*T
         d_dlnT_c_Rho(i_dS_dRho) = helm_res(h_dsdt)*T
         d_dlnT_c_Rho(i_mu) = 0
         d_dlnT_c_Rho(i_lnfree_e) = (helm_res(h_dxnet)*T/(avo*Rho))/free_e
         d_dlnT_c_Rho(i_gamma1) = helm_res(h_dgam1dt)*T
         d_dlnT_c_Rho(i_gamma3) = helm_res(h_dgam3dt)*T
         d_dlnT_c_Rho(i_eta) = helm_res(h_detat)*T

         if (.false.) then ! testing using mu
            ! pion ok.  pele ok for logRho large.  pcoul ok.
            ! egas = eion + eele + epos + ecoul
               ! ecoul ok. eion ok. eele bad.
            ! check sele
            res(i_mu) = helm_res(h_sele)
            d_dlnRho_c_T(i_mu) = helm_res(h_dsepd)*Rho
            d_dlnT_c_Rho(i_mu) = helm_res(h_dsept)*T
         end if
         
         d_dlnRho_c_T(i_lnE) = helm_res(h_ded)*Rho/energy
         d_dlnT_c_Rho(i_lnE) = helm_res(h_det)*T/energy         
         
         d_dlnRho_c_T(i_lnS) = helm_res(h_dsd)*Rho/entropy
         d_dlnT_c_Rho(i_lnS) = helm_res(h_dst)*T/entropy

         d_dlnRho_c_T(i_lnPgas) = helm_res(h_dpgasd)*Rho/Pgas
         d_dlnT_c_Rho(i_lnPgas) = helm_res(h_dpgast)*T/Pgas

         ! composition partials need set
         d_dabar_c_TRho(i_lnS) = helm_res(h_dsa)/entropy
         d_dabar_c_TRho(i_lnPgas) = helm_res(h_dpgasa)/Pgas
         d_dabar_c_TRho(i_lnE) = helm_res(h_dea)/energy
         
         d_dabar_c_TRho(i_grad_ad) = helm_res(h_dnabda)
         d_dabar_c_TRho(i_chiRho) = helm_res(h_dchidda)
         d_dabar_c_TRho(i_chiT) = helm_res(h_dchitda)
         d_dabar_c_TRho(i_Cp) = helm_res(h_dcpda)
         d_dabar_c_TRho(i_Cv) = helm_res(h_dcvda)
         d_dabar_c_TRho(i_dE_dRho) = helm_res(h_deda)
         d_dabar_c_TRho(i_dS_dT) = helm_res(h_dsta)
         d_dabar_c_TRho(i_dS_dRho) = helm_res(h_dsda)
         d_dabar_c_TRho(i_mu) = 0
         d_dabar_c_TRho(i_lnfree_e) = (helm_res(h_dxnea)/(avo*Rho))/free_e
         d_dabar_c_TRho(i_gamma1) = helm_res(h_dgam1da)
         d_dabar_c_TRho(i_gamma3) = helm_res(h_dgam3da)
         d_dabar_c_TRho(i_eta) = helm_res(h_detaa)

         d_dzbar_c_TRho(i_lnS) = helm_res(h_dsz)/entropy
         d_dzbar_c_TRho(i_lnPgas) = helm_res(h_dpgasz)/Pgas
         d_dzbar_c_TRho(i_lnE) = helm_res(h_dez)/energy
         
         d_dzbar_c_TRho(i_grad_ad) = helm_res(h_dnabdz)
         d_dzbar_c_TRho(i_chiRho) = helm_res(h_dchiddz)
         d_dzbar_c_TRho(i_chiT) = helm_res(h_dchitdz)
         d_dzbar_c_TRho(i_Cp) = helm_res(h_dcpdz)
         d_dzbar_c_TRho(i_Cv) = helm_res(h_dcvdz)
         d_dzbar_c_TRho(i_dE_dRho) = helm_res(h_dedz)
         d_dzbar_c_TRho(i_dS_dT) = helm_res(h_dstz)
         d_dzbar_c_TRho(i_dS_dRho) = helm_res(h_dsdz)
         d_dzbar_c_TRho(i_mu) = 0
         d_dzbar_c_TRho(i_lnfree_e) = (helm_res(h_dxnez)/(avo*Rho))/free_e
         d_dzbar_c_TRho(i_gamma1) = helm_res(h_dgam1dz)
         d_dzbar_c_TRho(i_gamma3) = helm_res(h_dgam3dz)
         d_dzbar_c_TRho(i_eta) = helm_res(h_detaz)

         
      end subroutine do_convert_helm_results

      
      real(dp) function eval_theta_e(eta, d_theta_e_deta)
         use interp_1d_lib, only : interp_values_and_slopes
         real(dp), intent(in) :: eta
         real(dp), intent(out) :: d_theta_e_deta
         integer, parameter :: nvar = 1
         integer :: ierr
         real(dp) :: x(nvar), vals(nvar), slopes(nvar)
         x(1) = eta
         call interp_values_and_slopes( &
                  x_theta_e, theta_e_nx, f_theta_e1, nvar, x, vals, slopes, ierr)
         eval_theta_e = vals(1)
         d_theta_e_deta = slopes(1)
      end function eval_theta_e 
      

      subroutine Get_HELM_Results( &
               X, abar, zbar, arho, alogrho, atemp, alogtemp, &
               coulomb_temp_cut, coulomb_den_cut, &
               include_radiation, always_skip_elec_pos, always_include_elec_pos, &
               logT_ion, logT_neutral, res, off_table, ierr)
         use const_def
         use helm

         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: X, abar, zbar
         real(dp), intent(in) :: arho, alogrho
         real(dp), intent(in) :: atemp, alogtemp 
         real(dp), intent(in) :: coulomb_temp_cut, coulomb_den_cut, & 
            logT_ion, logT_neutral
         logical, intent(in) :: include_radiation, &
            always_skip_elec_pos, always_include_elec_pos
         real(dp), intent(inout) :: res(:) ! (num_helm_results)
         logical, intent(out) :: off_table
         integer, intent(out) :: ierr ! 0 means AOK.

         real(dp) :: Rho, logRho, T, logT, dse, dpe, dsp
         
         logical, parameter :: clip_to_table_boundaries = .true.
         
         include 'formats.dek'
         
         ierr = 0
         off_table = .false.

         !..get temp and rho args
         T = atemp; logT = alogtemp
         if (atemp == arg_not_provided .and. alogtemp == arg_not_provided) then
            ierr = -2; return
         end if
         if (atemp == arg_not_provided) T = exp10(logT)
         
         Rho = arho; logrho = alogrho
         if (arho == arg_not_provided .and. alogrho == arg_not_provided) then
            ierr = -3; return
         end if
         if (arho == arg_not_provided) Rho = exp10(logRho)
         
         call helmeos2(T, logT, Rho, logRho, X, abar, zbar, &
                  coulomb_temp_cut, coulomb_den_cut, &
                  res, clip_to_table_boundaries, include_radiation, &
                  always_skip_elec_pos, always_include_elec_pos, &
                  logT_ion, logT_neutral, off_table, ierr)
         res(h_valid) = 1
         
      end subroutine Get_HELM_Results

         
      subroutine get_eosDT_HELMEOS_Results( &
               rq, Z_in, X_in, abar, zbar, &
               arho, alogrho, atemp, alogtemp, &
               include_radiation, always_skip_elec_pos, always_include_elec_pos, &
               logT_ion, logT_neutral, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, d_dabar_const_TRho, d_dzbar_const_TRho, &
               helm_res, off_table, ierr)
         use helm
         use const_def
         
         type (EoS_General_Info), pointer :: rq
         real(dp), intent(in) :: Z_in, X_in, abar, zbar
         real(dp), intent(in) :: arho, alogrho, atemp, alogtemp, &
            logT_ion, logT_neutral
         logical, intent(in) :: include_radiation, &
            always_skip_elec_pos, always_include_elec_pos
         real(dp), intent(inout) :: res(:) ! (nv)
         real(dp), intent(inout) :: d_dlnRho_c_T(:) ! (nv)
         real(dp), intent(inout) :: d_dlnT_c_Rho(:) ! (nv)
         real(dp), intent(inout) :: d_dabar_const_TRho(:) ! (num_eos_basic_results) 
         real(dp), intent(inout) :: d_dzbar_const_TRho(:) ! (num_eos_basic_results) 
         real(dp), intent(inout) :: helm_res(:) ! (num_helm_results)
         logical, intent(out) :: off_table
         integer, intent(out) :: ierr
         
         real(dp) :: Z, X, Rho, logRho, T, logT
         real(dp), parameter :: tiny = 1d-20
         integer :: ix, iz
         !logical, parameter :: clip_to_table_boundaries = .false.
         logical, parameter :: clip_to_table_boundaries = .true.
         
         include 'formats'
                  
         if (is_bad(X_in) .or. is_bad(Z_in)) then
            ierr = -1
            return
         end if
         
         X = X_in; Z = Z_in
         if (X < tiny) X = 0
         if (Z < tiny) Z = 0
         
         !..get temp and rho args
         T = atemp; logT = alogtemp
         if (atemp == arg_not_provided .and. alogtemp == arg_not_provided) then
            ierr = -2; return
         end if
         if (alogtemp == arg_not_provided) logT = log10(T)
         if (atemp == arg_not_provided) T = exp10(logT)
         
         if (T <= 0) then
            ierr = -1
            return
         end if
         
         Rho = arho; logrho = alogrho
         if (arho == arg_not_provided .and. alogrho == arg_not_provided) then
            ierr = -3; return
         end if
         if (alogrho == arg_not_provided) logRho = log10(Rho)
         if (arho == arg_not_provided) Rho = exp10(logRho)
         
         if (Rho <= 0) then
            ierr = -1
            return
         end if
         
         if (is_bad(Rho) .or. is_bad(T)) then
            ierr = -1
            return
         end if
         
         call helmeos2( &
               T, logT, Rho, logRho, X, abar, zbar, &
               rq% coulomb_temp_cut_HELM, rq% coulomb_den_cut_HELM, &
               helm_res, clip_to_table_boundaries, include_radiation, &
               always_skip_elec_pos, always_include_elec_pos, &
               logT_ion, logT_neutral, off_table, ierr)
         if (off_table) return
         if (ierr /= 0) then
            write(*,*) 'failed in helmeos2'
            return
         end if
         
         call do_convert_helm_results( &
               helm_res, Z, X, abar, zbar, Rho, T, &
               res, d_dlnRho_c_T, d_dlnT_c_Rho, d_dabar_const_TRho, d_dzbar_const_TRho, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in do_convert_helm_results'
            return
         end if
               
      end subroutine get_eosDT_HELMEOS_Results


      end module eos_HELM_eval
      
