! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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

module net_eval
   
   use const_def
   use math_lib
   use chem_def
   use chem_lib, only : get_mass_excess
   use net_def, only : Net_General_Info, Net_Info
   use utils_lib, only : fill_with_NaNs
   
   implicit none


contains
   
   
   subroutine eval_net(&
      n, g, rates_only, just_dxdt, &
      num_isos, num_reactions, num_wk_reactions, &
      x, temp, logtemp, rho, logrho, &
      abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
      rate_factors, weak_rate_factor, &
      reaction_Qs, reaction_neuQs, &
      eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, &
      dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, &
      screening_mode, &
      eps_nuc_categories, eps_neu_total, &
      actual_Qs, actual_neuQs, from_weaklib, symbolic, ierr)
      use net_initialize, only : &
         setup_net_info, set_ptrs_for_approx21
      use net_approx21, only : num_reactions_func => num_reactions
      use net_screen
      use net_derivs
      use net_def, only : &
         net_test_partials, &
         net_test_partials_val, net_test_partials_dval_dx, &
         net_test_partials_i, net_test_partials_iother
      
      type (Net_Info) :: n
      type (Net_General_Info), pointer :: g
      logical, intent(in) :: rates_only, just_dxdt
      integer, intent(in) :: num_isos
      integer, intent(in) :: num_reactions, num_wk_reactions
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: temp, logtemp
      real(dp), intent(in) :: rho, logrho
      real(dp), intent(in) :: abar  ! mean number of nucleons per nucleus
      real(dp), intent(in) :: zbar  ! mean charge per nucleus
      real(dp), intent(in) :: z2bar ! mean charge squared per nucleus
      real(dp), intent(in) :: ye
      real(dp), intent(in) :: eta, d_eta_dlnT, d_eta_dlnRho ! eta and derivatives
      real(dp), intent(in) :: rate_factors(:)
      real(dp), intent(in) :: weak_rate_factor
      real(dp), pointer, intent(in) :: reaction_Qs(:) ! (rates_reaction_id_max)
      real(dp), pointer, intent(in) :: reaction_neuQs(:) ! (rates_reaction_id_max)
      real(dp), intent(out) :: eps_nuc ! ergs/gram/second from burning
      real(dp), intent(out) :: d_eps_nuc_dT
      real(dp), intent(out) :: d_eps_nuc_dRho
      real(dp), intent(inout) :: d_eps_nuc_dx(:)
      real(dp), intent(inout) :: dxdt(:)
      real(dp), intent(inout) :: d_dxdt_dRho(:)
      real(dp), intent(inout) :: d_dxdt_dT(:)
      real(dp), intent(inout) :: d_dxdt_dx(:, :)
      real(dp), intent(inout) :: eps_nuc_categories(:)
      real(dp), intent(out) :: eps_neu_total
      integer, intent(in) :: screening_mode
      real(dp), pointer, dimension(:) :: actual_Qs, actual_neuQs ! ignore if null
      logical, pointer :: from_weaklib(:) ! ignore if null
      logical, intent(in) :: symbolic
      integer, intent(out) :: ierr
      
      integer, parameter :: max_z_for_cache = 14
      real(dp) :: enuc, T9, total, prev, curr, prev_T
      real(dp) :: eps_total, Ys, sum_dxdt, compare, Z_plus_N
      real(qp) :: eps_nuc_MeV(num_rvs)
      integer :: ci, i, j, ir, weak_id, h1, iwork
      integer(8) :: time0, time1
      logical :: doing_timing
      
      logical, parameter :: dbg = .false.
      !logical, parameter :: dbg = .true.
      
      include 'formats'
      
      if (dbg) write(*, *) 'enter eval_net'
      
      doing_timing = g% doing_timing
      if (doing_timing) then
         call system_clock(time0)
         g% doing_timing = .false.
      end if
      
      if (.not. g% net_has_been_defined) then
         ierr = -1
         if (dbg) write(*, *) 'failed (.not. g% net_has_been_defined)'
         return
      end if
      
      if (temp == arg_not_provided .or. logtemp == arg_not_provided .or. &
         rho == arg_not_provided .or. logrho == arg_not_provided) then
         write(*, *) "You must now eplxicity pass both the linear and log values of the temperature and density"
         ierr = -1
         return
      end if
      
      ierr = 0
      n% g => g
      
      if (dbg) write(*, *) 'call setup_net_info'
      call setup_net_info(n)
      
      n% reaction_Qs => reaction_Qs
      n% reaction_neuQs => reaction_neuQs
      n% weak_rate_factor = weak_rate_factor
      n% logT = logtemp
      n% temp = temp
      n% logRho = logrho
      n% rho = rho
      n% screening_mode = screening_mode
      n% x = x
      n% zbar = zbar
      n% abar = abar
      n% z2bar = z2bar
      n% ye = ye
      n% eta = eta
      n% d_eta_dlnT = d_eta_dlnT
      n% d_eta_dlnRho = d_eta_dlnRho
      n% rate_factors = rate_factors
      
      if (n% logT < rattab_tlo) then ! clip to table so can eval beta decays
         n% logT = rattab_tlo
         n% temp = rattab_temp_lo
      end if
      
      T9 = n% temp * 1d-9
      
      if (g% doing_approx21) then
         call set_ptrs_for_approx21(n)
      end if
      
      if (dbg) write(*, *) 'call set_molar_abundances'
      call set_molar_abundances(n, dbg, ierr)
      if (ierr /= 0) then
         if (dbg) write(*, *) 'failed in set_molar_abundances'
         return
      end if
      
      if (num_wk_reactions > 0) then
         if (dbg) write(*, *) 'call get_weaklib_rates'
         call get_weaklib_rates(n, ierr)
         if (ierr /= 0) then
            if (dbg) write(*, *) 'failed in get_weaklib_rates'
            return
         end if
      end if
      
      if (associated(actual_Qs) .and. associated(actual_neuQs)) then
         do i = 1, g% num_reactions
            ir = g% reaction_id(i)
            from_weaklib(i) = .false.
            actual_Qs(i) = n% reaction_Qs(ir)
            actual_neuQs(i) = n% reaction_neuQs(ir)
            weak_id = g% weak_reaction_index(i)
            if (weak_id > 0) then
               if (g% weaklib_ids(weak_id) > 0) then
                  from_weaklib(i) = .true.
                  actual_Qs(i) = n% Q(weak_id)
                  actual_neuQs(i) = n% Qneu(weak_id)
               end if
            end if
         end do
      end if
      
      if (doing_timing) then
         call system_clock(time1)
         g% clock_net_eval = g% clock_net_eval + (time1 - time0)
         time0 = time1
      end if
      
      if (dbg) write(*, *) 'call get_rates_with_screening'
      call get_rates_with_screening(n, ierr)
      if (dbg) write(*, *) 'done get_rates_with_screening'
      if (ierr /= 0) then
         if (dbg) write(*, *) 'failed in get_rates_with_screening'
         return
      end if
      
      if (rates_only) return
      
      ! n% d_eps_nuc_dT = 0
      ! n% d_eps_nuc_dRho = 0
      ! n% d_eps_nuc_dx(:) = 0
      
      ! n% dxdt(:) = 0
      ! n% d_dxdt_dRho(:) = 0
      ! n% d_dxdt_dT(:) = 0
      ! if (.not. just_dxdt) d_dxdt_dx(:,:) = 0
      n% eps_nuc_categories(:) = 0
      n% eps_neu_total = 0
      n% d_eps_nuc_dy = 0
      
      if (g% doing_approx21) then
         call eval_net_approx21_procs(n, just_dxdt, ierr)
         if (ierr /= 0) return
         
         if (net_test_partials) then
            net_test_partials_val = eps_nuc
            net_test_partials_dval_dx = d_eps_nuc_dx(net_test_partials_i)
            if (g% add_co56_to_approx21) then
               write(*, *) 'net: eval_net_approx21_plus_co56'
            else
               write(*, *) 'net: eval_net_approx21'
            end if
         end if
         
         call unpack_for_export(n, eps_nuc, d_eps_nuc_dT, d_eps_nuc_dRho, d_eps_nuc_dx, &
            eps_neu_total, &
            dxdt, d_dxdt_dT, d_dxdt_dRho, d_dxdt_dx, &
            eps_nuc_categories)
         
         return
      end if   ! End of approx21
      
      if (dbg) write(*, *) 'call get_derivs'
      call get_derivs(&
         n, n% dydt, eps_nuc_MeV(1:num_rvs), n% eta, n% ye, &
         n% logT, n% temp, n% rho, n% abar, n% zbar, &
         num_reactions, n% rate_factors, &
         symbolic, just_dxdt, ierr)
      if (ierr /= 0) then
         if (dbg) write(*, *) 'failed in get_derivs'
         return
      end if
      
      if (symbolic) then
         do j = 1, num_isos
            do i = 1, num_isos
               d_dxdt_dx(i, j) = n% d_dydt_dy(i, j)
            end do
         end do
         return
      end if
      
      if (doing_timing) then
         call system_clock(time1)
         g% clock_net_derivs = g% clock_net_derivs + (time1 - time0)
         time0 = time1
      end if
      
      ! convert the eps_nuc_categories
      do i = 1, num_categories
         n% eps_nuc_categories(i) = Qconv * n% eps_nuc_categories(i)
      end do
      
      ! store the results
      do i = 1, num_isos
         ci = g% chem_id(i)
         n% dxdt(i) = chem_isos% Z_plus_N(ci) * n% dydt(i_rate, i)
      end do
      
      if (.not. just_dxdt) call store_partials(n)
      
      n% eps_nuc = eps_nuc_MeV(i_rate) * Qconv
      n% d_eps_nuc_dT = eps_nuc_MeV(i_rate_dT) * Qconv
      n% d_eps_nuc_dRho = eps_nuc_MeV(i_rate_dRho) * Qconv
      
      n% eps_neu_total = n% eps_neu_total * Qconv
      
      if (doing_timing) then
         call system_clock(time1)
         g% clock_net_eval = g% clock_net_eval + (time1 - time0)
         g% doing_timing = .true.
      end if
      
      if (net_test_partials) then
         !net_test_partials_val = eps_nuc
         !net_test_partials_dval_dx = d_eps_nuc_dx(net_test_partials_i)
         net_test_partials_val = &
            n% rate_screened(g% net_reaction(irn14_to_c12)) / &
               n% rate_raw(g% net_reaction(irn14_to_c12))
         net_test_partials_dval_dx = 0d0
         write(*, *) 'net_test_partials'
      end if
      
      call unpack_for_export(n, eps_nuc, d_eps_nuc_dT, d_eps_nuc_dRho, d_eps_nuc_dx, &
         eps_neu_total, &
         dxdt, d_dxdt_dT, d_dxdt_dRho, d_dxdt_dx, &
         eps_nuc_categories)
   
   end subroutine eval_net
   
   subroutine unpack_for_export(n, eps_nuc, d_eps_nuc_dT, d_eps_nuc_dRho, d_eps_nuc_dx, &
      eps_neu_total, &
      dxdt, d_dxdt_dT, d_dxdt_dRho, d_dxdt_dx, &
      eps_nuc_categories)
      type(net_info) :: n
      real(dp), intent(out) :: eps_nuc ! ergs/gram/second from burning
      real(dp), intent(out) :: d_eps_nuc_dT
      real(dp), intent(out) :: d_eps_nuc_dRho
      real(dp), intent(out) :: d_eps_nuc_dx(:)
      real(dp), intent(out) :: dxdt(:)
      real(dp), intent(out) :: d_dxdt_dRho(:)
      real(dp), intent(out) :: d_dxdt_dT(:)
      real(dp), intent(out) :: d_dxdt_dx(:, :)
      real(dp), intent(out) :: eps_neu_total
      real(dp), intent(out) :: eps_nuc_categories(:)
      
      eps_nuc = n% eps_nuc
      d_eps_nuc_dT = n% d_eps_nuc_dT
      d_eps_nuc_dRho = n% d_eps_nuc_dRho
      d_eps_nuc_dx = n% d_eps_nuc_dx
      
      dxdt = n% dxdt
      d_dxdt_dT = n% d_dxdt_dT
      d_dxdt_dRho = n% d_dxdt_dRho
      d_dxdt_dx = n% d_dxdt_dx
      
      eps_neu_total = n% eps_neu_total
      
      eps_nuc_categories = n% eps_nuc_categories
   
   end subroutine unpack_for_export
   
   
   subroutine eval_net_approx21_procs(n, just_dxdt, ierr)
      use net_approx21
      use rates_def
      type(net_info) :: n
      type(net_general_info), pointer :: g => null()
      logical, intent(in) :: just_dxdt
      integer :: ierr
      
      integer :: ci, i, j, num_isos
      real(dp) :: Z_plus_N
      
      ierr = 0
      
      g => n% g
      
      num_isos = g% num_isos
      
      call approx21_special_reactions(&
         n% temp, n% rho, n% abar, n% zbar, n% y, &
         g% use_3a_fl87, Qconv * n% reaction_Qs(ir_he4_he4_he4_to_c12), &
         n% rate_screened, n% rate_screened_dT, n% rate_screened_dRho, &
         n% dratdumdy1, n% dratdumdy2, g% add_co56_to_approx21, ierr)
      if (ierr /= 0) return
      
      call approx21_dydt(&
         n% y, n% rate_screened, n% rate_screened, &
         n% dydt1, .false., g% fe56ec_fake_factor, g% min_T_for_fe56ec_fake_factor, &
         g% fe56ec_n_neut, n% temp, n% rho, g% add_co56_to_approx21, ierr)
      if (ierr /= 0) return
      
      n% fII = approx21_eval_PPII_fraction(n% y, n% rate_screened)
      
      call get_approx21_eps_info(n, &
         n% dydt1, n% rate_screened, .true., n% eps_total, n% eps_neu_total, &
         g% add_co56_to_approx21, ierr)
      
      if (ierr /= 0) return
      n% eps_nuc = n% eps_total - n% eps_neu_total
      
      do i = 1, num_isos
         n% dxdt(i) = chem_isos% Z_plus_N(g% chem_id(i)) * n% dydt1(i)
      end do
      
      if (just_dxdt) return
      
      call approx21_dfdy(&
         n% y, n% dfdy, &
         g% fe56ec_fake_factor, g% min_T_for_fe56ec_fake_factor, g% fe56ec_n_neut, &
         n% rate_screened, n% rate_screened_dT, n% rate_screened_dRho, &
         n% dratdumdy1, n% dratdumdy2, n% temp, g% add_co56_to_approx21, ierr)
      if (ierr /= 0) return
      
      call approx21_dfdT_dfdRho(&
         
         ! NOTE: currently this gives d_eps_total_dy -- should fix to account for neutrinos too
         
         n% y, g% mion, n% dfdy, n% rate_screened, n% rate_screened_dT, n% rate_screened_dRho, &
         g% fe56ec_fake_factor, g% min_T_for_fe56ec_fake_factor, &
         g% fe56ec_n_neut, n% temp, n% rho, n% dfdT, n% dfdRho, n% d_epsnuc_dy, g% add_co56_to_approx21, ierr)
      if (ierr /= 0) return
      
      call get_approx21_eps_info(n, &
         n% dfdT, n% rate_screened_dT, .false., n% deps_total_dT, n% deps_neu_dT, &
         g% add_co56_to_approx21, ierr)
      
      if (ierr /= 0) return
      n% d_eps_nuc_dT = n% deps_total_dT - n% deps_neu_dT
      
      call get_approx21_eps_info(n, &
         n% dfdRho, n% rate_screened_dRho, .false., n% deps_total_dRho, n% deps_neu_dRho, &
         g% add_co56_to_approx21, ierr)
      
      if (ierr /= 0) return
      n% d_eps_nuc_dRho = n% deps_total_dRho - n% deps_neu_dRho
      
      call approx21_d_epsneu_dy(&
         n% y, n% rate_screened, &
         n% reaction_neuQs(irpp_to_he3), &
         n% reaction_neuQs(ir34_pp2), &
         n% reaction_neuQs(ir34_pp3), &
         n% reaction_neuQs(irc12_to_n14), &
         n% reaction_neuQs(irn14_to_c12), &
         n% reaction_neuQs(iro16_to_n14), &
         n% d_epsneu_dy, &
         g% add_co56_to_approx21, ierr)
      if (ierr /= 0) return
      
      do i = 1, n% g%num_isos
         ci = g% chem_id(i)
         Z_plus_N = dble(chem_isos% Z_plus_N(ci))
         n% d_eps_nuc_dx(i) = (n% d_epsnuc_dy(i) - n% d_epsneu_dy(i)) / Z_plus_N
         n% d_dxdt_dRho(i) = Z_plus_N * n% dfdRho(i)
         n% d_dxdt_dT(i) = Z_plus_N * n% dfdT(i)
         do j = 1, num_isos
            n% d_dxdt_dx(i, j) = &
               n% dfdy(i, j) * Z_plus_N / chem_isos% Z_plus_N(g% chem_id(j))
         end do
      end do
   
   end subroutine eval_net_approx21_procs
   
   
   subroutine get_approx21_eps_info(n, &
      dydt1, rate_screened, do_eps_nuc_categories, eps_total, eps_neu_total, plus_co56, ierr)
      use net_approx21, only : approx21_eps_info
      use rates_def
      type(net_info) :: n
      type(net_general_info), pointer :: g => null()
      real(dp), intent(in), dimension(:) :: dydt1, rate_screened
      logical, intent(in) :: do_eps_nuc_categories
      real(dp), intent(out) :: eps_total, eps_neu_total
      logical, intent(in) :: plus_co56
      integer, intent(out) :: ierr
      real(dp) :: Qtotal_rfe56ec, Qneu_rfe56ec
      
      g => n% g
      
      ! Indexes into reaction_Qs and reaction_neuQs should be in terms of the
      ! normal rate ids not the approx21 rate ids (in net_approx21.f90)
      
      call get_Qs_rfe56ec(n, Qtotal_rfe56ec, Qneu_rfe56ec)
      
      call approx21_eps_info(&
         n, n% y, g% mion, dydt1, rate_screened, n% fII, &
         n% reaction_Qs(irpp_to_he3), n% reaction_neuQs(irpp_to_he3), &
         n% reaction_Qs(ir_he3_he3_to_h1_h1_he4), &
         n% reaction_Qs(ir34_pp2), n% reaction_neuQs(ir34_pp2), &
         n% reaction_Qs(ir34_pp3), n% reaction_neuQs(ir34_pp3), &
         n% reaction_Qs(irc12_to_n14), n% reaction_neuQs(irc12_to_n14), &
         n% reaction_Qs(irn14_to_c12), n% reaction_neuQs(irn14_to_c12), &
         n% reaction_Qs(iro16_to_n14), n% reaction_neuQs(iro16_to_n14), &
         n% reaction_Qs(irn14_to_o16), &
         
         n% reaction_Qs(irprot_to_neut), n% reaction_neuQs(irprot_to_neut), &
         n% reaction_Qs(irneut_to_prot), n% reaction_neuQs(irneut_to_prot), &
         n% reaction_Qs(irni56ec_to_co56), n% reaction_neuQs(irni56ec_to_co56), &
         n% reaction_Qs(irco56ec_to_fe56), n% reaction_neuQs(irco56ec_to_fe56), &
         Qtotal_rfe56ec, Qneu_rfe56ec, &
         
         n% reaction_Qs(irn14ag_lite), &
         n% reaction_Qs(ir_he4_he4_he4_to_c12), &
         n% reaction_Qs(ir_c12_ag_o16), n% reaction_Qs(ir_o16_ag_ne20), &
         n% reaction_Qs(ir1212), &
         n% reaction_Qs(ir1216_to_mg24), n% reaction_Qs(ir1216_to_si28), &
         n% reaction_Qs(ir1616a), n% reaction_Qs(ir1616g), &
         n% reaction_Qs(ir_ne20_ag_mg24), &
         n% reaction_Qs(ir_mg24_ag_si28), &
         n% reaction_Qs(ir_si28_ag_s32), &
         n% reaction_Qs(ir_s32_ag_ar36), &
         n% reaction_Qs(ir_ar36_ag_ca40), &
         n% reaction_Qs(ir_ca40_ag_ti44), &
         n% reaction_Qs(ir_ti44_ag_cr48), &
         n% reaction_Qs(ir_cr48_ag_fe52), &
         n% reaction_Qs(ir_fe52_ag_ni56), &
         n% reaction_Qs(ir_fe52_ng_fe53), &
         n% reaction_Qs(ir_fe53_ng_fe54), &
         n% reaction_Qs(ir_fe54_ng_fe55), &
         n% reaction_Qs(ir_fe55_ng_fe56), &
         n% reaction_Qs(irfe52neut_to_fe54), &
         n% reaction_Qs(irfe52aprot_to_fe54), &
         n% reaction_Qs(irfe54ng_to_fe56), &
         n% reaction_Qs(irfe54aprot_to_fe56), &
         n% reaction_Qs(irfe52aprot_to_ni56), &
         n% reaction_Qs(irfe54prot_to_ni56), &
         n% reaction_Qs(irhe4_breakup), &
         n% reaction_Qs(irhe4_rebuild), &
         eps_total, eps_neu_total, & ! Dont use n% here as we call this for both eps_neu and eps_neu_dt and drho
         do_eps_nuc_categories, n% eps_nuc_categories, &
         .false., plus_co56, ierr)
   
   end subroutine get_approx21_eps_info
   
   subroutine get_Qs_rfe56ec(n, Qtotal, Qneu)
      use chem_def
      use rates_def
      type(net_info) :: n
      real(dp), intent(out) :: Qtotal, Qneu
      integer :: id, ir
      include 'formats'
      id = n% g% approx21_ye_iso
      if (id == imn56) then
         ir = irfe56ec_fake_to_mn56
      else if (id == imn57) then
         ir = irfe56ec_fake_to_mn57
      else if (id == icr56) then
         ir = irfe56ec_fake_to_cr56
      else if (id == icr57) then
         ir = irfe56ec_fake_to_cr57
      else if (id == icr58) then
         ir = irfe56ec_fake_to_cr58
      else if (id == icr59) then
         ir = irfe56ec_fake_to_cr59
      else if (id == icr60) then
         ir = irfe56ec_fake_to_cr60
      else if (id == icr61) then
         ir = irfe56ec_fake_to_cr61
      else if (id == icr62) then
         ir = irfe56ec_fake_to_cr62
      else if (id == icr63) then
         ir = irfe56ec_fake_to_cr63
      else if (id == icr64) then
         ir = irfe56ec_fake_to_cr64
      else if (id == icr65) then
         ir = irfe56ec_fake_to_cr65
      else if (id == icr66) then
         ir = irfe56ec_fake_to_cr66
      else
         ir = irco56ec_to_fe56
      end if
      Qtotal = n% reaction_Qs(ir)
      Qneu = n% reaction_neuQs(ir)
   end subroutine get_Qs_rfe56ec
   
   subroutine store_partials(n)
      use rates_def, only : i_rate, i_rate_dT, i_rate_dRho
      type(net_info) :: n
      type(net_general_info), pointer :: g => null()
      integer :: i, j, ci
      real(dp) :: Z_plus_N
      include 'formats'
      g => n%g
      do i = 1, g% num_isos
         ci = g% chem_id(i)
         Z_plus_N = dble(chem_isos% Z_plus_N(ci))
         n% d_eps_nuc_dx(i) = Qconv * n% d_eps_nuc_dy(i) / Z_plus_N
         n% dxdt(i) = Z_plus_N * n% dydt(i_rate, i)
         n% d_dxdt_dRho(i) = Z_plus_N * n% dydt(i_rate_dRho, i)
         n% d_dxdt_dT(i) = Z_plus_N * n% dydt(i_rate_dT, i)
         do j = 1, g% num_isos
            n% d_dxdt_dx(i, j) = &
               n% d_dydt_dy(i, j) * Z_plus_N / chem_isos% Z_plus_N(g% chem_id(j))
         end do
      end do
   
   end subroutine store_partials
   
   subroutine get_rates_with_screening(n, ierr)
      use rates_def, only : reaction_inputs, nrattab
      use rates_lib, only : eval_using_rate_tables
      use net_approx21, only : num_reactions_func => num_reactions
      use net_screen
      
      type(net_info) :: n
      type(net_general_info), pointer :: g => null()
      
      integer, intent(out) :: ierr
      
      logical, parameter :: dbg = .false.
      integer(8) :: time0, time1
      
      integer :: i, num, num_reactions
      real(dp) :: f
      logical :: okay
      
      include 'formats'
      
      g => n% g
      
      do i = 1, g% num_reactions
         if (g% reaction_id(i) <= 0) then
            write(*, 2) 'g% reaction_id(i)', i, g% reaction_id(i)
            call mesa_error(__FILE__, __LINE__, 'get_rates_with_screening')
         end if
      end do
      
      if (dbg) write(*, *) 'call eval_using_rate_tables'
      call eval_using_rate_tables(&
         g% num_reactions, g% reaction_id, g% rate_table, g% rattab_f1, nrattab, &
         n% ye, n% logT, n% temp, n% rho, n% rate_factors, g% logttab, &
         n% rate_raw, n% rate_raw_dT, n% rate_raw_dRho, ierr)
      if (ierr /= 0) then
         if (dbg) write(*, *) 'ierr from eval_using_rate_tables'
         return
      end if
      
      if (g% doing_timing) then
         call system_clock(time0)
      end if
      
      if (g% doing_approx21) then
         call approx21_rates(n, g% add_co56_to_approx21, ierr)
         if (ierr /= 0) return
      end if
      
      ! get the reaction rates including screening factors
      if (dbg) write(*, *) 'call screen_net with init=.false.'
      call screen_net(&
         g, g% num_isos, n% y, n% temp, n% rho, n% logT, n% logrho, .false., &
         n% rate_raw, n% rate_raw_dT, n% rate_raw_dRho, &
         n% rate_screened, n% rate_screened_dT, n% rate_screened_dRho, &
         n% screening_mode, &
         n% zbar, n% abar, n% z2bar, n% ye, ierr)
      if (dbg) write(*, *) 'done screen_net with init=.false.'
      if (ierr /= 0) return
      if (g% doing_approx21) then
         num = num_reactions_func(g% add_co56_to_approx21)
         do i = g% num_reactions + 1, num
            n% rate_screened(i) = n% rate_raw(i)
            n% rate_screened_dT(i) = n% rate_raw_dT(i)
            n% rate_screened_dRho(i) = n% rate_raw_dRho(i)
         end do
         do i = 1, num
            n% dratdumdy1(i) = 0d0
            n% dratdumdy2(i) = 0d0
         end do
      end if
      
      if (g% doing_timing) then
         call system_clock(time1)
         g% clock_net_screen = g% clock_net_screen + (time1 - time0)
      end if
   
   end subroutine get_rates_with_screening
   
   subroutine approx21_rates(n, plus_co56, ierr)
      use net_approx21, only : &
         approx21_pa_pg_fractions, approx21_weak_rates
      type(net_info) :: n
      logical, intent(in) :: plus_co56
      integer, intent(out) :: ierr
      ierr = 0
      call approx21_pa_pg_fractions(&
         n% rate_raw, n% rate_raw_dT, n% rate_raw_dRho, ierr)
      if (ierr /= 0) return
      call approx21_weak_rates(&
         n% y, n% rate_raw, n% rate_raw_dT, n% rate_raw_dRho, &
         n% temp, n% rho, n% ye, n% eta, n% zbar, &
         n% weak_rate_factor, plus_co56, ierr)
      if (ierr /= 0) return
   end subroutine approx21_rates
   
   
   subroutine get_weaklib_rates(n, ierr)
      use rates_def, only : Coulomb_Info
      use rates_lib, only : eval_weak_reaction_info, coulomb_set_context
      use net_def, only : other_kind
      
      type (net_info) :: n
      type(net_general_info), pointer :: g
      
      type (Coulomb_Info), target :: cc_info
      type (Coulomb_Info), pointer :: cc
      
      integer, intent(out) :: ierr
      integer :: i, j, id, ir
      integer(8) :: time0, time1
      
      include 'formats'
      
      ! before getting the weaklib rates, the Coulomb_Info
      ! structure must be populated.  the ecapture routines need
      ! to know some local quantities (functions of the density,
      ! temperature, and composition), to calculate Coulomb
      ! corrections to the rates
      
      ierr = 0
      cc => cc_info
      
      g => n% g
      
      if(g% doing_timing) then
         call system_clock(time0)
      end if
      
      call coulomb_set_context(cc, n% temp, n% rho, n% logT, n% logRho, n% zbar, n% abar, n% z2bar)
      
      call eval_weak_reaction_info(&
         g% num_wk_reactions, &
         g% weaklib_ids(1:g% num_wk_reactions), &
         g% reaction_id_for_weak_reactions(1:g% num_wk_reactions), &
         cc, n% temp * 1d-9, n% ye * n% rho, n% eta, n% d_eta_dlnT, n% d_eta_dlnRho, &
         n% lambda, n% dlambda_dlnT, n% dlambda_dlnRho, &
         n% Q, n% dQ_dlnT, n% dQ_dlnRho, &
         n% Qneu, n% dQneu_dlnT, n% dQneu_dlnRho, &
         ierr)
      if (n% weak_rate_factor < 1d0) then
         do i = 1, g% num_wk_reactions
            n% lambda(i) = n% weak_rate_factor * n% lambda(i)
            n% dlambda_dlnT(i) = n% weak_rate_factor * n% dlambda_dlnT(i)
            n% dlambda_dlnRho(i) = n% weak_rate_factor * n% dlambda_dlnRho(i)
         end do
      end if
      if (g% doing_timing) then
         call system_clock(time1)
         g% clock_net_weak_rates = g% clock_net_weak_rates + (time1 - time0)
      end if
   
   end subroutine get_weaklib_rates
   
   
   subroutine get_T_limit_factor(&
      temp, lnT, T_lo, T_hi, lnT_lo, lnT_hi, &
      min_ln_factor, min_factor, &
      factor, d_factor_dT)
      real(dp), intent(in) :: &
         temp, lnT, T_lo, T_hi, lnT_lo, lnT_hi, &
         min_ln_factor, min_factor
      real(dp), intent(out) :: &
         factor, d_factor_dT
      real(dp) :: ln_factor, d_ln_factor_dlnT
      factor = 1d0
      d_factor_dT = 0d0
      if (temp <= T_lo) return
      if (temp >= T_hi) then
         factor = min_factor
         return
      end if
      ln_factor = min_ln_factor * (lnT - lnT_lo) / (lnT_hi - lnT_lo)
      d_ln_factor_dlnT = min_ln_factor / (lnT_hi - lnT_lo)
      factor = exp(ln_factor)
      d_factor_dT = d_ln_factor_dlnT * factor / temp
   end subroutine get_T_limit_factor
   
   
   subroutine set_molar_abundances(n, dbg, ierr)
      type (net_info) :: n
      type(net_general_info), pointer :: g
      logical, intent(in) :: dbg
      integer, intent(out) :: ierr
      
      real(dp) :: sum
      integer :: i, ci
      character (len = 256) :: message
      include 'formats'
      sum = 0
      g => n% g
      do i = 1, g% num_isos
         sum = sum + n% x(i)
         ci = g% chem_id(i)
         if (ci <= 0) then
            write(*, *) 'problem in set_molar_abundances'
            write(*, *) 'i', i
            write(*, *) 'g% num_isos', g% num_isos
            write(*, *) 'g% chem_id(i)', g% chem_id(i)
            call mesa_error(__FILE__, __LINE__, 'set_molar_abundances')
         end if
         n% y(i) = min(1d0, max(n% x(i), 0d0)) / chem_isos% Z_plus_N(ci)
      enddo
      
      return      ! let it go even with bad xsum.
      
      if (abs(sum - 1d0) > 1d-2) then
         ierr = -1
         if (dbg) then
            do i = 1, n% g% num_isos
               ci = g% chem_id(i)
               write(*, 2) chem_isos% name(ci), i, n% x(i)
            end do
            write(*, 1) 'abs(sum - 1d0)', abs(sum - 1d0)
         end if
         return
      end if
   
   end subroutine set_molar_abundances
   
   subroutine do_clean_up_fractions(nzlo, nzhi, species, nz, xa, max_sum_abs, xsum_tol, ierr)
      integer, intent(in) :: nzlo, nzhi, species, nz
      real(dp), intent(inout) :: xa(:, :) ! (species, nz)
      real(dp), intent(in) :: max_sum_abs, xsum_tol
      integer, intent(out) :: ierr
      integer :: k, op_err
      ierr = 0
      if (nzlo == nzhi) then
         call do_clean1(species, xa(1:species, nzlo), nzlo, max_sum_abs, xsum_tol, ierr)
         return
      end if
      !x$OMP  PARALLEL DO PRIVATE(k, op_err)
      do k = nzlo, nzhi
         op_err = 0
         call do_clean1(species, xa(1:species, k), k, max_sum_abs, xsum_tol, op_err)
         if (op_err /= 0) ierr = op_err
      end do
      !x$OMP  END PARALLEL DO
   end subroutine do_clean_up_fractions
   
   
   subroutine do_clean1(species, xa, k, max_sum_abs, xsum_tol, ierr)
      use utils_lib
      integer, intent(in) :: species, k
      real(dp), intent(inout) :: xa(:) ! (species)
      real(dp), intent(in) :: max_sum_abs, xsum_tol
      integer, intent(out) :: ierr
      integer :: j
      real(dp) :: xsum
      real(dp), parameter :: tiny_x = 1d-99
      character (len = 256) :: message
      if (max_sum_abs > 1) then ! check for crazy values
         xsum = sum(abs(xa(1:species)))
         if (is_bad(xsum) .or. xsum > max_sum_abs) then
            ierr = -1
            return
         end if
      end if
      ierr = 0
      do j = 1, species
         if (xa(j) < tiny_x) xa(j) = tiny_x
         if (xa(j) > 1) xa(j) = 1
      end do
      xsum = sum(xa(1:species))
      if (abs(xsum - 1) > xsum_tol) then
         ierr = -1
         return
      end if
      xa(1:species) = xa(1:species) / xsum
   end subroutine do_clean1


end module net_eval

