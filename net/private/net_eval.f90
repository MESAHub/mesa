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
      use chem_lib, only: get_mass_excess
      use net_def, only: Net_General_Info, Net_Info
      
      implicit none
      
      
      contains
      

      subroutine eval_net( &
            n, g, rates_only, just_dxdt, &
            num_isos, num_reactions, num_wk_reactions, &
            x, atemp, alogtemp, arho, alogrho, &
            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
            rate_factors, weak_rate_factor, &
            reaction_Qs, reaction_neuQs, &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
            screening_mode, &
            eps_nuc_categories, eps_neu_total, &
            lwork, work, actual_Qs, actual_neuQs, from_weaklib, symbolic, ierr)
         use net_initialize, only: &
            set_rate_ptrs, setup_net_info, set_ptrs_for_approx21
         use net_approx21, only: num_reactions_func => num_reactions
         use net_screen
         use net_derivs
         use net_def, only: &
            net_test_partials, &
            net_test_partials_val, net_test_partials_dval_dx, &
            net_test_partials_i, net_test_partials_iother

         type (Net_Info) :: n
         type (Net_General_Info), pointer :: g
         logical, intent(in) :: rates_only, just_dxdt
         integer, intent(in) :: num_isos
         integer, intent(in) :: num_reactions, num_wk_reactions
         real(dp), intent(in)  :: x(:)
         real(dp), intent(in)  :: atemp, alogtemp
         real(dp), intent(in)  :: arho, alogrho
         real(dp), intent(in)  :: abar  ! mean number of nucleons per nucleus
         real(dp), intent(in)  :: zbar  ! mean charge per nucleus
         real(dp), intent(in)  :: z2bar ! mean charge squared per nucleus
         real(dp), intent(in)  :: ye    
         real(dp), intent(in)  :: eta, d_eta_dlnT, d_eta_dlnRho ! eta and derivatives
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
         real(dp), intent(inout) :: d_dxdt_dx(:,:)
         real(dp), intent(inout) :: eps_nuc_categories(:)
         real(dp), intent(out) :: eps_neu_total
         integer, intent(in) :: screening_mode
         integer, intent(in) :: lwork
         real(dp), pointer :: work(:) ! (lwork)
         real(dp), pointer, dimension(:) :: actual_Qs, actual_neuQs ! ignore if null
         logical, pointer :: from_weaklib(:) ! ignore if null
         logical, intent(in) :: symbolic
         integer, intent(out) :: ierr

         real(dp), dimension(:), pointer :: &
            rate_raw, rate_raw_dT, rate_raw_dRho, &
            rate_screened, rate_screened_dT, rate_screened_dRho
         integer, parameter :: max_z_for_cache = 14
         real(qp), target :: dydt_a(num_rvs*num_isos)
         real(qp), pointer :: dydt(:,:) ! (num_rvs, num_isos)
         real(dp), target :: mion_a(num_isos)
         real(dp), pointer :: mion(:)
         real(dp) :: enuc, temp, logtemp, T9, rho, logrho, total, prev, curr, prev_T
         real(dp) :: btemp, bden, eps_total, Ys, sum_dxdt, compare, Z_plus_N
         real(qp) :: eps_nuc_MeV(num_rvs)
         integer :: ci, i, j, ir, weak_id, h1, iwork, approx21_num_rates
         integer, pointer :: chem_id(:)
         integer(8) :: time0, time1
         logical :: doing_timing
         
         ! for approx21
         real(dp), pointer :: dfdy(:,:)
         real(dp), dimension(:), pointer :: &
            dratdumdy1, dratdumdy2, d_epsnuc_dy, d_epsneu_dy, dydt1, dfdT, dfdRho
         real(dp) :: &
            deps_total_dRho, deps_total_dT, &
            deps_neu_dT, deps_neu_dRho, fII
               
         real(dp) :: mev2gr
         
         logical, parameter :: dbg = .false.
         !logical, parameter :: dbg = .true.
         
         include 'formats'

         if (dbg) write(*,*) 'enter eval_net'
         
         doing_timing = g% doing_timing
         if (doing_timing) then
            call system_clock(time0)
            g% doing_timing = .false.
         end if

         ierr = 0
         
         dydt(1:num_rvs,1:num_isos) => dydt_a(1:num_rvs*num_isos)
         chem_id => g% chem_id

         eps_nuc = 0

         temp = atemp; logtemp = alogtemp; rho = arho; logrho = alogrho
         call get_T_rho_args(temp, logtemp, rho, logrho, ierr)
         if (ierr /= 0) then
            if (dbg) write(*,*) 'failed in get_T_rho_args'
            return
         end if
         
         if (logtemp < rattab_tlo) then ! clip to table so can eval beta decays
            logtemp = rattab_tlo
            temp = exp10(logtemp)
         end if
         T9 = temp*1d-9
         
         n% reaction_Qs => reaction_Qs
         n% reaction_neuQs => reaction_neuQs
         n% eps_neu_total = 0
         n% weak_rate_factor = weak_rate_factor
         n% logT = logtemp
         n% temp = temp
         n% logRho = logrho
         n% rho = rho
         
         if (dbg) write(*,*) 'call set_rate_ptrs'
         call set_rate_ptrs(g, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho, lwork, work, &
            iwork, ierr) ! iwork is number of entries in work used for rates
         if (ierr /= 0) then
            if (dbg) write(*,*) 'failed in set_ptrs_in_work'
            return
         end if

         if (dbg) write(*,*) 'call setup_net_info'
         call setup_net_info( &
            g, n, eps_nuc_categories,  &
            screening_mode, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho, lwork, work, &
            iwork, ierr) ! iwork updated for amount now used in work
         if (ierr /= 0) then
            if (dbg) write(*,*) 'failed in setup_net_info'
            return
         end if

         if (g% doing_approx21) then
            approx21_num_rates = num_reactions_func(g%add_co56_to_approx21)
         else
            approx21_num_rates = -1
         end if
         
         if (g% doing_approx21) then
            call set_ptrs_for_approx21( &
               g% add_co56_to_approx21, &
               iwork, work, dfdy, dratdumdy1, dratdumdy2, &
               d_epsnuc_dy, d_epsneu_dy, dydt1, dfdT, dfdRho)
            mion => mion_a
            mev2gr = 1d6*ev2erg/(clight*clight)
            do i=1,num_isos
                mion(i) = get_mass_excess(chem_isos,g% chem_id(i))*mev2gr
            end do
         end if
         
         if (.not. g% net_has_been_defined) then
            ierr = -1
            if (dbg) write(*,*) 'failed (.not. g% net_has_been_defined)'
            return
         end if
         
         if (doing_timing) then
            call system_clock(time1)
            g% clock_net_eval = g% clock_net_eval + (time1 - time0)
            time0 = time1
         end if

         if (dbg) write(*,*) 'call set_molar_abundances'
         call set_molar_abundances(g, num_isos, x, n% y, dbg, ierr)
         if (ierr /= 0) then
            if (dbg) write(*,*) 'failed in set_molar_abundances'
            return
         end if
         
         if (num_wk_reactions > 0) then
            if (dbg) write(*,*) 'call get_weaklib_rates'
            call get_weaklib_rates(ierr)
            if (ierr /= 0) then
               if (dbg) write(*,*) 'failed in get_weaklib_rates'
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

         n% d_eps_nuc_dy(:) = 0

         ! limit range of temperatures and densities
         btemp = min(rattab_temp_hi, max(temp, rattab_temp_lo))
         bden = min(1d11, max(rho, 1d-10))
         
         if (doing_timing) then
            call system_clock(time1)
            g% clock_net_eval = g% clock_net_eval + (time1 - time0)
            time0 = time1
         end if

         if (dbg) write(*,*) 'call get_rates_with_screening'
         call get_rates_with_screening(ierr)
         if (dbg) write(*,*) 'done get_rates_with_screening'
         if (ierr /= 0) then
            if (dbg) write(*,*) 'failed in get_rates_with_screening'
            return
         end if
         
         if (rates_only) return

         d_eps_nuc_dT = 0
         d_eps_nuc_dRho = 0
         d_eps_nuc_dx(:) = 0
                  
         dxdt(:) = 0
         d_dxdt_dRho(:) = 0
         d_dxdt_dT(:) = 0
         if (.not. just_dxdt) d_dxdt_dx(:,:) = 0
         eps_nuc_categories(:) = 0
         eps_neu_total = 0
         
         if (g% doing_approx21) then
            call eval_net_approx21_procs()
            if (ierr /= 0) return                

            if (net_test_partials) then            
               net_test_partials_val = eps_nuc
               net_test_partials_dval_dx = d_eps_nuc_dx(net_test_partials_i)
               if (g% add_co56_to_approx21) then
                  write(*,*) 'net: eval_net_approx21_plus_co56'
               else
                  write(*,*) 'net: eval_net_approx21'
               end if
            end if

            return            
         end if         
     
         if (dbg) write(*,*) 'call get_derivs'
         call get_derivs(  &
             n, dydt, eps_nuc_MeV(1:num_rvs), eta, ye, &
             logtemp, btemp, bden, abar, zbar,  &
             num_reactions, rate_factors, &
             symbolic, just_dxdt, ierr)
         if (ierr /= 0) then
            if (dbg) write(*,*) 'failed in get_derivs'
            return
         end if
                  
         if (symbolic) then
            do j=1, num_isos
               do i=1, num_isos
                  d_dxdt_dx(i,j) = n% d_dydt_dy(i,j)
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
         do i=1,num_categories
            n% eps_nuc_categories(i) = Qconv*n% eps_nuc_categories(i)
         end do

         ! store the results
         do i=1,num_isos
            ci = chem_id(i)
            dxdt(i) = chem_isos% Z_plus_N(ci)*dydt(i_rate, i)
         end do
         
         if (.not. just_dxdt) call store_partials
   
         eps_nuc = eps_nuc_MeV(i_rate)*Qconv 
         d_eps_nuc_dT = eps_nuc_MeV(i_rate_dT)*Qconv 
         d_eps_nuc_dRho = eps_nuc_MeV(i_rate_dRho)*Qconv 

         eps_neu_total = n% eps_neu_total*Qconv

         if (doing_timing) then
            call system_clock(time1)
            g% clock_net_eval = g% clock_net_eval + (time1 - time0)
            g% doing_timing = .true.
         end if
         
         if (net_test_partials) then   
            !net_test_partials_val = eps_nuc
            !net_test_partials_dval_dx = d_eps_nuc_dx(net_test_partials_i)
            net_test_partials_val = &
               rate_screened(g% net_reaction(irn14_to_c12))/ &
               rate_raw(g% net_reaction(irn14_to_c12))
            net_test_partials_dval_dx = 0d0 
            write(*,*) 'net_test_partials'
         end if
         
         contains

         subroutine eval_net_approx21_procs()
            use net_approx21

            ierr = 0
            
            call approx21_special_reactions( &
               btemp, bden, abar, zbar, n% y, &
               g% use_3a_fl87, Qconv*reaction_Qs(ir_he4_he4_he4_to_c12), &
               rate_screened, rate_screened_dT, rate_screened_dRho, &
               dratdumdy1, dratdumdy2, g% add_co56_to_approx21, ierr)
            if (ierr /= 0) return            
            
            call approx21_dydt( &
               n% y, rate_screened, rate_screened, &
               dydt1, .false., g% fe56ec_fake_factor, g% min_T_for_fe56ec_fake_factor, &
               g% fe56ec_n_neut, btemp, bden, g% add_co56_to_approx21, ierr)
            if (ierr /= 0) return
               
            fII = approx21_eval_PPII_fraction(n% y, rate_screened)
            
            call get_approx21_eps_info( &
                  dydt1, rate_screened, .true., eps_total, eps_neu_total, &
                  g% add_co56_to_approx21,  ierr)
                  
            if (ierr /= 0) return               
            eps_nuc = eps_total - eps_neu_total
            
            do i=1,num_isos
               dxdt(i) = chem_isos% Z_plus_N(chem_id(i))*dydt1(i)
            end do

            if (just_dxdt) return

            call approx21_dfdy( &
               n% y, dfdy, &
               g% fe56ec_fake_factor, g% min_T_for_fe56ec_fake_factor, g% fe56ec_n_neut, &
               rate_screened, rate_screened_dT, rate_screened_dRho, &
               dratdumdy1, dratdumdy2, btemp,g% add_co56_to_approx21,  ierr)
            if (ierr /= 0) return

            call approx21_dfdT_dfdRho( & 
               
               ! NOTE: currently this gives d_eps_total_dy -- should fix to account for neutrinos too
               
               n% y, mion, dfdy, rate_screened, rate_screened_dT, rate_screened_dRho, &
               g% fe56ec_fake_factor, g% min_T_for_fe56ec_fake_factor, &
               g% fe56ec_n_neut, btemp, bden, dfdT, dfdRho, d_epsnuc_dy, g% add_co56_to_approx21,  ierr)
            if (ierr /= 0) return

            call get_approx21_eps_info( &
               dfdT, rate_screened_dT, .false., deps_total_dT, deps_neu_dT, &
               g% add_co56_to_approx21,  ierr)

            if (ierr /= 0) return
            d_eps_nuc_dT = deps_total_dT - deps_neu_dT
                              
            call get_approx21_eps_info( &
               dfdRho, rate_screened_dRho, .false., deps_total_dRho, deps_neu_dRho, &
               g% add_co56_to_approx21,  ierr)

            if (ierr /= 0) return             
            d_eps_nuc_dRho = deps_total_dRho - deps_neu_dRho
            
            call approx21_d_epsneu_dy( &
               n% y, rate_screened, &
               reaction_neuQs(irpp_to_he3), &   
               reaction_neuQs(ir34_pp2), &  
               reaction_neuQs(ir34_pp3), &  
               reaction_neuQs(irc12_to_n14), &  
               reaction_neuQs(irn14_to_c12), &  
               reaction_neuQs(iro16_to_n14), &  
               d_epsneu_dy, &
               g% add_co56_to_approx21,  ierr)
            if (ierr /= 0) return

            do i=1,num_isos
               ci = chem_id(i)
               Z_plus_N = dble(chem_isos% Z_plus_N(ci))
               d_eps_nuc_dx(i) = (d_epsnuc_dy(i) - d_epsneu_dy(i))/Z_plus_N 
               d_dxdt_dRho(i) = Z_plus_N*dfdRho(i)
               d_dxdt_dT(i) = Z_plus_N*dfdT(i)
               do j=1, num_isos
                  d_dxdt_dx(i,j) = &
                     dfdy(i,j)*Z_plus_N/chem_isos% Z_plus_N(chem_id(j))
               end do
            end do


         end subroutine eval_net_approx21_procs


         subroutine get_approx21_eps_info( &
               dydt1, rate_screened, do_eps_nuc_categories, eps_total, eps_neu_total, plus_co56, ierr)
            use net_approx21, only: approx21_eps_info
            real(dp), intent(in), dimension(:) :: dydt1, rate_screened
            logical, intent(in) :: do_eps_nuc_categories
            real(dp), intent(out) :: eps_total, eps_neu_total
            logical, intent(in) :: plus_co56
            integer, intent(out) :: ierr
            real(dp) :: Qtotal_rfe56ec, Qneu_rfe56ec
            
            ! Indexes into reaction_Qs and reaction_neuQs should be in terms of the
            ! normal rate ids not the approx21 rate ids (in net_approx21.f90)
            
            call get_Qs_rfe56ec(Qtotal_rfe56ec, Qneu_rfe56ec)

            call approx21_eps_info( &
               n, n% y, mion, dydt1, rate_screened, fII, &               
               reaction_Qs(irpp_to_he3), reaction_neuQs(irpp_to_he3), & 
               reaction_Qs(ir_he3_he3_to_h1_h1_he4), &
               reaction_Qs(ir34_pp2), reaction_neuQs(ir34_pp2), & 
               reaction_Qs(ir34_pp3), reaction_neuQs(ir34_pp3), & 
               reaction_Qs(irc12_to_n14), reaction_neuQs(irc12_to_n14), & 
               reaction_Qs(irn14_to_c12), reaction_neuQs(irn14_to_c12), & 
               reaction_Qs(iro16_to_n14), reaction_neuQs(iro16_to_n14), & 
               reaction_Qs(irn14_to_o16), &
               
               reaction_Qs(irprot_to_neut), reaction_neuQs(irprot_to_neut), & 
               reaction_Qs(irneut_to_prot), reaction_neuQs(irneut_to_prot), & 
               reaction_Qs(irni56ec_to_co56), reaction_neuQs(irni56ec_to_co56), & 
               reaction_Qs(irco56ec_to_fe56), reaction_neuQs(irco56ec_to_fe56), & 
               Qtotal_rfe56ec, Qneu_rfe56ec, &
               
               reaction_Qs(irn14ag_lite), &
               reaction_Qs(ir_he4_he4_he4_to_c12), &
               reaction_Qs(ir_c12_ag_o16), reaction_Qs(ir_o16_ag_ne20), &
               reaction_Qs(ir1212), &
               reaction_Qs(ir1216_to_mg24), reaction_Qs(ir1216_to_si28), &
               reaction_Qs(ir1616a), reaction_Qs(ir1616g), &
               reaction_Qs(ir_ne20_ag_mg24), &
               reaction_Qs(ir_mg24_ag_si28), &
               reaction_Qs(ir_si28_ag_s32), &
               reaction_Qs(ir_s32_ag_ar36), &
               reaction_Qs(ir_ar36_ag_ca40), &
               reaction_Qs(ir_ca40_ag_ti44), &
               reaction_Qs(ir_ti44_ag_cr48), &
               reaction_Qs(ir_cr48_ag_fe52), &
               reaction_Qs(ir_fe52_ag_ni56), &       
               reaction_Qs(ir_fe52_ng_fe53), &       
               reaction_Qs(ir_fe53_ng_fe54), &       
               reaction_Qs(ir_fe54_ng_fe55), &       
               reaction_Qs(ir_fe55_ng_fe56), &                              
               reaction_Qs(irfe52neut_to_fe54), &               
               reaction_Qs(irfe52aprot_to_fe54), &               
               reaction_Qs(irfe54ng_to_fe56), &               
               reaction_Qs(irfe54aprot_to_fe56), &               
               reaction_Qs(irfe52aprot_to_ni56), &               
               reaction_Qs(irfe54prot_to_ni56), &               
               reaction_Qs(irhe4_breakup), &
               reaction_Qs(irhe4_rebuild), &
               eps_total, eps_neu_total, &
               do_eps_nuc_categories, n% eps_nuc_categories, &
               .false., plus_co56, ierr)

         end subroutine get_approx21_eps_info
         
         subroutine get_Qs_rfe56ec(Qtotal, Qneu)
            use chem_def
            use rates_def, only: irfe56ec_fake_to_mn56, irfe56ec_fake_to_mn57
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
            Qtotal = reaction_Qs(ir)
            Qneu = reaction_neuQs(ir)
         end subroutine get_Qs_rfe56ec
                     
         subroutine store_partials
            integer :: i, j
            include 'formats'
            do i=1,num_isos
               ci = chem_id(i)
               Z_plus_N = dble(chem_isos% Z_plus_N(ci))
               d_eps_nuc_dx(i) = Qconv*n% d_eps_nuc_dy(i)/Z_plus_N
               dxdt(i) = Z_plus_N*dydt(i_rate, i)
               d_dxdt_dRho(i) = Z_plus_N*dydt(i_rate_dRho, i)
               d_dxdt_dT(i) = Z_plus_N*dydt(i_rate_dT, i)
               do j=1, num_isos
                  d_dxdt_dx(i,j) = &
                     n% d_dydt_dy(i,j)*Z_plus_N/chem_isos% Z_plus_N(chem_id(j))
               end do
            end do
         end subroutine store_partials
         
         subroutine get_rates_with_screening(ierr)
            use rates_def, only: reaction_inputs
            use rates_lib, only: eval_using_rate_tables
            use net_approx21, only: num_reactions_func => num_reactions
            
            integer, intent(out) :: ierr
            
            integer :: i, num
            real(dp) :: f
            logical :: okay
            
            include 'formats'

            do i=1,num_reactions
               if (g% reaction_id(i) <= 0) then
                  write(*,2) 'g% reaction_id(i)', i, g% reaction_id(i)
                  call mesa_error(__FILE__,__LINE__,'get_rates_with_screening')
               end if
            end do
            
            if (dbg) write(*,*) 'call eval_using_rate_tables'
            call eval_using_rate_tables( &
               num_reactions, g% reaction_id, g% rate_table, g% rattab_f1, nrattab,  &
               ye, logtemp, btemp, bden, rate_factors, g% logttab, &
               rate_raw, rate_raw_dT, rate_raw_dRho, ierr) 
            if (ierr /= 0) then
               if (dbg) write(*,*) 'ierr from eval_using_rate_tables'
               return
            end if
            
            if (doing_timing) then
               call system_clock(time1)
               g% clock_net_rate_tables = g% clock_net_rate_tables + (time1 - time0)
               time0 = time1
            end if

            if (g% doing_approx21) then
               call approx21_rates(g% add_co56_to_approx21,ierr)
               if (ierr /= 0) return            
            end if
            

            ! get the reaction rates including screening factors
            if (dbg) write(*,*) 'call screen_net with init=.false.'
            call screen_net( &
               g, num_isos, n% y, btemp, bden, logtemp, logrho, .false.,  &
               rate_raw, rate_raw_dT, rate_raw_dRho, &
               rate_screened, rate_screened_dT, rate_screened_dRho, &
               n% screening_mode, &
               zbar, abar, z2bar, ye, ierr)
            if (dbg) write(*,*) 'done screen_net with init=.false.'
            if (ierr /= 0) return
            if (g% doing_approx21) then
               num = num_reactions_func(g%add_co56_to_approx21)
               do i=num_reactions+1,num
                  rate_screened(i) = rate_raw(i)
                  rate_screened_dT(i) = rate_raw_dT(i)
                  rate_screened_dRho(i) = rate_raw_dRho(i)
               end do
               do i=1,num
                  dratdumdy1(i) = 0d0
                  dratdumdy2(i) = 0d0
               end do           
            end if

            
            if (doing_timing) then
               call system_clock(time1)
               g% clock_net_screen = g% clock_net_screen + (time1 - time0)
               time0 = time1
            end if
            
         end subroutine get_rates_with_screening 

         subroutine approx21_rates(plus_co56, ierr)
            use net_approx21, only: &
               approx21_pa_pg_fractions, approx21_weak_rates
            logical, intent(in) :: plus_co56
            integer, intent(out) :: ierr
            ierr = 0
            call approx21_pa_pg_fractions( &
               rate_raw, rate_raw_dT, rate_raw_dRho, ierr)
            if (ierr /= 0) return            
            call approx21_weak_rates( &
               n% y, rate_raw, rate_raw_dT, rate_raw_dRho, &
               btemp, bden, ye, eta, zbar, &
               weak_rate_factor, plus_co56, ierr)
            if (ierr /= 0) return            
         end subroutine approx21_rates


         subroutine get_weaklib_rates(ierr)
            use rates_def, only : Coulomb_Info
            use rates_lib, only: eval_weak_reaction_info, coulomb_set_context
            use net_def, only: other_kind

            type (Coulomb_Info), target :: cc_info
            type (Coulomb_Info), pointer :: cc

            integer, intent(out) :: ierr
            integer :: i, j, id, ir
            include 'formats'

            ! before getting the weaklib rates, the Coulomb_Info
            ! structure must be populated.  the ecapture routines need
            ! to know some local quantities (functions of the density,
            ! temperature, and composition), to calculate Coulomb
            ! corrections to the rates

            ierr = 0
            cc => cc_info

            call coulomb_set_context( &
               cc, temp, rho, logtemp, logrho, zbar, abar, z2bar)
            
            call eval_weak_reaction_info( &
               num_wk_reactions, &
               g% weaklib_ids(1:num_wk_reactions), &
               g% reaction_id_for_weak_reactions(1:num_wk_reactions), &
               cc, temp*1d-9, ye*rho, eta, d_eta_dlnT, d_eta_dlnRho, &
               n% lambda, n% dlambda_dlnT, n% dlambda_dlnRho, &
               n% Q, n% dQ_dlnT, n% dQ_dlnRho, &
               n% Qneu, n% dQneu_dlnT, n% dQneu_dlnRho, &
               ierr)
            if (weak_rate_factor < 1d0) then
               do i=1,num_wk_reactions
                  n% lambda(i) = weak_rate_factor*n% lambda(i)
                  n% dlambda_dlnT(i) = weak_rate_factor*n% dlambda_dlnT(i)
                  n% dlambda_dlnRho(i) = weak_rate_factor*n% dlambda_dlnRho(i)
               end do
            end if          
            if (doing_timing) then
               call system_clock(time1)
               g% clock_net_weak_rates = g% clock_net_weak_rates + (time1 - time0)
               time0 = time1
            end if
         end subroutine get_weaklib_rates
      
      end subroutine eval_net
         
         
      subroutine get_T_limit_factor( &
            temp, lnT, T_lo, T_hi, lnT_lo, lnT_hi, &
            min_ln_factor, min_factor, &
            factor, d_factor_dT)
         real(dp), intent(in) ::  &
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
         ln_factor = min_ln_factor*(lnT - lnT_lo)/(lnT_hi - lnT_lo)
         d_ln_factor_dlnT = min_ln_factor/(lnT_hi - lnT_lo)
         factor = exp(ln_factor)
         d_factor_dT = d_ln_factor_dlnT*factor/temp
      end subroutine get_T_limit_factor

         
      subroutine set_molar_abundances(g, num_isos, x, y, dbg, ierr)
         type (Net_General_Info), pointer :: g
         integer, intent(in) :: num_isos
         real(dp), intent(in) :: x(:)
         real(dp), intent(inout) :: y(:)
         logical, intent(in) :: dbg
         integer, intent(out) :: ierr
         
         real(dp) :: sum
         integer :: i, ci
         character (len=256) :: message
         include 'formats'
         sum = 0
         do i = 1, g% num_isos
            sum = sum + x(i)
            ci = g% chem_id(i)
            if (ci <= 0) then
               write(*,*) 'problem in set_molar_abundances'
               write(*,*) 'i', i
               write(*,*) 'g% num_isos', g% num_isos
               write(*,*) 'g% chem_id(i)', g% chem_id(i)
               call mesa_error(__FILE__,__LINE__,'set_molar_abundances') 
            end if
            y(i) = min(1d0, max(x(i), 0d0)) / chem_isos% Z_plus_N(ci)
         enddo
         
         
         
         return      ! let it go even with bad xsum.
         
         
         
   
         if (abs(sum - 1d0) > 1d-2) then
            ierr = -1
            if (dbg) then
               do i = 1, g% num_isos
                  ci = g% chem_id(i)
                  write(*,2) chem_isos% name(ci), i, x(i)
               end do
               write(*,1) 'abs(sum - 1d0)', abs(sum - 1d0)
            end if
            return
         end if
      
      end subroutine set_molar_abundances

      
      subroutine get_T_rho_args(temp, logtemp, rho, logrho, info)
         real(dp), intent(inout) :: temp, logtemp ! log10 of temp
         real(dp), intent(inout) :: rho, logrho ! log10 of rho
         integer, intent(out) :: info
         info = 0
         if (temp == arg_not_provided .and. logtemp == arg_not_provided) then
            info = -2
            return
         end if
         if (logtemp == arg_not_provided) logtemp = log10(temp)
         if (temp == arg_not_provided) temp = exp10(logtemp)
         if (temp <= 0) then
            info = -1
            return
         end if
         if (rho == arg_not_provided .and. logrho == arg_not_provided) then
            info = -3
            return
         end if
         if (logrho == arg_not_provided) logrho = log10(rho)
         if (rho == arg_not_provided) rho = exp10(logrho)
         if (rho <= 0) then
            info = -1
            return
         end if
      end subroutine get_T_rho_args
      
      
      subroutine do_clean_up_fractions(nzlo, nzhi, species, nz, xa, max_sum_abs, xsum_tol, ierr)
         integer, intent(in) :: nzlo, nzhi, species, nz
         real(dp), intent(inout) :: xa(:,:) ! (species, nz)
         real(dp), intent(in) :: max_sum_abs, xsum_tol
         integer, intent(out) :: ierr
         integer :: k, op_err
         ierr = 0
         if (nzlo == nzhi) then
            call do_clean1(species, xa(1: species, nzlo), nzlo, max_sum_abs, xsum_tol, ierr)
            return
         end if         
!x$OMP  PARALLEL DO PRIVATE(k, op_err)
         do k = nzlo, nzhi
            op_err = 0
            call do_clean1(species, xa(1: species, k), k, max_sum_abs, xsum_tol, op_err)
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
         character (len=256) :: message
         if (max_sum_abs > 1) then ! check for crazy values
            xsum = sum(abs(xa(1: species)))
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
         xsum = sum(xa(1: species))         
         if (abs(xsum-1) > xsum_tol) then
            ierr = -1
            return
         end if
         xa(1: species) = xa(1: species)/xsum
      end subroutine do_clean1
      

      end module net_eval

