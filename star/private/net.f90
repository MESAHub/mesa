! ***********************************************************************
!
!   Copyright (C) 2012-2019  The MESA Team
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

      module net

      use star_private_def
      use const_def
      use utils_lib, only: is_bad, mesa_error

      implicit none

      private
      public :: set_net, do_net, do1_net, do_micro_change_net, &
         get_screening_mode, default_set_which_rates, default_set_rate_factors, &
         default_set_op_mono_factors


      contains


      subroutine do_net(s, nzlo, nzhi, ierr)
         use star_utils, only: start_time, update_time
         use net_lib, only: net_work_size
         use rates_def, only: rates_other_screening
         use alloc
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         integer, intent(out) :: ierr

         logical, parameter :: use_omp = .true.
         integer :: k, op_err, net_lwork, j, jj, cnt, kmax
         integer(8) :: time0, clock_rate
         real(dp) :: total
         integer, pointer :: ks(:)
         logical, parameter :: only_dlnT = .false.
         logical :: okay, check_op_split_burn

         include 'formats'

         ierr = 0

         if (s% ctrl% eps_nuc_factor == 0d0 .and. s% ctrl% dxdt_nuc_factor == 0d0) then
            do k = nzlo, nzhi
               s% eps_nuc(k) = 0d0
               s% d_epsnuc_dlnd(k) = 0d0
               s% d_epsnuc_dlnT(k) = 0d0
               s% d_epsnuc_dx(:,k) = 0d0
               s% eps_nuc_categories(:,k) = 0d0
               s% dxdt_nuc(:,k) =  0d0
               s% d_dxdt_nuc_dRho(:,k) =  0d0
               s% d_dxdt_nuc_dT(:,k) =  0d0
               s% d_dxdt_nuc_dx(:,:,k) =  0d0
               s% eps_nuc_neu_total(k) = 0d0
               if (s% ctrl% op_split_burn) then
                  s% burn_avg_epsnuc(k) = 0d0
                  s% burn_num_iters(k) = 0
               end if
            end do
            return
         end if

         rates_other_screening => null()
         if(s% ctrl% use_other_screening) then
            rates_other_screening => s% other_screening
         end if

         net_lwork = net_work_size(s% net_handle, ierr)
         
         check_op_split_burn = s% ctrl% op_split_burn
         
         if (nzlo == nzhi) then
            call do1_net(s, nzlo, s% species, &
               s% num_reactions, &
               net_lwork, check_op_split_burn, ierr)
            return
         end if

         if (s% doing_timing) call start_time(s, time0, total)
         if (use_omp) then
            okay = .true.
!$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
            do k = nzlo, nzhi
               if (.not. okay) cycle
               op_err = 0
               call do1_net(s, k, s% species, &
                  s% num_reactions, &
                  net_lwork, check_op_split_burn, op_err)
               if (op_err /= 0) okay = .false.
            end do
!$OMP END PARALLEL DO
            if (.not. okay) ierr = -1
         else
            do k = nzlo, nzhi
               call do1_net(s, k, s% species, &
                  s% num_reactions, &
                  net_lwork, check_op_split_burn, ierr)
               if (ierr /= 0) exit
            end do
         end if
         if (s% doing_timing) call update_time(s, time0, total, s% time_nonburn_net)

      end subroutine do_net


      subroutine do1_net(s, k, species, &
            num_reactions, net_lwork, check_op_split_burn, ierr)
         use rates_def, only: std_reaction_Qs, std_reaction_neuQs, i_rate, &
            star_debugging_rates_flag, rates_test_partials_val, rates_test_partials_dval_dx
         use net_def, only: Net_Info, net_test_partials, &
            net_test_partials_val, net_test_partials_dval_dx, net_test_partials_i, &
            net_test_partials_iother
         use net_lib, only: net_get
         use star_utils, only: lookup_nameofvar
         use chem_def, only: chem_isos, category_name, i_ni56_co56, i_co56_fe56, &
            num_categories, iphoto, category_name
         use eos_def, only : i_eta
         use utils_lib,only: realloc_double, realloc_double3
         type (star_info), pointer :: s
         integer, intent(in) :: k, species, num_reactions, net_lwork
         logical, intent(in) :: check_op_split_burn
         integer, intent(out) :: ierr

         integer :: i, j, kk, screening_mode, sz, i_var, i_var_sink
         real(dp) :: log10_rho, log10_T, T, alfa, beta, eps_nuc_factor, &
            d_eps_nuc_dRho, d_eps_nuc_dT, cat_factor, tau_gamma, eps_cat_sum
         real(dp), target :: net_work_ary(net_lwork)
         real(dp), pointer :: net_work(:)
         type (Net_Info), target :: net_info_target
         type (Net_Info), pointer :: netinfo
         character (len=100) :: message
         real(dp), pointer :: reaction_neuQs(:)
         logical :: clipped_T

         logical, parameter :: dbg = .false.

         include 'formats'

         ierr = 0

         if (check_op_split_burn .and. &
             s% doing_struct_burn_mix .and. &
             s% T_start(k) >= s% ctrl% op_split_burn_min_T) then ! leave this to do_burn
            return
         end if

         net_work => net_work_ary
         netinfo => net_info_target

         s% eps_nuc(k) = 0d0
         s% d_epsnuc_dlnd(k) = 0d0
         s% d_epsnuc_dlnT(k) = 0d0
         s% d_epsnuc_dx(:,k) = 0d0
         s% eps_nuc_categories(:,k) = 0d0
         s% dxdt_nuc(:,k) =  0d0
         s% d_dxdt_nuc_dRho(:,k) =  0d0
         s% d_dxdt_nuc_dT(:,k) =  0d0
         s% d_dxdt_nuc_dx(:,:,k) =  0d0
         s% eps_nuc_neu_total(k) = 0d0

         if ((s% ctrl% eps_nuc_factor == 0d0 .and. s% ctrl% dxdt_nuc_factor == 0d0) .or. &
              s% abar(k) > s% ctrl% max_abar_for_burning) then
            return
         end if

         log10_rho = s% lnd(k)/ln10
         log10_T = s% lnT(k)/ln10
         T = s% T(k)
         
         clipped_T = (s% ctrl% max_logT_for_net > 0 .and. log10_T > s% ctrl% max_logT_for_net)
         if (clipped_T) then
            log10_T = s% ctrl% max_logT_for_net
            T = exp10(log10_T)
         end if

         screening_mode = get_screening_mode(s,ierr)
         if (ierr /= 0) then
            write(*,*) 'unknown string for screening_mode: ' // trim(s% ctrl% screening_mode)
            call mesa_error(__FILE__,__LINE__,'do1_net')
            return
         end if

         if (s% ctrl% reaction_neuQs_factor /= 1d0) then
            sz = size(std_reaction_neuQs,dim=1)
            allocate(reaction_neuQs(sz))
            do j=1,sz
               reaction_neuQs(j) = std_reaction_neuQs(j)*s% ctrl% reaction_neuQs_factor
            end do
         else
            reaction_neuQs => std_reaction_neuQs
         end if

         if (s% ctrl% solver_test_net_partials) then
             net_test_partials = (k == s% ctrl% solver_test_partials_k .and. &
               s% solver_call_number == s% ctrl% solver_test_partials_call_number .and. &
               s% solver_iter == s% ctrl% solver_test_partials_iter_number)
             ! if the test is for a partial wrt an abundance, do this
             ! in inlist set solver_test_partials_var_name and solver_test_partials_sink_name 
             ! set solver_test_partials_equ_name = ''
             i_var = lookup_nameofvar(s, s% ctrl% solver_test_partials_var_name)
             i_var_sink = lookup_nameofvar(s, s% ctrl% solver_test_partials_sink_name)
             s% solver_test_partials_var = i_var ! index in vars
             if (i_var > s% nvar_hydro) then ! index in xa for sink
                s% solver_test_partials_dx_sink = i_var_sink - s% nvar_hydro
             else
                s% solver_test_partials_dx_sink = 0
             end if
             net_test_partials_i = i_var - s% nvar_hydro ! index in xa for var
             net_test_partials_iother = i_var_sink - s% nvar_hydro ! index in xa for var
         end if
         
         if (s% ctrl% use_other_net_get) then
            call s% other_net_get( &
               s% id, k, &
               s% net_handle, .false., netinfo, species, num_reactions, s% xa(1:species,k), &
               T, log10_T, s% rho(k), log10_Rho, &
               s% abar(k), s% zbar(k), s% z2bar(k), s% ye(k), &
               s% eta(k), s% d_eos_dlnT(i_eta,k), s% d_eos_dlnd(i_eta,k), &
               s% rate_factors, s% ctrl% weak_rate_factor, &
               std_reaction_Qs, reaction_neuQs, &
               s% eps_nuc(k), d_eps_nuc_dRho, d_eps_nuc_dT, s% d_epsnuc_dx(:,k), &
               s% dxdt_nuc(:,k), s% d_dxdt_nuc_dRho(:,k), s% d_dxdt_nuc_dT(:,k), s% d_dxdt_nuc_dx(:,:,k), &
               screening_mode, s% eps_nuc_categories(:,k), &
               s% eps_nuc_neu_total(k), net_lwork, net_work, ierr)
         else
            call net_get( &
               s% net_handle, .false., netinfo, species, num_reactions, s% xa(1:species,k), &
               T, log10_T, s% rho(k), log10_Rho, &
               s% abar(k), s% zbar(k), s% z2bar(k), s% ye(k), &
               s% eta(k), s% d_eos_dlnT(i_eta,k), s% d_eos_dlnd(i_eta,k), &
               s% rate_factors, s% ctrl% weak_rate_factor, &
               std_reaction_Qs, reaction_neuQs, &
               s% eps_nuc(k), d_eps_nuc_dRho, d_eps_nuc_dT, s% d_epsnuc_dx(:,k), &
               s% dxdt_nuc(:,k), s% d_dxdt_nuc_dRho(:,k), s% d_dxdt_nuc_dT(:,k), s% d_dxdt_nuc_dx(:,:,k), &
               screening_mode, s% eps_nuc_categories(:,k), &
               s% eps_nuc_neu_total(k), net_lwork, net_work, ierr)
         end if

         if (is_bad(s% eps_nuc(k))) then
            ierr = -1
            if (s% ctrl% report_ierr) write(*,*) 'net_get returned bad eps_nuc', ierr
            if (s% ctrl% stop_for_bad_nums) call mesa_error(__FILE__,__LINE__,'do1_net')
            return
         end if

         if (-k == s% nz) then
            write(*,1) 'logT', log10_T
            call mesa_error(__FILE__,__LINE__,'net')
         end if
         !if (k == 864 .and. log10_T >= 7.522497408d0 .and. log10_T <= 7.5224974089d0) then
         if (-k == s% nz) then
            do j=1,num_categories
               write(*,2) trim(category_name(j)), j, s% eps_nuc_categories(j,k)
            end do
            write(*,'(A)')
            write(*,1) 'logRho', log10_Rho
            write(*,1) 'logT', log10_T
            write(*,1) 'eps_nuc', s% eps_nuc(k)
            write(*,1) 'sum(eps_nuc_categories)', sum(s% eps_nuc_categories(:,k))
            write(*,1) 'sum(eps_nuc_categories)/eps_nuc', &
               sum(s% eps_nuc_categories(:,k))/s% eps_nuc(k)
            write(*,2) trim(s% net_name), s% species
            call mesa_error(__FILE__,__LINE__,'after net_get in star')
         end if
         
         if (s% ctrl% solver_test_net_partials .and. net_test_partials) then
            s% solver_test_partials_val = net_test_partials_val
            s% solver_test_partials_dval_dx = net_test_partials_dval_dx
         end if
         
         if (ierr == 0) then
     
            if (clipped_T) then
               d_eps_nuc_dT = 0
               s% d_dxdt_nuc_dT(1:species,k) = 0
            end if

            if (s% ctrl% nonlocal_NiCo_kap_gamma > 0d0 .and. &
                  .not. s% ctrl% nonlocal_NiCo_decay_heat) then
               tau_gamma = 0
               do kk = 1, k
                  tau_gamma = tau_gamma + s% dm(kk)/(pi4*s% rmid(kk)*s% rmid(kk))
               end do
               tau_gamma = tau_gamma*s% ctrl% nonlocal_NiCo_kap_gamma
               s% eps_nuc(k) = s% eps_nuc(k)*(1d0 - exp(-tau_gamma))
            end if
         
            if (abs(s% eps_nuc(k)) > s% ctrl% max_abs_eps_nuc) then
               s% eps_nuc(k) = sign(s% ctrl% max_abs_eps_nuc, s% eps_nuc(k))
               d_eps_nuc_dRho = 0d0
               d_eps_nuc_dT = 0d0
               s% d_epsnuc_dx(:,k) = 0d0
            end if
            
            eps_cat_sum = sum(s% eps_nuc_categories(:,k))
            if (abs(eps_cat_sum) < 1d-10) then
               alfa = 1d0
            else
               alfa = s% eps_nuc(k)/eps_cat_sum
            end if
            if (.false. .and. abs(1d0 - alfa) > 1d-10) then
               write(*,3) 'do1_net: sum(categories) /= eps_nuc', k, s% model_number, &
                  log10(abs(alfa)), s% eps_nuc(k), eps_cat_sum, s% eps_nuc_categories(iphoto,k)
            else
               do i=1,num_categories
                  s% eps_nuc_categories(i,k) = s% eps_nuc_categories(i,k)*alfa
               end do
            end if

         end if

         if (s% ctrl% reaction_neuQs_factor /= 1d0) deallocate(reaction_neuQs)

         if (ierr /= 0) then
            if (s% ctrl% report_ierr) then
               write(*,'(A)')
               write(*,*) 'do1_net: net_get failure for cell ', k
               !return
               call show_stuff(s,k,net_lwork,net_work)
            end if
            if (is_bad_num(s% eps_nuc(k))) then
               if (s% ctrl% stop_for_bad_nums) then
                  write(*,2) 'eps_nuc', k, s% eps_nuc(k)
                  call mesa_error(__FILE__,__LINE__,'do1_net')
               end if
            end if
            return
         end if

         if (is_bad_num(s% eps_nuc(k))) then
            if (s% ctrl% stop_for_bad_nums) then
               write(*,2) 'eps_nuc', k, s% eps_nuc(k)
               call mesa_error(__FILE__,__LINE__,'do1_net')
            end if
            ierr = -1
            return
         end if

         s% d_epsnuc_dlnd(k) = d_eps_nuc_dRho*s% rho(k)
         s% d_epsnuc_dlnT(k) = d_eps_nuc_dT*s% T(k)

         eps_nuc_factor = s% ctrl% eps_nuc_factor
         if (eps_nuc_factor /= 1d0) then
            s% eps_nuc(k) = s% eps_nuc(k)*eps_nuc_factor
            s% d_epsnuc_dlnd(k) = s% d_epsnuc_dlnd(k)*eps_nuc_factor
            s% d_epsnuc_dlnT(k) = s% d_epsnuc_dlnT(k)*eps_nuc_factor
            s% d_epsnuc_dx(:,k) = s% d_epsnuc_dx(:,k)*eps_nuc_factor
            s% eps_nuc_categories(:,k) = s% eps_nuc_categories(:,k)*eps_nuc_factor
         end if

         if (s% ctrl% dxdt_nuc_factor /= 1d0) then
            s% dxdt_nuc(:,k) = s% dxdt_nuc(:,k)*s% ctrl% dxdt_nuc_factor
            s% d_dxdt_nuc_dRho(:,k) = s% d_dxdt_nuc_dRho(:,k)*s% ctrl% dxdt_nuc_factor
            s% d_dxdt_nuc_dT(:,k) = s% d_dxdt_nuc_dT(:,k)*s% ctrl% dxdt_nuc_factor
            s% d_dxdt_nuc_dx(:,:,k) = s% d_dxdt_nuc_dx(:,:,k)*s% ctrl% dxdt_nuc_factor
         end if

         if (is_bad_num(s% eps_nuc(k))) then
            write(*,*) 'k', k
            write(*,1) 's% eps_nuc(k)', s% eps_nuc(k)
            ierr = -1
            call show_stuff(s,k,net_lwork,net_work)
            write(*,*) '(is_bad_num(s% eps_nuc(k)))'
            write(*,*) 'failed in do1_net'
            if (s% ctrl% stop_for_bad_nums) call mesa_error(__FILE__,__LINE__,'do1_net')
            return
         end if
         
         if (k == -1) then
            write(*,'(A)')
            call show_stuff(s,k,net_lwork,net_work)
         end if

         if (s% model_number == -1) then
            write(*,5) 'eps_nuc', k, s% solver_iter, s% model_number, s% solver_adjust_iter, &
                        s% eps_nuc(k)
         end if
         
         if (.false.) then
            write(*,'(A)')
            call show_stuff(s,k,net_lwork,net_work)
            write(*,'(A)')
            write(*,'(A)')
            write(*,1) 's% eps_nuc(k)', s% eps_nuc(k)
            write(*,1) 's% d_epsnuc_dlnd(k)', s% d_epsnuc_dlnd(k)
            write(*,1) 's% d_epsnuc_dlnT(k)', s% d_epsnuc_dlnT(k)
            write(*,'(A)')
            write(*,*) 'do1_net'
            stop
            !ierr = -1
         end if

         if (.false.) call show_stuff(s,k,net_lwork,net_work)

      end subroutine do1_net



      subroutine show_stuff(s,k,lwork,work)
         use chem_def
         use eos_def, only : i_eta
         use rates_def
         use net_lib, only: get_reaction_id_table_ptr, get_net_rate_ptrs
         use num_lib, only: qsort
         use eos_def, only: i_eta
         type (star_info), pointer :: s
         integer, intent(in) :: k, lwork
         real(dp), pointer :: work(:)

         integer, pointer :: reaction_id(:) ! maps net reaction number to reaction id
         integer :: i, j, ierr, species, num_reactions
         real(dp) :: log10_Rho, log10_T
         real(dp), pointer :: v(:)
         integer, pointer :: index(:)
         real(dp), pointer, dimension(:) :: &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho

         include 'formats'

         logical, parameter :: do_sort = .true.

         ierr = 0
         species = s% species
         num_reactions = s% num_reactions
         log10_T = s% lnT(k)/ln10
         log10_Rho = s% lnd(k)/ln10

         call get_net_rate_ptrs(s% net_handle, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho, lwork, work, &
            ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in get_net_rate_ptrs'
            call mesa_error(__FILE__,__LINE__)
         end if

         call get_reaction_id_table_ptr(s% net_handle, reaction_id, ierr)
         if (ierr /= 0) return

         do i=1,num_reactions
            if (s% rate_factors(i) /= 1d0) then
               write(*,2) 'rate factor ' // trim(reaction_Name(reaction_id(i))), &
                  s% rate_factors(i)
            end if
         end do

         i = max(species, num_reactions)
         allocate(v(i), index(i))
         write(*,'(A)')
         if (.true.) then
            write(*, *)
            if (do_sort) then
               do j=1,num_reactions
                  v(j) = abs(rate_raw(j))
               end do
               call qsort(index, num_reactions, v)
            else
               do j=1,num_reactions
                  index(j) = j
               end do
            end if

            write(*,*) 'reaction rate_raw'
            do i=1,num_reactions
               j = index(num_reactions+1-i)
               write(*,2) trim(reaction_Name(reaction_id(j))), k, rate_raw(j)
            end do
         end if

         if (.false.) then
            write(*,'(A)')
            write(*,*) 'screened rates'
            do j=1,num_reactions
               write(*,3) 'screened rate ' // trim(reaction_Name(reaction_id(j))), &
                  j, k, rate_screened(j)
            end do
         end if

         if (.true.) then
            write(*,'(A)')
            do j=1,species
               write(*,2) 'dxdt ' // trim(chem_isos% name(s% chem_id(j))), k, s% dxdt_nuc(j, k)
            end do
            write(*,'(A)')
            do j=1,species
               write(*,2) 'd_epsnuc_dx ' // trim(chem_isos% name(s% chem_id(j))), k, s% d_epsnuc_dx(j, k)
            end do
         end if
         write(*,'(A)')

         if (.false.) then
            write(*,'(A)')
            do j=1,species
               write(*,2) 'dt*dxdt ' // trim(chem_isos% name(s% chem_id(j))), k, &
                  s% dt * s% dxdt_nuc(j, k)
            end do
         end if

         if (.true.) then
            if (do_sort) then
               do j=1,species
                  v(j) = s% xa(j,k)
               end do
               call qsort(index, species, v)
            else
               do j=1,num_reactions
                  index(j) = j
               end do
            end if
            write(*,'(A)')
            do i=1,species
               j = index(species+1-i)
               if (.true. .or. s% xa(j,k) > 1d-9) &
                  write(*,1) 'xin(net_iso(i' // &
                     trim(chem_isos% name(s% chem_id(j))) // '))= ', s% xa(j,k)
            end do
         end if
         
         if (.false.) then
            do i=1,species
               write(*,'(a,i3,a,d26.16)') 'values_for_Xinit(', i, ')= ', s% xa(i,k)
            end do
         end if

         write(*,2) 'k', k
         write(*,'(A)')
         write(*,*) 'net_name ', trim(s% net_name)
         write(*,*) 'species', species
         write(*,'(A)')
         write(*,1) 'logT =', log10_T
         write(*,1) 'T =', s% T(k)
         write(*,1) 'logRho =', log10_Rho
         write(*,1) 'rho =', s% rho(k)
         write(*,'(A)')
         write(*,1) 'eta =', s% eta(k)
         write(*,1) 'd_eta_lnT =', s% d_eos_dlnT(i_eta,k)
         write(*,1) 'd_eta_lnd =', s% d_eos_dlnd(i_eta,k)
         write(*,'(A)')
         write(*,1) 'abar =', s% abar(k)
         write(*,1) 'zbar =', s% zbar(k)
         write(*,1) 'z2bar =', s% z2bar(k)
         write(*,1) 'ye =', s% ye(k)
         write(*,*) 'screening_mode = ' // trim(s% ctrl% screening_mode)
         write(*,'(A)')

      end subroutine show_stuff



      integer function get_screening_mode(s,ierr)
         use rates_lib, only: screening_option
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         include 'formats'
         ierr = 0
         if (s% screening_mode_value >= 0) then
            get_screening_mode = s% screening_mode_value
            return
         end if
         get_screening_mode = screening_option(s% ctrl% screening_mode, ierr)
         if (ierr /= 0) return
         s% screening_mode_value = get_screening_mode
      end function get_screening_mode


      subroutine do_micro_change_net(s, new_net_name, ierr)
         use net_def
         type (star_info), pointer :: s
         character (len=*), intent(in) :: new_net_name
         integer, intent(out) :: ierr
         ierr = 0
         s% net_name = new_net_name
         call set_net(s, new_net_name, ierr)
      end subroutine do_micro_change_net


      subroutine set_net(s, new_net_name, ierr)
         use net_lib
         use utils_lib, only: realloc_double
         use alloc, only: update_nvar_allocs, set_chem_names
         use chem_def, only: ih1, ihe4
         use rates_def
         type (star_info), pointer :: s
         character (len=*), intent(in) :: new_net_name
         integer, intent(out) :: ierr

         integer :: i, ir
         integer :: old_num_reactions, old_nvar_chem, old_species
         integer, parameter :: num_lowT_rates = 10
         integer, pointer :: net_reaction_ptr(:)

         include 'formats'

         old_num_reactions = s% num_reactions

         if (s% net_handle /= 0) call free_net_handle(s% net_handle)

         s% net_handle = alloc_net_handle(ierr)
         if (ierr /= 0) then
            write(*,*) 'set_net failed in alloc_net_handle'
            return
         end if
         call net_ptr(s% net_handle, s% net_rq, ierr)
         if (ierr /= 0) then
            write(*,*) 'set_net failed in net_ptr'
            return
         end if

         call net_tables(s, ierr)
         if (ierr /= 0) then
            write(*,*) 'set_net failed in net_tables'
            return
         end if

         old_species = s% species
         s% species = net_num_isos(s% net_handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'set_net failed in net_num_isos'
            return
         end if

         s% num_reactions = net_num_reactions(s% net_handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'set_net failed in net_num_reactions'
            return
         end if

         old_nvar_chem = s% nvar_chem
         s% nvar_chem = s% species
         call update_nvar_allocs(s, s% nvar_hydro, old_nvar_chem, ierr)
         if (ierr /= 0) then
            write(*,*) 'set_net failed in update_nvar_allocs'
            return
         end if

         call get_chem_id_table_ptr(s% net_handle, s% chem_id, ierr)
         if (ierr /= 0) then
            write(*,*) 'set_net failed in get_chem_id_table_ptr'
            return
         end if

         call get_net_iso_table_ptr(s% net_handle, s% net_iso, ierr)
         if (ierr /= 0) then
            write(*,*) 'set_net failed in get_net_iso_table_ptr'
            return
         end if
         
         if (s% net_iso(ih1) == 0 .or. s% net_iso(ihe4) == 0) then
            write(*,*) 'mesa/star requires both h1 and he4 in net isotopes'
            write(*,*) 'but they are not included in ' // trim(new_net_name)
            ierr = -1
            return
         end if

         if (associated(s% xa_removed)) deallocate(s% xa_removed)
         allocate(s% xa_removed(s% species))

         call set_chem_names(s)

         call s% set_rate_factors(s% id, ierr)
         if (ierr /= 0) then
            write(*,*) 'set_net failed in s% set_rate_factors'
            return
         end if

         if (associated(s% op_mono_factors)) deallocate(s% op_mono_factors)
         allocate(s% op_mono_factors(s% species))

         call s% set_op_mono_factors(s% id, ierr)
         if (ierr /= 0) then
            write(*,*) 'set_net failed in s% set_op_mono_factors'
            return
         end if
         
         s% need_to_setvars = .true.

      end subroutine set_net


      subroutine net_tables(s, ierr)
         use net_lib ! setup net
         use rates_lib
         use rates_def, only: rates_reaction_id_max
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0

         call net_start_def(s% net_handle, ierr)
         if (ierr /= 0) then
            if (s% ctrl% report_ierr) write(*,*) 'failed in net_start_def'
            return
         end if

         if (len_trim(s% net_name) == 0) then
            write(*,*) 'missing net_name -- please set it and try again'
            ierr = -1
            return
         end if

         call read_net_file(s% net_name, s% net_handle, ierr)
         if (ierr /= 0) then
            if (s% ctrl% report_ierr) write(*,*) 'failed in read_net_file ' // trim(s% net_name)
            return
         end if

         call net_finish_def(s% net_handle, ierr)
         if (ierr /= 0) then
            if (s% ctrl% report_ierr) write(*,*) 'failed in net_finish_def'
            return
         end if

         if (associated(s% rate_factors)) deallocate(s% rate_factors)
         allocate(s% rate_factors(rates_reaction_id_max))

         call s% set_rate_factors(s% id, ierr)
         if (ierr /= 0) then
            if (s% ctrl% report_ierr) write(*,*) 'failed in set_rate_factors'
            return
         end if

         if (associated(s% which_rates)) deallocate(s% which_rates)
         allocate(s% which_rates(rates_reaction_id_max))

         call s% set_which_rates(s% id, ierr)
         if (ierr /= 0) then
            if (s% ctrl% report_ierr) write(*,*) 'failed in set_which_rates'
            return
         end if

         call net_set_which_rates(s% net_handle, s% which_rates, ierr)
         if (ierr /= 0) then
            if (s% ctrl% report_ierr) write(*,*) 'failed in net_set_which_rates'
            return
         end if

         call net_set_logTcut(s% net_handle, s% ctrl% net_logTcut_lo, s% ctrl% net_logTcut_lim, ierr)
         if (ierr /= 0) then
            if (s% ctrl% report_ierr) write(*,*) 'failed in net_set_logTcut'
            return
         end if

         call net_set_fe56ec_fake_factor( &
            s% net_handle, s% ctrl% fe56ec_fake_factor, s% ctrl% min_T_for_fe56ec_fake_factor, ierr)
         if (ierr /= 0) then
            if (s% ctrl% report_ierr) write(*,*) 'failed in net_set_fe56ec_fake_factor'
            return
         end if

         call net_setup_tables( &
            s% net_handle, rates_cache_suffix_for_star, ierr)
         if (ierr /= 0) then
            if (s% ctrl% report_ierr) write(*,*) 'failed in net_setup_tables'
            return
         end if

      end subroutine net_tables


      subroutine default_set_which_rates(id, ierr)
         use rates_def, only: rates_NACRE_if_available
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% which_rates(:) = rates_NACRE_if_available
      end subroutine default_set_which_rates


      subroutine default_set_rate_factors(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% rate_factors(:) = 1
      end subroutine default_set_rate_factors


      subroutine default_set_op_mono_factors(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% op_mono_factors(:) = 1
      end subroutine default_set_op_mono_factors


      end module net

