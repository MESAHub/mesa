! ***********************************************************************
!
!   Copyright (C) 2010  The MESA Team
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

      module net_def
      
      use const_def, only: dp, qp
      
      implicit none


      ! reaction_kind array in Net_General_Info
         integer, parameter :: neut_kind = 1 ! involves neut but no prot
         integer, parameter :: prot_kind = neut_kind + 1 ! involves prot and perhaps neut
         integer, parameter :: other_strong_kind = prot_kind + 1 ! strong, no neut or prot
         integer, parameter :: weak_kind = other_strong_kind + 1
         integer, parameter :: num_kinds = weak_kind

      ! for reaction_reaclib_kind array in Net_General_Info
         integer, parameter :: other_kind = 0 
            ! includes weak reactions and reactions that don't have reverse in net
            ! and one of each pair of 2 to 2 reactions (including np, pa, na)
         integer, parameter :: ng_kind = other_kind + 1
         integer, parameter :: pn_kind = ng_kind + 1
         integer, parameter :: pg_kind = pn_kind + 1
         integer, parameter :: ap_kind = pg_kind + 1
         integer, parameter :: an_kind = ap_kind + 1
         integer, parameter :: ag_kind = an_kind + 1
         integer, parameter :: general_one_one_kind = ag_kind + 1 ! 1 species in and 1 out (e.g., 3alfa)
         integer, parameter :: general_two_one_kind = general_one_one_kind + 1 ! 2 species in and 1 out
         integer, parameter :: general_two_two_kind = general_two_one_kind + 1 ! 2 species in and 2 out
         integer, parameter :: max_kind = general_two_two_kind
         





         
      type Net_General_Info ! things that are constant for the particular net
      ! it is okay to have multiple threads using the same instance of this simultaneously.

         integer :: num_isos ! total number in current net            
         integer :: num_reactions ! total number of reactions for current net
         
         logical :: doing_approx21, add_co56_to_approx21
         
         integer :: approx21_ye_iso ! e.g., icr56 for fake fe56ec
         integer :: fe56ec_n_neut ! number of neutrons consumed per fake fe56ec
         
         character (len=32) :: cache_suffix

         ! isotopes
         integer, pointer :: net_iso(:) ! maps chem id to net iso number
         ! index from 1 to num_chem_isos
         ! value is 0 if the iso is not in the current net
         ! else is value between 1 and num_isos in current net
         integer, pointer :: chem_id(:) ! maps net iso number to chem id
         ! index from 1 to num_isos in current net
         ! value is between 1 and num_chem_isos         

         ! reactions
                  
         integer, pointer :: net_reaction(:) ! maps reaction id to net reaction number
         ! index from 1 to rates_reaction_id_max (in rates_def)   
         ! value is 0 if the reaction is not in the current net
         ! else is value between 1 and num_reactions in current net
         integer, allocatable :: reaction_id(:) ! maps net reaction number to reaction id
         ! index from 1 to num_reactions in current net
         ! value is between 1 and rates_reaction_id_max (in rates_def)     

         integer, allocatable :: reaction_kind(:)

         integer, pointer :: reaction_reaclib_kind(:)
         integer, pointer :: reverse_id_for_kind_ne_other(:)
         
         integer, allocatable :: reaction_max_Z(:)
         integer, allocatable:: reaction_max_Z_plus_N_for_max_Z(:)
         
         ! extra info
         
         ! strong rates cutoff smoothly for logT < logTcut_lim
         real(dp) :: logTcut_lim
         ! strong rates are zero logT < logTcut_lo
         real(dp) :: logTcut_lo
         
         ! equilibrium eps_nuc cancelation for ng, pg, pn reactions
         ! at high T, these reactions are assumed in equilibrium with their reverses,
         ! so no net eps_nuc from the pair
         real(dp) :: logT_lo_eps_nuc_cancel ! no cancelation for logT <= this
         real(dp) :: logT_hi_eps_nuc_cancel ! full cancelation for logT >= this
         
         real(dp) :: fe56ec_fake_factor, min_T_for_fe56ec_fake_factor
   
         ! the following is private info for the implementation
         
         ! tables for screen5
         real(dp), allocatable :: zs13(:) ! (num_reactions) ! zs13 = (z1+z2)**(1./3.)
         real(dp), allocatable :: zhat(:) ! (num_reactions)
         real(dp), allocatable :: zhat2(:) ! (num_reactions)
         real(dp), allocatable :: lzav(:) ! (num_reactions)
         real(dp), allocatable :: aznut(:) ! (num_reactions)
         real(dp), allocatable :: zs13inv(:) ! (num_reactions) ! zs13inv = 1 / zs13
   
         ! info for evaluation of the raw reaction rates
         real(dp), pointer :: rate_table(:,:) ! (nrate_table,num_reactions)
         real(dp), pointer :: rattab_f1(:) ! =(4,nrattab,num_reactions) ! for interpolation
         real(dp), allocatable  :: ttab(:) ! (nrate_table)
         real(dp), allocatable  :: logttab(:) ! (nrate_table)

         ! Precomputed powers of Z
         real(dp), allocatable :: & ! (num_isos)
                           z158(:), & ! screen z**1.58
                           z52(:) ! columb z**5/2

         ! info for evaluation of weak rates
         integer :: num_wk_reactions ! number of weak reactions in the current net
         integer, pointer :: &
            weaklib_ids(:), & ! (1:num_wk_reactions) = num in 1:num_weak_reactions from rates_def
               ! get_weak_rate_id from rates_lib
               ! set for rates in weak_info_list file and weakreactions.tables
            weak_reaction_index(:), & ! (1:num_reactions) = num in 1:num_wk_reactions
            weak_reaction_num(:), & ! (1:num_wk_reactions) = num in 1:num_reactions
            reaction_id_for_weak_reactions(:) ! (1:num_wk_reactions) = rates reaction id
         
         ! top level file name for net
         character (len=256) :: net_filename
         
         ! timing
         logical :: doing_timing
         ! the following are sums of results from system_clock. 
         ! divide by clock_rate to get seconds.
         ! must set all of these to 0 before change doing_timing to true.
         integer(8) :: clock_net_eval
         integer(8) :: clock_net_weak_rates
         integer(8) :: clock_net_rate_tables
         integer(8) :: clock_net_screen
         integer(8) :: clock_net_derivs
         integer(8) :: clock_derivs_select
         integer(8) :: clock_derivs_setup
         integer(8) :: clock_derivs_general
         integer(8) :: clock_net_get
                  
         ! bookkeeping
         integer :: handle
         logical :: net_has_been_defined
         logical :: in_use

      end type Net_General_Info

      integer, parameter :: num_weak_info_arrays_in_Net_Info = 9 ! weaklib results
      
            
      type Net_Info
         ! this is working storage for the nuclear reaction calculations
         
         ! pointers to caller supplied arrays ----------------------------------

         real(dp), pointer :: reaction_Qs(:) ! if null, use standard values         
         real(dp), pointer :: reaction_neuQs(:) ! if null, use standard values

         real(dp), pointer :: eps_nuc_categories(:) ! (num_categories)
         ! eps_nuc subtotals for each reaction category
         
         real(dp), pointer, dimension(:) :: &
            rate_screened, rate_screened_dT, rate_screened_dRho ! (num_rates)
         ! the units here depend on the number of reactants.
         ! in all cases, the rate_screened times as many molar fractions as there are reactants
            ! gives a number with the same units as dy/dt.
         ! so for a 2-body reaction, there are 2 Y factors, each with units [moles/gram]
            ! and the rate_screened units for such a reaction are [grams/(mole-sec)], 
            ! which when multiplied by [moles/gram]^2 gives the same units as dydt.
         ! for a 1-body reaction (e.g., a decay),
         ! there is only 1 Y factor, so the units are [1/second].
         ! similarly, a 3 body reaction will have rate_screened
         ! with units of [gram^2/(mole^2-sec)].

         real(dp), pointer, dimension(:) :: &
            rate_raw, rate_raw_dT, rate_raw_dRho ! (num_rates)
         ! raw rates are unscreened (but include density factors)
                  
         ! pointers into work array ----------------------------------

         ! molar fractions and their rates of change
         real(dp), pointer :: y(:) ! units [moles/gram]     (num_isos)
         real(dp), pointer :: d_dydt_dy(:,:) ! units [1/second] (num_isos, num_isos)
         real(dp), pointer :: d_eps_nuc_dy(:) ! (num_isos)
         
         ! weaklib results
         real(dp), dimension(:), pointer :: &
            lambda, dlambda_dlnT, dlambda_dlnRho, &
            Q, dQ_dlnT, dQ_dlnRho, &
            Qneu, dQneu_dlnT, dQneu_dlnRho

         type (Net_General_Info), pointer  :: g

         integer :: screening_mode

         real(dp) :: temp, logT, rho, logRho

         real(dp) :: eps_neu_total
         real(dp) :: weak_rate_factor
      
      end type Net_Info


      ! Interface for net hooks
      interface
         subroutine other_net_derivs_interface( &
            n, dydt, eps_nuc_MeV, eta, ye, logtemp, temp, den, abar, zbar, &
            num_reactions, rate_factors, &
            symbolic, just_dydt, ierr)
         import dp, qp, Net_Info
         implicit none

         type(Net_Info), pointer :: n
         real(qp), pointer, intent(inout) :: dydt(:,:)
         real(qp), intent(out) :: eps_nuc_MeV(:)
         integer, intent(in) :: num_reactions
         real(dp), intent(in) ::eta, ye, logtemp, temp, den, abar, zbar, &
            rate_factors(:)
         logical, intent(in) :: symbolic, just_dydt
         integer, intent(out) :: ierr

         end subroutine other_net_derivs_interface

      end interface

      ! Other net_derivs handling
      procedure(other_net_derivs_interface), pointer  :: &
         net_other_net_derivs => null()

      
   ! private to the implementation
      integer, parameter :: max_net_handles = 10
      type (Net_General_Info), target :: net_handles(max_net_handles)
      
      character (len=256) :: net_dir

      integer :: weak_rate_id_for_ni56_ec, weak_rate_id_for_co56_ec


      ! parameters for net burn

      integer, parameter :: i_burn_caller_id = 1
      integer, parameter :: i_net_handle = 2
      integer, parameter :: i_screening_mode = 3
      integer, parameter :: i_net_lwork = 4
      integer, parameter :: i_eos_handle = 5
      integer, parameter :: i_sparse_format = 6
      integer, parameter :: i_clip = 7
      integer, parameter :: i_ntimes = 8
      
      integer, parameter :: burn_lipar = i_ntimes

      integer, parameter :: r_burn_temp = 1
      integer, parameter :: r_burn_lgT = 2
      integer, parameter :: r_burn_rho = 3
      integer, parameter :: r_burn_lgRho = 4
      integer, parameter :: r_burn_eta = 5
      integer, parameter :: r_burn_theta = 6
      integer, parameter :: r_burn_time_net = 7
      integer, parameter :: r_burn_prev_lgT = 8
      integer, parameter :: r_burn_prev_lgRho = 9
      integer, parameter :: r_burn_prev_eta = 10

      integer, parameter :: burn_lrpar = r_burn_prev_eta

      integer, parameter :: r_burn_const_P_rho = 1
      integer, parameter :: r_burn_const_P_pressure = 2
      integer, parameter :: r_burn_const_P_init_rho = 3
      integer, parameter :: r_burn_const_P_time_net = 4
      integer, parameter :: r_burn_const_P_time_eos = 5
      integer, parameter :: r_burn_const_P_temperature = 6
      integer, parameter :: r_burn_const_P_init_lnS = 7
      integer, parameter :: r_burn_const_P_lnS = 8
      
      integer, parameter :: burn_const_P_lrpar = r_burn_const_P_lnS

      logical :: net_test_partials
      real(dp) :: net_test_partials_val, net_test_partials_dval_dx
      integer :: net_test_partials_i, net_test_partials_iother


      contains

      
      subroutine do_net_def_init
         use const_def, only: mesa_data_dir
         use rates_lib, only: get_weak_rate_id
         integer :: i

         net_test_partials = .false.
         net_dir = trim(mesa_data_dir) // '/net_data'
         do i=1, max_net_handles
            net_handles(i)% handle = i
            net_handles(i)% in_use = .false.
            net_handles(i)% net_has_been_defined = .false.
            net_handles(i)% num_isos = 0
            net_handles(i)% num_reactions = 0
         end do
              
         weak_rate_id_for_ni56_ec = get_id('ni56','co56')
         weak_rate_id_for_co56_ec = get_id('co56','fe56')
         
         contains
         
         integer function get_id(iso1, iso2)
            character(len=*), intent(in) :: iso1, iso2
            include 'formats'
            get_id = get_weak_rate_id(iso1, iso2)
            if (get_id == 0) then
               write(*,2) 'failed to find weak reaction for ' // trim(iso1) &
                  // ' to ' // trim(iso2) 
            end if
         end function get_id
         
      end subroutine do_net_def_init


      integer function do_alloc_net(ierr)
         integer, intent(out) :: ierr
         integer :: i
         ierr = 0
         do_alloc_net = -1
!$omp critical (net_handle)
         do i = 1, max_net_handles
            if (.not. net_handles(i)% in_use) then
               net_handles(i)% in_use = .true.
               do_alloc_net = i
               exit
            end if
         end do
!$omp end critical (net_handle)
         if (do_alloc_net == -1) then
            ierr = -1
            return
         end if
         if (net_handles(do_alloc_net)% handle /= do_alloc_net) then
            ierr = -1
            return
         end if
         call init_net_handle_data(do_alloc_net)
      end function do_alloc_net
      
      
      subroutine init_net_handle_data(handle)
         use rates_def
         integer, intent(in) :: handle
         type (Net_General_Info), pointer :: g
         g => net_handles(handle)
         call do_free_net(handle)
         g% in_use = .true.
         g% doing_approx21 = .false.
         g% add_co56_to_approx21 = .false.
         g% approx21_ye_iso = -1
         g% doing_timing = .false.
         g% logTcut_lo = rattab_tlo
         g% logTcut_lim = rattab_tlo + 0.1d0
         g% logT_lo_eps_nuc_cancel = 9.4d0
         g% logT_hi_eps_nuc_cancel = 9.5d0
         g% fe56ec_fake_factor = 1d-4
         g% min_T_for_fe56ec_fake_factor = 3d9
         g% cache_suffix = '0'
      end subroutine init_net_handle_data

      subroutine do_free_net(handle)
         use rates_def
         integer, intent(in) :: handle
         type (Net_General_Info), pointer :: g
         if (handle >= 1 .and. handle <= max_net_handles) then
            g => net_handles(handle)
            if (associated(g% net_iso)) then
               deallocate(g% net_iso)
                  nullify(g% net_iso)
            end if
            if (associated(g% chem_id)) then
               deallocate(g% chem_id)
                  nullify(g% chem_id)
            end if
            if (associated(g% net_reaction)) then
               deallocate(g% net_reaction)
                  nullify(g% net_reaction)
            end if
            if (allocated(g% reaction_id)) then
               deallocate(g% reaction_id)
            end if
            if (allocated(g% reaction_kind)) then
               deallocate(g% reaction_kind)
            end if

            if (associated(g% reaction_reaclib_kind)) then
               deallocate(g% reaction_reaclib_kind)
                  nullify(g% reaction_reaclib_kind)
            end if
            if (associated(g% reaction_id_for_weak_reactions)) then
               deallocate(g% reaction_id_for_weak_reactions)
                  nullify(g% reaction_id_for_weak_reactions)
            end if
            if (associated(g% reverse_id_for_kind_ne_other)) then
               deallocate(g% reverse_id_for_kind_ne_other)
                  nullify(g% reverse_id_for_kind_ne_other)
            end if

            if (allocated(g% reaction_max_Z)) then
               deallocate(g% reaction_max_Z)
            end if
            if (allocated(g% reaction_max_Z_plus_N_for_max_Z)) then
               deallocate(g% reaction_max_Z_plus_N_for_max_Z)
            end if
            if (allocated(g% zs13)) then
               deallocate(g% zs13)
            end if
            if (allocated(g% zhat)) then
               deallocate(g% zhat)
            end if
            if (allocated(g% zhat2)) then
               deallocate(g% zhat2)
            end if
            if (allocated(g% lzav)) then
               deallocate(g% lzav)
            end if
            if (allocated(g% aznut)) then
               deallocate(g% aznut)
            end if
            if (allocated(g% zs13inv)) then
               deallocate(g% zs13inv)
            end if
            if (allocated(g% z158)) then
               deallocate(g% z158)
            end if
            if (allocated(g% z52)) then
               deallocate(g% z52)
            end if
            if (associated(g% rate_table)) then
               deallocate(g% rate_table)
                  nullify(g% rate_table)
            end if
            if (allocated(g% ttab)) then
               deallocate(g% ttab)
            end if
            if (allocated(g% logttab)) then
               deallocate(g% logttab)
            end if
            if (associated(g% rattab_f1)) then
               deallocate(g% rattab_f1)
                  nullify(g% rattab_f1)
            end if
            if (associated(g% weaklib_ids)) then
               deallocate(g% weaklib_ids)
                  nullify(g% weaklib_ids)
            end if
            if (associated(g% weak_reaction_num)) then
               deallocate(g% weak_reaction_num)
                  nullify(g% weak_reaction_num)
            end if
            if (associated(g% weak_reaction_index)) then
               deallocate(g% weak_reaction_index)
                  nullify(g% weak_reaction_index)
            end if
            g% in_use = .false.
            g% net_has_been_defined = .false.
            g% num_isos = 0
            g% num_reactions = 0
            g% num_wk_reactions = 0
         end if
         
         
      end subroutine do_free_net
      

      subroutine get_net_ptr(handle, g, ierr)
         integer, intent(in) :: handle
         type (Net_General_Info), pointer :: g
         integer, intent(out):: ierr         
         if (handle < 1 .or. handle > max_net_handles) then
            ierr = -1
            return
         end if
         g => net_handles(handle)
         ierr = 0
      end subroutine get_net_ptr


      integer function get_net_timing_total(g)
         type (Net_General_Info), pointer :: g
         get_net_timing_total = 0
         if (.not. g% doing_timing) return
         get_net_timing_total = &
            g% clock_net_eval + &
            g% clock_net_weak_rates + &
            g% clock_net_rate_tables + &
            g% clock_net_screen + &
            g% clock_net_derivs
      end function get_net_timing_total


      subroutine zero_net_timing(g)
         type (Net_General_Info), pointer :: g
         g% clock_net_eval = 0
         g% clock_net_weak_rates = 0
         g% clock_net_rate_tables = 0 
         g% clock_net_screen = 0 
         g% clock_net_derivs = 0
         
         g% clock_derivs_setup = 0
         g% clock_derivs_select = 0
         g% clock_derivs_general = 0
         g% clock_net_get = 0
      end subroutine zero_net_timing
      
      subroutine do_net_set_fe56ec_fake_factor( &
            handle, fe56ec_fake_factor, min_T_for_fe56ec_fake_factor, ierr)
         integer, intent(in) :: handle
         real(dp), intent(in) :: fe56ec_fake_factor, min_T_for_fe56ec_fake_factor
         integer, intent(out) :: ierr
         type (Net_General_Info), pointer :: g
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for do_net_set_fe56ec_fake_factor'
            return
         end if
         g% fe56ec_fake_factor = fe56ec_fake_factor
         g% min_T_for_fe56ec_fake_factor = min_T_for_fe56ec_fake_factor
      end subroutine do_net_set_fe56ec_fake_factor
      
      
      subroutine do_net_set_logTcut(handle, logTcut_lo, logTcut_lim, ierr)
         integer, intent(in) :: handle
         real(dp), intent(in) :: logTcut_lo 
         real(dp), intent(in) :: logTcut_lim 
         integer, intent(out) :: ierr
         type (Net_General_Info), pointer :: g
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for net_set_logTcut'
            return
         end if
         g% logTcut_lo = logTcut_lo
         g% logTcut_lim = logTcut_lim
      end subroutine do_net_set_logTcut
      
      
      subroutine do_net_set_eps_nuc_cancel( &
            handle, logT_lo_eps_nuc_cancel, logT_hi_eps_nuc_cancel, ierr)
         integer, intent(in) :: handle
         real(dp), intent(in) :: logT_lo_eps_nuc_cancel 
         real(dp), intent(in) :: logT_hi_eps_nuc_cancel 
         integer, intent(out) :: ierr
         type (Net_General_Info), pointer :: g
         call get_net_ptr(handle, g, ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for net_set_eps_nuc_cancel'
            return
         end if
         g% logT_lo_eps_nuc_cancel = logT_lo_eps_nuc_cancel
         g% logT_hi_eps_nuc_cancel = logT_hi_eps_nuc_cancel
      end subroutine do_net_set_eps_nuc_cancel



      end module net_def

