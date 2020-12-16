! ***********************************************************************
!
!  Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!  MESA is free software; you can use it and/or modify
!  it under the combined terms and restrictions of the MESA MANIFESTO
!  and the GNU General Library Public License as published
!  by the Free Software Foundation; either version 2 of the License,
!  or (at your option) any later version.
!
!  You should have received a copy of the MESA MANIFESTO along with
!  this software; if not, it is available at the mesa website:
!  http://mesa.sourceforge.net/
!
!  MESA is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!  See the GNU Library General Public License for more details.
!
!  You should have received a copy of the GNU Library General Public License
!  along with this software; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module rates_def
      use utils_def, only: integer_dict
      use const_def, only: dp
      use chem_def, only: iso_name_length, nuclide_data, npart
      
      implicit none

      ! weaklib
      
      character (len=256) :: weak_data_dir         

      ! ecapture

      character (len=1000) :: ecapture_states_file
      character (len=1000) :: ecapture_transitions_file
      
      ! reaclib


      integer, parameter :: max_nreaclib=100000
      integer, parameter :: max_species_per_reaction=6
      integer, parameter :: ncoefficients=7
      integer, parameter :: nchapters=11
      integer, parameter :: ninverse_coeff = 2
      integer, parameter :: max_terms_per_rate = 20
      integer, parameter :: max_id_length = 36


      ! i/o parameters
      integer, parameter :: internal_format = 0
      integer, parameter :: pretty_print_format = 1
      integer, parameter :: short_format = 2


      ! flags for reaction types -- values are chapter numbers
      integer, parameter :: &
         r_one_one = 1, &
         r_one_two = 2, &
         r_one_three = 3, &
         r_two_one = 4, &
         r_two_two = 5, &
         r_two_three = 6, &
         r_two_four = 7, &
         r_three_one = 8, &
         r_three_two = 9, &
         r_four_two = 10, &
         r_one_four = 11


      integer, dimension(nchapters) :: Nin = (/1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 1/)
      integer, dimension(nchapters) :: Nout = (/1, 2, 3, 1, 2, 3, 4, 1, 2, 2, 4/)


      type reaclib_data
         integer, dimension(:), pointer :: chapter=>NULL()
         character(len=iso_name_length), dimension(:,:), pointer :: species=>NULL()
         character(len=iso_name_length), dimension(:), pointer :: label=>NULL()
         character, dimension(:), pointer :: reaction_flag=>NULL()
         character, dimension(:), pointer :: reverse_flag=>NULL()
         real(dp), dimension(:), pointer :: Qvalue=>NULL()
         real(dp), dimension(:,:), pointer :: coefficients=>NULL()
         ! following are for 1D allocation of 2D arrays
         character(len=iso_name_length), dimension(:), pointer :: species1=>NULL()
         real(dp), dimension(:), pointer :: coefficients1=>NULL()
      end type reaclib_data


      type reaction_data
      
         integer :: nreactions
         integer :: nchapters_included
         integer, dimension(nchapters) :: chapters_present
         type(nuclide_data), pointer :: nuclides=>NULL()
         
         integer :: num_from_weaklib
         integer, dimension(:), pointer :: weaklib_ids=>NULL() ! (num_from_weaklib)
         logical, dimension(:), pointer :: also_in_reaclib=>NULL() ! (num_from_weaklib)
         
         integer, dimension(2, nchapters) :: bookmarks
         character(len=max_id_length), dimension(:), pointer :: reaction_handle=>NULL(), reverse_handle=>NULL()
         integer, dimension(:), pointer :: category=>NULL()
         integer, dimension(:), pointer :: chapter=>NULL()
         integer, dimension(:, :), pointer :: pspecies=>NULL()
         character, dimension(:), pointer :: reaction_flag=>NULL()
         real(dp), dimension(:), pointer :: weight=>NULL()
         real(dp), dimension(:), pointer :: weight_reverse=>NULL()
         real(dp), dimension(:, :), pointer :: coefficients=>NULL()
         real(dp), dimension(:), pointer :: weak_mask=>NULL()
         real(dp), dimension(:, :), pointer :: inverse_coefficients=>NULL()
         integer, dimension(:), pointer :: inverse_exp=>NULL()
         real(dp), dimension(:, :), pointer :: inverse_part=>NULL()
         real(dp), dimension(:), pointer :: Q=>NULL()
         real(dp), dimension(:), pointer :: Qneu=>NULL()
         type (integer_dict), pointer :: reaction_dict=>NULL(), reverse_dict=>NULL()
         
      end type reaction_data
      
      
      character (len=256) :: reaclib_dir, reaclib_filename


         ! reactions information
         
         integer, parameter :: maxlen_reaction_Name = 32
         character (len=maxlen_reaction_Name), pointer :: reaction_Name(:)=>NULL() ! (rates_reaction_id_max)
            
         integer, parameter :: maxlen_reaction_Info = 72
         character (len=maxlen_reaction_Info), pointer :: reaction_Info(:)=>NULL() ! (rates_reaction_id_max)
         
         real(dp), pointer :: std_reaction_Qs(:)=>NULL() ! (rates_reaction_id_max) 
            ! set at initialization; read-only afterwards.
            ! avg energy including neutrinos
         
         real(dp), pointer :: std_reaction_neuQs(:)=>NULL() ! (rates_reaction_id_max) 
            ! set at initialization; read-only afterwards.
            ! avg neutrino loss
         
         real(dp), pointer :: weak_lowT_rate(:)=>NULL() ! (rates_reaction_id_max) 
            ! these are from reaclib or weak_info.list
            ! set at initialization; read-only afterwards.
            
         integer, pointer :: reaction_screening_info(:,:)=>NULL() !(3,rates_reaction_id_max)
            ! reaction_screen_info(1:2,i) = [chem_id1, chem_id2] for screening.  0's if no screening.
            
         integer, pointer :: weak_reaction_info(:,:)=>NULL() ! (2,rates_reaction_id_max)
            ! weak_reaction_info(1:2,i) = [chem_id_in, chem_id_out].  0's if not a weak reaction.

         integer, pointer :: reaction_ye_rho_exponents(:,:)=>NULL() ! (2,rates_reaction_id_max)
            ! multiply T dependent rate by Ye^a(i) * Rho^b(i)
            ! reaction_ye_rho_coeffs(1,i) is a(i)
            ! reaction_ye_rho_coeffs(2,i) is b(i)
            ! (0,0) for photodisintegrations and decays
            ! (0,1) for standard 2 body reactions
            ! (0,2) for 3 body reactions such as triple alpha
            ! (1,1) for 2 body electron captures
            ! (1,2) for 3 body electron captures (e.g., pep)      
      
         integer, parameter :: max_num_reaction_inputs = 3
         integer, pointer :: reaction_inputs(:,:)=>NULL() ! (2*max_num_reaction_inputs,rates_reaction_id_max)
            ! up to max_num_reaction_inputs pairs of coefficients and chem id's, terminated by 0's.
            ! e.g.,  o16(p,g)f17 would be (/ 1, io16, 1, ih1, 0 /)
            ! triple alpha would be (/ 3, ihe4, 0 /)
            ! he3(he4, g)be7(e-,nu)li7(p,a)he4 would be (/ 1, ihe3, 1, ihe4, i, ih1, 0 /)

         integer, parameter :: max_num_reaction_outputs = 4
         integer, pointer :: reaction_outputs(:,:)=>NULL() ! (2*max_num_reaction_outputs,rates_reaction_id_max)
            ! up to max_num_reaction_outputs pairs of coefficients and chem id's, terminated by 0's.
            ! e.g.,  o16(p,g)f17 would be (/ 1, if17, 0 /)
            ! c12(a, p)n15 would be (/ 1, in15, 1, ih1, 0 /)
         
      ! weak_info_list
      
         integer :: num_weak_info_list_reactions
         real(dp), pointer :: weak_info_list_halflife(:)=>NULL(), weak_info_list_Qneu(:)=>NULL()
         type (integer_dict), pointer :: weak_info_list_dict=>NULL()

      ! weak
      
         real(dp) :: &
            T9_weaklib_full_off = 0.01d0, & ! use pure reaclib for T <= this
            T9_weaklib_full_on = 0.02d0 ! use pure weaklib for T >= this
            ! blend for intermediate temperatures

         ! for high Z elements, switch to reaclib at temp where no longer fully ionized
         ! as rough approximation for this, we switch at Fe to higher values of T9
         integer :: weaklib_blend_hi_Z = 26
         ! if input element has Z >= weaklib_blend_hi_Z, then use the following T9 limits
         real(dp) :: &
            T9_weaklib_full_off_hi_Z = 0.063d0, & ! use pure reaclib for T <= this
            T9_weaklib_full_on_hi_Z = 0.073d0 ! use pure weaklib for T >= this

         integer :: num_weak_reactions

         type, abstract :: weak_rate_table
            integer :: num_T9
            integer :: num_lYeRho

            real(dp), allocatable :: T9s(:)
            real(dp), allocatable :: lYeRhos(:)

            logical :: has_cc ! are there tabulated coulomb corrections?
            real(dp), allocatable :: data(:,:,:,:)
               ! (4, num_T9, num_lYeRho, ?)
               ! ? = 3 without coulomb corrections
               ! ? = 5 with coulomb corrections
            contains
            
              procedure(setup_weak_table), deferred :: setup
              procedure(interpolate_weak_table), deferred :: interpolate

         end type weak_rate_table

         abstract interface
            subroutine setup_weak_table(table, ierr)
              import weak_rate_table
              class(weak_rate_table), intent(inout) :: table
              integer, intent(out) :: ierr
            end subroutine setup_weak_table
         end interface

         abstract interface
            subroutine interpolate_weak_table(table, T9, lYeRho, &
                 lambda, dlambda_dlnT, dlambda_dlnRho, &
                 Qneu, dQneu_dlnT, dQneu_dlnRho, &
                 delta_Q, Vs, ierr)
              use const_def, only : dp
              import weak_rate_table
              class(weak_rate_table), intent(inout) :: table
              real(dp), intent(in) :: T9, lYeRho
              real(dp), intent(out) :: lambda, dlambda_dlnT, dlambda_dlnRho
              real(dp), intent(out) :: Qneu, dQneu_dlnT, dQneu_dlnRho
              real(dp), intent(out) :: delta_Q, Vs
              integer, intent(out) :: ierr
            end subroutine interpolate_weak_table
         end interface

         type :: table_c
            class(weak_rate_table), pointer :: t=>NULL()
         end type table_c

         type(table_c), dimension(:), allocatable :: weak_reactions_tables

         integer, pointer, dimension(:) :: & ! (num_weak_reactions)
            weak_lhs_nuclide_id=>NULL(), weak_rhs_nuclide_id=>NULL(), weak_reaclib_id=>NULL()
         character(len=iso_name_length), dimension(:), pointer :: &
            weak_lhs_nuclide_name=>NULL(), weak_rhs_nuclide_name=>NULL() ! (num_weak_reactions)
         type (integer_dict), pointer :: weak_reactions_dict=>NULL()

         logical :: weak_bicubic = .false.  
            ! true means do bicubic splines for interpolation
            ! false means just do bilinear
            ! bilinear is safe; bicubic can overshoot near jumps

         ! Suzuki et al. (2016)
         logical :: use_suzuki_tables = .false.
         
      ! ecapture
      
      logical :: do_ecapture = .false.

      type (integer_dict), pointer :: ecapture_states_number_dict=>NULL()
      type (integer_dict), pointer :: ecapture_states_offset_dict=>NULL()

      type (integer_dict), pointer :: ecapture_transitions_number_dict=>NULL()
      type (integer_dict), pointer :: ecapture_transitions_offset_dict=>NULL()

      integer, parameter :: max_num_ecapture_nuclei = 25
      integer, parameter :: max_num_ecapture_states = 25

      integer, parameter :: max_num_ecapture_reactions = 25
      integer, parameter :: max_num_ecapture_transitions = 25

      integer, parameter :: num_transitions_data = 2
      integer, parameter :: i_Si = 1, i_Sf = 2
      integer :: num_ecapture_transitions
      integer, pointer :: ecapture_transitions_data(:,:)=>NULL()
      real(dp), pointer :: ecapture_logft_data(:)=>NULL()

      integer, parameter :: num_states_data = 2
      integer, parameter :: i_E = 1, i_J = 2
      integer :: num_ecapture_states
      real(dp), pointer :: ecapture_states_data(:,:)=>NULL()

      integer :: num_ecapture_nuclei, num_ecapture_reactions
      integer, pointer, dimension(:) :: ecapture_nuclide_id=>NULL(), &
             ecapture_lhs_nuclide_id=>NULL(), ecapture_rhs_nuclide_id=>NULL() ! (num_ecapture_reactions)
      character(len=iso_name_length), dimension(:), pointer :: ecapture_nuclide_name=>NULL(), &
             ecapture_lhs_nuclide_name=>NULL(), ecapture_rhs_nuclide_name=>NULL() ! (num_ecapture_reactions)
      type (integer_dict), pointer :: ecapture_reactions_dict=>NULL()
      

      integer, pointer :: reaction_categories(:)=>NULL() ! (rates_reaction_id_max) set by net using reactions.list info
      

      integer, pointer,dimension(:) :: &
         reaction_is_reverse, reaction_reaclib_lo, reaction_reaclib_hi, reverse_reaction_id
         ! caches for get_reaclib_rate_and_dlnT (in raw_rates.f)
         ! for all of these, 0 means "cache entry not yet set -- don't have the information"



         
      ! for tabular evaluation of the raw reaction rates
      real(dp) :: rattab_thi != 10.301029995664d0 ! log10(highest temp = 2e10)
      real(dp) :: rattab_tlo != 5.30102999566398d0 ! log10(lowest temp = 2e5)
      real(dp) :: rattab_temp_hi != 10**rattab_thi
      real(dp) :: rattab_temp_lo != 10**rattab_tlo
      
      integer :: rattab_points_per_decade = 2000
      integer :: nrattab ! number of reaction rate table temperatures
         ! nrattab = <points per decade>*(rattab_thi - rattab_tlo) + 1
               
      real(dp) :: rattab_tstp != (rattab_thi-rattab_tlo)/(nrattab-1)! step size
               
      ! reactions for hardwired nets and reactions with multiple choices for rates
         integer, parameter :: ir1212 = 1
         integer, parameter :: ir1216 = ir1212+1
         integer, parameter :: ir1216_to_mg24 = ir1216+1
         integer, parameter :: ir1216_to_si28 = ir1216_to_mg24+1
         integer, parameter :: ir1616 = ir1216_to_si28+1
         integer, parameter :: ir1616a = ir1616+1
         integer, parameter :: ir1616g = ir1616a+1
         integer, parameter :: ir1616p_aux = ir1616g+1
         integer, parameter :: ir1616ppa = ir1616p_aux+1
         integer, parameter :: ir1616ppg = ir1616ppa+1
         integer, parameter :: ir_he3_ag_be7 = ir1616ppg+1
         integer, parameter :: ir34_pp2 = ir_he3_ag_be7+1
         integer, parameter :: ir34_pp3 = ir34_pp2+1
         integer, parameter :: ir_al27_pa_mg24 = ir34_pp3+1
         integer, parameter :: ir_ar36_ag_ca40 = ir_al27_pa_mg24+1
         integer, parameter :: ir_ar36_ga_s32 = ir_ar36_ag_ca40+1
         integer, parameter :: ir_b8_gp_be7 = ir_ar36_ga_s32+1
         integer, parameter :: ir_be7_pg_b8 = ir_b8_gp_be7+1
         integer, parameter :: ir_c12_ag_o16 = ir_be7_pg_b8+1
         integer, parameter :: ir_c12_ap_n15 = ir_c12_ag_o16+1
         integer, parameter :: ir_c12_pg_n13 = ir_c12_ap_n15+1
         integer, parameter :: ir_c12_to_he4_he4_he4 = ir_c12_pg_n13+1
         integer, parameter :: ir_c13_an_o16 = ir_c12_to_he4_he4_he4+1
         integer, parameter :: ir_c13_pg_n14 = ir_c13_an_o16+1
         integer, parameter :: ir_ca40_ag_ti44 = ir_c13_pg_n14+1
         integer, parameter :: ir_ca40_ga_ar36 = ir_ca40_ag_ti44+1
         integer, parameter :: ir_cr48_ag_fe52 = ir_ca40_ga_ar36+1
         integer, parameter :: ir_cr48_ga_ti44 = ir_cr48_ag_fe52+1
         integer, parameter :: ir_f17_ap_ne20 = ir_cr48_ga_ti44+1
         integer, parameter :: ir_f17_gp_o16 = ir_f17_ap_ne20+1
         integer, parameter :: ir_f17_pa_o14 = ir_f17_gp_o16+1
         integer, parameter :: ir_f18_gp_o17 = ir_f17_pa_o14+1
         integer, parameter :: ir_f18_pa_o15 = ir_f18_gp_o17+1

         integer, parameter :: ir_f17_pg_ne18 = ir_f18_pa_o15+1
         integer, parameter :: ir_f18_pg_ne19 = ir_f17_pg_ne18+1

         integer, parameter :: ir_f19_ap_ne22 = ir_f18_pg_ne19+1
         integer, parameter :: ir_f19_gp_o18 = ir_f19_ap_ne22+1
         integer, parameter :: ir_f19_pa_o16 = ir_f19_gp_o18+1
         integer, parameter :: ir_f19_pg_ne20 = ir_f19_pa_o16+1
         integer, parameter :: ir_fe52_ag_ni56 = ir_f19_pg_ne20+1
         integer, parameter :: ir_fe52_ga_cr48 = ir_fe52_ag_ni56+1
         integer, parameter :: ir_h2_be7_to_h1_he4_he4 = ir_fe52_ga_cr48+1
         integer, parameter :: ir_h2_h2_to_he4 = ir_h2_be7_to_h1_he4_he4+1
         integer, parameter :: ir_h2_he3_to_h1_he4 = ir_h2_h2_to_he4+1

         integer, parameter :: ir_he3_be7_to_h1_h1_he4_he4 = ir_h2_he3_to_h1_he4+1
         integer, parameter :: ir_he3_he3_to_h1_h1_he4 = ir_he3_be7_to_h1_h1_he4_he4+1
         integer, parameter :: ir_h1_h1_he4_to_he3_he3 = ir_he3_he3_to_h1_h1_he4+1
         integer, parameter :: ir_he4_he4_he4_to_c12 = ir_h1_h1_he4_to_he3_he3+1
         integer, parameter :: ir_li7_pa_he4 = ir_he4_he4_he4_to_c12+1
         integer, parameter :: ir_mg24_ag_si28 = ir_li7_pa_he4+1
         integer, parameter :: ir_mg24_ap_al27 = ir_mg24_ag_si28+1
         integer, parameter :: ir_mg24_ga_ne20 = ir_mg24_ap_al27+1
         integer, parameter :: ir_n13_ap_o16 = ir_mg24_ga_ne20+1
         integer, parameter :: ir_n13_gp_c12 = ir_n13_ap_o16+1
         integer, parameter :: ir_n13_pg_o14 = ir_n13_gp_c12+1
         integer, parameter :: ir_n14_ag_f18 = ir_n13_pg_o14+1
         integer, parameter :: ir_n14_ap_o17 = ir_n14_ag_f18+1
         integer, parameter :: ir_n14_gp_c13 = ir_n14_ap_o17+1
         integer, parameter :: ir_n14_pg_o15 = ir_n14_gp_c13+1
         integer, parameter :: ir_n15_ag_f19 = ir_n14_pg_o15+1
         integer, parameter :: ir_n15_ap_o18 = ir_n15_ag_f19+1
         integer, parameter :: ir_n15_pa_c12 = ir_n15_ap_o18+1
         integer, parameter :: ir_n15_pg_o16 = ir_n15_pa_c12+1
         integer, parameter :: ir_na23_pa_ne20 = ir_n15_pg_o16+1
         integer, parameter :: ir_ne18_gp_f17 = ir_na23_pa_ne20+1
         integer, parameter :: ir_ne19_ga_o15 = ir_ne18_gp_f17+1
         integer, parameter :: ir_ne19_gp_f18 = ir_ne19_ga_o15+1
         integer, parameter :: ir_ne20_ag_mg24 = ir_ne19_gp_f18+1
         integer, parameter :: ir_ne20_ap_na23 = ir_ne20_ag_mg24+1
         integer, parameter :: ir_ne20_ga_o16 = ir_ne20_ap_na23+1
         integer, parameter :: ir_ne20_gp_f19 = ir_ne20_ga_o16+1
         integer, parameter :: ir_ne22_ag_mg26 = ir_ne20_gp_f19+1
         integer, parameter :: ir_ne22_pg_na23 = ir_ne22_ag_mg26+1
         integer, parameter :: ir_ni56_ga_fe52 = ir_ne22_pg_na23+1
         integer, parameter :: ir_o14_ag_ne18 = ir_ni56_ga_fe52+1
         integer, parameter :: ir_o14_ap_f17 = ir_o14_ag_ne18+1
         integer, parameter :: ir_o14_gp_n13 = ir_o14_ap_f17+1
         integer, parameter :: ir_o15_ag_ne19 = ir_o14_gp_n13+1
         integer, parameter :: ir_o15_ap_f18 = ir_o15_ag_ne19+1
         integer, parameter :: ir_o15_gp_n14 = ir_o15_ap_f18+1
         integer, parameter :: ir_o16_ag_ne20 = ir_o15_gp_n14+1
         integer, parameter :: ir_o16_ap_f19 = ir_o16_ag_ne20+1
         integer, parameter :: ir_o16_ga_c12 = ir_o16_ap_f19+1
         integer, parameter :: ir_o16_gp_n15 = ir_o16_ga_c12+1
         integer, parameter :: ir_o16_pg_f17 = ir_o16_gp_n15+1
         integer, parameter :: ir_o17_pa_n14 = ir_o16_pg_f17+1
         integer, parameter :: ir_o17_pg_f18 = ir_o17_pa_n14+1
         integer, parameter :: ir_o18_ag_ne22 = ir_o17_pg_f18+1
         integer, parameter :: ir_o18_pa_n15 = ir_o18_ag_ne22+1
         integer, parameter :: ir_o18_pg_f19 = ir_o18_pa_n15+1
         integer, parameter :: ir_s32_ag_ar36 = ir_o18_pg_f19+1
         integer, parameter :: ir_s32_ga_si28 = ir_s32_ag_ar36+1
         integer, parameter :: ir_si28_ag_s32 = ir_s32_ga_si28+1
         integer, parameter :: ir_si28_ga_mg24 = ir_si28_ag_s32+1
         integer, parameter :: ir_ti44_ag_cr48 = ir_si28_ga_mg24+1
         integer, parameter :: ir_ti44_ga_ca40 = ir_ti44_ag_cr48+1
         integer, parameter :: iral27pa_aux = ir_ti44_ga_ca40+1
         integer, parameter :: iral27pg_aux = iral27pa_aux+1
         integer, parameter :: irar36ap_aux = iral27pg_aux+1
         integer, parameter :: irar36ap_to_ca40 = irar36ap_aux+1
         integer, parameter :: irar36gp_aux = irar36ap_to_ca40+1
         integer, parameter :: irar36gp_to_s32 = irar36gp_aux+1
         integer, parameter :: irbe7ec_li7_aux = irar36gp_to_s32+1
         integer, parameter :: irbe7pg_b8_aux = irbe7ec_li7_aux+1
         integer, parameter :: irc12_to_c13 = irbe7pg_b8_aux+1
         integer, parameter :: irc12_to_n14 = irc12_to_c13+1
         integer, parameter :: irc12ap_aux = irc12_to_n14+1
         integer, parameter :: irc12ap_to_o16 = irc12ap_aux+1
         integer, parameter :: irca40ap_aux = irc12ap_to_o16+1
         integer, parameter :: irca40ap_to_ti44 = irca40ap_aux+1
         integer, parameter :: irca40gp_aux = irca40ap_to_ti44+1
         integer, parameter :: irca40gp_to_ar36 = irca40gp_aux+1
         integer, parameter :: ircl35pa_aux = irca40gp_to_ar36+1
         integer, parameter :: ircl35pg_aux = ircl35pa_aux+1
         integer, parameter :: irco55gprot_aux = ircl35pg_aux+1
         integer, parameter :: irco55pg_aux = irco55gprot_aux+1
         integer, parameter :: irco55protg_aux = irco55pg_aux+1
         integer, parameter :: ircr48ap_aux = irco55protg_aux+1
         integer, parameter :: ircr48ap_to_fe52 = ircr48ap_aux+1
         integer, parameter :: ircr48gp_aux = ircr48ap_to_fe52+1
         integer, parameter :: ircr48gp_to_ti44 = ircr48gp_aux+1
         integer, parameter :: irf19pg_aux = ircr48gp_to_ti44+1
         integer, parameter :: irfe52ap_aux = irf19pg_aux+1
         integer, parameter :: irfe52ap_to_ni56 = irfe52ap_aux+1
         integer, parameter :: irfe52aprot_aux = irfe52ap_to_ni56+1
         integer, parameter :: irfe52aprot_to_fe54 = irfe52aprot_aux+1
         integer, parameter :: irfe52aprot_to_ni56 = irfe52aprot_to_fe54+1
         integer, parameter :: irfe52gp_aux = irfe52aprot_to_ni56+1
         integer, parameter :: irfe52gp_to_cr48 = irfe52gp_aux+1
         integer, parameter :: irfe52neut_to_fe54 = irfe52gp_to_cr48+1
         integer, parameter :: irfe52ng_aux = irfe52neut_to_fe54+1
         integer, parameter :: irfe53gn_aux = irfe52ng_aux+1
         integer, parameter :: irfe53ng_aux = irfe53gn_aux+1
         integer, parameter :: irfe54a_to_ni56 = irfe53ng_aux+1
         integer, parameter :: irfe54an_aux = irfe54a_to_ni56+1
         integer, parameter :: irfe54an_to_ni56 = irfe54an_aux+1
         integer, parameter :: irfe54aprot_to_fe56 = irfe54an_to_ni56+1
         integer, parameter :: irfe54g_to_fe52 = irfe54aprot_to_fe56+1
         integer, parameter :: irfe54ng_aux = irfe54g_to_fe52+1
         integer, parameter :: irfe54ng_to_fe56 = irfe54ng_aux+1
         integer, parameter :: irfe54prot_to_fe52 = irfe54ng_to_fe56+1
         integer, parameter :: irfe54prot_to_ni56 = irfe54prot_to_fe52+1
         integer, parameter :: irfe54protg_aux = irfe54prot_to_ni56+1
         integer, parameter :: irfe55gn_aux = irfe54protg_aux+1
         integer, parameter :: irfe55ng_aux = irfe55gn_aux+1
         
         integer, parameter :: irfe56ec_fake_to_mn56 = irfe55ng_aux+1
         integer, parameter :: irfe56ec_fake_to_mn57 = irfe56ec_fake_to_mn56+1
         integer, parameter :: irfe56ec_fake_to_cr56 = irfe56ec_fake_to_mn57+1
         integer, parameter :: irfe56ec_fake_to_cr57 = irfe56ec_fake_to_cr56+1
         integer, parameter :: irfe56ec_fake_to_cr58 = irfe56ec_fake_to_cr57+1
         integer, parameter :: irfe56ec_fake_to_cr59 = irfe56ec_fake_to_cr58+1
         integer, parameter :: irfe56ec_fake_to_cr60 = irfe56ec_fake_to_cr59+1
         integer, parameter :: irfe56ec_fake_to_cr61 = irfe56ec_fake_to_cr60+1
         integer, parameter :: irfe56ec_fake_to_cr62 = irfe56ec_fake_to_cr61+1
         integer, parameter :: irfe56ec_fake_to_cr63 = irfe56ec_fake_to_cr62+1
         integer, parameter :: irfe56ec_fake_to_cr64 = irfe56ec_fake_to_cr63+1
         integer, parameter :: irfe56ec_fake_to_cr65 = irfe56ec_fake_to_cr64+1
         integer, parameter :: irfe56ec_fake_to_cr66 = irfe56ec_fake_to_cr65+1
         
         integer, parameter :: irfe56ee_to_ni56 = irfe56ec_fake_to_cr66+1
         integer, parameter :: irfe56gn_aux = irfe56ee_to_ni56+1
         integer, parameter :: irfe56gn_to_fe54 = irfe56gn_aux+1
         integer, parameter :: irfe56prot_to_fe54 = irfe56gn_to_fe54+1
         integer, parameter :: irh2_protg_aux = irfe56prot_to_fe54+1
         integer, parameter :: irh2g_neut_aux = irh2_protg_aux+1
         integer, parameter :: irhe3_neutg_aux = irh2g_neut_aux+1
         integer, parameter :: irhe3gprot_aux = irhe3_neutg_aux+1
         integer, parameter :: irhe4_breakup = irhe3gprot_aux+1
         integer, parameter :: irhe4_rebuild = irhe4_breakup+1
         integer, parameter :: irhe4g_neut_aux = irhe4_rebuild+1
         integer, parameter :: irk39pa_aux = irhe4g_neut_aux+1
         integer, parameter :: irk39pg_aux = irk39pa_aux+1
         integer, parameter :: irmg24ap_aux = irk39pg_aux+1
         integer, parameter :: irmg24ap_to_si28 = irmg24ap_aux+1
         integer, parameter :: irmg24gp_aux = irmg24ap_to_si28+1
         integer, parameter :: irmg24gp_to_ne20 = irmg24gp_aux+1
         integer, parameter :: irmn51pg_aux = irmg24gp_to_ne20+1
         integer, parameter :: irn14_to_c12 = irmn51pg_aux+1
         integer, parameter :: irn14_to_n15 = irn14_to_c12+1
         integer, parameter :: irn14_to_o16 = irn14_to_n15+1
         integer, parameter :: irn14ag_lite = irn14_to_o16+1
         integer, parameter :: irn14gc12 = irn14ag_lite+1
         integer, parameter :: irn14pg_aux = irn14gc12+1
         integer, parameter :: irn15pa_aux = irn14pg_aux+1
         integer, parameter :: irn15pg_aux = irn15pa_aux+1
         integer, parameter :: irna23pa_aux = irn15pg_aux+1
         integer, parameter :: irna23pg_aux = irna23pa_aux+1
         integer, parameter :: irne18ag_to_mg24 = irna23pg_aux+1
         integer, parameter :: irne18ap_to_mg22 = irne18ag_to_mg24+1
         integer, parameter :: irne18ap_to_mg24 = irne18ap_to_mg22+1
         integer, parameter :: irne19pg_to_mg22 = irne18ap_to_mg24+1
         integer, parameter :: irne19pg_to_mg24 = irne19pg_to_mg22+1
         integer, parameter :: irne20ap_aux = irne19pg_to_mg24+1
         integer, parameter :: irne20ap_to_mg24 = irne20ap_aux+1
         integer, parameter :: irne20gp_aux = irne20ap_to_mg24+1
         integer, parameter :: irne20gp_to_o16 = irne20gp_aux+1
         integer, parameter :: irne20pg_to_mg22 = irne20gp_to_o16+1
         integer, parameter :: irne20pg_to_mg24 = irne20pg_to_mg22+1
         integer, parameter :: irneut_to_prot = irne20pg_to_mg24+1
         integer, parameter :: irni56ec_to_fe54 = irneut_to_prot+1
         integer, parameter :: irni56ec_to_fe56 = irni56ec_to_fe54+1
         integer, parameter :: irni56ec_to_co56 = irni56ec_to_fe56+1
         integer, parameter :: irco56ec_to_fe56 = irni56ec_to_co56+1
         integer, parameter :: irni56gp_aux = irco56ec_to_fe56+1
         integer, parameter :: irni56gp_to_fe52 = irni56gp_aux+1
         integer, parameter :: irni56gprot_aux = irni56gp_to_fe52+1
         integer, parameter :: irni56gprot_to_fe52 = irni56gprot_aux+1
         integer, parameter :: irni56gprot_to_fe54 = irni56gprot_to_fe52+1
         integer, parameter :: irni56ng_to_fe54 = irni56gprot_to_fe54+1
         integer, parameter :: irni57na_aux = irni56ng_to_fe54+1
         integer, parameter :: iro16_to_n14 = irni57na_aux+1
         integer, parameter :: iro16_to_o17 = iro16_to_n14+1
         integer, parameter :: iro16ap_aux = iro16_to_o17+1
         integer, parameter :: iro16ap_to_ne20 = iro16ap_aux+1
         integer, parameter :: iro16gp_aux = iro16ap_to_ne20+1
         integer, parameter :: iro16gp_to_c12 = iro16gp_aux+1
         integer, parameter :: iro17_to_o18 = iro16gp_to_c12+1
         integer, parameter :: irp31pa_aux = iro17_to_o18+1
         integer, parameter :: irp31pg_aux = irp31pa_aux+1
         integer, parameter :: irpep_to_he3 = irp31pg_aux+1
         integer, parameter :: irpp_to_he3 = irpep_to_he3+1
         integer, parameter :: irprot_neutg_aux = irpp_to_he3+1
         integer, parameter :: irprot_to_neut = irprot_neutg_aux+1
         integer, parameter :: irs32ap_aux = irprot_to_neut+1
         integer, parameter :: irs32ap_to_ar36 = irs32ap_aux+1
         integer, parameter :: irs32gp_aux = irs32ap_to_ar36+1
         integer, parameter :: irs32gp_to_si28 = irs32gp_aux+1
         integer, parameter :: irsc43pa_aux = irs32gp_to_si28+1
         integer, parameter :: irsc43pg_aux = irsc43pa_aux+1
         integer, parameter :: irsi28ap_aux = irsc43pg_aux+1
         integer, parameter :: irsi28ap_to_s32 = irsi28ap_aux+1
         integer, parameter :: irsi28gp_aux = irsi28ap_to_s32+1
         integer, parameter :: irsi28gp_to_mg24 = irsi28gp_aux+1
         integer, parameter :: irti44ap_aux = irsi28gp_to_mg24+1
         integer, parameter :: irti44ap_to_cr48 = irti44ap_aux+1
         integer, parameter :: irti44gp_aux = irti44ap_to_cr48+1
         integer, parameter :: irti44gp_to_ca40 = irti44gp_aux+1
         integer, parameter :: irv47pa_aux = irti44gp_to_ca40+1
         integer, parameter :: irv47pg_aux = irv47pa_aux+1
         integer, parameter :: ir_h1_h1_wk_h2 = irv47pg_aux+1
         integer, parameter :: ir_h1_h1_ec_h2 = ir_h1_h1_wk_h2+1
         integer, parameter :: irn14ag_to_ne22 = ir_h1_h1_ec_h2+1
         integer, parameter :: irf19pa_aux = irn14ag_to_ne22+1
         integer, parameter :: ir_b8_wk_he4_he4 = irf19pa_aux+1
         integer, parameter :: irmn51pa_aux = ir_b8_wk_he4_he4+1
         integer, parameter :: irfe54gn_aux = irmn51pa_aux+1
         integer, parameter :: irco55pa_aux = irfe54gn_aux+1
         integer, parameter :: irco55prota_aux = irco55pa_aux+1
         integer, parameter :: irn14ag_to_o18 = irco55prota_aux+1
         integer, parameter :: ir_h1_he3_wk_he4 = irn14ag_to_o18+1
         integer, parameter :: ir_be7_wk_li7 = ir_h1_he3_wk_he4+1

         ! rates added for approx21
         integer, parameter :: ir_al27_pg_si28 = ir_be7_wk_li7+1
         integer, parameter :: ir_si28_gp_al27 = ir_al27_pg_si28+1
         integer, parameter :: ir_si28_ap_p31 = ir_si28_gp_al27+1
         integer, parameter :: ir_p31_pa_si28 = ir_si28_ap_p31+1
         integer, parameter :: ir_p31_pg_s32 = ir_p31_pa_si28+1
         integer, parameter :: ir_s32_gp_p31 = ir_p31_pg_s32+1
         integer, parameter :: ir_s32_ap_cl35 = ir_s32_gp_p31+1
         integer, parameter :: ir_cl35_pa_s32 = ir_s32_ap_cl35+1
         integer, parameter :: ir_cl35_pg_ar36 = ir_cl35_pa_s32+1
         integer, parameter :: ir_ar36_gp_cl35 = ir_cl35_pg_ar36+1
         integer, parameter :: ir_ar36_ap_k39 = ir_ar36_gp_cl35+1
         integer, parameter :: ir_k39_pa_ar36 = ir_ar36_ap_k39+1
         integer, parameter :: ir_k39_pg_ca40 = ir_k39_pa_ar36+1
         integer, parameter :: ir_ca40_gp_k39 = ir_k39_pg_ca40+1
         integer, parameter :: ir_ca40_ap_sc43 = ir_ca40_gp_k39+1
         integer, parameter :: ir_sc43_pa_ca40 = ir_ca40_ap_sc43+1
         integer, parameter :: ir_sc43_pg_ti44 = ir_sc43_pa_ca40+1
         integer, parameter :: ir_ti44_gp_sc43 = ir_sc43_pg_ti44+1
         integer, parameter :: ir_ti44_ap_v47 = ir_ti44_gp_sc43+1
         integer, parameter :: ir_v47_pa_ti44 = ir_ti44_ap_v47+1
         integer, parameter :: ir_v47_pg_cr48 = ir_v47_pa_ti44+1
         integer, parameter :: ir_cr48_gp_v47 = ir_v47_pg_cr48+1
         integer, parameter :: ir_cr48_ap_mn51 = ir_cr48_gp_v47+1
         integer, parameter :: ir_mn51_pa_cr48 = ir_cr48_ap_mn51+1
         integer, parameter :: ir_mn51_pg_fe52 = ir_mn51_pa_cr48+1
         integer, parameter :: ir_fe52_gp_mn51 = ir_mn51_pg_fe52+1
         integer, parameter :: ir_fe52_ap_co55 = ir_fe52_gp_mn51+1
         integer, parameter :: ir_co55_pa_fe52 = ir_fe52_ap_co55+1
         integer, parameter :: ir_co55_pg_ni56 = ir_co55_pa_fe52+1
         integer, parameter :: ir_ni56_gp_co55 = ir_co55_pg_ni56+1
         integer, parameter :: ir_fe52_ng_fe53 = ir_ni56_gp_co55+1
         integer, parameter :: ir_fe53_gn_fe52 = ir_fe52_ng_fe53+1
         integer, parameter :: ir_fe53_ng_fe54 = ir_fe53_gn_fe52+1
         integer, parameter :: ir_fe54_gn_fe53 = ir_fe53_ng_fe54+1
         integer, parameter :: ir_fe54_pg_co55 = ir_fe54_gn_fe53+1
         integer, parameter :: ir_co55_gp_fe54 = ir_fe54_pg_co55+1
         integer, parameter :: ir_he3_ng_he4 = ir_co55_gp_fe54+1
         integer, parameter :: ir_he4_gn_he3 = ir_he3_ng_he4+1
         integer, parameter :: ir_h1_ng_h2 = ir_he4_gn_he3+1
         integer, parameter :: ir_h2_gn_h1 = ir_h1_ng_h2+1
         integer, parameter :: ir_h2_pg_he3 = ir_h2_gn_h1+1
         integer, parameter :: ir_he3_gp_h2 = ir_h2_pg_he3+1
         integer, parameter :: ir_fe54_ng_fe55 = ir_he3_gp_h2+1
         integer, parameter :: ir_fe55_gn_fe54 = ir_fe54_ng_fe55+1
         integer, parameter :: ir_fe55_ng_fe56 = ir_fe55_gn_fe54+1
         integer, parameter :: ir_fe56_gn_fe55 = ir_fe55_ng_fe56+1
         integer, parameter :: ir_fe54_ap_co57 = ir_fe56_gn_fe55+1
         integer, parameter :: ir_co57_pa_fe54 = ir_fe54_ap_co57+1
         integer, parameter :: ir_fe56_pg_co57 = ir_co57_pa_fe54+1
         integer, parameter :: ir_co57_gp_fe56 = ir_fe56_pg_co57+1
         
         integer, parameter :: ir_c12_c12_to_h1_na23 = ir_co57_gp_fe56+1
         integer, parameter :: ir_he4_ne20_to_c12_c12 = ir_c12_c12_to_h1_na23+1
         integer, parameter :: ir_c12_c12_to_he4_ne20 = ir_he4_ne20_to_c12_c12+1
         integer, parameter :: ir_he4_mg24_to_c12_o16 = ir_c12_c12_to_he4_ne20+1


! fxt for al26 isomers
         integer, parameter :: ir_al26_1_to_al26_2 = ir_he4_mg24_to_c12_o16 + 1
         integer, parameter :: ir_al26_2_to_al26_1 = ir_al26_1_to_al26_2 + 1


         integer, parameter :: num_predefined_reactions = ir_al26_2_to_al26_1

         integer :: rates_reaction_id_max


      
      ! for mazurek's ni56 electron capture rate interpolation
         real(dp) :: tv(7),rv(6),rfdm(4),rfd0(4),rfd1(4),rfd2(4),tfdm(5),tfd0(5),tfd1(5),tfd2(5)


      type T_Factors
         real(dp) :: lnT9
         real(dp) :: T9
         real(dp) :: T92
         real(dp) :: T93
         real(dp) :: T94
         real(dp) :: T95
         real(dp) :: T96
         real(dp) :: T912
         real(dp) :: T932
         real(dp) :: T952
         real(dp) :: T972
         real(dp) :: T913
         real(dp) :: T923
         real(dp) :: T943
         real(dp) :: T953
         real(dp) :: T973
         real(dp) :: T9113
         real(dp) :: T914
         real(dp) :: T934
         real(dp) :: T954
         real(dp) :: T974
         real(dp) :: T915
         real(dp) :: T935
         real(dp) :: T945
         real(dp) :: T965
         real(dp) :: T917
         real(dp) :: T927
         real(dp) :: T947
         real(dp) :: T918
         real(dp) :: T938
         real(dp) :: T958
         real(dp) :: T9i
         real(dp) :: T9i2
         real(dp) :: T9i3
         real(dp) :: T9i12
         real(dp) :: T9i32
         real(dp) :: T9i52
         real(dp) :: T9i72
         real(dp) :: T9i13
         real(dp) :: T9i23
         real(dp) :: T9i43
         real(dp) :: T9i53
         real(dp) :: T9i14
         real(dp) :: T9i34
         real(dp) :: T9i54
         real(dp) :: T9i15
         real(dp) :: T9i35
         real(dp) :: T9i45
         real(dp) :: T9i65
         real(dp) :: T9i17
         real(dp) :: T9i27
         real(dp) :: T9i47
         real(dp) :: T9i18
         real(dp) :: T9i38
         real(dp) :: T9i58
         real(dp) :: T916
         real(dp) :: T976
         real(dp) :: T9i76
      end type T_Factors
      
      
      
      ! rate results components
      
      integer, parameter :: i_rate = 1        
      integer, parameter :: i_rate_dT = 2  
      integer, parameter :: i_rate_dRho = 3 
      integer, parameter :: num_rvs = 3
      
      
      
      ! screening
      
      integer, parameter :: no_screening = 0
      integer, parameter :: classic_screening = 1
         ! DeWitt, Graboske, Cooper, "Screening Factors for Nuclear Reactions. 
         !    I. General Theory", ApJ, 181:439-456, 1973.
         ! Graboske, DeWitt, Grossman, Cooper, "Screening Factors for Nuclear Reactions. 
         !    II. Intermediate Screening and Astrophysical Applications", ApJ, 181:457-474, 1973.
      integer, parameter :: extended_screening = 2
         ! based on code from Frank Timmes
         ! extends the Graboske method using results from Alastuey and Jancovici (1978),
         ! along with plasma parameters from Itoh et al (1979) for strong screening.
      integer, parameter :: salpeter_screening = 3
         ! weak screening only.  following Salpeter (1954),
         ! with equations (4-215) and (4-221) of Clayton (1968).
      integer, parameter :: chugunov_screening = 4
        ! based on code from Sam Jones
        ! Implements screening from Chugunov et al (2007) 

      type Screen_Info
         real(dp) :: temp
         real(dp) :: den
         real(dp) :: logT
         real(dp) :: logRho
         real(dp) :: theta_e
         real(dp) :: zbar
         real(dp) :: abar
         real(dp) :: z2bar
         real(dp) :: zbar13
         real(dp) :: zbar0pt28
         real(dp) :: z1pt58bar
         real(dp) :: ztilda ! sqrt(z2bar + zbar*theta_e)  ! (Dewitt eqn 4)
         real(dp) :: ztilda0pt58
         real(dp) :: Lambda0 ! = 1.88d8*sqrt(rho/(abar*T**3)) ! (Graboske eqn 19; mu_I = abar)
         real(dp) :: Lambda0b ! Lambda0**0.86
         real(dp) :: Lambda0_23 ! Lambda0**(2d0/3d0)
         real(dp) :: ytot
         real(dp) :: rr
         real(dp) :: tempi
         real(dp) :: dtempi
         real(dp) :: deni
         real(dp) :: pp
         real(dp) :: dppdt
         real(dp) :: dppdd
         real(dp) :: qlam0z
         real(dp) :: qlam0zdt
         real(dp) :: qlam0zdd
         real(dp) :: taufac
         real(dp) :: taufacdt
         real(dp) :: xni
         real(dp) :: dxnidd
         real(dp) :: aa
         real(dp) :: daadt
         real(dp) :: daadd
         real(dp) :: ntot, a_e
         integer :: num_calls, num_cache_hits
      end type Screen_Info
      
      
      real(dp) :: reaclib_min_T9 ! for T9 < this, return 0 for reaclib strong rates
      

      ! integers to 1/3 (for graboske screening)
      integer, parameter :: num_one_thirds = 60
      real(dp) :: one_third_power(num_one_thirds)
      ! integers to 1.86 (for graboske screening)
      integer, parameter :: num_pow_186 = 60
      real(dp) :: pow_186(num_pow_186)


      type (integer_dict), pointer :: reaction_names_dict
      
      logical :: have_finished_initialization = .false.
      logical :: rates_use_cache = .true.

      


      ! choices for various rates
         ! NOTE: if change these, must edit raw_rates to match.
         
         ! NACRE = Angulo et al. 1999 Nucl. Phys. A, 656, 3
         ! JR = jina reaclib -- (Sakharuk et al. 2006)
         ! CF88 = Frank Timmes' version of 
            ! Caughlin, G. R. & Fowler, W. A. 1988, Atom. Data and Nuc. Data Tables, 40, 283
         
         
         ! when possible, NACRE is 1, jina reaclib is 2
         integer, parameter :: rates_NACRE_if_available = 1
         integer, parameter :: rates_JR_if_available = 2
         
         ! the jina reaclib rates are not valid below T = 10^7.
         ! if we have the option, we automatically blend over to nacre for low temperatures.
         ! the following values determine the blend region.
         real(dp) :: JR_T_full_off = 1.0d7 ! don't use JR below this
         real(dp) :: JR_T_full_on = 1.1d7 ! don't need to blend above this
         
                  
         ! triple alpha
         integer, parameter :: use_rate_3a_NACRE = 1
         integer, parameter :: use_rate_3a_JR = 2 
         integer, parameter :: use_rate_3a_CF88 = 3
         integer, parameter :: use_rate_3a_FL87 = 4 ! Fushiki and Lamb, Apj, 317, 368-388, 1987
            ! note: use_rate_3a_FL87 is a special case. see eval_FL_epsnuc_3alf in rate_lib
         
         ! c12(a,g)o16
         integer, parameter :: use_rate_c12ag_NACRE = 1
         integer, parameter :: use_rate_c12ag_JR = 2 
         integer, parameter :: use_rate_c12ag_Kunz = 3 ! Kunz et al. (2002)
         integer, parameter :: use_rate_c12ag_CF88 = 4
         
         ! c12 + c12
         integer, parameter :: use_rate_1212_CF88_multi = 1
            ! combines the rates for the n, p, and a channels.
            ! using neutron branching from dayras switkowski and woosley 1976
            ! and an estimate of proton branching.
         integer, parameter :: use_rate_1212_CF88_basic = 2
            ! the single rate approximation from CF88
      
         ! n14(p,g)o15
         integer, parameter :: use_rate_n14pg_NACRE = 1
         integer, parameter :: use_rate_n14pg_JR = 2
         integer, parameter :: use_rate_n14pg_CF88 = 3
         
         
         ! o16 + o16
         integer, parameter :: use_rate_1616_CF88 = 1
         integer, parameter :: use_rate_1616_reaclib = 2
         

         ! Here be dragons, higher temperatures
         ! will generate possibly un-physical values as the partition table cuts off at 1d10
         ! and the polynomial fits to the rates cuts off at 1d10. We truncate the rates whether
         ! the warn flag is on or off, to stop the truncation set a higher max_safe_logT_for_rates
         
         ! Warn if rates exceed the max usable temperature 
         logical :: warn_rates_for_high_temp = .true.
         real(dp) :: max_safe_logT_for_rates = 10d0
          
         ! Maximum sensible value for a reaction rate
         ! This is to try and catch rates that go bad
         real(dp),parameter :: max_safe_rate_for_any_temp = 1d40 

      ! info for rates being evaluated using tables (rate_list.txt)
      type rate_table_info
         logical :: use_rate_table
         logical :: need_to_read
         character (len=132) :: rate_fname
         integer :: nT8s
         real(dp), pointer :: T8s(:) ! (nT8s)
         real(dp), pointer :: f1(:) ! =(4,nT8s)
      end type rate_table_info
      
      type (rate_table_info), pointer :: raw_rates_records(:)
      character (len=1000) :: rates_table_dir

      type (integer_dict), pointer :: skip_warnings_dict
      
      type (reaction_data), target :: reaclib_rates
         
      character (len=1000) :: rates_dir, rates_cache_dir, rates_temp_cache_dir
      


      ! coulomb corrections for weak reactions
      integer :: which_vs_coulomb = 0
      integer :: which_mui_coulomb = 0

      type Coulomb_Info
         real(dp) :: temp
         real(dp) :: den
         real(dp) :: logT
         real(dp) :: logRho
         real(dp) :: theta_e
         real(dp) :: zbar
         real(dp) :: abar
         real(dp) :: z2bar
         real(dp) :: ye
         real(dp) :: z52bar
         real(dp) :: zbar13
         real(dp) :: abari
         real(dp) :: rr
         real(dp) :: tempi
         real(dp) :: dtempi
         real(dp) :: deni
         real(dp) :: pp
         real(dp) :: rs
         real(dp) :: gamma_e
      end type Coulomb_Info
      
      
      logical :: star_debugging_rates_flag
      real(dp) :: rates_test_partials_val, rates_test_partials_dval_dx
      real(dp) :: rates_test_partials_logT_lo, rates_test_partials_logT_hi
      real(dp) :: rates_test_partials_logRho_lo, rates_test_partials_logRho_hi

      
      contains
      
      
      subroutine set_rates_cache_dir(rates_cache_dir_in, ierr)
         use const_def, only: mesa_data_dir, mesa_caches_dir, mesa_temp_caches_dir, use_mesa_temp_cache
         use utils_lib, only : mkdir, switch_str
         character (len=*), intent(in) :: rates_cache_dir_in
         integer, intent(out) :: ierr
         ierr = 0    
         rates_dir = trim(mesa_data_dir) // '/rates_data'
         if (len_trim(rates_cache_dir_in) > 0) then
            rates_cache_dir = rates_cache_dir_in
         else if (len_trim(mesa_caches_dir) > 0) then
            rates_cache_dir = trim(mesa_caches_dir) // '/rates_cache'
         else
            rates_cache_dir = trim(rates_dir) // '/cache'
         end if
         if (rates_use_cache) call mkdir(rates_cache_dir)
         
         rates_temp_cache_dir=trim(mesa_temp_caches_dir)//'/rates_cache'
         if (use_mesa_temp_cache) call mkdir(rates_temp_cache_dir) 
                         
      end subroutine set_rates_cache_dir


      subroutine start_rates_def_init(ierr)
         use utils_lib, only: integer_dict_define
         use math_lib
         integer, intent(out) :: ierr
         
         integer :: i
         ierr = 0    
         star_debugging_rates_flag = .false.
         call create_skip_warnings_dict(ierr)  
         if (ierr /= 0) return
         nullify(reaction_names_dict)
         do i=1,rates_reaction_id_max
            call integer_dict_define(reaction_names_dict, reaction_Name(i), i, ierr)
            if (ierr /= 0) then
               write(*,*) 'FATAL ERROR: rates_def_init failed in integer_dict_define'
               return
            end if
         end do
         call do_start_rates_def_init(ierr)
         
      end subroutine start_rates_def_init
      
      
      subroutine create_skip_warnings_dict(ierr)
         use utils_lib
         use utils_def
         integer, intent(out) :: ierr
         
         integer :: iounit, n, i, t, id, read_int
         character (len=256) :: buffer, string, filename, list_filename
         
         ierr = 0
         
         nullify(skip_warnings_dict)

         list_filename = 'skip_warnings.list'
         ! first try the local directory
         filename = trim(list_filename)
         open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
         if (ierr /= 0) then ! if don't find that file, look in rates_data
            filename = trim(rates_dir) // '/' // trim(list_filename)
            ierr = 0
            open(newunit=iounit, file=trim(filename), action='read', status='old', iostat=ierr)
            if (ierr /= 0) then
               write(*,*) 'failed to open file ' // trim(list_filename)
               return
            end if
         end if
         
         n = 0
         i = 0
         
      reaction_loop: do
            t = token(iounit, n, i, buffer, string)
            if (t == eof_token) exit
            if (t /= name_token) then
               call error; return
            end if
            call integer_dict_define(skip_warnings_dict, string, 1, ierr)
            if (ierr /= 0) then
               write(*,*) 'FATAL ERROR: create_skip_warnings_dict failed in integer_dict_define'
               return
            end if
         end do reaction_loop
         
         close(iounit)

         call integer_dict_create_hash(skip_warnings_dict, ierr)
         if (ierr /= 0) then
            write(*,*) 'FATAL ERROR: create_skip_warnings_dict failed'
            return
         end if

         contains
         
         subroutine error
            ierr = -1
            close(iounit)
         end subroutine error

      end subroutine create_skip_warnings_dict
      
      
      integer function reaclib_index(handle) result(indx)
         use utils_lib, only: integer_dict_lookup
         character(len=*), intent(in) :: handle ! as in rates% reaction_handle
         integer :: ierr
         ierr = 0
         call integer_dict_lookup(reaclib_rates% reaction_dict, handle, indx, ierr)
         if (ierr /= 0) indx = 0
      end function reaclib_index
      
      
      integer function reaclib_reverse(handle) result(indx)
         use utils_lib, only: integer_dict_lookup
         character(len=*), intent(in) :: handle ! as in rates% reaction_handle
         integer :: ierr
         ierr = 0
         call integer_dict_lookup(reaclib_rates% reverse_dict, handle, indx, ierr)
         if (ierr /= 0) indx = 0
      end function reaclib_reverse
   
   
      subroutine weaklib_init(ierr)
         use const_def, only: mesa_data_dir
         integer, intent(out) :: ierr 
         integer :: i
         ierr = 0
         weak_data_dir = trim(mesa_data_dir) // '/rates_data'
         nullify(weak_reactions_dict)
         nullify(weak_info_list_dict)
      end subroutine weaklib_init
      
      
      subroutine free_weak_info

         use utils_lib, only: integer_dict_free

         integer :: i

         if (ALLOCATED(weak_reactions_tables)) then
            do i = 1, num_weak_reactions
               if (ASSOCIATED(weak_reactions_tables(i)%t)) deallocate(weak_reactions_tables(i)%t)
            end do
            deallocate(weak_reactions_tables)
         endif
            
         if (ASSOCIATED(weak_lhs_nuclide_id)) deallocate(weak_lhs_nuclide_id)
         if (ASSOCIATED(weak_rhs_nuclide_id)) deallocate(weak_rhs_nuclide_id)
         if (ASSOCIATED(weak_lhs_nuclide_name)) deallocate(weak_lhs_nuclide_name)
         if (ASSOCIATED(weak_rhs_nuclide_name)) deallocate(weak_rhs_nuclide_name)
         if (ASSOCIATED(weak_reaclib_id)) deallocate(weak_reaclib_id)

         if (ASSOCIATED(weak_info_list_halflife)) deallocate(weak_info_list_halflife)
         if (ASSOCIATED(weak_info_list_Qneu)) deallocate(weak_info_list_Qneu)

         if (ASSOCIATED(weak_reactions_dict)) call integer_dict_free(weak_reactions_dict)
         if (ASSOCIATED(weak_info_list_dict)) call integer_dict_free(weak_info_list_dict)

      end subroutine free_weak_info


      subroutine ecapture_init(ierr)

         use const_def, only: mesa_data_dir

         integer, intent(out) :: ierr 

         ierr = 0

         nullify(ecapture_reactions_dict)
         nullify(ecapture_transitions_number_dict)
         nullify(ecapture_transitions_offset_dict)
         nullify(ecapture_states_number_dict)
         nullify(ecapture_states_offset_dict)

      end subroutine ecapture_init


      subroutine free_ecapture_info

         use utils_lib, only: integer_dict_free

         if (ASSOCIATED(ecapture_transitions_data)) deallocate(ecapture_transitions_data)
         if (ASSOCIATED(ecapture_states_data)) deallocate(ecapture_states_data)
         if (ASSOCIATED(ecapture_logft_data)) deallocate(ecapture_logft_data)
         if (ASSOCIATED(ecapture_nuclide_id)) deallocate(ecapture_nuclide_id)
         if (ASSOCIATED(ecapture_lhs_nuclide_id)) deallocate(ecapture_lhs_nuclide_id)
         if (ASSOCIATED(ecapture_rhs_nuclide_id)) deallocate(ecapture_rhs_nuclide_id)
         if (ASSOCIATED(ecapture_nuclide_name)) deallocate(ecapture_nuclide_name)
         if (ASSOCIATED(ecapture_lhs_nuclide_name)) deallocate(ecapture_lhs_nuclide_name)
         if (ASSOCIATED(ecapture_rhs_nuclide_name)) deallocate(ecapture_rhs_nuclide_name)

         if (ASSOCIATED(ecapture_reactions_dict)) call integer_dict_free(ecapture_reactions_dict)
         if (ASSOCIATED(ecapture_transitions_number_dict)) call integer_dict_free(ecapture_transitions_number_dict)
         if (ASSOCIATED(ecapture_transitions_offset_dict)) call integer_dict_free(ecapture_transitions_offset_dict)
         if (ASSOCIATED(ecapture_states_number_dict)) call integer_dict_free(ecapture_states_number_dict)
         if (ASSOCIATED(ecapture_states_offset_dict)) call integer_dict_free(ecapture_states_offset_dict)

       end subroutine free_ecapture_info
      
      
      subroutine reaclib_init(jina_reaclib_filename)
         use const_def, only: mesa_data_dir
         character (len=*), intent(in) :: jina_reaclib_filename
         reaclib_dir = trim(mesa_data_dir) // '/rates_data'
         !reaclib_filename = 'jina_reaclib_results05301331'
         reaclib_filename = jina_reaclib_filename
         if (len_trim(reaclib_filename) == 0) &
            reaclib_filename = 'jina_reaclib_results_20171020_default'
      end subroutine reaclib_init
      
      
      subroutine allocate_reaclib_data(r, n, ierr)
         type(reaclib_data), intent(inout) :: r
         integer, intent(in) :: n
         integer, intent(out) :: ierr
         ierr = 0
         allocate( &
            r% chapter(n), r% species1(max_species_per_reaction*n), &
            r% label(n), r% reaction_flag(n), r% reverse_flag(n), &
            r% Qvalue(n), r% coefficients1(ncoefficients*n), stat=ierr)
         r% species(1:max_species_per_reaction,1:n) => r% species1(1:max_species_per_reaction*n)
         r% coefficients(1:ncoefficients,1:n) => r% coefficients1(1:ncoefficients*n)
      end subroutine allocate_reaclib_data
      

      subroutine free_reaclib_data(reaclib)
         type(reaclib_data), intent(inout) :: reaclib
         if (associated(reaclib% chapter)) & 
            deallocate( &
               reaclib% chapter, reaclib% species1, reaclib% label, reaclib% reaction_flag, &
               reaclib% reverse_flag, reaclib% Qvalue, reaclib% coefficients1)
      end subroutine free_reaclib_data
      

      subroutine allocate_reaction_data(r, n, nweak, ierr)
         type(reaction_data), intent(out) :: r
         integer, intent(in) :: n ! number of rates
         integer, intent(in) :: nweak ! number of weaklib rates
         integer, intent(out) :: ierr
         allocate( &
            r% reaction_handle(n), r% reverse_handle(n), r% category(n), r% chapter(n), &
            r% weaklib_ids(nweak), r% also_in_reaclib(nweak), &
            r% pspecies(max_species_per_reaction, n), r% reaction_flag(n), &
            r% weight(n), r% weight_reverse(n), r% coefficients(ncoefficients, n), &
            r% weak_mask(n), r% inverse_coefficients(ninverse_coeff, n), &
            r% inverse_exp(n), r% inverse_part(npart, n), r% Q(n), r% Qneu(n), stat=ierr)
         nullify(r% reaction_dict)
         nullify(r% reverse_dict)

      end subroutine allocate_reaction_data
      

      subroutine free_reaction_data(r)
         use utils_lib, only: integer_dict_free
         type(reaction_data), intent(inout) :: r
         if (associated(r% chapter)) then
            deallocate( &
               r% reaction_handle, r% reverse_handle, r% category, r% chapter, &
               r% weaklib_ids, r% also_in_reaclib, r% pspecies, &
               r% reaction_flag, r% weight, r% weight_reverse, r% coefficients, r% weak_mask, &
               r% inverse_coefficients, r% inverse_exp, r% inverse_part, r% Q, r% Qneu)
         end if
         if (associated(r% reaction_dict)) call integer_dict_free(r% reaction_dict)
         if (associated(r% reverse_dict)) call integer_dict_free(r% reverse_dict)
      end subroutine free_reaction_data


      subroutine do_start_rates_def_init(ierr)
         use math_lib
         integer, intent(out) :: ierr
         integer :: i         
         ierr = 0
         do i=1,num_one_thirds
            one_third_power(i) = pow(dble(i),1d0/3d0)
         end do
         do i=1,num_pow_186
            pow_186(i) = pow(dble(i),1.86d0)
         end do
         call set_rattab_range(5.30102999566398d0, 10.301029995664d0)
         
         reaclib_min_T9 = 1d-2 
            ! need <= 2d-3 for pre-ms li7 burning
            ! pre-ms deuterium burning needs much lower (4d-4)
            ! but that seems to cause problems during advanced burning.
                        
      end subroutine do_start_rates_def_init
      
      
      subroutine set_rattab_range(tlo, thi)
         use math_lib
         real(dp), intent(in) :: tlo, thi         
         if (abs(thi - tlo) < 1d-6) then
            rattab_tlo = tlo
            rattab_temp_lo = exp10(rattab_tlo)
            rattab_thi = rattab_tlo
            rattab_temp_hi = rattab_temp_lo
            nrattab = 1
            rattab_tstp = 0
            return
         end if
         rattab_thi = thi
         rattab_tlo = tlo
         rattab_temp_hi = exp10(rattab_thi)
         rattab_temp_lo = exp10(rattab_tlo)
         nrattab = rattab_points_per_decade*(rattab_thi - rattab_tlo) + 1
         if (nrattab <= 1) then
            rattab_thi = rattab_tlo
            nrattab = 1
            rattab_tstp = 0
         else
            rattab_tstp = (rattab_thi-rattab_tlo)/(nrattab-1)
         end if
      end subroutine set_rattab_range
      
      
      integer function get_rates_reaction_id(reaction_name) result(value)
         use utils_lib, only: integer_dict_lookup
         character (len=*), intent(in)  :: reaction_name 
         integer :: ierr
         integer :: indx
         ierr = 0
         call integer_dict_lookup(reaction_names_dict, reaction_name, value, ierr)
         if (ierr /= 0) value = 0
      end function get_rates_reaction_id


      subroutine create_ecapture_dict_key(ecapture_lhs, ecapture_rhs, key)
         character(len=iso_name_length), intent(in) :: ecapture_lhs, ecapture_rhs
         character(len=2*iso_name_length+1), intent(out) :: key
         key = trim(ecapture_lhs) // ' ' // trim(ecapture_rhs)
      end subroutine create_ecapture_dict_key
      
      
      subroutine create_weak_dict_key(weak_lhs, weak_rhs, key)
         character(len=iso_name_length), intent(in) :: weak_lhs, weak_rhs
         character(len=2*iso_name_length+1), intent(out) :: key
         key = trim(weak_lhs) // ' ' // trim(weak_rhs)
      end subroutine create_weak_dict_key
      

      integer function do_get_weak_rate_id(lhs, rhs) ! returns 0 if reaction not found
         use utils_lib
         character (len=*), intent(in)  :: lhs, rhs 
         integer :: ierr, i
         character(len=2*iso_name_length+1) :: key
         character (len=iso_name_length) :: lhs_name, rhs_name
         ierr = 0
         do_get_weak_rate_id = 0
         lhs_name = adjustl(lhs)
         rhs_name = adjustl(rhs)
         call create_weak_dict_key(lhs_name, rhs_name, key)
         call integer_dict_lookup(weak_reactions_dict, key, i, ierr)
         if (ierr /= 0) then
            !write(*,*) 'failed in integer_dict_lookup for key ' // trim(key)
            return
         end if
         do_get_weak_rate_id = i
      end function do_get_weak_rate_id
      

      integer function do_get_weak_info_list_id(lhs, rhs) ! returns 0 if reaction not found
         ! value can be used to index weak_info_list_halflife and weak_info_list_Qneu
         use utils_lib
         character (len=*), intent(in)  :: lhs, rhs ! names as in weak_info.list file
         integer :: ierr, i
         character(len=2*iso_name_length+1) :: key
         character (len=iso_name_length) :: lhs_name, rhs_name
         ierr = 0
         do_get_weak_info_list_id = 0
         lhs_name = adjustl(lhs)
         rhs_name = adjustl(rhs)
         call create_weak_dict_key(lhs_name, rhs_name, key)
         call integer_dict_lookup(weak_info_list_dict, key, i, ierr)
         if (ierr /= 0) then
            !write(*,'(a)') 'get_weak_info_list_id failed for ' // trim(key)
            return
         end if
         do_get_weak_info_list_id = i
      end function do_get_weak_info_list_id
      
      
      integer function get_num_reaction_inputs(ir)
         integer, intent(in) :: ir
         integer :: j
         include 'formats.dek'
         if (max_num_reaction_inputs == 3) then
            if (reaction_inputs(5,ir) /= 0) then
               get_num_reaction_inputs = 3
            else if (reaction_inputs(3,ir) /= 0) then
               get_num_reaction_inputs = 2
            else if (reaction_inputs(1,ir) /= 0) then
               get_num_reaction_inputs = 1
            else 
               get_num_reaction_inputs = 0
            end if
            return
         end if
         get_num_reaction_inputs = max_num_reaction_inputs
         do j = 1, 2*max_num_reaction_inputs-1, 2
            if (reaction_inputs(j,ir) == 0) then
               get_num_reaction_inputs = (j-1)/2
               exit
            end if
         end do
      end function get_num_reaction_inputs
      
      
      integer function get_num_reaction_outputs(ir)
         integer, intent(in) :: ir
         integer :: j
         if (max_num_reaction_outputs == 3) then
            if (reaction_outputs(5,ir) /= 0) then
               get_num_reaction_outputs = 3
            else if (reaction_outputs(3,ir) /= 0) then
               get_num_reaction_outputs = 2
            else if (reaction_outputs(1,ir) /= 0) then
               get_num_reaction_outputs = 1
            else 
               get_num_reaction_outputs = 0
            end if
            return
         end if
         get_num_reaction_outputs = max_num_reaction_outputs
         do j = 1, 2*max_num_reaction_outputs-1, 2
            if (reaction_outputs(j,ir) == 0) then
               get_num_reaction_outputs = (j-1)/2
               exit
            end if
         end do
      end function get_num_reaction_outputs


      end module rates_def

