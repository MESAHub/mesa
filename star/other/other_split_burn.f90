! ***********************************************************************
!
!   Copyright (C) 2022  The MESA Team
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
 
      module other_split_burn

      ! consult star/other/README for general usage instructions
      ! control name: use_other_split_burn = .true.
      ! procedure pointer: s% other_split_burn => my_routine
      ! This also require op_split_burn = .true. as well as setting the other op_split_burn options

      use star_def

      implicit none
      
            
      contains
      
      
      subroutine null_other_split_burn(  &
         id, k, net_handle, eos_handle, num_isos, num_reactions, t_start, t_end, starting_x, &
         num_times_for_interpolation, times, log10Ts_f1, log10Rhos_f1, etas_f1, &
         dxdt_source_term, rate_factors, &
         weak_rate_factor, reaction_Qs, reaction_neuQs, &
         screening_mode, &
         stptry, max_steps, eps, odescal, &
         use_pivoting, trace, dbg,  burner_finish_substep, &
         burn_lwork, burn_work_array, &
         net_lwork, net_work_array, &
         ! results
         ending_x, eps_nuc_categories, avg_eps_nuc, eps_neu_total, &
         nfcn, njac, nstep, naccpt, nrejct, ierr)
         use net_def
         use chem_def, only: num_categories
         
         integer, intent(in) :: net_handle, eos_handle, id
         integer, intent(in) :: num_isos
         integer, intent(in) :: num_reactions
         real(dp), intent(in) :: t_start, t_end, starting_x(:) ! (num_isos)
         
         integer, intent(in) :: num_times_for_interpolation 
            ! ending time is times(num_times); starting time is 0
         real(dp), pointer, intent(in) :: times(:) ! (num_times) 
         real(dp), pointer, intent(in) :: log10Ts_f1(:) ! =(4,numtimes) interpolant for log10T(time)
         real(dp), pointer, intent(in) :: log10Rhos_f1(:) ! =(4,numtimes) interpolant for log10Rho(time)
         real(dp), pointer, intent(in) :: etas_f1(:) ! =(4,numtimes) interpolant for eta(time)
         real(dp), pointer, intent(in) :: dxdt_source_term(:) ! (num_isos)  or null if no source term.
         real(dp), intent(in), pointer :: rate_factors(:) ! (num_reactions)
         real(dp), intent(in) :: weak_rate_factor
         real(dp), pointer, intent(in) :: reaction_Qs(:) ! (rates_reaction_id_max)
         real(dp), pointer, intent(in) :: reaction_neuQs(:) ! (rates_reaction_id_max)
         integer, intent(in) :: screening_mode ! see screen_def
         real(dp), intent(in) :: stptry ! try this for 1st step.  0 means try in 1 step.
         integer, intent(in) :: max_steps ! maximal number of allowed steps.
         real(dp), intent(in) :: eps, odescal ! tolerances.  e.g., set both to 1d-6
         logical, intent(in) :: use_pivoting ! for matrix solves
         logical, intent(in) :: trace, dbg
         interface
            include 'burner_finish_substep.inc'
         end interface
         integer, intent(in) :: net_lwork, burn_lwork
         real(dp), intent(inout), pointer :: burn_work_array(:) ! (burn_lwork)
         real(dp), intent(inout), pointer :: net_work_array(:) ! (net_lwork)

         ! These should be set for the output
         real(dp), intent(inout) :: ending_x(:) ! (num_isos)
         real(dp), intent(inout) :: eps_nuc_categories(:) ! (num_categories)
         real(dp), intent(out) :: avg_eps_nuc, eps_neu_total
         integer, intent(out) :: nfcn    ! number of function evaluations
         integer, intent(out) :: njac    ! number of jacobian evaluations
         integer, intent(out) :: nstep   ! number of computed steps
         integer, intent(out) :: naccpt  ! number of accepted steps
         integer, intent(out) :: nrejct  ! number of rejected steps
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         integer,intent(in) :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         contains


         subroutine my_burn_finish_substep(nstp, time, y, ierr)
            use chem_def, only: category_name
            integer,intent(in) :: nstp
            real(dp), intent(in) :: time, y(:)
            integer, intent(out) :: ierr
            real(dp) :: frac, step_time
            integer :: j, i
            include 'formats'
            ierr = 0
            ! This routine does nothing other than set ierr = 0,
            ! Leave this empty if you dont want to do anything otherwise
            ! Else replace the burn_finish_substep with this routine if
            ! passing to net_1_zone_burn
         end subroutine my_burn_finish_substep


      end subroutine null_other_split_burn


      end module other_split_burn
      
      
      
      
