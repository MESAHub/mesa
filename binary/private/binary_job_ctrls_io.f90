! ***********************************************************************
!
!   Copyright (C) 2013  Bill Paxton and Pablo Marchant
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

      module binary_job_ctrls_io

      use const_def, only: dp
      use binary_def

      implicit none

      include "binary_job_controls.inc"

      namelist /binary_job/ &
         show_binary_log_description_at_start, &
         binary_history_columns_file, &
         warn_binary_extra, &
         inlist_names, &
      ! extra files (Maybe overkill with so few inlist parameters)
         read_extra_binary_job_inlist1, extra_binary_job_inlist1_name, &
         read_extra_binary_job_inlist2, extra_binary_job_inlist2_name, &
         read_extra_binary_job_inlist3, extra_binary_job_inlist3_name, &
         read_extra_binary_job_inlist4, extra_binary_job_inlist4_name, &
         read_extra_binary_job_inlist5, extra_binary_job_inlist5_name, &
         evolve_both_stars, &
         relax_primary_to_th_eq, &
         log_Lnuc_div_L_for_relax_primary_to_th_eq, &
         min_age_for_relax_primary_to_th_eq, &
         max_steps_for_relax_primary_to_th_eq, &
         no_history_during_relax_primary_to_th_eq, &
         reset_age_for_relax_primary_to_th_eq, &
         tsync_for_relax_primary_to_th_eq, &
         change_ignore_rlof_flag, &
         change_initial_ignore_rlof_flag, &
         new_ignore_rlof_flag, &
         change_model_twins_flag, &
         change_initial_model_twins_flag, &
         new_model_twins_flag, &
         change_point_mass_i, &
         change_initial_point_mass_i, &
         new_point_mass_i, &
         change_m1, &
         change_initial_m1, &
         new_m1, &
         change_m2, &
         change_initial_m2, &
         new_m2, &
         change_separation_eccentricity, &
         change_initial_separation_eccentricity, &
         change_period_eccentricity, &
         change_initial_period_eccentricity, &
         new_separation, &
         new_period, &
         new_eccentricity

      contains


      subroutine do_read_binary_job(b, filename, ierr)
         use utils_lib
         type (binary_info), pointer :: b
         character(*), intent(in) :: filename
         integer, intent(out) :: ierr
         character (len=strlen) :: binary_job_namelist_name
         binary_job_namelist_name = ''
         ierr = 0
         call set_default_binary_job_controls
         call read_binary_job_file(b, filename, 1, ierr)
      end subroutine do_read_binary_job


      recursive subroutine read_binary_job_file(b, filename, level, ierr)
         use utils_lib
         character(*), intent(in) :: filename
         type (binary_info), pointer :: b
         integer, intent(in) :: level
         integer, intent(out) :: ierr
         logical :: read_extra1, read_extra2, read_extra3, read_extra4, read_extra5
         character (len=strlen) :: message, extra1, extra2, extra3, extra4, extra5
         integer :: unit

         ierr = 0

         if (level >= 10) then
            write(*,*) 'ERROR: too many levels of nested extra binary_job inlist files'
            ierr = -1
            return
         end if

         if (len_trim(filename) > 0) then
            open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
            if (ierr /= 0) then
               write(*, *) 'Failed to open control namelist file ', trim(filename)
               return
            end if
            read(unit, nml=binary_job, iostat=ierr)
            close(unit)
            if (ierr /= 0) then
               write(*, *)
               write(*, *)
               write(*, *)
               write(*, *)
               write(*, '(a)') &
                  'Failed while trying to read control namelist file: ' // trim(filename)
               write(*, '(a)') &
                  'Perhaps the following runtime error message will help you find the problem.'
               write(*, *)
               open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
               read(unit, nml=binary_job)
               close(unit)
               return
            end if
         end if

         call store_binary_job_controls(b, ierr)

         ! recursive calls to read other inlists

         read_extra1 = read_extra_binary_job_inlist1
         read_extra_binary_job_inlist1 = .false.
         extra1 = extra_binary_job_inlist1_name
         extra_binary_job_inlist1_name = 'undefined'

         read_extra2 = read_extra_binary_job_inlist2
         read_extra_binary_job_inlist2 = .false.
         extra2 = extra_binary_job_inlist2_name
         extra_binary_job_inlist2_name = 'undefined'

         read_extra3 = read_extra_binary_job_inlist3
         read_extra_binary_job_inlist3 = .false.
         extra3 = extra_binary_job_inlist3_name
         extra_binary_job_inlist3_name = 'undefined'

         read_extra4 = read_extra_binary_job_inlist4
         read_extra_binary_job_inlist4 = .false.
         extra4 = extra_binary_job_inlist4_name
         extra_binary_job_inlist4_name = 'undefined'

         read_extra5 = read_extra_binary_job_inlist5
         read_extra_binary_job_inlist5 = .false.
         extra5 = extra_binary_job_inlist5_name
         extra_binary_job_inlist5_name = 'undefined'

         if (read_extra1) then
            call read_binary_job_file(b, extra1, level+1, ierr)
            if (ierr /= 0) return
         end if

         if (read_extra2) then
            call read_binary_job_file(b, extra2, level+1, ierr)
            if (ierr /= 0) return
         end if

         if (read_extra3) then
            call read_binary_job_file(b, extra3, level+1, ierr)
            if (ierr /= 0) return
         end if

         if (read_extra4) then
            call read_binary_job_file(b, extra4, level+1, ierr)
            if (ierr /= 0) return
         end if

         if (read_extra5) then
            call read_binary_job_file(b, extra5, level+1, ierr)
            if (ierr /= 0) return
         end if

      end subroutine read_binary_job_file


      subroutine store_binary_job_controls(b, ierr)
         type (binary_info), pointer :: b
         integer, intent(out) :: ierr

         ierr = 0

         b% job% show_binary_log_description_at_start = show_binary_log_description_at_start
         b% job% binary_history_columns_file = binary_history_columns_file
         b% job% warn_binary_extra = warn_binary_extra
         b% job% inlist_names(:) = inlist_names(:)

         b% job% evolve_both_stars = evolve_both_stars
         b% job% relax_primary_to_th_eq = relax_primary_to_th_eq
         b% job% log_Lnuc_div_L_for_relax_primary_to_th_eq = log_Lnuc_div_L_for_relax_primary_to_th_eq
         b% job% min_age_for_relax_primary_to_th_eq = min_age_for_relax_primary_to_th_eq
         b% job% max_steps_for_relax_primary_to_th_eq = max_steps_for_relax_primary_to_th_eq
         b% job% no_history_during_relax_primary_to_th_eq = no_history_during_relax_primary_to_th_eq
         b% job% reset_age_for_relax_primary_to_th_eq = reset_age_for_relax_primary_to_th_eq
         b% job% tsync_for_relax_primary_to_th_eq = tsync_for_relax_primary_to_th_eq

         b% job% change_ignore_rlof_flag = change_ignore_rlof_flag
         b% job% change_initial_ignore_rlof_flag = change_initial_ignore_rlof_flag
         b% job% new_ignore_rlof_flag = new_ignore_rlof_flag
         b% job% change_model_twins_flag = change_model_twins_flag
         b% job% change_initial_model_twins_flag = change_initial_model_twins_flag
         b% job% new_model_twins_flag = new_model_twins_flag
         b% job% change_point_mass_i = change_point_mass_i
         b% job% change_initial_point_mass_i = change_initial_point_mass_i
         b% job% new_point_mass_i = new_point_mass_i
         b% job% change_m1 = change_m1
         b% job% change_initial_m1 = change_initial_m1
         b% job% new_m1 = new_m1
         b% job% change_m2 = change_m2
         b% job% change_initial_m2 = change_initial_m2
         b% job% new_m2 = new_m2
         b% job% change_separation_eccentricity = change_separation_eccentricity
         b% job% change_initial_separation_eccentricity = change_initial_separation_eccentricity
         b% job% change_period_eccentricity = change_period_eccentricity
         b% job% change_initial_period_eccentricity = change_initial_period_eccentricity
         b% job% new_separation = new_separation
         b% job% new_period = new_period
         b% job% new_eccentricity = new_eccentricity

      end subroutine store_binary_job_controls


      subroutine set_default_binary_job_controls
         include 'binary_job.defaults'
      end subroutine set_default_binary_job_controls


      subroutine set_binary_job_controls_for_writing(b, ierr)
         type (binary_info), pointer :: b
         integer, intent(out) :: ierr

         ierr = 0

         show_binary_log_description_at_start = b% job% show_binary_log_description_at_start
         binary_history_columns_file = b% job% binary_history_columns_file
         warn_binary_extra = b% job% warn_binary_extra
         inlist_names(:) = b% job% inlist_names(:)

         evolve_both_stars = b% job% evolve_both_stars
         evolve_both_stars = b% job% evolve_both_stars
         relax_primary_to_th_eq = b% job% relax_primary_to_th_eq
         log_Lnuc_div_L_for_relax_primary_to_th_eq = b% job% log_Lnuc_div_L_for_relax_primary_to_th_eq
         min_age_for_relax_primary_to_th_eq = b% job% min_age_for_relax_primary_to_th_eq
         max_steps_for_relax_primary_to_th_eq = b% job% max_steps_for_relax_primary_to_th_eq
         no_history_during_relax_primary_to_th_eq = b% job% no_history_during_relax_primary_to_th_eq
         reset_age_for_relax_primary_to_th_eq = b% job% reset_age_for_relax_primary_to_th_eq
         tsync_for_relax_primary_to_th_eq = b% job% tsync_for_relax_primary_to_th_eq

         change_ignore_rlof_flag = b% job% change_ignore_rlof_flag
         change_initial_ignore_rlof_flag = b% job% change_initial_ignore_rlof_flag
         new_ignore_rlof_flag = b% job% new_ignore_rlof_flag
         change_model_twins_flag = b% job% change_model_twins_flag
         change_initial_model_twins_flag = b% job% change_initial_model_twins_flag
         new_model_twins_flag = b% job% new_model_twins_flag
         change_point_mass_i = b% job% change_point_mass_i
         change_initial_point_mass_i = b% job% change_initial_point_mass_i
         new_point_mass_i = b% job% new_point_mass_i
         change_m1 = b% job% change_m1
         change_initial_m1 = b% job% change_initial_m1
         new_m1 = b% job% new_m1
         change_m2 = b% job% change_m2
         change_initial_m2 = b% job% change_initial_m2
         new_m2 = b% job% new_m2
         change_separation_eccentricity = b% job% change_separation_eccentricity
         change_initial_separation_eccentricity = b% job% change_initial_separation_eccentricity
         change_period_eccentricity = b% job% change_period_eccentricity
         change_initial_period_eccentricity = b% job% change_initial_period_eccentricity
         new_separation = b% job% new_separation
         new_period = b% job% new_period
         new_eccentricity = b% job% new_eccentricity

      end subroutine set_binary_job_controls_for_writing


      subroutine do_write_binary_job(b, filename, ierr)
         type (binary_info), pointer :: b
         character(*), intent(in) :: filename
         integer, intent(out) :: ierr
         integer :: io
         ierr = 0
         call set_binary_job_controls_for_writing(b, ierr)
         if (ierr /= 0) return
         open(newunit=io, file=trim(filename), action='write', status='replace', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(filename)
            return
         end if
         write(io, nml=binary_job, iostat=ierr)
         close(io)
      end subroutine do_write_binary_job


      end module binary_job_ctrls_io

