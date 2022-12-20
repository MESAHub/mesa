! ***********************************************************************
!
!   Copyright (C) 2013  Pablo Marchant
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
         read_extra_binary_job_inlist, extra_binary_job_inlist_name, &
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
         new_eccentricity, &
         pgbinary_flag

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
         logical, dimension(max_extra_inlists) :: read_extra
         character (len=strlen) :: message
         character (len=strlen), dimension(max_extra_inlists) :: extra
         integer :: unit, i

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
         do i=1, max_extra_inlists
            read_extra(i) = read_extra_binary_job_inlist(i)
            read_extra_binary_job_inlist(i) = .false.
            extra(i) = extra_binary_job_inlist_name(i)
            extra_binary_job_inlist_name(i) = 'undefined'
            
            if (read_extra(i)) then
               call read_binary_job_file(b, extra(i), level+1, ierr)
               if (ierr /= 0) return
            end if
         end do
         
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
         b% job% pgbinary_flag = pgbinary_flag

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
         pgbinary_flag = b% job% pgbinary_flag

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


      subroutine get_binary_job(b, name, val, ierr)
         use utils_lib, only: StrUpCase
         type (binary_info), pointer :: b
         character(len=*),intent(in) :: name
         character(len=*), intent(out) :: val
         integer, intent(out) :: ierr
   
         character(len(name)) :: upper_name
         character(len=512) :: str
         integer :: iounit,iostat,ind,i
   
   
         ! First save current controls
         call set_binary_job_controls_for_writing(b, ierr)
         if(ierr/=0) return
   
         ! Write namelist to temporay file
         open(newunit=iounit,status='scratch')
         write(iounit,nml=binary_job)
         rewind(iounit)
   
         ! Namelists get written in captials
         upper_name = StrUpCase(name)
         val = ''
         ! Search for name inside namelist
         do 
            read(iounit,'(A)',iostat=iostat) str
            ind = index(str,trim(upper_name))
            if( ind /= 0 ) then
               val = str(ind+len_trim(upper_name)+1:len_trim(str)-1) ! Remove final comma and starting =
               do i=1,len(val)
                  if(val(i:i)=='"') val(i:i) = ' '
               end do
               exit
            end if
            if(is_iostat_end(iostat)) exit
         end do   
   
         if(len_trim(val) == 0 .and. ind==0 ) ierr = -1
   
         close(iounit)
   
      end subroutine get_binary_job
   
      subroutine set_binary_job(b, name, val, ierr)
         type (binary_info), pointer :: b
         character(len=*), intent(in) :: name, val
         character(len=len(name)+len(val)+14) :: tmp
         integer, intent(out) :: ierr
   
         ! First save current controls
         call set_binary_job_controls_for_writing(b, ierr)
         if(ierr/=0) return
   
         tmp=''
         tmp = '&binary_job '//trim(name)//'='//trim(val)//' /'
   
         ! Load into namelist
         read(tmp, nml=binary_job)
   
         ! Add to star
         call store_binary_job_controls(b, ierr)
         if(ierr/=0) return
   
      end subroutine set_binary_job


      end module binary_job_ctrls_io

