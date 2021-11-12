! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use utils_lib
      
      implicit none

      include "test_suite_extras_def.inc"

      real(dp) :: target_Teff
      real(dp) :: target_L
      
      ! these routines are called by the standard run_star check_model
      contains

      include "test_suite_extras.inc"
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: iounit
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
         
         ! Uncomment these lines if you wish to use the functions in this file,
         ! otherwise we use a null_ version which does nothing.
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         ! Once you have set the function pointers you want,
         ! then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          s% job% warn_run_star_extras =.false.       

         if (s% x_integer_ctrl(1) == 2) then
            !load initial mass and target values
            iounit = alloc_iounit(ierr)
            if (ierr /= 0) then
               return
            end if
            open(unit=iounit, file='initial_mass.dat', action='read', iostat=ierr)
            if (ierr == 0) then
               read(iounit,*) s% initial_mass
            else
               return
            end if
            close(iounit)
            call free_iounit(iounit)

            iounit = alloc_iounit(ierr)
            if (ierr /= 0) then
               termination_code_str(t_xtra1) = "Error allocating iounit for target_vals.dat"
               return
            end if
            open(unit=iounit, file='target_vals.dat', action='read', iostat=ierr)
            if (ierr == 0) then
               read(iounit,*) target_Teff
               read(iounit,*) target_L
            else
               return
            end if
            close(iounit)
            call free_iounit(iounit)

            include 'formats'

            write(*,1) "Loading model with mass:", s% initial_mass
            write(*,1) "Target value for Teff and L:", target_Teff, target_L
         end if
            
      end subroutine extras_controls
      
      ! None of the following functions are called unless you set their
      ! function point in extras_control.
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
      end subroutine extras_startup
      

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         

         if (s% x_integer_ctrl(1) == 2) then
            s% termination_code = t_xtra1
            if(s% Teff > target_Teff*1.005d0 .or. s% Teff < target_Teff*0.995d0) then
               extras_check_model = terminate
            else if (s% L(1) > target_L*1.005d0 .or. s% L(1) < target_L*0.995d0) then
               extras_check_model = terminate
            end if
            if (extras_check_model == terminate) then
               termination_code_str(t_xtra1) = "Relaxed model does not match input"
            else
               extras_check_model = terminate
               termination_code_str(t_xtra1) = "Values for Teff and L are within tolerance"
            end if
            return
         end if

         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         !note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.
         

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         !note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) call mesa_error(__FILE__,__LINE__,'data_for_extra_profile_columns')
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr, k, iounit
         type (star_info), pointer :: s
         
 99      format(99(1pd26.18))
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         if (s% x_integer_ctrl(1) == 1) then
            if (s% Teff < 8000) then
               write(*,*) "Reached target Teff, generate input files for relax"
               extras_finish_step = terminate
               s% termination_code = t_xtra1

               iounit = alloc_iounit(ierr)
               if (ierr /= 0) then
                  termination_code_str(t_xtra1) = "Error allocating iounit for entropy.dat"
                  return
               end if
               open(unit=iounit, file='entropy.dat', action='write', status='replace', iostat=ierr)
               if (ierr == 0) then
                  write(iounit,*) s% nz
                  do k=1, s%nz
                     write(iounit, 99) 1-(s% q(k)-s% dq(k)/2), s% entropy(k)*avo*kerg
                  end do
               else
                  termination_code_str(t_xtra1) = "Error writing entropy.dat"
                  return
               end if
               close(iounit)
               call free_iounit(iounit)

               iounit = alloc_iounit(ierr)
               if (ierr /= 0) then
                  termination_code_str(t_xtra1) = "Error allocating iounit for angular_momentum.dat"
                  return
               end if
               open(unit=iounit, file='angular_momentum.dat', action='write', status='replace', iostat=ierr)
               if (ierr == 0) then
                  write(iounit,*) s% nz
                  do k=1, s%nz
                     write(iounit, 99) 1-s% q(k), s% j_rot(k)
                  end do
               else
                  termination_code_str(t_xtra1) = "Error writing angular_momentum.dat"
                  return
               end if
               close(iounit)
               call free_iounit(iounit)

               iounit = alloc_iounit(ierr)
               if (ierr /= 0) then
                  termination_code_str(t_xtra1) = "Error allocating iounit for composition.dat"
                  return
               end if
               open(unit=iounit, file='composition.dat', action='write', status='replace', iostat=ierr)
               if (ierr == 0) then
                  write(iounit,*) s% nz, s% species
                  do k=1, s%nz
                     write(iounit, 99) 1-(s% q(k)-s% dq(k)/2), s% xa(:,k)
                  end do
               else
                  termination_code_str(t_xtra1) = "Error writing composition.dat"
                  return
               end if
               close(iounit)
               call free_iounit(iounit)

               iounit = alloc_iounit(ierr)
               if (ierr /= 0) then
                  termination_code_str(t_xtra1) = "Error allocating iounit for initial_mass.dat"
                  return
               end if
               open(unit=iounit, file='initial_mass.dat', action='write', status='replace', iostat=ierr)
               if (ierr == 0) then
                  write(iounit, 99) s% star_mass
               else
                  termination_code_str(t_xtra1) = "Error writing initial_mass.dat"
                  return
               end if
               close(iounit)
               call free_iounit(iounit)

               iounit = alloc_iounit(ierr)
               if (ierr /= 0) then
                  termination_code_str(t_xtra1) = "Error allocating iounit for target_vals.dat"
                  return
               end if
               open(unit=iounit, file='target_vals.dat', action='write', status='replace', iostat=ierr)
               if (ierr == 0) then
                  write(iounit, 99) s% Teff
                  write(iounit, 99) s% L(1)
               else
                  termination_code_str(t_xtra1) = "Error writing target_vals.dat"
                  return
               end if
               close(iounit)
               call free_iounit(iounit)

               termination_code_str(t_xtra1) = "produced files for relax"
               return
            end if
         end if

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         real(dp) :: dt
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call test_suite_after_evolve(s, ierr)

      end subroutine extras_after_evolve
      
      

      end module run_star_extras
      
