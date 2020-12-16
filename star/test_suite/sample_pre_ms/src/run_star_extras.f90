! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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


      ! create zams basically just does a series of create pre-MS,
      ! run until it reaches zams, and save it to a file for future use.
      ! that can be efficient if you are going to do many runs
      ! using the same composition but different masses.

      ! however, making a zams file is a big job and probably only
      ! worth the effort if you are going to run a very large number
      ! of cases.   For most applications it is much easier just to
      ! create a single model using create_pre_main_sequence_model.


 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use run_star_support

      
      implicit none

      
      include "test_suite_extras_def.inc"
      
      
      contains

      include "test_suite_extras.inc"


      subroutine do_run
         
         integer :: id, ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         
         write(*,*) 'do create pre ms samples'

         ierr = 0
         call test_suite_startup(s, .false., ierr)

         call do_read_star_job('inlist', ierr)
         if (failed('do_read_star_job')) return

         id = id_from_read_star_job
         id_from_read_star_job = 0
         
         call star_ptr(id, s, ierr)
         if (failed('star_ptr')) return
      
         call starlib_init(s, ierr)
         if (failed('star_init')) return

         s% inlist_fname = 'inlist'
         
         call star_set_kap_and_eos_handles(id, ierr)
         if (failed('set_star_kap_and_eos_handles')) return
                  
         call star_setup(id, 'inlist', ierr)
         if (failed('star_setup')) return

         call extras_controls(s% id, ierr)
         if (failed('extras_controls')) return

         call do_star_job_controls_before(id, s, .false., ierr)
         if (failed('do_star_job_controls_before')) return
         
         call do_sample_pre_ms(s, ierr)
         
         call test_suite_after_evolve(s, ierr)

         
         contains
         

         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = (ierr /= 0)
            if (failed) then
               write(*, *) trim(str) // ' ierr', ierr
            end if
         end function failed


      end subroutine do_run    


      
      subroutine do_sample_pre_ms(s, ierr)
         use mtx_lib, only: lapack_decsol
         use num_def, only: square_matrix_type
         use num_lib
         type (star_info), pointer :: s
         integer, intent(out) :: ierr
         
         integer :: i, j, k, n, id, result, result_reason
         logical :: okay
         integer, parameter :: num_Zs = 3, num_Ms = 3
         real(dp) :: Zs(num_Zs), Ms(num_Ms)

         include 'formats'
         
         Zs = (/ 0.03d0, 0.02d0, 0.00d0 /)
         Ms = (/ 0.1d0, 2d0, 20d0 /)
         
         ierr = 0
         id = s% id

         okay = .true.

         ierr = 0         
         do i=1, num_Ms
            do j=1, num_Zs

               s% initial_z = Zs(j)
               s% initial_mass = Ms(i)
                     
               s% mesh_delta_coeff = 0.5d0
               if (s% initial_mass > 1) s% mesh_delta_coeff = 0.8d0
               if (s% initial_mass > 80) s% mesh_delta_coeff = 1d0
            
               write(*,*) 
               write(*,1) 's% initial_z', s% initial_z
               write(*,1) 's% initial_mass', s% initial_mass
               write(*,*) 
            
               call star_create_pre_ms_model( &
                  id, s% job% pre_ms_T_c, s% job% pre_ms_guess_rho_c, s% job% pre_ms_d_log10_P, &
                  s% job% pre_ms_logT_surf_limit, s% job% pre_ms_logP_surf_limit, &
                  s% job% initial_zfracs, &
                  s% job% dump_missing_metals_into_heaviest, &
                  .false., '', 0, ierr)
               if (failed('star_create_pre_ms_model')) exit
            
               call pre_ms_evolve(s, id, ierr)
               if (failed('pre_ms_evolve')) exit 
               
               do k = 1, 6
                  write(*,*)
               end do
               
            end do
         end do
         
         call free_star(id, ierr)
         if (failed('free_star')) return 
         
         call starlib_shutdown
         
         write(*, *)
         if (okay) then
            write(*, '(a)') 'finished sample pre-ms'
         else
            write(*, '(a)') 'failed during sample pre-ms'
         end if
         write(*, *)
         
         contains
         
         
         logical function failed(str)
            character (len=*), intent(in) :: str
            failed = (ierr /= 0)
            if (failed) then
               write(*, *)
               write(*, *) trim(str) // ' ierr', ierr
               okay = .false.
               stop 1
            end if
         end function failed

      end subroutine do_sample_pre_ms
      
      
      subroutine pre_ms_evolve(s, id, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         integer, parameter :: lipar=0, lrpar=0
         logical, parameter :: restore_at_end = .true.
         integer, target :: ipar_ary(lipar)
         real(dp), target :: rpar_ary(lrpar)
         integer, pointer :: ipar(:)
         real(dp), pointer :: rpar(:)
         ipar => ipar_ary
         rpar => rpar_ary
         call star_evolve_to_check_point( &
            id, before_pre_ms_evolve, pre_ms_evolve_adjust_model, pre_ms_evolve_check_model, &
            pre_ms_evolve_finish_step, restore_at_end, &
            lipar, ipar, lrpar, rpar, ierr)
      end subroutine pre_ms_evolve


      subroutine before_pre_ms_evolve(s, id, lipar, ipar, lrpar, rpar, ierr)
         use star_def, only:star_info
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         ierr = 0
      end subroutine before_pre_ms_evolve
      
      integer function pre_ms_evolve_adjust_model(s, id, lipar, ipar, lrpar, rpar)
         use star_def, only:star_info
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         pre_ms_evolve_adjust_model = keep_going
      end function pre_ms_evolve_adjust_model
      
      integer function pre_ms_evolve_check_model(s, id, lipar, ipar, lrpar, rpar)
         use star_def, only:star_info
         type (star_info), pointer :: s
         integer, intent(in) :: id, lipar, lrpar
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         pre_ms_evolve_check_model = bare_bones_check_model(id) 
         if (pre_ms_evolve_check_model /= keep_going) return
         if (s% model_number >= 100) then
            pre_ms_evolve_check_model = terminate
            s% termination_code = t_extras_check_model
         end if
      end function pre_ms_evolve_check_model
      
      
      integer function pre_ms_evolve_finish_step(s)
         type (star_info), pointer :: s
         pre_ms_evolve_finish_step = keep_going
      end function pre_ms_evolve_finish_step

      
      include 'standard_run_star_extras.inc'

      end module run_star_extras
      
