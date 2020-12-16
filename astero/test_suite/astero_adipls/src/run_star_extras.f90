! ***********************************************************************
!
!   Copyright (C) 2011  Bill Paxton
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
      
      implicit none
      
      include "test_suite_extras_def.inc"

      
      ! these routines are called by the standard run_star check_model
      contains

      include "test_suite_extras.inc"
      
      
      subroutine extras_controls(id, ierr)
         use astero_def, only: star_astero_procs
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
         include 'set_star_astero_procs.inc'
      end subroutine extras_controls

      
      subroutine set_my_vars(id, ierr) ! called from star_astero code
         !use astero_search_data, only: include_my_var1_in_chi2, my_var1
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ! my_var's are predefined in the simplex_search_data.
         ! this routine's job is to assign those variables to current value in the model.
         ! it is called whenever a new value of chi2 is calculated.
         ! only necessary to set the my_var's you are actually using.
         ierr = 0
         !if (include_my_var1_in_chi2) then
            call star_ptr(id, s, ierr)
            if (ierr /= 0) return
            !my_var1 = s% Teff
         !end if
      end subroutine set_my_vars
      
      
      subroutine will_set_my_param(id, i, new_value, ierr) ! called from star_astero code
         !use astero_search_data, only: vary_my_param1
         integer, intent(in) :: id
         integer, intent(in) :: i ! which of my_param's will be set
         real(dp), intent(in) :: new_value
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0

         ! old value has not yet been changed.
         ! do whatever is necessary for this new value.
         ! i.e. change whatever mesa params you need to adjust.
         ! as example, my_param1 is alpha_mlt
         ! if (i == 1) then
         !    call star_ptr(id, s, ierr)
         !    if (ierr /= 0) return
         !    s% mixing_length_alpha = new_value
         ! end if
      end subroutine will_set_my_param


      subroutine my_other_adipls_mode_info( &
            l, order, freq, inertia, x, y, aa, data, nn, iy, iaa, ispcpr, ierr)
         integer, intent(in) :: l, order
         real(dp), intent(in) :: freq, inertia
         real(dp), intent(in) :: x(1:nn), y(1:iy,1:nn), aa(1:iaa,1:nn), data(8)
         integer, intent(in) :: nn, iy, iaa, ispcpr
         integer, intent(out) :: ierr
         ierr = 0
         write(*,*) 'astero called my_other_adipls_mode_info'
      end subroutine my_other_adipls_mode_info
      
      
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
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         
         real(dp) :: dt, expected_freq, freq
         logical :: okay, store_for_adipls, save_mode_info
         integer :: l_to_match, order_to_match, order_to_save
         character (len=256) :: save_mode_filename
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         if (s% x_ctrl(1) > 0d0) then
         
            store_for_adipls = .true.
            l_to_match = 0
            order_to_match = 4
            expected_freq = s% x_ctrl(1)

            save_mode_info = .true.
            order_to_save = 5
            save_mode_filename = 'eigen.data'
         
            call get_adipls_frequency_info( &
               s, store_for_adipls, l_to_match, order_to_match, expected_freq, &
               save_mode_info, order_to_save, save_mode_filename, freq, okay, ierr)
            if (ierr /= 0) stop 1
            if (okay) then
               write(*,'(a,2f20.2)') 'got ok match for expected frequency', freq, expected_freq
            else
               write(*,'(a,2f20.2)') 'ERROR: bad match for expected frequency', freq, expected_freq
            end if
         
         end if

         call test_suite_after_evolve(s, ierr)
         
      end subroutine extras_after_evolve
      

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id

         type (star_info), pointer :: s
         real(dp) :: dt, expected_freq, freq
         logical :: okay, store_for_adipls, save_mode_info
         integer :: ierr, l_to_match, order_to_match, order_to_save
         character (len=256) :: save_mode_filename

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         extras_check_model = keep_going
         
         if (s% x_ctrl(1) > 0d0) then
         
            ! get frequencies for certain models
            if (mod(s% model_number,50) /= 0) return
         
            store_for_adipls = .true.
            l_to_match = 0
            order_to_match = 5
            expected_freq = -1
            save_mode_info = .false.
            order_to_save = 0
            save_mode_filename = ''
         
            call get_adipls_frequency_info( &
               s, store_for_adipls, l_to_match, order_to_match, expected_freq, &
               save_mode_info, order_to_save, save_mode_filename, freq, okay, ierr)
            if (ierr /= 0) extras_check_model = terminate
         
         end if
         
      end function extras_check_model
      
      
      subroutine set_my_param(s, i, new_value)
         type (star_info), pointer :: s
         integer, intent(in) :: i ! which of my_param's will be set
         real(dp), intent(in) :: new_value
         include 'formats'
         ! old value has not yet been changed.
         ! do whatever is necessary for this new value.
         ! i.e. change whatever mesa params you need to adjust.
         ! for example, my_param1 is mass
         if (i == 1) then
            s% job% new_mass = new_value
         end if
         
      end subroutine set_my_param


      subroutine get_adipls_frequency_info( &
            s, store_for_adipls, l_to_match, order_to_match, expected_freq, &
            save_mode_info, order_to_save, save_mode_filename, freq, okay, ierr)
         use astero_lib, only: adipls_get_one_el_info
         type (star_info), pointer :: s
         integer, intent(in) :: l_to_match, order_to_match, order_to_save
         logical, intent(in) :: store_for_adipls, save_mode_info
         character (len=*), intent(in) :: save_mode_filename
         real(dp), intent(in) :: expected_freq ! ignore if < 0
         real(dp), intent(out) :: freq
         logical, intent(out) :: okay ! true if expected_freq is okay
         integer, intent(out) :: ierr
         
         integer :: l, iscan, i, num
         real(dp) :: nu1, nu2, R, G, M
         real(dp), pointer, dimension(:) :: l_freq, l_inertia
         integer, pointer, dimension(:) :: l_order, l_em
         logical :: add_center_point, keep_surface_point, &
            add_atmosphere, do_restribute_mesh

         include 'formats'

         ierr = 0
         okay = .true.
         R = Rsun*s% photosphere_r
         G = s% cgrav(1)
         M = Msun*s% star_mass
         add_center_point = .true.
         keep_surface_point = .false.
         add_atmosphere = .true.
         do_restribute_mesh = .false.
         l = l_to_match
         
         nullify(l_freq)
         nullify(l_inertia)
         nullify(l_order)
         nullify(l_em)

         nu1 = 50
         nu2 = 1000
         iscan = 200
         
         !write(*,*) 'call adipls_get_one_el_info'
         call adipls_get_one_el_info( &
            s, l, nu1, nu2, iscan, R, G, M, &
            add_center_point, keep_surface_point, add_atmosphere, &
            do_restribute_mesh, store_for_adipls, &
            save_mode_info, order_to_save, save_mode_filename, &
            num, l_freq, l_inertia, l_order, l_em, ierr)
         !write(*,*) 'done adipls_get_one_el_info'
         if (ierr /= 0) then
            write(*,*) 'failed in adipls_get_one_el_info'
            stop 1
         end if
         write(*,*)
         write(*,'(2a8,99a20)') 'el', 'order', 'freq (microHz)', 'inertia'
         if (expected_freq > 0) okay = .false.
         do i = 1, num
            write(*,'(2i8,f20.10,e20.10,i20)') l, l_order(i), l_freq(i), l_inertia(i)
            if (expected_freq > 0 .and. l_order(i) == order_to_match) then
               freq = l_freq(i)
               okay = (abs(freq - expected_freq) < expected_freq*3d-2)
            end if
         end do
         deallocate(l_freq, l_inertia, l_order, l_em)
      
      end subroutine get_adipls_frequency_info


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
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
      end function extras_finish_step
      
      

      end module run_star_extras
      
