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
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use run_star_support
      
      implicit none
      
      include 'test_suite_extras_def.inc'
      include 'multi_stars_extras_def.inc'

      integer :: RSP2_num_periods
      real(dp) :: RSP2_period, time_started
            
      contains

      include 'test_suite_extras.inc'
      include 'multi_stars_extras.inc'


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         include 'formats'
         ierr = 0
         extras_start_step = keep_going
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end function extras_start_step
      
      
      ! returns either keep_going or terminate.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         real(dp) :: target_period, rel_run_E_err, time_ended
         type (star_info), pointer :: s, s_other
         integer :: id_other
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         if (.not. s% RSP2_flag) return
         if (s% x_integer_ctrl(1) <= 0) return
         ! check_cycle_completed when v(1) goes from positive to negative
         if (s% v(1)*s% v_start(1) > 0d0 .or. s% v(1) > 0d0) return
         ! at max radius
         ! either start of 1st cycle, or end of current
         if (time_started == 0) then
            time_started = s% time
            write(*,*) 'RSP2 first maximum radius, period calculations start at model, day', &
               s% model_number, s% time/(24*3600)
            return
         end if
         RSP2_num_periods = RSP2_num_periods + 1
         time_ended = s% time
         !if (abs(s% v(1)-s% v_start(1)).gt.1.0d-10) & ! tweak the end time
         !   time_ended = time_started + (s% time - time_started)*s% v_start(1)/(s% v_start(1) - s% v(1))
         RSP2_period = time_ended - time_started
         write(*,*) 'RSP2 period', RSP2_num_periods, RSP2_period/(24*3600)
         time_started = time_ended
         if (RSP2_num_periods < s% x_integer_ctrl(1)) return
         write(*,*)
         write(*,*)
         write(*,*)
         target_period = s% x_ctrl(1)
         rel_run_E_err = s% cumulative_energy_error/s% total_energy
         write(*,*) 'RSP2 rel_run_E_err', rel_run_E_err
         if (s% total_energy /= 0d0 .and. abs(rel_run_E_err) > 1d-5) then
            write(*,*) '*** RSP2 BAD rel_run_E_error ***', &
            s% cumulative_energy_error/s% total_energy
         else if (abs(RSP2_period/(24*3600) - target_period) > 1d-2) then
            write(*,*) '*** RSP2 BAD period ***', RSP2_period/(24*3600) - target_period, &
               RSP2_period/(24*3600), target_period
         else
            write(*,*) 'RSP2 good match for period', &
               RSP2_period/(24*3600), target_period
         end if
         write(*,*)
         write(*,*)
         write(*,*)
         extras_finish_step = terminate
      end function extras_finish_step
      
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_startup(s, restart, ierr)
         
         if (id == 1 .and. .not. s% RSP_flag) then
            write(*,*) 'star id==1, but not RSP_flag'
            stop 'extras_startup'
         end if
         
         if (id == 2 .and. .not. s% RSP2_flag) then
            write(*,*) 'star id==2, but not RSP2_flag'
            stop 'extras_startup'
         end if
         
         RSP2_num_periods = 0
         RSP2_period = 0
         time_started = 0

      end subroutine extras_startup
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: dt
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         call test_suite_after_evolve(s, ierr)
      end subroutine extras_after_evolve


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 4
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s, s_other
         integer :: id_other
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (n /= 4) then
            ierr = -1
            return
         end if
         if (id == 1) then
            id_other = 2
         else if (id == 2) then
            id_other = 1
         else
            ierr = -1
            return
         end if
         call star_ptr(id_other, s_other, ierr)
         if (ierr /= 0) return
         names(1) = 'r_Rsp'
         names(2) = 'v_Rsp'
         names(3) = 'Teff_Rsp'
         names(4) = 'L_Rsp'
         vals(1) = s_other% r(1)/Rsun
         vals(2) = s_other% v(1)/1d5 ! kms
         vals(3) = s_other% Teff
         vals(4) = s_other% L(1)/Lsun
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 21
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s, s_other
         integer :: k, id_other
         include 'formats'
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (id == 1) then
            id_other = 2
         else if (id == 2) then
            id_other = 1
         else
            ierr = -1
            return
         end if
         call star_ptr(id_other, s_other, ierr)
         if (ierr /= 0) return
         
         names(1) = 'v_drel' ! 'v_R'
         names(2) = 'Y_drel' ! 'Y_face_R'
         names(3) = 'w_drel' ! 'w_R'
         names(4) = 'Chi_drel' ! 'Lr_div_L_R'
         names(5) = 'logR_R'
         names(6) = 'logP_R'
         names(7) = 'logT_R'
         names(8) = 'logRho_R'
         names(9) = 'logL_R'
         names(10) = 'xCOUPL'
         names(11) = 'CPL_drel' ! 'COUPL_R'
         names(12) = 'xSOURCE'
         names(13) = 'SRC_drel' ! 'SRC_R'
         names(14) = 'xDAMP'
         names(15) = 'DAMP_drel' ! 'DAMP_R'
         names(16) = 'xEq'
         names(17) = 'Eq_drel' ! 'Eq_R'
         names(18) = 'xUq'
         names(19) = 'Uq_drel' ! 'Uq_R'
         names(20) = 'Uq_err'
         names(21) = 'Uq_diff'

         if (.not. associated(s_other% Y_face)) then
            vals(1:nz,:) = 0d0
         else if (s_other% nz /= nz) then
            vals(1:nz,:) = 0d0
         else
            do k=1,nz
               !vals(k,1) = s_other% v(k)*1d-5
               vals(k,1) = rel_diff(s_other% v(k), s% v(k)) ! v_drel
               
               !vals(k,2) = s_other% Y_face(k)
               vals(k,2) = rel_diff(s_other% Y_face(k), s% Y_face(k)) ! Y_drel
               
               if (s_other% RSP2_flag) then
                  vals(k,3) = s_other% w(k)
               else if (s_other% RSP_flag) then
                  vals(k,3) = s_other% RSP_w(k)
               else
                  vals(k,3) = 0d0
               end if
               !if (k == 194) then
               !   write(*,2) 'w RSP_w w_start RSP_w_start', k, &
               !      s% w(k), s_other% RSP_w(k), s% w_start(k), s_other% RSP_w_start(k)
               !end if
               vals(k,3) = rel_diff(vals(k,3), s% w(k)) ! w_drel    
                         
               !vals(k,4) = s_other% Lr(k)/s_other% L(k)
               vals(k,4) = rel_diff(s_other% Chi(k), s% Chi(k)) ! Chi_drel
               !if (k==109) then
               !   write(*,2) 'Chi_rsp Chi_rsp2 reldiff', k, &
               !      s_other% Chi(k), s% Chi(k), vals(k,4)
               !end if
                                
               vals(k,5) = safe_log10(s_other% r(k)/Rsun)
               vals(k,6) = s_other% lnPeos(k)/ln10
               vals(k,7) = s_other% lnT(k)/ln10
               vals(k,8) = s_other% lnd(k)/ln10
               vals(k,9) = safe_log10(s_other% L(k)/Lsun)
               vals(k,10) = s% COUPL(k)
               !vals(k,11) = s_other% COUPL(k)
               vals(k,11) = rel_diff(s_other% COUPL(k), s% COUPL(k)) ! CPL_drel

               vals(k,12) = s% SOURCE(k)
               !vals(k,13) = s_other% SOURCE(k)
               vals(k,13) = rel_diff(s_other% SOURCE(k), s% SOURCE(k)) ! SRC_drel
               
               vals(k,14) = s% DAMP(k)
               !vals(k,15) = s_other% DAMP(k)
               vals(k,15) = rel_diff(s_other% DAMP(k), s% DAMP(k)) ! DAMP_drel
               
               vals(k,16) = s% Eq(k)
               !vals(k,17) = s_other% Eq(k)
               vals(k,17) = rel_diff(s_other% Eq(k), s% Eq(k)) ! Eq_drel
               
               vals(k,18) = s% Uq(k)
               !vals(k,19) = s_other% Uq(k)
               vals(k,19) = rel_diff(s_other% Uq(k), s% Uq(k)) ! Uq_drel
               
                          
               vals(k,20) = s% Pvsc(k)
               vals(k,21) = s_other% Pvsc(k)
               vals(k,20) = (s% Uq(k) - s_other% Uq(k)) / &
                  (1d-6 + 1d-3*max(abs(s% Uq(k)),abs(s_other% Uq(k))))
               vals(k,21) = s% Uq(k) - s_other% Uq(k)
               
               if (.false.) then ! debugging xtra values
                  vals(k,10) = s% xtra1_array(k)
                  vals(k,11) = s_other% xtra1_array(k)
                  vals(k,12) = s% xtra2_array(k)
                  vals(k,13) = s_other% xtra2_array(k)
                  vals(k,14) = s% xtra3_array(k)
                  vals(k,15) = s_other% xtra3_array(k)
                  vals(k,16) = s% xtra4_array(k)
                  vals(k,17) = s_other% xtra4_array(k)
                  vals(k,18) = s% xtra5_array(k)
                  vals(k,19) = s_other% xtra5_array(k)
                  vals(k,20) = s% xtra6_array(k)
                  vals(k,21) = s_other% xtra6_array(k)
               end if
               
            end do
         end if
         
         contains
               
            real(dp) function rel_diff(a, b, atol, rtol) result(d)
               real(dp), intent(in) :: a, b
               real(dp), intent(in), optional :: atol, rtol
               real(dp) :: atl, rtl
               if (present(atol)) then
                  atl = atol
               else
                  atl = 1d-6
               end if
               if (present(rtol)) then
                  rtl = rtol
               else
                  rtl = 1d-3
               end if
               d = (a - b)/(atl + rtl*max(abs(a),abs(b)))
            end function rel_diff
         
            real(dp) function fix_if_bad(v)
               use utils_lib, only: is_bad
               real(dp), intent(in) :: v
               if (is_bad(v)) then
                  fix_if_bad = 100
               else
                  fix_if_bad = v
               end if
            end function fix_if_bad
            
      end subroutine data_for_extra_profile_columns

      end module run_star_extras
      
