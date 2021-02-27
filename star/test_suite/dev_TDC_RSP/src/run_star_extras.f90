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
      
      contains

      include 'test_suite_extras.inc'
      include 'multi_stars_extras.inc'
      
      
      subroutine report(s)
         type (star_info), pointer :: s
         character (len=128) :: filename
         integer :: id, ierr, k, iounit
         include 'formats'
         return
         if (s% model_number < 1) return
         id = s% id
         if (id == 1) then
            filename = 'l1'
         else
            filename = 'l2'
         end if
         if (s% TDC_flag) then
            write(*,*) 'write TDC results to ' // trim(filename), id
         else if (s% RSP_flag) then
            write(*,*) 'write RSP results to ' // trim(filename), id
         else
            return
         end if
         open(newunit=iounit, file=trim(filename), action='write', status='replace', iostat=ierr)
         if (ierr /= 0) stop 'failed in open'
         if (s% TDC_flag) then 
            write(iounit,*) 'TDC'
         else if (s% RSP_flag) then
            write(iounit,*) 'RSP'
         end if
         do k=1,s% nz
            write(iounit,2) 'm', k, s% m(k)
            write(iounit,2) 'dm', k, s% dm(k)
            write(iounit,2) 'dm_bar', k, s% dm_bar(k)
            write(iounit,2) 'logR', k, s% lnR(k)/ln10
            write(iounit,2) 'logT', k, s% lnT(k)/ln10
            write(iounit,2) 'logRho', k, s% lnd(k)/ln10
            write(iounit,2) 'v', k, s% v(k)
            write(iounit,2) 'L', k, s% L(k)
            write(iounit,2) 'Lr', k, s% Lr(k)
            write(iounit,2) 'Lc', k, s% Lc(k)
            write(iounit,2) 'Lt', k, s% Lt(k)
            write(iounit,2) 'COUPL', k, s% COUPL(k)
            write(iounit,2) 'SOURCE', k, s% SOURCE(k)
            write(iounit,2) 'DAMP', k, s% DAMP(k)
            write(iounit,2) 'DAMPR', k, s% DAMPR(k)
            write(iounit,2) 'Eq', k, s% Eq(k)
            write(iounit,2) 'Uq', k, s% Uq(k)
            write(iounit,2) 'Hp_face', k, s% Hp_face(k)
            write(iounit,2) 'Chi', k, s% Chi(k)
            write(iounit,2) 'Y_face', k, s% Y_face(k)
            write(iounit,2) 'PII', k, s% PII(k)
            if (s% TDC_flag) then 
               write(iounit,2) 'w', k, s% w(k)
               write(iounit,2) 'eturb', k, s% w(k)**2
            else if (s% RSP_flag) then
               write(iounit,2) 'w', k, s% RSP_w(k)
               write(iounit,2) 'eturb', k, s% RSP_Et(k)
            end if
         end do
         write(iounit,1) 'L_center', s% L_center
         write(iounit,1) 'm_center', s% m_center
         write(iounit,1) 'R_center', s% R_center
         write(iounit,1) 'v_center', s% v_center
         if (s% TDC_flag) then 
            write(iounit,*) 'TDC'
         else if (s% RSP_flag) then
            write(iounit,*) 'RSP'
         end if
         close(iounit)
         write(*,*)
         
         if (id == 2) stop 'terminate test after report for star2'
      end subroutine report


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
         real(dp) :: target_period, rel_run_E_err
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call report(s)
         if (s% x_integer_ctrl(1) <= 0) return
         if (s% rsp_num_periods < s% x_integer_ctrl(1)) return
         write(*,*)
         write(*,*)
         write(*,*)
         target_period = s% x_ctrl(1)
         rel_run_E_err = s% cumulative_energy_error/s% total_energy
         write(*,*) 'rel_run_E_err', rel_run_E_err
         if (s% total_energy /= 0d0 .and. abs(rel_run_E_err) > 1d-5) then
            write(*,*) '*** BAD rel_run_E_error ***', &
            s% cumulative_energy_error/s% total_energy
         else if (abs(s% rsp_period/(24*3600) - target_period) > 1d-2) then
            write(*,*) '*** BAD ***', s% rsp_period/(24*3600) - target_period, &
               s% rsp_period/(24*3600), target_period
         else
            write(*,*) 'good match for period', &
               s% rsp_period/(24*3600), target_period
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
         
         if (id == 1 .and. .not. s% TDC_flag) then
            write(*,*) 'star id==1, but not TDC_flag'
            stop 'extras_startup'
         end if
         
         if (id == 2 .and. .not. s% RSP_flag) then
            write(*,*) 'star id==2, but not RSP_flag'
            stop 'extras_startup'
         end if
         
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
         how_many_extra_profile_columns = 4
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
         names(1) = 'v_alt'
         names(2) = 'Y_face_alt'
         names(3) = 'w_alt'
         names(4) = 'Lr_div_L_alt'
         if (.not. associated(s_other% Y_face)) then
            vals(1:nz,:) = 0d0
         else if (s_other% nz /= nz) then
            vals(1:nz,:) = 0d0
         else
            do k=1,nz
               vals(k,1) = s_other% v(k)
               vals(k,2) = s_other% Y_face(k)
               if (s_other% TDC_flag) then
                  vals(k,3) = s_other% w(k)
               else if (s_other% RSP_flag) then
                  vals(k,3) = s_other% RSP_w(k)
               else
                  vals(k,3) = 0d0
               end if
               vals(k,4) = s_other% Lr(k)/s_other% L(k)
            end do
         end if
      end subroutine data_for_extra_profile_columns

!   Hp_face
!   Y_face
!   PII_face
!   Chi
!   COUPL
!   SOURCE
!   DAMP
!   DAMPR
!   Eq
!   Uq
!   Lr_div_L
!   Lc_div_L
!   Lr
!   Lc
!   w


      end module run_star_extras
      
