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
         
         if (id == 1 .and. .not. s% RSP2_flag) then
            write(*,*) 'star id==1, but not RSP2_flag'
            stop 'extras_startup'
         end if
         
         if (id == 2 .and. s% RSP2_flag) then
            write(*,*) 'star id==2, but RSP2_flag'
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
         names(1) = 'r_R'
         names(2) = 'v_R'
         names(3) = 'Teff_R'
         names(4) = 'L_R'
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
         logical, parameter :: do_drel = .false.
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
         names(1) = 'v_R'
         if (do_drel) names(1) = 'v_drel'
         names(2) = 'Y_face_R'
         if (do_drel) names(2) = 'Y_drel'
         names(3) = 'w_R'
         if (do_drel) names(3) = 'w_drel'
         names(4) = 'Lc_div_L_R'
         if (do_drel) names(4) = 'Chi_drel'
         names(5) = 'logR_R'
         names(6) = 'logP_R'
         names(7) = 'logT_R'
         names(8) = 'logRho_R'
         names(9) = 'logL_R'
         names(10) = 'xCOUPL'
         names(11) = 'COUPL_R'
         if (do_drel) names(11) = 'CPL_drel' 
         names(12) = 'xSOURCE'
         names(13) = 'SRC_R'
         if (do_drel) names(13) = 'SRC_drel'
         names(14) = 'xDAMP'
         names(15) = 'DAMP_R'
         if (do_drel) names(15) = 'DAMP_drel'
         names(16) = 'xEq'
         names(17) = 'Eq_R'
         if (do_drel) names(17) = 'Eq_drel'
         names(18) = 'xUq'
         names(19) = 'Uq_R'
         if (do_drel) names(19) = 'Uq_drel'
         names(20) = 'xL'
         names(21) = 'L_r'
         if (do_drel) names(19) = 'L_drel'

         if (.not. associated(s_other% Y_face)) then
            vals(1:nz,:) = 0d0
         else if (s_other% nz /= nz) then
            vals(1:nz,:) = 0d0
         else
            do k=1,nz
               vals(k,1) = s_other% v(k)*1d-5
               if (do_drel) vals(k,1) = rel_diff(s_other% v(k), s% v(k)) ! v_drel
               
               vals(k,2) = s_other% Y_face(k)
               if (do_drel) vals(k,2) = rel_diff(s_other% Y_face(k), s% Y_face(k)) ! Y_drel
               
               if (s_other% RSP2_flag) then
                  vals(k,3) = s_other% w(k)
               else if (s_other% RSP_flag) then
                  vals(k,3) = s_other% RSP_w(k)
               else
                  vals(k,3) = 0d0
               end if
               if (do_drel) vals(k,3) = rel_diff(vals(k,3), s% w(k)) ! w_drel    
                         
               vals(k,4) = s_other% Lc(k)/s_other% L(k)
               if (do_drel) vals(k,4) = rel_diff(s_other% Chi(k), s% Chi(k)) ! Chi_drel
                                
               vals(k,5) = safe_log10(s_other% r(k)/Rsun)
               vals(k,6) = s_other% lnPeos(k)/ln10
               vals(k,7) = s_other% lnT(k)/ln10
               vals(k,8) = s_other% lnd(k)/ln10
               vals(k,9) = safe_log10(s_other% L(k)/Lsun)
               vals(k,10) = s% COUPL(k)
               vals(k,11) = s_other% COUPL(k)
               if (do_drel) vals(k,11) = rel_diff(s_other% COUPL(k), s% COUPL(k)) ! CPL_drel

               vals(k,12) = s% SOURCE(k)
               vals(k,13) = s_other% SOURCE(k)
               if (do_drel) vals(k,13) = rel_diff(s_other% SOURCE(k), s% SOURCE(k)) ! SRC_drel
               
               vals(k,14) = s% DAMP(k)
               vals(k,15) = s_other% DAMP(k)
               if (do_drel) vals(k,15) = rel_diff(s_other% DAMP(k), s% DAMP(k)) ! DAMP_drel
               
               vals(k,16) = s% Eq(k)
               vals(k,17) = s_other% Eq(k)
               if (do_drel) vals(k,17) = rel_diff(s_other% Eq(k), s% Eq(k)) ! Eq_drel
               
               vals(k,18) = s% Uq(k)
               vals(k,19) = s_other% Uq(k)
               if (do_drel) vals(k,19) = rel_diff(s_other% Uq(k), s% Uq(k)) ! Uq_drel
                          
               vals(k,20) = s% L(k)
               vals(k,21) = s_other% L(k)
               if (do_drel) vals(k,21) = rel_diff(s_other% L(k), s% L(k)) ! L_drel
               
            end do
         end if
         
         contains
               
         real(dp) function rel_diff(b, a, atol, rtol) result(d)
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
      
