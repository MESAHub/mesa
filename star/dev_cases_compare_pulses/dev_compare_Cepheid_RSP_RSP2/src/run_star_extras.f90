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
      integer :: id_other, k
      include 'formats'
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_finish_step = keep_going
      if (s% id == 2 .and. s% RSP2_flag .and. s% x_logical_ctrl(1)) then
         !write(*,2) 'set star2 xh = star1 xh to resynchronize', s% model_number
         id_other = 1
         call star_ptr(id_other, s_other, ierr)
         if (ierr /= 0) return
         ! sanity check before start copying
         if (s% model_number /= s_other% model_number .or. &
            s% nz /= s_other% nz .or. &
            s% R_center /= s_other% R_center .or. &
            s% L_center /= s_other% L_center .or. &
            s% m_center /= s_other% m_center .or. &
            s% xmstar /= s_other% xmstar) then
            write(*, 2) 'something does not match', s% model_number
            call mesa_error(__FILE__, __LINE__, 'extras_finish_step')
         end if
         do k = 1, s% nz
            s% xh(s% i_lum, k) = s_other% L(k)
            s% xh(s% i_v, k) = s_other% v(k)
            s% xh(s% i_w, k) = s_other% RSP_w(k)
            s% xh(s% i_lnR, k) = s_other% lnR(k)
            s% xh(s% i_lnT, k) = s_other% lnT(k)
            s% xh(s% i_lnd, k) = s_other% lnd(k)
         end do
      end if
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
         write(*, *) 'star id==1, but not RSP_flag'
         call mesa_error(__FILE__, __LINE__, 'extras_startup')
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
      how_many_extra_history_columns = 6
   end function how_many_extra_history_columns
   
   
   subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len = maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s, s_other
      integer :: id_other
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
      names(1) = 'r_R'
      names(2) = 'v_R'
      names(3) = 'Teff_R'
      names(4) = 'L_R'
      names(5) = 'log_tot_KE_R'
      names(6) = 'num_periods'
      vals(1) = s_other% r(1) / Rsun
      vals(2) = s_other% v(1) / 1d5 ! kms
      vals(3) = s_other% Teff
      vals(4) = s_other% L(1) / Lsun
      vals(5) = safe_log10(s_other% total_radial_kinetic_energy_end)
      vals(6) = s_other% rsp_num_periods
   end subroutine data_for_extra_history_columns
   
   
   integer function how_many_extra_profile_columns(id)
      use star_def, only : star_info
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_profile_columns = 38
   end function how_many_extra_profile_columns
   
   
   subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
      use star_def, only : star_info, maxlen_profile_column_name
      use const_def, only : dp
      integer, intent(in) :: id, n, nz
      character (len = maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(nz, n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s, s_other
      real(dp) :: val
      integer :: i, k, id_other
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
      
      i = 1
      names(i) = 'v_R'; i = i + 1
      names(i) = 'v_diff'; i = i + 1
      names(i) = 'v_drel'; i = i + 1
      
      names(i) = 'Y_face_R'; i = i + 1
      names(i) = 'Y_diff'; i = i + 1
      names(i) = 'Y_drel'; i = i + 1
      
      names(i) = 'w_R'; i = i + 1
      names(i) = 'w_diff'; i = i + 1
      names(i) = 'w_drel'; i = i + 1
      
      names(i) = 'Lc_div_L_R'; i = i + 1
      names(i) = 'Lc_diff'; i = i + 1
      names(i) = 'Lc_drel'; i = i + 1
      
      names(i) = 'COUPL_R'; i = i + 1
      names(i) = 'CPL_diff'; i = i + 1
      names(i) = 'CPL_drel'; i = i + 1
      
      names(i) = 'SRC_R'; i = i + 1
      names(i) = 'SRC_diff'; i = i + 1
      names(i) = 'SRC_drel'; i = i + 1
      
      names(i) = 'DAMP_R'; i = i + 1
      names(i) = 'DAMP_diff'; i = i + 1
      names(i) = 'DAMP_drel'; i = i + 1
      
      names(i) = 'DAMPR_R'; i = i + 1
      names(i) = 'DAMPR_diff'; i = i + 1
      names(i) = 'DAMPR_drel'; i = i + 1
      
      names(i) = 'Eq_R'; i = i + 1
      names(i) = 'Eq_diff'; i = i + 1
      names(i) = 'Eq_drel'; i = i + 1
      
      names(i) = 'Uq_R'; i = i + 1
      names(i) = 'Uq_diff'; i = i + 1
      names(i) = 'Uq_drel'; i = i + 1
      
      names(i) = 'Pvsc_R'; i = i + 1
      names(i) = 'Pvsc_diff'; i = i + 1
      names(i) = 'Pvsc_drel'; i = i + 1
      
      names(i) = 'logR_R'; i = i + 1
      names(i) = 'logP_R'; i = i + 1
      names(i) = 'logT_R'; i = i + 1
      names(i) = 'logRho_R'; i = i + 1
      names(i) = 'logL_R'; i = i + 1
      
      if (.not. associated(s_other% Y_face)) then
         vals(1:nz, :) = 0d0
      else if (s_other% nz /= nz) then
         vals(1:nz, :) = 0d0
      else
         do k = 1, nz
            
            i = 1
            vals(k, i) = s_other% v(k) * 1d-5; i = i + 1
            vals(k, i) = s_other% v(k) - s% v(k); i = i + 1
            vals(k, i) = rel_diff(s_other% v(k), s% v(k)); i = i + 1
            
            vals(k, i) = s_other% Y_face(k); i = i + 1
            vals(k, i) = s_other% Y_face(k) - s% Y_face(k); i = i + 1
            vals(k, i) = rel_diff(s_other% Y_face(k), s% Y_face(k)); i = i + 1
            
            if (s_other% RSP2_flag) then
               val = s_other% w(k)
            else if (s_other% RSP_flag) then
               val = s_other% RSP_w(k)
            else
               val = 0d0
            end if
            vals(k, i) = val; i = i + 1
            vals(k, i) = val - s% w(k); i = i + 1
            if (.false. .and. abs(vals(k, i - 1)) > 1d0 .and. .not. s% doing_first_model_of_run) then
               write(*, 3) 'w diff', k, s% model_number, vals(k, i - 1), val, s% w(k)
               call mesa_error(__FILE__, __LINE__, 'data_for_extra_profile_columns')
            end if
            vals(k, i) = rel_diff(val, s% w(k)); i = i + 1
            
            val = s_other% Lc(k) / s_other% L(k)
            vals(k, i) = val; i = i + 1
            vals(k, i) = val - s% Lc(k) / s% L(k); i = i + 1
            vals(k, i) = rel_diff(val, s% Lc(k) / s% L(k)); i = i + 1
            
            vals(k, i) = s_other% COUPL(k); i = i + 1
            vals(k, i) = s_other% COUPL(k) - s% COUPL(k); i = i + 1
            vals(k, i) = rel_diff(s_other% COUPL(k), s% COUPL(k)); i = i + 1
            
            vals(k, i) = s_other% SOURCE(k); i = i + 1
            vals(k, i) = s_other% SOURCE(k) - s% SOURCE(k); i = i + 1
            vals(k, i) = rel_diff(s_other% SOURCE(k), s% SOURCE(k)); i = i + 1
            
            vals(k, i) = s_other% DAMP(k); i = i + 1
            vals(k, i) = s_other% DAMP(k) - s% DAMP(k); i = i + 1
            vals(k, i) = rel_diff(s_other% DAMP(k), s% DAMP(k)); i = i + 1
            
            vals(k, i) = s_other% DAMPR(k); i = i + 1
            vals(k, i) = s_other% DAMPR(k) - s% DAMPR(k); i = i + 1
            vals(k, i) = rel_diff(s_other% DAMPR(k), s% DAMPR(k)); i = i + 1
            
            vals(k, i) = s_other% Eq(k); i = i + 1
            vals(k, i) = s_other% Eq(k) - s% Eq(k); i = i + 1
            vals(k, i) = rel_diff(s_other% Eq(k), s% Eq(k)); i = i + 1
            
            vals(k, i) = s_other% Uq(k); i = i + 1
            vals(k, i) = s_other% Uq(k) - s% Uq(k); i = i + 1
            vals(k, i) = rel_diff(s_other% Uq(k), s% Uq(k)); i = i + 1
            
            vals(k, i) = s_other% Pvsc(k); i = i + 1
            vals(k, i) = s_other% Pvsc(k) - s% Pvsc(k); i = i + 1
            vals(k, i) = rel_diff(s_other% Pvsc(k), s% Pvsc(k)); i = i + 1
            
            vals(k, i) = safe_log10(s_other% r(k) / Rsun); i = i + 1
            vals(k, i) = s_other% lnPeos(k) / ln10; i = i + 1
            vals(k, i) = s_other% lnT(k) / ln10; i = i + 1
            vals(k, i) = s_other% lnd(k) / ln10; i = i + 1
            vals(k, i) = safe_log10(s_other% L(k) / Lsun); i = i + 1
         
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
            atl = 1d-9
         end if
         if (present(rtol)) then
            rtl = rtol
         else
            rtl = 1d0
         end if
         d = (a - b) / (atl + rtl * max(abs(a), abs(b)))
      end function rel_diff
      
      real(dp) function fix_if_bad(v)
         use utils_lib, only : is_bad
         real(dp), intent(in) :: v
         if (is_bad(v)) then
            fix_if_bad = 100
         else
            fix_if_bad = v
         end if
      end function fix_if_bad
   
   end subroutine data_for_extra_profile_columns

end module run_star_extras
      
