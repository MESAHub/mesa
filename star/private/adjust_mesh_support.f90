! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
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

module adjust_mesh_support
   
   use star_private_def
   use const_def
   
   implicit none
   
   private
   public :: get_gval_info, check_validity
   
   logical, parameter :: dbg = .false.

contains
   
   
   subroutine get_gval_info(&
      s, delta_gval_max, gvals1, nz, &
      num_gvals, gval_names, &
      gval_is_xa_function, gval_is_logT_function, ierr)
      use chem_def
      use num_lib, only : find0
      use rates_def
      use alloc
      use mesh_functions, only : set_mesh_function_data, max_allowed_gvals
      type (star_info), pointer :: s
      real(dp), dimension(:), pointer :: delta_gval_max
      real(dp), dimension(:), pointer :: gvals1
      integer, intent(in) :: nz, num_gvals
      integer, intent(out) :: ierr
      character (len = 32), intent(out) :: gval_names(max_allowed_gvals)
      logical, intent(out), dimension(max_allowed_gvals) :: &
         gval_is_xa_function, gval_is_logT_function
      
      integer :: j, k, other_ierr
      logical, parameter :: dbg = .false.
      real(dp), allocatable, dimension(:) :: src
      real(dp) :: eps_min_for_delta, &
         dlog_eps_dlogP_full_off, dlog_eps_dlogP_full_on, alfa_czb
      real(dp), dimension(:, :), pointer :: gvals
      
      gvals(1:nz, 1:num_gvals) => gvals1(1:nz * num_gvals)
      
      include 'formats'
      
      ierr = 0
      
      eps_min_for_delta = exp10(s% mesh_dlog_eps_min_for_extra)
      
      dlog_eps_dlogP_full_off = s% mesh_dlog_eps_dlogP_full_off
      dlog_eps_dlogP_full_on = s% mesh_dlog_eps_dlogP_full_on
      
      call set_mesh_function_data(&
         s, num_gvals, gval_names, &
         gval_is_xa_function, gval_is_logT_function, gvals1, ierr)
      if (ierr /= 0) return
      
      allocate(src(nz))
      
      call set_delta_gval_max(src, ierr)
      if (ierr /= 0) return
      
      call smooth_gvals(nz, src, num_gvals, gvals)
   
   
   contains
      
      
      subroutine set_delta_gval_max(src, ierr)
         real(dp) :: src(:)
         integer, intent(out) :: ierr
         real(dp) :: P_exp, beta
         
         logical, parameter :: dbg = .false.
         
         include 'formats'
         
         ierr = 0
         delta_gval_max(1:nz) = 1d0
         
         if (s% mesh_Pgas_div_P_exponent /= 0) then
            P_exp = s% mesh_Pgas_div_P_exponent
            do k = 1, nz
               beta = s% Pgas(k) / s% Peos(k)
               delta_gval_max(k) = delta_gval_max(k) * pow(beta, P_exp)
            end do
         end if
         
         if (s% use_other_mesh_delta_coeff_factor) then
            do k = 1, nz
               s% mesh_delta_coeff_factor(k) = delta_gval_max(k)
            end do
            call s% other_mesh_delta_coeff_factor(s% id, ierr)
            if (ierr /= 0) then
               write(*, *) 'other_mesh_delta_coeff_factor returned ierr', ierr
               return
            end if
            do k = 1, nz
               delta_gval_max(k) = s% mesh_delta_coeff_factor(k)
            end do
         end if
         
         do j = 1, num_mesh_logX
            if (s% mesh_dlogX_dlogP_extra(j) > 0 .and. s% mesh_dlogX_dlogP_extra(j) < 1) &
               call do_mesh_dlogX_dlogP_coef(s, j)
         end do
         
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_pp_dlogP_extra, ipp)
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_cno_dlogP_extra, icno)
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_3alf_dlogP_extra, i3alf)
         
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_burn_c_dlogP_extra, i_burn_c)
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_burn_n_dlogP_extra, i_burn_n)
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_burn_o_dlogP_extra, i_burn_o)
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_burn_ne_dlogP_extra, i_burn_ne)
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_burn_na_dlogP_extra, i_burn_na)
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_burn_mg_dlogP_extra, i_burn_mg)
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_burn_si_dlogP_extra, i_burn_si)
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_burn_s_dlogP_extra, i_burn_s)
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_burn_ar_dlogP_extra, i_burn_ar)
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_burn_ca_dlogP_extra, i_burn_ca)
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_burn_ti_dlogP_extra, i_burn_ti)
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_burn_cr_dlogP_extra, i_burn_cr)
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_burn_fe_dlogP_extra, i_burn_fe)
         
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_cc_dlogP_extra, icc)
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_co_dlogP_extra, ico)
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_oo_dlogP_extra, ioo)
         
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_pnhe4_dlogP_extra, ipnhe4)
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_photo_dlogP_extra, iphoto)
         call do1_dlog_eps_dlogP_coef(s% mesh_dlog_other_dlogP_extra, iother)
         
         if (s% mesh_delta_coeff_factor_smooth_iters > 0) then ! smooth delta_gval_max
            
            do k = 1, nz
               src(k) = delta_gval_max(k)
            end do
            
            do j = 1, s% mesh_delta_coeff_factor_smooth_iters
               delta_gval_max(1) = (2 * src(1) + src(2)) / 3
               do k = 2, nz - 1
                  delta_gval_max(k) = (src(k - 1) + src(k) + src(k + 1)) / 3
               end do
               delta_gval_max(nz) = (2 * src(nz) + src(nz - 1)) / 3
               if (j == 3) exit
               src(1) = (2 * delta_gval_max(1) + delta_gval_max(2)) / 3
               do k = 2, nz - 1
                  src(k) = &
                     (delta_gval_max(k - 1) + delta_gval_max(k) + delta_gval_max(k + 1)) / 3
               end do
               src(nz) = (2 * delta_gval_max(nz) + delta_gval_max(nz - 1)) / 3
            end do
         
         end if
      
      end subroutine set_delta_gval_max
      
      
      subroutine do_mesh_dlogX_dlogP_coef(s, which)
         use chem_lib, only : chem_get_iso_id
         type (star_info), pointer :: s
         integer, intent(in) :: which
         real(dp) :: &
            logX_min_for_extra, dlogX_dlogP_extra, dlogX_dlogP_full_on, dlogX_dlogP_full_off, &
            X_min_for_extra, dlogX, dlogP, dlogX_dlogP, coef
         integer :: k, cid, j
         include 'formats'
         if (len_trim(s% mesh_logX_species(which)) == 0) return
         cid = chem_get_iso_id(s% mesh_logX_species(which))
         if (cid <= 0) return
         j = s% net_iso(cid)
         if (j == 0) return
         logX_min_for_extra = s% mesh_logX_min_for_extra(which)
         dlogX_dlogP_extra = s% mesh_dlogX_dlogP_extra(which)
         dlogX_dlogP_full_on = s% mesh_dlogX_dlogP_full_on(which)
         dlogX_dlogP_full_off = s% mesh_dlogX_dlogP_full_off(which)
         X_min_for_extra = exp10(max(-50d0, logX_min_for_extra))
         do k = 2, nz
            if (s% xa(j, k) < X_min_for_extra .or. s% xa(j, k - 1) < X_min_for_extra) cycle
            dlogX = abs(log10(s% xa(j, k) / s% xa(j, k - 1)))
            dlogP = max(1d-10, abs(s% lnPeos(k) - s% lnPeos(k - 1)) * ln10)
            dlogX_dlogP = dlogX / dlogP
            if (dlogX_dlogP <= dlogX_dlogP_full_off) cycle
            if (dlogX_dlogP >= dlogX_dlogP_full_on) then
               coef = dlogX_dlogP_extra
            else
               coef = 1 - (1 - dlogX_dlogP_extra) * &
                  (dlogX_dlogP - dlogX_dlogP_full_off) / (dlogX_dlogP_full_on - dlogX_dlogP_full_off)
            end if
            if (coef < delta_gval_max(k)) delta_gval_max(k) = coef
            if (coef < delta_gval_max(k - 1)) delta_gval_max(k - 1) = coef
            if (k < nz) then
               if (coef < delta_gval_max(k + 1)) delta_gval_max(k + 1) = coef
            end if
         end do
      end subroutine do_mesh_dlogX_dlogP_coef
      
      
      subroutine do1_dlog_eps_dlogP_coef(dlog_eps_dlogP_extra, cat)
         real(dp), intent(in) :: dlog_eps_dlogP_extra
         integer, intent(in) :: cat
         integer :: k
         real(dp) :: eps, epsm1, dlog_eps, dlogP, dlog_eps_dlogP, &
            extra, new_max, maxv
         include 'formats'
         
         if (dlog_eps_dlogP_extra <= 0 .or. dlog_eps_dlogP_extra >=1) return
         maxv = maxval(s% eps_nuc_categories(cat, 1:nz))
         if (maxv < eps_min_for_delta) return
         
         do k = 2, nz
            
            eps = s% eps_nuc_categories(cat, k)
            if (eps < eps_min_for_delta) cycle
            
            epsm1 = s% eps_nuc_categories(cat, k - 1)
            if (epsm1 < eps_min_for_delta) cycle
            
            maxv = maxval(s% eps_nuc_categories(:, k))
            if (maxv /= eps) cycle
            
            dlogP = (s% lnPeos(k) - s% lnPeos(k - 1)) / ln10
            if (abs(dlogP) < 1d-50) cycle
            
            dlog_eps = abs(log10(eps / epsm1))
            
            if (is_bad(dlog_eps)) then
               write(*, 3) 'adjust mesh support ' // trim(category_name(cat)), &
                  cat, k, s% eps_nuc_categories(cat, k)
               stop
            end if
            
            dlog_eps_dlogP = dlog_eps / dlogP
            if (dlog_eps_dlogP <= dlog_eps_dlogP_full_off) cycle
            
            if (dlog_eps_dlogP >= dlog_eps_dlogP_full_on) then
               extra = dlog_eps_dlogP_extra
            else
               extra = 1 - (1 - dlog_eps_dlogP_extra) * &
                  (dlog_eps_dlogP - dlog_eps_dlogP_full_off) / &
                  (dlog_eps_dlogP_full_on - dlog_eps_dlogP_full_off)
            end if
            new_max = extra
            if (new_max < delta_gval_max(k)) then
               delta_gval_max(k) = new_max
            end if
            if (new_max < delta_gval_max(k - 1)) then
               delta_gval_max(k - 1) = new_max
            end if
            if (k < nz) then
               if (new_max < delta_gval_max(k + 1)) delta_gval_max(k + 1) = new_max
            end if
         
         end do
      end subroutine do1_dlog_eps_dlogP_coef
      
      
      function blend_coef(coef_start, alfa) result(coef)
         real(dp), intent(in) :: coef_start, alfa
         real(dp) :: coef
         real(dp) :: beta
         
         ! this implements the following piecewise blend
         ! 0 < alfa < 0.5 : constant (coeff_start)
         ! 0.5 < alfa < 1 : linear (coeff_start -> 1)
         
         beta = 2d0 * (alfa - 0.5d0)
         if (beta .lt. 0d0) then
            coef = coef_start
         else
            coef = coef_start * (1d0 - beta) + beta
         endif
      
      end function blend_coef
   
   
   end subroutine get_gval_info
   
   
   subroutine single_peak(nz, src)
      integer, intent(in) :: nz
      real(dp), pointer :: src(:)
      integer :: k, kmax
      real(dp) :: prev, val
      kmax = maxloc(src(1:nz), dim = 1)
      val = src(kmax)
      do k = kmax + 1, nz
         prev = val
         val = src(k)
         if (val > prev) val = prev
         src(k) = val
      end do
      val = src(kmax)
      do k = kmax - 1, 1, -1
         prev = val
         val = src(k)
         if (val > prev) val = prev
         src(k) = val
      end do
   end subroutine single_peak
   
   
   subroutine increasing_inward(nz, src)
      integer, intent(in) :: nz
      real(dp), pointer :: src(:)
      integer :: k
      real(dp) :: prev, val
      val = src(1)
      do k = 2, nz
         prev = val
         val = src(k)
         if (val < prev) val = prev
         src(k) = val
      end do
   end subroutine increasing_inward
   
   
   subroutine decreasing_outward(nz, src)
      integer, intent(in) :: nz
      real(dp), pointer :: src(:)
      integer :: k
      real(dp) :: prev, val
      val = src(nz)
      do k = nz - 1, 1, -1
         prev = val
         val = src(k)
         if (val > prev) val = prev
         src(k) = val
      end do
   end subroutine decreasing_outward
   
   
   subroutine increasing_outward(nz, src)
      integer, intent(in) :: nz
      real(dp), pointer :: src(:)
      integer :: k
      real(dp) :: prev, val
      val = src(nz)
      do k = nz - 1, 1, -1
         prev = val
         val = src(k)
         if (val < prev) val = prev
         src(k) = val
      end do
   end subroutine increasing_outward
   
   
   subroutine smooth_gvals(nz, src, num_gvals, gvals)
      integer, intent(in) :: nz, num_gvals
      real(dp) :: src(:), gvals(:, :)
      integer :: k, i
      
      do i = 1, num_gvals
         
         do k = 1, nz
            src(k) = gvals(k, i)
         end do
         
         gvals(1, i) = (2 * src(1) + src(2)) / 3
         do k = 2, nz - 1
            gvals(k, i) = (src(k - 1) + src(k) + src(k + 1)) / 3
         end do
         gvals(nz, i) = (2 * src(nz) + src(nz - 1)) / 3
         
         src(1) = (2 * gvals(1, i) + gvals(2, i)) / 3
         do k = 2, nz - 1
            src(k) = (gvals(k - 1, i) + gvals(k, i) + gvals(k + 1, i)) / 3
         end do
         src(nz) = (2 * gvals(nz, i) + gvals(nz - 1, i)) / 3
         
         gvals(1, i) = (2 * src(1) + src(2)) / 3
         do k = 2, nz - 1
            gvals(k, i) = (src(k - 1) + src(k) + src(k + 1)) / 3
         end do
         gvals(nz, i) = (2 * src(nz) + src(nz - 1)) / 3
      
      end do
   
   end subroutine smooth_gvals
   
   
   subroutine set_boundary_values(s, src, dest, j)
      type (star_info), pointer :: s
      real(dp), pointer :: src(:), dest(:, :)
      integer, intent(in) :: j
      integer :: k, nz
      nz = s% nz
      dest(1, j) = src(1)
      do k = 2, nz
         dest(k, j) = (src(k - 1) + src(k)) / 2
      end do
   end subroutine set_boundary_values
   
   
   subroutine check_validity(s, ierr)
      type (star_info), pointer :: s
      integer, intent(out) :: ierr
      
      integer :: k, nz
      
      include 'formats'
      
      ierr = 0
      nz = s% nz
      
      do k = 1, nz - 1
         if (s% xh(s% i_lnR, k) <= s% xh(s% i_lnR, k + 1)) then
            ierr = -1
            s% retry_message = 'at start of remesh: negative cell volume for cell'
            if (s% report_ierr) then
               write(*, *) 'at start of remesh: negative cell volume for cell', k
               write(*, *)
               write(*, 2) 's% xh(s% i_lnR, k)', k, s% xh(s% i_lnR, k)
               write(*, 2) 's% xh(s% i_lnR, k+1)', k + 1, s% xh(s% i_lnR, k + 1)
               write(*, *)
               write(*, *) 's% model_number', s% model_number
               write(*, *) 's% nz', s% nz
               write(*, *) 's% num_retries', s% num_retries
               write(*, *)
            end if
            return
         end if
         if (s% dq(k) <= 0) then
            ierr = -1
            s% retry_message = 'at start of remesh: non-positive cell mass for cell'
            if (s% report_ierr) then
               write(*, *) 'at start of remesh: non-positive cell mass for cell', k
               write(*, *) 's% model_number', s% model_number
               write(*, *) 's% nz', s% nz
               write(*, *) 's% num_retries', s% num_retries
               write(*, *)
            end if
            return
         end if
      end do
   
   end subroutine check_validity


end module adjust_mesh_support


