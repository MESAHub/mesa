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


! Routine eval_fp_ft for computing rotation corrections to the stellar structure equations.
! Following Endal & Sofia, 1976, ApJ 210:184.
! Based on code from Evert Glebbeek which in turn was based on code from Alex Heger.



      module hydro_rotation

      use const_def, only: pi, pi4, ln10, two_thirds, one_third
      use star_utils, only: get_r_from_xh

      use star_private_def

      implicit none

      ! Angular factors appearing in different integrals
      ! Precompute these for a quarter circle
      ! Read-only after initialization
      integer, parameter, private :: intmax = 101 ! divisions of 1/4 circle
      real(dp), private :: &
         dtheta, dtheta2, theta(intmax), cost(intmax),sint(intmax), &
         cost2(intmax),sint2(intmax), utheta(intmax)
      logical, save, private :: have_initialized = .false.

      logical, parameter :: dbg = .false.

      integer, parameter :: dbg_cell = -1


      contains

      ! compute w_div_w_roche for a known angular frequency omega, rphi, and Mphi
      real(dp) function w_div_w_roche_omega(rphi,Mphi,omega,cgrav, max_w, max_w2, w_div_wc_flag) result(w_roche)
         real(dp), intent(in) :: rphi,Mphi,omega,cgrav, max_w, max_w2
         logical, intent(in) :: w_div_wc_flag
         real(dp) :: wr, wr_high, wr_low, dimless_rphi, new_dimless_rphi, rphi_lim1, rphi_lim2

         if (omega == 0d0) then
            w_roche = 0d0
            return
         end if

         dimless_rphi = rphi*pow(abs(omega), two_thirds)/pow(cgrav*Mphi,one_third)
         if (.not. w_div_wc_flag) then
            ! verify if w_div_w_roche is not above max, otherwise limit it to that
            wr = max_w
            new_dimless_rphi = pow(wr,two_thirds)*(1-pow2(wr)/6d0+0.01726d0*pow4(wr)-0.03569d0*pow6(wr))
            if (dimless_rphi > new_dimless_rphi) then
               w_roche = wr
               return
            end if
         else
            ! smoothly cap to max_w to get a continuous function
            ! nothing is done when we are below max_w2, but between max_w2 and max_w we smoothly 
            ! produce an asymptote that would result in w_div_wc=max_w for jrot->infinity
            wr = max_w
            rphi_lim1 = pow(wr,two_thirds)*(1-pow2(wr)/6d0+0.01726d0*pow4(wr)-0.03569d0*pow6(wr))

            wr = max_w2
            rphi_lim2 = pow(wr,two_thirds)*(1-pow2(wr)/6d0+0.01726d0*pow4(wr)-0.03569d0*pow6(wr))

            if (abs(dimless_rphi) > rphi_lim1) then
               dimless_rphi = &
                 2*(rphi_lim2-rphi_lim1)/(1+exp(-2*(abs(dimless_rphi)-rphi_lim1)/(rphi_lim2-rphi_lim1)))-rphi_lim2+2*rphi_lim1
            end if
         end if

         ! otherwise, bisect result
         wr_high = wr
         wr_low = 0
         do while (wr_high-wr_low>1d-6)
            wr = 0.5d0*(wr_high+wr_low) 
            new_dimless_rphi = pow(wr,two_thirds)*(1-pow2(wr)/6d0+0.01726d0*pow4(wr)-0.03569d0*pow6(wr))
            if (dimless_rphi > new_dimless_rphi) then
               wr_low = wr
            else
               wr_high = wr
            end if
         end do
         w_roche = 0.5d0*(wr_high+wr_low)

         if (omega < 0d0) then
            w_roche = -w_roche
         end if

      end function w_div_w_roche_omega
      
      ! compute w_div_w_roche for a known specific angular momentum jrot, rphi, and Mphi
      real(dp) function w_div_w_roche_jrot(rphi,Mphi,jrot,cgrav, max_w, max_w2, w_div_wc_flag) result(w_roche)
         real(dp), intent(in) :: rphi,Mphi,jrot,cgrav, max_w, max_w2
         logical, intent(in) :: w_div_wc_flag
         real(dp) :: wr, wr_high, wr_low, dimless_factor, new_dimless_factor
         real(dp) :: w2, w4, w6, lg_one_sub_w4, jr_lim1, jr_lim2, A, C

         if (jrot == 0d0) then
            w_roche = 0d0
            return
         end if

         dimless_factor = abs(jrot)/sqrt(cgrav*Mphi*rphi)

         if (.not. w_div_wc_flag) then
            ! verify if w_div_w_roche is not above max, otherwise limit it to that
            wr = max_w
            w2 = pow2(wr)
            w4 = pow4(wr)
            w6 = pow6(wr)
            lg_one_sub_w4 = log(1d0-w4)
            new_dimless_factor = two_thirds*wr*(1+17d0/60d0*w2-0.3436d0*w4-0.4055d0*w6-0.9277d0*lg_one_sub_w4) &
                     /(1d0-0.1076d0*w4-0.2336d0*w6-0.5583d0*lg_one_sub_w4)
            if (dimless_factor > new_dimless_factor) then
               w_roche = wr
               return
            end if
         else
            ! smoothly cap to max_w to get a continuous function
            ! nothing is done when we are below max_w2, but between max_w2 and max_w we smoothly 
            ! produce an asymptote that would result in w_div_wc=max_w for jrot->infinity
            wr = max_w
            A = 1d0-0.1076d0*pow4(wr)-0.2336d0*pow6(wr)-0.5583d0*log(1d0-pow4(wr))
            C = 1d0+17d0/60d0*pow2(wr)-0.3436d0*pow4(wr)-0.4055d0*pow6(wr)-0.9277d0*log(1d0-pow4(wr))
            jr_lim1 = two_thirds*wr*C/A

            wr = max_w2
            A = 1d0-0.1076d0*pow4(wr)-0.2336d0*pow6(wr)-0.5583d0*log(1d0-pow4(wr))
            C = 1d0+17d0/60d0*pow2(wr)-0.3436d0*pow4(wr)-0.4055d0*pow6(wr)-0.9277d0*log(1d0-pow4(wr))
            jr_lim2 = two_thirds*wr*C/A

            if (abs(dimless_factor) > jr_lim1) then
               dimless_factor = 2*(jr_lim2-jr_lim1)/(1+exp(-2*(abs(dimless_factor)-jr_lim1)/(jr_lim2-jr_lim1)))-jr_lim2+2*jr_lim1
            end if
         end if

         ! otherwise, bisect result
         wr_high = wr
         wr_low = 0
         do while (wr_high-wr_low>1d-6)
            wr = 0.5d0*(wr_high+wr_low) 
            w2 = pow2(wr)
            w4 = pow4(wr)
            w6 = pow6(wr)
            lg_one_sub_w4 = log(1d0-w4)
            new_dimless_factor = two_thirds*wr*(1+17d0/60d0*w2-0.3436d0*w4-0.4055d0*w6-0.9277d0*lg_one_sub_w4) &
                     /(1d0-0.1076d0*w4-0.2336d0*w6-0.5583d0*lg_one_sub_w4)
            if (dimless_factor > new_dimless_factor) then
               wr_low = wr
            else
               wr_high = wr
            end if
         end do
         w_roche = 0.5d0*(wr_high+wr_low)

         if (jrot < 0d0) then
            w_roche = -w_roche
         end if

      end function w_div_w_roche_jrot

      subroutine eval_i_rot(s,k,r00,w_div_w_crit_roche, i_rot, di_rot_dlnr, di_rot_dw_div_wc)
         use auto_diff_support
         type (star_info), pointer :: s
         integer, intent(in) :: k ! just for debugging
         real(dp), intent(in) :: r00,w_div_w_crit_roche
         real(dp), intent(out) :: i_rot, di_rot_dlnr, di_rot_dw_div_wc

         type(auto_diff_real_2var_order1) :: ir, r, re, w, w2, w4, w6, lg_one_sub_w4, B, A

         include 'formats'
         if (s% use_other_eval_i_rot) then
            call s% other_eval_i_rot(s% id,k,r00,w_div_w_crit_roche, i_rot, di_rot_dlnr, di_rot_dw_div_wc)
         else if (s% simple_i_rot_flag) then
            i_rot = (2d0/3d0)*r00*r00
            di_rot_dlnr = 2*i_rot
            di_rot_dw_div_wc = 0d0
         else
            ! Compute i_rot following Paxton et al. 2019 (ApJs, 243, 10)
            w = w_div_w_crit_roche
            w%d1val1 = 1d0

            r = r00
            r%d1val2 = r00 ! Makes the independent variable lnR

            w2 = pow2(w)
            w4 = pow4(w)
            w6 = pow6(w)
            lg_one_sub_w4 = log(1d0-w4)
            re = r*(1d0+w2/6d0-0.0002507d0*w4+0.06075d0*w6)
            B = (1d0+w2/5d0-0.2735d0*w4-0.4327d0*w6-3d0/2d0*0.5583d0*lg_one_sub_w4)
            A = (1d0-0.1076d0*w4-0.2336d0*w6-0.5583d0*lg_one_sub_w4)
            
            ir =  two_thirds*pow2(re)*B/A

            i_rot = ir%val
            di_rot_dw_div_wc = ir%d1val1
            di_rot_dlnr = ir%d1val2
         end if

      end subroutine eval_i_rot


      subroutine set_i_rot(s, skip_w_div_w_crit_roche)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_w_div_w_crit_roche
         integer :: k, nz
         include 'formats'

!$OMP PARALLEL DO PRIVATE(k) SCHEDULE(dynamic,2)
         do k=1,s% nz
            if (.not. skip_w_div_w_crit_roche) then
               s% w_div_w_crit_roche(k) = &
                  w_div_w_roche_jrot(s% r(k),s% m(k),s% j_rot(k),s% cgrav(k), &
                     s% w_div_wcrit_max, s% w_div_wcrit_max2, s% w_div_wc_flag)
            end if
            call eval_i_rot(s, k, s% r(k), s% w_div_w_crit_roche(k), &
               s% i_rot(k), s% di_rot_dlnr(k), s% di_rot_dw_div_wc(k))
         end do
!$OMP END PARALLEL DO

      end subroutine set_i_rot

      subroutine set_i_rot_from_omega_and_j_rot(s)
         type (star_info), pointer :: s
         integer :: k
         include 'formats'
         do k=1,s% nz
            if (s% omega(k) /= 0d0) then
               ! we can directly compute i_rot using j_rot and omega
               s% i_rot(k) = s% j_rot(k)/s% omega(k)
            else
               call update1_i_rot_from_xh(s, k)
            end if
         end do
      end subroutine set_i_rot_from_omega_and_j_rot

      subroutine set_j_rot(s)
         type (star_info), pointer :: s
         integer :: k
         include 'formats'
         do k=1,s% nz
            s% j_rot(k) = s% i_rot(k)*s% omega(k)
         end do
      end subroutine set_j_rot


      subroutine set_omega(s, str)
         type (star_info), pointer :: s
         character (len=*) :: str
         integer :: k
         include 'formats'
         do k=1,s% nz
            s% omega(k) = s% j_rot(k)/s% i_rot(k)
         end do
      end subroutine set_omega


      subroutine check_omega(s, str)
         type (star_info), pointer :: s
         character (len=*) :: str
         integer :: k
         logical :: okay
         include 'formats'
         okay = .true.
         do k=1,s% nz
            if (abs(s% omega(k) - s% j_rot(k)/s% i_rot(k)) > 1d-14) then
               write(*,2) 'omega error', k, s% omega(k) - s% j_rot(k)/s% i_rot(k)
               okay = .false.
               exit
            end if
         end do
         if (okay) return
         write(*,*) trim(str)
         stop 'check_omega'
      end subroutine check_omega


      subroutine update1_i_rot_from_xh(s, k)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: r00, r003, rp1, rp13, rm1, rm13, r_in, r_out
         include 'formats'

         r00 = get_r_from_xh(s,k)

         call eval_i_rot(s, k, r00, s% w_div_w_crit_roche(k), &
            s% i_rot(k), s% di_rot_dlnr(k), s% di_rot_dw_div_wc(k))
      end subroutine update1_i_rot_from_xh

      subroutine use_xh_to_update_i_rot(s)
         type (star_info), pointer :: s
         integer :: k
         do k=1,s% nz
            if (s% j_rot(k) /= 0d0) then
               s% w_div_w_crit_roche(k) = &
                  w_div_w_roche_jrot(get_r_from_xh(s,k),s% m(k),s% j_rot(k),s% cgrav(k), &
                     s% w_div_wcrit_max, s% w_div_wcrit_max2, s% w_div_wc_flag)
            else
               s% w_div_w_crit_roche(k) = 0d0
            end if
         end do
         do k=1,s% nz
            call update1_i_rot_from_xh(s,k)
         end do
      end subroutine use_xh_to_update_i_rot

      subroutine use_xh_to_update_i_rot_and_j_rot(s)
         type (star_info), pointer :: s
         integer :: k
         do k=1,s% nz
            s% w_div_w_crit_roche(k) = &
               w_div_w_roche_omega(get_r_from_xh(s,k),s% m(k),s% omega(k),s% cgrav(k), &
                  s% w_div_wcrit_max, s% w_div_wcrit_max2, s% w_div_wc_flag)
         end do
         do k=1,s% nz
            call update1_i_rot_from_xh(s,k)
         end do
         call set_j_rot(s)
      end subroutine use_xh_to_update_i_rot_and_j_rot


      subroutine get_rotation_sigmas(s, nzlo, nzhi, dt, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         real(dp), intent(in) :: dt
         integer, intent(out) :: ierr

         integer :: k, nz
         real(dp), allocatable :: am_nu(:), am_sig(:)

         include 'formats'

         ierr = 0
         nz = s% nz
         
         allocate(am_nu(nz), am_sig(nz))

         call get1_am_sig(s, nzlo, nzhi, s% am_nu_j, s% am_sig_j, dt, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,1) 'failed in get_rotation_sigmas'
            return
         end if

         do k=1,nz
            am_nu(k) = s% am_nu_j(k) + s% am_nu_omega(k)
         end do
         ! do it this way so apply limit to sum; sum is used as diffusion coeff for omega
         call get1_am_sig(s, nzlo, nzhi, am_nu, am_sig, dt, ierr)
         if (ierr /= 0) then
            if (s% report_ierr) write(*,1) 'failed in get_rotation_sigmas'
            return
         end if

         do k=1,nz
            s% am_sig_omega(k) = max(0d0, am_sig(k) - s% am_sig_j(k))
         end do

      end subroutine get_rotation_sigmas


      subroutine get1_am_sig(s, nzlo, nzhi, am_nu, am_sig, dt, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nzlo, nzhi
         real(dp), intent(in) :: dt
         real(dp), dimension(:) :: am_nu, am_sig
         integer, intent(out) :: ierr

         integer :: k, nz, nz2
         real(dp) :: r, D, am_nu_E00, am_nu_Ep1, dmbar, s1, &
            sig_term_limit, xmstar, siglim

         include 'formats'

         ierr = 0
         xmstar = s% xmstar
         sig_term_limit = s% am_sig_term_limit
         nz = s% nz
         ! note: am_sig is cell centered, so combine adjacent am_nu face values.
         am_nu_E00 = 0; am_nu_Ep1 = 0
         nz2 = nzhi
         if (nzhi == nz) then
            k = nz
            D = am_nu_E00
            r = 0.5d0*s% r(k)
            s1 = pi4*r*r*s% rho(k)
            am_sig(k) = s1*s1*D/s% dm(k)
            nz2 = nz-1
         end if
         do k = nzlo, nz2
            am_nu_E00 = max(0d0, am_nu(k))
            am_nu_Ep1 = max(0d0, am_nu(k+1))
            ! Meynet, Maeder, & Mowlavi, A&A 416, 1023-1036, 2004, eqn 51 with f = 1/2.
            D = 2*(am_nu_E00*am_nu_Ep1)/max(1d-99, am_nu_E00 + am_nu_Ep1)
            r = 0.5d0*(s% r(k) + s% r(k+1)) ! consistent with f = 1/2
            s1 = pi4*r*r*s% rho(k)
            am_sig(k) = s1*s1*D/s% dm(k)
         end do

         ! can get numerical problems unless limit am_sig
         ! adjust am_sig to make sure am_sig*dt/dmbar is < allowed limit
         do k = nzlo, nzhi
            if (k < nz) then
               dmbar = xmstar*min(s% dq(k),s% dq(k+1))
               siglim = sig_term_limit*dmbar/dt
            else
               dmbar = xmstar*s% dq(k)
               siglim = sig_term_limit*dmbar/dt
            end if
            if (am_sig(k) > siglim) then
               am_sig(k) = siglim
            end if
         end do

      end subroutine get1_am_sig


      subroutine set_uniform_omega(id, omega, ierr)
         integer, intent(in) :: id
         real(dp), intent(in) :: omega
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, nz
         include 'formats'
         call get_star_ptr(id, s, ierr)
         if (ierr /= 0) return
         nz = s% nz
         do k=1, nz
            s% omega(k) = omega
         end do
         call use_xh_to_update_i_rot_and_j_rot(s)
         s% fp_rot(1:nz) = 0
         s% ft_rot(1:nz) = 0
         s% w_div_w_crit_roche(1:nz) = 0 !TODO: this should not be here, but removing it breaks a test case
         s% r_polar(1:nz) = 0
         s% r_equatorial(1:nz) = 0
         s% am_nu_rot(1:nz) = 0
         s% am_nu_non_rot(1:nz) = 0
         s% am_nu_omega(1:nz) = 0
         s% am_nu_j(1:nz) = 0
         s% am_sig_omega(1:nz) = 0
         s% am_sig_j(1:nz) = 0
         s% domega_dlnR(1:nz) = 0
         s% richardson_number(1:nz) = 0
         s% D_mix_non_rotation(1:nz) = 0
         s% D_visc(1:nz) = 0
         s% D_DSI(1:nz) = 0
         s% D_SH(1:nz) = 0
         s% D_SSI(1:nz) = 0
         s% D_ES(1:nz) = 0
         s% D_GSF(1:nz) = 0
         s% D_ST(1:nz) = 0
         s% nu_ST(1:nz) = 0
         s% omega_shear(1:nz) = 0
         s% dynamo_B_r(1:nz) = 0
         s% dynamo_B_phi(1:nz) = 0
         call set_rotation_info(s, .false., ierr)
         if (ierr /= 0) return
         s% need_to_setvars = .true.
      end subroutine set_uniform_omega


      subroutine set_rotation_info(s, skip_w_div_w_crit_roche, ierr)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_w_div_w_crit_roche
         integer, intent(out) :: ierr
         integer :: k
         include 'formats'
         ierr = 0

         if (.not. s% rotation_flag) return

         call set_i_rot(s, skip_w_div_w_crit_roche)
         call set_omega(s, 'set_rotation_info')

         if (.not. s% use_other_eval_fp_ft) then
            call eval_fp_ft( &
                  s% id, s% nz, s% m, s% r, s% rho, s% omega, s% ft_rot, s% fp_rot, &
                  s% r_polar, s% r_equatorial, s% report_ierr, ierr)
         else
            call s% other_eval_fp_ft( &
                  s% id, s% nz, s% m, s% r, s% rho, s% omega, s% ft_rot, s% fp_rot, &
                  s% r_polar, s% r_equatorial, s% report_ierr, ierr)
         end if
         if (ierr /= 0) then
            write(*,*) 'failed in eval_fp_ft'
         end if
      end subroutine set_rotation_info


      subroutine set_surf_avg_rotation_info(s)
         use star_utils, only: get_Lrad_div_Ledd
         type (star_info), pointer :: s
         real(dp) :: &
            dm, dmsum, omega_sum, omega_crit_sum, omega_div_omega_crit_sum, &
            v_rot_sum, v_crit_sum, v_div_v_crit_sum, Lrad_div_Ledd_sum, &
            kap_face, Ledd, gamma_factor, omega_crit, omega, kap_sum, &
            j_rot_sum, j_rot, v_rot, v_crit, Lrad_div_Ledd, dtau, tau, &
            cgrav, kap, mmid, Lmid, rmid, logT_sum, logRho_sum
         integer :: k, ierr
         logical, parameter :: dbg = .false.
         include 'formats'

         if (.not. s% rotation_flag) then
            s% omega_avg_surf = 0
            s% omega_crit_avg_surf = 0
            s% w_div_w_crit_avg_surf = 0
            s% j_rot_avg_surf = 0
            s% v_rot_avg_surf = 0
            s% v_crit_avg_surf = 0
            s% v_div_v_crit_avg_surf = 0
            s% Lrad_div_Ledd_avg_surf = 0
            s% opacity_avg_surf = 0
            s% logT_avg_surf = 0
            s% logRho_avg_surf = 0
            return
         end if
         
         ierr = 0
         call set_rotation_info(s,.true.,ierr)
         if (ierr /= 0) then
            write(*,*) 'got ierr from call set_rotation_info in set_surf_avg_rotation_info'
            write(*,*) 'just ignore it'
         end if

         tau = s% tau_factor*s% tau_base
         dmsum = 0d0
         Lrad_div_Ledd_sum = 0d0
         rmid = 0d0

         do k = 1, s% nz - 1
            kap = s% opacity(k)
            rmid = s% rmid(k)
            mmid = 0.5d0*(s% m_grav(k) + s% m_grav(k+1))
            Lmid = 0.5d0*(s% L(k) + s% L(k+1))
            cgrav = 0.5d0*(s% cgrav(k) + s% cgrav(k+1))
            dm = s% dm(k)
            dtau = dm*kap/(pi4*rmid*rmid)

            if (tau + dtau <= s% surf_avg_tau_min) then
               tau = tau + dtau
               cycle
            end if

            ! check for partial contribution from cell
            ! the tau < s% surf_avg_tau is meant for the case in which the surface tau is set
            ! equal or larger to surf_avg_tau. In that case we just use the values of the surface cell.
            if (tau < s% surf_avg_tau) then
               if (tau < s% surf_avg_tau_min) then ! only use part of this cell
                  dm = dm*(tau + dtau - s% surf_avg_tau_min)/dtau
               else if (tau + dtau > s% surf_avg_tau) then ! only use part of this cell
                  dm = dm*(s% surf_avg_tau - tau)/dtau
                  !write(*,2) 'tau limit', k, (s% surf_avg_tau - tau)/dtau
               end if
            end if
            dmsum = dmsum + dm
            Lrad_div_Ledd = get_Lrad_div_Ledd(s,k)
            Lrad_div_Ledd_sum = Lrad_div_Ledd_sum + dm*Lrad_div_Ledd
            tau = tau + dtau
            if (tau >= s% surf_avg_tau) exit
         end do

         s% Lrad_div_Ledd_avg_surf = Lrad_div_Ledd_sum/dmsum
         gamma_factor = 1d0 - min(s% Lrad_div_Ledd_avg_surf, 0.9999d0)

         tau = s% tau_factor*s% tau_base
         dmsum = 0
         j_rot_sum = 0
         omega_sum = 0
         omega_crit_sum = 0
         omega_div_omega_crit_sum = 0
         v_rot_sum = 0
         v_crit_sum = 0
         v_div_v_crit_sum = 0
         kap_sum = 0
         logT_sum = 0
         logRho_sum = 0

         do k = 1, s% nz - 1

            kap = s% opacity(k)
            ! TODO: better explain
            ! Use equatorial radius
            rmid = 0.5d0*(s% r_equatorial(k) + s% r_equatorial(k+1))
            dm = s% dm(k)
            dtau = dm*kap/(pi4*rmid*rmid)

            if (tau + dtau <= s% surf_avg_tau_min) then
               tau = tau + dtau
               cycle
            end if

            ! check for partial contribution from cell
            ! the tau < s% surf_avg_tau is meant for the case in which the surface tau is set
            ! equal or larger to surf_avg_tau. In this case we just use the values of the surface cell.
            if (tau < s% surf_avg_tau) then
               if (tau < s% surf_avg_tau_min) then ! only use part of this cell
                  dm = dm*(tau + dtau - s% surf_avg_tau_min)/dtau
               else if (tau + dtau > s% surf_avg_tau) then ! only use part of this cell
                  dm = dm*(s% surf_avg_tau - tau)/dtau
               end if
            end if

            dmsum = dmsum + dm
            cgrav = 0.5d0*(s% cgrav(k) + s% cgrav(k+1))
            mmid = 0.5d0*(s% m_grav(k) + s% m_grav(k+1))
            omega = 0.5d0*(s% omega(k) + s% omega(k+1))
            j_rot = 0.5d0*(s% j_rot(k) + s% j_rot(k+1))

            kap_sum = kap_sum + dm*kap
            j_rot_sum = j_rot_sum + dm*j_rot

            omega_crit = sqrt(gamma_factor*cgrav*mmid/pow3(rmid))
            omega_div_omega_crit_sum = omega_div_omega_crit_sum + dm*abs(omega/omega_crit)

            v_rot = omega*rmid
            v_crit = omega_crit*rmid
            omega_sum = omega_sum + dm*omega
            omega_crit_sum = omega_crit_sum + dm*omega_crit
            v_rot_sum = v_rot_sum + dm*v_rot
            v_crit_sum = v_crit_sum + dm*v_crit
            v_div_v_crit_sum = v_div_v_crit_sum + dm*abs(v_rot/v_crit)
            logT_sum = logT_sum + dm*s% lnT(k)/ln10
            logRho_sum = logRho_sum + dm*s% lnd(k)/ln10
            kap_sum = kap_sum + dm*kap
            tau = tau + dtau
            if (tau >= s% surf_avg_tau) exit

         end do

         s% logT_avg_surf = logT_sum/dmsum
         s% logRho_avg_surf = logRho_sum/dmsum
         s% opacity_avg_surf = kap_sum/dmsum
         s% j_rot_avg_surf = j_rot_sum/dmsum
         s% omega_avg_surf = omega_sum/dmsum
         s% omega_crit_avg_surf = omega_crit_sum/dmsum
         s% w_div_w_crit_avg_surf = omega_div_omega_crit_sum/dmsum
         s% v_rot_avg_surf = v_rot_sum/dmsum
         s% v_crit_avg_surf = v_crit_sum/dmsum
         s% v_div_v_crit_avg_surf = v_div_v_crit_sum/dmsum

      end subroutine set_surf_avg_rotation_info


      ! Input variables:
      !  N     Number of meshpoints used by the model (arrays are this size)
      !  XM    Mass coordinate [gram]
      !  R     Radius coordinate [cm]
      !  RHO   Density [gram/cm^3]
      !  AW    Angular velocity [rad/sec]
      ! Output variables:
      !  Correction factor FT at each meshpoint
      !  Correction factor FP at each meshpoint
      !  r_polar, r_equatorial at each meshpoint
      subroutine eval_fp_ft( &
            id, nz, xm, r, rho, aw, ft, fp, r_polar, r_equatorial, report_ierr, ierr)
         use num_lib
         integer, intent(in) :: id
         integer, intent(in) :: nz
         real(dp), intent(in) :: aw(:), r(:), rho(:), xm(:) ! (nz)
         real(dp), intent(inout) :: ft(:), fp(:), r_polar(:), r_equatorial(:) ! (nz)
         logical, intent(in) :: report_ierr
         integer, intent(out) :: ierr

         type (star_info), pointer :: s
         integer :: j, k, kmax, ift_in, ift_out, ifp_in, ifp_out, kanz, i_in, i_out
         real(dp) :: gur(nz),eta(nz),rnull(nz),geff(nz),spgp(nz),spgm(nz), &
            gp(intmax),gmh(intmax)
         real(dp) :: rae(intmax),rov,etanu,wr,xm_out,&
            a,ft_min,fp_min,tegra,rnorm,rom1,rho_out,rho_in,rnull_out,rnull_in, &
            r_in,r_out,lnr_in, lnr_out, dat, dot, dut, r1, r2, dfdx, rnuv1, rnuv2, rrs, rrs1

         integer, parameter :: lipar = 2, lrpar = 5, imax = 50
         real(dp), parameter :: epsx=1d-8, epsy=epsx
         real(dp) :: apot1, apot2

         integer :: liwork, lwork
         integer, pointer :: iwork(:)
         real(dp), pointer :: work(:)

         real(dp), pointer :: cgrav(:)

         integer, parameter :: out_io = 0, max_steps = 1000, itol = 0, lout = 0
         real(dp) :: rtol(1), atol(1), max_step_size, h
         integer :: idid

         real(dp), target :: veta_ary(1)
         real(dp), pointer :: veta(:)

         real(dp), target :: rpar_ary(lrpar)
         integer, target :: ipar_ary(lipar)
         real(dp), pointer :: rpar(:)
         integer, pointer :: ipar(:)

         logical :: dbg

         real(dp) :: A_omega,fp_numerator, ft_numerator, w, w2, w3, w4, w5, w6, lg_one_sub_w4, &
            d_A_omega_dw, d_fp_numerator_dw, d_ft_numerator_dw

         include 'formats'

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         rpar => rpar_ary
         ipar => ipar_ary
         cgrav => s% cgrav

         dbg = .false. ! (s% model_number >= 5)

!$OMP PARALLEL DO PRIVATE(j, A_omega, fp_numerator, ft_numerator, d_A_omega_dw, d_fp_numerator_dw, d_ft_numerator_dw, w, w2, w3, w4, w5, w6, lg_one_sub_w4) SCHEDULE(dynamic,2)
            do j=1, s% nz
               !Compute fp, ft, re and rp using fits to the Roche geometry of a single star.
               !by this point in the code, w_div_w_crit_roche is set
               w2 = pow2(s% w_div_w_crit_roche(j))
               w4 = pow4(s% w_div_w_crit_roche(j))
               w6 = pow6(s% w_div_w_crit_roche(j))
               lg_one_sub_w4 = log(1d0-w4)
               A_omega = (1d0-0.1076d0*w4-0.2336d0*w6-0.5583d0*lg_one_sub_w4)
               fp_numerator = (1d0-two_thirds*w2-0.06837d0*w4-0.2495d0*w6)
               ft_numerator = (1d0+0.2185d0*w4-0.1109d0*w6)
               !fits for fp, ft
               fp(j) = fp_numerator/A_omega
               ft(j) = ft_numerator/A_omega
               !re and rp can be derived analytically from w_div_wcrit
               r_equatorial(j) = r(j)*(1d0+w2/6d0-0.0002507d0*w4+0.06075d0*w6)
               r_polar(j) = r_equatorial(j)/(1d0+0.5d0*w2)
               ! Be sure they are consistent with r_Phi
               r_equatorial(j) = max(r_equatorial(j),r(j))
               r_polar(j) = min(r_polar(j),r(j))
               if (s% w_div_wc_flag) then
                  ! need to compute partials as well
                  w = s% w_div_w_crit_roche(j)
                  w3 = pow3(s% w_div_w_crit_roche(j))
                  w5 = pow5(s% w_div_w_crit_roche(j))
                  d_fp_numerator_dw = -2d0*two_thirds*w-4d0*0.06837d0*w3-6d0*0.2495d0*w5
                  d_ft_numerator_dw = 4d0*0.2185d0*w3-6d0*0.1109d0*w5
                  d_A_omega_dw = -4d0*0.1076d0*w3-6d0*0.2336d0*w5-0.5583d0/(1d0-w4)*(-4d0*w3)
                  s% dfp_rot_dw_div_wc(j) = (d_fp_numerator_dw/A_omega - fp_numerator*d_A_omega_dw/pow2(A_omega))
                  s% dft_rot_dw_div_wc(j) = (d_ft_numerator_dw/A_omega - ft_numerator*d_A_omega_dw/pow2(A_omega))
               end if
               !if (j == s% solver_test_partials_k) then
               !   s% solver_test_partials_val = fp(j)
               !   s% solver_test_partials_var = s% i_w_div_wc
               !   s% solver_test_partials_dval_dx = s% dfp_rot_dw_div_wc(j)
               !end if
            end do
!$OMP END PARALLEL DO

      end subroutine eval_fp_ft

      subroutine compute_j_fluxes_and_extra_jdot(id, ierr)
         use auto_diff_support
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         type(auto_diff_real_star_order1) :: omega00, omegap1, r00, rp1, i_rot00, i_rotp1, rho00, part1

         integer :: k
         real(dp) :: pi2_div4
         ierr = 0

         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (s% am_D_mix_factor==0d0) then
            s% am_nu_omega(:) = 0d0
         end if

         s% extra_jdot(:) = 0d0
         if (s% use_other_torque) then
            call s% other_torque(s% id, ierr)
            if (ierr /= 0) then
               if (s% report_ierr .or. dbg) &
                  write(*, *) 'solve_omega_mix: other_torque returned ierr', ierr
               return
            end if
         end if
         if (associated(s% binary_other_torque)) then
            call s% binary_other_torque(s% id, ierr)
            if (ierr /= 0) then
               if (s% report_ierr .or. dbg) &
                  write(*, *) 'solve_omega_mix: binary_other_torque returned ierr', ierr
               return
            end if
         end if

         pi2_div4 = pow2(pi)/4d0
         do k=1, s% nz-1

            r00 = wrap_r_00(s, k)
            rp1 = wrap_r_p1(s, k)
            rho00 = wrap_d_00(s, k)

            omega00 = 0d0
            omega00% val = s% omega(k)   
            omega00% d1Array(i_jrot_00) = 1d0/s% i_rot(k)
            omega00% d1Array(i_lnR_00) = -s% j_rot(k)/pow2(s% i_rot(k))*s% di_rot_dlnr(k)
            omega00% d1Array(i_w_div_wc_00) = -s% j_rot(k)/pow2(s% i_rot(k))*s% di_rot_dw_div_wc(k)
            omegap1 = 0d0
            omegap1% val = s% omega(k+1)   
            omegap1% d1Array(i_jrot_p1) = 1d0/s% i_rot(k+1)
            omegap1% d1Array(i_lnR_p1) = -s% j_rot(k+1)/pow2(s% i_rot(k+1))*s% di_rot_dlnr(k+1)
            omegap1% d1Array(i_w_div_wc_p1) = -s% j_rot(k+1)/pow2(s% i_rot(k+1))*s% di_rot_dw_div_wc(k+1)

            i_rot00 = 0d0
            i_rot00% val = s% i_rot(k)   
            i_rot00% d1Array(i_w_div_wc_00) = s% di_rot_dw_div_wc(k)
            i_rot00% d1Array(i_lnR_00) = s% di_rot_dlnr(k)
            i_rotp1 = 0d0
            i_rotp1% val = s% i_rot(k+1)
            i_rotp1% d1Array(i_w_div_wc_p1) = s% di_rot_dw_div_wc(k+1)
            i_rotp1% d1Array(i_lnR_p1) = s% di_rot_dlnr(k+1)

            part1 = -pi2_div4*pow4(r00+rp1)*pow2(rho00)*(i_rot00+i_rotp1)*(s% am_nu_omega(k)+s% am_nu_omega(k+1))
            s% j_flux(k) = part1*(omega00-omegap1)/s% dm(k)

            !! this is to test partials
            !if (k==188) then
            !   part1 = (omega00-omegap1)/s% dm(k)
            !   s% solver_test_partials_val = part1% val
            !   s% solver_test_partials_var = s% i_lnR !s% i_w_div_wc
            !   s% solver_test_partials_dval_dx =  part1% d1Array(i_lnR_00) !part1% d1Array(i_w_div_wc_00)
            !end if

         end do
         s% j_flux(s% nz) = 0d0
      end subroutine compute_j_fluxes_and_extra_jdot

      subroutine init_rotation(ierr)
         integer, intent(out) :: ierr
         integer :: i
         ierr = 0
         if (have_initialized) return
         dtheta=pi / (2.d0 * (intmax-1))
         dtheta2=0.5d0*dtheta
         do i=1,intmax
            theta(i)=(i-1)*dtheta
            cost(i)=cos(theta(i))
            cost2(i)=cost(i)*cost(i)
            sint(i)=sin(theta(i))
            sint2(i)=sint(i)*sint(i)
            utheta(i)=0.5d0*(3.d0*cost2(i)-1.d0)
         end do
         have_initialized = .true.
      end subroutine init_rotation


      end module hydro_rotation

