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

      ! inner radius of shell ri
      ! outer radius of shell ra
      subroutine eval_i_rot(s,k,ri,r00,ra,w_div_w_crit_roche, i_rot, di_rot_dlnr, di_rot_dw_div_wc)
         type (star_info), pointer :: s
         integer, intent(in) :: k ! just for debugging
         real(dp), intent(in) :: ri,r00,ra,w_div_w_crit_roche
         real(dp), intent(out) :: i_rot, di_rot_dlnr, di_rot_dw_div_wc
         real(dp) :: re, w, w2, w3, w4, w5, w6, lg_one_sub_w4, B, A
         real(dp) :: rai,ra2,ri2,rm2
         real(dp) :: dre_dw_div_wc, dB_dw_div_wc, dA_dw_div_wc
         include 'formats'
         if (s% use_other_eval_i_rot) then
            call s% other_eval_i_rot(s% id,ri,r00,ra,w_div_w_crit_roche, i_rot, di_rot_dlnr, di_rot_dw_div_wc)
         else if (s% fitted_fp_ft_i_rot) then
            !If fitted_fp_ft_i_rot is true, then compute i_rot following Paxton et al. 2019 (ApJs, 243, 10)
            w = w_div_w_crit_roche
            w2 = pow2(w_div_w_crit_roche)
            w3 = pow3(w_div_w_crit_roche)
            w4 = pow4(w_div_w_crit_roche)
            w5 = pow5(w_div_w_crit_roche)
            w6 = pow6(w_div_w_crit_roche)
            lg_one_sub_w4 = log(1d0-w4)
            re = r00*(1d0+w2/6d0-0.0002507d0*w4+0.06075d0*w6)
            dre_dw_div_wc = r00*(w/3d0-4d0*0.0002507d0*w3+6d0*0.06075d0*w5)
            B = (1d0+w2/5d0-0.2735d0*w4-0.4327d0*w6-3d0/2d0*0.5583d0*lg_one_sub_w4)
            dB_dw_div_wc = 2d0*w/5d0-4d0*0.2735d0*w3-6d0*0.4327d0*w5+3d0/2d0*0.5583d0/(1d0-w4)*4d0*w3
            A = (1d0-0.1076d0*w4-0.2336d0*w6-0.5583d0*lg_one_sub_w4)
            dA_dw_div_wc = -4d0*0.1076d0*w3-6d0*0.2336d0*w5+0.5583d0/(1d0-w4)*4d0*w3

            i_rot =  two_thirds*pow2(re)*B/A
            di_rot_dw_div_wc = i_rot*(2d0*dre_dw_div_wc/re + dB_dw_div_wc/B - dA_dw_div_wc/A)
            di_rot_dlnr = 2*i_rot
         else if (s% simple_i_rot_flag) then
            i_rot = (2d0/3d0)*r00*r00
            di_rot_dlnr = 2*i_rot
            di_rot_dw_div_wc = 0d0
         else
            ! expression for evaluation without subtraction from Langer code
            rai=ra*ri
            ra2=ra*ra
            ri2=ri*ri
            rm2=ri2+rai+ra2
            i_rot=0.4D0*(ri2*ri2+rai*rm2+ra2*ra2)/rm2
            di_rot_dlnr = 0d0 ! not supperted for implicit solver
            di_rot_dw_div_wc = 0d0 ! not supperted for implicit solver
         end if

      end subroutine eval_i_rot


      subroutine set_i_rot(s, skip_w_div_w_crit_roche)
         type (star_info), pointer :: s
         logical, intent(in) :: skip_w_div_w_crit_roche
         integer :: k, nz
         include 'formats'

!$OMP PARALLEL DO PRIVATE(k) SCHEDULE(dynamic,2)
         do k=1,s% nz
            if (s% fitted_fp_ft_i_rot .and. .not. skip_w_div_w_crit_roche) then
               s% w_div_w_crit_roche(k) = &
                  w_div_w_roche_jrot(s% r(k),s% m(k),s% j_rot(k),s% cgrav(k), &
                     s% w_div_wcrit_max, s% w_div_wcrit_max2, s% w_div_wc_flag)
            end if
            if (k==1) then
               call eval_i_rot(s, k, s% rmid(1), s% r(1), s% r(1), s% w_div_w_crit_roche(1), &
                  s% i_rot(1), s% di_rot_dlnr(1), s% di_rot_dw_div_wc(1))
            else
               call eval_i_rot(s, k, s% rmid(k), s% r(k), s% rmid(k-1), s% w_div_w_crit_roche(k), &
                  s% i_rot(k), s% di_rot_dlnr(k), s% di_rot_dw_div_wc(k))
            end if
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

         if (s% fitted_fp_ft_i_rot .or. s% simple_i_rot_flag) then
            call eval_i_rot(s, k, r00, r00, r00, s% w_div_w_crit_roche(k), &
               s% i_rot(k), s% di_rot_dlnr(k), s% di_rot_dw_div_wc(k))
            return
         end if

         ! Compute the moment of inertia of a thin shell, ignoring rotational deformation
         ! need to compute rmid at k+1 and k-1. These are respectively r_in and r_out
         r00 = get_r_from_xh(s,k)
         r003 = r00*r00*r00

         if (k == s% nz) then
            rp1 = s% R_center
         else
            rp1 = get_r_from_xh(s,k+1)
         end if
         rp13 = rp1*rp1*rp1
         r_in = pow(0.5*(r003 + rp13),one_third)

         if (k == 1) then
            r_out = r00
         else
            rm1 = get_r_from_xh(s,k-1)
            rm13 = rm1*rm1*rm1
            r_out = pow(0.5*(r003 + rm13),one_third)
         end if

         call eval_i_rot(s,k,r_in,r00,r_out,0d0,&
            s% i_rot(k), s% di_rot_dlnr(k), s% di_rot_dw_div_wc(k))

      end subroutine update1_i_rot_from_xh

      subroutine use_xh_to_update_i_rot(s)
         type (star_info), pointer :: s
         integer :: k
         if (s% fitted_fp_ft_i_rot) then
            do k=1,s% nz
               if (s% j_rot(k) /= 0d0) then
                  s% w_div_w_crit_roche(k) = &
                     w_div_w_roche_jrot(get_r_from_xh(s,k),s% m(k),s% j_rot(k),s% cgrav(k), &
                        s% w_div_wcrit_max, s% w_div_wcrit_max2, s% w_div_wc_flag)
               else
                  s% w_div_w_crit_roche(k) = 0d0
               end if
            end do
         end if
         do k=1,s% nz
            call update1_i_rot_from_xh(s,k)
         end do
      end subroutine use_xh_to_update_i_rot

      subroutine use_xh_to_update_i_rot_and_j_rot(s)
         type (star_info), pointer :: s
         integer :: k
         if (s% fitted_fp_ft_i_rot) then
            do k=1,s% nz
               s% w_div_w_crit_roche(k) = &
                  w_div_w_roche_omega(get_r_from_xh(s,k),s% m(k),s% omega(k),s% cgrav(k), &
                     s% w_div_wcrit_max, s% w_div_wcrit_max2, s% w_div_wc_flag)
            end do
         end if
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
            if (s% fitted_fp_ft_i_rot) then
               ! TODO: better explain
               ! Use equatorial radius
               rmid = 0.5d0*(s% r_equatorial(k) + s% r_equatorial(k+1))
            else
              rmid = s% rmid(k)
            end if
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

         if (s% fitted_fp_ft_i_rot) then
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
            return
         end if

         s% dfp_rot_dw_div_wc(:s% nz) = 0d0
         s% dft_rot_dw_div_wc(:s% nz) = 0d0

         veta => veta_ary
         if (.not. have_initialized) then
            write(*,*) 'must call init_rotation prior to getting rotation info'
            ierr = -1; return
         end if

         kmax = 0

         ft_min=1.0d0
         fp_min=1.0d0
         ift_in=0
         ift_out=0
         ifp_in=0
         ifp_out=0

         call cash_karp_work_sizes(1,liwork,lwork)
         allocate(work(lwork),iwork(liwork))

         ! main loop
         etanu = -1.00d-20
         tegra=0.d0
         rnull_out=0d0
         rho_out = rho(nz)
         r_out = 0
         lnr_out = 1d-10

         rnorm=1.d0
         do i_in=nz+1,2,-1
            i_out = i_in-1

            ! for each shell, compute the distortion (in terms of r0) of the shell. The
            ! inner edge of the shell is at r(i_in), the outer edge at r(i_out). The density is
            ! defined in the interior of each shell and needs to be interpolated when
            ! it is needed at the edge.


            rho_in = rho_out
            rho_out = rho(i_out)
            ! Compute mean density at this point
            r_in = r_out
            lnr_in = lnr_out
            r_out = r(i_out)
            lnr_out = log(r_out)
            rom1 = xm(i_out)*0.75d0/(pi*r_out*r_out*r_out)
            rov = rho_in/rom1
            veta(1) = etanu
            kanz = 1
            dat = 1.d-04
            dot = (lnr_out-lnr_in)/8.d0
            dut = (lnr_out-lnr_in)/300.d0

            ! E&S, A6 by integration of A7
            rtol = 1d-6
            atol = 1d-6
            h = lnr_out - lnr_in
            max_step_size = 0d0
            rpar(1) = rov
            call cash_karp( &
               1, rad_fcn, lnr_in, veta, lnr_out, &
               h, max_step_size, max_steps, &
               rtol, atol, itol, &
               null_solout, out_io, &
               work, lwork, iwork, liwork, &
               lrpar, rpar, lipar, ipar, &
               lout, idid)
            if (idid < 0) then
               if (report_ierr .or. dbg) write(*,2) 'eval_fp_ft failed in integration', i_out
               ierr = -1
               if (dbg) stop 'eval_fp_ft'
               exit
            end if

            eta(i_out)=veta(1)
            xm_out=xm(i_out)
            wr=aw(i_out)
            etanu=veta(1)
            r1=4.00d0*r_out
            r2=0.01d0*r_out
            do j = 1, 400 ! exit when have bracketed root
               ! needs etanu,wr,r_out,xm_out; sets a
               rpar(1) = etanu
               rpar(2) = wr
               rpar(3) = r_out
               rpar(4) = xm_out
               ipar(1) = i_out
               ipar(2) = 0
               ierr = 0
               apot1 = apot_fcn(r1, dfdx, lrpar, rpar, lipar, ipar, ierr)
               if (ierr /= 0) then
                  if (report_ierr .or. dbg) write(*,2) 'eval_fp_ft failed in calculation of apot1', i_out
                  if (dbg) stop 'eval_fp_ft'
                  exit
               end if
               apot2 = apot_fcn(r2, dfdx, lrpar, rpar, lipar, ipar, ierr)
               if (ierr /= 0) then
                  if (report_ierr .or. dbg) write(*,2) 'eval_fp_ft failed in calculation of apot2', i_out
                  if (dbg) stop 'eval_fp_ft'
                  exit
               end if
               if (apot1*apot2 <= 0) exit
               r1=4.00d0*r1
               r2=0.01d0*r2
            end do

            if (apot1*apot2 > 0) then
               ierr = -1
               if (report_ierr .or. dbg) &
                  write(*,2) 'eval_fp_ft failed in calculation of apot1, apot2', i_out, apot1, apot2


               !exit
               write(*,2) 'i_out', i_out
               write(*,1) 'r_out', r_out
               write(*,1) 'r1', r1
               write(*,1) 'r2', r2
               write(*,1) 'apot1', apot1
               write(*,1) 'apot2', apot2
               write(*,1) 'epsx', epsx
               write(*,1) 'epsy', epsy
               write(*,2) 'imax', imax
               stop 'hydro_rotation: eval_fp_ft'


               exit
            end if

            rnuv1 = safe_root_with_initial_guess( &
               apot_fcn, r_out, r1, r2, apot1, apot2, &
               imax, epsx, epsy, lrpar, rpar, lipar, ipar, ierr)
            if (ierr /= 0) then
               if (report_ierr .or. dbg) write(*,2) 'eval_fp_ft failed in calculation of apot', i_out
               if (.not. dbg) exit
               write(*,2) 'i_out', i_out
               write(*,1) 'r_out', r_out
               write(*,1) 'r1', r1
               write(*,1) 'r2', r2
               write(*,1) 'apot1', apot1
               write(*,1) 'apot2', apot2
               write(*,1) 'epsx', epsx
               write(*,1) 'epsy', epsy
               write(*,2) 'imax', imax
               stop 'hydro_rotation: eval_fp_ft'
            end if

            a = rpar(5)
            rnuv2=rnuv1
            rnull(i_out) =  rnuv2
            rnull_in = rnull_out
            rnull_out = rnuv2

            ! compute integral for the potential calculation. See Endal&Sofia for details

            tegra = psiint(0.5d0*(rho_out+rho_in),xm_out,wr,etanu,rnull_in,rnull_out) + tegra

            ! calculate g and 1/g on a quarter circle.
            gur(i_out)=0.0d0
            do k=1,intmax,1
               gp(k)=gpsi(cgrav,rnull(i_out), tegra, k, a, wr, xm_out, rae)
               gmh(k)=1.d0 / gp(k)
               gur(i_out)=gur(i_out)+gp(k)
            end do
            gur(i_out)=gur(i_out)*dtheta/pi*2.d0
            r_polar(i_out) = max(rae(1),r2)
            r_equatorial(i_out) = min(rae(intmax),r1)

            ! find spsi*<g> and spsi*<1/g>. spsi is the surface area of the equipotential
            rrs=rae(1)*rae(1)*sint(1)
            rrs1=rae(intmax)*rae(intmax)*sint(intmax)
            spgp(i_out)=(gp(1)*rrs+gp (intmax)*rrs1)*dtheta2
            spgm(i_out)=(gmh(1)*rrs+gmh(intmax)*rrs1)*dtheta2

            do k=2,intmax-1
               rrs=rae(k)*rae(k)*sint(k)
               spgp(i_out)=spgp(i_out)+gp(k)*rrs*dtheta
               spgm(i_out)=spgm(i_out)+gmh(k)*rrs*dtheta
            enddo
            spgp(i_out)=spgp(i_out)*pi4
            spgm(i_out)=spgm(i_out)*pi4

            !  Find fp and ft

            fp(i_out) =  pi4 * r_out*r_out*r_out*r_out     / (cgrav(i_out)*xm_out*spgm(i_out))
            ft(i_out) = pow2(pi4 * r_out*r_out) / (spgp(i_out)*spgm(i_out))

            if (ft(i_out) < s% ft_error_limit) then
               ierr = -1
               if (report_ierr .or. dbg) then
                  write(*,2) 'FT too small', i_out, ft(i_out), s% ft_error_limit
                  stop 'eval_fp_ft'
               end if
               exit
               if (ift_in == 0) ift_in=i_out
               ift_out=i_out
               ft_min=min(ft_min,ft(i_out))
            elseif (ift_in /= 0) then
               ift_in=0
               ift_out=0
               ft_min=1.0d0
            endif

            if (fp(i_out) < s% fp_error_limit) then
               ierr = -1
               if (report_ierr .or. dbg) then
                  write(*,2) 'FP too small', i_out, fp(i_out), s% fp_error_limit
                  stop 'eval_fp_ft'
               end if
               exit
               if (ifp_in == 0) ifp_in=i_out
               ifp_out=i_out
               fp_min=min(fp_min,fp(i_out))
            elseif (ifp_in /= 0) then
               ifp_in=0
               ifp_out=0
               fp_min=1.0d0
            endif

            ft(i_out)=max(s% ft_min, min(ft(i_out),1.d0))
            fp(i_out)=max(s% fp_min, min(fp(i_out),1.d0))

         end do

         deallocate(work, iwork)

         ft(1)=ft(2)
         fp(1)=fp(2)
         ft(nz)=ft(nz-1)
         fp(nz)=fp(nz-1)

         contains

         real(dp) function apot_fcn(rnu, df, lrpar, rpar, lipar, ipar, ierr) result(f)
            ! returns with ierr = 0 if was able to evaluate f and df/dx at x
            ! if df/dx not available, it is okay to set it to 0
            integer, intent(in) :: lrpar, lipar
            real(dp), intent(in) :: rnu
            real(dp), intent(out) :: df
            integer, intent(inout), pointer :: ipar(:) ! (lipar)
            real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
            integer, intent(out) :: ierr
            real(dp) :: wrnu2, a2, a3, dadrn, c1, rfk, drfk
            real(dp) :: etanu,wr,r,xm,a,r_psi
            real(dp), parameter :: c0=1.0d0/3.0d0
            real(dp), parameter :: c2=2.0d0/35.0d0
            real(dp), parameter :: c3=3.0d0*c2
            include 'formats'
            ierr = 0
            c1 = 0.6d0*cgrav(1)
            etanu = rpar(1)
            wr = rpar(2)
            r = rpar(3)
            xm = rpar(4)
            wrnu2=wr*wr*rnu*rnu/(c1*xm*(2.d0+etanu)) ! (1/3)*da/dr0
            a=rnu*wrnu2 ! E&S, A10  NOTE: typo in paper where it gives r0^2 instead of r0^3.
            a2=a*a
            a3=a2*a
            drfk= pow(max(1.0d-20,1.d0+0.6d0*a2-c2*a3),c0) ! dr_psi/dr0|A, E&S A13
            r_psi=rnu*drfk ! r_psi(r0)
            rfk=r_psi - r ! r_psi(r0) - true_r_psi
            dadrn=wrnu2*a ! (1/3)*a*da/dr0
            drfk=drfk+rnu/(drfk*drfk)*(1.2d0*dadrn-c3*a*dadrn) ! dr_psi/dr0
            ! divide by r to normalize
            f = rfk / r
            df = drfk / r
            if (rfk > 1.d30) f=1.d30
            if (df > 1.d30) df=1.d30
            rpar(5) = a
         end function apot_fcn

      end subroutine eval_fp_ft


      subroutine rad_fcn(n, x, h, y, dydx, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: n, lrpar, lipar
         real(dp), intent(in) :: x, h
         real(dp), intent(inout) :: y(:) ! (n)
         real(dp), intent(inout) :: dydx(:) ! (n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr ! nonzero means retry with smaller timestep.
         real(dp) :: rov
         ierr = 0
         rov = rpar(1)
         dydx(1)=(6.d0-6.d0*rov*(y(1)+1.d0)-y(1)*(y(1)-1.d0))
      end subroutine rad_fcn



      real(dp) function psiint(rho, xm, aw, eta, r0_in, r0_out)
         real(dp), intent(in) :: rho, xm, aw, eta, r0_in, r0_out
         psiint=rho*(pow8(r0_out)-pow8(r0_in))/xm
         psiint=psiint*(5.d0+eta)/(2.d0+eta)*0.125d0*aw*aw
      end function psiint


      real(dp) function gpsi(cgrav, rnu, tegra, kr, a, wr, xm_out, rae)
         ! rnu = r0; tegra = integral in A9; kr is the theta index
         real(dp), intent(in) :: cgrav(:), rnu, tegra, a, wr, xm_out
         real(dp), intent(inout) :: rae(:)
         integer, intent(in) :: kr
         real(dp) :: r, ri, ri2, wr2, dpdr, dpdthe
         r= rnu * (1.d0 - a*utheta(kr))
         ri=1.0d0/r
         ri2=ri*ri
         wr2=wr*wr
         rae(kr)=r
         dpdr=-cgrav(kr)*xm_out*ri2 +pi4*ri2*ri2*utheta(kr)*tegra +wr2*r*sint2(kr)
            ! dpsi/dr from A9
         dpdthe= (pi4*ri2*ri*tegra + (wr2*r*r))*cost(kr)*sint(kr)
            ! dpsi/dtheta from A9
         gpsi=dsqrt(dpdr*dpdr+(dpdthe*dpdthe*ri2)) !   E&S, A14
      end function gpsi

      subroutine compute_j_fluxes_and_extra_jdot(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         integer :: k
         real(dp) :: pi2_div4, part1
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
            part1 = -pi2_div4*pow4(s% r(k)+s% r(k+1))*pow2(s% rho(k))*(s% i_rot(k)+s% i_rot(k+1))*(s% am_nu_omega(k)+s% am_nu_omega(k+1))
            s% j_flux(k) = part1*(s% omega(k)-s% omega(k+1))/s% dm(k)

            s% dj_flux_dw00(k) = s% j_flux(k)/(s% i_rot(k)+s% i_rot(k+1))*s% di_rot_dw_div_wc(k) &
               -part1*s% j_rot(k)/pow2(s% i_rot(k))*s% di_rot_dw_div_wc(k)/s% dm(k)
            s% dj_flux_dwp1(k) = s% j_flux(k)/(s% i_rot(k)+s% i_rot(k+1))*s% di_rot_dw_div_wc(k+1) &
               +part1*s% j_rot(k+1)/pow2(s% i_rot(k+1))*s% di_rot_dw_div_wc(k+1)/s% dm(k)

            s% dj_flux_dj00(k) = part1/s% i_rot(k)/s% dm(k)
            s% dj_flux_djp1(k) = -part1/s% i_rot(k+1)/s% dm(k)

            s% dj_flux_dlnr00(k) = 4d0*s% j_flux(k)/(s% r(k)+s% r(k+1))*s% r(k) &
               + s% j_flux(k)/(s% i_rot(k)+s% i_rot(k+1))*s% di_rot_dlnr(k) &
               -part1*s% j_rot(k)/pow2(s% i_rot(k))*s% di_rot_dlnr(k)/s% dm(k)
            s% dj_flux_dlnrp1(k) = 4d0*s% j_flux(k)/(s% r(k)+s% r(k+1))*s% r(k+1) &
               + s% j_flux(k)/(s% i_rot(k)+s% i_rot(k+1))*s% di_rot_dlnr(k+1) &
               +part1*s% j_rot(k+1)/pow2(s% i_rot(k+1))*s% di_rot_dlnr(k+1)/s% dm(k)

            s% dj_flux_dlnd(k) = 2d0*s% j_flux(k)
         end do
         s% j_flux(s% nz) = 0d0
         s% dj_flux_dw00(s% nz) = 0d0
         s% dj_flux_dwp1(s% nz) = 0d0
         s% dj_flux_dj00(s% nz) = 0d0
         s% dj_flux_dlnr00(s% nz) = 0d0
         s% dj_flux_dlnrp1(s% nz) = 0d0
         s% dj_flux_dlnd(s% nz) = 0d0

         !! this is to test partials
         !s% solver_test_partials_val = s% j_flux(19)
         !s% solver_test_partials_var = s% i_j_rot
         !s% solver_test_partials_dval_dx =  s% dj_flux_dj00(19)
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

