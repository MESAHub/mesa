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


      module rotation_mix_info

      use const_def
      use num_lib
      use utils_lib
      use star_private_def
      use magnetic_diffusion

      implicit none

      private
      public :: set_rotation_mixing_info, smooth_for_rotation

      real(dp), parameter :: Ri_crit = 0.25d0 ! critical Richardson number
      real(dp), parameter :: R_crit = 2500d0 ! critical Reynolds number

      integer, parameter :: i_DSI = 1
      integer, parameter :: i_SH = i_DSI + 1
      integer, parameter :: i_SSI = i_SH + 1
      integer, parameter :: i_ES = i_SSI + 1
      integer, parameter :: i_GSF = i_ES + 1
      integer, parameter :: i_ST = i_GSF + 1
      integer, parameter :: num_instabilities = i_ST


      contains


      subroutine set_rotation_mixing_info(s, ierr)
         use star_utils, only: weighed_smoothing

         type (star_info), pointer :: s
         integer, intent(out) :: ierr

         real(dp) :: f_mu, q, D_ST_old, nu_ST_old
         ! the following are all defined at cell boundaries
         real(dp), dimension(:), pointer :: & ! just copies of pointers
            r, m, L, j_rot, gradT, grada, grav, visc, Ri
         real(dp), dimension(:), allocatable :: & ! allocated temporary storage
            csound, rho, T, P, cp, cv, chiRho, abar, zbar, gradT_sub_grada, &
            opacity, kap_cond, gamma1, mu_alt, eps_nuc, enu, L_neu, delta, &
            scale_height, omega, cell_dr, &
            dRho, dr, dPressure, domega, d_mu, d_j_rot, &
            dRho_dr, dRho_dr_ad, dr2omega, H_T, &
            domega_dlnR, Hj, dlnR_domega, &
            t_dyn, t_kh, &
            Ri_mu, Ri_T, &
            ve0, ve_mu, &
            v_ssi, h_ssi, Ris_1, Ris_2, &
            v_es, H_es, &
            v_gsf, H_gsf, &
            N2, N2_mu
         real(dp), allocatable :: smooth_work(:,:), saved(:,:)
         logical, allocatable :: unstable(:,:) ! (num_instabilities, nz)

         integer :: nz, i, j, k, k0, which, op_err
         real(dp) :: alfa, beta, growth_limit, age_fraction, &
            D_omega_source, max_replacement_factor, &
            dgtau, angsml, angsmt, out_q, prev_out_q, prev_out_q_m1, prev_q, prev_q_m1

         include 'formats'

         ierr = 0
         nz = s% nz

         s% D_visc(1:nz) = 0
         s% D_DSI(1:nz) = 0
         s% D_SH(1:nz) = 0
         s% D_SSI(1:nz) = 0
         s% D_ES(1:nz) = 0
         s% D_GSF(1:nz) = 0
         s% D_ST(1:nz) = 0
         s% nu_ST(1:nz) = 0
         s% dynamo_B_r(1:nz) = 0
         s% dynamo_B_phi(1:nz) = 0
         s% omega_shear(1:nz) = 0
         s% D_omega(1:nz) = 0

         if (all(s% omega(1:nz) == 0d0)) return

         call setup(ierr)
         if (failed('setup for set_rotation_mixing_info', ierr)) return

         unstable(:,1:nz) = .false.
         growth_limit = 1d-10

!$OMP PARALLEL DO PRIVATE(which, k, q, age_fraction, op_err) SCHEDULE(dynamic,2)
         do which = 1, num_instabilities

            if (ierr /= 0) cycle
            op_err = 0

            select case (which)

               case (i_DSI)

                  if (s% D_DSI_factor > 0 .or. s% am_nu_DSI_factor > 0) then
                     call set_D_DSI(op_err)
                     if (failed('set_D_DSI', op_err)) then
                        ierr = -1; cycle
                     end if
                     call smooth_for_rotation(s, s% D_DSI, s% smooth_D_DSI, smooth_work(1:nz,which))
                     if (s% skip_rotation_in_convection_zones) &
                        call zero_if_convective(nz, s% mixing_type, s% D_mix, s% D_DSI)
                     call zero_if_tiny(s,s% D_DSI)
                  end if

               case (i_SH)

                  if (s% D_SH_factor > 0 .or. s% am_nu_SH_factor > 0) then
                     call set_D_SH(op_err)
                     if (failed('set_D_SH', op_err)) then
                        ierr = -1; cycle
                     end if
                     call smooth_for_rotation(s, s% D_SH, s% smooth_D_SH, smooth_work(1:nz,which))
                     if (s% skip_rotation_in_convection_zones) &
                        call zero_if_convective(nz, s% mixing_type, s% D_mix, s% D_SH)
                     call zero_if_tiny(s,s% D_SH)
                  end if

               case (i_SSI)

                  if (s% D_SSI_factor > 0 .or. s% am_nu_SSI_factor > 0) then
                     call set_D_SSI(op_err)
                     if (failed('set_D_SSI', op_err)) then
                        ierr = -1; cycle
                     end if
                     call smooth_for_rotation(s, s% D_SSI, s% smooth_D_SSI, smooth_work(1:nz,which))
                     if (s% skip_rotation_in_convection_zones) &
                        call zero_if_convective(nz, s% mixing_type, s% D_mix, s% D_SSI)
                     call zero_if_tiny(s,s% D_SSI)
                  end if

               case (i_ES)

                  if (s% D_ES_factor > 0 .or. s% am_nu_ES_factor > 0) then
                     call set_D_ES(op_err)
                     if (failed('set_D_ES', op_err)) then
                        ierr = -1; cycle
                     end if
                     call smooth_for_rotation(s, s% D_ES, s% smooth_D_ES, smooth_work(1:nz,which))
                     if (s% skip_rotation_in_convection_zones) &
                        call zero_if_convective(nz, s% mixing_type, s% D_mix, s% D_ES)
                     call zero_if_tiny(s,s% D_ES)
                  end if

               case (i_GSF)

                  if (s% D_GSF_factor > 0 .or. s% am_nu_GSF_factor > 0) then
                     call set_D_GSF(op_err)
                     if (failed('set_D_GSF', op_err)) then
                        ierr = -1; cycle
                     end if
                     call smooth_for_rotation(s, s% D_GSF, s% smooth_D_GSF, smooth_work(1:nz,which))
                     if (s% skip_rotation_in_convection_zones) &
                        call zero_if_convective(nz, s% mixing_type, s% D_mix, s% D_GSF)
                     call zero_if_tiny(s,s% D_GSF)
                  end if

               case (i_ST)

                  if (s% D_ST_factor > 0 .or. s% am_nu_ST_factor > 0) then

                     call set_ST(s, &
                        rho, T, r, L, omega, Cp, abar, zbar, delta, grav, &
                        N2, N2_mu, opacity, kap_cond, scale_height, op_err)
                     if (failed('set_ST', op_err)) then
                        ierr = -1; cycle
                     end if

                     call smooth_for_rotation(s, s% D_ST, s% smooth_D_ST, smooth_work(1:nz,which))
                     call smooth_for_rotation(s, s% nu_ST, s% smooth_nu_ST, smooth_work(1:nz,which))

                     ! time smooth for ST only
                     if (s% prev_mesh_have_ST_start_info .and. .not. s% doing_relax &
                        .and. .not. s% doing_first_model_of_run .and. s% ST_angsmt>0d0) then

                        k0 = 1
                        angsml = s% ST_angsml
                        angsmt = s% ST_angsmt
                        prev_out_q = 0d0
                        prev_q = 1d0
                        do k=2, nz-1 ! ignore k=1 or nz edge case
                           if (s% m(k) > s% mstar_old) cycle
                           do while (k0 < s% prev_mesh_nz)
                              if (s% m(k)/s% mstar_old > prev_q) exit
                              prev_out_q = prev_out_q + s% prev_mesh_dq(k0)
                              prev_q = 1d0 - prev_out_q
                              k0 = k0+1
                           end do
                           if (s% m(k)/s% mstar_old < prev_q) exit
                           prev_out_q_m1 = prev_out_q - s% prev_mesh_dq(k0-1)
                           prev_q_m1 = 1 - prev_out_q_m1

                           ! linear interpolation
                           alfa = (s% m(k)/s% mstar_old - prev_q)/(prev_q_m1-prev_q)
                           D_ST_old = (1d0-alfa)*s% prev_mesh_D_ST_start(k0) + alfa*s% prev_mesh_D_ST_start(k0-1)
                           nu_ST_old = (1d0-alfa)*s% prev_mesh_nu_ST_start(k0) + alfa*s% prev_mesh_nu_ST_start(k0-1)
                           dgtau = angsml*(s% r(k)-s% r(k+1))*(s% r(k-1)*s% r(k))/s% dt
                           s% D_ST(k) = max(D_ST_old/(1d0 + angsmt), &
                              min(s% D_ST(k), max(D_ST_old*(1d0 + angsmt), D_ST_old + dgtau)))
                           s% nu_ST(k) = max(nu_ST_old/(1d0 + angsmt), &
                              min(s% nu_ST(k), max(nu_ST_old*(1d0 + angsmt), nu_ST_old + dgtau)))
                        end do
                     end if

                     ! calculate B_r and B_phi
                     do k = 1, nz
                        q = s% omega_shear(k)
                        s% dynamo_B_r(k) = & ! eqn 11, Heger et al. 2005
                           pow(pow2(pi4*rho(k)*s% nu_ST(k)*q/r(k))*abs(omega(k))*s% nu_ST(k),0.25D0)
                        s% dynamo_B_phi(k) = & ! eqn 12, Heger et al. 2005
                           pow(pow2(pi4*rho(k)*omega(k)*q*r(k))*abs(omega(k))*s% nu_ST(k),0.25d0)
                     end do

                     if (s% skip_rotation_in_convection_zones) &
                        call zero_if_convective(nz, s% mixing_type, s% D_mix, s% D_ST)
                     if (s% skip_rotation_in_convection_zones) &
                        call zero_if_convective(nz, s% mixing_type, s% D_mix, s% nu_ST)
                     if (s% skip_rotation_in_convection_zones) &
                        call zero_if_convective(nz, s% mixing_type, s% D_mix, s% dynamo_B_r)
                     if (s% skip_rotation_in_convection_zones) &
                        call zero_if_convective(nz, s% mixing_type, s% D_mix, s% dynamo_B_phi)
                     call zero_if_tiny(s,s% D_ST)
                     call zero_if_tiny(s,s% nu_ST)

                  end if

               case default
                  call mesa_error(__FILE__,__LINE__,'bad case for rotation_mix_info')

            end select

         end do
!$OMP END PARALLEL DO
         if (failed('set_rotation_mixing_info instabilities', ierr)) return
         
         if (s% D_omega_flag .and. s% doing_finish_load_model) then
            do k=1,nz
               s% D_omega(k) = 0d0
            end do
         else if (s% D_omega_flag) then
                     
            do k=1,nz
               if (s% q(k) <= s% max_q_for_D_omega_zero_in_convection_region .and. &
                   s% mixing_type(k) == convective_mixing) then
                  s% D_omega(k) = 0d0
                  cycle
               end if
               D_omega_source = &
                  s% D_DSI_factor  * s% D_DSI(k)  + &
                  s% D_SH_factor   * s% D_SH(k)   + &
                  s% D_SSI_factor  * s% D_SSI(k)  + &
                  s% D_ES_factor   * s% D_ES(k)   + &
                  s% D_GSF_factor  * s% D_GSF(k)  + &
                  s% D_ST_factor   * s% D_ST(k)
               if (is_bad(D_omega_source)) then
                  write(*,2) 'D_omega_source', k, D_omega_source
                  write(*,2) 's% D_DSI_factor  * s% D_DSI(k)', k, s% D_DSI_factor  * s% D_DSI(k), &
                     s% D_DSI_factor, s% D_DSI(k)
                  write(*,2) 's% D_SH_factor   * s% D_SH(k)', k, s% D_SH_factor   * s% D_SH(k), &
                     s% D_SH_factor, s% D_SH(k)
                  write(*,2) 's% D_SSI_factor  * s% D_SSI(k)', k, s% D_SSI_factor  * s% D_SSI(k), &
                     s% D_SSI_factor, s% D_SSI(k)
                  write(*,2) 's% D_ES_factor   * s% D_ES(k)', k, s% D_ES_factor   * s% D_ES(k), &
                     s% D_ES_factor, s% D_ES(k)
                  write(*,2) 's% D_GSF_factor  * s% D_GSF(k)', k, s% D_GSF_factor  * s% D_GSF(k), &
                     s% D_GSF_factor, s% D_GSF(k)
                  write(*,2) 's% D_ST_factor   * s% D_ST(k)', k, s% D_ST_factor   * s% D_ST(k), &
                     s% D_ST_factor, s% D_ST(k)
                  call mesa_error(__FILE__,__LINE__,'rotation mix')
               end if
               s% D_omega(k) = D_omega_source
               if (is_bad(s% D_omega(k))) then
                  write(*,2) 's% D_omega(k)', k, s% D_omega(k)
                  write(*,2) 'D_omega_source', k, D_omega_source
                  call mesa_error(__FILE__,__LINE__,'rotation mix')
               end if
            end do
            
            if (s% smooth_D_omega > 0) then
               call smooth_for_rotation(s, s% D_omega, s% smooth_D_omega, smooth_work(1:nz,1))
               do k=1,nz
                  if (is_bad(s% D_omega(k))) then
                     write(*,2) 'after smooth_for_rotation s% D_omega(k)', k, s% D_omega(k)
                     call mesa_error(__FILE__,__LINE__,'rotation mix')
                  end if
               end do
            end if
            
            if (s% D_omega_mixing_rate > 0d0 .and. s% dt > 0) then
               call mix_D_omega 
               do k=1,nz
                  if (is_bad(s% D_omega(k))) then
                     write(*,2) 'after mix_D_omega s% D_omega(k)', k, s% D_omega(k)
                     call mesa_error(__FILE__,__LINE__,'rotation mix')
                  end if
               end do
            end if
            
         end if
         
         if (s% D_omega_flag) then
            do k=1,nz
               if (is_bad(s% D_omega(k))) then
                  write(*,2) 'before return s% D_omega(k)', k, s% D_omega(k)
                  call mesa_error(__FILE__,__LINE__,'rotation mix')
               end if
               if (s% D_omega(k) < 0d0) s% D_omega(k) = 0d0
            end do
         end if


         contains
         
         
         subroutine mix_D_omega
            integer :: i, k, nz
            real(dp), dimension(:), allocatable :: & ! work vectors
               sig, rhs, d, du, dl, bp, vp, xp, x
            real(dp) :: &
               dt, rate, d_ddt_dm1, d_ddt_d00, d_ddt_dp1, m, &
               d_dt, d_dt_in, d_dt_out
            include 'formats'
            
            nz = s% nz
            dt = s% dt
            if (dt == 0) return
            
            allocate(sig(nz), rhs(nz), d(nz), du(nz), dl(nz), bp(nz), vp(nz), xp(nz), x(nz))
            
            rate = min(s% D_omega_mixing_rate, 1d0/dt)
            do k=2,nz-1
               if (s% D_omega(k) == 0 .or. s% D_omega(k+1) == 0) then
                  sig(k) = 0
               else if ((.not. s% D_omega_mixing_across_convection_boundary) .and. &
                  s% mixing_type(k) /= convective_mixing .and. &
                      (s% mixing_type(k-1) == convective_mixing .or. &
                       s% mixing_type(k+1) == convective_mixing)) then
                   sig(k) = 0
               else
                  sig(k) = rate*dt
               end if               
            end do
            sig(1) = 0
            sig(nz) = 0
            
            do k=1,nz
               if (k < nz) then
                  d_dt_in = sig(k)*(s% D_omega(k+1) - s% D_omega(k))
               else
                  d_dt_in = -sig(k)*s% D_omega(k)
               end if
               if (k > 1) then
                  d_dt_out = sig(k-1)*(s% D_omega(k) - s% D_omega(k-1))
                  d_ddt_dm1 = sig(k-1)
                  d_ddt_d00 = -(sig(k-1) + sig(k))
               else
                  d_dt_out = 0
                  d_ddt_dm1 = 0
                  d_ddt_d00 = -sig(k)
               end if
               d_dt = d_dt_in - d_dt_out
               d_ddt_dp1 = sig(k)
               rhs(k) = d_dt
               d(k) = 1d0 - d_ddt_d00
               if (k < nz) then
                  du(k) = -d_ddt_dp1
               else
                  du(k) = 0
               end if
               if (k > 1) dl(k-1) = -d_ddt_dm1               
            end do
            dl(nz) = 0
            
            ! solve tridiagonal
            bp(1) = d(1)
            vp(1) = rhs(1)
            do i = 2,nz
               m = dl(i-1)/bp(i-1)
               bp(i) = d(i) - m*du(i-1)
               vp(i) = rhs(i) - m*vp(i-1)
            end do
            xp(nz) = vp(nz)/bp(nz)
            x(nz) = xp(nz)
            do i = nz-1, 1, -1
               xp(i) = (vp(i) - du(i)*xp(i+1))/bp(i)
               x(i) = xp(i)
            end do
            
            do k=2,nz
               if (is_bad(x(k))) then
                  return
                  write(*,3) 's% D_omega(k) prev, x', k, s% model_number, s% D_omega(k), x(k), bp(i)
                  call mesa_error(__FILE__,__LINE__,'mix_D_omega')
               end if
            end do
            
            ! update D_omega
            
            do k=2,nz
               s% D_omega(k) = s% D_omega(k) + x(k)
               if (is_bad(s% D_omega(k))) then
                  write(*,3) 's% D_omega(k)', k, s% model_number, s% D_omega(k)
                  call mesa_error(__FILE__,__LINE__,'mix_D_omega')
               end if
               if (s% D_omega(k) < 0d0) s% D_omega(k) = 0d0
            end do
            s% D_omega(1) = 0d0
         
         end subroutine mix_D_omega


         subroutine setup(ierr)
            use kap_lib, only: kap_get_elect_cond_opacity
            integer, intent(out) :: ierr

            real(dp) :: &
               bracket_term, ri0, alfa, beta, enum1, enu00, &
               rho6, gamma, mu_e, rm23, ctmp, xi2, dynvisc, denom, &
               eps_nucm1, eps_nuc00, scale_height2, dlnRho_dlnP, dlnT_dlnP, &
               kap, dlnkap_dlnRho, dlnkap_dlnT

            integer :: i, k, j

            include 'formats'

            ierr = 0
            nz = s% nz

            f_mu = s% am_gradmu_factor

            ! copy some pointers
            r => s% r
            m => s% m
            L => s% L
            j_rot => s% j_rot
            gradT => s% gradT
            grada => s% grada_face
            grav => s% grav
            visc => s% D_visc
            Ri => s% richardson_number
            
            allocate( &
               csound(nz), rho(nz), T(nz), P(nz), cp(nz), cv(nz), chiRho(nz), abar(nz), zbar(nz), &
               opacity(nz), kap_cond(nz), gamma1(nz), mu_alt(nz), omega(nz), cell_dr(nz), eps_nuc(nz), enu(nz), L_neu(nz), &
               gradT_sub_grada(nz), delta(nz), scale_height(nz), &
               dRho(nz), dr(nz), dPressure(nz), domega(nz), d_mu(nz), d_j_rot(nz), &
               dRho_dr(nz), dRho_dr_ad(nz), dr2omega(nz), H_T(nz), &
               domega_dlnR(nz), Hj(nz), dlnR_domega(nz), t_dyn(nz), t_kh(nz), &
               Ri_mu(nz), Ri_T(nz), ve0(nz), ve_mu(nz), &
               v_ssi(nz), h_ssi(nz), Ris_1(nz), Ris_2(nz), &
               v_es(nz), H_es(nz), v_gsf(nz), H_gsf(nz), N2(nz), N2_mu(nz), &
               smooth_work(nz,num_instabilities), saved(nz,num_instabilities), &
               unstable(num_instabilities,nz))

            ! interpolate by mass to get values at cell boundaries
            enu00 = s% eps_nuc_neu_total(1) + s% non_nuc_neu(1)
            enu(1) = enu00
            eps_nuc00 = s% eps_nuc(1)
            eps_nuc(1) = eps_nuc00
            csound(1) = s% csound(1)
            rho(1) = s% rho(1)
            T(1) = s% T(1)
            P(1) = s% Peos(1)
            cp(1) = s% cp(1)
            cv(1) = s% Cv(1)
            chiRho(1) = s% chiRho(1)
            abar(1) = s% abar(1)
            zbar(1) = s% zbar(1)
            opacity(1) = s% opacity(1)
            gamma1(1) = s% gamma1(1)
            mu_alt(1) = s% abar(1)/(1 + s% zbar(1))
            delta(1) = s% chiT(1)/s% chiRho(1)
            L_neu(1) = enu00*s% dm(1)
            do k = 2, nz
               alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
               beta = 1 - alfa
               enum1 = enu00
               enu00 = s% eps_nuc_neu_total(k) + s% non_nuc_neu(k)
               enu(k) = alfa*enu00 + beta*enum1
               eps_nucm1 = eps_nuc00
               eps_nuc00 = s% eps_nuc(k)
               eps_nuc(k) = alfa*eps_nuc00 + beta*eps_nucm1
               csound(k) = alfa*s% csound(k) + beta*s% csound(k-1)
               rho(k) = alfa*s% rho(k) + beta*s% rho(k-1)
               T(k) = alfa*s% T(k) + beta*s% T(k-1)
               P(k) = alfa*s% Peos(k) + beta*s% Peos(k-1)
               cp(k) = alfa*s% cp(k) + beta*s% cp(k-1)
               cv(k) = alfa*s% Cv(k) + beta*s% Cv(k-1)
               chiRho(k) = alfa*s% chiRho(k) + beta*s% chiRho(k-1)
               abar(k) = alfa*s% abar(k) + beta*s% abar(k-1)
               zbar(k) = alfa*s% zbar(k) + beta*s% zbar(k-1)
               opacity(k) = alfa*s% opacity(k) + beta*s% opacity(k-1)
               gamma1(k) = alfa*s% gamma1(k) + beta*s% gamma1(k-1)
               mu_alt(k) = alfa*s% abar(k)/(1 + s% zbar(k)) + beta*s% abar(k-1)/(1 + s% zbar(k-1))
               delta(k) = alfa*s% chiT(k)/s% chiRho(k) + beta*s% chiT(k-1)/s% chiRho(k-1)
               L_neu(k) = enu00*s% dm(k) + L_neu(k-1)
               cell_dr(k-1) = s% rmid(k-1) - s% rmid(k)
            end do
            cell_dr(nz) = s% rmid(nz) - s% R_center

            do i = 1, nz
               gradT_sub_grada(i) = s% gradT(i) - s% grada_face(i)
               gradT_sub_grada(i) = & ! make sure it isn't too close to 0
                  sign(max(abs(gradT_sub_grada(i)),1d-99),gradT_sub_grada(i))
               scale_height(i) = P(i)*r(i)*r(i)/(s% cgrav(i)*m(i)*rho(i))
               scale_height2 = sqrt(P(i)/s% cgrav(i))/rho(i)
               if (scale_height2 < scale_height(i)) scale_height(i) = scale_height2
               omega(i) = s% omega(i)
            end do

            ! differences (at cell boundaries)
            do i = 2, nz-1
               dRho(i) = rho(i-1) - rho(i)
               dr(i) = max(1d0, 0.5D0*(r(i-1) - r(i+1)))
               dPressure(i) = min(-1d-10, P(i-1) - P(i))
               d_mu(i) = mu_alt(i-1) - mu_alt(i)
               d_j_rot(i) = j_rot(i-1) - j_rot(i)
               domega(i) = 0.5D0*(omega(i-1) - omega(i+1))
            end do

            dRho(1) = 0
            dRho(nz) = 0

            dr(1) = max(1d0,r(1)-r(2))
            dr(nz) = r(nz) - s% R_center

            dPressure(1) = 0
            dPressure(nz) = 0

            d_mu(1) = 0
            d_mu(nz) = 0

            d_j_rot(1) = 0
            d_j_rot(nz) = 0

            domega(1) = 0
            domega(nz) = 0

            do i = 2, nz-1
               dRho_dr(i) = dRho(i)/dr(i)
               dRho_dr_ad(i) = rho(i)*dPressure(i)/(P(i)*gamma1(i)*dr(i))
               dr2omega(i) = 4.5d0*j_rot(i)*d_j_rot(i)/dr(i) ! d(r^2 omega)^2/dr using i = (2/3)*r^2
               domega_dlnR(i) = domega(i)*r(i)/dr(i)
               if (gradT(i) > 1d-20) then
                  H_T(i) = scale_height(i)/gradT(i) ! -dr/dlnT, scale height of temperature
               else
                  H_T(i) = scale_height(i)
               endif
               if(abs(d_j_rot(i))>0.d0)then
                  Hj(i) = j_rot(i)*dr(i)/d_j_rot(i)
               else
                  Hj(i) =0d0
               end if
                  ! dr/dlnj, scale height of angular momentum
            end do

            dRho_dr(1) = 0; dRho_dr(nz) = 0
            dRho_dr_ad(1) = 0; dRho_dr_ad(nz) = 0
            dr2omega(1) = 0; dr2omega(nz) = 0
            domega_dlnR(1) = 0; domega_dlnR(nz) = 0
            H_T(1) = H_T(2); H_T(nz) = H_T(nz-1)
            Hj(1) = Hj(2); Hj(nz) = Hj(nz-1)

            do i = 1, nz
               dlnRho_dlnP = s% grad_density(i)
               dlnT_dlnP = s% grad_temperature(i)
               N2(i) = -grav(i)*(1/gamma1(i) - dlnRho_dlnP)/scale_height(i)
               N2_mu(i) = -grav(i)/scale_height(i)*(1/chiRho(i) - delta(i)*dlnT_dlnP - dlnRho_dlnP)
            end do

            do k=1,nz
               s% domega_dlnR(k) = domega_dlnR(k)
            end do

            ! safe inverse of domega/dlnR
            do i = 2, nz-1
               dlnR_domega(i) = sign(1d0/max(abs(domega_dlnR(i)),1d-30),domega_dlnR(i))
            end do
            dlnR_domega(1) = 0; dlnR_domega(nz) = 0

            do k = 2, nz - 1 ! shear
               s% omega_shear(k) = &
                  (omega(k-1) - omega(k+1))/(r(k-1) - r(k+1))*(r(k)/omega(k))
               s% omega_shear(k) = max(1d-30,min(1d30,abs(s% omega_shear(k))))
            end do
            s% omega_shear(1) = 0
            s% omega_shear(nz) = 0

            ! timescales
            do i = 1, nz
               t_dyn(i) = sqrt(r(i)*r(i)*r(i)/(s% cgrav(i)*m(i)))
               t_kh(i) = s% cgrav(i)*m(i)*m(i)/(r(i)*max(1d0,L(i)+L_neu(i)))
            end do

            ! Richardson numbers (Heger 2000, eqn 20)
            do i = 2, nz-1
               ri0 = (rho(i)*delta(i)/P(i))*pow2(dlnR_domega(i)*grav(i))
               Ri_T(i) = ri0*max(0d0,-gradT_sub_grada(i)) ! turn off Ri_T in convection zones
               Ri_mu(i) = ri0*f_mu*s% smoothed_brunt_B(i)
            end do
            Ri_T(1) = 0; Ri_T(nz) = 0
            Ri_mu(1) = 0; Ri_mu(nz) = 0
            do i=1,nz
               if (N2(i) < 0d0) then
                  Ri(i) = 1d0 ! disable in convective region
               else
                  Ri(i) = Ri_T(i) + Ri_mu(i)
               end if
            end do

            ! dynamic viscosity
            do i=1,nz
               rho6 = rho(i)*1d-6
               gamma = 0.2275d0*zbar(i)*zbar(i)*pow(rho6/abar(i),one_third)*1.d8/T(i)
                  ! gamma => eq (5) of Itoh et al 1987 ApJ 317,733
               ! electron viscosity according to Nandkumar & Pethick 1984 MNRAS
               mu_e = abar(i)/zbar(i)
               rm23 = pow(rho6/mu_e,2d0/3d0)
               ctmp = 1d0 + 1.018d0*rm23
               xi2 = sqrt(pi/3d0)*log(zbar(i))/3d0 + 2d0*log(1.32d0+2.33d0/sqrt(gamma)) - &
                     0.475d0*(1d0+2.036d0*rm23)/ctmp + 0.276d0*rm23/ctmp
               ! dynamic shear viscosity according to Wallenborn and Bauss (1978)
               !  and also Itoh et al 1987 ApJ 317,733
               ! fitting formula for eta* in Eq (12) of Itoh et al. 1987
               ctmp = -0.016321227d0+1.0198850d0*pow(gamma,-1.9217970d0) + &
                       0.024113535d0*pow(gamma,0.49999098d0)
               ! dynamic shear viscosity
               dynvisc = 5.53d3*zbar(i)*pow(rho6,5d0/6d0)*ctmp/pow(abar(i),one_third)
               ! add contibution of radiation
               dynvisc = dynvisc + 4.D0*crad*pow4(T(i))/(15.D0*clight*opacity(i)*rho(i))
               ! add contibution of electrons
               dynvisc = dynvisc + 1.893d6*pow(rm23,2.5d0)/(zbar(i)*ctmp*xi2)
               ! kinematic shear viscosity
               visc(i) = dynvisc/rho(i)
            end do

            ! velocities for ES and GSF
            if (s% D_ES_factor > 0 .or. s% D_GSF_factor > 0) then
               do i = 1, nz ! Heger 2000, eqns 35 and 36
                  ! the bracket_term blows up at center since r^2/L and r^2/m can both -> Inf
                  ! so bullet proof by including lower bounds
                  bracket_term = &
                     2*r(i)*r(i)*(eps_nuc(i)/max(1d-3*Lsun,abs(L(i))) - 1/max(1d-3*Msun,m(i))) - &
                     3/(pi4*rho(i)*max(1d-3*Rsun,r(i)))
                  if (abs(gradT_sub_grada(i)) < 1d-50) then
                     ve0(i) = 1d99
                     ve_mu(i) = 1d99
                  else
                     denom = (-gradT_sub_grada(i)*delta(i)*pow2(s% cgrav(i)*m(i)))
                     if (abs(denom) < 1d-50 .or. is_bad(denom)) then
                        if (s% stop_for_bad_nums) then
                           write(*,2) 'denom', i, denom
                           call mesa_error(__FILE__,__LINE__,'rotation mix info: velocities for ES and GSF')
                        end if
                        ve0(i) = 1d99
                     else
                        ve0(i) = grada(i)*omega(i)*omega(i)*r(i)*r(i)*r(i)*L(i)*bracket_term/denom
                     end if
                     ve_mu(i) = (scale_height(i)/t_kh(i))* &
                              (f_mu*s% smoothed_brunt_B(i))/(gradT_sub_grada(i))
                  end if
                  if (is_bad(ve0(i))) then
                     if (s% stop_for_bad_nums) then
                        write(*,2) 've0(i)', i, ve0(i)
                        call mesa_error(__FILE__,__LINE__,'rotation mix info')
                     end if
                     ve0(i) = 1d99
                  end if
                  if (is_bad(ve_mu(i))) then
                     if (s% stop_for_bad_nums) then
                        write(*,2) 've_mu(i)', i, ve_mu(i)
                        call mesa_error(__FILE__,__LINE__,'rotation mix info')
                     end if
                     ve_mu(i) = 1d99
                  end if

                  if (s% model_number == -1 .and. i == 316) then
                     write(*,2) 'grada(i)', i, grada(i)
                     write(*,2) 'gradT(i)', i, gradT(i)
                     write(*,2) 'gradT_sub_grada(i)', i, gradT_sub_grada(i)
                     write(*,2) 's% smoothed_brunt_B(i)', i, s% smoothed_brunt_B(i)
                     write(*,2) 'omega(i)', i, omega(i)
                     write(*,2) 's% omega(i)', i, s% omega(i)
                     write(*,2) '2*r**2*eps_nuc/L', i, 2*r(i)*r(i)*eps_nuc(i)/max(1d0,L(i))
                     write(*,2) '2*r**2/m', i, 2*r(i)*r(i)/m(i)
                     write(*,2) '3/(4*pi*rho*r)', i, 3/(4*pi*rho(i)*r(i))
                     write(*,2) 've0(i)', i, ve0(i)
                  end if

               end do

               ! conductive opacities for ST
               if (s% D_ST_factor > 0d0 .or. s% am_nu_factor > 0d0) then
                  do i = 1, nz
                     call kap_get_elect_cond_opacity( &
                        s% kap_handle, zbar(i), log10(rho(i)), log10(T(i)),  &
                        kap, dlnkap_dlnRho, dlnkap_dlnT, ierr)
                     if (ierr /= 0) return
                     kap_cond(i) = kap
                  end do
               end if

            end if

         end subroutine setup


         subroutine set_D_DSI(ierr)
            integer, intent(out) :: ierr
            integer :: i, k, kbot, ktop
            real(dp) :: instability_height, height, D
            logical, parameter :: dbg = .false.
            include 'formats'

            ierr = 0
            do i = 1, nz
               unstable(i_DSI,i) = (Ri(i) < Ri_crit) .and. (gradT_sub_grada(i) < 0)
               ! stable in convective regions where gradT >= grada
            end do

            kbot = nz
            do i = nz-1, 1, -1

               if (unstable(i_DSI,i) .and. .not. unstable(i_DSI,i+1)) kbot = i

               if (unstable(i_DSI,i+1) .and. &
                     (i == 1 .or. .not. unstable(i_DSI,i)) .and. kbot > 1) then

                  if (unstable(i_DSI,i)) then
                     ktop = i
                  else
                     ktop = i+1
                  end if

                  if (ktop >= kbot) cycle

                  instability_height = r(ktop) - r(kbot)
                  if (dbg) write(*,3) 'DSI: ktop, kbot', ktop, kbot
                  do k = ktop, kbot
                     if (.not. unstable(i_DSI,k)) then
                        write(*,2) 'D_DSI where stable?', k, s% q(k), &
                           s% D_DSI(k), D, scale_height(k)*csound(k), &
                           Ri(k), Ri_crit, height, t_dyn(k), &
                           instability_height, scale_height(k), csound(k)
                        call mesa_error(__FILE__,__LINE__,'set_D_DSI')
                     end if
                     height = min(instability_height, scale_height(k))
                     D = height*height/t_dyn(k)
                     s% D_DSI(k) = min(D, scale_height(k)*csound(k))
                     if (dbg) write(*,2) 'D_DSI', k, s% q(k), &
                        s% D_DSI(k), D, scale_height(k)*csound(k), &
                        Ri(k), Ri_crit, height, t_dyn(k), &
                        instability_height, scale_height(k), csound(k)
                  end do

                  if (dbg) write(*,*)

               end if

            end do
            if (dbg) call mesa_error(__FILE__,__LINE__,'set_D_DSI')
         end subroutine set_D_DSI


         subroutine set_D_SH(ierr) ! comment in Langer code says "DO NOT USE"
            integer, intent(out) :: ierr
            integer :: i, k, kbot, ktop
            real(dp) :: instability_height, height, D
            ierr = 0

            do i = 1, nz
               D = grav(i)/rho(i)*(dRho_dr_ad(i)-dRho_dr(i))+dr2omega(i)/(r(i)*r(i)*r(i))
               if (D < 0) then
                  unstable(i_SH,i) = .true.
                  s% D_SH(i) = D ! save for later
               else
                  s% D_SH(i) = 0
               end if
            end do

            kbot = nz
            do i = nz-1, 1, -1

               if (unstable(i_SH,i) .and. .not. unstable(i_SH,i+1)) kbot = i

               if (unstable(i_SH,i+1) .and. &
                     (i == 1 .or. .not. unstable(i_SH,i)) .and. kbot > 1) then
                  if (unstable(i_SH,i)) then
                     ktop = i
                  else
                     ktop = i+1
                  end if
                  if (ktop >= kbot) cycle
                  instability_height = r(ktop) - r(kbot)
                  do k = ktop, kbot
                     height = min(instability_height, scale_height(k))
                     ! use the previously calculated value saved in D_SH
                     D = pow2(height*s% D_SH(k)*r(k)/grav(k))/t_dyn(k)
                     s% D_SH(k) = min(D, scale_height(k)*csound(k))
                  end do
               end if

            end do
         end subroutine set_D_SH


         subroutine set_D_SSI(ierr)
            use chem_def
            integer, intent(out) :: ierr
            integer :: i, k, kbot, ktop
            real(dp) :: qe3, qe4, lambda, dynvisc, Prandtl, radcon, D
            include 'formats'

            ierr = 0
            qe3 = qe*qe*qe
            qe4 = qe3*qe

            do i=1,nz
               unstable(i_SSI,i) = .false.
               ! thermal conductivity
               radcon = 4.D0*crad*clight*T(i)*T(i)*T(i)/(3.D0*opacity(i)*rho(i)) ! erg / (K cm sec)
               ! Prandtl-number according to Tassoul
               dynvisc = visc(i)*rho(i)
               Prandtl = dynvisc*Cv(i)/radcon
               Ris_1(i) = 0.125D0*R_crit*Prandtl*Ri_T(i)
               if (Ris_1(i) <= Ri_crit) then
                  Ris_2(i) = Ri_mu(i)
                  if (Ris_2(i) <= Ri_crit) unstable(i_SSI,i) = .true.
               else
                  Ris_2(i) = 0
               end if
            end do

            kbot = nz
            do i = nz-1, 1, -1

               if (unstable(i_SSI,i) .and. .not. unstable(i_SSI,i+1)) kbot = i

               if (unstable(i_SSI,i+1) .and. &
                     (i == 1 .or. .not. unstable(i_SSI,i)) .and. kbot > 1) then
                  if (unstable(i_SSI,i)) then
                     ktop = i
                  else
                     ktop = i+1
                  end if

                  if (ktop >= kbot) cycle

                  do k = ktop, kbot ! Heger 2000, eqn 31
                     v_ssi(k)=sqrt(visc(k)/R_crit*abs(domega_dlnR(k)))
                  end do

                  H_ssi(kbot) = v_ssi(kbot)*(r(kbot-1) - r(kbot))/ &
                                 max(1d-99,abs(v_ssi(kbot-1) - v_ssi(kbot)))
                  do k = kbot-1, ktop+1, -1
                     H_ssi(k) = v_ssi(k)*dr(k)/ &
                                 max(1d-99,0.5d0*abs(v_ssi(k+1) - v_ssi(k-1)))
                  end do
                  H_ssi(ktop) = v_ssi(ktop)*(r(ktop) - r(ktop+1))/ &
                                 max(1d-99,abs(v_ssi(ktop) - v_ssi(ktop+1)))

                  do k = ktop, kbot ! Heger 2000, eqn 34
                     H_ssi(k) = min(H_ssi(k),scale_height(k))
                     v_ssi(k) = min(v_ssi(k),csound(k))
                     D = H_ssi(k)*v_ssi(k)* &
                           pow2(1d0-max(0d0,max(Ris_1(k),Ris_2(k))/Ri_crit))
                     s% D_SSI(k) = min(D, scale_height(k)*csound(k))

                     if (s% D_SSI(k) > 1d100) then ! bug
                        write(*,2) 's% D_SSI(k)', k, s% D_SSI(k)
                        write(*,2) 'D', k, D
                        write(*,2) 'scale_height(k)', k, scale_height(k)
                        write(*,2) 'csound(k)', k, csound(k)
                        write(*,2) 'H_ssi(k)', k, H_ssi(k)
                        write(*,2) 'v_ssi(k)', k, v_ssi(k)
                        write(*,2) 'Ris_1(k)', k, Ris_1(k)
                        write(*,2) 'Ris_2(k)', k, Ris_2(k)
                        write(*,2) 'dr(k)', k, dr(k)
                        write(*,2) 'visc(k)', k, visc(k)
                        write(*,2) 'domega_dlnR(k)', k, domega_dlnR(k)
                        call mesa_error(__FILE__,__LINE__,'set_D_SSI')
                     end if

                  end do

               end if

            end do

         end subroutine set_D_SSI


         subroutine set_D_ES(ierr)
            integer, intent(out) :: ierr
            integer :: i, k, kbot, ktop
            real(dp) :: instability_height, D, v, dln_v_es
            include 'formats'
            ierr = 0

            do i = 1, nz
               v = abs(ve0(i)) - abs(ve_mu(i)) ! heger 2000, eqn 38
               if (v > 0) then
                  unstable(i_ES,i) = .true.
                  v_es(i) = v
               else
                  v_es(i) = 0
               end if
            end do

            kbot = nz
            do i = nz-1, 1, -1

               if (unstable(i_ES,i) .and. .not. unstable(i_ES,i+1)) kbot = i

               if (unstable(i_ES,i+1) .and. &
                     (i == 1 .or. .not. unstable(i_ES,i)).and. kbot > 1) then
                  if (unstable(i_ES,i)) then
                     ktop = i
                  else
                     ktop = i+1
                  end if

                  if (ktop >= kbot) cycle

                  instability_height = r(ktop) - r(kbot)

                  ! heger 2000, eqn 39
                  H_es(kbot) = v_es(kbot)*(r(kbot-1) - r(kbot))/ &
                              max(1d-99,abs(v_es(kbot-1) - v_es(kbot)))
                  do k = kbot-1, ktop+1, -1
                     H_es(k) = v_es(k)*dr(k)/ &
                              max(1d-99,0.5d0*abs(v_es(k+1) - v_es(k-1)))
                  end do
                  H_es(ktop) = v_es(ktop)*(r(ktop) - r(ktop+1))/ &
                              max(1d-99,abs(v_es(ktop) - v_es(ktop+1)))

                  do k = ktop, kbot
                     H_es(k) = min(instability_height, H_es(k), scale_height(k))
                     v_es(k) = min(v_es(k), csound(k))
                     D = H_es(k)*v_es(k)
                     s% D_ES(k) = min(D, scale_height(k)*csound(k))
                  end do

               end if

            end do

         end subroutine set_D_ES


         subroutine set_D_GSF(ierr)
            integer, intent(out) :: ierr
            integer :: i, k, kbot, ktop
            real(dp) :: instability_height, D, v, v_diff
            include 'formats'
            ierr = 0

            do i = 1, nz

               ! heger 2000, eqn 42
               if(abs(Hj(i))>0d0) then
                  v = ve0(i)*2*H_T(i)*r(i)/(Hj(i)*Hj(i))/(1 + 2*omega(i)*dlnR_domega(i))
               else
                  v = 0d0
               end if
               if (is_bad(v)) then
                  write(*,2) 'bad v for GSF', i, v
                  write(*,2) 've0(i)', i, ve0(i)
                  write(*,2) 'H_T(i)', i, H_T(i)
                  write(*,2) 'r(i)', i, r(i)
                  write(*,2) 'Hj(i)', i, Hj(i)
                  write(*,2) 'omega(i)', i, omega(i)
                  write(*,2) 'dlnR_domega(i)', i, dlnR_domega(i)
                  call mesa_error(__FILE__,__LINE__,'set_D_GSF')
                  v = 0
               end if
               v_diff = abs(v) - abs(ve_mu(i)) ! heger 2000, eqn 43
               if (v_diff > 0) then
                  unstable(i_GSF,i) = .true.
                  v_gsf(i) = v_diff
               else
                  v_gsf(i) = 0
               end if

               if (s% model_number == -3 .and. i == -1) then
                  write(*,2) 've0(i)', i, ve0(i)
                  write(*,2) 'H_T(i)', i, H_T(i)
                  write(*,2) 'r(i)', i, r(i)
                  write(*,2) 'Hj(i)', i, Hj(i)
                  write(*,2) 'omega(i)', i, omega(i)
                  write(*,2) 'dlnR_domega(i)', i, dlnR_domega(i)
                  write(*,2) 'v', i, v
                  write(*,2) 've_mu(i)', i, ve_mu(i)
                  write(*,2) 'v_gsf(i)', i, v_gsf(i)
               end if

            end do

            kbot = nz
            do i = nz-1, 1, -1

               if (unstable(i_GSF,i) .and. .not. unstable(i_GSF,i+1)) kbot = i

               if (unstable(i_GSF,i+1) .and. &
                     (i == 1 .or. .not. unstable(i_GSF,i)) .and. kbot > 1) then
                  if (unstable(i_GSF,i)) then
                     ktop = i
                  else
                     ktop = i+1
                  end if

                  if (ktop >= kbot) cycle

                  instability_height = r(ktop) - r(kbot)

                  ! heger 2000, eqn 45
                  if (kbot == 1) then
                     H_gsf(kbot) = v_gsf(kbot)*(s% R_center - r(kbot))/ &
                              max(1d-99,abs(v_gsf(kbot)))
                  else
                     H_gsf(kbot) = v_gsf(kbot)*(r(kbot-1) - r(kbot))/ &
                              max(1d-99,abs(v_gsf(kbot-1) - v_gsf(kbot)))
                  end if

                  if (kbot == -1) then
                     write(*,2) 'r(kbot)', kbot, r(kbot)
                     write(*,2) 'r(kbot-1)', kbot-1, r(kbot-1)
                     write(*,2) 'v_gsf(kbot)', kbot, v_gsf(kbot)
                     write(*,2) 'v_gsf(kbot-1)', kbot-1, v_gsf(kbot-1)
                     write(*,2) 'H_gsf(kbot)', kbot, H_gsf(kbot)
                  end if

                  do k = max(2,kbot-1), min(nz-1,ktop+1), -1
                     H_gsf(k) = v_gsf(k)*dr(k)/ &
                              max(1d-99,0.5d0*abs(v_gsf(k+1) - v_gsf(k-1)))
                     if (k == -1) then
                        write(*,2) 'dr(k)', k, dr(k)
                        write(*,2) 'v_gsf(k-1)', k-1, v_gsf(k-1)
                        write(*,2) 'v_gsf(k)', k, v_gsf(k)
                        write(*,2) 'v_gsf(k+1)', k+1, v_gsf(k+1)
                        write(*,2) 'H_gsf(k)', k, H_gsf(k)
                     end if
                  end do

                  if (ktop == nz) then
                     H_gsf(ktop) = v_gsf(ktop)*(r(ktop) - s% R_center)/ &
                              max(1d-99,abs(v_gsf(ktop)))
                  else
                     H_gsf(ktop) = v_gsf(ktop)*(r(ktop) - r(ktop+1))/ &
                              max(1d-99,abs(v_gsf(ktop) - v_gsf(ktop+1)))
                  end if

                  if (ktop == -1) then
                     write(*,2) 'r(ktop)', ktop, r(ktop)
                     write(*,2) 'r(ktop+1)', ktop+1, r(ktop+1)
                     write(*,2) 'v_gsf(ktop)', ktop, v_gsf(ktop)
                     write(*,2) 'v_gsf(ktop+1)', ktop+1, v_gsf(ktop+1)
                     write(*,2) 'H_gsf(ktop)', ktop, H_gsf(ktop)
                  end if

                  do k = ktop, kbot
                     H_gsf(k) = min(instability_height, H_gsf(k), scale_height(k))
                     v_gsf(k) = min(v_gsf(k), csound(k))
                     D = H_gsf(k)*v_gsf(k)
                     s% D_GSF(k) = min(D, scale_height(k)*csound(k))
                  end do

               end if

            end do

            if (s% model_number == -1) then
               k = -1
               write(*,2) 've0(k)', k, ve0(k)
               write(*,2) 'H_T(k)', k, H_T(k)
               write(*,2) 'r(k)', k, r(k)
               write(*,2) 'Hj(k)', k, Hj(k)
               write(*,2) 'omega(k)', k, omega(k)
               write(*,2) 'dlnR_domega(k)', k, dlnR_domega(k)
               write(*,2) 'dr(k)', k, dr(k)
               write(*,2) 'H_gsf(k)', k, H_gsf(k)
               write(*,2) 'v_gsf(k)', k, v_gsf(k)
               write(*,2) 'csound(k)', k, csound(k)
               write(*,2) 'scale_height(k)', k, scale_height(k)
               write(*,2) 's% D_GSF(k)', k, s% D_GSF(k)
            end if

         end subroutine set_D_GSF


         logical function failed(str, ierr)
            character (len=*), intent(in) :: str
            integer, intent(in) :: ierr
            if (ierr == 0) then
               failed = .false.
               return
            end if
            write(s% retry_message,*) 'set_rotation_mixing_info failed in call to ' // trim(str)
            if (s% report_ierr) write(*, *) s% retry_message
            failed = .true.
         end function failed


      end subroutine set_rotation_mixing_info


      subroutine smooth_for_rotation(s, v, width, work)
         use star_utils, only: weighed_smoothing
         type (star_info), pointer :: s
         real(dp), dimension(:) :: v, work
         integer :: width
         logical, parameter :: preserve_sign = .false.
         if (width <= 0) return
         call weighed_smoothing(v, s% nz, width, preserve_sign, work)
      end subroutine smooth_for_rotation


      subroutine set_ST(s, &
            rho, T, r, L, omega, Cp, abar, zbar, delta, grav, &
            N2, N2_mu, opacity, kap_cond, scale_height, &
            ierr)
         ! with modifications by S.-C. Yoon, July 2003
         type (star_info), pointer :: s
         real(dp), dimension(:) :: & ! allocated temporary storage
            rho, T, r, L, omega, Cp, abar, zbar, delta, grav, &
            N2, N2_mu, opacity, kap_cond, scale_height
         integer, intent(out) :: ierr

         integer :: nz, k, j, kk
         real(dp) :: xmagfmu, xmagft, xmagfdif, xmagfnu, &
            xkap, xgamma, xlg, xsig1, xsig2, xsig3, xxx, ffff, xsig, &
            xeta, xmagn, xmagnmu, xmagnt, xmagw, xmagdn, xmagtn, xmagrn, xmag4pd, &
            dlnomega_dlnr, xmagq, xmager2w, &
            xmagnn, xmagwn, xmagq0, xmagwa0, xmags0, xmagbphi0, xmagbr0, xmageta0, &
            xmagkr2n, xmagq1, xmagwa1, xmags1a, xmags1b, xmags1, &
            xmagbphi1, xmagbr1, xmageta1a, xmageta1b, xmageta1, &
            xmagsm, xmagqm, xmagetam, xmagsf, xmags, xmagnu, xmagdif

         include 'formats'

         ierr = 0
         nz = s% nz

         s% D_ST(1:nz) = 0
         s% nu_ST(1:nz) = 0
         s% dynamo_B_r(1:nz) = 0
         s% dynamo_B_phi(1:nz) = 0

         xmagfmu = 1
         xmagft = 1
         xmagfdif = 1
         xmagfnu = 1

         xmageta0 = 0; xmageta1 = 0; xmagetam = 0
         xmagq0 = 0; xmagq1 = 0; xmagqm = 0
         xmags0 = 0; xmags1 = 0; xmagsm = 0;

         do k = 2, nz-1

            xkap = 16d0*boltz_sigma*T(k)*T(k)*T(k)/ &
               (3d0*opacity(k)*rho(k)*rho(k)*Cp(k)) ! thermal diffusivity

            xsig = calc_sige(abar(k), zbar(k), rho(k), T(k), Cp(k), kap_cond(k), opacity(k))
            xeta = calc_eta(xsig)  ! magnetic diffusivity

            xmagn = N2(k)
            xmagnmu = N2_mu(k)
            xmagnt = xmagn - xmagnmu ! N2_T
            xmagw = abs(omega(k))
            xmagdn = rho(k)
            xmagtn = T(k)
            xmagrn = r(k)
            xmag4pd = sqrt(pi4*xmagdn)

            dlnomega_dlnr = &
               (omega(k-1) - omega(k+1))/(r(k-1) - r(k+1))*(r(k)/omega(k))
            xmagq = max(1d-30,min(1d30,abs(dlnomega_dlnr))) ! shear
            xmager2w = xeta/(r(k)*r(k)*abs(omega(k)))

            ! magnetic quantities
            if (xmagnmu > 0.0D0) then
               xmagnn = xmagnmu*xmagfmu ! N^2
               xmagwn = xmagw/sqrt(xmagnn) ! omega/N
               xmagq0 = pow(xmagwn,-1.5D0)*pow(xmager2w,0.25D0) ! q0
               xmagwa0 = xmagq*xmagwn*xmagw ! omega_A
               xmags0 = xmagdn*xmagw*xmagw*xmagrn*xmagrn*pow3(xmagq)*pow4(xmagwn) ! S_0
               if (is_bad(xmags0)) then
                  do kk=1,s% nz
                     write(*,2) 'omega(kk)', kk, omega(kk)
                  end do
                  write(*,2) 'dlnomega_dlnr', k, dlnomega_dlnr
                  write(*,2) 'xmagfmu', k, xmagfmu
                  write(*,2) 'xmagnmu', k, xmagnmu
                  write(*,2) 'xmagnn', k, xmagnn
                  write(*,2) 'xmagwn', k, xmagwn
                  write(*,2) 'xmagq', k, xmagq
                  write(*,2) 'xmagrn', k, xmagrn
                  write(*,2) 'xmagw', k, xmagw
                  write(*,2) 'xmagdn', k, xmagdn
                  write(*,2) 'xmags0', k, xmags0
                  call mesa_error(__FILE__,__LINE__,'set_ST')
               end if
               xmagbphi0 = xmagwa0*xmag4pd*xmagrn ! B_\phi
               xmagbr0 = xmagbphi0*xmagq*xmagwn*xmagwn ! B_r
               xmageta0 = pow4(xmagq)*pow6(xmagwn)*xmagrn*xmagrn*xmagw ! eta_e
            end if

            if (xmagnt > 0.0D0) then
               xmagnn = xmagnt*xmagft ! N^2
               xmagwn = xmagw/sqrt(xmagnn) ! omega/N
               xmagkr2n = xkap/(xmagrn*xmagrn*sqrt(xmagnn)) ! kappa/(r^2 N)
               
               !xmagq1 = pow(xmager2w*xmagwn*pow3(xeta/xkap)/pow7(xmagwn),0.25D0) ! q_1
               if(xmagwn > 1d-42) then ! fix from rob
                  xmagq1 = pow(xmager2w*xmagwn*pow3(xeta/xkap)/pow7(xmagwn),0.25D0) ! q_1
               else
                  xmagq1 = 0d0
               end if
               
               xmagwa1 = sqrt(xmagq)*xmagw*pow(xmagwn*xmagkr2n,0.125D0) ! \omega_A
               xmags1a = xmagdn*pow2(xmagw*xmagrn)*xmagq*sqrt(xmagwn*xmagkr2n) ! S_1a
               xmags1b = xmagdn*pow2(xmagw)*pow2(xmagrn)*pow3(xmagq)*pow4(xmagwn) ! S_1b
               xmags1 = max(xmags1a,xmags1b) ! S_1
               xmagbphi1 = xmagwa1*xmag4pd*xmagrn ! B_\phi
               xmagbr1 = xmagbphi1*pow(xmagwn*xmagkr2n,0.25D0) ! B_r
               xmageta1a = xmagrn*xmagrn*xmagw*xmagq*pow(xmagwn*xmagkr2n,0.75D0) ! eta_e1a
               xmageta1b = pow4(xmagq)*pow6(xmagwn)*pow2(xmagrn)*xmagw ! eta_e1b
               xmageta1 = max(xmageta1a,xmageta1b) ! eta_e
            end if

            if ((xmagnt > 0.D0) .and. (xmagnmu > 0.D0) .and. xmags0+xmags1 /= 0d0) then
               xmagsm=xmags0*xmags1/(xmags0+xmags1) ! S_m
               if (is_bad(xmagsm)) then
                  write(*,2) 'xmags0', k, xmags0
                  write(*,2) 'xmags1', k, xmags1
                  write(*,2) 'xmags0+xmags1', k, xmags0+xmags1
                  write(*,2) 'xmags1', k, xmags1
                  write(*,2) 'xmags0', k, xmags0
                  write(*,2) 'xmagsm', k, xmagsm
                  call mesa_error(__FILE__,__LINE__,'set_ST')
               end if
               xmagqm=xmagq0+xmagq1 ! q_m
               xmagetam=xmageta0*xmageta1/(xmageta0+xmageta1) ! eta_m
            else if ((xmagnt <= 0.D0) .and. (xmagnmu > 0.D0)) then
               xmagsm=xmags0
               xmagqm=xmagq0
               xmagetam=xmageta0
            else if ((xmagnt > 0.D0) .and. (xmagnmu <= 0.D0)) then
               xmagsm=xmags1
               xmagqm=xmagq1
               xmagetam=xmageta1
            else if ((xmagnt <= 0.D0) .and. (xmagnmu <= 0.D0)) then
               xmagsm=0.0D0
               xmagqm=0.0D0
               xmagetam=0.0D0
            end if

            if (xmagn < 0.0D0) then
               xmagsm=0.0D0
               xmagqm=0.0D0
               xmagetam=0.0D0
            end if

            if (s% mixing_type(k) == convective_mixing .and. (xmagn > 0.D0)) then
               xmagsm=0.0D0
               xmagqm=0.0D0
               xmagetam=0.0D0
            end if

            if (s% mixing_type(k) == overshoot_mixing) then
               xmagsm=0.0D0
               xmagqm=0.0D0
               xmagetam=0.0D0
            end if

            xmagsf=1.D0-MIN(1.D0,xmagqm/xmagq)
            xmags=xmagsf*xmagsm ! S
            xmagnu=xmags/(xmagw*xmagq*xmagdn) ! \nu_e
            xmagdif=xmagsf*xmagetam ! \kappa_c (molecular diffusion)
               ! (set equal to eta_e and compute "mean")

            !..... Special treatment for semiconvective regions:
            !..... Geometric mean between semiconvective "effective viscosity"
            !..... and \nu_e as resulting from \mu-dominated case
            !..... (the T-dominated case is indefined in these regions).
            !..... Assume flux is dominated by convective flux and given by
            !..... the luminosity of the star at this point.
            !..... From Kippenhahn&Weigert, Eqs. (7.6) and (7.7), the convective
            !..... velocity can be computed.  Now assume \nu=1/3 l_mix v and assume
            !..... that l_mix = H_p.
            if ((xmagnt .LE. 0.D0) .AND. (xmagnmu .GT. 0.D0) .AND. (xmagn .GT. 0.D0)) &
               xmagnu = &
                  sqrt(xmagnu*scale_height(k)*(one_third)* &
                  pow(grav(k)*delta(k)*scale_height(k)*MAX(0.D0,L(k))/ &
                         (64.D0*pi*xmagdn*Cp(k)*xmagtn*xmagrn*xmagrn),one_third))

            s% D_ST(k) = min(xmagdif,xmagnu)*xmagfdif
            s% nu_ST(k) = xmagnu*xmagfnu

            if (is_bad(s% D_ST(k))) then
               write(*,4) 'set_ST mixing_type', k, s% model_number, s% mixing_type(k)
               write(*,3) 'set_ST xeta', k, s% model_number, xeta
               write(*,3) 'set_ST xsig', k, s% model_number, xsig
               write(*,3) 'set_ST dlnomega_dlnr', k, s% model_number, dlnomega_dlnr
               write(*,3) 'set_ST xmager2w', k, s% model_number, xmager2w
               write(*,3) 'set_ST xkap', k, s% model_number, xkap
               write(*,3) 'set_ST xgamma', k, s% model_number, xgamma
               write(*,3) 'set_ST xmagq', k, s% model_number, xmagq
               write(*,3) 'set_ST xmagqm', k, s% model_number, xmagqm
               write(*,3) 'set_ST xmagsf', k, s% model_number, xmagsf
               write(*,3) 'set_ST xmagsm', k, s% model_number, xmagsm
               write(*,3) 'set_ST xmags', k, s% model_number, xmags
               write(*,3) 'set_ST xmagnu', k, s% model_number, xmagnu
               write(*,3) 'set_ST xmagfdif', k, s% model_number, xmagfdif
               write(*,3) 'set_ST D_ST(k)', k, s% model_number, s% D_ST(k)
               call mesa_error(__FILE__,__LINE__,'set_ST')
            end if

         end do

         s% D_ST(1) = s% D_ST(2)
         s% D_ST(nz) = s% D_ST(nz-1)

         s% nu_ST(1) = s% nu_ST(2)
         s% nu_ST(nz) = s% nu_ST(nz-1)

      end subroutine set_ST


      subroutine zero_if_convective(nz, mixing_type, D_mix, dc)
         integer, intent(in) :: nz
         integer, dimension(:) :: mixing_type
         real(dp), dimension(:) :: D_mix, dc
         integer :: k
         do k=1,nz
            if (mixing_type(k) == convective_mixing) dc(k) = 0
         end do
      end subroutine zero_if_convective


      subroutine zero_if_tiny(s, dc)
         type (star_info), pointer :: s
         real(dp), dimension(:) :: dc
         integer :: k
         real(dp) :: tiny
         tiny = s% clip_D_limit
         do k=1,s% nz
            if (dc(k) < tiny) dc(k) = 0
         end do
      end subroutine zero_if_tiny


      end module rotation_mix_info
