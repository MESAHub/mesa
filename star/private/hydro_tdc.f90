! ***********************************************************************
!
!   Copyright (C) 2010-2020  Bill Paxton & The MESA Team
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

      module hydro_tdc

      use star_private_def
      use const_def
      use utils_lib, only: is_bad
      use auto_diff
      use auto_diff_support
      use accurate_sum_auto_diff_star_order1
      use star_utils, only: em1, e00, ep1

      implicit none

      private
      public :: do1_tdc_L_eqn, do1_turbulent_energy_eqn, compute_Eq_cell, &
         compute_Uq_face, set_TDC_vars, &
         set_using_TDC, set_etrb_start_vars, reset_etrb_using_L
      
      real(dp), parameter :: &
         x_ALFAP = 2.d0/3.d0, & ! Ptrb
         x_ALFAS = (1.d0/2.d0)*sqrt_2_div_3, & ! PII_face and Lc
         x_ALFAC = (1.d0/2.d0)*sqrt_2_div_3, & ! Lc
         x_CEDE  = (8.d0/3.d0)*sqrt_2_div_3, & ! DAMP
         x_GAMMAR = 2.d0*sqrt(3.d0) ! DAMPR

      contains
      
      
      subroutine set_using_TDC(s)
         type (star_info), pointer :: s      
         real(dp) :: alfa, beta
         call get_TDC_frac(s, alfa, beta)
         s% using_TDC = (alfa > 0d0)
      end subroutine set_using_TDC
      
      
      subroutine get_TDC_frac(s, alfa, beta)
         type (star_info), pointer :: s
         real(dp), intent(out) :: alfa, beta
         real(dp) :: dt, switch
         include 'formats'
         if (.not. s% TDC_flag) then
            alfa = 0d0; beta = 1d0
            return
         end if
         dt = s% dt
         if (s% TDC_min_dt_div_tau_conv_switch_to_MLT > 0) then
            switch = s% max_conv_time_scale*s% TDC_min_dt_div_tau_conv_switch_to_MLT
         else if (s% TDC_min_dt_years_switch_to_MLT > 0) then
            switch = s% TDC_min_dt_years_switch_to_MLT*secyer
         else
            switch = 0d0
         end if
         
         if (.false.) then ! blending
            !if (dt >= s% TDC_dt_seconds_etrb_all_MLT) then
            !   alfa = 0d0 ! all MLT, no TDC
            !else if (dt <= s% TDC_dt_seconds_etrb_no_MLT) then
            !   alfa = 1d0 ! no MLT, all TDC
            !else ! blend 0 < alfa < 1 = fraction TDC, beta = 1 - alfa = fraction MLT
            !   alfa = (s% TDC_dt_seconds_etrb_all_MLT - dt)/ &
            !       (s% TDC_dt_seconds_etrb_all_MLT - s% TDC_dt_seconds_etrb_no_MLT)
            !end if
         else
            if (dt >= switch) then
               alfa = 0d0 ! all MLT, no TDC
            else
               alfa = 1d0 ! no MLT, all TDC
            end if
         end if
         
         beta = 1d0 - alfa
      end subroutine get_TDC_frac
      
      
      subroutine set_TDC_vars(s,ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr    
         type(auto_diff_real_star_order1) :: x
         integer :: k
         include 'formats'
         ierr = 0
         if (s% need_to_reset_w) then
            write(*,2) 'reset_etrb_using_L', s% model_number
            call reset_etrb_using_L(s,ierr)
            if (ierr /= 0) then
               stop 'failed in reset_etrb_using_L'
               return
            end if
            s% need_to_reset_w = .false.
         end if
         do k=1,s%nz
            x = compute_Hp_face(s, k, ierr) ! sets Hp_face
            if (ierr /= 0) return
            x = compute_L_face(s, k, ierr) ! sets Lr, Lt, Lc
            if (ierr /= 0) return
         end do
         s% Y_face(1) = 0d0
         s% PII(1) = 0d0
         s% Lc(1) = 0d0
         s% Lt(1) = 0d0
      end subroutine set_TDC_vars
      

      subroutine do1_tdc_L_eqn(s, k, skip_partials, nvar, ierr)
         use star_utils, only: save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr         
         type(auto_diff_real_star_order1) ::  &
            L_expected, L_actual,resid
         real(dp) :: scale, residual, L_start_max
         logical :: test_partials
         include 'formats'

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         
         if (.not. s% using_TDC) then
            ierr = -1
            return
         end if

         ierr = 0
         L_expected = compute_L_face(s, k, ierr)
         if (ierr /= 0) return        
         L_actual = wrap_L_00(s, k)  
         L_start_max = maxval(s% L_start(1:s% nz))
         scale = 1d0/L_start_max
         if (is_bad(scale)) then
            write(*,2) 'do1_tdc_L_eqn scale', k, scale
            stop 'do1_tdc_L_eqn'
         end if
         resid = (L_expected - L_actual)*scale         
      
         residual = resid%val
         s% equ(s% i_equL, k) = residual         
         if (test_partials) then
            s% solver_test_partials_val = residual 
         end if
         
         if (skip_partials) return
         call save_eqn_residual_info(s, k, nvar, s% i_equL, resid, 'do1_tdc_L_eqn', ierr)
         if (ierr /= 0) return

         if (test_partials) then
            s% solver_test_partials_var = s% i_lnR
            s% solver_test_partials_dval_dx = resid%d1Array(i_lnR_00)
            write(*,4) 'do1_tdc_L_eqn', s% solver_test_partials_var
         end if      
      end subroutine do1_tdc_L_eqn


      subroutine eval_xis(s, k, xi0, xi1, xi2)
         ! Inputs
         type (star_info), pointer :: s
         integer, intent(in) :: k

         ! Intermediates
         integer :: ierr
         real(dp) :: gammar, alpha, POM, ALFAM_ALFA
         type(auto_diff_real_star_order1) :: tst, w_00, T_00, d_00, Cp_00, kap_00, Hp_cell, POM2
         type(auto_diff_real_star_order1) :: Peos_00, chiT_00, chiRho_00, QQ_00, grad_ad_00, fac
         type(auto_diff_real_star_order1) :: Hp_face_00, Hp_face_p1, PII_face_00, PII_face_p1, PII_div_Hp_cell
         type(auto_diff_real_star_order1) :: Chi_div_w_cell, d_v_div_r, r_00, r_p1, f, r6_cell
         type(auto_diff_real_star_order1) :: S0, Eq0, DR0, Pt0, D0, dV

         ! Outputs
         type(auto_diff_real_star_order1) :: xi0, xi1, xi2

         logical :: test_partials
         include 'formats'

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         ! Wrap and store basic variables
         alpha = s% TDC_alfa
         gammar = s% TDC_alfar*x_GAMMAR
         w_00 = wrap_w_00(s,k)
         T_00 = wrap_T_00(s,k)
         d_00 = wrap_d_00(s,k)
         Cp_00 = wrap_Cp_00(s,k)
         kap_00 = wrap_kap_00(s,k)
         Hp_cell = compute_Hp_cell(s,k,ierr)
         Peos_00 = wrap_Peos_00(s, k)         
         chiT_00 = wrap_chiT_00(s, k)
         chiRho_00 = wrap_chiRho_00(s, k)
         QQ_00 = chiT_00/(d_00*T_00*chiRho_00)
         Hp_face_00 = compute_Hp_face(s,k,ierr)
         PII_face_00 = compute_PII_face(s, k, ierr)
         grad_ad_00 = wrap_grad_ad_00(s, k)

         ! Evaluate S0
         if (k == s% nz) then
            PII_div_Hp_cell = PII_face_00/Hp_face_00
         else
            Hp_face_p1 = shift_p1(compute_Hp_face(s,k+1,ierr))
            if (ierr /= 0) return
            PII_face_p1 = shift_p1(compute_PII_face(s, k+1, ierr))
            if (ierr /= 0) return
            PII_div_Hp_cell = 0.5d0*(PII_face_00/Hp_face_00 + PII_face_p1/Hp_face_p1)
         end if
         S0 = PII_div_Hp_cell*T_00*grad_ad_00 ! = S / w (only equal if TDC_SEED=0)
         ! PII units same as Cp = erg g^-1 K^-1
         ! Source units = (erg g^-1 K^-1) cm^-1 K
         !     = (erg/g)/(cm) = (cm/s^2)

         ! Evaluate Eq0 = (Eq / w)
         if (k <= s% TDC_num_outermost_cells_forced_nonturbulent .or. &
             k > s% nz - s% TDC_num_innermost_cells_forced_nonturbulent) then
            Eq0 = 0d0
         else
            ALFAM_ALFA = s% TDC_alfam*s% TDC_alfa
            if (ALFAM_ALFA == 0d0 .or. &
                  k <= s% TDC_num_outermost_cells_forced_nonturbulent .or. &
                  k > s% nz - s% TDC_num_innermost_cells_forced_nonturbulent) then
               Chi_div_w_cell = 0d0
            else
               r_00 = wrap_r_00(s,k) ! Changed to be fully backwards Euler (not time centered)
               r_p1 = wrap_r_p1(s,k) ! Changed to be fully backwards Euler (not time centered)
               d_v_div_r = compute_d_v_div_r(s, k, ierr)
               f = (16d0/3d0)*pi*ALFAM_ALFA/s% dm(k)  
               r6_cell = 0.5d0*(pow6(r_00) + pow6(r_p1))
               Chi_div_w_cell = f*pow2(d_00)*r6_cell*d_v_div_r*Hp_cell
               ! units = (1/g)(g/cm^3)(g/cm^3)(cm^6)(cm/s/cm)(cm) = g cm/s = erg/(cm/s)
            end if

            if (ierr /= 0) return
            Eq0 = 4d0*pi*Chi_div_w_cell*d_v_div_r/s% dm(k)
            ! units = (erg/(cm/s)) (cm/s/cm) / g = (erg/g)/cm = cm/s^2
         end if

         ! Evaluate xi0
         xi0 = Eq0 + S0

         ! Evaluate DR0 = (Dr / etrb)
         POM = 4d0*boltz_sigma*(gammar/alpha)**2 ! erg cm^-2 K^-4 s^-1
         POM2 = pow3(T_00)/(pow2(d_00)*Cp_00*kap_00) 
         DR0 = POM*POM2/pow2(Hp_cell)

         ! Evaluate xi1
         Pt0 = s% TDC_alfap*d_00 ! Pturb / w^2
         dV = 1d0/d_00 - 1d0/s% rho_start(k)
         xi1 = -(DR0 + Pt0 * dV)

         ! Evaluate xi2
         D0 = (s% TDC_alfad*x_CEDE/alpha) / Hp_cell ! (D / w^3)
         xi2 = -D0
         ! units = (1/cm)

         if (k < -60 .and. k > 30) then
            write(*,*) D0%val, Eq0%val, S0%val
         end if
         
         if (test_partials) then
            tst = xi0
            s% solver_test_partials_val = tst%val
            s% solver_test_partials_var = s% i_lnd
            s% solver_test_partials_dval_dx = tst%d1Array(i_lnd_00)
            write(*,*) 'eval_xis', s% solver_test_partials_var, s% lnd(k)/ln10, tst%val
         end if

      end subroutine eval_xis
      
      !> Returns the smallest positive z such that tan(z) = y/x
      type(auto_diff_real_star_order1) function two_var_pos_atan(x,y) result(z)
         type(auto_diff_real_star_order1), intent(in) :: x,y
         type(auto_diff_real_star_order1) :: x1, y1

         x1 = abs(x) + 1d-50
         y1 = abs(y) + 1d-50
         z = atan(y1/x1)
         if (z < 0d0) then
            z = z + pi
         end if
      end function two_var_pos_atan

      subroutine eval_Af_from_A0(s, k, xi0, xi1, xi2, A0, Af, dt, i)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: dt
         integer, intent(in) :: i
         type(auto_diff_real_star_order1), intent(in) :: xi0, xi1, xi2, A0
         type(auto_diff_real_star_order1), intent(out) :: Af
         
         integer :: l
         type(auto_diff_real_star_order1) :: J2, J, Jt, em1, ep1, num, den, root, y, tst
         logical :: test_partials
         include 'formats'

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         J2 = pow2(xi1) - 4d0 * xi0 * xi2
         !if (test_partials) write(*,2) 'J2', k, J2%val, pow2(xi1%val), 4d0 * xi0%val * xi2%val

         if (J2 > 0d0) then
            ! Hyperbolic branch

            J = sqrt(J2)
            Jt = dt * J

            em1 = exp(min(30d0, 0.5d0 * Jt)) - 1d0
            ep1 = exp(min(30d0, 0.5d0 * Jt)) + 1d0

            num = em1 * (J2 - pow2(xi1)) - 2d0 * A0 * xi2 * (ep1 * J + em1 * xi1)
            den = 2d0 * xi2 * (-ep1 * J + em1 * (xi1 + 2d0 * A0 * xi2))

            Af = num / den
            if (Af < 0d0) then
               Af = -Af
            end if

            if (test_partials) then
               tst = Af
               s% solver_test_partials_val = tst%val
               s% solver_test_partials_var = s% i_lnd
               s% solver_test_partials_dval_dx = tst%d1Array(i_lnd_00)
               write(*,*) 'eval_Af_from_A0', s% solver_test_partials_var, &
                  s% lnd(k)/ln10, tst%val, J2%val
            end if

            if (i > 30 .and. i < -60) then
               write(*,*) 'k          ',i
               write(*,*) 'xis        ',xi0%val, xi1%val, xi2%val
               write(*,*) 'num,den    ',num%val, den%val
               write(*,*) 'J2,J       ',J2%val, J%val
               write(*,*) 'em1,ep1    ',em1%val, ep1%val
               write(*,*) 'A0,Af,xiR  ',A0%val,Af%val,sqrt(abs(xi0%val/xi2%val))
            end if

            do l=1,27
               if (is_bad(Af%d1Array(l))) then
                  !$omp critical (eval_Af_crit)
                  write(*,*) 'Bad partial in Af in Hyperbolic branch.'
                  write(*,*) 'Partial      ',l
                  write(*,*) 'Af           ',Af%val
                  write(*,*) 'Af Partial   ',Af%d1Array(l)
                  write(*,*) 'xi0          ',xi0%val
                  write(*,*) 'xi0 Partial  ',xi0%d1Array(l)
                  write(*,*) 'xi1          ',xi1%val
                  write(*,*) 'xi1 Partial  ',xi1%d1Array(l)
                  write(*,*) 'xi2          ',xi2%val
                  write(*,*) 'xi2 Partial  ',xi2%d1Array(l)
                  write(*,*) 'num          ',num%val
                  write(*,*) 'num Partial  ',num%d1Array(l)
                  write(*,*) 'den          ',den%val
                  write(*,*) 'den Partial  ',den%d1Array(l)
                  stop 'eval_Af_from_A0'
                  !$omp end critical (eval_Af_crit)
                  end if
            end do

         else
            ! Trigonometric branch

            J = sqrt(-J2)
            Jt = dt * J

            ! This branch contains decaying solutions that reach A = 0, at which point
            ! they switch onto the 'zero' branch. So we have to calculate the position of
            ! the first root to check it against dt.
            y = xi1 + 2d0 * A0 * xi2

            ! We had a choice above to pick which of +-I to use in switching branches.
            ! That choice has to be consistent with a decaying solution, which we check now.
            if ((J > 0d0 .and. y < 0d0) .or. (J < 0d0 .and. y > 0d0)) then
               J = -J
            end if

            root = two_var_pos_atan(J, y) - two_var_pos_atan(J, xi1)

            if (0.25d0 * Jt < root) then
               num = -xi1 + J * tan(0.25d0 * Jt + atan(y / J)) 
               den = 2d0 * xi2
               Af = num / den
            else
               Af = 0d0
            end if

            do l=1,27
               if (is_bad(Af%d1Array(l))) then
                  !$omp critical (eval_Af_crit)
                  write(*,*) 'Bad partial in Af in Trigonometric branch.'
                  write(*,*) 'Partial      ',l
                  write(*,*) 'Af           ',Af%val
                  write(*,*) 'Af Partial   ',Af%d1Array(l)
                  write(*,*) 'xi0          ',xi0%val
                  write(*,*) 'xi0 Partial  ',xi0%d1Array(l)
                  write(*,*) 'xi1          ',xi1%val
                  write(*,*) 'xi1 Partial  ',xi1%d1Array(l)
                  write(*,*) 'xi2          ',xi2%val
                  write(*,*) 'xi2 Partial  ',xi2%d1Array(l)
                  write(*,*) 'root         ',root%val
                  write(*,*) 'root Partial ',root%d1Array(l)
                  write(*,*) 'num          ',num%val
                  write(*,*) 'num Partial  ',num%d1Array(l)
                  write(*,*) 'den          ',den%val
                  write(*,*) 'den Partial  ',den%d1Array(l)
                  stop 'eval_Af_from_A0'
                  !$omp end critical (eval_Af_crit)
                  end if
            end do

            !write(*,*) 'trig',xi0%val,xi1%val,xi2%val,J%val,root%val,dt,num%val, den%val,Af%val

         end if

      end subroutine eval_Af_from_A0

      subroutine do1_turbulent_energy_eqn(s, k, skip_partials, nvar, ierr)
         use star_utils, only: set_energy_eqn_scal, save_eqn_residual_info
         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr         
         integer :: j
         type(auto_diff_real_star_order1) :: xi0, xi1, xi2, A0, Af, w_00, LHS, RHS
         type(auto_diff_real_star_order1) :: tst, resid_ad, &
            d_turbulent_energy_ad, Ptrb_dV_ad, dt_dLt_dm_ad, dt_C_ad, dt_Eq_ad, &
            etrb_mlt_cell, esum_mlt, esum_tdc
         type(accurate_auto_diff_real_star_order1) :: esum_ad
         logical :: non_turbulent_cell, test_partials
         real(dp) :: scal, residual, alfa, beta
         include 'formats'

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         
         ierr = 0
         
         non_turbulent_cell = &
            s% TDC_alfa == 0d0 .or. &
            k <= s% TDC_num_outermost_cells_forced_nonturbulent .or. &
            k > s% nz - s% TDC_num_innermost_cells_forced_nonturbulent
         
         if (.not. s% using_TDC) then
             
            resid_ad = wrap_w_00(s,k) - s% w_start(k)

         else if (non_turbulent_cell) then

            resid_ad = wrap_w_00(s,k)/(1d2*s% csound_start(k)) ! make w = 0
            
         else
         
            call get_TDC_frac(s, alfa, beta) ! alfa is frac TDC, beta is frac MLT

            ! always compute etrb_mlt_cell.  for debugging
            !if (beta > 0d0) then
               etrb_mlt_cell = compute_etrb_mlt_cell(s, k, ierr); if (ierr /= 0) return ! erg g^-1
               esum_mlt = etrb_mlt_cell - wrap_etrb_00(s,k)
            !else
            !   etrb_mlt_cell = 0d0
            !   esum_mlt = 0d0
            !end if
            
            !if (alfa > 0d0) then
               call setup_d_turbulent_energy(ierr); if (ierr /= 0) return ! erg g^-1 = cm^2 s^-2
               call setup_Ptrb_dV_ad(ierr); if (ierr /= 0) return ! erg g^-1
               call setup_dt_dLt_dm_ad(ierr); if (ierr /= 0) return ! erg g^-1
               call setup_dt_C_ad(ierr); if (ierr /= 0) return ! erg g^-1
               call setup_dt_Eq_ad(ierr); if (ierr /= 0) return ! erg g^-1
               call set_energy_eqn_scal(s, k, scal, ierr); if (ierr /= 0) return  ! 1/(erg g^-1 s^-1)         
               ! sum terms in esum_ad using accurate_auto_diff_real_star_order1
               esum_ad = d_turbulent_energy_ad + Ptrb_dV_ad + dt_dLt_dm_ad - dt_C_ad - dt_Eq_ad ! erg g^-1
               esum_tdc = esum_ad ! convert back to auto_diff_real_star_order1
            !else
            !   esum_tdc = 0d0
            !end if
            resid_ad = alfa*esum_tdc + beta*esum_mlt
            resid_ad = resid_ad*scal/s%dt ! to make residual unitless, must cancel out the dt in scal

       
            if (.true.) then ! NEW WAY
            
               call eval_xis(s, k, xi0, xi1, xi2)

               A0 = sqrt(max(0d0, get_etrb_start(s,k)))
               call eval_Af_from_A0(s, k, xi0, xi1, xi2, A0, Af, s%dt, k)

               ! sum terms in esum_ad using accurate_auto_diff_real_star_order1
               ! 2*w*(w - A) = 2*w*B*(xi1 + (2*A+B)*xi2) - 2*w*Lambda
               !
               ! B = (w-A) so
               !
               ! 2*w*(w - A) = 2*w*(w-A)*(xi1 + (w+A)*xi2) - 2*w*Lambda
               !
               ! -2*w*Lambda = -dt_dLt_dm_ad

               w_00 = wrap_w_00(s,k)
               if (.false. .and. s% TDC_alfat == 0d0) then
                  esum_ad = w_00 - Af
               else
                  esum_ad = (w_00 - Af) * (1d0 - s%dt * (xi1 + (w_00 + Af) * xi2)) + dt_dLt_dm_ad / (2d0 * (1d3 + w_00))
               end if
               resid_ad = esum_ad
            
            end if

            if (k > 30 .and. k < -60) then
               write(*,*) k,resid_ad%val, w_00%val, Af%val, xi0%val, xi1%val, xi2%val, sqrt(abs(xi0%val/xi2%val)), dt_dLt_dm_ad%val, A0%val
            end if

            scal = 1d0/s%csound_start(k)
            !scal = 1d0/(1d3 + s%csound_start(k) + maxval(s% w_start(1:s%nz)))
            !scal = sqrt(scal/s%dt) 
            !scal = 1d0/(1d4 + s% dt + pow(4d0,s% solver_iter))
               ! scal/dt -> 1/(erg g^-1 s^-1)*s^-1 = 1/(erg g^-1) = 1/(g cm^2 sec^-2 g^-1) = sec^2 cm^-2
               ! sqrt(scal/dt) => sec/cm = 1 / (cm/sec).   resid units are cm/sec.
            resid_ad = resid_ad*scal

            do j=1,27
               if (is_bad(resid_ad%d1Array(j))) then
                  write(*,*) 'Cell         ',k
                  write(*,*) 'Partial      ',j
                  write(*,*) 'Residual     ',resid_ad%val
                  write(*,*) 'Bad Partial  ',resid_ad%d1Array(j)
                  write(*,*) 'w_00         ',w_00%val
                  write(*,*) 'w_00 Partial ',w_00%d1Array(j)
                  write(*,*) 'Af           ',Af%val
                  write(*,*) 'Af Partial   ',Af%d1Array(j)
                  write(*,*) 'xi0          ',xi0%val
                  write(*,*) 'xi0 Partial  ',xi0%d1Array(j)
                  write(*,*) 'xi1          ',xi1%val
                  write(*,*) 'xi1 Partial  ',xi1%d1Array(j)
                  write(*,*) 'xi2          ',xi2%val
                  write(*,*) 'xi2 Partial  ',xi2%d1Array(j)
                  stop
               end if
            end do

         end if

         residual = resid_ad%val
         s% equ(s% i_detrb_dt, k) = residual

         if (test_partials) then
            tst = residual
            s% solver_test_partials_val = tst%val
            if (s% solver_iter == 12) &
               write(*,*) 'do1_turbulent_energy_eqn', s% solver_test_partials_var, s% lnd(k), tst%val
         end if
         
         if (skip_partials) return
         call save_eqn_residual_info(s, k, nvar, s% i_detrb_dt, resid_ad, 'do1_turbulent_energy_eqn', ierr)
         if (ierr /= 0) return

         if (test_partials) then
            s% solver_test_partials_var = s% i_lnd
            s% solver_test_partials_dval_dx = tst%d1Array(i_lnd_00)     ! xi0 good , xi1 partial 0, xi2 good.  Af horrible.'
            write(*,*) 'do1_turbulent_energy_eqn', s% solver_test_partials_var, s% lnd(k)/ln10, tst%val
         end if      

         contains
         
         subroutine setup_d_turbulent_energy(ierr) ! erg g^-1
            integer, intent(out) :: ierr
            ierr = 0
            d_turbulent_energy_ad = wrap_etrb_00(s,k) - get_etrb_start(s,k)
         end subroutine setup_d_turbulent_energy
         
         ! Ptrb_dV_ad = Ptrb_ad*dV_ad
         subroutine setup_Ptrb_dV_ad(ierr) ! erg g^-1
            use star_utils, only: calc_Ptrb_ad_tw
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: Ptrb_ad, dV_ad, d_00
            call calc_Ptrb_ad_tw(s, k, Ptrb_ad, ierr)
            if (ierr /= 0) return
            d_00 = wrap_d_00(s,k)
            dV_ad = 1d0/d_00 - 1d0/s% rho_start(k)
            Ptrb_dV_ad = Ptrb_ad*dV_ad ! erg cm^-3 cm^-3 g^-1 = erg g^-1
         end subroutine setup_Ptrb_dV_ad

         subroutine setup_dt_dLt_dm_ad(ierr) ! erg g^-1
            integer, intent(out) :: ierr            
            type(auto_diff_real_star_order1) :: Lt_00, Lt_p1, dLt_ad
            real(dp) :: Lt_00_start, Lt_p1_start
            logical :: time_centering
            include 'formats'
            ierr = 0
            time_centering = &
               s% using_velocity_time_centering .and. &
               s% include_L_in_velocity_time_centering
            Lt_00 = compute_Lt(s, k, ierr)
            if (ierr /= 0) return
            if (time_centering) Lt_00 = 0.5d0*(Lt_00 + s% Lt_start(k))
            if (k == s% nz) then
               Lt_p1 = 0d0
            else
               Lt_p1 = shift_p1(compute_Lt(s, k+1, ierr))
               if (ierr /= 0) return
               if (time_centering) Lt_p1 = 0.5d0*(Lt_p1 + s% Lt_start(k+1))
            end if
            dt_dLt_dm_ad = (Lt_00 - Lt_p1)*s%dt/s%dm(k)
         end subroutine setup_dt_dLt_dm_ad
         
         subroutine setup_dt_C_ad(ierr) ! erg g^-1
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: C
            C = compute_C(s, k, ierr) ! erg g^-1 s^-1
            if (ierr /= 0) return
            dt_C_ad = s%dt*C
         end subroutine setup_dt_C_ad
                  
         subroutine setup_dt_Eq_ad(ierr) ! erg g^-1
            integer, intent(out) :: ierr
            type(auto_diff_real_star_order1) :: Eq_cell
            Eq_cell = compute_Eq_cell(s, k, ierr) ! erg g^-1 s^-1
            if (ierr /= 0) return
            dt_Eq_ad = s%dt*Eq_cell
         end subroutine setup_dt_Eq_ad
      
      end subroutine do1_turbulent_energy_eqn


      function compute_etrb_mlt_cell(s, k, ierr) result(etrb_mlt) ! erg g^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: etrb_mlt
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: mlt_vc_00, mlt_vc_p1, vc_cell_mlt
         include 'formats'
         ierr = 0
         if (k > 1) then
            mlt_vc_00 = s% mlt_vc_ad(k)
         else
            mlt_vc_00 = 0d0
         end if
         if (k < s% nz) then
            mlt_vc_p1 = shift_p1(s% mlt_vc_ad(k+1))
         else
            mlt_vc_p1 = 0d0
         end if
         vc_cell_mlt = 0.5d0*(mlt_vc_00 + mlt_vc_p1)
         ! vc_mlt = sqrt(2/3)*w   
         etrb_mlt = 1.5d0*pow2(vc_cell_mlt)
      end function compute_etrb_mlt_cell
      
      
      function compute_Hp_cell(s, k, ierr) result(Hp_cell) ! cm
         ! instead of 0.5d0*(Hp_face(k) + Hp_face(k+1)) to keep block tridiagonal
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Hp_cell
         type(auto_diff_real_star_order1) :: r_mid, r_00, r_p1, &
            Peos_00, d_00, alt_Hp_cell, alfa
         real(dp) :: cgrav_00, cgrav_p1, cgrav_mid, m_00, m_p1, m_mid
         include 'formats'
         ierr = 0
         r_00 = wrap_opt_time_center_r_00(s, k)
         cgrav_00 = s% cgrav(k)
         m_00 = s% m(k)
         d_00 = wrap_d_00(s, k)
         Peos_00 = wrap_Peos_00(s, k)
         r_p1 = wrap_opt_time_center_r_p1(s, k)
         if (k < s% nz) then
            cgrav_p1 = s% cgrav(k+1)
            m_p1 = s% m(k+1)
         else
            cgrav_p1 = s% cgrav(k)
            m_p1 = s% m_center
         end if
         cgrav_mid = 0.5d0*(cgrav_00 + cgrav_p1)
         m_mid = 0.5d0*(m_00 + m_p1)
         r_mid = 0.5d0*(r_00 + r_p1)
         Hp_cell = pow2(r_mid)*Peos_00 / (d_00*cgrav_mid*m_mid)
         if (s% alt_scale_height_flag) then
            ! consider sound speed*hydro time scale as an alternative scale height
            alt_Hp_cell = sqrt(Peos_00/cgrav_mid)/d_00
            if (alt_Hp_cell%val < Hp_cell%val) then ! blend
               alfa = pow2(alt_Hp_cell/Hp_cell) ! 0 <= alfa%val < 1
               Hp_cell = alfa*Hp_cell + (1d0 - alfa)*alt_Hp_cell
            end if
         end if
      end function compute_Hp_cell
      
      
      subroutine get_TDC_alfa_beta_face_weights(s, k, alfa, beta)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: alfa, beta
         ! face_value = alfa*cell_value(k) + beta*cell_value(k-1)
         if (k == 1) stop 'bad k==1 for get_TDC_alfa_beta_face_weights'
         if (s% TDC_use_mass_interp_face_values) then
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            beta = 1d0 - alfa
         else
            alfa = 0.5d0
            beta = 0.5d0
         end if
      end subroutine get_TDC_alfa_beta_face_weights
      
      
      function compute_Hp_face(s, k, ierr) result(Hp_face) ! cm
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Hp_face
         type(auto_diff_real_star_order1) :: &
            rho_face, area, dlnPeos, &
            r_00, Peos_00, d_00, Peos_m1, d_m1, Peos_div_rho, &
            d_face, Peos_face, alt_Hp_face, A
         real(dp) :: alfa, beta
         integer :: j
         include 'formats'
         ierr = 0
         if (k > s% nz) then
            Hp_face = 1d0 ! not used
            s% Hp_face(k) = Hp_face%val
            return
         end if
         if (k > 1 .and. .not. s% TDC_assume_HSE) then
            call get_TDC_alfa_beta_face_weights(s, k, alfa, beta)
            rho_face = alfa*wrap_d_00(s,k) + beta*wrap_d_m1(s,k)
            area = 4d0*pi*pow2(wrap_r_00(s,k))
            dlnPeos = wrap_lnPeos_m1(s,k) - wrap_lnPeos_00(s,k)
            Hp_face = -s% dm_bar(k)/(area*rho_face*dlnPeos)
         else
            r_00 = wrap_opt_time_center_r_00(s, k)
            d_00 = wrap_d_00(s, k)
            Peos_00 = wrap_Peos_00(s, k)
            if (k == 1) then
               Peos_div_rho = Peos_00/d_00
               Hp_face = pow2(r_00)*Peos_div_rho/(s% cgrav(k)*s% m(k))
            else
               d_m1 = wrap_d_m1(s, k)
               Peos_m1 = wrap_Peos_m1(s, k)
               call get_TDC_alfa_beta_face_weights(s, k, alfa, beta)
               Peos_div_rho = alfa*Peos_00/d_00 + beta*Peos_m1/d_m1
               Hp_face = pow2(r_00)*Peos_div_rho/(s% cgrav(k)*s% m(k))
               if (s% alt_scale_height_flag) then
                  ! consider sound speed*hydro time scale as an alternative scale height
                  d_face = alfa*d_00 + beta*d_m1
                  Peos_face = alfa*Peos_00 + beta*Peos_m1
                  alt_Hp_face = sqrt(Peos_face/s% cgrav(k))/d_face
                  if (alt_Hp_face%val < Hp_face%val) then ! blend
                     A = pow2(alt_Hp_face/Hp_face) ! 0 <= A%val < 1
                     Hp_face = A*Hp_face + (1d0 - A)*alt_Hp_face
                  end if
               end if
            end if
         end if
         s% Hp_face(k) = Hp_face%val

      end function compute_Hp_face

      
      function compute_Y_face(s, k, ierr) result(Y_face) ! superadiabatic gradient [unitless]
         type (star_info), pointer :: s
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Y_face
         type(auto_diff_real_star_order1) :: Hp_face, Y1, Y2, QQ_div_Cp_face, &
            r_00, d_00, Peos_00, Cp_00, T_00, chiT_00, chiRho_00, QQ_00, lnT_00, &
            r_m1, d_m1, Peos_m1, Cp_m1, T_m1, chiT_m1, chiRho_m1, QQ_m1, lnT_m1, &
            dlnT_dlnP, grad_ad_00, grad_ad_m1, grad_ad_face, dlnT, dlnP, alt_Y_face
         real(dp) :: dm_bar, alfa, beta
         include 'formats'
         ierr = 0
         
         if (k > s% nz) then
            Y_face = 0d0
            return
         end if
         
         if (k == 1 .or. s% TDC_alfa == 0d0) then
            Y_face = 0d0
            s% Y_face(k) = 0d0
            return
         end if
         
         call get_TDC_alfa_beta_face_weights(s, k, alfa, beta)
         
         if (s% TDC_use_RSP_eqn_for_Y_face) then
      
            dm_bar = s% dm_bar(k)
            Hp_face = compute_Hp_face(s,k,ierr)
            if (ierr /= 0) return
      
            r_00 = wrap_opt_time_center_r_00(s, k)
            d_00 = wrap_d_00(s, k)
            Peos_00 = wrap_Peos_00(s, k)
            Cp_00 = wrap_Cp_00(s, k)
            T_00 = wrap_T_00(s, k)
            chiT_00 = wrap_chiT_00(s, k)
            chiRho_00 = wrap_chiRho_00(s, k)
            QQ_00 = chiT_00/(d_00*T_00*chiRho_00)
            lnT_00 = wrap_lnT_00(s,k)
      
            r_m1 = wrap_opt_time_center_r_m1(s, k)
            d_m1 = wrap_d_m1(s, k)
            Peos_m1 = wrap_Peos_m1(s, k)
            Cp_m1 = wrap_Cp_m1(s, k)
            T_m1 = wrap_T_m1(s, k)
            chiT_m1 = wrap_chiT_m1(s, k)
            chiRho_m1 = wrap_chiRho_m1(s, k)
            QQ_m1 = chiT_m1/(d_m1*T_m1*chiRho_m1)
            lnT_m1 = wrap_lnT_m1(s,k)
            QQ_div_Cp_face = alfa*QQ_00/Cp_00 + beta*QQ_m1/Cp_m1
            ! QQ units (g cm^-3 K)^-1 = g^-1 cm^3 K^-1
            ! Cp units erg g^-1 K^-1 = g cm^2 s^-2 g^-1 K^-1 = cm^2 s^-2 K^-1
            ! QQ/Cp units = (g^-1 cm^3 K^-1)/(cm^2 s^-2 K^-1)
            !  = g^-1 cm^3 K^-1 cm^-2 s^2 K
            !  = g^-1 cm s^2
            ! P units = erg cm^-3 = g cm^2 s^-2 cm^-3 = g cm^-1 s^-2
            ! QQ/Cp*P is unitless.
         
            Y1 = QQ_div_Cp_face*(Peos_m1 - Peos_00) - (lnT_m1 - lnT_00)
            ! Y1 unitless
         
            Y2 = 4d0*pi*pow2(r_00)*Hp_face*2d0/(1/d_00 + 1/d_m1)/dm_bar
            ! units = cm^2 cm / (cm^3 g^-1) / g
            !       = cm^2 cm cm^-3 g g^-1 = unitless
         
            Y_face = Y1*Y2 ! unitless
         
         else
         
            grad_ad_00 = wrap_grad_ad_00(s,k)
            grad_ad_m1 = wrap_grad_ad_m1(s,k)
            grad_ad_face = alfa*grad_ad_00 + beta*grad_ad_m1
            dlnT = wrap_lnT_m1(s,k) - wrap_lnT_00(s,k)
            dlnP = wrap_lnPeos_m1(s,k) - wrap_lnPeos_00(s,k)
            dlnT_dlnP = dlnT/dlnP
            if (is_bad(dlnT_dlnP%val)) then
               alt_Y_face = 0d0
            else if (s% use_Ledoux_criterion .and. s% calculate_Brunt_B) then
               ! gradL = grada + gradL_composition_term
               alt_Y_face = dlnT_dlnP - (grad_ad_face + s% gradL_composition_term(k))
            else
               alt_Y_face = dlnT_dlnP - grad_ad_face
            end if
            if (is_bad(alt_Y_face%val)) alt_Y_face = 0
            Y_face = alt_Y_face
            
         end if
         s% Y_face(k) = Y_face%val

      end function compute_Y_face
      
      
      function compute_PII_face(s, k, ierr) result(PII_face) ! ergs g^-1 K^-1 (like Cp)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: PII_face
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Cp_00, Cp_m1, Cp_face, Y_face
         real(dp) :: ALFAS_ALFA, alfa, beta
         include 'formats'
         ierr = 0
         if (k > s% nz) then
            PII_face = 0d0
            return
         end if
         if (k == 1 .or. s% TDC_alfa == 0d0 .or. &
               k > s% nz - s% TDC_num_innermost_cells_forced_nonturbulent) then
            PII_face = 0d0
            s% PII(k) = 0d0
            s% Y_face(k) = 0d0
            return
         end if
         Y_face = compute_Y_face(s, k, ierr)
         if (ierr /= 0) return
         Cp_00 = wrap_Cp_00(s, k)
         Cp_m1 = wrap_Cp_m1(s, k)
         call get_TDC_alfa_beta_face_weights(s, k, alfa, beta)
         Cp_face = alfa*Cp_00 + beta*Cp_m1 ! ergs g^-1 K^-1
         ALFAS_ALFA = x_ALFAS*s% TDC_alfa
         PII_face = ALFAS_ALFA*Cp_face*Y_face
         s% PII(k) = PII_face%val
         if (k == -2 .and. s% PII(k) < 0d0) then
            write(*,2) 's% PII(k)', k, s% PII(k)
            write(*,2) 'Cp_face', k, Cp_face%val
            write(*,2) 'Y_face', k, Y_face%val
            !write(*,2) 'PII_face%val', k, PII_face%val
            !write(*,2) 'T_rho_face%val', k, T_rho_face%val
            !write(*,2) '', k, 
            !write(*,2) '', k, 
            stop 'compute_PII_face'
         end if
      end function compute_PII_face
      
      
      function compute_d_v_div_r(s, k, ierr) result(d_v_div_r) ! s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: d_v_div_r
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: v_00, v_p1, r_00, r_p1
         include 'formats'
         ierr = 0
         v_00 = wrap_opt_time_center_v_00(s,k)
         v_p1 = wrap_opt_time_center_v_p1(s,k)
         r_00 = wrap_opt_time_center_r_00(s,k)
         r_p1 = wrap_opt_time_center_r_p1(s,k)
         if (r_p1%val == 0d0) r_p1 = 1d0
         d_v_div_r = v_00/r_00 - v_p1/r_p1 ! units s^-1
      end function compute_d_v_div_r
      
      
      function compute_Chi_cell(s, k, ierr) result(Chi_cell) ! eddy viscosity energy (Kuhfuss 1986) [erg]
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Chi_cell
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            w_rho2, r6_cell, d_v_div_r, Hp_cell, w_00, d_00, r_00, r_p1
         real(dp) :: f, ALFAM_ALFA
         include 'formats'
         ierr = 0
         ALFAM_ALFA = s% TDC_alfam*s% TDC_alfa
         if (ALFAM_ALFA == 0d0 .or. &
               k <= s% TDC_num_outermost_cells_forced_nonturbulent .or. &
               k > s% nz - s% TDC_num_innermost_cells_forced_nonturbulent) then
            Chi_cell = 0d0
         else
            Hp_cell = compute_Hp_cell(s, k, ierr)
            if (ierr /= 0) return
            d_v_div_r = compute_d_v_div_r(s, k, ierr)
            if (ierr /= 0) return
            w_00 = wrap_w_00(s,k)
            d_00 = wrap_d_00(s,k)
            f = (16d0/3d0)*pi*ALFAM_ALFA/s% dm(k)  
            w_rho2 = w_00*pow2(d_00)
            r_00 = wrap_opt_time_center_r_00(s,k)
            r_p1 = wrap_opt_time_center_r_p1(s,k)
            r6_cell = 0.5d0*(pow6(r_00) + pow6(r_p1))
            Chi_cell = f*w_rho2*r6_cell*d_v_div_r*Hp_cell
            ! units = g^-1 cm s^-1 g^2 cm^-6 cm^6 s^-1 cm
            !       = g cm^2 s^-2
            !       = erg
         end if
         s% Chi(k) = Chi_cell%val

      end function compute_Chi_cell
      
      
      function compute_dChi_dm_bar_face(s,k,ierr) result(dChi_dm_bar)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: dChi_dm_bar
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            d_v_div_r_00, w_rho2_00, f_00, Chi_00_div_r6_Hp, &
            d_v_div_r_m1, w_rho2_m1, f_m1, Chi_m1_div_r6_Hp, &
            Hp_face, r6_face, Chi_00, Chi_out
         real(dp) :: alfa, beta, ALFAM_ALFA
         if (k > 1 .and. .not. s% TDC_assume_HSE) then
            ! complexity needed to keep to block tridiagonal.
            ALFAM_ALFA = s% TDC_alfam*s% TDC_alfa
            call get_TDC_alfa_beta_face_weights(s, k, alfa, beta)
            d_v_div_r_00 = compute_d_v_div_r(s, k, ierr)
            if (ierr /= 0) return
            w_rho2_00 = wrap_w_00(s,k)*pow2(wrap_d_00(s,k))
            f_00 = (16d0/3d0)*pi*ALFAM_ALFA/s% dm(k)  
            Chi_00_div_r6_Hp = f_00*w_rho2_00*d_v_div_r_00
            d_v_div_r_m1 = shift_m1(compute_d_v_div_r(s, k-1, ierr))
            if (ierr /= 0) return
            w_rho2_m1 = wrap_w_m1(s,k)*pow2(wrap_d_m1(s,k))
            f_m1 = (16d0/3d0)*pi*ALFAM_ALFA/s% dm(k-1)  
            Chi_m1_div_r6_Hp = f_m1*w_rho2_m1*d_v_div_r_m1
            Hp_face = compute_Hp_face(s,k,ierr)
            if (ierr /= 0) return            
            r6_face = pow6(wrap_r_00(s,k))
            dChi_dm_bar = (Chi_m1_div_r6_Hp - Chi_00_div_r6_Hp)*r6_face*Hp_face
         else
            Chi_00 = compute_Chi_cell(s,k,ierr)
            if (k > 1) then
               Chi_out = shift_m1(compute_Chi_cell(s,k-1,ierr))
               if (ierr /= 0) return
            else
               Chi_out = 0d0
            end if
            dChi_dm_bar = (Chi_out - Chi_00)/s% dm_bar(k)
         end if
      end function compute_dChi_dm_bar_face

      
      function compute_Eq_cell(s, k, ierr) result(Eq_cell) ! erg g^-1 s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Eq_cell
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: d_v_div_r, Chi_cell
         include 'formats'
         ierr = 0
         if (k <= s% TDC_num_outermost_cells_forced_nonturbulent .or. &
             k > s% nz - s% TDC_num_innermost_cells_forced_nonturbulent) then
            Eq_cell = 0d0
         else
            Chi_cell = compute_Chi_cell(s,k,ierr)
            if (ierr /= 0) return
            d_v_div_r = compute_d_v_div_r(s, k, ierr)
            if (ierr /= 0) return
            Eq_cell = 4d0*pi*Chi_cell*d_v_div_r/s% dm(k) ! erg s^-1 g^-1
         end if
         s% Eq(k) = Eq_cell%val
      end function compute_Eq_cell


      function compute_Uq_face(s, k, ierr) result(Uq_face) ! cm s^-2, acceleration
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Uq_face
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: dChi_dm_bar, r_00
         include 'formats'
         ierr = 0         
         if (k <= s% TDC_num_outermost_cells_forced_nonturbulent .or. &
             k > s% nz - s% TDC_num_innermost_cells_forced_nonturbulent) then
            Uq_face = 0d0
         else
            r_00 = wrap_opt_time_center_r_00(s,k)
            dChi_dm_bar = compute_dChi_dm_bar_face(s,k,ierr)
            if (ierr /= 0) return
            Uq_face = 4d0*pi*dChi_dm_bar/r_00
         end if
         ! erg g^-1 cm^-1 = g cm^2 s^-2 g^-1 cm^-1 = cm s^-2, acceleration
         s% Uq(k) = Uq_face%val
      end function compute_Uq_face


      function compute_Source(s, k, ierr) result(Source) ! erg g^-1 s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Source
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            w_00, T_00, d_00, Peos_00, Cp_00, chiT_00, chiRho_00, QQ_00, &
            Hp_face_00, Hp_face_p1, PII_face_00, PII_face_p1, PII_div_Hp_cell, &
            grad_ad_00, fac
         include 'formats'
         ierr = 0
         w_00 = wrap_w_00(s, k)
         T_00 = wrap_T_00(s, k)                  
         d_00 = wrap_d_00(s, k)         
         Peos_00 = wrap_Peos_00(s, k)         
         Cp_00 = wrap_Cp_00(s, k)
         chiT_00 = wrap_chiT_00(s, k)
         chiRho_00 = wrap_chiRho_00(s, k)
         QQ_00 = chiT_00/(d_00*T_00*chiRho_00)
            
         Hp_face_00 = compute_Hp_face(s,k,ierr)
         if (ierr /= 0) return
         PII_face_00 = compute_PII_face(s, k, ierr)
         if (ierr /= 0) return
         
         if (k == s% nz) then
            PII_div_Hp_cell = PII_face_00/Hp_face_00
         else
            Hp_face_p1 = shift_p1(compute_Hp_face(s,k+1,ierr))
            if (ierr /= 0) return
            PII_face_p1 = shift_p1(compute_PII_face(s, k+1, ierr))
            if (ierr /= 0) return
            PII_div_Hp_cell = 0.5d0*(PII_face_00/Hp_face_00 + PII_face_p1/Hp_face_p1)
         end if
         
         !fac = w_00*T_00*(Peos_00*QQ_00/Cp_00)  ! create separate term just for debugging
         ! Peos_00*QQ_00/Cp_00 = grad_ad
         grad_ad_00 = wrap_grad_ad_00(s, k)
         fac = (w_00 + s% TDC_source_seed)*T_00*grad_ad_00

         Source = PII_div_Hp_cell*fac
         
         ! PII units same as Cp = erg g^-1 K^-1
         ! P*QQ/Cp is unitless (see Y_face)
         ! Source units = (erg g^-1 K^-1) cm^-1 cm s^-1 K
         !     = erg g^-1 s^-1
         
         if (k == -49) then
            write(*,2) 'Source%val', k, Source%val
            write(*,2) 'w_00%val', k, w_00%val
            write(*,2) 'grad_ad_00%val', k, grad_ad_00%val
            write(*,2) 'PII_face_00%val', k, PII_face_00%val
            write(*,2) 'PII_face_p1%val', k, PII_face_p1%val
            write(*,*)
            !stop 'compute_Source'
         end if
         s% SOURCE(k) = Source%val

      end function compute_Source


      function compute_D(s, k, ierr) result(D) ! erg g^-1 s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: D
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Hp_cell, w_00, dw3
         include 'formats'
         ierr = 0
         if (s% TDC_alfa == 0d0) then
            D = 0d0
         else
            Hp_cell = compute_Hp_cell(s,k,ierr)
            if (ierr /= 0) return
            w_00 = wrap_w_00(s,k)
            dw3 = pow3(w_00) - pow3(s% TDC_w_min_for_damping)
            D = (s% TDC_alfad*x_CEDE/s% TDC_alfa)*dw3/Hp_cell
            ! units cm^3 s^-3 cm^-1 = cm^2 s^-3 = erg g^-1 s^-1
         end if
         s% DAMP(k) = D%val
      end function compute_D


      function compute_Dr(s, k, ierr) result(Dr) ! erg g^-1 s^-1 = cm^2 s^-3
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Dr
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            w_00, T_00, d_00, Cp_00, kap_00, Hp_cell, POM2
         real(dp) :: gammar, alpha, POM
         include 'formats'
         ierr = 0
         alpha = s% TDC_alfa
         gammar = s% TDC_alfar*x_GAMMAR
         if (gammar == 0d0) then
            Dr = 0d0
            s% DAMPR(k) = 0d0
            return
         end if
         w_00 = wrap_w_00(s,k)
         T_00 = wrap_T_00(s,k)
         d_00 = wrap_d_00(s,k)
         Cp_00 = wrap_Cp_00(s,k)
         kap_00 = wrap_kap_00(s,k)
         Hp_cell = compute_Hp_cell(s,k,ierr)
         if (ierr /= 0) return
         POM = 4d0*boltz_sigma*(gammar/alpha)**2 ! erg cm^-2 K^-4 s^-1
         POM2 = pow3(T_00)/(pow2(d_00)*Cp_00*kap_00) 
            ! K^3 / ((g cm^-3)^2 (erg g^-1 K^-1) (cm^2 g^-1))
            ! K^3 / (cm^-4 erg K^-1) = K^4 cm^4 erg^-1
         Dr = POM*POM2*pow2(w_00)/pow2(Hp_cell)
         ! (erg cm^-2 K^-4 s^-1) (K^4 cm^4 erg^-1) cm^2 s^-2 cm^-2
         ! cm^2 s^-3 = erg g^-1 s^-1
         s% DAMPR(k) = Dr%val
      end function compute_Dr


      function compute_C(s, k, ierr) result(C) ! erg g^-1 s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: C
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Source, D, Dr
         if (k <= s% TDC_num_outermost_cells_forced_nonturbulent .or. &
             k > s% nz - s% TDC_num_innermost_cells_forced_nonturbulent) then
            s% SOURCE(k) = 0d0
            s% DAMP(k) = 0d0
            s% DAMPR(k) = 0d0
            s% COUPL(k) = 0d0
            C = 0d0
            return
         end if
         Source = compute_Source(s, k, ierr)
         if (ierr /= 0) return
         D = compute_D(s, k, ierr)
         if (ierr /= 0) return
         Dr = compute_Dr(s, k, ierr)
         if (ierr /= 0) return
         C = Source - D - Dr
         s% COUPL(k) = C%val
      end function compute_C


      function compute_L_face(s, k, ierr) result(L_face) ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: L_face
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Lr, Lc, Lt
         call compute_L_terms(s, k, L_face, Lr, Lc, Lt, ierr)
      end function compute_L_face


      subroutine compute_L_terms(s, k, L, Lr, Lc, Lt, ierr)
         type (star_info), pointer, intent(in) :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1), intent(out) :: L, Lr, Lc, Lt
         real(dp) :: L_val
         integer, intent(out) :: ierr         
         include 'formats'
         ierr = 0
         if (k > s% nz) then
            L = 0d0
            L%val = s% L_center
            Lr = 0d0
            Lc = 0d0
            Lt = 0d0
            return
         end if
         Lr = compute_Lr(s, k, ierr)
         if (ierr /= 0) return
         if (k == 1) then            
            Lc = 0d0
            Lt = 0d0
         else
            Lc = compute_Lc(s, k, ierr)
            if (ierr /= 0) return
            Lt = compute_Lt(s, k, ierr)
            if (ierr /= 0) return
         end if
         L = Lr + Lc + Lt
      end subroutine compute_L_terms


      function compute_Lr(s, k, ierr) result(Lr) ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lr
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: &
            r_00, area, T_00, T400, Erad, T_m1, T4m1, &
            kap_00, kap_m1, kap_face, diff_T4_div_kap, BW, BK
         real(dp) :: alfa
         include 'formats'
         ierr = 0
         if (k > s% nz) then
            Lr = s% L_center
         else
            r_00 = wrap_r_00(s,k) ! not time centered
            area = 4d0*pi*pow2(r_00)
            T_00 = wrap_T_00(s,k)
            T400 = pow4(T_00)
            if (k == 1) then ! Lr(1) proportional to Erad in cell(1)
               Erad = crad * T400
               Lr = s% TDC_Lsurf_factor * area * clight * Erad
               s% Lr(k) = Lr%val
               return
            end if
            T_m1 = wrap_T_m1(s,k)
            T4m1 = pow4(T_m1)            
            alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
            kap_00 = wrap_kap_00(s,k)
            kap_m1 = wrap_kap_m1(s,k)
            kap_face = alfa*kap_00 + (1d0 - alfa)*kap_m1
            diff_T4_div_kap = (T4m1 - T400)/kap_face

            if (s% TDC_use_Stellingwerf_Lr) then ! RSP style
               BW = log(T4m1/T400)
               if (abs(BW%val) > 1d-20) then
                  BK = log(kap_m1/kap_00)
                  if (abs(1d0 - BK%val/BW%val) > 1d-15 .and. abs(BW%val - BK%val) > 1d-15) then
                     diff_T4_div_kap = (T4m1/kap_m1 - T400/kap_00)/(1d0 - BK/BW)
                  end if
               end if
            end if
            Lr = -crad*clight/3d0*diff_T4_div_kap*pow2(area)/s% dm_bar(k)       
            ! units (erg cm^-3 K^-4) (cm s^-1) (K^4 cm^-2 g cm^4) g^-1 = erg s^-1  
         
            !s% xtra1_array(k) = s% T_start(k)
            !s% xtra2_array(k) = T4m1%val - T400%val
            !s% xtra3_array(k) = kap_face%val
            !s% xtra4_array(k) = diff_T4_div_kap%val
            !s% xtra5_array(k) = Lr%val/Lsun   
            !s% xtra6_array(k) = 1

         end if
         s% Lr(k) = Lr%val
      end function compute_Lr


      function compute_Lc(s, k, ierr) result(Lc) ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lc
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: Lc_w_face_factor
         Lc = compute_Lc_terms(s, k, Lc_w_face_factor, ierr)
         s% Lc(k) = Lc%val
      end function compute_Lc


      function compute_Lc_terms(s, k, Lc_w_face_factor, ierr) result(Lc)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lc, Lc_w_face_factor
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: r_00, area, &
            T_m1, T_00, d_m1, d_00, w_m1, w_00, T_rho_face, PII_face, w_face
         real(dp) :: ALFAC, ALFAS, alfa, beta
         include 'formats'
         ierr = 0
         if (k > s% nz .or. k == 1) then
            Lc = 0d0
            Lc_w_face_factor = 1
            return
         end if
         r_00 = wrap_r_00(s, k) ! not time centered
         area = 4d0*pi*pow2(r_00)
         T_m1 = wrap_T_m1(s, k)
         T_00 = wrap_T_00(s, k)         
         d_m1 = wrap_d_m1(s, k)
         d_00 = wrap_d_00(s, k)
         w_m1 = wrap_w_m1(s, k)
         w_00 = wrap_w_00(s, k)
         call get_TDC_alfa_beta_face_weights(s, k, alfa, beta)
         T_rho_face = alfa*T_00*d_00 + beta*T_m1*d_m1
         PII_face = compute_PII_face(s, k, ierr)
         w_face = alfa*w_00 + beta*w_m1
         ALFAC = x_ALFAC
         ALFAS = x_ALFAS
         Lc_w_face_factor = area*(ALFAC/ALFAS)*T_rho_face*PII_face
         ! units = cm^2 K g cm^-3 ergs g^-1 K^-1 = ergs cm^-1
         Lc = w_face*Lc_w_face_factor
         ! units = cm s^-1 ergs cm^-1 = ergs s^-1
         if (k == -458) then
            write(*,2) 'Lc%val', k, Lc%val
            write(*,2) 'w_face%val', k, w_face%val
            write(*,2) 'Lc_w_face_factor', k, Lc_w_face_factor%val
            write(*,2) 'PII_face%val', k, PII_face%val
            write(*,2) 'T_rho_face%val', k, T_rho_face%val
            !write(*,2) '', k, 
            !write(*,2) '', k, 
            stop 'compute_Lc_terms'
         end if
      end function compute_Lc_terms


      function compute_Lt(s, k, ierr) result(Lt) ! erg s^-1
         type (star_info), pointer :: s
         integer, intent(in) :: k
         type(auto_diff_real_star_order1) :: Lt
         integer, intent(out) :: ierr
         type(auto_diff_real_star_order1) :: r_00, area2, d_m1, d_00, &
            rho2_face, Hp_face, w_m1, w_00, w_face, etrb_m1, etrb_00
         real(dp) :: alpha_alpha_t, alfa, beta
         include 'formats'
         ierr = 0
         if (k > s% nz) then
            Lt = 0d0
            return
         end if 
         alpha_alpha_t = s% TDC_alfa*s% TDC_alfat
         if (k == 1 .or. alpha_alpha_t == 0d0) then
            Lt = 0d0
            s% Lt(k) = 0d0
            return
         end if
         r_00 = wrap_r_00(s,k) ! not time centered     
         area2 = (4d0*pi)**2*pow4(r_00)
         d_m1 = wrap_d_m1(s,k)
         d_00 = wrap_d_00(s,k)
         call get_TDC_alfa_beta_face_weights(s, k, alfa, beta)
         rho2_face = alfa*pow2(d_00) + beta*pow2(d_m1)
         w_m1 = wrap_w_m1(s,k)
         w_00 = wrap_w_00(s,k)
         w_face = alfa*w_00 + beta*w_m1
         etrb_m1 = wrap_etrb_m1(s,k)
         etrb_00 = wrap_etrb_00(s,k)
         Hp_face = compute_Hp_face(s,k,ierr)
         if (ierr /= 0) return         
         ! Ft = - alpha_t * rho_face * alpha * Hp_face * w_face * detrb/dr (thesis eqn 2.44)
         ! replace dr by dm_bar/(area*rho_face)
         ! Ft = - alpha_alpha_t * rho_face * Hp_face * w_face * (area*rho_face) * detrb/dm_bar
         ! Lt = area * Ft
         ! Lt = -alpha_alpha_t * (area*rho_face)**2 * Hp_face * w_face * (etrb(k-1) - etrb(k))/dm_bar
         Lt = - alpha_alpha_t * area2 * rho2_face * Hp_face * w_face * (etrb_m1 - etrb_00) / s% dm_bar(k)  
         ! this is slightly rewritten from the RSP form to use the TDC etrb variable
         ! units = (cm^4) (g^2 cm^-6) (cm) (cm s^-1) (ergs g^-1) g^-1 = erg s^-1
         s% Lt(k) = Lt%val      
      end function compute_Lt


      subroutine set_etrb_start_vars(s, ierr)
         type (star_info), pointer :: s
         integer, intent(out) :: ierr         
         integer :: k, op_err
         type(auto_diff_real_star_order1) :: Y_face, Lt
         include 'formats'
         ierr = 0
         do k=1,s%nz
            Y_face = compute_Y_face(s, k, ierr)
            if (ierr /= 0) return
            s% Y_face_start(k) = Y_face%val  
            Lt = compute_Lt(s, k, ierr)
            if (ierr /= 0) return
            s% Lt_start(k) = Lt%val  
            s% w_start(k) = s% w(k)
         end do         
      end subroutine set_etrb_start_vars
      
      
      subroutine reset_etrb_using_L(s, ierr)
         use star_utils, only: store_etrb_in_xh
         type (star_info), pointer :: s
         integer, intent(out) :: ierr   
         integer :: k, nz, j, k_maxerr
         real(dp) :: Lc_val, w_00, maxerr, dlnP, dlnT, gradT_actual, &
            super_ad_actual, super_ad_expected
         type(auto_diff_real_star_order1) :: &
            Lc_w_face_factor, L, Lr, Lc, Lt, Y_face
         real(dp), allocatable :: w_face(:), target_Lc(:)
         real(dp) :: alfa, beta
         real(dp), parameter :: atol = 1-6d0, rtol = 1d-9
         logical, parameter :: dbg = .false.
         include 'formats'
         ierr = 0
         if (s% TDC_alfa == 0d0) return ! no convection
         
         ! s% dt hasn't been set, so don't check it against TDC limits
         !call get_TDC_frac(s, alfa, beta)
         !s% using_TDC = (alfa > 0d0)
         !write(*,*) 'get_TDC_frac alfa TDC_flag', alfa, s% TDC_flag
         !stop 'reset_etrb_using_L'
         !if (.not. s% using_TDC) return
         
         nz = s% nz
         if (s% have_previous_conv_vel) then
            write(*,*) 'initial w_face set using conv_vel from file or from MLT'
         else
            write(*,*) 'need conv_vel from file or from MLT in order to set initial w_face'
            stop 'reset_etrb_using_L'
         end if
         allocate(w_face(nz), target_Lc(nz))

         if (.false.) then
            write(*,*) 'Compare standard and alternative forms of Y_face'
            write(*,'(a8,99(a18))') 'k', 'Y_face', 'old_Y_face/Y_face', 'new_Y_face/Y_face'
            do k=2,nz
               Y_face = compute_Y_face(s, k, ierr)
            end do
            stop 'reset_etrb_using_L'
         end if

         do k=1, nz
            Lr = compute_Lr(s, k, ierr)
            if (ierr /= 0) stop 'failed in compute_Lr'
            Lc = compute_Lc_terms(s, k, Lc_w_face_factor, ierr)
            if (ierr /= 0) stop 'failed in compute_Lc_terms'
            target_Lc(k) = s% L(k) - Lr%val
            Lc_val = target_Lc(k) ! assume Lt = 0 for this
            if (abs(Lc_w_face_factor%val) < 1d-20) then
               w_face(k) = 0d0
            else ! filter out negative guesses for w_face
               w_face(k) = max(0d0, Lc_val/Lc_w_face_factor%val)
            end if
            if (is_bad(w_face(k))) then
               write(*,2) 'bad w_face', k, w_face(k)
               stop 'reset_etrb_using_L'
               w_face(k) = 0d0
            end if
            ! CAUTION: using conv_vel is dangerous.
            !  better if can stick to L.  but problems with that too.
            if ((w_face(k) == 0d0 .and. s% conv_vel(k) > 0d0) .or. &
                (w_face(k) > 0d0 .and. s% conv_vel(k) == 0d0)) &
               w_face(k) = s% conv_vel(k)
         end do
         
         do k=1, nz
            if (k < nz) then
               w_00 = 0.5d0*(w_face(k) + w_face(k+1))
            else ! w_center = 0
               w_00 = 0.5d0*w_face(k)
            end if
            s% w(k) = w_00
            s% xh(s% i_w,k) = s% w(k)
         end do
         
         if (dbg) stop 'reset_etrb_using_L'
      end subroutine reset_etrb_using_L


      end module hydro_tdc

