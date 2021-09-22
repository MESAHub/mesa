module turb
   use const_def
   use num_lib
   use utils_lib
   use auto_diff

   implicit none

   private
   public :: set_thermohaline, set_mlt, set_tdc, set_semiconvection

   contains

   subroutine set_thermohaline(thermohaline_option, Lambda, grada, gradr, T, opacity, rho, Cp, gradL_composition_term, &
                              iso, XH1, thermohaline_coeff, &
                              D, gradT, Y_face, conv_vel, mixing_type, ierr)
      use thermohaline
      character(len=*), intent(in) :: thermohaline_option
      type(auto_diff_real_star_order1), intent(in) :: Lambda, grada, gradr, T, opacity, rho, Cp
      real(dp), intent(in) :: gradL_composition_term, XH1, thermohaline_coeff
      integer, intent(in) :: iso

      type(auto_diff_real_star_order1), intent(out) :: gradT, Y_face, conv_vel, D
      integer, intent(out) :: mixing_type, ierr

      real(dp) :: D_thrm

      call get_D_thermohaline(&
         thermohaline_option, grada%val, gradr%val, T%val, opacity%val, rho%val, &
         Cp%val, gradL_composition_term, &
         iso, XH1, thermohaline_coeff, D_thrm, ierr)

      D = D_thrm
      gradT = gradr
      Y_face = gradT - grada
      conv_vel = 3d0*D/Lambda
      mixing_type = thermohaline_mixing 
   end subroutine set_thermohaline

   subroutine set_TDC( &
            conv_vel_start, mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, cgrav, m, report, &
            mixing_type, scale, L, r, P, T, rho, dV, Cp, opacity, &
            scale_height, gradL, grada, conv_vel, D, Y_face, gradT, tdc_num_iters, ierr)
      use tdc
      real(dp), intent(in) :: conv_vel_start, mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, cgrav, m, scale
      type(auto_diff_real_star_order1), intent(in) :: &
         L, r, P, T, rho, dV, Cp, opacity, scale_height, gradL, grada
      logical, intent(in) :: report
      type(auto_diff_real_star_order1),intent(out) :: conv_vel, Y_face, gradT, D
      integer, intent(out) :: tdc_num_iters, mixing_type, ierr
      include 'formats'
      call get_TDC_solution( &
         conv_vel_start, mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, cgrav, m, report, &
         mixing_type, scale, L, r, P, T, rho, dV, Cp, opacity, &
         scale_height, gradL, grada, conv_vel, Y_face, tdc_num_iters, ierr)
      if (ierr /= 0) then
         write(*,*) 'get_TDC_solution failed in set_TDC'
         write(*,*) 'Repeating call with reporting on.'
         call get_TDC_solution( &
            conv_vel_start, mixing_length_alpha, alpha_TDC_DAMP, alpha_TDC_DAMPR, alpha_TDC_PtdVdt, dt, cgrav, m, .true., &
            mixing_type, scale, L, r, P, T, rho, dV, Cp, opacity, &
            scale_height, gradL, grada, conv_vel, Y_face, tdc_num_iters, ierr)
      end if

      gradT = Y_face + gradL
      D = conv_vel*scale_height*mixing_length_alpha/3d0     ! diffusion coefficient [cm^2/sec]

   end subroutine set_TDC

   subroutine set_semiconvection(L, Lambda, m, T, P, Pr, beta, opacity, rho, alpha_semiconvection, &
                                 semiconvection_option, cgrav, Cp, gradr, grada, gradL, &
                                 gradL_composition_term, &
                                 gradT, Y_face, conv_vel, D, mixing_type, ierr)
      use semiconvection
      type(auto_diff_real_star_order1), intent(in) :: L, Lambda, T, P, Pr, beta, opacity, rho
      type(auto_diff_real_star_order1), intent(in) :: Cp, gradr, grada, gradL
      character(len=*), intent(in) :: semiconvection_option
      real(dp), intent(in) :: alpha_semiconvection, cgrav, gradL_composition_term, m
      type(auto_diff_real_star_order1), intent(out) :: gradT, Y_face, conv_vel, D
      integer, intent(out) :: mixing_type, ierr

      call calc_semiconvection(L, Lambda, m, T, P, Pr, beta, opacity, rho, alpha_semiconvection, &
                                 semiconvection_option, cgrav, Cp, gradr, grada, gradL, &
                                 gradL_composition_term, &
                                 gradT, Y_face, conv_vel, D, mixing_type, ierr)
   end subroutine set_semiconvection

   subroutine set_MLT(MLT_option, mixing_length_alpha, Henyey_MLT_nu_param, Henyey_MLT_y_param, &
                     chiT, chiRho, Cp, grav, Lambda, rho, P, T, opacity, &
                     gradr, grada, gradL, &
                     Gamma, gradT, Y_face, conv_vel, D, mixing_type, ierr)
      use mlt
      type(auto_diff_real_star_order1), intent(in) :: chiT, chiRho, Cp, grav, Lambda, rho, P, T, opacity, gradr, grada, gradL
      character(len=*), intent(in) :: MLT_option
      real(dp), intent(in) :: mixing_length_alpha, Henyey_MLT_nu_param, Henyey_MLT_y_param
      type(auto_diff_real_star_order1), intent(out) :: Gamma, gradT, Y_face, conv_vel, D
      integer, intent(out) :: mixing_type, ierr

      call calc_MLT(MLT_option, mixing_length_alpha, Henyey_MLT_nu_param, Henyey_MLT_y_param, &
                     chiT, chiRho, Cp, grav, Lambda, rho, P, T, opacity, &
                     gradr, grada, gradL, &
                     Gamma, gradT, Y_face, conv_vel, D, mixing_type, ierr)
   end subroutine set_MLT

end module turb
