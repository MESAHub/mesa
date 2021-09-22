module turb
   use const_def
   use num_lib
   use utils_lib
   use auto_diff

   implicit none

   private
   public :: set_thermohaline, set_mlt, set_tdc, set_semiconvection

   contains

   !> Computes the diffusivity of thermohaline mixing when the
   !! thermal gradient is stable and the composition gradient is unstable.
   !!
   !! @param thermohaline_option A string specifying which thermohaline prescription to use.
   !! @param grada Adiabatic gradient dlnT/dlnP
   !! @param gradr Radiative temperature gradient dlnT/dlnP, equals the actual gradient because there's no convection
   !! @param T Temperature
   !! @param opacity opacity
   !! @param rho Density
   !! @param Cp Heat capacity at constant pressure
   !! @param gradL_composition_term dlnMu/dlnP where Mu is the mean molecular weight.
   !! @param iso The index of the species that drives thermohaline mixing.
   !! @param XH1 Mass fraction of H1.
   !! @param thermohaline_coeff Free parameter multiplying the thermohaline diffusivity.
   !! @param D_thrm Output, diffusivity.
   !! @param ierr Output, error index.
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

   !> Calculates the outputs of semiconvective mixing theory.
   !!
   !! @param L Luminosity across a face (erg/s).
   !! @param Lambda The mixing length (cm).
   !! @param m Mass inside the face (g).
   !! @param T temperature (K).
   !! @param P pressure (erg/cm^3).
   !! @param Pr radiation pressure (erg/cm^3).
   !! @param beta ratio of gas pressure to radiation pressure.
   !! @param opacity opacity (cm^2/g).
   !! @param rho density (g/cm^3).
   !! @param alpha_semiconvection The semiconvective alpha parameter.
   !! @param semiconvection_option A string specifying which semiconvection theory to use. Currently supported are 'Langer_85 mixing; gradT = gradr' and 'Langer_85'.
   !! @param cgrav gravitational constant (erg*cm/g^2).
   !! @param Cp Specific heat at constant pressure (erg/g/K).
   !! @param gradr The radiative temperature gradient dlnT/dlnP_{rad}
   !! @param grada The adiabatic temperature gradient dlnT/dlnP|s
   !! @param gradL The Ledoux temperature gradient dlnT/dlnP
   !! @param gradL_composition_term The contribution of composition gradients to the Ledoux temperature gradient.
   !! @param gradT The temperature gradient dlnT/dlnP (output).
   !! @param Y_face The superadiabaticity (dlnT/dlnP - grada, output).
   !! @param conv_vel The convection speed (cm/s).
   !! @param D The chemical diffusion coefficient (cm^2/s).
   !! @param mixing_type Set to semiconvective if convection operates (output).
   !! @param ierr Tracks errors (output).
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

   !> Calculates the outputs of convective mixing length theory.
   !!
   !! @param MLT_option A string specifying which MLT option to use. Currently supported are Cox, Henyey, ML1, ML2, Mihalas. Note that 'TDC' is also a valid input and will return the Cox result. This is for use when falling back from TDC -> MLT, as Cox is the most-similar prescription to TDC.
   !! @param mixing_length_alpha The mixing length parameter.
   !! @param Henyey_MLT_nu_param The nu parameter in Henyey's MLT prescription.
   !! @param Henyey_MLT_y_param The y parameter in Henyey's MLT prescription.
   !! @param chiT dlnP/dlnT|rho
   !! @param chiRho dlnP/dlnRho|T
   !! @param Cp Specific heat at constant pressure (erg/g/K).
   !! @param grav The acceleration due to gravity (cm/s^2).
   !! @param Lambda The mixing length (cm).
   !! @param rho density (g/cm^3).
   !! @param T temperature (K).
   !! @param opacity opacity (cm^2/g)
   !! @param gradr The radiative temperature gradient dlnT/dlnP_{rad}
   !! @param grada The adiabatic temperature gradient dlnT/dlnP|s
   !! @param gradL The Ledoux temperature gradient dlnT/dlnP
   !! @param Gamma The convective efficiency parameter (output).
   !! @param gradT The temperature gradient dlnT/dlnP (output).
   !! @param Y_face The superadiabaticity (dlnT/dlnP - grada, output).
   !! @param conv_vel The convection speed (cm/s).
   !! @param D The chemical diffusion coefficient (cm^2/s).
   !! @param mixing_type Set to convective if convection operates (output).
   !! @param ierr Tracks errors (output).
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
