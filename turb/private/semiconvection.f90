! ***********************************************************************
!
!   Copyright (C) 2010-2021  The MESA Team
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


module semiconvection

implicit none

private
public :: calc_semiconvection

contains

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
   subroutine calc_semiconvection(L, Lambda, m, T, P, Pr, beta, opacity, rho, alpha_semiconvection, &
                                 semiconvection_option, cgrav, Cp, gradr, grada, gradL, &
                                 gradL_composition_term, &
                                 gradT, Y_face, conv_vel, D, mixing_type, ierr) ! Langer 1983 & 1985
      use const_def
      use num_lib
      use utils_lib
      use auto_diff
      type(auto_diff_real_star_order1), intent(in) :: L, Lambda, T, P, Pr, beta, opacity, rho
      type(auto_diff_real_star_order1), intent(in) :: Cp, gradr, grada, gradL
      character(len=*), intent(in) :: semiconvection_option
      real(dp), intent(in) :: alpha_semiconvection, cgrav, gradL_composition_term, m
      type(auto_diff_real_star_order1), intent(out) :: gradT, Y_face, conv_vel, D
      integer, intent(out) :: mixing_type, ierr

      type(auto_diff_real_star_order1) :: bc, LG, &
         radiative_conductivity, a0, a1, a2, a3, a4, a5, a6, a, &
         b1, b2, b3, b4, b5, b6, b7, b, div, bsq    
      real(dp) :: alpha
      include 'formats'

      ! Pre-compute common pieces
      radiative_conductivity = &
         (4d0/3d0*crad*clight)*pow3(T)/(opacity*rho) ! erg / (K cm sec)
      D = alpha_semiconvection*radiative_conductivity/(6d0*Cp*rho) &
            *(gradr - grada)/(gradL - gradr)
      if (D%val <= 0) return         
      conv_vel = 3d0*D/Lambda             

      if (semiconvection_option == 'Langer_85 mixing; gradT = gradr') then
         gradT = gradr
         Y_face = gradT - grada
         mixing_type = semiconvective_mixing
         return
      else if (semiconvection_option == 'Langer_85') then
      !            Solve[{
      !                  L/Lrad - Lsc/Lrad - 1 == 0, 
      !                  Lrad == grad LG, 
      !                  gradMu == (4 - 3*beta)/beta*gradL_composition_term,
      !                  Lsc/Lrad == alpha (grad - gradA)/(2 grad (gradL - grad))
      !                              (grad - gradA - (beta (8 - 3 beta))/bc gradMu)}, 
      !                  grad, {Lsc, Lrad, gradMu}] // Simplify
         alpha = min(1d0, alpha_semiconvection)
         bc = 32d0 - 24d0*beta - beta*beta            
         LG = (16d0*pi*clight*m*cgrav*Pr)/(P*opacity)            
         a0 = alpha*gradL_composition_term*LG            
         a1 = -2d0*bc*L            
         a2 = 2d0*alpha*bc*grada*LG            
         a3 = -2d0*bc*gradL*LG            
         a4 = 32d0*a0            
         a5 = -36d0*beta*a0            
         a6 = 9d0*beta*beta*a0            
         a = a1 + a2 + a3 + a4 + a5 + a6                           
         b1 = 32d0 - 36d0*beta + 9d0*beta*beta            
         b2 = b1*a0            
         b3 = -2d0*gradL*L + alpha*grada*grada*LG            
         b4 = (-alpha*gradA + gradL)*LG            
         b5 = -b2 + 2d0*bc*(L + b4)            
         b6 = b2*grada + bc*b3            
         b7 = -4d0*(alpha - 2d0)*bc*LG*b6            
         b = b7 + b5*b5            
         div = 2d0*(alpha - 2d0)*bc*LG
         bsq = sqrt(b)
         gradT = (a + bsq)/div
         Y_face = gradT - grada
         conv_vel = 3d0*D/Lambda             
         mixing_type = semiconvective_mixing
      else
         write(*,*) 'turb: unknown values for semiconvection_option ' // &
            trim(semiconvection_option)
         ierr = -1
         return
      end if

   end subroutine calc_semiconvection

end module semiconvection
