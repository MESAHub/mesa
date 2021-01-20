! ***********************************************************************
!
!   Copyright (C) 2011  Bill Paxton
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
 
      module other_mlt

      ! consult star/other/README for general usage instructions
      ! control name: use_other_mlt = .true.
      ! procedure pointer: s% other_mlt => my_routine


      use star_def

      implicit none
      
            
      contains
      
      
      subroutine null_other_mlt( &
            id, k, cgrav, m, mstar, r, L, X, &            
            T_face, rho_face, P_face, &
            chiRho_face, chiT_face, &
            Cp_face, opacity_face, grada_face, &            
            alfa, beta, & ! f_face = alfa*f_00 + beta*f_m1
            T_00, T_m1, rho_00, rho_m1, P_00, P_m1, &
            chiRho_for_partials_00, chiT_for_partials_00, &
            chiRho_for_partials_m1, chiT_for_partials_m1, &
            chiRho_00, d_chiRho_00_dlnd, d_chiRho_00_dlnT, &
            chiRho_m1, d_chiRho_m1_dlnd, d_chiRho_m1_dlnT, &
            chiT_00, d_chiT_00_dlnd, d_chiT_00_dlnT, &
            chiT_m1, d_chiT_m1_dlnd, d_chiT_m1_dlnT, &
            Cp_00, d_Cp_00_dlnd, d_Cp_00_dlnT, &
            Cp_m1, d_Cp_m1_dlnd, d_Cp_m1_dlnT, &
            opacity_00, d_opacity_00_dlnd, d_opacity_00_dlnT, &
            opacity_m1, d_opacity_m1_dlnd, d_opacity_m1_dlnT, &
            grada_00, d_grada_00_dlnd, d_grada_00_dlnT, &
            grada_m1, d_grada_m1_dlnd, d_grada_m1_dlnT, &            
            gradr_factor, d_gradr_factor_dw, gradL_composition_term, &
            alpha_semiconvection, semiconvection_option, &
            thermohaline_coeff, thermohaline_option, &
            dominant_iso_for_thermohaline, &
            mixing_length_alpha, alt_scale_height, remove_small_D_limit, &
            MLT_option, Henyey_y_param, Henyey_nu_param, &
            normal_mlt_gradT_factor, &
            max_conv_vel, dt, tau, just_gradr, &
            mixing_type, mlt_basics, mlt_partials1, ierr)
         
! UNCOMMENT THIS
         !use star_lib, only: star_mlt_eval
         use eos_lib, only: Radiation_Pressure
         
         integer, intent(in) :: id ! id for star         
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         real(dp), intent(in) :: &
            cgrav, m, mstar, r, L, X, &            
            T_face, rho_face, P_face, &
            chiRho_face, chiT_face, &
            Cp_face, opacity_face, grada_face, &            
            alfa, beta, &
            T_00, T_m1, rho_00, rho_m1, P_00, P_m1, &
            chiRho_for_partials_00, chiT_for_partials_00, &
            chiRho_for_partials_m1, chiT_for_partials_m1, &
            chiRho_00, d_chiRho_00_dlnd, d_chiRho_00_dlnT, &
            chiRho_m1, d_chiRho_m1_dlnd, d_chiRho_m1_dlnT, &
            chiT_00, d_chiT_00_dlnd, d_chiT_00_dlnT, &
            chiT_m1, d_chiT_m1_dlnd, d_chiT_m1_dlnT, &
            Cp_00, d_Cp_00_dlnd, d_Cp_00_dlnT, &
            Cp_m1, d_Cp_m1_dlnd, d_Cp_m1_dlnT, &
            opacity_00, d_opacity_00_dlnd, d_opacity_00_dlnT, &
            opacity_m1, d_opacity_m1_dlnd, d_opacity_m1_dlnT, &
            grada_00, d_grada_00_dlnd, d_grada_00_dlnT, &
            grada_m1, d_grada_m1_dlnd, d_grada_m1_dlnT, &
            gradr_factor, d_gradr_factor_dw, gradL_composition_term, &
            alpha_semiconvection, thermohaline_coeff, mixing_length_alpha, &
            Henyey_y_param, Henyey_nu_param, &
            max_conv_vel, dt, tau, remove_small_D_limit, &
            normal_mlt_gradT_factor
         logical, intent(in) :: alt_scale_height
         character (len=*), intent(in) :: thermohaline_option, MLT_option, semiconvection_option
         integer, intent(in) :: dominant_iso_for_thermohaline
         logical, intent(in) :: just_gradr
         integer, intent(out) :: mixing_type
         real(dp), intent(inout) :: mlt_basics(:) ! (num_mlt_results)
         real(dp), intent(inout), pointer :: mlt_partials1(:) ! =(num_mlt_partials, num_mlt_results)
         integer, intent(out) :: ierr
         
         type (star_info), pointer :: s
         real(dp), pointer :: mlt_partials(:,:)
         real(dp) :: alfa_blend, beta_blend, factor, &
            T_start_face, P_start_face, Prad_start_face, min_ratio, max_ratio
         integer :: j
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         mixing_type = -1
         
! UNCOMMENT THIS
         !call star_mlt_eval(  &
            !id, k, cgrav, m, mstar, r, L, X, &            
            !T_face, rho_face, P_face, &
            !chiRho_face, chiT_face, &
            !Cp_face, opacity_face, grada_face, &            
            !alfa, beta, & ! f_face = alfa*f_00 + beta*f_m1
            !T_00, T_m1, rho_00, rho_m1, P_00, P_m1, &
            !chiRho_for_partials_00, chiT_for_partials_00, &
            !chiRho_for_partials_m1, chiT_for_partials_m1, &
            !chiRho_00, d_chiRho_00_dlnd, d_chiRho_00_dlnT, &
            !chiRho_m1, d_chiRho_m1_dlnd, d_chiRho_m1_dlnT, &
            !chiT_00, d_chiT_00_dlnd, d_chiT_00_dlnT, &
            !chiT_m1, d_chiT_m1_dlnd, d_chiT_m1_dlnT, &
            !Cp_00, d_Cp_00_dlnd, d_Cp_00_dlnT, &
            !Cp_m1, d_Cp_m1_dlnd, d_Cp_m1_dlnT, &
            !opacity_00, d_opacity_00_dlnd, d_opacity_00_dlnT, &
            !opacity_m1, d_opacity_m1_dlnd, d_opacity_m1_dlnT, &
            !grada_00, d_grada_00_dlnd, d_grada_00_dlnT, &
            !grada_m1, d_grada_m1_dlnd, d_grada_m1_dlnT, &            
            !gradr_factor, d_gradr_factor_dw, gradL_composition_term, &
            !alpha_semiconvection, semiconvection_option, &
            !thermohaline_coeff, thermohaline_option, &
            !dominant_iso_for_thermohaline, &
            !mixing_length_alpha, alt_scale_height, remove_small_D_limit, &
            !MLT_option, Henyey_y_param, Henyey_nu_param, &
            !normal_mlt_gradT_factor, &
            !max_conv_vel, dt, tau, just_gradr, &
            !mixing_type, mlt_basics, mlt_partials1, ierr)
         
         ! see star_data/public/star_data_def.inc for lists of mlt_basics and mlt_partials
         
         
         mlt_partials(1:num_mlt_partials,1:num_mlt_results) => &
            mlt_partials1(1:num_mlt_partials*num_mlt_results)

         ! factor = multiplier for gradr
         factor = 0.8d0 ! calculate according to local conditions
         
         ! blend from standard gradT to gradT = factor*gradr depending on ratio Prad_face/P_face
         ! alfa_blend = fraction factor*gradr, beta_blend = fraction standard gradT
         min_ratio = 0.2
         max_ratio = 0.8
         ! use start of step values for T and P so that the value of alfa_blend is constant during iterations
         ! that saves us from dealing with partials of alfa_blend
         if (k > 1) then
            T_start_face = alfa*s% T_start(k) + beta*s% T_start(k-1)
            P_start_face = alfa*s% P_start(k) + beta*s% P_start(k-1)
         else
            T_start_face = s% T_start(k)
            P_start_face = s% P_start(k)
         end if
         Prad_start_face = Radiation_Pressure(T_start_face)
         if (Prad_start_face/P_start_face <= min_ratio) then
            alfa_blend = 0d0
         else if (Prad_start_face/P_start_face >= max_ratio) then
            alfa_blend = 1d0
         else
            alfa_blend = (Prad_start_face/P_start_face - min_ratio)/(max_ratio - min_ratio)
         end if
         beta_blend = 1d0 - alfa_blend
                  
         ! combine factor*gradr with the standard gradT in mlt_basics and mlt_partials  
         mlt_basics(mlt_gradT) = alfa_blend*factor*mlt_basics(mlt_gradr) + beta_blend*mlt_basics(mlt_gradT)
         do j=1,num_mlt_partials
            mlt_partials(j,mlt_gradT) = alfa_blend*factor*mlt_partials(j,mlt_gradr) + beta_blend*mlt_partials(j,mlt_gradT)
         end do
         
         ! the caller will unpack all of this to set things like s% gradT(k) and s% d_gradT_dlnT00(k)
         
      end subroutine null_other_mlt


      end module other_mlt
      
      
      
      
