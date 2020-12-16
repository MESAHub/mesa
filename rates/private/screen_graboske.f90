! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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

      module screen_graboske
      use const_def
      use rates_def
      use math_lib
      
      implicit none

      contains
      
      subroutine graboske_init_z_info(zg1, zg2, zg3, zg4, z1, z2, ierr)
         use utils_lib, only: realloc_double
         !..compute and store things that only depend on z1 and z2
         real(dp), intent(out) :: zg1, zg2, zg3, zg4
         real(dp), intent(in) :: z1, z2
         integer, intent(out) :: ierr
         
         integer :: iz1, iz2, i, new_sz
         real(dp) :: z1_13, z2_13, z12_13, z1_23, z2_23, z12_23, b
         real(dp), parameter :: x13 = 1.0d0/3.0d0 
         real(dp), parameter :: x23 = 2.0d0/3.0d0
         real(dp), parameter :: x43 = 4.0d0/3.0d0
         real(dp), parameter :: x53 = 5.0d0/3.0d0
         real(dp), parameter :: x512  = 5.0d0/12.0d0
         
         logical, parameter :: debug = .false.
         
         include 'formats.dek'
         
         ierr = 0
         
         iz1 = int(z1)
         iz2 = int(z2)
         if (iz1 <= 0 .or. iz2 <= 0) then
            zg1 = 0d0
            zg2 = 0d0
            zg3 = 0d0
            zg4 = 0d0
            return
         end if
         
         if (iz1 > iz2) then ! switch so that iz1 <= iz2
            i = iz1; iz1 = iz2; iz2 = i
         end if
         
         if (iz1 <= num_one_thirds) then
            z1_13 = one_third_power(iz1)
         else
            z1_13 = pow(z1,x13)
         end if
         if (iz2 <= num_one_thirds) then
            z2_13 = one_third_power(iz2)
         else
            z2_13 = pow(z2,x13)
         end if
         if (iz1+iz2 <= num_one_thirds) then
            z12_13 = one_third_power(iz1+iz2)
         else
            z12_13 = pow(z1+z2,x13)
         end if
         z1_23 = z1_13*z1_13
         z2_23 = z2_13*z2_13
         z12_23 = z12_13*z12_13

         zg1 = (z1+z2)*z12_23  - z1*z1_23 - z2*z2_23
         zg2 = ((z1+z2)*z1_13 - z1*z1_13 - z2*z2_13)
         zg3 = (z12_23 - z1_23 - z2_23)

         b = 0.860d0
         zg4 = pow(z1+z2,1d0+b) - pow(z1,1d0+b) - pow(z2,1d0+b)
         if (debug) write(*,*) 'graboske_et_al_screening initization call'
      
      end subroutine graboske_init_z_info

      subroutine graboske_et_al_screening( &
               sc, zg1, zg2, zg3, zg4, theta_e, cache,  &
               a1, z1, a2, z2, scor, dscor_dT, dscor_dRho, ierr)
         use utils_lib, only:is_bad
         
         ! DeWitt, Graboske, Cooper, "Screening Factors for Nuclear Reactions. 
         !    I. General Theory", ApJ, 181:439-456, 1973.
         ! Graboske, DeWitt, Grossman, Cooper, "Screening Factors for Nuclear Reactions. 
         !    II. Intermediate Screening and Astrophysical Applications", ApJ, 181:457-474, 1973.

         type (Screen_Info), pointer :: sc
         real(dp), intent(in) :: zg1, zg2, zg3, zg4
         real(dp), intent(in) :: theta_e
         real(dp), pointer :: cache(:, :, :)
         real(dp), intent(in) :: a1, z1, a2, z2
         real(dp), intent(out) :: scor, dscor_dT, dscor_dRho
         integer, intent(out) :: ierr

         ! z1 and z2 are ion charge numbers
         ! T and rho are temperature and density in cgs units
         ! y = x / chem_A
         ! abar = 1 / sum(y(:))
         ! zbar = abar * sum(z(:) * y(:))
         ! z2bar = abar * sum(z(:)**2 * y(:))
         ! z1pt58bar = abar * sum(z(:)**1.58d0 * y(:))
         ! z1pt58bar is used in calculating intermediate screening

         real(dp) ::  &
               T, rho, abar, zbar, z2bar, ztilda, zg2_screen, zg3_screen, &
               H120, H120a, dH120_dT, dH120_dRho,  &
               Lambda0, Lambda12, b, k_b, eta_b, zeta_b, &
               d_Lambda0_dT, d_Lambda0_dRho, z1pt58bar, zbar13
         real(dp) :: alfa, H120_weak, H120_intermediate, H120_strong,  &
               dH120_dT_weak, dH120_dRho_weak, &
               dH120_dT_intermediate, dH120_dRho_intermediate, &
               dH120_dT_strong, dH120_dRho_strong
     
         real(dp), parameter :: H120_max = 300d0         
         real(dp), parameter :: Lam_1  = 0.1d0
         real(dp), parameter :: Lam_2  = 0.125d0
         real(dp), parameter :: Lam_3  = 2d0
         real(dp), parameter :: Lam_4  = 2.15d0
         real(dp), parameter :: Lam_5  = 4.85d0
         real(dp), parameter :: Lam_6  = 5d0
         
         integer :: max_z_for_cache, i, iz1, iz2
         
         logical, parameter :: use_cache = .true.
         logical :: debug
         
         include 'formats.dek'
         
         debug = .false.
         ierr = 0
         scor = 1
         dscor_dT = 0
         dscor_dRho = 0
         
         iz1 = int(z1)
         iz2 = int(z2)
         if (iz1 <= 0 .or. iz2 <= 0) return
         
         if (iz1 > iz2) then ! switch so that iz1 <= iz2
            i = iz1; iz1 = iz2; iz2 = i
         end if
         
         !debug = (iz1 == 1 .and. iz2 == 6)
         
         sc% num_calls = sc% num_calls + 1

         if (use_cache .and. associated(cache)) then ! check
            max_z_for_cache = max(size(cache, dim=2), size(cache, dim=3))
            if (iz1 <= max_z_for_cache .and. iz2 <= max_z_for_cache) then
               if (cache(1,iz1,iz2) /= 0) then
                  scor = cache(1,iz1,iz2)
                  dscor_dT = cache(2,iz1,iz2)
                  dscor_dRho = cache(3,iz1,iz2)
                  sc% num_cache_hits = sc% num_cache_hits + 1
                  if (debug) write(*,*) 'graboske_et_al_screening cache hit'
                  return
               end if
            end if
         else
            max_z_for_cache = -1
         end if

         if (debug) write(*,*) 'graboske_et_al_screening eval'
         zbar= sc% zbar
         abar= sc% abar
         z2bar = sc% z2bar
         z1pt58bar = sc% z1pt58bar
         zbar13 = sc% zbar13

         if (zbar <= 0) then
            ierr = -1
            write(*,*) 'bad zbar arg for screening', zbar
            return
         end if

         zg2_screen = 0.316d0*zbar13*zg2
         zg3_screen = (0.737d0/zbar)*zg3

         T = sc% temp
         rho = sc% den

         ztilda = sc% ztilda ! sqrt(z2bar + zbar*theta_e)  ! (Dewitt eqn 4)
                
         H120 = 0 
         ! Graboske 73, eqn 19. screening function H120 = k_b * eta_b * zeta_b * Lambda0^b
         ! with b, k_b, eta_b, and zeta_b given in Table 4.

         Lambda0 = sc% Lambda0 
            ! 0.88d8*sqrt(rho/(abar*T**3)) ! (Graboske eqn 19; mu_I = abar)
         
         d_Lambda0_dRho = 0.5d0*Lambda0/rho
         d_Lambda0_dT = -1.5d0*Lambda0/T
                  
         Lambda12 = z1*z2*ztilda*Lambda0  ! (Dewitt eqn 6)
         
         ! if Lambda12 < 0.1, weak screening
         ! else if Lambda12 < 2, intermediate screening
         ! else if Lambda12 < 5, min of intermediate or strong screening
         ! else strong screening
         
         if (is_bad(Lambda12)) then
            ierr = -1
            return
         end if
         
         ! calculate weak screening
         b = 1
         k_b = 1d0/2d0
         eta_b = ztilda
         zeta_b = 2*z1*z2
         H120 = Lambda12
         dH120_dT = d_Lambda0_dT * (H120 * b / Lambda0)
         dH120_dRho = d_Lambda0_dRho * (H120 * b / Lambda0)

         if (Lambda12 >= Lam_1) then

            H120_weak = H120
            dH120_dT_weak = dH120_dT
            dH120_dRho_weak = dH120_dRho
            
            ! calculate intermediate screening
            b = 0.860d0
            k_b = 0.380d0
            ! 3*b-1 = 1.58, z1pt58bar = <z^(3*b-1)>
            eta_b = z1pt58bar/(sc% ztilda0pt58*sc% zbar0pt28) 
            zeta_b = zg4                         
            H120 = k_b * eta_b * zeta_b * sc% Lambda0b ! intermediate screening value
            dH120_dT = d_Lambda0_dT * (H120 * b / Lambda0)
            dH120_dRho = d_Lambda0_dRho * (H120 * b / Lambda0)
            H120_intermediate = H120
            dH120_dT_intermediate = dH120_dT
            dH120_dRho_intermediate = dH120_dRho

            ! calculate strong screening
            b = 2d0/3d0
            k_b = 0.624d0
            eta_b = zbar13                                                          
            zeta_b = zg1 &
                      + zg2_screen &
                      + zg3_screen / sc% Lambda0_23
            H120_strong = k_b * eta_b * zeta_b * sc% Lambda0b ! strong screening value
            dH120_dT_strong = d_Lambda0_dT * (H120 * b / Lambda0)
            dH120_dRho_strong = d_Lambda0_dRho * (H120 * b / Lambda0)

            if (Lambda12 <= Lam_2) then ! blend intermediate with weak
            
               alfa = (Lambda12 - Lam_1)/(Lam_2 - Lam_1) ! fraction intermediate
               H120 = alfa*H120_intermediate + (1-alfa)*H120_weak
               dH120_dT = alfa*dH120_dT_intermediate + (1-alfa)*dH120_dT_weak
               dH120_dRho = alfa*dH120_dRho_intermediate + (1-alfa)*dH120_dRho_weak
               if (debug) write(*,*) 'blend intermediate with weak'
               
            else if (Lambda12 <= Lam_3) then
            
               H120 = H120_intermediate
               dH120_dT = dH120_dT_intermediate
               dH120_dRho = dH120_dRho_intermediate
               if (debug) write(*,*) 'intermediate'
               
            else if (Lambda12 >= Lam_6) then
            
               H120 = H120_strong
               dH120_dT = dH120_dT_strong
               dH120_dRho = dH120_dRho_strong
               if (debug) write(*,*) 'strong'
               
            else if (H120_intermediate <= H120_strong) then
            
               H120 = H120_intermediate
               dH120_dT = dH120_dT_intermediate
               dH120_dRho = dH120_dRho_intermediate
               if (Lambda12 >= Lam_5) then ! blend with strong
                  alfa = (Lambda12 - Lam_5) / (Lam_6 - Lam_5)
                  H120 = alfa*H120_strong + (1-alfa)*H120
                  dH120_dT = alfa*dH120_dT_strong + (1-alfa)*dH120_dT
                  dH120_dRho = alfa*dH120_dRho_strong + (1-alfa)*dH120_dRho
                  if (debug) write(*,*) 'blend intermediate with strong'
               end if
               
            else
            
               H120 = H120_strong
               dH120_dT = dH120_dT_strong
               dH120_dRho = dH120_dRho_strong
               if (Lambda12 <= Lam_4) then ! blend with intermediate
                  alfa = (Lambda12 - Lam_3) / (Lam_4 - Lam_3)
                  H120 = alfa*H120 + (1-alfa)*H120_intermediate
                  dH120_dT = alfa*dH120_dT + (1-alfa)*dH120_dT_intermediate
                  dH120_dRho = alfa*dH120_dRho + (1-alfa)*dH120_dRho_intermediate
                  if (debug) write(*,*) 'blend intermediate with strong'
               else
                  if (debug) write(*,*) 'strong'
               end if
               
            end if
         
         else
            if (debug) write(*,*) 'weak screening', Lambda12, Lam_1
               
         end if
         
         if (debug) write(*,*) 'H120', H120
                        
         !..machine limit the output
         H120 = max(min(H120, H120_max), 0.0d0) 
         scor = exp(H120) 
         if (H120 == H120_max) then
            dscor_dT = 0.0d0
            dscor_dRho = 0.0d0
         else 
            ! H120 = k_b * eta_b * zeta_b * Lambda0**b
            dscor_dT = scor * dH120_dT
            dscor_dRho = scor * dH120_dRho
         end if

         if (iz1 <= max_z_for_cache .and. iz2 <= max_z_for_cache) then ! store
            cache(1,iz1,iz2) = scor
            cache(2,iz1,iz2) = dscor_dT
            cache(3,iz1,iz2) = dscor_dRho
         end if
      
      if (.not. debug) return
      
111   format(a40,1pe26.16)
      write(*,111) 'ROK', rho
      write(*,111) 'TK', T
      write(*,111) 'L0', Lambda0
      write(*,111) 'z1', z1
      write(*,111) 'z2', z2
      write(*,111) 'ZSCHL', ztilda
      write(*,111) 'L12', Lambda12
      write(*,111) 'H120', H120
      write(*,111) 'SCREEN', scor
      write(*,111) 'zbar', zbar
      write(*,111) 'zbar/abar', zbar/abar
      write(*,111) 'z2bar', z2bar
      write(*,111) 'z2bar/abar', z2bar/abar
      write(*,111) 'z1pt58bar', z1pt58bar
      write(*,111) 'abar', abar
      write(*,111) '1/abar', 1d0/abar
      write(*,111) 'theta_e', theta_e
      stop 'graboske_et_al_screening'

      if (debug) write(*, *) 'scor', scor
      if (debug) write(*, *) 'dscor_dT', dscor_dT
      if (debug) write(*, *) 'dscor_dRho', dscor_dRho
      if (debug) write(*, *)
      if (debug) write(*, *)
      if (debug) stop 'graboske_et_al_screening'
         
      end subroutine graboske_et_al_screening                                                              
      

      end module screen_graboske

