! ***********************************************************************
!
!   Copyright (C) 2016  Bill Paxton
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

      module hydro_reconstruct

      use star_private_def
      use const_def
      use utils_lib

      implicit none

      private
      public :: do_uface_and_Pface, get_G, do1_uface_and_Pface


      contains
      

      subroutine do_uface_and_Pface(s, ierr)
         type (star_info), pointer :: s 
         integer, intent(out) :: ierr
         integer :: k, op_err
         include 'formats'
         ierr = 0
!$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
         do k = 1, s% nz
            op_err = 0
            call do1_uface_and_Pface(s, k, op_err)
            if (op_err /= 0) ierr = op_err
         end do
!$OMP END PARALLEL DO
      end subroutine do_uface_and_Pface


      subroutine get_G(s, k, G, dG_dlnR, dG_dw)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: G, dG_dlnR, dG_dw
         real(dp) :: cgrav, gr_factor, d_gr_factor_dlnR
         cgrav = s% cgrav(k)
         if (s% rotation_flag .and. s% use_gravity_rotation_correction) &
            cgrav = cgrav*s% fp_rot(k)
         G = cgrav
         dG_dlnR = 0d0
         if (s% rotation_flag .and. s% use_gravity_rotation_correction &
               .and. s% w_div_wc_flag) then
            dG_dw = G/s% fp_rot(k)*s% dfp_rot_dw_div_wc(k)
         else
            dG_dw = 0d0
         end if
      end subroutine get_G

      
      subroutine do1_uface_and_Pface(s, k, ierr)
         use eos_def, only: i_gamma1, i_lnfree_e, i_lnPgas
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         integer, intent(out) :: ierr
         logical :: test_partials
         integer :: iL, iR
         ! there are leftover here from various experiments that have gone away
         ! needs a sweep to remove unused locals
         real(dp) :: &
            numerator, denominator, Sl1, Sr1, Sl2, Sr2, Sl, Sr, &
            Ss, PL, PR, uL, uR, rhoL, rhoR
         real(dp) :: &
            G, dG_dlnR, dG_dw, r, A, csL, d_csL_dlndL, d_csL_dlnTL, d_dlnd, d_dlnT, &
            d_PL_dlndL, d_PL_dlnTL, d_PL_dlndR, d_PL_dlnTR, d_dPdm_grav_dlnR, grav, &
            csR, d_csR_dlndR, d_csR_dlnTR, &
            d_PR_dlndR, d_PR_dlndL, d_PR_dlnTR, d_PR_dlnTL, &
            dPdm_grav, &
            d_Sl1_duL, d_Sl1_duR, d_Sl1_dlndL, &
            d_Sl1_dlndR, d_Sl1_dlnTL, d_Sl1_dlnTR, &
            d_Sr1_duL, d_Sr1_duR, d_Sr1_dlndL, &
            d_Sr1_dlndR, d_Sr1_dlnTL, d_Sr1_dlnTR, &
            d_Sl2_duL, d_Sl2_duR, d_Sl2_dlndL, &
            d_Sl2_dlndR, d_Sl2_dlnTL, d_Sl2_dlnTR, &
            d_Sr2_duL, d_Sr2_duR, d_Sr2_dlndL, &
            d_Sr2_dlndR, d_Sr2_dlnTL, d_Sr2_dlnTR, &
            d_Sl_duL, d_Sl_duR, d_Sl_dlndL, &
            d_Sl_dlndR, d_Sl_dlnTL, d_Sl_dlnTR, &
            d_Sr_duL, d_Sr_duR, d_Sr_dlndL, &
            d_Sr_dlndR, d_Sr_dlnTL, d_Sr_dlnTR, &
            d_num_duL, d_num_duR, d_num_dlndL, &
            d_num_dlndR, d_num_dlnTL, d_num_dlnTR, &
            d_den_duL, d_den_duR, d_den_dlndL, &
            d_den_dlndR, d_den_dlnTL, d_den_dlnTR, &
            d_Ss_duL, d_Ss_duR, d_Ss_dlndL, d_Ss_dlnR, &
            d_Ss_dlndR, d_Ss_dlnTL, d_Ss_dlnTR, d_Ss_dv, d_Ss_dw, &
            P_face_L, d_PfaceL_duL, d_PfaceL_dlndL, d_PfaceL_dlnTL, &
            d_PfaceL_duR, d_PfaceL_dlndR, d_PfaceL_dlnTR, &
            P_face_R, d_PfaceR_duR, d_PfaceR_dlndR, d_PfaceR_dlnTR, &
            d_PfaceR_duL, d_PfaceR_dlndL,d_PfaceR_dlnTL, &
            rho_bar, d_rho_bar_dlndL, d_rho_bar_dlndR, &
            cs_bar, d_cs_bar_dlndL, d_cs_bar_dlnTL, &
            d_cs_bar_dlndR, d_cs_bar_dlnTR, &
            rho_cs_bar, d_rhocsbar_dlndL, d_rhocsbar_dlnTL, &
            d_rhocsbar_dlndR, d_rhocsbar_dlnTR, dP_div_rhocs, du, &
            duL_du, duR_du, drhoL_dlndL, drhoR_dlndR, &
            dP, dgamma1, dcsl, dlnd, d_num_dlnR, d_num_dw, &
            d_PL_dlnR, d_PR_dlnR, d_PL_dw, d_PR_dw, Et, &
            d_PfaceL_dlnR, d_PfaceR_dlnR, d_PfaceL_dw, d_PfaceR_dw, &
            kap_face, d_kap_face_dlndL, d_kap_face_dlnTL, &
            d_kap_face_dlndR, d_kap_face_dlnTR, &
            dPdm_rad, d_dPdm_rad_dlnR, delta_m, &
            alfa, beta, T4, prad, f, d_du_dlndR, d_du_dlndL, d_du_dlnR, &
            xh, xhe, Z, dPL_dPR_00, dPR_dPL_m1, duL_duR_00, duR_duL_m1, &
            drhoL_drhoR_00, drhoR_drhoL_m1, dcsL_dcsR_00, dcsR_dcsL_m1, &         
            d_Sl1_duR_00, d_Sl1_dcsR_00, d_Sl2_duL_m1, d_Sl2_dcsL_m1, &
            d_Sl_duR_00, d_Sl_dcsR_00, d_Sl_duL_m1, d_Sl_dcsL_m1, &
            d_Sr1_duL_m1, d_Sr1_dcsL_m1, d_Sr2_duR_00, d_Sr2_dcsR_00, &
            d_Sr_duL_m1, d_Sr_dcsL_m1, d_Sr_duR_00, d_Sr_dcsR_00, &
            d_num_dPR_00, d_num_dPL_m1, d_num_duR_00, d_num_duL_m1, &
            d_num_drhoR_00, d_num_drhoL_m1, d_num_dcsR_00, d_num_dcsL_m1, &
            d_den_duR_00, d_den_duL_m1, d_den_drhoR_00, d_den_drhoL_m1, &
            d_den_dcsR_00, d_den_dcsL_m1, &
            d_Ss_dPR_00, d_Ss_dPL_m1, d_Ss_duR_00, d_Ss_duL_m1, &
            d_Ss_drhoR_00, d_Ss_drhoL_m1, d_Ss_dcsR_00, d_Ss_dcsL_m1, &
            d_ufaceL_dPR_00,  d_ufaceR_dPR_00, &
            d_ufaceL_dPL_m1,  d_ufaceR_dPL_m1, &
            d_ufaceL_duR_00,  d_ufaceR_duR_00, &
            d_ufaceL_duL_m1,  d_ufaceR_duL_m1, &
            d_ufaceL_drhoR_00,  d_ufaceR_drhoR_00, &
            d_ufaceL_drhoL_m1,  d_ufaceR_drhoL_m1, &
            d_ufaceL_dcsR_00,  d_ufaceR_dcsR_00, &
            d_ufaceL_dcsL_m1,  d_ufaceR_dcsL_m1, &
            d_PfaceL_dPR_00,  d_PfaceR_dPR_00, &
            d_PfaceL_dPL_m1,  d_PfaceR_dPL_m1, &
            d_PfaceL_duR_00,  d_PfaceR_duR_00, &
            d_PfaceL_duL_m1,  d_PfaceR_duL_m1, &
            d_PfaceL_drhoR_00,  d_PfaceR_drhoR_00, &
            d_PfaceL_drhoL_m1,  d_PfaceR_drhoL_m1, &
            d_PfaceL_dcsR_00,  d_PfaceR_dcsR_00, &
            d_PfaceL_dcsL_m1,  d_PfaceR_dcsL_m1
            
         include 'formats'
         
         ierr = 0
         test_partials = .false.
         !test_partials = (k == s% solver_test_partials_k)
         
         s% RTI_du_diffusion_kick(k) = 0d0
         
         s% d_uface_du00(k) = 0
         s% d_uface_dum1(k) = 0
         s% d_uface_dlnd00(k) = 0
         s% d_uface_dlndm1(k) = 0
         s% d_uface_dlnT00(k) = 0
         s% d_uface_dlnTm1(k) = 0
         s% d_uface_dlnR(k) = 0
         s% d_uface_dw(k) = 0

         s% d_Pface_du00(k) = 0
         s% d_Pface_dum1(k) = 0
         s% d_Pface_dlnd00(k) = 0
         s% d_Pface_dlndm1(k) = 0
         s% d_Pface_dlnT00(k) = 0
         s% d_Pface_dlnTm1(k) = 0
         s% d_Pface_dlnR(k) = 0
         s% d_Pface_dL(k) = 0
         s% d_Pface_dw(k) = 0
                            
         if (k == 1) then         
            s% u_face(k) = s% u(k)
            s% d_uface_du00(k) = 1d0
            s% P_face(k) = s% P(k)
            s% d_Pface_dlnd00(k) = s% P(k)*s% chiRho_for_partials(k)
            s% d_Pface_dlnT00(k) = s% P(k)*s% chiT_for_partials(k)
            s% d_Pface_dlndm1(k) = 0d0
            s% d_Pface_dlnTm1(k) = 0d0
            return            
         end if

         d_PL_dlndL = 0
         d_PL_dlnTL = 0
         d_PR_dlndR = 0
         d_PR_dlnTR = 0           
         d_PL_dlndR = 0
         d_PL_dlnTR = 0
         d_PL_dlnR = 0
         d_PL_dw = 0
         d_PR_dlndL = 0
         d_PR_dlnTL = 0
         d_PR_dlnR = 0
         d_PR_dw = 0
         duL_du = 0         
         drhoL_dlndL = 0         
         d_csL_dlnTL = 0            
         d_csL_dlndL = 0
         duR_du = 0         
         drhoR_dlndR = 0         
         d_csR_dlnTR = 0          
         d_csR_dlndR = 0
         dPL_dPR_00 = 0
         dPR_dPL_m1 = 0
         duL_duR_00 = 0
         duR_duL_m1 = 0
         drhoL_drhoR_00 = 0
         drhoR_drhoL_m1 = 0
         dcsL_dcsR_00 = 0
         dcsR_dcsL_m1 = 0
      
         r = s% r(k)
         A = 4d0*pi*r*r
      
         iL = k
         iR = k-1
         
         
         PL = s% P(iL)
         PR = s% P(iR)
         d_PL_dlndL = s% chiRho_for_partials(iL)*PL
         d_PL_dlnTL = s% chiT_for_partials(iL)*PL
         d_PR_dlndR = s% chiRho_for_partials(iR)*PR
         d_PR_dlnTR = s% chiT_for_partials(iR)*PR

         uL = s% u(iL)
         duL_du = 1d0
      
         rhoL = s% rho(iL)
         drhoL_dlndL = rhoL
      
         csL = s% csound(iL) ! sqrt(gamma1*P/rho) 
         d_csL_dlnTL = s% P(iL)*(s% d_eos_dlnT(i_gamma1,iL) + &
            s% gamma1(iL)*s% chiT_for_partials(iL))/(2d0*rhoL*csL)               
         d_csL_dlndL = s% P(iL)*(s% d_eos_dlnd(i_gamma1,iL) + &
            s% gamma1(iL)*(s% chiRho_for_partials(iL) - 1d0))/(2d0*rhoL*csL)
      
         uR = s% u(iR)
         duR_du = 1d0
      
         rhoR = s% rho(iR)
         drhoR_dlndR = rhoR
      
         csR = s% csound(iR) ! sqrt(gamma1*P/rho)
         d_csR_dlnTR = s% P(iR)*(s% d_eos_dlnT(i_gamma1,iR) + &
            s% gamma1(iR)*s% chiT_for_partials(iR))/(2d0*rhoR*csR)               
         d_csR_dlndR = s% P(iR)*(s% d_eos_dlnd(i_gamma1,iR) + &
            s% gamma1(iR)*(s% chiRho_for_partials(iR) - 1d0))/(2d0*rhoR*csR)
         
         ! change PR and PL for gravity
         call get_G(s, k, G, dG_dlnR, dG_dw)
         dPdm_grav = -G*s% m_grav(k)/(r*r*A)  ! cm^-1 s^-2
         d_dPdm_grav_dlnR = -4d0*dPdm_grav - dG_dlnR*s% m_grav(k)/(r*r*A)            
         delta_m = 0.5d0*s% dm(iL) ! positive delta_m from left center to edge
         PL = PL + delta_m*dPdm_grav
         d_PL_dlnR = delta_m*d_dPdm_grav_dlnR            
         d_PL_dw = delta_m*dPdm_grav/G*dG_dw     
         delta_m = -0.5d0*s% dm(iR) ! negative delta_m from right center to edge
         PR = PR + delta_m*dPdm_grav
         d_PR_dlnR = delta_m*d_dPdm_grav_dlnR
         d_PR_dw = delta_m*dPdm_grav/G*dG_dw
            
         ! acoustic wavespeeds (eqn 2.38)
         Sl1 = uL - csL
         d_Sl1_duL = 1d0
         d_Sl1_duR = 0d0
         d_Sl1_dlndL = -d_csL_dlndL
         d_Sl1_dlndR = 0d0
         d_Sl1_dlnTL = -d_csL_dlnTL
         d_Sl1_dlnTR = 0d0
         
         Sl2 = uR - csR
         d_Sl2_duL = 0d0
         d_Sl2_duR = 1d0
         d_Sl2_dlndL = 0d0
         d_Sl2_dlndR = -d_csR_dlndR
         d_Sl2_dlnTL = 0d0
         d_Sl2_dlnTR = -d_csR_dlnTR

         ! take Sl = min(Sl1, Sl2)
         if (Sl1 < Sl2) then
            Sl = Sl1
            d_Sl_duL = d_Sl1_duL
            d_Sl_duR = d_Sl1_duR
            d_Sl_dlndL = d_Sl1_dlndL
            d_Sl_dlndR = d_Sl1_dlndR
            d_Sl_dlnTL = d_Sl1_dlnTL
            d_Sl_dlnTR = d_Sl1_dlnTR
         else
            Sl = Sl2
            d_Sl_duL = d_Sl2_duL
            d_Sl_duR = d_Sl2_duR
            d_Sl_dlndL = d_Sl2_dlndL
            d_Sl_dlndR = d_Sl2_dlndR
            d_Sl_dlnTL = d_Sl2_dlnTL
            d_Sl_dlnTR = d_Sl2_dlnTR
         end if

         Sr1 = uR + csR
         d_Sr1_duR = 1d0
         d_Sr1_duL = 0d0
         d_Sr1_dlndL = 0d0
         d_Sr1_dlndR = d_csR_dlndR
         d_Sr1_dlnTL = 0d0
         d_Sr1_dlnTR = d_csR_dlnTR
         
         Sr2 = uL + csL
         d_Sr2_duL = 1d0
         d_Sr2_duR = 0d0
         d_Sr2_dlndL = d_csL_dlndL
         d_Sr2_dlndR = 0d0
         d_Sr2_dlnTL = d_csL_dlnTL
         d_Sr2_dlnTR = 0d0
         
         ! take Sr = max(Sr1, Sr2)
         if (Sr1 > Sr2) then
            Sr = Sr1
            d_Sr_duL = d_Sr1_duL
            d_Sr_duR = d_Sr1_duR
            d_Sr_dlndL = d_Sr1_dlndL
            d_Sr_dlndR = d_Sr1_dlndR
            d_Sr_dlnTL = d_Sr1_dlnTL
            d_Sr_dlnTR = d_Sr1_dlnTR
         else
            Sr = Sr2
            d_Sr_duL = d_Sr2_duL
            d_Sr_duR = d_Sr2_duR
            d_Sr_dlndL = d_Sr2_dlndL
            d_Sr_dlndR = d_Sr2_dlndR
            d_Sr_dlnTL = d_Sr2_dlnTL
            d_Sr_dlnTR = d_Sr2_dlnTR
         end if
         
         ! contact velocity (eqn 2.20)
         numerator = uR*rhoR*(Sr - uR) + uL*rhoL*(uL - Sl) + (PL - PR)         
         d_num_dlnR = d_PL_dlnR - d_PR_dlnR         
         d_num_duR = rhoR*(sR + uR*(d_Sr_duR - 2d0)) - rhoL*uL*d_Sl_duR         
         d_num_duL = uR*rhoR*d_Sr_duL + uL*rhoL*(1d0 - d_Sl_duL) + rhoL*(uL - Sl)                  
         d_num_dlndR = &
            uR*rhoR*(sR - uR + d_Sr_dlndR) - uL*rhoL*d_Sl_dlndR + (d_PL_dlndR - d_PR_dlndR)
         d_num_dlnTR = &
            uR*rhoR*d_Sr_dlnTR - uL*rhoL*d_Sl_dlnTR + (d_PL_dlnTR - d_PR_dlnTR)          
         d_num_dlndL = &
            uR*rhoR*d_Sr_dlndL + uL*rhoL*(uL - Sl - d_Sl_dlndL) + (d_PL_dlndL - d_PR_dlndL)
         d_num_dlnTL = &
            uR*rhoR*d_Sr_dlnTL - uL*rhoL*d_Sl_dlnTL + (d_PL_dlnTL - d_PR_dlnTL)       
         d_num_dw = d_PL_dw - d_PR_dw
         
         denominator = rhoR*(Sr - uR) + rhoL*(uL - Sl)         
         d_den_duR = rhoR*(d_Sr_duR - 1d0) - rhoL*d_Sl_duR
         d_den_duL = rhoR*d_Sr_duL + rhoL*(1d0 - d_Sl_duL)         
         d_den_dlndR = rhoR*(Sr - uR + d_Sr_dlndR) - rhoL*d_Sl_dlndR
         d_den_dlndL = rhoR*d_Sr_dlndL + rhoL*(uL - Sl - d_Sl_dlndL)         
         d_den_dlnTR = rhoR*d_Sr_dlnTR - rhoL*d_Sl_dlnTR
         d_den_dlnTL = rhoR*d_Sr_dlnTL - rhoL*d_Sl_dlnTL

         if (denominator == 0d0 .or. is_bad(denominator)) then
            ierr = -1
            if (s% report_ierr) then
               write(*,2) 'u_face denominator bad', k, denominator
               write(*,2) 'rhoR', k, rhoR
               write(*,2) 'rhoL', k, rhoL
               write(*,2) 'Sr - uR', k, Sr - uR
               write(*,2) 'uL - Sl', k, uL - Sl
               write(*,2) 'Sr', k, Sr
               write(*,2) 'Sl', k, Sl
               write(*,2) 'uR', k, uR
               write(*,2) 'uL', k, uL
            end if
            return
         end if
         
         Ss = numerator/denominator
         d_Ss_dlnR = d_num_dlnR/denominator
         d_Ss_duR = d_num_duR/denominator - Ss*d_den_duR/denominator
         d_Ss_duL = d_num_duL/denominator - Ss*d_den_duL/denominator
         d_Ss_dlndR = d_num_dlndR/denominator - Ss*d_den_dlndR/denominator
         d_Ss_dlndL = d_num_dlndL/denominator - Ss*d_den_dlndL/denominator
         d_Ss_dlnTL = d_num_dlnTL/denominator - Ss*d_den_dlnTL/denominator
         d_Ss_dlnTR = d_num_dlnTR/denominator - Ss*d_den_dlnTR/denominator
         d_Ss_dv = 0
         d_Ss_dw = d_num_dw/denominator
     
         s% u_face(k) = Ss
         if (is_bad(s% u_face(k))) then
            write(*,2) 's% u_face(k)', k, s% u_face(k)
            stop 
         end if

         s% d_uface_dlnR(k) = d_Ss_dlnR
         s% d_uface_du00(k) = d_Ss_duL
         s% d_uface_dlnd00(k) = d_Ss_dlndL
         s% d_uface_dlnT00(k) = d_Ss_dlnTL
         s% d_uface_dum1(k) = d_Ss_duR
         s% d_uface_dlndm1(k) = d_Ss_dlndR
         s% d_uface_dlnTm1(k) = d_Ss_dlnTR
         s% d_uface_dw(k) = d_Ss_dw
         if (k==1) then
            write(*,*) "check k=1", s% d_uface_dw(k), d_Ss_dw
         end if

         ! contact pressure (eqn 2.19)
         P_face_L = rhoL*(uL-Sl)*(uL-Ss) + PL
         d_PfaceL_dlnR = d_PL_dlnR - d_Ss_dlnR*rhoL*(uL-Sl)
         d_PfaceL_duL = &
            rhoL*((1d0 - d_Sl_duL)*(uL-Ss) + (uL-Sl)*(1d0 - d_Ss_duL))
         d_PfaceL_dlndL = d_PL_dlndL + &
            rhoL*((uL-Sl)*(uL-Ss) - d_Sl_dlndL*(uL-Ss) - (uL-Sl)*d_Ss_dlndL)
         d_PfaceL_dlnTL = d_PL_dlnTL - &
            rhoL*(d_Sl_dlnTL*(uL-Ss) + (uL-Sl)*d_Ss_dlnTL)
         d_PfaceL_dlndR = d_PL_dlndR - & 
            rhoL*(d_Sl_dlndR*(uL-Ss) + (uL-Sl)*d_Ss_dlndR)
         d_PfaceL_dlnTR = d_PL_dlnTR - & 
            rhoL*(d_Sl_dlnTR*(uL-Ss) + (uL-Sl)*d_Ss_dlnTR)
         d_PfaceL_duR = &
            -rhoL*(d_Sl_duR*(uL-Ss) + (uL-Sl)*d_Ss_duR)
         d_PfaceL_dw = d_PL_dw - d_Ss_dw*rhoL*(uL-Sl)
         
         P_face_R = rhoR*(uR-Sr)*(uR-Ss) + PR
         d_PfaceR_dlnR = d_PR_dlnR - d_Ss_dlnR*rhoR*(uR-Sr)
         d_PfaceR_duR = &
            rhoR*((1d0-d_Sr_duR)*(uR-Ss) + (uR-Sr)*(1d0-d_Ss_duR))
         d_PfaceR_dlndR = d_PR_dlndR + &
            rhoR*((uR-Sr)*(uR-Ss) - d_Sr_dlndR*(uR-Ss) - (uR-Sr)*d_Ss_dlndR)
         d_PfaceR_dlnTR = d_PR_dlnTR - &
            rhoR*(d_Sr_dlnTR*(uR-Ss) + (uR-Sr)*d_Ss_dlnTR)               
         d_PfaceR_dlndL = d_PR_dlndL - &
            rhoR*(d_Sr_dlndL*(uR-Ss) + (uR-Sr)*d_Ss_dlndL)
         d_PfaceR_dlnTL = d_PR_dlnTL - &
            rhoR*(d_Sr_dlnTL*(uR-Ss) + (uR-Sr)*d_Ss_dlnTL)
         d_PfaceR_duL = &
            -rhoR*(d_Sr_duL*(uR-Ss) + (uR-Sr)*d_Ss_duL)
         d_PfaceR_dw = d_PR_dw - d_Ss_dw*rhoR*(uR-Sr)

         ! P is continuous across the face, so P_face_L == P_face_R if AOK
         ! skip this for low P as happens after end of plateau phase of ccsn IIp
         if (P_face_L > 1d3 .and. P_face_R > 1d3 .and. &
               rhoL > 1d-11 .and. rhoR > 1d-11 .and. &
               1d4*abs(P_face_L - P_face_R) > &
                  (1d0 + max(abs(P_face_L),abs(P_face_R),abs(PL),abs(PR)))) then
            write(*,2) 'bad reconstruction for P_face', k
            ierr = -1
            if (.not. s% report_ierr) return
            stop 'do1_uface_and_Pface'
         end if
         
         s% P_face(k) = 0.5d0*(P_face_L + P_face_R)
         s% d_Pface_du00(k) = 0.5d0*(d_PfaceL_duL + d_PfaceR_duL)
         s% d_Pface_dum1(k) = 0.5d0*(d_PfaceL_duR + d_PfaceR_duR)
         s% d_Pface_dlnd00(k) = 0.5d0*(d_PfaceL_dlndL + d_PfaceR_dlndL)
         s% d_Pface_dlnT00(k) = 0.5d0*(d_PfaceL_dlnTL + d_PfaceR_dlnTL)
         s% d_Pface_dlndm1(k) = 0.5d0*(d_PfaceL_dlndR + d_PfaceR_dlndR)
         s% d_Pface_dlnTm1(k) = 0.5d0*(d_PfaceL_dlnTR + d_PfaceR_dlnTR)         
         s% d_Pface_dlnR(k) = 0.5d0*(d_PfaceL_dlnR + d_PfaceR_dlnR)
         s% d_Pface_dL(k) = 0d0
         s% d_Pface_dw(k) = 0.5d0*(d_PfaceL_dw + d_PfaceR_dw)

         if(k < s% nz .and. s% RTI_flag) then
             if (s% eta_RTI(k) > 0d0 .and. &
                   s% dlnddt_RTI_diffusion_factor > 0d0 .and. s% dt > 0d0) then
                f = A*s% dlnddt_RTI_diffusion_factor*s% eta_RTI(k)/s% dm_bar(k)
                du = f*(rhoL - rhoR) ! bump uface in direction of lower density
                d_du_dlndL = f*rhoL
                d_du_dlndR = -f*rhoR
                d_du_dlnR = 2*du            
                s% RTI_du_diffusion_kick(k) = du            
                s% u_face(k) = s% u_face(k) + du
                s% d_uface_dlnd00(k) = s% d_uface_dlnd00(k) + d_du_dlndL
                s% d_uface_dlndm1(k) = s% d_uface_dlndm1(k) + d_du_dlndR
                s% d_uface_dlnR(k) = s% d_uface_dlnR(k) + d_du_dlnR
             end if
         end if

            !d_PfaceR_dlnTR =  - &
            !   rhoR*(d_Sr_dlnTR*(uR-Ss) + (uR-Sr)*d_Ss_dlnTR)               
         if (test_partials) then
            s% solver_test_partials_val = PL
            s% solver_test_partials_var = s% i_w_div_wc           
            s% solver_test_partials_dval_dx = d_PL_dw
            !stop
         end if
     
      end subroutine do1_uface_and_Pface


      end module hydro_reconstruct

