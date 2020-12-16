! ***********************************************************************
!
!   Copyright (C) 2015-2019  Bill Paxton & The MESA Team
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

      module hydro_riemann
      
      use star_private_def
      use const_def
      use star_utils, only: em1, e00, ep1
      use utils_lib
      
      implicit none

      ! Cheng, J, Shu, C-W, and Zeng, Q.,
      ! "A Conservative Lagrangian Scheme for Solving
      !  Compressible Fluid Flows with Multiple Internal Energy Equations",
      ! Commun. Comput. Phys., 12, pp 1307-1328, 2012.
      
      ! Cheng, J. and Shu, C-W, 
      ! "Positivity-preserving Lagrangian scheme for multi-material
      !  compressible flow", J. Comp. Phys., 257 (2014), 143-168.

      ! Kappeli, R. and Mishra, S.,
      ! "Well-balanced schemes for the Euler equations with gravitation",
      ! J. Comp. Phys., 259 (2014), 199-219.


      private
      public :: &
         do_surf_Riemann_dudt_eqn, do1_Riemann_momentum_eqn, &
         do1_Riemann_energy_eqn, do1_Riemann_dlnRdt_eqn

      contains
      
      ! ergs/s from cell(k) to cell(k-1)
      subroutine eval1_energy_flux(s, k, &
            eflux, d_eflux_dL, d_eflux_dlnR, &
            d_eflux_du00, d_eflux_dum1, &
            d_eflux_dlnd00, d_eflux_dlndm1, d_eflux_dlnT00, d_eflux_dlnTm1, d_eflux_dw, &
            ierr)
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         real(qp), intent(out) :: eflux
         real(dp), intent(out) :: &
            d_eflux_dL, d_eflux_dlnR, &
            d_eflux_du00, d_eflux_dum1, &
            d_eflux_dlnd00, d_eflux_dlndm1, d_eflux_dlnT00, d_eflux_dlnTm1, d_eflux_dw
         integer, intent(out) :: ierr
         real(qp) :: r, A, dA_dlnR, P_face, u_face, L, &
            eta, dr, vol, gam, ieL, ieR, rhoR, rhoL, vR, vL
         integer :: nz
         logical :: test_partials
         include 'formats'
         ierr = 0

         !test_partials = (k == s% solver_test_partials_k) 
         test_partials = .false.
         
         d_eflux_dL = 0 
         d_eflux_dlnR = 0 
         d_eflux_du00 = 0 
         d_eflux_dum1 = 0 
         d_eflux_dlnd00 = 0 
         d_eflux_dlndm1 = 0 
         d_eflux_dlnT00 = 0 
         d_eflux_dlnTm1 = 0 
         d_eflux_dw = 0 

         if (s% constant_L) then
            L = 0d0
         else if (k > s% nz) then
            L = s% L_center
         else
            L = s% L(k)
            d_eflux_dL = 1d0
         end if
         
         nz = s% nz

         if (k > nz) then
            eflux = L
            return    
         end if
         
         r = s% r(k)
         A = 4d0*pi*r*r
         dA_dlnR = 2d0*A
         u_face = s% u_face(k)
         P_face = s% P_face(k)

         eflux = A*P_face*u_face + L
         
         if (is_bad(dble(eflux))) then
            ierr = -1
            if (.not. s% report_ierr) return
            write(*,2) 'A', k, A
            write(*,2) 'P_face', k, P_face
            write(*,2) 'u_face', k, u_face
            write(*,2) 'L', k, L
            write(*,2) 'eflux', k, eflux
            stop 'eval1_energy_flux'
         end if
         
         d_eflux_dlnR = &
            dA_dlnR*P_face*u_face + &
            A*s% d_Pface_dlnR(k)*u_face + &
            A*P_face*s% d_uface_dlnR(k)
         d_eflux_dL = d_eflux_dL + A*s% d_Pface_dL(k)*u_face
         d_eflux_dlnd00 = &
            A*(s% d_Pface_dlnd00(k)*u_face + P_face*s% d_uface_dlnd00(k))
         d_eflux_dlndm1 = &
           A*(s% d_Pface_dlndm1(k)*u_face + P_face*s% d_uface_dlndm1(k))
         d_eflux_du00 = &
            A*(s% d_Pface_du00(k)*u_face + P_face*s% d_uface_du00(k))
         d_eflux_dum1 = &
            A*(s% d_Pface_dum1(k)*u_face + P_face*s% d_uface_dum1(k))
         d_eflux_dlnT00 = &
            A*(s% d_Pface_dlnT00(k)*u_face + P_face*s% d_uface_dlnT00(k))
         d_eflux_dlnTm1 = &
            A*(s% d_Pface_dlnTm1(k)*u_face + P_face*s% d_uface_dlnTm1(k))
         d_eflux_dw = &
            A*s% d_Pface_dw(k)*u_face + &
            A*P_face*s% d_uface_dw(k)
         
         if (test_partials) then
            s% solver_test_partials_val = u_face
            s% solver_test_partials_var = s% i_lnT
            s% solver_test_partials_dval_dx = s% d_uface_dlnT00(k)
         end if
         
      end subroutine eval1_energy_flux
      
      
      subroutine eval_surf_energy_flux(s, Psurf, &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT, &
            eflux, d_eflux_dL, d_eflux_dlnR, &
            d_eflux_du00, d_eflux_dum1, &
            d_eflux_dlnd00, d_eflux_dlndm1, d_eflux_dlnT00, d_eflux_dlnTm1, d_eflux_dw, &
            ierr)
         type (star_info), pointer :: s 
         real(dp), intent(in) :: &
            Psurf, dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT
         real(qp), intent(out) :: eflux
         real(dp), intent(out) :: &
            d_eflux_dL, d_eflux_dlnR, &
            d_eflux_du00, d_eflux_dum1, &
            d_eflux_dlnd00, d_eflux_dlndm1, d_eflux_dlnT00, d_eflux_dlnTm1, d_eflux_dw
         integer, intent(out) :: ierr
         real(qp) :: r, A, dA_dlnR, P_face, u_face, L, &
            eta, dr, vol, gam, ieL, ieR, rhoR, rhoL, vR, vL
         integer :: nz
         logical :: test_partials
         include 'formats'
         ierr = 0

         !test_partials = (k == s% solver_test_partials_k) 
         test_partials = .false.
         
         d_eflux_dL = 0 
         d_eflux_dlnR = 0 
         d_eflux_du00 = 0 
         d_eflux_dum1 = 0 
         d_eflux_dlnd00 = 0 
         d_eflux_dlndm1 = 0 
         d_eflux_dlnT00 = 0 
         d_eflux_dlnTm1 = 0 
         d_eflux_dw = 0 

         if (s% constant_L) then
            L = 0d0
         else
            L = s% L(1)
            d_eflux_dL = 1d0
         end if
         
         nz = s% nz
         
         r = s% r(1)
         A = 4d0*pi*r*r
         dA_dlnR = 2d0*A
         u_face = s% u_face(1)

         eflux = A*Psurf*u_face + L
         
         if (is_bad(dble(eflux))) then
            ierr = -1
            if (.not. s% report_ierr) return
            write(*,2) 'Psurf', Psurf
            write(*,2) 'u_face', u_face
            write(*,2) 'L', L
            write(*,2) 'eflux', eflux
            stop 'eval_surf_energy_flux'
         end if
         
         d_eflux_dlnR = &
            dA_dlnR*Psurf*u_face + &
            A*Psurf*dlnPsurf_dlnR*u_face + &
            A*Psurf*s% d_uface_dlnR(1)
         d_eflux_dL = d_eflux_dL + A*Psurf*dlnPsurf_dL*u_face
         d_eflux_dlnd00 = &
            A*(Psurf*dlnPsurf_dlnd*u_face + Psurf*s% d_uface_dlnd00(1))
         d_eflux_dlndm1 = A*Psurf*s% d_uface_dlndm1(1)
         d_eflux_du00 = A*Psurf*s% d_uface_du00(1)
         d_eflux_dum1 = A*Psurf*s% d_uface_dum1(1)
         d_eflux_dlnT00 = &
            A*(Psurf*dlnPsurf_dlnT*u_face + Psurf*s% d_uface_dlnT00(1))
         d_eflux_dlnTm1 = A*Psurf*s% d_uface_dlnTm1(1)
         d_eflux_dw = &
            A*Psurf*s% d_uface_dw(1)
         
         if (test_partials) then
            s% solver_test_partials_val = u_face
            s% solver_test_partials_var = s% i_lnT
            s% solver_test_partials_dval_dx = s% d_uface_dlnT00(1)
         end if
         
      end subroutine eval_surf_energy_flux
      
      ! (g cm/s)/s from cell(k) to cell(k-1)
      subroutine eval1_momentum_flux(s, k, &
            momflux, d_momflux_dlnR, d_momflux_dL, &
            d_momflux_du00, d_momflux_dum1, &
            d_momflux_dlnd00, d_momflux_dlndm1, &
            d_momflux_dlnT00, d_momflux_dlnTm1, &
            d_momflux_dw, ierr)
         type (star_info), pointer :: s 
         integer, intent(in) :: k
         real(qp), intent(out) :: momflux
         real(dp), intent(out) :: &
            d_momflux_dlnR, d_momflux_dL, &
            d_momflux_du00, d_momflux_dum1, &
            d_momflux_dlnd00, d_momflux_dlndm1, &
            d_momflux_dlnT00, d_momflux_dlnTm1, &
            d_momflux_dw
         integer, intent(out) :: ierr
         real(qp) :: r, A, Pface
         integer :: nz
         logical :: test_partials
         include 'formats'
         
         !test_partials = (k-1 == s% solver_test_partials_k)
         test_partials = .false.
         ierr = 0
         d_momflux_dL = 0 
         d_momflux_dlnR = 0 
         d_momflux_du00 = 0 
         d_momflux_dum1 = 0 
         d_momflux_dlnd00 = 0 
         d_momflux_dlndm1 = 0 
         d_momflux_dlnT00 = 0 
         d_momflux_dlnTm1 = 0 
         d_momflux_dw = 0 
         nz = s% nz
         if (k > nz) then
            momflux = 0
            return    
         end if
         r = s% r(k)
         A = 4d0*pi*r*r
         Pface = s% P_face(k)
         momflux = A*Pface
         
         d_momflux_dlnR = 2d0*momflux + A*s% d_Pface_dlnR(k)
         d_momflux_dL = s% d_Pface_dL(k)*A
         d_momflux_du00 = s% d_Pface_du00(k)*A
         d_momflux_dum1 = s% d_Pface_dum1(k)*A
         d_momflux_dlnd00 = s% d_Pface_dlnd00(k)*A
         d_momflux_dlndm1 = s% d_Pface_dlndm1(k)*A
         d_momflux_dlnT00 = s% d_Pface_dlnT00(k)*A
         d_momflux_dlnTm1 = s% d_Pface_dlnTm1(k)*A
         d_momflux_dw = A*s% d_Pface_dw(k)
              
         if (test_partials) then
            s% solver_test_partials_val = momflux
            s% solver_test_partials_var = s% i_lnT
            s% solver_test_partials_dval_dx = d_momflux_dlnTm1
         end if
         
      end subroutine eval1_momentum_flux
      
      
      subroutine eval_surf_momentum_flux(s, P_surf, &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT, &
            momflux, d_momflux_dlnR, d_momflux_dL, &
            d_momflux_du00, d_momflux_dum1, &
            d_momflux_dlnd00, d_momflux_dlndm1, &
            d_momflux_dlnT00, d_momflux_dlnTm1, &
            d_momflux_dw, ierr)
         type (star_info), pointer :: s 
         real(dp), intent(in) :: P_surf, &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT
         real(qp), intent(out) :: momflux
         real(dp), intent(out) :: &
            d_momflux_dlnR, d_momflux_dL, &
            d_momflux_du00, d_momflux_dum1, &
            d_momflux_dlnd00, d_momflux_dlndm1, &
            d_momflux_dlnT00, d_momflux_dlnTm1, &
            d_momflux_dw
         integer, intent(out) :: ierr
         integer :: k
         real(qp) :: r, A
         include 'formats'
         ierr = 0
         k = 1
         r = s% r(k)
         A = 4d0*pi*r*r
         momflux = A*P_surf
         d_momflux_dlnR = momflux*(2 + dlnPsurf_dlnR)
         d_momflux_dL = momflux*dlnPsurf_dL
         d_momflux_dlnd00 = momflux*dlnPsurf_dlnd
         d_momflux_dlnT00 = momflux*dlnPsurf_dlnT
         d_momflux_du00 = 0
         d_momflux_dum1 = 0
         d_momflux_dlndm1 = 0
         d_momflux_dlnTm1 = 0
         d_momflux_dw = 0
      end subroutine eval_surf_momentum_flux


      subroutine do_surf_Riemann_dudt_eqn( &
            s, P_surf, &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT, &
            xscale, equ, skip_partials, nvar, ierr)
         type (star_info), pointer :: s         
         real(dp), intent(in) :: P_surf, &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: skip_partials
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr
         include 'formats'
         call do1_dudt_eqn( &
            s, 1, P_surf, &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT, &
            xscale, equ, skip_partials, nvar, ierr)
      end subroutine do_surf_Riemann_dudt_eqn
      

      subroutine do1_Riemann_momentum_eqn( &
            s, k, P_surf, xscale, equ, skip_partials, nvar, ierr)
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         real(dp), intent(in) :: P_surf ! only used if k==1
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: skip_partials
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr
         call do1_dudt_eqn( &
            s, k, P_surf, 0d0, 0d0, 0d0, 0d0, &
            xscale, equ, skip_partials, nvar, ierr)
      end subroutine do1_Riemann_momentum_eqn
         

      subroutine do1_dudt_eqn( &
            s, k, P_surf, &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT, &
            xscale, equ, skip_partials, nvar, ierr)
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         real(dp), intent(in) :: P_surf ! only used if k==1
         real(dp), intent(in) :: &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: skip_partials
         integer, intent(in) :: nvar
         integer, intent(out) :: ierr
         
         integer :: i_du_dt, i_u, nz
         real(dp) :: &
            residual, &
            e00_u, em1_u, ep1_u, &
            e00_lum, ep1_lum, &
            e00_lnR, ep1_lnR, &
            e00_lnd, e00_lnT, e00_Et, &
            em1_lnd, em1_lnT, em1_Et, &
            ep1_lnd, ep1_lnT, ep1_Et, &  
            e00_PR, ep1_PR, &           
            e00_PL, em1_PL, &
            e00_uR, ep1_uR, &           
            e00_uL, em1_uL, &
            e00_rhoR, ep1_rhoR, &          
            e00_rhoL, em1_rhoL, &
            e00_csR, ep1_csR, &         
            e00_csL, em1_csL, &
            e00_w, ep1_w
         include 'formats'
         
         call do1_dudt_residual_and_partials( &
            s, k, P_surf, &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT, &
            skip_partials, &
            residual, &
            e00_u, em1_u, ep1_u, &
            e00_lum, ep1_lum, &
            e00_lnR, ep1_lnR, &
            e00_lnd, e00_lnT, e00_Et, &
            em1_lnd, em1_lnT, em1_Et, &
            ep1_lnd, ep1_lnT, ep1_Et, &
            e00_PR, ep1_PR, &           
            e00_PL, em1_PL, &
            e00_uR, ep1_uR, &           
            e00_uL, em1_uL, &
            e00_rhoR, ep1_rhoR, &          
            e00_rhoL, em1_rhoL, &
            e00_csR, ep1_csR, &         
            e00_csL, em1_csL, &            
            e00_w, ep1_w, &
            ierr)
         if (ierr /= 0) return
         
         call store1_dudt_residual_and_partials( &
            s, k, xscale, equ, skip_partials, nvar, &
            residual, &
            e00_u, em1_u, ep1_u, &
            e00_lum, ep1_lum, &
            e00_lnR, ep1_lnR, &
            e00_lnd, e00_lnT, e00_Et, &
            em1_lnd, em1_lnT, em1_Et, &
            ep1_lnd, ep1_lnT, ep1_Et, &
            e00_PR, ep1_PR, &           
            e00_PL, em1_PL, &
            e00_uR, ep1_uR, &           
            e00_uL, em1_uL, &
            e00_rhoR, ep1_rhoR, &          
            e00_rhoL, em1_rhoL, &
            e00_csR, ep1_csR, &         
            e00_csL, em1_csL, &            
            e00_w, ep1_w, &
            ierr)

      end subroutine do1_dudt_eqn
      
      
      subroutine do1_dudt_residual_and_partials( &
            s, k, P_surf, &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT, &
            skip_partials, &
            residual, &
            e00_u, em1_u, ep1_u, &
            e00_lum, ep1_lum, &
            e00_lnR, ep1_lnR, &
            e00_lnd, e00_lnT, e00_Et, &
            em1_lnd, em1_lnT, em1_Et, &
            ep1_lnd, ep1_lnT, ep1_Et, &
            e00_PR, ep1_PR, &           
            e00_PL, em1_PL, &
            e00_uR, ep1_uR, &           
            e00_uL, em1_uL, &
            e00_rhoR, ep1_rhoR, &          
            e00_rhoL, em1_rhoL, &
            e00_csR, ep1_csR, &         
            e00_csL, em1_csL, &            
            e00_w, ep1_w, &
            ierr)
         use hydro_reconstruct, only: get_G
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         real(dp), intent(in) :: P_surf ! only used if k==1
         real(dp), intent(in) :: &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT
         logical, intent(in) :: skip_partials
         real(dp), intent(out) :: &
            residual, &
            e00_u, em1_u, ep1_u, &
            e00_lum, ep1_lum, &
            e00_lnR, ep1_lnR, &
            e00_lnd, e00_lnT, e00_Et, &
            em1_lnd, em1_lnT, em1_Et, &
            ep1_lnd, ep1_lnT, ep1_Et, &
            e00_PR, ep1_PR, &           
            e00_PL, em1_PL, &
            e00_uR, ep1_uR, &           
            e00_uL, em1_uL, &
            e00_rhoR, ep1_rhoR, &          
            e00_rhoL, em1_rhoL, &
            e00_csR, ep1_csR, &         
            e00_csL, em1_csL, &      
            e00_w, ep1_w
         integer, intent(out) :: ierr
      
         integer :: j, nz, i_du_dt
         real(qp) :: &
            flux_out, flux_in, momflux00, momfluxp1, rR, rL, P, &
            geometry_source, gravity_source, diffusion_source, &
            sources, dm, mL, mR, dudt_ex_qp, &
            gsL, gsR, u_p1, u_00, u_m1, sigp1, sig00            
         real(dp) :: &
            G00, dG00_dlnR, dgsR_dlnRR, Gp1, dGp1_dlnR, dgsL_dlnRL, &
            dG00_dw, dGp1_dw, dgsR_dwR, dgsL_dwL, &
            dudt_expected, ds_dlnR00, ds_dlnRp1, ds_dlnd00, ds_dlnT00, ds_dw00, ds_dwp1, &
            d_momflux00_dlnR00, d_momflux00_dL00, &
            d_momflux00_du00, d_momflux00_dum1, &
            d_momflux00_dlnd00, d_momflux00_dlndm1, &
            d_momflux00_dlnT00, d_momflux00_dlnTm1, &
            d_momfluxp1_dlnRp1, d_momfluxp1_dLp1, &
            d_momfluxp1_dup1, d_momfluxp1_du00, &
            d_momfluxp1_dlndp1, d_momfluxp1_dlnd00, &
            d_momfluxp1_dlnTp1, d_momfluxp1_dlnT00, &
            d_momflux00_dPR_00, d_momflux00_dPL_m1, &
            d_momflux00_duR_00, d_momflux00_duL_m1, &
            d_momflux00_drhoR_00, d_momflux00_drhoL_m1, &
            d_momflux00_dcsR_00, d_momflux00_dcsL_m1, &
            d_momfluxp1_dPR_p1, d_momfluxp1_dPL_00, &
            d_momfluxp1_duR_p1, d_momfluxp1_duL_00, &
            d_momfluxp1_drhoR_p1, d_momfluxp1_drhoL_00, &
            d_momfluxp1_dcsR_p1, d_momfluxp1_dcsL_00, &
            d_momflux00_dw00, d_momfluxp1_dwp1, &
            Uq, dUq_dlnR00, dUq_dlnRp1, &
            dUq_dlnd00, dUq_dlnT00, dUq_dEt00, &
            dUq_dlndm1, dUq_dlnTm1, dUq_dEtm1, &
            dUq_dlndp1, dUq_dlnTp1, dUq_dEtp1, &
            dUq_dum1, dUq_du00, dUq_dup1, &
            d_dudt_dEt00, d_dudt_dEtp1, d_dudt_dEtm1, &
            d_dudt_dlnR00, d_dudt_dlnRp1, d_dudt_dL00, d_dudt_dLp1, &
            d_dudt_dw00, d_dudt_dwp1, &
            d_dudt_dum1, d_dudt_du00, d_dudt_dup1, &
            d_dudt_dlndm1, d_dudt_dlnd00, d_dudt_dlndp1, &
            d_dudt_dlnTm1, d_dudt_dlnT00, d_dudt_dlnTp1, &
            d_geom_source_dlnRL, d_geom_source_dlnRR, &
            d_geom_source_dlnd, d_geom_source_dlnT, &
            d_grav_source_dlnRL, d_grav_source_dlnRR, &
            d_grav_source_dwL, d_grav_source_dwR, &
            kap, F_avg, AL, AR, ds_dL00, ds_dLp1, &         
            FR, rad_force_R, d_rfR_dLR, d_rfR_dlnRR, &
            d_rfR_dlnd00, d_rfR_dlnT00, d_rfR_dlndm1, d_rfR_dlnTm1, &             
            kap_faceR, d_kap_faceR_dlnd00, d_kap_faceR_dlnT00, &
            d_kap_faceR_dlndm1, d_kap_faceR_dlnTm1, &         
            FL, rad_force_L, d_rfL_dLL, d_rfL_dlnRL, &
            d_rfL_dlnd00, d_rfL_dlnT00, d_rfL_dlndp1, d_rfL_dlnTp1, &             
            kap_faceL, d_kap_faceL_dlnd00, d_kap_faceL_dlnT00, &
            d_kap_faceL_dlndp1, d_kap_faceL_dlnTp1, &
            ds_dlndm1, ds_dlnTm1, ds_dlndp1, ds_dlnTp1, &            
            alfa, beta, ie_plus_ke, scal, dt, d_dlnd, d_dlnT, &
            d_diff_source_dum1, d_diff_source_du00, d_diff_source_dup1, &
            full_on, full_off, logRho, rmid, &
            dudt_actual, d_dudt_actual_du, dudt_factor, alpha, &
            v_drag, drag_factor, drag_fraction
         logical :: dbg, do_diffusion, test_partials

         include 'formats'
         
         if (s% use_other_momentum) &
            stop 'Riemann dudt does not support use_other_momentum'
         if (s% use_other_momentum_implicit) &
            stop 'Riemann dudt does not support use_other_momentum_implicit'
         if (s% use_mass_corrections) &
            stop 'Riemann dudt does not support use_mass_corrections'

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
            
         ierr = 0
         nz = s% nz
         i_du_dt = s% i_du_dt
         
         residual = 0
         e00_u = 0; em1_u = 0; ep1_u = 0
         e00_lum = 0; ep1_lum = 0
         e00_lnR = 0; ep1_lnR = 0
         e00_lnd = 0; e00_lnT = 0
         em1_lnd = 0; em1_lnT = 0
         ep1_lnd = 0; ep1_lnT = 0

         e00_uL = 0; em1_uL = 0
         e00_uR = 0; ep1_uR = 0
         e00_PL = 0; em1_PL = 0
         e00_PR = 0; ep1_PR = 0
         e00_rhoL = 0; em1_rhoL = 0
         e00_rhoR = 0; ep1_rhoR = 0
         e00_csL = 0; em1_csL = 0
         e00_csR = 0; ep1_csR = 0
         
         dbg = .false.
         if (k == 1) then
            call eval_surf_momentum_flux(s, P_surf, &
               dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT, &
               momflux00, d_momflux00_dlnR00, d_momflux00_dL00, &
               d_momflux00_du00, d_momflux00_dum1, &
               d_momflux00_dlnd00, d_momflux00_dlndm1, &
               d_momflux00_dlnT00, d_momflux00_dlnTm1, &
               d_momflux00_dw00,ierr)   
            if (ierr /= 0) return
         else
            call eval1_momentum_flux(s, k, &
               momflux00, d_momflux00_dlnR00, d_momflux00_dL00, &
               d_momflux00_du00, d_momflux00_dum1, &
               d_momflux00_dlnd00, d_momflux00_dlndm1, &
               d_momflux00_dlnT00, d_momflux00_dlnTm1, &
               d_momflux00_dw00,ierr)   
            if (ierr /= 0) return
         end if
         
         call eval1_momentum_flux(s, k+1, &
            momfluxp1, d_momfluxp1_dlnRp1, d_momfluxp1_dLp1, &
            d_momfluxp1_dup1, d_momfluxp1_du00, &
            d_momfluxp1_dlndp1, d_momfluxp1_dlnd00, &
            d_momfluxp1_dlnTp1, d_momfluxp1_dlnT00, &
            d_momfluxp1_dwp1, ierr)   
         if (ierr /= 0) return
                  
         dt = s% dt
         dm = s% dm(k)
         flux_out = momflux00
         rR = s% r(k)
         mR = s% m(k)
         if (k == nz) then
            rL = s% R_center
            mL = s% M_center
            flux_in = 0
         else
            rL = s% r(k+1)
            mL = s% m(k+1)
            flux_in = momfluxp1
         end if
         
         P = s% P(k)
                  
         if (k == nz) then 
            ! no flux in from left, so only have geometry source on right
            ! this matters for cases with R_center > 0.
            geometry_source = 4*pi*P*rR*rR
            d_geom_source_dlnRR = 8*pi*P*rR*rR
            d_geom_source_dlnRL = 0
         else
            geometry_source = 4*pi*P*(rR*rR - rL*rL)
            d_geom_source_dlnRR = 8*pi*P*rR*rR
            d_geom_source_dlnRL = -8*pi*P*rL*rL
         end if
         d_geom_source_dlnd = geometry_source*s% chiRho_for_partials(k)
         d_geom_source_dlnT = geometry_source*s% chiT_for_partials(k)
         
         if (s% zero_gravity) then
            gravity_source = 0
            d_grav_source_dlnRL = 0
            d_grav_source_dlnRR = 0
            d_grav_source_dwL = 0
            d_grav_source_dwR = 0
         else
            ! left 1/2 of dm gets gravity force at left face
            ! right 1/2 of dm gets gravity force at right face.
            ! this form is to match the gravity force equilibrium reconstruction.
            call get_G(s, k, G00, dG00_dlnR, dG00_dw)
            gsR = -G00*mR*0.5d0*dm/(rR*rR)
            dgsR_dlnRR = -2d0*gsR + gsR*dG00_dlnR/G00
            dgsR_dwR = gsR*dG00_dw/G00
            if (k == nz) then
               gsL = 0d0
               dgsL_dlnRL = 0
               dgsL_dwL = 0
            else
               call get_G(s, k+1, Gp1, dGp1_dlnR, dGp1_dw)
               gsL = -Gp1*mL*0.5d0*dm/(rL*rL)
               dgsL_dlnRL = -2d0*gsL + gsL*dGp1_dlnR/Gp1
               dgsL_dwL = gsL*dGp1_dw/Gp1
            end if
            gravity_source = gsL + gsR ! total gravitational force on cell
            d_grav_source_dlnRL = dgsL_dlnRL
            d_grav_source_dlnRR = dgsR_dlnRR
            d_grav_source_dwL = dgsL_dwL
            d_grav_source_dwR = dgsR_dwR
         end if

         sources = geometry_source + gravity_source

         do_diffusion = s% RTI_flag .and. s% dudt_RTI_diffusion_factor > 0d0
         if (do_diffusion) then ! add diffusion source term to dudt
            if (k < nz) then
               sigp1 = s% dudt_RTI_diffusion_factor*s% sig_RTI(k+1)
               u_p1 = s% u(k+1)
            else
               sigp1 = 0
               u_p1 = 0
            end if
            if (k > 1) then
               sig00 = s% dudt_RTI_diffusion_factor*s% sig_RTI(k)
               u_m1 = s% u(k-1)
            else
               sig00 = 0
               u_m1 = 0
            end if
            u_00 = s% u(k)
            diffusion_source = sig00*(u_m1 - u_00) - sigp1*(u_00 - u_p1)
            d_diff_source_dum1 = sig00
            d_diff_source_du00 = -(sig00 + sigp1)
            d_diff_source_dup1 = sigp1
         else
            diffusion_source = 0d0
            d_diff_source_dum1 = 0d0
            d_diff_source_du00 = 0d0
            d_diff_source_dup1 = 0d0
         end if
         sources = sources + diffusion_source
         s% dudt_RTI(k) = diffusion_source/dm
         
         dudt_ex_qp = (flux_in - flux_out + sources)/dm
         dudt_expected = dudt_ex_qp

         drag_factor = s% v_drag_factor
         v_drag = s% v_drag
         if (s% q(k) < s% q_for_v_drag_full_off) then
            drag_fraction = 0d0
         else if (s% q(k) > s% q_for_v_drag_full_on) then
            drag_fraction = 1d0
         else
            drag_fraction = (s% q(k) - s% q_for_v_drag_full_off)&
                               /(s% q_for_v_drag_full_on - s% q_for_v_drag_full_off)
         end if
         drag_factor = drag_factor*drag_fraction

         if (drag_factor > 0d0) then
            if (s% u(k) > v_drag) then
               dudt_expected = dudt_expected - drag_factor*pow2(s% u(k) - v_drag)/s% r(k)
            else if (s% u(k) < -v_drag) then
               dudt_expected = dudt_expected + drag_factor*pow2(s% u(k) + v_drag)/s% r(k)
            end if
         end if
         
         dudt_factor = 1d0
         
         ! make residual units be relative difference in energy
         ie_plus_ke = s% energy_start(k) + 0.5d0*s% u_start(k)*s% u_start(k)
         scal = dt*max(abs(s% u_start(k)),s% csound_start(k))/ie_plus_ke
         if (k == 1) scal = scal*1d-2
         
         
         dudt_actual = s% du_dt(k)
         d_dudt_actual_du = 1d0/dt
         
         residual = (dudt_expected - dudt_actual)*scal

         if (test_partials) then
            s% solver_test_partials_val = flux_in - flux_out
         end if
         
         if (is_bad(residual)) then
            ierr = -1
            return
!$omp critical (dudt_eqn)
            write(*,2) 'residual', k, residual
            write(*,2) 's% L(k+1)', k+1, s% L(k+1)
            write(*,2) 's% L(k)', k, s% L(k)
            write(*,2) 's% opacity(k)', k, s% opacity(k)
            write(*,2) 'rL', k, rL
            write(*,2) 'rR', k, rR
            stop 'do1_dudt_residual_and_partials'
!$omp end critical (dudt_eqn)
         end if
         
         if (skip_partials) return

         ! partials of geometry_source + gravity_source
         ds_dlnR00 = d_geom_source_dlnRR + d_grav_source_dlnRR
         ds_dlnRp1 = d_geom_source_dlnRL + d_grav_source_dlnRL
         ds_dlnd00 = d_geom_source_dlnd
         ds_dlnT00 = d_geom_source_dlnT
         ds_dw00 = d_grav_source_dwR
         ds_dwp1 = d_grav_source_dwL
         ds_dL00 = 0d0
         ds_dLp1 = 0d0
         ds_dlndm1 = 0d0
         ds_dlnTm1 = 0d0
         ds_dlndp1 = 0d0
         ds_dlnTp1 = 0d0

         d_dudt_dL00 = dudt_factor*(0 - d_momflux00_dL00 + ds_dL00)/dm
         d_dudt_dLp1 = dudt_factor*(d_momfluxp1_dLp1 - 0 + ds_dLp1)/dm
         
         d_dudt_dlnR00 = dudt_factor*(0 - d_momflux00_dlnR00 + ds_dlnR00)/dm
         d_dudt_dlnRp1 = dudt_factor*(d_momfluxp1_dlnRp1 - 0 + ds_dlnRp1)/dm  
         
         d_dudt_dw00 = dudt_factor*(0 - d_momflux00_dw00 + ds_dw00)/dm
         d_dudt_dwp1 = dudt_factor*(d_momfluxp1_dwp1 - 0 + ds_dwp1)/dm  

         if (drag_factor > 0d0) then
            if (s% u(k) > v_drag) then
               d_dudt_dlnR00 = d_dudt_dlnR00 + drag_factor*pow2(s% u(k) - v_drag)/s% r(k)
            else if (s% u(k) < -v_drag) then
               d_dudt_dlnR00 = d_dudt_dlnR00 - drag_factor*pow2(s% u(k) + v_drag)/s% r(k)
            end if
         end if
               
         d_dudt_dum1 = &
            dudt_factor*(0 - d_momflux00_dum1 + 0 + d_diff_source_dum1)/dm
         d_dudt_du00 = &
            dudt_factor*(d_momfluxp1_du00 - d_momflux00_du00 + 0 + d_diff_source_du00)/dm
         d_dudt_dup1 = &
            dudt_factor*(d_momfluxp1_dup1 - 0 + 0 + d_diff_source_dup1)/dm  

         if (drag_factor > 0d0) then
            if (s% u(k) > v_drag) then
               d_dudt_du00 = d_dudt_du00 - 2d0*drag_factor*(s% u(k) - v_drag)/s% r(k)
            else if (s% u(k) < -v_drag) then
               d_dudt_du00 = d_dudt_du00 + 2d0*drag_factor*(s% u(k) + v_drag)/s% r(k)
            end if
         end if
                   
         d_dudt_dlndm1 = dudt_factor*(0 - d_momflux00_dlndm1 + ds_dlndm1)/dm
         d_dudt_dlnd00 = &
            dudt_factor*(d_momfluxp1_dlnd00 - d_momflux00_dlnd00 + ds_dlnd00)/dm
         d_dudt_dlndp1 = dudt_factor*(d_momfluxp1_dlndp1 - 0 + ds_dlndp1)/dm         
         
         d_dudt_dlnTm1 = dudt_factor*(0 - d_momflux00_dlnTm1 + ds_dlnTm1)/dm
         d_dudt_dlnT00 = &
            dudt_factor*(d_momfluxp1_dlnT00 - d_momflux00_dlnT00 + ds_dlnT00)/dm
         d_dudt_dlnTp1 = dudt_factor*(d_momfluxp1_dlnTp1- 0 + ds_dlnTp1)/dm         
            
         d_dudt_dEtp1 = 0
         d_dudt_dEt00 = 0
         d_dudt_dEtm1 = 0

         e00_Et = d_dudt_dEt00*scal ! 
         if (k > 1) em1_Et = d_dudt_dEtm1*scal ! 
         if (k < nz) ep1_Et = d_dudt_dEtp1*scal ! 

         e00_u = d_dudt_du00*scal - d_dudt_actual_du*scal ! 
         if (k > 1) em1_u = d_dudt_dum1*scal ! 
         if (k < nz) ep1_u = d_dudt_dup1*scal ! 
      
         e00_lum = d_dudt_dL00*scal
         ep1_lum = d_dudt_dLp1*scal

         e00_lnR = d_dudt_dlnR00*scal ! 
         if (k < nz) ep1_lnR = d_dudt_dlnRp1*scal ! 
            
         e00_lnd = d_dudt_dlnd00*scal ! 
         e00_lnT = d_dudt_dlnT00*scal ! 

         e00_w = d_dudt_dw00*scal ! 
         if (k < nz) ep1_w = d_dudt_dwp1*scal ! 
         
         e00_PR = 0; ep1_PR = 0            
         e00_PL = 0; em1_PL = 0
         e00_uR = 0; ep1_uR = 0            
         e00_uL = 0; em1_uL = 0
         e00_rhoR = 0; ep1_rhoR = 0            
         e00_rhoL = 0; em1_rhoL = 0
         e00_csR = 0; ep1_csR = 0            
         e00_csL = 0; em1_csL = 0

         if (k > 1) then            
            em1_lnd = d_dudt_dlndm1*scal ! 
            em1_lnT = d_dudt_dlnTm1*scal ! 
         end if
         
         if (k < nz) then           
            ep1_lnd = d_dudt_dlndp1*scal ! 
            ep1_lnT = d_dudt_dlnTp1*scal ! 
         end if
   
         !e00_u = ok; em1_u = 0; ep1_u = 0
         !e00_lum = ok; ep1_lum = 0
         !e00_lnR = ok; ep1_lnR = 0
         !e00_lnd = ok; e00_lnT = BAD
         !em1_lnd = 0; em1_lnT = 0
         !ep1_lnd = 0; ep1_lnT = 0
         
         if (test_partials) then
            s% solver_test_partials_var = s% i_w_div_wc
            s% solver_test_partials_dval_dx = d_grav_source_dwR/dm
            !write(*,2) 'Uq', k, Uq
            !write(*,2) 'dudt_expected', k, dudt_expected
            !write(*,2) 'dudt_actual', k, dudt_actual
            !write(*,2) 'geometry_source', k, geometry_source
            !write(*,2) 'gravity_source', k, gravity_source
            !write(*,2) 'sources', k, sources
            !write(*,2) 'flux_in - flux_out', k, flux_in - flux_out
            !write(*,2) 'flux_in', k, flux_in
            !write(*,2) 'flux_out', k, flux_out
            !write(*,2) 'drag_factor', k, drag_factor
            !write(*,2) 'dudt_factor', k, dudt_factor
            !write(*,2) 'scal', k, scal
            !write(*,2) 's% u(k)', k, s% u(k)
            !write(*,2) 's% u_start(k)', k, s% u_start(k)
            !!write(*,2) 's% Et_Eq(k)', k, s% Et_Eq(k)
            !write(*,*)
         end if
  
      end subroutine do1_dudt_residual_and_partials
      
      
      subroutine store1_dudt_residual_and_partials( &
            s, k, xscale, equ, skip_partials, nvar, &
            residual, &
            e00_u, em1_u, ep1_u, &
            e00_lum, ep1_lum, &
            e00_lnR, ep1_lnR, &
            e00_lnd, e00_lnT, e00_Et, &
            em1_lnd, em1_lnT, em1_Et, &
            ep1_lnd, ep1_lnT, ep1_Et, &
            e00_PR, ep1_PR, &           
            e00_PL, em1_PL, &
            e00_uR, ep1_uR, &           
            e00_uL, em1_uL, &
            e00_rhoR, ep1_rhoR, &          
            e00_rhoL, em1_rhoL, &
            e00_csR, ep1_csR, &         
            e00_csL, em1_csL, &            
            e00_w, ep1_w, &
            ierr)
         type (star_info), pointer :: s         
         integer, intent(in) :: k
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: skip_partials
         integer, intent(in) :: nvar
         real(dp), intent(in) :: &
            residual, &
            e00_u, em1_u, ep1_u, &
            e00_lum, ep1_lum, &
            e00_lnR, ep1_lnR, &
            e00_lnd, e00_lnT, e00_Et, &
            em1_lnd, em1_lnT, em1_Et, &
            ep1_lnd, ep1_lnT, ep1_Et, &
            e00_PR, ep1_PR, &           
            e00_PL, em1_PL, &
            e00_uR, ep1_uR, &           
            e00_uL, em1_uL, &
            e00_rhoR, ep1_rhoR, &          
            e00_rhoL, em1_rhoL, &
            e00_csR, ep1_csR, &         
            e00_csL, em1_csL, &
            e00_w, ep1_w
         integer, intent(out) :: ierr
         
         integer :: nz, i_du_dt, i_u
         include 'formats'
         
         ierr = 0
         nz = s% nz
         i_du_dt = s% i_du_dt
         i_u = s% i_u

         equ(i_du_dt, k) = residual
         s% u_residual(k) = residual
         
         if (skip_partials) return

         call e00(s, xscale, i_du_dt, i_u, k, nvar, e00_u)
         if (k > 1) call em1(s, xscale, i_du_dt, i_u, k, nvar, em1_u)
         if (k < nz) call ep1(s, xscale, i_du_dt, i_u, k, nvar, ep1_u)
      
         call e00(s, xscale, i_du_dt, s% i_lum, k, nvar, e00_lum)
         if (k < nz) call ep1(s, xscale, i_du_dt, s% i_lum, k, nvar, ep1_lum)
         
         call e00(s, xscale, i_du_dt, s% i_lnR, k, nvar, e00_lnR)
         if (k < nz) call ep1(s, xscale, i_du_dt, s% i_lnR, k, nvar, ep1_lnR)

         call e00(s, xscale, i_du_dt, s% i_lnd, k, nvar, e00_lnd)
         if (s% do_struct_thermo) &
            call e00(s, xscale, i_du_dt, s% i_lnT, k, nvar, e00_lnT)
         
         if (s% w_div_wc_flag) then
            call e00(s, xscale, i_du_dt, s% i_w_div_wc, k, nvar, e00_w)
            if (k < nz) call ep1(s, xscale, s% i_du_dt, s% i_w_div_wc, k, nvar, ep1_w)
         end if
         
         if (k > 1) then            
            call em1(s, xscale, i_du_dt, s% i_lnd, k, nvar, em1_lnd)
            if (s% do_struct_thermo) &
               call em1(s, xscale, i_du_dt, s% i_lnT, k, nvar, em1_lnT)
         end if
         
         if (k < nz) then           
            call ep1(s, xscale, i_du_dt, s% i_lnd, k, nvar, ep1_lnd)
            if (s% do_struct_thermo) &
               call ep1(s, xscale, i_du_dt, s% i_lnT, k, nvar, ep1_lnT)
         end if
         
      end subroutine store1_dudt_residual_and_partials
      
      
      subroutine do1_Riemann_energy_eqn( &
            s, k, P_surf, dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT, &
            xscale, equ, skip_partials, do_chem, nvar, ierr)
         use eos_def, only: i_grad_ad
         use star_utils, only: cell_start_specific_PE, cell_specific_PE
         type (star_info), pointer :: s         
         integer, intent(in) :: k, nvar
         real(dp), intent(in) :: P_surf ! only used if k==1
         real(dp), intent(in) :: &
            dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: skip_partials, do_chem
         integer, intent(out) :: ierr

         real(qp) :: &
            eflux00, efluxp1, dt, dm, flux_out, flux_in, u_new, u_old, &
            dKE_dt, mC, rC_new, rC_start, G, dPE_dt, sources, dedt_expected_qp, &
            energy_qp, e0_qp, f_qp, de_div_e0_expected, de_div_e0_actual, &
            dEt_dt, sig00, sigp1, diffusion_eps, e_m1, e_p1
         real(dp) :: &
            d_eflux_dL00, d_eflux00_dlnR00, &
            d_eflux00_du00, d_eflux00_dum1, &
            d_eflux00_dlnd00, d_eflux00_dlndm1, &
            d_eflux00_dlnT00, d_eflux00_dlnTm1, &
            d_eflux_dLp1, d_efluxp1_dlnRp1, &
            d_efluxp1_dup1, d_efluxp1_du00, &
            d_efluxp1_dlndp1, d_efluxp1_dlnd00, &
            d_efluxp1_dlnTp1, d_efluxp1_dlnT00, &
            d_eflux00_dPR_00, d_eflux00_dPL_m1, &
            d_eflux00_duR_00, d_eflux00_duL_m1, &
            d_eflux00_drhoR_00, d_eflux00_drhoL_m1, &
            d_eflux00_dcsR_00, d_eflux00_dcsL_m1, &
            d_efluxp1_dPR_p1, d_efluxp1_dPL_00, &
            d_efluxp1_duR_p1, d_efluxp1_duL_00, &
            d_efluxp1_drhoR_p1, d_efluxp1_drhoL_00, &
            d_efluxp1_dcsR_p1, d_efluxp1_dcsL_00, &
            d_eflux00_dw00, d_efluxp1_dwp1, &
            e00_PR, ep1_PR, e00_PL, em1_PL, &
            e00_uR, ep1_uR, e00_uL, em1_uL, &
            e00_rhoR, ep1_rhoR, e00_rhoL, em1_rhoL, &
            e00_csR, ep1_csR, e00_csL, em1_csL, &         
            alpha, energy, e0, f, rR, rL, dedt_expected, d_dKEdt_du, E_src_sink, &
            d_rC_new_dlnR00, d_rC_new_dlnRp1, emin_start, vt0, dPE_dlnR00, dPE_dlnRp1, &
            d_dPEdt_dlnR00, d_dPEdt_dlnRp1, eps_nuc, non_nuc_neu, PE_start, PE_new, &
            e_div_e0, d_dedt_dL00, d_dedt_dLp1, d_dedt_dlnR00, d_dedt_dlnRp1, d_dedt_dw00, d_dedt_dwp1, &
            d_dedt_dum1, d_dedt_du00, d_dedt_dup1, Et0, Et, &
            d_dedt_dlndm1, d_dedt_dlnd00, d_dedt_dlndp1, d_dEtdt_dEt, &
            d_dedt_dlnTm1, d_dedt_dlnT00, d_dedt_dlnTp1, &
            d_dedt_dlnd_c_E, d_dedt_dE_c_Rho, d_dedt_dEt, &
            d_de_div_e0_dlnd_c_E, d_de_div_e0_dlnE_c_Rho, &
            num_dequ, num_dlnd, num_deflux00, num_defluxp1, &
            num_dP_face, num_du_face, dlgE, dequ, dlgT, dlgRho, &
            d_Etotal_dt, d_dlnd, d_dlnT, &
            d_diff_eps_dlnE_m1, d_diff_eps_dlnE_00, d_diff_eps_dlnE_p1, &
            diffusion_factor, diffusion_eps_in, diffusion_eps_out, scal, term
         integer :: nz, j, i, i_dlnE_dt, i_lnd, i_lnT, i_lnR, i_lum, i_u, i_w_div_wc
         logical :: do_diffusion, test_partials, doing_op_split_burn
         
         include 'formats'

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.
         
         ierr = 0
         i_dlnE_dt = s% i_dlnE_dt
         nz = s% nz
         E_src_sink = 0

         doing_op_split_burn = s% op_split_burn .and. &
            s% T_start(k) >= s% op_split_burn_min_T
         
         dt = 1d0/s% dVARdot_dVAR
         
         dm = s% dm(k)
         energy = s% energy(k)
         energy_qp = energy
         e0 = s% energy_start(k)
         e0_qp = e0
         e_div_e0 = energy/e0
         f_qp = dt/e0_qp ! factor to convert dedt to de_div_e0
         f = f_qp

         if (k==1) then
            call eval_surf_energy_flux(s, P_surf, &
               dlnPsurf_dL, dlnPsurf_dlnR, dlnPsurf_dlnd, dlnPsurf_dlnT, &
               eflux00, d_eflux_dL00, d_eflux00_dlnR00, &
               d_eflux00_du00, d_eflux00_dum1, &
               d_eflux00_dlnd00, d_eflux00_dlndm1, &
               d_eflux00_dlnT00, d_eflux00_dlnTm1, &
               d_eflux00_dw00, ierr)   
            if (ierr /= 0) return
         else
            call eval1_energy_flux(s, k, &
               eflux00, d_eflux_dL00, d_eflux00_dlnR00, &
               d_eflux00_du00, d_eflux00_dum1, &
               d_eflux00_dlnd00, d_eflux00_dlndm1, &
               d_eflux00_dlnT00, d_eflux00_dlnTm1, &
               d_eflux00_dw00, ierr)   
            if (ierr /= 0) return
         end if

         call eval1_energy_flux(s, k+1, &
            efluxp1, d_eflux_dLp1, d_efluxp1_dlnRp1, &
            d_efluxp1_dup1, d_efluxp1_du00, &
            d_efluxp1_dlndp1, d_efluxp1_dlnd00, &
            d_efluxp1_dlnTp1, d_efluxp1_dlnT00, &
            d_efluxp1_dwp1, ierr)   
         if (ierr /= 0) return
         
         flux_out = eflux00
         flux_in = efluxp1
         
         dKE_dt = 0d0
         dPE_dt = 0d0
         sources = 0d0
         do_diffusion = .false.
         
         rR = s% r(k)
         if (k == nz) then
            rL = s% R_center
         else
            rL = s% r(k+1)
         end if
         
         ! rate of change in specific KE (erg/g/s)
         u_new = s% u(k)
         u_old = s% u_start(k)
         dKE_dt = 0.5d0*(u_new*u_new - u_old*u_old)/dt ! erg/g/s
         d_dKEdt_du = u_new/dt
         
         ! rate of change in specific PE (erg/g/s)
         PE_start = cell_start_specific_PE(s,k)
         PE_new = cell_specific_PE(s,k,dPE_dlnR00,dPE_dlnRp1)
         dPE_dt = (PE_new - PE_start)/dt ! erg/g/s
         d_dPEdt_dlnR00 = dPE_dlnR00/dt
         d_dPEdt_dlnRp1 = dPE_dlnRp1/dt
         
         diffusion_factor = s% dedt_RTI_diffusion_factor
         do_diffusion = s% RTI_flag .and. diffusion_factor > 0d0
         if (do_diffusion) then ! add diffusion source term to dedt         
            if (k < s% nz) then
               if (s% alpha_RTI(k) > 1d-10 .and. k > 1) then
                  emin_start = min( &
                     s% energy_start(k+1), s% energy_start(k), s% energy_start(k-1))
                  if (emin_start < 5d0*s% RTI_energy_floor) then
                     diffusion_factor = diffusion_factor* &
                        (1d0 + (5d0*s% RTI_energy_floor - emin_start)/emin_start)
                     !write(*,2) 'boosted energy diff factor', k+1, diffusion_factor
                  end if
               end if
               sigp1 = diffusion_factor*s% sig_RTI(k+1)
               e_p1 = s% energy(k+1)
            else
               sigp1 = 0
               e_p1 = 0
            end if
            if (k > 1) then
               sig00 = diffusion_factor*s% sig_RTI(k)
               e_m1 = s% energy(k-1)
            else
               sig00 = 0
               e_m1 = 0
            end if
            diffusion_eps_in = sigp1*(e_p1 - energy_qp)/dm
            diffusion_eps_out = sig00*(energy_qp - e_m1)/dm
            diffusion_eps = diffusion_eps_in - diffusion_eps_out
            d_diff_eps_dlnE_m1 = sig00*e_m1/dm
            d_diff_eps_dlnE_00 = -(sig00 + sigp1)*energy_qp/dm
            d_diff_eps_dlnE_p1 = sigp1*e_p1/dm
         else
            diffusion_eps_in = 0d0
            diffusion_eps_out = 0d0
            diffusion_eps = 0d0
            d_diff_eps_dlnE_m1 = 0d0
            d_diff_eps_dlnE_00 = 0d0
            d_diff_eps_dlnE_p1 = 0d0
            e_p1 = 0
            e_m1 = 0
         end if
         sources = sources + diffusion_eps
         s% dedt_RTI(k) = diffusion_eps
      
         if (s% nonlocal_NiCo_decay_heat .or. &
               s% gamma_law_hydro > 0 .or. doing_op_split_burn) then
            eps_nuc = 0 ! get eps_nuc from extra_heat instead
         else
            eps_nuc = s% eps_nuc(k)
         end if

         if (s% gamma_law_hydro > 0) then
            non_nuc_neu = 0d0
         else
            non_nuc_neu = 0.5d0*(s% non_nuc_neu_start(k) + s% non_nuc_neu(k))
         end if
         
         s% eps_heat(k) = &
            eps_nuc - non_nuc_neu + s% extra_heat(k) + s% irradiation_heat(k)
         if (s% mstar_dot /= 0d0) &
            s% eps_heat(k) = s% eps_heat(k) + s% eps_mdot(k)
         if (s% do_element_diffusion .and. s% do_WD_sedimentation_heating) &
            s% eps_heat(k) = s% eps_heat(k) + s% eps_WD_sedimentation(k)
         if (s% do_element_diffusion .and. s% do_diffusion_heating) &
            s% eps_heat(k) = s% eps_heat(k) + s% eps_diffusion(k)
         if (s% do_conv_premix .and. s% do_premix_heating) &
            s% eps_heat(k) = s% eps_heat(k) + s% eps_pre_mix(k)
         sources = sources + s% eps_heat(k)
         
         E_src_sink = E_src_sink + s% eps_heat(k)
         dEt_dt = 0d0
         d_dEtdt_dEt = 0d0
         
         d_Etotal_dt = (flux_in - flux_out)/dm + sources
         ! Etotal = IE + KE + PE + Et,
         ! so energy = IE = Etotal - KE - PE - Et,
         ! d_internal_enegy = d_Etotal - dKE - dPE - dEt_dt
         dedt_expected_qp = d_Etotal_dt - dKE_dt - dPE_dt - dEt_dt
         dedt_expected = dedt_expected_qp ! erg/g/s     

         de_div_e0_expected = dedt_expected*f
         de_div_e0_actual = (energy_qp - e0_qp)/e0_qp
         
         if (k > 1) then
            scal = 1d0
         else
            scal = 1d-6
         end if
         
         equ(i_dlnE_dt, k) = scal*(de_div_e0_expected - de_div_e0_actual)
         s% E_residual(k) = equ(i_dlnE_dt, k)
         s% ergs_error(k) = dm*e0*(de_div_e0_actual - de_div_e0_expected) 
            ! this gives the right sign for reported errors
         
         if (is_bad(equ(i_dlnE_dt, k))) then
            ierr = -1
            if (.not. s% report_ierr) return
!$omp critical (hydro_reimann_crit1)
            write(*,2) 'equ(i_dlnE_dt, k)', k, equ(i_dlnE_dt, k)
            write(*,2) 'de_div_e0_expected', k, de_div_e0_expected
            write(*,2) 'de_div_e0_actual', k, de_div_e0_actual
            write(*,2) 'dedt_expected_qp', k, dedt_expected_qp
            write(*,2) 's% dVARdot_dVAR', k, s% dVARdot_dVAR
            write(*,2) 'dt', k, dt
            write(*,2) 'e0', k, e0
            write(*,2) 'f_qp', k, f_qp
            write(*,2) 'd_Etotal_dt', k, d_Etotal_dt
            write(*,2) 'dKE_dt', k, dKE_dt
            write(*,2) 'dPE_dt', k, dPE_dt
            write(*,2) 'dEt_dt', k, dEt_dt
            write(*,2) 'flux_in', k, flux_in
            write(*,2) 'flux_out', k, flux_out
            write(*,2) 'dm', k, dm
            write(*,2) 'sources', k, sources
            write(*,*)
            stop 'do1_Riemann_energy_eqn'
!$omp end critical (hydro_reimann_crit1)
         end if

         if (test_partials) then
            s% solver_test_partials_val = dedt_expected ! equ(i_dlnE_dt, k)
         end if
         
         if (skip_partials) return

         i_lnd = s% i_lnd
         i_lnT = s% i_lnT
         i_lnR = s% i_lnR
         i_lum = s% i_lum
         i_u = s% i_u
         i_w_div_wc = s% i_w_div_wc ! used when including rotation with implicit correction factors

         ! partial of -de_div_e0_actual
         call e00(s, xscale, i_dlnE_dt, i_lnd, k, nvar, -scal*s% dE_drho_for_partials(k)*s% rho(k)/e0)
         call e00(s, xscale, i_dlnE_dt, i_lnT, k, nvar, -scal*s% Cv_for_partials(k)*s% T(k)/e0)
         
         ! add partials of (flux_in - flux_out)/dm    
         d_dedt_dL00 = (0 - d_eflux_dL00)/dm
         d_dedt_dLp1 = (d_eflux_dLp1 - 0)/dm
         d_dedt_dlnR00 = (0 - d_eflux00_dlnR00)/dm
         d_dedt_dlnRp1 = (d_efluxp1_dlnRp1 - 0)/dm
         d_dedt_dum1 = (0 - d_eflux00_dum1)/dm
         d_dedt_du00 = (d_efluxp1_du00 - d_eflux00_du00)/dm
         d_dedt_dup1 = (d_efluxp1_dup1 - 0)/dm
         d_dedt_dlndm1 = (0 - d_eflux00_dlndm1)/dm
         d_dedt_dlnd00 = (d_efluxp1_dlnd00 - d_eflux00_dlnd00)/dm
         d_dedt_dlndp1 = (d_efluxp1_dlndp1 - 0)/dm
         d_dedt_dlnTm1 = (0 - d_eflux00_dlnTm1)/dm
         d_dedt_dlnT00 = (d_efluxp1_dlnT00 - d_eflux00_dlnT00)/dm
         d_dedt_dlnTp1 = (d_efluxp1_dlnTp1 - 0)/dm         
         d_dedt_dw00 = (0 - d_eflux00_dw00)/dm
         d_dedt_dwp1 = (d_efluxp1_dwp1 - 0)/dm
         
         ! add partials of extra heat
         d_dedt_dlnR00 = d_dedt_dlnR00 + s% d_extra_heat_dlnR00(k)
         d_dedt_dlnRp1 = d_dedt_dlnRp1 + s% d_extra_heat_dlnRp1(k)
         d_dedt_dlndm1 = d_dedt_dlndm1 + s% d_extra_heat_dlndm1(k)
         d_dedt_dlnd00 = d_dedt_dlnd00 + s% d_extra_heat_dlnd00(k)
         d_dedt_dlndp1 = d_dedt_dlndp1 + s% d_extra_heat_dlndp1(k)
         d_dedt_dlnTm1 = d_dedt_dlnTm1 + s% d_extra_heat_dlnTm1(k)
         d_dedt_dlnT00 = d_dedt_dlnT00 + s% d_extra_heat_dlnT00(k)
         d_dedt_dlnTp1 = d_dedt_dlnTp1 + s% d_extra_heat_dlnTp1(k)
                     
         ! add partials of eps_heat
         if (s% gamma_law_hydro == 0) then
            d_dedt_dlnT00 = d_dedt_dlnT00 - 0.5d0*s% d_nonnucneu_dlnT(k)
            d_dedt_dlnd00 = d_dedt_dlnd00 - 0.5d0*s% d_nonnucneu_dlnd(k)
            ! partials of eps_nuc
            if (.not. (s% nonlocal_NiCo_decay_heat .or. doing_op_split_burn)) then
               d_dedt_dlnT00 = d_dedt_dlnT00 + s% d_epsnuc_dlnT(k)
               d_dedt_dlnd00 = d_dedt_dlnd00 + s% d_epsnuc_dlnd(k)
               if (do_chem .and. s% dxdt_nuc_factor > 0d0) then
                  do j=1,s% nvar_chem
                     call e00(s, xscale, i_dlnE_dt, s% i_chem1+j-1, k, nvar, &
                              scal*s% d_epsnuc_dx(j,k)*f)
                  end do
               end if
            end if
         end if
         
         ! add partials of -dKE_dt
         d_dedt_du00 = d_dedt_du00 - d_dKEdt_du
         
         ! add partials of -dPE_dt
         d_dedt_dlnR00 = d_dedt_dlnR00 - d_dPEdt_dlnR00
         d_dedt_dlnRp1 = d_dedt_dlnRp1 - d_dPEdt_dlnRp1

         if (do_diffusion) then ! add partials of diffusion term
            d_dedt_dlnd00 = d_dedt_dlnd00 + d_diff_eps_dlnE_00/energy*s% dE_drho_for_partials(k)*s% rho(k)
            d_dedt_dlnT00 = d_dedt_dlnT00 + d_diff_eps_dlnE_00/energy*s% Cv_for_partials(k)*s% T(k)
            if (k > 1) then
               d_dedt_dlndm1 = d_dedt_dlndm1 + &
                  d_diff_eps_dlnE_m1/dble(e_m1)*s% dE_drho_for_partials(k-1)*s% rho(k-1)
               d_dedt_dlnTm1 = d_dedt_dlnTm1 + &
                  d_diff_eps_dlnE_m1/dble(e_m1)*s% Cv_for_partials(k-1)*s% T(k-1)
            end if
            if (k < nz) then
               d_dedt_dlndp1 = d_dedt_dlndp1 + &
                  d_diff_eps_dlnE_p1/dble(e_p1)*s% dE_drho_for_partials(k+1)*s% rho(k+1)
               d_dedt_dlnTp1 = d_dedt_dlnTp1 + &
                  d_diff_eps_dlnE_p1/dble(e_p1)*s% Cv_for_partials(k+1)*s% T(k+1)
            end if
         end if

         call e00(s, xscale, i_dlnE_dt, i_u, k, nvar, scal*d_dedt_du00*f)
         call e00(s, xscale, i_dlnE_dt, i_lnR, k, nvar, scal*d_dedt_dlnR00*f)
         if (i_lum /= 0) &
            call e00(s, xscale, i_dlnE_dt, i_lum, k, nvar, scal*d_dedt_dL00*f)
         call e00(s, xscale, i_dlnE_dt, i_lnd, k, nvar, scal*d_dedt_dlnd00*f)
         call e00(s, xscale, i_dlnE_dt, i_lnT, k, nvar, scal*d_dedt_dlnT00*f)
         if (s% w_div_wc_flag) &
            call e00(s, xscale, i_dlnE_dt, i_w_div_wc, k, nvar, scal*d_dedt_dw00*f)

         if (k > 1) then ! do k-1 partials   
            call em1(s, xscale, i_dlnE_dt, i_u, k, nvar, scal*d_dedt_dum1*f)            
            call em1(s, xscale, i_dlnE_dt, i_lnd, k, nvar, scal*d_dedt_dlndm1*f)              
            call em1(s, xscale, i_dlnE_dt, i_lnT, k, nvar, scal*d_dedt_dlnTm1*f) 
         end if

         if (k < nz) then ! do k+1 partials
            call ep1(s, xscale, i_dlnE_dt, i_u, k, nvar, scal*d_dedt_dup1*f)           
            call ep1(s, xscale, i_dlnE_dt, i_lnR, k, nvar, scal*d_dedt_dlnRp1*f)
            if (i_lum /= 0) &
               call ep1(s, xscale, i_dlnE_dt, i_lum, k, nvar, scal*d_dedt_dLp1*f)           
            call ep1(s, xscale, i_dlnE_dt, i_lnd, k, nvar, scal*d_dedt_dlndp1*f)              
            call ep1(s, xscale, i_dlnE_dt, i_lnT, k, nvar, scal*d_dedt_dlnTp1*f) 
            if (s% w_div_wc_flag) &
               call ep1(s, xscale, i_dlnE_dt, i_w_div_wc, k, nvar, scal*d_dedt_dwp1*f)
         end if



        ! d_dedt_dL00 
        ! d_dedt_dlnR00 
        ! d_dedt_du00 
        ! d_dedt_dlnd00 
        ! d_dedt_dlnT00 
         
        ! d_dedt_dum1 
        ! d_dedt_dlndm1 
        ! d_dedt_dlnTm1 

        ! d_dedt_dLp1 
        ! d_dedt_dlnRp1 
        ! d_dedt_dup1 
        ! d_dedt_dlndp1 
        ! d_dedt_dlnTp1 

         if (test_partials) then   
            s% solver_test_partials_var = i_lnT
            s% solver_test_partials_dval_dx = d_dedt_dlnT00
            write(*,2) 'logRho', k, s% lnd(k)/ln10
            write(*,2) 'logT', k, s% lnT(k)/ln10
            write(*,2) 'X', k, s% X(k)
            write(*,2) 'Z', k, s% Z(k)
         end if

      end subroutine do1_Riemann_energy_eqn
         
         
      subroutine do1_Riemann_dlnRdt_eqn( &
            s, k, xscale, equ, skip_partials, nvar, ierr)
         type (star_info), pointer :: s         
         integer, intent(in) :: k, nvar
         real(dp), pointer :: xscale(:,:)
         real(dp), pointer :: equ(:,:)
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr
         
         integer :: nz, i_dlnR_dt, i_u, i_lnR, i_w_div_wc
         real(dp) :: r, r0, r_div_r0, u_expected, u_factor, c_factor, &
            dr_div_r0_actual, dr_div_r0_expected, dt, alpha, uc_factor

         logical :: test_partials
         
         include 'formats'

         !test_partials = (k == s% solver_test_partials_k)
         test_partials = .false.

         ierr = 0
         dt = s% dt
         nz = s% nz
         i_dlnR_dt = s% i_dlnR_dt
         i_u = s% i_u
         i_lnR = s% i_lnR
         i_w_div_wc = s% i_w_div_wc
         
         r = s% r(k)
         r0 = s% r_start(k)
         r_div_r0 = r/r0

         u_expected = s% u_face(k)
         u_factor = 1d0

         ! dr = r - r0 = u_expected*dt
         ! eqn: dr/r0 = u_expected*dt/r0
         ! write dr/r0 using dlnR = dlnR_dt*dt 
         ! dr/r0 = (r - r0)/r0 = r/r0 - 1
         !        = exp(lnR)/exp(lnR_start) - 1
         !        = exp(lnR - lnR_start) - 1
         !        = exp(dlnR) - 1 = exp(dlnR_dt*dt) - 1
         ! eqn becomes: u_expected*dt/r0 = expm1(dlnR)
         dr_div_r0_actual = expm1(s% dlnR_dt(k)*dt) ! expm1(x) = E^x - 1

         c_factor = dt/r0 ! factor to convert dr_dt to dr_div_r0
         dr_div_r0_expected = u_expected*c_factor

         equ(i_dlnR_dt, k) = dr_div_r0_expected - dr_div_r0_actual         
         s% lnR_residual(k) = equ(i_dlnR_dt, k)
         
         if (is_bad(equ(i_dlnR_dt, k))) then
            ierr = -1
            if (.not. s% report_ierr) return
!$omp critical (do1_Riemann_dlnRdt_eqn_omp)
            write(*,2) 'equ(i_dlnR_dt, k)', k, equ(i_dlnR_dt, k)
            write(*,2) 'dr_div_r0_expected', k, dr_div_r0_expected
            write(*,2) 'dr_div_r0_actual', k, dr_div_r0_actual
            write(*,2) 'u_expected', k, u_expected
            write(*,2) 's% u_face(k)', k, s% u_face(k)
            write(*,2) 'dt/r0', k, dt/r0
            stop 'do1_Riemann_dlnRdt_eqn'
!$omp end critical (do1_Riemann_dlnRdt_eqn_omp)
         end if

         if (test_partials) then
            s% solver_test_partials_val = equ(i_dlnR_dt, k) 
         end if
     
         if (skip_partials) return
         
         ! partial of -dr_div_r0_actual
         call e00(s, xscale, i_dlnR_dt, i_lnR, k, nvar, -r_div_r0)
         
         ! partials of dr_div_r0_expected = u_face*f
         uc_factor = u_factor*c_factor
         call e00(s, xscale, i_dlnR_dt, i_lnR, k, nvar, uc_factor*s% d_uface_dlnR(k))
         call e00(s, xscale, i_dlnR_dt, i_u, k, nvar, uc_factor*s% d_uface_du00(k))
         
         call e00(s, xscale, i_dlnR_dt, s% i_lnd, k, nvar, uc_factor*s% d_uface_dlnd00(k))
         if (s% do_struct_thermo) &
            call e00(s, xscale, i_dlnR_dt, s% i_lnT, k, nvar, uc_factor*s% d_uface_dlnT00(k))
         
         if (k > 1) then
            call em1(s, xscale, i_dlnR_dt, i_u, k, nvar, uc_factor*s% d_uface_dum1(k))            
            call em1(s, xscale, i_dlnR_dt, s% i_lnd, k, nvar, uc_factor*s% d_uface_dlndm1(k))
            if (s% do_struct_thermo) &
               call em1(s, xscale, i_dlnR_dt, s% i_lnT, k, nvar, uc_factor*s% d_uface_dlnTm1(k))
         end if

         if (s% w_div_wc_flag) then
            call e00(s, xscale, i_dlnR_dt, i_w_div_wc, k, nvar, uc_factor*s% d_uface_dw(k))
         end if

         if (test_partials) then   
            s% solver_test_partials_var = i_w_div_wc
            s% solver_test_partials_dval_dx = uc_factor*s% d_uface_dw(k)
         end if
            
      end subroutine do1_Riemann_dlnRdt_eqn

         
      end module hydro_riemann

