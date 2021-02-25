! ***********************************************************************
!
!   Copyright (C) 2018-2019  Bill Paxton & The MESA Team
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful, 
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module rsp_eval_eos_and_kap
      use eos_def
      use eos_lib
      use chem_def
      use chem_lib, only: chem_Xsol, basic_composition_info
      use kap_lib
      use kap_def
      use const_lib
      use utils_lib
      use star_private_def
      use rsp_def, only: xa, X, Z, Y, &
         abar, zbar, z53bar, XC, XN, XO, Xne


      implicit none

      integer :: species
      integer, pointer, dimension(:) :: net_iso, chem_id
      integer :: eos_handle, kap_handle
      
      contains  
     
      subroutine mesa_eos_kap (s,k,G,H, &   !input: temp,volume  &
               P,PV,PT,E,EV,ET,CP,CPV,dCp_dT_00, &
               Q,QV,QT,OP,OPV,OPT,ierr) 
      implicit none
      type (star_info), pointer :: s
      integer, intent(out) :: ierr
      integer :: k, j
      real(8) :: G,H,P,PV,PT, &
         E,EV,ET,CP,CPV,dCp_dT_00, &
         Q,QV,QT,OP,OPV,OPT,cs, &
         Pgas,d_Pg_dV,d_Pg_dT,Prad,d_Pr_dT, &
         egas,d_egas_dV,d_egas_dT,erad,d_erad_dV,d_erad_dT
      include 'formats'
      if (k <= 0 .or. k > s% nz) then
         j = 0
      else
         j = s% nz + 1 - k
      end if
      if (is_bad(G+H)) then
         write(*,2) 'LINA mesa_eos_kap G H', k, G, H
         stop 'mesa_eos_kap'
      end if
      call eval1_mesa_eos_and_kap(s,j,.false.,G,H, &
               Pgas,d_Pg_dV,d_Pg_dT,Prad,d_Pr_dT,&
               egas,d_egas_dV,d_egas_dT,erad,d_erad_dV,d_erad_dT, &
               cs,CP,CPV,dCp_dT_00, &
               Q,QV,QT,OP,OPV,OPT,ierr) 
      if (ierr /= 0) return
      E = egas + erad
      EV = d_egas_dV + d_erad_dV
      ET = d_egas_dT + d_erad_dT
      P = Pgas + Prad
      PV = d_Pg_dV
      PT = d_Pg_dT + d_Pr_dT
      end subroutine mesa_eos_kap
      
      
      subroutine restart_rsp_eos_and_kap(s)
         type (star_info), pointer :: s
         integer :: ierr
         eos_handle = s% eos_handle
         kap_handle = s% kap_handle
         net_iso => s% net_iso
         chem_id => s% chem_id
         species = s% species
      end subroutine restart_rsp_eos_and_kap
      
      subroutine init_for_rsp_eos_and_kap(s)
         use adjust_xyz, only: get_xa_for_standard_metals
         type (star_info), pointer :: s
         
         integer :: i, k, j, iz, ierr
         real(dp) :: initial_z, initial_y, initial_x, &
            initial_h1, initial_h2, initial_he3, initial_he4, &
            xsol_he3, xsol_he4, z2bar, ye, mass_correction, sumx
         include 'formats'
         ierr = 0         
         eos_handle = s% eos_handle
         kap_handle = s% kap_handle
         net_iso => s% net_iso
         chem_id => s% chem_id
         species = s% species
         
         initial_x = max(0d0, min(1d0, s% RSP_X))
         initial_z = max(0d0, min(1d0, s% RSP_Z))
         initial_y = max(0d0,1d0 - (initial_x + initial_z))
         initial_h1 = initial_x
         initial_h2 = 0
         if (s% initial_he3 < 0d0) then
            xsol_he3 = chem_Xsol('he3')
            xsol_he4 = chem_Xsol('he4')
            initial_he3 = initial_y*xsol_he3/(xsol_he3 + xsol_he4)
            initial_he4 = initial_y*xsol_he4/(xsol_he3 + xsol_he4)
         else if (s% initial_he3 < initial_y) then
            initial_he3 = s% initial_he3
            initial_he4 = initial_y - s% initial_he3
         else
            write(*,*) 'ERROR: initial_he3 is larger than initial_y'
            ierr = -1
            return
         end if
         call get_xa_for_standard_metals(s, &
            species, chem_id, net_iso, &
            initial_h1, initial_h2, initial_he3, initial_he4, &
            s% job% initial_zfracs, s% initial_dump_missing_heaviest, &
            xa, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in get_xa_for_standard_metals'
            return
         end if
         if (abs(1-sum(xa(:))) > 1d-8) then
            write(*,1) 'initial_h1', initial_h1
            write(*,1) 'initial_h2', initial_h2
            write(*,1) 'initial_he3', initial_he3
            write(*,1) 'initial_he4', initial_he4
            write(*,1) 'initial_y', initial_y
            write(*,1) 'initial_z', initial_z
            write(*,1) 'initial_h1+h2+he3+he4+z', &
               initial_h1 + initial_h2 + initial_he3 + initial_he4 + initial_z
            write(*,1) 'sum(xa(:))', sum(xa(:))
            write(*,*) 'init_for_rsp_eos_and_kap'
            ierr = -1
            return
         end if
         call basic_composition_info( &
            species, s% chem_id, xa(:), X, Y, Z, abar, zbar, z2bar, z53bar, ye, &
            mass_correction, sumx)
         if (is_bad(zbar)) then
            write(*,1) 'basic_composition_info initial_x', initial_x
            write(*,1) 'basic_composition_info initial_z', initial_z
            write(*,1) 'basic_composition_info initial_y', initial_y
            write(*,1) 'basic_composition_info zbar', zbar
            do i=1,species
               write(*,2) 'xa', i, xa(i)
            end do
            stop
         end if
         xc = 0; xn = 0; xo = 0; xne = 0
         do i=1, species
            iz = chem_isos% Z(chem_id(i))
            select case(iz)
               case (6)
                  xc = xc + xa(i)
               case (7)
                  xn = xn + xa(i)
               case (8)
                  xo = xo + xa(i)
               case (10)
                  xne = xne + xa(i)
            end select
         end do

         return

         write(*,1) 'init_for_rsp_eos_and_kap X', X
         write(*,1) 'Y', Y
         write(*,1) 'Z', Z       
         write(*,1) 'abar', abar
         write(*,1) 'zbar', zbar
         write(*,1) 'XC', XC
         write(*,1) 'XN', XN
         write(*,1) 'XO', XO
         write(*,1) 'Xne', Xne
         do i=1,species
            write(*,2) 'xa', i, xa(i)
         end do
         write(*,*) trim(s% net_name)
         stop 'init_for_rsp_eos_and_kap'
      end subroutine init_for_rsp_eos_and_kap
      
      
      subroutine eval1_mesa_eos_and_kap( &
            s,k,skip_kap,T_in,V, &
            Pgas,d_Pg_dV,d_Pg_dT,Prad,d_Pr_dT, &
            egas,d_egas_dV,d_egas_dT,erad,d_erad_dV,d_erad_dT, &
            CSND,CP,CPV,CPT,Q,QV,QT,OP,OPV,OPT,ierr)
         use star_utils, only: write_eos_call_info
         use eos_support, only: get_eos
         use micro, only: store_eos_for_cell
         type (star_info), pointer :: s
         integer, intent(in) :: k
         logical, intent(in) :: skip_kap
         real(dp), intent(in) :: T_in,V
         real(dp), intent(out) :: &
            Pgas,d_Pg_dV,d_Pg_dT,Prad,d_Pr_dT, &
            egas,d_egas_dV,d_egas_dT,erad,d_erad_dV,d_erad_dT, &
            CSND,CP,CPV,CPT,Q,QV,QT,OP,OPV,OPT
         integer, intent(out) :: ierr
         
         integer :: j
         real(dp) :: logT, logRho, T, Rho, energy, entropy, &
            dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
            d_dlnRho_const_T, d_dlnT_const_Rho, E, dE_dV, dE_dT, &
            lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            frac_Type2, opacity, dlnkap_dlnd, dlnkap_dlnT, &
            dlnd_dV, chiRho, chiT, dchiRho_dlnd, dchiRho_dlnT, &
            dchiT_dlnd, dchiT_dlnT, dQ_dlnd, dQ_dlnT, opacity_factor
         real(dp), dimension(num_eos_basic_results) :: &
            res, d_dlnd, d_dlnT, d_dabar, d_dzbar
         real(dp) :: d_dxa(num_eos_d_dxa_results,s% species)

         include 'formats'
         
         rho = 1d0/V
         T = T_in
         logRho = log10(rho)
         dlnd_dV = -rho
         logT = log10(T)  
         
         if (k > 0 .and. k <= s% nz) then
            s% rho(k) = rho
            s% lnd(k) = logRho*ln10
            s% xh(s% i_lnd,k) = s% lnd(k)
            s% T(k) = T
            s% lnT(k) = logT*ln10
            s% xh(s% i_lnT,k) = s% lnT(k)
            s% abar(k) = abar
            s% zbar(k) = zbar
            s% z53bar(k) = z53bar
            call get_eos( &
               s, k, xa, &
               Rho, logRho, T, logT, &
               res, d_dlnd, d_dlnT, &
               d_dxa, ierr)
            if (ierr == 0) then
               s% lnPgas(k) = res(i_lnPgas)
               s% Pgas(k) = exp(s% lnPgas(k))
               do j=1,num_eos_basic_results
                  s% d_eos_dlnd(j,k) = d_dlnd(j)
                  s% d_eos_dlnT(j,k) = d_dlnT(j)
               end do
               call store_eos_for_cell(s, k, res, d_dlnd, d_dlnT, d_dxa, ierr)
               CSND = s% csound(k)
            end if
         else ! k <= 0 or k > nz
            call get_eos( &
               s, 0, xa, &
               Rho, logRho, T, logT, &
               res, d_dlnd, d_dlnT, &
               d_dxa, ierr)
            if (ierr == 0) then
               Prad = crad*T**4/3d0
               Pgas = exp(res(i_lnPgas))
               CSND = sqrt(max(0d0, res(i_gamma1)*(Prad+Pgas)/Rho)) 
            end if
         end if
                   
         if (ierr /= 0) then
            !$OMP critical
            if (k > 0 .and. k < s% nz) call write_eos_call_info(s,k)
            write(*,2) 'X', k, X
            write(*,2) 'Z', k, Z
            write(*,2) 'zbar', k, zbar
            write(*,2) 'abar', k, abar
            write(*,2) 'V', k, V
            write(*,2) 'rho', k, rho
            write(*,2) 'T', k, T
            write(*,2) 'logRho', k, logRho
            write(*,2) 'logT', k, logT
            if (s% stop_for_bad_nums .and. is_bad(logRho+logT)) stop 'do_eos_for_cell'
            !$OMP end critical
            !return
            stop 'RSP failed in get_eos'
         end if
         
         if (skip_kap) then
            OP=0; OPV=0; OPT=0
         else
            call eval1_kap(s, k, skip_kap, &
               V, logRho, dlnd_dV, T, logT, species, chem_id, net_iso, xa, &
               res, d_dlnd, d_dlnT, OP, OPV, OPT, ierr)
         end if
         
         lnfree_e = res(i_lnfree_e)
         d_lnfree_e_dlnRho = d_dlnd(i_lnfree_e)
         d_lnfree_e_dlnT = d_dlnT(i_lnfree_e)
         
         Prad = crad*T**4/3d0  ! erg/cm^2
         d_Pr_dT = 4d0*Prad/T

         erad = 3d0*Prad/rho ! 3*Prad*V   erg/gm
         d_erad_dT = 3d0*d_Pr_dT/rho
         d_erad_dV = 3d0*Prad
         
         E = exp(res(i_lnE))
         dE_dV = E*d_dlnd(i_lnE)*dlnd_dV
         dE_dT = E*d_dlnT(i_lnE)/T
         egas = E - erad
         d_egas_dV = dE_dV - d_erad_dV
         d_egas_dT = dE_dT - d_erad_dT
         
         Pgas = exp(res(i_lnPgas))
         d_Pg_dV = Pgas*d_dlnd(i_lnPgas)*dlnd_dV
         d_Pg_dT = Pgas*d_dlnT(i_lnPgas)/T
         
         CP = res(i_Cp)
         CPV = d_dlnd(i_Cp)*dlnd_dV
         CPT = d_dlnT(i_Cp)/T
         
         chiT = res(i_chiT)
         dchiT_dlnd = d_dlnd(i_chiT)
         dchiT_dlnT = d_dlnT(i_chiT)
         
         chiRho = res(i_chiRho)
         dchiRho_dlnd = d_dlnd(i_chiRho)
         dchiRho_dlnT = d_dlnT(i_chiRho)
         
         Q = chiT/(rho*T*chiRho) ! thermal expansion coefficient
         dQ_dlnd = Q*(dchiT_dlnd/chiT - dchiRho_dlnd/chiRho - 1d0)
         dQ_dlnT = Q*(dchiT_dlnT/chiT - dchiRho_dlnT/chiRho - 1d0)
         QV = dQ_dlnd*dlnd_dV
         QT = dQ_dlnT/T
         
         if (is_bad(egas) .or. egas <= 0d0) then
            ierr = -1
            return
         end if
          
      end subroutine eval1_mesa_eos_and_kap
      
      
      subroutine eval1_kap(s, k, skip_kap, &
            V, logRho, dlnd_dV, T, logT, species, chem_id, net_iso, xa, &
            res, d_dlnd, d_dlnT, OP, OPV, OPT, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         logical, intent(in) :: skip_kap
         integer, intent(in) :: species
         integer, pointer :: chem_id(:), net_iso(:)
         real(dp), intent(in) :: xa(:) 
         real(dp), intent(in) :: V, logRho, T, logT, dlnd_dV
         real(dp), intent(in), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT
         real(dp), intent(out) :: OP, OPV, OPT
         integer, intent(out) :: ierr
         
         real(dp) :: lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
            eta, d_eta_dlnRho, d_eta_dlnT, &
            frac_Type2, opacity, dlnkap_dlnd, dlnkap_dlnT, opacity_factor
         real(dp) :: kap_fracs(num_kap_fracs), dlnkap_dxa(s% species)
         
         include 'formats'
         ierr = 0
         
         if (s% RSP_kap_density_factor > 0d0) then
            OP = s% RSP_kap_density_factor*V
            OPV = s% RSP_kap_density_factor
            OPT = 0d0
            return
         end if
      
         lnfree_e = res(i_lnfree_e)
         d_lnfree_e_dlnRho = d_dlnd(i_lnfree_e)
         d_lnfree_e_dlnT = d_dlnT(i_lnfree_e)

         eta = res(i_eta)
         d_eta_dlnRho = d_dlnd(i_eta)
         d_eta_dlnT = d_dlnT(i_eta)

         if (s% use_other_kap) then
            call s% other_kap_get( &
               s% id, k, kap_handle, species, chem_id, net_iso, xa, &
               logRho, logT, &
               lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
               eta, d_eta_dlnRho, d_eta_dlnT, &
               kap_fracs, opacity, dlnkap_dlnd, dlnkap_dlnT, dlnkap_dxa, ierr)
         else
            call kap_get( &
               kap_handle, species, chem_id, net_iso, xa, &
               logRho, logT, &
               lnfree_e, d_lnfree_e_dlnRho, d_lnfree_e_dlnT, &
               eta, d_eta_dlnRho, d_eta_dlnT, &
               kap_fracs, opacity, dlnkap_dlnd, dlnkap_dlnT, dlnkap_dxa, ierr)
         end if
         frac_Type2 = kap_fracs(i_frac_Type2)
         if (ierr /= 0) then
!$OMP critical
            write(*,*) 'failed in eval1_mesa_eos_and_kap get kap'
            write(*,2) 'logRho', k, logRho
            write(*,2) 'logT', k, logT
            write(*,2) 'lnfree_e', k, lnfree_e
            write(*,2) 'zbar', k, zbar
            write(*,2) 'X', k, X
            write(*,2) 'Z', k, Z
            stop 'eval1_kap'
!$OMP end critical
            !return
            stop 1
         end if

         if (k > 0 .and. k <= s% nz) then
            opacity_factor = s% extra_opacity_factor(k)
         else
            opacity_factor = s% opacity_factor
         end if
         if (opacity_factor /= 1d0) then
            if (s% min_logT_for_opacity_factor_off > 0) then
               if (logT >= s% max_logT_for_opacity_factor_off .or. &
                   logT <= s% min_logT_for_opacity_factor_off) then
                  opacity_factor = 1
               else if (logT > s% max_logT_for_opacity_factor_on) then
                  opacity_factor = 1 + (opacity_factor-1)* &
                     (logT - s% max_logT_for_opacity_factor_off)/ &
                     (s% max_logT_for_opacity_factor_on - s% max_logT_for_opacity_factor_off)
               else if (logT < s% min_logT_for_opacity_factor_on) then
                  opacity_factor = 1 + (opacity_factor-1)* &
                     (logT - s% min_logT_for_opacity_factor_off)/ &
                     (s% min_logT_for_opacity_factor_on - s% min_logT_for_opacity_factor_off)
               end if
            end if
         end if

         OP = opacity*opacity_factor
         OPV = dlnkap_dlnd*opacity*dlnd_dV
         OPT = dlnkap_dlnT*opacity/T
         if (k > 0 .and. k <= s% nz) then
            s% opacity(k) = opacity
         end if
      
      end subroutine eval1_kap
      
            
      subroutine eval1_mesa_Rho_given_PT(s, k, P, T, rho_guess, rho, ierr)
         use eos_support, only: solve_eos_given_PT
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: P, T, rho_guess
         real(dp), intent(out) :: rho
         integer, intent(out) :: ierr           
                
         real(dp) :: logT, logP, logRho_guess, logRho, &
              logRho_tol, logP_tol
         real(dp), dimension(num_eos_basic_results) :: &
            res, d_dlnd, d_dlnT, d_dabar, d_dzbar
         integer :: iter
            
         include 'formats'

         ierr = 0

         logT = log10(T)
         logP = log10(P)
         logRho_guess = log10(rho_guess)

         logRho_tol = 5d-10
         logP_tol = 5d-13

         ! in some parts of (logRho,logT) cannot meet tight tolerances.
         ! so be prepared to relax them.
         do iter = 1, 9
            call solve_eos_given_PT( &
                 s, k, Z, X, abar, zbar, xa, &
                 logT, logP, logRho_guess, logRho_tol, logP_tol, &
                 logRho, res, d_dlnd, d_dlnT, d_dabar, d_dzbar, &
                 ierr)            
            if (ierr == 0) exit
            logRho_tol = logRho_tol*3d0
            logP_tol = logP_tol*3d0
         end do

         rho = exp10(logRho)

         if (ierr /= 0) then
            return
            write(*,2) 'Z', k, Z
            write(*,2) 'X', k, X
            write(*,2) 'abar', k, abar
            write(*,2) 'zbar', k, zbar
            write(*,2) 'logT', k, logT
            write(*,2) 'logP', k, logP
            write(*,2) 'logRho_guess', k, logRho_guess
         end if
               
      end subroutine eval1_mesa_Rho_given_PT
      
      
      real(dp) function eval1_gamma_PT_getRho(s, k, P, T, ierr)
         use eos_lib, only: eos_gamma_PT_get
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: P, T
         integer, intent(out) :: ierr                  
         real(dp), dimension(num_eos_basic_results) :: &
            res, d_dlnd, d_dlnT
         real(dp) :: logP, logT, logRho, rho, gamma
         include 'formats'
         logP = log10(P)
         logT = log10(T)         
         gamma = 5d0/3d0
         call eos_gamma_PT_get( &
            eos_handle, abar, P, logP, T, logT, gamma, &
            rho, logRho, res, d_dlnd, d_dlnT, &
            ierr)
         eval1_gamma_PT_getRho = rho
      end function eval1_gamma_PT_getRho


      subroutine get_surf_P_T_kap(s, &
            M, R, L, tau, kap_guess, &
            T, P, kap, Teff, ierr)

         use atm_support, only: get_atm_PT_legacy_grey_and_kap

         type (star_info), pointer :: s
         real(dp), intent(in) :: M, R, L, tau, kap_guess
         real(dp), intent(out) :: T, P, kap, Teff
         integer, intent(out) :: ierr                  
         
         real(dp) :: &
            lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
         integer :: iter
         
         ierr = 0

         call get_atm_PT_legacy_grey_and_kap( &
              s, tau, L, R, M, s% cgrav(1), .TRUE., & 
              Teff, kap, &
              lnT, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
              lnP, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, &
              ierr)

         T = exp(lnT)
         P = exp(lnP)
         
      end subroutine get_surf_P_T_kap

      end module rsp_eval_eos_and_kap
      

