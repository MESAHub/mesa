! ***********************************************************************
!
!   Copyright (C) 2009  Bill Paxton
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

      module eval_eosPT
      
      use eos_def
      use eos_lib
      use const_def
      use utils_lib
      use chem_def
      
      implicit none

      integer :: eos_handle

      integer, parameter :: species = 7
      integer, parameter :: h1=1, he4=2, c12=3, n14=4, o16=5, ne20=6, mg24=7
      integer, target :: chem_id_array(species)
      integer, pointer, dimension(:) :: chem_id, net_iso
      double precision :: xa(species)
      
      
      contains


      subroutine Init_Composition(X_in, Zinit_in, abar, zbar)
         use chem_lib
         double precision, intent(IN) :: X_in, Zinit_in
         double precision, intent(out) :: abar, zbar
         
         ! note that these don't matter for the tables we produce
         ! since OPAL has its own metals fractions, and SCVH doesn't have any.
         double precision, parameter :: Zfrac_C = 0.173312d0
         double precision, parameter :: Zfrac_N = 0.053177d0
         double precision, parameter :: Zfrac_O = 0.482398d0
         double precision, parameter :: Zfrac_Ne = 0.098675d0
         
         double precision :: frac, sumx, mass_correction
         double precision :: X, Y, Z, z2bar, z53bar, ye, Zinit, XC, XO
         
         chem_id => chem_id_array
         
         allocate(net_iso(num_chem_isos))         
         net_iso(:) = 0
         
         chem_id(h1) = ih1; net_iso(ih1) = h1
         chem_id(he4) = ihe4; net_iso(ihe4) = he4
         chem_id(c12) = ic12; net_iso(ic12) = c12
         chem_id(n14) = in14; net_iso(in14) = n14
         chem_id(o16) = io16; net_iso(io16) = o16
         chem_id(ne20) = ine20; net_iso(ine20) = ne20
         chem_id(mg24) = img24; net_iso(img24) = mg24
         
         X = X_in
         Zinit = Zinit_in
         XC = 0; XO = 0
         Y = 1 - (X + Zinit + XC + XO)
         if (Y < 0) then ! adjust XC and XO
            if (XC + XO <= 0) then
               write(*,*) 'bad args to Init_Composition'
               call mesa_error(__FILE__,__LINE__)
            end if
            frac = (1 - X - Zinit) / (XC + XO)
            if (frac <= 0) stop 'bad args to Init_Composition'
            XC = frac*XC; XO = frac*XO
            Y = 1 - (X+Zinit+XC+XO)
            if (Y < -1d-10) then
               write(*,*) 'screw up in Init_Composition'
               call mesa_error(__FILE__,__LINE__)
            end if
            if (Y < 0) Y = 0
         end if
      
         xa(h1) = X
         xa(he4) = Y
         xa(c12) = Zinit * Zfrac_C + XC
         xa(n14) = Zinit * Zfrac_N
         xa(o16) = Zinit * Zfrac_O + XO
         xa(ne20) = Zinit * Zfrac_Ne
         xa(species) = 1 - sum(xa(1:species-1))
         
         call basic_composition_info( &
               species, chem_id, xa, X, Y, Z, abar, zbar, z2bar, z53bar, ye, &
               mass_correction, sumx)
         
      end subroutine Init_Composition

      
      subroutine eos_4_Pgas_T( &
            logW, logT, T, logPgas, Pgas, &
            abar, zbar, X, Z, logRho_guess, &
            logRho, logE, logS, chiRho, chiT, &
            Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
            mu, log_free_e, gamma1, gamma3, grad_ad, eta, &
            expanded_scvh, call_number, ierr)
         double precision, intent(in) :: &
            logW, logT, T, logPgas, Pgas, &
            abar, zbar, X, Z, logRho_guess
         double precision, intent(out) :: &
            logRho, logE, logS, chiRho, chiT, &
            Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
            mu, log_free_e, gamma1, gamma3, grad_ad, eta
         logical, intent(in) :: expanded_scvh
         integer, intent(in) :: call_number
         integer, intent(out) :: ierr
         
         double precision :: &
            alfa, beta, logRho_scvh, logE_scvh, logS_scvh, chiRho_scvh, chiT_scvh, &
            Cp_scvh, Cv_scvh, dE_dRho_scvh, dS_dT_scvh, dS_dRho_scvh, &
            mu_scvh, log_free_e_scvh, gamma1_scvh, gamma3_scvh, grad_ad_scvh, eta_scvh, &
            logRho_mesa, logE_mesa, logS_mesa, chiRho_mesa, chiT_mesa, &
            Cp_mesa, Cv_mesa, dE_dRho_mesa, dS_dT_mesa, dS_dRho_mesa, &
            mu_mesa, log_free_e_mesa, gamma1_mesa, gamma3_mesa, grad_ad_mesa, eta_mesa
            
         include 'formats.dek'
         
         ierr = 0
         
         if (expanded_scvh) then
            alfa = eval_alfa_expanded_scvh()
         else
            alfa = eval_alfa()
         end if
         
         beta = 1d0 - alfa
         
         if (alfa > 0d0) then
            call scvh_4_PT( &
               logW, logT, T, logPgas, Pgas, &
               abar, zbar, X, Z, logRho_guess, &
               logRho_scvh, logE_scvh, logS_scvh, chiRho_scvh, chiT_scvh, &
               Cp_scvh, Cv_scvh, dE_dRho_scvh, dS_dT_scvh, dS_dRho_scvh, &
               mu_scvh, log_free_e_scvh, gamma1_scvh, gamma3_scvh, grad_ad_scvh, eta_scvh, &
               call_number, ierr)
            if (ierr /= 0 .or. is_bad(logE_scvh) .or.  &
                  is_bad(logRho_scvh) .or. is_bad(logS_scvh) .or. grad_ad_scvh <= 0) then
               ierr = 0
               alfa = 0d0
               beta = 1d0
            end if
         end if
         
         if (alfa == 0d0) then
            logRho_scvh = 0d0
            logE_scvh = 0d0
            logS_scvh = 0d0
            chiRho_scvh = 0d0
            chiT_scvh = 0d0
            Cp_scvh = 0d0
            Cv_scvh = 0d0
            dE_dRho_scvh = 0d0
            dS_dT_scvh = 0d0
            dS_dRho_scvh = 0d0
            mu_scvh = 0d0
            log_free_e_scvh = 0d0
            gamma1_scvh = 0d0
            gamma3_scvh = 0d0
            grad_ad_scvh = 0d0
            eta_scvh = 0d0
         end if
         
         if (beta > 0d0) then
            call mesa_4_PT( &
               logW, logT, T, logPgas, Pgas, &
               abar, zbar, X, Z, logRho_guess, &
               logRho_mesa, logE_mesa, logS_mesa, chiRho_mesa, chiT_mesa, &
               Cp_mesa, Cv_mesa, dE_dRho_mesa, dS_dT_mesa, dS_dRho_mesa, &
               mu_mesa, log_free_e_mesa, gamma1_mesa, gamma3_mesa, grad_ad_mesa, eta_mesa, &
               call_number, ierr)
         else
            logRho_mesa = 0d0
            logE_mesa = 0d0
            logS_mesa = 0d0
            chiRho_mesa = 0d0
            chiT_mesa = 0d0
            Cp_mesa = 0d0
            Cv_mesa = 0d0
            dE_dRho_mesa = 0d0
            dS_dT_mesa = 0d0
            dS_dRho_mesa = 0d0
            mu_mesa = 0d0
            log_free_e_mesa = 0d0
            gamma1_mesa = 0d0
            gamma3_mesa = 0d0
            grad_ad_mesa = 0d0
            eta_mesa = 0d0
         end if
         
         logRho = alfa*logRho_scvh + beta*logRho_mesa
         logE = alfa*logE_scvh + beta*logE_mesa
         logS = alfa*logS_scvh + beta*logS_mesa
         chiRho = alfa*chiRho_scvh + beta*chiRho_mesa
         chiT = alfa*chiT_scvh + beta*chiT_mesa
         Cp = alfa*Cp_scvh + beta*Cp_mesa
         Cv = alfa*Cv_scvh + beta*Cv_mesa
         dE_dRho = alfa*dE_dRho_scvh + beta*dE_dRho_mesa
         dS_dT = alfa*dS_dT_scvh + beta*dS_dT_mesa
         dS_dRho = alfa*dS_dRho_scvh + beta*dS_dRho_mesa
         mu = alfa*mu_scvh + beta*mu_mesa
         log_free_e = alfa*log_free_e_scvh + beta*log_free_e_mesa
         gamma1 = alfa*gamma1_scvh + beta*gamma1_mesa
         gamma3 = alfa*gamma3_scvh + beta*gamma3_mesa
         grad_ad = alfa*grad_ad_scvh + beta*grad_ad_mesa
         eta = alfa*eta_scvh + beta*eta_mesa

         if (is_bad(Cv)) then
            write(*,1) 'alfa', alfa
            write(*,1) 'beta', beta
            write(*,1) 'logW', logW
            write(*,1) 'logT', logT
            write(*,1) 'logPgas', logPgas
            write(*,1) 'logRho', logRho
            write(*,1) 'Cv_scvh', Cv_scvh
            write(*,1) 'Cv_mesa', Cv_mesa
            write(*,1) 'Cv', Cv
            write(*,1) 'X', X
            write(*,1) 'Z', Z
            write(*,*)
            stop 'eos_4_Pgas_T'
         end if
         
         contains
         
         integer function eval_alfa_expanded_scvh() result(alfa)
            double precision :: &
               dx, dy, logT_max, logT_max_minus, logT_min, logT_min_plus, &
               logPgas_max, logPgas_max_minus, logPgas_min, logPgas_min_plus, &
               logW_min, logW_min_plus, dlogT, dlogPgas, dlogW
            
            dlogT = 0.25d0
            logT_max = 7.0d0
            logT_max_minus = logT_max - dlogT
            logT_min = 2.1d0
            logT_min_plus = logT_min + dlogT
            
            dlogPgas = 1d0
            logPgas_max = 19d0
            logPgas_max_minus = logPgas_max - dlogPgas
            logPgas_min = 0.0d0
            logPgas_min_plus = logPgas_min + dlogPgas
            
            dlogW = 0.5d0
            logW_min = -12d0
            logW_min_plus = logW_min + dlogW

            ! alfa = fraction from scvh; beta = fraction from mesa
            if (logW <= logW_min) then
               alfa = 0 ! pure MESA
            else if (logW < logW_min_plus) then
               alfa = (logW - logW_min)/dlogW
            else if (logT >= logT_max) then
               alfa = 0 ! pure MESA
            else if (logT >= logT_max_minus) then
               if (logPgas >= logPgas_max) then
                  alfa = 0 ! pure MESA
               else if (logPgas >= logPgas_max_minus) then
                  dx = (logPgas - logPgas_max_minus)/dlogPgas
                  dy = (logT - logT_max_minus)/dlogT
                  alfa = 1d0 - min(1d0, sqrt(dx**2 + dy**2))
               else if (logPgas >= logPgas_min_plus) then
                  alfa = 1d0 - (logT - logT_max_minus)/dlogT              
               else if (logPgas >= logPgas_min) then
                  dx = (logPgas_min_plus - logPgas)/dlogPgas
                  dy = (logT - logT_max_minus)/dlogT
                  alfa = 1d0 - min(1d0, sqrt(dx**2 + dy**2))
               else ! logPgas < logPgas_min
                  alfa = 0 ! pure MESA
               end if
            else if (logT >= logT_min_plus) then
               if (logPgas >= logPgas_max) then
                  alfa = 0 ! pure MESA
               else if (logPgas >= logPgas_max_minus) then
                  alfa = 1d0 - (logPgas - logPgas_max_minus)/dlogPgas
               else if (logPgas >= logPgas_min_plus) then
                  alfa = 1d0 ! pure SCVH
               else if (logPgas >= logPgas_min) then
                  alfa = 1d0 - (logPgas_min_plus - logPgas)/dlogPgas
               else ! logPgas < logPgas2
                  alfa = 0 ! pure MESA
               end if
            else if (logT >= logT_min) then
               if (logPgas >= logPgas_max) then
                  alfa = 0 ! pure MESA
               else if (logPgas >= logPgas_max_minus) then
                  dx = (logPgas - logPgas_max_minus)/dlogPgas
                  dy = (logT_min_plus - logT)/dlogT
                  alfa = 1d0 - min(1d0, sqrt(dx**2 + dy**2))
               else if (logPgas >= logPgas_min_plus) then
                  alfa = 1d0 - (logT_min_plus - logT)/dlogT
               else if (logPgas >= logPgas_min) then
                  dx = (logPgas_min_plus - logPgas)/dlogPgas
                  dy = (logT_min_plus - logT)/dlogT
                  alfa = 1d0 - min(1d0, sqrt(dx**2 + dy**2))
               else ! logPgas < logPgas_min
                  alfa = 0 ! pure MESA
               end if
            else ! logT < logT_min
               alfa = 0 ! pure MESA
            end if
         end function eval_alfa_expanded_scvh
         
         integer function eval_alfa() result(alfa)
            double precision :: &
               logT1, logT2, logT3, logT4, logW1, logW2, dx, dy
            logT1 = 5.4d0
            logT2 = 5.2d0
            logT3 = 3.5d0
            logT4 = 3.3d0
            logW1 = -4.5d0
            logW2 = -7.5d0  
            ! alfa = fraction from scvh; beta = fraction from mesa
            if (logT >= logT1) then
               alfa = 0 ! pure MESA
            else if (logT >= logT2) then
               if (logW >= logW1) then
                  alfa = (logT - logT1) / (logT2 - logT1)
               else if (logW >= logW2) then
                  dx = (logW - logW1) / (logW2 - logW1)
                  dy = (logT - logT2) / (logT1 - logT2)
                  alfa = 1 - min(1d0, sqrt(dx**2 + dy**2))
               else ! logW < logW2
                  alfa = 0 ! pure MESA
               end if
            else if (logT >= logT3) then
               if (logW >= logW1) then
                  alfa = 1 ! pure SCVH
               else if (logW >= logW2) then
                  alfa = (logW - logW2) / (logW1 - logW2)
               else ! logW < logW2
                  alfa = 0 ! pure MESA
               end if
            else if (logT >= logT4) then
               if (logW >= logW1) then
                  alfa = 1 ! pure SCVH
               else if (logW >= logW2) then
                  dx = (logW - logW2) / (logW1 - logW2)
                  dy = (logT - logT3) / (logT4 - logT3)
                  alfa = min(1d0, sqrt(dx**2 + dy**2))
               else ! logW < logW2
                  alfa = (logT - logT3) / (logT4 - logT3)
               end if
            else ! logT < logT4
               alfa = 1 ! pure SCVH
            end if
         end function eval_alfa
         
         subroutine check(str, v_scvh, v_mesa)
            character (len=*) :: str
            double precision :: v_scvh, v_mesa
            double precision :: rel_diff
            double precision, parameter :: lim = 0.2d0
            include 'formats.dek'
            rel_diff = (v_scvh - v_mesa)/max(abs(v_scvh),abs(v_mesa),1d0)
            if (abs(rel_diff) < lim) return
            write(*,*) 'rel diff too large for ' // trim(str)
            write(*,*) 'call number', call_number
            write(*,1) 'v_scvh', v_scvh
            write(*,1) 'v_mesa', v_mesa
            write(*,1) 'rel_diff', rel_diff
            write(*,*)
            write(*,1) 'logW', logW
            write(*,1) 'logT', logT
            write(*,1) 'logPgas', logPgas
            write(*,1) 'logRho', logRho
            write(*,1) 'X', X
            write(*,1) 'Z', Z
            write(*,*)
            write(*,1) 'logE', logE
            write(*,1) 'logS', logS
            write(*,1) 'chiRho', chiRho
            write(*,1) 'chiT', chiT
            write(*,1) 'Cp', Cp
            write(*,1) 'Cv', Cv
            write(*,1) 'dE_dRho', dE_dRho
            write(*,1) 'dS_dT', dS_dT
            write(*,1) 'dS_dRho', dS_dRho
            write(*,1) 'mu', mu
            write(*,1) 'log_free_e', log_free_e
            write(*,1) 'gamma1', gamma1
            write(*,1) 'gamma3', gamma3
            write(*,1) 'grad_ad', grad_ad
            write(*,1) 'eta', eta
            write(*,*)
            write(*,1) 'alfa', alfa
            write(*,1) 'beta', beta
            write(*,*)
            stop
         end subroutine check
         
      end subroutine eos_4_Pgas_T


      subroutine scvh_4_PT( &
            logW, logT, T, logPgas, Pgas, &
            abar, zbar, X, Z, logRho_guess, &
            logRho, logE, logS, chiRho, chiT, &
            Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
            mu, log_free_e, gamma1, gamma3, grad_ad, eta, call_number, ierr)
         use scvh_eval, only: interpolate_scvh
         double precision, intent(in) :: &
            logW, logT, T, logPgas, Pgas, &
            abar, zbar, X, Z, logRho_guess
         double precision, intent(out) :: &
            logRho, logE, logS, chiRho, chiT, &
            Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
            mu, log_free_e, gamma1, gamma3, grad_ad, eta
         integer, intent(in) :: call_number
         integer, intent(out) :: ierr
         
         double precision :: &
            temp,den, &
            xnh2,xnh,xnhe,xnhep, &
            pres,dpdd,dpdt, &
            ener,dedd,gam3, &
            entr,dsdd,dsdt, &
            xtra,dxdd,dxdt, &
            xnhp, xnhepp, mu_M_scvh, Ne_scvh, logNe_scvh
         integer, parameter :: include_radiation = 1
         
         include 'formats.dek'
         
         ierr = 0
         call interpolate_scvh( &
            include_radiation,logT,T,logPgas,Pgas,X, &
            logRho, logE, logS, chiRho, chiT, &
            Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
            mu, log_free_e, gamma1, gamma3, grad_ad, eta, ierr)
         
         if (ierr == 1) then ! off table
            ierr = 0
            call mesa_4_PT( &
               logW, logT, T, logPgas, Pgas, &
               abar, zbar, X, Z, logRho_guess, &
               logRho, logE, logS, chiRho, chiT, &
               Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
               mu, log_free_e, gamma1, gamma3, grad_ad, eta, call_number, ierr)
         end if
         
         if (ierr /= 0) then
            write(*,*) 'scvh_4_PT failed'
            write(*,1) 'logT', logT
            write(*,1) 'logPgas', logPgas
            stop 'scvh_4_PT'
         end if

      end subroutine scvh_4_PT


      subroutine mesa_4_PT( &
            logW, logT, T, logPgas, Pgas, &
            abar, zbar, X, Z, logRho_guess, &
            logRho, logE, logS, chiRho, chiT, &
            Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
            mu, log_free_e, gamma1, gamma3, grad_ad, eta, &
            call_number, ierr)
         double precision, intent(in) :: &
            logW, logT, T, logPgas, Pgas, &
            abar, zbar, X, Z, logRho_guess
         double precision, intent(out) :: &
            logRho, logE, logS, chiRho, chiT, &
            Cp, Cv, dE_dRho, dS_dT, dS_dRho, &
            mu, log_free_e, gamma1, gamma3, grad_ad, eta
         integer, intent(in) :: call_number
         integer, intent(out) :: ierr
         
         double precision :: &
            P, logRho_bnd1, logRho_bnd2, other, other_tol, &
            logRho_tol, other_at_bnd1, other_at_bnd2, &
            logRho_result, lnd, res(num_eos_basic_results), &
            d_eos_dlnd(num_eos_basic_results), &
            d_dabar_const_TRho(num_eos_basic_results), &
            d_dzbar_const_TRho(num_eos_basic_results), &
            d_eos_dlnT(num_eos_basic_results)
         integer:: which_other, max_iter, eos_calls
         
         include 'formats.dek'

         ierr = 0
         
         which_other = i_lnPgas
         other = logPgas*ln10
         other_tol = 1d-8*ln10
      
         max_iter = 100
            
         logRho_tol = 1d-8

         logRho_bnd1 = arg_not_provided
         logRho_bnd2 = arg_not_provided
         other_at_bnd1 = arg_not_provided
         other_at_bnd2 = arg_not_provided

         call eosDT_get_Rho( &
            eos_handle, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            logT, which_other, other, &
            logRho_tol, other_tol, max_iter, logRho_guess, &
            logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
            logRho_result, res, d_eos_dlnd, d_eos_dlnT, &
            d_dabar_const_TRho, d_dzbar_const_TRho, eos_calls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in eos_get_Rho'
            write(*,1) 'Z = ', Z
            write(*,1) 'X = ', X
            write(*,1) 'abar = ', abar
            write(*,1) 'zbar = ', zbar
            write(*,1) 'logT = ', logT
            write(*,1) 'logPgas = ', logPgas
            write(*,1) 'logRho_tol = ', logRho_tol
            write(*,1) 'other_tol = ', other_tol
            write(*,1) 'logRho_guess = ', logRho_guess
            !return
            call mesa_error(__FILE__,__LINE__)
         end if
         
         logRho = logRho_result         
         logE = res(i_lnE)/ln10
         logS = res(i_lnS)/ln10
         chiRho = res(i_chiRho)
         chiT = res(i_chiT)
         Cp = res(i_Cp)
         Cv = res(i_Cv)
         dE_dRho = res(i_dE_dRho)
         dS_dT = res(i_dS_dT)
         dS_dRho = res(i_dS_dRho)
         mu = res(i_mu)
         log_free_e = res(i_lnfree_e)/ln10
         gamma1 = res(i_gamma1)
         gamma3 = res(i_gamma3)
         grad_ad = res(i_grad_ad)
         eta = res(i_eta)
      
      end subroutine mesa_4_PT

      
      end module eval_eosPT

