! ***********************************************************************
!
!   Copyright (C) 2018-2019  The MESA Team
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
      use star_utils, only: &
         store_rho_in_xh, store_lnd_in_xh, get_rho_and_lnd_from_xh, &
         store_T_in_xh, store_lnT_in_xh, get_T_and_lnT_from_xh
      use star_private_def
      use rsp_def, only: xa, X, Z, Y, &
         abar, zbar, z53bar, XC, XN, XO, Xne


      implicit none

      integer :: species
      integer, pointer, dimension(:) :: net_iso, chem_id
      integer :: eos_handle, kap_handle
      
      contains
      
      
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
         call mesa_error(__FILE__,__LINE__,'init_for_rsp_eos_and_kap')
      end subroutine init_for_rsp_eos_and_kap
      
      
      subroutine eval_mesa_eos_and_kap(&
            s,k,T_in,V, &
            Pg,d_Pg_dV,d_Pg_dT,Pr,d_Pr_dT, &
            egas,d_egas_dV,d_egas_dT,erad,d_erad_dV,d_erad_dT, &
            CSND,CP,CPV,CPT,Q,QV,QT,OP,OPV,OPT,ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: T_in,V
         real(dp), intent(out) :: &
            Pg,d_Pg_dV,d_Pg_dT,Pr,d_Pr_dT, &
            egas,d_egas_dV,d_egas_dT,erad,d_erad_dV,d_erad_dT, &
            CSND,CP,CPV,CPT,Q,QV,QT,OP,OPV,OPT
         integer, intent(out) :: ierr
         call eval1_mesa_eos_and_kap(&
            s,k,.false.,T_in,V, &
            Pg,d_Pg_dV,d_Pg_dT,Pr,d_Pr_dT, &
            egas,d_egas_dV,d_egas_dT,erad,d_erad_dV,d_erad_dT, &
            CSND,CP,CPV,CPT,Q,QV,QT,OP,OPV,OPT,ierr)
      end subroutine eval_mesa_eos_and_kap
      
      
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
            call store_rho_in_xh(s, k, rho)
            call get_rho_and_lnd_from_xh(s, k, s% rho(k), s% lnd(k))
            call store_T_in_xh(s, k, T)
            call get_T_and_lnT_from_xh(s, k, s% T(k), s% lnT(k))
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
            !$omp critical (rsp_eval_eos_and_kap_1)
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
            if (s% stop_for_bad_nums .and. is_bad(logRho+logT)) call mesa_error(__FILE__,__LINE__,'do_eos_for_cell')
            !$omp end critical (rsp_eval_eos_and_kap_1)
            !return
            call mesa_error(__FILE__,__LINE__,'RSP failed in get_eos')
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
      
      
      subroutine eval1_mesa_eosDEgas_and_kap( &
            s, k, skip_kap, egas, V, T, Pgas, CSND, CP, Q, OP, ierr)
         use star_utils, only: write_eos_call_info
         use eos_support, only: solve_eos_given_DEgas
         use micro, only: store_eos_for_cell
         type (star_info), pointer :: s
         integer, intent(in) :: k
         logical, intent(in) :: skip_kap
         real(dp), intent(in) :: egas, V
         real(dp), intent(out) :: T, Pgas, CSND, CP, Q, OP
         integer, intent(out) :: ierr
         
         integer :: j, eos_calls
         real(dp) :: rho, logRho, dlnd_dV, egas_tol, logT, &
            logT_guess, logT_tol, new_erad, new_egas, OPV, OPT
         real(dp), dimension(num_eos_basic_results) :: &
            res, d_dlnd, d_dlnT
         real(dp) :: d_dxa(num_eos_d_dxa_results, species)
         
         include 'formats'
         ierr = 0

         rho = 1d0/V
         logRho = log10(rho)
         dlnd_dV = -rho
         
         if (egas <= 0d0 .or. is_bad(egas)) then
            !$OMP critical (RSP_eosDEgas)
            write(*,2) 'egas', k, egas
            write(*,*) 'called eval1_mesa_eosDEgas_and_kap with bad value for egas'
            call mesa_error(__FILE__,__LINE__,'eval1_mesa_eosDEgas_and_kap')
            !$OMP end critical (RSP_eosDEgas)
            ierr = -1
            return
         end if

         if (k > 0 .and. k <= s% nz) then
            egas_tol = egas*1d-11
            s% rho(k) = rho
            s% lnd(k) = logRho*ln10
            s% xh(s% i_lnd,k) = s% lnd(k)
            s% abar(k) = abar
            s% zbar(k) = zbar
            logT_guess = s% lnT(k)/ln10
            logT_tol = 1d-11
            call solve_eos_given_DEgas( &
               s, k, xa, &
               logRho, egas, logT_guess, logT_tol, egas_tol, &
               logT, res, d_dlnd, d_dlnT, d_dxa, &
               ierr)
            if (ierr == 0) then
               call store_lnT_in_xh(s, k, logT*ln10)
               call get_T_and_lnT_from_xh(s, k, s% T(k), s% lnT(k))
               T = s% T(k)
               new_erad = crad*T**4/rho
               new_egas = exp(res(i_lnE)) - new_erad
               if (is_bad(new_egas) .or. new_egas <= 0d0 .or. &
                     abs(new_egas - egas) > 1d3*egas_tol) then               
                  !$OMP critical (RSP_eosDEgas)
                  write(*,1) 'logRho', s% lnd(k)/ln10
                  write(*,1) 'logT_guess', logT_guess
                  write(*,1) 'egas', egas
                  write(*,1) 'Z', Z
                  write(*,1) 'X', X
                  write(*,1) 'abar', abar
                  write(*,1) 'zbar', zbar
                  write(*,1) 'logT_tol', logT_tol
                  write(*,1) 'egas_tol', egas_tol
                  write(*,'(A)')
                  write(*,1) 'guess logT', logT_guess
                  write(*,1) 'found logT', logT
                  write(*,1) 'wanted egas', egas
                  write(*,1) 'got egas', new_egas
                  write(*,1) '(want - got)/got', (egas - new_egas)/new_egas
                  write(*,'(A)')
                  write(*,2) 'eos_calls', eos_calls
                  write(*,'(A)')
                  write(*,2) 'failed eval1_mesa_eosDEgas_and_kap', k
                  write(*,'(A)')
                  write(*,*) 'is_bad(new_egas)', is_bad(new_egas)
                  write(*,*) 'new_egas <= 0d0', new_egas <= 0d0
                  write(*,*) 'abs(new_egas - egas) > egas_tol', &
                     abs(new_egas - egas) > egas_tol, &
                     abs(new_egas - egas) - egas_tol, &
                     abs(new_egas - egas), egas_tol
                  call mesa_error(__FILE__,__LINE__,'eval1_mesa_eosDEgas_and_kap')
                  !$OMP end critical (RSP_eosDEgas)
               end if
               s% lnPgas(k) = res(i_lnPgas)
               s% Pgas(k) = exp(s% lnPgas(k))
               Pgas = s% Pgas(k)
               do j=1,num_eos_basic_results
                  s% d_eos_dlnd(j,k) = d_dlnd(j)
                  s% d_eos_dlnT(j,k) = d_dlnT(j)
               end do
               call store_eos_for_cell(s, k, res, d_dlnd, d_dlnT, d_dxa, ierr)
               CSND = s% csound(k)
            end if
         else ! k <= 0 or k > nz
            write(*,*) 'cannot call eval1_mesa_eosDEgas_and_kap with k <= 0 or k > nz'
            ierr = -1
            return
         end if
                   
         if (ierr /= 0) then
            !$OMP critical (RSP_eosDEgas)
            if (k > 0 .and. k < s% nz) call write_eos_call_info(s,k)
            write(*,2) 'X', k, X
            write(*,2) 'Z', k, Z
            write(*,2) 'zbar', k, zbar
            write(*,2) 'abar', k, abar
            write(*,2) 'V', k, V
            write(*,2) 'rho', k, rho
            write(*,2) 'egas', k, egas
            write(*,2) 'logRho', k, logRho
            write(*,2) 'T', k, T
            write(*,2) 'logT', k, logT
            if (s% stop_for_bad_nums .and. egas <= 0d0) call mesa_error(__FILE__,__LINE__,'do_eos_for_cell')
            !$OMP end critical (RSP_eosDEgas)
            return
            call mesa_error(__FILE__,__LINE__,'RSP failed in eval1_mesa_eosDEgas_and_kap')
         end if
         
         if (skip_kap) then
            OP=0
         else
            call eval1_kap(s, k, skip_kap, &
               V, logRho, dlnd_dV, T, logT, species, chem_id, net_iso, xa, &
               res, d_dlnd, d_dlnT, OP, OPV, OPT, ierr)
         end if
         
         Pgas = exp(res(i_lnPgas))
         CP = res(i_Cp)
         Q = res(i_chiT)/(rho*T*res(i_chiRho)) ! thermal expansion coefficient
         
      end subroutine eval1_mesa_eosDEgas_and_kap
      
      
      subroutine eval1_mesa_eosDE_and_kap( &  ! for eos, energy = egas + erad
            s, k, skip_kap, energy, V, T, Pgas, CSND, CP, Q, OP, ierr)
         use star_utils, only: write_eos_call_info
         use eos_support, only: solve_eos_given_DE
         use micro, only: store_eos_for_cell
         type (star_info), pointer :: s
         integer, intent(in) :: k
         logical, intent(in) :: skip_kap
         real(dp), intent(in) :: energy, V
         real(dp), intent(out) :: T, Pgas, CSND, CP, Q, OP
         integer, intent(out) :: ierr
         
         integer :: j, eos_calls
         real(dp) :: rho, logRho, dlnd_dV, logE, logE_want, logE_tol, logT, &
            logT_guess, logT_tol, new_erad, new_egas, OPV, OPT
         real(dp), dimension(num_eos_basic_results) :: &
            res, d_dlnd, d_dlnT
         real(dp) :: d_dxa(num_eos_d_dxa_results, species)
         
         include 'formats'
         ierr = 0

         rho = 1d0/V
         logRho = log10(rho)
         logE_want = log10(energy)
         dlnd_dV = -rho
         
         if (energy <= 0d0 .or. is_bad(energy)) then
            !$OMP critical (RSP_eosDE)
            write(*,2) 'energy', k, energy
            write(*,*) 'called eval1_mesa_eosDE_and_kap with bad value for energy'
            call mesa_error(__FILE__,__LINE__,'eval1_mesa_eosDE_and_kap')
            !$OMP end critical (RSP_eosDE)
            ierr = -1
            return
         end if

         if (k > 0 .and. k <= s% nz) then
            s% rho(k) = rho
            s% lnd(k) = logRho*ln10
            s% xh(s% i_lnd,k) = s% lnd(k)
            s% abar(k) = abar
            s% zbar(k) = zbar
            logT_guess = s% lnT(k)/ln10
            logT_tol = 1d-11
            logE_tol = 1d-11
            
            call solve_eos_given_DE( &
               s, k, xa, &
               logRho, logE_want, logT_guess, logT_tol, logE_tol, &
               logT, res, d_dlnd, d_dlnT, d_dxa, &
               ierr)
            if (ierr == 0) then
               call store_lnT_in_xh(s, k, logT*ln10)
               call get_T_and_lnT_from_xh(s, k, s% T(k), s% lnT(k))
               T = s% T(k)
               new_erad = crad*T**4/rho
               new_egas = exp(res(i_lnE)) - new_erad
               logE = res(i_lnE)/ln10
               if (is_bad(logE) .or. logE <= -20d0) then               
                  !$OMP critical (RSP_eosDE)
                  write(*,1) 'logRho', s% lnd(k)/ln10
                  write(*,1) 'Z', Z
                  write(*,1) 'X', X
                  write(*,1) 'abar', abar
                  write(*,1) 'zbar', zbar
                  write(*,1) 'logT_tol', logT_tol
                  write(*,1) 'logE_tol', logE_tol
                  write(*,'(A)')
                  write(*,1) 'guess logT', logT_guess
                  write(*,1) 'found logT', logT
                  write(*,'(A)')
                  write(*,1) 'wanted logE', logE_want
                  write(*,1) 'got logE', logE
                  write(*,'(A)')
                  write(*,2) 'eos_calls', eos_calls
                  write(*,'(A)')
                  write(*,2) 'failed eval1_mesa_eosDE_and_kap', k
                  write(*,'(A)')
                  call mesa_error(__FILE__,__LINE__,'eval1_mesa_eosDE_and_kap')
                  !$OMP end critical (RSP_eosDE)
               end if
               s% lnPgas(k) = res(i_lnPgas)
               s% Pgas(k) = exp(s% lnPgas(k))
               Pgas = s% Pgas(k)
               do j=1,num_eos_basic_results
                  s% d_eos_dlnd(j,k) = d_dlnd(j)
                  s% d_eos_dlnT(j,k) = d_dlnT(j)
               end do
               call store_eos_for_cell(s, k, res, d_dlnd, d_dlnT, d_dxa, ierr)
               CSND = s% csound(k)
            end if
         else ! k <= 0 or k > nz
            write(*,*) 'cannot call eval1_mesa_eosDE_and_kap with k <= 0 or k > nz'
            ierr = -1
            return
         end if
                   
         if (ierr /= 0) then
            !$OMP critical (RSP_eosDE)
            if (k > 0 .and. k < s% nz) call write_eos_call_info(s,k)
            write(*,2) 'X', k, X
            write(*,2) 'Z', k, Z
            write(*,2) 'zbar', k, zbar
            write(*,2) 'abar', k, abar
            write(*,2) 'V', k, V
            write(*,2) 'rho', k, rho
            write(*,2) 'logE', k, logE
            write(*,2) 'logRho', k, logRho
            write(*,2) 'T', k, T
            write(*,2) 'logT', k, logT
            !$OMP end critical (RSP_eosDE)
            return
            call mesa_error(__FILE__,__LINE__,'RSP failed in eval1_mesa_eosDE_and_kap')
         end if
         
         if (skip_kap) then
            OP=0
         else
            call eval1_kap(s, k, skip_kap, &
               V, logRho, dlnd_dV, T, logT, species, chem_id, net_iso, xa, &
               res, d_dlnd, d_dlnT, OP, OPV, OPT, ierr)
         end if
         
         Pgas = exp(res(i_lnPgas))
         CP = res(i_Cp)
         Q = res(i_chiT)/(rho*T*res(i_chiRho)) ! thermal expansion coefficient
         
      end subroutine eval1_mesa_eosDE_and_kap
      
      
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
!$omp critical (rsp_eval_eos_and_kap_2)
            write(*,*) 'failed in eval1_mesa_eos_and_kap get kap'
            write(*,2) 'logRho', k, logRho
            write(*,2) 'logT', k, logT
            write(*,2) 'lnfree_e', k, lnfree_e
            write(*,2) 'zbar', k, zbar
            write(*,2) 'X', k, X
            write(*,2) 'Z', k, Z
            call mesa_error(__FILE__,__LINE__,'eval1_kap')
!$omp end critical (rsp_eval_eos_and_kap_2)
            !return
            call mesa_error(__FILE__,__LINE__)
         end if

         if (k > 0 .and. k <= s% nz .and. s% use_other_opacity_factor) then
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
      
      
      subroutine eval1_mesa_T_given_DP(s, k, Vol, P, T_guess, T, ierr)
         use eos_support, only: solve_eos_given_DP
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(in) :: Vol, P, T_guess
         real(dp), intent(out) :: T
         integer, intent(out) :: ierr         
         real(dp) :: rho, logRho, logP, logT_guess, &
            logT_tol, logP_tol, logT
         real(dp), dimension(num_eos_basic_results) :: &
            res, d_dlnd, d_dlnT
         real(dp) :: d_dxa(num_eos_d_dxa_results, species)
         include 'formats'
         ierr = 0
         rho = 1d0/Vol
         logRho = log10(rho)
         logP = log10(P)
         logT_guess = log10(T_guess)
         logT_tol = 1d-11
         logP_tol = 1d-11
         call solve_eos_given_DP( &
            s, k, xa, &
            logRho, logP, logT_guess, logT_tol, logP_tol, &
            logT, res, d_dlnd, d_dlnT, d_dxa, &
            ierr)
         T = exp10(logT)
      end subroutine eval1_mesa_T_given_DP
      
      
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
            res, d_dlnd, d_dlnT
         real(dp) :: d_dxa(num_eos_d_dxa_results, species)
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
                 s, k, xa, &
                 logT, logP, logRho_guess, logRho_tol, logP_tol, &
                 logRho, res, d_dlnd, d_dlnT, d_dxa, &
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


      subroutine update_eos_and_kap(s,kk,ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: kk  
         integer, intent(out) :: ierr  
         real(dp) :: &
            T,V,Pgas,d_Pg_dV,d_Pg_dT,Prad,d_Pr_dT, &
            egas,d_egas_dV,d_egas_dT,erad,d_erad_dV,d_erad_dT, &
            CSND,CP,CPV,CPT,Q,QV,QT,OP,OPV,OPT
         T = s% T(kk)
         V = 1d0/s% rho(kk)
         call eval_mesa_eos_and_kap( &
            s,kk,T,V, &
            Pgas,d_Pg_dV,d_Pg_dT,Prad,d_Pr_dT, &
            egas,d_egas_dV,d_egas_dT,erad,d_erad_dV,d_erad_dT, &
            CSND,CP,CPV,CPT,Q,QV,QT,OP,OPV,OPT,ierr)       
      end subroutine update_eos_and_kap


      subroutine set_Rho_for_new_Pgas(s,kk,ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: kk  
         integer, intent(out) :: ierr  
         real(dp) :: &
            other_value, other_tol, logRho_bnd1, other_at_bnd1, &
            logRho_bnd2, other_at_bnd2, logRho_guess, logRho_result, logRho_tol
         real(dp), dimension(num_eos_basic_results) :: &
            res, d_dlnd, d_dlnT
         real(dp), dimension(num_eos_d_dxa_results, species) :: d_dxa
         integer :: max_iter, which_other, eos_calls, iter
            
         include 'formats'
         ierr = 0         
         max_iter = 100
         which_other = i_lnPgas
         other_value = log(s% Pgas(kk))
         other_tol = 1d-11
         logRho_tol = 1d-11
         logRho_bnd1 = arg_not_provided
         other_at_bnd1 = arg_not_provided
         logRho_bnd2 = arg_not_provided
         other_at_bnd2 = arg_not_provided
         logRho_guess = log10(s% rho(kk))
         call store_T_in_xh(s, kk, s% T(kk))
         call get_T_and_lnT_from_xh(s, kk, s% T(kk), s% lnT(kk))
         
         call eosDT_get_Rho( &
            eos_handle, &
            species, chem_id, net_iso, xa, &
            s% lnT(kk)/ln10, which_other, other_value, &
            logRho_tol, other_tol, max_iter, logRho_guess, &
            logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
            logRho_result, res, d_dlnd, d_dlnT, &
            d_dxa, eos_calls, ierr)
         if (ierr /= 0) return
               
         s% lnd(kk) = logRho_result*ln10
         s% xh(s% i_lnd,kk) = s% lnd(kk)
         s% rho(kk) = exp(s% lnd(kk))
         s% Vol(kk) = 1d0/s% rho(kk)
         
      end subroutine set_Rho_for_new_Pgas


      subroutine set_T_for_new_egas(s,kk,ierr)  ! uses s% T(kk), s% egas(kk) and s% lnd(kk)
         type (star_info), pointer :: s
         integer, intent(in) :: kk  
         integer, intent(out) :: ierr  

         real(dp) :: &
            egas_tol, logT_bnd1, egas_at_bnd1, new_egas, egas_want, &
            logT_bnd2, egas_at_bnd2, logT_guess, logT_result, logT_tol
         real(dp), dimension(num_eos_basic_results) :: &
            res, d_dlnd, d_dlnT
         real(dp) :: d_dxa(num_eos_d_dxa_results, species)
         integer :: max_iter, which_other, eos_calls, iter
            
         include 'formats'
         ierr = 0         
         max_iter = 100
         egas_want = s% egas(kk)
         egas_tol = egas_want*1d-12
         logT_tol = 1d-11
         logT_bnd1 = arg_not_provided
         egas_at_bnd1 = arg_not_provided
         logT_bnd2 = arg_not_provided
         egas_at_bnd2 = arg_not_provided
         logT_guess = log10(s% T(kk))
         
         call eosDT_get_T( &
            eos_handle, &
            species, chem_id, net_iso, xa, &
            s% lnd(kk)/ln10, i_egas, egas_want, &
            logT_tol, egas_tol, max_iter, logT_guess, &
            logT_bnd1, logT_bnd2, egas_at_bnd1, egas_at_bnd2, &
            logT_result, res, d_dlnd, d_dlnT, &
            d_dxa, eos_calls, ierr)
         if (ierr /= 0) return
               
         call store_lnT_in_xh(s, kk, logT_result*ln10)
         call get_T_and_lnT_from_xh(s, kk, s% T(kk), s% lnT(kk))
         
         new_egas = exp(res(i_lnE)) - crad*s% T(kk)**4/s% rho(kk)
         if (is_bad(new_egas) .or. new_egas <= 0d0 .or. &
               abs(new_egas - egas_want) > 1d3*egas_tol) then               
            write(*,1) 'logRho', s% lnd(kk)/ln10
            write(*,1) 'logT_guess', logT_guess
            write(*,1) 'egas_want', egas_want
            write(*,1) 'Z', Z
            write(*,1) 'X', X
            write(*,1) 'abar', abar
            write(*,1) 'zbar', zbar
            write(*,1) 'logT_tol', logT_tol
            write(*,1) 'egas_tol', egas_tol
            write(*,'(A)')
            write(*,1) 'guess logT', logT_guess
            write(*,1) 'found logT', logT_result
            write(*,1) 'wanted egas', egas_want
            write(*,1) 'got egas', new_egas
            write(*,1) '(want - got)/got', (egas_want - new_egas)/new_egas
            write(*,'(A)')
            write(*,*) 'eos_calls', eos_calls
            write(*,'(A)')
            write(*,2) 'failed set_T_for_new_egas', kk
            write(*,'(A)')
            write(*,*) 'is_bad(new_egas)', is_bad(new_egas)
            write(*,*) 'new_egas <= 0d0', new_egas <= 0d0
            write(*,*) 'abs(new_egas - egas_want) > egas_tol', &
               abs(new_egas - egas_want) > egas_tol, &
               abs(new_egas - egas_want) - egas_tol, &
               abs(new_egas - egas_want), egas_tol
            call mesa_error(__FILE__,__LINE__,'set_T_for_new_egas')
         end if
         
      end subroutine set_T_for_new_egas


      subroutine set_T_for_new_Pgas(s,kk,ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: kk  
         integer, intent(out) :: ierr  
         real(dp) :: &
            other_value, other_tol, logT_bnd1, other_at_bnd1, &
            logT_bnd2, other_at_bnd2, logT_guess, logT_result, logT_tol
         real(dp), dimension(num_eos_basic_results) :: &
            res, d_dlnd, d_dlnT
         real(dp), dimension(num_eos_d_dxa_results, species) :: d_dxa
         integer :: max_iter, which_other, eos_calls, iter
            
         include 'formats'
         ierr = 0         
         max_iter = 100
         which_other = i_lnPgas
         other_value = log(s% Pgas(kk))
         other_tol = 1d-11
         logT_tol = 1d-11
         logT_bnd1 = arg_not_provided
         other_at_bnd1 = arg_not_provided
         logT_bnd2 = arg_not_provided
         other_at_bnd2 = arg_not_provided
         logT_guess = log10(s% T(kk))
         
         call eosDT_get_T( &
            eos_handle, &
            species, chem_id, net_iso, xa, &
            s% lnd(kk)/ln10, which_other, other_value, &
            logT_tol, other_tol, max_iter, logT_guess, &
            logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
            logT_result, res, d_dlnd, d_dlnT, &
            d_dxa, eos_calls, ierr)
         if (ierr /= 0) return
               
         call store_lnT_in_xh(s, kk, logT_result*ln10)
         call get_T_and_lnT_from_xh(s, kk, s% T(kk), s% lnT(kk))

      end subroutine set_T_for_new_Pgas


      subroutine set_T_for_new_energy(s,kk,logT_tol,other_tol,ierr)
         use eos_lib, only: eosDT_get_T
         type (star_info), pointer :: s
         integer, intent(in) :: kk  
         real(dp), intent(in) :: logT_tol, other_tol  
         integer, intent(out) :: ierr  
         real(dp) :: &
            other_value, logT_bnd1, other_at_bnd1, &
            logT_bnd2, other_at_bnd2, logT_guess, logT_result
         real(dp), dimension(num_eos_basic_results) :: &
            res, d_dlnd, d_dlnT
         real(dp), dimension(num_eos_d_dxa_results, species) :: d_dxa
         integer :: max_iter, which_other, eos_calls, iter
         ierr = 0         
         max_iter = 100
         which_other = i_lnE
         other_value = log(s% energy(kk))
         !other_tol = 1d-12
         !logT_tol = 1d-12
         logT_bnd1 = arg_not_provided
         other_at_bnd1 = arg_not_provided
         logT_bnd2 = arg_not_provided
         other_at_bnd2 = arg_not_provided
         logT_guess = log10(s% T(kk))
         
         call eosDT_get_T( &
            eos_handle, &
            species, chem_id, net_iso, xa, &
            s% lnd(kk)/ln10, which_other, other_value, &
            logT_tol, other_tol, max_iter, logT_guess, &
            logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
            logT_result, res, d_dlnd, d_dlnT, &
            d_dxa, eos_calls, ierr)
         if (ierr /= 0) return
               
         call store_lnT_in_xh(s, kk, logT_result*ln10)
         call get_T_and_lnT_from_xh(s, kk, s% T(kk), s% lnT(kk))
      end subroutine set_T_for_new_energy                


      end module rsp_eval_eos_and_kap
      

