! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
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


      module hydro_chem_eqns

      use star_private_def
      use const_def
      use utils_lib

      implicit none

      private
      public :: do_chem_eqns, do1_chem_eqns


      contains


      subroutine do_chem_eqns( &
            s, nvar, skip_partials, ierr)
         type (star_info), pointer :: s
         integer, intent(in) :: nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr
         integer :: k, op_err
         include 'formats'
         ierr = 0
!$OMP PARALLEL DO PRIVATE(k,op_err) SCHEDULE(dynamic,2)
         do k = 1, s% nz
            if (ierr /= 0) cycle
            call do1_chem_eqns(s, k, nvar, skip_partials, op_err)
            if (op_err /= 0) ierr = op_err
         end do
!$OMP END PARALLEL DO
      end subroutine do_chem_eqns


      subroutine do1_chem_eqns(s, k, nvar, skip_partials, ierr)

         use chem_def
         use net_lib, only: show_net_reactions, show_net_params
         use rates_def, only: reaction_Name, i_rate
         use star_utils, only: em1, e00, ep1

         type (star_info), pointer :: s
         integer, intent(in) :: k, nvar
         logical, intent(in) :: skip_partials
         integer, intent(out) :: ierr

         integer, pointer :: reaction_id(:) ! maps net reaction number to reaction id
         integer :: nz, i_lnd, i_lnT, j, i, jj, ii, equchem1, species
         real(dp) :: &
            dxdt_expected_dxa, dxdt_expected, dxdt_actual, dVARdot_dVAR, &
            dxdt_expected_dlnd, dxdt_expected_dlnT, &
            dq, dm, dequ, dxdt_nuc, dxdt_mix, max_abs_residual, &
            sum_dxdt_nuc, dx_expected_dlnd, dx_expected_dlnT, &
            d_dxdt_mix_dx00, d_dxdt_mix_dxm1, d_dxdt_mix_dxp1, &
            sum_dx_burning, sum_dx_mixing, residual, &
            dxdt_factor, alpha, eqn_scale, d_dxdt_dx, &
            dequ_dlnd, dequ_dlnT, dequ_dlnPgas_const_T, dequ_dlnT_const_Pgas
         logical :: test_partials, doing_op_split_burn
         logical, parameter :: checking = .false.

         include 'formats'

         ierr = 0

         dVARdot_dVAR = s% dVARdot_dVAR
         equchem1 = s% equchem1
         species = s% species
         nz = s% nz
         i_lnd = s% i_lnd
         i_lnT = s% i_lnT

         dq = s% dq(k)
         dm = s% dm(k)
         
         max_abs_residual = 0
         sum_dxdt_nuc = 0

         if (s% do_mix) then

            d_dxdt_mix_dxm1 = s% d_dxdt_mix_dxm1(k)
            d_dxdt_mix_dx00 = s% d_dxdt_mix_dx00(k)
            d_dxdt_mix_dxp1 = s% d_dxdt_mix_dxp1(k)

         else

            d_dxdt_mix_dxm1 = 0
            d_dxdt_mix_dx00 = 0
            d_dxdt_mix_dxp1 = 0

         end if

         sum_dx_burning = 0
         sum_dx_mixing = 0
         
         do j=1,species ! composition equation for species j in cell k

            !test_partials = (k == s% solver_test_partials_k .and. s% net_iso(ihe4) == j)
            test_partials = .false.

            i = equchem1+j-1

            dxdt_actual = s% xa_sub_xa_start(j,k)*dVARdot_dVAR
            
            doing_op_split_burn = s% op_split_burn .and. s% T_start(k) >= s% op_split_burn_min_T
            if (s% do_burn .and. .not. doing_op_split_burn) then
               dxdt_nuc = s% dxdt_nuc(j,k)
            else
               dxdt_nuc = 0
            end if

            if (s% do_mix) then
               dxdt_mix = s% dxdt_mix(j,k)
            else
               dxdt_mix = 0
            end if

            dxdt_expected = dxdt_mix + dxdt_nuc

            dxdt_factor = 1d0

            eqn_scale = max(s% min_chem_eqn_scale, s% x_scale(i,k)*dVARdot_dVAR)
            residual = (dxdt_expected - dxdt_actual)/eqn_scale
            s% equ(i,k) = residual
            
            if (abs(residual) > max_abs_residual) &
               max_abs_residual = abs(s% equ(i,k))

            if (is_bad(s% equ(i,k))) then
               s% retry_message = 'bad residual for do1_chem_eqns'
!$OMP critical (star_chem_eqns_bad_num)
               if (s% report_ierr) then
                  write(*,3) 'do1_chem_eqns: equ ' // trim(s% nameofequ(i)), &
                        i, k, s% equ(i,k)
                  write(*,2) 'dxdt_expected', k, dxdt_expected
                  write(*,2) 'dxdt_actual', k, dxdt_actual
                  write(*,2) 'eqn_scale', k, eqn_scale
                  write(*,2) 'dxdt_mix', k, dxdt_mix
                  write(*,2) 'dxdt_nuc', k, dxdt_nuc
               end if
               if (s% stop_for_bad_nums) stop 'do1_chem_eqns'
!$OMP end critical (star_chem_eqns_bad_num)
            end if

            if (skip_partials) cycle

            call e00(s, i, i, k, nvar, -dVARdot_dVAR/eqn_scale)

            ! all the rest are jacobian terms for dxdt_expected/eqn_scale

            if (s% do_burn) then

               do jj=1,species
                  ii = equchem1+jj-1
                  dxdt_expected_dxa = s% d_dxdt_nuc_dx(j,jj,k)
                  dequ = dxdt_expected_dxa/eqn_scale
                  if (checking) call check_dequ(dequ,'d_dxdt_nuc_dx')
                  call e00(s, i, ii, k, nvar, dxdt_factor*dequ)
               end do

               dxdt_expected_dlnd = s% d_dxdt_nuc_drho(j,k)*s% rho(k)
               dequ_dlnd = dxdt_expected_dlnd/eqn_scale
               dxdt_expected_dlnT = s% d_dxdt_nuc_dT(j,k)*s% T(k)
               dequ_dlnT = dxdt_expected_dlnT/eqn_scale

               if (s% do_struct_hydro) then ! partial wrt lnd const lnT
                  call e00(s, i, i_lnd, k, nvar, dxdt_factor*dequ_dlnd)
               end if

               if (s% do_struct_thermo) then ! partial wrt lnT const lnd
                  call e00(s, i, i_lnT, k, nvar, dxdt_factor*dequ_dlnT)
               end if

            end if

            if (s% do_mix) then

               dxdt_expected_dxa = d_dxdt_mix_dx00
               dequ = dxdt_expected_dxa/eqn_scale
               if (checking) call check_dequ(dequ,'d_dxdt_mix_dx00')
               call e00(s, i, i, k, nvar, dxdt_factor*dequ)
               if (k > 1) then
                  dxdt_expected_dxa = d_dxdt_mix_dxm1
                  dequ = dxdt_expected_dxa/eqn_scale
                  if (checking) call check_dequ(dequ,'d_dxdt_mix_dxm1')
                  call em1(s, i, i, k, nvar, dxdt_factor*dequ)
               end if
               if (k < nz) then
                  dxdt_expected_dxa = d_dxdt_mix_dxp1
                  dequ = dxdt_expected_dxa/eqn_scale
                  if (checking) call check_dequ(dequ,'d_dxdt_mix_dxp1')
                  call ep1(s, i, i, k, nvar, dxdt_factor*dequ)
               end if

            end if

            if (test_partials) then   
               s% solver_test_partials_dx_sink = s% net_iso(img24)
               s% solver_test_partials_val = s% dxdt_nuc(j,k)
               s% solver_test_partials_var = s% nvar_hydro + j
               s% solver_test_partials_dval_dx = s% d_dxdt_nuc_dx(j,j,k)
               write(*,*) 'do1_chem_eqns', s% solver_test_partials_var
            end if

         end do
         
         s% max_abs_xa_residual(k) = max_abs_residual
         
         contains
         
         subroutine check_dequ(dequ, str)
            real(dp), intent(in) :: dequ
            character (len=*), intent(in) :: str
            include 'formats'
            if (is_bad(dequ)) then
               ierr = -1
               if (s% report_ierr) then
                  write(*,2) 'do1_chem_eqns: bad ' // trim(str), k
               end if
               if (s% stop_for_bad_nums) stop 'do1_chem_eqns'
               return
            end if
         end subroutine check_dequ

      end subroutine do1_chem_eqns


      end module hydro_chem_eqns

