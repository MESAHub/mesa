! ***********************************************************************
!
!   Copyright (C) 2010-2019  Josiah Schwab & The MESA Team
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

      module net_derivs_support
      use net_def, only: Net_General_Info, Net_Info, net_test_partials, &
         net_test_partials_val, net_test_partials_dval_dx, net_test_partials_i, &
         net_test_partials_iother
      use const_def
      use chem_def
      use rates_def

      implicit none


      real(dp), parameter :: r_min = 1d-99

      logical, parameter :: checking_deriv_flags = .false.


      logical, parameter :: show_rate = .false.
      logical, parameter :: show_jac = .false.
      logical, parameter :: show_neuQs = .false.

      real(dp), parameter :: show_dydt_y = 7.0097206738032283D-02
      logical, parameter :: show_dydt = .false.

      logical, parameter :: show_d_dydt_dRho = .false.
      logical, parameter :: show_d_dydt_dT = .false.

      logical, parameter :: show_eps_nuc = .false.
      logical, parameter :: show_d_eps_nuc_dy = .false.

      logical, parameter :: checkQs = .false.
      real(dp), parameter :: checkQ_frac = 1d-4


      contains


      real(dp) function isoB(ci)
         use chem_lib, only: get_Q
         integer, intent(in) :: ci

         isoB = get_Q(chem_isos,ci)

      end function isoB


      subroutine do_in_out( &
           n, dydt, eps_nuc_MeV, i, r_in, &
           n_in, i_in, c_in, n_out, i_out, c_out, &
           idr, dr, &
           deriv_flgs, symbolic, just_dydt)

        type (Net_Info) :: n
        real(qp) :: dydt(:,:) ! (num_rvs, num_isos)
        real(qp), intent(out) :: eps_nuc_MeV(num_rvs)
        integer, intent(in) :: i ! the reaction number
        real(dp), intent(in) :: r_in ! coefficient of rate for the reaction
        integer, intent(in) :: n_in, n_out ! number of inputs and outputs
        integer, dimension(3), intent(in) :: i_in, i_out ! net isotope numbers for the reaction
        real(dp), dimension(3), intent(in) :: c_in, c_out ! isotope coefficients in reaction equation
        integer, dimension(3), intent(in) :: idr ! isotope number for dr
        real(dp), dimension(3), intent(in) :: dr ! coefficient for Jacobian entries d_dydt_dy(idr)
        logical, pointer :: deriv_flgs(:)
        logical, intent(in) :: symbolic, just_dydt

        ! the purpose of this wrapper is to automatically set the Q values associated with the reaction

        integer :: j
        j = n% g% reaction_id(i)
        call do_in_out_neu( &
             n, dydt, eps_nuc_MeV, i, r_in, &
             n_in, i_in, c_in, n_out, i_out, c_out, &
             idr, dr, &
             n% reaction_Qs(j), n% reaction_neuQs(j), 0d0, 0d0, &
             deriv_flgs, symbolic, just_dydt)

      end subroutine do_in_out


      subroutine do_in_out_neu( &
           n, dydt, eps_nuc_MeV, i, r_in, &
           n_in, i_in, c_in, n_out, i_out, c_out, &
           idr, dr, Q, Qneu, dQneu_dT, dQneu_dRho, &
           deriv_flgs, symbolic, just_dydt)

        ! this function handles reactions with 1-3 inputs going to 1-3 outputs

        type (Net_Info) :: n
        real(qp) :: dydt(:,:) ! (num_rvs, num_isos)
        real(qp), intent(out) :: eps_nuc_MeV(num_rvs)
        integer, intent(in) :: i ! the reaction number
        real(dp), intent(in) :: r_in ! coefficient of rate for the reaction
        integer, intent(in) :: n_in, n_out ! number of inputs and outputs
        integer, dimension(3), intent(in) :: i_in, i_out ! net isotope numbers for the reaction
        real(dp), dimension(3), intent(in) :: c_in, c_out ! isotope coefficients in reaction equation
        integer, dimension(3), intent(in) :: idr ! isotope number for dr
        real(dp), dimension(3), intent(in) :: dr ! coefficient for Jacobian entries d_dydt_dy(idr)
        real(dp), intent(in) :: Q, Qneu, dQneu_dT, dQneu_dRho
        logical, pointer :: deriv_flgs(:)
        logical, intent(in) :: symbolic, just_dydt

        real(dp) :: rvs(num_rvs), d, d1, d2, d3, lhs, rhs, r, checkQ
        type (Net_General_Info), pointer  :: g
        integer, pointer :: chem_id(:)
        integer :: j, cid, icat, reaction_id

        logical :: condition ! for debugging output

        include 'formats'

        ! enforce minimum rate threshold
        r = r_in
        if (r < r_min .or. n% rate_screened(i) < r_min) r = 0

        ! identify reaction category
        icat = reaction_categories(n% g% reaction_id(i))

        g => n% g
        chem_id => g% chem_id

        d = n% rate_screened(i)
        d1  = dr(1) * d
        d2  = dr(2) * d
        d3  = dr(3) * d
        rvs(i_rate) = r * n% rate_screened(i)
        rvs(i_rate_dT) = r * n% rate_screened_dT(i)
        rvs(i_rate_dRho) = r * n% rate_screened_dRho(i)

        n% raw_rate(i) = n% rate_raw(i) * r * avo
        n% screened_rate(i) = n% rate_screened(i) * r * avo
        n% eps_nuc_rate(i) = n% rate_screened(i) * r * (Q - Qneu)  * Qconv
        n% eps_neu_rate(i) = n% rate_screened(i) * r * Qneu  * Qconv

        ! evaluate left hand side (inputs)
        lhs = 0
        do j = 1, n_in
           call check(i_in(j), 'input', j)
           cid = chem_id(i_in(j))
           call do_lhs_iso(n, dydt, i, c_in(j), i_in(j), rvs, idr(1), d1, idr(2), d2, idr(3), d3, &
                symbolic, just_dydt)
           lhs = lhs + c_in(j)*(chem_isos% Z(cid) + chem_isos% N(cid))
        end do

        ! evaluate right hand side (outputs)
        rhs = 0
        do j = 1, n_out
           call check(i_out(j), 'output', j)
           cid = chem_id(i_out(j))
           call do_rhs_iso(n, dydt, i, c_out(j), i_out(j), rvs, idr(1), d1, idr(2), d2, idr(3), d3, &
                symbolic, just_dydt)
           rhs = rhs + c_out(j)*(chem_isos% Z(cid) + chem_isos% N(cid))
        end do

        ! construct eps_nuc and its Rho & T derivatives
        eps_nuc_MeV(i_rate) = eps_nuc_MeV(i_rate) + rvs(i_rate)*(Q-Qneu)
        eps_nuc_MeV(i_rate_dT) = eps_nuc_MeV(i_rate_dT) + rvs(i_rate_dT)*(Q-Qneu) - rvs(i_rate)*dQneu_dT
        eps_nuc_MeV(i_rate_dRho) = eps_nuc_MeV(i_rate_dRho) + rvs(i_rate_dRho)*(Q-Qneu) - rvs(i_rate)*dQneu_dRho


        ! for debugging
        condition = abs(rvs(1)*Q) > 1d2
        if (show_eps_nuc .and. condition) &
             write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' eps_nuc',  rvs(1)*Q


        ! set categories and neutrinos
        n% eps_nuc_categories(icat) = n% eps_nuc_categories(icat) + rvs(i_rate)*(Q-Qneu)
        n% eps_neu_total = n% eps_neu_total + rvs(i_rate) * Qneu

        ! for debugging
        condition = n% reaction_neuQs(g% reaction_id(i))*rvs(i_rate) /= 0 .and. &
             abs(n% y(g% net_iso(ihe4)) - show_dydt_y) < 1d-20
        if (show_neuQs .and. condition)  &
             write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' neu',  &
             n% reaction_neuQs(reaction_id)*rvs(:)


        ! set composition derivatives
        !write(*,*) trim(reaction_Name(g% reaction_id(i)))
        !write(*,*) idr(1), idr(2), idr(3), lbound(n% d_eps_nuc_dy),ubound(n% d_eps_nuc_dy),d1,d2,d3
        if(idr(1)>0) n% d_eps_nuc_dy(idr(1)) = n% d_eps_nuc_dy(idr(1)) + d1*(Q-Qneu)
        if(idr(2)>0) n% d_eps_nuc_dy(idr(2)) = n% d_eps_nuc_dy(idr(2)) + d2*(Q-Qneu)
        if(idr(3)>0) n% d_eps_nuc_dy(idr(3)) = n% d_eps_nuc_dy(idr(3)) + d3*(Q-Qneu)
    
        ! for debugging
        if(idr(1)>0) then
            condition = chem_id(idr(1)) == ic12
            if (show_d_eps_nuc_dy .and. d1 > 0 .and. condition)  &
                write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' d_epsnuc_dy',  &
                d, dr(1), d1, Q, d1*Q, n% d_eps_nuc_dy(idr(1)), n% reaction_Qs(reaction_id), &
                n% reaction_neuQs(reaction_id)
        end if
    
        if(idr(2)>0) then
            condition = chem_id(idr(2)) == ic12
            if (show_d_eps_nuc_dy .and. d2 > 0 .and. condition)  &
                write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' d_epsnuc_dy',  &
                d, dr(2), d2, Q, d2*Q, n% d_eps_nuc_dy(idr(2)), n% reaction_Qs(reaction_id), &
                n% reaction_neuQs(reaction_id)
        end if
    
        if(idr(3)>0) then
            condition = chem_id(idr(3)) == ic12
            if (show_d_eps_nuc_dy .and. d3 > 0 .and. condition)  &
                write(*,1) trim(reaction_Name(g% reaction_id(i))) // ' d_epsnuc_dy',  &
                d, dr(3), d3, Q, d3*Q, n% d_eps_nuc_dy(idr(3)), n% reaction_Qs(reaction_id), &
                n% reaction_neuQs(reaction_id)
        end if

        call check_balance(n, i, lhs, rhs)
        if (checking_deriv_flags) deriv_flgs(i) = .true.

        if (checkQs) then
           checkQ = 0
           do j = 1, n_out
              cid = chem_id(i_out(j))
              checkQ = checkQ + c_out(j)*isoB(cid)
           end do
           do j = 1, n_in
              cid = chem_id(i_in(j))
              checkQ = checkQ - c_in(j)*isoB(cid)
           end do
           if (abs(Q - checkQ) > checkQ_frac*abs(checkQ)) then
              write(*,1) 'do_in_out checkQ ' // trim(reaction_Name(g% reaction_id(i))),  &
                   Q, checkQ
              !stop
           end if
        end if

      contains

        subroutine check(ii, str, jj)
          integer, intent(in) :: ii, jj
          character (len=*), intent(in) :: str
          if (ii <= 0) then
             write(*,*)  &
                  'do_in_out: bad iso num for ' // trim(str), jj, &
                  ' in ' // trim(reaction_Name(g% reaction_id(i)))
             call mesa_error(__FILE__,__LINE__)
          end if
        end subroutine check

      end subroutine do_in_out_neu


      subroutine do_lhs_iso( &
            n, dydt, i, c, i1, rvs, i2, dr2, i3, dr3, i4, dr4, symbolic, just_dydt)
         type (Net_Info) :: n
         real(qp) :: dydt(:,:) ! (num_rvs, num_isos)
         integer, intent(in) :: i, i1, i2, i3, i4
         real(dp), intent(in) :: c, rvs(:), dr2, dr3, dr4
         logical, intent(in) :: symbolic, just_dydt

         ! i1, i2, i3, and 14 are isotope numbers
         ! dr2 = dr/dy(i2); dr3 = dr/dy(i3); dr4 = dr/dy(i4)
         ! -c * r   = dydt(i1)
         ! -c * dr2 = d_dydt(i1)_dy(i2)
         ! -c * dr3 = d_dydt(i1)_dy(i3)
         ! -c * dr4 = d_dydt(i1)_dy(i4)

         integer :: j
         integer, pointer :: chem_id(:)
         chem_id => n% g% chem_id

         include 'formats'

         if (symbolic) then
            n% d_dydt_dy(i1, i2) = 1
            if (i3 <= 0) return
            n% d_dydt_dy(i1, i3) = 1
            if (i4 <= 0) return
            n% d_dydt_dy(i1, i4) = 1
            return
         end if

         ! update the dydt terms for i1
         do j=1,num_rvs
            dydt(j,i1) = dydt(j,i1) - c * rvs(j)
         end do
         if (chem_id(i1) == io16 .and. show_dydt .and. &
                  abs(n% y(i1) - show_dydt_y) < 1d-20) &
               write(*,1) 'lhs ' // trim(reaction_Name(n% g% reaction_id(i))), &
                  -c * rvs(i_rate), dydt(i_rate,i1), &
                  n% rate_screened(i), &
                  n% rate_raw(i), &
                  n% y(i1)
         if (chem_id(i1) == ihe4 .and. show_d_dydt_dRho .and. &
                  abs(n% y(i1) - show_dydt_y) < 1d-20) &
               write(*,1) trim(reaction_Name(n% g% reaction_id(i))) // ' d_dydt_dRho ' &
                  // trim(chem_isos% name(chem_id(i1))), &
                  -c * rvs(i_rate_dRho), n% rate_screened(i), &
                  n% rate_screened_dRho(i), n% y(i1)
         if (chem_id(i1) == ihe4 .and. show_d_dydt_dT .and. &
                  abs(n% y(i1) - show_dydt_y) < 1d-20) &
               write(*,1) trim(reaction_Name(n% g% reaction_id(i))) // ' d_dydt_dT ' &
                  // trim(chem_isos% name(chem_id(i1))), &
                  -c * rvs(i_rate_dT), n% rate_screened(i), &
                  n% rate_screened_dT(i), n% y(i1)

         if (just_dydt) return

         ! update the Jacobian for d_dydt(i1)_dy(i2)
         n% d_dydt_dy(i1, i2) = n% d_dydt_dy(i1, i2)  - c * dr2

         if (chem_id(i1) == ini56 .and. show_jac .and. c * dr2 /= 0) &
               write(*,1) trim(reaction_Name(n% g% reaction_id(i))) // ' d_dydt_dy l2 dr2 '  &
                  // trim(chem_isos% name(chem_id(i1))) // ' ' // trim(chem_isos% name(chem_id(i2))),  &
                  - c * dr2, n% d_dydt_dy(i1, i2)

         if (i3 <= 0) return

         ! update the Jacobian for d_dydt(i1)_dy(i3)
         n% d_dydt_dy(i1, i3) = n% d_dydt_dy(i1, i3)  - c * dr3

         if (chem_id(i1) == ini56 .and. show_jac .and. c * dr3 /= 0) &
               write(*,1) trim(reaction_Name(n% g% reaction_id(i))) // ' d_dydt_dy l3 dr3 '  &
           // trim(chem_isos% name(chem_id(i1))) // ' ' // trim(chem_isos% name(chem_id(i3))),  &
           - c * dr3, n% d_dydt_dy(i1, i3)

         if (i4 <= 0) return

         !write(*,4) trim(reaction_Name(n% g% reaction_id(i))) // ' do_lhs_iso', i, i1, i4
         ! update the Jacobian for d_dydt(i1)_dy(i4)
         n% d_dydt_dy(i1, i4) = n% d_dydt_dy(i1, i4)  - c * dr4

         if (chem_id(i1) == ini56 .and. show_jac .and. c * dr4 /= 0) &
               write(*,1) trim(reaction_Name(n% g% reaction_id(i))) // ' d_dydt_dy l4 dr4 '  &
           // trim(chem_isos% name(chem_id(i1))) // ' ' // &
                  trim(chem_isos% name(chem_id(i4))),  &
                  - c * dr4, n% d_dydt_dy(i1, i4)

      end subroutine do_lhs_iso


      subroutine do_rhs_iso( &
            n, dydt, i, c, i1, rvs, i2, dr2, i3, dr3, i4, dr4, symbolic, just_dydt)
         type (Net_Info) :: n
         real(qp) :: dydt(:,:) ! (num_rvs, num_isos)
         integer, intent(in) :: i, i1, i2, i3, i4
         real(dp), intent(in) :: c, rvs(:), dr2, dr3, dr4
         logical, intent(in) :: symbolic, just_dydt

         ! i1, i2, i3, and 14 are isotope numbers
         ! dr2 = dr/dy(i2); dr3 = dr/dy(i3); dr4 = dr/dy(i4)
         ! c * r   = dydt(i1)
         ! c * dr2 = d_dydt(i1)_dy(i2)
         ! c * dr3 = d_dydt(i1)_dy(i3)
         ! c * dr4 = d_dydt(i1)_dy(i4)

         integer :: j
         integer, pointer :: chem_id(:)
         chem_id => n% g% chem_id

         include 'formats'

         if (symbolic) then
            n% d_dydt_dy(i1, i2) = 1
            if (i3 <= 0) return
            n% d_dydt_dy(i1, i3) = 1
            if (i4 <= 0) return
            n% d_dydt_dy(i1, i4) = 1
            return
         end if

         ! update the dydt terms for i1
         do j=1,num_rvs
            dydt(j,i1) = dydt(j,i1) + c * rvs(j)
         end do
         if (chem_id(i1) == io16 .and. show_dydt .and. &
                  abs(n% y(i1) - show_dydt_y) < 1d-20) &
               write(*,1) 'rhs ' // trim(reaction_Name(n% g% reaction_id(i))), &
               c * rvs(i_rate), dydt(i_rate,i1), &
               n% rate_screened(i), n% rate_raw(i), n% y(i1)
         if (chem_id(i1) == ihe4 .and. show_d_dydt_dRho .and. &
                  abs(n% y(i1) - show_dydt_y) < 1d-20) &
               write(*,1) trim(reaction_Name(n% g% reaction_id(i))) // ' d_dydt_dRho ' &
                  // trim(chem_isos% name(chem_id(i1))), &
                  c * rvs(i_rate_dRho), n% rate_screened(i), &
                  n% rate_screened_dRho(i), n% y(i1)
         if (chem_id(i1) == ihe4 .and. show_d_dydt_dT .and. &
                  abs(n% y(i1) - show_dydt_y) < 1d-20) &
               write(*,1) trim(reaction_Name(n% g% reaction_id(i))) // ' d_dydt_dT ' &
                  // trim(chem_isos% name(chem_id(i1))), &
                  c * rvs(i_rate_dT), n% rate_screened(i), &
                  n% rate_screened_dT(i), n% y(i1)

         if (just_dydt) return

         ! update the Jacobian for d_dydt(i1)_dy(i2)
         n% d_dydt_dy(i1, i2) = n% d_dydt_dy(i1, i2)  + c * dr2

         if (chem_id(i1) == ini56 .and. show_jac .and. c * dr2 /= 0) &
               write(*,1) trim(reaction_Name(n% g% reaction_id(i))) // ' d_dydt_dy r2 dr2 ' &
                  // trim(chem_isos% name(chem_id(i1))) // ' ' // &
                  trim(chem_isos% name(chem_id(i2))),  &
                  c * dr2, n% d_dydt_dy(i1, i2)

         if (i3 <= 0) return

         ! update the Jacobian for d_dydt(i1)_dy(i3)
         n% d_dydt_dy(i1, i3) = n% d_dydt_dy(i1, i3)  + c * dr3

         if (chem_id(i1) == ini56 .and. show_jac .and. c * dr3 /= 0) &
               write(*,1) trim(reaction_Name(n% g% reaction_id(i))) // ' d_dydt_dy r3 dr3 ' &
            // trim(chem_isos% name(chem_id(i1))) // ' ' // trim(chem_isos% name(chem_id(i3))),  &
           c * dr3, n% d_dydt_dy(i1, i3)

         if (i4 <= 0) return

         ! update the Jacobian for d_dydt(i1)_dy(i4)
         n% d_dydt_dy(i1, i4) = n% d_dydt_dy(i1, i4)  + c * dr4

         if (chem_id(i1) == ini56 .and. show_jac .and. c * dr4 /= 0) &
               write(*,1) trim(reaction_Name(n% g% reaction_id(i))) // ' d_dydt_dy r4 dr4 ' &
            // trim(chem_isos% name(chem_id(i1))) // ' ' // trim(chem_isos% name(chem_id(i4))),  &
           c * dr4, n% d_dydt_dy(i1, i4)

      end subroutine do_rhs_iso


      subroutine check_balance(n, i, lhs, rhs) ! check conservation of nucleons
         type (Net_Info) :: n
         integer, intent(in) :: i
         real(dp), intent(in) :: lhs, rhs
         if (lhs == rhs) return
         if (abs(lhs-rhs) < 1d-6) return
         write(*,'(2a)') 'non-conservation of nucleons in reaction ',  &
               reaction_Name(n% g% reaction_id(i))
         write(*,*) 'lhs aion sum', lhs
         write(*,*) 'rhs aion sum', rhs
         stop
      end subroutine check_balance



      end module net_derivs_support
