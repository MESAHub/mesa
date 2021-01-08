! ***********************************************************************
!
!   Copyright (C) 2009-2019  Bill Paxton, Frank Timmes & The MESA Team
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

      module test_net_do_one
      use chem_def
      use chem_lib
      use net_def
      use net_lib
      use const_def
      use rates_def
      use utils_lib, only: mesa_error
      
      implicit none
         
      logical, parameter :: extended_set = .false.
      logical, parameter :: sorted = .true.

      logical :: qt
      character (len=64) :: net_file
      integer :: handle
      integer :: species
      real(dp) :: z, abar, zbar, z2bar, z53bar, ye, &
         eta, d_eta_dlnT, d_eta_dlnRho, eps_neu_total
      integer :: screening_mode
      real(dp) :: test_logT, test_logRho     
      logical :: reuse_rate_raw, reuse_rate_screened = .false.
      integer, pointer :: reaction_id(:)
      real(dp), dimension(:), pointer ::  &
         xin, xin_copy, d_eps_nuc_dx, dxdt, d_dxdt_dRho, d_dxdt_dT
      real(dp), pointer :: d_dxdt_dx(:,:)  
      type (Net_Info), target :: net_info_target
      type (Net_Info), pointer :: n
      

      contains
      
      
      subroutine do1_net(symbolic)
         use chem_lib, only:composition_info
         logical, intent(in) :: symbolic

         real(dp) :: logRho, logT, Rho, T, sum, mass_correction, weak_rate_factor, &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, xh, xhe, rate_limit

         integer :: info, i, j, k, lwork, chem_id(species), num_reactions
         real(dp), dimension(species) :: dabar_dx, dzbar_dx, dmc_dx
         real(dp), dimension(:), pointer :: eps_nuc_categories
         real(dp), pointer :: work(:), rate_factors(:)
         type (Net_General_Info), pointer  :: g
         logical :: skip_jacobian
         
         info = 0
         call get_net_ptr(handle, g, info)
         if (info /= 0) call mesa_error(__FILE__,__LINE__)

         num_reactions = g% num_reactions         

         logRho = test_logRho
         logT   = test_logT

         if (.not. qt) write(*,*)
         
         lwork = net_work_size(handle, info) 
         if (info /= 0) call mesa_error(__FILE__,__LINE__)
         
         allocate(work(lwork),  &
            rate_factors(num_reactions), &
            eps_nuc_categories(num_categories), &
            stat=info)
         if (info /= 0) call mesa_error(__FILE__,__LINE__)
         
         call get_chem_id_table(handle, species, chem_id, info)
         if (info /= 0) call mesa_error(__FILE__,__LINE__)

         call composition_info( &
            species, chem_id, xin, xh, xhe, z, abar, zbar, z2bar, z53bar, ye, &
            mass_correction, sum, dabar_dx, dzbar_dx, dmc_dx)

         Rho = exp10(logRho)
         T   = exp10(logT)
         
         if (net_file == '19_to_ni56.net') then
            logT = 9D+00
            logRho = 8D+00
            eta = 3D+00
         end if
         
         if (net_file == 'approx21_cr60_plus_co56.net') then
            logT =    4.6233007922659333D+00
            logRho =   -1.0746410107891649D+01
            eta =   -2.2590260158215202D+01
         end if

         rate_factors(:) = 1
         weak_rate_factor = 1
         rate_limit = 0d0
         skip_jacobian = .false.
         d_eta_dlnT = 0d0
         d_eta_dlnRho = 0d0
         
         if (symbolic) then
            call net_get_symbolic_d_dxdt_dx(handle, n, species, num_reactions,  &
                  xin, T, logT, Rho, logRho,  &
                  abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
                  rate_factors, weak_rate_factor, &
                  std_reaction_Qs, std_reaction_neuQs, &
                  eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
                  dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
                  screening_mode,      &
                  eps_nuc_categories, eps_neu_total, &
                  lwork, work, info)
         else
            call net_get(handle, skip_jacobian, n, species, num_reactions,  &
                  xin, T, logT, Rho, logRho,  &
                  abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
                  rate_factors, weak_rate_factor, &
                  std_reaction_Qs, std_reaction_neuQs, reuse_rate_raw, reuse_rate_screened, &
                  eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
                  dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
                  screening_mode,    &
                  eps_nuc_categories, eps_neu_total, &
                  lwork, work, info)
         end if
         if (info /= 0) then
            write(*, *) 'bad return from net_get'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         if (symbolic .and..not. qt) then
            write(*,*) 'nonzero d_dxdt_dx entries'
            k = 0
            do j=1,species
               do i=1,species
                  if (d_dxdt_dx(i,j) /= 0) then
                     k = k + 1
                     write(*,'(a50,2i5)')  &
                           trim(chem_isos% name(chem_id(i))) //  &
                           ' ' // trim(chem_isos% name(chem_id(j))), i, j
                  end if
               end do
            end do
            write(*,*)
            write(*,'(a50,i5)') 'num non zeros', k
            write(*,*)
         else if (.not. qt) then
            call show_results( &
                  g, lwork, work, logT, logRho, species, num_reactions, xin,  &
                  eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
                  dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
                  eps_nuc_categories, extended_set, sorted)         
         end if

         deallocate(work, rate_factors, eps_nuc_categories)
 
         return

         write(*, *)
 1       format(a40, 6x, e25.10)
 2       format(a40, a6, e25.10)
         write(*, 1) 'abar', abar
         do i=1, species
            write(*, 2) 'dabar_dx', trim(chem_isos% name(chem_id(i))), dabar_dx(i)
         end do
         write(*, *)
         write(*, 1) 'zbar', zbar
         do i=1, species
            write(*, 2) 'dzbar_dx', trim(chem_isos% name(chem_id(i))), dzbar_dx(i)
         end do
         write(*, *)

      end subroutine do1_net



      subroutine show_results( &
            g, lwork, work, logT, logRho, species, num_reactions, xin,  &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx,  &
            eps_nuc_categories,  &
            extended_set, sorted)
         type (Net_General_Info), pointer  :: g
         integer, intent(in) :: lwork
         real(dp), pointer :: work(:)
         real(dp), intent(in) :: logT, logRho
         integer, intent(in) :: species, num_reactions
         real(dp), intent(in) :: xin(species)
         real(dp), intent(in) :: eps_nuc
         real(dp), intent(in) :: d_eps_nuc_dT
         real(dp), intent(in) :: d_eps_nuc_dRho
         real(dp), intent(in) :: d_eps_nuc_dx(species) 
         real(dp), intent(in) :: dxdt(species) 
         real(dp), intent(in) :: d_dxdt_dRho(species) 
         real(dp), intent(in) :: d_dxdt_dT(species) 
         real(dp), intent(in) :: d_dxdt_dx(species, species) 
         real(dp), intent(in), dimension(num_categories) :: eps_nuc_categories 
         logical, intent(in) :: extended_set
         logical, intent(in) :: sorted
         
         integer :: ierr, j
         real(dp), dimension(:), pointer :: &
            rate_raw, rate_raw_dT, rate_raw_dRho, &
            rate_screened, rate_screened_dT, rate_screened_dRho
            
         include 'formats'
         
         call get_net_rate_ptrs(g% handle, &
            rate_screened, rate_screened_dT, rate_screened_dRho, &
            rate_raw, rate_raw_dT, rate_raw_dRho, lwork, work, &
            ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in get_net_rate_ptrs'
            call mesa_error(__FILE__,__LINE__)
         end if
                  
         if (net_file == 'approx21_cr60_plus_co56.net') then
            write(*, *)
            write(*, '(a40, f20.9)') 'log temp', logT
            write(*, '(a40, f20.9)') 'log rho', logRho
            write(*, *)
            j = irco56ec_to_fe56
            write(*,1) 'rate_raw rco56ec_to_fe56', &
               rate_raw(g% net_reaction(j))
            j = irni56ec_to_co56
            write(*,1) 'rate_raw rni56ec_to_co56', &
               rate_raw(g% net_reaction(j))
            return
         end if
         
         write(*, *)
         call show_summary_results(logT, logRho, &
               eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx) 
         
         if (extended_set) then
            write(*, *)
            call show_all_rates( &
                g, num_reactions, rate_raw, rate_raw_dT, rate_raw_dRho, &
                rate_screened, rate_screened_dT, rate_screened_dRho, &
                logT, logRho, sorted)
         end if
         
         write(*, *)
         call show_by_category( &
                  g, num_reactions, eps_nuc_categories,  &
                  eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
                  sorted)
         
         if (.not. extended_set) return
         
         write(*, *)
         call show_dx_dt(g, species, xin, dxdt, sorted)
         
         write(*, *)
         call show_d_eps_nuc_dx(g, species, xin, d_eps_nuc_dx, sorted)

         write(*, *)

      end subroutine show_results

      
      subroutine show_summary_results(logT, logRho,  &
               eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx) 
         real(dp), intent(in) :: logT, logRho
         real(dp), intent(in) :: eps_nuc
         real(dp), intent(in) :: d_eps_nuc_dT
         real(dp), intent(in) :: d_eps_nuc_dRho
         real(dp), intent(in) :: d_eps_nuc_dx(species) 

         real(dp) :: T, Rho, eps, d_eps_dt, d_eps_dd
         T = exp10(logT); Rho = exp10(logRho)

         write(*, *)
         write(*, '(a40, f20.9)') 'log temp', logT
         write(*, '(a40, f20.9)') 'log rho', logRho
         eps = eps_nuc
         d_eps_dt = d_eps_nuc_dT
         d_eps_dd = d_eps_nuc_dRho
         write(*, *)
         write(*, '(a40, f20.9)') 'log(eps_nuc)', safe_log10(eps_nuc)
         write(*, *)
         write(*, '(a40, e20.9)') 'eps_nuc', eps_nuc
         write(*, *)
         write(*, '(a40, f20.9)') 'd_lneps_dlnT', d_eps_dt * T / eps
         write(*, '(a40, f20.9)') 'd_lneps_dlnRho', d_eps_dd * Rho / eps
      
      end subroutine show_summary_results

      
      subroutine show_all_rates( &
            g, num_reactions, rate_raw, rate_raw_dT, rate_raw_dRho, &
            rate_screened, rate_screened_dT, rate_screened_dRho, logT, logRho, sorted)
         type (Net_General_Info), pointer  :: g
         integer, intent(in) :: num_reactions
         real(dp), intent(in) :: logT, logRho
         real(dp), dimension(num_reactions), intent(in) :: &
            rate_raw, rate_raw_dT, rate_raw_dRho, &
            rate_screened, rate_screened_dT, rate_screened_dRho
         logical, intent(in) :: sorted
         
         real(dp), dimension(num_reactions) :: rfact
         integer :: i
         real(dp) :: T, Rho
         T = exp10(logT); Rho = exp10(logRho)

         write(*, *)
         write(*, *) 'summary of log raw rates'
         write(*, *)
         call show_log_rates(g, rate_raw, T, Rho, sorted)
         write(*, *)
         write(*, *) 'summary of screening factors'
         write(*, *)
         do i=1,num_reactions
            if (rate_raw(i) > 1d-50) then
               rfact(i) = rate_screened(i) / rate_raw(i)
            else
               rfact(i) = 1
            end if
         end do
         call show_rates(g, rfact, T, Rho, sorted)
         write(*, *)
         write(*, *) 'summary of log screened rates (reactions/gm/sec)'
         write(*, *)
         do i=1,num_reactions
            rfact(i) = rate_screened(i)
         end do
         call show_log_rates(g, rfact, T, Rho, sorted)
         write(*, *)

      end subroutine show_all_rates

      subroutine show_by_category( &
               g, num_reactions, eps_nuc_categories, &
               eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx,  &
               sorted)
         type (Net_General_Info), pointer  :: g
         integer, intent(in) :: num_reactions
         real(dp), intent(in), dimension(num_categories) :: eps_nuc_categories
         real(dp), intent(in) ::  &
               eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx(species)
         logical, intent(in) :: sorted

         real(dp) :: mx
         integer :: k, j, jmx
         logical :: flgs(rates_reaction_id_max)
   
         integer :: info
         
         write(*, *)
         write(*, *) 'energy generation by category'
         write(*, *)
         write(*, '(a40, 3x, a20)') 'category', 'log rate (erg/g/sec)'
         write(*, *)
         flgs = .false.
         do k = 1, num_categories
            if (.not. sorted) then
               jmx = k
            else
               mx = -99d99; jmx = -1
               do j = 1, num_categories
                  if ((.not. flgs(j)) .and. eps_nuc_categories(j) > mx) then
                     mx = eps_nuc_categories(j); jmx = j
                  end if
               end do
               if (jmx <= 0) exit
               if (mx < 1) exit ! FOR TEST OUTPUT
               flgs(jmx) = .true.
            end if
            write(*, '(a40, 2x, f15.6, e15.6)')  &
                  trim(category_name(jmx)), safe_log10(eps_nuc_categories(jmx)), &
                  eps_nuc_categories(jmx)     
         end do
         write(*, *)
         !write(*, '(a40, 2x, f15.6, e15.6)')  &
         !         'log10(-photodisintegration)', safe_log10(-eps_nuc_categories(iphoto)), &
         !         -eps_nuc_categories(iphoto)
         
         write(*, *)
         
      end subroutine show_by_category


      subroutine show_rates(g, rts, T, Rho, sorted)
         type (Net_General_Info), pointer  :: g
         real(dp), intent(in) :: rts(rates_reaction_id_max), T, Rho
         logical, intent(in) :: sorted
         
         logical :: flgs(rates_reaction_id_max)
         real(dp) :: mx
         integer :: k, j, jmx
         
         flgs = .false.
         
         do k = 1, g% num_reactions
            if (.not. sorted) then
               jmx = k; mx = rts(jmx)
            else
               mx = -99d99; jmx = -1
               do j = 1, g% num_reactions
                  if ((.not. flgs(j)) .and. rts(j) > mx) then
                     mx = rts(j); jmx = j
                  end if
               end do
               if (jmx <= 0) exit
               if (mx < 1d-60) exit
               flgs(jmx) = .true.
            end if
            if (mx == 1) cycle
            write(*, '(a40, e20.9, 2e17.6)') trim(reaction_name(reaction_id(jmx))), mx
         end do
         
      end subroutine show_rates


      subroutine show_log_rates(g, rts, T, Rho, sorted)
         type (Net_General_Info), pointer  :: g
         real(dp), intent(in) :: rts(rates_reaction_id_max), T, Rho
         logical, intent(in) :: sorted
         
         logical :: flgs(rates_reaction_id_max)
         real(dp) :: mx
         integer :: k, j, jmx
         
         flgs = .false.
         
         do k = 1, g% num_reactions
            if (.not. sorted) then
               jmx = k; mx = rts(jmx)
            else
               mx = -99d99; jmx = -1
               do j = 1, g% num_reactions
                  if ((.not. flgs(j)) .and. rts(j) > mx) then
                     mx = rts(j); jmx = j
                  end if
               end do
               if (jmx <= 0) exit
               if (mx < 1d-40) exit
               flgs(jmx) = .true.
            end if
            if (mx == 1) cycle
            write(*, '(a40, f20.9, 2e17.6)') trim(reaction_name(reaction_id(jmx))), safe_log10(mx)
         end do
         
      end subroutine show_log_rates
      
      
      subroutine show_dx_dt(g, species, xin, dxdt, sorted)
         type (Net_General_Info), pointer  :: g
         integer, intent(in) :: species
         real(dp), intent(in) :: xin(species)
         real(dp), intent(in) :: dxdt(species)
         logical, intent(in) :: sorted

         write(*, *)
         write(*, *) 'summary of isotope mass abundance changes'
         write(*, *)
         write(*, '(a40, 2(a17))') 'isotope', 'x initial', 'dx_dt   '
         call show_partials(g, species, xin, dxdt, .true., sorted)
         
      end subroutine show_dx_dt
      
      
      subroutine show_d_eps_nuc_dx(g, species, xin, d_eps_nuc_dx, sorted)
         type (Net_General_Info), pointer  :: g
         integer, intent(in) :: species
         real(dp), intent(in) :: xin(species)
         real(dp), intent(in) :: d_eps_nuc_dx(species)
         logical, intent(in) :: sorted

         write(*, *)
         write(*, *) 'summary of d_eps_nuc_dx'
         write(*, *)
         write(*, '(a40, a17)') 'isotope', 'd_eps_nuc_dx'
         call show_partials(g, species, xin, d_eps_nuc_dx, .false., sorted)
         
      end subroutine show_d_eps_nuc_dx
      
      
      subroutine show_partials(g, species, xin, derivs, initX_flag, sorted)
         type (Net_General_Info), pointer  :: g
         integer, intent(in) :: species
         real(dp), intent(in) :: xin(species)
         real(dp), intent(in) :: derivs(species)
         logical, intent(in) :: initX_flag, sorted

         real(dp) :: mx
         integer :: k, j, jmx
         integer, pointer :: chem_id(:)
         logical :: iflgs(species)
         chem_id => g% chem_id
         write(*, *)
         iflgs = .false.
         do k = 1, species
            if (.not. sorted) then
               jmx = k
            else
               mx = -99d99; jmx = -1
               do j = 1, species
                  if ((.not. iflgs(j)) .and. abs(derivs(j)) > mx) then
                     mx = abs(derivs(j)); jmx = j
                  end if
               end do
               if (jmx <= 0) exit
               if (mx < 1d-40) exit
            end if
            if (initX_flag) then
               write(*, '(a40, 2e17.6)') trim(chem_isos% name(chem_id(jmx))), xin(jmx), derivs(jmx)
            else
               write(*, '(a40, e25.14)') trim(chem_isos% name(chem_id(jmx))), derivs(jmx)
            end if
            iflgs(jmx) = .true.
         end do
         
      end subroutine show_partials

      
      end module test_net_do_one




