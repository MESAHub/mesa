! ***********************************************************************
!
!   Copyright (C) 2011-2019  Bill Paxton & The MESA Team

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

      program sample_net
      use utils_lib, only: mesa_error
      implicit none
      
      call test
      
      
      contains
      
      
      
      subroutine test
         use rates_def, only: rates_NACRE_if_available
         use chem_def, only: num_categories
         
         integer :: ierr, handle, which_rates_choice, species, &
            num_reactions, lwork
         integer, pointer :: which_rates(:), chem_id(:), net_iso(:)
         character (len=100) :: net_file
         character (len=64) :: mesa_dir
         
         mesa_dir = '../..'         
         net_file = 'basic.net'
         which_rates_choice = rates_NACRE_if_available

         ierr = 0
         call initialize(mesa_dir, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         call setup_net( &
            net_file, handle, which_rates, which_rates_choice, &
            species, chem_id, net_iso, num_reactions, lwork, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         call do1_net_eval( &
            handle, species, num_reactions, &
            chem_id, net_iso, lwork, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)         
         
      end subroutine test
      
      
      subroutine initialize(mesa_dir, ierr)
         use const_lib, only: const_init
         use math_lib
         use chem_lib, only: chem_init
         use rates_lib, only: rates_init, rates_warning_init
         use net_lib, only : net_init
         character (len=*), intent(in) :: mesa_dir
         integer, intent(out) :: ierr
         ierr = 0

         call math_init()
         
         call const_init(mesa_dir,ierr)     
         if (ierr /= 0) then
            write(*,*) 'const_init failed'
            call mesa_error(__FILE__,__LINE__)
         end if     
            
         call chem_init('isotopes.data', ierr)
         if (ierr /= 0) then
            write(*,*) 'chem_init failed'
            return
         end if


        call rates_init('reactions.list', '', 'rate_tables', .false., .false.,&
                     '', '', '',  ierr)
         if (ierr /= 0) then
            write(*,*) 'rates_init failed'
            return
         end if
         
         call rates_warning_init(.true., 10d0)
         
         call net_init(ierr)
         if (ierr /= 0) then
            write(*,*) 'net_init failed'
            return
         end if        
          
      end subroutine initialize
      
      
      subroutine setup_net( &
            net_file, handle, which_rates, which_rates_choice, &
            species, chem_id, net_iso, num_reactions, lwork, ierr)
         use net_lib
         use rates_def, only: rates_reaction_id_max
         
         character (len=*), intent(in) :: net_file
         integer, intent(in) :: which_rates_choice
         integer, pointer :: which_rates(:) ! will be allocated
         integer, pointer :: chem_id(:), net_iso(:) ! set, but not allocated
         integer, intent(out) :: handle, species, num_reactions, lwork, ierr
         
         ierr = 0
         handle = alloc_net_handle(ierr)
         if (ierr /= 0) then
            write(*,*) 'alloc_net_handle failed'
            return
         end if
         
         call net_start_def(handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'net_start_def failed'
            return
         end if
         
         write(*,*) 'load ' // trim(net_file)
         call read_net_file(net_file, handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'read_net_file failed ', trim(net_file)
            return
         end if
         
         call net_finish_def(handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'net_finish_def failed'
            return
         end if
      
         allocate(which_rates(rates_reaction_id_max))
         which_rates(:) = which_rates_choice

         call net_set_which_rates(handle, which_rates, ierr)
         if (ierr /= 0) then
            write(*,*) 'net_set_which_rate_f17pg failed'
            return
         end if
         
         call net_setup_tables(handle, '', ierr)
         if (ierr /= 0) then
            write(*,*) 'net_setup_tables failed'
            return
         end if
         
         species = net_num_isos(handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in net_num_isos'
            return
         end if
         
         call get_chem_id_table_ptr(handle, chem_id, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in get_chem_id_table_ptr'
            return
         end if
         
         call get_net_iso_table_ptr(handle, net_iso, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in get_net_iso_table_ptr'
            return
         end if
         
         num_reactions = net_num_reactions(handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in net_num_reactions'
            return
         end if

         lwork = net_work_size(handle, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in net_work_size'
            return
         end if
         
      end subroutine setup_net
      
      
      subroutine do1_net_eval( &
            handle, species, num_reactions, chem_id, net_iso, lwork, ierr)
            
         use rates_def
         use chem_def
         use net_def
         use net_lib
         use chem_lib
         
         integer, intent(in) :: handle, species, num_reactions, &
            chem_id(:), net_iso(:), lwork
         integer, intent(out) :: ierr
         
         integer :: screening_mode
         real(dp) :: xa(species), T, logT, Rho, logRho, eta, d_eta_dlnT, d_eta_dlnRho, &
            d_eps_nuc_dx(species), dabar_dx(species), dzbar_dx(species), dmc_dx(species), &
            weak_rate_factor, xh, xhe, z, abar, zbar, z2bar, z53bar, ye, mass_correction, xsum, &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, &
            eps_nuc_categories(num_categories), eps_neu_total, &
            dxdt(species), d_dxdt_dRho(species), d_dxdt_dT(species), &
            d_dxdt_dx(species,species)
         real(dp), target :: work_ary(lwork), rate_factors_ary(num_reactions)
         real(dp), pointer, dimension(:) :: work, rate_factors
         logical :: skip_jacobian
         type (Net_Info), target :: netinfo_target
         type (Net_Info), pointer :: netinfo
         
         include "formats"
         
         ierr = 0
         work => work_ary
         rate_factors => rate_factors_ary
         netinfo => netinfo_target
         
         ! set mass fractions -- must add to 1.0
         xa = 0
         xa(net_iso(ih1)) = 7.5876644280605066d-01
         xa(net_iso(ihe4)) = 2.3952230737160904d-01
         xa(net_iso(img24)) = 1 - sum(xa(:))
         
         call composition_info( &
            species, chem_id, xa, xh, xhe, z, abar, zbar, z2bar, z53bar, &
            ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)
            
         logT = 8
         T = 1d8
         logRho = 6
         Rho = 1d6
         
         eta = 0
         rate_factors(:) = 1
         weak_rate_factor = 1

         screening_mode = extended_screening
         
         skip_jacobian = .false.
         
         call net_get( &
            handle, skip_jacobian, netinfo, species, num_reactions, &
            xa, T, logT, Rho, logRho, & 
            abar, zbar, z2bar, ye, eta, d_eta_dlnT, d_eta_dlnRho, &
            rate_factors, weak_rate_factor, & 
            std_reaction_Qs, std_reaction_neuQs, &
            eps_nuc, d_eps_nuc_dRho, d_eps_nuc_dT, d_eps_nuc_dx, & 
            dxdt, d_dxdt_dRho, d_dxdt_dT, d_dxdt_dx, & 
            screening_mode, &     
            eps_nuc_categories, eps_neu_total, & 
            lwork, work, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in net_get'
            return
         end if
         
         write(*,1) 'logT', logT
         write(*,1) 'logRho', logRho
         write(*,1) 'eps_nuc', eps_nuc
         write(*,*)

      end subroutine do1_net_eval
      
      
      
      end program sample_net




