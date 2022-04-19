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

! this program shows (a) how to setup a network in mesa, and (b) return the
! instantaneous energy generation rate and mass fraction rates of change.
! this program does not do a time integration of a reaction network;
! for that you want to see $MESA_DIR/net/test/one_zone_burn.f90.

      
      call test
      
      contains
      
      
      
      subroutine test
         use chem_def, only: num_categories
         
         integer :: ierr, handle, species, &
            num_reactions, lwork
         integer, pointer :: chem_id(:), net_iso(:)
         character (len=100) :: net_file
         character (len=64) :: mesa_dir


! explicitly set my_mesa_dir to your $MESA_DIR, or use a blank string, in which case your $MESA_DIR is automagically used         

         mesa_dir = '../..'         
!         mesa_dir = ''         

! choose the network to use

         net_file = 'approx21.net'


! initialize
         ierr = 0
         call initialize(mesa_dir, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

! set up the network         
         call setup_net( &
            net_file, handle, &
            species, chem_id, net_iso, num_reactions, lwork, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)


! call the burner         
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
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
            
         call chem_init('isotopes.data', ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)


        call rates_init('reactions.list', '', 'rate_tables', .false., .false.,&
                     '', '', '',  ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         call rates_warning_init(.true., 10d0)
         
         call net_init(ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
          
      end subroutine initialize
      
      
      subroutine setup_net( &
            net_file, handle, &
            species, chem_id, net_iso, num_reactions, lwork, ierr)
         use net_lib
         use rates_def, only: rates_reaction_id_max
         
         character (len=*), intent(in) :: net_file
         integer, pointer :: chem_id(:), net_iso(:) ! set, but not allocated
         integer, intent(out) :: handle, species, num_reactions, lwork, ierr
         
         ierr = 0
         handle = alloc_net_handle(ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         call net_start_def(handle, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         write(*,*) 'load ' // trim(net_file)
         call read_net_file(net_file, handle, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         call net_finish_def(handle, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
               
         call net_setup_tables(handle, '', ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         species = net_num_isos(handle, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         call get_chem_id_table_ptr(handle, chem_id, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         call get_net_iso_table_ptr(handle, net_iso, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
         num_reactions = net_num_reactions(handle, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)

         lwork = net_work_size(handle, ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         
      end subroutine setup_net
      
      
      subroutine do1_net_eval( &
            handle, species, num_reactions, chem_id, net_iso, lwork, ierr)
            
         use rates_def
         use chem_def
         use net_def
         use net_lib
         use chem_lib
      
! declare the pass   
         integer, intent(in) :: handle, species, num_reactions, &
            chem_id(:), net_iso(:), lwork
         integer, intent(out) :: ierr
         

! locals
         integer :: screening_mode, i
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
         character (len=80) :: string
         
         include "formats"


! popular format statements
21    format(1x,t2,a,t16,a,t31,a,t46,a,t59,a,t74,a,t89,a)
22    format(1x,t2,a,1p7e15.6)
23    format(1x,t2,a7,1pe14.6,t24,a7,1pe14.6,t46,a7,1pe14.6,t68,a7,1pe14.6)
24    format(1x,t2,a12,1pe14.6,t30,a12,1pe14.6,t60,a12,1pe14.6,t90,a12,1pe14.6,t120,a12,1pe14.6)
         

! set some pointers and options
         ierr = 0
         work => work_ary
         rate_factors => rate_factors_ary
         netinfo => netinfo_target

         eta = 0
         rate_factors(:) = 1
         weak_rate_factor = 1
         screening_mode = extended_screening
         skip_jacobian = .false.




! main loop, keep returning here
100   xa(:) = 0

      write(6,*)  
      write(6,*) 'give the temperature, density, and mass fractions (h1, he4, c12, n14, o16, ne20, mg24) =>'
      write(6,*) 'hit return for T = 1e9 K, Rho = 1e4 g/cc, x(c12) = 1 ; enter -1 to stop'
      write(6,*)  
      read(5,'(a)') string

! stop
      if (string(1:2) .eq. '-1') then
       call mesa_error(__FILE__,__LINE__,'normal termination')

! read the conditions
      else
       if (string(1:6) .ne. '      ') then
        read(string,*) T,Rho, xa(net_iso(ih1)), xa(net_iso(ihe4)), xa(net_iso(ic12)), & 
                       xa(net_iso(in14)), xa(net_iso(io16)), xa(net_iso(ine20)), xa(net_iso(img24))
! or set some defaults
       else
        T = 1.0d9 ; Rho = 1.0d4 ; xa(net_iso(ic12)) = 1.0d0
       end if
      end if

      logT   = log10(T)
      logRho = log10(Rho)



! get some composition variables
         call composition_info( &
            species, chem_id, xa, xh, xhe, z, abar, zbar, z2bar, z53bar, &
            ye, mass_correction, xsum, dabar_dx, dzbar_dx, dmc_dx)

         
! this is the instantaneous eps_nuc only

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
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__)
         

! say the initial conditions
         write(6,23) 'T     =',T,       'Rho   =',Rho,      'abar  =',abar,   'zbar  =',zbar
         write(6,23) 'h1    =',xa(net_iso(ih1)),  'he4   =',xa(net_iso(ihe4)),  'c12   =',xa(net_iso(ic12)), 'n14   =',xa(net_iso(in14))
         write(6,23) 'o16   =',xa(net_iso(io16)), 'ne20  =',xa(net_iso(ine20)), 'mg24  =',xa(net_iso(img24))
 

! write out the mass fraction changes
         write(6,'(A)')
         write(6,24) 'd(h1)/dt   =',dxdt(net_iso(ih1)), 'd(he4)/dt  =',dxdt(net_iso(ihe4)), 'd(c12)/dt  =',dxdt(net_iso(ic12)), 'd(n14)/dt  =',dxdt(net_iso(in14))
         write(6,24) 'd(o16)/dt  =',dxdt(net_iso(io16)), 'd(ne20)/dt =',dxdt(net_iso(ine20)), 'd(mg24)/dt =',dxdt(net_iso(img24))


! check non-conservation
         xsum =   dxdt(net_iso(ih1))  + dxdt(net_iso(ihe4))  + dxdt(net_iso(ic12)) + dxdt(net_iso(in14)) &
                + dxdt(net_iso(io16)) + dxdt(net_iso(ine20)) + dxdt(net_iso(img24))
         write(6,24) '1 - sum    =',xsum


! instantaneous net energy generation rate
         write(6,24) 'eps_nuc    =',eps_nuc, 'erg/g/sec'


! back for more
         goto 100


      end subroutine do1_net_eval
      
      
      
      end program sample_net




