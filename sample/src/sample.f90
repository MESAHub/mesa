! ***********************************************************************
!
!   Copyright (C) 2011-2019  The MESA Team
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

      program sample_eos
      use eos_def
      use eos_lib
      use chem_def
      use chem_lib
      use const_lib
      use math_lib

      implicit none

      real(dp) :: X, Z, Y, abar, zbar, z2bar, z53bar, ye
      integer, parameter :: species = 7
      integer, parameter :: h1=1, he4=2, c12=3, n14=4, o16=5, ne20=6, mg24=7
      integer, pointer, dimension(:) :: net_iso, chem_id
      real(dp) :: xa(species)
      character (len=256) :: my_mesa_dir


      call Sample
      
      contains
      
      subroutine Sample
         
         integer :: handle
         real(dp) :: Rho, T, Pgas, log10Rho, log10T, Prad, energy, entropy
         real(dp) :: dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, d_dlnRho_const_T, d_dlnT_const_Rho
         real(dp), dimension(num_eos_basic_results) :: res, d_dlnd, d_dlnT, d_dabar, d_dzbar
         real(dp), dimension(num_eos_d_dxa_results, species) :: d_dxa
         integer :: ierr
         
         ierr = 0

         my_mesa_dir = '..' ! if empty string, uses environment variable MESA_DIR
         call const_init(my_mesa_dir,ierr)     
         if (ierr /= 0) then
            write(*,*) 'const_init failed'
            call mesa_error(__FILE__,__LINE__)
         end if        
         
         call math_init()
         
         call chem_init('isotopes.data', ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in chem_init'
            call mesa_error(__FILE__,__LINE__)
         end if

         ! allocate and initialize the eos tables
         call Setup_eos(handle)
         
         allocate(net_iso(num_chem_isos), chem_id(species), stat=ierr)
         if (ierr /= 0) call mesa_error(__FILE__,__LINE__,'allocate failed')
         X = 0.70_dp
         Z = 0.02_dp
         call Init_Composition
         
         Rho = 1.3519d2
         log10T = 8.15_dp
         T = exp10(log10T)
         
         ! get a set of results for given temperature and density
         call eosDT_get( &
               handle, &
               species, chem_id, net_iso, xa, &
               Rho, log10(Rho), T, log10T,  &
               res, d_dlnd, d_dlnT, &
               d_dxa, ierr)
         
        
 1       format(a20,3x,e20.12)

         Pgas = exp(res(i_lnPgas))

         ! the indices for the results are defined in eos_def.f
         write(*,*)
         write(*,1) 'temperature', T
         write(*,1) 'density', Rho
         write(*,1) 'logT', log10T
         write(*,1) 'logRho', log10(Rho)
         write(*,*)
         write(*,1) 'Z', Z
         write(*,1) 'X', X
         write(*,1) 'abar', abar
         write(*,1) 'zbar', zbar
         write(*,*)
         write(*,1) 'Pgas', Pgas
         write(*,1) 'logPgas', res(i_lnPgas)/ln10
         write(*,1) 'grad_ad', res(i_grad_ad)
         write(*,1) 'c_P', res(i_Cp)
         write(*,*)

         call eosPT_get(  &
               handle, &
               species, chem_id, net_iso, xa, &
               Pgas, log10(Pgas), T, log10(T),  &
               Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas,  &
               res, d_dlnd, d_dlnT, d_dxa, ierr)
      
         ! the indices for the results are defined in eos_def.f
         write(*,*)
         write(*,1) 'temperature', T
         write(*,1) 'Pgas', Pgas
         write(*,1) 'logT', log10(T)
         write(*,1) 'logPgas', res(i_lnPgas)/ln10
         write(*,*)
         write(*,1) 'Z', Z
         write(*,1) 'X', X
         write(*,1) 'abar', abar
         write(*,1) 'zbar', zbar
         write(*,*)
         write(*,1) 'density', Rho
         write(*,1) 'logRho', log10Rho
         write(*,1) 'grad_ad', res(i_grad_ad)
         write(*,1) 'c_P', res(i_Cp)
         write(*,*)

         ! deallocate the eos tables
         call Shutdown_eos(handle)
         
         deallocate(net_iso, chem_id)
         
         if (ierr /= 0) then
            write(*,*) 'bad result from eos_get'
            call mesa_error(__FILE__,__LINE__)
         end if

      end subroutine Sample
      

      subroutine Setup_eos(handle)
         ! allocate and load the eos tables
         use eos_def
         use eos_lib
         integer, intent(out) :: handle

         integer :: ierr
         logical, parameter :: use_cache = .true.

         call eos_init('', use_cache, ierr)
         if (ierr /= 0) then
            write(*,*) 'eos_init failed in Setup_eos'
            call mesa_error(__FILE__,__LINE__)
         end if
         
         write(*,*) 'loading eos tables'
         
         handle = alloc_eos_handle(ierr)
         if (ierr /= 0) then
            write(*,*) 'failed trying to allocate eos handle'
            call mesa_error(__FILE__,__LINE__)
         end if
      
      end subroutine Setup_eos
      
      
      subroutine Shutdown_eos(handle)
         use eos_def
         use eos_lib
         integer, intent(in) :: handle
         call free_eos_handle(handle)
         call eos_shutdown
      end subroutine Shutdown_eos


      subroutine Init_Composition
         use chem_lib

         real(dp), parameter :: Zfrac_C = 0.173312d0
         real(dp), parameter :: Zfrac_N = 0.053177d0
         real(dp), parameter :: Zfrac_O = 0.482398d0
         real(dp), parameter :: Zfrac_Ne = 0.098675d0
         
         real(dp) :: frac, dabar_dx(species), dzbar_dx(species),  &
               sumx, xh, xhe, xz, mass_correction, dmc_dx(species)
         
         net_iso(:) = 0
         
         chem_id(h1) = ih1; net_iso(ih1) = h1
         chem_id(he4) = ihe4; net_iso(ihe4) = he4
         chem_id(c12) = ic12; net_iso(ic12) = c12
         chem_id(n14) = in14; net_iso(in14) = n14
         chem_id(o16) = io16; net_iso(io16) = o16
         chem_id(ne20) = ine20; net_iso(ine20) = ne20
         chem_id(mg24) = img24; net_iso(img24) = mg24
         
         Y = 1 - (X + Z)
               
         xa(h1) = X
         xa(he4) = Y
         xa(c12) = Z * Zfrac_C
         xa(n14) = Z * Zfrac_N
         xa(o16) = Z * Zfrac_O
         xa(ne20) = Z * Zfrac_Ne
         xa(species) = 1 - sum(xa(1:species-1))
         
         call composition_info( &
               species, chem_id, xa, xh, xhe, xz, abar, zbar, z2bar, z53bar, ye,  &
               mass_correction, sumx, dabar_dx, dzbar_dx, dmc_dx)
         
      end subroutine Init_Composition


      end   

