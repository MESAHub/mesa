! ***********************************************************************
!
!   Copyright (C) 2011  Bill Paxton
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
!
! ***********************************************************************

      module ionization_def
      
      use const_def, only: dp, use_mesa_temp_cache
      
      implicit none
      
      ! ionization results

      integer, parameter :: ion_ilogPgas = 1 ! log10 Pgas
      ! log10 of partial pressures
      integer, parameter :: ion_ilogpp_H = ion_ilogPgas + 1
      integer, parameter :: ion_ilogpp_He = ion_ilogpp_H + 1
      integer, parameter :: ion_ilogpp_C = ion_ilogpp_He + 1
      integer, parameter :: ion_ilogpp_N = ion_ilogpp_C + 1
      integer, parameter :: ion_ilogpp_O = ion_ilogpp_N + 1
      integer, parameter :: ion_ilogpp_Ne = ion_ilogpp_O + 1
      integer, parameter :: ion_ilogpp_Mg = ion_ilogpp_Ne + 1
      integer, parameter :: ion_ilogpp_Si = ion_ilogpp_Mg + 1
      integer, parameter :: ion_ilogpp_Fe = ion_ilogpp_Si + 1
      ! charge
      integer, parameter :: ion_iZ_H = ion_ilogpp_Fe + 1
      integer, parameter :: ion_iZ_He = ion_iZ_H + 1
      integer, parameter :: ion_iZ_C = ion_iZ_He + 1
      integer, parameter :: ion_iZ_N = ion_iZ_C + 1
      integer, parameter :: ion_iZ_O = ion_iZ_N + 1
      integer, parameter :: ion_iZ_Ne = ion_iZ_O + 1
      integer, parameter :: ion_iZ_Mg = ion_iZ_Ne + 1
      integer, parameter :: ion_iZ_Si = ion_iZ_Mg + 1
      integer, parameter :: ion_iZ_Fe = ion_iZ_Si + 1
      ! fraction neutral
      integer, parameter :: ion_ifneut_H = ion_iZ_Fe + 1
      integer, parameter :: ion_ifneut_He = ion_ifneut_H + 1
      integer, parameter :: ion_ifneut_C = ion_ifneut_He + 1
      integer, parameter :: ion_ifneut_N = ion_ifneut_C + 1
      integer, parameter :: ion_ifneut_O = ion_ifneut_N + 1
      integer, parameter :: ion_ifneut_Ne = ion_ifneut_O + 1
      integer, parameter :: ion_ifneut_Mg = ion_ifneut_Ne + 1
      integer, parameter :: ion_ifneut_Si = ion_ifneut_Mg + 1
      integer, parameter :: ion_ifneut_Fe = ion_ifneut_Si + 1

      integer, parameter :: num_ion_vals = ion_ifneut_Fe
      
      character (len=20) :: ion_result_names(num_ion_vals)




      
      ! based on scheme for eos tables

      
      integer, parameter :: num_ion_Zs = 5
      real(dp), parameter :: ion_Zs(num_ion_Zs) = (/ 0.00d0, 0.02d0, 0.04d0, 0.20d0, 1.00d0 /)
      integer, parameter :: num_ion_Xs_for_Z(num_ion_Zs) = (/ 6, 5, 5, 5, 1 /)
      
      integer, parameter :: num_ion_Xs = 6
      real(dp), parameter :: ion_Xs(num_ion_Xs) = &
            (/ 0.0d0, 0.2d0, 0.4d0, 0.6d0, 0.8d0, 1.0d0 /)

      integer, parameter :: sz_per_ion_point = 4 ! for bicubic spline interpolation

      real(dp) :: ion_logQ_min ! logQ = logRho - 2*logT + 12
      real(dp) :: ion_logQ_max
      real(dp) :: ion_del_logQ ! spacing for the logQs
      integer :: ion_num_logQs
      real(dp) :: ion_logT_min
      real(dp) :: ion_logT_max
      real(dp) :: ion_del_logT ! spacing for the logTs
      integer :: ion_num_logTs
      real(dp), dimension(:), pointer :: ion_logQs, ion_logTs
      real(dp), dimension(:), pointer :: ion_tbl1
      real(dp), dimension(:,:,:,:,:,:), pointer :: ion_tbl
         ! dimension(sz_per_ion_point, num_ion_vals, num_logQs, num_logTs, num_ion_Xs, num_ion_Zs)
      integer :: ion_version

      integer, parameter :: min_version = 49

      character(len=32) :: ion_file_prefix, ion_Z1_suffix
      character(len=1000) :: ionization_cache_dir, ionization_temp_cache_dir
      
      logical :: use_cache_for_ion = .true.
      logical :: ion_root_is_initialized = .false.
      logical :: ion_is_initialized = .false.




      integer, parameter :: table_version = 1
      
      type Ionization_Info
         integer :: num_log_ne, num_logT
         real(dp), pointer :: log_ne(:), logT(:)
         real(dp), pointer :: f1(:)
         logical :: have_interpolation_info
         integer :: ilinx, iliny
      end type Ionization_Info
      
      
      type (Ionization_Info), target :: fe_he_info
      type (Ionization_Info), pointer :: fe_he_ptr


      logical :: table_is_initialized = .false.






      contains

      
      subroutine ion_def_init(ionization_cache_dir_in)
         use const_def, only: mesa_data_dir, mesa_caches_dir, mesa_temp_caches_dir
         use utils_lib, only : mkdir
         character (len=*), intent(in) :: ionization_cache_dir_in
         
         ion_is_initialized = .false.
         use_cache_for_ion = .true.

         if (len_trim(ionization_cache_dir_in) > 0) then
            ionization_cache_dir = ionization_cache_dir_in
         else if (len_trim(mesa_caches_dir) > 0) then
            ionization_cache_dir = trim(mesa_caches_dir) // '/ionization_cache'
         else
            ionization_cache_dir = trim(mesa_data_dir) // '/ionization_data/cache'
         end if
         call mkdir(ionization_cache_dir)
         
         ionization_temp_cache_dir = trim(mesa_temp_caches_dir)//'/ionization_cache'
         if(use_mesa_temp_cache) call mkdir(ionization_temp_cache_dir)

         ion_file_prefix = 'ion'
         ion_Z1_suffix = '_CO_1'
         ion_result_names(ion_ilogPgas) = 'logPgas'
         ion_result_names(ion_ilogpp_H) = 'logpp_H'
         ion_result_names(ion_ilogpp_He) = 'logpp_He'
         ion_result_names(ion_ilogpp_C) = 'logpp_C'
         ion_result_names(ion_ilogpp_N) = 'logpp_N'
         ion_result_names(ion_ilogpp_O) = 'logpp_O'
         ion_result_names(ion_ilogpp_Ne) = 'logpp_Ne'
         ion_result_names(ion_ilogpp_Mg) = 'logpp_Mg'
         ion_result_names(ion_ilogpp_Si) = 'logpp_Si'
         ion_result_names(ion_ilogpp_Fe) = 'logpp_Fe'
         ion_result_names(ion_iZ_H) = 'Z_H'
         ion_result_names(ion_iZ_He) = 'Z_He'
         ion_result_names(ion_iZ_C) = 'Z_C'
         ion_result_names(ion_iZ_N) = 'Z_N'
         ion_result_names(ion_iZ_O) = 'Z_O'
         ion_result_names(ion_iZ_Ne) = 'Z_Ne'
         ion_result_names(ion_iZ_Mg) = 'Z_Mg'
         ion_result_names(ion_iZ_Si) = 'Z_Si'
         ion_result_names(ion_iZ_Fe) = 'Z_Fe'
         ion_result_names(ion_ifneut_H) = 'fneut_H'
         ion_result_names(ion_ifneut_He) = 'fneut_He'
         ion_result_names(ion_ifneut_C) = 'fneut_C'
         ion_result_names(ion_ifneut_N) = 'fneut_N'
         ion_result_names(ion_ifneut_O) = 'fneut_O'
         ion_result_names(ion_ifneut_Ne) = 'fneut_Ne'
         ion_result_names(ion_ifneut_Mg) = 'fneut_Mg'
         ion_result_names(ion_ifneut_Si) = 'fneut_Si'
         ion_result_names(ion_ifneut_Fe) = 'fneut_Fe'
      end subroutine ion_def_init


      end module ionization_def
      
