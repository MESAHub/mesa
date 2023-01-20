! ***********************************************************************
!
! Copyright (C) 2020 The MESA Team
!
! MESA is free software; you can use it and/or modify
! it under the combined terms and restrictions of the MESA MANIFESTO
! and the GNU General Library Public License as published
! by the Free Software Foundation; either version 2 of the License,
! or (at your option) any later version.
!
! You should have received a copy of the MESA MANIFESTO along with
! this software; if not, it is available at the mesa website:
! http://mesa.sourceforge.net/
!
! MESA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
!
! You should have received a copy of the GNU Library General Public License
! along with this software; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

   module eos_ctrls_io

   use const_def
   use eos_def
   use math_lib
   use utils_lib, only: mesa_error

   implicit none

   public :: read_namelist, write_namelist, get_eos_controls, set_eos_controls
   private

   ! controls for HELM
   real(dp) :: Z_all_HELM ! all HELM for Z >= this unless use_FreeEOS
   real(dp) :: logT_all_HELM ! all HELM for lgT >= this
   real(dp) :: logT_low_all_HELM ! all HELM for lgT <= this
   real(dp) :: coulomb_temp_cut_HELM, coulomb_den_cut_HELM

   ! controls for OPAL_SCVH
   logical :: use_OPAL_SCVH
   real(dp) :: logT_low_all_SCVH ! SCVH for lgT >= this
   real(dp) :: logT_all_OPAL ! OPAL for lgT <= this
   real(dp) :: logRho1_OPAL_SCVH_limit ! don't use OPAL_SCVH for logRho > this
   real(dp) :: logRho2_OPAL_SCVH_limit ! full OPAL_SCVH okay for logRho < this
   real(dp) :: logRho_min_OPAL_SCVH_limit ! no OPAL/SCVH for logRho < this
   real(dp) :: logQ_max_OPAL_SCVH ! no OPAL/SCVH for logQ > this
   real(dp) :: logQ_min_OPAL_SCVH ! no OPAL/SCVH for logQ <= this.
   real(dp) :: Z_all_OPAL ! all OPAL for Z <= this

   ! controls for FreeEOS
   logical :: use_FreeEOS
   real(dp) :: logQ_max_FreeEOS_hi
   real(dp) :: logQ_max_FreeEOS_lo
   real(dp) :: logQ_min_FreeEOS_hi
   real(dp) :: logQ_min_FreeEOS_lo
   real(dp) :: logRho_min_FreeEOS_hi
   real(dp) :: logRho_min_FreeEOS_lo
   real(dp) :: logRho_max_FreeEOS_hi
   real(dp) :: logRho_max_FreeEOS_lo
   real(dp) :: logT_min_FreeEOS_hi
   real(dp) :: logT_min_FreeEOS_lo
   real(dp) :: logT_max_FreeEOS_hi
   real(dp) :: logT_max_FreeEOS_lo
   real(dp) :: logQ_cut_FreeEOS_lo_Z_max
   real(dp) :: logQ_cut_lo_Z_FreeEOS_hi
   real(dp) :: logQ_cut_lo_Z_FreeEOS_lo
   real(dp) :: logQ_cut_hi_Z_FreeEOS_hi
   real(dp) :: logQ_cut_hi_Z_FreeEOS_lo
   real(dp) :: logRho_cut_FreeEOS_hi
   real(dp) :: logRho_cut_FreeEOS_lo
   real(dp) :: logT_cut_FreeEOS_hi
   real(dp) :: logT_cut_FreeEOS_lo
   character (len=30) :: suffix_for_FreeEOS_Z(num_FreeEOS_Zs)
         
   ! controls for CMS
   logical :: use_CMS, CMS_use_fixed_composition
   integer :: CMS_fixed_composition_index
   real(dp) :: max_Z_for_any_CMS, max_Z_for_all_CMS ! set to -1 to disable CMS
   real(dp) :: logQ_max_for_any_CMS, logQ_max_for_all_CMS      ! for upper blend zone in logQ = logRho - 2*logT + 12
   real(dp) :: logQ_min_for_any_CMS, logQ_min_for_all_CMS      ! for lower blend zone in logQ
   real(dp) :: logRho_max_for_all_CMS, logRho_max_for_any_CMS  ! for upper blend zone in logRho
   real(dp) :: logRho_min_for_all_CMS, logRho_min_for_any_CMS  ! for lower blend zone in logRho
   real(dp) :: logT_max_for_all_CMS, logT_max_for_any_CMS      ! for upper blend zone in logT
   real(dp) :: logT_min_for_all_CMS, logT_min_for_any_CMS      ! for lower blend zone in logT      
   real(dp) :: logT_max_for_all_CMS_pure_He, logT_max_for_any_CMS_pure_He ! upper logT blend zone is different for pure He

   ! controls for PC
   logical :: use_PC
   real(dp) :: mass_fraction_limit_for_PC ! skip any species with abundance < this
   real(dp) :: logRho1_PC_limit ! okay for pure PC for logRho > this
   real(dp) :: logRho2_PC_limit ! don't use PC for logRho < this (>= 2.8)
   logical :: PC_use_Gamma_limit_instead_of_T
   real(dp) :: logT1_PC_limit ! okay for pure PC for logT < this (like logT_all_OPAL)
   real(dp) :: logT2_PC_limit ! don't use PC for logT > this (like logT_all_HELM)
   real(dp) :: log_Gamma_e_all_HELM ! HELM for log_Gamma_e <= this
   real(dp) :: Gamma_e_all_HELM ! 10**log_Gamma_e_all_HELM
   real(dp) :: log_Gamma_e_all_PC ! PC for log_Gamma_e >= this
   ! crystallization boundaries
   real(dp) :: PC_Gamma_start_crystal ! Begin releasing latent heat of crystallization
   real(dp) :: PC_Gamma_full_crystal ! Fully into the solid phase

   ! limits for Skye
   logical :: use_Skye
   logical :: Skye_use_ion_offsets
   real(dp) :: mass_fraction_limit_for_Skye
   real(dp) :: Skye_min_gamma_for_solid ! The minimum Gamma_i at which to use the solid free energy fit (below this, extrapolate).
   real(dp) :: Skye_max_gamma_for_liquid ! The maximum Gamma_i at which to use the liquid free energy fit (above this, extrapolate).
   character(len=128) :: Skye_solid_mixing_rule ! Currently support 'Ogata' or 'PC'

   logical :: use_simple_Skye_blends
   real(dp) :: logRho_min_for_any_Skye, logRho_min_for_all_Skye
   real(dp) :: logT_min_for_any_Skye, logT_min_for_all_Skye

   ! misc
   logical :: include_radiation, include_elec_pos
   logical :: eosDT_use_linear_interp_for_X
   logical :: eosDT_use_linear_interp_to_HELM
   character(len=128) :: eosDT_file_prefix
   logical :: okay_to_convert_ierr_to_skip
   real(dp) :: tiny_fuzz

   ! other eos
   logical :: use_other_eos_component, use_other_eos_results
   
   ! debugging
   logical :: dbg
   real(dp) :: logT_lo, logT_hi
   real(dp) :: logRho_lo, logRho_hi
   real(dp) :: X_lo, X_hi
   real(dp) :: Z_lo, Z_hi

   logical :: read_extra_eos_inlist1
   character (len=128) :: extra_eos_inlist1_name

   logical :: read_extra_eos_inlist2
   character (len=128) :: extra_eos_inlist2_name

   logical :: read_extra_eos_inlist3
   character (len=128) :: extra_eos_inlist3_name

   logical :: read_extra_eos_inlist4
   character (len=128) :: extra_eos_inlist4_name

   logical :: read_extra_eos_inlist5
   character (len=128) :: extra_eos_inlist5_name


   namelist /eos/ &
      use_FreeEOS, &
      
      ! controls for HELM
      Z_all_HELM, & ! all HELM for Z >= this unless use_FreeEOS
      logT_all_HELM, & ! all HELM for lgT >= this
      logT_low_all_HELM, & ! all HELM for lgT <= this
      coulomb_temp_cut_HELM, &
      coulomb_den_cut_HELM, &
      
      ! controls for OPAL_SCVH
      use_OPAL_SCVH, &
      logT_low_all_SCVH, & ! SCVH for lgT >= this
      logT_all_OPAL, & ! OPAL for lgT <= this
      logRho1_OPAL_SCVH_limit, & ! don't use OPAL_SCVH for logRho > this
      logRho2_OPAL_SCVH_limit, & ! full OPAL_SCVH okay for logRho < this
      logRho_min_OPAL_SCVH_limit, & ! no OPAL/SCVH for logRho < this
      logQ_max_OPAL_SCVH, & ! no OPAL/SCVH for logQ > this
      logQ_min_OPAL_SCVH, & ! no OPAL/SCVH for logQ <= this.
      Z_all_OPAL, & ! all OPAL for Z <= this
      
      ! controls for FreeEOS
      use_FreeEOS, &
      logQ_max_FreeEOS_hi, &
      logQ_max_FreeEOS_lo, &
      logQ_min_FreeEOS_hi, &
      logQ_min_FreeEOS_lo, &
      logRho_min_FreeEOS_hi, &
      logRho_min_FreeEOS_lo, &
      logRho_max_FreeEOS_hi, &
      logRho_max_FreeEOS_lo, &
      logT_min_FreeEOS_hi, &
      logT_min_FreeEOS_lo, &
      logT_max_FreeEOS_hi, &
      logT_max_FreeEOS_lo, &
      logQ_cut_FreeEOS_lo_Z_max, &
      logQ_cut_lo_Z_FreeEOS_hi, &
      logQ_cut_lo_Z_FreeEOS_lo, &
      logQ_cut_hi_Z_FreeEOS_hi, &
      logQ_cut_hi_Z_FreeEOS_lo, &
      logRho_cut_FreeEOS_hi, &
      logRho_cut_FreeEOS_lo, &
      logT_cut_FreeEOS_hi, &
      logT_cut_FreeEOS_lo, &
      suffix_for_FreeEOS_Z, &
      
      ! controls for CMS
      use_CMS, CMS_use_fixed_composition, &
      CMS_fixed_composition_index, &
      max_Z_for_any_CMS, &
      max_Z_for_all_CMS, & ! set to -1 to disable CMS
      logQ_max_for_any_CMS, &
      logQ_max_for_all_CMS, &      ! for upper blend zone in logQ = logRho - 2*logT + 12
      logQ_min_for_any_CMS, &
      logQ_min_for_all_CMS, &      ! for lower blend zone in logQ
      logRho_max_for_all_CMS, &
      logRho_max_for_any_CMS, &  ! for upper blend zone in logRho
      logRho_min_for_all_CMS, &
      logRho_min_for_any_CMS, &  ! for lower blend zone in logRho
      logT_max_for_all_CMS, &
      logT_max_for_any_CMS, &      ! for upper blend zone in logT
      logT_min_for_all_CMS, &
      logT_min_for_any_CMS, &      ! for lower blend zone in logT      
      logT_max_for_all_CMS_pure_He, &
      logT_max_for_any_CMS_pure_He, & ! upper logT blend zone is different for pure He
      
      ! controls for PC
      use_PC, &
      mass_fraction_limit_for_PC, & ! skip any species with abundance < this
      logRho1_PC_limit, & ! okay for pure PC for logRho > this
      logRho2_PC_limit, & ! don't use PC for logRho < this (>= 2.8)
      PC_use_Gamma_limit_instead_of_T, &
      logT1_PC_limit, & ! okay for pure PC for logT < this (like logT_all_OPAL)
      logT2_PC_limit, & ! don't use PC for logT > this (like logT_all_HELM)
      log_Gamma_e_all_HELM, & ! HELM for log_Gamma_e <= this
      log_Gamma_e_all_PC, & ! PC for log_Gamma_e >= this
      ! crystallization boundaries
      PC_Gamma_start_crystal, & ! Begin releasing latent heat of crystallization
      PC_Gamma_full_crystal, & ! Fully into the solid phase

      ! controls for Skye
      use_Skye, &
      Skye_use_ion_offsets, &
      mass_fraction_limit_for_Skye, &
      Skye_min_gamma_for_solid, & ! The minimum Gamma_i at which to use the solid free energy fit (below this, extrapolate).
      Skye_max_gamma_for_liquid, & ! The maximum Gamma_i at which to use the liquid free energy fit (above this, extrapolate).
      Skye_solid_mixing_rule, &

      use_simple_Skye_blends, &
      logRho_min_for_any_Skye, &
      logRho_min_for_all_Skye, &
      logT_min_for_any_Skye, &
      logT_min_for_all_Skye, &

      ! misc
      include_radiation, &
      include_elec_pos, &
      eosDT_use_linear_interp_for_X, &
      eosDT_use_linear_interp_to_HELM, &
      eosDT_file_prefix, &
      
      okay_to_convert_ierr_to_skip, &
      tiny_fuzz, &

      ! other eos
      use_other_eos_component, use_other_eos_results, &

      ! debugging
      dbg, &
      logT_lo, logT_hi, &
      logRho_lo, logRho_hi, &
      X_lo, X_hi, &
      Z_lo, Z_hi, &
      
      read_extra_eos_inlist1, extra_eos_inlist1_name, &
      read_extra_eos_inlist2, extra_eos_inlist2_name, &
      read_extra_eos_inlist3, extra_eos_inlist3_name, &
      read_extra_eos_inlist4, extra_eos_inlist4_name, &
      read_extra_eos_inlist5, extra_eos_inlist5_name


   contains


   ! read a "namelist" file and set parameters
   subroutine read_namelist(handle, inlist, ierr)
      integer, intent(in) :: handle
      character (len=*), intent(in) :: inlist
      integer, intent(out) :: ierr ! 0 means AOK.
      type (EoS_General_Info), pointer :: rq
      integer :: iz, j
      include 'formats'
      call get_eos_ptr(handle,rq,ierr)
      if (ierr /= 0) return
      call set_default_controls
      call read_controls_file(rq, inlist, 1, ierr)
      if (ierr /= 0) return
      rq% Gamma_e_all_HELM = exp10(rq% log_Gamma_e_all_HELM)
      if (FreeEOS_XZ_struct% Zs(num_FreeEOS_Zs) /= 1d0) then
         write(*,*) 'ERROR: expect FreeEOS_XZ_struct% Zs(num_FreeEOS_Zs) == 1d0'
         call mesa_error(__FILE__,__LINE__,'init_eos_handle_data')
      end if
   end subroutine read_namelist


   recursive subroutine read_controls_file(rq, filename, level, ierr)
      use ISO_FORTRAN_ENV, only: IOSTAT_END
      character(*), intent(in) :: filename
      type (EoS_General_Info), pointer :: rq
      integer, intent(in) :: level
      integer, intent(out) :: ierr
      logical :: read_extra1, read_extra2, read_extra3, read_extra4, read_extra5
      character (len=128) :: message, extra1, extra2, extra3, extra4, extra5
      integer :: unit

      ierr = 0
      if (level >= 10) then
         write(*,*) 'ERROR: too many levels of nested extra controls inlist files'
         ierr = -1
         return
      end if
      
      if (len_trim(filename) > 0) then
         open(newunit=unit, file=trim(filename), &
            action='read', delim='quote', status='old', iostat=ierr)
         if (ierr /= 0) then
            if (level == 1) then
               ierr = 0 ! no inlist file so just use defaults
               call store_controls(rq)
            else
               write(*, *) 'Failed to open eos namelist file ', trim(filename)
            end if
            return
         end if
         read(unit, nml=eos, iostat=ierr)
         close(unit)
         if (ierr == IOSTAT_END) then ! end-of-file means didn't find an &eos namelist
            ierr = 0
            write(*, *) 'WARNING: Failed to find eos namelist in file: ', trim(filename)
            call store_controls(rq)
            close(unit)
            return
         else if (ierr /= 0) then
            write(*, *)
            write(*, *)
            write(*, *)
            write(*, *)
            write(*, '(a)') 'Failed while trying to read eos namelist file: ' // trim(filename)
            write(*, '(a)') 'Perhaps the following runtime error message will help you find the problem.'
            write(*, *)
            open(newunit=unit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
            read(unit, nml=eos)
            close(unit)
            return
         end if
      end if

      call store_controls(rq)
      
      if (len_trim(filename) == 0) return

      ! recursive calls to read other inlists

      read_extra1 = read_extra_eos_inlist1
      read_extra_eos_inlist1 = .false.
      extra1 = extra_eos_inlist1_name
      extra_eos_inlist1_name = 'undefined'

      read_extra2 = read_extra_eos_inlist2
      read_extra_eos_inlist2 = .false.
      extra2 = extra_eos_inlist2_name
      extra_eos_inlist2_name = 'undefined'

      read_extra3 = read_extra_eos_inlist3
      read_extra_eos_inlist3 = .false.
      extra3 = extra_eos_inlist3_name
      extra_eos_inlist3_name = 'undefined'

      read_extra4 = read_extra_eos_inlist4
      read_extra_eos_inlist4 = .false.
      extra4 = extra_eos_inlist4_name
      extra_eos_inlist4_name = 'undefined'

      read_extra5 = read_extra_eos_inlist5
      read_extra_eos_inlist5 = .false.
      extra5 = extra_eos_inlist5_name
      extra_eos_inlist5_name = 'undefined'

      if (read_extra1) then
         call read_controls_file(rq, extra1, level+1, ierr)
         if (ierr /= 0) return
      end if

      if (read_extra2) then
         call read_controls_file(rq, extra2, level+1, ierr)
         if (ierr /= 0) return
      end if

      if (read_extra3) then
         call read_controls_file(rq, extra3, level+1, ierr)
         if (ierr /= 0) return
      end if

      if (read_extra4) then
         call read_controls_file(rq, extra4, level+1, ierr)
         if (ierr /= 0) return
      end if

      if (read_extra5) then
         call read_controls_file(rq, extra5, level+1, ierr)
         if (ierr /= 0) return
      end if

   end subroutine read_controls_file


   subroutine set_default_controls
      include 'eos.defaults'
   end subroutine set_default_controls


   subroutine store_controls(rq)
      type (EoS_General_Info), pointer :: rq
      ! controls for HELM
      rq% Z_all_HELM = Z_all_HELM
      rq% logT_all_HELM = logT_all_HELM
      rq% logT_low_all_HELM = logT_low_all_HELM
      rq% coulomb_temp_cut_HELM = coulomb_temp_cut_HELM
      rq% coulomb_den_cut_HELM = coulomb_den_cut_HELM      
      ! controls for OPAL_SCVH
      rq% use_OPAL_SCVH = use_OPAL_SCVH
      rq% logT_low_all_SCVH = logT_low_all_SCVH
      rq% logT_all_OPAL = logT_all_OPAL
      rq% logRho1_OPAL_SCVH_limit = logRho1_OPAL_SCVH_limit
      rq% logRho2_OPAL_SCVH_limit = logRho2_OPAL_SCVH_limit
      rq% logRho_min_OPAL_SCVH_limit = logRho_min_OPAL_SCVH_limit
      rq% logQ_max_OPAL_SCVH = logQ_max_OPAL_SCVH
      rq% logQ_min_OPAL_SCVH = logQ_min_OPAL_SCVH
      rq% Z_all_OPAL = Z_all_OPAL      
      ! controls for FreeEOS
      rq% use_FreeEOS = use_FreeEOS
      rq% logQ_max_FreeEOS_hi = logQ_max_FreeEOS_hi
      rq% logQ_max_FreeEOS_lo = logQ_max_FreeEOS_lo
      rq% logQ_min_FreeEOS_hi = logQ_min_FreeEOS_hi
      rq% logQ_min_FreeEOS_lo = logQ_min_FreeEOS_lo
      rq% logRho_min_FreeEOS_hi = logRho_min_FreeEOS_hi
      rq% logRho_min_FreeEOS_lo = logRho_min_FreeEOS_lo
      rq% logRho_max_FreeEOS_hi = logRho_max_FreeEOS_hi
      rq% logRho_max_FreeEOS_lo = logRho_max_FreeEOS_lo
      rq% logT_min_FreeEOS_hi = logT_min_FreeEOS_hi
      rq% logT_min_FreeEOS_lo = logT_min_FreeEOS_lo
      rq% logT_max_FreeEOS_hi = logT_max_FreeEOS_hi
      rq% logT_max_FreeEOS_lo = logT_max_FreeEOS_lo
      rq% logQ_cut_FreeEOS_lo_Z_max = logQ_cut_FreeEOS_lo_Z_max
      rq% logQ_cut_lo_Z_FreeEOS_hi = logQ_cut_lo_Z_FreeEOS_hi
      rq% logQ_cut_lo_Z_FreeEOS_lo = logQ_cut_lo_Z_FreeEOS_lo
      rq% logQ_cut_hi_Z_FreeEOS_hi = logQ_cut_hi_Z_FreeEOS_hi
      rq% logQ_cut_hi_Z_FreeEOS_lo = logQ_cut_hi_Z_FreeEOS_lo
      rq% logRho_cut_FreeEOS_hi = logRho_cut_FreeEOS_hi
      rq% logRho_cut_FreeEOS_lo = logRho_cut_FreeEOS_lo
      rq% logT_cut_FreeEOS_hi = logT_cut_FreeEOS_hi
      rq% logT_cut_FreeEOS_lo = logT_cut_FreeEOS_lo
      rq% suffix_for_FreeEOS_Z(1:num_FreeEOS_Zs) = &
         suffix_for_FreeEOS_Z(1:num_FreeEOS_Zs)      
      ! controls for CMS
      rq% use_CMS = use_CMS
      rq% CMS_use_fixed_composition = CMS_use_fixed_composition
      rq% CMS_fixed_composition_index = CMS_fixed_composition_index
      rq% max_Z_for_any_CMS = max_Z_for_any_CMS
      rq% max_Z_for_all_CMS = max_Z_for_all_CMS
      rq% logQ_max_for_any_CMS = logQ_max_for_any_CMS
      rq% logQ_max_for_all_CMS = logQ_max_for_all_CMS
      rq% logQ_min_for_any_CMS = logQ_min_for_any_CMS
      rq% logQ_min_for_all_CMS = logQ_min_for_all_CMS
      rq% logRho_max_for_all_CMS = logRho_max_for_all_CMS
      rq% logRho_max_for_any_CMS = logRho_max_for_any_CMS
      rq% logRho_min_for_all_CMS = logRho_min_for_all_CMS
      rq% logRho_min_for_any_CMS = logRho_min_for_any_CMS
      rq% logT_max_for_all_CMS = logT_max_for_all_CMS
      rq% logT_max_for_any_CMS = logT_max_for_any_CMS
      rq% logT_min_for_all_CMS = logT_min_for_all_CMS
      rq% logT_min_for_any_CMS = logT_min_for_any_CMS
      rq% logT_max_for_all_CMS_pure_He = logT_max_for_all_CMS_pure_He
      rq% logT_max_for_any_CMS_pure_He = logT_max_for_any_CMS_pure_He      
      ! controls for PC
      rq% use_PC = use_PC
      rq% mass_fraction_limit_for_PC = mass_fraction_limit_for_PC
      rq% logRho1_PC_limit = logRho1_PC_limit
      rq% logRho2_PC_limit = logRho2_PC_limit
      rq% PC_use_Gamma_limit_instead_of_T = PC_use_Gamma_limit_instead_of_T
      rq% logT1_PC_limit = logT1_PC_limit
      rq% logT2_PC_limit = logT2_PC_limit
      rq% log_Gamma_e_all_HELM = log_Gamma_e_all_HELM
      rq% log_Gamma_e_all_PC = log_Gamma_e_all_PC
      rq% PC_Gamma_start_crystal = PC_Gamma_start_crystal
      rq% PC_Gamma_full_crystal = PC_Gamma_full_crystal
      ! controls for Skye
      rq% use_Skye = use_Skye
      rq% Skye_use_ion_offsets = Skye_use_ion_offsets
      rq% mass_fraction_limit_for_Skye = mass_fraction_limit_for_Skye
      rq%Skye_min_gamma_for_solid = Skye_min_gamma_for_solid
      rq%Skye_max_gamma_for_liquid = Skye_max_gamma_for_liquid
      rq%Skye_solid_mixing_rule = Skye_solid_mixing_rule
      rq% use_simple_Skye_blends = use_simple_Skye_blends
      rq% logRho_min_for_any_Skye = logRho_min_for_any_Skye
      rq% logRho_min_for_all_Skye = logRho_min_for_all_Skye
      rq% logT_min_for_any_Skye = logT_min_for_any_Skye
      rq% logT_min_for_all_Skye = logT_min_for_all_Skye

      ! misc
      rq% include_radiation = include_radiation
      rq% include_elec_pos = include_elec_pos
      rq% eosDT_use_linear_interp_for_X = eosDT_use_linear_interp_for_X
      rq% eosDT_use_linear_interp_to_HELM = eosDT_use_linear_interp_to_HELM      
      rq% eosDT_file_prefix = eosDT_file_prefix      
      rq% okay_to_convert_ierr_to_skip = okay_to_convert_ierr_to_skip
      rq% tiny_fuzz = tiny_fuzz

      ! other eos
      rq% use_other_eos_component = use_other_eos_component
      rq% use_other_eos_results = use_other_eos_results

      ! debugging
      rq% dbg = dbg
      rq% logT_lo = logT_lo
      rq% logT_hi = logT_hi
      rq% logRho_lo = logRho_lo
      rq% logRho_hi = logRho_hi
      rq% X_lo = X_lo
      rq% X_hi = X_hi
      rq% Z_lo = Z_lo
      rq% Z_hi = Z_hi
   end subroutine store_controls


   subroutine write_namelist(handle, filename, ierr)
      integer, intent(in) :: handle
      character(*), intent(in) :: filename
      integer, intent(out) :: ierr 
      type (EoS_General_Info), pointer :: rq
      integer :: iounit
      open(newunit=iounit, file=trim(filename), &
         action='write', status='replace', iostat=ierr)
      if (ierr /= 0) then
         write(*,*) 'failed to open ' // trim(filename)
         return
      endif
      call get_eos_ptr(handle,rq,ierr)
      if (ierr /= 0) then
         close(iounit)
         return
      end if      
      call set_controls_for_writing(rq)      
      write(iounit, nml=eos, iostat=ierr)
      close(iounit)
   end subroutine write_namelist


   subroutine set_controls_for_writing(rq)
      type (EoS_General_Info), pointer :: rq
      ! controls for HELM
      Z_all_HELM = rq% Z_all_HELM
      logT_all_HELM = rq% logT_all_HELM
      logT_low_all_HELM = rq% logT_low_all_HELM
      coulomb_temp_cut_HELM = rq% coulomb_temp_cut_HELM
      coulomb_den_cut_HELM = rq% coulomb_den_cut_HELM      
      ! controls for OPAL_SCVH
      use_OPAL_SCVH = rq% use_OPAL_SCVH
      logT_low_all_SCVH = rq% logT_low_all_SCVH
      logT_all_OPAL = rq% logT_all_OPAL
      logRho1_OPAL_SCVH_limit = rq% logRho1_OPAL_SCVH_limit
      logRho2_OPAL_SCVH_limit = rq% logRho2_OPAL_SCVH_limit
      logRho_min_OPAL_SCVH_limit = rq% logRho_min_OPAL_SCVH_limit
      logQ_max_OPAL_SCVH = rq% logQ_max_OPAL_SCVH
      logQ_min_OPAL_SCVH = rq% logQ_min_OPAL_SCVH
      Z_all_OPAL = rq% Z_all_OPAL      
      ! controls for FreeEOS
      use_FreeEOS = rq% use_FreeEOS
      logQ_max_FreeEOS_hi = rq% logQ_max_FreeEOS_hi
      logQ_max_FreeEOS_lo = rq% logQ_max_FreeEOS_lo
      logQ_min_FreeEOS_hi = rq% logQ_min_FreeEOS_hi
      logQ_min_FreeEOS_lo = rq% logQ_min_FreeEOS_lo
      logRho_min_FreeEOS_hi = rq% logRho_min_FreeEOS_hi
      logRho_min_FreeEOS_lo = rq% logRho_min_FreeEOS_lo
      logRho_max_FreeEOS_hi = rq% logRho_max_FreeEOS_hi
      logRho_max_FreeEOS_lo = rq% logRho_max_FreeEOS_lo
      logT_min_FreeEOS_hi = rq% logT_min_FreeEOS_hi
      logT_min_FreeEOS_lo = rq% logT_min_FreeEOS_lo
      logT_max_FreeEOS_hi = rq% logT_max_FreeEOS_hi
      logT_max_FreeEOS_lo = rq% logT_max_FreeEOS_lo
      logQ_cut_FreeEOS_lo_Z_max = rq% logQ_cut_FreeEOS_lo_Z_max
      logQ_cut_lo_Z_FreeEOS_hi = rq% logQ_cut_lo_Z_FreeEOS_hi
      logQ_cut_lo_Z_FreeEOS_lo = rq% logQ_cut_lo_Z_FreeEOS_lo
      logQ_cut_hi_Z_FreeEOS_hi = rq% logQ_cut_hi_Z_FreeEOS_hi
      logQ_cut_hi_Z_FreeEOS_lo = rq% logQ_cut_hi_Z_FreeEOS_lo
      logRho_cut_FreeEOS_hi = rq% logRho_cut_FreeEOS_hi
      logRho_cut_FreeEOS_lo = rq% logRho_cut_FreeEOS_lo
      logT_cut_FreeEOS_hi = rq% logT_cut_FreeEOS_hi
      logT_cut_FreeEOS_lo = rq% logT_cut_FreeEOS_lo
      suffix_for_FreeEOS_Z(1:num_FreeEOS_Zs) = &
         rq% suffix_for_FreeEOS_Z(1:num_FreeEOS_Zs)      
      ! controls for CMS
      use_CMS = rq% use_CMS
      CMS_use_fixed_composition = rq% CMS_use_fixed_composition
      CMS_fixed_composition_index = rq% CMS_fixed_composition_index
      max_Z_for_any_CMS = rq% max_Z_for_any_CMS
      max_Z_for_all_CMS = rq% max_Z_for_all_CMS
      logQ_max_for_any_CMS = rq% logQ_max_for_any_CMS
      logQ_max_for_all_CMS = rq% logQ_max_for_all_CMS
      logQ_min_for_any_CMS = rq% logQ_min_for_any_CMS
      logQ_min_for_all_CMS = rq% logQ_min_for_all_CMS
      logRho_max_for_all_CMS = rq% logRho_max_for_all_CMS
      logRho_max_for_any_CMS = rq% logRho_max_for_any_CMS
      logRho_min_for_all_CMS = rq% logRho_min_for_all_CMS
      logRho_min_for_any_CMS = rq% logRho_min_for_any_CMS
      logT_max_for_all_CMS = rq% logT_max_for_all_CMS
      logT_max_for_any_CMS = rq% logT_max_for_any_CMS
      logT_min_for_all_CMS = rq% logT_min_for_all_CMS
      logT_min_for_any_CMS = rq% logT_min_for_any_CMS
      logT_max_for_all_CMS_pure_He = rq% logT_max_for_all_CMS_pure_He
      logT_max_for_any_CMS_pure_He = rq% logT_max_for_any_CMS_pure_He      
      ! controls for PC
      use_PC = rq% use_PC
      mass_fraction_limit_for_PC = rq% mass_fraction_limit_for_PC
      logRho1_PC_limit = rq% logRho1_PC_limit
      logRho2_PC_limit = rq% logRho2_PC_limit
      PC_use_Gamma_limit_instead_of_T = rq% PC_use_Gamma_limit_instead_of_T
      logT1_PC_limit = rq% logT1_PC_limit
      logT2_PC_limit = rq% logT2_PC_limit
      log_Gamma_e_all_HELM = rq% log_Gamma_e_all_HELM
      log_Gamma_e_all_PC = rq% log_Gamma_e_all_PC
      PC_Gamma_start_crystal = rq% PC_Gamma_start_crystal
      PC_Gamma_full_crystal = rq% PC_Gamma_full_crystal
      ! controls for Skye
      use_Skye = rq% use_Skye
      Skye_use_ion_offsets = rq% Skye_use_ion_offsets
      mass_fraction_limit_for_Skye = rq% mass_fraction_limit_for_Skye   
      Skye_min_gamma_for_solid = rq% Skye_min_gamma_for_solid
      Skye_max_gamma_for_liquid = rq% Skye_max_gamma_for_liquid  
      Skye_solid_mixing_rule = rq% Skye_solid_mixing_rule
      use_simple_Skye_blends = rq% use_simple_Skye_blends
      logRho_min_for_any_Skye = rq% logRho_min_for_any_Skye
      logRho_min_for_all_Skye = rq% logRho_min_for_all_Skye
      logT_min_for_any_Skye = rq% logT_min_for_any_Skye
      logT_min_for_all_Skye = rq% logT_min_for_all_Skye

      ! misc
      include_radiation = rq% include_radiation
      include_elec_pos = rq% include_elec_pos
      eosDT_use_linear_interp_for_X = rq% eosDT_use_linear_interp_for_X
      eosDT_use_linear_interp_to_HELM = rq% eosDT_use_linear_interp_to_HELM      
      eosDT_file_prefix = rq% eosDT_file_prefix      
      okay_to_convert_ierr_to_skip = rq% okay_to_convert_ierr_to_skip
      tiny_fuzz = rq% tiny_fuzz

      ! other eos
      use_other_eos_component = rq% use_other_eos_component
      use_other_eos_results = rq% use_other_eos_results

      ! debugging
      dbg = rq% dbg
      logT_lo = rq% logT_lo
      logT_hi = rq% logT_hi
      logRho_lo = rq% logRho_lo
      logRho_hi = rq% logRho_hi
      X_lo = rq% X_lo
      X_hi = rq% X_hi
      Z_lo = rq% Z_lo
      Z_hi = rq% Z_hi
   end subroutine set_controls_for_writing
   

   subroutine get_eos_controls(rq, name, val, ierr)
      use utils_lib, only: StrUpCase
      type (EoS_General_Info), pointer :: rq
      character(len=*),intent(in) :: name
      character(len=*), intent(out) :: val
      integer, intent(out) :: ierr

      character(len(name)+1) :: upper_name
      character(len=512) :: str
      integer :: iounit,iostat,ind,i

      ierr = 0


      ! First save current controls
      call set_controls_for_writing(rq)

      ! Write namelist to temporay file
      open(newunit=iounit,status='scratch')
      write(iounit,nml=eos)
      rewind(iounit)

      ! Namelists get written in captials
      upper_name = trim(StrUpCase(name))//'='
      val = ''
      ! Search for name inside namelist
      do 
         read(iounit,'(A)',iostat=iostat) str
         ind = index(trim(str),trim(upper_name))
         if( ind /= 0 ) then
            val = str(ind+len_trim(upper_name):len_trim(str)-1) ! Remove final comma and starting =
            do i=1,len(val)
               if(val(i:i)=='"') val(i:i) = ' '
            end do
            exit
         end if
         if(is_iostat_end(iostat)) exit
      end do   

      if(len_trim(val) == 0 .and. ind==0 ) ierr = -1

      close(iounit)

   end subroutine get_eos_controls

   subroutine set_eos_controls(rq, name, val, ierr)
      type (EoS_General_Info), pointer :: rq
      character(len=*), intent(in) :: name, val
      character(len=len(name)+len(val)+8) :: tmp
      integer, intent(out) :: ierr

      ierr = 0

      ! First save current eos_controls
      call set_controls_for_writing(rq)

      tmp=''
      tmp = '&eos '//trim(name)//'='//trim(val)//' /'

      ! Load into namelist
      read(tmp, nml=eos)

      ! Add to eos
      call store_controls(rq)
      if(ierr/=0) return

   end subroutine set_eos_controls


   end module eos_ctrls_io

