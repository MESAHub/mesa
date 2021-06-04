! ***********************************************************************
!
!   Copyright (C) 2010-2021  The MESA Team
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

      module eos_def
      
      use const_def, only: dp, use_mesa_temp_cache
      use chem_def, only: max_el_z
      
      implicit none
  
      ! interfaces for procedure pointers
      abstract interface

         subroutine other_eos_frac_interface( &
            handle, &
            species, chem_id, net_iso, xa, &
            Rho, logRho, T, logT, &
            frac, dfrac_dlogRho, dfrac_dlogT, &
            ierr)
            use const_def, only: dp
            integer, intent(in) :: handle ! eos handle; from star, pass s% eos_handle
            integer, intent(in) :: species
            integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos
            integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
            real(dp), intent(in) :: xa(:) ! mass fractions

            real(dp), intent(in) :: Rho, logRho ! the density
            real(dp), intent(in) :: T, logT ! the temperature

            real(dp), intent(out) :: frac ! fraction of other_eos to use
            real(dp), intent(out) :: dfrac_dlogRho ! its partial derivative at constant T
            real(dp), intent(out) :: dfrac_dlogT   ! its partial derivative at constant Rho

            integer, intent(out) :: ierr ! 0 means AOK.
         end subroutine other_eos_frac_interface

         subroutine other_eos_interface( &
            handle, &
            species, chem_id, net_iso, xa, &
            Rho, logRho, T, logT, &
            res, d_dlnd, d_dlnT, d_dxa, ierr)
            use const_def, only: dp
            integer, intent(in) :: handle ! eos handle; from star, pass s% eos_handle

            integer, intent(in) :: species
            integer, pointer :: chem_id(:) ! maps species to chem id
            integer, pointer :: net_iso(:) ! maps chem id to species number
            real(dp), intent(in) :: xa(:) ! mass fractions

            real(dp), intent(in) :: Rho, logRho ! the density
            real(dp), intent(in) :: T, logT ! the temperature

            real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
            real(dp), intent(inout) :: d_dlnd(:) ! (num_eos_basic_results)
            real(dp), intent(inout) :: d_dlnT(:) ! (num_eos_basic_results)
            real(dp), intent(inout) :: d_dxa(:,:) ! (num_eos_basic_results,species)

            integer, intent(out) :: ierr ! 0 means AOK.
         end subroutine other_eos_interface

      end interface


      logical, parameter :: show_allocations = .false.  ! for debugging memory usage
      integer, parameter :: eos_name_length = 20 ! String length for storing EOS variable names
 
        
      ! cgs units

      ! the basic eos results
      
      integer, parameter :: i_lnPgas = 1
            ! gas pressure (total pressure minus radiation pressure)
      integer, parameter :: i_lnE = i_lnPgas+1 
            ! internal energy per gram
      integer, parameter :: i_lnS = i_lnE+1 
            ! entropy per gram
      integer, parameter :: i_mu = i_lnS+1 
            ! mean molecular weight per gas particle (ions + free electrons)
      integer, parameter :: i_lnfree_e = i_mu+1
            ! free_e := total combined number per nucleon of free electrons
      integer, parameter :: i_eta = i_lnfree_e+1 
            ! electron degeneracy parameter (eta > 1 for significant degeneracy)
            ! eta = ratio of electron chemical potential to kT
      integer, parameter :: i_grad_ad = i_eta+1
            ! dlnT_dlnP at constant S
      integer, parameter :: i_chiRho = i_grad_ad+1
            ! dlnP_dlnRho at constant T
      integer, parameter :: i_chiT = i_chiRho+1 
            ! dlnP_dlnT at constant Rho
      integer, parameter :: i_Cp = i_chiT+1
            ! dh_dT at constant P, specific heat at constant total pressure
            ! where h is enthalpy, h = E + P/Rho
      integer, parameter :: i_Cv = i_Cp+1 
            ! dE_dT at constant Rho, specific heat at constant volume
      integer, parameter :: i_dE_dRho = i_Cv+1
            ! at constant T
      integer, parameter :: i_dS_dT = i_dE_dRho+1 
            ! at constant Rho
      integer, parameter :: i_dS_dRho = i_dS_dT+1 
            ! at constant T
      integer, parameter :: i_gamma1 = i_dS_dRho+1 
            ! dlnP_dlnRho at constant S
      integer, parameter :: i_gamma3 = i_gamma1+1 
            ! gamma3 - 1 = dlnT_dlnRho at constant S
      integer, parameter :: i_phase = i_gamma3+1
            ! phase: 1 for solid, 0 for liquid, in-between for blend
      integer, parameter :: i_latent_ddlnT = i_phase+1
            ! T*dS/dlnT from phase transition
      integer, parameter :: i_latent_ddlnRho = i_latent_ddlnT+1
            ! T*dS/dlnRho from phase transition

      ! blend information
      integer, parameter :: i_eos_HELM = 1
      integer, parameter :: i_frac_HELM = i_latent_ddlnRho+1
            ! fraction HELM
      integer, parameter :: i_eos_OPAL_SCVH = 2
      integer, parameter :: i_frac_OPAL_SCVH = i_frac_HELM+1
            ! fraction OPAL/SCVH
      integer, parameter :: i_eos_FreeEOS = 3
      integer, parameter :: i_frac_FreeEOS = i_frac_OPAL_SCVH+1
            ! fraction FreeEOS
      integer, parameter :: i_eos_PC = 4
      integer, parameter :: i_frac_PC = i_frac_FreeEOS+1
            ! fraction PC
      integer, parameter :: i_eos_Skye = 5
      integer, parameter :: i_frac_Skye = i_frac_PC+1
            ! fraction Skye
      integer, parameter :: i_eos_CMS = 6
      integer, parameter :: i_frac_CMS = i_frac_Skye+1
            ! fraction CMS

      ! eos frac entries correspond to the slice
      ! i_frac:i_frac+num_eos_frac_results-1
      integer, parameter :: i_frac = i_frac_HELM ! first frac entry
      integer, parameter :: num_eos_frac_results = 6
      
      integer, parameter :: num_eos_basic_results = i_frac_CMS
      integer, parameter :: nv = num_eos_basic_results

      ! only return d_dxa of lnE and lnPgas to star
      integer, parameter :: num_eos_d_dxa_results = 2

      character (len=eos_name_length) :: eosDT_result_names(nv)

      ! non-positive indices for special EOS quantities
      ! these are supported in the root-finds (eosDT_get_{Rho,T})
      ! even though their values are not basic EOS results

      integer, parameter :: i_logPtot = 0  ! log10 total pressure (gas + radiation)
      integer, parameter :: i_egas = -1 ! gas specific energy density (no radiation)

      
      ! NOTE: the calculation of eta is based on the following equation for ne, 
      ! the mean number of free electrons per cm^3,
      ! assuming non-relativistic electrons (okay for T < 10^9 at least)
      !
      !  ne = 4 Pi / h^3 (2 me k T)^1.5 F[1/2,eta]   -- see, for example, Clayton, eqn 2-57
      !        where F is the fermi-dirac integral: 
      !        F[beta,eta] := Integrate[(u^beta)/(1+E^(u-eta)),{u,0,Infinity}]
      ! 
      ! CAVEAT: when free_e, the mean number of free electrons per nucleon gets really small, 
      ! eta isn't very interesting because there aren't a lot of free electrons to be degenerate!
      ! our calculation of eta gets flaky at this point as well.
      ! we sweep this problem under the rug by making eta tend to a fairly large negative value 
      ! when free_e < 0.01 or so. this roughly corresponds to T < 10^4 or less.
      
      integer, parameter :: eosdt_OPAL_SCVH = 1
      integer, parameter :: eosdt_max_FreeEOS = 2

      integer, parameter :: max_num_DT_Zs = 15
      integer, parameter :: max_num_DT_Xs = 12

      integer, parameter :: num_eosDT_Zs = 3
      integer, parameter :: num_eosDT_Xs = 6
      
      integer, parameter :: num_FreeEOS_Zs = 15
      integer, parameter :: num_FreeEOS_Xs = 12

      type DT_XZ_Info
         integer :: nZs
         real(dp) :: Zs(max_num_DT_Zs)
         integer :: nXs_for_Z(max_num_DT_Zs)
         real(dp) :: Xs_for_Z(max_num_DT_Xs, max_num_DT_Zs)
      end type DT_XZ_Info
      
      type (DT_XZ_Info), target :: eosDT_XZ_struct, FreeEOS_XZ_struct

      integer, parameter :: sz_per_eos_point = 4 ! for bicubic spline interpolation

      type EoS_General_Info

         ! limits for HELM
         real(dp) :: Z_all_HELM ! all HELM for Z >= this unless eos_use_FreeEOS
         real(dp) :: logT_all_HELM ! all HELM for lgT >= this
         real(dp) :: logT_low_all_HELM ! all HELM for lgT <= this
         real(dp) :: logT_ion_HELM, logT_neutral_HELM, max_logRho_neutral_HELM
         real(dp) :: coulomb_temp_cut_HELM, coulomb_den_cut_HELM
         
         ! limits for OPAL_SCVH
         logical :: use_OPAL_SCVH
         real(dp) :: logT_low_all_SCVH ! SCVH for lgT >= this
         real(dp) :: logT_all_OPAL ! OPAL for lgT <= this
         real(dp) :: logRho1_OPAL_SCVH_limit ! don't use OPAL_SCVH for logRho > this
         real(dp) :: logRho2_OPAL_SCVH_limit ! full OPAL_SCVH okay for logRho < this
         real(dp) :: logRho_min_OPAL_SCVH_limit ! no OPAL/SCVH for logRho < this
         real(dp) :: logQ_max_OPAL_SCVH ! no OPAL/SCVH for logQ > this
         real(dp) :: logQ_min_OPAL_SCVH ! no OPAL/SCVH for logQ <= this.
         real(dp) :: Z_all_OPAL ! all OPAL for Z <= this
         
         ! limits for FreeEOS
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
         
         ! limits for CMS
         logical :: use_CMS, CMS_use_fixed_composition
         integer :: CMS_fixed_composition_index ! in [0,10]
         real(dp) :: max_Z_for_any_CMS, max_Z_for_all_CMS ! set to -1 to disable CMS
         real(dp) :: logQ_max_for_any_CMS, logQ_max_for_all_CMS      ! for upper blend zone in logQ = logRho - 2*logT + 12
         real(dp) :: logQ_min_for_any_CMS, logQ_min_for_all_CMS      ! for lower blend zone in logQ
         real(dp) :: logRho_max_for_all_CMS, logRho_max_for_any_CMS  ! for upper blend zone in logRho
         real(dp) :: logRho_min_for_all_CMS, logRho_min_for_any_CMS  ! for lower blend zone in logRho
         real(dp) :: logT_max_for_all_CMS, logT_max_for_any_CMS      ! for upper blend zone in logT
         real(dp) :: logT_min_for_all_CMS, logT_min_for_any_CMS      ! for lower blend zone in logT      
         real(dp) :: logT_max_for_all_CMS_pure_He, logT_max_for_any_CMS_pure_He ! upper logT blend zone is different for pure He
         
         ! limits for PC
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
         real(dp) :: tiny_fuzz
         ! crystallization boundaries
         real(dp) :: PC_Gamma_start_crystal ! Begin releasing latent heat of crystallization
         real(dp) :: PC_Gamma_full_crystal ! Fully into the solid phase

         ! limits for Skye
         logical :: use_Skye
         real(dp) :: mass_fraction_limit_for_Skye
         real(dp) :: Skye_min_gamma_for_solid ! The minimum Gamma_i at which to use the solid free energy fit (below this, extrapolate).
         real(dp) :: Skye_max_gamma_for_liquid ! The maximum Gamma_i at which to use the liquid free energy fit (above this, extrapolate).
         character(len=128) :: Skye_solid_mixing_rule ! Currently support 'Ogata' or 'PC'

         logical :: use_simple_Skye_blends
         real(dp) :: logRho_min_for_any_Skye, logRho_min_for_all_Skye
         real(dp) :: logT_min_for_any_Skye, logT_min_for_all_Skye

         ! misc
         logical :: include_radiation, always_skip_elec_pos, always_include_elec_pos
         logical :: eosDT_use_linear_interp_for_X
         logical :: eosDT_use_linear_interp_to_HELM
      
         character(len=128) :: eosDT_file_prefix

         logical :: okay_to_convert_ierr_to_skip
         
         ! other eos
         logical :: use_other_eos_component
         procedure (other_eos_frac_interface), pointer, nopass :: &
         other_eos_frac => null()
         procedure (other_eos_interface), pointer, nopass :: &
            other_eos_component => null()
         logical :: use_other_eos_results
         procedure (other_eos_interface), pointer, nopass :: &
            other_eos_results => null()
         
         ! debugging
         logical :: dbg
         real(dp) :: logT_lo, logT_hi
         real(dp) :: logRho_lo, logRho_hi
         real(dp) :: X_lo, X_hi
         real(dp) :: Z_lo, Z_hi
         
         ! bookkeeping
         integer :: handle
         logical :: in_use

      end type EoS_General_Info

      
      include 'helm_def.dek'


      ! THE FOLLOWING ARE PRIVATE DEFS -- NOT FOR USE BY CLIENTS
      
      
      ! data table types
      
      type (HELM_Table), pointer :: eos_ht
      
      ! for mesa (logQ,logT) tables
      type EosDT_XZ_Info
         real(dp) :: logQ_min ! logQ = logRho - 2*logT + 12
         real(dp) :: logQ_max
         real(dp) :: del_logQ ! spacing for the logQs
         integer :: num_logQs
         real(dp) :: logT_min
         real(dp) :: logT_max
         real(dp) :: del_logT ! spacing for the logTs
         integer :: num_logTs
         real(dp), pointer :: logQs(:), logTs(:)
         real(dp), pointer :: tbl1(:)
         integer :: version
      end type EosDT_XZ_Info


      ! data table variables
      type (EosDT_XZ_Info), dimension(num_eosDT_Xs, num_eosDT_Zs), target :: &
         eosDT_XZ_data, eosSCVH_XZ_data, eosCMS_XZ_data
      type (EosDT_XZ_Info), dimension(num_FreeEOS_Xs, num_FreeEOS_Zs), target :: &
         FreeEOS_XZ_data

      logical, dimension(num_eosDT_Xs, num_eosDT_Zs) :: &
         eosDT_XZ_loaded, eosSCVH_XZ_loaded, eosCMS_XZ_loaded
      logical, dimension(num_FreeEOS_Xs, num_FreeEOS_Zs) :: FreeEOS_XZ_loaded
      

      ! interpolation info for eosPC support tables FITION9
      type FITION_Info
         real(dp), pointer :: f1(:), f(:,:)
      end type FITION_Info
      integer, parameter :: FITION_vals = 6
      real(dp), pointer :: FITION9_lnGAMIs(:)
      type (FITION_Info), target :: FITION_data(FITION_vals)
      logical :: FITION9_loaded ! all get created at once

      ! interpolation info for eosPC support tables, FSCRliq8 and EXCOR7
      type eosPC_Support_Info
         integer :: Zion, nlnRS, nlnGAME, nvals
         real(dp) :: lnRS_min, lnRS_max, dlnRS
         real(dp) :: lnGAME_min, lnGAME_max, dlnGAME
         real(dp), pointer :: lnRSs(:) ! (nlnRS)
         real(dp), pointer :: lnGAMEs(:) ! (nlnGAME)
         real(dp), pointer :: tbl1(:) ! (4,nvals,nlnRS,nlnGAME) 
         real(dp), pointer :: tbl(:,:,:,:) ! => tbl1(:)
      end type eosPC_Support_Info
      
      integer, parameter :: max_FSCRliq8_Zion = max_el_z
      type (eosPC_Support_Info), target :: FSCRliq8_data(max_FSCRliq8_Zion)
      logical, dimension(max_FSCRliq8_Zion) :: FSCRliq8_Zion_loaded

      type (eosPC_Support_Info), target :: EXCOR7_data
      logical :: EXCOR7_table_loaded

      integer, parameter :: max_eos_handles = 10
      type (EoS_General_Info), target :: eos_handles(max_eos_handles)
      
      logical :: use_cache_for_eos = .true.
      logical :: eos_root_is_initialized = .false.

      character(len=1000) :: eosDT_cache_dir
      character(len=1000) :: eosDT_temp_cache_dir

      logical :: eos_test_partials
      real(dp) :: eos_test_partials_val, eos_test_partials_dval_dx ! for dfridr from star
      
      
      contains
      
      
      subroutine eos_def_init
         integer :: i
         type (DT_XZ_Info), pointer :: eosDT_XZ_ptr, FreeEOS_XZ_ptr
         include 'formats'
         
         use_cache_for_eos = .true.
         eos_root_is_initialized = .false.
         eos_test_partials = .false.

         EXCOR7_table_loaded = .false.
         FSCRliq8_Zion_loaded(:) = .false.
         FITION9_loaded = .false.

         do i=1,max_eos_handles
            eos_handles(i)% handle = i
            eos_handles(i)% in_use = .false.
         end do

         eosDT_XZ_ptr => eosDT_XZ_struct
         eosDT_XZ_ptr% nZs = num_eosDT_Zs
         eosDT_XZ_ptr% Zs(1:num_eosDT_Zs) = (/ 0.00d0, 0.02d0, 0.04d0 /)
         eosDT_XZ_ptr% nXs_for_Z(1:num_eosDT_Zs) = (/ 6, 5, 5 /)
         do i=1,num_eosDT_Zs
            eosDT_XZ_ptr% Xs_for_Z(1:num_eosDT_Xs,i) = &
               (/ 0.0d0, 0.2d0, 0.4d0, 0.6d0, 0.8d0, 1.0d0 /)
         end do

         FreeEOS_XZ_ptr => FreeEOS_XZ_struct
         FreeEOS_XZ_ptr% nZs = num_FreeEOS_Zs
         FreeEOS_XZ_ptr% Zs(1:num_FreeEOS_Zs) = &
            (/ 0.00d0, 0.02d0, 0.04d0, 0.06d0, 0.08d0, 0.1d0, 0.2d0, &
               0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0, 0.9d0, 1.0d0 /)
         FreeEOS_XZ_ptr% nXs_for_Z(1:num_FreeEOS_Zs) = &
            (/ 11, 11, 11, 11, 11, 4, 3, 3, 3, 3, 3, 3, 3, 2, 1 /)
         do i=1,5 ! 0.0 to 0.08
            FreeEOS_XZ_ptr% Xs_for_Z(1:11,i) = &
               (/ 0.0d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, &
                  0.6d0, 0.7d0, 0.8d0, 0.9d0, 1.0d0-FreeEOS_XZ_ptr% Zs(i) /)
         end do
         FreeEOS_XZ_ptr% Xs_for_Z(1:4,6) = (/ 0.0d0, 0.1d0, 0.4d0, 0.9d0 /) ! 0.1
         do i=7,13 ! 0.2, 0.8
            FreeEOS_XZ_ptr% Xs_for_Z(1:3,i) = &
               (/ 0.0d0, 0.1d0, 1.0d0-FreeEOS_XZ_ptr% Zs(i) /)
         end do
         FreeEOS_XZ_ptr% Xs_for_Z(1:2,14) = (/ 0.0d0, 0.1d0 /) ! 0.9
         FreeEOS_XZ_ptr% Xs_for_Z(1,15) = 0.0d0 ! 1.0
         
         eosDT_XZ_loaded(:,:) = .false.
         eosSCVH_XZ_loaded(:,:)=.false.
         FreeEOS_XZ_loaded(:,:)=.false.
         eosCMS_XZ_loaded(:,:)=.false.
         
         eosDT_result_names(i_lnPgas) = 'lnPgas'
         eosDT_result_names(i_lnE) = 'lnE'
         eosDT_result_names(i_lnS) = 'lnS'
         eosDT_result_names(i_grad_ad) = 'grad_ad'
         eosDT_result_names(i_chiRho) = 'chiRho'
         eosDT_result_names(i_chiT) = 'chiT'
         eosDT_result_names(i_Cp) = 'Cp'
         eosDT_result_names(i_Cv) = 'Cv'
         eosDT_result_names(i_dE_dRho) = 'dE_dRho'
         eosDT_result_names(i_dS_dT) = 'dS_dT'
         eosDT_result_names(i_dS_dRho) = 'dS_dRho'
         eosDT_result_names(i_mu) = 'mu'
         eosDT_result_names(i_lnfree_e) = 'lnfree_e'
         eosDT_result_names(i_gamma1) = 'gamma1'
         eosDT_result_names(i_gamma3) = 'gamma3'
         eosDT_result_names(i_eta) = 'eta'
         eosDT_result_names(i_phase) = 'phase'
         eosDT_result_names(i_latent_ddlnT) = 'latent_ddlnT'
         eosDT_result_names(i_latent_ddlnRho) = 'latent_ddlnRho'
         eosDT_result_names(i_frac_OPAL_SCVH) = 'OPAL/SCVH'
         eosDT_result_names(i_frac_HELM) = 'HELM'
         eosDT_result_names(i_frac_Skye) = 'Skye'
         eosDT_result_names(i_frac_PC) = 'PC'
         eosDT_result_names(i_frac_FreeEOS) = 'FreeEOS'
         eosDT_result_names(i_frac_CMS) = 'CMS'

      end subroutine eos_def_init

      
      integer function do_alloc_eos(ierr) result(handle)
         integer, intent(out) :: ierr
         integer :: i
         ierr = 0
         handle = -1
!$omp critical (eos_handle)
         do i = 1, max_eos_handles
            if (.not. eos_handles(i)% in_use) then
               eos_handles(i)% in_use = .true.
               handle = i
               exit
            end if
         end do
!$omp end critical (eos_handle)
         if (handle == -1) then
            ierr = -1
            return
         end if
         if (eos_handles(handle)% handle /= handle) then
            ierr = -1
            return
         end if
         call init_eos_handle_data(handle)
      end function do_alloc_eos


      subroutine init_eos_handle_data(handle)
         use other_eos
         use math_lib
         integer, intent(in) :: handle
         type (EoS_General_Info), pointer :: rq
         integer :: iz
         rq => eos_handles(handle)
         rq% in_use = .true.
         rq% handle = handle

         rq% other_eos_frac => null_other_eos_frac
         rq% other_eos_component => null_other_eos_component
         rq% other_eos_results => null_other_eos_results
         
      end subroutine init_eos_handle_data
            
      
      subroutine do_free_eos_handle(handle)
         integer, intent(in) :: handle
         type (EoS_General_Info), pointer :: rq
         if (handle >= 1 .and. handle <= max_eos_handles) then
            rq => eos_handles(handle)
            eos_handles(handle)% in_use = .false.
         end if
      end subroutine do_free_eos_handle
      

      subroutine get_eos_ptr(handle,rq,ierr)
         integer, intent(in) :: handle
         type (EoS_General_Info), pointer :: rq
         integer, intent(out):: ierr         
         if (handle < 1 .or. handle > max_eos_handles) then
            ierr = -1
            return
         end if
         rq => eos_handles(handle)
         ierr = 0
      end subroutine get_eos_ptr
      
      
      subroutine eos_def_shutdown
         type (eosPC_Support_Info), pointer :: fq
         type (FITION_Info), pointer :: fi
         integer :: ix, iz

         ! DT tables
         call free_EosDT_XZ_Info( &
            eosDT_XZ_data, eosDT_XZ_loaded, num_eosDT_Xs, num_eosDT_Zs)
         call free_EosDT_XZ_Info( &
            eosSCVH_XZ_data, eosSCVH_XZ_loaded, num_eosDT_Xs, num_eosDT_Zs)
         call free_EosDT_XZ_Info( &
            eosCMS_XZ_data, eosCMS_XZ_loaded, num_eosDT_Xs, num_eosDT_Zs)
         call free_EosDT_XZ_Info( &
            FreeEOS_XZ_data, FreeEOS_XZ_loaded, num_FreeEOS_Xs, num_FreeEOS_Zs)

         ! Misc. stuff

         if (FITION9_loaded) then
            do iz = 1, FITION_vals
               fi => FITION_data(iz)
               if (ASSOCIATED(fi% f1)) deallocate(fi% f1)
               nullify(fi% f)
            end do
            if (ASSOCIATED(FITION9_lnGAMIs)) deallocate(FITION9_lnGAMIs)
            FITION9_loaded = .false.
         end if

         fq => EXCOR7_data
         call free_eosPC_support_Info(fq)
         EXCOR7_table_loaded = .false.
         
         do iz = 1, max_FSCRliq8_Zion
            if (.not. FSCRliq8_Zion_loaded(iz)) cycle
            fq => FSCRliq8_data(iz)
            call free_eosPC_support_Info(fq)
         end do
         FSCRliq8_Zion_loaded(:) = .false.
         
         eos_root_is_initialized = .false.
         
         
         contains
         
         subroutine free_eosPC_support_Info(fq)
            type (eosPC_Support_Info), pointer :: fq
            if (ASSOCIATED(fq% lnRSs)) deallocate(fq% lnRSs)
            if (ASSOCIATED(fq% lnGAMEs)) deallocate(fq% lnGAMEs)
            if (ASSOCIATED(fq% tbl1)) deallocate(fq% tbl1)
            nullify(fq% tbl)
         end subroutine free_eosPC_support_Info
         
         subroutine free_EosDT_XZ_Info(d, flgs, numXs, numZs)
            integer, intent(in) :: numXs, numZs
            type (EosDT_XZ_Info), dimension(numXs, numZs) :: d
            logical, dimension(numXs, numZs) :: flgs
            integer :: ix, iz
            do ix = 1, numXs
               do iz = 1, numZs
                  if (ASSOCIATED(d(ix,iz)% tbl1)) then
                     deallocate(d(ix,iz)% tbl1)
                  end if
                  if (ASSOCIATED(d(ix,iz)% logQs)) then
                     deallocate(d(ix,iz)% logQs)
                  end if
                  if (ASSOCIATED(d(ix,iz)% logTs)) then
                     deallocate(d(ix,iz)% logTs)
                  end if
               end do
            end do
            flgs(:,:) = .false.
         end subroutine free_EosDT_XZ_Info

      end subroutine eos_def_shutdown


      end module eos_def
      
