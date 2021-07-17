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


      module eos_lib

      use const_def, only: dp, crad
      use math_lib

      implicit none

      contains ! the procedure interface for the library
      ! client programs should only call these routines.
            
            
      subroutine eos_init( &
            eosDT_cache_dir, use_cache, info)      
         use eos_initialize, only : Init_eos
         character(*), intent(in) :: eosDT_cache_dir ! blank string means use default
         logical, intent(in) :: use_cache
         integer, intent(out) :: info ! 0 means AOK. 
         info = 0
         call Init_eos( &
            eosDT_cache_dir, use_cache, info)
         if (info /= 0) return      
      end subroutine eos_init      
      
      subroutine eos_shutdown
         use eos_def
         use helm_alloc,only:free_helm_table
         if (associated(eos_ht)) call free_helm_table(eos_ht)
         call eos_def_shutdown
      end subroutine eos_shutdown
      
      
      ! after eos_init has finished, you can allocate a "handle"
      ! and set control parameter values using an inlist
      
      integer function alloc_eos_handle(ierr) result(handle)
         integer, intent(out) :: ierr ! 0 means AOK.  
         character (len=0) :: inlist 
         handle = alloc_eos_handle_using_inlist(inlist, ierr)
      end function alloc_eos_handle      
      
      integer function alloc_eos_handle_using_inlist(inlist,ierr) result(handle)
         use eos_def, only:do_alloc_eos
         use eos_ctrls_io, only:read_namelist
         character (len=*), intent(in) :: inlist ! empty means just use defaults.
         integer, intent(out) :: ierr ! 0 means AOK.   
         ierr = 0
         handle = do_alloc_eos(ierr)
         if (ierr /= 0) return
         call read_namelist(handle, inlist, ierr)
      end function alloc_eos_handle_using_inlist      
      
      subroutine free_eos_handle(handle)
         ! frees the handle and all associated data
         use eos_def,only:do_free_eos_handle
         integer, intent(in) :: handle
         call do_free_eos_handle(handle)
      end subroutine free_eos_handle      


      subroutine eos_ptr(handle,rq,ierr)
         use eos_def,only:EoS_General_Info,get_eos_ptr
         integer, intent(in) :: handle ! from alloc_eos_handle
         type (EoS_General_Info), pointer :: rq
         integer, intent(out):: ierr
         call get_eos_ptr(handle,rq,ierr)
      end subroutine eos_ptr

      
      ! as a convenience
      
      elemental real(dp) function Radiation_Pressure(T)
         use const_def, only: crad
         real(dp), intent(in) :: T
         Radiation_Pressure = crad*T*T*T*T/3d0
      end function Radiation_Pressure
      
      elemental real(dp) function Radiation_Energy(T)
         use const_def, only: crad
         real(dp), intent(in) :: T
         Radiation_Energy = crad*T*T*T*T
      end function Radiation_Energy
      
      
      
      ! eos evaluation
      
      ! you can call these routines after you've allocated a handle.
      ! NOTE: the information referenced via the handle is read-only,
      ! so you can do multiple evaluations in parallel using the same handle.      
      

      subroutine eosDT_get( &
               handle, species, chem_id, net_iso, xa, &
               Rho, logRho, T, logT, &
               res, d_dlnd, d_dlnT, d_dxa, ierr)
         use eos_def
         use eosDT_eval, only: Get_eosDT_Results
         use chem_lib, only: basic_composition_info
         integer, intent(in) :: handle ! eos handle; from star, pass s% eos_handle
         integer, intent(in) :: species ! number of species
         integer, pointer :: chem_id(:) ! maps species to chem id
         integer, pointer :: net_iso(:) ! maps chem id to species number
         real(dp), intent(in) :: xa(:) ! mass fractions         
         real(dp), intent(in) :: Rho, logRho ! the density
         real(dp), intent(in) :: T, logT ! the temperature         
         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)         
         real(dp), intent(inout) :: d_dlnd(:) ! (num_eos_basic_results) 
         real(dp), intent(inout) :: d_dlnT(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dxa(:,:) ! (num_eos_d_dxa_results,species)
         integer, intent(out) :: ierr ! 0 means AOK.
         real(dp), allocatable :: d_dxa_eos(:,:) ! eos internally returns derivs of all quantities
         type (EoS_General_Info), pointer :: rq
         real(dp) :: X, Y, Z, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx
         call get_eos_ptr(handle,rq,ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for eos_get -- did you call alloc_eos_handle?'
            return
         end if
         call basic_composition_info( &
            species, chem_id, xa, X, Y, Z, &
            abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)
         allocate(d_dxa_eos(num_eos_basic_results, species))
         call Get_eosDT_Results( &
            rq, Z, X, abar, zbar, &
            species, chem_id, net_iso, xa, &
            Rho, logRho, T, logT, &
            res, d_dlnd, d_dlnT, d_dxa_eos, ierr)
         ! only return 1st two d_dxa results (lnE and lnPgas) to star
         d_dxa(1:num_eos_d_dxa_results,1:species) = d_dxa_eos(1:num_eos_d_dxa_results, 1:species)
      end subroutine eosDT_get
      

      subroutine eosDT_get_component( &
               handle, which_eos, &
               species, chem_id, net_iso, xa, &
               Rho, log10Rho, T, log10T, &
               res, d_dlnRho_const_T, d_dlnT_const_Rho, d_dxa_const_TRho, &
               ierr)

         use chem_lib, only: basic_composition_info
         use eos_def
         use eosDT_eval, only: Test_one_eosDT_component

         ! INPUT
         
         integer, intent(in) :: handle ! eos handle; from star, pass s% eos_handle
         
         integer, intent(in) :: which_eos ! see eos_def: i_eos_<component>

         integer, intent(in) :: species ! number of species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos         
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions
         
         real(dp), intent(in) :: Rho, log10Rho ! the density
            ! provide both if you have them.  else pass one and set the other to arg_not_provided
            ! "arg_not_provided" is defined in mesa const_def
            
         real(dp), intent(in) :: T, log10T ! the temperature
            ! provide both if you have them.  else pass one and set the other to arg_not_provided
                     
         ! OUTPUT
         
         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
         ! partial derivatives of the basic results wrt lnd and lnT
         
         real(dp), intent(inout) :: d_dlnRho_const_T(:) ! (num_eos_basic_results) 
         ! d_dlnRho_const_T(i) = d(res(i))/dlnd|T,X where X = composition
         real(dp), intent(inout) :: d_dlnT_const_Rho(:) ! (num_eos_basic_results) 
         ! d_dlnT(i) = d(res(i))/dlnT|Rho,X where X = composition
         real(dp), intent(inout) :: d_dxa_const_TRho(:,:) ! (num_eos_d_dxa_results,species)
         
         integer, intent(out) :: ierr ! 0 means AOK.

         real(dp), allocatable :: d_dxa_eos(:,:) ! eos internally returns derivs of all quantities

         type (EoS_General_Info), pointer :: rq

         real(dp) :: X, Y, Z, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx
         
         call get_eos_ptr(handle,rq,ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for eos_get -- did you call alloc_eos_handle?'
            return
         end if

         call basic_composition_info( &
            species, chem_id, xa, X, Y, Z, &
            abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)

         allocate(d_dxa_eos(num_eos_basic_results,species))

         call Test_one_eosDT_component( &
               rq, which_eos, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               Rho, log10Rho, T, log10T, &
               res, d_dlnRho_const_T, d_dlnT_const_Rho, d_dxa_eos, ierr)

         ! only return 1st two d_dxa results (lnE and lnPgas)
         d_dxa_const_TRho(1:num_eos_d_dxa_results,1:species) = d_dxa_eos(1:num_eos_d_dxa_results, 1:species)
         
      end subroutine eosDT_get_component


      subroutine helmeos2_eval( &
            T, logT, Rho, logRho, abar, zbar, &
            coulomb_temp_cut, coulomb_den_cut, helm_res, &
            clip_to_table_boundaries, include_radiation, &
            include_elec_pos, &
            off_table, ierr)         
         use helm
         real(dp), intent(in) :: T, logT, Rho, logRho, abar, zbar, &
            coulomb_temp_cut, coulomb_den_cut
         real(dp), intent(inout) :: helm_res(:) ! (num_helm_results)
         logical, intent(in) :: clip_to_table_boundaries, include_radiation, &
            include_elec_pos
         logical, intent(out) :: off_table
         integer, intent(out) :: ierr ! 0 means AOK.
         call helmeos2( &
            T, logT, Rho, logRho, abar, zbar, coulomb_temp_cut, coulomb_den_cut, &
            helm_res, clip_to_table_boundaries, include_radiation, include_elec_pos, &
            off_table, ierr)
      end subroutine helmeos2_eval
      
      
      ! the following routine uses gas pressure and temperature as input variables
      subroutine eosPT_get( &
               handle, &
               species, chem_id, net_iso, xa, &
               Pgas, log10Pgas, T, log10T, &
               Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               d_dxa_const_TRho, ierr)

         use chem_lib, only: basic_composition_info
         use eos_def
         use eosPT_eval, only: Get_eosPT_Results

         ! INPUT

         integer, intent(in) :: handle ! eos handle; from star, pass s% eos_handle

         integer, intent(in) :: species ! number of species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions

         real(dp), intent(in) :: Pgas, log10Pgas ! the gas pressure
            ! provide both if you have them.  else pass one and set the other to arg_not_provided
            ! "arg_not_provided" is defined in mesa const_def

         real(dp), intent(in) :: T, log10T ! the temperature
            ! provide both if you have them.  else pass one and set the other to arg_not_provided

         ! OUTPUT

         real(dp), intent(out) :: Rho, log10Rho ! density
         real(dp), intent(out) :: dlnRho_dlnPgas_const_T
         real(dp), intent(out) :: dlnRho_dlnT_const_Pgas

         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)

         ! partial derivatives of the basic results

         real(dp), intent(inout) :: d_dlnRho_const_T(:) ! (num_eos_basic_results)
         ! d_dlnRho_const_T(i) = d(res(i))/dlnd|T,X where X = composition
         real(dp), intent(inout) :: d_dlnT_const_Rho(:) ! (num_eos_basic_results)
         ! d_dlnT_const_Rho(i) = d(res(i))/dlnT|Rho,X where X = composition
         real(dp), intent(inout) :: d_dxa_const_TRho(:,:) ! (num_eos_d_dxa_results, species)
         ! d_dxa_const_TRho(i) = d(res(i))/X|T,Rho,X where X = composition

         real(dp), allocatable :: d_dxa_eos(:,:) ! eos internally returns derivs of all quantities

         integer, intent(out) :: ierr ! 0 means AOK.

         type (EoS_General_Info), pointer :: rq

         real(dp) :: X, Y, Z, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx

         call get_eos_ptr(handle,rq,ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for eos -- did you call alloc_eos_handle?'
            return
         end if

         call basic_composition_info( &
            species, chem_id, xa, X, Y, Z, &
            abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)

         allocate(d_dxa_eos(num_eos_basic_results, species))

         call Get_eosPT_Results( &
                  rq, Z, X, abar, zbar, &
                  species, chem_id, net_iso, xa, &
                  Pgas, log10Pgas, T, log10T, &
                  Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
                  res, d_dlnRho_const_T, d_dlnT_const_Rho, d_dxa_eos, &
                  ierr)
         ! only return 1st two d_dxa results (lnE and lnPgas) to star
         d_dxa_const_TRho(1:num_eos_d_dxa_results,1:species) = d_dxa_eos(1:num_eos_d_dxa_results, 1:species)

      end subroutine eosPT_get
      
      
      ! gamma law eos routines -- ignores radiation (e.g.., P = Pgas only)
      
      
      subroutine eos_gamma_DP_get_ET( &
            abar, rho, P, gamma, energy, T, ierr)
         use const_def, only: avo, kerg
         real(dp), intent(in) :: abar, rho, P, gamma
         real(dp), intent(out) :: energy, T
         integer, intent(out) :: ierr
         ierr = 0      
         energy = (P/rho)/(gamma - 1d0)
         T = (gamma - 1d0)*energy*abar/(avo*kerg)
      end subroutine eos_gamma_DP_get_ET
      
      
      subroutine eos_gamma_DE_get_PT( &
            abar, rho, energy, gamma, P, T, ierr)
         use const_def, only: avo, kerg
         real(dp), intent(in) :: abar, rho, energy, gamma
         real(dp), intent(out) :: P, T
         integer, intent(out) :: ierr
         ierr = 0         
         P = (gamma - 1d0)*energy*rho
         T = (gamma - 1d0)*energy*abar/(avo*kerg)
      end subroutine eos_gamma_DE_get_PT
      
      
      subroutine eos_gamma_DT_get_P_energy( &
            abar, rho, T, gamma, P, energy, ierr)
         use const_def, only: avo, kerg
         real(dp), intent(in) :: abar, rho, T, gamma
         real(dp), intent(out) :: P, energy
         integer, intent(out) :: ierr
         ierr = 0
         P = avo*kerg*rho*T/abar
         energy = (P/rho)/(gamma - 1d0)
      end subroutine eos_gamma_DT_get_P_energy
      
      
      subroutine eos_gamma_PRho_get_T_energy( &
            abar, P, rho, gamma, T, energy, ierr)
         use const_def, only: avo, kerg
         real(dp), intent(in) :: abar, P, rho, gamma
         real(dp), intent(out) :: T, energy
         integer, intent(out) :: ierr
         ierr = 0
         energy = (P/rho)/(gamma - 1d0)
         T = (gamma - 1d0)*energy*abar/(avo*kerg)
      end subroutine eos_gamma_PRho_get_T_energy
      
      
      subroutine eos_gamma_PT_get_rho_energy( &
            abar, P, T, gamma, rho, energy, ierr)
         use const_def, only: avo, kerg
         real(dp), intent(in) :: abar, P, T, gamma
         real(dp), intent(out) :: rho, energy
         integer, intent(out) :: ierr
         ierr = 0
         rho = (P/T)*abar/(avo*kerg)
         energy = (P/rho)/(gamma - 1d0)
      end subroutine eos_gamma_PT_get_rho_energy
      
      
      subroutine eos_gamma_DE_get( &
            handle, abar, energy, log10E, rho, log10Rho, gamma, &
            T, log10T, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
            dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
            dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E, ierr)
         use eos_def
         use eosDE_eval, only: Get_eos_gamma_DE_Results
         integer, intent(in) :: handle            
         real(dp), intent(in) :: abar, energy, log10E, Rho, log10Rho, gamma
         real(dp), intent(out) :: T, log10T
         real(dp), intent(inout), dimension(:) :: &
            res, d_dlnRho_const_T, d_dlnT_const_Rho
         real(dp), intent(out) :: & 
            dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
            dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E
         integer, intent(out) :: ierr
         type (EoS_General_Info), pointer :: rq
         call get_eos_ptr(handle,rq,ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for eos -- did you call alloc_eos_handle?'
            return
         end if
         ! require positive density and energy
         if ((rho .le. 0) .or. (energy .le. 0)) then
            ierr = -1
            return
         endif
         call Get_eos_gamma_DE_Results( &
            rq, abar, energy, log10E, rho, log10Rho, gamma, &
            T, log10T, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
            dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
            dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E, ierr)
      end subroutine eos_gamma_DE_get


      subroutine eos_gamma_PT_get( &
            handle, abar, P, log10P, T, log10T, gamma, &
            rho, log10Rho, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
            ierr)
         use eos_def
         use eosDE_eval, only: Get_eos_gamma_DE_Results
         use math_lib
         integer, intent(in) :: handle            
         real(dp), intent(in) :: abar, P, log10P, T, log10T, gamma
         real(dp), intent(out) :: rho, log10Rho
         real(dp), intent(inout), dimension(:) :: &
            res, d_dlnRho_const_T, d_dlnT_const_Rho
         integer, intent(out) :: ierr
         type (EoS_General_Info), pointer :: rq
         real(dp) :: energy, temp, log10temp, & 
            dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
            dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E
         call get_eos_ptr(handle,rq,ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for eos -- did you call alloc_eos_handle?'
            return
         end if
         ! require positive pressure and temperature
         if ((P .le. 0) .or. (T .le. 0)) then
            ierr = -1
            return
         endif
         call eos_gamma_PT_get_rho_energy( &
            abar, P, T, gamma, rho, energy, ierr)
         log10Rho = log10(rho)
         if (ierr /= 0) return
         call eos_gamma_DE_get( &
            handle, abar, energy, log10(energy), rho, log10Rho, gamma, &
            temp, log10temp, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
            dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
            dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E, ierr)
      end subroutine eos_gamma_PT_get


      subroutine eos_gamma_DT_get( &
            handle, abar, rho, log10Rho, T, log10T, gamma, &
            res, d_dlnRho_const_T, d_dlnT_const_Rho, &
            Pgas, Prad, energy, entropy, ierr)
         use eos_def
         use eosDE_eval, only: Get_eos_gamma_DE_Results
         use math_lib
         integer, intent(in) :: handle            
         real(dp), intent(in) :: abar, rho, log10Rho, T, log10T, gamma
         real(dp), intent(inout), dimension(:) :: &
            res, d_dlnRho_const_T, d_dlnT_const_Rho
         real(dp), intent(out) :: Pgas, Prad, energy, entropy
         integer, intent(out) :: ierr
         type (EoS_General_Info), pointer :: rq
         real(dp) :: P, temp, log10temp, & 
            dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
            dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E
         call get_eos_ptr(handle,rq,ierr)
         if (ierr /= 0) then
            write(*,*) 'invalid handle for eos -- did you call alloc_eos_handle?'
            return
         end if
         ! require positive density and temperature
         if ((rho .le. 0) .or. (T .le. 0)) then
            ierr = -1
            return
         endif
         call eos_gamma_DT_get_P_energy( &
            abar, rho, T, gamma, P, energy, ierr)
         if (ierr /= 0) return
         call eos_gamma_DE_get( &
            handle, abar, energy, log10(energy), rho, log10Rho, gamma, &
            temp, log10temp, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
            dlnT_dlnE_c_Rho, dlnT_dlnd_c_E, &
            dlnPgas_dlnE_c_Rho, dlnPgas_dlnd_c_E, ierr)
         Pgas = exp(res(i_lnPgas))
         Prad = crad*T*T*T*T/3d0
         energy = exp(res(i_lnE))
         entropy = exp(res(i_lnS))
      end subroutine eos_gamma_DT_get


      ! misc
      

      subroutine eos_fermi_dirac_integral(dk, eta, theta, fd, fdeta, fdtheta)
         !..from Frank Timmes' site, http://www.cococubed.com/code_pages/fermi_dirac.shtml
         !..this routine computes the fermi-dirac integrals of 
         !..index dk, with degeneracy parameter eta and relativity parameter theta.
         !..input is dk the real(dp) index of the fermi-dirac function,
         !..eta the degeneracy parameter, and theta the relativity parameter.
         !..theta = (k * T)/(mass_electron * c^2), k = Boltzmann const.
         !..the output is fd is computed by applying three 10-point 
         !..gauss-legendre and one 10-point gauss-laguerre rules over
         !..four appropriate subintervals. the derivative with respect to eta is
         !..output in fdeta, and the derivative with respct to theta is in fdtheta.
         !..within each subinterval the fd kernel.
         !..
         !..this routine delivers 14 significant figure accuracy
         !..
         !..reference: j.m. aparicio, apjs 117, 632 1998
         use gauss_fermi, only: dfermi
        
         real(dp), intent(in) :: dk
         real(dp), intent(in) :: eta
         real(dp), intent(in) :: theta
         real(dp), intent(out) :: fd
         real(dp), intent(out) :: fdeta
         real(dp), intent(out) :: fdtheta
         call dfermi(dk, eta, theta, fd, fdeta, fdtheta)
      end subroutine eos_fermi_dirac_integral
      

      subroutine eos_get_helm_results( &
               X, abar, zbar, Rho, log10Rho, T, log10T, &
               coulomb_temp_cut, coulomb_den_cut, &
               include_radiation, include_elec_pos, &
               res, ierr)
         ! direct call to the helm eos.  
         ! returns much more info than the standard

         use eos_def
         use eosDT_eval, only: Get_HELM_Results

         ! INPUT

         real(dp), intent(in) :: X ! the hydrogen mass fraction
            
         real(dp), intent(in) :: abar
            ! mean atomic number (nucleons per nucleus; grams per mole)
         real(dp), intent(in) :: zbar ! mean charge per nucleus
         
         real(dp), intent(in) :: Rho, log10Rho ! the density
            ! provide both if you have them.  
            ! else pass one and set the other to arg_not_provided
            
         real(dp), intent(in) :: T, log10T ! the temperature
            ! provide both if you have them.  
            ! else pass one and set the other to arg_not_provided

         real(dp), intent(in) :: coulomb_temp_cut, coulomb_den_cut

         logical, intent(in) :: include_radiation, include_elec_pos
         
         ! OUTPUT
         
         real(dp), intent(inout) :: res(:) ! (num_helm_results)
            ! array to hold the results
         integer, intent(out) :: ierr ! 0 means AOK.     
         
         logical :: off_table

         call Get_HELM_Results( &
            abar, zbar, Rho, log10Rho, T, log10T, &
            coulomb_temp_cut, coulomb_den_cut, &
            include_radiation, include_elec_pos, &
            res, off_table, ierr)
              
      end subroutine eos_get_helm_results
      
      
      !subroutine eos_convert_helm_results( &
      !      helm_res, Z, X, abar, zbar, Rho, T, res, &
      !      d_dlnRho_const_T, d_dlnT_const_Rho, ierr)      
      subroutine eos_convert_helm_results( &
            helm_res, Z, X, abar, zbar, Rho, T, basic_flag, res,  &
            d_dlnRho_const_T, d_dlnT_const_Rho,  &
            d_dabar_const_TRho, d_dzbar_const_TRho, ierr)      
         use eos_def
         use eos_helm_eval, only: do_convert_helm_results
         real(dp), intent(in) :: helm_res(:) ! (num_helm_results)
         real(dp), intent(in) :: Z, X, abar, zbar, Rho, T
         logical, intent(in) :: basic_flag ! if true, then only want basic results
         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnRho_const_T(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnT_const_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dabar_const_TRho(:) ! (num_eos_basic_results) 
         real(dp), intent(inout) :: d_dzbar_const_TRho(:) ! (num_eos_basic_results) 
         !real(dp), intent(inout), dimension(:) :: d2_dlnd2, d2_dlnd_dlnT, d2_dlnT2
         integer, intent(out) :: ierr      
         d_dabar_const_TRho = 0
         d_dzbar_const_TRho = 0
         call do_convert_helm_results( &
                  helm_res, Z, abar, zbar, Rho, T, &
                  res, d_dlnRho_const_T, d_dlnT_const_Rho, &
                  d_dabar_const_TRho, d_dzbar_const_TRho, &
                  ierr)      
      end subroutine eos_convert_helm_results

      
      ! eosDT search routines.  these use num_lib safe_root to find T or Rho.
      
      subroutine eosDT_get_T( &
               handle, &
               species, chem_id, net_iso, xa, &
               logRho, which_other, other_value, &
               logT_tol, other_tol, max_iter, logT_guess, &
               logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
               logT_result, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               d_dxa_const_TRho, eos_calls, ierr)

         ! finds log10 T given values for density and 'other', and initial guess for temperature.
         ! does up to max_iter attempts until logT changes by less than tol.

         ! 'other' can be any of the basic result variables for the eos
         ! specify 'which_other' by means of the definitions in eos_def (e.g., i_lnE)

         use chem_lib, only: basic_composition_info
         use eos_def
         use eosDT_eval, only : get_T

         integer, intent(in) :: handle ! eos handle; from star, pass s% eos_handle

         integer, intent(in) :: species ! number of species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions

         real(dp), intent(in) :: logRho ! log10 of density
         integer, intent(in) :: which_other ! from eos_def.  e.g., i_lnE
         real(dp), intent(in) :: other_value ! desired value for the other variable
         real(dp), intent(in) :: other_tol

         real(dp), intent(in) :: logT_tol
         integer, intent(in) :: max_iter ! max number of iterations

         real(dp), intent(in) :: logT_guess ! log10 of temperature
         real(dp), intent(in) :: logT_bnd1, logT_bnd2 ! bounds for logT
            ! if don't know bounds, just set to arg_not_provided (defined in const_def)
         real(dp), intent(in) :: other_at_bnd1, other_at_bnd2 ! values at bounds
            ! if don't know these values, just set to arg_not_provided (defined in const_def)

         real(dp), intent(inout) :: logT_result ! log10 of temperature
         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnRho_const_T(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnT_const_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dxa_const_TRho(:,:) ! (num_eos_d_dxa_results, species)
         real(dp), allocatable :: d_dxa_eos(:,:) ! eos internally returns derivs of all quantities

         integer, intent(out) :: eos_calls
         integer, intent(out) :: ierr ! 0 means AOK.

         ! compute composition info
         real(dp) :: Y, Z, X, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx

         call basic_composition_info( &
            species, chem_id, xa, X, Y, Z, &
            abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)

         allocate(d_dxa_eos(num_eos_basic_results, species))

         call get_T( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logRho, which_other, other_value, &
               logT_tol, other_tol, max_iter, logT_guess, &
               logT_bnd1, logT_bnd2,  other_at_bnd1, other_at_bnd2, &
               logT_result, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               d_dxa_eos, eos_calls, ierr)
         ! only return 1st two d_dxa results (lnE and lnPgas) to star
         d_dxa_const_TRho(1:num_eos_d_dxa_results,1:species) = d_dxa_eos(1:num_eos_d_dxa_results, 1:species)

         deallocate(d_dxa_eos)

      end subroutine eosDT_get_T


      subroutine eosDT_get_Rho( &
               handle, &
               species, chem_id, net_iso, xa, &
               logT, which_other, other_value, &
               logRho_tol, other_tol, max_iter, logRho_guess, &
               logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
               logRho_result, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               d_dxa_const_TRho, eos_calls, ierr)

         ! finds log10 Rho given values for temperature and 'other', and initial guess for density.
         ! does up to max_iter attempts until logRho changes by less than tol.

         ! 'other' can be any of the basic result variables for the eos
         ! specify 'which_other' by means of the definitions in eos_def (e.g., i_lnE)

         use chem_lib, only: basic_composition_info
         use eos_def
         use eosDT_eval, only : get_Rho

         integer, intent(in) :: handle ! eos handle; from star, pass s% eos_handle

         integer, intent(in) :: species ! number of species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions

         real(dp), intent(in) :: logT ! log10 of temperature

         integer, intent(in) :: which_other ! from eos_def.  e.g., i_lnE
         real(dp), intent(in) :: other_value ! desired value for the other variable
         real(dp), intent(in) :: other_tol

         real(dp), intent(in) :: logRho_tol

         integer, intent(in) :: max_iter ! max number of Newton iterations

         real(dp), intent(in) :: logRho_guess ! log10 of density
         real(dp), intent(in) :: logRho_bnd1, logRho_bnd2 ! bounds for logRho
            ! if don't know bounds, just set to arg_not_provided (defined in const_def)
         real(dp), intent(in) :: other_at_bnd1, other_at_bnd2 ! values at bounds
            ! if don't know these values, just set to arg_not_provided (defined in const_def)

         real(dp), intent(out) :: logRho_result ! log10 of density

         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnRho_const_T(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnT_const_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dxa_const_TRho(:,:) ! (num_eos_d_dxa_results, species)
         real(dp), allocatable :: d_dxa_eos(:,:) ! eos internally returns derivs of all quantities

         integer, intent(out) :: eos_calls
         integer, intent(out) :: ierr ! 0 means AOK.

         ! compute composition info
         real(dp) :: Y, Z, X, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx

         call basic_composition_info( &
            species, chem_id, xa, X, Y, Z, &
            abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)

         allocate(d_dxa_eos(num_eos_basic_results, species))

         call get_Rho( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logT, which_other, other_value, &
               logRho_tol, other_tol, max_iter, logRho_guess, &
               logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
               logRho_result, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               d_dxa_eos, eos_calls, ierr)
         ! only return 1st two d_dxa results (lnE and lnPgas) to star
         d_dxa_const_TRho(1:num_eos_d_dxa_results,1:species) = d_dxa_eos(1:num_eos_d_dxa_results, 1:species)

         deallocate(d_dxa_eos)

      end subroutine eosDT_get_Rho


      ! eosPT search routines.  these use num_lib safe_root to find T or Pgas.

      subroutine eosPT_get_T( &
               handle, &
               species, chem_id, net_iso, xa, &
               logPgas, which_other, other_value, &
               logT_tol, other_tol, max_iter, logT_guess, &
               logT_bnd1, logT_bnd2, other_at_bnd1, other_at_bnd2, &
               logT_result, Rho, log10Rho, &
               dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               d_dxa_const_TRho, &
               eos_calls, ierr)

         ! finds log10 T given values for gas pressure and 'other',
         ! and initial guess for temperature.
         ! does up to max_iter attempts until logT changes by less than tol.

         ! 'other' can be any of the basic result variables for the eos
         ! specify 'which_other' by means of the definitions in eos_def (e.g., i_lnE)

         use chem_lib, only: basic_composition_info
         use eos_def
         use eosPT_eval, only : get_T

         integer, intent(in) :: handle ! eos handle; from star, pass s% eos_handle

         integer, intent(in) :: species ! number of species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions

         real(dp), intent(in) :: logPgas ! log10 of gas pressure
         integer, intent(in) :: which_other ! from eos_def.  e.g., i_lnE
         real(dp), intent(in) :: other_value ! desired value for the other variable
         real(dp), intent(in) :: other_tol

         real(dp), intent(in) :: logT_tol
         integer, intent(in) :: max_iter ! max number of iterations

         real(dp), intent(in) :: logT_guess ! log10 of temperature
         real(dp), intent(in) :: logT_bnd1, logT_bnd2 ! bounds for logT
            ! if don't know bounds, just set to arg_not_provided (defined in const_def)
         real(dp), intent(in) :: other_at_bnd1, other_at_bnd2 ! values at bounds
            ! if don't know these values, just set to arg_not_provided (defined in const_def)

         real(dp), intent(out) :: logT_result ! log10 of temperature
         real(dp), intent(out) :: Rho, log10Rho ! density
         real(dp), intent(out) :: dlnRho_dlnPgas_const_T
         real(dp), intent(out) :: dlnRho_dlnT_const_Pgas

         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnRho_const_T(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnT_const_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dxa_const_TRho(:,:) ! (num_eos_d_dxa_results, species)
         real(dp), allocatable :: d_dxa_eos(:,:) ! eos internally returns derivs of all quantities

         integer, intent(out) :: eos_calls
         integer, intent(out) :: ierr ! 0 means AOK.

         ! compute composition info
         real(dp) :: Y, Z, X, abar, zbar, z2bar, z53bar, ye, mass_correction, sumx

         call basic_composition_info( &
            species, chem_id, xa, X, Y, Z, &
            abar, zbar, z2bar, z53bar, ye, mass_correction, sumx)

         allocate(d_dxa_eos(num_eos_basic_results, species))

         call get_T( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               logPgas, which_other, other_value, &
               logT_tol, other_tol, max_iter, logT_guess, &
               logT_bnd1, logT_bnd2,  other_at_bnd1, other_at_bnd2, &
               logT_result, Rho, log10Rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas, &
               res, d_dlnRho_const_T, d_dlnT_const_Rho, d_dxa_eos, &
               eos_calls, ierr)
         ! only return 1st two d_dxa results (lnE and lnPgas) to star
         d_dxa_const_TRho(1:num_eos_d_dxa_results,1:species) = d_dxa_eos(1:num_eos_d_dxa_results, 1:species)

         deallocate(d_dxa_eos)

      end subroutine eosPT_get_T
      

      subroutine num_eos_files_loaded( &
            num_DT, num_FreeEOS)
         use eos_def
         integer, intent(out) :: &
            num_DT, num_FreeEOS
         num_DT = count(eosDT_XZ_loaded)
         num_FreeEOS = count(FreeEOS_XZ_loaded)
      end subroutine num_eos_files_loaded


      subroutine eos_get_control_namelist(handle, name, val, ierr)
         use eos_def
         use eos_ctrls_io, only: get_eos_controls
         integer, intent(in) :: handle ! eos handle; from star, pass s% eos_handle
         character(len=*),intent(in) :: name
         character(len=*),intent(out) :: val
         integer, intent(out) :: ierr
         type (EoS_General_Info), pointer :: rq
         ierr = 0
         call get_eos_ptr(handle,rq,ierr)
         if(ierr/=0) return
         call get_eos_controls(rq, name, val, ierr)

      end subroutine eos_get_control_namelist

      subroutine eos_set_control_namelist(handle, name, val, ierr)
         use eos_def
         use eos_ctrls_io, only: set_eos_controls
         integer, intent(in) :: handle ! eos handle; from star, pass s% eos_handle
         character(len=*),intent(in) :: name
         character(len=*),intent(in) :: val
         integer, intent(out) :: ierr
         type (EoS_General_Info), pointer :: rq
         ierr = 0
         call get_eos_ptr(handle,rq,ierr)
         if(ierr/=0) return
         call set_eos_controls(rq, name, val, ierr)

      end subroutine eos_set_control_namelist


      end module eos_lib
