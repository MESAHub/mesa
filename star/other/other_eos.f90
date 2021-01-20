! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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
 
      module other_eos

      ! set use_other_eos = .true. in your controls namelist
        
      ! edit the extras_controls routine to set the procedure pointers
      ! e.g.,
         ! s% other_eosDT_get => my_eosDT_get
         ! s% other_eosDT_get_T => my_eosDT_get_T
         ! s% other_eosDT_get_Rho => my_eosDT_get_Rho


      use star_def
      use eos_def

      implicit none
      
            
      contains
      
      
      subroutine null_other_eosDT_get( &
              id, k, handle, &
              species, chem_id, net_iso, xa, &
              Rho, log10Rho, T, log10T, & 
              res, d_dlnRho_const_T, d_dlnT_const_Rho, d_dxa_const_TRho, ierr)

         ! INPUT
         use chem_def, only: num_chem_isos
         
         integer, intent(in) :: id ! star id if available; 0 otherwise
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         integer, intent(in) :: handle ! eos handle

         integer, intent(in) :: species
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
         ! d_dlnRho(i) = d(res(i))/dlnd|T
         real(dp), intent(inout) :: d_dlnT_const_Rho(:) ! (num_eos_basic_results) 
         ! d_dlnT(i) = d(res(i))/dlnT|Rho
         real(dp), intent(inout) :: d_dxa_const_TRho(:,:) ! (num_eos_d_dxa_results,species)
         ! d_dxa(i,j) = d(res(i))/dxa(j)|T,Rho
         
         integer, intent(out) :: ierr ! 0 means AOK.
         
         res = 0
         d_dlnRho_const_T = 0
         d_dlnT_const_Rho = 0
         d_dxa_const_TRho = 0

         write(*,*) 'no implementation for other_eosDT_get'
         ierr = -1
         
      end subroutine null_other_eosDT_get
      
      
      subroutine null_other_eosDE_get( &
               handle, Z, X, abar, zbar, &
               species, chem_id, net_iso, xa, &
               energy, log10E, rho, log10Rho, log10T_guess, &
               T, log10T, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               d_dabar_const_TRho, d_dzbar_const_TRho, ierr)

         use eos_def

         ! INPUT
         
         integer, intent(in) :: handle

         real(dp), intent(in) :: Z ! the metals mass fraction
         real(dp), intent(in) :: X ! the hydrogen mass fraction
            
         real(dp), intent(in) :: abar
            ! mean atomic number (nucleons per nucleus; grams per mole)
         real(dp), intent(in) :: zbar ! mean charge per nucleus
         
         integer, intent(in) :: species
         integer, pointer :: chem_id(:) ! maps species to chem id
            ! index from 1 to species
            ! value is between 1 and num_chem_isos         
         integer, pointer :: net_iso(:) ! maps chem id to species number
            ! index from 1 to num_chem_isos (defined in chem_def)
            ! value is 0 if the iso is not in the current net
            ! else is value between 1 and number of species in current net
         real(dp), intent(in) :: xa(:) ! mass fractions
         
         real(dp), intent(in) :: energy, log10E ! the internal energy
            ! provide both if you have them.  else pass one and set the other to arg_not_provided
            ! "arg_not_provided" is defined in mesa const_def
            
         real(dp), intent(in) :: Rho, log10Rho ! the density
            ! provide both if you have them.  else pass one and set the other to arg_not_provided
            
         real(dp), intent(in) :: log10T_guess ! guess for logT to use if off table
                     
         ! OUTPUT
         
         real(dp), intent(out) :: T, log10T
         
         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
         
         ! partial derivatives of the basic results
         
         real(dp), intent(inout) :: d_dlnRho_const_T(:) ! (num_eos_basic_results) 
         ! d_dlnRho_const_T(i) = d(res(i))/dlnd|T,X where X = composition
         real(dp), intent(inout) :: d_dlnT_const_Rho(:) ! (num_eos_basic_results) 
         ! d_dlnT(i) = d(res(i))/dlnT|Rho,X where X = composition

         real(dp), intent(inout) :: d_dabar_const_TRho(:) ! (num_eos_basic_results) 
         ! d_dabar(i) = d(res(i))/dabar|TRho,zbar
         real(dp), intent(inout) :: d_dzbar_const_TRho(:) ! (num_eos_basic_results) 
         ! d_dzbar(i) = d(res(i))/dzbar|TRho,abar
         
         integer, intent(out) :: ierr ! 0 means AOK.

         ! NOTE: when converting partials for f = f(lnd,lnT(lnd,lnE)), 
         ! df_dlnE_const_Rho = df_dlnT_const_Rho*dlnT_dlnE_const_Rho
         !     dlnT_dlnE_const_Rho = E/(Cv*T)
         ! df_dlnd_const_E = df_dlnd_const_T + df_dlnT_const_Rho*dlnT_dlnd_const_E
         !     dlnT_dlnd_const_E = -Rho*dE_dRho/(Cv*T)
         
         T = 0
         log10T = 0
         res = 0
         d_dlnRho_const_T = 0
         d_dlnT_const_Rho = 0
         d_dabar_const_TRho = 0
         d_dzbar_const_TRho = 0

         write(*,*) 'no implementation for other_eosDE_get'
         ierr = -1
         
      end subroutine null_other_eosDE_get
      

      ! eosDT search routines.
      
      subroutine null_other_eosDT_get_T( &
               id, k, handle, &
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
         
         use chem_def, only: num_chem_isos

         integer, intent(in) :: id ! star id if available; 0 otherwise
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         integer, intent(in) :: handle

         integer, intent(in) :: species
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

         real(dp), intent(out) :: logT_result ! log10 of temperature
         real(dp), intent(inout) :: res(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dlnRho_const_T(:) ! (num_eos_basic_results) 
         real(dp), intent(inout) :: d_dlnT_const_Rho(:) ! (num_eos_basic_results)
         real(dp), intent(inout) :: d_dxa_const_TRho(:,:) ! (num_eos_d_dxa_results, species)
         
         integer, intent(out) :: eos_calls
         integer, intent(out) :: ierr ! 0 means AOK.
         
         logT_result = 0
         res = 0
         d_dlnRho_const_T = 0
         d_dlnT_const_Rho = 0
         d_dxa_const_TRho = 0
         eos_calls = 0

         write(*,*) 'no implementation for other_eosDT_get_T'
         ierr = -1

      end subroutine null_other_eosDT_get_T
      

      subroutine null_other_eosDT_get_Rho( &
               id, k, handle, &
               species, chem_id, net_iso, xa, &
               logT, which_other, other_value, &
               logRho_tol, other_tol, max_iter, logRho_guess,  &
               logRho_bnd1, logRho_bnd2, other_at_bnd1, other_at_bnd2, &
               logRho_result, res, d_dlnRho_const_T, d_dlnT_const_Rho, &
               d_dxa_const_TRho, eos_calls, ierr)
     
         ! finds log10 Rho given values for temperature and 'other', and initial guess for density.
         ! does up to max_iter attempts until logRho changes by less than tol.
         
         ! 'other' can be any of the basic result variables for the eos
         ! specify 'which_other' by means of the definitions in eos_def (e.g., i_lnE)
         
         use chem_def, only: num_chem_isos
         
         integer, intent(in) :: id ! star id if available; 0 otherwise
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         integer, intent(in) :: handle

         integer, intent(in) :: species
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

         integer, intent(out) :: eos_calls
         integer, intent(out) :: ierr ! 0 means AOK.
         
         logRho_result = 0
         res = 0
         d_dlnRho_const_T = 0
         d_dlnT_const_Rho = 0
         d_dxa_const_TRho = 0
         eos_calls = 0

         write(*,*) 'no implementation for other_eosDT_get_Rho'
         ierr = -1

      end subroutine null_other_eosDT_get_Rho
      

      end module other_eos
      
      
      
      
